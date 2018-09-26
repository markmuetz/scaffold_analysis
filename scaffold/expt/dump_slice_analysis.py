import os
from logging import getLogger

import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib import colors

from omnium import Analyser, OmniumError
from omnium.consts import Re, L, cp, g
from omnium.utils import get_cube

from scaffold.vertlev import VertLev

logger = getLogger('scaf.dump_slice')

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


class DumpSliceAnalyser(Analyser):
    analysis_name = 'dump_slice_analysis'
    single_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da???.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.dump_slice_analysis.dummy']
    uses_runid = True
    runid_pattern = 'atmosa_da(?P<runid>\d{3}).nc'

    def load(self):
        self.load_cubes()

    def run(self):
        dump = self.cubes
        self.rho = get_cube(dump, 0, 253) / Re ** 2
        self.rho_d = get_cube(dump, 0, 389)

        self.th = get_cube(dump, 0, 4)
        self.ep = get_cube(dump, 0, 255)

        self.q = get_cube(dump, 0, 10)
        self.qcl = get_cube(dump, 0, 254)
        self.qcf = get_cube(dump, 0, 12)
        self.qrain = get_cube(dump, 0, 272)
        self.qgraup = get_cube(dump, 0, 273)
        self.w = get_cube(dump, 0, 150)

        try:
            qcf2 = get_cube(dump, 0, 271)
            self.qcf2 = qcf2
        except OmniumError:
            logger.info('dump has no qcf2')

    def display_results(self):
        os.makedirs(self.file_path('/xy'), exist_ok=True)
        os.makedirs(self.file_path('/xz'), exist_ok=True)
        os.makedirs(self.file_path('/yz'), exist_ok=True)
        self.vertlevs = VertLev(self.suite.suite_dir)
        self._plot(self.task.expt)

    def _plot(self, expt):
        # self._w_slices_plots(expt)
        # self._w_mean_y_plots(expt)
        self._th_slices_plots(expt)
        self._u_slices_plots(expt)
        # self._qvar_mean_y_plots(expt)

    def _w_slices_plots(self, expt):
        wcube = self.w
        z = wcube.coord('atmosphere_hybrid_height_coordinate').points / 1000

        for i in range(wcube.shape[0]):
            fig, ax = plt.subplots(dpi=100)
            data = wcube.data[i]
            # Coords are model_level, y, x or model_level, lat, lon
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data, norm=norm, origin='lower', cmap='bwr')

            ax.set_title('w xy slice at z={:.2f} km'.format(z[i]))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xy/{}_{}_w_slice_{}.png'.format(expt,
                                                                         self.task.runid,
                                                                         i)))
            plt.close('all')

        for i in range(wcube.shape[1]):
            fig, ax = plt.subplots(dpi=100)
            data = wcube.data[:, i]
            # Coords are model_level, y, x or model_level, lat, lon
            Nx = wcube.shape[2]
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Nx),
                                                             data)
            data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)

            ax.set_title('w xz slice at y={} gridbox'.format(i))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('height (100 m)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xz/{}_{}_w_slice_{}.png'.format(expt,
                                                                         self.task.runid,
                                                                         i)))
            plt.close('all')

        for i in range(wcube.shape[2]):
            fig, ax = plt.subplots(dpi=100)
            data = wcube.data[:, :, i]
            # Coords are model_level, y, x or model_level, lat, lon
            Ny = wcube.shape[1]
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Ny),
                                                             data)
            data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Ny - 1, Ny))
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)

            ax.set_title('w yz slice at x={} gridbox'.format(i))
            ax.set_xlabel('y (km)')
            ax.set_ylabel('height (100 m)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/yz/{}_{}_w_slice_{}.png'.format(expt,
                                                                         self.task.runid,
                                                                         i)))
            plt.close('all')

    def _u_slices_plots(self, expt):
        ucube = self.u
        z = ucube.coord('atmosphere_hybrid_height_coordinate').points / 1000

        height_level = 18
        mean_wind = ucube.data[height_level].mean()

        for i in range(ucube.shape[0]):
            fig, ax = plt.subplots(dpi=100)
            data = ucube.data[i] - mean_wind
            # Coords are model_level, y, x or model_level, lat, lon
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data, norm=norm, origin='lower', cmap='bwr')

            ax.set_title('u - u_mean={:.2f} m/s at z={:.2f} xy slice at z={:.2f} km'.format(mean_wind,
                                                                                            z[height_level],
                                                                                            z[i]))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xy/{}_{}_u_slice_{}.png'.format(expt,
                                                                         self.task.runid,
                                                                         i)))
            plt.close('all')

        for i in range(ucube.shape[1]):
            fig, ax = plt.subplots(dpi=100)
            data = ucube.data[:, i] - mean_wind
            # Coords are model_level, y, x or model_level, lat, lon
            Nx = ucube.shape[2]
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Nx),
                                                             data)
            data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)

            ax.set_title('u - u_mean={:.2f} at z={:.2f} xz slice at y={} gridbox'.format(mean_wind,
                                                                                         z[height_level],
                                                                                         i))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('height (100 m)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xz/{}_{}_u_slice_{}.png'.format(expt,
                                                                         self.task.runid,
                                                                         i)))
            plt.close('all')

    def _th_slices_plots(self, expt):
        thcube = self.th
        z = thcube.coord('atmosphere_hybrid_height_coordinate').points / 1000

        for i in range(thcube.shape[0]):
            fig, ax = plt.subplots(dpi=100)
            # N.B. subtract mean.
            data = thcube.data[i] - thcube.data[i].mean()
            # Coords are model_level, y, x or model_level, lat, lon
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data, norm=norm, origin='lower', cmap='bwr')

            ax.set_title('theta xy slice at z={:.2f} km'.format(z[i]))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xy/{}_{}_theta_slice_{}.png'.format(expt,
                                                                             self.task.runid,
                                                                             i)))
            plt.close('all')

    def _w_mean_y_plots(self, expt):
        wcube = self.w
        fig, ax = plt.subplots(dpi=100)
        data = wcube.data
        # Coords are model_level, y, x or model_level, lat, lon
        data_mean = data.mean(axis=1)
        Nx = data.shape[2]
        data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Nx),
                                                         data_mean)
        data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
        # Only go up to 20 km and use aspect ratio to plot equal aspect
        # (allowing for diff in coords).
        norm = MidpointNormalize(midpoint=0,
                                 vmin=data_interp.min(),
                                 vmax=data_interp.max())
        im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)

        ax.set_title('w mean over y')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('height (100 m)')
        plt.colorbar(im)
        plt.savefig(self.file_path('/xz/{}_{}_w_mean_over_y.png'.format(expt,
                                                                        self.task.runid)))
        plt.close('all')

    def _qvar_mean_y_plots(self, expt):
        qvars = ['qcl', 'qcf', 'qrain', 'qgraup']
        for qvar in qvars:
            if not hasattr(self, qvar):
                continue
            qcube = getattr(self, qvar)

            fig, ax = plt.subplots(dpi=100)
            data = qcube.data
            # Coords are model_level, y, x or model_level, lat, lon
            data_mean = data.mean(axis=1)
            Nx = data.shape[2]
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Nx),
                                                             data_mean)
            data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
            # Only go up to 20 km and use aspect ratio to plot equal aspect
            # (allowing for diff in coords).
            im = ax.imshow(data_interp[:200], origin='lower', cmap='Blues', aspect=0.1)

            ax.set_title('{} mean over y'.format(qvar))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('height (100 m)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xz/{}_{}_{}_mean_over_y.png'.format(expt,
                                                                             self.task.runid,
                                                                             qvar)))
            plt.close('all')

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
