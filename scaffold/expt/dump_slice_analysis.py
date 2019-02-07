import os
from logging import getLogger

import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib import colors
import iris

has_metpy = False
try:
    import metpy
    import metpy.calc as mpcalc
    from metpy.units import units
    has_metpy = True
except ImportError:
    pass

from omnium import Analyser, OmniumError
from omnium.consts import Re, p_ref, kappa
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
    input_filename_glob = '{input_dir}/atmosa_da480.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.{runid:03}.dump_slice_analysis.dummy']
    uses_runid = True
    runid_pattern = 'atmosa_da(?P<runid>\d{3}).nc'

    def _create_cube(self, archetype, data, name, units):
        cube = archetype.copy()
        cube.rename(name)
        cube.units = units
        cube.data = data
        return cube

    def load(self):
        self.load_cubes()

    def run(self):
        dump = self.cubes
        self.rho = get_cube(dump, 0, 253) / Re ** 2
        self.rho_d = get_cube(dump, 0, 389)

        self.theta = get_cube(dump, 0, 4)
        self.exnerp = get_cube(dump, 0, 255)

        self.q = get_cube(dump, 0, 10)
        self.qcl = get_cube(dump, 0, 254)
        self.qcf = get_cube(dump, 0, 12)
        self.qrain = get_cube(dump, 0, 272)
        self.qgraup = get_cube(dump, 0, 273)

        self.u = get_cube(dump, 0, 2)
        self.w = get_cube(dump, 0, 150)

        try:
            qcf2 = get_cube(dump, 0, 271)
            self.qcf2 = qcf2
        except OmniumError:
            logger.info('dump has no qcf2')

        qvdata = get_cube(dump, 0, 10).data
        Tdata = self.theta.data * self.exnerp.data
        pdata = self.exnerp.data ** (1 / kappa) * p_ref

        p = pdata * units('Pa')
        qv = qvdata * units('kg/kg')
        T = Tdata * units('K')
        Td = mpcalc.dewpoint_from_specific_humidity(qv, T, p)
        theta_e = mpcalc.equivalent_potential_temperature(p, T, Td)

        self.theta_e = self._create_cube(self.theta, theta_e.magnitude, 'theta_e', 'K')

    def display_results(self):
        self.vertlevs = VertLev(self.suite.suite_dir)
        self._plot(self.task.expt)

    def _create_dir_savefig(self, filename):
        full_filename = self.file_path(filename)
        os.makedirs(os.path.dirname(full_filename), exist_ok=True)
        plt.savefig(full_filename)

    def _plot(self, expt):
        # self._w_slices_plots(expt)
        # self._w_mean_y_plots(expt)
        # self._th_slices_plots(expt)
        # self._u_slices_plots(expt)
        # self._qvar_slices_y_plots(expt)
        # self._qvar_mean_y_plots(expt)
        # self._th_e_slices_plots(expt)
#
        self._plot_slices(expt, 'w', use_norm=True)
        self._plot_slices(expt, 'theta', use_norm=True, anomaly=True)
        self._plot_slices(expt, 'theta_e', use_norm=True, anomaly=True, cmap='RdBu_r')
        qvars = ['qcl', 'qcf', 'qcf2', 'qrain', 'qgraup']
        for qvar in qvars:
            if not hasattr(self, qvar):
                continue
            self._plot_slices(expt, qvar, cmap='Blues')

    def _plot_slices(self, expt, var, **kwargs):
        cube = getattr(self, var)
        z = cube.coord('atmosphere_hybrid_height_coordinate').points / 1000
        self._plot_xy_slices(expt, var, cube, z, **kwargs)
        self._plot_vert_slices(expt, var, cube, 'xz', **kwargs)
        self._plot_vert_slices(expt, var, cube, 'yz', **kwargs)

        self._plot_vert_mean_slice(expt, var, cube, 'xz', **kwargs)

    def _plot_xy_slices(self, expt, var, cube, z, use_norm=False, cmap='bwr',
                        anomaly=False):
        for i in range(cube.shape[0]):
            fig, ax = plt.subplots(dpi=100)
            data = cube.data[i]
            title = '{} xy slice at z={:.2f} km'.format(var, z[i])
            if anomaly:
                data = data - data.mean()
                title = 'anom ' + title
            ax.set_title(title)
            # Coords are model_level, y, x or model_level, lat, lon
            kwargs = {}
            if use_norm:
                kwargs['norm'] = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            if cmap:
                kwargs['cmap'] = cmap
            im = ax.imshow(data, origin='lower', **kwargs)

            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            plt.colorbar(im)
            self._create_dir_savefig('/xy/{2}/{0}_{1}_{2}_slice_{3}.png'
                                     .format(expt, self.task.runid, var, i))
            plt.close('all')

    def _plot_vert_slices(self, expt, var, cube, plane='xy', use_norm=False, cmap='bwr',
                          anomaly=False):
        # Coords are model_level, y, x or model_level, lat, lon
        if plane == 'xz':
            N = cube.shape[1]
        elif plane == 'yz':
            N = cube.shape[2]
        else:
            raise ValueError('plane=[xz|yz]')

        for i in range(N):
            fig, ax = plt.subplots(dpi=100)
            ax.set_ylabel('height (km)')
            filename = ('/{0}/{3}/{1}_{2}_{3}_slice_{4}.png'
                        .format(plane, expt, self.task.runid, var, i))
            if plane == 'xz':
                ax.set_xlabel('x (km)')
                title = '{} xz slice at y={} gridbox'.format(var, i)
                data = cube.data[:, i]
            elif plane == 'yz':
                ax.set_xlabel('y (km)')
                title = '{} yz slice at x={} gridbox'.format(var, i)
                data = cube.data[:, :, i]
            if anomaly:
                data = data - data.mean(axis=1)[:, None]
                title = 'anom ' + title
            ax.set_title(title)

            self._plot_vert_slice(ax, data, N, use_norm, cmap)
            self._create_dir_savefig(filename)
            plt.close('all')

    def _plot_vert_slice(self, ax, data, N, use_norm, cmap='bwr'):
        data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(N),
                                                         data)
        data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, N - 1, N))
        kwargs = {}
        if use_norm:
            kwargs['norm'] = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
        if cmap:
            kwargs['cmap'] = cmap

        im = ax.imshow(data_interp[:200], origin='lower', aspect=10, extent=(0, 256, 0, 20),
                       **kwargs)

        plt.colorbar(im)

    def _plot_vert_mean_slice(self, expt, var, cube, plane='xy', use_norm=False, cmap='',
                              anomaly=False):
        # Coords are model_level, y, x or model_level, lat, lon
        if plane == 'xz':
            cube_index = 1
        elif plane == 'yz':
            cube_index = 2
        else:
            raise ValueError('plane=[xz|yz]')
        N = cube.shape[cube_index]

        fig, ax = plt.subplots(dpi=100)
        ax.set_ylabel('height (km)')
        filename = ('/{0}/{3}/{1}_{2}_{3}_mean.png'
                    .format(plane, expt, self.task.runid, var))
        data = cube.data.mean(axis=cube_index)
        if plane == 'xz':
            ax.set_xlabel('x (km)')
            title = '{} xz mean'.format(var)
        elif plane == 'yz':
            ax.set_xlabel('y (km)')
            title = '{} yz mean'.format(var)
        if anomaly:
            data = data - data.mean(axis=1)[:, None]
            title = 'anom ' + title
        ax.set_title(title)

        # Coords are model_level, y, x or model_level, lat, lon
        self._plot_vert_slice(ax, data, N, use_norm, cmap)
        self._create_dir_savefig(filename)
        plt.close('all')

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
            # Use for equal aspect ratio.
            # im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)
            im = ax.imshow(data_interp[:200],
                           norm=norm, origin='lower', cmap='bwr',
                           aspect=10,
                           extent=(0, 256, 0, 20))

            ax.set_title('w xz slice at y={} gridbox'.format(i))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('height (km)')
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

        mean_wind_min_index = ucube.data.mean(axis=(1, 2)).argmin()
        mean_wind = ucube.data[mean_wind_min_index].mean()

        for i in range(ucube.shape[0]):
            fig, ax = plt.subplots(dpi=100)
            data = ucube.data[i] - mean_wind
            # Coords are model_level, y, x or model_level, lat, lon
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data, norm=norm, origin='lower', cmap='bwr')

            ax.set_title('u - u_mean={:.2f} m/s at z={:.2f} km xy slice at z={:.2f} km'.format(mean_wind,
                                                                                               z[mean_wind_min_index],
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
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_rho, np.arange(Nx),
                                                             data)
            data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)

            ax.set_title('u - u_mean={:.2f} at z={:.2f} km xz slice at y={} gridbox'.format(mean_wind,
                                                                                            z[mean_wind_min_index],
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

    def _th_e_slices_plots(self, expt):
        thcube = self.th
        z = thcube.coord('atmosphere_hybrid_height_coordinate').points / 1000

        exnerp = self.ep
        theta = self.th

        qvdata = get_cube(self.cubes, 0, 10).data
        Tdata = theta.data * exnerp.data
        pdata = exnerp.data ** (1 / kappa) * p_ref

        p = pdata * units('Pa')
        qv = qvdata * units('kg/kg')
        T = Tdata * units('K')
        Td = mpcalc.dewpoint_from_specific_humidity(qv, T, p)
        theta_e = mpcalc.equivalent_potential_temperature(p, T, Td)

        for i in range(theta_e.shape[0]):
            fig, ax = plt.subplots(dpi=100)
            # N.B. subtract mean.
            data = theta_e[i] - theta_e[i].mean()
            # Coords are model_level, y, x or model_level, lat, lon
            norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            im = ax.imshow(data, norm=norm, origin='lower', cmap='bwr')

            ax.set_title('theta xy slice at z={:.2f} km'.format(z[i]))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xy/{}_{}_theta_e_slice_{}.png'.format(expt,
                                                                               self.task.runid,
                                                                               i)))
            plt.close('all')

        for i in range(theta_e.shape[1]):
            fig, ax = plt.subplots(dpi=100)
            data = theta_e[:, i]
            data[np.isnan(data)] = 273 * units('K')
            # Coords are model_level, y, x or model_level, lat, lon
            Nx = theta_e.shape[2]
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Nx),
                                                             data)
            data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
            # norm = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
            # Use for equal aspect ratio.
            # im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)
            im = ax.imshow(data_interp[:200],
                           origin='lower', cmap='RdBu_r',
                           aspect=10,
                           vmin=325, vmax=370,
                           extent=(0, 256, 0, 20))

            ax.set_title('theta_e xz slice at y={} gridbox'.format(i))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('height (km)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xz/{}_{}_theta_e_slice_{}.png'.format(expt,
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
        im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr',
                       extent=(0, 256, 0, 20), aspect=10)

        ax.set_title('w mean over y')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('height (km)')
        plt.colorbar(im)
        plt.savefig(self.file_path('/xz/{}_{}_w_mean_over_y.png'.format(expt,
                                                                        self.task.runid)))
        plt.close('all')

    def _qvar_slices_y_plots(self, expt):
        qvars = ['qcl', 'qcf', 'qcf2', 'qrain', 'qgraup']
        for qvar in qvars:
            if not hasattr(self, qvar):
                continue
            qcube = getattr(self, qvar)
            z = qcube.coord('atmosphere_hybrid_height_coordinate').points / 1000

            for i in range(qcube.shape[0]):
                fig, ax = plt.subplots(dpi=100)
                data = qcube.data[i]
                # Coords are model_level, y, x or model_level, lat, lon
                im = ax.imshow(data, origin='lower', cmap='Blues')

                ax.set_title('q xy slice at z={:.2f} km'.format(z[i]))
                ax.set_xlabel('x (km)')
                ax.set_ylabel('y (km)')
                plt.colorbar(im)
                plt.savefig(self.file_path('/xy/{}_{}_{}_slice_{}.png'.format(expt,
                                                                              self.task.runid,
                                                                              qvar,
                                                                              i)))
                plt.close('all')

            for i in range(qcube.shape[1]):
                fig, ax = plt.subplots(dpi=100)
                data = qcube.data[:, i]
                # Coords are model_level, y, x or model_level, lat, lon
                Nx = qcube.shape[2]
                data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(Nx),
                                                                 data)
                data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
                # Use for equal aspect ratio.
                # im = ax.imshow(data_interp[:200], norm=norm, origin='lower', cmap='bwr', aspect=0.1)
                im = ax.imshow(data_interp[:200],
                               origin='lower', cmap='Blues',
                               aspect=10,
                               extent=(0, 256, 0, 20))

                ax.set_title('q xz slice at y={} gridbox'.format(i))
                ax.set_xlabel('x (km)')
                ax.set_ylabel('height (km)')
                plt.colorbar(im)
                plt.savefig(self.file_path('/xz/{}_{}_{}_slice_{}.png'.format(expt,
                                                                              self.task.runid,
                                                                              qvar,
                                                                              i)))
                plt.close('all')

    def _qvar_mean_y_plots(self, expt):
        qvars = ['qcl', 'qcf', 'qcf2', 'qrain', 'qgraup']
        for qvar in qvars:
            if not hasattr(self, qvar):
                continue
            qcube = getattr(self, qvar)
            z = qcube.coord('atmosphere_hybrid_height_coordinate').points / 1000

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
            im = ax.imshow(data_interp[:200], origin='lower', cmap='Blues',
                           extent=(0, 256, 0, 20), aspect=10)

            ax.set_title('{} mean over y'.format(qvar))
            ax.set_xlabel('x (km)')
            ax.set_ylabel('height (km)')
            plt.colorbar(im)
            plt.savefig(self.file_path('/xz/{}_{}_{}_mean_over_y.png'.format(expt,
                                                                             self.task.runid,
                                                                             qvar)))
            plt.close('all')

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
