import os
from logging import getLogger

import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib import colors

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
    input_filename_glob = '{input_dir}/atmosa_da4??.nc'
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
        # self._plot_slices(expt, 'u', use_norm=True, anomaly=True, vlev='rho')
        # self._plot_slices(expt, 'u', use_norm=True, use_mean_wind=True, vlev='rho')
        # self._plot_slices(expt, 'w', use_norm=True)
        # self._plot_slices(expt, 'theta', use_norm=True, anomaly=True)
        # self._plot_slices(expt, 'theta_e', cmap='RdBu_r')
        # qvars = ['qcl', 'qcf', 'qcf2', 'qrain', 'qgraup']
        # for qvar in qvars:
        #     if not hasattr(self, qvar):
        #         continue
        #     self._plot_slices(expt, qvar, cmap='Blues')

        expt_slices = {
            'S0W0Forced': {480: [
                (self._plot_indiv_cross_sections, (expt, 'w', 10, 125, 50), {'use_norm': True}),
                (self._plot_indiv_cross_sections, (expt, 'qcl', 10, 125, 50), {'cmap': 'Blues'}),
                (self._plot_indiv_cross_sections, (expt, 'qcf', 10, 125, 50), {'cmap': 'Blues'}),
                (self._plot_indiv_cross_sections, (expt, 'theta_e', 10, 125, 50),
                 {'cmap': 'RdBu_r'}),
                (self._plot_zoom, (expt, 'w', 10, 5, 15, 125, 120, 130, [50]), {'use_norm': True}),
                (self._plot_zoom, (expt, 'theta_e', 10, 5, 15, 125, 120, 130, [50]), {'cmap': 'RdBu_r'}),
                (self._plot_zoom, (expt, 'qcl', 10, 5, 15, 125, 120, 130, [50]), {'cmap': 'Blues'}),
                (self._plot_zoom, (expt, 'qcf', 10, 5, 15, 125, 120, 130, [50]), {'cmap': 'Blues'}),
                (self._plot_zoom, (expt, 'qcf2', 10, 5, 15, 125, 120, 130, [50]), {'cmap': 'Blues'}),
                (self._plot_zoom, (expt, 'qrain', 10, 5, 15, 125, 120, 130, [50]), {'cmap': 'Blues'}),
                (self._plot_zoom, (expt, 'qgraup', 10, 5, 15, 125, 120, 130, [50]), {'cmap': 'Blues'}),
            ]},
            'S4W5rel2_S4W5Forced': {480: [
                (self._plot_zoom, (expt, 'w', 0, 100, 185, 225, [25, 50]),
                 {'use_norm': True, 'aspect':4}),
                (self._plot_zoom, (expt, 'qcl', 0, 100, 185, 225, [25, 50]),
                 {'cmap': 'Blues', 'aspect':4}),
                (self._plot_zoom, (expt, 'qcf', 0, 100, 185, 225, [25, 50]),
                 {'cmap': 'Blues', 'aspect':4}),
                (self._plot_zoom, (expt, 'qcf2', 0, 100, 185, 225, [25, 50]),
                 {'cmap': 'Blues', 'aspect':4}),
                (self._plot_zoom, (expt, 'qrain', 0, 100, 185, 225, [25, 50]),
                 {'cmap': 'Blues', 'aspect':4}),
                (self._plot_zoom, (expt, 'qgraup', 0, 100, 185, 225, [25, 50]),
                 {'cmap': 'Blues', 'aspect':4}),
                (self._plot_zoom, (expt, 'theta', 0, 100, 185, 225, [25, 50]),
                 {'use_norm': True, 'anomaly': True, 'aspect':4}),
            ]},
        }
        if expt in expt_slices:
            if self.task.runid in expt_slices[expt]:
                for fn, args, kwargs in expt_slices[expt][self.task.runid]:
                    fn(*args, **kwargs)

    def _plot_zoom(self, expt, var, i, imin, imax, j, jmin, jmax, ks, **kwargs):
        extent = (imin, imax, jmin, jmax)

        for k in ks:
            # self._plot_zoom_hor_slice(expt, var, imin, imax, jmin, jmax, k, **kwargs)
            self._plot_indiv_hor_slice(expt, var, i, j, k, mode='zoom', extent=extent, **kwargs)
        self._plot_indiv_vert_slice(expt, var, 'xz', i, j, None, mode='zoom', extent=extent, **kwargs)
        self._plot_indiv_vert_slice(expt, var, 'yz', i, j, None, mode='zoom', extent=extent, **kwargs)


    def _plot_indiv_cross_sections(self, expt, var, i, j, k, **kwargs):
        self._plot_indiv_hor_slice(expt, var, i, j, k, mode='indiv', **kwargs)
        self._plot_indiv_vert_slice(expt, var, 'xz', i, j, k, mode='indiv', **kwargs)
        self._plot_indiv_vert_slice(expt, var, 'yz', i, j, k, mode='indiv', **kwargs)

    def _plot_zoom_vert_slice(self, expt, var, plane, imin, imax, jmin, jmax, ks,
                              use_norm=False, cmap='bwr', anomaly=False, aspect=1):
        cube = getattr(self, var)
        z = cube.coord('atmosphere_hybrid_height_coordinate').points / 1000
        # Coords are model_level, y, x or model_level, lat, lon
        if plane == 'xz':
            N = imax - imin
            index = int((jmax + jmin) / 2)
            #vline_index = i
        elif plane == 'yz':
            N = jmax - jmin
            index = int((imax + imin) / 2)
            #vline_index = j
        else:
            raise ValueError('plane=[xz|yz]')

        if anomaly:
            var = var + '_anom'

        fig, ax = plt.subplots(dpi=100)
        ax.set_ylabel('height (km)')
        filename = ('/{2}/{0}/zoom_{3}/{1}_{2}_{3}_slice_{4}.png'
                    .format(plane, expt, self.task.runid, var, index))
        if plane == 'xz':
            ax.set_xlabel('x (km)')
            title = '{} xz slice at y={} gridbox'.format(var, index)
            data = cube.data[:, index, imin:imax]
            extent = (imin, imax, 0, 20)
        elif plane == 'yz':
            ax.set_xlabel('y (km)')
            title = '{} yz slice at x={} gridbox'.format(var, index)
            data = cube.data[:, jmin:jmax, index]
            extent = (jmin, jmax, 0, 20)
        if anomaly:
            data = data - data.mean(axis=1)[:, None]
        ax.set_title(title)

        self._plot_vert_slice(ax, data, N, use_norm, cmap, extent=extent, aspect=aspect)
        self._create_dir_savefig(filename)
        plt.close('all')

    def _plot_indiv_hor_slice(self, expt, var, i, j, k, use_norm=False, cmap='bwr',
                              anomaly=False, use_mean_wind=False, mode='', extent=None,
                              **kwargs):
        cube = getattr(self, var)
        z = cube.coord('atmosphere_hybrid_height_coordinate').points / 1000
        if anomaly:
            var = var + '_anom'
        if use_mean_wind:
            var = var + '_mean_wind'
            mean_wind_min_index = cube.data.mean(axis=(1, 2)).argmin()
            mean_wind = cube.data[mean_wind_min_index].mean()
        if mode:
            var = mode + '_' + var

        fig, ax = plt.subplots(dpi=100)
        data = cube.data[k]
        title = '{} xy slice at z={:.2f} km'.format(var, z[k])
        if anomaly:
            data = data - data.mean()
        if use_mean_wind:
            data = data - mean_wind
        if extent is None:
            extent = (0, 256, 0, 256)
        else:
            imin, imax, jmin, jmax = extent
            data = data[imin:imax, jmin:jmax]

        ax.set_title(title)
        # Coords are model_level, y, x or model_level, lat, lon
        kwargs = {}
        if use_norm:
            kwargs['norm'] = MidpointNormalize(midpoint=0, vmin=data.min(), vmax=data.max())
        if cmap:
            kwargs['cmap'] = cmap
        im = ax.imshow(data, origin='lower', extent=extent, **kwargs)

        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        plt.colorbar(im)
        if i:
            ax.vlines(x=i, ymin=0, ymax=256, color='k', linestyles='--')
        if j:
            ax.hlines(y=j, xmin=0, xmax=256, color='k', linestyles='--')
        self._create_dir_savefig('/{1}/xy/{2}/{0}_{1}_{2}_slice_{3}.png'
                                 .format(expt, self.task.runid, var, k))
        plt.close('all')

    def _plot_indiv_vert_slice(self, expt, var, plane, i, j, k, extent=None,
                               use_norm=False, cmap='bwr', anomaly=False, use_mean_wind=False,
                               vlev='theta', mode=''):
        cube = getattr(self, var)
        z = cube.coord('atmosphere_hybrid_height_coordinate').points / 1000
        # Coords are model_level, y, x or model_level, lat, lon
        if plane == 'xz':
            N = cube.shape[1]
            index = j
            vline_index = i
        elif plane == 'yz':
            N = cube.shape[2]
            index = i
            vline_index = j
        else:
            raise ValueError('plane=[xz|yz]')

        if anomaly:
            var = var + '_anom'
        if use_mean_wind:
            var = var + '_mean_wind'
            mean_wind_min_index = cube.data.mean(axis=(1, 2)).argmin()
            mean_wind = cube.data[mean_wind_min_index].mean()
        if mode:
            var = mode + '_' + var

        fig, ax = plt.subplots(dpi=100)
        ax.set_ylabel('height (km)')
        filename = ('/{2}/{0}/{3}/{1}_{2}_{3}_slice_{4}.png'
                    .format(plane, expt, self.task.runid, var, index))
        if plane == 'xz':
            ax.set_xlabel('x (km)')
            title = '{} xz slice at y={} gridbox'.format(var, index)
            data = cube.data[:, index]
        elif plane == 'yz':
            ax.set_xlabel('y (km)')
            title = '{} yz slice at x={} gridbox'.format(var, index)
            data = cube.data[:, :, index]
        if anomaly:
            data = data - data.mean(axis=1)[:, None]
        if use_mean_wind:
            data = data - mean_wind

        if extent is None:
            extent = (0, 256, 0, 20)
        else:
            imin, imax, _, _ = extent
            data = data[:, imin:imax]
            N = imax - imin

        ax.set_title(title)

        self._plot_vert_slice(ax, data, N, use_norm, cmap, vlev=vlev, extent=extent)
        if vline_index:
            ax.vlines(x=vline_index, ymin=0, ymax=20, color='k', linestyles='--')
        if k:
            ax.hlines(y=z[k], xmin=0, xmax=256, color='k', linestyles='--')
        self._create_dir_savefig(filename)
        plt.close('all')

    def _plot_slices(self, expt, var, **kwargs):
        cube = getattr(self, var)
        self._plot_xy_slices(expt, var, cube, **kwargs)
        self._plot_vert_slices(expt, var, cube, 'xz', **kwargs)
        self._plot_vert_slices(expt, var, cube, 'yz', **kwargs)

        self._plot_vert_mean_slice(expt, var, cube, 'xz', **kwargs)

    def _plot_xy_slices(self, expt, var, cube, **kwargs):
        for k in range(cube.shape[0]):
            self._plot_indiv_hor_slice(expt, var, None, None, k, **kwargs)

    def _plot_vert_slices(self, expt, var, cube, plane, **kwargs):
        # Coords are model_level, y, x or model_level, lat, lon
        if plane == 'xz':
            N = cube.shape[1]
        elif plane == 'yz':
            N = cube.shape[2]
        else:
            raise ValueError('plane=[xz|yz]')

        for i in range(N):
            if plane == 'xz':
                self._plot_indiv_vert_slice(expt, var, plane, None, i, None, **kwargs)
            elif plane == 'yz':
                self._plot_indiv_vert_slice(expt, var, plane, i, None, None, **kwargs)

    def _plot_vert_slice(self, ax, data, N, use_norm, cmap='bwr',
                         extent=(0, 256, 0, 20), aspect=10, vlev='theta'):
        if vlev == 'theta':
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_theta, np.arange(N),
                                                             data)
        elif vlev == 'rho':
            data_rbs = scipy.interpolate.RectBivariateSpline(self.vertlevs.z_rho, np.arange(N),
                                                             data)
        data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, N - 1, N))
        kwargs = {}
        if use_norm:
            kwargs['norm'] = MidpointNormalize(midpoint=0,
                                               vmin=data_interp[:200].min(),
                                               vmax=data_interp[:200].max())
        if cmap:
            kwargs['cmap'] = cmap

        im = ax.imshow(data_interp[:200], origin='lower', aspect=aspect, extent=extent,
                       **kwargs)

        plt.colorbar(im)

    def _plot_vert_mean_slice(self, expt, var, cube, plane='xy', use_norm=False, cmap='bwr',
                              anomaly=False, use_mean_wind=False, vlev='theta'):
        # Coords are model_level, y, x or model_level, lat, lon
        if plane == 'xz':
            cube_index = 1
        elif plane == 'yz':
            cube_index = 2
        else:
            raise ValueError('plane=[xz|yz]')
        N = cube.shape[cube_index]
        if anomaly:
            var = var + '_anom'
        if use_mean_wind:
            var = var + '_mean_wind'
            mean_wind_min_index = cube.data.mean(axis=(1, 2)).argmin()
            mean_wind = cube.data[mean_wind_min_index].mean()

        fig, ax = plt.subplots(dpi=100)
        ax.set_ylabel('height (km)')
        filename = ('/{2}/{0}/{3}/{1}_{2}_{3}_mean.png'
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
        if use_mean_wind:
            data = data - mean_wind
        ax.set_title(title)

        # Coords are model_level, y, x or model_level, lat, lon
        self._plot_vert_slice(ax, data, N, use_norm, cmap, vlev=vlev)
        self._create_dir_savefig(filename)
        plt.close('all')

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
