from logging import getLogger

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt
import iris

from omnium import Analyser, ExptList
from omnium.utils import get_cube
from omnium.consts import Re, cp, g
from cloud_tracking.utils import label_clds

logger = getLogger('scaf.prof_an')


def mean(a, axis=None):
    """Workaround for old versions of numpy, where using axis=tuple on masked_arrays doesn't work.

    https://stackoverflow.com/questions/30209624/numpy-mean-used-with-a-tuple-as-axis-argument-not-working-with-a-masked-arr
    """
    sums = a.sum(axis=axis)
    counts = np.logical_not(a.mask).sum(axis=axis)
    result = sums * 1. / counts
    return result


class ProfileAnalyser(Analyser):
    """Calculates profiles for u, v, w, theta, qcl (+qi), in-cloud theta, qcl, mass-flux and # cld.

    Calculates energy lost using pressure from model.
    Calculates *average* profile for all times in the loaded cubes.
    """
    analysis_name = 'profile_analysis'
    single_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.pp2.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.{runid:03}.profile_analysis.nc']

    uses_runid = True
    runid_pattern = 'atmos.(?P<runid>\d{3}).pp2.nc'

    def load(self):
        self.load_cubes()

    def run(self):
        cubes = self.cubes
        logger.debug('cubes: {}'.format(cubes))

        expts = ExptList(self.suite)
        expts.find([self.task.expt])
        expt_obj = expts.get(self.task.expt)
        dx, dy = expt_obj.dx, expt_obj.dy

        # u/v profile.
        u = get_cube(cubes, 0, 2)
        v = get_cube(cubes, 0, 3)
        u_profile = u.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        v_profile = v.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.u_heights = u.coord('level_height').points

        self.results['u_profile'] = u_profile
        self.results['v_profile'] = v_profile

        # momentum flux profile.
        u_inc = get_cube(cubes, 53, 185)
        v_inc = get_cube(cubes, 53, 186)
        rho = (get_cube(cubes, 0, 253) / Re ** 2)

        u_inc_ts = u_inc.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        v_inc_ts = v_inc.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        rho_ts = rho.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

        theta = get_cube(cubes, 0, 4)
        qcl = get_cube(cubes, 0, 254)
        qgr = get_cube(cubes, 0, 273)
        qcf = get_cube(cubes, 0, 12)

        w = get_cube(cubes, 0, 150)
        rho = get_cube(cubes, 0, 253)
        pressure = get_cube(cubes, 0, 408)

        logger.debug('got cubes of interest')

        w_mask = w.data > self.settings.w_thresh
        qcl_mask = qcl.data > self.settings.qcl_thresh
        cloud_mask = w_mask & qcl_mask

        self.results['theta_profile'] = theta.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.results['qcl_profile'] = qcl.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.results['qgr_profile'] = qgr.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.results['qcf_profile'] = qcf.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.results['pressure_profile'] = pressure.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

        logger.debug('performed avging')

        # Make masked arrays to grab subsets of data I want.
        # N.B. masked where data will be ignored - hence ~cloud_mask.
        theta_cloud_data = np.ma.array(theta.data, mask=~cloud_mask)
        theta_not_cloud_data = np.ma.array(theta.data, mask=cloud_mask)
        qcl_cloud_data = np.ma.array(qcl.data, mask=~cloud_mask)
        qcl_not_cloud_data = np.ma.array(qcl.data, mask=cloud_mask)

        logger.debug('masked data')

        # Work out theta profiles in/out of clouds.
        theta_cloud_profile = self.results['theta_profile'].copy()
        theta_not_cloud_profile = self.results['theta_profile'].copy()
        # This line caused problems under numpy 1.9.2 - due to axis being tuple.
        theta_cloud_profile.data = mean(theta_cloud_data, axis=(0, 2, 3)).data
        theta_not_cloud_profile.data = mean(theta_not_cloud_data, axis=(0, 2, 3)).data
        self.results['theta_cloud_profile'] = theta_cloud_profile
        self.results['theta_not_cloud_profile'] = theta_not_cloud_profile

        qcl_cloud_profile = self.results['qcl_profile'].copy()
        qcl_not_cloud_profile = self.results['qcl_profile'].copy()
        qcl_cloud_profile.data = mean(qcl_cloud_data, axis=(0, 2, 3)).data
        qcl_not_cloud_profile.data = mean(qcl_not_cloud_data, axis=(0, 2, 3)).data
        self.results['qcl_cloud_profile'] = qcl_cloud_profile
        self.results['qcl_not_cloud_profile'] = qcl_not_cloud_profile

        logger.debug('created masked arrays')

        z = theta.coord('level_height').points
        dz = z[1:] - z[:-1]

        dz_ts = dz.repeat(u_inc_ts.shape[0]).reshape(*u_inc_ts.shape)
        delta_t = 30 # 30s.
        u_mom_flux_ts = rho_ts * u_inc_ts * dz_ts / delta_t
        v_mom_flux_ts = rho_ts * v_inc_ts * dz_ts / delta_t

        start_time = u_inc_ts.coord('time').points[0]
        self.times = u_inc_ts.coord('time').points - start_time

        self.results['u_mom_flux_ts'] = u_mom_flux_ts
        self.results['v_mom_flux_ts'] = v_mom_flux_ts

        # mass flux profile.
        rho.data = rho.data / Re**2
        rho.units = 'kg m-3'
        height = rho.coord('level_height').points

        # Interp theta grid vars onto rho grid.
        w_rho_grid = (w[:, :-1].data + w[:, 1:].data) / 2
        qcl_rho_grid = (qcl[:, :-1].data + qcl[:, 1:].data) / 2

        # Threshold.
        w_rho_mask = w_rho_grid > self.settings.w_thresh
        qcl_rho_mask = qcl_rho_grid > self.settings.qcl_thresh
        cloud_rho_mask = w_rho_mask & qcl_rho_mask

        # Total mass-flux across domain.
        mf = w_rho_grid * rho.data
        # Calc. profiles.
        # N.B. mf[:, i][cloud_rho_mask[:, i]].sum() is the conv. mass-flux at level i.
        #mf_cloud_profile_data = np.array([mf[:, i][cloud_rho_mask[:, i]].sum() for i in range(mf.shape[1])])
        #mf_w_profile_data = np.array([mf[:, i][w_mask[:, i]].sum() for i in range(mf.shape[1])])
        #mf_qcl_profile_data = np.array([mf[:, i][qcl_rho_mask[:, i]].sum() for i in range(mf.shape[1])])

        height_coord = iris.coords.DimCoord(height, long_name='level_height', units='m')
        for name, mask in [('cloud_profile', cloud_rho_mask),
                           ('w_profile', w_rho_mask),
                           ('qcl_profile', qcl_rho_mask)]:
            logger.debug('running {}'.format(name))
            profile_data = []
            cloud_data = []

            # Loop over height.
            for i in range(mf.shape[1]):
                mf_conv = mf[:, i][mask[:, i]].sum()
                num_clouds = 0
                for itime in range(mf.shape[0]):
                    num_clouds += label_clds(mask[itime, i], diagonal=True)[0]
                profile_data.append(mf_conv/num_clouds * dx * dy)
                cloud_data.append(num_clouds)

            profile_data = np.array(profile_data)
            profile_data[np.isnan(profile_data)] = 0
            mf_name = 'mf_' + name
            cloud_name = 'clouds_' + name
            self.results[mf_name] = iris.cube.Cube(profile_data,
                                                   long_name=mf_name,
                                                   dim_coords_and_dims=[(height_coord, 0)])
            self.results[cloud_name] = iris.cube.Cube(cloud_data,
                                                      long_name=cloud_name,
                                                      dim_coords_and_dims=[(height_coord, 0)])

        self.calc_energy_loss_rate()

    def save(self, state, suite):
        self.save_results_cubes(state, suite)

    def display_results(self):
        self.mfpercloud_profile_xlim = None
        self.mftotal_profile_xlim = None
        self.cloud_profile_xlim = None
        self._plot_uv()
        self._plot_theta_qcl()
        self._plot_momentum_flux()
        self._plot_mass_flux()
        plt.close('all')

    def _plot_uv(self):
        u_profile = self.results['u_profile']
        v_profile = self.results['v_profile']
        height = u_profile.coord('level_height').points
        plt.figure(self.output_filename + '_uv_profile')
        plt.clf()
        plt.title(self.output_filename + '_uv_profile')

        plt.plot(u_profile.data, height, 'g-', label='u')
        plt.plot(v_profile.data, height, 'b--', label='v')

        plt.ylabel('height (m)')
        plt.xlabel('(m s$^{-1}$)')

        plt.legend()
        plt.savefig(self.file_path('uv_profile.png'))

    def _plot_theta_qcl(self):
        plt.figure('theta')
        theta_profile = self.results['theta_profile']
        theta_cloud_profile = self.results['theta_cloud_profile']
        theta_not_cloud_profile = self.results['theta_not_cloud_profile']
        height = theta_profile.coord('level_height').points
        plt.title(self.task.expt)
        plt.plot(theta_profile.data, height, 'g-', label='theta')
        plt.plot(theta_cloud_profile.data, height, 'g--', label='theta cloud')
        plt.plot(theta_not_cloud_profile.data, height, 'g.', label='theta not_cloud')
        plt.xlim((250, 500))
        plt.ylim((0, 20000))
        plt.xlabel('theta (K)')
        plt.ylabel('height (m)')
        plt.legend()
        plt.savefig(self.file_path('theta_profile.png'))

        plt.figure('qcl')
        qcl_profile = self.results['qcl_profile']
        qcl_cloud_profile = self.results['qcl_cloud_profile']
        qcl_not_cloud_profile = self.results['qcl_not_cloud_profile']
        height = qcl_profile.coord('level_height').points
        plt.title(self.task.expt)
        plt.plot(qcl_profile.data * 1000, height, 'g-', label='qcl')
        plt.plot(qcl_cloud_profile.data * 1000, height, 'g--', label='qcl cloud')
        plt.plot(qcl_not_cloud_profile.data * 1000, height, 'g.', label='qcl not_cloud')
        plt.ylim((0, 20000))
        plt.xlabel('qcl (g kg$^{-1}$)')
        plt.ylabel('height (m)')
        plt.legend()
        plt.savefig(self.file_path('qcl_profile.png'))

        #plt.figure('qgr')
        qgr_profile = self.results['qgr_profile']
        plt.title(self.task.expt)
        plt.plot(qgr_profile.data * 1000, height, 'r-', label='qgr')
        plt.ylim((0, 20000))
        plt.xlabel('qgr (g kg$^{-1}$)')
        plt.ylabel('height (m)')
        plt.legend()
        #plt.savefig(self.file_path('qgr_profile.png'))

        #plt.figure('qcf')
        qcf_profile = self.results['qcf_profile']
        plt.title(self.task.expt)
        plt.plot(qcf_profile.data * 1000, height, 'b-', label='qcf')
        plt.xlim((0, 0.1))
        plt.ylim((0, 20000))
        plt.xlabel('qcf (g kg$^{-1}$)')
        plt.ylabel('height (m)')
        plt.legend()
        plt.savefig(self.file_path('hydrom_profile.png'))

    def _plot_momentum_flux(self):
        u_mom_flux_ts = self.results['u_mom_flux_ts']
        v_mom_flux_ts = self.results['v_mom_flux_ts']
        z = u_mom_flux_ts.coord('level_height').points

        plt.figure(self.output_filename + '_momentum_flux_profile')
        plt.clf()
        plt.title(self.output_filename + '_momentum_flux_profile')

        plt.plot(u_mom_flux_ts.data.mean(axis=0), z, 'g-', label='u')
        plt.plot(v_mom_flux_ts.data.mean(axis=0), z, 'b--', label='v')
        plt.ylabel('height (m)')
        plt.xlabel('mom flux (kg m$^{-1}$ s$^{-2}$)')
        plt.legend()
        plt.savefig(self.file_path('momentum_flux_profile.png'))

    def _plot_mass_flux(self):
        mf_cloud_profile = self.results['mf_cloud_profile']
        mf_w_profile = self.results['mf_w_profile']
        mf_qcl_profile = self.results['mf_qcl_profile']
        height = mf_cloud_profile.coord('level_height').points

        clouds_cloud_profile = self.results['clouds_cloud_profile']
        clouds_w_profile = self.results['clouds_w_profile']
        clouds_qcl_profile = self.results['clouds_qcl_profile']

        plt.clf()
        plt.title(self.task.expt + ': mf/cloud')
        plt.plot(mf_cloud_profile.data, height, label='cloud mask')
        plt.plot(mf_w_profile.data, height, label='w mask')
        plt.plot(mf_qcl_profile.data, height, label='qcl mask')
        plt.ylim((0, 20000))
        if self.mfpercloud_profile_xlim:
            plt.xlim(self.mfpercloud_profile_xlim)
        plt.legend()
        plt.savefig(self.file_path('mfpercloud_profile.png'))

        plt.clf()
        plt.title(self.task.expt + ': #clouds')
        plt.plot(clouds_cloud_profile.data, height, label='cloud mask')
        plt.plot(clouds_w_profile.data, height, label='w mask')
        plt.plot(clouds_qcl_profile.data, height, label='qcl mask')
        plt.ylim((0, 20000))
        if self.cloud_profile_xlim:
            plt.xlim(self.cloud_profile_xlim)
        plt.legend()
        plt.savefig(self.file_path('numclouds_profile.png'))

        plt.clf()
        plt.title(self.task.expt + ': total mf')
        plt.plot(clouds_cloud_profile.data * mf_cloud_profile.data, height, label='cloud mask')
        plt.plot(clouds_w_profile.data * mf_w_profile.data, height, label='w mask')
        plt.plot(clouds_qcl_profile.data * mf_qcl_profile.data, height, label='qcl mask')
        plt.ylim((0, 20000))
        if self.mftotal_profile_xlim:
            plt.xlim(self.mftotal_profile_xlim)
        plt.legend()
        plt.savefig(self.file_path('totalmf_profile.png'))

    def calc_energy_loss_rate(self):
        # Integ(-cp / g * C(z), p_0, p_TOA, dp)
        # 1e2: convert hPa to Pa.
        p_profile = self.results['pressure_profile']  # Pa
        pdata = p_profile.data

        p0 = pdata[0]
        # N.B. these are all -ve.
        dp = (pdata[2:] - pdata[:-2]) / 2

        # Calc analytically (C(z) const)
        # loss, therefore want +ve.
        C_max = 2. / 86400  # [K/s]
        E_loss_up_to_200hPa = cp / g * C_max * (p0 - 200 * 1e2)

        # N.B slice means pdata[1:-1] has same length as dp.
        # lin_region can be used on pdata[1:-1] or dp.
        lin_region = (pdata[1:-1] < 200 * 1e2) & (pdata[1:-1] > 100 * 1e2)
        # -ve.
        cooling = -(1 - (200 * 1e2 - pdata[1:-1][lin_region])/(200 * 1e2 - 100 * 1e2)) * C_max

        # +ve.
        E_loss_lin_region = (cp / g * cooling * dp[lin_region]).sum()
        E_loss = E_loss_up_to_200hPa + E_loss_lin_region
        self.save_text('energy_loss.txt', 'Energy loss rate [W m-2]: {}\n'.format(E_loss))
        p_profile.attributes['Energy loss rate [W m-2]'] = E_loss
