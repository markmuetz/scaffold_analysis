import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import iris

from omnium.analyzer import Analyzer
from omnium.utils import get_cube, count_blobs_mask
from omnium.consts import Re, L, cp, g


class ProfileAnalyzer(Analyzer):
    analysis_name = 'profile_analysis'

    def set_config(self, config):
	super(ProfileAnalyzer, self).set_config(config)
	self.w_thresh = config.getfloat('w_thresh', 1)
	self.qcl_thresh = config.getfloat('qcl_thresh', 0.0001)
        if 'mf_profile_xlim' in config:
            self.mf_profile_xlim = [float(v) for v in config['mf_profile_xlim'].split(',')]
        else:
            self.mf_profile_xlim = None
	    
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
        plt.savefig(self.figpath('_uv_profile.png'))

    def _plot_theta_qcl(self):
	plt.figure('theta')
        theta_profile = self.results['theta_profile']
        height = theta_profile.coord('level_height').points
	plt.title(self.expt)
        plt.plot(theta_profile.data, height, 'g-', label='theta')
	plt.xlim((250, 500))
	plt.ylim((0, 20000))
        plt.xlabel('theta (K)')
        plt.ylabel('height (m)')
        plt.savefig(self.figpath('_theta_profile.png'))

	plt.figure('qcl')
        qcl_profile = self.results['qcl_profile']
        height = qcl_profile.coord('level_height').points
	plt.title(self.expt)
        plt.plot(qcl_profile.data * 1000, height, 'g-', label='qcl')
	plt.ylim((0, 20000))
        plt.xlabel('qcl (g kg$^{-1}$)')
        plt.ylabel('height (m)')
        plt.savefig(self.figpath('_qcl_profile.png'))

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
        plt.savefig(self.figpath('_momentum_flux_profile.png'))

    def _plot_mass_flux(self):
        mf_cloud_profile = self.results['mf_cloud_profile']
        mf_w_profile = self.results['mf_w_profile']
        mf_qcl_profile = self.results['mf_qcl_profile']
        height = mf_cloud_profile.coord('level_height').points

        plt.clf()
	plt.title(self.expt)
        plt.plot(mf_cloud_profile.data, height, label='cloud mask')
        plt.plot(mf_w_profile.data, height, label='w mask')
        plt.plot(mf_qcl_profile.data, height, label='qcl mask')
	plt.ylim((0, 20000))
	if self.mf_profile_xlim:
	    plt.xlim(self.mf_profile_xlim)
        plt.legend()
        plt.savefig(self.figpath('_mass_flux_profile.png'))

    def run_analysis(self):
        cubes = self.cubes

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

        self.results['theta_profile'] = theta.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.results['qcl_profile'] = qcl.collapsed(['time', 'grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

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
        w = get_cube(cubes, 0, 150)
        rho = get_cube(cubes, 0, 253)

        rho.data = rho.data / Re**2
        rho.units = 'kg m-3'
        height = rho.coord('level_height').points

        # Interp theta grid vars onto rho grid.
        w_rho_grid = (w[:, :-1].data + w[:, 1:].data) / 2
        qcl_rho_grid = (qcl[:, :-1].data + qcl[:, 1:].data) / 2

        # Threshold.
        w_mask = w_rho_grid > self.w_thresh
        qcl_mask = qcl_rho_grid > self.qcl_thresh
        cloud_mask = w_mask & qcl_mask

        # Total mass-flux across domain.
        mf = w_rho_grid * rho.data
        # Calc. profiles.
        # N.B. mf[:, i][cloud_mask[:, i]].sum() is the conv. mass-flux at level i.
        #mf_cloud_profile_data = np.array([mf[:, i][cloud_mask[:, i]].sum() for i in range(mf.shape[1])])
        #mf_w_profile_data = np.array([mf[:, i][w_mask[:, i]].sum() for i in range(mf.shape[1])])
        #mf_qcl_profile_data = np.array([mf[:, i][qcl_mask[:, i]].sum() for i in range(mf.shape[1])])
        mf_cloud_profile_data = []
        mf_w_profile_data = []
        mf_qcl_profile_data = []

        for mask, profile_data in [(cloud_mask, mf_cloud_profile_data),
                                   (w_mask, mf_w_profile_data), 
                                   (qcl_mask, mf_qcl_profile_data)]:
            for i in range(mf.shape[1]):
                mf_conv = mf[:, i][mask[:, i]].sum()
                num_clouds = 0
                for itime in range(mf.shape[0]):
                    num_clouds += count_blobs_mask(mask[itime, i], diagonal=True)[0]
                profile_data.append(mf_conv/num_clouds * 4e6)

        mf_cloud_profile_data = np.array(mf_cloud_profile_data)
        mf_w_profile_data = np.array(mf_w_profile_data)
        mf_qcl_profile_data = np.array(mf_qcl_profile_data)

        mf_cloud_profile_data[np.isnan(mf_cloud_profile_data)] = 0
        mf_w_profile_data[np.isnan(mf_w_profile_data)] = 0
        mf_qcl_profile_data[np.isnan(mf_qcl_profile_data)] = 0

	height_coord = iris.coords.DimCoord(height, long_name='level_height', units='m')
        self.results['mf_cloud_profile'] = iris.cube.Cube(mf_cloud_profile_data, 
                                                          long_name='mass_flux_profile_cloud_mask',
                                                          dim_coords_and_dims=[(height_coord, 0)])
        self.results['mf_w_profile'] = iris.cube.Cube(mf_w_profile_data, 
                                                      long_name='mass_flux_profile_w_mask',
                                                      dim_coords_and_dims=[(height_coord, 0)])
        self.results['mf_qcl_profile'] = iris.cube.Cube(mf_qcl_profile_data, 
                                                        long_name='mass_flux_profile_qcl_mask',
                                                        dim_coords_and_dims=[(height_coord, 0)])

    def display_results(self):
        self._plot_uv()
        self._plot_theta_qcl()
        self._plot_momentum_flux()
        self._plot_mass_flux()
        plt.close('all')
