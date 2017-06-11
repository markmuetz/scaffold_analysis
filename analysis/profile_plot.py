import os
from collections import OrderedDict
from itertools import groupby

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class ProfilePlotter(Analyzer):
    analysis_name = 'profile_plot'
    multi_expt = True
    base_u_profile = np.array([(0, -2), (1e3, -3), (12e3, 2.5), (14.5e3, 0), (40e3, 0)])

    def run_analysis(self):
        pass

    def _plot_uv_profile(self):
        plt.figure('uv_profile')
	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            u_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'u_profile')
            # v_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'v_profile')
            height = u_profile.coord('level_height').points
            shear_factor = int(expt[1]) # i.e. 0-5.
            shear_u_profile = self.base_u_profile.copy()
            shear_u_profile[:, 1] *= shear_factor
	    # N.B. convert m->km.
            plot = plt.plot(shear_u_profile[:, 1], shear_u_profile[:, 0] / 1e3, label=expt)
            colour = plot[0].get_color()
            plt.plot(u_profile.data, height / 1e3, color=colour, linestyle='--')
        plt.xlim((-15, 15))
        plt.ylim((0, 20))
        plt.ylabel('height (km)')
        plt.xlabel('u profile (m s$^{-1}$)')
        plt.legend()
        plt.savefig(self.figpath('uv_profile.png'))

    def _plot_thermodynamic_profile(self):
        plt.figure('thermodynamic_profile')
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            theta_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta_profile')
            theta_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta_cloud_profile')
            theta_not_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta_not_cloud_profile')
            height = theta_profile.coord('level_height').points / 1e3  # Convert m->km.

            cloud_profile_mask = theta_cloud_profile.data > 1
            plot = ax1.plot(theta_profile.data, height, label=expt)
            colour = plot[0].get_color()
            ax1.plot(theta_cloud_profile.data[cloud_profile_mask], height[cloud_profile_mask], 
                     color=colour, linestyle='--')
            ax1.plot(theta_cloud_profile.data[cloud_profile_mask][-1],
                     height[cloud_profile_mask][-1], 
                     color=colour, marker='o')
            # This is the TD profile outside of diagnosed clouds.
            # It almost exactly overlies the domain mean.
            # plt.plot(theta_not_cloud_profile.data, height, color=colour, marker='+')

        ax1.set_xlim((290, 360))
        ax1.set_ylim((0, 18))
        ax1.set_xlabel('theta (K)')
        ax1.set_ylabel('height (km)')

	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            qcl_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcl_profile')
            qcl_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcl_cloud_profile')
            qcl_not_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcl_not_cloud_profile')

            height = qcl_profile.coord('level_height').points / 1e3  # Convert m->km.

            cloud_profile_mask = qcl_cloud_profile.data > 0.0000001
            # kg kg-1 -> g kg-1
            plot = ax2.plot(qcl_profile.data * 1e3, height, label=expt)
            colour = plot[0].get_color()
            #ax2.plot(qcl_cloud_profile.data[cloud_profile_mask] * 1e3, height[cloud_profile_mask], 
            #         color=colour, linestyle='--')
            #ax2.plot(qcl_cloud_profile.data[cloud_profile_mask][-1],
            #         height[cloud_profile_mask][-1], 
            #         color=colour, marker='o')
            #ax2.plot(theta_not_cloud_profile.data, height, color=colour, marker='+')

        ax2.set_xlabel('cloud liquid water (g kg$^{-1}$')
	plt.legend(loc='upper right')
        plt.savefig(self.figpath('thermodynamic_profile.png'))

    def _plot_mf_profile(self):
        plt.figure('mf_profile')
        f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
	for expt in self.expts:
	    cubes = self.expt_cubes[expt]

            mf_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mf_cloud_profile')
            mf_w_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mf_w_profile')
            mf_qcl_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mf_qcl_profile')
            height = mf_cloud_profile.coord('level_height').points / 1e3  # m -> km

            clouds_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'clouds_cloud_profile')
            clouds_w_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'clouds_w_profile')
            clouds_qcl_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'clouds_qcl_profile')

            plot = ax1.plot(mf_cloud_profile.data / 1e7, height)
            colour = plot[0].get_color()
            ax2.plot(clouds_cloud_profile.data, height, color=colour)
            ax3.plot(clouds_cloud_profile.data * mf_cloud_profile.data / 1e8, height, 
                     color=colour, label=expt)

        ax1.set_ylim((0, 18))
        ax1.set_ylabel('height (km)')

        ax1.set_xlabel('MF per cloud ($\\times 10^7$ kg s$^{-1}$)')
        ax2.set_xlabel('Number of clouds')
        ax3.set_xlabel('Total MF ($\\times 10^8$ kg s$^{-1}$)')
        ax3.set_xlim((0, 100))

	plt.legend(loc='upper right')
        plt.savefig(self.figpath('mf_profile.png'))

    def display_results(self):
        self._plot_uv_profile()
        self._plot_thermodynamic_profile()
        self._plot_mf_profile()
        plt.close('all')
