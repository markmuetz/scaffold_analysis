import os
from logging import getLogger

import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib import rcParams
import pylab as plt

from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr

from scaffold.vertlev import VertLev
from scaffold.utils import cm_to_inch

logger = getLogger('scaf.prof_plot')


class ProfilePlotter(Analyser):
    analysis_name = 'profile_plot'
    multi_expt = True
    base_u_profile = np.array([(0, -2), (1e3, -3), (12e3, 2.5), (14.5e3, 0), (40e3, 0)])

    def run_analysis(self):
        pass

    def _plot_poster_shear_profiles(self):
        f, ax1 = plt.subplots(1, 1, sharey=True)
        f.set_size_inches(*cm_to_inch(8, 12))
        vertlev = VertLev(self.suite.suite_dir)
        ax1.plot(vertlev.dz_theta, vertlev.z_rho / 1e3)
        ax1.set_xlabel('$\\Delta z$ (m)')
        ax1.set_ylabel('height (km)')

        ax1.set_ylim((0, 25))

        for expt in self.expts:
            cubes = self.expt_cubes[expt]
            u_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'u_profile')
            # v_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'v_profile')
            height = u_profile.coord('level_height').points
            shear_factor = int(expt[1]) # i.e. 0-5.
            shear_u_profile = self.base_u_profile.copy()
            shear_u_profile[:, 1] *= shear_factor
            # N.B. convert m->km.
            plot = ax1.plot(shear_u_profile[:, 1], shear_u_profile[:, 0] / 1e3, label=expt)
            colour = plot[0].get_color()
            #import ipdb; ipdb.set_trace()
            ax1.plot(u_profile.data, height / 1e3, color=colour, linestyle='--')

        ax1.set_xlim((-15, 15))
        #ax1.set_ylim((0, 20))
        ax1.set_xlabel('u profile (m s$^{-1}$)')
        ax1.legend(bbox_to_anchor=(1.05, 1.05))

        #plt.tight_layout()
        #plt.subplots_adjust(wspace=0)
        plt.savefig(self.figpath('poster_shear_profiles.png'))

    def _plot_input_profiles(self):
        f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
        vertlev = VertLev(self.suite.suite_dir)
        ax1.plot(vertlev.dz_theta, vertlev.z_rho / 1e3)
        ax1.set_xlabel('$\\Delta z$ (m)')
        ax1.set_ylabel('height (km)')
        ax1.set_ylim((0, 25))

        # Cooling profile is very similar in height coords for different expts.
        # Just use first.
        expt = self.expts[0]
        cubes = self.expt_cubes[expt]
        # in Pa.
        pressure_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'pressure_profile')
        z = pressure_profile.coord('level_height').points
        pressure_data = pressure_profile.data / 100  # Convert to hPa.

        cooling = np.zeros_like(pressure_profile.data)

        cooling[pressure_data > 200] = -2
        lin_region = (pressure_data < 200) & (pressure_data > 100)
        cooling[lin_region] = -(1 - (200 - pressure_data[lin_region])/(200 - 100)) * 2

        ax2.plot(cooling, z / 1e3, 'b-')
        ax2.set_xlim((-3, 3))
        #plt.ylabel('height (km)')
        ax2.set_xlabel('prescribed heating (K day$^{-1}$)')
        ax2.axvline(x=0, color='k', linestyle='--')

        for expt in self.expts:
            cubes = self.expt_cubes[expt]
            u_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'u_profile')
            # v_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'v_profile')
            height = u_profile.coord('level_height').points
            shear_factor = int(expt[1]) # i.e. 0-5.
            shear_u_profile = self.base_u_profile.copy()
            shear_u_profile[:, 1] *= shear_factor
            # N.B. convert m->km.
            plot = ax3.plot(shear_u_profile[:, 1], shear_u_profile[:, 0] / 1e3, label=expt)
            colour = plot[0].get_color()
            #import ipdb; ipdb.set_trace()
            ax3.plot(u_profile.data, height / 1e3, color=colour, linestyle='--')

        ax3.set_xlim((-15, 15))
        #ax3.set_ylim((0, 20))
        ax3.set_xlabel('u profile (m s$^{-1}$)')
        ax3.legend(loc='upper left')

        #plt.tight_layout()
        #plt.subplots_adjust(wspace=0)
        plt.savefig(self.figpath('input_profiles.png'))

    def _plot_uv_profile(self):
        fig = plt.figure('uv_profile', figsize=(3.5, 4.5), dpi=1200)
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
        plt.legend(loc='upper left')
        plt.tight_layout()
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
            # Plot a marker at each end of the cloud profile.
            ax1.plot(theta_cloud_profile.data[cloud_profile_mask][0],
                     height[cloud_profile_mask][0],
                     color=colour, marker='+')
            ax1.plot(theta_cloud_profile.data[cloud_profile_mask][-1],
                     height[cloud_profile_mask][-1],
                     color=colour, marker='o')
            # This is the TD profile outside of diagnosed clouds.
            # It almost exactly overlies the domain mean.
            # plt.plot(theta_not_cloud_profile.data, height, color=colour, marker='+')

        ax1.set_xlim((290, 360))
        ax1.set_ylim((0, 18))
        ax1.set_xlabel('$\\theta$ (K)')
        ax1.set_ylabel('height (km)')
        ax1.legend(loc='upper left')

        for expt in self.expts:
            cubes = self.expt_cubes[expt]

            qcl_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcl_profile')
            qgr_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qgr_profile')
            qcf_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcf_profile')

            # qcl_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcl_cloud_profile')
            # qcl_not_cloud_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'qcl_not_cloud_profile')

            height = qcl_profile.coord('level_height').points / 1e3  # Convert m->km.

            # cloud_profile_mask = qcl_cloud_profile.data > 0.0000001
            # kg kg-1 -> g kg-1
            plot = ax2.plot(qcl_profile.data * 1e3, height, label=expt)
            colour = plot[0].get_color()
            ax2.plot(qgr_profile.data * 1e3, height, color=colour, linestyle=':')
            ax2.plot(qcf_profile.data * 1e3, height, color=colour, linestyle='-.')
            ax2.set_xlim((1e-3, 1e-1))
            ax2.set_xscale('log')
            #ax2.plot(qcl_cloud_profile.data[cloud_profile_mask] * 1e3, height[cloud_profile_mask],
            #         color=colour, linestyle='--')
            #ax2.plot(qcl_cloud_profile.data[cloud_profile_mask][-1],
            #         height[cloud_profile_mask][-1],
            #         color=colour, marker='o')
            #ax2.plot(theta_not_cloud_profile.data, height, color=colour, marker='+')

        ax2.set_xlabel('hydrometeors (g kg$^{-1}$)')
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

            plot = ax1.plot(clouds_cloud_profile.data, height)
            colour = plot[0].get_color()
            ax2.plot(mf_cloud_profile.data / 1e8, height, color=colour)
            ax3.plot(clouds_cloud_profile.data * mf_cloud_profile.data / 1e8, height,
                     color=colour, label=expt)

        ax1.set_ylim((0, 18))
        ax1.set_ylabel('height (km)')

        ax1.set_xlabel('Number of clouds')
        ax2.set_xlabel('MF per cloud ($\\times 10^8$ kg s$^{-1}$)')
        ax3.set_xlabel('Total MF ($\\times 10^8$ kg s$^{-1}$)')
        ax3.set_xlim((0, 100))

        plt.legend(loc='upper right')
        plt.savefig(self.figpath('mf_profile.png'))

    def _plot_dz_profile(self):
        try:
            import f90nml
        except:
            logger.warn('f90nml not installed')
            return

        plt.figure('dz_profile')
        vertlevs_filename = os.path.join(self.suite.suite_dir, 'app/um/file/rce_vertlevs.nml')

        vertlevs = f90nml.read(vertlevs_filename)['vertlevs']
        eta_theta = np.array(vertlevs['eta_theta'])
        eta_rho = np.array(vertlevs['eta_rho'])
        z_top = vertlevs['z_top_of_model']
        assert len(eta_theta) == len(eta_rho) + 1

        z_theta = eta_theta * z_top
        z_rho = eta_rho * z_top
        dz_theta = z_theta[1:] - z_theta[:-1]

        plt.plot(dz_theta, z_rho)
        plt.xlabel('$\\Delta z$ theta (m)')
        plt.ylabel('height (m)')
        plt.savefig(self.figpath('dz_profile.png'))

    def _plot_momentum_flux(self):
        fig = plt.figure('momf_profile', figsize=(3.5, 4.5), dpi=1200)
        for expt in self.expts:
            cubes = self.expt_cubes[expt]
            u_mom_flux_ts = get_cube_from_attr(cubes, 'omnium_cube_id', 'u_mom_flux_ts')
            v_mom_flux_ts = get_cube_from_attr(cubes, 'omnium_cube_id', 'v_mom_flux_ts')
            z = u_mom_flux_ts.coord('level_height').points

            plt.plot(u_mom_flux_ts.data.mean(axis=0) * 1e3, z / 1e3, label=expt)
            #plt.plot(v_mom_flux_ts.data.mean(axis=0), z, 'b--', label='v')
            plt.ylabel('height (km)')
            plt.xlabel('mom flux ($\\times 10^{-3}$ kg m$^{-1}$ s$^{-2}$)')

        plt.ylim((0, 20))
        plt.legend(loc='upper left')
        plt.savefig(self.figpath('momentum_flux_profile.png'))

    def _plot_cooling(self):
        for expt in self.expts:
            cubes = self.expt_cubes[expt]
            # in Pa.
            pressure_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'pressure_profile')
            z = pressure_profile.coord('level_height').points
            pressure_data = pressure_profile.data / 100  # Convert to hPa.

            plt.figure(self.output_filename + '_pressure_profile')
            plt.clf()

            cooling = np.zeros_like(pressure_profile.data)

            cooling[pressure_data > 200] = -2
            lin_region = (pressure_data < 200) & (pressure_data > 100)
            cooling[lin_region] = -(1 - (200 - pressure_data[lin_region])/(200 - 100)) * 2

            plt.plot(cooling, z / 1e3, 'b-')
            plt.xlim((-3, 3))
            plt.ylabel('height (km)')
            plt.xlabel('prescribed heating (K day$^{-1}$)')

            plt.savefig(self.figpath('{}.pressure_profile.png'.format(expt)))

    def display_results(self):
        rcParams.update({'figure.autolayout': True})
        self._plot_input_profiles()
        self._plot_uv_profile()
        self._plot_thermodynamic_profile()
        self._plot_mf_profile()
        self._plot_dz_profile()
        self._plot_momentum_flux()
        self._plot_cooling()
        self._plot_poster_shear_profiles()
        rcParams.update({'figure.autolayout': False})
        plt.close('all')
