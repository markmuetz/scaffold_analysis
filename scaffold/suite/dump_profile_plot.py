from logging import getLogger

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt

has_metpy = False
try:
    import metpy
    import metpy.calc as mpcalc
    from metpy.plots import SkewT
    from metpy.units import units
    has_metpy = True
except ImportError:
    pass

from omnium import Analyser
from omnium.utils import get_cube
from omnium.consts import p_ref, kappa

from scaffold.utils import cm_to_inch, find_intersections
from scaffold.colour import EXPT_DETAILS

logger = getLogger('scaf.dump_prof_plot')
if not has_metpy:
    logger.warning('metpy not available')

# name, specific density, number conc., colour.
VARS = [
    ('qcl', (0, 254), (0, 75), 'b'),
    ('qcf', (0, 12), (0, 79), 'c'),
    ('qcf2', (0, 271), (0, 78), 'g'),
    ('qgr', (0, 273), (0, 81), 'r'),
    ('qrn', (0, 272), (0, 76), 'k'),
]


def plot_hydrometeor_profile(da, expt, ax1, ax2):
    for var in VARS:
        logger.debug('plotting profile of {} for {}', var, expt)
        varname = var[0]
        qvar = get_cube(da, *var[1])
        nvar = get_cube(da, *var[2])

        z = qvar.coord('atmosphere_hybrid_height_coordinate').points / 1000
        qvar_profile = qvar.data.mean(axis=(1, 2))
        nvar_profile = nvar.data.mean(axis=(1, 2))

        ax1.plot(qvar_profile * 1000, z, var[3], label=varname)
        ax2.plot(nvar_profile * 1000, z, var[3], label=varname)

        ax1.set_xscale('log')
        ax1.set_ylim((0, 20))
        ax1.set_xlim((10e-7, 10e-1))
        # ax1.set_title(expt)

        ax2.set_xscale('log')
        ax2.set_ylim((0, 20))
        ax2.set_xlim((10e-1, 10e8))

def plot_skewT(fig, name, p_profile, T_profile, Td_profile):
    skew = SkewT(fig, rotation=55)

    skew.plot(p_profile.to('hPa'), T_profile.to('degC'), 'r-')
    skew.plot(p_profile.to('hPa'), Td_profile.to('degC'), 'r--')

    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-30, 30)

    # skew.plot_dry_adiabats(t0=np.linspace(253.15, 303.15, 6) * units('K'))
    t0 = np.linspace(-30, 30, 7) * units('degC')
    skew.plot_dry_adiabats(t0=t0)
    skew.plot_moist_adiabats(t0=t0)
    skew.plot_mixing_lines()

    Tparcel_profile = mpcalc.parcel_profile(p_profile, T_profile[0], Td_profile[0]).to('degC')
    skew.plot(p_profile.to('hPa'), Tparcel_profile, 'k-')
    skew.shade_cape(p_profile, T_profile, Tparcel_profile)
    skew.shade_cin(p_profile, T_profile, Tparcel_profile)

    try:
        cape, cin = mpcalc.surface_based_cape_cin(p_profile, T_profile, Td_profile)
        skew.ax.set_title('{}\n' 
                          'CAPE = {:.2f} J kg$^{{-1}}$\n'
                          'CIN = {:.2f} J kg$^{{-1}}$'
                          .format(name, cape.magnitude, cin.magnitude))
    except Exception as e:
        logger.debug(e)
    skew.ax.set_xlabel('temperature ($^\circ$C)')
    skew.ax.set_ylabel('pressure (hPa)')

    fig.tight_layout()

class DumpProfilePlotter(Analyser):
    analysis_name = 'dump_profile_plot'
    multi_expt = True

    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da480.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.dump_profile_plot.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        # self._plot_hydrometeors()
        self._plot_skewT()
        self._plot_theta_profiles()
        plt.close('all')

    def _plot_hydrometeors(self):
        fig = plt.figure(dpi=100, figsize=(10, 10))
        axes = []
        num_expts = len(self.task.expts)
        axL = fig.add_subplot(num_expts, 2, 1)
        axR = fig.add_subplot(num_expts, 2, 2, sharey=axL)

        for i, expt in enumerate(self.task.expts):
            da = self.expt_cubes[expt]
            print(i)
            if i == 0:
                ax1, ax2 = axL, axR
            else:
                ax1 = fig.add_subplot(num_expts, 2, i * 2 + 1, sharey=axL, sharex=axL)
                ax2 = fig.add_subplot(num_expts, 2, i * 2 + 2, sharey=axL, sharex=axR)
            axes.append([ax1, ax2])

            plot_hydrometeor_profile(da, expt, ax1, ax2)

            indiv_fig = plt.figure(dpi=100)
            indiv_ax1 = indiv_fig.add_subplot(1, 2, 1)
            indiv_ax2 = indiv_fig.add_subplot(1, 2, 2, sharey=axL)
            plot_hydrometeor_profile(da, expt, indiv_ax1, indiv_ax2)
            indiv_ax1.set_xlabel('mass fraction (g kg$^{-1}$)')
            indiv_ax2.set_xlabel('number (# kg$^{-1}$)')
            indiv_ax2.legend(bbox_to_anchor=(0.8, 0.75))
            indiv_fig.savefig(self.file_path('hydrometeors_{}.png'.format(expt)))

        axR.legend(bbox_to_anchor=(0.8, 0.95))
        for i, expt in enumerate(self.task.expts):
            axes[i][0].set_ylabel('{}\nheight (km)'.format(expt))

            if i != num_expts - 1:
                plt.setp(axes[i][0].get_xticklabels(), visible=False)
                plt.setp(axes[i][1].get_xticklabels(), visible=False)
                plt.setp(axes[i][1].get_yticklabels(), visible=False)

        axes[-1][0].set_xlabel('mass fraction (g kg$^{-1}$)')
        axes[-1][1].set_xlabel('number (# kg$^{-1}$)')

        # plt.tight_layout()
        plt.savefig(self.file_path('hydrometeors.png'))
        plt.show()

    def _plot_skewT(self):
        for i, expt in enumerate(self.task.expts):
            da = self.expt_cubes[expt]

            exnerp = get_cube(da, 0, 255)
            theta = get_cube(da, 0, 4)

            qvdata = get_cube(da, 0, 10).data
            Tdata = theta.data * exnerp.data
            pdata = exnerp.data ** (1 / kappa) * p_ref

            p = pdata * units('Pa')
            qv = qvdata * units('kg/kg')
            T = Tdata * units('K')
            Td = mpcalc.dewpoint_from_specific_humidity(qv, T, p)

            fig = plt.figure(dpi=100, figsize=cm_to_inch(10, 12))

            plot_skewT(fig, expt,
                       p.mean(axis=(1, 2)),
                       T.mean(axis=(1, 2)),
                       Td.mean(axis=(1, 2)))
            plt.savefig(self.file_path('skewT_{}.png'.format(expt)))

            fig = plt.figure(dpi=100, figsize=cm_to_inch(10, 12))

            if expt in EXPT_DETAILS:
                ucp_kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
                plot_skewT(fig, ucp_kwargs['label'],
                           p.mean(axis=(1, 2)),
                           T.mean(axis=(1, 2)),
                           Td.mean(axis=(1, 2)))
                plt.savefig(self.file_path('UCP_skewT_{}.png'.format(expt)))

    def _plot_theta_profiles(self):
        for i, expt in enumerate(self.task.expts):
            logger.debug('plot thetas for {}', expt)
            da = self.expt_cubes[expt]

            theta_cube = get_cube(da, 0, 4)
            pi_cube = get_cube(da, 0, 255)
            q_cube = get_cube(da, 0, 10)
            q = q_cube.data

            theta = theta_cube.data
            p = pi_cube.data**(1 / kappa) * (p_ref / 100)
            T = theta * pi_cube.data

            z = theta_cube.coord('atmosphere_hybrid_height_coordinate').points

            Td = mpcalc.dewpoint_from_specific_humidity(q * units('kg/kg'), T * units.degK, p * units('hPa'))
            theta_e = mpcalc.equivalent_potential_temperature(p * units('hPa'), T * units.degK, Td)
            theta_es = mpcalc.saturation_equivalent_potential_temperature(p * units('hPa'), T * units.degK)
            title = 'thetas_' + expt

            plt.figure(title)
            plt.clf()
            plt.title(title)

            theta_profile = theta.mean(axis=(1, 2))
            theta_e_profile = theta_e.mean(axis=(1, 2))
            theta_es_profile = theta_es.mean(axis=(1, 2))
            p_profile = p.mean(axis=(1, 2))

            plt.plot(theta_profile, z, 'r-', label='$\\theta$')
            plt.plot(theta_e_profile, z, 'g-', label='$\\theta_{e}$')
            plt.plot(theta_es_profile, z, 'b-', label='$\\theta_{es}$')
            plt.axvline(x=theta_e_profile[0], linestyle='--', color='k', label='parcel ascent')

            indices, weights = find_intersections(theta_es_profile.magnitude,
                                                  np.ones_like(theta_e_profile) * theta_e_profile[0].magnitude)

            i, w = indices[0], weights[0]
            z_lfc = z[i] + w * (z[i + 1] - z[i])
            p_lfc = p_profile[i] + w * (p_profile[i + 1] - p_profile[i])
            plt.plot((theta_e_profile.magnitude[0] - 2, theta_e_profile.magnitude[0] + 2),
                     (z_lfc, z_lfc),
                     linestyle='--', color='grey', label='LFC ({:.0f} m)'.format(z_lfc))

            i, w = indices[1], weights[1]
            z_lnb = z[i] + w * (z[i + 1] - z[i])
            p_lnb = p_profile[i] + w * (p_profile[i + 1] - p_profile[i])
            plt.plot((theta_e_profile.magnitude[0] - 2, theta_e_profile.magnitude[0] + 2),
                     (z_lnb, z_lnb),
                     linestyle='--', color='brown', label='LNB ({:.0f} m)'.format(z_lnb))


            logger.debug('LFC: {:.0f} m, {:.0f} hPa'.format(z_lfc, p_lfc))
            logger.debug('LNB: {:.0f} m, {:.0f} hPa'.format(z_lnb, p_lnb))

            plt.legend()

            plt.ylabel('height (m)')
            plt.xlim((290, 360))
            plt.ylim((0, 15000))
            plt.savefig(self.file_path(title))
