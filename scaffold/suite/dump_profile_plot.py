from logging import getLogger

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt

import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units

from omnium import Analyser
from omnium.utils import get_cube
from omnium.consts import p_ref, kappa

logger = getLogger('scaf.dump_prof_plot')

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
    skew = SkewT(fig, rotation=45)

    skew.plot(p_profile.to('hPa'), T_profile.to('degC'), 'r-')
    skew.plot(p_profile.to('hPa'), Td_profile.to('degC'), 'r--')

    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-30, 40)

    Tparcel_profile = mpcalc.parcel_profile(p_profile, T_profile[0], Td_profile[0]).to('degC')
    skew.plot(p_profile.to('hPa'), Tparcel_profile, 'k-')
    skew.shade_cape(p_profile, T_profile, Tparcel_profile)
    skew.shade_cin(p_profile, T_profile, Tparcel_profile)

    try:
        cape, cin = mpcalc.surface_based_cape_cin(p_profile, T_profile, Td_profile)
        skew.ax.set_title('CAPE = {:.2f} J kg$^{{-1}}$\n'
                          'CIN = {:.2f} J kg$^{{-1}}$'
                          .format(cape.magnitude, cin.magnitude))
    except:
        pass

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
        self._plot_hydrometeors()
        self._plot_skewT()
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

            fig = plt.figure(dpi=100, figsize=(10, 10))

            plot_skewT(fig, expt,
                       p.mean(axis=(1, 2)),
                       T.mean(axis=(1, 2)),
                       Td.mean(axis=(1, 2)))
            plt.savefig(self.file_path('skewT_{}.png'.format(expt)))
