from logging import getLogger

import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium import Analyser
from omnium.utils import get_cube
from omnium.consts import p_ref, kappa, L, cp, g

from scaffold.expt_settings import EXPT_DETAILS

logger = getLogger('scaf.rel_plot')


class RelaxationPlot(Analyser):
    analysis_name = 'relaxation_plot'
    multi_expt = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmos.456.pp2.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.relaxation_plot.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self._plot_theta_mv_incs()
        self._plot_UCP_theta_mv_incs()
        self._calc_total_heating()

    def _calc_total_heating(self):
        total_heating = ['Expt,TrelFE [W m-2],QrelFE [W m-2]']
        for expt in self.task.expts:
            cube = self.expt_cubes[expt]

            # incs on theta levels.
            T_inc = get_cube(cube, 53, 181)
            mv_inc = get_cube(cube, 53, 182)

            # exnerp on rho levels.
            exnerp = get_cube(cube, 0, 255)
            p = exnerp.data ** (1 / kappa) * p_ref

            # N.B. this is the correct way round if you want the right signs.
            # TODO: this calculation is not exact. Midpoints of rho levels are not on theta
            # TODO: this calculation also assumes hydrostatic equilibrium - this is not right
            # TODO: this calculation also doesn't take into account theta-level 0 - not sure about?
            # levels so, working out dp should take this into account.
            dp = (p[:, :-1] - p[:, 1:])

            # z = T_inc.coord('level_height').points
            # z_theta = exnerp.coord('level_height').points

            # Note, I am ignoring top level!
            T_inc_fe = ((cp / g) * (T_inc.data[:, :-1] / 30 * dp)).sum(axis=1).mean(axis=(1, 2))
            q_inc = mv_inc.data / (1 - mv_inc.data)
            q_inc_fe = ((L / g) * (q_inc[:, :-1] / 30 * dp)).sum(axis=1).mean(axis=(1, 2))

            logger.debug('{} - T_inc flux equiv: {} W m-2', expt, T_inc_fe)
            logger.debug('{} - q_inc flux equiv: {} W m-2', expt, q_inc_fe)
            total_heating.append('{},{},{}'.format(expt, T_inc_fe.mean(), q_inc_fe.mean()))
        self.save_text('final_day_relaxation_energy_flux.csv', '\n'.join(total_heating) + '\n')
        self.save_text('final_day_relaxation_energy_flux.csv.done', 'done')

    def _plot_theta_mv_incs(self):
        plt.clf()
        plt.xlabel('T_inc (K day$^{-1}$)')
        plt.ylabel('height (km)')
        for expt in self.task.expts:
            T_inc = get_cube(self.expt_cubes[expt], 53, 181)

            z = T_inc.coord('level_height').points
            colour = None

            for i in range(8):
                print(i)
                # 2880 steps_per_perdiodm: converts to K/day
                T_inc_profile = T_inc[i].data.mean(axis=(1, 2)) * 2880
                if i == 0:
                    plot = plt.plot(T_inc_profile, z / 1000, label=expt)
                    colour = plot[0].get_color()
                else:
                    plt.plot(T_inc_profile, z / 1000, color=colour)

        plt.legend()
        plt.savefig(self.file_path('T_incs.png'))

        plt.clf()
        plt.xlabel('mv_inc (g kg$^{-1}$ day$^{-1}$)')
        plt.ylabel('height (km)')
        for expt in self.task.expts:
            mv_inc = get_cube(self.expt_cubes[expt], 53, 182)
            z = mv_inc.coord('level_height').points
            colour = None

            for i in range(8):
                print(i)
                # 2880 steps_per_perdiodm: converts to g/kg/day
                mv_inc_profile = mv_inc[i].data.mean(axis=(1, 2)) * 2880 * 1000
                if i == 0:
                    plot = plt.plot(mv_inc_profile, z / 1000, label=expt)
                    colour = plot[0].get_color()
                else:
                    plt.plot(mv_inc_profile, z / 1000, color=colour)

        plt.legend()
        plt.savefig(self.file_path('mv_incs.png'))

    def _plot_UCP_theta_mv_incs(self):
        plt.clf()
        plt.xlabel('T_inc (K day$^{-1}$)')
        plt.ylabel('height (km)')
        for expt in self.task.expts:
            ucp_kwargs = {}
            if expt in EXPT_DETAILS:
                ucp_kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
            T_inc = get_cube(self.expt_cubes[expt], 53, 181)

            z = T_inc.coord('level_height').points

            # 2880 steps_per_perdiodm: converts to K/day
            T_inc_profile = T_inc.data.mean(axis=(0, 2, 3)) * 2880
            plt.plot(T_inc_profile, z / 1000, **ucp_kwargs)

        plt.xlim((-3, 1))
        plt.ylim((0, 15))
        plt.axvline(x=0, ls='--', color='k')
        plt.legend(loc='lower right')
        plt.savefig(self.file_path('UCP_T_incs.png'))

        plt.clf()
        plt.xlabel('mv_inc (g kg$^{-1}$ day$^{-1}$)')
        plt.ylabel('height (km)')
        for expt in self.task.expts:
            ucp_kwargs = {}
            if expt in EXPT_DETAILS:
                ucp_kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
            mv_inc = get_cube(self.expt_cubes[expt], 53, 182)
            z = mv_inc.coord('level_height').points

            # 2880 steps_per_perdiodm: converts to g/kg/day
            mv_inc_profile = mv_inc.data.mean(axis=(0, 2, 3)) * 2880 * 1000
            plt.plot(mv_inc_profile, z / 1000, **ucp_kwargs)

        plt.ylim((0, 15))
        plt.legend(loc='upper right')
        plt.savefig(self.file_path('UCP_mv_incs.png'))
