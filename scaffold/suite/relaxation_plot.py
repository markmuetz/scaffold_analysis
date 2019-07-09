from logging import getLogger
import warnings

import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium import Analyser, ExptList
from omnium.utils import get_cube, cm_to_inch
from omnium.consts import p_ref, kappa, L, cp, g, Re

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
        self._plot_T_mv_incs()
        self._plot_UCP_T_mv_incs()
        self._plot_combined_T_mv_incs()
        self._calc_total_heating()

    def _calc_total_heating(self):
        total_heating = ['Expt,relFE [W m-2],QrelFE [W m-2]']
        expts = ExptList(self.suite)
        expts.find(self.task.expts)
        for expt in self.task.expts:
            expt_obj = expts.get(expt)
            dt = expt_obj.dt

            cube = self.expt_cubes[expt]

            # incs on theta levels.
            T_inc = get_cube(cube, 53, 181)
            mv_inc = get_cube(cube, 53, 182)

            self._calc_vertical_integral(expt, cube, T_inc, mv_inc, dt, total_heating)
            # self._calc_vertical_integral_hydrostatic(expt, cube, T_inc, mv_inc, dt, total_heating)

        self.save_text('final_day_relaxation_energy_flux.csv', '\n'.join(total_heating) + '\n')
        self.save_text('final_day_relaxation_energy_flux.csv.done', 'done')

    def _calc_vertical_integral(self, expt, cube, T_inc, mv_inc, dt, total_heating):
        rho = get_cube(cube, 0, 253)
        rho.data = rho.data / Re**2
        for coord in T_inc.coords():
            if coord.name() == 'level_height':
                height_name = 'level_height'
                break
            elif coord.name() == 'atmosphere_hybrid_height_coordinate':
                height_name = 'atmosphere_hybrid_height_coordinate'
                break
        z_coord = T_inc.coord(height_name)
        dz = z_coord.bounds[:, 1] - z_coord.bounds[:, 0]
        T_inc_fe = (cp * (rho.data * T_inc.data / dt * dz[None, :, None, None])).sum(axis=1).mean(axis=(1, 2))
        q_inc = mv_inc.data / (1 - mv_inc.data)
        q_inc_fe = (L * (rho.data * q_inc / dt * dz[None, :, None, None])).sum(axis=1).mean(axis=(1, 2))
        logger.debug('{} - T_inc flux equiv: {} W m-2', expt, T_inc_fe)
        logger.debug('{} - q_inc flux equiv: {} W m-2', expt, q_inc_fe)
        # total_heating.append('{},nonhydrostatic,{},{}'.format(expt, T_inc_fe.mean(), q_inc_fe.mean()))
        total_heating.append('{},{},{}'.format(expt, T_inc_fe.mean(), q_inc_fe.mean()))

    def _calc_vertical_integral_hydrostatic(self, expt, cube, T_inc, mv_inc, dt, total_heating):
        msg = ("This method has some incorrect assumptions in it\n"
               "Please use _calc_vertical_integral instead")
        warnings.warn(msg, DeprecationWarning)
        # See ASSUMPTION
        # N.B. difference is small, but using the correct calc leads to better cons. of energy in
        # budget calcs.

        # exnerp on rho levels.
        exnerp = get_cube(cube, 0, 255)
        p = exnerp.data ** (1 / kappa) * p_ref

        # ASSUMPTION: this calculation is not exact. Midpoints of rho levels are not on theta
        # ASSUMPTION: this calculation also assumes hydrostatic equilibrium - this is not right
        # ASSUMPTION: this calculation also doesn't take into account theta-level 0 - not sure about?
        # N.B. this is the correct way round if you want the right signs.
        dp = (p[:, :-1] - p[:, 1:])
        # z = T_inc.coord('level_height').points
        # z_theta = exnerp.coord('level_height').points
        # Note, I am ignoring top level!
        T_inc_fe = ((cp / g) * (T_inc.data[:, :-1] / dt * dp)).sum(axis=1).mean(axis=(1, 2))
        q_inc = mv_inc.data / (1 - mv_inc.data)
        q_inc_fe = ((L / g) * (q_inc[:, :-1] / dt * dp)).sum(axis=1).mean(axis=(1, 2))
        logger.debug('{} - T_inc flux equiv: {} W m-2', expt, T_inc_fe)
        logger.debug('{} - q_inc flux equiv: {} W m-2', expt, q_inc_fe)
        total_heating.append('{},{},{}'.format(expt, T_inc_fe.mean(), q_inc_fe.mean()))

    def _plot_T_mv_incs(self):
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

    def _plot_UCP_T_mv_incs(self):
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

        plt.xlim((-5, 1))
        plt.ylim((0, 15))
        plt.axvline(x=0, ls='--', color='k')
        plt.legend(loc='lower left')
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

    def _plot_combined_T_mv_incs(self):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=cm_to_inch(20, 12))

        ax1.set_xlabel('T_inc (K day$^{-1}$)')
        ax1.set_ylabel('height (km)')
        for expt in self.task.expts:
            kwargs = {}
            if expt in EXPT_DETAILS:
                kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
            T_inc = get_cube(self.expt_cubes[expt], 53, 181)

            z = T_inc.coord('level_height').points

            # 2880 steps_per_perdiodm: converts to K/day
            T_inc_profile = T_inc.data.mean(axis=(0, 2, 3)) * 2880
            ax1.plot(T_inc_profile, z / 1000, **kwargs)

        ax1.set_xlim((-4, 1))
        ax1.set_ylim((0, 15))
        ax1.axvline(x=0, ls='--', color='k')

        ax2.set_xlabel('mv_inc (g kg$^{-1}$ day$^{-1}$)')
        for expt in self.task.expts:
            kwargs = {}
            if expt in EXPT_DETAILS:
                kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
            mv_inc = get_cube(self.expt_cubes[expt], 53, 182)
            z = mv_inc.coord('level_height').points

            # 2880 steps_per_perdiodm: converts to g/kg/day
            mv_inc_profile = mv_inc.data.mean(axis=(0, 2, 3)) * 2880 * 1000
            ax2.plot(mv_inc_profile, z / 1000, **kwargs)
            logger.info('{}: mv_inc max;min: {};{}',
                        expt, mv_inc_profile.max(), mv_inc_profile.min())

        ax2.set_xlim((-1, 0.4))
        ax2.set_ylim((0, 15))
        ax2.axvline(x=0, ls='--', color='k')
        ax2.legend(loc='upper left')

        ax3.set_xlabel('combined_inc (K day$^{-1}$)')
        for expt in self.task.expts:
            kwargs = {}
            if expt in EXPT_DETAILS:
                kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))
            T_inc = get_cube(self.expt_cubes[expt], 53, 181)
            mv_inc = get_cube(self.expt_cubes[expt], 53, 182)
            z = mv_inc.coord('level_height').points

            # 2880 steps_per_perdiodm: converts to K/day
            T_inc_profile = T_inc.data.mean(axis=(0, 2, 3)) * 2880
            # 2880 steps_per_perdiodm: converts to g/kg/day
            mv_inc_profile = mv_inc.data.mean(axis=(0, 2, 3)) * 2880 * 1000
            ax3.plot(T_inc_profile + L / cp * mv_inc_profile / 1000, z / 1000, **kwargs)

        ax3.set_xlim((-4, 1))
        ax3.set_ylim((0, 15))
        ax3.axvline(x=0, ls='--', color='k')

        plt.savefig(self.file_path('combined_T_mv_incs.png'))
