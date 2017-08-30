from logging import getLogger

import matplotlib

matplotlib.use('Agg')
import numpy as np
import pylab as plt
import iris

from omnium.analyzer import Analyzer
from omnium.utils import get_cube
from omnium.consts import L

from analysis.colour import EXPT_COLOUR

logger = getLogger('om.prof_an')


class SurfFluxPlot(Analyzer):
    analysis_name = 'surf_flux_plot'
    multi_expt = True
    expts_to_plot = ['S0', 'S4']

    def _plot(self):
        lines = ['Expt,PFE [W m-2],LHF [W m-2],SHF [W m-2]']
        for expt in self.expts:
            cubes = self.expt_cubes[expt]

            precip = get_cube(cubes, 4, 203)
            lhf = get_cube(cubes, 3, 234)
            shf = get_cube(cubes, 3, 217)

            precip_ts = precip.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
            lhf_ts = lhf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
            shf_ts = shf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

            if expt in self.expts_to_plot:
                start_time = precip_ts.coord('time').points[0]
                times = (precip_ts.coord('time').points - start_time) / 24

                colour = EXPT_COLOUR[expt]
                plot = plt.plot(times, lhf_ts.data, label=expt, color=colour, linestyle='-')
                #colour = plot[0].get_color()
                colour = plot[0].get_color()
                plt.plot(times, shf_ts.data, color=colour, linestyle='-.')
                # TODO: fix for any time delta.
                # daily smoothing with 15 min means.
                precip_ts_smoothed = np.convolve(precip_ts.data, np.ones((96, )) / 96., mode='same')
                plt.plot(times[96:-96], precip_ts_smoothed[96:-96] * L, color=colour, linestyle='--')

            day20index = len(precip_ts.data) / 2
            mean_pfe = precip_ts.data[day20index:].mean() * L
            mean_lhf = lhf_ts.data[day20index:].mean()
            mean_shf = shf_ts.data[day20index:].mean()
            lines.append('{},{},{},{}'.format(expt, mean_pfe, mean_lhf, mean_shf))

        self.save_text('energy_flux.csv', '\n'.join(lines) + '\n')
        plt.ylim((-50, 300))
        plt.xlim((0, 40))
        plt.legend()

        plt.ylabel('flux (W m$^{-2}$)')
        plt.xlabel('time (day)')
        plt.axvline(x=20, linestyle='--', color='k')
        plt.savefig(self.figpath('energy_fluxes.png'))

    def run_analysis(self):
        pass

    def display_results(self):
        """Save all results for surf flux analysis."""
        self._plot()

