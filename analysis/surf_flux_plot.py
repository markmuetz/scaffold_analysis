import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
import iris
from matplotlib.colors import LogNorm

from omnium.analyzer import Analyzer
from omnium.utils import get_cube
from omnium.consts import L

from analysis.colour import EXPT_COLOUR


class SurfFluxPlot(Analyzer):
    analysis_name = 'surf_flux_plot'
    multi_expt = True

    def _plot(self):
	for expt in self.expts:
	    cubes = self.expt_cubes[expt]

            precip = get_cube(cubes, 4, 203)
            lhf = get_cube(cubes, 3, 234)
            shf = get_cube(cubes, 3, 217)

            precip_ts = precip.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
            lhf_ts = lhf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
            shf_ts = shf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

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

        plt.ylim((0, 350))
        plt.xlim((0, 20))
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


