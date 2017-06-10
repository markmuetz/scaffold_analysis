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

    def _plot(self):
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


    def display_results(self):
        self._plot()
        plt.close('all')
