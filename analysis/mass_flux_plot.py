import os
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MassFlugPlotter(Analyzer):
    analysis_name = 'mass_flux_plot'
    multi_expt = True

    def run_analysis(self):
	self.expt_mses = OrderedDict()
	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            for cube in cubes:
                print(cube.name())

    def _plot_mass_flux_hist(self):
	self.append_log('plotting MSEs')

	self.expt_mses = OrderedDict()
	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            for cube in cubes:
                name = self.output_filename + cube.name() + '_mass_flux_hist'
                plt.figure(name)
                plt.hist(cube.data, bins=50)

                plot_filename = os.path.join(self.results_dir, name + '.png')
                plt.savefig(plot_filename)
                self.append_log('Saved to {}'.format(plot_filename))

    def save_analysis(self):
        self._plot_mass_flux_hist()
        plt.close('all')
