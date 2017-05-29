import os
from collections import OrderedDict
from itertools import groupby

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MassFluxPlotter(Analyzer):
    analysis_name = 'mass_flux_plot'
    multi_expt = True

    def run_analysis(self):
        pass

    def _plot_mass_flux_hist(self):
	self.append_log('plotting MSEs')

	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            sorted_cubes = []

            for cube in cubes:
                (height_level_index, thresh_index) = cube.attributes['mass_flux_key']
                mf_key = (height_level_index, thresh_index)
                sorted_cubes.append((mf_key, cube))

            # Each element is a tuple like: ((1, 3), cube)
            # Sorting will put in correct order, sorting on initial tuple.
            sorted_cubes.sort()

            # Group on first element of tuple, i.e. on 1 for ((1, 3), cube)
            for group, cubes in groupby(sorted_cubes, lambda x: x[0][0]):
                hist_data = []
                dmax = 0
                for i, item in enumerate(cubes):
                    cube = item[1]
                    hist_data.append(cube)
                    dmax = max(cube.data.max(), dmax)

                assert len(hist_data) == 3
                name = '{}.{}.z{}.mass_flux_hist'.format(self.output_filename, expt, group)
                plt.figure(name)

                # TODO: User configurable.
                y_min, bin_edges = np.histogram(hist_data[2].data, bins=50, range=(0, dmax))
                y, bin_edges = np.histogram(hist_data[1].data, bins=50, range=(0, dmax))
                y_max, bin_edges = np.histogram(hist_data[0].data, bins=50, range=(0, dmax))

                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

                plot_filename = os.path.join(self.results_dir, name + '.png')
                # yerr is a rel, not abs, value.
                plt.bar(bin_centers, y, yerr=[y - y_min, y_max - y])
                plt.savefig(plot_filename)
                self.append_log('Saved to {}'.format(plot_filename))


    def save_analysis(self):
        self._plot_mass_flux_hist()
        plt.close('all')
