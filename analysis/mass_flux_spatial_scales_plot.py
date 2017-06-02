import os
from collections import OrderedDict
from itertools import groupby

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MassFluxSpatialScalesPlotter(Analyzer):
    analysis_name = 'mass_flux_spatial_scales_plot'
    multi_expt = True

    def read_lim(self, lim):
        return [float(v) for v in lim.split(',')]

    def set_config(self, config):
	super(MassFluxSpatialScalesPlotter, self).set_config(config)
        self.nbins = config.getint('nbins', None)

    def run_analysis(self):
        pass

    def _plot_mass_flux_spatial(self):
	self.append_log('plotting mass_flux_spatial')

        heights = []
        ns = []

	for expt in self.expts:
	    cubes = self.expt_cubes[expt]
            sorted_cubes = []

            for cube in cubes:
                (height_level_index, thresh_index, n) = cube.attributes['mass_flux_spatial_key']
                mf_key = (height_level_index, thresh_index, n)
                sorted_cubes.append((mf_key, cube))

            # Each element is a tuple like: ((1, 2, 32), cube)
            # Sorting will put in correct order, sorting on initial tuple.
            sorted_cubes.sort()

            # Group on first element of tuple, i.e. on 1 for ((1, 2, 32), cube)
            for height_index, key_cubes in groupby(sorted_cubes, lambda x: x[0][0]):
                if height_index not in heights:
                    heights.append(height_index)
                hist_data = []
                dmax = 0
                for i, key_cube in enumerate(key_cubes):
                    # middle cube is the one with the middle thresh_index.
                    mf_key = key_cube[0]
                    cube = key_cube[1]
                    # Pick out middle element, i.e. thresh_index == 1.
                    if mf_key[1] == 1:
                        hist_data.append((mf_key, cube))
                        dmax = max(cube.data.max(), dmax)

                # assert len(hist_data) == 3
                for mf_key, hist_datum in hist_data:
                    (height_index, thresh_index, n) = mf_key
                    if n not in ns:
                        ns.append(n)
                    name = '{}.{}.z{}.n{}.mass_flux_spatial_hist'.format(self.output_filename, 
			                                                 expt, height_index, n)
                    plt.figure(name)
                    plt.clf()
                    plt.title('{} z{} n{} mass_flux_spatial_hist'.format(expt, height_index, n))

                    hist_kwargs = {}
		    xlim_key = 'z{}_n{}_xlim'.format(height_index, n)
		    ylim_key = 'z{}_n{}_ylim'.format(height_index, n)
		    xlim = None
		    ylim = None
		    if xlim_key in self._config:
			xlim = self.read_lim(self._config[xlim_key])
                        hist_kwargs['range'] = xlim
                    else:
                        #hist_kwargs['range'] = (0, 0.1)
			pass

		    if ylim_key in self._config:
			ylim = self.read_lim(self._config[ylim_key])

                    if self.nbins:
                        hist_kwargs['bins'] = self.nbins
		    filtered_data = hist_datum.data[hist_datum.data > 0.01]
                    y, bin_edges = np.histogram(filtered_data, **hist_kwargs)
                    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

                    # N.B. full width bins.
                    width = bin_edges[1:] - bin_edges[:-1]
                    plt.bar(bin_centers, y, width=width)

		    plt.xlim((0, 0.1))
                    if xlim:
                        plt.xlim(xlim)
                    if ylim:
                        plt.ylim(ylim)
                    plt.savefig(self.figpath('.png'))

                    plt.figure('combined_expt_z{}_n{}'.format(height_index, n))
                    plt.plot(bin_centers, y, label=expt)

        for height_index in heights:
	    for n in ns:
		plt.figure('combined_expt_z{}_n{}'.format(height_index, n))
		plt.title('combined_expt_z{}_n{}'.format(height_index, n))
		plt.legend()
		plt.xlim((0, 0.1))
		plt.savefig(self.figpath('_z{}_n{}_combined.png'.format(height_index, n)))

    def display_results(self):
        self._plot_mass_flux_spatial()
        plt.close('all')
