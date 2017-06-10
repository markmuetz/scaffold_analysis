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

    def set_config(self, config):
	super(MassFluxPlotter, self).set_config(config)
        if 'xlim' in config:
            self.xlim = [float(v) for v in config['xlim'].split(',')]
        else:
            self.xlim = None

        if 'ylim' in config:
            self.ylim = [float(v) for v in config['ylim'].split(',')]
        else:
            self.ylim = None
        self.nbins = config.getint('nbins', None)

    def run_analysis(self):
        pass

    def _plot_mass_flux_hist(self):
	self.append_log('plotting mass_flux')

        groups = []

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
                if group not in groups:
                    groups.append(group)
                hist_data = []
                dmax = 0
                for i, item in enumerate(cubes):
                    cube = item[1]
                    hist_data.append(cube)
                    dmax = max(cube.data.max(), dmax)

                assert len(hist_data) == 3
                name = '{}.z{}.mass_flux_hist'.format(expt, group)
                plt.figure(name)
                plt.clf()
                plt.title(name)

                hist_kwargs = {}
                if self.xlim:
                    hist_kwargs['range'] = self.xlim
                else:
                    hist_kwargs['range'] = (0, dmax)

                if self.nbins:
                    hist_kwargs['bins'] = self.nbins
                #y_min, bin_edges = np.histogram(hist_data[2].data, bins=50, range=(0, dmax))
                #y_max, bin_edges = np.histogram(hist_data[0].data, bins=50, range=(0, dmax))
                y, bin_edges = np.histogram(hist_data[1].data, **hist_kwargs)
                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
		y2 = bin_centers * y

                # yerr is a rel, not abs, value.
                # N.B. full width bins.
                width = bin_edges[1:] - bin_edges[:-1]
                plt.bar(bin_centers, y, width=width)
                #plt.bar(bin_centers, y, width=width, yerr=[y - y_min, y_max - y])

                if self.xlim:
                    plt.xlim(self.xlim)
                plt.yscale('log')
                if self.ylim:
                    plt.ylim(ymax=self.ylim[1])

                plt.yscale('linear')
                if self.ylim:
                    plt.ylim(self.ylim)
                plt.savefig(self.figpath(name + '.png'))

		plt.figure(name + 'mf_wieghted_plot_filename')
		plt.clf()
		plt.plot(bin_centers, y2)
                plt.savefig(self.figpath(name + '.mf_weighted.png'))

                plt.figure('combined_expt_z{}'.format(group))
                plt.plot(bin_centers, y, label=expt)

                plt.figure('combined_expt_mf_weighted_z{}'.format(group))
		plt.plot(bin_centers, y2, label=expt)

        for group in groups:
            plt.figure('combined_expt_z{}'.format(group))
            plt.title('combined_expt_z{}'.format(group))
            plt.legend()
            plt.yscale('log')
            plt.savefig(self.figpath('z{}_combined.png'.format(group)))

	    plt.figure('combined_expt_mf_weighted_z{}'.format(group))
            plt.legend()
            plt.savefig(self.figpath('z{}_mf_weighted_comb.png'.format(group)))

    def display_results(self):
        self._plot_mass_flux_hist()
        plt.close('all')
