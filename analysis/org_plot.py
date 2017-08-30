from itertools import groupby

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt

from omnium.analyzer import Analyzer

from analysis.utils import cm_to_inch


LX = 256000
LY = 256000


class OrgPlotter(Analyzer):
    analysis_name = 'org_plot'
    multi_expt = True

    def set_config(self, config):
        super(OrgPlotter, self).set_config(config)
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

    def _plot_org_hist(self):
        self.append_log('plotting org')

        groups = []

        for expt in self.expts:
            cubes = self.expt_cubes[expt]
            sorted_cubes = []

            for cube in cubes:
                (height_level_index, thresh_index) = cube.attributes['dist_key']
                dist_key = (height_level_index, thresh_index)
                sorted_cubes.append((dist_key, cube))

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
                name = '{}.z{}.dist_hist'.format(expt, group)
                plt.figure(name)
                plt.clf()
                #plt.title(name)

                hist_kwargs = {}
                if self.xlim:
                    hist_kwargs['range'] = self.xlim
                else:
                    hist_kwargs['range'] = (0, dmax)

                if self.nbins:
                    hist_kwargs['bins'] = self.nbins
                #y, bin_edges = np.histogram(hist_data[1].data, **hist_kwargs)

                #n, bins = np.histogram(hist_data[1].data, **hist_kwargs)
                plt.figure('Not_used')
                n, bins, patch = plt.hist(hist_data[1].data, 700)
                plt.figure('combined_expt_z{}'.format(group))

                areas = np.pi * (bins[1:]**2 - bins[:-1]**2)
                cloud_densities = n / areas

                # Normalize based on mean density over domain.
                # WRONG WAY TO DO IT!:
                # mean = cloud_densities[bins < LX / 2.].mean()
                # self.plt.plot((bins[:-1] + bins[1:]) / 2, cloud_densities / mean)
                # calculates the mean of the densities, not the mean density.

                # Correct way to normalize:
                # Divide the total number in a circle by the circle's area.
                imax = np.argmax(bins[1:] > (LX / 2))
                mean_density = n[:imax].sum() / (np.pi * bins[imax]**2)
                xpoints = (bins[:-1] + bins[1:]) / 2

                plt.figure('combined_expt_z{}'.format(group))
                plt.plot(xpoints / 1000, cloud_densities / mean_density, label=expt)
                plt.xlabel('Distance (km)')
                plt.ylabel('Normalized cloud number density')
                plt.axhline(y=1, ls='--')
                #print(bins[:21] / 1000)
                #print(n[:20])
                #print(cloud_densities[:20])
                plt.xlim((0, 120))
                plt.ylim((0, 20))

                plt.figure('combined_expt_z{}_log'.format(group))
                plt.plot(xpoints / 1000, cloud_densities / mean_density, label=expt)
                plt.yscale('log')
                plt.xlabel('Distance (km)')
                plt.ylabel('Normalized cloud number density')
                plt.axhline(y=1, ls='--')

                plt.xlim((0, 256))
                plt.ylim((1e-1, 2e1))

                plt.figure('poster_combined_expt_z{}_log'.format(group))
                plt.plot(xpoints / 1000, cloud_densities / mean_density, label=expt)
                plt.yscale('log')
                plt.xlabel('Distance (km)')
                plt.ylabel('Normalized cloud\nnumber density')
                plt.axhline(y=1, ls='--')

                plt.xlim((0, 256))
                plt.ylim((1e-1, 2e1))

        for group in groups:
            plt.figure('combined_expt_z{}'.format(group))
            #plt.title('combined_expt_z{}'.format(group))
            plt.legend(loc='upper right')
            plt.savefig(self.figpath('z{}_combined.png'.format(group)))

            plt.figure('combined_expt_z{}_log'.format(group))
            plt.legend(loc='upper right')
            plt.savefig(self.figpath('z{}_combined_log.png'.format(group)))

            fig = plt.figure('poster_combined_expt_z{}_log'.format(group))
            fig.set_size_inches(*cm_to_inch(25, 7))
            plt.legend(loc='upper center', ncol=5)
            plt.tight_layout()
            plt.savefig(self.figpath('poster_z{}_combined_log.png'.format(group)))

    def display_results(self):
        self._plot_org_hist()
        plt.close('all')

