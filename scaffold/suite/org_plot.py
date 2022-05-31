import os
from itertools import groupby
from logging import getLogger

import matplotlib
import numpy as np
import scipy.interpolate as interpolate
import scipy.integrate as integrate

matplotlib.use('Agg')
import pylab as plt

from omnium import Analyser, ExptList

from scaffold.utils import cm_to_inch, find_intersections
from scaffold.expt_settings import EXPT_DETAILS

logger = getLogger('scaf.org_plot')


class OrgPlotter(Analyser):
    """Performs organization analysis similar to that in Cohen & Craig 2006.

    Calcs and plots a 'histogram' of cloud-cloud distances.
    Uses cloud-cloud distances, divided by area of annulus for each bin, and normalalizes by the
    area avg cloud density (to get a value that would be 1 if cloud field randomly distributed.
    """
    analysis_name = 'org_plot'
    multi_expt = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename = '{input_dir}/atmos.org_combined.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.org_plot.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        self.xlim = None
        self.ylim = None
        self.nbins = 128

        self.groups = []
        self.org_data = {}
        self.thresh = 0  # 0, 1 or 2 -- 10% lower, actual, 10% higher.

        expts = ExptList(self.suite)
        expts.find(self.task.expts)
        for expt in self.task.expts:
            expt_obj = expts.get(expt)
            cubes = self.expt_cubes[expt]
            sorted_cubes = []
            lx = expt_obj.lx

            for cube in cubes:
                (height_level_index, thresh_index) = cube.attributes['dist_key']
                dist_key = (height_level_index, thresh_index)
                sorted_cubes.append((dist_key, cube))

            # Each element is a tuple like: ((1, 3), cube)
            # Sorting will put in correct order, sorting on initial tuple.
            sorted_cubes.sort()

            # Group on first element of tuple, i.e. on 1 for ((1, 3), cube)
            for group, cubes in groupby(sorted_cubes, lambda x: x[0][0]):
                if group not in self.groups:
                    self.groups.append(group)
                hist_data = []
                dmax = 0
                for i, item in enumerate(cubes):
                    cube = item[1]
                    hist_data.append(cube)
                    dmax = max(cube.data.max(), dmax)

                assert len(hist_data) == 3

                hist_kwargs = {}
                if self.xlim:
                    hist_kwargs['range'] = self.xlim
                else:
                    hist_kwargs['range'] = (0, dmax)

                if self.nbins:
                    hist_kwargs['bins'] = self.nbins

                plt.figure('Not_used')
                n, bins, patch = plt.hist(hist_data[self.thresh].data, **hist_kwargs)
                plt.figure('combined_expt_z{}'.format(group))

                areas = np.pi * (bins[1:]**2 - bins[:-1]**2)
                cloud_densities = n / areas

                # Correct way to normalize:
                # Divide the total number in a circle by the circle's area.
                imax = np.argmax(bins[1:] > (lx / 2))
                mean_density = n[:imax].sum() / (np.pi * bins[imax]**2)
                xpoints = (bins[:-1] + bins[1:]) / 2

                org_data = cloud_densities / mean_density
                self.org_data[(expt, group)] = {'xpoints': xpoints, 'org_data': org_data}
                ones = np.ones_like(org_data)
                indices, weights = find_intersections(org_data, ones)
                assert np.all(org_data[:indices[0]] > 0), 'Expected init org_data to be +ve'

                org_data_fn = interpolate.interp1d(xpoints, org_data - ones)
                crosses = []
                for i, index in enumerate(indices[:-1]):
                    cross = (xpoints[index] + weights[i] * (xpoints[index + 1] - xpoints[index]))

                    assert np.isclose(org_data_fn(cross), 0), 'Cross point {} not 0'.format(index)
                    crosses.append(cross)

                logger.info(crosses)
                cluster_radius = crosses[0]
                suppr_radius = crosses[1]
                logger.info('cluster_radius: {}', cluster_radius)
                logger.info('suppr_radius: {}', suppr_radius)

                x1 = list(xpoints[:indices[0] + 1]) + [crosses[0]]
                y1 = org_data_fn(x1)
                cluster_index = integrate.trapz(y1, x1)
                logger.debug('cluster_index: {}', cluster_index)

                x2 = list(xpoints[indices[0] + 1:indices[1] + 1]) + [crosses[1]]
                y2 = org_data_fn(x2)

                # Make +ve.
                suppr_index = np.abs(integrate.trapz(y2, x2))
                logger.debug('suppr_index: {}', suppr_index)

                self.org_data[(expt, group)]['crosses'] = crosses
                self.org_data[(expt, group)]['cluster_radius'] = cluster_radius
                self.org_data[(expt, group)]['suppr_radius'] = suppr_radius
                self.org_data[(expt, group)]['cluster_index'] = cluster_index
                self.org_data[(expt, group)]['suppr_index'] = suppr_index

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        csv_lines = ['expt,group,cluster_radius,suppr_radius,cluster_index,suppr_index']
        for (expt, group), org_data_item in self.org_data.items():
            org_data = org_data_item['org_data']
            xpoints = org_data_item['xpoints']
            self._plot_org_hist(expt, group, xpoints, org_data)

            org_data_filename = os.path.join(os.path.dirname(self.task.output_filenames[0]),
                                             'org_data_{}_z{}.np'.format(expt, group))
            xpoints_filename = os.path.join(os.path.dirname(self.task.output_filenames[0]),
                                            'org_data_xpoints_{}_z{}.np'.format(expt, group))

            org_data.dump(org_data_filename)
            xpoints.dump(xpoints_filename)
            csv_lines.append('{expt},{group},{cluster_radius},{suppr_radius},{cluster_index},{suppr_index}'
                             .format(expt=expt, group=group, **org_data_item))
        org_csv_filename = os.path.join(os.path.dirname(self.task.output_filenames[0]),
                                        'org_data.csv')
        with open(org_csv_filename, 'w') as f:
            f.write('\n'.join(csv_lines) + '\n')

        for group in self.groups:
            plt.figure('combined_expt_z{}'.format(group))
            plt.legend(loc='upper right')
            plt.savefig(self.file_path('z{}_t{}_combined.png'.format(group, self.thresh)))

            plt.figure('combined_expt_z{}_log'.format(group))
            plt.legend(loc='upper right')
            plt.savefig(self.file_path('z{}_t{}_combined_log.png'.format(group, self.thresh)))

            fig = plt.figure('poster_combined_expt_z{}_log'.format(group))
            fig.set_size_inches(*cm_to_inch(25, 7))
            plt.legend(loc='upper center', ncol=5)
            plt.tight_layout()
            plt.savefig(self.file_path('poster_z{}_t{}_combined_log.png'.format(group, self.thresh)))

            fig = plt.figure('fig_combined_expt_z{}_log'.format(group))
            fig.set_size_inches(*cm_to_inch(20, 7))
            plt.legend(loc='upper right', ncol=2)
            plt.tight_layout()
            plt.savefig(self.file_path('fig_z{}_t{}_combined_log.png'.format(group, self.thresh)))

        plt.close('all')

    def _plot_org_hist(self, expt, group, xpoints, org_data):
        self.append_log('plotting org')

        ucp_kwargs = {}
        if expt in EXPT_DETAILS:
            ucp_kwargs = dict(zip(['label', 'color', 'linestyle'], EXPT_DETAILS[expt]))

        plt.figure('combined_expt_z{}'.format(group))
        plt.plot(xpoints / 1000, org_data, label=expt)
        plt.xlabel('Distance (km)')
        # plt.ylabel('Normalized cloud number density')
        plt.ylabel('RDF')
        plt.axhline(y=1, color='k', ls='--')

        plt.xlim((0, 120))
        plt.ylim((0, 20))

        plt.figure('combined_expt_z{}_log'.format(group))
        plt.plot(xpoints / 1000, org_data, label=expt)
        plt.yscale('log')
        plt.xlabel('Distance (km)')
        # plt.ylabel('Normalized cloud number density')
        plt.ylabel('RDF')
        plt.axhline(y=1, ls='--')

        plt.xlim((0, 256))
        plt.ylim((1e-1, 2e1))

        plt.figure('poster_combined_expt_z{}_log'.format(group))
        plt.plot(xpoints / 1000, org_data, label=expt)
        plt.yscale('log')
        plt.xlabel('Distance (km)')
        plt.ylabel('Normalized cloud\nnumber density')
        plt.axhline(y=1, ls='--')

        plt.xlim((0, 256))
        plt.ylim((1e-1, 2e1))

        plt.figure('fig_combined_expt_z{}_log'.format(group))
        plt.plot(xpoints / 1000, org_data, **ucp_kwargs)
        plt.yscale('log')
        plt.xlabel('Distance (km)')
        # plt.ylabel('Normalized cloud\nnumber density')
        plt.ylabel('RDF')
        plt.axhline(y=1, ls='--', color='k')

        plt.xlim((0, 128))
        if len(self.task.expts) == 10:
            # RWPs.
            plt.ylim((2e-1, 2e1))
        else:
            plt.ylim((6e-1, 1e1))
