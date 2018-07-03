from itertools import groupby

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt
from scipy.stats import linregress

from omnium.analyser import Analyser

from scaffold.utils import cm_to_inch
from scaffold.suite_settings import dx
from scaffold.scaffold_settings import settings


class MassFluxPlotter(Analyser):
    """Plots mass flux histograms for various combinations of height_level, thresh_index.

    Fits a linear regression to the histograms (fitting a straight line to a log plot).
    """
    analysis_name = 'mass_flux_plot'
    multi_expt = True
    input_dir = 'omnium_output/{version_dir}'
    input_filename = '{input_dir}/{expt}/atmos.mass_flux_combined.nc'
    output_dir = 'omnium_output/{version_dir}/suite'
    output_filenames = ['{output_dir}/atmos.mass_flux_plot.dummy']

    settings = settings

    # Input values are in kg m-2 s-1, i.e. MF/cloud is an average over the cloud's area.
    # I want total MF/cloud though: multiply by the area of a grid cell or dx**2
    mass_flux_scaling = 1e8

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self.xlim = None
        self.ylim = None
        self.nbins = None
        self._plot_mass_flux_hist()
        plt.close('all')

    def _plot_mass_flux_hist(self):
        self.append_log('plotting mass_flux')

        groups = []

        linregress_details = ['z, expt, m, c, rval, pval, stderr']

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
                    dmax = max(cube.data.max() * dx**2 / self.mass_flux_scaling, dmax)

                assert len(hist_data) == 3
                name = '{}.z{}.mass_flux_hist'.format(expt, group)
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
                #y_min, bin_edges = np.histogram(hist_data[2].data, bins=50, range=(0, dmax))
                #y_max, bin_edges = np.histogram(hist_data[0].data, bins=50, range=(0, dmax))
                y, bin_edges = np.histogram(hist_data[1].data * dx**2 / self.mass_flux_scaling,
                                            **hist_kwargs)
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
                plt.savefig(self.file_path(name + '.png'))

                plt.figure(name + 'mf_wieghted_plot_filename')
                plt.clf()
                plt.plot(bin_centers, y2)
                plt.savefig(self.file_path(name + '.mf_weighted.png'))

                plt.figure('combined_expt_z{}'.format(group))
                plot = plt.plot(bin_centers, y, label=expt)
                colour = plot[0].get_color()
                # Rem. y = m * x + c
                #def fn(x, A, B):
                #    return A * np.exp(B * x)

                #popt, pcov = curve_fit(fn, bin_centers, y)

                log_y = np.log(y[y > 10])
                x = bin_centers[:len(log_y)]
                # TODO: Not a good idea to fit a straight line to exp data.
                # Would be better to do fit to explonential function surely.
                m, c, rval, pval, stderr = linregress(x[1:], log_y[1:])
                #import ipdb; ipdb.set_trace()
                #plt.plot(bin_centers, fn(bin_centers, *popt), color=colour, linestyle='--')
                plt.plot(x, np.exp(m * x + c), color=colour, linestyle='--')
                linregress_details.append('{},{},{},{},{},{},{}'.format(group, expt, m, c, rval, pval, stderr))

                plt.figure('combined_expt_mf_weighted_z{}'.format(group))
                plt.plot(bin_centers, y2, label=expt)

                if plt.fignum_exists('both_z{}'.format(group)):
                    f = plt.figure('both_z{}'.format(group))
                    ax1, ax2 = f.axes

                    # poster.
                    f_p = plt.figure('poster_z{}'.format(group))
                    ax1_p = f_p.axes[0]
                else:
                    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, num='both_z{}'.format(group))
                    if self.xlim:
                        ax1.set_xlim(self.xlim)
                    if self.ylim:
                        ax1.set_ylim(self.ylim)
                    ax1.set_yscale('log')
                    ax1.set_ylabel('Number of clouds')

                    # MUST MATCH self.mass_flux_scaling.
                    assert self.mass_flux_scaling == 1e8
                    ax2.set_ylabel('Mass flux contrib. ($\\times 10^8$ kg s$^{-1}$)')
                    #ax1.set_xlabel('MF per cloud ($\\times 10^7$ kg s$^{-1}$)')
                    ax2.set_xlabel('Mass flux per cloud ($\\times 10^8$ kg s$^{-1}$)')

                    f_p, ax1_p = plt.subplots(1, 1, num='poster_z{}'.format(group))
                    f_p.set_size_inches(*cm_to_inch(25, 9))
                    if self.xlim:
                        ax1_p.set_xlim(self.xlim)
                    if self.ylim:
                        ax1_p.set_ylim(self.ylim)
                    ax1_p.set_yscale('log')
                    ax1_p.set_ylabel('Number of clouds')
                    ax1_p.set_xlabel('Mass flux per cloud ($\\times 10^8$ kg s$^{-1}$)')

                plot = ax1.plot(bin_centers, y, label=expt)
                colour = plot[0].get_color()
                ax1.plot(x, np.exp(m * x + c), color=colour, linestyle='--')
                ax2.plot(bin_centers, y2, label=expt)

                ax1_p.plot(bin_centers, y, color=colour, label=expt)
                ax1_p.plot(x, np.exp(m * x + c), color=colour, linestyle='--')

        self.save_text('mf_linregress.csv', '\n'.join(linregress_details) + '\n')

        for group in groups:
            plt.figure('combined_expt_z{}'.format(group))
            #plt.title('combined_expt_z{}'.format(group))
            plt.legend()
            plt.yscale('log')
            plt.savefig(self.file_path('z{}_combined.png'.format(group)))

            plt.figure('combined_expt_mf_weighted_z{}'.format(group))
            plt.legend()
            plt.savefig(self.file_path('z{}_mf_weighted_comb.png'.format(group)))

            plt.figure('both_z{}'.format(group))
            plt.legend()
            plt.savefig(self.file_path('z{}_both.png'.format(group)))

            plt.figure('poster_z{}'.format(group))
            plt.tight_layout()
            plt.legend()
            plt.savefig(self.file_path('poster_z{}.png'.format(group)))