from itertools import groupby
from logging import getLogger

import matplotlib
import numpy as np

matplotlib.use('Agg')
import pylab as plt
from scipy import stats
from scipy import optimize

from omnium import Analyser, ExptList

from scaffold.utils import cm_to_inch

logger = getLogger('scaf.mfp')


class PlotDensityMfContrib:
    """Plots density (top) and MF contrib (bottom) for all heights/expts."""
    def __init__(self, height_level_index, xlim, ylim, mass_flux_scaling=1e8):
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                       num='both_z{}'.format(height_level_index))

        if xlim:
            ax1.set_xlim(xlim)
        if ylim:
            ax1.set_ylim(ylim)
        ax1.set_yscale('log')
        ax1.set_ylabel('Cloud number\nprobability density')

        # MUST MATCH self.mass_flux_scaling.
        assert mass_flux_scaling == 1e8
        ax2.set_ylabel('Mass flux contrib. ($\\times 10^8$ kg s$^{-1}$)')
        #ax1.set_xlabel('MF per cloud ($\\times 10^7$ kg s$^{-1}$)')
        ax2.set_xlabel('Mass flux per cloud ($\\times 10^8$ kg s$^{-1}$)')
        self.fig, self.ax1, self.ax2 = fig, ax1, ax2

    def plot(self, expt_obj, bin_centers, width, y_density, y2, linregress, **kwargs):
        (x, m, c, rval, pval, stderr) = linregress

        plot = self.ax1.plot(bin_centers, y_density * width, label=expt_obj.name)
        colour = plot[0].get_color()
        # ax1.plot(x, np.exp(m * x + c), color=colour, linestyle='--')
        self.ax2.plot(bin_centers, y2, label=expt_obj.name)

    def finish(self, filename):
        self.ax1.legend()
        self.fig.savefig(filename)


class PlotDensityPoster:
    """Plots large fig of density for all heights/expts."""
    def __init__(self, height_level_index, xlim, ylim):
        fig, ax1 = plt.subplots(1, 1, num='poster_z{}'.format(height_level_index))
        fig.set_size_inches(*cm_to_inch(25, 9))
        if xlim:
            ax1.set_xlim(xlim)
        if ylim:
            ax1.set_ylim(ylim)
        ax1.set_yscale('log')
        ax1.set_ylabel('Cloud number probability density')
        ax1.set_xlabel('Mass flux per cloud ($\\times 10^8$ kg s$^{-1}$)')
        ax1.set_ylim((1e-4, 1))

        self.fig, self.ax1 = fig, ax1

    def plot(self, expt_obj, bin_centers, width, y_density, num_clds, **kwargs):
        self.ax1.plot(bin_centers, y_density * width,
                      label='{} - {} clouds'.format(expt_obj.name, num_clds))
        logger.debug('Sum y_density={}'.format((y_density * width).sum()))
        # ax1_p.fill_between(bin_centers, y_hist + np.sqrt(y_hist), y_hist - np.sqrt(y_hist),
        # color=colour, alpha=0.3)
        # ax1_p.plot(x, np.exp(m * x + c), color=colour, linestyle='--')

    def finish(self, filename):
        self.ax1.legend()
        self.fig.savefig(filename)


class MassFluxPlotter(Analyser):
    """Plots mass flux histograms for various combinations of height_level, thresh_index.

    Fits a linear regression to the histograms (fitting a straight line to a log plot).
    """
    analysis_name = 'mass_flux_plot'
    multi_expt = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename = '{input_dir}/atmos.mass_flux_combined.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.mass_flux_plot.dummy']

    # Input values are in kg m-2 s-1, i.e. MF/cloud is an average over the cloud's area.
    # I want total MF/cloud though: multiply by the area of a grid cell or dx**2
    mass_flux_scaling = 1e8

    def load(self):
        self.load_cubes()

    def run(self):
        self.xlim = (0, 3)
        self.ylim = None
        self.nbins = 100

        self._do_group_data()
        self._calc_histograms()
        self._calc_fit_linregress()
        self._calc_fit_expt()

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self._plot_mass_flux_hist()
        plt.close('all')

    def _do_group_data(self):
        """Group data into easier to use lists.

        self.height_levels contains indices of used height_levels
        self.all_expt_data has one item for each expt/height_level
        """
        self.height_levels = []
        self.all_expt_data = []

        expts = ExptList(self.suite)
        expts.find(self.task.expts)
        for expt in self.task.expts:
            expt_obj = expts.get(expt)
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
            # i.e. group on height_level_index.
            for height_level_index, cubes in groupby(sorted_cubes, lambda x: x[0][0]):
                if height_level_index not in self.height_levels:
                    self.height_levels.append(height_level_index)
                expt_data = {
                    'cubes': list(cubes),
                    'expt_obj': expt_obj,
                    'height_level_index': height_level_index
                }
                self.all_expt_data.append(expt_data)

    def _calc_histograms(self):
        for expt_data in self.all_expt_data:
            cubes = expt_data['cubes']
            expt_obj = expt_data['expt_obj']
            dx, dy = expt_obj.dx, expt_obj.dy

            hist_data = []
            dmax = 0
            for i, item in enumerate(cubes):
                cube = item[1]
                hist_data.append(cube)
                dmax = max(cube.data.max() * dx * dy / self.mass_flux_scaling, dmax)

            assert len(hist_data) == 3

            hist_kwargs = {}
            if self.xlim:
                hist_kwargs['range'] = self.xlim
            else:
                hist_kwargs['range'] = (0, dmax)

            if self.nbins:
                hist_kwargs['bins'] = self.nbins

            #y_min, bin_edges = np.histogram(hist_data[2].data, bins=50, range=(0, dmax))
            #y_max, bin_edges = np.histogram(hist_data[0].data, bins=50, range=(0, dmax))

            y, bin_edges = np.histogram(hist_data[1].data * dx * dy / self.mass_flux_scaling,
                                        **hist_kwargs)
            y_density, bin_edges = np.histogram(hist_data[1].data * dx * dy / self.mass_flux_scaling,
                                                density=True,
                                                **hist_kwargs)
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            width = bin_edges[1:] - bin_edges[:-1]
            y2 = bin_centers * y

            expt_data['y'] = y
            expt_data['y_density'] = y_density
            expt_data['y2'] = y2
            expt_data['bin_centers'] = bin_centers
            expt_data['width'] = width
            expt_data['num_clds'] = len(hist_data[1].data)

    def _calc_fit_linregress(self):
        linregress_details = ['z, expt, m, c, rval, pval, stderr']
        for expt_data in self.all_expt_data:
            logger.debug('Calc fitted linregress on log y for {}, {}',
                         expt_data['height_level_index'],
                         expt_data['expt_obj'].name)
            y = expt_data['y']
            bin_centers = expt_data['bin_centers']

            # TODO: Fix! Dodgy logic here.
            # y[y > 10] selects all values where y > 10
            # :len(log_y) selects values up to len - the two are not nec. the same!
            # Only calc. logarithm where y > 10.
            log_y = np.log(y[y > 10])
            x = bin_centers[:len(log_y)]

            m, c, rval, pval, stderr = stats.linregress(x[1:], log_y[1:])
            linregress_details.append('{},{},{},{},{},{},{}'
                                      .format(expt_data['height_level_index'],
                                              expt_data['expt_obj'].name,
                                              m, c, rval, pval, stderr))
            expt_data['linregress'] = (x, m, c, rval, pval, stderr)
        self.save_text('mf_linregress.csv', '\n'.join(linregress_details) + '\n')

    def _calc_fit_expt(self):
        for expt_data in self.all_expt_data:
            logger.debug('Calc fitted exp for {}, {}',
                         expt_data['height_level_index'],
                         expt_data['expt_obj'].name)
            y_density = expt_data['y_density']
            bin_centers = expt_data['bin_centers']
            width = expt_data['width']

            def exp_dn(x, lmbda):
                return lmbda * np.exp(-lmbda * x)

            beta = (bin_centers * y_density * width).sum()
            lmbda_guess = 1 / beta
            logger.debug('lambda guess: {}', lmbda_guess)
            popt, pcov = optimize.curve_fit(exp_dn, bin_centers, y_density, p0=(lmbda_guess,))
            expt_data['lambda_guess'] = lmbda_guess
            expt_data['lambda_fit'] = popt[0]

    def _plot_mass_flux_hist(self):
        combined_plotters = {}
        for height_level_index in self.height_levels:
            # Set up the plotters that will combine plots of multiple expts on one fig.
            combined_plotters[height_level_index] = {
                'density_mf_contrib': PlotDensityMfContrib(height_level_index,
                                                           self.xlim, self.ylim,
                                                           self.mass_flux_scaling),
                'density_poster': PlotDensityPoster(height_level_index,
                                                    self.xlim, self.ylim),
            }

        for expt_data in self.all_expt_data:
            logger.debug('plot mass flux hist for {}, {}',
                         expt_data['height_level_index'],
                         expt_data['expt_obj'].name)
            height_level_index = expt_data['height_level_index']

            self._plot_indiv_mass_flux(**expt_data)
            # Do multiple plots - expt_data dict deref'ed into keyword args.
            combined_plotters[height_level_index]['density_mf_contrib'].plot(**expt_data)
            combined_plotters[height_level_index]['density_poster'].plot(**expt_data)

        for height_level_index in self.height_levels:
            # Finish up (save and add legends) for multiple plots.
            filepath = self.file_path('density_mf_contrib_z{}'.format(height_level_index))
            combined_plotters[height_level_index]['density_mf_contrib'].finish(filepath)
            filepath = self.file_path('density_poster_z{}'.format(height_level_index))
            combined_plotters[height_level_index]['density_poster'].finish(filepath)

    def _plot_indiv_mass_flux(self, expt_obj, height_level_index, bin_centers, width, y, **kwargs):
        name = '{}.z{}.mass_flux_hist'.format(expt_obj.name, height_level_index)

        # Plot individual mass flux plots.
        for log in [False, True]:
            name = 'log_' + name if log else name
            plt.figure(name)
            plt.bar(bin_centers, y, width=width)
            #plt.bar(bin_centers, y, width=width, yerr=[y - y_min, y_max - y])

            if self.xlim:
                plt.xlim(self.xlim)

            if log:
                if self.ylim:
                    plt.ylim(ymax=self.ylim[1])
                plt.yscale('log')
            else:
                plt.yscale('linear')
                if self.ylim:
                    plt.ylim(self.ylim)
            plt.savefig(self.file_path(name + '.png'))
