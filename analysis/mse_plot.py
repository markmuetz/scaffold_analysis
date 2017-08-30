from collections import OrderedDict

import matplotlib

matplotlib.use('Agg')
import pylab as plt

from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MsePlotter(Analyzer):
    analysis_name = 'mse_plot'

    def run_analysis(self):
        self.expt_mses = OrderedDict()
        for expt in self.expts:
            cubes = self.expt_cubes[expt]

            mse_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mse_profile')
            mse_profile.rename(expt + ': ' + mse_profile.name())
            print(mse_profile.name())
            self.results[expt + '_mse_profile'] = mse_profile
            th = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta')
            # WTF? Why has the name of these coords changed?
            # z = th.coord('level_height').points
            z = th.coord('atmosphere_hybrid_height_coordinate').points
            dz = z[1:] - z[:-1]

            self.expt_mses[expt] = [(msep.data * dz).sum() for msep in mse_profile.slices_over('time')]

    def _plot_mses(self):
        self.append_log('plotting MSEs')
        plt.figure(self.output_filename + '_mse_timeseries')
        for expt in self.expts:
            mses = self.expt_mses[expt]
            plt.title(self.output_filename + '_mse_timeseries')
            plt.plot(mses, label=expt)
        plt.legend()
        plt.savefig(self.figpath('_mse_timeseries.png'))

    def display_results(self):
        self._plot_mses()
        plt.close('all')
