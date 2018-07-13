from collections import OrderedDict

import matplotlib

matplotlib.use('Agg')
import pylab as plt

from omnium import Analyser
from omnium.utils import get_cube_from_attr


class MsePlotter(Analyser):
    """Plots timeseries of MSE."""
    analysis_name = 'mse_plot'
    multi_expt = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.mse_combined.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.mse_plot.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        self.expt_mses = OrderedDict()
        for expt in self.task.expts:
            cubes = self.expt_cubes[expt]

            mse_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mse_profile')
            mse_profile.rename(expt + ': ' + mse_profile.name())
            print(mse_profile.name())
            self.results[expt + '_mse_profile'] = mse_profile
            th = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta')
            # Old:
            # This is from fields file.
            # z = th.coord('level_height').points
            # This is after iris conv to nc.
            # z = th.coord('atmosphere_hybrid_height_coordinate').points
            # TODO: check this is right.
            z = th.coord('altitude').points[:, 0, 0]
            dz = z[1:] - z[:-1]
            self.expt_mses[expt] = [(msep.data * dz).sum() for msep in mse_profile.slices_over('time')]

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        self._plot_mses()
        plt.close('all')

    def _plot_mses(self):
        self.append_log('plotting MSEs')
        plt.figure(self.output_filename + '_mse_timeseries')
        for expt in self.task.expts:
            mses = self.expt_mses[expt]
            plt.title(self.output_filename + '_mse_timeseries')
            plt.plot(mses, label=expt)
        plt.legend()
        plt.savefig(self.file_path('_mse_timeseries.png'))
