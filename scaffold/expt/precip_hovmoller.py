import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import numpy as np
import pylab as plt

from omnium.analyser import Analyser
from omnium.utils import get_cube


class PrecipHovmollerAnalyser(Analyser):
    """Produce a Hovmoller of precipitation."""
    analysis_name = 'precip_hovmoller_analysis'
    single_file = True

    input_dir = 'work/20000101T0000Z/{expt}_atmos'
    input_filename_glob = '{input_dir}/atmos.pp3.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.precip_hovmoller_analysis.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        """Save all results for surf flux analysis."""
        self._plot(self.task.expt)

    def _plot(self, expt):
        precip = get_cube(self.cubes, 4, 203)

        start_time = precip.coord('time').points[0]
        times = precip.coord('time').points - start_time
        timelength_hours = times[-1] - times[0]
        timelength_days = timelength_hours / 24

        # Convert m to km.
        x = precip.coord('grid_longitude').points / 1e3
        y = precip.coord('grid_latitude').points / 1e3
        x_res = x[1] - x[0]
        y_res = y[1] - y[0]
        nx = len(x)
        ny = len(y)

        hov_x = precip.data.mean(axis=1)
        hov_y = precip.data.mean(axis=2)

        plt.clf()
        plt.title('{} Hovmöller (x)'.format(expt))
        # Mask out very small values of precip.
        mask = hov_x < 0.00001
        plt.imshow(np.ma.masked_array(hov_x, mask), interpolation='nearest',
                   extent=[0, x_res * nx, timelength_days, 0], aspect=100, cmap=plt.get_cmap('Blues'),
                   norm=LogNorm(), vmax=2e-3)
        cb = plt.colorbar()
        cb.set_label('precip. (mm hr$^{-1}$)')
        plt.xlabel('x (km)')
        plt.ylabel('time (day)')
        plt.savefig(self.file_path(expt + '_hovmoller_x.png'))
        plt.clf()

        plt.title('{} Hovmöller (y)'.format(expt))
        # Mask out very small values of precip.
        mask = hov_y < 0.00001
        plt.imshow(np.ma.masked_array(hov_y, mask), interpolation='nearest',
                   extent=[0, y_res * ny, timelength_days, 0], aspect=100, cmap=plt.get_cmap('Blues'),
                   norm=LogNorm(), vmax=2e-3)
        cb = plt.colorbar()
        cb.set_label('precip. (mm hr$^{-1}$)')
        plt.xlabel('y (km)')
        plt.ylabel('time (day)')
        plt.savefig(self.file_path(expt + '_hovmoller_y.png'))
