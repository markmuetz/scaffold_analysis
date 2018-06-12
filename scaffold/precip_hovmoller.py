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

    def _plot(self):
        precip = self.precip

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
        plt.title('{} Hovmöller (x)'.format(self.expt))
        # Mask out very small values of precip.
        mask = hov_x < 0.00001
        plt.imshow(np.ma.masked_array(hov_x, mask), interpolation='nearest',
                   extent=[0, x_res * nx, timelength_days, 0], aspect=100, cmap=plt.get_cmap('Blues'),
                   norm=LogNorm(), vmax=2e-3)
        cb = plt.colorbar()
        cb.set_label('precip. (mm hr$^{-1}$)')
        plt.xlabel('x (km)')
        plt.ylabel('time (day)')
        plt.savefig(self.figpath('hovmoller_x.png'))
        plt.clf()

        plt.title('{} Hovmöller (y)'.format(self.expt))
        # Mask out very small values of precip.
        mask = hov_y < 0.00001
        plt.imshow(np.ma.masked_array(hov_y, mask), interpolation='nearest',
                   extent=[0, y_res * ny, timelength_days, 0], aspect=100, cmap=plt.get_cmap('Blues'),
                   norm=LogNorm(), vmax=2e-3)
        cb = plt.colorbar()
        cb.set_label('precip. (mm hr$^{-1}$)')
        plt.xlabel('y (km)')
        plt.ylabel('time (day)')
        plt.savefig(self.figpath('hovmoller_y.png'))

    def run_analysis(self):
        cubes = self.cubes

        self.precip = get_cube(cubes, 4, 203)

    def display_results(self):
        """Save all results for surf flux analysis."""
        self._plot()
