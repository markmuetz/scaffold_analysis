from logging import getLogger

import matplotlib

matplotlib.use('Agg')
import numpy as np
import pylab as plt
import iris

from omnium import Analyser
from omnium.utils import get_cube
from omnium.consts import L

from scaffold.colour import EXPT_COLOUR
logger = getLogger('scaf.prof_an')


class SurfFluxPlot(Analyser):
    """Plots surf fluxes and precip timeseries for all expts."""
    analysis_name = 'surf_flux_plot'
    multi_expt = True

    input_dir = 'work/19700101T0000Z/{expt}_atmos'
    input_filename = '{input_dir}/atmos.pp3.nc'
    output_dir = 'omnium_output/{version_dir}/suite_{expts}'
    output_filenames = ['{output_dir}/atmos.surf_flux_plot.dummy']

    def load(self):
        self.load_cubes()

    def run(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')

    def display_results(self):
        """Save all results for surf flux analysis."""
        self.expts_to_plot = self.task.expts
        self._plot()

    def _plot(self):
        lines = ['Expt,PFE [W m-2],LHF [W m-2],SHF [W m-2]']
        for expt in self.task.expts:
            cubes = self.expt_cubes[expt]

            precip = get_cube(cubes, 4, 203)
            lhf = get_cube(cubes, 3, 234)
            shf = get_cube(cubes, 3, 217)

            precip_ts = precip.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
            lhf_ts = lhf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
            shf_ts = shf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

            if expt in self.expts_to_plot:
                start_time = precip_ts.coord('time').points[0]
                times = (precip_ts.coord('time').points - start_time) / 24

                kwargs = {}
                if expt in EXPT_COLOUR:
                    kwargs['color'] = EXPT_COLOUR[expt]
                plot = plt.plot(times, lhf_ts.data, label=expt, linestyle='-', **kwargs)
                #colour = plot[0].get_color()
                colour = plot[0].get_color()
                plt.plot(times, shf_ts.data, color=colour, linestyle='-.')
                # TODO: fix for any time delta.
                # daily smoothing with 15 min means.
                precip_ts_smoothed = np.convolve(precip_ts.data, np.ones((96, )) / 96., mode='same')
                plt.plot(times[96:-96], precip_ts_smoothed[96:-96] * L, color=colour, linestyle='--')

            half_way_index = int(len(precip_ts.data) / 2)
            mean_pfe = precip_ts.data[half_way_index:].mean() * L
            mean_lhf = lhf_ts.data[half_way_index:].mean()
            mean_shf = shf_ts.data[half_way_index:].mean()
            lines.append('{},{},{},{}'.format(expt, mean_pfe, mean_lhf, mean_shf))

        self.save_text('energy_flux.csv', '\n'.join(lines) + '\n')
        plt.ylim((-50, 300))
        plt.xlim((0, 20))
        plt.legend()

        plt.ylabel('flux (W m$^{-2}$)')
        plt.xlabel('time (day)')
        plt.axvline(x=10, linestyle='--', color='k')
        plt.savefig(self.file_path('energy_fluxes.png'))
