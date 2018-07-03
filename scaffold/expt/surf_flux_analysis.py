import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
import iris

from omnium.analyser import Analyser
from omnium.utils import get_cube
from omnium.consts import L

from scaffold.scaffold_settings import settings


class SurfFluxAnalyser(Analyser):
    """Analyse surface fluxes, plot graphs of energy/moisture fluxes."""
    analysis_name = 'surf_flux_analysis'
    single_file = True
    input_dir = 'work/20000101T0000Z/{expt}_atmos'
    input_filename = '{input_dir}/atmos.pp3.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.surf_flux_analysis.nc']

    settings = settings

    def load(self):
        self.load_cubes()

    def save(self, state, suite):
        self.save_results_cubes(state, suite)

    def run(self):
        cubes = self.cubes

        precip = get_cube(cubes, 4, 203)
        lhf = get_cube(cubes, 3, 234)
        shf = get_cube(cubes, 3, 217)

        self.precip_ts = precip.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.lhf_ts = lhf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        self.shf_ts = shf.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)

        self.results['precip_ts'] = self.precip_ts
        self.results['lhf_ts'] = self.lhf_ts
        self.results['shf_ts'] = self.shf_ts

    def display_results(self):
        """Save all results for surf flux analysis."""
        self._plot()

    def _plot(self):
        precip_ts = self.results['precip_ts']
        lhf_ts = self.results['lhf_ts']
        shf_ts = self.results['shf_ts']

        start_time = precip_ts.coord('time').points[0]
        times = precip_ts.coord('time').points - start_time

        plt.figure(self.output_filename + '_energy_fluxes')
        plt.clf()
        plt.title(self.output_filename + '_energy_fluxes')

        plt.plot(times, lhf_ts.data, 'g-', label='LHF')
        plt.plot(times, shf_ts.data, 'r-', label='SHF')
        # TODO: fix for any time delta.
        # daily smoothing with 15 min means.
        precip_ts_smoothed = np.convolve(precip_ts.data, np.ones((96, )) / 96., mode='same')
        plt.plot(times[96:-96], precip_ts_smoothed[96:-96] * L, 'b-', label='PFE (smoothed)')
        plt.ylim((0, 400))
        plt.legend()
        plt.ylabel('flux (W m$^{-2}$)')
        plt.xlabel('time (hrs)')
        plt.savefig(self.file_path('energy_fluxes.png'))

        plt.figure(self.output_filename + '_water_fluxes')
        plt.clf()
        plt.title(self.output_filename + '_water_fluxes')
        plt.plot(times, lhf_ts.data / L * 3600, 'g-', label='UWF')
        plt.plot(times[96:-96], precip_ts_smoothed[96:-96] * 3600, 'b-', label='Precip. (smoothed)')
        plt.plot(times, precip_ts.data * 3600, 'b--', label='Precip.')
        plt.ylim((0, 1))
        plt.legend()
        plt.ylabel('water flux (mm hr$^{-1}$)')
        plt.xlabel('time (hrs)')
        plt.savefig(self.file_path('water_fluxes.png'))
