from logging import getLogger

import numpy as np
import scipy
import matplotlib.pyplot as plt

from omnium import Analyser
from omnium.consts import Re, L, cp, g
from omnium.utils import get_cube

from scaffold.vertlev import VertLev

logger = getLogger('scaf.dump_slice')


class DumpSliceAnalyser(Analyser):
    analysis_name = 'dump_slice_analysis'
    single_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da4??.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.dump_slice_analysis.dummy']
    uses_runid = True
    runid_pattern = 'atmosa_da(?P<runid>\d{3}).nc'

    def load(self):
        self.load_cubes()

    def run(self):
        dump = self.cubes
        self.rho = get_cube(dump, 0, 253) / Re ** 2
        self.rho_d = get_cube(dump, 0, 389)

        self.th = get_cube(dump, 0, 4)
        self.ep = get_cube(dump, 0, 255)

        self.q = get_cube(dump, 0, 10)
        self.qcl = get_cube(dump, 0, 254)
        self.qcf = get_cube(dump, 0, 12)
        self.qrain = get_cube(dump, 0, 272)
        self.qgraup = get_cube(dump, 0, 273)

    def display_results(self):
        self._plot(self.task.expt)

    def _plot(self, expt):
        vertlevs = VertLev(self.suite.suite_dir)
        fig, ax = plt.subplots(dpi=100)
        data = self.qcl.data
        # Coords are model_level, y, x or model_level, lat, lon
        data_mean = data.mean(axis=1)
        Nx = data.shape[2]
        data_rbs = scipy.interpolate.RectBivariateSpline(vertlevs.z_theta, np.arange(Nx),
                                                         data_mean)
        data_interp = data_rbs(np.linspace(0, 40000, 400), np.linspace(0, Nx - 1, Nx))
        # Only go up to 20 km and use aspect ratio to plot equal aspect (allowing for diff in coords).
        im = ax.imshow(data_interp[:200], origin='lower', cmap='Blues', aspect=0.1)

        ax.set_title('')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('height (100 m)')
        plt.colorbar(im)
        plt.savefig(self.file_path(expt + '_slice_mean_over_y.png'))
        plt.close('all')

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
