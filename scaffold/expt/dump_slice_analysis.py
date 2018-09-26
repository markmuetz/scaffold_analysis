from logging import getLogger

import iris
import numpy as np

from omnium import Analyser
from omnium.consts import Re, L, cp, g
from omnium.utils import get_cube

logger = getLogger('scaf.dump_slice')


class DumpSliceAnalyser(Analyser):
    analysis_name = 'dump_slice_analysis'
    single_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da4??.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.precip_hovmoller_analysis.dummy']
    uses_runid = True
    runid_pattern = 'atmosa_da(?P<runid>\d{3}).nc'

    def load(self):
        self.load_cubes()

    def run(self):
        """Get useful cubes from self.dump, perform sanity chacks and calc MSE, TCW."""
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
        """Save all results for surf flux analysis."""
        self._plot(self.task.expt)

    def _plot(self, expt):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
