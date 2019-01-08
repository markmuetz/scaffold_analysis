from logging import getLogger

import iris

from omnium import Analyser
from omnium.utils import get_cube

logger = getLogger('scaf.dump_extract')


class DumpExtractState(Analyser):
    analysis_name = 'dump_extract_state'
    multi_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmosa_da???.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.dump_extract_state.dummy']
    uses_runid = True
    runid_pattern = 'atmosa_da(?P<runid>\d{3}).nc'

    def load(self):
        import ipdb; ipdb.set_trace()
        filenames = self.task.filenames
        theta = []
        mv = []
        for dump_filename in filenames:
            da = iris.load(dump_filename)
            da_theta = get_cube(da, 0, 4)
            da_mv = get_cube(da, 0, 391)
            theta.append(da_theta)
            mv.append(da_mv)

        self.theta = iris.cube.CubeList.merge(theta)
        self.mv = iris.cube.CubeList.merge(mv)

    def run(self):
        pass

    def display_results(self):
        pass

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
