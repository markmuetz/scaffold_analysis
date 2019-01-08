from logging import getLogger
import configparser as cp

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
            if self.suite.check_filename_missing(dump_filename):
                logger.debug('filename {} missing, skipping', dump_filename)
                continue
            da = iris.load(dump_filename)
            da_theta = get_cube(da, 0, 4)
            da_mv = get_cube(da, 0, 391)
            theta.append(da_theta)
            mv.append(da_mv)

        self.theta = iris.cube.CubeList.merge(theta)[0]
        self.mv = iris.cube.CubeList.merge(mv)[0]
        self.z = self.theta.coord('atmosphere_hybrid_height_coordinate').points

    def run(self):
        self.theta_profile = self.theta.data.mean(axis=(0, 2, 3))
        self.mv_profile = self.mv.data.mean(axis=(0, 2, 3))

    def display_results(self):
        opt = cp.ConfigParser()
        opt.add_section('namelist:idealised')

        opt.set('namelist:idealised', 'mv_relax_data', ','.join(['{:.5f}'.format(v) for v in self.mv_profile]))
        opt.set('namelist:idealised', 'mv_relax_height', ','.join(['{:.2f}'.format(v) for v in self.z]))
        opt.set('namelist:idealised', 'mv_relax_time', '0')
        opt.set('namelist:idealised', 'mv_relax_timescale', '25600.0')

        opt.set('namelist:idealised', 'num_mv_relax_heights', str(len(self.z)))
        opt.set('namelist:idealised', 'num_mv_relax_times', '1')
        opt.set('namelist:idealised', 'num_theta_relax_heights', str(len(self.z)))
        opt.set('namelist:idealised', 'num_theta_relax_times', '1')

        opt.set('namelist:idealised', 'theta_relax_data', ','.join(['{:.2f}'.format(v) for v in self.theta_profile]))
        opt.set('namelist:idealised', 'theta_relax_height', ','.join(['{:.2f}'.format(v) for v in self.z]))
        opt.set('namelist:idealised', 'theta_relax_time', '0')
        opt.set('namelist:idealised', 'theta_relax_timescale', '25600.0')

        opt.write(sys.stdout)
        with open('rose-app-S0_relax.conf', 'w') as f:
            opt.write(f)

    def save(self, state, suite):
        with open(self.task.output_filenames[0], 'w') as f:
            f.write('done')
