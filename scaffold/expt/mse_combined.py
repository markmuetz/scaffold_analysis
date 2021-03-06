from omnium import Analyser
from omnium.utils import get_cube_from_attr


class MseCombinedAnalysis(Analyser):
    analysis_name = 'mse_combined'
    multi_file = True
    input_dir = 'omnium_output/{version_dir}/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.restart_dump_analysis.nc'
    output_dir = 'omnium_output/{version_dir}/{expt}'
    output_filenames = ['{output_dir}/atmos.mse_combined.nc']

    def load(self):
        self.load_cubes()

    def run(self):
        cubes = self.cubes

        mse_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mse_profile')
        theta = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta')

        self.results['mse_profile'] = mse_profile
        self.results['theta'] = theta

    def save(self, state, suite):
        self.save_results_cubes(state, suite)
