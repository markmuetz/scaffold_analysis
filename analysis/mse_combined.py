from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MseCombinedAnalysis(Analyzer):
    analysis_name = 'mse_combined'
    multi_file = True

    def run_analysis(self):
        cubes = self.cubes

        mse_profile = get_cube_from_attr(cubes, 'omnium_cube_id', 'mse_profile')
        theta = get_cube_from_attr(cubes, 'omnium_cube_id', 'theta')[0] # Just want z from it.

        self.results['mse_profile'] = mse_profile
        self.results['theta'] = theta
