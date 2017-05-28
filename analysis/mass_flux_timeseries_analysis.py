from omnium.analyzer import Analyzer
from omnium.utils import get_cube_from_attr


class MassFluxTimeSeriesAnalyzer(Analyzer):
    analysis_name = 'mass_flux_timeseries_analysis'

    def run_analysis(self):
        cubes = self.cubes

        w = get_cube_from_attr(cubes, 'id', 'w')

	self.results['w'] = w
