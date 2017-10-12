from omnium.analyser import Analyser
from omnium.utils import get_cube_from_attr


class MassFluxTimeSeriesAnalyser(Analyser):
    analysis_name = 'mass_flux_timeseries_analysis'
    single_file = True

    def run_analysis(self):
        cubes = self.cubes

        w = get_cube_from_attr(cubes, 'id', 'w')

        self.results['w'] = w
