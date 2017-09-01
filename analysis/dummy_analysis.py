from omnium.analyser import Analyser


class DummyAnalyser(Analyser):
    analysis_name = 'dummy_analysis'

    def run_analysis(self):
        pass
