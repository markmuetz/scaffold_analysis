from omnium.analyzer import Analyzer


class DummyAnalyzer(Analyzer):
    analysis_name = 'dummy_analysis'

    def run_analysis(self):
        pass
