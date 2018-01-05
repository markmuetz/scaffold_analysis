import omnium
import scaffold


def test_omnium_version():
    omnium.omnium_main(['omnium', 'version'], 'scaffold_testing')

def test_scaffold_analysis_classes():
    assert len(scaffold.analysis_classes) > 0
