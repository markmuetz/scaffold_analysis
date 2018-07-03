from mock import Mock, patch

import omnium
from omnium.setup_logging import setup_logger
import scaffold


def test_omnium_version():
    omnium.omnium_main(['omnium', 'version'], 'scaffold_testing')

@patch('os.path.exists')
@patch('os.makedirs')
def test_scaffold_analysis_classes(mock_makedirs, mock_exists):
    assert len(scaffold.analysis_classes) > 0
    setup_logger(debug=True, colour=False)
    mock_exists.return_value = False
    suite = Mock()
    task = Mock()
    suite.check_filename_missing.return_value = False
    task.filenames = ['/a/b/c.txt']
    task.output_filenames = ['/a/b/c.out']
    for analysis_class in scaffold.analysis_classes:
        analyser = analysis_class(suite, task)
        assert analyser
