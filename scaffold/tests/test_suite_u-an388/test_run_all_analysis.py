import os
import shutil
import omnium


def setup_module():
    suite_dir = os.getenv('SCAFFOLD_SUITE_UAN388_DIR', None)
    assert suite_dir, 'Env set SCAFFOLD_SUITE_UAN388_DIR to your u-an388 test suite'
    assert os.path.exists(suite_dir)
    if os.path.exists(os.path.join(suite_dir, 'omnium_output')):
        shutil.rmtree(os.path.join(suite_dir, 'omnium_output'))

def test_01_run_cycle_analysis():
    suite_dir = os.getenv('SCAFFOLD_SUITE_UAN388_DIR')
    os.chdir(suite_dir)
    omnium.omnium_main(['omnium', '-D', '-b', 'run', '-s', 'test',
                        '-t', 'cycle', '--all', 'S0', 'S4'])

def test_02_run_expt_analysis():
    suite_dir = os.getenv('SCAFFOLD_SUITE_UAN388_DIR')
    os.chdir(suite_dir)
    omnium.omnium_main(['omnium', '-D', '-b', 'run', '-s', 'test',
                        '-t', 'expt', '--all', 'S0', 'S4'])

def test_03_run_suite_analysis():
    suite_dir = os.getenv('SCAFFOLD_SUITE_UAN388_DIR')
    os.chdir(suite_dir)
    omnium.omnium_main(['omnium', '-D', '-b', 'run', '-s', 'test',
                        '-t', 'suite', '--all', 'S0', 'S4'])
