"""Minimal, omnium-free replacement for the omnium Analyser framework.

The scaffold analysers keep their omnium-style interface (`self.settings`, `self.task`,
`self.suite`, `load_cubes`, `save_results_cubes`, `file_path`, ...), but no longer depend on
omnium. Instead of omnium's task discovery (globs, `.done` files, version_dir hashing), the
caller (a remake rule) gives each analyser its exact inputs and outputs and calls `run_all`.
"""
import os
import re
from collections import OrderedDict
from logging import getLogger
from types import SimpleNamespace

import iris

from scaffold.cube_utils import ScaffoldError
from scaffold.um_config import UMConfig

logger = getLogger('scaffold.analysis')


class Settings:
    """Attribute access to a plain settings dict (replaces omnium AnalysisSettings)."""
    def __init__(self, settings):
        self._settings = dict(settings)

    def __getattr__(self, item):
        if item.startswith('__'):
            raise AttributeError(item)
        try:
            return self._settings[item]
        except KeyError:
            raise AttributeError(item)

    def has_item(self, item):
        return item in self._settings


class Suite:
    """The bits of an omnium Suite that the analysers use."""
    def __init__(self, suite_dir, expt_dataw_dir='work/19700101T0000Z/{expt}_atmos'):
        self.suite_dir = str(suite_dir)
        self.expt_dataw_dir = expt_dataw_dir

    def expt_config(self, expt):
        conf = os.path.join(self.suite_dir, self.expt_dataw_dir.format(expt=expt),
                            'rose-app-run.conf')
        config = UMConfig(conf)
        config.name = expt
        return config

    @staticmethod
    def check_filename_missing(filename):
        # omnium symlinked missing files to a placeholder; inputs are now declared explicitly.
        return False


class ExptList(list):
    """omnium ExptList lookalike: `find(names)` then `get(name)` returns the expt's UMConfig."""
    def __init__(self, suite):
        super().__init__()
        self._suite = suite

    def find(self, expt_names):
        for expt_name in expt_names:
            self.append(self._suite.expt_config(expt_name))

    def get(self, expt_name):
        for expt in self:
            if expt.name == expt_name:
                return expt
        raise ScaffoldError('No expt with name {}'.format(expt_name))


class Analyser:
    analysis_name = None
    # One of these must be overridden by the subclass.
    single_file = False
    multi_file = False
    multi_expt = False

    # Kept for reference: these described omnium's file discovery. Inputs and outputs are now
    # passed in explicitly, but some analysers still use `runid_pattern`/`min_runid`.
    input_dir = None
    input_filename_glob = None
    input_filename = None
    input_filenames = None
    output_dir = None
    output_filenames = None

    uses_runid = False
    min_runid = 0
    max_runid = int(1e10)
    runid_pattern = None

    force = False

    @classmethod
    def get_runid_filename_vars(cls, filename):
        assert cls.uses_runid and cls.runid_pattern
        match = re.match(cls.runid_pattern, os.path.basename(str(filename)))
        if match:
            return int(match.group('runid')), match.groupdict()
        raise ScaffoldError('Could not find runid in {} using {}'.format(filename, cls.runid_pattern))

    def __init__(self, suite, task, settings):
        assert sum([self.single_file, self.multi_file, self.multi_expt]) == 1
        self.suite = suite
        self.task = task
        self.settings = settings
        self.output_filename = task.output_filenames[0]
        self.results = OrderedDict()
        for output_dir in {os.path.dirname(fn) for fn in task.output_filenames}:
            os.makedirs(output_dir, exist_ok=True)

    def append_log(self, message):
        logger.debug('{}: {}'.format(self.analysis_name, message))

    def load_cubes(self):
        if self.multi_file:
            cubes = iris.load(self.task.filenames)
            self.cubes = iris.cube.CubeList.concatenate(cubes)
        elif self.multi_expt:
            self.expt_cubes = OrderedDict()
            for expt, filename in zip(self.task.expts, self.task.filenames):
                self.expt_cubes[expt] = iris.load(filename)
        else:
            self.cubes = iris.load(self.task.filenames[0])

    def save_results_cubes(self, state=None, suite=None):
        for cube_id, cube in self.results.items():
            cube.attributes['omnium_cube_id'] = cube_id
            cube.attributes['omnium_process'] = self.analysis_name
        cubelist = iris.cube.CubeList(list(self.results.values()))
        if not len(cubelist):
            logger.warning('No results to save')
            return
        save_kwargs = {'zlib': True}
        if self.settings.has_item('zlib'):
            save_kwargs['zlib'] = self.settings.zlib
        if self.settings.has_item('complevel'):
            save_kwargs['complevel'] = self.settings.complevel
        if os.path.islink(self.output_filename):
            os.remove(self.output_filename)
        iris.save(cubelist, self.output_filename, **save_kwargs)

    def file_path(self, name):
        """Path for an extra output next to the task's main output."""
        filename = os.path.join(os.path.dirname(self.output_filename),
                                self.analysis_name + '_' + name)
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        return filename

    def save_text(self, name, text):
        with open(self.file_path(name), 'w') as f:
            f.write(text)

    def run_all(self):
        """load -> run -> display_results (if any) -> save; checks declared outputs exist."""
        self.load()
        self.run()
        if hasattr(self, 'display_results'):
            self.display_results()
        self.save(None, self.suite)
        missing = [fn for fn in self.task.output_filenames if not os.path.exists(fn)]
        if missing:
            raise ScaffoldError('{} did not produce: {}'.format(self.analysis_name, missing))


def run_analyser(analyser_cls, suite_dir, inputs, outputs, settings, expt=None, expts=None,
                 runid=None):
    """Run one analyser task with explicit inputs/outputs.

    inputs: list of input filenames (for multi_expt, one per expt, in `expts` order).
    outputs: list of output filenames (the analyser's declared outputs).
    """
    inputs = [str(fn) for fn in inputs]
    if analyser_cls.uses_runid and runid is None and analyser_cls.single_file:
        runid, _ = analyser_cls.get_runid_filename_vars(inputs[0])
    task = SimpleNamespace(expt=expt if not analyser_cls.multi_expt else expts,
                           expts=expts, runid=runid, filenames=inputs,
                           output_filenames=[str(fn) for fn in outputs])
    analyser = analyser_cls(Suite(suite_dir), task, Settings(settings))
    analyser.run_all()
    return analyser
