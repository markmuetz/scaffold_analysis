"""Run a shard of omnium's expt tasks, using omnium's own task generation and run_task.

Used to produce the reference output_omnium tree in parallel on SLURM. Equivalent to
`omnium run -s production -t expt --all <expt>`, restricted to the given analyses and to
tasks with index % nshards == shard. Must be run from the omnium suite (workspace) dir.

    python omnium_shard_run.py <expt> <analysis1,analysis2,...> <shard> <nshards> [--dry-run]
"""
import os
import socket
import sys
import traceback

import omnium
from omnium.analysis import AnalysisPkgs
from omnium.run_control import RunControl
from omnium.setup_logging import add_file_logging, setup_logger
from omnium.suite import Suite


def main(expt, analyses, shard, nshards, dry_run=False):
    logger = setup_logger(debug=False, colour=False, warn_stderr=True)
    suite = Suite(os.getcwd(), False)
    assert suite.is_in_suite, 'not in an omnium suite'
    log_dir = os.path.join(suite.suite_dir, '.omnium/log')
    add_file_logging(os.path.join(log_dir, '{}_shard{}of{}_{}.log'
                                  .format(expt, shard, nshards, socket.gethostname())))
    suite.set_analysis_pkgs(AnalysisPkgs(['scaffold'], suite))

    run_control = RunControl(suite, 'expt', [expt], 'production')
    run_control.gen_analysis_workflow()
    run_control.gen_tasks()
    all_tasks = run_control._task_master.all_tasks
    unknown = set(analyses) - {t.analysis_name for t in all_tasks}
    if unknown:
        logger.warning('no tasks for analyses: {}', unknown)
    tasks = [t for t in all_tasks
             if t.analysis_name in analyses and t.index % nshards == shard]
    logger.info('omnium {}: {} {} shard {}/{}: {} tasks',
                omnium.__file__, expt, analyses, shard, nshards, len(tasks))

    failed = []
    for task in tasks:
        if dry_run:
            print(task)
            continue
        try:
            run_control.run_task(task)
        except Exception:
            logger.error('FAILED task {}: {}', task.index, task)
            logger.error(traceback.format_exc())
            failed.append(task)

    if failed:
        logger.error('{} of {} tasks failed', len(failed), len(tasks))
        for task in failed:
            logger.error('  {}', task)
        return 1
    return 0


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if a != '--dry-run']
    sys.exit(main(args[0], args[1].split(','), int(args[2]), int(args[3]),
                  dry_run='--dry-run' in sys.argv))
