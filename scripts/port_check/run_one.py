"""Run one scaffold analyser with explicit inputs/outputs (no omnium), for checking the port.

    python run_one.py <analysis_name> <suite_dir> <output> <input> [<input> ...]
        [--expt E] [--expts E1,E2,...]
"""
import argparse

from scaffold import ANALYSERS
from scaffold.analysis_base import run_analyser
from scaffold.config import SETTINGS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('analysis_name')
    parser.add_argument('suite_dir')
    parser.add_argument('outputs', help='comma-separated output filenames')
    parser.add_argument('inputs', nargs='+')
    parser.add_argument('--expt')
    parser.add_argument('--expts')
    args = parser.parse_args()
    run_analyser(ANALYSERS[args.analysis_name], args.suite_dir, args.inputs,
                 args.outputs.split(','), SETTINGS, expt=args.expt,
                 expts=args.expts.split(',') if args.expts else None)


if __name__ == '__main__':
    main()
