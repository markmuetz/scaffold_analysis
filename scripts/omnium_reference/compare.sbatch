#!/bin/bash
#SBATCH --account=afesp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH -o slurm_output/%x_%A_%a.out
# Usage: sbatch --array=0-N -J <name> compare.sbatch <archive_version_dir> <subdirs_file> <report_dir>
# Compares one top-level subdir (line SLURM_ARRAY_TASK_ID+1 of subdirs_file) of the archived
# omnium output against the reference rerun.
source /home/users/mmuetz/projects/local/scaffold_analysis/scripts/omnium_reference/common.sh
A=$1
SUB=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $2)
R=$WORKSPACE/omnium_output/om_v0.11.3.0_scaffold_v0.6.5.0_2abba69afd
mkdir -p $3
pyrun $REPO/scripts/compare_outputs.py $A $R --include "^$SUB/" --report $3/$SUB.csv
