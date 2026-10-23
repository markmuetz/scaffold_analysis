#!/bin/bash
#SBATCH --account=afesp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH -o %x_%A_%a.out
# Compare one top-level subdir of two output trees (one array element per subdir).
# Usage: sbatch --array=0-N -J <name> compare_trees.sbatch <tree_a> <tree_b> <subdirs_file> <report_dir>
P=/home/users/mmuetz/projects/local/scaffold_analysis_port
SUB=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $3)
mkdir -p $4
cd /tmp && ~/.pixi/bin/pixi run --frozen --manifest-path $P/../scaffold_analysis/pixi.toml -e remake \
    python $P/scripts/compare_outputs.py $1 $2 --include "^$SUB/" --report $4/$SUB.csv
