#!/bin/bash
#SBATCH --account=afesp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH -o slurm_output/%x_%j.out
# Run omnium suite analyses for the rce_org1 figures, using the rce_org1_figs worktree.
# Usage: sbatch -J <name> run_suite.sbatch "<analysis> <GROUP>" ["<analysis> <GROUP>" ...]
# Runs each pair in order (later ones may read earlier outputs).
set -e
FIGS=/home/users/mmuetz/projects/local/scaffold_analysis_figs
MAIN=/home/users/mmuetz/projects/local/scaffold_analysis
cd /gws/ssde/j25b/afesp/users/mmuetz/scaffold_rce_org1_figs/u-be530
source ./expts.sh
export MPLBACKEND=Agg
for pair in "$@"; do
    read ANALYSIS GROUP <<< "$pair"
    echo "$(date) $ANALYSIS on $GROUP"
    ~/.pixi/bin/pixi run --frozen --manifest-path $MAIN/pixi.toml -e omnium \
        env PYTHONPATH=$FIGS:$MAIN/../omnium:$MAIN/../cloud_tracking \
        python $MAIN/../omnium/bin/omnium run -s production -t suite --force -a $ANALYSIS ${!GROUP}
done
