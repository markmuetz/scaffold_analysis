#!/bin/bash
#SBATCH --account=afesp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH -o slurm_output/%x_%A_%a.out
# Usage: sbatch --array=0-$((22*NSHARDS-1)) -J <stage> expt_stage.sbatch <analyses,comma,sep> <nshards>
ANALYSES=$1
NSHARDS=$2
source /home/users/mmuetz/projects/local/scaffold_analysis/scripts/omnium_reference/common.sh
EXPT=${ALL_EXPTS[$((SLURM_ARRAY_TASK_ID / NSHARDS))]}
SHARD=$((SLURM_ARRAY_TASK_ID % NSHARDS))
echo "$(date) $SLURM_JOB_NAME $EXPT shard $SHARD/$NSHARDS: $ANALYSES"
pyrun $REPO/scripts/omnium_reference/omnium_shard_run.py $EXPT $ANALYSES $SHARD $NSHARDS
