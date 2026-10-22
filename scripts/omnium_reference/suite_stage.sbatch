#!/bin/bash
#SBATCH --account=afesp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH -o slurm_output/%x_%A_%a.out
# Usage: sbatch --array=0-N -J <stage> suite_stage.sbatch <cmd_file>
# Each line of cmd_file is: <analysis> <expt group var name>, e.g. "precip_plot FORCED".
source /home/users/mmuetz/projects/local/scaffold_analysis/scripts/omnium_reference/common.sh
read ANALYSIS GROUP <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $REPO/scripts/omnium_reference/$1)"
EXPTS=${!GROUP}
echo "$(date) $SLURM_JOB_NAME: $ANALYSIS on $GROUP ($EXPTS)"
pyrun $OMNIUM_BIN run -s production -t suite -a $ANALYSIS $EXPTS
