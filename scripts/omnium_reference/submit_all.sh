#!/bin/bash
# Submit the full omnium reference run (all 22 expts) as 5 dependent SLURM stages.
set -euo pipefail
source /home/users/mmuetz/projects/local/scaffold_analysis/scripts/omnium_reference/common.sh
D=$REPO/scripts/omnium_reference
mkdir -p slurm_output
NE=${#ALL_EXPTS[@]}
S1=cloud_analysis,profile_analysis,restart_dump_analysis,surf_flux_analysis,precip_hovmoller_analysis,dump_extract_state,dump_slice_analysis,restart_dump_summary_info
S2=mass_flux_analysis,mass_flux_spatial_scales_analysis,org_analysis,cloud_track_analysis,mse_combined
S3=mass_flux_combined,org_combined,mass_flux_spatial_scales_combined
# summary_correlations takes a hand-made input file (it ignores omnium's file loading).
VERSION_DIR=omnium_output/om_v0.11.3.0_scaffold_v0.6.5.0_2abba69afd
mkdir -p $VERSION_DIR
touch $VERSION_DIR/summary_correlations.dummy $VERSION_DIR/summary_correlations.dummy.done

j1=$(sbatch --parsable -J om_s1 --array=0-$((NE*10-1)) $D/expt_stage.sbatch $S1 10)
j2=$(sbatch --parsable -J om_s2 --dependency=afterany:$j1 --array=0-$((NE*10-1)) $D/expt_stage.sbatch $S2 10)
j3=$(sbatch --parsable -J om_s3 --dependency=afterany:$j2 --time=04:00:00 --mem=16G --array=0-$((NE-1)) $D/expt_stage.sbatch $S3 1)
j4=$(sbatch --parsable -J om_s4 --dependency=afterany:$j3 --array=0-$(($(wc -l < $D/suite_cmds.txt)-1)) $D/suite_stage.sbatch suite_cmds.txt)
j5=$(sbatch --parsable -J om_s5 --dependency=afterany:$j4 --mem=16G --array=0-$(($(wc -l < $D/final_cmds.txt)-1)) $D/suite_stage.sbatch final_cmds.txt)
echo "submitted: s1=$j1 s2=$j2 s3=$j3 s4=$j4 s5=$j5" | tee -a slurm_output/submitted.txt
