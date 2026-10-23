#!/bin/bash
# Rerun the mass flux and thermodynamic summary analyses with the rho-level and theta/Exner fixes,
# in the scaffold_rce_org1_mf_thermo_fix workspace. That workspace is a copy of
# scaffold_rce_org1_centroid_fix (without its per-runid org_analysis files), with the
# centroid-fixed org_combined files linked in, the archive's dumps (runids 240-480) linked into
# share/data/history, and runcontrol_expt restricted to the analyses rerun here.
# Stages: mass_flux_analysis (14 expts x NSHARDS) and restart_dump_summary_info (14)
#         -> mass_flux_combined (14) -> suite mass_flux_plot/summary_correlations.
set -e
NSHARDS=4
WORKSPACE=/gws/ssde/j25b/afesp/users/mmuetz/scaffold_rce_org1_mf_thermo_fix/u-be530
BRANCH=/home/users/mmuetz/projects/local/scaffold_analysis_mf_thermo
MAIN=/home/users/mmuetz/projects/local/scaffold_analysis
PYPATH=$BRANCH:$MAIN/../omnium:$MAIN/../cloud_tracking
PYRUN="$HOME/.pixi/bin/pixi run --frozen --manifest-path $MAIN/pixi.toml -e omnium env PYTHONPATH=$PYPATH MPLBACKEND=Agg python"
SHARD_RUN=$MAIN/scripts/omnium_reference/omnium_shard_run.py
SB="sbatch --parsable --account=afesp --partition=standard --qos=standard"
EXPTS="S0W0Forced S4W0Forced S0W5Forced S4W5Forced RWP_C1 RWP_C2 RWP_C3 RWP_C4 RWP_C5 RWP_C6 RWP_C7 RWP_C8 RWP_C9 RWP_C10"

cd $WORKSPACE
mkdir -p slurm_output

j1=$($SB -J mt_mf_an --time=12:00:00 --mem=16G --array=0-$((14 * NSHARDS - 1)) \
    -o slurm_output/%x_%A_%a.out --wrap "
E=($EXPTS); EXPT=\${E[\$((SLURM_ARRAY_TASK_ID / $NSHARDS))]}; SHARD=\$((SLURM_ARRAY_TASK_ID % $NSHARDS))
cd $WORKSPACE && $PYRUN $SHARD_RUN \$EXPT mass_flux_analysis \$SHARD $NSHARDS")

j2=$($SB -J mt_rdsi --time=04:00:00 --mem=64G --array=0-13 \
    -o slurm_output/%x_%A_%a.out --wrap "
E=($EXPTS); EXPT=\${E[\$SLURM_ARRAY_TASK_ID]}
cd $WORKSPACE && $PYRUN $SHARD_RUN \$EXPT restart_dump_summary_info 0 1")

j3=$($SB -J mt_mf_comb --time=02:00:00 --mem=16G --array=0-13 --dependency=afterok:$j1 \
    -o slurm_output/%x_%A_%a.out --wrap "
E=($EXPTS); EXPT=\${E[\$SLURM_ARRAY_TASK_ID]}
cd $WORKSPACE && $PYRUN $SHARD_RUN \$EXPT mass_flux_combined 0 1")

j4=$($SB -J mt_suite --time=04:00:00 --mem=64G --dependency=afterok:$j2:$j3 \
    -o slurm_output/%x_%j.out --wrap "
set -e; cd $WORKSPACE && source ./expts.sh
for pair in 'mass_flux_plot FORCED' 'mass_flux_plot RWPs' 'summary_correlations FORCED'; do
    read ANALYSIS GROUP <<< \"\$pair\"; echo \"\$(date) \$ANALYSIS on \$GROUP\"
    $PYRUN $MAIN/../omnium/bin/omnium run -s production -t suite --force -a \$ANALYSIS \${!GROUP}
done")

echo "mass_flux_analysis $j1, restart_dump_summary_info $j2, mass_flux_combined $j3, suite $j4"
