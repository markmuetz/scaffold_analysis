# Sourced by the sbatch scripts. Runs from the omnium reference workspace.
WORKSPACE=/gws/ssde/j25b/afesp/users/mmuetz/scaffold_omnium_remake/u-be530
REPO=/home/users/mmuetz/projects/local/scaffold_analysis
OMNIUM_BIN=/home/users/mmuetz/projects/local/omnium/bin/omnium
export MPLBACKEND=Agg
cd $WORKSPACE
source $WORKSPACE/expts.sh
ALL_EXPTS=($FORCED $RWPs $W0 $W5)
pyrun() { ~/.pixi/bin/pixi run --frozen --manifest-path $REPO/pixi.toml -e omnium python "$@"; }
