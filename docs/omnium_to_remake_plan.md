# Plan: port scaffold_analysis from omnium to remake

Branch: `omnium_remake_v0.8.5` (remake v0.8.5+, `../remake` at `04462f3`).

Goal: rewrite every omnium workflow in `scaffold/` as remake rules, and prove the port is
correct by producing two output trees from the same u-be530 inputs, `output_omnium` and
`output_remake`, and diffing them.

Data (read-only): `/gws/ssde/j25b/afesp/users/mmuetz/phd/rdf_mirror/um11.0_runs/archive/u-be530`.

## 1. What omnium does today

omnium's job here is the bookkeeping around ~33 `Analyser` subclasses in `scaffold/`:

| omnium concept | what it does | used by |
|---|---|---|
| `runcontrol_{cycle,expt,suite}` in `app/omnium/rose-app.conf` | ordered list of enabled analysers per run type | the whole workflow |
| `single_file` + `uses_runid` + `runid_pattern` | one task per input file; `runid` parsed from the filename with a regex | cloud, mass_flux, mass_flux_spatial_scales, org, profile, restart_dump, dump_slice, converters |
| `multi_file` | one task per expt over a glob of files | *_combined, cloud_track, restart_dump_summary_info, dump_extract_state |
| `multi_expt` | one task over the expts given on the command line, output in `suite_{expt1_expt2_...}` | all `scaffold/suite/*` plots |
| "virtual drive" | globs input dirs, then adds each task's planned outputs, so `input_filename_glob` can match files that don't exist yet; a file only counts if `<file>.done` exists | DAG construction |
| `TaskMaster` | links tasks by matching output to input filenames | DAG construction |
| `AnalysisSettings` + `version_dir` | settings dict → JSON → sha1; output dir is `om_v<omnium>_<pkg>_v<ver>_<sha1[:10]>` | every output path |
| `.done` files, `--no-run-if-started`, `force` | rerun control | every task |
| `production` mode | refuses to run unless omnium and scaffold git trees are clean | provenance |
| `save_results_cubes` | saves `self.results` with `omnium_cube_id`, version and git attrs on each cube | most expt analysers |
| `file_path(name)` / `save_text` | writes extra outputs (figures, csv, txt) next to the declared output | all plots and several analysers |
| `suite.expts.get(expt).dx/lx` | parses UM namelists from `work/<cycle>/<expt>_atmos/rose-app-run.conf` | org, cloud_track, mass_flux_spatial_scales |
| `VertLev(suite_dir)` | reads `app/um/file/rce_vertlevs.nml` | mass_flux, profile, dump_* |
| MPI master/worker | parallel `cycle` runs | converters only |

## 2. Target design in remake

- **Remakefile(s)**: `remakefiles/rce_org1.py` (split into `expt.py` and `suite.py` using
  the `ex5_multifile` pattern if it gets large). Rules call plain functions in `scaffold/`; no
  `Analyser` base class.
- **Config, replacing `AnalysisSettings`**: a plain module `scaffold/config.py` with:
  - `SUITE_DIR` (input root) and `OUTPUT_DIR` (`.../output_remake`);
  - expt groups (`FORCED`, `RWPS`, `W0`, `W5`, from `expts.sh`) and their runid ranges (forced
    and RWP 0–596, rel2 0–476);
  - `SETTINGS`, an ordinary dict with the production values.

  No JSON, no sha1, no `version_dir`, no `settings.json`/`.name` files. Each rule declares only
  the settings it reads with `uses=`, so changing `w_thresh` reruns only the rules that read
  it. The `production`/`test` switch becomes a query (`-Q`) or a `TEST` flag that shrinks the
  matrix.
- **Output layout**: `output_remake/<expt>/atmos.<runid>.<analysis>.nc` and
  `output_remake/suite_<group>/...`. The same relative layout is used for the omnium reference
  once the version-hash directory level is stripped, so the two trees map 1:1.
- **Mapping**:
  - `single_file`+runid → `matrix=[{'expt': e, 'runid': r} ...]`. This is a `list[dict]`
    because rel2 expts have fewer runids. Runids come from known ranges, not globs.
  - `multi_file` → fan-in rule with `matrix={'expt': ...}` and a callable `inputs(expt)` over
    the runids.
  - `multi_expt` → `matrix={'group': [...]}` with `inputs(group)` over the group's expts.
  - `min_runid`/`max_runid` and `settings.start_runid/end_runid` → filter the runid lists
    when building inputs, so each task declares only what it reads.
  - `.done`, `already_analysed`, `force` → remake's DB (`remake run`, `set-state`, `why`).
  - MPI / the old `slurm_array_run_omnium.sh` → `remake run -E slurm` (LOTUS2) or `-E multiproc`.
  - `Analyser` helpers → `scaffold/io.py`: `load_cubes`, `save_cubes(results, path)`, which
    keeps the `omnium_cube_id` attribute because downstream code selects cubes by it, but
    drops the `omnium_vn`/git attributes. `omnium.utils.get_cube*` is vendored into
    `scaffold/utils.py`.
  - `suite.expts.get(expt).dx` → `scaffold/um_config.py`, which reads `rose-app-run.conf`.
    The conf and `rce_vertlevs.nml` become declared inputs.
- **Converters** (`cycle_converter`, `cycle_dump_converter`, `expt_converter`): not ported. All
  inputs already exist as `.nc`. Keep a one-off script only if the raw FF files ever need
  re-converting.

## 3. Validation: `output_omnium` vs `output_remake`

1. **Reference env** (`pixi -e omnium`): re-creates the old `omnium/envs/omnium_env.yml`
   (python 3.7.6, iris 2.4.0, numpy 1.18.1, scipy 1.4.1, matplotlib 3.2.1, metpy 0.12.0,
   pandas 1.0.3). omnium, cloud_tracking and scaffold come from PYTHONPATH (sibling checkouts).
2. **Reference workspace** on the afesp GWS (not in the archive): a suite-shaped directory
   with `rose-suite.info`, `.omnium/suite.conf`, `app/omnium/rose-app.conf`,
   `app/um/file/rce_vertlevs.nml`, and file-level symlinks to the archive's
   `share/data/history/<expt>/*` and `work/*/<expt>_atmos/*`. omnium writes
   `omnium_output/<version_dir>/`, which becomes `output_omnium`. The archive is never written.
3. **Sanity check**: diff the fresh `output_omnium` against the archived
   `omnium_output/om_v0.11.2.0_scaffold_v0.6.5.0_2abba69afd`. Any difference here is caused by
   the env or code drift since 2019, not by the port.
4. **Port in two steps, to keep diffs attributable**:
   - (a) Refactor the scaffold science into plain functions and run them in the **omnium env**
     through a tiny driver (no omnium, no remake). The result must be **bit-identical** to
     `output_omnium`.
   - (b) Run the remake rules in the **remake env** (python ≥3.10, iris 3, numpy 2, current
     matplotlib/metpy) to produce `output_remake`. Remaining diffs come from the env alone,
     and each one gets explained and signed off.
5. **`scripts/compare_outputs.py`**, by output type:
   - `.nc`: per cube/variable values exact, then `allclose`; coords and units; ignore
     `omnium_*`/git attributes.
   - `.hdf`/`.csv`: `pandas.testing.assert_frame_equal`.
   - `.txt`: text diff.
   - `.pkl`: unpickle, then compare stats dicts and tracker summaries.
   - `.png`/`.pdf`: list files present/missing; pixel diffs only as a flag (matplotlib
     versions differ). The comparison relies on the data files behind the figures.
6. **Order of work**: smoke test on one expt with a handful of runids, then all 14 forced+RWP
   expts (the paper set), then the rel2 expts and the remaining suite plots.

## 4. What I need to learn

- **remake** (done for basics): rules, `matrix` forms incl. `list[dict]`, callable inputs, fan-in,
  `uses` one level deep, `check_outputs`/`set-state --success` for adoption,
  `-E slurm` + per-rule `config={'slurm': ...}`, `remake why/info/ls-tasks`, the
  `remake-slurm` monitoring workflow. Still to read: `docs/guide/slurm.md` (job env on LOTUS2),
  `debugging.md`, `ex5_multifile`, and how imported-module code is treated by change detection.
- **Each scaffold analyser in detail** (6.6k lines): which inputs it really reads (vs the
  glob), every side output it writes, and any hidden state.
- **cloud_tracking** package: `Tracker`, `generate_stats`, `output_stats_to_file`, what is pickled.
- **iris 2.4 → 3.x** migration (load/concatenate behaviour, netCDF save attributes, coord
  APIs), **metpy 0.12 → 1.x** (units/function signatures), numpy 1.18 → 2.x deprecations,
  pandas 1.0 → 2.x (`to_hdf`).
- **JASMIN LOTUS2**: partitions, QoS and accounts. The old `short-serial-4hr` script is LOTUS1.

## 5. Pain points and risks

1. **Python split**: remake needs Python ≥3.10, while the reference env is 3.7/iris 2.4. One env
   can't run both, hence the two-step validation in §3.4.
2. **Hidden outputs**: plots declare a `.dummy` sentinel and write the real figures/csv through
   `file_path`. remake should declare real outputs where practical. Otherwise a manifest output
   listing what was written.
3. **Hidden state**: `precip_plot` writes `*_cached_hov_{x,y}.nc` and reads them back on later
   runs. This becomes an explicit upstream rule.
4. **Undeclared cross-rule reads**: `energy_balance` reads csvs that are side outputs of
   `surf_flux_plot` and `relaxation_plot`; `summary_correlations` reads
   `omnium_output/{version_dir}/summary_correlations.dummy` and gathers across suite dirs. These
   need real edges.
5. **Inputs from globs + `.done`** become static lists. The irregular selections need care:
   `restart_dump_analysis` uses `atmosa_da4??`, `dump_slice` uses `da4[58][60]`, `profile_plot`
   uses runid 456 only, `cloud_track` uses `min_runid=480`, `dump_extract_state` uses
   `min_runid=264`. Rel2 runs are shorter.
6. **Settings read in two places**: some runid limits are class attributes, others come from
   settings inside `run()`, after all files were "inputs". Unify them in config so each task
   declares exactly what it reads.
7. **Change detection**: remake hashes the rule body plus `uses` one level deep. Edits inside
   imported `scaffold/*.py` functions won't trigger reruns unless the function is listed in
   `uses` (or the task is reset with `set-state --pending`). Decide on a convention: list the
   science function in `uses`.
8. **Environment-induced diffs** (iris save conventions, float ordering in scipy/numpy labelling,
   metpy units, matplotlib rendering): expected in step 3.4(b). Each one needs explaining, not
   just tolerating.
9. **Pickles**: cloud-track outputs pickle `cloud_tracking` objects, so the comparison must load
   them with matching class definitions.
10. **Cost**: cloud_analysis etc. over 150 files × 14 expts of 256×256×71 fields (48 steps per
    file after day 20) needs SLURM. Both trees have to be produced, so it costs twice. The
    smoke test comes first.
11. **Archive safety**: the reference workspace only symlinks to the archive. omnium's
    converters would write next to their inputs, so they are disabled in the reference
    `rose-app.conf` (inputs are already `.nc` with `.done`).
12. **netCDF variable names**: iris/netCDF versions differ in sanitising variable names
    (`w1.0` vs `w1_0`). Compare by normalised name. Downstream code selects cubes by
    `omnium_cube_id`, not by variable name, so the analysis itself is unaffected.
13. **Paper reproducibility**: figure names/paths used by the rce_org1 paper must be kept or
    mapped (`org_plot` figure names were fixed in `c02d325`).

## 6. Status

- [x] Branch created; plan written.
- [x] `pixi -e omnium` env installed (`pixi.toml`/`pixi.lock`): conda-forge still solves
      py3.7.6 + iris 2.4.0. cartopy has to be pinned to 0.17.0, because metpy 0.12 breaks on
      newer cartopy.
- [x] Reference workspace: `/gws/ssde/j25b/afesp/users/mmuetz/scaffold_omnium_remake/u-be530`.
      It holds copies of `rose-suite.info`, `.omnium/suite.conf`, `rce_vertlevs.nml` and
      `expts.sh`, plus file-level symlinks into the archive (16,536 links). Its
      `app/omnium/rose-app.conf` has the converters removed and `commit=omnium_remake_v0.8.5`,
      and `.omnium/analysis_pkgs/scaffold_analysis` links to this repo. omnium writes only
      under `<workspace>/omnium_output/`.
- [x] omnium runs: `-t expt --all -o S0W0Forced` plans 640 tasks. Smoke test with one
      `cloud_analysis` on `S0W0Forced/atmos.240.pp1.nc`: 2.6 min, 1.7 GB peak RSS. All 24
      variables are bit-identical to the archived `om_v0.11.3.0_scaffold_v0.6.5.0_2abba69afd`
      output, and the settings hash is the same (`2abba69afd`).
      Caveat: the archived files keep `.`/`-` in netCDF variable names
      (`labelled_clouds_z16_w1.0_qcl5e-06`), while this env writes `_`
      (`..._w1_0_qcl5e_06`), so comparisons match names after `re.sub(r'[.\-]', '_', name)`.
- [~] Full `output_omnium` for all 22 expts submitted on 2026-10-22 (SLURM jobs
      61547971–75, `--account=afesp --partition=standard --qos=standard`, as in kasbex-essd and
      wescon-tools). Scripts are in `scripts/omnium_reference/`:
      - `omnium_shard_run.py` uses omnium's own RunControl/run_task but runs one shard of an
        expt's tasks (index % nshards).
      - `submit_all.sh` chains five stages with `afterany`:
        - s1: raw-input analyses, 22×10 shards
        - s2: cloud-dependent analyses + `mse_combined`, 22×10
        - s3: `*_combined`, 22 jobs
        - s4: suite plots per group (`suite_cmds.txt`), mirroring the archived `om_v0.11.2.0`
          suite dirs
        - s5: `energy_balance` + `summary_correlations` (`final_cmds.txt`)
      Logs are in `<workspace>/slurm_output/`, and omnium's per-shard logs are in
      `<workspace>/.omnium/log/`. Next: sanity diff vs the archived outputs.
      - s1 result: everything ran except `restart_dump_summary_info`, which failed for all 22
        expts. iris 2.4 returns masked arrays (nothing masked), and `masked * units(...)`
        drops the pint unit, so metpy raised. Fixed in `46ca3fe` (`np.ma.getdata`) and rerun
        (job 61562209); s5 waits for it. Compared with the archived 2019 HDFs (14 expts):
        `dynamic_summary` is bit-identical; in `thermodynamic_summary`, Tsurf is identical and
        the metpy-derived LCL/LFC/LNB/CAPE/CIN differ by ≤3e-4 relative (LNB ≤3.7 Pa, CAPE
        ≤0.008 J/kg). That points to a different metpy version in 2019; that env can't be
        recovered, since the archived `conda_list.txt` was overwritten in 2024. Treat as a
        known, explained env difference.
      - Old pickles (2019) reference `cloud_tracking.cloud_tracking`, renamed to
        `cloud_tracking.tracking` in cloud_tracking `8b240fa` (Mar 2020). Loading them needs
        `sys.modules['cloud_tracking.cloud_tracking'] = cloud_tracking.tracking`.
- [ ] Step (a) refactor, bit-identical in omnium env.
- [ ] `remake` env + remakefile + `output_remake`; `compare_outputs.py` report.

How to run omnium in the workspace (pixi tasks run from the project root, so call python
directly):

```sh
cd /gws/ssde/j25b/afesp/users/mmuetz/scaffold_omnium_remake/u-be530
~/.pixi/bin/pixi run --manifest-path ~/projects/local/scaffold_analysis/pixi.toml -e omnium \
    python ~/projects/local/omnium/bin/omnium run -s production -t expt --all S0W0Forced
```

## 7. Port design as built (branch `remake_port`, worktree `../scaffold_analysis_port`)

**Change from §2:** the analysers are not rewritten as plain functions. Instead,
`scaffold/analysis_base.py` is a ~200-line omnium-free stand-in for the omnium `Analyser`
framework, with the same interface the analysers use (`self.settings`, `self.task`,
`self.suite`, `ExptList`, `load_cubes`, `save_results_cubes`, `file_path`, `save_text`). The
caller passes explicit inputs and outputs to `run_analyser`. The ~6.6k lines of science code
are unchanged apart from imports, so results stay as close as possible to the reference.

The goals of §2 are kept:
- **No settings hashing:** `scaffold/config.py` holds a plain `SETTINGS` dict, and each rule
  declares the settings keys it uses in `uses=` (with the analyser class).
- **Explicit task lists:** inputs are static lists from `config` (runid ranges per expt), with
  the subsets each analyser actually reads (`240 <= r < 480` for `*_combined`,
  `240 < r <= 480` for `restart_dump_summary_info`, `r >= 480` for `cloud_track_analysis`, ...).
- **Real dependencies:** the hidden ones are declared, including vertlevs, the UM conf, the csv
  side outputs read by `energy_balance` and the `org_data.csv` read by `summary_correlations`.
- **Output tree:** omnium's layout minus the version-hash level:
  `output_remake/<expt>/...` and `output_remake/suite_<expts>/...`.

Removed: FF->nc converters, omnium settings, and the omnium-driven tests. `dump_extract_state`
now writes `local_conf/` next to its output, not into the suite dir.

Remakefile: `remakefiles/u_be530.py` (13,113 tasks). Run with:

```sh
cd ~/projects/local/scaffold_analysis_port
~/.pixi/bin/pixi run --manifest-path ../scaffold_analysis/pixi.toml -e remake \
    env PYTHONPATH=$PWD:$PWD/../cloud_tracking MPLBACKEND=Agg \
    remake run remakefiles/u_be530.py -E slurm
```

Checks so far: `cloud_analysis` (S0W0Forced atmos.240) through `run_analyser` is identical
to the omnium reference in the omnium env. In the remake env it is close: the only difference
is the time `calendar` attribute (`standard` vs `gregorian`, iris 3).
