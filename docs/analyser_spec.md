# Analyser spec: what each scaffold analyser really reads and writes

The rule-by-rule spec for the remake port. It is derived from the code at `7dd82b2`/`46ca3fe`
and checked against the omnium task list for u-be530 (640 expt tasks for a 25-day expt).

Conventions:
- `H` = `share/data/history/<expt>`, `W` = `work/19700101T0000Z/<expt>_atmos`,
  `O` = output dir for `<expt>`, `SO` = suite output dir for `<group>`.
- Settings are the production values in `scaffold/scaffold_settings.py`:
  height_levels=[15,17,19], w_threshs=[0.9,1,1.1], qcl_threshs=[4.5e-6,5e-6,5.5e-6], npow=4,
  qcl_thresh=5e-6, w_thresh=1, start_runid=240, end_runid=480, mf_profile_xlim=(0,10000).
- "UM conf" = `W/rose-app-run.conf` (dx, lx via omnium `Expt`).
  "vertlevs" = `app/um/file/rce_vertlevs.nml` (via `VertLev`).

## Runids

| stream | 25-day expts (FORCED, RWPs) | 20-day expts (rel2: W0, W5) |
|---|---|---|
| pp1 (every 4 h, 30-min steps to day 20, 5-min after) | 000–596 step 4 (150) | 000–476 step 4 (120) |
| pp2 (daily) | 000–576 step 24 (25) | 000–456 step 24 (20) |
| dumps `atmosa_daNNN` | 024–600 step 24 (25) | 024–480 step 24 (20) |

## Expt analysers

| analyser | matrix | real inputs | settings | declared outputs | side outputs / notes |
|---|---|---|---|---|---|
| cloud_analysis | expt × pp1 runid | `H/atmos.R.pp1.nc` | height_levels, w_threshs, qcl_threshs | `O/atmos.R.cloud_analysis.nc` | — |
| mass_flux_analysis | expt × pp1 runid | `O/atmos.R.cloud_analysis.nc`, vertlevs | — | `O/atmos.R.mass_flux_analysis.nc` | vertlevs not declared |
| mass_flux_spatial_scales_analysis | expt × pp1 runid | cloud_analysis R, vertlevs | npow | `O/atmos.R.mass_flux_spatial_scales_analysis.nc` | |
| org_analysis | expt × pp1 runid | cloud_analysis R, UM conf (lx) | — | `O/atmos.R.org_analysis.nc` | `org_analysis_dists_hist.pdf` |
| profile_analysis | expt × pp2 runid | `H/atmos.R.pp2.nc`, UM conf | w_thresh, qcl_thresh | `O/atmos.R.profile_analysis.nc` | 9 `profile_analysis_*.png`, `profile_analysis_energy_loss.txt` (overwritten by every task: last runid wins) |
| restart_dump_analysis | expt × dumps 408,432,456,480 (glob `da4??`) | `H/atmosa_daR.nc` | — | `O/atmos.R.restart_dump_analysis.nc` | |
| dump_slice_analysis | expt × dumps 456, 480 (glob `da4[58][60]`) | `H/atmosa_daR.nc`, vertlevs | — | `O/atmos.R.dump_slice_analysis.dummy` | many PDFs under `dump_slice_analysis_/R/{figs,xy,...}/`; hard-coded slices per expt |
| surf_flux_analysis | expt | `W/atmos.pp3.nc` | — | `O/atmos.surf_flux_analysis.nc` | 2 png |
| precip_hovmoller_analysis | expt | `W/atmos.pp3.nc` | — | `O/atmos.precip_hovmoller_analysis.dummy` | 2 png |
| dump_extract_state | expt | dumps with runid ≥ 264 | — (class `min_runid=264`) | `O/atmos.dump_extract_state.dummy` | **writes `<suite_dir>/local_conf/rose-app-relax_to_<expt>.conf`, outside the output tree** |
| restart_dump_summary_info | expt | dumps with `240 < R <= 480` (note `<`/`<=`) | start/end_runid | `O/atmos.restart_dump_summary_info.hdf` (2 tables) | metpy; fixed in `46ca3fe` |
| mass_flux_combined | expt | mass_flux_analysis with `240 <= R < 480` | start/end_runid | `O/atmos.mass_flux_combined.nc` | all 150 declared as inputs, 60 used |
| mass_flux_spatial_scales_combined | expt | mfss_analysis with `240 <= R < 480` | start/end_runid | `O/atmos.mass_flux_spatial_scales_combined.nc` | as above |
| org_combined | expt | org_analysis with `240 <= R < 480` | start/end_runid | `O/atmos.org_combined.nc` | as above |
| mse_combined | expt | all restart_dump_analysis (4) | — | `O/atmos.mse_combined.nc` | |
| cloud_track_analysis | expt (25-day only) | cloud_analysis with R ≥ 480 (class `min_runid`), UM conf (dx) | — | `O/atmos.cloud_track_analysis.{all_stats,trackers}.pkl` | `O/atmos.cloud_tracking_z1_t1.{txt,...}` via `file_path('cloud_tracking')`; only key (1,1) |

## Suite analysers

Groups: FORCED (4), RWPs (10), FORCED_PLUS_RWPS (14), W0 (4 rel2), W5 (4 rel2).
The output dir is `suite_<expt1>_<expt2>...` (in remake: `suite_<group>`).

| analyser | groups run | real inputs per expt | declared output | side outputs / notes |
|---|---|---|---|---|
| mass_flux_plot | FORCED, RWPs, W0, W5 | `O/atmos.mass_flux_combined.nc`, UM conf | `.dummy` | many pdf/png, `mass_flux_plot_mean_num_clds.hdf`, `mass_flux_plot_mf_linregress.csv` |
| mass_flux_spatial_scales_plot | FORCED, RWPs, W0, W5 | `O/atmos.mass_flux_spatial_scales_combined.nc` | `.dummy` | png |
| profile_plot | FORCED, RWPs, W0, W5 | `O/atmos.456.profile_analysis.nc`, vertlevs | `.dummy` | 9+ png incl. `UCP_shear_profiles.png` |
| precip_plot | FORCED, RWPs, W0, W5 | `W/atmos.pp3.nc` | `.dummy` | `time_index*.pdf`, `hovmollers.pdf`, `max_precip.csv`; **writes `<expt>_cached_hov_{x,y}.nc` (numpy pickles despite `.nc`)**; the cache read-back is commented out, so it is write-only |
| org_plot | FORCED, RWPs, W0, W5, FORCED_PLUS_RWPS | `O/atmos.org_combined.nc`, UM conf (lx) | `.dummy` | `org_data.csv`, `org_data[_xpoints]_<expt>_z<g>.np` (pickles), pdfs. Output changed after 2019 (`c02d325`, `7dd82b2`): matches archive `v0.11.3.0`, not `v0.11.2.0` |
| surf_flux_plot | FORCED, RWPs, W0, W5 | `W/atmos.pp3.nc` | `.dummy` | `surf_flux_plot_{final_day,half_sim}_energy_flux.csv`, png |
| mse_plot | FORCED, RWPs, W0, W5 | `O/atmos.mse_combined.nc` | `.dummy` | png |
| dump_profile_plot | FORCED, W0, W5 | `H/atmosa_da480.nc` | `.dummy` | `melting_level.csv`, `density_at_12p5km.csv`, png |
| relaxation_plot | FORCED, W0, W5 | `H/atmos.456.pp2.nc` | `.dummy` | `relaxation_plot_final_day_relaxation_energy_flux.csv`, png |
| cloud_track_plot | FORCED | `O/atmos.cloud_track_analysis.all_stats.pkl` | `.dummy` | png |
| cloud_track_lifecycle_plot | FORCED | `O/atmos.cloud_track_analysis.trackers.pkl` | `.dummy` | png |
| energy_balance | FORCED, W0, W5 | **`SO/surf_flux_plot_final_day_energy_flux.csv`, `SO/relaxation_plot_final_day_relaxation_energy_flux.csv`** (side outputs of two other plots) | `SO/energy_balance.hdf` | `energy_balance.tex` |
| summary_correlations | — (fixed list of 14 expts) | **`suite_<FORCED_PLUS_RWPS>/org_data.csv`**; per expt `atmos.restart_dump_summary_info.hdf`, `atmos.mass_flux_combined.nc`, `atmos.cloud_tracking_z1_t1.txt` | `suite_summary_correlations/{atmos.summary_correlations.dummy, all_col_values.hdf}` | `summary_info_*.pdf/png` (paper Figs 7, 9), `ordered_pvals.csv`. Its declared input is a hand-made `summary_correlations.dummy`; the real inputs are found via `cd` |

Not in the u-be530 workflow: `dump_entr` (needs da480 for a group, no runcontrol entry),
`mass_flux_spatial_scales_point`, `CycleConverter`/`CycleDumpConverter`/`ExptConverter`
(FF → nc, already done).

## Hidden dependencies to make explicit in remake

1. `energy_balance` ← csv side outputs of `surf_flux_plot` and `relaxation_plot`.
2. `summary_correlations` ← `org_plot[FORCED_PLUS_RWPS]/org_data.csv`,
   `restart_dump_summary_info`, `mass_flux_combined`, and the `cloud_track_analysis` txt for
   all 14 expts.
3. Undeclared reads of vertlevs and the UM conf (mass_flux, mfss, org, profile, dump_slice,
   cloud_track, several plots).
4. `*_combined` and `restart_dump_summary_info` declare all runids but use a subset. Declare
   only the subset.
5. `dump_extract_state` writes into the suite dir. In remake it writes into the output tree.
6. `profile_analysis_energy_loss.txt` and its pngs are rewritten by every runid task, so the
   result depends on the order tasks run in. In remake, write per-runid files, or only for a
   chosen runid.
