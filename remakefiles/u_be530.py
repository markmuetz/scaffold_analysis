"""remake pipeline for the u-be530 RCE analysis (replaces the omnium workflow).

Rules follow docs/analyser_spec.md: each task gets exactly the inputs its analyser reads,
including the dependencies omnium left implicit (vertlevs, UM conf, csv side outputs, the
org_data.csv read by summary_correlations).

Output layout matches omnium's, minus the version_dir level:
    OUT/<expt>/atmos.<runid>.<analysis>.nc
    OUT/suite_<expt1>_<expt2>.../...       (same suite dir names as omnium)
    OUT/suite_summary_correlations/...

Run (from the scaffold_analysis_port worktree):
    PYTHONPATH=$PWD:$PWD/../cloud_tracking pixi run --manifest-path ../scaffold_analysis/pixi.toml \
        -e remake remake run remakefiles/u_be530.py -E slurm
"""
from remake import Remake, rule

from scaffold import ANALYSERS
from scaffold.config import (ALL_EXPTS, FORCED, RWPS, GROUPS, SETTINGS, pp1_runids, pp2_runids,
                             dump_runids)

SUITE_DIR = '/gws/ssde/j25b/afesp/users/mmuetz/scaffold_omnium_remake/u-be530'
OUT = '/gws/ssde/j25b/afesp/users/mmuetz/scaffold_omnium_remake/output_remake'

HIST = SUITE_DIR + '/share/data/history'
WORK = SUITE_DIR + '/work/19700101T0000Z'
VERTLEVS = SUITE_DIR + '/app/um/file/rce_vertlevs.nml'

rmk = Remake(config={
    'slurm': {
        'account': 'afesp',
        'partition': 'standard',
        'qos': 'standard',
        'time': '12:00:00',
        'mem': '32G',
    },
})

# Expts with 25-day runs (the only ones with 5-min pp1 output after day 20 -> cloud tracking).
LONG_EXPTS = FORCED + RWPS


def rid(runid):
    return '{:03d}'.format(runid)


def suite_dir(group):
    return '{}/suite_{}'.format(OUT, '_'.join(GROUPS[group]))


def um_conf(expt):
    return '{}/{}_atmos/rose-app-run.conf'.format(WORK, expt)


def numbered(paths):
    """{'in000': path0, 'in001': path1, ...}: analyser inputs, in order."""
    return {'in{:03d}'.format(i): str(p) for i, p in enumerate(paths)}


def analyser_files(inputs):
    return [inputs[k] for k in sorted(inputs) if k.startswith('in')]


def declared_outputs(outputs):
    return [outputs[k] for k in sorted(outputs) if k.startswith('out')]


def settings_subset(keys):
    return {k: SETTINGS[k] for k in keys}


# ---------------------------------------------------------------------------
# Expt rules: matrix (expt, runid)
# ---------------------------------------------------------------------------

def _per_runid_rule(name, matrix, inputs, outputs, settings_keys=(), depends_on=None, config=None):
    @rule(inputs=inputs, outputs=outputs, matrix=matrix, depends_on=depends_on,
          uses={'analyser': ANALYSERS[name], 'settings': settings_subset(settings_keys),
                'suite_root': SUITE_DIR, 'analyser_files': analyser_files,
                'declared_outputs': declared_outputs},
          config=config, name=name)
    def run(inputs, outputs, expt, runid):
        from scaffold.analysis_base import run_analyser
        run_analyser(analyser, suite_root, analyser_files(inputs), declared_outputs(outputs),
                     settings, expt=expt, runid=int(runid))
    return run


PP1_MATRIX = [{'expt': e, 'runid': rid(r)} for e in ALL_EXPTS for r in pp1_runids(e)]
PP2_MATRIX = [{'expt': e, 'runid': rid(r)} for e in ALL_EXPTS for r in pp2_runids(e)]
RESTART_DUMP_MATRIX = [{'expt': e, 'runid': rid(r)} for e in ALL_EXPTS
                       for r in dump_runids(e) if 400 <= r < 500]   # omnium glob atmosa_da4??
DUMP_SLICE_MATRIX = [{'expt': e, 'runid': rid(r)} for e in ALL_EXPTS
                     for r in dump_runids(e) if r in (456, 480)]  # glob atmosa_da4[58][60]


def _out(analysis, ext='nc'):
    return lambda expt, runid: {'out0': '{}/{}/atmos.{}.{}.{}'.format(OUT, expt, runid,
                                                                      analysis, ext)}


cloud_analysis = _per_runid_rule(
    'cloud_analysis', PP1_MATRIX,
    inputs=lambda expt, runid: numbered(['{}/{}/atmos.{}.pp1.nc'.format(HIST, expt, runid)]),
    outputs=_out('cloud_analysis'),
    settings_keys=('height_levels', 'w_threshs', 'qcl_threshs'))


def _cloud_analysis_input(expt, runid, extra=None):
    inputs = numbered(['{}/{}/atmos.{}.cloud_analysis.nc'.format(OUT, expt, runid)])
    inputs.update(extra or {})
    return inputs


mass_flux_analysis = _per_runid_rule(
    'mass_flux_analysis', PP1_MATRIX,
    inputs=lambda expt, runid: _cloud_analysis_input(expt, runid, {'vertlevs': VERTLEVS}),
    outputs=_out('mass_flux_analysis'), depends_on=[cloud_analysis])

mass_flux_spatial_scales_analysis = _per_runid_rule(
    'mass_flux_spatial_scales_analysis', PP1_MATRIX,
    inputs=lambda expt, runid: _cloud_analysis_input(expt, runid, {'vertlevs': VERTLEVS}),
    outputs=_out('mass_flux_spatial_scales_analysis'), settings_keys=('npow',),
    depends_on=[cloud_analysis])

org_analysis = _per_runid_rule(
    'org_analysis', PP1_MATRIX,
    inputs=lambda expt, runid: _cloud_analysis_input(expt, runid, {'um_conf': um_conf(expt)}),
    outputs=_out('org_analysis'), depends_on=[cloud_analysis])

profile_analysis = _per_runid_rule(
    'profile_analysis', PP2_MATRIX,
    inputs=lambda expt, runid: dict(numbered(['{}/{}/atmos.{}.pp2.nc'.format(HIST, expt, runid)]),
                                    um_conf=um_conf(expt)),
    outputs=_out('profile_analysis'), settings_keys=('w_thresh', 'qcl_thresh'))

restart_dump_analysis = _per_runid_rule(
    'restart_dump_analysis', RESTART_DUMP_MATRIX,
    inputs=lambda expt, runid: numbered(['{}/{}/atmosa_da{}.nc'.format(HIST, expt, runid)]),
    outputs=_out('restart_dump_analysis'))

dump_slice_analysis = _per_runid_rule(
    'dump_slice_analysis', DUMP_SLICE_MATRIX,
    inputs=lambda expt, runid: dict(numbered(['{}/{}/atmosa_da{}.nc'.format(HIST, expt, runid)]),
                                    vertlevs=VERTLEVS),
    outputs=_out('dump_slice_analysis', 'dummy'))


# ---------------------------------------------------------------------------
# Expt rules: matrix (expt), usually fan-in over runids
# ---------------------------------------------------------------------------

def _per_expt_rule(name, expts, inputs, outputs, settings_keys=(), depends_on=None, config=None):
    @rule(inputs=inputs, outputs=outputs, matrix={'expt': list(expts)}, depends_on=depends_on,
          uses={'analyser': ANALYSERS[name], 'settings': settings_subset(settings_keys),
                'suite_root': SUITE_DIR, 'analyser_files': analyser_files,
                'declared_outputs': declared_outputs},
          config=config, name=name)
    def run(inputs, outputs, expt):
        from scaffold.analysis_base import run_analyser
        run_analyser(analyser, suite_root, analyser_files(inputs), declared_outputs(outputs),
                     settings, expt=expt)
    return run


def _expt_out(*names):
    return lambda expt: {'out{}'.format(i): '{}/{}/{}'.format(OUT, expt, n)
                         for i, n in enumerate(names)}


def _pp3(expt):
    return numbered(['{}/{}_atmos/atmos.pp3.nc'.format(WORK, expt)])


def _runid_window(runids, lo, hi, lo_incl=True, hi_incl=False):
    return [r for r in runids
            if (lo <= r if lo_incl else lo < r) and (r <= hi if hi_incl else r < hi)]


START, END = SETTINGS['start_runid'], SETTINGS['end_runid']

surf_flux_analysis = _per_expt_rule(
    'surf_flux_analysis', ALL_EXPTS, inputs=_pp3, outputs=_expt_out('atmos.surf_flux_analysis.nc'))

precip_hovmoller_analysis = _per_expt_rule(
    'precip_hovmoller_analysis', ALL_EXPTS, inputs=_pp3,
    outputs=_expt_out('atmos.precip_hovmoller_analysis.dummy'))

dump_extract_state = _per_expt_rule(
    'dump_extract_state', ALL_EXPTS,
    inputs=lambda expt: numbered(['{}/{}/atmosa_da{}.nc'.format(HIST, expt, rid(r))
                                  for r in dump_runids(expt) if r >= 264]),
    outputs=lambda expt: {
        'out0': '{}/{}/atmos.dump_extract_state.dummy'.format(OUT, expt),
        'side_conf': '{}/{}/local_conf/rose-app-relax_to_{}.conf'.format(OUT, expt, expt)})

restart_dump_summary_info = _per_expt_rule(
    'restart_dump_summary_info', ALL_EXPTS,
    inputs=lambda expt: numbered(['{}/{}/atmosa_da{}.nc'.format(HIST, expt, rid(r))
                                  for r in _runid_window(dump_runids(expt), START, END,
                                                         lo_incl=False, hi_incl=True)]),
    outputs=_expt_out('atmos.restart_dump_summary_info.hdf'),
    settings_keys=('start_runid', 'end_runid'))


def _combined_inputs(analysis):
    def inputs(expt):
        return numbered(['{}/{}/atmos.{}.{}.nc'.format(OUT, expt, rid(r), analysis)
                         for r in _runid_window(pp1_runids(expt), START, END)])
    return inputs


mass_flux_combined = _per_expt_rule(
    'mass_flux_combined', ALL_EXPTS, inputs=_combined_inputs('mass_flux_analysis'),
    outputs=_expt_out('atmos.mass_flux_combined.nc'), settings_keys=('start_runid', 'end_runid'),
    depends_on=[mass_flux_analysis])

mass_flux_spatial_scales_combined = _per_expt_rule(
    'mass_flux_spatial_scales_combined', ALL_EXPTS,
    inputs=_combined_inputs('mass_flux_spatial_scales_analysis'),
    outputs=_expt_out('atmos.mass_flux_spatial_scales_combined.nc'),
    settings_keys=('start_runid', 'end_runid'), depends_on=[mass_flux_spatial_scales_analysis])

org_combined = _per_expt_rule(
    'org_combined', ALL_EXPTS, inputs=_combined_inputs('org_analysis'),
    outputs=_expt_out('atmos.org_combined.nc'), settings_keys=('start_runid', 'end_runid'),
    depends_on=[org_analysis])

mse_combined = _per_expt_rule(
    'mse_combined', ALL_EXPTS,
    inputs=lambda expt: numbered(['{}/{}/atmos.{}.restart_dump_analysis.nc'.format(OUT, expt, rid(r))
                                  for r in dump_runids(expt) if 400 <= r < 500]),
    outputs=_expt_out('atmos.mse_combined.nc'), depends_on=[restart_dump_analysis])

cloud_track_analysis = _per_expt_rule(
    'cloud_track_analysis', LONG_EXPTS,
    inputs=lambda expt: dict(numbered(['{}/{}/atmos.{}.cloud_analysis.nc'.format(OUT, expt, rid(r))
                                       for r in pp1_runids(expt) if r >= 480]),
                             um_conf=um_conf(expt)),
    outputs=lambda expt: {
        'out0': '{}/{}/atmos.cloud_track_analysis.all_stats.pkl'.format(OUT, expt),
        'out1': '{}/{}/atmos.cloud_track_analysis.trackers.pkl'.format(OUT, expt),
        'side_txt': '{}/{}/atmos.cloud_tracking_z1_t1.txt'.format(OUT, expt)},
    depends_on=[cloud_analysis])


# ---------------------------------------------------------------------------
# Suite rules: matrix (group), fan-in over the group's expts
# ---------------------------------------------------------------------------

def _per_group_rule(name, groups, expt_input, side_outputs=(), extra_inputs=None,
                    depends_on=None, config=None):
    """One task per group; analyser input is one file per expt (expt_input(expt))."""
    def inputs(group):
        files = numbered([expt_input(expt) for expt in GROUPS[group]])
        if extra_inputs:
            files.update(extra_inputs(group))
        return files

    def outputs(group):
        outs = {'out0': '{}/atmos.{}.dummy'.format(suite_dir(group), name)}
        for i, side in enumerate(side_outputs):
            outs['side{}'.format(i)] = '{}/{}'.format(suite_dir(group), side)
        return outs

    @rule(inputs=inputs, outputs=outputs, matrix={'group': list(groups)}, depends_on=depends_on,
          uses={'analyser': ANALYSERS[name], 'suite_root': SUITE_DIR, 'analyser_files': analyser_files,
                'declared_outputs': declared_outputs, 'groups': GROUPS},
          config=config, name=name)
    def run(inputs, outputs, group):
        from scaffold.analysis_base import run_analyser
        run_analyser(analyser, suite_root, analyser_files(inputs), declared_outputs(outputs),
                     {}, expts=groups[group])
    return run


SIMPLE_GROUPS = ['FORCED', 'RWPs', 'W0', 'W5']
DUMP_GROUPS = ['FORCED', 'W0', 'W5']


def _um_confs(group):
    return {'um_conf_{}'.format(e): um_conf(e) for e in GROUPS[group]}


def _vertlevs(group):
    return {'vertlevs': VERTLEVS}


mass_flux_plot = _per_group_rule(
    'mass_flux_plot', SIMPLE_GROUPS,
    lambda expt: '{}/{}/atmos.mass_flux_combined.nc'.format(OUT, expt),
    extra_inputs=_um_confs, depends_on=[mass_flux_combined])

mass_flux_spatial_scales_plot = _per_group_rule(
    'mass_flux_spatial_scales_plot', SIMPLE_GROUPS,
    lambda expt: '{}/{}/atmos.mass_flux_spatial_scales_combined.nc'.format(OUT, expt),
    depends_on=[mass_flux_spatial_scales_combined])

profile_plot = _per_group_rule(
    'profile_plot', SIMPLE_GROUPS,
    lambda expt: '{}/{}/atmos.456.profile_analysis.nc'.format(OUT, expt),
    extra_inputs=_vertlevs, depends_on=[profile_analysis])

precip_plot = _per_group_rule(
    'precip_plot', SIMPLE_GROUPS, lambda expt: '{}/{}_atmos/atmos.pp3.nc'.format(WORK, expt),
    config={'slurm': {'mem': '64G'}})

org_plot = _per_group_rule(
    'org_plot', SIMPLE_GROUPS + ['FORCED_PLUS_RWPS'],
    lambda expt: '{}/{}/atmos.org_combined.nc'.format(OUT, expt),
    side_outputs=('org_data.csv',), extra_inputs=_um_confs, depends_on=[org_combined])

surf_flux_plot = _per_group_rule(
    'surf_flux_plot', SIMPLE_GROUPS, lambda expt: '{}/{}_atmos/atmos.pp3.nc'.format(WORK, expt),
    side_outputs=('surf_flux_plot_final_day_energy_flux.csv',))

mse_plot = _per_group_rule(
    'mse_plot', SIMPLE_GROUPS, lambda expt: '{}/{}/atmos.mse_combined.nc'.format(OUT, expt),
    depends_on=[mse_combined])

dump_profile_plot = _per_group_rule(
    'dump_profile_plot', DUMP_GROUPS, lambda expt: '{}/{}/atmosa_da480.nc'.format(HIST, expt),
    config={'slurm': {'mem': '64G'}})

relaxation_plot = _per_group_rule(
    'relaxation_plot', DUMP_GROUPS, lambda expt: '{}/{}/atmos.456.pp2.nc'.format(HIST, expt),
    side_outputs=('relaxation_plot_final_day_relaxation_energy_flux.csv',), extra_inputs=_um_confs,
    config={'slurm': {'mem': '64G'}})

cloud_track_plot = _per_group_rule(
    'cloud_track_plot', ['FORCED'],
    lambda expt: '{}/{}/atmos.cloud_track_analysis.all_stats.pkl'.format(OUT, expt),
    depends_on=[cloud_track_analysis])

cloud_track_lifecycle_plot = _per_group_rule(
    'cloud_track_lifecycle_plot', ['FORCED'],
    lambda expt: '{}/{}/atmos.cloud_track_analysis.trackers.pkl'.format(OUT, expt),
    depends_on=[cloud_track_analysis])


def _energy_balance_inputs(group):
    return numbered(['{}/surf_flux_plot_final_day_energy_flux.csv'.format(suite_dir(group)),
                     '{}/relaxation_plot_final_day_relaxation_energy_flux.csv'.format(suite_dir(group))])


@rule(inputs=_energy_balance_inputs,
      outputs=lambda group: {'out0': '{}/energy_balance.hdf'.format(suite_dir(group))},
      matrix={'group': DUMP_GROUPS}, depends_on=[surf_flux_plot, relaxation_plot],
      uses={'analyser': ANALYSERS['energy_balance'], 'suite_root': SUITE_DIR, 'analyser_files': analyser_files,
                'declared_outputs': declared_outputs, 'groups': GROUPS})
def energy_balance(inputs, outputs, group):
    from scaffold.analysis_base import run_analyser
    run_analyser(analyser, suite_root, analyser_files(inputs), declared_outputs(outputs), {},
                 expts=groups[group])


SC_EXPTS = FORCED + RWPS


def _summary_correlations_inputs():
    inputs = {'org_data': '{}/org_data.csv'.format(suite_dir('FORCED_PLUS_RWPS'))}
    for expt in SC_EXPTS:
        inputs['{}_summary_info'.format(expt)] = '{}/{}/atmos.restart_dump_summary_info.hdf'.format(OUT, expt)
        inputs['{}_mass_flux'.format(expt)] = '{}/{}/atmos.mass_flux_combined.nc'.format(OUT, expt)
        inputs['{}_cloud_tracking'.format(expt)] = '{}/{}/atmos.cloud_tracking_z1_t1.txt'.format(OUT, expt)
    return inputs


@rule(inputs=_summary_correlations_inputs(),
      outputs={'out0': OUT + '/suite_summary_correlations/atmos.summary_correlations.dummy',
               'out1': OUT + '/suite_summary_correlations/all_col_values.hdf'},
      depends_on=[org_plot, restart_dump_summary_info, mass_flux_combined, cloud_track_analysis],
      uses={'analyser': ANALYSERS['summary_correlations'], 'suite_root': SUITE_DIR, 'analyser_files': analyser_files,
                'declared_outputs': declared_outputs, 'out': OUT})
def summary_correlations(inputs, outputs):
    from scaffold.analysis_base import run_analyser
    # summary_correlations finds its real inputs relative to the dir of its first (nominal)
    # input file, so pass a path in OUT; the file itself is never opened.
    run_analyser(analyser, suite_root, [out + '/summary_correlations.dummy'],
                 declared_outputs(outputs), {})


rmk.rules_from_current_module()
