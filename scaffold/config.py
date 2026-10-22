"""u-be530 experiment layout and analysis settings (replaces omnium AnalysisSettings).

Settings are a plain dict: no JSON, no hashing, no version_dir. Rules declare the
settings they read with remake's `uses=`.
"""
FORCED = ['S0W0Forced', 'S4W0Forced', 'S0W5Forced', 'S4W5Forced']
RWPS = ['RWP_C{}'.format(i) for i in range(1, 11)]
W0 = ['S0W0rel2_S0W0Forced', 'S0W0rel2_S4W0Forced', 'S4W0rel2_S0W0Forced', 'S4W0rel2_S4W0Forced']
W5 = ['S0W5rel2_S0W5Forced', 'S0W5rel2_S4W5Forced', 'S4W5rel2_S0W5Forced', 'S4W5rel2_S4W5Forced']
FORCED_PLUS_RWPS = FORCED + RWPS
ALL_EXPTS = FORCED + RWPS + W0 + W5

GROUPS = {
    'FORCED': FORCED,
    'RWPs': RWPS,
    'FORCED_PLUS_RWPS': FORCED_PLUS_RWPS,
    'W0': W0,
    'W5': W5,
}

# FORCED and RWP runs are 25 days long; the rel2 runs (W0, W5) are 20 days.
NDAYS = {expt: 25 for expt in FORCED + RWPS}
NDAYS.update({expt: 20 for expt in W0 + W5})


def pp1_runids(expt):
    """pp1 output files: one every 4 h (30-min output to day 20, 5-min after)."""
    return list(range(0, NDAYS[expt] * 24, 4))


def pp2_runids(expt):
    """pp2 output files: one per day."""
    return list(range(0, NDAYS[expt] * 24, 24))


def dump_runids(expt):
    """Restart dumps atmosa_daNNN: one per day, from day 1."""
    return list(range(24, NDAYS[expt] * 24 + 1, 24))


SETTINGS = dict(
    # cloud_analysis
    height_levels=[15, 17, 19],
    qcl_threshs=[4.5e-6, 5e-6, 5.5e-6],
    w_threshs=[0.9, 1., 1.1],
    # mass_flux_spatial_scales_analysis
    npow=4,
    # profile_analysis
    mf_profile_xlim=(0, 10000),
    qcl_thresh=5e-6,
    w_thresh=1.,
    # multiple.
    start_runid=240,
    end_runid=480,
)
