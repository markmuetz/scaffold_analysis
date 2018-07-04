from omnium import AnalyserSetting

import scaffold


production_settings = AnalyserSetting(scaffold, dict(
    # cloud_analysis
    height_levels = [15, 17, 19],
    qcl_threshs=[4.5e-6,5e-6,5.5e-6],
    w_threshs=[0.9,1.,1.1],
    # mass_flux_spatial_scales_analysis
    npow=4,
    # profile_analysis
    mf_profile_xlim=(0,10000),
    qcl_thresh=5e-6,
    w_thresh=1.,
    # multiple.
    start_runid=120,
))
