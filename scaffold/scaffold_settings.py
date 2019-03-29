from omnium import AnalysisSettings


production_settings = AnalysisSettings(dict(
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
    start_runid=240,
    end_runid=480,
    # TODO: See if I can get this working in settings.
    # cloud_track_touching=True,
    # cloud_track_touching_diagonal=True,
))

test_settings = AnalysisSettings(dict(
    version = 1,
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
    start_runid=0,
    end_runid=int(1e6),
))
