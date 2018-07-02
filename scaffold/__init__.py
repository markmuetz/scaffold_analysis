from scaffold.version import VERSION
from scaffold.cycle.cloud_analysis import CloudAnalyser
from scaffold.cycle.mass_flux_analysis import MassFluxAnalyser
from scaffold.cycle.mass_flux_spatial_scales_analysis import MassFluxSpatialScalesAnalyser
from scaffold.cycle.org_analysis import OrgAnalyser
from scaffold.cycle.profile_analysis import ProfileAnalyser
from scaffold.cycle.restart_dump_analysis import RestartDumpAnalyser
from scaffold.cloud_track_analysis import CloudTrackAnalyser
from scaffold.dummy_analysis import DummyAnalyser
from scaffold.expt.mass_flux_combined import MassFluxCombinedAnalysis
from scaffold.mass_flux_plot import MassFluxPlotter
from scaffold.mass_flux_spatial_scales_combined import MassFluxSpatialScalesCombined
from scaffold.mass_flux_spatial_scales_plot import MassFluxSpatialScalesPlotter
from scaffold.expt.mse_combined import MseCombinedAnalysis
from scaffold.mse_plot import MsePlotter
from scaffold.expt.org_combined import OrgCombined
from scaffold.org_plot import OrgPlotter
from scaffold.precip_plot import PrecipPlot
from scaffold.precip_hovmoller import PrecipHovmollerAnalyser
from scaffold.profile_plot import ProfilePlotter
from scaffold.surf_flux_analysis import SurfFluxAnalyser
from scaffold.surf_flux_plot import SurfFluxPlot


__version__ = VERSION

analysis_classes = [
    CloudAnalyser,
    CloudTrackAnalyser,
    DummyAnalyser,
    MassFluxAnalyser,
    MassFluxCombinedAnalysis,
    MassFluxPlotter,
    MassFluxSpatialScalesAnalyser,
    MassFluxSpatialScalesCombined,
    MassFluxSpatialScalesPlotter,
    MseCombinedAnalysis,
    MsePlotter,
    OrgAnalyser,
    OrgCombined,
    OrgPlotter,
    PrecipPlot,
    PrecipHovmollerAnalyser,
    ProfileAnalyser,
    ProfilePlotter,
    RestartDumpAnalyser,
    SurfFluxAnalyser,
    SurfFluxPlot,
]
