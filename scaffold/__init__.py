from scaffold.version import VERSION
from scaffold.scaffold_settings import production_settings, test_settings

from scaffold.cycle.cycle_converter import CycleConverter
from scaffold.cycle.cloud_analysis import CloudAnalyser
from scaffold.cycle.mass_flux_analysis import MassFluxAnalyser
from scaffold.cycle.mass_flux_spatial_scales_analysis import MassFluxSpatialScalesAnalyser
from scaffold.cycle.org_analysis import OrgAnalyser
from scaffold.cycle.profile_analysis import ProfileAnalyser
from scaffold.cycle.restart_dump_analysis import RestartDumpAnalyser

from scaffold.expt.expt_converter import ExptConverter
from scaffold.expt.cloud_track_analysis import CloudTrackAnalyser
from scaffold.expt.mass_flux_combined import MassFluxCombinedAnalysis
from scaffold.expt.mass_flux_spatial_scales_combined import MassFluxSpatialScalesCombined
from scaffold.expt.mse_combined import MseCombinedAnalysis
from scaffold.expt.org_combined import OrgCombined
from scaffold.expt.surf_flux_analysis import SurfFluxAnalyser
from scaffold.expt.precip_hovmoller import PrecipHovmollerAnalyser

from scaffold.suite.mass_flux_plot import MassFluxPlotter
from scaffold.suite.mass_flux_spatial_scales_plot import MassFluxSpatialScalesPlotter
from scaffold.suite.mse_plot import MsePlotter
from scaffold.suite.org_plot import OrgPlotter
from scaffold.suite.precip_plot import PrecipPlot
from scaffold.suite.profile_plot import ProfilePlotter
from scaffold.suite.surf_flux_plot import SurfFluxPlot

__version__ = VERSION

analysis_settings = {
    'default': production_settings,
    'production': production_settings,
    'test': test_settings,
}

analysis_settings_filename = 'omnium_output/{version_dir}/settings.json'

analysis_classes = [
    CycleConverter,
    CloudAnalyser,
    CloudTrackAnalyser,
    ExptConverter,
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
