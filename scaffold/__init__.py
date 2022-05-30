import iris

from scaffold.cycle.cycle_converter import CycleConverter, CycleDumpConverter
from scaffold.expt.cloud_analysis import CloudAnalyser
from scaffold.expt.cloud_track_analysis import CloudTrackAnalyser
from scaffold.expt.dump_extract_state import DumpExtractState
from scaffold.expt.dump_slice_analysis import DumpSliceAnalyser
from scaffold.expt.expt_converter import ExptConverter
from scaffold.expt.mass_flux_analysis import MassFluxAnalyser
from scaffold.expt.mass_flux_combined import MassFluxCombinedAnalysis
from scaffold.expt.mass_flux_spatial_scales_analysis import MassFluxSpatialScalesAnalyser
from scaffold.expt.mass_flux_spatial_scales_combined import MassFluxSpatialScalesCombined
from scaffold.expt.mse_combined import MseCombinedAnalysis
from scaffold.expt.org_analysis import OrgAnalyser
from scaffold.expt.org_combined import OrgCombined
from scaffold.expt.precip_hovmoller import PrecipHovmollerAnalyser
from scaffold.expt.profile_analysis import ProfileAnalyser
from scaffold.expt.restart_dump_analysis import RestartDumpAnalyser
from scaffold.expt.restart_dump_summary_info import RestartDumpSummaryInfo
from scaffold.expt.surf_flux_analysis import SurfFluxAnalyser
from scaffold.scaffold_settings import production_settings, test_settings
from scaffold.suite.cloud_track_plot import CloudTrackPlot
from scaffold.suite.cloud_track_lifecycles_plot import CloudTrackLifecyclesPlot
from scaffold.suite.dump_entr import DumpEntr
from scaffold.suite.dump_profile_plot import DumpProfilePlotter
from scaffold.suite.energy_balance import EnergyBalance
from scaffold.suite.mass_flux_plot import MassFluxPlotter
from scaffold.suite.mass_flux_spatial_scales_plot import MassFluxSpatialScalesPlotter
from scaffold.suite.mse_plot import MsePlotter
from scaffold.suite.org_plot import OrgPlotter
from scaffold.suite.precip_plot import PrecipPlot
from scaffold.suite.profile_plot import ProfilePlotter
from scaffold.suite.relaxation_plot import RelaxationPlot
from scaffold.suite.summary_correlations import SummaryCorrelations
from scaffold.suite.surf_flux_plot import SurfFluxPlot
from scaffold.version import VERSION

# Check and apply LBHEM fix.
# DISABLED due to not being able to pip install iris.
# assert iris.HAS_LBHEM_FIX
# iris.site_configuration['fix_lbhem_eq_4'] = 3

__version__ = VERSION

analysis_settings = {
    'default': production_settings,
    'production': production_settings,
    'test': test_settings,
}

analysis_settings_filename = 'omnium_output/{version_dir}/settings.json'

analyser_classes = [
    CycleConverter,
    CycleDumpConverter,
    CloudAnalyser,
    CloudTrackAnalyser,
    CloudTrackLifecyclesPlot,
    CloudTrackPlot,
    DumpEntr,
    DumpExtractState,
    DumpProfilePlotter,
    DumpSliceAnalyser,
    EnergyBalance,
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
    RelaxationPlot,
    RestartDumpAnalyser,
    RestartDumpSummaryInfo,
    SummaryCorrelations,
    SurfFluxAnalyser,
    SurfFluxPlot,
]
