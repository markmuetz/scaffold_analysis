import iris

from scaffold.expt.cloud_analysis import CloudAnalyser
from scaffold.expt.cloud_track_analysis import CloudTrackAnalyser
from scaffold.expt.dump_extract_state import DumpExtractState
from scaffold.expt.dump_slice_analysis import DumpSliceAnalyser
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

# All analysers, by analysis_name. Run one with scaffold.analysis_base.run_analyser.
ANALYSERS = {cls.analysis_name: cls for cls in [
    CloudAnalyser,
    CloudTrackAnalyser,
    CloudTrackLifecyclesPlot,
    CloudTrackPlot,
    DumpEntr,
    DumpExtractState,
    DumpProfilePlotter,
    DumpSliceAnalyser,
    EnergyBalance,
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
]}
