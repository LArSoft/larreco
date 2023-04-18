cet_enable_asserts()

cet_make_library(LIBRARY_NAME CandidateHitFinderTool INTERFACE
  SOURCE ICandidateHitFinder.h
  LIBRARIES INTERFACE
  lardataobj::RecoBase
)

cet_write_plugin_builder(lar::CandidateHitFinderTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larreco::CandidateHitFinderTool
)

cet_make_library(LIBRARY_NAME PeakFitterTool INTERFACE
  SOURCE IPeakFitter.h
  LIBRARIES INTERFACE
  larreco::CandidateHitFinderTool
)

cet_write_plugin_builder(lar::PeakFitterTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larreco::PeakFitterTool)

cet_make_library(LIBRARY_NAME WaveformTool INTERFACE
  SOURCE IWaveformTool.h
)

cet_write_plugin_builder(lar::WaveformTool art::tool Modules
  INSTALL_BUILDER
)

cet_collect_plugin_builders(Modules lar::HitFinderToolBuilders
  LIST
  lar::CandidateHitFinderTool
  lar::PeakFitterTool
  lar::WaveformTool
)

include(lar::HitFinderToolBuilders)

cet_build_plugin(CandHitDerivative lar::CandidateHitFinderTool
  LIBRARIES PRIVATE
  larreco::WaveformTool
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art_plugin_support::toolMaker
  cetlib_except::cetlib_except
  ROOT::Hist
)

cet_build_plugin(CandHitMorphological lar::CandidateHitFinderTool
  LIBRARIES PRIVATE
  larreco::WaveformTool
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Utilities
  art_plugin_support::toolMaker
  cetlib_except::cetlib_except
  ROOT::Hist
)

cet_build_plugin(CandHitStandard lar::CandidateHitFinderTool
  LIBRARIES PRIVATE
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  art_plugin_support::toolMaker
)

cet_build_plugin(PeakFitterGaussElimination lar::PeakFitterTool
  LIBRARIES PRIVATE
  larreco::HitFinder
)

cet_build_plugin(PeakFitterGaussian lar::PeakFitterTool
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Hist
)

cet_build_plugin(PeakFitterMrqdt lar::PeakFitterTool
  LIBRARIES PRIVATE
  larreco::CandidateHitFinderTool  
  larreco::RecoAlg
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  larvecutils::MarqFitAlg
  cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger
)

cet_build_plugin(WaveformTools lar::WaveformTool
  LIBRARIES PRIVATE
  ROOT::Core
)

install_headers()
install_source()
install_fhicl()
