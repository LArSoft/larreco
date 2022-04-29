cet_make_library(SOURCE
  GaussianEliminationAlg.cxx
  HitAnaAlg.cxx
  HitFilterAlg.cxx
  RFFHitFinderAlg.cxx
  RFFHitFitter.cxx
  RegionAboveThresholdFinder.cxx
  LIBRARIES
  PUBLIC
  lardataobj::MCBase
  lardataobj::RecoBase 
  larcoreobj::SimpleTypesAndConstants
  PRIVATE
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Tree
)

add_subdirectory(HitFinderTools)

cet_build_plugin(APAHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RawData
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
)

cet_build_plugin(DPRawHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Services_System_TriggerNamesService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(DisambigCheater art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(FFTHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::MathCore
  ROOT::Matrix
)

cet_build_plugin(GausHitFinderAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(GausHitFinder art::SharedProducer
  LIBRARIES PRIVATE
  larreco::HitFinder
  larreco::CandidateHitFinderTool
  larreco::PeakFitterTool
  lardata::ArtDataHelper
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Utilities
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::MathCore
  TBB::tbb
)

cet_build_plugin(HitAnaModule art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::HitFinder
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::MCBase
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(HitCheater art::EDProducer
  LIBRARIES PRIVATE
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataalg::UtilitiesHeaders
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(HitFinderAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(HitFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::ArtDataHelper
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MCHitAnaExample art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Core
  ROOT::Hist
)

cet_build_plugin(MCHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::Simulation
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(RFFHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::HitFinder
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(RawHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  lardataobj::RawData
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(TTHitFinder art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  lardataobj::RawData
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
