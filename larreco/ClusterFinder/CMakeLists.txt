add_subdirectory(TrackShowerSplitter)

cet_make_library(SOURCE ClusterCreator.cxx
  LIBRARIES
  PUBLIC
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  PRIVATE
  larreco::RecoAlg_ClusterRecoUtil
)

cet_build_plugin(BlurredClustering art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg_ClusterRecoUtil
  larreco::RecoAlg
  larreco::ClusterFinder
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardata::Utilities
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(Cluster3D art::EDProducer
  LIBRARIES PRIVATE
  larreco::ClusterFinder
  larreco::RecoAlg_ClusterRecoUtil
  larreco::RecoAlg_Cluster3DAlgs
  larreco::ClusterParamsImportWrapper
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Principal
  art_plugin_support::toolMaker
  art::Persistency_Common
  cetlib::cetlib
  ROOT::Physics
  ROOT::Tree
)

cet_build_plugin(ClusterAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardataobj::RecoBase
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(ClusterCheater art::EDProducer
  LIBRARIES PRIVATE
  larreco::ClusterFinder
  larreco::ClusterParamsImportWrapper
  larreco::RecoAlg_ClusterRecoUtil
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::Utilities
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RecoBase
  nug4::ParticleNavigation
  art::Framework_Services_Registry
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(ClusterCrawler art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  larcoreobj::SimpleTypesAndConstants
  lardataobj::RawData
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(ClusterPCA art::EDAnalyzer
  LIBRARIES PRIVATE
  larcoreobj::SimpleTypesAndConstants
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

cet_build_plugin(ClusterTrackAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
)

cet_build_plugin(ClusteringValidation art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Gpad
  ROOT::Hist
)

cet_build_plugin(DBCluster3D art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(DBclusterAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  larcorealg::Geometry
  lardataobj::RecoBase
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(DBcluster art::EDProducer
  LIBRARIES PRIVATE
  larreco::ClusterFinder
  larreco::RecoAlg
  larreco::RecoAlg_ClusterRecoUtil
  larreco::ClusterParamsImportWrapper
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  lardata::Utilities
  lardata::AssociationUtil
  lardata::DetectorClocksService
  larcorealg::Geometry
  lardataobj::RecoBase
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(HoughLineFinderAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(HoughLineFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardataobj::RecoBase
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(LineCluster art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::ArtDataHelper
  lardata::AssociationUtil
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(LineMerger art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg_ClusterRecoUtil
  larreco::ClusterParamsImportWrapper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::Utilities
  lardata::AssociationUtil
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(SmallClusterFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::ClusterFinder
  larreco::RecoAlg
  larreco::RecoAlg_ClusterRecoUtil
  larreco::ClusterParamsImportWrapper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::Utilities
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(TrajCluster art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larreco::RecoAlg_TCAlg
  larsim::MCCheater_ParticleInventoryService_service
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::Utilities
  lardata::AssociationUtil
  larcore::ServiceUtil
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
