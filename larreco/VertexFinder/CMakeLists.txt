cet_build_plugin(AggregateVertexAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Hist
)

cet_build_plugin(FeatureVertexFinderAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Hist
)

cet_build_plugin(FeatureVertexFinder art::EDProducer
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Hist
)

cet_build_plugin(PrimaryVertexFinder art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Hist
  ROOT::Physics
)

cet_build_plugin(VertexCheater art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  lardata::AssociationUtil
  lardataobj::RecoBase
  nug4::ParticleNavigation
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Hist
  ROOT::Physics
)

cet_build_plugin(VertexFinder2D art::EDProducer
  LIBRARIES PRIVATE
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Hist
  ROOT::MathCore
)

cet_build_plugin(VertexFitter art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::types
  canvas::canvas
)

install_fhicl()
install_source()
