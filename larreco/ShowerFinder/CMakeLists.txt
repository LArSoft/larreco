add_subdirectory(ShowerReco3D)
add_subdirectory(job)

cet_build_plugin(EMShower art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::ArtDataHelper
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  range-v3::range-v3
  ROOT::Physics
)

cet_build_plugin(NeutrinoShowerEff art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardata::ArtDataHelper
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(NuShowerEff art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardata::ArtDataHelper
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(ShowerCheater art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(ShowerFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::MathCore
)

cet_build_plugin(ShowerReco art::EDProducer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::ArtDataHelper
  lardata::AssociationUtil
  lardata::Utilities
  larcorealg::Geometry
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::MathCore
  ROOT::Tree
  range-v3::range-v3
)

cet_build_plugin(TCShowerAnalysis art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardataobj::RecoBase
  lardataobj::Simulation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(TCShowerElectronLikelihood art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(TCShowerTemplateMaker art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::MCBase
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(TCShower art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  range-v3::range-v3
)

install_headers()
install_fhicl()
install_source()
