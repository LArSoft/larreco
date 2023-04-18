include(lar::TrackMakerTool)

cet_build_plugin(CCTrackMaker art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(CosmicTracker art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::MathCore
  ROOT::Physics
)

cet_build_plugin(FeatureTracker art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Physics
)

cet_build_plugin(KalmanFilterFinalTrackFitter art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larreco::TrackMaker
  lardata::DetectorPropertiesService
  lardata::RecoObjects
  lardataobj::AnalysisBase
  lardataobj::MCBase
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::types
  fhiclcpp::fhiclcpp
)

cet_build_plugin(KalmanFilterFitTrackMaker lar::TrackMakerTool
  LIBRARIES PRIVATE
  larreco::TrackCreationBookKeeper
  larreco::TrackMaker
  larreco::RecoAlg
  lardata::RecoObjects
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  fhiclcpp::types
)

cet_build_plugin(KalmanFilterTrajectoryFitter art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larreco::TrackMaker
  lardata::DetectorPropertiesService
  lardata::RecoObjects
  lardataobj::MCBase
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MCSFitProducer art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardataobj::RecoBase
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::types
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MagDriftAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  ######################################################################
  # FIXME: there should be an INTERFACE target, e.g. nug4::MagneticFieldService
  nug4::MagneticField
  Geant4::G4global
  ######################################################################
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Graf
  ROOT::Hist
)

cet_build_plugin(MuonTrackingEff art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::RecoBase
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Graf
  ROOT::Hist
  ROOT::Physics
)

cet_build_plugin(NeutrinoTrackingEff art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::RecoBase
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(PFPAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
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

cet_build_plugin(PMAlgTrackMaker art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::ArtDataHelper
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  larcorealg::Geometry
  art_root_io::TFileService_service
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(PMAlgTrajFitter art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
)

cet_build_plugin(SeedAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  ROOT::Hist
  ROOT::Matrix
)

cet_build_plugin(SeedFinderModule art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
)

cet_build_plugin(SpacePointAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larsim::MCCheater_BackTrackerService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  ROOT::Hist
  range-v3::range-v3
)

cet_build_plugin(SpacePointCheater art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
)

cet_build_plugin(SpacePointFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
)

cet_build_plugin(SpacePts art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::MathCore
  ROOT::Physics
)

cet_build_plugin(TCTrack art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(Track3DKalman art::EDProducer
  LIBRARIES PRIVATE
  larreco::Genfit
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  CLHEP::Random
  ROOT::Matrix
  ROOT::Physics
  ROOT::Tree
)

cet_build_plugin(Track3DKalmanHit art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::RecoObjects
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(Track3DKalmanSPS art::EDProducer
  LIBRARIES PRIVATE
  larreco::Genfit
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::EG
  ROOT::Hist
  ROOT::MathCore
  ROOT::Matrix
  ROOT::Physics
  ROOT::Tree
)

cet_build_plugin(Track3Dreco art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::MathCore
  ROOT::Physics
)

cet_build_plugin(TrackAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::RecoBase
  lardataobj::Simulation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  ROOT::Hist
  ROOT::Matrix
)

cet_build_plugin(TrackCheater art::EDProducer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Physics
)

cet_build_plugin(TrackKalmanCheater art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larsim::MCCheater_ParticleInventoryService_service
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  lardata::RecoObjects
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Hist
)

cet_build_plugin(TrackProducerFromPFParticle art::EDProducer
  LIBRARIES PRIVATE
  larreco::TrackMaker
  lardata::DetectorPropertiesService
  lardata::ForEachAssociatedGroup
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(TrackProducerFromTrack art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorPropertiesService
  lardata::ForEachAssociatedGroup
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(TrackProducerFromTrackTrajectory art::EDProducer
  LIBRARIES PRIVATE
  larreco::TrackMaker
  lardata::DetectorPropertiesService
  lardata::ForEachAssociatedGroup
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  art::Framework_Principal
  art::Persistency_Common
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

cet_build_plugin(TrackStitcher art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
