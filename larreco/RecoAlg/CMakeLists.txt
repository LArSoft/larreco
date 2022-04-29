add_subdirectory(ClusterRecoUtil)
add_subdirectory(CMTool)
add_subdirectory(Cluster3DAlgs)
add_subdirectory(ImagePatternAlgs)
add_subdirectory(PMAlg)
add_subdirectory(TCAlg)
add_subdirectory(xml)

cet_make_library(LIBRARY_NAME ClusterParamsImportWrapper INTERFACE
  SOURCE ClusterParamsImportWrapper.h
  LIBRARIES INTERFACE
  lardata::headers
)

cet_make_library(LIBRARY_NAME RootMathFunctor INTERFACE
  SOURCE RootMathFunctor.h
  LIBRARIES INTERFACE
  larcorealg::Math_Functor
)

cet_make_library(LIBRARY_NAME Track3DKalmanHit INTERFACE
  SOURCE Track3DKalmanHit.h
  LIBRARIES INTERFACE
  lardata::RecoObjects
  lardataobj::RecoBase
  canvas::canvas
)

cet_make_library(LIBRARY_NAME TrackTrajectoryCreationBookKeeper INTERFACE
  SOURCE TrackTrajectoryCreationBookKeeper.h
  LIBRARIES INTERFACE
  canvas::canvas
  lardataobj::RecoBase
  lardataobj::TrackingTypes
)

cet_make_library(LIBRARY_NAME TrackCreationBookKeeper INTERFACE
  SOURCE TrackCreationBookKeeper.h
  LIBRARIES INTERFACE
  larreco::TrackTrajectoryCreationBookKeeper
  larreco::TrackMaker
  lardataobj::RecoBase
  lardataobj::TrackFitHitInfo
  lardataobj::TrackingTypes
)

cet_make_library(LIBRARY_NAME VertexFitMinuitStruct INTERFACE
  SOURCE VertexFitMinuitStruct.h
  LIBRARIES INTERFACE
  ROOT::Physics
)

cet_make_library(LIBRARY_NAME VertexWrapper INTERFACE
  SOURCE VertexWrapper.h
  LIBRARIES INTERFACE
  lardataobj::RecoBase
)

cet_make_library(SOURCE
  APAGeometryAlg.cxx
  BlurredClusteringAlg.cxx
  CCHitFinderAlg.cxx
  ClusterCrawlerAlg.cxx
  ClusterMatchAlg.cxx
  ClusterMatchTQ.cxx
  ClusterMergeHelper.cxx
  CornerFinderAlg.cxx
  CosmicTrackerAlg.cxx
  DBScan3DAlg.cxx
  DBScanAlg.cxx
  DisambigAlg.cxx
  EMShowerAlg.cxx
  EndPointAlg.cxx
  GausFitCache.cxx
  Geometric3DVertexFitter.cxx
  HoughBaseAlg.cxx
  KalmanFilterAlg.cxx
  LinFitAlg.cxx
  MergeClusterAlg.cxx
  PMAlgCosmicTagger.cxx
  PMAlgStitching.cxx
  PMAlgTracking.cxx
  PMAlgVertexing.cxx
  ProjectionMatchingAlg.cxx
  SeedFinderAlgorithm.cxx
  ShowerEnergyAlg.cxx
  SmallClusterFinderAlg.cxx
  SpacePointAlg.cxx
  SpacePointAlg_TimeSort.cxx
  StitchAlg.cxx
  TCShowerAlg.cxx
  Track3DKalmanHitAlg.cxx
  TrackKalmanFitter.cxx
  TrackLineFitAlg.cxx
  TrackMomentumCalculator.cxx
  TrackShowerSeparationAlg.cxx
  TrackTrajectoryAlg.cxx
  TrajClusterAlg.cxx
  TrajectoryMCSFitter.cxx
  VertexFitAlg.cxx
  LIBRARIES
  PUBLIC
  larreco::RecoAlg_CMTool_CMToolBase
  larreco::RecoAlg_ClusterRecoUtil
  larreco::RecoAlg_ImagePatternAlgs_DataProvider
  larreco::RecoAlg_PMAlg
  larreco::RecoAlg_TCAlg
  larreco::VertexWrapper
  larreco::VertexFitMinuitStruct
  larreco::Calorimetry
  larsim::MCCheater_BackTrackerService_service
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  lardata::RecoObjects
  lardata::Utilities
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Core
  ROOT::Gpad
  ROOT::Hist
  ROOT::MathCore
  ROOT::Physics
  ROOT::TMVA
  Boost::headers
  RStarTree::RStarTree
  PRIVATE
  larreco::RecoAlg_ImagePatternAlgs_DataProvider
  larreco::TrackMaker
  larreco::TrackCreationBookKeeper
  larevt::ChannelStatusProvider
  larcore::ServiceUtil
  larcorealg::Geometry
  lardata::ArtDataHelper
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardata::RecoObjects
  lardata::SimpleFits
  lardataalg::DetectorInfo
  lardataalg::UtilitiesHeaders
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib::container_algorithms
  range-v3::range-v3
  ROOT::GenVector
  ROOT::Graf
  ROOT::Matrix
  ROOT::Minuit
  ROOT::Minuit2
  ROOT::RIO
  ROOT::Tree
  CLHEP::Random
)

install_headers()
install_fhicl()
install_source()
