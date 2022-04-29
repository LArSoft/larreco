add_subdirectory(HitReaders)

cet_make_library(SOURCE
  HashTuple.h
  QuadExpr.cxx
  Solver.cxx
  TripletFinder.cxx
  LIBRARIES
  PUBLIC
  canvas::canvas
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataalg::DetectorInfo
  art::Framework_Services_Registry
  ROOT::Physics
)

cet_build_plugin(PlotSpacePoints art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_BackTrackerService_service
  lardata::DetectorClocksService
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Core
  ROOT::Hist
)

cet_build_plugin(SpacePointSolver art::EDProducer
  LIBRARIES PRIVATE
  larreco::SpacePointSolver
  larreco::HitReaderTool
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  lardata::ArtDataHelper
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_plugin_support::toolMaker
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
)


install_headers()
install_fhicl()
install_source()
