cet_build_plugin(CellTree art::EDAnalyzer
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larcorealg::GeoAlgo
  lardataalg::DetectorInfo
  lardataobj::RawData
  lardataobj::RecoBase
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  canvas::canvas
  ROOT::Core
  ROOT::EG
  ROOT::Hist
  ROOT::Physics
  ROOT::RIO
  ROOT::Tree
)

install_fhicl()
install_source()
