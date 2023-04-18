cet_build_plugin(CellTree art::EDAnalyzer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  larcorealg::Geometry
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
