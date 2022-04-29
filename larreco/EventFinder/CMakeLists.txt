cet_build_plugin(EventCheater art::EDProducer
  LIBRARIES PRIVATE
  lardata::AssociationUtil
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
