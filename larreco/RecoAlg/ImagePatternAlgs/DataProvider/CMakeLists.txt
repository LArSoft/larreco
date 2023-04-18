cet_make_library(SOURCE DataProviderAlg.cxx
  LIBRARIES
  PUBLIC
  larreco::Calorimetry
  lardataobj::RecoBase
  fhiclcpp::types
  cetlib_except::cetlib_except
  CLHEP::Random # For testing on noise, not reco.
  PRIVATE
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()
