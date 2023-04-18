cet_build_plugin(SCECorrection art::EDProducer
  LIBRARIES PRIVATE
  larevt::SpaceChargeService
  larevt::SpaceCharge
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  larcorealg::Geometry
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Core
  ROOT::Physics
)

install_headers()
install_fhicl()
install_source()
