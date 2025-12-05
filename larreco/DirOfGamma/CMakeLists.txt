cet_make_library(SOURCE DirOfGamma.cxx
  LIBRARIES
  PUBLIC
  lardataobj::RecoBase
  canvas::canvas
  ROOT::Hist
  PRIVATE
  larreco::RecoAlg_PMAlg
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataalg::DetectorInfo
  messagefacility::MF_MessageLogger
  ROOT::MathCore
)

cet_build_plugin(EMShower3D art::EDProducer
  LIBRARIES PRIVATE
  larreco::DirOfGamma
  larreco::RecoAlg
  larreco::RecoAlg_PMAlg
  larreco::Calorimetry
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  larcore::ServiceUtil
  lardataalg::UtilitiesHeaders
  lardataobj::RecoBase
  art::Framework_Services_Registry
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(MultiEMShowers art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::RecoAlg_PMAlg
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Physics
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
