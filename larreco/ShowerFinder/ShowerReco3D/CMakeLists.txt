add_subdirectory(job)

cet_make_library(SOURCE
  ShowerCalo.cxx
  ShowerRecoAlg.cxx
  ShowerRecoAlgBase.cxx
  ShowerRecoException.cxx
  ShowerRecoManager.cxx
  LIBRARIES
  PUBLIC
  lardata::headers
  lardataobj::RecoBase
  PRIVATE
  larreco::RecoAlg_CMTool_CMToolBase
  larreco::RecoAlg_ClusterRecoUtil
  larreco::Calorimetry
  lardata::Utilities
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Services_Registry
  ROOT::Core
  ROOT::Hist
  ROOT::MathCore
)

cet_build_plugin(ShowerQuality art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::MCComp
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(ShowerReco3D art::EDProducer
  LIBRARIES PRIVATE
  larreco::ShowerFinder_ShowerReco3D
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larreco::RecoAlg_CMTool_CMTAlgMatch
  larreco::RecoAlg_CMTool_CMTAlgPriority
  larreco::RecoAlg_CMTool_CMToolBase
  larreco::Calorimetry
  lardata::AssociationUtil
  lardata::Utilities
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
