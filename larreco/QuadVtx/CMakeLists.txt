cet_make_library(SOURCE HeatMap.cxx
  LIBRARIES PRIVATE
  ROOT::Hist
)

cet_build_plugin(EvalVtx art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(QuadVtx art::EDProducer
  LIBRARIES PRIVATE
  larreco::QuadVtx
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::Matrix
  ROOT::Physics
)

install_headers()
install_fhicl()
install_source()
