cet_make_library(SOURCE
  MCBTAlgConstants.h
  MCBTAlg.cxx
  MCBTException.cxx
  MCMatchAlg.cxx
  LIBRARIES
  PUBLIC
  lardataobj::RecoBase
  lardataobj::Simulation
  canvas::canvas
  PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::DetectorClocksService
  lardataalg::DetectorInfo
  art::Framework_Services_Registry
  canvas::canvas
  ROOT::Core
)

cet_build_plugin(MCBTDemo art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::MCComp
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardataobj::RecoBase
  lardataobj::MCBase
  lardataobj::Simulation
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
)

install_headers()
install_fhicl()
install_source()
