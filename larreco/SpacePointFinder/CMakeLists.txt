cet_build_plugin(TTSpacePointFinder art::EDProducer
  LIBRARIES PRIVATE
  larreco::RecoAlg
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
)

install_headers()
install_fhicl()
install_source()
