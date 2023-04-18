add_subdirectory(Segmentation2D)

cet_build_plugin(TrackShowerHits art::EDProducer
  LIBRARIES PRIVATE
  larreco::ClusterFinder_TrackShowerSplitter_Segmentation2D
  larcore::Geometry_Geometry_service
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()
install_scripts()
