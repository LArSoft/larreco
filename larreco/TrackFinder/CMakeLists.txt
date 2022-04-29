include(lar::TrackMakerTool)

cet_build_plugin(KalmanFilterFitTrackMaker lar::TrackMakerTool
  LIBRARIES PRIVATE
  larreco::TrackCreationBookKeeper
  larreco::TrackMaker
  larreco::RecoAlg
  lardata::RecoObjects
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  fhiclcpp::types
)

install_headers()
install_fhicl()
install_source()
