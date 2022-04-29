cet_enable_asserts()

cet_make_plugin_builder(lar::PathFinderTool art::tool Modules
  LIBRARIES CONDITIONAL
  larreco::RecoAlg_Cluster3DAlgs
)

include(lar::PathFinderTool)

cet_build_plugin(ClusterPathFinder lar::PathFinderTool
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
  larreco::RecoAlg_Cluster3DAlgs_Voronoi
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  Eigen3::Eigen
)

cet_build_plugin(ConvexHullPathFinder lar::PathFinderTool
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
  art_root_io::tfile_support
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  Eigen3::Eigen
)

cet_build_plugin(MSTPathFinder lar::PathFinderTool
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardata::AssociationUtil
  lardataobj::RecoBase
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  Eigen3::Eigen
)

cet_build_plugin(VoronoiPathFinder art::tool
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
  larreco::RecoAlg_Cluster3DAlgs_Voronoi
  art_root_io::tfile_support
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  Eigen3::Eigen
)

install_headers()
install_fhicl()
install_source()
