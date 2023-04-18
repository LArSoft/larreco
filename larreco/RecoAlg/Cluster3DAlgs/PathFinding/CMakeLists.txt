cet_enable_asserts()

cet_build_plugin(ClusterPathFinder lar::ClusterModAlg
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
  larreco::RecoAlg_Cluster3DAlgs_Voronoi
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  Eigen3::Eigen
)

cet_build_plugin(ConvexHullPathFinder lar::ClusterModAlg
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
  art_root_io::tfile_support
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  Eigen3::Eigen
)

cet_build_plugin(MSTPathFinder lar::ClusterModAlg
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

cet_build_plugin(VoronoiPathFinder lar::ClusterModAlg
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
