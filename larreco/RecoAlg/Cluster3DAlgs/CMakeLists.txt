cet_enable_asserts()

# FIXME: Voronoi/DCEL.h should be moved to its own subdirectory to avoid
#        circular dependency issues.
#        2022-04-15 CHG
cet_make_library(LIBRARY_NAME Voronoi_DCEL INTERFACE
  SOURCE Voronoi/DCEL.h
  LIBRARIES INTERFACE
  Eigen3::Eigen
)

cet_make_library(SOURCE
  Cluster3D.cxx
  HoughSeedFinderAlg.cxx
  PCASeedFinderAlg.cxx
  ParallelHitsSeedFinderAlg.cxx
  PrincipalComponentsAlg.cxx
  SeedFinderAlgBase.h
  SkeletonAlg.cxx
  kdTree.cxx
  LIBRARIES
  PUBLIC
  larreco::Voronoi_DCEL
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  ROOT::Gpad
  Eigen3::Eigen
  PRIVATE
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataalg::DetectorInfo
  larcorealg::Geometry
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::Matrix
  ROOT::Physics
)

add_subdirectory(ConvexHull)
add_subdirectory(Voronoi)
add_subdirectory(PathFinding)


install_headers()
install_fhicl()
install_source()
