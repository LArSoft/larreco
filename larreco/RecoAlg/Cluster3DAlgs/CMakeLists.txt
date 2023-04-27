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

cet_make_library(LIBRARY_NAME ClusterAlg INTERFACE
  SOURCE IClusterAlg.h
  LIBRARIES INTERFACE larreco::RecoAlg_Cluster3DAlgs)

cet_write_plugin_builder(lar::ClusterAlg art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larreco::ClusterAlg)

cet_make_library(LIBRARY_NAME ClusterModAlg INTERFACE
  SOURCE IClusterModAlg.h
  LIBRARIES INTERFACE larreco::RecoAlg_Cluster3DAlgs)

cet_write_plugin_builder(lar::ClusterModAlg art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larreco::ClusterModAlg art::Utilities)

cet_make_library(LIBRARY_NAME ClusterParamsBuilder INTERFACE
  SOURCE IClusterParamsBuilder.h
  LIBRARIES INTERFACE larreco::RecoAlg_Cluster3DAlgs)

cet_write_plugin_builder(lar::ClusterParamsBuilder art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larreco::ClusterParamsBuilder)

cet_make_library(LIBRARY_NAME Hit3DBuilder INTERFACE
  SOURCE IHit3DBuilder.h
  LIBRARIES INTERFACE
  larreco::RecoAlg_Cluster3DAlgs
  lardataobj::RecoBase
  canvas::canvas
)

cet_write_plugin_builder(lar::Hit3DBuilder art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larreco::Hit3DBuilder)

cet_collect_plugin_builders(Modules lar::Cluster3DToolBuilders
  LIST
  lar::ClusterAlg
  lar::ClusterModAlg
  lar::ClusterParamsBuilder
  lar::Hit3DBuilder
)

include(lar::Cluster3DToolBuilders)

add_subdirectory(ConvexHull)
add_subdirectory(Voronoi)
add_subdirectory(PathFinding)

cet_build_plugin(ClusterMergeAlg lar::ClusterModAlg
  LIBRARIES PRIVATE
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Core
  ROOT::Hist
)

cet_build_plugin(ClusterParamsBuilder lar::ClusterParamsBuilder
  LIBRARIES PRIVATE
  larcoreobj::SimpleTypesAndConstants
)

cet_build_plugin(DBScanAlg lar::ClusterAlg
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  art::Framework_Services_Registry
  art_plugin_support::toolMaker
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
)

cet_build_plugin(MinSpanTreeAlg lar::ClusterAlg
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Services_Registry
  art_plugin_support::toolMaker
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Physics
  Eigen3::Eigen
)

cet_build_plugin(SnippetHit3DBuilder lar::Hit3DBuilder
  LIBRARIES PRIVATE
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::Tree
  Eigen3::Eigen
)

cet_build_plugin(SpacePointHit3DBuilder lar::Hit3DBuilder
  LIBRARIES PRIVATE
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Tree
)

cet_build_plugin(StandardHit3DBuilder lar::Hit3DBuilder
  LIBRARIES PRIVATE
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
