cet_enable_asserts()

cet_make_library(LIBRARY_NAME Voronoi_IEvent INTERFACE
  SOURCE IEvent.h
  LIBRARIES INTERFACE
  larreco::Voronoi_DCEL
)

cet_make_library(LIBRARY_NAME Voronoi_SweepEvent INTERFACE
  SOURCE SweepEvent.h
  LIBRARIES INTERFACE
  larreco::Voronoi_IEvent
  larreco::Voronoi_DCEL
)

cet_make_library(SOURCE
  BeachLine.cxx
  EventUtilities.cxx
  Voronoi.cxx
  LIBRARIES
  PUBLIC
  larreco::Voronoi_IEvent
  larreco::Voronoi_DCEL
  Boost::boost
  Eigen3::Eigen
  PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_ConvexHull
)

install_headers()
install_fhicl()
install_source()
