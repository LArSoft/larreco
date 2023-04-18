cet_enable_asserts()

cet_make_library(SOURCE ConvexHull.cxx
  LIBRARIES
  PUBLIC
  larreco::RecoAlg_Cluster3DAlgs
  PRIVATE
  Boost::headers
  Eigen3::Eigen
)

install_headers()
install_fhicl()
install_source()
