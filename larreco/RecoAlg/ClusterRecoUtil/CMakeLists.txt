cet_make_library(LIBRARY_NAME ClusterParamsAlgBase INTERFACE
  SOURCE ClusterParamsAlgBase.h
  LIBRARIES INTERFACE
  lardataobj::RecoBase
)

cet_make_library(LIBRARY_NAME OverriddenClusterParamsAlg INTERFACE
  SOURCE OverriddenClusterParamsAlg.h
  LIBRARIES INTERFACE
  larreco::ClusterParamsAlgBase
)

cet_make_library(SOURCE
  ClusterParams.h
  ClusterParamsAlg.cxx
  LazyClusterParamsAlg.cxx
  Polygon2D.cxx
  StandardClusterParamsAlg.cxx
  LIBRARIES
  PUBLIC
  larreco::ClusterParamsAlgBase
  lardata::headers
  lardataobj::RecoBase
  PRIVATE
  lardata::Utilities
  lardata::SimpleFits
  lardataalg::UtilitiesHeaders
  larcoreobj::headers
  cetlib::cetlib
  ROOT::Core
  ROOT::Gpad
  ROOT::Hist
  ROOT::MathCore
)
install_headers()
install_fhicl()
install_source()
