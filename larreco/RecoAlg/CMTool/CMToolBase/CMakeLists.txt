cet_make_library(SOURCE
  CMManagerBase.cxx
  CMatchBookKeeper.cxx
  CMatchManager.cxx
  CMergeBookKeeper.cxx
  CMergeManager.cxx
  LIBRARIES
  PUBLIC
  larreco::RecoAlg_ClusterRecoUtil
  lardata::headers
  ROOT::Core
)

install_headers()
install_fhicl()
install_source()
