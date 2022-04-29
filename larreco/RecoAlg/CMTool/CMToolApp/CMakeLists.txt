cet_make_library(SOURCE CMergeHelper.cxx
  LIBRARIES PUBLIC
  larreco::RecoAlg_CMTool_CMToolBase
)

install_headers()
install_fhicl()
install_source()
