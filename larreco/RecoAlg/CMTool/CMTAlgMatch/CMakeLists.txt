cet_make_library(SOURCE
  CFAlgoQRatio.cxx
  CFAlgoShowerCompat.cxx
  CFAlgoTimeOverlap.cxx
  CFAlgoZOverlap.cxx
  LIBRARIES
  PUBLIC
  larreco::RecoAlg_CMTool_CMToolBase
  PRIVATE
  ROOT::RIO
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
