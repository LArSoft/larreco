cet_make_library(SOURCE
  CPAlgoArray.cxx
  CPAlgoIgnoreTracks.cxx
  CPAlgoNHits.cxx
  CPAlgoPolyArea.cxx
  CPAlgoQSum.cxx
  LIBRARIES
  PUBLIC
  larreco::RecoAlg_CMTool_CMToolBase
  ROOT::Core
  PRIVATE
  ROOT::MathCore
)

install_headers()
install_fhicl()
install_source()
