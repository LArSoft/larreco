cet_enable_asserts()

cet_make_library(LIBRARY_NAME HitReaderTool INTERFACE
  SOURCE IHitReader.h
  LIBRARIES INTERFACE
  lardataobj::RecoBase
  canvas::canvas
)

cet_write_plugin_builder(lar::HitReaderTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larreco::HitReaderTool
)

include(lar::HitReaderTool)

cet_build_plugin(HitsICARUS lar::HitReaderTool
  LIBRARIES PRIVATE
  messagefacility::MF_MessageLogger
)

cet_build_plugin(HitsStandard lar::HitReaderTool
  LIBRARIES PRIVATE
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()
