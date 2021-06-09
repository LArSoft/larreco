art_make(MODULE_LIBRARIES
  lardataalg_DetectorInfo
  larcorealg_Geometry
  lardataobj_RecoBase
  art::Framework_Services_Registry
  ${ROOT_CORE}
  ${ROOT_PHYSICS}
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  )

install_headers()
install_fhicl()
install_source()
