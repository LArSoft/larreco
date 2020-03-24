art_make(MODULE_LIBRARIES
  larcorealg_Geometry
  lardataobj_RecoBase
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ROOT_CORE}
  ${ROOT_PHYSICS}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  )

install_headers()
install_fhicl()
install_source()
