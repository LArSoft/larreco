art_make(MODULE_LIBRARIES
  lardataalg_DetectorInfo
  larcorealg_Geometry
  lardataobj_RecoBase
  art::Framework_Services_Registry
  ROOT::Core
  ROOT::Physics
  art_root_io::TFileService_service
  messagefacility::MF_MessageLogger
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
