add_subdirectory(ShowerReco3D)
add_subdirectory(job)

art_make(MODULE_LIBRARIES
          larreco_RecoAlg
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
