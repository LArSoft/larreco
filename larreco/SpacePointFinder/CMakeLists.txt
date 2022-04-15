include_directories ( $ENV{G4INCLUDE} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make(MODULE_LIBRARIES
           larreco_RecoAlg
           messagefacility::MF_MessageLogger
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types
install_headers()
install_fhicl()
install_source()
