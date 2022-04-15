cet_enable_asserts()

art_make(TOOL_LIBRARIES
           canvas::canvas
           messagefacility::MF_MessageLogger
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
