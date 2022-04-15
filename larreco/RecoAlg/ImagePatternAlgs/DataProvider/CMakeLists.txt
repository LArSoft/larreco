art_make(LIB_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           art::Framework_Services_Registry
           ROOT::Core
           canvas::canvas
           messagefacility::MF_MessageLogger
           cetlib_except::cetlib_except
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
