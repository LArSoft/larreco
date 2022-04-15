art_make(LIB_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           art::Framework_Services_Registry
           art::Utilities
           canvas::canvas
           cetlib_except::cetlib_except
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Matrix
           ROOT::Physics
           ROOT::Minuit2
           range-v3::range-v3
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
