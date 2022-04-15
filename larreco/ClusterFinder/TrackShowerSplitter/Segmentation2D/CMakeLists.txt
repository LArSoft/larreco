art_make(LIB_LIBRARIES
           larreco_RecoAlg_PMAlg
           ROOT::Core
           ROOT::Physics
           canvas::canvas
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           cetlib_except::cetlib_except
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
