art_make(LIB_LIBRARIES
           larreco_RecoAlg_PMAlg
           ROOT::Core
           ROOT::Physics
           canvas::canvas
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           cetlib_except::cetlib_except
        )

install_headers()
install_fhicl()
install_source()
