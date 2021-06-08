art_make(LIB_LIBRARIES
           larreco_RecoAlg_PMAlg
           ROOT::Core
           ROOT::Physics
           canvas
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp fhiclcpp::types
           cetlib_except
        )

install_headers()
install_fhicl()
install_source()
