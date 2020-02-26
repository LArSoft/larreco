art_make(LIB_LIBRARIES
           larreco_RecoAlg_PMAlg
           ROOT::Core
           ROOT::Physics
           canvas
           ${MF_MESSAGELOGGER}
           ${FHICLCPP}
           cetlib_except
        )

install_headers()
install_fhicl()
install_source()
