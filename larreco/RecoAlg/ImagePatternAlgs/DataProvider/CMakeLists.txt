art_make(LIB_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
           canvas
           ${MF_MESSAGELOGGER}
           cetlib_except
           ${CLHEP}
        )

install_headers()
install_fhicl()
install_source()
