art_make(LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
           ROOT::Matrix
           ROOT::Physics
           canvas
           ${MF_MESSAGELOGGER}
           cetlib_except
           ROOT::Minuit2
           ${ART_UTILITIES}
        )

install_headers()
install_fhicl()
install_source()
