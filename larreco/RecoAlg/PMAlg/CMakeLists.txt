art_make(LIB_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_UTILITIES}
           ${CANVAS}
           ${CETLIB_EXCEPT}
           ${MF_MESSAGELOGGER}
           ROOT::Core
           ROOT::Matrix
           ROOT::Physics
           ROOT::Minuit2
        )

install_headers()
install_fhicl()
install_source()
