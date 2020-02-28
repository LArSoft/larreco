art_make(LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ROOT_CORE}
           ${ROOT_MATRIX}
           ${ROOT_PHYSICS}
           canvas
           ${MF_MESSAGELOGGER}
           cetlib_except
           ${ROOT_MINUIT2}
        )

install_headers()
install_fhicl()
install_source()
