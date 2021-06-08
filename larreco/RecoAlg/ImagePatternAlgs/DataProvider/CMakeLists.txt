art_make(LIB_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
           canvas
           messagefacility::MFMessageLogger
           cetlib_except
           ${CLHEP_LIBRARIES}
        )

install_headers()
install_fhicl()
install_source()
