art_make(LIB_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           art::Framework_Services_Registry
           art::Utilities
           ${CANVAS}
           ${CETLIB_EXCEPT}
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Matrix
           ROOT::Physics
           ROOT::Minuit2
        )

install_headers()
install_fhicl()
install_source()
