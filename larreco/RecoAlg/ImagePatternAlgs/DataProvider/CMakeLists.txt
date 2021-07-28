art_make(LIB_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           art::Framework_Services_Registry
           ROOT::Core
           canvas::canvas
           messagefacility::MF_MessageLogger
           cetlib_except::cetlib_except
           ${CLHEP}
        )

install_headers()
install_fhicl()
install_source()
