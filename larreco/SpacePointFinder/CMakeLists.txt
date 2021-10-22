include_directories ( $ENV{G4INCLUDE} )

art_make(MODULE_LIBRARIES
           larreco_RecoAlg
           messagefacility::MF_MessageLogger
         )
install_headers()
install_fhicl()
install_source()
