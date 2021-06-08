include_directories ( $ENV{G4INCLUDE} )

art_make(MODULE_LIBRARIES
           larreco_RecoAlg
           ${MF_MESSAGELOGGER}
         )
install_headers()
install_fhicl()
install_source()
