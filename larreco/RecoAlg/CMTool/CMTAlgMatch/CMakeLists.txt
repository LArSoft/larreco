art_make(LIB_LIBRARIES
           larreco_RecoAlg_CMTool_CMToolBase
           larcorealg_Geometry
           art::Framework_Services_Registry
           ${ART_UTILITIES}
           ROOT::Core
           ROOT::RIO
           ROOT::Tree
           canvas::canvas
           cetlib_except::cetlib_except
        )

install_headers()
install_fhicl()
install_source()
