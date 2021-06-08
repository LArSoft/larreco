cet_enable_asserts()

art_make(TOOL_LIBRARIES
           larreco_RecoAlg_Cluster3DAlgs
           larreco_RecoAlg_Cluster3DAlgs_Voronoi
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ROOT::Core
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           canvas
           messagefacility::MF_MessageLogger
        )

install_headers()
install_fhicl()
install_source()
