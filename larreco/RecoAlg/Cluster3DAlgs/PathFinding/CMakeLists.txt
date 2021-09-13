cet_enable_asserts()

art_make(TOOL_LIBRARIES
           larreco_RecoAlg_Cluster3DAlgs
           larreco_RecoAlg_Cluster3DAlgs_Voronoi
           art_root_io::tfile_support
           ROOT::Core
           art_root_io::TFileService_service
           canvas::canvas
           messagefacility::MF_MessageLogger
        )

install_headers()
install_fhicl()
install_source()
