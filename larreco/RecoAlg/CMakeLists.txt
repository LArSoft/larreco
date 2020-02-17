art_make(LIB_LIBRARIES
           larcorealg_Geometry
           larreco_RecoAlg_TCAlg
           larreco_RecoAlg_PMAlg
           larreco_RecoAlg_ClusterRecoUtil
           larreco_RecoAlg_CMTool_CMToolBase
           larreco_RecoAlg_ImagePatternAlgs_DataProvider
           ROOT::Core
           ROOT::Physics
           ROOT::Matrix
           ROOT::Minuit
           canvas
           ${FHICLCPP}
           cetlib_except
        )

add_subdirectory(CMTool)
add_subdirectory(ClusterRecoUtil)
add_subdirectory(Cluster3DAlgs)
add_subdirectory(ImagePatternAlgs)
add_subdirectory(PMAlg)
add_subdirectory(TCAlg)
add_subdirectory(xml)

install_headers()
install_fhicl()
install_source()
