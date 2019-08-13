art_make(LIB_LIBRARIES
           larreco_RecoAlg_CMTool_CMToolBase
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ROOT_CORE}
           ${ROOT_RIO}
           ${ROOT_TREE}
           canvas
           cetlib_except
        )

install_headers()
install_fhicl()
install_source()
