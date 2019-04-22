cet_enable_asserts()

add_subdirectory(ConvexHull)
add_subdirectory(Voronoi)
add_subdirectory(PathFinding)

art_make( LIB_LIBRARIES  lardataobj_RecoBase
                         larcorealg_Geometry
                         larcore_Geometry_Geometry_service
                         lardata_ArtDataHelper
                         lardata_RecoObjects
                         lardata_Utilities
                         larreco_RecoAlg_Cluster3DAlgs_ConvexHull
                         larreco_RecoAlg_Cluster3DAlgs_Voronoi
                         ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_BASIC}
                         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                         art_Persistency_Common
                         art_Persistency_Provenance
                         art_Utilities
                         canvas
                         ${MF_MESSAGELOGGER}

                         ${FHICLCPP}
                         ${CETLIB}
                         cetlib_except
                         ${Boost_FILESYSTEM_LIBRARY}
                         ${Boost_SYSTEM_LIBRARY}
                         ${CLHEP}
                         ${ROOT_BASIC_LIB_LIST}
                         ${ROOT_GEOM}
                         ${ROOT_MINUIT}
          TOOL_LIBRARIES larreco_RecoAlg_Cluster3DAlgs
        )

install_headers()
install_fhicl()
install_source()
