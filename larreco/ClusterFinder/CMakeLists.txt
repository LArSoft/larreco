
art_make( 
          
          LIB_LIBRARIES larcore_Geometry
                        lardataobj_RecoBase
                        larreco_RecoAlg_ClusterRecoUtil
          MODULE_LIBRARIES larreco_RecoAlg
                           larreco_RecoAlg_Cluster3DAlgs
                           larreco_RecoAlg_ClusterRecoUtil
                           larreco_ClusterFinder
                           larsim_MCCheater_BackTracker_service
                           larevt_Filters
                           lardataobj_RecoBase
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           nutools_RandomUtils_NuRandomService_service
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common canvas_Persistency_Common
                           art_Utilities canvas_Utilities
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(RStarTree)
add_subdirectory(TrackShowerSplitter)
