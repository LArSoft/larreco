
art_make( 
          LIBRARY_NAME ClusterFinder
          LIB_LIBRARIES Geometry
                        RecoBase
                        ClusterRecoUtil
          MODULE_LIBRARIES RecoAlg
                           Cluster3DAlgs
                           ClusterRecoUtil
                           ClusterFinder
                           larsim_MCCheater_BackTracker_service
                           Filters
                           RecoBase
                           Geometry
                           larcore_Geometry_Geometry_service
                           SeedService_service # artextensions
                           ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_PERSISTENCY_COMMON}
                           ${ART_UTILITIES}
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
