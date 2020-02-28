include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(HitReaders)

art_make(BASENAME_ONLY
         LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ROOT_CORE}
           ${ROOT_PHYSICS}
           canvas
           cetlib_except
         MODULE_LIBRARIES
           larcorealg_Geometry
           lardataobj_RecoBase
           larsim_MCCheater_BackTrackerService_service
           lardata_ArtDataHelper
           larreco_SpacePointSolver
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ROOT_CORE}
           ${ROOT_HIST}
        )

install_headers()
install_fhicl()
install_source()
