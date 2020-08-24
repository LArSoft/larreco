include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(HitReaders)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_UTILITIES}
           ${CANVAS}
           ${CETLIB_EXCEPT}
           ROOT::Core
           ROOT::Physics
         MODULE_LIBRARIES
           larcorealg_Geometry
           lardataobj_RecoBase
           larsim_MCCheater_BackTrackerService_service
           lardata_ArtDataHelper
           larreco_SpacePointSolver
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
           ROOT::Hist
        )

install_headers()
install_fhicl()
install_source()
