include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(HitReaders)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           art::Framework_Services_Registry
           art::Utilities
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
           art::Framework_Services_Registry
           ROOT::Core
           ROOT::Hist
        )

install_headers()
install_fhicl()
install_source()
