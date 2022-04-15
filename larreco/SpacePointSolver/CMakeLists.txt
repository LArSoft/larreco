include_directories(${CMAKE_CURRENT_SOURCE_DIR}) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible
add_subdirectory(HitReaders)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           art::Framework_Services_Registry
           art::Utilities
           canvas::canvas
           cetlib_except::cetlib_except
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
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
