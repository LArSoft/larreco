include_directories(${CMAKE_CURRENT_SOURCE_DIR}) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make(LIB_LIBRARIES
           ROOT::Core
           ROOT::Hist
         MODULE_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           larreco_QuadVtx
           nusimdata::SimulationBase
           lardata_ArtDataHelper
           art_root_io::TFileService_service
           art_root_io::tfile_support
           art::Framework_Services_Registry
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Hist
           ROOT::Matrix
           ROOT::Physics
           ROOT::Tree
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
