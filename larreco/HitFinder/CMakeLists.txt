add_subdirectory(HitFinderTools)

art_make(LIB_LIBRARIES
           lardataobj_RawData
           lardataobj_RecoBase
           larcorealg_Geometry
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           cetlib_except::cetlib_except
           ROOT::Core
           art::Utilities
         MODULE_LIBRARIES
           larreco_HitFinder
           larsim_MCCheater_BackTrackerService_service
           larreco_RecoAlg
           lardataobj_RecoBase
           larcorealg_Geometry
           lardata_ArtDataHelper
           art::Framework_Services_Registry
           art_root_io::tfile_support
           ROOT::Core
           art_root_io::TFileService_service
           art_Framework_Services_System_TriggerNamesService_service
           messagefacility::MF_MessageLogger
           ROOT::Hist
           ROOT::MathCore
           ROOT::Physics
           ${TBB}
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
