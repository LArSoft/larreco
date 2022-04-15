add_subdirectory(job)

art_make(LIB_LIBRARIES
           lardata_Utilities
           larreco_Calorimetry
           larreco_RecoAlg_CMTool_CMToolBase
           larreco_RecoAlg_CMTool_CMTAlgPriority
           larreco_RecoAlg_CMTool_CMTAlgMatch
           art::Framework_Services_Registry
           ROOT::Core
           canvas::canvas
           cetlib_except::cetlib_except
         MODULE_LIBRARIES
           larreco_MCComp
           larreco_ShowerFinder_ShowerReco3D
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           ROOT::Core
           ROOT::Physics
           messagefacility::MF_MessageLogger
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
