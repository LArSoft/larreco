art_make(LIB_LIBRARIES
           larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
           nug4::ParticleNavigation
           nusimdata::SimulationBase
           lardata_RecoObjects
           lardataalg_DetectorInfo
           larcorealg_Geometry
           larreco_Calorimetry
           art::Framework_Principal
           art::Framework_Services_Registry
           art::Persistency_Provenance
           art_root_io::TFileService_service
           art_root_io::tfile_support
           art::Utilities
           canvas::canvas
           cetlib_except::cetlib_except
           cetlib::cetlib
           ROOT::Core
           ROOT::GenVector
           ROOT::Gpad
           ROOT::Graf
           ROOT::Graf3d
           ROOT::Hist
           ROOT::MathCore
           ROOT::Matrix
           ROOT::RIO
           ROOT::Tree
           ROOT::TMVA) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

add_subdirectory(TCDebugTools)

install_headers()
install_source()
