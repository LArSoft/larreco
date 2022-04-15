art_make(LIB_LIBRARIES
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art::Utilities
           canvas::canvas
           messagefacility::MF_MessageLogger
           cetlib_except::cetlib_except
           ROOT::Geom
           ROOT::Graf3d
           ROOT::Core
           ROOT::MathCore
           ROOT::Matrix
           ROOT::Physics
           ROOT::EG
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
