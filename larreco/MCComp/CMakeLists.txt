#add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           art::Framework_Services_Registry
           canvas::canvas
           cetlib_except::cetlib_except
           ROOT::Core
           art::Utilities
         MODULE_LIBRARIES
           larreco_MCComp
           larcorealg_Geometry
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           ROOT::Core
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

# install_headers()
install_fhicl()
install_source()
