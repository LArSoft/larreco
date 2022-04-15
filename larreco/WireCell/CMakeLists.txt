art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           lardataobj_RawData
           lardataobj_RecoBase
           nusimdata::SimulationBase
           art::Framework_Services_Registry
           ROOT::Core
           ROOT::RIO
           ROOT::Tree
           messagefacility::MF_MessageLogger
           ROOT::EG
           ROOT::Hist
           ROOT::Physics
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

# install_headers()
install_fhicl()
install_source()
