art_make(EXCLUDE FeatureVertexFinder_module.cc
         MODULE_LIBRARIES
           larsim_MCCheater_ParticleInventoryService_service
           lardataobj_RecoBase
           larcorealg_Geometry
           larreco_RecoAlg
           nusimdata::SimulationBase
           lardata_RecoObjects
           art::Framework_Services_Registry
           art_root_io::TFileService_service
           art_root_io::tfile_support
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           messagefacility::MF_MessageLogger
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(FeatureVertexFinder "module"
                messagefacility::MF_MessageLogger
                art::Framework_Services_Registry
                ROOT::Core
                ROOT::Hist
                larcorealg_Geometry
                lardataalg_DetectorInfo
                lardataobj_RecoBase
              ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

install_fhicl()
install_source()
