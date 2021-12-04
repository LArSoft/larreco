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
         )

simple_plugin(FeatureVertexFinder "module"
                messagefacility::MF_MessageLogger
                art::Framework_Services_Registry
                ROOT::Core
                ROOT::Hist
                larcorealg_Geometry
                lardataalg_DetectorInfo
                lardataobj_RecoBase
              )

install_fhicl()
install_source()
