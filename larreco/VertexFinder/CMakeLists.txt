art_make(EXCLUDE FeatureVertexFinder_module.cc
         MODULE_LIBRARIES
           larsim_MCCheater_ParticleInventoryService_service
           lardataobj_RecoBase
           larcorealg_Geometry
           larreco_RecoAlg
           nusimdata_SimulationBase
           lardata_RecoObjects
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ROOT_CORE}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
           ${MF_MESSAGELOGGER}
         )


simple_plugin(FeatureVertexFinder "module"
                ${MF_MESSAGELOGGER}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ROOT_CORE}
                ${ROOT_HIST}
                larcorealg_Geometry
                lardataobj_RecoBase
              )

install_fhicl()
install_source()
