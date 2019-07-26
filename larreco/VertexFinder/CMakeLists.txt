art_make(EXCLUDE FeatureVertexFinder_module.cc
         MODULE_LIBRARIES
           larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
           larevt_Filters
           lardataobj_RecoBase
           larcorealg_Geometry
           larcore_Geometry_Geometry_service
           lardata_Utilities
           larsim_Simulation
           nug4_ParticleNavigation
           lardataobj_Simulation
           nusimdata_SimulationBase
           lardata_RecoObjects
           ${ART_FRAMEWORK_CORE}
           ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_FRAMEWORK_SERVICES_BASIC}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ROOT_CORE}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           art_Persistency_Common
           ${MF_MESSAGELOGGER}
           ${ROOT_BASIC_LIB_LIST}
           ${ROOT_GENVECTOR}
           larreco_RecoAlg
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
