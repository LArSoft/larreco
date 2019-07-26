art_make(MODULE_LIBRARIES
           larcorealg_Geometry
           lardataobj_RawData
           lardataobj_RecoBase
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ROOT_CORE}
           ${ROOT_RIO}
           ${ROOT_TREE}
           ${MF_MESSAGELOGGER}
           ${ROOT_EG}
           ${ROOT_HIST}
           ${ROOT_PHYSICS}
        )

# install_headers()
install_fhicl()
install_source()
