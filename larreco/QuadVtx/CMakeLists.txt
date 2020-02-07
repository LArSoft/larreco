include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make(LIB_LIBRARIES
           ${ROOT_CORE}
           ${ROOT_HIST}
         MODULE_LIBRARIES
           larcorealg_Geometry
           lardataobj_RecoBase
           larreco_QuadVtx
           nusimdata_SimulationBase
           lardata_ArtDataHelper
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${MF_MESSAGELOGGER}
           ${ROOT_BASIC_LIB_LIST}
           ${ROOT_CORE}
           ${ROOT_HIST}
        )

install_headers()
install_fhicl()
install_source()
