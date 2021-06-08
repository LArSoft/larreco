include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make(LIB_LIBRARIES
           ROOT::Core
           ROOT::Hist
         MODULE_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           larreco_QuadVtx
           nusimdata_SimulationBase
           lardata_ArtDataHelper
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           messagefacility::MFMessageLogger
           ROOT::Core
           ROOT::Hist
           ROOT::Matrix
           ROOT::Physics
           ROOT::Tree
        )

install_headers()
install_fhicl()
install_source()
