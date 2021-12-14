include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make(LIB_LIBRARIES
           ROOT::Core
           ROOT::Hist
         MODULE_LIBRARIES
           larcorealg_Geometry
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           larreco_QuadVtx
           nusimdata::SimulationBase
           lardata_ArtDataHelper
           art_root_io::TFileService_service
           art_root_io::tfile_support
           art::Framework_Services_Registry
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Hist
           ROOT::Matrix
           ROOT::Physics
           ROOT::Tree
        )

install_headers()
install_fhicl()
install_source()
