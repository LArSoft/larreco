add_subdirectory(HitFinderTools)

art_make(LIB_LIBRARIES
           lardataobj_RawData
           lardataobj_RecoBase
           larcorealg_Geometry
           messagefacility::MFMessageLogger
           fhiclcpp::fhiclcpp fhiclcpp::types
           cetlib_except
           ROOT::Core
           ${ART_UTILITIES}
         MODULE_LIBRARIES
           larreco_HitFinder
           larsim_MCCheater_BackTrackerService_service
           larreco_RecoAlg
           lardataobj_RecoBase
           larcorealg_Geometry
           lardata_ArtDataHelper
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ROOT::Core
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           art_Framework_Services_System_TriggerNamesService_service
           messagefacility::MFMessageLogger
           ROOT::Hist
           ROOT::MathCore
           ROOT::Physics
           ${TBB}
         )

install_headers()
install_fhicl()
install_source()
