add_subdirectory(job)

art_make(LIB_LIBRARIES
           lardata_Utilities
           larreco_Calorimetry
           larreco_RecoAlg_CMTool_CMToolBase
           larreco_RecoAlg_CMTool_CMTAlgPriority
           larreco_RecoAlg_CMTool_CMTAlgMatch
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
           canvas
           cetlib_except
         MODULE_LIBRARIES
           larreco_MCComp
           larreco_ShowerFinder_ShowerReco3D
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ROOT::Core
           ROOT::Physics
           messagefacility::MFMessageLogger
         )

install_headers()
install_fhicl()
install_source()
