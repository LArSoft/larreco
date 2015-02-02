
add_subdirectory(job)

art_make(
	  LIBRARY_NAME  ShowerReco3D
          LIB_LIBRARIES
                        RecoBase
                        Utilities
			AnalysisAlg
                        ClusterRecoUtil
			CMToolBase
	  		CMTAlgPriority
			CMTAlgMerge
			CMTAlgMatch
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_MINUIT}
          MODULE_LIBRARIES 
	  		   ShowerReco3D
	  		   RecoAlg
                           BackTracker_service
                           Filters
                           RecoBase
                           Geometry
                           Geometry_service
                           ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   ${ART_PERSISTENCY_COMMON}
			   ${ART_UTILITIES}
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	                   ${FHICLCPP}
			   ${CETLIB}
         )

install_headers()
install_fhicl()
install_source()

