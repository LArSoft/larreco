
art_make( LIBRARY_NAME  RecoAlg
          LIB_LIBRARIES larsim_MCCheater_BackTracker_service
                        RecoBase
                        RecoBaseArt
                        Simulation
                        Filters
                        Geometry
                        larcore_Geometry_Geometry_service
                        RecoObjects
                        AnalysisBase
			AnalysisAlg
                        Utilities
                        ClusterRecoUtil
			CMToolBase
			CMTAlgMerge
			CMTAlgPriority
			CMTAlgMatch
			CMToolApp
			PMAlg
                        ${SIMULATIONBASE}
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
			${ROOT_MINUIT2}	
        )

add_subdirectory(CMTool)
add_subdirectory(ClusterRecoUtil)
add_subdirectory(Cluster3DAlgs)
add_subdirectory(PMAlg)

install_headers()
install_fhicl()
install_source()

