
cet_remove_compiler_flag( CXX -Werror )

art_make( LIBRARY_NAME  RecoAlg
          LIB_LIBRARIES HarrisVertexFinder_module
                        BackTracker_service
                        RecoBase
                        Simulation
                        Filters
                        Geometry
                        Geometry_service
                        RecoObjects
                        AnalysisBase
                        LArProperties_service
                        DetectorProperties_service
                        Utilities
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
        )

install_headers()
install_fhicl()
install_source()
