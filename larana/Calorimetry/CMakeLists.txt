include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( BASENAME_ONLY
	  LIBRARY_NAME Calorimetry
	  LIB_LIBRARIES	Filters
		        RecoBase
		        Geometry
		        Geometry_service
		        RecoAlg
                        AnalysisAlg
                        AnalysisBase
                        LArProperties_service
                        DetectorProperties_service
                        Utilities
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        ${CETLIB}
		        ${CLHEP}
                        ${Boost_SYSTEM_LIBRARY}
		        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
	  		Calorimetry
                        Filters
                        RecoObjects
                        RecoBase
                        Geometry
                        Geometry_service
                        Simulation
                        RawData
                        AnalysisAlg
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
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
