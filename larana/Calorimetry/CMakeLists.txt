include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( 
	  
	  LIB_LIBRARIES	larevt_Filters
		        lardata_RecoBase
		        larcore_Geometry
		        larcore_Geometry_Geometry_service
		        larreco_RecoAlg
                        lardata_AnalysisAlg
                        lardata_AnalysisBase
                        lardata_Utilities
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        ${CETLIB}
		        ${CLHEP}
                        ${Boost_SYSTEM_LIBRARY}
		        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
	  		larana_Calorimetry
                        larevt_Filters
                        lardata_RecoObjects
                        lardata_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation
                        lardata_RawData
                        lardata_AnalysisAlg
                        lardata_AnalysisBase
                        lardata_Utilities
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
