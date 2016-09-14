
art_make( 
          LIB_LIBRARIES lardataobj_RecoBase
                        larsim_Simulation lardataobj_Simulation
                        larevt_Filters
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
                        lardataobj_AnalysisBase
                        lardata_Utilities
           		${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            		${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	        	${FHICLCPP}
			${CETLIB}
            		${CLHEP}
			${ROOT_BASIC_LIB_LIST}
            		${ROOT_MINUIT}
			${ROOT_MINUIT2}	
			${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()
