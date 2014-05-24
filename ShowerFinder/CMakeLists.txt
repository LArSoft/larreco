
art_make( BASENAME_ONLY
#          EXCLUDE ShowerSelectorFilter_module.cc
          MODULE_LIBRARIES 
                        LArProperties_service
                        DetectorProperties_service
                        BackTracker_service
                        AnalysisBase
                        RecoBase
                        Geometry
                        Geometry_service
                        AnalysisAlg
                        AnalysisBase
                        Simulation
                        Utilities
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
                        ${ROOT_BASIC_LIB_LIST}
         )

#simple_plugin( ShowerSelectorFilter "module"
#         		  Filters
# 			  RecoAlg
# 		    	  RawData
#                          LArProperties_service
# 			  ${SIMULATIONBASE}
#               BASENAME_ONLY 
# 	     )     

install_headers()
install_fhicl()
install_source()
