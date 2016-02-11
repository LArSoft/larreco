
add_subdirectory(ShowerReco3D)
add_subdirectory(job)

art_make( 
#          EXCLUDE ShowerSelectorFilter_module.cc
          MODULE_LIBRARIES
			larreco_RecoAlg
                        lardata_Utilities_LArProperties_service
                        lardata_Utilities_DetectorProperties_service
                        larsim_MCCheater_BackTracker_service
                        lardata_AnalysisBase
                        lardata_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_AnalysisAlg
                        lardata_AnalysisBase
                        larsim_Simulation
                        lardata_Utilities
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
#         		  larevt_Filters
# 			  larreco_RecoAlg
# 		    	  lardata_RawData
#                          lardata_Utilities_LArProperties_service
# 			  ${SIMULATIONBASE}
#                
# 	     )     

install_headers()
install_fhicl()
install_source()
