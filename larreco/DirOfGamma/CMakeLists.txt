
art_make( 
          LIBRARY_NAME  DirOfGamma
          LIB_LIBRARIES Geometry
                        larcore_Geometry_Geometry_service
                        Simulation
                        RawData
                        RecoBase
		        RecoAlg
			PMAlg
                        Utilities
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
 			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
		        ${Boost_SYSTEM_LIBRARY}
	  MODULE_LIBRARIES RecoAlg
			   PMAlg
			   DirOfGamma
			   larsim_MCCheater_BackTracker_service
      			   RecoBase
			   Simulation
			   Geometry
	    		   larcore_Geometry_Geometry_service
      			   Utilities
			   ${SIMULATIONBASE}
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   ${ART_PERSISTENCY_COMMON}
			   ${ART_PERSISTENCY_PROVENANCE}
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
