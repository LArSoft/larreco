
art_make( BASENAME_ONLY
	  MODULE_LIBRARIES RecoAlg
			   PMAlg
			   DirOfGamma
			   BackTracker_service
      			   RecoBase
			   Simulation
			   Geometry
	    		   Geometry_service
      			   Utilities
			   DetectorProperties_service
			   LArProperties_service  
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

add_subdirectory(DirOfGamma)

install_headers()
install_fhicl()
install_source()
