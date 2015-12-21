art_make_library( LIBRARY_NAME RFFHitFitter
                             SOURCE RegionAboveThresholdFinder.cxx
					   GaussianEliminationAlg.cxx
					   RFFHitFitter.cxx
                             LIBRARIES
                                                     ${ROOT_BASIC_LIB_LIST}
						     ${Boost_SYSTEM_LIBRARY}
)

art_make( 
	  LIBRARY_NAME  HitFinder
          LIB_LIBRARIES
			RawData
			RecoBase 
			RFFHitFitter
                        Geometry
                        larcore_Geometry_Geometry_service
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
	  		HitFinder
                        larsim_MCCheater_BackTracker_service
			RecoAlg
                        RecoBase
			MCBase
                        Geometry
                        larcore_Geometry_Geometry_service
                        Utilities
                        Simulation
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

install_headers()
install_fhicl()
install_source()
