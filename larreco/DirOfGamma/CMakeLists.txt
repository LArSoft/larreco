
art_make( 
          
          LIB_LIBRARIES larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation
                        lardata_RawData
                        lardata_RecoBase
		        larreco_RecoAlg
			larreco_RecoAlg_PMAlg
                        lardata_Utilities
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
 			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
		        ${Boost_SYSTEM_LIBRARY}
	  MODULE_LIBRARIES larreco_RecoAlg
			   larreco_RecoAlg_PMAlg
			   larreco_DirOfGamma
			   larsim_MCCheater_BackTracker_service
      			   lardata_RecoBase
			   larsim_Simulation
			   larcore_Geometry
	    		   larcore_Geometry_Geometry_service
      			   lardata_Utilities
			   nusimdata_SimulationBase
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   art_Persistency_Common canvas_Persistency_Common
			   art_Persistency_Provenance canvas_Persistency_Provenance
			   art_Utilities canvas_Utilities
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	      		   ${FHICLCPP}
			   ${CETLIB}
			   ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
