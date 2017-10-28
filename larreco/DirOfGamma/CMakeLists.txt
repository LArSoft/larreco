
art_make( 
          
          LIB_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
		        larreco_RecoAlg
			larreco_RecoAlg_PMAlg
                        lardata_Utilities
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        art_Persistency_Common
			art_Persistency_Provenance
                        art_Utilities
			canvas
 			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
		        ${Boost_SYSTEM_LIBRARY}
	  MODULE_LIBRARIES larreco_RecoAlg
			   larreco_RecoAlg_PMAlg
			   larreco_DirOfGamma
			   larsim_MCCheater_BackTracker_service
      			   lardataobj_RecoBase
			   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			   larcorealg_Geometry
	    		   larcore_Geometry_Geometry_service
      			   lardata_Utilities
			   nusimdata_SimulationBase
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   art_Persistency_Common
			   art_Persistency_Provenance
			   art_Utilities
			   canvas
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	      		   ${FHICLCPP}
			   cetlib cetlib_except
			   ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
