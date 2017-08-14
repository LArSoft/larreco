add_subdirectory(HitFinderTools)

art_make( 
          LIB_LIBRARIES
			lardataobj_RawData
			lardataobj_RecoBase 
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES 
	  		larreco_HitFinder
                        larsim_MCCheater_BackTracker_service
			larreco_RecoAlg
                        lardataobj_RecoBase
			lardataobj_MCBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
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
