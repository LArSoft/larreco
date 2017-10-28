
add_subdirectory(ShowerReco3D)
add_subdirectory(job)

art_make( 
#          EXCLUDE ShowerSelectorFilter_module.cc
          MODULE_LIBRARIES
			larreco_RecoAlg
                        larsim_MCCheater_BackTracker_service
                        lardataobj_AnalysisBase
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larreco_Calorimetry
                        lardataobj_AnalysisBase
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardata_Utilities
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

#simple_plugin( ShowerSelectorFilter "module"
#         		  larevt_Filters
# 			  larreco_RecoAlg
# 		    	  lardataobj_RawData
# 			  nusimdata_SimulationBase
#                
# 	     )     

install_headers()
install_fhicl()
install_source()
