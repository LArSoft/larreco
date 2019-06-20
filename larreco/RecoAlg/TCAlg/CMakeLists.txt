
art_make( 
      LIB_LIBRARIES larsim_MCCheater_BackTrackerService_service
	            larsim_MCCheater_ParticleInventoryService_service
      larsim_Simulation
      nutools_ParticleNavigation
      lardataobj_Simulation
      nusimdata_SimulationBase
			lardataobj_RecoBase
			lardata_RecoObjects
			larcorealg_Geometry
      larcore_Geometry_Geometry_service
      larreco_Calorimetry
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
      ${ROOT_BASIC_LIB_LIST}
      ${ROOT_GEOM}
			${Boost_SYSTEM_LIBRARY}
			${ROOT_GENVECTOR}
      ${ROOT_TMVA}
      art_Framework_Services_Optional
        )

add_subdirectory(TCDebugTools)

install_headers()
install_source()
