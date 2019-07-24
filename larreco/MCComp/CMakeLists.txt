#add_subdirectory(job)

art_make( 
	  
          LIB_LIBRARIES    larcorealg_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			   lardataobj_RecoBase
			   lardata_Utilities
			   ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
			   canvas
                           ${MF_MESSAGELOGGER}
                           
                           ${FHICLCPP}
                           cetlib cetlib_except
			   ${Boost_SYSTEM_LIBRARY}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larreco_MCComp
			   larcorealg_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			   ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
			   canvas
                           ${MF_MESSAGELOGGER}
                           
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
			)

# install_headers()
install_fhicl()
install_source()
