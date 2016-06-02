#add_subdirectory(job)

art_make( 
	  
          LIB_LIBRARIES    larcore_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation
			   lardata_RecoBase
			   lardata_Utilities
			   ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common canvas_Persistency_Common
                           art_Persistency_Provenance canvas_Persistency_Provenance
                           art_Utilities canvas_Utilities
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
			   ${Boost_SYSTEM_LIBRARY}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larreco_MCComp
			   larcore_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation
			   ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common canvas_Persistency_Common
                           art_Persistency_Provenance canvas_Persistency_Provenance
                           art_Utilities canvas_Utilities
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
			)

# install_headers()
install_fhicl()
install_source()
