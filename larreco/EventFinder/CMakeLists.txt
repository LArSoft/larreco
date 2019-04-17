
art_make( 
          MODULE_LIBRARIES 
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common
                        art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
         )

install_headers()
install_fhicl()
install_source()
