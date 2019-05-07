
art_make( 
          LIB_LIBRARIES larreco_Calorimetry
                        lardataobj_AnalysisBase
                        lardataobj_RecoBase
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
        )


install_headers()
install_fhicl()
install_source()


