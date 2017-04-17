
art_make( 
          LIB_LIBRARIES lardataobj_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
                        art_Persistency_Common
                        art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_source()
