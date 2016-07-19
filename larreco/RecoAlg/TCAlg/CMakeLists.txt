
art_make( 
          LIB_LIBRARIES lardataobj_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
#                        lardata_RecoObjects
      ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
      art_Persistency_Common canvas_Persistency_Common
      art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
      ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_source()
