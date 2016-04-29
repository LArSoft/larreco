
art_make( 
          LIB_LIBRARIES lardata_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
#                        lardata_RecoObjects
      ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
      ${ART_PERSISTENCY_COMMON}
      ${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
      ${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_source()
