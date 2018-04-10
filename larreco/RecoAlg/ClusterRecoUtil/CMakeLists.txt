art_make( 
          LIB_LIBRARIES
                        lardata_Utilities
                        larcorealg_Geometry
			${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
