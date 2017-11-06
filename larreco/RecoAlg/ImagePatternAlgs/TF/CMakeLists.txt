art_make( 
          LIB_LIBRARIES 
			${FHICLCPP}
			cetlib cetlib_except
			${TENSORFLOW}
			${PROTOBUF}
			rt # this is required for glibc earlier than 2.17 (i.e. SLF6)
        )

install_headers()
install_source()

