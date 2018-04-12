art_make( 
          LIB_LIBRARIES 
			${FHICLCPP}
			cetlib cetlib_except
			${TENSORFLOW}
			${PROTOBUF}
			${FOR_TENSORFLOW}
      #rt #workaround for larsoft issue 17960
        )

install_headers()
install_source()

