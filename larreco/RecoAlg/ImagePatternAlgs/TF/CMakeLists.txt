art_make( 
          LIB_LIBRARIES 
			${FHICLCPP}
			cetlib cetlib_except
			${TENSORFLOW}
			${PROTOBUF}
      rt #WORKAROUND FOR LARSOFT ISSUE 17960
        )

install_headers()
install_source()

