add_subdirectory(Cheating)

art_make(TOOL_LIBRARIES 
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			larreco_RecoAlg
			${ART_FRAMEWORK_CORE}  
)

install_headers()
install_source()
install_fhicl()
