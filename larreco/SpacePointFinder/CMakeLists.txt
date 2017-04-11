
include_directories ( $ENV{G4INCLUDE} )

art_make(  
          MODULE_LIBRARIES  larreco_RecoAlg
			    lardata_RecoObjects larreco_Deprecated
			    lardataobj_RecoBase
			    larcore_Geometry
			    larcore_Geometry_Geometry_service
			    lardata_Utilities
			    nutools_MagneticField_MagneticField_service
			    nusimdata_SimulationBase
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    art_Persistency_Common
			    art_Utilities
			    canvas
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
			    ${FHICLCPP}
			    cetlib cetlib_except
			    ${CLHEP}
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
         )
install_headers()
install_fhicl()
install_source()


