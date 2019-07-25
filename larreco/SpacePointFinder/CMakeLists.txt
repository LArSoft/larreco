
include_directories ( $ENV{G4INCLUDE} )

art_make(  
          MODULE_LIBRARIES  larreco_RecoAlg
			    lardata_RecoObjects
			    lardataobj_RecoBase
			    larcorealg_Geometry
			    larcore_Geometry_Geometry_service
			    lardata_Utilities
			    nug4_MagneticField_MagneticField_service
			    nusimdata_SimulationBase
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    art_Persistency_Common
			    art_Utilities
			    canvas
			    ${MF_MESSAGELOGGER}
			    
			    ${FHICLCPP}
			    cetlib cetlib_except
			    ${CLHEP}
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
         )
install_headers()
install_fhicl()
install_source()


