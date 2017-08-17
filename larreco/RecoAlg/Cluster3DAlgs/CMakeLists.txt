
art_make( EXCLUDE        MinSpanTreeAlg_tool.cc
          LIB_LIBRARIES  lardataobj_RecoBase
                         larcorealg_Geometry
                         larcore_Geometry_Geometry_service
                         lardata_RecoObjects
                         lardata_Utilities
                         ${ART_FRAMEWORK_CORE}
			             ${ART_FRAMEWORK_PRINCIPAL}
			             ${ART_FRAMEWORK_SERVICES_REGISTRY}
			             ${ART_FRAMEWORK_SERVICES_BASIC}
			             ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			             ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			             art_Persistency_Common
			             art_Persistency_Provenance
			             art_Utilities
			             canvas
			             ${MF_MESSAGELOGGER}
			             ${MF_UTILITIES}
	                     ${FHICLCPP}
                         ${CETLIB}
                         cetlib_except
                         ${Boost_FILESYSTEM_LIBRARY}
                         ${Boost_SYSTEM_LIBRARY}
                         ${CLHEP}
			             ${ROOT_BASIC_LIB_LIST}
                         ${ROOT_MINUIT}
        )

simple_plugin( MinSpanTreeAlg  "tool"
                         larreco_RecoAlg_Cluster3DAlgs
             )

simple_plugin( DBScanAlg  "tool"
                         larreco_RecoAlg_Cluster3DAlgs
             )

install_headers()
install_fhicl()
install_source()
