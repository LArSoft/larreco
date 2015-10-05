include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  DirOfGamma
          LIB_LIBRARIES DetectorProperties_service
                        LArProperties_service
                        Geometry
                        Geometry_service
                        Simulation
                        RawData
                        RecoBase
		        RecoAlg
			PMAlg
                        Utilities
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
 			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
		        ${Boost_SYSTEM_LIBRARY}
          )

install_headers()
install_fhicl()
install_source()
