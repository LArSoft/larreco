if (DEFINED ENV{TENSORFLOW_DIR})
set (LARRECO_RECOALG_IMAGEPATTERNALGS_TF larreco_RecoAlg_ImagePatternAlgs_TF)
else(DEFINED ENV{TENSORFLOW_DIR})
set (EXCLUDE_TF PointIdAlg.cxx)
endif (DEFINED ENV{TENSORFLOW_DIR})

art_make(
          EXCLUDE ${EXCLUDE_TF}
          LIB_LIBRARIES larreco_RecoAlg_ImagePatternAlgs_Keras
                        ${LARRECO_RECOALG_IMAGEPATTERNALGS_TF}
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			nusimdata_SimulationBase
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
			cetlib cetlib_except
			${CLHEP}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_MINUIT}
			${ROOT_MINUIT2}	
			${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()
