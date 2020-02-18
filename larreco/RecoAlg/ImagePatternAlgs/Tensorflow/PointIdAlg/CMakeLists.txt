include_directories($ENV{TRTIS_CLIENTS_INC})
cet_find_library(TRTIS_CLIENTS_LIBRARY NAMES request PATHS $ENV{TRTIS_CLIENTS_LIB})

art_make(
          LIB_LIBRARIES larreco_RecoAlg_ImagePatternAlgs_DataProvider 
			larreco_RecoAlg_ImagePatternAlgs_Keras
                        larreco_RecoAlg_ImagePatternAlgs_Tensorflow_TF
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${TRTIS_CLIENTS_LIBRARY}
			
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
