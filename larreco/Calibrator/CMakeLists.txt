#include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( LIB_LIBRARIES
		        lardataobj_RecoBase 
		        larcoreobj_SimpleTypesAndConstants
		        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common
			art_Utilities
                        art_Persistency_Provenance
			canvas
                        ${MF_MESSAGELOGGER}
		        
		        ${FHICLCPP}
		        cetlib cetlib_except
		        ${CLHEP}
                        ${Boost_SYSTEM_LIBRARY}
		        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
                        lardataobj_RecoBase
		        larcoreobj_SimpleTypesAndConstants
		        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
