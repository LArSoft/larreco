#include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( LIB_LIBRARIES	larevt_Filters
		        lardataobj_RecoBase
		        lardata_RecoBaseArt
		        larcore_Geometry
		        larcore_Geometry_Geometry_service
                        lardataobj_AnalysisBase
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common canvas_Persistency_Common
			art_Utilities canvas_Utilities
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        ${MF_MESSAGELOGGER}
		        ${MF_UTILITIES}
		        ${FHICLCPP}
		        ${CETLIB}
		        ${CLHEP}
                        ${Boost_SYSTEM_LIBRARY}
		        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
                        larevt_Filters
                        lardata_RecoObjects larreco_Deprecated
                        lardataobj_RecoBase
		        lardata_RecoBaseArt
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_AnalysisBase
                        lardata_Utilities
		        larreco_RecoAlg
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
