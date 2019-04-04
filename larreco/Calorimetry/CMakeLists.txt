#include_directories ( ${PROJECT_SOURCE_DIR} )
art_make( LIB_LIBRARIES	larevt_Filters
		        lardataobj_RecoBase
		        lardata_ArtDataHelper
		        larcorealg_Geometry
		        larcore_Geometry_Geometry_service
                        lardataobj_AnalysisBase
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
                        larevt_Filters
                        lardata_RecoObjects larreco_Deprecated
                        lardataobj_RecoBase
		        lardata_ArtDataHelper
                        larcorealg_Geometry
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
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
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
