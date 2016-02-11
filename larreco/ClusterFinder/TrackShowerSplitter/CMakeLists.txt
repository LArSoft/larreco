add_subdirectory(Segmentation2D)

simple_plugin(TrackShowerHits   "module"
			lardata_RawData
		 	lardata_RecoBase
			larreco_RecoAlg
			larreco_ClusterFinder_TrackShowerSplitter_Segmentation2D
		 	larcore_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			lardata_Utilities_DetectorProperties_service
			lardata_Utilities_LArProperties_service
               ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
         		${ART_FRAMEWORK_SERVICES_OPTIONAL}
    			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
              
)

install_headers()
install_fhicl()
install_source()
install_scripts()
