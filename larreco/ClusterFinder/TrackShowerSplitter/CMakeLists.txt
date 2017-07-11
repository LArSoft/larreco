add_subdirectory(Segmentation2D)

simple_plugin(TrackShowerHits   "module"
			lardataobj_RawData
		 	lardataobj_RecoBase
			larreco_RecoAlg
			larreco_ClusterFinder_TrackShowerSplitter_Segmentation2D
		 	larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
               ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
         		${ART_FRAMEWORK_SERVICES_OPTIONAL}
    			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
              
)

install_headers()
install_fhicl()
install_source()
install_scripts()
