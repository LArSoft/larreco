art_make( MODULE_LIBRARIES
                larcore_Geometry_Geometry_service
		larcorealg_Geometry
		lardataobj_RecoBase
		lardata_ArtDataHelper
		larreco_Calorimetry
		larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
		nusimdata_SimulationBase
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		${ART_ROOT_IO_TFILESERVICE_SERVICE}
                art_Framework_Services_System_TriggerNamesService_service
 		art_Persistency_Common
		art_Persistency_Provenance
		art_Utilities
		canvas
		${MF_MESSAGELOGGER}	
		cetlib cetlib_except
		${ROOT_BASIC_LIB_LIST}
)

install_fhicl()
