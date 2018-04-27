simple_plugin(PointIdEffTest "module"
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

simple_plugin(PointIdTrainingData "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
			nusimdata_SimulationBase
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

simple_plugin(EmTrackClusterId2out "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
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

simple_plugin(EmTrackClusterId3out "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
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

simple_plugin(EmTrackMichelId "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
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

simple_plugin(ParticleDecayId "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
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

install_fhicl()
