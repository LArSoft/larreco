
simple_plugin(PointIdTrainingData "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
)

simple_plugin(PointIdTrainingNuevent "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			larreco_RecoAlg_PMAlg
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
)

simple_plugin(PointIdEffTest "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
)

simple_plugin(EmTrackClusterId "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
)

simple_plugin(ParticleDecayId "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardataobj_RecoBase
			lardata_ArtDataHelper
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
)

add_subdirectory(Keras)
add_subdirectory(MLP)
add_subdirectory(PointIdAlg)

install_headers()
install_fhicl()
install_source()
