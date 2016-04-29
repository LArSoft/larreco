
simple_plugin(PointIdTrainingData "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardata_RecoBase
			lardata_AnalysisAlg
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			${SIMULATIONBASE}
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

simple_plugin(PointIdEffTest "module"
			larcore_Geometry_Geometry_service
			larcore_Geometry
			lardata_RecoBase
			lardata_AnalysisAlg
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			${SIMULATIONBASE}
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

add_subdirectory(PointIdAlg)

install_headers()
install_fhicl()
install_source()
