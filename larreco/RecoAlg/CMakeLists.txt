
art_make( 
          LIB_LIBRARIES larsim_MCCheater_BackTracker_service
                        lardataobj_RecoBase
                        lardata_ArtDataHelper
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larevt_Filters
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects larreco_Deprecated
                        lardataobj_AnalysisBase
			larreco_Calorimetry
                        lardata_Utilities
                        larreco_RecoAlg_ClusterRecoUtil
			larreco_RecoAlg_CMTool_CMToolBase
			larreco_RecoAlg_CMTool_CMTAlgMerge
			larreco_RecoAlg_CMTool_CMTAlgPriority
			larreco_RecoAlg_CMTool_CMTAlgMatch
			larreco_RecoAlg_CMTool_CMToolApp
			${ROOT_MATHCORE}
			${ROOT_GENVECTOR}
                        ${ROOT_GEOM}
			${ROOT_BASIC_LIB_LIST}
      larreco_RecoAlg_TCAlg
			larreco_RecoAlg_PMAlg
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${CLHEP}
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_MINUIT}
			${ROOT_MINUIT2}	
        )

add_subdirectory(CMTool)
add_subdirectory(ClusterRecoUtil)
add_subdirectory(Cluster3DAlgs)
add_subdirectory(ImagePatternAlgs)
add_subdirectory(PMAlg)
add_subdirectory(TCAlg)

install_headers()
install_fhicl()
install_source()

