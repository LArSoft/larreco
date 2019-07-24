
art_make( 
          LIB_LIBRARIES larsim_MCCheater_BackTrackerService_service
	                larsim_MCCheater_ParticleInventoryService_service
                        lardataobj_RecoBase
                        lardata_ArtDataHelper
                        larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
                        larevt_Filters
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
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
			larreco_RecoAlg_ImagePatternAlgs_DataProvider
			larreco_RecoAlg_PMAlg
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
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
add_subdirectory(xml)

install_headers()
install_fhicl()
install_source()

