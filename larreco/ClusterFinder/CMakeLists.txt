art_make(LIB_LIBRARIES
           lardataobj_RecoBase
         MODULE_LIBRARIES
           larreco_Calorimetry
           larreco_RecoAlg_Cluster3DAlgs
           larreco_RecoAlg_TCAlg
           larreco_RecoAlg
           larreco_ClusterFinder
           larsim_MCCheater_ParticleInventoryService_service
           lardataobj_AnalysisBase
           lardataobj_RecoBase
           nurandom_RandomUtils_NuRandomService_service
           nusimdata_SimulationBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ROOT::Core
           ROOT::Physics
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           messagefacility::MF_MessageLogger
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(RStarTree)
add_subdirectory(TrackShowerSplitter)
