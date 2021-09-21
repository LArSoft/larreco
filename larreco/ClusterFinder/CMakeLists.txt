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
           art::Framework_Services_Registry
           art_root_io::tfile_support
           ROOT::Core
           ROOT::Physics
           art_root_io::TFileService_service
           messagefacility::MF_MessageLogger
         )

install_headers()
install_fhicl()
install_source()

add_subdirectory(RStarTree)
add_subdirectory(TrackShowerSplitter)
