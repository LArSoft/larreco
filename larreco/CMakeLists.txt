add_subdirectory(Calorimetry)

# FIXME: TrackFinder/TrackMaker.h should be moved to its own
#        subdirectory to avoid circular dependency issues.
#        2022-04-15 CHG
cet_make_library(LIBRARY_NAME TrackMaker INTERFACE
  SOURCE TrackFinder/TrackMaker.h
  LIBRARIES INTERFACE
  lardataobj::RecoBase
  lardataobj::TrackingTypes
  canvas::canvas
)
cet_write_plugin_builder(lar::TrackMakerTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larreco::TrackMaker
)

add_subdirectory(RecoAlg)

add_subdirectory(Calibrator)
add_subdirectory(ClusterFinder)
add_subdirectory(DirOfGamma)
add_subdirectory(EventFinder)
add_subdirectory(Genfit)
add_subdirectory(HitFinder)
add_subdirectory(MCComp)
add_subdirectory(QuadVtx)
add_subdirectory(SCECorrections)
add_subdirectory(ShowerFinder)
add_subdirectory(SpacePointFinder)
add_subdirectory(SpacePointSolver)
add_subdirectory(TrackFinder)
add_subdirectory(VertexFinder)
add_subdirectory(WireCell)
