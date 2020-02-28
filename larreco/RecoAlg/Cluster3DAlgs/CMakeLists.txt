cet_enable_asserts()

add_subdirectory(ConvexHull)
add_subdirectory(Voronoi)
add_subdirectory(PathFinding)

art_make(LIB_LIBRARIES
           lardataobj_RecoBase
           larcorealg_Geometry
           lardata_ArtDataHelper
           lardata_RecoObjects
           ${ART_FRAMEWORK_CORE}
           ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_ROOT_IO_TFILE_SUPPORT}
           ${ART_UTILITIES}
           ROOT::Core
           ROOT::Gpad
           ROOT::Hist
           ROOT::Matrix
           ROOT::Physics
           ROOT::Tree
           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
           art_Persistency_Provenance
           canvas
           ${MF_MESSAGELOGGER}
           ${FHICLCPP}
           ${CETLIB}
           cetlib_except
          TOOL_LIBRARIES larreco_RecoAlg_Cluster3DAlgs
        )

install_headers()
install_fhicl()
install_source()
