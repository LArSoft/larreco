# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_test(GausFitCache_test USE_BOOST_UNIT
                           LIBRARIES larreco_RecoAlg
        )

cet_test(VoronoiDiagram_test LIBRARIES larreco_RecoAlg_Cluster3DAlgs_Voronoi
                                       larreco_RecoAlg_Cluster3DAlgs)
