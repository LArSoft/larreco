# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_test(GausFitCache_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larreco::RecoAlg
  cetlib::cetlib
  ROOT::Gpad
  ROOT::Hist  
)

cet_test(VoronoiDiagram_test
  LIBRARIES PRIVATE
  larreco::RecoAlg_Cluster3DAlgs_Voronoi
  larreco::RecoAlg_Cluster3DAlgs
  messagefacility::MF_MessageLogger
)
