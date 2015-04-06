# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_test(GausFitCache_test USE_BOOST_UNIT
                           LIBRARIES RecoAlg
        )
