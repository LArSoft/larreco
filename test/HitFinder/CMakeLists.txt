# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_test(HitAnaAlg_test USE_BOOST_UNIT
                        LIBRARIES larreco_HitFinder
)
