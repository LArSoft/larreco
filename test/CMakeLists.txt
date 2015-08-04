# various integration tests

include(CetTest)

cet_enable_asserts()

install_fhicl()

add_subdirectory(RecoAlg)
add_subdirectory(HitFinder)
