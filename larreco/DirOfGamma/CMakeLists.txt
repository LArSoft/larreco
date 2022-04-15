art_make(LIB_LIBRARIES
           larcorealg_Geometry
           larreco_RecoAlg_PMAlg
           art::Framework_Services_Registry
           ROOT::Core
           ROOT::Physics
           canvas::canvas
           messagefacility::MF_MessageLogger
           cetlib_except::cetlib_except
         MODULE_LIBRARIES
           larreco_RecoAlg
           larreco_RecoAlg_PMAlg
           larreco_Calorimetry
           larreco_DirOfGamma
           art::Framework_Services_Registry
           ROOT::Core
           ROOT::Physics
           messagefacility::MF_MessageLogger
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
