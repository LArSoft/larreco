include_directories ($ENV{G4INCLUDE}) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make(EXCLUDE
           Track3DKalman_module.cc
           Track3DKalmanSPS_module.cc
           DumpTracks_module.cc
           KalmanFilterFitTrackMaker_tool.cc
         MODULE_LIBRARIES
           larreco_RecoAlg
           nug4::MagneticFieldServices_MagneticFieldServiceStandard_service
           ROOT::Core
           messagefacility::MF_MessageLogger
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(Track3DKalman "module"
                larreco_Genfit
                lardataobj_RecoBase
                larcorealg_Geometry
                nurandom::RandomUtils_NuRandomService_service
                nusimdata::SimulationBase
                art::Framework_Services_Registry
                art_root_io::TFileService_service
                art_root_io::tfile_support
                ROOT::Core
                ROOT::Matrix
                ROOT::Physics
                ROOT::Tree
                messagefacility::MF_MessageLogger
              ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(Track3DKalmanSPS "module"
                larreco_Genfit
                lardataobj_RecoBase
                larcorealg_Geometry
                nusimdata::SimulationBase
                art::Framework_Services_Registry
                art_root_io::TFileService_service
                art_root_io::tfile_support
                messagefacility::MF_MessageLogger
                ROOT::Core
                ROOT::Hist
                ROOT::Matrix
                ROOT::Physics
                ROOT::Tree
              ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(KalmanFilterFitTrackMaker "tool"
                canvas::canvas
                lardataobj_RecoBase
                larreco_RecoAlg
              ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
