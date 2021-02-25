include_directories ($ENV{G4INCLUDE})

art_make(EXCLUDE
           Track3DKalman_module.cc
           Track3DKalmanSPS_module.cc
           DumpTracks_module.cc
           KalmanFilterFitTrackMaker_tool.cc
         MODULE_LIBRARIES
           larreco_RecoAlg
           nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
           ROOT::Core
           ${MF_MESSAGELOGGER}
         )

simple_plugin(Track3DKalman "module"
                larreco_Genfit
                lardataobj_RecoBase
                larcorealg_Geometry
                nurandom_RandomUtils_NuRandomService_service
                nusimdata_SimulationBase
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                ${ART_ROOT_IO_TFILE_SUPPORT}
                ${CLHEP}
                ROOT::Core
                ROOT::Matrix
                ROOT::Physics
                ROOT::Tree
                ${MF_MESSAGELOGGER}
              )

simple_plugin(Track3DKalmanSPS "module"
                larreco_Genfit
                lardataobj_RecoBase
                larcorealg_Geometry
                nusimdata_SimulationBase
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                ${ART_ROOT_IO_TFILE_SUPPORT}
                ${MF_MESSAGELOGGER}
                ROOT::Core
                ROOT::Hist
                ROOT::Matrix
                ROOT::Physics
                ROOT::Tree
              )

simple_plugin(KalmanFilterFitTrackMaker "tool"
                canvas
                lardataobj_RecoBase
                larreco_RecoAlg
              )

install_headers()
install_fhicl()
install_source()
