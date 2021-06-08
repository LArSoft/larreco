add_subdirectory(ShowerReco3D)
add_subdirectory(job)

art_make(MODULE_LIBRARIES
          larreco_RecoAlg
          lardataobj_RecoBase
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ROOT::Core
          ROOT::Physics
          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
          ${MF_MESSAGELOGGER}
        )

install_headers()
install_fhicl()
install_source()
