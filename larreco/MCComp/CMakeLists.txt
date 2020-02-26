#add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           canvas
           cetlib_except
           ROOT::Core
           ${ART_UTILITIES}
         MODULE_LIBRARIES
           larreco_MCComp
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
         )

# install_headers()
install_fhicl()
install_source()
