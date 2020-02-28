#add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           canvas
           cetlib_except
           ${ROOT_CORE}
         MODULE_LIBRARIES
           larreco_MCComp
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ROOT_CORE}
         )

# install_headers()
install_fhicl()
install_source()
