#add_subdirectory(job)

art_make(LIB_LIBRARIES
           larcorealg_Geometry
           art::Framework_Services_Registry
           canvas::canvas
           cetlib_except::cetlib_except
           ROOT::Core
           ${ART_UTILITIES}
         MODULE_LIBRARIES
           larreco_MCComp
           larcorealg_Geometry
           art::Framework_Services_Registry
           ROOT::Core
         )

# install_headers()
install_fhicl()
install_source()
