install_fhicl()

# copy in the source tree
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )
