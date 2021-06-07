# ======================================================================
#  larreco main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larreco/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(larreco LANGUAGES CXX)

find_package(cetmodules REQUIRED)
find_package(canvas_root_io)
find_package(art)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( nutools )
find_ups_product( nug4 )
find_ups_product( nurandom )
find_ups_product( art )
find_ups_product(art_root_io)
find_ups_product( postgresql )
find_ups_product( range )
find_ups_product( eigen )
find_ups_root()

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(larreco)

# tests
add_subdirectory(test)
add_subdirectory(cmakeConfig)
# ups - table and config files
#add_subdirectory(ups)

# packaging utility
#include(UseCPack)
