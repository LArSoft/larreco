# ======================================================================
#  larreco main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larreco/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larreco VERSION 09.08.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package(Range-v3 REQUIRED PUBLIC)
find_package (Eigen3 3.3...3.3.9 REQUIRED NO_MODULE)
find_package(messagefacility REQUIRED PUBLIC)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package(nusimdata REQUIRED PUBLIC)
find_package(nug4 REQUIRED PUBLIC)
find_package(nurandom REQUIRED PUBLIC)
find_package(nutools REQUIRED PUBLIC)
find_package(PostgreSQL REQUIRED PUBLIC)
find_package(Boost COMPONENTS filesystem REQUIRED PUBLIC)
find_package(ROOT COMPONENTS Core RIO Matrix Physics Tree Hist EG REQUIRED PUBLIC)
find_package(CLHEP COMPONENTS Random REQUIRED PUBLIC)

find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

include_directories( $ENV{EIGEN_INC} )

# source
add_subdirectory(larreco)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
