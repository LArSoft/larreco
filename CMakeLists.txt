# ======================================================================
#  larreco main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larreco/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19) ### MIGRATE-ACTION-COMPATIBILITY (migrate-3.10.00) - Uncomment to preserve compatibility with older CMake versions

find_package(cetmodules REQUIRED)
project(larreco VERSION 09.11.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - add args: REPORT_THRESHOLD VERBOSE

find_package(Range-v3 REQUIRED EXPORT)
find_package (Eigen3 3.3...3.3.9 REQUIRED NO_MODULE)
find_package(messagefacility REQUIRED EXPORT)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED EXPORT)
find_package(art_root_io REQUIRED EXPORT)
find_package(nusimdata REQUIRED EXPORT)
find_package(nug4 REQUIRED EXPORT)
find_package(nurandom REQUIRED EXPORT)
find_package(nutools REQUIRED EXPORT)
find_package(PostgreSQL REQUIRED EXPORT)
find_package(Boost COMPONENTS filesystem REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core RIO Matrix Physics Tree Hist EG REQUIRED EXPORT)
find_package(CLHEP COMPONENTS Random REQUIRED EXPORT)

find_package( larcoreobj REQUIRED )
find_package( larcorealg REQUIRED )
find_package( larcore REQUIRED )
find_package( lardataobj REQUIRED )
find_package( lardataalg REQUIRED )
find_package( lardata REQUIRED )
find_package( larevt REQUIRED )
find_package( larsim REQUIRED )
find_package( larvecutils REQUIRED )
find_package( RStarTree REQUIRED )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

include_directories( $ENV{EIGEN_INC} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.10.00) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

# source
add_subdirectory(larreco)

# tests
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
