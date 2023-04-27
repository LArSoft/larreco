# ======================================================================
#  larreco main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larreco/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.12.00 REQUIRED)
project(larreco VERSION 09.19.00 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

cet_cmake_module_directories(Modules BINARY)

find_package(larcore REQUIRED EXPORT)
find_package(larcorealg REQUIRED EXPORT)
find_package(larcoreobj REQUIRED EXPORT)
find_package(lardata REQUIRED EXPORT)
find_package(lardataalg REQUIRED EXPORT)
find_package(lardataobj REQUIRED EXPORT)
find_package(larevt REQUIRED EXPORT)
find_package(larsim REQUIRED EXPORT)
find_package(larvecutils REQUIRED EXPORT)

find_package(nug4 REQUIRED EXPORT)
find_package(nurandom REQUIRED EXPORT)
find_package(nusimdata REQUIRED EXPORT)

find_package(art_root_io REQUIRED EXPORT)
find_package(canvas_root_io REQUIRED EXPORT)
find_package(art REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)

find_package(Boost COMPONENTS headers REQUIRED EXPORT)
find_package(CLHEP COMPONENTS Random REQUIRED EXPORT)
find_package(Eigen3 REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core EG GenVector Geom Gpad Graf Graf3d Hist MathCore Matrix Minuit Minuit2 Physics RIO TMVA Tree REQUIRED EXPORT)
find_package(Range-v3 REQUIRED EXPORT)
find_package(RStarTree REQUIRED EXPORT)
find_package(TBB REQUIRED EXPORT)

include(BuildDictionary)
include(BuildPlugins)

# source
add_subdirectory(Modules)
add_subdirectory(larreco)

# tests
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
