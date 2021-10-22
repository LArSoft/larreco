cet_enable_asserts()

set( hitfinder_tool_lib_list
                        larreco_RecoAlg
                        larcorealg_Geometry
                        lardataobj_RecoBase
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        nurandom_RandomUtils_NuRandomService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art_root_io::tfile_support ROOT::Core
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        art_root_io::TFileService_service
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        Boost::filesystem
    )

include(FindOpenMP)
if(OPENMP_FOUND)
  # even if OpenMP is found on a SLF6 machine, it cannot be used.
  # see the discussion in https://cdcvs.fnal.gov/redmine/issues/23001
  if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6.32.*" )
    message(STATUS "OpenMP is not supported on SLF6")
  else()
    add_definitions(-DWITH_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif(OPENMP_FOUND)

art_make(TOOL_LIBRARIES ${hitfinder_tool_lib_list} )

install_headers()
install_source()
install_fhicl()
