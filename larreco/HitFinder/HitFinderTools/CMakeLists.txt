cet_enable_asserts()

add_library (gshf_MarqFitAlg STATIC MarqFitAlg.cxx)
set_property(TARGET gshf_MarqFitAlg PROPERTY POSITION_INDEPENDENT_CODE ON)

set( hitfinder_tool_lib_list
                        larreco_RecoAlg
			larcorealg_Geometry
			lardataobj_RecoBase
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nurandom_RandomUtils_NuRandomService_service
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common
			art_Persistency_Provenance
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		#	gshf_MarqFitAlg
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			
			${FHICLCPP}
			${CETLIB}
			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}
            ${CLHEP}
			${ROOT_BASIC_LIB_LIST}
    )

#set(${CMAKE_CXX_FLAGS}
#		${CMAKE_CXXFLAGS} 
#		-O3 
#		-fopenmp		
#  )

include(FindOpenMP)
#if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")


art_make(TOOL_LIBRARIES ${hitfinder_tool_lib_list}
			gshf_MarqFitAlg
        )

install_headers()
install_source()
install_fhicl()
