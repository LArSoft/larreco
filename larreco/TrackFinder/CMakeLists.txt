
include_directories ( $ENV{G4INCLUDE} )

art_make(  
          EXCLUDE Track3DKalman_module.cc Track3DKalmanSPS_module.cc DumpTracks_module.cc KalmanFilterFitTrackMaker_tool.cc
          MODULE_LIBRARIES  larreco_RecoAlg
			    larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			    lardata_RecoObjects larreco_Deprecated
			    lardataobj_RecoBase
			    larcorealg_Geometry
			    larcore_Geometry_Geometry_service
			    larreco_Calorimetry
			    lardataobj_AnalysisBase
			    lardata_Utilities
			    nutools_MagneticField_MagneticField_service
			    nusimdata_SimulationBase
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    art_Persistency_Common
			    art_Utilities
			    canvas
			    ${MF_MESSAGELOGGER}
			    
			    ${FHICLCPP}
			    cetlib cetlib_except
			    ${CLHEP}
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
         )

simple_plugin( Track3DKalman "module"
  			 larreco_Genfit 
			 larreco_RecoAlg
			 larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			 lardata_RecoObjects larreco_Deprecated
			 lardataobj_RecoBase
			 larcorealg_Geometry
			 larcore_Geometry_Geometry_service
			 larreco_Calorimetry
			 lardataobj_AnalysisBase
			 lardata_Utilities
			 nurandom_RandomUtils_NuRandomService_service
			 nutools_MagneticField_MagneticField_service
			 nusimdata_SimulationBase
			 ${ART_FRAMEWORK_CORE}
			 ${ART_FRAMEWORK_PRINCIPAL}
			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
			 ${ART_FRAMEWORK_SERVICES_BASIC}
			 ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			 ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			 art_Persistency_Common
			 art_Utilities
			 canvas
			 ${MF_MESSAGELOGGER}
			 
			 ${FHICLCPP}
			 cetlib cetlib_except
	   
        )

simple_plugin( Track3DKalmanSPS "module"
			 larreco_Genfit 
			 larreco_RecoAlg
			 larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			 lardata_RecoObjects larreco_Deprecated
			 lardataobj_RecoBase
			 larcorealg_Geometry
			 larcore_Geometry_Geometry_service
			 larreco_Calorimetry
			 lardataobj_AnalysisBase
			 lardata_Utilities
			 nutools_MagneticField_MagneticField_service
			 nusimdata_SimulationBase
			 ${ART_FRAMEWORK_CORE}
			 ${ART_FRAMEWORK_PRINCIPAL}
			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
			 ${ART_FRAMEWORK_SERVICES_BASIC}
			 ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			 ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			 art_Persistency_Common
			 art_Utilities
			 canvas
			 ${MF_MESSAGELOGGER}
			 
			 ${FHICLCPP}
			 cetlib cetlib_except
 	   
       )

simple_plugin( KalmanFilterFitTrackMaker "tool"
			    ${Boost_SYSTEM_LIBRARY}
			    ${Boost_FILESYSTEM_LIBRARY}
			    ${ROOT_CORE}
			    canvas
			    lardataobj_RecoBase
			    larreco_RecoAlg
       )

install_headers()
install_fhicl()
install_source()


