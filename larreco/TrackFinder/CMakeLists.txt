
include_directories ( $ENV{G4INCLUDE} )

art_make(  
          EXCLUDE Track3DKalman_module.cc Track3DKalmanSPS_module.cc DumpTracks_module.cc
          MODULE_LIBRARIES  RecoAlg
			    Simulation
			    RecoObjects
			    RecoBase
			    Geometry
			    larcore_Geometry_Geometry_service
			    AnalysisAlg
			    AnalysisBase
			    Utilities
			    ${MAGNETICFIELD_SERVICE}
			    ${SIMULATIONBASE}
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    ${ART_PERSISTENCY_COMMON}
			    ${ART_UTILITIES}
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
			    ${FHICLCPP}
			    ${CETLIB}
			    ${CLHEP}
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
         )

simple_plugin( Track3DKalman "module"
  			 Genfit 
			 RecoAlg
			 Simulation
			 RecoObjects
			 RecoBase
			 Geometry
			 larcore_Geometry_Geometry_service
			 AnalysisAlg
			 AnalysisBase
			 Utilities
			 SeedService_service # artextensions
			 ${MAGNETICFIELD_SERVICE}
			 ${SIMULATIONBASE}
			 ${ART_FRAMEWORK_CORE}
			 ${ART_FRAMEWORK_PRINCIPAL}
			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
			 ${ART_FRAMEWORK_SERVICES_BASIC}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			 ${ART_PERSISTENCY_COMMON}
			 ${ART_UTILITIES}
			 ${MF_MESSAGELOGGER}
			 ${MF_UTILITIES}
			 ${FHICLCPP}
			 ${CETLIB}
	   
        )

simple_plugin( Track3DKalmanSPS "module"
			 Genfit 
			 RecoAlg
			 Simulation
			 RecoObjects
			 RecoBase
			 Geometry
			 larcore_Geometry_Geometry_service
			 AnalysisAlg
			 AnalysisBase
			 Utilities
			 ${MAGNETICFIELD_SERVICE}
			 ${SIMULATIONBASE}
			 ${ART_FRAMEWORK_CORE}
			 ${ART_FRAMEWORK_PRINCIPAL}
			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
			 ${ART_FRAMEWORK_SERVICES_BASIC}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			 ${ART_PERSISTENCY_COMMON}
			 ${ART_UTILITIES}
			 ${MF_MESSAGELOGGER}
			 ${MF_UTILITIES}
			 ${FHICLCPP}
			 ${CETLIB}
 	   
       )


simple_plugin( DumpTracks "module"
               RecoBase
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${MF_MESSAGELOGGER}
               ${MF_UTILITIES}
               ${FHICLCPP}
               ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
               
        )


install_headers()
install_fhicl()
install_source()


