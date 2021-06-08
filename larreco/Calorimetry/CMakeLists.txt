art_make(LIB_LIBRARIES
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           lardata_ArtDataHelper
           larcorealg_Geometry
           lardataobj_AnalysisBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_UTILITIES}
           ${CANVAS}
           ${CETLIB_EXCEPT}
           fhiclcpp::fhiclcpp fhiclcpp::types
           messagefacility::MF_MessageLogger
           ROOT::Core
           ROOT::Physics
         MODULE_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           messagefacility::MF_MessageLogger)

install_headers()
install_fhicl()
install_source()
