cet_make_library(SOURCE CalorimetryAlg.cxx TrackCalorimetryAlg.cxx
  LIBRARIES
  PUBLIC
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RecoBase
  art::Framework_Services_Registry
  fhiclcpp::types
  ROOT::Physics
  PRIVATE
  larevt::ElectronLifetimeProvider
  larevt::ElectronLifetimeService
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataalg::DetectorInfo
  lardataobj::AnalysisBase
  lardataobj::TrackingTypes
  lardataobj::RecoBase
)

cet_build_plugin(CaloChecker art::EDAnalyzer
  LIBRARIES PRIVATE
  larevt::SpaceCharge
  larevt::SpaceChargeService
  lardata::ArtDataHelper
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcorealg::Geometry
  larcoreobj::headers
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  canvas::canvas
  cetlib::cetlib
  ROOT::Hist
  ROOT::MathCore
)

cet_build_plugin(Calorimetry art::EDProducer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  larevt::SpaceCharge
  larevt::SpaceChargeService
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardataalg::UtilitiesHeaders
  larcorealg::Geometry
  larcoreobj::headers
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::MathCore
  ROOT::Physics
)

cet_build_plugin(GnocchiCalorimetry art::EDProducer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larevt::ChannelStatusService
  larevt::SpaceCharge
  larevt::SpaceChargeService
  lardata::ArtDataHelper
  lardata::AssociationUtil
  lardata::DetectorPropertiesService
  larcore::ServiceUtil
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  ROOT::Hist
  ROOT::MathCore
)

cet_build_plugin(PrintCalorimetry art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(ShowerCalorimetry art::EDProducer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larevt::SpaceCharge
  larevt::SpaceChargeService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardata::AssociationUtil
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardataalg::UtilitiesHeaders
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  canvas::canvas
  ROOT::Physics
)

cet_build_plugin(TrackCalorimetry art::EDProducer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::LArPropertiesService
  lardata::ServicePack
  lardataalg::UtilitiesHeaders
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
