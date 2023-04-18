cet_make_library(LIBRARY_NAME PhotonCalibrator INTERFACE
  SOURCE IPhotonCalibrator.h
  LIBRARIES INTERFACE
  lardataobj::RecoBase
)

cet_make_library(LIBRARY_NAME PhotonCalibratorService INTERFACE
  SOURCE IPhotonCalibratorService.h
  LIBRARIES INTERFACE
  larreco::PhotonCalibrator
  art_plugin_types::serviceDeclaration
)

cet_make_library(LIBRARY_NAME PhotonCalibratorStandard INTERFACE
  SOURCE PhotonCalibratorStandard.h
  LIBRARIES INTERFACE
  larreco::PhotonCalibrator
  fhiclcpp::fhiclcpp
)

cet_write_plugin_builder(lar::PhotonCalibratorService
  art::service
  Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larreco::PhotonCalibratorService
)

include(lar::PhotonCalibratorService)

cet_build_plugin(PhotonCalibratorServiceStandard
  lar::PhotonCalibratorService
  IMPL_SOURCE # Empty: export a linkable INTERFACE target.
  LIBRARIES INTERFACE
  larreco::PhotonCalibratorStandard
  art::Framework_Principal
  fhiclcpp::types
  fhiclcpp::fhiclcpp
)

install_headers()
install_fhicl()
install_source()
