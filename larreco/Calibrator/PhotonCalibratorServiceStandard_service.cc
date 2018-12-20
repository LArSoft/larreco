/**
 * Required minimal implementation file for calibrator service
 * which only returns a provider.
 */
#include "larreco/Calibrator/PhotonCalibratorServiceStandard.h"

DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceStandard, 
                                  calib::IPhotonCalibratorService)
