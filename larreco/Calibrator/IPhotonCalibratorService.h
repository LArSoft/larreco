////////////////////////////////////////////////////////////////////////
// \file IPhotonCalibratorSerice.h
//
// \brief Generic framework interface to IPhotonCalibrator
//
// \author ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef IPHOTONCALIBRATORSERVICE_H
#define IPHOTONCALIBRATORSERVICE_H

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"

namespace calib {
  class IPhotonCalibratorService {
  public:
    using provider_type = calib::IPhotonCalibrator;

    virtual ~IPhotonCalibratorService() = default;
    virtual provider_type const* provider() const = 0;
  };
}

DECLARE_ART_SERVICE_INTERFACE(calib::IPhotonCalibratorService, SHARED)

#endif
