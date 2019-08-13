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

// LArSoft includes
#include "larreco/Calibrator/IPhotonCalibrator.h"

// ART includes
#include "art/Framework/Services/Registry/ServiceMacros.h"


namespace calib {
  class IPhotonCalibratorService {

  public:
    using provider_type = calib::IPhotonCalibrator;

  public:
    virtual ~IPhotonCalibratorService() = default;

    virtual provider_type const* provider() const = 0;

  }; // class IPhotonCalibratorService
}

DECLARE_ART_SERVICE_INTERFACE(calib::IPhotonCalibratorService, LEGACY)

#endif
