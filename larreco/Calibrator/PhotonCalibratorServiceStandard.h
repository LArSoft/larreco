////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorServiceStandard.h
//
// \brief Framework interface to PhotonCalibratorStandard
//
// \author ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORSERVICESTANDARD
#define PHOTONCALIBRATORSERVICESTANDARD

// LArSoft Includes
#include "larreco/Calibrator/PhotonCalibratorStandard.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"


// ART Includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "art/Framework/Principal/Run.h"

namespace calib {


  class PhotonCalibratorServiceStandard : public IPhotonCalibratorService
  {
  public:
    using provider_type = PhotonCalibratorStandard;

    struct ServiceConfiguration_t
    {
      fhicl::Atom<float> SPESize  { fhicl::Name("SPESize")  };
      fhicl::Atom<float> SPEShift { fhicl::Name("SPEShift") };
      fhicl::Atom<float> UseArea  { fhicl::Name("UseArea")  };
    };

    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

  public:
    PhotonCalibratorServiceStandard(Parameters const & config,
                                    art::ActivityRegistry& aReg)
      : fProvider( new PhotonCalibratorStandard(config.get_PSet(), aReg) )
      { }

    provider_type const* provider() const override { return fProvider.get(); }

  private:
    std::unique_ptr<PhotonCalibratorStandard> fProvider;
  };

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceStandard,
                                   calib::IPhotonCalibratorService,
                                   LEGACY)


#endif // PHOTONCALIBRATORSERVICESTANDARD
