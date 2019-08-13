////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorStandard.h
//
// \brief Standard service provider applying a flat scale factor to all optical hits.
//
// \author ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORDEFAULT_H
#define PHOTONCALIBRATORDEFAULT_H

#include "larreco/Calibrator/IPhotonCalibrator.h"

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"


namespace calib {

  class PhotonCalibratorStandard : public IPhotonCalibrator
  {
  public:
    PhotonCalibratorStandard(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
      : fSPESize  ( pset.get< float >("SPESize")     ),
        fSPEShift ( pset.get< float >("SPEShift", 0.)),
        fUseArea  ( pset.get< bool  >("UseArea")     )
      {}


    PhotonCalibratorStandard(float size, float shift, bool useArea)
      : fSPESize  ( size ),
        fSPEShift ( shift ),
        fUseArea  ( useArea )
      {}

    // Override base class functions
    double PE(double adcs, int opchannel) const override { return adcs/fSPESize + fSPEShift; }
    bool   UseArea() const override                      { return fUseArea; }

    // Setters for this implementation
    void SetSPESize (float size)    { fSPESize  = size; }
    void SetSPEShift(float shift)   { fSPEShift = shift; }
    void SetUseArea (bool  useArea) { fUseArea  = useArea; }

    /// Need a 3D position because result depends on position along length of
    /// bar. This is going to be pretty imprecise even so.
    // virtual double GeV(double PE, int opchannel, TVector3 pos) override;

  private:
    float  fSPESize;
    float  fSPEShift;
    bool   fUseArea;


  }; // class PhotonCalibratorStandard
}



#endif
