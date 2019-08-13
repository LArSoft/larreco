////////////////////////////////////////////////////////////////////////
// \file IPhotonCalibrator.h
//
// \brief Generic interface for service provider for calibrating optical hits
//
// \author ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef IPHOTONCALIBRATOR_H
#define IPHOTONCALIBRATOR_H

// LArSoft includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"




namespace calib
{
  /// May want to swap in dummy charge and photon calibrations in various
  /// combinations.
  class IPhotonCalibrator
  {
  public:
    // Make non-copiable
    IPhotonCalibrator(IPhotonCalibrator const &) = delete;
    IPhotonCalibrator(IPhotonCalibrator &&) = delete;
    IPhotonCalibrator& operator= (IPhotonCalibrator const &) = delete;
    IPhotonCalibrator& operator= (IPhotonCalibrator &&) = delete;

    IPhotonCalibrator() {};
    virtual ~IPhotonCalibrator() = default;

    virtual double PE(double adcs, int opchannel) const = 0;
    virtual bool   UseArea() const = 0;

    /// Need a 3D position because result depends on position along length of
    /// bar. This is going to be pretty imprecise even so.
    // virtual double GeV(double PE, int opchannel, TVector3 pos) = 0;

    /// Convenience
    double PE(const recob::OpHit& oh) const
    {
      return oh.PE();
    }

    double PE(const recob::OpFlash& of) const
    {
      return of.TotalPE();
    }

    // double GeV(const recob::OpHit& oh, TVector3 pos)
    //{
    //  return GeV(oh.PE(), oh.OpChannel(), pos);
    //}

    //double GeV(const OpFlash& of, TVector3 pos)
    //{
    //  // This function would be in the .cxx in practice
    //  const std::vector<double>& pes = of.PEs();
    //  double ret = 0;
    //  for(int chan = 0; chan < pes.size(); ++chan)
    //    ret += GeV(of.PE(chan), chan, pos);
    //  return ret;
    //}


  };

}


#endif
