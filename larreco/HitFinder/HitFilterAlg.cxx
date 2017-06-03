#include "HitFilterAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace hit{

  HitFilterAlg::HitFilterAlg(fhicl::ParameterSet const & p) {   
    this->reconfigure(p);
  }
  
  void HitFilterAlg::reconfigure(fhicl::ParameterSet const & p) {
    // Implementation of optional member function here.
    fMinPulseHeight   = p.get< std::vector<float> >("MinPulseHeight");
    fMinPulseSigma    = p.get< std::vector<float> >("MinPulseSigma");
    fMaxIntegralToADC = p.get< std::vector<float> >("MaxIntegralToADC");  
  }

  bool HitFilterAlg::IsGoodHit(const recob::Hit& hit) const {

    float hitPH    = hit.PeakAmplitude();
    float hitSigma = hit.RMS();
    
    float hitSummedADC = hit.SummedADC();
    float hitIntegral  = hit.Integral();
    float val = (hitSummedADC - hitIntegral)/(hitIntegral + hitSummedADC);

    const geo::WireID& wireID   = hit.WireID();
    size_t             view     = wireID.Plane;

    if (view >= fMinPulseSigma.size() || view >= fMinPulseHeight.size() || view >= fMaxIntegralToADC.size()) {
      mf::LogError("HitFilterAlg") << "Filtering settings not configured for all views! Will not filter hits in unconfigured views!";
      return true;
    }
    
    if ( hitPH    > fMinPulseHeight[view] &&
	 hitSigma > fMinPulseSigma[view] &&
	 val < fMaxIntegralToADC[view] ) {
      return true;
    }
    else return false;
  }
}//end namespace hit
