#include "HitFilterAlg.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace hit{

  HitFilterAlg::HitFilterAlg(fhicl::ParameterSet const & p) :  
    fMinPulseHeight(p.get< std::vector<float> >("MinPulseHeight")),
    fMinPulseSigma(p.get< std::vector<float> >("MinPulseSigma"))
{
}

  bool HitFilterAlg::IsGoodHit(const recob::Hit& hit) const {

    const float              hitPH    = hit.PeakAmplitude();
    const float              hitSigma = hit.RMS();

    const geo::WireID& wireID   = hit.WireID();
    const size_t             view     = wireID.Plane;

    if (view >= fMinPulseSigma.size() || view >= fMinPulseHeight.size()) {
      mf::LogError("HitFilterAlg") << "Filtering settings not configured for all views! Will not filter hits in unconfigured views!";
      return true;
    }

    if ( hitPH    > fMinPulseHeight[view] &&
	 hitSigma > fMinPulseSigma[view] ) {
      return true;
    }
    else return false;
  }
}//end namespace hit
