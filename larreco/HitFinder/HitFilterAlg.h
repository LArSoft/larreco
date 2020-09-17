////////////////////////////////////////////////////////////////////////
// Class:       HitFilterAlg
// Purpose:     Provide configurable hit filtering functions for identifying
//              noise and other undesireable hits
//
// Original code by Tracy Usher, converted to a larsoft algorithm by Brandon Eberly
////////////////////////////////////////////////////////////////////////

#ifndef HITFILTERALG_H
#define HITFILTERALG_H

#include <vector>

namespace fhicl { class ParameterSet; }
namespace recob { class Hit; }

namespace hit{

  class HitFilterAlg {
    public:
      explicit HitFilterAlg(fhicl::ParameterSet const & p);
      virtual ~HitFilterAlg() {}

      bool IsGoodHit(const recob::Hit& hit) const;

    private:

      const std::vector<float> fMinPulseHeight;
      const std::vector<float> fMinPulseSigma;
  };

}//end namespace hit
#endif
