#ifndef RFFHITFINDERALG_H
#define RFFHITFINDERALG_H

/*!
 * Title:   RFFHitFinderAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Class that runs the RFF HitFinder. Implements an RFFHitFitter, and takes 
 * the result and stores it in recob::Hit objects. 
 *
 * Input:  recob::Wire
 * Output: recob::Hit
*/

#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "lardata/RecoBase/Wire.h"
#include "lardata/RecoBase/Hit.h"

#include "RFFHitFitter.h"

namespace hit{

  class RFFHitFinderAlg{

  const float SQRT_TWO_PI = 2.506628;

  public:
    RFFHitFinderAlg(fhicl::ParameterSet const&);

    void SetFitterParamsVectors(geo::Geometry const&);
    void Run(std::vector<recob::Wire> const&,
	     std::vector<recob::Hit>&,
	     geo::Geometry const&);

  private:

    std::vector<float> fMatchThresholdVec;
    std::vector<unsigned int> fMergeMultiplicityVec;
    std::vector<float> fAmpThresholdVec;

    void SetFitterParams(unsigned int);

    void EmplaceHit(std::vector<recob::Hit>&,
		    recob::Wire const&,
		    float const&,
		    raw::TDCtick_t const&, raw::TDCtick_t const&,
		    geo::SigType_t const&, geo::WireID const&);

    
    RFFHitFitter fFitter;
    
  };

}


#endif
