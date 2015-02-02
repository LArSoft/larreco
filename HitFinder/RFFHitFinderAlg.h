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

#include "Geometry/Geometry.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"

#include "RFFHitFitter.h"

namespace hit{

  class RFFHitFinderAlg{

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

    RFFHitFitter fFitter;

  };

}


#endif
