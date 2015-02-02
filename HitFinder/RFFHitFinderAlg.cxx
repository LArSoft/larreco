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

#include "RFFHitFinderAlg.h"

hit::RFFHitFinderAlg::RFFHitFinderAlg(fhicl::ParameterSet const& p)
{
  fMatchThresholdVec = p.get< std::vector<float> >("MeanMatchThreshold");
  fMergeMultiplicityVec = p.get< std::vector<unsigned int> >("MinMergeMultiplicity");
  fAmpThresholdVec = p.get< std::vector<float> >("AmplitudeThreshold",std::vector<float>(1,0.0));
}

void hit::RFFHitFinderAlg::SetFitterParamsVectors(geo::Geometry const& geo)
{
  const unsigned int n_planes = geo.Nplanes();

  //If size zero, throw. If size one, assume same for all planes.
  //If size > 1 but < n_planes, throw. If size = n_plane, good.

  if(fMatchThresholdVec.size()==0 || 
     fMergeMultiplicityVec.size()==0 ||
     fAmpThresholdVec.size()==0)
    throw std::runtime_error("Error in RFFHitFinderAlg: Configured with zero planes.");

  if( (fMatchThresholdVec.size()>1 && fMatchThresholdVec.size()<n_planes) ||
      (fMergeMultiplicityVec.size()>1 && fMergeMultiplicityVec.size()<n_planes) ||
      (fAmpThresholdVec.size()>1 && fAmpThresholdVec.size()<n_planes) )
    throw std::runtime_error("Error in RFFHitFinderAlg: Configured with incorrect n_planes.");
    
  if(fMatchThresholdVec.size()==1)
    fMatchThresholdVec.resize(n_planes,fMatchThresholdVec[0]);

  if(fMergeMultiplicityVec.size()==1)
    fMergeMultiplicityVec.resize(n_planes,fMergeMultiplicityVec[0]);

  if(fAmpThresholdVec.size()==1)
    fAmpThresholdVec.resize(n_planes,fAmpThresholdVec[0]);
}

void hit::RFFHitFinderAlg::SetFitterParams(unsigned int p)
{
  fFitter.SetFitterParams(fMatchThresholdVec[p],fMergeMultiplicityVec[p],fAmpThresholdVec[p]);
}

void hit::RFFHitFinderAlg::Run(std::vector<recob::Wire> const& wireVector,
			       std::vector<recob::Hit>& hitVector,
			       geo::Geometry const& geo)
{
  hitVector.reserve(wireVector.size());
  for(auto const& wire : wireVector){

    SetFitterParams(wire.View());
    for(auto const& roi : wire.SignalROI().get_ranges()){
      fFitter.RunFitter(roi.data());
    }

  }
}
