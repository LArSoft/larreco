#ifndef RFFHITFITTER_H
#define RFFHITFITTER_H

/*!
 * Title:   RFFHitFitter Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Class that does the base RFF algorithm. RFF works by simplifiying a Gaussian
 * fit by dividing a pulse by its derivative. for a Guassian, the result is a
 * line, with the slope and intercept related to the sigma and mean of the 
 * Gaussian. 
 *
 * Input:  Signal (vector of floats)
 * Output: Guassian means and sigmas
*/

#include <vector>
#include <set>

#include "GaussianEliminationAlg.h"

namespace hit{

  struct SignalSetComp{
    bool operator() (const std::pair<float,float>& lhs,
		     const std::pair<float,float>& rhs) const
    { return lhs.first < rhs.first; } 		     
  };

  class RFFHitFitter {

    typedef std::pair<float,float> MeanSigmaPair;
    
  public:
    RFFHitFitter(float,unsigned int,float,float step=0.1,float max=5.0);
    RFFHitFitter(float step=0.1,float max=5.0);

    void SetFitterParams(float,unsigned int,float);

    void RunFitter(const std::vector<float>& signal);

    const std::vector<float>& MeanVector() { return fMeanVector; }
    const std::vector<float>& SigmaVector() { return fSigmaVector; }
    const std::vector<float>& MeanErrorVector() { return fMeanErrorVector; }
    const std::vector<float>& SigmaErrorVector() { return fSigmaErrorVector; }
    const std::vector<float>& AmplitudeVector() { return fAmpVector; }
    const std::vector<float>& AmplitudeErrorVector() { return fAmpErrorVector; }
    unsigned int NHits() { return fMeanVector.size(); }

    void ClearResults();

    void PrintResults();
			 
  private:
    float fMeanMatchThreshold;
    unsigned int fMinMergeMultiplicity;
    float fFinalAmpThreshold;

    float fGEAlgStepSize;
    float fGEAlgMax;
    util::GaussianEliminationAlg fGEAlg;
    
    std::vector<float> fMeanVector;
    std::vector<float> fSigmaVector;
    std::vector<float> fMeanErrorVector;
    std::vector<float> fSigmaErrorVector;
    std::vector<float> fAmpVector;
    std::vector<float> fAmpErrorVector;

    std::multiset< MeanSigmaPair, SignalSetComp > fSignalSet;
    std::vector< std::vector< std::multiset<MeanSigmaPair>::iterator > >
      fMergeVector;

    void CalculateAllMeansAndSigmas(const std::vector<float>& signal);
    void CalculateMergedMeansAndSigmas();
    void CalculateAmplitudes(const std::vector<float>& signal);
    void CreateMergeVector();

    bool HitsBelowThreshold();

  };

}

#endif
