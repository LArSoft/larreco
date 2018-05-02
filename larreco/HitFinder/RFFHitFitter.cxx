/*!
 * Title:   RFFHitFitter Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Class that does the base RFF algorithm. RFF works by simplifiying a Gaussian
 * fit by dividing a pulse by its derivative. for a Gaussian, the result is a
 * line, with the slope and intercept related to the sigma and mean of the 
 * Gaussian. 
 *
 * Input:  Signal (vector of floats)
 * Output: Guassian means and sigmas
*/

#include "RFFHitFitter.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "cetlib_except/exception.h"

hit::RFFHitFitter::RFFHitFitter(float step, float max):
  fGEAlg(step,max)
{}

hit::RFFHitFitter::RFFHitFitter(float        max_mean,
                                unsigned int min_multi,
                                float        threshold,
                                float        step,
                                float        max):
                                fGEAlg(step,max)
{
    SetFitterParams(max_mean,min_multi,threshold);
}

void hit::RFFHitFitter::SetFitterParams(float        max_mean,
                                        unsigned int min_multi,
                                        float        threshold)
{
    fMeanMatchThreshold = max_mean;
    fMinMergeMultiplicity = min_multi;
  
    if(fMinMergeMultiplicity==0) fMinMergeMultiplicity=1;

    fFinalAmpThreshold = threshold;

    ClearResults();
}				

void hit::RFFHitFitter::RunFitter(const std::vector<float>& signal)
{
    ClearResults();
    CalculateAllMeansAndSigmas(signal);
    CreateMergeVector();
    CalculateMergedMeansAndSigmas(signal.size());
    CalculateAmplitudes(signal);
}

void hit::RFFHitFitter::CalculateAllMeansAndSigmas(const std::vector<float>& signal)
{
    if(signal.size()<=2) return;
  
    float prev_dev=0,this_dev=0;;
    float slope=0; float sigma=0;
    float intercept=0; float mean=0;
    
    for(size_t i_tick=1; i_tick < signal.size()-1; i_tick++)
    {
        this_dev = 0.5*(signal[i_tick+1]-signal[i_tick-1])/signal[i_tick];
        slope = this_dev - prev_dev;
    
        prev_dev = this_dev;
    
        if(slope>=0) continue;
    
        sigma = std::sqrt(-1/slope);
        intercept = 0.5*(signal[i_tick+1]-signal[i_tick-1])/signal[i_tick] - slope*i_tick;
        mean = -1*intercept/slope;
    
        fSignalSet.insert(std::make_pair(mean,sigma));
    }
}

void hit::RFFHitFitter::CreateMergeVector()
{
    fMergeVector.clear(); fMergeVector.reserve( fSignalSet.size() );
  
    float prev_mean=-9e6;
    for(std::multiset<MeanSigmaPair>::iterator it=fSignalSet.begin(); it!=fSignalSet.end(); it++)
    {
        if( std::abs(it->first - prev_mean) > fMeanMatchThreshold || fMergeVector.size()==0 )
            fMergeVector.push_back( std::vector< std::multiset<MeanSigmaPair>::iterator >(1,it) );
        else
            fMergeVector.back().push_back(it);
        prev_mean = it->first;
    }
}

void hit::RFFHitFitter::CalculateMergedMeansAndSigmas(size_t signal_size)
{
    fMeanVector.reserve(fMergeVector.size());
    fSigmaVector.reserve(fMergeVector.size());
    fMeanErrorVector.reserve(fMergeVector.size());
    fSigmaErrorVector.reserve(fMergeVector.size());

    for(size_t i_col=0; i_col<fMergeVector.size(); i_col++)
    {
        if(fMergeVector[i_col].size()<fMinMergeMultiplicity) continue;

        fMeanVector.push_back(0.0);
        fSigmaVector.push_back(0.0);

        for(auto const& sigpair : fMergeVector[i_col])
        {
            fMeanVector.back() += sigpair->first;
            fSigmaVector.back() += sigpair->second;
        }

        fMeanVector.back() /= fMergeVector[i_col].size();
        fSigmaVector.back() /= fMergeVector[i_col].size();

        if(fMeanVector.back() < 0 || fMeanVector.back()>signal_size-1)
        {
            fMeanVector.pop_back();
            fSigmaVector.pop_back();
            continue;
        }

        fMeanErrorVector.push_back(0.0);
        fSigmaErrorVector.push_back(0.0);

        for(auto const& sigpair : fMergeVector[i_col])
        {
            fMeanErrorVector.back() +=
                (sigpair->first-fMeanVector.back())*(sigpair->first-fMeanVector.back());
            fSigmaErrorVector.back() +=
                (sigpair->second-fSigmaVector.back())*(sigpair->second-fSigmaVector.back());
        }

        fMeanErrorVector.back() = std::sqrt(fMeanErrorVector.back()) / fMergeVector[i_col].size();
        fSigmaErrorVector.back() = std::sqrt(fSigmaErrorVector.back()) / fMergeVector[i_col].size();

    }

}

void hit::RFFHitFitter::CalculateAmplitudes(const std::vector<float>& signal)
{
    std::vector<float> heightVector(fMeanVector.size());
    size_t bin=0;
  
    for(size_t i=0; i<fMeanVector.size(); i++)
    {
        if     (fMeanVector[i]<0)                 bin=0;
        else if(fMeanVector[i]+1 > signal.size()) bin=signal.size()-2;
        else                                      bin = std::floor(fMeanVector[i]);

        if(bin >= signal.size()-1)
            throw cet::exception("RFFHitFitter") << "Error in CalculatAmplitudes! bin is out of range!\n"
					   << "\tFor element " << i << " bin is " << bin << "(" << fMeanVector[i] << ")"
					   << " but size is " << signal.size() << ".\n";
    
        heightVector[i] = signal[bin] - (fMeanVector[i]-(float)bin)*(signal[bin]-signal[bin+1]);
    }
  
    fAmpVector = fGEAlg.SolveEquations(fMeanVector,fSigmaVector,heightVector);

    while(HitsBelowThreshold())
    {
        for(size_t i=0; i<fAmpVector.size(); i++)
        {
            if(fAmpVector[i] < fFinalAmpThreshold)
            {
                fMeanVector.erase(fMeanVector.begin()+i);
                fMeanErrorVector.erase(fMeanErrorVector.begin()+i);
                fSigmaVector.erase(fSigmaVector.begin()+i);
                fSigmaErrorVector.erase(fSigmaErrorVector.begin()+i);
                fAmpVector.erase(fAmpVector.begin()+i);
                heightVector.erase(heightVector.begin()+i);
            }
        }
        fAmpVector = fGEAlg.SolveEquations(fMeanVector,fSigmaVector,heightVector);
    }

    fAmpErrorVector.resize(fAmpVector.size(),0.0);
}

bool hit::RFFHitFitter::HitsBelowThreshold()
{
    for(auto const& amp : fAmpVector)
        if(amp < fFinalAmpThreshold) return true;
    return false;
}

void hit::RFFHitFitter::ClearResults()
{
    fMeanVector.clear();
    fSigmaVector.clear();
    fMeanErrorVector.clear();
    fSigmaErrorVector.clear();
    fAmpVector.clear();
    fAmpErrorVector.clear();
    fSignalSet.clear();
    fMergeVector.clear();
}

void hit::RFFHitFitter::PrintResults()
{
    std::cout << "InitialSignalSet" << std::endl;
  
    for(auto const& sigpair : fSignalSet)
        std::cout << "\t" << sigpair.first << " / " << sigpair.second << std::endl;
  
    std::cout << "\nNHits = " << NHits() << std::endl;
    std::cout << "\tMean / Sigma / Amp" << std::endl;
    for(size_t i=0; i<NHits(); i++)
        std::cout << "\t"
            << fMeanVector[i] <<  " +- " << fMeanErrorVector[i] << " / "
            << fSigmaVector[i] << " +- " << fSigmaErrorVector[i] << " / "
            << fAmpVector[i] << " +- " << fAmpErrorVector[i]
            << std::endl;
}
