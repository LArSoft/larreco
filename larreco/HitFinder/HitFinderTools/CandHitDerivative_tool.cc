////////////////////////////////////////////////////////////////////////
/// \file   CandHitDerivative.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IWaveformAlgs.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <fstream>

namespace reco_tool
{

class CandHitDerivative : ICandidateHitFinder
{
public:
    explicit CandHitDerivative(const fhicl::ParameterSet& pset);
    
    ~CandHitDerivative();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void findHitCandidates(const std::vector<float>&,
                           size_t,
                           double,
                           HitCandidateVec&) const override;
    
    void findHitCandidates(std::vector<float>::const_iterator,
                           std::vector<float>::const_iterator,
                           size_t,
                           double,
                           HitCandidateVec&) const override;
    
    void MergeHitCandidates(const std::vector<float>&,
                            const HitCandidateVec&,
                            MergeHitCandidateVec&) const override;
    
private:
    // Member variables from the fhicl file
    std::unique_ptr<reco_tool::IWaveformAlgs> fWaveformAlgs;
    
    const geo::GeometryCore*  fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
CandHitDerivative::CandHitDerivative(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
CandHitDerivative::~CandHitDerivative()
{
}
    
void CandHitDerivative::configure(const fhicl::ParameterSet& pset)
{
    // Recover the baseline tool
    fWaveformAlgs = art::make_tool<reco_tool::IWaveformAlgs> (pset.get<fhicl::ParameterSet>("WaveformAlgs"));
    
    return;
}
    
void CandHitDerivative::findHitCandidates(const std::vector<float>& waveform,
                                          size_t                    roiStartTick,
                                          double                    roiThreshold,
                                          HitCandidateVec&          hitCandidateVec) const
{
    // In this case we want to find hit candidates based on the derivative of of the input waveform
    // We get this from our waveform algs too...
    std::vector<float> derivativeVec;
    
    fWaveformAlgs->getSmoothDerivativeVec(waveform, derivativeVec);
    
    // Now find the hits
    findHitCandidates(derivativeVec.begin(),derivativeVec.end(),roiStartTick,roiThreshold,hitCandidateVec);
    
    // Reset the hit height from the input waveform
    for(auto& hitCandidate : hitCandidateVec)
    {
        size_t centerIdx = hitCandidate.hitCenter;
        
        hitCandidate.hitHeight = waveform.at(centerIdx);
    }
    
    return;
}
    
void CandHitDerivative::findHitCandidates(std::vector<float>::const_iterator startItr,
                                          std::vector<float>::const_iterator stopItr,
                                          size_t                             roiStartTick,
                                          double                             roiThreshold,
                                          HitCandidateVec&                   hitCandidateVec) const
{
    // Require a minimum length
    size_t roiLength = std::distance(startItr,stopItr);
    
    if (roiLength > 2)
    {
        // Now search for candidate hits... our first is to simply find the range between the maximum and minimum bins in the
        // input derivative vector. We'll consider special cases further down.
        // So, start by looking for the peak bin, peak positive over baseline, then search from there for the largest negative bin
        std::vector<float>::const_iterator maxItr = std::max_element(startItr,stopItr);
        std::vector<float>::const_iterator minItr = std::min_element(maxItr,  stopItr);
        
        float range = *maxItr - *minItr;
        
        // At some point small rolling oscillations on the waveform need to be ignored...
        if (range > 0.1)
        {
            // Need to back up to find zero crossing
            std::vector<float>::const_iterator newEndItr = maxItr;
            
            while(newEndItr != startItr && *newEndItr > 0.) newEndItr--;
            
            size_t startTick = std::distance(startItr,newEndItr);
            
            // Find hits in the section of the waveform leading up to this candidate hit
            if (startTick > 2) findHitCandidates(startItr,newEndItr,roiStartTick,roiThreshold,hitCandidateVec);
            
            // Now need to go forward to again get close to zero
            std::vector<float>::const_iterator newStartItr = minItr;
            
            while(newStartItr != stopItr && *newStartItr < 0.) newStartItr++;
            
            size_t stopTick = std::distance(startItr,newStartItr);
            
            // The range from the max derivative to the min derivative found by a simple
            // search above may contain several intervening max/min. This can occur when
            // a delta ray is ejected and you have two hits overlapping. So we need to
            // step through from max to min and look for these special cases.
            std::vector<std::vector<float>::const_iterator> maxItrVec = {maxItr};
            std::vector<std::vector<float>::const_iterator> minItrVec;
            
            std::vector<float>::const_iterator adcItr  = maxItr;
            std::vector<float>::const_iterator lastItr = maxItr;
            bool                               foundMin(false);
            
            // We are gauranteed that we start at the maximum and end at the minimum. If
            // we find an intervening minimum then there MUST be a corresponding maximum
            // to get the derivative turning down again to reach the absolute minimum
            // This loop is constructed with that assumption in mind
            while(++adcItr != minItr)
            {
                // If we found an intervening minimum then switch modes to looking for next max
                if (foundMin)
                {
                    // Peak condition is that the waveform slope is decreasing again
                    if (*adcItr < *lastItr)
                    {
                        maxItrVec.push_back(lastItr);
                        foundMin = false;
                    }
                }
                // Condition for a minimum is that the waveform slope starts increasing
                else if (*adcItr > *lastItr)
                {
                    minItrVec.push_back(lastItr);
                    foundMin = true;
                }
                lastItr = adcItr;
            }
            
            minItrVec.push_back(minItr);
            
            // Fill candidate hit vector
            for(size_t idx = 0; idx < maxItrVec.size(); idx++)
            {
                // get a new hit object and fill it
                HitCandidate_t hitCandidate;
                
                std::vector<float>::const_iterator peakItr = std::min_element(maxItrVec.at(idx),minItrVec.at(idx),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
                
                hitCandidate.startTick     = roiStartTick + startTick;
                hitCandidate.stopTick      = roiStartTick + stopTick;
                hitCandidate.maxTick       = roiStartTick + std::distance(startItr,maxItrVec.at(idx));
                hitCandidate.minTick       = roiStartTick + std::distance(startItr,minItrVec.at(idx));
                hitCandidate.maxDerivative = maxItrVec.at(idx) != stopItr ? *maxItrVec.at(idx) : 0.;
                hitCandidate.minDerivative = minItrVec.at(idx) != stopItr ? *minItrVec.at(idx) : 0.;
                hitCandidate.hitCenter     = roiStartTick + std::distance(startItr,peakItr); //0.5 * float(hitCandidate.maxTick + hitCandidate.minTick);
                hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
                hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130; //0.6065;
                
                hitCandidateVec.push_back(hitCandidate);
            }
            
            // Finally, search the section of the waveform following this candidate for more hits
            findHitCandidates(newStartItr,stopItr,roiStartTick + stopTick,roiThreshold,hitCandidateVec);
        }
    }
    
    return;
}
    
void CandHitDerivative::MergeHitCandidates(const std::vector<float>& signalVec,
                                           const HitCandidateVec&    hitCandidateVec,
                                           MergeHitCandidateVec&     mergedHitsVec) const
{
    // If nothing on the input end then nothing to do
    if (hitCandidateVec.empty()) return;
    
    // The idea is to group hits that "touch" so they can be part of common fit, those that
    // don't "touch" are fit independently. So here we build the output vector to achieve that
    HitCandidateVec groupedHitVec;
    int             lastTick = hitCandidateVec.front().stopTick;
    
    // Step through the input hit candidates and group them by proximity
    for(const auto& hitCandidate : hitCandidateVec)
    {
        // Weed out the little people
        if (hitCandidate.hitHeight < 2.0) continue;
        
        // Check condition that we have a new grouping
        if (int(hitCandidate.startTick) - lastTick > 1 && !groupedHitVec.empty())
        {
            // Weed out the little people
            if (groupedHitVec.size() > 1 || groupedHitVec.front().hitHeight > 2.0)
                mergedHitsVec.emplace_back(groupedHitVec);
            
            groupedHitVec.clear();
        }

        // Add the current hit to the current group
        groupedHitVec.emplace_back(hitCandidate);
        
        lastTick = hitCandidate.stopTick;
    }
    
    // Check end condition
    if (!groupedHitVec.empty()) mergedHitsVec.emplace_back(groupedHitVec);
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(CandHitDerivative)
}
