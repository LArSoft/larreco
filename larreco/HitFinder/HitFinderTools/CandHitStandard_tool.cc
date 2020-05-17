////////////////////////////////////////////////////////////////////////
/// \file   CandHitStandard.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

#include "art/Utilities/ToolMacros.h"
#include "larcore/Geometry/Geometry.h"

#include <algorithm>

namespace reco_tool
{

class CandHitStandard : ICandidateHitFinder
{
public:
    explicit CandHitStandard(const fhicl::ParameterSet& pset);

    void findHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t&,
                           const size_t,
                           const size_t,
                           const size_t,
                           HitCandidateVec&) const override;

    void MergeHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t&,
                            const HitCandidateVec&,
                            MergeHitCandidateVec&) const override;

private:

    void findHitCandidates(std::vector<float>::const_iterator,
                           std::vector<float>::const_iterator,
                           const size_t,
                           const size_t,
                           HitCandidateVec&) const;

    // Member variables from the fhicl file
    const float                    fRoiThreshold;               ///< minimum maximum to minimum peak distance

    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};

//----------------------------------------------------------------------
// Constructor.
CandHitStandard::CandHitStandard(const fhicl::ParameterSet& pset):
    fRoiThreshold(pset.get< float >("RoiThreshold", 5.))
    {
}

void CandHitStandard::findHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t& dataRange,
                                        const size_t                                         roiStartTick,
                                        const size_t                                         channel,
                                        const size_t                                         eventCount,
                                        HitCandidateVec&                                     hitCandidateVec) const
{
    // Recover the actual waveform
    const Waveform& waveform = dataRange.data();

    // Recover the plane index for this method
    std::vector<geo::WireID> wids  = fGeometry->ChannelToWire(channel);
    const size_t             plane = wids[0].Plane;

    // Use the recursive version to find the candidate hits
    findHitCandidates(waveform.begin(),waveform.end(),roiStartTick,plane,hitCandidateVec);

    return;
}

void CandHitStandard::findHitCandidates(std::vector<float>::const_iterator startItr,
                                        std::vector<float>::const_iterator stopItr,
                                        const size_t                       roiStartTick,
                                        const size_t                       planeIdx,
                                        HitCandidateVec&                   hitCandidateVec) const
{
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        auto maxItr = std::max_element(startItr, stopItr);

        float maxValue = *maxItr;
        int   maxTime  = std::distance(startItr,maxItr);

        if (maxValue > fRoiThreshold)
        {
            // backwards to find first bin for this candidate hit
            auto firstItr = std::distance(startItr,maxItr) > 2 ? maxItr - 1 : startItr;

            while(firstItr != startItr)
            {
                // Check for pathology where waveform goes too negative
                if (*firstItr < -fRoiThreshold) break;

                // Check both sides of firstItr and look for min/inflection point
                if (*firstItr < *(firstItr+1) && *firstItr <= *(firstItr-1)) break;

                firstItr--;
            }

            int firstTime = std::distance(startItr,firstItr);

            // Recursive call to find all candidate hits earlier than this peak
            findHitCandidates(startItr, firstItr + 1, roiStartTick, planeIdx, hitCandidateVec);

            // forwards to find last bin for this candidate hit
            auto lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;

            while(lastItr != stopItr - 1)
            {
                // Check for pathology where waveform goes too negative
                if (*lastItr < -fRoiThreshold) break;

                // Check both sides of firstItr and look for min/inflection point
                if (*lastItr <= *(lastItr+1) && *lastItr < *(lastItr-1)) break;

                lastItr++;
            }

            int lastTime = std::distance(startItr,lastItr);

            // Now save this candidate's start and max time info
            HitCandidate hitCandidate;
            hitCandidate.startTick     = roiStartTick + firstTime;
            hitCandidate.stopTick      = roiStartTick + lastTime;
            hitCandidate.maxTick       = roiStartTick + firstTime;
            hitCandidate.minTick       = roiStartTick + lastTime;
            hitCandidate.maxDerivative = *(startItr + firstTime);
            hitCandidate.minDerivative = *(startItr + lastTime);
            hitCandidate.hitCenter     = roiStartTick + maxTime;
            hitCandidate.hitSigma      = std::max(2.,float(lastTime - firstTime) / 6.);
            hitCandidate.hitHeight     = maxValue;

            hitCandidateVec.push_back(hitCandidate);

            // Recursive call to find all candidate hits later than this peak
            findHitCandidates(lastItr + 1, stopItr, roiStartTick + std::distance(startItr,lastItr + 1), planeIdx, hitCandidateVec);
        }
    }

    return;
}

void CandHitStandard::MergeHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t& rangeData,
                                         const HitCandidateVec&                               hitCandidateVec,
                                         MergeHitCandidateVec&                                mergedHitsVec) const
{
    // If no hits then nothing to do here
    if (hitCandidateVec.empty()) return;

    // The idea is to group hits that "touch" so they can be part of common fit, those that
    // don't "touch" are fit independently. So here we build the output vector to achieve that
    HitCandidateVec groupedHitVec;
    int             lastTick = hitCandidateVec.front().stopTick;

    // Step through the input hit candidates and group them by proximity
    for(const auto& hitCandidate : hitCandidateVec)
    {
        // Check condition that we have a new grouping
        if (int(hitCandidate.startTick) - lastTick > 1)
        {
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

DEFINE_ART_CLASS_TOOL(CandHitStandard)
}
