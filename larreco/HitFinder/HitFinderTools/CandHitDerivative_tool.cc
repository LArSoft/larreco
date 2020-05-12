////////////////////////////////////////////////////////////////////////
/// \file   CandHitDerivative.cc
/// \author T. Usher
//note for MT: this implementation is not thread-safe 
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IWaveformTool.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "larcore/Geometry/Geometry.h"

#include "TProfile.h"

#include <cmath>

namespace reco_tool
{

class CandHitDerivative : ICandidateHitFinder
{
public:
    explicit CandHitDerivative(const fhicl::ParameterSet& pset);

    void findHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t&,
                           const size_t,
                           const size_t,
                           const size_t,
                           HitCandidateVec&) const override;

    void MergeHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t&,
                            const HitCandidateVec&,
                            MergeHitCandidateVec&) const override;

private:
    // Internal functions
    void findHitCandidates(Waveform::const_iterator,
                           Waveform::const_iterator,
                           const size_t,
                           int,
                           float,
                           HitCandidateVec&) const;

    // Finding the nearest maximum/minimum from current point
    Waveform::const_iterator findNearestMax(Waveform::const_iterator, Waveform::const_iterator) const;
    Waveform::const_iterator findNearestMin(Waveform::const_iterator, Waveform::const_iterator) const;

    // handle finding the "start" and "stop" of a candidate hit
    Waveform::const_iterator findStartTick(Waveform::const_iterator, Waveform::const_iterator) const;
    Waveform::const_iterator findStopTick(Waveform::const_iterator, Waveform::const_iterator)  const;

    // some control variables
    size_t               fPlane;               //< Identifies the plane this tool is meant to operate on
    int                  fMinDeltaTicks;       //< minimum ticks from max to min to consider
    int                  fMaxDeltaTicks;       //< maximum ticks from max to min to consider
    float                fMinDeltaPeaks;       //< minimum maximum to minimum peak distance
    float                fMinHitHeight;        //< Drop candidate hits with height less than this
    size_t               fNumInterveningTicks; //< Number ticks between candidate hits to merge
    bool                 fOutputHistograms;    //< If true will generate a very large file of hists!

    art::TFileDirectory* fHistDirectory;

    // Global histograms
    TH1F*                fDStopStartHist;
    TH1F*                fDMaxTickMinTickHist;
    TH1F*                fDMaxDerivMinDerivHist;

    mutable std::map<size_t,int> fChannelCntMap;

    // Member variables from the fhicl file
    std::unique_ptr<reco_tool::IWaveformTool> fWaveformTool;

    const geo::GeometryCore*  fGeometry = lar::providerFrom<geo::Geometry>();
};

//----------------------------------------------------------------------
// Constructor.
CandHitDerivative::CandHitDerivative(const fhicl::ParameterSet& pset)
{
    fPlane               = pset.get< size_t >("Plane",               0);
    fMinDeltaTicks       = pset.get< int    >("MinDeltaTicks",       0);
    fMaxDeltaTicks       = pset.get< int    >("MaxDeltaTicks",       30);
    fMinDeltaPeaks       = pset.get< float  >("MinDeltaPeaks",       0.025);
    fMinHitHeight        = pset.get< float  >("MinHitHeight",        2.0);
    fNumInterveningTicks = pset.get< size_t >("NumInterveningTicks", 6);
    fOutputHistograms    = pset.get< bool   >("OutputHistograms",    false);

    // Recover the baseline tool
    fWaveformTool = art::make_tool<reco_tool::IWaveformTool> (pset.get<fhicl::ParameterSet>("WaveformAlgs"));

    // If asked, define the global histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;

        fHistDirectory = tfs.get();

        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("HitPlane_%1zu",fPlane));

        fDStopStartHist        = dir.make<TH1F>(Form("DStopStart_%1zu", fPlane), ";Delta Stop/Start;",    200, 0., 200.);
        fDMaxTickMinTickHist   = dir.make<TH1F>(Form("DMaxTMinT_%1zu",  fPlane), ";Delta Max/Min Tick;",  200, 0., 200.);
        fDMaxDerivMinDerivHist = dir.make<TH1F>(Form("DMaxDMinD_%1zu",  fPlane), ";Delta Max/Min Deriv;", 200, 0., 200.);
    }

    return;
}

void CandHitDerivative::findHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t& dataRange,
                                          const size_t                                         roiStartTick,
                                          const size_t                                         channel,
                                          const size_t                                         eventCount,
                                          HitCandidateVec&                                     hitCandidateVec) const
{
    // In this case we want to find hit candidates based on the derivative of of the input waveform
    // We get this from our waveform algs too...
    Waveform rawDerivativeVec;
    Waveform derivativeVec;

    // Recover the actual waveform
    const Waveform& waveform = dataRange.data();

    fWaveformTool->firstDerivative(waveform, rawDerivativeVec);
    fWaveformTool->triangleSmooth(rawDerivativeVec, derivativeVec);

    std::vector<geo::WireID> wids  = fGeometry->ChannelToWire(channel);
    size_t                   plane = wids[0].Plane;
    size_t                   cryo  = wids[0].Cryostat;
    size_t                   tpc   = wids[0].TPC;
    size_t                   wire  = wids[0].Wire;

    // Just make sure the input candidate hit vector has been cleared
    hitCandidateVec.clear();

    // Now find the hits
    findHitCandidates(derivativeVec.begin(),derivativeVec.end(),roiStartTick,fMinDeltaTicks,fMinDeltaPeaks,hitCandidateVec);

    if (hitCandidateVec.empty())
    {
        if (plane == 0)
        {
            std::cout << "** C/T/P: " << cryo << "/" << tpc << "/" << plane << ", wire: " << wire << " has not hits with input size: " << waveform.size() << std::endl;
        }
    }

    // Reset the hit height from the input waveform
    for(auto& hitCandidate : hitCandidateVec)
    {
        size_t centerIdx = hitCandidate.hitCenter;

        hitCandidate.hitHeight = waveform.at(centerIdx);
    }

    // Keep track of histograms if requested
    if (fOutputHistograms)
    {
        // Recover the details...
//        std::vector<geo::WireID> wids  = fGeometry->ChannelToWire(channel);
//        size_t                   plane = wids[0].Plane;
//        size_t                   cryo  = wids[0].Cryostat;
//        size_t                   tpc   = wids[0].TPC;
//        size_t                   wire  = wids[0].Wire;

        size_t channelCnt = fChannelCntMap[channel]++;

        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("HitPlane_%1zu/ev%04zu/c%1zut%1zuwire_%05zu",plane,eventCount,cryo,tpc,wire));

        size_t waveformSize = waveform.size();
        int    waveStart    = roiStartTick;
        int    waveStop     = waveStart + waveformSize;

        TProfile* waveHist     = dir.make<TProfile>(Form("HWfm_%03zu_ctw%01zu-%01zu-%01zu-%05zu",channelCnt,cryo,tpc,plane,wire), "Waveform",   waveformSize, waveStart, waveStop, -500., 500.);
        TProfile* derivHist    = dir.make<TProfile>(Form("HDer_%03zu_ctw%01zu-%01zu-%01zu-%05zu",channelCnt,cryo,tpc,plane,wire), "Derivative", waveformSize, waveStart, waveStop, -500., 500.);
        TProfile* candHitHist  = dir.make<TProfile>(Form("HCan_%03zu_ctw%01zu-%01zu-%01zu-%05zu",channelCnt,cryo,tpc,plane,wire), "Cand Hits",  waveformSize, waveStart, waveStop, -500., 500.);
        TProfile* maxDerivHist = dir.make<TProfile>(Form("HMax_%03zu_ctw%01zu-%01zu-%01zu-%05zu",channelCnt,cryo,tpc,plane,wire), "Maxima",     waveformSize, waveStart, waveStop, -500., 500.);

        // Fill wave/derivative
        for(size_t idx = 0; idx < waveform.size(); idx++)
        {
            waveHist->Fill(roiStartTick + idx, waveform.at(idx));
            derivHist->Fill(roiStartTick + idx, derivativeVec.at(idx));
        }

        // Fill hits
        for(const auto& hitCandidate : hitCandidateVec)
        {
            candHitHist->Fill(hitCandidate.hitCenter, hitCandidate.hitHeight);
            maxDerivHist->Fill(hitCandidate.maxTick, hitCandidate.maxDerivative);
            maxDerivHist->Fill(hitCandidate.minTick, hitCandidate.minDerivative);

            fDStopStartHist->Fill(hitCandidate.stopTick - hitCandidate.startTick, 1.);
            fDMaxTickMinTickHist->Fill(hitCandidate.minTick - hitCandidate.maxTick, 1.);
            fDMaxDerivMinDerivHist->Fill(hitCandidate.maxDerivative - hitCandidate.minDerivative, 1.);
        }
    }

    return;
}

void CandHitDerivative::findHitCandidates(Waveform::const_iterator startItr,
                                          Waveform::const_iterator stopItr,
                                          const size_t             roiStartTick,
                                          int                      dTicksThreshold,
                                          float                    dPeakThreshold,
                                          HitCandidateVec&         hitCandidateVec) const
{
    // Search for candidate hits...
    // The idea will be to find the largest deviation in the input derivative waveform as the starting point. Depending
    // on if a maximum or minimum, we search forward or backward to find the minimum or maximum that our extremum
    // corresponds to.
    std::pair<Waveform::const_iterator, Waveform::const_iterator> minMaxPair = std::minmax_element(startItr, stopItr);

    Waveform::const_iterator maxItr = minMaxPair.second;
    Waveform::const_iterator minItr = minMaxPair.first;

    // Use the larger of the two as the starting point and recover the nearest max or min
    if (std::fabs(*maxItr) > std::fabs(*minItr)) minItr = findNearestMin(maxItr, stopItr);
    else                                         maxItr = findNearestMax(minItr,startItr);

    int   deltaTicks = std::distance(maxItr,minItr);
    float range      = *maxItr - *minItr;

    // At some point small rolling oscillations on the waveform need to be ignored...
    if (deltaTicks >= dTicksThreshold && range > dPeakThreshold)
    {
        // Need to back up to find zero crossing, this will be the starting point of our
        // candidate hit but also the endpoint of the pre sub-waveform we'll search next
        Waveform::const_iterator newEndItr = findStartTick(maxItr, startItr);

        int startTick = std::distance(startItr,newEndItr);

        // Now need to go forward to again get close to zero, this will then be the end point
        // of our candidate hit and the starting point for the post sub-waveform to search
        Waveform::const_iterator newStartItr = findStopTick(minItr, stopItr);

        int stopTick = std::distance(startItr,newStartItr);

        // Find hits in the section of the waveform leading up to this candidate hit
        if (startTick > dTicksThreshold)
        {
            // Special handling for merged hits
            if (*(newEndItr-1) > 0.) {dTicksThreshold = 2;              dPeakThreshold = 0.;            }
            else                     {dTicksThreshold = fMinDeltaTicks; dPeakThreshold = fMinDeltaPeaks;}

            findHitCandidates(startItr,newEndItr+1,roiStartTick,dTicksThreshold,dPeakThreshold,hitCandidateVec);
        }

        // Create a new hit candidate and store away
        HitCandidate hitCandidate;

        Waveform::const_iterator peakItr = std::min_element(maxItr,minItr,[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});

        // Check balance
        if      (2 * std::distance(peakItr,minItr) < std::distance(maxItr,peakItr)) peakItr--;
        else if (2 * std::distance(maxItr,peakItr) < std::distance(peakItr,minItr)) peakItr++;

        hitCandidate.startTick     = roiStartTick + startTick;
        hitCandidate.stopTick      = roiStartTick + stopTick;
        hitCandidate.maxTick       = roiStartTick + std::distance(startItr,maxItr);
        hitCandidate.minTick       = roiStartTick + std::distance(startItr,minItr);
        hitCandidate.maxDerivative = maxItr != stopItr ? *maxItr : 0.;
        hitCandidate.minDerivative = minItr != stopItr ? *minItr : 0.;
        hitCandidate.hitCenter     = roiStartTick + std::distance(startItr,peakItr) + 0.5;
        hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
        hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130;

        hitCandidateVec.push_back(hitCandidate);

        // Finally, search the section of the waveform following this candidate for more hits
        if (std::distance(newStartItr,stopItr) > dTicksThreshold)
        {
            // Special handling for merged hits
            if (*(newStartItr+1) < 0.) {dTicksThreshold = 2;              dPeakThreshold = 0.;            }
            else                       {dTicksThreshold = fMinDeltaTicks; dPeakThreshold = fMinDeltaPeaks;}

            findHitCandidates(newStartItr,stopItr,roiStartTick + stopTick,dTicksThreshold,dPeakThreshold,hitCandidateVec);
        }
    }

    return;
}

void CandHitDerivative::MergeHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t& rangeData,
                                           const HitCandidateVec&                               hitCandidateVec,
                                           MergeHitCandidateVec&                                mergedHitsVec) const
{
    // If nothing on the input end then nothing to do
    if (hitCandidateVec.empty()) return;

    // The idea is to group hits that "touch" so they can be part of common fit, those that
    // don't "touch" are fit independently. So here we build the output vector to achieve that
    // Get a container for the hits...
    HitCandidateVec groupedHitVec;

    // Initialize the end of the last hit which we'll set to the first input hit's stop
    size_t lastStopTick = hitCandidateVec.front().stopTick;

    // Step through the input hit candidates and group them by proximity
    for(const auto& hitCandidate : hitCandidateVec)
    {
        // Small pulse height hits should not be considered?
        if (hitCandidate.hitHeight > fMinHitHeight)
        {
            // Check condition that we have a new grouping
            if (hitCandidate.startTick > lastStopTick + fNumInterveningTicks && !groupedHitVec.empty())
            {
                mergedHitsVec.emplace_back(groupedHitVec);

                groupedHitVec.clear();
            }

            // Add the current hit to the current group
            groupedHitVec.emplace_back(hitCandidate);

            lastStopTick = hitCandidate.stopTick;
        }
    }

    // Check end condition
    if (!groupedHitVec.empty()) mergedHitsVec.emplace_back(groupedHitVec);

    return;
}

ICandidateHitFinder::Waveform::const_iterator CandHitDerivative::findNearestMin(Waveform::const_iterator maxItr, Waveform::const_iterator stopItr) const
{
    // reset the min iterator and search forward to find the nearest minimum
    Waveform::const_iterator lastItr = maxItr;

    // The strategy is simple...
    // We are at a maximum so we search forward until we find the lowest negative point
    while((lastItr + 1) != stopItr)
    {
        if (*(lastItr + 1) > *lastItr) break;

        lastItr++;
    }

    // The minimum will be the last iterator value...
    return lastItr;
}

ICandidateHitFinder::Waveform::const_iterator CandHitDerivative::findNearestMax(Waveform::const_iterator minItr, Waveform::const_iterator startItr) const
{
    // Set the internal loop variable...
    Waveform::const_iterator lastItr = minItr;

    // One extra condition to watch for here, make sure we can actually back up!
    if (std::distance(startItr,minItr) > 0)
    {
        // Similar to searching for a maximum, we loop backward over ticks looking for the waveform to start decreasing
        while((lastItr - 1) != startItr)
        {
            if (*(lastItr - 1) < *lastItr) break;

            lastItr--;
        }
    }

    return lastItr;
}

ICandidateHitFinder::Waveform::const_iterator CandHitDerivative::findStartTick(Waveform::const_iterator maxItr, Waveform::const_iterator startItr) const
{
    Waveform::const_iterator lastItr = maxItr;

    // If we can't back up then there is nothing to do
    if (std::distance(startItr,lastItr) > 0)
    {
        // In theory, we are starting at a maximum and want to find the "start" of the candidate peak
        // Ideally we would look to search backward to the point where the (derivative) waveform crosses zero again.
        // However, the complication is that we need to watch for the case where two peaks are merged together and
        // we might run through another peak before crossing zero...
        // So... loop until we hit the startItr...
        Waveform::const_iterator loopItr = lastItr - 1;

        while(loopItr != startItr)
        {
            // Ideal world case, we cross zero... but we might encounter a minimum... or an inflection point
            if (*loopItr < 0. || !(*loopItr < *lastItr)) break;

            lastItr = loopItr--;
        }
    }
    else lastItr = startItr;

    return lastItr;
}

ICandidateHitFinder::Waveform::const_iterator CandHitDerivative::findStopTick(Waveform::const_iterator minItr, Waveform::const_iterator stopItr)   const
{
    Waveform::const_iterator lastItr = minItr;

    // If we can't go forward then there is really nothing to do
    if (std::distance(minItr,stopItr) > 1)
    {
        // Pretty much the same strategy as for finding the start tick...
        Waveform::const_iterator loopItr = lastItr + 1;

        while(loopItr != stopItr)
        {
            // Ideal case that we have crossed zero coming from a minimum... but watch for a maximum as well
            if (*loopItr > 0. || !(*loopItr > *lastItr)) break;

            lastItr = loopItr++;
        }
    }

    return lastItr;
}

DEFINE_ART_CLASS_TOOL(CandHitDerivative)
}
