////////////////////////////////////////////////////////////////////////
/// \file   CandHitMorphological.cc
/// \author T. Usher
// MT note: This implementation is not thread-safe. 
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

class CandHitMorphological : ICandidateHitFinder
{
public:
    explicit CandHitMorphological(const fhicl::ParameterSet& pset);

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
    //< Top level hit finding using erosion/dilation vectors
    void findHitCandidates(Waveform::const_iterator, Waveform::const_iterator,   //< derivative
                           Waveform::const_iterator, Waveform::const_iterator,   //< erosion
                           Waveform::const_iterator, Waveform::const_iterator,   //< dilation
                           const size_t,
                           float,
                           HitCandidateVec&) const;

    //< Fine grain hit finding within candidate peak regions using derivative method
    void findHitCandidates(Waveform::const_iterator,
                           Waveform::const_iterator,
                           const size_t,
                           int,
                           float,
                           HitCandidateVec&) const;

    //< For a given range, return the list of max/min pairs
    using MaxMinPair       = std::pair<Waveform::const_iterator,Waveform::const_iterator>;
    using CandHitParams    = std::tuple<Waveform::const_iterator,Waveform::const_iterator,Waveform::const_iterator,Waveform::const_iterator>;
    using CandHitParamsVec = std::vector<CandHitParams>;

    bool getListOfHitCandidates(Waveform::const_iterator,
                                Waveform::const_iterator,
                                int,
                                float,
                                CandHitParamsVec&) const;

    // Finding the nearest maximum/minimum from current point
    Waveform::const_iterator findNearestMax(Waveform::const_iterator, Waveform::const_iterator) const;
    Waveform::const_iterator findNearestMin(Waveform::const_iterator, Waveform::const_iterator) const;

    // handle finding the "start" and "stop" of a candidate hit
    Waveform::const_iterator findStartTick(Waveform::const_iterator, Waveform::const_iterator) const;
    Waveform::const_iterator findStopTick(Waveform::const_iterator, Waveform::const_iterator)  const;

    // some fhicl control variables
    const size_t               fPlane;                //< Identifies the plane this tool is meant to operate on
    const float                fDilationThreshold;    //< Dilation threshold
    const float                fDilationFraction;     //< Fraction of max dilation to set for min dilation
    const float                fErosionFraction;      //< Fraction of max dilation value to set min erosion
    const int                  fMinDeltaTicks;        //< minimum ticks from max to min to consider
    const float                fMinDeltaPeaks;        //< minimum maximum to minimum peak distance
    const float                fMinHitHeight;         //< Drop candidate hits with height less than this
    const size_t               fNumInterveningTicks;  //< Number ticks between candidate hits to merge
    const int                  fStructuringElement;   //< Window size for morphologcial filter
    const bool                 fOutputHistograms;     //< If true will generate summary style histograms
    const bool                 fOutputWaveforms;      //< If true will output waveform related info <<< very big output file!
    const float                fFitNSigmaFromCenter;  //< Limit ticks to fit to NSigma from hit center; not applied if zero or negative

    art::TFileDirectory* fHistDirectory;

    // Global histograms
    TH1F*                fDStopStartHist;        //< Basically keeps track of the length of hit regions
    TH1F*                fDMaxTickMinTickHist;   //< This will be a measure of the width of candidate hits
    TH1F*                fDMaxDerivMinDerivHist; //< This is the difference peak to peak of derivative for cand hit
    TH1F*                fMaxErosionHist;        //< Keep track of the maximum erosion
    TH1F*                fMaxDilationHist;       //< Keep track of the maximum dilation
    TH1F*                fMaxDilEroRatHist;      //< Ratio of the maxima of the two

    //MT note: The mutable data members are only used in the histogram filling functions
    //and histogram filling can only be done in single-threaded mode.
    //Will need to consider design changes if this behavior changes. 
    mutable size_t       fLastChannel;           //< Kludge to keep track of last channel when histogramming in effect
    mutable size_t       fChannelCnt;            //< Counts the number of times a channel is used (assumed in order)

    //< All of the real work is done in the waveform tool
    std::unique_ptr<reco_tool::IWaveformTool> fWaveformTool;

    const geo::GeometryCore*  fGeometry = lar::providerFrom<geo::Geometry>();
};

//----------------------------------------------------------------------
// Constructor.
CandHitMorphological::CandHitMorphological(const fhicl::ParameterSet& pset):
    fPlane              (pset.get< size_t >("Plane",               0)),
    fDilationThreshold  (pset.get< float  >("DilationThreshold",   4.)),
    fDilationFraction   (pset.get< float  >("DilationFraction",    0.75)),
    fErosionFraction    (pset.get< float  >("ErosionFraction",     0.2)),
    fMinDeltaTicks      (pset.get< int    >("MinDeltaTicks",       0)),
    fMinDeltaPeaks      (pset.get< float  >("MinDeltaPeaks",       0.025)),
    fMinHitHeight       (pset.get< float  >("MinHitHeight",        1.0)),
    fNumInterveningTicks(pset.get< size_t >("NumInterveningTicks", 6)),
    fStructuringElement (pset.get< int    >("StructuringElement",  20)),
    fOutputHistograms   (pset.get< bool   >("OutputHistograms",    false)),
    fOutputWaveforms    (pset.get< bool   >("OutputWaveforms",     false)),
    fFitNSigmaFromCenter(pset.get< float  >("FitNSigmaFromCenter", 5.))
{

if (art::Globals::instance()->nthreads() > 1u) {
  if (fOutputHistograms) {
         throw art::Exception(art::errors::Configuration) << "Cannot fill histograms when multiple threads configured, please set fOutputHistograms to false or change number of threads to 1\n";
     }
  
  if (fOutputWaveforms) {
         throw art::Exception(art::errors::Configuration) << "Cannot write output waveforms when multiple threads configured, please set fOutputHistograms to false or change number of threads to 1\n";
  }
  }
  // Recover the baseline tool
    fWaveformTool = art::make_tool<reco_tool::IWaveformTool> (pset.get<fhicl::ParameterSet>("WaveformAlgs"));

    // Set the last channel to some nonsensical value
    fLastChannel = std::numeric_limits<size_t>::max();
    fChannelCnt  = 0;

    // If asked, define the global histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;

        fHistDirectory = tfs.get();

        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("HitPlane_%1zu",fPlane));

        fDStopStartHist        = dir.make<TH1F>(Form("DStopStart_%1zu",   fPlane), ";Delta Stop/Start;",    100,   0., 100.);
        fDMaxTickMinTickHist   = dir.make<TH1F>(Form("DMaxTMinT_%1zu",    fPlane), ";Delta Max/Min Tick;",  100,   0., 100.);
        fDMaxDerivMinDerivHist = dir.make<TH1F>(Form("DMaxDMinD_%1zu",    fPlane), ";Delta Max/Min Deriv;", 200,   0., 100.);
        fMaxErosionHist        = dir.make<TH1F>(Form("MaxErosion_%1zu",   fPlane), ";Max Erosion;",         200, -50., 150.);
        fMaxDilationHist       = dir.make<TH1F>(Form("MaxDilation_%1zu",  fPlane), ";Max Dilation;",        200, -50., 150.);
        fMaxDilEroRatHist      = dir.make<TH1F>(Form("MaxDilEroRat_%1zu", fPlane), ";Max Dil/Ero;",         200,  -1.,   1.);
    }

    return;
}

void CandHitMorphological::findHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t& dataRange,
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

    // Now we get the erosion/dilation vectors
    Waveform erosionVec;
    Waveform dilationVec;
    Waveform averageVec;
    Waveform differenceVec;

    reco_tool::HistogramMap histogramMap;

    // Compute the morphological filter vectors
    fWaveformTool->getErosionDilationAverageDifference(waveform, fStructuringElement, histogramMap, erosionVec, dilationVec, averageVec, differenceVec);

    // Now find the hits
    findHitCandidates(derivativeVec.begin(), derivativeVec.end(),
                      erosionVec.begin(),    erosionVec.end(),
                      dilationVec.begin(),   dilationVec.end(),
                      roiStartTick,
                      fDilationThreshold,
                      hitCandidateVec);

    // Limit start and stop tick to the neighborhood of the peak
    if (fFitNSigmaFromCenter>0) {
      for (auto& hc: hitCandidateVec) {
	auto startCand = hc.hitCenter-fFitNSigmaFromCenter*hc.hitSigma;
        if (startCand>=0) hc.startTick = std::max(hc.startTick, size_t(startCand));
        hc.stopTick = std::min(hc.stopTick, size_t(hc.hitCenter+fFitNSigmaFromCenter*hc.hitSigma));
      }
    }

    // Reset the hit height from the input waveform
    for(auto& hitCandidate : hitCandidateVec)
    {
        size_t centerIdx = hitCandidate.hitCenter;

        hitCandidate.hitHeight = waveform.at(centerIdx);
    }

    // Keep track of histograms if requested
    if (fOutputWaveforms)
    {
        // Recover the details...
        std::vector<geo::WireID> wids  = fGeometry->ChannelToWire(channel);
        size_t                   plane = wids[0].Plane;
        size_t                   cryo  = wids[0].Cryostat;
        size_t                   tpc   = wids[0].TPC;
        size_t                   wire  = wids[0].Wire;

        if (channel != fLastChannel) fChannelCnt = 0;

        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("Event%04zu/c%1zuT%1zuP%1zu/Wire_%05zu",eventCount,cryo,tpc,plane,wire));

        size_t waveformSize = waveform.size();
        size_t waveStart    = dataRange.begin_index();

        TProfile* waveHist     = dir.make<TProfile>(Form("HWfm_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Waveform",   waveformSize, 0, waveformSize, -500., 500.);
        TProfile* derivHist    = dir.make<TProfile>(Form("HDer_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Derivative", waveformSize, 0, waveformSize, -500., 500.);
        TProfile* erosionHist  = dir.make<TProfile>(Form("HEro_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Erosion",    waveformSize, 0, waveformSize, -500., 500.);
        TProfile* dilationHist = dir.make<TProfile>(Form("HDil_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Dilation",   waveformSize, 0, waveformSize, -500., 500.);
        TProfile* candHitHist  = dir.make<TProfile>(Form("HCan_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Cand Hits",  waveformSize, 0, waveformSize, -500., 500.);
        TProfile* maxDerivHist = dir.make<TProfile>(Form("HMax_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Maxima",     waveformSize, 0, waveformSize, -500., 500.);
        TProfile* strtStopHist = dir.make<TProfile>(Form("HSSS_%03zu_roiStart-%05zu",fChannelCnt,waveStart), "Start/Stop", waveformSize, 0, waveformSize, -500., 500.);

        // Fill wave/derivative
        for(size_t idx = 0; idx < waveform.size(); idx++)
        {
            waveHist->Fill(roiStartTick + idx, waveform.at(idx));
            derivHist->Fill(roiStartTick + idx, derivativeVec.at(idx));
            erosionHist->Fill(roiStartTick + idx, erosionVec.at(idx));
            dilationHist->Fill(roiStartTick + idx, dilationVec.at(idx));
        }

        // Fill hits
        for(const auto& hitCandidate : hitCandidateVec)
        {
            candHitHist->Fill(hitCandidate.hitCenter,  hitCandidate.hitHeight);
            maxDerivHist->Fill(hitCandidate.maxTick,   hitCandidate.maxDerivative);
            maxDerivHist->Fill(hitCandidate.minTick,   hitCandidate.minDerivative);
            strtStopHist->Fill(hitCandidate.startTick, waveform.at(hitCandidate.startTick));
            strtStopHist->Fill(hitCandidate.stopTick,  waveform.at(hitCandidate.stopTick));
        }

        fLastChannel = channel;
        fChannelCnt++;
    }

    if (fOutputHistograms)
    {
        // Fill hits
        for(const auto& hitCandidate : hitCandidateVec)
        {
            fDStopStartHist->Fill(hitCandidate.stopTick - hitCandidate.startTick, 1.);
            fDMaxTickMinTickHist->Fill(hitCandidate.minTick - hitCandidate.maxTick, 1.);
            fDMaxDerivMinDerivHist->Fill(hitCandidate.maxDerivative - hitCandidate.minDerivative, 1.);
        }

        // Get the max dilation/erosion
        Waveform::const_iterator maxDilationItr = std::max_element(dilationVec.begin(), dilationVec.end());
        Waveform::const_iterator maxErosionItr  = std::max_element(erosionVec.begin(),  erosionVec.end());

        float dilEroRat(1.);

        if (std::abs(*maxDilationItr) > 0.) dilEroRat = *maxErosionItr / *maxDilationItr;

        fMaxErosionHist->Fill(*maxErosionItr,  1.);
        fMaxDilationHist->Fill(*maxDilationItr, 1.);
        fMaxDilEroRatHist->Fill(dilEroRat, 1.);
    }

    return;
}

void CandHitMorphological::findHitCandidates(Waveform::const_iterator derivStartItr,    Waveform::const_iterator derivStopItr,
                                             Waveform::const_iterator erosionStartItr,  Waveform::const_iterator erosionStopItr,
                                             Waveform::const_iterator dilationStartItr, Waveform::const_iterator dilationStopItr,
                                             const size_t                   roiStartTick,
                                             float                    dilationThreshold,
                                             HitCandidateVec&         hitCandidateVec) const
{
    // This function aims to use the erosion/dilation vectors to find candidate hit regions
    // Once armed with a region then the "standard" differential approach is used to return the candidate peaks

    // Don't do anything if not enough ticks
    int ticksInInputWaveform = std::distance(derivStartItr, derivStopItr);

    if (ticksInInputWaveform < fMinDeltaTicks) return;

    // This function is recursive, we start by finding the largest element of the dilation vector
    Waveform::const_iterator maxItr = std::max_element(dilationStartItr,dilationStopItr);
    float                    maxVal = *maxItr;

    // Check that the peak is of reasonable height...
    if (maxVal < dilationThreshold) return;

    int maxBin = std::distance(dilationStartItr,maxItr);

    // The candidate hit region we want lies between the two nearest minima to the maximum we just found
    // subject to the condition that the erosion vector has return to less than zero
    Waveform::const_iterator firstDerItr = derivStartItr   + maxBin;
    Waveform::const_iterator erosionItr  = erosionStartItr + maxBin;

    float firstDerivValue = -1.;
    float erosionCutValue = fErosionFraction * maxVal;

    // Search for starting point
    while(firstDerItr != derivStartItr)
    {
        // Look for the turnover point
        if (*erosionItr-- < erosionCutValue)
        {
            // We are looking for the zero crossing signifying a minimum value in the waveform
            // (so the previous derivative < 0 while current is > 0)
            // We are moving "backwards" so the current value <= 0, the previous value > 0
            if (*firstDerItr * firstDerivValue <= 0. && firstDerivValue > 0.) break;
        }

        firstDerivValue = *firstDerItr--;
    }

    // Set the start bin
    int hitRegionStart = std::distance(derivStartItr,firstDerItr);

    // Now go the other way
    Waveform::const_iterator lastDerItr = derivStartItr + maxBin;

    // Reset the local variables
    float lastDerivValue = 1.;

    erosionItr = erosionStartItr + maxBin;

    // Search for starting point
    while(lastDerItr != derivStopItr)
    {
        if (*erosionItr++ <= erosionCutValue)
        {
            // We are again looking for the zero crossing signifying a minimum value in the waveform
            // This time we are moving forward, so test is that previous value < 0, new value >= 0
            if (*lastDerItr * lastDerivValue <= 0. && lastDerivValue < 0.) break;
        }

        lastDerivValue = *lastDerItr++;
    }

    // Set the stop bin
    int hitRegionStop = std::distance(derivStartItr,lastDerItr);

    // Recursive call to find any hits in front of where we are now
    if (hitRegionStart > fMinDeltaTicks)
        findHitCandidates(derivStartItr,    derivStartItr    + hitRegionStart,
                          erosionStartItr,  erosionStartItr  + hitRegionStart,
                          dilationStartItr, dilationStartItr + hitRegionStart,
                          roiStartTick,
                          fDilationThreshold,
                          hitCandidateVec);

    // Call the differential hit finding to get the actual hits within the region
    findHitCandidates(derivStartItr + hitRegionStart, derivStartItr + hitRegionStop,
                      roiStartTick  + hitRegionStart,
                      fMinDeltaTicks,
                      fMinDeltaPeaks,
                      hitCandidateVec);

    // Now call ourselves again to find any hits trailing the region we just identified
    if (std::distance(lastDerItr,derivStopItr) > fMinDeltaTicks)
        findHitCandidates(derivStartItr    + hitRegionStop,    derivStopItr,
                          erosionStartItr  + hitRegionStop,    erosionStopItr,
                          dilationStartItr + hitRegionStop,    dilationStopItr,
                          roiStartTick     + hitRegionStop,
                          fDilationThreshold,
                          hitCandidateVec);

    return;
}

void CandHitMorphological::findHitCandidates(Waveform::const_iterator startItr,
                                             Waveform::const_iterator stopItr,
                                             const size_t             roiStartTick,
                                             int                      dTicksThreshold,
                                             float                    dPeakThreshold,
                                             HitCandidateVec&         hitCandidateVec) const
{
    // Search for candidate hits...
    // Strategy is to get the list of all possible max/min pairs of the input derivative vector and then
    // look for candidate hits in that list
    CandHitParamsVec candHitParamsVec;

    if (getListOfHitCandidates(startItr, stopItr, dTicksThreshold, dPeakThreshold, candHitParamsVec))
    {
        // We've been given a list of candidate hits so now convert to hits
        // Version one... simply convert all the candidates
        for(const auto& tuple : candHitParamsVec)
        {
            // Create a new hit candidate and store away
            HitCandidate hitCandidate;

            Waveform::const_iterator candStartItr = std::get<0>(tuple);
            Waveform::const_iterator maxItr       = std::get<1>(tuple);
            Waveform::const_iterator minItr       = std::get<2>(tuple);
            Waveform::const_iterator candStopItr  = std::get<3>(tuple);

            Waveform::const_iterator peakItr = std::min_element(maxItr,minItr,[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});

            // Check balance
            if      (2 * std::distance(peakItr,minItr) < std::distance(maxItr,peakItr)) peakItr--;
            else if (2 * std::distance(maxItr,peakItr) < std::distance(peakItr,minItr)) peakItr++;

            // Special handling of the start tick for multiple hits
            size_t hitCandidateStartTick = roiStartTick + std::distance(startItr,candStartItr);

            if (!hitCandidateVec.empty())
            {
                int deltaTicks = hitCandidateStartTick - hitCandidateVec.back().stopTick;

                if (deltaTicks > 0)
                {
                    hitCandidateStartTick           -= deltaTicks / 2;
                    hitCandidateVec.back().stopTick += deltaTicks / 2;
                }
            }

            hitCandidate.startTick     = hitCandidateStartTick;
            hitCandidate.stopTick      = roiStartTick + std::distance(startItr,candStopItr);
            hitCandidate.maxTick       = roiStartTick + std::distance(startItr,maxItr);
            hitCandidate.minTick       = roiStartTick + std::distance(startItr,minItr);
            hitCandidate.maxDerivative = maxItr != stopItr ? *maxItr : 0.;
            hitCandidate.minDerivative = minItr != stopItr ? *minItr : 0.;
            hitCandidate.hitCenter     = roiStartTick + std::distance(startItr,peakItr) + 0.5;
            hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
            hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130;

            hitCandidateVec.push_back(hitCandidate);
        }
    }

//    // The idea will be to find the largest deviation in the input derivative waveform as the starting point. Depending
//    // on if a maximum or minimum, we search forward or backward to find the minimum or maximum that our extremum
//    // corresponds to.
//    std::pair<Waveform::const_iterator, Waveform::const_iterator> minMaxPair = std::minmax_element(startItr, stopItr);
//
//    Waveform::const_iterator maxItr = minMaxPair.second;
//    Waveform::const_iterator minItr = minMaxPair.first;
//
//    // Use the larger of the two as the starting point and recover the nearest max or min
//    if (std::fabs(*maxItr) > std::fabs(*minItr)) minItr = findNearestMin(maxItr, stopItr);
//    else                                         maxItr = findNearestMax(minItr,startItr);
//
//    int   deltaTicks = std::distance(maxItr,minItr);
//    float range      = *maxItr - *minItr;
//
//    // At some point small rolling oscillations on the waveform need to be ignored...
//    if (deltaTicks >= dTicksThreshold && range > dPeakThreshold)
//    {
//        // Need to back up to find zero crossing, this will be the starting point of our
//        // candidate hit but also the endpoint of the pre sub-waveform we'll search next
//        Waveform::const_iterator newEndItr = findStartTick(maxItr, startItr);
//
//        int startTick = std::distance(startItr,newEndItr);
//
//        // Now need to go forward to again get close to zero, this will then be the end point
//        // of our candidate hit and the starting point for the post sub-waveform to search
//        Waveform::const_iterator newStartItr = findStopTick(minItr, stopItr);
//
//        int stopTick = std::distance(startItr,newStartItr);
//
//        // Find hits in the section of the waveform leading up to this candidate hit
//        if (startTick > 2)
//        {
//            // Special handling for merged hits
//            if (*(newEndItr-1) > 0.) {dTicksThreshold = 2;              dPeakThreshold = 0.;            }
//            else                     {dTicksThreshold = fMinDeltaTicks; dPeakThreshold = fMinDeltaPeaks;}
//
//            findHitCandidates(startItr,newEndItr+1,roiStartTick,dTicksThreshold,dPeakThreshold,hitCandidateVec);
//        }
//
//        // Create a new hit candidate and store away
//        HitCandidate_t hitCandidate;
//
//        Waveform::const_iterator peakItr = std::min_element(maxItr,minItr,[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
//
//        // Check balance
//        if      (2 * std::distance(peakItr,minItr) < std::distance(maxItr,peakItr)) peakItr--;
//        else if (2 * std::distance(maxItr,peakItr) < std::distance(peakItr,minItr)) peakItr++;
//
//        // Special handling of the start tick for multiple hits
//        size_t hitCandidateStartTick = roiStartTick + startTick;
//
//        if (!hitCandidateVec.empty())
//        {
//            int deltaTicks = hitCandidateStartTick - hitCandidateVec.back().stopTick;
//
//            if (deltaTicks > 0)
//            {
//                hitCandidateStartTick           -= deltaTicks / 2;
//                hitCandidateVec.back().stopTick += deltaTicks / 2;
//            }
//        }
//
//        hitCandidate.startTick     = hitCandidateStartTick;
//        hitCandidate.stopTick      = roiStartTick + stopTick;
//        hitCandidate.maxTick       = roiStartTick + std::distance(startItr,maxItr);
//        hitCandidate.minTick       = roiStartTick + std::distance(startItr,minItr);
//        hitCandidate.maxDerivative = maxItr != stopItr ? *maxItr : 0.;
//        hitCandidate.minDerivative = minItr != stopItr ? *minItr : 0.;
//        hitCandidate.hitCenter     = roiStartTick + std::distance(startItr,peakItr) + 0.5;
//        hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
//        hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130;
//
//        hitCandidateVec.push_back(hitCandidate);
//
//        // Finally, search the section of the waveform following this candidate for more hits
//        if (std::distance(newStartItr,stopItr) > 2)
//        {
//            // Special handling for merged hits
//            if (*(newStartItr+1) < 0.) {dTicksThreshold = 2;              dPeakThreshold = 0.;            }
//            else                       {dTicksThreshold = fMinDeltaTicks; dPeakThreshold = fMinDeltaPeaks;}
//
//            findHitCandidates(newStartItr,stopItr,roiStartTick + stopTick,dTicksThreshold,dPeakThreshold,hitCandidateVec);
//        }
//    }

    return;
}

bool CandHitMorphological::getListOfHitCandidates(Waveform::const_iterator startItr,
                                                  Waveform::const_iterator stopItr,
                                                  int                      dTicksThreshold,
                                                  float                    dPeakThreshold,
                                                  CandHitParamsVec&        candHitParamsVec) const
{
    // We'll check if any of our candidates meet the requirements so declare the result here
    bool foundCandidate(false);

    int dTicks = std::distance(startItr,stopItr);

    // Search for candidate hits...
    // But only if enough ticks
    if (dTicks < fMinDeltaTicks) return foundCandidate;

    // Generally, the mission is simple... the goal is to find all possible combinations of maximum/minimum pairs in
    // the input (presumed) derivative waveform. We can do this with a divice and conquer approach where we start by
    // finding the largerst max or min and start from there
    MaxMinPair minMaxPair = std::minmax_element(startItr, stopItr);

    Waveform::const_iterator maxItr = minMaxPair.second;
    Waveform::const_iterator minItr = minMaxPair.first;

    // Use the larger of the two as the starting point and recover the nearest max or min
    if (std::fabs(*maxItr) > std::fabs(*minItr)) minItr = findNearestMin(maxItr, stopItr);
    else                                         maxItr = findNearestMax(minItr,startItr);

    int   deltaTicks = std::distance(maxItr,minItr);
    float range      = *maxItr - *minItr;

    if (deltaTicks < 2) return foundCandidate;

    // Check if this particular max/min pair would meet the requirements...
    if (deltaTicks >= dTicksThreshold && range > dPeakThreshold) foundCandidate = true;

    // Need to back up to find zero crossing, this will be the starting point of our
    // candidate hit but also the endpoint of the pre sub-waveform we'll search next
    Waveform::const_iterator candStartItr = findStartTick(maxItr, startItr);

    // Now need to go forward to again get close to zero, this will then be the end point
    // of our candidate hit and the starting point for the post sub-waveform to search
    Waveform::const_iterator candStopItr = findStopTick(minItr, stopItr);

    // Call ourself to find hit candidates preceding this one
    bool prevTicks = getListOfHitCandidates(startItr, candStartItr, dTicksThreshold, dPeakThreshold, candHitParamsVec);

    // The above call will have populated the list of candidate max/min pairs preceding this one, so now add our contribution
    candHitParamsVec.emplace_back(candStartItr, maxItr, minItr, candStopItr);

    // Now catch any that might follow this one
    bool postTicks = getListOfHitCandidates(candStopItr, stopItr, dTicksThreshold, dPeakThreshold, candHitParamsVec);

    return foundCandidate || prevTicks || postTicks;
}


void CandHitMorphological::MergeHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t& rangeData,
                                              const HitCandidateVec& hitCandidateVec,
                                              MergeHitCandidateVec&  mergedHitsVec) const
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

ICandidateHitFinder::Waveform::const_iterator CandHitMorphological::findNearestMin(Waveform::const_iterator maxItr, Waveform::const_iterator stopItr) const
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

ICandidateHitFinder::Waveform::const_iterator CandHitMorphological::findNearestMax(Waveform::const_iterator minItr, Waveform::const_iterator startItr) const
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

ICandidateHitFinder::Waveform::const_iterator CandHitMorphological::findStartTick(Waveform::const_iterator maxItr, Waveform::const_iterator startItr) const
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

ICandidateHitFinder::Waveform::const_iterator CandHitMorphological::findStopTick(Waveform::const_iterator minItr, Waveform::const_iterator stopItr)   const
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

DEFINE_ART_CLASS_TOOL(CandHitMorphological)
}
