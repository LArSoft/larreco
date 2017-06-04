///////////////////////////////////////////////////////////////////////
///
/// \file   ICandidateHitFinder.h
///
/// \brief  This provides an interface for tools which are tasked with
///         finding candidate hits on input waveforms
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef ICandidateHitFinder_H
#define ICandidateHitFinder_H

#include "fhiclcpp/ParameterSet.h"

namespace reco_tool
{
    class ICandidateHitFinder
    {
    public:
        virtual ~ICandidateHitFinder() noexcept = default;
        
        // Define standard art tool interface
        virtual void configure(const fhicl::ParameterSet& pset) = 0;
        
        // Define a structure to contain hits
        using HitCandidate_t = struct HitCandidate
        {
            size_t startTick;
            size_t stopTick;
            size_t maxTick;
            size_t minTick;
            float  maxDerivative;
            float  minDerivative;
            float  hitCenter;
            float  hitSigma;
            float  hitHeight;
        };
        
        using HitCandidateVec      = std::vector<HitCandidate_t>;
        using MergeHitCandidateVec = std::vector<HitCandidateVec>;
        
        // Search for candidate hits on the input waveform
        virtual void findHitCandidates(const std::vector<float>&,           // Waveform to analyze
                                       size_t,                              // waveform start tick
                                       double,                              // threshold
                                       HitCandidateVec&) const = 0;         // output candidate hits
        
        // Search for candidate hits on the input waveform
        virtual void findHitCandidates(std::vector<float>::const_iterator,  // Start of waveform
                                       std::vector<float>::const_iterator,  // end of waveform
                                       size_t,                              // waveform start tick
                                       double,                              // threshold
                                       HitCandidateVec&) const = 0;         // output candidate hits
        
        virtual void MergeHitCandidates(const std::vector<float>&,
                                        const HitCandidateVec&,
                                        MergeHitCandidateVec&) const = 0;
    };
}

#endif
