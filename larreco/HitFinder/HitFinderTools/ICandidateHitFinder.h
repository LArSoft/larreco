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
#include "lardataobj/RecoBase/Wire.h"

namespace reco_tool
{
    class ICandidateHitFinder
    {
    public:
        virtual ~ICandidateHitFinder() noexcept = default;

        // Define a structure to contain hits
        struct HitCandidate
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

        using HitCandidateVec      = std::vector<HitCandidate>;
        using MergeHitCandidateVec = std::vector<HitCandidateVec>;

        using Waveform = std::vector<float>;

        // Search for candidate hits on the input waveform
        virtual void findHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t&, // Waveform (with range info) to analyze
                                       const size_t,                                               // waveform start tick
                                       const size_t,                                               // channel #
                                       const size_t,                                               // Event count (for histograms)
                                       HitCandidateVec&) const = 0;                          // output candidate hits

        virtual void MergeHitCandidates(const recob::Wire::RegionsOfInterest_t::datarange_t&,
                                        const HitCandidateVec&,
                                        MergeHitCandidateVec&) const = 0;
    };
}

#endif
