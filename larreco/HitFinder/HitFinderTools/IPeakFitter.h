///////////////////////////////////////////////////////////////////////
///
/// \file   IPeakFitter.h
///
/// \brief  This provides an interface for tools which are tasked with
///         fitting peaks on input waveforms
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IPeakFitter_H
#define IPeakFitter_H

#include "fhiclcpp/ParameterSet.h"
#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

namespace reco_tool
{
    class IPeakFitter
    {
    public:
        // Define standard art tool interface
        // virtual void configure(const fhicl::ParameterSet& pset) = 0;

        // Define a structure to contain hits
        struct PeakFitParams_t
        {
            float peakCenter;
            float peakCenterError;
            float peakSigma;
            float peakSigmaError;
            float peakAmplitude;
            float peakAmplitudeError;
        };

        using PeakParamsVec = std::vector<PeakFitParams_t>;
        virtual ~IPeakFitter() = default;
        // Get parameters for input candidate peaks
        virtual void findPeakParameters(const std::vector<float>&,
                                        const ICandidateHitFinder::HitCandidateVec&,
                                        PeakParamsVec&,
                                        double&,
                                        int&) const = 0;
    };
}

#endif
