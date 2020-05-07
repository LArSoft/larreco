////////////////////////////////////////////////////////////////////////
/// \file   PeakFitterGaussElimination.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"

#include "art/Utilities/ToolMacros.h"
#include "larreco/HitFinder/GaussianEliminationAlg.h"

#include <algorithm>

namespace reco_tool
{

class PeakFitterGaussElimination : IPeakFitter
{
public:
    explicit PeakFitterGaussElimination(const fhicl::ParameterSet& pset);

    void findPeakParameters(const std::vector<float>&,
                            const ICandidateHitFinder::HitCandidateVec&,
                            PeakParamsVec&,
                            double&,
                            int&) const override;

private:
    // Member variables from the fhicl file
    float fStepSize;     ///< Step size used by gaussian elim alg
    float fMax;          ///< Max

    std::unique_ptr<util::GaussianEliminationAlg> fGEAlg;
};

//----------------------------------------------------------------------
// Constructor.
PeakFitterGaussElimination::PeakFitterGaussElimination(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fStepSize = pset.get<float>("StepSize", 0.1);
    fMax      = pset.get<float>("Max",      0.5);

//    fGEAlg    = std::make_unique<util::GaussianEliminationAlg>(fStepSize, fMax);

    return;
}

// --------------------------------------------------------------------------------------------
void PeakFitterGaussElimination::findPeakParameters(const std::vector<float>&                   roiSignalVec,
                                                    const ICandidateHitFinder::HitCandidateVec& hitCandidateVec,
                                                    PeakParamsVec&                              peakParamsVec,
                                                    double&                                     chi2PerNDF,
                                                    int&                                        NDF) const
{
    // This module tries to use the method for fitting hits found in the RRFHitFinder
    // from Wes Ketchum. It uses the gaussian elimation algorithm he set up.
    //
    // *** NOTE: this algorithm assumes the reference time for input hit candidates is to
    //           the first tick of the input waveform (ie 0)
    //
    if (hitCandidateVec.empty()) return;

    std::vector<float> meanVec;
    std::vector<float> sigmaVec;
    std::vector<float> heightVec;

    for(const auto& hitCandidate : hitCandidateVec)
    {
        float  candMean   = hitCandidate.hitCenter;
        float  candSigma  = hitCandidate.hitSigma;
        size_t bin        = std::floor(candMean);

        bin = std::min(bin, roiSignalVec.size() - 1);

        float  candHeight = roiSignalVec[bin] - (candMean-(float)bin)*(roiSignalVec[bin]-roiSignalVec[bin+1]);

        meanVec.push_back(candMean);
        sigmaVec.push_back(candSigma);
        heightVec.push_back(candHeight);
    }

    return;
}

DEFINE_ART_CLASS_TOOL(PeakFitterGaussElimination)
}
