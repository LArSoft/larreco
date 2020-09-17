////////////////////////////////////////////////////////////////////////
/// \file   PeakFitterGaussian.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"
#include "larreco/RecoAlg/GausFitCache.h" // hit::GausFitCache

#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/Geometry.h"

#include <cassert>
#include <fstream>

#include "TF1.h"
#include "TH1F.h"

namespace reco_tool
{

/// Customized function cache for Gaussians with a baseline.
///
/// The baseline parameter is always the last one.
class BaselinedGausFitCache: public hit::GausFitCache {

public:
    /// Constructor (see base class constructor).
    BaselinedGausFitCache(std::string const& new_name="BaselinedGausFitCache")
      : hit::GausFitCache(new_name)
      {}

protected:

    /// Creates and returns the function with specified number of Gaussians.
    ///
    /// The formula is `gaus(0) + gaus(3) + ... + gaus(3*(nFunc-1)) + [nFunc*3]`.
    virtual TF1* CreateFunction(size_t nFunc) const
    {
         // add the Gaussians first
         std::string formula;
         std::size_t iGaus = 0;
         while (iGaus < nFunc) formula += "gaus(" + std::to_string(3 * (iGaus++)) + ") + ";
         formula += "[" + std::to_string(3 * nFunc) + "]";

         std::string const func_name = FunctionName(nFunc);
         auto* pF = new TF1(func_name.c_str(), formula.c_str());
         pF->SetParName(iGaus * 3, "baseline");
         return pF;
    } // CreateFunction()

}; // BaselinedGausFitCache


class PeakFitterGaussian : IPeakFitter
{
public:
    explicit PeakFitterGaussian(const fhicl::ParameterSet& pset);

    void findPeakParameters(const std::vector<float>&,
                            const ICandidateHitFinder::HitCandidateVec&,
                            PeakParamsVec&,
                            double&,
                            int&) const override;

private:
    // Member variables from the fhicl file
    const double                   fMinWidth;          ///< minimum initial width for gaussian fit
    const double                   fMaxWidthMult;      ///< multiplier for max width for gaussian fit
    const double                   fPeakRange;         ///< set range limits for peak center
    const double                   fAmpRange;          ///< set range limit for peak amplitude
    const bool                     fFloatBaseline;     ///< Allow baseline to "float" away from zero
    const bool                     fOutputHistograms;  ///< If true will generate summary style histograms

    TH1F*                    fNumCandHitsHist;
    TH1F*                    fROISizeHist;
    TH1F*                    fCandPeakPositionHist;
    TH1F*                    fCandPeakWidHist;
    TH1F*                    fCandPeakAmpitudeHist;
    TH1F*                    fCandBaselineHist;
    TH1F*                    fFitPeakPositionHist;
    TH1F*                    fFitPeakWidHist;
    TH1F*                    fFitPeakAmpitudeHist;
    TH1F*                    fFitBaselineHist;

    mutable BaselinedGausFitCache fFitCache; ///< Preallocated ROOT functions for the fits.

    mutable TH1F             fHistogram;

    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};

//----------------------------------------------------------------------
// Constructor.
PeakFitterGaussian::PeakFitterGaussian(const fhicl::ParameterSet& pset):
    fMinWidth(pset.get<double>("MinWidth",         0.5)),
    fMaxWidthMult (pset.get<double>("MaxWidthMult",     3.)),
    fPeakRange(pset.get<double>("PeakRangeFact",    2.)),
    fAmpRange(pset.get<double>("PeakAmpRange",     2.)),
    fFloatBaseline(pset.get< bool >("FloatBaseline",    false)),
    fOutputHistograms(pset.get< bool >("OutputHistograms", false))
{
    fHistogram = TH1F("PeakFitterHitSignal","",500,0.,500.);

    fHistogram.Sumw2();

    std::string function = "Gaus(0)";

    // If asked, define the global histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;

        // Make a directory for these histograms
        art::TFileDirectory dir = tfs->mkdir("PeakFit");

        fNumCandHitsHist      = dir.make<TH1F>("NumCandHits",    "# Candidate Hits", 100,   0., 100.);
        fROISizeHist          = dir.make<TH1F>("ROISize",        "ROI Size",         400,   0., 400.);
        fCandPeakPositionHist = dir.make<TH1F>("CPeakPosition",  "Peak Position",    200,   0., 400.);
        fCandPeakWidHist      = dir.make<TH1F>("CPeadWidth",     "Peak Width",       100,   0.,  25.);
        fCandPeakAmpitudeHist = dir.make<TH1F>("CPeakAmplitude", "Peak Amplitude",   100,   0., 200.);
        fCandBaselineHist     = dir.make<TH1F>("CBaseline",      "Baseline",         200, -25.,  25.);
        fFitPeakPositionHist  = dir.make<TH1F>("FPeakPosition",  "Peak Position",    200,   0., 400.);
        fFitPeakWidHist       = dir.make<TH1F>("FPeadWidth",     "Peak Width",       100,   0.,  25.);
        fFitPeakAmpitudeHist  = dir.make<TH1F>("FPeakAmplitude", "Peak Amplitude",   100,   0., 200.);
        fFitBaselineHist      = dir.make<TH1F>("FBaseline",      "Baseline",         200, -25.,  25.);
    }

    return;
}

// --------------------------------------------------------------------------------------------
void PeakFitterGaussian::findPeakParameters(const std::vector<float>&                   roiSignalVec,
                                            const ICandidateHitFinder::HitCandidateVec& hitCandidateVec,
                                            PeakParamsVec&                              peakParamsVec,
                                            double&                                     chi2PerNDF,
                                            int&                                        NDF) const
{
    // The following is a translation of the original FitGaussians function in the original
    // GausHitFinder module originally authored by Jonathan Asaadi
    //
    // *** NOTE: this algorithm assumes the reference time for input hit candidates is to
    //           the first tick of the input waveform (ie 0)
    //
    if (hitCandidateVec.empty()) return;

    // in case of a fit failure, set the chi-square to infinity
    chi2PerNDF = std::numeric_limits<double>::infinity();

    int startTime = hitCandidateVec.front().startTick;
    int endTime   = hitCandidateVec.back().stopTick;
    int roiSize   = endTime - startTime;

    // Check to see if we need a bigger histogram for fitting
    if (roiSize > fHistogram.GetNbinsX())
    {
        std::string histName = "PeakFitterHitSignal_" + std::to_string(roiSize);
        fHistogram = TH1F(histName.c_str(),"",roiSize,0.,roiSize);
        fHistogram.Sumw2();
    }

    for(int idx = 0; idx < roiSize; idx++) fHistogram.SetBinContent(idx+1,roiSignalVec[startTime+idx]);

    // Build the string to describe the fit formula
#if 0
    std::string equation = "gaus(0)";

    for(size_t idx = 1; idx < hitCandidateVec.size(); idx++) equation += "+gaus(" + std::to_string(3*idx) + ")";

    // Set the baseline if so desired
    float baseline(0.);

    if (fFloatBaseline)
    {
        baseline = roiSignalVec[startTime];

        equation += "+" + std::to_string(baseline);
    }

    // Now define the complete function to fit
    TF1 Gaus("Gaus",equation.c_str(),0,roiSize,TF1::EAddToList::kNo);
#else
    unsigned int const nGaus = hitCandidateVec.size();
    assert(fFitCache.Get(nGaus));
    TF1& Gaus = *(fFitCache.Get(nGaus));

    // Set the baseline if so desired
    float baseline(0.);

    if (fFloatBaseline)
    {
        baseline = roiSignalVec[startTime];
        Gaus.SetParameter(nGaus * 3, baseline);
        Gaus.SetParLimits( nGaus * 3, baseline - 12., baseline + 12.);
    }
    else Gaus.FixParameter(nGaus * 3, baseline); // last parameter is the baseline

#endif // 0

    if (fOutputHistograms)
    {
        fNumCandHitsHist->Fill(hitCandidateVec.size(), 1.);
        fROISizeHist->Fill(roiSize, 1.);
        fCandBaselineHist->Fill(baseline, 1.);
    }

    // ### Setting the parameters for the Gaussian Fit ###
    int parIdx{0};
    for(auto const& candidateHit : hitCandidateVec)
    {
        double const peakMean   = candidateHit.hitCenter - float(startTime);
        double const peakWidth  = candidateHit.hitSigma;
        double const amplitude  = candidateHit.hitHeight - baseline;
        double const meanLowLim = std::max(peakMean - fPeakRange * peakWidth,              0.);
        double const meanHiLim  = std::min(peakMean + fPeakRange * peakWidth, double(roiSize));

        if (fOutputHistograms)
        {
            fCandPeakPositionHist->Fill(peakMean, 1.);
            fCandPeakWidHist->Fill(peakWidth, 1.);
            fCandPeakAmpitudeHist->Fill(amplitude, 1.);
        }

        Gaus.SetParameter(  parIdx, amplitude);
        Gaus.SetParameter(1+parIdx, peakMean);
        Gaus.SetParameter(2+parIdx, peakWidth);
        Gaus.SetParLimits(  parIdx, 0.1 * amplitude,  fAmpRange * amplitude);
        Gaus.SetParLimits(1+parIdx, meanLowLim,       meanHiLim);
        Gaus.SetParLimits(2+parIdx, std::max(fMinWidth, 0.1 * peakWidth), fMaxWidthMult * peakWidth);

        parIdx += 3;
    }

    int fitResult{-1};

    try
    { fitResult = fHistogram.Fit(&Gaus,"QNWB","", 0., roiSize);}
    catch(...)
    {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}

    // If the fit result is not zero there was an error
    if (!fitResult)
    {
        // ##################################################
        // ### Getting the fitted parameters from the fit ###
        // ##################################################
        chi2PerNDF = (Gaus.GetChisquare() / Gaus.GetNDF());
        NDF        = Gaus.GetNDF();

        int parIdx { 0 };
        for(size_t idx = 0; idx < hitCandidateVec.size(); idx++)
        {
            PeakFitParams_t peakParams;

            peakParams.peakAmplitude      = Gaus.GetParameter(parIdx);
            peakParams.peakAmplitudeError = Gaus.GetParError( parIdx);
            peakParams.peakCenter         = Gaus.GetParameter(parIdx + 1) + float(startTime);
            peakParams.peakCenterError    = Gaus.GetParError( parIdx + 1);
            peakParams.peakSigma          = Gaus.GetParameter(parIdx + 2);
            peakParams.peakSigmaError     = Gaus.GetParError( parIdx + 2);

            if (fOutputHistograms)
            {
                fFitPeakPositionHist->Fill(peakParams.peakCenter, 1.);
                fFitPeakWidHist->Fill(peakParams.peakSigma, 1.);
                fFitPeakAmpitudeHist->Fill(peakParams.peakAmplitude, 1.);
            }

            peakParamsVec.emplace_back(peakParams);

            parIdx += 3;
        }

        if (fOutputHistograms) fFitBaselineHist->Fill(Gaus.GetParameter(3*nGaus), 1.);
    }
#if 0
    Gaus.Delete();
#endif // 0
    return;
}

DEFINE_ART_CLASS_TOOL(PeakFitterGaussian)
}
