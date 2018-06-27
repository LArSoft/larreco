////////////////////////////////////////////////////////////////////////
/// \file   PeakFitterGaussian.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <fstream>
#include "TH1F.h"
#include "TF1.h"

namespace reco_tool
{

class PeakFitterGaussian : IPeakFitter
{
public:
    explicit PeakFitterGaussian(const fhicl::ParameterSet& pset);
    
    ~PeakFitterGaussian();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void findPeakParameters(const std::vector<float>&,
                            const ICandidateHitFinder::HitCandidateVec&,
                            PeakParamsVec&,
                            double&,
                            int&) const override;
    
private:
    // Member variables from the fhicl file
    double                   fMinWidth;     ///< minimum initial width for gaussian fit
    double                   fMaxWidthMult; ///< multiplier for max width for gaussian fit
    double                   fPeakRange;    ///< set range limits for peak center
    double                   fAmpRange;     ///< set range limit for peak amplitude
    
    mutable TH1F             fHistogram;
    
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
PeakFitterGaussian::PeakFitterGaussian(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
PeakFitterGaussian::~PeakFitterGaussian()
{
}
    
void PeakFitterGaussian::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fMinWidth     = pset.get<double>("MinWidth",      0.5);
    fMaxWidthMult = pset.get<double>("MaxWidthMult",  3.);
    fPeakRange    = pset.get<double>("PeakRangeFact", 2.);
    fAmpRange     = pset.get<double>("PeakAmpRange",  2.);
    
    fHistogram    = TH1F("PeakFitterHitSignal","",500,0.,500.);
    
    fHistogram.Sumw2();
    
    std::string function = "Gaus(0)";
    
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
    
    for(int idx = 0; idx < roiSize; idx++) fHistogram.SetBinContent(idx+1,roiSignalVec.at(startTime+idx));
    
    // Build the string to describe the fit formula
    std::string equation = "gaus(0)";
    
    for(size_t idx = 1; idx < hitCandidateVec.size(); idx++) equation += "+gaus(" + std::to_string(3*idx) + ")";
    
    // Set the baseline
    float baseline = roiSignalVec.at(startTime);
    
    equation += "+" + std::to_string(baseline);

    // Now define the complete function to fit
    TF1 Gaus("Gaus",equation.c_str(),0,roiSize);
    
    // ### Setting the parameters for the Gaussian Fit ###
    int parIdx(0);
    for(auto& candidateHit : hitCandidateVec)
    {
        double peakMean   = candidateHit.hitCenter - float(startTime);
        double peakWidth  = candidateHit.hitSigma;
        double amplitude  = candidateHit.hitHeight - baseline;
        double meanLowLim = std::max(peakMean - fPeakRange * peakWidth,              0.);
        double meanHiLim  = std::min(peakMean + fPeakRange * peakWidth, double(roiSize));
        
        Gaus.SetParameter(  parIdx, amplitude);
        Gaus.SetParameter(1+parIdx, peakMean);
        Gaus.SetParameter(2+parIdx, peakWidth);
        Gaus.SetParLimits(  parIdx, 0.1 * amplitude,  fAmpRange * amplitude);
        Gaus.SetParLimits(1+parIdx, meanLowLim,       meanHiLim);
        Gaus.SetParLimits(2+parIdx, std::max(fMinWidth, 0.1 * peakWidth), fMaxWidthMult * peakWidth);
        
        parIdx += 3;
    }
    
    int fitResult(-1);
    
    try
    { fitResult = fHistogram.Fit(&Gaus,"QNRWB","", 0., roiSize);}
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
    
        parIdx = 0;
        for(size_t idx = 0; idx < hitCandidateVec.size(); idx++)
        {
            PeakFitParams_t peakParams;
        
            peakParams.peakAmplitude      = Gaus.GetParameter(parIdx);
            peakParams.peakAmplitudeError = Gaus.GetParError( parIdx);
            peakParams.peakCenter         = Gaus.GetParameter(parIdx + 1) + float(startTime);
            peakParams.peakCenterError    = Gaus.GetParError( parIdx + 1);
            peakParams.peakSigma          = Gaus.GetParameter(parIdx + 2);
            peakParams.peakSigmaError     = Gaus.GetParError( parIdx + 2);
        
            peakParamsVec.emplace_back(peakParams);
        
            parIdx += 3;
        }
    }
    
    Gaus.Delete();
    
    return;
}

DEFINE_ART_CLASS_TOOL(PeakFitterGaussian)
}
