#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"
#include "larreco/RecoAlg/GausFitCache.h" // hit::GausFitCache
//#include "larreco/HitFinder/MarqFitAlg.h"//marqfit functions
#include "MarqFitAlg.h"//marqfit functions

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <cassert>
#include <fstream>

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
//#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"

namespace reco_tool
{
  /*
  class BaselinedGauseFitCache: public hit::GauseFitCache {

  public:
    BaselinedGausFitCache(std::string const & new_name="BaselinedGausFitCache")
      : hit::GausFitCache(new_name)
    {}

  protected:

  }; // BaselinedGausFitCache
  */
  class PeakFitterMrqdt : IPeakFitter
  {
  public:
    explicit PeakFitterMrqdt(const fhicl::ParameterSet& pset);

    ~PeakFitterMrqdt();

    void configure(const fhicl::ParameterSet& pset) override;

    void findPeakParameters(const std::vector<float>&,
			    const ICandidateHitFinder::HitCandidateVec&,
			    PeakParamsVec&,
			    double&,
			    int&) const override;

  private:
    //Variables from the fhicl file
    double fMinWidth;
    double fMaxWidthMult;
    double fPeakRange;
    double fAmpRange;
    bool fFloatBaseline;
    bool fOutputHistograms;

    // std::unique_ptr<marqfitgaus::MarqFitAlg> fMarqFitAlg;
    //    mutable BaselinedGausFitCache  fFitCache;

    std::unique_ptr<gshf::MarqFitAlg> fMarqFitAlg;

    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
  };
  
  //--------------------------
  //constructor

  PeakFitterMrqdt::PeakFitterMrqdt(const fhicl::ParameterSet & pset)
  {
    configure(pset);
  }

  PeakFitterMrqdt::~PeakFitterMrqdt()
  {
  }

  void PeakFitterMrqdt::configure(const fhicl::ParameterSet& pset)
  {
    //get parameters
    fMinWidth = pset.get<double>("MinWidth", 0.5);
    fMaxWidthMult = pset.get<double>("MaxWidthMult", 3.);
    fPeakRange = pset.get<double>("PeakRangeFact", 2.);
    fAmpRange = pset.get<double>("PeakAmpRange", 2.);
    fFloatBaseline = pset.get< bool >("FloatBaseline", false);
    //    fOutputHistograms = pset.get< bool >("OutputHistograms", false);

    return;
  }

  //------------------------
  void PeakFitterMrqdt::findPeakParameters(const std::vector<float>&                      signal, 
					   const ICandidateHitFinder::HitCandidateVec& fhc_vec, 
					   PeakParamsVec&                              mhpp_vec, 
					   double&                                     chi2PerNDF, 
					   int&                                        NDF) const
  {

    if(fhc_vec.empty()) return;

    const float chiCut   = 1e-3;
    float lambda   = 0.001;      /* Marquardt damping parameter */
    float chiSqr = std::numeric_limits<float>::max(), dchiSqr = std::numeric_limits<float>::max();
    int nParams=0;
    float y[1000],p[15],perr[15];//pmin[15],pmax[15],perr[15];

    int startTime = fhc_vec.front().startTick;
    int endTime = fhc_vec.back().stopTick;

    int roiSize = endTime - startTime;

    /* choose the fit function and set the parameters */
    nParams = 0;

    //    consider  moving to  auto loop  like most other things in larsoft? or is this less efficient?
    //    for(auto const& candidateHit : hitCandidateVec)

    for(size_t ih=0;ih<fhc_vec.size();ih++){
      float const peakMean   = fhc_vec[ih].hitCenter - (float)startTime;
      float const peakWidth  = fhc_vec[ih].hitSigma;
      float const amplitude  = fhc_vec[ih].hitHeight;
      //float const meanLowLim = fmax(peakMean - fPeakRange * peakWidth,       0.);
      //float const meanHiLim  = fmin(peakMean + fPeakRange * peakWidth, (float)roiSize);
      p[0+nParams]=amplitude;
      p[1+nParams]=peakMean;
      p[2+nParams]=peakWidth;
      //pmin[0+nParams]=0.1 * amplitude;
      //pmax[0+nParams]=fAmpRange * amplitude;
      //pmin[1+nParams]=meanLowLim;
      //pmax[1+nParams]=meanHiLim;
      //pmin[2+nParams]=fmax(fMinWidth, 0.1 * peakWidth);
      //pmax[2+nParams]=fMaxWidthMult * peakWidth;
      nParams += 3;
    }
    int fitResult=-1;

    //recast signal -> y to keep the  catch about negative adc values
    for(size_t idx=0; idx<signal.size(); idx++){

      if(signal[idx]<=0.) y[idx]=0;
      else y[idx]=signal[idx];

    }

    int trial=0;
    lambda=-1.;   /* initialize lambda on first call */
    do{
      fitResult=fMarqFitAlg->gshf::MarqFitAlg::mrqdtfit(lambda, p, y, nParams, roiSize, chiSqr, dchiSqr);
      //      fitResult=fMarqFitAlg->mrqdtfit(lambda, p, y, nParams, roiSize, chiSqr, dchiSqr);
      trial++;
      if(fitResult||(trial>100))break;
    }
    while (fabs(dchiSqr) >= chiCut);

    //    int fitStat=-1;
    if (!fitResult){
      int fitResult2=fMarqFitAlg->gshf::MarqFitAlg::cal_perr(p,y,nParams,roiSize,perr);
      //int fitResult2=fMarqFitAlg->cal_perr(p,y,nParams,roiSize,perr);
      if (!fitResult2){
	float NDF = roiSize - nParams;
	chi2PerNDF = chiSqr / NDF;
	int parIdx = 0;
	for(size_t i=0;i<fhc_vec.size();i++){
	  
	  PeakFitParams_t mhpp;	  
	  
	  mhpp.peakAmplitude      = p[parIdx + 0];
	  mhpp.peakAmplitudeError = perr[parIdx + 0];
	  mhpp.peakCenter         = p[parIdx + 1] + 0.5 + float(startTime);
	  mhpp.peakCenterError    = perr[parIdx + 1];
	  mhpp.peakSigma          = p[parIdx + 2];
	  mhpp.peakSigmaError     = perr[parIdx + 2];
	  
	  mhpp_vec.emplace_back(mhpp);
	  
	  parIdx += 3;

	}

	//	fitStat=0;
      }
    }
    return;
  }

  DEFINE_ART_CLASS_TOOL(PeakFitterMrqdt)
}
