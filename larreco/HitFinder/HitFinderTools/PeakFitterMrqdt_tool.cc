#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"
#include "larreco/RecoAlg/GausFitCache.h" // hit::GausFitCache
#include "lardata/Utilities/MarqFitAlg.h"//marqfit functions

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <cassert>
#include <fstream>

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

namespace reco_tool
{

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
    //bool fOutputHistograms; // unused

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

    int startTime = fhc_vec.front().startTick;
    int endTime = fhc_vec.back().stopTick;

    int roiSize = endTime - startTime;

    std::vector<float> y(roiSize);
    std::vector<float> p(3*fhc_vec.size());
    std::vector<float> perr(3*fhc_vec.size());

    /* choose the fit function and set the parameters */
    nParams = 0;

    for(size_t ih=0; ih<fhc_vec.size(); ih++){
      float const peakMean   = fhc_vec[ih].hitCenter - (float)startTime;
      float const peakWidth  = fhc_vec[ih].hitSigma;
      float const amplitude  = fhc_vec[ih].hitHeight;
      //float const meanLowLim = fmax(peakMean - fPeakRange * peakWidth,       0.);
      //float const meanHiLim  = fmin(peakMean + fPeakRange * peakWidth, (float)roiSize);
      p[0+nParams]=amplitude;
      p[1+nParams]=peakMean;
      p[2+nParams]=peakWidth;

      nParams += 3;
    }
    int fitResult=-1;

    //recast signal -> y to keep the  catch about negative adc values
    for(size_t idx=0; idx<size_t(roiSize); idx++){
      // float adc=signal[startTime+idx];
      //if(adc<=0.) adc=0;
      //y[idx]=adc;

      if(signal[startTime+idx]<=0.) y[idx]=0;
      else y[idx]=signal[startTime+idx];
    }

    int trial=0;
    lambda=-1.;   /* initialize lambda on first call */
    do{
      fitResult=fMarqFitAlg->gshf::MarqFitAlg::mrqdtfit(lambda, &p[0], &y[0], nParams, roiSize, chiSqr, dchiSqr);
      trial++;
      if(fitResult||(trial>100))break;
    }
    while (fabs(dchiSqr) >= chiCut);

    if (!fitResult){
      int fitResult2=fMarqFitAlg->gshf::MarqFitAlg::cal_perr(&p[0],&y[0],nParams,roiSize,&perr[0]);
      if (!fitResult2){
	float NDF = roiSize - nParams;
	chi2PerNDF = chiSqr / NDF;
	int parIdx = 0;
	for(size_t i=0;i<fhc_vec.size();i++){

	  /* stand alone method  
	  mhpp_vec[i].peakAmplitude      = p[parIdx + 0];
	  mhpp_vec[i].peakAmplitudeError = perr[parIdx + 0];
	  mhpp_vec[i].peakCenter         = p[parIdx + 1] + 0.5 + float(startTime);
	  mhpp_vec[i].peakCenterError    = perr[parIdx + 1];
	  mhpp_vec[i].peakSigma          = p[parIdx + 2];
	  mhpp_vec[i].peakSigmaError     = perr[parIdx + 2];
	  */
	  
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
