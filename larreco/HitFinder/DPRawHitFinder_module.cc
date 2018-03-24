#ifndef DPRAWHITFINDER_H
#define DPRAWHITFINDER_H

////////////////////////////////////////////////////////////////////////
//
// DPRawHitFinder class
//
// christoph.alt@cern.ch
//
// This algorithm is designed to find hits on raw waveforms in collection planes (dual phase/single phase)
// It is based on GausHitFinder.
// -----------------------------------
// 
// 1. The algorithm walks along the wire and looks for pulses above threshold.
// 2. Depending on distance and minimum ADC count between peaks inside the same ROI,
// peaks can be grouped. Grouped peaks are fitted together (see later).
// 3. Inside each group of peaks, look for pairs of peaks that are very close, with
// one peak being much smaller than the other one (in integral and amplitude).
// If such a pair of peaks is found, merge the two peaks to one peak.
//
// For pulse trains with #peaks <= MaxMultiHit and width < MaxGroupLength:
// 4. Fit n double exponentials to each group of peaks, where n is the number
// of peaks inside this group. 
// 5. If the Chi2/NDF returned > Chi2NDFRetry, attempt to fit n+1 double exponentials
// to the group of peaks by adding a peak close to the maximum deviation between
// fit and waveform. If this is a better fit it then uses the parameters of this
// fit to characterize the "hit" object. If not, try n+2 exponentials and so on.
// Stop when Chi2/NDF is good or 3 times the number of inital exponentials is reached.
//
// If Chi2/NDF is still bad or if #peaks > MaxMultiHit or width > MaxGroupLength:
// 6. Split pulse into equally long hits.
//
// The parameters of the fit are saved in a feature vector by using MVAWriter to
// draw the fitted function in the event display.
//
////////////////////////////////////////////////////////////////////////


// C/C++ standard library
#include <algorithm> // std::accumulate()
#include <vector>
#include <string>
#include <memory> // std::unique_ptr()
#include <utility> // std::move()
#include <cmath>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "fhiclcpp/ParameterSet.h"


// LArSoft Includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/ArtDataHelper/MVAWriter.h"

// ROOT Includes
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TStopwatch.h"

namespace hit{
  class DPRawHitFinder : public art::EDProducer {
    
  public:
    
    explicit DPRawHitFinder(fhicl::ParameterSet const& pset); 
         
    void produce(art::Event& evt) override;
    void beginJob() override;
    void endJob() override;
    void reconfigure(fhicl::ParameterSet const& p) ;

  private:

    using TimeValsVec      = std::vector<std::tuple<int,int,int>>;
    using PeakTimeWidVec   = std::vector<std::tuple<int,int,int,int>>;
    using MergedTimeWidVec = std::vector<std::tuple<int,int,PeakTimeWidVec>>;
    using PeakDevVec	   = std::vector<std::tuple<double,int,int,int>>;
    using ParameterVec 	   = std::vector<std::pair<double,double>>;  //< parameter/error vec 
   
    void findCandidatePeaks(std::vector<float>::const_iterator startItr,
                            std::vector<float>::const_iterator stopItr,
                            TimeValsVec&                       timeValsVec,
                            float&                             PeakMin,
                            int                                firstTick) const;
    
    void mergeCandidatePeaks(const std::vector<float>&, TimeValsVec&, MergedTimeWidVec&) const;

    // ### This function will fit N-Exponentials to a TH1D where N is set ###
    // ###            by the number of peaks found in the pulse         ###
      
    void FitExponentials(const std::vector<float> fSignalVector,
                         const PeakTimeWidVec     fPeakVals,
                         int                      fStartTime,
                         int                      fEndTime,
                         ParameterVec&            fparamVec,
                         double&                  fchi2PerNDF,
                         int&                     fNDF);
    
    void FindPeakWithMaxDeviation(const std::vector<float> fSignalVector,
			  	  int			   fNPeaks,
                          	  int                      fStartTime,
                          	  int                      fEndTime,
                          	  ParameterVec             fparamVec,
                         	  PeakTimeWidVec           fpeakVals,
			  	  PeakDevVec& 		   fPeakDev);

    std::string CreateFitFunction(int fNPeaks);

    void AddPeak(std::tuple<double,int,int,int> fPeakDevCand,
		 PeakTimeWidVec& fpeakValsTemp);

    void SplitPeak(std::tuple<double,int,int,int> fPeakDevCand,
		   PeakTimeWidVec& fpeakValsTemp);

    double WidthFunc(double fPeakMean,
		   double fPeakAmp,
		   double fPeakTau1,
		   double fPeakTau2,
		   double fStartTime,
		   double fEndTime,
		   double fPeakMeanTrue);

    double ChargeFunc(double fPeakMean,
		      double fPeakAmp,
		      double fPeakTau1,
		      double fPeakTau2,
		      double fChargeNormFactor,
		      double fPeakMeanTrue);

    void FillOutHitParameterVector(const std::vector<double>& input,
				   std::vector<double>& output);
      
    void doBinAverage(const std::vector<float>& inputVec,
                      std::vector<float>&       outputVec,
                      size_t                    binsToAverage) const;
      
    void reBin(const std::vector<float>& inputVec,
               std::vector<float>&       outputVec,
               size_t                    nBinsToCombine) const;


    struct Comp {
    // important: we need two overloads, because the comparison
    // needs to be done both ways to check for equality

    bool operator()(std::tuple<int,int,int,int> p, int s) const
    { return std::get<0>(p) < s; }
    bool operator()(int s, std::tuple<int,int,int,int> p) const
    { return s < std::get<0>(p); }
    };
    
    std::string      fCalDataModuleLabel;

    //FHiCL parameter (see "hitfindermodules.fcl" for details)
    int    fLogLevel;
    float  fMinSig;
    int    fTicksToStopPeakFinder;
    int    fMinWidth;
    double fMinADCSum;
    double fMinADCSumOverWidth;
    int    fMaxMultiHit;
    int    fMaxGroupLength;
    double fChargeNorm;
    bool   fTryNplus1Fits;
    double fChi2NDFRetry;
    double fChi2NDFRetryFactorMultiHits;
    double fChi2NDFMax;
    double fChi2NDFMaxFactorMultiHits;
    size_t fNumBinsToAverage;
    double fMinTau;
    double fMaxTau;      
    double fFitPeakMeanRange;
    int    fGroupMaxDistance;
    double fGroupMinADC;
    bool   fDoMergePeaks;
    double fMergeADCSumThreshold;
    double fMergeMaxADCThreshold;
    double fWidthNormalization;
    int    fLongMaxHits;
    int    fLongPulseWidth;

    art::InputTag fNewHitsTag;              // tag of hits produced by this module, need to have it for fit parameter data products 
    anab::FVectorWriter<3> fHitParamWriter; // helper for saving hit fit parameters in data products
    
    TH1F* fFirstChi2;
    TH1F* fChi2;
		
  protected: 
    
  
  }; // class DPRawHitFinder
  

//-------------------------------------------------
//-------------------------------------------------
DPRawHitFinder::DPRawHitFinder(fhicl::ParameterSet const& pset) :
	fNewHitsTag(
	    pset.get<std::string>("module_label"), "",
	    art::ServiceHandle<art::TriggerNamesService>()->getProcessName()),
	fHitParamWriter(this)
{
    this->reconfigure(pset);
  
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);

    // declare that data products with feature vectors describing
    // hits is going to be produced
    fHitParamWriter.produces_using< recob::Hit >();

} // DPRawHitFinder::DPRawHitFinder()


//-------------------------------------------------
//-------------------------------------------------
void DPRawHitFinder::FillOutHitParameterVector(const std::vector<double>& input,
                                              std::vector<double>&       output)
{
    if(input.size()==0)
        throw std::runtime_error("DPRawHitFinder::FillOutHitParameterVector ERROR! Input config vector has zero size.");

    art::ServiceHandle<geo::Geometry> geom;
    const unsigned int N_PLANES = geom->Nplanes();

    if(input.size()==1)
        output.resize(N_PLANES,input[0]);
    else if(input.size()==N_PLANES)
        output = input;
    else
        throw std::runtime_error("DPRawHitFinder::FillOutHitParameterVector ERROR! Input config vector size !=1 and !=N_PLANES.");
      
}
						
  
//-------------------------------------------------
//-------------------------------------------------
void DPRawHitFinder::reconfigure(fhicl::ParameterSet const& p)
{
    // Implementation of optional member function here.
    fLogLevel              	 = p.get< int >("LogLevel");
    fCalDataModuleLabel    	 = p.get< std::string  >("CalDataModuleLabel");
    fMaxMultiHit      	   	 = p.get< int    >("MaxMultiHit");
    fMaxGroupLength	   	 = p.get< int    >("MaxGroupLength");
    fTryNplus1Fits    	   	 = p.get< bool   >("TryNplus1Fits");
    fChi2NDFRetry     	   	 = p.get< double >("Chi2NDFRetry");
    fChi2NDFRetryFactorMultiHits = p.get< double >("Chi2NDFRetryFactorMultiHits");
    fChi2NDFMax        	   	 = p.get< double >("Chi2NDFMax");
    fChi2NDFMaxFactorMultiHits 	 = p.get< double >("Chi2NDFMaxFactorMultiHits");
    fNumBinsToAverage 	   	 = p.get< size_t >("NumBinsToAverage");
    fMinSig           	   	 = p.get< float    >("MinSig");
    fMinWidth         	   	 = p.get< int >("MinWidth");
    fMinADCSum		   	 = p.get< double >("MinADCSum");
    fMinADCSumOverWidth    	 = p.get< double >("MinADCSumOverWidth");
    fChargeNorm       	   	 = p.get< double >("ChargeNorm");
    fTicksToStopPeakFinder 	 = p.get< double >("TicksToStopPeakFinder");
    fMinTau          	   	 = p.get< double >("MinTau");
    fMaxTau           	   	 = p.get< double >("MaxTau");
    fFitPeakMeanRange	   	 = p.get< double >("FitPeakMeanRange");
    fGroupMaxDistance      	 = p.get< int >("GroupMaxDistance");
    fGroupMinADC           	 = p.get< int >("GroupMinADC");
    fDoMergePeaks	   	 = p.get< bool   >("DoMergePeaks");
    fMergeADCSumThreshold  	 = p.get< double >("MergeADCSumThreshold");
    fMergeMaxADCThreshold  	 = p.get< double >("MergeMaxADCThreshold");
    fWidthNormalization    	 = p.get< double >("WidthNormalization");
    fLongMaxHits           	 = p.get< double >("LongMaxHits");
    fLongPulseWidth        	 = p.get< double >("LongPulseWidth");
}  

//-------------------------------------------------
//-------------------------------------------------
void DPRawHitFinder::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
   
    // ======================================
    // === Hit Information for Histograms ===
    fFirstChi2 = tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
    fChi2      = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
}

//-------------------------------------------------
//-------------------------------------------------
void DPRawHitFinder::endJob()
{

}

//-------------------------------------------------
void DPRawHitFinder::produce(art::Event& evt)
{
  //==================================================================================================
  TH1::AddDirectory(kFALSE);
   
  //Instantiate and Reset a stop watch
  //TStopwatch StopWatch;
  //StopWatch.Reset();
   
  // ################################
  // ### Calling Geometry service ###
  // ################################
  art::ServiceHandle<geo::Geometry> geom;

  // ###############################################
  // ### Making a ptr vector to put on the event ###
  // ###############################################
  // this contains the hit collection
  // and its associations to wires and raw digits
  recob::HitCollectionCreator hcol(*this, evt);
    
  // start collection of fit parameters, initialize metadata describing it
  auto hitID = fHitParamWriter.initOutputs<recob::Hit>(fNewHitsTag, { "t0", "tau1", "tau2" });
   
  // ##########################################
  // ### Reading in the Wire List object(s) ###
  // ##########################################
  art::Handle< std::vector<recob::Wire> > wireVecHandle;
  evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
   
  // #################################################################
  // ### Reading in the RawDigit associated with these wires, too  ###
  // #################################################################
  art::FindOneP<raw::RawDigit> RawDigits(wireVecHandle, evt, fCalDataModuleLabel);
  // Channel Number
  raw::ChannelID_t channel = raw::InvalidChannelID;
    
  //##############################
  //### Looping over the wires ###
  //##############################
  for(size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++)
  {
    // ####################################
    // ### Getting this particular wire ###
    // ####################################
    art::Ptr<recob::Wire>   wire(wireVecHandle, wireIter);
    art::Ptr<raw::RawDigit> rawdigits = RawDigits.at(wireIter);
    // --- Setting Channel Number and Signal type ---
    channel = wire->Channel();        
    // get the WireID for this hit
    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
    // for now, just take the first option returned from ChannelToWire
    geo::WireID wid  = wids[0];

      if(fLogLevel >= 1)
      {
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "-----------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "Channel: " << channel << std::endl;
	std::cout << "Cryostat: " << wid.Cryostat << std::endl;
	std::cout << "TPC: " << wid.TPC << std::endl;
	std::cout << "Plane: " << wid.Plane << std::endl;
	std::cout << "Wire: " << wid.Wire << std::endl;
      }


      // #################################################
      // ### Set up to loop over ROI's for this wire   ###
      // #################################################
      const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();
       
      int CountROI=0;
      for(const auto& range : signalROI.get_ranges())
      {
        // #################################################
        // ### Getting a vector of signals for this wire ###
        // #################################################
        const std::vector<float>& signal = range.data();

        // ##########################################################
        // ### Making an iterator for the time ticks of this wire ###
        // ##########################################################
        std::vector<float>::const_iterator timeIter;  	    // iterator for time bins
           
        // ROI start time
        raw::TDCtick_t roiFirstBinTick = range.begin_index();
        MergedTimeWidVec mergedVec;

        // ###########################################################
        // ### If option set do bin averaging before finding peaks ###
        // ###########################################################
            
        if (fNumBinsToAverage > 1)
        {
          std::vector<float> timeAve;
          doBinAverage(signal, timeAve, fNumBinsToAverage);
            
          // ###################################################################
          // ### Search current averaged ROI for candidate peaks and widths  ###
          // ###################################################################
          TimeValsVec timeValsVec;
          findCandidatePeaks(timeAve.begin(),timeAve.end(),timeValsVec,fMinSig,0);
                
          // ####################################################
          // ### If no startTime hit was found skip this wire ###
          // ####################################################
          if (timeValsVec.empty()) continue;
                
          // #############################################################
          // ### Merge potentially overlapping peaks and do multi fit  ###
          // #############################################################
          mergeCandidatePeaks(timeAve, timeValsVec, mergedVec);
        }
            
        // ###########################################################
        // ### Otherwise, operate directonly on signal vector      ###
        // ###########################################################
        else
        {
          // ##########################################################
          // ### Search current ROI for candidate peaks and widths  ###
          // ##########################################################
          TimeValsVec timeValsVec;
          findCandidatePeaks(signal.begin(),signal.end(),timeValsVec,fMinSig,0);

	  if(fLogLevel >=1)
	  {
	    std::cout << std::endl;
	    std::cout << std::endl;
	    std::cout << "-------------------- ROI #" << CountROI << " -------------------- " << std::endl;
	    if(timeValsVec.size() == 1) std::cout << "ROI #" << CountROI << " (" << timeValsVec.size() << " peak):   ROIStartTick: " << range.offset << "    ROIEndTick:" << range.offset+range.size() << std::endl;
	    else std::cout << "ROI #" << CountROI << " (" << timeValsVec.size() << " peaks):   ROIStartTick: " << range.offset << "    ROIEndTick:" << range.offset+range.size() << std::endl;
            CountROI++;
	  }

 	  if(fLogLevel >=2)
	  {
	    int CountPeak=0;
            for( auto const& timeValsTmp : timeValsVec )
	    {
	      std::cout << "Peak #" << CountPeak << ":   PeakStartTick: " << range.offset + std::get<0>(timeValsTmp) << "    PeakMaxTick: " << range.offset + std::get<1>(timeValsTmp) << "    PeakEndTick: " << range.offset + std::get<2>(timeValsTmp) << std::endl;
	      CountPeak++; 
	    }
	  }
          // ####################################################
          // ### If no startTime hit was found skip this wire ###
          // ####################################################
          if (timeValsVec.empty()) continue;
            
          // #############################################################
          // ### Merge potentially overlapping peaks and do multi fit  ###
          // #############################################################
          mergeCandidatePeaks(signal, timeValsVec, mergedVec);

        }
            
        // #######################################################
        // ### Creating the parameter vector for the new pulse ###
        // #######################################################
        ParameterVec paramVec;

        // === Number of Exponentials to try ===
	int NumberOfPeaksBeforeFit=0;
        unsigned int nExponentialsForFit=0;
        double       chi2PerNDF=0.;
        int          NDF=0;

	unsigned int NumberOfMergedVecs = mergedVec.size();

        // ################################################################
        // ### Lets loop over the groups of peaks we found on this wire ###
        // ################################################################

        for(unsigned int j=0; j < NumberOfMergedVecs; j++)
        {
          int startT = std::get<0>(mergedVec.at(j));
          int endT   = std::get<1>(mergedVec.at(j));
	  int width  = endT + 1 - startT;
          PeakTimeWidVec& peakVals = std::get<2>(mergedVec.at(j));

 	  if(fLogLevel >=3)
	  {
	    std::cout << std::endl;
	    if(peakVals.size() == 1) std::cout << "- Group #" << j << " (" << peakVals.size() << " peak):  GroupStartTick: " << range.offset + startT << "    GroupEndTick: " << range.offset + endT << std::endl;
	    else std::cout << "- Group #" << j << " (" << peakVals.size() << " peaks):  GroupStartTick: " << range.offset + startT << "    GroupEndTick: " << range.offset + endT << std::endl;

	    int CountPeakInGroup=0;
	    for( auto const& peakValsTmp : peakVals )
	    {
	      std::cout << "Peak #" << CountPeakInGroup << " in group #" << j << ":  PeakInGroupStartTick: " << range.offset + std::get<2>(peakValsTmp) << "    PeakInGroupMaxTick: " <<  range.offset + std::get<0>(peakValsTmp) << "    PeakInGroupEndTick: " << range.offset + std::get<3>(peakValsTmp) << std::endl;
	      CountPeakInGroup++;
	    }
	  }

          // ### Getting rid of noise hits ###
          if (width < fMinWidth || (double)std::accumulate(signal.begin()+startT, signal.begin()+endT+1, 0) < fMinADCSum || (double)std::accumulate(signal.begin()+startT, signal.begin()+endT+1, 0)/width < fMinADCSumOverWidth)
	  {
	    if(fLogLevel >=3)
	    {
	      std::cout << "Delete this group of peaks because width, integral or width/intergral is too small." << std::endl;
	    }
	    continue;
	  }


          // #####################################################################################################
          // ### Only attempt to fit if number of peaks <= fMaxMultiHit and if group length <= fMaxGroupLength ###
          // #####################################################################################################
	  NumberOfPeaksBeforeFit = peakVals.size();
	  nExponentialsForFit = peakVals.size();
	  chi2PerNDF = 0.;
	  NDF = 0;
	  if(NumberOfPeaksBeforeFit <= fMaxMultiHit && width <= fMaxGroupLength)
	  {
	    // #####################################################
            // ### Calling the function for fitting Exponentials ###
            // #####################################################
	    paramVec.clear();
	    FitExponentials(signal, peakVals, startT, endT, paramVec, chi2PerNDF, NDF);

	    if(fLogLevel >=4)
	    {
	      std::cout << std::endl;
	      std::cout << "--- First fit ---" << std::endl;
	      if (nExponentialsForFit == 1) std::cout << "- Fitted " << nExponentialsForFit << " peak in group #"  << j << ":" << std::endl;
	      else std::cout << "- Fitted " << nExponentialsForFit << " peaks in group #"  << j << ":" << std::endl;
	      std::cout << "chi2/ndf = " << std::setprecision(2) << std::fixed << chi2PerNDF << std::endl;
	      std::cout << "tau1 [mus] = " << std::setprecision(3) << std::fixed << paramVec[0].first << std::endl;
	      std::cout << "tau2 [mus] = " << std::setprecision(3) << std::fixed << paramVec[1].first << std::endl;

	      for(unsigned int i = 0; i < nExponentialsForFit; i++)
	      {
		std::cout << "Peak #" << i << ": A [ADC] = " << std::setprecision(1) << std::fixed << paramVec[2*(i+1)].first << std::endl;
		std::cout << "Peak #" << i << ": t0 [ticks] = " << std::setprecision(1) << std::fixed << range.offset + paramVec[2*(i+1)+1].first << std::endl;
	      }
	    }

	    // If the chi2 is infinite then there is a real problem so we bail
	    if (!(chi2PerNDF < std::numeric_limits<double>::infinity())) continue;
                   
	    fFirstChi2->Fill(chi2PerNDF);
                
	    // ########################################################
	    // ### Trying extra Exponentials for an initial bad fit ###
	    // ########################################################

	    if( (fTryNplus1Fits && nExponentialsForFit == 1 && chi2PerNDF > fChi2NDFRetry) ||
	        (fTryNplus1Fits && nExponentialsForFit > 1 && chi2PerNDF > fChi2NDFRetryFactorMultiHits*fChi2NDFRetry) )
	    {
	      unsigned int nExponentialsBeforeRefit=nExponentialsForFit;  
	      unsigned int nExponentialsAfterRefit=nExponentialsForFit; 
	      double oldChi2PerNDF = chi2PerNDF;
	      double chi2PerNDF2;
	      int    NDF2;
	      bool   RefitSuccess;
	      PeakTimeWidVec peakValsTemp;
	      while( (nExponentialsForFit == 1 && nExponentialsAfterRefit < 3*nExponentialsBeforeRefit && chi2PerNDF > fChi2NDFRetry) ||
		     (nExponentialsForFit > 1 && nExponentialsAfterRefit < 3*nExponentialsBeforeRefit && chi2PerNDF > fChi2NDFRetryFactorMultiHits*fChi2NDFRetry) )
	      {
		RefitSuccess = false;
		PeakDevVec PeakDev;
	 	FindPeakWithMaxDeviation(signal, nExponentialsForFit, startT, endT, paramVec, peakVals, PeakDev);

		//Add peak and re-fit
		for(auto& PeakDevCand : PeakDev)
		{
		  chi2PerNDF2=0.;
		  NDF2=0.;
		  ParameterVec paramVecRefit;
		  peakValsTemp = peakVals;

		  AddPeak(PeakDevCand, peakValsTemp);
		  FitExponentials(signal, peakValsTemp, startT, endT, paramVecRefit, chi2PerNDF2, NDF2);

		  if (chi2PerNDF2 < chi2PerNDF)
		  {
		    paramVec 	    = paramVecRefit;
		    peakVals	    = peakValsTemp;
		    nExponentialsForFit = peakVals.size();
		    chi2PerNDF  	    = chi2PerNDF2;
		    NDF         	    = NDF2;
		    nExponentialsAfterRefit++;
		    RefitSuccess = true;
		    break;
		  }
		}
			
		//Split peak and re-fit
		if(RefitSuccess == false)
		{
		  for(auto& PeakDevCand : PeakDev)
		  {
		    chi2PerNDF2=0.;
		    NDF2=0.;
		    ParameterVec paramVecRefit;
		    peakValsTemp=peakVals;

		    SplitPeak(PeakDevCand, peakValsTemp);
		    FitExponentials(signal, peakValsTemp, startT, endT, paramVecRefit, chi2PerNDF2, NDF2);

		    if (chi2PerNDF2 < chi2PerNDF)
		    {
		      paramVec 	        = paramVecRefit;
		      peakVals	        = peakValsTemp;
		      nExponentialsForFit = peakVals.size();
		      chi2PerNDF  	= chi2PerNDF2;
		      NDF         	= NDF2;
		      nExponentialsAfterRefit++;
		      RefitSuccess = true;
		      break;
		    }
		  }
		}

		if(RefitSuccess == false)
		{
		  break;
		}	
	      }

	      if(fLogLevel >=5)
	      {
		std::cout << std::endl;
		std::cout << "--- Refit ---" << std::endl;
		if( chi2PerNDF == oldChi2PerNDF) std::cout << "chi2/ndf didn't improve. Keep first fit." << std::endl;
		else
		{
		  std::cout << "- Added peaks to group #" << j << ". This group now has " << nExponentialsForFit << " peaks:" << std::endl;
		  std::cout << "- Group #" << j << " (" << peakVals.size() << " peaks):  GroupStartTick: " << range.offset + startT << "    GroupEndTick: " << range.offset + endT << std::endl;

		  int CountPeakInGroup=0;
		  for( auto const& peakValsTmp : peakVals )
		  {
		    std::cout << "Peak #" << CountPeakInGroup << " in group #" << j << ":  PeakInGroupStartTick: " << range.offset + std::get<2>(peakValsTmp) << "    PeakInGroupMaxTick: " <<  range.offset + std::get<0>(peakValsTmp) << "    PeakInGroupEndTick: " << range.offset + std::get<3>(peakValsTmp) << std::endl;
		    CountPeakInGroup++;
		  }

		  std::cout << "chi2/ndf = " << std::setprecision(2) << std::fixed << chi2PerNDF << std::endl;
		  std::cout << "tau1 [mus] = " << std::setprecision(3) << std::fixed << paramVec[0].first << std::endl;
		  std::cout << "tau2 [mus] = " << std::setprecision(3) << std::fixed << paramVec[1].first << std::endl;

		  for(unsigned int i = 0; i < nExponentialsForFit; i++)
		  {
		    std::cout << "Peak #" << i << ": A [ADC] = " << std::setprecision(1) << std::fixed << paramVec[2*(i+1)].first << std::endl;
		    std::cout << "Peak #" << i << ": t0 [ticks] = " << std::setprecision(1) << std::fixed << range.offset + paramVec[2*(i+1)+1].first << std::endl;
		  }
		}
	      }
	    }

            // #######################################################
            // ### Loop through returned peaks and make recob hits ###
            // #######################################################
                
            int numHits(0);
 	    //Check chi2PerNDF
	    if( ( nExponentialsForFit == 1 && chi2PerNDF <= fChi2NDFMax ) || ( nExponentialsForFit >= 2 && chi2PerNDF <= fChi2NDFMaxFactorMultiHits*fChi2NDFMax ) )
	    {
              for(unsigned int i = 0; i < nExponentialsForFit; i++)
              {
                //Extract fit parameters for this hit
		double peakTau1 = paramVec[0].first;
	        double peakTau2 = paramVec[1].first;
                double peakAmp   = paramVec[2*(i+1)].first;
                double peakMean  = paramVec[2*(i+1)+1].first;

	 	//Calculate mean
		TF1 Exponentials("Exponentials","( [0] * exp(0.4*(x-[1])/[2]) / ( 1 + exp(0.4*(x-[1])/[3]) ) )",startT,endT);
        	Exponentials.SetParameter(0, peakAmp);
        	Exponentials.SetParameter(1, peakMean);
        	Exponentials.SetParameter(2, peakTau1);
        	Exponentials.SetParameter(3, peakTau2);
		double peakMeanTrue = Exponentials.GetMaximumX(startT,endT);
		Exponentials.Delete();

		//Calculate width (=FWHM)
		double peakWidth = WidthFunc(peakMean, peakAmp, peakTau1, peakTau2, startT, endT, peakMeanTrue);
		peakWidth /= fWidthNormalization; //from FWHM to "standard deviation": standard deviation = FWHM/(2*sqrt(2*ln(2)))

                // Extract fit parameters errors
                double peakAmpErr   = paramVec[2*(i+1)].second;
                double peakMeanErr  = paramVec[2*(i+1)+1].second;
                double peakWidthErr = 0.1*peakWidth;

                // ### Charge ###
                double charge = ChargeFunc(peakMean, peakAmp, peakTau1, peakTau2, fChargeNorm, peakMeanTrue);
                double chargeErr = std::sqrt(TMath::Pi()) * (peakAmpErr*peakWidthErr + peakWidthErr*peakAmpErr);    

                // ### limits for getting sum of ADC counts
	        int startTthisHit = std::get<2>(peakVals.at(i));
	        int endTthisHit = std::get<3>(peakVals.at(i));
                std::vector<float>::const_iterator sumStartItr = signal.begin() + startTthisHit;
                std::vector<float>::const_iterator sumEndItr   = signal.begin() + endTthisHit;

                // ### Sum of ADC counts
                double sumADC = std::accumulate(sumStartItr, sumEndItr, 0.);


		//Check if fit returns reasonable values
		if(peakWidth <= 0 || charge <= 0. || charge != charge)
		{
		  if(fLogLevel >= 1)
		  {
		    std::cout << std::endl;
		    std::cout << "WARNING: For peak #" << i << " in this group:" << std::endl;
		    std::cout << "Fit function returned width < 0 or charge < 0 or charge = nan." << std::endl;
		    std::cout << "---> DO NOT create hit object from fit parameters but use peak values instead." << std::endl;
		    std::cout << "---> Set fit parameter so that a sharp peak with a width of 1 tick is shown in the event display. This indicates that the fit failed." << std::endl;
		  }
		  peakWidth = ( ( (double)endTthisHit - (double)startTthisHit )/4. ) / fWidthNormalization; //~4 is the factor between FWHM and full width of the hit (last bin - first bin). no drift: 4.4, 6m drift: 3.7
                  peakAmp   = 0.3989 * sumADC / peakWidth;  // Use gaussian formulation
		  peakAmpErr = 0.1*peakAmp;
		  peakMeanErr=peakWidth/2;
		  charge = sumADC;
		  peakMeanTrue = std::get<0>(peakVals.at(i));
		  peakWidth *= 2;	//double the peak width again (overestimating the width is safer than underestimating it)

		  //set the fit values to make it visible in the event display that this fit failed
                  peakMean = peakMeanTrue-2;
                  peakTau1 = 0.008;
                  peakTau2 = 0.0065;
		}

                // Create the hit
		recob::HitCreator hitcreator(*wire,                            // wire reference
                                             wid,                              // wire ID
                                             startT+roiFirstBinTick,           // start_tick TODO check
                                             endT+roiFirstBinTick,             // end_tick TODO check
                                             peakWidth,                        // rms
                                             peakMeanTrue+roiFirstBinTick,     // peak_time
                                             peakMeanErr,                      // sigma_peak_time
                                             peakAmp,                          // peak_amplitude
                                             peakAmpErr,                       // sigma_peak_amplitude
                                             charge,                           // hit_integral
                                             chargeErr,                        // hit_sigma_integral
                                             sumADC,                           // summedADC FIXME
                                             nExponentialsForFit,              // multiplicity
                                             numHits,                          // local_index TODO check that the order is correct
                                             chi2PerNDF,                       // goodness_of_fit
                                             NDF                               // dof
                                             );

		if(fLogLevel >=6)
	    	{
		  std::cout << std::endl;
		  std::cout << "- Created hit object for peak #" << i << " in this group with the following parameters (obtained from fit):" << std::endl;
		  std::cout << "HitStartTick: " << startT+roiFirstBinTick << std::endl;
		  std::cout << "HitEndTick: " << endT+roiFirstBinTick << std::endl;
		  std::cout << "HitWidthTicks: " << std::setprecision(2) << std::fixed << peakWidth << std::endl;
		  std::cout << "HitMeanTick: " << std::setprecision(2) << std::fixed << peakMeanTrue+roiFirstBinTick << " +- " << peakMeanErr << std::endl;
		  std::cout << "HitAmplitude [ADC]: " << std::setprecision(1) << std::fixed << peakAmp << " +- " << peakAmpErr << std::endl;
		  std::cout << "HitIntegral [ADC*ticks]: " << std::setprecision(1) << std::fixed << charge << " +- " << chargeErr << std::endl;
		  std::cout << "HitADCSum [ADC*ticks]: " << std::setprecision(1) << std::fixed << sumADC << std::endl;
		  std::cout << "HitMultiplicity: " << nExponentialsForFit << std::endl;
		  std::cout << "HitIndex in group: " << numHits << std::endl;
		  std::cout << "Hitchi2/ndf: " << std::setprecision(2) << std::fixed << chi2PerNDF << std::endl;
		  std::cout << "HitNDF: " << NDF << std::endl;
		}

                const recob::Hit hit(hitcreator.move());
		    
                hcol.emplace_back(std::move(hit), wire, rawdigits);
                // add fit parameters associated to the hit just pushed to the collection
                std::array<float, 3> fitParams;
                fitParams[0] = peakMean+roiFirstBinTick;
                fitParams[1] = peakTau1;
                fitParams[2] = peakTau2;
                fHitParamWriter.addVector(hitID, fitParams);
                numHits++;
              } // <---End loop over Exponentials
            } // <---End if chi2 <= chi2Max
	  } // <---End if(NumberOfPeaksBeforeFit <= fMaxMultiHit && width <= fMaxGroupLength), then fit

          // #######################################################
          // ### If too large then force alternate solution      ###
          // ### - Make n hits from pulse train where n will     ###
          // ###   depend on the fhicl parameter fLongPulseWidth ###
          // ### Also do this if chi^2 is too large              ###
          // #######################################################
          if( NumberOfPeaksBeforeFit > fMaxMultiHit || (width > fMaxGroupLength) || ( nExponentialsForFit == 1 && chi2PerNDF > fChi2NDFMax ) || ( nExponentialsForFit >= 2 && chi2PerNDF > fChi2NDFMaxFactorMultiHits*fChi2NDFMax ) )
          {

            int nHitsInThisGroup = (endT - startT) / fLongPulseWidth;
                    
            if (nHitsInThisGroup > fLongMaxHits)
            {
              nHitsInThisGroup = fLongMaxHits;
              fLongPulseWidth = (endT - startT) / nHitsInThisGroup;
            }
                    
            if (nHitsInThisGroup * fLongPulseWidth < endT - startT) nHitsInThisGroup++;
                    
            int firstTick = startT;
            int lastTick  = firstTick + std::min(endT,fLongPulseWidth) -1;

	    if(fLogLevel >= 1)
	    {
	      if( NumberOfPeaksBeforeFit > fMaxMultiHit)
	      {
		std::cout << std::endl;
		std::cout << "WARNING: Number of peaks in this group (" << NumberOfPeaksBeforeFit << ") is higher than threshold (" <<  fMaxMultiHit << ")." << std::endl;
		std::cout << "---> DO NOT fit. Split group of peaks into hits with equal length instead." << std::endl;
	      }
	      if( width > fMaxGroupLength)
	      {
		std::cout << std::endl;
		std::cout << "WARNING: group of peak is longer (" << width << " ticks) than threshold (" <<  fMaxGroupLength << " ticks)." << std::endl;
		std::cout << "---> DO NOT fit. Split group of peaks into hits with equal length instead." << std::endl;
	      }

	      if( ( nExponentialsForFit == 1 && chi2PerNDF > fChi2NDFMax ) || ( nExponentialsForFit >= 2 && chi2PerNDF > fChi2NDFMaxFactorMultiHits*fChi2NDFMax ) )
	      {
		std::cout << std::endl;
	      	std::cout << "WARNING: For fit of this group (" <<  NumberOfPeaksBeforeFit << " peaks before refit, " << nExponentialsForFit << " peaks after refit): " << std::endl;
	      	if ( nExponentialsForFit == 1 && chi2PerNDF > fChi2NDFMax ) std::cout << "chi2/ndf of this fit (" << chi2PerNDF << ") is higher than threshold (" << fChi2NDFMax << ")." << std::endl;
	      	if ( nExponentialsForFit >= 2 && chi2PerNDF > fChi2NDFMaxFactorMultiHits*fChi2NDFMax ) std::cout << "chi2/ndf of this fit (" << chi2PerNDF << ") is higher than threshold (" << fChi2NDFMaxFactorMultiHits*fChi2NDFMax << ")." << std::endl;
	        std::cout << "---> DO NOT create hit object but split group of peaks into hits with equal length instead." << std::endl;
	      }
	      std::cout << "---> Group goes from tick " << roiFirstBinTick+startT << " to " << roiFirstBinTick+endT << ". Split group into (" << roiFirstBinTick+endT << " - " << roiFirstBinTick+startT << ")/" << fLongPulseWidth << " = " <<  (endT - startT) << "/" << fLongPulseWidth << " = " << nHitsInThisGroup << " peaks (" << fLongPulseWidth << " = LongPulseWidth), or maximum LongMaxHits = " << fLongMaxHits << " peaks." << std::endl;
	    }

                   
            for(int hitIdx = 0; hitIdx < nHitsInThisGroup; hitIdx++)
            {
              // This hit parameters
              double peakWidth = (lastTick - firstTick) / 3.; 
              double peakMeanTrue  = (firstTick + lastTick) / 2.;
	      double peakMeanErr = (lastTick - firstTick) / 2.;
              double sumADC    = std::accumulate(signal.begin() + firstTick, signal.begin() + lastTick, 0.);
	      double charge = sumADC;
	      double chargeErr = 0.1*sumADC;
              double peakAmp   = 0.3989 * sumADC / peakWidth;  // Use gaussian formulation
	      double peakAmpErr = 0.1*peakAmp;
	      nExponentialsForFit = nHitsInThisGroup;
              NDF         = 1;
              chi2PerNDF  =  chi2PerNDF > fChi2NDFRetryFactorMultiHits*fChi2NDFRetry ? chi2PerNDF : -1.;
              //chi2PerNDF  = -1.;

	      //set the fit values to make it visible in the event display that this fit failed
              double peakMean = peakMeanTrue-2;
              double peakTau1 = 0.008;
              double peakTau2 = 0.0065;
                        
              recob::HitCreator hitcreator(*wire,                            // wire reference
                                           wid,                              // wire ID
                                           firstTick+roiFirstBinTick,        // start_tick TODO check
                                           lastTick+roiFirstBinTick,         // end_tick TODO check
                                           peakWidth,                        // rms
                                           peakMeanTrue+roiFirstBinTick,     // peak_time
                                           peakMeanErr,                      // sigma_peak_time
                                           peakAmp,                          // peak_amplitude
                                           peakAmpErr,                       // sigma_peak_amplitude
                                           charge,                           // hit_integral
                                           chargeErr,                        // hit_sigma_integral
                                           sumADC,                           // summedADC FIXME
                                           nExponentialsForFit,              // multiplicity
                                           hitIdx,                          // local_index TODO check that the order is correct
                                           chi2PerNDF,                       // goodness_of_fit
                                           NDF                               // dof
                                           );


	      if(fLogLevel >=6)
	      {
	        std::cout << std::endl;
	        std::cout << "- Created hit object for peak #" << hitIdx << " in this group with the following parameters (obtained from waveform):" << std::endl;
	        std::cout << "HitStartTick: " << firstTick+roiFirstBinTick << std::endl;
	        std::cout << "HitEndTick: " << lastTick+roiFirstBinTick << std::endl;
	        std::cout << "HitWidthTicks: " << std::setprecision(2) << std::fixed << peakWidth << std::endl;
	        std::cout << "HitMeanTick: " << std::setprecision(2) << std::fixed << peakMeanTrue+roiFirstBinTick << " +- " << peakMeanErr << std::endl;
	        std::cout << "HitAmplitude [ADC]: " << std::setprecision(1) << std::fixed << peakAmp << " +- " << peakAmpErr << std::endl;
	        std::cout << "HitIntegral [ADC*ticks]: " << std::setprecision(1) << std::fixed << charge << " +- " << chargeErr << std::endl;
	        std::cout << "HitADCSum [ADC*ticks]: " << std::setprecision(1) << std::fixed << sumADC << std::endl;
	        std::cout << "HitMultiplicity: " << nExponentialsForFit << std::endl;
	        std::cout << "HitIndex in group: " << hitIdx << std::endl;
	        std::cout << "Hitchi2/ndf: " << std::setprecision(2) << std::fixed << chi2PerNDF << std::endl;
	        std::cout << "HitNDF: " << NDF << std::endl;
	      }   
              const recob::Hit hit(hitcreator.move());
              hcol.emplace_back(std::move(hit), wire, rawdigits);

              std::array<float, 3> fitParams;
              fitParams[0] = peakMean+roiFirstBinTick;
              fitParams[1] = peakTau1;
              fitParams[2] = peakTau2;
              fHitParamWriter.addVector(hitID, fitParams);
    
              // set for next loop
              firstTick = lastTick+1;
              lastTick  = std::min(lastTick + fLongPulseWidth, endT);
            }//<---Hits in this group
	  }//<---End if #peaks > MaxMultiHit
          fChi2->Fill(chi2PerNDF);
         }//<---End loop over merged candidate hits
       } //<---End looping over ROI's
     }//<---End looping over all the wires

    //==================================================================================================
    // End of the event
   
    // move the hit collection and the associations into the event
    hcol.put_into(evt);

    // and put hit fit parameters together with metadata into the event
    fHitParamWriter.saveOutputs(evt);

} // End of produce() 
    
// --------------------------------------------------------------------------------------------
// Initial finding of candidate peaks
// --------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::findCandidatePeaks(std::vector<float>::const_iterator   startItr,
                                            std::vector<float>::const_iterator    stopItr,
                                            std::vector<std::tuple<int,int,int>>& timeValsVec,
                                            float&                                PeakMin,
                                            int                                   firstTick) const
{
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        auto maxItr = std::max_element(startItr, stopItr);
        
        float maxValue = *maxItr;
        int   maxTime  = std::distance(startItr,maxItr);
        
        if (maxValue >= PeakMin)
        {
            // backwards to find first bin for this candidate hit
            auto firstItr = std::distance(startItr,maxItr) > 2 ? maxItr - 1 : startItr;
	    bool PeakStartIsHere = true;

            while(firstItr != startItr)
            {    
                //Check for inflection point & ADC count <= 0
		PeakStartIsHere = true;
		for(int i=1; i <= fTicksToStopPeakFinder; i++)
		{
		    if( *firstItr >= *(firstItr-i) )
		    {
		    PeakStartIsHere = false;
		    break;
		    }
	
		}
                if (*firstItr <= 0 || PeakStartIsHere) break;
                firstItr--;
            }

            int firstTime = std::distance(startItr,firstItr);
            
            // Recursive call to find all candidate hits earlier than this peak
            findCandidatePeaks(startItr, firstItr - 1, timeValsVec, PeakMin, firstTick);
            
            // forwards to find last bin for this candidate hit
            auto lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;
	    bool PeakEndIsHere = true;

            while(lastItr != stopItr)
            {     
                //Check for inflection point & ADC count <= 0
		PeakEndIsHere = true;
		for(int i=1; i <= fTicksToStopPeakFinder; i++)
		{
		    if( *lastItr >= *(lastItr+i) )
		    {
		    PeakEndIsHere = false;
		    break;
		    }
		}
                if (*lastItr <= 0 || PeakEndIsHere) break;
                lastItr++;
            }

            int lastTime = std::distance(startItr,lastItr);
            
            // Now save this candidate's start and max time info
            timeValsVec.push_back(std::make_tuple(firstTick+firstTime,firstTick+maxTime,firstTick+lastTime));
            
            // Recursive call to find all candidate hits later than this peak
            findCandidatePeaks(lastItr + 1, stopItr, timeValsVec, PeakMin, firstTick + std::distance(startItr,lastItr + 1));
        }
    }
    
    return;
}
    
// --------------------------------------------------------------------------------------------
// Merging of nearby candidate peaks
// --------------------------------------------------------------------------------------------
    
void hit::DPRawHitFinder::mergeCandidatePeaks(const std::vector<float>& signalVec, TimeValsVec& timeValsVec, MergedTimeWidVec& mergedVec) const
{
    // ################################################################
    // ### Lets loop over the candidate pulses we found in this ROI ###
    // ################################################################
    auto timeValsVecItr = timeValsVec.begin();
    unsigned int PeaksInThisMergedPeak = 0;
    //int EndTickOfPreviousMergedPeak=0;
    while(timeValsVecItr != timeValsVec.end())
    {
        PeakTimeWidVec peakVals;
        
        // Setting the start, peak, and end time of the pulse
        auto& timeVal = *timeValsVecItr++;
        int startT = std::get<0>(timeVal);
        int maxT   = std::get<1>(timeVal);
        int endT   = std::get<2>(timeVal);
        int widT   = endT - startT;
	int FinalStartT=startT;
	int FinalEndT=endT;

        peakVals.emplace_back(maxT,widT,startT,endT);
        
        // See if we want to merge pulses together
        // First check if we have more than one pulse on the wire
        bool checkNextHit = timeValsVecItr != timeValsVec.end();

        // Loop until no more merged pulses (or candidates in this ROI)
        while(checkNextHit)
        {
            // group hits if start time of the next pulse is < end time + fGroupMaxDistance of current pulse and if intermediate signal between two pulses doesn't go below fMinBinToGroup
            int NextStartT = std::get<0>(*timeValsVecItr);
	    
	    double MinADC = signalVec[endT];
	    for(int i = endT; i <= NextStartT; i++)
	    {
		if(signalVec[i]<MinADC)
		{
		MinADC = signalVec[i];
		}
	    }
	    
	    // Group peaks (grouped peaks are fitted together and can be merged)
            if( MinADC >= fGroupMinADC && NextStartT - endT <= fGroupMaxDistance )
            {
	    	int PrevStartT=startT;
	    	int PrevMaxT=maxT;
	    	int PrevEndT=endT;
	    	//int PrevWidT=widT;

                timeVal = *timeValsVecItr++;
                int NextMaxT = std::get<1>(timeVal);
                int NextEndT = std::get<2>(timeVal);
                int NextWidT = NextEndT - NextStartT;
	    	FinalEndT=NextEndT;

		int PrevSumADC = 0;
		for(int i = PrevStartT; i<= PrevEndT; i++)
		{
		    PrevSumADC+=signalVec[i];
		} 

		int NextSumADC = 0;
		for (int i = NextStartT; i<= NextEndT; i++)
		{
		    NextSumADC+=signalVec[i];
		} 
		
		//Merge peaks within a group
		if(fDoMergePeaks)
		{
		    if( signalVec[NextMaxT] <= signalVec[PrevMaxT] && signalVec[NextMaxT] < fMergeMaxADCThreshold*signalVec[PrevMaxT] &&  NextSumADC < fMergeADCSumThreshold*PrevSumADC )
		    {
		    	maxT=PrevMaxT;
		    	startT=PrevStartT;
		    	endT=NextEndT;
		    	widT=endT - startT;
		    	peakVals.pop_back();
                    	peakVals.emplace_back(maxT,widT,startT,endT);
		    }
		    else if( signalVec[NextMaxT] > signalVec[PrevMaxT] && signalVec[PrevMaxT] < fMergeMaxADCThreshold*signalVec[NextMaxT] &&  PrevSumADC < fMergeADCSumThreshold*NextSumADC  )
		    {
		    	maxT=NextMaxT;
		    	startT=PrevStartT;
		    	endT=NextEndT;
		    	widT=endT - startT;
		    	peakVals.pop_back();
                    	peakVals.emplace_back(maxT,widT,startT,endT);
		    }
		    else
		    {
		    	maxT=NextMaxT;
		    	startT=NextStartT;
		    	endT=NextEndT;
		    	widT=NextWidT;
                    	peakVals.emplace_back(maxT,widT,startT,endT);
		    	PeaksInThisMergedPeak++;
                    }
		}
		else
		{
		    maxT=NextMaxT;
		    startT=NextStartT;
		    endT=NextEndT;
		    widT=NextWidT;
                    peakVals.emplace_back(maxT,widT,startT,endT);
		    PeaksInThisMergedPeak++;
		}
                checkNextHit = timeValsVecItr != timeValsVec.end();
            }//<---Checking adjacent pulses
            else 
	    {
		checkNextHit = false;
		PeaksInThisMergedPeak = 0;
	    }
            
        }//<---End checking if there is more than one pulse on the wire   
        // Add these to our merged vector
        mergedVec.emplace_back(FinalStartT, FinalEndT, peakVals);
    }
    return;
}

// --------------------------------------------------------------------------------------------
// Fit Exponentials
// --------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::FitExponentials(const std::vector<float>  fSignalVector,
                                          const PeakTimeWidVec      fPeakVals,
                                          int                       fStartTime,
                                          int                       fEndTime,
                                          ParameterVec&             fparamVec,
                                          double&                   fchi2PerNDF,
                                          int&                      fNDF)
{
    int size = fEndTime - fStartTime + 1;
    int NPeaks = fPeakVals.size();

    // #############################################
    // ### If size < 0 then set the size to zero ###
    // #############################################
    if(fEndTime - fStartTime < 0){size = 0;}
   
    // --- TH1D HitSignal ---
    TH1F hitSignal("hitSignal","",std::max(size,1),fStartTime,fEndTime+1);
    hitSignal.Sumw2();
   
    // #############################
    // ### Filling the histogram ###
    // #############################
    for(int i = fStartTime; i < fEndTime+1; i++)
    {
	hitSignal.Fill(i,fSignalVector[i]);
	hitSignal.SetBinError(i,0.288675); //1/sqrt(12)
    }

    // ################################################
    // ### Building TFormula for basic Exponentials ###
    // ################################################
    std::string eqn = CreateFitFunction(NPeaks);

    // --------------------------------------
    // --- TF1 function for Exponentials  ---
    // --------------------------------------

    TF1 Exponentials("Exponentials",eqn.c_str(),fStartTime,fEndTime+1);

    Exponentials.SetParameter(0, 0.5);
    Exponentials.SetParameter(1, 0.5);
    Exponentials.SetParLimits(0, fMinTau, fMaxTau);
    Exponentials.SetParLimits(1, fMinTau, fMaxTau);
    double amplitude=0;
    double peakMean=0;

    double peakMeanShift=2;
    double peakMeanSeed=0;
    double peakMeanRangeLow=0;
    double peakMeanRangeHi=0;
    double peakStart=0;
    double peakEnd=0;

    	for(int i = 0; i < NPeaks; i++)
    	{
        peakMean = std::get<0>(fPeakVals.at(i));
	peakStart = std::get<2>(fPeakVals.at(i));
	peakEnd = std::get<3>(fPeakVals.at(i));
        peakMeanSeed=peakMean-peakMeanShift;
        peakMeanRangeLow = std::max(peakStart-peakMeanShift, peakMeanSeed-fFitPeakMeanRange);
        peakMeanRangeHi = std::min(peakEnd, peakMeanSeed+fFitPeakMeanRange);
        amplitude = fSignalVector[peakMean];

	Exponentials.SetParameter(2*(i+1), 1.65*amplitude);
	Exponentials.SetParLimits(2*(i+1), 0.1*amplitude, 2*1.65*amplitude);
	Exponentials.SetParameter(2*(i+1)+1, peakMeanSeed);

	    if(NPeaks == 1)
	    {
	    Exponentials.SetParLimits(2*(i+1)+1, peakMeanRangeLow, peakMeanRangeHi);
	    }
	    else if(NPeaks >= 2 && i == 0)
	    {
	    double HalfDistanceToNextMean = 0.5*(std::get<0>(fPeakVals.at(i+1)) - peakMean);
	    Exponentials.SetParLimits( 2*(i+1)+1, peakMeanRangeLow, std::min(peakMeanRangeHi, peakMeanSeed+HalfDistanceToNextMean) );
	    }
	    else if(NPeaks >= 2 && i == NPeaks-1)
	    {
	    double HalfDistanceToPrevMean = 0.5*(peakMean - std::get<0>(fPeakVals.at(i-1)));
	    Exponentials.SetParLimits(2*(i+1)+1, std::max(peakMeanRangeLow, peakMeanSeed-HalfDistanceToPrevMean), peakMeanRangeHi );
	    }
	    else
	    {
	    double HalfDistanceToNextMean = 0.5*(std::get<0>(fPeakVals.at(i+1)) - peakMean);
	    double HalfDistanceToPrevMean = 0.5*(peakMean - std::get<0>(fPeakVals.at(i-1)));
	    Exponentials.SetParLimits(2*(i+1)+1, std::max(peakMeanRangeLow, peakMeanSeed-HalfDistanceToPrevMean), std::min(peakMeanRangeHi, peakMeanSeed+HalfDistanceToNextMean) );
	    }	
	}

    // ###########################################
    // ### PERFORMING THE TOTAL FIT OF THE HIT ###
    // ###########################################
    try
      { hitSignal.Fit(&Exponentials,"QNRWM","", fStartTime, fEndTime+1);}
    catch(...)
      {mf::LogWarning("DPRawHitFinder") << "Fitter failed finding a hit";}
   
    // ##################################################
    // ### Getting the fitted parameters from the fit ###
    // ##################################################
    fchi2PerNDF = (Exponentials.GetChisquare() / Exponentials.GetNDF());
    fNDF        = Exponentials.GetNDF();

    fparamVec.emplace_back(Exponentials.GetParameter(0),Exponentials.GetParError(0));
    fparamVec.emplace_back(Exponentials.GetParameter(1),Exponentials.GetParError(1));
 
    for(int i = 0; i < NPeaks; i++)
    {
        fparamVec.emplace_back(Exponentials.GetParameter(2*(i+1)),Exponentials.GetParError(2*(i+1)));
        fparamVec.emplace_back(Exponentials.GetParameter(2*(i+1)+1),Exponentials.GetParError(2*(i+1)+1));
    }
    Exponentials.Delete();
    hitSignal.Delete();
}//<----End FitExponentials


//---------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::FindPeakWithMaxDeviation(const std::vector<float> fSignalVector,
			  	  		   int			    fNPeaks,
                          	  		   int                      fStartTime,
                          	  		   int                      fEndTime,
                          	  		   ParameterVec             fparamVec,
                         	  		   PeakTimeWidVec           fpeakVals,
			  	 		   PeakDevVec& 		    fPeakDev)
{
//   int size = fEndTime - fStartTime + 1;
//    if(fEndTime - fStartTime < 0){size = 0;}

    std::string eqn = CreateFitFunction(fNPeaks);  // string for equation of Exponentials fit

    TF1 Exponentials("Exponentials",eqn.c_str(),fStartTime,fEndTime+1);

 	for(size_t i=0; i < fparamVec.size(); i++)
    	{
    	Exponentials.SetParameter(i, fparamVec[i].first);
	}

    // ##########################################################################
    // ### Finding the peak with the max chi2 fit and signal ###
    // ##########################################################################
    double Chi2PerNDFPeak;
    double MaxPosDeviation;
    double MaxNegDeviation;
    int BinMaxPosDeviation;
    int BinMaxNegDeviation;
    for(int i = 0; i < fNPeaks; i++)
    {
	Chi2PerNDFPeak = 0.;
	MaxPosDeviation=0.;
	MaxNegDeviation=0.;
	BinMaxPosDeviation=0;
	BinMaxNegDeviation=0;

    	for(int j = std::get<2>(fpeakVals.at(i)); j < std::get<3>(fpeakVals.at(i))+1; j++)
	{
            if( (Exponentials(j+0.5)-fSignalVector[j]) > MaxPosDeviation && j != std::get<0>(fpeakVals.at(i)) )
	    {
	    MaxPosDeviation = Exponentials(j+0.5)-fSignalVector[j];
	    BinMaxPosDeviation = j;
	    }
	    if( (Exponentials(j+0.5)-fSignalVector[j]) < MaxNegDeviation && j != std::get<0>(fpeakVals.at(i)) )
	    {
	    MaxNegDeviation = Exponentials(j+0.5)-fSignalVector[j];
	    BinMaxNegDeviation = j;
	    }
	Chi2PerNDFPeak += pow((Exponentials(j+0.5)-fSignalVector[j])/sqrt(fSignalVector[j]),2);
	}

	if(BinMaxNegDeviation != 0)
	{
	Chi2PerNDFPeak /= static_cast<double>((std::get<3>(fpeakVals.at(i))-std::get<2>(fpeakVals.at(i))));
	fPeakDev.emplace_back(Chi2PerNDFPeak,i,BinMaxNegDeviation,BinMaxPosDeviation);
	}
    }

std::sort(fPeakDev.begin(),fPeakDev.end(), [](std::tuple<double,int,int,int> const &t1, std::tuple<double,int,int,int> const &t2) {return std::get<0>(t1) > std::get<0>(t2);} );	
Exponentials.Delete();
}

//---------------------------------------------------------------------------------------------
std::string hit::DPRawHitFinder::CreateFitFunction(int fNPeaks)
{
    std::string feqn = "";  // string for equation for Exponentials fit
    std::stringstream numConv;
    
    for(int i = 0; i < fNPeaks; i++)
    {
        feqn.append("+( [");
        numConv.str("");
        numConv << 2*(i+1);
        feqn.append(numConv.str());
        feqn.append("] * exp(0.4*(x-[");
        numConv.str("");
        numConv << 2*(i+1)+1;
        feqn.append(numConv.str());
        feqn.append("])/[");
        numConv.str("");
        numConv << 0;
        feqn.append(numConv.str());    
        feqn.append("]) / ( 1 + exp(0.4*(x-[");
        numConv.str("");
        numConv << 2*(i+1)+1;
        feqn.append(numConv.str()); 
        feqn.append("])/[");       
        numConv.str("");
        numConv << 1;
        feqn.append(numConv.str()); 
    feqn.append("]) ) )");
    }

return feqn;
}


//---------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::AddPeak(std::tuple<double,int,int,int> fPeakDevCand,
				  PeakTimeWidVec& fpeakValsTemp)
{
int PeakNumberWithNewPeak = std::get<1>(fPeakDevCand);
int NewPeakMax = std::get<2>(fPeakDevCand);
int OldPeakMax = std::get<0>(fpeakValsTemp.at(PeakNumberWithNewPeak));
int OldPeakOldStart = std::get<2>(fpeakValsTemp.at(PeakNumberWithNewPeak));
int OldPeakOldEnd = std::get<3>(fpeakValsTemp.at(PeakNumberWithNewPeak));
			    
int NewPeakStart=0;
int NewPeakEnd=0;
int OldPeakNewStart=0;
int OldPeakNewEnd=0;
int DistanceBwOldAndNewPeak=0;

    if(NewPeakMax<OldPeakMax)
    {
    NewPeakStart = OldPeakOldStart;
    OldPeakNewEnd = OldPeakOldEnd;
    DistanceBwOldAndNewPeak = OldPeakMax - NewPeakMax;
    NewPeakEnd = NewPeakMax+0.5*(DistanceBwOldAndNewPeak-(DistanceBwOldAndNewPeak%2));
    if(DistanceBwOldAndNewPeak%2 == 0) NewPeakEnd -= 1; 
    OldPeakNewStart = NewPeakEnd+1;
    }
    else if(OldPeakMax<NewPeakMax)
    {
    NewPeakEnd = OldPeakOldEnd;
    OldPeakNewStart = OldPeakOldStart;
    DistanceBwOldAndNewPeak = NewPeakMax - OldPeakMax;
    OldPeakNewEnd = OldPeakMax+0.5*(DistanceBwOldAndNewPeak-(DistanceBwOldAndNewPeak%2));
    if(DistanceBwOldAndNewPeak%2 == 0) OldPeakNewEnd -= 1;
    NewPeakStart = OldPeakNewEnd+1;
    }
    else if(OldPeakMax==NewPeakMax){return;} //This shouldn't happen

fpeakValsTemp.at(PeakNumberWithNewPeak) = std::make_tuple(OldPeakMax,0,OldPeakNewStart,OldPeakNewEnd);
fpeakValsTemp.emplace_back(NewPeakMax,0,NewPeakStart,NewPeakEnd);
std::sort(fpeakValsTemp.begin(),fpeakValsTemp.end(), [](std::tuple<int,int,int,int> const &t1, std::tuple<int,int,int,int> const &t2) {return std::get<0>(t1) < std::get<0>(t2);} );

return;
}


//---------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::SplitPeak(std::tuple<double,int,int,int> fPeakDevCand,
				    PeakTimeWidVec& fpeakValsTemp)
{
int PeakNumberWithNewPeak = std::get<1>(fPeakDevCand);
int OldPeakOldStart = std::get<2>(fpeakValsTemp.at(PeakNumberWithNewPeak));
int OldPeakOldEnd = std::get<3>(fpeakValsTemp.at(PeakNumberWithNewPeak));
int WidthOldPeakOld = OldPeakOldEnd - OldPeakOldStart;

if(WidthOldPeakOld<3) {return;} //can't split peaks with widths < 3

int NewPeakMax = 0;
int NewPeakStart = 0;
int NewPeakEnd = 0;
int OldPeakNewMax = 0; 
int OldPeakNewStart = 0;
int OldPeakNewEnd = 0;


OldPeakNewStart = OldPeakOldStart;
NewPeakEnd = OldPeakOldEnd;

OldPeakNewEnd = OldPeakNewStart + 0.5*(WidthOldPeakOld + (WidthOldPeakOld%2));
NewPeakStart = OldPeakNewEnd+1;

int WidthOldPeakNew = OldPeakNewEnd-OldPeakNewStart;
int WidthNewPeak = NewPeakEnd-NewPeakStart;

OldPeakNewMax = OldPeakNewStart + 0.5*(WidthOldPeakNew - (WidthOldPeakNew%2));
NewPeakMax = NewPeakStart + 0.5*(WidthNewPeak - (WidthNewPeak%2));

fpeakValsTemp.at(PeakNumberWithNewPeak) = std::make_tuple(OldPeakNewMax,0,OldPeakNewStart,OldPeakNewEnd);
fpeakValsTemp.emplace_back(NewPeakMax,0,NewPeakStart,NewPeakEnd);
std::sort(fpeakValsTemp.begin(),fpeakValsTemp.end(), [](std::tuple<int,int,int,int> const &t1, std::tuple<int,int,int,int> const &t2) {return std::get<0>(t1) < std::get<0>(t2);} );

return;
}

//---------------------------------------------------------------------------------------------
double hit::DPRawHitFinder::WidthFunc(double fPeakMean,
		    		      double fPeakAmp,
		    		      double fPeakTau1,
		    		      double fPeakTau2,
				      double fStartTime,
				      double fEndTime,
			    	      double fPeakMeanTrue)
{
double MaxValue = ( fPeakAmp * exp(0.4*(fPeakMeanTrue-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(fPeakMeanTrue-fPeakMean)/fPeakTau2) );
double FuncValue = 0.;
double HalfMaxLeftTime = 0.;
double HalfMaxRightTime = 0.;
double ReBin=10.;

    //First iteration (+- 1 bin)
    for(double x = fPeakMeanTrue; x > fStartTime-1000.; x--)
    {
    FuncValue = ( fPeakAmp * exp(0.4*(x-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x-fPeakMean)/fPeakTau2) ); 
	if( FuncValue < 0.5*MaxValue )
	{
	HalfMaxLeftTime = x;
	break;
	}
    }

    for(double x = fPeakMeanTrue; x < fEndTime+1000.; x++)
    {
    FuncValue = ( fPeakAmp * exp(0.4*(x-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x-fPeakMean)/fPeakTau2) ); 
	if( FuncValue < 0.5*MaxValue )
	{
	HalfMaxRightTime = x;
	break;
	}
    }

    //Second iteration (+- 0.1 bin)
    for(double x = HalfMaxLeftTime+1; x > HalfMaxLeftTime; x-=(1/ReBin))
    {
    FuncValue = ( fPeakAmp * exp(0.4*(x-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x-fPeakMean)/fPeakTau2) ); 
	if( FuncValue < 0.5*MaxValue )
	{
	HalfMaxLeftTime = x;
	break;
	}
    }

    for(double x = HalfMaxRightTime-1; x < HalfMaxRightTime; x+=(1/ReBin))
    {
    FuncValue = ( fPeakAmp * exp(0.4*(x-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x-fPeakMean)/fPeakTau2) ); 
	if( FuncValue < 0.5*MaxValue )
	{
	HalfMaxRightTime = x;
	break;
	}
    }

    //Third iteration (+- 0.01 bin)
    for(double x = HalfMaxLeftTime+1/ReBin; x > HalfMaxLeftTime; x-=1/(ReBin*ReBin))
    {
    FuncValue = ( fPeakAmp * exp(0.4*(x-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x-fPeakMean)/fPeakTau2) ); 
	if( FuncValue < 0.5*MaxValue )
	{
	HalfMaxLeftTime = x;
	break;
	}
    }

    for(double x = HalfMaxRightTime-1/ReBin; x < HalfMaxRightTime; x+=1/(ReBin*ReBin))
    {
    FuncValue = ( fPeakAmp * exp(0.4*(x-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x-fPeakMean)/fPeakTau2) ); 
	if( FuncValue < 0.5*MaxValue )
	{
	HalfMaxRightTime = x;
	break;
	}
    }

return HalfMaxRightTime-HalfMaxLeftTime;
}

//---------------------------------------------------------------------------------------------
double hit::DPRawHitFinder::ChargeFunc(double fPeakMean,
		      		       double fPeakAmp,
		      		       double fPeakTau1,
		      		       double fPeakTau2,
				       double fChargeNormFactor,
				       double fPeakMeanTrue)

{
double ChargeSum = 0.;
double Charge = 0.;
double ReBin = 10.;

bool ChargeBigEnough=true;
    for(double x = fPeakMeanTrue - 1/ReBin; ChargeBigEnough && x > fPeakMeanTrue-1000.; x-=1.)
    {
	for(double i=0.; i > -1.; i-=(1/ReBin))
	{
    	Charge = ( fPeakAmp * exp(0.4*(x+i-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x+i-fPeakMean)/fPeakTau2) );
	ChargeSum += Charge;
    	}
	if(Charge < 0.01) ChargeBigEnough = false;
    }

ChargeBigEnough=true;
    for(double x = fPeakMeanTrue; ChargeBigEnough && x < fPeakMeanTrue+1000.; x+=1.)
    {
	for(double i=0.; i < 1.; i+=(1/ReBin))
	{
    	Charge = ( fPeakAmp * exp(0.4*(x+i-fPeakMean)/fPeakTau1)) / ( 1 + exp(0.4*(x+i-fPeakMean)/fPeakTau2) );
	ChargeSum += Charge;
    	}
	if(Charge < 0.01) ChargeBigEnough = false;
    }


return ChargeSum*fChargeNormFactor/ReBin;
}

//---------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::doBinAverage(const std::vector<float>& inputVec,
                                      std::vector<float>&       outputVec,
                                      size_t                    binsToAverage) const
{
    size_t halfBinsToAverage(binsToAverage/2);
    
    float runningSum(0.);
    
    for(size_t idx = 0; idx < halfBinsToAverage; idx++) runningSum += inputVec[idx];
    
    outputVec.resize(inputVec.size());
    std::vector<float>::iterator outputVecItr = outputVec.begin();
    
    // First pass through to build the erosion vector
    for(std::vector<float>::const_iterator inputItr = inputVec.begin(); inputItr != inputVec.end(); inputItr++)
    {
        size_t startOffset = std::distance(inputVec.begin(),inputItr);
        size_t stopOffset  = std::distance(inputItr,inputVec.end());
        size_t count       = std::min(2 * halfBinsToAverage, std::min(startOffset + halfBinsToAverage + 1, halfBinsToAverage + stopOffset - 1));
        
        if (startOffset >= halfBinsToAverage) runningSum -= *(inputItr - halfBinsToAverage);
        if (stopOffset  >  halfBinsToAverage) runningSum += *(inputItr + halfBinsToAverage);
        
        *outputVecItr++ = runningSum / float(count);
    }
    
    return;
}
    
//---------------------------------------------------------------------------------------------    
void hit::DPRawHitFinder::reBin(const std::vector<float>& inputVec,
                               std::vector<float>&       outputVec,
                               size_t                    nBinsToCombine) const
{
    size_t nNewBins = inputVec.size() / nBinsToCombine;
    
    if (inputVec.size() % nBinsToCombine > 0) nNewBins++;
    
    outputVec.resize(nNewBins, 0.);
    
    size_t outputBin = 0;
    
    for(size_t inputIdx = 0; inputIdx < inputVec.size();)
    {
        outputVec[outputBin] += inputVec[inputIdx++];
        
        if (inputIdx % nBinsToCombine == 0) outputBin++;
        
        if (outputBin > outputVec.size())
        {
            std::cout << "***** DISASTER!!! ****** outputBin: " << outputBin << ", inputIdx = " << inputIdx << std::endl;
            break;
        }
    }
    
    return;
}


  DEFINE_ART_MODULE(DPRawHitFinder)

} // end of hit namespace
#endif // DPRawHitFinder_H
