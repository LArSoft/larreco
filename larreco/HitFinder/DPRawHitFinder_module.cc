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
// The algorithm walks along the wire and looks for pulses above threshold.
// The algorithm then attempts to fit n double exponentials to these 
// pulses where n is set by the number of peaks found in the pulse.
// If the Chi2/NDF returned is "bad" it attempts to fit n+1 double exponentials
// to the pulse. If this is a better fit it then uses the parameters of this
// fit to characterize the "hit" object. The parameters of the fit are saved
// in a feature vector by using MVAWriter.
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
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "HitFilterAlg.h"

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
    void reconfigure(fhicl::ParameterSet const& p) override;

  private:

    using TimeValsVec      = std::vector<std::tuple<int,int,int>>;
    using PeakTimeWidVec   = std::vector<std::tuple<int,int,int,int>>;
    using MergedTimeWidVec = std::vector<std::tuple<int,int,PeakTimeWidVec>>;
    
    void findCandidatePeaks(std::vector<float>::const_iterator startItr,
                            std::vector<float>::const_iterator stopItr,
                            TimeValsVec&                       timeValsVec,
                            int&                               PeakMin,
                            int                                firstTick) const;
    
    void mergeCandidatePeaks(const std::vector<float>&, TimeValsVec&, MergedTimeWidVec&) const;
  
    // ### This function will fit N-Exponentials to a TH1D where N is set ###
    // ###            by the number of peaks found in the pulse         ###
      
    using ParameterVec = std::vector<std::pair<double,double>>;  //< parameter/error vec
    
    void FitExponentials(const std::vector<float>& SignalVector,
                      const PeakTimeWidVec&     PeakVals,
                      int                       StartTime,
                      int                       EndTime,
                      double                    ampScaleFctr,
                      ParameterVec&             paramVec,
                      double&                   chi2PerNDF,
                      int&                      NDF);
    
    void FillOutHitParameterVector(const std::vector<double>& input,
				   std::vector<double>& output);
      
    void doBinAverage(const std::vector<float>& inputVec,
                      std::vector<float>&       outputVec,
                      size_t                    binsToAverage) const;
      
    void reBin(const std::vector<float>& inputVec,
               std::vector<float>&       outputVec,
               size_t                    nBinsToCombine) const;
    
    std::string      fCalDataModuleLabel;

    int    fMinSig;                   ///<signal height threshold
    double fInitTau;                ///<Initial value for tau1 and tau2
    int    fMinWidth;                 ///<Minimum hit width
    int    fLongMaxHits;              ///<Maximum number hits on a really long pulse train
    int    fLongPulseWidth;           ///<Sets width of hits used to describe long pulses
    
    unsigned int fMaxMultiHit;              ///<maximum hits for multi fit
    int          fAreaMethod;               ///<Type of area calculation
    double	 fAreaNorm;                ///<factors for converting area to same units as peak height
    bool	 fTryNplus1Fits;            ///<whether we will (true) or won't (false) try n+1 fits
    double	 fChi2NDFRetry;             ///<Value at which a second n+1 Fit will be tried
    double	 fChi2NDF;                  ///maximum Chisquared / NDF allowed for a hit to be saved
    size_t       fNumBinsToAverage;         ///< If bin averaging for peak finding, number bins to average
    double	 fNoiseRMS = 0;
    double	 fMIPSignal = 0;
    
    std::unique_ptr<HitFilterAlg> fHitFilterAlg;   ///algorithm used to filter out noise hits

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
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
  
    
    bool const doHitFiltering = p.get<bool>("FilterHits", false);
    if (doHitFiltering) {
      if (fHitFilterAlg) { // create a new algorithm instance
        fHitFilterAlg->reconfigure(p.get<fhicl::ParameterSet>("HitFilterAlg"));
      }
      else { // reconfigure the existing instance
        fHitFilterAlg = std::make_unique<HitFilterAlg>
          (p.get<fhicl::ParameterSet>("HitFilterAlg"));
      }
    }

    fLongMaxHits      = p.get< int    >("LongMaxHits");
    fLongPulseWidth   = p.get< int    >("LongPulseWidth");
    fMaxMultiHit      = p.get< int    >("MaxMultiHit");
    fAreaMethod       = p.get< int    >("AreaMethod");
    fTryNplus1Fits    = p.get< bool   >("TryNplus1Fits");
    fChi2NDFRetry     = p.get< double >("Chi2NDFRetry");
    fChi2NDF          = p.get< double >("Chi2NDF");
    fNumBinsToAverage = p.get< size_t >("NumBinsToAverage", 0);
    fNoiseRMS         = p.get< double >("NoiseRMS");
    fMIPSignal        = p.get< double >("MIPSignal");
    fMinSig           = p.get< int    >("MinSig");
    fMinWidth         = p.get< double >("MinWidth");
    fInitTau	      = p.get< double >("InitTau");
    fAreaNorm         = p.get< double >("AreaNorm");
}  

//-------------------------------------------------
//-------------------------------------------------
void DPRawHitFinder::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
   
    
    // ======================================
    // === Hit Information for Histograms ===
    fFirstChi2	= tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
    fChi2	        = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
}

//-------------------------------------------------
//-------------------------------------------------
void DPRawHitFinder::endJob()
{

}

//  This algorithm uses the fact that deconvolved signals are very smooth 
//  and looks for hits as areas between local minima that have signal above 
//  threshold.
//-------------------------------------------------
void DPRawHitFinder::produce(art::Event& evt)
{
    //==================================================================================================
    TH1::AddDirectory(kFALSE);
   
    // Instantiate and Reset a stop watch
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
    art::FindOneP<raw::RawDigit> RawDigits
        (wireVecHandle, evt, fCalDataModuleLabel);
    // Channel Number
    raw::ChannelID_t channel = raw::InvalidChannelID;
    
    //#################################################
    //###    Set the charge determination method    ###
    //### Default is to compute the normalized area ###
    //#################################################
    //std::function<double (double,double,double,double,int,int)> chargeFunc = [](double peakMean, double peakAmp, double peaktau1, double peaktau2, double areaNorm, int low, int hi){return std::sqrt(2*TMath::Pi())*peakAmp*peakWidth/areaNorm;};
    
    //##############################################
    //### Alternative is to integrate over pulse ###
    //##############################################
    //if (fAreaMethod == 0)
        std::function<double (double,double,double,double,double,int,int)> chargeFunc = [](double peakMean, double peakAmp, double peaktau1, double peaktau2, double areaNorm, int low, int hi)
                        {
                            double charge(0);
			    double ReBin = 10.;
                            for(int x = low; x < hi*ReBin; x++)
                                charge += ( peakAmp * exp(0.4*((x/ReBin)-peakMean)/peaktau1)) / ( 1 + exp(0.4*((x/ReBin)-peakMean)/peaktau2) );
                            return charge/ReBin;
                        };
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

        // #################################################
        // ### Set up to loop over ROI's for this wire   ###
        // #################################################
        const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();
       
        for(const auto& range : signalROI.get_ranges())
        {
            // #################################################
            // ### Getting a vector of signals for this wire ###
            // #################################################
            //std::vector<float> signal(wire->Signal());

            const std::vector<float>& signal = range.data();
      
            // ##########################################################
            // ### Making an iterator for the time ticks of this wire ###
            // ##########################################################
            std::vector<float>::const_iterator timeIter;  	    // iterator for time bins
           
            // ROI start time
            raw::TDCtick_t roiFirstBinTick = range.begin_index();
 	    //std::cout << "roiFirstBinTick: " << roiFirstBinTick << std::endl;
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
            // ### Lets loop over the pulses we found on this wire ###
            // #######################################################
            
            for(auto& mergedCands : mergedVec)
            {
                int             startT   = std::get<0>(mergedCands);
                int             endT     = std::get<1>(mergedCands);
                PeakTimeWidVec& peakVals = std::get<2>(mergedCands);

                // ### Putting in a protection in case things went wrong ###
                // ### In the end, this primarily catches the case where ###
                // ### a fake pulse is at the start of the ROI           ###
                if (endT - startT < fMinWidth) continue;
	 
                // #######################################################
                // ### Clearing the parameter vector for the new pulse ###
                // #######################################################
                // === Setting the number of Exponential functions to try ===
                int nExponentialsForFit = peakVals.size();
	 
                // #####################################################
                // ### Calling the function for fitting Exponentials ###
                // #####################################################
                double       chi2PerNDF(0.);
                int          NDF(0);
                ParameterVec paramVec;
                
                // ##########################################################
                // ### If # requested Exponentials is too large then punt ###
                // ##########################################################

                if (peakVals.size() <= fMaxMultiHit)
                {
                    FitExponentials(signal, peakVals, startT, endT, 1.0, paramVec, chi2PerNDF, NDF);
               
                    // If the chi2 is infinite then there is a real problem so we bail
                    if (!(chi2PerNDF < std::numeric_limits<double>::infinity())) continue;
                   
                    fFirstChi2->Fill(chi2PerNDF);

                    // #######################################################
                    // ### Clearing the parameter vector for the new pulse ###
                    // #######################################################
              /*    double       chi2PerNDF2(0.);
                    int          NDF2(0);
                    ParameterVec paramVec2;
                
                    // ########################################################
                    // ### Trying extra Exponentials for an initial bad fit ###
                    // ########################################################
                    if( (chi2PerNDF > (fChi2NDFRetry) && fTryNplus1Fits && nExponentialsForFit == 1)||
                        (chi2PerNDF > (fChi2NDFRetry)   && fTryNplus1Fits && nExponentialsForFit > 1))
                    {
                        // ###############################################################
                        // ### Modify input parameters for re-fitting n+1 Exponentials ###
                        // ###############################################################
           		for(auto& peakValsIter : peakVals)
			{
			
			}
             		peakVals.emplace_back(maxT,widT,nextStartT,endT);
                    
                        // We need to make sure we are not out of range and new peak amplitude is non-negative
                        if (newPeakTime < endT - 1 && signal[newPeakTime] > 0.)
                        {
                            peakVals.emplace_back(newPeakTime,999,newPeakTimeLowB,newPeakTimeUpB,0,0);
	    
                            // #########################################################
                            // ### Calling the function for re-fitting n+1 Exponentials ###
                            // #########################################################
                            FitExponentials(signal, peakVals, startT, endT, 0.5, paramVec2, chi2PerNDF2, NDF2);
	    
                            // #########################################################
                            // ### Getting the appropriate parameter into the vector ###
                            // #########################################################

                            if (chi2PerNDF2 < chi2PerNDF)
                            {
                                nExponentialsForFit = peakVals.size();
                                chi2PerNDF  = chi2PerNDF2;
                                NDF         = NDF2;
                                paramVec    = paramVec2;
                            }
                        }
                    } */
                }
                else{std::cout << "hit too long" << std::endl;}
                // #######################################################
                // ### If too large then force alternate solution      ###
                // ### - Make n hits from pulse train where n will     ###
                // ###   depend on the fhicl parameter fLongPulseWidth ###
                // ### Also do this if chi^2 is too large              ###
                // #######################################################
     /*           if (peakVals.size() > fMaxMultiHit || chi2PerNDF > fChi2NDF)
                {
                    int longPulseWidth = fLongPulseWidth.at(view);
                    int nHitsThisPulse = (endT - startT) / longPulseWidth;
                    
                    if (nHitsThisPulse > fLongMaxHits.at(view))
                    {
                        nHitsThisPulse = fLongMaxHits.at(view);
                        longPulseWidth = (endT - startT) / nHitsThisPulse;
                    }
                    
                    if (nHitsThisPulse * longPulseWidth < endT - startT) nHitsThisPulse++;
                    
                    int firstTick = startT;
                    int lastTick  = firstTick + std::min(endT,longPulseWidth);
                    
                    paramVec.clear();
                    nExponentialsForFit = nHitsThisPulse;
                    NDF         = 1.;
                    chi2PerNDF  =  chi2PerNDF > fChi2NDF ? chi2PerNDF : -1.;
                    
                    for(int hitIdx = 0; hitIdx < nHitsThisPulse; hitIdx++)
                    {
                        // This hit parameters
                        double sumADC    = std::accumulate(signal.begin() + firstTick, signal.begin() + lastTick, 0.);
                        double peakSigma = (lastTick - firstTick) / 3.;  // Set the width...
                        double peakAmp   = 0.3989 * sumADC / peakSigma;  // Use Exponentials formulation
                        double peakMean  = (firstTick + lastTick) / 2.;
                    
                        // Store hit params
                        paramVec.emplace_back(peakAmp,   0.1 * peakAmp);
                        paramVec.emplace_back(peakMean,  0.1 * peakMean);
                        paramVec.emplace_back(peakSigma, 0.1 * peakSigma);
                        
                        // set for next loop
                        firstTick = lastTick;
                        lastTick  = std::min(lastTick  + longPulseWidth, endT);
                    }
                } */
	    
                // #######################################################
                // ### Loop through returned peaks and make recob hits ###
                // #######################################################
                
                int numHits(0);

                for(int hitIdx = 0; hitIdx < nExponentialsForFit; hitIdx++)
                {
                    // Extract values for this hit
                    double peakAmp   = paramVec[4*hitIdx].first;
                    double peakMean  = paramVec[4*hitIdx+1].first;
		    double peakTau1 = paramVec[4*hitIdx+2].first;
		    double peakTau2 = paramVec[4*hitIdx+3].first;
                    double peakWidth = std::max(8.175, 5.61959 + 7.93513*peakTau1 - 1.45595*peakTau1*peakTau1)/(2*sqrt(2*log(2))); //Conversion from tau1 to "width at 0.5*amplitude" to standard deviatiopn of a Gaussian. Width of elec. response (= min. width of a hit) is 8.175 at 0.5*amplitude.
                    // Extract errors
                    double peakAmpErr   = paramVec[4*hitIdx].second;
                    double peakMeanErr  = paramVec[4*hitIdx+1].second;
                    double peakWidthErr = sqrt(pow(paramVec[4*hitIdx+2].second,2) + pow(paramVec[4*hitIdx+3].second,2));



                    // ### Charge ###
                    double charge    = chargeFunc(peakMean, peakAmp, peakTau1, peakTau2, fAreaNorm,startT-10,endT+20);
                    double chargeErr = std::sqrt(TMath::Pi()) * (peakAmpErr*peakWidthErr + peakWidthErr*peakAmpErr);    

                    // ### limits for getting sums
	            int startTthisHit = std::get<2>(peakVals.at(hitIdx));
	            int endTthisHit = std::get<3>(peakVals.at(hitIdx));
                    std::vector<float>::const_iterator sumStartItr = signal.begin() + startTthisHit;
                    std::vector<float>::const_iterator sumEndItr   = signal.begin() + endTthisHit;

                    // ### Sum of ADC counts
                    double sumADC = std::accumulate(sumStartItr, sumEndItr, 0.);

		    //Mean
		    TF1 Exponentials("Exponentials","( [0] * exp(0.4*(x-[1])/[2]) / ( 1 + exp(0.4*(x-[1])/[3]) ) )",startTthisHit,endTthisHit);
        	    Exponentials.SetParameter(0, peakAmp);
        	    Exponentials.SetParameter(1, peakMean);
        	    Exponentials.SetParameter(2, peakTau1);
        	    Exponentials.SetParameter(3, peakTau2);
		    double peakMeanTrue = Exponentials.GetMaximumX(startTthisHit,endTthisHit);
		    Exponentials.Delete();

		    //std::cout << "peakMean: " << peakMean << std::endl;
		    //std::cout << "peakMeanTrue: " << peakMeanTrue << std::endl;
                    // ok, now create the hit
                    recob::HitCreator hitcreator(*wire,                            // wire reference
                                                 wid,                              // wire ID
                                                 startT+roiFirstBinTick,           // start_tick TODO check
                                                 endT+roiFirstBinTick,             // end_tick TODO check
                                                 peakWidth,                        // rms
                                                 peakMeanTrue+roiFirstBinTick,         // peak_time
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
                    
                    const recob::Hit hit(hitcreator.move());
		    
                    if (!fHitFilterAlg || fHitFilterAlg->IsGoodHit(hit)) {
                        hcol.emplace_back(std::move(hit), wire, rawdigits);
                        // add fit parameters associated to the hit just pushed to the collection
                        std::array<float, 3> fitParams;
                        fitParams[0] = peakMean+roiFirstBinTick;
                        fitParams[1] = paramVec[4*hitIdx+2].first;
                        fitParams[2] = paramVec[4*hitIdx+3].first;
                        fHitParamWriter.addVector(hitID, fitParams);
                        numHits++;
                    }
                } // <---End loop over Exponentials
                
                fChi2->Fill(chi2PerNDF);
	    
           }//<---End loop over merged candidate hits

       } //<---End looping over ROI's
	//std::cout << std::endl;	 

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
                                            int&                                  PeakMin,
                                            int                                   firstTick) const
{
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        auto maxItr = std::max_element(startItr, stopItr);
        
        float maxValue = *maxItr;
        int   maxTime  = std::distance(startItr,maxItr);
        
        if (maxValue > PeakMin)
        {
            // backwards to find first bin for this candidate hit
            auto firstItr = std::distance(startItr,maxItr) > 2 ? maxItr - 1 : startItr;

            while(firstItr != startItr)
            {    
                // Check both sides of firstItr and look for min/inflection point
		if ( *firstItr <= *(firstItr+1) && *firstItr < (*(firstItr-1)) && *firstItr < (*(firstItr-2)) && *firstItr < (*(firstItr-3)) ) break;            
                if (*firstItr < 0) break;

                firstItr--;
            }

            int firstTime = std::distance(startItr,firstItr);
            
            // Recursive call to find all candidate hits earlier than this peak
            findCandidatePeaks(startItr, firstItr - 1, timeValsVec, PeakMin, firstTick);
            
            // forwards to find last bin for this candidate hit
            auto lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;

            while(lastItr != stopItr)
            {        
                if ( *lastItr <= *(lastItr-1) && *lastItr < (*(lastItr+1)) && *lastItr < (*(lastItr+2)) && *lastItr < (*(lastItr+3)) ) break;       
                if (*lastItr < 0) break;
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

    while(timeValsVecItr != timeValsVec.end())
    {
        PeakTimeWidVec peakVals;
        
        // Setting the start, peak, and end time of the pulse
        auto& timeVal = *timeValsVecItr++;
        int startT = std::get<0>(timeVal);
        int maxT   = std::get<1>(timeVal);
        int endT   = std::get<2>(timeVal);
        int widT   = std::max(2,(endT - startT) / 6);

        peakVals.emplace_back(maxT,widT,startT,endT);
        
	std::cout << "startT: " << startT << std::endl;
	std::cout << "endT: " << endT << std::endl;

        // See if we want to merge pulses together
        // First check if we have more than one pulse on the wire
        bool checkNextHit = timeValsVecItr != timeValsVec.end();

        // Loop until no more merged pulses (or candidates in this ROI)
        while(checkNextHit)
        {
            // If the start time of the next pulse is the end time of the current pulse then merge
            // Alternatively, if the start time of the next pulses is one tick away then
            // merge if the intervening signal is above 0.
            int nextStartT = std::get<0>(*timeValsVecItr);
	    int PrevStartT=startT;
	    int PrevMaxT=maxT;
	    int PrevEndT=endT;
	    int PrevWidT=widT;

	    
	    bool NegBinBetweenPeaks = 0;
	    for(int i = endT; i <= nextStartT; i++)
	    {
		if(signalVec[i]<0)
		{
		NegBinBetweenPeaks = 1;
		break;
		}
	    }

            //if( ( ( nextStartT - endT <= 1 ) || ( signalVec[endT] >= 0 && signalVec[nextStartT] >= 0 ) ) && PeaksInThisMergedPeak < fMaxMultiHit-1 )
            if( ( ( nextStartT - endT <= 1 ) || NegBinBetweenPeaks == 0 ) && PeaksInThisMergedPeak < fMaxMultiHit-1 )
            {
                timeVal = *timeValsVecItr++;
                maxT = std::get<1>(timeVal);
                endT = std::get<2>(timeVal);
                widT = std::max(2,(endT - nextStartT) / 6);

		std::cout << "nextStartT: " << nextStartT << std::endl;
		std::cout << "nextEndT: " << endT << std::endl;

		if(  maxT-PrevEndT <= 1 && (signalVec[maxT]-signalVec[PrevEndT]) < 0.2*(signalVec[PrevMaxT]-signalVec[PrevEndT]) )
		{
		    peakVals.pop_back();  
                    peakVals.emplace_back(PrevMaxT,PrevWidT,PrevStartT,endT);
		}
		else
		{
                    peakVals.emplace_back(maxT,widT,nextStartT,endT);
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
        mergedVec.emplace_back(startT, endT, peakVals);
    }
    
    return;
}

// --------------------------------------------------------------------------------------------
// Fit Exponentialssians
// --------------------------------------------------------------------------------------------
void hit::DPRawHitFinder::FitExponentials(const std::vector<float>& SignalVector,
                                      const PeakTimeWidVec&     PeakVals,
                                      int                       StartTime,
                                      int                       EndTime,
                                      double                    ampScaleFctr,
                                      ParameterVec&             paramVec,
                                      double&                   chi2PerNDF,
                                      int&                      NDF)
{
    int size = (EndTime) - (StartTime);
    // #############################################
    // ### If size < 0 then set the size to zero ###
    // #############################################
    if(EndTime - StartTime < 0){size = 0;}
   
    // --- TH1D HitSignal ---
    TH1F hitSignal("hitSignal","",std::max(size,1),StartTime,EndTime);
    hitSignal.Sumw2();
   
    // #############################
    // ### Filling the histogram ###
    // #############################
    for(int aa = StartTime; aa < EndTime; aa++)
    {
        hitSignal.Fill(aa,SignalVector[aa]);
        if(EndTime > 10000){break;} // FIXME why?
    }//<---End aa loop
    // ############################################
    // ### Building TFormula for basic Exponentials ###
    // ############################################
    std::string eqn = "( [0] * exp(0.4*(x-[1])/[2]) / ( 1 + exp(0.4*(x-[1])/[3]) ) )";  // string for equation for Exponentials fit
    std::stringstream numConv;
    
    int NPeaks = PeakVals.size();

    for(int i = 4; i < 4 + (NPeaks-1)*2; i+=2)
    {
        eqn.append("+( [");
        numConv.str("");
        numConv << i;
        eqn.append(numConv.str());
        eqn.append("] * exp(0.4*(x-[");
        numConv.str("");
        numConv << i+1;
        eqn.append(numConv.str());
        eqn.append("])/[");
        numConv.str("");
        numConv << 2;
        eqn.append(numConv.str());    
        eqn.append("]) / ( 1 + exp(0.4*(x-[");
        numConv.str("");
        numConv << i+1;
        eqn.append(numConv.str()); 
        eqn.append("])/[");       
        numConv.str("");
        numConv << 3;
        eqn.append(numConv.str()); 
    eqn.append("]) ) )");
    }
    
    // ---------------------------------
    // --- TF1 function for ExponentialsHit  ---
    // ---------------------------------
    TF1 Exponentials("Exponentials",eqn.c_str(),0,std::max(size,1));
   

        double peakMean = std::get<0>(PeakVals.at(0));
        double amplitude  = SignalVector[peakMean];
	double peakMeanShift=2.;

        Exponentials.SetParameter(0, 1.65*amplitude);
        Exponentials.SetParameter(1, peakMean-peakMeanShift);
        Exponentials.SetParameter(2, 0.5);
        Exponentials.SetParameter(3, 0.5);
        Exponentials.SetParLimits(0, 0.001*amplitude, 10*1.65*amplitude);
        Exponentials.SetParLimits(1, peakMean-peakMeanShift-5., peakMean-peakMeanShift+5.);

		for(int k = 4, j=1; k < 4 + (NPeaks-1)*2; k+=2, j++)
    		{
        	peakMean = std::get<0>(PeakVals.at(j));
        	amplitude = SignalVector[peakMean];
		Exponentials.SetParameter(k,1.65*amplitude);
		Exponentials.SetParameter(k+1,peakMean-peakMeanShift);
       		
		Exponentials.SetParLimits(k, 0.5*amplitude, 2*1.65*amplitude);
		Exponentials.SetParLimits(k+1, peakMean-peakMeanShift-5., peakMean-peakMeanShift+5.);
			if(std::get<1>(PeakVals.at(j))==999)
			{
				Exponentials.SetParLimits(k+1, std::get<2>(PeakVals.at(j)), std::get<3>(PeakVals.at(j)));
				Exponentials.SetParLimits(k, 0.1*amplitude, 2*1.65*amplitude);
			}
		}


    // ####################################################
    // ### PERFORMING THE TOTAL Exponentials FIT OF THE HIT ###
    // ####################################################
    try
      { hitSignal.Fit(&Exponentials,"QNRWM","", StartTime, EndTime);}
    catch(...)
      {mf::LogWarning("DPRawHitFinder") << "Fitter failed finding a hit";}
   
    // ##################################################
    // ### Getting the fitted parameters from the fit ###
    // ##################################################
    chi2PerNDF = (Exponentials.GetChisquare() / Exponentials.GetNDF());
    NDF        = Exponentials.GetNDF();
  

    //std::cout << "chi2PerNDF: " << chi2PerNDF << "\t" << "NPeaks: " << NPeaks <<  std::endl;
    //std::cout << "chi2PerNDF/NPeaks: " << chi2PerNDF/NPeaks << std::endl;


    paramVec.emplace_back(Exponentials.GetParameter(0),Exponentials.GetParError(0));
    paramVec.emplace_back(Exponentials.GetParameter(1),Exponentials.GetParError(1));
    paramVec.emplace_back(Exponentials.GetParameter(2),Exponentials.GetParError(2));
    paramVec.emplace_back(Exponentials.GetParameter(3),Exponentials.GetParError(3));
    
    for(int i = 4; i < (4 + 2*(NPeaks-1)); i+=2)
    {
        paramVec.emplace_back(Exponentials.GetParameter(i),Exponentials.GetParError(i));
        paramVec.emplace_back(Exponentials.GetParameter(i+1),Exponentials.GetParError(i));
        paramVec.emplace_back(Exponentials.GetParameter(2),Exponentials.GetParError(2));
        paramVec.emplace_back(Exponentials.GetParameter(3),Exponentials.GetParError(3));
    }
    Exponentials.Delete();
    hitSignal.Delete();
}//<----End FitExponentials

    
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
