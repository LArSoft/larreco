////////////////////////////////////////////////////////////////////////
//
// GaussHitFinder class
//
// jaasaadi@syr.edu
//
//  This algorithm is designed to find hits on wires after deconvolution.
// -----------------------------------
// This algorithm is based on the FFTHitFinder written by Brian Page,
// Michigan State University, for the ArgoNeuT experiment.
//
//
// The algorithm walks along the wire and looks for pulses above threshold
// The algorithm then attempts to fit n-gaussians to these pulses where n
// is set by the number of peaks found in the pulse
// If the Chi2/NDF returned is "bad" it attempts to fit n+1 gaussians to
// the pulse. If this is a better fit it then uses the parameters of the
// Gaussian fit to characterize the "hit" object
//
// To use this simply include the following in your producers:
// gaushit:     @local::microboone_gaushitfinder
// gaushit:	@local::argoneut_gaushitfinder
////////////////////////////////////////////////////////////////////////


// C/C++ standard library
#include <algorithm> // std::accumulate()
#include <string>
#include <memory> // std::unique_ptr()
#include <utility> // std::move()

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"

// LArSoft Includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "HitFilterAlg.h"

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"

// ROOT Includes
#include "TH1F.h"
#include "TMath.h"

namespace hit{
class GausHitFinder : public art::EDProducer {

public:

    explicit GausHitFinder(fhicl::ParameterSet const& pset);

private:

    void produce(art::Event& evt) override;
    void beginJob() override;

    void FillOutHitParameterVector(const std::vector<double>& input, std::vector<double>& output);

    bool                fFilterHits;
    bool                fMakeRawDigitAssns;

    std::string         fCalDataModuleLabel;
    std::string         fAllHitsInstanceName;

    std::vector<int>    fLongMaxHitsVec;           ///<Maximum number hits on a really long pulse train
    std::vector<int>    fLongPulseWidthVec;        ///<Sets width of hits used to describe long pulses

    size_t              fMaxMultiHit;              ///<maximum hits for multi fit
    int                 fAreaMethod;               ///<Type of area calculation
    std::vector<double> fAreaNormsVec;             ///<factors for converting area to same units as peak height
    double	            fChi2NDF;                  ///maximum Chisquared / NDF allowed for a hit to be saved

    std::vector<float>  fPulseHeightCuts;
    std::vector<float>  fPulseWidthCuts;
    std::vector<float>  fPulseRatioCuts;

    size_t              fEventCount;

    std::vector<std::unique_ptr<reco_tool::ICandidateHitFinder>> fHitFinderToolVec;  ///< For finding candidate hits
    std::unique_ptr<reco_tool::IPeakFitter>                      fPeakFitterTool;    ///< Perform fit to candidate peaks
    std::unique_ptr<HitFilterAlg>                                fHitFilterAlg;      ///< algorithm used to filter out noise hits

    TH1F* fFirstChi2;
    TH1F* fChi2;

}; // class GausHitFinder


//-------------------------------------------------
//-------------------------------------------------
GausHitFinder::GausHitFinder(fhicl::ParameterSet const& pset)
  : EDProducer{pset},
    fEventCount(0)
{
    fCalDataModuleLabel  = pset.get< std::string >("CalDataModuleLabel");
    fAllHitsInstanceName = pset.get< std::string >("AllHitsInstanceName","");
    fFilterHits          = pset.get< bool        >("FilterHits",false);
    fMakeRawDigitAssns   = pset.get< bool        >("MakeRawDigitAssns",true);

    if (fFilterHits) {
        fHitFilterAlg = std::make_unique<HitFilterAlg>(pset.get<fhicl::ParameterSet>("HitFilterAlg"));
    }

    FillOutHitParameterVector(pset.get< std::vector<double> >("AreaNorms"), fAreaNormsVec);

    fLongMaxHitsVec    = pset.get< std::vector<int>>("LongMaxHits",    std::vector<int>() = {25,25,25});
    fLongPulseWidthVec = pset.get< std::vector<int>>("LongPulseWidth", std::vector<int>() = {16,16,16});
    fMaxMultiHit       = pset.get< int             >("MaxMultiHit");
    fAreaMethod        = pset.get< int             >("AreaMethod");
    fChi2NDF           = pset.get< double          >("Chi2NDF");

    fPulseHeightCuts   = pset.get< std::vector<float>>("PulseHeightCuts", std::vector<float>() = {3.0,  3.0,  3.0});
    fPulseWidthCuts    = pset.get< std::vector<float>>("PulseWidthCuts",  std::vector<float>() = {2.0,  1.5,  1.0});
    fPulseRatioCuts    = pset.get< std::vector<float>>("PulseRatioCuts",  std::vector<float>() = {0.35, 0.40, 0.20});

    // recover the tool to do the candidate hit finding
    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& hitFinderTools = pset.get<fhicl::ParameterSet>("HitFinderToolVec");

    fHitFinderToolVec.resize(hitFinderTools.get_pset_names().size());

    for(const std::string& hitFinderTool : hitFinderTools.get_pset_names())
    {
        const fhicl::ParameterSet& hitFinderToolParamSet = hitFinderTools.get<fhicl::ParameterSet>(hitFinderTool);
        size_t                     planeIdx              = hitFinderToolParamSet.get<size_t>("Plane");

        fHitFinderToolVec.at(planeIdx) = art::make_tool<reco_tool::ICandidateHitFinder>(hitFinderToolParamSet);
    }

    // Recover the peak fitting tool
    fPeakFitterTool = art::make_tool<reco_tool::IPeakFitter>(pset.get<fhicl::ParameterSet>("PeakFitter"));

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // We want the option to output two hit collections, one filtered
    // and one with all hits. The key to doing this will be a non-null
    // instance name for the second collection
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this,fAllHitsInstanceName,true,fMakeRawDigitAssns);

    // and now the filtered hits...
    if (fAllHitsInstanceName != "") recob::HitCollectionCreator::declare_products(*this,"",true,fMakeRawDigitAssns);

    return;
} // GausHitFinder::GausHitFinder()


//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::FillOutHitParameterVector(const std::vector<double>& input,
                                              std::vector<double>&       output)
{
    if(input.size()==0)
        throw std::runtime_error("GausHitFinder::FillOutHitParameterVector ERROR! Input config vector has zero size.");

    art::ServiceHandle<geo::Geometry const> geom;
    const unsigned int N_PLANES = geom->Nplanes();

    if(input.size()==1)
        output.resize(N_PLANES,input[0]);
    else if(input.size()==N_PLANES)
        output = input;
    else
        throw std::runtime_error("GausHitFinder::FillOutHitParameterVector ERROR! Input config vector size !=1 and !=N_PLANES.");

}



//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;


    // ======================================
    // === Hit Information for Histograms ===
    fFirstChi2	= tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
    fChi2	        = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
}

//  This algorithm uses the fact that deconvolved signals are very smooth
//  and looks for hits as areas between local minima that have signal above
//  threshold.
//-------------------------------------------------
void GausHitFinder::produce(art::Event& evt)
{
    //==================================================================================================

    TH1::AddDirectory(kFALSE);

    // Instantiate and Reset a stop watch
    //TStopwatch StopWatch;
    //StopWatch.Reset();

    // ################################
    // ### Calling Geometry service ###
    // ################################
    art::ServiceHandle<geo::Geometry const> geom;

    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // this contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator allHitCol(*this, evt, fAllHitsInstanceName, true, fMakeRawDigitAssns);

    // Handle the filtered hits collection...
    recob::HitCollectionCreator  hcol(*this, evt, "", true, fMakeRawDigitAssns);
    recob::HitCollectionCreator* filteredHitCol = 0;

    if( fFilterHits ) filteredHitCol = &hcol;
//    if (fAllHitsInstanceName != "") filteredHitCol = &hcol;

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
    std::function<double (double,double,double,double,int,int)> chargeFunc = [](double peakMean, double peakAmp, double peakWidth, double areaNorm, int low, int hi){return std::sqrt(2*TMath::Pi())*peakAmp*peakWidth/areaNorm;};

    //##############################################
    //### Alternative is to integrate over pulse ###
    //##############################################
    if (fAreaMethod == 0)
        chargeFunc = [](double peakMean, double peakAmp, double peakWidth, double areaNorm, int low, int hi)
                        {
                            double charge(0);
                            for(int sigPos = low; sigPos < hi; sigPos++)
                                charge += peakAmp * TMath::Gaus(sigPos,peakMean,peakWidth);
                            return charge;
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

        // --- Setting Channel Number and Signal type ---
        channel = wire->Channel();

        // get the WireID for this hit
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        // for now, just take the first option returned from ChannelToWire
        geo::WireID wid  = wids[0];
        // We need to know the plane to look up parameters
        geo::PlaneID::PlaneID_t plane = wid.Plane;

        // ----------------------------------------------------------
        // -- Setting the appropriate signal widths and thresholds --
        // --    for the right plane.      --
        // ----------------------------------------------------------

        // #################################################
        // ### Set up to loop over ROI's for this wire   ###
        // #################################################
        const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();

        for(const auto& range : signalROI.get_ranges())
        {
            // ##########################################################
            // ### Making an iterator for the time ticks of this wire ###
            // ##########################################################
            std::vector<float>::const_iterator timeIter;  	    // iterator for time bins

            // ROI start time
            raw::TDCtick_t roiFirstBinTick = range.begin_index();

            // ###########################################################
            // ### Scan the waveform and find candidate peaks + merge  ###
            // ###########################################################

            reco_tool::ICandidateHitFinder::HitCandidateVec      hitCandidateVec;
            reco_tool::ICandidateHitFinder::MergeHitCandidateVec mergedCandidateHitVec;

            fHitFinderToolVec.at(plane)->findHitCandidates(range, 0, channel, fEventCount, hitCandidateVec);
            fHitFinderToolVec.at(plane)->MergeHitCandidates(range, hitCandidateVec, mergedCandidateHitVec);

            // #######################################################
            // ### Lets loop over the pulses we found on this wire ###
            // #######################################################

            for(auto& mergedCands : mergedCandidateHitVec)
            {
                int startT= mergedCands.front().startTick;
                int endT  = mergedCands.back().stopTick;

                // ### Putting in a protection in case things went wrong ###
                // ### In the end, this primarily catches the case where ###
                // ### a fake pulse is at the start of the ROI           ###
                if (endT - startT < 5) continue;

                // #######################################################
                // ### Clearing the parameter vector for the new pulse ###
                // #######################################################

                // === Setting the number of Gaussians to try ===
                int nGausForFit = mergedCands.size();

                // ##################################################
                // ### Calling the function for fitting Gaussians ###
                // ##################################################
                double                                chi2PerNDF(0.);
                int                                   NDF(1);
		/*stand alone
                reco_tool::IPeakFitter::PeakParamsVec peakParamsVec(nGausForFit);
		*/
                reco_tool::IPeakFitter::PeakParamsVec peakParamsVec;

                // #######################################################
                // ### If # requested Gaussians is too large then punt ###
                // #######################################################
                if (mergedCands.size() <= fMaxMultiHit)
                {
                    fPeakFitterTool->findPeakParameters(range.data(), mergedCands, peakParamsVec, chi2PerNDF, NDF);

                    // If the chi2 is infinite then there is a real problem so we bail
                    if (!(chi2PerNDF < std::numeric_limits<double>::infinity()))
                    {
                        chi2PerNDF = 2.*fChi2NDF;
                        NDF        = 2;
                    }

                    fFirstChi2->Fill(chi2PerNDF);
                }

                // #######################################################
                // ### If too large then force alternate solution      ###
                // ### - Make n hits from pulse train where n will     ###
                // ###   depend on the fhicl parameter fLongPulseWidth ###
                // ### Also do this if chi^2 is too large              ###
                // #######################################################
                if (mergedCands.size() > fMaxMultiHit || nGausForFit * chi2PerNDF > fChi2NDF)
                {
                    int longPulseWidth = fLongPulseWidthVec.at(plane);
                    int nHitsThisPulse = (endT - startT) / longPulseWidth;

                    if (nHitsThisPulse > fLongMaxHitsVec.at(plane))
                    {
                        nHitsThisPulse = fLongMaxHitsVec.at(plane);
                        longPulseWidth = (endT - startT) / nHitsThisPulse;
                    }

                    if (nHitsThisPulse * longPulseWidth < endT - startT) nHitsThisPulse++;

                    int firstTick = startT;
                    int lastTick  = std::min(firstTick + longPulseWidth, endT);

                    peakParamsVec.clear();
                    nGausForFit = nHitsThisPulse;
                    NDF         = 1.;
                    chi2PerNDF  =  chi2PerNDF > fChi2NDF ? chi2PerNDF : -1.;

                    for(int hitIdx = 0; hitIdx < nHitsThisPulse; hitIdx++)
                    {
                        // This hit parameters
                        double sumADC    = std::accumulate(range.begin() + firstTick, range.begin() + lastTick, 0.);
                        double peakSigma = (lastTick - firstTick) / 3.;  // Set the width...
                        double peakAmp   = 0.3989 * sumADC / peakSigma;  // Use gaussian formulation
                        double peakMean  = (firstTick + lastTick) / 2.;

                        // Store hit params
                        reco_tool::IPeakFitter::PeakFitParams_t peakParams;

                        peakParams.peakCenter         = peakMean;
                        peakParams.peakCenterError    = 0.1 * peakMean;
                        peakParams.peakSigma          = peakSigma;
                        peakParams.peakSigmaError     = 0.1 * peakSigma;
                        peakParams.peakAmplitude      = peakAmp;
                        peakParams.peakAmplitudeError = 0.1 * peakAmp;

                        peakParamsVec.push_back(peakParams);

                        // set for next loop
                        firstTick = lastTick;
                        lastTick  = std::min(lastTick  + longPulseWidth, endT);
                    }
                }

                // #######################################################
                // ### Loop through returned peaks and make recob hits ###
                // #######################################################

                int numHits(0);

                // Make a container for what will be the filtered collection
                std::vector<recob::Hit> filteredHitVec;

                for(const auto& peakParams : peakParamsVec)
                {
                    // Extract values for this hit
                    float peakAmp   = peakParams.peakAmplitude;
                    float peakMean  = peakParams.peakCenter;
                    float peakWidth = peakParams.peakSigma;

                    // Place one bit of protection here
                    if (std::isnan(peakAmp))
                    {
                        std::cout << "**** hit peak amplitude is a nan! Channel: " << channel << ", start tick: " << startT << std::endl;
                        continue;
                    }

                    // Extract errors
                    float peakAmpErr   = peakParams.peakAmplitudeError;
                    float peakMeanErr  = peakParams.peakCenterError;
                    float peakWidthErr = peakParams.peakSigmaError;

                    // ### Charge ###
                    float charge    = chargeFunc(peakMean, peakAmp, peakWidth, fAreaNormsVec[plane],startT,endT);;
                    float chargeErr = std::sqrt(TMath::Pi()) * (peakAmpErr*peakWidthErr + peakWidthErr*peakAmpErr);

                    // ### limits for getting sums
                    std::vector<float>::const_iterator sumStartItr = range.begin() + startT;
                    std::vector<float>::const_iterator sumEndItr   = range.begin() + endT;

                    // ### Sum of ADC counts
                    double sumADC = std::accumulate(sumStartItr, sumEndItr, 0.);

                    // ok, now create the hit
                    recob::HitCreator hitcreator(*wire,                            // wire reference
                                                 wid,                              // wire ID
                                                 startT+roiFirstBinTick,           // start_tick TODO check
                                                 endT+roiFirstBinTick,             // end_tick TODO check
                                                 peakWidth,                        // rms
                                                 peakMean+roiFirstBinTick,         // peak_time
                                                 peakMeanErr,                      // sigma_peak_time
                                                 peakAmp,                          // peak_amplitude
                                                 peakAmpErr,                       // sigma_peak_amplitude
                                                 charge,                           // hit_integral
                                                 chargeErr,                        // hit_sigma_integral
                                                 sumADC,                           // summedADC FIXME
                                                 nGausForFit,                      // multiplicity
                                                 numHits,                          // local_index TODO check that the order is correct
                                                 chi2PerNDF,                       // goodness_of_fit
                                                 NDF                               // dof
                                                 );

                    filteredHitVec.push_back(hitcreator.copy());

                    const recob::Hit hit(hitcreator.move());

                    // This loop will store ALL hits
                    if (fMakeRawDigitAssns) allHitCol.emplace_back(std::move(hit), wire, RawDigits.at(wireIter));
                    else                    allHitCol.emplace_back(std::move(hit), wire);
                    numHits++;
                } // <---End loop over gaussians

                // Should we filter hits?
                if (filteredHitCol && !filteredHitVec.empty())
                {
                    // #######################################################################
                    // Is all this sorting really necessary?  Would it be faster to just loop
                    // through the hits and perform simple cuts on amplitude and width on a
                    // hit-by-hit basis, either here in the module (using fPulseHeightCuts and
                    // fPulseWidthCuts) or in HitFilterAlg?
                    // #######################################################################

                    // Sort in ascending peak height
                    std::sort(filteredHitVec.begin(),filteredHitVec.end(),[](const auto& left, const auto& right){return left.PeakAmplitude() > right.PeakAmplitude();});

                    // Reject if the first hit fails the PH/wid cuts
                    if (filteredHitVec.front().PeakAmplitude() < fPulseHeightCuts.at(plane) || filteredHitVec.front().RMS() < fPulseWidthCuts.at(plane)) filteredHitVec.clear();

                    // Now check other hits in the snippet
                    if (filteredHitVec.size() > 1)
                    {
                        // The largest pulse height will now be at the front...
                        float largestPH = filteredHitVec.front().PeakAmplitude();

                        // Find where the pulse heights drop below threshold
                        float threshold(fPulseRatioCuts.at(plane));

                        std::vector<recob::Hit>::iterator smallHitItr = std::find_if(filteredHitVec.begin(),filteredHitVec.end(),[largestPH,threshold](const auto& hit){return hit.PeakAmplitude() < 8. && hit.PeakAmplitude() / largestPH < threshold;});

                        // Shrink to fit
                        if (smallHitItr != filteredHitVec.end()) filteredHitVec.resize(std::distance(filteredHitVec.begin(),smallHitItr));

                        // Resort in time order
                        std::sort(filteredHitVec.begin(),filteredHitVec.end(),[](const auto& left, const auto& right){return left.PeakTime() < right.PeakTime();});
                    }

                    // Copy the hits we want to keep to the filtered hit collection
                    for(const auto& filteredHit : filteredHitVec)
                        if (!fHitFilterAlg || fHitFilterAlg->IsGoodHit(filteredHit))
                        {
                            if (fMakeRawDigitAssns) filteredHitCol->emplace_back(filteredHit, wire, RawDigits.at(wireIter));
                            else                    filteredHitCol->emplace_back(filteredHit, wire);
                        }
                }

                fChi2->Fill(chi2PerNDF);

            }//<---End loop over merged candidate hits

        } //<---End looping over ROI's

   }//<---End looping over all the wires


    //==================================================================================================
    // End of the event -- move the hit collection and the associations into the event

    if ( filteredHitCol ){

      // If we filtered hits but no instance name was
      // specified for the "all hits" collection, then
      // only save the filtered hits to the event
      if( fAllHitsInstanceName == "" ) {
        filteredHitCol->put_into(evt);

      // otherwise, save both
      } else {
        filteredHitCol->put_into(evt);
        allHitCol.put_into(evt);
      }

    } else {
      allHitCol.put_into(evt);
    }

    // Keep track of events processed
    fEventCount++;

} // End of produce()


  DEFINE_ART_MODULE(GausHitFinder)

} // end of hit namespace
