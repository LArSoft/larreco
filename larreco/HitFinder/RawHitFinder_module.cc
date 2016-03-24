#ifndef RAWHITFINDER_H
#define RAWHITFINDER_H

///////////////////////////////
//Hit finder that runs on raw signals instead of deconvolutes ones //
// Written initially for dune 35t online filter //
///////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <fstream>

//Framework
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Utilities/InputTag.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBaseArt/HitCreator.h"

//LArSoft From FFT

#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"

//ROOT from CalData

#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

//ROOT From Gauss

#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"

namespace hit {

  class RawHitFinder : public art::EDProducer {

    public:

      explicit RawHitFinder(fhicl::ParameterSet const& pset); 
      virtual ~RawHitFinder();

      //    std::vector<float>  Signal() const;
      void produce(art::Event& evt); 
      void beginJob(); 
      void endJob(); 
      void reconfigure(fhicl::ParameterSet const& p);


    private:
      unsigned int  fDataSize;               ///< size of raw data on one wire
      // int        fPostsample;            ///< number of postsample bins
      art::InputTag fDigitModuleLabel;     ///< module that made digits (constants)
      std::string   fSpillName;           ///< nominal spill is an empty string

      //FFT Copied Variables

      std::string         fCalDataModuleLabel;
      std::string         fHitLabelName;
      double              fMinSigInd;            ///<Induction signal height threshold 
      double              fMinSigCol;           ///<Collection signal height threshold 
      double              fIndWidth;           ///<Initial width for induction fit
      double              fColWidth;          ///<Initial width for collection fit
      double              fIndMinWidth;      ///<Minimum induction hit width
      double              fColMinWidth;     ///<Minimum collection hit width
      double              fIncludeMoreTail;
      int                 fMaxMultiHit;    ///<maximum hits for multi fit
      int                 fAreaMethod;    ///<Type of area calculation  
      std::vector<double> fAreaNorms;    ///<factors for converting area to same units as peak height 
      bool                fUncompressWithPed;                       ///< Option to uncompress with pedestal.

    protected: 

  }; // class RawHitFinder


  //-------------------------------------------------
  RawHitFinder::RawHitFinder(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations to raw digits BUT NOT associations to wires
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this, 
        /*instance_name*/"", 
        /*doWireAssns*/false, 
        /*doRawDigitAssns*/true);
  }

  //-------------------------------------------------
  RawHitFinder::~RawHitFinder()
  {
  }

  //////////////////////////////////////////////////////
  void RawHitFinder::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< art::InputTag >("DigitModuleLabel", "daq");
    //  fPostsample       = p.get< int >        ("PostsampleBins");
    // fHitLabelName = p.get< std::string >("HitLabelName", "hit");
    // size_t pos = fDigitModuleLabel.find(":");
    // if( pos!=std::string::npos ) {
    //   fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    // }

    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = p.get< double       >("MinSigInd");
    fMinSigCol          = p.get< double       >("MinSigCol"); 
    fIncludeMoreTail    = p.get< double       >("IncludeMoreTail", 0.);
    fIndWidth           = p.get< double       >("IndWidth");  
    fColWidth           = p.get< double       >("ColWidth");
    fIndMinWidth        = p.get< double       >("IndMinWidth");
    fColMinWidth        = p.get< double       >("ColMinWidth"); 	  	
    fMaxMultiHit        = p.get< int          >("MaxMultiHit");
    fAreaMethod         = p.get< int          >("AreaMethod");
    fAreaNorms          = p.get< std::vector< double > >("AreaNorms");
    fUncompressWithPed  = p.get< bool         >("UncompressWithPed", true);
    mf::LogInfo("RawHitFinder_module") << "fDigitModuleLabel: " << fDigitModuleLabel << std::endl;

  }

  //-------------------------------------------------
  void RawHitFinder::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void RawHitFinder::endJob()
  {  
  }

  void RawHitFinder::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;

    bool retVal = evt.getByLabel(fDigitModuleLabel, digitVecHandle);
    if(retVal==true) 
      mf::LogInfo("RawHitFinder_module") << "I got fDigitModuleLabel: " << fDigitModuleLabel << std::endl;
    else 
      mf::LogWarning("RawHitFinder_module") << "Could not get fDigitModuleLabel: " << fDigitModuleLabel << std::endl;

    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc;                // vector holding uncompressed adc values

    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // this contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator hcol
      (*this, evt, false /* doWireAssns */, true /* doRawDigitAssns */);


    std::vector<float> startTimes;             // stores time of window start
    std::vector<float> maxTimes;    	     // stores time of local maximum    
    std::vector<float> endTimes;    	     // stores time of window end
    std::vector<float> peakHeight;    	     // stores adc counts at maximum
    std::vector<float> hitrms;    	     // stores charge-weighted rms of time across hit
    std::vector<double> charge;         // stores the total charge assoc. with the hit
    //   calculated as the sum of the ADC counts over the window    

    uint32_t channel      = 0;              // channel number

    double threshold       = 0.;             // minimum signal size for id'ing a hit
    double totSig = 0;
    double myrms = 0;
    double mynorm = 0;
    // double fitWidth        = 0.;             // hit fit width initial value
    // double minWidth        = 0.;             // minimum hit width
    geo::SigType_t sigType = geo::kInduction;// type of plane we are looking at
    std::stringstream numConv;

    // loop over all channels    
    hcol.reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();

      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel   = digitVec->Channel();
      fDataSize = digitVec->Samples();
      //      std::cout << " Channel " << channel << std::endl;

      rawadc.resize(fDataSize);
      holder.resize(fDataSize);

      // uncompress the data
      if (fUncompressWithPed){
        int pedestal = (int)digitVec->GetPedestal();
        raw::Uncompress(digitVec->ADCs(), rawadc, pedestal, digitVec->Compression());
      }
      else{
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      }

      for(unsigned int bin = 0; bin < fDataSize; ++bin){ 
        holder[bin]=(rawadc[bin]-digitVec->GetPedestal());
      }
      //holder[bin]=rawadc[bin];
      //now holder and rawadc should be filled correctly

      sigType       = geom->SignalType(channel);

      peakHeight.clear();
      endTimes.clear();
      startTimes.clear();
      maxTimes.clear();
      charge.clear();
      hitrms.clear();


      // ###############################################
      // ###             Induction Planes            ###
      // ###############################################

      if(sigType == geo::kInduction){
        threshold = fMinSigInd;
        //	std::cout<< "Threshold is " << threshold << std::endl;
        // fitWidth = fIndWidth;
        // minWidth = fIndMinWidth;
        //	continue;
        float negthr=-1.0*threshold;
        unsigned int bin =1;
        float minadc=0;

        // find the dips
        while (bin<fDataSize) {  // loop over ticks
          float thisadc = holder[bin];
          if (thisadc<negthr) { // new region
            //	    std::cout << "new region" << bin << " " << thisadc << std::endl;
            // step back to find zero crossing
            unsigned int place = bin;
            while (thisadc<=0 && bin>0) {
              //	      std::cout << bin << " " << thisadc << std::endl;
              bin--;
              thisadc=holder[bin];
            }
            float hittime = bin+thisadc/(thisadc-holder[bin+1]);
            maxTimes.push_back(hittime);

            // step back more to find the hit start time
            while (thisadc<threshold && bin>0) {
              //	      std::cout << bin << " " << thisadc << std::endl;
              bin--;
              thisadc=holder[bin];
            }
            if (bin>=2) bin-=2;
            while (thisadc>threshold && bin>0) {

              //	        std::cout << bin << " " << thisadc << std::endl;
              bin--;
              thisadc=holder[bin];
            }
            startTimes.push_back(bin+1);
            // now step forward from hit time to find end time
            bin=place; 	      
            thisadc=holder[bin];
            minadc=thisadc;

            totSig = 0;
            while (thisadc<negthr && bin<fDataSize) {
              //	        std::cout << bin << " " << thisadc << std::endl;
              totSig += fabs(thisadc); 
              bin++;
              thisadc=holder[bin];
              //	          std::cout << "ADC VALUE INDUCTION" << thisadc << std::endl;
              if (thisadc<minadc) minadc=thisadc;		
            }
            endTimes.push_back(bin-1);
            peakHeight.push_back(-1.0*minadc);
            charge.push_back(totSig);
            hitrms.push_back(5.0);
            //	    std::cout << "TOTAL SIGNAL INDUCTION " << totSig << "  5.0" << std::endl; 
            // std::cout << "filled end times " << bin-1 << "peak height vector size " << peakHeight.size() << std::endl;

            // don't look for a new hit until it returns to baseline
            while (thisadc<0 && bin<fDataSize) {
              //	      std::cout << bin << " " << thisadc << std::endl;
              bin++;
              thisadc=holder[bin];
            }
          } // end region
          bin++;	  
        }// loop over ticks
      }

      // ###############################################
      // ###             Collection Plane            ###
      // ###############################################    

      else if(sigType == geo::kCollection){
        threshold = fMinSigCol;
        // fitWidth  = fColWidth;
        // minWidth  = fColMinWidth;

        //find local maxima
        float madc=threshold;
        int ibin=0;
        //	int nohits = 0;
        unsigned int bin =0;
        int start = 0;
        int end = 0;

        while (bin<fDataSize) {
          float thisadc = holder[bin];
          madc=threshold;

          if (thisadc>madc) { // new region
            startTimes.push_back(bin);
            start = bin;
            bin++;

            while (thisadc>threshold && bin<fDataSize) {
              thisadc=holder[bin];
              if (thisadc>madc) {ibin=bin; madc=thisadc;}
              bin++;
            }
            maxTimes.push_back(ibin-1);
            peakHeight.push_back(madc);	    
            endTimes.push_back(bin-1);
            end = bin-1;

            totSig = 0;
            myrms  = 0;
            mynorm = 0;

            for(int i = start-std::ceil(fIncludeMoreTail*(end-start)); i <= end+std::ceil(fIncludeMoreTail*(end-start)); i++){
              totSig += holder[i];
              float temp2 = holder[i]*holder[i];
              mynorm += temp2;
              float temp = ibin-1.0-i;
              myrms += temp*temp*temp2;
            }
            charge.push_back(totSig);
            myrms/=mynorm;
            if((end-start+2*std::ceil(fIncludeMoreTail*(end-start)+1))!=0)
            {
              myrms/=(float)(end-start+2*std::ceil(fIncludeMoreTail*(end-start)+1));
              hitrms.push_back(sqrt(myrms));
            }
            else
            {
              hitrms.push_back(sqrt(myrms));
            }
            //   std::cout << "CHARGE ON ADC####################### " << totSig 
            // << " RMS " << myrms << std::endl;
            //  nohits++;
          }// end region
          //nohits = 0;

          start = 0;
          end = 0;
          bin++;
        }
      }
      //      std::cout << "channel " << channel << "hits found  " <<  maxTimes.size() << std::endl;

      int numHits(0);   // number of consecutive hits being fitted
      int hitIndex(0);  // index of current hit in sequence
      //      double amplitude(0), position(0), width(0);  //fit parameters
      double amplitude(0), position(0);  //fit parameters
      double start(0), end(0);
      double amplitudeErr(0), positionErr(0);  //fit errors
      double goodnessOfFit(0), chargeErr(0);  //Chi2/NDF and error on charge
      double hrms(0);

      numHits = maxTimes.size();
      for (int i=0;i<numHits;++i) {
        //	int index = int(maxTimes[i]+0.5);
        //	 std::cout << " Channel " << channel << " Hit " << i+1 << " Max Time " << maxTimes[i] << " Start//  Time " <<
        // startTimes[i] << " End Time " << endTimes[i] << " Peak Height " << peakHeight[i] << " Charge " << charge[i] << "  RMS  "  << hitrms[i] << std::endl;

        //	amplitude     = holder[index];

        amplitude     = peakHeight[i];
        position      = maxTimes[i];
        start         = startTimes[i];
        end           = endTimes[i];
        hrms = hitrms[i];
        amplitudeErr  = -1;
        positionErr   = 1.0;
        goodnessOfFit = -1;
        chargeErr = -1;
        totSig = charge[i];

        // get the WireID for this hit
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        geo::WireID wid = wids[0];

        // make the hit
        if (start>=end){
          mf::LogWarning("RawHitFinder_module") << "Hit start "<<start<<" is >= hit end "<<end;
          continue;
        }

        recob::HitCreator hit(
            *digitVec,        // raw digit reference
            wid,              // wire ID
            start,            // start_tick FIXME
            end,              // end_tick FIXME
            hrms,              // rms FIXME
            position,         // peak_time
            positionErr,      // sigma_peak_time
            amplitude,        // peak_amplitude
            amplitudeErr,     // sigma_peak_amplitude
            totSig,           // hit_integral
            chargeErr,        // hit_sigma_integral
            std::accumulate   // summedADC FIXME
            (holder.begin() + (int) start, holder.begin() + (int) end, 0.), 
            1,                // multiplicity FIXME
            -1,               // local_index FIXME
            goodnessOfFit,    // goodness_of_fit
            int(end - start)  // dof
            );
        hcol.emplace_back(hit.move(), digitVec);

        /*
        if(evt.event()==21&&channel==2022)
        {
          std::cout << position << std::endl;
        }*/

        ++hitIndex;
      }//end loop over hits
      //hitIndex += numHits;	
    } // end loop over channels

    //    std::cerr << "hcol.size(): " << hcol.size() << std::endl;//jpd

    // std::cerr << "I produced fHitLabelName: " << fHitLabelName << std::endl;


    // move the hit collection and the associations into the event

    hcol.put_into(evt);
  }

  DEFINE_ART_MODULE(RawHitFinder)   

} // end namespace hit


#endif //RAWHITFINDER_H
