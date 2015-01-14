#ifndef RAWHITFINDER_H
#define RAWHITFINDER_H

///////////////////////////////
//Hit finder that runs on raw signals instead of deconvolutes ones //
// Written initially for lbne 35t online filter //
///////////////////////////////

#include <string>
#include <vector>
#include <stdint.h>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
//Framework
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

//LArSoft

#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "Utilities/LArFFT.h"
#include "RecoBase/Hit.h"

//LArSoft From FFT

#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"

//ROOT from CalData

#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

//ROOT From Gauss

#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"

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
    
    unsigned int          fDataSize;               ///< size of raw data on one wire
    // int          fPostsample;            ///< number of postsample bins
    std::string  fDigitModuleLabel;     ///< module that made digits (constants)
    std::string  fSpillName;           ///< nominal spill is an empty string
    ///< it is set by the DigitModuleLabel
    ///< ex.:  "daq:preSpill" for prespill data

    unsigned int            fMaxSamples; ///< max number of ADC samples possible on the wire
    
    //FFT Copied Variables
    
    std::string         fCalDataModuleLabel;
    std::string         fHitLabelName;
    double              fMinSigInd;            ///<Induction signal height threshold 
    double              fMinSigCol;           ///<Collection signal height threshold 
    double              fIndWidth;           ///<Initial width for induction fit
    double              fColWidth;          ///<Initial width for collection fit
    double              fIndMinWidth;      ///<Minimum induction hit width
    double              fColMinWidth;     ///<Minimum collection hit width
    int                 fMaxMultiHit;    ///<maximum hits for multi fit
    int                 fAreaMethod;    ///<Type of area calculation  
    std::vector<double> fAreaNorms;    ///<factors for converting area to same units as peak height 

  protected: 
    
  }; // class RawHitFinder

  
  //-------------------------------------------------
  RawHitFinder::RawHitFinder(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    //    produces< std::vector<recob::Hit> > (fHitLabelName);
    produces< std::vector<recob::Hit> > ();
  }
  
  //-------------------------------------------------
  RawHitFinder::~RawHitFinder()
  {
  }

  //////////////////////////////////////////////////////
  void RawHitFinder::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    //  fPostsample       = p.get< int >        ("PostsampleBins");
    // fHitLabelName = p.get< std::string >("HitLabelName", "hit");
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
 
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = p.get< double       >("MinSigInd");
    fMinSigCol          = p.get< double       >("MinSigCol"); 
    fIndWidth           = p.get< double       >("IndWidth");  
    fColWidth           = p.get< double       >("ColWidth");
    fIndMinWidth        = p.get< double       >("IndMinWidth");
    fColMinWidth        = p.get< double       >("ColMinWidth"); 	  	
    fMaxMultiHit        = p.get< int          >("MaxMultiHit");
    fAreaMethod         = p.get< int          >("AreaMethod");
    fAreaNorms          = p.get< std::vector< double > >("AreaNorms");
    
    //    std::cerr << "fDigitModuleLabel: " << fDigitModuleLabel << std::endl;//jpd

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
    if(retVal==true) std::cerr << "I got fDigitModuleLabel: " << fDigitModuleLabel << std::endl;
    else std::cerr << "I did not get fDigitModuleLabel: " << fDigitModuleLabel << std::endl;

    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc;                // vector holding uncompressed adc values

    // Making a ptr vector to put on the event
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);

    
    std::vector<float> startTimes;             // stores time of 1st local minimum
    std::vector<float> maxTimes;    	     // stores time of local maximum    
    std::vector<float> endTimes;    	     // stores time of 2nd local minimum
    std::vector<float> peakHeight;    	     // stores time of 2nd local minimum
    uint32_t channel      = 0;              // channel number


    double threshold       = 0.;             // minimum signal size for id'ing a hit
    // double fitWidth        = 0.;             // hit fit width initial value
    // double minWidth        = 0.;             // minimum hit width
    geo::SigType_t sigType = geo::kInduction;// type of plane we are looking at
    std::stringstream numConv;

 
    // loop over all channels    
    hcol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
	// get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel =  digitVec->Channel();
      fMaxSamples = digitVec->Samples();
      fDataSize = digitVec->fADC.size();
      //      std::cout << " Channel " << channel << std::endl;


      // uncompress the data
      raw::Uncompress(digitVec->fADC, rawadc, digitVec->Compression());

      fDataSize = rawadc.size();
      rawadc.resize(fDataSize);
      holder.resize(fDataSize);

       for(unsigned int bin = 0; bin < fDataSize; ++bin) 
      	holder[bin]=(rawadc[bin]-digitVec->GetPedestal());
      //now holder and rawadc should be filled correctly
      
      channel       = digitVec->Channel();
      sigType       = geom->SignalType(channel);
	
      peakHeight.clear();
      endTimes.clear();
      startTimes.clear();
      maxTimes.clear();


      //Set the appropriate signal widths and thresholds
      if(sigType == geo::kInduction){
	threshold = fMinSigInd;
	// fitWidth = fIndWidth;
	// minWidth = fIndMinWidth;
	//	continue;
	float negthr=-1.0*threshold;
	unsigned int bin =0;
	float minadc=0;
      	
	// find the dips
	while (bin<fDataSize) {  // loop over ticks
	  float thisadc = holder[bin];
	  if (thisadc<negthr) { // new region
	    //std::cout << "new region" << bin << " " << thisadc << std::endl;
	    // step back to find zero crossing
	    unsigned int place = bin;
	    while (thisadc<0 && bin>0) {
	      //std::cout << bin << " " << thisadc << std::endl;
	      bin--;
	      thisadc=holder[bin];
	    }
	    float hittime = bin+thisadc/(thisadc-holder[bin+1]);
	    maxTimes.push_back(hittime);
	    // step back more to find the hit start time
	    while (thisadc<threshold && bin>0) {
	      //std::cout << bin << " " << thisadc << std::endl;
	      bin--;
	      thisadc=holder[bin];
	    }
	    if (bin>=2) bin-=2;
	    //	    std::cout << bin << " " << thisadc << std::endl;	    bin--;  
	    //	    std::cout << bin << " " << thisadc << std::endl;	    bin--;  
	    while (thisadc>threshold && bin>0) {
	      //std::cout << bin << " " << thisadc << std::endl;
	      bin--;
	      thisadc=holder[bin];
	    }
	    startTimes.push_back(bin+1);
	    // now step forward from hit time to find end time
	    bin=place; 	      
	    thisadc=holder[bin];
	    minadc=thisadc;
	    while (thisadc<negthr && bin<fDataSize) {
	      //	      std::cout << bin << " " << thisadc << std::endl;
	      bin++;
	      thisadc=holder[bin];
	      if (thisadc<minadc) minadc=thisadc;		
	    }
	    endTimes.push_back(bin-1);
	    peakHeight.push_back(-1.0*minadc);
	    while (thisadc<0 && bin<fDataSize) {
	      //	      std::cout << bin << " " << thisadc << std::endl;
	      bin++;
	    }
	  }// end region
	  bin++;	  
	}// loop over ticks
      }
      else if(sigType == geo::kCollection){
	threshold = fMinSigCol;
	// fitWidth  = fColWidth;
	// minWidth  = fColMinWidth;
	
	//find local maxima
	float madc=threshold;
	int ibin=0;
	unsigned int bin =0;
      	
	while (bin<fDataSize) {
	  float thisadc = holder[bin];
	  madc=threshold;
	  if (thisadc>madc) { // new region
	    startTimes.push_back(bin);
	    while (thisadc>madc && bin<fDataSize) {
	      madc=thisadc; ibin=bin;
	      bin++;
	      thisadc=holder[bin];
	    }
	    maxTimes.push_back(ibin-1);
	    peakHeight.push_back(madc);
	    while (thisadc>threshold && bin<fDataSize) {
	      bin++;
	      thisadc=holder[bin];
	    }
	    endTimes.push_back(bin-1);
	  }// end region
	  bin++;
	}
      }
      //      std::cout << "channel " << channel << "hits found  " <<  maxTimes.size() << std::endl;
      double totSig(0); // stores the total hit signal
      int numHits(0);   // number of consecutive hits being fitted
      int hitIndex(0);  // index of current hit in sequence
      //      double amplitude(0), position(0), width(0);  //fit parameters
      double amplitude(0), position(0);  //fit parameters
      double start(0), end(0);
      double amplitudeErr(0), positionErr(0), widthErr(0);  //fit errors
      double goodnessOfFit(0), chargeErr(0);  //Chi2/NDF and error on charge
      std::vector<double>  hitSig;
      
      numHits = maxTimes.size();
      for (int i=0;i<numHits;++i) {
	//	int index = int(maxTimes[i]+0.5);
	// std::cout << channel << " " << i+1 << " " << maxTimes[i] << " " <<
	//   startTimes[i] << " " << endTimes[i] 
	// 	  << " " << peakHeight[i] << std::endl;
	
	//	amplitude     = holder[index];
	amplitude     = peakHeight[i];
	position      = maxTimes[i];
	start         = startTimes[i];
	end         = endTimes[i];
	amplitudeErr  = -1;
	positionErr   = 1.0;
	widthErr      = 0.5;
	goodnessOfFit = -1;
	chargeErr = -1;
	
	// hitSig.resize(3);
	
	// for(int sigPos = 0; sigPos < 3; ++sigPos){
	//   hitSig[sigPos] = holder[position-1+sigPos];
	// }              	    
	
	// totSig = -1;
	
	// get the WireID for this hit
	std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
	geo::WireID wid = wids[0];
	geo::PlaneID plid(wid.Cryostat,wid.TPC,wid.Plane);
	geo::View_t view = geom->View(plid);	      
	geo::SigType_t st = geom->SignalType(plid);
	
	// make the hit
	
	art::Ptr<raw::RawDigit> dv = digitVec;
	recob::Hit hit(dv,
		       view,st,wid,
		       start,
		       widthErr,
		       end,
		       widthErr,
		       position,
		       positionErr,
		       totSig,         
		       chargeErr,                
		       amplitude,
		       amplitudeErr,
		       1,                  /// \todo - mulitplicity has to be determined
		       goodnessOfFit);               	    
	hcol->push_back(hit);
	
      }//end loop over hits
      hitIndex += numHits;	
    } // end loop over channels
    

    //    std::cerr << "hcol->size(): " << hcol->size() << std::endl;//jpd

    // std::cerr << "I produced fHitLabelName: " << fHitLabelName << std::endl;


    //    evt.put(std::move(hcol), fHitLabelName);
    evt.put(std::move(hcol));

  }

  DEFINE_ART_MODULE(RawHitFinder)   

} // end namespace hit


#endif //RAWHITFINDER_H
