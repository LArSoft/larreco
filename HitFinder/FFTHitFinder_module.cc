#ifndef FFTHITFINDER_H
#define FFTHITFINDER_H
////////////////////////////////////////////////////////////////////////
//
// FFTHitFinder class
//
// pagebri3@msu.edu
//
//  This algorithm is designed to find hits on wires after deconvolution
//  with an average shape used as the input response.
////////////////////////////////////////////////////////////////////////
#include <string>
#include <stdint.h>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Core/EDProducer.h" 

// LArSoft Includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Hit.h"

// ROOT Includes 
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"

namespace hit{

  class FFTHitFinder : public art::EDProducer {
    
  public:
    
    explicit FFTHitFinder(fhicl::ParameterSet const& pset); 
    virtual ~FFTHitFinder();
         
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                

  private:
        
    std::string     fCalDataModuleLabel;
    double          fMinSigInd;     ///<Induction signal height threshold 
    double          fMinSigCol;     ///<Collection signal height threshold 
    double          fIndWidth;      ///<Initial width for induction fit
    double          fColWidth;      ///<Initial width for collection fit
    double          fIndMinWidth;   ///<Minimum induction hit width
    double          fColMinWidth;   ///<Minimum collection hit width
    int             fMaxMultiHit;   ///<maximum hits for multi fit
    int             fAreaMethod;    ///<Type of area calculation  
    std::vector<double> fAreaNorms; ///<factors for converting area to same units as peak height 
  protected: 
    
  
  }; // class FFTHitFinder  
  
  //-------------------------------------------------
  FFTHitFinder::FFTHitFinder(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Hit> >();
  }


  //-------------------------------------------------
  FFTHitFinder::~FFTHitFinder()
  {
  }
  
  //-------------------------------------------------
  void FFTHitFinder::reconfigure(fhicl::ParameterSet const& p)
  {
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
  }  

  //-------------------------------------------------
  void FFTHitFinder::beginJob()
  { 
  }

  //-------------------------------------------------
  void FFTHitFinder::endJob()
  {
  }

  //  This algorithm uses the fact that deconvolved signals are very smooth 
  //  and looks for hits as areas between local minima that have signal above 
  //  threshold.
  //-------------------------------------------------
  void FFTHitFinder::produce(art::Event& evt)
  { 
    
    std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
    // Read in the wire List object(s).
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
    art::ServiceHandle<geo::Geometry> geom;
   
    std::vector<int> startTimes;             // stores time of 1st local minimum
    std::vector<int> maxTimes;    	     // stores time of local maximum    
    std::vector<int> endTimes;    	     // stores time of 2nd local minimum
    int time               = 0;              // current time bin
    int minTimeHolder      = 0;              // current start time
    uint32_t channel       = 0;              // channel number
    bool maxFound          = false;          // Flag for whether a peak > threshold has been found
    double threshold       = 0.;             // minimum signal size for id'ing a hit
    double fitWidth        = 0.;             // hit fit width initial value
    double minWidth        = 0.;             // minimum hit width
    std::string eqn        = "gaus(0)";      // string for equation to fit
    geo::SigType_t sigType = geo::kInduction;// type of plane we are looking at
    std::stringstream numConv;


    //loop over wires
    for(unsigned int wireIter = 0; wireIter < wireVecHandle->size(); wireIter++) {

      art::Ptr<recob::Wire> wire(wireVecHandle, wireIter);
      startTimes.clear();
      maxTimes.clear();
      endTimes.clear();
      std::vector<float> signal(wire->Signal());
      std::vector<float>::iterator timeIter;   // iterator for time bins
      time          = 0;
      minTimeHolder = 0;
      maxFound      = false;
      channel       = wire->RawDigit()->Channel();
      sigType       = geom->SignalType(channel);

      //Set the appropriate signal widths and thresholds
      if(sigType == geo::kInduction){
	threshold     = fMinSigInd;
	fitWidth      = fIndWidth;
	minWidth      = fIndMinWidth;
      }
      else if(sigType == geo::kCollection){
	threshold = fMinSigCol;
	fitWidth  = fColWidth;
	minWidth  = fColMinWidth;
      }
      // loop over signal
      for(timeIter = signal.begin(); timeIter+2 < signal.end(); timeIter++){    
	//test if timeIter+1 is a local minimum
	if(*timeIter > *(timeIter+1) && *(timeIter+1) < *(timeIter+2)){
	  //only add points if already found a local max above threshold.
	  if(maxFound) {
	    endTimes.push_back(time+1);
	    maxFound = false;
	    //keep these in case new hit starts right away
	    minTimeHolder = time+2; 
	  }
	  else minTimeHolder = time+1; 
	}
	//if not a minimum, test if we are at a local maximum 
	//if so, and the max value is above threshold, add it and proceed.
	else if(*timeIter < *(timeIter+1) && 
		*(timeIter+1) > *(timeIter+2) && 
		*(timeIter+1) > threshold){ 
	  maxFound = true;
	  maxTimes.push_back(time+1);
	  startTimes.push_back(minTimeHolder);          
	}
	time++;
      }//end loop over signal vec

      
      //if no inflection found before end, but peak found add end point
      while(maxTimes.size()>endTimes.size()) 
	endTimes.push_back(signal.size()-1); 
      if(startTimes.size() == 0) continue;
      
      //All code below does the fitting, adding of hits
      //to the hit vector and when all wires are complete 
      //saving them 
      double totSig(0); // stores the total hit signal
      double startT(0); // stores the start time
      double endT(0);   // stores the end time
      int numHits(0);   // number of consecutive hits being fitted
      int size(0);      // size of data vector for fit
      int hitIndex(0);  // index of current hit in sequence
      double amplitude(0), position(0), width(0);  //fit parameters
      double amplitudeErr(0), positionErr(0), widthErr(0);  //fit errors
      double goodnessOfFit(0), chargeErr(0);  //Chi2/NDF and error on charge
      double minPeakHeight(0);  //lowest peak height in multi-hit fit
     
      //stores gaussian paramters first index is the hit number
      //the second refers to height, position, and width respectively
      std::vector<double>  hitSig;

      //add found hits to hit vector      
      while(hitIndex < (signed)startTimes.size()) {
	
	startT = endT = 0;
	numHits = 1;
        minPeakHeight = signal[maxTimes[hitIndex]];

	//consider adding pulse to group of consecutive hits if:
        //1 less than max consecutive hits
        //2 we are not at the last point in the signal vector
        //3 the height of the dip between the two is greater than threshold/2
        //4 and there is no gap between them
        while(numHits < fMaxMultiHit &&
	      numHits+hitIndex < (signed)endTimes.size() && 
	      signal[endTimes[hitIndex+numHits-1]] >threshold/2.0 &&  
	      startTimes[hitIndex+numHits] - endTimes[hitIndex+numHits-1] < 2){

	  if(signal[maxTimes[hitIndex+numHits]] < minPeakHeight) 
	    minPeakHeight = signal[maxTimes[hitIndex+numHits]];

	  ++numHits;
	}

	//finds the first point > 1/2 the smallest peak
	startT = startTimes[hitIndex];

	while(signal[(int)startT] < minPeakHeight/2.0) ++startT;

	//finds the first point from the end > 1/2 the smallest peak
	endT = endTimes[hitIndex+numHits-1];

	while(signal[(int)endT] <minPeakHeight/2.0) --endT;
	size = (int)(endT-startT);
	TH1D hitSignal("hitSignal","",size,startT,endT);
	for(int i = (int)startT; i < (int)endT; ++i)
	  hitSignal.Fill(i,signal[i]);

        //build the TFormula
        eqn = "gaus(0)";	

	for(int i = 3; i < numHits*3; i+=3) {
	  eqn.append("+gaus(");
	  numConv.str("");
	  numConv << i;
	  eqn.append(numConv.str());
	  eqn.append(")");
	}

	TF1 gSum("gSum",eqn.c_str(),0,size);

	if(numHits > 1) {
	  TArrayD data(numHits*numHits);
	  TVectorD amps(numHits); 
	  for(int i = 0; i < numHits; ++i) {
	    amps[i] = signal[maxTimes[hitIndex+i]];
	    for(int j = 0; j < numHits;j++) 
	      data[i+numHits*j] = TMath::Gaus(maxTimes[hitIndex+j],
					      maxTimes[hitIndex+i],
					      fitWidth);
	  }//end loop over hits
	  
          //This section uses a linear approximation in order to get an
	  //initial value of the individual hit amplitudes 
	  try{
	    TMatrixD h(numHits,numHits);
	    h.Use(numHits,numHits,data.GetArray());
	    TDecompSVD a(h);
	    a.Solve(amps);
	  }
	  catch(...){
	    mf::LogInfo("FFTHitFinder")<<"TDcompSVD failed";
	    hitIndex += numHits;
	    continue;
	  }
      
	  for(int i = 0; i < numHits; ++i) {
	    //if the approximation makes a peak vanish
            //set initial height as average of threshold and
            //raw peak height
            if(amps[i] > 0 ) amplitude = amps[i];
            else amplitude = 0.5*(threshold+signal[maxTimes[hitIndex+i]]);
            gSum.SetParameter(3*i,amplitude);
	    gSum.SetParameter(1+3*i, maxTimes[hitIndex+i]);
	    gSum.SetParameter(2+3*i, fitWidth);
	    gSum.SetParLimits(3*i, 0.0, 3.0*amplitude);
	    gSum.SetParLimits(1+3*i, startT , endT);
	    gSum.SetParLimits(2+3*i, 0.0, 10.0*fitWidth);
	  }//end loop over hits
	}//end if numHits > 1
	else {
	  gSum.SetParameters(signal[maxTimes[hitIndex]],maxTimes[hitIndex],fitWidth);
	  gSum.SetParLimits(0,0.0,1.5*signal[maxTimes[hitIndex]]);
	  gSum.SetParLimits(1, startT , endT);
	  gSum.SetParLimits(2,0.0,10.0*fitWidth);
	}

	/// \todo - just get the integral from the fit for totSig
        hitSignal.Fit(&gSum,"QNRW","", startT, endT);
	for(int hitNumber = 0; hitNumber < numHits; ++hitNumber) {
          totSig = 0;
	  if(gSum.GetParameter(3*hitNumber)   > threshold/2.0 && 
	     gSum.GetParameter(3*hitNumber+2) > minWidth) { 
	    amplitude     = gSum.GetParameter(3*hitNumber);
	    position      = gSum.GetParameter(3*hitNumber+1);
	    width         = gSum.GetParameter(3*hitNumber+2);
            amplitudeErr  = gSum.GetParError(3*hitNumber);
	    positionErr   = gSum.GetParError(3*hitNumber+1);
	    widthErr      = gSum.GetParError(3*hitNumber+2);
            goodnessOfFit = gSum.GetChisquare()/(double)gSum.GetNDF();

	    //estimate error from area of Gaussian
            chargeErr = std::sqrt(TMath::Pi())*(amplitudeErr*width+widthErr*amplitude);   

	    hitSig.resize(size);

	    for(int sigPos = 0; sigPos < size; ++sigPos){
	      hitSig[sigPos] = amplitude*TMath::Gaus(sigPos+startT,position, width);
	      totSig += hitSig[(int)sigPos];              
	    }              	    

            if(fAreaMethod) 
              totSig = std::sqrt(2*TMath::Pi())*amplitude*width/fAreaNorms[(size_t)sigType];              

	    // get the WireID for this hit
	    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
	    ///\todo need to have a disambiguation algorithm somewhere in here
	    // for now, just take the first option returned from ChannelToWire
	    geo::WireID wid = wids[0];

	    // make the hit
	    recob::Hit hit(wire, 
			   wid,
			   position - width, 
                           widthErr,
			   position + width, 
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
	    
	  }//end if over threshold
	}//end loop over hits
	hitIndex += numHits;	
      } // end while on hitIndex<(signed)startTimes.size()

    } // while on Wires
    
    evt.put(std::move(hcol));

  } // End of produce()  
  

  
  DEFINE_ART_MODULE(FFTHitFinder)

} // end of hit namespace


#endif // FFTHITFINDER_H
