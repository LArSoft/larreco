////////////////////////////////////////////////////////////////////////
// Class:       HitCheater
// Module Type: producer
// File:        HitCheater_module.cc
//
// Generated at Tue Nov  8 09:41:20 2011 by Brian Rebel using artmod
// from art v1_00_02.
////////////////////////////////////////////////////////////////////////
#ifndef HitCheater_h
#define HitCheater_h

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib/exception.h"

#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "Utilities/DetectorProperties.h"
#include "SimpleTypesAndConstants/geo_types.h"

#include "TH1.h"
#include "TString.h"

class TH1D;

namespace hit {
  class HitCheater;
}

namespace recob {
  class Wire;
  class Hit;
}

namespace sim {
  class IDE;
}

class hit::HitCheater : public art::EDProducer {
public:
  explicit HitCheater(fhicl::ParameterSet const & p);
  virtual ~HitCheater();

  virtual void produce(art::Event & e);

  virtual void beginJob();
  virtual void reconfigure(fhicl::ParameterSet const & p);

private:

  void FindHitsOnChannel(const sim::SimChannel*   sc,
			 std::vector<recob::Hit>& hits,
  			 art::Ptr<recob::Wire>&   wire, 
			 int                      spill);


  std::string         fG4ModuleLabel;              ///< label name for module making sim::SimChannels		
  std::string         fWireModuleLabel;      	   ///< label name for module making recob::Wires		
  double              fMinCharge;            	   ///< Minimum charge required to make a hit                 
  double              fElectronsToADC;             ///< Conversion factor of electrons to ADC counts
  std::string         fCalDataProductInstanceName; ///< label name for module making recob::Wires
  int                 fReadOutWindowSize;          ///< Number of samples in a readout window; NOT total samples
  int                 fNumberTimeSamples;          ///< Number of total time samples (N*readoutwindowsize)
  double              fSamplingRate;               ///< from util::DetectorProperties
  int                 fTriggerOffset;              ///< from util::DetectorProperties
  int                 fNewHitTDCGap;               ///< gap allowed in tdcs without charge before making a new hit
};

//-------------------------------------------------------------------
hit::HitCheater::HitCheater(fhicl::ParameterSet const & p)
{
  this->reconfigure(p);

  produces< std::vector<recob::Hit> >();
}

//-------------------------------------------------------------------
hit::HitCheater::~HitCheater()
{
}

//-------------------------------------------------------------------
void hit::HitCheater::produce(art::Event & e)
{
  // make the unique_ptr for the hits
  std::unique_ptr< std::vector<recob::Hit> > hits(new std::vector<recob::Hit>);

  // Read in the wire List object(s).
  art::Handle< std::vector<recob::Wire> > wHandle;
  int whatSpill = 1;
  if( fCalDataProductInstanceName.size()>0 ) {
    e.getByLabel(fWireModuleLabel,fCalDataProductInstanceName,wHandle);
    if( fCalDataProductInstanceName.find("ost") != std::string::npos) whatSpill=2;
    else whatSpill=0;
  }
  else {
    e.getByLabel(fWireModuleLabel,wHandle);
  }
  
  // make a map of wires to channel numbers
  std::map<unsigned int, art::Ptr<recob::Wire> > wireMap;
  
  for(size_t wc = 0; wc < wHandle->size(); ++wc){
    art::Ptr<recob::Wire> wire(wHandle, wc);
    wireMap[wire->RawDigit()->Channel()] = wire;
  }

  // get the sim::SimChannels out of the event
  std::vector<const sim::SimChannel*> sccol;
  e.getView(fG4ModuleLabel, sccol);

  // find the hits on each channel
  for(size_t sc = 0; sc < sccol.size(); ++sc){

    FindHitsOnChannel(sccol.at(sc), *hits, wireMap.find(sccol.at(sc)->Channel())->second, whatSpill);

  }// end loop over SimChannels

  // put the cheated hits into the event
  LOG_DEBUG("HitCheater") << "putting " << hits->size() << " hits into the event";
  e.put(std::move(hits));

  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::FindHitsOnChannel(const sim::SimChannel*   sc,
					std::vector<recob::Hit>& hits,
					art::Ptr<recob::Wire>&   wire, 
					int                      spill)
{
  art::ServiceHandle<geo::Geometry> geo;

  // determine the possible geo::WireIDs for this particular channel
  // then make a map of tdc to electrons for each one of those geo::WireIDs
  // then find hits on each geo::WireID
  std::vector<geo::WireID> wireids = geo->ChannelToWire(wire->Channel());
  
  std::map<geo::WireID, std::map< unsigned int, double> > wireIDSignals;

  auto const& idemap = sc->TDCIDEMap();

  for(auto const& mapitr : idemap){
    unsigned short tdc = mapitr.first;

    if( fReadOutWindowSize != fNumberTimeSamples ) {
      if( tdc < spill*fReadOutWindowSize || 
	  tdc > (spill+1)*fReadOutWindowSize )  continue;
    }
    
    // figure out which TPC we are in for each voxel
    std::vector<double> pos(3, 0.);
    unsigned int tpc   = 0;
    unsigned int cstat = 0;
    float        edep  = 0.;

    for(auto const& ideitr : mapitr.second){

      edep = ideitr.numElectrons;
      
      pos[0] = ideitr.x;
      pos[1] = ideitr.y;
      pos[2] = ideitr.z;

      try{
	geo->PositionToTPC(&pos[0], tpc, cstat);
      }
      catch(cet::exception &e){
	mf::LogWarning("HitCheater") << "caught exception \n"
				     << e
				     << "when attempting to find TPC for position "
				     << "move on to the next sim::IDE";
	continue;
      }

      for( auto const& wid : wireids){
	if(wid.TPC == tpc && wid.Cryostat == cstat){
	  // in the right TPC, now figure out which wire we want
	  // this works because there is only one plane option for 
	  // each WireID in each TPC
	  if(wid.Wire == geo->NearestWire(pos, wid.Plane, wid.TPC, wid.Cryostat))
	    wireIDSignals[wid][tdc] += edep;
	}// end if in the correct TPC and Cryostat
      }// end loop over wireids for this channel
    }// end loop over ides for this channel
  }// end loop over tdcs for this channel

  // now loop over each wire ID and determine where the hits are
  for( auto widitr : wireIDSignals){

    // get the first tdc in the 
    unsigned short prev         = widitr.second.begin()->first;
    double         startTime    = prev;
    double         totCharge    = 0.;
    double         maxCharge    = -1.;
    double         peakTime     = 0.;
    int            multiplicity =  1 ;

    // loop over all the tdcs for this geo::WireID
    for( auto tdcitr : widitr.second){
      unsigned short tdc = tdcitr.first;
      
      double adc = fElectronsToADC*tdcitr.second;

      // more than a one tdc gap between times with 
      // signal, start a new hit
      if(tdc - prev > fNewHitTDCGap){
	
	if(totCharge > fMinCharge){
	  hits.push_back(recob::Hit(wire, 
				    widitr.first,
				    startTime, 1.,
				    prev,      1.,
				    peakTime,  1.,
				    totCharge, std::sqrt(totCharge),
				    maxCharge, std::sqrt(maxCharge),
				    multiplicity,
				    1.)
			 );
	  
	  LOG_DEBUG("HitCheater") << "new hit is " << hits.back();
	  
	}// end if charge is large enough

	// reset the variables for each hit
	startTime = tdc;
	peakTime  = tdc;
	totCharge = 0.;
	maxCharge = -1.;

      }// end if need to start a new hit

      totCharge += adc;
      if(adc > maxCharge){
	maxCharge = adc;
	peakTime = tdc;
      }

      prev = tdc;

    }// end loop over tdc values for the current geo::WireID


    // We might have missed the last hit, so do it now
    if(totCharge > fMinCharge){
      hits.push_back(recob::Hit(wire, 
				widitr.first,
				startTime, 1.,
				prev,      1.,
				peakTime,  1.,
				totCharge, std::sqrt(totCharge),
				maxCharge, std::sqrt(maxCharge),
				multiplicity,
				1.)
		     );
      
      LOG_DEBUG("HitCheater") << "new hit is " << hits.back();
      
    }// end if charge is large enough
    
  }// end loop over map of geo::WireID to map<tdc,electrons>
  

  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::beginJob()
{
  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::reconfigure(fhicl::ParameterSet const & p)
{
  fG4ModuleLabel   = p.get< std::string >("G4ModuleLabel",   "largeant");
  fWireModuleLabel = p.get< std::string >("WireModuleLabel", "caldata" );
  fMinCharge       = p.get< double      >("MinimumCharge",   5.        );
  fNewHitTDCGap    = p.get< int         >("NewHitTDCGap",    1         );

  art::ServiceHandle<util::DetectorProperties> detprop;
  fElectronsToADC = detprop->ElectronsToADC();
  fSamplingRate   = detprop->SamplingRate();
  fTriggerOffset  = detprop->TriggerOffset();

   fCalDataProductInstanceName="";
   size_t pos = fWireModuleLabel.find(":");
   if( pos!=std::string::npos ) {
     fCalDataProductInstanceName = fWireModuleLabel.substr( pos+1 );
     fWireModuleLabel = fWireModuleLabel.substr( 0, pos );
   }

   fReadOutWindowSize  = detprop->ReadOutWindowSize();
   fNumberTimeSamples  = detprop->NumberTimeSamples();

  return;
}

#endif /* HitCheater_h */

DEFINE_ART_MODULE(hit::HitCheater)
