////////////////////////////////////////////////////////////////////////
// Class:       DisambigCheater
// Module Type: producer
// File:        DisambigCheater_module.cc
//
// tylerdalion@gmail.com
//
////////////////////////////////////////////////////////////////////////
#ifndef DisambigCheater_h
#define DisambigCheater_h

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Simulation/SimChannel.h"
#include "RecoBase/Hit.h"
#include "Utilities/DetectorProperties.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "MCCheater/BackTracker.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "TH1.h"
#include "TString.h"


namespace recob {
  class Hit;
}

namespace sim {
  class IDE;
}

namespace hit{

  class DisambigCheater : public art::EDProducer {
  public:
    explicit DisambigCheater(fhicl::ParameterSet const & p);
    virtual ~DisambigCheater();
    
    void produce(art::Event & e);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const & p);
    
  private:
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<cheat::BackTracker> bt;
    art::ServiceHandle<art::TFileService> tfs;
    
    std::string fChanHitLabel;
    std::string fWidHitLabel;
    
    std::map< std::pair<double,double>, std::vector<geo::WireID> > fHitToWids;
    void InitHitToWids( const std::vector< art::Ptr<recob::Hit> >& ChHits );

    void MakeDisambigHit(art::Ptr<recob::Hit> hit, geo::WireID wid, 
			 std::unique_ptr< std::vector<recob::Hit> >& hcol);
    
    unsigned int fFalseChanHits = 0; // hits with no IDEs - could be noise or other technical problems
    unsigned int fBadIDENearestWire = 0; // Just to count the inconsistent IDE wire returns
    
  };
  
  //-------------------------------------------------------------------
  DisambigCheater::DisambigCheater(fhicl::ParameterSet const & p)
  {
    this->reconfigure(p);
    produces< std::vector<recob::Hit> >();
  }
  
  //-------------------------------------------------------------------
  hit::DisambigCheater::~DisambigCheater()
  {
  }
  
  //-------------------------------------------------------------------
  void DisambigCheater::produce(art::Event & evt)
  {
    
    // make the unique_ptr for the hits
    std::unique_ptr< std::vector<recob::Hit> > hits(new std::vector<recob::Hit>);
    
    // get hits on channels
    art::Handle< std::vector<recob::Hit> > ChanHits;
    evt.getByLabel(fChanHitLabel, ChanHits);
    std::vector< art::Ptr<recob::Hit> > ChHits;
    art::fill_ptr_vector(ChHits, ChanHits);
    
    // find the wireIDs each hit is on
    this->InitHitToWids( ChHits );
    
    
    // make all of the hits
    for(size_t h=0; h<ChHits.size(); h++){
      
      // the trivial Z hits
      if(ChHits[h]->View()==geo::kZ){
	this->MakeDisambigHit(ChHits[h], ChHits[h]->WireID(), hits);
	continue;
      }
      
      // make U/V hits if any wire IDs are associated
      // count hits without a wireID
      /// \todo: Decide how to handle the hits with multiple-wid activity. For now, randomly choose.
      std::pair<double,double> ChanTime(ChHits[h]->Channel()*1., ChHits[h]->PeakTime()*1.);
      if( fHitToWids[ChanTime].size() == 1 ) 
	this->MakeDisambigHit(ChHits[h], fHitToWids[ChanTime][0], hits);
      else if( fHitToWids[ChanTime].size() == 0 ) 
	fFalseChanHits++;
      else if( fHitToWids[ChanTime].size() > 1 ) 
	this->MakeDisambigHit(ChHits[h], fHitToWids[ChanTime][0], hits); // same thing for now
      
    }
    
    
    evt.put(std::move(hits));
    fHitToWids.clear();
    return;
    
  }
  
  
  //-------------------------------------------------
  void DisambigCheater::MakeDisambigHit(art::Ptr<recob::Hit> hit, geo::WireID wid, 
					std::unique_ptr< std::vector<recob::Hit> >& hcol)
  {
    
    if( !wid.isValid ){ 
      mf::LogWarning("InvalidWireID") << "wid is invalid, hit not being made\n";
      return; }
    
    art::Ptr<recob::Wire> wire = hit->Wire();
    recob::Hit WidHit( wire, 			     wid,
		       hit->StartTime(),     	     hit->SigmaStartTime(),
		       hit->EndTime(),       	     hit->SigmaEndTime(),
		       hit->PeakTime(),      	     hit->SigmaPeakTime(),
		       hit->Charge(),        	     hit->SigmaCharge(),                
		       hit->Charge(true),   	     hit->SigmaCharge(true),
		       hit->Multiplicity(),  	     hit->GoodnessOfFit() );     
    
    hcol->push_back(WidHit);
    return;
    
  }
  
  
  //-------------------------------------------------------------------
  void DisambigCheater::InitHitToWids( const std::vector< art::Ptr<recob::Hit> >& ChHits )
  {
    
    unsigned int Ucount(0), Vcount(0);
    
    for( size_t h = 0; h < ChHits.size(); h++ ){
      if( ChHits[h]->View() == geo::kZ ) continue;
      if( ChHits[h]->View() == geo::kU ) Ucount++;
      if( ChHits[h]->View() == geo::kV ) Vcount++;
      art::Ptr<recob::Hit> chit = ChHits[h];
      std::vector<geo::WireID> cwids = geom->ChannelToWire(chit->Channel());
      std::pair<double,double> ChanTime( chit->Channel()*1., chit->PeakTime()*1.); // hit key value 
      
      // get hit IDEs
      std::vector< sim::IDE > ides;
      bt->HitToSimIDEs( chit, ides );
      
      // catch the hits that have no IDEs
      std::vector<double> xyzpos;
      try{	xyzpos = bt->SimIDEsToXYZ(ides);      }
      catch(...){
	mf::LogVerbatim("DisambigCheat") << "Hit on channel " << chit->Channel()
					 << " between " << chit->StartTime()
					 << " and " << chit->EndTime() << " matches no IDEs.";
	std::vector< geo::WireID > emptyWids;
	fHitToWids[ChanTime] = emptyWids;
	  continue;
      }
      
      // see what wire ID(s) the hit-IDEs are on
      //  if none: hit maps to empty vector
      //  if one: we have a unique hit, vector of size 1
      //  if more than one: make a vector of all wids
      std::vector<geo::WireID> widsWithIdes; 
      for( size_t i=0; i<ides.size(); i++ ){
	double xyzIde[] = { ides[i].x, ides[i].y, ides[i].z };
	unsigned int tpc, cryo;

	// Occasionally, ide position is not in TPC
	/// \todo: Why would an IDE xyz position be outside of a TPC?
	try{	geom->PositionToTPC( xyzIde, tpc, cryo );     }
	catch(...){   mf::LogWarning("DisambigCheat") << "IDE at x = " << xyzIde[0] 
						      << ", y = " << xyzIde[1]
						      << ", z = " << xyzIde[2]
						      << " does not correspond to a TPC.";
	              continue;    }
	
	// NearestWire seems to be missing some correction that is applied in LArG4, and is
	// sometimes off by one or two wires. Use it and then correct it given channel
	/// \todo: Why would an IDE ossociated with a hit return a nearest wire not on the hit's channel? Usually only one off.
	geo::WireID IdeWid = geom->NearestWireID(  xyzIde,
						   cwids[0].Plane,
						   tpc,   cryo  );
	geo::WireID storethis = IdeWid; // default...
	bool foundmatch(false);	
	for( size_t w=0; w<cwids.size(); w++ ){	 
	  if( std::abs((int)(IdeWid.Wire) - (int)(cwids[w].Wire)) <= 1 ){
	    storethis = cwids[w]; // ...apply correction
	    foundmatch = true;
	    break;
	  } 	  
	}
	if( !foundmatch ){
	  mf::LogWarning("DisambigCheat") << "IDE NearestWire return more than 1 off from channel wids: wire " 
					    << IdeWid.Wire; 
	  fBadIDENearestWire++;
	}
	
	bool alreadyStored(false);
	for( size_t wid=0; wid<widsWithIdes.size(); wid++ ) if( storethis == widsWithIdes[wid] ) alreadyStored = true;
	if( !alreadyStored ) widsWithIdes.push_back(storethis);
      } // end loop through ides from HitToSimIdes
      
      fHitToWids[ChanTime] = widsWithIdes;
      
    } // end U/V channel hit loop
    
    
    if(fHitToWids.size() != Ucount + Vcount) throw cet::exception("DisambigCheat");
    
    return;
}


//-------------------------------------------------------------------
void DisambigCheater::beginJob()
{
  return;
}


//-------------------------------------------------------------------
void DisambigCheater::endJob()
{

  if(fFalseChanHits > 0)
    mf::LogWarning("DisambigCheater") << fFalseChanHits << " hits had no associated IDE or WireIDs";

  return;
}


//-------------------------------------------------------------------
void DisambigCheater::reconfigure(fhicl::ParameterSet const & p)
{
  fChanHitLabel =  p.get< std::string >("ChanHitLabel");
  return;
}

#endif // DisambigCheater_h

DEFINE_ART_MODULE(DisambigCheater)

} // end hit namespace
