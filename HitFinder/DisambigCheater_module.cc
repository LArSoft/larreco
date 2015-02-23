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

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindOneP.h"
#include "art/Framework/Principal/Event.h"

#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Simulation/SimChannel.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBaseArt/HitCreator.h"
#include "Utilities/DetectorProperties.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "MCCheater/BackTracker.h"


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
    
    std::string fChanHitLabel;
    std::string fWidHitLabel;
    
    std::map< std::pair<double,double>, std::vector<geo::WireID> > fHitToWids;
    void InitHitToWids( const std::vector< art::Ptr<recob::Hit> >& ChHits );

    void MakeDisambigHit(
      art::Ptr<recob::Hit> const& hit,
      geo::WireID const& wid,
      art::Ptr<recob::Wire> const& wire,
      art::Ptr<raw::RawDigit> const& rawdigits,
      recob::HitCollectionCreator& hcol
      );
    
    unsigned int fFalseChanHits = 0; // hits with no IDEs - could be noise or other technical problems
    unsigned int fBadIDENearestWire = 0; // Just to count the inconsistent IDE wire returns
    
  };
  
  //-------------------------------------------------------------------
  DisambigCheater::DisambigCheater(fhicl::ParameterSet const & p)
  {
    this->reconfigure(p);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
  }
  
  //-------------------------------------------------------------------
  hit::DisambigCheater::~DisambigCheater()
  {
  }
  
  //-------------------------------------------------------------------
  void DisambigCheater::produce(art::Event & evt)
  {
    
    // get hits on channels
    art::Handle< std::vector<recob::Hit> > ChanHits;
    evt.getByLabel(fChanHitLabel, ChanHits);
    std::vector< art::Ptr<recob::Hit> > ChHits;
    art::fill_ptr_vector(ChHits, ChanHits);
  
    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> ChannelHitRawDigits
      (ChanHits, evt, fChanHitLabel);
    const bool doRawDigitAssns = ChannelHitRawDigits.isValid();
    
    art::FindOneP<recob::Wire> ChannelHitWires(ChanHits, evt, fChanHitLabel);
    const bool doWireAssns = ChannelHitWires.isValid();
    
    // this object contains the hit collection
    // and its associations to wires and raw digits
    // (if the original objects have them):
    recob::HitCollectionCreator hits(*this, evt, doWireAssns, doRawDigitAssns);
    

    // find the wireIDs each hit is on
    this->InitHitToWids( ChHits );
    
    
    // make all of the hits
    for(size_t h=0; h<ChHits.size(); h++){
      
      // get the objects associated with this hit
      art::Ptr<recob::Wire> wire;
      if (doWireAssns) wire = ChannelHitWires.at(h);
      art::Ptr<raw::RawDigit> rawdigits;
      if (doRawDigitAssns) rawdigits = ChannelHitRawDigits.at(h);
      
      // the trivial Z hits
      if(ChHits[h]->View()==geo::kZ){
        this->MakeDisambigHit(ChHits[h], ChHits[h]->WireID(), wire, rawdigits, hits);
        continue;
      }

      
      // make U/V hits if any wire IDs are associated
      // count hits without a wireID
      /// \todo: Decide how to handle the hits with multiple-wid activity. For now, randomly choose.
      std::pair<double,double> ChanTime(ChHits[h]->Channel()*1., ChHits[h]->PeakTime()*1.);
      if( fHitToWids[ChanTime].size() == 1 )
        this->MakeDisambigHit(ChHits[h], fHitToWids[ChanTime][0], wire, rawdigits, hits);
      else if( fHitToWids[ChanTime].size() == 0 )
        fFalseChanHits++;
      else if( fHitToWids[ChanTime].size() > 1 )
        this->MakeDisambigHit(ChHits[h], fHitToWids[ChanTime][0], wire, rawdigits, hits); // same thing for now
      
    } // for
    
    
    // put the hit collection and associations into the event
    hits.put_into(evt);
    
    fHitToWids.clear();
    return;
    
  }
  
  
  //-------------------------------------------------
  void DisambigCheater::MakeDisambigHit(
    art::Ptr<recob::Hit> const& original_hit, geo::WireID const& wid,
    art::Ptr<recob::Wire> const& wire, art::Ptr<raw::RawDigit> const& rawdigits,
    recob::HitCollectionCreator& hcol
    )
  {
    
    if( !wid.isValid ){ 
      mf::LogWarning("InvalidWireID") << "wid is invalid, hit not being made\n";
      return;
    }
    
    // create a hit copy of the original one, but with a different wire ID
    recob::HitCreator hit(*original_hit, wid);
    hcol.emplace_back(hit.move(), wire, rawdigits);
    
  }
  
  
  //-------------------------------------------------------------------
  void DisambigCheater::InitHitToWids( const std::vector< art::Ptr<recob::Hit> >& ChHits )
  {
    
    unsigned int Ucount(0), Vcount(0);
    for( size_t h = 0; h < ChHits.size(); h++ ){
      recob::Hit const& chit = *(ChHits[h]);
      if( ChHits[h]->View() == geo::kZ ) continue;
      if( ChHits[h]->View() == geo::kU ) Ucount++;
      else if( ChHits[h]->View() == geo::kV ) Vcount++;
      std::vector<geo::WireID> cwids = geom->ChannelToWire(chit.Channel());
      std::pair<double,double> ChanTime( (double) chit.Channel(), (double) chit.PeakTime()); // hit key value 
      
      // get hit IDEs
      std::vector< sim::IDE > ides;
      bt->HitToSimIDEs( chit, ides );
      
      // catch the hits that have no IDEs
      bool hasIDEs = !ides.empty();
      if (hasIDEs) {
        try        { bt->SimIDEsToXYZ(ides); }
        catch(...) { hasIDEs = false; }
      } // if
      if (!hasIDEs) {
        mf::LogVerbatim("DisambigCheat") << "Hit on channel " << chit.Channel()
                                         << " between " << chit.PeakTimeMinusRMS()
                                         << " and " << chit.PeakTimePlusRMS() << " matches no IDEs.";
        fHitToWids[ChanTime] = std::vector< geo::WireID >();
        continue;
      } // if no IDEs
      
      // see what wire ID(s) the hit-IDEs are on
      //  if none: hit maps to empty vector
      //  if one: we have a unique hit, vector of size 1
      //  if more than one: make a vector of all wids
      std::vector<geo::WireID> widsWithIdes; 
      for( size_t i=0; i<ides.size(); i++ ){
	const double xyzIde[] = { ides[i].x, ides[i].y, ides[i].z };

	// Occasionally, ide position is not in TPC
	/// \todo: Why would an IDE xyz position be outside of a TPC?
	geo::TPCID tpcID = geom->FindTPCAtPosition(xyzIde);
	if (!tpcID.isValid) {
	  mf::LogWarning("DisambigCheat") << "IDE at x = " << xyzIde[0] 
	                                  << ", y = " << xyzIde[1]
	                                  << ", z = " << xyzIde[2]
	                                  << " does not correspond to a TPC.";
	  continue;
	}
	unsigned int tpc = tpcID.TPC, cryo = tpcID.Cryostat;
	
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
    
    if(fHitToWids.size() != Ucount + Vcount){
      //throw cet::exception("DisambigCheat");
      mf::LogWarning("DisambigCheat")<<"Nhits mismatch: "<<fHitToWids.size()<<" "<<Ucount+Vcount;
    }
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
