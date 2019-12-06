////////////////////////////////////////////////////////////////////////
// Class:       DisambigCheater
// Module Type: producer
// File:        DisambigCheater_module.cc
//
// tylerdalion@gmail.com
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"


namespace hit{

  class DisambigCheater : public art::EDProducer {
    public:
      explicit DisambigCheater(fhicl::ParameterSet const & p);

    private:
      void produce(art::Event & e) override;
      void endJob() override;

      art::ServiceHandle<geo::Geometry const> geom;
      art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

      std::string fChanHitLabel;
      std::string fWidHitLabel;

      std::map< std::pair<double,double>, std::vector<geo::WireID> > fHitToWids;
      void InitHitToWids( const std::vector< art::Ptr<recob::Hit> >& ChHits );

      void MakeDisambigHit(
          art::Ptr<recob::Hit> const& hit,
          geo::WireID const& wid,
          art::Ptr<recob::Wire> const& wire,
          //art::Ptr<raw::RawDigit> const& rawdigits,
          recob::HitCollectionCreator& hcol
          );

      unsigned int fFalseChanHits = 0; // hits with no IDEs - could be noise or other technical problems
      unsigned int fBadIDENearestWire = 0; // Just to count the inconsistent IDE wire returns

      std::vector<unsigned int> fMaxWireShift; // shift to account for space charge and diffusion

  };

  //-------------------------------------------------------------------
  DisambigCheater::DisambigCheater(fhicl::ParameterSet const & p)
    : EDProducer{p}
  {
    fChanHitLabel =  p.get< std::string >("ChanHitLabel");

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(producesCollector(), "", true, false);


    // Space charge can shift true IDE postiion to far-off channels.
    // Calculate maximum number of wires to shift hit in order to be on correct channel.
    // Shift no more than half of the number of channels, as beyond there will
    //   have been a closer wire segment.
    if( geom->Ncryostats()!=1 || geom->NTPC()<1 ){
      fMaxWireShift.resize(3);
      fMaxWireShift[0]=1;
      fMaxWireShift[1]=1;
      fMaxWireShift[2]=1;
    } else {
      // assume TPC 0 is typical of all in terms of number of channels
      unsigned int np = geom->Cryostat(0).TPC(0).Nplanes();
      fMaxWireShift.resize(np);
      for(unsigned int p = 0; p<np; ++p){
        double xyz[3] = {0.};
        double xyz_next[3] = {0.};
        unsigned int nw = geom->Cryostat(0).TPC(0).Plane(p).Nwires();
        for(unsigned int w = 0; w<nw; ++w){

          // for vertical planes
          if(geom->Cryostat(0).TPC(0).Plane(p).View()==geo::kZ)   {
            fMaxWireShift[2] = geom->Cryostat(0).TPC(0).Plane(p).Nwires();
            break;
          }

          geom->Cryostat(0).TPC(0).Plane(p).Wire(w).GetCenter(xyz);
          geom->Cryostat(0).TPC(0).Plane(p).Wire(w+1).GetCenter(xyz_next);

          if(xyz[2]==xyz_next[2]){
            fMaxWireShift[p] = w;
            break;
          }
        }// end wire loop
      }// end plane loop

      for(unsigned int i=0; i<np; i++)
        fMaxWireShift[i] = std::floor(fMaxWireShift[i]/2);
    }

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
//    art::FindOneP<raw::RawDigit> ChannelHitRawDigits
//      (ChanHits, evt, fChanHitLabel);
//    const bool doRawDigitAssns = ChannelHitRawDigits.isValid();
    const bool doRawDigitAssns = false;

    art::FindOneP<recob::Wire> ChannelHitWires(ChanHits, evt, fChanHitLabel);
    const bool doWireAssns = ChannelHitWires.isValid();

    // this object contains the hit collection
    // and its associations to wires and raw digits
    // (if the original objects have them):
    recob::HitCollectionCreator hits(evt, doWireAssns, doRawDigitAssns);


    // find the wireIDs each hit is on
    this->InitHitToWids( ChHits );


    // make all of the hits
    for(size_t h=0; h<ChHits.size(); h++){

      // get the objects associated with this hit
      art::Ptr<recob::Wire> wire;
      if (doWireAssns) wire = ChannelHitWires.at(h);
//      art::Ptr<raw::RawDigit> rawdigits;
//      if (doRawDigitAssns) rawdigits = ChannelHitRawDigits.at(h);

      // the trivial Z hits
      if(ChHits[h]->View()==geo::kZ){
        this->MakeDisambigHit(ChHits[h], ChHits[h]->WireID(), wire, hits);
        continue;
      }


      // make U/V hits if any wire IDs are associated
      // count hits without a wireID
      /// \todo: Decide how to handle the hits with multiple-wid activity. For now, randomly choose.
      std::pair<double,double> ChanTime(ChHits[h]->Channel()*1., ChHits[h]->PeakTime()*1.);
      if( fHitToWids[ChanTime].size() == 1 )
        this->MakeDisambigHit(ChHits[h], fHitToWids[ChanTime][0], wire, hits);
      else if( fHitToWids[ChanTime].size() == 0 )
        fFalseChanHits++;
      else if( fHitToWids[ChanTime].size() > 1 )
        this->MakeDisambigHit(ChHits[h], fHitToWids[ChanTime][0], wire, hits); // same thing for now

    } // for


    // put the hit collection and associations into the event
    hits.put_into(evt);

    fHitToWids.clear();
  }


  //-------------------------------------------------
  void DisambigCheater::MakeDisambigHit(
      art::Ptr<recob::Hit> const& original_hit, geo::WireID const& wid,
      art::Ptr<recob::Wire> const& wire,// art::Ptr<raw::RawDigit> const& rawdigits,
      recob::HitCollectionCreator& hcol
      )
  {

    if( !wid.isValid ){
      mf::LogWarning("InvalidWireID") << "wid is invalid, hit not being made\n";
      return;
    }

    // create a hit copy of the original one, but with a different wire ID
    recob::HitCreator hit(*original_hit, wid);
    hcol.emplace_back(hit.move(), wire);

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
      std::vector<const sim::IDE* > ides;
      try{
        ides=bt_serv->HitToSimIDEs_Ps( chit);
      }
      catch(...){};

      // catch the hits that have no IDEs
      bool hasIDEs = !ides.empty();
      if (hasIDEs) {
        try        { bt_serv->SimIDEsToXYZ(ides); } //What is this supposed to be doing? It get a vector, but never assigns it anywhere. There has to be a better way to do this check
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
        const double xyzIde[] = { ides[i]->x, ides[i]->y, ides[i]->z };

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
        geo::WireID IdeWid;
        try {
          IdeWid = geom->NearestWireID(xyzIde, cwids[0].Plane, tpc, cryo);
        }
        catch (geo::InvalidWireError const& e) { // adopt suggestion if possible
          if (!e.hasSuggestedWire()) throw;
          IdeWid = e.suggestedWireID();
          mf::LogError("DisambigCheat") << "Detected a point out of its wire plane:\n"
            << e.what() << "\nUsing suggested wire " << IdeWid << "\n";
        }
        geo::WireID storethis = IdeWid; // default...
        bool foundmatch(false);
        for( size_t w=0; w<cwids.size(); w++ ){
          if (cwids[w].TPC!=tpc || cwids[w].Cryostat!=cryo) continue;
          if( (unsigned int)std::abs((int)(IdeWid.Wire) - (int)(cwids[w].Wire)) <= fMaxWireShift[cwids[0].Plane] ){
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
  }


  //-------------------------------------------------------------------
  void DisambigCheater::endJob()
  {

    if(fFalseChanHits > 0)
      mf::LogWarning("DisambigCheater") << fFalseChanHits << " hits had no associated IDE or WireIDs";
  }



  DEFINE_ART_MODULE(DisambigCheater)

} // end hit namespace
