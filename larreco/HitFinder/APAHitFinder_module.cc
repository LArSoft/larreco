////////////////////////////////////////////////////////////////////////
//
// APAHitFinder class
//
// talion@gmail.com
//
//  This algorithm is designed to find hits on APA channels after
//  deconvolution, and then disambiguate those hits, attempting to
//  localize the hit to one segment, on one side of the APA.
//
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr()
#include <utility> // std::move()

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"


// LArSoft Includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larreco/RecoAlg/DisambigAlg.h"

namespace apa{
  class APAHitFinder : public art::EDProducer {

  public:

    explicit APAHitFinder(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    apa::DisambigAlg    fDisambigAlg;
    art::ServiceHandle<geo::Geometry const> fGeom;

    std::string fChanHitLabel;

  }; // class APAHitFinder


//-------------------------------------------------
//-------------------------------------------------
APAHitFinder::APAHitFinder(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  , fDisambigAlg(pset.get< fhicl::ParameterSet >("DisambigAlg"))
{
  fChanHitLabel =  pset.get< std::string >("ChanHitLabel");

  // let HitCollectionCreator declare that we are going to produce
  // hits and associations with wires and raw digits
  // (with no particular product label)
  recob::HitCollectionCreator::declare_products(producesCollector());
}


//-------------------------------------------------
void APAHitFinder::produce(art::Event& evt)
{
  // this object contains the hit collection
  // and its associations to wires and raw digits:
  recob::HitCollectionCreator hcol(evt);

  art::Handle< std::vector<recob::Hit> > ChannelHits;
  evt.getByLabel(fChanHitLabel, ChannelHits);

  // also get the associated wires and raw digits;
  // we assume they have been created by the same module as the hits
  art::FindOneP<raw::RawDigit> ChannelHitRawDigits
    (ChannelHits, evt, fChanHitLabel);
  art::FindOneP<recob::Wire> ChannelHitWires
    (ChannelHits, evt, fChanHitLabel);

  // Make unambiguous collection hits
  std::vector< art::Ptr<recob::Hit> >  ChHits;
  art::fill_ptr_vector(ChHits, ChannelHits);
  for( size_t h = 0; h < ChHits.size(); h++ ){
    if( ChHits[h]->View() != geo::kZ ) continue;

    art::Ptr<recob::Wire> wire = ChannelHitWires.at(h);
    art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(h);

    // just copy it
    hcol.emplace_back(*ChHits[h], wire, rawdigits);
  }


  // Run alg on all APAs
  fDisambigAlg.RunDisambig(ChannelHits);


  for( size_t t=0; t < fDisambigAlg.fDisambigHits.size(); t++ ){
    art::Ptr<recob::Hit>  hit = fDisambigAlg.fDisambigHits[t].first;
    geo::WireID           wid = fDisambigAlg.fDisambigHits[t].second;

    // create a new hit copy of the original one, but with new wire ID
    recob::HitCreator disambiguous_hit(*hit, wid);

    // get the objects associated with the original hit;
    // since hit comes from ChannelHits, its key is the index in that collection
    // and also the index for the query of associated objects
    art::Ptr<recob::Hit>::key_type hit_index = hit.key();
    art::Ptr<recob::Wire> wire = ChannelHitWires.at(hit_index);
    art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(hit_index);

    hcol.emplace_back(disambiguous_hit.move(), wire, rawdigits);
  } // for

  // put the hit collection and associations into the event
  hcol.put_into(evt);

}


DEFINE_ART_MODULE(APAHitFinder)

} // end of apa namespace
