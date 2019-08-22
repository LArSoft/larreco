/*!
 * Title:   TTSpacePointFinder class
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Hit
 * Outputs: recob::SpacePoint
 *
 * Description:
 * This module, TimeTickSpacePointFinder (or TTSpacePointFinder for short) was
 * originally designed to produce a spacepoint object based on hits from
 * TTHitFinder, with the  intention to allow for a significant number of
 * ghost spacepoints, where some downstream application would deal with the
 * results.
 *
 * In reality, it is just a generic SpacePointFinder implementing the
 * SpacePointAlg_TimeSort algorithm.
 *
 * This code is totally microboone specific, btw.
 */

#include <string>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"

// LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/SpacePointAlg_TimeSort.h"

namespace sppt{

  class TTSpacePointFinder : public art::EDProducer {

  public:

    explicit TTSpacePointFinder(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt);
    void beginRun(art::Run& run);

  private:

    std::string    fHitModuleLabel;     /// Input hit module name
    std::string    fUHitsInstanceLabel; /// Input U hits instance name
    std::string    fVHitsInstanceLabel; /// Input V hits instance name
    std::string    fYHitsInstanceLabel; /// Input Y hits instance name

    SpacePointAlg_TimeSort fSpptAlg;

  protected:

  }; // class TTSpacePointFinder

  //-------------------------------------------------
  TTSpacePointFinder::TTSpacePointFinder(fhicl::ParameterSet const& pset):
    EDProducer{pset},
    fSpptAlg( pset.get< fhicl::ParameterSet >("SpacePointAlgParams") )
  {
    fHitModuleLabel     = pset.get< std::string >("HitModuleLabel","tthit");
    fUHitsInstanceLabel = pset.get< std::string >("UHitsInstaceLabel","uhits");
    fVHitsInstanceLabel = pset.get< std::string >("VHitsInstaceLabel","vhits");
    fYHitsInstanceLabel = pset.get< std::string >("YHitsInstaceLabel","yhits");

    produces< std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::SpacePoint, recob::Hit>       >();
  }

  //-------------------------------------------------
  void TTSpacePointFinder::beginRun(art::Run &run){
    fSpptAlg.setTimeOffsets();
    fSpptAlg.fillCoordinatesArrays();
  }

  //-------------------------------------------------
  void TTSpacePointFinder::produce(art::Event& evt)
  {

    //initialize our spacepoint collection
    std::unique_ptr<std::vector<recob::SpacePoint> > spptCollection(new std::vector<recob::SpacePoint>);
    std::unique_ptr<std::vector<std::vector<art::Ptr<recob::Hit> > > >
      spptAssociatedHits(new std::vector<std::vector<art::Ptr<recob::Hit> > >);

    // Read in the hits. Note, we will reorder hit vector, so we do in fact need a copy.
    art::Handle< std::vector<recob::Hit> > hitHandle_U;
    evt.getByLabel(fHitModuleLabel,fUHitsInstanceLabel,hitHandle_U);
    std::vector< art::Ptr<recob::Hit> > hitVec_U;
    art::fill_ptr_vector(hitVec_U,hitHandle_U);

    art::Handle< std::vector<recob::Hit> > hitHandle_V;
    evt.getByLabel(fHitModuleLabel,fVHitsInstanceLabel,hitHandle_V);
    std::vector< art::Ptr<recob::Hit> > hitVec_V;
    art::fill_ptr_vector(hitVec_V,hitHandle_V);

    art::Handle< std::vector<recob::Hit> > hitHandle_Y;
    evt.getByLabel(fHitModuleLabel,fYHitsInstanceLabel,hitHandle_Y);
    std::vector< art::Ptr<recob::Hit> > hitVec_Y;
    art::fill_ptr_vector(hitVec_Y,hitHandle_Y);

    MF_LOG_DEBUG("TTSpacePointFinder")
      << "Got handles to hits:\n"
      << hitVec_U.size() << " u hits, "
      << hitVec_V.size() << " v hits, "
      << hitVec_Y.size() << " y hits.";

    //now, call the space point alg
    fSpptAlg.createSpacePoints(hitVec_U,hitVec_V,hitVec_Y,
                               spptCollection,spptAssociatedHits);

    mf::LogInfo("TTSpacePointFinder") << "Finished spacepoint alg. Created " << spptCollection->size() << " spacepoints.";


    //check that spptCollection and spptAssociatedHits have same size
    if(spptCollection->size() != spptAssociatedHits->size()){
      mf::LogError("TTSpacePointFinder") << "ERROR in space point alg.\n"
                                         << "Spacepoint and associated hit collections have different sizes!\n"
                                         << spptCollection->size() << " != " << spptAssociatedHits->size()
                                         << "\nWill return with no space points put on event.";
      return;
    }

    //Now, need to fill in the sppt<-->hit associations
    std::unique_ptr<art::Assns<recob::SpacePoint,recob::Hit> > spptAssns(new art::Assns<recob::SpacePoint,recob::Hit>);
    for(unsigned int isppt=0; isppt<spptCollection->size(); isppt++){
      util::CreateAssn(*this,evt,*spptCollection,spptAssociatedHits->at(isppt),*spptAssns,isppt);
    }

    //finally, put things on the event
    evt.put(std::move(spptCollection));
    evt.put(std::move(spptAssns));

  } // End of produce()

  DEFINE_ART_MODULE(TTSpacePointFinder)

} // end of sppt namespace
