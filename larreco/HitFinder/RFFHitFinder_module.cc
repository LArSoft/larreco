////////////////////////////////////////////////////////////////////////
// Class:       RFFHitFinder
// Module Type: producer
// File:        RFFHitFinder_module.cc
//
// Generated at Fri Jan 30 17:27:31 2015 by Wesley Ketchum using artmod
// from cetpkgsupport v1_08_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include <memory>

#include "RFFHitFinderAlg.h"

namespace hit{
  class RFFHitFinder;
}

namespace hit{

  class RFFHitFinder : public art::EDProducer {
  public:
    explicit RFFHitFinder(fhicl::ParameterSet const & p);

    // Plugins should not be copied or assigned.
    RFFHitFinder(RFFHitFinder const &) = delete;
    RFFHitFinder(RFFHitFinder &&) = delete;
    RFFHitFinder & operator = (RFFHitFinder const &) = delete;
    RFFHitFinder & operator = (RFFHitFinder &&) = delete;

  private:
    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    art::InputTag   fWireModuleLabel;
    RFFHitFinderAlg fAlg;
  };


  RFFHitFinder::RFFHitFinder(fhicl::ParameterSet const & p)
    :
    EDProducer{p},
    fWireModuleLabel(p.get<std::string>("WireModuleLabel")),
    fAlg(p.get<fhicl::ParameterSet>("RFFHitFinderAlgParams"))
  {
    //calls the produces stuff for me!
    recob::HitCollectionCreator::declare_products(producesCollector());
  }

  void RFFHitFinder::produce(art::Event & e)
  {
    art::ServiceHandle<geo::Geometry const> geoHandle;

    art::Handle< std::vector<recob::Wire> > wireHandle;
    e.getByLabel(fWireModuleLabel,wireHandle);

    std::unique_ptr< std::vector<recob::Hit> > hitCollection(new std::vector<recob::Hit>);
    fAlg.Run(*wireHandle,*hitCollection,*geoHandle);

    recob::HitCollectionAssociator hcol(e,fWireModuleLabel,true);
    hcol.use_hits(std::move(hitCollection));
    hcol.put_into(e);

    //e.put(std::move(hitCollection));
  }

  void RFFHitFinder::beginJob()
  {
    art::ServiceHandle<geo::Geometry const> geoHandle;
    geo::Geometry const& geo(*geoHandle);
    fAlg.SetFitterParamsVectors(geo);
  }

}

DEFINE_ART_MODULE(hit::RFFHitFinder)
