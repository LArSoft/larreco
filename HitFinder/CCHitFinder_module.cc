/**
 * @file   CCHitFinder_module.cc
 * @brief  Hit finder for cluster crawler algorithm
 * @author Bruce Baller (bballer@fnal.gov)
 * 
 * Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod 
 * from cetpkgsupport v1_02_00.
 */


// C/C++ standard libraries
#include <string>

// Framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/InputTag.h"

//LArSoft includes
#include "RecoAlg/CCHitFinderAlg.h"

// ... more includes in the implementation section


namespace hit {
  class CCHitFinder;
}

class hit::CCHitFinder: public art::EDProducer {

  public:
    explicit CCHitFinder(fhicl::ParameterSet const & pset);
    virtual ~CCHitFinder() = default;

    void reconfigure(fhicl::ParameterSet const & pset) override;
    void produce(art::Event & evt) override;

  private:
    std::string fCalDataModuleLabel; ///< label of module producing input wires
    CCHitFinderAlg fCCHFAlg; // define CCHitFinderAlg object
}; // hit::CCHitFinder()


//******************************************************************************
//***  implementation
//***

// C/C++ standard libraries
#include <vector>
#include <memory> // std::move()
#include <utility> // std::unique_ptr<>

// Framework libraries
#include "art/Framework/Principal/Handle.h"

//LArSoft includes
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBaseArt/HitCreator.h" // recob::HitCollectionAssociator


namespace hit {
  
  
  //----------------------------------------------------------------------------
  CCHitFinder::CCHitFinder(fhicl::ParameterSet const& pset) :
    // TODO replace with art::InputTag when LArSoft uses art >=1_13_00
    fCalDataModuleLabel(pset.get<std::string>("CalDataModuleLabel")),
    fCCHFAlg           (pset.get<fhicl::ParameterSet>("CCHitFinderAlg"))
  {
    this->reconfigure(pset);
    
    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label);
    // TODO this should be marked as transient when art will implement issue #8018
    recob::HitCollectionAssociator::declare_products(*this);
    
  } // CCHitFinder::CCHitFinder()
  
  
  //----------------------------------------------------------------------------
  void CCHitFinder::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCCHFAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitFinderAlg"));
  }
  
  
  //----------------------------------------------------------------------------
  void CCHitFinder::produce(art::Event & evt)
  {
    // fetch the wires needed by CCHitFinder

    // make this accessible to ClusterCrawler_module
    art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
     = evt.getValidHandle<std::vector<recob::Wire>>(fCalDataModuleLabel);

    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(*wireVecHandle);
    
    // extract the result of the algorithm (it's moved)
    std::unique_ptr<std::vector<recob::Hit>> Hits
      (new std::vector<recob::Hit>(std::move(fCCHFAlg.YieldHits())));
    
    mf::LogInfo("CCHitFinder") << Hits->size() << " hits produced.";
    
    // shcol contains the hit collection
    // and its associations to wires and raw digits;
    // we get the association to raw digits through wire associations
    recob::HitCollectionAssociator shcol(*this, evt, fCalDataModuleLabel, true);

    shcol.use_hits(std::move(Hits));
    
    // move the hit collection and the associations into the event:
    shcol.put_into(evt);

  } // produce()
  

  DEFINE_ART_MODULE(CCHitFinder)
  
} // namespace hit
