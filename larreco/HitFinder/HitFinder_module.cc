/**
 * @file   HitFinder_module.cc
 * @brief  Hit finder (originating for cluster crawler algorithm)
 * @author Bruce Baller (bballer@fnal.gov)
 * 
 * Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod 
 * from cetpkgsupport v1_02_00.
 */


// C/C++ standard libraries
#include <string>
#include <utility> // std::unique_ptr<>

// Framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

//LArSoft includes
#include "larreco/RecoAlg/CCHitFinderAlg.h"

// ... more includes in the implementation section


namespace hit {
  
  class HitFinder: public art::EDProducer {
  
    public:
      explicit HitFinder(fhicl::ParameterSet const & pset);
      virtual ~HitFinder() = default;
  
      virtual void reconfigure(fhicl::ParameterSet const & pset) override;
      virtual void produce(art::Event & evt) override;
      
      virtual void endJob() override;
  
    private:
      art::InputTag fCalDataModuleLabel; ///< label of module producing input wires
      std::unique_ptr<CCHitFinderAlg> fCCHFAlg; // define CCHitFinderAlg object
      
  }; // hit::HitFinder()
  
} // namespace hit

//******************************************************************************
//***  implementation
//***

// C/C++ standard libraries
#include <vector>
#include <memory> // std::move()

// Framework libraries
#include "art/Framework/Principal/Handle.h"

//LArSoft includes
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h" // recob::HitCollectionAssociator


namespace hit {
  
  
  //----------------------------------------------------------------------------
  HitFinder::HitFinder(fhicl::ParameterSet const& pset)
  {
    reconfigure(pset);
    
    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label);
    // TODO this should be marked as transient when art will implement issue #8018
    recob::HitCollectionAssociator::declare_products(*this);
    
  } // HitFinder::HitFinder()
  
  
  //----------------------------------------------------------------------------
  void HitFinder::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCalDataModuleLabel = pset.get<art::InputTag>("CalDataModuleLabel");
    
    // this trick avoids double configuration on construction
    if (fCCHFAlg)
      fCCHFAlg->reconfigure(pset.get<fhicl::ParameterSet>("CCHitFinderAlg"));
    else {
      fCCHFAlg.reset
        (new CCHitFinderAlg(pset.get<fhicl::ParameterSet>("CCHitFinderAlg")));
    }
  } // HitFinder::reconfigure()
  
  
  //----------------------------------------------------------------------------
  void HitFinder::produce(art::Event & evt)
  {
    // fetch the wires needed by HitFinder

    // make this accessible to ClusterCrawler_module
    art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
      = evt.getValidHandle<std::vector<recob::Wire>>(fCalDataModuleLabel);

    // find hits in all planes
    fCCHFAlg->RunCCHitFinder(*wireVecHandle);
    
    // extract the result of the algorithm (it's moved)
    std::unique_ptr<std::vector<recob::Hit>> Hits
      (new std::vector<recob::Hit>(std::move(fCCHFAlg->YieldHits())));
    
    mf::LogInfo("HitFinder") << Hits->size() << " hits produced.";
    
    // shcol contains the hit collection
    // and its associations to wires and raw digits;
    // we get the association to raw digits through wire associations
    recob::HitCollectionAssociator shcol(*this, evt, fCalDataModuleLabel, true);

    shcol.use_hits(std::move(Hits));
    
    // move the hit collection and the associations into the event:
    shcol.put_into(evt);
    
  } // produce()
  
  
  //----------------------------------------------------------------------------
  void HitFinder::endJob() {
    // print the statistics about fits
    mf::LogInfo log("HitFinder"); // messages are printed on "log" destruction
    fCCHFAlg->PrintStats(log);
  } // HitFinder::endJob()
  
  
  DEFINE_ART_MODULE(HitFinder)
  
} // namespace hit
