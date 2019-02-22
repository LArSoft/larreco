/**
 * @file   TCTrack_module.cc
 * @brief  Create seeds, spacepoints and tracks from 3D matched clusters produced by TrajCluster_module
 * @author Bruce Baller (baller@fnal.gov)
 * 
*
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
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "lardata/Utilities/AssociationUtil.h"

// ... more includes in the implementation section

namespace trkf {

  class TCTrack: public art::EDProducer {
    
  public:
    explicit TCTrack(fhicl::ParameterSet const & pset);
    
    void reconfigure(fhicl::ParameterSet const & pset) ;
    void produce(art::Event & evt) override;
    
  private:

    SpacePointAlg fSptalg;

    std::string            fPFPModuleLabel;
//    std::string            fHitModuleLabel;

    
  }; // class TCTrack
  
  //----------------------------------------------------------------------------
  void TCTrack::reconfigure(fhicl::ParameterSet const & pset)
  {
//    fHitModuleLabel         = pset.get< std::string >("HitModuleLabel");
    fPFPModuleLabel     = pset.get< std::string >("PFPModuleLabel");
    
    fSptalg.reconfigure(pset.get<fhicl::ParameterSet>("SpacePointAlg"));

  }
    
  //----------------------------------------------------------------------------
  TCTrack::TCTrack(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fSptalg(pset.get<fhicl::ParameterSet>("SpacePointAlg"))
  {
    
    reconfigure(pset);
    
    produces<std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::SpacePoint, recob::Hit>       >();
    
  } // TCTrack::TCTrack()
  
  //----------------------------------------------------------------------------
  void TCTrack::produce(art::Event & evt)
  {
    
    // all data products are assumed to be produced by the same module that produced the PFParticles -> TrajCluster_module
    art::InputTag DataInputTag(fPFPModuleLabel);
    art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(DataInputTag);

    art::ValidHandle<std::vector<recob::Cluster>> clsHandle = evt.getValidHandle<std::vector<recob::Cluster>>(DataInputTag);

    art::FindManyP<recob::Cluster> pfp_cls(pfpHandle, evt, fPFPModuleLabel);
    art::FindManyP<recob::Hit> cls_hit(clsHandle, evt, fPFPModuleLabel);
    
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> > sphitassn(new art::Assns<recob::SpacePoint, recob::Hit>);
    std::unique_ptr<std::vector<recob::SpacePoint> > spts(new std::vector<recob::SpacePoint>);
    
    art::PtrVector<recob::Hit> hits;
    for(unsigned short ipfp = 0; ipfp < pfpHandle->size(); ++ipfp) {
      // Get the clusters associated with this PFParticle - there should be one in each plane
      std::vector<art::Ptr<recob::Cluster> > clsList;
      pfp_cls.get(ipfp, clsList);
      hits.clear();
      std::cout<<"PFP "<<ipfp<<"\n";
      for(unsigned short icl = 0; icl < clsList.size(); ++icl) {
        std::vector<art::Ptr<recob::Hit> > hitList;
        unsigned int clsIndex = clsList[icl]->ID() - 1;
        cls_hit.get(clsIndex, hitList);
        std::cout<<" cls index "<<clsIndex<<" hits size "<<hitList.size()<<"  "<<(int)clsList[icl]->StartWire()<<":"<<(int)clsList[icl]->StartTick()<<"  EndWire "<<(int)clsList[icl]->EndWire()<<":"<<(int)clsList[icl]->EndTick()<<"\n";
        hits.reserve(hits.size() + hitList.size());
        for(std::vector< art::Ptr<recob::Hit> >::const_iterator i = hitList.begin(); i != hitList.end(); ++i) hits.push_back(*i);
      } // icl
      // make new space points using these hits
      std::vector<recob::SpacePoint> new_spts;
      fSptalg.makeSpacePoints(hits, new_spts);
      if(new_spts.empty()) continue;
/*
      if(ipfp == 0) {
        std::cout<<"new_spts size "<<new_spts.size()<<"\n";
        for(auto& spt : new_spts) {
          std::cout<<" "<<spt<<"\n";
        } // spt
        std::cout<<"\n";
      } // ipfp == 0
*/
      int nspt = spts->size();
      spts->insert(spts->end(), new_spts.begin(), new_spts.end());
      // associate the hits with the spacepoints
      art::PtrVector<recob::SpacePoint> sptvec;
      for(unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
        const recob::SpacePoint& spt = (*spts)[ispt];
        const art::PtrVector<recob::Hit>& hits = fSptalg.getAssociatedHits(spt);
        util::CreateAssn(*this, evt, *spts, hits, *sphitassn, ispt);
      } // ispt
    } // ipfp
    
    evt.put(std::move(spts));
    evt.put(std::move(sphitassn));
      
  } // TCTrack::produce()
    
  DEFINE_ART_MODULE(TCTrack)

} // namespace trkf
