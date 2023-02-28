/**
 * @file   TCTrack_module.cc
 * @brief  Create seeds, spacepoints and tracks from 3D matched clusters produced by TrajCluster_module
 * @author Bruce Baller (baller@fnal.gov)
 */

// C/C++ standard libraries
#include <string>

// Framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/SpacePointAlg.h"

// ... more includes in the implementation section

namespace trkf {

  class TCTrack : public art::EDProducer {

  public:
    explicit TCTrack(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;

    SpacePointAlg fSptalg;

    art::InputTag fPFPtag;
  }; // class TCTrack

  //----------------------------------------------------------------------------
  TCTrack::TCTrack(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fSptalg{pset.get<fhicl::ParameterSet>("SpacePointAlg")}
    , fPFPtag{pset.get<std::string>("PFPModuleLabel")}
  {
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
  }

  //----------------------------------------------------------------------------
  void TCTrack::produce(art::Event& evt)
  {
    // all data products are assumed to be produced by the same module that produced the PFParticles -> TrajCluster_module
    auto pfpHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFPtag);
    auto clsHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fPFPtag);

    art::FindManyP<recob::Cluster> pfp_cls(pfpHandle, evt, fPFPtag);
    art::FindManyP<recob::Hit> cls_hit(clsHandle, evt, fPFPtag);

    auto sphitassn = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();
    auto spts = std::make_unique<std::vector<recob::SpacePoint>>();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    art::PtrVector<recob::Hit> hits;
    for (unsigned short ipfp = 0; ipfp < pfpHandle->size(); ++ipfp) {
      // Get the clusters associated with this PFParticle - there should be one in each plane
      std::vector<art::Ptr<recob::Cluster>> clsList;
      pfp_cls.get(ipfp, clsList);
      hits.clear();
      std::cout << "PFP " << ipfp << "\n";
      for (unsigned short icl = 0; icl < clsList.size(); ++icl) {
        std::vector<art::Ptr<recob::Hit>> hitList;
        unsigned int clsIndex = clsList[icl]->ID() - 1;
        cls_hit.get(clsIndex, hitList);
        std::cout << " cls index " << clsIndex << " hits size " << hitList.size() << "  "
                  << (int)clsList[icl]->StartWire() << ":" << (int)clsList[icl]->StartTick()
                  << "  EndWire " << (int)clsList[icl]->EndWire() << ":"
                  << (int)clsList[icl]->EndTick() << "\n";
        hits.reserve(hits.size() + hitList.size());
        hits.insert(hits.end(), hitList.begin(), hitList.end());
      } // icl
      // make new space points using these hits
      std::vector<recob::SpacePoint> new_spts;
      fSptalg.makeSpacePoints(clockData, detProp, hits, new_spts);
      if (new_spts.empty()) continue;

      int nspt = spts->size();
      spts->insert(spts->end(), new_spts.begin(), new_spts.end());
      // associate the hits with the spacepoints
      art::PtrVector<recob::SpacePoint> sptvec;
      for (unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
        const recob::SpacePoint& spt = (*spts)[ispt];
        const art::PtrVector<recob::Hit>& hits = fSptalg.getAssociatedHits(spt);
        util::CreateAssn(evt, *spts, hits, *sphitassn, ispt);
      } // ispt
    }   // ipfp

    evt.put(std::move(spts));
    evt.put(std::move(sphitassn));

  } // TCTrack::produce()

  DEFINE_ART_MODULE(TCTrack)

} // namespace trkf
