////////////////////////////////////////////////////////////////////////
// Class:       TrackShowerHits
// Module Type: producer
// File:        TrackShowerHits_module.cc
// Authors: dorota.stefan@cern.ch robert.sulej@cern.ch
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larreco/ClusterFinder/TrackShowerSplitter/Segmentation2D/Segmentation2D.h"

#include <memory>

namespace tss {

  typedef std::map<unsigned int, std::vector<tss::Hit2D>> view_hitmap;
  typedef std::map<unsigned int, view_hitmap> tpc_view_hitmap;
  typedef std::map<unsigned int, tpc_view_hitmap> cryo_tpc_view_hitmap;

  class TrackShowerHits : public art::EDProducer {
  public:
    explicit TrackShowerHits(fhicl::ParameterSet const& p);

    TrackShowerHits(TrackShowerHits const&) = delete;
    TrackShowerHits(TrackShowerHits&&) = delete;
    TrackShowerHits& operator=(TrackShowerHits const&) = delete;
    TrackShowerHits& operator=(TrackShowerHits&&) = delete;

  private:
    void produce(art::Event& e) override;

    cryo_tpc_view_hitmap fHitMap;
    bool sortHits(const art::Event& evt);

    art::ServiceHandle<geo::Geometry const> fGeom;

    bool fHugeShowers, fShowersBySeg2D;

    tss::SimpleClustering fSimpleClustering;
    tss::Segmentation2D fSegmentation2D;

    std::string fHitModuleLabel;
  };
  // ------------------------------------------------------

  TrackShowerHits::TrackShowerHits(fhicl::ParameterSet const& p)
    : EDProducer{p}, fSegmentation2D(p.get<fhicl::ParameterSet>("Segmentation2DAlg"))
  {
    fHitModuleLabel = p.get<std::string>("HitModuleLabel");
    fHugeShowers = p.get<bool>("FindHugeShowers");
    fShowersBySeg2D = p.get<bool>("FindMoreShowers");

    produces<std::vector<recob::Cluster>>();
    produces<art::Assns<recob::Cluster, recob::Hit>>();
  }
  // ------------------------------------------------------

  bool
  TrackShowerHits::sortHits(const art::Event& evt)
  {
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    fHitMap.clear();

    art::Handle<std::vector<recob::Hit>> hitListHandle;
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if (evt.getByLabel(fHitModuleLabel, hitListHandle)) {
      art::fill_ptr_vector(hitlist, hitListHandle);

      unsigned int cryo, tpc, view;
      for (auto const& h : hitlist) {
        cryo = h->WireID().Cryostat;
        tpc = h->WireID().TPC;
        view = h->WireID().Plane;

        fHitMap[cryo][tpc][view].emplace_back(detProp, h);
      }
      return true;
    }
    else
      return false;
  }
  // ------------------------------------------------------

  void
  TrackShowerHits::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<recob::Cluster>> clusters(new std::vector<recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>> clu2hit(
      new art::Assns<recob::Cluster, recob::Hit>);

    if (sortHits(evt)) {
      unsigned int cidx = 0;
      const unsigned int emTag = 0x10000;

      for (auto tpc_iter = fGeom->begin_TPC_id(); tpc_iter != fGeom->end_TPC_id(); tpc_iter++) {
        for (const auto& v : fHitMap[tpc_iter->Cryostat][tpc_iter->TPC]) {
          auto cls = fSimpleClustering.run(v.second);

          if (fHugeShowers) {
            mf::LogVerbatim("TrackShowerHits") << "Find huge EM showers (cores).";

            int c = 0, clsSize = cls.size();
            while (c < clsSize) {
              if (cls[c].hits().size() < 2) {
                c++;
                continue;
              }

              std::vector<const tss::Hit2D*> trks, ems;
              fSegmentation2D.splitHitsNaive(cls[c], trks, ems);
              cls.erase(cls.begin() + c);
              clsSize--;

              cls.emplace_back(Cluster2D(trks));
              cls.emplace_back(Cluster2D(ems));
              cls.back().tagEM(true);
            }
          }

          if (fShowersBySeg2D) {
            mf::LogVerbatim("TrackShowerHits") << "Find EM showers by density of vtxs.";

            int c = 0, clsSize = cls.size();
            while (c < clsSize) {
              if (cls[c].isEM() || (cls[c].hits().size() < 2)) {
                c++;
                continue;
              }

              auto segs = fSegmentation2D.run(cls[c]);

              for (const auto& s : segs)
                cls.emplace_back(Cluster2D(s));

              cls.erase(cls.begin() + c);
              clsSize--;
            }
          }

          for (auto& c : cls) {
            if (!c.hits().size()) continue; // skip 0-size clusters

            if (!c.isEM()) continue; // create clusters only for em parts now

            clusters->emplace_back(
              recob::Cluster(0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             0.0F,
                             c.hits().size(),
                             0.0F,
                             0.0F,
                             cidx + emTag,
                             (geo::View_t)c.hits().front()->View(),
                             c.hits().front()->Hit2DPtr()->WireID().planeID()));

            std::vector<art::Ptr<recob::Hit>> hits2d;
            hits2d.reserve(c.hits().size());

            for (auto h2d : c.hits())
              hits2d.push_back(h2d->Hit2DPtr());

            if (hits2d.size()) util::CreateAssn(*this, evt, *clusters, hits2d, *clu2hit);

            ++cidx;
          }
        }
      }
    }
    else
      mf::LogWarning("TrackShowerHits") << "Hits not found in the event.";

    evt.put(std::move(clusters));
    evt.put(std::move(clu2hit));
  }
  // ------------------------------------------------------

  DEFINE_ART_MODULE(TrackShowerHits)

}
