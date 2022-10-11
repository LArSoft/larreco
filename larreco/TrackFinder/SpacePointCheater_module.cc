//
// Name: SpacePointCheater_module.cc
//
// Purpose: Module SpacePointCheater.
//
// Configuration parameters.
//
// ClusterModuleLabel;  // Cluster module label (e.g. "dbcluster").
// MinHits:             // Ignore clusters with fewer than this number of hits.
// ClusterAssns:        // If true, make associations between space points and clusters.
// SpacePointAlg:       // Configuration for SpacePointtAlg algoriithm.
//
// Created: 15-Dec-2011  H. Greenlee
//

#include <cassert>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/SpacePointAlg.h"

namespace trkf {

  class SpacePointCheater : public art::EDProducer {
  public:
    explicit SpacePointCheater(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;
    void endJob() override;

    // Fcl Attributes.

    SpacePointAlg fSptalg; // Algorithm object.
    std::string fClusterModuleLabel;
    unsigned int fMinHits; // Minimum number of hits per cluster.
    bool fClusterAssns;    // Make Cluster-SpacePoint associations.

    // Statistics.

    int fNumEvent; // Number of events.
    int fNumSpt2;  // Number of 2-view space points.
    int fNumSpt3;  // Number of 3-view space points.
  };

  DEFINE_ART_MODULE(SpacePointCheater)

  //----------------------------------------------------------------------------
  SpacePointCheater::SpacePointCheater(const fhicl::ParameterSet& pset)
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    : EDProducer{pset}
    , fSptalg(pset.get<fhicl::ParameterSet>("SpacePointAlg"))
    , fMinHits(0)
    , fClusterAssns(false)
    , fNumEvent(0)
    , fNumSpt2(0)
    , fNumSpt3(0)
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fMinHits = pset.get<unsigned int>("MinHits");
    fClusterAssns = pset.get<bool>("ClusterAssns");

    produces<std::vector<art::PtrVector<recob::SpacePoint>>>();
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
    if (fClusterAssns) produces<art::Assns<recob::SpacePoint, recob::Cluster>>();

    mf::LogInfo("SpacePointCheater")
      << "SpacePointCheater configured with the following parameters:\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  Minimum Hits per Cluster = " << fMinHits << "\n"
      << "  Cluster associations = " << fClusterAssns;
  }

  //----------------------------------------------------------------------------
  void SpacePointCheater::produce(art::Event& evt)
  //
  // Purpose: Produce method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Get Services.

    art::ServiceHandle<geo::Geometry const> geom;

    // Get clusters.

    art::Handle<std::vector<recob::Cluster>> clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointAlg).

    if (clusterh.isValid()) {

      // Make a collection of space points that will be inserted into the event.

      auto sptvecs = std::make_unique<std::vector<art::PtrVector<recob::SpacePoint>>>();
      auto spts = std::make_unique<std::vector<recob::SpacePoint>>();
      auto sphitassn = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();
      auto spclassn = std::make_unique<art::Assns<recob::SpacePoint, recob::Cluster>>();

      // Make a hit vector which will be used to store hits to be passed
      // to SpacePointAlg.

      art::PtrVector<recob::Hit> hits;
      art::FindManyP<recob::Hit> fm(clusterh, evt, fClusterModuleLabel);

      // Loop over first cluster.

      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

      int nclus = clusterh->size();
      for (int iclus = 0; iclus < nclus; ++iclus) {
        art::Ptr<recob::Cluster> piclus(clusterh, iclus);
        geo::View_t iview = piclus->View();

        std::vector<art::Ptr<recob::Hit>> ihits = fm.at(iclus);

        // Test first view.

        if (ihits.size() >= fMinHits &&
            ((iview == geo::kU && fSptalg.enableU()) || (iview == geo::kV && fSptalg.enableV()) ||
             (iview == geo::kZ && fSptalg.enableW()))) {

          // Store hits from first view into hit vector.

          unsigned int nihits = ihits.size();
          hits.clear();
          hits.reserve(nihits);
          for (std::vector<art::Ptr<recob::Hit>>::const_iterator i = ihits.begin();
               i != ihits.end();
               ++i)
            hits.push_back(*i);

          // Loop over second cluster.

          for (int jclus = 0; jclus < iclus; ++jclus) {
            art::Ptr<recob::Cluster> pjclus(clusterh, jclus);
            geo::View_t jview = pjclus->View();

            std::vector<art::Ptr<recob::Hit>> jhits = fm.at(jclus);

            // Test second view.

            if (jhits.size() >= fMinHits &&
                ((jview == geo::kU && fSptalg.enableU()) ||
                 (jview == geo::kV && fSptalg.enableV()) ||
                 (jview == geo::kZ && fSptalg.enableW())) &&
                jview != iview) {

              // Store hits from second view into hit vector.

              unsigned int njhits = jhits.size();
              assert(hits.size() >= nihits);
              //hits.resize(nihits);
              while (hits.size() > nihits)
                hits.pop_back();
              assert(hits.size() == nihits);
              hits.reserve(nihits + njhits);
              for (std::vector<art::Ptr<recob::Hit>>::const_iterator j = jhits.begin();
                   j != jhits.end();
                   ++j)
                hits.push_back(*j);

              // If two-view space points are allowed, make them here.

              if (fSptalg.minViews() <= 2) {
                std::vector<recob::SpacePoint> new_spts;
                fSptalg.makeMCTruthSpacePoints(clockData, detProp, hits, new_spts);

                // If we found some space points, insert them into the event.

                if (new_spts.size() > 0) {
                  fNumSpt2 += new_spts.size();
                  art::PtrVector<recob::Cluster> clusters;
                  clusters.reserve(2);
                  clusters.push_back(piclus);
                  clusters.push_back(pjclus);

                  // Insert newly found space points into event collection.

                  int nspt = spts->size();
                  spts->insert(spts->end(), new_spts.begin(), new_spts.end());

                  // Associate space points with hits and clusters.

                  art::PtrVector<recob::SpacePoint> sptvec;
                  for (unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
                    const recob::SpacePoint& spt = (*spts)[ispt];
                    const art::PtrVector<recob::Hit>& hits = fSptalg.getAssociatedHits(spt);
                    util::CreateAssn(evt, *spts, hits, *sphitassn, ispt);
                    if (fClusterAssns) util::CreateAssn(evt, *spts, clusters, *spclassn, ispt);

                    // make the PtrVector for this collection of space points
                    // Do not reproduce the following lines
                    // Contact brebel@fnal.gov if you think you need to reproduce these lines.
                    art::ProductID spid = evt.getProductID<std::vector<recob::SpacePoint>>();
                    art::Ptr<recob::SpacePoint> spptr(spid, ispt, evt.productGetter(spid));
                    sptvec.push_back(spptr);
                  }
                  sptvecs->push_back(sptvec);
                }
              }

              // Loop over third cluster.

              for (int kclus = 0; kclus < jclus; ++kclus) {
                art::Ptr<recob::Cluster> pkclus(clusterh, kclus);
                geo::View_t kview = pkclus->View();

                std::vector<art::Ptr<recob::Hit>> khits = fm.at(kclus);

                // Test third view.

                if (khits.size() >= fMinHits &&
                    ((kview == geo::kU && fSptalg.enableU()) ||
                     (kview == geo::kV && fSptalg.enableV()) ||
                     (kview == geo::kZ && fSptalg.enableW())) &&
                    kview != iview && kview != jview) {

                  // Store hits from third view into hit vector.

                  unsigned int nkhits = khits.size();
                  assert(hits.size() >= nihits + njhits);
                  //hits.resize(nihits + njhits);
                  while (hits.size() > nihits + njhits)
                    hits.pop_back();
                  assert(hits.size() == nihits + njhits);
                  hits.reserve(nihits + njhits + nkhits);
                  for (std::vector<art::Ptr<recob::Hit>>::const_iterator k = khits.begin();
                       k != khits.end();
                       ++k)
                    hits.push_back(*k);

                  // Make three-view space points.

                  std::vector<recob::SpacePoint> new_spts;
                  fSptalg.makeMCTruthSpacePoints(clockData, detProp, hits, new_spts);

                  // If we found some space points, insert them into the event.

                  if (new_spts.size() > 0) {
                    fNumSpt3 += new_spts.size();
                    art::PtrVector<recob::Cluster> clusters;
                    clusters.reserve(3);
                    clusters.push_back(piclus);
                    clusters.push_back(pjclus);
                    clusters.push_back(pkclus);

                    // Insert newly found space points into event collection.

                    int nspt = spts->size();
                    spts->insert(spts->end(), new_spts.begin(), new_spts.end());

                    // Associate space points with hits and clusters.

                    art::PtrVector<recob::SpacePoint> sptvec;
                    for (unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
                      const recob::SpacePoint& spt = (*spts)[ispt];
                      const art::PtrVector<recob::Hit>& hits = fSptalg.getAssociatedHits(spt);
                      util::CreateAssn(evt, *spts, hits, *sphitassn, ispt);
                      if (fClusterAssns) util::CreateAssn(evt, *spts, clusters, *spclassn, ispt);

                      // make the PtrVector for this collection of space points
                      // Do not reproduce the following lines
                      // Contact brebel@fnal.gov if you think you need to reproduce these lines.
                      art::ProductID spid = evt.getProductID<std::vector<recob::SpacePoint>>();
                      art::Ptr<recob::SpacePoint> spptr(spid, ispt, evt.productGetter(spid));
                      sptvec.push_back(spptr);
                    }
                    sptvecs->push_back(sptvec);
                  }
                }
              }
            }
          }
        }
      }

      // Add space points and associations to event.

      evt.put(std::move(spts));
      evt.put(std::move(sptvecs));
      evt.put(std::move(sphitassn));
      if (fClusterAssns) evt.put(std::move(spclassn));
    }
  }

  //----------------------------------------------------------------------------
  void SpacePointCheater::endJob()
  //
  // Purpose: Print summary.
  //
  {
    mf::LogInfo("SpacePointCheater")
      << "SpacePointCheater statistics:\n"
      << "  Number of events = " << fNumEvent << "\n"
      << "  Number of 2-view space points created = " << fNumSpt2 << "\n"
      << "  Number of 3-view space points created = " << fNumSpt3;
  }
}
