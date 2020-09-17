////////////////////////////////////////////////////////////////////////
// Class:       TrackKalmanCheater
// Module Type: producer
// File:        TrackKalmanCheater.h
//
// This class produces RecoBase/Track objects using KalmanFilterService.
// MC truth information is used to associate Hits used as input to
// the Kalman filter.
//
// Configuration parameters:
//
// Hist               - Histogram flag (generate histograms if true).
// UseClusterHits     - Use clustered hits if true, use all hits if false.
// HitModuleLabel     - Module label for unclustered Hits.
// ClusterModuleLabel - Module label for Clusters.
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// KalmanFilterAlg    - Parameter set for KalmanFilterAlg.
// SpacePointAlg      - Parmaeter set for space points.
//
////////////////////////////////////////////////////////////////////////

#include <deque>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoObjects/KHit.h"
#include "lardata/RecoObjects/KHitContainerWireX.h"
#include "lardata/RecoObjects/PropYZPlane.h"
#include "lardata/RecoObjects/SurfYZPlane.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/RecoAlg/KalmanFilterAlg.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TH1F.h"

namespace {
  bool
  accepted_particle(int apdg)
  {
    return apdg == 13 ||  // Muon
           apdg == 211 || // Charged pion
           apdg == 321 || // Charged kaon
           apdg == 2212;  // (Anti)proton
  }
}

namespace trkf {

  class TrackKalmanCheater : public art::EDProducer {
  public:
    explicit TrackKalmanCheater(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& e) override;
    void beginJob() override;
    void endJob() override;

    // Fcl parameters.

    bool fHist;                      ///< Make histograms.
    KalmanFilterAlg fKFAlg;          ///< Kalman filter algorithm.
    SpacePointAlg fSpacePointAlg;    ///< Space point algorithm.
    bool fUseClusterHits;            ///< Use cluster hits or all hits?
    std::string fHitModuleLabel;     ///< Unclustered Hits.
    std::string fClusterModuleLabel; ///< Clustered Hits.
    double fMaxTcut;                 ///< Maximum delta ray energy in MeV for restricted dE/dx.

    // Histograms.

    TH1F* fHIncChisq{nullptr}; ///< Incremental chisquare.
    TH1F* fHPull{nullptr};     ///< Hit pull.

    // Statistics.

    int fNumEvent{0}; ///< Number of events seen.
    int fNumTrack{0}; ///< Number of tracks produced.
  };

  DEFINE_ART_MODULE(TrackKalmanCheater)

} // namespace trkf

//------------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// p - Fcl parameters.
///
trkf::TrackKalmanCheater::TrackKalmanCheater(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  , fHist(pset.get<bool>("Hist"))
  , fKFAlg(pset.get<fhicl::ParameterSet>("KalmanFilterAlg"))
  , fSpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"))
  , fUseClusterHits(pset.get<bool>("UseClusterHits"))
  , fHitModuleLabel{pset.get<std::string>("HitModuleLabel")}
  , fClusterModuleLabel{pset.get<std::string>("ClusterModuleLabel")}
  , fMaxTcut(pset.get<double>("MaxTcut"))
{
  produces<std::vector<recob::Track>>();
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Track, recob::Hit>>();
  produces<art::Assns<recob::Track, recob::SpacePoint>>();
  produces<art::Assns<recob::SpacePoint, recob::Hit>>();

  // Report.

  mf::LogInfo("TrackKalmanCheater")
    << "TrackKalmanCheater configured with the following parameters:\n"
    << "  UseClusterHits = " << fUseClusterHits << "\n"
    << "  HitModuleLabel = " << fHitModuleLabel << "\n"
    << "  ClusterModuleLabel = " << fClusterModuleLabel;
}

//------------------------------------------------------------------------------
/// Begin job method.
void
trkf::TrackKalmanCheater::beginJob()
{
  if (fHist) {

    // Book histograms.

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("hitkalman", "TrackKalmanCheater histograms");

    fHIncChisq = dir.make<TH1F>("IncChisq", "Incremental Chisquare", 100, 0., 20.);
    fHPull = dir.make<TH1F>("Pull", "Hit Pull", 100, -10., 10.);
  }
}

//------------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// e - Art event.
///
/// This method extracts Hit from the event and produces and adds
/// Track objects.
///
void
trkf::TrackKalmanCheater::produce(art::Event& evt)
{
  ++fNumEvent;

  // Make a collection of tracks, plus associations, that will
  // eventually be inserted into the event.

  std::unique_ptr<std::vector<recob::Track>> tracks(new std::vector<recob::Track>);
  std::unique_ptr<art::Assns<recob::Track, recob::Hit>> th_assn(
    new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint>> tsp_assn(
    new art::Assns<recob::Track, recob::SpacePoint>);

  // Make a collection of space points, plus associations, that will
  // be inserted into the event.

  std::unique_ptr<std::vector<recob::SpacePoint>> spts(new std::vector<recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit>> sph_assn(
    new art::Assns<recob::SpacePoint, recob::Hit>);

  // Make a collection of KGTracks where we will save our results.

  std::deque<KGTrack> kalman_tracks;

  // Get Services.

  art::ServiceHandle<geo::Geometry const> geom;
  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  // Reset space point algorithm.

  fSpacePointAlg.clearHitMap();

  // Get Hits.

  art::PtrVector<recob::Hit> hits;

  if (fUseClusterHits) {

    // Get clusters.

    art::Handle<std::vector<recob::Cluster>> clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    // Get hits from all clusters.
    art::FindManyP<recob::Hit> fm(clusterh, evt, fClusterModuleLabel);

    if (clusterh.isValid()) {
      int nclus = clusterh->size();

      for (int i = 0; i < nclus; ++i) {
        art::Ptr<recob::Cluster> pclus(clusterh, i);
        std::vector<art::Ptr<recob::Hit>> clushits = fm.at(i);
        int nhits = clushits.size();
        hits.reserve(hits.size() + nhits);

        for (std::vector<art::Ptr<recob::Hit>>::const_iterator ihit = clushits.begin();
             ihit != clushits.end();
             ++ihit) {
          hits.push_back(*ihit);
        }
      }
    }
  }
  else {

    // Get unclustered hits.

    art::Handle<std::vector<recob::Hit>> hith;
    evt.getByLabel(fHitModuleLabel, hith);
    if (hith.isValid()) {
      int nhits = hith->size();
      hits.reserve(nhits);

      for (int i = 0; i < nhits; ++i)
        hits.push_back(art::Ptr<recob::Hit>(hith, i));
    }
  }

  // Sort hits into separate PtrVectors based on track id.

  std::map<int, art::PtrVector<recob::Hit>> hitmap;

  // Loop over hits.

  for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end(); ++ihit) {
    //const recob::Hit& hit = **ihit;

    // Get track ids for this hit.

    std::vector<sim::TrackIDE> tids = bt_serv->HitToTrackIDEs(clockData, *ihit);

    // Loop over track ids.

    for (std::vector<sim::TrackIDE>::const_iterator itid = tids.begin(); itid != tids.end();
         ++itid) {
      int trackID = itid->trackID;

      // Add hit to PtrVector corresponding to this track id.

      hitmap[trackID].push_back(*ihit);
    }
  }

  // Extract geant mc particles.

  sim::ParticleList const& plist = pi_serv->ParticleList();

  auto const det_prop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

  PropYZPlane const prop{det_prop, fMaxTcut, true};

  // Loop over geant particles.

  {
    // use mf::LogDebug instead of MF_LOG_DEBUG because we reuse it in many lines;
    // insertions are protected by mf::isDebugEnabled()
    mf::LogDebug log("TrackKalmanCheater");
    for (sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      const simb::MCParticle* part = (*ipart).second;
      if (!part) throw cet::exception("TrackKalmanCheater") << "no particle!\n";
      int pdg = part->PdgCode();

      // Ignore everything except stable charged nonshowering particles.

      int apdg = std::abs(pdg);
      if (not accepted_particle(apdg)) { continue; }

      int trackid = part->TrackId();
      int stat = part->StatusCode();
      int nhit = 0;
      if (hitmap.count(trackid) != 0) nhit = hitmap[trackid].size();
      const TLorentzVector& pos = part->Position();
      const TLorentzVector& mom = part->Momentum();

      if (mf::isDebugEnabled()) {
        log << "Trackid=" << trackid << ", pdgid=" << pdg << ", status=" << stat
            << ", NHit=" << nhit << "\n"
            << "  x = " << pos.X() << ", y = " << pos.Y() << ", z = " << pos.Z() << "\n"
            << "  px= " << mom.Px() << ", py= " << mom.Py() << ", pz= " << mom.Pz() << "\n";
      } // if debugging

      double x = pos.X();
      double y = pos.Y();
      double z = pos.Z();
      double px = mom.Px();
      double py = mom.Py();
      double pz = mom.Pz();
      double p = std::hypot(px, py, pz);

      if (nhit > 0 && pz != 0.) {
        const art::PtrVector<recob::Hit>& trackhits = hitmap[trackid];
        if (trackhits.empty()) throw cet::exception("TrackKalmanCheater") << "No hits in track\n";

        // Make a seed track (KTrack).

        std::shared_ptr<const Surface> psurf(new SurfYZPlane(0., y, z, 0.));
        TrackVector vec(5);
        vec(0) = x;
        vec(1) = 0.;
        vec(2) = px / pz;
        vec(3) = py / pz;
        vec(4) = 1. / p;
        Surface::TrackDirection dir = (pz > 0. ? Surface::FORWARD : Surface::BACKWARD);
        KTrack trk(psurf, vec, dir, pdg);

        // Fill KHitContainer with hits.

        KHitContainerWireX cont;
        cont.fill(det_prop, trackhits, -1);

        // Count hits in each plane.  Set the preferred plane to be
        // the one with the most hits.

        std::vector<unsigned int> planehits(3, 0);
        for (art::PtrVector<recob::Hit>::const_iterator ih = trackhits.begin();
             ih != trackhits.end();
             ++ih) {
          const recob::Hit& hit = **ih;
          unsigned int plane = hit.WireID().Plane;

          if (plane >= planehits.size()) {
            throw cet::exception("TrackKalmanCheater") << "plane " << plane << "...\n";
          }
          ++planehits[plane];
        }
        unsigned int prefplane = 0;
        for (unsigned int i = 0; i < planehits.size(); ++i) {
          if (planehits[i] > planehits[prefplane]) prefplane = i;
        }
        if (mf::isDebugEnabled()) log << "Preferred plane = " << prefplane << "\n";

        // Build and smooth track.

        KGTrack trg(prefplane);
        fKFAlg.setPlane(prefplane);
        if (fKFAlg.buildTrack(trk, trg, prop, Propagator::FORWARD, cont, false)) {
          if (fKFAlg.smoothTrackIter(5, trg, prop)) {
            KETrack tremom;
            if (fKFAlg.fitMomentum(trg, prop, tremom)) { fKFAlg.updateMomentum(tremom, prop, trg); }
            ++fNumTrack;
            kalman_tracks.push_back(trg);
          }
        }
      }
    }
  }

  // Fill histograms.

  if (fHist) {

    // First loop over tracks.

    for (auto const& trg : kalman_tracks) {

      // Loop over measurements in this track.

      const std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();
      for (std::multimap<double, KHitTrack>::const_iterator ih = trackmap.begin();
           ih != trackmap.end();
           ++ih) {
        const KHitTrack& trh = (*ih).second;
        const std::shared_ptr<const KHitBase>& hit = trh.getHit();
        double chisq = hit->getChisq();
        fHIncChisq->Fill(chisq);
        const KHit<1>* ph1 = dynamic_cast<const KHit<1>*>(&*hit);
        if (ph1 != 0) {
          double pull = ph1->getResVector()(0) / std::sqrt(ph1->getResError()(0, 0));
          fHPull->Fill(pull);
        }
      }
    }
  }

  // Process Kalman filter tracks into persistent objects.

  tracks->reserve(kalman_tracks.size());
  for (auto const& kalman_track : kalman_tracks) {
    tracks->push_back(recob::Track());
    kalman_track.fillTrack(det_prop, tracks->back(), tracks->size() - 1);

    // Make Track to Hit associations.

    art::PtrVector<recob::Hit> trhits;
    std::vector<unsigned int> hittpindex;
    kalman_track.fillHits(hits, hittpindex);
    util::CreateAssn(evt, *tracks, trhits, *th_assn, tracks->size() - 1);

    // Make space points from this track.

    int nspt = spts->size();
    fSpacePointAlg.fillSpacePoints(det_prop, *spts, kalman_track.TrackMap());

    // Associate newly created space points with hits.
    // Also associate track with newly created space points.

    art::PtrVector<recob::SpacePoint> sptvec;

    // Loop over newly created space points.

    for (unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
      const recob::SpacePoint& spt = (*spts)[ispt];
      art::ProductID sptid = evt.getProductID<std::vector<recob::SpacePoint>>();
      art::Ptr<recob::SpacePoint> sptptr(sptid, ispt, evt.productGetter(sptid));
      sptvec.push_back(sptptr);

      // Make space point to hit associations.

      const art::PtrVector<recob::Hit>& sphits = fSpacePointAlg.getAssociatedHits(spt);
      util::CreateAssn(evt, *spts, sphits, *sph_assn, ispt);
    }

    // Make track to space point associations.

    util::CreateAssn(evt, *tracks, sptvec, *tsp_assn, tracks->size() - 1);
  }

  // Add tracks and associations to event.

  evt.put(std::move(tracks));
  evt.put(std::move(spts));
  evt.put(std::move(th_assn));
  evt.put(std::move(tsp_assn));
  evt.put(std::move(sph_assn));
}

//------------------------------------------------------------------------------
/// End job method.
void
trkf::TrackKalmanCheater::endJob()
{
  mf::LogInfo("TrackKalmanCheater") << "TrackKalmanCheater statistics:\n"
                                    << "  Number of events = " << fNumEvent << "\n"
                                    << "  Number of tracks created = " << fNumTrack;
}
