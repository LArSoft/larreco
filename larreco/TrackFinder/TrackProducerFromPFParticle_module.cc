#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//
#include <memory>
//
#include "art/Utilities/make_tool.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larreco/TrackFinder/TrackMaker.h"
//
/**
   * @file  larreco/TrackFinder/TrackProducerFromPFParticle_module.cc
   * @class TrackProducerFromPFParticle
   *
   * @brief Produce a reco::Track collection, as a result of the fit of an existing recob::PFParticle collection.
   *
   * This producer takes an input an existing recob::PFParticle collection and refits the associated tracks;
   * it can make track fits of the associated showers as well, but this is experimental - do it at your own risk.
   * The mandatory outputs are: the resulting recob::Track collection, the associated hits, and the association
   * between the input PFParticle and the output Track.
   * Optional outputs are recob::TrackFitHitInfo and recob::SpacePoint collections, plus the Assns of SpacePoints to Hits.
   * An option is provided to create SpacePoints from the TrajectoryPoints in the Track.
   * The fit is performed by an user-defined tool, which must inherit from larreco/TrackFinder/TrackMaker.
   *
   * Parameters: trackMaker (fhicl::ParameterSet for the trkmkr::TrackMaker tool used to do the fit), inputCollection (art::InputTag of the input recob::Track collection),
   * doTrackFitHitInfo (bool to decide whether to produce recob::TrackFitHitInfo's), doSpacePoints (bool to decide whether to produce recob::SpacePoint's),
   * spacePointsFromTrajP (bool to decide whether the produced recob::SpacePoint's are taken from the recob::tracking::TrajectoryPoint_t's of the fitted recob::Track),
   * trackFromPF (bool to decide whether to fit the recob::Track associated to the recob::PFParticle), and
   * showerFromPF (bool to decide whether to fit the recob::Shower associated to the recob::PFParticle - this option is intended to mitigate possible problems due to tracks being mis-identified as showers)
   * seedFromPF (bool to decide whether to fit the recob::PFParticle using the associated seed)
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */
//
//
class TrackProducerFromPFParticle : public art::EDProducer {
public:
  explicit TrackProducerFromPFParticle(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  //
  // Plugins should not be copied or assigned.
  TrackProducerFromPFParticle(TrackProducerFromPFParticle const&) = delete;
  TrackProducerFromPFParticle(TrackProducerFromPFParticle&&) = delete;
  TrackProducerFromPFParticle& operator=(TrackProducerFromPFParticle const&) = delete;
  TrackProducerFromPFParticle& operator=(TrackProducerFromPFParticle&&) = delete;

private:
  // Required functions.
  void produce(art::Event& e) override;
  std::unique_ptr<trkmkr::TrackMaker> trackMaker_;
  art::InputTag pfpInputTag;
  art::InputTag trkInputTag;
  art::InputTag shwInputTag;
  art::InputTag clsInputTag;
  bool doTrackFitHitInfo_;
  bool doSpacePoints_;
  bool spacePointsFromTrajP_;
  bool trackFromPF_;
  bool showerFromPF_;
  bool seedFromPF_;
};
//
TrackProducerFromPFParticle::TrackProducerFromPFParticle(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , trackMaker_{art::make_tool<trkmkr::TrackMaker>(p.get<fhicl::ParameterSet>("trackMaker"))}
  , pfpInputTag{p.get<art::InputTag>("inputCollection")}
  , doTrackFitHitInfo_{p.get<bool>("doTrackFitHitInfo")}
  , doSpacePoints_{p.get<bool>("doSpacePoints")}
  , spacePointsFromTrajP_{p.get<bool>("spacePointsFromTrajP")}
  , trackFromPF_{p.get<bool>("trackFromPF")}
  , showerFromPF_{p.get<bool>("showerFromPF")}
  , seedFromPF_{p.get<bool>("seedFromPF")}
{
  //
  if (p.has_key("trackInputTag"))
    trkInputTag = p.get<art::InputTag>("trackInputTag");
  else
    trkInputTag = pfpInputTag;
  if (p.has_key("showerInputTag"))
    shwInputTag = p.get<art::InputTag>("showerInputTag");
  else
    shwInputTag = pfpInputTag;
  if (p.has_key("clusterInputTag"))
    clsInputTag = p.get<art::InputTag>("clusterInputTag");
  else
    clsInputTag = pfpInputTag;
  produces<std::vector<recob::Track>>();
  produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  produces<art::Assns<recob::PFParticle, recob::Track>>();
  if (doTrackFitHitInfo_) produces<std::vector<std::vector<recob::TrackFitHitInfo>>>();
  if (doSpacePoints_) {
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  }
}
//
void
TrackProducerFromPFParticle::produce(art::Event& e)
{
  // Output collections
  auto outputTracks = std::make_unique<std::vector<recob::Track>>();
  auto outputHits = std::make_unique<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  auto outputPfpTAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Track>>();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo>>>();
  auto outputSpacePoints = std::make_unique<std::vector<recob::SpacePoint>>();
  auto outputHitSpacePointAssn = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
  //
  // PtrMakers for Assns
  art::PtrMaker<recob::Track> trackPtrMaker(e);
  art::PtrMaker<recob::SpacePoint>* spacePointPtrMaker = nullptr;
  if (doSpacePoints_) spacePointPtrMaker = new art::PtrMaker<recob::SpacePoint>(e);
  //
  // Input from event
  art::ValidHandle<std::vector<recob::PFParticle>> inputPfps =
    e.getValidHandle<std::vector<recob::PFParticle>>(pfpInputTag);
  std::unique_ptr<art::FindManyP<recob::Track>> assocTracks;
  art::Assns<recob::Track, recob::Hit> tkHitsAssn;
  std::unique_ptr<art::FindManyP<recob::Shower>> assocShowers;
  std::unique_ptr<art::FindManyP<recob::Seed>> assocSeeds;
  if (trackFromPF_) {
    assocTracks = std::unique_ptr<art::FindManyP<recob::Track>>(
      new art::FindManyP<recob::Track>(inputPfps, e, trkInputTag));
    tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit>>(trkInputTag);
  }
  if (showerFromPF_)
    assocShowers = std::unique_ptr<art::FindManyP<recob::Shower>>(
      new art::FindManyP<recob::Shower>(inputPfps, e, shwInputTag));
  if (seedFromPF_)
    assocSeeds = std::unique_ptr<art::FindManyP<recob::Seed>>(
      new art::FindManyP<recob::Seed>(inputPfps, e, pfpInputTag));
  const auto& trackHitsGroups = util::associated_groups(tkHitsAssn);
  //
  std::unique_ptr<art::FindManyP<recob::Cluster>> assocClusters =
    std::unique_ptr<art::FindManyP<recob::Cluster>>(
      new art::FindManyP<recob::Cluster>(inputPfps, e, pfpInputTag));
  auto const& clHitsAssn = *e.getValidHandle<art::Assns<recob::Cluster, recob::Hit>>(clsInputTag);
  const auto& clusterHitsGroups = util::associated_groups(clHitsAssn);
  //
  // Initialize tool for this event
  trackMaker_->initEvent(e);
  //
  // Loop over pfps to fit
  for (unsigned int iPfp = 0; iPfp < inputPfps->size(); ++iPfp) {
    //
    const art::Ptr<recob::PFParticle> pfp(inputPfps, iPfp);
    //
    if (trackFromPF_) {
      // Tracks associated to PFParticles
      const std::vector<art::Ptr<recob::Track>>& tracks = assocTracks->at(iPfp);
      // Loop over tracks to refit
      for (art::Ptr<recob::Track> const& track : tracks) {
        //
        // Get track and its hits
        std::vector<art::Ptr<recob::Hit>> inHits;
        decltype(auto) hitsRange = util::groupByIndex(trackHitsGroups, track.key());
        for (art::Ptr<recob::Hit> const& hit : hitsRange)
          inHits.push_back(hit);
        //
        // Declare output objects
        recob::Track outTrack;
        std::vector<art::Ptr<recob::Hit>> outHits;
        trkmkr::OptionalOutputs optionals;
        if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
        if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
        //
        // Invoke tool to fit track and fill output objects
        bool fitok = trackMaker_->makeTrack(track, inHits, outTrack, outHits, optionals);
        if (!fitok) continue;
        //
        // Check that the requirement Nhits == Npoints is satisfied
        // We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
        if (outTrack.NumberTrajectoryPoints() != outHits.size()) {
          throw cet::exception("TrackProducerFromPFParticle")
            << "Produced recob::Track required to have 1-1 correspondance between hits and "
               "points.\n";
        }
        //
        // Fill output collections, including Assns
        outputTracks->emplace_back(std::move(outTrack));
        const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size() - 1);
        outputPfpTAssn->addSingle(pfp, aptr);
        unsigned int ip = 0;
        for (auto const& trhit : outHits) {
          recob::TrackHitMeta metadata(
            outputTracks->back().HasValidPoint(ip) ? ip : std::numeric_limits<int>::max(),
            -std::numeric_limits<double>::max());
          outputHits->addSingle(aptr, trhit, metadata);
          //
          if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
            auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
            const double fXYZ[3] = {tp.X(), tp.Y(), tp.Z()};
            const double fErrXYZ[6] = {0};
            recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
            outputSpacePoints->emplace_back(std::move(sp));
            const art::Ptr<recob::SpacePoint> apsp =
              (*spacePointPtrMaker)(outputSpacePoints->size() - 1);
            outputHitSpacePointAssn->addSingle(trhit, apsp);
          }
          ip++;
        }
        if (doSpacePoints_ && !spacePointsFromTrajP_) {
          auto osp = optionals.spacePointHitPairs();
          for (auto it = osp.begin(); it != osp.end(); ++it) {
            outputSpacePoints->emplace_back(std::move(it->first));
            const art::Ptr<recob::SpacePoint> apsp =
              (*spacePointPtrMaker)(outputSpacePoints->size() - 1);
            outputHitSpacePointAssn->addSingle(it->second, apsp);
          }
        }
        if (doTrackFitHitInfo_) { outputHitInfo->emplace_back(optionals.trackFitHitInfos()); }
      }
    }
    //
    if (showerFromPF_) {
      //
      // Showers associated to PFParticles
      const std::vector<art::Ptr<recob::Shower>>& showers = assocShowers->at(iPfp);
      // if there is more than one shower the logic below to get the hits does not work! this works, at least for uboone
      if (showers.size() != 1) continue;
      //
      // Get hits for shower (through the chain pfp->clusters->hits)
      std::vector<art::Ptr<recob::Hit>> inHits;
      const std::vector<art::Ptr<recob::Cluster>> clustersRange = assocClusters->at(iPfp);
      for (art::Ptr<recob::Cluster> const& cluster : clustersRange) {
        // for hits we use groupByIndex since it preserves the order (and we can use it since each cluster must have associated hits)
        decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, cluster.key());
        for (art::Ptr<recob::Hit> const& hit : hitsRange)
          inHits.push_back(hit);
      }
      // Loop over showers to refit (should be only one)
      for (unsigned int iShower = 0; iShower < showers.size(); ++iShower) {
        //
        // Get the shower and convert/hack it into a trajectory so that the fit is initialized
        art::Ptr<recob::Shower> shower = showers[iShower];
        recob::tracking::Point_t pos(
          shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
        recob::tracking::Vector_t dir(
          shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());
        float mom = 1.;
        if (shower->Energy().size() == 3) mom = shower->Energy()[2] * 0.001;
        std::vector<recob::tracking::Point_t> p;
        std::vector<recob::tracking::Vector_t> d;
        for (unsigned int i = 0; i < inHits.size(); ++i) {
          p.push_back(pos);
          d.push_back(mom * dir);
        }
        recob::TrackTrajectory traj(
          std::move(p), std::move(d), recob::TrackTrajectory::Flags_t(p.size()), false);
        //
        // Declare output objects
        recob::Track outTrack;
        std::vector<art::Ptr<recob::Hit>> outHits;
        trkmkr::OptionalOutputs optionals;
        if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
        if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
        //
        // Invoke tool to fit track and fill output objects
        bool fitok = trackMaker_->makeTrack(traj, iPfp, inHits, outTrack, outHits, optionals);
        if (!fitok) continue;
        //
        // Check that the requirement Nhits == Npoints is satisfied
        // We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
        if (outTrack.NumberTrajectoryPoints() != outHits.size()) {
          throw cet::exception("TrackProducerFromPFParticle")
            << "Produced recob::Track required to have 1-1 correspondance between hits and "
               "points.\n";
        }
        //
        // Fill output collections, including Assns
        outputTracks->emplace_back(std::move(outTrack));
        const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size() - 1);
        outputPfpTAssn->addSingle(pfp, aptr);
        unsigned int ip = 0;
        for (auto const& trhit : outHits) {
          recob::TrackHitMeta metadata(
            outputTracks->back().HasValidPoint(ip) ? ip : std::numeric_limits<int>::max(),
            -std::numeric_limits<double>::max());
          outputHits->addSingle(aptr, trhit, metadata);
          //
          if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
            auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
            const double fXYZ[3] = {tp.X(), tp.Y(), tp.Z()};
            const double fErrXYZ[6] = {0};
            recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
            outputSpacePoints->emplace_back(std::move(sp));
            const art::Ptr<recob::SpacePoint> apsp =
              (*spacePointPtrMaker)(outputSpacePoints->size() - 1);
            outputHitSpacePointAssn->addSingle(trhit, apsp);
          }
          ip++;
        }
        if (doSpacePoints_ && !spacePointsFromTrajP_) {
          auto osp = optionals.spacePointHitPairs();
          for (auto it = osp.begin(); it != osp.end(); ++it) {
            outputSpacePoints->emplace_back(std::move(it->first));
            const art::Ptr<recob::SpacePoint> apsp =
              (*spacePointPtrMaker)(outputSpacePoints->size() - 1);
            outputHitSpacePointAssn->addSingle(it->second, apsp);
          }
        }
        if (doTrackFitHitInfo_) { outputHitInfo->emplace_back(optionals.trackFitHitInfos()); }
      }
    }
    //
    //
    if (seedFromPF_) {
      //
      // Seeds associated to PFParticles
      const std::vector<art::Ptr<recob::Seed>>& seeds = assocSeeds->at(iPfp);
      // if there is more than one seed the logic below to get the hits does not work! this works, at least for uboone
      if (seeds.size() != 1) continue;
      //
      // Get hits for pfp (through the chain pfp->clusters->hits)
      std::vector<art::Ptr<recob::Hit>> inHits;
      const std::vector<art::Ptr<recob::Cluster>> clustersRange = assocClusters->at(iPfp);
      for (art::Ptr<recob::Cluster> const& cluster : clustersRange) {
        // for hits we use groupByIndex since it preserves the order (and we can use it since each cluster must have associated hits)
        decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, cluster.key());
        for (art::Ptr<recob::Hit> const& hit : hitsRange)
          inHits.push_back(hit);
      }
      if (inHits.size() < 4) continue;
      // Loop over seeds should be only one)
      for (unsigned int iS = 0; iS < seeds.size(); ++iS) {
        //
        // Get the seed and convert/hack it into a trajectory so that the fit is initialized
        art::Ptr<recob::Seed> seed = seeds[iS];
        double p0[3], pe[3];
        seed->GetPoint(p0, pe);
        double d0[3], de[3];
        seed->GetDirection(d0, de);
        recob::tracking::Point_t pos(p0[0], p0[1], p0[2]);
        recob::tracking::Vector_t dir(d0[0], d0[1], d0[2]);
        std::vector<recob::tracking::Point_t> p;
        std::vector<recob::tracking::Vector_t> d;
        for (unsigned int i = 0; i < inHits.size(); ++i) {
          p.push_back(pos);
          d.push_back(dir);
        }
        recob::TrackTrajectory traj(
          std::move(p), std::move(d), recob::TrackTrajectory::Flags_t(p.size()), false);
        //
        // Declare output objects
        recob::Track outTrack;
        std::vector<art::Ptr<recob::Hit>> outHits;
        trkmkr::OptionalOutputs optionals;
        if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
        if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
        //
        // Invoke tool to fit track and fill output objects
        bool fitok = trackMaker_->makeTrack(traj, iPfp, inHits, outTrack, outHits, optionals);
        if (!fitok) continue;
        //
        // Check that the requirement Nhits == Npoints is satisfied
        // We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
        if (outTrack.NumberTrajectoryPoints() != outHits.size()) {
          throw cet::exception("TrackProducerFromPFParticle")
            << "Produced recob::Track required to have 1-1 correspondance between hits and "
               "points.\n";
        }
        //
        // Fill output collections, including Assns
        outputTracks->emplace_back(std::move(outTrack));
        const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size() - 1);
        outputPfpTAssn->addSingle(pfp, aptr);
        unsigned int ip = 0;
        for (auto const& trhit : outHits) {
          recob::TrackHitMeta metadata(
            outputTracks->back().HasValidPoint(ip) ? ip : std::numeric_limits<int>::max(),
            -std::numeric_limits<double>::max());
          outputHits->addSingle(aptr, trhit, metadata);
          //
          if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
            auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
            const double fXYZ[3] = {tp.X(), tp.Y(), tp.Z()};
            const double fErrXYZ[6] = {0};
            recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
            outputSpacePoints->emplace_back(std::move(sp));
            const art::Ptr<recob::SpacePoint> apsp =
              (*spacePointPtrMaker)(outputSpacePoints->size() - 1);
            outputHitSpacePointAssn->addSingle(trhit, apsp);
          }
          ip++;
        }
        if (doSpacePoints_ && !spacePointsFromTrajP_) {
          auto osp = optionals.spacePointHitPairs();
          for (auto it = osp.begin(); it != osp.end(); ++it) {
            outputSpacePoints->emplace_back(std::move(it->first));
            const art::Ptr<recob::SpacePoint> apsp =
              (*spacePointPtrMaker)(outputSpacePoints->size() - 1);
            outputHitSpacePointAssn->addSingle(it->second, apsp);
          }
        }
        if (doTrackFitHitInfo_) { outputHitInfo->emplace_back(optionals.trackFitHitInfos()); }
      }
    }
    //
  }
  //
  // Put collections in the event
  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));
  e.put(std::move(outputPfpTAssn));
  if (doTrackFitHitInfo_) { e.put(std::move(outputHitInfo)); }
  if (doSpacePoints_) {
    e.put(std::move(outputSpacePoints));
    e.put(std::move(outputHitSpacePointAssn));
  }
  if (doSpacePoints_) delete spacePointPtrMaker;
}
//
DEFINE_ART_MODULE(TrackProducerFromPFParticle)
