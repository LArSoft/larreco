#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "art/Utilities/make_tool.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"
#include "larreco/TrackFinder/TrackMaker.h"

/**
   * @file  larreco/TrackFinder/TrackProducerFromTrack_module.cc
   * @class TrackProducerFromTrack
   *
   * @brief Produce a reco::Track collection, as a result of the fit of an existing recob::Track collection.
   *
   * This producer takes an input an existing recob::Track collection (and the associated hits) and re-fits them.
   * The mandatory output are: the resulting recob::Track collection and the associated hits.
   * Optional outputs are recob::TrackFitHitInfo and recob::SpacePoint collections, plus the Assns of SpacePoints to Hits.
   * An option is provided to create SpacePoints from the TrajectoryPoints in the Track.
   * The fit is performed by an user-defined tool, which must inherit from larreco/TrackFinder/TrackMaker.
   *
   * Parameters: trackMaker (fhicl::ParameterSet for the trkmkr::TrackMaker tool used to do the fit), inputCollection (art::InputTag of the input recob::Track collection),
   * doTrackFitHitInfo (bool to decide whether to produce recob::TrackFitHitInfo's), doSpacePoints (bool to decide whether to produce recob::SpacePoint's), and
   * spacePointsFromTrajP (bool to decide whether the produced recob::SpacePoint's are taken from the recob::tracking::TrajectoryPoint_t's of the fitted recob::Track).
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */
//
//
class TrackProducerFromTrack : public art::EDProducer {
public:
  explicit TrackProducerFromTrack(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  TrackProducerFromTrack(TrackProducerFromTrack const&) = delete;
  TrackProducerFromTrack(TrackProducerFromTrack&&) = delete;
  TrackProducerFromTrack& operator=(TrackProducerFromTrack const&) = delete;
  TrackProducerFromTrack& operator=(TrackProducerFromTrack&&) = delete;

private:
  void produce(art::Event& e) override;
  std::unique_ptr<trkmkr::TrackMaker> trackMaker_;
  art::InputTag trackInputTag;
  bool doTrackFitHitInfo_;
  bool doSpacePoints_;
  bool spacePointsFromTrajP_;
};

TrackProducerFromTrack::TrackProducerFromTrack(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , trackMaker_{art::make_tool<trkmkr::TrackMaker>(p.get<fhicl::ParameterSet>("trackMaker"))}
  , trackInputTag{p.get<art::InputTag>("inputCollection")}
  , doTrackFitHitInfo_{p.get<bool>("doTrackFitHitInfo")}
  , doSpacePoints_{p.get<bool>("doSpacePoints")}
  , spacePointsFromTrajP_{p.get<bool>("spacePointsFromTrajP")}
{
  produces<std::vector<recob::Track>>();
  produces<art::Assns<recob::Track, recob::Hit>>();
  if (doTrackFitHitInfo_) produces<std::vector<std::vector<recob::TrackFitHitInfo>>>();
  if (doSpacePoints_) {
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  }
}

void
TrackProducerFromTrack::produce(art::Event& e)
{
  // Output collections
  auto outputTracks = std::make_unique<std::vector<recob::Track>>();
  auto outputHits = std::make_unique<art::Assns<recob::Track, recob::Hit>>();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo>>>();
  auto outputSpacePoints = std::make_unique<std::vector<recob::SpacePoint>>();
  auto outputHitSpacePointAssn = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();

  // PtrMakers for Assns
  art::PtrMaker<recob::Track> trackPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> spacePointPtrMaker(e);

  // Input from event
  auto const inputTracks = e.getValidHandle<std::vector<recob::Track>>(trackInputTag);
  auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit>>(trackInputTag);
  const auto& tracksWithHits = util::associated_groups(tkHitsAssn);

  // Initialize tool for this event
  trackMaker_->initEvent(e);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  // Loop over tracks to refit
  unsigned int iTrack = 0;
  for (auto hitsRange : tracksWithHits) {

    // Get track and its hits
    art::Ptr<recob::Track> track(inputTracks, iTrack++);
    std::vector<art::Ptr<recob::Hit>> inHits;
    for (art::Ptr<recob::Hit> const& hit : hitsRange)
      inHits.push_back(hit);

    // Declare output objects
    recob::Track outTrack;
    std::vector<art::Ptr<recob::Hit>> outHits;
    trkmkr::OptionalOutputs optionals;
    if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
    if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();

    // Invoke tool to fit track and fill output objects
    bool fitok = trackMaker_->makeTrack(detProp, track, inHits, outTrack, outHits, optionals);
    if (!fitok) continue;

    // Check that the requirement Nhits == Npoints is satisfied
    // We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
    if (outTrack.NumberTrajectoryPoints() != outHits.size()) {
      throw cet::exception("TrackProducerFromTrack")
        << "Produced recob::Track required to have 1-1 correspondance between hits and "
           "points.\noutTrack.NumberTrajectoryPoints()="
        << outTrack.NumberTrajectoryPoints() << " outHits.size()=" << outHits.size() << "\n";
    }

    // Fill output collections, including Assns
    outputTracks->emplace_back(std::move(outTrack));
    const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size() - 1);
    unsigned int ip = 0;
    for (auto const& trhit : outHits) {
      outputHits->addSingle(aptr, trhit);

      if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
        auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
        const double fXYZ[3] = {tp.X(), tp.Y(), tp.Z()};
        const double fErrXYZ[6] = {0};
        recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
        outputSpacePoints->emplace_back(std::move(sp));
        const art::Ptr<recob::SpacePoint> apsp = spacePointPtrMaker(outputSpacePoints->size() - 1);
        outputHitSpacePointAssn->addSingle(trhit, apsp);
      }
      ip++;
    }
    if (doSpacePoints_ && !spacePointsFromTrajP_) {
      auto osp = optionals.spacePointHitPairs();
      for (auto it = osp.begin(); it != osp.end(); ++it) {
        outputSpacePoints->emplace_back(std::move(it->first));
        const art::Ptr<recob::SpacePoint> apsp = spacePointPtrMaker(outputSpacePoints->size() - 1);
        outputHitSpacePointAssn->addSingle(it->second, apsp);
      }
    }
    if (doTrackFitHitInfo_) { outputHitInfo->emplace_back(optionals.trackFitHitInfos()); }
  }

  // Put collections in the event
  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));
  if (doTrackFitHitInfo_) { e.put(std::move(outputHitInfo)); }
  if (doSpacePoints_) {
    e.put(std::move(outputSpacePoints));
    e.put(std::move(outputHitSpacePointAssn));
  }
}

DEFINE_ART_MODULE(TrackProducerFromTrack)
