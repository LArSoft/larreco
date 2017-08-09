////////////////////////////////////////////////////////////////////////
// Class:       TrackProducerFromTrack
// Plugin Type: producer (art v2_07_03)
// File:        TrackProducerFromTrack_module.cc
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//
#include <memory>
//
#include "art/Utilities/make_tool.h" 
#include "larreco/TrackFinder/TrackMaker.h"
//
  /**
   * @brief Produce a reco::Track collection, as a result of the fit of an existing recob::Track collection.
   *
   * This producer takes an input an existing recob::Track collection (and the associated hits) and re-fits them.
   * The mandatory output are: the resulting recob::Track collection and the associated hits. 
   * Optional outputs are recob::TrackFitHitInfo and recob::SpacePoint collections, plus the Assns of SpacePoints to Hits.
   * An option is provided to create SpacePoints from the TrajectoryPoints in the Track.
   * Note: SpacePoints should not be used and will be soon deprecated as their functionality is covered by TrajectoryPoints.
   * The fit is performed by an user-defined tool, which must inherit from larreco/TrackFinder/TrackMaker.
   */
//
//
class TrackProducerFromTrack : public art::EDProducer {
public:
  explicit TrackProducerFromTrack(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  //
  // Plugins should not be copied or assigned.
  TrackProducerFromTrack(TrackProducerFromTrack const &) = delete;
  TrackProducerFromTrack(TrackProducerFromTrack &&) = delete;
  TrackProducerFromTrack & operator = (TrackProducerFromTrack const &) = delete;
  TrackProducerFromTrack & operator = (TrackProducerFromTrack &&) = delete;
  // Required functions.
  void produce(art::Event & e) override;
private:
  std::unique_ptr<trkmkr::TrackMaker> trackMaker_;
  art::InputTag trackInputTag;
  bool doTrackFitHitInfo_;
  bool doSpacePoints_;
  bool spacePointsFromTrajP_;
};
//
TrackProducerFromTrack::TrackProducerFromTrack(fhicl::ParameterSet const & p)
  : trackMaker_{art::make_tool<trkmkr::TrackMaker>(p.get<fhicl::ParameterSet>("trackMaker"))}
  , trackInputTag{p.get<art::InputTag>("inputCollection")}
  , doTrackFitHitInfo_{p.get<bool>("doTrackFitHitInfo")}
  , doSpacePoints_{p.get<bool>("doSpacePoints")}
  , spacePointsFromTrajP_{p.get<bool>("spacePointsFromTrajP")}
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  if (doTrackFitHitInfo_) produces<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  if (doSpacePoints_) {
    produces<std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::Hit, recob::SpacePoint> >();
  }
}
//
void TrackProducerFromTrack::produce(art::Event & e)
{
  // Output collections
  auto outputTracks  = std::make_unique<std::vector<recob::Track> >();
  auto outputHits    = std::make_unique<art::Assns<recob::Track, recob::Hit> >();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  auto outputSpacePoints  = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputHitSpacePointAssn = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint> >();
  //
  // Product ids for assns
  auto const tid = getProductID<std::vector<recob::Track> >(e);
  auto const tidgetter = e.productGetter(tid);
  auto const spid = getProductID<std::vector<recob::SpacePoint> >(e);
  auto const spidgetter = e.productGetter(spid);
  //
  // Input from event
  art::ValidHandle<std::vector<recob::Track> > inputTracks = e.getValidHandle<std::vector<recob::Track> >(trackInputTag);
  auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(trackInputTag);
  //
  // Loop over tracks to refit
  for (unsigned int iTrack = 0; iTrack < inputTracks->size(); ++iTrack) {
    //
    // Get track and its hits
    const recob::Track& track = inputTracks->at(iTrack);
    //not computationally optimal, but preserves the order (unlike FindManyP)
    std::vector<art::Ptr<recob::Hit> > inHits;
    art::Ptr<recob::Track> ptrack(inputTracks, iTrack);
    for (auto it = tkHitsAssn.begin(); it!=tkHitsAssn.end(); ++it) {
      if (it->first == ptrack) inHits.push_back(it->second);
      else if (inHits.size()>0) break;
    }
    //
    // Declare output objects
    recob::Track outTrack;
    art::PtrVector<recob::Hit> outHits;
    std::vector<recob::TrackFitHitInfo> outTrackFitHitInfos;
    std::vector<recob::SpacePoint> outSpacePoints;
    art::Assns<recob::Hit, recob::SpacePoint> outHitSpacePointAssn;
    //
    // Invoke tool to fit track and fill output objects
    trackMaker_->makeTrack(track, inHits, outTrack, outHits, outTrackFitHitInfos, outSpacePoints, outHitSpacePointAssn, e);
    //
    // Check that the requirement Nhits == Npoints is satisfied
    // We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
    if (outTrack.NumberTrajectoryPoints()!=outHits.size()) {
      throw std::runtime_error("TrackProducerFromTrack_module.cc - produced recob::Track required to have 1-1 correspondance between hits and points.");
    }
    //
    // Fill output collections, including Assns
    outputTracks->emplace_back(std::move(outTrack));
    art::Ptr<recob::Track> aptr(tid, outputTracks->size()-1, tidgetter);
    unsigned int ip = 0;
    for (auto const& trhit: outHits) {
      outputHits->addSingle(aptr, trhit);
      //
      if (spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
	auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
	const double fXYZ[3] = {tp.X(),tp.Y(),tp.Z()};
	const double fErrXYZ[6] = {0};
	recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
	outputSpacePoints->emplace_back(std::move(sp));
	art::Ptr<recob::SpacePoint> apsp(spid, outputSpacePoints->size()-1, spidgetter);
	outputHitSpacePointAssn->addSingle(trhit, apsp);
      }
      ip++;
    }
    if (doSpacePoints_ && !spacePointsFromTrajP_) {
      outputSpacePoints->insert(outputSpacePoints->end(),outSpacePoints.begin(),outSpacePoints.end());
      for (auto it = outHitSpacePointAssn.begin(); it!=outHitSpacePointAssn.end(); ++it ) {
	outputHitSpacePointAssn->addSingle(it->first,it->second);
      } 
    }
    if (doTrackFitHitInfo_) {
      outputHitInfo->emplace_back(std::move(outTrackFitHitInfos));
    }
  }
  //
  // Put collections in the event
  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));
  if (doTrackFitHitInfo_) {
    e.put(std::move(outputHitInfo));
  }
  if (doSpacePoints_) {
    e.put(std::move(outputSpacePoints));
    e.put(std::move(outputHitSpacePointAssn));
  }
}
//
DEFINE_ART_MODULE(TrackProducerFromTrack)
