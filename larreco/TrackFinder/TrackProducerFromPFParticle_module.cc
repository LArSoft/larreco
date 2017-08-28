////////////////////////////////////////////////////////////////////////
// Class:       TrackProducerFromPFParticle
// Plugin Type: producer (art v2_07_03)
// File:        TrackProducerFromPFParticle_module.cc
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
#include "lardata/Utilities/PtrMaker.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"
#include "cetlib/exception.h"
//
#include <memory>
//
#include "art/Utilities/make_tool.h"
#include "larreco/TrackFinder/TrackMaker.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
//
  /**
   * @brief Produce a reco::Track collection, as a result of the fit of an existing recob::PFParticle collection.
   *
   * This producer takes an input an existing recob::PFParticle collection and refits the associated tracks;
   * it can make track fits of the associated showers as well, but this is experimental - do it at your own risk.
   * The mandatory outputs are: the resulting recob::Track collection, the associated hits, and the association
   * between the input PFParticle and the output Track.
   * Optional outputs are recob::TrackFitHitInfo and recob::SpacePoint collections, plus the Assns of SpacePoints to Hits.
   * An option is provided to create SpacePoints from the TrajectoryPoints in the Track.
   * Note: SpacePoints should not be used and will be soon deprecated as their functionality is covered by TrajectoryPoints.
   * The fit is performed by an user-defined tool, which must inherit from larreco/TrackFinder/TrackMaker.
   */
//
//
class TrackProducerFromPFParticle : public art::EDProducer {
public:
  explicit TrackProducerFromPFParticle(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  //
  // Plugins should not be copied or assigned.
  TrackProducerFromPFParticle(TrackProducerFromPFParticle const &) = delete;
  TrackProducerFromPFParticle(TrackProducerFromPFParticle &&) = delete;
  TrackProducerFromPFParticle & operator = (TrackProducerFromPFParticle const &) = delete;
  TrackProducerFromPFParticle & operator = (TrackProducerFromPFParticle &&) = delete;
  // Required functions.
  void produce(art::Event & e) override;
private:
  std::unique_ptr<trkmkr::TrackMaker> trackMaker_;
  art::InputTag pfpInputTag;
  art::InputTag trkInputTag;
  art::InputTag shwInputTag;
  bool doTrackFitHitInfo_;
  bool doSpacePoints_;
  bool spacePointsFromTrajP_;
  bool trackFromPF_;
  bool showerFromPF_;
};
//
TrackProducerFromPFParticle::TrackProducerFromPFParticle(fhicl::ParameterSet const & p)
  : trackMaker_{art::make_tool<trkmkr::TrackMaker>(p.get<fhicl::ParameterSet>("trackMaker"))}
  , pfpInputTag{p.get<art::InputTag>("inputCollection")}
  , doTrackFitHitInfo_{p.get<bool>("doTrackFitHitInfo")}
  , doSpacePoints_{p.get<bool>("doSpacePoints")}
  , spacePointsFromTrajP_{p.get<bool>("spacePointsFromTrajP")}
  , trackFromPF_{p.get<bool>("trackFromPF")}
  , showerFromPF_{p.get<bool>("showerFromPF")}
{
  //
  if (p.has_key("trackInputTag")) trkInputTag = p.get<art::InputTag>("trackInputTag");
  else trkInputTag = pfpInputTag;
  if (p.has_key("showerInputTag")) shwInputTag = p.get<art::InputTag>("showerInputTag");
  else shwInputTag = pfpInputTag;
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  produces<art::Assns<recob::PFParticle, recob::Track> >();
  if (doTrackFitHitInfo_) produces<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  if (doSpacePoints_) {
    produces<std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::Hit, recob::SpacePoint> >();
  }
}
//
void TrackProducerFromPFParticle::produce(art::Event & e)
{
  // Output collections
  auto outputTracks  = std::make_unique<std::vector<recob::Track> >();
  auto outputHits    = std::make_unique<art::Assns<recob::Track, recob::Hit> >();
  auto outputPfpTAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Track> >();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  auto outputSpacePoints  = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputHitSpacePointAssn = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint> >();
  //
  // PtrMakers for Assns
  lar::PtrMaker<recob::Track> trackPtrMaker(e, *this);
  lar::PtrMaker<recob::SpacePoint> spacePointPtrMaker(e, *this);
  //
  // Input from event
  art::ValidHandle<std::vector<recob::PFParticle> > inputPfps = e.getValidHandle<std::vector<recob::PFParticle> >(pfpInputTag);
  const auto assocTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPfps, e, pfpInputTag));
  const auto assocShowers = std::unique_ptr<art::FindManyP<recob::Shower> >(new art::FindManyP<recob::Shower>(inputPfps, e, pfpInputTag));
  auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(trkInputTag);
  const auto& trackHitsGroups = util::associated_groups(tkHitsAssn);
  //
  auto const& pfClustersAssn = *e.getValidHandle<art::Assns<recob::PFParticle, recob::Cluster> >(pfpInputTag);
  const auto& pfpClusterGroups = util::associated_groups(pfClustersAssn);
  auto const& clHitsAssn = *e.getValidHandle<art::Assns<recob::Cluster, recob::Hit> >(shwInputTag);
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
      const std::vector<art::Ptr<recob::Track> >& tracks = assocTracks->at(iPfp);
      // Loop over tracks to refit
      for (art::Ptr<recob::Track> const& track: tracks) {
	//
	// Get track and its hits
	std::vector<art::Ptr<recob::Hit> > inHits;
	decltype(auto) hitsRange = util::groupByIndex(trackHitsGroups, track.key());
	for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
	//
	// Declare output objects
	recob::Track outTrack;
	std::vector<art::Ptr<recob::Hit> > outHits;
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
	if (outTrack.NumberTrajectoryPoints()!=outHits.size()) {
	  throw cet::exception("TrackProducerFromPFParticle") << "Produced recob::Track required to have 1-1 correspondance between hits and points.\n";
	}
	//
	// Fill output collections, including Assns
	outputTracks->emplace_back(std::move(outTrack));
	const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size()-1);
	outputPfpTAssn->addSingle(pfp, aptr);
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
	    const art::Ptr<recob::SpacePoint> apsp = spacePointPtrMaker(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(trhit, apsp);
	  }
	  ip++;
	}
	if (doSpacePoints_ && !spacePointsFromTrajP_) {
	  auto osp = optionals.spacePointHitPairs();
	  for (auto it = osp.begin(); it!=osp.end(); ++it ) {
	    outputSpacePoints->emplace_back(std::move(it->first));
	    const art::Ptr<recob::SpacePoint> apsp = spacePointPtrMaker(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(it->second,apsp);
	  }
	}
	if (doTrackFitHitInfo_) {
	  outputHitInfo->emplace_back(std::move(optionals.trackFitHitInfos()));
	}
      }
    }
    //
    if (showerFromPF_) {
      //
      // Showers associated to PFParticles
      const std::vector<art::Ptr<recob::Shower> >& showers = assocShowers->at(iPfp);
      // if there is more than one shower the logic below to get the hits does not work! this works, at least for uboone
      if (showers.size()!=1) continue;
      //
      // Get hits for shower (through the chain pfp->clusters->hits)
      std::vector<art::Ptr<recob::Hit> > inHits;
      decltype(auto) clustersRange = util::groupByIndex(pfpClusterGroups, pfp.key());
      for (art::Ptr<recob::Cluster> const& cluster: clustersRange) {
	decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, cluster.key());
	for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
      }
      // Loop over showers to refit (should be only one)
      for (unsigned int iShower = 0; iShower < showers.size(); ++iShower) {
	//
	// Get the shower and convert/hack it into a trajectory so that the fit is initialized
	art::Ptr<recob::Shower> shower = showers[iShower++];
	recob::tracking::Point_t pos(shower->ShowerStart().X(),shower->ShowerStart().Y(),shower->ShowerStart().Z());
	recob::tracking::Vector_t dir(shower->Direction().X(),shower->Direction().Y(),shower->Direction().Z());
	std::vector<recob::tracking::Point_t> p;
	std::vector<recob::tracking::Vector_t> d;
	for (unsigned int i=0; i<inHits.size(); ++i) {
	  p.push_back(pos);
	  d.push_back(dir);
	}
	recob::TrackTrajectory traj(std::move(p), std::move(d), recob::TrackTrajectory::Flags_t(p.size()), false);
	//
	// Declare output objects
	recob::Track outTrack;
	std::vector<art::Ptr<recob::Hit> > outHits;
	trkmkr::OptionalOutputs optionals;
	if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
	if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
	//
	// Invoke tool to fit track and fill output objects
	bool fitok = trackMaker_->makeTrack(traj, iShower, inHits, outTrack, outHits, optionals);
	if (!fitok) continue;
	//
	// Check that the requirement Nhits == Npoints is satisfied
	// We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
	if (outTrack.NumberTrajectoryPoints()!=outHits.size()) {
	  throw cet::exception("TrackProducerFromPFParticle") << "Produced recob::Track required to have 1-1 correspondance between hits and points.\n";
	}
	//
	// Fill output collections, including Assns
	outputTracks->emplace_back(std::move(outTrack));
	const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size()-1);
	outputPfpTAssn->addSingle(pfp, aptr);
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
	    const art::Ptr<recob::SpacePoint> apsp = spacePointPtrMaker(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(trhit, apsp);
	  }
	  ip++;
	}
	if (doSpacePoints_ && !spacePointsFromTrajP_) {
	  auto osp = optionals.spacePointHitPairs();
	  for (auto it = osp.begin(); it!=osp.end(); ++it ) {
	    outputSpacePoints->emplace_back(std::move(it->first));
	    const art::Ptr<recob::SpacePoint> apsp = spacePointPtrMaker(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(it->second,apsp);
	  }
	}
	if (doTrackFitHitInfo_) {
	  outputHitInfo->emplace_back(std::move(optionals.trackFitHitInfos()));
	}
      }
    }
    //
  }
  //
  // Put collections in the event
  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));
  e.put(std::move(outputPfpTAssn));
  if (doTrackFitHitInfo_) {
    e.put(std::move(outputHitInfo));
  }
  if (doSpacePoints_) {
    e.put(std::move(outputSpacePoints));
    e.put(std::move(outputHitSpacePointAssn));
  }
}
//
DEFINE_ART_MODULE(TrackProducerFromPFParticle)
