#ifndef TRACKMAKER_H
#define TRACKMAKER_H

////////////////////////////////////////////////////////////////////////
// Class:       TrackMaker
// File:        TrackMaker.h
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "canvas/Persistency/Common/Assns.h"

namespace trkmkr {

  /**
   * @brief Struct holding point-by-point elements used in OptionalOutputs.
   */

  struct OptionalPointElement {
  public:
    void setTrackFitHitInfo(recob::TrackFitHitInfo&& aTrackFitHitInfo) { trackFitHitInfo = std::make_unique< recob::TrackFitHitInfo >(aTrackFitHitInfo); }
    bool isTrackFitInfoSet() { return bool(trackFitHitInfo); }
    recob::TrackFitHitInfo&& moveTrackFitHitInfo() { return std::move( *(trackFitHitInfo.release()) ); }
    //
    void setSpacePoint(recob::SpacePoint&& aSpacePoint) { spacePoint = std::make_unique< recob::SpacePoint >(aSpacePoint); }
    bool isSpacePointSet() { return bool(spacePoint); }
    recob::SpacePoint&& moveSpacePoint() { return std::move( *(spacePoint.release()) ); }
  private:
    std::unique_ptr< recob::TrackFitHitInfo > trackFitHitInfo;
    std::unique_ptr<recob::SpacePoint> spacePoint;
  };

  /**
   * @brief Struct holding optional TrackMaker outputs.
   */

  struct OptionalOutputs {
  public:
    typedef std::pair<recob::SpacePoint, art::Ptr<recob::Hit> > SpHitPair;
    void addPoint(OptionalPointElement& ope) {
      if (isTrackFitInfosInit() && ope.isTrackFitInfoSet()) {
	outTrackFitHitInfos->push_back( std::move(ope.moveTrackFitHitInfo()) );
      }
    }
    void addPoint(OptionalPointElement& ope, art::Ptr<recob::Hit> hptr) {
      if (isSpacePointsInit() && ope.isSpacePointSet()) {
	outSpacePointHitPairs->push_back( std::move( SpHitPair(ope.moveSpacePoint(), hptr) ) );
      }
      addPoint(ope);
    }
    void reset() {
      if (isTrackFitInfosInit()) {
	outTrackFitHitInfos.reset();
	initTrackFitInfos();
      }
      if (isSpacePointsInit()) {
	outSpacePointHitPairs.reset();
	initSpacePoints();
      }
    }
    void initTrackFitInfos() {
      outTrackFitHitInfos = std::make_unique< std::vector<recob::TrackFitHitInfo> >();
    }
    void initSpacePoints() {
      outSpacePointHitPairs = std::make_unique< std::vector<SpHitPair> >();
    }
    bool isTrackFitInfosInit() { return bool(outTrackFitHitInfos); }
    bool isSpacePointsInit() { return bool(outSpacePointHitPairs); }
    std::vector<recob::TrackFitHitInfo> trackFitHitInfos() {
      if (!isTrackFitInfosInit()) throw std::logic_error("outTrackFitHitInfos is not available (any more?).");
      return std::move(*(outTrackFitHitInfos.release() ));
    }
    std::vector<SpHitPair> spacePointHitPairs() {
      if (!isSpacePointsInit()) throw std::logic_error("outSpacePointHitPairs is not available (any more?).");
      return std::move(*(outSpacePointHitPairs.release() ));
    }
  private:
    std::unique_ptr< std::vector<recob::TrackFitHitInfo> > outTrackFitHitInfos;
    std::unique_ptr< std::vector<SpHitPair> > outSpacePointHitPairs;
  };

  /**
   * @brief Base abstract class for tools used to fit tracks.
   *
   * Base abstract class for tools used to fit tracks. The virtual function makeTrack comes in three flavours,
   * one for each possible input (Trajectory, TrackTrajectory, Track).
   * The only purely virtual function is the one for input Trajectories (by default the other two just forward the call to it).
   * Its arguments are the const inputs (Trajectory, TrajectoryPointFlags, track ID)
   * and the non-const ouputs (mandatory: outTrack and outHits; optional outputs stored in OptionalOutputs).
   * In case other products are needed from the event (e.g. associations to the input), they can be retrieved overloading the initEvent function.
   * The tool is not meant to put collections in the event.
   * Requirements are that a Track has at least 2 points, that it has the same number of Points and Momenta,
   * that TrajectoryPoints and Hit have a 1-1 correspondance (same number and  same order).
   * The functions return a bool corresponding to the success or failure status of the fit.
   */

  class TrackMaker {
  public:
    virtual ~TrackMaker() noexcept = default;

    virtual void initEvent(const art::Event& e) { return; }

    virtual bool makeTrack(const recob::Trajectory& traj, const std::vector<recob::TrajectoryPointFlags>& flags,
			   const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
			   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const
    {
      return makeTrack(recob::TrackTrajectory(traj,std::move(recob::TrackTrajectory::Flags_t(traj.NPoints()))), tkID, inHits, outTrack, outHits, optionals);
    }

    virtual bool makeTrack(const art::Ptr<recob::Trajectory> traj, const std::vector<recob::TrajectoryPointFlags>& flags,
			   const std::vector<art::Ptr<recob::Hit> >& inHits,
			   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const
    {
      return makeTrack(recob::TrackTrajectory(*traj,std::move(recob::TrackTrajectory::Flags_t(traj->NPoints()))), traj.key(), inHits, outTrack, outHits, optionals);
    }

    virtual bool makeTrack(const art::Ptr<recob::TrackTrajectory> ttraj, const std::vector<art::Ptr<recob::Hit> >& inHits,
			   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const
    {
      return makeTrack(*ttraj, ttraj.key(), inHits, outTrack, outHits, optionals);
    }

    virtual bool makeTrack(const recob::TrackTrajectory& ttraj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
			   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const = 0;

    virtual bool makeTrack(const art::Ptr<recob::Track> track, const std::vector<art::Ptr<recob::Hit> >& inHits,
			   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const
    {
      return makeTrack(track->Trajectory(), track.key(), inHits, outTrack, outHits, optionals);
    }

    virtual bool makeTrack(const recob::Track& track, const std::vector<art::Ptr<recob::Hit> >& inHits,
			   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const
    {
      return makeTrack(track.Trajectory(), track.ID(), inHits, outTrack, outHits, optionals);
    }
  };
}

#endif
