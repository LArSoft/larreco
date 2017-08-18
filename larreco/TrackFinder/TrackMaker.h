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
 * @brief Base abstract class for tools used to fit tracks.
 *
 * Base abstract class for tools used to fit tracks. The virtual function makeTrack comes in three flavours,
 * one for each possible input (Trajectory, TrackTrajectory, Track).
 * The only purely virtual function is the one for input Trajectories (by default the other two just forward the call to it).
 * Its arguments are the const inputs (Trajectory, TrajectoryPointFlags, track ID)
 * and the non-const ouputs (mandatory: outTrack and outHits; optional: TrackFitHitInfos, SpacePoints, SpacePoints-Hits Assns);
 * a const reference to the Event is also provided in case the tool needs to access other data products.
 * The tool is not meant to put collections in the event.
 * Requirements are that a Track has at least 2 points, that it has the same number of Points and Momenta,
 * that TrajectoryPoints and Hit have a 1-1 correspondance (same number and  same order).
 * The functions return a bool corresponding to the success or failure status of the fit.
 */

  struct OptionalOutputs {
  public:
    OptionalOutputs() { outTrackFitHitInfos = nullptr; outSpacePoints = nullptr; outHitSpacePointAssn = nullptr; }
    ~OptionalOutputs() {
      if (outTrackFitHitInfos) delete outTrackFitHitInfos;
      if (outSpacePoints) delete outSpacePoints;
      if (outHitSpacePointAssn) delete outHitSpacePointAssn;
    }
    void initTrackFitInfos() { outTrackFitHitInfos = new std::vector<recob::TrackFitHitInfo>(); }
    void initSpacePoints() {
      outSpacePoints = new std::vector<recob::SpacePoint>();
      outHitSpacePointAssn = new art::Assns<recob::Hit, recob::SpacePoint>();
    }
    bool isTrackFitInfosInit() { return outTrackFitHitInfos != nullptr; }
    bool isSpacePointssInit() { return outSpacePoints != nullptr; }
    std::vector<recob::TrackFitHitInfo>* trackFitHitInfos() { return outTrackFitHitInfos; }
    std::vector<recob::SpacePoint>* spacePoints() { return outSpacePoints; }
    art::Assns<recob::Hit, recob::SpacePoint>* hitSpacePointAssn() { return outHitSpacePointAssn; }
  private:
    std::vector<recob::TrackFitHitInfo>* outTrackFitHitInfos;
    std::vector<recob::SpacePoint>* outSpacePoints;
    art::Assns<recob::Hit, recob::SpacePoint>* outHitSpacePointAssn;
  };

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
      return makeTrack(ttraj->Trajectory(), ttraj->Flags(), ttraj.key(), inHits, outTrack, outHits, optionals);
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
