#ifndef TRACKMAKER_H
#define TRACKMAKER_H

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
namespace detinfo {
  class DetectorPropertiesData;
}

namespace trkmkr {

  /**
   * @file  larreco/TrackFinder/TrackMaker.h
   * @struct trkmkr::OptionalPointElement
   *
   * @brief Struct holding point-by-point elements used in OptionalOutputs.
   *
   * This struct holds the elements of OptionalOutputs that are added for each point (i.e. each hit).
   *
   * It stores a unique_ptr to each optional output object element.
   * Functions are provided to set the unique_ptr and to check if it set.
   * When the elements are returned, the unique_ptr is reset.
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  struct OptionalPointElement {
  public:
    /// set the recob::TrackFitHitInfo unique_ptr
    void setTrackFitHitInfo(recob::TrackFitHitInfo&& aTrackFitHitInfo)
    {
      trackFitHitInfo = std::make_unique<recob::TrackFitHitInfo>(std::move(aTrackFitHitInfo));
    }
    void setTrackFitHitInfo(const recob::TrackFitHitInfo& aTrackFitHitInfo)
    {
      trackFitHitInfo = std::make_unique<recob::TrackFitHitInfo>(aTrackFitHitInfo);
    }
    /// check if the recob::TrackFitHitInfo unique_ptr is set
    bool isTrackFitInfoSet() { return bool(trackFitHitInfo); }
    /// get the recob::TrackFitHitInfo object, and reset the unique_ptr
    recob::TrackFitHitInfo getTrackFitHitInfo()
    {
      auto tmp = *trackFitHitInfo;
      trackFitHitInfo.reset();
      return tmp;
    }
    //
    /// set the recob::SpacePoint unique_ptr
    void setSpacePoint(recob::SpacePoint&& aSpacePoint)
    {
      spacePoint = std::make_unique<recob::SpacePoint>(aSpacePoint);
    }
    void setSpacePoint(const recob::SpacePoint& aSpacePoint)
    {
      spacePoint = std::make_unique<recob::SpacePoint>(aSpacePoint);
    }
    /// check if the recob::SpacePoint unique_ptr is set
    bool isSpacePointSet() { return bool(spacePoint); }
    /// get the recob::SpacePoint object, and release the unique_ptr
    recob::SpacePoint getSpacePoint()
    {
      auto tmp = *spacePoint;
      spacePoint.reset();
      return tmp;
    }

  private:
    std::unique_ptr<recob::TrackFitHitInfo> trackFitHitInfo;
    std::unique_ptr<recob::SpacePoint> spacePoint;
  };

  /**
   * @file  larreco/TrackFinder/TrackMaker.h
   * @struct trkmkr::OptionalOutputs
   *
   * @brief Struct holding optional TrackMaker outputs.
   *
   * This struct holds the optional outputs of track making and hides their details to the actual track making tools.
   * In this way, adding a new optional output will affect only those tools that produce such new ouput.
   *
   * It stores a unique_ptr to the vector of each optional output object (meant to be per-track).
   * Track making tools need to init the outional outputs they will produce, so that only the unique_ptrs that are needed are actually created.
   * Functions are provided (called addPoint) to add point-by-point elements (see OptionalPointElement).
   * When the output objects are returned, the unique_ptr is reset, so that no new elements should be added and a new initialization is needed.
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  struct OptionalOutputs {
  public:
    typedef std::pair<recob::SpacePoint, art::Ptr<recob::Hit>> SpHitPair;

    /// add one OptionalPointElement
    void addPoint(OptionalPointElement& ope)
    {
      if (isTrackFitInfosInit() && ope.isTrackFitInfoSet()) {
        outTrackFitHitInfos->push_back(ope.getTrackFitHitInfo());
      }
    }
    /// add one OptionalPointElement and the corresponding hit
    void addPoint(OptionalPointElement& ope, art::Ptr<recob::Hit> hptr)
    {
      if (isSpacePointsInit() && ope.isSpacePointSet()) {
        outSpacePointHitPairs->emplace_back(ope.getSpacePoint(), hptr);
      }
      addPoint(ope);
    }
    /// reset the stored vectors
    void reset()
    {
      if (isTrackFitInfosInit()) {
        outTrackFitHitInfos.reset();
        initTrackFitInfos();
      }
      if (isSpacePointsInit()) {
        outSpacePointHitPairs.reset();
        initSpacePoints();
      }
    }
    /// initialize the output vector of TrackFitHitInfos
    void initTrackFitInfos()
    {
      outTrackFitHitInfos = std::make_unique<std::vector<recob::TrackFitHitInfo>>();
    }
    /// initialize the output vector of SpHitPair
    void initSpacePoints() { outSpacePointHitPairs = std::make_unique<std::vector<SpHitPair>>(); }
    /// check initialization of the output vector of TrackFitHitInfos
    bool isTrackFitInfosInit() { return bool(outTrackFitHitInfos); }
    /// check initialization of the output vector of SpHitPair
    bool isSpacePointsInit() { return bool(outSpacePointHitPairs); }
    /// get the output vector of TrackFitHitInfos by releasing and moving
    std::vector<recob::TrackFitHitInfo> trackFitHitInfos()
    {
      if (!isTrackFitInfosInit())
        throw std::logic_error("outTrackFitHitInfos is not available (any more?).");
      auto tmp = *outTrackFitHitInfos;
      outTrackFitHitInfos.reset();
      return tmp;
    }
    /// get the output vector of SpHitPair by releasing and moving
    std::vector<SpHitPair> spacePointHitPairs()
    {
      if (!isSpacePointsInit())
        throw std::logic_error("outSpacePointHitPairs is not available (any more?).");
      auto tmp = *outSpacePointHitPairs;
      outSpacePointHitPairs.reset();
      return tmp;
    }

  private:
    std::unique_ptr<std::vector<recob::TrackFitHitInfo>> outTrackFitHitInfos;
    std::unique_ptr<std::vector<SpHitPair>> outSpacePointHitPairs;
  };

  /**
   * @file  larreco/TrackFinder/TrackMaker.h
   * @class trkmkr::TrackMaker
   *
   * @brief Base abstract class for tools used to fit tracks.
   *
   * The virtual function makeTrack comes in different versions, one for each possible input (Trajectory, TrackTrajectory, Track), both using const references and art pointers as input).
   * The functions return a bool corresponding to the success or failure status of the fit.
   *
   * The only purely virtual function is the one for input Trajectories (by default the other two just forward the call to it).
   * Its arguments are the const inputs (Trajectory, TrajectoryPointFlags, track ID)
   * and the non-const ouputs (mandatory: outTrack and outHits; optional outputs stored in OptionalOutputs).
   *
   * In case other products are needed from the event (e.g. associations to the input), they can be retrieved overriding the initEvent function.
   *
   * The tool is not meant to put collections in the event.
   *
   * Requirements are that a Track has at least 2 points, that it has the same number of Points and Momenta,
   * that TrajectoryPoints and Hit have a 1-1 correspondance (same number and  same order).
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  class TrackMaker {
  public:
    virtual ~TrackMaker() noexcept = default;

    /// per-event initialization; concrete classes may override this function to retrieve other products or associations from the event.
    virtual void initEvent(const art::Event& e) {}

    //@{
    /// makeTrack functions with recob::Trajectory as argument; calls the version with recob::TrackTrajectory using a dummy flags vector.
    virtual bool makeTrack(const detinfo::DetectorPropertiesData& detProp,
                           const recob::Trajectory& traj,
                           const std::vector<recob::TrajectoryPointFlags>& flags,
                           const int tkID,
                           const std::vector<art::Ptr<recob::Hit>>& inHits,
                           recob::Track& outTrack,
                           std::vector<art::Ptr<recob::Hit>>& outHits,
                           OptionalOutputs& optionals) const
    {
      return makeTrack(
        detProp,
        recob::TrackTrajectory(traj, recob::TrackTrajectory::Flags_t(traj.NPoints())),
        tkID,
        inHits,
        outTrack,
        outHits,
        optionals);
    }
    virtual bool makeTrack(const detinfo::DetectorPropertiesData& detProp,
                           const art::Ptr<recob::Trajectory> traj,
                           const std::vector<recob::TrajectoryPointFlags>& flags,
                           const std::vector<art::Ptr<recob::Hit>>& inHits,
                           recob::Track& outTrack,
                           std::vector<art::Ptr<recob::Hit>>& outHits,
                           OptionalOutputs& optionals) const
    {
      return makeTrack(
        detProp,
        recob::TrackTrajectory(*traj, recob::TrackTrajectory::Flags_t(traj->NPoints())),
        traj.key(),
        inHits,
        outTrack,
        outHits,
        optionals);
    }
    //@}

    /// makeTrack functions with art::Ptr<recob::TrackTrajectory>; calls the purely virtual version with const recob::TrackTrajectory reference as argument.
    virtual bool makeTrack(const detinfo::DetectorPropertiesData& detProp,
                           const art::Ptr<recob::TrackTrajectory> ttraj,
                           const std::vector<art::Ptr<recob::Hit>>& inHits,
                           recob::Track& outTrack,
                           std::vector<art::Ptr<recob::Hit>>& outHits,
                           OptionalOutputs& optionals) const
    {
      return makeTrack(detProp, *ttraj, ttraj.key(), inHits, outTrack, outHits, optionals);
    }

    /// makeTrack functions with const recob::TrackTrajectory reference as
    /// argument: purely virtual, to be implemented in concrete classes.
    virtual bool makeTrack(const detinfo::DetectorPropertiesData& detProp,
                           const recob::TrackTrajectory& ttraj,
                           const int tkID,
                           const std::vector<art::Ptr<recob::Hit>>& inHits,
                           recob::Track& outTrack,
                           std::vector<art::Ptr<recob::Hit>>& outHits,
                           OptionalOutputs& optionals) const = 0;

    //@{
    /// makeTrack functions with recob::Track as argument; calls the version
    /// with recob::TrackTrajectory.
    virtual bool makeTrack(const detinfo::DetectorPropertiesData& detProp,
                           const art::Ptr<recob::Track> track,
                           const std::vector<art::Ptr<recob::Hit>>& inHits,
                           recob::Track& outTrack,
                           std::vector<art::Ptr<recob::Hit>>& outHits,
                           OptionalOutputs& optionals) const
    {
      return makeTrack(
        detProp, track->Trajectory(), track.key(), inHits, outTrack, outHits, optionals);
    }

    virtual bool makeTrack(const detinfo::DetectorPropertiesData& detProp,
                           const recob::Track& track,
                           const std::vector<art::Ptr<recob::Hit>>& inHits,
                           recob::Track& outTrack,
                           std::vector<art::Ptr<recob::Hit>>& outHits,
                           OptionalOutputs& optionals) const
    {
      return makeTrack(
        detProp, track.Trajectory(), track.ID(), inHits, outTrack, outHits, optionals);
    }
    //@}
  };
}

#endif
