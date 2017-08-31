#ifndef TRACKKALMANFITTER_H
#define TRACKKALMANFITTER_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardata/RecoObjects/KFTrackState.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

namespace recob {
  class Hit;
}

class TVector3;
namespace trkmkr {
  class OptionalOutputs;
}

namespace trkf {

  /**
   * @brief Fit tracks using Kalman Filter fit+smooth.
   *
   * This algorithm fits tracks using a Kalman Filter forward fit followed by a backward smoothing.
   *
   * Inputs are: track object, associated hits, and momentum estimate.
   * Output are: resulting track and associated hits. The resulting track will feature covariance matrices at start and end positions.
   *
   */

  class TrackKalmanFitter {

  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> useRMS {
	Name("useRMSError"),
	Comment("Flag to replace the default hit error SigmaPeakTime() with RMS().")
      };
      fhicl::Atom<bool> sortHitsByPlane {
        Name("sortHitsByPlane"),
        Comment("Flag to sort hits along the forward fit. The hit order in each plane is preserved, the next hit to process in 3D is chosen as the one with shorter 3D propagation distance among the next hit in all planes.")
      };
      fhicl::Atom<bool> sortOutputHitsMinLength {
        Name("sortOutputHitsMinLength"),
        Comment("Flag to decide whether the hits are sorted before creating the output track in order to avoid tracks with huge length.")
      };
      fhicl::Atom<bool> skipNegProp {
        Name("skipNegProp"),
        Comment("Flag to decide whether, during the forward fit, the hits corresponding to a negative propagation distance should be dropped. Also, if sortOutputHitsMinLength is true, during sorting it rejects hits at a negative distance with respect to the previous.")
      };
      fhicl::Atom<bool> cleanZigzag {
        Name("cleanZigzag"),
        Comment("Flag to decide whether hits with a zigzag pattern should be iteratively removed. Zigzag identified as negative dot product of segments connecting a point to the points before and after it.")
      };
      fhicl::Atom<bool> rejectHighMultHits {
        Name("rejectHighMultHits"),
        Comment("Flag to rejects hits with Multiplicity()>1.")
      };
      fhicl::Atom<bool> rejectHitsNegativeGOF {
        Name("rejectHitsNegativeGOF"),
        Comment("Flag to rejects hits with GoodnessOfFit<0.")
      };
      fhicl::Atom<float> hitErr2ScaleFact {
        Name("hitErr2ScaleFact"),
	Comment("Scale the hit error squared by this factor.")
      };
      fhicl::Atom<bool> tryNoSkipWhenFails {
        Name("tryNoSkipWhenFails"),
        Comment("In case skipNegProp is true and the track fit fails, make a second attempt to fit the track with skipNegProp=false in order to attempt to avoid losing efficiency.")
      };
      fhicl::Atom<int> dumpLevel {
        Name("dumpLevel"),
	Comment("0 for no debug printouts, 1 for moderate, 2 for maximum.")
      };
    };
    using Parameters = fhicl::Table<Config>;

    TrackKalmanFitter(const TrackStatePropagator* prop, bool useRMS, bool sortHitsByPlane, bool sortOutputHitsMinLength, bool skipNegProp, bool cleanZigzag,
		      bool rejectHighMultHits, bool rejectHitsNegativeGOF, float hitErr2ScaleFact, bool tryNoSkipWhenFails, int dumpLevel){
      propagator=prop;
      useRMS_=useRMS;
      sortHitsByPlane_=sortHitsByPlane;
      sortOutputHitsMinLength_=sortOutputHitsMinLength;
      skipNegProp_=skipNegProp;
      cleanZigzag_=cleanZigzag;
      rejectHighMultHits_=rejectHighMultHits;
      rejectHitsNegativeGOF_=rejectHitsNegativeGOF;
      hitErr2ScaleFact_=hitErr2ScaleFact;
      tryNoSkipWhenFails_=tryNoSkipWhenFails;
      dumpLevel_=dumpLevel;
      detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    }
    explicit TrackKalmanFitter(const TrackStatePropagator* prop, Parameters const & p)
      : TrackKalmanFitter(prop,p().useRMS(),p().sortHitsByPlane(),p().sortOutputHitsMinLength(),p().skipNegProp(),p().cleanZigzag(),
			  p().rejectHighMultHits(),p().rejectHitsNegativeGOF(),p().hitErr2ScaleFact(),p().tryNoSkipWhenFails(),p().dumpLevel()) {}

    bool fitTrack(const recob::TrackTrajectory& traj, int tkID,  const SMatrixSym55& covVtx, const SMatrixSym55& covEnd,
		  const std::vector<art::Ptr<recob::Hit> >& hits, const double pval, const int pdgid, const bool flipDirection,
		  recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, trkmkr::OptionalOutputs& optionals) const;

    bool fitTrack(const Point_t& position, const Vector_t& direction, SMatrixSym55& trackStateCov,
		  const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags,
		  const int tkID, const double pval, const int pdgid,
		  recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, trkmkr::OptionalOutputs& optionals) const;

    KFTrackState setupInitialTrackState(const Point_t& position, const Vector_t& direction, SMatrixSym55& trackStateCov, const double pval, const int pdgid) const;

    bool setupInputStates(const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags,
			  const KFTrackState& trackState, bool& reverseHits,
			  std::vector<HitState>& hitstatev, std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv) const;

    bool doFitWork(KFTrackState& trackState, std::vector<HitState>& hitstatev, std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
		   std::vector<KFTrackState>& fwdPrdTkState, std::vector<KFTrackState>& fwdUpdTkState,
		   std::vector<unsigned int>& hitstateidx, std::vector<unsigned int>& rejectedhsidx, std::vector<unsigned int>& sortedtksidx, 
		   bool applySkipClean = true) const;

    void sortOutput(std::vector<HitState>& hitstatev, std::vector<KFTrackState>& fwdUpdTkState,
		    std::vector<unsigned int>& hitstateidx, std::vector<unsigned int>& rejectedhsidx,
		    std::vector<unsigned int>& sortedtksidx, bool applySkipClean = true) const;

    bool fillResult(const std::vector<art::Ptr<recob::Hit> >& inHits, const int tkID, const int pdgid, const bool reverseHits,
		    std::vector<HitState>& hitstatev, std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
		    std::vector<KFTrackState>& fwdPrdTkState, std::vector<KFTrackState>& fwdUpdTkState,
		    std::vector<unsigned int>& hitstateidx, std::vector<unsigned int>& rejectedhsidx, std::vector<unsigned int>& sortedtksidx,
		    recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, trkmkr::OptionalOutputs& optionals) const;

  private:
    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::DetectorProperties* detprop;
    const TrackStatePropagator* propagator;
    bool useRMS_;
    bool sortHitsByPlane_;
    bool sortOutputHitsMinLength_;
    bool skipNegProp_;
    bool cleanZigzag_;
    bool rejectHighMultHits_;
    bool rejectHitsNegativeGOF_;
    float hitErr2ScaleFact_;
    bool tryNoSkipWhenFails_;
    int dumpLevel_;
  };

}

#endif
