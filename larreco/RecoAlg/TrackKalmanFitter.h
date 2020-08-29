#ifndef TRACKKALMANFITTER_H
#define TRACKKALMANFITTER_H

#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoObjects/KFTrackState.h"
#include "lardataobj/RecoBase/TrajectoryPointFlags.h"

namespace detinfo {
  class DetectorPropertiesData;
}

namespace recob {
  class Hit;
  class Track;
  class TrackTrajectory;
}

namespace trkmkr {
  struct OptionalOutputs;
}

namespace trkf {

  class TrackStatePropagator;

  /**
   * @file  larreco/RecoAlg/TrackKalmanFitter.h
   * @class trkf::TrackKalmanFitter
   *
   * @brief Fit tracks using Kalman Filter fit+smooth.
   *
   * This algorithm fits tracks using a Kalman Filter forward fit followed by a backward smoothing. The resulting track will feature covariance matrices at start and end positions.
   * Fundamental components of the fit are trkf::KFTrackState and trkf::TrackStatePropagator.
   *
   * Inputs are: recob::TrackTrajectory, initial covariance, associated hits, momentum estimate, particle id hypothesis. Alternatively, instead of a TrackTrajectory the fit can be input with initial position and direction, and vector of recob::TrajectoryPointFlags.
   *
   * Outputs are: resulting recob::Track, associated hits, and trkmkr::OptionalOutputs.
   *
   * For configuration options see TrackKalmanFitter#Config
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  class TrackKalmanFitter {

  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> useRMS{Name("useRMSError"),
                               Comment("Flag to replace the default hit error "
                                       "recob::Hit::SigmaPeakTime() with recob::Hit::RMS()."),
                               true};
      fhicl::Atom<bool> sortHitsByPlane{
        Name("sortHitsByPlane"),
        Comment("Flag to sort hits along the forward fit. The hit order in each plane is preserved "
                "(unless sortHitsByWire is true), the next hit to process in 3D is chosen as the "
                "one with shorter 3D propagation distance among the next hit in all planes."),
        true};
      fhicl::Atom<bool> sortHitsByWire{
        Name("sortHitsByWire"),
        Comment("Set to true if, instead of keeping the hit sorting in each plane from the pattern "
                "recognition stage, the hits need to be sorted by wire number. Ignored if "
                "sortHitsByPlane = false."),
        false};
      fhicl::Atom<bool> sortOutputHitsMinLength{
        Name("sortOutputHitsMinLength"),
        Comment("Flag to decide whether the hits are sorted before creating the output track in "
                "order to avoid tracks with huge length."),
        true};
      fhicl::Atom<bool> skipNegProp{
        Name("skipNegProp"),
        Comment(
          "Flag to decide whether, during the forward fit, the hits corresponding to a negative "
          "propagation distance should be dropped. Also, if sortOutputHitsMinLength is true, "
          "during sorting hits at a negative distance with respect to the previous are rejected."),
        true};
      fhicl::Atom<bool> cleanZigzag{
        Name("cleanZigzag"),
        Comment("Flag to decide whether hits with a zigzag pattern should be iteratively removed. "
                "Zigzag identified as negative dot product of segments connecting a point to the "
                "points before and after it."),
        false};
      fhicl::Atom<bool> rejectHighMultHits{
        Name("rejectHighMultHits"),
        Comment("Flag to rejects hits with recob::Hit::Multiplicity()>1."),
        false};
      fhicl::Atom<bool> rejectHitsNegativeGOF{
        Name("rejectHitsNegativeGOF"),
        Comment("Flag to rejects hits with recob::Hit::GoodnessOfFit<0."),
        true};
      fhicl::Atom<float> hitErr2ScaleFact{Name("hitErr2ScaleFact"),
                                          Comment("Scale the hit error squared by this factor."),
                                          1.0};
      fhicl::Atom<bool> tryNoSkipWhenFails{
        Name("tryNoSkipWhenFails"),
        Comment("In case skipNegProp is true and the track fit fails, make a second attempt to fit "
                "the track with skipNegProp=false in order to attempt to avoid losing efficiency."),
        true};
      fhicl::Atom<bool> tryBothDirs{
        Name("tryBothDirs"),
        Comment("Try fit in both with default and reversed direction, choose the track with "
                "highest score=CountValidPoints/(Length*Chi2PerNdof)."),
        false};
      fhicl::Atom<bool> pickBestHitOnWire{
        Name("pickBestHitOnWire"),
        Comment("If there is >1 consecutive hit on the same wire, choose the one with best chi2 "
                "and exclude the others."),
        false};
      fhicl::Atom<float> maxResidue{
        Name("maxResidue"),
        Comment("Reject hits with residue > maxResidue [cm]. If negative, it is set to "
                "std::numeric_limits<float>::max()."),
        -1.};
      fhicl::Atom<float> maxResidueFirstHit{
        Name("maxResidueFirstHit"),
        Comment("Reject firt hit if has residue > maxResidueFirstHit [cm]. If negative, it is set "
                "to std::numeric_limits<float>::max()."),
        -1.};
      fhicl::Atom<float> maxChi2{Name("maxChi2"),
                                 Comment("Reject hits with chi2 > maxChi2. If negative, it is set "
                                         "to std::numeric_limits<float>::max()."),
                                 -1.};
      fhicl::Atom<float> maxDist{
        Name("maxDist"),
        Comment("Reject hits with propagation distance > maxDist [cm]. If negative, it is set to "
                "std::numeric_limits<float>::max()."),
        -1.};
      fhicl::Atom<float> negDistTolerance{
        Name("negDistTolerance"),
        Comment("Tolerance for negative propagation distance to avoid hit rejection (so this is "
                "expected to be a small negative number)."),
        0.};
      fhicl::Atom<int> dumpLevel{
        Name("dumpLevel"),
        Comment("0 for no debug printouts, 1 for moderate, 2 for maximum."),
        0};
    };
    using Parameters = fhicl::Table<Config>;

    /// Constructor from TrackStatePropagator and values of configuration parameters
    TrackKalmanFitter(const TrackStatePropagator* prop,
                      bool useRMS,
                      bool sortHitsByPlane,
                      bool sortHitsByWire,
                      bool sortOutputHitsMinLength,
                      bool skipNegProp,
                      bool cleanZigzag,
                      bool rejectHighMultHits,
                      bool rejectHitsNegativeGOF,
                      float hitErr2ScaleFact,
                      bool tryNoSkipWhenFails,
                      bool tryBothDirs,
                      bool pickBestHitOnWire,
                      float maxResidue,
                      float maxResidueFirstHit,
                      float maxChi2,
                      float maxDist,
                      float negDistTolerance,
                      int dumpLevel)
    {
      propagator = prop;
      useRMS_ = useRMS;
      sortHitsByPlane_ = sortHitsByPlane;
      sortHitsByWire_ = sortHitsByWire;
      sortOutputHitsMinLength_ = sortOutputHitsMinLength;
      skipNegProp_ = skipNegProp;
      cleanZigzag_ = cleanZigzag;
      rejectHighMultHits_ = rejectHighMultHits;
      rejectHitsNegativeGOF_ = rejectHitsNegativeGOF;
      hitErr2ScaleFact_ = hitErr2ScaleFact;
      tryNoSkipWhenFails_ = tryNoSkipWhenFails;
      tryBothDirs_ = tryBothDirs;
      pickBestHitOnWire_ = pickBestHitOnWire;
      maxResidue_ = (maxResidue > 0 ? maxResidue : std::numeric_limits<float>::max());
      maxResidueFirstHit_ =
        (maxResidueFirstHit > 0 ? maxResidueFirstHit : std::numeric_limits<float>::max());
      maxChi2_ = (maxChi2 > 0 ? maxChi2 : std::numeric_limits<float>::max());
      maxDist_ = (maxDist > 0 ? maxDist : std::numeric_limits<float>::max());
      negDistTolerance_ = negDistTolerance;
      dumpLevel_ = dumpLevel;
    }

    /// Constructor from TrackStatePropagator and Parameters table
    explicit TrackKalmanFitter(const TrackStatePropagator* prop, Parameters const& p)
      : TrackKalmanFitter(prop,
                          p().useRMS(),
                          p().sortHitsByPlane(),
                          p().sortHitsByWire(),
                          p().sortOutputHitsMinLength(),
                          p().skipNegProp(),
                          p().cleanZigzag(),
                          p().rejectHighMultHits(),
                          p().rejectHitsNegativeGOF(),
                          p().hitErr2ScaleFact(),
                          p().tryNoSkipWhenFails(),
                          p().tryBothDirs(),
                          p().pickBestHitOnWire(),
                          p().maxResidue(),
                          p().maxResidueFirstHit(),
                          p().maxChi2(),
                          p().maxDist(),
                          p().negDistTolerance(),
                          p().dumpLevel())
    {}

    /// Fit track starting from TrackTrajectory
    bool fitTrack(detinfo::DetectorPropertiesData const& detProp,
                  const recob::TrackTrajectory& traj,
                  int tkID,
                  const SMatrixSym55& covVtx,
                  const SMatrixSym55& covEnd,
                  const std::vector<art::Ptr<recob::Hit>>& hits,
                  const double pval,
                  const int pdgid,
                  const bool flipDirection,
                  recob::Track& outTrack,
                  std::vector<art::Ptr<recob::Hit>>& outHits,
                  trkmkr::OptionalOutputs& optionals) const;

    /// Fit track starting from intial position, direction, and flags
    bool fitTrack(detinfo::DetectorPropertiesData const& detProp,
                  const Point_t& position,
                  const Vector_t& direction,
                  SMatrixSym55& trackStateCov,
                  const std::vector<art::Ptr<recob::Hit>>& hits,
                  const std::vector<recob::TrajectoryPointFlags>& flags,
                  const int tkID,
                  const double pval,
                  const int pdgid,
                  recob::Track& outTrack,
                  std::vector<art::Ptr<recob::Hit>>& outHits,
                  trkmkr::OptionalOutputs& optionals) const;

    /// Function where the core of the fit is performed
    bool doFitWork(KFTrackState& trackState,
                   detinfo::DetectorPropertiesData const& detProp,
                   std::vector<HitState>& hitstatev,
                   std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
                   std::vector<KFTrackState>& fwdPrdTkState,
                   std::vector<KFTrackState>& fwdUpdTkState,
                   std::vector<unsigned int>& hitstateidx,
                   std::vector<unsigned int>& rejectedhsidx,
                   std::vector<unsigned int>& sortedtksidx,
                   bool applySkipClean = true) const;

  private:
    /// Return track state from intial position, direction, and covariance
    KFTrackState setupInitialTrackState(const Point_t& position,
                                        const Vector_t& direction,
                                        SMatrixSym55& trackStateCov,
                                        const double pval,
                                        const int pdgid) const;

    /// Setup vectors of HitState and Masks to be used during the fit
    bool setupInputStates(detinfo::DetectorPropertiesData const& detProp,
                          const std::vector<art::Ptr<recob::Hit>>& hits,
                          const std::vector<recob::TrajectoryPointFlags>& flags,
                          const KFTrackState& trackState,
                          std::vector<HitState>& hitstatev,
                          std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv) const;

    /// Sort the output states
    void sortOutput(std::vector<HitState>& hitstatev,
                    std::vector<KFTrackState>& fwdUpdTkState,
                    std::vector<unsigned int>& hitstateidx,
                    std::vector<unsigned int>& rejectedhsidx,
                    std::vector<unsigned int>& sortedtksidx,
                    std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
                    bool applySkipClean = true) const;

    /// Fill the output objects
    bool fillResult(const std::vector<art::Ptr<recob::Hit>>& inHits,
                    const int tkID,
                    const int pdgid,
                    std::vector<HitState>& hitstatev,
                    std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
                    std::vector<KFTrackState>& fwdPrdTkState,
                    std::vector<KFTrackState>& fwdUpdTkState,
                    std::vector<unsigned int>& hitstateidx,
                    std::vector<unsigned int>& rejectedhsidx,
                    std::vector<unsigned int>& sortedtksidx,
                    recob::Track& outTrack,
                    std::vector<art::Ptr<recob::Hit>>& outHits,
                    trkmkr::OptionalOutputs& optionals) const;

    art::ServiceHandle<geo::Geometry const> geom;
    const TrackStatePropagator* propagator;
    bool useRMS_;
    bool sortHitsByPlane_;
    bool sortHitsByWire_;
    bool sortOutputHitsMinLength_;
    bool skipNegProp_;
    bool cleanZigzag_;
    bool rejectHighMultHits_;
    bool rejectHitsNegativeGOF_;
    float hitErr2ScaleFact_;
    bool tryNoSkipWhenFails_;
    bool tryBothDirs_;
    bool pickBestHitOnWire_;
    float maxResidue_;
    float maxResidueFirstHit_;
    float maxChi2_;
    float maxDist_;
    float negDistTolerance_;
    int dumpLevel_;
  };

}

#endif
