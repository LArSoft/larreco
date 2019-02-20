#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/RecoAlg/TrackCreationBookKeeper.h"
#include "larreco/RecoAlg/TrackKalmanFitter.h"
#include "larreco/TrackFinder/TrackMaker.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace trkmkr {

  /**
   * @file larreco/TrackFinder/KalmanFilterFitTrackMaker_tool.cc
   * @class trkmkr::KalmanFilterFitTrackMaker
   *
   * @brief Concrete implementation of a tool to fit tracks with
   * TrackKalmanFitter.
   *
   * Concrete implementation of a tool to fit tracks with
   * trkf::TrackKalmanFitter; inherits from abstract class TrackMaker. It
   * prepares the input needed by the fitter (momentum, particleId, direction),
   * and returns a track with all outputs filled. If the flag
   * keepInputTrajectoryPoints is set to true, the tracjetory points from the
   * input track are copied into the output, so that only the covariance
   * matrices, the chi2 and the ndof in the output track are resulting from the
   * fit.
   *
   * For configuration options see KalmanFilterFitTrackMaker#Options and
   * KalmanFilterFitTrackMaker#Config.
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  class KalmanFilterFitTrackMaker : public TrackMaker {

  public:
    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      //
      fhicl::Atom<double> defaultMomInGeV{
        Name("defaultMomInGeV"),
        Comment("Default momentum estimate value (all other options are set to "
                "false, or if the estimate is not available)."),
        1.0};
      fhicl::Atom<bool> momFromMCSCollection{
        Name("momFromMCSCollection"),
        Comment("Flag used to get initial momentum estimate from MCSFitResult "
                "collection specified by mcsInputTag."),
        false};
      fhicl::Atom<art::InputTag> mcsInputTag{
        Name("mcsInputTag"),
        Comment("InputTag of MCSFitResult collection.")};
      fhicl::Atom<bool> momFromCombAndPid{
        Name("momFromCombAndPid"),
        Comment("Flag used to get initial momentum estimate from either range "
                "or mcs fit, based on particle id and containement (from "
                "contInputTag collection)."),
        false};
      fhicl::Atom<art::InputTag> contInputTag{
        Name("contInputTag"),
        Comment("InputTag of CosmicTag collection for containement.")};
      fhicl::Atom<bool> pidFromCollection{
        Name("pidFromCollection"),
        Comment("Flag used to get initial particle id estimate from ParticleID "
                "collection specified by pidInputTag."),
        false};
      fhicl::Atom<art::InputTag> pidInputTag{
        Name("pidInputTag"),
        Comment("InputTag of ParticleID collection.")};
      fhicl::Atom<double> pidFromLengthCut{
        Name("pidFromLengthCut"),
        Comment("Particle ID based on length: if shorted than cut is assumed "
                "to be a proton, if longer a muon; disabled if negative."),
        -1.};
      fhicl::Atom<int> defaultPdgId{
        Name("defaultPdgId"),
        Comment("Default particle id hypothesis (all other options are set to "
                "false, or if the estimate is not available)."),
        13};
      fhicl::Atom<bool> dirFromVec{
        Name("dirFromVec"),
        Comment("Assume track direction as the one giving positive dot product "
                "with vector specified by dirVec."),
        false};
      fhicl::Sequence<float, 3u> dirVec{
        Name("dirVec"),
        Comment("Fhicl sequence defining the vector used when dirFromVec=true. "
                "It must have 3 elements.")};
      fhicl::Atom<bool> alwaysInvertDir{
        Name("alwaysInvertDir"),
        Comment("If true, fit all tracks from end to vertex assuming inverted "
                "direction."),
        false};
      fhicl::Atom<bool> keepInputTrajectoryPoints{
        Name("keepInputTrajectoryPoints"),
        Comment(
          "Option to keep positions and directions from input "
          "trajectory/track. The fit will provide only covariance matrices, "
          "chi2, ndof, particle Id and absolute momentum. It may also modify "
          "the trajectory point flags. In order to avoid inconsistencies, it "
          "has to be used with the following fitter options all set to false: "
          "sortHitsByPlane, sortOutputHitsMinLength, skipNegProp."),
        false};
      //
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<KalmanFilterFitTrackMaker::Options> options{Name("options")};
      fhicl::Table<trkf::TrackStatePropagator::Config> propagator{
        Name("propagator")};
      fhicl::Table<trkf::TrackKalmanFitter::Config> fitter{Name("fitter")};
      fhicl::Table<trkf::TrajectoryMCSFitter::Config> mcsfit{Name("mcsfit")};
    };
    using Parameters = art::ToolConfigTable<Config>;

    /// Constructor from Parameters
    explicit KalmanFilterFitTrackMaker(Parameters const& p) :
      p_(p),
      prop{p_().propagator},
      kalmanFitter{&prop, p_().fitter},
      mcsfitter{p_().mcsfit},
      mom_def_{p_().options().defaultMomInGeV()},
      momFromMCSColl_{p_().options().momFromMCSCollection()},
      momFromCombAndPid_{p_().options().momFromCombAndPid()},
      pidFromColl_{p_().options().pidFromCollection()},
      mom_len_cut_{p_().options().pidFromLengthCut()},
      pid_def_{p_().options().defaultPdgId()},
      alwaysFlip_{p_().options().alwaysInvertDir()},
      dirFromVec_{p_().options().dirFromVec()}
    {
      if (momFromMCSColl_) {
        mcsInputTag_ = p_().options().mcsInputTag();
      }
      if (momFromCombAndPid_) {
        contInputTag_ = p_().options().contInputTag();
      }
      if (pidFromColl_) {
        pidInputTag_ = p_().options().pidInputTag();
      }
      if (dirFromVec_) {
        auto d = p_().options().dirVec();
        dirVec = recob::tracking::Vector_t{d[0], d[1], d[2]};
      }
      //
      if (p_().options().keepInputTrajectoryPoints() &&
          (p_().fitter().sortHitsByPlane() ||
           p_().fitter().sortOutputHitsMinLength() ||
           p_().fitter().skipNegProp())) {
        throw cet::exception("KalmanFilterFitTrackMaker")
          << "Incompatible configuration parameters: keepInputTrajectoryPoints "
             "needs the following fitter options all set to false: "
             "sortHitsByPlane, sortOutputHitsMinLength, skipNegProp."
          << "\n";
      }
      if (momFromMCSColl_ && momFromCombAndPid_) {
        throw cet::exception("KalmanFilterFitTrackMaker")
          << "Incompatible configuration parameters: momFromMCSCollection and "
             "momFromCombAndPid cannot be both true at the same time."
          << "\n";
      }
      if (pidFromColl_ && mom_len_cut_ > 0) {
        throw cet::exception("KalmanFilterFitTrackMaker")
          << "Incompatible configuration parameters: pidFromCollection and "
             "pidFromLengthCut cannot be respectively true and >0. at the same "
             "time."
          << "\n";
      }
      if (alwaysFlip_ && dirFromVec_) {
        throw cet::exception("KalmanFilterFitTrackMaker")
          << "Incompatible configuration parameters: alwaysInvertDir and "
             "dirFromVec cannot be both true at the same time."
          << "\n";
      }
      //
    }

    /// initialize event: get collection of recob::MCSFitResult
    virtual void
    initEvent(const art::Event& e) override
    {
      if (momFromMCSColl_)
        mcs = e.getValidHandle<std::vector<recob::MCSFitResult>>(mcsInputTag_).product();
      if (momFromCombAndPid_) {
        cont = e.getValidHandle<std::vector<anab::CosmicTag>>(contInputTag_).product();
      }
      if (pidFromColl_) {
        pid = e.getValidHandle<std::vector<anab::ParticleID>>(pidInputTag_).product();
      }
      return;
    }

    /// function that actually calls the fitter
    bool makeTrackImpl(const recob::TrackTrajectory& traj,
                       const int tkID,
                       const std::vector<art::Ptr<recob::Hit>>& inHits,
                       const SMatrixSym55& covVtx,
                       const SMatrixSym55& covEnd,
                       recob::Track& outTrack,
                       std::vector<art::Ptr<recob::Hit>>& outHits,
                       OptionalOutputs& optionals) const;

    /// override of TrackMaker purely virtual function with
    /// recob::TrackTrajectory as argument
    bool
    makeTrack(const recob::TrackTrajectory& traj,
              const int tkID,
              const std::vector<art::Ptr<recob::Hit>>& inHits,
              recob::Track& outTrack,
              std::vector<art::Ptr<recob::Hit>>& outHits,
              OptionalOutputs& optionals) const override
    {
      return makeTrackImpl(traj,
                           tkID,
                           inHits,
                           trkf::SMatrixSym55{},
                           trkf::SMatrixSym55{},
                           outTrack,
                           outHits,
                           optionals);
    }

    /// override of TrackMaker virtual function with recob::Track as argument
    bool
    makeTrack(const recob::Track& track,
              const std::vector<art::Ptr<recob::Hit>>& inHits,
              recob::Track& outTrack,
              std::vector<art::Ptr<recob::Hit>>& outHits,
              OptionalOutputs& optionals) const override
    {
      auto covs = track.Covariances();
      return makeTrackImpl(track.Trajectory(),
                           track.ID(),
                           inHits,
                           covs.first,
                           covs.second,
                           outTrack,
                           outHits,
                           optionals);
    }

    /// set the particle ID hypothesis
    int getParticleID(const recob::TrackTrajectory& traj, const int tkID) const;
    /// set the initial momentum estimate
    double getMomentum(const recob::TrackTrajectory& traj,
                       const int pid,
                       const bool isFlip,
                       const int tkID) const;
    /// decide whether to flip the direction or not
    bool isFlipDirection(const recob::TrackTrajectory& traj,
                         const int tkID) const;

    /// restore the TrajectoryPoints in the Track to be the same as those in the
    /// input TrackTrajectory (but keep covariance matrices and chi2 from fit).
    void restoreInputPoints(const recob::TrackTrajectory& traj,
                            const std::vector<art::Ptr<recob::Hit>>& inHits,
                            recob::Track& outTrack,
                            std::vector<art::Ptr<recob::Hit>>& outHits,
                            OptionalOutputs& optionals) const;

  private:
    Parameters p_;
    const trkf::TrackStatePropagator prop;
    const trkf::TrackKalmanFitter kalmanFitter;
    const trkf::TrajectoryMCSFitter mcsfitter;
    double mom_def_;
    bool momFromMCSColl_;
    art::InputTag mcsInputTag_;
    bool momFromCombAndPid_;
    art::InputTag contInputTag_;
    bool pidFromColl_;
    art::InputTag pidInputTag_;
    double mom_len_cut_;
    int pid_def_;
    bool alwaysFlip_;
    bool dirFromVec_;
    recob::tracking::Vector_t dirVec;
    const std::vector<recob::MCSFitResult>* mcs = nullptr;
    const std::vector<anab::CosmicTag>* cont = nullptr;
    const std::vector<anab::ParticleID>* pid = nullptr;
    trkf::TrackMomentumCalculator tmc;
  };
}

bool
trkmkr::KalmanFilterFitTrackMaker::makeTrackImpl(const recob::TrackTrajectory& traj,
                                                 const int tkID,
                                                 const std::vector<art::Ptr<recob::Hit>>& inHits,
                                                 const SMatrixSym55& covVtx,
                                                 const SMatrixSym55& covEnd,
                                                 recob::Track& outTrack,
                                                 std::vector<art::Ptr<recob::Hit>>& outHits,
                                                 OptionalOutputs& optionals) const
{
  //
  const int pid = getParticleID(traj, tkID);
  const bool flipDirection = isFlipDirection(traj, tkID);
  const double mom = getMomentum(traj, pid, flipDirection, tkID); // what about mom uncertainty?
  // std::cout << "fitting track with mom=" << mom << " pid=" << pid << " flip="
  // << flipDirection << " start pos=" << traj.Start() << " dir=" <<
  // traj.StartDirection() << std::endl;
  bool fitok = kalmanFitter.fitTrack(traj,
                                     tkID,
                                     covVtx,
                                     covEnd,
                                     inHits,
                                     mom,
                                     pid,
                                     flipDirection,
                                     outTrack,
                                     outHits,
                                     optionals);
  if (!fitok)
    return false;
  //
  // std::cout << "fitted track with mom=" << outTrack.StartMomentum() << "
  // pid=" << outTrack.ParticleId() << " flip=" << flipDirection << " start pos="
  // << outTrack.Start() << " dir=" << outTrack.StartDirection() << " nchi2=" <<
  // outTrack.Chi2PerNdof() << std::endl;
  //
  if (p_().options().keepInputTrajectoryPoints()) {
    restoreInputPoints(traj, inHits, outTrack, outHits, optionals);
  }
  //
  return true;
}

double
trkmkr::KalmanFilterFitTrackMaker::getMomentum(const recob::TrackTrajectory& traj,
                                               const int pid,
                                               const bool isFlip,
                                               const int tkID) const
{
  double mom = (mom_def_ >0 ? mom_def_ : traj.StartMomentum());
  if (momFromMCSColl_) {
    double mcsmom =
      (isFlip ? mcs->at(tkID).bwdMomentum() : mcs->at(tkID).fwdMomentum());
    // make sure the mcs fit converged
    if (mcsmom > 0.01 && mcsmom < 7.)
      mom = mcsmom;
  }
  if (momFromCombAndPid_) {
    bool isContained = cont->at(tkID).CosmicType() == anab::kNotTagged;
    // for now momentum from range implemented only for muons and protons
    // so treat pions as muons (~MIPs) and kaons as protons
    int pidtmp = pid;
    if (pidtmp == 211 || pidtmp == 13)
      pidtmp = 13;
    else
      pidtmp = 2212;
    mom = tmc.GetTrackMomentum(traj.Length(), pidtmp);
    if (isContained == false) {
      auto mcsresult = mcsfitter.fitMcs(traj, pid);
      double mcsmom =
        (isFlip ? mcsresult.bwdMomentum() : mcsresult.fwdMomentum());
      // make sure the mcs fit converged, also the mcsmom should not be less
      // than the range!
      if (mcsmom > 0.01 && mcsmom < 7. && mcsmom > mom)
        mom = mcsmom;
    }
  }
  return mom;
}

int
trkmkr::KalmanFilterFitTrackMaker::getParticleID(const recob::TrackTrajectory& traj,
                                                 const int tkID) const
{
  if (pidFromColl_) {
    return -1; //pid->at(tkID).Pdg();
  }
  if (mom_len_cut_ > 0.) {
    return (traj.Length() < mom_len_cut_ ? 2212 : 13);
  }
  return pid_def_;
}

bool
trkmkr::KalmanFilterFitTrackMaker::isFlipDirection(const recob::TrackTrajectory& traj,
                                                   const int tkID) const
{
  if (alwaysFlip_) {
    return true;
  } else if (dirFromVec_) {
    auto tdir = traj.VertexDirection();
    if (tdir.Dot(dirVec) < 0.)
      return true;
  }
  return false;
}

void
trkmkr::KalmanFilterFitTrackMaker::restoreInputPoints(const recob::TrackTrajectory& traj,
                                                      const std::vector<art::Ptr<recob::Hit>>& inHits,
                                                      recob::Track& outTrack,
                                                      std::vector<art::Ptr<recob::Hit>>& outHits,
                                                      OptionalOutputs& optionals) const
{
  if (optionals.isTrackFitInfosInit()) {
    throw cet::exception("KalmanFilterFitTrackMaker")
      << "Option keepInputTrajectoryPoints not compatible with "
         "doTrackFitHitInfo, please set doTrackFitHitInfo to false in the "
         "track producer.\n";
  }
  const auto np = outTrack.NumberTrajectoryPoints();
  trkmkr::TrackCreationBookKeeper tcbk(outHits,
                                       optionals,
                                       outTrack.ID(),
                                       outTrack.ParticleId(),
                                       traj.HasMomentum());
  //
  std::vector<unsigned int> flagsmap(np, -1);
  for (unsigned int i = 0; i < np; ++i)
    flagsmap[outTrack.FlagsAtPoint(i).fromHit()] = i;
  //
  for (unsigned int p = 0; p < np; ++p) {
    auto mask = outTrack.FlagsAtPoint(flagsmap[p]).mask();
    if (mask.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
      mask.unset(recob::TrajectoryPointFlagTraits::NoPoint);
    tcbk.addPoint(traj.LocationAtPoint(p),
                  traj.MomentumVectorAtPoint(p),
                  inHits[p],
                  recob::TrajectoryPointFlags(p, mask),
                  0);
  }
  auto covs = outTrack.Covariances();
  tcbk.setTotChi2(outTrack.Chi2());
  outTrack = tcbk.finalizeTrack(std::move(covs.first), std::move(covs.second));
  //
}

DEFINE_ART_CLASS_TOOL(trkmkr::KalmanFilterFitTrackMaker)
