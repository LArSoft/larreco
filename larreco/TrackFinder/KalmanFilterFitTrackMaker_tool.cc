////////////////////////////////////////////////////////////////////////
// Class:       KalmanFilterFitTrackMaker
// File:        KalmanFilterFitTrackMaker_tool.cc
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/TrackFinder/TrackMaker.h"
#include "larreco/RecoAlg/TrackKalmanFitter.h"
#include "larreco/RecoAlg/TrackCreationBookKeeper.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

namespace trkmkr {

  /**
   * @brief Concrete implementation of a tool to fit tracks with TrackKalmanFitter.
   *
   * Concrete implementation of a tool to fit tracks with TrackKalmanFitter; inherits from abstract class TrackMaker.
   * It prepares the input needed by the fitter (momentum, particleId, direction), and returns a track with all outputs filled.
   * If the flag keepInputTrajectoryPoints is set to true, the tracjetory points from the input track are copied into the output,
   * so that only the covariance matrices, the chi2 and the ndof in the output track are resulting from the fit.
   */

  class KalmanFilterFitTrackMaker : public TrackMaker {

  public:

    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      //
      fhicl::Atom<double> defaultMomInGeV {
	Name("defaultMomInGeV"),
	Comment("Default momentum estimate value (all other options are set to false, or if the estimate is not available).")
      };
      fhicl::Atom<bool> momFromMCS {
        Name("momFromMCS"),
        Comment("Flag used to get initial momentum estimate from MCSFitResult.")
      };
      fhicl::Atom<art::InputTag> mcsInputTag {
        Name("mcsInputTag"),
        Comment("InputTag of MCSFitResult collection.")
      };
      fhicl::Atom<int> defaultPdgId {
	Name("defaultPdgId"),
        Comment("Default particle id hypothesis (all other options are set to false, or if the estimate is not available).")
      };
      fhicl::Atom<bool> dirFromVec {
        Name("dirFromVec"),
        Comment("Assume track direction as the one giving positive dot product with vector specified by dirVec.")
      };
      fhicl::Sequence<float,3u> dirVec {
        Name("dirVec"),
        Comment("Fhicl sequence defining the vector used when dirFromVec=true. It must have 3 elements.")
      };
      fhicl::Atom<bool> alwaysInvertDir {
	Name("alwaysInvertDir"),
	Comment("If true, fit all tracks from end to vertex assuming inverted direction.")
      };
      fhicl::Atom<bool> keepInputTrajectoryPoints {
        Name("keepInputTrajectoryPoints"),
        Comment("Option to keep positions and directions from input trajectory/track. The fit will provide only covariance matrices, chi2, ndof, particle Id and absolute momentum. It may also modify the trajectory point flags. In order to avoid inconsistencies, it has to be used with the following fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp.")
      };
      //
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<KalmanFilterFitTrackMaker::Options> options {
	Name("options")
      };
      fhicl::Table<trkf::TrackStatePropagator::Config> propagator {
	Name("propagator")
      };
      fhicl::Table<trkf::TrackKalmanFitter::Config> fitter {
	Name("fitter")
      };
    };
    using Parameters = art::ToolConfigTable<Config>;

    explicit KalmanFilterFitTrackMaker(Parameters const& p) : p_(p)
    {
      prop = new trkf::TrackStatePropagator(p_().propagator);
      kalmanFitter = new trkf::TrackKalmanFitter(prop,p_().fitter);
      mom_def_ = p_().options().defaultMomInGeV();
      momFromMCS_ = p_().options().momFromMCS();
      mcsInputTag_ = p_().options().mcsInputTag();
      pid_def_ = p_().options().defaultPdgId();
      alwaysFlip_ = p_().options().alwaysInvertDir();
      //
      if ( p_().options().keepInputTrajectoryPoints() && (p_().fitter().sortHitsByPlane() || p_().fitter().sortOutputHitsMinLength() || p_().fitter().skipNegProp()) ) {
	throw cet::exception("KalmanFilterFitTrackMaker")
	  << "Incompatible configuration parameters: keepInputTrajectoryPoints needs the following fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp." << "\n";
      }
      //
    }

    ~KalmanFilterFitTrackMaker() {
      delete prop;
      delete kalmanFitter;
    }

    virtual void initEvent(const art::Event& e) override {
      if (momFromMCS_) mcs = e.getValidHandle<std::vector<recob::MCSFitResult> >(mcsInputTag_).product();
      return;
    }

    bool makeTrackImpl(const recob::TrackTrajectory& traj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
		       const SMatrixSym55& covVtx, const SMatrixSym55& covEnd,
		       recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const;

    bool makeTrack(const recob::TrackTrajectory& traj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
		   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const override {
      return makeTrackImpl(traj, tkID, inHits, trkf::SMatrixSym55(), trkf::SMatrixSym55(), outTrack, outHits, optionals);
    }

    bool makeTrack(const recob::Track& track, const std::vector<art::Ptr<recob::Hit> >& inHits,
		   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const override {
      auto covs = track.Covariances();
      return makeTrackImpl(track.Trajectory(), track.ID(), inHits, covs.first, covs.second, outTrack, outHits, optionals);
    }

    double setMomentum  (const recob::TrackTrajectory& traj, const int tkID) const;
    int    setParticleID(const recob::TrackTrajectory& traj, const int tkID) const;
    bool   setDirection (const recob::TrackTrajectory& traj, const int tkID) const;

    void restoreInputPoints(const recob::TrackTrajectory& traj, const std::vector<art::Ptr<recob::Hit> >& inHits,
			    recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const;

  private:
    Parameters p_;
    trkf::TrackKalmanFitter* kalmanFitter;
    trkf::TrackStatePropagator* prop;
    double mom_def_;
    bool   momFromMCS_;
    art::InputTag mcsInputTag_;
    int    pid_def_;
    bool   alwaysFlip_;
    const std::vector<recob::MCSFitResult>* mcs = nullptr;
  };

}

bool trkmkr::KalmanFilterFitTrackMaker::makeTrackImpl(const recob::TrackTrajectory& traj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
						      const SMatrixSym55& covVtx, const SMatrixSym55& covEnd,
						      recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const {
  //
  const double mom = setMomentum(traj, tkID);// what about uncertainty?
  const int    pid = setParticleID(traj, tkID);
  const bool   flipDirection = setDirection(traj, tkID);
  bool fitok = kalmanFitter->fitTrack(traj, tkID, covVtx, covEnd, inHits, mom, pid, flipDirection, outTrack, outHits, optionals);
  if (!fitok) return false;
  //
  if (p_().options().keepInputTrajectoryPoints()) {
    restoreInputPoints(traj,inHits,outTrack,outHits,optionals);
  }
  //
  return true;
}

double trkmkr::KalmanFilterFitTrackMaker::setMomentum(const recob::TrackTrajectory& traj, const int tkID) const {
  double mom = mom_def_;
  if (momFromMCS_) {
    double mcsmom = mcs->at(tkID).bestMomentum();
    //make sure the mcs fit converged
    if (mcsmom>0.01 && mcsmom<7.) mom = mcsmom;
  }
  return mom;
}

int trkmkr::KalmanFilterFitTrackMaker::setParticleID(const recob::TrackTrajectory& traj, const int tkID) const {
  return pid_def_;
}

bool trkmkr::KalmanFilterFitTrackMaker::setDirection(const recob::TrackTrajectory& traj, const int tkID) const {
  if (alwaysFlip_) {
    return true;
  } else if (p_().options().dirFromVec()) {
    std::array<float, 3> dir = p_().options().dirVec();
    auto tdir =  traj.VertexDirection();
    if ( (dir[0]*tdir.X() + dir[1]*tdir.Y() + dir[2]*tdir.Z())<0. ) return true;
  }
  return false;
}

void trkmkr::KalmanFilterFitTrackMaker::restoreInputPoints(const recob::TrackTrajectory& traj, const std::vector<art::Ptr<recob::Hit> >& inHits,
							   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const {
  if (optionals.isTrackFitInfosInit()) {
    throw cet::exception("KalmanFilterFitTrackMaker")
      << "Option keepInputTrajectoryPoints not compatible with doTrackFitHitInfo, please set doTrackFitHitInfo to false in the track producer.\n";
  }
  const auto np = outTrack.NumberTrajectoryPoints();
  trkmkr::TrackCreationBookKeeper tcbk(outHits, optionals, outTrack.ID(), outTrack.ParticleId(), traj.HasMomentum());
  //
  std::vector<unsigned int> flagsmap(np,-1);
  for (unsigned int i=0; i<np; ++i) flagsmap[outTrack.FlagsAtPoint(i).fromHit()] = i;
  //
  for (unsigned int p=0; p<np; ++p) {
    auto mask = outTrack.FlagsAtPoint(flagsmap[p]).mask();
    if (mask.isSet(recob::TrajectoryPointFlagTraits::NoPoint)) mask.unset(recob::TrajectoryPointFlagTraits::NoPoint);
    tcbk.addPoint(Point_t(traj.LocationAtPoint(p)), Vector_t(traj.MomentumVectorAtPoint(p)), inHits[p], recob::TrajectoryPointFlags(p,mask), 0);
  }
  auto covs = outTrack.Covariances();
  tcbk.setTotChi2(outTrack.Chi2());
  outTrack = tcbk.finalizeTrack(std::move(covs.first),std::move(covs.second));
  //
}

DEFINE_ART_CLASS_TOOL(trkmkr::KalmanFilterFitTrackMaker)
