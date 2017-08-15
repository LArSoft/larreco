#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/TrackFinder/TrackMaker.h"
#include "larreco/RecoAlg/TrackKalmanFitter.h"
#include "larreco/RecoAlg/TrackCreationBookKeeper.h"
//#include "lardataobj/RecoBase/MCSFitResult.h"

namespace trkmkr {

  class KalmanFilterFitTrackMaker : public TrackMaker {

  public:

    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      // fhicl::Atom<bool> pFromMSChi2 {
      //   Name("momFromMSChi2"),
      //   Comment("Flag used to get initial momentum estimate from trkf::TrackMomentumCalculator::GetMomentumMultiScatterChi2().")
      // };
      // fhicl::Atom<bool> pFromLength {
      //   Name("momFromLength"),
      //   Comment("Flag used to get initial momentum estimate from trkf::TrackMomentumCalculator::GetTrackMomentum().")
      // };
      // fhicl::Atom<bool> pFromCalo {
      // 	Name("momFromCalo"),
      // 	Comment("Flag used to get initial momentum estimate from inputCaloLabel collection.")
      // };
      fhicl::Atom<double> defaultMomInGeV {
	Name("defaultMomInGeV"),
	Comment("Default momentum estimate value (all other options are set to false, or if the estimate is not available).")
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
      fhicl::Atom<bool> tryNoSkipWhenFails {
        Name("tryNoSkipWhenFails"),
        Comment("In case skipNegProp is true and the track fit fails, make a second attempt to fit the track with skipNegProp=false in order to attempt to avoid losing efficiency.")
      };
      fhicl::Atom<bool> keepInputTrajectoryPoints {
        Name("keepInputTrajectoryPoints"),
        Comment("Option to keep positions and directions from input trajectory/track. The fit will provide only covariance matrices, chi2, ndof, particle Id and absolute momentum. It may also modify the trajectory point flags. In order to avoid inconsistencies, it has to be used with the following fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp.")
      };
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
      pid_def_ = p_().options().defaultPdgId();
      alwaysFlip_ = p_().options().alwaysInvertDir();
      //
      if (p_().options().keepInputTrajectoryPoints()) {
	if (p_().fitter().sortHitsByPlane() || p_().fitter().sortOutputHitsMinLength() || p_().fitter().skipNegProp()) {
	  throw cet::exception("KalmanFilterFitTrackMaker")
	    << "Incompatible configuration parameters: keepInputTrajectoryPoints needs the following fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp." << "\n";
	}
      }
      //
    }

    ~KalmanFilterFitTrackMaker() {
      delete prop;
      delete kalmanFitter;
    }

    bool makeTrack(const recob::TrackTrajectory& traj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
		   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals,
		   const art::Event& e) const override;

    bool setDirFlip(const recob::TrackTrajectory& traj) const;

    void restoreInputPoints(const recob::TrackTrajectory& traj, const std::vector<art::Ptr<recob::Hit> >& inHits,
			    recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals) const;

  private:
    Parameters p_;
    trkf::TrackKalmanFitter* kalmanFitter;
    trkf::TrackStatePropagator* prop;
    double mom_def_;
    int    pid_def_;
    bool   alwaysFlip_;
  };

}

bool trkmkr::KalmanFilterFitTrackMaker::makeTrack(const recob::TrackTrajectory& traj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
						  recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals,
						  const art::Event& e) const {
  //
  // art::ValidHandle<std::vector<recob::MCSFitResult> > mcs = e.getValidHandle<std::vector<recob::MCSFitResult> >("mcsproducer");
  // std::cout << "bestMom=" << mcs->at(tkID).bestMomentum() << " isBestFwd=" << mcs->at(tkID).isBestFwd() << std::endl;
  //
  const double mom = mom_def_;//mcs->at(tkID).bestMomentum();//mom_def_;//what about uncertainty? what about fit failures?
  const int    pid = pid_def_;
  const bool   flipDirection = setDirFlip(traj);//(mcs->at(tkID).isBestFwd()==false);//setDirFlip(traj);
  bool fitok = kalmanFitter->fitTrack(traj, tkID, trkf::SMatrixSym55(), trkf::SMatrixSym55(), inHits, mom, pid, flipDirection, outTrack, outHits, optionals);
  //
  if (!fitok && (kalmanFitter->getSkipNegProp() || kalmanFitter->getCleanZigzag()) && p_().options().tryNoSkipWhenFails()) {
    //ok try once more without skipping hits
    mf::LogWarning("KalmanFilterFitTrackMaker") << "Try to recover with skipNegProp = false and cleanZigzag = false\n";
    kalmanFitter->setSkipNegProp(false);
    kalmanFitter->setCleanZigzag(false);
    fitok = kalmanFitter->fitTrack(traj, tkID, trkf::SMatrixSym55(), trkf::SMatrixSym55(), inHits, mom, pid, flipDirection, outTrack, outHits, optionals);
    kalmanFitter->setSkipNegProp(p_().fitter().skipNegProp());
    kalmanFitter->setCleanZigzag(p_().fitter().cleanZigzag());
  }
  if (!fitok) {
    mf::LogWarning("KalmanFilterFitTrackMaker") << "Fit failed for track with ID=" << tkID << "\n";
    return false;
  }
  //
  if (p_().options().keepInputTrajectoryPoints()) {
    restoreInputPoints(traj,inHits,outTrack,outHits,optionals);
  }
  //
  return true;
}

bool trkmkr::KalmanFilterFitTrackMaker::setDirFlip(const recob::TrackTrajectory& traj) const {
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
  trkmkr::TrackCreationBookKeeper tcbk(outTrack, outHits, optionals, outTrack.ID(), outTrack.ParticleId(), traj.HasMomentum());
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
  tcbk.finalizeTrack(std::move(covs.first),std::move(covs.second));
  //
}

DEFINE_ART_CLASS_TOOL(trkmkr::KalmanFilterFitTrackMaker)
