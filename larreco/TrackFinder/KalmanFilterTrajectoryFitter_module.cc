/// \class KalmanFilterTrajectoryFitter
///
/// \brief Producer for fitting Trajectories and TrackTrajectories using TrackKalmanFitter.
///
/// \author G. Cerati
///

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larreco/RecoAlg/TrackKalmanFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/TrackFinder/TrackMaker.h"

#include "lardataobj/MCBase/MCTrack.h"

#include <memory>

namespace trkf {

  class KalmanFilterTrajectoryFitter : public art::EDProducer {
  public:

    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputTrajectoryLabel {
	Name("inputTrajectoryLabel"),
	Comment("Label of recob::Trajectory or recob::TrackTrajectory Collection to be fit")
      };
      fhicl::Atom<bool> isTrackTrajectory {
        Name("isTrackTrajectory"),
        Comment("If true, we assume the input collection is made of recob::TrackTrajectory objects, otherwise of recob::Trajectory objects.")
      };
      fhicl::Atom<art::InputTag> inputMCLabel {
	Name("inputMCLabel"),
	Comment("Label of sim::MCTrack Collection to be used for initial momentum estimate. Used only if momFromMC is set to true.")
      };
    };

    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      // fhicl::Atom<bool> pFromMSChi2 {
      //   Name("momFromMSChi2"),
      //   Comment("Flag used to get initial momentum estimate from trkf::TrackMomentumCalculator::GetMomentumMultiScatterChi2().")
      // };
      fhicl::Atom<bool> pFromLength {
        Name("momFromLength"),
        Comment("Flag used to get initial momentum estimate from trkf::TrackMomentumCalculator::GetTrackMomentum().")
      };
      fhicl::Atom<bool> pFromMC {
	Name("momFromMC"),
	Comment("Flag used to get initial momentum estimate from inputMCLabel collection.")
      };
      fhicl::Atom<double> pval {
	Name("momentumInGeV"),
	  Comment("Fixed momentum estimate value, to be used when momFromCalo, momFromLength and momFromMC are all false, or if the estimate is not available.")//momFromMSChi2,
      };
      fhicl::Atom<int> pdgId {
	Name("pdgId"),
        Comment("Default particle id hypothesis in case no valid id is provided either via PFParticle or in the ParticleId collection.")
      };
      fhicl::Atom<bool> dirFromMC {
        Name("dirFromMC"),
        Comment("Assume track direction from MC.")
      };
      fhicl::Atom<bool> dirFromVec {
        Name("dirFromVec"),
        Comment("Assume track direction from as the one giving positive dot product with vector specified by dirVec.")
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
      fhicl::Atom<bool> produceTrackFitHitInfo {
        Name("produceTrackFitHitInfo"),
        Comment("Option to produce (or not) the detailed TrackFitHitInfo.")
      };
      fhicl::Atom<bool> keepInputTrajectoryPoints {
        Name("keepInputTrajectoryPoints"),
        Comment("Option to keep positions and directions from input trajectory. The fit will provide only covariance matrices, chi2, ndof, particle Id and absolute momentum. It may also modify the trajectory point flags. In order to avoid inconsistencies, it has to be used with the following fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp.")
      };
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<KalmanFilterTrajectoryFitter::Inputs> inputs {
	Name("inputs"),
      };
      fhicl::Table<KalmanFilterTrajectoryFitter::Options> options {
	Name("options")
      };
      fhicl::Table<TrackStatePropagator::Config> propagator {
	Name("propagator")
      };
      fhicl::Table<TrackKalmanFitter::Config> fitter {
	Name("fitter")
      };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit KalmanFilterTrajectoryFitter(Parameters const & p);
    ~KalmanFilterTrajectoryFitter();

    // Plugins should not be copied or assigned.
    KalmanFilterTrajectoryFitter(KalmanFilterTrajectoryFitter const &) = delete;
    KalmanFilterTrajectoryFitter(KalmanFilterTrajectoryFitter &&) = delete;
    KalmanFilterTrajectoryFitter & operator = (KalmanFilterTrajectoryFitter const &) = delete;
    KalmanFilterTrajectoryFitter & operator = (KalmanFilterTrajectoryFitter &&) = delete;

    void produce(art::Event & e) override;

  private:
    Parameters p_;
    trkf::TrackKalmanFitter* kalmanFitter;
    TrackStatePropagator* prop;
    trkf::TrackMomentumCalculator* tmc;

    art::InputTag trajectoryInputTag;
    art::InputTag simTrackInputTag;

    bool isTT;

    double setMomValue(const recob::TrackTrajectory* ptraj, const double pMC, const int pId) const;
    int    setPId() const;
    bool   setDirFlip(const recob::TrackTrajectory* ptraj, TVector3& mcdir) const;

    void restoreInputPoints(const recob::TrackTrajectory& track,const std::vector<art::Ptr<recob::Hit> >& inHits,recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits) const;
  };
}

trkf::KalmanFilterTrajectoryFitter::KalmanFilterTrajectoryFitter(trkf::KalmanFilterTrajectoryFitter::Parameters const & p)
  : p_(p)
{

  prop = new TrackStatePropagator(p_().propagator);
  kalmanFitter = new trkf::TrackKalmanFitter(prop,p_().fitter);
  tmc = new trkf::TrackMomentumCalculator();

  trajectoryInputTag = art::InputTag(p_().inputs().inputTrajectoryLabel());
  if (p_().options().pFromMC() || p_().options().dirFromMC()) simTrackInputTag = art::InputTag(p_().inputs().inputMCLabel());

  isTT = p_().inputs().isTrackTrajectory();

  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  // produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();
  if (isTT) {
    produces<art::Assns<recob::TrackTrajectory, recob::Track> >();
  } else {
    produces<art::Assns<recob::Trajectory, recob::Track> >();
  }
  if (p_().options().produceTrackFitHitInfo()) {
    produces<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  }

  //throw expections to avoid possible silent failures due to incompatible configuration options

  unsigned int nDirs = 0;
  if (p_().options().dirFromMC())       nDirs++;
  if (p_().options().dirFromVec())      nDirs++;
  if (p_().options().alwaysInvertDir()) nDirs++;
  if (nDirs>1) {
    throw cet::exception("KalmanFilterTrajectoryFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among dirFromMC, dirFromVec, and alwaysInvertDir." << "\n";
  }

  unsigned int nPFroms = 0;
  // if (p_().options().pFromMSChi2()) nPFroms++;
  if (p_().options().pFromLength()) nPFroms++;
  if (p_().options().pFromMC())     nPFroms++;
  if (nPFroms>1) {
    throw cet::exception("KalmanFilterTrajectoryFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among pFromLength, and pFromMC." << "\n";//pFromMSChi2,
  }

  if (p_().options().keepInputTrajectoryPoints()) {
    if (p_().fitter().sortHitsByPlane() || p_().fitter().sortOutputHitsMinLength() || p_().fitter().skipNegProp()) {
      throw cet::exception("KalmanFilterTrajectoryFitter")
	<< "Incompatible configuration parameters: keepInputTrajectoryPoints needs the following fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp." << "\n";
    }
  }
}

trkf::KalmanFilterTrajectoryFitter::~KalmanFilterTrajectoryFitter() {
  delete prop;
  delete kalmanFitter;
  delete tmc;
}

void trkf::KalmanFilterTrajectoryFitter::produce(art::Event & e)
{

  auto outputTracks  = std::make_unique<std::vector<recob::Track> >();
  // auto outputHits    = std::make_unique<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();
  auto outputHits    = std::make_unique<art::Assns<recob::Track, recob::Hit> >();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo> > >();

  //only one will be filled and pushed into the event:
  auto outputTTjTAssn = std::make_unique<art::Assns<recob::TrackTrajectory, recob::Track> >();
  auto outputTjTAssn  = std::make_unique<art::Assns<recob::Trajectory     , recob::Track> >();

  auto const tid = getProductID<std::vector<recob::Track> >(e);
  auto const tidgetter = e.productGetter(tid);

  //FIXME, eventually remove this (ok only for single particle MC)
  double pMC = -1.;
  TVector3 mcdir;
  if (p_().options().pFromMC() || p_().options().dirFromMC()) {
    art::ValidHandle<std::vector<sim::MCTrack> > simTracks = e.getValidHandle<std::vector<sim::MCTrack> >(simTrackInputTag);
    for (unsigned int iMC = 0; iMC < simTracks->size(); ++iMC) {
      const sim::MCTrack& mctrack = simTracks->at(iMC);
      //fiducial cuts on MC tracks
      if (mctrack.PdgCode()!=13)   continue;
      if (mctrack.Process()!="primary")   continue;
      pMC = mctrack.Start().Momentum().P()*0.001;
      mcdir = TVector3(mctrack.Start().Momentum().X()*0.001/pMC,mctrack.Start().Momentum().Y()*0.001/pMC,mctrack.Start().Momentum().Z()*0.001/pMC);
      break;
    }
    //std::cout << "mc momentum value = " << pval << " GeV" << std::endl;
  }

  unsigned int nTrajs = 0;

  art::Handle<std::vector<recob::TrackTrajectory> > inputTrackTrajectoryH;
  art::Handle<std::vector<recob::Trajectory     > > inputTrajectoryH;
  const std::vector<recob::TrackTrajectory>* trackTrajectoryVec = 0;
  const std::vector<recob::Trajectory     >* trajectoryVec      = 0;
  const art::Assns<recob::TrackTrajectory, recob::Hit>* trackTrajectoryHitsAssn = 0;
  const art::Assns<recob::Trajectory     , recob::Hit>* trajectoryHitsAssn      = 0;
  if (isTT) {
    bool ok = e.getByLabel(trajectoryInputTag,inputTrackTrajectoryH);
    if (!ok) throw cet::exception("KalmanFilterTrajectoryFitter") << "Cannot find recob::TrackTrajectory art::Handle with inputTag " << trajectoryInputTag;
    trackTrajectoryVec = inputTrackTrajectoryH.product();
    trackTrajectoryHitsAssn = e.getValidHandle<art::Assns<recob::TrackTrajectory, recob::Hit> >(trajectoryInputTag).product();
    nTrajs = trackTrajectoryVec->size();
  } else {
    bool ok = e.getByLabel(trajectoryInputTag,inputTrajectoryH);
    if (!ok) throw cet::exception("KalmanFilterTrajectoryFitter") << "Cannot find recob::Trajectory art::Handle with inputTag " << trajectoryInputTag;
    trajectoryVec = inputTrajectoryH.product();
    trajectoryHitsAssn = e.getValidHandle<art::Assns<recob::Trajectory, recob::Hit> >(trajectoryInputTag).product();
    nTrajs = trajectoryVec->size();
  }

  for (unsigned int iTraj = 0; iTraj < nTrajs; ++iTraj) {

    const recob::TrackTrajectory& inTraj = (isTT ? trackTrajectoryVec->at(iTraj) : recob::TrackTrajectory(trajectoryVec->at(iTraj), std::vector<recob::TrajectoryPointFlags>()) );
    //const std::vector<recob::TrajectoryPointFlags>& inFlags = (isTT ? trackTrajectoryVec->at(iTraj).Flags() : std::vector<recob::TrajectoryPointFlags>());
    //this is not computationally optimal, but at least preserves the order unlike FindManyP
    std::vector<art::Ptr<recob::Hit> > inHits;
    if (isTT) {
      for (auto it = trackTrajectoryHitsAssn->begin(); it!=trackTrajectoryHitsAssn->end(); ++it) {
	if (it->first.key() == iTraj) inHits.push_back(it->second);
	else if (inHits.size()>0) break;
      }
    } else {
      for (auto it = trajectoryHitsAssn->begin(); it!=trajectoryHitsAssn->end(); ++it) {
	if (it->first.key() == iTraj) inHits.push_back(it->second);
	else if (inHits.size()>0) break;
      }
    }
    const int pId = setPId();
    const double mom = setMomValue(&inTraj, pMC, pId);
    const bool flipDir = setDirFlip(&inTraj, mcdir);

    recob::Track outTrack;
    std::vector<art::Ptr<recob::Hit> > outHits;
    trkmkr::OptionalOutputs optionals;
    if (p_().options().produceTrackFitHitInfo()) optionals.initTrackFitInfos();
    bool fitok = kalmanFitter->fitTrack(inTraj,iTraj,
					SMatrixSym55(),SMatrixSym55(),
					inHits,//inFlags,
					mom, pId, flipDir, outTrack, outHits, optionals);
    if (!fitok && (kalmanFitter->getSkipNegProp() || kalmanFitter->getCleanZigzag()) && p_().options().tryNoSkipWhenFails()) {
      //ok try once more without skipping hits
      mf::LogWarning("KalmanFilterTrajectoryFitter") << "Try to recover with skipNegProp = false and cleanZigzag = false\n";
      kalmanFitter->setSkipNegProp(false);
      kalmanFitter->setCleanZigzag(false);
      fitok = kalmanFitter->fitTrack(inTraj,iTraj,
				     SMatrixSym55(),SMatrixSym55(),
				     inHits,//inFlags,
				     mom, pId, flipDir, outTrack, outHits, optionals);
      kalmanFitter->setSkipNegProp(p_().fitter().skipNegProp());
      kalmanFitter->setCleanZigzag(p_().fitter().cleanZigzag());
    }
    if (!fitok) {
      mf::LogWarning("KalmanFilterTrajectoryFitter") << "Fit failed for track #" << iTraj << "\n";
      continue;
    }

    if (p_().options().keepInputTrajectoryPoints()) {
      restoreInputPoints(inTraj,inHits,outTrack,outHits);
    }

    outputTracks->emplace_back(std::move(outTrack));
    art::Ptr<recob::Track> aptr(tid, outputTracks->size()-1, tidgetter);
    unsigned int ip = 0;
    for (auto const& trhit: outHits) {
      //the fitter produces collections with 1-1 match between hits and point
      // recob::TrackHitMeta metadata(ip,-1);
      // outputHits->addSingle(aptr, trhit, metadata);
      outputHits->addSingle(aptr, trhit);
      ip++;
    }
    outputHitInfo->emplace_back(std::move(*optionals.trackFitHitInfos()));
    if (isTT) {
      outputTTjTAssn->addSingle(art::Ptr<recob::TrackTrajectory>(inputTrackTrajectoryH, iTraj),aptr);
    } else {
      outputTjTAssn->addSingle(art::Ptr<recob::Trajectory>(inputTrajectoryH, iTraj),aptr);
    }
  }
  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));
  if (p_().options().produceTrackFitHitInfo()) {
    e.put(std::move(outputHitInfo));
  }
  if (isTT) e.put(std::move(outputTTjTAssn));
  else e.put(std::move(outputTjTAssn));
}

void trkf::KalmanFilterTrajectoryFitter::restoreInputPoints(const recob::TrackTrajectory& track,const std::vector<art::Ptr<recob::Hit> >& inHits,recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits) const {
  const auto np = outTrack.NumberTrajectoryPoints();
  std::vector<Point_t>                     positions(np);
  std::vector<Vector_t>                    momenta(np);
  std::vector<recob::TrajectoryPointFlags> outFlags(np);
  //
  for (unsigned int p=0; p<np; ++p) {
    auto flag = outTrack.FlagsAtPoint(p);
    auto mom = outTrack.VertexMomentum();
    auto op = flag.fromHit();
    positions[op] = track.LocationAtPoint(op);
    momenta[op] = mom*track.DirectionAtPoint(op);
    auto mask = flag.mask();
    if (mask.isSet(recob::TrajectoryPointFlagTraits::NoPoint)) mask.unset(recob::TrajectoryPointFlagTraits::NoPoint);
    outFlags[op] = recob::TrajectoryPointFlags(op,mask);
  }
  auto covs = outTrack.Covariances();
  outTrack = recob::Track(recob::TrackTrajectory(std::move(positions),std::move(momenta),std::move(outFlags),true),
			  outTrack.ParticleId(),outTrack.Chi2(),outTrack.Ndof(),std::move(covs.first),std::move(covs.second),outTrack.ID());
  //
  outHits.clear();
  for (auto h : inHits) outHits.push_back(h);
}

double trkf::KalmanFilterTrajectoryFitter::setMomValue(const recob::TrackTrajectory* ptraj, const double pMC, const int pId) const {
  double result = p_().options().pval();
  // if (p_().options().pFromMSChi2()) {
  //   result = tmc->GetMomentumMultiScatterChi2(ptrack);
  // } else
  if (p_().options().pFromLength()) {
    result = tmc->GetTrackMomentum(ptraj->Length(), pId);
  } else if (p_().options().pFromMC() && pMC>0.) {
    result = pMC;
  }
  return result;
}

int trkf::KalmanFilterTrajectoryFitter::setPId() const {
  int result = p_().options().pdgId();
  return result;
}

bool trkf::KalmanFilterTrajectoryFitter::setDirFlip(const recob::TrackTrajectory* ptraj, TVector3& mcdir) const {
  bool result = false;
  if (p_().options().alwaysInvertDir()) {
    return true;
  } else if (p_().options().dirFromMC()) {
    auto tdir =  ptraj->VertexDirection();
    if ( (mcdir.X()*tdir.X() + mcdir.Y()*tdir.Y() + mcdir.Z()*tdir.Z())<0. ) result = true;
  } else if (p_().options().dirFromVec()) {
    std::array<float, 3> dir = p_().options().dirVec();
    auto tdir =  ptraj->VertexDirection();
    if ( (dir[0]*tdir.X() + dir[1]*tdir.Y() + dir[2]*tdir.Z())<0. ) result = true;
  }
  return result;
}

DEFINE_ART_MODULE(trkf::KalmanFilterTrajectoryFitter)
