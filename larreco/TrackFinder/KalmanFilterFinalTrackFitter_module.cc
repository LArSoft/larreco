////////////////////////////////////////////////////////////////////////
// Class:       KalmanFilterFinalTrackFitter
// Module Type: producer
// File:        KalmanFilterFinalTrackFitter_module.cc
//
// Generated at Fri Sep  2 10:48:46 2016 by Giuseppe Cerati using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larreco/RecoAlg/TrackKalmanFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardata/RecoObjects/PropYZPlane.h"

#include "lardataobj/MCBase/MCTrack.h"

#include <memory>

namespace trkf {

  class KalmanFilterFinalTrackFitter : public art::EDProducer {
  public:

    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputPFParticleLabel {
       Name("inputPFParticleLabel"),
       Comment("Label of recob::PFParticle Collection to be fit")
      };
      fhicl::Atom<art::InputTag> inputTracksLabel {
	Name("inputTracksLabel"),
	Comment("Label of recob::Track Collection to be fit")
      };
      fhicl::Atom<art::InputTag> inputCaloLabel {
	Name("inputCaloLabel"),
	Comment("Label of anab::Calorimetry Collection, matching inputTracksLabel, to be used for initial momentum estimate. Used only if momFromCalo is set to true.")
      };
      fhicl::Atom<art::InputTag> inputMCLabel {
	Name("inputMCLabel"),
	Comment("Label of sim::MCTrack Collection to be used for initial momentum estimate. Used only if momFromMC is set to true.")
      };
      fhicl::Atom<art::InputTag> inputPidLabel {
       Name("inputPidLabel"),
       Comment("Label of anab::ParticleID Collection, matching inputTracksLabel, to be used for particle Id.")
      };
    };

    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> trackFromPF {
        Name("trackFromPF"),
        Comment("If true extract tracks from inputPFParticleLabel collection, if false from inputTracksLabel.")
      };
      fhicl::Atom<bool> pFromMSChi2 {
        Name("momFromMSChi2"),
        Comment("Flag used to get initial momentum estimate from trkf::TrackMomentumCalculator::GetMomentumMultiScatterChi2().")
      };
      fhicl::Atom<bool> pFromLength {
        Name("momFromLength"),
        Comment("Flag used to get initial momentum estimate from trkf::TrackMomentumCalculator::GetTrackMomentum().")
      };
      fhicl::Atom<bool> pFromCalo {
	Name("momFromCalo"),
	Comment("Flag used to get initial momentum estimate from inputCaloLabel collection.")
      };
      fhicl::Atom<bool> pFromMC {
	Name("momFromMC"),
	Comment("Flag used to get initial momentum estimate from inputMCLabel collection.")
      };
      fhicl::Atom<double> pval {
	Name("momentumInGeV"),
	Comment("Fixed momentum estimate value, to be used when momFromCalo, momFromMSChi2, momFromLength and momFromMC are all false, or if the estimate is not available.")
      };
      fhicl::Atom<bool> idFromPF {
        Name("idFromPF"),
        Comment("Flag used to get particle ID estimate from corresponding PFParticle. Needs trackFromPF=true.")
      };
      fhicl::Atom<bool> idFromCollection {
        Name("idFromCollection"),
        Comment("Flag used to get particle ID estimate from inputPidLabel collection.")
      };
      fhicl::Atom<int> pdgId {
	Name("pdgId"),
        Comment("Default particle id hypothesis in case no valid id is provided either via PFParticle or in the ParticleId collection.")
      };
      fhicl::Atom<bool> dirFromVtxPF {
        Name("dirFromVtxPF"),
        Comment("Assume track direction from Vertex in PFParticle. Needs trackFromPF=true.")
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
        Comment("Flag to decide whether, during the forward fit, the hits corresponding to a negative propagation distance should be dropped.")
      };
      fhicl::Atom<float> hitErrScaleFact {
        Name("hitErrScaleFact"),
        Comment(".")
      };
  };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<KalmanFilterFinalTrackFitter::Inputs> inputs {
	Name("inputs"),
      };
      fhicl::Table<KalmanFilterFinalTrackFitter::Options> options {
	Name("options")
      };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit KalmanFilterFinalTrackFitter(Parameters const & p);
    ~KalmanFilterFinalTrackFitter();

    // Plugins should not be copied or assigned.
    KalmanFilterFinalTrackFitter(KalmanFilterFinalTrackFitter const &) = delete;
    KalmanFilterFinalTrackFitter(KalmanFilterFinalTrackFitter &&) = delete;
    KalmanFilterFinalTrackFitter & operator = (KalmanFilterFinalTrackFitter const &) = delete;
    KalmanFilterFinalTrackFitter & operator = (KalmanFilterFinalTrackFitter &&) = delete;

    void produce(art::Event & e) override;

  private:
    Parameters p_;
    trkf::TrackKalmanFitter* kalmanFitter;
    trkf::Propagator* prop;
    trkf::TrackMomentumCalculator* tmc;

    art::InputTag pfParticleInputTag;
    art::InputTag trackInputTag;
    art::InputTag caloInputTag;
    art::InputTag pidInputTag;
    art::InputTag simTrackInputTag;

    std::unique_ptr<art::FindManyP<anab::Calorimetry> > trackCalo;
    std::unique_ptr<art::FindManyP<anab::ParticleID> > trackId;
    std::unique_ptr<art::FindManyP<recob::Track> > assocTracks;
    std::unique_ptr<art::FindManyP<recob::Vertex> > assocVertices;

    double setMomValue(art::Ptr<recob::Track> ptrack, const std::unique_ptr<art::FindManyP<anab::Calorimetry> >& trackCalo, const double pMC, const int pId) const;
    int    setPId(const unsigned int iTrack, const std::unique_ptr<art::FindManyP<anab::ParticleID> >& trackId, const int pfPid = 0) const;
    bool   setDirFlip(const recob::Track& track, TVector3& mcdir, const std::vector<art::Ptr<recob::Vertex> >* vertices = 0) const;
  };
}

trkf::KalmanFilterFinalTrackFitter::KalmanFilterFinalTrackFitter(trkf::KalmanFilterFinalTrackFitter::Parameters const & p)
  : p_(p)
{

  prop = new trkf::PropYZPlane(0., true);
  kalmanFitter = new trkf::TrackKalmanFitter(prop,p_().options().useRMS(),
					     p_().options().sortHitsByPlane(),
					     p_().options().sortOutputHitsMinLength(),
					     p_().options().skipNegProp(),
					     p_().options().hitErrScaleFact());
  tmc = new trkf::TrackMomentumCalculator();

  if (p_().options().trackFromPF()) pfParticleInputTag = art::InputTag(p_().inputs().inputPFParticleLabel());
  else {
    trackInputTag = art::InputTag(p_().inputs().inputTracksLabel());
    if (p_().options().idFromCollection()) pidInputTag = art::InputTag(p_().inputs().inputPidLabel());
  }
  if (p_().options().pFromCalo()) caloInputTag = art::InputTag(p_().inputs().inputCaloLabel());
  if (p_().options().pFromMC() || p_().options().dirFromMC()) simTrackInputTag = art::InputTag(p_().inputs().inputMCLabel());

  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  if (p_().options().trackFromPF()) {
    produces<art::Assns<recob::PFParticle, recob::Track> >();
  }

  //throw expections to avoid possible silent failures due to incompatible configuration options
  if (p_().options().trackFromPF()==0 && p_().options().idFromPF())
    throw cet::exception("KalmanFilterFinalTrackFitter") << "Incompatible configuration parameters: cannot use idFromPF=true with trackFromPF=false." << "\n";
  if (p_().options().trackFromPF()==0 && p_().options().dirFromVtxPF())
    throw cet::exception("KalmanFilterFinalTrackFitter") << "Incompatible configuration parameters: cannot use dirFromVtxPF=true with trackFromPF=false." << "\n";

  unsigned int nIds = 0;
  if (p_().options().idFromPF())         nIds++;
  if (p_().options().idFromCollection()) nIds++;
  if (nIds>1) {
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among idFromPF and idFromCollection." << "\n";
  }

  unsigned int nDirs = 0;
  if (p_().options().dirFromVtxPF())    nDirs++;
  if (p_().options().dirFromMC())       nDirs++;
  if (p_().options().dirFromVec())      nDirs++;
  if (p_().options().alwaysInvertDir()) nDirs++;
  if (nDirs>1) {
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among dirFromVtxPF, dirFromMC, dirFromVec, and alwaysInvertDir." << "\n";
  }

  unsigned int nPFroms = 0;
  if (p_().options().pFromCalo())   nPFroms++;
  if (p_().options().pFromMSChi2()) nPFroms++;
  if (p_().options().pFromLength()) nPFroms++;
  if (p_().options().pFromMC())     nPFroms++;
  if (nPFroms>1) {
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among pFromCalo, pFromMSChi2, pFromLength, and pFromMC." << "\n";
  }
}

trkf::KalmanFilterFinalTrackFitter::~KalmanFilterFinalTrackFitter() {
  delete prop;
  delete kalmanFitter;
  delete tmc;
}

void trkf::KalmanFilterFinalTrackFitter::produce(art::Event & e)
{

  auto outputTracks = std::make_unique<std::vector<recob::Track> >();
  auto outputHits   = std::make_unique<art::Assns<recob::Track, recob::Hit> >();

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

  if (p_().options().trackFromPF()) {

    auto outputPFAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Track> >();

    art::ValidHandle<std::vector<recob::PFParticle> > inputPFParticle = e.getValidHandle<std::vector<recob::PFParticle> >(pfParticleInputTag);
    assocTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPFParticle, e, pfParticleInputTag));
    assocVertices = std::unique_ptr<art::FindManyP<recob::Vertex> >(new art::FindManyP<recob::Vertex>(inputPFParticle, e, pfParticleInputTag));
    std::vector<art::Ptr<recob::Track> > tracks;
    std::vector<art::Ptr<recob::Vertex> > vertices;

    for (unsigned int iPF = 0; iPF < inputPFParticle->size(); ++iPF) {

      tracks = assocTracks->at(iPF);
      auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(pfParticleInputTag);
      vertices = assocVertices->at(iPF);

      if (p_().options().pFromCalo()) {
	trackCalo = std::unique_ptr<art::FindManyP<anab::Calorimetry> >(new art::FindManyP<anab::Calorimetry>(tracks, e, caloInputTag));
      }

      for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack) {

	const recob::Track& track = *tracks[iTrack];
	art::Ptr<recob::Track> ptrack = tracks[iTrack];
	const int pId = setPId(iTrack, trackId, inputPFParticle->at(iPF).PdgCode());
	const double mom = setMomValue(ptrack, trackCalo, pMC, pId);
	const bool flipDir = setDirFlip(track, mcdir, &vertices);

	//this is not computationally optimal, but at least preserves the order unlike FindManyP
	std::vector<art::Ptr<recob::Hit> > inHits;
	for (auto it = tkHitsAssn.begin(); it!=tkHitsAssn.end(); ++it) {
	  if (it->first == ptrack) inHits.push_back(it->second);
	  else if (inHits.size()>0) break;
	}

	recob::Track outTrack;
	art::PtrVector<recob::Hit> outHits;
	bool fitok = kalmanFitter->fitTrack(track, inHits, mom, pId, flipDir, outTrack, outHits);

	if (!fitok) {
	  mf::LogWarning("KalmanFilterFinalTrackFitter") << "Fit failed for PFP # " << iPF << " track #" << iTrack << "\n";
	  continue;
	}

	outputTracks->emplace_back(std::move(outTrack));
	art::Ptr<recob::Track> aptr(tid, outputTracks->size()-1, tidgetter);
	for (auto const& trhit: outHits) {
	  outputHits->addSingle(aptr, trhit);
	}
	outputPFAssn->addSingle(art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr);
      }
    }
    e.put(std::move(outputTracks));
    e.put(std::move(outputHits));
    e.put(std::move(outputPFAssn));
  } else {

    art::ValidHandle<std::vector<recob::Track> > inputTracks = e.getValidHandle<std::vector<recob::Track> >(trackInputTag);
    auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(trackInputTag);

    if (p_().options().pFromCalo()) {
      trackCalo = std::unique_ptr<art::FindManyP<anab::Calorimetry> >(new art::FindManyP<anab::Calorimetry>(inputTracks, e, caloInputTag));
    }

    if (p_().options().idFromCollection()) {
      trackId = std::unique_ptr<art::FindManyP<anab::ParticleID> >(new art::FindManyP<anab::ParticleID>(inputTracks, e, pidInputTag));
    }

    for (unsigned int iTrack = 0; iTrack < inputTracks->size(); ++iTrack) {

      const recob::Track& track = inputTracks->at(iTrack);
      art::Ptr<recob::Track> ptrack(inputTracks, iTrack);
      const int pId = setPId(iTrack, trackId);
      const double mom = setMomValue(ptrack, trackCalo, pMC, pId);
      const bool flipDir = setDirFlip(track, mcdir);

      //this is not computationally optimal, but at least preserves the order unlike FindManyP
      std::vector<art::Ptr<recob::Hit> > inHits;
      for (auto it = tkHitsAssn.begin(); it!=tkHitsAssn.end(); ++it) {
	if (it->first == ptrack) inHits.push_back(it->second);
	else if (inHits.size()>0) break;
      }

      recob::Track outTrack;
      art::PtrVector<recob::Hit> outHits;
      bool fitok = kalmanFitter->fitTrack(track, inHits, mom, pId, flipDir, outTrack, outHits);

      if (!fitok) {
	mf::LogWarning("KalmanFilterFinalTrackFitter") << "Fit failed for track #" << iTrack << "\n";
	continue;
      }

      outputTracks->emplace_back(std::move(outTrack));
      art::Ptr<recob::Track> aptr(tid, outputTracks->size()-1, tidgetter);
      for (auto const& trhit: outHits) {
	outputHits->addSingle(aptr, trhit);
      }
    }
    e.put(std::move(outputTracks));
    e.put(std::move(outputHits));
  }
}

double trkf::KalmanFilterFinalTrackFitter::setMomValue(art::Ptr<recob::Track> ptrack, const std::unique_ptr<art::FindManyP<anab::Calorimetry> >& trackCalo, const double pMC, const int pId) const {
  double result = p_().options().pval();
  if (p_().options().pFromMSChi2()) {
    result = tmc->GetMomentumMultiScatterChi2(ptrack);
  } else if (p_().options().pFromLength()) {
    result = tmc->GetTrackMomentum(ptrack->Length(), pId);
  } else if (p_().options().pFromCalo()) {
    //take average energy from available views
    const std::vector<art::Ptr<anab::Calorimetry> >& calo = trackCalo->at(ptrack.key());
    double sumenergy = 0.;
    int nviews = 0.;
    for (auto caloit : calo) {
      if (caloit->KineticEnergy()>0.) {
	sumenergy+=caloit->KineticEnergy();
	nviews+=1;
      }
    }
    if (nviews!=0 && sumenergy!=0.) {
      //protect against problematic cases
      result = sumenergy/(nviews*1000.);
    }
  } else if (p_().options().pFromMC() && pMC>0.) {
    result = pMC;
  }
  return result;
}

int trkf::KalmanFilterFinalTrackFitter::setPId(const unsigned int iTrack, const std::unique_ptr<art::FindManyP<anab::ParticleID> >& trackId, const int pfPid) const {
  int result = p_().options().pdgId();
  if (p_().options().trackFromPF() && p_().options().idFromPF()) {
    result = pfPid;
  } else if (p_().options().idFromCollection()) {
    //take the pdgId corresponding to the minimum chi2 (should we give preference to the majority? fixme)
    double minChi2 = -1.;
    for (auto idit : trackId->at(iTrack)) {
      if ( idit->MinChi2()>0. && (minChi2<0. || idit->MinChi2()<minChi2) ) {
	result = idit->Pdg();
	minChi2 = idit->MinChi2();
      }
    }
  }
  return result;
}

bool trkf::KalmanFilterFinalTrackFitter::setDirFlip(const recob::Track& track, TVector3& mcdir, const std::vector<art::Ptr<recob::Vertex> >* vertices) const {
  bool result = false;
  if (p_().options().alwaysInvertDir()) {
    return true;
  } else if (p_().options().dirFromMC()) {
    auto tdir =  track.VertexDirection();
    if ( (mcdir.X()*tdir.X() + mcdir.Y()*tdir.Y() + mcdir.Z()*tdir.Z())<0. ) result = true;
  } else if (p_().options().dirFromVec()) {
    std::array<float, 3> dir = p_().options().dirVec();
    auto tdir =  track.VertexDirection();
    if ( (dir[0]*tdir.X() + dir[1]*tdir.Y() + dir[2]*tdir.Z())<0. ) result = true;
  } else if (p_().options().trackFromPF() && p_().options().dirFromVtxPF() && vertices->size()>0) {
    //if track end is closer to first vertex then track vertex, flip direction
    double xyz[3];
    (*vertices)[0]->XYZ(xyz);
    auto& tv = track.Trajectory().Vertex();
    auto& te = track.Trajectory().End();
    if ( ((xyz[0]-te.X())*(xyz[0]-te.X()) + (xyz[1]-te.Y())*(xyz[1]-te.Y()) + (xyz[2]-te.Z())*(xyz[2]-te.Z())) >
	 ((xyz[0]-tv.X())*(xyz[0]-tv.X()) + (xyz[1]-tv.Y())*(xyz[1]-tv.Y()) + (xyz[2]-tv.Z())*(xyz[2]-tv.Z())) ) result = true;
  }
  return result;
}

DEFINE_ART_MODULE(trkf::KalmanFilterFinalTrackFitter)
