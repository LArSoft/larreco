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

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "lardata/RecoObjects/KHitWireX.h"
#include "lardata/RecoObjects/KHitContainerWireX.h"
#include "lardata/RecoObjects/Surface.h"
#include "lardata/RecoObjects/KGTrack.h"
#include "lardata/RecoObjects/PropAny.h"
#include "lardata/RecoObjects/SurfYZPlane.h"
#include "lardata/RecoObjects/PropYZPlane.h"

#include "larreco/RecoAlg/TrackKalmanFitter.h"

#include "lardataobj/MCBase/MCTrack.h"

#include <memory>

namespace trkf {

  class KalmanFilterFinalTrackFitter : public art::EDProducer {
  public:

    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> inputTracksLabel {
	Name("inputTracksLabel"),
	Comment("Label of reco::Track Collection to be fit")
      };
      fhicl::Atom<art::InputTag> inputCaloLabel {
	Name("inputCaloLabel"),
	Comment("Label of anab::Calorimetry Collection, matching inputTracksLabel, to be used for initial momentum estimate. Used only if momFromCalo is set to true.")
      };
      fhicl::Atom<art::InputTag> inputMCLabel {
	Name("inputMCLabel"),
	Comment("Label of sim::MCTrack Collection to be used for initial momentum estimate. Used only if momFromMC is set to true and momFromCalo is set to false.")
      };
    };

    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<bool> pFromCalo {
	Name("momFromCalo"),
	Comment("Flag used to get initial momentum estimate from inputCaloLabel collection. Has precedence over momFromMC.")
      };
      fhicl::Atom<bool> pFromMC {
	Name("momFromMC"),
	Comment("Flag used to get initial momentum estimate from inputMCLabel collection. Bypassed if momFromCalo is true.")
      };
      fhicl::Atom<double> pval {
	Name("momentumInGeV"),
	Comment("Fixed momentum estimate value, to be used when momFromCalo and momFromMC are both false, or if the estimate is not available.")
      };
      fhicl::Atom<int> pdgId {
	Name("pdgId"),
	Comment("PDG id hypothesis to be used during the fit.")
      };
      fhicl::Atom<bool> useRMS {
	Name("useRMSError"),
	Comment("Flag to replace the default hit error SigmaPeakTime() with RMS().")
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
    double pval;
    trkf::TrackKalmanFitter* kalmanFitter;
    trkf::Propagator* prop;
  };
}

trkf::KalmanFilterFinalTrackFitter::KalmanFilterFinalTrackFitter(trkf::KalmanFilterFinalTrackFitter::Parameters const & p)
  : p_(p)
{
  pval = p_().options().pval();

  prop = new trkf::PropYZPlane(0., false);
  kalmanFitter = new trkf::TrackKalmanFitter(prop,p_().options().useRMS());

  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
}

trkf::KalmanFilterFinalTrackFitter::~KalmanFilterFinalTrackFitter() {
  delete prop;
  delete kalmanFitter;
}

void trkf::KalmanFilterFinalTrackFitter::produce(art::Event & e)
{

  auto outputTracks = std::make_unique<std::vector<recob::Track> >();
  auto outputHits   = std::make_unique<art::Assns<recob::Track, recob::Hit> >();

  auto const tid = getProductID<std::vector<recob::Track> >(e);
  auto const tidgetter = e.productGetter(tid);

  art::InputTag TrackInputTag(p_().inputs().inputTracksLabel());
  art::ValidHandle<std::vector<recob::Track> > inputTracks = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag);
  std::unique_ptr<art::FindManyP<recob::Hit> > trackHits(new art::FindManyP<recob::Hit>(inputTracks, e, TrackInputTag));

  art::InputTag CaloInputTag(p_().inputs().inputCaloLabel());
  std::unique_ptr<art::FindManyP<anab::Calorimetry> > trackCalo;
  if (p_().options().pFromCalo()) {
    trackCalo = std::unique_ptr<art::FindManyP<anab::Calorimetry> >(new art::FindManyP<anab::Calorimetry>(inputTracks, e, CaloInputTag));
  } else if (p_().options().pFromMC()) {
    //FIXME, eventually remove this (ok only for single particle MC)
    art::InputTag SimTrackInputTag(p_().inputs().inputMCLabel());
    art::ValidHandle<std::vector<sim::MCTrack>> simTracks = e.getValidHandle<std::vector<sim::MCTrack>>(SimTrackInputTag);
    for (unsigned int iMC = 0; iMC < simTracks->size(); ++iMC) {
      const sim::MCTrack& mctrack = simTracks->at(iMC);
      //fiducial cuts on MC tracks
      if (mctrack.Start().Position().X()< 80 || mctrack.Start().Position().X()>176) continue;
      if (mctrack.Start().Position().Y()<-30 || mctrack.Start().Position().Y()> 30) continue;
      if (mctrack.Start().Position().Z()< 10 || mctrack.Start().Position().Z()>100) continue;
      if (mctrack.Start().Momentum().P()<500) continue;
      pval = mctrack.Start().Momentum().P()*0.001;
      break;
    }
    //std::cout << "mc momentum value = " << pval << " GeV" << std::endl;
  }

  for (unsigned int iTrack = 0; iTrack < inputTracks->size(); ++iTrack) {

    const recob::Track& track = inputTracks->at(iTrack);
    if (p_().options().pFromCalo()) {
      const std::vector<art::Ptr<anab::Calorimetry> >& calo = trackCalo->at(iTrack);
      double sumenergy = 0.;
      int nviews = 0.;
      for (auto caloit :  calo) {
	if (caloit->KineticEnergy()>0.) {
	  sumenergy+=caloit->KineticEnergy();
	  nviews+=1;
	}
      }
      if (nviews==0 || sumenergy==0.) pval = 1.;//protect against problematic cases
      else pval = sumenergy/(nviews*1000.);
    }

    recob::Track outTrack;
    art::PtrVector<recob::Hit> outHits;
    bool fitok = kalmanFitter->fitTrack(track, trackHits->at(iTrack), pval, p_().options().pdgId(), outTrack, outHits);

    if (!fitok) continue;

    outputTracks->emplace_back(std::move(outTrack));
    art::Ptr<recob::Track> aptr(tid, outputTracks->size()-1, tidgetter);
    for (auto const& trhit: outHits) {
      outputHits->addSingle(aptr, trhit);
    }
  }

  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));

}

DEFINE_ART_MODULE(trkf::KalmanFilterFinalTrackFitter)
