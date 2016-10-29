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

class KalmanFilterFinalTrackFitter;

class KalmanFilterFinalTrackFitter : public art::EDProducer {
public:
  explicit KalmanFilterFinalTrackFitter(fhicl::ParameterSet const & p);
  ~KalmanFilterFinalTrackFitter();

  // Plugins should not be copied or assigned.
  KalmanFilterFinalTrackFitter(KalmanFilterFinalTrackFitter const &) = delete;
  KalmanFilterFinalTrackFitter(KalmanFilterFinalTrackFitter &&) = delete;
  KalmanFilterFinalTrackFitter & operator = (KalmanFilterFinalTrackFitter const &) = delete;
  KalmanFilterFinalTrackFitter & operator = (KalmanFilterFinalTrackFitter &&) = delete;

  void produce(art::Event & e) override;

private:
  std::string inputTracksLabel;
  bool pFromCalo, pFromMC;
  std::string inputCaloLabel;
  std::string inputMCLabel;
  double pval;
  int pdgId;
  trkf::TrackKalmanFitter* kalmanFitter;
  trkf::Propagator* prop;
};

KalmanFilterFinalTrackFitter::KalmanFilterFinalTrackFitter(fhicl::ParameterSet const & p)
  : inputTracksLabel(p.get<std::string>("inputTracksLabel"))
{
  pFromCalo = p.get<bool>("momFromCalo");
  pFromMC = p.get<bool>("momFromMC");
  if (pFromCalo) {
    inputCaloLabel = p.get<std::string>("inputCaloLabel");
  } else if (pFromMC) {
    inputMCLabel = p.get<std::string>("inputMCLabel");
  } else {
    pval = p.get<double>("momentumInGeV");
  }
  pdgId = p.get<int>("pdgId");
  
  prop = new trkf::PropYZPlane(0., false);
  kalmanFitter = new trkf::TrackKalmanFitter(prop);

  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
}

KalmanFilterFinalTrackFitter::~KalmanFilterFinalTrackFitter() {
  delete prop;
  delete kalmanFitter;
}

void KalmanFilterFinalTrackFitter::produce(art::Event & e)
{
  
  auto outputTracks = std::make_unique<std::vector<recob::Track> >();
  auto outputHits   = std::make_unique<art::Assns<recob::Track, recob::Hit> >();

  auto const tid = getProductID<std::vector<recob::Track> >(e);
  auto const tidgetter = e.productGetter(tid);
  
  art::InputTag TrackInputTag(inputTracksLabel);
  art::ValidHandle<std::vector<recob::Track> > inputTracks = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag);
  std::unique_ptr<art::FindManyP<recob::Hit> > trackHits(new art::FindManyP<recob::Hit>(inputTracks, e, TrackInputTag));

  art::InputTag CaloInputTag(inputCaloLabel);
  std::unique_ptr<art::FindManyP<anab::Calorimetry> > trackCalo;
  if (pFromCalo) {
    trackCalo = std::unique_ptr<art::FindManyP<anab::Calorimetry> >(new art::FindManyP<anab::Calorimetry>(inputTracks, e, CaloInputTag));
  } else if (pFromMC) {
    //FIXME, eventually remove this (ok only for single particle MC)
    art::InputTag SimTrackInputTag(inputMCLabel);
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
    if (pFromCalo) {
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
    bool fitok = kalmanFitter->fitTrack(track, trackHits->at(iTrack), pval, pdgId, outTrack, outHits);

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

DEFINE_ART_MODULE(KalmanFilterFinalTrackFitter)
