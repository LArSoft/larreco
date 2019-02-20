/// \class KalmanFilterFinalTrackFitter
///
/// \brief Producer for fitting Trajectories and TrackTrajectories using TrackKalmanFitter.
///
/// \author G. Cerati
///

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/TrackKalmanFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/TrackFinder/TrackMaker.h"

#include <memory>

namespace trkf {

  class KalmanFilterFinalTrackFitter : public art::EDProducer {
  public:
    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputPFParticleLabel{
        Name("inputPFParticleLabel"),
        Comment("Label of recob::PFParticle Collection to be fit")};
      fhicl::Atom<art::InputTag> inputTracksLabel{
        Name("inputTracksLabel"),
        Comment("Label of recob::Track Collection to be fit")};
      fhicl::Atom<art::InputTag> inputShowersLabel{
        Name("inputShowersLabel"),
        Comment("Label of recob::Shower Collection (associated to PFParticles) to be fit")};
      fhicl::Atom<art::InputTag> inputCaloLabel{
        Name("inputCaloLabel"),
        Comment("Label of anab::Calorimetry Collection, matching inputTracksLabel, to be used for "
                "initial momentum estimate. Used only if momFromCalo is set to true.")};
      fhicl::Atom<art::InputTag> inputMCLabel{
        Name("inputMCLabel"),
        Comment("Label of sim::MCTrack Collection to be used for initial momentum estimate. Used "
                "only if momFromMC is set to true.")};
      fhicl::Atom<art::InputTag> inputPidLabel{
        Name("inputPidLabel"),
        Comment("Label of anab::ParticleID Collection, matching inputTracksLabel, to be used for "
                "particle Id.")};
    };

    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> trackFromPF{Name("trackFromPF"),
                                    Comment("If true extract tracks from inputPFParticleLabel "
                                            "collection, if false from inputTracksLabel.")};
      fhicl::Atom<bool> showerFromPF{
        Name("showerFromPF"),
        Comment("If true extract showers from inputPFParticleLabel collection.")};
      fhicl::Atom<bool> pFromMSChi2{
        Name("momFromMSChi2"),
        Comment("Flag used to get initial momentum estimate from "
                "trkf::TrackMomentumCalculator::GetMomentumMultiScatterChi2().")};
      fhicl::Atom<bool> pFromLength{Name("momFromLength"),
                                    Comment("Flag used to get initial momentum estimate from "
                                            "trkf::TrackMomentumCalculator::GetTrackMomentum().")};
      fhicl::Atom<bool> pFromCalo{
        Name("momFromCalo"),
        Comment("Flag used to get initial momentum estimate from inputCaloLabel collection.")};
      fhicl::Atom<bool> pFromMC{
        Name("momFromMC"),
        Comment("Flag used to get initial momentum estimate from inputMCLabel collection.")};
      fhicl::Atom<double> pval{
        Name("momentumInGeV"),
        Comment("Fixed momentum estimate value, to be used when momFromCalo, momFromMSChi2, "
                "momFromLength and momFromMC are all false, or if the estimate is not available.")};
      fhicl::Atom<bool> idFromPF{Name("idFromPF"),
                                 Comment("Flag used to get particle ID estimate from corresponding "
                                         "PFParticle. Needs trackFromPF=true.")};
      fhicl::Atom<bool> idFromCollection{
        Name("idFromCollection"),
        Comment("Flag used to get particle ID estimate from inputPidLabel collection.")};
      fhicl::Atom<int> pdgId{
        Name("pdgId"),
        Comment("Default particle id hypothesis in case no valid id is provided either via "
                "PFParticle or in the ParticleId collection.")};
      fhicl::Atom<bool> dirFromVtxPF{
        Name("dirFromVtxPF"),
        Comment("Assume track direction from Vertex in PFParticle. Needs trackFromPF=true.")};
      fhicl::Atom<bool> dirFromMC{Name("dirFromMC"), Comment("Assume track direction from MC.")};
      fhicl::Atom<bool> dirFromVec{Name("dirFromVec"),
                                   Comment("Assume track direction from as the one giving positive "
                                           "dot product with vector specified by dirVec.")};
      fhicl::Sequence<float, 3u> dirVec{Name("dirVec"),
                                        Comment("Fhicl sequence defining the vector used when "
                                                "dirFromVec=true. It must have 3 elements.")};
      fhicl::Atom<bool> alwaysInvertDir{
        Name("alwaysInvertDir"),
        Comment("If true, fit all tracks from end to vertex assuming inverted direction.")};
      fhicl::Atom<bool> produceTrackFitHitInfo{
        Name("produceTrackFitHitInfo"),
        Comment("Option to produce (or not) the detailed TrackFitHitInfo.")};
      fhicl::Atom<bool> produceSpacePoints{
        Name("produceSpacePoints"),
        Comment("Option to produce (or not) the associated SpacePoints.")};
      fhicl::Atom<bool> keepInputTrajectoryPoints{
        Name("keepInputTrajectoryPoints"),
        Comment("Option to keep positions and directions from input trajectory/track. The fit will "
                "provide only covariance matrices, chi2, ndof, particle Id and absolute momentum. "
                "It may also modify the trajectory point flags. In order to avoid inconsistencies, "
                "it has to be used with the following fitter options all set to false: "
                "sortHitsByPlane, sortOutputHitsMinLength, skipNegProp.")};
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<KalmanFilterFinalTrackFitter::Inputs> inputs{
        Name("inputs"),
      };
      fhicl::Table<KalmanFilterFinalTrackFitter::Options> options{Name("options")};
      fhicl::Table<TrackStatePropagator::Config> propagator{Name("propagator")};
      fhicl::Table<TrackKalmanFitter::Config> fitter{Name("fitter")};
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit KalmanFilterFinalTrackFitter(Parameters const& p);

    // Plugins should not be copied or assigned.
    KalmanFilterFinalTrackFitter(KalmanFilterFinalTrackFitter const&) = delete;
    KalmanFilterFinalTrackFitter(KalmanFilterFinalTrackFitter&&) = delete;
    KalmanFilterFinalTrackFitter& operator=(KalmanFilterFinalTrackFitter const&) = delete;
    KalmanFilterFinalTrackFitter& operator=(KalmanFilterFinalTrackFitter&&) = delete;

  private:
    void produce(art::Event& e) override;

    Parameters p_;
    TrackStatePropagator prop;
    trkf::TrackKalmanFitter kalmanFitter;
    mutable trkf::TrackMomentumCalculator tmc{};
    bool inputFromPF;

    art::InputTag pfParticleInputTag;
    art::InputTag trackInputTag;
    art::InputTag showerInputTag;
    art::InputTag caloInputTag;
    art::InputTag pidInputTag;
    art::InputTag simTrackInputTag;

    std::unique_ptr<art::FindManyP<anab::Calorimetry>> trackCalo;
    std::unique_ptr<art::FindManyP<anab::ParticleID>> trackId;
    std::unique_ptr<art::FindManyP<recob::Track>> assocTracks;
    std::unique_ptr<art::FindManyP<recob::Shower>> assocShowers;
    std::unique_ptr<art::FindManyP<recob::Vertex>> assocVertices;

    double setMomValue(art::Ptr<recob::Track> ptrack,
                       const std::unique_ptr<art::FindManyP<anab::Calorimetry>>& trackCalo,
                       const double pMC,
                       const int pId) const;
    int setPId(const unsigned int iTrack,
               const std::unique_ptr<art::FindManyP<anab::ParticleID>>& trackId,
               const int pfPid = 0) const;
    bool setDirFlip(const recob::Track& track,
                    TVector3& mcdir,
                    const std::vector<art::Ptr<recob::Vertex>>* vertices = 0) const;

    void restoreInputPoints(const recob::Trajectory& track,
                            const std::vector<art::Ptr<recob::Hit>>& inHits,
                            recob::Track& outTrack,
                            std::vector<art::Ptr<recob::Hit>>& outHits) const;
  };
}

trkf::KalmanFilterFinalTrackFitter::KalmanFilterFinalTrackFitter(
  trkf::KalmanFilterFinalTrackFitter::Parameters const& p)
  : EDProducer{p}
  , p_(p)
  , prop{p_().propagator}
  , kalmanFitter{&prop, p_().fitter}
  , inputFromPF{p_().options().trackFromPF() || p_().options().showerFromPF()}
{

  if (inputFromPF) {
    pfParticleInputTag = art::InputTag(p_().inputs().inputPFParticleLabel());
    if (p_().options().showerFromPF())
      showerInputTag = art::InputTag(p_().inputs().inputShowersLabel());
  }
  else {
    trackInputTag = art::InputTag(p_().inputs().inputTracksLabel());
    if (p_().options().idFromCollection())
      pidInputTag = art::InputTag(p_().inputs().inputPidLabel());
  }
  if (p_().options().pFromCalo()) caloInputTag = art::InputTag(p_().inputs().inputCaloLabel());
  if (p_().options().pFromMC() || p_().options().dirFromMC())
    simTrackInputTag = art::InputTag(p_().inputs().inputMCLabel());

  produces<std::vector<recob::Track>>();
  produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  produces<art::Assns<recob::Track, recob::Hit>>();
  if (inputFromPF) { produces<art::Assns<recob::PFParticle, recob::Track>>(); }
  if (p_().options().produceTrackFitHitInfo()) {
    produces<std::vector<std::vector<recob::TrackFitHitInfo>>>();
  }
  if (p_().options().produceSpacePoints()) {
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  }

  //throw expections to avoid possible silent failures due to incompatible configuration options
  if (p_().options().trackFromPF() == 0 && p_().options().idFromPF())
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: cannot use idFromPF=true with trackFromPF=false."
      << "\n";
  if (p_().options().trackFromPF() == 0 && p_().options().dirFromVtxPF())
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: cannot use dirFromVtxPF=true with "
         "trackFromPF=false."
      << "\n";

  unsigned int nIds = 0;
  if (p_().options().idFromPF()) nIds++;
  if (p_().options().idFromCollection()) nIds++;
  if (nIds > 1) {
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among "
         "idFromPF and idFromCollection."
      << "\n";
  }

  unsigned int nDirs = 0;
  if (p_().options().dirFromVtxPF()) nDirs++;
  if (p_().options().dirFromMC()) nDirs++;
  if (p_().options().dirFromVec()) nDirs++;
  if (p_().options().alwaysInvertDir()) nDirs++;
  if (nDirs > 1) {
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among "
         "dirFromVtxPF, dirFromMC, dirFromVec, and alwaysInvertDir."
      << "\n";
  }

  unsigned int nPFroms = 0;
  if (p_().options().pFromCalo()) nPFroms++;
  if (p_().options().pFromMSChi2()) nPFroms++;
  if (p_().options().pFromLength()) nPFroms++;
  if (p_().options().pFromMC()) nPFroms++;
  if (nPFroms > 1) {
    throw cet::exception("KalmanFilterFinalTrackFitter")
      << "Incompatible configuration parameters: only at most one can be set to true among "
         "pFromCalo, pFromMSChi2, pFromLength, and pFromMC."
      << "\n";
  }

  if (p_().options().keepInputTrajectoryPoints()) {
    if (p_().fitter().sortHitsByPlane() || p_().fitter().sortOutputHitsMinLength() ||
        p_().fitter().skipNegProp()) {
      throw cet::exception("KalmanFilterTrajectoryFitter")
        << "Incompatible configuration parameters: keepInputTrajectoryPoints needs the following "
           "fitter options all set to false: sortHitsByPlane, sortOutputHitsMinLength, skipNegProp."
        << "\n";
    }
  }

  if (p_().options().showerFromPF()) {
    if (nPFroms > 0 || nIds > 0 || nDirs > 0) {
      throw cet::exception("KalmanFilterTrajectoryFitter")
        << "Incompatible configuration parameters: showerFromPF currently does not support "
           "optional momentum values, particle hypotheses and directions."
        << "\n";
    }
  }
}

void
trkf::KalmanFilterFinalTrackFitter::produce(art::Event& e)
{
  auto outputTracks = std::make_unique<std::vector<recob::Track>>();
  auto outputHitsMeta =
    std::make_unique<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  auto outputHits = std::make_unique<art::Assns<recob::Track, recob::Hit>>();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo>>>();

  auto const tid = e.getProductID<std::vector<recob::Track>>();
  auto const tidgetter = e.productGetter(tid);

  auto outputSpacePoints = std::make_unique<std::vector<recob::SpacePoint>>();
  auto outputHitSpacePointAssn = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
  auto const spid = e.getProductID<std::vector<recob::SpacePoint>>();
  auto const spidgetter = e.productGetter(spid);

  //FIXME, eventually remove this (ok only for single particle MC)
  double pMC = -1.;
  TVector3 mcdir;
  if (p_().options().pFromMC() || p_().options().dirFromMC()) {
    art::ValidHandle<std::vector<sim::MCTrack>> simTracks =
      e.getValidHandle<std::vector<sim::MCTrack>>(simTrackInputTag);
    for (unsigned int iMC = 0; iMC < simTracks->size(); ++iMC) {
      const sim::MCTrack& mctrack = simTracks->at(iMC);
      //fiducial cuts on MC tracks
      if (mctrack.PdgCode() != 13) continue;
      if (mctrack.Process() != "primary") continue;
      pMC = mctrack.Start().Momentum().P() * 0.001;
      mcdir = TVector3(mctrack.Start().Momentum().X() * 0.001 / pMC,
                       mctrack.Start().Momentum().Y() * 0.001 / pMC,
                       mctrack.Start().Momentum().Z() * 0.001 / pMC);
      break;
    }
  }

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  if (inputFromPF) {

    auto outputPFAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Track>>();

    auto inputPFParticle = e.getValidHandle<std::vector<recob::PFParticle>>(pfParticleInputTag);
    if (p_().options().trackFromPF())
      assocTracks =
        std::make_unique<art::FindManyP<recob::Track>>(inputPFParticle, e, pfParticleInputTag);
    if (p_().options().showerFromPF())
      assocShowers =
        std::make_unique<art::FindManyP<recob::Shower>>(inputPFParticle, e, showerInputTag);
    assocVertices =
      std::make_unique<art::FindManyP<recob::Vertex>>(inputPFParticle, e, pfParticleInputTag);

    for (unsigned int iPF = 0; iPF < inputPFParticle->size(); ++iPF) {

      if (p_().options().trackFromPF()) {
        const std::vector<art::Ptr<recob::Track>>& tracks = assocTracks->at(iPF);
        auto const& tkHitsAssn =
          *e.getValidHandle<art::Assns<recob::Track, recob::Hit>>(pfParticleInputTag);
        const std::vector<art::Ptr<recob::Vertex>>& vertices = assocVertices->at(iPF);

        if (p_().options().pFromCalo()) {
          trackCalo = std::make_unique<art::FindManyP<anab::Calorimetry>>(tracks, e, caloInputTag);
        }

        for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack) {

          const recob::Track& track = *tracks[iTrack];
          art::Ptr<recob::Track> ptrack = tracks[iTrack];
          const int pId = setPId(iTrack, trackId, inputPFParticle->at(iPF).PdgCode());
          const double mom = setMomValue(ptrack, trackCalo, pMC, pId);
          const bool flipDir = setDirFlip(track, mcdir, &vertices);

          //this is not computationally optimal, but at least preserves the order unlike FindManyP
          std::vector<art::Ptr<recob::Hit>> inHits;
          for (auto it = tkHitsAssn.begin(); it != tkHitsAssn.end(); ++it) {
            if (it->first == ptrack)
              inHits.push_back(it->second);
            else if (inHits.size() > 0)
              break;
          }

          recob::Track outTrack;
          std::vector<art::Ptr<recob::Hit>> outHits;
          trkmkr::OptionalOutputs optionals;
          if (p_().options().produceTrackFitHitInfo()) optionals.initTrackFitInfos();
          bool fitok = kalmanFitter.fitTrack(detProp,
                                             track.Trajectory(),
                                             track.ID(),
                                             track.VertexCovarianceLocal5D(),
                                             track.EndCovarianceLocal5D(),
                                             inHits,
                                             mom,
                                             pId,
                                             flipDir,
                                             outTrack,
                                             outHits,
                                             optionals);
          if (!fitok) continue;

          if (p_().options().keepInputTrajectoryPoints()) {
            restoreInputPoints(track.Trajectory().Trajectory(), inHits, outTrack, outHits);
          }

          outputTracks->emplace_back(std::move(outTrack));
          art::Ptr<recob::Track> aptr(tid, outputTracks->size() - 1, tidgetter);
          unsigned int ip = 0;
          for (auto const& trhit : outHits) {
            //the fitter produces collections with 1-1 match between hits and point
            recob::TrackHitMeta metadata(ip, -1);
            outputHitsMeta->addSingle(aptr, trhit, metadata);
            outputHits->addSingle(aptr, trhit);
            ip++;
          }
          outputPFAssn->addSingle(art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr);
          outputHitInfo->emplace_back(optionals.trackFitHitInfos());
        }
      }

      if (p_().options().showerFromPF()) {
        art::Ptr<recob::PFParticle> pPF(inputPFParticle, iPF);
        const std::vector<art::Ptr<recob::Shower>>& showers = assocShowers->at(iPF);
        if (showers.size() == 0) continue;
        auto const& pfClustersAssn =
          *e.getValidHandle<art::Assns<recob::PFParticle, recob::Cluster>>(showerInputTag);
        auto const& clHitsAssn =
          *e.getValidHandle<art::Assns<recob::Cluster, recob::Hit>>(showerInputTag);
        std::vector<art::Ptr<recob::Hit>> inHits;
        for (auto itpf = pfClustersAssn.begin(); itpf != pfClustersAssn.end(); ++itpf) {
          if (itpf->first == pPF) {
            art::Ptr<recob::Cluster> clust = itpf->second;
            for (auto it = clHitsAssn.begin(); it != clHitsAssn.end(); ++it) {
              if (it->first == clust) inHits.push_back(it->second);
            }
          }
          else if (inHits.size() > 0)
            break;
        }
        for (unsigned int iShower = 0; iShower < showers.size(); ++iShower) {
          const recob::Shower& shower = *showers[iShower];
          recob::Track outTrack;
          std::vector<art::Ptr<recob::Hit>> outHits;
          trkmkr::OptionalOutputs optionals;
          if (p_().options().produceTrackFitHitInfo()) optionals.initTrackFitInfos();
          Point_t pos(shower.ShowerStart().X(), shower.ShowerStart().Y(), shower.ShowerStart().Z());
          Vector_t dir(shower.Direction().X(), shower.Direction().Y(), shower.Direction().Z());
          auto cov = SMatrixSym55();
          auto pid = p_().options().pdgId();
          auto mom = p_().options().pval();
          bool fitok = kalmanFitter.fitTrack(detProp,
                                             pos,
                                             dir,
                                             cov,
                                             inHits,
                                             std::vector<recob::TrajectoryPointFlags>(),
                                             shower.ID(),
                                             mom,
                                             pid,
                                             outTrack,
                                             outHits,
                                             optionals);
          if (!fitok) continue;

          outputTracks->emplace_back(std::move(outTrack));
          art::Ptr<recob::Track> aptr(tid, outputTracks->size() - 1, tidgetter);
          unsigned int ip = 0;
          for (auto const& trhit : outHits) {
            // the fitter produces collections with 1-1 match between hits and point
            recob::TrackHitMeta metadata(ip, -1);
            outputHitsMeta->addSingle(aptr, trhit, metadata);
            outputHits->addSingle(aptr, trhit);
            if (p_().options().produceSpacePoints() && outputTracks->back().HasValidPoint(ip)) {
              auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
              double fXYZ[3] = {tp.X(), tp.Y(), tp.Z()};
              double fErrXYZ[6] = {0};
              recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
              outputSpacePoints->emplace_back(std::move(sp));
              art::Ptr<recob::SpacePoint> apsp(spid, outputSpacePoints->size() - 1, spidgetter);
              outputHitSpacePointAssn->addSingle(trhit, apsp);
            }
            ip++;
          }
          outputPFAssn->addSingle(art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr);
          outputHitInfo->emplace_back(optionals.trackFitHitInfos());
        }
      }
    }
    e.put(std::move(outputTracks));
    e.put(std::move(outputHitsMeta));
    e.put(std::move(outputHits));
    e.put(std::move(outputPFAssn));
    if (p_().options().produceTrackFitHitInfo()) { e.put(std::move(outputHitInfo)); }
    if (p_().options().produceSpacePoints()) {
      e.put(std::move(outputSpacePoints));
      e.put(std::move(outputHitSpacePointAssn));
    }
  }
  else {

    art::ValidHandle<std::vector<recob::Track>> inputTracks =
      e.getValidHandle<std::vector<recob::Track>>(trackInputTag);
    auto const& tkHitsAssn = *e.getValidHandle<art::Assns<recob::Track, recob::Hit>>(trackInputTag);

    if (p_().options().pFromCalo()) {
      trackCalo = std::make_unique<art::FindManyP<anab::Calorimetry>>(inputTracks, e, caloInputTag);
    }

    if (p_().options().idFromCollection()) {
      trackId = std::make_unique<art::FindManyP<anab::ParticleID>>(inputTracks, e, pidInputTag);
    }

    for (unsigned int iTrack = 0; iTrack < inputTracks->size(); ++iTrack) {

      const recob::Track& track = inputTracks->at(iTrack);
      art::Ptr<recob::Track> ptrack(inputTracks, iTrack);
      const int pId = setPId(iTrack, trackId);
      const double mom = setMomValue(ptrack, trackCalo, pMC, pId);
      const bool flipDir = setDirFlip(track, mcdir);

      //this is not computationally optimal, but at least preserves the order unlike FindManyP
      std::vector<art::Ptr<recob::Hit>> inHits;
      for (auto it = tkHitsAssn.begin(); it != tkHitsAssn.end(); ++it) {
        if (it->first == ptrack)
          inHits.push_back(it->second);
        else if (inHits.size() > 0)
          break;
      }

      recob::Track outTrack;
      std::vector<art::Ptr<recob::Hit>> outHits;
      trkmkr::OptionalOutputs optionals;
      if (p_().options().produceTrackFitHitInfo()) optionals.initTrackFitInfos();
      bool fitok = kalmanFitter.fitTrack(detProp,
                                         track.Trajectory(),
                                         track.ID(),
                                         track.VertexCovarianceLocal5D(),
                                         track.EndCovarianceLocal5D(),
                                         inHits,
                                         mom,
                                         pId,
                                         flipDir,
                                         outTrack,
                                         outHits,
                                         optionals);
      if (!fitok) continue;

      if (p_().options().keepInputTrajectoryPoints()) {
        restoreInputPoints(track.Trajectory().Trajectory(), inHits, outTrack, outHits);
      }

      outputTracks->emplace_back(std::move(outTrack));
      art::Ptr<recob::Track> aptr(tid, outputTracks->size() - 1, tidgetter);
      unsigned int ip = 0;
      for (auto const& trhit : outHits) {
        //the fitter produces collections with 1-1 match between hits and point
        recob::TrackHitMeta metadata(ip, -1);
        outputHitsMeta->addSingle(aptr, trhit, metadata);
        outputHits->addSingle(aptr, trhit);
        if (p_().options().produceSpacePoints() && outputTracks->back().HasValidPoint(ip)) {
          auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
          double fXYZ[3] = {tp.X(), tp.Y(), tp.Z()};
          double fErrXYZ[6] = {0};
          recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
          outputSpacePoints->emplace_back(std::move(sp));
          art::Ptr<recob::SpacePoint> apsp(spid, outputSpacePoints->size() - 1, spidgetter);
          outputHitSpacePointAssn->addSingle(trhit, apsp);
        }
        ip++;
      }
      outputHitInfo->emplace_back(optionals.trackFitHitInfos());
    }
    e.put(std::move(outputTracks));
    e.put(std::move(outputHitsMeta));
    e.put(std::move(outputHits));
    if (p_().options().produceTrackFitHitInfo()) { e.put(std::move(outputHitInfo)); }
    if (p_().options().produceSpacePoints()) {
      e.put(std::move(outputSpacePoints));
      e.put(std::move(outputHitSpacePointAssn));
    }
  }
}

void
trkf::KalmanFilterFinalTrackFitter::restoreInputPoints(
  const recob::Trajectory& track,
  const std::vector<art::Ptr<recob::Hit>>& inHits,
  recob::Track& outTrack,
  std::vector<art::Ptr<recob::Hit>>& outHits) const
{
  const auto np = outTrack.NumberTrajectoryPoints();
  std::vector<Point_t> positions(np);
  std::vector<Vector_t> momenta(np);
  std::vector<recob::TrajectoryPointFlags> outFlags(np);
  //
  for (unsigned int p = 0; p < np; ++p) {
    auto flag = outTrack.FlagsAtPoint(p);
    auto mom = outTrack.VertexMomentum();
    auto op = flag.fromHit();
    positions[op] = track.LocationAtPoint(op);
    momenta[op] = mom * track.DirectionAtPoint(op);
    auto mask = flag.mask();
    if (mask.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
      mask.unset(recob::TrajectoryPointFlagTraits::NoPoint);
    outFlags[op] = recob::TrajectoryPointFlags(op, mask);
  }
  auto covs = outTrack.Covariances();
  outTrack = recob::Track(
    recob::TrackTrajectory(std::move(positions), std::move(momenta), std::move(outFlags), true),
    outTrack.ParticleId(),
    outTrack.Chi2(),
    outTrack.Ndof(),
    std::move(covs.first),
    std::move(covs.second),
    outTrack.ID());
  //
  outHits.clear();
  for (auto h : inHits)
    outHits.push_back(h);
}

double
trkf::KalmanFilterFinalTrackFitter::setMomValue(
  art::Ptr<recob::Track> ptrack,
  const std::unique_ptr<art::FindManyP<anab::Calorimetry>>& trackCalo,
  const double pMC,
  const int pId) const
{
  double result = p_().options().pval();
  if (p_().options().pFromMSChi2()) { result = tmc.GetMomentumMultiScatterChi2(ptrack); }
  else if (p_().options().pFromLength()) {
    result = tmc.GetTrackMomentum(ptrack->Length(), pId);
  }
  else if (p_().options().pFromCalo()) {
    //take average energy from available views
    const std::vector<art::Ptr<anab::Calorimetry>>& calo = trackCalo->at(ptrack.key());
    double sumenergy = 0.;
    int nviews = 0.;
    for (auto caloit : calo) {
      if (caloit->KineticEnergy() > 0.) {
        sumenergy += caloit->KineticEnergy();
        nviews += 1;
      }
    }
    if (nviews != 0 && sumenergy != 0.) {
      //protect against problematic cases
      result = sumenergy / (nviews * 1000.);
    }
  }
  else if (p_().options().pFromMC() && pMC > 0.) {
    result = pMC;
  }
  return result;
}

int
trkf::KalmanFilterFinalTrackFitter::setPId(
  const unsigned int iTrack,
  const std::unique_ptr<art::FindManyP<anab::ParticleID>>& trackId,
  const int pfPid) const
{
  /*
  int result = p_().options().pdgId();
  if (p_().options().trackFromPF() && p_().options().idFromPF()) { result = pfPid; }
  else if (p_().options().idFromCollection()) {
    //take the pdgId corresponding to the minimum chi2 (should we give preference to the majority? fixme)
    double minChi2 = -1.;
    for (auto idit : trackId->at(iTrack)) {
      if (idit->MinChi2() > 0. && (minChi2 < 0. || idit->MinChi2() < minChi2)) {
        result = idit->Pdg();
        minChi2 = idit->MinChi2();
      }
    }
  }
  */
  return -1;
}

bool
trkf::KalmanFilterFinalTrackFitter::setDirFlip(
  const recob::Track& track,
  TVector3& mcdir,
  const std::vector<art::Ptr<recob::Vertex>>* vertices) const
{
  bool result = false;
  if (p_().options().alwaysInvertDir()) { return true; }
  else if (p_().options().dirFromMC()) {
    auto tdir = track.VertexDirection();
    if ((mcdir.X() * tdir.X() + mcdir.Y() * tdir.Y() + mcdir.Z() * tdir.Z()) < 0.) result = true;
  }
  else if (p_().options().dirFromVec()) {
    std::array<float, 3> dir = p_().options().dirVec();
    auto tdir = track.VertexDirection();
    if ((dir[0] * tdir.X() + dir[1] * tdir.Y() + dir[2] * tdir.Z()) < 0.) result = true;
  }
  else if (p_().options().trackFromPF() && p_().options().dirFromVtxPF() && vertices->size() > 0) {
    //if track end is closer to first vertex then track vertex, flip direction
    double xyz[3];
    (*vertices)[0]->XYZ(xyz);
    auto& tv = track.Trajectory().Vertex();
    auto& te = track.Trajectory().End();
    if (((xyz[0] - te.X()) * (xyz[0] - te.X()) + (xyz[1] - te.Y()) * (xyz[1] - te.Y()) +
         (xyz[2] - te.Z()) * (xyz[2] - te.Z())) >
        ((xyz[0] - tv.X()) * (xyz[0] - tv.X()) + (xyz[1] - tv.Y()) * (xyz[1] - tv.Y()) +
         (xyz[2] - tv.Z()) * (xyz[2] - tv.Z())))
      result = true;
  }
  return result;
}

DEFINE_ART_MODULE(trkf::KalmanFilterFinalTrackFitter)
