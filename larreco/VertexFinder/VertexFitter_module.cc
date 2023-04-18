#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

#include <memory>

namespace trkf {

  /**
   * @file  larreco/VertexFinder/VertexFitter_module.cc
   * @class trkf::VertexFitter
   *
   * @brief Module for fitting a vertex using the Geometric3DVertexFitter.
   *
   * It selects primary PFParticles, and then collects all tracks associated to its daughters;
   * if at least 2 tracks are found, they are passed to the vertex fitter.
   *
   * Inputs are: a PFParticle collection, and the associated tracks.
   *
   * Outputs are: vector of recob::Vertex, Assns of (neutrino) recob::PFParticle to recob::Vertex,
   * Assns of recob::Vertex and recob::Track with recob::VertexAssnMeta.
   *
   * For configuration options see Geometric3DVertexFitter#Parameters
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  class VertexFitter : public art::EDProducer {
  public:
    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputPFParticleLabel{
        Name("inputPFParticleLabel"),
        Comment("Label of recob::PFParticle Collection to be fit")};
      fhicl::Atom<art::InputTag> inputTracksLabel{
        Name("inputTracksLabel"),
        Comment("Label of recob::Track Collection associated to PFParticles")};
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<VertexFitter::Inputs> inputs{
        Name("inputs"),
      };
      fhicl::Table<Geometric3DVertexFitter::Config> geom3dvtxfit{Name("geom3dvtxfit")};
      fhicl::Table<TrackStatePropagator::Config> propagator{Name("propagator")};
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit VertexFitter(Parameters const& p);

    // Plugins should not be copied or assigned.
    VertexFitter(VertexFitter const&) = delete;
    VertexFitter(VertexFitter&&) = delete;
    VertexFitter& operator=(VertexFitter const&) = delete;
    VertexFitter& operator=(VertexFitter&&) = delete;

  private:
    void produce(art::Event& e) override;

    art::InputTag pfParticleInputTag;
    art::InputTag trackInputTag;
    Geometric3DVertexFitter fitter;
  };
}

trkf::VertexFitter::VertexFitter(Parameters const& p)
  : EDProducer{p}
  , pfParticleInputTag(p().inputs().inputPFParticleLabel())
  , trackInputTag(p().inputs().inputTracksLabel())
  , fitter(p().geom3dvtxfit, p().propagator)
{
  produces<std::vector<recob::Vertex>>();
  produces<art::Assns<recob::PFParticle, recob::Vertex>>();
  produces<art::Assns<recob::Vertex, recob::Track, recob::VertexAssnMeta>>();
}

void trkf::VertexFitter::produce(art::Event& e)
{
  using namespace std;

  auto outputVertices = make_unique<vector<recob::Vertex>>();
  auto outputPFVxAssn = make_unique<art::Assns<recob::PFParticle, recob::Vertex>>();
  auto outputVxTkMtAssn =
    make_unique<art::Assns<recob::Vertex, recob::Track, recob::VertexAssnMeta>>();

  const auto& inputPFParticle = e.getValidHandle<vector<recob::PFParticle>>(pfParticleInputTag);
  art::FindManyP<recob::Track> const assocTracks{inputPFParticle, e, trackInputTag};

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  // PtrMakers for Assns
  art::PtrMaker<recob::Vertex> vtxPtrMaker(e);

  for (size_t iPF = 0; iPF < inputPFParticle->size(); ++iPF) {

    art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPF);
    if (pfp->IsPrimary() == false || pfp->NumDaughters() < 2) continue;
    vector<art::Ptr<recob::Track>> tracks;
    auto& pfd = pfp->Daughters();
    for (auto ipfd : pfd) {
      // Daugthers returns the id as in pfp->Self() not the key
      // so need to find the key for the pfp with Self==ipfd
      for (size_t jPF = 0; jPF < inputPFParticle->size(); ++jPF) {
        art::Ptr<recob::PFParticle> pfpd(inputPFParticle, jPF);
        if (pfpd->Self() != ipfd) continue;
        vector<art::Ptr<recob::Track>> pftracks = assocTracks.at(jPF);
        for (auto t : pftracks) {
          tracks.push_back(t);
        }
        break;
      }
    }
    if (tracks.size() < 2) continue;

    VertexWrapper vtx = fitter.fitTracks(detProp, tracks);
    if (vtx.isValid() == false) continue;
    vtx.setVertexId(outputVertices->size());

    auto meta = fitter.computeMeta(detProp, vtx, tracks);

    // Fill the output collections

    outputVertices->emplace_back(vtx.vertex());
    const art::Ptr<recob::Vertex> aptr = vtxPtrMaker(outputVertices->size() - 1);
    outputPFVxAssn->addSingle(art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr);

    size_t itt = 0;
    for (auto t : tracks) {
      outputVxTkMtAssn->addSingle(aptr, t, meta[itt]);
      itt++;
    }
  }

  e.put(std::move(outputVertices));
  e.put(std::move(outputPFVxAssn));
  e.put(std::move(outputVxTkMtAssn));
}

DEFINE_ART_MODULE(trkf::VertexFitter)
