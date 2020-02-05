// Chris Backhouse - c.backhouse@ucl.ac.uk - Oct 2019

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

namespace{
  int gEvt;
  double gTrueX, gTrueY, gTrueZ;
  std::vector<double> gRecoX, gRecoY, gRecoZ;
}

namespace quad
{

class EvalVtx: public art::EDAnalyzer
{
public:
  explicit EvalVtx(const fhicl::ParameterSet& pset);

  void analyze(const art::Event& evt) override;

protected:
  std::string fTruthLabel;
  std::vector<std::string> fVertexLabels;

  TTree* fTree;
};

DEFINE_ART_MODULE(EvalVtx)

// ---------------------------------------------------------------------------
EvalVtx::EvalVtx(const fhicl::ParameterSet& pset) :
  EDAnalyzer(pset),
  fTruthLabel(pset.get<std::string>("TruthLabel")),
  fVertexLabels(pset.get<std::vector<std::string>>("VertexLabels"))
{
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("vtxs", "vtxs");

  fTree->Branch("evt", &gEvt);
  fTree->Branch("true_x", &gTrueX);
  fTree->Branch("true_y", &gTrueY);
  fTree->Branch("true_z", &gTrueZ);

  gRecoX.resize(fVertexLabels.size());
  gRecoY.resize(fVertexLabels.size());
  gRecoZ.resize(fVertexLabels.size());

  for(unsigned int i = 0; i < fVertexLabels.size(); ++i){
    const std::string& l = fVertexLabels[i];
    fTree->Branch((l+"_x").c_str(), &gRecoX[i]);
    fTree->Branch((l+"_y").c_str(), &gRecoY[i]);
    fTree->Branch((l+"_z").c_str(), &gRecoZ[i]);
  }
}

// ---------------------------------------------------------------------------
recob::Vertex GetFirstVertex(const std::string& label, const art::Event& evt)
{
  art::Handle<std::vector<recob::Vertex>> vtxs;
  evt.getByLabel(label, vtxs);

  if(vtxs->empty()) return recob::Vertex();

  return (*vtxs)[0];
}

// ---------------------------------------------------------------------------
recob::Vertex GetVtxByAssns(const std::string& label, const art::Event& evt)
{
  art::Handle<std::vector<recob::Vertex>> vtxs;
  evt.getByLabel(label, vtxs);

  if(vtxs->empty()) return recob::Vertex();

  art::Handle<std::vector<recob::PFParticle>> parts;
  evt.getByLabel(label, parts);

  art::FindManyP<recob::Vertex> fm(parts, evt, label);

  for(unsigned int i = 0; i < parts->size(); ++i){
    const int pdg = abs((*parts)[i].PdgCode());
    if(pdg == 12 || pdg == 14 || pdg == 16){
      const std::vector<art::Ptr<recob::Vertex>> vtxs = fm.at(i);
      if(vtxs.size() == 1) return *vtxs[0];

      if(vtxs.empty()) std::cout << "Warning: vertex list empty!" << std::endl;
      if(!vtxs.empty()) std::cout << "Warning: " << vtxs.size() << " vertices for daughter?" << std::endl;
    }
  }

  return recob::Vertex();
}

// ---------------------------------------------------------------------------
void EvalVtx::analyze(const art::Event& evt)
{
  art::Handle<std::vector<simb::MCTruth>> truths;
  evt.getByLabel(fTruthLabel, truths);
  if(truths->empty()) return;

  const simb::MCParticle& nu = (*truths)[0].GetNeutrino().Nu();

  gEvt = evt.event();

  gTrueX = nu.Vx();
  gTrueY = nu.Vy();
  gTrueZ = nu.Vz();

  for(unsigned int i = 0; i < fVertexLabels.size(); ++i){
    const recob::Vertex vtx = (fVertexLabels[i] == "pandora") ? GetVtxByAssns(fVertexLabels[i], evt) : GetFirstVertex(fVertexLabels[i], evt);
    const recob::Vertex::Point_t reco_vtx = vtx.position();

    gRecoX[i] = reco_vtx.x();
    gRecoY[i] = reco_vtx.y();
    gRecoZ[i] = reco_vtx.z();
  } // end for i

  fTree->Fill();
}

} // namespace
