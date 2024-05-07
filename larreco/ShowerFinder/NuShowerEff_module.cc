////////////////////////////////////////////////////////////////////////
// Class:       NuShowerEff
// Plugin Type: analyzer (art v2_11_03)
// File:        NuShowerEff_module.cc
//
// Generated at Tue Sep 18 14:20:57 2018 by Wanwei Wu using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArsoft includes
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TTree.h"

// user's defined constant
#define MAX_SHOWERS 1000

using namespace std;

class NuShowerEff;

class NuShowerEff : public art::EDAnalyzer {
public:
  explicit NuShowerEff(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuShowerEff(NuShowerEff const&) = delete;
  NuShowerEff(NuShowerEff&&) = delete;
  NuShowerEff& operator=(NuShowerEff const&) = delete;
  NuShowerEff& operator=(NuShowerEff&&) = delete;

private:
  // Required functions.
  void analyze(art::Event const& e) override;

  // user's defined functions
  void beginJob() override;
  void endJob() override;
  void beginRun(art::Run const& run) override;
  void doEfficiencies();
  bool insideFV(double vertex[4]);
  void reset();

  // Declare member data here.

  // read from fchicl: Input Tag
  std::string fMCTruthModuleLabel;
  std::string fHitModuleLabel;
  std::string fShowerModuleLabel;
  std::string fTruthMatchDataModuleLabel;

  // read from fchicl: user's defined parameters
  int fNeutrinoPDGcode;
  int fLeptonPDGcode;
  bool fSaveMCTree;
  bool fHitShowerThroughPFParticle; // an option for getting hits associated with shower
  double fMinPurity;                // for reference
  double fMinCompleteness;          // for reference
  float fFidVolCutX;                // Fiducail Volume cut [cm]
  float fFidVolCutY;
  float fFidVolCutZ;

  // Fiducial Volume parameters
  float fFidVolXmin;
  float fFidVolXmax;
  float fFidVolYmin;
  float fFidVolYmax;
  float fFidVolZmin;
  float fFidVolZmax;

  // Event
  int Event;
  int Run;
  int SubRun;

  // MC Truth: Generator
  int MC_incoming_PDG;
  int MC_lepton_PDG;
  int MC_isCC;
  int MC_channel;
  int MC_target;
  double MC_Q2;
  double MC_W;
  double MC_vertex[4];
  double MC_incoming_P[4];
  double MC_lepton_startMomentum[4];
  double MC_lepton_endMomentum[40];
  double MC_lepton_startXYZT[4];
  double MC_lepton_endXYZT[4];
  double MC_lepton_theta;
  int MC_leptonID;
  int MC_LeptonTrack;

  double MC_incoming_E; // incoming neutrino energy
  double MC_lepton_startEnergy;

  // recob::Shower
  int n_recoShowers;
  double sh_direction_X[MAX_SHOWERS];
  double sh_direction_Y[MAX_SHOWERS];
  double sh_direction_Z[MAX_SHOWERS];
  double sh_start_X[MAX_SHOWERS];
  double sh_start_Y[MAX_SHOWERS];
  double sh_start_Z[MAX_SHOWERS];
  double sh_length[MAX_SHOWERS];
  double sh_ehit_Q[MAX_SHOWERS];
  int sh_TrackId[MAX_SHOWERS];
  int sh_hasPrimary_e[MAX_SHOWERS];
  double sh_allhit_Q[MAX_SHOWERS];
  double sh_purity[MAX_SHOWERS];
  double sh_completeness[MAX_SHOWERS];
  double esh_1_purity; // largest shower in a CC event that contains electron contributions
  double esh_1_completeness;
  double
    esh_each_purity[MAX_SHOWERS]; // each shower in a CC event that contains electron contributions
  double esh_each_completeness[MAX_SHOWERS];
  double esh_purity; // average over all electron showers in a CC event
  double esh_completeness;
  int count_primary_e_in_Event; // number of showers containing electron contribution in a CC event

  // TTree
  TTree* fEventTree;

  TH1D* h_Ev_den; // incoming neutrino energy from MC. den: denominator
  TH1D* h_Ev_num; // recon. num: numerator

  TH1D* h_Ee_den; // primary electron energy from MC
  TH1D* h_Ee_num;

  TEfficiency* h_Eff_Ev = 0;
  TEfficiency* h_Eff_Ee = 0;
};

// =====================================================================================
NuShowerEff::NuShowerEff(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
  , fMCTruthModuleLabel(p.get<std::string>("MCTruthModuleLabel"))
  , fHitModuleLabel(p.get<std::string>("HitModuleLabel"))
  , fShowerModuleLabel(p.get<std::string>("ShowerModuleLabel"))
  , fTruthMatchDataModuleLabel(p.get<std::string>("TruthMatchDataModuleLabel"))
  , fNeutrinoPDGcode(p.get<int>("NeutrinoPDGcode"))
  , fLeptonPDGcode(p.get<int>("LeptonPDGcode"))
  , fSaveMCTree(p.get<bool>("SaveMCTree"))
  , fHitShowerThroughPFParticle(p.get<bool>("HitShowerThroughPFParticle"))
  , fMinPurity(p.get<double>("MinPurity"))
  , fMinCompleteness(p.get<double>("MinCompleteness"))
  , fFidVolCutX(p.get<float>("FidVolCutX"))
  , fFidVolCutY(p.get<float>("FidVolCutY"))
  , fFidVolCutZ(p.get<float>("FidVolCutZ"))
{}

//============================================================================
void NuShowerEff::beginJob()
{
  // Get geometry: Fiducial Volume
  auto const* geo = lar::providerFrom<geo::Geometry>();
  double minx = 1e9; // [cm]
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (auto const& tpc : geo->Iterate<geo::TPCGeo>()) {
    auto const world = tpc.GetCenter();

    if (minx > world.X() - tpc.HalfWidth()) minx = world.X() - tpc.HalfWidth();
    if (maxx < world.X() + tpc.HalfWidth()) maxx = world.X() + tpc.HalfWidth();
    if (miny > world.Y() - tpc.HalfHeight()) miny = world.Y() - tpc.HalfHeight();
    if (maxy < world.Y() + tpc.HalfHeight()) maxy = world.Y() + tpc.HalfHeight();
    if (minz > world.Z() - tpc.Length() / 2.) minz = world.Z() - tpc.Length() / 2.;
    if (maxz < world.Z() + tpc.Length() / 2.) maxz = world.Z() + tpc.Length() / 2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  art::ServiceHandle<art::TFileService const> tfs;

  double E_bins[21] = {0,   0.5, 1.0, 1.5, 2.0,  2.5,  3.0,  3.5,  4,    4.5, 5.0,
                       5.5, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 17.0, 20.0, 25.0};

  h_Ev_den =
    tfs->make<TH1D>("h_Ev_den",
                    "Neutrino Energy; Neutrino Energy (GeV); Shower Reconstruction Efficiency",
                    20,
                    E_bins);
  h_Ev_den->Sumw2();
  h_Ev_num =
    tfs->make<TH1D>("h_Ev_num",
                    "Neutrino Energy; Neutrino Energy (GeV); Shower Reconstruction Efficiency",
                    20,
                    E_bins);
  h_Ev_num->Sumw2();

  h_Ee_den =
    tfs->make<TH1D>("h_Ee_den",
                    "Electron Energy; Electron Energy (GeV); Shower Reconstruction Efficiency",
                    20,
                    E_bins);
  h_Ee_den->Sumw2();
  h_Ee_num =
    tfs->make<TH1D>("h_Ee_num",
                    "Electron Energy; Electron Energy (GeV); Shower Reconstruction Efficiency",
                    20,
                    E_bins);
  h_Ee_num->Sumw2();

  if (fSaveMCTree) {
    fEventTree = new TTree("Event", "Event Tree from Sim & Reco");
    fEventTree->Branch("eventNo", &Event); // only select event in FV
    fEventTree->Branch("runNo", &Run);
    fEventTree->Branch("subRunNo", &SubRun);

    fEventTree->Branch("MC_incoming_E", &MC_incoming_E);
    fEventTree->Branch("MC_lepton_startEnergy", &MC_lepton_startEnergy);

    fEventTree->Branch("n_showers", &n_recoShowers);
    fEventTree->Branch("sh_hasPrimary_e", &sh_hasPrimary_e, "sh_hasPrimary_e[n_showers]/I");
    fEventTree->Branch("count_primary_e_in_Event", &count_primary_e_in_Event);
    fEventTree->Branch("esh_1_purity", &esh_1_purity);
    fEventTree->Branch("esh_1_completeness", &esh_1_completeness);
    fEventTree->Branch(
      "esh_each_purity", &esh_each_purity, "esh_each_purity[count_primary_e_in_Event]/D");
    fEventTree->Branch("esh_each_completeness",
                       &esh_each_completeness,
                       "esh_each_completeness[count_primary_e_in_Event]/D");
    fEventTree->Branch("esh_purity", &esh_purity);
    fEventTree->Branch("esh_completeness", &esh_completeness);
  }
}

//============================================================================
void NuShowerEff::endJob()
{
  doEfficiencies();
}

//============================================================================
void NuShowerEff::beginRun(art::Run const&)
{
  mf::LogInfo("ShowerEff") << "==== begin run ... ====" << endl;
}

//============================================================================
void NuShowerEff::analyze(art::Event const& e)
{
  reset(); // for some variables

  Event = e.id().event();
  Run = e.run();
  SubRun = e.subRun();

  // -------- find Geant4 TrackId that corresponds to e+/e- from neutrino interaction ----------
  // MCTruth: Generator
  art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
  e.getByLabel(fMCTruthModuleLabel, MCtruthHandle);
  std::vector<art::Ptr<simb::MCTruth>> MCtruthlist;
  art::fill_ptr_vector(MCtruthlist, MCtruthHandle);
  art::Ptr<simb::MCTruth> MCtruth;
  // MC (neutrino) interaction
  int MCinteractions = MCtruthlist.size();
  for (int i = 0; i < MCinteractions; i++) {
    MCtruth = MCtruthlist[i];
    if (MCtruth->NeutrinoSet()) { // NeutrinoSet(): whether the neutrino information has been set
      simb::MCNeutrino nu = MCtruth->GetNeutrino(); // GetNeutrino(): reference to neutrino info
      if (nu.CCNC() == 0)
        MC_isCC = 1;
      else if (nu.CCNC() == 1)
        MC_isCC = 0;
      simb::MCParticle neutrino = nu.Nu();           // Nu(): the incoming neutrino
      MC_target = nu.Target();                       // Target(): nuclear target, as PDG code
      MC_incoming_PDG = std::abs(nu.Nu().PdgCode()); // here not use std::abs()
      MC_Q2 = nu.QSqr();                             // QSqr(): momentum transfer Q^2, in GeV^2
      MC_channel = nu.InteractionType();             // 1001: CCQE
      MC_W = nu.W();                                 // W(): hadronic invariant mass
      const TLorentzVector& nu_momentum = nu.Nu().Momentum(0);
      nu_momentum.GetXYZT(MC_incoming_P);
      const TLorentzVector& vertex = neutrino.Position(0);
      vertex.GetXYZT(MC_vertex);
      simb::MCParticle lepton = nu.Lepton(); // Lepton(): the outgoing lepton
      MC_lepton_PDG = lepton.PdgCode();

      MC_incoming_E = MC_incoming_P[3];
    }
  }

  // Geant4: MCParticle -> lepton (e)
  // Note: generator level MCPartilceList is different from Geant4 level MCParticleList.  MCParticleList(Geant4) contains all particles in MCParticleList(generator) but their the indexes (TrackIds) are different.
  simb::MCParticle* MClepton = NULL; //Geant4 level
  art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
  simb::MCParticle* particle = 0;

  for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar) {
    particle = ipar->second; // first is index(TrackId), second is value (point address)

    auto& truth = pi_serv->ParticleToMCTruth_P(particle); // beam neutrino only
    if (truth->Origin() == simb::kBeamNeutrino && std::abs(particle->PdgCode()) == fLeptonPDGcode &&
        particle->Mother() == 0) { // primary lepton; Mother() = 0 means e^{-} for v+n=e^{-}+p
      const TLorentzVector& lepton_momentum = particle->Momentum(0);
      const TLorentzVector& lepton_position = particle->Position(0);
      const TLorentzVector& lepton_positionEnd = particle->EndPosition();
      const TLorentzVector& lepton_momentumEnd = particle->EndMomentum();
      lepton_momentum.GetXYZT(MC_lepton_startMomentum);
      lepton_position.GetXYZT(MC_lepton_startXYZT);
      lepton_positionEnd.GetXYZT(MC_lepton_endXYZT);
      lepton_momentumEnd.GetXYZT(MC_lepton_endMomentum);
      MC_leptonID = particle->TrackId();

      MC_lepton_startEnergy = MC_lepton_startMomentum[3];
      MClepton = particle;
    }
  }

  // check if the interaction is inside the Fiducial Volume
  bool isFiducial = false;
  isFiducial = insideFV(MC_vertex);
  if (!isFiducial) { return; }

  if (MC_isCC && (fNeutrinoPDGcode == std::abs(MC_incoming_PDG))) {
    if (MClepton) {
      h_Ev_den->Fill(MC_incoming_P[3]);
      h_Ee_den->Fill(MC_lepton_startMomentum[3]);
    }
  }

  // recob::Hit
  // ---- build a map for total hit charges corresponding to MC_leptonID (Geant4) ----

  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> all_hits;
  if (e.getByLabel(fHitModuleLabel, hitHandle)) { art::fill_ptr_vector(all_hits, hitHandle); }

  double ehit_total = 0.0;

  art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> fmhitmc(
    hitHandle, e, fTruthMatchDataModuleLabel);

  std::map<int, double> all_hits_trk_Q; //Q for charge: Integral()
  for (size_t i = 0; i < all_hits.size(); ++i) {
    art::Ptr<recob::Hit> hit = all_hits[i];
    auto particles =
      fmhitmc.at(hit.key()); // particles here is a pointer. A hit may come from multiple particles.
    auto hitmatch = fmhitmc.data(hit.key()); // metadata
    double maxenergy = -1e9;
    int hit_TrackId = 0;
    for (size_t j = 0; j < particles.size(); ++j) {
      if (!particles[j]) continue;
      if (!pi_serv->TrackIdToMotherParticle_P(particles[j]->TrackId())) continue;
      size_t trkid = (pi_serv->TrackIdToMotherParticle_P(particles[j]->TrackId()))->TrackId();

      if ((hitmatch[j]->energy) > maxenergy) {
        maxenergy = hitmatch[j]->energy;
        hit_TrackId = trkid;
      }
    }
    if (hit_TrackId == MC_leptonID) ehit_total += hit->Integral();
    all_hits_trk_Q[hit_TrackId] +=
      hit
        ->Integral(); // Integral(): the integral under the calibrated signal waveform of the hit, in tick x ADC units
  }

  //--------- Loop over all showers: find the associated hits for each shower ------------
  const simb::MCParticle* MClepton_reco = NULL;

  double temp_sh_ehit_Q = 0.0;
  double temp_sh_allHit_Q = 0.0;
  int temp_sh_TrackId = -999;
  count_primary_e_in_Event = 0;

  art::Handle<std::vector<recob::Shower>> showerHandle;
  std::vector<art::Ptr<recob::Shower>> showerlist;
  if (e.getByLabel(fShowerModuleLabel, showerHandle)) {
    art::fill_ptr_vector(showerlist, showerHandle);
  }
  n_recoShowers = showerlist.size();
  art::FindManyP<recob::Hit> sh_hitsAll(showerHandle, e, fShowerModuleLabel);

  std::map<int, double> showers_trk_Q;
  for (int i = 0; i < n_recoShowers; i++) { // loop over showers
    sh_hasPrimary_e[i] = 0;

    std::map<int, double> sh_hits_trk_Q; //Q for charge: Integral()
    sh_hits_trk_Q.clear();

    art::Ptr<recob::Shower> shower = showerlist[i];
    sh_direction_X[i] =
      shower->Direction().X(); // Direction(): direction cosines at start of the shower
    sh_direction_Y[i] = shower->Direction().Y();
    sh_direction_Z[i] = shower->Direction().Z();
    sh_start_X[i] = shower->ShowerStart().X();
    sh_start_Y[i] = shower->ShowerStart().Y();
    sh_start_Z[i] = shower->ShowerStart().Z();
    sh_length[i] = shower->Length(); // shower length in [cm]

    std::vector<art::Ptr<recob::Hit>> sh_hits; // associated hits for ith shower
    //In mcc8, if we get hits associated with the shower through shower->hits association directly for pandora, the hit list is incomplete. The recommended way of getting hits is through association with pfparticle:
    //shower->pfparticle->clusters->hits
    //----------getting hits through PFParticle (an option here)-------------------
    if (fHitShowerThroughPFParticle) {

      // PFParticle
      art::Handle<std::vector<recob::PFParticle>> pfpHandle;
      std::vector<art::Ptr<recob::PFParticle>> pfps;
      if (e.getByLabel(fShowerModuleLabel, pfpHandle)) art::fill_ptr_vector(pfps, pfpHandle);
      // Clusters
      art::Handle<std::vector<recob::Cluster>> clusterHandle;
      std::vector<art::Ptr<recob::Cluster>> clusters;
      if (e.getByLabel(fShowerModuleLabel, clusterHandle))
        art::fill_ptr_vector(clusters, clusterHandle);
      art::FindManyP<recob::PFParticle> fmps(
        showerHandle, e, fShowerModuleLabel); // PFParticle in Shower
      art::FindManyP<recob::Cluster> fmcp(
        pfpHandle, e, fShowerModuleLabel); // Cluster vs. PFParticle
      art::FindManyP<recob::Hit> fmhc(clusterHandle, e, fShowerModuleLabel); // Hit in Shower
      if (fmps.isValid()) {
        std::vector<art::Ptr<recob::PFParticle>> pfs = fmps.at(i);
        for (size_t ipf = 0; ipf < pfs.size(); ++ipf) {
          if (fmcp.isValid()) {
            std::vector<art::Ptr<recob::Cluster>> clus = fmcp.at(pfs[ipf].key());
            for (size_t iclu = 0; iclu < clus.size(); ++iclu) {
              if (fmhc.isValid()) {
                std::vector<art::Ptr<recob::Hit>> hits = fmhc.at(clus[iclu].key());
                for (size_t ihit = 0; ihit < hits.size(); ++ihit) {
                  sh_hits.push_back(hits[ihit]);
                }
              }
            }
          }
        }
      }
    }
    else {

      // ----- using shower->hit association for getting hits associated with shower-----
      sh_hits =
        sh_hitsAll.at(i); // associated hits for ith shower (using association of hits and shower)
    }                     // we only choose one method for hits associated with shower here.

    for (size_t k = 0; k < sh_hits.size(); ++k) {
      art::Ptr<recob::Hit> hit = sh_hits[k];
      auto particles = fmhitmc.at(hit.key());
      auto hitmatch = fmhitmc.data(hit.key());
      double maxenergy = -1e9;
      int hit_TrackId = 0;
      for (size_t j = 0; j < particles.size(); ++j) {
        if (!particles[j]) continue;
        if (!pi_serv->TrackIdToMotherParticle_P(particles[j]->TrackId())) continue;
        size_t trkid = (pi_serv->TrackIdToMotherParticle_P(particles[j]->TrackId()))->TrackId();

        if ((hitmatch[j]->energy) > maxenergy) {
          maxenergy = hitmatch[j]->energy;
          hit_TrackId = trkid;
        }
      }
      if (hit_TrackId == MC_leptonID) { sh_ehit_Q[i] += hit->Integral(); }
      sh_allhit_Q[i] += hit->Integral();
      sh_hits_trk_Q[hit_TrackId] += hit->Integral(); // Q from the same hit_TrackId
    }

    // get TrackId for a shower
    double maxshowerQ = -1.0e9;
    for (std::map<int, double>::iterator k = sh_hits_trk_Q.begin(); k != sh_hits_trk_Q.end(); k++) {
      if (k->second > maxshowerQ) {
        maxshowerQ = k->second;
        sh_TrackId[i] =
          k->first; //assign a sh_TrackId with TrackId from hit(particle) that contributing largest to the shower.
      }
    }

    //---------------------------------------------------------------------------------

    if (sh_TrackId[i] == MC_leptonID && sh_ehit_Q[i] > 0) {
      temp_sh_TrackId = sh_TrackId[i];
      sh_hasPrimary_e[i] = 1;
      count_primary_e_in_Event += 1;
      MClepton_reco = pi_serv->TrackIdToParticle_P(sh_TrackId[i]);

      if (sh_allhit_Q[i] > 0 && sh_ehit_Q[i] <= sh_allhit_Q[i]) {
        sh_purity[i] = sh_ehit_Q[i] / sh_allhit_Q[i];
      }
      else {
        sh_purity[i] = 0;
      }
      if (ehit_total > 0 && sh_ehit_Q[i] <= sh_allhit_Q[i]) {
        sh_completeness[i] = sh_ehit_Q[i] / ehit_total;
      }
      else {
        sh_completeness[i] = 0;
      }
      temp_sh_ehit_Q += sh_ehit_Q[i];
      temp_sh_allHit_Q += sh_allhit_Q[i];
    }

    showers_trk_Q[sh_TrackId[i]] += maxshowerQ;

  } // end: for loop over showers

  // ---------------------------------------------------------------
  if (MClepton_reco && MClepton) {
    if ((temp_sh_TrackId == MC_leptonID) && MC_isCC &&
        (fNeutrinoPDGcode == std::abs(MC_incoming_PDG))) {
      if ((temp_sh_allHit_Q >= temp_sh_ehit_Q) && (temp_sh_ehit_Q > 0.0)) {
        esh_purity = temp_sh_ehit_Q / temp_sh_allHit_Q;
        esh_completeness = temp_sh_ehit_Q / ehit_total;

        if (esh_purity > fMinPurity && esh_completeness > fMinCompleteness) {
          h_Ev_num->Fill(MC_incoming_P[3]);
          h_Ee_num->Fill(MC_lepton_startMomentum[3]);
        }
      }
    }
  }
  // --------------------------------------------------------------

  if ((MClepton_reco && MClepton) && MC_isCC && (fNeutrinoPDGcode == std::abs(MC_incoming_PDG))) {
    for (int j = 0; j < count_primary_e_in_Event; j++) {
      esh_each_purity[j] = 0.0;
    }

    double temp_esh_1_allhit = -1.0e9;
    int temp_shower_index = -999;
    int temp_esh_index = 0;
    for (int i = 0; i < n_recoShowers; i++) {
      if (sh_TrackId[i] == MC_leptonID) {

        // for each electron shower
        if (sh_ehit_Q[i] > 0) {
          esh_each_purity[temp_esh_index] = sh_purity[i];
          esh_each_completeness[temp_esh_index] = sh_completeness[i];
          temp_esh_index += 1;
        }

        // find largest shower
        if ((sh_allhit_Q[i] > temp_esh_1_allhit) && (sh_ehit_Q[i] > 0.0)) {
          temp_esh_1_allhit = sh_allhit_Q[i];
          temp_shower_index = i;
        }
      }
    }
    // largest shower with electron contribution
    esh_1_purity = sh_purity[temp_shower_index];
    esh_1_completeness = sh_completeness[temp_shower_index];
  }

  if (count_primary_e_in_Event > 0 && MC_isCC && fSaveMCTree) {
    fEventTree->Fill();
  } // so far, only save CC event
}

// ====================================================================================
void NuShowerEff::doEfficiencies()
{
  art::ServiceHandle<art::TFileService const> tfs;

  if (TEfficiency::CheckConsistency(*h_Ev_num, *h_Ev_den)) {
    h_Eff_Ev = tfs->make<TEfficiency>(*h_Ev_num, *h_Ev_den);
    TGraphAsymmErrors* grEff_Ev = h_Eff_Ev->CreateGraph();
    grEff_Ev->Write("grEff_Ev");
    h_Eff_Ev->Write("h_Eff_Ev");
  }

  if (TEfficiency::CheckConsistency(*h_Ee_num, *h_Ee_den)) {
    h_Eff_Ee = tfs->make<TEfficiency>(*h_Ee_num, *h_Ee_den);
    TGraphAsymmErrors* grEff_Ee = h_Eff_Ee->CreateGraph();
    grEff_Ee->Write("grEff_Ee");
    h_Eff_Ee->Write("h_Eff_Ee");
  }
}

// ====================================================================================
bool NuShowerEff::insideFV(double vertex[4])
{
  double x = vertex[0];
  double y = vertex[1];
  double z = vertex[2];

  if (x > fFidVolXmin && x < fFidVolXmax && y > fFidVolYmin && y < fFidVolYmax && z > fFidVolZmin &&
      z < fFidVolZmax) {
    return true;
  }
  else {
    return false;
  }
}

// ====================================================================================
void NuShowerEff::reset()
{
  MC_incoming_PDG = -999;
  MC_lepton_PDG = -999;
  MC_isCC = -999;
  MC_channel = -999;
  MC_target = -999;
  MC_Q2 = -999.0;
  MC_W = -999.0;
  MC_lepton_theta = -999.0;
  MC_leptonID = -999;
  MC_LeptonTrack = -999;

  for (int i = 0; i < MAX_SHOWERS; i++) {
    sh_direction_X[i] = -999.0;
    sh_direction_Y[i] = -999.0;
    sh_direction_Z[i] = -999.0;
    sh_start_X[i] = -999.0;
    sh_start_Y[i] = -999.0;
    sh_start_Z[i] = -999.0;
    sh_length[i] = -999.0;
    sh_hasPrimary_e[i] = -999;
    sh_ehit_Q[i] = 0.0;
    sh_allhit_Q[i] = 0.0;
    sh_purity[i] = 0.0;
    sh_completeness[i] = 0.0;
  }
}
DEFINE_ART_MODULE(NuShowerEff)
