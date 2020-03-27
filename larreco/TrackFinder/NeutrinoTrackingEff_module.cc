//
//**Tracking Efficiency module***
//The basic idea is to loop over the hits from a given track and call BackTracker
//then look at std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackID(hit);
//then associete the hits to a G4 track ID (particle) that generate those hits(electrons)
//It was developed for CC neutrio interactions, it also can handle proton decay events p->k+nu_bar
//And protons, pion and muons from particle cannon by using the option isNeutrinoInt = false;
//All histograms are as a function of truth info (momentum, length)
//
// A. Higuera
// ahiguera@central.uh.edu

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"

// C/C++ standard libraries

#define MAX_TRACKS 1000
using namespace std;

//========================================================================

namespace DUNE {

  class NeutrinoTrackingEff : public art::EDAnalyzer {
  public:
    explicit NeutrinoTrackingEff(fhicl::ParameterSet const& pset);

  private:
    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void processEff(const art::Event& evt);
    void truthMatcher(detinfo::DetectorClocksData const& clockData,
                      std::vector<art::Ptr<recob::Hit>> all_hits,
                      std::vector<art::Ptr<recob::Hit>> track_hits,
                      const simb::MCParticle*& MCparticle,
                      double& Efrac,
                      double& Ecomplet);
    double truthLength(const detinfo::DetectorClocksData& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       const simb::MCParticle* MCparticle);
    bool insideFV(double vertex[4]);
    void doEfficiencies();

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    int fNeutrinoPDGcode;
    int fLeptonPDGcode;
    double fMaxNeutrinoE;
    double fMaxLeptonP;
    bool fisNeutrinoInt;

    int MC_isCC;
    int MC_incoming_PDG;
    double MC_incoming_P[4];
    double MC_vertex[4];
    double MC_lepton_startMomentum[4];

    int MC_leading_protonID;
    int MC_leading_PionPlusID;
    int MC_leading_PionMinusID;
    int MC_leptonID;
    int MC_kaonID;
    int MC_michelID;

    double MC_leptonP;
    double MC_leading_PionPlusP;
    double MC_leading_ProtonP;
    double MC_leading_PionMinusP;
    double MC_kaonP;
    double MC_michelP;

    TH1D* h_Ev_den;
    TH1D* h_Ev_num;
    TH1D* h_Pmu_den;
    TH1D* h_Pmu_num;
    TH1D* h_theta_den;
    TH1D* h_theta_num;
    TH1D* h_Pproton_den;
    TH1D* h_Pproton_num;
    TH1D* h_Ppion_plus_den;
    TH1D* h_Ppion_plus_num;
    TH1D* h_Ppion_minus_den;
    TH1D* h_Ppion_minus_num;

    TH1D* h_Efrac_lepton;
    TH1D* h_Ecomplet_lepton;
    TH1D* h_Efrac_proton;
    TH1D* h_Ecomplet_proton;
    TH1D* h_Efrac_pion_plus;
    TH1D* h_Ecomplet_pion_plus;
    TH1D* h_Efrac_pion_minus;
    TH1D* h_Ecomplet_pion_minus;
    TH1D* h_trackRes_lepton;
    TH1D* h_trackRes_proton;
    TH1D* h_trackRes_pion_plus;
    TH1D* h_trackRes_pion_minus;

    TH1D* h_muon_length;
    TH1D* h_proton_length;
    TH1D* h_pionp_length;
    TH1D* h_pionm_length;
    TH1D* h_muonwtrk_length;
    TH1D* h_protonwtrk_length;
    TH1D* h_pionpwtrk_length;
    TH1D* h_pionmwtrk_length;

    TEfficiency* h_Eff_Ev = 0;
    TEfficiency* h_Eff_Pmu = 0;
    TEfficiency* h_Eff_theta = 0;
    TEfficiency* h_Eff_Pproton = 0;
    TEfficiency* h_Eff_Ppion_plus = 0;
    TEfficiency* h_Eff_Ppion_minus = 0;

    TEfficiency* h_Eff_Lmuon = 0;
    TEfficiency* h_Eff_Lproton = 0;
    TEfficiency* h_Eff_Lpion_plus = 0;
    TEfficiency* h_Eff_Lpion_minus = 0;

    //nucleon decay histograms
    TH1D* h_Pkaon_den;
    TH1D* h_Pkaon_num;
    TH1D* h_Pmichel_e_den;
    TH1D* h_Pmichel_e_num;
    TH1D* h_Efrac_kaon;
    TH1D* h_Ecomplet_kaon;
    TH1D* h_trackRes_kaon;
    TH1D* h_Efrac_michel;
    TH1D* h_Ecomplet_michel;
    TH1D* h_trackRes_michel;
    TH1D* h_kaon_length;
    TH1D* h_michel_length;
    TH1D* h_kaonwtrk_length;
    TH1D* h_michelwtrk_length;
    TEfficiency* h_Eff_Pkaon = 0;
    TEfficiency* h_Eff_Pmichel = 0;
    TEfficiency* h_Eff_Lkaon = 0;
    TEfficiency* h_Eff_Lmichel = 0;

    float fFidVolCutX;
    float fFidVolCutY;
    float fFidVolCutZ;

    float fFidVolXmin;
    float fFidVolXmax;
    float fFidVolYmin;
    float fFidVolYmax;
    float fFidVolZmin;
    float fFidVolZmax;

    double fDriftVelocity; // in cm/ns
    art::ServiceHandle<geo::Geometry const> geom;

  }; // class NeutrinoTrackingEff

  //========================================================================
  NeutrinoTrackingEff::NeutrinoTrackingEff(fhicl::ParameterSet const& p) : EDAnalyzer(p)
  {
    fMCTruthModuleLabel = p.get<std::string>("MCTruthModuleLabel");
    fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
    fisNeutrinoInt = p.get<bool>("isNeutrinoInt");
    fLeptonPDGcode = p.get<int>("LeptonPDGcode");
    fNeutrinoPDGcode = p.get<int>("NeutrinoPDGcode");
    fMaxNeutrinoE = p.get<double>("MaxNeutrinoE");
    fMaxLeptonP = p.get<double>("MaxLeptonP");
    fFidVolCutX = p.get<float>("FidVolCutX");
    fFidVolCutY = p.get<float>("FidVolCutY");
    fFidVolCutZ = p.get<float>("FidVolCutZ");

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    fDriftVelocity = detProp.DriftVelocity() * 1e-3; // cm/ns
  }
  //========================================================================
  void
  NeutrinoTrackingEff::beginJob()
  {
    std::cout << "job begin..." << std::endl;
    auto const* geo = lar::providerFrom<geo::Geometry>();
    // Define histogram boundaries (cm).
    // For now only draw cryostat=0.
    double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;
    for (size_t i = 0; i < geo->NTPC(); ++i) {
      double local[3] = {0., 0., 0.};
      double world[3] = {0., 0., 0.};
      const geo::TPCGeo& tpc = geo->TPC(i);
      tpc.LocalToWorld(local, world);
      if (minx > world[0] - geo->DetHalfWidth(i)) minx = world[0] - geo->DetHalfWidth(i);
      if (maxx < world[0] + geo->DetHalfWidth(i)) maxx = world[0] + geo->DetHalfWidth(i);
      if (miny > world[1] - geo->DetHalfHeight(i)) miny = world[1] - geo->DetHalfHeight(i);
      if (maxy < world[1] + geo->DetHalfHeight(i)) maxy = world[1] + geo->DetHalfHeight(i);
      if (minz > world[2] - geo->DetLength(i) / 2.) minz = world[2] - geo->DetLength(i) / 2.;
      if (maxz < world[2] + geo->DetLength(i) / 2.) maxz = world[2] + geo->DetLength(i) / 2.;
    }

    fFidVolXmin = minx + fFidVolCutX;
    fFidVolXmax = maxx - fFidVolCutX;
    fFidVolYmin = miny + fFidVolCutY;
    fFidVolYmax = maxy - fFidVolCutY;
    fFidVolZmin = minz + fFidVolCutZ;
    fFidVolZmax = maxz - fFidVolCutZ;

    std::cout << "Fiducial volume:"
              << "\n"
              << fFidVolXmin << "\t< x <\t" << fFidVolXmax << "\n"
              << fFidVolYmin << "\t< y <\t" << fFidVolYmax << "\n"
              << fFidVolZmin << "\t< z <\t" << fFidVolZmax << "\n";

    art::ServiceHandle<art::TFileService const> tfs;

    double E_bins[21] = {0,   0.5, 1.0, 1.5, 2.0,  2.5,  3.0,  3.5,  4,    4.5, 5.0,
                         5.5, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 17.0, 20.0, 25.0};
    double theta_bin[44] = {0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10.,
                            11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22.,
                            24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44.,
                            46., 48., 50., 55., 60., 65., 70., 75., 80., 85., 90.};
    double Pbins[18] = {
      0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0};

    for (int i = 0; i < 21; ++i)
      E_bins[i] *= fMaxNeutrinoE / 25.;
    for (int i = 0; i < 18; ++i)
      Pbins[i] *= fMaxLeptonP / 3.0;

    h_Ev_den = tfs->make<TH1D>(
      "h_Ev_den", "Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency", 20, E_bins);
    h_Ev_num = tfs->make<TH1D>(
      "h_Ev_num", "Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency", 20, E_bins);
    h_Pmu_den = tfs->make<TH1D>(
      "h_Pmu_den", "Muon Momentum; Muon Momentum (GeV); Tracking Efficiency", 20, E_bins);
    h_Pmu_num = tfs->make<TH1D>(
      "h_Pmu_num", "Muon Momentum; Muon Momentum (GeV); Tracking Efficiency", 20, E_bins);
    h_theta_den =
      tfs->make<TH1D>("h_theta_den",
                      "Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",
                      43,
                      theta_bin);
    h_theta_num =
      tfs->make<TH1D>("h_theta_num",
                      "Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",
                      43,
                      theta_bin);
    h_Pproton_den = tfs->make<TH1D>(
      "h_Pproton_den", "Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
    h_Pproton_num = tfs->make<TH1D>(
      "h_Pproton_num", "Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
    h_Ppion_plus_den = tfs->make<TH1D>(
      "h_Ppion_plus_den", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ppion_plus_num = tfs->make<TH1D>(
      "h_Ppion_plus_num", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ppion_minus_den = tfs->make<TH1D>(
      "h_Ppion_minus_den", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ppion_minus_num = tfs->make<TH1D>(
      "h_Ppion_minus_num", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ev_den->Sumw2();
    h_Ev_num->Sumw2();
    h_Pmu_den->Sumw2();
    h_Pmu_num->Sumw2();
    h_theta_den->Sumw2();
    h_theta_num->Sumw2();
    h_Pproton_den->Sumw2();
    h_Pproton_num->Sumw2();
    h_Ppion_plus_den->Sumw2();
    h_Ppion_plus_num->Sumw2();
    h_Ppion_minus_den->Sumw2();
    h_Ppion_minus_num->Sumw2();

    h_Efrac_lepton = tfs->make<TH1D>("h_Efrac_lepton", "Efrac Lepton; Track Purity;", 60, 0, 1.2);
    h_Ecomplet_lepton =
      tfs->make<TH1D>("h_Ecomplet_lepton", "Ecomplet Lepton; Track Completeness;", 60, 0, 1.2);
    h_Efrac_proton = tfs->make<TH1D>("h_Efrac_proton", "Efrac Proton; Track Purity;", 60, 0, 1.2);
    h_Ecomplet_proton =
      tfs->make<TH1D>("h_Ecomplet_proton", "Ecomplet Proton; Track Completeness;", 60, 0, 1.2);
    h_Efrac_pion_plus =
      tfs->make<TH1D>("h_Efrac_pion_plus", "Efrac Pion +; Track Purity;", 60, 0, 1.2);
    h_Ecomplet_pion_plus =
      tfs->make<TH1D>("h_Ecomplet_pion_plus", "Ecomplet Pion +; Track Completeness;", 60, 0, 1.2);
    h_Efrac_pion_minus =
      tfs->make<TH1D>("h_Efrac_pion_minus", "Efrac Pion -; Track Purity;", 60, 0, 1.2);
    h_Ecomplet_pion_minus =
      tfs->make<TH1D>("h_Ecomplet_pion_minus", "Ecomplet Pion -; Track Completeness;", 60, 0, 1.2);
    h_trackRes_lepton = tfs->make<TH1D>(
      "h_trackRes_lepton", "Muon Residual; Truth length - Reco length (cm);", 200, -100, 100);
    h_trackRes_proton = tfs->make<TH1D>(
      "h_trackRes_proton", "Proton Residual; Truth length - Reco length (cm);", 200, -100, 100);
    h_trackRes_pion_plus = tfs->make<TH1D>(
      "h_trackRes_pion_plus", "Pion + Residual; Truth length - Reco length (cm);", 200, -100, 100);
    h_trackRes_pion_minus = tfs->make<TH1D>(
      "h_trackRes_pion_minus", "Pion - Residual; Truth length - Reco length (cm);", 200, -100, 100);
    h_Efrac_lepton->Sumw2();
    h_Ecomplet_lepton->Sumw2();
    h_Efrac_proton->Sumw2();
    h_Ecomplet_proton->Sumw2();
    h_Efrac_pion_plus->Sumw2();
    h_Ecomplet_pion_plus->Sumw2();
    h_Efrac_pion_minus->Sumw2();
    h_Ecomplet_pion_minus->Sumw2();
    h_trackRes_lepton->Sumw2();
    h_trackRes_proton->Sumw2();
    h_trackRes_pion_plus->Sumw2();
    h_trackRes_pion_minus->Sumw2();

    h_muon_length =
      tfs->make<TH1D>("h_muon_length", "Muon Length;  Muon Truth Length (cm)", 40, 0, 100);
    h_proton_length =
      tfs->make<TH1D>("h_proton_length", "Proton Length; Proton Truth Length (cm)", 40, 0, 100);
    h_pionp_length =
      tfs->make<TH1D>("h_pionp_length", "Pion + Length; Pion^{+} Truth Length (cm)", 40, 0, 100);
    h_pionm_length =
      tfs->make<TH1D>("h_pionm_length", "Pion - Length; Pion^{-} Truth Length (cm)", 40, 0, 100);

    h_muonwtrk_length =
      tfs->make<TH1D>("h_muonwtrk_length", "Muon Length; Muon Truth Length (cm)", 40, 0, 100);
    h_protonwtrk_length =
      tfs->make<TH1D>("h_protonwtrk_length", "Proton Length; Proton Truth Length (cm)", 40, 0, 100);
    h_pionpwtrk_length = tfs->make<TH1D>(
      "h_pionpwtrk_length", "Pion + Length; Pion^{+} Truth Length (cm)", 40, 0, 100);
    h_pionmwtrk_length = tfs->make<TH1D>(
      "h_pionmwtrk_length", "Pion - Length; Pion^{-} Truth Length (cm)", 40, 0, 100);

    h_muon_length->Sumw2();
    h_muonwtrk_length->Sumw2();
    h_proton_length->Sumw2();
    h_protonwtrk_length->Sumw2();
    h_pionp_length->Sumw2();
    h_pionpwtrk_length->Sumw2();
    h_pionm_length->Sumw2();
    h_pionmwtrk_length->Sumw2();

    h_Pkaon_den =
      tfs->make<TH1D>("h_Pkaon_den", "Kaon; Kaon Momentum (GeV); Tracking Efficiency", 17, Pbins);
    h_Pkaon_num =
      tfs->make<TH1D>("h_Pkaon_num", "Kaon; Kaon Momentum (GeV); Tracking Efficiency", 17, Pbins);
    h_Pmichel_e_den =
      tfs->make<TH1D>("h_Pmichel_e_den",
                      "Michel Electron; Michele e Momentum (GeV); Tracking Efficiency",
                      17,
                      Pbins);
    h_Pmichel_e_num =
      tfs->make<TH1D>("h_Pmichel_e_num",
                      "Michel Electron; Michele e Momentum (GeV); Tracking Efficiency",
                      17,
                      Pbins);
    h_Pkaon_den->Sumw2();
    h_Pkaon_num->Sumw2();
    h_Pmichel_e_den->Sumw2();
    h_Pmichel_e_num->Sumw2();
    h_Efrac_kaon = tfs->make<TH1D>("h_Efrac_kaon", "Efrac Kaon; Track Purity;", 60, 0, 1.2);
    h_Ecomplet_kaon =
      tfs->make<TH1D>("h_Ecomplet_kaon", "Ecomplet Kaon; Track Completeness;", 60, 0, 1.2);
    h_trackRes_kaon = tfs->make<TH1D>(
      "h_trackRes_kaon", "Kaon Residual; Truth length - Reco length (cm);", 200, -100, 100);
    h_Efrac_michel =
      tfs->make<TH1D>("h_Efrac_michel", "Efrac Michel; Track Energy fraction;", 60, 0, 1.2);
    h_Ecomplet_michel =
      tfs->make<TH1D>("h_Ecomplet_michel", "Ecomplet Michel; Track Completeness;", 60, 0, 1.2);
    h_trackRes_michel = tfs->make<TH1D>(
      "h_trackRes_michel", "Michel Residual; Truth length - Reco length (cm);", 200, -100, 100);
    h_kaon_length =
      tfs->make<TH1D>("h_kaon_length", "Kaon Length; Kaon Truth Length (cm)", 40, 0, 100);
    h_kaonwtrk_length =
      tfs->make<TH1D>("h_kaonwtrk_length", "Kaon Length; Kaon Truth Length (cm)", 40, 0, 100);
    h_michel_length =
      tfs->make<TH1D>("h_michel_length", "Michel Length; Michel e Truth Length (cm)", 40, 0, 100);
    h_michelwtrk_length = tfs->make<TH1D>(
      "h_michelwtrk_length", "Michel Length; Michel e Truth Length (cm)", 40, 0, 100);

    h_Efrac_kaon->Sumw2();
    h_Ecomplet_kaon->Sumw2();
    h_trackRes_kaon->Sumw2();
    h_Efrac_michel->Sumw2();
    h_Ecomplet_michel->Sumw2();
    h_trackRes_michel->Sumw2();
    h_kaon_length->Sumw2();
    h_kaonwtrk_length->Sumw2();
    h_michel_length->Sumw2();
    h_michelwtrk_length->Sumw2();
  }
  //========================================================================
  void
  NeutrinoTrackingEff::endJob()
  {
    doEfficiencies();
  }
  //========================================================================
  void
  NeutrinoTrackingEff::beginRun(const art::Run& /*run*/)
  {
    mf::LogInfo("NeutrinoTrackingEff") << "begin run..." << std::endl;
  }
  //========================================================================
  void
  NeutrinoTrackingEff::analyze(const art::Event& event)
  {
    if (event.isRealData()) return;

    processEff(event);
  }
  //========================================================================
  void
  NeutrinoTrackingEff::processEff(const art::Event& event)
  {
    // Save neutrino's interaction info
    art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
    event.getByLabel(fMCTruthModuleLabel, MCtruthHandle);
    std::vector<art::Ptr<simb::MCTruth>> MCtruthlist;
    art::fill_ptr_vector(MCtruthlist, MCtruthHandle);
    art::Ptr<simb::MCTruth> MCtruth;

    //For now assume that there is only one neutrino interaction...
    MCtruth = MCtruthlist[0];
    if (MCtruth->NeutrinoSet()) {
      simb::MCNeutrino nu = MCtruth->GetNeutrino();
      if (nu.CCNC() == 0)
        MC_isCC = 1;
      else if (nu.CCNC() == 1)
        MC_isCC = 0;
      simb::MCParticle neutrino = nu.Nu();
      MC_incoming_PDG = nu.Nu().PdgCode();
      const TLorentzVector& nu_momentum = nu.Nu().Momentum(0);
      nu_momentum.GetXYZT(MC_incoming_P);
      const TLorentzVector& vertex = neutrino.Position(0);
      vertex.GetXYZT(MC_vertex);
    }

    //!save FS particles

    double tmp_leadingPionPlusE = 0.0;
    double tmp_leadingPionMinusE = 0.0;
    double tmp_leadingProtonE = 0.0;

    simb::MCParticle* MClepton = nullptr;
    simb::MCParticle* MCproton = nullptr;
    simb::MCParticle* MCpion_plus = nullptr;
    simb::MCParticle* MCpion_minus = nullptr;
    simb::MCParticle* MCkaon = nullptr;
    simb::MCParticle* MCmichel = nullptr;

    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    simb::MCParticle* particle = 0;
    int i = 0; // particle index

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);

    for (sim::ParticleList::const_iterator ipar = plist.begin(); ipar != plist.end(); ++ipar) {
      particle = ipar->second;
      if (particle->PdgCode() == fLeptonPDGcode && particle->Mother() == 0) { //primary lepton
        const TLorentzVector& lepton_momentum = particle->Momentum(0);
        lepton_momentum.GetXYZT(MC_lepton_startMomentum);
        MC_leptonID = particle->TrackId();
        MC_leptonP = sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                          pow(particle->Momentum().Pz(), 2));
        MClepton = particle;
      }
      if (particle->Mother() == 0) { //save primary particle i.e. from the neutrino interaction
        //save leading pion and proton
        if (particle->PdgCode() == 2212) {
          if (particle->Momentum().E() > tmp_leadingProtonE) {
            tmp_leadingProtonE = particle->Momentum().E();
            MC_leading_protonID = particle->TrackId();
            MC_leading_ProtonP =
              sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                   pow(particle->Momentum().Pz(), 2));
            MCproton = particle;
          }
        }
        else if (particle->PdgCode() == 211) {
          if (particle->Momentum().E() > tmp_leadingPionPlusE) {
            tmp_leadingPionPlusE = particle->Momentum().E();
            MC_leading_PionPlusID = particle->TrackId();
            MC_leading_PionPlusP =
              sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                   pow(particle->Momentum().Pz(), 2));
            MCpion_plus = particle;
          }
        }
        else if (particle->PdgCode() == -211) {
          if (particle->Momentum().E() > tmp_leadingPionMinusE) {
            tmp_leadingPionMinusE = particle->Momentum().E();
            MC_leading_PionMinusID = particle->TrackId();
            MC_leading_PionMinusP =
              sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                   pow(particle->Momentum().Pz(), 2));
            MCpion_minus = particle;
          }
        }
        i++; //paticle index
      }

      //=======================================================================================
      //add Nucleon decay stuff and particle cannon
      //=======================================================================================
      if (!fisNeutrinoInt) {
        if (particle->Mother() == 0) {
          const TLorentzVector& positionStart = particle->Position(0);
          positionStart.GetXYZT(MC_vertex);
        }
        if (particle->PdgCode() == 321) { //save primary Kaon
          MC_kaonID = particle->TrackId();
          MC_kaonP = sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                          pow(particle->Momentum().Pz(), 2));
          MCkaon = particle;
        }
        else if (particle->PdgCode() == fLeptonPDGcode) { // Particle cannon muon
          MC_leptonID = particle->TrackId();
          MC_leptonP = sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                            pow(particle->Momentum().Pz(), 2));
          MClepton = particle;
        }
        else if (particle->PdgCode() == 2212) {
          if (particle->Momentum().E() > tmp_leadingProtonE) {
            tmp_leadingProtonE = particle->Momentum().E();
            MC_leading_protonID = particle->TrackId();
            MC_leading_ProtonP =
              sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                   pow(particle->Momentum().Pz(), 2));
            MCproton = particle;
          }
        }
        else if (particle->PdgCode() == 211) {
          if (particle->Momentum().E() > tmp_leadingPionPlusE) {
            tmp_leadingPionPlusE = particle->Momentum().E();
            MC_leading_PionPlusID = particle->TrackId();
            MC_leading_PionPlusP =
              sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                   pow(particle->Momentum().Pz(), 2));
            MCpion_plus = particle;
          }
        }
        else if (particle->PdgCode() == -211) {
          if (particle->Momentum().E() > tmp_leadingPionMinusE) {
            tmp_leadingPionMinusE = particle->Momentum().E();
            MC_leading_PionMinusID = particle->TrackId();
            MC_leading_PionMinusP =
              sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                   pow(particle->Momentum().Pz(), 2));
            MCpion_minus = particle;
          }
        }
        else if (particle->Process() == "Decay" &&
                 particle->PdgCode() == -11) { // michel electron from muon decay
          MC_michelID = particle->TrackId();
          MC_michelP = sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                            pow(particle->Momentum().Pz(), 2));
          MCmichel = particle;
        }
        else if (TMath::Abs(particle->PdgCode() == 321)) { //save primary Kaon
          MC_kaonID = particle->TrackId();
          MC_kaonP = sqrt(pow(particle->Momentum().Px(), 2) + pow(particle->Momentum().Py(), 2) +
                          pow(particle->Momentum().Pz(), 2));
          MCkaon = particle;
        }
      }
    }
    //===================================================================
    //Saving denominator histograms
    //===================================================================
    if (not insideFV(MC_vertex)) return;
    double Pv =
      sqrt(pow(MC_incoming_P[0], 2) + pow(MC_incoming_P[1], 2) + pow(MC_incoming_P[2], 2));
    double theta_mu = acos((MC_incoming_P[0] * MC_lepton_startMomentum[0] +
                            MC_incoming_P[1] * MC_lepton_startMomentum[1] +
                            MC_incoming_P[2] * MC_lepton_startMomentum[2]) /
                           (Pv * MC_leptonP));
    theta_mu *= (180.0 / 3.14159);
    double truth_lengthLepton = truthLength(clockData, detProp, MClepton);
    double proton_length = truthLength(clockData, detProp, MCproton);
    double pion_plus_length = truthLength(clockData, detProp, MCpion_plus);
    double pion_minus_length = truthLength(clockData, detProp, MCpion_minus);
    double kaonLength = truthLength(clockData, detProp, MCkaon);
    double michelLength = truthLength(clockData, detProp, MCmichel);

    // save CC events within the fiducial volume with the favorite neutrino
    // flavor
    if (MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE)) {
      if (MClepton) {
        h_Ev_den->Fill(MC_incoming_P[3]);
        h_Pmu_den->Fill(MC_leptonP);
        h_theta_den->Fill(theta_mu);
        h_muon_length->Fill(truth_lengthLepton);
      }
      if (MCproton) {
        h_Pproton_den->Fill(MC_leading_ProtonP);
        h_proton_length->Fill(proton_length);
      }
      if (MCpion_plus) {
        h_Ppion_plus_den->Fill(MC_leading_PionPlusP);
        h_pionp_length->Fill(pion_plus_length);
      }
      if (MCpion_minus) {
        h_Ppion_minus_den->Fill(MC_leading_PionMinusP);
        h_pionm_length->Fill(pion_minus_length);
      }
      if (MCkaon) {
        h_Pkaon_den->Fill(MC_kaonP);
        h_kaon_length->Fill(kaonLength);
      }
    }

    //save events for Nucleon decay and particle cannon
    if (!fisNeutrinoInt) {
      if (MClepton) {
        h_Pmu_den->Fill(MC_leptonP);
        h_muon_length->Fill(truth_lengthLepton);
      }
      if (MCkaon) {
        h_Pkaon_den->Fill(MC_kaonP);
        h_kaon_length->Fill(kaonLength);
      }
      if (MCproton) {
        h_Pproton_den->Fill(MC_leading_ProtonP);
        h_proton_length->Fill(proton_length);
      }
      if (MCpion_plus) {
        h_Ppion_plus_den->Fill(MC_leading_PionPlusP);
        h_pionp_length->Fill(pion_plus_length);
      }
      if (MCpion_minus) {
        h_Ppion_minus_den->Fill(MC_leading_PionMinusP);
        h_pionm_length->Fill(pion_minus_length);
      }
      if (MCmichel) {
        h_Pmichel_e_den->Fill(MC_michelP);
        h_michel_length->Fill(michelLength);
      }
    }

    //========================================================================
    // Reco stuff, once we have selected a MC Particle let's find out if there is a track associated
    //========================================================================
    art::Handle<std::vector<recob::Track>> trackListHandle;
    if (!event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track>> tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);
    int n_recoTrack = tracklist.size();

    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    if (n_recoTrack == 0) {
      MF_LOG_DEBUG("NeutrinoTrackingEff") << "There are no reco tracks... bye";
      return;
    }
    MF_LOG_DEBUG("NeutrinoTrackingEff") << "Found this many reco tracks " << n_recoTrack;

    double Efrac_lepton = 0.0;
    double Ecomplet_lepton = 0.0;
    double Efrac_proton = 0.0;
    double Ecomplet_proton = 0.0;
    double Efrac_pionplus = 0.0;
    double Ecomplet_pionplus = 0.0;
    double Efrac_pionminus = 0.0;
    double Ecomplet_pionminus = 0.0;
    double Efrac_kaon = 0.0;
    double Ecomplet_kaon = 0.0;
    double Efrac_michel = 0.0;
    double Ecomplet_michel = 0.0;
    double trackLength_lepton = 0.0;
    double trackLength_proton = 0.0;
    double trackLength_pion_plus = 0.0;
    double trackLength_pion_minus = 0.0;
    double trackLength_kaon = 0.0;
    double trackLength_michel = 0.0;
    const simb::MCParticle* MClepton_reco = nullptr;
    const simb::MCParticle* MCproton_reco = nullptr;
    const simb::MCParticle* MCpion_plus_reco = nullptr;
    const simb::MCParticle* MCpion_minus_reco = nullptr;
    const simb::MCParticle* MCkaon_reco = nullptr;
    const simb::MCParticle* MCmichel_reco = nullptr;

    std::vector<art::Ptr<recob::Hit>> tmp_all_trackHits = track_hits.at(0);
    std::vector<art::Ptr<recob::Hit>> all_hits;
    art::Handle<std::vector<recob::Hit>> hithandle;
    auto const pd = event.getProductDescription(tmp_all_trackHits[0].id());
    if (pd && event.getByLabel(pd->inputTag(), hithandle)) {
      art::fill_ptr_vector(all_hits, hithandle);
    }

    for (int i = 0; i < n_recoTrack; i++) {
      art::Ptr<recob::Track> track = tracklist[i];
      std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);
      double tmpEfrac = 0;
      double tmpEcomplet = 0;
      const simb::MCParticle* particle;
      truthMatcher(clockData, all_hits, all_trackHits, particle, tmpEfrac, tmpEcomplet);
      if (!particle) continue;
      if ((particle->PdgCode() == fLeptonPDGcode) && (particle->TrackId() == MC_leptonID)) {
        // save the best track ... based on completeness if there is more than
        // one track if( tmpEfrac > Efrac_lepton ){ ///this was base on purity
        if (tmpEcomplet > Ecomplet_lepton) {
          Ecomplet_lepton = tmpEcomplet;
          Efrac_lepton = tmpEfrac;
          trackLength_lepton = track->Length();
          MClepton_reco = particle;
        }
      }
      else if ((particle->PdgCode() == 2212) && (particle->TrackId() == MC_leading_protonID)) {
        //save the best track ... based on completeness if there is more than one track
        if (tmpEcomplet > Ecomplet_proton) {
          Ecomplet_proton = tmpEcomplet;
          Efrac_proton = tmpEfrac;
          trackLength_proton = track->Length();
          MCproton_reco = particle;
        }
      }
      else if ((particle->PdgCode() == 211) && (particle->TrackId() == MC_leading_PionPlusID)) {
        //save the best track ... based on completeness if there is more than one track
        if (tmpEcomplet > Ecomplet_pionplus) {
          Ecomplet_pionplus = tmpEcomplet;
          Efrac_pionplus = tmpEfrac;
          trackLength_pion_plus = track->Length();
          MCpion_plus_reco = particle;
        }
      }
      else if ((particle->PdgCode() == -211) && (particle->TrackId() == MC_leading_PionMinusID)) {
        //save the best track ... based on completeness if there is more than one track
        if (tmpEcomplet > Ecomplet_pionminus) {
          Ecomplet_pionminus = tmpEcomplet;
          Efrac_pionminus = tmpEfrac;
          trackLength_pion_minus = track->Length();
          MCpion_minus_reco = particle;
        }
      }
      //kaon from nucleon decay
      else if ((TMath::Abs(particle->PdgCode()) == 321) && (particle->TrackId() == MC_kaonID)) {
        //save the best track ... based on completeness if there is more than one track
        if (tmpEcomplet > Ecomplet_kaon) {
          Ecomplet_kaon = tmpEcomplet;
          Efrac_kaon = tmpEfrac;
          trackLength_kaon = track->Length();
          MCkaon_reco = particle;
        }
      }
      //michel from nucleon decay
      else if ((particle->PdgCode() == -11) && (particle->TrackId() == MC_michelID)) {
        //save the best track ... based on completeness if there is more than one track
        if (tmpEcomplet > Ecomplet_michel) {
          Ecomplet_michel = tmpEcomplet;
          Efrac_michel = tmpEfrac;
          trackLength_michel = track->Length();
          MCmichel_reco = particle;
        }
      }
    }

    double Reco_LengthRes = truth_lengthLepton - trackLength_lepton;
    double Reco_LengthResProton = proton_length - trackLength_proton;
    double Reco_LengthResPionPlus = pion_plus_length - trackLength_pion_plus;
    double Reco_LengthResPionMinus = pion_minus_length - trackLength_pion_minus;

    if (MClepton_reco && MClepton) {
      if (MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE)) {
        h_Pmu_num->Fill(MC_leptonP);
        h_Ev_num->Fill(MC_incoming_P[3]);
        h_theta_num->Fill(theta_mu);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_Ecomplet_lepton->Fill(Ecomplet_lepton);
        h_trackRes_lepton->Fill(Reco_LengthRes);
        h_muonwtrk_length->Fill(truth_lengthLepton);
      }
    }
    if (MCproton_reco && MCproton) {
      if (MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE)) {
        h_Pproton_num->Fill(MC_leading_ProtonP);
        h_Efrac_proton->Fill(Efrac_proton);
        h_Ecomplet_proton->Fill(Ecomplet_proton);
        h_trackRes_proton->Fill(Reco_LengthResProton);
        h_protonwtrk_length->Fill(proton_length);
      }
    }
    if (MCpion_plus_reco && MCpion_plus) {
      if (MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE)) {
        h_Ppion_plus_num->Fill(MC_leading_PionPlusP);
        h_Efrac_pion_plus->Fill(Efrac_pionplus);
        h_Ecomplet_pion_plus->Fill(Ecomplet_pionplus);
        h_trackRes_pion_plus->Fill(Reco_LengthResPionPlus);
        h_pionpwtrk_length->Fill(pion_plus_length);
      }
    }
    if (MCpion_minus_reco && MCpion_minus) {
      if (MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE)) {
        h_Ppion_minus_num->Fill(MC_leading_PionMinusP);
        h_Efrac_pion_minus->Fill(Efrac_pionminus);
        h_Ecomplet_pion_minus->Fill(Ecomplet_pionminus);
        h_trackRes_pion_minus->Fill(Reco_LengthResPionMinus);
        h_pionmwtrk_length->Fill(pion_minus_length);
      }
    }
    if (MCkaon_reco && MCkaon) {
      if (MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE)) {
        h_Pkaon_num->Fill(MC_kaonP);
        h_Efrac_kaon->Fill(Efrac_kaon);
        h_Ecomplet_kaon->Fill(Ecomplet_kaon);
        h_trackRes_kaon->Fill(kaonLength - trackLength_kaon);
        h_kaonwtrk_length->Fill(kaonLength);
      }
    }
    //Non neutrino events
    //=========================================================
    if (!fisNeutrinoInt) {
      if (MClepton_reco && MClepton) {
        h_Pmu_num->Fill(MC_leptonP);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_Ecomplet_lepton->Fill(Ecomplet_lepton);
        h_trackRes_lepton->Fill(Reco_LengthRes);
        h_muonwtrk_length->Fill(truth_lengthLepton);
      }
      if (MCkaon_reco && MCkaon) {
        h_Pkaon_num->Fill(MC_kaonP);
        h_Efrac_kaon->Fill(Efrac_kaon);
        h_Ecomplet_kaon->Fill(Ecomplet_kaon);
        h_trackRes_kaon->Fill(kaonLength - trackLength_kaon);
        h_kaonwtrk_length->Fill(kaonLength);
      }
      if (MCproton_reco && MCproton) {
        h_Pproton_num->Fill(MC_leading_ProtonP);
        h_Efrac_proton->Fill(Efrac_proton);
        h_Ecomplet_proton->Fill(Ecomplet_proton);
        h_trackRes_proton->Fill(Reco_LengthResProton);
        h_protonwtrk_length->Fill(proton_length);
      }
      if (MCpion_plus_reco && MCpion_plus) {
        h_Ppion_plus_num->Fill(MC_leading_PionPlusP);
        h_Efrac_pion_plus->Fill(Efrac_pionplus);
        h_Ecomplet_pion_plus->Fill(Ecomplet_pionplus);
        h_trackRes_pion_plus->Fill(Reco_LengthResPionPlus);
        h_pionpwtrk_length->Fill(pion_plus_length);
      }
      if (MCpion_minus_reco && MCpion_minus) {
        h_Ppion_minus_num->Fill(MC_leading_PionMinusP);
        h_Efrac_pion_minus->Fill(Efrac_pionminus);
        h_Ecomplet_pion_minus->Fill(Ecomplet_pionminus);
        h_trackRes_pion_minus->Fill(Reco_LengthResPionMinus);
        h_pionmwtrk_length->Fill(pion_minus_length);
      }
      if (MCmichel_reco && MCmichel) {
        h_Pmichel_e_num->Fill(MC_michelP);
        h_Efrac_michel->Fill(Efrac_michel);
        h_Ecomplet_michel->Fill(Ecomplet_michel);
        h_trackRes_michel->Fill(michelLength - trackLength_michel);
        h_michelwtrk_length->Fill(michelLength);
      }
    }
  }
  //========================================================================
  void
  NeutrinoTrackingEff::truthMatcher(detinfo::DetectorClocksData const& clockData,
                                    std::vector<art::Ptr<recob::Hit>> all_hits,
                                    std::vector<art::Ptr<recob::Hit>> track_hits,
                                    const simb::MCParticle*& MCparticle,
                                    double& Efrac,
                                    double& Ecomplet)
  {
    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
    std::map<int, double> trkID_E;
    for (size_t j = 0; j < track_hits.size(); ++j) {
      art::Ptr<recob::Hit> hit = track_hits[j];
      std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
      for (size_t k = 0; k < TrackIDs.size(); k++) {
        trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
      }
    }
    double E_em = 0.0;
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =
      0.0; // amount of energy deposited by the particle that deposited more energy...

    // If the collection of hits have more than one particle associate
    // save the particle w/ the highest energy deposition since we are
    // looking for muons/pions/protons this should be enough
    if (!trkID_E.size()) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for (std::map<int, double>::iterator ii = trkID_E.begin(); ii != trkID_E.end(); ++ii) {
      total_E += ii->second;
      if ((ii->second) > max_E) {
        partial_E = ii->second;
        max_E = ii->second;
        TrackID = ii->first;
        if (TrackID < 0) E_em += ii->second;
      }
    }
    MCparticle = pi_serv->TrackIdToParticle_P(TrackID);

    // In the current simulation, we do not save EM Shower daughters
    // in GEANT. But we do save the energy deposition in TrackIDEs. If
    // the energy deposition is from a particle that is the daughter
    // of an EM particle, the negative of the parent track ID is saved
    // in TrackIDE for the daughter particle we don't want to track
    // gammas or any other EM activity
    if (TrackID < 0) return;

    Efrac = (partial_E) / total_E;

    // Completeness
    double totenergy = 0;
    for (size_t k = 0; k < all_hits.size(); ++k) {
      art::Ptr<recob::Hit> hit = all_hits[k];
      std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
      for (size_t l = 0; l < TrackIDs.size(); ++l) {
        if (TrackIDs[l].trackID == TrackID) totenergy += TrackIDs[l].energy;
      }
    }
    Ecomplet = partial_E / totenergy;
  }
  //========================================================================
  double
  NeutrinoTrackingEff::truthLength(const detinfo::DetectorClocksData& clockData,
                                   const detinfo::DetectorPropertiesData& detProp,
                                   const simb::MCParticle* MCparticle)
  {
    // Calculate the truth length considering only the part that is
    // inside the TPC Base on a peace of code from
    // dune/TrackingAna/TrackingEfficiency_module.cc

    if (!MCparticle) return -999.0;
    int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
    std::vector<double> TPCLengthHits(numberTrajectoryPoints, 0.0);
    int FirstHit = 0, LastHit = 0;
    double TPCLength = 0.0;
    bool BeenInVolume = false;

    double const WindowSize = detProp.NumberTimeSamples() * clockData.TPCClock().TickPeriod() * 1e3;

    for (unsigned int MCHit = 0; MCHit < TPCLengthHits.size(); ++MCHit) {
      const TLorentzVector& tmpPosition = MCparticle->Position(MCHit);
      double const tmpPosArray[] = {tmpPosition[0], tmpPosition[1], tmpPosition[2]};
      if (MCHit != 0)
        TPCLengthHits[MCHit] = sqrt(pow((MCparticle->Vx(MCHit - 1) - MCparticle->Vx(MCHit)), 2) +
                                    pow((MCparticle->Vy(MCHit - 1) - MCparticle->Vy(MCHit)), 2) +
                                    pow((MCparticle->Vz(MCHit - 1) - MCparticle->Vz(MCHit)), 2));
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if (tpcid.isValid) {
        // -- Check if hit is within drift window...
        geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
        geo::TPCGeo const& tpc = cryo.TPC(tpcid.TPC);
        double XPlanePosition = tpc.PlaneLocation(0)[0];
        double DriftTimeCorrection = fabs(tmpPosition[0] - XPlanePosition) / fDriftVelocity;
        double TimeAtPlane = MCparticle->T() + DriftTimeCorrection;

        if (TimeAtPlane < trigger_offset(clockData) ||
            TimeAtPlane > trigger_offset(clockData) + WindowSize)
          continue;
        LastHit = MCHit;
        if (!BeenInVolume) {
          BeenInVolume = true;
          FirstHit = MCHit;
        }
      }
    }
    for (int Hit = FirstHit + 1; Hit <= LastHit; ++Hit)
      TPCLength += TPCLengthHits[Hit];
    return TPCLength;
  }
  //========================================================================
  bool
  NeutrinoTrackingEff::insideFV(double vertex[4])
  {
    double const x = vertex[0];
    double const y = vertex[1];
    double const z = vertex[2];

    return x > fFidVolXmin && x < fFidVolXmax && y > fFidVolYmin && y < fFidVolYmax &&
           z > fFidVolZmin && z < fFidVolZmax;
  }
  //========================================================================
  void
  NeutrinoTrackingEff::doEfficiencies()
  {
    art::ServiceHandle<art::TFileService const> tfs;

    if (TEfficiency::CheckConsistency(*h_Ev_num, *h_Ev_den)) {
      h_Eff_Ev = tfs->make<TEfficiency>(*h_Ev_num, *h_Ev_den);
      TGraphAsymmErrors* grEff_Ev = h_Eff_Ev->CreateGraph();
      grEff_Ev->Write("grEff_Ev");
      h_Eff_Ev->Write("h_Eff_Ev");
    }
    if (TEfficiency::CheckConsistency(*h_Pmu_num, *h_Pmu_den)) {
      h_Eff_Pmu = tfs->make<TEfficiency>(*h_Pmu_num, *h_Pmu_den);
      TGraphAsymmErrors* grEff_Pmu = h_Eff_Pmu->CreateGraph();
      grEff_Pmu->Write("grEff_Pmu");
      h_Eff_Pmu->Write("h_Eff_Pmu");
    }
    if (TEfficiency::CheckConsistency(*h_theta_num, *h_theta_den)) {
      h_Eff_theta = tfs->make<TEfficiency>(*h_theta_num, *h_theta_den);
      TGraphAsymmErrors* grEff_theta = h_Eff_theta->CreateGraph();
      grEff_theta->Write("grEff_theta");
      h_Eff_theta->Write("h_Eff_theta");
    }
    if (TEfficiency::CheckConsistency(*h_Pproton_num, *h_Pproton_den)) {
      h_Eff_Pproton = tfs->make<TEfficiency>(*h_Pproton_num, *h_Pproton_den);
      TGraphAsymmErrors* grEff_Pproton = h_Eff_Pproton->CreateGraph();
      grEff_Pproton->Write("grEff_Pproton");
      h_Eff_Pproton->Write("h_Eff_Pproton");
    }
    if (TEfficiency::CheckConsistency(*h_Ppion_plus_num, *h_Ppion_plus_den)) {
      h_Eff_Ppion_plus = tfs->make<TEfficiency>(*h_Ppion_plus_num, *h_Ppion_plus_den);
      TGraphAsymmErrors* grEff_Ppion_plus = h_Eff_Ppion_plus->CreateGraph();
      grEff_Ppion_plus->Write("grEff_Ppion_plus");
      h_Eff_Ppion_plus->Write("h_Eff_Ppion_plus");
    }
    if (TEfficiency::CheckConsistency(*h_Ppion_minus_num, *h_Ppion_minus_den)) {
      h_Eff_Ppion_minus = tfs->make<TEfficiency>(*h_Ppion_minus_num, *h_Ppion_minus_den);
      TGraphAsymmErrors* grEff_Ppion_minus = h_Eff_Ppion_minus->CreateGraph();
      grEff_Ppion_minus->Write("grEff_Ppion_minus");
      h_Eff_Ppion_minus->Write("h_Eff_Ppion_minus");
    }
    if (TEfficiency::CheckConsistency(*h_muonwtrk_length, *h_muon_length)) {
      h_Eff_Lmuon = tfs->make<TEfficiency>(*h_muonwtrk_length, *h_muon_length);
      TGraphAsymmErrors* grEff_Lmuon = h_Eff_Lmuon->CreateGraph();
      grEff_Lmuon->Write("grEff_Lmuon");
      h_Eff_Lmuon->Write("h_Eff_Lmuon");
    }
    if (TEfficiency::CheckConsistency(*h_protonwtrk_length, *h_proton_length)) {
      h_Eff_Lproton = tfs->make<TEfficiency>(*h_protonwtrk_length, *h_proton_length);
      TGraphAsymmErrors* grEff_Lproton = h_Eff_Lproton->CreateGraph();
      grEff_Lproton->Write("grEff_Lproton");
      h_Eff_Lproton->Write("h_Eff_Lproton");
    }
    if (TEfficiency::CheckConsistency(*h_pionpwtrk_length, *h_pionp_length)) {
      h_Eff_Lpion_plus = tfs->make<TEfficiency>(*h_pionpwtrk_length, *h_pionp_length);
      TGraphAsymmErrors* grEff_Lpion_plus = h_Eff_Lpion_plus->CreateGraph();
      grEff_Lpion_plus->Write("grEff_Lpion_plus");
      h_Eff_Lpion_plus->Write("h_Eff_Lpion_plus");
    }
    if (TEfficiency::CheckConsistency(*h_pionpwtrk_length, *h_pionp_length)) {
      h_Eff_Lpion_minus = tfs->make<TEfficiency>(*h_pionmwtrk_length, *h_pionm_length);
      TGraphAsymmErrors* grEff_Lpion_minus = h_Eff_Lpion_minus->CreateGraph();
      grEff_Lpion_minus->Write("grEff_Lpion_minus");
      h_Eff_Lpion_minus->Write("h_Eff_Lpion_minus");
    }
    if (TEfficiency::CheckConsistency(*h_Pkaon_num, *h_Pkaon_den)) {
      h_Eff_Pkaon = tfs->make<TEfficiency>(*h_Pkaon_num, *h_Pkaon_den);
      TGraphAsymmErrors* grEff_Pkaon = h_Eff_Pkaon->CreateGraph();
      grEff_Pkaon->Write("grEff_Pkaon");
      h_Eff_Pkaon->Write("h_Eff_Pkaon");
    }
    if (TEfficiency::CheckConsistency(*h_kaonwtrk_length, *h_kaon_length)) {
      h_Eff_Lkaon = tfs->make<TEfficiency>(*h_kaonwtrk_length, *h_kaon_length);
      TGraphAsymmErrors* grEff_Lkaon = h_Eff_Lkaon->CreateGraph();
      grEff_Lkaon->Write("grEff_Lkaon");
      h_Eff_Lkaon->Write("h_Eff_Lkaon");
    }
    if (TEfficiency::CheckConsistency(*h_Pmichel_e_num, *h_Pmichel_e_den)) {
      h_Eff_Pmichel = tfs->make<TEfficiency>(*h_Pmichel_e_num, *h_Pmichel_e_den);
      TGraphAsymmErrors* grEff_Pmichel = h_Eff_Pmichel->CreateGraph();
      grEff_Pmichel->Write("grEff_Pmichel");
      h_Eff_Pmichel->Write("h_Eff_Pmichel");
    }
    if (TEfficiency::CheckConsistency(*h_michelwtrk_length, *h_michel_length)) {
      h_Eff_Lmichel = tfs->make<TEfficiency>(*h_michelwtrk_length, *h_michel_length);
      TGraphAsymmErrors* grEff_Lmichel = h_Eff_Lmichel->CreateGraph();
      grEff_Lmichel->Write("grEff_Lmichel");
      h_Eff_Lmichel->Write("h_Eff_Lmichel");
    }
  }
  //========================================================================
  DEFINE_ART_MODULE(NeutrinoTrackingEff)

}
