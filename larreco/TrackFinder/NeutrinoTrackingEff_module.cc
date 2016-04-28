#ifndef NeutrinoTrackingEff_Module
#define NeutrinoTrackingEff_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Track.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/RawData/ExternalTrigger.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#define MAX_TRACKS 1000
using namespace std;

//========================================================================

namespace DUNE{

class NeutrinoTrackingEff : public art::EDAnalyzer {
public:

    explicit NeutrinoTrackingEff(fhicl::ParameterSet const& pset);
    virtual ~NeutrinoTrackingEff();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

    void processEff(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac);
    double truthLength( const simb::MCParticle *MCparticle );
    bool insideFV(double vertex[4]);
    void doEfficiencies();
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    int         fNeutrinoPDGcode;
    int		fLeptonPDGcode;
    double      fMaxNeutrinoE;
    double      fMaxLeptonP;
    bool	fSaveMCTree; 
    bool    	fisNeutrinoInt;
 
    TTree *fEventTree;
    TTree *fHitsTree;

    TH1D *h_Ev_den;
    TH1D *h_Ev_num;
    TH1D *h_Pmu_den;
    TH1D *h_Pmu_num;
    TH1D *h_theta_den;
    TH1D *h_theta_num;
    TH1D *h_Pproton_den;
    TH1D *h_Pproton_num;
    TH1D *h_Ppion_plus_den; 
    TH1D *h_Ppion_plus_num; 
    TH1D *h_Ppion_minus_den; 
    TH1D *h_Ppion_minus_num; 

    TH1D *h_Efrac_lepton;     
    TH1D *h_Efrac_proton;     
    TH1D *h_Efrac_pion_plus;     
    TH1D *h_Efrac_pion_minus;     
    TH1D *h_trackRes_lepton;
    TH1D *h_trackRes_proton;
    TH1D *h_trackRes_pion_plus;
    TH1D *h_trackRes_pion_minus;

    TEfficiency* h_Eff_Ev = 0;
    TEfficiency* h_Eff_Pmu = 0;
    TEfficiency* h_Eff_theta = 0;
    TEfficiency* h_Eff_Pproton = 0;
    TEfficiency* h_Eff_Ppion_plus = 0;
    TEfficiency* h_Eff_Ppion_minus = 0;

    //nucleon decay histograms
    TH1D *h_Pkaon_den;
    TH1D *h_Pkaon_num; 
    TH1D *h_Pmichel_e_den;
    TH1D *h_Pmichel_e_num;
    TH1D *h_Efrac_kaon;
    TH1D *h_trackRes_kaon; 
    TH1D *h_Efrac_michel;
    TH1D *h_trackRes_michel;
    TEfficiency* h_Eff_Pkaon =0;
    TEfficiency* h_Eff_Pmichel =0;


    // Event 
    int Event;
    int Run;
    int SubRun;

    //MC truth
    int    MC_incoming_PDG;
    int    MC_lepton_PDG;
    int    MC_isCC;
    int    MC_channel;
    int    MC_target;
    double MC_Q2;
    double MC_W;
    double MC_vertex[4];
    double MC_incoming_P[4];
    double MC_lepton_startMomentum[4]; 
    double MC_lepton_endMomentum[4]; 
    double MC_lepton_startXYZT[4];
    double MC_lepton_endXYZT[4];
    double MC_lepton_theta; 
    double MC_leptonP;
    int    MC_leptonID;
    int    MC_kaonID;
    double MC_kaonP;
    double MC_kaon_startMomentum[4]; 
    double MC_kaon_endMomentum[4]; 
    double MC_kaon_startXYZT[4];
    double MC_kaon_endXYZT[4];
    int    MC_kaonTrack;
    int    MC_michelID;
    double MC_michelP;
    double MC_michel_startMomentum[4]; 
    double MC_michel_endMomentum[4]; 
    double MC_michel_startXYZT[4];
    double MC_michel_endXYZT[4];
    int    MC_michelTrack;
    int    MC_leading_protonID;
    double MC_leading_ProtonP;
    int    MC_leading_PionPlusID;
    int    MC_leading_PionMinusID;
    double MC_leading_PionPlusP;
    double MC_leading_PionMinusP;
    int    MC_LeptonTrack;
    int    MC_ProtonTrack;
    int    MC_PionPlusTrack; 
    int    MC_PionMinusTrack;
 
    int    MC_Ntrack;  
    int    MC_id[MAX_TRACKS];  
    int    MC_pdg[MAX_TRACKS]; 
    int    MC_mother[MAX_TRACKS];  
    double MC_startXYZT[MAX_TRACKS][4]; 
    double MC_endXYZT[MAX_TRACKS][4];  
    double MC_startMomentum[MAX_TRACKS][4]; 
    double MC_endMomentum[MAX_TRACKS][4];  

    double Reco_EfracLepton;
    double Reco_LengthRes; 
    double Reco_LengthResProton;
    double Reco_LengthResPionPlus;
    double Reco_LengthResPionMinus;
    double Reco_LengthReskaon;
    double Reco_Efrac_kaon;
    double Reco_Efrac_michel;
    double Reco_LengthResmichel;
    int    n_recoTrack;

    float fFidVolCutX;
    float fFidVolCutY;
    float fFidVolCutZ;

    float fFidVolXmin;
    float fFidVolXmax;
    float fFidVolYmin;
    float fFidVolYmax;
    float fFidVolZmin;
    float fFidVolZmax;

    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;
    art::ServiceHandle<geo::Geometry> geom;

 }; // class NeutrinoTrackingEff


//========================================================================
NeutrinoTrackingEff::NeutrinoTrackingEff(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
NeutrinoTrackingEff::~NeutrinoTrackingEff(){
  //destructor
}
//========================================================================
void NeutrinoTrackingEff::reconfigure(fhicl::ParameterSet const& p){

    fMCTruthModuleLabel  = p.get<std::string>("MCTruthModuleLabel");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fisNeutrinoInt	 = p.get<bool>("isNeutrinoInt");
    fLeptonPDGcode       = p.get<int>("LeptonPDGcode");
    fNeutrinoPDGcode     = p.get<int>("NeutrinoPDGcode");
    fMaxNeutrinoE	 = p.get<double>("MaxNeutrinoE");
    fMaxLeptonP          = p.get<double>("MaxLeptonP");
    fSaveMCTree		 = p.get<bool>("SaveMCTree");
    fFidVolCutX          = p.get<float>("FidVolCutX");
    fFidVolCutY          = p.get<float>("FidVolCutY");
    fFidVolCutZ          = p.get<float>("FidVolCutZ");
}
//========================================================================
void NeutrinoTrackingEff::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  // Get geometry.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  std::cout<<"Fiducial volume:"<<"\n"
	   <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
	   <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
	   <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";

  art::ServiceHandle<art::TFileService> tfs;

  double E_bins[21] ={0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,14.0,17.0,20.0,25.0};
  double theta_bin[44]= { 0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,55.,60.,65.,70.,75.,80.,85.,90.};
  double Pbins[18] ={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0};
  
  for (int i = 0; i<21; ++i) E_bins[i] *= fMaxNeutrinoE/25.;
  for (int i = 0; i<18; ++i) Pbins[i] *= fMaxLeptonP/3.0;

  h_Ev_den = tfs->make<TH1D>("h_Ev_den","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
  h_Ev_num = tfs->make<TH1D>("h_Ev_num","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
  h_Pmu_den = tfs->make<TH1D>("h_Pmu_den","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,E_bins);
  h_Pmu_num = tfs->make<TH1D>("h_Pmu_num","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,E_bins);
  h_theta_den = tfs->make<TH1D>("h_theta_den","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
  h_theta_num = tfs->make<TH1D>("h_theta_num","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
  h_Pproton_den = tfs->make<TH1D>("h_Pproton_den","Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Pproton_num = tfs->make<TH1D>("h_Pproton_num","Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
  h_Ppion_plus_den = tfs->make<TH1D>("h_Ppion_plus_den", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ppion_plus_num = tfs->make<TH1D>("h_Ppion_plus_num", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ppion_minus_den = tfs->make<TH1D>("h_Ppion_minus_den", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
  h_Ppion_minus_num = tfs->make<TH1D>("h_Ppion_minus_num", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
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

  h_Efrac_lepton = tfs->make<TH1D>("h_Efrac_lepton","Efrac Lepton; Track Energy fraction;",60,0,1.2);
  h_Efrac_proton = tfs->make<TH1D>("h_Efrac_proton","Efrac Proton; Track Energy fraction;",60,0,1.2);
  h_Efrac_pion_plus = tfs->make<TH1D>("h_Efrac_pion_plus","Efrac Pion +; Track Energy fraction;",60,0,1.2);
  h_Efrac_pion_minus = tfs->make<TH1D>("h_Efrac_pion_minus","Efrac Pion -; Track Energy fraction;",60,0,1.2);
  h_trackRes_lepton = tfs->make<TH1D>("h_trackRes_lepton", "Muon Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_proton = tfs->make<TH1D>("h_trackRes_proton", "Proton Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_pion_plus = tfs->make<TH1D>("h_trackRes_pion_plus", "Pion + Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_pion_minus = tfs->make<TH1D>("h_trackRes_pion_minus", "Pion - Residual; Truth length - Reco length (cm);",200,-100,100);
  h_Efrac_lepton->Sumw2();
  h_Efrac_proton->Sumw2();
  h_Efrac_pion_plus->Sumw2();
  h_Efrac_pion_minus->Sumw2();
  h_trackRes_lepton->Sumw2();
  h_trackRes_proton->Sumw2();
  h_trackRes_pion_plus->Sumw2();
  h_trackRes_pion_minus->Sumw2();

  if(!fisNeutrinoInt ){
    double Pbins[21] ={0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0};
    h_Pmu_den = tfs->make<TH1D>("h_lepton_den","Lepton Momentum; Lepton Momentum (GeV); Tracking Efficiency",20,Pbins);
    h_Pmu_num = tfs->make<TH1D>("h_lepton_num","Lepton Momentum; Lepton Momentum (GeV); Tracking Efficiency",20,Pbins);
    h_Pkaon_den = tfs->make<TH1D>("h_Pkaon_den","Kaon; Kaon Momentum (GeV); Tracking Efficiency", 20, Pbins);
    h_Pkaon_num = tfs->make<TH1D>("h_Pkaon_num","Kaon; Kaon Momentum (GeV); Tracking Efficiency", 20, Pbins);
    h_Pmichel_e_den = tfs->make<TH1D>("h_Pmichel_e_den","Michel Electron; Michele e Momentum (GeV); Tracking Efficiency", 20, Pbins);
    h_Pmichel_e_num = tfs->make<TH1D>("h_Pmichel_e_num","Michel Electron; Michele e Momentum (GeV); Tracking Efficiency", 20, Pbins);
    //h_Pmu_den->Sumw2();
    //h_Pmu_num->Sumw2();
    h_Pkaon_den->Sumw2();
    h_Pkaon_num->Sumw2();
    h_Pmichel_e_den->Sumw2(); 
    h_Pmichel_e_num->Sumw2(); 
    h_Efrac_kaon = tfs->make<TH1D>("h_Efrac_kaon","Efrac Kaon; Track Energy fraction;",60,0,1.2);
    h_trackRes_kaon = tfs->make<TH1D>("h_trackRes_kaon","Kaon Residual; Truth length - Reco length (cm);",200,-100,100);
    h_Efrac_michel = tfs->make<TH1D>("h_Efrac_michel","Efrac Michel; Track Energy fraction;",60,0,1.2);
    h_trackRes_michel = tfs->make<TH1D>("h_trackRes_michel","Michel Residual; Truth length - Reco length (cm);",200,-100,100);
    h_Efrac_kaon->Sumw2();
    h_trackRes_kaon->Sumw2();
    h_Efrac_michel->Sumw2();
    h_trackRes_michel->Sumw2();
  }


  if( fSaveMCTree ){
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");
    fEventTree->Branch("eventNo", &Event);
    fEventTree->Branch("runNo", &Run);
    fEventTree->Branch("subRunNo", &SubRun);
    fEventTree->Branch("mc_incoming_PDG", &MC_incoming_PDG);
    fEventTree->Branch("mc_lepton_PDG", &MC_lepton_PDG);
    fEventTree->Branch("mc_isCC", &MC_isCC);
    fEventTree->Branch("mc_target", &MC_target);
    fEventTree->Branch("mc_channel", &MC_channel);
    fEventTree->Branch("mc_Q2", &MC_Q2);
    fEventTree->Branch("mc_W", &MC_W);
    fEventTree->Branch("mc_vertex", &MC_vertex, "mc_vertex[4]/D");
    fEventTree->Branch("mc_incoming_P", &MC_incoming_P, "mc_incoming_P[4]/D");
    fEventTree->Branch("mc_lepton_startMomentum", &MC_lepton_startMomentum, "mc_lepton_startMomentum[4]/D");
    fEventTree->Branch("mc_lepton_endMomentum", &MC_lepton_endMomentum, "mc_lepton_endMomentum[4]/D");
    fEventTree->Branch("mc_lepton_startXYZT", &MC_lepton_startXYZT, "mc_lepton_startXYZT[4]/D");
    fEventTree->Branch("mc_lepton_endXYZT", &MC_lepton_endXYZT, "mc_lepton_endXYZT[4]/D");
    fEventTree->Branch("mc_lepton_theta", &MC_lepton_theta, "mc_lepton_theta/D");
    fEventTree->Branch("mc_leptonID", &MC_leptonID);
    fEventTree->Branch("mc_leptonP", &MC_leptonP, "mc_leptonP/D");
    fEventTree->Branch("mc_kaonID", &MC_kaonID);
    fEventTree->Branch("mc_kaonTrack", &MC_kaonTrack);
    fEventTree->Branch("mc_michelID", &MC_michelID);
    fEventTree->Branch("mc_michelTrack", &MC_michelTrack);
    fEventTree->Branch("mc_leadingProtonID", &MC_leading_protonID);
    fEventTree->Branch("mc_leadingProtonID", &MC_leading_ProtonP);
    fEventTree->Branch("mc_leadingPionPlusID", &MC_leading_PionPlusID);
    fEventTree->Branch("mc_leadingPionPlusP", &MC_leading_PionPlusP);
    fEventTree->Branch("mc_leadingPionMinusID", &MC_leading_PionMinusID);
    fEventTree->Branch("mc_leadingPionMinusP", &MC_leading_PionMinusP);
    fEventTree->Branch("mc_kaon_startMomentum", &MC_kaon_startMomentum, "mc_kaon_startMomentum[4]/D");
    fEventTree->Branch("mc_kaon_endMomentum", &MC_kaon_endMomentum, "mc_kaon_endMomentum[4]/D");
    fEventTree->Branch("mc_kaon_startXYZT", &MC_kaon_startXYZT, "mc_kaon_startXYZT[4]/D");
    fEventTree->Branch("mc_kaon_endXYZT", &MC_kaon_endXYZT, "mc_kaon_endXYZT[4]/D");
    fEventTree->Branch("mc_michel_startMomentum", &MC_michel_startMomentum, "mc_michel_startMomentum[4]/D");
    fEventTree->Branch("mc_michel_endMomentum", &MC_michel_endMomentum, "mc_michel_endMomentum[4]/D");
    fEventTree->Branch("mc_michel_startXYZT", &MC_michel_startXYZT, "mc_michel_startXYZT[4]/D");
    fEventTree->Branch("mc_leptonTrack", &MC_LeptonTrack);
    fEventTree->Branch("mc_protonTrack", &MC_ProtonTrack);
    fEventTree->Branch("mc_kaonTrack", &MC_kaonTrack);
    fEventTree->Branch("mc_michelTrack", &MC_michelTrack);
    fEventTree->Branch("mc_pionPlusTrack", &MC_PionPlusTrack);
    fEventTree->Branch("mc_pionMinusTrack", &MC_PionMinusTrack);
    fEventTree->Branch("mc_Ntrack", &MC_Ntrack);  // number of particles 
    fEventTree->Branch("mc_id", &MC_id, "mc_id[mc_Ntrack]/I");  
    fEventTree->Branch("mc_pdg", &MC_pdg, "mc_pdg[mc_Ntrack]/I"); 
    fEventTree->Branch("mc_mother", &MC_mother, "mc_mother[mc_Ntrack]/I"); 
    fEventTree->Branch("mc_startXYZT", &MC_startXYZT, "mc_startXYZT[mc_Ntrack][4]/D");  
    fEventTree->Branch("mc_endXYZT", &MC_endXYZT, "mc_endXYZT[mc_Ntrack][4]/D"); 
    fEventTree->Branch("mc_startMomentum", &MC_startMomentum, "mc_startMomentum[mc_Ntrack][4]/D");  
    fEventTree->Branch("mc_endMomentum", &MC_endMomentum, "mc_endMomentum[mc_Ntrack][4]/D"); 
    
    fEventTree->Branch("reco_efracLepton", &Reco_EfracLepton);
    fEventTree->Branch("reco_lengthResProton", &Reco_LengthResProton);
    fEventTree->Branch("reco_lengthResPionPlus", &Reco_LengthResPionPlus);
    fEventTree->Branch("reco_lengthResPionMinus", &Reco_LengthResPionMinus);
    fEventTree->Branch("reco_MC_LeptonTrack", &MC_LeptonTrack);
    fEventTree->Branch("reco_MC_ProtonTrack", &MC_ProtonTrack);
    fEventTree->Branch("reco_MC_PionPlusTrack", &MC_PionPlusTrack);
    fEventTree->Branch("reco_MC_PionMinusTrack", &MC_PionMinusTrack);
    fEventTree->Branch("reco_MC_kaonTrack", &MC_kaonTrack);
    fEventTree->Branch("reco_MC_michelTrack", &MC_michelTrack);
    fEventTree->Branch("reco_efrackaon", &Reco_Efrac_kaon);
    fEventTree->Branch("reco_efracmichel", &Reco_Efrac_michel);
    
  }


}
//========================================================================
void NeutrinoTrackingEff::endJob(){
     
  doEfficiencies();

}
//========================================================================
void NeutrinoTrackingEff::beginRun(const art::Run& /*run*/){
  mf::LogInfo("NeutrinoTrackingEff")<<"begin run..."<<std::endl;
}
//========================================================================
void NeutrinoTrackingEff::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    reset();

    Event  = event.id().event(); 
    Run    = event.run();
    SubRun = event.subRun();
    bool isFiducial = false;
    processEff(event, isFiducial);
    if( fSaveMCTree ){
      if(isFiducial && MC_isCC == 1) fEventTree->Fill();
    }
}
//========================================================================
void NeutrinoTrackingEff::processEff( const art::Event& event, bool &isFiducial){

    //!save neutrino's interaction info 
    art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
    event.getByLabel(fMCTruthModuleLabel, MCtruthHandle);
    std::vector<art::Ptr<simb::MCTruth>> MCtruthlist;
    art::fill_ptr_vector(MCtruthlist, MCtruthHandle);
    art::Ptr<simb::MCTruth> MCtruth; 
    //For now assume that there is only one neutrino interaction...
    int MCinteractions = MCtruthlist.size();
    for( int i =0; i<MCinteractions; i++){
       MCtruth = MCtruthlist[i];
       if( MCtruth->NeutrinoSet() ){
         simb::MCNeutrino nu = MCtruth->GetNeutrino();
         if( nu.CCNC() == 0 ) MC_isCC = 1;
         else if ( nu.CCNC() == 1 ) MC_isCC = 0; 
         simb::MCParticle neutrino = nu.Nu();
         MC_target = nu.Target();
         MC_incoming_PDG = nu.Nu().PdgCode();
         MC_Q2 = nu.QSqr();
         MC_channel = nu.InteractionType();
         MC_W = nu.W();
         const TLorentzVector& nu_momentum = nu.Nu().Momentum(0); 
         nu_momentum.GetXYZT(MC_incoming_P); 
         const TLorentzVector& vertex =neutrino.Position(0); 
         vertex.GetXYZT(MC_vertex);
         simb::MCParticle lepton = nu.Lepton();
         MC_lepton_PDG = lepton.PdgCode();
          LOG_DEBUG("NeutrinoTrackingEff")<<"Incoming E "<<MC_incoming_P[3]<<" is CC? "<<MC_isCC<<" nuPDG "<<MC_incoming_PDG<<" target "<<MC_target<<" vtx "<<MC_vertex[0]<<" "<<MC_vertex[1]<<" "<<MC_vertex[2]<<" "<<MC_vertex[3];
       }
    }

    //!save FS particles
    double tmp_leadingPionPlusE = 0.0;
    double tmp_leadingPionMinusE = 0.0;
    double tmp_leadingProtonE  = 0.0; 
  
    simb::MCParticle *MClepton = NULL; 
    simb::MCParticle *MCproton = NULL;
    simb::MCParticle *MCpion_plus = NULL;
    simb::MCParticle *MCpion_minus = NULL;
    simb::MCParticle *MCkaon = NULL;
    simb::MCParticle *MCmichel = NULL;
 
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index
    MC_Ntrack = plist.size();

    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
       particle = ipar->second;
       if( particle->PdgCode() == fLeptonPDGcode && particle->Mother() == 0 ){  //primary lepton
         const TLorentzVector& lepton_momentum =particle->Momentum(0); 
         const TLorentzVector& lepton_position =particle->Position(0); 
         const TLorentzVector& lepton_positionEnd   = particle->EndPosition();
         const TLorentzVector& lepton_momentumEnd   = particle->EndMomentum();
         lepton_momentum.GetXYZT(MC_lepton_startMomentum);
         lepton_position.GetXYZT(MC_lepton_startXYZT);
         lepton_positionEnd.GetXYZT(MC_lepton_endXYZT);
         lepton_momentumEnd.GetXYZT(MC_lepton_endMomentum);
         MC_leptonID = particle->TrackId();
         MC_leptonP = sqrt(pow(MC_lepton_startMomentum[0],2)+pow(MC_lepton_startMomentum[1],2)+pow(MC_lepton_startMomentum[2],2));
         MClepton = particle;
       }
       if( particle->Mother() == 0 ){   //save primary particle i.e. from the neutrino interaction
         MC_id[i] = particle->TrackId();
         MC_pdg[i] = particle->PdgCode();
         MC_mother[i] = particle->Mother();
         const TLorentzVector& positionStart = particle->Position(0);
         const TLorentzVector& positionEnd   = particle->EndPosition();
         const TLorentzVector& momentumStart = particle->Momentum(0);
         const TLorentzVector& momentumEnd   = particle->EndMomentum();
         positionStart.GetXYZT(MC_startXYZT[i]);
         positionEnd.GetXYZT(MC_endXYZT[i]);
         momentumStart.GetXYZT(MC_startMomentum[i]);
         momentumEnd.GetXYZT(MC_endMomentum[i]);
         //save leading pion and proton
         if( particle->PdgCode() == 2212 ){
           if(particle->Momentum().E() > tmp_leadingProtonE){
             tmp_leadingProtonE = particle->Momentum().E();
             MC_leading_protonID = particle->TrackId();          
             MC_leading_ProtonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             //const TLorentzVector& protonEnd   = particle->EndPosition();
             //const TLorentzVector& protonStart = particle->Position();
             //proton_length = sqrt( pow((protonEnd.X()-protonStart.X()),2)+pow((protonEnd.Y()-protonStart.Y()),2)+pow((protonEnd.Z()-protonStart.Z()),2)); 
             MCproton = particle;
           } 
         }
         else if( particle->PdgCode() == 211 ){
           if(particle->Momentum().E() > tmp_leadingPionPlusE){
             tmp_leadingPionPlusE = particle->Momentum().E();
             MC_leading_PionPlusID = particle->TrackId();          
             MC_leading_PionPlusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             //const TLorentzVector& pionEnd = particle->EndPosition();
             //const TLorentzVector& pionStart = particle->Position();
             //pion_plus_length = sqrt( pow((pionEnd.X()-pionStart.X()),2)+pow((pionEnd.Y()-pionStart.Y()),2)+pow((pionEnd.Z()-pionStart.Z()),2));
             MCpion_plus = particle;
           } 
         }
         else if( particle->PdgCode() == -211 ){
           if(particle->Momentum().E() > tmp_leadingPionMinusE){
             tmp_leadingPionMinusE = particle->Momentum().E();
             MC_leading_PionMinusID = particle->TrackId();          
             MC_leading_PionMinusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             //const TLorentzVector& pionEnd = particle->EndPosition();
             //const TLorentzVector& pionStart = particle->Position();
             //pion_minus_length = sqrt( pow((pionEnd.X()-pionStart.X()),2)+pow((pionEnd.Y()-pionStart.Y()),2)+pow((pionEnd.Z()-pionStart.Z()),2));
             MCpion_minus = particle;
           } 
         }
         i++; //paticle index
       }

       //add Nucleon decay stuff
       if(!fisNeutrinoInt ){
         if( particle->Mother() == 0 && particle->PdgCode() == 321 ){   //save primary Kaon
           const TLorentzVector& kaon_momentum =particle->Momentum(0); 
           const TLorentzVector& kaon_position =particle->Position(0); 
           const TLorentzVector& kaon_positionEnd   = particle->EndPosition();
           const TLorentzVector& kaon_momentumEnd   = particle->EndMomentum();
           kaon_momentum.GetXYZT(MC_lepton_startMomentum);
           kaon_position.GetXYZT(MC_lepton_startXYZT);
           kaon_positionEnd.GetXYZT(MC_lepton_endXYZT);
           kaon_momentumEnd.GetXYZT(MC_lepton_endMomentum);
           kaon_position.GetXYZT(MC_vertex); //Primary vertex
           MC_kaonID = particle->TrackId();    
           MC_kaonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
           MCkaon = particle;
         }   
         else if( particle->Mother() ==2 && particle->Process() =="Decay" && particle->PdgCode() == fLeptonPDGcode ){  // lepton
           const TLorentzVector& lepton_momentum =particle->Momentum(0); 
           const TLorentzVector& lepton_position =particle->Position(0); 
           const TLorentzVector& lepton_positionEnd   = particle->EndPosition();
           const TLorentzVector& lepton_momentumEnd   = particle->EndMomentum();
           lepton_momentum.GetXYZT(MC_lepton_startMomentum);
           lepton_position.GetXYZT(MC_lepton_startXYZT);
           lepton_positionEnd.GetXYZT(MC_lepton_endXYZT);
           lepton_momentumEnd.GetXYZT(MC_lepton_endMomentum);
           MC_leptonID = particle->TrackId();
           MC_leptonP = sqrt(pow(MC_lepton_startMomentum[0],2)+pow(MC_lepton_startMomentum[1],2)+pow(MC_lepton_startMomentum[2],2));
           MClepton = particle;
         }
         else if( particle->Process() =="Decay" && particle->PdgCode() == -11){  // michel electron
           const TLorentzVector& michel_momentum =particle->Momentum(0); 
           const TLorentzVector& michel_position =particle->Position(0); 
           const TLorentzVector& michel_positionEnd   = particle->EndPosition();
           const TLorentzVector& michel_momentumEnd   = particle->EndMomentum();
           michel_momentum.GetXYZT(MC_michel_startMomentum);
           michel_position.GetXYZT(MC_michel_startXYZT);
           michel_positionEnd.GetXYZT(MC_michel_endXYZT);
           michel_momentumEnd.GetXYZT(MC_michel_endMomentum);
           MC_michelID = particle->TrackId();
           MC_michelP = sqrt(pow(MC_michel_startMomentum[0],2)+pow(MC_michel_startMomentum[1],2)+pow(MC_michel_startMomentum[2],2));
           MCmichel = particle;
         }

       }
    } 
    //Saving denominator histograms for lepton pions and protons 
    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
    double Pv  = sqrt(pow(MC_incoming_P[0],2)+pow(MC_incoming_P[1],2)+pow(MC_incoming_P[2],2));
    double theta_mu = acos((MC_incoming_P[0]*MC_lepton_startMomentum[0] + MC_incoming_P[1]*MC_lepton_startMomentum[1] +MC_incoming_P[2]*MC_lepton_startMomentum[2])/(Pv*MC_leptonP) );
    theta_mu *= (180.0/3.14159);
    double truth_lengthLepton = truthLength(MClepton); 
    double proton_length = truthLength(MCproton);
    double pion_plus_length = truthLength(MCpion_plus);
    double pion_minus_length = truthLength(MCpion_minus);
    double kaonLength = truthLength(MCkaon);
    double michelLength = truthLength(MCmichel);

    //save CC events within the fiducial volume with the favorite neutrino flavor 
    if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
       if( MClepton ){
         h_Ev_den->Fill(MC_incoming_P[3]);
         h_Pmu_den->Fill(MC_leptonP);
         h_theta_den->Fill(theta_mu);
       }
       if( MCproton ){
         h_Pproton_den->Fill(MC_leading_ProtonP);
       }
       if( MCpion_plus ){
         h_Ppion_plus_den->Fill( MC_leading_PionPlusP);
       }
       if( MCpion_minus ){
         h_Ppion_minus_den->Fill( MC_leading_PionMinusP);
       }
    } 
  
    //save events for Nucleon decay studies
    if(!fisNeutrinoInt ){
      if( MClepton ){
         h_Pmu_den->Fill(MC_leptonP);
       }
       if( MCkaon ){
         h_Pkaon_den->Fill(MC_kaonP);
       }
       if( MCmichel ){
         h_Pmichel_e_den->Fill(MC_michelP);
       }
 
    }
 
    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================
    art::Handle< std::vector<recob::Track> > trackListHandle;
    if(! event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track> > tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);
    n_recoTrack = tracklist.size();

    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    if( n_recoTrack == 0 ){
      LOG_DEBUG("NeutrinoTrackingEff")<<"There are no reco tracks... bye";
      return; 
    }
    LOG_DEBUG("NeutrinoTrackingEff")<<"Found this many reco tracks "<<n_recoTrack;
    
    double Efrac_lepton =0.0;
    double Efrac_proton =0.0;
    double Efrac_pionplus =0.0;
    double Efrac_pionminus =0.0;
    double Efrac_kaon =0.0;
    double Efrac_michel =0.0;
    double trackLength_lepton =0.0;
    double trackLength_proton =0.0;
    double trackLength_pion_plus =0.0;
    double trackLength_pion_minus =0.0;
    double trackLength_kaon =0.0;
    double trackLength_michel =0.0;
    const simb::MCParticle *MClepton_reco = NULL; 
    const simb::MCParticle *MCproton_reco = NULL;
    const simb::MCParticle *MCpion_plus_reco = NULL;
    const simb::MCParticle *MCpion_minus_reco = NULL;
    const simb::MCParticle *MCkaon_reco = NULL;
    const simb::MCParticle *MCmichel_reco = NULL;
    for(int i=0; i<n_recoTrack; i++) {
       art::Ptr<recob::Track> track = tracklist[i];
//       const TVector3 tmp_track_vtx = track->Vertex();
//       double track_vtx[4] ={tmp_track_vtx[0], tmp_track_vtx[1], tmp_track_vtx[2], -999};
//       bool track_isFiducial = insideFV( track_vtx );
//       if( !track_isFiducial ) continue;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       truthMatcher( all_trackHits, particle, tmpEfrac );
       if (!particle) continue;
       //std::cout<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac<<std::endl;
       if(  (particle->PdgCode() == fLeptonPDGcode) && (particle->TrackId() == MC_leptonID) ){
         //save the best track ... based on Efrac if there is more than one track 
         if( tmpEfrac > Efrac_lepton ){
           Efrac_lepton = tmpEfrac;
           trackLength_lepton = track->Length(); 
           MClepton_reco = particle;
         }
       }
       else if( (particle->PdgCode() == 2212) && (particle->TrackId() == MC_leading_protonID) ){
         //save the best track ... based on Efrac if there is more than one track 
         if( tmpEfrac > Efrac_proton ){
           Efrac_proton = tmpEfrac;
           trackLength_proton = track->Length();
           MCproton_reco = particle;
         }
       }
       else if( (particle->PdgCode() == 211) && (particle->TrackId() == MC_leading_PionPlusID) ){
         //save the best track ... based on Efrac if there is more than one track 
         if( tmpEfrac > Efrac_pionplus ){
           Efrac_pionplus = tmpEfrac;
           trackLength_pion_plus = track->Length();
           MCpion_plus_reco = particle;
         }
       }
       else if( (particle->PdgCode() == -211) && (particle->TrackId() == MC_leading_PionMinusID) ){
         //save the best track ... based on Efrac if there is more than one track 
         if( tmpEfrac > Efrac_pionminus ){
           Efrac_pionminus = tmpEfrac;
           trackLength_pion_minus = track->Length();
           MCpion_minus_reco = particle;
         }
       }
       //kaon from nucleon decay
       else if( (particle->PdgCode() == 321) && (particle->TrackId() == MC_kaonID) ){
         if( tmpEfrac > Efrac_kaon ){
           Efrac_kaon = tmpEfrac;
           trackLength_kaon = track->Length();
           MCkaon_reco = particle;
         }
       }
       //michel from nucleon decay
       else if( (particle->PdgCode() == -11) && (particle->TrackId() == MC_michelID) ){
         if( tmpEfrac > Efrac_michel ){
           Efrac_michel = tmpEfrac;
           trackLength_michel = track->Length();
           MCmichel_reco = particle;
         }
       } 

    }

    Reco_LengthRes =  truth_lengthLepton-trackLength_lepton;
    Reco_LengthResProton = proton_length-trackLength_proton; 
    Reco_LengthResPionPlus = pion_plus_length-trackLength_pion_plus; 
    Reco_LengthResPionMinus = pion_minus_length-trackLength_pion_minus;

    if( MClepton_reco && MClepton  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){ 
        MC_LeptonTrack = 1;
        h_Pmu_num->Fill(MC_leptonP);
        h_Ev_num->Fill(MC_incoming_P[3]);
        h_theta_num->Fill(theta_mu);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_trackRes_lepton->Fill(Reco_LengthRes);  
      }
    }
    if( MCproton_reco && MCproton ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
        MC_ProtonTrack = 1;
        h_Pproton_num->Fill(MC_leading_ProtonP);     
        h_Efrac_proton->Fill(Efrac_proton);
        h_trackRes_proton->Fill(Reco_LengthResProton);       
      }
    }
    if( MCpion_plus_reco && MCpion_plus ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
        MC_PionPlusTrack = 1;
        h_Ppion_plus_num->Fill(MC_leading_PionPlusP);     
        h_Efrac_pion_plus->Fill(Efrac_pionplus);
        h_trackRes_pion_plus->Fill(Reco_LengthResPionPlus);
      }
    }
    if( MCpion_minus_reco && MCpion_minus  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ) {
        MC_PionMinusTrack = 1;
        h_Ppion_minus_num->Fill(MC_leading_PionMinusP);     
        h_Efrac_pion_minus->Fill(Efrac_pionminus);
        h_trackRes_pion_minus->Fill(Reco_LengthResPionMinus);
      }
    }
    if(!fisNeutrinoInt ){
      if( MClepton_reco && MClepton  ){
        MC_LeptonTrack = 1;
        h_Pmu_num->Fill(MC_leptonP);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_trackRes_lepton->Fill(Reco_LengthRes);
      }
      if( MCkaon_reco && MCkaon ){
        MC_kaonTrack = 1;
        h_Pkaon_num->Fill(MC_kaonP);
        h_Efrac_kaon->Fill(Efrac_kaon);
        h_trackRes_kaon->Fill(kaonLength-trackLength_kaon);
      }
      if( MCmichel_reco && MCmichel ){
        MC_michelTrack = 1;
        h_Pmichel_e_num->Fill(MC_michelP);
        h_Efrac_michel->Fill(Efrac_michel);
        h_trackRes_michel->Fill(michelLength-trackLength_michel);
      }

    }
    
}
//========================================================================
void NeutrinoTrackingEff::truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac){

    //std::cout<<"truthMatcher..."<<std::endl;
    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < track_hits.size(); ++j){
       art::Ptr<recob::Hit> hit = track_hits[j];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t k = 0; k < TrackIDs.size(); k++){
          trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
       }            
    }
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 
    //!since we are looking for muons/pions/protons this should be enough 
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         partial_E = ii->second;
         max_E = ii->second;
         TrackID = ii->first;
       }
    } 

    MCparticle = bt->TrackIDToParticle(TrackID);
    Efrac = partial_E/total_E;
    //std::cout<<"total "<<total_E<<" frac "<<Efrac<<" which particle "<<MCparticle->PdgCode()<<std::endl;
}
//========================================================================
double NeutrinoTrackingEff::truthLength( const simb::MCParticle *MCparticle ){
   //calculate the truth length considering only the part that is inside the TPC
   //Base on a peace of code from dune/TrackingAna/TrackingEfficiency_module.cc

   if( !MCparticle ) return -999.0;
   int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
   double TPCLengthHits[numberTrajectoryPoints];
   int FirstHit=0, LastHit=0;
   double TPCLength = 0.0;
   bool BeenInVolume = false;

   for(int MCHit=0; MCHit < numberTrajectoryPoints; ++MCHit) {
      const TLorentzVector& tmpPosition= MCparticle->Position(MCHit);
      double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      if (MCHit!=0) TPCLengthHits[MCHit] = sqrt( pow( (MCparticle->Vx(MCHit-1)-MCparticle->Vx(MCHit)),2)+ pow( (MCparticle->Vy(MCHit-1)-MCparticle->Vy(MCHit)),2)+ pow( (MCparticle->Vz(MCHit-1)-MCparticle->Vz(MCHit)),2));
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if(tpcid.isValid) {
        // -- Check if hit is within drift window...
        geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
        geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
        double XPlanePosition      = tpc.PlaneLocation(0)[0];
        double DriftTimeCorrection = fabs( tmpPosition[0] - XPlanePosition ) / XDriftVelocity;
        double TimeAtPlane         = MCparticle->T() + DriftTimeCorrection; 
        if( TimeAtPlane < detprop->TriggerOffset() || TimeAtPlane > detprop->TriggerOffset() + WindowSize ) continue;
        LastHit = MCHit;
        if( !BeenInVolume ) {
	  BeenInVolume = true;
          FirstHit = MCHit;
	}
      }		
   }
   for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) TPCLength += TPCLengthHits[Hit];
   return TPCLength;
}
//========================================================================
bool NeutrinoTrackingEff::insideFV( double vertex[4]){ 

     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     if (x>fFidVolXmin && x<fFidVolXmax&&
	 y>fFidVolYmin && y<fFidVolYmax&&
	 z>fFidVolZmin && z<fFidVolZmax)
       return true;
     else
       return false;
}
//========================================================================
void NeutrinoTrackingEff::doEfficiencies(){

   art::ServiceHandle<art::TFileService> tfs;

   if(TEfficiency::CheckConsistency(*h_Ev_num,*h_Ev_den)){
      h_Eff_Ev = tfs->make<TEfficiency>(*h_Ev_num,*h_Ev_den);
      TGraphAsymmErrors *grEff_Ev = h_Eff_Ev->CreateGraph();
      grEff_Ev->Write("grEff_Ev");
      h_Eff_Ev->Write("h_Eff_Ev");
   }
   if(TEfficiency::CheckConsistency(*h_Pmu_num,*h_Pmu_den)){ 
     h_Eff_Pmu = tfs->make<TEfficiency>(*h_Pmu_num,*h_Pmu_den);
     TGraphAsymmErrors *grEff_Pmu = h_Eff_Pmu->CreateGraph();
     grEff_Pmu->Write("grEff_Pmu");
     h_Eff_Pmu->Write("h_Eff_Pmu");
   }
   if(TEfficiency::CheckConsistency(*h_theta_num,*h_theta_den)){
     h_Eff_theta = tfs->make<TEfficiency>(*h_theta_num,*h_theta_den);
     TGraphAsymmErrors *grEff_theta = h_Eff_theta->CreateGraph();
     grEff_theta->Write("grEff_theta");
     h_Eff_theta->Write("h_Eff_theta");
   }
   if(TEfficiency::CheckConsistency(*h_Pproton_num,*h_Pproton_den)){
     h_Eff_Pproton = tfs->make<TEfficiency>(*h_Pproton_num,*h_Pproton_den);
     TGraphAsymmErrors *grEff_Pproton = h_Eff_Pproton->CreateGraph();
     grEff_Pproton->Write("grEff_Pproton");
     h_Eff_Pproton->Write("h_Eff_Pproton");
   }
   if(TEfficiency::CheckConsistency(*h_Ppion_plus_num,*h_Ppion_plus_den)){
     h_Eff_Ppion_plus = tfs->make<TEfficiency>(*h_Ppion_plus_num,*h_Ppion_plus_den);
     TGraphAsymmErrors *grEff_Ppion_plus = h_Eff_Ppion_plus->CreateGraph();
     grEff_Ppion_plus->Write("grEff_Ppion_plus");
     h_Eff_Ppion_plus->Write("h_Eff_Ppion_plus");
   }
   if(TEfficiency::CheckConsistency(*h_Ppion_minus_num,*h_Ppion_minus_den)){
     h_Eff_Ppion_minus = tfs->make<TEfficiency>(*h_Ppion_minus_num,*h_Ppion_minus_den);
     TGraphAsymmErrors *grEff_Ppion_minus = h_Eff_Ppion_minus->CreateGraph();
     grEff_Ppion_minus->Write("grEff_Ppion_minus");
     h_Eff_Ppion_minus->Write("h_Eff_Ppion_minus");
   }

   if(!fisNeutrinoInt ){ 
     if(TEfficiency::CheckConsistency(*h_Pkaon_num,*h_Pkaon_den)){
       h_Eff_Pkaon = tfs->make<TEfficiency>(*h_Pkaon_num,*h_Pkaon_den);
       TGraphAsymmErrors *grEff_Pkaon = h_Eff_Pkaon->CreateGraph();
       grEff_Pkaon->Write("grEff_Pkaon");
       h_Eff_Pkaon->Write("h_Eff_Pkaonn");
     }
     if(TEfficiency::CheckConsistency(*h_Pmichel_e_num,*h_Pmichel_e_den)){
       h_Eff_Pmichel = tfs->make<TEfficiency>(*h_Pmichel_e_num,*h_Pmichel_e_den);
       TGraphAsymmErrors *grEff_Pmichel = h_Eff_Pmichel->CreateGraph();
       grEff_Pmichel->Write("grEff_Pmichel");
       h_Eff_Pmichel->Write("h_Eff_Pmichel");
     }

   }

}
//========================================================================
void NeutrinoTrackingEff::reset(){

   MC_incoming_PDG = -999;
   MC_lepton_PDG =-999;
   MC_isCC =-999;
   MC_channel =-999;
   MC_target =-999;
   MC_Q2 =-999.0;
   MC_W =-999.0;
   MC_lepton_theta = -999.0;
   MC_leptonID = -999;
   MC_leading_protonID = -999;
   MC_leading_ProtonP = -999.0;
   MC_leading_PionPlusID = -999;
   MC_leading_PionMinusID = -999;
   MC_leading_PionPlusP = -999.0;
   MC_leading_PionMinusP = -999.0;
   MC_LeptonTrack = -999;
   MC_ProtonTrack =-999;
   MC_PionPlusTrack =-999; 
   MC_PionMinusTrack =-999;
   MC_kaonID = -999;
   MC_kaonP = -999.0;
   MC_LeptonTrack = -999;
   MC_kaonTrack =-999;
   MC_michelID = -999;
   MC_michelP = -999.0;
   MC_michelTrack =-999;
 
   for(int i = 0; i<4; i++){
      MC_vertex[i] = 0;
      MC_incoming_P[i] =  0;
      MC_lepton_startMomentum[i] = 0;
      MC_lepton_endMomentum[i] = 0;
      MC_lepton_startXYZT[i] = 0;
      MC_lepton_endXYZT[i] = 0;
      MC_kaon_startMomentum[i] = 0;
      MC_kaon_endMomentum[i] = 0;
      MC_kaon_startXYZT[i] = 0;
      MC_kaon_endXYZT[i] = 0;
      MC_michel_startMomentum[i] = 0;
      MC_michel_endMomentum[i] = 0;
      MC_michel_startXYZT[i] = 0;
      MC_michel_endXYZT[i] = 0;
  
   }

   Reco_EfracLepton = -999.0; 
   Reco_LengthRes = -999.0;
   Reco_LengthResProton = -999.0;
   Reco_LengthResPionPlus = -999.0;
   Reco_LengthResPionMinus = -999.0;
   Reco_LengthReskaon = -999.0;
   Reco_LengthResmichel = -999.0;
  
   for(int i=0; i<MAX_TRACKS; i++) {
       MC_id[i] = 0;
       MC_pdg[i] = 0;
       MC_mother[i] = 0;
       for(int j=0; j<4; j++) {
          MC_startXYZT[i][j]      = 0;
          MC_endXYZT[i][j]        = 0;
          MC_startMomentum[i][j] = 0;
          MC_endMomentum[i][j]   = 0;
       }
    }

}
//========================================================================
DEFINE_ART_MODULE(NeutrinoTrackingEff)

} 

#endif // NeutrinoTrackingEff_Module
