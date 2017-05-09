#ifndef NeutrinoShowerEff_Module
#define NeutrinoShowerEff_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/ArtDataHelper/MVAReader.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#define MAX_SHOWERS 1000
using namespace std;


//========================================================================

namespace DUNE{

  class NeutrinoShowerEff : public art::EDAnalyzer {
  public:

    explicit NeutrinoShowerEff(fhicl::ParameterSet const& pset);
    virtual ~NeutrinoShowerEff();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    
    void processEff(const art::Event& evt, bool &isFiducial);
    void truthMatcher(std::vector<art::Ptr<recob::Hit>>all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
    template <size_t N> void checkCNNtrkshw(const art::Event& evt, std::vector<art::Ptr<recob::Hit>>all_hits);
    bool insideFV(double vertex[4]);
    void doEfficiencies();
    void reset();

  private:

    // the parameters we'll read from the .fcl
    
    art::InputTag fMCTruthModuleLabel;
    art::InputTag fHitModuleLabel;
    art::InputTag fShowerModuleLabel;
    art::InputTag fCNNEMModuleLabel;
    int           fNeutrinoPDGcode;
    int	       	  fLeptonPDGcode;
    double        fMaxNeutrinoE;
    bool	  fSaveMCTree; 
    double 	  fMaxEfrac;
    double        fMinCompleteness;

 
    TTree *fEventTree;
    TTree *fHitsTree;

    TH1D *h_Ev_den;
    TH1D *h_Ev_num;

    TH1D *h_Ee_den;
    TH1D *h_Ee_num;

    TH1D *h_Pe_den;
    TH1D *h_Pe_num;

    TH1D *h_theta_den;
    TH1D *h_theta_num;

    TH1D *h_Efrac_shContamination;
    TH1D *h_Efrac_shPurity;     
    TH1D *h_Ecomplet_lepton;     

    TH1D *h_Efrac_NueCCPurity;     //Signal
    TH1D *h_Ecomplet_NueCC;     
    TH1D *h_HighestHitsProducedParticlePDG_NueCC;     
    TH1D *h_HighestHitsProducedParticlePDG_bkg;     
    TH1D *h_Efrac_bkgPurity;     //Background
    TH1D *h_Ecomplet_bkg;     


    TEfficiency* h_Eff_Ev = 0;
    TEfficiency* h_Eff_Ee = 0;  
    TEfficiency* h_Eff_Pe = 0;
    TEfficiency* h_Eff_theta = 0;


    //Electron shower Best plane number
    TH1D *h_esh_bestplane_NueCC;
    TH1D *h_esh_bestplane_NC;
    TH1D *h_dEdX_electronorpositron_NueCC;
    TH1D *h_dEdX_electronorpositron_NC;
    TH1D *h_dEdX_photon_NueCC;
    TH1D *h_dEdX_photon_NC;
    TH1D *h_dEdX_proton_NueCC;
    TH1D *h_dEdX_proton_NC;
    TH1D *h_dEdX_neutron_NueCC;
    TH1D *h_dEdX_neutron_NC;
    TH1D *h_dEdX_chargedpion_NueCC;
    TH1D *h_dEdX_chargedpion_NC;
    TH1D *h_dEdX_neutralpion_NueCC;
    TH1D *h_dEdX_neutralpion_NC;
    TH1D *h_dEdX_everythingelse_NueCC;
    TH1D *h_dEdX_everythingelse_NC;

    //Study CNN track/shower id
    TH1D *h_trklike_em;
    TH1D *h_trklike_nonem;

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
    int    MC_leptonID;
    int    MC_LeptonTrack;
 
    double sh_direction_X[MAX_SHOWERS];
    double sh_direction_Y[MAX_SHOWERS];
    double sh_direction_Z[MAX_SHOWERS];
    double sh_start_X[MAX_SHOWERS];
    double sh_start_Y[MAX_SHOWERS];
    double sh_start_Z[MAX_SHOWERS];
    double sh_energy[MAX_SHOWERS][3];
    double sh_MIPenergy[MAX_SHOWERS][3];
    double sh_dEdx[MAX_SHOWERS][3];
    int    sh_bestplane[MAX_SHOWERS];
    double sh_length[MAX_SHOWERS];
    int    sh_hasPrimary_e[MAX_SHOWERS];
    double sh_Efrac_contamination[MAX_SHOWERS];
    int    sh_nHits[MAX_SHOWERS];
    int    n_recoShowers;
    double sh_Efrac_best;



    float fFidVolCutX;
    float fFidVolCutY;
    float fFidVolCutZ;

    float fFidVolXmin;
    float fFidVolXmax;
    float fFidVolYmin;
    float fFidVolYmax;
    float fFidVolZmin;
    float fFidVolZmax;

    art::ServiceHandle<geo::Geometry> geom;





  }; // class NeutrinoShowerEff


  //========================================================================
  NeutrinoShowerEff::NeutrinoShowerEff(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
    
  }
  //========================================================================
  NeutrinoShowerEff::~NeutrinoShowerEff(){
    //destructor
  }
  //========================================================================
  void NeutrinoShowerEff::reconfigure(fhicl::ParameterSet const& p){

 
    fMCTruthModuleLabel  = p.get<art::InputTag>("MCTruthModuleLabel");
    fHitModuleLabel      = p.get<art::InputTag>("HitModuleLabel");
    fShowerModuleLabel   = p.get<art::InputTag>("ShowerModuleLabel");
    fCNNEMModuleLabel    = p.get<art::InputTag>("CNNEMModuleLabel","");
    fLeptonPDGcode       = p.get<int>("LeptonPDGcode");
    fNeutrinoPDGcode     = p.get<int>("NeutrinoPDGcode");
    fMaxNeutrinoE	 = p.get<double>("MaxNeutrinoE");
    fMaxEfrac		 = p.get<double>("MaxEfrac");
    fMinCompleteness     = p.get<double>("MinCompleteness");
    fSaveMCTree		 = p.get<bool>("SaveMCTree");
    fFidVolCutX          = p.get<float>("FidVolCutX");
    fFidVolCutY          = p.get<float>("FidVolCutY");
    fFidVolCutZ          = p.get<float>("FidVolCutZ");}
  //========================================================================
  //========================================================================
  void NeutrinoShowerEff::beginJob(){
    cout<<"job begin..."<<endl;

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
    //  double Pbins[18] ={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0};
  
  
    h_Ev_den = tfs->make<TH1D>("h_Ev_den","Neutrino Energy; Neutrino Energy (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Ev_den->Sumw2();
    h_Ev_num = tfs->make<TH1D>("h_Ev_num","Neutrino Energy; Neutrino Energy (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Ev_num->Sumw2();

    h_Ee_den = tfs->make<TH1D>("h_Ee_den","Electron Energy; Electron Energy (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Ee_den->Sumw2();
    h_Ee_num = tfs->make<TH1D>("h_Ee_num","Electron Energy; Electron Energy (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Ee_num->Sumw2();

    h_Pe_den = tfs->make<TH1D>("h_Pe_den","Electron Momentum; Electron Momentum (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Pe_den->Sumw2();
    h_Pe_num = tfs->make<TH1D>("h_Pe_num","Electron Momentum; Electron Momentum (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Pe_num->Sumw2();

    h_theta_den = tfs->make<TH1D>("h_theta_den","Theta; Theta w.r.t beam direction (Degrees); Shower reconstruction Efficiency",43,theta_bin);
    h_theta_den->Sumw2();
    h_theta_num = tfs->make<TH1D>("h_theta_num","Theta; Theta w.r.t beam direction (Degrees); Shower reconstruction Efficiency",43,theta_bin);
    h_theta_num->Sumw2();

    h_Efrac_shContamination = tfs->make<TH1D>("h_Efrac_shContamination","Efrac Lepton; Energy fraction (contamination);",60,0,1.2);
    h_Efrac_shContamination->Sumw2();
    h_Efrac_shPurity = tfs->make<TH1D>("h_Efrac_shPurity","Efrac Lepton; Energy fraction (Purity);",60,0,1.2);
    h_Efrac_shPurity->Sumw2();
    h_Ecomplet_lepton = tfs->make<TH1D>("h_Ecomplet_lepton","Ecomplet Lepton; Track Completeness;",60,0,1.2);
    h_Ecomplet_lepton->Sumw2();

    h_HighestHitsProducedParticlePDG_NueCC= tfs->make<TH1D>("h_HighestHitsProducedParticlePDG_NueCC","PDG Code; PDG Code;",4,-0.5,3.5);//0 for undefined, 1=electron, 2=photon, 3=anything else     //Signal
    h_HighestHitsProducedParticlePDG_NueCC->Sumw2();     
    h_HighestHitsProducedParticlePDG_bkg= tfs->make<TH1D>("h_HighestHitsProducedParticlePDG_bkg","PDG Code; PDG Code;",4,-0.5,3.5);//0 for undefined, 1=electron, 2=photon, 3=anything else     //bkg    
    h_HighestHitsProducedParticlePDG_bkg->Sumw2();


    h_Efrac_NueCCPurity= tfs->make<TH1D>("h_Efrac_NueCCPurity","Efrac NueCC; Energy fraction (Purity);",60,0,1.2);     //Signal
    h_Efrac_NueCCPurity->Sumw2();
    h_Ecomplet_NueCC= tfs->make<TH1D>("h_Ecomplet_NueCC","Ecomplet NueCC; Track Completeness;",60,0,1.2);     
    h_Ecomplet_NueCC->Sumw2();

    
    h_Efrac_bkgPurity= tfs->make<TH1D>("h_Efrac_bkgPurity","Efrac bkg; Energy fraction (Purity);",60,0,1.2);     //Background
    h_Efrac_bkgPurity->Sumw2();
    h_Ecomplet_bkg= tfs->make<TH1D>("h_Ecomplet_bkg","Ecomplet bkg; Track Completeness;",60,0,1.2);     
    h_Ecomplet_bkg->Sumw2();



    h_esh_bestplane_NueCC=tfs->make<TH1D>("h_esh_bestplane_NueCC","Best plane; Best plane;",4,-0.5,3.5); 
    h_esh_bestplane_NC=tfs->make<TH1D>("h_esh_bestplane_NC","Best plane; Best plane;",4,-0.5,3.5); 
    h_dEdX_electronorpositron_NueCC=tfs->make<TH1D>("h_dEdX_electronorpositron_NueCC","dE/dX; Electron or Positron dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_electronorpositron_NC=tfs->make<TH1D>("h_dEdX_electronorpositron_NC","dE/dX; Electron or Positron dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_photon_NueCC=tfs->make<TH1D>("h_dEdX_photon_NueCC","dE/dX; photon dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_photon_NC=tfs->make<TH1D>("h_dEdX_photon_NC","dE/dX; photon dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_proton_NueCC=tfs->make<TH1D>("h_dEdX_proton_NueCC","dE/dX; proton dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_proton_NC=tfs->make<TH1D>("h_dEdX_proton_NC","dE/dX; proton dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_neutron_NueCC=tfs->make<TH1D>("h_dEdX_neutron_NueCC","dE/dX; neutron dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_neutron_NC=tfs->make<TH1D>("h_dEdX_neutron_NC","dE/dX; neutron dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_chargedpion_NueCC=tfs->make<TH1D>("h_dEdX_chargedpion_NueCC","dE/dX; charged pion dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_chargedpion_NC=tfs->make<TH1D>("h_dEdX_chargedpion_NC","dE/dX; charged pion dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_neutralpion_NueCC=tfs->make<TH1D>("h_dEdX_neutralpion_NueCC","dE/dX; neutral pion dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_neutralpion_NC=tfs->make<TH1D>("h_dEdX_neutralpion_NC","dE/dX; neutral pion dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_everythingelse_NueCC=tfs->make<TH1D>("h_dEdX_everythingelse_NueCC","dE/dX; everythingelse dE/dX (MeV/cm);",100,0.0,15.0); 
    h_dEdX_everythingelse_NC=tfs->make<TH1D>("h_dEdX_everythingelse_NC","dE/dX; everythingelse dE/dX (MeV/cm);",100,0.0,15.0); 
    

    h_esh_bestplane_NueCC->Sumw2();
    h_esh_bestplane_NC->Sumw2();
    h_dEdX_electronorpositron_NueCC->Sumw2();
    h_dEdX_electronorpositron_NC->Sumw2();
    h_dEdX_photon_NueCC->Sumw2();
    h_dEdX_photon_NC->Sumw2();
    h_dEdX_proton_NueCC->Sumw2();
    h_dEdX_proton_NC->Sumw2();
    h_dEdX_neutron_NueCC->Sumw2();
    h_dEdX_neutron_NC->Sumw2();
    h_dEdX_chargedpion_NueCC->Sumw2();
    h_dEdX_chargedpion_NC->Sumw2();
    h_dEdX_neutralpion_NueCC->Sumw2();
    h_dEdX_neutralpion_NC->Sumw2();
    h_dEdX_everythingelse_NueCC->Sumw2();
    h_dEdX_everythingelse_NC->Sumw2();
    

    h_trklike_em = tfs->make<TH1D>("h_trklike_em","EM hits; Track-like Score;",100,0,1);
    h_trklike_nonem = tfs->make<TH1D>("h_trklike_nonem","Non-EM hits; Track-like Score;",100,0,1);
    
    if( fSaveMCTree ){
      fEventTree = new TTree("Event", "Event Tree from Sim & Reco");
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
      fEventTree->Branch("mc_leptonID", &MC_leptonID, "mc_leptonID/I");

      fEventTree->Branch("n_showers", &n_recoShowers);
      fEventTree->Branch("sh_direction_X", &sh_direction_X, "sh_direction_X[n_showers]/D");
      fEventTree->Branch("sh_direction_Y", &sh_direction_Y, "sh_direction_Y[n_showers]/D");
      fEventTree->Branch("sh_direction_Z", &sh_direction_Z, "sh_direction_Z[n_showers]/D");
      fEventTree->Branch("sh_start_X", &sh_start_X, "sh_start_X[n_showers]/D");
      fEventTree->Branch("sh_start_Y", &sh_start_Y, "sh_start_Y[n_showers]/D");
      fEventTree->Branch("sh_start_Z", &sh_start_Z, "sh_start_Z[n_showers]/D");
      fEventTree->Branch("sh_energy", &sh_energy, "sh_energy[n_showers][3]/D");
      fEventTree->Branch("sh_MIPenergy", &sh_MIPenergy, "sh_MIPenergy[n_showers][3]/D");
      fEventTree->Branch("sh_dEdx", &sh_dEdx, "sh_dEdx[n_showers][3]/D");
      fEventTree->Branch("sh_bestplane", &sh_bestplane, "sh_bestplane[n_showers]/I");
      fEventTree->Branch("sh_length", &sh_length, "sh_length[n_showers]/D");
      fEventTree->Branch("sh_hasPrimary_e", &sh_hasPrimary_e, "sh_hasPrimary_e[n_showers]/I");
      fEventTree->Branch("sh_Efrac_contamination", &sh_Efrac_contamination, "sh_Efrac_contamination[n_showers]/D");
      fEventTree->Branch("sh_Efrac_best", &sh_Efrac_best, "sh_Efrac_best/D");
      fEventTree->Branch("sh_nHits",&sh_nHits, "sh_nHits[n_showers]/I");
      
    }

  }
  //========================================================================
  void NeutrinoShowerEff::endJob(){
    doEfficiencies();
  }
  //========================================================================
  void NeutrinoShowerEff::beginRun(const art::Run& /*run*/){
    mf::LogInfo("NeutrinoShowerEff")<<"begin run..."<<endl;
  }
  //========================================================================
  void NeutrinoShowerEff::analyze( const art::Event& event ){

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
  void NeutrinoShowerEff::processEff( const art::Event& event, bool &isFiducial){

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
        MC_incoming_PDG = std::abs(nu.Nu().PdgCode());
        MC_Q2 = nu.QSqr();
        MC_channel = nu.InteractionType();
        MC_W = nu.W();
        const TLorentzVector& nu_momentum = nu.Nu().Momentum(0); 
        nu_momentum.GetXYZT(MC_incoming_P); 
        const TLorentzVector& vertex =neutrino.Position(0); 
        vertex.GetXYZT(MC_vertex);
        simb::MCParticle lepton = nu.Lepton();
        MC_lepton_PDG = lepton.PdgCode();
        //cout<<"Incoming E "<<MC_incoming_P[3]<<" is CC? "<<MC_isCC<<" nuPDG "<<MC_incoming_PDG<<" target "<<MC_target<<" vtx "<<MC_vertex[0]<<" "<<MC_vertex[1]<<" "<<MC_vertex[2]<<" "<<MC_vertex[3]<<endl;
      }
    }

    //!save lepton 
    simb::MCParticle *MClepton = NULL; 
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;

    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
      particle = ipar->second;
      if( std::abs(particle->PdgCode()) == fLeptonPDGcode && particle->Mother() == 0 ){  //primary lepton
        const TLorentzVector& lepton_momentum =particle->Momentum(0); 
        const TLorentzVector& lepton_position =particle->Position(0); 
        const TLorentzVector& lepton_positionEnd   = particle->EndPosition();
        const TLorentzVector& lepton_momentumEnd   = particle->EndMomentum();
        lepton_momentum.GetXYZT(MC_lepton_startMomentum);
        lepton_position.GetXYZT(MC_lepton_startXYZT);
        lepton_positionEnd.GetXYZT(MC_lepton_endXYZT);
        lepton_momentumEnd.GetXYZT(MC_lepton_endMomentum);
        MC_leptonID = particle->TrackId();
        MClepton = particle;
      }
    } 
    //Saving denominator histograms for lepton pions and protons 
    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
    double Pe = sqrt(pow(MC_lepton_startMomentum[0],2)+pow(MC_lepton_startMomentum[1],2)+pow(MC_lepton_startMomentum[2],2));
    double Pv  = sqrt(pow(MC_incoming_P[0],2)+pow(MC_incoming_P[1],2)+pow(MC_incoming_P[2],2));
    double theta_e = acos((MC_incoming_P[0]*MC_lepton_startMomentum[0] + MC_incoming_P[1]*MC_lepton_startMomentum[1] +MC_incoming_P[2]*MC_lepton_startMomentum[2])/(Pv*Pe) );
    theta_e *= (180.0/3.14159);

    //save CC events within the fiducial volume with the favorite neutrino flavor 
    if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) ){
      if( MClepton ){
        h_Ev_den->Fill(MC_incoming_P[3]);
        h_Ee_den->Fill(MC_lepton_startMomentum[3]);
        h_Pe_den->Fill(Pe);
        h_theta_den->Fill(theta_e);
      }
    } 
   
    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================
    art::Handle<std::vector<recob::Shower>> showerHandle;
    if(!event.getByLabel(fShowerModuleLabel,showerHandle)){
      mf::LogError("NeutrinoShowerEff")<<"Could not find shower with label "<<fShowerModuleLabel.encode();
      return;
    }
    std::vector<art::Ptr<recob::Shower>> showerlist;
    art::fill_ptr_vector(showerlist, showerHandle);

    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> all_hits;
    if(event.getByLabel(fHitModuleLabel,hitHandle)){
      art::fill_ptr_vector(all_hits, hitHandle);
    }

    n_recoShowers= showerlist.size();
    //if ( n_recoShowers == 0 || n_recoShowers> MAX_SHOWERS ) return;
    art::FindManyP<recob::Hit> sh_hitsAll(showerHandle, event, fShowerModuleLabel);
    cout<<"Found this many showers "<<n_recoShowers<<endl; 
    double Efrac_contamination= 999.0;
    double Efrac_contaminationNueCC= 999.0;
    

    double Ecomplet_lepton =0.0;
    double Ecomplet_NueCC =0.0;
    int ParticlePGD_HighestShHits=0;//undefined
    int shower_bestplane=0;
    double Showerparticlededx_inbestplane=0.0;
    int showerPDGwithHighestHitsforFillingdEdX=0;//0=undefined,1=electronorpositronshower,2=photonshower,3=protonshower,4=neutronshower,5=chargedpionshower,6=neutralpionshower,7=everythingelseshower

    const simb::MCParticle *MClepton_reco = NULL; 
    int nHits =0;
    for(int i=0; i<n_recoShowers; i++){
      
      art::Ptr<recob::Shower> shower = showerlist[i];
      sh_direction_X[i] = shower->Direction().X();  
      sh_direction_Y[i] = shower->Direction().Y();  
      sh_direction_Z[i] = shower->Direction().Z();  
      sh_start_X[i] = shower->ShowerStart().X();
      sh_start_Y[i] = shower->ShowerStart().Y();
      sh_start_Z[i] = shower->ShowerStart().Z();
      sh_bestplane[i] = shower->best_plane();
      sh_length[i] = shower->Length();
      for( size_t j =0; j<shower->Energy().size(); j ++) sh_energy[i][j] = shower->Energy()[j];
      for( size_t j =0; j<shower->MIPEnergy().size(); j++) sh_MIPenergy[i][j] = shower->MIPEnergy()[j];
      for( size_t j =0; j<shower->dEdx().size(); j++) sh_dEdx[i][j] = shower->dEdx()[j];

      std::vector<art::Ptr<recob::Hit>> sh_hits = sh_hitsAll.at(i);  

      //  std::cout<<" shower best plane:"<<shower->best_plane()<<" shower dEdx size:"<<shower->dEdx().size()<<std::endl;
      //for( size_t j =0; j<shower->dEdx().size(); j++) std::cout<<shower->dEdx()[j]<<" ";

      const simb::MCParticle *particle;
      double tmpEfrac_contamination = 0.0;  //fraction of non EM energy contatiminatio (see truthMatcher for definition)
  
              
      double tmpEcomplet =0;
  
      int tmp_nHits = sh_hits.size();
      truthMatcher( all_hits, sh_hits, particle, tmpEfrac_contamination,tmpEcomplet);
      //truthMatcher( all_hits, sh_hits, particle, tmpEfrac_contaminationNueCC,tmpEcompletNueCC );
       


      sh_Efrac_contamination[i] = tmpEfrac_contamination;
      sh_nHits[i] = tmp_nHits; 
      sh_hasPrimary_e[i] = 0;
       

      //Shower with highest hits       
      if( tmp_nHits > nHits ){
        nHits = tmp_nHits;
        Ecomplet_NueCC =tmpEcomplet;
        Efrac_contaminationNueCC = tmpEfrac_contamination; 
	 
        if(std::abs(particle->PdgCode())==11){
          ParticlePGD_HighestShHits=1;
        }else if(particle->PdgCode()==22){
          ParticlePGD_HighestShHits=2;
        }else{
          ParticlePGD_HighestShHits=3;
        }
	 
	 

        //dedx for different showers
        //Highest hits shower pdg for the dEdx study 0=undefined,1=electronorpositronshower,2=photonshower,3=protonshower,4=neutronshower,5=chargedpionshower,6=neutralpionshower,7=everythingelseshower
        shower_bestplane=shower->best_plane();
        if (shower_bestplane<0 || shower_bestplane>=int(shower->dEdx().size())){
          //bestplane is not set properly, just pick the first plane that has dEdx
          for (size_t i = 0; i<shower->dEdx().size(); ++i){
            if (shower->dEdx()[i]){
              shower_bestplane = i;
              break;
            }
          }
        }
        if (shower_bestplane<0 || shower_bestplane>=int(shower->dEdx().size())){
          //still a problem? just set it to 0
          shower_bestplane = 0;
        }
          
        if (shower_bestplane>=0 and shower_bestplane<int(shower->dEdx().size()))
          Showerparticlededx_inbestplane=shower->dEdx()[shower_bestplane]; 
	   
        if(std::abs(particle->PdgCode())==11){//lepton shower
          showerPDGwithHighestHitsforFillingdEdX=1;
        }else if(particle->PdgCode()==22){//photon shower
          showerPDGwithHighestHitsforFillingdEdX=2;
        }else if(particle->PdgCode()==2212){//proton shower
          showerPDGwithHighestHitsforFillingdEdX=3;
        }else if(particle->PdgCode()==2112){//neutron shower
          showerPDGwithHighestHitsforFillingdEdX=4;
        }else if(std::abs(particle->PdgCode())==211){//charged pion shower
          showerPDGwithHighestHitsforFillingdEdX=5;
        }else if(particle->PdgCode()==111){//neutral pion shower
          showerPDGwithHighestHitsforFillingdEdX=6;
        }else{//everythingelse shower
          showerPDGwithHighestHitsforFillingdEdX=7;
        }
	 

        //Efrac_contamination = tmpEfrac_contamination;
        //MClepton_reco = particle;
        //sh_Efrac_best =Efrac_contamination; 
        //cout<<"this is the best shower "<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac_contamination<<" "<<sh_hits.size()<<endl;
      } 



      if( particle->PdgCode()  == fLeptonPDGcode && particle->TrackId() == MC_leptonID ) sh_hasPrimary_e[i] = 1;
      //cout<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac_contamination<<" "<<sh_hits.size()<<" "<<particle->TrackId()<<" "<<MC_leptonID<<endl;
      //save the best shower based on non EM and number of hits

      if( particle->PdgCode()  == fLeptonPDGcode && particle->TrackId() == MC_leptonID ){

        if(tmpEcomplet>Ecomplet_lepton){

          Ecomplet_lepton = tmpEcomplet;
	   
          Efrac_contamination = tmpEfrac_contamination;
          MClepton_reco = particle;
          sh_Efrac_best =Efrac_contamination; 
           
        }
      }          
    }//end of looping all the showers
   
    if( MClepton_reco && MClepton  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) ){ 
        h_Efrac_shContamination->Fill(Efrac_contamination);
        h_Efrac_shPurity->Fill(1-Efrac_contamination);
	h_Ecomplet_lepton->Fill(Ecomplet_lepton);
	
	// Selecting good showers requires completeness of gretaer than 70 % and Purity > 70 %	
	if( Efrac_contamination < fMaxEfrac && Ecomplet_lepton> fMinCompleteness ){
        
	  h_Pe_num->Fill(Pe);
          h_Ev_num->Fill(MC_incoming_P[3]);
	  h_Ee_num->Fill(MC_lepton_startMomentum[3]);
	  h_theta_num->Fill(theta_e);
	}
      }
    }

    //NueCC SIgnal and background Completeness
    if(MC_isCC==1
       &&(fNeutrinoPDGcode == std::abs(MC_incoming_PDG))
       &&isFiducial){
        h_HighestHitsProducedParticlePDG_NueCC->Fill(ParticlePGD_HighestShHits);
	
        if(ParticlePGD_HighestShHits>0){// atleat one shower is reconstructed
          h_Ecomplet_NueCC->Fill(Ecomplet_NueCC);
          h_Efrac_NueCCPurity->Fill(1-Efrac_contaminationNueCC);    
	  
          h_esh_bestplane_NueCC->Fill(shower_bestplane);
          if(showerPDGwithHighestHitsforFillingdEdX==1)//electron or positron shower
            {
              h_dEdX_electronorpositron_NueCC->Fill(Showerparticlededx_inbestplane);
            }else if(showerPDGwithHighestHitsforFillingdEdX==2)//photon shower
            {
              h_dEdX_photon_NueCC->Fill(Showerparticlededx_inbestplane);
            }else if(showerPDGwithHighestHitsforFillingdEdX==3)//proton shower
            {
              h_dEdX_proton_NueCC->Fill(Showerparticlededx_inbestplane);
            }else if(showerPDGwithHighestHitsforFillingdEdX==4)//neutron shower
            {
              h_dEdX_neutron_NueCC->Fill(Showerparticlededx_inbestplane);
            }else if(showerPDGwithHighestHitsforFillingdEdX==5)//charged pion shower
            {
              h_dEdX_chargedpion_NueCC->Fill(Showerparticlededx_inbestplane);
            }else if(showerPDGwithHighestHitsforFillingdEdX==6)//neutral pion shower
            {
              h_dEdX_neutralpion_NueCC->Fill(Showerparticlededx_inbestplane);
            }else if(showerPDGwithHighestHitsforFillingdEdX==7)//everythingelse shower
            {
              h_dEdX_everythingelse_NueCC->Fill(Showerparticlededx_inbestplane);
            }
        }
    }
    else if(!MC_isCC&&
	    isFiducial){
      h_HighestHitsProducedParticlePDG_bkg->Fill(ParticlePGD_HighestShHits);
      
      
      if(ParticlePGD_HighestShHits>0){
	h_Ecomplet_bkg->Fill(Ecomplet_NueCC);
	h_Efrac_bkgPurity->Fill(1-Efrac_contaminationNueCC);	
        
        
	h_esh_bestplane_NC->Fill(shower_bestplane);
        if(showerPDGwithHighestHitsforFillingdEdX==1)//electron or positron shower
          {
            h_dEdX_electronorpositron_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==2)//photon shower
          {
            h_dEdX_photon_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==3)//proton shower
          {
            h_dEdX_proton_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==4)//neutron shower
          {
            h_dEdX_neutron_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==5)//charged pion shower
          {
            h_dEdX_chargedpion_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==6)//neutral pion shower
          {
            h_dEdX_neutralpion_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==7)//everythingelse shower
          {
            h_dEdX_everythingelse_NC->Fill(Showerparticlededx_inbestplane);
          }
      }//if(ParticlePGD_HighestShHits>0)
    }//else if(!MC_isCC&&isFiducial)

    checkCNNtrkshw<4>(event, all_hits);
  }

  //========================================================================
  void NeutrinoShowerEff::truthMatcher(std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

    MCparticle=0;
    Efrac=1.0;
    Ecomplet=0;

    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < shower_hits.size(); ++j){
      art::Ptr<recob::Hit> hit = shower_hits[j];
      //For know let's use collection plane to look at the shower reconstruction
      //if( hit->View() != 2) continue;
      std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
      for(size_t k = 0; k < TrackIDs.size(); k++){
        if (trkID_E.find(std::abs(TrackIDs[k].trackID))==trkID_E.end()) trkID_E[std::abs(TrackIDs[k].trackID)] = 0;
        trkID_E[std::abs(TrackIDs[k].trackID)] += TrackIDs[k].energy;
      }            
    }
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E=0.0;
    //double noEM_E = 0.0;  //non electromagnetic energy is defined as energy from charged pion and protons 
    if( !trkID_E.size() ) return; //Ghost shower???
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
      total_E += ii->second;
      if((ii->second)>max_E){
        partial_E = ii->second;
        max_E = ii->second;
        TrackID = ii->first;
      }
      //int ID = ii->first;
      // const simb::MCParticle *particle = bt->TrackIDToParticle(ID);
      //if( abs(particle->PdgCode()) == 211 || particle->PdgCode() == 2212 ){
      //if( particle->PdgCode() != 22 && abs(particle->PdgCode()) != 11){
      //noEM_E += ii->second;
      //}

    } 
    MCparticle = bt->TrackIDToParticle(TrackID);
    
    
    Efrac = 1-(partial_E/total_E);
    
    //completeness
    double totenergy =0;
    for(size_t k = 0; k < all_hits.size(); ++k){
      art::Ptr<recob::Hit> hit = all_hits[k];
      std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
      for(size_t l = 0; l < TrackIDs.size(); ++l){
        if(std::abs(TrackIDs[l].trackID)==TrackID) {
          totenergy += TrackIDs[l].energy;
        }
      }
    } 
    Ecomplet = partial_E/totenergy;


  }
  //========================================================================
  bool NeutrinoShowerEff::insideFV( double vertex[4]){ 

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    //This is temporarily we should define a common FV    
    double x = vertex[0];
    double y = vertex[1];
    double z = vertex[2];

    /*   if( fabs(x) > 350.0 ) return false;
         else if( fabs(y) > 550.0 ) return false;
         else if( z< 0 || z> 400.0 ) return false;
         else return true;
    */

    if (x>fFidVolXmin && x<fFidVolXmax&&
        y>fFidVolYmin && y<fFidVolYmax&&
        z>fFidVolZmin && z<fFidVolZmax)
      return true;
    else
      return false;



  }
  //========================================================================
  void NeutrinoShowerEff::doEfficiencies(){

    art::ServiceHandle<art::TFileService> tfs;

    if(TEfficiency::CheckConsistency(*h_Ev_num,*h_Ev_den)){
      h_Eff_Ev = tfs->make<TEfficiency>(*h_Ev_num,*h_Ev_den);
      TGraphAsymmErrors *grEff_Ev = h_Eff_Ev->CreateGraph();
      grEff_Ev->Write("grEff_Ev");
      h_Eff_Ev->Write("h_Eff_Ev");
    }
   
    if(TEfficiency::CheckConsistency(*h_Ee_num,*h_Ee_den)){
      h_Eff_Ee = tfs->make<TEfficiency>(*h_Ee_num,*h_Ee_den);
      TGraphAsymmErrors *grEff_Ee = h_Eff_Ee->CreateGraph();
      grEff_Ee->Write("grEff_Ee");
      h_Eff_Ee->Write("h_Eff_Ee");
    }
   
    if(TEfficiency::CheckConsistency(*h_Pe_num,*h_Pe_den)){ 
      h_Eff_Pe = tfs->make<TEfficiency>(*h_Pe_num,*h_Pe_den);
      TGraphAsymmErrors *grEff_Pe = h_Eff_Pe->CreateGraph();
      grEff_Pe->Write("grEff_Pe");
      h_Eff_Pe->Write("h_Eff_Pe");
    }
    if(TEfficiency::CheckConsistency(*h_theta_num,*h_theta_den)){
      h_Eff_theta = tfs->make<TEfficiency>(*h_theta_num,*h_theta_den);
      TGraphAsymmErrors *grEff_theta = h_Eff_theta->CreateGraph();
      grEff_theta->Write("grEff_theta");
      h_Eff_theta->Write("h_Eff_theta");
    }
   
  }

  //============================================
  //Check CNN track/shower ID
  //============================================
  template <size_t N> void NeutrinoShowerEff::checkCNNtrkshw(const art::Event& evt, std::vector<art::Ptr<recob::Hit>>all_hits){
    if (fCNNEMModuleLabel=="") return;

    art::ServiceHandle<cheat::BackTracker> bt;
    //auto const* geo = lar::providerFrom<geo::Geometry>();

    auto hitResults = anab::MVAReader<recob::Hit, N>::create(evt, fCNNEMModuleLabel);
    if (hitResults){
      int trkLikeIdx = hitResults->getIndex("track");
      int emLikeIdx = hitResults->getIndex("em");
      if ((trkLikeIdx < 0) || (emLikeIdx < 0)){
        throw cet::exception("NeutrinoShowerEff") << "No em/track labeled columns in MVA data products." << std::endl;
      }
      //std::cout<<all_hits.size()<<std::endl;
      for (size_t i = 0; i<all_hits.size(); ++i){
        //find out if the hit was generated by an EM particle
        bool isEMparticle = false;
        int  pdg = INT_MAX;
        std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(all_hits[i]);
        if (!TrackIDs.size()) continue;
//        raw::ChannelID_t channel = all_hits[i]->Channel();
//        bool firstwire = false;
//        std::vector<geo::WireID> wires = geo->ChannelToWire(channel);
//        for (auto &w : wires){
//          if (w.TPC == all_hits[i]->WireID().TPC){
//            if (w==all_hits[i]->WireID()) firstwire = true;
//            break;
//          }
//        }
        int trkid = INT_MAX;
        double maxE = -1;
        for(size_t k = 0; k < TrackIDs.size(); k++){
          if (TrackIDs[k].energy>maxE){
            maxE = TrackIDs[k].energy;
            trkid = TrackIDs[k].trackID;
          }
        }
        if (trkid!=INT_MAX){
          auto *particle = bt->TrackIDToParticle(trkid);
          if (particle){
            pdg = particle->PdgCode();
            if (std::abs(pdg)==11||//electron/positron
                pdg == 22 ||//photon
                pdg == 111){//pi0
              isEMparticle = true;
            }
          }
        }
        auto vout = hitResults->getOutput(all_hits[i]);        
        //std::cout<<i<<" "<<all_hits[i]->View()<<" "<<vout[0]<<" "<<vout[1]<<" "<<vout[2]<<" "<<vout[3]<<" "<<firstwire<<std::endl;
        double trk_like = -1, trk_or_em = vout[trkLikeIdx] + vout[emLikeIdx];
        if (trk_or_em > 0){
          trk_like = vout[trkLikeIdx] / trk_or_em;
          //std::cout<<"trk_like "<<trk_like<<std::endl;
          if (isEMparticle){
            h_trklike_em->Fill(trk_like);
//            if (trk_like>0.4&&trk_like<0.41){
//              std::cout<<std::string(all_hits[i]->WireID())<<" "<<all_hits[i]->PeakTime()<<std::endl;
//              std::cout<<vout[trkLikeIdx]<<" "<<vout[emLikeIdx]<<" "<<trk_like<<std::endl;
//            }
          }
          else{
            h_trklike_nonem->Fill(trk_like);
          }
        }
      }
    }
    else{
      std::cout<<"Couldn't get hitResults."<<std::endl;
    }
  }


  //========================================================================
  void NeutrinoShowerEff::reset(){

    MC_incoming_PDG = -999;
    MC_lepton_PDG =-999;
    MC_isCC =-999;
    MC_channel =-999;
    MC_target =-999;
    MC_Q2 =-999.0;
    MC_W =-999.0;
    MC_lepton_theta = -999.0;
    MC_leptonID = -999;
    MC_LeptonTrack = -999;
 
    for(int i=0; i<MAX_SHOWERS; i++){
      sh_direction_X[i] = -999.0;
      sh_direction_Y[i] = -999.0;
      sh_direction_Z[i] = -999.0;
      sh_start_X[i] = -999.0;
      sh_start_Y[i] = -999.0;
      sh_start_Z[i] = -999.0;
      sh_bestplane[i] = -999.0;
      sh_length[i] = -999.0;
      sh_hasPrimary_e[i] = -999.0;
      sh_Efrac_contamination[i] = -999.0;
      sh_nHits[i] = -999.0;
      for( int j=0; j<3; j++){
        sh_energy[i][j] = -999.0;
        sh_MIPenergy[i][j] = -999.0;
        sh_dEdx[i][j] = -999.0;
      }
    } 
  }
  //========================================================================
  DEFINE_ART_MODULE(NeutrinoShowerEff)

} 

#endif // NeutrinoShowerEff_Module
