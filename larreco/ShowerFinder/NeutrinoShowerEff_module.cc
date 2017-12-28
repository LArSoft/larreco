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
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

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
    TH1D *h_Ev_num_dEdx;

    TH1D *h_Ee_den;
    TH1D *h_Ee_num;
    TH1D *h_Ee_num_dEdx;

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
    TEfficiency* h_Eff_Ev_dEdx = 0;
    TEfficiency* h_Eff_Ee = 0;  
    TEfficiency* h_Eff_Ee_dEdx = 0;  
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
    
    TH1D *h_dEdXasymm_electronorpositron_NueCC;
    TH1D *h_dEdXasymm_photon_NC;

    TH1D *h_mpi0_electronorpositron_NueCC;
    TH1D *h_mpi0_photon_NC;

    //Study CNN track/shower id
    TH1D *h_trklike_em;
    TH1D *h_trklike_nonem;


    //Study the angle between the reconstructed shower direction w.r.t MC true particle direction
    TH1D *h_CosThetaShDirwrtTrueparticle_electronorpositron_NueCC;
    TH1D *h_CosThetaShDirwrtTrueparticle_electronorpositron_NC;
    TH1D *h_CosThetaShDirwrtTrueparticle_photon_NueCC;
    TH1D *h_CosThetaShDirwrtTrueparticle_photon_NC;
    TH1D *h_CosThetaShDirwrtTrueparticle_proton_NueCC;
    TH1D *h_CosThetaShDirwrtTrueparticle_proton_NC;
    TH1D *h_CosThetaShDirwrtTrueparticle_chargedpion_NueCC;
    TH1D *h_CosThetaShDirwrtTrueparticle_chargedpion_NC;
    
    //Study the reconstructed shower start position (x,y,z) w.r.t MC true particle start position
    TH1D *h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NueCC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NueCC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NueCC;

    TH1D *h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NC;

    TH1D *h_ShStartXwrtTrueparticleStartXDiff_photon_NueCC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_photon_NueCC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_photon_NueCC;

    TH1D *h_ShStartXwrtTrueparticleStartXDiff_photon_NC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_photon_NC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_photon_NC;

    //True photon end position comparison with the reconstructed shower start position
    TH1D *h_ShStartXwrtTrueparticleEndXDiff_photon_NueCC;
    TH1D *h_ShStartYwrtTrueparticleEndYDiff_photon_NueCC;
    TH1D *h_ShStartZwrtTrueparticleEndZDiff_photon_NueCC;

    TH1D *h_ShStartXwrtTrueparticleEndXDiff_photon_NC;
    TH1D *h_ShStartYwrtTrueparticleEndYDiff_photon_NC;
    TH1D *h_ShStartZwrtTrueparticleEndZDiff_photon_NC;

    TH1D *h_ShStartXwrtTrueparticleStartXDiff_proton_NueCC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_proton_NueCC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_proton_NueCC;

    TH1D *h_ShStartXwrtTrueparticleStartXDiff_proton_NC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_proton_NC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_proton_NC;
    
    TH1D *h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NueCC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NueCC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NueCC;
    
    TH1D *h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NC;
    TH1D *h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NC;
    TH1D *h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NC;
    


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
    double sh_purity[MAX_SHOWERS];
    double sh_completeness[MAX_SHOWERS];
    int    sh_nHits[MAX_SHOWERS];
    int    sh_largest;
    int    sh_pdg[MAX_SHOWERS];
    double sh_dEdxasymm[MAX_SHOWERS];
    double sh_mpi0;

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
  
  
    h_Ev_den = tfs->make<TH1D>("h_Ev_den","Neutrino Energy; Neutrino Energy (GeV); Shower Reconstruction Efficiency",20,E_bins);
    h_Ev_den->Sumw2();
    h_Ev_num = tfs->make<TH1D>("h_Ev_num","Neutrino Energy; Neutrino Energy (GeV); Shower Reconstruction Efficiency",20,E_bins);
    h_Ev_num->Sumw2();
    h_Ev_num_dEdx = tfs->make<TH1D>("h_Ev_num_dEdx","Neutrino Energy; Neutrino Energy (GeV); Shower Reconstruction Efficiency",20,E_bins);
    h_Ev_num_dEdx->Sumw2();

    h_Ee_den = tfs->make<TH1D>("h_Ee_den","Electron Energy; Electron Energy (GeV); Shower Reconstruction Efficiency",20,E_bins);
    h_Ee_den->Sumw2();
    h_Ee_num = tfs->make<TH1D>("h_Ee_num","Electron Energy; Electron Energy (GeV); Shower Reconstruction Efficiency",20,E_bins);
    h_Ee_num->Sumw2();
    h_Ee_num_dEdx = tfs->make<TH1D>("h_Ee_num_dEdx","Electron Energy; Electron Energy (GeV); Shower Reconstruction Efficiency",20,E_bins);
    h_Ee_num_dEdx->Sumw2();

    h_Pe_den = tfs->make<TH1D>("h_Pe_den","Electron Momentum; Electron Momentum (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Pe_den->Sumw2();
    h_Pe_num = tfs->make<TH1D>("h_Pe_num","Electron Momentum; Electron Momentum (GeV); Shower reconstruction Efficiency",20,E_bins);
    h_Pe_num->Sumw2();

    h_theta_den = tfs->make<TH1D>("h_theta_den","CosTheta; CosTheta w.r.t beam direction (Degrees); Shower reconstruction Efficiency",43,theta_bin);
    h_theta_den->Sumw2();
    h_theta_num = tfs->make<TH1D>("h_theta_num","CosTheta; CosTheta w.r.t beam direction (Degrees); Shower reconstruction Efficiency",43,theta_bin);
    h_theta_num->Sumw2();

    h_Efrac_shContamination = tfs->make<TH1D>("h_Efrac_shContamination","Efrac Lepton; Energy fraction (contamination);",60,0,1.2);
    h_Efrac_shContamination->Sumw2();
    h_Efrac_shPurity = tfs->make<TH1D>("h_Efrac_shPurity","Efrac Lepton; Energy fraction (Purity);",60,0,1.2);
    h_Efrac_shPurity->Sumw2();
    h_Ecomplet_lepton = tfs->make<TH1D>("h_Ecomplet_lepton","Ecomplet Lepton; Shower Completeness;",60,0,1.2);
    h_Ecomplet_lepton->Sumw2();

    h_HighestHitsProducedParticlePDG_NueCC= tfs->make<TH1D>("h_HighestHitsProducedParticlePDG_NueCC","PDG Code; PDG Code;",4,-0.5,3.5);//0 for undefined, 1=electron, 2=photon, 3=anything else     //Signal
    h_HighestHitsProducedParticlePDG_NueCC->Sumw2();     
    h_HighestHitsProducedParticlePDG_bkg= tfs->make<TH1D>("h_HighestHitsProducedParticlePDG_bkg","PDG Code; PDG Code;",4,-0.5,3.5);//0 for undefined, 1=electron, 2=photon, 3=anything else     //bkg    
    h_HighestHitsProducedParticlePDG_bkg->Sumw2();


    h_Efrac_NueCCPurity= tfs->make<TH1D>("h_Efrac_NueCCPurity","Efrac NueCC; Energy fraction (Purity);",60,0,1.2);     //Signal
    h_Efrac_NueCCPurity->Sumw2();
    h_Ecomplet_NueCC= tfs->make<TH1D>("h_Ecomplet_NueCC","Ecomplet NueCC; Shower Completeness;",60,0,1.2);     
    h_Ecomplet_NueCC->Sumw2();

    
    h_Efrac_bkgPurity= tfs->make<TH1D>("h_Efrac_bkgPurity","Efrac bkg; Energy fraction (Purity);",60,0,1.2);     //Background
    h_Efrac_bkgPurity->Sumw2();
    h_Ecomplet_bkg= tfs->make<TH1D>("h_Ecomplet_bkg","Ecomplet bkg; Shower Completeness;",60,0,1.2);     
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
    
    h_dEdXasymm_electronorpositron_NueCC=tfs->make<TH1D>("h_dEdXasymm_electronorpositron_NueCC","dE/dX asymmetry; Electron or Positron dE/dX asymmetry;",60,0.0,1.2); 
    h_dEdXasymm_photon_NC=tfs->make<TH1D>("h_dEdXasymm_photon_NC","dE/dX asymmetry; photon dE/dx asymmetry;",60,0.0,1.2); 

    h_mpi0_electronorpositron_NueCC=tfs->make<TH1D>("h_mpi0_electronorpositron_NueCC","m(#gamma#gamma); Electron or Positron dE/dX (MeV/cm);",100,0,1); 
    h_mpi0_photon_NC=tfs->make<TH1D>("h_mpi0_photon_NC","m(#gamma#gamma); Electron or Positron dE/dX (MeV/cm);",100,0,1);

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
    
    h_dEdXasymm_electronorpositron_NueCC->Sumw2();
    h_dEdXasymm_photon_NC->Sumw2();

    h_mpi0_electronorpositron_NueCC->Sumw2();
    h_mpi0_photon_NC->Sumw2();

    h_trklike_em = tfs->make<TH1D>("h_trklike_em","EM hits; Track-like Score;",100,0,1);
    h_trklike_nonem = tfs->make<TH1D>("h_trklike_nonem","Non-EM hits; Track-like Score;",100,0,1);
    

    //Study the constheta angle between the reconstructed shower direction w.r.t MC true particle direction
    h_CosThetaShDirwrtTrueparticle_electronorpositron_NueCC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_electronorpositron_NueCC","CosTheta; cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_electronorpositron_NC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_electronorpositron_NC","CosTheta;cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_photon_NueCC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_photon_NueCC","CosTheta;cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_photon_NC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_photon_NC","CosTheta;cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_proton_NueCC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_proton_NueCC","CosTheta;cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_proton_NC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_proton_NC","CosTheta;cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_chargedpion_NueCC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_chargedpion_NueCC","CosTheta;cos#theta;",110,-1.1,1.1); 
    h_CosThetaShDirwrtTrueparticle_chargedpion_NC=tfs->make<TH1D>("h_CosThetaShDirwrtTrueparticle_chargedpion_NC","CosTheta;cos#theta;",110,-1.1,1.1); 
    
    //Study the reconstructed shower start position (x,y,z) w.r.t MC true particle start position
    h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NueCC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NueCC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NueCC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NueCC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NueCC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NueCC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 

    h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 


    h_ShStartXwrtTrueparticleStartXDiff_photon_NueCC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_photon_NueCC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_photon_NueCC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_photon_NueCC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_photon_NueCC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_photon_NueCC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 

    h_ShStartXwrtTrueparticleStartXDiff_photon_NC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_photon_NC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_photon_NC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_photon_NC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_photon_NC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_photon_NC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 




    h_ShStartXwrtTrueparticleEndXDiff_photon_NueCC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleEndXDiff_photon_NueCC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleEndYDiff_photon_NueCC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleEndYDiff_photon_NueCC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleEndZDiff_photon_NueCC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleEndZDiff_photon_NueCC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 

    h_ShStartXwrtTrueparticleEndXDiff_photon_NC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleEndXDiff_photon_NC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleEndYDiff_photon_NC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleEndYDiff_photon_NC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleEndZDiff_photon_NC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleEndZDiff_photon_NC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0);


    h_ShStartXwrtTrueparticleStartXDiff_proton_NueCC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_proton_NueCC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_proton_NueCC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_proton_NueCC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_proton_NueCC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_proton_NueCC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 

    h_ShStartXwrtTrueparticleStartXDiff_proton_NC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_proton_NC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_proton_NC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_proton_NC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_proton_NC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_proton_NC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 



    h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NueCC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NueCC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NueCC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NueCC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NueCC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NueCC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 

    h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NC=tfs->make<TH1D>("h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NC","ShVx-TrueParticleVx; ShVx-TrueParticleVx (cm);",100,-5.0,5.0); 
    h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NC=tfs->make<TH1D>("h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NC","ShVy-TrueParticleVy; ShVy-TrueParticleVy (cm);",100,-5.0,5.0); 
    h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NC=tfs->make<TH1D>("h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NC","ShVz-TrueParticleVz; ShVz-TrueParticleVz (cm);",100,-5.0,5.0); 


    //Study the constheta angle between the reconstructed shower direction w.r.t MC true particle direction
    h_CosThetaShDirwrtTrueparticle_electronorpositron_NueCC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_electronorpositron_NC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_photon_NueCC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_photon_NC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_proton_NueCC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_proton_NC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_chargedpion_NueCC->Sumw2();
    h_CosThetaShDirwrtTrueparticle_chargedpion_NC->Sumw2();
    
    //Study the reconstructed shower start position (x,y,z) w.r.t MC true particle start position
    h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NueCC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NueCC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NueCC->Sumw2();

    h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NC->Sumw2();

    h_ShStartXwrtTrueparticleStartXDiff_photon_NueCC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_photon_NueCC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_photon_NueCC->Sumw2();

    h_ShStartXwrtTrueparticleStartXDiff_photon_NC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_photon_NC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_photon_NC->Sumw2();


    h_ShStartXwrtTrueparticleEndXDiff_photon_NueCC->Sumw2();
    h_ShStartYwrtTrueparticleEndYDiff_photon_NueCC->Sumw2();
    h_ShStartZwrtTrueparticleEndZDiff_photon_NueCC->Sumw2();

    h_ShStartXwrtTrueparticleEndXDiff_photon_NC->Sumw2();
    h_ShStartYwrtTrueparticleEndYDiff_photon_NC->Sumw2();
    h_ShStartZwrtTrueparticleEndZDiff_photon_NC->Sumw2();


    h_ShStartXwrtTrueparticleStartXDiff_proton_NueCC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_proton_NueCC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_proton_NueCC->Sumw2();

    h_ShStartXwrtTrueparticleStartXDiff_proton_NC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_proton_NC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_proton_NC->Sumw2();
    
    h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NueCC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NueCC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NueCC->Sumw2();
    
    h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NC->Sumw2();
    h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NC->Sumw2();
    h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NC->Sumw2();



    



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
      fEventTree->Branch("sh_purity",&sh_purity,"sh_purity[n_showers]/D");
      fEventTree->Branch("sh_completeness",&sh_completeness,"sh_completeness[n_showers]/D");
      fEventTree->Branch("sh_Efrac_best", &sh_Efrac_best, "sh_Efrac_best/D");
      fEventTree->Branch("sh_nHits",&sh_nHits, "sh_nHits[n_showers]/I");
      fEventTree->Branch("sh_largest",&sh_largest,"sh_largest/I");
      fEventTree->Branch("sh_pdg",&sh_pdg,"sh_pdg[n_showers]/I");
      fEventTree->Branch("sh_dEdxasymm", &sh_dEdxasymm, "sh_dEdxasymm[n_showers]/D");
      fEventTree->Branch("sh_mpi0",&sh_mpi0,"sh_mpi0/D");
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
      if(isFiducial) fEventTree->Fill();
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
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
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
    int ParticlePDG_HighestShHits=0;//undefined
    int shower_bestplane=0;
    double Showerparticlededx_inbestplane=0.0;
    double dEdxasymm_largestshw = -1.;

    int showerPDGwithHighestHitsforFillingdEdX=0;//0=undefined,1=electronorpositronshower,2=photonshower,3=protonshower,4=neutronshower,5=chargedpionshower,6=neutralpionshower,7=everythingelseshower

    
    double ShAngle=-9999.0,ShVxTrueParticleVxDiff=-9999.0,ShVyTrueParticleVyDiff=-9999.0,ShVzTrueParticleVzDiff=-9999.0, ShStartVxTrueParticleEndVxDiff=-9999.0,ShStartVyTrueParticleEndVyDiff=-9999.0,ShStartVzTrueParticleEndVzDiff=-9999.0;

    const simb::MCParticle *MClepton_reco = NULL; 
    int nHits =0;
    
    TVector3 p1, p2;
    double E1st = 0;
    double E2nd = 0;

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

      double dEdxasymm = -1;
      double dEdx0 = 0;
      if (shower->best_plane()>=0&&shower->best_plane()<int(shower->dEdx().size())){
        dEdx0 = shower->dEdx()[shower->best_plane()];
      }
      double dEdx1 = 0;
      double maxE = 0;
      for (int j = 0; j<3; ++j){
        if (j==shower->best_plane()) continue;
        if (j>=int(shower->Energy().size())) continue;
        if (shower->Energy()[j]>maxE){
          maxE = shower->Energy()[j];
          dEdx1 = shower->dEdx()[j];
        }
      }
      if (dEdx0||dEdx1){
        dEdxasymm = std::abs(dEdx0-dEdx1)/(dEdx0+dEdx1);
      }
      sh_dEdxasymm[i] = dEdxasymm;

      if (shower->best_plane()>=0 && shower->best_plane()<int(shower->Energy().size())){
        if (shower->Energy()[shower->best_plane()]>E1st){
          if (p1.Mag()){
            E2nd = E1st;
            p2 = p1;
          }
          E1st = shower->Energy()[shower->best_plane()];
          p1[0] = E1st * shower->Direction().X();
          p1[1] = E1st * shower->Direction().Y();
          p1[2] = E1st * shower->Direction().Z();
        }
        else{
          if (shower->Energy()[shower->best_plane()]>E2nd){
            E2nd = shower->Energy()[shower->best_plane()];
            p2[0] = E2nd * shower->Direction().X();
            p2[1] = E2nd * shower->Direction().Y();
            p2[2] = E2nd * shower->Direction().Z();
          }
        }
      }

      std::vector<art::Ptr<recob::Hit>> sh_hits = sh_hitsAll.at(i);  

      if (!sh_hits.size()){
        //no shower hits found, try pfparticle
        // PFParticles
        art::Handle<std::vector<recob::PFParticle> > pfpHandle;
        std::vector<art::Ptr<recob::PFParticle> > pfps;
        if (event.getByLabel(fShowerModuleLabel, pfpHandle))
          art::fill_ptr_vector(pfps, pfpHandle);
        // Clusters
        art::Handle<std::vector<recob::Cluster> > clusterHandle;
        std::vector<art::Ptr<recob::Cluster> > clusters;
        if (event.getByLabel(fShowerModuleLabel, clusterHandle))
          art::fill_ptr_vector(clusters, clusterHandle);
        art::FindManyP<recob::PFParticle> fmps(showerHandle, event, fShowerModuleLabel);
        art::FindManyP<recob::Cluster> fmcp(pfpHandle, event, fShowerModuleLabel);
        art::FindManyP<recob::Hit> fmhc(clusterHandle, event, fShowerModuleLabel);
        if (fmps.isValid()){
          std::vector<art::Ptr<recob::PFParticle>> pfs = fmps.at(i);
          for (size_t ipf = 0; ipf<pfs.size(); ++ipf){
            if (fmcp.isValid()){
              std::vector<art::Ptr<recob::Cluster>> clus = fmcp.at(pfs[ipf].key());
              for (size_t iclu = 0; iclu<clus.size(); ++iclu){
                if (fmhc.isValid()){
                  std::vector<art::Ptr<recob::Hit>> hits = fmhc.at(clus[iclu].key());
                  for (size_t ihit = 0; ihit<hits.size(); ++ihit){
                    sh_hits.push_back(hits[ihit]);
                  }
                }
              }
            }
          }
        }
      }
      //  std::cout<<" shower best plane:"<<shower->best_plane()<<" shower dEdx size:"<<shower->dEdx().size()<<std::endl;
      //for( size_t j =0; j<shower->dEdx().size(); j++) std::cout<<shower->dEdx()[j]<<" ";

      const simb::MCParticle *particle;
      double tmpEfrac_contamination = 0.0;  //fraction of non EM energy contatiminatio (see truthMatcher for definition)              
      double tmpEcomplet =0;
  
      int tmp_nHits = sh_hits.size();
     

      truthMatcher( all_hits, sh_hits, particle, tmpEfrac_contamination,tmpEcomplet);
      //truthMatcher( all_hits, sh_hits, particle, tmpEfrac_contaminationNueCC,tmpEcompletNueCC );
      if (!particle) continue;

      sh_Efrac_contamination[i] = tmpEfrac_contamination;
      sh_purity[i] = 1 - tmpEfrac_contamination;
      sh_completeness[i] = tmpEcomplet;
      sh_nHits[i] = tmp_nHits; 
      sh_hasPrimary_e[i] = 0;
      sh_pdg[i] = particle->PdgCode();

      //Shower with highest hits       
      if( tmp_nHits > nHits ){
        sh_largest = i;
        dEdxasymm_largestshw = dEdxasymm;
        nHits = tmp_nHits;
        Ecomplet_NueCC =tmpEcomplet;
        Efrac_contaminationNueCC = tmpEfrac_contamination; 
	//Calculate Shower anagle w.r.t True particle
	double ShDirMag  = sqrt(pow(sh_direction_X[i],2)+pow(sh_direction_Y[i],2)+pow(sh_direction_Z[i],2));  
	ShAngle = (sh_direction_X[i]*particle->Px() + sh_direction_Y[i]*particle->Py() +sh_direction_Z[i]*particle->Pz())/(ShDirMag*particle->P()) ;


	ShVxTrueParticleVxDiff=sh_start_X[i]-particle->Vx();
      	ShVyTrueParticleVyDiff=sh_start_Y[i]-particle->Vy();
	ShVzTrueParticleVzDiff=sh_start_Z[i]-particle->Vz();

	//put overflow and underflow at top and bottom bins:
	if (ShVxTrueParticleVxDiff > 5) ShVxTrueParticleVxDiff = 4.99;
	else if (ShVxTrueParticleVxDiff < -5) ShVxTrueParticleVxDiff = -5;
	if (ShVyTrueParticleVyDiff > 5) ShVyTrueParticleVyDiff = 4.99;
	else if (ShVyTrueParticleVyDiff < -5) ShVyTrueParticleVyDiff = -5;
	if (ShVzTrueParticleVzDiff > 5) ShVzTrueParticleVzDiff = 4.99;
	else if (ShVzTrueParticleVzDiff < -5) ShVzTrueParticleVzDiff = -5;


	ShStartVxTrueParticleEndVxDiff=sh_start_X[i]-particle->EndX();
	ShStartVyTrueParticleEndVyDiff=sh_start_Y[i]-particle->EndY();
	ShStartVzTrueParticleEndVzDiff=sh_start_Z[i]-particle->EndZ();
 

        if(std::abs(particle->PdgCode())==11){
          ParticlePDG_HighestShHits=1;
        }else if(particle->PdgCode()==22){
          ParticlePDG_HighestShHits=2;
        }else{
          ParticlePDG_HighestShHits=3;
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

      if( std::abs(particle->PdgCode())  == fLeptonPDGcode && particle->TrackId() == MC_leptonID ){

        if(tmpEcomplet>Ecomplet_lepton){

          Ecomplet_lepton = tmpEcomplet;
	   
          Efrac_contamination = tmpEfrac_contamination;
          MClepton_reco = particle;
          sh_Efrac_best =Efrac_contamination; 
           
        }
      }          
    }//end of looping all the showers

    if (p1.Mag()&&p2.Mag()){
      sh_mpi0 = sqrt(pow(p1.Mag()+p2.Mag(),2)-(p1+p2).Mag2());
    }
    else sh_mpi0 = 0;
   
    if( MClepton_reco && MClepton  ){
      if( MC_isCC && (fNeutrinoPDGcode == std::abs(MC_incoming_PDG)) ){ 
        h_Efrac_shContamination->Fill(Efrac_contamination);
        h_Efrac_shPurity->Fill(1-Efrac_contamination);
	h_Ecomplet_lepton->Fill(Ecomplet_lepton);
	
	// Selecting good showers requires completeness of gretaer than 70 % and Purity > 70 %	
	if( Efrac_contamination < fMaxEfrac && Ecomplet_lepton> fMinCompleteness ){
        
	  h_Pe_num->Fill(Pe);
          h_Ev_num->Fill(MC_incoming_P[3]);
	  h_Ee_num->Fill(MC_lepton_startMomentum[3]);
	  h_theta_num->Fill(theta_e);

	  if (Showerparticlededx_inbestplane > 1 && Showerparticlededx_inbestplane < 3) {
	    h_Ev_num_dEdx->Fill(MC_incoming_P[3]);
	    h_Ee_num_dEdx->Fill(MC_lepton_startMomentum[3]);
	  }
	}
      }
    }

    //NueCC SIgnal and background Completeness
    if(MC_isCC==1
       &&(fNeutrinoPDGcode == std::abs(MC_incoming_PDG))
       &&isFiducial){
        h_HighestHitsProducedParticlePDG_NueCC->Fill(ParticlePDG_HighestShHits);
	
        if(ParticlePDG_HighestShHits>0){// atleat one shower is reconstructed
          h_Ecomplet_NueCC->Fill(Ecomplet_NueCC);
          h_Efrac_NueCCPurity->Fill(1-Efrac_contaminationNueCC);    
	  
          h_esh_bestplane_NueCC->Fill(shower_bestplane);
          if(showerPDGwithHighestHitsforFillingdEdX==1)//electron or positron shower
            {
              h_dEdX_electronorpositron_NueCC->Fill(Showerparticlededx_inbestplane);
	      //Study the angle between the reconstructed shower direction w.r.t MC true particle direction
	      h_CosThetaShDirwrtTrueparticle_electronorpositron_NueCC->Fill(ShAngle);
	       
	      //Study the reconstructed shower start position (x,y,z) w.r.t MC true particle start position
	      h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NueCC->Fill(ShVxTrueParticleVxDiff);
	      h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NueCC->Fill(ShVyTrueParticleVyDiff);
	      h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NueCC->Fill(ShVzTrueParticleVzDiff);

	      h_dEdXasymm_electronorpositron_NueCC->Fill(dEdxasymm_largestshw);

              h_mpi0_electronorpositron_NueCC->Fill(sh_mpi0);

	    }else if(showerPDGwithHighestHitsforFillingdEdX==2)//photon shower
            {
              h_dEdX_photon_NueCC->Fill(Showerparticlededx_inbestplane);
	      h_CosThetaShDirwrtTrueparticle_photon_NueCC->Fill(ShAngle);
	      h_ShStartXwrtTrueparticleStartXDiff_photon_NueCC->Fill(ShVxTrueParticleVxDiff);
	      h_ShStartYwrtTrueparticleStartYDiff_photon_NueCC->Fill(ShVyTrueParticleVyDiff);
	      h_ShStartZwrtTrueparticleStartZDiff_photon_NueCC->Fill(ShVzTrueParticleVzDiff);


	      h_ShStartXwrtTrueparticleEndXDiff_photon_NueCC->Fill(ShStartVxTrueParticleEndVxDiff);
	      h_ShStartYwrtTrueparticleEndYDiff_photon_NueCC->Fill(ShStartVyTrueParticleEndVyDiff);
	      h_ShStartZwrtTrueparticleEndZDiff_photon_NueCC->Fill(ShStartVzTrueParticleEndVzDiff);
 

	      
	    }else if(showerPDGwithHighestHitsforFillingdEdX==3)//proton shower
            {
              h_dEdX_proton_NueCC->Fill(Showerparticlededx_inbestplane);
	      h_CosThetaShDirwrtTrueparticle_proton_NueCC->Fill(ShAngle);

	      h_ShStartXwrtTrueparticleStartXDiff_proton_NueCC->Fill(ShVxTrueParticleVxDiff);
	      h_ShStartYwrtTrueparticleStartYDiff_proton_NueCC->Fill(ShVyTrueParticleVyDiff);
	      h_ShStartZwrtTrueparticleStartZDiff_proton_NueCC->Fill(ShVzTrueParticleVzDiff);


   	    }else if(showerPDGwithHighestHitsforFillingdEdX==4)//neutron shower
            {
              h_dEdX_neutron_NueCC->Fill(Showerparticlededx_inbestplane);
            
	    }else if(showerPDGwithHighestHitsforFillingdEdX==5)//charged pion shower
            {
              h_dEdX_chargedpion_NueCC->Fill(Showerparticlededx_inbestplane);
	      h_CosThetaShDirwrtTrueparticle_chargedpion_NueCC->Fill(ShAngle);
	      h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NueCC->Fill(ShVxTrueParticleVxDiff);
	      h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NueCC->Fill(ShVyTrueParticleVyDiff);
	      h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NueCC->Fill(ShVzTrueParticleVzDiff);
   
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
      h_HighestHitsProducedParticlePDG_bkg->Fill(ParticlePDG_HighestShHits);
      
      
      if(ParticlePDG_HighestShHits>0){
	h_Ecomplet_bkg->Fill(Ecomplet_NueCC);
	h_Efrac_bkgPurity->Fill(1-Efrac_contaminationNueCC);	
        
        
	h_esh_bestplane_NC->Fill(shower_bestplane);
        if(showerPDGwithHighestHitsforFillingdEdX==1)//electron or positron shower
          {
            h_dEdX_electronorpositron_NC->Fill(Showerparticlededx_inbestplane);
	    h_CosThetaShDirwrtTrueparticle_electronorpositron_NC->Fill(ShAngle);


	    h_ShStartXwrtTrueparticleStartXDiff_electronorpositron_NC->Fill(ShVxTrueParticleVxDiff);
	    h_ShStartYwrtTrueparticleStartYDiff_electronorpositron_NC->Fill(ShVyTrueParticleVyDiff);
	    h_ShStartZwrtTrueparticleStartZDiff_electronorpositron_NC->Fill(ShVzTrueParticleVzDiff);

	    
	  }else if(showerPDGwithHighestHitsforFillingdEdX==2)//photon shower
          {
            h_dEdX_photon_NC->Fill(Showerparticlededx_inbestplane);
	    h_CosThetaShDirwrtTrueparticle_photon_NC->Fill(ShAngle);
	    

	    h_ShStartXwrtTrueparticleStartXDiff_photon_NC->Fill(ShVxTrueParticleVxDiff);
	    h_ShStartYwrtTrueparticleStartYDiff_photon_NC->Fill(ShVyTrueParticleVyDiff);
	    h_ShStartZwrtTrueparticleStartZDiff_photon_NC->Fill(ShVzTrueParticleVzDiff);
	    
	    h_ShStartXwrtTrueparticleEndXDiff_photon_NC->Fill(ShStartVxTrueParticleEndVxDiff);
	    h_ShStartYwrtTrueparticleEndYDiff_photon_NC->Fill(ShStartVyTrueParticleEndVyDiff);
	    h_ShStartZwrtTrueparticleEndZDiff_photon_NC->Fill(ShStartVzTrueParticleEndVzDiff);

            h_dEdXasymm_photon_NC->Fill(dEdxasymm_largestshw);

            h_mpi0_photon_NC->Fill(sh_mpi0);

	  }else if(showerPDGwithHighestHitsforFillingdEdX==3)//proton shower
          {
            h_dEdX_proton_NC->Fill(Showerparticlededx_inbestplane);
            h_CosThetaShDirwrtTrueparticle_proton_NC->Fill(ShAngle);

	    h_ShStartXwrtTrueparticleStartXDiff_proton_NC->Fill(ShVxTrueParticleVxDiff);
	    h_ShStartYwrtTrueparticleStartYDiff_proton_NC->Fill(ShVyTrueParticleVyDiff);
	    h_ShStartZwrtTrueparticleStartZDiff_proton_NC->Fill(ShVzTrueParticleVzDiff);
    
    
	    
	    
	  }else if(showerPDGwithHighestHitsforFillingdEdX==4)//neutron shower
          {
            h_dEdX_neutron_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==5)//charged pion shower
          {
            h_dEdX_chargedpion_NC->Fill(Showerparticlededx_inbestplane);
	    h_CosThetaShDirwrtTrueparticle_chargedpion_NC->Fill(ShAngle);
	    h_ShStartXwrtTrueparticleStartXDiff_chargedpion_NC->Fill(ShVxTrueParticleVxDiff);
	    h_ShStartYwrtTrueparticleStartYDiff_chargedpion_NC->Fill(ShVyTrueParticleVyDiff);
	    h_ShStartZwrtTrueparticleStartZDiff_chargedpion_NC->Fill(ShVzTrueParticleVzDiff);


	  }else if(showerPDGwithHighestHitsforFillingdEdX==6)//neutral pion shower
          {
            h_dEdX_neutralpion_NC->Fill(Showerparticlededx_inbestplane);
          }else if(showerPDGwithHighestHitsforFillingdEdX==7)//everythingelse shower
          {
            h_dEdX_everythingelse_NC->Fill(Showerparticlededx_inbestplane);
          }
      }//if(ParticlePDG_HighestShHits>0)
    }//else if(!MC_isCC&&isFiducial)

    checkCNNtrkshw<4>(event, all_hits);
  }

  //========================================================================
  void NeutrinoShowerEff::truthMatcher(std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

    MCparticle=0;
    Efrac=1.0;
    Ecomplet=0;
    
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < shower_hits.size(); ++j){
      art::Ptr<recob::Hit> hit = shower_hits[j];
      //For know let's use collection plane to look at the shower reconstruction
      //if( hit->View() != 2) continue;
      std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
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
      // const simb::MCParticle *particle = pi_serv->TrackIDToParticle(ID);
      //if( abs(particle->PdgCode()) == 211 || particle->PdgCode() == 2212 ){
      //if( particle->PdgCode() != 22 && abs(particle->PdgCode()) != 11){
      //noEM_E += ii->second;
      //}

    } 
    MCparticle = pi_serv->TrackIdToParticle_P(TrackID);
    
    
    Efrac = 1-(partial_E/total_E);
    
    //completeness
    double totenergy =0;
    for(size_t k = 0; k < all_hits.size(); ++k){
      art::Ptr<recob::Hit> hit = all_hits[k];
      std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
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

    if(TEfficiency::CheckConsistency(*h_Ev_num_dEdx,*h_Ev_den)){
      h_Eff_Ev_dEdx = tfs->make<TEfficiency>(*h_Ev_num_dEdx,*h_Ev_den);
      TGraphAsymmErrors *grEff_Ev_dEdx = h_Eff_Ev_dEdx->CreateGraph();
      grEff_Ev_dEdx->Write("grEff_Ev_dEdx");
      h_Eff_Ev_dEdx->Write("h_Eff_Ev_dEdx");
    }
   
    if(TEfficiency::CheckConsistency(*h_Ee_num,*h_Ee_den)){
      h_Eff_Ee = tfs->make<TEfficiency>(*h_Ee_num,*h_Ee_den);
      TGraphAsymmErrors *grEff_Ee = h_Eff_Ee->CreateGraph();
      grEff_Ee->Write("grEff_Ee");
      h_Eff_Ee->Write("h_Eff_Ee");
    }

    if(TEfficiency::CheckConsistency(*h_Ee_num_dEdx,*h_Ee_den)){
      h_Eff_Ee_dEdx = tfs->make<TEfficiency>(*h_Ee_num_dEdx,*h_Ee_den);
      TGraphAsymmErrors *grEff_Ee_dEdx = h_Eff_Ee_dEdx->CreateGraph();
      grEff_Ee_dEdx->Write("grEff_Ee_dEdx");
      h_Eff_Ee_dEdx->Write("h_Eff_Ee_dEdx");
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

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
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
        std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(all_hits[i]);
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
          auto *particle = pi_serv->TrackIdToParticle_P(trkid);
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
      sh_purity[i] = -999.0;
      sh_completeness[i] = -999.0;
      sh_nHits[i] = -999.0;
      for( int j=0; j<3; j++){
        sh_energy[i][j] = -999.0;
        sh_MIPenergy[i][j] = -999.0;
        sh_dEdx[i][j] = -999.0;
      }
      sh_pdg[i] = -999;
      sh_dEdxasymm[i] = -999;
    }
    sh_largest = -999;
    sh_mpi0 = -999;
  }
  //========================================================================
  DEFINE_ART_MODULE(NeutrinoShowerEff)

} 

#endif // NeutrinoShowerEff_Module
