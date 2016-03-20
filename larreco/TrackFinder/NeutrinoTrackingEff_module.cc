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
//using namespace std;

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
    void initOutput();

    void processEff(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac);
    bool insideFV(double vertex[4]);
    void doEfficiencies();
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fOutFileName;
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    int         fNeutrinoPDGcode;
    int		fLeptonPDGcode;
    double      fMaxNeutrinoE;
    bool	fSaveMCTree; 

    TFile *fOutFile;
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


 }; // class NeutrinoTrackingEff


//========================================================================
NeutrinoTrackingEff::NeutrinoTrackingEff(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
    initOutput();
}
//========================================================================
NeutrinoTrackingEff::~NeutrinoTrackingEff(){
  //destructor
}
//========================================================================
void NeutrinoTrackingEff::reconfigure(fhicl::ParameterSet const& p){

    fOutFileName         = p.get<std::string>("outFile");
    fMCTruthModuleLabel  = p.get<std::string>("MCTruthModuleLabel");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fLeptonPDGcode       = p.get<int>("LeptonPDGcode");
    fNeutrinoPDGcode     = p.get<int>("NeutrinoPDGcode");
    fMaxNeutrinoE	 = p.get<double>("MaxNeutrinoE");
    fSaveMCTree		 = p.get<bool>("SaveMCTree");
    fFidVolCutX          = p.get<float>("FidVolCutX");
    fFidVolCutY          = p.get<float>("FidVolCutY");
    fFidVolCutZ          = p.get<float>("FidVolCutZ");
}
//========================================================================
void NeutrinoTrackingEff::initOutput(){
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");
    double E_bins[21] ={0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,14.0,17.0,20.0,25.0};
    double theta_bin[44]= { 0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,55.,60.,65.,70.,75.,80.,85.,90.};
    double Pbins[18] ={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0};

    TDirectory* subDir = fOutFile->mkdir("Histograms");
    subDir->cd();
    h_Ev_den = new TH1D("h_Ev_den","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
    h_Ev_num = new TH1D("h_Ev_num","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
    h_Pmu_den = new TH1D("h_Pmu_den","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,E_bins);
    h_Pmu_num = new TH1D("h_Pmu_num","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,E_bins);
    h_theta_den = new TH1D("h_theta_den","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
    h_theta_num = new TH1D("h_theta_num","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
    h_Pproton_den = new TH1D("h_Pproton_den","Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
    h_Pproton_num = new TH1D("h_Pproton_num","Protons; Proton Momentum (GeV); Tracking Efficiency", 17, Pbins);
    h_Ppion_plus_den = new TH1D("h_Ppion_plus_den", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ppion_plus_num = new TH1D("h_Ppion_plus_num", "Pions Plus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ppion_minus_den = new TH1D("h_Ppion_minus_den", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);
    h_Ppion_minus_num = new TH1D("h_Ppion_minus_num", "Pions Minus; Pion Momentum (GeV);  Tracking Efficiency", 17, Pbins);

    h_Efrac_lepton = new TH1D("h_Efrac_lepton","Efrac Lepton; Track Energy fraction;",60,0,1.2);
    h_Efrac_proton = new TH1D("h_Efrac_proton","Efrac Proton; Track Energy fraction;",60,0,1.2);
    h_Efrac_pion_plus = new TH1D("h_Efrac_pion_plus","Efrac Pion +; Track Energy fraction;",60,0,1.2);
    h_Efrac_pion_minus = new TH1D("h_Efrac_pion_minus","Efrac Pion -; Track Energy fraction;",60,0,1.2);
    h_trackRes_lepton = new TH1D("h_trackRes_lepton", "Residual; Truth lenght - Reco length (cm);",200,-100,100);
    h_trackRes_proton = new TH1D("h_trackRes_proton", "Proton Residual; Truth lenght - Reco length (cm);",200,-100,100);
    h_trackRes_pion_plus = new TH1D("h_trackRes_pion_plus", "Pion + Residual; Truth lenght - Reco length (cm);",200,-100,100);
    h_trackRes_pion_minus = new TH1D("h_trackRes_pion_minus", "Pion - Residual; Truth lenght - Reco length (cm);",200,-100,100);
 
    if( fSaveMCTree ){
      TDirectory* subDirTree = fOutFile->mkdir("Events");
      subDirTree->cd();
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
      fEventTree->Branch("mc_leptonID", &MC_leptonID);
      fEventTree->Branch("mc_leadingProtonID", &MC_leading_protonID);
      fEventTree->Branch("mc_leadingProtonID", &MC_leading_ProtonP);
      fEventTree->Branch("mc_leadingPionPlusID", &MC_leading_PionPlusID);
      fEventTree->Branch("mc_leadingPionPlusP", &MC_leading_PionPlusP);
      fEventTree->Branch("mc_leadingPionMinusID", &MC_leading_PionMinusID);
      fEventTree->Branch("mc_leadingPionMinusP", &MC_leading_PionMinusP);
      fEventTree->Branch("mc_leptonTrack", &MC_LeptonTrack);
      fEventTree->Branch("mc_protonTrack", &MC_ProtonTrack);
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

      gDirectory = tmpDir;
    }

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

}
//========================================================================
void NeutrinoTrackingEff::endJob(){

    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Histograms");
    h_Efrac_lepton->Write();
    h_Efrac_proton->Write(); 
    h_Efrac_pion_plus->Write();
    h_Efrac_pion_minus->Write();
    h_trackRes_lepton->Write();
    h_trackRes_proton->Write();
    h_trackRes_pion_plus->Write();
    h_trackRes_pion_minus->Write();

    h_Ev_den->Write();
    h_Ev_num->Write();
    h_theta_den->Write();
    h_theta_num->Write();
    h_Pproton_den->Write();
    h_Pproton_num->Write();
    h_Ppion_plus_num->Write();
    h_Ppion_minus_num->Write();
    h_Ppion_plus_den->Write();
    h_Ppion_minus_den->Write();
     
    doEfficiencies();

    if( fSaveMCTree ){
      fOutFile->cd("/Events");
      fEventTree->Write();
      gDirectory = tmpDir;
    }
    fOutFile->Close();
}
//========================================================================
void NeutrinoTrackingEff::beginRun(const art::Run& /*run*/){
  mf::LogInfo("NeutrinoTrackingEff")<<"begin run..."<<std::endl;
}
//========================================================================
void NeutrinoTrackingEff::analyze( const art::Event& event ){

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
         std::cout<<"Incoming E "<<MC_incoming_P[3]<<" is CC? "<<MC_isCC<<" nuPDG "<<MC_incoming_PDG<<" target "<<MC_target<<" vtx "<<MC_vertex[0]<<" "<<MC_vertex[1]<<" "<<MC_vertex[2]<<" "<<MC_vertex[3]<<std::endl;
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
  
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index
    MC_Ntrack = plist.size();

    double proton_length =0.0;
    double pion_plus_length =0.0;
    double pion_minus_length =0.0;
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
             const TLorentzVector& protonEnd   = particle->EndPosition();
             const TLorentzVector& protonStart = particle->Position();
             proton_length = sqrt( pow((protonEnd.X()-protonStart.X()),2)+pow((protonEnd.Y()-protonStart.Y()),2)+pow((protonEnd.Z()-protonStart.Z()),2)); 
             MCproton = particle;
           } 
         }
         else if( particle->PdgCode() == 211 ){
           if(particle->Momentum().E() > tmp_leadingPionPlusE){
             tmp_leadingPionPlusE = particle->Momentum().E();
             MC_leading_PionPlusID = particle->TrackId();          
             MC_leading_PionPlusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             const TLorentzVector& pionEnd = particle->EndPosition();
             const TLorentzVector& pionStart = particle->Position();
             pion_plus_length = sqrt( pow((pionEnd.X()-pionStart.X()),2)+pow((pionEnd.Y()-pionStart.Y()),2)+pow((pionEnd.Z()-pionStart.Z()),2));
             MCpion_plus = particle;
           } 
         }
         else if( particle->PdgCode() == -211 ){
           if(particle->Momentum().E() > tmp_leadingPionMinusE){
             tmp_leadingPionMinusE = particle->Momentum().E();
             MC_leading_PionMinusID = particle->TrackId();          
             MC_leading_PionMinusP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
             const TLorentzVector& pionEnd = particle->EndPosition();
             const TLorentzVector& pionStart = particle->Position();
             pion_minus_length = sqrt( pow((pionEnd.X()-pionStart.X()),2)+pow((pionEnd.Y()-pionStart.Y()),2)+pow((pionEnd.Z()-pionStart.Z()),2));
             MCpion_minus = particle;
           } 
         }
         i++; //paticle index
       }
    } 
    //Saving denominator histograms for lepton pions and protons 
    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
    double Pmu = sqrt(pow(MC_lepton_startMomentum[0],2)+pow(MC_lepton_startMomentum[1],2)+pow(MC_lepton_startMomentum[2],2));
    double Pv  = sqrt(pow(MC_incoming_P[0],2)+pow(MC_incoming_P[1],2)+pow(MC_incoming_P[2],2));
    double theta_mu = acos((MC_incoming_P[0]*MC_lepton_startMomentum[0] + MC_incoming_P[1]*MC_lepton_startMomentum[1] +MC_incoming_P[2]*MC_lepton_startMomentum[2])/(Pv*Pmu) );
    theta_mu *= (180.0/3.14159);
    double truth_lengthLepton =sqrt(pow((MC_lepton_endXYZT[0]-MC_lepton_startXYZT[0]),2) + pow((MC_lepton_endXYZT[1]-MC_lepton_startXYZT[1]),2) + pow((MC_lepton_endXYZT[2]-MC_lepton_startXYZT[2]),2));

    //save CC events within the fiducial volume with the favorite neutrino flavor 
    if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
       if( MClepton ){
         h_Ev_den->Fill(MC_incoming_P[3]);
         h_Pmu_den->Fill(Pmu);
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
      std::cout<<"There are no reco tracks... bye"<<std::endl;
      return; 
    }
    std::cout<<"Found this many reco tracks "<<n_recoTrack<<std::endl;

    double Efrac_lepton =0.0;
    double Efrac_proton =0.0;
    double Efrac_pionplus =0.0;
    double Efrac_pionminus =0.0;
    double trackLength_lepton =0.0;
    double trackLength_proton =0.0;
    double trackLength_pion_plus =0.0;
    double trackLength_pion_minus =0.0;
    const simb::MCParticle *MClepton_reco = NULL; 
    const simb::MCParticle *MCproton_reco = NULL;
    const simb::MCParticle *MCpion_plus_reco = NULL;
    const simb::MCParticle *MCpion_minus_reco = NULL;
    for(int i=0; i<n_recoTrack; i++) {
       art::Ptr<recob::Track> track = tracklist[i];
       const TVector3 tmp_track_vtx = track->Vertex();
       double track_vtx[4] ={tmp_track_vtx[0], tmp_track_vtx[1], tmp_track_vtx[2], -999};
       bool track_isFiducial = insideFV( track_vtx );
       if( !track_isFiducial ) continue;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       truthMatcher( all_trackHits, particle, tmpEfrac );
       //std::cout<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac<<std::endl;
       if (!particle) continue;
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
 
    }

    Reco_LengthRes =  truth_lengthLepton-trackLength_lepton;
    Reco_LengthResProton = proton_length-trackLength_proton; 
    Reco_LengthResPionPlus = pion_plus_length-trackLength_pion_plus; 
    Reco_LengthResPionMinus = pion_minus_length-trackLength_pion_minus;
    if( MClepton_reco && MClepton  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){ 
        MC_LeptonTrack = 1;
        h_Pmu_num->Fill(Pmu);
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

   if(TEfficiency::CheckConsistency(*h_Ev_num,*h_Ev_den)){
      h_Eff_Ev = new TEfficiency(*h_Ev_num,*h_Ev_den);
      TGraphAsymmErrors *grEff_Ev = h_Eff_Ev->CreateGraph();
      grEff_Ev->Write("grEff_Ev");
   }
   if(TEfficiency::CheckConsistency(*h_Pmu_num,*h_Pmu_den)){ 
     h_Eff_Pmu = new TEfficiency(*h_Pmu_num,*h_Pmu_den);
     TGraphAsymmErrors *grEff_Pmu = h_Eff_Pmu->CreateGraph();
     grEff_Pmu->Write("grEff_Pmu");
   }
   if(TEfficiency::CheckConsistency(*h_theta_num,*h_theta_den)){
     h_Eff_theta = new TEfficiency(*h_theta_num,*h_theta_den);
     TGraphAsymmErrors *grEff_theta = h_Eff_theta->CreateGraph();
     grEff_theta->Write("grEff_theta");
   }
   if(TEfficiency::CheckConsistency(*h_Pproton_num,*h_Pproton_den)){
     h_Eff_Pproton = new TEfficiency(*h_Pproton_num,*h_Pproton_den);
     TGraphAsymmErrors *grEff_Pproton = h_Eff_Pproton->CreateGraph();
     grEff_Pproton->Write("grEff_Pproton");
   }
   if(TEfficiency::CheckConsistency(*h_Ppion_plus_num,*h_Ppion_plus_den)){
     h_Eff_Ppion_plus = new TEfficiency(*h_Ppion_plus_num,*h_Ppion_plus_den);
     TGraphAsymmErrors *grEff_Ppion_plus = h_Eff_Ppion_plus->CreateGraph();
     grEff_Ppion_plus->Write("grEff_Ppion_plus");
   }
   if(TEfficiency::CheckConsistency(*h_Ppion_minus_num,*h_Ppion_minus_den)){
     h_Eff_Ppion_minus = new TEfficiency(*h_Ppion_minus_num,*h_Ppion_minus_den);
     TGraphAsymmErrors *grEff_Ppion_minus = h_Eff_Ppion_minus->CreateGraph();
     grEff_Ppion_minus->Write("grEff_Ppion_minus");
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
 
   for(int i = 0; i<4; i++){
      MC_vertex[i] = 0;
      MC_incoming_P[i] =  0;
      MC_lepton_startMomentum[i] = 0;
      MC_lepton_endMomentum[i] = 0;
      MC_lepton_startXYZT[i] = 0;
      MC_lepton_endXYZT[i] = 0;
   }

   Reco_EfracLepton = -999.0; 
   Reco_LengthRes = -999.0;
   Reco_LengthResProton = -999.0;
   Reco_LengthResPionPlus = -999.0;
   Reco_LengthResPionMinus = -999.0;

  
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
