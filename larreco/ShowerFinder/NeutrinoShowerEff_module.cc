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
    void initOutput();

    void processEff(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac);
    bool insideFV(double vertex[4]);
    void doEfficiencies();
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fOutFileName;
    std::string fMCTruthModuleLabel;
    std::string fShowerModuleLabel;
    int         fNeutrinoPDGcode;
    int		fLeptonPDGcode;
    double      fMaxNeutrinoE;
    bool	fSaveMCTree; 
    double 	fMaxEfrac;

    TFile *fOutFile;
    TTree *fEventTree;
    TTree *fHitsTree;

    TH1D *h_Ev_den;
    TH1D *h_Ev_num;
    TH1D *h_Pe_den;
    TH1D *h_Pe_num;
    TH1D *h_theta_den;
    TH1D *h_theta_num;

    TH1D *h_Efrac_shContamination;
    TH1D *h_Efrac_lepton;     

    TEfficiency* h_Eff_Ev = 0;
    TEfficiency* h_Eff_Pe = 0;
    TEfficiency* h_Eff_theta = 0;

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

 }; // class NeutrinoShowerEff


//========================================================================
NeutrinoShowerEff::NeutrinoShowerEff(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
    initOutput();
}
//========================================================================
NeutrinoShowerEff::~NeutrinoShowerEff(){
  //destructor
}
//========================================================================
void NeutrinoShowerEff::reconfigure(fhicl::ParameterSet const& p){

    fOutFileName         = p.get<std::string>("outFile");
    fMCTruthModuleLabel  = p.get<std::string>("MCTruthModuleLabel");
    fShowerModuleLabel   = p.get<std::string>("ShowerModuleLabel");
    fLeptonPDGcode       = p.get<int>("LeptonPDGcode");
    fNeutrinoPDGcode     = p.get<int>("NeutrinoPDGcode");
    fMaxNeutrinoE	 = p.get<double>("MaxNeutrinoE");
    fMaxEfrac		 = p.get<double>("MaxEfrac");
    fSaveMCTree		 = p.get<bool>("SaveMCTree");
}
//========================================================================
void NeutrinoShowerEff::initOutput(){
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");
    double E_bins[21] ={0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,14.0,17.0,20.0,25.0};
    double theta_bin[44]= { 0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,55.,60.,65.,70.,75.,80.,85.,90.};

    TDirectory* subDir = fOutFile->mkdir("Histograms");
    subDir->cd();
    h_Ev_den = new TH1D("h_Ev_den","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
    h_Ev_num = new TH1D("h_Ev_num","Neutrino Energy; Neutrino Energy (GeV); Tracking Efficiency",20,E_bins);
    h_Pe_den = new TH1D("h_Pe_den","Electron Momentum; Electron Momentum (GeV); Tracking Efficiency",20,E_bins);
    h_Pe_num = new TH1D("h_Pe_num","Electron Momentum; Electron Momentum (GeV); Tracking Efficiency",20,E_bins);
    h_theta_den = new TH1D("h_theta_den","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);
    h_theta_num = new TH1D("h_theta_num","Theta; Theta w.r.t beam direction (Degrees); Tracking Efficiency",43,theta_bin);

    h_Efrac_shContamination = new TH1D("h_Efrac_shContamination","Efrac Lepton; Energy fraction (contamination);",60,0,1.2);
 
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

      gDirectory = tmpDir;
    }

}
//========================================================================
void NeutrinoShowerEff::beginJob(){
  cout<<"job begin..."<<endl;
}
//========================================================================
void NeutrinoShowerEff::endJob(){

    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Histograms");
    h_Efrac_shContamination->Write();
    h_Ev_den->Write();
    h_Ev_num->Write();
    h_theta_den->Write();
    h_theta_num->Write();
     
    doEfficiencies();

    if( fSaveMCTree ){
      fOutFile->cd("/Events");
      fEventTree->Write();
      gDirectory = tmpDir;
    }
    fOutFile->Close();
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
         cout<<"Incoming E "<<MC_incoming_P[3]<<" is CC? "<<MC_isCC<<" nuPDG "<<MC_incoming_PDG<<" target "<<MC_target<<" vtx "<<MC_vertex[0]<<" "<<MC_vertex[1]<<" "<<MC_vertex[2]<<" "<<MC_vertex[3]<<endl;
       }
    }

    //!save lepton 
    simb::MCParticle *MClepton = NULL; 
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;

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
    } 
    //Saving denominator histograms for lepton pions and protons 
    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
    double Pe = sqrt(pow(MC_lepton_startMomentum[0],2)+pow(MC_lepton_startMomentum[1],2)+pow(MC_lepton_startMomentum[2],2));
    double Pv  = sqrt(pow(MC_incoming_P[0],2)+pow(MC_incoming_P[1],2)+pow(MC_incoming_P[2],2));
    double theta_e = acos((MC_incoming_P[0]*MC_lepton_startMomentum[0] + MC_incoming_P[1]*MC_lepton_startMomentum[1] +MC_incoming_P[2]*MC_lepton_startMomentum[2])/(Pv*Pe) );
    theta_e *= (180.0/3.14159);

    //save CC events within the fiducial volume with the favorite neutrino flavor 
    if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){
       if( MClepton ){
         h_Ev_den->Fill(MC_incoming_P[3]);
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
    if(!event.getByLabel(fShowerModuleLabel,showerHandle)) return;
    std::vector<art::Ptr<recob::Shower>> showerlist;
    art::fill_ptr_vector(showerlist, showerHandle);

    n_recoShowers= showerlist.size();
    if ( n_recoShowers == 0 || n_recoShowers> MAX_SHOWERS ) return;
    art::FindManyP<recob::Hit> sh_hitsAll(showerHandle, event, fShowerModuleLabel);
    cout<<"Found this many showers "<<n_recoShowers<<endl; 
    double Efrac_contamination= 999.0;
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
       const simb::MCParticle *particle;
       double tmpEfrac_contamination = 0.0;  //fraction of non EM energy contatiminatio (see truthMatcher for definition)
       int tmp_nHits = sh_hits.size();
       truthMatcher( sh_hits, particle, tmpEfrac_contamination );
       sh_Efrac_contamination[i] = tmpEfrac_contamination;
       sh_nHits[i] = tmp_nHits; 
       sh_hasPrimary_e[i] = 0;
       if( particle->PdgCode()  == fLeptonPDGcode && particle->TrackId() == MC_leptonID ) sh_hasPrimary_e[i] = 1;
       //cout<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac_contamination<<" "<<sh_hits.size()<<" "<<particle->TrackId()<<" "<<MC_leptonID<<endl;
       //save the best shower based on non EM and number of hits
      
       if( particle->PdgCode()  == fLeptonPDGcode && particle->TrackId() == MC_leptonID ){
         if( tmp_nHits > nHits ){
            nHits = tmp_nHits;
            Efrac_contamination = tmpEfrac_contamination;
            MClepton_reco = particle;
            sh_Efrac_best =Efrac_contamination; 
            //cout<<"this is the best shower "<<particle->PdgCode()<<" "<<particle->TrackId()<<" Efrac "<<tmpEfrac_contamination<<" "<<sh_hits.size()<<endl;
         } 
       }          
    }
   
    if( MClepton_reco && MClepton  ){
      if( MC_isCC && (fNeutrinoPDGcode == MC_incoming_PDG) && (MC_incoming_P[3] <= fMaxNeutrinoE) ){ 
        h_Efrac_shContamination->Fill(Efrac_contamination);
        if( Efrac_contamination < fMaxEfrac ){
          h_Pe_num->Fill(Pe);
          h_Ev_num->Fill(MC_incoming_P[3]);
          h_theta_num->Fill(theta_e);
        }
      }
    }
}
//========================================================================
void NeutrinoShowerEff::truthMatcher( std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac){

    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < shower_hits.size(); ++j){
       art::Ptr<recob::Hit> hit = shower_hits[j];
       //For know let's use collection plane to look at the shower reconstruction
       //if( hit->View() != 2) continue;
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t k = 0; k < TrackIDs.size(); k++){
          trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
       }            
    }
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double noEM_E = 0.0;  //non electromagnetic energy is defined as energy from charged pion and protons 
    if( !trkID_E.size() ) return; //Ghost shower???
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         max_E = ii->second;
         TrackID = ii->first;
       }
       int ID = ii->first;
       const simb::MCParticle *particle = bt->TrackIDToParticle(ID);
       //if( abs(particle->PdgCode()) == 211 || particle->PdgCode() == 2212 ){
       if( particle->PdgCode() != 22 && abs(particle->PdgCode()) != 11){
         noEM_E += ii->second;
       }
    } 

    MCparticle = bt->TrackIDToParticle(TrackID);
    Efrac = noEM_E/total_E;
}
//========================================================================
bool NeutrinoShowerEff::insideFV( double vertex[4]){ 

     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     //This is temporarily we should define a common FV    
     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     if( fabs(x) > 350.0 ) return false;
     else if( fabs(y) > 550.0 ) return false;
     else if( z< 0 || z> 400.0 ) return false;
     else return true;

}
//========================================================================
void NeutrinoShowerEff::doEfficiencies(){

   if(TEfficiency::CheckConsistency(*h_Ev_num,*h_Ev_den)){
      h_Eff_Ev = new TEfficiency(*h_Ev_num,*h_Ev_den);
      TGraphAsymmErrors *grEff_Ev = h_Eff_Ev->CreateGraph();
      grEff_Ev->Write("grEff_Ev");
   }
   if(TEfficiency::CheckConsistency(*h_Pe_num,*h_Pe_den)){ 
     h_Eff_Pe = new TEfficiency(*h_Pe_num,*h_Pe_den);
     TGraphAsymmErrors *grEff_Pe = h_Eff_Pe->CreateGraph();
     grEff_Pe->Write("grEff_Pe");
   }
   if(TEfficiency::CheckConsistency(*h_theta_num,*h_theta_den)){
     h_Eff_theta = new TEfficiency(*h_theta_num,*h_theta_den);
     TGraphAsymmErrors *grEff_theta = h_Eff_theta->CreateGraph();
     grEff_theta->Write("grEff_theta");
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
