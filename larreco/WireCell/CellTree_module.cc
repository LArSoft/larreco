// Read data from MC raw files and convert it into ROOT tree
// Chao Zhang (chao@bnl.gov) 5/13/2014

#ifndef CELLTREE_MODULE
#define CELLTREE_MODULE

// LArSoft includes
#include "lardata/Utilities/GeometryUtilities.h"

#include "larsimobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
// #include "MCCheater/BackTracker.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"


// ROOT includes.
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TUnixSystem.h"
#include "TDatabasePDG.h"
#include "TObjArray.h"

// C++ Includes
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>

#define MAX_TRACKS 30000

using namespace std;

namespace wc {

class CellTree : public art::EDAnalyzer {
public:

    explicit CellTree(fhicl::ParameterSet const& pset);
    virtual ~CellTree();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void initOutput();
    void printEvent();
    void print_vector(ostream& out, vector<double>& v, TString desc, bool end=false);


    void processRaw(const art::Event& evt);
    void processCalib(const art::Event& evt);
    void processSpacePoint( const art::Event& event, TString option, ostream& out=cout);
    void processSimChannel(const art::Event& evt);
    void processMC(const art::Event& evt);
    void processMCTracks();

    void reset();
    void InitProcessMap();

    bool IsPrimary(int i) { return mc_mother[i] == 0 ; }
    bool KeepMC(int i);
    double KE(float* momentum);  // KE
    TString PDGName(int pdg);
    bool DumpMCJSON(int id, ostream& out);
    void DumpMCJSON(ostream& out=cout);


private:

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    std::string fCalibLabel;
    std::vector<std::string> fSpacePointLabels;
    std::string fOutFileName;
    std::string mcOption;
    bool fSaveMCTrackPoints;
    bool fSaveSimChannel;
    bool fSaveRaw;
    bool fSaveCalib;
    bool fSaveMC;
    bool fSaveJSON;
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

    // art::ServiceHandle<geo::Geometry> fGeom;
    // // auto const* larp = lar::providerFrom<detinfo::LArPropertiesService>();

    TFile *fOutFile;
    TTree *fEventTree;
    std::map<std::string, int> processMap;
    std::map<int, int> savedMCTrackIdMap;  // key: id; value: pdg;

    int entryNo;

    // Event Tree Leafs
    int fEvent;
    int fRun;
    int fSubRun;


    int fCalib_nChannel;
    // int fCalib_channelId[MAX_CHANNEL];  // hit channel id; size == fCalib_Nhit
    // // FIXEME:: cannot save e.g std::vector<std::vector<float> > in ttree
    std::vector<int> fCalib_channelId;
    // std::vector<std::vector<float> > fCalib_wf;
    TClonesArray *fCalib_wf;
    // std::vector<std::vector<int> > fCalib_wfTDC;


    int fRaw_nChannel;
    std::vector<int> fRaw_channelId;
    TClonesArray *fRaw_wf;

    int fSIMIDE_size;
    vector<int> fSIMIDE_channelIdY;
    vector<int> fSIMIDE_trackId;
    vector<unsigned short> fSIMIDE_tdc;
    vector<float> fSIMIDE_x;
    vector<float> fSIMIDE_y;
    vector<float> fSIMIDE_z;
    vector<float> fSIMIDE_numElectrons;

    int mc_Ntrack;  // number of tracks in MC
    int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
    int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
    int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
    int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
    float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
    float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
    float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
    float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > mc_daughters;  // daughters id of this track; vector
    TObjArray *fMC_trackPosition;

    int mc_isnu; // is neutrino interaction
    int mc_nGeniePrimaries; // number of Genie primaries
    int mc_nu_pdg; // pdg code of neutrino
    int mc_nu_ccnc; // cc or nc
    int mc_nu_mode; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
    int mc_nu_intType; // interaction type
    int mc_nu_target; // target interaction
    int mc_hitnuc; // hit nucleon
    int mc_hitquark; // hit quark
    double mc_nu_Q2; // Q^2
    double mc_nu_W; // W
    double mc_nu_X; // X
    double mc_nu_Y; // Y
    double mc_nu_Pt; // Pt
    double mc_nu_Theta; // angle relative to lepton
    float mc_nu_pos[4];  // interaction position of nu
    float mc_nu_mom[4];  // interaction momentum of nu

    // ----- derived ---
    std::map<int, int> trackIndex;
    std::vector<std::vector<int> > trackParents;
    std::vector<std::vector<int> > trackChildren;
    std::vector<std::vector<int> > trackSiblings;
    TDatabasePDG *dbPDG;

}; // class CellTree


//-----------------------------------------------------------------------
CellTree::CellTree(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    dbPDG = new TDatabasePDG();

    reconfigure(parameterSet);
    InitProcessMap();
    initOutput();
    entryNo = 0;
}

//-----------------------------------------------------------------------
CellTree::~CellTree()
{
}

//-----------------------------------------------------------------------
void CellTree::reconfigure(fhicl::ParameterSet const& p){
    fRawDigitLabel   = p.get<std::string>("RawDigitLabel");
    fCalibLabel      = p.get<std::string>("CalibLabel");
    fSpacePointLabels= p.get<std::vector<std::string> >("SpacePointLabels");
    fOutFileName     = p.get<std::string>("outFile");
    mcOption        = p.get<std::string>("mcOption");
    fSaveMCTrackPoints = p.get<bool>("saveMCTrackPoints");
    fSaveRaw         = p.get<bool>("saveRaw");
    fSaveCalib       = p.get<bool>("saveCalib");
    fSaveMC          = p.get<bool>("saveMC");
    fSaveSimChannel  = p.get<bool>("saveSimChannel");
    fSaveJSON        = p.get<bool>("saveJSON");
}

//-----------------------------------------------------------------------
void CellTree::initOutput()
{
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");

    // 3.1: add mc_trackPosition
    TNamed version("version", "3.1");
    version.Write();

    // init Event TTree
    TDirectory* subDir = fOutFile->mkdir("Event");
    subDir->cd();
    fEventTree = new TTree("Sim", "Event Tree from Simulation");
    fEventTree->Branch("eventNo", &fEvent);
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);

    fEventTree->Branch("raw_nChannel", &fRaw_nChannel);  // number of hit channels above threshold
    fEventTree->Branch("raw_channelId" , &fRaw_channelId); // hit channel id; size == raw_nChannel
    fRaw_wf = new TClonesArray("TH1F");
    fEventTree->Branch("raw_wf", &fRaw_wf, 256000, 0);  // raw waveform adc of each channel


    fEventTree->Branch("calib_nChannel", &fCalib_nChannel);  // number of hit channels above threshold
    fEventTree->Branch("calib_channelId" , &fCalib_channelId); // hit channel id; size == calib_Nhit
    fCalib_wf = new TClonesArray("TH1F");
    fEventTree->Branch("calib_wf", &fCalib_wf, 256000, 0);  // calib waveform adc of each channel
    // fCalib_wf->BypassStreamer();
    // fEventTree->Branch("calib_wfTDC", &fCalib_wfTDC);  // calib waveform tdc of each channel

    fEventTree->Branch("simide_size", &fSIMIDE_size);  // size of stored sim:IDE
    fEventTree->Branch("simide_channelIdY", &fSIMIDE_channelIdY);
    fEventTree->Branch("simide_trackId", &fSIMIDE_trackId);
    fEventTree->Branch("simide_tdc", &fSIMIDE_tdc);
    fEventTree->Branch("simide_x", &fSIMIDE_x);
    fEventTree->Branch("simide_y", &fSIMIDE_y);
    fEventTree->Branch("simide_z", &fSIMIDE_z);
    fEventTree->Branch("simide_numElectrons", &fSIMIDE_numElectrons);

    fEventTree->Branch("mc_Ntrack", &mc_Ntrack);  // number of tracks in MC
    fEventTree->Branch("mc_id", &mc_id, "mc_id[mc_Ntrack]/I");  // track id; size == mc_Ntrack
    fEventTree->Branch("mc_pdg", &mc_pdg, "mc_id[mc_Ntrack]/I");  // track particle pdg; size == mc_Ntrack
    fEventTree->Branch("mc_process", &mc_process, "mc_process[mc_Ntrack]/I");  // track generation process code; size == mc_Ntrack
    fEventTree->Branch("mc_mother", &mc_mother, "mc_mother[mc_Ntrack]/I");  // mother id of this track; size == mc_Ntrack
    fEventTree->Branch("mc_daughters", &mc_daughters);  // daughters id of this track; vector
    fEventTree->Branch("mc_startXYZT", &mc_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
    fEventTree->Branch("mc_endXYZT", &mc_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
    fEventTree->Branch("mc_startMomentum", &mc_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
    fEventTree->Branch("mc_endMomentum", &mc_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
    fMC_trackPosition = new TObjArray();
    fMC_trackPosition->SetOwner(kTRUE);
    fEventTree->Branch("mc_trackPosition", &fMC_trackPosition);

    fEventTree->Branch("mc_isnu", &mc_isnu);
    fEventTree->Branch("mc_nGeniePrimaries", &mc_nGeniePrimaries);
    fEventTree->Branch("mc_nu_pdg", &mc_nu_pdg);
    fEventTree->Branch("mc_nu_ccnc", &mc_nu_ccnc);
    fEventTree->Branch("mc_nu_mode", &mc_nu_mode);
    fEventTree->Branch("mc_nu_intType", &mc_nu_intType);
    fEventTree->Branch("mc_nu_target", &mc_nu_target);
    fEventTree->Branch("mc_hitnuc", &mc_hitnuc);
    fEventTree->Branch("mc_hitquark", &mc_hitquark);
    fEventTree->Branch("mc_nu_Q2", &mc_nu_Q2);
    fEventTree->Branch("mc_nu_W", &mc_nu_W);
    fEventTree->Branch("mc_nu_X", &mc_nu_X);
    fEventTree->Branch("mc_nu_Y", &mc_nu_Y);
    fEventTree->Branch("mc_nu_Pt", &mc_nu_Pt);
    fEventTree->Branch("mc_nu_Theta", &mc_nu_Theta);
    fEventTree->Branch("mc_nu_pos", &mc_nu_pos, "mc_nu_pos[4]/F");
    fEventTree->Branch("mc_nu_mom", &mc_nu_mom, "mc_nu_mom[4]/F");

    gDirectory = tmpDir;

    if (fSaveJSON) {
        system("rm -rf bee");
        gSystem->MakeDirectory("bee");
        gSystem->ChangeDirectory("bee");
        gSystem->MakeDirectory("data");
    }

}

//-----------------------------------------------------------------------
void CellTree::beginJob()
{


}


//-----------------------------------------------------------------------
void CellTree::endJob()
{
    // Write fEventTree to file
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Event");

    fEventTree->Write();

    gDirectory = tmpDir;

    fOutFile->Close();

    if (fSaveJSON) {
        system("zip -r bee_upload data");
        gSystem->ChangeDirectory("..");
    }
}

//-----------------------------------------------------------------------
void CellTree::beginRun(const art::Run& /*run*/)
{
  mf::LogInfo("CellTree") << "begin run";
}

//-----------------------------------------------------------------------
void CellTree::analyze( const art::Event& event )
{
    reset();
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    if (fSaveRaw) processRaw(event);
    if (fSaveCalib) processCalib(event);
    if (fSaveSimChannel) processSimChannel(event);
    if (fSaveMC) processMC(event);

    if (fSaveJSON) {
        gSystem->MakeDirectory(TString::Format("data/%i", entryNo).Data());
        int nSp = fSpacePointLabels.size();
        for (int i=0; i<nSp; i++) {
            TString jsonfile;
            jsonfile.Form("data/%i/%i-%s.json", entryNo, entryNo, fSpacePointLabels[i].c_str());
            std::ofstream out(jsonfile.Data());
            processSpacePoint(event, fSpacePointLabels[i], out);
            out.close();
        }

        if(fSaveMC) {
            processMCTracks();
            TString jsonfile;
            jsonfile.Form("data/%i/%i-mc.json", entryNo, entryNo);
            std::ofstream out(jsonfile.Data());
            DumpMCJSON(out);
            out.close();
        }
    }

    // printEvent();
    fEventTree->Fill();

    entryNo++;
}

//-----------------------------------------------------------------------
void CellTree::reset()
{

    fRaw_channelId.clear();
    // fRaw_wf->Clear();
    fRaw_wf->Delete();

    fCalib_channelId.clear();
    fCalib_wf->Clear();

    fSIMIDE_channelIdY.clear();
    fSIMIDE_trackId.clear();
    fSIMIDE_tdc.clear();
    fSIMIDE_x.clear();
    fSIMIDE_y.clear();
    fSIMIDE_z.clear();
    fSIMIDE_numElectrons.clear();

    for (int i=0; i<MAX_TRACKS; i++) {
        mc_id[i] = 0;
        mc_pdg[i] = 0;
        mc_mother[i] = 0;
        for (int j=0; j<4; j++) {
            mc_startXYZT[i][j]      = 0;
            mc_endXYZT[i][j]        = 0;
            mc_startMomentum[i][j] = 0;
            mc_endMomentum[i][j]   = 0;
        }
    }
    mc_daughters.clear();
    savedMCTrackIdMap.clear();
    fMC_trackPosition->Clear();

    mc_isnu = 0;
    mc_nGeniePrimaries = -1;
    mc_nu_pdg = -1;
    mc_nu_ccnc = -1;
    mc_nu_mode = -1;
    mc_nu_intType = -1;
    mc_nu_target = -1;
    mc_hitnuc = -1;
    mc_hitquark = -1;
    mc_nu_Q2 = -1;
    mc_nu_W = -1;
    mc_nu_X = -1;
    mc_nu_Y = -1;
    mc_nu_Pt = -1;
    mc_nu_Theta = -1;
    for (int i=0; i<4; i++) {
        mc_nu_pos[i] = 0;
        mc_nu_mom[i] = 0;
    }

    trackIndex.clear();
    trackParents.clear();
    trackChildren.clear();
    trackSiblings.clear();
}

//-----------------------------------------------------------------------
void CellTree::processRaw( const art::Event& event )
{
    art::Handle< std::vector<raw::RawDigit> > rawdigit;
    if (! event.getByLabel(fRawDigitLabel, rawdigit)) {
        cout << "WARNING: no label " << fRawDigitLabel << endl;
        return;
    }
    std::vector< art::Ptr<raw::RawDigit> >  wires;
    art::fill_ptr_vector(wires, rawdigit);

    fRaw_nChannel = wires.size();

    int i=0;
    for (auto const& wire: wires) {
        int chanId = wire->Channel();
        fRaw_channelId.push_back(chanId);

        int nSamples = wire->Samples();
        std::vector<short> uncompressed(nSamples);
        raw::Uncompress(wire->ADCs(), uncompressed, wire->Compression());

        TH1F *h = new((*fRaw_wf)[i]) TH1F("", "", 9600, 0, 9600);
        for (int j=1; j<=nSamples; j++) {
            h->SetBinContent(j, uncompressed[j-1]);
        }
        i++;
    }

}

//-----------------------------------------------------------------------
void CellTree::processCalib( const art::Event& event )
{

    art::Handle< std::vector<recob::Wire> > wires_handle;
    if (! event.getByLabel(fCalibLabel, wires_handle)) {
        cout << "WARNING: no label " << fCalibLabel << endl;
        return;
    }
    std::vector< art::Ptr<recob::Wire> >  wires;
    art::fill_ptr_vector(wires, wires_handle);

    // wires size should == Nchannels == 1992; (no hit channel has a flat 0-waveform)
    // cout << "\n wires size: " << wires.size() << endl;
    fCalib_nChannel = wires.size();

    int i=0;
    for (auto const& wire: wires) {
        std::vector<float> calibwf = wire->Signal();
        int chanId = wire->Channel();
        fCalib_channelId.push_back(chanId);
        TH1F *h = new((*fCalib_wf)[i]) TH1F("", "", 9600, 0, 9600);
        for (int j=1; j<=9600; j++) {
            h->SetBinContent(j, calibwf[j]);
        }
        // fCalib_wf.push_back(calibwf);
        // cout << chanId << ", " << nSamples << endl;
        i++;
    }

}

//-----------------------------------------------------------------------
void CellTree::processSimChannel( const art::Event& event )
{
    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel("largeant", simChannelHandle);
    // cout << "total simChannel: " << (*simChannelHandle).size() << endl;
    fSIMIDE_size = 0;
    for ( auto const& channel : (*simChannelHandle) ) {
        auto channelNumber = channel.Channel();
        // cout << channelNumber << endl;
        // if (! (fGeometry->SignalType( channelNumber ) == geo::kCollection) ) {
        //     continue;
        // }
        auto const& timeSlices = channel.TDCIDEMap();
        for ( auto const& timeSlice : timeSlices ) {
            auto const& energyDeposits = timeSlice.second;
            for ( auto const& energyDeposit : energyDeposits ) {
                fSIMIDE_size++;
                fSIMIDE_channelIdY.push_back(channelNumber);
                fSIMIDE_tdc.push_back(timeSlice.first);
                fSIMIDE_trackId.push_back(energyDeposit.trackID);
                fSIMIDE_x.push_back(energyDeposit.x);
                fSIMIDE_y.push_back(energyDeposit.y);
                fSIMIDE_z.push_back(energyDeposit.z);
                fSIMIDE_numElectrons.push_back(energyDeposit.numElectrons);
                // cout << channelNumber << ": " << energyDeposit.trackID << ": " << timeSlice.first << ": "
                //      << energyDeposit.x << ", " << energyDeposit.y << ", " << energyDeposit.z << ", "
                //      << energyDeposit.numElectrons << endl;
            }
        }

    }
    cout << "total IDEs: " << fSIMIDE_size << endl;


}

//-----------------------------------------------------------------------
void CellTree::processMC( const art::Event& event )
{
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    if (! event.getByLabel("largeant", particleHandle)) return;
    std::vector< art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);

    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel("largeant", simChannelHandle);

    // art::ServiceHandle<cheat::BackTracker> bt;
    art::FindOneP<simb::MCTruth> fo(particleHandle, event, "largeant");

    int i=0; // track index in saved MCParticles;
    int i_all=0; // track index in all MCParticles;
    for (auto const& particle: particles ) {
        art::Ptr<simb::MCTruth> mctruth = fo.at(i_all);
        i_all++;

        if (mcOption == "nuOnly") {
            if ( !(mctruth->Origin() == 1 && particle->Mother() == 0) ) {
                continue;
            }
        }

        // if (mctruth->Origin() == 1 || mc_mother[i] == 0) {
        //     cout << "process: " << particle->Process()
        //          << ", id: " << mc_id[i]
        //          << ", pdg: " <<  mc_pdg[i]
        //          << ", mother: " << mc_mother[i]
        //          << ", nDaughter: " << (particle->NumberDaughters())
        //          << ", truth: " << mctruth->Origin()
        //          << endl;
        // }
        // const art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(mc_id[i]);

        mc_process[i] = processMap[particle->Process()];
        if (mc_process[i] == 0) cout << "unknown process: " << particle->Process() << endl;
        mc_id[i] = particle->TrackId();
        mc_pdg[i] = particle->PdgCode();
        mc_mother[i] = particle->Mother();
        savedMCTrackIdMap[mc_id[i]] = mc_pdg[i];

        int Ndaughters = particle->NumberDaughters();
        vector<int> daughters;
        for (int i=0; i<Ndaughters; i++) {
            daughters.push_back(particle->Daughter(i));
        }
        mc_daughters.push_back(daughters);
        size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
        int last = numberTrajectoryPoints - 1;
        const TLorentzVector& positionStart = particle->Position(0);
        const TLorentzVector& positionEnd   = particle->Position(last);
        const TLorentzVector& momentumStart = particle->Momentum(0);
        const TLorentzVector& momentumEnd   = particle->Momentum(last);
        positionStart.GetXYZT(mc_startXYZT[i]);
        positionEnd.GetXYZT(mc_endXYZT[i]);
        momentumStart.GetXYZT(mc_startMomentum[i]);
        momentumEnd.GetXYZT(mc_endMomentum[i]);

        if (fSaveMCTrackPoints) {
            TClonesArray *Lposition = new TClonesArray("TLorentzVector", numberTrajectoryPoints);
            // Read the position and momentum along this particle track
            for(unsigned int j=0; j<numberTrajectoryPoints; j++) {
                new ((*Lposition)[j]) TLorentzVector(particle->Position(j));
            }
            fMC_trackPosition->Add(Lposition);
        }

        i++;
        if (i==MAX_TRACKS) {
            cout << "WARNING:: # tracks exceeded MAX_TRACKS " << MAX_TRACKS << endl;
            break;
        }
    } // particle loop done
    mc_Ntrack = i;
    // cout << "MC_Ntracks:" << mc_Ntrack << endl;

    // Generator Info
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    event.getByLabel("generator",mctruthListHandle);
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mctruthListHandle);
    art::Ptr<simb::MCTruth> mctruth;

    if (mclist.size()>0) {
        mctruth = mclist.at(0);
        if (mctruth->NeutrinoSet()) {
            simb::MCNeutrino nu = mctruth->GetNeutrino();
            mc_isnu = 1;
            mc_nGeniePrimaries = mctruth->NParticles();
            mc_nu_pdg = nu.Nu().PdgCode();
            mc_nu_ccnc = nu.CCNC();
            mc_nu_mode = nu.Mode();
            mc_nu_intType = nu.InteractionType();
            mc_nu_target = nu.Target();
            mc_hitnuc = nu.HitNuc();
            mc_hitquark = nu.HitQuark();
            mc_nu_Q2 = nu.QSqr();
            mc_nu_W = nu.W();
            mc_nu_X = nu.X();
            mc_nu_Y = nu.Y();
            mc_nu_Pt = nu.Pt();
            mc_nu_Theta = nu.Theta();

            const TLorentzVector& position = nu.Nu().Position(0);
            const TLorentzVector& momentum = nu.Nu().Momentum(0);
            position.GetXYZT(mc_nu_pos);
            momentum.GetXYZT(mc_nu_mom);
            // cout << "nu: " << mc_nu_pdg << ", nPrim: " << mc_nGeniePrimaries
            //      << ", ccnc: " << mc_nu_ccnc << endl;
            // for (int i=0; i<mc_nGeniePrimaries; i++) {
            //     simb::MCParticle particle = mctruth->GetParticle(i);
            //     cout << "id: " << particle.TrackId()
            //          << ", pdg: " << particle.PdgCode()
            //          << endl;
            // }
        }
    }

}


//-----------------------------------------------------------------------
void CellTree::processSpacePoint( const art::Event& event, TString option, ostream& out)
{

    art::Handle< std::vector<recob::SpacePoint> > handle;
    if (! event.getByLabel(option.Data(), handle)) {
        cout << "WARNING: no label " << option << endl;
        return;
    }
    std::vector< art::Ptr<recob::SpacePoint> >  sps;
    art::fill_ptr_vector(sps, handle);

    double x=0, y=0, z=0, q=0, nq=1;
    vector<double> vx, vy, vz, vq, vnq;

    for (auto const& sp: sps ) {
        // cout << sp->XYZ()[0] << ", " << sp->XYZ()[1] << ", " << sp->XYZ()[2] << endl;
        x = sp->XYZ()[0];
        y = sp->XYZ()[1];
        z = sp->XYZ()[2];
        vx.push_back(x);
        vy.push_back(y);
        vz.push_back(z);
        vq.push_back(q);
        vnq.push_back(nq);
    }

    out << fixed << setprecision(1);
    out << "{" << endl;

    out << '"' << "runNo" << '"' << ":" << '"' << fRun << '"' << "," << endl;
    out << '"' << "subRunNo" << '"' << ":" << '"' << fSubRun << '"' << "," << endl;
    out << '"' << "eventNo" << '"' << ":" << '"' << fEvent << '"' << "," << endl;

    TString geomName(fGeometry->DetectorName().c_str());
    if (geomName.Contains("35t")) { geomName = "dune35t"; }
    else if (geomName.Contains("protodune")) { geomName = "protodune"; }
    else if (geomName.Contains("workspace")) { geomName = "dune10kt_workspace"; }
    else { geomName = "uboone"; } // use uboone as default
    out << '"' << "geom" << '"' << ":" << '"' << geomName << '"' << "," << endl;


    print_vector(out, vx, "x");
    print_vector(out, vy, "y");
    print_vector(out, vz, "z");

    out << fixed << setprecision(0);
    print_vector(out, vq, "q");
    print_vector(out, vnq, "nq");

    out << '"' << "type" << '"' << ":" << '"' << option << '"' << endl;
    out << "}" << endl;
}

//-----------------------------------------------------------------------
void CellTree::print_vector(ostream& out, vector<double>& v, TString desc, bool end)
{
    int N = v.size();

    out << '"' << desc << '"' << ":[";
    for (int i=0; i<N; i++) {
        out << v[i];
        if (i!=N-1) {
            out << ",";
        }
    }
    out << "]";
    if (!end) out << ",";
    out << endl;
}

//-----------------------------------------------------------------------
void CellTree::processMCTracks()
{
    // map track id to track index in the array
    for (int i=0; i<mc_Ntrack; i++) {
        trackIndex[mc_id[i]] = i;
    }

    // in trackParents, trackChildren, trackSiblings vectors, store track index (not track id)
    for (int i=0; i<mc_Ntrack; i++) {
        // currently, parent size == 1;
        // for primary particle, parent id = 0;
        vector<int> parents;
        if ( !IsPrimary(i) ) {
            parents.push_back(trackIndex[mc_mother[i]]);
        }
        trackParents.push_back(parents); // primary track will have 0 parents

        vector<int> children;
        int nChildren = mc_daughters.at(i).size();
        for (int j=0; j<nChildren; j++) {
            children.push_back(trackIndex[mc_daughters.at(i).at(j)]);
        }
        trackChildren.push_back(children);

    }

    // siblings
    for (int i=0; i<mc_Ntrack; i++) {
        vector<int> siblings;
        if ( IsPrimary(i) ) {
            for (int j=0; j<mc_Ntrack; j++) {
                if( IsPrimary(j) ) {
                    siblings.push_back(j);
                }
            }
        }
        else {
            // siblings are simply children of the mother
            int mother = trackIndex[mc_mother[i]];
            int nSiblings = trackChildren.at(mother).size();
            for (int j=0; j<nSiblings; j++) {
                siblings.push_back(trackChildren.at(mother).at(j));
            }
        }
        trackSiblings.push_back(siblings);
    }

}

//-----------------------------------------------------------------------
bool CellTree::DumpMCJSON(int id, ostream& out)
{
    int i = trackIndex[id];
    if (!KeepMC(i)) return false;

    int e = KE(mc_startMomentum[i])*1000;

    int nDaughter = mc_daughters.at(i).size();
    vector<int> saved_daughters;
    for (int j=0; j<nDaughter; j++) {
        int daughter_id = mc_daughters.at(i).at(j);
        // int e_daughter = KE(mc_startMomentum[ trackIndex[daughter_id] ])*1000;
        // if (e_daughter >= thresh_KE) {
        if ( KeepMC(trackIndex[daughter_id]) ) {
            saved_daughters.push_back(daughter_id);
        }
    }

    out << fixed << setprecision(1);
    out << "{";

    out << "\"id\":" << id << ",";
    out << "\"text\":" << "\"" << PDGName(mc_pdg[i]) << "  " << e << " MeV\",";
    out << "\"data\":{";
    out << "\"start\":[" << mc_startXYZT[i][0] << ", " <<  mc_startXYZT[i][1] << ", " << mc_startXYZT[i][2] << "],";
    out << "\"end\":[" << mc_endXYZT[i][0] << ", " <<  mc_endXYZT[i][1] << ", " << mc_endXYZT[i][2] << "]";
    out << "},";
    out << "\"children\":[";
    int nSavedDaughter = saved_daughters.size();
    if (nSavedDaughter == 0) {
        out << "],";
        out << "\"icon\":" << "\"jstree-file\"";
        out << "}";
        return true;
    }
    else {
        for (int j=0; j<nSavedDaughter; j++) {
            DumpMCJSON(saved_daughters.at(j), out);
            if (j!=nSavedDaughter-1) {
                out << ",";
            }
        }
        out << "]";
        out << "}";
        return true;
    }
}

//-----------------------------------------------------------------------
void CellTree::DumpMCJSON(ostream& out)
{
    out << "[";
    vector<int> primaries;
    for (int i=0; i<mc_Ntrack; i++) {
        if (IsPrimary(i)) {
            // int e = KE(mc_startMomentum[i])*1000;
            // if (e<thresh_KE) continue;
            if (KeepMC(i)) {
                primaries.push_back(i);
            }
        }
    }
    int size = primaries.size();
    // cout << size << endl;
    for (int i=0; i<size; i++) {
        if (DumpMCJSON(mc_id[primaries[i]], out) && i!=size-1) {
            out << ", ";
        }
    }

    out << "]";
}

//-----------------------------------------------------------------------
double CellTree::KE(float* momentum)
{
    TLorentzVector particle(momentum);
    return particle.E()-particle.M();
}

//-----------------------------------------------------------------------
bool CellTree::KeepMC(int i)
{
    double e = KE(mc_startMomentum[i])*1000;
    double thresh_KE_em = 5.; // MeV
    double thresh_KE_np = 10; // MeV
    if (mc_pdg[i]==22 || mc_pdg[i]==11 || mc_pdg[i]==-11) {
        if (e>=thresh_KE_em) return true;
        else return false;
    }
    else if (mc_pdg[i]==2112 || mc_pdg[i]==2212 || mc_pdg[i]>1e9) {
        if (e>=thresh_KE_np) return true;
        else return false;
    }
    return true;
}

//-----------------------------------------------------------------------
TString CellTree::PDGName(int pdg)
{
    TParticlePDG *p = dbPDG->GetParticle(pdg);
    if (p == 0) {
        if (pdg>1e9) {
            int z = (pdg - 1e9) / 10000;
            int a = (pdg - 1e9 - z*1e4) / 10;
            TString name;
            if (z == 18) name = "Ar";

            else if (z == 17) name = "Cl";
            else if (z == 19) name = "Ca";
            else if (z == 16) name = "S";
            else if (z == 15) name = "P";
            else if (z == 14) name = "Si";
            else if (z == 1) name = "H";
            else if (z == 2) name = "He";

            else return Form("%i", pdg);
            return Form("%s-%i", name.Data(), a);
        }
        return Form("%i", pdg);
    }
    else {
        return p->GetName();
    }
}

//-----------------------------------------------------------------------
void CellTree::printEvent()
{
    cout << " Run/SubRun/Event: " << fRun << "/" << fSubRun << "/" << fEvent << endl;
    cout << "      Ntracks:" << mc_Ntrack << endl;

    for (int i=0; i<mc_Ntrack; i++) {
        cout << "\n              id: " << mc_id[i];
        cout << "\n             pdg: " << mc_pdg[i];
        cout << "\n          mother: " << mc_mother[i];
        cout << "\n      Ndaughters: " << mc_daughters.at(i).size();
        cout << "\n      start XYZT: (" << mc_startXYZT[i][0] << ", " << mc_startXYZT[i][1] << ", " << mc_startXYZT[i][2] << ", " << mc_startXYZT[i][3] << ")";
        cout << "\n        end XYZT: (" << mc_endXYZT[i][0] << ", " << mc_endXYZT[i][1] << ", " << mc_endXYZT[i][2] << ", " << mc_endXYZT[i][3] << ")";
        cout << "\n  start momentum: (" << mc_startMomentum[i][0] << ", " << mc_startMomentum[i][1] << ", " << mc_startMomentum[i][2] << ", " << mc_startMomentum[i][3] << ")";
        cout << "\n    end momentum: (" << mc_endMomentum[i][0] << ", " << mc_endMomentum[i][1] << ", " << mc_endMomentum[i][2] << ", " << mc_endMomentum[i][3] << ")";

        cout << endl;
    }
}

//-----------------------------------------------------------------------
void CellTree::InitProcessMap()
{
    processMap["unknown"]              = 0;
    processMap["primary"]              = 1;
    processMap["compt"]                = 2;
    processMap["phot"]                 = 3;
    processMap["annihil"]              = 4;
    processMap["eIoni"]                = 5;
    processMap["eBrem"]                = 6;
    processMap["conv"]                 = 7;
    processMap["muIoni"]               = 8;
    processMap["muMinusCaptureAtRest"] = 9;
    processMap["NeutronInelastic"]     = 10;
    processMap["nCapture"]             = 11;
    processMap["hadElastic"]           = 12;
    processMap["Decay"]                = 13;
    processMap["CoulombScat"]          = 14;
    processMap["muPairProd"]           = 15;
    processMap["muBrems"]              = 16;
    processMap["muPairProd"]           = 17;
    processMap["PhotonInelastic"]      = 18;
    processMap["hIoni"]                = 19;
    processMap["ProtonInelastic"]      = 20;
    processMap["PionPlusInelastic"]    = 21;
    processMap["CHIPSNuclearCaptureAtRest"] = 22;
    processMap["PionMinusInelastic"]   = 23;
}

//-----------------------------------------------------------------------
DEFINE_ART_MODULE(CellTree)
} // namespace microboone


#endif
