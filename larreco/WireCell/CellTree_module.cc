// Read data from MC raw files and convert it into ROOT tree
// Chao Zhang (chao@bnl.gov) 5/13/2014
// Added optical info --- Brooke Russell (brussell@yale.edu) 1/31/2017

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes.
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TTree.h"

// C++ Includes
#include <cstdio>
#include <fstream>
#include <map>

#define MAX_TRACKS 30000

using namespace std;

namespace wc {

  class CellTree : public art::EDAnalyzer {
  public:
    explicit CellTree(fhicl::ParameterSet const& pset);

  private:
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void initOutput();
    void printEvent();
    void print_vector(ostream& out, vector<double>& v, TString desc, bool end = false);

    void processRaw(const art::Event& evt);
    void processCalib(const art::Event& evt);
    void processOpHit(const art::Event& evt);
    void processOpFlash(const art::Event& evt);
    void processSpacePoint(const art::Event& event, TString option, ostream& out = cout);
    void processSpacePointTruthDepo(const art::Event& event,
                                    TString option,
                                    ostream& out = cout,
                                    bool t0_corrected = true);
    void processSimChannel(const art::Event& evt);
    void processMC(const art::Event& evt);
    void processMCTracks();
    void processTrigger(const art::Event& evt);

    void reset();
    void InitProcessMap();

    bool IsPrimary(int i) { return mc_mother[i] == 0; }
    bool KeepMC(int i);
    double KE(float* momentum); // KE
    TString PDGName(int pdg);
    bool DumpMCJSON(int id, ostream& out);
    void DumpMCJSON(ostream& out = cout);

    // the parameters we'll read from the .fcl
    std::string fRawDigitLabel;
    std::string fCalibLabel;
    std::string fOpHitLabel;
    std::string fOpFlashLabel;
    std::string fTriggerLabel;
    std::string fSimEnergyDepositLabel;
    std::vector<std::string> fSpacePointLabels;
    std::string fSimChannelLabel;
    std::string fOutFileName;
    std::string mcOption;
    int nRawSamples;
    float opMultPEThresh;
    float drift_speed;
    bool fSaveMCTrackPoints;
    bool fSaveSimChannel;
    bool fSaveRaw;
    bool fSaveCalib;
    bool fSaveOpHit;
    bool fSaveOpFlash;
    bool fSaveMC;
    bool fSaveTrigger;
    bool fSaveJSON;
    bool fT0_corrected;
    art::ServiceHandle<geo::Geometry const> fGeometry; // pointer to Geometry service

    // art::ServiceHandle<geo::Geometry const> fGeom;
    // // auto const* larp = lar::providerFrom<detinfo::LArPropertiesService>();

    TFile* fOutFile;
    TTree* fEventTree;
    std::map<std::string, int> processMap;
    std::map<int, int> savedMCTrackIdMap; // key: id; value: pdg;

    int entryNo;

    // Event Tree Leafs
    int fEvent;
    int fRun;
    int fSubRun;
    double fEventTime;

    unsigned int fTriggernumber; //trigger counter
    double fTriggertime;         //trigger time w.r.t. electronics clock T0
    double fBeamgatetime;        //beamgate time w.r.t. electronics clock T0
    unsigned int fTriggerbits;   //trigger bits

    int fCalib_nChannel;
    // int fCalib_channelId[MAX_CHANNEL];  // hit channel id; size == fCalib_Nhit
    // // FIXEME:: cannot save e.g std::vector<std::vector<float> > in ttree
    std::vector<int> fCalib_channelId;
    // std::vector<std::vector<float> > fCalib_wf;
    TClonesArray* fCalib_wf;
    // std::vector<std::vector<int> > fCalib_wfTDC;

    int oh_nHits;
    vector<int> oh_channel;
    vector<double> oh_bgtime;
    vector<double> oh_trigtime;
    vector<double> oh_pe;

    int of_nFlash;
    vector<float> of_t;
    vector<float> of_peTotal;
    vector<int> of_multiplicity;
    TClonesArray* fPEperOpDet;

    int fRaw_nChannel;
    std::vector<int> fRaw_channelId;
    TClonesArray* fRaw_wf;

    int fSIMIDE_size;
    vector<int> fSIMIDE_channelIdY;
    vector<int> fSIMIDE_trackId;
    vector<unsigned short> fSIMIDE_tdc;
    vector<float> fSIMIDE_x;
    vector<float> fSIMIDE_y;
    vector<float> fSIMIDE_z;
    vector<float> fSIMIDE_numElectrons;

    int mc_Ntrack;                              // number of tracks in MC
    int mc_id[MAX_TRACKS];                      // track id; size == mc_Ntrack
    int mc_pdg[MAX_TRACKS];                     // track particle pdg; size == mc_Ntrack
    int mc_process[MAX_TRACKS];                 // track generation process code; size == mc_Ntrack
    int mc_mother[MAX_TRACKS];                  // mother id of this track; size == mc_Ntrack
    float mc_startXYZT[MAX_TRACKS][4];          // start position of this track; size == mc_Ntrack
    float mc_endXYZT[MAX_TRACKS][4];            // end position of this track; size == mc_Ntrack
    float mc_startMomentum[MAX_TRACKS][4];      // start momentum of this track; size == mc_Ntrack
    float mc_endMomentum[MAX_TRACKS][4];        // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int>> mc_daughters; // daughters id of this track; vector
    TObjArray* fMC_trackPosition;

    int mc_isnu;            // is neutrino interaction
    int mc_nGeniePrimaries; // number of Genie primaries
    int mc_nu_pdg;          // pdg code of neutrino
    int mc_nu_ccnc;         // cc or nc
    int mc_nu_mode;    // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
    int mc_nu_intType; // interaction type
    int mc_nu_target;  // target interaction
    int mc_hitnuc;     // hit nucleon
    int mc_hitquark;   // hit quark
    double mc_nu_Q2;   // Q^2
    double mc_nu_W;    // W
    double mc_nu_X;    // X
    double mc_nu_Y;    // Y
    double mc_nu_Pt;   // Pt
    double mc_nu_Theta; // angle relative to lepton
    float mc_nu_pos[4]; // interaction position of nu
    float mc_nu_mom[4]; // interaction momentum of nu

    // ----- derived ---
    std::map<int, int> trackIndex;
    std::vector<std::vector<int>> trackParents;
    std::vector<std::vector<int>> trackChildren;
    std::vector<std::vector<int>> trackSiblings;
    TDatabasePDG* dbPDG;

  }; // class CellTree

  //-----------------------------------------------------------------------
  CellTree::CellTree(fhicl::ParameterSet const& p) : EDAnalyzer(p)
  {
    dbPDG = new TDatabasePDG();
    entryNo = 0;

    fRawDigitLabel = p.get<std::string>("RawDigitLabel");
    fCalibLabel = p.get<std::string>("CalibLabel");
    fOpHitLabel = p.get<std::string>("OpHitLabel");
    fOpFlashLabel = p.get<std::string>("OpFlashLabel");
    fTriggerLabel = p.get<std::string>("TriggerLabel");
    fSimEnergyDepositLabel = p.get<std::string>("SimEnergyDepositLabel");
    fSpacePointLabels = p.get<std::vector<std::string>>("SpacePointLabels");
    fSimChannelLabel = p.get<std::string>("SimChannelLabel");
    fOutFileName = p.get<std::string>("outFile");
    mcOption = p.get<std::string>("mcOption");
    fSaveMCTrackPoints = p.get<bool>("saveMCTrackPoints");
    fSaveRaw = p.get<bool>("saveRaw");
    fSaveCalib = p.get<bool>("saveCalib");
    fSaveOpHit = p.get<bool>("saveOpHit");
    fSaveOpFlash = p.get<bool>("saveOpFlash");
    fSaveMC = p.get<bool>("saveMC");
    fSaveSimChannel = p.get<bool>("saveSimChannel");
    fSaveTrigger = p.get<bool>("saveTrigger");
    fSaveJSON = p.get<bool>("saveJSON");
    fT0_corrected = p.get<bool>("t0_corrected");
    opMultPEThresh = p.get<float>("opMultPEThresh");
    drift_speed = p.get<float>("drift_speed"); // mm/us
    nRawSamples = p.get<int>("nRawSamples");

    InitProcessMap();
    initOutput();
  }

  //-----------------------------------------------------------------------
  void CellTree::initOutput()
  {
    TDirectory* tmpDir = gDirectory;

    fOutFile = new TFile(fOutFileName.c_str(), "recreate");

    // 3.1: add mc_trackPosition
    TNamed version("version", "4.0");
    version.Write();

    // init Event TTree
    TDirectory* subDir = fOutFile->mkdir("Event");
    subDir->cd();
    fEventTree = new TTree("Sim", "Event Tree from Simulation");
    fEventTree->Branch("eventNo", &fEvent);
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);
    fEventTree->Branch("eventTime", &fEventTime); // timestamp

    fEventTree->Branch("triggerNo", &fTriggernumber);   // timestamp
    fEventTree->Branch("triggerTime", &fTriggertime);   // timestamp
    fEventTree->Branch("beamgateTime", &fBeamgatetime); // timestamp
    fEventTree->Branch("triggerBits", &fTriggerbits);   // timestamp

    fEventTree->Branch("raw_nChannel", &fRaw_nChannel);   // number of hit channels above threshold
    fEventTree->Branch("raw_channelId", &fRaw_channelId); // hit channel id; size == raw_nChannel
    fRaw_wf = new TClonesArray("TH1F");
    fEventTree->Branch("raw_wf", &fRaw_wf, 256000, 0); // raw waveform adc of each channel

    fEventTree->Branch("calib_nChannel",
                       &fCalib_nChannel); // number of hit channels above threshold
    fEventTree->Branch("calib_channelId", &fCalib_channelId); // hit channel id; size == calib_Nhit
    fCalib_wf = new TClonesArray("TH1F");
    fEventTree->Branch("calib_wf", &fCalib_wf, 256000, 0); // calib waveform adc of each channel
    // fCalib_wf->BypassStreamer();
    // fEventTree->Branch("calib_wfTDC", &fCalib_wfTDC);  // calib waveform tdc of each channel

    fEventTree->Branch("oh_nHits", &oh_nHits);     // number of op hits
    fEventTree->Branch("oh_channel", &oh_channel); //opchannel id; size == no ophits
    fEventTree->Branch("oh_bgtime",
                       &oh_bgtime); // optical pulse peak time w.r.t. start of beam gate
    fEventTree->Branch("oh_trigtime", &oh_trigtime); // optical pulse peak time w.r.t. trigger
    fEventTree->Branch("oh_pe", &oh_pe);             // pulse PE

    fEventTree->Branch("of_nFlash", &of_nFlash);
    fEventTree->Branch("of_t", &of_t);             // time in us w.r.t. the trigger for each flash
    fEventTree->Branch("of_peTotal", &of_peTotal); // total PE (sum of all PMTs) for each flash
    fEventTree->Branch("of_multiplicity",
                       &of_multiplicity); // total number of PMTs above threshold for each flash
    fPEperOpDet = new TClonesArray("TH1F");
    fEventTree->Branch("pe_opdet", &fPEperOpDet, 256000, 0);

    fEventTree->Branch("simide_size", &fSIMIDE_size); // size of stored sim:IDE
    fEventTree->Branch("simide_channelIdY", &fSIMIDE_channelIdY);
    fEventTree->Branch("simide_trackId", &fSIMIDE_trackId);
    fEventTree->Branch("simide_tdc", &fSIMIDE_tdc);
    fEventTree->Branch("simide_x", &fSIMIDE_x);
    fEventTree->Branch("simide_y", &fSIMIDE_y);
    fEventTree->Branch("simide_z", &fSIMIDE_z);
    fEventTree->Branch("simide_numElectrons", &fSIMIDE_numElectrons);

    fEventTree->Branch("mc_Ntrack", &mc_Ntrack);               // number of tracks in MC
    fEventTree->Branch("mc_id", &mc_id, "mc_id[mc_Ntrack]/I"); // track id; size == mc_Ntrack
    fEventTree->Branch(
      "mc_pdg", &mc_pdg, "mc_pdg[mc_Ntrack]/I"); // track particle pdg; size == mc_Ntrack
    fEventTree->Branch(
      "mc_process",
      &mc_process,
      "mc_process[mc_Ntrack]/I"); // track generation process code; size == mc_Ntrack
    fEventTree->Branch("mc_mother",
                       &mc_mother,
                       "mc_mother[mc_Ntrack]/I");      // mother id of this track; size == mc_Ntrack
    fEventTree->Branch("mc_daughters", &mc_daughters); // daughters id of this track; vector
    fEventTree->Branch(
      "mc_startXYZT",
      &mc_startXYZT,
      "mc_startXYZT[mc_Ntrack][4]/F"); // start position of this track; size == mc_Ntrack
    fEventTree->Branch(
      "mc_endXYZT",
      &mc_endXYZT,
      "mc_endXYZT[mc_Ntrack][4]/F"); // start position of this track; size == mc_Ntrack
    fEventTree->Branch(
      "mc_startMomentum",
      &mc_startMomentum,
      "mc_startMomentum[mc_Ntrack][4]/F"); // start momentum of this track; size == mc_Ntrack
    fEventTree->Branch(
      "mc_endMomentum",
      &mc_endMomentum,
      "mc_endMomentum[mc_Ntrack][4]/F"); // start momentum of this track; size == mc_Ntrack
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
      // gSystem->ChangeDirectory("bee");
      gSystem->MakeDirectory("bee/data");
    }
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
      gSystem->ChangeDirectory("bee");
      system("zip -r bee_upload data");
      gSystem->ChangeDirectory("..");
    }
  }

  //-----------------------------------------------------------------------
  void CellTree::beginRun(const art::Run& /*run*/) { mf::LogInfo("CellTree") << "begin run"; }

  //-----------------------------------------------------------------------
  void CellTree::analyze(const art::Event& event)
  {
    reset();
    fEvent = event.id().event();
    fRun = event.run();
    fSubRun = event.subRun();
    art::Timestamp ts = event.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    fEventTime = tts.AsDouble();

    if (fSaveRaw) processRaw(event);
    if (fSaveCalib) processCalib(event);
    if (fSaveOpHit) processOpHit(event);
    if (fSaveOpFlash) processOpFlash(event);
    if (fSaveSimChannel) processSimChannel(event);
    if (fSaveMC) processMC(event);
    if (fSaveTrigger) processTrigger(event);

    if (fSaveJSON) {
      gSystem->MakeDirectory(TString::Format("bee/data/%i", entryNo).Data());
      int nSp = fSpacePointLabels.size();
      for (int i = 0; i < nSp; i++) {
        TString jsonfile;
        jsonfile.Form("bee/data/%i/%i-%s.json", entryNo, entryNo, fSpacePointLabels[i].c_str());
        std::ofstream out(jsonfile.Data());
        if (fSpacePointLabels[i] == "truthDepo") {
          processSpacePointTruthDepo(event, fSpacePointLabels[i], out, fT0_corrected);
        }
        else {
          processSpacePoint(event, fSpacePointLabels[i], out);
        }
        out.close();
      }

      if (fSaveMC) {
        processMCTracks();
        TString jsonfile;
        jsonfile.Form("bee/data/%i/%i-mc.json", entryNo, entryNo);
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

    oh_channel.clear();
    oh_bgtime.clear();
    oh_trigtime.clear();
    oh_pe.clear();

    of_t.clear();
    of_peTotal.clear();
    of_multiplicity.clear();
    fPEperOpDet->Delete();

    fSIMIDE_channelIdY.clear();
    fSIMIDE_trackId.clear();
    fSIMIDE_tdc.clear();
    fSIMIDE_x.clear();
    fSIMIDE_y.clear();
    fSIMIDE_z.clear();
    fSIMIDE_numElectrons.clear();

    mc_Ntrack = 0;
    for (int i = 0; i < MAX_TRACKS; i++) {
      mc_id[i] = 0;
      mc_pdg[i] = 0;
      mc_mother[i] = 0;
      for (int j = 0; j < 4; j++) {
        mc_startXYZT[i][j] = 0;
        mc_endXYZT[i][j] = 0;
        mc_startMomentum[i][j] = 0;
        mc_endMomentum[i][j] = 0;
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
    for (int i = 0; i < 4; i++) {
      mc_nu_pos[i] = 0;
      mc_nu_mom[i] = 0;
    }

    trackIndex.clear();
    trackParents.clear();
    trackChildren.clear();
    trackSiblings.clear();
  }

  //-----------------------------------------------------------------------
  void CellTree::processRaw(const art::Event& event)
  {
    art::Handle<std::vector<raw::RawDigit>> rawdigit;
    if (!event.getByLabel(fRawDigitLabel, rawdigit)) {
      cout << "WARNING: no label " << fRawDigitLabel << endl;
      return;
    }
    std::vector<art::Ptr<raw::RawDigit>> wires;
    art::fill_ptr_vector(wires, rawdigit);

    fRaw_nChannel = wires.size();

    int i = 0;
    for (auto const& wire : wires) {
      int chanId = wire->Channel();
      fRaw_channelId.push_back(chanId);

      int nSamples = wire->Samples();
      std::vector<short> uncompressed(nSamples);
      raw::Uncompress(wire->ADCs(), uncompressed, wire->Compression());

      TH1F* h = new ((*fRaw_wf)[i]) TH1F("", "", nRawSamples, 0, nRawSamples);
      for (int j = 1; j <= nSamples; j++) {
        h->SetBinContent(j, uncompressed[j - 1]);
      }
      i++;
      if (i == 1) { cout << nSamples << " samples expanding to " << nRawSamples << endl; }
    }
  }

  //-----------------------------------------------------------------------
  void CellTree::processCalib(const art::Event& event)
  {

    art::Handle<std::vector<recob::Wire>> wires_handle;
    if (!event.getByLabel(fCalibLabel, wires_handle)) {
      cout << "WARNING: no label " << fCalibLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::Wire>> wires;
    art::fill_ptr_vector(wires, wires_handle);

    // wires size should == Nchannels == 1992; (no hit channel has a flat 0-waveform)
    // cout << "\n wires size: " << wires.size() << endl;
    fCalib_nChannel = wires.size();

    int i = 0;
    for (auto const& wire : wires) {
      std::vector<float> calibwf = wire->Signal();
      int chanId = wire->Channel();
      fCalib_channelId.push_back(chanId);
      TH1F* h = new ((*fCalib_wf)[i]) TH1F("", "", nRawSamples, 0, nRawSamples);
      for (int j = 1; j <= nRawSamples; j++) {
        h->SetBinContent(j, calibwf[j]);
      }
      // fCalib_wf.push_back(calibwf);
      // cout << chanId << ", " << nSamples << endl;
      i++;
    }
  }

  //----------------------------------------------------------------------
  void CellTree::processOpHit(const art::Event& event)
  {
    art::Handle<std::vector<recob::OpHit>> ophit_handle;
    if (!event.getByLabel(fOpHitLabel, ophit_handle)) {
      cout << "WARNING: no label " << fOpHitLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::OpHit>> ophits;
    art::fill_ptr_vector(ophits, ophit_handle);
    oh_nHits = (int)ophits.size();

    for (auto const& oh : ophits) {
      oh_channel.push_back(oh->OpChannel());
      oh_bgtime.push_back(oh->PeakTime());
      oh_trigtime.push_back(oh->PeakTimeAbs());
      oh_pe.push_back(oh->PE());
    }
  }

  //----------------------------------------------------------------------
  void CellTree::processOpFlash(const art::Event& event)
  {
    art::Handle<std::vector<recob::OpFlash>> flash_handle;
    if (!event.getByLabel(fOpFlashLabel, flash_handle)) {
      cout << "WARNING: no label " << fOpFlashLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::OpFlash>> flashes;
    art::fill_ptr_vector(flashes, flash_handle);
    of_nFlash = (int)flashes.size();

    int a = 0;
    int nOpDet = fGeometry->NOpDets();

    for (auto const& flash : flashes) {
      of_t.push_back(flash->Time());
      of_peTotal.push_back(flash->TotalPE());
      TH1F* h = new ((*fPEperOpDet)[a]) TH1F("", "", nOpDet, 0, nOpDet);

      int mult = 0;
      for (int i = 0; i < nOpDet; ++i) {
        if (flash->PE(i) >= opMultPEThresh) { mult++; }
        h->SetBinContent(i, flash->PE(i));
      }
      of_multiplicity.push_back(mult);
      a++;
    }
  }

  //-----------------------------------------------------------------------
  void CellTree::processSimChannel(const art::Event& event)
  {
    art::Handle<std::vector<sim::SimChannel>> simChannelHandle;
    // event.getByLabel("largeant", simChannelHandle);
    if (!event.getByLabel(fSimChannelLabel, simChannelHandle)) {
      cout << "WARNING: no label " << fSimChannelLabel << endl;
      return;
    }

    // cout << "total simChannel: " << (*simChannelHandle).size() << endl;
    fSIMIDE_size = 0;
    for (auto const& channel : (*simChannelHandle)) {
      auto channelNumber = channel.Channel();
      // cout << channelNumber << endl;
      // if (! (fGeometry->SignalType( channelNumber ) == geo::kCollection) ) {
      //     continue;
      // }
      auto const& timeSlices = channel.TDCIDEMap();
      for (auto const& timeSlice : timeSlices) {
        auto const& energyDeposits = timeSlice.second;
        for (auto const& energyDeposit : energyDeposits) {
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
  void CellTree::processMC(const art::Event& event)
  {
    art::Handle<std::vector<simb::MCParticle>> particleHandle;
    if (!event.getByLabel("largeant", particleHandle)) return;
    std::vector<art::Ptr<simb::MCParticle>> particles;
    art::fill_ptr_vector(particles, particleHandle);

    art::Handle<std::vector<sim::SimChannel>> simChannelHandle;
    event.getByLabel("largeant", simChannelHandle);

    // art::ServiceHandle<cheat::BackTracker const> bt;
    art::FindOneP<simb::MCTruth> fo(particleHandle, event, "largeant");

    int i = 0;     // track index in saved MCParticles;
    int i_all = 0; // track index in all MCParticles;
    for (auto const& particle : particles) {
      art::Ptr<simb::MCTruth> mctruth = fo.at(i_all);
      i_all++;

      if (mcOption == "nuOnly") {
        if (!(mctruth->Origin() == 1 && particle->Mother() == 0)) { continue; }
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
      // const art::Ptr<simb::MCTruth> mctruth = bt_serv->TrackIDToMCTruth(mc_id[i]);

      mc_process[i] = processMap[particle->Process()];
      if (mc_process[i] == 0) cout << "unknown process: " << particle->Process() << endl;
      mc_id[i] = particle->TrackId();
      mc_pdg[i] = particle->PdgCode();
      mc_mother[i] = particle->Mother();
      savedMCTrackIdMap[mc_id[i]] = mc_pdg[i];

      int Ndaughters = particle->NumberDaughters();
      vector<int> daughters;
      for (int i = 0; i < Ndaughters; i++) {
        daughters.push_back(particle->Daughter(i));
      }
      mc_daughters.push_back(daughters);
      size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
      int last = numberTrajectoryPoints - 1;
      const TLorentzVector& positionStart = particle->Position(0);
      const TLorentzVector& positionEnd = particle->Position(last);
      const TLorentzVector& momentumStart = particle->Momentum(0);
      const TLorentzVector& momentumEnd = particle->Momentum(last);
      positionStart.GetXYZT(mc_startXYZT[i]);
      positionEnd.GetXYZT(mc_endXYZT[i]);
      momentumStart.GetXYZT(mc_startMomentum[i]);
      momentumEnd.GetXYZT(mc_endMomentum[i]);

      if (fSaveMCTrackPoints) {
        TClonesArray* Lposition = new TClonesArray("TLorentzVector", numberTrajectoryPoints);
        // Read the position and momentum along this particle track
        for (unsigned int j = 0; j < numberTrajectoryPoints; j++) {
          new ((*Lposition)[j]) TLorentzVector(particle->Position(j));
        }
        fMC_trackPosition->Add(Lposition);
      }

      i++;
      if (i == MAX_TRACKS) {
        cout << "WARNING:: # tracks exceeded MAX_TRACKS " << MAX_TRACKS << endl;
        break;
      }
    } // particle loop done
    mc_Ntrack = i;
    // cout << "MC_Ntracks:" << mc_Ntrack << endl;

    // Generator Info
    art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
    event.getByLabel("generator", mctruthListHandle);
    std::vector<art::Ptr<simb::MCTruth>> mclist;
    art::fill_ptr_vector(mclist, mctruthListHandle);
    art::Ptr<simb::MCTruth> mctruth;

    if (mclist.size() > 0) {
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
  void CellTree::processSpacePoint(const art::Event& event, TString option, ostream& out)
  {

    art::Handle<std::vector<recob::SpacePoint>> sp_handle;
    art::Handle<std::vector<recob::PointCharge>> pc_handle;
    bool sp_exists = event.getByLabel(option.Data(), sp_handle);
    bool pc_exists = event.getByLabel(option.Data(), pc_handle);
    if (!sp_exists) {
      cout << "WARNING: no label " << option << endl;
      return;
    }
    std::vector<art::Ptr<recob::SpacePoint>> sps;
    std::vector<art::Ptr<recob::PointCharge>> pcs;
    art::fill_ptr_vector(sps, sp_handle);
    if (pc_exists) {
      art::fill_ptr_vector(pcs, pc_handle);
      if (sps.size() != pcs.size()) {
        cout << "WARNING: SpacePoint and PointCharge length mismatch" << endl;
        return;
      }
    }
    double x = 0, y = 0, z = 0, q = 0, nq = 1;
    vector<double> vx, vy, vz, vq, vnq;

    for (uint i = 0; i < sps.size(); i++) {
      // cout << sp->XYZ()[0] << ", " << sp->XYZ()[1] << ", " << sp->XYZ()[2] << endl;
      x = sps[i]->XYZ()[0];
      y = sps[i]->XYZ()[1];
      z = sps[i]->XYZ()[2];
      if (pc_exists && pcs[i]->hasCharge()) { q = pcs[i]->charge(); }
      else {
        q = 0;
      }
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
    else if (geomName.Contains("protodunevd")) {
      geomName = "protodunevd";
    }
    else if (geomName.Contains("protodune")) {
      geomName = "protodune";
    }
    else if (geomName.Contains("workspace")) {
      geomName = "dune10kt_workspace";
    }
    else if (geomName.Contains("icarus")) {
      geomName = "icarus";
    }
    else if (geomName.Contains("sbnd")) {
      geomName = "sbnd";
    }
    else {
      geomName = "uboone";
    } // use uboone as default
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

  //---- The X-axis position along drift changes to wire plane readout view without t0 correction ----
  void CellTree::processSpacePointTruthDepo(const art::Event& event,
                                            TString option,
                                            ostream& out,
                                            bool t0_corrected)
  {

    art::Handle<std::vector<sim::SimEnergyDeposit>> sed_handle;
    if (!event.getByLabel(fSimEnergyDepositLabel, sed_handle)) {
      cout << "WARNING: no label " << fSimEnergyDepositLabel << " for SimEnergyDeposit" << endl;
      return;
    }
    std::vector<art::Ptr<sim::SimEnergyDeposit>> sed;
    art::fill_ptr_vector(sed, sed_handle);
    int size = sed.size();
    double t = 0, x = 0, y = 0, z = 0, q = 0, nq = 1, cluster_id = 1;
    vector<double> vx, vy, vz, vq, vnq, vcluster;

    TString geomName(fGeometry->DetectorName().c_str());
    if (geomName.Contains("35t")) { geomName = "dune35t"; }
    else if (geomName.Contains("protodunevd")) {
      geomName = "protodunevd";
    }
    else if (geomName.Contains("protodune")) {
      geomName = "protodune";
    }
    else if (geomName.Contains("workspace")) {
      geomName = "dune10kt_workspace";
    }
    else if (geomName.Contains("icarus")) {
      geomName = "icarus";
    }
    else if (geomName.Contains("sbnd")) {
      geomName = "sbnd";
    }
    else {
      geomName = "uboone";
    } // use uboone as default

    for (int i = 0; i < size; i++) {
      // cout << sp->XYZ()[0] << ", " << sp->XYZ()[1] << ", " << sp->XYZ()[2] << endl;
      x = sed[i]->MidPointX(); // unit: cm
      t = sed[i]->Time();      // unit: ns
      y = sed[i]->MidPointY();
      z = sed[i]->MidPointZ();
      q = sed[i]->NumElectrons();
      if (q < 0) q = sed[i]->Energy() * 25000; // approx. #electrons
      // cout << q << ", " << sed[i]->Energy()*25000 << endl;
      if (q < 1000) continue; // skip small dots to reduce size
      if (!t0_corrected) {    // t0 ''enters'' in position along drift
        if (geomName == "sbnd") {
          if (x < 0) {
            x = x + t * 1e-3 * drift_speed * 0.1;
            cluster_id = 1;
          }
          else if (x > 0) {
            x = x - t * 1e-3 * drift_speed * 0.1;
            cluster_id = 2;
          }
          else {
            cluster_id = 3;
          }
          if (t * 1e-3 > 0 && t * 1e-3 < 5) cluster_id = 0;
        }
        else if (geomName == "uboone") {
          x = x + t * 1e-3 * drift_speed * 0.1;
        }
        else {
          cout << "t0 uncorrection for drift volume(s) yet to be added for " << geomName << endl;
        }
      }
      vx.push_back(x);
      vy.push_back(y);
      vz.push_back(z);
      vq.push_back(q);
      vnq.push_back(nq);
      vcluster.push_back(cluster_id);
    }

    out << fixed << setprecision(1);
    out << "{" << endl;

    out << '"' << "runNo" << '"' << ":" << '"' << fRun << '"' << "," << endl;
    out << '"' << "subRunNo" << '"' << ":" << '"' << fSubRun << '"' << "," << endl;
    out << '"' << "eventNo" << '"' << ":" << '"' << fEvent << '"' << "," << endl;

    out << '"' << "geom" << '"' << ":" << '"' << geomName << '"' << "," << endl;

    print_vector(out, vx, "x");
    print_vector(out, vy, "y");
    print_vector(out, vz, "z");

    out << fixed << setprecision(0);
    print_vector(out, vq, "q");
    print_vector(out, vnq, "nq");
    print_vector(out, vcluster, "cluster_id");

    out << '"' << "type" << '"' << ":" << '"' << option << '"' << endl;
    out << "}" << endl;
  }

  //-----------------------------------------------------------------------
  void CellTree::print_vector(ostream& out, vector<double>& v, TString desc, bool end)
  {
    int N = v.size();

    out << '"' << desc << '"' << ":[";
    for (int i = 0; i < N; i++) {
      out << v[i];
      if (i != N - 1) { out << ","; }
    }
    out << "]";
    if (!end) out << ",";
    out << endl;
  }

  //-----------------------------------------------------------------------
  void CellTree::processMCTracks()
  {
    // map track id to track index in the array
    for (int i = 0; i < mc_Ntrack; i++) {
      trackIndex[mc_id[i]] = i;
    }

    // in trackParents, trackChildren, trackSiblings vectors, store track index (not track id)
    for (int i = 0; i < mc_Ntrack; i++) {
      // currently, parent size == 1;
      // for primary particle, parent id = 0;
      vector<int> parents;
      if (!IsPrimary(i)) { parents.push_back(trackIndex[mc_mother[i]]); }
      trackParents.push_back(parents); // primary track will have 0 parents

      vector<int> children;
      int nChildren = mc_daughters.at(i).size();
      for (int j = 0; j < nChildren; j++) {
        children.push_back(trackIndex[mc_daughters.at(i).at(j)]);
      }
      trackChildren.push_back(children);
    }

    // siblings
    for (int i = 0; i < mc_Ntrack; i++) {
      vector<int> siblings;
      if (IsPrimary(i)) {
        for (int j = 0; j < mc_Ntrack; j++) {
          if (IsPrimary(j)) { siblings.push_back(j); }
        }
      }
      else {
        // siblings are simply children of the mother
        int mother = trackIndex[mc_mother[i]];
        int nSiblings = trackChildren.at(mother).size();
        for (int j = 0; j < nSiblings; j++) {
          siblings.push_back(trackChildren.at(mother).at(j));
        }
      }
      trackSiblings.push_back(siblings);
    }
  }

  //-----------------------------------------------------------------------
  void CellTree::processTrigger(const art::Event& event)
  {
    art::Handle<std::vector<raw::Trigger>> triggerListHandle;
    std::vector<art::Ptr<raw::Trigger>> triggerlist;
    if (event.getByLabel(fTriggerLabel, triggerListHandle)) {
      art::fill_ptr_vector(triggerlist, triggerListHandle);
    }
    else {
      cout << "WARNING: no trigger label " << fTriggerLabel << endl;
    }
    if (triggerlist.size()) {
      fTriggernumber = triggerlist[0]->TriggerNumber();
      fTriggertime = triggerlist[0]->TriggerTime();
      fBeamgatetime = triggerlist[0]->BeamGateTime();
      fTriggerbits = triggerlist[0]->TriggerBits();
    }
    else {
      fTriggernumber = 0;
      fTriggertime = 0;
      fBeamgatetime = 0;
      fTriggerbits = 0;
    }
  }

  //-----------------------------------------------------------------------
  bool CellTree::DumpMCJSON(int id, ostream& out)
  {
    int i = trackIndex[id];
    if (!KeepMC(i)) return false;

    int e = KE(mc_startMomentum[i]) * 1000;

    int nDaughter = mc_daughters.at(i).size();
    vector<int> saved_daughters;
    for (int j = 0; j < nDaughter; j++) {
      int daughter_id = mc_daughters.at(i).at(j);
      // int e_daughter = KE(mc_startMomentum[ trackIndex[daughter_id] ])*1000;
      // if (e_daughter >= thresh_KE) {
      if (KeepMC(trackIndex[daughter_id])) { saved_daughters.push_back(daughter_id); }
    }

    vector<double> vx, vy, vz;
    if (fSaveMCTrackPoints) {
      // fMC_trackPosition->Print();
      TClonesArray* traj = (TClonesArray*)(*fMC_trackPosition)[i];
      int nPoints = traj->GetEntries();
      // cout << "traj points: " << nPoints << endl;
      for (int j = 0; j < nPoints; j++) {
        TLorentzVector* pos = (TLorentzVector*)(*traj)[j];
        vx.push_back(pos->X());
        vy.push_back(pos->Y());
        vz.push_back(pos->Z());
      }
    }

    out << fixed << setprecision(1);
    out << "{";

    out << "\"id\":" << id << ",";
    out << "\"text\":"
        << "\"" << PDGName(mc_pdg[i]) << "  " << e << " MeV\",";
    out << "\"data\":{";
    print_vector(out, vx, "traj_x");
    print_vector(out, vy, "traj_y");
    print_vector(out, vz, "traj_z");
    out << "\"start\":[" << mc_startXYZT[i][0] << ", " << mc_startXYZT[i][1] << ", "
        << mc_startXYZT[i][2] << "],";
    out << "\"end\":[" << mc_endXYZT[i][0] << ", " << mc_endXYZT[i][1] << ", " << mc_endXYZT[i][2]
        << "]";
    out << "},";
    out << "\"children\":[";
    int nSavedDaughter = saved_daughters.size();
    if (nSavedDaughter == 0) {
      out << "],";
      out << "\"icon\":"
          << "\"jstree-file\"";
      out << "}";
      return true;
    }
    else {
      for (int j = 0; j < nSavedDaughter; j++) {
        DumpMCJSON(saved_daughters.at(j), out);
        if (j != nSavedDaughter - 1) { out << ","; }
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
    for (int i = 0; i < mc_Ntrack; i++) {
      if (IsPrimary(i)) {
        // int e = KE(mc_startMomentum[i])*1000;
        // if (e<thresh_KE) continue;
        if (KeepMC(i)) { primaries.push_back(i); }
      }
    }
    int size = primaries.size();
    // cout << size << endl;
    for (int i = 0; i < size; i++) {
      if (DumpMCJSON(mc_id[primaries[i]], out) && i != size - 1) { out << ", "; }
    }

    out << "]";
  }

  //-----------------------------------------------------------------------
  double CellTree::KE(float* momentum)
  {
    TLorentzVector particle(momentum);
    return particle.E() - particle.M();
  }

  //-----------------------------------------------------------------------
  bool CellTree::KeepMC(int i)
  {
    double e = KE(mc_startMomentum[i]) * 1000;
    double thresh_KE_em = 5.; // MeV
    double thresh_KE_np = 50; // MeV
    // cout << "pdg: " << mc_pdg[i] << ", KE: " << e << " MeV, process: " << mc_process[i] << endl;
    if (mc_process[i] == 8    // muIoni
        || mc_process[i] == 6 // eBrem
        || mc_process[i] == 5 // eIoni
    ) {
      return false; // skip those ionization and radiation electrons as there are too many to show.
    }

    if (mc_pdg[i] == 22 || mc_pdg[i] == 11 || mc_pdg[i] == -11) {
      if (e >= thresh_KE_em)
        return true;
      else
        return false;
    }
    else if (mc_pdg[i] == 2112 || mc_pdg[i] == 2212 || mc_pdg[i] > 1e9) {
      if (e >= thresh_KE_np)
        return true;
      else
        return false;
    }
    return true;
  }

  //-----------------------------------------------------------------------
  TString CellTree::PDGName(int pdg)
  {
    TParticlePDG* p = dbPDG->GetParticle(pdg);
    if (p == 0) {
      if (pdg > 1e9) {
        int z = (pdg - 1e9) / 10000;
        int a = (pdg - 1e9 - z * 1e4) / 10;
        TString name;
        if (z == 18)
          name = "Ar";

        else if (z == 17)
          name = "Cl";
        else if (z == 19)
          name = "Ca";
        else if (z == 16)
          name = "S";
        else if (z == 15)
          name = "P";
        else if (z == 14)
          name = "Si";
        else if (z == 1)
          name = "H";
        else if (z == 2)
          name = "He";

        else
          return Form("%i", pdg);
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

    for (int i = 0; i < mc_Ntrack; i++) {
      cout << "\n              id: " << mc_id[i];
      cout << "\n             pdg: " << mc_pdg[i];
      cout << "\n          mother: " << mc_mother[i];
      cout << "\n      Ndaughters: " << mc_daughters.at(i).size();
      cout << "\n      start XYZT: (" << mc_startXYZT[i][0] << ", " << mc_startXYZT[i][1] << ", "
           << mc_startXYZT[i][2] << ", " << mc_startXYZT[i][3] << ")";
      cout << "\n        end XYZT: (" << mc_endXYZT[i][0] << ", " << mc_endXYZT[i][1] << ", "
           << mc_endXYZT[i][2] << ", " << mc_endXYZT[i][3] << ")";
      cout << "\n  start momentum: (" << mc_startMomentum[i][0] << ", " << mc_startMomentum[i][1]
           << ", " << mc_startMomentum[i][2] << ", " << mc_startMomentum[i][3] << ")";
      cout << "\n    end momentum: (" << mc_endMomentum[i][0] << ", " << mc_endMomentum[i][1]
           << ", " << mc_endMomentum[i][2] << ", " << mc_endMomentum[i][3] << ")";

      cout << endl;
    }
  }

  //-----------------------------------------------------------------------
  void CellTree::InitProcessMap()
  {
    processMap["unknown"] = 0;
    processMap["primary"] = 1;
    processMap["compt"] = 2;
    processMap["phot"] = 3;
    processMap["annihil"] = 4;
    processMap["eIoni"] = 5;
    processMap["eBrem"] = 6;
    processMap["conv"] = 7;
    processMap["muIoni"] = 8;
    processMap["muMinusCaptureAtRest"] = 9;
    processMap["neutronInelastic"] = 10;
    processMap["nCapture"] = 11;
    processMap["hadElastic"] = 12;
    processMap["Decay"] = 13;
    processMap["CoulombScat"] = 14;
    processMap["muPairProd"] = 15;
    processMap["muBrems"] = 16;
    processMap["muPairProd"] = 17;
    processMap["PhotonInelastic"] = 18;
    processMap["hIoni"] = 19;
    processMap["protonInelastic"] = 20;
    processMap["pi+Inelastic"] = 21;
    processMap["CHIPSNuclearCaptureAtRest"] = 22;
    processMap["pi-Inelastic"] = 23;
  }

  //-----------------------------------------------------------------------
  DEFINE_ART_MODULE(CellTree)
} // namespace microboone
