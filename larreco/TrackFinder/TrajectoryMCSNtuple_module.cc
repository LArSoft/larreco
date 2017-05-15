////////////////////////////////////////////////////////////////////////
// Class:       TrajectoryMCSNtuple
// Plugin Type: analyzer (art v2_05_00)
// File:        TrajectoryMCSNtuple_module.cc
//
// Generated at Mon Feb  6 10:06:04 2017 by Giuseppe Cerati using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "lardata/RecoObjects/TrackState.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class TrajectoryMCSNtuple;

class TrajectoryMCSNtuple : public art::EDAnalyzer {
public:

  struct Inputs {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<std::string> inputLabel {
      Name("inputLabel"),
      Comment("Label of recob::TrackTrajectory Collection to be fit")
     };
  };

  struct Config {
    using Name = fhicl::Name;
    fhicl::Table<TrajectoryMCSNtuple::Inputs> inputs {
      Name("inputs"),
    };
    fhicl::Table<trkf::TrajectoryMCSFitter::Config> mcsfitter {
      Name("mcsfitter")
    };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit TrajectoryMCSNtuple(Parameters const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  ~TrajectoryMCSNtuple();

  // Plugins should not be copied or assigned.
  TrajectoryMCSNtuple(TrajectoryMCSNtuple const &) = delete;
  TrajectoryMCSNtuple(TrajectoryMCSNtuple &&) = delete;
  TrajectoryMCSNtuple & operator = (TrajectoryMCSNtuple const &) = delete;
  TrajectoryMCSNtuple & operator = (TrajectoryMCSNtuple &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  void beginJob() override;

  void resetTree();

private:

  std::string inputTracksLabel;
  trkf::TrackMomentumCalculator tmc;
  trkf::TrajectoryMCSFitter mcsfitter;
  TTree* tree;
  //
  int    run, subrun, eventid;
  double trkLength;
  double trkMom_MuFwd, trkMomErr_MuFwd, trkLL_MuFwd;
  double trkMom_MuBwd, trkMomErr_MuBwd, trkLL_MuBwd;
  double trkMom_Mu, trkMomErr_Mu, trkLL_Mu;
  double trkDeltaLL_Mu;
  int    trkIsBestFwd_Mu;
  double trkMom_PFwd, trkMomErr_PFwd, trkLL_PFwd;
  double trkMom_PBwd, trkMomErr_PBwd, trkLL_PBwd;
  double trkMom_P, trkMomErr_P, trkLL_P;
  double trkDeltaLL_P;
  int    trkIsBestFwd_P;
  std::vector<double> trkSegRadLengths;
  std::vector<double> trkScattAngles;
  //
  int    trkNHits, trkIsContained;
  double trkMom_RangeMu, trkMom_RangeP;
  int    trkID_MCSRange;
  double trkStartPosX, trkStartPosY, trkStartPosZ;
  double trkStartDirX, trkStartDirY, trkStartDirZ;
  double trkEndPosX, trkEndPosY, trkEndPosZ;
  double trkEndDirX, trkEndDirY, trkEndDirZ;
  std::vector<int> trkSegNHits;
  std::vector<double> trkSegStartPosX, trkSegStartPosY, trkSegStartPosZ;
  std::vector<double> trkSegEndPosX, trkSegEndPosY, trkSegEndPosZ;
  std::vector<double> trkSegDirX, trkSegDirY, trkSegDirZ;
  //
  int    trkid0;
  int    trkid1;
  int    trkid2;
  double trkpida0;
  double trkpida1;
  double trkpida2;
  //
  double simMom;
  double simLength;
  double simStartPosX, simStartPosY, simStartPosZ;
  double simStartDirX, simStartDirY, simStartDirZ;
  double simEndPosX, simEndPosY, simEndPosZ;
  double simEndDirX, simEndDirY, simEndDirZ;
  int    simID;
  std::string simProc;
  int    simIsContained;
  int    simAndTrkSameDir;
  //
};

void TrajectoryMCSNtuple::resetTree() {
  run = -999;
  subrun = -999;
  eventid = -999;
  trkLength = -999;
  trkMom_MuFwd = -999;
  trkMomErr_MuFwd = -999;
  trkLL_MuFwd = -999;
  trkMom_MuBwd = -999;
  trkMomErr_MuBwd = -999;
  trkLL_MuBwd = -999;
  trkMom_Mu = -999;
  trkMomErr_Mu = -999;
  trkLL_Mu = -999;
  trkDeltaLL_Mu = -999;
  trkIsBestFwd_Mu = -999;
  trkMom_PFwd = -999;
  trkMomErr_PFwd = -999;
  trkLL_PFwd = -999;
  trkMom_PBwd = -999;
  trkMomErr_PBwd = -999;
  trkLL_PBwd = -999;
  trkMom_P = -999;
  trkMomErr_P = -999;
  trkLL_P = -999;
  trkDeltaLL_P = -999;
  trkIsBestFwd_P = -999;
  trkSegRadLengths.clear();
  trkScattAngles.clear();
  //
  trkNHits = -999;
  trkIsContained = -999;
  trkMom_RangeMu = -999;
  trkMom_RangeP = -999;
  trkID_MCSRange = -999;
  trkStartPosX = -999;
  trkStartPosY = -999;
  trkStartPosZ = -999;
  trkStartDirX = -999;
  trkStartDirY = -999;
  trkStartDirZ = -999;
  trkEndPosX = -999;
  trkEndPosY = -999;
  trkEndPosZ = -999;
  trkEndDirX = -999;
  trkEndDirY = -999;
  trkEndDirZ = -999;
  trkSegNHits.clear();
  trkSegStartPosX.clear();
  trkSegStartPosY.clear();
  trkSegStartPosZ.clear();
  trkSegEndPosX.clear();
  trkSegEndPosY.clear();
  trkSegEndPosZ.clear();
  trkSegDirX.clear();
  trkSegDirY.clear();
  trkSegDirZ.clear();
  //
  trkid0 = -999;
  trkid1 = -999;
  trkid2 = -999;
  trkpida0 = -999;
  trkpida1 = -999;
  trkpida2 = -999;
  //
  simMom = -999;
  simLength = -999;
  simStartPosX = -999;
  simStartPosY = -999;
  simStartPosZ = -999;
  simStartDirX = -999;
  simStartDirY = -999;
  simStartDirZ = -999;
  simEndPosX = -999;
  simEndPosY = -999;
  simEndPosZ = -999;
  simEndDirX = -999;
  simEndDirY = -999;
  simEndDirZ = -999;
  simID = -999;
  simProc = "";
  simIsContained = -999;
  simAndTrkSameDir = -999;
}

void TrajectoryMCSNtuple::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //
  tree = tfs->make<TTree>("tree", "tree");
  //
  tree->Branch("run", &run,"run/I");
  tree->Branch("subrun", &subrun, "subrun/I");
  tree->Branch("eventid", &eventid, "eventid/I");
  //
  tree->Branch("trkLength", &trkLength, "trkLength/D");
  tree->Branch("trkMom_MuFwd"   , &trkMom_MuFwd   , "trkMom_MuFwd/D"   );
  tree->Branch("trkMomErr_MuFwd", &trkMomErr_MuFwd, "trkMomErr_MuFwd/D");
  tree->Branch("trkLL_MuFwd"    , &trkLL_MuFwd    , "trkLL_MuFwd/D"    );
  tree->Branch("trkMom_MuBwd"   , &trkMom_MuBwd   , "trkMom_MuBwd/D"   );
  tree->Branch("trkMomErr_MuBwd", &trkMomErr_MuBwd, "trkMomErr_MuBwd/D");
  tree->Branch("trkLL_MuBwd"    , &trkLL_MuBwd    , "trkLL_MuBwd/D"    );
  tree->Branch("trkMom_Mu"      , &trkMom_Mu      , "trkMom_Mu/D"      );
  tree->Branch("trkMomErr_Mu"   , &trkMomErr_Mu   , "trkMomErr_Mu/D"   );
  tree->Branch("trkLL_Mu"       , &trkLL_Mu       , "trkLL_Mu/D"       );
  tree->Branch("trkDeltaLL_Mu"  , &trkDeltaLL_Mu  , "trkDeltaLL_Mu/D"  );
  tree->Branch("trkIsBestFwd_Mu", &trkIsBestFwd_Mu, "trkIsBestFwd_Mu/I");
  //
  tree->Branch("trkMom_PFwd"   , &trkMom_PFwd   , "trkMom_PFwd/D"   );
  tree->Branch("trkMomErr_PFwd", &trkMomErr_PFwd, "trkMomErr_PFwd/D");
  tree->Branch("trkLL_PFwd"    , &trkLL_PFwd    , "trkLL_PFwd/D"    );
  tree->Branch("trkMom_PBwd"   , &trkMom_PBwd   , "trkMom_PBwd/D"   );
  tree->Branch("trkMomErr_PBwd", &trkMomErr_PBwd, "trkMomErr_PBwd/D");
  tree->Branch("trkLL_PBwd"    , &trkLL_PBwd    , "trkLL_PBwd/D"    );
  tree->Branch("trkMom_P"      , &trkMom_P      , "trkMom_P/D"      );
  tree->Branch("trkMomErr_P"   , &trkMomErr_P   , "trkMomErr_P/D"   );
  tree->Branch("trkLL_P"       , &trkLL_P       , "trkLL_P/D"       );
  tree->Branch("trkDeltaLL_P"  , &trkDeltaLL_P  , "trkDeltaLL_P/D"  );
  tree->Branch("trkIsBestFwd_P", &trkIsBestFwd_P, "trkIsBestFwd_P/I");
  //
  tree->Branch("trkSegRadLengths", &trkSegRadLengths);
  tree->Branch("trkScattAngles"  , &trkScattAngles  );
  //
  tree->Branch("trkNHits", &trkNHits, "trkNHits/I");
  tree->Branch("trkIsContained", &trkIsContained, "trkIsContained/I");
  tree->Branch("trkMom_RangeMu", &trkMom_RangeMu, "trkMom_RangeMu/D");
  tree->Branch("trkMom_RangeP" , &trkMom_RangeP , "trkMom_RangeP/D" );
  tree->Branch("trkID_MCSRange", &trkID_MCSRange, "trkID_MCSRange/I");
  tree->Branch("trkStartPosX", &trkStartPosX, "trkStartPosX/D");
  tree->Branch("trkStartPosY", &trkStartPosY, "trkStartPosY/D");
  tree->Branch("trkStartPosZ", &trkStartPosZ, "trkStartPosZ/D");
  tree->Branch("trkStartDirX", &trkStartDirX, "trkStartDirX/D");
  tree->Branch("trkStartDirY", &trkStartDirY, "trkStartDirY/D");
  tree->Branch("trkStartDirZ", &trkStartDirZ, "trkStartDirZ/D");
  tree->Branch("trkEndPosX"  , &trkEndPosX  , "trkEndPosX/D"  );
  tree->Branch("trkEndPosY"  , &trkEndPosY  , "trkEndPosY/D"  );
  tree->Branch("trkEndPosZ"  , &trkEndPosZ  , "trkEndPosZ/D"  );
  tree->Branch("trkEndDirX"  , &trkEndDirX  , "trkEndDirX/D"  );
  tree->Branch("trkEndDirY"  , &trkEndDirY  , "trkEndDirY/D"  );
  tree->Branch("trkEndDirZ"  , &trkEndDirZ  , "trkEndDirZ/D"  );
  //
  tree->Branch("trkSegNHits"    , &trkSegNHits    );
  tree->Branch("trkSegStartPosX", &trkSegStartPosX);
  tree->Branch("trkSegStartPosY", &trkSegStartPosY);
  tree->Branch("trkSegStartPosZ", &trkSegStartPosZ);
  tree->Branch("trkSegEndPosX"  , &trkSegEndPosX  );
  tree->Branch("trkSegEndPosY"  , &trkSegEndPosY  );
  tree->Branch("trkSegEndPosZ"  , &trkSegEndPosZ  );
  tree->Branch("trkSegDirX"     , &trkSegDirX	  );
  tree->Branch("trkSegDirY"     , &trkSegDirY	  );
  tree->Branch("trkSegDirZ"     , &trkSegDirZ     );
  //
  tree->Branch("trkid0"  , &trkid0  , "trkid0/I"  );
  tree->Branch("trkid1"  , &trkid1  , "trkid1/I"  );
  tree->Branch("trkid2"  , &trkid2  , "trkid2/I"  );
  tree->Branch("trkpida0", &trkpida0, "trkpida0/D");
  tree->Branch("trkpida1", &trkpida1, "trkpida1/D");
  tree->Branch("trkpida2", &trkpida2, "trkpida2/D");
  //
  tree->Branch("simMom"   , &simMom   , "simMom/D"   );
  tree->Branch("simLength", &simLength, "simLength/D");
  tree->Branch("simStartPosX", &simStartPosX, "simStartPosX/D");
  tree->Branch("simStartPosY", &simStartPosY, "simStartPosY/D");
  tree->Branch("simStartPosZ", &simStartPosZ, "simStartPosZ/D");
  tree->Branch("simStartDirX", &simStartDirX, "simStartDirX/D");
  tree->Branch("simStartDirY", &simStartDirY, "simStartDirY/D");
  tree->Branch("simStartDirZ", &simStartDirZ, "simStartDirZ/D");
  tree->Branch("simEndPosX"  , &simEndPosX  , "simEndPosX/D"  );
  tree->Branch("simEndPosY"  , &simEndPosY  , "simEndPosY/D"  );
  tree->Branch("simEndPosZ"  , &simEndPosZ  , "simEndPosZ/D"  );
  tree->Branch("simEndDirX"  , &simEndDirX  , "simEndDirX/D"  );
  tree->Branch("simEndDirY"  , &simEndDirY  , "simEndDirY/D"  );
  tree->Branch("simEndDirZ"  , &simEndDirZ  , "simEndDirZ/D"  );
  tree->Branch("simID", &simID, "simID/I");
  tree->Branch("simProc", &simProc);
  tree->Branch("simIsContained", &simIsContained, "simIsContained/I");
  tree->Branch("simAndTrkSameDir", &simAndTrkSameDir, "simAndTrkSameDir/I");
  //
}

TrajectoryMCSNtuple::TrajectoryMCSNtuple(Parameters const & p)
  : EDAnalyzer(p), inputTracksLabel(p().inputs().inputLabel()), mcsfitter(p().mcsfitter) {}

TrajectoryMCSNtuple::~TrajectoryMCSNtuple() {}

void TrajectoryMCSNtuple::analyze(art::Event const & e)
{

  using namespace std;
  using namespace trkf;
  using namespace recob;
  using namespace recob::tracking;

  run = e.run();
  subrun = e.subRun();
  eventid = e.event();

  art::ValidHandle<std::vector<sim::MCTrack> >* simTracks = 0;
  if (e.isRealData()==0) {
    art::InputTag SimTrackInputTag("mcreco");
    simTracks = new art::ValidHandle<std::vector<sim::MCTrack> >(e.getValidHandle<std::vector<sim::MCTrack> >(SimTrackInputTag));
  }

  art::InputTag TrackInputTag(inputTracksLabel);
  art::ValidHandle<std::vector<recob::Track> > Tracks = e.getValidHandle<std::vector<recob::Track> >(TrackInputTag);
  //
  art::InputTag MuMCSInputTag(inputTracksLabel+"MCSFitMu");
  art::ValidHandle<std::vector<recob::MCSFitResult> > MCSMu = e.getValidHandle<std::vector<recob::MCSFitResult> >(MuMCSInputTag);
  //
  art::InputTag PMCSInputTag(inputTracksLabel+"MCSFitP");
  art::ValidHandle<std::vector<recob::MCSFitResult> > MCSP = e.getValidHandle<std::vector<recob::MCSFitResult> >(PMCSInputTag);
  //
  art::InputTag PidInputTag(inputTracksLabel+"pid");
  art::FindMany<anab::ParticleID> AssPid(Tracks, e, PidInputTag);
  //
  assert(Tracks->size()==MCSMu->size());
  //
  for (unsigned int iTrack = 0; iTrack < Tracks->size(); ++iTrack) {
    const recob::Track& track = Tracks->at(iTrack);
    const recob::MCSFitResult& mcsMu = MCSMu->at(iTrack);
    const recob::MCSFitResult& mcsP = MCSP->at(iTrack);
    std::vector<const anab::ParticleID*> pids = AssPid.at(iTrack);
    //
    resetTree();
    //
    trkLength = track.Length();
    //
    trkMom_MuFwd    = mcsMu.fwdMomentum();
    trkMomErr_MuFwd = mcsMu.fwdMomUncertainty();
    trkLL_MuFwd     = mcsMu.fwdLogLikelihood();
    trkMom_MuBwd    = mcsMu.bwdMomentum();
    trkMomErr_MuBwd = mcsMu.bwdMomUncertainty();
    trkLL_MuBwd     = mcsMu.bwdLogLikelihood();
    trkMom_Mu       = mcsMu.bestMomentum();
    trkMomErr_Mu    = mcsMu.bestMomUncertainty();
    trkLL_Mu        = mcsMu.bestLogLikelihood();
    trkDeltaLL_Mu   = mcsMu.deltaLogLikelihood();
    trkIsBestFwd_Mu = mcsMu.isBestFwd();
    //
    trkMom_PFwd    = mcsP.fwdMomentum();
    trkMomErr_PFwd = mcsP.fwdMomUncertainty();
    trkLL_PFwd     = mcsP.fwdLogLikelihood();
    trkMom_PBwd    = mcsP.bwdMomentum();
    trkMomErr_PBwd = mcsP.bwdMomUncertainty();
    trkLL_PBwd     = mcsP.bwdLogLikelihood();
    trkMom_P       = mcsP.bestMomentum();
    trkMomErr_P    = mcsP.bestMomUncertainty();
    trkLL_P        = mcsP.bestLogLikelihood();
    trkDeltaLL_P   = mcsP.deltaLogLikelihood();
    trkIsBestFwd_P = mcsP.isBestFwd();
    //
    trkSegRadLengths = mcsMu.segmentRadLengths();
    trkScattAngles   = mcsMu.scatterAngles();
    //
    trkNHits = track.NPoints();
    //
    Point_t start(track.Start().X(),track.Start().Y(),track.Start().Z());
    Point_t end(track.End().X(),track.End().Y(),track.End().Z());
    trkIsContained = (start.X()>30.  && start.X()<230.  && end.X()>30.  && end.X()<230. &&
		      start.Y()>-85. && start.Y()<85.   && end.Y()>-85. && end.Y()<85.  &&
		      start.Z()>30.  && start.Z()<1010. && end.Z()>30.  && end.Z()<1010.);
    trkMom_RangeMu = tmc.GetTrackMomentum(trkLength,13  );
    trkMom_RangeP  = tmc.GetTrackMomentum(trkLength,2212);
    trkID_MCSRange = ( ( fabs(trkMom_RangeMu-trkMom_Mu)<fabs(trkMom_RangeP-trkMom_P) ) ? 13 : 2212);
    trkStartPosX = start.X();
    trkStartPosY = start.Y();
    trkStartPosZ = start.Z();
    trkStartDirX = track.StartDirection().X();
    trkStartDirY = track.StartDirection().Y();
    trkStartDirZ = track.StartDirection().Z();
    trkEndPosX = end.X();
    trkEndPosY = end.Y();
    trkEndPosZ = end.Z();
    trkEndDirX = track.EndDirection().X();
    trkEndDirY = track.EndDirection().Y();
    trkEndDirZ = track.EndDirection().Z();
    vector<size_t> breakpoints;
    vector<double> segradlengths;
    vector<double> cumseglens;
    mcsfitter.breakTrajInSegments(track.Trajectory(), breakpoints, segradlengths, cumseglens);
    //
    std::vector<int>    trkSegNHitsTmp;
    std::vector<double> trkSegStartPosXtmp;
    std::vector<double> trkSegStartPosYtmp;
    std::vector<double> trkSegStartPosZtmp;
    std::vector<double> trkSegEndPosXtmp;
    std::vector<double> trkSegEndPosYtmp;
    std::vector<double> trkSegEndPosZtmp;
    std::vector<double> trkSegDirXtmp;
    std::vector<double> trkSegDirYtmp;
    std::vector<double> trkSegDirZtmp;
    Vector_t pcdir;
    for (unsigned int p = 0; p<segradlengths.size(); p++) {
      mcsfitter.linearRegression(track.Trajectory(), breakpoints[p], breakpoints[p+1], pcdir);
      //
      trkSegNHitsTmp.push_back(breakpoints[p+1]-breakpoints[p]);
      trkSegStartPosXtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p]).X() );
      trkSegStartPosYtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p]).Y() );
      trkSegStartPosZtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p]).Z() );
      trkSegEndPosXtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p+1]-1).X());
      trkSegEndPosYtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p+1]-1).Y());
      trkSegEndPosZtmp.push_back( track.Trajectory().LocationAtPoint(breakpoints[p+1]-1).Z());
      //
      trkSegDirXtmp.push_back(pcdir.X());
      trkSegDirYtmp.push_back(pcdir.Y());
      trkSegDirZtmp.push_back(pcdir.Z());
    }
    trkSegNHits     = trkSegNHitsTmp;
    trkSegStartPosX = trkSegStartPosXtmp;
    trkSegStartPosY = trkSegStartPosYtmp;
    trkSegStartPosZ = trkSegStartPosZtmp;
    trkSegEndPosX   = trkSegEndPosXtmp;
    trkSegEndPosY   = trkSegEndPosYtmp;
    trkSegEndPosZ   = trkSegEndPosZtmp;
    trkSegDirX      = trkSegDirXtmp;
    trkSegDirY      = trkSegDirYtmp;
    trkSegDirZ      = trkSegDirZtmp;
    //
    for (size_t ipid = 0; ipid < pids.size(); ++ipid){
      if (!pids[ipid]->PlaneID().isValid) continue;
      int planenum = pids[ipid]->PlaneID().Plane;
      if (planenum<0||planenum>2) continue;
      if (planenum==0) {
	trkid0   = pids[ipid]->Pdg();
	trkpida0 = pids[ipid]->PIDA();
      }
      if (planenum==1) {
	trkid1   = pids[ipid]->Pdg();
	trkpida1 = pids[ipid]->PIDA();
      }
      if (planenum==2) {
	trkid2   = pids[ipid]->Pdg();
	trkpida2 = pids[ipid]->PIDA();
      }
    }
    //
    if (e.isRealData()==0) {
      for (unsigned int iMC = 0; iMC < (*simTracks)->size(); ++iMC) {
	const sim::MCTrack& mctrack = (*simTracks)->at(iMC);
	//
	Vector_t mcstartmom(mctrack.Start().Momentum().X()*0.001,mctrack.Start().Momentum().Y()*0.001,mctrack.Start().Momentum().Z()*0.001);
	Vector_t mcstartdir = mcstartmom.Unit();
	Vector_t mcendmom(mctrack.End().Momentum().X()*0.001,mctrack.End().Momentum().Y()*0.001,mctrack.End().Momentum().Z()*0.001);
	Vector_t mcenddir = mcendmom.Unit();
	double dotvtx = track.VertexDirection().X()*mcstartdir.X()+track.VertexDirection().Y()*mcstartdir.Y()+track.VertexDirection().Z()*mcstartdir.Z();
	double dotend = track.EndDirection().X()*mcstartdir.X()+track.EndDirection().Y()*mcstartdir.Y()+track.EndDirection().Z()*mcstartdir.Z();
	//
	Point_t mcstart(mctrack.Start().Position().X(),mctrack.Start().Position().Y(),mctrack.Start().Position().Z());
	Point_t mcend(mctrack.End().Position().X(),mctrack.End().Position().Y(),mctrack.End().Position().Z());
	//
	bool match = false;
	if ( (mcstart-start).R()<5 && dotvtx>0.9 )  match = true;
	if ( (mcstart-end).R()<5   && dotend<-0.9 ) match = true;
	if (!match) continue;
	//
	bool mccontained = (mctrack.Start().Position().X()>30.  && mctrack.Start().Position().X()<230.  && mctrack.End().Position().X()>30.  && mctrack.End().Position().X()<230. &&
			    mctrack.Start().Position().Y()>-85. && mctrack.Start().Position().Y()<85.   && mctrack.End().Position().Y()>-85. && mctrack.End().Position().Y()<85.  &&
			    mctrack.Start().Position().Z()>30.  && mctrack.Start().Position().Z()<1010. && mctrack.End().Position().Z()>30.  && mctrack.End().Position().Z()<1010.);
	//
	double mclen = 0.;
	for (unsigned int imc=0; imc<mctrack.size(); ++imc) {
	  if (imc>0) {
	    mclen+=sqrt( (mctrack[imc].X()-mctrack[imc-1].X())*(mctrack[imc].X()-mctrack[imc-1].X()) +
			 (mctrack[imc].Y()-mctrack[imc-1].Y())*(mctrack[imc].Y()-mctrack[imc-1].Y()) +
			 (mctrack[imc].Z()-mctrack[imc-1].Z())*(mctrack[imc].Z()-mctrack[imc-1].Z()) );
	  }
	}
	//
	simMom = mcstartmom.R();
	simLength = mclen;
	simStartPosX = mctrack.Start().Position().X();
	simStartPosY = mctrack.Start().Position().Y();
	simStartPosZ = mctrack.Start().Position().Z();
	simStartDirX = mcstartdir.X();
	simStartDirY = mcstartdir.Y();
	simStartDirZ = mcstartdir.Z();
	simEndPosX = mctrack.End().Position().X();
	simEndPosY = mctrack.End().Position().Y();
	simEndPosZ = mctrack.End().Position().Z();
	simEndDirX = mcenddir.X();
	simEndDirY = mcenddir.Y();
	simEndDirZ = mcenddir.Z();
	simID = std::abs(mctrack.PdgCode());
	simProc = mctrack.Process();
	simIsContained = mccontained;
	simAndTrkSameDir = dotvtx>0;
	break;
      }
    }
    //
    tree->Fill();
    //
  }
  if (e.isRealData()==0) delete simTracks;
}

DEFINE_ART_MODULE(TrajectoryMCSNtuple)
