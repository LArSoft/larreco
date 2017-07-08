//////////////////////////////////////////////////////////////////////
///
/// Step crawling code used by TrajClusterAlg
///
/// Bruce Baller, baller@fnal.gov
///
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrajClusterAlg.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"

class TH1F;
class TH2F;
class TProfile;

struct SortEntry{
  unsigned int index;
  float length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}


namespace tca {

  //------------------------------------------------------------------------------

  TrajClusterAlg::TrajClusterAlg(fhicl::ParameterSet const& pset):fCaloAlg(pset.get<fhicl::ParameterSet>("CaloAlg"))
  {
    reconfigure(pset);

    art::ServiceHandle<art::TFileService> tfs;
    
    // True - Reco vertex difference
    fNuVtx_dx = tfs->make<TH1F>("Vtx dx","Vtx dx",80,-10,10);
    fNuVtx_dy = tfs->make<TH1F>("Vtx dy","Vtx dy",80,-10,10);
    fNuVtx_dz = tfs->make<TH1F>("Vtx dz","Vtx dz",80,-10,10);
    fNuVtx_Score = tfs->make<TH1F>("Reco-True Vtx Score","Vtx Score",80, 0, 80);
    fNuVtx_Enu_Score_p = tfs->make<TProfile>("NuVtx_E_Score_p","Score vs Enu (MeV)", 20, 0, 2000);
    
    fVx2_Score = tfs->make<TH1F>("Vx2_Score","Vx2 Score",80, 0, 80);
    fVx3_Score = tfs->make<TH1F>("Vx3_Score","Vx2 Score",100, 0, 100);
    
    
    fdWire[0] = tfs->make<TH1F>("dWireEl","dWire - Electrons",21,-10,10);
    fdWire[1] = tfs->make<TH1F>("dWireMu","dWire - Muons",21,-10,10);
    fdWire[2] = tfs->make<TH1F>("dWirePi","dWire - Pions",21,-10,10);
    fdWire[3] = tfs->make<TH1F>("dWireKa","dWire - Kaons",21,-10,10);
    fdWire[4] = tfs->make<TH1F>("dWirePr","dWire - Protons",21,-10,10);
    
    fEP_T[0] = tfs->make<TProfile>("EP_T_El","EP vs T(MeV) - Electrons", 20, 0, 100);
    fEP_T[1] = tfs->make<TProfile>("EP_T_Mu","EP vs T(MeV) - Muons", 20, 0, 1000);
    fEP_T[2] = tfs->make<TProfile>("EP_T_Pi","EP vs T(MeV) - Pions", 20, 0, 1000);
    fEP_T[3] = tfs->make<TProfile>("EP_T_Ka","EP vs T(MeV) - Kaons", 20, 0, 1000);
    fEP_T[4] = tfs->make<TProfile>("EP_T_Pr","EP vs T(MeV) - Protons", 20, 0, 1000);
    

    fMCSMom_TruMom_e = tfs->make<TH2F>("MCSMom_TruMom_e","MCSMom vs Tru Mom electrons", 50, 0, 100, 50, 0, 1000);
    fMCSMom_TruMom_mu = tfs->make<TH2F>("MCSMom_TruMom_mu","MCSMom vs Tru Mom electrons", 50, 0, 1000, 50, 0, 1000);
    fMCSMom_TruMom_pi = tfs->make<TH2F>("MCSMom_TruMom_pi","MCSMom vs Tru Mom electrons", 50, 0, 1000, 50, 0, 1000);
    fMCSMom_TruMom_p = tfs->make<TH2F>("MCSMom_TruMom_p","MCSMom vs Tru Mom electrons", 50, 0, 1000, 50, 0, 1000);

    // Same as above but with good Efficiency * Purity
    fMCSMomEP_TruMom_e = tfs->make<TH2F>("MCSMomEP_TruMom_e","MCSMom vs Tru Mom electrons", 50, 0, 100, 50, 0, 1000);

    // Initialize the variables used to calculate Efficiency * Purity (aka EP) for matching to truth
    for(unsigned short pdgIndex = 0; pdgIndex < 6; ++pdgIndex) {
      EPTSums[pdgIndex] = 0;
      EPSums[pdgIndex] = 0;
      EPCnts[pdgIndex] = 0;
    }
    fEventsProcessed = 0;
    
    nTruPrimaryVtxOK = 0;
    nTruPrimaryVtxReco = 0;
    
  }
  
  bool TrajClusterAlg::SortByMultiplet(TCHit const& a, TCHit const& b)
  {
    // compare the wire IDs first:
    int cmp_res = a.WireID.cmp(b.WireID);
    if (cmp_res != 0) return cmp_res < 0; // order is decided, unless equal
    // decide by start time
    if (a.StartTick != b.StartTick) return a.StartTick < b.StartTick;
    // if still undecided, resolve by local index
    return a.LocalIndex < b.LocalIndex; // if still unresolved, it's a bug!
  } // SortByMultiplet

  //------------------------------------------------------------------------------
  void TrajClusterAlg::reconfigure(fhicl::ParameterSet const& pset)
  {

    bool badinput = false;
    fHitFinderModuleLabel = pset.get<art::InputTag>("HitFinderModuleLabel");
    fMakeNewHits          = pset.get< bool >("MakeNewHits", true);
    fMode                 = pset.get< short >("Mode", 1);
    fHitErrFac            = pset.get< float >("HitErrFac", 0.4);
    fMinAmp               = pset.get< float >("MinAmp", 0.1);
    tjs.AngleRanges       = pset.get< std::vector<float>>("AngleRanges");
    fNPtsAve              = pset.get< short >("NPtsAve", 20);
    fMinPtsFit            = pset.get< std::vector<unsigned short >>("MinPtsFit");
    fMinPts               = pset.get< std::vector<unsigned short >>("MinPts");
    if(pset.has_key("MaxAngleRange")) {
      fMaxAngleCode         = pset.get< std::vector<unsigned short>>("MaxAngleRange");
      std::cout<<"MaxAngleRange has been re-named MaxAngleCode: 0 = small angle (first, second range), 1 = large angle (2nd to last range), 2 = Very large angle (last range)\n";
    } else {
      fMaxAngleCode         = pset.get< std::vector<unsigned short>>("MaxAngleCode");
    }
    fMinMCSMom             = pset.get< std::vector<short >>("MinMCSMom");

    fMaxChi               = pset.get< float >("MaxChi", 10);
    fChargeCuts           = pset.get< std::vector<float >>("ChargeCuts", {3, 0.15, 0.25});
    fMultHitSep           = pset.get< float >("MultHitSep", 2.5);
    fKinkCuts             = pset.get< std::vector<float >>("KinkCuts", {0.4, 1.5, 4});
    fQualityCuts          = pset.get< std::vector<float >>("QualityCuts", {0.8, 3});
    fMaxWireSkipNoSignal  = pset.get< float >("MaxWireSkipNoSignal", 1);
    fMaxWireSkipWithSignal= pset.get< float >("MaxWireSkipWithSignal", 100);
    fProjectionErrFactor  = pset.get< float >("ProjectionErrFactor", 2);
    // Step size for very large angle trajectories
    if(pset.has_key("VLAStepSize")) {
      fVLAStepSize        = pset.get< float >("VLAStepSize", 1.5);
    } else {
      fVLAStepSize = 1.5;
    }
    fJTMaxHitSep2         = pset.get< float >("JTMaxHitSep", 2);
    
    std::vector<std::string> skipAlgsVec = pset.get< std::vector<std::string>  >("SkipAlgs");
    
    std::vector<std::string> specialAlgsVec;
    if(pset.has_key("SpecialAlgs")) specialAlgsVec = pset.get< std::vector<std::string>  >("SpecialAlgs");
    
    tjs.DeltaRayTag       = pset.get< std::vector<short>>("DeltaRayTag", {-1, -1, -1});
    tjs.MuonTag           = pset.get< std::vector<short>>("MuonTag", {-1, -1, -1, - 1});
    tjs.ShowerTag         = pset.get< std::vector<float>>("ShowerTag", {-1, -1, -1, -1, -1, -1});
    fChkStopCuts          = pset.get< std::vector<float>>("ChkStopCuts", {-1, -1, -1});
    fMaxTrajSep           = pset.get< float >("MaxTrajSep", 4);
    
    fStudyMode            = pset.get< bool  >("StudyMode", false);
    fMatchTruth           = pset.get< std::vector<float> >("MatchTruth", {-1, -1, -1, -1});
    tjs.Vertex2DCuts      = pset.get< std::vector<float >>("Vertex2DCuts", {-1, -1, -1, -1, -1, -1, -1});
    if(pset.has_key("VertexScoreWeights")) tjs.VertexScoreWeights = pset.get< std::vector<float> >("VertexScoreWeights");
    tjs.Vertex3DChiCut    = pset.get< float >("Vertex3DChiCut", -1);
    fMaxVertexTrajSep     = pset.get< std::vector<float>>("MaxVertexTrajSep");
    tjs.Match3DCuts       = pset.get< std::vector<float >>("Match3DCuts", {-1, -1, -1, -1});
    
    debug.Cryostat = 0;
    debug.TPC = 0;
    if(pset.has_key("DebugCryostat")) debug.Cryostat = pset.get< int >("DebugCryostat");
    if(pset.has_key("DebugTPC"))      debug.TPC = pset.get< int >("DebugTPC");
    debug.Plane           = pset.get< int >("DebugPlane", -1);
    debug.Wire            = pset.get< int >("DebugWire", -1);
    debug.Tick            = pset.get< int >("DebugTick", -1);
    debug.WorkID          = pset.get< int >("DebugWorkID", 0);

    // convert the max traj separation into a separation^2
    fMaxTrajSep *= fMaxTrajSep;
    if(fJTMaxHitSep2 > 0) fJTMaxHitSep2 *= fJTMaxHitSep2;
    
    // in the following section we ensure that the fcl vectors are appropriately sized so that later references are valid
    if(fMinPtsFit.size() != fMinPts.size()) badinput = true;
    if(fMaxVertexTrajSep.size() != fMinPts.size()) badinput = true;
    if(fMaxAngleCode.size() != fMinPts.size()) badinput = true;
    if(fMinMCSMom.size() != fMinPts.size()) badinput = true;
    if(badinput) throw art::Exception(art::errors::Configuration)<< "Bad input from fcl file. Vector lengths for MinPtsFit, MaxVertexTrajSep, MaxAngleRange and MinMCSMom should be defined for each reconstruction pass";
    
    if(tjs.Vertex2DCuts.size() < 7) throw art::Exception(art::errors::Configuration)<<"Vertex2DCuts must be size 7\n 0 = Max length definition for short TJs\n 1 = Max vtx-TJ sep short TJs\n 2 = Max vtx-TJ sep long TJs\n 3 = Max position pull for >2 TJs\n 4 = Max vtx position error\n 5 = Min MCSMom for one of two TJs\n 6 = Min fraction of wires hit btw vtx and Tjs\n 7 = Min Score\n 8 = ID of a 2D vertex to print";
    // resize for a debug element of the vector
    if(tjs.Vertex2DCuts.size() < 9) {
      tjs.Vertex2DCuts.resize(9);
      // Set a default minimum Score
      tjs.Vertex2DCuts[7] = 3;
    }
    if(fKinkCuts.size() != 3) throw art::Exception(art::errors::Configuration)<<"KinkCuts must be size 2\n 0 = Hard kink angle cut\n 1 = Kink angle significance\n 2 = nPts fit";
    if(fChargeCuts.size() != 3) throw art::Exception(art::errors::Configuration)<<"ChargeCuts must be size 3\n 0 = Charge pull cut\n 1 = Min allowed fractional chg RMS\n 2 = Max allowed fractional chg RMS";
    
    if(tjs.MuonTag.size() != 4) throw art::Exception(art::errors::Configuration)<<"MuonTag must be size 4\n 0 = minPtsFit\n 1 = minMCSMom\n 2= maxWireSkipNoSignal\n 3 = min delta ray length for tagging";
    if(tjs.DeltaRayTag.size() != 3) throw art::Exception(art::errors::Configuration)<<"DeltaRayTag must be size 3\n 0 = Max endpoint sep\n 1 = min MCSMom\n 2 = max MCSMom";
    if(fChkStopCuts.size() != 3) throw art::Exception(art::errors::Configuration)<<"ChkStopCuts must be size 3\n 0 = Min Charge ratio\n 1 = Charge slope pull cut\n 2 = Charge fit chisq cut";
    if(tjs.ShowerTag.size() < 13) {
      std::cout<< "ShowerTag must be size 13\n 0 = Mode\n 1 = max MCSMom\n 2 = max separation (WSE units)\n 3 = Max angle diff\n 4 = Factor * rms width\n 5 = Min half width\n 6 = min total Tps\n 7 = Min Tjs\n 8 = max parent FOM\n 9 = max direction FOM 10 = max AspectRatio\n 11 = min Score to preserve a vertex\n 12 = Debug showers in CTP\n";
      std::cout<<" Fixing this problem...";
      tjs.ShowerTag.resize(13);
      // set the min score to 0
      tjs.ShowerTag[11] = 0;
      // turn off printing
      tjs.ShowerTag[12] = -1;
    }
    if(tjs.Match3DCuts.size() != 4) throw art::Exception(art::errors::Configuration)<< "Match3DCuts must be size 3\n 0 = dx(cm) match\n 1 = MinMCSMom\n 2 = Min length for 2-view match\n 3 = 2-view match require dead region in 3rd view?";
    
    // check the angle ranges and convert from degrees to radians
    if(tjs.AngleRanges.back() < 90) {
      mf::LogVerbatim("TC")<<"Last element of AngleRange != 90 degrees. Fixing it\n";
      tjs.AngleRanges.back() = 90;
    }
    
    fExpectNarrowHits = (fMode == 4);
    
    // decide whether debug information should be printed
    bool validCTP = debug.Cryostat >= 0 && debug.TPC >= 0 && debug.Plane >= 0;
    if(validCTP) debug.CTP = EncodeCTP((unsigned int)debug.Cryostat, (unsigned int)debug.TPC, (unsigned int)debug.Plane);
    fDebugMode = validCTP || debug.WorkID < 0;
    if(fDebugMode) {
      std::cout<<"**************** Debug mode: debug.CTP "<<debug.CTP<<" ****************\n";
      std::cout<<"Cryostat "<<debug.Cryostat<<" TPC "<<debug.TPC<<" Plane "<<debug.Plane<<"\n";
      std::cout<<"Pass MinPts  MinPtsFit Max Angle\n";
      for(unsigned short pass = 0; pass < fMinPts.size(); ++pass) {
        unsigned short ir = fMaxAngleCode[pass];
        if(ir > tjs.AngleRanges.size() - 1) ir = tjs.AngleRanges.size() - 1;
        std::cout<<std::setw(3)<<pass;
        std::cout<<std::setw(7)<<fMinPts[pass];
        std::cout<<std::setw(7)<<fMinPtsFit[pass];
        std::cout<<std::setw(12)<<(int)tjs.AngleRanges[ir]<<"\n";
      }
      std::cout<<"Max angle ranges\n range  degrees  radians\n";
      for(unsigned short ir = 0; ir < tjs.AngleRanges.size(); ++ir) {
        std::cout<<std::setw(3)<<ir;
        std::cout<<std::setw(8)<<(int)tjs.AngleRanges[ir];
        std::cout<<std::setw(8)<<std::setprecision(3)<<tjs.AngleRanges[ir] * M_PI / 180;
        std::cout<<"\n";
      } // ir
    }
    
    for(auto& range : tjs.AngleRanges) {
      if(range < 0 || range > 90) throw art::Exception(art::errors::Configuration)<< "Invalid angle range "<<range<<" Must be 0 - 90 degrees";
      range *= M_PI / 180;
    } // range
    
    // Ensure that the size of the AlgBitNames vector is consistent with the AlgBit typedef enum
    if(kAlgBitSize != AlgBitNames.size()) throw art::Exception(art::errors::Configuration)<<"kAlgBitSize "<<kAlgBitSize<<" != AlgBitNames size "<<AlgBitNames.size();
    fAlgModCount.resize(kAlgBitSize);

    if(kFlagBitSize != StopFlagNames.size()) throw art::Exception(art::errors::Configuration)<<"kFlagBitSize "<<kFlagBitSize<<" != StopFlagNames size "<<StopFlagNames.size();
    
    if(kFlagBitSize > 64) throw art::Exception(art::errors::Configuration)<<"Increase the size of UseAlgs to at least "<<kFlagBitSize;
    
    bool gotit, printHelp = false;
    for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) tjs.UseAlg[ib] = true;
    
    // turn off the special algs
    // A lightly tested algorithm that should only be turned on for testing
    tjs.UseAlg[kStopAtTj] = false;
    // Do an exhaustive (and slow) check of the hit -> trajectory associations
    tjs.UseAlg[kChkInTraj] = false;
    
    for(auto strng : skipAlgsVec) {
      gotit = false;
      if(strng == "All") {
        // turn everything off
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) tjs.UseAlg[ib] = false;
        gotit = true;
        break;
      } // All off
      for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
        if(strng == AlgBitNames[ib]) {
          tjs.UseAlg[ib] = false;
          gotit = true;
          break;
        }
      } // ib
      if(!gotit) {
        std::cout<<"******* Unknown SkipAlgs input string '"<<strng<<"'\n";
        printHelp = true;
      }
    } // strng
    if(printHelp) {
      std::cout<<"Valid AlgNames:";
      for(auto strng : AlgBitNames) std::cout<<" "<<strng;
      std::cout<<"\n";
      std::cout<<"Or specify All to turn all algs off\n";
      throw art::Exception(art::errors::Configuration)<< "Invalid SkipAlgs specification";
    }
    // overwrite any settings above with special algs
    for(auto strng : specialAlgsVec) {
      gotit = false;
      for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
        if(strng == AlgBitNames[ib]) {
          tjs.UseAlg[ib] = true;
          gotit = true;
          break;
        }
      } // ib
      if(!gotit) {
        std::cout<<"******* Unknown SpecialAlgs input string '"<<strng<<"'\n";
        printHelp = true;
      }
    } // strng
    if(printHelp) {
      std::cout<<"Valid AlgNames:";
      for(auto strng : AlgBitNames) std::cout<<" "<<strng;
      std::cout<<"\n";
      std::cout<<"Or specify All to turn all algs off\n";
      throw art::Exception(art::errors::Configuration)<< "Invalid SkipAlgs specification";
    }
    
    if(fDebugMode) {
      std::cout<<"Using algs:";
      for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
        if(ib % 10 == 0) std::cout<<"\n";
        if(tjs.UseAlg[ib] && ib != kKilled) std::cout<<" "<<AlgBitNames[ib];
      }
      std::cout<<"\n";
      std::cout<<"Skipping algs:";
      for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(!tjs.UseAlg[ib] && ib != kKilled) std::cout<<" "<<AlgBitNames[ib];
      std::cout<<"\n";
      if(fExpectNarrowHits) std::cout<<"Configured to expect narrow hits produced by gaushit with LongMaxHits set large.\n";
      if(tjs.IgnoreNegChiHits) std::cout<<"Configured to ignore hits with GoodnessOfFit < 0.\n";
    }
   
  } // reconfigure
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ClearResults() {
    tjs.vtx3.clear();
    tjs.vtx.clear();
    tjs.tcl.clear();
    tjs.inClus.clear();
    tjs.matchVec.clear();
    tjs.matchVecPFPList.clear();
    tjs.WireHitRange.clear();
    tjs.fHits.clear();
    tjs.allTraj.clear();
    tjs.cots.clear();
    tjs.showers.clear();
    tjs.MCPartList.clear();
  } // ClearResults()

  ////////////////////////////////////////////////
  void TrajClusterAlg::RunTrajClusterAlg(art::Event & evt)
  {

    if(fMode == 0) return;
    
    // Get the hits
    art::ValidHandle< std::vector<recob::Hit>> hitVecHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitFinderModuleLabel);

    ++fEventsProcessed;
    if(hitVecHandle->size() < 3) return;
    
    // a gratuitous clearing of everything before we start
    ClearResults();
 
    tjs.detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    tjs.geom = lar::providerFrom<geo::Geometry>();
    
    tjs.fHits.reserve(hitVecHandle->size());

    // transfer the hits into the local vector so we can modify them
    float minAmp = fMinAmp;
    if(fMinAmp < 0) minAmp = 0;
    for(unsigned int iht = 0; iht < hitVecHandle->size(); ++iht) {
      art::Ptr<recob::Hit> hit = art::Ptr<recob::Hit>(hitVecHandle, iht);
      // Look for a hit with negative amplitude
      if(hit->PeakAmplitude() < minAmp) continue;
      TCHit localHit;
      localHit.StartTick = hit->StartTick();
      localHit.EndTick = hit->EndTick();
      localHit.PeakTime = hit->PeakTime();
      localHit.SigmaPeakTime = hit->SigmaPeakTime();
      localHit.PeakAmplitude = hit->PeakAmplitude();
      localHit.SigmaPeakAmp = hit->SigmaPeakAmplitude();
      localHit.Integral = hit->Integral();
      localHit.SigmaIntegral = hit->SigmaIntegral();
      localHit.RMS = hit->RMS();
      localHit.GoodnessOfFit = hit->GoodnessOfFit();
      localHit.NDOF = hit->DegreesOfFreedom();
      localHit.Multiplicity = hit->Multiplicity();
      localHit.LocalIndex = hit->LocalIndex();
      localHit.WireID = hit->WireID();
      tjs.fHits.push_back(localHit);
    } // iht    

    // sort it as needed;
    // that is, sorted by wire ID number,
    // then by start of the region of interest in time, then by the multiplet
    std::sort(tjs.fHits.begin(), tjs.fHits.end(), &SortByMultiplet);
    
    // Match these hits to MC tracks
    MatchTrueHits();

    // check for debugging mode triggered by Plane, Wire, Tick
    debug.Hit = UINT_MAX;
    if(fDebugMode) {
      std::cout<<"Look for debug hit "<<debug.Plane<<":"<<debug.Wire<<":"<<debug.Tick;
      for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
        if((int)tjs.fHits[iht].WireID.Plane != debug.Plane) continue;
        if((int)tjs.fHits[iht].WireID.Wire != debug.Wire) continue;
        if(tjs.fHits[iht].PeakTime < debug.Tick - 5) continue;
        if(tjs.fHits[iht].PeakTime > debug.Tick + 5) continue;
        debug.Hit = iht;
        std::cout<<" iht "<<iht<<" "<<debug.Cryostat<<":"<<debug.TPC<<":"<<debug.Plane<<":"<<PrintHit(tjs.fHits[iht]);
        std::cout<<" Amp "<<(int)tjs.fHits[iht].PeakAmplitude;
        std::cout<<" RMS "<<std::fixed<<std::setprecision(1)<<tjs.fHits[iht].RMS;
        std::cout<<" Chisq "<<std::fixed<<std::setprecision(1)<<tjs.fHits[iht].GoodnessOfFit;
        std::cout<<" Mult "<<tjs.fHits[iht].Multiplicity;
        std::cout<<"\n";
        break;
      } // iht
      if(debug.Hit == UINT_MAX) std::cout<<" not found\n";
    } // debugging mode

    
    fRun = evt.run();
    fSubRun  = evt.subRun();
    fEvent = evt.event();
    fWorkID = 0;
    
    // Set true if a truly bad situation occurs
    fQuitAlg = false;
    fIsRealData = evt.isRealData();
    didPrt = false;
    
    if(fMode > 0) {
      // Step from upstream (small wire number) to downstream (large wire number)
      fStepDir = 1;
    } else {
      // or step in the other direction
      fStepDir = -1;
    }
    InitializeAllTraj();
    for (const geo::TPCID& tpcid: tjs.geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
      FillWireHitRange(tpcid);
      if(fQuitAlg) return;
      for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
        // special mode for only reconstructing the collection plane
        if(fMode == 2 && fPlane != TPC.Nplanes() - 1) continue;
        // no hits on this plane?
        if(tjs.FirstWire[fPlane] > tjs.LastWire[fPlane]) continue;
        // Set the CTP code to ensure objects are compared within the same plane
        fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
        fCstat = tpcid.Cryostat;
        fTpc = tpcid.TPC;
        // save some processing time if in debug mode. This may cause confusion when debugging
        // subtle failures
//        if(debug.CTP != UINT_MAX && debug.CTP != fCTP) continue;
        // reconstruct all trajectories in the current plane
        ReconstructAllTraj();
        if(fQuitAlg) {
          mf::LogVerbatim("TC")<<"Found fQuitAlg after ReconstructAllTraj";
          ClearResults();
          return;
        }
      } // fPlane
      // No sense taking muon direction if delta ray tagging is disabled
      if(tjs.DeltaRayTag[0] >= 0) TagMuonDirections(tjs, debug.WorkID);
      if(tjs.Vertex3DChiCut > 0) Find3DVertices(tjs, debug, tpcid);
      for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
        fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
        ChkVtxAssociations(tjs, fCTP);
      }
      if(tjs.ShowerTag[0] > 0) {
        // find showers after 3D vertex finding - which may split trajectories
        for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
          // no hits on this plane?
          if(tjs.FirstWire[fPlane] > tjs.LastWire[fPlane]) continue;
          // Set the CTP code to ensure objects are compared within the same plane
          fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
          FindShowers(tjs, fCTP);
        }
      } // make showers
      // Match3D should be the last thing called for this tpcid
      Match3D(tpcid);
    } // tpcid

    MatchTruth();
    
    if(fStudyMode) {
      // output MC-reco stuff to optimize the vertex weights
      for (const geo::TPCID& tpcid: tjs.geom->IterateTPCIDs()) {
        geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
        for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
          fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
          ChkVtxAssociations(tjs, fCTP);
        }
      } // tpcid
    }
    
    FillPFPInfo();
    // convert the cots vector into recob::Shower
    MakeShowers(tjs, fCaloAlg);
    KillPoorVertices(tjs);
    // Convert trajectories in allTraj into clusters
    MakeAllTrajClusters();
    if(fQuitAlg) {
      mf::LogVerbatim("TC")<<"RunTrajCluster failed in MakeAllTrajClusters";
      ClearResults();
      return;
    }
    if(!CheckHitClusterAssociations(tjs)) {
      ClearResults();
      mf::LogVerbatim("TC")<<"RunTrajCluster failed in CheckHitClusterAssociations";
      return;
    }
/*
    if(fDebugMode) {
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        if(tjs.allTraj[itj].WorkID == TJPrt) {
          PrintAllTraj("DBG", tjs, debug, itj, USHRT_MAX);
          break;
        }
      } // itj
    } // TJPrt > 0
*/
    // print trajectory summary report?
    if(tjs.ShowerTag[0] >= 0) debug.Plane = tjs.ShowerTag[11];
    if(fDebugMode) {
      mf::LogVerbatim("TC")<<"Done in RunTrajClusterAlg";
      PrintAllTraj("RTC", tjs, debug, USHRT_MAX, 0);
    }

    unsigned short ntj = 0;
    unsigned short nsh = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      ++ntj;
      if(tj.AlgMod[kShowerTj]) ++nsh;
    } // tj
    if(fDebugMode) std::cout<<"RTC done ntj "<<ntj<<" nsh "<<nsh<<" events processed "<<fEventsProcessed<<"\n";

    if(fStudyMode) {
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.AlgMod[kKilled]) continue;
        // reco MCSMom vs reco range
        float len = TrajLength(tj);
        if(len > 99) len = 99;
        // ignore really short Tjs
        if(len < 2) continue;
        if(tj.MCPartListIndex == USHRT_MAX) continue;
        auto& mcp = tjs.MCPartList[tj.MCPartListIndex];
        int truKE = 1000 * (mcp->E() - mcp->Mass());
        int pdg = abs(mcp->PdgCode());
        double mass = 0.511;
        if(pdg == 13) mass = 105.7;
        if(pdg == 211) mass = 139.6;
        if(pdg == 2212) mass = 938.3;
        double tPlusM = truKE + mass;
        double truMom = sqrt(tPlusM * tPlusM - mass * mass);
        if(pdg == 11) fMCSMom_TruMom_e->Fill(truMom, tj.MCSMom);
        if(pdg == 13) fMCSMom_TruMom_mu->Fill(truMom, tj.MCSMom);
        if(pdg == 211) fMCSMom_TruMom_pi->Fill(truMom, tj.MCSMom);
        if(pdg == 2212) fMCSMom_TruMom_p->Fill(truMom, tj.MCSMom);
        // See if a parameterization of expected MCSMom(Length) 
        if(pdg == 11 && tj.EffPur > 0.7) fMCSMomEP_TruMom_e->Fill(truMom, tj.MCSMom);
        // check charge rms for electrons vs others
        if(tj.EffPur > 0.7) {
          float rms = 0;
          float cnt = 0;
          for(unsigned short ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1] - 1; ++ipt) {
            TrajPoint& tpm = tj.Pts[ipt - 1];
            if(tpm.Chg == 0) continue;
            TrajPoint& tp = tj.Pts[ipt];
            if(tp.Chg == 0) continue;
            TrajPoint& tpp = tj.Pts[ipt + 1];
            if(tpp.Chg == 0) continue;
            if(std::abs(tpp.Pos[0] - tpm.Pos[0]) > 4) continue;
            ++cnt;
            float expectChg = 0.5 * (tpm.Chg + tpp.Chg);
            float chgRat = std::abs(tp.Chg - expectChg) / expectChg;
            rms += chgRat;
          } // tp
          if(cnt == 0) continue;
          rms /= cnt;
          mf::LogVerbatim("TC")<<"NTP "<<pdg<<" "<<(int)cnt<<" "<<truKE<<" "<<std::fixed<<std::setprecision(3)<<rms;
        } // long matched trajectory
      } // itj
      
      for(auto& vx2 : tjs.vtx) {
        if(vx2.Stat[kVtxKilled]) continue;
        fVx2_Score->Fill((float)vx2.Score);
      } // vx2
      for(auto& vx3 : tjs.vtx3) {
        if(vx3.ID == 0) continue;
        float score = 0;
        for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
          if(vx3.Vtx2ID[ipl] == 0) continue;
          unsigned short iv2 = vx3.Vtx2ID[ipl] - 1;
          score += tjs.vtx[iv2].Score;
        } // ipl
        fVx3_Score->Fill(score);
      } //  vx3
    } // studymode
    
    if(fMatchTruth[0] >= 0) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Event "<<evt.event();
      float sum = 0;
      float sumt = 0;
      for(unsigned short pdgIndex = 0; pdgIndex < EPSums.size(); ++pdgIndex) {
        if(EPSums[pdgIndex] == 0) continue;
        if(pdgIndex == 0) myprt<<" Electron";
        if(pdgIndex == 1) myprt<<" Muon";
        if(pdgIndex == 2) myprt<<" Pion";
        if(pdgIndex == 3) myprt<<" Kaon";
        if(pdgIndex == 4) myprt<<" Proton";
        float ave = EPTSums[pdgIndex] / (float)EPSums[pdgIndex];
        myprt<<" "<<std::fixed<<std::setprecision(2)<<ave;
        myprt<<" "<<EPCnts[pdgIndex];
        if(pdgIndex > 0) {
          sum  += EPSums[pdgIndex];
          sumt += EPTSums[pdgIndex];
        }
      } // pdgIndex
      if(sum > 0) myprt<<" MuPiKP "<<std::fixed<<std::setprecision(2)<<sumt / sum;
    } // fMatchTruth[0] >= 0

    // convert vertex time from WSE to ticks
    for(auto& avtx : tjs.vtx) avtx.Pos[1] /= tjs.UnitsPerTick;
    
    if(fDebugMode) mf::LogVerbatim("TC")<<"RunTrajCluster success run "<<fRun<<" event "<<fEvent<<" allTraj size "<<tjs.allTraj.size()<<" events processed "<<fEventsProcessed;
    
  } // RunTrajClusterAlg

  ////////////////////////////////////////////////
  void TrajClusterAlg::InitializeAllTraj()
  {
    tjs.allTraj.clear();
    tjs.vtx.clear();
    tjs.vtx3.clear();
  } // InitializeAllTraj
  
  ////////////////////////////////////////////////
  bool TrajClusterAlg::MergeAndStore(unsigned short itj1, unsigned short itj2, bool doPrt)
  {
    // Merge the two trajectories in allTraj and store them. Returns true if it was successfull.
    // Merging is done between the end (end = 1) of tj1 and the beginning (end = 0) of tj2. This function preserves the
    // AlgMod state of itj1
    
    if(itj1 > tjs.allTraj.size() - 1) return false;
    if(itj2 > tjs.allTraj.size() - 1) return false;
    if(tjs.allTraj[itj1].AlgMod[kKilled] || tjs.allTraj[itj2].AlgMod[kKilled]) return false;
    
    // Merging shower Tjs requires merging the showers as well.
    if(tjs.allTraj[itj1].AlgMod[kShowerTj] || tjs.allTraj[itj2].AlgMod[kShowerTj]) return MergeShowersAndStore(tjs, itj1, itj2, prt);
    
    // make a copy so they can be trimmed as needed
    Trajectory tj1 = tjs.allTraj[itj1];
    Trajectory tj2 = tjs.allTraj[itj2];
    
    // ensure that the tj1 - tj2 order is correct
    if(std::abs(tj1.Pts[tj1.EndPt[1]].Pos[0] - tj2.Pts[tj2.EndPt[1]].Pos[0]) < 
       std::abs(tj1.Pts[tj1.EndPt[1]].Pos[0] - tj2.Pts[tj2.EndPt[0]].Pos[0])) {
      return false;
    }

    if(doPrt) {
      mf::LogVerbatim("TC")<<"MergeAndStore: tj1.ID "<<tj1.ID<<" tj2.ID "<<tj2.ID;
    }
    
    // ensure that these are in the same step order
    if(tj1.StepDir != tj2.StepDir) ReverseTraj(tjs, tj2);
    
    if(tj1.StopFlag[1][kBragg]) std::cout<<"MergeAndStore: You are merging the end of a trajectory "<<tj1.ID<<" with a Bragg peak. Is this wise?\n";
    
    // assume that everything will succeed
    fQuitAlg = false;
    
    // remove any points at the end of tj1 that don't have used hits
    tj1.Pts.resize(tj1.EndPt[1] + 1);
    
    // determine if they overlap by finding the point on tj2 that is closest
    // to the end point of tj1.
    TrajPoint& endtj1TP = tj1.Pts[tj1.EndPt[1]];
    // Set minSep large so that dead wire regions are accounted for
    float minSep = 1000;
    unsigned short tj2ClosePt = 0;
    // Note that TrajPointTrajDOCA only considers TPs that have charge
    TrajPointTrajDOCA(tjs, endtj1TP, tj2, tj2ClosePt, minSep);
    if(doPrt) mf::LogVerbatim("TC")<<" Merge point tj1 "<<PrintPos(tjs, endtj1TP)<<" tj2ClosePt "<<tj2ClosePt<<" Pos "<<PrintPos(tjs, tj2.Pts[tj2ClosePt]);
    // check for full overlap
    if(tj2ClosePt > tj2.EndPt[1]) return false;
    // check for the following possibilities, where - indicates a TP with charge
    // tj1:  --------------
    // tj2:                  -------------
    // Another possibility with overlap
    // tj1:  --------------
    // tj2:               ---------------
    
    // The approach is to append tj2 to tj1, store tj1 as a new trajectory,
    // and re-assign all hits to the new trajectory
    
    // First ensure that any hit will appear only once in the merged trajectory in the overlap region
    // whether it is used or unused. The point on tj2 where the merge will begin, tj2ClosePt, will be
    // increased until this condition is met.
    // Make a temporary vector of tj1 hits in the end points for simpler searching
    std::vector<unsigned int> tj1Hits;
    for(unsigned short ii = 0; ii < tj1.Pts.size(); ++ii) {
      // only go back a few points in tj1
      if(ii > 10) break;
      unsigned short ipt = tj1.Pts.size() - 1 - ii;
      tj1Hits.insert(tj1Hits.end(), tj1.Pts[ipt].Hits.begin(), tj1.Pts[ipt].Hits.end());
      if(ipt == 0) break;
    } // ii
    
    bool bumpedPt = true;
    while(bumpedPt) {
      bumpedPt = false;
      for(unsigned short ii = 0; ii < tj2.Pts[tj2ClosePt].Hits.size(); ++ii) {
        unsigned int iht = tj2.Pts[tj2ClosePt].Hits[ii];
        if(std::find(tj1Hits.begin(), tj1Hits.end(), iht) != tj1Hits.end()) bumpedPt = true;
      } // ii
      if(bumpedPt && tj2ClosePt < tj2.EndPt[1]) {
        ++tj2ClosePt;
      } else {
        break;
      }
    } // bumpedPt
    if(doPrt) mf::LogVerbatim("TC")<<" revised tj2ClosePt "<<tj2ClosePt;
    // append tj2 hits to tj1

    tj1.Pts.insert(tj1.Pts.end(), tj2.Pts.begin() + tj2ClosePt, tj2.Pts.end());
    // re-define the end points
    SetEndPoints(tjs, tj1);
    
    // A more exhaustive check that hits only appear once
    if(HasDuplicateHits(tjs, tj1, doPrt)) {
      if(doPrt) {
        mf::LogVerbatim("TC")<<"MergeAndStore found duplicate hits. Coding error";
        PrintTrajectory("MAS", tjs, tj1, USHRT_MAX);
      }
      return false;
    }
    // kill the original trajectories
    MakeTrajectoryObsolete(tjs, itj1);
    MakeTrajectoryObsolete(tjs, itj2);
    // Do this so that StoreTraj keeps the correct WorkID (of itj1)
    tj1.ID = tj1.WorkID;
    SetPDGCode(tjs, tj1);
    ChkStop(tj1);
    return StoreTraj(tjs, tj1);
    if(tjs.UseAlg[kChkInTraj]) {
      fQuitAlg = !InTrajOK(tjs, "MAS");
      if(fQuitAlg) mf::LogVerbatim("TC")<<"InTrajOK failed in MergeAndStore";
      return false;
    }
    return true;
    
  } // MergeAndStore
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ReconstructAllTraj()
  {
    // Reconstruct clusters in fPlane and put them in allTraj
    
    unsigned int ii, iwire, jwire, iht, jht;
    
    if(fPlane > tjs.FirstWire.size() - 1) {
      mf::LogWarning("TC")<<"ReconstructAllTraj called with invalid fPlane "<<fPlane;
      fQuitAlg = true;
      return;
    }
    
    unsigned int nwires = tjs.LastWire[fPlane] - tjs.FirstWire[fPlane] - 1;
    std::vector<unsigned int> iHitsInMultiplet, jHitsInMultiplet;
    unsigned short ihtIndex, jhtIndex;
    
    // turn of trajectory printing
    TJPrt = 0;
    didPrt = false;
    
    // Make several passes through the hits with user-specified cuts for each
    // pass. In general these are to not reconstruct large angle trajectories on
    // the first pass
    float maxHitsRMS = 4 * fAveHitRMS[fPlane];
    for(unsigned short pass = 0; pass < fMinPtsFit.size(); ++pass) {
      for(ii = 0; ii < nwires; ++ii) {
        // decide which way to step given the sign of fStepDir
        if(fStepDir > 0) {
          // step DS
          iwire = tjs.FirstWire[fPlane] + ii;
          jwire = iwire + 1;
        } else {
          // step US
          iwire = tjs.LastWire[fPlane] - ii - 1;
          jwire = iwire - 1;
        }
        // skip bad wires or no hits on the wire
        if(tjs.WireHitRange[fPlane][iwire].first < 0) continue;
        if(tjs.WireHitRange[fPlane][jwire].first < 0) continue;
        unsigned int ifirsthit = (unsigned int)tjs.WireHitRange[fPlane][iwire].first;
        unsigned int ilasthit = (unsigned int)tjs.WireHitRange[fPlane][iwire].second;
        unsigned int jfirsthit = (unsigned int)tjs.WireHitRange[fPlane][jwire].first;
        unsigned int jlasthit = (unsigned int)tjs.WireHitRange[fPlane][jwire].second;
        for(iht = ifirsthit; iht < ilasthit; ++iht) {
          prt = (iht == debug.Hit);
          if(prt)  {
            mf::LogVerbatim("TC")<<"+++++++ Pass "<<pass<<" Found debug hit "<<fPlane<<":"<<PrintHit(tjs.fHits[iht]);
            didPrt = true;
          }
          // Only consider hits that are available
          if(tjs.fHits[iht].InTraj != 0) continue;
          if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
          // We hope to make a trajectory point at the hit position of iht in WSE units
          // with a direction pointing to jht
          unsigned int fromWire = tjs.fHits[iht].WireID.Wire;
          float fromTick = tjs.fHits[iht].PeakTime;
          float iqtot = tjs.fHits[iht].Integral;
          float hitsRMSTick = tjs.fHits[iht].RMS;
          std::vector<unsigned int> iHitsInMultiplet;
          if(pass == 0) {
            // only use the hit on the first pass
            iHitsInMultiplet.resize(1);
            iHitsInMultiplet[0] = iht;
          } else {
            GetHitMultiplet(iht, iHitsInMultiplet, ihtIndex);
          }
          if(iHitsInMultiplet.size() > 1) {
            fromTick = HitsPosTick(tjs, iHitsInMultiplet, iqtot, kUnusedHits);
            hitsRMSTick = HitsRMSTick(tjs, iHitsInMultiplet, kUnusedHits);
          }
          bool fatIHit = (hitsRMSTick > maxHitsRMS);
          if(prt) mf::LogVerbatim("TC")<<" hit RMS "<<tjs.fHits[iht].RMS<<" BB Multiplicity "<<iHitsInMultiplet.size()<<" AveHitRMS["<<fPlane<<"] "<<fAveHitRMS[fPlane]<<" HitsRMSTick "<<hitsRMSTick<<" fatIHit "<<fatIHit;
          for(jht = jfirsthit; jht < jlasthit; ++jht) {
            // Only consider hits that are available
            if(tjs.fHits[iht].InTraj != 0) continue;
            if(tjs.fHits[jht].InTraj != 0) continue;
            if(tjs.IgnoreNegChiHits && tjs.fHits[jht].GoodnessOfFit < 0) continue;
            // clear out any leftover work inTraj's that weren't cleaned up properly
/*
            for(unsigned short oht = jfirsthit; oht < jlasthit; ++oht) {
              if(tjs.fHits[oht].InTraj < 0) {
                mf::LogVerbatim("TC")<<"Bad cleanup "<<PrintHit(tjs.fHits[oht])<<" "<<tjs.fHits[oht].InTraj<<" events processed "<<fEventsProcessed;
                std::cout<<"Bad cleanup "<<PrintHit(tjs.fHits[oht])<<" "<<tjs.fHits[oht].InTraj<<" events processed "<<fEventsProcessed<<" fWorkID "<<fWorkID<<"\n";
                tjs.fHits[oht].InTraj = 0;
              }
            }
*/
            unsigned int toWire = jwire;
            float toTick = tjs.fHits[jht].PeakTime;
            float jqtot = tjs.fHits[jht].Integral;
            if(jqtot < 1) continue;
            std::vector<unsigned int> jHitsInMultiplet;
            if(pass == 0) {
              // only use the hit on the first pass
              jHitsInMultiplet.resize(1);
              jHitsInMultiplet[0] = jht;
            } else {
              GetHitMultiplet(jht, jHitsInMultiplet, jhtIndex);
            }
            hitsRMSTick = HitsRMSTick(tjs, jHitsInMultiplet, kUnusedHits);
            bool fatJHit = (hitsRMSTick > maxHitsRMS);
            if(pass == 0) {
              // require both hits to be consistent
              if((fatIHit && !fatJHit) || (!fatIHit && fatJHit)) {
                if(prt) mf::LogVerbatim("TC")<<" jht fails "<<PrintHit(tjs.fHits[jht])<<" hit RMS "<<tjs.fHits[jht].RMS<<" mRMS "<<hitsRMSTick<<" fatJhit "<<fatJHit<<" max RMS "<<maxHitsRMS;
                continue;
              }
            } else {
              // pass > 0
              if(jHitsInMultiplet.size() > 1) toTick = HitsPosTick(tjs, jHitsInMultiplet, jqtot, kUnusedHits);
            }
            if(prt) {
              mf::LogVerbatim("TC")<<"+++++++ checking TrajHitsOK with jht "<<fPlane<<":"<<PrintHit(tjs.fHits[jht])<<" BB Multiplicity "<<jHitsInMultiplet.size()<<" HitsRMSTick "<<HitsRMSTick(tjs, jHitsInMultiplet, kUnusedHits)<<" fatJhit "<<fatJHit<<" setting toTick to "<<(int)toTick<<" TrajHitsOK "<<TrajHitsOK(tjs, iHitsInMultiplet, jHitsInMultiplet);
            }
            // Ensure that the hits StartTick and EndTick have the proper overlap
            if(!TrajHitsOK(tjs, iHitsInMultiplet, jHitsInMultiplet)) continue;
            // start a trajectory with direction from iht -> jht
            Trajectory work;
            if(!StartTraj(work, fromWire, fromTick, toWire, toTick, fCTP, pass)) continue;
            if(didPrt) TJPrt = work.WorkID;
            // check for a major failure
            if(fQuitAlg) return;
            if(work.Pts.empty()) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: StartTraj failed";
              ReleaseHits(tjs, work);
              continue;
            }
            // check for a large angle crawl
            if(work.Pts[0].AngleCode > fMaxAngleCode[work.Pass]) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: Wrong angle code "<<work.Pts[0].AngleCode<<" for this pass "<<work.Pass;
              ReleaseHits(tjs, work);
              continue;
            }
            work.Pts[0].DeltaRMS = fHitErrFac * tjs.UnitsPerTick * hitsRMSTick;
            // don't include the charge of iht since it will be low if this
            // is a starting/ending track
            work.AveChg = jqtot;
            // try to add close hits
            bool sigOK;
            AddHits(work, 0, sigOK);
            // check for a major failure
            if(fQuitAlg) return;
            if(!sigOK || work.Pts[0].Chg == 0) {
              if(prt) mf::LogVerbatim("TC")<<" No hits at initial trajectory point ";
              ReleaseHits(tjs, work);
              continue;
            }
            // move the TP position to the hit position but don't mess with the direction
            work.Pts[0].Pos = work.Pts[0].HitPos;
            // print the header and the first TP
            if(prt) PrintTrajectory("RAT", tjs, work, USHRT_MAX);
            // We can't update the trajectory yet because there is only one TP.
            work.EndPt[0] = 0;
            // now try stepping away
            StepCrawl(work);
            // check for a major failure
            if(fQuitAlg) return;
            if(prt) mf::LogVerbatim("TC")<<" After first StepCrawl. fGoodTraj "<<fGoodTraj<<" fTryWithNextPass "<<fTryWithNextPass;
            if(!fGoodTraj && fTryWithNextPass) {
              StepCrawl(work);
              if(!fGoodTraj || !fUpdateTrajOK) {
                if(prt) mf::LogVerbatim("TC")<<" xxxxxxx StepCrawl failed AGAIN. fTryWithNextPass "<<fTryWithNextPass;
                ReleaseHits(tjs, work);
                continue;
              } // Failed again
            }
            // Check the quality of the work trajectory
            CheckTraj(work);
            // check for a major failure
            if(fQuitAlg) return;
            if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: After CheckTraj EndPt "<<work.EndPt[0]<<"-"<<work.EndPt[1]<<" fGoodTraj "<<fGoodTraj<<" fTryWithNextPass "<<fTryWithNextPass;
            if(fTryWithNextPass) {
              // Most likely, the first part of the trajectory was good but the latter part
              // had too many unused hits. The work vector was
              // truncated and pass incremented, so give it another try
              work.AlgMod[kTryWithNextPass] = true;
              StepCrawl(work);
              // check for a major failure
              if(fQuitAlg) return;
              if(!fGoodTraj) {
                if(prt) mf::LogVerbatim("TC")<<" xxxxxxx StepCrawl failed AGAIN after CheckTraj";
                ReleaseHits(tjs, work);
                continue;
              } // Failed again
            } // fTryWithNextPass
            if(prt) mf::LogVerbatim("TC")<<"StepCrawl done: fGoodTraj "<<fGoodTraj<<" NumPtsWithCharge "<<NumPtsWithCharge(tjs, work, true)<<" cut "<<fMinPts[work.Pass];
            // decide if the trajectory is long enough for this pass
            if(!fGoodTraj || NumPtsWithCharge(tjs, work, true) < fMinPts[work.Pass]) {
              if(prt) mf::LogVerbatim("TC")<<" xxxxxxx Not enough points "<<NumPtsWithCharge(tjs, work, true)<<" minimum "<<fMinPts[work.Pass]<<" or !fGoodTraj";
              ReleaseHits(tjs, work);
              continue;
            }
            if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: calling StoreTraj with npts "<<work.EndPt[1];
            fQuitAlg = !StoreTraj(tjs, work);
            // check for a major failure
            if(fQuitAlg) return;
            if(tjs.UseAlg[kChkInTraj]) {
              fQuitAlg = !InTrajOK(tjs, "RAT");
              if(fQuitAlg) {
                mf::LogVerbatim("TC")<<"InTrajOK failed in ReconstructAllTraj";
                return;
              }
            }
            break;
          } // jht
        } // iht
      } // iwire
      // Try to merge trajectories before making vertices
      EndMerge();
      if(fQuitAlg) return;
      
      // Tag delta rays before making vertices
      TagDeltaRays(tjs, fCTP, debug.WorkID);

      // TY: Split high charge hits near the trajectory end
      ChkHiChgHits();

      Find2DVertices(tjs, debug, fCTP);
      FindVtxTjs();
      if(fQuitAlg) return;

    } // pass
    
    // Use unused hits in all trajectories
    UseUnusedHits();
    
    // make junk trajectories using nearby un-assigned hits
    if(fJTMaxHitSep2 > 0) {
      FindJunkTraj();
      if(fQuitAlg) return;
    }
    TagDeltaRays(tjs, fCTP, debug.WorkID);
    Find2DVertices(tjs, debug, fCTP);
    // check for a major failure
    if(fQuitAlg) return;

    // last attempt to attach Tjs to vertices
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) if(tjs.vtx[ivx].NTraj > 0) AttachAnyTrajToVertex(tjs, ivx, vtxPrt);
    
    // Check the Tj <-> vtx associations and define the vertex quality
    ChkVtxAssociations(tjs, fCTP);

    // TY: Improve hit assignments near vertex 
    VtxHitsSwap(tjs, debug, fCTP);

    // Refine vertices, trajectories and nearby hits
//    Refine2DVertices();
    
  } // ReconstructAllTraj

  //////////////////////////////////////////
  void TrajClusterAlg::UseUnusedHits()
  {
    if(tjs.allTraj.size() == 0) return;
    if(!tjs.UseAlg[kUseUnusedHits]) return;
    
    // max change in position allowed after adding all unused hits in a multiplet 
    float maxPosSep2 = 0.25;
    // Relax this cut for special tracking mode
    if(fMode == 2) maxPosSep2 = 1;
    
    std::vector<unsigned int> hitsInMultiplet;
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      // Find the max delta
      unsigned short firstPt = tj.EndPt[0];
      unsigned short lastPt = tj.EndPt[1];
      for(unsigned short ipt = firstPt; ipt <= lastPt; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(AngleRange(tjs, tp) == 0) continue;
        if(tp.Hits.empty()) continue;
        bool hitsAdded = false;
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          GetHitMultiplet(iht, hitsInMultiplet);
          if(hitsInMultiplet.size() > 1) {
            for(unsigned short jj = ii + 1; jj < tp.Hits.size(); ++jj) {
              if(!tp.UseHit[jj]) continue;
              if(std::find(hitsInMultiplet.begin(), hitsInMultiplet.end(), tp.Hits[jj]) != hitsInMultiplet.end()) {
                tp.UseHit[jj] = true;
                tjs.fHits[tp.Hits[jj]].InTraj = tj.ID;
                hitsAdded = true;
              }
            } // jj
          }
        } // ii
        if(hitsAdded) {
          // save the hit position
          std::array<float, 2> oldHitPos = tp.HitPos;
          DefineHitPos(tp);
          // keep it if 
          if(PosSep2(tj.Pts[ipt].HitPos, oldHitPos) < maxPosSep2) {
            tj.AlgMod[kUseUnusedHits] = true;
          } else {
            UnsetUsedHits(tjs, tj.Pts[ipt]);
          }
        }
      } // ipt
      if(tj.AlgMod[kUseUnusedHits]) SetEndPoints(tjs, tj);
    } // itj
    
  } // UseUnusedHits
  
  //////////////////////////////////////////
  void TrajClusterAlg::ReversePropagate(Trajectory& tj)
  {
    // Reverse the trajectory and step in the opposite direction. The
    // updated trajectory is returned if this process is successful
    
    if(!tjs.UseAlg[kRevProp]) return;
    
    if(tj.Pts.size() < 6) return;
    // only do this once
    if(tj.AlgMod[kRevProp]) return;
    
    // This code can't handle VLA trajectories
    if(tj.Pts[tj.EndPt[0]].AngleCode == 2) return;
    
    if(tj.EndPt[0] > 0) {
      tj.Pts.erase(tj.Pts.begin(), tj.Pts.begin() + tj.EndPt[0]);
      SetEndPoints(tjs, tj);
    }
    
    if(prt) mf::LogVerbatim("TC")<<"ReversePropagate: Prepping Tj "<<tj.ID<<" incoming StepDir "<<tj.StepDir;
    
    short stepDir = tj.StepDir;

    // find the wire on which EndPt resides
    unsigned int wire0 = std::nearbyint(tj.Pts[0].Pos[0]);
    unsigned int nextWire = wire0 - tj.StepDir;
    
    // check for dead wires
    geo::PlaneID planeID = DecodeCTP(tj.CTP);
    unsigned short ipl = planeID.Plane;
    while(nextWire > tjs.FirstWire[ipl] && nextWire < tjs.LastWire[ipl]) {
      if(tjs.WireHitRange[ipl][nextWire].first >= 0) break;
      nextWire -= tj.StepDir;
    }
    if(nextWire == tjs.LastWire[ipl] - 1) return;
    // clone the first point
    TrajPoint tp = tj.Pts[0];
    // strip off the hits
    tp.Hits.clear(); tp.UseHit.reset();
    // move it to the next wire
    MoveTPToWire(tp, (float)nextWire);
    // find close unused hits near this position
    float maxDelta = 10 * tj.Pts[tj.EndPt[1]].DeltaRMS;
    if(!FindCloseHits(tjs, tp, maxDelta, kUnusedHits)) return;
    if(prt) mf::LogVerbatim("TC")<<" nUnused hits "<<tp.Hits.size()<<" at Pos "<<PrintPos(tjs, tp);
    if(tp.Hits.empty()) return;
     // There are hits on the next wire. Make a copy, reverse it and try
    // to extend it with StepCrawl
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" tp.Hits ";
      for(auto& iht : tp.Hits) myprt<<" "<<PrintHit(tjs.fHits[iht])<<"_"<<tjs.fHits[iht].InTraj;
    } // prt
    //
    // Make a working copy of tj
    Trajectory tjWork = tj;
    // So the first shall be last and the last shall be first
    ReverseTraj(tjs, tjWork);
    // Flag it to use special cuts in StepCrawl
    tjWork.AlgMod[kRevProp] = true;
    // We are doing this probably because the trajectory is stopping.
    // Reduce the number of fitted points to a small number
    unsigned short lastPt = tjWork.Pts.size() - 1;
    if(lastPt < 4) return;
    // update the charge
    float chg = 0;
    float cnt = 0;
    for(unsigned short ii = 0; ii < 4; ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tjWork.Pts[ipt].Chg == 0) continue;
      chg += tjWork.Pts[ipt].Chg;
      ++cnt;
    } // ii
    if(cnt == 0) return;
    if(cnt > 1) tjWork.Pts[lastPt].AveChg = chg / cnt;
    StepCrawl(tjWork);
    if(!fGoodTraj) {
      if(prt) mf::LogVerbatim("TC")<<" ReversePropagate StepCrawl failed";
      return;
    }
    // restore the original direction
    if(tjWork.StepDir != stepDir) ReverseTraj(tjs, tjWork);
    tj = tjWork;
    // re-check for a stopping track
    ChkStop(tj);
    if(prt) {
      mf::LogVerbatim("TC")<<" ReversePropagate success. Outgoing StepDir "<<tj.StepDir;
      if(tj.Pts.size() < 50) PrintTrajectory("RP", tjs, tj, USHRT_MAX);
    }

  } // ReversePropagate
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindJunkTraj()
  {
    // Makes junk trajectories using unassigned hits
    
    if(fPlane > tjs.FirstWire.size() - 1) {
      mf::LogWarning("TC")<<"FindJunkTraj called with invalid fPlane "<<fPlane;
      fQuitAlg = true;
      return;
    }
    
    unsigned int iwire, ifirsthit, ilasthit, iht;
    unsigned int jwire, jfirsthit, jlasthit, jht;
    unsigned int kwire, kfirsthit, klasthit, kht;
    unsigned int fromIndex;
    unsigned int loWire, hiWire;
    float loTime, hiTime;
    bool hitsAdded;
    unsigned short nit, tht, newTjIndex;
    
    // shouldn't have to do this but...
    for(iht = 0; iht < tjs.fHits.size(); ++iht) if(tjs.fHits[iht].InTraj < 0) tjs.fHits[iht].InTraj = 0;
    
    std::vector<unsigned int> tHits;
    // vectors for checking hit consistency
    std::vector<unsigned int> iHit(1), jHit(1);
    bool jtPrt = false;
    for(iwire = tjs.FirstWire[fPlane]; iwire < tjs.LastWire[fPlane] - 1; ++iwire) {
      // skip bad wires or no hits on the wire
      if(tjs.WireHitRange[fPlane][iwire].first < 0) continue;
      jwire = iwire + 1;
      if(tjs.WireHitRange[fPlane][jwire].first < 0) continue;
      ifirsthit = (unsigned int)tjs.WireHitRange[fPlane][iwire].first;
      ilasthit = (unsigned int)tjs.WireHitRange[fPlane][iwire].second;
      jfirsthit = (unsigned int)tjs.WireHitRange[fPlane][jwire].first;
      jlasthit = (unsigned int)tjs.WireHitRange[fPlane][jwire].second;
      for(iht = ifirsthit; iht < ilasthit; ++iht) {
        jtPrt = (iht == debug.Hit);
        if(jtPrt) {
          mf::LogVerbatim("TC")<<"FindJunkTraj: Found debug hit "<<PrintHit(tjs.fHits[iht])<<" fJTMaxHitSep2 "<<fJTMaxHitSep2<<" iht "<<iht<<" jfirsthit "<<jfirsthit<<" jlasthit "<<jlasthit;
        }
        if(tjs.fHits[iht].InTraj != 0) continue;
        if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
        for(jht = jfirsthit; jht < jlasthit; ++jht) {
          if(tjs.fHits[jht].InTraj != 0) continue;
          if(tjs.IgnoreNegChiHits && tjs.fHits[jht].GoodnessOfFit < 0) continue;
          if(prt && HitSep2(tjs, iht, jht) < 100) mf::LogVerbatim("TC")<<" use "<<PrintHit(tjs.fHits[jht])<<" HitSep2 "<<HitSep2(tjs, iht, jht);
          if(HitSep2(tjs, iht, jht) > fJTMaxHitSep2) continue;
          jHit[0] = jht;
          // check for hit width consistency
          if(!TrajHitsOK(tjs, iHit, jHit)) continue;
          tHits.clear();
          // add all hits and flag them
          fromIndex = iht - tjs.fHits[iht].LocalIndex;
          for(kht = fromIndex; kht < fromIndex + tjs.fHits[iht].Multiplicity; ++kht) {
            if(tjs.fHits[kht].InTraj != 0) continue;
            if(HitSep2(tjs, iht, kht) > fJTMaxHitSep2) continue;
            tHits.push_back(kht);
            tjs.fHits[kht].InTraj = -4;
          } // kht
          fromIndex = jht - tjs.fHits[jht].LocalIndex;
          for(kht = fromIndex; kht < fromIndex + tjs.fHits[jht].Multiplicity; ++kht) {
            if(tjs.fHits[kht].InTraj != 0) continue;
            if(HitSep2(tjs, jht, kht) > fJTMaxHitSep2) continue;
            tHits.push_back(kht);
            tjs.fHits[kht].InTraj = -4;
          } // kht
          if(iwire != 0) { loWire = iwire - 1; } else { loWire = 0; }
          if(jwire < tjs.NumWires[fPlane] - 3) { hiWire = jwire + 2; } else { hiWire = tjs.NumWires[fPlane] - 1; }
          hitsAdded = true;
          nit = 0;
          while(hitsAdded && nit < 100) {
            hitsAdded = false;
            for(kwire = loWire; kwire < hiWire + 1; ++kwire) {
              if(tjs.WireHitRange[fPlane][kwire].first < 0) continue;
              kfirsthit = (unsigned int)tjs.WireHitRange[fPlane][kwire].first;
              klasthit = (unsigned int)tjs.WireHitRange[fPlane][kwire].second;
              for(kht = kfirsthit; kht < klasthit; ++kht) {
                if(tjs.fHits[kht].InTraj != 0) continue;
                if(tjs.IgnoreNegChiHits && tjs.fHits[kht].GoodnessOfFit < 0) continue;
                // this shouldn't be needed but do it anyway
                if(std::find(tHits.begin(), tHits.end(), kht) != tHits.end()) continue;
                // re-purpose jHit and check for consistency
                jHit[0] = kht;
                if(!TrajHitsOK(tjs, tHits, jHit)) continue;
                // check w every hit in tHit
                for(tht = 0; tht < tHits.size(); ++tht) {
                  if(jtPrt && HitSep2(tjs, kht, tHits[tht]) < 100) mf::LogVerbatim("TC")<<" kht "<<PrintHit(tjs.fHits[kht])<<" tht "<<PrintHit(tjs.fHits[tHits[tht]])<<" HitSep2 "<<HitSep2(tjs, kht, tHits[tht])<<" cut "<<fJTMaxHitSep2;
                  if(HitSep2(tjs, kht, tHits[tht]) > fJTMaxHitSep2) continue;
                  hitsAdded = true;
                  tHits.push_back(kht);
                  tjs.fHits[kht].InTraj = -4;
                  if(tHits.size() > 50) break;
                  if(kwire > hiWire) hiWire = kwire;
                  if(kwire < loWire) loWire = kwire;
                  break;
                } // tht
              } // kht
              if(jtPrt) mf::LogVerbatim("TC")<<" kwire "<<kwire<<" thits size "<<tHits.size();
            } // kwire
            ++nit;
          } // hitsAdded && nit < 100
          loTime = 1E6; hiTime = 0;
          loWire = USHRT_MAX; hiWire = 0;
          for(tht = 0; tht < tHits.size(); ++tht) {
            if(tjs.fHits[tHits[tht]].WireID.Wire < loWire) loWire = tjs.fHits[tHits[tht]].WireID.Wire;
            if(tjs.fHits[tHits[tht]].WireID.Wire > hiWire) hiWire = tjs.fHits[tHits[tht]].WireID.Wire;
            if(tjs.fHits[tHits[tht]].PeakTime < loTime) loTime = tjs.fHits[tHits[tht]].PeakTime;
            if(tjs.fHits[tHits[tht]].PeakTime > hiTime) hiTime = tjs.fHits[tHits[tht]].PeakTime;
          }
          if(jtPrt) {
            mf::LogVerbatim myprt("TC");
            myprt<<" tHits";
            for(auto tht : tHits) myprt<<" "<<PrintHit(tjs.fHits[tht]);
            myprt<<"\n";
          } // prt
          // See if this is a ghost trajectory
          unsigned short ofTraj = USHRT_MAX;
          if(IsGhost(tHits, ofTraj)) continue;
          MakeJunkTraj(tHits, newTjIndex);
          if(fQuitAlg) return;
          // release any hits that weren't included in a trajectory
          for(auto iht : tHits) if(tjs.fHits[iht].InTraj == -4) tjs.fHits[iht].InTraj = 0;
          if(hitsAdded) break;
        } // jht
      } // iht
    } // iwire
  } // FindJunkTraj

  //////////////////////////////////////////
  void TrajClusterAlg::MakeJunkTraj(std::vector<unsigned int> tHits, unsigned short& newTjIndex)
  {
    
    if(!tjs.UseAlg[kJunkTj]) return;
     // Make a crummy trajectory using the provided hits
    newTjIndex = USHRT_MAX;
    
    if(tHits.size() < 2) return;

    std::vector<std::vector<unsigned int>> tpHits;
    unsigned short ii, iht, ipt;
    
    // Start the trajectory using the first and last hits to
    // define a starting direction. Use the last pass settings
    Trajectory work;
    unsigned short pass = fMinPts.size() - 1;
    if(!StartTraj(work, tHits[0], tHits[tHits.size()-1], pass)) return;
    if(work.ID == debug.WorkID) {
      mf::LogVerbatim("TC")<<" Turning on debug mode in MakeJunkTraj";
      debug.CTP = work.CTP;
      prt = true;
    }
    
    // Do a more complicated specification of TP hits if there
    // are enough of them
    if(tHits.size() > 6) {
      // fit all of the hits to a line
      std::vector<float> x(tHits.size()), y(tHits.size()), yerr2(tHits.size());
      float intcpt, slope, intcpterr, slopeerr, chidof, qtot = 0;
      
      for(ii = 0; ii < tHits.size(); ++ii) {
        iht = tHits[ii];
        x[ii] = tjs.fHits[iht].WireID.Wire;
        y[ii] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
        qtot += tjs.fHits[iht].Integral;
        yerr2[ii] = tjs.fHits[iht].Integral;
      } // ii
      fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
      
      if(prt) mf::LogVerbatim("TC")<<" tHits line fit chidof "<<chidof<<" Angle "<<atan(slope);
      // return without making a junk trajectory
      if(chidof > 900) return;
      // A rough estimate of the trajectory angle
      work.Pts[0].Ang = atan(slope);
      work.Pts[0].Dir[0] = cos(work.Pts[0].Ang);
      work.Pts[0].Dir[1] = sin(work.Pts[0].Ang);
      SetAngleCode(tjs, work.Pts[0]);
      // Rotate the hits into this coordinate system to find the start and end
      // points and general direction
      double cs = cos(-work.Pts[0].Ang);
      double sn = sin(-work.Pts[0].Ang);
      float tAlong, minAlong = 1E6, maxAlong = -1E6;
      float pointSize = 2.1;
      // sort the hits by the distance along the general direction
      std::vector<SortEntry> sortVec(tHits.size());
      SortEntry sortEntry;
      for(ii = 0; ii < x.size(); ++ii) {
        tAlong = cs * x[ii] - sn * y[ii];
        if(tAlong < minAlong) minAlong = tAlong;
        if(tAlong > maxAlong) maxAlong = tAlong;
        sortEntry.index = ii;
        sortEntry.length = tAlong;
        sortVec[ii] = sortEntry;
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), lessThan);
      // make a temp vector
      std::vector<unsigned int> tmp(sortVec.size());
      // overwrite with the sorted values
      for(ii = 0; ii < sortVec.size(); ++ii) tmp[ii] = tHits[sortVec[ii].index];
      tHits = tmp;
      // create a trajectory point at each WSE unit (if there are hits at that point)
      unsigned short npts = (unsigned short)((maxAlong - minAlong) / pointSize);
      // rotate back into normal coordinate system
      if(prt) mf::LogVerbatim("TC")<<" minAlong "<<minAlong<<" maxAlong "<<maxAlong<<" work.Pts[0].Ang "<<work.Pts[0].Ang<<" npts "<<npts;
      if(npts < 2) npts = 2;
      tpHits.resize(npts);
      for(ii = 0; ii < tHits.size(); ++ii) {
        ipt = (unsigned short)((sortVec[ii].length - minAlong) / pointSize);
        if(ipt > npts - 1) ipt = npts - 1;
        if(prt) mf::LogVerbatim("TC")<<"tHit "<<PrintHit(tjs.fHits[tHits[ii]])<<" length "<<sortVec[ii].length<<" ipt "<<ipt<<" Chg "<<(int)tjs.fHits[tHits[ii]].Integral;
        tpHits[ipt].push_back(tHits[ii]);
      }
    }  else {
      // just a few hits. Put each one at a TP in the order that
      // they were found
      tpHits.resize(tHits.size());
      for(ii = 0; ii < tHits.size(); ++ii) {
        tpHits[ii].push_back(tHits[ii]);
      }
    } // tHits.size()
    // make the TPs
    // work.Pts[0] is already defined but it needs hits added
    work.Pts[0].Hits = tpHits[0];
    // set all bits true
    work.Pts[0].UseHit.set();
    DefineHitPos(work.Pts[0]);
    work.Pts[0].Pos = work.Pts[0].HitPos;
    if(prt) PrintTrajectory("MJT", tjs, work, USHRT_MAX);
    // another TP to get the direction
    TrajPoint tpd;
    // make the rest of the TPs
    for(ipt = 1; ipt < tpHits.size(); ++ipt) {
      if(tpHits[ipt].empty()) continue;
      // Use the previous TP as a starting point
      unsigned short lastPt = work.Pts.size() - 1;
      TrajPoint tp = work.Pts[lastPt];
      tp.Step = ipt;
      tp.Hits = tpHits[ipt];
      // use all hits
      tp.UseHit.set();
      DefineHitPos(tp);
      // Just use the hit position as the traj position
      tp.Pos = tp.HitPos;
      if(TrajPointSeparation(work.Pts[ipt-1], tp) < 0.5) continue;
/* This sometimes gives bogus angles. Just use the original TP angle
      // define the direction
      if(!MakeBareTrajPoint(tjs, work.Pts[ipt-1], tp, tpd)) continue;
      tp.Dir = tpd.Dir;
      tp.Ang = tpd.Ang;
      SetAngleCode(tjs, tp);
      if(ipt == 1) {
        work.Pts[0].Dir = tpd.Dir;
        work.Pts[0].Ang = tpd.Ang;
        work.Pts[0].AngleCode = tpd.AngleCode;
      }
*/
      work.Pts.push_back(tp);
      SetEndPoints(tjs, work);
    }
    if(prt) {
      PrintTrajectory("MJT", tjs, work, USHRT_MAX);
    }
    work.AlgMod[kJunkTj] = true;
    fGoodTraj = true;
    // Finally push it onto tjs.allTraj
    fQuitAlg = !StoreTraj(tjs, work);
    if(fQuitAlg) return;
    if(tjs.UseAlg[kChkInTraj]) {
      fQuitAlg = !InTrajOK(tjs, "MJT");
      if(fQuitAlg) {
        mf::LogVerbatim("TC")<<"InTrajOK failed in MakeJunkTraj";
        return;
      }
    }
    // return with a valid index for the new trajectory
    newTjIndex = tjs.allTraj.size() - 1;
  } // MakeJunkTraj

  //////////////////////////////////////////
  void TrajClusterAlg::MatchTrueHits()
  {
    // Matches reco hits to MC true tracks and puts the association into
    // TCHit TruTrkID. This code is almost identical to the first part of MatchTruth.
    
    if(fIsRealData) return;
    if(fMatchTruth[0] < 0) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    tjs.MCPartList.clear();
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;
    
    simb::Origin_t sourceOrigin = simb::kUnknown;
    // partList is the vector of MC particles that we want to use
    tjs.MCPartList.reserve(plist.size());
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* part = (*ipart).second;
      int trackID = part->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(sourcePtclTrackID < 0) {
        if(fMatchTruth[0] == 1) {
          // Look for beam neutrino or single particle
          if(theTruth->Origin() == simb::kBeamNeutrino) {
            fSourceParticleEnergy = 1000 * part->E();
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
            if(fMatchTruth[1] > 2) std::cout<<"Found beam neutrino sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<"\n";
          } // beam neutrino
          if(theTruth->Origin() == simb::kSingleParticle) {
            fSourceParticleEnergy = 1000 * part->E();
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
            if(fMatchTruth[1] > 0) {
              TVector3 dir;
              dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
              dir.SetMag(1);
              std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<" Vx "<<(int)part->Vx()<<" Vy "<<(int)part->Vy()<<" Vz "<<(int)part->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            }
          } // single particle
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            if(fMatchTruth[1] > 2) std::cout<<"True vertex position "<<(int)part->Vx()<<" "<<(int)part->Vy()<<" "<<(int)part->Vz()<<"\n";
          } // sourceOrigin != simb::kUnknown
        } else {
          // look for cosmic rays
          if(theTruth->Origin() == simb::kCosmicRay) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kCosmicRay;
          }
        }
      }
      // ignore anything that has the incorrect origin
      if(theTruth->Origin() != sourceOrigin) continue;
      // ignore processes that aren't a stable final state particle
      if(part->Process() == "neutronInelastic") continue;
      if(part->Process() == "hadElastic") continue;
      // ignore anything that isn't charged
      unsigned short pdg = abs(part->PdgCode());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      tjs.MCPartList.push_back(part);
    } // ipart
    
    tjs.MCPartList.shrink_to_fit();
    
    // vector of (mother, daughter) pairs of partList indices
    std::vector<std::pair<unsigned short, unsigned short>> moda;
    
    if(sourcePtclTrackID > 0) {
      // enter grandmother-daughter pairs for primary electrons if we are using shower finding code
      if(tjs.ShowerTag[0] > 1) {
        for(unsigned short ii = 0; ii < tjs.MCPartList.size(); ++ii) {
          // check for the end of the primary particles
          if(tjs.MCPartList[ii]->Mother() != 0) break;
          // check for an electron
          if(abs(tjs.MCPartList[ii]->PdgCode()) != 11) continue;
          int primElectronTrackID = tjs.MCPartList[ii]->TrackId();
          for(unsigned short jj = ii + 1; jj < tjs.MCPartList.size(); ++jj) {
            int trackID = tjs.MCPartList[jj]->TrackId();
            const simb::MCParticle* gmom = bt->TrackIDToMotherParticle(trackID);
            if(gmom == 0 || gmom->TrackId() != primElectronTrackID) continue;
            moda.push_back(std::make_pair(ii, jj));
          } // jj
        } // ii
      } // Shower finding mode
      // Now enter mother-daughter pairs for soft interactions
      // daughters appear later in the list so reverse iterate
      for(unsigned short ii = 0; ii < tjs.MCPartList.size(); ++ii) {
        unsigned short dpl = tjs.MCPartList.size() - 1 - ii;
        if(tjs.MCPartList[dpl]->Mother() == 0) continue;
        // ignore previous entries
        int trackID = tjs.MCPartList[ii]->TrackId();
        bool skipit = false;
        for(auto& md : moda) if(md.second == trackID) skipit = true;
        if(skipit) continue;
        int motherID = tjs.MCPartList[dpl]->Mother() + sourcePtclTrackID - 1;
        // count the number of daughters
        unsigned short ndtr = 0;
        for(unsigned short jj = 0; jj < tjs.MCPartList.size(); ++jj) {
          // some processes to ignore
          if(tjs.MCPartList[jj]->Process() == "hIoni") continue;
          if(tjs.MCPartList[jj]->Process() == "eIoni") continue;
          if(tjs.MCPartList[jj]->Mother() == tjs.MCPartList[dpl]->Mother()) ++ndtr;
        } // jj
        // require only one daughter
        if(ndtr != 1) continue;
        // Then find the mother index
        unsigned short momIndex = USHRT_MAX;
        for(unsigned short jj = 0; jj < tjs.MCPartList.size(); ++jj) {
          if(tjs.MCPartList[jj]->TrackId() == motherID) {
            momIndex = jj;
            break;
          }
        } // jj
        // Mother not found for some reason
        if(momIndex == USHRT_MAX) continue;
        // ensure that mother and daughter have the same PDG code
        if(tjs.MCPartList[momIndex]->PdgCode() != tjs.MCPartList[dpl]->PdgCode()) continue;
        moda.push_back(std::make_pair(ii, dpl));
      } // ii
    } // sourcePtclTrackID >= 0

    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      TCHit& hit = tjs.fHits[iht];
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)hit.WireID.Plane, (int)hit.WireID.Wire, (int)hit.WireID.TPC, (int)hit.WireID.Cryostat);
      double startTick = hit.PeakTime - hit.RMS;
      double endTick = hit.PeakTime + hit.RMS;
      unsigned short hitTruTrkID = 0;
      // get a list of track IDEs that are close to this hit
      std::vector<sim::TrackIDE> tides;
      bt->ChannelToTrackIDEs(tides, channel, startTick, endTick);
      // Declare a match to the one which has an energy fraction > 0.5
      for(auto itide = tides.begin(); itide != tides.end(); ++itide) {
        if(itide->energyFrac > 0.5) {
          hitTruTrkID = itide->trackID;
          break;
        }
      } // itid
      // not matched (confidently) to a MC track
      if(hitTruTrkID == 0) continue;
      // find out which partList entry corresponds to this track ID
      unsigned short partListIndex;
      for(partListIndex = 0; partListIndex < tjs.MCPartList.size(); ++partListIndex) if(hitTruTrkID == tjs.MCPartList[partListIndex]->TrackId()) break;
      if(partListIndex == tjs.MCPartList.size()) {
//        std::cout<<"MatchTrueHits: Didn't find partList entry for MC Track ID "<<hitTruTrkID<<"\n";
        continue;
      }
      
      // Try to re-assign it to a mother. Note that the mother-daughter pairs
      // are in reverse order, so this loop will transfer all-generation daughters
      // to the (grand) mother
      for(auto& md : moda) if(md.second == partListIndex) partListIndex = md.first;
      hit.MCPartListIndex = partListIndex;
    } // iht

    
    if(fMatchTruth[1] > 1) {
      mf::LogVerbatim myprt("TC");
      myprt<<"part   PDG TrkID MomID KE(MeV)   Process         Trajectory_extent_in_plane \n";
      for(unsigned short ipl = 0; ipl < tjs.MCPartList.size(); ++ipl) {
        unsigned short pdg = abs(tjs.MCPartList[ipl]->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (tjs.MCPartList[ipl]->E() - tjs.MCPartList[ipl]->Mass());
        int motherID = tjs.MCPartList[ipl]->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipl;
        myprt<<std::setw(6)<<tjs.MCPartList[ipl]->PdgCode();
        myprt<<std::setw(6)<<tjs.MCPartList[ipl]->TrackId();
        myprt<<std::setw(6)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<tjs.MCPartList[ipl]->Process();
        // print the extent of the particle in each plane
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          unsigned int fht = UINT_MAX;
          unsigned int lht = 0;
          for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
            if(tjs.fHits[iht].WireID.Plane != plane) continue;
            unsigned short partListIndex = ipl;
            // Look for the real mother
            for(auto& md : moda) if(md.second == partListIndex) partListIndex = md.first;
            if(tjs.fHits[iht].MCPartListIndex != partListIndex) continue;
            if(fht == UINT_MAX) fht = iht;
            lht = iht;
          } // iht
          if(fht == UINT_MAX) continue;
          myprt<<" "<<PrintHitShort(tjs.fHits[fht])<<"-"<<PrintHitShort(tjs.fHits[lht]);
        } // plane
        myprt<<"\n";
      } // ipl
    }

  } // MatchTrueHits

  //////////////////////////////////////////
  void TrajClusterAlg::MatchTruth()
  {
    
    if(fIsRealData) return;
    if(fMatchTruth[0] < 0) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    // set true if there is a reconstructed 3D vertex within 1 cm of the true vertex
    bool nuVtxRecoOK = false;

    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
    fSourceParticleEnergy = -1;
    
    simb::Origin_t sourceOrigin = simb::kUnknown;
    std::vector<simb::MCParticle*> partList;
    // partList is the vector of MC particles that we want to use
    partList.reserve(plist.size());
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      simb::MCParticle* part = (*ipart).second;
      int trackID = part->TrackId();
      art::Ptr<simb::MCTruth> theTruth = bt->TrackIDToMCTruth(trackID);
      if(sourcePtclTrackID < 0) {
        if(fMatchTruth[0] == 1) {
          // Look for beam neutrino or single particle
          if(theTruth->Origin() == simb::kBeamNeutrino) {
            fSourceParticleEnergy = 1000 * part->E(); // in MeV
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
            fNeutrinoEnergy = 1000 * theTruth->GetNeutrino().Nu().E();
            if(fMatchTruth[1] > 2) std::cout<<"Found beam neutrino E = "<<fNeutrinoEnergy<<" sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<"\n";
          }
          if(theTruth->Origin() == simb::kSingleParticle) {
            fSourceParticleEnergy = 1000 * part->E(); // in MeV
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
            if(fMatchTruth[1] > 0) {
              TVector3 dir;
              dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
              dir.SetMag(1);
              std::cout<<"Found single particle sourcePtclTrackID "<<trackID<<" PDG code "<<part->PdgCode()<<" Vx "<<(int)part->Vx()<<" Vy "<<(int)part->Vy()<<" Vz "<<(int)part->Vz()<<" dir "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
            }
          }
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            if(fMatchTruth[1] > 2) std::cout<<" True vertex position "<<(int)part->Vx()<<" "<<(int)part->Vy()<<" "<<(int)part->Vz()<<" energy "<<(int)(1000*part->E())<<"\n";
            for(auto& aVtx3 : tjs.vtx3) {
              fNuVtx_dx->Fill(part->Vx() - aVtx3.X);
              fNuVtx_dy->Fill(part->Vy() - aVtx3.Y);
              fNuVtx_dz->Fill(part->Vz() - aVtx3.Z);
              if(std::abs(part->Vx()-aVtx3.X) < 1 && std::abs(part->Vy()-aVtx3.Y) < 1 && std::abs(part->Vz()-aVtx3.Z) < 1) {
                nuVtxRecoOK = true;
                float score = 0;
                for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
                  if(aVtx3.Vtx2ID[ipl] == 0) continue;
                  unsigned short iv2 = aVtx3.Vtx2ID[ipl] - 1;
                  score += tjs.vtx[iv2].Score;
                } // ipl
                fNuVtx_Score->Fill(score);
                fNuVtx_Enu_Score_p->Fill(fNeutrinoEnergy, score);
              }
            } // aVtx3
          } // sourceOrigin != simb::kUnknown
        } else {
          // look for cosmic rays
          if(theTruth->Origin() == simb::kCosmicRay) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kCosmicRay;
          }
        }
      }
      // ignore anything that has the incorrect origin
      if(theTruth->Origin() != sourceOrigin) continue;
      // ignore processes that aren't a stable final state particle
      if(part->Process() == "neutronInelastic") continue;
      if(part->Process() == "hadElastic") continue;
      // ignore anything that isn't charged
      unsigned short pdg = abs(part->PdgCode());
      bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
      if(!isCharged) continue;
      partList.push_back(part);
    } // ipart
    
    if(fMatchTruth[1] > 2) {
      for(unsigned int ii = 0; ii < partList.size(); ++ii) {
        int trackID = partList[ii]->TrackId();
        const simb::MCParticle* gmom = bt->TrackIDToMotherParticle(trackID);
        std::cout<<ii<<" PDG Code  "<<partList[ii]->PdgCode()<<" TrackId "<<trackID<<" Mother  "<<partList[ii]->Mother()<<" Grandmother "<<gmom->TrackId()<<" Process "<<partList[ii]->Process()<<"\n";
      } // ii
    }
    
    // vector of (mother, daughter) pairs of TrackIds
    std::vector<std::pair<int, int>> moda;
    
    if(sourcePtclTrackID > 0) {
      // enter grandmother-daughter pairs for primary electrons if we are using shower finding code
      if(tjs.ShowerTag[0] > 1) {
        for(unsigned short ii = 0; ii < partList.size(); ++ii) {
          // check for the end of the primary particles
          if(partList[ii]->Mother() != 0) break;
          // check for an electron
          if(abs(partList[ii]->PdgCode()) != 11) continue;
          int primElectronTrackID = partList[ii]->TrackId();
          for(unsigned short jj = ii + 1; jj < partList.size(); ++jj) {
            int trackID = partList[jj]->TrackId();
            const simb::MCParticle* gmom = bt->TrackIDToMotherParticle(trackID);
            if(gmom == 0 || gmom->TrackId() != primElectronTrackID) continue;
            moda.push_back(std::make_pair(primElectronTrackID, trackID));
          } // jj
        } // ii
      } // Shower finding mode
      // Now enter mother-daughter pairs for soft interactions
      // daughters appear later in the list so reverse iterate
      for(unsigned short ii = 0; ii < partList.size(); ++ii) {
        unsigned short dpl = partList.size() - 1 - ii;
        if(partList[dpl]->Mother() == 0) continue;
        // ignore previous entries
        int trackID = partList[ii]->TrackId();
        bool skipit = false;
        for(auto& md : moda) if(md.second == trackID) skipit = true;
        if(skipit) continue;
        int motherID = partList[dpl]->Mother() + sourcePtclTrackID - 1;
        // count the number of daughters
        unsigned short ndtr = 0;
        for(unsigned short jj = 0; jj < partList.size(); ++jj) {
          // some processes to ignore
          if(partList[jj]->Process() == "hIoni") continue;
          if(partList[jj]->Process() == "eIoni") continue;
          if(partList[jj]->Mother() == partList[dpl]->Mother()) ++ndtr;
        } // jj
        // require only one daughter
        if(ndtr != 1) continue;
        // Then find the mother index
        unsigned short momIndex = USHRT_MAX;
        for(unsigned short jj = 0; jj < partList.size(); ++jj) {
          if(partList[jj]->TrackId() == motherID) {
            momIndex = jj;
            break;
          }
        } // jj
        // Mother not found for some reason
        if(momIndex == USHRT_MAX) continue;
        // ensure that mother and daughter have the same PDG code
        if(partList[momIndex]->PdgCode() != partList[dpl]->PdgCode()) continue;
        moda.push_back(std::make_pair(partList[momIndex]->TrackId(), partList[dpl]->TrackId()));
      } // ii
    } // sourcePtclTrackID >= 0

    if(fMatchTruth[1] > 2 && !moda.empty()) {
      std::cout<<"Mother-Daughter track IDs\n";
      unsigned short cnt = 0;
      for(auto& md : moda) {
        std::cout<<" "<<md.first<<"-"<<md.second;
        ++cnt;
        if(!(cnt % 20)) std::cout<<"\n";
      } // md
      std::cout<<"\n";
    } // fMatchTruth[1] > 2

    // Match all hits to the truth. Put the MC track ID in a temp vector
    std::vector<int> hitTruTrkID(tjs.fHits.size());
    // Prepare to count of the number of hits matched to each MC Track in each plane
    std::vector<std::vector<unsigned short>> nMatchedHitsInPartList(plist.size());
    for(unsigned short ipl = 0; ipl < plist.size(); ++ipl) nMatchedHitsInPartList[ipl].resize(tjs.NumPlanes);
    // and make a list of the TJs and hit count for each MC Track
    std::vector<std::vector<std::array<unsigned short, 2>>> nMatchedHitsInTj(partList.size());

    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      TCHit& hit = tjs.fHits[iht];
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)hit.WireID.Plane, (int)hit.WireID.Wire, (int)hit.WireID.TPC, (int)hit.WireID.Cryostat);
      double startTick = hit.PeakTime - hit.RMS;
      double endTick = hit.PeakTime + hit.RMS;
      unsigned short plane = tjs.fHits[iht].WireID.Plane;
      // get a list of track IDEs that are close to this hit
      std::vector<sim::TrackIDE> tides;
      bt->ChannelToTrackIDEs(tides, channel, startTick, endTick);
      // Declare a match to the one which has an energy fraction > 0.5
      for(auto itide = tides.begin(); itide != tides.end(); ++itide) {
        if(itide->energyFrac > 0.5) {
          hitTruTrkID[iht] = itide->trackID;
          break;
        }
      } // itid
      // not matched (confidently) to a MC track
      if(hitTruTrkID[iht] == 0) continue;

      // Try to re-assign it to a mother. Note that the mother-daughter pairs
      // are in reverse order, so this loop will transfer all-generation daughters
      // to the (grand) mother
      for(auto& md : moda) if(md.second == hitTruTrkID[iht]) hitTruTrkID[iht] = md.first;

      // count the number of matched hits for each MC track in each plane
      for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
        if(hitTruTrkID[iht] == partList[ipl]->TrackId()) {
          ++nMatchedHitsInPartList[ipl][plane];
          if(tjs.fHits[iht].InTraj > 0) {
            unsigned short itj = tjs.fHits[iht].InTraj - 1;
            bool gotit = false;
            for(auto& hitInTj : nMatchedHitsInTj[ipl]) {
              if(hitInTj[0] == itj) {
                ++hitInTj[1];
                gotit = true;
              }
            } //  hitInTj
            if(!gotit) {
              std::array<unsigned short, 2> tmp {itj, 1};
              nMatchedHitsInTj[ipl].push_back(tmp);
            }
          } // inTraj > 0
        } // hit matched to partList
      } // ipl
    } // iht
    
    // remove partList elements that have no matched hits
    std::vector<simb::MCParticle*> newPartList;
    std::vector<std::vector<unsigned short>> newnMatchedHitsInPartList;
    std::vector<std::vector<std::array<unsigned short, 2>>> newnMatchedHitsInTj;
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      unsigned short nht = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) nht += nMatchedHitsInPartList[ipl][plane];
      if(nht == 0) continue;
      newPartList.push_back(partList[ipl]);
      newnMatchedHitsInPartList.push_back(nMatchedHitsInPartList[ipl]);
      newnMatchedHitsInTj.push_back(nMatchedHitsInTj[ipl]);
    } // ipl
    partList = newPartList;
    nMatchedHitsInPartList = newnMatchedHitsInPartList;
    nMatchedHitsInTj = newnMatchedHitsInTj;
    
    // count the number of primary tracks that have at least 3 hits in at least 2 planes
    nTruPrimary = 0;
    nTruPrimaryOK = 0;
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      if(partList[ipl]->Mother() != 0) continue;
      ++nTruPrimary;
      unsigned short nInPln = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(nMatchedHitsInPartList[ipl][plane] > 2) ++nInPln;
      } // plane
      if(nInPln > 1) ++nTruPrimaryOK;
    } // ipl
    
    if(nTruPrimaryOK > 1) {
      // More than one reconstructable primaries so there must a reconstructable neutrino vertex
      ++nTruPrimaryVtxOK;
      // was it reconstructed?
      if(nuVtxRecoOK) ++nTruPrimaryVtxReco;
      if(fSourceParticleEnergy > 0 && !nuVtxRecoOK) mf::LogVerbatim("TC")<<"BadVtx fSourceParticleEnergy "<<std::fixed<<std::setprecision(2)<<fSourceParticleEnergy<<" events processed "<<fEventsProcessed;
    }
    
    if(fMatchTruth[1] > 1) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Number of primary particles "<<nTruPrimary<<" Number reconstructable "<<nTruPrimaryOK<<" Found neutrino vertex? "<<nuVtxRecoOK<<"\n";
      myprt<<"part   PDG TrkID MomID KE(MeV)   Process         Trajectory_extent_in_plane \n";
      for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
        unsigned short pdg = abs(partList[ipl]->PdgCode());
        bool isCharged = (pdg == 11) || (pdg == 13) || (pdg == 211) || (pdg == 321) || (pdg == 2212);
        if(!isCharged) continue;
        // Kinetic energy in MeV
        int TMeV = 1000 * (partList[ipl]->E() - partList[ipl]->Mass());
        int motherID = partList[ipl]->Mother() + sourcePtclTrackID - 1;
        myprt<<std::setw(4)<<ipl;
        myprt<<std::setw(6)<<partList[ipl]->PdgCode();
        myprt<<std::setw(6)<<partList[ipl]->TrackId();
        myprt<<std::setw(6)<<motherID;
        myprt<<std::setw(6)<<TMeV;
        myprt<<std::setw(20)<<partList[ipl]->Process();
        // print the extent of the particle in each plane
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          unsigned int fht = UINT_MAX;
          unsigned int lht = 0;
          for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
            if(tjs.fHits[iht].WireID.Plane != plane) continue;
            unsigned short momTrackID = partList[ipl]->TrackId();
            // Look for the real mother
            for(auto& md : moda) if(md.second == momTrackID) momTrackID = md.first;
            if(hitTruTrkID[iht] != momTrackID) continue;
            if(fht == UINT_MAX) fht = iht;
            lht = iht;
          } // iht
          if(fht == UINT_MAX) continue;
          myprt<<" "<<PrintHitShort(tjs.fHits[fht])<<"-"<<PrintHitShort(tjs.fHits[lht]);
        } // plane
        myprt<<"\n";
      } // ipl
    }

    // Declare a TJ - partlist match for the trajectory which has the most true hits.
    // another temp vector for the one-to-one match
    std::vector<std::vector<unsigned short>> partListToTjID(partList.size());
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) partListToTjID[ipl].resize(tjs.NumPlanes);
    
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(nMatchedHitsInPartList[ipl][plane] < 2) continue;
        unsigned short mostHits = 0;
        unsigned short tjWithMostHits = USHRT_MAX;
        for(unsigned short ii = 0; ii < nMatchedHitsInTj[ipl].size(); ++ii) {
          unsigned short itj = nMatchedHitsInTj[ipl][ii][0];
          geo::PlaneID planeID = DecodeCTP(tjs.allTraj[itj].CTP);
          // ensure we only check Tjs in the correct plane
          if(planeID.Plane != plane) continue;
          unsigned short nMatHits = nMatchedHitsInTj[ipl][ii][1];
          if(nMatHits > mostHits) {
            mostHits = nMatHits;
            tjWithMostHits = itj;
          }
        } // ii
        if(tjWithMostHits == USHRT_MAX) continue;
        // Count the total number of used hits in the TJ
        auto tmp = PutTrajHitsInVector(tjs.allTraj[tjWithMostHits], kUsedHits);
        if(tjs.allTraj[tjWithMostHits].ParentTrajID > 0) {
          // This is a daughter trajectory that has more truth matched hits than the parent
          unsigned short ptj = tjs.allTraj[tjWithMostHits].ParentTrajID - 1;
          // add the parent hits to the vector of daughter hits
          auto ptmp = PutTrajHitsInVector(tjs.allTraj[ptj], kUsedHits);
          tmp.insert(tmp.end(), ptmp.begin(), ptmp.end());
          // revise the mostHits count. To do this we need to find the parent tj index in nMatchedHitsInTj
          for(unsigned short ii = 0; ii < nMatchedHitsInTj[ipl].size(); ++ii) {
            unsigned short itj = nMatchedHitsInTj[ipl][ii][0];
            if(itj == ptj) {
              mostHits += nMatchedHitsInTj[ipl][ii][1];
              // re-direct the calculation to the parent
              tjWithMostHits = ptj;
              break;
            } // found the parent tj
          } // ii
        } // deal with daughters
        // count the number matched to a true particle
        float nTjHits = 0;
        for(auto& iht : tmp) if(hitTruTrkID[iht] > 0) ++nTjHits;
        float nTruHits = nMatchedHitsInPartList[ipl][plane];
        float nTjTruRecHits = mostHits;
        float eff = nTjTruRecHits / nTruHits;
        float pur = nTjTruRecHits / nTjHits;
        float effpur = eff * pur;
        // This overwrites any previous match that has poorer efficiency * purity
        if(effpur > tjs.allTraj[tjWithMostHits].EffPur) {
          tjs.allTraj[tjWithMostHits].MCPartListIndex = ipl;
          tjs.allTraj[tjWithMostHits].EffPur = effpur;
          partListToTjID[ipl][plane] = tjs.allTraj[tjWithMostHits].ID;
        }
      } // plane
    } // ipl
    
    // Update the EP sums
    for(unsigned short ipl = 0; ipl < partList.size(); ++ipl) {
      float TMeV = 1000 * (partList[ipl]->E() - partList[ipl]->Mass());
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        // require at least 2 matched hits
        if(nMatchedHitsInPartList[ipl][plane] < 2) continue;
        unsigned short pdgIndex = PDGCodeIndex(tjs, partList[ipl]->PdgCode());
        // count the number of EP sums for this PDG code
        EPSums[pdgIndex] += TMeV;
        ++EPCnts[pdgIndex];
        // find the first and last matched hit in this plane
        unsigned int fht = UINT_MAX;
        unsigned int lht = 0;
        // find the first and last matched hit in this plane
        for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
          if(tjs.fHits[iht].WireID.Plane != plane) continue;
          unsigned short momTrackID = partList[ipl]->TrackId();
          // Look for the real mother
          for(auto& md : moda) if(md.second == momTrackID) momTrackID = md.first;
          if(hitTruTrkID[iht] != momTrackID) continue;
          if(fht == UINT_MAX) fht = iht;
          lht = iht;
        } // iht
        if(fht == UINT_MAX) continue;
        if(partListToTjID[ipl][plane] == 0) {
          // Enter 0 in the profile histogram
          fEP_T[pdgIndex]->Fill(TMeV, 0);
          if(nMatchedHitsInPartList[ipl][plane] > fMatchTruth[3]) {
            mf::LogVerbatim myprt("TC");
            myprt<<"pdgIndex "<<pdgIndex<<" BadEP TMeV "<<(int)TMeV<<" No matched trajectory to partList["<<ipl<<"]";
            myprt<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipl][plane];
            myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht])<<" events processed "<<fEventsProcessed;
          }
          continue;
        }
        unsigned short itj = partListToTjID[ipl][plane] - 1;
        EPTSums[pdgIndex] += TMeV * tjs.allTraj[itj].EffPur;
        fEP_T[pdgIndex]->Fill(TMeV, tjs.allTraj[itj].EffPur);
        // print out some debugging information if the EP was pitiful and the number of matched hits is large
        if(tjs.allTraj[itj].EffPur < fMatchTruth[2] && nMatchedHitsInPartList[ipl][plane] > fMatchTruth[3]) {
          mf::LogVerbatim myprt("TC");
          myprt<<"pdgIndex "<<pdgIndex<<" BadEP "<<std::fixed<<std::setprecision(2)<<tjs.allTraj[itj].EffPur;
          myprt<<" TMeV "<<(int)TMeV<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipl][plane];
          myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht])<<" events processed "<<fEventsProcessed;
          // print alg names
          for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tjs.allTraj[itj].AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        }
        // check for a bad match to a primary electron shower
        if(tjs.ShowerTag[0] > 1 && pdgIndex == 0 && partList[ipl]->Mother() == 0 && fSourceParticleEnergy > 50) {
          Trajectory& ptj = tjs.allTraj[itj];
          // determine if this is a parent of a shower Tj
          int dtrID = 0;
          for(auto& tj : tjs.allTraj) {
            if(!tj.AlgMod[kShowerTj]) continue;
            if(tj.ParentTrajID == ptj.ID) {
              dtrID = tj.ID;
              break;
            }
          } // tj
          if(dtrID == 0) mf::LogVerbatim("TC")<<"BadShower Primary electron -> Traj "<<ptj.ID<<"_"<<ptj.CTP<<" Wrong shower parent. Events processed "<<fEventsProcessed;
        } // check primary electron shower
        // histogram the MC-reco stopping point difference
        unsigned short endPt = tjs.allTraj[itj].EndPt[0];
        float recoWire0 = tjs.allTraj[itj].Pts[endPt].Pos[0];
        endPt = tjs.allTraj[itj].EndPt[1];
        float recoWire1 = tjs.allTraj[itj].Pts[endPt].Pos[0];
        float trueFirstWire = tjs.fHits[fht].WireID.Wire;
        float trueLastWire = tjs.fHits[lht].WireID.Wire;
        // decide which ends should be compared
        if(std::abs(recoWire0 - trueFirstWire) < std::abs(recoWire1 - trueFirstWire)) {
          fdWire[pdgIndex]->Fill(recoWire0 - trueFirstWire);
          fdWire[pdgIndex]->Fill(recoWire1 - trueLastWire);
        } else {
          fdWire[pdgIndex]->Fill(recoWire1 - trueFirstWire);
          fdWire[pdgIndex]->Fill(recoWire0 - trueLastWire);
        }
      } // plane
    } // ipl

    // match 2D vertices (crudely)
    for(auto& tj : tjs.allTraj) {
      // obsolete vertex
      if(tj.AlgMod[kKilled]) continue;
      // require a truth match
      if(tj.MCPartListIndex == USHRT_MAX) continue;
      // ignore electrons unless it is a primary electron
      auto& mcp = tjs.MCPartList[tj.MCPartListIndex];
      int pdg = abs(mcp->PdgCode());
      if(pdg == 11 && mcp->Mother() != 0) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == 0) continue;
        VtxStore& vx2 = tjs.vtx[tj.VtxID[end]-1];
        vx2.Stat[kVtxTruMatch] = true;
      } // end
    } // vx2

  } // MatchTruth

  ////////////////////////////////////////////////
  void TrajClusterAlg::AddLAHits(Trajectory& tj, unsigned short ipt, bool& sigOK)
  {
    // Very Large Angle version of AddHits to be called for the last angle range

    if(ipt > tj.Pts.size() - 1) return;
    TrajPoint& tp = tj.Pts[ipt];
    tp.Hits.clear();
    tp.UseHit.reset();
    sigOK = false;

    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short ipl = planeID.Plane;

    // look at adjacent wires for larger angle trajectories
    // We will check the most likely wire first
    std::vector<int> wires(1);
    wires[0] = std::nearbyint(tp.Pos[0]);
    if(wires[0] < 0 || wires[0] > (int)tjs.LastWire[ipl]-1) return;
    
    if(tp.AngleCode != 2) {
      mf::LogVerbatim("TC")<<"AddLAHits called with a bad angle code. "<<tp.AngleCode<<" Don't do this";
      return;
    }
    // and the adjacent wires next in the most likely order only
    // after the first two points have been defined
    if(ipt > 1) {
      if(tp.Dir[0] > 0) {
        if(wires[0] < (int)tjs.LastWire[ipl]-1) wires.push_back(wires[0] + 1);
        if(wires[0] > 0) wires.push_back(wires[0] - 1);
       } else {
        if(wires[0] > 0) wires.push_back(wires[0] - 1);
        if(wires[0] < (int)tjs.LastWire[ipl]-1) wires.push_back(wires[0] + 1);
      }
    } // ipt > 0 ...
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" AddLAHits: Pos "<<PrintPos(tjs, tp)<<" tp.AngleCode "<<tp.AngleCode<<" Wires under consideration";
      for(auto& wire : wires) myprt<<" "<<wire;
    }
    
    // a temporary tp that we can move around
    TrajPoint ltp = tp;
    // do this while testing
    sigOK = false;
    
    tp.Hits.clear();
    std::array<int, 2> wireWindow;
    std::array<float, 2> timeWindow;
    float pos1Window = fVLAStepSize/2;
    timeWindow[0] = ltp.Pos[1] - pos1Window;
    timeWindow[1] = ltp.Pos[1] + pos1Window;
    // Put the existing hits in to a vector so we can ensure that they aren't added again
    std::vector<unsigned int> oldHits = PutTrajHitsInVector(tj, kAllHits);
    
    for(unsigned short ii = 0; ii < wires.size(); ++ii) {
      int wire = wires[ii];
      if(wire < 0 || wire > (int)tjs.LastWire[ipl]) continue;
      // Assume a signal exists on a dead wire
      if(tjs.WireHitRange[fPlane][wire].first == -1) sigOK = true;
      if(tjs.WireHitRange[fPlane][wire].first < 0) continue;
      wireWindow[0] = wire;
      wireWindow[1] = wire;
      bool hitsNear;
      // Look for hits using the requirement that the timeWindow overlaps with the hit StartTick and EndTick
      std::vector<unsigned int> closeHits = FindCloseHits(tjs, wireWindow, timeWindow, ipl, kAllHits, fExpectNarrowHits, hitsNear);
      if(hitsNear) sigOK = true;
      for(auto& iht : closeHits) {
        // Ensure that none of these hits are already used by this trajectory
        if(tjs.fHits[iht].InTraj == tj.ID) continue;
        // or in another trajectory in any previously added point
        if(std::find(oldHits.begin(), oldHits.end(), iht) != oldHits.end()) continue;
        if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
        tp.Hits.push_back(iht);
      }
    } // ii
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" LAPos "<<PrintPos(tjs, ltp)<<" Tick window "<<(int)(timeWindow[0]/tjs.UnitsPerTick)<<" to "<<(int)(timeWindow[1]/tjs.UnitsPerTick);
      for(auto& iht : tp.Hits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
    } // prt
    
    // no hits found
    if(tp.Hits.empty()) return;
    
    if(tp.Hits.size() > 16) tp.Hits.resize(16);
    
    tp.UseHit.reset();
    
    if(tjs.UseAlg[kStopAtTj]) {
      // don't continue if we have run into another trajectory that has a similar angle
      unsigned short nAvailable = 0;
      unsigned int otherTjHit = INT_MAX;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tjs.fHits[tp.Hits[ii]].InTraj > 0) {
          otherTjHit = tp.Hits[ii];
          continue;
        }
        ++nAvailable;
      } // ii
      if(nAvailable == 0 && otherTjHit != UINT_MAX) {
        // get the trajectory index
        unsigned short otherTj = tjs.fHits[otherTjHit].InTraj - 1;
        Trajectory& otj = tjs.allTraj[otherTj];
        // find out what point the hit is in
        unsigned short atPt = USHRT_MAX;
        for(unsigned short ipt = 0; ipt < otj.Pts.size(); ++ipt) {
          for(auto& iht : otj.Pts[ipt].Hits) {
            if(iht == otherTjHit) {
              atPt = ipt;
              break;
            } // iht == otherTjHit
          } // iht
          if(atPt != USHRT_MAX) break;
        } // ipt
        if(atPt != USHRT_MAX && DeltaAngle(tp.Ang, otj.Pts[atPt].Ang) < 0.1) {
          if(prt) mf::LogVerbatim("TC")<<" Found a VLA merge candidate trajectory "<<otj.ID<<". Set StopFlag[kAtTj] and stop stepping";
          tj.StopFlag[1][kAtTj] = true;
          return;
        } // atPt is valid
      } // nAvailable == 0 &&
    } // stop at Tj
    
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      unsigned int iht = tp.Hits[ii];
      if(tjs.fHits[iht].InTraj != 0) continue;
      tp.UseHit[ii] = true;
      tjs.fHits[iht].InTraj = tj.ID;
    } // ii
    DefineHitPos(tp);
    SetEndPoints(tjs, tj);
    UpdateAveChg(tj);
 
  } // AddLAHits

  ////////////////////////////////////////////////
  void TrajClusterAlg::AddHits(Trajectory& tj, unsigned short ipt, bool& sigOK)
  {
    // Try to add hits to the trajectory point ipt on the supplied
    // trajectory

    // assume failure
    sigOK = false;

    if(fPlane > tjs.FirstWire.size() - 1) {
      mf::LogWarning("TC")<<"AddHits called with invalid fPlane "<<fPlane;
      fQuitAlg = true;
      return;
    }

    if(tj.Pts.empty()) return;
    if(ipt > tj.Pts.size() - 1) return;
    
    // Call large angle hit finding if the last tp is large angle
    if(tj.Pts[ipt].AngleCode == 2) {
      AddLAHits(tj, ipt, sigOK);
      return;
    }
    
    std::vector<unsigned int> closeHits;

    unsigned int lastPtWithUsedHits = tj.EndPt[1];
    TrajPoint& tp = tj.Pts[ipt];

    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if(wire < tjs.FirstWire[fPlane] || wire > tjs.LastWire[fPlane]-1) return;
    // Move the TP to this wire
    MoveTPToWire(tp, (float)wire);
    
    // find the projection error to this point. Note that if this is the first
    // TP, lastPtWithUsedHits = 0, so the projection error is 0
    float dw = tp.Pos[0] - tj.Pts[lastPtWithUsedHits].Pos[0];
    float dt = tp.Pos[1] - tj.Pts[lastPtWithUsedHits].Pos[1];
    float dpos = sqrt(dw * dw + dt * dt);
    float projErr = dpos * tj.Pts[lastPtWithUsedHits].AngErr;
    // Add this to the Delta RMS factor and construct a cut
    float deltaCut = 3 * (projErr + tp.DeltaRMS);
    
    deltaCut *= fProjectionErrFactor;
    if(prt) mf::LogVerbatim("TC")<<" AddHits: calculated deltaCut "<<deltaCut;
    
    if(deltaCut < 1) deltaCut = 1;
    if(deltaCut > 4) deltaCut = 4;

    // TY: open it up for RevProp, since we might be following a stopping track
    if(tj.AlgMod[kRevProp]) deltaCut *= 2;
    
    // loosen up a bit if we just passed a block of dead wires
    if(abs(dw) > 20 && DeadWireCount(tjs, tp.Pos[0], tj.Pts[lastPtWithUsedHits].Pos[0], tj.CTP) > 10) deltaCut *= 2;
    
    // Create a larger cut to use in case there is nothing close
    float bigDelta = 2 * deltaCut;
    unsigned int imBig = UINT_MAX;
    tp.Delta = deltaCut;
    
    // projected time in ticks for testing the existence of a hit signal
    raw::TDCtick_t rawProjTick = (float)(tp.Pos[1] / tjs.UnitsPerTick);
    if(prt) {
      mf::LogVerbatim("TC")<<" AddHits: wire "<<wire<<" tp.Pos[0] "<<tp.Pos[0]<<" projTick "<<rawProjTick<<" deltaRMS "<<tp.DeltaRMS<<" tp.Dir[0] "<<tp.Dir[0]<<" deltaCut "<<deltaCut<<" dpos "<<dpos<<" projErr "<<projErr<<" ExpectedHitsRMS "<<ExpectedHitsRMS(tp);
    }
    
    std::vector<unsigned int> hitsInMultiplet;
    unsigned short localIndex;
    
    geo::PlaneID planeID = DecodeCTP(tj.CTP);
    unsigned int ipl = planeID.Plane;
    if(wire > tjs.LastWire[ipl]) return;
    // Assume a signal exists on a dead wire
    if(tjs.WireHitRange[ipl][wire].first == -1) sigOK = true;
    if(tjs.WireHitRange[ipl][wire].first < 0) return;
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
    float fwire = wire;
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      if(tjs.fHits[iht].InTraj == tj.ID) continue;
      if(rawProjTick > tjs.fHits[iht].StartTick && rawProjTick < tjs.fHits[iht].EndTick) sigOK = true;
      if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
      float ftime = tjs.UnitsPerTick * tjs.fHits[iht].PeakTime;
      float delta = PointTrajDOCA(tjs, fwire, ftime, tp);
      float dt = std::abs(ftime - tp.Pos[1]);
      GetHitMultiplet(iht, hitsInMultiplet, localIndex);
      if(prt && delta < 100 && dt < 100) {
        mf::LogVerbatim myprt("TC");
        myprt<<"  iht "<<iht;
        myprt<<" "<<tjs.fHits[iht].WireID.Plane<<":"<<PrintHit(tjs.fHits[iht]);
        myprt<<" delta "<<std::fixed<<std::setprecision(2)<<delta<<" deltaCut "<<deltaCut<<" dt "<<dt;
        myprt<<" BB Mult "<<hitsInMultiplet.size()<<" localIndex "<<localIndex<<" RMS "<<std::setprecision(1)<<tjs.fHits[iht].RMS;
        myprt<<" Chi "<<std::setprecision(1)<<tjs.fHits[iht].GoodnessOfFit;
        myprt<<" InTraj "<<tjs.fHits[iht].InTraj;
        myprt<<" Chg "<<(int)tjs.fHits[iht].Integral;
        myprt<<" Signal? "<<sigOK;
      }
      if(tjs.fHits[iht].InTraj == 0 && delta < bigDelta && hitsInMultiplet.size() < 3 && !tj.AlgMod[kRevProp]) {
        // An available hit that is just outside the window that is not part of a large multiplet
        bigDelta = delta;
        imBig = iht;
      }
      if(delta > deltaCut) continue;
      if(std::find(closeHits.begin(), closeHits.end(), iht) != closeHits.end()) continue;
      closeHits.push_back(iht);
      if(hitsInMultiplet.size() > 1) {
        // include all the hits in a multiplet
        for(auto& jht : hitsInMultiplet) {
          if(tjs.fHits[jht].InTraj == tj.ID) continue;
          if(std::find(closeHits.begin(), closeHits.end(), jht) != closeHits.end()) continue;
          closeHits.push_back(jht);
        } // jht
      } // multiplicity > 1
    } // iht

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"closeHits ";
      for(auto iht : closeHits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
      if(imBig < tjs.fHits.size()) {
        myprt<<" imBig "<<PrintHit(tjs.fHits[imBig]);
      } else {
        myprt<<" imBig "<<imBig;
      }
    }
    if(closeHits.empty() && imBig == UINT_MAX) {
      if(prt) mf::LogVerbatim("TC")<<" no signal on any wire at tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" tick "<<(int)tp.Pos[1]/tjs.UnitsPerTick<<" closeHits size "<<closeHits.size();
      return;
    }
    if(imBig < tjs.fHits.size() && closeHits.empty()) {
      closeHits.push_back(imBig);
      if(prt) mf::LogVerbatim("TC")<<" Added bigDelta hit "<<PrintHit(tjs.fHits[imBig])<<" w delta = "<<bigDelta;
    }
    if(!closeHits.empty()) sigOK = true;
    if(!sigOK) return;
    tp.Hits.insert(tp.Hits.end(), closeHits.begin(), closeHits.end());
    if(tp.Hits.size() > 16) {
      // Actually this is a hopelessly messy region that we should ignore
      tp.Hits.clear();
      tp.Chg = 0;
      return;
    }
    // reset UseHit and assume that none of these hits will be used (yet)
    tp.UseHit.reset();
    // decide which of these hits should be used in the fit. Use a generous maximum delta
    // and require a charge check if we'not just starting out
    bool useChg = true;
    FindUseHits(tj, ipt, 10, useChg);
    DefineHitPos(tp);
    SetEndPoints(tjs, tj);
    if(prt) mf::LogVerbatim("TC")<<" number of close hits "<<closeHits.size()<<" used hits "<<NumHitsInTP(tp, kUsedHits);
  } // AddHits
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindUseHits(Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg)
  {
    // Hits have been associated with trajectory point ipt but none are used. Here we will
    // decide which hits to use.
    
    if(ipt > tj.Pts.size() - 1) return;
    TrajPoint& tp = tj.Pts[ipt];
    
    if(tp.Hits.empty()) return;

    // don't check charge when starting out
    if(ipt < 5) useChg = false; 
    float chgPullCut = 1000;
    if(useChg) chgPullCut = fChargeCuts[0];
    // open it up for RevProp, since we might be following a stopping track
    if(tj.AlgMod[kRevProp]) chgPullCut *= 2;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"FUH:  maxDelta "<<maxDelta<<" useChg requested "<<useChg<<" Norm AveChg "<<(int)tp.AveChg<<" tj.ChgRMS "<<tj.ChgRMS<<" chgPullCut "<<chgPullCut<<" TPHitsRMS "<<(int)TPHitsRMSTick(tjs, tp, kUnusedHits)<<" ExpectedHitsRMS "<<(int)ExpectedHitsRMS(tp)<<" AngCode "<<tp.AngleCode;
    }

    // inverse of the path length for normalizing hit charge to 1 WSE unit
    float pathInv = std::abs(tp.Dir[0]);
    if(pathInv < 0.05) pathInv = 0.05;
   
    // Find the hit that has the smallest delta and the number of available hits
    tp.Delta = maxDelta;
    float delta;
    unsigned short imbest = USHRT_MAX;
    std::vector<float> deltas(tp.Hits.size(), 100);
    // keep track of the best delta - even if it is used
    float bestDelta = maxDelta;
    unsigned short nAvailable = 0;
    unsigned short imBadRecoHit = USHRT_MAX;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      tp.UseHit[ii] = false;
      unsigned int iht = tp.Hits[ii];
      delta = PointTrajDOCA(tjs, iht, tp);
      if(delta < bestDelta) bestDelta = delta;
      if(tjs.fHits[iht].InTraj > 0) continue;
      if(tjs.fHits[iht].GoodnessOfFit < 0 || tjs.fHits[iht].GoodnessOfFit > 100) imBadRecoHit = ii;
      ++nAvailable;
      if(prt) {
        if(useChg) {
          if(prt) mf::LogVerbatim("TC")<<" "<<ii<<"  "<<PrintHit(tjs.fHits[iht])<<" delta "<<delta<<" Norm Chg "<<(int)(tjs.fHits[iht].Integral * pathInv);
        } else {
          if(prt) mf::LogVerbatim("TC")<<" "<<ii<<"  "<<PrintHit(tjs.fHits[iht])<<" delta "<<delta;
        }
      } // prt
      deltas[ii] = delta;
      if(delta < tp.Delta) {
        tp.Delta = delta;
        imbest = ii;
      }
    } // ii
    
    float chgWght = 0.5;
    
    if(prt) mf::LogVerbatim("TC")<<" nAvailable "<<nAvailable<<" imbest "<<imbest<<" single hit. tp.Delta "<<tp.Delta<<" bestDelta "<<bestDelta<<" path length "<<1 / pathInv<<" imBadRecoHit "<<imBadRecoHit;
    if(imbest == USHRT_MAX || nAvailable == 0) return;
    unsigned int bestDeltaHit = tp.Hits[imbest];
    if(tp.AngleCode == 1) {
      // Get the hits that are in the same multiplet as bestDeltaHit
      std::vector<unsigned int> hitsInMultiplet;
      unsigned short localIndex;
      GetHitMultiplet(bestDeltaHit, hitsInMultiplet, localIndex);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" bestDeltaHit "<<PrintHit(tjs.fHits[bestDeltaHit]);
        myprt<<" in multiplet:";
        for(auto& iht : hitsInMultiplet) myprt<<" "<<PrintHit(tjs.fHits[iht]);
      }
      // Consider the case where a bad reco hit might be better. It is probably wider and
      // has more charge
      if(imBadRecoHit != USHRT_MAX) {
        unsigned int iht = tp.Hits[imBadRecoHit];
        if(tjs.fHits[iht].RMS > HitsRMSTick(tjs, hitsInMultiplet, kUnusedHits)) {
          if(prt) mf::LogVerbatim("TC")<<" Using imBadRecoHit "<<PrintHit(tjs.fHits[iht]);
          tp.UseHit[imBadRecoHit] = true;
          tjs.fHits[iht].InTraj = tj.ID;
          return;
        }
      } // bad reco hit
      // Use the hits in the muliplet instead
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(tjs.fHits[iht].InTraj > 0) continue;
        if(std::find(hitsInMultiplet.begin(), hitsInMultiplet.end(), iht) == hitsInMultiplet.end()) continue;
        tp.UseHit[ii] = true;
        tjs.fHits[iht].InTraj = tj.ID;
      } // ii
      return;
    } // isLA

    // don't use the best UNUSED hit if the best delta is for a USED hit and it is much better
    // TY: ignore for PevProg
    if(bestDelta < 0.5 * tp.Delta && !tj.AlgMod[kRevProp]) return;
    
    if(!useChg || (useChg && (tj.AveChg <= 0 || tj.ChgRMS <= 0))) {
      // necessary quantities aren't available for more careful checking
      if(prt) mf::LogVerbatim("TC")<<" tj.AveChg "<<tj.AveChg<<" or tj.ChgRMS "<<tj.ChgRMS<<". Use the best hit";
      tp.UseHit[imbest] = true;
      tjs.fHits[bestDeltaHit].InTraj = tj.ID;
      return;
    }
    
    // Don't try to get fancy if we are tracking a long muon
    if(tj.PDGCode == 13 && bestDelta < 0.5) {
      if(prt) mf::LogVerbatim("TC")<<" Tracking muon. Use the best hit";
      tp.UseHit[imbest] = true;
      tjs.fHits[bestDeltaHit].InTraj = tj.ID;
      return;
    }
    
    // The best hit is the only one available or this is a small angle trajectory
    if(nAvailable == 1 || tp.AngleCode == 0) {
      float bestDeltaHitChgPull = std::abs(tjs.fHits[bestDeltaHit].Integral * pathInv / tp.AveChg - 1) / tj.ChgRMS;
      if(prt) mf::LogVerbatim("TC")<<" bestDeltaHitChgPull "<<bestDeltaHitChgPull<<" chgPullCut "<<chgPullCut;
      if(bestDeltaHitChgPull < chgPullCut || tp.Delta < 0.1) {
        tp.UseHit[imbest] = true;
        tjs.fHits[bestDeltaHit].InTraj = tj.ID;
      } // good charge or very good delta
      return;
    } // bestDeltaHitMultiplicity == 1
    
    // Find the expected width for the angle of this TP (ticks)
    float expectedWidth = ExpectedHitsRMS(tp);
    
    // Handle two available hits
    if(nAvailable == 2) {
      // See if these two are in the same multiplet and both are available
      std::vector<unsigned int> tHits;
      unsigned short localIndex;
      GetHitMultiplet(bestDeltaHit, tHits, localIndex);
      // ombest is the index of the other hit in tp.Hits that is in the same multiplet as bestDeltaHit
      // if we find it
      unsigned short ombest = USHRT_MAX;
      unsigned int otherHit = INT_MAX;
      if(tHits.size() == 2) {
        otherHit = tHits[1 - localIndex];
        // get the index of this hit in tp.Hits
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tjs.fHits[tp.Hits[ii]].InTraj > 0) continue;
          if(tp.Hits[ii] == otherHit) {
            ombest = ii;
            break;
          }
        } // ii
      } // tHits.size() == 2
      if(prt) {
        mf::LogVerbatim("TC")<<" Doublet: imbest "<<imbest<<" ombest "<<ombest;
      }
      // The other hit exists in the tp and it is available
      if(ombest < tp.Hits.size()) {
        // compare the best delta hit and the other hit separately and the doublet as a merged pair
        float bestHitDeltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(bestDeltaHit);
        // Construct a FOM starting with the delta pull
        float bestDeltaHitFOM = deltas[imbest] /  bestHitDeltaErr;
        if(bestDeltaHitFOM < 0.5) bestDeltaHitFOM = 0.5;
        // multiply by the charge pull if it is significant
        float bestDeltaHitChgPull = std::abs(tjs.fHits[bestDeltaHit].Integral * pathInv / tp.AveChg - 1) / tj.ChgRMS;
        if(bestDeltaHitChgPull > 1) bestDeltaHitFOM *= chgWght * bestDeltaHitChgPull;
        // scale by the ratio
        float rmsRat = tjs.fHits[bestDeltaHit].RMS / expectedWidth;
        if(rmsRat < 1) rmsRat = 1 / rmsRat;
        bestDeltaHitFOM *= rmsRat;
        if(prt) mf::LogVerbatim("TC")<<" bestDeltaHit FOM "<<deltas[imbest]/bestHitDeltaErr<<" bestDeltaHitChgPull "<<bestDeltaHitChgPull<<" rmsRat "<<rmsRat<<" bestDeltaHitFOM "<<bestDeltaHitFOM;
        // Now do the same for the other hit
        float otherHitDeltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(otherHit);
        float otherHitFOM = deltas[ombest] /  otherHitDeltaErr;
        if(otherHitFOM < 0.5) otherHitFOM = 0.5;
        float otherHitChgPull = std::abs(tjs.fHits[otherHit].Integral * pathInv / tp.AveChg - 1) / tj.ChgRMS;
        if(otherHitChgPull > 1) otherHitFOM *= chgWght * otherHitChgPull;
        rmsRat = tjs.fHits[otherHit].RMS / expectedWidth;
        if(rmsRat < 1) rmsRat = 1 / rmsRat;
        otherHitFOM *= rmsRat;
        if(prt) mf::LogVerbatim("TC")<<" otherHit FOM "<<deltas[ombest]/otherHitDeltaErr<<" otherHitChgPull "<<otherHitChgPull<<" rmsRat "<<rmsRat<<" otherHitFOM "<<otherHitFOM;
        // And for the doublet
        float doubletChg = tjs.fHits[bestDeltaHit].Integral + tjs.fHits[otherHit].Integral;
        float doubletTime = (tjs.fHits[bestDeltaHit].Integral * tjs.fHits[bestDeltaHit].PeakTime + tjs.fHits[otherHit].Integral * tjs.fHits[otherHit].PeakTime) / doubletChg;
        doubletChg *= pathInv;
        doubletTime *= tjs.UnitsPerTick;
        float doubletWidthTick = TPHitsRMSTick(tjs, tp, kUnusedHits);
        float doubletRMSTimeErr = doubletWidthTick * tjs.UnitsPerTick;
        if(prt) mf::LogVerbatim("TC")<<" doublet Chg "<<doubletChg<<" doubletTime "<<doubletTime<<" doubletRMSTimeErr "<<doubletRMSTimeErr;
        float doubletFOM = PointTrajDOCA(tjs, tp.Pos[0], doubletTime, tp) / doubletRMSTimeErr;
        if(doubletFOM < 0.5) doubletFOM = 0.5;
        float doubletChgPull = std::abs(doubletChg * pathInv / tp.AveChg - 1) / tj.ChgRMS;
        if(doubletChgPull > 1) doubletFOM *= chgWght * doubletChgPull;
        rmsRat = doubletWidthTick / expectedWidth;
        if(rmsRat < 1) rmsRat = 1 / rmsRat;
        doubletFOM *= rmsRat;
        if(prt) mf::LogVerbatim("TC")<<" doublet FOM "<<PointTrajDOCA(tjs, tp.Pos[0], doubletTime, tp)/doubletRMSTimeErr<<" doubletChgPull "<<doubletChgPull<<" rmsRat "<<rmsRat<<" doubletFOM "<<doubletFOM;
        // Assume the doublet is best
        if(tjs.fHits[bestDeltaHit].InTraj > 0 || tjs.fHits[otherHit].InTraj > 0) {
          std::cout<<"Crap \n";
          exit(1);
        }
        if(doubletFOM < bestDeltaHitFOM && doubletFOM < otherHitFOM) {
          tp.UseHit[imbest] = true;
          tjs.fHits[bestDeltaHit].InTraj = tj.ID;
          tp.UseHit[ombest] = true;
          tjs.fHits[otherHit].InTraj = tj.ID;
        } else {
          // the doublet is not the best
          if(bestDeltaHitFOM < otherHitFOM) {
            tp.UseHit[imbest] = true;
            tjs.fHits[bestDeltaHit].InTraj = tj.ID;
          } else {
            tp.UseHit[ombest] = true;
            tjs.fHits[otherHit].InTraj = tj.ID;
          } // otherHit is the best
        } // doublet is not the best
      } else {
        // the other hit isn't available. Just use the singlet
        tp.UseHit[imbest] = true;
        tjs.fHits[bestDeltaHit].InTraj = tj.ID;
      }
      return;
    } // nAvailable == 2
    
    // we are left with nAvailable > 2 

    // Use all of the hits if they are all available and are in the same multiplet
    if(nAvailable == tp.Hits.size()) {
      // the first hit in the multiplet
      unsigned int hit0 = tp.Hits[0] - tjs.fHits[tp.Hits[0]].LocalIndex;
      unsigned short cnt = 0;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(iht - tjs.fHits[iht].LocalIndex == hit0) ++cnt;
      } // ii
      // all in the same multiplet
      if(cnt == tp.Hits.size()) {
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          unsigned int iht = tp.Hits[ii];
          if(tjs.fHits[iht].InTraj > 0) continue;
          tp.UseHit[ii] = true;
          tjs.fHits[iht].InTraj = tj.ID;
        } // ii
        return;
      } // all in the same multiplet
    } // nAvailable == tp.Hits.size()

    float hitsWidth = TPHitsRMSTick(tjs, tp, kUnusedHits);
    float maxTick = tp.Pos[1] / tjs.UnitsPerTick + 0.6 * expectedWidth;
    float minTick = tp.Pos[1] / tjs.UnitsPerTick - 0.6 * expectedWidth;
    if(prt) mf::LogVerbatim("TC")<<" Multiplet: hitsWidth "<<hitsWidth<<" expectedWidth "<<expectedWidth<<" tick range "<<(int)minTick<<" "<<(int)maxTick;
    // use all of the hits in the tick window
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      unsigned int iht = tp.Hits[ii];
      if(tjs.fHits[iht].InTraj > 0) continue;
      if(tjs.fHits[iht].PeakTime < minTick) continue;
      if(tjs.fHits[iht].PeakTime > maxTick) continue;
      tp.UseHit[ii] = true;
      tjs.fHits[iht].InTraj = tj.ID;
    }

  } // FindUseHits
  
  //////////////////////////////////////////
  void TrajClusterAlg::DefineHitPos(TrajPoint& tp)
  {
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    
    tp.Chg = 0;
    if(tp.Hits.empty()) return;
    
    unsigned short nused = 0;
    unsigned int iht = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(tp.UseHit[ii]) {
        iht = tp.Hits[ii];
        ++nused;
      }
    }
    if(nused == 0) return;
    
    // don't bother with rest of this if there is only one hit since it can
    // only reside on one wire
    if(nused == 1) {
      tp.Chg = tjs.fHits[iht].Integral;
      // Normalize to 1 WSE path length
      float pathInv = std::abs(tp.Dir[0]);
      if(pathInv < 0.05) pathInv = 0.05;
      tp.Chg *= pathInv;
      tp.HitPos[0] = tjs.fHits[iht].WireID.Wire;
      tp.HitPos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
      float wireErr = tp.Dir[1] * 0.289;
      float timeErr = tp.Dir[0] * HitTimeErr(iht);
      tp.HitPosErr2 = wireErr * wireErr + timeErr * timeErr;
      if(prt) mf::LogVerbatim("TC")<<"DefineHitPos: singlet "<<std::fixed<<std::setprecision(1)<<tp.HitPos[0]<<":"<<(int)(tp.HitPos[1]/tjs.UnitsPerTick)<<" ticks. HitPosErr "<<sqrt(tp.HitPosErr2);
      return;
    } // nused == 1
    
    // multiple hits possibly on different wires
    std::vector<unsigned int> hitVec;
    tp.Chg = 0;
    std::array<float, 2> newpos;
    float chg;
    newpos[0] = 0;
    newpos[1] = 0;
    // Find the wire range for hits used in the TP
    unsigned int loWire = INT_MAX;
    unsigned int hiWire = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      unsigned int iht = tp.Hits[ii];
      chg = tjs.fHits[iht].Integral;
      unsigned int wire = tjs.fHits[iht].WireID.Wire;
      if(wire < loWire) loWire = wire;
      if(wire > hiWire) hiWire = wire;
      newpos[0] += chg * wire;
      newpos[1] += chg * tjs.fHits[iht].PeakTime;
      tp.Chg += chg;
      hitVec.push_back(iht);
    } // ii
 
    if(tp.Chg == 0) return;
    
    tp.HitPos[0] = newpos[0] / tp.Chg;
    tp.HitPos[1] = newpos[1] * tjs.UnitsPerTick / tp.Chg;
    
    // Normalize to 1 WSE path length
    float pathInv = std::abs(tp.Dir[0]);
    if(pathInv < 0.05) pathInv = 0.05;
    tp.Chg *= pathInv;
    
    // Error is the wire error (1/sqrt(12))^2 if all hits are on one wire.
    // Scale it by the wire range
    float dWire = 1 + hiWire - loWire;
    float wireErr = tp.Dir[1] * dWire * 0.289;
    float timeErr2 = tp.Dir[0] * tp.Dir[0] * HitsTimeErr2(hitVec);
    tp.HitPosErr2 = wireErr * wireErr + timeErr2;
    if(prt) mf::LogVerbatim("TC")<<"DefineHitPos: multiplet "<<std::fixed<<std::setprecision(1)<<tp.HitPos[0]<<":"<<(int)(tp.HitPos[1]/tjs.UnitsPerTick)<<" ticks. HitPosErr "<<sqrt(tp.HitPosErr2);

  } // HitPosErr2

  //////////////////////////////////////////
  float TrajClusterAlg::HitTimeErr(const unsigned int iht)
  {
    return tjs.fHits[iht].RMS * tjs.UnitsPerTick * fHitErrFac * tjs.fHits[iht].Multiplicity;
  } // HitTimeErr
  
  //////////////////////////////////////////
  float TrajClusterAlg::HitsTimeErr2(const std::vector<unsigned int>& hitVec)
  {
    // Estimates the error^2 of the time using all hits in hitVec
    if(hitVec.empty()) return 0;
    float err = fHitErrFac * HitsRMSTime(tjs, hitVec, kUnusedHits);
    return err * err;
  } // HitsTimeErr2

  ////////////////////////////////////////////////
  void TrajClusterAlg::EndMerge()
  {
    // Merges trajectories end-to-end or makes vertices
    
    if(tjs.allTraj.size() < 2) return;
    if(!tjs.UseAlg[kEndMerge]) return;
    
    mrgPrt = (debug.Plane == (int)fPlane && debug.Wire < 0);
    if(mrgPrt) mf::LogVerbatim("TC")<<"inside EndMerge2 on plane "<<fPlane;
    
    // Ensure that all tjs are in the same order
    bool first = true;
    unsigned short firstTj = 0;
    for(unsigned short it1 = 1; it1 < tjs.allTraj.size() - 1; ++it1) {
      Trajectory& tj = tjs.allTraj[it1];
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != fCTP) continue;
      if(first) {
        firstTj = it1;
        first = false;
      }
      if(tj.StepDir != tjs.allTraj[firstTj].StepDir) {
        mf::LogVerbatim("TC")<<"EndMerge2: Trajectory "<<tj.ID<<" is out of order. Fixing it.";
        ReverseTraj(tjs, tj);
      }
    } // it1
    
    unsigned short maxShortTjLen = tjs.Vertex2DCuts[0];
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size() - 1; ++it1) {
      if(tjs.allTraj[it1].AlgMod[kKilled]) continue;
      if(tjs.allTraj[it1].CTP != fCTP) continue;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // no merge if there is a vertex at the end
        if(tjs.allTraj[it1].VtxID[end1] > 0) continue;
        // make a copy of tp1 so we can mess with it
        TrajPoint tp1 = tjs.allTraj[it1].Pts[tjs.allTraj[it1].EndPt[end1]];
        bool isVLA = (tp1.AngleCode == 2);
        float bestFOM = 5;
        if(isVLA) bestFOM = 20;
        float bestDOCA;
        unsigned short imbest = USHRT_MAX;
        for(unsigned short it2 = 0; it2 < tjs.allTraj.size(); ++it2) {
          if(it1 == it2) continue;
          if(tjs.allTraj[it2].AlgMod[kKilled]) continue;
          if(tjs.allTraj[it2].CTP != fCTP) continue;
          // ensure that MCSMom isn't wildly different if they are both long
          if(tjs.allTraj[it1].Pts.size() > 8 && tjs.allTraj[it2].Pts.size() > 8) {
            if(tjs.allTraj[it2].MCSMom > tjs.allTraj[it1].MCSMom) {
              if(tjs.allTraj[it1].MCSMom < 0.5 * tjs.allTraj[it2].MCSMom) continue;
            } else {
              if(tjs.allTraj[it2].MCSMom < 0.5 * tjs.allTraj[it1].MCSMom) continue;
            }
          }
          unsigned short end2 = 1 - end1;
          // check for a vertex at this end
          if(tjs.allTraj[it2].VtxID[end2] > 0) continue;
          TrajPoint& tp2 = tjs.allTraj[it2].Pts[tjs.allTraj[it2].EndPt[end2]];
          TrajPoint& tp2OtherEnd = tjs.allTraj[it2].Pts[tjs.allTraj[it2].EndPt[end1]];
          // ensure that the other end isn't closer
          if(std::abs(tp2OtherEnd.Pos[0] - tp1.Pos[0]) < std::abs(tp2.Pos[0] - tp1.Pos[0])) continue;
          // ensure that the order is correct
          if(tjs.allTraj[it1].StepDir > 0) {
            if(tp2.Pos[0] < tp1.Pos[0] - 2) continue;
          } else {
            if(tp2.Pos[0] > tp1.Pos[0] + 2) continue;
          }
          // ensure that there is a signal on most of the wires between these points
          if(!SignalBetween(tjs, tp1, tp2, 0.8, mrgPrt)) {
//            if(mrgPrt) mf::LogVerbatim("TC")<<" no signal between these points "<<PrintPos(tjs, tp1.Pos)<<" "<<PrintPos(tjs, tp2.Pos);
            continue;
          }
          // Find the distance of closest approach for small angle merging
          // Inflate the doca cut if we are bridging a block of dead wires
          float dang = DeltaAngle(tp1.Ang, tp2.Ang);
          float doca = 5;
          if(isVLA) {
            // compare the minimum separation between Large Angle trajectories using a generous cut
            unsigned short ipt1, ipt2;
            TrajTrajDOCA(tjs, tjs.allTraj[it1], tjs.allTraj[it2], ipt1, ipt2, doca);
            if(mrgPrt) mf::LogVerbatim("TC")<<" isVLA check ipt1 "<<ipt1<<" ipt2 "<<ipt2<<" doca "<<doca;
          } else {
            // small angle
            doca = PointTrajDOCA(tjs, tp1.Pos[0], tp1.Pos[1], tp2);
          }
          float fom = dang * doca;
          if(fom < bestFOM) {
            bestFOM = fom;
            bestDOCA = doca;
            imbest = it2;
          }
        } // it2
        // No merge/vertex candidates
        if(imbest == USHRT_MAX) continue;
        
        // Make angle adjustments to tp1.
        unsigned short it2 = imbest;
        unsigned short end2 = 1 - end1;
        bool loMCSMom = (tjs.allTraj[it1].MCSMom + tjs.allTraj[it2].MCSMom) < 150;
        // Don't use the angle at the end Pt for high momentum long trajectories in case there is a little kink at the end
        if(tjs.allTraj[it1].Pts.size() > 50 && tjs.allTraj[it1].MCSMom > 100) {
          if(end1 == 0) {
            tp1.Ang = tjs.allTraj[it1].Pts[tjs.allTraj[it1].EndPt[0] + 2].Ang;
          } else {
            tp1.Ang = tjs.allTraj[it1].Pts[tjs.allTraj[it1].EndPt[1] - 2].Ang;
          }
        } else if(loMCSMom) {
          // Low momentum - calculate the angle using the two Pts at the end
          unsigned short pt1, pt2;
          if(end1 == 0) {
            pt1 = tjs.allTraj[it1].EndPt[0];
            pt2 = pt1 + 1;
          } else {
            pt2 = tjs.allTraj[it1].EndPt[1];
            pt1 = pt2 - 1;
          }
          TrajPoint tpdir;
          if(MakeBareTrajPoint(tjs, tjs.allTraj[it1].Pts[pt1], tjs.allTraj[it1].Pts[pt2], tpdir)) tp1.Ang = tpdir.Ang;
        } // low MCSMom
        // Now do the same for tj2
        TrajPoint tp2 = tjs.allTraj[it2].Pts[tjs.allTraj[it2].EndPt[end2]];
        if(tjs.allTraj[it2].Pts.size() > 50 && tjs.allTraj[it2].MCSMom > 100) {
          if(end1 == 0) {
            tp2.Ang = tjs.allTraj[it2].Pts[tjs.allTraj[it2].EndPt[0] + 2].Ang;
          } else {
            tp2.Ang = tjs.allTraj[it2].Pts[tjs.allTraj[it2].EndPt[1] - 2].Ang;
          }
        } else if(loMCSMom) {
          // Low momentum - calculate the angle using the two Pts at the end
          unsigned short pt1, pt2;
          if(end2 == 0) {
            pt1 = tjs.allTraj[it2].EndPt[0];
            pt2 = pt1 + 1;
          } else {
            pt2 = tjs.allTraj[it2].EndPt[1];
            pt1 = pt2 - 1;
          }
          TrajPoint tpdir;
          if(MakeBareTrajPoint(tjs, tjs.allTraj[it2].Pts[pt1], tjs.allTraj[it2].Pts[pt2], tpdir)) tp2.Ang = tpdir.Ang;
        } // low MCSMom
        
        // decide whether to merge or make a vertex
        float dang = DeltaAngle(tp1.Ang, tp2.Ang);

        float dangCut;
        float docaCut;
        float chgPull = 0;
        float minChgRMS = tjs.allTraj[it1].ChgRMS;
        if(tjs.allTraj[it2].ChgRMS < minChgRMS) minChgRMS = tjs.allTraj[it2].ChgRMS;
        if(tp1.Chg > tp2.Chg) {
          chgPull = (tp1.Chg / tp2.Chg - 1) / minChgRMS;
        } else {
          chgPull = (tp2.Chg / tp1.Chg - 1) / minChgRMS;
        }
        if(loMCSMom) {
          // increase dangCut dramatically for low MCSMom tjs
          dangCut = 1.0;
          // and the doca cut
          docaCut = 2;
        } else {
          // do a more careful calculation of the angle cut
          unsigned short e0 = tjs.allTraj[it1].EndPt[0];
          unsigned short e1 = tjs.allTraj[it1].EndPt[1];
          float tj1len = TrajPointSeparation(tjs.allTraj[it1].Pts[e0], tjs.allTraj[it1].Pts[e1]);
          float thetaRMS1 = MCSThetaRMS(tjs, tjs.allTraj[it1]);
          // calculate (thetaRMS / sqrt(length) )^2
          thetaRMS1 *= thetaRMS1 / tj1len;
          // and now tj2
          e0 = tjs.allTraj[it2].EndPt[0];
          e1 = tjs.allTraj[it2].EndPt[1];
          float tj2len = TrajPointSeparation(tjs.allTraj[it2].Pts[e0], tjs.allTraj[it2].Pts[e1]);
          float thetaRMS2 = MCSThetaRMS(tjs, tjs.allTraj[it2]);
          thetaRMS2 *= thetaRMS2 / tj2len;
          float dangErr = 0.5 * sqrt(thetaRMS1 + thetaRMS2);
          dangCut = fKinkCuts[0] + fKinkCuts[1] * dangErr;
          docaCut = 1;
          if(isVLA) docaCut = 15;
        }
        
        // check the merge cuts. Start with doca and dang requirements
        bool doMerge = bestDOCA < docaCut && dang < dangCut;
        bool showerTjs = tjs.allTraj[it1].PDGCode == 11 || tjs.allTraj[it2].PDGCode == 11;
        bool hiMCSMom = tjs.allTraj[it1].MCSMom > 200 || tjs.allTraj[it2].MCSMom > 200;
        // add a charge similarity requirement if not shower-like or low momentum or not LA
        if(doMerge && !showerTjs && hiMCSMom && chgPull > fChargeCuts[0] && !isVLA) doMerge = false;
        // ignore the charge pull cut if both are high momentum and dang is really small
        if(!doMerge && tjs.allTraj[it1].MCSMom > 900 && tjs.allTraj[it2].MCSMom > 900 && dang < 0.1 && bestDOCA < docaCut) doMerge = true;

        // do not merge if chgPull is really high
        if(doMerge && chgPull > 2*fChargeCuts[0]) doMerge = false;

        bool signalBetween = true;
        if(!isVLA) signalBetween = SignalBetween(tjs, tp1, tp2, 0.99, mrgPrt);
        doMerge = doMerge && signalBetween;
        
        if(mrgPrt) mf::LogVerbatim("TC")<<" EM2 "<<tjs.allTraj[it1].ID<<"_"<<end1<<"-"<<tjs.allTraj[it2].ID<<"_"<<end2<<" tp1-tp2 "<<PrintPos(tjs, tp1)<<"-"<<PrintPos(tjs, tp2)<<" bestFOM "<<bestFOM<<" bestDOCA "<<bestDOCA<<" docaCut "<<docaCut<<" isVLA? "<<isVLA<<" dang "<<dang<<" dangCut "<<dangCut<<" chgPull "<<chgPull<<" doMerge? "<<doMerge<<" loMCSMom? "<<loMCSMom<<" hiMCSMom? "<<hiMCSMom<<" signalBetween? "<<signalBetween;

        if(doMerge) {
          if(mrgPrt) mf::LogVerbatim("TC")<<"  Merge ";
          bool didMerge = false;
          if(end1 == 1) {
            didMerge = MergeAndStore(it1, it2, mrgPrt);
          } else {
            didMerge = MergeAndStore(it2, it1, mrgPrt);
          }
          if(didMerge) {
            // wipe out the AlgMods for the new Trajectory
            unsigned short newTjIndex = tjs.allTraj.size()-1;
            tjs.allTraj[newTjIndex].AlgMod.reset();
            // and set the EndMerge bit
            tjs.allTraj[newTjIndex].AlgMod[kEndMerge] = true;
            // and maybe the RevProp bit
            if(tjs.allTraj[it1].AlgMod[kRevProp] || tjs.allTraj[it2].AlgMod[kRevProp]) tjs.allTraj[newTjIndex].AlgMod[kRevProp] = true;
            // Set the end merge flag for the killed trajectories to aid tracing merges
            tjs.allTraj[it1].AlgMod[kEndMerge] = true;
            tjs.allTraj[it2].AlgMod[kEndMerge] = true;
          } // Merge and store successfull
          else {
            if(mrgPrt) mf::LogVerbatim("TC")<<"  MergeAndStore failed ";
          }
        } else {
          // create a vertex instead if it passes the vertex cuts
          VtxStore aVtx;
          aVtx.CTP = tjs.allTraj[it1].CTP;
          aVtx.ID = tjs.vtx.size() + 1;
          // keep it simple if tp1 and tp2 are very close
          if(std::abs(tp1.Pos[0] - tp2.Pos[0]) < 2 && std::abs(tp1.Pos[1] - tp2.Pos[1]) < 2) {
            aVtx.Pos[0] = 0.5 * (tp1.Pos[0] + tp2.Pos[0]);
            aVtx.Pos[1] = 0.5 * (tp1.Pos[1] + tp2.Pos[1]);
          } else {
            // Tps not so close
            float sepCut = tjs.Vertex2DCuts[2];
            bool tj1Short = (tjs.allTraj[it1].EndPt[1] - tjs.allTraj[it1].EndPt[0] < maxShortTjLen);
            bool tj2Short = (tjs.allTraj[it2].EndPt[1] - tjs.allTraj[it2].EndPt[0] < maxShortTjLen);
            if(tj1Short || tj2Short) sepCut = tjs.Vertex2DCuts[1];
            TrajIntersection(tp1, tp2, aVtx.Pos);
            float dw = aVtx.Pos[0] - tp1.Pos[0];
            if(std::abs(dw) > sepCut) continue;
            float dt = aVtx.Pos[1] - tp1.Pos[1];
            if(std::abs(dt) > sepCut) continue;
            dw = aVtx.Pos[0] - tp2.Pos[0];
            if(std::abs(dw) > sepCut) continue;
            dt = aVtx.Pos[1] - tp2.Pos[1];
            if(std::abs(dt) > sepCut) continue;
            // ensure that the vertex is not closer to the other end if the tj is short
            if(tj1Short) {
              TrajPoint otp1 = tjs.allTraj[it1].Pts[tjs.allTraj[it1].EndPt[1-end1]];
              if(PosSep2(otp1.Pos, aVtx.Pos) < PosSep2(tp1.Pos, aVtx.Pos)) continue;
            }
            if(tj2Short) {
              TrajPoint otp2 = tjs.allTraj[it2].Pts[tjs.allTraj[it2].EndPt[1-end2]];
              if(PosSep2(otp2.Pos, aVtx.Pos) < PosSep2(tp2.Pos, aVtx.Pos)) continue;
            }
            // we expect the vertex to be between tp1 and tp2
            if(aVtx.Pos[0] < tp1.Pos[0] && aVtx.Pos[0] < tp2.Pos[0]) {
              aVtx.Pos[0] = std::min(tp1.Pos[0], tp2.Pos[0]);
              aVtx.Stat[kFixed] = true;
            }
            if(aVtx.Pos[0] > tp1.Pos[0] && aVtx.Pos[0] > tp2.Pos[0]) {
              aVtx.Pos[0] = std::max(tp1.Pos[0], tp2.Pos[0]);
              aVtx.Stat[kFixed] = true;
            }
          } // Tps not so close
          // We got this far. Try a vertex fit to ensure that the errors are reasonable
          tjs.allTraj[it1].VtxID[end1] = aVtx.ID;
          tjs.allTraj[it2].VtxID[end2] = aVtx.ID;
          // do a fit
          if(!FitVertex(tjs, aVtx, mrgPrt)) {
            // back out
            tjs.allTraj[it1].VtxID[end1] = 0;
            tjs.allTraj[it2].VtxID[end2] = 0;
            if(mrgPrt) mf::LogVerbatim("TC")<<" Vertex fit failed ";
            continue;
          }
          aVtx.NTraj = 2;
          aVtx.Pass = tjs.allTraj[it1].Pass;
          aVtx.Topo = end1 + end2;
          tjs.allTraj[it1].AlgMod[kEndMerge] = true;
          tjs.allTraj[it2].AlgMod[kEndMerge] = true;
          if(mrgPrt) mf::LogVerbatim("TC")<<"  Make vertex ID "<<aVtx.ID<<" at "<<(int)aVtx.Pos[0]<<":"<<(int)(aVtx.Pos[1]/tjs.UnitsPerTick);
          if(!StoreVertex(tjs, aVtx)) continue;
        } // create a vertex
        if(tjs.allTraj[it1].AlgMod[kKilled]) break;
      } // end1
    } // it1
    
  } // EndMerge
  
  //////////////////////////////////////////
  void TrajClusterAlg::Match3D(const geo::TPCID& tpcid)
  {
    // Match trajectories in each view using hits. This function loops through hits that are used in trajectories
    // in all (3) planes in the current tpcid. The triplet of Tj IDs is stored in a MatchStruct, ms, and pushed onto
    // tjs.matchVec. A count is made of all subsequent hit triplets that have the same Tj ID triplet.  The matchVec vector
    // is sorted by the decreasing number of hit triplet counts. The matchVec entries with the most hit triplet counts are
    // identified as to-be-made PFParticles and put in matchVecPFPList. The function Find3DEndPoints defines the sXYZ and eXYZ
    // elements of the Match Struct using a convention that sXYZ[0] is at low X end of the matched trajectories.
    
    if(!tjs.UseAlg[kMatch3D]) return;
    if(tjs.Match3DCuts[0] < 0) return;
    
    prt = (debug.Plane >= 0) && (debug.Tick == 3333);
    
    if(prt) mf::LogVerbatim("TC")<<"inside Match3D. dX (cm) cut "<<tjs.Match3DCuts[0];
    
    // vector of X positions of all hits
    std::vector<float> xx(tjs.fHits.size());
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      TCHit& hit = tjs.fHits[iht];
      xx[iht] = tjs.detprop->ConvertTicksToX(hit.PeakTime, hit.WireID.Plane, hit.WireID.TPC, hit.WireID.Cryostat);
    } // iht
    
    if(tjs.NumPlanes == 3) {
      Match3Views(tpcid, xx);
      // Now try to find 3D matches in 2 views using the unmatched hits
      Match2Views(tpcid, xx);
    } else {
      Match2Views(tpcid, xx);
    }
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"matchVec\n";
      for(unsigned int ii = 0; ii < tjs.matchVec.size(); ++ii) {
        myprt<<ii<<" Count "<<tjs.matchVec[ii].Count<<" TjIDs: ";
        unsigned short nstj = 0;
        for(auto& tjID : tjs.matchVec[ii].TjIDs) {
          myprt<<" "<<tjID;
          // count shower Tjs
          if(tjs.allTraj[tjID - 1].AlgMod[kShowerTj]) ++nstj;
        }
        myprt<<" nShowerTj "<<nstj<<"\n";
      } // ii
/*
      myprt<<"matchVecPFPList";
      for(auto& pfp : tjs.matchVecPFPList) {
        myprt<<" "<<pfp<<" TjIDs";
        for(auto& tjID : tjs.matchVec[pfp].TjIDs) myprt<<" "<<tjID;
      }
*/
    } // prt
    
    // ensure that a Tj will only used by one PFParticle
    std::vector<unsigned short> tjInPFParticle(tjs.allTraj.size(), USHRT_MAX);
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned short imv = tjs.matchVecPFPList[ipfp];
      auto& ims = tjs.matchVec[imv];
      for(unsigned short ii = 0; ii < ims.TjIDs.size(); ++ii) {
        unsigned short itj = ims.TjIDs[ii] - 1;
        if(tjInPFParticle[itj] != USHRT_MAX) {
          std::cout<<"Error! itj "<<itj<<" tj.ID "<<tjs.allTraj[itj].ID<<" TjIDs["<<ii<<"] = "<<ims.TjIDs[ii]<<" already set.";
          std::cout<<" Trying to set ipfp "<<ipfp<<" imv "<<imv<<"\n";
        }
        tjInPFParticle[itj] = ipfp;
      } // ii
    } // ipfp
    
    // Define the parent (j) - daughter (i) relationship and the PDGCode
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned short imv = tjs.matchVecPFPList[ipfp];
      auto& ims = tjs.matchVec[imv];
      // assume that this is its own parent
      ims.ParentMSIndex = ipfp;
      unsigned short n11 = 0;
      unsigned short n13 = 0;
      unsigned short nsh = 0;
      // look for a parent (j) in the list of trajectories
      for(unsigned short ii = 0; ii < ims.TjIDs.size(); ++ii) {
        unsigned short itj = ims.TjIDs[ii] - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.PDGCode == 11) ++n11;
        if(tj.PDGCode == 13) ++n13;
        if(tj.AlgMod[kShowerTj]) ++nsh;
        // Look for a parent trajectory that is matched in 3D
        if(tj.ParentTrajID > 0 && tjs.allTraj[tj.ParentTrajID - 1].AlgMod[kMatch3D]) {
          // Look for this parent in matchVecPFPList
          for(unsigned short jpfp = 0; jpfp < tjs.matchVecPFPList.size(); ++jpfp) {
            if(jpfp == ipfp) continue;
            unsigned short jmv = tjs.matchVecPFPList[jpfp];
            auto& jms = tjs.matchVec[jmv];
            if(std::find(jms.TjIDs.begin(), jms.TjIDs.end(), tj.ParentTrajID) != jms.TjIDs.end()) {
              ims.ParentMSIndex = jpfp;
              if(std::find(jms.DtrIndices.begin(), jms.DtrIndices.end(), ipfp) == jms.DtrIndices.end()) jms.DtrIndices.push_back(ipfp);
              break;
            }
          } // jpfp
        } // ParentTrajID > 0
      } // ii
      ims.PDGCode = 13;
      if(nsh > 1) {
        // use PDGCode = 1111 for a shower Tj
        ims.PDGCode = 1111;
      } else if(n11 > n13) {
        // a generic shower-like Tj
        ims.PDGCode = 11;
      }
      if(prt) mf::LogVerbatim("TC")<<"ipfp "<<ipfp<<" n11 "<<n11<<" n13 "<<n13<<" nsh "<<nsh<<" PDGCode "<<ims.PDGCode<<" Parent "<<ims.ParentMSIndex;
    } // ipfp
    
    Find3DEndPoints(tpcid);
    Find3DShowerEndPoints(tjs, tpcid);

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"matchVec out\n";
      for(unsigned int ii = 0; ii < tjs.matchVec.size(); ++ii) {
        myprt<<ii<<" Count "<<tjs.matchVec[ii].Count<<" TjIDs: ";
        for(auto& tjID : tjs.matchVec[ii].TjIDs) myprt<<" "<<tjID;
        myprt<<" PDGCode "<<tjs.matchVec[ii].PDGCode;
        myprt<<"\n";
      } // ii
      myprt<<"matchVecPFPList";
      for(auto& pfp : tjs.matchVecPFPList) myprt<<" "<<pfp;
    } // prt

/*
    // Ensure that there is no more than 1 trajectory in each plane in each ms entry
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned short imv = tjs.matchVecPFPList[ipfp];
      auto& ims = tjs.matchVec[imv];
      std::array<unsigned short, 3> nInPln = {0};
      for(auto& tjID : ims.TjIDs) {
        geo::PlaneID plnID = DecodeCTP(tjs.allTraj[tjID - 1].CTP);
        ++nInPln[plnID.Plane];
      }
      if(nInPln[0] > 1 ||nInPln[1] > 1 ||nInPln[2] > 1) std::cout<<">>>>>>>>>>> Didn't merge broken Tjs\n";
    } // ipfp check
*/
  } // Match3D
  
  //////////////////////////////////////////
  void TrajClusterAlg::Match3Views(const geo::TPCID& tpcid, const std::vector<float>& xx)
  {
    // Match Tjs in 3 views

    
    // decide if we should take the time to check the MCSMom cut
    bool chkMCSMom = (tjs.Match3DCuts[1] > 0);
    short momCut = tjs.Match3DCuts[1];
    
    // Define the order of the planes in which hits will be considered. 
    unsigned short ipl = 2;
    unsigned short jpl = 0;
    unsigned short kpl = 1;
    unsigned int klastWire = tjs.LastWire[kpl] - 1;
    
    for(unsigned int iwire = tjs.FirstWire[ipl]; iwire < tjs.LastWire[ipl]; ++iwire) {
      if(tjs.WireHitRange[ipl][iwire].first < 0) continue;
      for(unsigned int jwire = tjs.FirstWire[jpl]; jwire < tjs.LastWire[jpl]; ++jwire) {
        if(tjs.WireHitRange[jpl][jwire].first < 0) continue;
        double yp, zp;
        tjs.geom->IntersectionPoint(iwire, jwire, ipl, jpl, tpcid.Cryostat, tpcid.TPC, yp, zp);
        // ensure this is inside the TPC
        if(yp < tjs.YLo || yp > tjs.YHi || zp < tjs.ZLo || zp > tjs.ZHi) continue;
        float fkwire = tjs.geom->WireCoordinate(yp, zp, kpl, tpcid.TPC, tpcid.Cryostat);
        if(fkwire < 0 || fkwire > tjs.MaxPos0[kpl]) continue;
        unsigned int kwire = std::nearbyint(fkwire);
        if(kwire > klastWire) continue;
        if(tjs.WireHitRange[kpl][kwire].first < 0) continue;
        // Have intersecting wire matches. Now look for time matches
        unsigned int ifirsthit = (unsigned int)tjs.WireHitRange[ipl][iwire].first;
        unsigned int ilasthit = (unsigned int)tjs.WireHitRange[ipl][iwire].second;
        unsigned int jfirsthit = (unsigned int)tjs.WireHitRange[jpl][jwire].first;
        unsigned int jlasthit = (unsigned int)tjs.WireHitRange[jpl][jwire].second;
        unsigned int kfirsthit = (unsigned int)tjs.WireHitRange[kpl][kwire].first;
        unsigned int klasthit = (unsigned int)tjs.WireHitRange[kpl][kwire].second;
        for(unsigned int iht = ifirsthit; iht < ilasthit; ++iht) {
          // Only consider hits that are associated with a trajectory
          if(tjs.fHits[iht].InTraj <= 0) continue;
          // see if this trajectory already has a 3D match
          unsigned short itj = tjs.fHits[iht].InTraj - 1;
          if(tjs.allTraj[itj].AlgMod[kMatch3D]) continue;
          // don't try to match Tjs in showers
          if(tjs.allTraj[itj].AlgMod[kInShower]) continue;
          for(unsigned int jht = jfirsthit; jht < jlasthit; ++jht) {
            if(tjs.fHits[jht].InTraj <= 0) continue;
            // see if this trajectory already has a 3D match
            unsigned short jtj = tjs.fHits[jht].InTraj - 1;
            if(tjs.allTraj[jtj].AlgMod[kMatch3D]) continue;
            // don't try to match Tjs in showers
            if(tjs.allTraj[jtj].AlgMod[kInShower]) continue;
            if(xx[jht] < xx[iht] - tjs.Match3DCuts[0]) continue;
            if(xx[jht] > xx[iht] + tjs.Match3DCuts[0]) break;
            for(unsigned int kht = kfirsthit; kht < klasthit; ++kht) {
              if(tjs.fHits[kht].InTraj <= 0) continue;
              unsigned short ktj = tjs.fHits[kht].InTraj - 1;
              if(tjs.allTraj[ktj].AlgMod[kMatch3D]) continue;
              if(xx[kht] < xx[iht] - tjs.Match3DCuts[0]) continue;
              if(xx[kht] > xx[iht] + tjs.Match3DCuts[0]) break;
              // we have a time match
              // check MCSMom?
              if(chkMCSMom) {
                unsigned short npass = 0;
                unsigned short itj = tjs.fHits[iht].InTraj - 1;
                if(tjs.allTraj[itj].MCSMom > momCut) ++npass;
                itj = tjs.fHits[jht].InTraj - 1;
                if(tjs.allTraj[itj].MCSMom > momCut) ++npass;
                itj = tjs.fHits[kht].InTraj - 1;
                if(tjs.allTraj[itj].MCSMom > momCut) ++npass;
                // require at least two tjs to pass the cut
                if(npass < 2) continue;
              }
              // next see if the Tj IDs are in the match list
              unsigned short indx = 0;
              for(indx = 0; indx < tjs.matchVec.size(); ++indx) {
                if(tjs.fHits[iht].InTraj != tjs.matchVec[indx].TjIDs[ipl]) continue;
                if(tjs.fHits[jht].InTraj != tjs.matchVec[indx].TjIDs[jpl]) continue;
                if(tjs.fHits[kht].InTraj != tjs.matchVec[indx].TjIDs[kpl]) continue;
                ++tjs.matchVec[indx].Count;
                break;
              } // indx
              if(indx == tjs.matchVec.size()) {
                // not found in the match vector so add it
                MatchStruct ms;
                ms.TjIDs.resize(3);
                ms.TjIDs[ipl] = tjs.fHits[iht].InTraj;
                ms.TjIDs[jpl] = tjs.fHits[jht].InTraj;
                ms.TjIDs[kpl] = tjs.fHits[kht].InTraj;
                ms.Count = 1;
                tjs.matchVec.push_back(ms);
              } // not found in the list
            } // kht
          } // jht
        } // iht
      } // jwire
    } // iwire
    
    if(tjs.matchVec.empty()) {
      if(prt) mf::LogVerbatim("TC")<<"Match3Views: no 3-plane matches found";
      return;
    }
    
    // sort by decreasing number of matches
    SortEntry se;
    std::vector<SortEntry> sortVec(tjs.matchVec.size());
    for(unsigned int indx = 0; indx < tjs.matchVec.size(); ++indx) {
      se.index = indx;
      se.length = tjs.matchVec[indx].Count;
      sortVec[indx] = se;
    }
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), greaterThan);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Match3Views matchVec sorted\n";
      for(unsigned int ii = 0; ii < tjs.matchVec.size(); ++ii) {
        unsigned int indx = sortVec[ii].index;
        myprt<<ii<<" Count "<<tjs.matchVec[indx].Count<<" TjIDs: ";
        for(auto& tjID : tjs.matchVec[indx].TjIDs) myprt<<" "<<tjID;
        myprt<<"\n";
      } // ii
    } // prt
    
    for(unsigned int ii = 0; ii < tjs.matchVec.size(); ++ii) {
      unsigned int indx = sortVec[ii].index;
      // skip this match if any of the trajectories is already matched or merged and killed
      bool skipit = false;
      // count the number of shower Tjs
      unsigned short nstj = 0;
      for(unsigned short ipl = 0; ipl < tjs.matchVec[indx].TjIDs.size(); ++ipl) {
        unsigned short itj = tjs.matchVec[indx].TjIDs[ipl] - 1;
        if(tjs.allTraj[itj].AlgMod[kMatch3D]) skipit = true;
        if(tjs.allTraj[itj].AlgMod[kKilled]) skipit = true;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) ++nstj;
      }
      if(skipit) continue;
      // Require 0 or 3 shower Tjs
      if(nstj != 0 && nstj != tjs.matchVec[indx].TjIDs.size()) continue;
      auto& imv = tjs.matchVec[indx];
      if(prt && imv.Count > 20) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Match3D indx "<<indx<<" Count "<<imv.Count;
      } // prt
      if(imv.TjIDs.empty()) continue;
      // Set the AlgMod bit
      for(unsigned short ipl = 0; ipl < imv.TjIDs.size(); ++ipl) {
        unsigned short itj = imv.TjIDs[ipl] - 1;
        Trajectory& tj = tjs.allTraj[itj];
        tj.AlgMod[kMatch3D] = true;
      } // ipl
      // Look for fragments of broken trajectories.
      // Allow one fragment for each trajectory (assumes that there is one Tj per plane...)
      std::vector<int> fragID(imv.TjIDs.size(), INT_MAX);
      // Stop searching when the triplet count gets low
      int minCount = 0.3 * imv.Count;
      for(unsigned int jj = ii + 1; jj < tjs.matchVec.size(); ++jj) {
        unsigned int jndx = sortVec[jj].index;
        if(tjs.matchVec[jndx].Count < minCount) break;
        auto& jmv = tjs.matchVec[jndx];
        // Count the number of ID matches in each matchVec entry
        unsigned short nmat = 0;
        // Keep track of the ID of an entry that is not matched - the presumed fragment
        int fID = INT_MAX;
        // and the index of fID in the matchVec entry
        unsigned short fIDIndex = 0;
        for(unsigned short jjj = 0; jjj < jmv.TjIDs.size(); ++jjj) {
          int jID = tjs.matchVec[jndx].TjIDs[jjj];
          if(jjj < imv.TjIDs.size() && jID == imv.TjIDs[jjj]) {
            ++nmat;
          } else if(!tjs.allTraj[jID - 1].AlgMod[kMatch3D]) {
            // jID wasn't found. Set fragID in case this is a fragment
            fID = jID;
            fIDIndex = jjj;
          }
        } // jjj
        // Look for two (out of three) matched ID
        if(nmat != 2) continue;
        if(fID == INT_MAX) continue;
        // The fragment should be fragID
        // Don't overwrite the fragID if one already exists
        if(fragID[fIDIndex] == INT_MAX) fragID[fIDIndex] = fID;
      } // jj
      // try to merge
      for(unsigned short jfg = 0; jfg < fragID.size(); ++jfg) {
        if(fragID[jfg] == INT_MAX) continue;
        if(prt) mf::LogVerbatim("TC")<<"matchVec index "<<indx<<" Try to merge "<<imv.TjIDs[jfg]<<" and "<<fragID[jfg];
        unsigned short ftj1 = imv.TjIDs[jfg] - 1;
        unsigned short ftj2 = fragID[jfg] - 1;
        // Put ftj1 first so that its AlgMod state is preserved - specifically the kMatch3D bit
        if(MergeAndStore(ftj1, ftj2, prt)) {
          // merge successfull. 
          int newTjID = tjs.allTraj[tjs.allTraj.size()-1].ID;
          if(prt) mf::LogVerbatim("TC")<<" merge successfull. newTjID "<<newTjID;
          std::replace(imv.TjIDs.begin(), imv.TjIDs.end(), imv.TjIDs[jfg], newTjID);
          // Set the count to 0 for all other entries that have these fragments
          // TODO Maybe it should be set to minCount instead...
          for(auto& jmv : tjs.matchVec) {
            if(std::find(jmv.TjIDs.begin(), jmv.TjIDs.end(), imv.TjIDs[jfg]) != jmv.TjIDs.end()) jmv.Count = 0;
            if(std::find(jmv.TjIDs.begin(), jmv.TjIDs.end(), fragID[jfg]) != jmv.TjIDs.end()) jmv.Count = 0;
          }
        } // successfull merge
        else if(prt) { mf::LogVerbatim("TC")<<" merge failed"; }
      } // jfg
      //  Index the matches that will become PFParticles into matchVecPFPList
      tjs.matchVecPFPList.push_back(indx);
    } // ii (indx)
  } // Match3Views

  //////////////////////////////////////////
  void TrajClusterAlg::Match2Views(const geo::TPCID& tpcid, const std::vector<float>& xx)
  {
    // Try to recover failed 3D matched trajectories using 2 views or match in 3D in a 2 view TPC. 
    // The xx vector of size tjs.fHits.size() has been pre-loaded by the calling routine
    
    // 2-view matching not requested
    if(tjs.Match3DCuts[2] < 0) return;
    
    unsigned int nTriple = tjs.matchVec.size();
    // make a temporary list of triple Tj match IDs to simplify searching
    std::vector<unsigned short> tripleTjList;
    for(auto& ms : tjs.matchVec) tripleTjList.insert(tripleTjList.end(), ms.TjIDs.begin(), ms.TjIDs.end());
    
    unsigned short minTjLen = tjs.Match3DCuts[2];
    
    // Three plane TPC - require a 2-plane match reside in a dead region of the 3rd plane?
    bool require3rdPlnDeadRegion = (tjs.NumPlanes == 3 && tjs.Match3DCuts[3] > 0);
    
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes - 1; ++ipl) {
      for(unsigned int iwire = tjs.FirstWire[ipl]; iwire < tjs.LastWire[ipl]; ++iwire) {
        if(tjs.WireHitRange[ipl][iwire].first < 0) continue;
        for(unsigned short jpl = ipl + 1; jpl < tjs.NumPlanes; ++jpl) {
          for(unsigned int jwire = tjs.FirstWire[jpl]; jwire < tjs.LastWire[jpl]; ++jwire) {
            if(tjs.WireHitRange[jpl][jwire].first < 0) continue;
            double yp, zp;
            tjs.geom->IntersectionPoint(iwire, jwire, ipl, jpl, tpcid.Cryostat, tpcid.TPC, yp, zp);
            // ensure this is inside the TPC
            if(yp < tjs.YLo || yp > tjs.YHi || zp < tjs.ZLo || zp > tjs.ZHi) continue;
            if(require3rdPlnDeadRegion) {
              unsigned short kpl = 3 - ipl - jpl;
              float fkwire = tjs.geom->WireCoordinate(yp, zp, kpl, tpcid.TPC, tpcid.Cryostat);
              if(fkwire < 0 || fkwire > tjs.MaxPos0[kpl]) continue;
              unsigned int kwire = std::nearbyint(fkwire);
              unsigned int lastWire = tjs.LastWire[2] - 1;
              if(kwire > lastWire) continue;
              if(tjs.WireHitRange[kpl][kwire].first == -1) continue;
            }
            // Have intersecting wire matches. Now look for time matches
            unsigned int ifirsthit = (unsigned int)tjs.WireHitRange[ipl][iwire].first;
            unsigned int ilasthit = (unsigned int)tjs.WireHitRange[ipl][iwire].second;
            unsigned int jfirsthit = (unsigned int)tjs.WireHitRange[jpl][jwire].first;
            unsigned int jlasthit = (unsigned int)tjs.WireHitRange[jpl][jwire].second;
            for(unsigned int iht = ifirsthit; iht < ilasthit; ++iht) {
              // Only consider hits that are associated with a trajectory
              if(tjs.fHits[iht].InTraj <= 0) continue;
              // see if this trajectory already has a 3D match
              unsigned short itj = tjs.fHits[iht].InTraj - 1;
              // or is in the triple match list
              if(std::find(tripleTjList.begin(), tripleTjList.end(), tjs.fHits[iht].InTraj) != tripleTjList.end()) continue;
              if(tjs.allTraj[itj].AlgMod[kMatch3D]) continue;
              // don't try to match Tjs in showers
              if(tjs.allTraj[itj].AlgMod[kInShower]) continue;
              if(tjs.allTraj[itj].Pts.size() < minTjLen) continue;
              for(unsigned int jht = jfirsthit; jht < jlasthit; ++jht) {
                // Only consider hits that are associated with a trajectory
                if(tjs.fHits[jht].InTraj <= 0) continue;
                // see if this trajectory already has a 3D match
                unsigned short jtj = tjs.fHits[jht].InTraj - 1;
                if(tjs.allTraj[jtj].AlgMod[kMatch3D]) continue;
                // or is in the triple match list
                if(std::find(tripleTjList.begin(), tripleTjList.end(), tjs.fHits[jht].InTraj) != tripleTjList.end()) continue;
                // don't try to match Tjs in showers
                if(tjs.allTraj[jtj].AlgMod[kInShower]) continue;
                if(tjs.allTraj[jtj].Pts.size() < minTjLen) continue;
                if(xx[jht] < xx[iht] - tjs.Match3DCuts[0]) continue;
                if(xx[jht] > xx[iht] + tjs.Match3DCuts[0]) break;
                // next see if the Tj IDs are in the match list
                unsigned short indx = nTriple;
                for(indx = nTriple; indx < tjs.matchVec.size(); ++indx) {
                  unsigned short nm = 0;
                  for(unsigned short ii = 0; ii < tjs.matchVec[indx].TjIDs.size(); ++ii) {
                    if(tjs.fHits[iht].InTraj == tjs.matchVec[indx].TjIDs[ii]) ++nm;
                    if(tjs.fHits[jht].InTraj == tjs.matchVec[indx].TjIDs[ii]) ++nm;
                  } // ii
                  if(nm == 2) {
                    ++tjs.matchVec[indx].Count;
                    break;
                  }
                } // indx
                if(indx == tjs.matchVec.size()) {
                  MatchStruct ms;
                  ms.TjIDs.resize(2);
                  ms.TjIDs[0] = tjs.fHits[iht].InTraj;
                  ms.TjIDs[1] = tjs.fHits[jht].InTraj;
                  ms.Count = 1;
                  tjs.matchVec.push_back(ms);
                } // not found in the list
              } // jht
            } // iht
          } // jwire
        } // jpl
      } // iwire
    } // ipl
    
    if(tjs.matchVec.size() == nTriple) {
      if(prt) mf::LogVerbatim("TC")<<"Match2Views: no 2-view matches found";
      return;
    }
    
    std::vector<SortEntry> sortVec;
    SortEntry se;
    for(unsigned int ii = nTriple; ii < tjs.matchVec.size(); ++ii) {
      se.index = ii;
      se.length = tjs.matchVec[ii].Count;
      sortVec.push_back(se);
    }
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), greaterThan);
    
    if(prt) {
      unsigned int nDouble = tjs.matchVec.size() - nTriple;
      mf::LogVerbatim myprt("TC");
      myprt<<"Match2Views: nTriple "<<nTriple<<" nDouble "<<nDouble<<"\n";
      for(unsigned int ii = 0; ii < sortVec.size(); ++ii) {
        myprt<<ii<<" indx "<<sortVec[ii].index;
        unsigned int indx = sortVec[ii].index;
        myprt<<" Count "<<tjs.matchVec[indx].Count;
        for(auto& tjID : tjs.matchVec[indx].TjIDs) myprt<<" "<<tjID;
        myprt<<"\n";
      } // ii
    } // prt
    
    // Index the matches in matchVecPFPList
    for(unsigned int ii = 0; ii < sortVec.size(); ++ii) {
      unsigned int indx = sortVec[ii].index;
      // skip this match if any of the trajectories is already matched
      bool skipit = false;
      unsigned short nstj = 0;
      for(unsigned short ipl = 0; ipl < tjs.matchVec[indx].TjIDs.size(); ++ipl) {
        unsigned short itj = tjs.matchVec[indx].TjIDs[ipl] - 1;
        if(tjs.allTraj[itj].AlgMod[kMatch3D]) skipit = true;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) ++nstj;
      }
      if(skipit) continue;
      if(nstj != 0 && nstj != tjs.matchVec[indx].TjIDs.size()) continue;
      for(unsigned short ipl = 0; ipl < tjs.matchVec[indx].TjIDs.size(); ++ipl) {
        int mtid = tjs.matchVec[indx].TjIDs[ipl];
        // Set the match flag
        tjs.allTraj[mtid - 1].AlgMod[kMatch3D] = true;
      } // ipl
      //  Put the subset of matches that will become PFParticles into a list
      tjs.matchVecPFPList.push_back(indx);
    } // ii (indx)
    
  } // Match2Views
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::Find3DEndPoints(const geo::TPCID& tpcid)
  {
    // Define the sXYZ and eXYZ endpoints of the matched trajectories. This is done by matching the
    // endpoints of the two longest trajectories in all planes. The endpoints and trajectories are then
    // reversed if necessary to put sXYZ[0] > eXYZ[0] and EndPt[0] => sXYZ 
    
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    
    for(auto& im : tjs.matchVecPFPList) {
      // a reference to a set of 3D matched trajectories
      auto& ms = tjs.matchVec[im];
      if(ms.TjIDs.empty()) continue;
      // ensure we are in the correct tpcid using the first Tj CTP
      unsigned short it1 = ms.TjIDs[0] - 1;
      geo::PlaneID plane1ID = DecodeCTP(tjs.allTraj[it1].CTP);
      if(plane1ID.Cryostat != cstat) continue;
      if(plane1ID.TPC != tpc) continue;
      // ignore shower Tjs. These are handled in FindShowerEndPoints
      if(ms.PDGCode == 1111) continue;
      // Check for the existence of a 3D vertex with these trajectories and if so define the
      // 3D vertex start and end indices. Note that in it's current state MatchHas3DVertex may not
      // always define eXYZ.
      if(Matched3DVtx(tjs, im) == 2) continue;
      // find the longest Tj
      unsigned short itjLong = USHRT_MAX;
      unsigned short ilen = 0;
      for(unsigned short ii = 0; ii < ms.TjIDs.size(); ++ii) {
        unsigned short itj = ms.TjIDs[ii] - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.Pts.size() > ilen) {
          ilen = tj.Pts.size();
          itjLong = itj;
        }
      } // ii
      if (itjLong == USHRT_MAX){
        mf::LogWarning("TC")<<"In Find3DEndPoints, itjLong = "<<itjLong<<", will return";
        return;
      }
      Trajectory& iTj = tjs.allTraj[itjLong];
      // find the 2nd longest tj
      unsigned short jtjLong = USHRT_MAX;
      unsigned short jlen = 0;
      for(unsigned short jj = 0; jj < ms.TjIDs.size(); ++jj) {
        unsigned short jtj = ms.TjIDs[jj] - 1;
        if(jtj == itjLong) continue;
        Trajectory& tj = tjs.allTraj[jtj];
        // ensure that this Tj is in a different plane than the longest one
        if(tj.CTP == iTj.CTP) continue;
        if(tj.Pts.size() > jlen) {
          jlen = tj.Pts.size();
          jtjLong = jtj;
        }
      } // ii
      if (jtjLong == USHRT_MAX){
        mf::LogWarning("TC")<<"In Find3DEndPoints, jtjLong = "<<jtjLong<<", will return";
        return;
      }
      // consider both ends of each Tj
      unsigned short iPln = DecodeCTP(iTj.CTP).Plane;
      Trajectory& jTj = tjs.allTraj[jtjLong];
      unsigned short jPln = DecodeCTP(jTj.CTP).Plane;
      unsigned short iEnd = 0, jEnd = 0;
      float bestDx = 1E6;
      float bestDxX = 0;
      for(unsigned short ie = 0; ie < 2; ++ie) {
        unsigned short iEndPt = iTj.EndPt[ie];
        TrajPoint& iTp = iTj.Pts[iEndPt];
        unsigned short iPln = DecodeCTP(iTp.CTP).Plane;
        float iX = tjs.detprop->ConvertTicksToX(iTp.Pos[1]/tjs.UnitsPerTick, iPln, tpc, cstat);
        for(unsigned short je = 0; je < 2; ++je) {
          unsigned short jEndPt = jTj.EndPt[je];
          TrajPoint& jTp = jTj.Pts[jEndPt];
          unsigned short jPln = DecodeCTP(jTp.CTP).Plane;
          float jX = tjs.detprop->ConvertTicksToX(jTp.Pos[1]/tjs.UnitsPerTick, jPln, tpc, cstat);
          float dx = std::abs(iX - jX);
          if(dx < bestDx) {
            bestDx = dx;
            iEnd = ie;
            jEnd = je;
            bestDxX = iX;
          }
        } // e2
      } // e1
      // Now decide which end should be the start, taking into account the possible existence of a 3D vertex.
      // Reference the Tps at the end with the best match
      unsigned short iEndPt = iTj.EndPt[iEnd];
      TrajPoint& iTp = iTj.Pts[iEndPt];
      unsigned short jEndPt = jTj.EndPt[jEnd];
      TrajPoint& jTp = jTj.Pts[jEndPt];
      std::array<float, 3> xyz;
      double yp, zp;
      unsigned int iwire = std::nearbyint(iTp.Pos[0]);
      if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, iPln, iwire))) continue;
      unsigned int jwire = std::nearbyint(jTp.Pos[0]);
      if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, jPln, jwire))) continue;
      tjs.geom->IntersectionPoint(iwire, jwire, iPln, jPln, cstat, tpc, yp, zp);
      // ensure this is inside the TPC
      if(yp < tjs.YLo) yp = tjs.YLo;
      if(yp > tjs.YHi) yp = tjs.YLo;
      if(zp < tjs.ZLo) zp = tjs.ZLo;
      if(zp > tjs.ZHi) zp = tjs.ZLo;
      xyz[0] = bestDxX;
      xyz[1] = yp;
      xyz[2] = zp;
      // find the xyz position of the other end
      iEndPt = iTj.EndPt[1 - iEnd];
      TrajPoint& ioTp = iTj.Pts[iEndPt];
      jEndPt = jTj.EndPt[1 - jEnd];
      TrajPoint& joTp = jTj.Pts[jEndPt];
      std::array<float, 3> oxyz;
      iwire = std::nearbyint(ioTp.Pos[0]);
      if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, iPln, iwire))) continue;
      jwire = std::nearbyint(joTp.Pos[0]);
      if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, jPln, jwire))) continue;
      tjs.geom->IntersectionPoint(iwire, jwire, iPln, jPln, cstat, tpc, yp, zp);
      // ensure this is inside the TPC
      if(yp < tjs.YLo) yp = tjs.YLo;
      if(yp > tjs.YHi) yp = tjs.YLo;
      if(zp < tjs.ZLo) zp = tjs.ZLo;
      if(zp > tjs.ZHi) zp = tjs.ZLo;
      geo::PlaneID iPlaneID = DecodeCTP(ioTp.CTP);
      float xi = tjs.detprop->ConvertTicksToX(ioTp.Pos[1]/tjs.UnitsPerTick, iPlaneID.Plane, tpc, cstat);
      geo::PlaneID jPlaneID = DecodeCTP(joTp.CTP);
      float xj = tjs.detprop->ConvertTicksToX(joTp.Pos[1]/tjs.UnitsPerTick, jPlaneID.Plane, tpc, cstat);
      oxyz[0] = 0.5 * (xi + xj);
      oxyz[1] = yp;
      oxyz[2] = zp;
      
      if(ms.sVtx3DIndex != USHRT_MAX) {
        // 3D vertex assignment exists at the start. Define the other end
        std::array<float, 3> vpos;
        vpos[0] = tjs.vtx3[ms.sVtx3DIndex].X;
        vpos[1] = tjs.vtx3[ms.sVtx3DIndex].Y;
        vpos[2] = tjs.vtx3[ms.sVtx3DIndex].Z;
        if(PosSep2(xyz, vpos) < PosSep2(oxyz, vpos)) {
          ms.eXYZ = oxyz;
        } else {
          ms.eXYZ = xyz;
        }
      } else {
        // No 3D vertex assignment exists at either end so define the end positions. A new 3D vertex will be created in FillPFPInfo
        if(xyz[0] > oxyz[0]) {
          ms.sXYZ = xyz;
          ms.eXYZ = oxyz;
        } else {
          ms.sXYZ = oxyz;
          ms.eXYZ = xyz;
        }
      } // No 3D vertex assignment exists at either end
      
      // Reverse trajectories as necessary so that EndPt[0] is at the sXYZ end
      std::array<float, 2> vpos;
      for(unsigned short ii = 0; ii < ms.TjIDs.size(); ++ii) {
        unsigned short itj = ms.TjIDs[ii] - 1;
        Trajectory& tj = tjs.allTraj[itj];
        unsigned short endPt0 = tj.EndPt[0];
        unsigned short endPt1 = tj.EndPt[1];
        // Project sXYZ to this plane coordinate system
        geo::PlaneID planeID = DecodeCTP(tj.CTP);
        vpos[0] = tjs.geom->WireCoordinate(ms.sXYZ[1], ms.sXYZ[2], planeID) + 0.5;
        if(vpos[0] < 0 || (unsigned int)vpos[0] > tjs.NumWires[planeID.Plane]) continue;
        vpos[1] = tjs.detprop->ConvertXToTicks(ms.sXYZ[0], planeID.Plane, tpc, cstat) * tjs.UnitsPerTick;
        // Reverse if end 0 is further away from vpos than end 1
        if(PosSep2(tj.Pts[endPt0].Pos, vpos) > PosSep2(tj.Pts[endPt1].Pos, vpos)) ReverseTraj(tjs, tj);
      } // ii
     if(prt)  mf::LogVerbatim("TC")<<"F3DEP: im "<<im<<" itjLong "<<itjLong<<" iEnd "<<iEnd<<" jtjLong "<<jtjLong<<" jEnd "<<jEnd;

    } // im (ms)
        
  } // Find3DEndPoints
  
  //////////////////////////////////////////
  void TrajClusterAlg::FillPFPInfo()
  {
    // Fills the PFParticle info in matchVec. Each entry contains a list of trajectories and has a defined
    // start and end XYZ position, possibly with the index of a 3D vertex. 
    
    unsigned short pfpCount = 0;
    for(auto& im : tjs.matchVecPFPList) {
      auto& ms = tjs.matchVec[im];
      // ms contains a list of tjs that were matched in 3D. 
      if(ms.Count == 0) continue;
      // Set the parent index to itself for now
      ms.ParentMSIndex = pfpCount;
      ++pfpCount;
      // Trajectories have been reversed so that the Tj and hit order are consistent. The
      // next step is to decide whether the start of the PFParticle (where the 3D vertex association will be made) is indeed the
      // start of the trajectory from a physics standpoint, e.g. dQ/dx, muon delta-ray tag, cosmic rays entering the detector, etc.
      // Make this decision in a utility function.
      if(Reverse3DMatchTjs(tjs, im, prt)) {
        // All of the trajectories should be reversed
        for(auto& tjID : ms.TjIDs) {
          unsigned short itj = tjID - 1;
          ReverseTraj(tjs, tjs.allTraj[itj]);
        } // tjID
        if(prt)  mf::LogVerbatim("TC")<<" Reversed all Tjs using advice from Reverse3DMatchTjs";
        // swap the matchVec end info also
        std::swap(ms.sXYZ, ms.eXYZ);
        std::swap(ms.sVtx3DIndex, ms.eVtx3DIndex);
      } // Reverse3DMatchTjs
      // we may need to clobber 2D vertices so start a list and count them
      std::vector<unsigned short> vtxIDs, vtxIDCnt;
      for(auto& tjID : ms.TjIDs) {
        unsigned short itj = tjID - 1;
        Trajectory& tj = tjs.allTraj[itj];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] > 0) {
            // Count the number of times this 2D vertex is used.
            unsigned short indx;
            for(indx = 0; indx < vtxIDs.size(); ++indx) if(tj.VtxID[end] == vtxIDs[indx]) break;
            if(indx == vtxIDs.size()) {
              // not found so add it
              vtxIDs.push_back(tj.VtxID[end]);
              vtxIDCnt.push_back(1);
            } else {
              ++vtxIDCnt[indx];
            }
          } // 2D vertex exists
        } // end
      } // tjID
      // clobber any 2D vertices that exist between trajectories in the match list
      for(unsigned short ii = 0; ii < vtxIDs.size(); ++ii) {
        if(vtxIDCnt[ii] > 1) MakeVertexObsolete(tjs, vtxIDs[ii]);
      } // ii
      // Make a 3D vertex at the start of the PFParticle if one doesn't already exist
      if(ms.sVtx3DIndex == USHRT_MAX) {
        // get a reference to one of the Tjs to define the cryostat and TPC
        unsigned short itj = ms.TjIDs[0] - 1;
        geo::PlaneID planeID = DecodeCTP(tjs.allTraj[itj].CTP);
        Vtx3Store newVx3;
        newVx3.CStat = planeID.Cryostat;
        newVx3.TPC = planeID.TPC;
        // Set Wire < 0 as a flag that this is a "complete" 3D vertex even though no 2D vertices have been made.
        newVx3.Wire = -2;
        newVx3.X = ms.sXYZ[0];
        newVx3.Y = ms.sXYZ[1];
        newVx3.Z = ms.sXYZ[2];
        tjs.vtx3.push_back(newVx3);
        ms.sVtx3DIndex = tjs.vtx3.size() - 1;
        if(prt) mf::LogVerbatim("TC")<<" Made 3D start vertex "<<tjs.vtx3.size() - 1<<" at "<<newVx3.X<<" "<<newVx3.Y<<" "<<newVx3.Z;
      }
    } // im (ms)

  } // FillPFPInfo
  
  //////////////////////////////////////////
  void TrajClusterAlg::StepCrawl(Trajectory& tj)
  {
    // Crawl along the direction specified in the traj vector in steps of size step
    // (wire spacing equivalents). Find hits between the last trajectory point and
    // the last trajectory point + step. A new trajectory point is added if hits are
    // found. Crawling continues until no signal is found for two consecutive steps
    // or until a wire or time boundary is reached.
    
    fGoodTraj = false;
    fTryWithNextPass = false;
    if(tj.Pts.empty()) return;
    
    if(fCTP != tj.CTP || !WireHitRangeOK(tjs, tj.CTP)) {
      std::cout<<"StepCrawl: Warning fCTP != tj.CTP or invalid WireHitRange.\n";
      fQuitAlg = true;
      return;
    }
 
    unsigned short lastPt;
    unsigned short lastPtWithUsedHits = tj.EndPt[1];
    unsigned short lastPtWithHits;

    lastPt = lastPtWithUsedHits;
    // Construct a local TP from the last TP that will be moved on each step.
    // Only the Pos and Dir variables will be used
    TrajPoint ltp;
    ltp.CTP = tj.CTP;
    ltp.Pos = tj.Pts[lastPt].Pos;
    ltp.Dir = tj.Pts[lastPt].Dir;
    // A second TP is cloned from the leading TP of tj, updated with hits, fit
    // parameters,etc and possibly pushed onto tj as the next TP
    TrajPoint tp ;
    
    // assume it is good from here on
    fGoodTraj = true;
    
    unsigned short nMissedSteps = 0;

    bool sigOK, keepGoing;
    unsigned short killPts;
    for(unsigned short step = 1; step < 10000; ++step) {
      // make a copy of the previous TP
      lastPt = tj.Pts.size() - 1;
      tp = tj.Pts[lastPt];
      ++tp.Step;
      double stepSize = fVLAStepSize;
      if(tp.AngleCode < 2) stepSize = std::abs(1/ltp.Dir[0]);
      // move the local TP position by one step in the right direction
      for(unsigned short iwt = 0; iwt < 2; ++iwt) ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;

      unsigned short ivx = TPNearVertex(tjs, ltp);
      if(ivx != USHRT_MAX) {
        // Trajectory stops near a vertex so make the assignment
        tj.StopFlag[1][kAtVtx] = true;
        tj.VtxID[1] = tjs.vtx[ivx].ID;
        break;
      }

      SetPDGCode(tjs, tj);

      // copy this position into tp
      tp.Pos = ltp.Pos;
      tp.Dir = ltp.Dir;
      if(prt) {
        mf::LogVerbatim("TC")<<"StepCrawl "<<step<<" Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" Dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" stepSize "<<stepSize<<" AngCode "<<tp.AngleCode;
      }
      // hit the boundary of the TPC?
      if(tp.Pos[0] < 0 || tp.Pos[0] > tjs.MaxPos0[fPlane] ||
         tp.Pos[1] < 0 || tp.Pos[1] > tjs.MaxPos1[fPlane]) break;
      // remove the old hits and other stuff
      tp.Hits.clear();
      tp.UseHit.reset();
      tp.FitChi = 0; tp.Chg = 0;
      // append to the trajectory
      tj.Pts.push_back(tp);
      // update the index of the last TP
      lastPt = tj.Pts.size() - 1;
      // look for hits
      AddHits(tj, lastPt, sigOK);
      // Check the stop flag
      if(tj.StopFlag[1][kAtTj]) break;
      // If successfull, AddHits has defined UseHit for this TP,
      // set the trajectory endpoints, and define HitPos.
      if(tj.Pts[lastPt].Hits.empty()) {
        // Require three points with charge on adjacent wires for small angle
        // stepping.
        if(tj.Pts[lastPt].AngleCode == 0 && lastPt == 2) return;
        // No close hits added.
        ++nMissedSteps;
        // First check for no signal in the vicinity
        if(lastPt > 0) {
          // break if this is a reverse propagate activity and there was no signal (not on a dead wire)
          if(!sigOK && tj.AlgMod[kRevProp]) break;
          // Ensure that there is a signal here after missing a number of steps on a LA trajectory
          if(tj.Pts[lastPt].AngleCode > 0 && nMissedSteps > 4 && !SignalAtTp(tjs, ltp)) {
            tj.StopFlag[1][kSignal] = false;
            break;
          }
          // the last point with hits (used or not) is the previous point
          lastPtWithHits = lastPt - 1;
          float tps = TrajPointSeparation(tj.Pts[lastPtWithHits], ltp);
          float dwc = DeadWireCount(tjs, ltp, tj.Pts[lastPtWithHits]);
          float nMissedWires = tps * std::abs(ltp.Dir[0]) - dwc;
          float maxWireSkip = fMaxWireSkipNoSignal;
          if(tj.PDGCode == 13) maxWireSkip = tjs.MuonTag[2];
          if(prt) mf::LogVerbatim("TC")<<" StepCrawl: no signal at ltp "<<PrintPos(tjs, ltp)<<" nMissedWires "<<std::fixed<<std::setprecision(1)<<nMissedWires<<" dead wire count "<<dwc<<" maxWireSkip "<<maxWireSkip<<" tj.PGDCode "<<tj.PDGCode;
          if(nMissedWires > maxWireSkip) {
            // We passed a number of wires without adding hits and are ready to quit.
            // First see if there is one good unused hit on the end TP and if so use it
            // lastPtWithHits + 1 == lastPt && tj.Pts[lastPtWithHits].Chg == 0 && tj.Pts[lastPtWithHits].Hits.size() == 1
            if(tj.EndPt[1] < tj.Pts.size() - 1 && tj.Pts[tj.EndPt[1]+1].Hits.size() == 1) {
              unsigned short lastLonelyPoint = tj.EndPt[1] + 1;
              unsigned int iht = tj.Pts[lastLonelyPoint].Hits[0];
              if(tjs.fHits[iht].InTraj == 0 && tj.Pts[lastLonelyPoint].Delta < 3 * tj.Pts[lastLonelyPoint].DeltaRMS) {
                tjs.fHits[iht].InTraj = tj.ID;
                tj.Pts[lastLonelyPoint].UseHit[0] = true;
                DefineHitPos(tj.Pts[lastLonelyPoint]);
                SetEndPoints(tjs, tj);
                if(prt) {
                  mf::LogVerbatim("TC")<<" Added a Last Lonely Hit before breaking ";
                  PrintTrajPoint("LLH", tjs, lastPt, tj.StepDir, tj.Pass, tj.Pts[lastLonelyPoint]);
                }
              }
            }
            break;
          }
        } // lastPt > 0
        // no sense keeping this TP on tj if no hits were added
        tj.Pts.pop_back();
        continue;
      } // tj.Pts[lastPt].Hits.empty()
      // ensure that we actually moved
      if(lastPt > 0 && PosSep2(tj.Pts[lastPt].Pos, tj.Pts[lastPt-1].Pos) < 0.1) return;
      // Found hits at this location so reset the missed steps counter
      nMissedSteps = 0;
      // Update the last point fit, etc using the just added hit(s)
      UpdateTraj(tj);
      // a failure occurred
      if(!fUpdateTrajOK) return;
      if(tj.Pts[lastPt].Chg == 0) {
        // There are points on the trajectory by none used in the last step. See
        // how long this has been going on
        float tps = TrajPointSeparation(tj.Pts[tj.EndPt[1]], ltp);
        float dwc = DeadWireCount(tjs, ltp, tj.Pts[tj.EndPt[1]]);
        float nMissedWires = tps * std::abs(ltp.Dir[0]) - dwc;
        if(prt)  mf::LogVerbatim("TC")<<" Hits exist on the trajectory but are not used. Missed wires "<<nMissedWires<<" dead wire count "<<(int)dwc;
        // break if this is a reverse propagate activity with no dead wires
        if(tj.AlgMod[kRevProp] && dwc == 0) break;
        if(nMissedWires > fMaxWireSkipWithSignal) break;
        // try this out
        if(!MaskedHitsOK(tj)) {
          return;
        }
        // check for a series of bad fits and stop stepping
        if(tjs.UseAlg[kStopBadFits] && nMissedWires > 4 && StopIfBadFits(tj)) break;
        // Keep stepping
        if(prt) {
          if(tj.AlgMod[kRevProp]) {
            PrintTrajectory("RP", tjs, tj, lastPt);
          } else {
            PrintTrajectory("SC", tjs, tj, lastPt);
          }
        }
        continue;
      } // tp.Hits.empty()
      if(tj.Pts.size() == 3) {
        // ensure that the last hit added is in the same direction as the first two.
        // This is a simple way of doing it
        if(PosSep2(tj.Pts[0].HitPos, tj.Pts[2].HitPos) < PosSep2(tj.Pts[0].HitPos, tj.Pts[1].HitPos)) return;
        // ensure that this didn't start as a small angle trajectory and immediately turn
        // into a large angle one
        if(tj.Pts[lastPt].AngleCode > fMaxAngleCode[tj.Pass]) {
          if(prt) mf::LogVerbatim("TC")<<" Wandered into an invalid angle range. Quit stepping.";
          fGoodTraj = false;
          return;
        }
      } // tj.Pts.size() == 3
      // Update the local TP with the updated position and direction
      ltp.Pos = tj.Pts[lastPt].Pos;
      ltp.Dir = tj.Pts[lastPt].Dir;
      if(fMaskedLastTP) {
        // see if TPs have been masked off many times and if the
        // environment is clean. If so, return and try with next pass
        // cuts
        if(!MaskedHitsOK(tj)) {
          if(prt) {
            if(tj.AlgMod[kRevProp]) {
              PrintTrajectory("RP", tjs, tj, lastPt);
            } else {
              PrintTrajectory("SC", tjs, tj, lastPt);
            }
          }
          return;
        }
        // Don't bother with the rest of the checking below if we
        // set all hits not used on this TP
        if(prt) {
          if(tj.AlgMod[kRevProp]) {
            PrintTrajectory("RP", tjs, tj, lastPt);
          } else {
            PrintTrajectory("SC", tjs, tj, lastPt);
          }
        }
        continue;
      }
      // We have added a TP with hits
      // assume that we aren't going to kill the point we just added, or any
      // of the previous points...
      killPts = 0;
      // assume that we should keep going after killing points
      keepGoing = true;
      // check for a kink. Stop crawling if one is found
      GottaKink(tj, killPts);
      if(tj.StopFlag[1][kAtKink]) keepGoing = false;
      // See if the Chisq/DOF exceeds the maximum.
      // UpdateTraj should have reduced the number of points fit
      // as much as possible for this pass, so this trajectory is in trouble.
      if(killPts == 0 &&  tj.Pts[lastPt].FitChi > fMaxChi && tj.PDGCode != 13) {
        if(prt) mf::LogVerbatim("TC")<<"   bad FitChi "<<tj.Pts[lastPt].FitChi<<" cut "<<fMaxChi;
        fGoodTraj = (NumPtsWithCharge(tjs, tj, true) > fMinPtsFit[tj.Pass]);
        return;
      }
      // print the local tp unless we have killing to do
      if(killPts == 0) {
        if(prt) {
          if(tj.AlgMod[kRevProp]) {
            PrintTrajectory("RP", tjs, tj, lastPt);
          } else {
            PrintTrajectory("SC", tjs, tj, lastPt);
          }
        }
      } else {
        MaskTrajEndPoints(tj, killPts);
        if(!fGoodTraj) return;
        unsigned int onWire = (float)(std::nearbyint(tj.Pts[lastPt].Pos[0]));
        float nSteps = (float)(step - tj.Pts[lastPt - killPts].Step);
        if(prt) mf::LogVerbatim("TC")<<"TRP   killing "<<killPts<<" after "<<nSteps<<" steps from prev TP.  Current tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1];
        // move the position
        tj.Pts[lastPt].Pos[0] += nSteps * tj.Pts[lastPt].Dir[0];
        tj.Pts[lastPt].Pos[1] += nSteps * tj.Pts[lastPt].Dir[1];
        if(tj.Pts[lastPt].AngleCode == 0) {
          // put the TP at the wire position prior to the move
          float dw = onWire - tj.Pts[lastPt].Pos[0];
          tj.Pts[lastPt].Pos[0] = onWire;
          tj.Pts[lastPt].Pos[1] += dw * tj.Pts[lastPt].Dir[1] / tj.Pts[lastPt].Dir[0];
        }
        // check the MCSMom after we going
        if(tj.Pts.size() > 20 && tj.Pass < fMinMCSMom.size() && tj.MCSMom < fMinMCSMom[tj.Pass]) break;
        // copy to the local trajectory point
        ltp.Pos = tj.Pts[lastPt].Pos;
        ltp.Dir = tj.Pts[lastPt].Dir;
        if(prt) mf::LogVerbatim("TC")<<"  New ltp.Pos     "<<ltp.Pos[0]<<" "<<ltp.Pos[1]<<" ticks "<<(int)ltp.Pos[1]/tjs.UnitsPerTick;
        if(!keepGoing) break;
      }
    } // step
    
    if(prt) mf::LogVerbatim("TC")<<"End StepCrawl with tj size "<<tj.Pts.size()<<" fGoodTraj = "<<fGoodTraj<<" with fTryWithNextPass "<<fTryWithNextPass;

  } // StepCrawl
  
  ////////////////////////////////////////////////
  bool TrajClusterAlg::IsGhost(Trajectory& tj)
  {
    // Sees if trajectory tj shares many hits with another trajectory and if so merges them.
    
    if(!tjs.UseAlg[kUseGhostHits]) return false;
    // ensure that tj is not a saved trajectory
    if(tj.ID > 0) return true;
    // or an already killed trajectory
    if(tj.AlgMod[kKilled]) return true;
    if(tj.Pts.size() < 3) return false;
    
    // vectors of traj IDs, and the occurrence count
    std::vector<unsigned short> tID, tCnt;
    
    unsigned short hitCnt = 0;
    unsigned short nAvailable = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        // ignore hits used by this trajectory
        if(tj.Pts[ipt].UseHit[ii]) {
          ++hitCnt;
          continue;
        }
        unsigned int iht = tj.Pts[ipt].Hits[ii];
        if(tjs.fHits[iht].InTraj > 0) {
          unsigned short itj = tjs.fHits[iht].InTraj;
          unsigned short indx;
          for(indx = 0; indx < tID.size(); ++indx) if(tID[indx] == itj) break;
          if(indx == tID.size()) {
            tID.push_back(itj);
            tCnt.push_back(1);
          } else {
            ++tCnt[indx];
          }
        } else {
          ++nAvailable;
        }
      } // ii
    } // ipt
    
    // Call it a ghost if > 1/3 of the hits are used by another trajectory
    hitCnt /= 3;
    int oldTjID = INT_MAX;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"IsGhost tj hits size cut "<<hitCnt<<" tID_tCnt";
      for(unsigned short ii = 0; ii < tCnt.size(); ++ii) myprt<<" "<<tID[ii]<<"_"<<tCnt[ii];
      myprt<<"\nAvailable hits "<<nAvailable;
    } // prt
    
    for(unsigned short ii = 0; ii < tCnt.size(); ++ii) {
      if(tCnt[ii] > hitCnt) {
        oldTjID = tID[ii];
        hitCnt = tCnt[ii];
      }
    } // ii
    if(oldTjID == INT_MAX) return false;
    int oldTjIndex = oldTjID - 1;
    
    // See if this looks like a short delta-ray on a long muon
    Trajectory& oTj = tjs.allTraj[oldTjIndex];
    if(oTj.PDGCode == 13 && hitCnt < 0.1 * oTj.Pts.size()) return false;
    
    // See if there are gaps in this trajectory indicating that it is really a ghost and not
    // just a crossing trajectory 
    // find the range of wires spanned by oTj
    int wire0 = INT_MAX;
    int wire1 = 0;
    for(auto& otp : oTj.Pts) {
      int wire = std::nearbyint(otp.Pos[0]);
      if(wire < wire0) wire0 = wire;
      if(wire > wire1) wire1 = wire;
    } // tp
    
    unsigned short nwires = wire1 - wire0 + 1;
    std::vector<float> oTjPos1(nwires, -1);
    unsigned short nMissedWires = 0;
    for(unsigned short ipt = oTj.EndPt[0]; ipt <= oTj.EndPt[1]; ++ipt) {
      if(oTj.Pts[ipt].Chg == 0) continue;
      int wire = std::nearbyint(oTj.Pts[ipt].Pos[0]);
      int indx = wire - wire0;
      if(indx < 0 || indx > nwires - 1) {
        std::cout<<"Ooops "<<ipt<<" "<<indx<<" tj.ID "<<tj.ID<<"\n";
        PrintTrajectory("tj", tjs, tj, USHRT_MAX);
        PrintTrajectory("oTj", tjs, oTj, USHRT_MAX);
        continue;
      }
      oTjPos1[indx] = oTj.Pts[ipt].Pos[1];
      ++nMissedWires;
    } // ipt
    // count the number of ghost TPs
    unsigned short ngh = 0;
    // and the number with Delta > 0 relative to oTj
    unsigned short nghPlus = 0;
    // keep track of the first point and last point appearance of oTj
    unsigned short firstPtInoTj = USHRT_MAX;
    unsigned short lastPtInoTj = 0;
    TrajPoint tp = tj.Pts[tj.EndPt[0]];
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) {
        tp = tj.Pts[ipt];
        continue;
      }
      int wire = std::nearbyint(tj.Pts[ipt].Pos[0]);
      int indx = wire - wire0;
      if(indx < 0 || indx > nwires - 1) continue;
      if(oTjPos1[indx] > 0) {
        // ensure that the hits in this tp are used in oTj
        bool HitInoTj = false;
        for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(tjs.fHits[iht].InTraj ==  oldTjID) HitInoTj = true;
        } // ii
        if(HitInoTj) {
          ++ngh;
          MoveTPToWire(tp, tj.Pts[ipt].Pos[0]);
          if(tp.Pos[1] > oTjPos1[indx]) ++nghPlus;
          if(firstPtInoTj == USHRT_MAX) firstPtInoTj = ipt;
          lastPtInoTj = ipt;
        }
      } // oTjHasChg[indx]
    } // ipt
    
    if(prt) mf::LogVerbatim("TC")<<" Number of missed wires in oTj gaps "<<nMissedWires<<" Number of ghost hits in these gaps "<<ngh<<" nghPlus "<<nghPlus<<" cut "<<0.2 * nMissedWires;
    
    if(ngh < 0.2 * nMissedWires) return false;
    
    // require all of the tj TPs to be on either the + or - side of the oTj trajectory
    if(!(nghPlus > 0.8 * ngh || nghPlus < 0.2 * ngh) ) return false;
    
    if(prt) mf::LogVerbatim("TC")<<" Trajectory is a ghost of "<<oldTjID<<" first point in oTj "<<firstPtInoTj<<" last point "<<lastPtInoTj;
    
    // unset all of the shared hits
    for(unsigned short ipt = firstPtInoTj; ipt <= lastPtInoTj; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      UnsetUsedHits(tjs, tj.Pts[ipt]);
      if(prt) PrintTrajectory("IG", tjs, tj, ipt);
    }
    // see how many points are left at the end
    ngh = 0;
    for(unsigned short ipt = lastPtInoTj; ipt <= tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ngh;
    } // ipt
    // clobber those too?
    if(ngh > 0 && ngh < fMinPts[tj.Pass]) {
      for(unsigned short ipt = lastPtInoTj; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg > 0) UnsetUsedHits(tjs, tj.Pts[ipt]);
      } // ipt
    }
    SetEndPoints(tjs, tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    tjs.allTraj[oldTjIndex].AlgMod[kUseGhostHits] = true;
    TrimEndPts(tjs, tj, fQualityCuts, prt);
    if(tj.AlgMod[kKilled]) {
      fGoodTraj = false;
      return true;
    }
    tj.MCSMom = MCSMom(tjs, tj);
    if(prt)  mf::LogVerbatim("TC")<<" New tj size "<<tj.Pts.size();
    return true;
    
  } // IsGhost
  
  ////////////////////////////////////////////////
  bool TrajClusterAlg::IsGhost(std::vector<unsigned int>& tHits, unsigned short& ofTraj)
  {
    // Called by FindJunkTraj to see if the passed hits are close to an existing
    // trajectory and if so, they will be used in that other trajectory
    
    ofTraj = USHRT_MAX;
    
    if(!tjs.UseAlg[kUseGhostHits]) return false;
    
    if(tHits.size() < 2) return false;
    // find all nearby hits
    std::vector<unsigned int> hitsInMuliplet, nearbyHits;
    for(auto iht : tHits) {
      GetHitMultiplet(iht, hitsInMuliplet);
      // prevent double counting
      for(auto mht : hitsInMuliplet) {
        if(std::find(nearbyHits.begin(), nearbyHits.end(), mht) == nearbyHits.end()) {
          nearbyHits.push_back(mht);
        }
      } // mht
    } // iht
    
    // vectors of traj IDs, and the occurrence count
    std::vector<unsigned short> tID, tCnt;
    unsigned short itj, indx;
    for(auto iht : nearbyHits) {
      if(tjs.fHits[iht].InTraj <= 0) continue;
      itj = tjs.fHits[iht].InTraj;
      for(indx = 0; indx < tID.size(); ++indx) if(tID[indx] == itj) break;
      if(indx == tID.size()) {
        tID.push_back(itj);
        tCnt.push_back(1);
      }  else {
        ++tCnt[indx];
      }
    } // iht
    if(tCnt.empty()) return false;
    
    // Call it a ghost if > 50% of the hits are used by another trajectory
    unsigned short tCut = 0.5 * tHits.size();
    unsigned short ii, jj;
    itj = USHRT_MAX;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"IsGhost tHits size "<<tHits.size()<<" cut fraction "<<tCut<<" tID_tCnt";
      for(ii = 0; ii < tCnt.size(); ++ii) myprt<<" "<<tID[ii]<<"_"<<tCnt[ii];
    } // prt
    
    for(ii = 0; ii < tCnt.size(); ++ii) {
      if(tCnt[ii] > tCut) {
        itj = tID[ii] - 1;
        break;
      }
    } // ii
    if(itj > tjs.allTraj.size() - 1) return false;
    
    if(prt) mf::LogVerbatim("TC")<<"is ghost of trajectory "<<tjs.allTraj[itj].ID;

    // Use all hits in tHits that are found in itj
    unsigned int iht, tht;
    for(auto& tp : tjs.allTraj[itj].Pts) {
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        iht = tp.Hits[ii];
        if(tjs.fHits[iht].InTraj != 0) continue;
        for(jj = 0; jj < tHits.size(); ++jj) {
          tht = tHits[jj];
          if(tht != iht) continue;
          tp.UseHit[ii] = true;
          tjs.fHits[iht].InTraj = tjs.allTraj[itj].ID;
          break;
        } // jj
      } // ii
    } // tp
    tjs.allTraj[itj].AlgMod[kUseGhostHits] = true;
    ofTraj = itj;
    return true;
    
  } // IsGhost

  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckTraj(Trajectory& tj)
  {
    // Check the quality of the trajectory and possibly trim it or flag it for deletion
    
    if(!fGoodTraj) return;
    
    fTryWithNextPass = false;

    // ensure that the end points are defined
    SetEndPoints(tjs, tj);
    if(tj.EndPt[0] == tj.EndPt[1]) return;
    
    tj.MCSMom = MCSMom(tjs, tj);
    
    if(prt) {
      mf::LogVerbatim("TC")<<"inside CheckTraj with tj.Pts.size = "<<tj.Pts.size()<<" MCSMom "<<tj.MCSMom;
    }
    
    // remove any points at the end that don't have charge
    tj.Pts.resize(tj.EndPt[1] + 1);

    // Ensure that a hit only appears once in the TJ
    if(HasDuplicateHits(tjs, tj, prt)) {
      if(prt) mf::LogVerbatim("TC")<<" HasDuplicateHits ";
       fGoodTraj = false;
      return;
    }
    
    // See if this is a ghost trajectory
    if(IsGhost(tj)) {
      if(prt) mf::LogVerbatim("TC")<<" CT: Ghost trajectory - trimmed hits ";
    }
    
     // checks are different for Very Large Angle trajectories
    bool isVLA = (tj.Pts[tj.EndPt[1]].AngleCode == 2);
    // The last two ranges are Large Angle and Very Large Angle. Determine if the TJ is Small Angle
    bool isSA = (tj.Pts[tj.EndPt[1]].AngleCode == 0);
    
    // First remove any TPs at the end that have no hits after
    // setting the StopFlag. Assume that there are no hits on TPs after the end
    tj.StopFlag[1][kSignal] = false;
    if(tj.EndPt[1] < tj.Pts.size() - 1) {
      // There must be hits at the end so set the kSignal StopFlag
      if(!tj.Pts[tj.EndPt[1]+1].Hits.empty()) tj.StopFlag[1][kSignal] = true;
    }
    tj.Pts.resize(tj.EndPt[1] + 1);

    // Fill in any gaps with hits that were skipped, most likely delta rays on muon tracks
    if(!isVLA) FillGaps(tj);
    
    if(prt) mf::LogVerbatim("TC")<<" CheckTraj MCSMom "<<tj.MCSMom<<" isVLA? "<<isVLA<<" NumPtsWithCharge "<<NumPtsWithCharge(tjs, tj, false)<<" Min Req'd "<<fMinPts[tj.Pass];
    
    if(NumPtsWithCharge(tjs, tj, false) < fMinPts[tj.Pass]) {
      fGoodTraj = false;
      return;
    }
    
    // Check for hit width consistency on short trajectories
    if(tj.Pts.size() < 10) {
      float maxWidth = 0;
      float minWidth = 999;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        if(tj.Pts[ipt].HitPosErr2 > maxWidth) maxWidth = tj.Pts[ipt].HitPosErr2;
        if(tj.Pts[ipt].HitPosErr2 < minWidth) minWidth = tj.Pts[ipt].HitPosErr2;
      } // ipt
      // Require less than a 3X difference in the hit width or 10X for HitPosErr2
      if(maxWidth > 10 * minWidth) {
        if(prt) mf::LogVerbatim("TC")<<" TP width variation too large: minWidth "<<minWidth<<" maxWidth "<<maxWidth;
        fGoodTraj = false;
        return;
      }
    } // short trajectory

    // Trim the end points until the TJ meets the quality cuts
    TrimEndPts(tjs, tj, fQualityCuts, prt);
    if(tj.AlgMod[kKilled]) {
      fGoodTraj = false;
      return;
    }
    
    // Check for a Bragg peak at both ends. This may be used by FixTrajBegin.
    ChkStop(tj);

    // Update the trajectory parameters at the beginning of the trajectory
    FixTrajBegin(tj);

    // ignore short trajectories
    if(tj.EndPt[1] < 4) return;
    
    if(isSA && !tj.StopFlag[1][kBragg]) {
      // Small angle checks

      if(tjs.UseAlg[kCTKink] && tj.EndPt[1] > 8 && !tj.StopFlag[1][kAtKink] && tj.MCSMom > 50) {
        // look for the signature of a kink near the end of the trajectory.
        // These are: Increasing delta for the last few hits
        unsigned short newSize = USHRT_MAX;
        unsigned short lastPtToChk = tj.EndPt[1] - 4;
        float deltaCut = 2 * tj.Pts[lastPtToChk].DeltaRMS;
        for(unsigned short ipt = tj.EndPt[1]; ipt > lastPtToChk; --ipt) {
          // Stop checking if delta is good
          if(tj.Pts[ipt].Delta < deltaCut) break;
          float drat = tj.Pts[ipt].Delta / tj.Pts[ipt-1].Delta;
          if(drat > 1.2) newSize = ipt;
        } // ipt
        if(newSize != USHRT_MAX) {
          if(prt) mf::LogVerbatim("TC")<<"CTKink: Masking end points to newSize "<<newSize;
          for(unsigned short ipt = newSize; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
          SetEndPoints(tjs, tj);
          tj.AlgMod[kCTKink] = true;
        }
      } // tjs.UseAlg[kCTKink]

      if(tjs.UseAlg[kCTStepChk] && !tj.AlgMod[kRevProp]) {
        // Compare the number of steps taken per TP near the beginning and
        // at the end. This will get confused if RevProp is used
        short nStepBegin = tj.Pts[2].Step - tj.Pts[1].Step;
        short nStepEnd;
        unsigned short lastPt = tj.Pts.size() - 1;
        unsigned short newSize = tj.Pts.size();
        for(unsigned short ipt = lastPt; ipt > lastPt - 2; --ipt) {
          nStepEnd = tj.Pts[ipt].Step - tj.Pts[ipt - 1].Step;
          if(nStepEnd > 3 * nStepBegin) newSize = ipt;
        }
        if(prt) mf::LogVerbatim("TC")<<"CTStepChk: check number of steps. newSize "<<newSize<<" tj.Pts.size() "<<tj.Pts.size();
        if(newSize < tj.Pts.size()) {
          for(unsigned short ipt = newSize; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
          SetEndPoints(tjs, tj);
          tj.AlgMod[kCTStepChk] = true;
          tj.Pts.resize(newSize);
          return;
        } // newSize < tj.Pts.size()
      } // tjs.UseAlg[kCTStepChk]
    } // isSA
    
    FindSoftKink(tj);
    
    HiEndDelta(tj);
    
    CheckHiMultUnusedHits(tj);
    if(!fGoodTraj || fQuitAlg) return;
    
    // lop off high multiplicity hits at the end
    CheckHiMultEndHits(tj);
    
    // Check for a Bragg peak at both ends. This may be used by FixTrajBegin.
    ChkStop(tj);

    if(prt && tj.Pts.size() < 100) PrintTrajectory("CTo", tjs, tj, USHRT_MAX);
    
  } // CheckTraj

  //////////////////////////////////////////
  void TrajClusterAlg::FindSoftKink(Trajectory& tj)
  {
    // Looks for a soft kink in the trajectory and truncates it if one is found.
    // This is best done after FixTrajBegin has been called.
    
    if(!tjs.UseAlg[kSoftKink]) return;
    if(tj.Pts.size() < 15) return;
    if(tj.MCSMom < 100) return;
    
    float dang = DeltaAngle(tj.Pts[tj.EndPt[0]].Ang, tj.Pts[tj.EndPt[1]].Ang);
    
    if(prt) {
      mf::LogVerbatim("TC")<<"FindSoftKink: "<<tj.ID<<" dang "<<dang<<" cut "<<0.5 * fKinkCuts[0];
    }
    if(dang < 0.5 * fKinkCuts[0]) return;
    // require at least 5 points fitted at the end of the trajectory
    unsigned short endPt = tj.EndPt[1];
    if(tj.Pts[endPt].NTPsFit < 5) return;
    if(tj.Pts[endPt].NTPsFit > endPt) return;
    // Estimate where where the kink would be
    unsigned short kinkPt = endPt - tj.Pts[endPt].NTPsFit;
    // Require at least 5 points in the trajectory before the kink
    if(prt) mf::LogVerbatim("TC")<<" kinkPt "<<kinkPt<<" NTPsFit at kinkPt "<<tj.Pts[kinkPt].NTPsFit<<" max "<<0.5 * kinkPt;
    if(kinkPt < 5) return;
    // require fewer points fitted in this region compared the number of points prior to it
    if(tj.Pts[kinkPt].NTPsFit > 0.5 * kinkPt) return;
    // scan back until we find the maximum number of points fitted
    unsigned short maxPtsFit = tj.Pts[kinkPt].NTPsFit;
    unsigned short atPt = kinkPt;
    for(unsigned short ipt = kinkPt; kinkPt > tj.EndPt[0] + 5; --ipt) {
      if(tj.Pts[ipt].NTPsFit > maxPtsFit) {
        maxPtsFit = tj.Pts[ipt].NTPsFit;
        atPt = ipt;
      }
      // stop scanning when the max starts falling
      if(tj.Pts[ipt].NTPsFit < maxPtsFit) break;
      if(ipt == 0) break;
    } // ipt
    if(atPt < 5) return;
    // require the trajectory be straight before the kink - the section we are going to keep
    if(MCSMom(tjs, tj, tj.EndPt[0], atPt) < 500) return;
    // release the hits in TPs after this point
    for(unsigned short ipt = atPt; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
    // Truncate the trajectory at this point
    tj.Pts.resize(atPt + 1);
    SetEndPoints(tjs, tj);
    tj.AlgMod[kSoftKink] = true;
    if(prt) mf::LogVerbatim("TC")<<" truncated trajectory at "<<PrintPos(tjs, tj.Pts[tj.Pts.size()-1]);
    
  } // FindSoftKinks

  ////////////////////////////////////////////////
  void TrajClusterAlg::FixTrajBegin(Trajectory& tj)
  {
    // Update the parameters at the beginning of the trajectory. The first
    // points may not belong to this trajectory since they were added when there was
    // little information. This information may be updated later if ReversePropagate is used
    
    if(!tjs.UseAlg[kFixBegin]) return;
    
    // don't do anything if this tj has been modified by ReversePropagate
    if(tj.AlgMod[kRevProp]) return;
    
    // don't bother with really short tjs
    if(tj.Pts.size() < 3) return;
    
    unsigned short lastPtToChk = 10;
    if(tjs.UseAlg[kFTBRevProp]) lastPtToChk = tj.EndPt[1];

    unsigned short atPt = tj.EndPt[1];
    unsigned short maxPtsFit = 0;
    for(unsigned short ipt = 3; ipt < lastPtToChk; ++ipt) {
      if(ipt == tj.Pts.size()) break;
      if(tj.Pts[ipt].Chg == 0) continue;
      if(tj.Pts[ipt].NTPsFit >= maxPtsFit) {
        maxPtsFit = tj.Pts[ipt].NTPsFit;
        atPt = ipt;
        // no reason to continue if there are a good number of points fitted
        if(maxPtsFit > 20) break;
      }
    } // ipt
    // find the first point that is in this fit
    unsigned short firstPtFit = tj.EndPt[0];
    unsigned short cnt = 0;
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      if(ii > atPt) break;
      unsigned short ipt = atPt - ii;
      if(tj.Pts[ipt].Chg == 0) continue;
      ++cnt;
      if(cnt == maxPtsFit) {
        firstPtFit = ipt;
        break;
      } // full count
    } // ii
    
    bool needsRevProp = firstPtFit > 3;
    
    if(!needsRevProp) {
      // check one wire on the other side of EndPt[0] to see if there are hits that are available which could
      // be picked up by reverse propagation
      TrajPoint tp = tj.Pts[0];
      tp.Hits.clear();
      tp.UseHit.reset();
      // Move the TP "backwards"
      double stepSize = fVLAStepSize;
      if(tp.AngleCode < 2) stepSize = std::abs(1/tp.Dir[0]);
      tp.Pos[0] -= tp.Dir[0] * stepSize * tj.StepDir;
      tp.Pos[1] -= tp.Dir[1] * stepSize * tj.StepDir;
      float maxDelta = 3 * tp.DeltaRMS;
      if(FindCloseHits(tjs, tp, maxDelta, kUnusedHits) && !tp.Hits.empty()) {
        needsRevProp = true;
        if(prt) {
          mf::LogVerbatim("TC")<<"FTB: Close unused hits found near EndPt[0] "<<tp.Hits.size()<<" or dead wire. Call ReversePropagate";
          PrintTrajPoint("FTB", tjs, 0, tj.StepDir, tj.Pass, tp);
        }
      }
    } // !needsRevProp
    
    if(prt) mf::LogVerbatim("TC")<<"FTB: maxPtsFit "<<maxPtsFit<<" at point "<<atPt<<" firstPtFit "<<firstPtFit<<" Needs ReversePropagate? "<<needsRevProp;

    if(tjs.UseAlg[kFTBRevProp] && needsRevProp) {
      // lop off the points before firstPtFit and reverse propagate
      if(prt) mf::LogVerbatim("TC")<<"  FTB call ReversePropagate ";
      for(unsigned short ipt = 0; ipt < firstPtFit; ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
      SetEndPoints(tjs, tj);
      tj.AlgMod[kFTBRevProp] = true;
      // update the trajectory 
      for(unsigned short ipt = tj.EndPt[0]; ipt < atPt; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        tp.Dir = tj.Pts[atPt].Dir;
        tp.Ang = tj.Pts[atPt].Ang;
        tp.AngErr = tj.Pts[atPt].AngErr;
        tp.AngleCode = tj.Pts[atPt].AngleCode;
        // Correct the projected time to the wire
        float dw = tp.Pos[0] - tj.Pts[atPt].Pos[0];
        if(tp.Dir[0] != 0) tp.Pos[1] = tj.Pts[atPt].Pos[1] + dw * tp.Dir[1] / tp.Dir[0];
        tp.Delta = PointTrajDOCA(tjs, tp.HitPos[0], tp.HitPos[1], tp);
        tp.DeltaRMS = tj.Pts[atPt].DeltaRMS;
        tp.NTPsFit = tj.Pts[atPt].NTPsFit;
        tp.FitChi = tj.Pts[atPt].FitChi;
        tp.AveChg = tj.Pts[firstPtFit].AveChg;
        tp.ChgPull = (tj.Pts[ipt].Chg / tj.AveChg - 1) / tj.ChgRMS;
        if(prt) PrintTrajectory("ftbPrep", tjs, tj, ipt);
      } // ii
      // Check for quality and trim if necessary
      TrimEndPts(tjs, tj, fQualityCuts, prt);
      if(tj.AlgMod[kKilled]) {
        fGoodTraj = false;
        return;
      }
      ReversePropagate(tj);
    } else if(firstPtFit > 0) {
      FixTrajBegin(tj, firstPtFit);
    } else {
      // The first points were in the fit but the angle may not be well defined
      for(unsigned short ipt = tj.EndPt[0]; ipt < atPt; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        tp.Dir = tj.Pts[atPt].Dir;
        tp.Ang = tj.Pts[atPt].Ang;
        tp.AngErr = tj.Pts[atPt].AngErr;
        tp.AngleCode = tj.Pts[atPt].AngleCode;
      } // ipt
    }

  } // FixTrajBegin
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FixTrajBegin(Trajectory& tj, unsigned short atPt)
  {
    // Update the parameters at the beginning of the trajectory starting at point atPt
    
    if(!tjs.UseAlg[kFixBegin]) return;
    // ignore short trajectories
    unsigned short npwc = NumPtsWithCharge(tjs, tj, false);
    if(npwc < 6) return;
    // ignore somewhat longer trajectories that are curly
    if(npwc < 10 && tj.MCSMom < 100) return;
    // ignore shower-like trajectories
    if(tj.PDGCode == 11) return;
    // ignore junk trajectories
    if(tj.AlgMod[kJunkTj]) return;
    // ignore stopping trajectories
    if(tj.StopFlag[0][kBragg]) return;
    
    
    unsigned short firstPt = tj.EndPt[0];
    if(prt) {
      mf::LogVerbatim("TC")<<"FixTrajBegin: atPt "<<atPt<<" firstPt "<<firstPt<<" Stops at end 0? "<<PrintStopFlag(tj, 0);
    }
    
    if(atPt == tj.EndPt[0]) return;
    
    float maxDelta = 4 * tj.Pts[tj.EndPt[1]].DeltaRMS;
    
    // update the trajectory for all the points up to atPt
    // assume that we will use all of these points
    bool maskPts = false;
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      if(ii > atPt) break;
      unsigned int ipt = atPt - ii;
      TrajPoint& tp = tj.Pts[ipt];
      tp.Dir = tj.Pts[atPt].Dir;
      tp.Ang = tj.Pts[atPt].Ang;
      tp.AngErr = tj.Pts[atPt].AngErr;
      tp.AngleCode = tj.Pts[atPt].AngleCode;
      // Correct the projected time to the wire
      float dw = tp.Pos[0] - tj.Pts[atPt].Pos[0];
      if(tp.Dir[0] != 0) tp.Pos[1] = tj.Pts[atPt].Pos[1] + dw * tp.Dir[1] / tp.Dir[0];
      bool newHits = false;
      tj.Pts[ipt].Delta = PointTrajDOCA(tjs, tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tj.Pts[ipt]);
      tj.Pts[ipt].DeltaRMS = tj.Pts[atPt].DeltaRMS;
      tj.Pts[ipt].NTPsFit = tj.Pts[atPt].NTPsFit;
      tj.Pts[ipt].FitChi = tj.Pts[atPt].FitChi;
      tj.Pts[ipt].AveChg = tj.Pts[atPt].AveChg;
      tj.Pts[ipt].ChgPull = (tj.Pts[ipt].Chg / tj.AveChg - 1) / tj.ChgRMS;
      if(tj.Pts[ipt].Delta > maxDelta) maskPts = true;
      if(maskPts) UnsetUsedHits(tjs, tp);
      if(prt) {
        if(newHits) {
          PrintTrajectory("FTB", tjs, tj, ipt);
        } else {
          PrintTrajectory("ftb", tjs, tj, ipt);
        }
      }
      if(ipt == 0) break;
    } // ii
    if(maskPts) SetEndPoints(tjs, tj);
    tj.AlgMod[kFixBegin] = true;
    
  } // FixTrajBegin
  
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FixTrajEnd(Trajectory& tj, unsigned short atPt)
  {
    // Update the parameters at the end of the trajectory starting at point atPt
    
    if(!tjs.UseAlg[kFixEnd]) return;
    // ignore short trajectories
    unsigned short npwc = NumPtsWithCharge(tjs, tj, false);
    if(npwc < 6) return;
    // ignore somewhat longer trajectories that are curly
    if(npwc < 10 && tj.MCSMom < 100) return;
    // ignore shower-like trajectories
    if(tj.PDGCode == 11) return;
    // ignore junk trajectories
    if(tj.AlgMod[kJunkTj]) return;
    // ingore stopping trajectories
    if(tj.StopFlag[1][kBragg]) return;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"FixTrajEnd: atPt "<<atPt;
    }
    
    if(atPt == tj.EndPt[1]) return;

    // update the trajectory for all the intervening points
    for(unsigned short ipt = atPt + 1; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      tp.Dir = tj.Pts[atPt].Dir;
      tp.Ang = tj.Pts[atPt].Ang;
      tp.AngErr = tj.Pts[atPt].AngErr;
      tp.AngleCode = tj.Pts[atPt].AngleCode;
      // Correct the projected time to the wire
      float dw = tp.Pos[0] - tj.Pts[atPt].Pos[0];
      if(tp.Dir[0] != 0) tp.Pos[1] = tj.Pts[atPt].Pos[1] + dw * tp.Dir[1] / tp.Dir[0];
      tj.Pts[ipt].Delta = PointTrajDOCA(tjs, tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tj.Pts[ipt]);
      tj.Pts[ipt].DeltaRMS = tj.Pts[atPt].DeltaRMS;
      tj.Pts[ipt].NTPsFit = tj.Pts[atPt].NTPsFit;
      tj.Pts[ipt].FitChi = tj.Pts[atPt].FitChi;
      tj.Pts[ipt].AveChg = tj.Pts[atPt].AveChg;
      tj.Pts[ipt].ChgPull = (tj.Pts[ipt].Chg / tj.AveChg - 1) / tj.ChgRMS;
      if(prt) {
        PrintTrajectory("FTE", tjs, tj, ipt);
      }
    } // ipt
    tj.AlgMod[kFixEnd] = true;
    
  } // FixTrajEnd
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FillGaps(Trajectory& tj)
  {
    // Fill in any gaps in the trajectory with close hits regardless of charge (well maybe not quite that)
   
    if(!tjs.UseAlg[kFillGap]) return;
    
    // start with the first point that has charge
    short firstPtWithChg = tj.EndPt[0];
    bool first = true;
    float maxDelta = 1;
    // don't let MCSMom suffer too much while filling gaps
    short minMCSMom = 0.9 * tj.MCSMom;
    while(firstPtWithChg < tj.EndPt[1]) {
      short nextPtWithChg = firstPtWithChg + 1;
      // find the next point with charge
      for(nextPtWithChg = firstPtWithChg + 1; nextPtWithChg < tj.EndPt[1]; ++nextPtWithChg) {
        if(tj.Pts[nextPtWithChg].Chg > 0) break;
      } // nextPtWithChg
      if(nextPtWithChg == firstPtWithChg + 1) {
        // the next point has charge
        ++firstPtWithChg;
        continue;
      }
      // Found a gap. Make a bare trajectory point at firstPtWithChg that points to nextPtWithChg
      TrajPoint tp;
      if(!MakeBareTrajPoint(tjs, tj.Pts[firstPtWithChg], tj.Pts[nextPtWithChg], tp)) {
        mf::LogVerbatim("TC")<<"FillGaps: Failure from MakeBareTrajPoint ";
        PrintTrajectory("FG", tjs, tj, USHRT_MAX);
        fGoodTraj = false;
        return;
      }
      // Find the maximum delta between hits and the trajectory Pos for all
      // hits on this trajectory
      if(first) {
        maxDelta = 2.5 * MaxHitDelta(tjs, tj);
        first = false;
      } // first
      // define a loose charge cut using the average charge at the first point with charge
      float maxChg = tj.Pts[firstPtWithChg].AveChg * (1 + 2 * fChargeCuts[0] * tj.ChgRMS);
      // Eliminate the charge cut altogether if we are close to an end
      if(tj.Pts.size() < 10) {
        maxChg = 1E6;
      } else {
        short chgCutPt = tj.EndPt[0] + 5;
        if(firstPtWithChg < chgCutPt) {
          // gap is near end 0
          maxChg = 1E6;
        } else {
          // check for gap near end 1
          chgCutPt = tj.EndPt[1] - 5;
          if(chgCutPt < tj.EndPt[0]) chgCutPt = tj.EndPt[0];
          if(nextPtWithChg > chgCutPt) maxChg = 1E6;
        }
      }

      // fill in the gap
      for(unsigned short mpt = firstPtWithChg + 1; mpt < nextPtWithChg; ++mpt) {
        if(tj.Pts[mpt].Chg > 0) {
          mf::LogWarning("TC")<<"FillGaps coding error: firstPtWithChg "<<firstPtWithChg<<" mpt "<<mpt<<" nextPtWithChg "<<nextPtWithChg;
          fQuitAlg = true;
          return;
        }
        bool filled = false;
        float chg = 0;
        for(unsigned short ii = 0; ii < tj.Pts[mpt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[mpt].Hits[ii];
          if(tjs.fHits[iht].InTraj > 0) continue;
          float delta = PointTrajDOCA(tjs, iht, tp);
          if(prt) mf::LogVerbatim("TC")<<" FG "<<PrintPos(tjs,tj.Pts[mpt])<<" hit "<<PrintHit(tjs.fHits[iht])<<" delta "<<delta<<" maxDelta "<<maxDelta<<" Chg "<<tjs.fHits[iht].Integral<<" maxChg "<<maxChg;
          if(delta > maxDelta) continue;
          if(tj.Pts[mpt].UseHit[ii]) {
            mf::LogWarning("TC")<<"FillGaps: Found UseHit true on TP with no charge "<<tj.ID<<" mpt "<<mpt<<" hit "<<PrintHit(tjs.fHits[iht]);
            fQuitAlg = true;
            return;
          }
          tj.Pts[mpt].UseHit[ii] = true;
          tjs.fHits[iht].InTraj = tj.ID;
          chg += tjs.fHits[iht].Integral;
          filled = true;
        } // ii
        if(chg > maxChg || MCSMom(tjs, tj) < minMCSMom) {
          // don't use these hits after all
          UnsetUsedHits(tjs, tj.Pts[mpt]);
          filled = false;
        }
        if(filled) {
          DefineHitPos(tj.Pts[mpt]);
          tj.AlgMod[kFillGap] = true;
          if(prt) {
            PrintTrajPoint("FG", tjs, mpt, tj.StepDir, tj.Pass, tj.Pts[mpt]);
            mf::LogVerbatim("TC")<<"Check MCSMom "<<MCSMom(tjs, tj);
          }
        } // filled
      } // mpt
      firstPtWithChg = nextPtWithChg;
    } // firstPtWithChg
    
    if(tj.AlgMod[kFillGap]) tj.MCSMom = MCSMom(tjs, tj);
    
  } // FillGaps 
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::HiEndDelta(Trajectory& tj)
  {
    // Modify the trajectory at the end if there is a consistent increase in delta. It
    // is called from CheckTraj.
    // This needs to be done carefully...
    
    if(!tjs.UseAlg[kHED]) return;
    if(tj.StopFlag[1][kBragg]) return;
    // Only consider long high momentum.
    if(tj.MCSMom < 100) return;
    if(tj.Pts.size() < 100) return;

    unsigned short ept = tj.EndPt[1];

    TrajPoint& lastTp = tj.Pts[ept];

    if(lastTp.AngleCode > 1) return;
    if(lastTp.FitChi < 1) return;
    
    unsigned short npts = USHRT_MAX;
    float lastDelta = lastTp.Delta;
    // check the last 20 points on the trajectory for a systematic increase in Delta and FitChi
    for(unsigned short ii = 1; ii < 20; ++ii) {
      unsigned short ipt = ept - ii;
      TrajPoint& tp = tj.Pts[ipt];
      if(prt) std::cout<<ii<<" Chk "<<ipt<<" "<<tp.FitChi<<" "<<tp.Delta<<"\n";
      if(tp.Chg == 0) continue;
      if(tp.FitChi < 1 || tp.Delta > lastDelta) {
        npts = ii;
        break;
      }
      lastDelta = tp.Delta;
    } // ii
    
    if(prt) mf::LogVerbatim("TC")<<"HED: last point FitChi "<<lastTp.FitChi<<" NTPsFit "<<lastTp.NTPsFit<<" new npts "<<npts;
    
    // something bad happened
    if(npts == USHRT_MAX) return;
    // The Tj end has some other problem
    if(npts < 4) return;
    
    // re-fit the end of the trajectory yyy
    lastTp.NTPsFit = npts;
    FitTraj(tjs, tj);
    if(prt) PrintTrajPoint("HED", tjs, ept, tj.StepDir, tj.Pass, lastTp);
    // update the last points
    for(unsigned short ii = 1; ii <= npts; ++ii) {
      unsigned short ipt = ept - ii;
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      tp.Dir = tj.Pts[ept].Dir;
      tp.Ang = tj.Pts[ept].Ang;
      tp.AngErr = tj.Pts[ept].AngErr;
      tp.AngleCode = tj.Pts[ept].AngleCode;
      // Correct the projected time to the wire
      float dw = tp.Pos[0] - tj.Pts[ept].Pos[0];
      if(tp.Dir[0] != 0) tp.Pos[1] = tj.Pts[ept].Pos[1] + dw * tp.Dir[1] / tp.Dir[0];
      tp.Delta = PointTrajDOCA(tjs, tp.HitPos[0], tp.HitPos[1], tp);
      tp.DeltaRMS = tj.Pts[ept].DeltaRMS;
      tp.NTPsFit = tj.Pts[ept].NTPsFit;
      tp.FitChi = tj.Pts[ept].FitChi;
      if(prt) PrintTrajPoint("HED", tjs, ipt, tj.StepDir, tj.Pass, tp);
    } // ii

    tj.AlgMod[kHED] = true;
    
  } // HiEndDelta
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiMultUnusedHits(Trajectory& tj)
  {
    // Check for many unused hits in high multiplicity TPs in work and try to use them
    
    if(!tjs.UseAlg[kChkHiMultHits]) return;
    
    // This code might do bad things to short trajectories
    if(NumPtsWithCharge(tjs, tj, true) < 6) return;
    if(tj.EndPt[0] == tj.EndPt[1]) return;
    
    // count the number of unused hits multiplicity > 1 hits and decide
    // if the unused hits should be used. This may trigger another
    // call to StepCrawl
    unsigned short ii, stopPt;
    // Use this to see if the high multiplicity Pts are mostly near
    // the end or further upstream
    unsigned short lastMult1Pt = USHRT_MAX;
    // the number of TPs with > 1 hit (HiMult)
    unsigned short nHiMultPt = 0;
    // the total number of hits associated with HiMult TPs
    unsigned short nHiMultPtHits = 0;
    // the total number of used hits associated with HiMult TPs
    unsigned short nHiMultPtUsedHits = 0;
    unsigned int iht;
    // start counting at the leading edge and break if a hit
    // is found that is used in a trajectory
    bool doBreak = false;
    unsigned short jj;
    for(ii = 1; ii < tj.Pts.size(); ++ii) {
      stopPt = tj.EndPt[1] - ii;
      for(jj = 0; jj < tj.Pts[stopPt].Hits.size(); ++jj) {
        iht = tj.Pts[stopPt].Hits[jj];
        if(tjs.fHits[iht].InTraj > 0) {
          doBreak = true;
          break;
        }
      } // jj
      if(doBreak) break;
      // require 2 consecutive multiplicity = 1 points
      if(lastMult1Pt == USHRT_MAX && tj.Pts[stopPt].Hits.size() == 1 && tj.Pts[stopPt-1].Hits.size() == 1) lastMult1Pt = stopPt;
      if(tj.Pts[stopPt].Hits.size() > 1) {
        ++nHiMultPt;
        nHiMultPtHits += tj.Pts[stopPt].Hits.size();
        nHiMultPtUsedHits += NumHitsInTP(tj.Pts[stopPt], kUsedHits);
      } // high multiplicity TP
      // stop looking when two consecutive single multiplicity TPs are found
      if(lastMult1Pt != USHRT_MAX) break;
      if(stopPt == 1) break;
    } // ii
    // Don't do this if there aren't a lot of high multiplicity TPs
    float fracHiMult = (float)nHiMultPt / (float)ii;
    if(lastMult1Pt != USHRT_MAX) {
      float nchk = tj.EndPt[1] - lastMult1Pt + 1;
      fracHiMult = (float)nHiMultPt / nchk;
    } else {
      fracHiMult = (float)nHiMultPt / (float)ii;
    }
    float fracHitsUsed = 0;
    if(nHiMultPt > 0 && nHiMultPtHits > 0) fracHitsUsed = (float)nHiMultPtUsedHits / (float)nHiMultPtHits;
    // Use this to limit the number of points fit for trajectories that
    // are close the LA tracking cut
    ii = tj.EndPt[1];
    bool sortaLargeAngle = (AngleRange(tjs, tj.Pts[ii]) == 1);

    if(prt) mf::LogVerbatim("TC")<<"CHMUH: First InTraj stopPt "<<stopPt<<" fracHiMult "<<fracHiMult<<" fracHitsUsed "<<fracHitsUsed<<" lastMult1Pt "<<lastMult1Pt<<" sortaLargeAngle "<<sortaLargeAngle;
    if(fracHiMult < 0.3) return;
    if(fracHitsUsed > 0.98) return;
    
    float maxDelta = 2.5 * MaxHitDelta(tjs, tj);

    if(prt) {
      mf::LogVerbatim("TC")<<" Pts size "<<tj.Pts.size()<<" nHiMultPt "<<nHiMultPt<<" nHiMultPtHits "<<nHiMultPtHits<<" nHiMultPtUsedHits "<<nHiMultPtUsedHits<<" sortaLargeAngle "<<sortaLargeAngle<<" maxHitDelta "<<maxDelta;
    }

    // Use next pass cuts if available
    if(sortaLargeAngle && tj.Pass < fMinPtsFit.size()-1) ++tj.Pass;

    // Make a copy of tj in case something bad happens
    Trajectory TjCopy = tj;
    // and the list of used hits
    auto inTrajHits = PutTrajHitsInVector(tj, kUsedHits);
    unsigned short ipt;

    // unset the used hits from stopPt + 1 to the end
    for(ipt = stopPt + 1; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
    SetEndPoints(tjs, tj);
    unsigned short killPts;
    float delta;
    bool added;
    for(ipt = stopPt + 1; ipt < tj.Pts.size(); ++ipt) {
      // add hits that are within maxDelta and re-fit at each point
      added = false;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        iht = tj.Pts[ipt].Hits[ii];
        if(prt) mf::LogVerbatim("TC")<<" ipt "<<ipt<<" hit "<<PrintHit(tjs.fHits[iht])<<" inTraj "<<tjs.fHits[iht].InTraj<<" delta "<<PointTrajDOCA(tjs, iht, tj.Pts[ipt]);
        if(tjs.fHits[iht].InTraj != 0) continue;
        delta = PointTrajDOCA(tjs, iht, tj.Pts[ipt]);
        if(delta > maxDelta) continue;
        if (!NumHitsInTP(TjCopy.Pts[ipt], kUsedHits)||TjCopy.Pts[ipt].UseHit[ii]){
          tj.Pts[ipt].UseHit[ii] = true;
          tjs.fHits[iht].InTraj = tj.ID;
          added = true;
        }
      } // ii
      if(added) DefineHitPos(tj.Pts[ipt]);
      if(tj.Pts[ipt].Chg == 0) continue;
      tj.EndPt[1] = ipt;
      // This will be incremented by one in UpdateTraj
      if(sortaLargeAngle) tj.Pts[ipt].NTPsFit = 2;
      UpdateTraj(tj);
      if(!fUpdateTrajOK) {
        if(prt) mf::LogVerbatim("TC")<<"UpdateTraj failed on point "<<ipt;
        // Clobber the used hits from the corrupted points in tj
        for(unsigned short jpt = stopPt + 1; jpt <= ipt; ++jpt) {
          for(unsigned short jj = 0; jj < tj.Pts[jpt].Hits.size(); ++jj) {
            if(tj.Pts[jpt].UseHit[jj]) tjs.fHits[tj.Pts[jpt].Hits[jj]].InTraj = 0;
          } // jj
        } // jpt
        // restore the original trajectory
        tj = TjCopy;
        // restore the hits
        for(unsigned short jpt = stopPt + 1; jpt <= ipt; ++jpt) {
          for(unsigned short jj = 0; jj < tj.Pts[jpt].Hits.size(); ++jj) {
            if(tj.Pts[jpt].UseHit[jj]) tjs.fHits[tj.Pts[jpt].Hits[jj]].InTraj = tj.ID;
          } // jj
        } // jpt
        return;
      }
      GottaKink(tj, killPts);
      if(killPts > 0) {
        MaskTrajEndPoints(tj, killPts);
        if(!fGoodTraj) return;
        break;
      }
      if(prt) PrintTrajectory("CHMUH", tjs, tj, ipt);
    } // ipt
    // if we made it here it must be OK
    SetEndPoints(tjs, tj);
    // Try to extend it, unless there was a kink
    if(tj.StopFlag[1][kAtKink]) return;
    // trim the end points although this shouldn't happen
    if(tj.EndPt[1] != tj.Pts.size() - 1) tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kChkHiMultHits] = true;
  } // CheckHiMultUnusedHits

  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiMultEndHits(Trajectory& tj)
  {
    // mask off high multiplicity TPs at the end
    if(!tjs.UseAlg[kChkHiMultEndHits]) return;
    if(tj.StopFlag[1][kBragg]) return;
    if(tj.Pts.size() < 10) return;
    // find the average multiplicity in the first half
    unsigned short aveMult= 0;
    unsigned short ipt, nhalf = tj.Pts.size() / 2;
    unsigned short cnt = 0;
    for(auto& tp : tj.Pts) {
      if(tp.Chg == 0) continue;
      aveMult += tp.Hits.size();
      ++cnt;
      if(cnt == nhalf) break;
    } //  pt
    if(cnt == 0) return;
    aveMult /= cnt;
    if(aveMult == 0) aveMult = 1;
    // convert this into a cut
    aveMult *= 3;
    cnt = 0;
    for(ipt = tj.EndPt[1]; ipt > tj.EndPt[0]; --ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      if(tj.Pts[ipt].Hits.size() > aveMult) {
        UnsetUsedHits(tjs, tj.Pts[ipt]);
        ++cnt;
        continue;
      }
      break;
    } // ipt
    if(prt) mf::LogVerbatim("TC")<<"CHMEH multiplicity cut "<<aveMult<<" number of TPs masked off "<<cnt;
    if(cnt > 0) {
      tj.AlgMod[kChkHiMultEndHits] = true;
      SetEndPoints(tjs, tj);
    }
  } // CheckHiMultEndHits

  ////////////////////////////////////////////////
  bool TrajClusterAlg::StopIfBadFits(Trajectory& tj)
  {
    // Returns true if there are a number of Tps that were not used in the trajectory because the fit was poor and the
    // charge pull is not really high. This 
    
    // don't consider muons
    if(tj.PDGCode == 13) return false;
    // or long straight Tjs
    if(tj.Pts.size() > 40 && tj.MCSMom > 200) return false;
    
    unsigned short nBadFit = 0;
    unsigned short nHiChg = 0;
    unsigned short cnt = 0;
    for(unsigned short ipt = tj.Pts.size() - 1; ipt > tj.EndPt[1]; --ipt ) {
      if(tj.Pts[ipt].FitChi > 2) ++nBadFit;
      if(tj.Pts[ipt].ChgPull > 3) ++nHiChg;
      ++cnt;
      if(cnt == 5) break;
    } // ipt
    
    if(prt) mf::LogVerbatim("TC")<<"StopIfBadFits: nBadFit "<<nBadFit<<" nHiChg "<<nHiChg;
    if(nBadFit > 3 && nHiChg == 0) return true;
    return false;
    
  } // StopIfBadFits

  ////////////////////////////////////////////////
  bool TrajClusterAlg::MaskedHitsOK(Trajectory& tj)
  {
    // Version 2 of MaskedHitsOK.
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration (true) or whether to stop and possibly try with the next pass settings (false)
    
    unsigned short lastPt = tj.Pts.size() - 1;
    if(tj.Pts[lastPt].Chg > 0) return true;
    unsigned short endPt = tj.EndPt[1];
    
    // count the number of points w/o used hits and the number with one unused hit
    unsigned short nMasked = 0;
    unsigned short nOneHit = 0;
    unsigned short nOKChg = 0;
    unsigned short nOKDelta = 0;
    float maxOKDelta = 10 * tj.Pts[endPt].DeltaRMS;
    float maxOKChg = 0;
    // find the maximum charge point on the trajectory
    // TODO: Might have to be more careful if we use 1 WSEU normalized charge
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) if(tj.Pts[ipt].Chg > maxOKChg) maxOKChg = tj.Pts[ipt].Chg;
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - ii;
      if(tj.Pts[ipt].Chg > 0) break;
      unsigned short nUnusedHits = 0;
      float chg = 0;
      for(unsigned short jj = 0; jj < tj.Pts[ipt].Hits.size(); ++jj) {
        unsigned int iht = tj.Pts[ipt].Hits[jj];
        if(tjs.fHits[iht].InTraj != 0) continue;
        ++nUnusedHits;
        chg += tjs.fHits[iht].Integral;
      } // jj
      if(chg < maxOKChg) ++nOKChg;
      if(nUnusedHits == 1) ++nOneHit;
      if(tj.Pts[ipt].Delta < maxOKDelta) ++nOKDelta;
      ++nMasked;
    } // ii
    
    if(prt) {
      mf::LogVerbatim("TC")<<"MaskedHitsOK2:  nMasked "<<nMasked<<" nOneHit "<<nOneHit<<" nOKChg "<<nOKChg<<" nOKDelta "<<nOKDelta;
    }

    if(nMasked < 10 || nOneHit < 10) return true;
    
    if(nOKDelta != nMasked) return true;
    if(nOKChg != nMasked) return true;
    // we would like to reduce the number of fitted points to a minimum and include
    // the masked hits, but we can only do that if there are enough points
    if(tj.Pts[endPt].NTPsFit <= fMinPtsFit[tj.Pass]) {
      // stop stepping if we have masked off more points than are in the fit
      if(nMasked > tj.Pts[endPt].NTPsFit) return false;
      return true;
    }
    // Reduce the number of points fit and try to include the points
    unsigned short newNTPSFit = tj.Pts[endPt].NTPsFit / 2;
    if(newNTPSFit < fMinPtsFit[tj.Pass]) newNTPSFit = fMinPtsFit[tj.Pass];
    for(unsigned ipt = endPt + 1; ipt < tj.Pts.size(); ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(tjs.fHits[iht].InTraj == 0) {
          tp.UseHit[ii] = true;
          tjs.fHits[iht].InTraj = tj.ID;
          break;
        }
      } // ii
      DefineHitPos(tp);
      SetEndPoints(tjs, tj);
      tp.NTPsFit = newNTPSFit;
      FitTraj(tjs, tj);
      if(prt) PrintTrajectory("MHOK2", tjs, tj, ipt);
    } // ipt
    
    tj.AlgMod[kMaskHits] = true;
    UpdateAveChg(tj);
    return true;
    
  } // MaskedHitsOK

  ////////////////////////////////////////////////
  void TrajClusterAlg::PrepareForNextPass(Trajectory& tj)
  {
    // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
    // fitted points to get FitCHi < 2
    
    fTryWithNextPass = false;

    // See if there is another pass available
    if(tj.Pass > fMinPtsFit.size()-2) return;
    ++tj.Pass;
    
    unsigned short lastPt = tj.Pts.size() - 1;
    // Return if the last fit chisq is OK
    if(tj.Pts[lastPt].FitChi < 1.5) {
      fTryWithNextPass = true;
      return;
    }
    TrajPoint& lastTP = tj.Pts[lastPt];
    unsigned short newNTPSFit = lastTP.NTPsFit;
    // only give it a few tries before giving up
    unsigned short nit = 0;

     while(lastTP.FitChi > 1.5 && lastTP.NTPsFit > 2) {
      if(lastTP.NTPsFit > 3) newNTPSFit -= 2;
      else if(lastTP.NTPsFit == 3) newNTPSFit = 2;
      lastTP.NTPsFit = newNTPSFit;
      FitTraj(tjs, tj);
      if(prt) mf::LogVerbatim("TC")<<"PrepareForNextPass: FitChi is > 1.5 "<<lastTP.FitChi<<" Reduced NTPsFit to "<<lastTP.NTPsFit<<" tj.Pass "<<tj.Pass;
      if(lastTP.NTPsFit <= fMinPtsFit[tj.Pass]) break;
      ++nit;
      if(nit == 3) break;
    }
    // decide if the next pass should indeed be attempted
    if(lastTP.FitChi > 2) return;
    fTryWithNextPass = true;
    
  } // PrepareForNextPass
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::GottaKink(Trajectory& tj, unsigned short& killPts)
  {
    // Checks the last few points on the trajectory and returns with the number of
    // points (killPts) that should be killed (aka masked) at the end
    // fKinkCuts
    // 0 = kink angle cut (radians)
    // 1 = kink angle significance cut
    // 2 = nPts fit at the end of the tj
    // Kink angle cut = fKinkCuts[0] + fKinkCuts[1] * MCSThetaRMS
    
    killPts = 0;
    
    // decide whether to turn kink checking back on
    if(fKinkCuts[0] > 0 && tj.EndPt[1] == 20) {
      if(MCSMom(tjs, tj, 10, 19) > 50) tj.AlgMod[kNoKinkChk] = false;
      if(prt) mf::LogVerbatim("TC")<<"GottaKink turn kink checking back on? "<<tj.AlgMod[kNoKinkChk]<<" with MCSMom "<<MCSMom(tjs, tj, 10, 19);
    }
    if(tj.AlgMod[kNoKinkChk]) return;

    unsigned short lastPt = tj.EndPt[1];
    if(lastPt < 5) return;
    if(tj.Pts[lastPt].Chg == 0) return;
    
    // MCSThetaRMS is the scattering angle for the entire length of the trajectory. Convert
    // this to the scattering angle for one WSE unit
    float thetaRMS = MCSThetaRMS(tjs, tj, tj.EndPt[0], tj.EndPt[1]) / sqrt(TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[lastPt]));
    float kinkAngCut = fKinkCuts[0] + fKinkCuts[1] * thetaRMS;
    // relax this a bit when doing RevProp
    if(tj.AlgMod[kRevProp]) kinkAngCut *= 1.3;
    
    // A simple check when there are few points being fit and the TJ is short. MCSMom isn't well known at this point so don't use it
    if(tj.Pts[lastPt].NTPsFit < 6 && tj.Pts.size() < 20) {
      unsigned short ii, prevPtWithHits = USHRT_MAX;
      unsigned short ipt;
      for(ii = 1; ii < tj.Pts.size(); ++ii) {
        ipt = lastPt - ii;
        if(tj.Pts[ipt].Chg > 0) {
          prevPtWithHits = ipt;
          break;
        }
        if(ipt == 0) break;
      } // ii
      if(prevPtWithHits == USHRT_MAX) return;
      float dang = DeltaAngle(tj.Pts[lastPt].Ang, tj.Pts[prevPtWithHits].Ang);
      kinkAngCut = 1.2 * fKinkCuts[0];
      if(prt) mf::LogVerbatim("TC")<<"GottaKink Simple check lastPt "<<PrintPos(tjs,tj.Pts[lastPt])<<" dang "<<dang<<" cut "<<kinkAngCut;
      if(dang > fKinkCuts[0]) {
        killPts = 1;
        tj.StopFlag[1][kAtKink] = true;
      }
      // Another case where there are few hits fit just prior to a dead wire
      // section or there were no hits added for several steps or due to a large
      // value of fMaxWireSkipNoSignal. We just added a bogus hit just after this section
      // so the trajectory angle change will be small. Find the angle between the previous
      // point fitted angle and the angle formed by the last two TPs
      if(std::abs(tj.Pts[lastPt-1].Pos[0] - tj.Pts[lastPt].Pos[0]) > 3) {
        TrajPoint tmp;
        if(!MakeBareTrajPoint(tjs, tj.Pts[lastPt-1], tj.Pts[lastPt], tmp)) {
          mf::LogVerbatim("TC")<<"GottaKink failure from MakeBareTrajPoint ";
          PrintTrajectory("GK", tjs, tj, USHRT_MAX);
          fGoodTraj = false;
          return;
        }
        dang = DeltaAngle(tmp.Ang, tj.Pts[prevPtWithHits].Ang);
        if(prt) mf::LogVerbatim("TC")<<"GottaKink Simple check after gap lastPt "<<lastPt<<" prevPtWithHits "<<prevPtWithHits<<" dang "<<dang<<" cut "<<kinkAngCut;
        if(dang > 1.5 * kinkAngCut) {
          killPts = 1;
          tj.StopFlag[1][kAtKink] = true;
        }
      }
      return;
    } // tj.Pts[lastPt].NTPsFit < 6 && tj.Pts.size() < 20

    if(tj.EndPt[1] < 10) return;
    
    unsigned short kinkPt = USHRT_MAX;
    
    // Find the kinkPt which is fKinkCuts[2] from the end that has charge
    unsigned short cnt = 0;
    unsigned short nPtsFit = fKinkCuts[2];
    unsigned short nHiMultPt = 0;
    unsigned short nHiChg = 0;
    
    for(unsigned short ii = 1; ii < lastPt; ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg == 0) continue;
      ++cnt;
      if(tj.Pts[ipt].Hits.size() > 1) ++nHiMultPt;
      if(tj.Pts[ipt].ChgPull > 1.5) ++nHiChg;
      if(cnt == nPtsFit) {
        kinkPt = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(kinkPt == USHRT_MAX) return;

    TrajPoint tpFit;
    unsigned short npts = 4;
    unsigned short fitDir = -1;
    FitTraj(tjs, tj, lastPt, npts, fitDir, tpFit);
    if(tpFit.FitChi > 1) return;
 
    float dang = DeltaAngle(tj.Pts[kinkPt].Ang, tpFit.Ang);
    
    if(dang > kinkAngCut) {
      killPts = nPtsFit;
      tj.StopFlag[1][kAtKink] = true;
    }
    
    if(killPts > 0) {
      // See if we are tracking a low momentum particle in which case we should just
      // turn off kink checking
      if(tjs.UseAlg[kNoKinkChk] && tj.EndPt[1] < 20) {
        // Find MCSMom if it hasn't been done
        if(tj.MCSMom < 0) tj.MCSMom = MCSMom(tjs, tj);
        if(tj.MCSMom < 50) {
          killPts = 0;
          tj.StopFlag[1][kAtKink] = false;
          tj.AlgMod[kNoKinkChk] = true;
          if(prt) mf::LogVerbatim("TC")<<"GottaKink turning off kink checking. MCSMom "<<tj.MCSMom;
        }
      } // turn off kink check
      // Don't stop if the last few points had high charge pull and we are tracking a muon, but do mask off the hits
      if(killPts > 0 && tj.PDGCode == 13 && tj.Pts[lastPt].ChgPull > 2  && tj.Pts[lastPt-1].ChgPull > 2) tj.StopFlag[1][kAtKink] = false;
      // Don't keep stepping or mask off any TPs if we hit a kink while doing RevProp
      if(tj.AlgMod[kRevProp]) killPts = 0;
    }
    
    if(prt) mf::LogVerbatim("TC")<<"GottaKink "<<kinkPt<<" Pos "<<PrintPos(tjs, tj.Pts[kinkPt])<<" dang "<<std::fixed<<std::setprecision(2)<<dang<<" cut "<<kinkAngCut<<" tpFit chi "<<tpFit.FitChi<<" killPts "<<killPts<<" GottaKink? "<<tj.StopFlag[1][kAtKink]<<" MCSMom "<<tj.MCSMom<<" thetaRMS "<<thetaRMS;
    
  } // GottaKink

  //////////////////////////////////////////
  void TrajClusterAlg::UpdateTraj(Trajectory& tj)
  {
    // Updates the last added trajectory point fit, average hit rms, etc.

    fUpdateTrajOK = false;
    fMaskedLastTP = false;
    
    if(tj.EndPt[1] < 1) return;
    unsigned int lastPt = tj.EndPt[1];
    TrajPoint& lastTP = tj.Pts[lastPt];

    // find the previous TP that has hits (and was therefore in the fit)
    unsigned short prevPtWithHits = USHRT_MAX;
    unsigned short firstFitPt = tj.EndPt[0];
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg > 0) {
        prevPtWithHits = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(prevPtWithHits == USHRT_MAX) return;
    
    // define the FitChi threshold above which something will be done
    float maxChi = 2;
    unsigned short minPtsFit = fMinPtsFit[tj.Pass];
    // just starting out?
    if(lastPt < 4) minPtsFit = 2;
    // was !TrajIsClean...
    if(tj.PDGCode == 13 && TrajIsClean(tjs, tj, prt)) {
      // Fitting a clean muon
      maxChi = fMaxChi;
      minPtsFit = lastPt / 3;
    }
    
    // Set the lastPT delta before doing the fit
    lastTP.Delta = PointTrajDOCA(tjs, lastTP.HitPos[0], lastTP.HitPos[1], lastTP);
    
    // update MCSMom. First ensure that nothing bad has happened
    float newMCSMom = MCSMom(tjs, tj);
    if(lastPt > 5 && newMCSMom < 0.6 * tj.MCSMom) {
      if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: MCSMom took a nose-dive "<<newMCSMom;
      UnsetUsedHits(tjs, lastTP);
      DefineHitPos(lastTP);
      SetEndPoints(tjs, tj);
      fUpdateTrajOK = true;
      return;
    }
    tj.MCSMom = newMCSMom;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"UpdateTraj: lastPt "<<lastPt<<" lastTP.Delta "<<lastTP.Delta<<" previous point with hits "<<prevPtWithHits<<" tj.Pts size "<<tj.Pts.size()<<" AngleCode "<<lastTP.AngleCode<<" PDGCode "<<tj.PDGCode<<" maxChi "<<maxChi<<" minPtsFit "<<minPtsFit<<" MCSMom "<<tj.MCSMom;
    }
    
    UpdateAveChg(tj);

    if(lastPt == 1) {
      // Handle the second trajectory point. No error calculation. Just update
      // the position and direction
      lastTP.NTPsFit = 2;
      FitTraj(tjs, tj);
      lastTP.FitChi = 0.01;
      lastTP.AngErr = tj.Pts[0].AngErr;
      if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: Second traj point pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1];
      fUpdateTrajOK = true;
      SetAngleCode(tjs, lastTP);
      return;
    }
    
    if(lastPt == 2) {
      // Third trajectory point. Keep it simple
      lastTP.NTPsFit = 3;
      FitTraj(tjs, tj);
      fUpdateTrajOK = true;
      if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: Third traj point fit "<<lastTP.FitChi;
      SetAngleCode(tjs, lastTP);
      return;
    }
    
    // Fit with > 2 TPs
    // Keep adding hits until Chi/DOF exceeds 1
    if(tj.Pts[prevPtWithHits].FitChi < 1) lastTP.NTPsFit += 1;
    // Reduce the number of points fit if the trajectory is long and chisq is getting a bit larger
    if(lastPt > 20 && tj.Pts[prevPtWithHits].FitChi > 1.5 && lastTP.NTPsFit > minPtsFit) lastTP.NTPsFit -= 2;
    if(tj.MCSMom < 100) {
      float localMCSMom = tj.MCSMom;
      if(NumPtsWithCharge(tjs, tj, false) > 10) {
        // long trajectory - only use the last 10 points
        localMCSMom = MCSMom(tjs, tj, lastPt - 10, lastPt);
        if(localMCSMom < 100) lastTP.NTPsFit = fMinPtsFit[tj.Pass];
      } else {
        // short trajectory
        lastTP.NTPsFit = fMinPtsFit[tj.Pass];
      }
      if(prt) mf::LogVerbatim("TC")<<" Low MCSMom? "<<localMCSMom<<" NTPsFit "<<lastTP.NTPsFit;
    } // tj.MCSMom < 100
    
    FitTraj(tjs, tj);
    
    // don't get too fancy when we are starting out
    if(lastPt < 6) {
      fUpdateTrajOK = true;
      UpdateDeltaRMS(tj);
      SetAngleCode(tjs, lastTP);
      if(prt) mf::LogVerbatim("TC")<<" Return with lastTP.FitChi "<<lastTP.FitChi<<" Chg "<<lastTP.Chg;
      return;
    }
    
    // find the first point that was fit.
    unsigned short cnt = 0;
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg > 0) {
        firstFitPt = ipt;
        ++cnt;
      }
      if(cnt == lastTP.NTPsFit) break;
      if(ipt == 0) break;
    }
    
    unsigned short ndead = DeadWireCount(tjs, lastTP.HitPos[0], tj.Pts[firstFitPt].HitPos[0], fCTP);
    
    if(lastTP.FitChi > 1.5 && tj.Pts.size() > 6) {
      // A large chisq jump can occur if we just jumped a large block of dead wires. In
      // this case we don't want to mask off the last TP but reduce the number of fitted points
      // This count will be off if there a lot of dead or missing wires...
      // reduce the number of points significantly
      if(ndead > 5) {
        if(lastTP.NTPsFit > 5) lastTP.NTPsFit = 5;
      } else {
        // Have a longish trajectory and chisq was a bit large.
        // Was this a sudden occurrence and the fraction of TPs are included
        // in the fit? If so, we should mask off this
        // TP and keep going. If these conditions aren't met, we
        // should reduce the number of fitted points
        float chirat = 0;
        if(prevPtWithHits != USHRT_MAX && tj.Pts[prevPtWithHits].FitChi > 0) chirat = lastTP.FitChi / tj.Pts[prevPtWithHits].FitChi;
        // Don't mask hits when doing RevProp. Reduce NTPSFit instead
        fMaskedLastTP = (chirat > 1.5 && lastTP.NTPsFit > 0.3 * NumPtsWithCharge(tjs, tj, false) && !tj.AlgMod[kRevProp]);
        if(prt) {
          mf::LogVerbatim("TC")<<" First fit chisq too large "<<lastTP.FitChi<<" prevPtWithHits chisq "<<tj.Pts[prevPtWithHits].FitChi<<" chirat "<<chirat<<" NumPtsWithCharge "<<NumPtsWithCharge(tjs, tj, false)<<" fMaskedLastTP "<<fMaskedLastTP;
        }
        // we should also mask off the last TP if there aren't enough hits
        // to satisfy the minPtsFit constraint
        if(!fMaskedLastTP && NumPtsWithCharge(tjs, tj, true) < minPtsFit) fMaskedLastTP = true;
      } // few dead wires
    } // lastTP.FitChi > 2 ...
    
    // Deal with a really long trajectory that is in trouble (uB cosmic).
    if(tj.PDGCode == 13 && lastTP.FitChi > fMaxChi) {
      if(lastTP.NTPsFit > 1.3 * tjs.MuonTag[0]) {
        lastTP.NTPsFit *= 0.8;
        if(prt) mf::LogVerbatim("TC")<<" Muon - Reduce NTPsFit "<<lastPt;
      } else {
        fMaskedLastTP = true;
        if(prt) mf::LogVerbatim("TC")<<" Muon - mask last point "<<lastPt;
      }
    }
    
    if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: First fit "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" FitChi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit<<" ndead wires "<<ndead<<" fMaskedLastTP "<<fMaskedLastTP;
    if(fMaskedLastTP) {
      UnsetUsedHits(tjs, lastTP);
      DefineHitPos(lastTP);
      SetEndPoints(tjs, tj);
      lastPt = tj.EndPt[1];
      lastTP.NTPsFit -= 1;
      FitTraj(tjs, tj);
      fUpdateTrajOK = true;
      SetAngleCode(tjs, lastTP);
      return;
    }  else {
      // a more gradual increase in chisq. Reduce the number of points
      unsigned short newNTPSFit = lastTP.NTPsFit;
      // reduce the number of points fit to keep Chisq/DOF < 2 adhering to the pass constraint
      // and also a minimum number of points fit requirement for long muons
      float prevChi = lastTP.FitChi;
      unsigned short ntry = 0;
      while(lastTP.FitChi > 1.5 && lastTP.NTPsFit > minPtsFit) {
        if(lastTP.NTPsFit > 15) {
          newNTPSFit = 0.7 * newNTPSFit;
        } else if(lastTP.NTPsFit > 4) {
          newNTPSFit -= 2;
        } else {
          newNTPSFit -= 1;
        }
        if(lastTP.NTPsFit < 3) newNTPSFit = 2;
        if(newNTPSFit < minPtsFit) newNTPSFit = minPtsFit;
        lastTP.NTPsFit = newNTPSFit;
        if(prt) mf::LogVerbatim("TC")<<"  Bad FitChi "<<lastTP.FitChi<<" Reduced NTPsFit to "<<lastTP.NTPsFit<<" Pass "<<tj.Pass;
        FitTraj(tjs, tj);
        if(lastTP.FitChi > prevChi) {
          if(prt) mf::LogVerbatim("TC")<<"  Chisq is increasing "<<lastTP.FitChi<<"  Try to remove an earlier bad hit";
          MaskBadTPs(tj, 1.5);
          ++ntry;
          if(ntry == 2) break;
        }
        prevChi = lastTP.FitChi;
        if(lastTP.NTPsFit == minPtsFit) break;
      } // lastTP.FitChi > 2 && lastTP.NTPsFit > 2
    }
    // last ditch attempt. Drop the last hit
    if(tj.Pts.size() > fMinPtsFit[tj.Pass] && lastTP.FitChi > maxChi) {
      if(prt) mf::LogVerbatim("TC")<<"  Last try. Drop last TP "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;
      UnsetUsedHits(tjs, lastTP);
      DefineHitPos(lastTP);
      SetEndPoints(tjs, tj);
      lastPt = tj.EndPt[1];
      FitTraj(tjs, tj);
      fUpdateTrajOK = true;
      fMaskedLastTP = true;
    }
    if(prt) mf::LogVerbatim("TC")<<"  Fit done. Chi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;

    if(tj.EndPt[0] == tj.EndPt[1]) {
      fUpdateTrajOK = false;
      return;
    }
    
    // Don't let the angle error get too small too soon. Stepping would stop if the first
    // few hits on a low momentum wandering track happen to have a very good fit to a straight line.
    // We will do this by averaging the default starting value of AngErr of the first TP with the current
    // value from FitTraj.
    if(lastPt < 14) {
      float defFrac = 1 / (float)(tj.EndPt[1]);
      lastTP.AngErr = defFrac * tj.Pts[0].AngErr + (1 - defFrac) * lastTP.AngErr;
    }

    UpdateDeltaRMS(tj);
    SetAngleCode(tjs, lastTP);

    fUpdateTrajOK = true;
    return;

  } // UpdateTraj

  //////////////////////////////////////////
  void TrajClusterAlg::UpdateDeltaRMS(Trajectory& tj)
  {
    // Estimate the Delta RMS of the TPs on the end of tj.
    
    unsigned int lastPt = tj.EndPt[1];
    TrajPoint& lastTP = tj.Pts[lastPt];
    
    if(lastTP.Chg == 0) return;
    if(lastPt < 6) return;

    unsigned short ii, ipt, cnt = 0;
    float sum = 0;
    for(ii = 1; ii < tj.Pts.size(); ++ii) {
      ipt = lastPt - ii;
      if(ipt > tj.Pts.size() - 1) break;
      if(tj.Pts[ipt].Chg == 0) continue;
      sum += PointTrajDOCA(tjs, tj.Pts[ipt].Pos[0], tj.Pts[ipt].Pos[1], lastTP);
      ++cnt;
      if(cnt == lastTP.NTPsFit) break;
      if(ipt == 0) break;
    }
    if(cnt < 3) return;
    // RMS of Gaussian distribution is ~1.2 x the average
    // of a one-sided Gaussian distribution (since Delta is > 0)
    lastTP.DeltaRMS = 1.2 * sum / (float)cnt;
    if(lastTP.DeltaRMS < 0.02) lastTP.DeltaRMS = 0.02;

  } // UpdateDeltaRMS
  
  //////////////////////////////////////////
  void TrajClusterAlg::UpdateAveChg(Trajectory& tj)
  {
    
    if(tj.EndPt[1] == 0) return;
    unsigned short lastPt = tj.EndPt[1];
    tj.AveChg = 0;
    tj.Pts[lastPt].AveChg = 0;

    // calculate ave charge and charge RMS using hits in the trajectory
    unsigned short ii, ipt, cnt = 0;
    float fcnt, sum = 0;
    float sum2 = 0;
    // Don't include the first point in the average. It will be too
    // low if this is a stopping/starting particle
    for(ii = 0; ii < tj.Pts.size(); ++ii) {
      ipt = tj.EndPt[1] - ii;
      if(ipt == 0) break;
      if(tj.Pts[ipt].Chg == 0) continue;
      ++cnt;
      sum += tj.Pts[ipt].Chg;
      sum2 += tj.Pts[ipt].Chg * tj.Pts[ipt].Chg;
      if(cnt == fNPtsAve) break;
    } // iii
    if(cnt == 0) return;
    fcnt = cnt;
    sum /= fcnt;
    tj.AveChg = sum;
    tj.Pts[lastPt].AveChg = sum;
    // define the first point average charge if necessary
    if(tj.Pts[tj.EndPt[0]].AveChg <= 0) tj.Pts[tj.EndPt[0]].AveChg = sum;
    if(cnt > 3) {
      float arg = sum2 - fcnt * sum * sum;
      if(arg < 0) arg = 0;
      float rms = sqrt(arg / (fcnt - 1));
      // convert this to a normalized RMS
      rms /= sum;
      // don't let the calculated charge RMS dominate the default
      // RMS until it is well known. Start with 50% error on the
      // charge RMS
      float defFrac = 1 / (float)(tj.EndPt[1]);
      tj.ChgRMS = defFrac * 0.5 + (1 - defFrac) * rms;
      if(tj.EndPt[1] > 10) {
        // don't let it get crazy small
        if(tj.ChgRMS < fChargeCuts[1]) tj.ChgRMS = fChargeCuts[1];
        // or crazy large
        if(tj.ChgRMS > fChargeCuts[2]) tj.ChgRMS = fChargeCuts[2];
      }
      tj.Pts[lastPt].ChgPull = (tj.Pts[lastPt].Chg / tj.AveChg - 1) / tj.ChgRMS;
    } // cnt > 3

  } // UpdateAveChg

  ////////////////////////////////////////////////
  bool TrajClusterAlg::StartTraj(Trajectory& tj, const unsigned int& fromHit, const unsigned int& toHit, const unsigned short& pass)
  {
    // Start a trajectory located at fromHit with direction pointing to toHit
    float fromWire = tjs.fHits[fromHit].WireID.Wire;
    float fromTick = tjs.fHits[fromHit].PeakTime;
    float toWire = tjs.fHits[toHit].WireID.Wire;
    float toTick = tjs.fHits[toHit].PeakTime;
    CTP_t tCTP = EncodeCTP(tjs.fHits[fromHit].WireID);
    return StartTraj(tj, fromWire, fromTick, toWire, toTick, tCTP, pass);
  } // StartTraj

  ////////////////////////////////////////////////
  bool TrajClusterAlg::StartTraj(Trajectory& tj, const float& fromWire, const float& fromTick, const float& toWire, const float& toTick, const CTP_t& tCTP, const unsigned short& pass)
  {
    // Start a simple (seed) trajectory going from (fromWire, toTick) to (toWire, toTick).
    
    // decrement the work ID so we can use it for debugging problems
    --fWorkID;
    if(fWorkID == INT_MIN) fWorkID = -1;
    tj.ID = fWorkID;
    tj.Pass = pass;
    // Assume we are stepping in the positive WSE units direction
    short stepdir = 1;
    int fWire = std::nearbyint(fromWire);
    int tWire = std::nearbyint(toWire);
    if(tWire < fWire) {
      stepdir = -1;
    } else if(tWire == fWire) {
      // on the same wire
      if(toTick < fromTick) stepdir = -1;
    }
    tj.StepDir = stepdir;
    tj.CTP = tCTP;
    
    // create a trajectory point
    TrajPoint tp;
    if(!MakeBareTrajPoint(tjs, fromWire, fromTick, toWire, toTick, tCTP, tp)) {
      mf::LogVerbatim("TC")<<"StartTraj: Failure from MakeBareTrajPoint fromWire "<<fromWire<<" fromTick "<<fromTick<<" toWire "<<toWire<<" toTick "<<toTick;
      return false;
    }
    SetAngleCode(tjs, tp);
    tp.AngErr = 0.1;
    if(tj.ID == debug.WorkID) { prt = true; didPrt = true; debug.Plane = fPlane; TJPrt = tj.ID; debug.WorkID = tj.ID; }
    if(prt) mf::LogVerbatim("TC")<<"StartTraj "<<(int)fromWire<<":"<<(int)fromTick<<" -> "<<(int)toWire<<":"<<(int)toTick<<" StepDir "<<tj.StepDir<<" dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" ang "<<tp.Ang<<" AngleCode "<<tp.AngleCode<<" angErr "<<tp.AngErr<<" ExpectedHitsRMS "<<ExpectedHitsRMS(tp);
    tj.Pts.push_back(tp);
    return true;
    
  } // StartTraj
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ChkInTraj(std::string someText)
  {
    // Check tjs.allTraj -> InTraj associations
    
    if(!tjs.UseAlg[kChkInTraj]) return;
    
    ++fAlgModCount[kChkInTraj];
    
    int tID;
    unsigned int iht;
    unsigned short itj = 0;
    std::vector<unsigned int> tHits;
    std::vector<unsigned int> atHits;
    for(auto& tj : tjs.allTraj) {
      // ignore abandoned trajectories
      if(tj.AlgMod[kKilled]) continue;
      tID = tj.ID;
      for(auto& tp : tj.Pts) {
        if(tp.Hits.size() > 16) {
          tj.AlgMod[kKilled] = true;
          mf::LogWarning("TC")<<"ChkInTraj: More than 16 hits created a UseHit bitset overflow\n";
          fQuitAlg = true;
          return;
        }
      } // tp
      if(tj.AlgMod[kKilled]) {
        std::cout<<someText<<" ChkInTraj hit size mis-match in tj ID "<<tj.ID<<" AlgBitNames";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) std::cout<<" "<<AlgBitNames[ib];
        std::cout<<"\n";
        continue;
      }
      tHits = PutTrajHitsInVector(tj, kUsedHits);
      if(tHits.size() < 2) {
        mf::LogVerbatim("TC")<<someText<<" ChkInTraj: Insufficient hits in traj "<<tj.ID<<" Killing it";
        PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
        tj.AlgMod[kKilled] = true;
        continue;
      }
      std::sort(tHits.begin(), tHits.end());
      atHits.clear();
      for(iht = 0; iht < tjs.fHits.size(); ++iht) {
        if(tjs.fHits[iht].InTraj == tID) atHits.push_back(iht);
      } // iht
      if(atHits.size() < 2) {
        mf::LogVerbatim("TC")<<someText<<" ChkInTraj: Insufficient hits in atHits in traj "<<tj.ID<<" Killing it";
        tj.AlgMod[kKilled] = true;
        continue;
      }
      if(!std::equal(tHits.begin(), tHits.end(), atHits.begin())) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ChkInTraj failed: inTraj - UseHit mis-match for tj ID "<<tID<<" tj.WorkID "<<tj.WorkID<<" atHits size "<<atHits.size()<<" tHits size "<<tHits.size()<<" in CTP "<<tj.CTP<<"\n";
        myprt<<"AlgMods: ";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
        myprt<<"index     inTraj     UseHit \n";
        for(iht = 0; iht < atHits.size(); ++iht) {
          myprt<<"iht "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]]);
          if(iht < tHits.size()) myprt<<" "<<PrintHit(tjs.fHits[tHits[iht]]);
          if(atHits[iht] != tHits[iht]) myprt<<" <<< "<<atHits[iht]<<" != "<<tHits[iht];
          myprt<<"\n";
          fQuitAlg = true;
        } // iht
        if(tHits.size() > atHits.size()) {
          for(iht = atHits.size(); iht < atHits.size(); ++iht) {
            myprt<<"atHits "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]])<<"\n";
          } // iht
          PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
        } // tHit.size > atHits.size()
      }
      // check the VtxID
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] > tjs.vtx.size()) {
          mf::LogVerbatim("TC")<<someText<<" ChkInTraj: Bad VtxID "<<tj.ID;
          std::cout<<someText<<" ChkInTraj: Bad VtxID "<<tj.ID<<" vtx size "<<tjs.vtx.size()<<"\n";
          tj.AlgMod[kKilled] = true;
          PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
          fQuitAlg = true;
          return;
        }
      } // end
      ++itj;
      if(fQuitAlg) return;
    } // tj
    
  } // ChkInTraj

 ////////////////////////////////////////////////
  void TrajClusterAlg::MakeAllTrajClusters()
  {
    // Make clusters from all trajectories in tjs.allTraj
    
    // Merge hits in trajectory points?
    if(fMakeNewHits) MergeTPHits();
    
    ClusterStore cls;
    tjs.tcl.clear();
    tjs.inClus.resize(tjs.fHits.size());
    unsigned int iht;
    for(iht = 0; iht < tjs.inClus.size(); ++iht) tjs.inClus[iht] = 0;
    
    if(prt) mf::LogVerbatim("TC")<<"MakeAllTrajClusters: tjs.allTraj size "<<tjs.allTraj.size();
    
    if(tjs.UseAlg[kChkInTraj]) {
      fQuitAlg = !InTrajOK(tjs, "MATC");
      if(fQuitAlg) {
        mf::LogVerbatim("TC")<<"InTrajOK failed in MakeAllTrajClusters";
        return;
      }
    }
    
    unsigned short itj, endPt0, endPt1;
    
    // Make one cluster for each trajectory. The indexing of trajectory parents
    // should map directly to cluster parents
    short clID = 0;
    for(itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      // big temp
      if(fMode == 3 && !tj.AlgMod[kShowerTj]) continue;
      // ensure that the endPts are correct
      SetEndPoints(tjs, tj);
      auto tHits = PutTrajHitsInVector(tj, kUsedHits);
      if(tHits.empty()) {
        mf::LogWarning("TC")<<"MakeAllTrajClusters: No hits found in trajectory "<<itj<<" so skip it";
        PrintTrajectory("MATC", tjs, tj, USHRT_MAX);
        continue;
      } // error
      // special handling for shower Tjs: Sort the hits by distance from the start position
      if(tj.AlgMod[kShowerTj]) {
        std::vector<SortEntry> sortVec(tHits.size());
        SortEntry sortEntry;
        std::array<float, 2> hpos;
        for(unsigned short ii = 0; ii < tHits.size(); ++ii) {
          unsigned int iht = tHits[ii];
          hpos[0] = tjs.fHits[iht].WireID.Wire;
          hpos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
          sortEntry.index = ii;
          sortEntry.length = PosSep2(hpos, tj.Pts[0].Pos);
          sortVec[ii] = sortEntry;
        } // ii
        std::sort(sortVec.begin(), sortVec.end(), lessThan);
        // make a temp vector
        std::vector<unsigned int> tmp(sortVec.size());
        // enter the sorted hits
        for(unsigned short ii = 0; ii < sortVec.size(); ++ii) tmp[ii] = tHits[sortVec[ii].index];
        tHits = tmp;
      } // showerTj
      // count AlgMod bits
      for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) ++fAlgModCount[ib];
      ++clID;
      cls.ID = clID;
      cls.CTP = tj.CTP;
      cls.PDGCode = tj.PDGCode;
      cls.ParentCluster = tj.ParentTrajID - 1;
      endPt0 = tj.EndPt[0];
      cls.BeginWir = tj.Pts[endPt0].Pos[0];
      cls.BeginTim = tj.Pts[endPt0].Pos[1] / tjs.UnitsPerTick;
      cls.BeginAng = tj.Pts[endPt0].Ang;
      cls.BeginChg = tj.Pts[endPt0].Chg;
      cls.BeginVtx = tj.VtxID[0]-1;
      endPt1 = tj.EndPt[1];
      cls.EndWir = tj.Pts[endPt1].Pos[0];
      cls.EndTim = tj.Pts[endPt1].Pos[1] / tjs.UnitsPerTick;
      cls.EndAng = tj.Pts[endPt1].Ang;
      cls.EndChg = tj.Pts[endPt1].Chg;
      cls.EndVtx = tj.VtxID[1]-1;
      cls.tclhits = tHits;
      // Set the traj info
      tj.ClusterIndex = tjs.tcl.size();
      tjs.tcl.push_back(cls);
      // do some checking and define tjs.inClus
      geo::PlaneID planeID = DecodeCTP(cls.CTP);
      for(unsigned short ii = 0; ii < cls.tclhits.size(); ++ii) {
        unsigned int iht = cls.tclhits[ii];
        if(tjs.fHits[iht].WireID.Plane != planeID.Plane ||
           tjs.fHits[iht].WireID.Cryostat != planeID.Cryostat ||
           tjs.fHits[iht].WireID.TPC != planeID.TPC) {
          mf::LogWarning("TC")<<"MakeAllTrajClusters: Bad OLD hit CTP in itj "<<itj<<" hit "<<PrintHit(tjs.fHits[iht])<<" WorkID "<<tjs.allTraj[itj].WorkID<<" Plane "<<tjs.fHits[iht].WireID.Plane<<" vs "<<planeID.Plane<<" Cstat "<<tjs.fHits[iht].WireID.Cryostat<<" vs "<<planeID.Cryostat<<" TPC "<<tjs.fHits[iht].WireID.TPC<<" vs "<<planeID.TPC;
          fQuitAlg = true;
          return;
        }
        if(tjs.inClus[iht] != 0) {
          mf::LogWarning("TC")<<"MakeAllTrajClusters: Trying to assign tj.ID "<<tj.ID<<" hit "<<iht<<" "<<PrintHit(tjs.fHits[iht])<<" to already-assigned cluster "<<tjs.inClus[iht]<<" workID "<<tj.WorkID;
          fQuitAlg = true;
          return;
        }
        tjs.inClus[iht] = clID;
      } //iht
    } // itj
    
    // Ensure that all PFParticles have associated clusters
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      MatchStruct& ms = tjs.matchVec[tjs.matchVecPFPList[ipfp]];
      for(auto& tjid : ms.TjIDs) {
        Trajectory& tj = tjs.allTraj[tjid - 1];
        if(tj.ID != tjid || tj.AlgMod[kKilled] || !tj.AlgMod[kMatch3D] || tj.ClusterIndex > tjs.tcl.size() - 1) {
          std::cout<<"MATC bad PFP -> cluster association "<<ipfp<<" TjID  "<<tjid<<" ClusterIndex "<<tj.ClusterIndex<<"\n";
          PrintAllTraj("Bad", tjs, debug, tjid - 1, USHRT_MAX, false);
        }
      } // tjid
    } // ipfp

    // Define the Hits vector for showers
    for(unsigned short ish = 0; ish < tjs.showers.size(); ++ish) {
      for(auto& TjID : tjs.showers[ish].TjIDs) {
        Trajectory& tj = tjs.allTraj[TjID - 1];
        auto tjHits = PutTrajHitsInVector(tj, kUsedHits);
        tjs.showers[ish].Hits.insert(tjs.showers[ish].Hits.begin(), tjHits.begin(), tjHits.end());
      } // TjID
    } // ish
    
  } // MakeAllTrajClusters
  
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindVtxTjs()
  {
    // Look for vertex trajectories in all vertices in the current fCTP
    if(!tjs.UseAlg[kVtxTj]) return;
    
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      if(tjs.vtx[ivx].Stat[kVtxKilled] == 0) continue;
      if(tjs.vtx[ivx].CTP != fCTP) continue;
     if(tjs.vtx[ivx].Stat[kVtxTrjTried]) continue;
      FindVtxTraj(ivx);
    } // ivx
    
  } // FindVtxTjs
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindVtxTraj(unsigned short ivx)
  {
    // Look for available hits in the vicinity of this vertex and try to make
    // a vertex trajectory from them
    
    if(!tjs.UseAlg[kVtxTj]) return;
    
    if(ivx > tjs.vtx.size() - 1) return;
    if(tjs.vtx[ivx].Stat[kVtxTrjTried]) return;
    VtxStore& theVtx = tjs.vtx[ivx];
    
    std::array<int, 2> wireWindow;
    std::array<float, 2> timeWindow;
    
    // on the first try we look for small angle trajectories which will have hits
    // with a large wire window and a small time window
    // Vertex2DCuts fcl input usage
    // 0 User definition of a short Tj => max number of Tj points
    // 1 max separation between a vertex and the start of a trajectory for a user-defined short Tj
    // 2 max separation for a user-defined long Tj
    // 3 max position pull when attaching a Tj to a vertex
    // 4 max position error for creating a Tj or attaching Tjs to an existing vertex
    // 5 Min MCSMom of Tjs that can be used to create a vertex
    // 6 min frac of Points/Wire between a vtx and a Tj. Ideally one if the efficiency is good
    // 7 min Score
    // 8 ID of a vertex for printing special debugging information
    wireWindow[0] = std::nearbyint(theVtx.Pos[0] - tjs.Vertex2DCuts[2]);
    wireWindow[1] = std::nearbyint(theVtx.Pos[0] + tjs.Vertex2DCuts[2]);
    timeWindow[0] = theVtx.Pos[1] - 5;
    timeWindow[1] = theVtx.Pos[1] + 5;
    
    geo::PlaneID planeID = DecodeCTP(theVtx.CTP);
    unsigned short ipl = planeID.Plane;
    
    if(prt) mf::LogVerbatim("TC")<<"inside FindVtxTraj "<<theVtx.ID<<" Window "<<wireWindow[0]<<" "<<wireWindow[1]<<" "<<timeWindow[0]<<" "<<timeWindow[1]<<" in plane "<<ipl;
    
    // find nearby available hits
    bool hitsNear;
    std::vector<unsigned int> closeHits = FindCloseHits(tjs, wireWindow, timeWindow, ipl, kUnusedHits, true, hitsNear);
    if(closeHits.empty()) return;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"closeHits";
      for(auto& iht : closeHits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
    }
    // sort by distance from the vertex
    std::vector<SortEntry> sortVec(closeHits.size());
    SortEntry sortEntry;
    for(unsigned short ii = 0; ii < closeHits.size(); ++ii) {
      unsigned int iht = closeHits[ii];
      float dw = tjs.fHits[iht].WireID.Wire - theVtx.Pos[0];
      float dt = tjs.UnitsPerTick * tjs.fHits[iht].PeakTime - theVtx.Pos[1];
      float d2 = dw * dw + dt * dt;
      sortEntry.index = ii;
      sortEntry.length = d2;
      sortVec[ii] = sortEntry;
    } // ii
    std::sort(sortVec.begin(), sortVec.end(), lessThan);
    unsigned int vWire = std::nearbyint(theVtx.Pos[0]);
    int vTick = theVtx.Pos[1]/tjs.UnitsPerTick;
    if(prt) PrintHeader("FVT");
    for(unsigned short ii = 0; ii < closeHits.size(); ++ii) {
      unsigned int iht = closeHits[sortVec[ii].index];
      if(tjs.fHits[iht].InTraj > 0) continue;
      // the direction will be poorly defined if a hit is very close to the vertex and it is in this list.
      // Ignore these hits
      if(tjs.fHits[iht].WireID.Wire == vWire) {
        // on the vertex wire. Check for a close time
        if(abs(tjs.fHits[iht].PeakTime - vTick) < 10) continue;
      } // hit on vtx wire
      float toWire = tjs.fHits[iht].WireID.Wire;
      float toTick = tjs.fHits[iht].PeakTime;
      // assume the last pass and fix it later after the angle is calculated
      unsigned short pass = fMinPts.size() - 1;
      Trajectory tj;
      if(!StartTraj(tj, theVtx.Pos[0], theVtx.Pos[1]/tjs.UnitsPerTick, toWire, toTick, theVtx.CTP, pass)) continue;
      // ensure that the first TP is good
      if(tj.Pts[0].Pos[0] < 0) continue;
      //      std::cout<<"fvt "<<theVtx.ID<<" "<<tj.ID<<" vtx0 "<<theVtx.Pos[0]<<" hit "<<PrintHit(tjs.fHits[iht])<<" StepDir "<<tj.StepDir<<"\n";
      tj.VtxID[0] = theVtx.ID;
      TrajPoint& tp = tj.Pts[0];
      // Move the Pt to the hit
      MoveTPToWire(tp, toWire);
      // attach the hit
      tp.Hits.push_back(iht);
      tp.UseHit[tp.Hits.size()-1] = true;
      tjs.fHits[iht].InTraj = tj.ID;
      tp.UseHit[tp.Hits.size()-1] = false;
      if(prt) PrintTrajPoint("FVT", tjs, 0, tj.StepDir, tj.Pass, tp);
      // Step away and see what happens
      prt = prt;
      StepCrawl(tj);
      // check for a major failure
      if(fQuitAlg) return;
      // Check the quality of the trajectory
      CheckTraj(tj);
      if(!fGoodTraj || NumPtsWithCharge(tjs, tj, true) < fMinPts[tj.Pass]) {
        if(prt) mf::LogVerbatim("TC")<<" xxxxxxx Not enough points "<<NumPtsWithCharge(tjs, tj, true)<<" minimum "<<fMinPts[tj.Pass]<<" or !fGoodTraj";
        ReleaseHits(tjs, tj);
        continue;
      }
      if(prt) prt = false;
      tj.AlgMod[kVtxTj] = true;
      fQuitAlg = !StoreTraj(tjs, tj);
      if(tjs.UseAlg[kChkInTraj]) {
        fQuitAlg = !InTrajOK(tjs, "FVT");
        if(fQuitAlg) {
          mf::LogVerbatim("TC")<<"InTrajOK failed in FindVtxTraj";
          return;
        }
      }
      if(prt) mf::LogVerbatim("TC")<<"FindVtxTraj: calling StoreTraj with npts "<<tj.EndPt[1];
    } // ii
    
    // Flag this as tried so we don't try again
    tjs.vtx[ivx].Stat[kVtxTrjTried] = true;
  } // FindVtxTraj
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet)
  {
    unsigned short localIndex;
    GetHitMultiplet(theHit, hitsInMultiplet, localIndex);
  } // GetHitMultiplet
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet, unsigned short& localIndex)
  {
    hitsInMultiplet.clear();
    if(theHit > tjs.fHits.size() - 1) return;
    hitsInMultiplet.resize(1);
    hitsInMultiplet[0] = theHit;
    
    float hitSep;
    unsigned int theWire = tjs.fHits[theHit].WireID.Wire;
    unsigned short ipl = tjs.fHits[theHit].WireID.Plane;
    float theTime = tjs.fHits[theHit].PeakTime;
    float theRMS = tjs.fHits[theHit].RMS;
    float narrowHitCut = 1.5 * fAveHitRMS[ipl];
    bool theHitIsNarrow = (theRMS < narrowHitCut);
    float maxPeak = tjs.fHits[theHit].PeakAmplitude;
    unsigned short imTall = theHit;
    unsigned short nNarrow = 0;
    if(theHitIsNarrow) nNarrow = 1;
//    if(prt) mf::LogVerbatim("TC")<<"GetHitMultiplet theHit "<<theHit<<" "<<PrintHit(tjs.fHits[theHit])<<" RMS "<<tjs.fHits[theHit].RMS<<" aveRMS "<<fAveHitRMS[ipl]<<" Amp "<<(int)tjs.fHits[theHit].PeakAmplitude;
    // look for hits < theTime but within hitSep
    if(theHit > 0) {
      for(unsigned int iht = theHit - 1; iht != 0; --iht) {
        if(tjs.fHits[iht].WireID.Wire != theWire) break;
        if(tjs.fHits[iht].WireID.Plane != ipl) break;
        if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
        if(tjs.fHits[iht].RMS > theRMS) {
          hitSep = fMultHitSep * tjs.fHits[iht].RMS;
          theRMS = tjs.fHits[iht].RMS;
        } else {
          hitSep = fMultHitSep * theRMS;
        }
        if(theTime - tjs.fHits[iht].PeakTime > hitSep) break;
//        if(prt) mf::LogVerbatim("TC")<<" iht- "<<iht<<" "<<tjs.fHits[iht].WireID.Plane<<":"<<PrintHit(tjs.fHits[iht])<<" RMS "<<tjs.fHits[iht].RMS<<" dt "<<theTime - tjs.fHits[iht].PeakTime<<" "<<hitSep<<" Amp "<<(int)tjs.fHits[iht].PeakAmplitude;
         hitsInMultiplet.push_back(iht);
        if(tjs.fHits[iht].RMS < narrowHitCut) ++nNarrow;
        if(tjs.fHits[iht].PeakAmplitude > maxPeak) {
          maxPeak = tjs.fHits[iht].PeakAmplitude;
          imTall = iht;
        }
        theTime = tjs.fHits[iht].PeakTime;
        if(iht == 0) break;
      } // iht
    } // iht > 0
    localIndex = hitsInMultiplet.size() - 1;
    // reverse the order so that hitsInMuliplet will be
    // returned in increasing time order
    if(hitsInMultiplet.size() > 1) std::reverse(hitsInMultiplet.begin(), hitsInMultiplet.end());
    // look for hits > theTime but within hitSep
    theTime = tjs.fHits[theHit].PeakTime;
    theRMS = tjs.fHits[theHit].RMS;
    for(unsigned int iht = theHit + 1; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].WireID.Wire != theWire) break;
      if(tjs.fHits[iht].WireID.Plane != ipl) break;
      if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
      if(tjs.fHits[iht].RMS > theRMS) {
        hitSep = fMultHitSep * tjs.fHits[iht].RMS;
        theRMS = tjs.fHits[iht].RMS;
      } else {
        hitSep = fMultHitSep * theRMS;
      }
      if(tjs.fHits[iht].PeakTime - theTime > hitSep) break;
//      if(prt) mf::LogVerbatim("TC")<<" iht+ "<<iht<<" "<<PrintHit(tjs.fHits[iht])<<" dt "<<(theTime - tjs.fHits[iht].PeakTime)<<" RMS "<<tjs.fHits[iht].RMS<<" "<<hitSep<<" Amp "<<(int)tjs.fHits[iht].PeakAmplitude;
       hitsInMultiplet.push_back(iht);
      if(tjs.fHits[iht].RMS < narrowHitCut) ++nNarrow;
      if(tjs.fHits[iht].PeakAmplitude > maxPeak) {
        maxPeak = tjs.fHits[iht].PeakAmplitude;
        imTall = iht;
      }
      theTime = tjs.fHits[iht].PeakTime;
    } // iht

    if(hitsInMultiplet.size() == 1) return;
    
    if(hitsInMultiplet.size() > 16) {
      // Found > 16 hits in a multiplet which would be bad for UseHit. Truncate it
      hitsInMultiplet.resize(16);
      return;
    }
    
    // Don't make a multiplet that includes a tall narrow hit with short fat hits
    if(nNarrow == hitsInMultiplet.size()) return;
    if(nNarrow == 0) return;
    
    if(theHitIsNarrow && theHit == imTall) {
//      if(prt) mf::LogVerbatim("TC")<<" theHit is narrow and tall. Use only it";
      // theHit is narrow and it is the highest amplitude hit in the multiplet. Ignore any
      // others that are short and fat
      auto tmp = hitsInMultiplet;
      tmp.resize(1);
      tmp[0] = theHit;
      hitsInMultiplet = tmp;
    } else {
      // theHit is not narrow and it is not the tallest. Ignore a single hit if it is
      // the tallest and narrow
//      if(prt) mf::LogVerbatim("TC")<<" theHit  is not narrow or tall";
      if(tjs.fHits[imTall].RMS < narrowHitCut) {
        unsigned short killMe = 0;
        for(unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
          if(hitsInMultiplet[ii] == imTall) {
            killMe = ii;
            break;
          }
        } // ii
        hitsInMultiplet.erase(hitsInMultiplet.begin() + killMe);
      } // tjs.fHits[imTall].RMS < narrowHitCut
    } // narrow / tall test

  } // GetHitMultiplet
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FillWireHitRange(const geo::TPCID& tpcid)
  {
    // fills the WireHitRange vector. Slightly modified version of the one in ClusterCrawlerAlg.
    
    // determine the number of planes
    geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    unsigned short nplanes = TPC.Nplanes();
    tjs.NumPlanes = nplanes;
    tjs.TPCID = tpcid;
    
    // Y,Z limits of the detector
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &thetpc = tjs.geom->TPC(tpc, cstat);
    thetpc.LocalToWorld(local,world);
    // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
    tjs.XLo = world[0]-tjs.geom->DetHalfWidth(tpc,cstat) + 5;
    tjs.XHi = world[0]+tjs.geom->DetHalfWidth(tpc,cstat) - 5;
    tjs.YLo = world[1]-tjs.geom->DetHalfHeight(tpc,cstat) + 5;
    tjs.YHi = world[1]+tjs.geom->DetHalfHeight(tpc,cstat) - 5;
    tjs.ZLo = world[2]-tjs.geom->DetLength(tpc,cstat)/2 + 5;
    tjs.ZHi = world[2]+tjs.geom->DetLength(tpc,cstat)/2 - 5;
    
    lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    if(!tjs.WireHitRange.empty()) tjs.WireHitRange.clear();
    
    // initialize everything
    tjs.WireHitRange.resize(nplanes);
    tjs.WireHitRangeCstat = cstat;
    tjs.WireHitRangeTPC = tpc;
    tjs.FirstWire.resize(nplanes);
    tjs.LastWire.resize(nplanes);
    tjs.NumWires.resize(nplanes);
    tjs.MaxPos0.resize(nplanes);
    tjs.MaxPos1.resize(nplanes);
    fAveHitRMS.resize(nplanes, nplanes);
    
    std::pair<int, int> flag;
    flag.first = -2; flag.second = -2;
    
    // Calculate tjs.UnitsPerTick, the scale factor to convert a tick into
    // Wire Spacing Equivalent (WSE) units where the wire spacing in this plane = 1.
    // Strictly speaking this factor should be calculated for each plane to handle the
    // case where the wire spacing is different in each plane. Deal with this later if
    // the approximation used here fails.
    
    raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel(0, 0, (int)tpc, (int)cstat);
    float wirePitch = tjs.geom->WirePitch(tjs.geom->View(channel));
    float tickToDist = tjs.detprop->DriftVelocity(tjs.detprop->Efield(),tjs.detprop->Temperature());
    tickToDist *= 1.e-3 * tjs.detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    tjs.UnitsPerTick = tickToDist / wirePitch;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      tjs.FirstWire[ipl] = INT_MAX;
      tjs.LastWire[ipl] = 0;
      tjs.NumWires[ipl] = tjs.geom->Nwires(ipl, tpc, cstat);
      tjs.WireHitRange[ipl].resize(tjs.NumWires[ipl], flag);
      tjs.MaxPos0[ipl] = (float)tjs.NumWires[ipl] - 0.5;
      tjs.MaxPos1[ipl] = (float)tjs.detprop->NumberTimeSamples() * tjs.UnitsPerTick;
    }
    
    // overwrite with the "dead wires" condition
    flag.first = -1; flag.second = -1;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        raw::ChannelID_t chan = tjs.geom->PlaneWireToChannel((int)ipl, (int)wire, (int)tpc, (int)cstat);
        if(!channelStatus.IsGood(chan)) tjs.WireHitRange[ipl][wire] = flag;
      } // wire
    } // ipl

    unsigned int lastwire = 0, lastipl = 0;
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].WireID.Cryostat != cstat) continue;
      if(tjs.fHits[iht].WireID.TPC != tpc) continue;
      unsigned short ipl = tjs.fHits[iht].WireID.Plane;
      unsigned int wire = tjs.fHits[iht].WireID.Wire;
      if(wire > tjs.NumWires[ipl] - 1) {
        mf::LogWarning("TC")<<"FillWireHitRange: Invalid wire number "<<wire<<" > "<<tjs.NumWires[ipl] - 1<<" in plane "<<ipl<<" Quitting";
        fQuitAlg = true;
        return;
      } // too large wire number
      if(ipl == lastipl && wire < lastwire) {
        mf::LogWarning("TC")<<"FillWireHitRange: Hits are not in increasing wire order. Quitting ";
        fQuitAlg = true;
        return;
      } // hits out of order
      lastwire = wire;
      lastipl = ipl;
      if(tjs.FirstWire[ipl] == INT_MAX) tjs.FirstWire[ipl] = wire;
      if(tjs.WireHitRange[ipl][wire].first < 0) tjs.WireHitRange[ipl][wire].first = iht;
      tjs.WireHitRange[ipl][wire].second = iht + 1;
      tjs.LastWire[ipl] = wire + 1;
    } // iht
    
    if(!CheckWireHitRange()) {
      fQuitAlg = true;
      return;
    }

    // Find the average multiplicity 1 hit RMS and calculate the expected max RMS for each range
    if(fDebugMode && (int)tpc == debug.TPC) std::cout<<"tjs.UnitsPerTick "<<std::setprecision(3)<<tjs.UnitsPerTick<<"\n";
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      float sumRMS = 0;
      float sumAmp = 0;
      unsigned int cnt = 0;
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        // don't let noisy wires screw up the calculation
        if(lastHit - firstHit > 100) continue;
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].Multiplicity != 1) continue;
          if(tjs.fHits[iht].GoodnessOfFit < 0 || tjs.fHits[iht].GoodnessOfFit > 100) continue;
          // don't let a lot of runt hits screw up the calculation
          if(tjs.fHits[iht].PeakAmplitude < 1) continue;
          ++cnt;
          sumRMS += tjs.fHits[iht].RMS;
          sumAmp += tjs.fHits[iht].PeakAmplitude;
        } // iht
      } // wire
      if(cnt < 4) continue;
      fAveHitRMS[ipl] = sumRMS/(float)cnt;
      sumAmp  /= (float)cnt;
      if(fDebugMode) std::cout<<"Pln "<<ipl<<" fAveHitRMS "<<fAveHitRMS[ipl]<<" Ave PeakAmplitude "<<sumAmp<<"\n";
    } // ipl
    
  } // FillWireHitRange
  
  ////////////////////////////////////////////////
  std::vector<recob::Hit> TrajClusterAlg::YieldHits()
  {
    // Create the final recob::hits and return them
    std::vector<recob::Hit> tmp;
    tmp.reserve(tjs.fHits.size());
    for(auto& tcHit : tjs.fHits) {
      geo::PlaneID planeID = geo::PlaneID(tcHit.WireID.Cryostat, tcHit.WireID.TPC, tcHit.WireID.Plane);
      raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel((int)tcHit.WireID.Plane, (int)tcHit.WireID.Wire, (int)tcHit.WireID.TPC, (int)tcHit.WireID.Cryostat);
      tmp.emplace_back(channel,
                       tcHit.StartTick, tcHit.EndTick,
                       tcHit.PeakTime, tcHit.SigmaPeakTime,
                       tcHit.RMS,
                       tcHit.PeakAmplitude, tcHit.SigmaPeakAmp,
                       tcHit.Integral, tcHit.Integral, tcHit.SigmaIntegral,
                       tcHit.Multiplicity, tcHit.LocalIndex,
                       tcHit.GoodnessOfFit, tcHit.NDOF,
                       tjs.geom->View(channel),
                       tjs.geom->SignalType(planeID),
                       tcHit.WireID
                       );
    } // tcHit
     return tmp;
  } // YieldHits
  
  ////////////////////////////////////////////////
  float TrajClusterAlg::ExpectedHitsRMS(TrajPoint const& tp)
  {
    // returns the expected RMS of hits for the trajectory point in ticks
    if(std::abs(tp.Dir[0]) > 0.001) {
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      return 1.5 * fAveHitRMS[planeID.Plane] + 2 * std::abs(tp.Dir[1]/tp.Dir[0])/tjs.UnitsPerTick;
    } else {
      return 500;
    }
  } // ExpectedHitsRMS
  
  ////////////////////////////////////////////////
  bool TrajClusterAlg::CheckWireHitRange()
  {
    // do a QC check
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        // No hits or dead wire
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        if(lastHit > tjs.fHits.size()) {
          mf::LogWarning("TC")<<"CheckWireHitRange: Invalid lastHit "<<lastHit<<" > fHits.size "<<tjs.fHits.size()<<" in plane "<<ipl;
          std::cout<<"CheckWireHitRange: Invalid lastHit "<<lastHit<<" > fHits.size "<<tjs.fHits.size()<<" in plane "<<ipl<<"\n";
          return false;
        }
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].WireID.Plane != ipl) {
            mf::LogWarning("TC")<<"CheckWireHitRange: Invalid plane "<<tjs.fHits[iht].WireID.Plane<<" != "<<ipl;
            std::cout<<"CheckWireHitRange: Invalid plane "<<tjs.fHits[iht].WireID.Plane<<" != "<<ipl<<"\n";
            return false;
          }
          if(tjs.fHits[iht].WireID.Wire != wire) {
            mf::LogWarning("TC")<<"CheckWireHitRange: Invalid wire "<<tjs.fHits[iht].WireID.Wire<<" != "<<wire<<" in plane "<<ipl;
            std::cout<<"CheckWireHitRange: Invalid wire "<<tjs.fHits[iht].WireID.Wire<<" != "<<wire<<" in plane "<<ipl<<"\n";
            return false;
          }
        } // iht
      } // wire
    } // ipl
    
    return true;

  } // CheckWireHitRange
  
  ////////////////////////////////////////////////
  bool TrajClusterAlg::EraseHit(const unsigned int& delHit)
  {
    // Erases delHit and makes corrections to allTraj and WireHitRange
    if(delHit > tjs.fHits.size() - 1) {
      mf::LogWarning("TC")<<"Trying to erase an invalid hit";
      return false;
    }
    // erase the hit
    tjs.fHits.erase(tjs.fHits.begin()+delHit);
    // Correct WireHitRange
    int idelHit = delHit;
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      for(unsigned int wire = tjs.FirstWire[ipl]; wire < tjs.LastWire[ipl];  ++wire) {
        // ignore wires with no hits or dead
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        if(idelHit > 0 && tjs.WireHitRange[ipl][wire].first > idelHit) --tjs.WireHitRange[ipl][wire].first;
        if(tjs.WireHitRange[ipl][wire].second > idelHit) --tjs.WireHitRange[ipl][wire].second;
        // Deal with the situation where this is the only hit on a wire
        int firstHit = tjs.WireHitRange[ipl][wire].first;
        int lastHit = tjs.WireHitRange[ipl][wire].second - 1;
        if(lastHit <= firstHit) {
          // erasing the only hit on this wire
          tjs.WireHitRange[ipl][wire].first = -2;
          tjs.WireHitRange[ipl][wire].second = -2;
          // skip checking
          continue;
        }
        // check the first hit
        if(tjs.fHits[firstHit].WireID.Plane != ipl || tjs.fHits[firstHit].WireID.Wire != wire) {
          std::cout<<"WireHitRange screwup on firstHit "<<tjs.fHits[firstHit].WireID.Plane<<":"<<tjs.fHits[firstHit].WireID.Wire;
          std::cout<<" != "<<ipl<<":"<<wire<<"\n";
          exit(1);
        } // and the last hit
        if(tjs.fHits[lastHit].WireID.Plane != ipl || tjs.fHits[lastHit].WireID.Wire != wire) {
          std::cout<<"WireHitRange screwup on lastHit "<<tjs.fHits[lastHit].WireID.Plane<<":"<<tjs.fHits[lastHit].WireID.Wire;
          std::cout<<" != "<<ipl<<":"<<wire<<"\n";
          exit(1);
        } // error checking
      } // wire
    } // ipl
    
    // do another sanity check
    if(!CheckWireHitRange()) return false;
    
    // now fix the Trajectory point hit -> tjs.fHits.InTraj association. The first step is to
    // remove any use of delHit in all trajectory points. The second is to remove any trajectory point that
    // uses only delHit to define the hit position and is therefore no longer valid
    for(auto& tj : tjs.allTraj) {
      unsigned short killPt = USHRT_MAX;
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        unsigned short killii = USHRT_MAX;
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tp.Hits[ii] == delHit) {
            // delHit is used in this TP so we need to remove it
            killii = ii;
          } else if(tp.Hits[ii] > delHit) {
            // delHit comes later in the hit collection so we need to simply correct it
            --tp.Hits[ii];
          }
        } // ii
        if(killii != USHRT_MAX) {
          // We need to erase the reference to this hit in the TP
          tp.Hits.erase(tp.Hits.begin() + killii);
          // shift UseHit being careful not to go outside the bounds
          unsigned short maxSize = tp.Hits.size();
          if(maxSize == 16) maxSize = 15;
          for(unsigned short ii = killii; ii < maxSize; ++ii) tp.UseHit[ii] = tp.UseHit[ii + 1];
          // Flag this TP for deletion if there are no other hits
          if(tp.Hits.empty()) killPt = ipt;
        } // killii != USHRT_MAX
      } // ipt
      // delHit was the only hit used in this TP and it was erased so delete the TP
      if(killPt != USHRT_MAX) {
        tj.Pts.erase(tj.Pts.begin() + killPt);
        SetEndPoints(tjs, tj);
      }
    } // tj

    return true;
  } // EraseHit
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::DefineHit(TCHit& tcHit, CTP_t& hitCTP, unsigned int& hitWire)
  {
    // Defines the hit WireID, channel, etc using hitCTP and hitWire
    geo::PlaneID planeID = DecodeCTP(hitCTP);
    tcHit.WireID = geo::WireID(planeID, hitWire);
//    hitCTP.Channel = tjs.geom->PlaneWireToChannel((int)planeID.Plane,(int)hitWire,(int)planeID.TPC,(int)planeID.Cryostat);
  } // DefineHit
  
  ////////////////////////////////////////////////
  unsigned int TrajClusterAlg::CreateHit(TCHit tcHit)
  {
    // Creates a hit in tjs.fHits using the supplied information. Returns UINT_MAX if there is failure.
    // Returns the index of the newly created hit.
    unsigned short newHitPlane = tcHit.WireID.Plane;
    unsigned int newHitWire = tcHit.WireID.Wire;
    // don't try to create a hit on a dead wire
    if(tjs.WireHitRange[newHitPlane][newHitWire].first == -1) return UINT_MAX;
    
    // Figure out where to put it
    unsigned int newHitIndex = UINT_MAX;
    if(tjs.WireHitRange[newHitPlane][newHitWire].first == -2) {
      // We want to put this hit on a wire that currently has none. Find the next wire that has a hit.
      // First look in the plane in which we want to put it
      for(unsigned int wire = newHitWire + 1; wire < tjs.NumWires[newHitPlane]; ++wire) {
        if(tjs.WireHitRange[newHitPlane][wire].first >= 0) {
          newHitIndex = tjs.WireHitRange[newHitPlane][wire].first;
          break;
        }
      } // wire
      // if not found in this plane look in the rest of the planes
      if(newHitIndex == UINT_MAX) {
        for(unsigned short ipl = newHitPlane + 1; ipl < tjs.NumPlanes; ++ipl) {
          for(unsigned int wire = tjs.FirstWire[ipl]; wire < tjs.LastWire[ipl]; ++wire) {
            if(tjs.WireHitRange[ipl][wire].first >= 0) {
              newHitIndex = tjs.WireHitRange[ipl][wire].first;
              break;
            }
          } // wire
          if(newHitIndex != UINT_MAX) break;
        } // ipl
      } // newHitIndex == UINT_MAX
    } else {
      // Hits exist on this wire
      unsigned int firstHit = tjs.WireHitRange[newHitPlane][newHitWire].first;
      unsigned int lastHit = tjs.WireHitRange[newHitPlane][newHitWire].second - 1;
      if(tcHit.PeakTime <= tjs.fHits[firstHit].PeakTime) {
        // new hit is earlier in time so it should be inserted before firstHit
        newHitIndex = firstHit;
      } else if(tcHit.PeakTime > tjs.fHits[lastHit].PeakTime) {
        // new hit is later so it should inserted after lastHit
        newHitIndex = lastHit + 1;
      } else {
        // new hit is somewhere in the middle
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tcHit.PeakTime > tjs.fHits[iht].PeakTime && tcHit.PeakTime <= tjs.fHits[iht+1].PeakTime) {
            // found it
            newHitIndex = iht + 1;
            break;
          }
        } // iht
      } // new hit in the middle
    } // Hits exist on this wire
    
    // this shouldn't be possible
    if(newHitIndex == UINT_MAX) {
      std::cout<<"CreateHit: Failed to find newHitIndex for new hit "<<PrintHit(tcHit)<<"\n";
      return newHitIndex;
    }
    
    // insert the hit
    tjs.fHits.insert(tjs.fHits.begin() + newHitIndex, tcHit);
    
    // Correct WireHitRange
    
    // Put the hit on a wire with no existing hits
    if(tjs.WireHitRange[newHitPlane][newHitWire].first == -2) {
      tjs.WireHitRange[newHitPlane][newHitWire].first = newHitIndex;
      tjs.WireHitRange[newHitPlane][newHitWire].second = newHitIndex + 1;
    } else {
      // This wire has hits, one of which is the new hits, so only correct the last hit
      ++tjs.WireHitRange[newHitPlane][newHitWire].second;
    }
    
    // correct the hit ranges in newHitPlane on wires after newHitWire
    for(unsigned int wire = newHitWire + 1; wire <  tjs.LastWire[newHitPlane]; ++wire) {
      // dead wire
      if(tjs.WireHitRange[newHitPlane][wire].first < 0) continue;
      ++tjs.WireHitRange[newHitPlane][wire].first;
      ++tjs.WireHitRange[newHitPlane][wire].second;
      // check the hits
      int firstHit = tjs.WireHitRange[newHitPlane][wire].first;
      int lastHit = tjs.WireHitRange[newHitPlane][wire].second - 1;
      if(tjs.fHits[firstHit].WireID.Plane != newHitPlane || tjs.fHits[firstHit].WireID.Wire != wire) {
        std::cout<<"WireHitRange1 screwup on firstHit "<<tjs.fHits[firstHit].WireID.Plane<<":"<<tjs.fHits[firstHit].WireID.Wire;
        std::cout<<" != "<<newHitPlane<<":"<<wire<<"\n";
        exit(1);
      } // error checking
      if(tjs.fHits[lastHit].WireID.Plane != newHitPlane || tjs.fHits[lastHit].WireID.Wire != wire) {
        std::cout<<"WireHitRange1 screwup on lastHit "<<tjs.fHits[lastHit].WireID.Plane<<":"<<tjs.fHits[lastHit].WireID.Wire;
        std::cout<<" != "<<newHitPlane<<":"<<wire<<"\n";
        exit(1);
      } // error checking
    } // wire
    
    // correct the hit ranges for the later planes
    for(unsigned short ipl = newHitPlane + 1; ipl < tjs.NumPlanes; ++ipl) {
      for(unsigned int wire = tjs.FirstWire[ipl]; wire < tjs.LastWire[ipl]; ++wire) {
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        ++tjs.WireHitRange[ipl][wire].first;
        ++tjs.WireHitRange[ipl][wire].second;
        // check the hits
        int firstHit = tjs.WireHitRange[ipl][wire].first;
        int lastHit = tjs.WireHitRange[ipl][wire].second - 1;
        if(tjs.fHits[firstHit].WireID.Plane != ipl || tjs.fHits[firstHit].WireID.Wire != wire) {
          std::cout<<"WireHitRange2 screwup on firstHit "<<tjs.fHits[firstHit].WireID.Plane<<":"<<tjs.fHits[firstHit].WireID.Wire;
          std::cout<<" != "<<ipl<<":"<<wire<<"\n";
          exit(1);
        } // error checking
        if(tjs.fHits[lastHit].WireID.Plane != ipl || tjs.fHits[lastHit].WireID.Wire != wire) {
          std::cout<<"WireHitRange2 screwup on lastHit "<<tjs.fHits[lastHit].WireID.Plane<<":"<<tjs.fHits[lastHit].WireID.Wire;
          std::cout<<" != "<<ipl<<":"<<wire<<"\n";
          exit(1);
        } // error checking
      } // wire
    } // ipl
    
    
    if(!CheckWireHitRange()) return UINT_MAX;
    
    // now correct the hit indices in the trajectories
    for(auto& tj : tjs.allTraj) {
      for(auto& tp : tj.Pts) {
        for(unsigned short iht = 0; iht < tp.Hits.size(); ++iht) {
          if(tp.Hits[iht] >= newHitIndex) ++tp.Hits[iht];
          
          if(tp.Hits.size() == 1 && tp.Chg > 0 && tjs.fHits[tp.Hits[iht]].WireID.Wire != std::nearbyint(tp.Pos[0])) {
            std::cout<<"  Create index problem tj.ID "<<tj.ID<<" iht "<<iht<<" newHitIndex "<<newHitIndex;
            std::cout<<" hit "<<PrintHit(tjs.fHits[tp.Hits[iht]])<<" Pos "<<PrintPos(tjs, tp)<<"\n";
            exit(1);
          }
          
        } // iht
      } // tp
    }

    return newHitIndex;
    
  } // CreateHit
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::MergeTPHits()
  {
    
    // Merge all hits that are used in one TP into a single hit
    // Make a list of hits that are slated for deletion
    std::vector<unsigned int> delHits;
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
      Trajectory& tj = tjs.allTraj[itj];
      // ignore shower Tj hits
      if(tj.AlgMod[kShowerTj]) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(NumHitsInTP(tp, kUsedHits) < 2) continue;
        // Make a list of the old hits on this TP before doing anything invasive
        std::vector<unsigned int> oldHits;
        // get some info so we can calculate the RMS
        raw::TDCtick_t loTick = INT_MAX;
        raw::TDCtick_t hiTick = 0;
        float mChg = 0;
        float mTick = 0;
        // estimate the uncertainties
        float mSigmaPeakAmp = 0;
        float mSigmaPeakTime = 0;
        float mSigmaIntegral = 0;
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          oldHits.push_back(iht);
          if(tjs.fHits[iht].StartTick < loTick) loTick = tjs.fHits[iht].StartTick;
          if(tjs.fHits[iht].EndTick > hiTick) hiTick = tjs.fHits[iht].EndTick;
          mChg += tjs.fHits[iht].Integral;
          mTick += tjs.fHits[iht].Integral * tjs.fHits[iht].PeakTime;
          mSigmaPeakAmp += tjs.fHits[iht].Integral * tjs.fHits[iht].SigmaPeakAmp;
          mSigmaPeakTime += tjs.fHits[iht].Integral * tjs.fHits[iht].SigmaPeakTime;
          mSigmaIntegral += tjs.fHits[iht].Integral * tjs.fHits[iht].SigmaIntegral;
        } // ii
        mTick /= mChg;
        if(mTick < 0) mTick = 0;
        mSigmaPeakAmp /= mChg;
        mSigmaPeakTime /= mChg;
        mSigmaIntegral /= mChg;
        // make a temporary signal waveform vector
        std::vector<float> signal(hiTick - loTick, 0);
        // fill it with the hit shapes
        for(auto& iht : oldHits) {
          float& peakTime = tjs.fHits[iht].PeakTime;
          float& amp = tjs.fHits[iht].PeakAmplitude;
          float& rms = tjs.fHits[iht].RMS;
          // add charge in the range +/- 3 sigma
          short loTime = (short)(peakTime - 3 * rms);
          if(loTime < loTick) loTime = loTick;
          short hiTime = (short)(peakTime + 3 * rms);
          if(hiTime > hiTick) hiTime = hiTick;
          for(short time = loTime; time < hiTime; ++time) {
            unsigned short indx = time - loTick;
            if(indx > signal.size() - 1) continue;
            float arg = (time - peakTime) / rms;
            signal[indx] += amp * exp(-0.5 * arg * arg);
          } // time
        } // iht
        // aveIndx is the index of the charge-weighted average in the signal vector
        float aveIndx = (mTick - loTick);
        // find the merged hit RMS
        float mRMS = 0;
        for(unsigned short indx = 0; indx < signal.size(); ++indx) {
          float dindx = indx - aveIndx;
          mRMS += signal[indx] * dindx * dindx;
        } // indx
        mRMS = std::sqrt(mRMS / mChg);
        // Modify the first hit in the list
        unsigned int mht = oldHits[0];
        tjs.fHits[mht].PeakTime = mTick;
        tjs.fHits[mht].SigmaPeakTime = mSigmaPeakTime;
        tjs.fHits[mht].PeakAmplitude = mChg / (2.5066 * mRMS);
        tjs.fHits[mht].SigmaPeakAmp = mSigmaPeakAmp;
        tjs.fHits[mht].Integral = mChg;
        tjs.fHits[mht].SigmaIntegral = mSigmaIntegral;
        tjs.fHits[mht].RMS = mRMS;
        tjs.fHits[mht].Multiplicity = 1;
        tjs.fHits[mht].LocalIndex = 0;
        tjs.fHits[mht].GoodnessOfFit = 1; // flag?
        tjs.fHits[mht].NDOF = 0;
        // then flag the other hits for erasing
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          for(unsigned short jj = 1; jj < oldHits.size(); ++jj) {
            if (tp.Hits[ii]==oldHits[jj]){
              tp.UseHit[ii] = false;
              // put it in the removal list
              delHits.push_back(tp.Hits[ii]);
              // Flag this hit
              tjs.fHits[tp.Hits[ii]].InTraj = SHRT_MAX;
              tp.Hits[ii] = INT_MAX;
            }
          }
        } // ii
      } // ipt
    } // itj
    
    // Erase the hits. Start by sorting them in decreasing order so that
    // the local delHits vector doesn't need to be modified when a hit is deleted
    if(delHits.size() > 1) std::sort(delHits.begin(), delHits.end(), std::greater<unsigned int>());

    for(auto& delHit : delHits) EraseHit(delHit);

  } // MergeTPHits
  
  //////////////////////////////////////////
  void TrajClusterAlg::MaskBadTPs(Trajectory& tj, float const& maxChi)
  {
    // Remove TPs that have the worst values of delta until the fit chisq < maxChi
    
    if(!tjs.UseAlg[kMaskBadTPs]) return;
    //don't use this function for reverse propagation
    if(!tjs.UseAlg[kRevProp]) return;
    
    if(tj.Pts.size() < 3) {
      mf::LogError("TC")<<"MaskBadTPs: Trajectory ID "<<tj.ID<<" too short to mask hits ";
      fGoodTraj = false;
      return;
    }
    unsigned short nit = 0;
    TrajPoint& lastTP = tj.Pts[tj.Pts.size() - 1];
    while(lastTP.FitChi > maxChi && nit < 3) {
      float maxDelta = 0;
      unsigned short imBad = USHRT_MAX;
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.Pts.size() - 1 - ii;
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        if(tp.Delta > maxDelta) {
          maxDelta = tp.Delta;
          imBad = ipt;
        }
        ++cnt;
        if(cnt == tp.NTPsFit) break;
      } // ii
      if(imBad == USHRT_MAX) return;
      if(prt) mf::LogVerbatim("TC")<<"MaskBadTPs: lastTP.FitChi "<<lastTP.FitChi<<"  Mask point "<<imBad;
      // mask the point
      UnsetUsedHits(tjs, tj.Pts[imBad]);
      FitTraj(tjs, tj);
      if(prt) mf::LogVerbatim("TC")<<"  after FitTraj "<<lastTP.FitChi;
      tj.AlgMod[kMaskBadTPs] = true;
      ++nit;
    } // lastTP.FItChi > maxChi && nit < 3
    
  } // MaskBadTPs
  
  //////////////////////////////////////////
  void TrajClusterAlg::MaskTrajEndPoints(Trajectory& tj, unsigned short nPts)
  {
    //PrintTrajectory("MTEP", tjs, tj, USHRT_MAX);

    // Masks off (sets all hits not-Used) nPts trajectory points at the leading edge of the
    // trajectory, presumably because the fit including this points is poor. The position, direction
    // and Delta of the last nPts points is updated as well
    
    if(tj.Pts.size() < 3) {
      mf::LogError("TC")<<"MaskTrajEndPoints: Trajectory ID "<<tj.ID<<" too short to mask hits ";
      fGoodTraj = false;
      return;
    }
    if(nPts > tj.Pts.size() - 2) {
      mf::LogError("TC")<<"MaskTrajEndPoints: Trying to mask too many points "<<nPts<<" Pts.size "<<tj.Pts.size();
      fGoodTraj = false;
      return;
    }
    
    // find the last good point (with charge)
    unsigned short lastGoodPt = USHRT_MAX ;

    if (!ChkMichel(tj, lastGoodPt)){ //did not find michel electron
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.EndPt[1] - nPts - ii;
        if(tj.Pts[ipt].Chg > 0) {
          lastGoodPt = ipt;
          break;
        }
        if(ipt == 0) break;
      } // ii
    }
    if(prt) {
      mf::LogVerbatim("TC")<<"MTEP: lastGoodPt "<<lastGoodPt<<" Pts size "<<tj.Pts.size()<<" fGoodTraj "<<fGoodTraj;
    }
    if(lastGoodPt == USHRT_MAX) return;
    tj.EndPt[1] = lastGoodPt;
    
    //for(unsigned short ii = 0; ii < nPts; ++ii) {
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      if (ipt==lastGoodPt) break;
      UnsetUsedHits(tjs, tj.Pts[ipt]);
      // Reset the position and direction of the masked off points
      tj.Pts[ipt].Dir = tj.Pts[lastGoodPt].Dir;
      if(tj.Pts[lastGoodPt].AngleCode == 2) {
        // Very large angle: Move by path length
        float path = TrajPointSeparation(tj.Pts[lastGoodPt], tj.Pts[ipt]);
        tj.Pts[ipt].Pos[0] = tj.Pts[lastGoodPt].Pos[0] + path * tj.Pts[ipt].Dir[0];
        tj.Pts[ipt].Pos[1] = tj.Pts[lastGoodPt].Pos[1] + path * tj.Pts[ipt].Dir[1];
      } else {
        // Not large angle: Move by wire
        float dw = tj.Pts[ipt].Pos[0] - tj.Pts[lastGoodPt].Pos[0];
        // Correct the projected time to the wire
        float newpos = tj.Pts[lastGoodPt].Pos[1] + dw * tj.Pts[ipt].Dir[1] / tj.Pts[ipt].Dir[0];
        if(prt) mf::LogVerbatim("TC")<<"MTEP: ipt "<<ipt<<" Pos[0] "<<tj.Pts[ipt].Pos[0]<<". Move Pos[1] from "<<tj.Pts[ipt].Pos[1]<<" to "<<newpos;
        tj.Pts[ipt].Pos[1] = tj.Pts[lastGoodPt].Pos[1] + dw * tj.Pts[ipt].Dir[1] / tj.Pts[ipt].Dir[0];
      }
      tj.Pts[ipt].Delta = PointTrajDOCA(tjs, tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tj.Pts[ipt]);
      if(prt) mf::LogVerbatim("TC")<<" masked ipt "<<ipt<<" Pos "<<PrintPos(tjs, tj.Pts[ipt])<<" Chg "<<tj.Pts[ipt].Chg;
    } // ii
    SetEndPoints(tjs, tj);
    
  } // MaskTrajEndPoints
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ChkStop(Trajectory& tj)
  {
    // Sets the StopFlag[kBragg] bits on the trajectory by identifying the Bragg peak
    // at each end. This function checks both ends, finding the point with the highest charge nearest the
    // end and considering the first (when end = 0) 4 points or last 4 points (when end = 1). The next
    // 5 - 10 points (fChkStop[0]) are fitted to a line, Q(x - x0) = Qo + (x - x0) * slope where x0 is the
    // wire position of the highest charge point. A large negative slope indicates that there is a Bragg
    // peak at the end.
    
    tj.StopFlag[0][kBragg] = false;
    tj.StopFlag[1][kBragg] = false;
    
    if(fChkStopCuts[0] < 0) return;
    
    // don't attempt with low momentum trajectories
    if(tj.MCSMom < 50) return;
    
    // ignore trajectories that are very large angle at both ends
    if(tj.Pts[tj.EndPt[0]].AngleCode == 2 || tj.Pts[tj.EndPt[1]].AngleCode == 2) return;

    unsigned short nPtsToCheck = fChkStopCuts[1];
    if(tj.Pts.size() < nPtsToCheck) return;
    
    if(prt) mf::LogVerbatim("TC")<<"ChkStop: requiring "<<nPtsToCheck<<" points with charge slope > "<<fChkStopCuts[0]<<" Chg/WSEU";
    
    // find the highest charge hit in the first 3 points at each end
    for(unsigned short end = 0; end < 2; ++end) {
      short dir = 1 - 2 * end;
      // find the point with the highest charge considering the first 3 points
      float big = 0;
      unsigned short hiPt = 0;
      float x0 = 0;
      for(unsigned short ii = 0; ii < 4; ++ii) {
        short ipt = tj.EndPt[end] + ii * dir;
        if(ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) break;
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg > big) {
          big = tp.Chg;
          x0 = tp.Pos[0];
          hiPt = ipt;
        }
      } // ii
      if(prt) mf::LogVerbatim("TC")<<" end "<<end<<" hiPt "<<hiPt<<" big "<<big;
      std::vector<float> x, y, yerr2;
      float intcpt, intcpterr;
      float slope, slopeerr, chidof;
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        short ipt = hiPt + ii * dir;
        if(ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) break;
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        x.push_back(std::abs(tp.Pos[0] - x0));
        y.push_back(tp.Chg);
        // Assume 10% point-to-point charge fluctuations
        float err = 0.1 * tp.Chg;
        if(prt) mf::LogVerbatim("TC")<<ipt<<"  "<<PrintPos(tjs, tp.Pos)<<" "<<x[x.size()-1]<<" Chg "<<(int)tp.Chg;
        yerr2.push_back(err * err);
        if(x.size() == nPtsToCheck) break;
      } // ii
      if(x.size() < nPtsToCheck) continue;
      fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
      // check for really bad chidof indicating a major failure
      if(chidof > 100) continue;
      // The charge slope is negative for a stopping track in the way that the fit was constructed.
      // Flip the sign so we can make a cut against fChkStopCuts[0] which is positive.
      slope = -slope;
      if(slope > fChkStopCuts[0] && chidof < fChkStopCuts[2] && slope > 2 * slopeerr) {
        tj.StopFlag[end][kBragg] = true;
        tj.AlgMod[kChkStop] = true;
        if(prt) mf::LogVerbatim("TC")<<" end "<<end<<" fit chidof "<<chidof<<" slope "<<slope<<" +/- "<<slopeerr<<" Stopping ";
      } else {
        if(prt) mf::LogVerbatim("TC")<<" end "<<end<<" fit chidof "<<chidof<<" slope "<<slope<<" +/- "<<slopeerr<<" Not stopping";
      }
   } // end

  } // ChkStop

  //////////////////////TY://////////////////////////
  bool TrajClusterAlg::ChkMichel(Trajectory& tj, unsigned short& lastGoodPt){

    //find number of hits that are consistent with Michel electron
    unsigned short nmichelhits = 0;
    //find number of hits that are consistent with Bragg peak
    unsigned short nbragghits = 0;
    float lastChg = 0;

    bool isfirsthit = true;
    unsigned short braggpeak = 0;

    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      if (ii>tj.EndPt[1]) continue;
      unsigned short ipt = tj.EndPt[1] - ii;
      if (tj.Pts[ipt].Chg>0){
        if (isfirsthit){
          isfirsthit = false;
          if (tj.Pts[ipt].ChgPull<0){
            ++nmichelhits;
          }
        }
        else{
          if (tj.Pts[ipt].ChgPull<0&&nmichelhits&&!nbragghits){//still Michel
            ++nmichelhits;
          }
          else{
            if (!nbragghits){
              ++nbragghits; //Last Bragg peak hit
              lastChg  = tj.Pts[ipt].Chg;
              braggpeak = ipt;
            }
            else if (tj.Pts[ipt].Chg<lastChg){ //still Bragg peak
              ++nbragghits;
              lastChg  = tj.Pts[ipt].Chg;
            }
            else break;
          }
        }
      }
    }
    if(prt) mf::LogVerbatim("TC")<<"ChkMichel Michel hits: "<<nmichelhits<<" Bragg peak hits: "<<nbragghits;
    if (nmichelhits>0&&nbragghits>2){//find Michel topology
      lastGoodPt = braggpeak;
      return true;
    }
    else{
      return false;
    }
  }

  ////////////////////////////////////////////////
  void TrajClusterAlg::ChkHiChgHits()
  {
    // Check allTraj trajectories in the current CTP to see if they are stopping
    if(!tjs.UseAlg[kSplitHiChgHits]) return;
    
    for(size_t i = 0; i< tjs.allTraj.size(); ++i) {
      auto & tj = tjs.allTraj[i];
      if(tj.CTP != fCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      SplitHiChgHits(tj);
    } // tj

  } // ChkHiChgHits

  /////////////////////TY:///////////////////////////
  void TrajClusterAlg::SplitHiChgHits(Trajectory& tj){
    
    // Check allTraj trajectories in the current CTP and split high charge hits 
    if(!tjs.UseAlg[kSplitHiChgHits]) return;

    // Only do it once
    if (tj.AlgMod[kSplitHiChgHits]) return;

    if(tj.CTP != fCTP) return;
    if(tj.AlgMod[kKilled]) return;
    //Ignore short trajectories
    if (tj.EndPt[1]<10) return;
    for(unsigned short end = 0; end < 2; ++end) {
      if(prt) mf::LogVerbatim("TC")<<"SplitHiChghits "<<end<<" "<<tj.VtxID[end];
      float hichg = 0;
      unsigned short tp = tj.EndPt[end];
      unsigned short nlohits = 0;
      unsigned short lastHiTP = USHRT_MAX;
      while (tp != tj.EndPt[1-end]){
        float ptchg = TpSumHitChg(tjs, tj.Pts[tp]);
        if (prt) mf::LogVerbatim("TC")<<"SplitHiChgHits "<<tp<<" "<<ptchg<<" "<<PrintPos(tjs, tj.Pts[tp]);
        if (ptchg){
          if (tp == tj.EndPt[end]){
            hichg = ptchg;
            lastHiTP = tp;
          }
          else if (ptchg>0.4*hichg){
            if (!nlohits){
              hichg = ptchg;
              lastHiTP = tp;
            }
            else{
              break;
            }
          }
          else ++nlohits;
        }
        if (end==0){
          ++tp;
        }
        else{
          --tp;
        }
      }
      //if (prt) mf::LogVerbatim("TC")<<"SplitHiChgHits "<<end<<" "<<nlohits;
      if (nlohits>4&&lastHiTP!=USHRT_MAX){
        //Create new vertex
        VtxStore aVtx;
        aVtx.Pos = tj.Pts[lastHiTP].Pos;
        aVtx.NTraj = 2;
        aVtx.Pass = tj.Pass;
        aVtx.Topo = 7;
        aVtx.ChiDOF = 0;
        aVtx.CTP = fCTP;
        aVtx.ID = tjs.vtx.size() + 1;
        if(!StoreVertex(tjs, aVtx)) {
          if(prt) mf::LogVerbatim("TC")<<" Failed storing vertex "<<tj.VtxID[end];
          return;
        }

        // make a copy
        Trajectory newTj = tj;
        newTj.ID = tjs.allTraj.size() + 1;

        // keep high charge hits, reassign other hits to the new trajectory
        unsigned short tp1 = lastHiTP+1;
        if (end==1) tp1 = lastHiTP-1;
        for (unsigned short ipt = std::min(tj.EndPt[1-end], tp1); ipt <= std::max(tj.EndPt[1-end], tp1); ++ipt){
          tj.Pts[ipt].Chg = 0;
          for (unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
            if(!tj.Pts[ipt].UseHit[ii]) continue;
            unsigned int iht = tj.Pts[ipt].Hits[ii];
            // This shouldn't happen but check anyway
            if(tjs.fHits[iht].InTraj != tj.ID) continue;
            tjs.fHits[iht].InTraj = newTj.ID;
            tj.Pts[ipt].UseHit[ii] = false;
          }//ii
        }//ipt
        SetEndPoints(tjs, tj);
        tj.VtxID[1-end] = aVtx.ID;
        tj.AlgMod[kSplitHiChgHits] = true;
        if(prt) {
          mf::LogVerbatim("TC")<<"Splitting trajectory ID "<<tj.ID<<" new EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
        }

        for (unsigned short ipt = std::min(newTj.EndPt[end], lastHiTP); ipt <= std::max(newTj.EndPt[end], lastHiTP); ++ipt){
          newTj.Pts[ipt].Chg = 0;
          for (unsigned short ii = 0; ii < newTj.Pts[ipt].Hits.size(); ++ii) {
            newTj.Pts[ipt].UseHit[ii] = false;
          }//ii
        }//ipt
        SetEndPoints(tjs, newTj);
        newTj.VtxID[end] = aVtx.ID;
        newTj.AlgMod[kSplitHiChgHits] = true;
        tjs.allTraj.push_back(newTj);
        
        break;     
      }
    }
  }


} // namespace cluster
