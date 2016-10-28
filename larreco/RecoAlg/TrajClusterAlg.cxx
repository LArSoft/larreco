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


// TEMP for MatchTruth
#include "larsim/MCCheater/BackTracker.h"

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
  TrajClusterAlg::TrajClusterAlg(fhicl::ParameterSet const& pset)
  {
    reconfigure(pset);
    

    art::ServiceHandle<art::TFileService> tfs;
    
    // True - Reco vertex difference
    fNuVtx_dx = tfs->make<TH1F>("Vtx dx","Vtx dx",80,-10,10);
    fNuVtx_dy = tfs->make<TH1F>("Vtx dy","Vtx dy",80,-10,10);
    fNuVtx_dz = tfs->make<TH1F>("Vtx dz","Vtx dz",80,-10,10);
    
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
    
    fDeltaN[0] = tfs->make<TH1F>("DeltaN0","Normalized Delta Pln 0", 50, 0, 4);
    fDeltaN[1] = tfs->make<TH1F>("DeltaN1","Normalized Delta Pln 1", 50, 0, 4);
    fDeltaN[2] = tfs->make<TH1F>("DeltaN2","Normalized Delta Pln 2", 50, 0, 4);

    fHitRMS[0] = tfs->make<TH1F>("hitrms0","Hit RMS Pln 0", 80, 0, 20);
    fHitRMS[1] = tfs->make<TH1F>("hitrms1","Hit RMS Pln 1", 80, 0, 20);
    fHitRMS[2] = tfs->make<TH1F>("hitrms2","Hit RMS Pln 2", 80, 0, 20);

    fTPWidth_Angle[0] = tfs->make<TH2F>("tpwidth_angle0","TP hit width vs Angle Pln 0", 20, 0, M_PI/2, 20, 0, 200);
    fTPWidth_Angle[1] = tfs->make<TH2F>("tpwidth_angle1","TP hit width vs Angle Pln 1", 20, 0, M_PI/2, 20, 0, 200);
    fTPWidth_Angle[2] = tfs->make<TH2F>("tpwidth_angle2","TP hit width vs Angle Pln 2", 20, 0, M_PI/2, 20, 0, 200);

    fTPWidth_AngleP[0] = tfs->make<TProfile>("tpwidth_anglep0","TP hit width vs Angle Pln 0", 10, 0, M_PI/2, "S");
    fTPWidth_AngleP[1] = tfs->make<TProfile>("tpwidth_anglep1","TP hit width vs Angle Pln 1", 10, 0, M_PI/2, "S");
    fTPWidth_AngleP[2] = tfs->make<TProfile>("tpwidth_anglep2","TP hit width vs Angle Pln 2", 10, 0, M_PI/2, "S");

    fExpect_Angle[0] = tfs->make<TProfile>("expect_angle0","Expected width vs Angle Pln 0", 11, 0, M_PI/2, "S");
    fExpect_Angle[1] = tfs->make<TProfile>("expect_angle1","Expected width vs Angle Pln 1", 11, 0, M_PI/2, "S");
    fExpect_Angle[2] = tfs->make<TProfile>("expect_angle2","Expected width vs Angle Pln 2", 11, 0, M_PI/2, "S");
    
    for(unsigned short pdgIndex = 0; pdgIndex < 6; ++pdgIndex) {
      EPCounts[pdgIndex] = 0;
      EPSums[pdgIndex] = 0;
    }
    fEventsProcessed = 0;
//    if(fStudyMode) outFile.open("quality.txt");
    
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
    fAngleRanges          = pset.get< std::vector<float>>("AngleRanges");
    fNPtsAve              = pset.get< short >("NPtsAve", 20);
    fMinPtsFit            = pset.get< std::vector<unsigned short >>("MinPtsFit");
    fMinPts               = pset.get< std::vector<unsigned short >>("MinPts");
    fMaxAngleRange        = pset.get< std::vector<unsigned short>>("MaxAngleRange");
    fMaxChi               = pset.get< float >("MaxChi", 10);
    fChgPullCut           = pset.get< float >("ChgPullCut", 3);
    fMultHitSep           = pset.get< float >("MultHitSep", 2.5);
    fKinkAngCut           = pset.get< float >("KinkAngCut", 0.4);
    fMaxWireSkipNoSignal  = pset.get< float >("MaxWireSkipNoSignal", 1);
    fMaxWireSkipWithSignal= pset.get< float >("MaxWireSkipWithSignal", 100);
    fProjectionErrFactor  = pset.get< float >("ProjectionErrFactor", 2);
    fJTMaxHitSep2         = pset.get< float >("JTMaxHitSep", 2);
    
    std::vector<std::string> skipAlgsVec = pset.get< std::vector<std::string>  >("SkipAlgs");
    
    fTagAllTraj           = pset.get< bool  >("TagAllTraj", false);
    fDeltaRayTag          = pset.get< std::vector<short>>("DeltaRayTag", {-1, -1, -1});
    fMuonTag              = pset.get< std::vector<short>>("MuonTag", {-1, -1, -1, - 1});
    fShowerTag            = pset.get< std::vector<short>>("ShowerTag", {-1, -1});
    fChkStopCuts          = pset.get< std::vector<float>>("ChkStopCuts", {-1, -1, -1});
    fMaxTrajSep           = pset.get< float >("MaxTrajSep", 4);
    
    fStudyMode            = pset.get< bool  >("StudyMode", false);
    fMatchTruth           = pset.get< std::vector<float> >("MatchTruth", {-1, -1, -1, -1});
    fVertex2DCuts         = pset.get< std::vector<float >>("Vertex2DCuts", {-1, -1, -1, -1, -1});
    fVertex3DChiCut       = pset.get< float >("Vertex3DChiCut", -1);
    fMaxVertexTrajSep     = pset.get< std::vector<float>>("MaxVertexTrajSep");
    
    debug.Plane           = pset.get< int >("DebugPlane", -1);
    debug.Wire            = pset.get< int >("DebugWire", -1);
    debug.Tick            = pset.get< int >("DebugTick", -1);
    debug.WorkID          = pset.get< short>("DebugWorkID", 0);
    
    // check the angle ranges and convert from degrees to radians
    if(fAngleRanges.back() < 90) {
      std::cout<<"Last element of AngleRange != 90 degrees. Fixing it\n";
      fAngleRanges.back() = 90;
    }
    for(auto& range : fAngleRanges) {
      if(range < 0 || range > 90) throw art::Exception(art::errors::Configuration)<< "Invalid angle range "<<range<<" Must be 0 - 90 degrees";
      range *= M_PI / 180;
    } // range
    // size this and set it later
    fAngleRangesMaxHitsRMS.resize(fAngleRanges.size());

    // convert the max traj separation into a separation^2
    fMaxTrajSep *= fMaxTrajSep;
    if(fJTMaxHitSep2 > 0) fJTMaxHitSep2 *= fJTMaxHitSep2;
    
    if(fMinPtsFit.size() != fMinPts.size()) badinput = true;
    if(fMaxVertexTrajSep.size() != fMinPts.size()) badinput = true;
    if(fMaxAngleRange.size() != fMinPts.size()) badinput = true;
    if(badinput) throw art::Exception(art::errors::Configuration)<< "Bad input from fcl file. Vector lengths for MinPtsFit, MaxVertexTrajSep and MaxAngleRange are not the same";
    
    if(fVertex2DCuts.size() < 5) throw art::Exception(art::errors::Configuration)<<"Vertex2DCuts must be size 5";
    
    if(fMuonTag.size() < 4) throw art::Exception(art::errors::Configuration)<<"MuonTag must be size 4 [minPtsFit, minMCSMom, maxWireSkipNoSignal, minDeltaLen]";
    if(fDeltaRayTag.size() < 3) throw art::Exception(art::errors::Configuration)<<"DeltaRayTag must be size 3 [max endpoint sep, min MCSMom, max MCSMom]";
    if(fChkStopCuts.size() < 3) throw art::Exception(art::errors::Configuration)<<"ChkStopCuts must be size 3 [Min Chg ratio, Chg slope pull cut, Chg fit chi cut]";
    if(fShowerTag.size() < 3) throw art::Exception(art::errors::Configuration)<< "ShowerTag must be size 3";
    
    if(kAlgBitSize != AlgBitNames.size())
      throw art::Exception(art::errors::Configuration)<<"kAlgBitSize "<<kAlgBitSize<<" != AlgBitNames size "<<AlgBitNames.size();
    fAlgModCount.resize(kAlgBitSize);
    
    unsigned short ib;
    bool gotit, printHelp = false;
    for(ib = 0; ib < AlgBitNames.size(); ++ib) fUseAlg[ib] = true;
    for(auto strng : skipAlgsVec) {
      gotit = false;
      if(strng == "All") {
        // turn everything off
        for(ib = 0; ib < AlgBitNames.size(); ++ib) fUseAlg[ib] = false;
        std::cout<<"Turning all algs off\n";
        gotit = true;
        break;
      } // All off
      for(ib = 0; ib < AlgBitNames.size(); ++ib) {
        if(strng == AlgBitNames[ib]) {
          fUseAlg[ib] = false;
          gotit = true;
          break;
        }
      } // ib
      if(!gotit) {
        std::cout<<"Unknown SkipAlgs input string '"<<strng<<"'\n";
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
    // Change the polarity of ChkInTraj
    if(fUseAlg[kChkInTraj]) { fUseAlg[kChkInTraj] = false; } else { fUseAlg[kChkInTraj] = true; std::cout<<"Note: ChkInTraj will be slow...\n"; }
   
  } // reconfigure
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ClearResults() {
    // clear everything in reverse order in which it appears in tjs DataStructs.h
    // so that vectors that appear later in the list are not unnecessarily shifted
    // before they are cleared.
    tjs.vtx3.clear();
    tjs.vtx.clear();
    tjs.tcl.clear();
    tjs.inClus.clear();
    tjs.WireHitRange.clear();
    tjs.fHits.clear();
    tjs.allTraj.clear();
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
 
    larprop = lar::providerFrom<detinfo::LArPropertiesService>();
    detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    tjs.fHits.reserve(hitVecHandle->size());

    for(unsigned int iht = 0; iht < hitVecHandle->size(); ++iht) {
      art::Ptr<recob::Hit> hit = art::Ptr<recob::Hit>(hitVecHandle, iht);
      if(hit->PeakAmplitude() < fMinAmp) continue;
      TCHit localHit;
      localHit.StartTick = hit->StartTick();
      localHit.EndTick = hit->EndTick();
      localHit.PeakTime = hit->PeakTime(); // This will be converted to Time later
      localHit.PeakAmplitude = hit->PeakAmplitude(); // This will be converted to Time later
      localHit.Integral = hit->Integral();
      localHit.RMS = hit->RMS();
      localHit.GoodnessOfFit = hit->GoodnessOfFit();
      localHit.NDOF = hit->DegreesOfFreedom();
      localHit.Multiplicity = hit->Multiplicity();
      localHit.LocalIndex = hit->LocalIndex();
      localHit.WireID = hit->WireID();
      tjs.fHits.push_back(localHit);
    } // iht
    // set a flag so that we do the tick -> time conversion later and only once
//    tjs.ConvertTicksToTime = true;

    // sort it as needed;
    // that is, sorted by wire ID number,
    // then by start of the region of interest in time, then by the multiplet
    std::sort(tjs.fHits.begin(), tjs.fHits.end(), &SortByMultiplet);
    
    // check for debugging mode triggered by Plane, Wire, Tick
    if(debug.Plane >= 0 && debug.Plane < 3 && debug.WorkID >= 0 && debug.Wire > 0 && debug.Tick > 0) {
      std::cout<<"Looking for debug hit "<<debug.Plane<<":"<<debug.Wire<<":"<<debug.Tick;
      for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
        if((int)tjs.fHits[iht].WireID.Plane != debug.Plane) continue;
        if((int)tjs.fHits[iht].WireID.Wire != debug.Wire) continue;
        if(tjs.fHits[iht].PeakTime < debug.Tick - 5) continue;
        if(tjs.fHits[iht].PeakTime > debug.Tick + 5) continue;
        debug.Hit = iht;
        std::cout<<" iht "<<iht<<" "<<debug.Plane<<":"<<PrintHit(tjs.fHits[iht]);
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
    
    fStepDir = fMode;
    InitializeAllTraj();
    for (geo::TPCID const& tpcid: geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      FillWireHitRange(tpcid);
      if(fQuitAlg) return;
      for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
        // no hits on this plane?
        if(tjs.FirstWire[fPlane] > tjs.LastWire[fPlane]) continue;
        // Set the CTP code to ensure objects are compared within the same plane
        fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
        fCstat = tpcid.Cryostat;
        fTpc = tpcid.TPC;
        // reconstruct all trajectories in the current plane
        ReconstructAllTraj();
        if(fQuitAlg) {
          std::cout<<"RunTrajCluster failed in ReconstructAllTraj Run "<<fRun<<" Event "<<fEvent<<" EventsProcessed "<<fEventsProcessed<<"\n";
          mf::LogVerbatim("TC")<<"RunTrajCluster failed after ReconstructAllTraj";
          ClearResults();
          return;
        }
      } // fPlane
      // No sense taking muon direction if delta ray tagging is disabled
      if(fDeltaRayTag[0] >= 0) TagMuonDirections(tjs, fMuonTag[3], debug.WorkID);
      if(fVertex3DChiCut > 0) Find3DVertices(tpcid);
    } // tpcid

    MatchTruth();
 
    // Convert trajectories in allTraj into clusters
    MakeAllTrajClusters();
    if(fQuitAlg) {
      std::cout<<"RunTrajCluster failed in MakeAllTrajClusters\n";
      ClearResults();
      return;
    }
    CheckHitClusterAssociations();
    if(fQuitAlg) {
      ClearResults();
      std::cout<<"RunTrajCluster failed in CheckHitClusterAssociations\n";
      return;
    }

    if(TJPrt > 0 || debug.Plane >= 0) {
      mf::LogVerbatim("TC")<<"Done in RunTrajClusterAlg";
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        if(tjs.allTraj[itj].WorkID == TJPrt) {
          PrintAllTraj("DBG", tjs, debug, itj, USHRT_MAX);
          break;
        }
      } // itj
      // Print summary of all TJs
      PrintAllTraj("RTC", tjs, debug, USHRT_MAX, 0);
    }

    // temp
    if(fStudyMode) {
      for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
        if(tjs.fHits[iht].Multiplicity != 1) continue;
        if(tjs.fHits[iht].GoodnessOfFit < 0) continue;
        if(tjs.fHits[iht].GoodnessOfFit > 50) continue;
        unsigned short ipl = tjs.fHits[iht].WireID.Plane;
        fHitRMS[ipl]->Fill(tjs.fHits[iht].RMS);
      } // iht
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.AlgMod[kKilled]) continue;
        if(tj.MCSMom == 0) continue;
/*
        // TP hit width plots
        unsigned short ipl = tj.CTP;
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          TrajPoint& tp = tj.Pts[ipt];
          if(tp.Chg == 0) continue;
          float dang = tp.Ang;
          if(dang > M_PI) dang = M_PI;
          if(dang < -M_PI) dang = M_PI;
          if(dang < 0) dang = -dang;
          if(dang > M_PI/2) dang = M_PI - dang;
          // width of all used hits in this tp
          float hitWid = TPHitsRMSTick(tjs, tp, kAllHits);
          fTPWidth_Angle[ipl]->Fill(dang, hitWid);
          fTPWidth_AngleP[ipl]->Fill(dang, hitWid);
          float expect = fAveHitRMS[ipl];
          if(std::abs(tp.Dir[0]) > 0.001) expect += std::abs(tp.Dir[1]/tp.Dir[0])/tjs.UnitsPerTick;
          fExpect_Angle[ipl]->Fill(dang, expect);
//          std::cout<<dang<<" "<<hitWid<<" "<<expect<<" Dir0 "<<tp.Dir[0]<<"\n";
          float dn = tp.Delta / fHitErrFac;
          if(dn > 0) fDeltaN[ipl]->Fill(dn);
        } // ipt
        if(tj.TruKE == 0) continue;
        unsigned short pdg = std::abs(tj.TruPDG);
        double mass = 0.511;
        if(pdg == 13) mass = 105.7;
        if(pdg == 211) mass = 139.6;
        if(pdg == 2212) mass = 938.3;
        double tPlusM = tjs.allTraj[itj].TruKE + mass;
        double truMom = sqrt(tPlusM * tPlusM - mass * mass);
        std::cout<<tj.CTP<<":"<<PrintPos(tjs, tj.Pts[tj.EndPt[0]])<<"-"<<tj.CTP<<":"<<PrintPos(tjs, tj.Pts[tj.EndPt[1]]);
        std::cout<<" MCS TruKE "<<tj.TruKE<<" pdg "<<tj.TruPDG<<" MCSMom "<<(int)tj.MCSMom<<" EffPur "<<std::setprecision(2)<<tj.EffPur<<"\n";
        if(pdg == 11) fMCSMom_KE_e->Fill(tj.TruKE, tj.MCSMom);
        if(pdg == 13) fMCSMom_KE_mu->Fill(tj.TruKE, tj.MCSMom);
        if(pdg == 211) fMCSMom_KE_pi->Fill(tj.TruKE, tj.MCSMom);
        if(pdg == 2212) fMCSMom_KE_p->Fill(tj.TruKE, tj.MCSMom);
 */
      } // itj
    } // studymode
    
    if(fMatchTruth[0] >= 0) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Event "<<evt.event();
      float sum = 0;
      float cnt = 0;
      for(unsigned short pdgIndex = 0; pdgIndex < 6; ++pdgIndex) {
        if(EPCounts[pdgIndex] == 0) continue;
        if(pdgIndex == 0) myprt<<" Electron";
        if(pdgIndex == 1) myprt<<" Muon";
        if(pdgIndex == 2) myprt<<" Pion";
        if(pdgIndex == 3) myprt<<" Kaon";
        if(pdgIndex == 4) myprt<<" Proton";
        float ave = EPSums[pdgIndex] / (float)EPCounts[pdgIndex];
        myprt<<" ave "<<std::fixed<<std::setprecision(2)<<ave;
        myprt<<" cnt "<<EPCounts[pdgIndex];
        sum += EPSums[pdgIndex];
        cnt += EPCounts[pdgIndex];
      } // pdgIndex
      if(cnt > 0) myprt<<" combined "<<std::fixed<<std::setprecision(2)<<sum / cnt;
    } // fMatchTruth[0] >= 0

    // convert vertex time from WSE to ticks
    for(auto& avtx : tjs.vtx) avtx.Pos[1] /= tjs.UnitsPerTick;
    
    std::cout<<"RunTrajCluster success run "<<fRun<<" event "<<fEvent<<" allTraj size "<<tjs.allTraj.size()<<" events processed "<<fEventsProcessed<<"\n";
    
  } // RunTrajClusterAlg

  ////////////////////////////////////////////////
  void TrajClusterAlg::InitializeAllTraj()
  {
    tjs.allTraj.clear();
    tjs.vtx.clear();
    tjs.vtx3.clear();
  } // InitializeAllTraj
  ////////////////////////////////////////////////
  bool TrajClusterAlg::MergeAndStore(unsigned short itj1, unsigned short itj2)
  {
    // Merge the two trajectories and store them. Returns true if it was successfull.
    // Merging is done between the end of tj1 and the beginning of tj2
    // First check for major failures
    fQuitAlg = false;
    if(itj1 > tjs.allTraj.size() - 1) fQuitAlg = true;
    if(itj2 > tjs.allTraj.size() - 1) fQuitAlg = true;
    if(tjs.allTraj[itj1].AlgMod[kKilled] || tjs.allTraj[itj2].AlgMod[kKilled]) {
      mf::LogWarning("TC")<<"MergeAndStore: Trying to merge a killed trajectory. Here they are ";
      PrintAllTraj("tj1", tjs, debug, itj1, USHRT_MAX);
      PrintAllTraj("tj1", tjs, debug, itj2, USHRT_MAX);
      fQuitAlg = true;
    }
    
    if(fQuitAlg) {
      mf::LogError("TC")<<"Failed in MergeAndStore";
      return false;
    }
    
    // make a copy so they can be trimmed as needed
    Trajectory tj1 = tjs.allTraj[itj1];
    Trajectory tj2 = tjs.allTraj[itj2];
    
    if(mrgPrt) {
      mf::LogVerbatim("TC")<<"MergeAngStore: tj1.ID "<<tj1.ID<<" tj2.ID "<<tj2.ID;
    }
    
    // ensure that these are in the same step order
    if(tj1.StepDir != tj2.StepDir) return false;
    
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
    
    // determine if some re-fitting of one of the trajectories at one end is required
    // mergePt will be the first point in the new trajectory that belonged to tj2
    unsigned short mergePt = tj1.Pts.size();
    bool fixTj1 = (tj1.Pts.size() < 10 && tj2.Pts.size() > 10);
    bool fixTj2 = (tj1.Pts.size() > 10 && tj2.Pts.size() < 10);
    // append tj2 hits to tj1

    tj1.Pts.insert(tj1.Pts.end(), tj2.Pts.begin() + tj2ClosePt, tj2.Pts.end());
    // re-define the end points
    SetEndPoints(tjs, tj1);
    
    // A more exhaustive check that hits only appear once
    if(HasDuplicateHits(tj1)) return false;
    
    // Correct the trajectory in the merging region if one of the tjs is short
    if(fixTj1) {
      // tj1 is short so extrapolate tj2 into tj1
      // We can use FixTrajBegin for this
      prt = mrgPrt;
      if(prt) mf::LogVerbatim("TC")<<"MergeAngStore: Fix Begin of tj1.ID "<<tj1.ID;
      FixTrajBegin(tj1, mergePt);
      if(mrgPrt) prt = false;
    } else if (fixTj2) {
      // tj2 is short so extrapolate tj1 into tj2
      --mergePt;
      prt = mrgPrt;
      if(prt) mf::LogVerbatim("TC")<<"MergeAngStore: Fix End of tj1.ID "<<tj1.ID;
      FixTrajEnd(tj1, mergePt);
      if(mrgPrt) prt = false;
    }

    // kill the original trajectories
    MakeTrajectoryObsolete(tjs, itj1);
    MakeTrajectoryObsolete(tjs, itj2);
    // Do this so that StoreTraj keeps the correct WorkID (of itj1)
    tj1.ID = tj1.WorkID;
    StoreTraj(tj1);
    return true;
    
  } // MergeAndStore

  ////////////////////////////////////////////////
  void TrajClusterAlg::CountSameHits(std::vector<unsigned int>& iHitVec, std::vector<unsigned int>& jHitVec, unsigned short& nSameHits)
  {
    // Counts the number of hits that are used in two different vectors of hits
    nSameHits = 0;
    for(unsigned short ii = 0; ii < iHitVec.size(); ++ii) {
      if(std::find(jHitVec.begin(), jHitVec.end(), iHitVec[ii]) != jHitVec.end()) ++nSameHits;
    }
  } // CountSameHits
  
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
          if(prt) didPrt = true;
          if(prt)  mf::LogVerbatim("TC")<<"+++++++ Pass "<<pass<<" Found debug hit "<<fPlane<<":"<<PrintHit(tjs.fHits[iht]);
          // ignore below-threshold hits
          // clear out any leftover work tjs.inTraj's that weren't cleaned up properly
//          for(oht = ifirsthit; oht < ilasthit; ++oht) if(tjs.inTraj[oht] < 0) tjs.inTraj[oht] = 0;
          // Only consider hits that are available
          if(tjs.fHits[iht].InTraj != 0) continue;
          // We hope to make a trajectory point at the hit position of iht in WSE units
          // with a direction pointing to jht
          unsigned int fromWire = tjs.fHits[iht].WireID.Wire;
          float fromTick = tjs.fHits[iht].PeakTime;
          float iqtot = tjs.fHits[iht].Integral;
          float hitsRMSTick = tjs.fHits[iht].RMS;
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
          if(prt) mf::LogVerbatim("TC")<<" RMS "<<tjs.fHits[iht].RMS<<" BB Multiplicity "<<iHitsInMultiplet.size()<<" AveHitRMS["<<fPlane<<"] "<<fAveHitRMS[fPlane]<<" HitsRMSTick "<<hitsRMSTick<<" fatIHit "<<fatIHit;
          for(jht = jfirsthit; jht < jlasthit; ++jht) {
            // Only consider hits that are available
            if(tjs.fHits[iht].InTraj != 0) continue;
            if(tjs.fHits[jht].InTraj != 0) continue;
            // clear out any leftover work inTraj's that weren't cleaned up properly
            for(unsigned short oht = jfirsthit; oht < jlasthit; ++oht) {
              if(tjs.fHits[oht].InTraj < 0) {
                mf::LogVerbatim("TC")<<"Bad cleanup "<<PrintHit(tjs.fHits[oht])<<" "<<tjs.fHits[oht].InTraj<<" events processed "<<fEventsProcessed;
                std::cout<<"Bad cleanup "<<PrintHit(tjs.fHits[oht])<<" "<<tjs.fHits[oht].InTraj<<" events processed "<<fEventsProcessed<<" fWorkID "<<fWorkID<<"\n";
                tjs.fHits[oht].InTraj = 0;
//                fQuitAlg = true;
//                return;
              }
            }
            unsigned int toWire = jwire;
            float toTick = tjs.fHits[jht].PeakTime;
            float jqtot = tjs.fHits[jht].Integral;
            if(jqtot < 1) continue;
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
//              HitMultipletPosition(jht, toTick, deltaRms, jqtot);
            }
            if(prt) mf::LogVerbatim("TC")<<"+++++++ checking ClusterHitsOK with jht "<<fPlane<<":"<<PrintHit(tjs.fHits[jht])<<" BB Multiplicity "<<jHitsInMultiplet.size()<<" HitsRMSTick "<<HitsRMSTick(tjs, jHitsInMultiplet, kUnusedHits)<<" fatJhit "<<fatJHit<<" setting toTick to "<<(int)toTick;
            // Ensure that the hits StartTick and EndTick have the proper overlap
            if(!TrajHitsOK(iHitsInMultiplet, jHitsInMultiplet)) continue;
            // start a trajectory with direction from iht -> jht
            Trajectory work;
            if(!StartTraj(work, fromWire, fromTick, toWire, toTick, fCTP, pass)) continue;
            if(didPrt) TJPrt = work.WorkID;
            // check for a major failure
            if(fQuitAlg) return;
            if(work.Pts.empty()) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: StartTraj failed";
              prt = false;
              ReleaseHits(work);
              continue;
            }
            unsigned short angRange = AngleRange(work.Pts[0]);
            // check for a large angle crawl
            if(angRange > fMaxAngleRange[work.Pass]) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: Wrong angle range "<<angRange<<" for this pass "<<work.Pass;
              prt = false;
              ReleaseHits(work);
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
            if(!sigOK || NumUsedHits(work.Pts[0]) == 0) {
              if(prt) mf::LogVerbatim("TC")<<" No hits at initial trajectory point ";
              prt = false;
              ReleaseHits(work);
              continue;
            }
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
                prt = false;
                ReleaseHits(work);
                continue;
              } // Failed again
            }
            // Check the quality of the work trajectory
            CheckTraj(work);
            // check for a major failure
            if(fQuitAlg) return;
            if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: After CheckWork EndPt "<<work.EndPt[0]<<"-"<<work.EndPt[1]<<" fGoodTraj "<<fGoodTraj<<" fTryWithNextPass "<<fTryWithNextPass;
            if(fTryWithNextPass) {
              // Most likely, the first part of the trajectory was good but the latter part
              // had too many unused hits. The work vector was
              // truncated and pass incremented, so give it another try
              work.AlgMod[kTryWithNextPass] = true;
              StepCrawl(work);
              // check for a major failure
              if(fQuitAlg) return;
              if(!fGoodTraj) {
                if(prt) mf::LogVerbatim("TC")<<" xxxxxxx StepCrawl failed AGAIN after CheckWork";
                ReleaseHits(work);
                continue;
              } // Failed again
            } // fTryWithNextPass
            if(prt) mf::LogVerbatim("TC")<<"StepCrawl done: fGoodTraj "<<fGoodTraj<<" NumPtsWithCharge "<<NumPtsWithCharge(work, true)<<" cut "<<fMinPts[work.Pass];
            // decide if the trajectory is long enough for this pass
            if(!fGoodTraj || NumPtsWithCharge(work, true) < fMinPts[work.Pass]) {
              if(prt) mf::LogVerbatim("TC")<<" xxxxxxx Not enough points "<<NumPtsWithCharge(work, true)<<" minimum "<<fMinPts[work.Pass]<<" or !fGoodTraj";
              ReleaseHits(work);
              continue;
            }
            ReversePropagate(work);
            if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: calling StoreTraj with npts "<<work.EndPt[1];
            StoreTraj(work);
            // check for a major failure
            if(fQuitAlg) return;
            break;
          } // jht
        } // iht
      } // iwire
      EndMerge();
      if(fQuitAlg) return;
      TagDeltaRays(tjs, fCTP, fDeltaRayTag, debug.WorkID);
      Find2DVertices();
      if(fQuitAlg) return;
    } // pass
    
    // Use unused hits in all trajectories
    UseUnusedHits();
    
    // make junk trajectories using nearby un-assigned hits
    if(fJTMaxHitSep2 > 0) {
      FindJunkTraj();
      if(fQuitAlg) return;
    }
    TagDeltaRays(tjs, fCTP, fDeltaRayTag, debug.WorkID);
    TagShowerTraj(tjs, fCTP, fShowerTag, debug.WorkID);
    Find2DVertices();
     // check for a major failure
    if(fQuitAlg) return;
    
    // last attempt to attach Tjs to vertices
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) if(tjs.vtx[ivx].NTraj > 0) AttachAnyTrajToVertex(tjs, ivx, fVertex2DCuts, vtxPrt);
    
    // Refine vertices, trajectories and nearby hits
    Refine2DVertices();
    
  } // ReconstructAllTraj

  //////////////////////////////////////////
  void TrajClusterAlg::UseUnusedHits()
  {
    if(tjs.allTraj.size() == 0) return;
    if(!fUseAlg[kUseUnusedHits]) return;
    
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
        if(AngleRange(tj.Pts[ipt]) == 0) continue;
//        if(!IsLargeAngle(tj.Pts[ipt])) continue;
        bool hitsAdded = false;
        for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          // hit is associated with this point and it is not used
          if(tj.Pts[ipt].UseHit[ii]) continue;
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          // and it is not used in any other trajectory
          if(tjs.fHits[iht].InTraj != 0) continue;
          tj.Pts[ipt].UseHit[ii] = true;
          tjs.fHits[iht].InTraj = tj.ID;
          hitsAdded = true;
        } // ii
        if(hitsAdded) {
          DefineHitPos(tj.Pts[ipt]);
          tj.AlgMod[kUseUnusedHits] = true;
          if(prt) mf::LogVerbatim("TC")<<"UseUnusedHits: Using hits on ipt "<<ipt;
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
    
    if(!fUseAlg[kRevProp]) return;
    
    if(tj.Pts.size() < 6) return;
    
    // Only consider trajectories that have had their beginning trajectory points
    // updated by FixTrajBegin
    if(!tj.AlgMod[kFixEnd]) return;
    
    // decide how complicated this should get. If the delta values at the beginning are not
    // too bad, we can just mask them off. If this is not the case, then the trajectory should
    // be truly reverse propagated, which necessitates re-fitting the trajectory points at the beginning

    if(prt) {
      mf::LogVerbatim("TC")<<"ReversePropagate: ";
    }
    
    // deal with this case later
    if(tj.EndPt[0] != 0) return;
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
 //   std::cout<<"chk "<<ipl<<":"<<nextWire<<" "<<tjs.LastWire[ipl]<<" "<<theLastWire<<"\n";
    // clone the first point
    TrajPoint tp = tj.Pts[0];
    // strip off the hits
    tp.Hits.clear(); tp.UseHit.reset();
    // move it to the next wire
    MoveTPToWire(tp, (float)nextWire);
    // find close unused hits near this position
    float maxDelta = 10 * tj.Pts[tj.EndPt[1]].DeltaRMS;
    if(!FindCloseHits(tjs, tp, maxDelta, kUnusedHits)) return;
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
    tjWork.Pts[lastPt].NTPsFit = 3;
    // update the charge
    float chg = 0;
    float cnt = 0;
    for(unsigned short ii = 0; ii < 4; ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tjWork.Pts[ipt].Chg == 0) continue;
      chg += tjWork.Pts[ipt].Chg;
      ++cnt;
    } // ii
    if(cnt > 1) tjWork.Pts[lastPt].AveChg = chg / cnt;
    StepCrawl(tjWork);
    if(!fGoodTraj) {
      if(prt) mf::LogVerbatim("TC")<<" ReversePropagate StepCrawl failed";
      return;
    }
    // restore the original direction
    ReverseTraj(tjs, tjWork);
    tj = tjWork;
    if(prt) mf::LogVerbatim("TC")<<" ReversePropagate success";

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
        prt = (iht == debug.Hit);
        if(prt) {
          mf::LogVerbatim("TC")<<"FindJunkTraj: Found debug hit "<<PrintHit(tjs.fHits[iht])<<" InTraj "<<tjs.fHits[iht].InTraj<<" fJTMaxHitSep2 "<<fJTMaxHitSep2;
        }
        if(tjs.fHits[iht].InTraj != 0) continue;
        for(jht = jfirsthit; jht < jlasthit; ++jht) {
          if(tjs.fHits[jht].InTraj != 0) continue;
          if(prt && HitSep2(tjs, iht, jht) < 100) mf::LogVerbatim("TC")<<" use "<<PrintHit(tjs.fHits[jht])<<" HitSep2 "<<HitSep2(tjs, iht, jht);
          if(HitSep2(tjs, iht, jht) > fJTMaxHitSep2) continue;
          tHits.clear();
          // add all hits and flag them
          fromIndex = iht - tjs.fHits[iht].LocalIndex;
          for(kht = fromIndex; kht < fromIndex + tjs.fHits[iht].Multiplicity; ++kht) {
            if(tjs.fHits[kht].InTraj != 0) continue;
            tHits.push_back(kht);
            tjs.fHits[kht].InTraj = -4;
          } // kht
          fromIndex = jht - tjs.fHits[jht].LocalIndex;
          for(kht = fromIndex; kht < fromIndex + tjs.fHits[jht].Multiplicity; ++kht) {
            if(tjs.fHits[kht].InTraj != 0) continue;
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
                // this shouldn't be needed but do it anyway
                if(std::find(tHits.begin(), tHits.end(), kht) != tHits.end()) continue;
                // check w every hit in tHit
                for(tht = 0; tht < tHits.size(); ++tht) {
//                  if(prt && HitSep2(kht, tHits[tht]) < 100) mf::LogVerbatim("TC")<<" kht "<<PrintHit(tjs.fHits[kht])<<" tht "<<PrintHit(tjs.fHits[tHits[tht]])<<" HitSep2 "<<HitSep2(kht, tHits[tht])<<" cut "<<fJTMaxHitSep2;
                  if(HitSep2(tjs, kht, tHits[tht]) > fJTMaxHitSep2) continue;
                  hitsAdded = true;
                  tHits.push_back(kht);
                  tjs.fHits[kht].InTraj = -4;
                  if(tHits.size() > 50) {
                    mf::LogWarning("TC")<<"FindJunkTraj: More than 50 hits found in junk trajectory. Stop looking";
                    break;
                  }
                  if(kwire > hiWire) hiWire = kwire;
                  if(kwire < loWire) loWire = kwire;
                  break;
                } // tht
              } // kht
//              if(prt) mf::LogVerbatim("TC")<<" kwire "<<kwire<<" thits size "<<tHits.size();
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
          if(prt) {
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
    
    if(!fUseAlg[kJunkTj]) return;
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
      // Rotate the hits into this coordinate system to find the start and end
      // points and general direction
      float cs = cos(-work.Pts[0].Ang);
      float sn = sin(-work.Pts[0].Ang);
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
//    work.Pts[0].UseHit.resize(work.Pts[0].Hits.size(), true);
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
//      tp.UseHit.resize(tp.Hits.size(), true);
      DefineHitPos(tp);
      // Just use the hit position as the traj position
      tp.Pos = tp.HitPos;
      if(TrajPointSeparation(work.Pts[ipt-1], tp) < 0.5) continue;
      // define the direction
      MakeBareTrajPoint(tjs, work.Pts[ipt-1], tp, tpd);
      if(tpd.Pos[0] < 0) {
        // bad direction
//        mf::LogError("TC")<<"Bad direction in MakeJunkTraj using points on ipt "<<ipt;
        continue;
      }
      tp.Dir = tpd.Dir;
      tp.Ang = tpd.Ang;
      if(ipt == 1) {
        work.Pts[0].Dir = tpd.Dir;
        work.Pts[0].Ang = tpd.Ang;
      }
      work.Pts.push_back(tp);
      SetEndPoints(tjs, work);
    }
    if(prt) {
      PrintTrajectory("MJT", tjs, work, USHRT_MAX);
    }
    work.AlgMod[kJunkTj] = true;
    // See if this works or just does damage
    fGoodTraj = true;
//    CheckTraj(work);
//    if(!fGoodTraj) return;
    // Finally push it onto tjs.allTraj
    StoreTraj(work);
    if(fQuitAlg) return;
    // return with a valid index for the new trajectory
    newTjIndex = tjs.allTraj.size() - 1;
  } // MakeJunkTraj

  //////////////////////////////////////////
  void TrajClusterAlg::MatchTruth()
  {
    
    if(fIsRealData) return;
    if(fMatchTruth[0] < 0) return;
    
    art::ServiceHandle<cheat::BackTracker> bt;
    // list of all true particles
    sim::ParticleList const& plist = bt->ParticleList();
    if(plist.empty()) return;
    
    // MC Particles for the desired true particles
    int sourcePtclTrackID = -1;
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
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kBeamNeutrino;
          }
          if(theTruth->Origin() == simb::kSingleParticle) {
            sourcePtclTrackID = trackID;
            sourceOrigin = simb::kSingleParticle;
          }
          if(sourceOrigin == simb::kBeamNeutrino) {
            // histogram the vertex position difference
            for(auto& aVtx3 : tjs.vtx3) {
              fNuVtx_dx->Fill(part->Vx() - aVtx3.X);
              fNuVtx_dy->Fill(part->Vy() - aVtx3.Y);
              fNuVtx_dz->Fill(part->Vz() - aVtx3.Z);
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
      if(fMatchTruth[1] > 2) std::cout<<partList.size()<<" PDG Code  "<<part->PdgCode()<<" TrackId "<<part->TrackId()<<" sourcePtclTrackID  "<<sourcePtclTrackID<<" Origin "<<theTruth->Origin()<<" Process "<<part->Process()<<"\n";
      partList.push_back(part);
    } // ipart
    // Match all hits to the truth. Put the MC track ID in a temp vector
    std::vector<int> hitTruTrkID(tjs.fHits.size());
    // Prepare to count of the number of hits matched to each MC Track in each plane
    std::vector<std::vector<unsigned short>> nMatchedHitsInPartList(plist.size());
    for(unsigned short ipl = 0; ipl < plist.size(); ++ipl) nMatchedHitsInPartList[ipl].resize(tjs.NumPlanes);
    // and make a list of the TJs and hit count for each MC Track
    std::vector<std::vector<std::array<unsigned short, 2>>> nMatchedHitsInTj(partList.size());
    
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      TCHit& hit = tjs.fHits[iht];
      raw::ChannelID_t channel = geom->PlaneWireToChannel((int)hit.WireID.Plane, (int)hit.WireID.Wire, (int)hit.WireID.TPC, (int)hit.WireID.Cryostat);
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
    
    if(fMatchTruth[1] > 1) {
      mf::LogVerbatim myprt("TC");
      myprt<<"part   PDG TrkID MomID KE(MeV)   Process nMatchedHits[plane]\n";
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
        myprt<<std::setw(12)<<partList[ipl]->Process();
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) myprt<<std::setw(4)<<nMatchedHitsInPartList[ipl][plane];
        myprt<<"  tjID_nTruHits";
        for(auto& tmp : nMatchedHitsInTj[ipl]) myprt<<" "<<tmp[0]+1<<"_"<<tmp[1];
        myprt<<"\n";
      } // ipl
    }

    // Declare a TJ - partlist match for the trajectory which has the most true hits
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
        if(tjWithMostHits > tjs.allTraj.size() - 1) continue;
        // the total number of hits used in the TJ
        auto tmp = PutTrajHitsInVector(tjs.allTraj[tjWithMostHits], kUsedHits);
        float nTjHits = tmp.size();
        float nTruHits = nMatchedHitsInPartList[ipl][plane];
        float nTjTruRecHits = mostHits;
        float eff = nTjTruRecHits / nTruHits;
        float pur = nTjTruRecHits / nTjHits;
        float effpur = eff * pur;
        // This overwrites any previous match that has poorer efficiency * purity
        if(effpur > tjs.allTraj[tjWithMostHits].EffPur) {
          tjs.allTraj[tjWithMostHits].TruPDG = partList[ipl]->PdgCode();
          tjs.allTraj[tjWithMostHits].TruKE = 1000 * (partList[ipl]->E() - partList[ipl]->Mass());
          tjs.allTraj[tjWithMostHits].EffPur = effpur;
          partListToTjID[ipl][plane] = tjs.allTraj[tjWithMostHits].ID;
//          std::cout<<ipl<<" plane "<<plane<<" nTruHits "<<nTruHits<<" Tj ID "<<tjs.allTraj[tjWithMostHits].ID<<" nTjHits "<<nTjHits<<" nTjTruRecHits "<<nTjTruRecHits<<" eff "<<eff<<" pur "<<pur<<" effpur "<<eff*pur<<"\n";
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
        ++EPCounts[pdgIndex];
        // find the first and last matched hit in this plane
        unsigned int fht = UINT_MAX;
        unsigned int lht = 0;
        // find the first and last matched hit in this plane
        for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
          if(tjs.fHits[iht].WireID.Plane != plane) continue;
          if(hitTruTrkID[iht] != partList[ipl]->TrackId()) continue;
          if(fht == UINT_MAX) fht = iht;
          lht = iht;
        } // iht
        if(partListToTjID[ipl][plane] == 0) {
          // Enter 0 in the profile histogram
          fEP_T[pdgIndex]->Fill(TMeV, 0);
          if(nMatchedHitsInPartList[ipl][plane] > fMatchTruth[3]) {
            mf::LogVerbatim myprt("TC");
            myprt<<"pdgIndex "<<pdgIndex<<" No matched trajectory to partList["<<ipl<<"]";
            myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht]);
          }
          continue;
        }
        unsigned short itj = partListToTjID[ipl][plane] - 1;
        EPSums[pdgIndex] += tjs.allTraj[itj].EffPur;
        fEP_T[pdgIndex]->Fill(TMeV, tjs.allTraj[itj].EffPur);
        // print out some debugging information if the EP was pitiful and the number of matched hits is large
        if(tjs.allTraj[itj].EffPur < fMatchTruth[2] && nMatchedHitsInPartList[ipl][plane] > fMatchTruth[3]) {
          mf::LogVerbatim myprt("TC");
          myprt<<"pdgIndex "<<pdgIndex<<" BadEP "<<tjs.allTraj[itj].EffPur<<" nMatchedHitsInPartList "<<nMatchedHitsInPartList[ipl][plane];
          myprt<<" from true hit "<<PrintHit(tjs.fHits[fht])<<" to "<<PrintHit(tjs.fHits[lht]);
        }
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

  } // MatchTruth

  ////////////////////////////////////////////////
  void TrajClusterAlg::AddLAHits(Trajectory& tj, unsigned short ipt, bool& sigOK)
  {
    // Large Angle version of AddHits

    TrajPoint& tp = tj.Pts[ipt];
    tp.Hits.clear();
    tp.UseHit.reset();
    sigOK = false;

    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short ipl = planeID.Plane;

    // look at adjacent wires for larger angle trajectories
    // We will check the most likely wire first
    std::vector<unsigned int> wires(1);
    wires[0] = std::nearbyint(tp.Pos[0]);
    if(wires[0] > tjs.LastWire[ipl]-1) return;
    
    unsigned short angRange = AngleRange(tp);
    if(angRange == 0) {
      mf::LogVerbatim("TC")<<"AddLAHits called with AngleRange = 0. Don't do this";
      return;
    }
    // and the adjacent wires next in the most likely order only
    // after the first point has been defined and for the largest angle range
    if(ipt > 0 && angRange == fAngleRanges.size() - 1) {
      if(wires[0] > 0) wires.push_back(wires[0] - 1);
      if(wires[0] < tjs.LastWire[ipl]-1) wires.push_back(wires[0] + 1);
    } // ipt > 0 ...
    
    // Determine the time (tick) window expected for a TP with this angle
//    raw::TDCtick_t tickWindow = ExpectedHitsRMS(tp);
    // Modify this if necessary to account for the step size of 2 WSE units for LA trajectories
//    raw::TDCtick_t stepWindow = 2 / tjs.UnitsPerTick;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" AddLAHits: AngleRange "<<angRange<<" Wires under consideration";
      for(auto& wire : wires) myprt<<" "<<wire;
    }
    
    // a temporary tp that we can move around
    TrajPoint ltp = tp;
    // do this while testing
    sigOK = false;
    
    tp.Hits.clear();
    for(unsigned short ii = 0; ii < wires.size(); ++ii) {
      unsigned int wire = wires[ii];
      if(wire > tjs.LastWire[ipl]) continue;
      // Assume a signal exists on a dead wire
      if(tjs.WireHitRange[fPlane][wire].first == -1) sigOK = true;
      if(tjs.WireHitRange[fPlane][wire].first < 0) continue;
      std::array<unsigned int, 2> wireWindow {wire, wire};
      MoveTPToWire(ltp, wire);
      std::array<float, 2> timeWindow {ltp.Pos[1] - 3, ltp.Pos[1] + 3};
      bool hitsNear;
      // Look for hits using the requirement that the timeWindow overlaps with the hit StartTick and EndTick
      auto closeHits = FindCloseHits(tjs, wireWindow, timeWindow, ipl, kUnusedHits, false, hitsNear);
      if(hitsNear) sigOK = true;
      tp.Hits.insert(tp.Hits.end(), closeHits.begin(), closeHits.end());
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" LAPos "<<tp.CTP<<":"<<PrintPos(tjs, ltp)<<" Tick window "<<(int)(timeWindow[0]/tjs.UnitsPerTick)<<" to "<<(int)(timeWindow[1]/tjs.UnitsPerTick);
        for(auto& iht : closeHits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
      } // prt
    } // ii
/*
    std::vector<unsigned int> hitsInMultiplet;
    
    for(unsigned short ii = 0; ii < wires.size(); ++ii) {
      unsigned int wire = wires[ii];
      if(wire > tjs.LastWire[fPlane]) continue;
      // Assume a signal exists on a dead wire
      if(tjs.WireHitRange[fPlane][wire].first == -1) sigOK = true;
      if(tjs.WireHitRange[fPlane][wire].first < 0) continue;
      raw::TDCtick_t rawProjTick = tp.Pos[1];
      if(std::abs(tp.Dir[0]) > 1E-4) rawProjTick = (raw::TDCtick_t)(tp.Pos[1] + ((float)wire - tp.Pos[0]) * tp.Dir[1] / tp.Dir[0]);
      rawProjTick /= tjs.UnitsPerTick;
      
      raw::TDCtick_t minTick = 0;
      if(rawProjTick > tickWindow) minTick = rawProjTick - tickWindow;
      raw::TDCtick_t maxTick = rawProjTick + tickWindow;
      if(prt) mf::LogVerbatim("TC")<<"Wire "<<wire<<" minTick "<<minTick<<" rawProjTick "<<rawProjTick<<" maxTick "<<maxTick;
      unsigned int firstHit = (unsigned int)tjs.WireHitRange[fPlane][wire].first;
      unsigned int lastHit = (unsigned int)tjs.WireHitRange[fPlane][wire].second;
      for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
        // ensure that it isn't already in this tj
        if(std::find(tjHits.begin(), tjHits.end(), iht) != tjHits.end()) {
          sigOK = true;
          continue;
        }
        // ignore this hit if it was part of the multiplet found in the previous iteration
        if(hitsInMultiplet.size() > 1 && std::find(hitsInMultiplet.begin(), hitsInMultiplet.end(), iht) != hitsInMultiplet.end()) {
          sigOK = true;
          continue;
        }
        GetHitMultiplet(iht, hitsInMultiplet);
        auto closeHits = FindCloseHits();
        if(prt) {
          std::cout<<" mHits ";
          // only print the hits that will be used
          for(auto& jht : hitsInMultiplet) std::cout<<" "<<PrintHit(tjs.fHits[jht]);
          std::cout<<"\n";
        }
        // find the position and charge of the multiplet
        float mChg;
        float mTime = HitsPosTime(tjs, hitsInMultiplet, mChg, kUnusedHits);
        // find the rms of these hits in time
        float mRMS = HitsRMSTime(tjs, hitsInMultiplet, kUnusedHits);
        // estimate the position error transverse to the direction in WSE units
        if(prt) std::cout<<" mTime (ticks) "<<(int)mTime/tjs.UnitsPerTick<<" mChg "<<(int)mChg<<" mRMS (ticks) "<<(int)mRMS/tjs.UnitsPerTick<<"\n";
        float mTimeRMS = fHitErrFac * mRMS * std::abs(tp.Dir[0]);
        if(mTimeRMS < 0.1) mTimeRMS = 0.1;
        float delta = PointTrajDOCA(tjs, wire, mTime, tp);
        // we shouldn't use charge unless the hit finder handles long signal really well
        float mPull = delta / mTimeRMS;
        if(prt && abs(mTime/tjs.UnitsPerTick - rawProjTick) < 500) {
          mf::LogVerbatim myprt("TC");
          myprt<<" LAPos "<<tjs.fHits[iht].WireID.Plane<<":"<<wire<<":"<<rawProjTick<<"_"<<tjs.fHits[iht].InTraj;
          myprt<<" Mult "<<hitsInMultiplet.size();
          myprt<<" mTime (ticks) "<<(int)(mTime/tjs.UnitsPerTick);
          myprt<<" mChg "<<(int)mChg;
          myprt<<" mPos1RMS (ticks) "<<(int)(mTimeRMS/tjs.UnitsPerTick);
          myprt<<" delta "<<std::fixed<<std::setprecision(1)<<delta;
          myprt<<" mPull "<<std::setprecision(1)<<mPull;
          myprt<<" Signal? "<<sigOK;
          myprt<<" mHits ";
          // only print the hits that will be used
          for(auto& jht : hitsInMultiplet) myprt<<" "<<PrintHit(tjs.fHits[jht]);
//          for(auto& jht : hitsInMultiplet) if(tjs.fHits[jht].InTraj == 0) myprt<<" "<<PrintHit(tjs.fHits[jht])<<"_"<<tjs.fHits[jht].InTraj;
        }
        if(mPull > 5) continue;
        sigOK = true;
        // Associate all hits in the multiplet with the TP
        for(auto& jht : hitsInMultiplet) {
          // ensure that it isn't already in this tj
          if(std::find(tjHits.begin(), tjHits.end(), jht) != tjHits.end()) continue;
          tp.Hits.push_back(jht);
          tjHits.push_back(jht);
        } // jht
      } // iht
    } // ii
*/
    // no hits found
    if(tp.Hits.empty()) {
      if(prt) mf::LogVerbatim("TC")<<" AddLAHits: No hits found ";
     return;
    }
    
    if(tp.Hits.size() > 16) {
//      mf::LogWarning("TC")<<"AddLAHits: More than 16 hits added to TP. Truncating to the max for UseHit";
      tp.Hits.resize(16);
    }
    // Use all of the hits that are available
    tp.UseHit.reset();
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      unsigned int iht = tp.Hits[ii];
      if(tjs.fHits[iht].InTraj == 0) {
        tp.UseHit[ii] = true;
        tjs.fHits[iht].InTraj = tj.ID;
      }
    } // ii
    DefineHitPos(tp);
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" HitPos "<<std::fixed<<std::setprecision(1)<<tp.HitPos[0]<<" "<<(int)tp.HitPos[1]/tjs.UnitsPerTick;
      myprt<<" using hits";
      for(auto& iht : tp.Hits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
    }
    SetEndPoints(tjs, tj);
 
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
    
    // Call large angle hit finding only on the last 2 ranges
    unsigned short lastRange = fAngleRanges.size() - 1;
    if(lastRange > 1) --lastRange;
    if(AngleRange(tj.Pts[ipt]) > lastRange) {
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
    // loosen up a bit if we just passed a block of dead wires
    if(abs(dw) > 20 && DeadWireCount(tp.Pos[0], tj.Pts[lastPtWithUsedHits].Pos[0], tj.CTP) > 10) deltaCut *= 2;
    
    deltaCut *= fProjectionErrFactor;
    if(deltaCut > 3) deltaCut = 3;
    if(deltaCut < 1) deltaCut = 1;
    
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
      if(tjs.fHits[iht].InTraj == 0 && delta < bigDelta) {
        // An available hit that is just outside the window
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
      GetHitMultiplet(imBig, hitsInMultiplet, localIndex);
      for(auto jht : hitsInMultiplet) {
        if(std::find(closeHits.begin(), closeHits.end(), jht) != closeHits.end()) continue;
        if(tjs.fHits[jht].InTraj == tj.ID) continue;
        closeHits.push_back(jht);
      } // jht
      if(prt) mf::LogVerbatim("TC")<<" Added bigDelta hit "<<PrintHit(tjs.fHits[imBig])<<" w delta = "<<bigDelta;
    }
    if(!closeHits.empty()) sigOK = true;
    if(!sigOK) return;
    tp.Hits.insert(tp.Hits.end(), closeHits.begin(), closeHits.end());
    if(tp.Hits.size() > 16) {
//      mf::LogWarning("TC")<<"AddHits: More than 16 hits added to TP. Truncating to the max for UseHit";
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
    if(ipt < 4) useChg = false;
    // Don't use charge for reverse propagation since we are probably trying to attach large charge
    // hits at the beginning of a trajectory
    if(tj.AlgMod[kRevProp]) useChg = false;
    FindUseHits(tj, ipt, 10, useChg);
    DefineHitPos(tp);
    SetEndPoints(tjs, tj);
    if(prt) mf::LogVerbatim("TC")<<" number of close hits "<<closeHits.size()<<" used hits "<<NumUsedHits(tp);
  } // AddHits
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ReleaseHits(Trajectory& tj)
  {
    // Sets InTraj[] = 0 and UseHit false for all TPs in work. Called when abandoning work
//    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt)  UnsetUsedHits(tj.Pts[ipt]);
    for(auto& tp : tj.Pts) {
      for(auto& iht : tp.Hits) {
        if(tjs.fHits[iht].InTraj == tj.ID) tjs.fHits[iht].InTraj = 0;
      }
    } // tp
    
  } // ReleaseWorkHits

  //////////////////////////////////////////
  void TrajClusterAlg::UnsetUsedHits(TrajPoint& tp)
  {
    // Sets InTraj = 0 and UseHit false for all used hits in tp
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(tp.UseHit[ii]) {
        tjs.fHits[tp.Hits[ii]].InTraj = 0;
        tp.UseHit[ii] = false;
      } // UseHit
    } // ii
    tp.Chg = 0;
  } // UnsetUsedHits

  //////////////////////////////////////////
  void TrajClusterAlg::FindUseHits(Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg)
  {
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting InTraj < 0.
    
    if(ipt > tj.Pts.size() - 1) return;
    TrajPoint& tp = tj.Pts[ipt];
    
    if(tp.Hits.empty()) return;

    // Use everything (unused) for large angle TPs as long
    // as the multiplet doesn't include a hit used in another
    // trajectory

    if(prt) {
      mf::LogVerbatim("TC")<<"FUH:  maxDelta "<<maxDelta<<" useChg requested "<<useChg<<" TPHitsRMS "<<TPHitsRMSTick(tjs, tp, kUnusedHits)<<" AngleRange "<<AngleRange(tp);
    }
    float chgPullCut = 1000;
    if(useChg) chgPullCut = fChgPullCut;
    
    // large angle or maybe starting out a large angle trajectory
    geo::PlaneID iplID = DecodeCTP(tj.CTP);
    unsigned short ipl = iplID.Plane;
    bool fatHit = (TPHitsRMSTick(tjs, tp, kUnusedHits) > 4 * fAveHitRMS[ipl]);
    if(AngleRange(tp) > 0 || (fatHit && tj.Pts.size() < 4)) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(tjs.fHits[iht].InTraj > 0) continue;
        if(tjs.fHits[iht].InTraj == tj.ID) continue;
        tp.UseHit[ii] = true;
        tjs.fHits[iht].InTraj = tj.ID;
      } // ii
      if(prt) mf::LogVerbatim("TC")<<"FUH: isLA or short and fat. Using all hits ";
      return;
    } // IsLargeAngle

    // Find the hit that has the smallest delta
    tp.Delta = maxDelta;
    float delta;
    unsigned short imBest = USHRT_MAX;
    std::vector<float> deltas(tp.Hits.size(), 100);
    // keep track of the best delta - even if it is used
    float bestDelta = maxDelta;
    unsigned short nAvailable = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      tp.UseHit[ii] = false;
      unsigned int iht = tp.Hits[ii];
      delta = PointTrajDOCA(tjs, iht, tp);
      if(delta < bestDelta) bestDelta = delta;
      if(tjs.fHits[iht].InTraj > 0) continue;
      ++nAvailable;
      if(prt) {
        if(useChg) {
          if(prt) mf::LogVerbatim("TC")<<" "<<ii<<"  "<<PrintHit(tjs.fHits[iht])<<" delta "<<delta<<" Chg "<<(int)tjs.fHits[iht].Integral;
        } else {
          if(prt) mf::LogVerbatim("TC")<<" "<<ii<<"  "<<PrintHit(tjs.fHits[iht])<<" delta "<<delta;
        }
      } // prt
      deltas[ii] = delta;
      if(delta < tp.Delta) {
        tp.Delta = delta;
        imBest = ii;
      }
    } // ii
    
    if(prt) mf::LogVerbatim("TC")<<" imBest available "<<imBest<<" single hit. tp.Delta "<<tp.Delta<<" bestDelta "<<bestDelta;
    if(imBest == USHRT_MAX) return;
    unsigned int bestHit = tp.Hits[imBest];
    
    // don't use the best UNUSED hit if the best delta is for a USED hit and it is much better
    if(bestDelta < 0.5 * tp.Delta) return;
    
    if(!useChg || (useChg && (tj.AveChg == 0 || tj.ChgRMS == 0))) {
      // necessary quantities aren't available for more carefull checking
      if(prt) mf::LogVerbatim("TC")<<" tj.AveChg "<<tj.AveChg<<" or tj.ChgRMS "<<tj.ChgRMS<<" not defined yet. Use the best hit";
      tp.UseHit[imBest] = true;
      tjs.fHits[bestHit].InTraj = tj.ID;
      return;
    }
    
    // calculate the charge pull
    
    float bestHitChgPull = (tjs.fHits[bestHit].Integral / tj.AveChg - 1) / tj.ChgRMS;
    
    if(prt) mf::LogVerbatim("TC")<<" bestHit "<<PrintHit(tjs.fHits[bestHit])<<" Delta "<<tp.Delta<<" Charge "<<(int)tjs.fHits[bestHit].Integral<<" ChgPull "<<bestHitChgPull<<" nAvailable "<<nAvailable<<" tj.AveChg "<<tj.AveChg<<" tj.ChgRMS "<<tj.ChgRMS;
    
    // always use the best hit if the charge pull is OK
    if(bestHitChgPull > -chgPullCut && bestHitChgPull < chgPullCut) {
      tp.UseHit[imBest] = true;
      tjs.fHits[bestHit].InTraj = tj.ID;
    } // good charge
     else if(nAvailable == 1 && tj.PDGCode == 13 && tp.Delta < 2 * tp.DeltaRMS) {
       // special handling for muons. Allow higher or lower charge if the delta is very good
       if(bestHitChgPull > -2 * chgPullCut && bestHitChgPull < 2 * chgPullCut) {
         tp.UseHit[imBest] = true;
         tjs.fHits[bestHit].InTraj = tj.ID;
       }
    }
      
    // nothing fancy if there is only one hit available or if we are just starting out
    if(nAvailable == 1 || tj.Pts.size() < 4) {
      if(prt) mf::LogVerbatim("TC")<<" One hit available. Use it? "<<tp.UseHit[imBest];
      return;
    } // single hit
    
    // calculate some quantities for selecting the best hit
    float bestHitDeltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(bestHit);
    float arg = deltas[imBest] /  bestHitDeltaErr;
    float tpDRMS2 = tp.DeltaRMS * tp.DeltaRMS;
    float bestHitFOM = sqrt(arg * arg + tpDRMS2);
    // scale by charge pull if > 1
    if(useChg && bestHitChgPull > 1) bestHitFOM *= bestHitChgPull;
    
    // find the closest unused neighbor - the second best available hit
    unsigned short secondBest = USHRT_MAX;
    float secondBestDelta = 5;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(ii == imBest) continue;
      unsigned int iht = tp.Hits[ii];
      if(tjs.fHits[iht].InTraj > 0) continue;
      if(deltas[ii] < secondBestDelta) {
        secondBestDelta = deltas[ii];
        secondBest = ii;
      }
    } // ii
    if(secondBest == USHRT_MAX) return;
    
    // determine if the second best hit should be considered with the
    // first as a multiplet. Find the hit separation significance.
    unsigned int secondBestHit = tp.Hits[secondBest];
    float dtick = std::abs(tjs.fHits[bestHit].PeakTime - tjs.fHits[secondBestHit].PeakTime);
    float rms = tjs.fHits[bestHit].RMS;
    if(tjs.fHits[secondBestHit].RMS > rms) rms = tjs.fHits[secondBestHit].RMS;
    if(dtick / rms > fMultHitSep) {
      if(prt)  mf::LogVerbatim("TC")<<" secondBestHit separation too large. Use the single hit ";
      return;
    }
    
    // calculate some quantities for selecting the second best hit
    float secondBestHitChgPull = (tjs.fHits[secondBestHit].Integral / tj.AveChg - 1) / tj.ChgRMS;
    float secondBestHitDeltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(secondBestHit);
    arg = deltas[secondBest] / secondBestHitDeltaErr;
    float secondBestHitFOM = sqrt(arg * arg + tpDRMS2);
    // scale by charge pull if > 1
    if(useChg && secondBestHitChgPull > 1) secondBestHitFOM *= secondBestHitChgPull;
    
    // calculate doublet quantities
    float doubletChg = tjs.fHits[bestHit].Integral + tjs.fHits[secondBestHit].Integral;
    float doubletChgPull = (doubletChg / tj.AveChg - 1) / tj.ChgRMS;
    float doubletTime = (tjs.fHits[bestHit].Integral * tjs.fHits[bestHit].PeakTime + tjs.fHits[secondBestHit].Integral * tjs.fHits[secondBestHit].PeakTime) / doubletChg;
    doubletTime *= tjs.UnitsPerTick;
    // Square of delta here
    float doubletDelta2 = PointTrajDOCA2(tjs, tp.Pos[0], doubletTime, tp);
    float doubletTimeErr = (tjs.fHits[bestHit].Integral * tjs.fHits[bestHit].RMS + tjs.fHits[secondBestHit].Integral * tjs.fHits[secondBestHit].RMS) / doubletChg;
    doubletTimeErr *= tjs.UnitsPerTick;
    arg = doubletDelta2 / (doubletTimeErr * doubletTimeErr);
    // note that arg is already squared
    float doubletFOM = sqrt(arg + tpDRMS2);
    // scale by charge pull if > 1
    if(useChg && doubletChgPull > 1) doubletFOM *= doubletChgPull;
    
    if(prt) mf::LogVerbatim("TC")<<" bestHit FOM "<<bestHitFOM<<" secondBestHit "<<PrintHit(tjs.fHits[secondBestHit])<<" FOM "<<secondBestHitFOM<<" Doublet FOM "<<doubletFOM;
    
    // Only consider 1 or 2 of the available hits
    // don't use either hit if both are inconsistent with the average charge
    if(useChg && bestHitChgPull > chgPullCut && secondBestHitChgPull > chgPullCut) return;
    // Have 3 choices: the best hit, the second best hit or both combined (doublet)
    // calculate the combined charge pull and delta
    if(doubletChgPull > chgPullCut) {
      // doublet doesn't work. choose the best or second best
      if(bestHitFOM < secondBestHitFOM) {
        tp.UseHit[imBest] = true;
        tjs.fHits[bestHit].InTraj = tj.ID;
        return;
      } else {
        tp.UseHit[secondBest] = true;
        tjs.fHits[secondBestHit].InTraj = tj.ID;
        return;
      }
    } else {
      // consider the doublet. Assume the bestHit is the best
      float fom = bestHitFOM;
      unsigned short jj = imBest;
      unsigned int jhit = bestHit;
      if(secondBestHitFOM < fom) {
        // second best is better
        fom = secondBestHitFOM;
        jj = secondBest;
        jhit = secondBestHit;
      } // secondBestFOM < fom
      if(doubletFOM < fom) {
        // doublet is best
        tp.UseHit[imBest] = true;
        tjs.fHits[bestHit].InTraj = tj.ID;
        tp.UseHit[secondBest] = true;
        tjs.fHits[secondBestHit].InTraj = tj.ID;
        return;
      } // doubletFOM < fom
      // if we got here the jj and jhit variables make sense
      tp.UseHit[jj] = true;
      tjs.fHits[jhit].InTraj = tj.ID;
    } // consider the doublet
 
  } //  FindUseHits

  //////////////////////////////////////////
  void TrajClusterAlg::SetPoorUsedHits(Trajectory& tj, unsigned short ipt)
  {
    // Try to use the hits on this TP by reducing the number of points fitted. This
    // should only be done for reasonably long TJ's
    if(tj.Pts.size() < 2 * fMinPtsFit[tj.Pass]) return;
    
    // Restrict to TPs with only one hit for now
    if(tj.Pts[ipt].Hits.size() != 1) return;
    
    // Count the number of previously added TPs that were not used in the fit starting at ipt
    unsigned short jpt, nNotUsed = 0;
    for(jpt = ipt; jpt > 1; --jpt) {
      if(tj.Pts[jpt].Chg > 0) break;
      ++nNotUsed;
    }
    // Missed 0 or 1 TP. Don't mess with the fit
    if(nNotUsed < 1) return;
    
    // put the hit into the fit and see how it goes
    TrajPoint& tp = tj.Pts[ipt];
    unsigned int iht = tp.Hits[0];
    tp.UseHit[0] = true;
    tjs.fHits[iht].InTraj = tj.ID;
    unsigned short nptsFit = tp.NTPsFit;
    
    // Temp until we get this sorted out. Make sure that tj == work.
    if(tj.Pts.size() != tj.Pts.size() || ipt != tj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"SetPoorUsedHits: tj != work. Fix this code";
      fQuitAlg = true;
      return;
    }

    // increment the number of points fit ala UpdateWork
    ++tp.NTPsFit;
    FitTraj(tj);
    if(tp.FitChi < 2) {
      // Looks like this will work. Return victorious after
      // decrementing NTPsFit
      --tp.NTPsFit;
      return;
    }
    while(tp.FitChi > 2 && tp.NTPsFit > fMinPtsFit[tj.Pass]) {
      if(tp.NTPsFit > 15) {
        nptsFit = 0.7 * nptsFit;
      } else if(tp.NTPsFit > 4) {
        nptsFit -= 2;
      } else {
        nptsFit -= 1;
      }
      if(nptsFit < 3) nptsFit = 2;
      tp.NTPsFit = nptsFit;
      FitTraj(tj);
      if(prt) mf::LogVerbatim("TC")<<"  SetPoorUsedHits: FitChi "<<tp.FitChi<<" after reducing NTPsFit to "<<tp.NTPsFit;
      if(tp.NTPsFit <= fMinPtsFit[tj.Pass]) break;
    } // reducing number of fitted points
    
    if(tp.FitChi > 2) {
      // Well that didn't work. Undo UseHit and return
      tp.UseHit[0] = false;
      tjs.fHits[iht].InTraj = 0;
    }
    
  } // SetPoorUsedHits

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
      tp.HitPos[0] = tjs.fHits[iht].WireID.Wire;
      tp.HitPos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
      float wireErr = tp.Dir[1] * 0.289;
      float timeErr = tp.Dir[0] * HitTimeErr(iht);
      tp.HitPosErr2 = wireErr * wireErr + timeErr * timeErr;
      if(prt) mf::LogVerbatim("TC")<<"DefineHitPos: singlet "<<std::fixed<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]<<" HitPosErr "<<sqrt(tp.HitPosErr2);
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
    
    // Error is the wire error (1/sqrt(12))^2 if all hits are on one wire.
    // Scale it by the wire range
    float dWire = 1 + hiWire - loWire;
    float wireErr = tp.Dir[1] * dWire * 0.289;
    float timeErr2 = tp.Dir[0] * tp.Dir[0] * HitsTimeErr2(hitVec);
    tp.HitPosErr2 = wireErr * wireErr + timeErr2;
    if(prt) mf::LogVerbatim("TC")<<"DefineHitPos: multiplet "<<std::fixed<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]<<" HitPosErr "<<sqrt(tp.HitPosErr2);

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
  
  //////////////////////////////////////////
  void TrajClusterAlg::SplitTrajCrossingVertices()
  {
    // This is kind of self-explanatory...
    
    if(tjs.vtx.empty()) return;
    if(tjs.allTraj.empty()) return;
    
    vtxPrt = (debug.Plane == (int)fPlane && debug.Tick < 0);
    if(vtxPrt) mf::LogVerbatim("TC")<<"vtxPrt set for plane "<<fPlane<<" in SplitTrajCrossingVertices";
//    if(vtxPrt) PrintAllTraj("STCV", tjs, debug, USHRT_MAX, 999);
    
    // Splits trajectories in tjs.allTraj that cross a vertex
    unsigned short itj, iv, nTraj = tjs.allTraj.size();
    unsigned short tPass, closePt;
    float doca;
    for(itj = 0; itj < nTraj; ++itj) {
      // obsolete trajectory
      if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
      tPass = tjs.allTraj[itj].Pass;
      for(iv = 0; iv < tjs.vtx.size(); ++iv) {
        // obsolete vertex
        if(tjs.vtx[iv].NTraj == 0) continue;
        // not in the cryostat/tpc/plane
        if(tjs.allTraj[itj].CTP != tjs.vtx[iv].CTP) continue;
        // too short
        if(tjs.allTraj[itj].EndPt[1] < 6) continue;
        TrajClosestApproach(tjs.allTraj[itj], tjs.vtx[iv].Pos[0], tjs.vtx[iv].Pos[1], closePt, doca);
        if(vtxPrt)  mf::LogVerbatim("TC")<<" doca "<<doca<<" btw traj "<<itj<<" and tjs.vtx "<<iv<<" closePt "<<closePt<<" in plane "<<fPlane<<" CTP "<<tjs.vtx[iv].CTP;
        if(doca > fMaxVertexTrajSep[tPass]) continue;
        if(vtxPrt)  {
          mf::LogVerbatim("TC")<<"Good doca "<<doca<<" btw traj "<<itj<<" and tjs.vtx "<<iv<<" closePt "<<closePt<<" in plane "<<fPlane<<" CTP "<<tjs.vtx[iv].CTP;
          PrintTrajPoint("STCV", tjs, closePt, 1, tPass, tjs.allTraj[itj].Pts[closePt]);
        }
      } // iv
    } // itj
    
  } // SplitTrajCrossingVertices

  //////////////////////////////////////////
  void TrajClusterAlg::EndMerge()
  {
    // Merge Trajectorys in the order tj1 end = 1 to tj2 end = 0
    if(tjs.allTraj.size() < 2) return;
    if(!fUseAlg[kEndMerge]) return;

    unsigned short tj1, tj2, ipt;
    float chg1rms, chg2rms, chgpull;
    bool notJunk;

    mrgPrt = (debug.Plane == (int)fPlane && debug.Wire < 0);
    if(mrgPrt) mf::LogVerbatim("TC")<<"inside EndMerge on plane "<<fPlane;
    
    // Try several passes if there was a merge
    for(unsigned short nit = 0; nit < 4; ++nit) {
      bool didMerge = false;
      unsigned short tjsize = tjs.allTraj.size();
      for(tj1 = 0; tj1 < tjsize; ++tj1) {
        if(tjs.allTraj[tj1].AlgMod[kKilled]) continue;
        if(tjs.allTraj[tj1].CTP != fCTP) continue;
        // no merge if there is a vertex at the end
        if(tjs.allTraj[tj1].VtxID[1] > 0) continue;
        // no merge if it stops at end 1
        if(tjs.allTraj[tj1].StopsAtEnd[1]) continue;
        // ignore bad trajectories
        if(tjs.allTraj[tj1].MCSMom < 5) continue;
        ipt = tjs.allTraj[tj1].EndPt[1];
        TrajPoint& tp1 = tjs.allTraj[tj1].Pts[ipt];
        for(tj2 = 0; tj2 < tjsize; ++tj2) {
          if(tjs.allTraj[tj1].AlgMod[kKilled]) continue;
          if(tj1 == tj2) continue;
          if(tjs.allTraj[tj2].StepDir != tjs.allTraj[tj1].StepDir) continue;
          if(tjs.allTraj[tj2].AlgMod[kKilled]) continue;
          if(tjs.allTraj[tj2].CTP != fCTP) continue;
          if(tjs.allTraj[tj2].VtxID[0] > 0) continue;
          // no merge if it stops at end 1
          if(tjs.allTraj[tj1].StopsAtEnd[0]) continue;
          // ignore bad trajectories
          if(tjs.allTraj[tj2].MCSMom < 5) continue;
          ipt = tjs.allTraj[tj2].EndPt[0];
          TrajPoint& tp2 = tjs.allTraj[tj2].Pts[ipt];
          float doca = PointTrajDOCA(tjs, tp1.Pos[0], tp1.Pos[1], tp2);
          if(doca > 50) continue;
          float dang = DeltaAngle(tp1.Ang, tp2.Ang);
          float dwc = DeadWireCount(tp1, tp2);
          float ptSep0;
          if(tjs.allTraj[tj1].StepDir > 0) {
            ptSep0 = tp2.Pos[0] - tp1.Pos[0] - dwc;
          } else {
            ptSep0 = tp1.Pos[0] - tp2.Pos[0] - dwc;
          }
          bool twoMuons = tjs.allTraj[tj1].PDGCode == 13 && tjs.allTraj[tj2].PDGCode == 13;
          float maxWireSkip = fMaxWireSkipNoSignal;
          float tpSep2 = PosSep2(tp1.Pos, tp2.Pos);
          if(mrgPrt) {
            mf::LogVerbatim("TC")<<"Candidate tj1-tj2 "<<tjs.allTraj[tj1].ID<<"_1"<<"-"<<tjs.allTraj[tj2].ID<<"_0"<<" dang "<<dang<<" DeadWireCount "<<(int)dwc<<" ptSep0 "<<ptSep0<<" maxWireSkip "<<maxWireSkip<<" twoMuons? "<<twoMuons<<" doca "<<doca<<" tp1 "<<PrintPos(tjs, tp1)<<" tp2 "<<PrintPos(tjs, tp2)<<" tpSep2 "<<tpSep2;
          }
          // maximum overlap && maximum missed (live) wires
          if(ptSep0 < -2 || ptSep0 > maxWireSkip) continue;
          // Merge short trajectories, regardless of charge as long as the separation is small
          bool isClose = tpSep2 < 15 && dang < 1 && doca < 2;
          bool shortTraj = NumPtsWithCharge(tjs.allTraj[tj1], false) < 15 && NumPtsWithCharge(tjs.allTraj[tj2], false) < 15;
          bool OneIsDeltaRay = tjs.allTraj[tj1].PDGCode == 12 || tjs.allTraj[tj2].PDGCode == 12;
          if(!twoMuons && isClose && (shortTraj || OneIsDeltaRay)) {
            if(mrgPrt) mf::LogVerbatim("TC")<<" EM short tj1-tj2 "<<tjs.allTraj[tj1].ID<<"_1"<<"-"<<tjs.allTraj[tj2].ID<<"_0"<<" tp1 "<<PrintPos(tjs, tp1)<<" tp2 "<<PrintPos(tjs, tp2)<<" doca "<<doca;
            didMerge = MergeAndStore(tj1, tj2);
            if(fQuitAlg) return;
            if(!didMerge) continue;
            tjsize = tjs.allTraj.size();
            tjs.allTraj[tjsize-1].AlgMod[kEndMerge] = true;
            if(OneIsDeltaRay) tjs.allTraj[tjsize-1].PDGCode = 12;
            continue;
          }
          // handle muons here + uB problem with uncharacterized dead wires
          if(twoMuons && dang < 0.4 && doca < 10) {
            if(mrgPrt) mf::LogVerbatim("TC")<<" merge muon tj1-tj2 "<<tjs.allTraj[tj1].ID<<"_1"<<"-"<<tjs.allTraj[tj2].ID<<"_0"<<" tp1 "<<PrintPos(tjs, tp1)<<" tp2 "<<PrintPos(tjs, tp2)<<" doca "<<doca;
            didMerge = MergeAndStore(tj1, tj2);
            if(fQuitAlg) return;
            if(!didMerge) continue;
            tjsize = tjs.allTraj.size();
            tjs.allTraj[tjsize-1].AlgMod[kEndMerge] = true;
            continue;
          }
          if(dang > 0.4) continue;
          // Inflate the doca cut if we are bridging a block of dead wires
          float docaCut = 1;
          if(dwc > 10) docaCut *= 2;
          //        doca = PointTrajDOCA(tjs, tp1.Pos[0], tp1.Pos[1], tp2);
          if(mrgPrt) mf::LogVerbatim("TC")<<" EM "<<tjs.allTraj[tj1].ID<<"_1"<<"-"<<tjs.allTraj[tj2].ID<<"_0 DOCA "<<doca<<" docaCut "<<docaCut<<" dang "<<dang;
          if(doca > docaCut) continue;
          notJunk = (!tjs.allTraj[tj1].AlgMod[kJunkTj] && !tjs.allTraj[tj1].AlgMod[kJunkTj]);
          if(doca > 10) continue;
          // check the charge?
          chgpull = 0;
          if(notJunk && tjs.allTraj[tj1].AveChg > 0 && tjs.allTraj[tj2].AveChg > 0) {
            // calculate charge pull using the trajectory with the lowest charge RMS
            chg1rms = tjs.allTraj[tj1].ChgRMS * tjs.allTraj[tj1].AveChg;
            chg2rms = tjs.allTraj[tj2].ChgRMS * tjs.allTraj[tj2].AveChg;
            if(chg2rms < chg1rms) chg1rms = chg2rms;
            // What is this for?
//            if(chg1rms < 1) chg1rms = 0.15 * (tjs.allTraj[tj1].AveChg + tjs.allTraj[tj2].AveChg);
            chgpull = std::abs(tjs.allTraj[tj1].AveChg - tjs.allTraj[tj2].AveChg) / chg1rms;
            if(mrgPrt) mf::LogVerbatim("TC")<<"   chk tp1 chg "<<tp1.AveChg<<" tj1 chg "<<tjs.allTraj[tj1].AveChg<<" rms "<<tjs.allTraj[tj1].ChgRMS<<" tp2 chg "<<tp2.AveChg<<" tj2 chg "<<tjs.allTraj[tj2].AveChg<<" rms "<<tjs.allTraj[tj2].ChgRMS;
          } // charge is known.
          if(mrgPrt) mf::LogVerbatim("TC")<<" chgpull "<<chgpull<<" require < 4";
          if(chgpull > 4) continue;
          // time to merge them
          didMerge = MergeAndStore(tj1, tj2);
          if(fQuitAlg) return;
          if(!didMerge) continue;
          tjsize = tjs.allTraj.size();
          tjs.allTraj[tjsize-1].AlgMod[kEndMerge] = true;
          // flag the killed trajectories as well to aid in debugging
          tjs.allTraj[tj1].AlgMod[kEndMerge] = true;
          tjs.allTraj[tj2].AlgMod[kEndMerge] = true;
          if(mrgPrt) mf::LogVerbatim("TC")<<" merged "<<tjs.allTraj[tj1].ID<<"_1"<<"-"<<tjs.allTraj[tj2].ID<<"_0"<<" dang "<<dang<<" ptSep0 "<<ptSep0<<" doca "<<doca<<" tp1 "<<PrintPos(tjs, tp1)<<" tp2 "<<PrintPos(tjs, tp2)<<" . tjsize "<<tjsize;
        } // tj2
      } // tj1
      if(!didMerge) break;
    } // nit

  }  // EndMerge

  //////////////////////////////////////////
  void TrajClusterAlg::Refine2DVertices()
  {
    // Improve trajectories near vertices in the current plane
    if(tjs.vtx.empty()) return;
    
    if(!fUseAlg[kRefineVtx]) return;
    
    geo::PlaneID planeID = DecodeCTP(fCTP);
    unsigned short ipl = planeID.Plane;
    if(ipl != 0) return;
    
    // store a copy of everything so that we can recover gracefully if there is a major failure
    auto lHits = tjs.fHits;
    auto lWireHitRange = tjs.WireHitRange;
    auto lAllTraj = tjs.allTraj;
    bool majorFailure = false;
    
    if(vtxPrt) PrintHeader("R2D");
    
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      VtxStore& rvx = tjs.vtx[ivx];
      if(rvx.CTP != fCTP) continue;
      if(rvx.NTraj < 2) continue;
      // ensure that it is within the active volume of the TPC
      if(rvx.Pos[0] < 0 || rvx.Pos[0] > tjs.MaxPos0[ipl]) continue;
      if(rvx.Pos[1] < 0 || rvx.Pos[1] > tjs.MaxPos1[ipl]) continue;
      // make a list of TJs attached at each end and find the Region Of Confusion
      // wire and time ranges
      std::array<float, 2> wROC = {rvx.Pos[0], rvx.Pos[0]};
      std::array<float, 2> tROC = {rvx.Pos[1], rvx.Pos[1]};
      std::array<std::vector<unsigned short>, 2> tjlist;
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
        if(tjs.allTraj[itj].CTP != rvx.CTP) continue;
        Trajectory& tj = tjs.allTraj[itj];
        // ensure that the ID is OK so the code below doesn't choke
        if(tj.ID != itj + 1) {
          std::cout<<"Refine2DVertices allTraj ID "<<tj.ID<<" != itj "<<itj<<" + 1\n";
          fQuitAlg = true;
          return;
        }
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == rvx.ID) {
            tjlist[end].push_back(itj);
            unsigned short endPt = tj.EndPt[end];
            if(vtxPrt) PrintTrajectory("R2D", tjs, tj, endPt);
            // Find the lo/hi wire/time
            float arg = tj.Pts[endPt].Pos[0];
            if(arg < wROC[0]) wROC[0] = arg;
            if(arg > wROC[1]) wROC[1] = arg;
            arg = tj.Pts[endPt].Pos[1];
            if(arg < tROC[0]) tROC[0] = arg;
            if(arg > tROC[1]) tROC[1] = arg;
          }
        } // end
      } // itj
      // round to the nearest integer WSE unit
      wROC[0] = std::floor(wROC[0]);
      wROC[1] = std::ceil(wROC[1]);
      tROC[0] = std::floor(tROC[0]);
      tROC[1] = std::ceil(tROC[1]);
      std::cout<<"vtx "<<rvx.ID<<" tjlist[0] ";
      for(auto itj : tjlist[0]) std::cout<<" "<<itj+1;
      std::cout<<" tjlist[1] ";
      for(auto itj : tjlist[1]) std::cout<<" "<<itj+1;
      std::cout<<"\n";
      std::cout<<"wROC "<<wROC[0]<<" "<<wROC[1]<<" tROC "<<tROC[0]/tjs.UnitsPerTick<<" "<<tROC[1]/tjs.UnitsPerTick<<"\n";
      // no sense continuing unless there are 2 or more Tjs at at least one end
      if(tjlist[0].size() < 2 && tjlist[1].size() < 2) continue;
      // create a list of temporary hits in this region
      // Note that the ROC includes loWire AND hiWire
      unsigned int loWire = std::nearbyint(wROC[0]);
      unsigned int hiWire = std::nearbyint(wROC[1]);
      unsigned short ROCsize = hiWire - loWire + 1;
      // the wire that the vertex is on
      std::vector<TCHit> wireHits;
      std::cout<<"ROCsize "<<ROCsize<<"\n";

      // create hits on the ROC boundary for testing
      TCHit boxHit;
      boxHit.Integral = 100;
      boxHit.RMS = 1;
      boxHit.PeakAmplitude = 5;
      boxHit.InTraj = 0;
      for(unsigned int wire = loWire; wire <= hiWire; ++wire) {
        for(unsigned short tb = 0; tb < 2; ++tb) {
          DefineHit(boxHit, rvx.CTP, wire);
          boxHit.PeakTime = tROC[tb] / tjs.UnitsPerTick;
          CreateHit(boxHit);
        } // tb
      } // wire

      // Make a vector of ALL fHits that are inside the ROC so that we can erase them later
      std::array<unsigned int, 2> iwROC {loWire, hiWire};
      bool hitsNear;
      std::vector<unsigned int> fHitsInROC = FindCloseHits(tjs, iwROC, tROC, ipl, kAllHits, true, hitsNear);
      // sort by decreasing index so that hits that are later in fhits will be erased
      // before the earlier hits, obviating the need to correct fHitsInROC
      std::sort(fHitsInROC.begin(), fHitsInROC.end(), std::greater<unsigned int>());
      std::cout<<"fHitsInROC";
      for(auto& iht : fHitsInROC) std::cout<<" "<<iht<<"_"<<PrintHit(tjs.fHits[iht]);
      std::cout<<"\n";
      
      // Look for a trajectory that has a hit in the ROC but is not in tjlist
      bool skipVtx = false;
      for(auto& iht : fHitsInROC) {
        unsigned short inTj = tjs.fHits[iht].InTraj;
        if(inTj == 0) continue;
        unsigned short itj = inTj - 1;
        if(std::find(tjlist[0].begin(), tjlist[0].end(), itj) == tjlist[0].end() &&
           std::find(tjlist[1].begin(), tjlist[1].end(), itj) == tjlist[1].end()) {
          std::cout<<"Traj ID "<<inTj<<" not found in tjlist . Kill or keep?\n";
          std::array<float, 2> pos0 = tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[0]].Pos;
          std::array<float, 2> pos1 = tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[1]].Pos;
          if(pos0[0] > wROC[0] && pos0[0] < wROC[1] && pos1[0] > wROC[0] && pos1[0] < wROC[1] &&
             pos0[1] > tROC[0] && pos0[1] < tROC[1] && pos1[1] > tROC[0] && pos1[1] < tROC[1]) {
            // completely contained - kill it
            std::cout<<"Traj ID "<<inTj<<" completely contained in the ROC. Killing it\n";
            MakeTrajectoryObsolete(tjs, itj);
          } else {
            std::cout<<"Traj ID "<<inTj<<" is in the ROC but isn't attached to the vertex. Skip this vertex \n";
            skipVtx = true;
            break;
          }
        } // find
        if(skipVtx) break;
      } // iht
      if(skipVtx) break;

      // matching vectors of points outside the boundary of the ROC
      std::array<std::vector<unsigned short>, 2> edgePts;
      for(unsigned short end = 0; end < 2; ++end) {
        edgePts[end].resize(tjlist[end].size());

        // We now have a number of trajectories in VtxTraj that enter the ROC. The hits in fHits are still assigned to the
        // original trajectories in allTraj. Now create a set of vtx TCHits associated with VtxTraj within the ROC
        for(unsigned short iitj = 0; iitj < tjlist[end].size(); ++iitj) {
          unsigned short itj = tjlist[end][iitj];
          Trajectory& vtj = tjs.allTraj[itj];
          // reverse the trajectory to make changes easier
          if(end == 0)  ReverseTraj(tjs, vtj);
          if(vtj.ID == 1) PrintTrajectory("chk1", tjs, vtj, USHRT_MAX);
          // find the TP that is just outside the ROC. First assume that the end is inside.
          unsigned short edgePt = vtj.EndPt[1];
          // loWire   vtx      hiWire
          //     |     V          |
          // tj  |       E--------|-     end = 0, StepDir =  1 OR end = 1, StepDir = -1 (typical)
          // tj  |  E-------------|----  end = 0, StepDir =  1 OR end = 1, StepDir = -1 (not typical but happens)
          // tj  |       E------- |      end = 0, StepDir =  1 OR end = 1, StepDir = -1 (short tj inside the ROC)
          // tj -|---E            |      end = 0, StepDir = -1 OR end = 1, StepDir =  1
          for(unsigned short ii = 0; ii < ROCsize; ++ii) {
            edgePt = vtj.EndPt[1] - 1 - ii;
            if(edgePt == 0) break;
            unsigned int tWire = std::nearbyint(vtj.Pts[edgePt].Pos[0]);
            // keep going if there is a hit on this tp that is in fHitsInROC
            bool hitInROC = false;
            for(auto& iht : vtj.Pts[edgePt].Hits) {
              if(std::find(fHitsInROC.begin(), fHitsInROC.end(), iht) != fHitsInROC.end()) {
                hitInROC = true;
                break;
              }
            } // iht
            if(hitInROC) continue;
            // hit the wire boundary
            if(tWire < loWire || tWire > hiWire) break;
            // hit the time boundary
            if(vtj.Pts[edgePt].Pos[1] < tROC[0] || vtj.Pts[edgePt].Pos[1] > tROC[1]) break;
            // don't allow the trajectory to have < 2 points
            if(edgePt == 2) break;
          } // ii
          
          if(edgePt < 2) {
            std::cout<<"Not enough points left on vtxTraj "<<vtj.ID<<"\n";
            majorFailure = true;
            break;
          }
          
          edgePts[end][itj] = edgePt;
          // make a local TP that we can move around
          TrajPoint ltp = vtj.Pts[edgePt];
          
          std::cout<<"end "<<end<<" vtj.ID "<<vtj.ID<<" edgePt "<<edgePt<<" pos "<<PrintPos(tjs, vtj.Pts[edgePt])<<"\n";
          // find the first used hit in the tp and use it to characterize the
          // Charge and RMS of VtxHits inside the ROC
          float chg = vtj.Pts[edgePt].AveChg;
          float rms = TPHitsRMSTick(tjs, vtj.Pts[edgePt], kUsedHits);
          float amp = chg / (2.5066 * rms);
          // Modify the existing hits inside the ROC.
          // Form a list of hits that should be erased when we are done
          std::vector<unsigned int> killMe;
          for(unsigned short ipt = edgePt + 1; ipt < vtj.Pts.size(); ++ipt) {
            MoveTPToWire(ltp, vtj.Pts[ipt].Pos[0]);
            unsigned int nused = 0;
            for(unsigned short ii = 0; ii < vtj.Pts[ipt].Hits.size(); ++ii) {
              if(!vtj.Pts[ipt].UseHit[ii]) continue;
              unsigned int iht = vtj.Pts[ipt].Hits[ii];
              std::cout<<" tweak hit "<<PrintHit(tjs.fHits[iht]);
              ++nused;
              if(nused == 1) {
                tjs.fHits[iht].PeakTime = ltp.Pos[1] / tjs.UnitsPerTick;
                tjs.fHits[iht].PeakAmplitude = amp;
                tjs.fHits[iht].Integral = chg;
                tjs.fHits[iht].RMS = rms;
                std::cout<<" to "<<PrintHit(tjs.fHits[iht])<<"\n";
              } else {
                std::cout<<" erase this hit\n";
                killMe.push_back(iht);
              }
            } // ii
          } // ipt
          // erase hits?
          if(!killMe.empty()) {
            if(killMe.size() > 1) std::sort(killMe.begin(), killMe.end(), std::greater<unsigned int>());
            for(auto& iht : killMe) {
              tjs.fHits[iht].InTraj = 0;
              EraseHit(iht);
            }
          } // killMe not empty
        } // itj
      } // end
      if(majorFailure) break;
    } // ivx

    if(majorFailure) {
      // recover after a major failure
      tjs.fHits = lHits;
      tjs.WireHitRange = lWireHitRange;
      tjs.allTraj = lAllTraj;
    }

  } // Refine2DVertices

  //////////////////////////////////////////
  void TrajClusterAlg::FindVtxTraj(unsigned short ivx)
  {
    // Look for available hits in the vicinity of this vertex and try to make
    // a vertex trajectory from them
    
    if(!fUseAlg[kVtxTj]) return;
    
    if(ivx > tjs.vtx.size() - 1) return;
    if(tjs.vtx[ivx].Stat[kVtxTrjTried]) return;
    VtxStore& theVtx = tjs.vtx[ivx];
    
    std::array<unsigned int, 2> wireWindow;
    std::array<float, 2> timeWindow;
    
    // on the first try we look for small angle trajectories which will have hits
    // with a large wire window and a small time window
    // fVertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    wireWindow[0] = std::nearbyint(theVtx.Pos[0] - fVertex2DCuts[2]);
    wireWindow[1] = std::nearbyint(theVtx.Pos[0] + fVertex2DCuts[2]);
    timeWindow[0] = theVtx.Pos[1] - 5;
    timeWindow[1] = theVtx.Pos[1] + 5;
    
    geo::PlaneID planeID = DecodeCTP(theVtx.CTP);
    unsigned short ipl = planeID.Plane;

    if(vtxPrt) mf::LogVerbatim("TC")<<"inside FindVtxTraj "<<theVtx.ID<<" Window "<<wireWindow[0]<<" "<<wireWindow[1]<<" "<<timeWindow[0]<<" "<<timeWindow[1]<<" in plane "<<ipl;

    // find nearby available hits
    bool hitsNear;
    std::vector<unsigned int> closeHits = FindCloseHits(tjs, wireWindow, timeWindow, ipl, kUnusedHits, true, hitsNear);
    if(closeHits.empty()) return;
    if(vtxPrt) {
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
    if(vtxPrt) PrintHeader("FVT");
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
      if(vtxPrt) PrintTrajPoint("FVT", tjs, 0, tj.StepDir, tj.Pass, tp);
      // Step away and see what happens
      prt = vtxPrt;
      StepCrawl(tj);
      // check for a major failure
      if(fQuitAlg) return;
      // Check the quality of the trajectory
      CheckTraj(tj);
      if(!fGoodTraj || NumPtsWithCharge(tj, true) < fMinPts[tj.Pass]) {
        if(prt) mf::LogVerbatim("TC")<<" xxxxxxx Not enough points "<<NumPtsWithCharge(tj, true)<<" minimum "<<fMinPts[tj.Pass]<<" or !fGoodTraj";
        ReleaseHits(tj);
        continue;
      }
      if(vtxPrt) prt = false;
      tj.AlgMod[kVtxTj] = true;
      if(vtxPrt) mf::LogVerbatim("TC")<<"FindVtxTraj: calling StoreTraj with npts "<<tj.EndPt[1];
      StoreTraj(tj);
    } // ii
    
    // Flag this as tried so we don't try again
    tjs.vtx[ivx].Stat[kVtxTrjTried] = true;
  } // FindVtxTraj
  
  //////////////////////////////////////////
  void TrajClusterAlg::Find2DVertices()
  {
    // Self-explanatory
    // fVertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    
    if(fVertex2DCuts[0] <= 0) return;

    if(tjs.allTraj.size() < 2) return;
    
    vtxPrt = (debug.Plane == (int)fPlane && debug.Tick < 0);
    if(vtxPrt) {
      mf::LogVerbatim("TC")<<"vtxPrt set for plane "<<fPlane<<" in Find2DVertices";
      PrintAllTraj("F2DVi", tjs, debug, USHRT_MAX, tjs.allTraj.size());
    }
    
    unsigned short maxShortTjLen = fVertex2DCuts[0];

    for(unsigned short tj1 = 0; tj1 < tjs.allTraj.size() - 1; ++tj1) {
      if(tjs.allTraj[tj1].AlgMod[kKilled]) continue;
      if(tjs.allTraj[tj1].CTP != fCTP) continue;
      if(tjs.allTraj[tj1].PDGCode == 12) continue;
      bool tj1Short = (tjs.allTraj[tj1].EndPt[1] - tjs.allTraj[tj1].EndPt[0] < maxShortTjLen);
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tjs.allTraj[tj1].VtxID[end1] > 0) continue;
        unsigned short  endPt1 = tjs.allTraj[tj1].EndPt[end1];
        unsigned short oendPt1 = tjs.allTraj[tj1].EndPt[1-end1];
        TrajPoint& tp1 = tjs.allTraj[tj1].Pts[endPt1];
        for(unsigned short tj2 = tj1 + 1; tj2 < tjs.allTraj.size(); ++tj2) {
          if(tjs.allTraj[tj2].AlgMod[kKilled]) continue;
          if(tjs.allTraj[tj2].CTP != fCTP) continue;
          if(tjs.allTraj[tj2].PDGCode == 12) continue;
          bool tj2Short = (tjs.allTraj[tj2].EndPt[1] - tjs.allTraj[tj2].EndPt[0] < maxShortTjLen);
          for(unsigned short end2 = 0; end2 < 2; ++end2) {
            unsigned short  endPt2 = tjs.allTraj[tj2].EndPt[end2];
            unsigned short oendPt2 = tjs.allTraj[tj2].EndPt[1-end2];
            if(tjs.allTraj[tj1].VtxID[end1] > 0) continue;
            if(tjs.allTraj[tj2].VtxID[end2] > 0) continue;
            TrajPoint& tp2 = tjs.allTraj[tj2].Pts[endPt2];
            // Rough first cut on the separation between the end points of the
            // two trajectories
            float sepCut = 20;
            if(std::abs(tp1.Pos[0] - tp2.Pos[0]) > sepCut) continue;
            if(std::abs(tp1.Pos[1] - tp2.Pos[1]) > sepCut) continue;
            float wint, tint;
            TrajIntersection(tp1, tp2, wint, tint);
            // make sure this is inside the TPC
            if(wint < 0 || wint > tjs.MaxPos0[fPlane]) continue;
            if(tint < 0 || tint > tjs.MaxPos1[fPlane]) continue;
            // Next cut on separation between the TPs and the intersection point
            if(tj1Short || tj2Short) { sepCut = fVertex2DCuts[1]; } else { sepCut = fVertex2DCuts[2]; }
            float dw1 = wint - tp1.Pos[0];
            if(std::abs(dw1) > sepCut) continue;
            float dt1 = tint - tp1.Pos[1];
            if(std::abs(dt1) > sepCut) continue;
            float dw2 = wint - tp2.Pos[0];
            if(std::abs(dw2) > sepCut) continue;
            float dt2 = tint - tp2.Pos[1];
            if(std::abs(dt2) > sepCut) continue;
            if(vtxPrt) mf::LogVerbatim("TC")<<"F2DV candidate tj1-tj2 "<<tjs.allTraj[tj1].ID<<"_"<<end1<<"-"<<tjs.allTraj[tj2].ID<<"_"<<end2<<" tjs.vtx pos "<<(int)wint<<":"<<(int)(tint/tjs.UnitsPerTick)<<" tp1 "<<PrintPos(tjs, tp1)<<" tp2 "<<PrintPos(tjs, tp2);
            // make sure that the other end isn't closer
            if(PointTrajDOCA2(tjs, wint, tint, tp1) > PointTrajDOCA2(tjs, wint, tint, tjs.allTraj[tj1].Pts[oendPt1])) continue;
            if(PointTrajDOCA2(tjs, wint, tint, tp2) > PointTrajDOCA2(tjs, wint, tint, tjs.allTraj[tj2].Pts[oendPt2])) continue;
            float dang = DeltaAngle(tp1.Ang, tp2.Ang);
            if(vtxPrt) mf::LogVerbatim("TC")<<"  dang "<<dang;
            if(dang < 0.15) continue;
            if(!SignalAtPos(wint, tint, fCTP)) {
              // try nudging it by +/- 1 wire
              wint += 1;
              if(!SignalAtPos(wint, tint, fCTP)) {
                wint -= 2;
                if(!SignalAtPos(wint, tint, fCTP)) continue;
              } // nudge +1
            } // SignalAtPos
            // Ensure that the vertex position is close to an end
            unsigned short closePt1;
            float doca;
            TrajClosestApproach(tjs.allTraj[tj1], wint, tint, closePt1, doca);
            short dpt = abs((short)endPt1 - (short)closePt1);
            if(tjs.allTraj[tj1].EndPt[1] > 4) {
              if(dpt > 3) continue;
            } else {
              // tighter cut for short trajectories
              if(dpt > 0) continue;
            }
            unsigned short closePt2;
            TrajClosestApproach(tjs.allTraj[tj2], wint, tint, closePt2, doca);
            dpt = abs((short)endPt2 - (short)closePt2);
            if(tjs.allTraj[tj2].EndPt[1] > 4) {
              if(dpt > 3) continue;
            } else {
              // tighter cut for short trajectories
              if(dpt > 0) continue;
            }
            if(vtxPrt) mf::LogVerbatim("TC")<<" wint:tint "<<(int)wint<<":"<<(int)tint<<" ticks "<<(int)(tint/tjs.UnitsPerTick)<<" dang "<<dang;
            // make a new temporary vertex
            VtxStore aVtx;
            aVtx.Pos[0] = wint;
            aVtx.Pos[1] = tint;
            aVtx.NTraj = 0;
            aVtx.Pass = tjs.allTraj[tj1].Pass;
            aVtx.Topo = end1 + end2;
            aVtx.ChiDOF = 0;
            aVtx.CTP = fCTP;
//            aVtx.VtxMod[kFixed] = false;
            // try to fit it. We need to give it an ID to do that. Take the next
            // available ID
            unsigned short newVtxID = tjs.vtx.size() + 1;
            aVtx.ID = newVtxID;
            tjs.allTraj[tj1].VtxID[end1] = newVtxID;
            tjs.allTraj[tj2].VtxID[end2] = newVtxID;
            if(!FitVertex(tjs, aVtx, fVertex2DCuts, vtxPrt)) {
              tjs.allTraj[tj1].VtxID[end1] = 0;
              tjs.allTraj[tj2].VtxID[end2] = 0;
              continue;
            }
            // Save it
            tjs.vtx.push_back(aVtx);
            // Try to attach other tjs to it
            unsigned short ivx = tjs.vtx.size() - 1;
            if(vtxPrt) {
              mf::LogVerbatim myprt("TC");
              myprt<<" New vtx "<<aVtx.ID;
              myprt<<" Tjs "<<tjs.allTraj[tj1].ID<<"_"<<end1<<"-"<<tjs.allTraj[tj2].ID<<"_"<<end2;
              myprt<<" at "<<std::fixed<<std::setprecision(1)<<aVtx.Pos[0]<<":"<<aVtx.Pos[1]/tjs.UnitsPerTick<<" call AttachAnyTrajToVertex ";
            }
            AttachAnyTrajToVertex(tjs, ivx, fVertex2DCuts, vtxPrt);
            // try stepping away from this vertex
            FindVtxTraj(ivx);
          } // end2
        } // tj2
      } // end1
    } // tj1

    // Split trajectories that cross a vertex
    SplitTrajCrossingVertices();
    FindHammerVertices();
    FindHammerVertices2();

    if(vtxPrt) PrintAllTraj("F2DVo", tjs, debug, USHRT_MAX, USHRT_MAX);

  } // Find2DVertices
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindHammerVertices2()
  {
    // Variant of FindHammerVertices with slightly different requirements:
    // 1) tj1 is a straight trajectory with most of the points fit
    // 2) No angle requirement between tj1 and tj2
    // 3) Large charge near the intersection point X on tj2
    // tj2       ---X---
    // tj1         /
    // tj1        /
    // tj1       /
    // minimum^2 DOCA of tj1 endpoint with tj2
    
    if(!fUseAlg[kHammerVx2]) return;
   
    unsigned short itj1, end1, endPt1, itj2, closePt2, ipt, ivx;
    unsigned short tjSize = tjs.allTraj.size();
    unsigned short maxPtsFit = 0;
    unsigned short numPtsWithCharge1, numPtsWithCharge2;
    float doca, minDOCA = 3, nPtsFitFrac, intChg, cnt, chgPull;
    bool didaSplit;
    float wint, tint;
    for(itj1 = 0; itj1 < tjSize; ++itj1) {
      if(tjs.allTraj[itj1].AlgMod[kKilled]) continue;
      Trajectory& tj1 = tjs.allTraj[itj1];
      numPtsWithCharge1 = NumPtsWithCharge(tj1, false);
      if(numPtsWithCharge1 < 6) continue;
      // Require that most of the points in tj1 are fitted to a straight line
      maxPtsFit = 0;
      for(ipt = 0; ipt < tj1.Pts.size(); ++ipt) if(tj1.Pts[ipt].Chg > 0 && tj1.Pts[ipt].NTPsFit > maxPtsFit) maxPtsFit = tj1.Pts[ipt].NTPsFit;
      nPtsFitFrac = (float)maxPtsFit / (float)numPtsWithCharge1;
//      std::cout<<"FHV2 "<<tjs.allTraj[itj1].ID<<" "<<maxPtsFit<<" "<<numPtsWithCharge1<<" frac "<<nPtsFitFrac<<"\n";
      if(nPtsFitFrac < 0.8) continue;
      if(vtxPrt) mf::LogVerbatim("TC")<<"FindHammerVertices2: tj1 "<<tj1.ID;
      // Check each end of tj1
      didaSplit = false;
      for(end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tj1.VtxID[end1] > 0) continue;
        if(end1 == 0) {
          endPt1 = tj1.EndPt[0];
        } else {
          endPt1 = tj1.EndPt[1];
        }
        for(itj2 = 0; itj2 < tjSize; ++itj2) {
          if(itj1 == itj2) continue;
          Trajectory& tj2 = tjs.allTraj[itj2];
          if(tj2.AlgMod[kKilled]) continue;
          // require that both be in the same CTP
          if(tj2.CTP != tj1.CTP) continue;
          numPtsWithCharge2 = NumPtsWithCharge(tj2, true);
          if(numPtsWithCharge2 < 6) continue;
          // ignore if tj1 is a lot shorter than tj2
          // ignore if ChgRMS isn't known
          if(tjs.allTraj[itj2].ChgRMS == 0) continue;
          if(vtxPrt) mf::LogVerbatim("TC")<<" Candidate? "<<tj1.ID<<"  NPWC1 "<<numPtsWithCharge1<<" NPWC2 "<<numPtsWithCharge2;
          if(numPtsWithCharge1 < 0.2 * numPtsWithCharge2) continue;
          // Find the minimum separation between tj1 and tj2
          doca = 5;
          TrajPointTrajDOCA(tjs, tj1.Pts[endPt1], tj2, closePt2, doca);
          if(vtxPrt) mf::LogVerbatim("TC")<<" Candidate? "<<tj1.ID<<"  "<<tj2.ID<<" doca "<<doca<<" tj2.EndPt[0] "<<tj2.EndPt[0]<<" closePt2 "<<closePt2<<" tj2.EndPt[1] "<<tj2.EndPt[1];
          if(doca == minDOCA) continue;
          // ensure that the closest point is not near an end
          if(closePt2 < tj2.EndPt[0] + 3) continue;
          if(closePt2 > tj2.EndPt[1] - 3) continue;
          // Find the intersection point between the tj1 end and tj2 closest Pt
          TrajIntersection(tj1.Pts[endPt1], tj2.Pts[closePt2], wint, tint);
          if(vtxPrt) mf::LogVerbatim("TC")<<" intersection W:T "<<wint<<"  "<<tint<<" ticks "<<tint/tjs.UnitsPerTick;
          // Find the point on tj2 that is closest to this point, overwriting closePt
          doca = minDOCA;
          TrajClosestApproach(tj2, wint, tint, closePt2, doca);
          if(vtxPrt) mf::LogVerbatim("TC")<<" closePt2 on tj2 "<<closePt2<<" doca "<<doca;
          if(doca == minDOCA) continue;
          // require a lot of charge in tj2 in this vicinity compared with the average.
          intChg = 0;
          cnt = 0;
          for(ipt = closePt2 - 2; ipt < closePt2 + 2; ++ipt) {
            if(tjs.allTraj[itj2].Pts[ipt].Chg == 0) continue;
            intChg += tjs.allTraj[itj2].Pts[ipt].Chg;
            ++cnt;
          } // ipt
          if(cnt > 0) {
            intChg /= cnt;
            chgPull = (intChg - tjs.allTraj[itj2].AveChg) / tjs.allTraj[itj2].ChgRMS;
            if(vtxPrt) mf::LogVerbatim("TC")<<" chgPull at intersection point "<<chgPull;
            if(chgPull < 10) continue;
          }
          // we have a winner
          // create a new vertex
          VtxStore aVtx;
          aVtx.Pos = tj2.Pts[closePt2].Pos;
          aVtx.NTraj = 3;
          aVtx.Pass = tj2.Pass;
          aVtx.Topo = 6;
          aVtx.ChiDOF = 0;
          aVtx.CTP = fCTP;
          tjs.vtx.push_back(aVtx);
          ivx = tjs.vtx.size() - 1;
          tjs.vtx[ivx].ID = ivx + 1;
          if(!SplitAllTraj(tjs, itj2, closePt2, ivx, vtxPrt)) {
            if(vtxPrt) mf::LogVerbatim("TC")<<"FindHammerVertices2: Failed to split trajectory";
            tjs.vtx.pop_back();
            continue;
          }
          tjs.allTraj[itj1].VtxID[end1] = tjs.vtx[ivx].ID;
          tjs.allTraj[itj1].AlgMod[kHammerVx2] = true;
          tj2.AlgMod[kHammerVx] = true;
          tjs.allTraj[tjs.allTraj.size()-1].AlgMod[kHammerVx2] = true;
          didaSplit = true;
          break;
        } // itj2
        if(didaSplit) break;
      } // end1
    } // itj1
  } // FindHammerVertices2
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindHammerVertices()
  {
    // Look for a trajectory that intersects another. Split
    // the trajectory and make a vertex. The convention used
    // is shown pictorially here. Trajectory tj1 must be longer
    // than tj2
    // tj2       ------
    // tj1         /
    // tj1        /
    // tj1       /
    
    if(!fUseAlg[kHammerVx]) return;
    
    unsigned short itj1, end1, endPt1;
    unsigned short itj2, closePt2, end20, end21, ivx;
    unsigned short tjSize = tjs.allTraj.size();
    
    // minimum^2 DOCA of tj1 endpoint with tj2
    float doca, minDOCA = 3, dang;
    bool didaSplit;
    unsigned short tj1len, tj2len;
    
    for(itj1 = 0; itj1 < tjSize; ++itj1) {
      if(tjs.allTraj[itj1].AlgMod[kKilled]) continue;
      // minimum length requirements
      tj1len = tjs.allTraj[itj1].EndPt[1] - tjs.allTraj[itj1].EndPt[0];
      if(tj1len < 6) continue;
      // Check each end of tj1
      didaSplit = false;
      for(end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tjs.allTraj[itj1].VtxID[end1] > 0) continue;
        if(end1 == 0) {
          endPt1 = tjs.allTraj[itj1].EndPt[0];
        } else {
          endPt1 = tjs.allTraj[itj1].EndPt[1];
        }
        for(itj2 = 0; itj2 < tjSize; ++itj2) {
          if(itj1 == itj2) continue;
          Trajectory& tj2 = tjs.allTraj[itj2];
          if(tj2.AlgMod[kKilled]) continue;
          // require that both be in the same CTP
          if(tj2.CTP != tjs.allTraj[itj1].CTP) continue;
          // length of tj2 cut
          tj2len = tj2.EndPt[1] - tj2.EndPt[0];
          if(tj2len < 6) continue;
          // ignore if tj1 is a lot shorter than tj2
          if(tj1len < 0.5 * tj2len) continue;
          // ignore very long straight trajectories (probably a cosmic muon)
          end20 = tj2.EndPt[0];
          end21 = tj2.EndPt[1];
          if(tj2len > 100 && DeltaAngle(tj2.Pts[end20].Ang, tj2.Pts[end21].Ang) < 0.2) continue;
          // Require no vertex associated with itj2
          if(tj2.VtxID[0] > 0 || tj2.VtxID[1] > 0) continue;
          doca = minDOCA;
          TrajPointTrajDOCA(tjs, tjs.allTraj[itj1].Pts[endPt1], tj2, closePt2, doca);
//          std::cout<<"FHV "<<tjs.allTraj[itj1].ID<<"  "<<tj2.ID<<" doca "<<doca<<" tj2.EndPt[0] "<<tj2.EndPt[0]<<" closePt2 "<<closePt2<<" tj2.EndPt[1] "<<tj2.EndPt[1]<<" doca "<<doca<<"\n";
          if(doca == minDOCA) continue;
          // ensure that the closest point is not near an end
          if(vtxPrt) mf::LogVerbatim("TC")<<"FindHammerVertices: Candidate "<<tjs.allTraj[itj1].ID<<"  "<<tj2.ID<<" doca "<<doca<<" tj2.EndPt[0] "<<tj2.EndPt[0]<<" closePt2 "<<closePt2<<" tj2.EndPt[1] "<<tj2.EndPt[1];
          if(closePt2 < tj2.EndPt[0] + 3) continue;
          if(closePt2 > tj2.EndPt[1] - 3) continue;
          // make an angle cut
          dang = DeltaAngle(tjs.allTraj[itj1].Pts[endPt1].Ang, tj2.Pts[closePt2].Ang);
          if(vtxPrt) mf::LogVerbatim("TC")<<" dang "<<dang<<" fKinkAngCut "<<fKinkAngCut;
          if(dang < fKinkAngCut) continue;
          // we have a winner
          // create a new vertex
          VtxStore aVtx;
          aVtx.Pos = tj2.Pts[closePt2].Pos;
          aVtx.NTraj = 3;
          aVtx.Pass = tj2.Pass;
          aVtx.Topo = 6;
          aVtx.ChiDOF = 0;
          aVtx.CTP = fCTP;
          tjs.vtx.push_back(aVtx);
          ivx = tjs.vtx.size() - 1;
          tjs.vtx[ivx].ID = ivx + 1;
          if(!SplitAllTraj(tjs, itj2, closePt2, ivx, vtxPrt)) {
            if(vtxPrt) mf::LogVerbatim("TC")<<"FindHammerVertices: Failed to split trajectory";
            tjs.vtx.pop_back();
            continue;
          }
          tjs.allTraj[itj1].VtxID[end1] = tjs.vtx[ivx].ID;
          tjs.allTraj[itj1].AlgMod[kHammerVx] = true;
          tj2.AlgMod[kHammerVx] = true;
          tjs.allTraj[tjs.allTraj.size()-1].AlgMod[kHammerVx] = true;
          didaSplit = true;
          break;
        } // tj2
        if(didaSplit) break;
      } // end1
    } // tj1
    
  } // FindHammerVertices
  
  //////////////////////////////////////
  void TrajClusterAlg::Find3DVertices(geo::TPCID const& tpcid)
  {
    // Create 3D vertices from 2D vertices. 3D vertices that are matched
    // in all three planes have Ptr2D >= 0 for all planes
    
    
    if(fVertex3DChiCut < 0) return;
    if(tjs.vtx.size() < 2) return;
    
    geo::TPCGeo const& TPC = geom->TPC(tpcid);
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    const unsigned int cstat = tpcid.Cryostat;
    const unsigned int tpc = tpcid.TPC;

    // create a array/vector of 2D vertex indices in each plane
    std::vector<std::vector<unsigned short>> vIndex(3);
    unsigned short ipl;
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      if(tjs.vtx[ivx].NTraj == 0) continue;
      geo::PlaneID iplID = DecodeCTP(tjs.vtx[ivx].CTP);
      if(iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
      ipl = iplID.Plane;
      if(ipl > 2) continue;
      vIndex[ipl].push_back(ivx);
    }
    
    unsigned short vtxInPln = 0;
    for(unsigned short ipl = 0; ipl < TPC.Nplanes(); ++ipl) if(vIndex[ipl].size() > 0) ++vtxInPln;
    if(vtxInPln < 2) return;

    // Y,Z limits of the detector
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    
    const geo::TPCGeo &thetpc = geom->TPC(tpc, cstat);
    thetpc.LocalToWorld(local,world);
    // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
    float YLo = world[1]-geom->DetHalfHeight(tpc,cstat) + 1;
    float YHi = world[1]+geom->DetHalfHeight(tpc,cstat) - 1;
    float ZLo = world[2]-geom->DetLength(tpc,cstat)/2 + 1;
    float ZHi = world[2]+geom->DetLength(tpc,cstat)/2 - 1;
    
    vtxPrt = (debug.Plane >= 0) && (debug.Tick == 6666);
    
    if(vtxPrt) {
      mf::LogVerbatim("CC")<<"Inside Find3DVertices";
      PrintAllTraj("F3DV", tjs, debug, USHRT_MAX, tjs.allTraj.size());
    }
    
    // wire spacing in cm
    float wirePitch = geom->WirePitch(0, 1, 0, tpcid.TPC, tpcid.Cryostat);
    
    size_t vsize = tjs.vtx.size();
    // vector of 2D vertices -> 3D vertices.
    std::vector<short> vPtr(vsize, -1);
    // fill temp vectors of 2D vertex X and X errors
    std::vector<float> vX(vsize);
    std::vector<float> vXsigma(vsize);
    float vXp;
    
    double ticks;
    for(unsigned short ivx = 0; ivx < vsize; ++ivx) {
      if(tjs.vtx[ivx].NTraj == 0) continue;
      geo::PlaneID iplID = DecodeCTP(tjs.vtx[ivx].CTP);
      if(iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
      // Convert 2D vertex time error to X error
      ticks = tjs.vtx[ivx].Pos[1] / tjs.UnitsPerTick;
      vX[ivx]  = detprop->ConvertTicksToX(ticks, (int)iplID.Plane, (int)tpc, (int)cstat);
      ticks = (tjs.vtx[ivx].Pos[1] + tjs.vtx[ivx].PosErr[1]) / tjs.UnitsPerTick;
      vXp = detprop->ConvertTicksToX(ticks, (int)iplID.Plane, (int)tpc, (int)cstat);
      vXsigma[ivx] = fabs(vXp - vX[ivx]);
    } // ivx
    
    // temp vector of all 2D vertex matches
    std::vector<Vtx3Store> v3temp;
    
    double y = 0, z = 0;
    TVector3 WPos = {0, 0, 0};
    TrajPoint tp;
    // i, j, k indicates 3 different wire planes
    unsigned short ii, jpl, jj, kpl, kk, ivx, jvx, kvx, i3t;
    float kX, kWire, kChi, dX, dXChi, dXSigma, dW;
    bool gotit, sigOK;
    // compare vertices in each view
    for(ipl = 0; ipl < 2; ++ipl) {
      for(ii = 0; ii < vIndex[ipl].size(); ++ii) {
        ivx = vIndex[ipl][ii];
        if(ivx > tjs.vtx.size() - 1) {
          mf::LogError("CC")<<"Find3DVertices: bad ivx "<<ivx;
          return;
        }
        // vertex has been matched already
        if(vPtr[ivx] >= 0) continue;
        unsigned int iWire = std::nearbyint(tjs.vtx[ivx].Pos[0]);
        for(jpl = ipl + 1; jpl < 3; ++jpl) {
          for(jj = 0; jj < vIndex[jpl].size(); ++jj) {
            jvx = vIndex[jpl][jj];
            if(jvx > tjs.vtx.size() - 1) {
              mf::LogError("CC")<<"Find3DVertices: bad jvx "<<jvx;
              return;
            }
            // vertex has been matched already
            if(vPtr[jvx] >= 0) continue;
            unsigned int jWire = std::nearbyint(tjs.vtx[jvx].Pos[0]);
            // new stuff
            dX = fabs(vX[ivx] - vX[jvx]);
            dXSigma = sqrt(vXsigma[ivx] * vXsigma[ivx] + vXsigma[jvx] * vXsigma[jvx]);
            dXChi = dX / dXSigma;
            
            if(vtxPrt) mf::LogVerbatim("CC")<<"Find3DVertices: ipl "<<ipl<<" ivx "<<ivx<<" ivX "<<vX[ivx]
              <<" jpl "<<jpl<<" jvx "<<jvx<<" jvX "<<vX[jvx]<<" W:T "<<(int)tjs.vtx[jvx].Pos[0]<<":"<<(int)tjs.vtx[jvx].Pos[1]<<" dXChi "<<dXChi<<" fVertex3DChiCut "<<fVertex3DChiCut;
            
            if(dXChi > fVertex3DChiCut) continue;
            if (geom->HasWire(geo::WireID(cstat, tpc, ipl, iWire))&&
                geom->HasWire(geo::WireID(cstat, tpc, jpl, jWire))){
              geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
            }
            else continue;
            if(y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
            WPos[1] = y;
            WPos[2] = z;
            kpl = 3 - ipl - jpl;
            kX = 0.5 * (vX[ivx] + vX[jvx]);
            kWire = -1;
            if(TPC.Nplanes() > 2) {
              kWire = geom->NearestWire(WPos, kpl, tpc, cstat);
              tp.Pos[0] = kWire;
              // See if there is a wire signal nearby in kpl
              tp.Pos[1] = detprop->ConvertXToTicks(kX, kpl, fTpc, fCstat) * tjs.UnitsPerTick;
              tp.CTP = EncodeCTP(fCstat, fTpc, kpl);
              sigOK = SignalAtTp(tp);
              if(vtxPrt) mf::LogVerbatim("TC")<<" 3rd plane pos "<<kWire<<" "<<tp.Pos[1]<<" sig present? "<<sigOK;
              if(!sigOK) continue;
            }
            kpl = 3 - ipl - jpl;
            // save this incomplete 3D vertex
            Vtx3Store v3d;
            v3d.Ptr2D[ipl] = ivx;
            v3d.Ptr2D[jpl] = jvx;
            v3d.Ptr2D[kpl] = -1;
            // see if this is already in the list
            gotit = false;
            for(i3t = 0; i3t < v3temp.size(); ++i3t) {
              if(v3temp[i3t].Ptr2D[0] == v3d.Ptr2D[0] && v3temp[i3t].Ptr2D[1] == v3d.Ptr2D[1] && v3temp[i3t].Ptr2D[2] == v3d.Ptr2D[2]) {
                gotit = true;
                break;
              }
            } // i3t
            if(gotit) continue;
            v3d.X = kX;
            // Use XErr to store dXChi
            v3d.XErr = dXChi;
            v3d.Y = y;
            float yzSigma = wirePitch * sqrt(tjs.vtx[ivx].PosErr[0] * tjs.vtx[ivx].PosErr[0] + tjs.vtx[jvx].PosErr[0] * tjs.vtx[jvx].PosErr[0]);
            v3d.YErr = yzSigma;
            v3d.Z = z;
            v3d.ZErr = yzSigma;
            v3d.Wire = kWire;
            v3d.CStat = cstat;
            v3d.TPC = tpc;
            // push the incomplete vertex onto the list
            v3temp.push_back(v3d);
            
            if(vtxPrt) mf::LogVerbatim("CC")<<"Find3DVertices: 2 Plane match ivx "<<ivx<<" P:W:T "<<ipl<<":"<<(int)tjs.vtx[ivx].Pos[0]<<":"<<(int)tjs.vtx[ivx].Pos[1]<<" jvx "<<jvx<<" P:W:T "<<jpl<<":"<<(int)tjs.vtx[jvx].Pos[0]<<":"<<(int)tjs.vtx[jvx].Pos[1]<<" dXChi "<<dXChi<<" yzSigma "<<yzSigma;
            
            if(TPC.Nplanes() == 2) continue;

            // See if the expected position of the vertex in the 3rd view
            // look for a 3 plane match
            for(kk = 0; kk < vIndex[kpl].size(); ++kk) {
              kvx = vIndex[kpl][kk];
              if(vPtr[kvx] >= 0) continue;
              // Wire difference error
              dW = wirePitch * (tjs.vtx[kvx].Pos[0] - kWire) / yzSigma;
              // X difference error
              dX = (vX[kvx] - kX) / dXSigma;
              kChi = 0.5 * sqrt(dW * dW + dX * dX);
              if(kChi < fVertex3DChiCut) {
                // push all complete vertices onto the list
                v3d.X = (vX[kvx] + 2 * kX) / 3;
                v3d.XErr = kChi;
                v3d.Ptr2D[kpl] = kvx;
                // see if this is already in the list
                gotit = false;
                for(i3t = 0; i3t < v3temp.size(); ++i3t) {
                  if(v3temp[i3t].Ptr2D[0] == v3d.Ptr2D[0] && v3temp[i3t].Ptr2D[1] == v3d.Ptr2D[1] && v3temp[i3t].Ptr2D[2] == v3d.Ptr2D[2]) {
                    gotit = true;
                    break;
                  }
                } // i3t
                if(gotit) continue;
                v3temp.push_back(v3d);
                if(vtxPrt) mf::LogVerbatim("CC")<<" kvx "<<kvx<<" kpl "<<kpl
                  <<" wire "<<(int)tjs.vtx[kvx].Pos[0]<<" kTime "<<(int)tjs.vtx[kvx].Pos[1]<<" kChi "<<kChi<<" dW "<<tjs.vtx[kvx].Pos[0] - kWire;
              } // kChi < best
            } // kk
          } // jj
        } // jpl
      } // ii
    } // ipl
    
    if(v3temp.empty()) return;
    
    if(vtxPrt) {
      mf::LogVerbatim("TC")<<"v3temp list";
      for(auto& v3d : v3temp) {
        mf::LogVerbatim("TC")<<v3d.Ptr2D[0]<<" "<<v3d.Ptr2D[1]<<" "<<v3d.Ptr2D[2]<<" wire "<<v3d.Wire<<" "<<v3d.XErr;
      } // v3d
    }
    
    // find the best 3 plane matches
    if(TPC.Nplanes() > 2) {
      float best;
      unsigned short imbest;
      bool gotone = true;
      std::vector<Vtx3Store> v3temp2;
      while(gotone) {
        gotone = false;
        best = fVertex3DChiCut;
        imbest = USHRT_MAX;
        for(ivx = 0; ivx < v3temp.size(); ++ivx) {
          if(v3temp[ivx].Wire < 0) continue;
          if(v3temp[ivx].XErr > best) continue;
          best = v3temp[ivx].XErr;
          imbest = ivx;
        } // ivx
        if(imbest == USHRT_MAX) break;
        // stash the best one temporarily in a second vector
        v3temp2.push_back(v3temp[imbest]);
        // set the XErr large so it isn't checked again
        v3temp[imbest].XErr = 9999;
        // look for these trajectory IDs in other vertices (complete or incomplete)
        for(ivx = 0; ivx < v3temp.size(); ++ivx) {
          if(ivx == imbest) continue;
          for(ipl = 0; ipl < 3; ++ipl) {
            if(v3temp[ivx].Ptr2D[ipl] == v3temp[imbest].Ptr2D[ipl]) {
              v3temp[ivx].XErr = 9999;
              gotone = true;
            }
          } // ipl
        } // ivx
      } // gotone
      // concatenate the two vectors
      v3temp.insert(v3temp.end(), v3temp2.begin(), v3temp2.end());
    } // TPC.Nplanes() > 2

    if(vtxPrt) {
      mf::LogVerbatim("TC")<<"v3temp list after";
      for(auto& v3d : v3temp) {
        if(v3d.XErr > fVertex3DChiCut) continue;
        mf::LogVerbatim("TC")<<v3d.Ptr2D[0]<<" "<<v3d.Ptr2D[1]<<" "<<v3d.Ptr2D[2]<<" wire "<<v3d.Wire<<" "<<v3d.XErr;
      } // v3d
    }

    // Store the vertices
    unsigned short ninc = 0;
    for(auto& v3d : v3temp) {
      if(v3d.XErr > fVertex3DChiCut) continue;
      if(TPC.Nplanes() == 2) {
        v3d.Ptr2D[2] = 666;
      } else {
        if(v3d.Wire >= 0) ++ninc;
      }
      if(vtxPrt) mf::LogVerbatim("CC")<<"3D vtx "<<tjs.vtx3.size()<<" Ptr2D "<<v3d.Ptr2D[0]<<" "<<v3d.Ptr2D[1]<<" "<<v3d.Ptr2D[2]
        <<" wire "<<v3d.Wire;
      tjs.vtx3.push_back(v3d);
    } // ivx

    // Try to complete incomplete vertices
//    if(ninc > 0) CompleteIncomplete3DVertices(tpcid);
    
  } // Find3DVertices
  
  //////////////////////////////////////////
  void TrajClusterAlg::CompleteIncomplete3DVertices(geo::TPCID const& tpcid)
  {
    // Look for trajectories in a plane that lack a 2D vertex as listed in
    // Ptr2D that are near the projected wire. This may trigger splitting trajectories,
    // assigning them to a new 2D vertex and completing 3D vertices

    geo::TPCGeo const& TPC = geom->TPC(tpcid);
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    unsigned short ipl, itj, mPlane, closePt;
    unsigned short nEndPt, ii, end, aVtxIndx;
    unsigned short dpt;
    CTP_t mCTP;
    // A TP for the missing 2D vertex
    TrajPoint tp;
    // 2D vertex
    VtxStore aVtx;
    float doca, maxdoca = 6;
    // list of candidate trajectory indices and ipt index
    std::vector<std::pair<unsigned short, unsigned short>> mTjs;
    for(auto& vx3 : tjs.vtx3) {
      // check for a completed 3D vertex
      if(vx3.Wire < 0) continue;
      mPlane = USHRT_MAX;
      for(ipl = 0; ipl < TPC.Nplanes(); ++ipl) {
        if(vx3.Ptr2D[ipl] >= 0) continue;
        mPlane = ipl;
        break;
      } // ipl
      if(mPlane == USHRT_MAX) continue;
      mCTP = EncodeCTP(vx3.CStat, vx3.TPC, mPlane);
      // X position of the purported missing vertex
      tp.Pos[0] = vx3.Wire;
      tp.Pos[1] = detprop->ConvertXToTicks(vx3.X, mPlane, vx3.TPC, vx3.CStat) * tjs.UnitsPerTick;
      for(itj = 0; itj < tjs.allTraj.size(); ++itj) {
        if(tjs.allTraj[itj].CTP != mCTP) continue;
        if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
        doca = maxdoca;
        // find the closest distance between the vertex and the trajectory
        TrajPointTrajDOCA(tjs, tp, tjs.allTraj[itj], closePt, doca);
        if(prt) mf::LogVerbatim("TC")<<"CI3DV itj ID "<<tjs.allTraj[itj].ID<<" closePT "<<closePt<<" doca "<<doca;
        if(closePt > tjs.allTraj[itj].EndPt[1]) continue;
        mTjs.push_back(std::make_pair(itj, closePt));
      } // itj
      // handle the case where there are one or more TJs with TPs near the ends
      // that make a vertex (a failure by Find2DVertices)
      if(mTjs.empty()) continue;
      aVtxIndx = tjs.vtx.size();
      aVtx.CTP = mCTP;
      aVtx.Topo = 6;
      aVtx.NTraj = mTjs.size();
      // Give it a bogus pass to indicate it wasn't created while stepping
      aVtx.Pass = 9;
      aVtx.Pos = tp.Pos;
      nEndPt = 0;
      for(ii = 0; ii < mTjs.size(); ++ii) {
        itj = mTjs[ii].first;
        closePt = mTjs[ii].second;
        // determine which end is the closest
        if(fabs(closePt - tjs.allTraj[itj].EndPt[0]) < fabs(closePt - tjs.allTraj[itj].EndPt[1])) {
          // closest to the beginning
          end = 0;
        } else {
          // closest to the end
          end = 1;
        }// closest to the end
        dpt = fabs(closePt - tjs.allTraj[itj].EndPt[end]);
        if(dpt < 4) {
          // close to an end
          tjs.allTraj[itj].VtxID[end] = tjs.vtx[aVtxIndx].ID;
          tjs.allTraj[itj].AlgMod[kComp3DVx] = true;
          ++nEndPt;
        } else {
          // closePt is not near an end, so split the trajectory
          if(!SplitAllTraj(tjs, itj, closePt, aVtxIndx, vtxPrt)) {
            mf::LogVerbatim("TC")<<"SplitAllTraj failed. Traj ID "<<tjs.allTraj[itj].ID;
            // failure occurred. Recover
            for(auto mTj : mTjs) {
              unsigned short jtj = mTj.first;
              if(tjs.allTraj[jtj].VtxID[0] == tjs.vtx[aVtxIndx].ID) tjs.allTraj[jtj].VtxID[0] = 0;
              if(tjs.allTraj[jtj].VtxID[1] == tjs.vtx[aVtxIndx].ID) tjs.allTraj[jtj].VtxID[1] = 0;
            } // itj
            continue;
          } // !SplitAllTraj
        } // closePt is not near an end, so split the trajectory
        tjs.allTraj[itj].AlgMod[kComp3DVx] = true;
        itj = tjs.allTraj.size() - 1;
        tjs.allTraj[itj].AlgMod[kComp3DVx] = true;
      } // ii
      tjs.vtx.push_back(aVtx);
      unsigned short ivx = tjs.vtx.size() - 1;
      tjs.vtx[ivx].ID = ivx + 1;
      std::cout<<"CIC new vtx "<<tjs.vtx[ivx].ID<<"\n";
      vx3.Ptr2D[mPlane] = aVtxIndx;
      vx3.Wire = -1;
      if(prt) mf::LogVerbatim("TC")<<"CompleteIncomplete3DVertices: new 2D tjs.vtx "<<aVtxIndx<<" points to 3D tjs.vtx ";
    } // vx3
  } // CompleteIncomplete3DVertices

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
    
    unsigned int step;
    float stepSize;
    unsigned short nMissedSteps = 0;

    bool sigOK, keepGoing;
    unsigned short killPts;
    for(step = 1; step < 10000; ++step) {
      unsigned short angRange = AngleRange(ltp);
      if(angRange > 1) { stepSize = 2; } else { stepSize = std::abs(1/ltp.Dir[0]); }
      // make a copy of the previous TP
      lastPt = tj.Pts.size() - 1;
      tp = tj.Pts[lastPt];
      ++tp.Step;
      // move the local TP position by one step in the right direction
      for(unsigned short iwt = 0; iwt < 2; ++iwt) ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;

      unsigned short ivx = TPNearVertex(tjs, ltp);
      if(ivx != USHRT_MAX) {
        // Trajectory stops near a vertex so make the assignment
        tj.AlgMod[kStopAtVtx] = true;
        tj.VtxID[1] = tjs.vtx[ivx].ID;
        break;
      }

      // Special handling of very long straight trajectories, e.g. uB cosmic rays
      if(fMuonTag[0] >= 0 && tj.PDGCode == 0 && tj.Pts.size() > (unsigned short)fMuonTag[0]) {
        tj.MCSMom = MCSMom(tjs, tj);
        if(tj.MCSMom > fMuonTag[1]) tj.PDGCode = 13;
        if(prt) mf::LogVerbatim("TC")<<" Check MCSMom "<<tj.MCSMom<<" PDGCode "<<tj.PDGCode;
      }

      // anything really really long must be a muon
      if(tj.Pts.size() > 200) tj.PDGCode = 13;
      // copy this position into tp
      tp.Pos = ltp.Pos;
      tp.Dir = ltp.Dir;
      if(prt) {
        mf::LogVerbatim("TC")<<"StepCrawl "<<step<<" Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" Dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" stepSize "<<stepSize<<" AngleRange "<<angRange;
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
      // If successfull, AddHits has defined UseHit for this TP,
      // set the trajectory endpoints, and defined HitPos.
      lastPtWithUsedHits = tj.EndPt[1];
      if(tj.Pts[lastPt].Hits.empty()) {
        // Require three points with charge on adjacent wires for small angle
        // stepping.
        if(angRange == 0 && lastPt == 2) return;
        // No close hits added.
        ++nMissedSteps;
        // First check for no signal in the vicinity
        if(lastPt > 0) {
          // break if this is a reverse propagate activity and there was no signal (not on a dead wire)
          if(!sigOK && tj.AlgMod[kRevProp]) break;
          // Ensure that there is a signal here after missing a number of steps on a LA trajectory
          if(angRange > 0 && nMissedSteps > 4 && !SignalAtTp(ltp)) break;
          // the last point with hits (used or not) is the previous point
          lastPtWithHits = lastPt - 1;
          float tps = TrajPointSeparation(tj.Pts[lastPtWithHits], ltp);
          float dwc = DeadWireCount(ltp, tj.Pts[lastPtWithHits]);
          float nMissedWires = tps * std::abs(ltp.Dir[0]) - dwc;
          float maxWireSkip = fMaxWireSkipNoSignal;
          if(tj.PDGCode == 13) maxWireSkip = fMuonTag[2];
          if(prt) mf::LogVerbatim("TC")<<" StepCrawl: no signal at ltp "<<PrintPos(tjs, ltp)<<" nMissedWires "<<std::fixed<<std::setprecision(1)<<nMissedWires<<" dead wire count "<<dwc<<" maxWireSkip "<<maxWireSkip<<" tj.PGDCode "<<tj.PDGCode;
          if(nMissedWires > maxWireSkip) break;
        }
        // no sense keeping this TP on tj if no hits were added
        tj.Pts.pop_back();
        continue;
      } // tj.Pts[lastPt].Hits.empty()
      // Found hits at this location so reset the missed steps counter
      nMissedSteps = 0;
      // Update the last point fit, etc using the just added hit(s)
      UpdateTraj(tj);
      // a failure occurred
      if(!fUpdateTrajOK) return;
      if(tj.Pts[lastPt].Chg == 0) {
        // There are points on the trajectory by none used in the last step. See
        // how long this has been going on
        float tps = TrajPointSeparation(tj.Pts[lastPtWithUsedHits], ltp);
        float dwc = DeadWireCount(ltp, tj.Pts[lastPtWithUsedHits]);
        float nMissedWires = tps * std::abs(ltp.Dir[0]) - dwc;
        if(prt)  mf::LogVerbatim("TC")<<" Hits exist on the trajectory but are not used. Missed wires "<<nMissedWires<<" dead wire count "<<(int)dwc;
        // break if this is a reverse propagate activity with no dead wires
        if(tj.AlgMod[kRevProp] && dwc == 0) break;
        if(nMissedWires > fMaxWireSkipWithSignal) break;
        // try this out
        if(!MaskedHitsOK(tj)) {
//          if(prt) PrintTrajectory("MHOK", tjs, tj, lastPt);
          return;
        }
        // Keep stepping
        if(prt) PrintTrajectory("SC", tjs, tj, lastPt);
        continue;
      } // tp.Hits.empty()
      if(tj.Pts.size() == 3) {
        // ensure that the last hit added is in the same direction as the first two.
        // This is a simple way of doing it
        if(PosSep2(tj.Pts[0].HitPos, tj.Pts[2].HitPos) < PosSep2(tj.Pts[0].HitPos, tj.Pts[1].HitPos)) return;
        // ensure that this didn't start as a small angle trajectory and immediately turn
        // into a large angle one
        if(angRange > fMaxAngleRange[tj.Pass]) {
          if(prt) mf::LogVerbatim("TC")<<" Wandered into an invalid angle range. Quit stepping.";
          fGoodTraj = false;
          return;
        }
      } // tj.Pts.size() == 3
      // Ensure that the trajectory meets the angle requirements now that the direction is well known
      if(tj.Pts.size() == 5 && AngleRange(tj.Pts[lastPt]) > fMaxAngleRange[tj.Pass]) {
        if(prt) mf::LogVerbatim("TC")<<" Wandered into an invalid angle range"<<AngleRange(tj.Pts[lastPt])<<" for this pass. Quit stepping.";
        return;
      }
      // Update the local TP with the updated position and direction
      ltp.Pos = tj.Pts[lastPt].Pos;
      ltp.Dir = tj.Pts[lastPt].Dir;
      if(fMaskedLastTP) {
        // see if TPs have been masked off many times and if the
        // environment is clean. If so, return and try with next pass
        // cuts
        if(!MaskedHitsOK(tj)) {
          if(prt) PrintTrajectory("SC", tjs, tj, lastPt);
          return;
        }
        // Don't bother with the rest of the checking below if we
        // set all hits not used on this TP
        if(prt) PrintTrajectory("SC", tjs, tj, lastPt);
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
      if(tj.AlgMod[kGottaKink]) keepGoing = false;
      // See if the Chisq/DOF exceeds the maximum.
      // UpdateTraj should have reduced the number of points fit
      // as much as possible for this pass, so this trajectory is in trouble.
      if(killPts == 0 &&  tj.Pts[lastPt].FitChi > fMaxChi && tj.PDGCode != 13) {
        if(prt) mf::LogVerbatim("TC")<<"   bad FitChi "<<tj.Pts[lastPt].FitChi<<" cut "<<fMaxChi;
        fGoodTraj = (NumPtsWithCharge(tj, true) > fMinPtsFit[tj.Pass]);
        return;
      }
      // print the local tp unless we have killing to do
      if(killPts == 0) {
        if(prt) PrintTrajectory("SC", tjs, tj, lastPt);
      } else {
        MaskTrajEndPoints(tj, killPts);
        if(!fGoodTraj) return;
        unsigned int onWire = (float)(std::nearbyint(tj.Pts[lastPt].Pos[0]));
        float nSteps = (float)(step - tj.Pts[lastPt - killPts].Step);
        if(prt) mf::LogVerbatim("TC")<<"TRP   killing "<<killPts<<" after "<<nSteps<<" steps from prev TP.  Current tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1];
        // move the position
        tj.Pts[lastPt].Pos[0] += nSteps * tj.Pts[lastPt].Dir[0];
        tj.Pts[lastPt].Pos[1] += nSteps * tj.Pts[lastPt].Dir[1];
        if(angRange == 0) {
          // put the TP at the wire position prior to the move
          float dw = onWire - tj.Pts[lastPt].Pos[0];
          tj.Pts[lastPt].Pos[0] = onWire;
          tj.Pts[lastPt].Pos[1] += dw * tj.Pts[lastPt].Dir[1] / tj.Pts[lastPt].Dir[0];
        }
        // copy to the local trajectory point
        ltp.Pos = tj.Pts[lastPt].Pos;
        ltp.Dir = tj.Pts[lastPt].Dir;
        if(prt) mf::LogVerbatim("TC")<<"  New ltp.Pos     "<<ltp.Pos[0]<<" "<<ltp.Pos[1]<<" ticks "<<(int)ltp.Pos[1]/tjs.UnitsPerTick;
        if(!keepGoing) break;
      }
    } // step
    
    if(prt) mf::LogVerbatim("TC")<<"End StepCrawl with "<<step<<" steps. tj size "<<tj.Pts.size()<<" fGoodTraj = "<<fGoodTraj<<" with fTryWithNextPass "<<fTryWithNextPass;

    if(fGoodTraj && fTryWithNextPass) {
      mf::LogVerbatim("TC")<<"StepCrawl: Have fGoodTraj && fTryWithNextPass true. This shouldn't happen. Fixing it.";
      fTryWithNextPass = false;
    }

  } // StepCrawl

  ////////////////////////////////////////////////
  bool TrajClusterAlg::IsGhost(std::vector<unsigned int>& tHits, unsigned short& ofTraj)
  {
    // Called by FindJunkTraj to see if the passed hits are close to an existing
    // trajectory and if so, they will be used in that other trajectory
    
    ofTraj = USHRT_MAX;
    
    if(!fUseAlg[kUseGhostHits]) return false;
    
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
    // Check the quality of the trajectory and possibly trim it.
    
    if(!fGoodTraj) return;
    
    fTryWithNextPass = false;

    // ensure that the end points are defined
    SetEndPoints(tjs, tj);
    if(tj.EndPt[0] == tj.EndPt[1]) return;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"inside CheckTraj";
    }
    
    // Ensure that a hit only appears once in the TJ
    if(HasDuplicateHits(tj)) {
      if(prt) mf::LogVerbatim("TC")<<" HasDuplicateHits ";
      fGoodTraj = false;
      return;
    }
    
    unsigned short angRange = AngleRange(tj.Pts[tj.EndPt[1]]);
    // checks are different for Very Large Angle trajectories
    bool isVLA = (angRange == fAngleRanges.size() - 1);
    // The last two ranges are Large Angle and Very Large Angle. Determine if the TJ is Small Angle
    bool isSA = (angRange < fAngleRanges.size() - 2);
    
    // First remove any TPs at the end that have no hits
    // TODO This shouldn't be done but first check to see what code will break
    // if we don't do it.
    tj.Pts.resize(tj.EndPt[1] + 1);

    if(isVLA && HasDuplicateHits(tj)) {
      fGoodTraj = false;
      return;
    }
    
    // Fill in any gaps with hits that were skipped, most likely delta rays on muon tracks
    
    CalculateQuality(tj);
    if(!isVLA) FillGaps(tj);
    
    if(prt) mf::LogVerbatim("TC")<<" CheckTraj MCSMom "<<tj.MCSMom<<" isVLA? "<<isVLA<<" NumPtsWithCharge "<<NumPtsWithCharge(tj, false)<<" Min Req'd "<<fMinPts[tj.Pass];
    
    if(NumPtsWithCharge(tj, false) < fMinPts[tj.Pass]) {
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
      // Require less than a 2X difference in the hit width or 4X for HitPosErr2
      if(maxWidth > 4 * minWidth) {
        if(prt) mf::LogVerbatim("TC")<<" TP width variation too large: minWidth "<<minWidth<<" maxWidth "<<maxWidth;
      fGoodTraj = false;
        return;
      }
    } // short trajectory
    
    // Set the StopsAtEnd flags. This will be used by FixTrajBegin
    ChkStop(tj);

    // Update the trajectory parameters at the beginning of the trajectory
    FixTrajBegin(tj);

    // check the fraction of the trajectory points that have hits
    if(fUseAlg[kTrimHits]) {
      // First ensure that there are at least two points with charge at the end
      while(NumPtsWithCharge(tj, false) > fMinPts[tj.Pass]) {
        unsigned short lastPt = tj.EndPt[1];
        if(tj.Pts[lastPt - 1].Chg == 0) {
          UnsetUsedHits(tj.Pts[lastPt]);
          SetEndPoints(tjs, tj);
          tj.AlgMod[kTrimHits] = true;
        } else {
          break;
        }
      } // tjSize > fMinPts[tj.Pass]
      // Next ensure that there are at least 2 good TPs at the end after a dead wire section.
      unsigned short dwc = DeadWireCount(tj.Pts[tj.EndPt[1]-1], tj.Pts[tj.EndPt[1]]);
      if(prt) mf::LogVerbatim("TC")<<"kTrimHits: dead wire count "<<dwc;
      while(dwc > 0 && NumPtsWithCharge(tj, false) > fMinPts[tj.Pass]) {
        // clobber the TPs (most of which are presumed to be on dead wires) after
        // the last TP that has used hits
        UnsetUsedHits(tj.Pts[tj.EndPt[1]]);
        tj.Pts.resize(tj.EndPt[1]+1);
        SetEndPoints(tjs, tj);
        dwc = DeadWireCount(tj.Pts[tj.EndPt[1]-1], tj.Pts[tj.EndPt[1]]);
        if(prt) mf::LogVerbatim("TC")<<"kTrimHits: trimmed single hit after dead wire section. new dead wire count "<<dwc;
        tj.AlgMod[kTrimHits] = true;
      }
      // Ensure that there are points on three consecutive wires at the end
      if(!isVLA) {
        unsigned short tjSize = tj.Pts.size();
        while(tjSize > fMinPts[tj.Pass]) {
          unsigned short lastPt = tj.Pts.size() - 1;
          unsigned int wire1 = std::nearbyint(tj.Pts[lastPt].Pos[0]);
          unsigned int wire2 = std::nearbyint(tj.Pts[lastPt - 1].Pos[0]);
          unsigned int wire3 = std::nearbyint(tj.Pts[lastPt - 2].Pos[0]);
          if(prt) mf::LogVerbatim("TC")<<" consecutive end wire check "<<wire1<<" "<<wire2<<" "<<wire3;
          if(abs(wire1 - wire2) == 1 && abs(wire2 - wire3) == 1) break;
          UnsetUsedHits(tj.Pts[tj.EndPt[1]]);
          tj.Pts.resize(tj.EndPt[1]+1);
          SetEndPoints(tjs, tj);
          --tjSize;
          tj.AlgMod[kTrimHits] = true;
        } // tj.Pts.size() > fMinPts[tj.Pass]
        if(prt) PrintTrajectory("CT", tjs, tj, USHRT_MAX);
      } // not isVLA

      // impose the requirement that 70% of the trajectory points should have hits with charge.
      float nPts = tj.EndPt[1] - tj.EndPt[0] + 1;
      dwc = DeadWireCount(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
      float nPtsWithCharge = NumPtsWithCharge(tj, false);
      float ptFrac = (nPtsWithCharge + dwc) /(nPts + dwc);
      if(prt) mf::LogVerbatim("TC")<<"kTrimHits: nPts "<<(int)nPts<<" DeadWireCount "<<(int)dwc<<" nPtsWithCharge "<<(int)nPtsWithCharge<<" ptFrac "<<ptFrac;
      while(ptFrac < 0.7 && nPts > 1) {
        // mask off points until this condition is satisfied
        UnsetUsedHits(tj.Pts[tj.EndPt[1]]);
        SetEndPoints(tjs, tj);
        nPts = tj.EndPt[1] - tj.EndPt[0] + 1;
        dwc = DeadWireCount(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
        nPtsWithCharge = NumPtsWithCharge(tj, false);
        if(nPtsWithCharge < fMinPts[tj.Pass]) return;
        ptFrac = (nPtsWithCharge + dwc) /(nPts + dwc);
        tj.AlgMod[kTrimHits] = true;
      } // ptFrac < 0.7 && nPts > 1
      if(prt) mf::LogVerbatim("TC")<<" after trim nPts "<<(int)nPts<<" DeadWireCount "<<(int)dwc<<" nPtsWithCharge "<<(int)nPtsWithCharge<<" ptFrac "<<ptFrac;
    } // fUseAlg[kTrimHits]

    // See if this looks like a ghost trajectory and if so, merge the
    // hits and kill this one
    if(fUseAlg[kUseGhostHits]) {
      auto tHits = PutTrajHitsInVector(tj, kUsedHits);
      unsigned short ofTraj = USHRT_MAX;
      if(IsGhost(tHits, ofTraj)) {
        fGoodTraj = false;
        return;
      }
    } // fUseAlg[kUseGhostHits]
    // ignore short trajectories
    if(tj.EndPt[1] < 4) return;
    
    if(isSA) {
      // Small angle checks

      if(fUseAlg[kCWKink] && tj.EndPt[1] > 8) {
        // look for the signature of a kink near the end of the trajectory.
        // These are: Increasing chisq for the last few hits. Presence of
        // a removed hit near the end. A sudden decrease in the number of
        // TPs in the fit. A change in the average charge of hits. These may
        // not all be present in every situation.
        unsigned short tpGap = USHRT_MAX;
        unsigned short nBigRat = 0;
        float chirat;
        for(unsigned short ii = 1; ii < 5; ++ii) {
          unsigned short ipt = tj.EndPt[1] - 1 - ii;
          if(tj.Pts[ipt].Chg == 0) {
            tpGap = ipt;
            break;
          }
          chirat = tj.Pts[ipt+1].FitChi / tj.Pts[ipt].FitChi;
          if(chirat > 1.5) ++nBigRat;
          if(prt) mf::LogVerbatim("TC")<<"CheckTraj: chirat "<<ipt<<" chirat "<<chirat<<" "<<tj.Pts[ipt-1].FitChi<<" "<<tj.Pts[ipt].FitChi<<" "<<tj.Pts[ipt+1].FitChi;
        } // ii
        if(prt) mf::LogVerbatim("TC")<<"CheckTraj: nBigRat "<<nBigRat<<" tpGap "<<tpGap;
        if(tpGap != USHRT_MAX && nBigRat > 0) {
          unsigned short newSize;
          if(tpGap != USHRT_MAX) {
            newSize = tpGap;
          } else {
            newSize = tj.Pts.size() - 3;
          }
          if(prt) mf::LogVerbatim("TC")<<"  Setting tj UseHits from "<<newSize<<" to "<<tj.Pts.size()-1<<" false";
          for(unsigned short ipt = newSize; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tj.Pts[ipt]);
          SetEndPoints(tjs, tj);
          tj.Pts.resize(newSize);
          tj.AlgMod[kCWKink] = true;
          return;
        }
      } // fUseAlg[kCWKink]

      if(fUseAlg[kCWStepChk]) {
        // Compare the number of steps taken per TP near the beginning and
        // at the end
        short nStepBegin = tj.Pts[2].Step - tj.Pts[1].Step;
        short nStepEnd;
        unsigned short lastPt = tj.Pts.size() - 1;
        unsigned short newSize = tj.Pts.size();
        for(unsigned short ipt = lastPt; ipt > lastPt - 2; --ipt) {
          nStepEnd = tj.Pts[ipt].Step - tj.Pts[ipt - 1].Step;
          if(nStepEnd > 3 * nStepBegin) newSize = ipt;
        }
        if(prt) mf::LogVerbatim("TC")<<"CheckTraj: check number of steps. newSize "<<newSize<<" tj.Pts.size() "<<tj.Pts.size();
        if(newSize < tj.Pts.size()) {
          for(unsigned short ipt = newSize; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tj.Pts[ipt]);
          SetEndPoints(tjs, tj);
          tj.AlgMod[kCWStepChk] = true;
          tj.Pts.resize(newSize);
          return;
        } // newSize < tj.Pts.size()
      } // fUseAlg[kCWStepChk]
    } // isSA
    FindSoftKink(tj);
    
    // Check either large angle or not-large angles
    CheckHiDeltas(tj);
    
    CheckHiMultUnusedHits(tj);
    if(!fGoodTraj || fQuitAlg) return;
    
    // lop off high multiplicity hits at the end
    CheckHiMultEndHits(tj);
    
  } // CheckTraj
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindSoftKink(Trajectory& tj)
  {
    // Looks for a soft kink in the trajectory and truncates it if one is found.
    // This is best done after FixTrajBegin has been called.
    
    if(!fUseAlg[kSoftKink]) return;
    if(tj.Pts.size() < 15) return;
    float dang = DeltaAngle(tj.Pts[tj.EndPt[0]].Ang, tj.Pts[tj.EndPt[1]].Ang);
    
    if(prt) {
      mf::LogVerbatim("TC")<<"FindSoftKink "<<tj.ID<<" dang "<<dang<<" cut "<<0.5 * fKinkAngCut;
    }
    if(dang < 0.5 * fKinkAngCut) return;
    // require at least 5 points fitted at the end of the trajectory
    unsigned short endPt = tj.EndPt[1];
    if(tj.Pts[endPt].NTPsFit < 5) return;
    if(tj.Pts[endPt].NTPsFit > endPt) return;
    // Estimate where where the kink would be
    unsigned short kinkPt = endPt - tj.Pts[endPt].NTPsFit;
    // Require at least 5 points in the trajectory before the kink
    if(kinkPt < 5) return;
    // require very few points fitted in this region compared the number of points prior to it
    if(tj.Pts[kinkPt].NTPsFit > 0.3 * kinkPt) return;
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
    for(unsigned short ipt = atPt; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tj.Pts[ipt]);
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
    
    // assume that all points in the trajectory were fitted to a line
    unsigned short lastPtFit = tj.EndPt[1];
    // Find the last point that includes the first point in the fit
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.EndPt[1] - 1 - ii;
      unsigned short firstPtFit = ipt + 1 - tj.Pts[ipt].NTPsFit;
      if(ipt == 0) break;
      if(firstPtFit > 1) continue;
      lastPtFit = ipt;
      break;
    }
    FixTrajBegin(tj, lastPtFit);
    
  } // FixTrajBegin
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FixTrajBegin(Trajectory& tj, unsigned short atPt)
  {
    // Update the parameters at the beginning of the trajectory starting at point atPt
    
    if(!fUseAlg[kFixEnd]) return;
    // ignore short trajectories
    unsigned short npwc = NumPtsWithCharge(tj, false);
    if(npwc < 6) return;
    // ignore somewhat longer trajectories that are curly
    if(npwc < 10 && tj.MCSMom < 100) return;
    // ignore shower-like trajectories
    if(tj.PDGCode == 12) return;
    // ignore junk trajectories
    if(tj.AlgMod[kJunkTj]) return;
    // ignore stopping trajectories
    if(tj.StopsAtEnd[0]) return;
    
    
    unsigned short firstPt = tj.EndPt[0];
    if(prt) {
      mf::LogVerbatim("TC")<<"FixTrajBegin: atPt "<<atPt<<" firstPt "<<firstPt<<" Stops at end 0? "<<tj.StopsAtEnd[0];
    }
    
    if(atPt == tj.EndPt[0]) return;
    
    // update the trajectory for all the intervening points
    for(unsigned short ipt = firstPt; ipt < atPt; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      tp.Dir = tj.Pts[atPt].Dir;
      tp.Ang = tj.Pts[atPt].Ang;
      tp.AngErr = tj.Pts[atPt].AngErr;
      // Correct the projected time to the wire
      float dw = tp.Pos[0] - tj.Pts[atPt].Pos[0];
      if(tp.Dir[0] != 0) tp.Pos[1] = tj.Pts[atPt].Pos[1] + dw * tp.Dir[1] / tp.Dir[0];
      // Try this out. Drop the TPs with charge and re-find them unless they are large angle
      bool newHits = false;
      if(tp.Chg > 0 && AngleRange(tp) == 0) {
        float chgIn = tp.Chg;
        if(tj.StopsAtEnd[0] && ipt < firstPt + 6 && tp.Hits.size() > 1) {
          // Pick up all of the hits near the end of a stopping TJ. Start by finding
          // a hit that is used in this point
          unsigned int myht = INT_MAX;
          for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
            if(!tp.UseHit[ii]) continue;
            myht = tp.Hits[ii];
            break;
          } // ii
          if(myht == INT_MAX) {
            mf::LogWarning("TC")<<"FixTrajBegin: Didn't find myht";
            fQuitAlg = true;
            return;
          }
          // next find all hits in the multiplet in which myht resides
          std::vector<unsigned int> hitsInMuliplet;
          GetHitMultiplet(myht, hitsInMuliplet);
          // Use the hits in the multiplet if they are associated with the TP and are available
          for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
            unsigned int iht = tp.Hits[ii];
            if(tjs.fHits[iht].InTraj > 0) continue;
            if(std::find(hitsInMuliplet.begin(), hitsInMuliplet.end(), iht) == hitsInMuliplet.end()) continue;
            tp.UseHit[ii] = true;
            tjs.fHits[iht].InTraj = tj.ID;
          } // ii
        } else {
          float maxDelta = 5 * tj.Pts[tj.EndPt[1]].DeltaRMS;
          if(tj.Pts[ipt].Hits.size() == 1) maxDelta *= 1.5;
          UnsetUsedHits(tp);
          bool useChg = true;
          if(tj.StopsAtEnd[0]) useChg = false;
          FindUseHits(tj, ipt, maxDelta, useChg);
        }
        DefineHitPos(tj.Pts[ipt]);
        if(tp.Chg != chgIn) newHits = true;
      }
      tj.Pts[ipt].Delta = PointTrajDOCA(tjs, tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tj.Pts[ipt]);
      tj.Pts[ipt].DeltaRMS = tj.Pts[atPt].DeltaRMS;
      tj.Pts[ipt].NTPsFit = tj.Pts[atPt].NTPsFit;
      tj.Pts[ipt].FitChi = tj.Pts[atPt].FitChi;
      tj.Pts[ipt].AveChg = tj.Pts[atPt].AveChg;
      tj.Pts[ipt].ChgPull = (tj.Pts[ipt].Chg / tj.AveChg - 1) / tj.ChgRMS;
      if(prt) {
        if(newHits) {
          PrintTrajectory("FTB", tjs, tj, ipt);
        } else {
          PrintTrajectory("ftb", tjs, tj, ipt);
        }
      }
    } // ipt
    tj.AlgMod[kFixEnd] = true;
    
  } // FixTrajBegin
  
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FixTrajEnd(Trajectory& tj, unsigned short atPt)
  {
    // Update the parameters at the end of the trajectory starting at point atPt
    
    if(!fUseAlg[kFixEnd]) return;
    // ignore short trajectories
    unsigned short npwc = NumPtsWithCharge(tj, false);
    if(npwc < 6) return;
    // ignore somewhat longer trajectories that are curly
    if(npwc < 10 && tj.MCSMom < 100) return;
    // ignore shower-like trajectories
    if(tj.PDGCode == 12) return;
    // ignore junk trajectories
    if(tj.AlgMod[kJunkTj]) return;
    // ingore stopping trajectories
    if(tj.StopsAtEnd[1]) return;
    
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
   
    if(!fUseAlg[kFillGap]) return;
    
    // Find the min (max) charge that we should allow = (half) double the ChgPullCut
    float minChg = tj.AveChg * (1 - 2 * fChgPullCut * tj.ChgRMS);
    float maxChg = tj.AveChg * (1 + 2 * fChgPullCut * tj.ChgRMS);
    
    // start with the first point that has charge
    unsigned short firstPtWithChg = tj.EndPt[0];
    bool first = true;
    float maxDelta = 1;
    while(firstPtWithChg < tj.EndPt[1]) {
      unsigned short nextPtWithChg = firstPtWithChg + 1;
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
      MakeBareTrajPoint(tjs, tj.Pts[firstPtWithChg], tj.Pts[nextPtWithChg], tp);
      // Find the maximum delta between hits and the trajectory Pos for all
      // hits on this trajectory
      if(first) {
        maxDelta = MaxHitDelta(tj);
        first = false;
      } // first
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
        if(chg < minChg  || chg > maxChg) {
          // don't use these hits after all
          UnsetUsedHits(tj.Pts[mpt]);
          filled = false;
        }
        if(filled) {
          DefineHitPos(tj.Pts[mpt]);
          tj.AlgMod[kFillGap] = true;
          if(prt) PrintTrajPoint("FG", tjs, mpt, tj.StepDir, tj.Pass, tj.Pts[mpt]);
        } // filled
      } // mpt
      firstPtWithChg = nextPtWithChg;
    } // firstPtWithChg
    
  } // FillGaps

  ////////////////////////////////////////////////
  float TrajClusterAlg::MaxHitDelta(Trajectory& tj)
  {
    float delta, md = 0;
    unsigned short ii;
    unsigned int iht;
    for(auto& tp : tj.Pts) {
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        iht = tp.Hits[ii];
        delta = PointTrajDOCA(tjs, iht, tp);
        if(delta > md) md = delta;
      } // ii
    } // pts
    return md;
  } // MaxHitDelta
 
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiDeltas(Trajectory& tj)
  {
    // Mask off hits at the beginning and end of the trajectory if
    // Delta is too high
    
    if(!fUseAlg[kHiEndDelta]) return;

    // don't bother with LA trajectories
    if(AngleRange(tj.Pts[tj.EndPt[1]]) > 0) return;
    
    float drms, pull;
    bool didit;
    unsigned short ipt, usePt;

    // Don't check the end if this appears to be a stopping particle since
    // there can be a lot of scatter near the stopping point.

    // Check the beginning first.
    unsigned short endPt = tj.EndPt[0];
    bool checkEnd = (endPt > 0 && !tj.Pts[0].Hits.empty());
    if(checkEnd) {
      // find the first value of delta RMS that is not the default value
      didit = false;
      usePt = USHRT_MAX;
      for(ipt = tj.EndPt[0] + fNPtsAve + 3; ipt < tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        // ignore the default value
        if(tj.Pts[ipt].DeltaRMS == 0.02) continue;
        usePt = ipt;
        break;
      } // ipt
      if(usePt < tj.EndPt[1]) {
        // Scan from usePt back to the beginning. Stop if delta looks suspiciously large
        // and remove all points from the beginning to this point
        unsigned short ii;
        drms = tj.Pts[usePt].DeltaRMS;
        if(drms < 0.02) drms = 0.02;
        for(ii = 0; ii < tj.Pts.size(); ++ii) {
          ipt = usePt - ii;
          pull = tj.Pts[ipt].Delta / drms;
          if(prt) mf::LogVerbatim("TC")<<"CHD begin "<<ipt<<" "<<fPlane<<":"<<PrintPos(tjs, tj.Pts[ipt])<<" Delta "<<tj.Pts[ipt].Delta<<" pull "<<pull<<" usePt "<<usePt;
          if(pull > 5) {
            // unset all TPs from here to the beginning
            for(unsigned short jpt = 0; jpt < ipt + 1; ++jpt) UnsetUsedHits(tj.Pts[jpt]);
            tj.AlgMod[kHiEndDelta] = true;
            didit = true;
            break;
          } // pull > 5
          if(ipt == 0) break;
          if(didit) break;
        } // ii
        if(tj.AlgMod[kHiEndDelta] == true) SetEndPoints(tjs, tj);
      } // usePt != USHRT_MAX
    }

    // now check the other end using a similar procedure
    endPt = tj.EndPt[1];
    checkEnd = (endPt < tj.Pts.size() - 1 && !tj.Pts[endPt + 1].Hits.empty());
    if(checkEnd) {
      usePt = USHRT_MAX;
      didit = false;
      unsigned short cnt = 0;
      for(ipt = tj.EndPt[1]; ipt > tj.EndPt[0]; --ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        if(tj.Pts[ipt].DeltaRMS < 0.02) continue;
        usePt = ipt;
        ++cnt;
        if(ipt == 0) break;
        if(cnt > 8) break;
      } // ipt
      if(usePt == USHRT_MAX) return;
      drms = tj.Pts[usePt].DeltaRMS;
      if(drms < 0.02) drms = 0.02;
      if(prt) mf::LogVerbatim("TC")<<"CHD end usePt "<<usePt<<" drms "<<drms;
      if(usePt == USHRT_MAX) return;
      for(ipt = usePt; ipt < tj.EndPt[1] + 1; ++ipt) {
        pull = tj.Pts[ipt].Delta / drms;
        if(prt) mf::LogVerbatim("TC")<<"CHD end "<<ipt<<" "<<fPlane<<":"<<PrintPos(tjs, tj.Pts[ipt])<<" Delta "<<tj.Pts[ipt].Delta<<" pull "<<pull;
        if(pull > 5) {
          // unset all TPs from here to the end
          for(unsigned short jpt = ipt; jpt < tj.EndPt[1] + 1; ++jpt) UnsetUsedHits(tj.Pts[jpt]);
          didit = true;
          tj.AlgMod[kHiEndDelta] = true;
          break;
        } // pull > 5
        if(didit) break;
      } // ipt
    }

    if(tj.AlgMod[kHiEndDelta] == true) SetEndPoints(tjs, tj);
    
  } // CheckHiDeltas
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiMultUnusedHits(Trajectory& tj)
  {
    // Check for many unused hits in high multiplicity TPs in work and try to use them
    
    
    if(!fUseAlg[kChkHiMultHits]) return;
    
    // This code might do bad things to short trajectories
    if(NumPtsWithCharge(tj, true) < 6) return;
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
        nHiMultPtUsedHits += NumUsedHits(tj.Pts[stopPt]);
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
    bool sortaLargeAngle = (AngleRange(tj.Pts[ii]) == 1);

    if(prt) mf::LogVerbatim("TC")<<"CHMUH: First InTraj stopPt "<<stopPt<<" fracHiMult "<<fracHiMult<<" fracHitsUsed "<<fracHitsUsed<<" lastMult1Pt "<<lastMult1Pt<<" sortaLargeAngle "<<sortaLargeAngle;
    if(fracHiMult < 0.3) return;
    if(fracHitsUsed > 0.98) return;
    
    float maxDelta = 2.5 * MaxHitDelta(tj);

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
    for(ipt = stopPt + 1; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tj.Pts[ipt]);
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
        if(tjs.fHits[iht].InTraj > 0) continue;
        delta = PointTrajDOCA(tjs, iht, tj.Pts[ipt]);
        if(delta > maxDelta) continue;
        tj.Pts[ipt].UseHit[ii] = true;
        tjs.fHits[iht].InTraj = tj.ID;
        added = true;
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
    if(tj.AlgMod[kGottaKink]) return;
    // trim the end points although this shouldn't happen
    if(tj.EndPt[1] != tj.Pts.size() - 1) tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kChkHiMultHits] = true;
    // Don't try to extend trajectories that were trimmed
    if(tj.AlgMod[kTrimHits]) return;
    if(prt) mf::LogVerbatim("TC")<<"TRP CheckHiMultUnusedHits successfull. Calling StepCrawl to extend it";
    StepCrawl(tj);
    
  } // CheckHiMultUnusedHits

  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiMultEndHits(Trajectory& tj)
  {
    // mask off high multiplicity TPs at the end
    if(!fUseAlg[kChkHiMultEndHits]) return;
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
        UnsetUsedHits(tj.Pts[ipt]);
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
  bool TrajClusterAlg::MaskedHitsOK(Trajectory& tj)
  {
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration (true) or whether to stop and possibly try with the next pass settings (false)
    
    if(!fUseAlg[kUnMaskHits]) return true;
    
    unsigned short lastPt = tj.Pts.size() - 1;
    if(tj.Pts[lastPt].Chg > 0) return true;
    
    // count the number of points w/o used hits and the number with one unused hit
    unsigned short nMasked = 0;
    unsigned short nOneHit = 0;
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - ii;
      if(tj.Pts[ipt].Chg > 0) break;
      unsigned short nUnusedHits = 0;
      for(unsigned short jj = 0; jj < tj.Pts[ipt].Hits.size(); ++jj) {
        unsigned int iht = tj.Pts[ipt].Hits[jj];
        if(tjs.fHits[iht].InTraj == 0) ++nUnusedHits;
      } // jj
      if(nUnusedHits == 1) ++nOneHit;
      ++nMasked;
    } // ii
    
    if(prt) {
      mf::LogVerbatim("TC")<<"MaskedHitsOK:  nMasked "<<nMasked<<" nOneHit "<<nOneHit<<" fMaxWireSkipWithSignal "<<fMaxWireSkipWithSignal;
    }

    if(nMasked < 3 || nOneHit < 3) return true;
    
    unsigned short endPt = tj.EndPt[1];
    if(tj.Pts[endPt].NTPsFit > fMinPtsFit[tj.Pass]) {
      // We missed a number of points. See if the charge is OK on these points and delta isn't too bad
      // and there is only one hit on the tp
      unsigned short nOKChg = 0;
      unsigned short nOKDelta = 0;
      for(unsigned ipt = endPt + 1; ipt < tj.Pts.size(); ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        unsigned short nUnusedHits = 0;
        unsigned short iiIndex = 0;
        for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(tjs.fHits[iht].InTraj == 0) {
            ++nUnusedHits;
            iiIndex = ii;
          }
        } // ii
        if(nUnusedHits > 1) break;
        unsigned int iht = tp.Hits[iiIndex];
        float chgPull = std::abs(tjs.fHits[iht].Integral / tj.Pts[endPt].Chg - 1) / tj.ChgRMS;
        if(chgPull < fChgPullCut) ++nOKChg;
        if(tp.Delta < 1.5 * tj.Pts[endPt].Delta) ++nOKDelta;
      } // ipt
      if(prt) mf::LogVerbatim("TC")<<" nOKChg "<<nOKChg<<" nOKDelta "<<nOKDelta;
      if(nOKChg != nMasked || nOKDelta != nMasked) return true;
      // Reduce the number of points fit to minimum for this pass and include the points
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
        tp.NTPsFit = fMinPtsFit[tj.Pass];
        FitTraj(tj);
        if(prt) PrintTrajectory("MHOK", tjs, tj, ipt);
        tj.AlgMod[kUnMaskHits] = true;
        if(tp.FitChi > 2) return false;
      } // ipt
    }
    
    return true;
/*
    // Check for many masked hits and nearby hits to see if they can be included
    // Only allow this to happen once
    if(!tj.AlgMod[kUnMaskHits] && lastPt > 10 && nMasked > 4 && nMasked == nClose) {
      std::cout<<"Trying it "<<tj.ID<<"\n";
      float maxDelta = 4 * tj.Pts[lastPt].DeltaRMS;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.Pts.size() - ii;
        if(tj.Pts[ipt].Chg > 0) break;
        if(tj.Pts[ipt].Delta > maxDelta) continue;
        if(tj.Pts[ipt].Hits.size() != 1) continue;
        unsigned int iht = tj.Pts[ipt].Hits[0];
        if(tjs.fHits[iht].InTraj > 0) continue;
        tj.Pts[ipt].UseHit[0] = true;
        tjs.fHits[iht].InTraj = tj.ID;
        DefineHitPos(tj.Pts[ipt]);
      } // ii
      SetEndPoints(tjs, tj);
      // try to fit them all
      tj.Pts[lastPt].NTPsFit = NumPtsWithCharge(tj, false);
      if(prt) mf::LogVerbatim("TC")<<" Add close hits and try UpdateTraj ";
      UpdateTraj(tj);
      tj.AlgMod[kUnMaskHits] = true;
      return true;
    }

    // Be a bit more lenient with short trajectories on the first pass if
    // the FitChi is not terribly bad and there is ony one hit associated with the last TP
    if(tj.Pass < (fMinPtsFit.size()-1) && tj.Pts.size() > 5 && tj.Pts.size() < 15 && nMasked < 4
       && tj.Pts[lastPt].FitChi < 2 * fMaxChi && tj.Pts[lastPt].Hits.size() == 1) {
      // set this hit used if it is available
      unsigned int iht = tj.Pts[lastPt].Hits[0];
      if(tjs.fHits[iht].InTraj <= 0) {
        tj.Pts[lastPt].UseHit[0] = true;
        tjs.fHits[iht].InTraj = tj.ID;
      }
      SetEndPoints(tjs, tj);
      PrepareForNextPass(tj);
      if(prt) mf::LogVerbatim("TC")<<"MaskedWorkHitsOK: Try next pass with kinda short, kinda bad trajectory. "<<fTryWithNextPass;
      return false;
    }

    // OK if we haven't exceeded the user cut
    if(nMasked < fMaxWireSkipWithSignal) return true;
    
    // not a lot of close hits try with the next pass settings
    if(nClose < 1.5 * nMasked) {
      // trim the trajectory
      unsigned short newSize = tj.Pts.size() - nMasked;
      if(prt) mf::LogVerbatim("TC")<<"MaskedHitsOK:  Trimming  to size "<<newSize;
      tj.Pts.resize(newSize);
      PrepareForNextPass(tj);
    }
    return false;
*/
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
      FitTraj(tj);
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
  void TrajClusterAlg::FindHit(std::string someText, unsigned int iht)
  {

    // look in tjs.allTraj
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      for(unsigned short ipt = 0; ipt < tjs.allTraj[itj].Pts.size(); ++ipt) {
        TrajPoint& tp = tjs.allTraj[itj].Pts[ipt];
        if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) {
          mf::LogVerbatim("TC")<<"FindHit: found hit from "<<someText<<" "<<tjs.allTraj[itj].CTP<<" "<<PrintHit(tjs.fHits[iht])<<" InTraj "<<tjs.fHits[iht].InTraj<<" tjs.allTraj ID "<<tjs.allTraj[itj].ID<<" trajectory below ";
          PrintTrajectory("FH", tjs, tjs.allTraj[itj], USHRT_MAX);
          return;
        }
      } // ipt
    } // itj
  } // FindHit
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::GottaKink(Trajectory& tj, unsigned short& killPts)
  {
    // This routine requires that EndPt is defined
    killPts = 0;
    
    if(tj.AlgMod[kNoKinkChk]) return;

    unsigned short lastPt = tj.EndPt[1];
    if(lastPt < 6) return;
    if(tj.Pts[lastPt].Chg == 0) return;
    
    // A simple check when there are few points being fit and the TJ is short
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
      if(prt) mf::LogVerbatim("TC")<<"GottaKink Simple check lastPt "<<lastPt<<" prevPtWithHits "<<prevPtWithHits<<" dang "<<dang<<" cut "<<fKinkAngCut;
      // need to be more generous since the angle error isn't being considered - BAD IDEA
      if(dang > fKinkAngCut) {
        killPts = 1;
        tj.AlgMod[kGottaKink] = true;
      }
      // Another case where there are few hits fit just prior to a dead wire
      // section or there were no hits added for several steps or due to a large
      // value of fMaxWireSkipNoSignal. We just added a bogus hit just after this section
      // so the trajectory angle change will be small. Find the angle between the previous
      // point fitted angle and the angle formed by the last two TPs
      if(std::abs(tj.Pts[lastPt-1].Pos[0] - tj.Pts[lastPt].Pos[0]) > 3) {
        TrajPoint tmp;
        MakeBareTrajPoint(tjs, tj.Pts[lastPt-1], tj.Pts[lastPt], tmp);
        dang = DeltaAngle(tmp.Ang, tj.Pts[prevPtWithHits].Ang);
        if(prt) mf::LogVerbatim("TC")<<"GottaKink Simple check after gap lastPt "<<lastPt<<" prevPtWithHits "<<prevPtWithHits<<" dang "<<dang<<" cut "<<fKinkAngCut;
        if(dang > 1.5 * fKinkAngCut) {
          killPts = 1;
          tj.AlgMod[kGottaKink] = true;
        }
      }
      // Few points fit in a long trajectory, indicating that chisq went to pot
      // maybe because we just passed through a bad section. Go back a few points
      // and check
      if(tj.Pts.size() < 10) return;
      if(dang < 0.5 * fKinkAngCut) return;
      for(ii = 1; ii < 6; ++ii) {
        ipt = lastPt - ii;
        if(tj.Pts[ipt].Chg == 0) {
          killPts = ii;
          break;
        }
        if(ipt == 0) break;
      } // ii
      if(killPts > 0) tj.AlgMod[kGottaKink] = true;
      return;
    } // tj.NTPSFit < 4

    
    if(tj.EndPt[1] < 10) return;
    
    unsigned short kinkPt = USHRT_MAX;
    
    // Find the kinkPt which is the third Pt from the end that has charge
    unsigned short cnt = 0;
    unsigned short nHiMultPt = 0;
    unsigned short nHiChg = 0;
    
    for(unsigned short ii = 1; ii < lastPt; ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg == 0) continue;
      ++cnt;
      if(tj.Pts[ipt].Hits.size() > 1) ++nHiMultPt;
      if(tj.Pts[ipt].ChgPull > 1.5) ++nHiChg;
      if(cnt == 3) {
        kinkPt = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(kinkPt == USHRT_MAX) return;

    TrajPoint tpFit;
    unsigned short npts = 4;
    unsigned short fitDir = -1;
    FitTraj(tj, lastPt, npts, fitDir, tpFit);
    if(tpFit.FitChi > 900) return;
 
    float dang = DeltaAngle(tj.Pts[kinkPt].Ang, tpFit.Ang);

    float kinkSig = dang / tpFit.AngErr;
    
    if(dang > fKinkAngCut && kinkSig > 3) {
      // Kink larger than hard user cut
      killPts = 3;
    } else if(dang > 0.8 * fKinkAngCut && kinkSig > 10 && tpFit.FitChi < 2) {
      // Found a very significant smaller angle kink
      killPts = 3;
    }
    
    if(killPts > 0) {
      // This could be two crossing trajectories, in which case calling
      // this is a kink would be wrong. Instead we should kill some points
      // and keep going. The signature of crossing trajectories is a multiplicity > 1
      // near the purported kink point and/or increased charge and a long trajectory
      if(tj.Pts.size() < 100 && nHiMultPt == 0) tj.AlgMod[kGottaKink] = true;
      // See if the points we want to kill have higher than expected charge, in which case
      // we want to mask them off but keep stepping
      if(nHiChg > 1) tj.AlgMod[kGottaKink] = false;
      // See if we are tracking a low momentum particle in which case we should just
      // turn off kink checking
      if(tj.EndPt[1] < 20) {
        // Find MCSMom if it hasn't been done
        if(tj.MCSMom == USHRT_MAX) tj.MCSMom = MCSMom(tjs, tj);
        if(tj.MCSMom < 50) {
          killPts = 0;
          tj.AlgMod[kGottaKink] = false;
          tj.AlgMod[kNoKinkChk] = true;
          if(prt) mf::LogVerbatim("TC")<<"GottaKink turning off kink checking. MCSMom "<<tj.MCSMom;
        }
      } // turn off kink check
    }
    if(prt) mf::LogVerbatim("TC")<<"GottaKink kinkPt "<<kinkPt<<" Pos "<<PrintPos(tjs, tj.Pts[kinkPt])<<" dang "<<dang<<" cut "<<fKinkAngCut<<" kinkSig "<<kinkSig<<" tpFit chi "<<tpFit.FitChi<<" nHiMultPt "<<nHiMultPt<<" killPts "<<killPts<<" GottaKink? "<<tj.AlgMod[kGottaKink];
    
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
    unsigned short ii, prevPtWithHits = USHRT_MAX;
    unsigned short firstFitPt = tj.EndPt[0];
    unsigned short ipt, ndead;
    for(ii = 1; ii < tj.Pts.size(); ++ii) {
      ipt = lastPt - ii;
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
    if(lastPt < 6) minPtsFit = 2;
    if(tj.PDGCode == 13) {
      // Fitting a muon
      maxChi = fMaxChi;
      minPtsFit = lastPt / 3;
    }
    
    // Set the lastPT delta before doing the fit
    lastTP.Delta = PointTrajDOCA(tjs, lastTP.HitPos[0], lastTP.HitPos[1], lastTP);
    
    if(prt) {
      mf::LogVerbatim("TC")<<"UpdateTraj: lastPt "<<lastPt<<" lastTP.Delta "<<lastTP.Delta<<" previous point with hits "<<prevPtWithHits<<" tj.Pts size "<<tj.Pts.size()<<" AngleRange "<<AngleRange(lastTP)<<" PDGCode "<<tj.PDGCode<<" maxChi "<<maxChi<<" minPtsFit "<<minPtsFit;
    }
    
    UpdateAveChg(tj);

    if(lastPt == 1) {
      // Handle the second trajectory point. No error calculation. Just update
      // the position and direction
      lastTP.NTPsFit = 2;
      FitTraj(tj);
      lastTP.FitChi = 0.01;
      lastTP.AngErr = tj.Pts[0].AngErr;
      if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: Second traj point pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1];
      fUpdateTrajOK = true;
      return;
    }
    
    if(lastPt == 2) {
      // Third trajectory point. Keep it simple
      lastTP.NTPsFit = 3;
      FitTraj(tj);
      fUpdateTrajOK = true;
      if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: Third traj point fit "<<lastTP.FitChi;
      return;
    }
    
    // Fit with > 2 TPs
    // Keep adding hits until Chi/DOF exceeds 1
    if(tj.Pts[prevPtWithHits].FitChi < 1) lastTP.NTPsFit += 1;
    // Reduce the number of points fit if the trajectory is long and chisq is getting a bit larger
    if(lastPt > 20 && tj.Pts[prevPtWithHits].FitChi > 1.5 && lastTP.NTPsFit > minPtsFit) lastTP.NTPsFit -= 2;
    // Restrict the number of fitted hits if this is a low momentum trajectory
    if(tj.MCSMom < 100 && lastTP.NTPsFit > fMinPtsFit[tj.Pass]) lastTP.NTPsFit = fMinPtsFit[tj.Pass];
    // also restrict for VLA trajectories
    if(AngleRange(lastTP) == fAngleRanges.size() - 1) lastTP.NTPsFit = fMinPtsFit[tj.Pass];
    
    FitTraj(tj);
    
    // don't get too fancy when we are starting out
    if(lastPt < 6) {
      fUpdateTrajOK = true;
      UpdateDeltaRMS(tj);
      if(prt) mf::LogVerbatim("TC")<<" Return with lastTP.FitChi "<<lastTP.FitChi<<" Chg "<<lastTP.Chg;
      return;
    }
    
    // find the first point that was fit.
    unsigned short cnt = 0;
    for(ii = 0; ii < tj.Pts.size(); ++ii) {
      ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg > 0) {
        firstFitPt = ipt;
        ++cnt;
      }
      if(cnt == lastTP.NTPsFit) break;
      if(ipt == 0) break;
    }
    
    // Maybe maxChi should be replaced with 1.5 here to make this decision independent of the user settings
    if(lastTP.FitChi > 1.5 && tj.Pts.size() > 6) {
      // A large chisq jump can occur if we just jumped a large block of dead wires. In
      // this case we don't want to mask off the last TP but reduce the number of fitted points
      // This count will be off if there a lot of dead or missing wires...
      ndead = 0;
      if(lastTP.FitChi > 3 && firstFitPt < lastPt) ndead = DeadWireCount(lastTP.HitPos[0], tj.Pts[firstFitPt].HitPos[0], fCTP);
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
        if(prevPtWithHits != USHRT_MAX) chirat = lastTP.FitChi / tj.Pts[prevPtWithHits].FitChi;
        fMaskedLastTP = (chirat > 1.5 && lastTP.NTPsFit > 0.3 * NumPtsWithCharge(tj, false));
        if(prt) {
          mf::LogVerbatim("TC")<<" First fit chisq too large "<<lastTP.FitChi<<" prevPtWithHits chisq "<<tj.Pts[prevPtWithHits].FitChi<<" chirat "<<chirat<<" NumPtsWithCharge "<<NumPtsWithCharge(tj, false)<<" fMaskedLastTP "<<fMaskedLastTP;
        }
        // we should also mask off the last TP if there aren't enough hits
        // to satisfy the minPtsFit constraint
        if(!fMaskedLastTP && NumPtsWithCharge(tj, true) < minPtsFit) fMaskedLastTP = true;
      } // few dead wires
    } // lastTP.FitChi > 2 ...
    
    // Deal with a really long trajectory that is in trouble (uB cosmic).
    if(tj.PDGCode == 13 && lastTP.FitChi > fMaxChi) {
      if(lastTP.NTPsFit > 1.3 * fMuonTag[0]) {
        lastTP.NTPsFit *= 0.8;
        if(prt) mf::LogVerbatim("TC")<<" Muon - Reduce NTPsFit "<<lastPt;
      } else {
        fMaskedLastTP = true;
        if(prt) mf::LogVerbatim("TC")<<" Muon - mask last point "<<lastPt;
      }
    }
    
    if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: First fit "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" FitChi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit<<" fMaskedLastTP "<<fMaskedLastTP;
    if(fMaskedLastTP) {
      UnsetUsedHits(lastTP);
      DefineHitPos(lastTP);
      SetEndPoints(tjs, tj);
      lastPt = tj.EndPt[1];
      lastTP.NTPsFit -= 1;
      FitTraj(tj);
      fUpdateTrajOK = true;
      return;
    }  else {
      // a more gradual increase in chisq. Reduce the number of points
      unsigned short newNTPSFit = lastTP.NTPsFit;
      // reduce the number of points fit to keep Chisq/DOF < 2 adhering to the pass constraint
      // and also a minimum number of points fit requirement for long muons
      float prevChi = lastTP.FitChi;
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
        FitTraj(tj);
        if(lastTP.FitChi > prevChi) {
          if(prt) mf::LogVerbatim("TC")<<"  Chisq is increasing "<<lastTP.FitChi<<"  Try to remove an earlier bad hit";
          MaskBadTPs(tj, 1.5);
        }
        prevChi = lastTP.FitChi;
        if(lastTP.NTPsFit == minPtsFit) break;
      } // lastTP.FitChi > 2 && lastTP.NTPsFit > 2
    }
    // last ditch attempt. Drop the last hit
    if(tj.Pts.size() > fMinPtsFit[tj.Pass] && lastTP.FitChi > maxChi) {
      if(prt) mf::LogVerbatim("TC")<<"  Last try. Drop last TP "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;
      UnsetUsedHits(lastTP);
      DefineHitPos(lastTP);
      SetEndPoints(tjs, tj);
      lastPt = tj.EndPt[1];
      FitTraj(tj);
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
  void TrajClusterAlg::FitTraj(Trajectory& tj)
  {
    // Jacket around FitTraj to fit the leading edge of the supplied trajectory
    unsigned short originPt = tj.EndPt[1];
    unsigned short npts = tj.Pts[originPt].NTPsFit;
    TrajPoint tpFit;
    unsigned short fitDir = -1;
    FitTraj(tj, originPt, npts, fitDir, tpFit);
    tj.Pts[originPt] = tpFit;
    
  } // FitTraj

  //////////////////////////////////////////
  void TrajClusterAlg::FitTraj(Trajectory& tj, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint& tpFit)
  {
    // Fit the supplied trajectory using HitPos positions with the origin at originPt.
    // The npts is interpreted as the number of points on each side of the origin
    // The allowed modes are as follows, where i denotes a TP that is included, . denotes
    // a TP with no hits, and x denotes a TP that is not included
    //TP 012345678  fitDir  originPt npts
    //   Oiiixxxxx   1        0       4 << npts in the fit
    //   xi.iiOxxx  -1        5       4
    //   xiiiOiiix   0        4       4 << 2 * npts + 1 points in the fit
    //   xxxiO.ixx   0        4       1
    //   0iiixxxxx   0        0       4
   // This routine puts the results into tp if the fit is successfull. The
    // fit "direction" is in increasing order along the trajectory from 0 to tj.Pts.size() - 1.
    
//    static const float twoPi = 2 * M_PI;
    
    if(originPt > tj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"FitTraj: Requesting fit of invalid TP "<<originPt;
      return;
    }
    
    // copy the origin TP into the fit TP
    tpFit = tj.Pts[originPt];
    // Assume that the fit will fail
    tpFit.FitChi = 999;
    if(fitDir < -1 || fitDir > 1) return;
    
    std::vector<double> x, y;
    std::array<float, 2> origin = tj.Pts[originPt].HitPos;
    // Use TP position if there aren't any hits on it
    if(tj.Pts[originPt].Chg == 0) origin = tj.Pts[originPt].Pos;
    
    // simple two point case
    if(NumPtsWithCharge(tj, false) == 2) {
      for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        double xx = tj.Pts[ipt].HitPos[0] - origin[0];
        double yy = tj.Pts[ipt].HitPos[1] - origin[1];
        x.push_back(xx);
        y.push_back(yy);
      } // ii
      if(x.size() != 2) return;
      if(x[0] == x[1]) {
        // Either + or - pi/2
        tpFit.Ang = M_PI/2;
        if(y[1] < y[0]) tpFit.Ang = -tpFit.Ang;
      } else {
        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        tpFit.Ang = atan2(dy, dx);
      }
      tpFit.Dir[0] = cos(tpFit.Ang);
      tpFit.Dir[1] = sin(tpFit.Ang);
      tpFit.Pos[0] += origin[0];
      tpFit.Pos[1] += origin[1];
      tpFit.AngErr = 0.01;
      tpFit.FitChi = 0.01;
      return;
    } // two points

    std::vector<double> w, q;
    std::array<float, 2> dir;
    double xx, yy, xr, yr;
    double chgWt;

    // Rotate the traj hit position into the coordinate system defined by the
    // originPt traj point, where x = along the trajectory, y = transverse
    double rotAngle = tj.Pts[originPt].Ang;
    double cs = cos(-rotAngle);
    double sn = sin(-rotAngle);

    // enter the originPT hit info if it exists
    if(tj.Pts[originPt].Chg > 0) {
      xx = tj.Pts[originPt].HitPos[0] - origin[0];
      yy = tj.Pts[originPt].HitPos[1] - origin[1];
//      if(prt) std::cout<<" originPT "<<originPt<<" xx "<<xx<<" "<<yy<<" chg "<<tj.Pts[originPt].Chg<<"\n";
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      x.push_back(xr);
      y.push_back(yr);
      chgWt = tj.Pts[originPt].ChgPull;
      if(chgWt < 1) chgWt = 1;
      chgWt *= chgWt;
      w.push_back(chgWt * tj.Pts[originPt].HitPosErr2);
    }
    
    // correct npts to account for the origin point
    if(fitDir != 0) --npts;
    
    // step in the + direction first
    if(fitDir != -1) {
      unsigned short cnt = 0;
      for(unsigned short ipt = originPt + 1; ipt < tj.Pts.size(); ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
//        if(prt) std::cout<<"ipt "<<ipt<<" xx "<<xx<<" yy "<<yy<<" chg "<<tj.Pts[ipt].Chg<<"\n";
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if(chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if(cnt == npts) break;
      } // ipt
    } // fitDir != -1
    
    // step in the - direction next
    if(fitDir != 1 && originPt > 0) {
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = originPt - ii;
        if(ipt > tj.Pts.size() - 1) continue;
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
//        if(prt) std::cout<<"ipt "<<ipt<<" xx "<<xx<<" "<<yy<<" chg "<<tj.Pts[ipt].Chg<<"\n";
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if(chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if(cnt == npts) break;
        if(ipt == 0) break;
      } // ipt
    } // fitDir != -1
    
    // Not enough points to define a line?
    if(x.size() < 2) return;
    
//    if(prt) std::cout<<"FitTraj: npts "<<npts<<" origin "<<origin[0]<<" "<<origin[1]<<" ticks "<<origin[1]/tjs.UnitsPerTick<<" rotAngle "<<rotAngle<<"\n";
    
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;

    // weight by the charge ratio and accumulate sums
    double wght;
    for(unsigned short ipt = 0; ipt < x.size(); ++ipt) {
      if(w[ipt] < 0.00001) w[ipt] = 0.00001;
      wght = 1 / w[ipt];
      sum   += wght;
      sumx  += wght * x[ipt];
      sumy  += wght * y[ipt];
      sumx2 += wght * x[ipt] * x[ipt];
      sumy2 += wght * y[ipt] * y[ipt];
      sumxy += wght * x[ipt] * y[ipt];
    }
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0) return;
    // A is the intercept
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // B is the slope
    double B = (sumxy * sum  - sumx * sumy) / delta;
    
    // The chisq will be set below if there are enough points. Don't allow it to be 0
    // so we can take Chisq ratios later
    tpFit.FitChi = 0.01;
    double newang = atan(B);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAngle);
    sn = sin(rotAngle);
    tpFit.Dir[0] = cs * dir[0] - sn * dir[1];
    tpFit.Dir[1] = sn * dir[0] + cs * dir[1];
    // ensure that the direction is consistent with the originPt direction
    bool flipDir = false;
    if(AngleRange(tj.Pts[originPt]) > 0) {
      flipDir = std::signbit(tpFit.Dir[1]) != std::signbit(tj.Pts[originPt].Dir[1]);
    } else {
      flipDir = std::signbit(tpFit.Dir[0]) != std::signbit(tj.Pts[originPt].Dir[0]);
    }
    if(flipDir) {
      tpFit.Dir[0] = -tpFit.Dir[0];
      tpFit.Dir[1] = -tpFit.Dir[1];
    }
    tpFit.Ang = atan2(tpFit.Dir[1], tpFit.Dir[0]);
//    if(prt) mf::LogVerbatim("TC")<<"FitTraj "<<originPt<<" originPt Dir "<<tj.Pts[originPt].Dir[0]<<" "<<tj.Pts[originPt].Dir[1]<<" rotAngle "<<rotAngle<<" tpFit.Dir "<<tpFit.Dir[0]<<" "<<tpFit.Dir[1]<<" Ang "<<tpFit.Ang<<" flipDir "<<flipDir<<" fit vector size "<<x.size();

    // rotate (0, intcpt) into (W,T) coordinates
    tpFit.Pos[0] = -sn * A + origin[0];
    tpFit.Pos[1] =  cs * A + origin[1];
    // force the origin to be at origin[0]
    MoveTPToWire(tpFit, origin[0]);
    
    if(x.size() < 3) return;
    
    // Calculate chisq/DOF
    double ndof = x.size() - 2;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    if(varnce > 0.) {
      // Intercept error is not used
//      InterceptError = sqrt(varnce * sumx2 / delta);
      double slopeError = sqrt(varnce * sum / delta);
      tpFit.AngErr = std::abs(atan(slopeError));
    } else {
      tpFit.AngErr = 0.01;
    }
    sum = 0;
    // calculate chisq
    double arg;
    for(unsigned short ii = 0; ii < y.size(); ++ii) {
      arg = y[ii] - A - B * x[ii];
      sum += arg * arg / w[ii];
    }
    tpFit.FitChi = sum / ndof;
  
  } // FitTraj
  
  //////////////////////////////////////////
  void TrajClusterAlg::UpdateAveChg(Trajectory& tj)
  {
    if(tj.EndPt[1] == 0) return;
    unsigned short lastPt = tj.EndPt[1];
    tj.AveChg = 0;
    tj.Pts[lastPt].AveChg = 0;

    // calculate ave charge and charge RMS using ALL hits in the trajectory
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
      // RMS until it is well known. Start with 100% error on the
      // charge RMS
      float defFrac = 1 / (float)(tj.EndPt[1]);
      tj.ChgRMS = defFrac + (1 - defFrac) * rms;
      // don't let it get crazy small
      if(tj.ChgRMS < 0.15) tj.ChgRMS = 0.15;
      tj.Pts[lastPt].ChgPull = (tj.Pts[lastPt].Chg / tj.AveChg - 1) / tj.ChgRMS;
    } // cnt > 3

  } // UpdateAveChg

  ////////////////////////////////////////////////
  bool TrajClusterAlg::StartTraj(Trajectory& tj, const unsigned int& fromHit, const unsigned int& toHit, const unsigned short& pass)
  {
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
    // Start a simple (seed) trajectory going from a hit to a position (toWire, toTick).
    
    // decrement the work ID so we can use it for debugging problems
    --fWorkID;
    if(fWorkID == SHRT_MIN) fWorkID = -1;
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
    MakeBareTrajPoint(tjs, fromWire, fromTick, toWire, toTick, tCTP, tp);
    if(tp.Pos[0] < 0) return false;

    tp.AngErr = 0.1;
    if(tj.ID == debug.WorkID) { prt = true; didPrt = true; debug.Plane = fPlane; TJPrt = tj.ID; }
    if(prt) mf::LogVerbatim("TC")<<"StartTraj "<<(int)fromWire<<":"<<(int)fromTick<<" -> "<<(int)toWire<<":"<<(int)toTick<<" dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" ang "<<tp.Ang<<" AngleRange "<<AngleRange(tp)<<" angErr "<<tp.AngErr<<" ExpectedHitsRMS "<<ExpectedHitsRMS(tp);
    tj.Pts.push_back(tp);
    return true;
    
  } // StartTraj
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ChkInTraj(std::string someText)
  {
    // Check tjs.allTraj -> InTraj associations
    
    if(!fUseAlg[kChkInTraj]) return;
    
    unsigned short tID;
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
        myprt<<someText<<" ChkInTraj: inTraj - UseHit mis-match for tj ID "<<tID<<" tj.WorkID "<<tj.WorkID<<" atHits size "<<atHits.size()<<" tHits size "<<tHits.size()<<" in CTP "<<tj.CTP<<"\n";
        myprt<<"AlgMods: ";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
        myprt<<"     inTraj     UseHit \n";
        for(iht = 0; iht < atHits.size(); ++iht) {
          myprt<<"iht "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]])<<"_"<<tjs.fHits[atHits[iht]].InTraj;
          if(iht < tHits.size()) myprt<<" "<<PrintHit(tjs.fHits[tHits[iht]])<<"_"<<tjs.fHits[tHits[iht]].InTraj;
          if(atHits[iht] != tHits[iht]) myprt<<" <<< ";
          myprt<<"\n";
          fQuitAlg = true;
        } // iht
        if(tHits.size() > atHits.size()) {
          for(iht = atHits.size(); iht < atHits.size(); ++iht) {
            myprt<<"atHits "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]])<<"\n";
          } // iht
//          PrintAllTraj(tjs, debug, USHRT_MAX, 0);
        } // tHit.size > atHits.size()
      }
      ++itj;
      if(fQuitAlg) return;
    } // tj
    
  } // ChkInTraj

  ////////////////////////////////////////////////
  void TrajClusterAlg::StoreTraj(Trajectory& tj)
  {

    if(tj.EndPt[1] <= tj.EndPt[0]) return;
    if(tj.AlgMod[kKilled]) {
      mf::LogWarning("TC")<<"StoreTraj: Trying to store a killed trajectory. tj ID "<<tj.ID;
      return;
    }

    if(!(tj.StepDir == 1 || tj.StepDir == -1)) {
      mf::LogError("TC")<<"StoreTraj: Invalid StepDir "<<tj.StepDir;
      fQuitAlg = true;
      return;
    }
    // put trajectories in order of US -> DS
    if(tj.StepDir < 0) ReverseTraj(tjs, tj);
    // This shouldn't be necessary but do it anyway
    SetEndPoints(tjs, tj);
    
    // Calculate the charge near the end and beginning if necessary. This must be a short
    // trajectory. Find the average using 4 points
    if(tj.Pts[tj.EndPt[0]].AveChg <= 0) {
      unsigned short cnt = 0;
      float sum = 0;
      for(unsigned short ipt = tj.EndPt[0] + 1; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == 4) break;
       }
      tj.Pts[tj.EndPt[0]].AveChg = sum / (float)cnt;
    }
    if(tj.Pts[tj.EndPt[1]].AveChg <= 0) {
      float sum = 0;
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.EndPt[1] - ii;
        if(tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == 4) break;
        if(ipt == 0) break;
      } // ii
      tj.Pts[tj.EndPt[1]].AveChg = sum / (float)cnt;
    } // begin charge == end charge
    
    CalculateQuality(tj);
    
    short trID = tjs.allTraj.size() + 1;
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(tj.Pts[ipt].UseHit[ii]) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(tjs.fHits[iht].InTraj > 0) {
            mf::LogWarning("TC")<<"StoreTraj: Failed trying to store hit "<<PrintHit(tjs.fHits[iht])<<" in new tjs.allTraj "<<trID<<" but it is used in traj ID = "<<tjs.fHits[iht].InTraj<<" print and quit";
            PrintTrajectory("SW", tjs, tj, USHRT_MAX);
            ReleaseHits(tj);
            fQuitAlg = true;
            return;
          } // error
          tjs.fHits[iht].InTraj = trID;
        }
      } // ii
    } // ipt
    
    // ensure that inTraj is clean for the ID
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].InTraj == tj.ID) {
        mf::LogWarning("TC")<<"StoreTraj: Hit "<<PrintHit(tjs.fHits[iht])<<" thinks it belongs to traj ID "<<tj.ID<<" but it wasn't stored\n";
        PrintTrajectory("SW", tjs, tj, USHRT_MAX);
        fQuitAlg = true;
        return;
      }
    } // iht

    tj.WorkID = tj.ID;
    tj.ID = trID;
    tjs.allTraj.push_back(tj);
    if(prt) mf::LogVerbatim("TC")<<"StoreTraj trID "<<trID<<" CTP "<<tj.CTP<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
    if(debug.Hit != UINT_MAX) {
      // print out some debug info
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
        for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(iht == debug.Hit) std::cout<<"Debug hit appears in trajectory w WorkID "<<tj.WorkID<<" UseHit "<<tj.Pts[ipt].UseHit[ii]<<"\n";
        } // ii
      } // ipt
    } // debug.Hit ...
    ChkInTraj("StoreTraj");
    
  } // StoreTraj
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CalculateQuality(Trajectory& tj)
  {
    // Calculate a quality metric using the deviations of all hits used in the trajectory
    
    tj.MCSMom = MCSMom(tjs, tj);
    
  } // CalculateQuality
  
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
    
    ChkInTraj("MATC");
    if(fQuitAlg) return;
    
    unsigned short itj, endPt0, endPt1, ii;
    
    // Make one cluster for each trajectory. The indexing of trajectory parents
    // should map directly to cluster parents
    short clID = 0;
    for(itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      if(tj.StepDir > 0) ReverseTraj(tjs, tj);
      // ensure that the endPts are correct
      SetEndPoints(tjs, tj);
      // some sort of error occurred
      if(tj.EndPt[0] >= tj.EndPt[1]) {
        mf::LogWarning("TC")<<"MakeAllTrajClusters failed in SetEndPoints "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        continue;
      }
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
      auto tHits = PutTrajHitsInVector(tj, kUsedHits);
      if(tHits.empty()) {
        mf::LogWarning("TC")<<"MakeAllTrajClusters: No hits found in trajectory "<<itj<<" so skip it";
        continue;
      } // error
      cls.tclhits = tHits;
      // Set the traj info
      tj.ClusterIndex = tjs.tcl.size();
      tjs.tcl.push_back(cls);
      // do some checking and define tjs.inClus
      geo::PlaneID planeID = DecodeCTP(cls.CTP);
      for(ii = 0; ii < cls.tclhits.size(); ++ii) {
        iht = cls.tclhits[ii];
        if(tjs.fHits[iht].WireID.Plane != planeID.Plane ||
           tjs.fHits[iht].WireID.Cryostat != planeID.Cryostat ||
           tjs.fHits[iht].WireID.TPC != planeID.TPC) {
          mf::LogWarning("TC")<<"MakeAllTrajClusters: Bad OLD hit CTP in itj "<<itj<<" hit "<<PrintHit(tjs.fHits[iht])<<" WorkID "<<tjs.allTraj[itj].WorkID<<" Plane "<<tjs.fHits[iht].WireID.Plane<<" vs "<<planeID.Plane<<" Cstat "<<tjs.fHits[iht].WireID.Cryostat<<" vs "<<planeID.Cryostat<<" TPC "<<tjs.fHits[iht].WireID.TPC<<" vs "<<planeID.TPC;
          fQuitAlg = true;
          return;
        }
        if(tjs.inClus[iht] != 0) {
          mf::LogWarning("TC")<<"MakeAllTrajClusters: Trying to assign tj.ID "<<tj.ID<<" hit "<<iht<<"_"<<PrintHit(tjs.fHits[iht])<<" to already-assigned cluster "<<tjs.inClus[iht];
          fQuitAlg = true;
          return;
        }
        tjs.inClus[iht] = clID;
      } //iht
    } // itj

  } // MakeAllTrajClusters

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
      // theHit is narrow and it is the highest amplitude hit in the multiplet. Ignore any
      // others that are short and fat
      auto tmp = hitsInMultiplet;
      tmp.resize(1);
      tmp[0] = theHit;
      float shortCut = 0.6 * maxPeak;
      for(auto& iht : hitsInMultiplet) {
        // reject short fat hits
        if(tjs.fHits[iht].PeakAmplitude < shortCut && tjs.fHits[iht].RMS > narrowHitCut) continue;
        tmp.push_back(iht);
      } // iht
      hitsInMultiplet = tmp;
    } else {
      // theHit is not narrow and it is not the tallest. Ignore a single hit if it is
      // the tallest and narrow
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
  bool TrajClusterAlg::TrajHitsOK(const std::vector<unsigned int>& iHitsInMultiplet, const std::vector<unsigned int>& jHitsInMultiplet)
  {
    // Hits (assume to be on adjacent wires have an acceptable signal overlap
    
    if(iHitsInMultiplet.empty() || jHitsInMultiplet.empty()) return false;
    
    float sum;
    float cvI = HitsPosTick(tjs, iHitsInMultiplet, sum, kAllHits);
    float minI = 1E6;
    float maxI = 0;
    for(auto& iht : iHitsInMultiplet) {
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - 3 * rms;
      if(arg < minI) minI = arg;
      arg = cv + 3 * rms;
      if(arg > maxI) maxI = arg;
    }
    
    float cvJ = HitsPosTick(tjs, jHitsInMultiplet, sum, kAllHits);
    float minJ = 1E6;
    float maxJ = 0;
    for(auto& jht : jHitsInMultiplet) {
      float cv = tjs.fHits[jht].PeakTime;
      float rms = tjs.fHits[jht].RMS;
      float arg = cv - 3 * rms;
      if(arg < minJ) minJ = arg;
      arg = cv + 3 * rms;
      if(arg > maxJ) maxJ = arg;
    }
    
    if(cvI < cvJ) {
      if(maxI > minJ) return true;
    } else {
      if(minI < maxJ) return true;
    }
    return false;
  } // TrajHitsOK
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FillWireHitRange(geo::TPCID const& tpcid)
  {
    // fills the WireHitRange vector. Slightly modified version of the one in ClusterCrawlerAlg.
    // Also fills the WirePtr vector
    
    // determine the number of planes
    art::ServiceHandle<geo::Geometry> geom;
    geo::TPCGeo const& TPC = geom->TPC(tpcid);
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    unsigned short nplanes = TPC.Nplanes();
    tjs.NumPlanes = nplanes;
    
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
    
    raw::ChannelID_t channel = geom->PlaneWireToChannel(0, 0, (int)tpc, (int)cstat);
    float wirePitch = geom->WirePitch(geom->View(channel));
    float tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    tjs.UnitsPerTick = tickToDist / wirePitch;
/* do this later after we make sure that things work
    // convert hit ticks to time if it hasn't already been done
    if(!tjs.ConvertTicksToTime) {
      for(auto& hit : tjs.fHits) hit.PeakTime *= tjs.UnitsPerTick;
      tjs.ConvertTicksToTime = false;
    }
*/
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      tjs.FirstWire[ipl] = INT_MAX;
      tjs.LastWire[ipl] = 0;
      tjs.NumWires[ipl] = geom->Nwires(ipl, tpc, cstat);
      tjs.WireHitRange[ipl].resize(tjs.NumWires[ipl], flag);
      tjs.MaxPos0[ipl] = (float)(tjs.NumWires[ipl] - 0.5);
      tjs.MaxPos1[ipl] = (float)detprop->NumberTimeSamples() * tjs.UnitsPerTick;
    }
    
    // overwrite with the "dead wires" condition
    flag.first = -1; flag.second = -1;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        raw::ChannelID_t chan = geom->PlaneWireToChannel((int)ipl, (int)wire, (int)tpc, (int)cstat);
        if(!channelStatus.IsGood(chan)) tjs.WireHitRange[ipl][wire] = flag;
//        if(!channelStatus.IsGood(chan)) std::cout<<"chan "<<chan<<" is not good "<<ipl<<":"<<wire<<"\n";
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
/*
      if(tjs.fHits[iht].StartTick > (int)detprop->NumberTimeSamples() || tjs.fHits[iht].EndTick > (int)detprop->NumberTimeSamples()) {
        std::cout<<"Bad StartTick "<<tjs.fHits[iht].StartTick<<" or EndTick "<<tjs.fHits[iht].EndTick<<" NumberTimeSamples "<<detprop->NumberTimeSamples()<<"\n";
      }
*/
    } // iht
    
    if(!CheckWireHitRange()) {
      fQuitAlg = true;
      return;
    }
    
    // Find the average multiplicity 1 hit RMS and calculate the expected max RMS for each range
//    std::cout<<"pln  AngleRange  AngleRangeMaxHitsRMS\n";
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      float sumRMS = 0;
      unsigned int cnt = 0;
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].Multiplicity != 1) continue;
          ++cnt;
          sumRMS += tjs.fHits[iht].RMS;
        } // iht
      } // wire
      if(cnt < 4) continue;
      fAveHitRMS[ipl] = sumRMS/(float)cnt;
      // calculate the max RMS expected for each angle range
      for(unsigned short ii = 0; ii < fAngleRanges.size(); ++ii) {
        float angle = fAngleRanges[ii];
        if(angle < M_PI/2) {
          fAngleRangesMaxHitsRMS[ii] = 1.5 * fAveHitRMS[ipl] + tan(angle) / tjs.UnitsPerTick;
        } else {
          fAngleRangesMaxHitsRMS[ii] = 1000;
        }
//        std::cout<<"Pln "<<ipl<<" MaxAngle "<<(int)(angle*180/M_PI)<<" "<<std::fixed<<std::setprecision(1)<<fAngleRangesMaxHitsRMS[ii]<<"\n";
      } // ii
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
      raw::ChannelID_t channel = geom->PlaneWireToChannel((int)tcHit.WireID.Plane, (int)tcHit.WireID.Wire, (int)tcHit.WireID.TPC, (int)tcHit.WireID.Cryostat);
      tmp.emplace_back(channel,
                       tcHit.StartTick, tcHit.EndTick,
                       tcHit.PeakTime, 0,
                       tcHit.RMS,
                       tcHit.PeakAmplitude, 0,
                       tcHit.Integral, tcHit.Integral, 0,
                       tcHit.Multiplicity, tcHit.LocalIndex,
                       tcHit.GoodnessOfFit, tcHit.NDOF,
                       geom->View(channel),
                       geom->SignalType(planeID),
                       tcHit.WireID
                       );
    } // tcHit
     return tmp;
  } // YieldHits
  
  ////////////////////////////////////////////////
  float TrajClusterAlg::ExpectedHitsRMS(TrajPoint const& tp)
  {
    // returns the expected RMS of hits for the trajectory point
    if(std::abs(tp.Dir[0]) > 0.001) {
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      return 1.5 * fAveHitRMS[planeID.Plane] + std::abs(tp.Dir[1]/tp.Dir[0])/tjs.UnitsPerTick;
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
//    std::cout<<"EraseHit "<<delHit<<"  "<<PrintHit(tjs.fHits[delHit])<<"\n";
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
//    hitCTP.Channel = geom->PlaneWireToChannel((int)planeID.Plane,(int)hitWire,(int)planeID.TPC,(int)planeID.Cryostat);
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
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(NumUsedHits(tp) < 2) continue;
        // Make a list of the old hits on this TP before doing anything invasive
        std::vector<unsigned int> oldHits;
        // get some info so we can calculate the RMS
        raw::TDCtick_t loTick = INT_MAX;
        raw::TDCtick_t hiTick = 0;
        float mChg = 0;
        float mTick = 0;
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          oldHits.push_back(iht);
          if(tjs.fHits[iht].StartTick < loTick) loTick = tjs.fHits[iht].StartTick;
          if(tjs.fHits[iht].EndTick > hiTick) hiTick = tjs.fHits[iht].EndTick;
          mChg += tjs.fHits[iht].Integral;
          mTick += tjs.fHits[iht].Integral * tjs.fHits[iht].PeakTime;
        } // ii
        mTick /= mChg;
        if(mTick < 0) mTick = 0;
        // make a temporary signal waveform vector
        std::vector<float> signal(hiTick - loTick, 0);
        // fill it with the hit shapes
        for(auto& iht : oldHits) {
          float peakTime = tjs.fHits[iht].PeakTime;
          float amp = tjs.fHits[iht].PeakAmplitude;
          float rms = tjs.fHits[iht].RMS;
//          std::cout<<"tj.ID "<<tj.ID<<" oldHit "<<iht<<" "<<PrintHit(tjs.fHits[iht])<<" "<<(int)tjs.fHits[iht].Integral<<"\n";
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
//        std::cout<<"mHit tick "<<mTick<<" mRMS "<<mRMS<<" mChg "<<(int)mChg<<"\n";
        // Modify the first hit in the list
        unsigned int mht = oldHits[0];
        tjs.fHits[mht].PeakTime = mTick;
        tjs.fHits[mht].PeakAmplitude = mChg / (2.5066 * mRMS);
        tjs.fHits[mht].Integral = mChg;
        tjs.fHits[mht].RMS = mRMS;
        tjs.fHits[mht].Multiplicity = 1;
        tjs.fHits[mht].LocalIndex = 0;
        tjs.fHits[mht].GoodnessOfFit = 1; // flag?
        tjs.fHits[mht].NDOF = 0;
        // then flag the other hits for erasing
        for(unsigned short ii = 1; ii < oldHits.size(); ++ii) {
          tp.UseHit[ii] = false;
          // put it in the removal list
          delHits.push_back(tp.Hits[ii]);
          // Flag this hit
          tp.Hits[ii] = INT_MAX;
          tjs.fHits[ii].InTraj = SHRT_MAX;
        } // ii
      } // ipt
    } // itj
    
    // Erase the hits. Start by sorting them in decreasing order so that
    // the local delHits vector doesn't need to be modified when a hit is deleted
    if(delHits.size() > 1) std::sort(delHits.begin(), delHits.end(), std::greater<unsigned int>());

    for(auto& delHit : delHits) EraseHit(delHit);

  } // MergeTPHits
  
  /////////////////////////////////////////
  unsigned short TrajClusterAlg::AngleRange(TrajPoint const& tp)
  {
    // returns the index of the angle range
    float dang = tp.Ang;
    if(dang > M_PI) dang = M_PI;
    if(dang < -M_PI) dang = M_PI;
    if(dang < 0) dang = -dang;
    if(dang > M_PI/2) dang = M_PI - dang;
    for(unsigned short ir = 0; ir < fAngleRanges.size(); ++ir) {
      if(dang < fAngleRanges[ir]) return ir;
    }
    return fAngleRanges.size() - 1;
  } // AngleRange

  /////////////////////////////////////////
  bool TrajClusterAlg::SignalAtTp(TrajPoint const& tp)
  {
    return SignalAtPos(tp.Pos[0], tp.Pos[1], tp.CTP);
  } // SignalAtTp
  
  
  /////////////////////////////////////////
  bool TrajClusterAlg::SignalAtPos(float pos0, float pos1, CTP_t tCTP)
  {
    // Returns true if the TP is within the TPC
    
    if(pos0 < 0) return false;
    if(pos1 < 0) return false;
    unsigned int wire = std::nearbyint(pos0);
    geo::PlaneID planeID = DecodeCTP(tCTP);
    unsigned int ipl = planeID.Plane;
    if(wire >= tjs.NumWires[ipl]) return false;
    if(pos1 > tjs.MaxPos1[ipl]) return false;
    // Assume dead wires have a signal
    if(tjs.WireHitRange[ipl][wire].first == -1) return true;
    raw::TDCtick_t rawProjTick = (float)(pos1 / tjs.UnitsPerTick);
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      if(rawProjTick > tjs.fHits[iht].StartTick && rawProjTick < tjs.fHits[iht].EndTick) return true;
    } // iht
    return false;
  } // SignalAtPos
  
  //////////////////////////////////////////
  unsigned short TrajClusterAlg::NumPtsWithCharge(Trajectory& tj, bool includeDeadWires)
  {
    unsigned short ntp = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]+1; ++ipt) if(tj.Pts[ipt].Chg > 0) ++ntp;
    // Add the count of deadwires
    unsigned short ipt0 = tj.EndPt[0];
    unsigned short ipt1 = tj.EndPt[1];
    if(ipt1 > tj.Pts.size() - 1) return 0;
    if(includeDeadWires) ntp += DeadWireCount(tj.Pts[ipt0], tj.Pts[ipt1]);
    return ntp;
  } // NumTPsWithCharge
  
  //////////////////////////////////////////
  unsigned short TrajClusterAlg::NumUsedHits(TrajPoint& tp)
  {
    // Counts the number of used hits in tp
    unsigned short nused = 0;
    if(tp.Hits.empty()) return nused;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if(tp.UseHit[ii]) ++nused;
    return nused;
  } // NumUsedHits
  
  //////////////////////////////////////////
  float TrajClusterAlg::DeadWireCount(TrajPoint& tp1, TrajPoint& tp2)
  {
    return DeadWireCount(tp1.Pos[0], tp2.Pos[0], tp1.CTP);
  } // DeadWireCount
  
  //////////////////////////////////////////
  float TrajClusterAlg::DeadWireCount(float inWirePos1, float inWirePos2, CTP_t tCTP)
  {
    if(inWirePos1 < -0.4 || inWirePos2 < -0.4) return 0;
    unsigned int inWire1 = std::nearbyint(inWirePos1);
    unsigned int inWire2 = std::nearbyint(inWirePos2);
    geo::PlaneID planeID = DecodeCTP(tCTP);
    unsigned short plane = planeID.Plane;
    if(inWire1 > tjs.NumWires[plane] || inWire2 > tjs.NumWires[plane]) return 0;
    if(inWire1 > inWire2) {
      // put in increasing order
      unsigned int tmp = inWire1;
      inWire1 = inWire2;
      inWire2 = tmp;
    } // inWire1 > inWire2
    ++inWire2;
    unsigned int wire, ndead = 0;
    for(wire = inWire1; wire < inWire2; ++wire) if(tjs.WireHitRange[plane][wire].first == -1) ++ndead;
    return ndead;
  } // DeadWireCount

  //////////////////////////////////////////
  void TrajClusterAlg::CheckHitClusterAssociations()
  {
    // check hit - cluster associations
    
    if(tjs.fHits.size() != tjs.inClus.size()) {
      mf::LogWarning("TC")<<"CHCA: Sizes wrong "<<tjs.fHits.size()<<" "<<tjs.inClus.size();
      fQuitAlg = true;
      return;
    }
    
    unsigned int iht;
    short clID;
    
    // check cluster -> hit association
    for(unsigned short icl = 0; icl < tjs.tcl.size(); ++icl) {
      if(tjs.tcl[icl].ID < 0) continue;
      clID = tjs.tcl[icl].ID;
      for(unsigned short ii = 0; ii < tjs.tcl[icl].tclhits.size(); ++ii) {
        iht = tjs.tcl[icl].tclhits[ii];
        if(iht > tjs.fHits.size() - 1) {
          mf::LogWarning("CC")<<"CHCA: Bad tclhits index "<<iht<<" tjs.fHits size "<<tjs.fHits.size();
          fQuitAlg = true;
          return;
        } // iht > tjs.fHits.size() - 1
        if(tjs.inClus[iht] != clID) {
          mf::LogError("TC")<<"CHCA: Bad cluster -> hit association. clID "<<clID<<" hit "<<PrintHit(tjs.fHits[iht])<<" tjs.inClus "<<tjs.inClus[iht]<<" CTP "<<tjs.tcl[icl].CTP;
          FindHit("CHCA ", iht);
          fQuitAlg = true;
          return;
        }
      } // ii
    } // icl
    
    // check hit -> cluster association
    unsigned short icl;
    for(iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.inClus[iht] <= 0) continue;
      icl = tjs.inClus[iht] - 1;
      // see if the cluster is obsolete
      if(tjs.tcl[icl].ID < 0) {
        mf::LogError("TC")<<"CHCA: Hit "<<PrintHit(tjs.fHits[iht])<<" associated with an obsolete cluster tjs.tcl[icl].ID "<<tjs.tcl[icl].ID;
        fQuitAlg = true;
        return;
      }
      if (std::find(tjs.tcl[icl].tclhits.begin(), tjs.tcl[icl].tclhits.end(), iht) == tjs.tcl[icl].tclhits.end()) {
        mf::LogError("TC")<<"CHCA: Hit "<<tjs.tcl[icl].CTP<<":"<<PrintHit(tjs.fHits[iht])<<" -> tjs.inClus "<<tjs.inClus[iht]<<" but isn't in tjs.tcl[icl].ID "<<tjs.tcl[icl].ID<<" list of hits. icl "<<icl<<" iht "<<iht;
        for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
          if(tjs.allTraj[itj].ClusterIndex == icl) mf::LogError("TC")<<"CHCA: Cluster index "<<icl<<" found in traj ID "<<tjs.allTraj[itj].ID;
        } // itj
        PrintAllTraj("CHCA", tjs, debug, USHRT_MAX, USHRT_MAX);
        fQuitAlg = true;
        return;
      }
    } // iht
    
  } // CheckHitClusterAssociations()
  
  //////////////////////////////////////////
  void TrajClusterAlg::MaskBadTPs(Trajectory& tj, float const& maxChi)
  {
    // Remove TPs that have the worst values of delta until the fit chisq < maxChi
    
    if(!fUseAlg[kMaskBadTPs]) return;
    
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
      UnsetUsedHits(tj.Pts[imBad]);
      FitTraj(tj);
      if(prt) mf::LogVerbatim("TC")<<"  after FitTraj "<<lastTP.FitChi;
      tj.AlgMod[kMaskBadTPs] = true;
      ++nit;
    } // lastTP.FItChi > maxChi && nit < 3
    
  } // MaskBadTPs
  
  //////////////////////////////////////////
  void TrajClusterAlg::MaskTrajEndPoints(Trajectory& tj, unsigned short nPts)
  {
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
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.EndPt[1] - nPts - ii;
      if(tj.Pts[ipt].Chg > 0) {
        lastGoodPt = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(prt) {
      mf::LogVerbatim("TC")<<"MTEP: lastGoodPt "<<lastGoodPt<<" Pts size "<<tj.Pts.size()<<" fGoodTraj "<<fGoodTraj;
    }
    if(lastGoodPt == USHRT_MAX) return;
    tj.EndPt[1] = lastGoodPt;
    bool isVLA = (AngleRange(tj.Pts[lastGoodPt]) == fAngleRanges.size() - 1);
    
    for(unsigned short ii = 0; ii < nPts; ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      UnsetUsedHits(tj.Pts[ipt]);
      // Reset the position and direction of the masked off points
      tj.Pts[ipt].Dir = tj.Pts[lastGoodPt].Dir;
      if(isVLA) {
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
    // Sets the StopsAtEnd bits on the trajectory by identifying the Bragg peak
    // at each end. This is done by inspecting points near both ends for long trajectories
    // and all points for short trajectories and fitting them to an exponential charge pattern hypothesis
    
    tj.StopsAtEnd[0] = false;
    tj.StopsAtEnd[1] = false;
    
    if(!fUseAlg[kChkStop]) return;
    if(fChkStopCuts[0] < 0) return;
    
    // ignore trajectories that are very large angle at both ends
    bool isVLA0 = (AngleRange(tj.Pts[tj.EndPt[0]]) == fAngleRanges.size() - 1);
    bool isVLA1 = (AngleRange(tj.Pts[tj.EndPt[1]]) == fAngleRanges.size() - 1);
    if(isVLA0 && isVLA1) return;
    
    // find the minimum and maximum charged points
    float minChg = 1E6;
    float maxChg = 0;
    unsigned short PtWithMaxChg = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      if(tj.Pts[ipt].Chg > maxChg) {
        maxChg = tj.Pts[ipt].Chg;
        PtWithMaxChg = ipt;
      }
      if(tj.Pts[ipt].Chg < minChg) minChg = tj.Pts[ipt].Chg;
    } // ipt
    // no sense continuing if there is little variation
    if(prt) mf::LogVerbatim("TC")<<"ChkStop: Min Chg "<<(int)minChg<<" maxChg "<<(int)maxChg<<" at point "<<PtWithMaxChg<<" Pos "<<PrintPos(tjs, tj.Pts[PtWithMaxChg])<<" charge ratio "<<maxChg / minChg;
    if(maxChg / minChg < fChkStopCuts[0]) return;
    
    // Try to fit the points whose charge > minChgCut to an exponential
    float minChgCut = fChkStopCuts[0] * minChg;
    
    // check short trajectories using all points.
    if(tj.Pts.size() < 7) {
      // determine which end is most likely the stopping point
      unsigned short end = 0;
      short dir = 1;
      if(abs(PtWithMaxChg - tj.EndPt[1]) < abs(PtWithMaxChg - tj.EndPt[0])) {
        end = 1;
        dir = -1;
      } // end check
      // Fit the charge pattern to an exponential hypothesis
      std::vector<float> x, y, yerr2;
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = PtWithMaxChg + ii * dir;
        if(tj.Pts[ipt].Chg > 0) {
          // See if we are below the minimum charge cut
          if(tj.Pts[ipt].Chg < minChgCut) break;
          float rr = log(0.5 + TrajPointSeparation(tj.Pts[PtWithMaxChg], tj.Pts[ipt]));
          x.push_back(rr);
          y.push_back(tj.Pts[ipt].Chg);
          // Expected charge error is ~10%
          float err = 0.1 * tj.Pts[ipt].Chg;
          yerr2.push_back(err * err);
          if(prt) mf::LogVerbatim("TC")<<"ChkStop "<<PrintPos(tjs, tj.Pts[ipt])<<" rr "<<rr<<" Chg "<<tj.Pts[ipt].Chg;
        } //  Chg > 0
        // hit the end of the trajectory?
        if(ipt == tj.EndPt[1 - end]) break;
      } // ii
      if(x.size() < 3) return;
      float intcpt, slope, intcpterr, slopeerr, chidof;
      fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
      if(prt) mf::LogVerbatim("TC")<<"ChkStop: Short Trajectory Fit intcpt "<<intcpt<<" slope "<<slope<<" err "<<slopeerr<<" chidof "<<chidof;
      // should be negative slope
      if(slope < -fChkStopCuts[1] * slopeerr && chidof < fChkStopCuts[2]) {
        // looks like it stops
        tj.StopsAtEnd[end] = true;
        tj.AlgMod[kChkStop] = true;
        if(prt) mf::LogVerbatim("TC")<<" Stops at end "<<end;
        // mask off the points before/after atPtEnd
        if(PtWithMaxChg != tj.EndPt[end]) {
          for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
            unsigned short ipt = PtWithMaxChg + ii * dir;
            UnsetUsedHits(tj.Pts[ipt]);
            if(ipt == tj.EndPt[1 - end]) break;
          } // ii
          SetEndPoints(tjs, tj);
        } // not at EndPt
      }
      return;
    } // short trajectory
    
    // Long trajectory >= 7 points. Check both ends
    for(unsigned short end = 0; end < 2; ++end) {
      short dir = 1;
      if(end == 1) dir = -1;
      unsigned short nPtsToCheck = tj.Pts.size() / 2;
      float maxChg = 0;
      float PtWithMaxChg = 0;
      unsigned short cnt = 0;
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.EndPt[end] + ii * dir;
        if(tj.Pts[ipt].Chg > maxChg) {
          maxChg = tj.Pts[ipt].Chg;
          PtWithMaxChg = ipt;
          ++cnt;
        }
        if(cnt == nPtsToCheck) break;
        if(ipt == tj.EndPt[1 - end]) break;
      } // ii
      // Fit the charge pattern to an exponential hypothesis
      std::vector<float> x, y, yerr2;
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = PtWithMaxChg + ii * dir;
        if(tj.Pts[ipt].Chg > 0) {
          // See if we are below the minimum charge cut
          if(tj.Pts[ipt].Chg < minChgCut) break;
          float rr = log(0.5 + TrajPointSeparation(tj.Pts[PtWithMaxChg], tj.Pts[ipt]));
          x.push_back(rr);
          y.push_back(tj.Pts[ipt].Chg);
          // Expected charge error is ~10%
          float err = 0.1 * tj.Pts[ipt].Chg;
          yerr2.push_back(err * err);
          if(prt) mf::LogVerbatim("TC")<<"ChkStop "<<PrintPos(tjs, tj.Pts[ipt])<<" rr "<<rr<<" Chg "<<tj.Pts[ipt].Chg;
        }
        if(x.size() == nPtsToCheck) break;
        if(ipt == tj.EndPt[1 - end]) break;
      } // ii
      if(x.size() < 4) continue;
      float intcpt, slope, intcpterr, slopeerr, chidof;
      fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
      if(prt) mf::LogVerbatim("TC")<<"ChkStop: Long trajectory end "<<end<<" Fit intcpt "<<intcpt<<" slope "<<slope<<" err "<<slopeerr<<" chidof "<<chidof<<" nPts fit "<<x.size()<<" PtWithMaxChg "<<PrintPos(tjs, tj.Pts[PtWithMaxChg]);
      // should be negative slope
      if(slope < -fChkStopCuts[1] * slopeerr && chidof < fChkStopCuts[2]) {
        // looks like it stops at atPtEnd1
        tj.StopsAtEnd[end] = true;
        tj.AlgMod[kChkStop] = true;
        if(prt) mf::LogVerbatim("TC")<<" Stops at end "<<end;
        if(PtWithMaxChg != tj.EndPt[end]) {
          for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
            unsigned short ipt = tj.EndPt[end] + ii * dir;
            if(ipt == PtWithMaxChg) break;
            UnsetUsedHits(tj.Pts[ipt]);
          } // ii
          SetEndPoints(tjs, tj);
          if(prt) mf::LogVerbatim("TC")<<" new end points "<<PrintPos(tjs, tj.Pts[tj.EndPt[0]])<<" "<<PrintPos(tjs, tj.Pts[tj.EndPt[1]]);
        } // not at EndPt
      }
    } // end
    
  } // StopsAtEnd

} // namespace cluster
