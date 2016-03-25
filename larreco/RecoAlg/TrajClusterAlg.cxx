//////////////////////////////////////////////////////////////////////
///
/// Step crawling code used by TrajClusterAlg
///
/// Bruce Baller, baller@fnal.gov
///
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrajClusterAlg.h"


// TEMP for TagAllTraj
#include "larsim/MCCheater/BackTracker.h"

class TH1F;
class TH2F;

struct SortEntry{
  unsigned int index;
  float length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}

namespace cluster {
  
  //------------------------------------------------------------------------------
  TrajClusterAlg::TrajClusterAlg(fhicl::ParameterSet const& pset)
  {
    reconfigure(pset);
    
    // define some histograms
    
    art::ServiceHandle<art::TFileService> tfs;
    if(fStudyMode) {
      
      fnHitsPerTP[0] = tfs->make<TH1F>("nHitsPerTP0","nHits / TP Pln 0", 10, 0, 10);
      fnHitsPerTP[1] = tfs->make<TH1F>("nHitsPerTP1","nHits / TP Pln 1", 10, 0, 10);
      fnHitsPerTP[2] = tfs->make<TH1F>("nHitsPerTP2","nHits / TP Pln 2", 10, 0, 10);
      
      fnHitsFitPerTP[0] = tfs->make<TH1F>("nHitsFitPerTP0","nHits fit / TP Pln 0", 10, 0, 10);
      fnHitsFitPerTP[1] = tfs->make<TH1F>("nHitsFitPerTP1","nHits fit / TP Pln 1", 10, 0, 10);
      fnHitsFitPerTP[2] = tfs->make<TH1F>("nHitsFitPerTP2","nHits fit / TP Pln 2", 10, 0, 10);
      
      fDelta[0] = tfs->make<TH1F>("Delta0","Delta Pln 0", 50, 0, 2);
      fDelta[1] = tfs->make<TH1F>("Delta1","Delta Pln 1", 50, 0, 2);
      fDelta[2] = tfs->make<TH1F>("Delta2","Delta Pln 2", 50, 0, 2);
      
      fDeltaN[0] = tfs->make<TH1F>("DeltaN0","Normalized Delta Pln 0", 50, 0, 4);
      fDeltaN[1] = tfs->make<TH1F>("DeltaN1","Normalized Delta Pln 1", 50, 0, 4);
      fDeltaN[2] = tfs->make<TH1F>("DeltaN2","Normalized Delta Pln 2", 50, 0, 4);

      fCharge[0] = tfs->make<TH1F>("Charge0","Charge/Pt Pln 0", 100, 0, 500);
      fCharge[1] = tfs->make<TH1F>("Charge1","Charge/Pt Pln 1", 100, 0, 500);
      fCharge[2] = tfs->make<TH1F>("Charge2","Charge/Pt Pln 2", 100, 0, 500);

      fnHitsPerTP_Angle[0] = tfs->make<TH2F>("nhtpertp_angle0","Hits/TP vs Angle Pln 0", 10, 0 , M_PI/2, 9, 1, 10);
      fnHitsPerTP_Angle[1] = tfs->make<TH2F>("nhtpertp_angle1","Hits/TP vs Angle Pln 1", 10, 0 , M_PI/2, 9, 1, 10);
      fnHitsPerTP_Angle[2] = tfs->make<TH2F>("nhtpertp_angle2","Hits/TP vs Angle Pln 2", 10, 0 , M_PI/2, 9, 1, 10);
      
      fnHitsPerTP_AngleP[0] = tfs->make<TProfile>("nhtpertp_anglep0","Hits/TP vs Angle Pln 0", 10, 0 , M_PI/2, "S");
      fnHitsPerTP_AngleP[1] = tfs->make<TProfile>("nhtpertp_anglep1","Hits/TP vs Angle Pln 1", 10, 0 , M_PI/2, "S");
      fnHitsPerTP_AngleP[2] = tfs->make<TProfile>("nhtpertp_anglep2","Hits/TP vs Angle Pln 2", 10, 0 , M_PI/2, "S");

    }
    
    if(fShowerStudy) {
      fShowerNumTrjint = tfs->make<TH1F>("showernumtrjint","Shower Num Traj Intersections",100, 0, 200);
      fShowerDVtx = tfs->make<TH1F>("showerdvtx","Shower dVtx",100, 0, 50);
      fShowerTheta_Sep = tfs->make<TH2F>("showertheta_sep","Shower dTheta vs Sep",40, 0, 4, 10, 0, 1);
      fShowerDVtx_Sep = tfs->make<TH2F>("showerdvtx_sep","Shower dVtx vs Sep",40, 0, 4, 10, 0, 1);
    }
    
  }
  
  //------------------------------------------------------------------------------
  void TrajClusterAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
 
    bool badinput = false;
    fHitFinderModuleLabel = pset.get<art::InputTag>("HitFinderModuleLabel");
    fMode                 = pset.get< short >("Mode", 0); // Default is don't use it
    fHitErrFac            = pset.get< float >("HitErrFac", 0.4);
    fMinAmp               = pset.get< float >("MinAmp", 5);
    fLargeAngle           = pset.get< float >("LargeAngle", 80);
    fNPtsAve              = pset.get< short >("NPtsAve", 20);
    fMinPtsFit            = pset.get< std::vector<unsigned short >>("MinPtsFit");
    fMinPts               = pset.get< std::vector<unsigned short >>("MinPts");
    fLAStep               = pset.get< std::vector<bool>>("LAStep");
    fMaxChi               = pset.get< float >("MaxChi", 10);
    fMultHitSep           = pset.get< float >("MultHitSep", 2.5);
    fRecoveryAlgs         = pset.get< float >("RecoveryAlgs", 0);
    fReversePropagate     = pset.get< float >("ReversePropagate", 0);
    fHitFOMCut            = pset.get< float >("HitFOMCut", 3);
    fChgRatCut            = pset.get< float >("ChgRatCut", 0.5);
    fKinkAngCut           = pset.get< float >("KinkAngCut", 0.4);
    fMaxWireSkipNoSignal  = pset.get< float >("MaxWireSkipNoSignal", 1);
    fMaxWireSkipWithSignal= pset.get< float >("MaxWireSkipWithSignal", 100);
    fProjectionErrFactor  = pset.get< float >("ProjectionErrFactor", 2);
    fJTMaxHitSep2         = pset.get< float >("JTMaxHitSep", 2);
    
    fStudyMode            = pset.get< bool  >("StudyMode", false);
    fShowerStudy          = pset.get< bool  >("ShowerStudy", false);
    fTagAllTraj           = pset.get< bool  >("TagAllTraj", false);
    fMaxTrajSep           = pset.get< float >("MaxTrajSep", 4);
    fShowerPrtPlane       = pset.get< short >("ShowerPrtPlane", -1);
    fVertex2DIPCut        = pset.get< float  >("Vertex2DIPCut", -1);
    fVertex3DChiCut       = pset.get< float  >("Vertex3DChiCut", -1);
    fMaxVertexTrajSep     = pset.get< std::vector<float>>("MaxVertexTrajSep");
    
    fDebugPlane         = pset.get< int  >("DebugPlane", -1);
    fDebugWire          = pset.get< int  >("DebugWire", -1);
    fDebugHit           = pset.get< int  >("DebugHit", -1);
    
    // convert angle (degrees) into a direction cosine cut in the wire coordinate
    // It should be in the range 0 < fLargeAngle < 90
    fLargeAngle = cos(fLargeAngle * M_PI / 180);
    
    // convert the max traj separation into a separation^2
    fMaxTrajSep *= fMaxTrajSep;
    if(fJTMaxHitSep2 > 0) fJTMaxHitSep2 *= fJTMaxHitSep2;
    
    if(fMinPtsFit.size() != fMinPts.size()) badinput = true;
    if(fMaxVertexTrajSep.size() != fMinPts.size()) badinput = true;
    if(fLAStep.size() != fMinPts.size()) badinput = true;
    if(badinput) throw art::Exception(art::errors::Configuration)<< "TrajClusterAlg: Bad input from fcl file. Vector lengths are not the same";
    
    if(kAlgBitSize != AlgBitNames.size())
      throw art::Exception(art::errors::Configuration)<<"kAlgBitSize "<<kAlgBitSize<<" != AlgBitNames size "<<AlgBitNames.size()<<"\n";

    
  } // reconfigure
  
  // used for sorting hits on wires
//  bool SortByLowHit(unsigned int i, unsigned int j) {return ((i > j));}
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ClearResults() {
    allTraj.clear();
    inTraj.clear();
    trial.clear();
    tjphs.clear();
    inClus.clear();
    tcl.clear();
    vtx.clear();
    vtx3.clear();
  } // ClearResults()

  ////////////////////////////////////////////////
  void TrajClusterAlg::RunTrajClusterAlg(art::Event & evt)
  {

    if(fMode == 0) return;
    
    //Get the hits for this event:
//    art::Handle< std::vector<recob::Hit> > hitVecHandle;
//    evt.getByLabel(fHitFinderModuleLabel, hitVecHandle);
    
    art::ValidHandle< std::vector<recob::Hit>> hitVecHandle
    = evt.getValidHandle<std::vector<recob::Hit>>(fHitFinderModuleLabel);

    fHits.resize(hitVecHandle->size());
    if(fHits.size() == 0) return;
 
    larprop = lar::providerFrom<detinfo::LArPropertiesService>();
    detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    
    for (unsigned int iht = 0; iht < fHits.size(); iht++) fHits[iht] = art::Ptr< recob::Hit>(hitVecHandle, iht);
    
    ClearResults();
    // set all hits to the available state
    inTraj.resize(fHits.size(), 0);

    // TODO Insert hit sorting code here
    
    std::cout<<"Event "<<evt.event()<<"\n";

    fQuitAlg = false;
    fIsRealData = evt.isRealData();
    
    for (geo::TPCID const& tpcid: geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
        didPrt = false;
        WireHitRange.clear();
        // define a code to ensure clusters are compared within the same plane
        fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
        fCstat = tpcid.Cryostat;
        fTpc = tpcid.TPC;
        // Calculate fScaleF, the scale factor to convert a tick into
        // Wire Spacing Equivalent (WSE) units where the wire spacing in this plane = 1
        raw::ChannelID_t channel = fHits[0]->Channel();
        float wirePitch = geom->WirePitch(geom->View(channel));
        float tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
        tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
        fScaleF = tickToDist / wirePitch;
        // fill the WireHitRange vector with first/last hit on each wire
        // dead wires and wires with no hits are flagged < 0
        GetHitRange();
        // no hits on this plane?
        if(fFirstWire == fLastWire) continue;

        if(fMode < 2) {
          fStepDir = fMode;
          // make all trajectories with this setting
          ReconstructAllTraj();
          if(fQuitAlg) {
            ClearResults();
            return;
          }
          // tag shower-like clusters/hits
          if(fTagAllTraj) TagAllTraj();
          if(fQuitAlg) {
            ClearResults();
            return;
          }
        } else {
          // here is where we get fancy
          // Step in the + direction
          fStepDir = 1;
          ReconstructAllTraj();
          if(fQuitAlg) {
            ClearResults();
            return;
          }
          // tag shower-like clusters/hits after the first call
          if(fTagAllTraj) TagAllTraj();
          if(fQuitAlg) {
            ClearResults();
            return;
          }
          // store the results
          trial.push_back(allTraj);
          inTrialTraj.push_back(inTraj);
          // Kill the trajectories and reconstruct again with different parameters
          KillAllTrajInCTP(fCTP);
//          for(auto& intj : inTraj) intj = 0;
          // Step in the - direction
          fStepDir = -1;
          ReconstructAllTraj();
          if(fQuitAlg) {
            ClearResults();
            return;
          }
          // store the results
          trial.push_back(allTraj);
          inTrialTraj.push_back(inTraj);
          AnalyzeTrials();
//          AdjudicateTrials(reAnalyze);
        }
      } // fPlane
      
      if(fVertex3DChiCut > 0) Find3DVertices(tpcid);
      trial.clear();
      inTrialTraj.clear();
    } // tpcid
    
    //    bool reAnalyze = false;
 
    // Convert trajectories in allTraj into clusters
    MakeAllTrajClusters();
    if(fQuitAlg) {
      ClearResults();
      return;
    }
    CheckHitClusterAssociations();
    if(fQuitAlg) {
      ClearResults();
      return;
    }
    
    if(didPrt || fDebugPlane >= 0) {
      mf::LogVerbatim("TC")<<"Done in ReconstructAllTraj";
      PrintAllTraj(USHRT_MAX, 0);
    }
    
    // temp
    if(fStudyMode) {
      unsigned short ipl, itj, nInPln;
      float ang, ave;
      std::vector<unsigned int> hitVec;
      for(ipl = 0; ipl < 3; ++ipl) {
        nInPln = 0;
        for(itj = 0; itj < allTraj.size(); ++itj) if(allTraj[itj].CTP == ipl) ++nInPln;
        if(nInPln != 1) continue;
        for(itj = 0; itj < allTraj.size(); ++itj) {
          if(allTraj[itj].CTP != ipl) continue;
//          PrintTrajectory(allTraj[itj], USHRT_MAX);
          for(auto& tp : allTraj[itj].Pts) {
            if(tp.Chg == 0) continue;
            fnHitsPerTP[ipl]->Fill(tp.Hits.size());
            fnHitsFitPerTP[ipl]->Fill(NumUsedHits(tp));
            fDelta[ipl]->Fill(tp.Delta);
            if(tp.HitPosErr2 > 0) fDeltaN[ipl]->Fill(tp.Delta/ sqrt(tp.HitPosErr2));
            fCharge[ipl]->Fill(tp.Chg);
          }
          ang = std::abs(work.Pts[0].Ang);
          if(ang > M_PI/2) ang = M_PI - ang;
          // find average number of hits / TP
          PutTrajHitsInVector(allTraj[itj], true, hitVec);
          ave = (float)hitVec.size() / (float)(allTraj[itj].EndPt[1] - allTraj[itj].EndPt[0] + 1);
          fnHitsPerTP_Angle[ipl]->Fill(ang, ave);
          fnHitsPerTP_AngleP[ipl]->Fill(ang, ave);
        } // itj
      } // ipl
    } // studymode

    // convert vertex time from WSE to ticks
    for(auto& avtx : vtx) avtx.Time /= fScaleF;
    
  } // RunTrajClusterAlg
  
  //////////////////////////////////
  void TrajClusterAlg::KillAllTrajInCTP(CTP_t tCTP)
  {
    // Kills all trajectories in the provided tCTP
    // and releases the hits
    unsigned short itj;
    std::vector<unsigned int> thits;
    for(itj = 0; itj < allTraj.size(); ++itj) {
      if(allTraj[itj].AlgMod[kKilled]) continue;
      if(allTraj[itj].CTP != tCTP) continue;
      allTraj[itj].AlgMod[kKilled] = true;
      if(tCTP == fCTP) {
        std::cout<<"KillAllTrajInCTP: write some code\n";
        exit(1);
      } else {
        std::cout<<"KillAllTrajInCTP: write some code\n";
        exit(1);
      }
    } // itj
  } // KillAllTrajInCTP

  ////////////////////////////////////////////////
  void TrajClusterAlg::AnalyzeTrials()
  {
    // Analyze the Set of All Trajectories to construct a single set of trajectories
    // that is the best. The results are stored in allTraj.
    
    // This shouldn't happen but do it anyway
    if(trial.size() == 1) {
      allTraj = trial[0];
      return;
    }
    mf::LogVerbatim myprt("TC");

    unsigned short itrial, itj, jtrial, jtj, nSameHits = 0;
/*
    unsigned short nHitsInTraj;
    for(itrial = 0; itrial < trial.size(); ++itrial) {
      nHitsInTraj = 0;
      for(itj = 0; itj < trial[itrial].size(); ++itj) nHitsInTraj += trial[itrial][itj].NHits;
      myprt<<"Plane "<<fPlane<<" trial "<<itrial<<" NHits in all trajectories "<<nHitsInTraj<<"\n";
    } // is
*/
    std::vector<unsigned int> iHitVec, jHitVec;
    tjphs.clear();
    TjPairHitShare tmp;
    for(itrial = 0; itrial < trial.size() - 1; ++itrial) {
      for(itj = 0; itj < trial[itrial].size(); ++itj) {
        // ignore obsolete trajectories
        if(trial[itrial][itj].AlgMod[kKilled]) continue;
        // look at long trajectories for testing
        if(trial[itrial][itj].EndPt[1] < 5) continue;
        PutTrajHitsInVector(trial[itrial][itj], true, iHitVec);
//        std::cout<<"itr "<<itr<<" itj "<<itj<<" hit size "<<iHitVec.size()<<" nSameHits ";
        for(jtrial = itrial + 1; jtrial < trial.size(); ++jtrial) {
          for(jtj = 0; jtj < trial[jtrial].size(); ++jtj) {
            if(trial[jtrial][jtj].AlgMod[kKilled]) continue;
            PutTrajHitsInVector(trial[jtrial][jtj], true, jHitVec);
            CountSameHits(iHitVec, jHitVec, nSameHits);
//            std::cout<<" "<<nSameHits;
            if(nSameHits == 0) continue;
            tmp.iTrial = itrial; tmp.iTj = itj;
            tmp.jTrial = jtrial; tmp.jTj = jtj;
            tmp.nSameHits = nSameHits;
            tjphs.push_back(tmp);
            // keep track of the
          } // jtj
        } // jtr
//        std::cout<<"\n";
      } // itj
    } // itr
    
    myprt<<"  itr  itj    jtr  jtj   nSameHits\n";
    for(auto& tmp : tjphs) {
      itj = tmp.iTj; jtj = tmp.jTj;
      myprt<<std::setw(5)<<tmp.iTrial<<std::setw(5)<<itj;
      myprt<<std::setw(5)<<tmp.jTrial<<std::setw(5)<<jtj;
      myprt<<std::setw(6)<<tmp.nSameHits<<"\n";
      PrintTrajectory(trial[tmp.iTrial][itj], USHRT_MAX);
      PrintTrajectory(trial[tmp.jTrial][jtj], USHRT_MAX);
    }

  } // AnalyzeTrials
/*
  ////////////////////////////////////////////////
  void TrajClusterAlg::AdjudicateTrials(bool& reAnalyze)
  {
    // returns redo true if AnalyzeTrials needs to be called again
    
    if(tjphs.size() == 0) return;
    
    unsigned short ipr, ii, ndup, itr, itj;
    for(ipr = 0; ipr < tjphs.size(); ++ipr) {
      // look for a trajectory that only appears once in iTrial
      // ensure that itj hasn't been clobbered
      itr = tjphs[ipr].iTrial;
      itj = tjphs[ipr].iTj;
      if(trial[itr][itj].ProcCode == USHRT_MAX) continue;
      ndup = 0;
      for(ii = 0; ii < tjphs.size(); ++ii) {
        if(tjphs[ii].iTrial != tjphs[ipr].iTrial) continue;
        if(tjphs[ii].iTj != tjphs[ipr].iTj) continue;
        ++ndup;
      } // ii
      std::cout<<"Adjudicate "<<ipr<<" ndup "<<ndup<<"\n";
      if(ndup != 1) continue;
      // TODO Should we ensure that jtj is unique also?
      MergeTrajPair(ipr, reAnalyze);
      if(reAnalyze) return;
    } // ipr
    
  } // AdjudicateTrials

  ////////////////////////////////////////////////
  void TrajClusterAlg::MergeTrajPair(unsigned short ipr, bool& reAnalyze)
  {
    // merge the trajectories and put the results into allTraj. Returns
    // reAnalyze true if AnalyzeTrials needs to be called again
    
    Trajectory iTj = trial[tjphs[ipr].iTrial][tjphs[ipr].iTj];
    Trajectory jTj = trial[tjphs[ipr].jTrial][tjphs[ipr].jTj];
    
    mf::LogVerbatim("TC")<<"MergeTrajPair "<<tjphs[ipr].iTrial<<"-"<<tjphs[ipr].iTj<<" and "<<tjphs[ipr].jTrial<<"-"<<tjphs[ipr].jTj;
    PrintTrajectory(iTj, USHRT_MAX);
    PrintTrajectory(jTj, USHRT_MAX);
    std::vector<float> tSep;
    TrajSeparation(iTj, jTj, tSep);
    
  } // MergeTrajPair
*/
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
  void TrajClusterAlg::PutTrajHitsInVector(Trajectory const& tj, bool onlyUsedHits, std::vector<unsigned int>& hitVec)
  {
    // Put hits in each trajectory point into a flat vector. Only hits with UseHit if onlyUsedHits == true
    hitVec.clear();
    unsigned short ipt, iht;
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      for(iht = 0; iht < tj.Pts[ipt].Hits.size(); ++iht) {
        if(onlyUsedHits) {
          if(tj.Pts[ipt].UseHit[iht]) hitVec.push_back(tj.Pts[ipt].Hits[iht]);
        } else {
          hitVec.push_back(tj.Pts[ipt].Hits[iht]);
        }
      } // iht
    } // ipt
  } // PutTrajHitsInVector
   
  ////////////////////////////////////////////////
  void TrajClusterAlg::ReconstructAllTraj()
  {
    // Reconstruct clusters in fPlane and put them in allTraj
    
    unsigned int ii, iwire, jwire, iht, jht, oht;
    
    unsigned int nwires = fLastWire - fFirstWire - 1;
    unsigned int ifirsthit, ilasthit, jfirsthit, jlasthit;
    float fromWire, fromTick, toWire, toTick, deltaRms, qtot;
    bool sigOK;

    for(unsigned short pass = 0; pass < fMinPtsFit.size(); ++pass) {
      fPass = pass;
      for(ii = 0; ii < nwires; ++ii) {
        // decide which way to step given the sign of fStepDir
        if(fStepDir > 0) {
          // step DS
          iwire = fFirstWire + ii;
          jwire = iwire + 1;
        } else {
          // step US
          iwire = fLastWire - ii - 1;
          jwire = iwire - 1;
        }
        // skip bad wires or no hits on the wire
        if(WireHitRange[iwire].first < 0) continue;
        if(WireHitRange[jwire].first < 0) continue;
        ifirsthit = (unsigned int)WireHitRange[iwire].first;
        ilasthit = (unsigned int)WireHitRange[iwire].second;
        jfirsthit = (unsigned int)WireHitRange[jwire].first;
        jlasthit = (unsigned int)WireHitRange[jwire].second;
        for(iht = ifirsthit; iht < ilasthit; ++iht) {
          // clear out any leftover work inTraj's that weren't cleaned up properly
          for(oht = ifirsthit; oht < ilasthit; ++oht) if(inTraj[oht] < 0) inTraj[oht] = 0;
          prt = (fDebugPlane == (int)fPlane && (int)iwire == fDebugWire && std::abs((int)fHits[iht]->PeakTime() - fDebugHit) < 10);
          if(prt) didPrt = true;
//          if(!prt) continue;
          if(prt) mf::LogVerbatim("TC")<<"+++++++ Pass "<<fPass<<" Found debug hit "<<fPlane<<":"<<PrintHit(iht)<<" inTraj "<<inTraj[iht]<<" RMS "<<fHits[iht]->RMS()<<" Multiplicity "<<fHits[iht]->Multiplicity()<<" LocalIndex "<<fHits[iht]->LocalIndex();
          if(inTraj[iht] != 0) continue;
          fromWire = fHits[iht]->WireID().Wire;
          fromTick = fHits[iht]->PeakTime();
          qtot = fHits[iht]->Integral();
          // decide whether to combine two hits in a multiplet into one HitPos
          if(fHits[iht]->Multiplicity() == 2) {
            // get the index of the other hit in the multiplet
            oht = 1 - fHits[iht]->LocalIndex();
            // Make sure it isn't used and check the separation
            if(inTraj[oht] == 0 && std::abs(fHits[oht]->PeakTime() - fromTick) < fMultHitSep * fHits[oht]->RMS()) {
              HitMultipletPosition(iht, fromTick, deltaRms, qtot);
              if(qtot == 0) continue;
              if(prt) mf::LogVerbatim("TC")<<" Hit doublet: back from HitMultipletPosition ";
            }
          }
          for(jht = jfirsthit; jht < jlasthit; ++jht) {
            // clear out any leftover work inTraj's that weren't cleaned up properly
            for(oht = jfirsthit; oht < jlasthit; ++oht) if(inTraj[oht] < 0) inTraj[oht] = 0;
            if(inTraj[iht] != 0) continue;
            fHitDoublet = false;
//            if(inTraj[jht] == -3) mf::LogVerbatim("TC")<<"ReconstructAllTraj bad jht flag "<<fPlane<<":"<<PrintHit(jht)<<" using iht "<<fPlane<<":"<<PrintHit(iht);
            if(inTraj[jht] != 0) continue;
            if(prt) mf::LogVerbatim("TC")<<"+++++++ checking ClusterHitsOK with jht "<<fPlane<<":"<<PrintHit(jht)<<" RMS "<<fHits[jht]->RMS()<<" Multiplicity "<<fHits[jht]->Multiplicity()<<" LocalIndex "<<fHits[jht]->LocalIndex();
            // Ensure that the hits StartTick and EndTick have the proper overlap
            if(!TrajHitsOK(iht, jht)) continue;
            // start a trajectory in the direction from iht -> jht
            toWire = jwire;
            toTick = fHits[jht]->PeakTime();
            if(fHits[jht]->Multiplicity() > 1) {
              HitMultipletPosition(jht, toTick, deltaRms, qtot);
              if(qtot == 0) continue;
              if(fHits[jht]->Multiplicity() > 2 && SkipHighMultHitCombo(iht, jht)) continue;
            }
            if(prt) mf::LogVerbatim("TC")<<"  Starting trajectory from "<<(int)fromWire<<":"<<(int)fromTick<<" to "<<(int)toWire<<":"<<(int)toTick;
            StartWork(fromWire, fromTick, toWire, toTick, fCTP);
            // check for a major failure
            if(fQuitAlg) return;
            // check for a failure
            if(work.Pts.size() == 0) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: StartWork failed";
              ReleaseWorkHits();
              continue;
            }
            // check for a large angle crawl
            if(IsLargeAngle(work.Pts[0]) && !fLAStep[fPass]) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: No LA stepping on this pass";
              ReleaseWorkHits();
              continue;
            }
            fHitDoublet = (fHits[iht]->Multiplicity() == 2) && (fHits[jht]->Multiplicity() == 2) && IsLargeAngle(work.Pts[0]);
            if(fHitDoublet) {
              // Adjust the TP delta RMS
              work.Pts[0].DeltaRMS = deltaRms;
            } // fHitDoublet
            // try to add close hits
            AddHits(work, 0, sigOK);
            // check for a major failure
            if(fQuitAlg) return;
            if(!sigOK || NumUsedHits(work.Pts[0]) == 0) {
              if(prt) mf::LogVerbatim("TC")<<" No hits at initial trajectory point ";
              ReleaseWorkHits();
              continue;
            }
            // print the header and the first TP
            if(prt) PrintTrajectory(work, USHRT_MAX);
            // We can't update the trajectory yet because there is only one TP.
            // Define the "average charge" so we can at least get going
//            work.Pts[0].AveChg = work.Pts[0].Chg;
            work.EndPt[0] = 0;
            // now try stepping away
            StepCrawl();
            // check for a major failure
            if(fQuitAlg) return;
            if(prt) mf::LogVerbatim("TC")<<" After first StepCrawl. fGoodWork "<<fGoodWork<<" fTryWithNextPass "<<fTryWithNextPass;
            if(!fGoodWork && fTryWithNextPass) {
              StepCrawl();
              if(!fUpdateTrajOK) {
                if(prt) mf::LogVerbatim("TC")<<" xxxxxxx StepCrawl failed AGAIN. fTryWithNextPass "<<fTryWithNextPass;
                ReleaseWorkHits();
                continue;
              } // Failed again
            }
            // Check the quality of the work trajectory
            CheckWork();
            // check for a major failure
            if(fQuitAlg) return;
            if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: After CheckWork EndPt "<<work.EndPt[0]<<"-"<<work.EndPt[1]<<" fGoodWork "<<fGoodWork<<" fTryWithNextPass "<<fTryWithNextPass;
            if(fTryWithNextPass) {
              // The first part of the trajectory was good but the latter part
              // had too many unused hits in the work TPs. The work vector was
              // truncated and fPass incremented, so give it another try
              work.AlgMod[kTryWithNextPass] = true;
              StepCrawl();
              // check for a major failure
              if(fQuitAlg) return;
              if(!fGoodWork) {
                if(prt) mf::LogVerbatim("TC")<<" xxxxxxx StepCrawl failed AGAIN after CheckWork";
                ReleaseWorkHits();
                continue;
              } // Failed again
            } // fTryWithNextPass
            if(prt) mf::LogVerbatim("TC")<<"StepCrawl done: work.EndPt[1] "<<work.EndPt[1]<<" cut "<<fMinPts[work.Pass];
            if(fRecoveryAlgs > 0) TryRecoveryAlgs();
            if(!fGoodWork) continue;
            // decide if the trajectory is long enough
            if(NumPtsWithCharge(work) < fMinPts[work.Pass]) {
              if(prt) mf::LogVerbatim("TC")<<" xxxxxxx Not enough points "<<NumPtsWithCharge(work)<<" minimum "<<fMinPts[work.Pass];
              ReleaseWorkHits();
              continue;
            }
            if(fReversePropagate > 0) ReversePropagate(work);
            if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: calling StoreWork with npts "<<work.EndPt[1];
            StoreWork();
            // check for a major failure
            if(fQuitAlg) return;
            break;
          } // jht
          if(inTraj[iht] > 0) break;
        } // iht
      } // iwire
      Find2DVertices();
    } // fPass
    
    if(fQuitAlg) return;
    
    // make junk trajectories using nearby un-assigned hits
    if(fJTMaxHitSep2 > 0) {
      FindJunkTraj();
      // check for a major failure
      if(fQuitAlg) return;
      Find2DVertices();
      // check for a major failure
      if(fQuitAlg) return;
    }
    
    FindHammerVertices();
    
    // last attempt to attach Tjs to vertices
    unsigned short lastPass = fMaxVertexTrajSep.size() - 1;
    for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) if(vtx[ivx].NTraj > 0) AttachAnyTrajToVertex(ivx, fMaxVertexTrajSep[lastPass], false );
    
    prt = false;
    // put MC info into the trajectory struct
//    if(!fIsRealData) FillTrajTruth();
    work.Pts.clear();
    
    if(didPrt) PrintAllTraj(USHRT_MAX, USHRT_MAX);

    CheckInTraj("RCAT");
    if(fQuitAlg) return;
    
    
  } // ReconstructAllTraj

  
  //////////////////////////////////////////
  void TrajClusterAlg::TryRecoveryAlgs()
  {
    // Check the work trajectory and try to recover it if it has poor quality
    

    if(prt) mf::LogVerbatim("TC")<<"inside TryRecoveryAlgs ";

    if(work.Pts.size() < fMinPts[work.Pass]) return;
    
    // a short shallow angle trajectory which has hit multiplets but the
    // hit multiplets haven't been "merged". Just try using every associated
    // hit and see if we can get back on the road again
    unsigned short ii, ipt;
    unsigned int iht;
    bool hitAdded = false;
    if(!fGoodWork && work.Pts.size() < 10) {
      for(auto& tp : work.Pts) {
        hitAdded = false;
        for(ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tp.UseHit[ii]) continue;
          iht = tp.Hits[ii];
          if(inTraj[iht] > 0) continue;
          hitAdded = true;
          inTraj[iht] = work.ID;
          tp.UseHit[ii] = true;
        } // ii
        if(hitAdded) DefineHitPos(tp);
      } // tp
      SetEndPoints(work);
      if(NumPtsWithCharge(work) > 2) {
        FitWork();
        ipt = work.Pts.size() - 1;
        if(prt) mf::LogVerbatim("TC")<<" Use all hits in short tj. fit chisq "<<work.Pts[ipt].FitChi;
        if(work.Pts[ipt].FitChi > fMaxChi) return;
        work.AlgMod[kRecovery1] = true;
        // keep on stepping
        StepCrawl();
        return;
      }
    } // short traj with unmerged hits
    
    if(work.Pts.size() < 10) return;
    
    // The most common situation is that the trajectory starts out OK but as
    // hits are added, they are close but not used. Start by counting the
    // number of close but unused hits going backwards through the trajectory.
    unsigned short nCloseHits = 0, nBadPts = 0, nSingleCloseHits = 0;
//    unsigned short firstHighMultPt = USHRT_MAX;
    for(ii = 0; ii < work.Pts.size(); ++ii) {
      ipt = work.Pts.size() - 1 - ii;
      if(work.Pts[ipt].Chg > 0) break;
      nCloseHits += work.Pts[ipt].Hits.size();
      if(work.Pts[ipt].Hits.size() == 1) {
        // only count available close hits
        iht = work.Pts[ipt].Hits[0];
        if(inTraj[iht] == 0) ++nSingleCloseHits;
      }
//      if(work.Pts[ipt].Hits.size() > 1 && firstHighMultPt == USHRT_MAX) firstHighMultPt = ipt;
      ++nBadPts;
    } // ii
    
    if(nBadPts == 0) return;
    if(prt) mf::LogVerbatim("TC")<<"TryRecoveryAlgs: nBadPts "<<nBadPts<<" nCloseHits "<<nCloseHits<<" nSingleCloseHits "<<nSingleCloseHits;
    
    // Require a minimum number of TPs that could be used to step in the opposite direction
    if(nSingleCloseHits < 5) return;
    
    // Use a copy of work in case something bad happens
    Trajectory twork = work;

    // Use single hit TPs at the end
    unsigned short nused = 0, iptBreak = 0;
    for(ii = 0; ii < twork.Pts.size(); ++ii) {
      ipt = twork.Pts.size() - 1 - ii;
      if(twork.Pts[ipt].Chg > 0) break;
      if(twork.Pts[ipt].Hits.size() != 1) continue;
      iht = twork.Pts[ipt].Hits[0];
      if(inTraj[iht] > 0) continue;
      inTraj[iht] = work.ID;
      twork.Pts[ipt].UseHit[0] = true;
      DefineHitPos(twork.Pts[ipt]);
      if(prt) mf::LogVerbatim("TC")<<" Use hit "<<PrintHit(iht)<<" HitPos "<<twork.Pts[ipt].HitPos[0]<<" "<<twork.Pts[ipt].HitPos[1];
      ++nused;
      iptBreak = ipt;
      if(nused == 4) break;
    } // ii
    // release the points from the beginning to iptBreak
    for(ipt = 0; ipt < iptBreak; ++ipt) UnsetUsedHits(twork.Pts[ipt]);
    work.AlgMod[kRecovery2] = true;
    
    // Here is where we might slot in additional algs
    
  } // TryRecoveryAlgs

  //////////////////////////////////////////
  void TrajClusterAlg::ReversePropagate(Trajectory& tj)
  {
    // Reverse the trajectory and propagate the fit from the end
    // to the beginning, possibly masking off hits that were associated
    // with TPs along the way
    
    // assume we aren't going to do this
    fRevProp = false;
    if(!prt) return;
    
    if(NumPtsWithCharge(tj) < 5) return;
    
    if(prt) mf::LogVerbatim("TC")<<"inside ReversePropagate";
    PrintTrajectory(tj, USHRT_MAX);
    
    // Find 4 TPs at the end of to start the process
    unsigned short ii, ipt, iptStart = USHRT_MAX, cnt = 0;
    for(ii = 0; ii < tj.Pts.size() - 1; ++ii) {
      ipt = tj.Pts.size() - 1 - ii;
      if(NumUsedHits(tj.Pts[ipt]) == 0) continue;
      ++cnt;
      if(cnt == 4) {
        iptStart = ipt;
        break;
      }
    } // ii
    if(iptStart == USHRT_MAX) return;
    
    // TODO will need to reverse the order of hits for large angle tjs
    
    // copy tj into a work tj
    Trajectory workTj = tj;
    // change the ID to keep track of what has been changed
    workTj.ID -= 1;
    // resize the number of trajectory points
    mf::LogVerbatim("TC")<<"iptStart "<<iptStart;
    workTj.Pts.erase(workTj.Pts.begin(), workTj.Pts.begin() + iptStart);
    // then reverse it
    ReverseTraj(workTj);
    unsigned short lastPt = workTj.Pts.size() - 1;
    // Decrement the number of points fit since it will be incremented in UpdateTraj
    // using the previous TP NTPsFit
    workTj.Pts[lastPt - 1].NTPsFit = cnt - 1;
    // These lines are just the tail end of AddHits after hits were inserted in Pts
    FindUseHits(workTj, lastPt);
    SetEndPoints(workTj);
//    DefineHitPos(workTj.Pts[iptStart]);
    // This line is ala StepCrawl
    UpdateTraj(workTj);
    PrintTrajectory(workTj, USHRT_MAX);
    mf::LogVerbatim("TC")<<"Start propagating";
    if(workTj.Pts[lastPt].FitChi > fMaxChi) return;
    fRevProp = true;
    // Append TPs from tj to workTj, find hits to use and fit
    for(ii = 1; ii < tj.Pts.size(); ++ii) {
      ipt = iptStart - ii;
      workTj.Pts.push_back(tj.Pts[ipt]);
      lastPt = workTj.Pts.size() - 1;
      UnsetUsedHits(workTj.Pts[lastPt]);
      FindUseHits(workTj, lastPt);
      SetEndPoints(workTj);
      UpdateTraj(workTj);
      PrintTrajectory(workTj, lastPt);
      if(ipt == 0) break;
    }
    fRevProp = false;
//    if(prt) PrintTrajectory(workTj, USHRT_MAX);
    // TODO be smarter about this maybe
    tj = workTj;
    tj.AlgMod[kRevProp] = true;
    
  } // ReversePropagate
  
/*
  //////////////////////////////////////////
  void TrajClusterAlg::ReversePropagate(Trajectory& tj)
  {
    // This is a StepCrawl-like routine but it only uses hits that are
    // already associated with the trajectory. The calling routine should have
    // set UseHit and inTraj appropriately and called DefineHitPos. The trajectory
    // EndPt information is ignored and not reliable unless this algorithm is
    // successful. Note that fGoodWork is used as a success/fail flag even though
    // we may not be reverse propagating the work trajectory
 
    fGoodWork = false;
    unsigned short endPt = tj.Pts.size() - 1;
    unsigned short startPt;
    for(startPt = 0; startPt < endPt; ++startPt) if(tj.Pts[startPt].Chg > 0) break;
    // we need at least two points to get started
    if(startPt == endPt) return;
    
    // do a first fit
    unsigned short npts = endPt - startPt;
    TrajPoint tpFit;
    FitTraj(tj, startPt, npts, 1, tpFit);
    if(prt) mf::LogVerbatim("TC")<<"ReversePropagate: first fit "<<tpFit.FitChi<<" origin "<<startPt<<" npts "<<npts;
    if(tpFit.FitChi > fMaxChi) return;
    tpFit.NTPsFit = npts;
    // overwrite the startPt
    tj.Pts[startPt] = tpFit;
    // update the other fitted TPs with the same information
    unsigned short ipt;
    for(ipt = startPt + 1; ipt <= endPt; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      tj.Pts[ipt].Delta = PointTrajDOCA(tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[0], tpFit);
      tj.Pts[ipt].FitChi = tpFit.FitChi;
      tj.Pts[ipt].NTPsFit = tpFit.NTPsFit;
      tj.Pts[ipt].Pos = tpFit.Pos;
      tj.Pts[ipt].Dir = tpFit.Dir;
    } // ipt
    // Now back propagate. The trajectory is projected to point ipt
    // using the trajectory stored in projPt.
    unsigned short projPt = startPt;
    if(prt) {
      mf::LogVerbatim("TC")<<"start reverse propagation";
      PrintTrajectory(tj, startPt);
    }
    // temp TP
    TrajPoint tp;
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      ipt = startPt - ii;
      // Update the position and direction at TP ipt
      tp = tj.Pts[projPt];
      MoveTPToWire(tp, tj.Pts[ipt].Pos[0]);
      tj.Pts[ipt].Pos = tp.Pos;
      tj.Pts[ipt].Dir = tp.Dir;
      if(prt) PrintTrajectory(tj, ipt);
      // update the used hits
      FindUseHits(tj, ipt);
      if(tj.Pts[ipt].Chg == 0) continue;
      SetEndPoints(tj);
      DefineHitPos(tj.Pts[ipt]);
      UpdateTraj(tj, ipt, projPt, 1);
      if(!fGoodWork) return;
      projPt = ipt;
      if(ipt == 0) break;
    } // jj
    
    if(prt) {
      mf::LogVerbatim("TC")<<"Reverse propagation done";
      PrintTrajectory(tj, USHRT_MAX);

    }
  } // ReversePropagate
*/
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindJunkTraj()
  {
    // Makes junk trajectories using unassigned hits
    
    unsigned int iwire, ifirsthit, ilasthit, iht;
    unsigned int jwire, jfirsthit, jlasthit, jht;
    unsigned int kwire, kfirsthit, klasthit, kht;
    unsigned int fromIndex;
    unsigned int loWire, hiWire;
    float loTime, hiTime;
    bool hitsAdded;
    unsigned short nit, tht;
    
    // shouldn't have to do this but...
    for(iht = 0; iht < fHits.size(); ++iht) if(inTraj[iht] < 0) inTraj[iht] = 0;
    
    std::vector<unsigned int> tHits;
    for(iwire = fFirstWire; iwire < fLastWire; ++iwire) {
      // skip bad wires or no hits on the wire
      if(WireHitRange[iwire].first < 0) continue;
      jwire = iwire + 1;
      if(WireHitRange[jwire].first < 0) continue;
      ifirsthit = (unsigned int)WireHitRange[iwire].first;
      ilasthit = (unsigned int)WireHitRange[iwire].second;
      jfirsthit = (unsigned int)WireHitRange[jwire].first;
      jlasthit = (unsigned int)WireHitRange[jwire].second;
      for(iht = ifirsthit; iht < ilasthit; ++iht) {
        prt = true;
        prt = (fDebugPlane == (int)fPlane && (int)iwire == fDebugWire && std::abs((int)fHits[iht]->PeakTime() - fDebugHit) < 10);
        if(prt) {
          mf::LogVerbatim("TC")<<"FindJunkTraj: Found debug hit "<<PrintHit(iht)<<" inTraj "<<inTraj[iht]<<" fJTMaxHitSep2 "<<fJTMaxHitSep2;
        }
        if(inTraj[iht] != 0) continue;
        for(jht = jfirsthit; jht < jlasthit; ++jht) {
          if(inTraj[jht] != 0) continue;
          if(prt && HitSep2(iht, jht) < 100) mf::LogVerbatim("TC")<<" use "<<PrintHit(jht)<<" HitSep2 "<<HitSep2(iht, jht);
          if(HitSep2(iht, jht) > fJTMaxHitSep2) continue;
          tHits.clear();
          // add all hits and flag them
          fromIndex = iht - fHits[iht]->LocalIndex();
          for(kht = fromIndex; kht < fromIndex + fHits[iht]->Multiplicity(); ++kht) {
            if(inTraj[kht] != 0) continue;
            tHits.push_back(kht);
            inTraj[kht] = -4;
          } // kht
          fromIndex = jht - fHits[jht]->LocalIndex();
          for(kht = fromIndex; kht < fromIndex + fHits[jht]->Multiplicity(); ++kht) {
            if(inTraj[kht] != 0) continue;
            tHits.push_back(kht);
            inTraj[kht] = -4;
          } // kht
          if(iwire != 0) { loWire = iwire - 1; } else { loWire = 0; }
          if(jwire < fNumWires - 1) { hiWire = jwire + 2; } else { hiWire = fNumWires; }
          hitsAdded = true;
          nit = 0;
          while(hitsAdded && nit < 100) {
            hitsAdded = false;
            for(kwire = loWire; kwire < hiWire + 1; ++kwire) {
              if(WireHitRange[kwire].first < 0) continue;
              kfirsthit = (unsigned int)WireHitRange[kwire].first;
              klasthit = (unsigned int)WireHitRange[kwire].second;
              for(kht = kfirsthit; kht < klasthit; ++kht) {
                if(inTraj[kht] != 0) continue;
                // this shouldn't be needed but do it anyway
                if(std::find(tHits.begin(), tHits.end(), kht) != tHits.end()) continue;
                // check w every hit in tHit
                for(tht = 0; tht < tHits.size(); ++tht) {
//                  if(prt && HitSep2(kht, tHits[tht]) < 100) mf::LogVerbatim("TC")<<" kht "<<PrintHit(kht)<<" tht "<<PrintHit(tHits[tht])<<" HitSep2 "<<HitSep2(kht, tHits[tht])<<" cut "<<fJTMaxHitSep2;
                  if(HitSep2(kht, tHits[tht]) > fJTMaxHitSep2) continue;
                  hitsAdded = true;
                  tHits.push_back(kht);
                  inTraj[kht] = -4;
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
            if(fHits[tHits[tht]]->WireID().Wire < loWire) loWire = fHits[tHits[tht]]->WireID().Wire;
            if(fHits[tHits[tht]]->WireID().Wire > hiWire) hiWire = fHits[tHits[tht]]->WireID().Wire;
            if(fHits[tHits[tht]]->PeakTime() < loTime) loTime = fHits[tHits[tht]]->PeakTime();
            if(fHits[tHits[tht]]->PeakTime() > hiTime) hiTime = fHits[tHits[tht]]->PeakTime();
          }
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<" tHits";
            for(auto tht : tHits) myprt<<" "<<PrintHit(tht);
            myprt<<"\n";
          } // prt
          MakeJunkTraj(tHits);
          // release any hits that weren't included in a trajectory
          for(auto iht : tHits) if(inTraj[iht] == -4) inTraj[iht] = 0;
          if(hitsAdded) break;
        } // jht
      } // iht
    } // iwire
  } // FindJunkTraj

  //////////////////////////////////////////
  void TrajClusterAlg::MakeJunkTraj(std::vector<unsigned int> tHits)
  {
    
     // Make a crummy trajectory using the provided hits
    if(tHits.size() < 2) return;

    std::vector<std::vector<unsigned int>> tpHits;
    unsigned short ii, iht, ipt;
    
    // Start the work trajectory using the first and last hits to
    // define a starting direction
    StartWork(tHits[0], tHits[tHits.size()-1]);
    
    // Do a more complicated specification of TP hits if there
    // are a large number of them
    if(tHits.size() > 8) {
      // fit all of the hits to a line
      std::vector<float> x(tHits.size()), y(tHits.size()), yerr2(tHits.size());
      float intcpt, slope, intcpterr, slopeerr, chidof, qtot = 0;
      
      for(ii = 0; ii < tHits.size(); ++ii) {
        iht = tHits[ii];
        x[ii] = fHits[iht]->WireID().Wire;
        y[ii] = fHits[iht]->PeakTime() * fScaleF;
        qtot += fHits[iht]->Integral();
        yerr2[ii] = 1;
      } // ii
      fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
//      std::cout<<" chidof "<<chidof<<" slope "<<slope<<"\n";
      // return without making a junk trajectory
      if(chidof > 900) return;
      // A rough estimate of the trajectory angle
      work.Pts[0].Ang = atan(slope);
      // Rotate the hits into this coordinate system to find the start and end
      // points and general direction
      float cs = cos(-work.Pts[0].Ang);
      float sn = sin(-work.Pts[0].Ang);
      float tAlong, minAlong = 1E6, maxAlong = -1E6;
      // sort the hits by the distance along the general direction
//      std::cout<<" x size "<<x.size()<<" slope "<<slope<<"\n";
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
      unsigned short npts = (unsigned short)(maxAlong - minAlong);
      // rotate back into normal coordinate system
      if(prt) mf::LogVerbatim("TC")<<" minAlong "<<minAlong<<" maxAlong "<<maxAlong<<" work.Pts[0].Ang "<<work.Pts[0].Ang<<" npts "<<npts;
      if(npts < 2) npts = 2;
      tpHits.resize(npts);
      for(ii = 0; ii < tHits.size(); ++ii) {
        ipt = (unsigned short)(sortVec[ii].length - minAlong);
        if(ipt > npts - 1) ipt = npts - 1;
        if(prt) mf::LogVerbatim("TC")<<"tHit "<<PrintHit(tHits[ii])<<" length "<<sortVec[ii].length<<" ipt "<<ipt;
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
    work.Pts[0].UseHit.resize(work.Pts[0].Hits.size(), true);
    DefineHitPos(work.Pts[0]);
    PrintTrajectory(work, USHRT_MAX);
    // make the rest of the TPs
    for(ipt = 1; ipt < tpHits.size(); ++ipt) {
      if(tpHits[ipt].size() == 0) continue;
      // Use the first TP as a starting point
      TrajPoint tp = work.Pts[0];
      tp.Step = ipt;
      tp.Hits = tpHits[ipt];
      // use all hits
      tp.UseHit.resize(tp.Hits.size(), true);
      DefineHitPos(tp);
      // Just use the hit position as the traj position
      tp.Pos = tp.HitPos;
      work.Pts.push_back(tp);
      SetEndPoints(work);
//      if(work.Pts.size() > 1) UpdateTraj(work);
    }
    if(prt) PrintTrajectory(work, USHRT_MAX);
    work.AlgMod[kJunkTj] = true;
    // Finally push it onto allTraj
    StoreWork();
    
  } // MakeJunkTraj

  
  //////////////////////////////////////////
  void TrajClusterAlg::FillTrajTruth()
  {
    
    sim::ParticleList plist;
    art::ServiceHandle<cheat::BackTracker> bt;
    plist = bt->ParticleList();
    std::vector<const simb::MCParticle*> plist2;

    unsigned int trackID;
    float KE;
    std::vector<unsigned int> tidlist;
 
    for(sim::ParticleList::const_iterator ipart = plist.begin(); ipart != plist.end(); ++ipart) {
      const simb::MCParticle* part = (*ipart).second;
      assert(part != 0);
      trackID = part->TrackId();
      // KE in MeV
      KE = 1000 * (part->E() - part->Mass());
      if(KE < 10) continue;
//      mf::LogVerbatim("TC")<<"TrackID "<<trackID<<" pdg "<<part->PdgCode()<<" E "<<part->E()<<" mass "<<part->Mass()<<" KE "<<KE<<" Mother "<<part->Mother()<<" Proc "<<part->Process();
      tidlist.push_back(trackID);
      plist2.push_back(part);
    }
    
    if(tidlist.size() == 0) return;
    
    // count of the number of Geant tracks used in each trajectory hit
    std::vector<unsigned short> tidcnt(tidlist.size());
    unsigned short itj, ipt, ii, iht, jj, maxCnt;
    std::vector<sim::IDE> sIDEs;
    for(itj = 0; itj < allTraj.size(); ++itj) {
      if(allTraj[itj].AlgMod[kKilled]) continue;
      sIDEs.clear();
      for(ii = 0; ii < tidcnt.size(); ++ii) tidcnt[ii] = 0;
      for(ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
        if(allTraj[itj].Pts[ipt].Chg == 0) continue;
        for(ii = 0; ii < allTraj[itj].Pts[ipt].Hits.size(); ++ii) {
          iht = allTraj[itj].Pts[ipt].Hits[ii];
          bt->HitToSimIDEs(fHits[iht], sIDEs);
          if(sIDEs.size() == 0) continue;
          for(jj = 0; jj < sIDEs.size(); ++jj) {
            trackID = sIDEs[jj].trackID;
            for(unsigned short kk = 0; kk < tidlist.size(); ++kk) {
              if(trackID == tidlist[kk]) {
                ++tidcnt[kk];
                break;
              }
            } // kk
          } // jj
        } // ii
      } // ipt
      // find the Geant track ID with the most entries. Require at
      // least 50% of the hits be associated with the track
      maxCnt = allTraj[itj].EndPt[1] / 2;
      trackID = tidcnt.size();
      for(ii = 0; ii < tidcnt.size(); ++ii) {
        if(tidcnt[ii] > maxCnt) {
          maxCnt = tidcnt[ii];
          trackID = ii;
        }
      } // ii
      if(trackID > tidlist.size()-1) continue;
//      mf::LogVerbatim("TC")<<"Traj "<<itj<<" match with trackID "<<trackID;
      allTraj[itj].TruPDG = plist2[trackID]->PdgCode();
      allTraj[itj].TruKE = 1000 * (plist2[trackID]->E() - plist2[trackID]->Mass());
      if(plist2[trackID]->Process() == "primary") allTraj[itj].IsPrimary = true;
    } // itj

    
  } // FillTrajTruth

  //////////////////////////////////////////
  void TrajClusterAlg::TagAllTraj()
  {
    // try to tag as shower-like or track-like
    
    if(allTraj.size() < 2) return;
    
    shPrt = (fShowerPrtPlane == (short)fPlane);
    if(shPrt) didPrt = true;
    
    if(shPrt) {
      mf::LogVerbatim("TC")<<"Inside TagAllTraj: plane "<<fPlane;
      // hijack fDebugPlane for a bit
      int tmp = fDebugPlane;
      fDebugPlane = fShowerPrtPlane;
      PrintAllTraj(USHRT_MAX, 0);
      fDebugPlane = tmp;
    }
    
    trjint.clear();
    TrjInt aTrjInt;
    ClsOfTrj.clear();
    
    // maximum separation^2
    float maxSep2 = fMaxTrajSep;
    float minSep2;
    float dang, vw, vt, dw, dt, dvtx2;
    std::vector<std::array<unsigned short, 3>> nCloseEnd(allTraj.size());
    unsigned short i1, i2, ipt1, ipt2;
    unsigned short endPt1, endPt2;
    unsigned short bin1, bin2;
    
    for(i1 = 0; i1 < allTraj.size() - 1; ++i1) {
      Trajectory& tj1 = allTraj[i1];
      if(tj1.AlgMod[kKilled]) continue;
       if(tj1.CTP != fPlane) continue;
      for(i2 = i1 + 1; i2 < allTraj.size(); ++i2) {
        Trajectory& tj2 = allTraj[i2];
        if(tj2.AlgMod[kKilled]) continue;
        if(tj2.CTP != tj1.CTP) continue;
        // find the closest approach
        minSep2 = maxSep2;
        TrajTrajDOCA(tj1, tj2, ipt1, ipt2, minSep2);
        if(minSep2 == maxSep2) continue;
        // Count the number at each end and in the middle
        bin1 = (unsigned short)(3 * (float)ipt1 / (float)tj1.EndPt[1]);
        if(bin1 > 2) bin1 = 2;
        // only count if this end doesn't have a vertex
        if(tj1.Vtx[bin1] < 0) ++nCloseEnd[i1][bin1];
        bin2 = (unsigned short)(3 * (float)ipt2 / (float)tj2.EndPt[1]);
        if(bin2 > 2) bin2 = 2;
        if(tj2.Vtx[bin2] < 0) ++nCloseEnd[i2][bin2];
        // find the angle between the TPs at the intersection
        dang = std::abs(tj1.Pts[ipt1].Ang - tj2.Pts[ipt2].Ang);
        if(bin1 != 1  && bin2 != 1) {
          // the DOCA point is at the ends of the two TJs.
          // Find the intersection using the appropriate end points
          endPt1 = allTraj[i1].EndPt[0];
          if(bin1 == 2) endPt1 = allTraj[i1].EndPt[1];
          endPt2 = allTraj[i2].EndPt[0];
          if(bin2 == 2) endPt2 = allTraj[i2].EndPt[1];
          TrajIntersection(allTraj[i1].Pts[endPt1], allTraj[i2].Pts[endPt2], vw, vt);
//          if(shPrt) mf::LogVerbatim("TC")<<"TI check i1 "<<i1<<" endPt1 "<<endPt1<<" Ang "<<allTraj[i1].Pts[endPt1].Ang<<" i2 "<<i2<<" endPt2 "<<endPt2<<" Ang "<<allTraj[i2].Pts[endPt2].Ang<<" W:T "<<(int)vw<<":"<<(int)vt/fScaleF;
        } else {
          vw = -1;
          vt = -1;
        }
        // distance between the vertex position and the closest separation
        dw = vw - allTraj[i1].Pts[ipt1].Pos[0];
        dt = vt - allTraj[i1].Pts[ipt1].Pos[1];
        dvtx2 = dw * dw + dt * dt;
        float dangErr = allTraj[i1].Pts[ipt1].AngErr;
        if(allTraj[i2].Pts[ipt2].AngErr > dangErr) dangErr = allTraj[i2].Pts[ipt2].AngErr;
//        float dangSig = dang / dangErr;
//        if(shPrt) mf::LogVerbatim("TC")<<"Close "<<i1<<" bin1 "<<bin1<<" i2 "<<i2<<" bin2 "<<bin2<<" minSep2 "<<minSep2<<" dang "<<dang<<" fKinkAngCut "<<fKinkAngCut<<" dangSig "<<dangSig<<" dvtx2 "<<dvtx2;
        // save it
        aTrjInt.itj1 = i1;
        aTrjInt.ipt1 = ipt1;
        aTrjInt.itj2 = i2;
        aTrjInt.ipt2 = ipt2;
        aTrjInt.sep2 = minSep2;
        aTrjInt.dang = dang;
        aTrjInt.vw = vw;
        aTrjInt.vt = vt;
        trjint.push_back(aTrjInt);
        if(fShowerStudy) {
          float sep = sqrt(minSep2);
          fShowerTheta_Sep->Fill(sep, dang);
          float dvtx = sqrt(dvtx2);
          fShowerDVtx->Fill(dvtx);
          fShowerDVtx_Sep->Fill(sep, dvtx);
        }
      } // jj
    } // ii
    
    if(fShowerStudy) {
      fShowerNumTrjint->Fill((float)trjint.size());
    }
    
    if(trjint.size() == 0) return;
    if(shPrt) {
      for(i1 = 0; i1 < trjint.size(); ++i1) {
        mf::LogVerbatim("TC")<<i1<<" trjint "<<" "<<trjint[i1].itj1<<" "<<trjint[i1].itj2<<" sep2 "<<trjint[i1].sep2<<" dang "<<trjint[i1].dang<<" vtx "<<fPlane<<":"<<(int)trjint[i1].vw<<":"<<(int)(trjint[i1].vt / fScaleF);
      }
    }
    
    std::vector<std::vector<unsigned short>> trjintIndices;
    FindClustersOfTrajectories(trjintIndices);
    
    unsigned short icot;
    // Deal with each COT (Note that ClsOfTrj and trjintIndices have the same size)
    for(icot = 0; icot < ClsOfTrj.size(); ++icot) {
      // try fitting the DOCA points to a line to get the shower axis if there are enough points
      if(trjintIndices[icot].size() > 3) {
        DefineShowerTraj(icot, trjintIndices);
      } // trjintIndices[icot].size() > 3
      else {
        mf::LogVerbatim("TC")<<" no shower angle and endpoint information available... Write some code";
      }
    } // icot
    
    // clean up
    trjint.clear();
    ClsOfTrj.clear();
    
    if(shPrt) {
      mf::LogVerbatim("TC")<<"Inside TagAllTraj: plane "<<fPlane;
      // hijack fDebugPlane for a bit
      int tmp = fDebugPlane;
      fDebugPlane = fShowerPrtPlane;
      PrintAllTraj(USHRT_MAX, 0);
      fDebugPlane = tmp;
    }
    
  } // TagAllTraj
  
  //////////////////////////////////////////
  void TrajClusterAlg::DefineShowerTraj(unsigned short icot, std::vector<std::vector<unsigned short>> trjintIndices)
  {
    // Try to define a shower trajectory consisting of a single track-like trajectory (the electron/photon)
    // plus a group of shower-like trajectories, using the vector of trjintIndices
    
    unsigned short ii, iTji, i1, ipt1, i2, ipt2, nEndInside, ipt0;
    std::vector<float> x, y, yerr2;
    float arg, fom, pos0, pos1, minsep, dang;
    float showerAngle, cs, sn;
    float intcpt, slope, intcpterr, slopeerr, chidof, shlong;
    unsigned short primTraj, primTrajEnd;
    for(ii = 0; ii < trjintIndices[icot].size(); ++ii) {
      iTji = trjintIndices[icot][ii];
      i1 = trjint[iTji].itj1;
      ipt1 = trjint[iTji].ipt1;
      x.push_back(allTraj[i1].Pts[ipt1].Pos[0]);
      y.push_back(allTraj[i1].Pts[ipt1].Pos[1]);
      // weight by the length of both trajectories
      i2 = trjint[iTji].itj2;
      arg = allTraj[i1].EndPt[1] + allTraj[i2].EndPt[1];
      yerr2.push_back(arg);
    } // ii
    fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    showerAngle = atan(slope);
    if(shPrt) mf::LogVerbatim("TC")<<"Fit intcpt "<<intcpt<<" slope "<<slope<<" angle "<<showerAngle<<" chidof "<<chidof;
    // Rotate the intersection points into a coordinate system along the
    // shower direction. We will use this to find the (restricted) longitudinal and
    // transverse extent of the shower
    cs = cos(-showerAngle);
    sn = sin(-showerAngle);
    // min and max of traj intersections along the shower angle direction
    float minshlong = 1E6, maxshlong = -1;
    unsigned short miniTji = 0, maxiTji = 0;
    for(ii = 0; ii < trjintIndices[icot].size(); ++ii) {
      iTji = trjintIndices[icot][ii];
      i1 = trjint[iTji].itj1;
      ipt1 = trjint[iTji].ipt1;
      shlong = cs * allTraj[i1].Pts[ipt1].Pos[0] - sn * allTraj[i1].Pts[ipt1].Pos[1];
      //          std::cout<<"i1 "<<i1<<" ipt1 "<<ipt1<<" pos "<<(int)allTraj[i1].Pts[ipt1].Pos[0]<<":"<<(int)allTraj[i1].Pts[ipt1].Pos[1]<<" shlong "<<shlong<<"\n";
      if(shlong < minshlong) { minshlong = shlong; miniTji = iTji; } // shlong < minshlong
      if(shlong > maxshlong) { maxshlong = shlong; maxiTji = iTji; } // shlong > maxshlong
      i2 = trjint[iTji].itj2;
      ipt2 = trjint[iTji].ipt2;
      shlong = cs * allTraj[i2].Pts[ipt2].Pos[0] - sn * allTraj[i2].Pts[ipt2].Pos[1];
      //          std::cout<<"i2 "<<i2<<" ipt2 "<<ipt2<<" pos "<<(int)allTraj[i2].Pts[ipt2].Pos[0]<<":"<<(int)allTraj[i2].Pts[ipt2].Pos[1]<<" shlong "<<shlong<<"\n";
      if(shlong < minshlong) { minshlong = shlong; miniTji = iTji; } // shlong < minshlong
      if(shlong > maxshlong) { maxshlong = shlong; maxiTji = iTji; } // shlong > maxshlong
    } // ii
    // reduce the extend by 1 WSE at each end to allow some tolerance
    minshlong += 1;
    maxshlong -= 1;
    if(shPrt) {
      mf::LogVerbatim("TC")<<"minshlong "<<minshlong<<" maxshlong "<<maxshlong;
      mf::LogVerbatim("TC")<<"Primary traj is probably in "<<miniTji<<" or "<<maxiTji;
    }
    // Form a list of candidates
    std::vector<unsigned short> primCandidates(4);
    primCandidates[0] = trjint[miniTji].itj1;
    primCandidates[1] = trjint[miniTji].itj2;
    primCandidates[2] = trjint[maxiTji].itj1;
    primCandidates[3] = trjint[maxiTji].itj2;
    fom = 999;
    for(ii = 0; ii < 4; ++ii) {
      i1 = primCandidates[ii];
      ipt1 = allTraj[i1].EndPt[1];
      // reject if both end points are inside the shower boundaries.
      // Rotate the endpoints into the shower coordinate system
      nEndInside = 0;
      shlong = cs * allTraj[i1].Pts[0].Pos[0] - sn * allTraj[i1].Pts[0].Pos[1];
      if(shlong > minshlong && shlong < maxshlong) ++nEndInside;
      if(shPrt) mf::LogVerbatim("TC")<<" Candidate min "<<i1<<" shlong "<<shlong<<" nEndInside "<<nEndInside;
      // rotate the other end and test
      shlong = cs * allTraj[i1].Pts[ipt1].Pos[0] - sn * allTraj[i1].Pts[ipt1].Pos[1];
      if(shlong > minshlong && shlong < maxshlong) ++nEndInside;
      if(shPrt) mf::LogVerbatim("TC")<<" Candidate max "<<i1<<" shlong "<<shlong<<" nEndInside "<<nEndInside;
      if(nEndInside > 1) continue;
      // average the angle of the TP end points and find the difference
      // wrt the shower angle
      dang = std::abs(0.5 * (allTraj[i1].Pts[0].Ang + allTraj[i1].Pts[ipt1].Ang) - showerAngle);
      arg = 10 * (1 + allTraj[i1].Pass) * dang * allTraj[i1].EndPt[1];
//      mf::LogVerbatim("TC")<<"Candidate "<<i1<<" nCloseEnd "<<nCloseEnd[i1][0]<<" "<<nCloseEnd[i1][1]<<" "<<nCloseEnd[i1][2]<<" pass "<<allTraj[i1].Pass<<" dang "<<dang<<" arg "<<arg;
      if(arg < fom) {
        fom = arg;
        primTraj = i1;
      }
    } // ii
    // determine which end is closest to the shower end
    // distance between TP0 and the shower min position
//    ipt1 = allTraj[primTraj].Pts.size()-1;
    ipt0 = allTraj[primTraj].EndPt[0];
    ipt1 = allTraj[primTraj].EndPt[1];
    pos0 = cs * allTraj[primTraj].Pts[ipt0].Pos[0] - sn * allTraj[primTraj].Pts[ipt0].Pos[1];
    pos1 = cs * allTraj[primTraj].Pts[ipt1].Pos[0] - sn * allTraj[primTraj].Pts[ipt1].Pos[1];
    minsep = std::abs(pos0 - minshlong);
    primTrajEnd = 0;
    //        mf::LogVerbatim("TC")<<"0min "<<pos0<<" minsep "<<minsep<<" end "<<primTrajEnd;
    arg = std::abs(pos1 - minshlong);
    if(arg < minsep) { minsep = arg;  primTrajEnd = ipt1; }
    //        mf::LogVerbatim("TC")<<"1min "<<pos1<<" arg "<<arg<<" end "<<primTrajEnd;
    arg = std::abs(pos0 - maxshlong);
    if(arg < minsep) { minsep = arg;  primTrajEnd = ipt0; }
    //        mf::LogVerbatim("TC")<<"0max "<<pos0<<" minsep "<<minsep<<" end "<<primTrajEnd;
    arg = std::abs(pos1 - maxshlong);
    if(arg < minsep) { minsep = arg;  primTrajEnd = ipt1; }
    //        mf::LogVerbatim("TC")<<"1max "<<pos1<<" arg "<<arg<<" end "<<primTrajEnd;
    if(shPrt) mf::LogVerbatim("TC")<<" set Primary traj = "<<primTraj<<" end "<<primTrajEnd;
    TagShowerTraj(icot, primTraj, primTrajEnd, showerAngle);

  }

  //////////////////////////////////////////
  void TrajClusterAlg::TagShowerTraj(unsigned short icot, unsigned short primTraj, unsigned short primTrajEnd, float showerAngle)
  {
    // Tag a track-like trajectory using primTraj and a shower-like cluster consisting
    // of all other trajectories in ClsOfTrj[icot]
    
    // set the PDG for all TJ's to shower-like
    unsigned short ii, itj;
    for(ii = 0; ii < ClsOfTrj[icot].size(); ++ii) {
      itj = ClsOfTrj[icot][ii];
      allTraj[itj].PDG = 12;
      allTraj[itj].ParentTraj = primTraj;
    }
    
    allTraj[primTraj].PDG = 13;
    allTraj[primTraj].ParentTraj = USHRT_MAX;
    
  } // MakeShowerClusters
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindClustersOfTrajectories(std::vector<std::vector<unsigned short>>& trjintIndices)
  {
    // Associate trajectories that are close to each into Clusters Of Trajectories (COTs)
    
    mf::LogVerbatim myprt("TC");

    if(trjint.size() == 0) return;
    ClsOfTrj.clear();
    trjintIndices.clear();
    std::vector<unsigned short> tmp;
    std::vector<unsigned short> imp;
    std::vector<bool> inCOT(allTraj.size(), false);
    bool tjAdded;
    unsigned short iTji, itj1, itj2;
    for(unsigned short nit = 0; nit < 100; ++nit) {
      tmp.clear();
      imp.clear();
      // look for trajectories that haven't been assigned to a COT
      for(iTji = 0; iTji < trjint.size(); ++iTji) {
        itj1 = trjint[iTji].itj1;
        itj2 = trjint[iTji].itj2;
        if(!inCOT[itj1]) {
          tmp.push_back(itj1);
          imp.push_back(iTji);
          inCOT[itj1] = true;
          break;
        }
        if(!inCOT[itj2]) {
          tmp.push_back(itj2);
          imp.push_back(iTji);
          inCOT[itj2] = true;
          break;
        }
      } // itji
      // no new un-assigned trajectories?
      if(tmp.size() == 0) break;
       // try to add more TJ's to the COT
      tjAdded = true;
      while(tjAdded) {
        tjAdded = false;
        for(iTji = 0; iTji < trjint.size(); ++iTji) {
          itj1 = trjint[iTji].itj1;
          itj2 = trjint[iTji].itj2;
          // itj2 available &&  found itj1 in tmp
          if(!inCOT[itj2] && (std::find(tmp.begin(), tmp.end(), itj1) != tmp.end())) {
            // add itj2 to tmp
            tmp.push_back(itj2);
            inCOT[itj2] = true;
            imp.push_back(iTji);
            tjAdded = true;
          } // !inCOT[itj1]
          // itj1 available &&  found itj2 in tmp
          if(!inCOT[itj1] && (std::find(tmp.begin(), tmp.end(), itj2) != tmp.end())) {
            // add itj1 to tmp
            tmp.push_back(itj1);
            inCOT[itj1] = true;
            imp.push_back(iTji);
            tjAdded = true;
          } // !inCOT[itj1]
        } // iTji
      } // tjAdded
      ClsOfTrj.push_back(tmp);
      trjintIndices.push_back(imp);
    } // nit
    
    if(shPrt) {
      myprt<<"FindClustersOfTrajectories list of COTs\n";
      for(unsigned short ii = 0; ii < ClsOfTrj.size(); ++ii) {
        myprt<<ii<<" COT ";
        for(unsigned short jj = 0; jj < ClsOfTrj[ii].size(); ++jj) myprt<<" "<<ClsOfTrj[ii][jj];
        myprt<<"\n";
      } // ii
      myprt<<"not in a COT ";
      for(itj1 = 0; itj1 < inCOT.size(); ++itj1) {
        if(!inCOT[itj1] && allTraj[itj1].CTP == fPlane) myprt<<" "<<itj1;
      }
    } // shPrt
    
  } // FindClustersOfTrajectories
 
  ////////////////////////////////////////////////
  void TrajClusterAlg::AddHits(Trajectory& tj, unsigned short ipt, bool& sigOK)
  {
    // Try to add hits to the trajectory point ipt on the supplied
    // trajectory
    
    fAddedBigDeltaHit = false;

    if(tj.Pts.size() == 0) return;
    if(ipt > tj.Pts.size() - 1) return;
    
    std::vector<unsigned int> closeHits;
    unsigned int wire, loWire, hiWire, iht, firstHit, lastHit;

    unsigned int lastPtWithHits = tj.EndPt[1];
    unsigned int prevPt = 0;
    // This is going to fail if at some point we decide to use this code
    // to add hits to a point that is not at the leading edge of the TJ
    if(tj.Pts.size() > 0) prevPt = tj.Pts.size() - 1;
    TrajPoint& tp = tj.Pts[ipt];

    // figure out which wires to consider
    // On the first entry only consider the wire the TP is on
    if(tj.Pts.size() == 1) {
      loWire = std::nearbyint(tp.Pos[0]);
      hiWire = loWire + 1;
    } else  {
      if(IsLargeAngle(tp)) {
        // look at adjacent wires for larger angle trajectories
        loWire = std::nearbyint(tp.Pos[0] - 1);
        hiWire = loWire + 3;
      } else {
        // not large angle
        loWire = std::nearbyint(tp.Pos[0]);
        hiWire = loWire + 1;
        // Move the TP to this wire
        MoveTPToWire(tp, (float)loWire);
      }
    } // tj.Pts.size > 1
    
    float fwire, ftime, delta;
    
    // find the projection error to this point. Note that if this is the first
    // TP, lastPtWithHits = 0, so the projection error is 0
    float dw = tp.Pos[0] - tj.Pts[lastPtWithHits].Pos[0];
    float dt = tp.Pos[1] - tj.Pts[lastPtWithHits].Pos[1];
    float dpos = sqrt(dw * dw + dt * dt);
    float projErr = dpos * tj.Pts[lastPtWithHits].AngErr;
    // Add this to the Delta RMS factor and construct a cut
    float deltaCut = 3 * (projErr + tp.DeltaRMS);
    
    bool isLA = IsLargeAngle(tp);
    // Very Large Angle
    bool isVLA = std::abs(tp.Dir[0]) < fLargeAngle - 0.1;
    if(isVLA) {
      // VLA deltaCut is delta / hit RMS
      deltaCut = 2;
    } else if(isLA) {
      if(deltaCut < 0.6) deltaCut = 0.6;
    } else {
      if(deltaCut < 0.1) deltaCut = 0.1;
    }
    deltaCut *= fProjectionErrFactor;
    
    // TEMP
    deltaCut = fProjectionErrFactor + 0.2 * std::abs(dw);
//    float bigDelta = 5 * deltaCut;
//    unsigned int imBig = UINT_MAX;
    
    // projected time in ticks for testing the existence of a hit signal
    raw::TDCtick_t rawProjTick = (float)(tp.Pos[1] / fScaleF);
    
    // put the close hits from previous TPs in to a vector to facilitate searching
    // to ensure that hits are not associated with multiple TPs. This is only required
    // for large angle trajectories
    std::vector<unsigned int> pclHits;
    // put all hits in the vector
    PutTrajHitsInVector(tj, false, pclHits);
/*
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"pclHits in ";
      for(auto iht : pclHits) myprt<<" "<<PrintHit(iht);
    }
*/
    // assume failure
    sigOK = false;
    if(prt) mf::LogVerbatim("TC")<<" AddHits: loWire "<<loWire<<" tp.Pos[0] "<<tp.Pos[0]<<" hiWire "<<hiWire<<" projTick "<<rawProjTick<<" deltaRMS "<<tp.DeltaRMS<<" tp.Dir[0] "<<tp.Dir[0]<<" isLA "<<isLA<<" isVLA "<<isVLA<<" deltaCut "<<deltaCut;
    
    for(wire = loWire; wire < hiWire; ++wire) {
      // Assume a signal exists on a dead wire
      if(WireHitRange[wire].first == -1) sigOK = true;
      if(WireHitRange[wire].first < 0) continue;
      firstHit = (unsigned int)WireHitRange[wire].first;
      lastHit = (unsigned int)WireHitRange[wire].second;
      fwire = wire;
      for(iht = firstHit; iht < lastHit; ++iht) {
        if(rawProjTick > fHits[iht]->StartTick() && rawProjTick < fHits[iht]->EndTick()) sigOK = true;
        ftime = fScaleF * fHits[iht]->PeakTime();
        delta = PointTrajDOCA(fwire, ftime, tp);
        float dt = std::abs(ftime - tp.Pos[1]);
        if(prt && delta < 100 && dt < 200) {
          mf::LogVerbatim myprt("TC");
          myprt<<"  chk "<<fHits[iht]->WireID().Plane<<":"<<PrintHit(iht);
          myprt<<" delta "<<std::fixed<<std::setprecision(2)<<delta<<" deltaCut "<<deltaCut<<" dt "<<dt;
          myprt<<" Mult "<<fHits[iht]->Multiplicity()<<" localIndex "<<fHits[iht]->LocalIndex()<<" RMS "<<std::setprecision(1)<<fHits[iht]->RMS();
          myprt<<" Chi "<<std::setprecision(1)<<fHits[iht]->GoodnessOfFit();
          myprt<<" inTraj "<<inTraj[iht];
          myprt<<" Chg "<<(int)fHits[iht]->Integral();
          myprt<<" Signal? "<<sigOK;
        }
        if(isVLA) {
          // Very Large Angle
          // Cut on dt using the RMS of a crude hit
          if(fHits[iht]->GoodnessOfFit() < 0) {
            if(dt > 2 * fHits[iht]->RMS()) continue;
          } else {
            if(dt > 3) continue;
          }
        } else if(isLA) {
          // Large Angle
          // The impact parameter delta may be good but we may be projecting
          // the trajectory too far away (in time) from the current position.
          // The LA step size is 1 so make the cut a bit larger than that.
          if(dt > 1.1) continue;
          if(delta > deltaCut) continue;
        } else {
          // Not large angle
          if(delta > deltaCut) continue;
        }
        sigOK = true;
        // No requirement is made at this point that a hit cannot be associated with
        // more than 1 TP in different trajectories, but make sure that it isn't assigned
        // more than once to this trajectory
        if(std::find(pclHits.begin(), pclHits.end(), iht) != pclHits.end()) continue;
        closeHits.push_back(iht);
        pclHits.push_back(iht);
        if(fHits[iht]->Multiplicity() > 1) {
          // include all the hits in a multiplet for not large angle TPs
          unsigned int loHit = iht - fHits[iht]->LocalIndex();
          unsigned int hiHit = loHit + fHits[iht]->Multiplicity();
          for(unsigned int jht = loHit; jht < hiHit; ++jht) {
            if(std::find(pclHits.begin(), pclHits.end(), jht) != pclHits.end()) continue;
            closeHits.push_back(jht);
            pclHits.push_back(jht);
          } // jht
        } // multiplicity > 1
      } // iht
    } // wire
    if(!sigOK) {
      if(prt) mf::LogVerbatim("TC")<<" no signal on any wire at tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" tick "<<(int)tp.Pos[1]/fScaleF<<" closeHits size "<<closeHits.size();
      return;
    }
/*
    if(imBig == UINT_MAX) return;
    if(closeHits.size() == 0 && std::find(pclHits.begin(), pclHits.end(), imBig) == pclHits.end()) {
      // use the closest hit available
      closeHits.push_back(imBig);
      if(!isLA && fHits[imBig]->Multiplicity() > 1) {
        // include all the hits in a multiplet for not large angle TPs
        unsigned int loHit = imBig - fHits[imBig]->LocalIndex();
        unsigned int hiHit = loHit + fHits[imBig]->Multiplicity();
        for(unsigned int jht = loHit; jht < hiHit; ++jht) {
          if(std::find(pclHits.begin(), pclHits.end(), jht) != pclHits.end()) continue;
          closeHits.push_back(jht);
          pclHits.push_back(jht);
        } // jht
      } // multiplicity > 1
      if(prt) mf::LogVerbatim("TC")<<" Added bigDelta hit "<<PrintHit(imBig)<<" w delta = "<<bigDelta;
    }
*/
    tp.Hits.insert(tp.Hits.end(), closeHits.begin(), closeHits.end());
    // sort the close hits by distance from the previous traj point
    if(tp.Hits.size() > 1) {
      std::vector<SortEntry> sortVec;
      SortEntry sortEntry;
      unsigned short ii;
      float dw, dt;
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        sortEntry.index = ii;
        iht = tp.Hits[ii];
        dw = fHits[iht]->WireID().Wire - tj.Pts[prevPt].Pos[0];
        dt = fHits[iht]->PeakTime() * fScaleF - tj.Pts[prevPt].Pos[1];
        sortEntry.length = dw * dw + dt * dt;
        sortVec.push_back(sortEntry);
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), lessThan);
      // make a temp vector
      std::vector<unsigned int> tmp(sortVec.size());
      // overwrite with the sorted values
      for(ii = 0; ii < sortVec.size(); ++ii) tmp[ii] = tp.Hits[sortVec[ii].index];
      // replace
      tp.Hits = tmp;
    }
    // resize the UseHit vector and assume that none of these hits will be used (yet)
    tp.UseHit.resize(tp.Hits.size(), false);
    // decide which of these hits should be used in the fit
    FindUseHits(tj, ipt);
    SetEndPoints(tj);
    DefineHitPos(tp);
    if(prt) mf::LogVerbatim("TC")<<" number of close hits "<<closeHits.size()<<" used hits "<<NumUsedHits(tp);
/*
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"pclHits out";
      for(auto iht : pclHits) myprt<<" "<<PrintHit(iht);
    }
*/
  } // AddHits
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::ReleaseWorkHits()
  {
    // Sets inTraj[] = 0 and UseHit false for all TPs in work. Called when abandoning work
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt)  UnsetUsedHits(work.Pts[ipt]);
  } // ReleaseWorkHits

  //////////////////////////////////////////
  void TrajClusterAlg::SetAllHitsUsed(TrajPoint& tp)
  {
    // Sets UseHit true for all availalable hits in tp. Ensure
    // that hits on the opposite side of a used hit (belonging to
    // a finished trajectory) are not set used
    unsigned short ii, lo = 0, hi = tp.Hits.size();
    for(ii = 0; ii < tp.Hits.size(); ++ii) {
      if(inTraj[tp.Hits[ii]] > 0) {
        // found a used hit. See which side the TP position (hopefully
        // from a good previous fit) is. Let's hope that there is only
        // one used hit...
        float hitpos = fHits[tp.Hits[ii]]->PeakTime() * fScaleF;
        if(tp.Pos[1] < hitpos) {
          hi = ii;
        } else {
          lo = ii + 1;
        }
        break;
      } // inTraj[tp.Hits[ii]] > 0
    } // ii
    for(ii = lo; ii < hi; ++ii) {
      if(inTraj[tp.Hits[ii]] > 0) continue;
      tp.UseHit[ii] = true;
    } // ii
    DefineHitPos(tp);
  } // SetAllHitsUsed

  //////////////////////////////////////////
  void TrajClusterAlg::UnsetUsedHits(TrajPoint& tp)
  {
    // Sets inTraj = 0 and UseHit false for all used hits in tp
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(tp.UseHit[ii] && inTraj[tp.Hits[ii]] < 0) {
        inTraj[tp.Hits[ii]] = 0;
        tp.UseHit[ii] = false;
      } // UseHit
    } // ii
    tp.Chg = 0;
  } // UnsetUsedHits

  //////////////////////////////////////////
  void TrajClusterAlg::FindUseHits(Trajectory& tj, unsigned short ipt)
  {
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting inTraj < 0.
    
    if(ipt > tj.Pts.size() - 1) return;
    TrajPoint& tp = tj.Pts[ipt];
    
    if(tp.Hits.size() == 0) return;
    // some error checking
    if(tp.Hits.size() != tp.UseHit.size()) {
      mf::LogWarning("TC")<<"FindUseHits: Hits and UseHits vectors not the same size "<<tp.Hits.size()<<" "<<tp.UseHit.size();
      tp.Hits.clear();
      tp.UseHit.clear();
      return;
    }
    unsigned short ii;
    unsigned int iht;
    
    // Use everything (unused) for large angle TPs as long
    // as the multiplet doesn't include a hit used in another
    // trajectory
    if(IsLargeAngle(tp)) {
      bool useAll = true;
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        iht = tp.Hits[ii];
        if(inTraj[iht] > 0) useAll = false;
      } // ii
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        iht = tp.Hits[ii];
        if(inTraj[iht] > 0) continue;
        if(useAll) {
          tp.UseHit[ii] = true;
          inTraj[iht] = tj.ID;
        } else {
          tp.UseHit[ii] = false;
        }
      } // ii
      return;
    } // IsLargeAngle
    
    // Find the hit that has the smallest delta
    tp.Delta = 10;
    float delta;
    unsigned short imbest = USHRT_MAX;
    std::vector<float> deltas(tp.Hits.size(), 100);
    for(ii = 0; ii < tp.Hits.size(); ++ii) {
      tp.UseHit[ii] = false;
      iht = tp.Hits[ii];
      if(inTraj[iht] > 0) continue;
      delta = PointTrajDOCA(iht, tp);
      deltas[ii] = delta;
      if(delta < tp.Delta) {
        tp.Delta = delta;
        imbest = ii;
      }
    } // ii
    
    if(imbest == USHRT_MAX) return;
    iht = tp.Hits[imbest];
    
    // return victorious if the hit multiplicity is 1 and the charge is reasonable,
    // or if this is a high multiplicity hit
    if(fHits[iht]->Multiplicity() == 1 || fHits[iht]->Multiplicity() > 2) {
      if(fAveChg > 0) {
        float chgrat = (fHits[iht]->Integral() - fAveChg) / fAveChg;
        if(chgrat > 2) return;
      }
      tp.UseHit[imbest] = true;
      inTraj[iht] = tj.ID;
      return;
    }
    
    // Handle hit muliplets here. Construct a figure of merit for
    // the single hit, a hit doublet and for a hit multiplet
    // Calculate sfom
    float deltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(iht);
    float chgrat;
    deltaErr *= deltaErr;
    deltaErr += tp.DeltaRMS * tp.DeltaRMS;
    deltaErr = sqrt(deltaErr);
    float sfom = tp.Delta / deltaErr;
    if(prt) mf::LogVerbatim("TC")<<"  FindUseHits: single hit delta "<<tp.Delta<<" deltaErr "<<deltaErr<<" tp.DeltaRMS "<<tp.DeltaRMS<<" sfom "<<sfom;
    // Compare doublets using charge
    if(fAveChg > 0 && fHits[iht]->Multiplicity() == 2 && deltas.size() == 2) {
      // index of the other hit in tp.Hits
      unsigned short ombest = 1 - imbest;
      if(prt) mf::LogVerbatim("TC")<<"  doublet deltas "<<deltas[imbest]<<" "<<deltas[ombest];
      if(deltas[ombest] < 2 * deltas[imbest]) {
        unsigned int oht = tp.Hits[ombest];
        // make a temporary charge weighted sfom
        chgrat = std::abs(fHits[iht]->Integral() - fAveChg) / fAveChg;
        if(chgrat < 0.3) chgrat = 0.3;
        float scfom = sfom * chgrat;
        // make a fom for the other hit
        deltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(oht);
        deltaErr += tp.DeltaRMS * tp.DeltaRMS;
        deltaErr = sqrt(deltaErr);
        float ofom = deltas[ombest] / deltaErr;
        chgrat  = std::abs(fHits[oht]->Integral() - fAveChg) / fAveChg;
        if(chgrat < 0.3) chgrat = 0.3;
        ofom *= chgrat;
        if(prt) mf::LogVerbatim("TC")<<"  scfom "<<scfom<<"  ofom "<<ofom;
        if(ofom < scfom) {
          imbest = ombest;
          iht = oht;
        }
      } // deltas[ombest] < 2 * deltas[imbest]
    } // a doublet
    // Calculate mfom
    float fwire, ftime, dRms, qtot = 0;
    fwire = fHits[iht]->WireID().Wire;
    HitMultipletPosition(iht, ftime, dRms, qtot);
    float mfom = 100;
    if(qtot > 0) {
      ftime *= fScaleF;
      float mdelta = PointTrajDOCA(fwire, ftime, tp);
      deltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * dRms;
      deltaErr *= deltaErr;
      deltaErr += tp.DeltaRMS * tp.DeltaRMS;
      deltaErr = sqrt(deltaErr);
      mfom = mdelta / deltaErr;
      if(prt) mf::LogVerbatim("TC")<<"  FindUseHits: multi-hit delta "<<mdelta<<" deltaErr "<<deltaErr<<" tp.DeltaRMS "<<tp.DeltaRMS<<" mfom "<<mfom;
    }
    // Use the charge ratio if it is known
    if(fAveChg > 0 && qtot > 0) {
      // Weight by the charge difference
      chgrat = std::abs(fHits[iht]->Integral() - fAveChg) / fAveChg;
      // We expect the charge ratio to be within 30% so don't let a small
      // charge ratio dominate the decision
      if(chgrat < 0.3) chgrat = 0.3;
      sfom *= chgrat;
      if(prt) mf::LogVerbatim("TC")<<"  single hit chgrat "<<chgrat<<" sfom "<<sfom;
      chgrat = (qtot - fAveChg) /fAveChg;
      if(chgrat < 0.3) chgrat = 0.3;
      mfom *= chgrat;
      if(prt) mf::LogVerbatim("TC")<<"  multi-hit chgrat "<<chgrat<<" mfom "<<mfom;
    }
    // require the sfom to be significantly better than mfom
    // Use the multiplet if both are good
    if(mfom > 1 && mfom > 1.5 * sfom) {
      // use the single hit
      tp.UseHit[imbest] = true;
      inTraj[iht] = tj.ID;
    } else {
      unsigned int jht;
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        jht = tp.Hits[ii];
        if(inTraj[jht] > 0) continue;
//        if(fHits[jht]->StartTick() != fHits[iht]->StartTick()) continue;
        tp.UseHit[ii] = true;
        inTraj[jht] = tj.ID;
      } // ii
    } // use the multiplet
  } //  FindUseHits

  //////////////////////////////////////////
  void TrajClusterAlg::SetPoorUsedHits(Trajectory& tj, unsigned short ipt)
  {
    // Try to use the hits on this TP by reducing the number of points fitted. This
    // should only be done for reasonably long TJ's
    if(tj.Pts.size() < 2 * fMinPtsFit[fPass]) return;
    
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
    inTraj[iht] = tj.ID;
    unsigned short nptsFit = tp.NTPsFit;
    
    // Temp until we get this sorted out. Make sure that tj == work.
    if(tj.Pts.size() != work.Pts.size() || ipt != work.Pts.size() - 1) {
      mf::LogWarning("TC")<<"SetPoorUsedHits: tj != work. Fix this code";
      fQuitAlg = true;
      return;
    }

    // increment the number of points fit ala UpdateWork
    ++tp.NTPsFit;
    FitWork();
    if(tp.FitChi < 2) {
      // Looks like this will work. Return victorious after
      // decrementing NTPsFit
      --tp.NTPsFit;
      return;
    }
    while(tp.FitChi > 2 && tp.NTPsFit > fMinPtsFit[fPass]) {
      if(tp.NTPsFit > 15) {
        nptsFit = 0.7 * nptsFit;
      } else if(tp.NTPsFit > 4) {
        nptsFit -= 2;
      } else {
        nptsFit -= 1;
      }
      if(nptsFit < 3) nptsFit = 2;
      tp.NTPsFit = nptsFit;
      FitWork();
      if(prt) mf::LogVerbatim("TC")<<"  SetPoorUsedHits: FitChi "<<tp.FitChi<<" after reducing NTPsFit to "<<tp.NTPsFit;
      if(tp.NTPsFit <= fMinPtsFit[work.Pass]) break;
    } // reducing number of fitted points
    
    if(tp.FitChi > 2) {
      // Well that didn't work. Undo UseHit and return
      tp.UseHit[0] = false;
      inTraj[iht] = 0;
    }
    
  } // SetPoorUsedHits

  //////////////////////////////////////////
  bool TrajClusterAlg::HitChargeOK(Trajectory& tj, unsigned short ipt, unsigned short iht)
  {
    // returns true if the charge of the hit pointed to by ipt and iht is consistent
    // with the average charge fAveChg
    
    if(fAveChg == 0) return true;
    if(tj.Pts.size() < 3) return true;
    
    TrajPoint& tp = tj.Pts[ipt];
    unsigned int hit = tp.Hits[iht];
    float chgrat = (fHits[hit]->Integral() - fAveChg) / fAveChg;
    // Ignore this hit if it is really high charge
    if(chgrat > 3 * fChgRatCut) return false;
    // Accept if it is reasonable
    if(chgrat < fChgRatCut) return true;
    // Somewhat high charge. Check the charge of the previous point if it exists
    if(ipt == 0) return true;
    chgrat = (tj.Pts[ipt - 1].Chg - fAveChg) / fAveChg;
    // The previous point had high charge as well
    if(chgrat > fChgRatCut) return false;
    return true;
  } // HitChargeOK
  
  //////////////////////////////////////////
  void TrajClusterAlg::DefineHitPos(TrajPoint& tp)
  {
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    
    tp.HitPosErr2 = -1;
    if(tp.Hits.size() == 0) return;
    if(tp.Hits.size() != tp.UseHit.size()) {
      mf::LogWarning("TC")<<" Hits - UseHit size mis-match";
      fQuitAlg = true;
      return;
    }

    std::vector<unsigned int> hitVec;
    tp.Chg = 0;
    std::array<float, 2> newpos;
    newpos[0] = 0;
    newpos[1] = 0;
    unsigned short ii, iht;
    for(ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      iht = tp.Hits[ii];
      newpos[0] += fHits[iht]->Integral() * fHits[iht]->WireID().Wire;
      newpos[1] += fHits[iht]->Integral() * fHits[iht]->PeakTime() * fScaleF;
      tp.Chg += fHits[iht]->Integral();
      hitVec.push_back(iht);
    } // ii
 
    if(tp.Chg == 0) return;
    
    tp.HitPos[0] = newpos[0] / tp.Chg;
    tp.HitPos[1] = newpos[1] / tp.Chg;
    
    // Error is the wire error (2/sqrt(12))^2 and time error
    tp.HitPosErr2 = std::abs(tp.Dir[1]) * 0.16 + std::abs(tp.Dir[0]) * HitsTimeErr2(hitVec);

  } // HitPosErr2

  //////////////////////////////////////////
  float TrajClusterAlg::HitTimeErr(unsigned int iht)
  {
    return fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
  } // HitTimeErr
  
  //////////////////////////////////////////
  float TrajClusterAlg::HitsTimeErr2(std::vector<unsigned int> const& hitVec)
  {
    // Estimates the error^2 of the time using all hits in hitVec
    
    if(hitVec.size() == 0) return 0;
    float err;
    if(hitVec.size() == 1) {
      err = HitTimeErr(hitVec[0]);
      return err * err;
    } // hitVec.size() == 1
    
    // This approximation works for two hits of roughly similar RMS
    // and with roughly similar (within 2X) amplitude. TODO deal with the
    // case of more than 2 hits if the need arises
    float averms = 0.5 * (fHits[hitVec[0]]->RMS() + fHits[hitVec[1]]->RMS());
    float hitsep = (fHits[hitVec[0]]->PeakTime() - fHits[hitVec[1]]->PeakTime()) / averms;
    // This will estimate the RMS of two hits separated by 1 sigma by 1.5 * RMS of one hit
    err = averms * (1 + 0.14 * hitsep * hitsep) * fScaleF * fHitErrFac;
    // inflate this further for high multiplicity hits
    err *= fHits[hitVec[0]]->Multiplicity();
    return err * err;
    
  } // HitsTimeErr2
  
  //////////////////////////////////////////
  void TrajClusterAlg::SplitTrajCrossingVertices()
  {
    // This is kind of self-explanatory...
    
    if(vtx.size() == 0) return;
    if(allTraj.size() == 0) return;
    
    vtxPrt = (fDebugPlane == (int)fPlane && fDebugHit < 0);
    if(vtxPrt) mf::LogVerbatim("TC")<<"vtxPrt set for plane "<<fPlane<<" in SplitTrajCrossingVertices";
    if(vtxPrt) PrintAllTraj(USHRT_MAX, 999);
    
    // Splits trajectories in allTraj that cross a vertex
    unsigned short itj, iv, nTraj = allTraj.size();
    unsigned short tPass, closePt;
    float doca;
    for(itj = 0; itj < nTraj; ++itj) {
      // obsolete trajectory
      if(allTraj[itj].AlgMod[kKilled]) continue;
      tPass = allTraj[itj].Pass;
      for(iv = 0; iv < vtx.size(); ++iv) {
        // obsolete vertex
        if(vtx[iv].NTraj == 0) continue;
        // not in the cryostat/tpc/plane
        if(allTraj[itj].CTP != vtx[iv].CTP) continue;
        // already assigned to this vertex
//        if(allTraj[itj].Vtx[0] == iv) continue;
//        if(allTraj[itj].Vtx[1] == iv) continue;
        // too short
        if(allTraj[itj].EndPt[1] < 6) continue;
        TrajClosestApproach(allTraj[itj], vtx[iv].Wire, vtx[iv].Time, closePt, doca);
        if(vtxPrt)  mf::LogVerbatim("TC")<<" doca "<<doca<<" btw traj "<<itj<<" and vtx "<<iv<<" closePt "<<closePt<<" in plane "<<fPlane<<" CTP "<<vtx[iv].CTP;
        if(doca > fMaxVertexTrajSep[tPass]) continue;
        if(vtxPrt)  {
          mf::LogVerbatim("TC")<<"Good doca "<<doca<<" btw traj "<<itj<<" and vtx "<<iv<<" closePt "<<closePt<<" in plane "<<fPlane<<" CTP "<<vtx[iv].CTP;
          PrintTrajPoint(closePt, 1, tPass, allTraj[itj].Pts[closePt]);
        }
      } // iv
    } // itj
    
  } // SplitTrajCrossingVertices


  //////////////////////////////////////////
  void TrajClusterAlg::Find2DVertices()
  {
    
    if(fVertex2DIPCut <= 0) return;

    if(allTraj.size() < 2) return;
    
    vtxPrt = (fDebugPlane == (int)fPlane && fDebugHit < 0);
    if(vtxPrt) {
      mf::LogVerbatim("TC")<<"vtxPrt set for plane "<<fPlane<<" in Find2DVertices";
      PrintAllTraj(USHRT_MAX, allTraj.size());
    }
    
    unsigned short tj1, end1, endPt1, oendPt1, ivx;
    unsigned short tj2, end2, endPt2, oendPt2;
    unsigned short closePt1, closePt2;
    short dpt;
    float wint, tint, dw1, dt1, dw2, dt2, dang, doca;
    bool sigOK;
    float maxWire = fNumWires + 1;
    
    // Rough cut on dWire and dTime in WSE units
    float firstCut = 50;
    
    unsigned short tjSize = allTraj.size();

    for(tj1 = 0; tj1 < tjSize; ++tj1) {
      if(allTraj[tj1].AlgMod[kKilled]) continue;
      if(allTraj[tj1].CTP != fCTP) continue;
      if(allTraj[tj1].PDG != 13) continue;
      for(end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(allTraj[tj1].Vtx[end1] >= 0) continue;
        if(end1 == 0) {
          endPt1 = allTraj[tj1].EndPt[0];
          oendPt1 = allTraj[tj1].EndPt[1];
        } else {
          endPt1 = allTraj[tj1].EndPt[1];
          oendPt1 = allTraj[tj1].EndPt[0];
        }
        TrajPoint& tp1 = allTraj[tj1].Pts[endPt1];
        for(tj2 = tj1 + 1; tj2 < tjSize; ++tj2) {
          if(allTraj[tj2].AlgMod[kKilled]) continue;
          if(allTraj[tj2].CTP != fCTP) continue;
          if(allTraj[tj2].PDG != 13) continue;
          for(end2 = 0; end2 < 2; ++end2) {
            if(end2 == 0) {
              endPt2 = allTraj[tj2].EndPt[0];
              oendPt2 = allTraj[tj2].EndPt[1];
            } else {
              endPt2 = allTraj[tj2].EndPt[1];
              oendPt2 = allTraj[tj2].EndPt[0];
            }
            if(allTraj[tj2].Vtx[end2] >= 0) continue;
            TrajPoint& tp2 = allTraj[tj2].Pts[endPt2];
            TrajIntersection(tp1, tp2, wint, tint);
            // make sure this is inside the TPC
            if(wint < 0 || wint > maxWire) continue;
            if(tint < 0 || tint > fMaxTime) continue;
            if(vtxPrt) mf::LogVerbatim("TC")<<"Vtx candidate tj1-tj2 "<<allTraj[tj1].ID<<"_"<<end1<<"-"<<allTraj[tj2].ID<<"_"<<end2<<" vtx pos "<<(int)wint<<":"<<(int)(tint/fScaleF);
            dw1 = wint - tp1.Pos[0];
            if(std::abs(dw1) > firstCut) continue;
            dt1 = tint - tp1.Pos[1];
            if(std::abs(dt1) > firstCut) continue;
            dw2 = wint - tp2.Pos[0];
            if(std::abs(dw2) > firstCut) continue;
            dt2 = tint - tp2.Pos[1];
            if(std::abs(dt2) > firstCut) continue;
            // make sure that the other end isn't closer
            if(PointTrajDOCA2(wint, tint, tp1) > PointTrajDOCA2(wint, tint, allTraj[tj1].Pts[oendPt1])) continue;
            if(PointTrajDOCA2(wint, tint, tp2) > PointTrajDOCA2(wint, tint, allTraj[tj2].Pts[oendPt2])) continue;
            dang = std::abs(tp1.Ang - tp2.Ang);
            // A dang cut would go here
            // see if there is a wire signal at this point
            if(!SignalPresent(wint, tint, wint, tint) && end1 == end2) {
              // try nudging it by +/- 1 wire
              wint += 1;
              if(!SignalPresent(wint, tint, wint, tint)) {
                wint -= 2;
                if(!SignalPresent(wint, tint, wint, tint)) continue;
              }
            } // !SignalPresent
            // Ensure that the vertex position is close to an end
            TrajClosestApproach(allTraj[tj1], wint, tint, closePt1, doca);
            dpt = abs((short)endPt1 - (short)closePt1);
            sigOK = SignalPresent(wint, tint, allTraj[tj1].Pts[closePt1].Pos[0], allTraj[tj1].Pts[closePt1].Pos[1]);
            if(vtxPrt) mf::LogVerbatim("TC")<<" tj1 closest approach "<<doca<<" closePt1 "<<closePt1<<" dpt "<<dpt<<" signal? "<<sigOK;
            if(!sigOK) continue;
            if(allTraj[tj1].EndPt[1] > 4) {
              if(dpt > 3) continue;
            } else {
              // tighter cut for short trajectories
              if(dpt > 0) continue;
            }
            TrajClosestApproach(allTraj[tj2], wint, tint, closePt2, doca);
            dpt = abs((short)endPt2 - (short)closePt2);
            sigOK = SignalPresent(wint, tint, allTraj[tj2].Pts[closePt2].Pos[0], allTraj[tj2].Pts[closePt2].Pos[1]);
            if(vtxPrt) mf::LogVerbatim("TC")<<" tj2 closest approach "<<doca<<" closePt2 "<<closePt2<<" dpt "<<dpt<<" signal? "<<sigOK;
            if(!sigOK) continue;
            if(allTraj[tj2].EndPt[1] > 4) {
              if(dpt > 3) continue;
            } else {
              // tighter cut for short trajectories
              if(dpt > 0) continue;
            }
            if(vtxPrt) mf::LogVerbatim("TC")<<" wint:tint "<<(int)wint<<":"<<(int)tint<<" ticks "<<(int)(tint/fScaleF)<<" dang "<<dang;
            // Make sure this is a distinct vertex
            ivx = USHRT_MAX;
            for(unsigned short ii = 0; ii < vtx.size(); ++ii) {
              if(vtx[ii].NTraj == 0) continue;
              if(vtx[ii].CTP != fCTP) continue;
//              mf::LogVerbatim("TC")<<"chk "<<ivx<<" dw "<<std::abs(vtx[ii].Wire - wint)<<" dt "<<std::abs(vtx[ii].Time - tint);
              if(std::abs(vtx[ii].Wire - wint) > fVertex2DIPCut) continue;
              if(std::abs(vtx[ii].Time - tint) > fVertex2DIPCut) continue;
              ivx = ii;
              break;
            } // ivx
            if(ivx == USHRT_MAX) {
              // Found a new vertex
              VtxStore aVtx;
              aVtx.Wire = wint;
              aVtx.WireErr = 2;
              aVtx.Time = tint;
              aVtx.TimeErr = 2 * fScaleF;
              aVtx.NTraj = 0;
              aVtx.Topo = end1 + end2;
              aVtx.ChiDOF = 0;
              aVtx.CTP = fCTP;
              aVtx.Fixed = false;
              vtx.push_back(aVtx);
              ivx = vtx.size() - 1;
            }
            allTraj[tj1].Vtx[end1] = ivx;
            allTraj[tj2].Vtx[end2] = ivx;
            vtx[ivx].NTraj += 2;
            if(vtxPrt) mf::LogVerbatim("TC")<<" New vtx "<<ivx<<" CTP "<<fCTP<<" Pos "<<(int)wint<<":"<<(int)(tint/fScaleF)<<" ticks. NTraj "<<vtx[ivx].NTraj;
            // try to attach other TJs to this vertex if it new
            if(ivx == vtx.size()-1) AttachAnyTrajToVertex(ivx, fMaxVertexTrajSep[fPass], false);
          } // end2
        } // tj2
      } // end1
    } // tj1

    // Split trajectories that cross a vertex
    SplitTrajCrossingVertices();
    
    // Delete any vertex that has lost the required number of trajectories.
    // Update NClusters while we are here
    unsigned short tjVtx = 0, tjVtxEnd = 0;
    for(ivx = 0; ivx < vtx.size(); ++ivx) {
      if(vtx[ivx].CTP != fCTP) continue;
      vtx[ivx].NTraj = 0;
      for(tj1 = 0; tj1 < tjSize; ++tj1) {
        if(allTraj[tj1].AlgMod[kKilled]) continue;
        if(allTraj[tj1].CTP != fCTP) continue;
        for(end1 = 0; end1 < 2; ++end1) {
          if(allTraj[tj1].Vtx[end1] == ivx) {
            ++vtx[ivx].NTraj;
            tjVtx = tj1;
            tjVtxEnd = end1;
          }
        } // end1
      } // tj1
      // clobber any vertex that has only one trajectory
      if(vtx[ivx].NTraj == 1) {
        vtx[ivx].NTraj = 0;
        allTraj[tjVtx].Vtx[tjVtxEnd] = -1;
      }
    } // ivx

    if(vtxPrt) PrintAllTraj(USHRT_MAX, USHRT_MAX);

  } // Find2DVertices
  
  //////////////////////////////////////////
  void TrajClusterAlg::FindHammerVertices()
  {
    // Look for a trajectory that intersects another. Split
    // the trajectory and make a vertex. The convention used
    // is shown pictorially here
    // tj2  ------------------
    // tj1         /
    // tj1        /
    // tj1       /
    
    unsigned short itj1, end1, endPt1;
    unsigned short itj2, closePt2, ivx;
    unsigned short tjSize = allTraj.size();
    
    // minimum^2 DOCA of tj1 endpoint with tj2
    float doca, minDOCA = 3, dang;
    bool didaSplit;
    
    for(itj1 = 0; itj1 < tjSize; ++itj1) {
      if(allTraj[itj1].AlgMod[kKilled]) continue;
      // minimum length requirements
      if(NumPtsWithCharge(allTraj[itj1]) < 6) continue;
      // Track-like
      if(allTraj[itj1].PDG != 13) continue;
      // Check each end of tj1
      didaSplit = false;
      for(end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(allTraj[itj1].Vtx[end1] >= 0) continue;
        if(end1 == 0) {
          endPt1 = allTraj[itj1].EndPt[0];
        } else {
          endPt1 = allTraj[itj1].EndPt[1];
        }
        for(itj2 = 0; itj2 < tjSize; ++itj2) {
          if(itj1 == itj2) continue;
          Trajectory& tj2 = allTraj[itj2];
          if(tj2.AlgMod[kKilled]) continue;
          // require that both be in the same CTP
          if(tj2.CTP != allTraj[itj1].CTP) continue;
          // minimum length requirements
          if(NumPtsWithCharge(tj2) < 6) continue;
          // Require no vertex associated with itj2
          if(tj2.Vtx[0] >= 0 || tj2.Vtx[1] >= 0) continue;
          // track-like
          if(tj2.PDG != 13) continue;
          doca = minDOCA;
          TrajPointTrajDOCA(allTraj[itj1].Pts[endPt1], tj2, closePt2, doca);
          if(doca == minDOCA) continue;
          // ensure that the closest point is not near an end
          if(closePt2 < tj2.EndPt[0] + 3) continue;
          if(closePt2 > tj2.EndPt[1] - 3) continue;
          // make an angle cut
          dang = std::abs(allTraj[itj1].Pts[endPt1].Ang - tj2.Pts[closePt2].Ang);
          if(dang > M_PI/2) dang = M_PI - dang;
          if(vtxPrt) mf::LogVerbatim("TC")<<"FindHammerVertices: Candidate "<<allTraj[itj1].ID<<"  "<<tj2.ID<<" closePt2 "<<closePt2<<" doca "<<doca<<" dang "<<dang;
          if(dang < fKinkAngCut) continue;
          // we have a winner
          // create a new vertex
          VtxStore aVtx;
          aVtx.Wire = tj2.Pts[closePt2].Pos[0];
          aVtx.WireErr = 2;
          aVtx.Time = tj2.Pts[closePt2].Pos[1];
          aVtx.TimeErr = 2 * fScaleF;
          aVtx.NTraj = 3;
          aVtx.Topo = 6;
          aVtx.ChiDOF = 0;
          aVtx.CTP = fCTP;
          aVtx.Fixed = false;
          vtx.push_back(aVtx);
          ivx = vtx.size() - 1;
          SplitAllTraj(itj2, closePt2, ivx);
          if(!fSplitTrajOK) {
            if(vtxPrt) mf::LogVerbatim("TC")<<"FindHammerVertices: Failed to split trajectory";
            vtx.pop_back();
            continue;
          }
          allTraj[itj1].Vtx[end1] = ivx;
          allTraj[itj1].AlgMod[kHammer2DVx] = true;
          tj2.AlgMod[kHammer2DVx] = true;
          allTraj[allTraj.size()-1].AlgMod[kHammer2DVx] = true;
          didaSplit = true;
          break;
        } // tj2
        if(didaSplit) break;
      } // end1
    } // tj1
    
  } // FindHammerVertices

  //////////////////////////////////////////
  void TrajClusterAlg::AttachAnyTrajToVertex(unsigned short ivx, float docaCut, bool requireSignal)
  {
    // try to attach to existing vertices
    float doca, oldDoca;
    unsigned short itj, end, endPt, oendPt, closePt, oldVtx;
    short dpt;
    
    if(vtx[ivx].NTraj == 0) {
      mf::LogWarning("TC")<<"AttachAnyTrajToVertex: Trying to attach to an abandoned vertex";
      return;
    }
    unsigned short tjSize = allTraj.size();
    for(itj = 0; itj < tjSize; ++itj) {
      if(allTraj[itj].AlgMod[kKilled]) continue;
      if(allTraj[itj].CTP != vtx[ivx].CTP) continue;
      if(allTraj[itj].Vtx[0] == ivx || allTraj[itj].Vtx[1] == ivx) continue;
      for(end = 0; end < 2; ++end) {
        if(allTraj[itj].Vtx[end] >= 0) continue;
        endPt = allTraj[itj].EndPt[end];
        oendPt = allTraj[itj].EndPt[1-end];
        doca = PointTrajDOCA(vtx[ivx].Wire, vtx[ivx].Time, allTraj[itj].Pts[endPt]);
        if(doca > docaCut) continue;
        // ensure that the other end isn't closer
        if(doca > PointTrajDOCA(vtx[ivx].Wire, vtx[ivx].Time, allTraj[itj].Pts[oendPt])) continue;
        // See if the vertex position is close to an end
        TrajClosestApproach(allTraj[itj], vtx[ivx].Wire, vtx[ivx].Time, closePt, doca);
        dpt = abs((short)endPt - (short)closePt);
        if(vtxPrt) mf::LogVerbatim("TC")<<"AttachAnyTrajToVertex: tj "<<itj<<" closest approach "<<doca<<" closePt "<<closePt<<" endPt "<<endPt<<" oendPt "<<oendPt<<" dpt "<<dpt;
        if(doca > fVertex2DIPCut) continue;
// require a wire signal between these points? TODO write this code
//        if(requireSignal && !SignalPresent(vtx[ivx].Wire, vtx[ivx].Time, tj.Pts[oendPt])) continue;
        // See if the existing vertex assignment is better
        if(allTraj[itj].Vtx[end] >= 0) {
          oldVtx = allTraj[itj].Vtx[end];
          oldDoca = PointTrajDOCA(vtx[oldVtx].Wire, vtx[oldVtx].Time, allTraj[itj].Pts[endPt]);
          if(oldDoca < doca) continue;
        }
        if(dpt < 3) {
          allTraj[itj].Vtx[end] = ivx;
          ++vtx[ivx].NTraj;
          if(vtxPrt) mf::LogVerbatim("TC")<<" Attach tj "<<itj<<" to vtx "<<ivx;
        } else {
          // split the trajectory at this point A BAD IDEA
//          SplitAllTraj(itj, closePt, ivx);
        }
      } // end
    } // itj
  } // AttachTrajToAVertex

  
  //////////////////////////////////////
  void TrajClusterAlg::Find3DVertices(geo::TPCID const& tpcid)
  {
    // Create 3D vertices from 2D vertices. 3D vertices that are matched
    // in all three planes have Ptr2D >= 0 for all planes
    
    
    geo::TPCGeo const& TPC = geom->TPC(tpcid);
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    const unsigned int cstat = tpcid.Cryostat;
    const unsigned int tpc = tpcid.TPC;
    
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
    
    vtxPrt = (fDebugPlane >= 0) && (fDebugHit == 6666);
    
    if(vtxPrt) {
      mf::LogVerbatim("CC")<<"Inside Find3DVertices";
      PrintAllTraj(USHRT_MAX, allTraj.size());
    }
    
    // wire spacing in cm
    float wirePitch = geom->WirePitch(0, 1, 0, tpcid.TPC, tpcid.Cryostat);
    
    // fill temp vectors of 2D vertex X and X errors
    std::vector<float> vX(vtx.size());
    std::vector<float> vXsigma(vtx.size());
    float vXp;
    double ticks;
    for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if(vtx[ivx].NTraj == 0) continue;
      geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
      if(iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
      // Convert 2D vertex time error to X error
      ticks = vtx[ivx].Time / fScaleF;
      vX[ivx]  = detprop->ConvertTicksToX(ticks, (int)iplID.Plane, (int)tpc, (int)cstat);
      ticks = (vtx[ivx].Time + vtx[ivx].TimeErr) / fScaleF;
      vXp = detprop->ConvertTicksToX(ticks, (int)iplID.Plane, (int)tpc, (int)cstat);
      vXsigma[ivx] = fabs(vXp - vX[ivx]);
    } // ivx
    
    // create a array/vector of 2D vertex indices in each plane
    std::array<std::vector<unsigned short>, 3> vIndex;
    unsigned short indx, ipl;
    for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if(vtx[ivx].NTraj == 0) continue;
      geo::PlaneID iplID = DecodeCTP(vtx[ivx].CTP);
      if(iplID.TPC != tpc || iplID.Cryostat != cstat) continue;
      ipl = iplID.Plane;
      if(ipl > 2) continue;
      indx = vIndex[ipl].size();
      vIndex[ipl].resize(indx + 1);
      vIndex[ipl][indx] = ivx;
    }
    
    // vector of 2D vertices -> 3D vertices.
    std::vector<short> vPtr;
    for(unsigned short ii = 0; ii < vtx.size(); ++ii) vPtr.push_back(-1);
    
    // temp vector of all 2D vertex matches
    std::vector<Vtx3Store> v3temp;
    
    double y = 0, z = 0;
    TVector3 WPos = {0, 0, 0};
    // i, j, k indicates 3 different wire planes
    unsigned short ii, jpl, jj, kpl, kk, ivx, jvx, kvx;
    unsigned int iWire, jWire;
    unsigned short v3dBest = 0;
    float xbest = 0, ybest = 0, zbest = 0;
    float kX, kWire;
    // compare vertices in each view
    bool gotit = false;
    for(ipl = 0; ipl < 2; ++ipl) {
      for(ii = 0; ii < vIndex[ipl].size(); ++ii) {
        ivx = vIndex[ipl][ii];
        if(ivx > vtx.size() - 1) {
          mf::LogError("CC")<<"Find3DVertices: bad ivx "<<ivx;
          return;
        }
        // vertex has been matched already
        if(vPtr[ivx] >= 0) continue;
        iWire = vtx[ivx].Wire;
        float best = fVertex3DChiCut;
        // temp array of 2D vertex indices in each plane
        std::array<short, 3> t2dIndex = {-1, -1, -1};
        std::array<short, 3> tmpIndex = {-1, -1, -1};
        for(jpl = ipl + 1; jpl < 3; ++jpl) {
          for(jj = 0; jj < vIndex[jpl].size(); ++jj) {
            jvx = vIndex[jpl][jj];
            if(jvx > vtx.size() - 1) {
              mf::LogError("CC")<<"Find3DVertices: bad jvx "<<jvx;
              return;
            }
            // vertex has been matched already
            if(vPtr[jvx] >= 0) continue;
            jWire = vtx[jvx].Wire;
            // new stuff
            float dX = fabs(vX[ivx] - vX[jvx]);
            float dXSigma = sqrt(vXsigma[ivx] * vXsigma[ivx] + vXsigma[jvx] * vXsigma[jvx]);
            float dXChi = dX / dXSigma;
            
            if(vtxPrt) mf::LogVerbatim("CC")<<"Find3DVertices: ipl "<<ipl<<" ivx "<<ivx<<" ivX "<<vX[ivx]
              <<" jpl "<<jpl<<" jvx "<<jvx<<" jvX "<<vX[jvx]<<" W:T "<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time<<" dXChi "<<dXChi<<" fVertex3DChiCut "<<fVertex3DChiCut;
            
            if(dXChi > fVertex3DChiCut) continue;
            geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
            if(y < YLo || y > YHi || z < ZLo || z > ZHi) continue;
            WPos[1] = y;
            WPos[2] = z;
            kpl = 3 - ipl - jpl;
            kX = 0.5 * (vX[ivx] + vX[jvx]);
            kWire = -1;
            if(TPC.Nplanes() > 2) kWire = geom->NearestWire(WPos, kpl, tpc, cstat);
            kpl = 3 - ipl - jpl;
            // save this incomplete 3D vertex
            Vtx3Store v3d;
            v3d.ProcCode = 1;
            tmpIndex[ipl] = ivx;
            tmpIndex[jpl] = jvx;
            tmpIndex[kpl] = -1;
            v3d.Ptr2D = tmpIndex;
            v3d.X = kX;
            v3d.XErr = dXSigma;
            v3d.Y = y;
            float yzSigma = wirePitch * sqrt(vtx[ivx].WireErr * vtx[ivx].WireErr + vtx[jvx].WireErr * vtx[jvx].WireErr);
            v3d.YErr = yzSigma;
            v3d.Z = z;
            v3d.ZErr = yzSigma;
            v3d.Wire = kWire;
            v3d.CStat = cstat;
            v3d.TPC = tpc;
            v3temp.push_back(v3d);
            
            if(vtxPrt) mf::LogVerbatim("CC")<<"Find3DVertices: 2 Plane match ivx "<<ivx<<" P:W:T "<<ipl<<":"<<(int)vtx[ivx].Wire<<":"<<(int)vtx[ivx].Time<<" jvx "<<jvx<<" P:W:T "<<jpl<<":"<<(int)vtx[jvx].Wire<<":"<<(int)vtx[jvx].Time<<" dXChi "<<dXChi<<" yzSigma "<<yzSigma;
            
            if(TPC.Nplanes() == 2) continue;
            // look for a 3 plane match
            best = fVertex3DChiCut;
            for(kk = 0; kk < vIndex[kpl].size(); ++kk) {
              kvx = vIndex[kpl][kk];
              if(vPtr[kvx] >= 0) continue;
              // Wire difference error
              float dW = wirePitch * (vtx[kvx].Wire - kWire) / yzSigma;
              // X difference error
              float dX = (vX[kvx] - kX) / dXSigma;
              float kChi = 0.5 * sqrt(dW * dW + dX * dX);
              if(kChi < best) {
                best = kChi;
                xbest = (vX[kvx] + 2 * kX) / 3;
                ybest = y;
                zbest = z;
                t2dIndex[ipl] = ivx;
                t2dIndex[jpl] = jvx;
                t2dIndex[kpl] = kvx;
                v3dBest = v3temp.size() - 1;
              }
              
              if(vtxPrt) mf::LogVerbatim("CC")<<" kvx "<<kvx<<" kpl "<<kpl
                <<" wire "<<(int)vtx[kvx].Wire<<" kTime "<<(int)vtx[kvx].Time<<" kChi "<<kChi<<" best "<<best<<" dW "<<vtx[kvx].Wire - kWire;
              
            } // kk
            if(vtxPrt) mf::LogVerbatim("CC")<<" done best = "<<best<<" fVertex3DChiCut "<<fVertex3DChiCut;
            if(TPC.Nplanes() > 2 && best < fVertex3DChiCut) {
              // create a real 3D vertex using the previously entered incomplete 3D vertex as a template
              if(v3dBest > v3temp.size() - 1) {
                mf::LogError("CC")<<"Find3DVertices: bad v3dBest "<<v3dBest;
                return;
              }
              Vtx3Store v3d = v3temp[v3dBest];
              v3d.Ptr2D = t2dIndex;
              v3d.Wire = -1;
              // TODO need to average ybest and zbest here with error weighting
              v3d.X = xbest;
              v3d.Y = ybest;
              v3d.Z = zbest;
              vtx3.push_back(v3d);
              gotit = true;
              // mark the 2D vertices as used
              for(unsigned short jj = 0; jj < 3; ++jj) if(t2dIndex[jj] >= 0) vPtr[t2dIndex[jj]] = vtx3.size() - 1;
              
              if(vtxPrt) mf::LogVerbatim("CC")<<"New 3D vtx "<<vtx3.size()<<" X "<<v3d.X<<" Y "<<v3d.Y<<" Z "<<v3d.Z
                <<" t2dIndex "<<t2dIndex[0]<<" "<<t2dIndex[1]<<" "<<t2dIndex[2]<<" best Chi "<<best;
              
            } // best < dRCut
            if(gotit) break;
          } // jj
          if(gotit) break;
        } // jpl
        if(gotit) break;
      } // ii
    } // ipl
    
    // Store incomplete 3D vertices but ignore those that are part of a complete 3D vertex
    unsigned short vsize = vtx3.size();
    unsigned short ninc = 0;
    for(unsigned short it = 0; it < v3temp.size(); ++it) {
      bool keepit = true;
      for(unsigned short i3d = 0; i3d < vsize; ++i3d) {
        for(unsigned short plane = 0; plane < 3; ++plane) {
          if(v3temp[it].Ptr2D[plane] == vtx3[i3d].Ptr2D[plane]) {
            keepit = false;
            break;
          }
        } // plane
        if(!keepit) break;
      } // i3d
      if(keepit) {
        vtx3.push_back(v3temp[it]);
        ++ninc;
      }
    } // it
    
    // Modify Ptr2D for 2-plane detector
    if(TPC.Nplanes() == 2) {
      for(unsigned short iv3 = 0; iv3 < vtx3.size(); ++iv3) {
        vtx3[iv3].Ptr2D[2] = 666;
      } //iv3
    } // 2 planes
    
    if(vtxPrt) {
       for(unsigned short it = 0; it < vtx3.size(); ++it) {
        mf::LogVerbatim("CC")<<"vtx3 "<<it<<" Ptr2D "<<vtx3[it].Ptr2D[0]<<" "<<vtx3[it].Ptr2D[1]<<" "<<vtx3[it].Ptr2D[2]
        <<" wire "<<vtx3[it].Wire;
      }
    }
    
    // Try to complete incomplete vertices
    if(ninc > 0) CompleteIncomplete3DVertices(tpcid);
    
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
    for(auto& vx3 : vtx3) {
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
      tp.Pos[1] = detprop->ConvertXToTicks(vx3.X, mPlane, vx3.TPC, vx3.CStat) * fScaleF;
      for(itj = 0; itj < allTraj.size(); ++itj) {
        if(allTraj[itj].CTP != mCTP) continue;
        if(allTraj[itj].AlgMod[kKilled]) continue;
        if(allTraj[itj].PDG != 13) continue;
        doca = maxdoca;
        // find the closest distance between the vertex and the trajectory
        TrajPointTrajDOCA(tp, allTraj[itj], closePt, doca);
        mf::LogVerbatim("TC")<<"CI3DV itj ID "<<allTraj[itj].ID<<" closePT "<<closePt<<" doca "<<doca;
        if(doca == maxdoca) continue;
        mTjs.push_back(std::make_pair(itj, closePt));
      } // itj
      // handle the case where there are one or more TJs with TPs near the ends
      // that make a vertex (a failure by Find2DVertices)
      if(mTjs.size() == 0) continue;
      aVtxIndx = vtx.size();
      aVtx.CTP = mCTP;
      aVtx.Topo = 6;
      aVtx.NTraj = mTjs.size();
      aVtx.Wire = tp.Pos[0];
      aVtx.Time = tp.Pos[1];
      nEndPt = 0;
      for(ii = 0; ii < mTjs.size(); ++ii) {
        itj = mTjs[ii].first;
        closePt = mTjs[ii].second;
        // determine which end is the closest
        if(fabs(closePt - allTraj[itj].EndPt[0]) < fabs(closePt - allTraj[itj].EndPt[1])) {
          // closest to the beginning
          end = 0;
        } else {
          // closest to the end
          end = 1;
        }// closest to the end
        dpt = fabs(closePt - allTraj[itj].EndPt[end]);
        if(dpt < 4) {
          // close to an end
          allTraj[itj].Vtx[end] = aVtxIndx;
          allTraj[itj].AlgMod[kComp3DVx] = true;
          ++nEndPt;
        } else {
          // closePt is not near an end, so split the trajectory
          SplitAllTraj(itj, closePt, aVtxIndx);
          if(!fSplitTrajOK) {
            mf::LogVerbatim("TC")<<"SplitAllTraj failed ";
            // failure occurred. Recover
            for(auto mTj : mTjs) {
              unsigned short jtj = mTj.first;
              if(allTraj[jtj].Vtx[0] == aVtxIndx) allTraj[jtj].Vtx[0] = -1;
              if(allTraj[jtj].Vtx[1] == aVtxIndx) allTraj[jtj].Vtx[1] = -1;
            } // itj
            continue;
          } // !fSplitTrajOK
        } // closePt is not near an end, so split the trajectory
        allTraj[itj].AlgMod[kComp3DVx] = true;
        itj = allTraj.size() - 1;
        allTraj[itj].AlgMod[kComp3DVx] = true;
      } // ii
      vtx.push_back(aVtx);
      vx3.Ptr2D[mPlane] = aVtxIndx;
      vx3.Wire = -1;
       mf::LogVerbatim("TC")<<"CompleteIncomplete3DVertices: new 2D vtx "<<aVtxIndx<<" points to 3D vtx ";
    } // vx3
  } // CompleteIncomplete3DVertices
  
  //////////////////////////////////////////
  short TrajClusterAlg::TPNearVertex(const TrajPoint& tp)
  {
    // Returns the index of a vertex if tp is heading towards it
    for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
      if(vtx[ivx].NTraj == 0) continue;
      if(vtx[ivx].CTP != tp.CTP) continue;
      if(std::abs(vtx[ivx].Wire - tp.Pos[0]) > 3) continue;
      if(std::abs(vtx[ivx].Time - tp.Pos[1]) > 3) continue;
      if(PointTrajDOCA2(vtx[ivx].Wire, vtx[ivx].Time, tp) > 3) continue;
      // see if the TP points to the vertex and not away from it
      if(tp.Pos[0] < vtx[ivx].Wire && tp.Dir[0] > 0) return (short)ivx;
      if(tp.Pos[0] > vtx[ivx].Wire && tp.Dir[0] < 0) return (short)ivx;
    } // ivx
    return -1;
  } // TPNearVertex


  //////////////////////////////////////////
  void TrajClusterAlg::StepCrawl()
  {
    // Crawl along the direction specified in the traj vector in steps of size step
    // (wire spacing equivalents). Find hits between the last trajectory point and
    // the last trajectory point + step. A new trajectory point is added if hits are
    // found. Crawling continues until no signal is found for two consecutive steps
    // or until a wire or time boundary is reached.
    
    fGoodWork = false;
    fTryWithNextPass = false;
    if(work.Pts.size() == 0) return;
 
    unsigned short iwt, lastPt;
    unsigned short lastPtWithHits = work.EndPt[1];
    if(lastPtWithHits != work.Pts.size() - 1) {
      mf::LogWarning("TC")<<"StepCrawl: Starting trajectory has no hits on the leading edge. Assume this is an error and quit.";
      PrintTrajectory(work, USHRT_MAX);
      return;
    }
    lastPt = lastPtWithHits;
    // Construct a local TP from the last work TP that will be moved on each step.
    // Only the Pos and Dir variables will be used
    TrajPoint ltp;
    ltp.Pos = work.Pts[lastPt].Pos;
    ltp.Dir = work.Pts[lastPt].Dir;
    // A second TP is cloned from the leading TP of work, updated with hits, fit
    // parameters,etc and possibly pushed onto work as the next TP
    TrajPoint tp ;
    
    unsigned int step;
    float stepSize;
    float nMissedWires;
    unsigned short nNoSignalSteps = 0;
    
    // decide on the number of steps expected without adding a hit
    unsigned short missedStepCut = SetMissedStepCut(ltp);
    
    if(prt) mf::LogVerbatim("TC")<<"Start StepCrawl with missedStepCut "<<missedStepCut<<" fScaleF "<<fScaleF;

    bool sigOK, keepGoing;
    unsigned short killPts;
    bool isLA;
    for(step = 1; step < 10000; ++step) {
      isLA = IsLargeAngle(ltp);
      if(isLA) { stepSize = 1; } else { stepSize = std::abs(1/ltp.Dir[0]); }
      // make a copy of the previous TP
      lastPt = work.Pts.size() - 1;
      tp = work.Pts[lastPt];
      ++tp.Step;
      // move the local TP position by one step in the right direction
      for(iwt = 0; iwt < 2; ++iwt) ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;
      if(TPNearVertex(ltp) >= 0) {
        work.AlgMod[kStopAtVtx] = true;
        break;
      }
      // copy this position into tp
      tp.Pos = ltp.Pos;
      tp.Dir = ltp.Dir;
      if(prt) mf::LogVerbatim("TC")<<"StepCrawl "<<step<<" Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" Dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" stepSize "<<stepSize;
      // hit the boundary of the TPC?
      if(tp.Pos[0] < 0 || tp.Pos[0] > fMaxWire) break;
      if(tp.Pos[1] < 0 || tp.Pos[1] > fMaxTime) break;
      // remove the old hits and other stuff
      tp.Hits.clear();
      tp.UseHit.clear();
      tp.FitChi = 0; tp.Chg = 0;
      // append to the work trajectory
      work.Pts.push_back(tp);
      // update the index of the last TP
      lastPt = work.Pts.size() - 1;
      // look for hits
      AddHits(work, lastPt, sigOK);
      // If successfull, AddHits has defined UseHit for this TP,
      // set the trajectory endpoints, and defined HitPos
      // overwrite the local TP with the updated work TP. The trajectory
      // fit and average values will be updated further down by UpdateWork
      if(work.Pts[lastPt].Hits.size() == 0) {
        // No close hits added.
        // Check not-large angle trajectories
        if(!SignalAtTp(ltp)) {
          ++nNoSignalSteps;
          if(isLA) {
            if(nNoSignalSteps > 20) {
              if(prt) mf::LogVerbatim("TC")<<" StepCrawl: no signal at ltp "<<PrintPos(ltp)<<" after "<<nNoSignalSteps<<" max = 20";
              break;
            }
          } else {
            if(nNoSignalSteps > fMaxWireSkipNoSignal) {
              if(prt) mf::LogVerbatim("TC")<<" StepCrawl: no signal at ltp "<<PrintPos(ltp)<<" after "<<nNoSignalSteps<<" max = "<<fMaxWireSkipNoSignal;
              break;
            }
          }
          // no sense keeping this TP on work if no hits were added
          work.Pts.pop_back();
          continue;
        }
        lastPtWithHits = work.EndPt[1];
        // TODO: TrajPointSeparation measures the distance along the trajectory but DeadWireCount
        // is the distance in the wire direction...
        nMissedWires = TrajPointSeparation(work.Pts[lastPtWithHits], ltp) - DeadWireCount(ltp.Pos[0], work.Pts[lastPtWithHits].Pos[0]);
        if(prt) mf::LogVerbatim("TC")<<" StepCrawl: nMissedWires "<<nMissedWires<<" user cut "<<fMaxWireSkipNoSignal;
        if(nMissedWires > fMaxWireSkipNoSignal) break;
        // no sense keeping this TP on work if no hits were added
        work.Pts.pop_back();
        continue;
      } // work.Pts[lastPt].Hits.size() == 0
      nNoSignalSteps = 0;
      if(work.Pts[lastPt].Chg == 0) {
        // There are points on the trajectory by none used in the last step. See
        // how long this has been going on
        lastPtWithHits = work.EndPt[1];
        nMissedWires = std::abs(ltp.Pos[0] - work.Pts[lastPtWithHits].Pos[0])
        - DeadWireCount(ltp.Pos[0], work.Pts[lastPtWithHits].Pos[0]);
        if(prt)  mf::LogVerbatim("TC")<<" Hits exist on the trajectory but are not used. Missed wires "<<nMissedWires;
        if(nMissedWires > 0.5 * work.Pts.size()) {
          if(MaybeDeltaRay(work)) {
            if(prt) mf::LogVerbatim("TC")<<" Tagged as a possible delta ray";
            work.PDG = 12;
            break;
          }
        } // nMissedWires > 0.5 * work.Pts.size()
        // Keep stepping
        if(prt) PrintTrajectory(work, lastPt);
        continue;
      } // tp.Hits.size() == 0
      // Update the last point fit, etc using the just added hit(s)
      UpdateTraj(work);
      // a failure occurred
      if(!fUpdateTrajOK) return;
      // Quit if we are starting out poorly. This can happen on the 2nd trajectory point
      // when a hit is picked up in the wrong direction
      if(work.Pts.size() == 2 && std::signbit(work.Pts[1].Dir[0]) != std::signbit(work.Pts[0].Dir[0])) {
        if(prt) mf::LogVerbatim("TC")<<" Picked wrong hit on 2nd traj point. Drop trajectory.";
        return;
      }
      if(work.Pts.size() == 3) {
        // ensure that the last hit added is in the same direction as the first two.
        // This is a simple way of doing it
        if(prt) mf::LogVerbatim("TC")<<"StepCrawl: Third TP. 0-2 sep "<<TrajPointHitSep2(work.Pts[0], work.Pts[2])<<" 0-1 sep "<<TrajPointHitSep2(work.Pts[0], work.Pts[1]);
        if(TrajPointHitSep2(work.Pts[0], work.Pts[2]) < TrajPointHitSep2(work.Pts[0], work.Pts[1])) return;
      } // work.Pts.size() == 3
      // Update the local TP with the updated position and direction
      ltp.Pos = work.Pts[lastPt].Pos;
      ltp.Dir = work.Pts[lastPt].Dir;
      if(fMaskedLastTP) {
        // see if TPs have been masked off many times and if the
        // environment is clean. If so, return and try with next pass
        // cuts
        if(!MaskedWorkHitsOK()) {
          if(prt) PrintTrajectory(work, lastPt);
          return;
        }
        // Don't bother with the rest of the checking below if we
        // set all hits not used on this TP
        if(prt) PrintTrajectory(work, lastPt);
        continue;
      }
      // Try to use all hits in doublets
      // TODO something wrong in this code
//      ModifyShortTraj(work);
      // We have added a TP with hits
      // assume that we aren't going to kill the point we just added, or any
      // of the previous points...
      killPts = 0;
      // assume that we should keep going after killing points
      keepGoing = true;
      // check for a kink. Stop crawling if one is found
      GottaKink(work, killPts);
      if(killPts > 0) keepGoing = false;
/*
      // testing
      if(work.Pts.size() == 6 && MaybeDeltaRay(work)) {
        mf::LogVerbatim("TC")<<"Might be\n";
        PrintTrajectory(work, USHRT_MAX);
       }
*/
      // See if the Chisq/DOF exceeds the maximum.
      // UpdateTraj should have reduced the number of points fit
      // as much as possible for this pass, so this trajectory is in trouble.
      if(killPts == 0 &&  work.Pts[lastPt].FitChi > fMaxChi) {
        if(prt) mf::LogVerbatim("TC")<<"   bad FitChi "<<work.Pts[lastPt].FitChi<<" cut "<<fMaxChi;
        fGoodWork = (NumPtsWithCharge(work) > fMinPtsFit[fPass]);
        return;
      }
      // print the local tp unless we have killing to do
      if(killPts == 0) {
        if(prt) PrintTrajectory(work, lastPt);
      } else {
        MaskTrajEndPoints(work, killPts);
/*
        // Kill some trajectory points. Keep track of the wire we are on at the current TP.
        // The trajectory has been corrupted by these hits, so we need to first remove them
        // and then recalculate the trajectory position of the current step.
        onWire = (float)(std::nearbyint(work.Pts[lastPt].Pos[0]));
        // don't remove points but simply set UseHit false
        unsigned short ii, indx, iht;
        for(unsigned short kill = 0; kill < killPts; ++kill) {
          indx = work.Pts.size() - 1 - kill;
          if(prt) mf::LogVerbatim("TC")<<"TRP   killing hits in TP "<<indx;
          for(ii = 0; ii < work.Pts[indx].Hits.size(); ++ii) {
            iht = work.Pts[indx].Hits[ii];
            work.Pts[indx].UseHit[ii] = false;
            if(inTraj[iht] == work.ID) inTraj[iht] = 0;
          } // ii
        } // kill
*/
        unsigned int onWire = (float)(std::nearbyint(work.Pts[lastPt].Pos[0]));
        float nSteps = (float)(step - work.Pts[lastPt - killPts].Step);
//        if(prt) mf::LogVerbatim("TC")<<"TRP   killing "<<killPts<<" after "<<nSteps<<" steps from prev TP.  Current tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1];
        // move the position
        work.Pts[lastPt].Pos[0] += nSteps * work.Pts[lastPt].Dir[0];
        work.Pts[lastPt].Pos[1] += nSteps * work.Pts[lastPt].Dir[1];
        if(!isLA) {
          // put the TP at the wire position prior to the move
          float dw = onWire - work.Pts[lastPt].Pos[0];
          work.Pts[lastPt].Pos[0] = onWire;
          work.Pts[lastPt].Pos[1] += dw * work.Pts[lastPt].Dir[1] / work.Pts[lastPt].Dir[0];
        }
        // copy to the local trajectory point
        ltp.Pos = work.Pts[lastPt].Pos;
        ltp.Dir = work.Pts[lastPt].Dir;
        if(prt) mf::LogVerbatim("TC")<<"  New ltp.Pos     "<<ltp.Pos[0]<<" "<<ltp.Pos[1]<<" ticks "<<(int)ltp.Pos[1]/fScaleF;
        if(!keepGoing) break;
      }
    } // step
    
    fGoodWork = true;
    if(prt) mf::LogVerbatim("TC")<<"End StepCrawl with "<<step<<" steps. work size "<<work.Pts.size()<<" fGoodWork = "<<fGoodWork<<" with fTryWithNextPass "<<fTryWithNextPass;

    if(fGoodWork && fTryWithNextPass) {
      mf::LogVerbatim("TC")<<"StepCrawl: Have fGoodWork && fTryWithNextPass true. This shouldn't happen. Fixing it.";
      fTryWithNextPass = false;
    }

  } // StepCrawl
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::MaskTrajEndPoints(Trajectory& tj, unsigned short nPts)
  {
    // Mask off some trajectory points at the leading edge. Keep track of the wire we are on at the current TP.
    // The trajectory has been corrupted by these hits, so we need to first remove them
    // and then recalculate the trajectory position of the current step.
    // don't remove points but simply set UseHit false
    unsigned short ii, ipt, iht;
    for(unsigned short kill = 0; kill < nPts; ++kill) {
      ipt = tj.Pts.size() - 1 - kill;
      if(prt) mf::LogVerbatim("TC")<<"TRP   masking hits in TP "<<ipt;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        iht = tj.Pts[ipt].Hits[ii];
        tj.Pts[ipt].UseHit[ii] = false;
        if(inTraj[iht] == tj.ID) inTraj[iht] = 0;
      } // ii
    } // kill
    
  } // MaskTrajEndPoints

  ////////////////////////////////////////////////
  void TrajClusterAlg::ModifyShortTraj(Trajectory& tj)
  {
    // See if most of the hits in the short trajectory are doublets and
    // we are only using single hits in each TP. If so, see if we would do
    // better by using all of the hits

    if(tj.Pts.size() != 4) return;
    unsigned short ndouble = 0, nsingle = 0;
    for(auto& tp : tj.Pts) {
      if(tp.Hits.size() == 2) {
        ++ndouble;
        if(NumUsedHits(tp) == 1) ++nsingle;
      }
    } // tp
    
    if(ndouble < 2) return;
    if(nsingle < 2) return;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"Inside ModifyShortTraj. ndouble "<<ndouble<<" nsingle "<<nsingle<<" Starting trajectory";
      PrintTrajectory(tj, USHRT_MAX);
    }
    
    // clone the trajectory, use all the hits in doublets and check the fit
    Trajectory ntj = tj;
    unsigned short ii;
    unsigned int iht;
    std::vector<unsigned int> inTrajHits;
    for(auto& tp : ntj.Pts) {
      if(tp.Hits.size() == 2 && NumUsedHits(tp) == 1) {
        for(ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tp.UseHit[ii]) continue;
          tp.UseHit[ii] = true;
          iht = tp.Hits[ii];
          inTraj[iht] = ntj.ID;
          // save the hit indices in case we need to undo this
          inTrajHits.push_back(iht);
          mf::LogVerbatim("TC")<<"Use hit "<<PrintHit(iht);
        } // ii
        mf::LogVerbatim("TC")<<" HitPos1 "<<tp.HitPos[1];
        DefineHitPos(tp);
        mf::LogVerbatim("TC")<<" new HitPos1 "<<tp.HitPos[1];
      } // tp.Hits.size() == 2 && NumUsedHits(tp) == 1
    } // tp
    // do the fit
    unsigned short originPt = ntj.Pts.size() - 1;
    unsigned short npts = ntj.Pts[originPt].NTPsFit;
    TrajPoint tpFit;
    unsigned short fitDir = -1;
    FitTraj(ntj, originPt, npts, fitDir, tpFit);
    if(prt) {
      mf::LogVerbatim("TC")<<"tpFit trajectory";
      PrintTrajectory(ntj, USHRT_MAX);
      mf::LogVerbatim("TC")<<"compare "<<tj.Pts[originPt].FitChi<<" "<<ntj.Pts[originPt].FitChi<<"\n";
    }

    if(ntj.Pts[originPt].FitChi < 1.2 * tj.Pts[originPt].FitChi) {
      // Use the modified trajectory
      tj = ntj;
      tj.AlgMod[kModifyShortTraj] = true;
      if(prt) mf::LogVerbatim("TC")<<"ModifyShortTraj: use modified trajectory";
    } else {
      // Use the original trajectory
      for(auto iht : inTrajHits) inTraj[iht] = 0;
      if(prt) mf::LogVerbatim("TC")<<"ModifyShortTraj: use original trajectory";
    }
    
  } // ModifyShortTraj

  ////////////////////////////////////////////////
  bool TrajClusterAlg::MaybeDeltaRay(Trajectory& tj)
  {
    // See if the trajectory appears to be a delta ray. This is characterized by a significant fraction of hits
    // in the trajectory belonging to an existing trajectory. This may also flag ghost trajectories...
    
    // vectors of traj IDs, and the occurrence count
    std::vector<unsigned short> tID, tCnt;
    unsigned short itj, indx;
    unsigned short tCut = 0.5 * tj.Pts.size();
    for(auto& tp : tj.Pts) {
      for(auto iht : tp.Hits) {
        if(inTraj[iht] <= 0) continue;
        itj = inTraj[iht];
        for(indx = 0; indx < tID.size(); ++indx) if(tID[indx] == itj) break;
        if(indx == tID.size()) {
          tID.push_back(itj);
          tCnt.push_back(1);
        }  else {
          ++tCnt[indx];
        }
      } // iht
    } // tp
    if(tCnt.size() == 0) return false;
/*
    mf::LogVerbatim myprt("TC");
    myprt<<"MaybeDeltaRay "<<tj.ID<<" size "<<tj.Pts.size()<<" tCnt";
    for(unsigned short ii = 0; ii < tID.size(); ++ii) myprt<<" "<<tID[ii]<<"_"<<tCnt[ii];
*/
    for(indx = 0; indx < tCnt.size(); ++indx) if(tCnt[indx] > tCut) {
      tj.AlgMod[kMaybeDeltaRay] = true;
      return true;
    }
    return false;
    
  } // MaybeDeltaRay

  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckWork()
  {
    // Check the quality of the work trajectory and possibly trim it.
   
    if(!fGoodWork) return;
    
    fTryWithNextPass = false;
    fCheckWorkModified = false;
    
    if(prt) mf::LogVerbatim("TC")<<"inside CheckWork";
    
    unsigned short ipt, ii, newSize;
    
    // check the fraction of wires that have hits
    float nTjWires = work.EndPt[1] - work.EndPt[0] + 1;
    float fracWiresWithGoodPts = NumPtsWithCharge(work) / nTjWires;
    if(fracWiresWithGoodPts < 0.6) {
      if(prt) mf::LogVerbatim("TC")<<"CW Poor trajectory. "<<fracWiresWithGoodPts<<" Trying to truncate and recover";
      // trim Pts off the end until this condition is satisfied
      while(fracWiresWithGoodPts > 0.6 && work.EndPt[1] > work.EndPt[1]) {
        work.Pts.pop_back();
        SetEndPoints(work);
        nTjWires = work.EndPt[1] - work.EndPt[0] + 1;
        fracWiresWithGoodPts = NumPtsWithCharge(work) / nTjWires;
        if(prt) mf::LogVerbatim("TC")<<" fracWiresWithGoodPts "<<fracWiresWithGoodPts;
      }
    } // fracWiresWithGoodPts < 0.6

    // ignore short trajectories
    if(work.EndPt[1] < 4) return;
    
    ipt = work.EndPt[1];
    if(!IsLargeAngle(work.Pts[ipt])) {
      // Not large angle checks
      // look for the signature of a kink near the end of the trajectory.
      // These are: Increasing chisq for the last few hits. Presence of
      // a removed hit near the end. A sudden decrease in the number of
      // TPs in the fit. A change in the average charge of hits. These may
      // not all be present in every situation.
      if(work.EndPt[1] > 8) {
        if(prt) PrintTrajectory(work, USHRT_MAX);
        unsigned short tpGap = USHRT_MAX;
        unsigned short nBigRat = 0;
        float chirat;
        for(ii = 1; ii < 3; ++ii) {
          ipt = work.EndPt[1] - 1 - ii;
          if(work.Pts[ipt].Chg == 0) {
            tpGap = ipt;
            break;
          }
          chirat = work.Pts[ipt+1].FitChi / work.Pts[ipt].FitChi;
          if(chirat > 1.5) ++nBigRat;
          if(prt) mf::LogVerbatim("TC")<<"CheckWork: chirat "<<ipt<<" chirat "<<chirat<<" "<<work.Pts[ipt-1].FitChi<<" "<<work.Pts[ipt].FitChi<<" "<<work.Pts[ipt+1].FitChi;
        } // ii
        if(prt) mf::LogVerbatim("TC")<<"CheckWork: nBigRat "<<nBigRat<<" tpGap "<<tpGap;
        if(tpGap != USHRT_MAX && nBigRat > 0) {
          if(tpGap != USHRT_MAX) {
            newSize = tpGap;
          } else {
            newSize = work.Pts.size() - 3;
          }
          if(prt) mf::LogVerbatim("TC")<<"  Setting work UseHits from "<<newSize<<" to "<<work.Pts.size()-1<<" false";
          for(ipt = newSize; ipt < work.Pts.size(); ++ipt) UnsetUsedHits(work.Pts[ipt]);
          SetEndPoints(work);
          fCheckWorkModified = true;
          work.Pts.resize(newSize);
          work.AlgMod[kCWKink] = true;
          return;
        }
      } // work.Pts.size() > 9
      
      // Compare the number of steps taken per TP near the beginning and
      // at the end
      short nStepBegin = work.Pts[2].Step - work.Pts[1].Step;
      short nStepEnd;
      unsigned short lastPt = work.Pts.size() - 1;
      newSize = work.Pts.size();
      for(ipt = lastPt; ipt > lastPt - 2; --ipt) {
        nStepEnd = work.Pts[ipt].Step - work.Pts[ipt - 1].Step;
        if(nStepEnd > 3 * nStepBegin) newSize = ipt;
      }
      if(prt) mf::LogVerbatim("TC")<<"CheckWork: check number of steps. newSize "<<newSize<<" work.Pts.size() "<<work.Pts.size();
      if(newSize < work.Pts.size()) {
        for(ipt = newSize; ipt < work.Pts.size(); ++ipt) UnsetUsedHits(work.Pts[ipt]);
        SetEndPoints(work);
        fCheckWorkModified = true;
        work.AlgMod[kCWStepChk] = true;
        work.Pts.resize(newSize);
        return;
      } // newSize < work.Pts.size()
    }
    
    // Check either large angle or not-large angles
    CheckHiDeltas();
    
    CheckHiMultUnusedHits();
    
  } // CheckWork
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiDeltas()
  {
    // Mask off hits at the beginning and end of the work trajectory if
    // Delta is too high
    
    // don't bother with LA trajectories
    if(IsLargeAngle(work.Pts[work.EndPt[1]])) return;
    
    // Check the beginning
    // find the first value of delta RMS that is not the default value
    unsigned short ipt, usePt = USHRT_MAX;
    for(ipt = work.EndPt[0] + fNPtsAve + 3; ipt < work.EndPt[1]; ++ipt) {
      if(work.Pts[ipt].Chg > 0 && work.Pts[ipt].DeltaRMS > 0 && work.Pts[ipt].DeltaRMS != 0.02) {
        usePt = ipt;
        break;
      }
    } // ipt
    float pull;
    if(usePt != USHRT_MAX) {
      // Scan from usePt back to the beginning. Stop if delta looks suspiciously large
      // and remove all points from the beginning to this point
      unsigned short ii;
      for(ii = 0; ii < work.Pts.size(); ++ii) {
        ipt = usePt - ii;
        pull = work.Pts[ipt].Delta / work.Pts[usePt].DeltaRMS;
        if(prt) mf::LogVerbatim("TC")<<"CHD begin "<<ipt<<" "<<fPlane<<":"<<PrintPos(work.Pts[ipt])<<" Delta "<<work.Pts[ipt].Delta<<" pull "<<pull<<" usePt "<<usePt;
        if(pull > 5) {
          // unset all TPs from here to the beginning
          for(unsigned short jpt = 0; jpt < ipt + 1; ++jpt) UnsetUsedHits(work.Pts[jpt]);
          work.AlgMod[kHiEndDelta] = true;
          break;
        } // pull > 5
        if(ipt == 0) break;
      } // ii
/*
      for(ipt = work.EndPt[0]; ipt < work.EndPt[0] + fNPtsAve; ++ipt) {
        pull = work.Pts[ipt].Delta / work.Pts[usePt].DeltaRMS;
        if(prt) mf::LogVerbatim("TC")<<"CHD begin "<<ipt<<" "<<fPlane<<":"<<PrintPos(work.Pts[ipt])<<" Delta "<<work.Pts[ipt].Delta<<" pull "<<pull<<" usePt "<<usePt;
        if(pull > 5) {
          UnsetUsedHits(work.Pts[ipt]);
          work.AlgMod[kHiEndDelta] = true;
        } // pull > 5
      } // ipt
*/
      if(work.AlgMod[kHiEndDelta] == true) SetEndPoints(work);
    } // usePt != USHRT_MAX
    
    // now check the other end using the same procedure
    usePt = USHRT_MAX;
    for(ipt = work.EndPt[1]; ipt > work.EndPt[0]; --ipt) {
      if(work.Pts[ipt].Chg > 0 && work.Pts[ipt].DeltaRMS > 0 && work.Pts[ipt].DeltaRMS != 0.02) {
        usePt = ipt;
        break;
      }
    } // ipt
    if(usePt == USHRT_MAX) return;
    for(ipt = usePt; ipt < work.EndPt[1] + 1; ++ipt) {
      pull = work.Pts[ipt].Delta / work.Pts[usePt].DeltaRMS;
      if(prt) mf::LogVerbatim("TC")<<"CHD end "<<ipt<<" "<<fPlane<<":"<<PrintPos(work.Pts[ipt])<<" Delta "<<work.Pts[ipt].Delta<<" pull "<<pull<<" usePt "<<usePt;
      if(pull > 5) {
        // unset all TPs from here to the beginning
        for(unsigned short jpt = ipt; jpt < work.EndPt[1] + 1; ++jpt) UnsetUsedHits(work.Pts[jpt]);
        work.AlgMod[kHiEndDelta] = true;
        break;
      } // pull > 5
    } // ipt
/*
    for(ipt = work.EndPt[1]; ipt > work.EndPt[0]; --ipt) {
      pull = work.Pts[ipt].Delta / work.Pts[usePt].DeltaRMS;
      if(prt) mf::LogVerbatim("TC")<<"CHD end "<<ipt<<" "<<fPlane<<":"<<PrintPos(work.Pts[ipt])<<" Delta "<<work.Pts[ipt].Delta<<" pull "<<pull;
      if(pull > 5) {
        UnsetUsedHits(work.Pts[ipt]);
        work.AlgMod[kHiEndDelta] = true;
      } // pull > 5
    } // ipt
*/
    if(work.AlgMod[kHiEndDelta] == true) SetEndPoints(work);
    
  } // CheckHiDeltas
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckHiMultUnusedHits()
  {
    // Check for many unused hits in high multiplicity TPs in work and try to use them
    
    // This code might do bad things to short trajectories
    if(work.Pts.size() < 10) return;
    
    // count the number of unused hits multiplicity > 1 hits and decide
    // if the unused hits should be used. This may trigger another
    // call to StepCrawl
    unsigned short ii, stopPt;
    // Use this to see if the high multiplicity Pts are mostly near
    // the end or di
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
    for(ii = 1; ii < work.Pts.size(); ++ii) {
      stopPt = work.EndPt[1] - ii;
      for(jj = 0; jj < work.Pts[stopPt].Hits.size(); ++jj) {
        iht = work.Pts[stopPt].Hits[jj];
        if(inTraj[iht] > 0) {
          doBreak = true;
          break;
        }
      } // jj
      if(doBreak) break;
      if(work.Pts[stopPt].Hits.size() == 1 && lastMult1Pt == USHRT_MAX) lastMult1Pt = stopPt;
      if(work.Pts[stopPt].Hits.size() > 1) {
        ++nHiMultPt;
        nHiMultPtHits += work.Pts[stopPt].Hits.size();
        nHiMultPtUsedHits += NumUsedHits(work.Pts[stopPt]);
      } // high multiplicity TP
      if(stopPt == 0) break;
    } // ii
    // Don't do this if there aren't a lot of high multiplicity TPs
    float fracHiMult = (float)nHiMultPt / (float)ii;
    float fracHitsUsed = 0;
    if(nHiMultPt > 0 && nHiMultPtHits > 0) fracHitsUsed = (float)nHiMultPtUsedHits / (float)nHiMultPtHits;
    if(prt) mf::LogVerbatim("TC")<<"CW First inTraj stopPt "<<stopPt<<" fracHiMult "<<fracHiMult<<" fracHitsUsed "<<fracHitsUsed;
    if(fracHiMult < 0.3) return;
    if(fracHitsUsed > 0.98) return;
    // Use this to limit the number of points fit for trajectories that
    // are close the LA tracking cut
    ii = work.EndPt[1];
    bool sortaLargeAngle = (std::abs(work.Pts[ii].Dir[0]) < fLargeAngle + 0.1);

    if(prt) mf::LogVerbatim("TC")<<" Pts size "<<work.Pts.size()<<" nHiMultPt "<<nHiMultPt<<" nHiMultPtHits "<<nHiMultPtHits<<" nHiMultPtUsedHits "<<nHiMultPtUsedHits<<" sortaLargeAngle "<<sortaLargeAngle;
    // Most of the hits are used so don't try too hard
    if(fracHitsUsed < 0.2) return;

    // Use next pass cuts if available
    if(sortaLargeAngle && work.Pass < fMinPtsFit.size()-1) ++work.Pass;

    // Make a copy of work and try to include the unused hits
    Trajectory newTj = work;
    unsigned short ipt;
    unsigned short iptStart = newTj.Pts.size() - 1;
    // remove any TPs at the end that have no hits
    for(ipt = iptStart; ipt > newTj.EndPt[1]; --ipt) if(newTj.Pts[ipt].Chg == 0) newTj.Pts.pop_back();
    // unset all of the used hits
    for(auto& tp : newTj.Pts) UnsetUsedHits(tp);
    if(stopPt < 2) {
      for(ipt = 0; ipt < 3; ++ipt) SetAllHitsUsed(newTj.Pts[ipt]);
      SetEndPoints(newTj);
      UpdateTraj(newTj);
     } // stopPt < 2
    unsigned short killPts;
    for(ipt = newTj.EndPt[1] + 1; ipt < newTj.Pts.size(); ++ipt) {
      SetAllHitsUsed(newTj.Pts[ipt]);
      SetEndPoints(newTj);
      // This will be incremented by one in UpdateTraj
      if(sortaLargeAngle) newTj.Pts[ipt].NTPsFit = 2;
      UpdateTraj(newTj);
      if(!fUpdateTrajOK) {
        if(prt) mf::LogVerbatim("TC")<<"UpdateTraj failed on point "<<ipt<<"\n";
        return;
      }
      GottaKink(newTj, killPts);
      if(killPts > 0) {
        MaskTrajEndPoints(newTj, killPts);
        break;
      }
      if(prt) PrintTrajectory(newTj, ipt);
    } // ipt
    // if we made it here it must be OK
    SetEndPoints(newTj);
    if(newTj.EndPt[0] == newTj.EndPt[1]) return;
    // Try to extend it, unless there was a kink
    if(work.AlgMod[kGottaKink]) return;
    // trim the end points although this shouldn't happen
    if(newTj.EndPt[1] != newTj.Pts.size() - 1) newTj.Pts.resize(newTj.EndPt[1] + 1);
    work = newTj;
    work.AlgMod[kManyHitsAdded] = true;
    if(prt) mf::LogVerbatim("TC")<<"TRP CheckHiMultUnusedHits successfull. Calling StepCrawl to extend it";
    StepCrawl();
    
  } // CheckHiMultUnusedHits

  ////////////////////////////////////////////////
  bool TrajClusterAlg::MaskedWorkHitsOK()
  {
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration or whether to stop and possibly try with the next pass settings
    unsigned short lastPt = work.Pts.size() - 1;
    if(NumUsedHits(work.Pts[lastPt]) > 0) return true;
    
    // count the number of points w/o used hits and the number of close hits
    unsigned short nClose = 0, nMasked = 0;
    unsigned short ii, indx;
    for(ii = 0; ii < work.Pts.size(); ++ii) {
      indx = work.Pts.size() - 1 - ii;
      if(NumUsedHits(work.Pts[indx]) > 0) break;
      ++nMasked;
      nClose += work.Pts[indx].Hits.size();
    } // ii
    
    if(nMasked == 0) return true;
    
    if(prt) mf::LogVerbatim("TC")<<"MaskedWorkHitsOK:  nMasked "<<nMasked<<" nClose "<<nClose<<" fMaxWireSkipWithSignal "<<fMaxWireSkipWithSignal;

    // See if we have masked off several TPs and there are a number of close hits that haven't
    // been added. This indicates that maybe we should be using the next pass cuts. The approach
    // here is to add the (previously excluded) hits in the fit and re-fit with next pass settings
    if(work.Pass < (fMinPtsFit.size()-1) && work.Pts.size() > 10 && nMasked > 1 && nClose < 3 * nMasked) {
      // set all of the close hits used for the last TPs. This is a bit scary since we
      // aren't sure that all of these hits are the right ones
      for(ii = 0; ii < work.Pts.size(); ++ii) {
        indx = work.Pts.size() - 1 - ii;
        if(NumUsedHits(work.Pts[indx]) > 0) break;
        unsigned int iht;
        for(unsigned short jj = 0; jj < work.Pts[indx].Hits.size(); ++jj) {
          iht = work.Pts[indx].Hits[jj];
          if(inTraj[iht] > 0) continue;
          work.Pts[indx].UseHit[jj] = true;
          inTraj[iht] = work.ID;
        } // jj
      } // ii
      SetEndPoints(work);
      PrepareWorkForNextPass();
      if(prt) mf::LogVerbatim("TC")<<"MaskedWorkHitsOK: Try next pass with too many missed close hits "<<fTryWithNextPass;
      work.AlgMod[kMaskedWorkHits] = true;
      return false;
    }

    // Be a bit more lenient with short trajectories on the first pass if
    // the FitChi is not terribly bad and there is ony one hit associated with the last TP
    if(work.Pass < (fMinPtsFit.size()-1) && work.Pts.size() > 5 && work.Pts.size() < 15 && nMasked < 4
       && work.Pts[lastPt].FitChi < 2 * fMaxChi && work.Pts[lastPt].Hits.size() == 1) {
      // set this hit used if it is available
      unsigned int iht = work.Pts[lastPt].Hits[0];
      if(inTraj[iht] <= 0) {
        work.Pts[lastPt].UseHit[0] = true;
        inTraj[iht] = work.ID;
      }
      SetEndPoints(work);
      PrepareWorkForNextPass();
      if(prt) mf::LogVerbatim("TC")<<"MaskedWorkHitsOK: Try next pass with kinda short, kinda bad trajectory. "<<fTryWithNextPass;
      return false;
    }

    // OK if we haven't exceeded the user cut
    if(nMasked < fMaxWireSkipWithSignal) return true;
    
    // not a lot of close hits try with the next pass settings
    if(nClose < 1.5 * nMasked) {
      // trim the trajectory
      unsigned short newSize = work.Pts.size() - nMasked;
      if(prt) mf::LogVerbatim("TC")<<"MaskedWorkHitsOK:  Trimming work to size "<<newSize;
      work.Pts.resize(newSize);
      PrepareWorkForNextPass();
    }
    return false;
    
  } // MaskedWorkHitsOK
/*
  ////////////////////////////////////////////////
  unsigned short TrajClusterAlg::NumGoodWorkTPs()
  {
    unsigned short cnt = 0;
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      if(work.Pts[ipt].Chg == 0) continue;
      ++cnt;
    }
    return cnt;
  } // NumGoodWorkTPs
*/
  ////////////////////////////////////////////////
  void TrajClusterAlg::PrepareWorkForNextPass()
  {
    // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
    // fitted points to get FitCHi < 2
    
    fTryWithNextPass = false;

    // See if there is another pass available
    if(work.Pass > fMinPtsFit.size()-2) return;
    ++work.Pass;
    
    unsigned short lastPt = work.Pts.size() - 1;
    // Return if the last fit chisq is OK
    if(work.Pts[lastPt].FitChi < 1.5) {
      fTryWithNextPass = true;
      return;
    }
    TrajPoint& lastTP = work.Pts[lastPt];
    unsigned short newNTPSFit = lastTP.NTPsFit;
    // only give it a few tries before giving up
    unsigned short nit = 0;

    while(lastTP.FitChi > 1.5 && lastTP.NTPsFit > 2) {
      if(lastTP.NTPsFit > 3) newNTPSFit -= 2;
      else if(lastTP.NTPsFit == 3) newNTPSFit = 2;
      lastTP.NTPsFit = newNTPSFit;
      FitWork();
      if(prt) mf::LogVerbatim("TC")<<"PrepareWorkForNextPass: FitChi is > 1.5 "<<lastTP.FitChi<<" Reduced NTPsFit to "<<lastTP.NTPsFit<<" work.Pass "<<work.Pass;
      if(lastTP.NTPsFit <= fMinPtsFit[work.Pass]) break;
      ++nit;
      if(nit == 3) break;
    }
    // decide if the next pass should indeed be attempted
    if(lastTP.FitChi > 2) return;
    fTryWithNextPass = true;
    fAveChg = lastTP.AveChg;
    
  } // PrepareWorkForNextPass
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FindHit(std::string someText, unsigned int iht)
  {
    // finds a hit in work
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      TrajPoint& tp = work.Pts[ipt];
      if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) {
        mf::LogVerbatim("TC")<<someText<<" Found hit "<<work.CTP<<":"<<PrintHit(iht)<<" in work. inTraj = "<<inTraj[iht]<<" work trajectory below ";
        PrintTrajectory(work, USHRT_MAX);
        return;
      }
    } // ipt
    // look in allTraj
    for(unsigned short itj = 0; itj < allTraj.size(); ++itj) {
      for(unsigned short ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
        TrajPoint& tp = allTraj[itj].Pts[ipt];
        if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) {
          mf::LogVerbatim("TC")<<someText<<" Found hit "<<allTraj[itj].CTP<<" "<<PrintHit(iht)<<" inTraj "<<inTraj[iht]<<" allTraj ID "<<allTraj[itj].ID<<" trajectory below ";
          PrintTrajectory(allTraj[itj], USHRT_MAX);
          return;
        }
      } // ipt
    } // itj
  } // FindHit

  ////////////////////////////////////////////////
  bool TrajClusterAlg::CheckAllHitsFlag()
  {
    // returns true if all hits have flag != -3
    unsigned int firsthit = (unsigned int)WireHitRange[fFirstWire].first;
    unsigned int lasthit = (unsigned int)WireHitRange[fLastWire-1].second;
    bool itsBad = false;
    for(unsigned int iht = firsthit; iht < lasthit; ++iht) {
      if(inTraj[iht] < 0) {
//        mf::LogVerbatim("TC")<<"CheckAllHitsFlag bad "<<PrintHit(iht);
        inTraj[iht] = 0;
        itsBad = true;
//        break;
      }
    } // iht
    
    if(!itsBad) return true;
    
    // print out the information
    mf::LogVerbatim myprt("TC");
    myprt<<"Not released hits ";
    for(unsigned int iht = firsthit; iht < lasthit; ++iht) {
      if(inTraj[iht] < 0) {
        myprt<<" "<<iht<<"_"<<fPlane<<":"<<PrintHit(iht)<<" inTraj set to "<<inTraj[iht];
        // look for it in a trajectory
        for(unsigned short itj = 0; itj < allTraj.size(); ++itj) {
          for(unsigned short ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
            TrajPoint& tp = allTraj[itj].Pts[ipt];
            if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) {
              myprt<<" in traj "<<itj;
            } // found it
          } // ip
        } // itj
        inTraj[iht] = 0;
      }
    } // iht
    return false;
    
  } // CheckAllHitsFlag

   
  ////////////////////////////////////////////////
  void TrajClusterAlg::GottaKink(Trajectory& tj, unsigned short& killPts)
  {
    // This routine requires that EndPt is defined
    killPts = 0;

    unsigned short lastPt = tj.EndPt[1];
    if(lastPt < 2) return;
    if(tj.Pts[lastPt].Chg == 0) return;
    
    float dang;
    
    // A simple check when there are few points being fit
    if(tj.Pts[lastPt].NTPsFit < 6) {
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
      dang = std::abs(tj.Pts[lastPt].Ang - tj.Pts[prevPtWithHits].Ang);
      if(prt) mf::LogVerbatim("TC")<<"GottaKink Simple check lastPt "<<lastPt<<" prevPtWithHits "<<prevPtWithHits<<" dang "<<dang<<" cut "<<fKinkAngCut;
      if(dang > fKinkAngCut) {
        killPts = 1;
        tj.AlgMod[kGottaKink] = true;
        tj.Pts[prevPtWithHits].KinkAng = dang;
      }
      return;
    } // tj.NTPSFit < 4

    
    if(tj.EndPt[1] < 10) return;
    
    unsigned short kinkPt = USHRT_MAX;
    unsigned short ii, ipt, cnt;
    
    // Find the kinkPt which is the third Pt from the end that has charge
    cnt = 0;
    for(ii = 1; ii < lastPt; ++ii) {
      ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg == 0) continue;
      ++cnt;
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
 
    dang = std::abs(tj.Pts[kinkPt].Ang - tpFit.Ang);
    if(tpFit.FitChi > 900) return;
    
    tj.Pts[kinkPt].KinkAng = dang;

    float kinkSig = dang / tpFit.AngErr;
    if(prt) mf::LogVerbatim("TC")<<"GottaKink kinkPt "<<kinkPt<<" Pos "<<PrintPos(tj.Pts[kinkPt])<<" dang "<<dang<<" cut "<<fKinkAngCut<<" kinkSig "<<kinkSig<<" tpFit chi "<<tpFit.FitChi;
    
    if(dang > fKinkAngCut) {
      // Kink larger than hard user cut
      killPts = 3;
    } else if(dang > 0.5 * fKinkAngCut && kinkSig > 10 && tpFit.FitChi < 2) {
      // Found a very significant smaller angle kink
      killPts = 3;
    }
    
    if(killPts > 0) tj.AlgMod[kGottaKink] = true;
    
  } // GottaKink

 
  //////////////////////////////////////////
  unsigned short TrajClusterAlg::SetMissedStepCut(TrajPoint const& tp)
  {
    // maximum of 100 WSE's
    unsigned short itmp = 100;
    if(std::abs(tp.Dir[0]) > 0.01) itmp = 1 + 4 * std::abs(1/tp.Dir[0]);
    if(itmp < 2) itmp = 2;
    return itmp;
  }

  //////////////////////////////////////////
  void TrajClusterAlg::UpdateTraj(Trajectory& tj)
  {
    // Updates the last added trajectory point fit, average hit rms, etc.

    fUpdateTrajOK = false;
    fMaskedLastTP = false;
    
    if(tj.EndPt[1] < 1) return;
    unsigned int lastPt = tj.EndPt[1];
    TrajPoint& lastTP = tj.Pts[lastPt];

    // find the previous TP that was has hits (and was therefore in the fit)
    unsigned short ii, prevPtWithHits = USHRT_MAX;
    unsigned short ipt;
    for(ii = 1; ii < tj.Pts.size(); ++ii) {
//      ipt = tj.Pts.size() - 1 - ii;
      ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg > 0) {
        prevPtWithHits = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(prevPtWithHits == USHRT_MAX) return;

    if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: lastPt "<<lastPt<<" previous point with hits "<<prevPtWithHits<<" tj.Pts size "<<tj.Pts.size()<<" LargeAngle? "<<IsLargeAngle(lastTP);
    
//    if(NumUsedHits(lastTP) == 0) return;
    
    // Set the lastPT delta before doing the fit
    lastTP.Delta = PointTrajDOCA(lastTP.HitPos[0], lastTP.HitPos[1], lastTP);
    
    UpdateAveChg(tj, lastPt, -1);

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
    } else if(lastPt == 2) {
       // Third trajectory point. Keep it simple
      lastTP.NTPsFit = 3;
      FitTraj(tj);
      fUpdateTrajOK = true;
      return;
    } else {
      // Fit with > 2 TPs
      lastTP.NTPsFit += 1;
      FitTraj(tj);
      if(lastTP.FitChi > 2 && tj.Pts.size() > 10) {
        // Have a longish trajectory and chisq was a bit large.
        // Was this a sudden occurrence and the fraction of TPs are included
        // in the fit? If so, we should mask off this
        // TP and keep going. If these conditions aren't met, we
        // should reduce the number of fitted points
        float chirat = 0;
        if(prevPtWithHits != USHRT_MAX) chirat = lastTP.FitChi / tj.Pts[prevPtWithHits].FitChi;
        fMaskedLastTP = (chirat > 1.5 && lastTP.NTPsFit > 0.5 * tj.Pts.size());
        if(prt) mf::LogVerbatim("TC")<<" First fit chisq too large "<<lastTP.FitChi<<" prevPtWithHits chisq "<<tj.Pts[prevPtWithHits].FitChi<<" chirat "<<chirat<<" fMaskedLastTP "<<fMaskedLastTP;
        // we should also mask off the last TP if there aren't enough hits
        // to satisfy the minPtsFit constraint
        if(!fMaskedLastTP && NumPtsWithCharge(tj) < fMinPtsFit[tj.Pass]) fMaskedLastTP = true;
        // A large chisq jump can occur if we just jumped a large block of dead wires. In
        // this case we don't want to mask off the last TP but reduce the number of fitted points
        if(fMaskedLastTP) {
          float nMissedWires = std::abs(lastTP.HitPos[0] - tj.Pts[prevPtWithHits].HitPos[0]);
          if(nMissedWires > 5) {
            nMissedWires -= DeadWireCount(lastTP.HitPos[0], tj.Pts[prevPtWithHits].HitPos[0]);
            fMaskedLastTP = (nMissedWires > 5);
          } // nMissedWires > 5
        } // fMaskedLastTP
      } // lastTP.FitChi > 2 ...
      
      if(fMaskedLastTP && tj.Pass == fMinPtsFit.size() - 1 && lastTP.NTPsFit > 2 * fMinPts[tj.Pass]) {
        // Turn off the TP mask on the last pass. Try to reduce the number of TPs fit instead
        fMaskedLastTP = false;
      }
      
      if(prt) mf::LogVerbatim("TC")<<"UpdateTraj: First fit "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" FitChi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit<<" fMaskedLastTP "<<fMaskedLastTP;
      if(fMaskedLastTP) {
        UnsetUsedHits(lastTP);
        DefineHitPos(lastTP);
        SetEndPoints(tj);
        lastPt = tj.EndPt[1];
        lastTP.NTPsFit -= 1;
        FitTraj(tj);
        fUpdateTrajOK = true;
        return;
      }  else {
        // a more gradual increase in chisq. Reduce the number of points
        unsigned short newNTPSFit = lastTP.NTPsFit;
        // reduce the number of points fit to keep Chisq/DOF < 2 adhering to the pass constraint
        while(lastTP.FitChi > 2 && lastTP.NTPsFit > 2) {
          if(lastTP.NTPsFit > 15) {
            newNTPSFit = 0.7 * newNTPSFit;
          } else if(lastTP.NTPsFit > 4) {
            newNTPSFit -= 2;
          } else {
            newNTPSFit -= 1;
          }
          if(lastTP.NTPsFit < 3) newNTPSFit = 2;
          lastTP.NTPsFit = newNTPSFit;
          if(prt) mf::LogVerbatim("TC")<<"  Bad FitChi "<<lastTP.FitChi<<" Reduced NTPsFit to "<<lastTP.NTPsFit<<" tj "<<tj.Pass;
          FitTraj(tj);
          if(lastTP.NTPsFit <= fMinPtsFit[tj.Pass]) break;
        } // lastTP.FitChi > 2 && lastTP.NTPsFit > 2
      }
      // last ditch attempt. Drop the last hit
      if(tj.Pts.size() > 2 && lastTP.FitChi > 2) {
        if(prt) mf::LogVerbatim("TC")<<"  Last try. Drop last TP "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;
        UnsetUsedHits(lastTP);
        DefineHitPos(lastTP);
        SetEndPoints(tj);
        lastPt = tj.EndPt[1];
        FitTraj(tj);
        fUpdateTrajOK = true;
        fMaskedLastTP = true;
      }
      if(prt) mf::LogVerbatim("TC")<<"  Fit done. Chi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;
    } // tj.Pts.size
    
    if(tj.EndPt[0] == tj.EndPt[1]) {
      mf::LogWarning("TC")<<"UpdateTraj: bad endpoints "<<tj.Pts.size();
      PrintTrajectory(tj, USHRT_MAX);
      fUpdateTrajOK = false;
      return;
    }
    
    // Don't let the angle error get too small too soon. Crawling would stop if the first
    // few hits on a low momentum wandering track happen to have a very good fit to a straight line.
    // We will do this by averaging the default starting value of AngErr of the first TP with the current
    // value from FitTraj.
    if(lastPt < 14) {
      float defFrac = 1 / (float)(tj.EndPt[1]);
      lastTP.AngErr = defFrac * tj.Pts[0].AngErr + (1 - defFrac) * lastTP.AngErr;
    }
    
    // put the fit result into all previous points that
    // were in the fit if it was really good
    if(lastTP.FitChi < 1.5) {
      unsigned short cnt = 0;
      TrajPoint tp;
      for(ii = 1; ii < tj.Pts.size(); ++ii) {
        ipt = lastPt - ii;
        if(tj.Pts[ipt].Chg > 0) {
          tj.Pts[ipt].Dir = lastTP.Dir;
          tj.Pts[ipt].Ang = lastTP.Ang;
          tj.Pts[ipt].AngErr = lastTP.AngErr;
          tp = lastTP;
          MoveTPToWire(tp, tj.Pts[ipt].Pos[0]);
          tj.Pts[ipt].Pos[1] = tp.Pos[1];
          tj.Pts[ipt].Delta = PointTrajDOCA(tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tj.Pts[ipt]);
          tj.Pts[ipt].NTPsFit = lastTP.NTPsFit;
          tj.Pts[ipt].FitChi = lastTP.FitChi;
          ++cnt;
        }
        if(cnt == lastTP.NTPsFit) break;
        if(ipt == 0) break;
      } // ii
    }

    UpdateDeltaRMS(tj);

    fUpdateTrajOK = true;
    return;

  } // UpdateWork

  //////////////////////////////////////////
  void TrajClusterAlg::MoveTPToWire(TrajPoint& tp, float wire)
  {
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    if(tp.Dir[0] == 0) return;
    float dw = wire - tp.Pos[0];
    if(dw == 0) return;
    tp.Pos[0] = wire;
    tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
  } // MoveTPToWire

  //////////////////////////////////////////
  void TrajClusterAlg::UpdateDeltaRMS(Trajectory& tj)
  {
    // Estimate the Delta RMS of the TPs on the end of tj.
    
    unsigned int lastPt = tj.EndPt[1];
    TrajPoint& lastTP = tj.Pts[lastPt];
    // put in reasonable guess in case something doesn't work out
    lastTP.DeltaRMS = 0.02;
    
    if(lastTP.Chg == 0) return;
    if(lastPt < 6) return;
    
    float sum = 0;
    unsigned short ii, ipt, cnt = 0;
     for(ii = tj.EndPt[0]; ii < lastPt; ++ii) {
      ipt = lastPt - ii;
      if(ipt > tj.Pts.size() - 1) {
        std::cout<<"UpdateDeltaRMS oops "<<ii<<" ipt "<<lastPt<<" size "<<tj.Pts.size();
        exit(1);
      }
       // Delta for the first two points is meaningless
       if(ipt < 2) break;
      if(tj.Pts[ipt].Chg == 0) continue;
      sum += tj.Pts[ipt].Delta;
      ++cnt;
      if(cnt == fNPtsAve) break;
    }
    if(cnt < 2) return;
    // RMS of Gaussian distribution is ~1.2 x the average
    // of a one-sided Gaussian distribution (since Delta is > 0)
    lastTP.DeltaRMS = 1.2 * sum / (float)cnt;
    if(lastTP.DeltaRMS < 0.02) lastTP.DeltaRMS = 0.02;

  } // UpdateDeltaRMS
  
  //////////////////////////////////////////
  void TrajClusterAlg::FitWork()
  {
    // Jacket around FitTraj to fit the leading edge of the supplied trajectory
    unsigned short originPt = work.Pts.size() - 1;
    unsigned short npts = work.Pts[originPt].NTPsFit;
    TrajPoint tpFit;
    unsigned short fitDir = -1;
    FitTraj(work, originPt, npts, fitDir, tpFit);
    work.Pts[originPt] = tpFit;
    
  } // FitWork
  
  //////////////////////////////////////////
  void TrajClusterAlg::FitTraj(Trajectory& tj)
  {
    // Jacket around FitTraj to fit the leading edge of the supplied trajectory
    unsigned short originPt = tj.Pts.size() - 1;
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
    
    if(originPt > tj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"FitTraj: Requesting fit of invalid TP "<<originPt;
      return;
    }
    
    // copy the origin TP into the fit TP
    tpFit = tj.Pts[originPt];
    // Assume that the fit will fail
    tpFit.FitChi = 999;
    if(fitDir < -1 || fitDir > 1) return;
    
    std::vector<double> x, y, w, q;
    std::array<float, 2> dir, origin = tj.Pts[originPt].HitPos;
    // Use TP position if there aren't any hits on it
    if(tj.Pts[originPt].Chg == 0) origin = tj.Pts[originPt].Pos;
    double xx, yy, xr, yr;

    // Rotate the traj hit position into the coordinate system defined by the
    // originPt traj point, where x = along the trajectory, y = transverse
    double rotAngle = tj.Pts[originPt].Ang;
    double cs = cos(-rotAngle);
    double sn = sin(-rotAngle);
    
    unsigned short ipt, cnt;
    double aveChg = 0;
    // enter the originPT hit info if it exists
    if(tj.Pts[originPt].Chg > 0) {
      xx = tj.Pts[originPt].HitPos[0] - origin[0];
      yy = tj.Pts[originPt].HitPos[1] - origin[1];
//      if(prt) std::cout<<" originPT "<<originPt<<" xx "<<xx<<" "<<yy<<" chg "<<tj.Pts[originPt].Chg<<"\n";
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      x.push_back(xr);
      y.push_back(yr);
      w.push_back(tj.Pts[originPt].HitPosErr2);
      q.push_back(tj.Pts[originPt].Chg);
      aveChg += tj.Pts[originPt].Chg;
    }
    
    // correct npts to account for the origin point
    if(fitDir != 0) --npts;
    
    // step in the + direction first
    if(fitDir != -1) {
      cnt = 0;
      for(ipt = originPt + 1; ipt < tj.Pts.size(); ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
//        if(prt) std::cout<<"ipt "<<ipt<<" xx "<<xx<<" yy "<<yy<<" chg "<<tj.Pts[ipt].Chg<<"\n";
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        w.push_back(tj.Pts[ipt].HitPosErr2);
        q.push_back(tj.Pts[ipt].Chg);
        aveChg += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == npts) break;
      } // ipt
    } // fitDir != -1
    
    // step in the - direction next
    if(fitDir != 1 && originPt > 0) {
      cnt = 0;
      for(ipt = originPt - 1; ipt > 0; --ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
//        if(prt) std::cout<<"ipt "<<ipt<<" xx "<<xx<<" "<<yy<<" chg "<<tj.Pts[ipt].Chg<<"\n";
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        w.push_back(tj.Pts[ipt].HitPosErr2);
        q.push_back(tj.Pts[ipt].Chg);
        aveChg += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == npts) break;
        if(ipt == 0) break;
      } // ipt
    } // fitDir != -1
    
    // Not enough points to define a line?
    if(x.size() < 2) return;
    
//    if(prt) std::cout<<"FitTraj: npts "<<npts<<" origin "<<origin[0]<<" "<<origin[1]<<" ticks "<<origin[1]/fScaleF<<" rotAngle "<<rotAngle<<"\n";
    
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;

    // weight by the charge ratio and accumulate sums
    aveChg /= (double)cnt;
    double chgrat, wght;
    for(ipt = 0; ipt < x.size(); ++ipt) {
      chgrat = std::abs(q[ipt] - aveChg) / aveChg;
      if(chgrat < 0.3) chgrat = 0.3;
      w[ipt] *= chgrat;
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
    
    // The chisq will be set below if there are enough points
    tpFit.FitChi = 0.01;
    double newang = atan(B);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAngle);
    sn = sin(rotAngle);
    tpFit.Dir[0] = cs * dir[0] - sn * dir[1];
    tpFit.Dir[1] = sn * dir[0] + cs * dir[1];
    // Reverse the direction?
    bool flipDir = false;
    if(std::abs(tpFit.Dir[0]) > std::abs(tpFit.Dir[1])) {
      flipDir = std::signbit(tpFit.Dir[0]) != std::signbit(tj.Pts[originPt].Dir[0]);
    } else {
      flipDir = std::signbit(tpFit.Dir[1]) != std::signbit(tj.Pts[originPt].Dir[1]);
    }
    if(flipDir) {
      tpFit.Dir[0] = -tpFit.Dir[0];
      tpFit.Dir[1] = -tpFit.Dir[1];
    }
    tpFit.Ang = atan2(tpFit.Dir[1], tpFit.Dir[0]);
    // rotate (0, intcpt) into W,T
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
//      InterceptError = 0.;
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
  void TrajClusterAlg::UpdateAveChg(Trajectory& tj, unsigned short updatePt, short dir)
  {
    
    if(updatePt < tj.EndPt[0]) return;
    if(updatePt > tj.EndPt[1]) return;
    
    // Find the average charge using fNPtsAve values of TP Chg. Result put in fAveChg and Pts[updatePt].Chg
    unsigned short ii, ipt, cnt = 0;
    
    float sum = 0;
    if(dir > 0) {
      for(ipt = updatePt; ipt < tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg > 0) {
          sum += tj.Pts[ipt].Chg;
          ++cnt;
        }
        if(cnt == fNPtsAve) break;
      }
    } else {
      // dir < 0
      for(ii = 0; ii < tj.Pts.size(); ++ii) {
        ipt = updatePt - ii;
        if(tj.Pts[ipt].Chg > 0) {
          sum += tj.Pts[ipt].Chg;
          ++cnt;
        }
        if(cnt == fNPtsAve) break;
        if(ipt == tj.EndPt[0]) break;
      }
    } // dir < 0
    
    if(cnt < fNPtsAve) return;
    fAveChg = sum / (float)cnt;
    if(tj.Pts[0].AveChg == 0) {
      // update all TPs to the current one
      for(unsigned short jpt = 0; jpt < updatePt; ++jpt) {
        tj.Pts[jpt].AveChg = fAveChg;
        if(tj.Pts[jpt].Chg > 0) tj.Pts[jpt].ChgRat = (tj.Pts[jpt].Chg - fAveChg) / fAveChg;
      } // jpt
    } // fAveChg == 0
     // update the point at the update point
    tj.Pts[updatePt].AveChg = fAveChg;
    tj.Pts[updatePt].ChgRat = (tj.Pts[updatePt].Chg - fAveChg) / fAveChg;
    
  } // UpdateAveChg

  ////////////////////////////////////////////////
  void TrajClusterAlg::StartWork(unsigned int fromHit, unsigned int toHit)
  {
    float fromWire = fHits[fromHit]->WireID().Wire;
    float fromTick = fHits[fromHit]->PeakTime();
    float toWire = fHits[toHit]->WireID().Wire;
    float toTick = fHits[toHit]->PeakTime();
    CTP_t tCTP = EncodeCTP(fHits[fromHit]->WireID());
    StartWork(fromWire, fromTick, toWire, toTick, tCTP);
  } // StartWork

  ////////////////////////////////////////////////
  void TrajClusterAlg::StartWork(float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP)
  {
    // Start a simple (seed) trajectory going from a hit to a position (toWire, toTick).
    // The traj vector is cleared if an error occurs
    
    // construct a default trajectory
    Trajectory tj;
    // and use it to blow out work
    work = tj;
    work.ID = -1;
    work.Pass = fPass;
    work.StepDir = fStepDir;
    work.CTP = tCTP;
    
    // create a trajectory point
    TrajPoint tp;
    MakeBareTrajPoint(fromWire, fromTick, toWire, toTick, tCTP, tp);

    tp.AngErr = 0.1;
    if(prt) mf::LogVerbatim("TC")<<"StartWork "<<(int)fromWire<<":"<<(int)fromTick<<" -> "<<(int)toWire<<":"<<(int)toTick<<" dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" ang "<<tp.Ang<<" angErr "<<tp.AngErr;
    work.Pts.push_back(tp);
    
    fAveChg = 0;
    
  } // StartWork
  
  //////////////////////////////////////////
  void TrajClusterAlg::ReverseTraj(Trajectory& tj)
  {
    // reverse the trajectory
    if(tj.Pts.size() == 0) return;
    // reverse the crawling direction flag
    tj.StepDir = -tj.StepDir;
    // Vertices
    short tmp = tj.Vtx[0];
    tj.Vtx[0] = tj.Vtx[1];
    tj.Vtx[0] = tmp;
    // trajectory points
    std::reverse(tj.Pts.begin(), tj.Pts.end());
    // reverse the direction vector on all points
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Dir[0] != 0) tj.Pts[ipt].Dir[0] = -tj.Pts[ipt].Dir[0];
      if(tj.Pts[ipt].Dir[1] != 0) tj.Pts[ipt].Dir[1] = -tj.Pts[ipt].Dir[1];
      tj.Pts[ipt].Ang = atan2(tj.Pts[ipt].Dir[1], tj.Pts[ipt].Dir[0]);
      // and the order of hits if more than 1
      if(tj.Pts[ipt].Hits.size() > 1) {
        std::reverse(tj.Pts[ipt].Hits.begin(), tj.Pts[ipt].Hits.end());
        std::reverse(tj.Pts[ipt].UseHit.begin(), tj.Pts[ipt].UseHit.end());
      }
    } // ipt
    SetEndPoints(tj);
  }
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckInTraj(std::string someText)
  {
    // Check allTraj -> inTraj associations
    unsigned short tID;
    unsigned int iht;
    unsigned short itj = 0;
    std::vector<unsigned int> tHits;
    std::vector<unsigned int> atHits;
    for(auto& tj : allTraj) {
      // ignore abandoned trajectories
      if(tj.AlgMod[kKilled]) continue;
      tID = tj.ID;
      PutTrajHitsInVector(tj, true,  tHits);
      if(tHits.size() < 2) {
        mf::LogVerbatim("TC")<<"CheckInTraj: Insufficient hits in tHits";
        PrintTrajectory(tj, USHRT_MAX);
        continue;
      }
      std::sort(tHits.begin(), tHits.end());
      atHits.clear();
      for(iht = 0; iht < inTraj.size(); ++iht) {
        if(inTraj[iht] == tID) atHits.push_back(iht);
      } // iht
      if(atHits.size() < 2) {
        mf::LogVerbatim("TC")<<"CheckInTraj: Insufficient hits in atHits";
        PrintTrajectory(tj, USHRT_MAX);
        continue;
      }
      if(!std::equal(tHits.begin(), tHits.end(), atHits.begin())) {
        mf::LogVerbatim myprt("TC");
        myprt<<"CheckIntTraj: allTraj - inTraj mis-match for tj ID "<<tID<<" atHits size "<<atHits.size()<<" tHits size "<<tHits.size()<<"\n";
        for(iht = 0; iht < atHits.size(); ++iht) {
          myprt<<"iht "<<iht<<" "<<PrintHit(atHits[iht]);
          if(iht < tHits.size()) myprt<<" "<<PrintHit(tHits[iht]);
          if(atHits[iht] != tHits[iht]) myprt<<" <<< ";
          myprt<<"\n";
          fQuitAlg = true;
        } // iht
        if(tHits.size() > atHits.size()) {
          for(iht = atHits.size(); iht < atHits.size(); ++iht) {
            myprt<<"atHits "<<iht<<" "<<PrintHit(atHits[iht])<<"\n";
          } // iht
        } // tHit.size > atHits.size()
      }
      ++itj;
      if(fQuitAlg) return;
    } // tj
    
  } // CheckInTraj

  ////////////////////////////////////////////////
  void TrajClusterAlg::StoreWork()
  {

    // put trajectories in order of US -> DS
    if(work.StepDir < 0) ReverseTraj(work);
    // This shouldn't be necessary but do it anyway
    SetEndPoints(work);
    if(work.EndPt[1] <= work.EndPt[0]) {
      mf::LogWarning("TC")<<"StoreWork: No meaningfull points on work trajectory";
      std::cout<<"StoreWork: No meaningfull points on work trajectory\n";
      return;
    }
    short trID = allTraj.size() + 1;
    unsigned short ii, iht;
    for(unsigned short ipt = work.EndPt[0]; ipt < work.EndPt[1] + 1; ++ipt) {
      for(ii = 0; ii < work.Pts[ipt].Hits.size(); ++ii) {
        if(work.Pts[ipt].UseHit[ii]) {
          iht = work.Pts[ipt].Hits[ii];
          if(inTraj[iht] > 0) {
            mf::LogWarning("TC")<<"StoreWork: Trying to store hit "<<PrintHit(iht)<<" in new allTraj "<<trID<<" but it is used in traj ID = "<<inTraj[iht]<<" print work and quit";
            PrintTrajectory(work, USHRT_MAX);
            ReleaseWorkHits();
            fQuitAlg = true;
            return;
          } // error
          inTraj[iht] = trID;
        }
      } // ii
    } // ipt
    work.ID = trID;
    allTraj.push_back(work);
    if(prt) mf::LogVerbatim("TC")<<"StoreWork trID "<<trID<<" CTP "<<work.CTP<<" EndPts "<<work.EndPt[0]<<" "<<work.EndPt[1];
    
  } // StoreWork

  ////////////////////////////////////////////////
  void TrajClusterAlg::SetEndPoints(Trajectory& tj)
  {
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
    
    tj.EndPt[0] = 0; tj.EndPt[1] = 0;
    if(tj.Pts.size() == 0) return;
    
    // check the end point pointers
    unsigned short ipt, ii;
    // make sure that the Chg is set correctly
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].UseHit.size() != tj.Pts[ipt].Hits.size()) {
        std::cout<<"UseHit size != Hits size\n";
        fQuitAlg = true;
        return;
      }
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].UseHit.size(); ++ii)
        if(tj.Pts[ipt].UseHit[ii]) tj.Pts[ipt].Chg += fHits[tj.Pts[ipt].Hits[ii]]->Integral();
    } // ipt
    
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[0] = ipt;
        break;
      }
    }
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      ipt = tj.Pts.size() - 1 - ii;
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[1] = ipt;
        break;
      }
    }
  } // SetEndPoints


  ////////////////////////////////////////////////
  void TrajClusterAlg::MakeAllTrajClusters()
  {
    // Make clusters from all trajectories in allTraj
    
    ClusterStore cls;
    tcl.clear();
    inClus.resize(fHits.size());
    unsigned int iht;
    for(iht = 0; iht < inClus.size(); ++iht) inClus[iht] = 0;
    
//    if(prt) mf::LogVerbatim("TC")<<"MakeAllTrajClusters: allTraj size "<<allTraj.size();
    mf::LogVerbatim("TC")<<"MakeAllTrajClusters: allTraj size "<<allTraj.size();
    
    CheckInTraj("MATC");
    if(fQuitAlg) return;
    
    unsigned short itj, endPt0, endPt1;
    
    std::vector<unsigned int> tHits;
    
    // Make one cluster for each trajectory. The indexing of trajectory parents
    // should map directly to cluster parents
    for(itj = 0; itj < allTraj.size(); ++itj) {
      if(allTraj[itj].AlgMod[kKilled]) continue;
      if(allTraj[itj].StepDir > 0) ReverseTraj(allTraj[itj]);
      // ensure that the endPts are correct
      SetEndPoints(allTraj[itj]);
      // some sort of error occurred
      if(allTraj[itj].EndPt[0] >= allTraj[itj].EndPt[1]) continue;
      cls.ID = allTraj[itj].ID;
      cls.CTP = allTraj[itj].CTP;
      cls.PDG = allTraj[itj].PDG;
      cls.ParentCluster = allTraj[itj].ParentTraj;
      endPt0 = allTraj[itj].EndPt[0];
      cls.BeginWir = allTraj[itj].Pts[endPt0].Pos[0];
      cls.BeginTim = allTraj[itj].Pts[endPt0].Pos[1] / fScaleF;
      cls.BeginAng = allTraj[itj].Pts[endPt0].Ang;
      cls.BeginChg = allTraj[itj].Pts[endPt0].Chg;
      cls.BeginVtx = allTraj[itj].Vtx[0];
      endPt1 = allTraj[itj].EndPt[1];
      cls.EndWir = allTraj[itj].Pts[endPt1].Pos[0];
      cls.EndTim = allTraj[itj].Pts[endPt1].Pos[1] / fScaleF;
      cls.EndAng = allTraj[itj].Pts[endPt1].Ang;
      cls.EndChg = allTraj[itj].Pts[endPt1].Chg;
      cls.EndVtx = allTraj[itj].Vtx[1];
      PutTrajHitsInVector(allTraj[itj], true, tHits);
      if(tHits.size() == 0) {
        mf::LogWarning("TC")<<"MakeAllTrajClusters: No hits found in trajectory "<<itj<<" so skip it";
        continue;
      } // error
      cls.tclhits = tHits;
      // Set the traj info
      allTraj[itj].ClusterIndex = tcl.size();
      tcl.push_back(cls);
      // do some checking
      for(iht = 0; iht < cls.tclhits.size(); ++iht) {
        unsigned int hit = cls.tclhits[iht];
        if(fHits[hit]->WireID().Plane == fPlane) continue;
        if(fHits[hit]->WireID().Cryostat == fCstat) continue;
        if(fHits[hit]->WireID().TPC == fTpc) continue;
        mf::LogError("TC")<<"MakeAllTrajClusters: Bad hit CTP in itj "<<itj;
        fQuitAlg = true;
        return;
      } //iht
      for(iht = 0; iht < cls.tclhits.size(); ++iht) inClus[cls.tclhits[iht]] = cls.ID;
    } // itj

    PrintAllTraj(USHRT_MAX, 0);
    PrintClusters();
    
  } // MakeAllTrajClusters

  //////////////////////////////////////////
  void TrajClusterAlg::PrintAllTraj(unsigned short itj, unsigned short ipt)
  {
    
    mf::LogVerbatim myprt("TC");

    if(vtx3.size() > 0) {
      // print out 3D vertices
      myprt<<"****** 3D vertices ******************************************__2DVtx_Indx__*******\n";
      myprt<<"Vtx  Cstat  TPC Proc     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for(unsigned short iv = 0; iv < vtx3.size(); ++iv) {
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(7)<<vtx3[iv].CStat;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].TPC;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].ProcCode;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].X;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].Y;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].Z;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].XErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].YErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].ZErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[0];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[1];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[2];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Wire;
        if(vtx3[iv].Wire < 0) {
          myprt<<"    Matched in all planes";
        } else {
          myprt<<"    Incomplete";
        }
        myprt<<"\n";
      }
    } // vtx3.size
    
    if(vtx.size() > 0) {
      // print out 2D vertices
      myprt<<"************ 2D vertices ************\n";
      myprt<<"Vtx   CTP    wire     error   tick     error  ChiDOF  NCl  topo  traj IDs\n";
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(fDebugPlane < 3 && fDebugPlane != (int)vtx[iv].CTP) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<vtx[iv].CTP;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Wire<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].WireErr;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Time/fScaleF<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].TimeErr/fScaleF;
        myprt<<std::right<<std::setw(8)<<vtx[iv].ChiDOF;
        myprt<<std::right<<std::setw(5)<<vtx[iv].NTraj;
        myprt<<std::right<<std::setw(6)<<vtx[iv].Topo;
        myprt<<"    ";
        // display the traj indices
        for(unsigned short ii = 0; ii < allTraj.size(); ++ii) {
          if(fDebugPlane < 3 && fDebugPlane != (int)allTraj[ii].CTP) continue;
          if(allTraj[ii].AlgMod[kKilled]) continue;
          for(unsigned short end = 0; end < 2; ++end)
            if(allTraj[ii].Vtx[end] == (short)iv) myprt<<std::right<<std::setw(4)<<allTraj[ii].ID<<"_"<<end;
        }
        myprt<<"\n";
      } // iv
    } // vtx.size
    
     // Print all trajectories in allTraj if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    unsigned short endPt;
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
       std::vector<unsigned int> tmp;
      myprt<<"TRJ  ID CTP Pass Pts frm  to     W:Tick     Ang   AveQ     W:T        Ang   AveQ  Hits/TP __Vtx__ PDG Parent TRuPDG  Prnt   KE  \n";
      for(unsigned short ii = 0; ii < allTraj.size(); ++ii) {
        if(fDebugPlane >=0 && fDebugPlane < 3 && (unsigned short)fDebugPlane != allTraj[ii].CTP) continue;
        myprt<<"TRJ"<<std::fixed;
        Trajectory tj = allTraj[ii];
        myprt<<std::setw(4)<<tj.ID;
        myprt<<std::setw(3)<<tj.CTP;
        myprt<<std::setw(5)<<tj.Pass;
        myprt<<std::setw(5)<<tj.Pts.size();
        myprt<<std::setw(4)<<tj.EndPt[0];
        myprt<<std::setw(4)<<tj.EndPt[1];
        endPt = tj.EndPt[0];
        TrajPoint tp = tj.Pts[endPt];
        unsigned short itick = tp.Pos[1]/fScaleF;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.ChgRat;
        endPt = tj.EndPt[1];
        tp = tj.Pts[endPt];
        itick = tp.Pos[1]/fScaleF;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.ChgRat;
        // find average number of used hits / TP
        PutTrajHitsInVector(allTraj[ii], true, tmp);
        float ave = (float)tmp.size() / (float)tj.Pts.size();
        myprt<<std::setw(8)<<std::setprecision(2)<<ave;
        myprt<<std::setw(4)<<tj.Vtx[0];
        myprt<<std::setw(4)<<tj.Vtx[1];
        myprt<<std::setw(6)<<tj.PDG;
        myprt<<std::setw(6)<<tj.ParentTraj;
        myprt<<std::setw(6)<<tj.TruPDG;
        myprt<<std::setw(6)<<tj.IsPrimary;
        myprt<<std::setw(7)<<(int)tj.TruKE;
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
      } // ii
      return;
    } // itj > allTraj.size()-1
    
    if(itj > allTraj.size()-1) return;
    
    Trajectory tj = allTraj[itj];
    
    mf::LogVerbatim("TC")<<"Print allTraj["<<itj<<"]: ClusterIndex "<<tj.ClusterIndex<<" Vtx[0] "<<tj.Vtx[0]<<" Vtx[1] "<<tj.Vtx[1];
    
    PrintHeader();
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) PrintTrajPoint(ii, tj.StepDir, tj.Pass, tj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
    }
  } // PrintAllTraj

  
  //////////////////////////////////////////
  void TrajClusterAlg::PrintTrajectory(Trajectory const& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    unsigned short first = 0;
    unsigned short last = tj.Pts.size();
    if(tPoint == USHRT_MAX) {
      if(tj.ID == 0) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Work:    ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDG<<" TruPDG "<<tj.TruPDG<<" vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[0];
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      } else {
        mf::LogVerbatim myprt("TC");
        myprt<<"allTraj: ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDG<<" TruPDG "<<tj.TruPDG<<" vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[0];
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      }
      PrintHeader();
      for(unsigned short ipt = first; ipt < last; ++ipt) PrintTrajPoint(ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
    } else {
      // just print one traj point
      if(tPoint > tj.Pts.size() -1) {
        mf::LogVerbatim("TC")<<"Can't print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(tPoint, tj.StepDir, tj.Pass, tj.Pts[tPoint]);
    }
  } // PrintTrajectory
  
  //////////////////////////////////////////
  void TrajClusterAlg::PrintHeader()
  {
    mf::LogVerbatim("TC")<<"TRP  CTP Ind  Stp      W:Tick   Delta  RMS    Ang   Err   Kink  Dir0  Dir1      Q  QRat    AveQ FitChi NTPF  Hits ";
  } // PrintHeader

  ////////////////////////////////////////////////
  void TrajClusterAlg::PrintTrajPoint(unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<"TRP"<<std::fixed;
    myprt<<pass;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(3)<<tp.CTP;
    myprt<<std::setw(4)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(7)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]/fScaleF; // W:T
    if(tp.Pos[1] < 10) myprt<<"  "; if(tp.Pos[1] < 100) myprt<<" "; if(tp.Pos[1] < 1000) myprt<<" ";
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta;
/*
    int itick = (int)(tp.Delta/fScaleF);
//    if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
    myprt<<std::setw(6)<<std::setprecision(1)<<itick; // dTick
*/
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.DeltaRMS;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(7)<<std::setprecision(2)<<tp.KinkAng;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[0];
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[1];
    myprt<<std::setw(7)<<(int)tp.Chg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgRat;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(7)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NTPsFit;
    // print the hits associated with this traj point
    for(unsigned short iht = 0; iht < tp.Hits.size(); ++iht) {
      myprt<<" "<<fHits[tp.Hits[iht]]->WireID().Plane;
      myprt<<":"<<fHits[tp.Hits[iht]]->WireID().Wire;
      if(tp.UseHit[iht]) {
        // Distinguish used hits from nearby hits
        myprt<<":";
      } else {
        myprt<<"x";
      }
      myprt<<(int)fHits[tp.Hits[iht]]->PeakTime();
      myprt<<"_"<<inTraj[tp.Hits[iht]];
    } // iht
  } // PrintTrajPoint

  ////////////////////////////////////////////////
  void TrajClusterAlg::HitMultipletPosition(unsigned int theHit, float& hitTick, float& deltaRms, float& qtot)
  {
    // returns the charge weighted wire, time position of hits in the multiplet which are within
    // fMultHitSep of iht
    
    // ignore multiplets. Just use hit separation
    float hitSep = fMultHitSep * fHits[theHit]->RMS();
    unsigned int iht;
    unsigned int wire = fHits[theHit]->WireID().Wire;
    unsigned int firsthit = (unsigned int)WireHitRange[wire].first;
    unsigned int lasthit = (unsigned int)WireHitRange[wire].second;
    qtot = 0;
    hitTick = 0;
    std::vector<unsigned int> closeHits;
    for(iht = theHit; iht < lasthit; ++iht) {
      if(fHits[iht]->PeakTime() - fHits[theHit]->PeakTime() > hitSep) break;
      // break if a hit is found that belongs to a trajectory
      if(inTraj[iht] > 0) break;
      qtot += fHits[iht]->Integral();
      hitTick += fHits[iht]->Integral() * fHits[iht]->PeakTime();
      closeHits.push_back(iht);
    } // iht
    for(iht = theHit - 1; iht >= firsthit; --iht) {
      if(fHits[theHit]->PeakTime() - fHits[iht]->PeakTime() > hitSep) break;
      // break if a hit is found that belongs to a trajectory
      if(inTraj[iht] > 0) break;
      qtot += fHits[iht]->Integral();
      hitTick += fHits[iht]->Integral() * fHits[iht]->PeakTime();
      closeHits.push_back(iht);
      if(iht == 0) break;
    } // iht
    if(qtot == 0) return;
    if(closeHits.size() == 1) {
      qtot = 0;
      return;
    }
    
    hitTick /= qtot;
    
    deltaRms = sqrt(HitsTimeErr2(closeHits));
    if(prt) mf::LogVerbatim("TC")<<" hitTick "<<hitTick<<" qtot "<<qtot;

  } // HitMultipletPosition
/*
  ////////////////////////////////////////////////
  void TrajClusterAlg::TrajSeparation(Trajectory& iTj, Trajectory& jTj, std::vector<float>& tSep)
  {
    // Returns a vector with the size of the number of trajectory points on iTj
    // with the separation distance between the closest traj points
    
    tSep.resize(iTj.Pts.size());
    unsigned short ipt, jpt, ibest;
    float dw, dt, sep, best;
    for(ipt = 0; ipt < iTj.Pts.size(); ++ipt) {
 deal with EndPt here
      best = maxSep2;
      ibest = USHRT_MAX;
      for(jpt = 0; jpt < jTj.Pts.size(); ++jpt) {
        dw = jTj.Pts[jpt].Pos[0] - iTj.Pts[ipt].Pos[0];
        dt = jTj.Pts[jpt].Pos[1] - iTj.Pts[ipt].Pos[1];
        sep = dw * dw + dt * dt;
        if(sep < best) {
          best = sep;
          ibest = jpt;
        }
        // decide whether to break out
        if(ibest != USHRT_MAX && sep > best) break;
      } // jpt
      tSep[ipt] = sqrt(best);
//      std::cout<<"Sep "<<ipt<<" "<<best<<"\n";
    } // ipt
    
  } // TrajectorySeparation
*/

  ////////////////////////////////////////////////
  bool TrajClusterAlg::TrajHitsOK(unsigned int iht, unsigned int jht)
  {
    // Hits (assume to be on adjacent wires with wire A > wire B) have an acceptable signal overlap
    
    if(iht > fHits.size() - 1) return false;
    if(jht > fHits.size() - 1) return false;
    
    raw::TDCtick_t hiStartTick = fHits[iht]->StartTick();
    if(fHits[jht]->StartTick() > hiStartTick) hiStartTick = fHits[jht]->StartTick();
    raw::TDCtick_t loEndTick = fHits[iht]->EndTick();
    if(fHits[jht]->EndTick() < loEndTick) loEndTick = fHits[jht]->EndTick();
    // add a tolerance to the StartTick - EndTick overlap
    raw::TDCtick_t tol = 30;
    // expand the tolerance for induction planes
    if(fPlane < geom->Cryostat(fCstat).TPC(fTpc).Nplanes()-1) tol = 40;

    if(fHits[jht]->PeakTime() > fHits[iht]->PeakTime()) {
      // positive slope
      if(loEndTick + tol < hiStartTick) {
//          if(prt) mf::LogVerbatim("TC")<<" bad overlap pos Slope "<<loEndTick<<" > "<<hiStartTick;
        return false;
      }
    } else {
      // negative slope
      if(loEndTick + tol < hiStartTick) {
//          if(prt) mf::LogVerbatim("TC")<<" bad overlap neg Slope "<<loEndTick<<" < "<<hiStartTick;
        return false;
      }
    }
    
    return true;
    
  } // TrajHitsOK
  
  /////////////////////////////////////////
  void TrajClusterAlg::PrintClusters()
  {
    
    // prints clusters to the screen for code development
    mf::LogVerbatim myprt("TC");
    
    if(vtx3.size() > 0) {
      // print out 3D vertices
      myprt<<"****** 3D vertices ******************************************__2DVtx_Indx__*******\n";
      myprt<<"Vtx  Cstat  TPC Proc     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for(unsigned short iv = 0; iv < vtx3.size(); ++iv) {
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(7)<<vtx3[iv].CStat;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].TPC;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].ProcCode;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].X;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].Y;
        myprt<<std::right<<std::setw(8)<<vtx3[iv].Z;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].XErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].YErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].ZErr;
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[0];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[1];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Ptr2D[2];
        myprt<<std::right<<std::setw(5)<<vtx3[iv].Wire;
        if(vtx3[iv].Wire < 0) {
          myprt<<"    Matched in all planes";
        } else {
          myprt<<"    Incomplete";
        }
        myprt<<"\n";
      }
    } // vtx3.size
    
    if(vtx.size() > 0) {
      // print out 2D vertices
      myprt<<"************ 2D vertices ************\n";
      myprt<<"Vtx   CTP    wire     error   tick     error  ChiDOF  NCl  topo  cluster IDs\n";
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(fDebugPlane < 3 && fDebugPlane != (int)vtx[iv].CTP) continue;
        if(vtx[iv].NTraj == 0) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<vtx[iv].CTP;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Wire<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].WireErr;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Time/fScaleF<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].TimeErr/fScaleF;
        myprt<<std::right<<std::setw(8)<<vtx[iv].ChiDOF;
        myprt<<std::right<<std::setw(5)<<vtx[iv].NTraj;
        myprt<<std::right<<std::setw(6)<<vtx[iv].Topo;
        myprt<<"    ";
        // display the cluster IDs
        for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
          if(fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
          if(tcl[ii].ID < 0) continue;
          if(tcl[ii].BeginVtx == (short)iv) myprt<<std::right<<std::setw(4)<<tcl[ii].ID<<"_B";
          if(tcl[ii].EndVtx == (short)iv) myprt<<std::right<<std::setw(4)<<tcl[ii].ID<<"_E";
        }
        myprt<<"\n";
      } // iv
    } // vtx.size
    
    myprt<<"Total number of clusters "<<tcl.size()<<"\n";
    
    float aveRMS, aveRes;
    myprt<<"*************************************** Clusters *********************************************************************\n";
    myprt<<"  ID CTP nht beg_W:T      bAng   bChg end_W:T      eAng   eChg  bVx  eVx aveRMS Qual cnt\n";
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      // print clusters in all planes (fDebugPlane = 3) or in a selected plane
//      if(fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
      myprt<<std::right<<std::setw(4)<<tcl[ii].ID;
      myprt<<std::right<<std::setw(3)<<tcl[ii].CTP;
      myprt<<std::right<<std::setw(5)<<tcl[ii].tclhits.size();
      unsigned short iTime = tcl[ii].BeginTim;
      myprt<<std::right<<std::setw(6)<<(int)(tcl[ii].BeginWir+0.5)<<":"<<iTime;
      if(iTime < 10) {
        myprt<<"   ";
      } else if(iTime < 100) {
        myprt<<"  ";
      } else if(iTime < 1000) myprt<<" ";
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].BeginAng;
      myprt<<std::right<<std::setw(7)<<(int)tcl[ii].BeginChg;
      iTime = tcl[ii].EndTim;
      myprt<<std::right<<std::setw(6)<<(int)(tcl[ii].EndWir+0.5)<<":"<<iTime;
      if(iTime < 10) {
        myprt<<"   ";
      } else if(iTime < 100) {
        myprt<<"  ";
      } else if(iTime < 1000) myprt<<" ";
      myprt<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<tcl[ii].EndAng;
      myprt<<std::right<<std::setw(7)<<(int)tcl[ii].EndChg;
      myprt<<std::right<<std::setw(5)<<tcl[ii].BeginVtx;
      myprt<<std::right<<std::setw(5)<<tcl[ii].EndVtx;
      aveRMS = 0;
      unsigned int iht = 0;
      for(unsigned short jj = 0; jj < tcl[ii].tclhits.size(); ++jj) {
        iht = tcl[ii].tclhits[jj];
        aveRMS += fHits[iht]->RMS();
      }
      aveRMS /= (float)tcl[ii].tclhits.size();
      myprt<<std::right<<std::setw(5)<<std::fixed<<std::setprecision(1)<<aveRMS;
      aveRes = 0;
      // find cluster tracking resolution
      unsigned int hit0, hit1, hit2, cnt = 0;
      float arg;
      for(unsigned short iht = 1; iht < tcl[ii].tclhits.size()-1; ++iht) {
        hit1 = tcl[ii].tclhits[iht];
        hit0 = tcl[ii].tclhits[iht-1];
        hit2 = tcl[ii].tclhits[iht+1];
        // require hits on adjacent wires
        if(fHits[hit1]->WireID().Wire + 1 != fHits[hit0]->WireID().Wire) continue;
        if(fHits[hit2]->WireID().Wire + 1 != fHits[hit1]->WireID().Wire) continue;
        arg = (fHits[hit0]->PeakTime() + fHits[hit2]->PeakTime())/2 - fHits[hit1]->PeakTime();
        aveRes += arg * arg;
        ++cnt;
      }
      if(cnt > 1) {
        aveRes /= (float)cnt;
        aveRes = sqrt(aveRes);
        // convert to a quality factor
        aveRes /= (aveRMS * fHitErrFac);
        myprt<<std::right<<std::setw(6)<<std::fixed<<std::setprecision(1)<<aveRes;
        myprt<<std::right<<std::setw(5)<<std::fixed<<cnt;
      } else {
        myprt<<"    NA";
        myprt<<std::right<<std::setw(5)<<std::fixed<<cnt;
      }
      myprt<<"\n";
    } // ii
    
  } // PrintClusters()
  
  /////////////////////////////////////////
  void TrajClusterAlg::MakeBareTrajPoint(unsigned int fromHit, unsigned int toHit, TrajPoint& tp)
  {
    CTP_t tCTP = EncodeCTP(fHits[fromHit]->WireID());
    MakeBareTrajPoint((float)fHits[fromHit]->WireID().Wire, fHits[fromHit]->PeakTime(),
                      (float)fHits[toHit]->WireID().Wire,   fHits[toHit]->PeakTime(), tCTP, tp);
    
  } // MakeBareTrajPoint

  /////////////////////////////////////////
  void TrajClusterAlg::MakeBareTrajPoint(float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp)
  {
    tp.CTP = tCTP;
    tp.Pos[0] = fromWire;
    tp.Pos[1] = fScaleF * fromTick;
    tp.Dir[0] = toWire - fromWire;
    tp.Dir[1] = fScaleF * (toTick - fromTick);
    float norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    
  } // MakeBareTrajPoint

  /////////////////////////////////////////
  float TrajClusterAlg::TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Calculates the angle of a line between two TPs
    float dw = tp2.Pos[0] - tp1.Pos[0];
    float dt = tp2.Pos[1] - tp1.Pos[1];
    return atan2(dw, dt);
  } // TwoTPAngle

  /////////////////////////////////////////
  bool TrajClusterAlg::SkipHighMultHitCombo(unsigned int iht, unsigned int jht)
  {
    // Return true if iht and jht are both in a multiplet but have the wrong local index to start a trajectory
    if(fHits[iht]->Multiplicity() < 3) return false;
    if(fHits[jht]->Multiplicity() < 3) return false;
    
    if(jht > iht && fHits[jht]->StartTick() > fHits[iht]->StartTick()) {
      // "positive slope" as visualized in the event display
      // ^    -
      // |    -
      // t    -
      // i   --
      // m   -
      // e   -
      //     ij
      // wire ->
      if(fHits[iht]->LocalIndex() != 0) return true;
      if(fHits[jht]->LocalIndex() != 0) return true;
    } else {
      // "negative slope"
      if(fHits[iht]->LocalIndex() != fHits[iht]->Multiplicity()-1) return true;
      if(fHits[jht]->LocalIndex() != fHits[jht]->Multiplicity()-1) return true;
    }
    return false;
  } // SkipHighMultHitCombo
  
  /////////////////////////////////////////
  bool TrajClusterAlg::LargeHitSep(unsigned int iht, unsigned int jht)
  {
    // Returns true if the time separation of two hits exceeds fMultHitSep
    if(fHits[iht]->WireID().Wire != fHits[jht]->WireID().Wire) return true;
    float minrms = fHits[iht]->RMS();
    if(fHits[jht]->RMS() < minrms) minrms = fHits[jht]->RMS();
    float hitsep = std::abs(fHits[iht]->RMS()) - fHits[jht]->RMS() / minrms;
    return (hitsep > fMultHitSep);
    
  } // LargeHitSep
 
  
  /////////////////////////////////////////
  bool TrajClusterAlg::IsLargeAngle(TrajPoint const& tp)
  {
    // standard criterion for using large angle cuts
    return (std::abs(tp.Dir[0]) < fLargeAngle);
  } // IsLargeAngle
  
  /////////////////////////////////////////
  bool TrajClusterAlg::SignalAtTp(TrajPoint const& tp)
  {
    // Returns true if there a is wire signal at tp
    unsigned int wire = tp.Pos[0] + 0.5;
    // Assume dead wires have a signal
    if(WireHitRange[wire].first == -1) return true;
    raw::TDCtick_t rawProjTick = (float)(tp.Pos[1] / fScaleF);
    unsigned int firsthit = (unsigned int)WireHitRange[wire].first;
    unsigned int lasthit = (unsigned int)WireHitRange[wire].second;
    for(unsigned int iht = firsthit; iht < lasthit; ++iht) {
      if(rawProjTick > fHits[iht]->StartTick() && rawProjTick < fHits[iht]->EndTick()) return true;
    } // iht
    return false;
  } // SignalAtTp
  
  //////////////////////////////////////////
  unsigned short TrajClusterAlg::NumPtsWithCharge(Trajectory& tj)
  {
    unsigned short ntp = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]+1; ++ipt) if(tj.Pts[ipt].Chg > 0) ++ntp;
    return ntp;
  } // NumTPsWithCharge
  
  //////////////////////////////////////////
  unsigned short TrajClusterAlg::NumUsedHits(TrajPoint& tp)
  {
    // Counts the number of used hits in tp
    unsigned short nused = 0;
    if(tp.Hits.size() == 0) return nused;
    for(unsigned short ii = 0; ii < tp.UseHit.size(); ++ii) if(tp.UseHit[ii]) ++nused;
    return nused;
  } // NumUsedHits
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::HitSanityCheck()
  {
    // checks to see if the hits are properly sorted
    
    std::cout<<"Plane "<<fPlane<<" sanity check. Wire range "<<fFirstWire<<" "<<fLastWire;
    unsigned int nhts = 0;
    for(unsigned int wire = fFirstWire; wire < fLastWire; ++wire) {
      if(WireHitRange[wire].first < 0) continue;
      unsigned int fhit = WireHitRange[wire].first;
      unsigned int lhit = WireHitRange[wire].second;
      for(unsigned int hit = fhit; hit < lhit; ++hit) {
        ++nhts;
        if(fHits[hit]->WireID().Wire != wire) {
          std::cout<<"Bad wire "<<hit<<" "<<fHits[hit]->WireID().Wire<<" "<<wire<<"\n";
          return;
        } // check wire
        if(fHits[hit]->WireID().Plane != fPlane) {
          std::cout<<"Bad plane "<<hit<<" "<<fHits[hit]->WireID().Plane<<" "<<fPlane<<"\n";
          return;
        } // check plane
        if(fHits[hit]->WireID().TPC != fTpc) {
          std::cout<<"Bad tpc "<<hit<<" "<<fHits[hit]->WireID().TPC<<" "<<fTpc<<"\n";
          return;
        } // check tpc
      } // hit
    } // wire
    std::cout<<" is OK. nhits "<<nhts<<"\n";
    
  } // HitSanityCheck

  /////////////////////////////////////////
  std::string TrajClusterAlg::PrintHit(unsigned int iht)
  {
    if(iht > fHits.size() - 1) return "Bad Hit";
    return std::to_string(fHits[iht]->WireID().Wire) + ":" + std::to_string((int)fHits[iht]->PeakTime());
  } // PrintHit
  
  /////////////////////////////////////////
  std::string TrajClusterAlg::PrintPos(TrajPoint const& tp)
  {
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    unsigned int time = std::nearbyint(tp.Pos[1]/fScaleF);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

  /////////////////////////////////////////
  bool TrajClusterAlg::SignalPresent(TrajPoint const& tp)
  {
    return SignalPresent(tp.Pos[1], tp.Pos[0], tp.Pos[1], tp.Pos[0]);
  }

  /////////////////////////////////////////
  bool TrajClusterAlg::SignalPresent(float wire1, float time1, TrajPoint const& tp)
  {
    unsigned int w1 = std::nearbyint(wire1);
    unsigned int w2 = std::nearbyint(tp.Pos[0]);
    return SignalPresent(w1, time1, w2, tp.Pos[1]);
  }

  /////////////////////////////////////////
  bool TrajClusterAlg::SignalPresent(float wire1, float time1, float wire2, float time2)
  {
    unsigned int w1 = std::nearbyint(wire1);
    unsigned int w2 = std::nearbyint(wire2);
    return SignalPresent(w1, time1, w2, time2);
  } // SignalPresent

  /////////////////////////////////////////
  bool TrajClusterAlg::SignalPresent(unsigned int wire1, float time1, unsigned int wire2, float time2)
  {
    // returns  true if there is a signal on the line between (wire1, time1) and (wire2, time2).
    
    // Gaussian amplitude in bins of size 0.15
    const float gausAmp[20] = {1, 0.99, 0.96, 0.90, 0.84, 0.75, 0.67, 0.58, 0.49, 0.40, 0.32, 0.26, 0.20, 0.15, 0.11, 0.08, 0.06, 0.04, 0.03, 0.02};
    
    // convert time to tick
    time1 /= fScaleF;
    time2 /= fScaleF;
//    mf::LogVerbatim("TC")<<"SignalPresent: check "<<wire1<<":"<<(int)time1<<" to "<<wire2<<":"<<(int)time2;
    
    // get the begin and end right
    unsigned int wireb = wire1;
    float timeb = time1;
    unsigned int wiree = wire2;
    float timee = time2;
    // swap them?
    if(wiree > wireb) {
      wireb = wire2;
      timeb = time2;
      wiree = wire1;
      timee = time1;
    }
    if(wiree < fFirstWire || wiree > fLastWire) return false;
    if(wireb < fFirstWire || wireb > fLastWire) return false;
    
    int wire0 = wiree;
    // checking a single wire?
    float slope = 0;
    bool oneWire = false;
    float prTime, prTimeLo = 0, prTimeHi = 0;
    if(wireb == wiree) {
      oneWire = true;
      if(time1 < time2) {
        prTimeLo = time1;
        prTimeHi = time2;
      } else {
        prTimeLo = time2;
        prTimeHi = time1;
      }
    } else {
      slope = (timeb - timee) / (wireb - wiree);
    }
    
    int bin;
    for(unsigned int wire = wiree; wire < wireb + 1; ++wire) {
      if(oneWire) {
        prTime = (prTimeLo + prTimeHi) / 2;
      } else {
        prTime = timee + (wire - wire0) * slope;
      }
      // skip dead wires
      if(WireHitRange[wire].first == -1) continue;
      // no hits on this wire
      if(WireHitRange[wire].first == -2) return false;
      unsigned int firsthit = WireHitRange[wire].first;
      unsigned int lasthit = WireHitRange[wire].second;
//      mf::LogVerbatim("TC")<<" wire "<<wire<<" Hit range "<<firsthit<<" "<<lasthit<<" prTime "<<prTime;
      float amp = 0;
      for(unsigned int khit = firsthit; khit < lasthit; ++khit) {
//        mf::LogVerbatim("TC")<<"  hit "<<PrintHit(khit)<<" rms "<<fHits[khit]->RMS()<<" amp "<<(int)fHits[khit]->PeakAmplitude()<<" StartTick "<<fHits[khit]->StartTick()<<" EndTick "<<fHits[khit]->EndTick();
        if(oneWire) {
          // TODO: This sometimes doesn't work with overlapping hits
//            if(prTimeHi > fHits[khit].EndTick()) continue;
//            if(prTimeLo < fHits[khit].StartTick()) continue;
          // A not totally satisfactory solution
          if(prTime < fHits[khit]->StartTick()) continue;
          if(prTime > fHits[khit]->EndTick()) continue;
          return true;
        } else {
          // skip checking if we are far away from prTime on the positive side
          if(fHits[khit]->PeakTime() - prTime > 500) continue;
          bin = std::abs(fHits[khit]->PeakTime() - prTime) / fHits[khit]->RMS();
          bin /= 0.15;
          if(bin > 19) continue;
          if(bin < 0) continue;
//          mf::LogVerbatim("CC")<<"  bin "<<bin<<" add "<<fHits[khit]->PeakAmplitude() * gausAmp[bin]<<" to amp "<<amp;
          // add amplitude from all hits
          amp += fHits[khit]->PeakAmplitude() * gausAmp[bin];
        }
      } // khit
//      mf::LogVerbatim("TC")<<"Amp "<<amp<<" fMinAmp "<<fMinAmp;
      if(amp < fMinAmp) return false;
    } // wire
    return true;
    
  } // SignalPresent
  
  //////////////////////////////////
  void TrajClusterAlg::GetHitRange()
  {
    // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
    // Hits must have been sorted by increasing wire number
    fFirstHit = 0;
    geo::PlaneID planeID = DecodeCTP(fCTP);
    fNumWires = geom->Nwires(planeID.Plane, planeID.TPC, planeID.Cryostat);
    fMaxWire = fNumWires + 0.5;
//    fMaxTick = (float)detprop->NumberTimeSamples();
    fMaxTime = (float)detprop->NumberTimeSamples() * fScaleF;
    WireHitRange.resize(fNumWires + 1);
    
    // These will be re-defined later
    fFirstWire = 0;
    fLastWire = 0;
    
    unsigned int wire, iht;
    unsigned int nHitInPlane;
    std::pair<int, int> flag;
    
    // Define the "no hits on wire" condition
    flag.first = -2; flag.second = -2;
    for(auto& apair : WireHitRange) apair = flag;
    
    nHitInPlane = 0;
    
    std::vector<bool> firsthit;
    firsthit.resize(fNumWires+1, true);
    bool firstwire = true;
    for(iht = 0; iht < fHits.size(); ++iht) {
      if(fHits[iht]->WireID().TPC != planeID.TPC) continue;
      if(fHits[iht]->WireID().Cryostat != planeID.Cryostat) continue;
      if(fHits[iht]->WireID().Plane != planeID.Plane) continue;
      wire = fHits[iht]->WireID().Wire;
      // define the first hit start index in this TPC, Plane
      if(firsthit[wire]) {
        WireHitRange[wire].first = iht;
        firsthit[wire] = false;
      }
      if(firstwire){
        fFirstWire = wire;
        firstwire = false;
      }
      WireHitRange[wire].second = iht+1;
      fLastWire = wire+1;
      ++nHitInPlane;
    }
    // overwrite with the "dead wires" condition
    lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    flag.first = -1; flag.second = -1;
    for(wire = 0; wire < fNumWires; ++wire) {
      raw::ChannelID_t chan = geom->PlaneWireToChannel
      ((int)planeID.Plane,(int)wire,(int)planeID.TPC,(int)planeID.Cryostat);
      if(channelStatus.IsBad(chan)) WireHitRange[wire] = flag;
    }
    
    unsigned int firstHit, lastHit;
    for(wire = 0; wire < fNumWires; ++wire) {
      // ignore dead wires and wires with no hits
      if(WireHitRange[wire].first < 0) continue;
      firstHit = WireHitRange[wire].first;
      lastHit = WireHitRange[wire].second;
      for(iht = firstHit; iht < lastHit; ++iht) {
        if(fHits[iht]->WireID().Wire != wire) {
          mf::LogWarning("TC")<<"Bad WireHitRange wire "<<fHits[iht]->WireID().Wire<<" != "<<wire;
          fQuitAlg = true;
          return;
        }
      } // iht
    } // wire
    
  } // GetHitRange()
  
  //////////////////////////////////////////
  float TrajClusterAlg::DeadWireCount(float inWirePos1, float inWirePos2)
  {
    unsigned int inWire1 = std::nearbyint(inWirePos1);
    unsigned int inWire2 = std::nearbyint(inWirePos2);
    if(inWire1 > inWire2) {
      // put in increasing order
      unsigned int tmp = inWire1;
      inWire1 = inWire2;
      inWire2 = tmp;
    } // inWire1 > inWire2
    ++inWire2;
    unsigned int wire, ndead = 0;
    for(wire = inWire1; wire < inWire2; ++wire) if(WireHitRange[wire].first == -1) ++ndead;
    return ndead;
  } // DeadWireCount

  //////////////////////////////////////////
  void TrajClusterAlg::CheckHitClusterAssociations()
  {
    // check hit - cluster associations
    
    if(fHits.size() != inClus.size()) {
      mf::LogError("TC")<<"CHCA: Sizes wrong "<<fHits.size()<<" "<<inClus.size();
      fQuitAlg = true;
      return;
    }
    
    unsigned int iht;
    short clID;
    
    // check cluster -> hit association
    for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
      if(tcl[icl].ID < 0) continue;
      clID = tcl[icl].ID;
      for(unsigned short ii = 0; ii < tcl[icl].tclhits.size(); ++ii) {
        iht = tcl[icl].tclhits[ii];
        if(iht > fHits.size() - 1) {
          mf::LogError("CC")<<"CHCA: Bad tclhits index "<<iht<<" fHits size "<<fHits.size();
          fQuitAlg = true;
          return;
        } // iht > fHits.size() - 1
        if(inClus[iht] != clID) {
          mf::LogError("TC")<<"CHCA: Bad cluster -> hit association. clID "<<clID<<" hit "<<PrintHit(iht)<<" inClus "<<inClus[iht]<<" CTP "<<tcl[icl].CTP;
          FindHit("CHCA ", iht);
          fQuitAlg = true;
          return;
        }
      } // ii
    } // icl
    
    // check hit -> cluster association
    unsigned short icl;
    for(iht = 0; iht < fHits.size(); ++iht) {
      if(inClus[iht] <= 0) continue;
      icl = inClus[iht] - 1;
      // see if the cluster is obsolete
      if(tcl[icl].ID < 0) {
        mf::LogError("TC")<<"CHCA: Hit "<<PrintHit(iht)<<" associated with an obsolete cluster tcl[icl].ID "<<tcl[icl].ID;
        fQuitAlg = true;
        return;
      }
      if (std::find(tcl[icl].tclhits.begin(), tcl[icl].tclhits.end(), iht) == tcl[icl].tclhits.end()) {
        mf::LogError("TC")<<"CHCA: Bad hit -> cluster association. hit "<<PrintHit(iht)<<" inClus "<<inClus[iht];
        fQuitAlg = true;
        return;
      }
    } // iht
    
  } // CheckHitClusterAssociations()
  
  
  //////////////////////////////////////////
  void TrajClusterAlg::SplitAllTraj(unsigned short itj, unsigned short pos, unsigned short ivx)
  {
    // Splits the trajectory itj in the allTraj vector into two trajectories at position pos. Splits
    // the trajectory and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
    
    fSplitTrajOK = false;
    if(itj > allTraj.size()-1) return;
    if(pos < allTraj[itj].EndPt[0] + 1 || pos > allTraj[itj].EndPt[1] - 1) return;
    
    Trajectory& tj = allTraj[itj];
    
    // ensure that there will be at least 3 TPs on each trajectory
    unsigned short ipt, ii, ntp = 0;
    for(ipt = 0; ipt < pos; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) return;
    ntp = 0;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) return;
    
    // make a copy
    Trajectory newTj = allTraj[itj];
    newTj.ID = allTraj.size() + 1;
    
    // Leave the first section of tj in place. Re-assign the hits
    // to the new trajectory
    unsigned int iht;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(!tj.Pts[ipt].UseHit[ii]) continue;
        iht = tj.Pts[ipt].Hits[ii];
        // This shouldn't happen but check anyway
        if(inTraj[iht] != tj.ID) continue;
        inTraj[iht] = newTj.ID;
        tj.Pts[ipt].UseHit[ii] = false;
      } // ii
    } // ipt
    SetEndPoints(tj);
    if(ivx != USHRT_MAX) tj.Vtx[1] = ivx;
    tj.AlgMod[kSplitTraj] = true;
    if(vtxPrt) {
      mf::LogVerbatim("TC")<<"SplitAllTraj: itj "<<tj.ID<<" EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
      PrintTrajectory(allTraj[itj], USHRT_MAX);
    }
    
    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if(eraseSize > newTj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"SplitAllTraj: Bad erase size ";
      fSplitTrajOK = false;
      return;
    }
    
    // erase the TPs at the beginning of the new trajectory
    //    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + pos + 1);
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for(ipt = 0; ipt < 3; ++ipt) {
      for(ii = 0; ii < newTj.Pts[ipt].UseHit.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
    } // ipt
    SetEndPoints(newTj);
    if(ivx != USHRT_MAX) newTj.Vtx[0] = ivx;
    newTj.AlgMod[kSplitTraj] = true;
    allTraj.push_back(newTj);
    if(vtxPrt) {
      mf::LogVerbatim("TC")<<"SplitAllTraj: NewTj "<<newTj.ID<<" EndPts "<<newTj.EndPt[0]<<" to "<<newTj.EndPt[1];
      PrintTrajectory(newTj, USHRT_MAX);
    }
    fSplitTrajOK = true;
    
  } // SplitAllTraj
  
  //////////////////////////////////////////
  void TrajClusterAlg::TrajPointTrajDOCA(TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep)
  {
    // Finds the point, ipt, on trajectory tj that is closest to trajpoint tp
    float best = minSep * minSep;
    closePt = 0;
    float dw, dt, dp2;
    unsigned short ipt;
    for(ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
      dw = tj.Pts[ipt].Pos[0] - tp.Pos[0];
      dt = tj.Pts[ipt].Pos[1] - tp.Pos[1];
      dp2 = dw * dw + dt * dt;
      if(dp2 < best) {
        best = dp2;
        closePt = ipt;
      }
    } // ipt
    minSep = sqrt(best);
  } // TrajPointTrajDOCA
  
  //////////////////////////////////////////
  void TrajClusterAlg::TrajTrajDOCA(Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
  {
    // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
    float best = minSep * minSep;
    ipt1 = 0; ipt2 = 0;
    float dw, dt, dp2;
    unsigned short i1, i2;
    for(i1 = tj1.EndPt[0]; i1 < tj1.EndPt[1] + 1; ++i1) {
      for(i2 = tj2.EndPt[0]; i2 < tj2.EndPt[1] + 1; ++i2) {
        // TODO: What about TPs with UseHit = false?
        dw = tj1.Pts[i1].Pos[0] - tj2.Pts[i2].Pos[0];
        dt = tj1.Pts[i1].Pos[1] - tj2.Pts[i2].Pos[1];
        dp2 = dw * dw + dt * dt;
        if(dp2 < best) {
          best = dp2;
          ipt1 = i1;
          ipt2 = i2;
        }
      } // i2
    } // i1
    minSep = sqrt(best);
  } // TrajTrajDOCA
  
  //////////////////////////////////////////
  float TrajClusterAlg::TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2)
  {
    // returns the separation^2 between the hit position of two trajectory points
    float dw = tp1.HitPos[0] - tp2.HitPos[0];
    float dt = tp1.HitPos[1] - tp2.HitPos[1];
    return dw * dw + dt * dt;
  } // TrajPointHitSep2
  
  //////////////////////////////////////////
  float TrajClusterAlg::HitSep2(unsigned int iht, unsigned int jht)
  {
    // returns the separation^2 between two hits in WSE units
    if(iht > fHits.size()-1 || jht > fHits.size()-1) return 1E6;
    float dw = (float)fHits[iht]->WireID().Wire - (float)fHits[jht]->WireID().Wire;
    float dt = (fHits[iht]->PeakTime() - fHits[jht]->PeakTime()) * fScaleF;
    return dw * dw + dt * dt;
  } // TrajPointHitSep2
  
  //////////////////////////////////////////
  float TrajClusterAlg::PointTrajSep2(float wire, float time, TrajPoint const& tp)
  {
    float dw = wire - tp.Pos[0];
    float dt = time - tp.Pos[1];
    return dw * dw + dt * dt;
  }
  
  //////////////////////////////////////////
  float TrajClusterAlg::PointTrajDOCA(unsigned int iht, TrajPoint const& tp)
  {
    float wire = fHits[iht]->WireID().Wire;
    float time = fHits[iht]->PeakTime() * fScaleF;
    return sqrt(PointTrajDOCA2(wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float TrajClusterAlg::PointTrajDOCA(float wire, float time, TrajPoint const& tp)
  {
    return sqrt(PointTrajDOCA2(wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float TrajClusterAlg::PointTrajDOCA2(float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach squared between a (wire, time(WSE)) point
    // and a trajectory point
    
    float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
    float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    float dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return (dw * dw + dt * dt);
    
  } // PointTrajDOCA2
  
  //////////////////////////////////////////
  void TrajClusterAlg::TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y)
  {
    // returns the intersection position, intPos, of two trajectory points
    
    x = -9999; y = -9999;
    
    double arg1 = tp1.Pos[0] * tp1.Dir[1] - tp1.Pos[1] * tp1.Dir[0];
    double arg2 = tp2.Pos[0] * tp1.Dir[1] - tp2.Pos[1] * tp1.Dir[0];
    double arg3 = tp2.Dir[0] * tp1.Dir[1] - tp2.Dir[1] * tp1.Dir[0];
    if(arg3 == 0) return;
    double s = (arg1 - arg2) / arg3;
    
    x = (float)(tp2.Pos[0] + s * tp2.Dir[0]);
    y = (float)(tp2.Pos[1] + s * tp2.Dir[1]);
    
  } // TrajIntersection
  
  //////////////////////////////////////////
  float TrajClusterAlg::TrajLength(Trajectory& tj)
  {
    float len = 0, dx, dy;
    unsigned short ipt;
    unsigned short prevPt = tj.EndPt[0];
    for(ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1] + 1; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dx = tj.Pts[ipt].Pos[0] - tj.Pts[prevPt].Pos[0];
      dy = tj.Pts[ipt].Pos[1] - tj.Pts[prevPt].Pos[1];
      len += sqrt(dx * dx + dy * dy);
      prevPt = ipt;
    }
    return len;
  } // TrajLength
  
  //////////////////////////////////////////
  float TrajClusterAlg::TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Returns the separation distance between two trajectory points
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation
  
  //////////////////////////////////////////
  void TrajClusterAlg::TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& Distance)
  {
    // find the closest approach between a trajectory tj and a point (x,y). Returns
    // the index of the closest trajectory point and the distance
    
    float dx, dy, dist, best = 1E6;
    closePt = 0;
    Distance = best;
    
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dx = tj.Pts[ipt].Pos[0] - x;
      dy = tj.Pts[ipt].Pos[1] - y;
      dist = dx * dx + dy * dy;
      if(dist < best) {
        best = dist;
        closePt = ipt;
      }
      // TODO is this wise?
      //      if(dist > best) break;
    } // ipt
    
    Distance = sqrt(best);
    
  } // TrajClosestApproach



} // namespace cluster