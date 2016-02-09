//////////////////////////////////////////////////////////////////////
///
/// Step crawling code used by TrajClusterAlg
///
/// Bruce Baller, baller@fnal.gov
///
///
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"

// LArSoft libraries
#include "RecoBase/Hit.h"
#include "RecoAlg/TrajClusterAlg.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

// TEMP for TagAllTraj
#include "MCCheater/BackTracker.h"

class TH1F;
class TH2F;

struct SortEntry{
  int index;
  int length;
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
 
    fHitFinderModuleLabel = pset.get<art::InputTag>("HitFinderModuleLabel");
    fMode                 = pset.get< short >("Mode", 0); // Default is don't use it
    fHitErrFac            = pset.get< float >("HitErrFac", 0.4);
    fLargeAngle           = pset.get< float >("LargeAngle", 80);
    fNPtsAve              = pset.get< short >("NPtsAve", 20);
    fMinNPtsFit           = pset.get< std::vector<unsigned short >>("MinNPtsFit");
    fMinPts               = pset.get< unsigned short >("MinPts", 20);
    fMaxChi               = pset.get< float >("MaxChi", 10);
    fMultHitSep           = pset.get< float >("MultHitSep", 2.5);
    fTP3ChiCut            = pset.get< float >("TP3ChiCut", 2.);
    fChgDiffCut           = pset.get< float >("ChgDiffCut", 10);
    fKinkAngCut           = pset.get< float >("KinkAngCut", 0.4);
    fMaxWireSkip          = pset.get< float >("MaxWireSkip", 1);
    fMaxDeltaJump         = pset.get< float >("MaxDeltaJump", 10);
    fProjectionErrFactor  = pset.get< float >("ProjectionErrFactor", 2);
    fStudyMode            = pset.get< bool  >("StudyMode", false);
    fShowerStudy          = pset.get< bool  >("ShowerStudy", false);
    fTagAllTraj           = pset.get< bool  >("TagAllTraj", false);
    fMaxTrajSep           = pset.get< float >("MaxTrajSep", 4);
    fShowerPrtPlane       = pset.get< short >("ShowerPrtPlane", -1);
    fFindTrajVertices     = pset.get< bool  >("FindTrajVertices", false);
    
    fDebugPlane         = pset.get< int  >("DebugPlane", -1);
    fDebugWire          = pset.get< int  >("DebugWire", -1);
    fDebugHit           = pset.get< int  >("DebugHit", -1);
    
    // convert angle (degrees) into a direction cosine cut in the wire coordinate
    // It should be in the range 0 < fLargeAngle < 90
    fLargeAngle = cos(fLargeAngle * M_PI / 180);
    
    // convert the max traj separation into a separation^2
    fMaxTrajSep *= fMaxTrajSep;
    
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
    art::Handle< std::vector<recob::Hit> > hitVecHandle;
    evt.getByLabel(fHitFinderModuleLabel, hitVecHandle);

    fHits.resize(hitVecHandle->size());
    if(fHits.size() == 0) return;
    
    for (unsigned int iht = 0; iht < fHits.size(); iht++) fHits[iht] = art::Ptr< recob::Hit>(hitVecHandle, iht);
    
    ClearResults();
    // set all hits to the available state
    inTraj.resize(fHits.size(), 0);

    // TODO Insert hit sorting code here
    
    std::cout<<"Event "<<evt.event()<<"\n";
    
    for (geo::TPCID const& tpcid: geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      for(fPlane = 0; fPlane < TPC.Nplanes(); ++fPlane) {
        WireHitRange.clear();
        // define a code to ensure clusters are compared within the same plane
        fCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, fPlane);
        fCstat = tpcid.Cryostat;
        fTpc = tpcid.TPC;
        // fill the WireHitRange vector with first/last hit on each wire
        // dead wires and wires with no hits are flagged < 0
        GetHitRange();
        // no hits on this plane?
        if(fFirstWire == fLastWire) continue;
        // Calculate some default values for hits
        FindDefaults();
        // This assumes (reasonably...) that all planes have the same wire pitch, sampling rate, etc
        raw::ChannelID_t channel = fHits[0]->Channel();
        // get the scale factor to convert dTick/dWire to dX/dU. This is used
        // to make the kink and merging cuts
        float wirePitch = geom->WirePitch(geom->View(channel));
        float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
        tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
        fScaleF = tickToDist / wirePitch;
//        HitSanityCheck();
        if(fMode < 2) {
          fStepDir = fMode;
          // make all trajectories with this setting
          ReconstructAllTraj();
        } else {
          // TEMP
          if(fPlane != 1) return;
          // here is where we get fancy
          // Step in the + direction
          fStepDir = 1;
          ReconstructAllTraj();
          // store the results
          trial.push_back(allTraj);
          inTrialTraj.push_back(inTraj);
          // THIS IS WRONG. IT WILL CLOBBER TRAJ'S IN OTHER PLANES
//          allTraj.clear();
          for(auto& intj : inTraj) intj = 0;
          // Step in the - direction
          fStepDir = -1;
          ReconstructAllTraj();
          // store the results
          trial.push_back(allTraj);
          inTrialTraj.push_back(inTraj);
//          allTraj.clear();
          AnalyzeTrials();
//          AdjudicateTrials(reAnalyze);
        }
      } // fPlane
    } // tpcid
    
    //    bool reAnalyze = false;

 
    // Convert trajectories in allTraj into clusters
    MakeAllTrajClusters();
//    CheckHitClusterAssociations();

    
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
          PrintTrajectory(allTraj[itj], USHRT_MAX);
          ang = std::abs(work.Pts[0].Ang);
          if(ang > M_PI/2) ang = M_PI - ang;
          // find average number of hits / TP
          PutTrajHitsInVector(allTraj[itj], hitVec);
          ave = (float)hitVec.size() / (float)allTraj[itj].Pts.size();
          fnHitsPerTP_Angle[ipl]->Fill(ang, ave);
          fnHitsPerTP_AngleP[ipl]->Fill(ang, ave);
        } // itj
      } // ipl
/*
      unsigned short lastPt = work.Pts.size() - 1;
      float ang = std::abs(work.Pts[lastPt].Ang);
      if(ang > M_PI/2) ang = M_PI - ang;
      mf::LogVerbatim("TC")<<"ANG "<<fPlane<<" "<<ang<<" "<<work.Pts[lastPt].AveHitRes<<" "<<work.Pts[lastPt].TP3Chi;
*/
    } // studymode

    
  } // RunStepCrawl
  
  //////////////////////////////////
  void TrajClusterAlg::FindDefaults()
  {
    // find default values of some tracking variables
    // TODO: This may not work well for events with few hits
    float avechg = 0;
    unsigned int wire, fhit, lhit, iht;
    float fcnt = 0;
    for(wire = fFirstWire; wire < fLastWire; ++wire) {
      fhit = WireHitRange[wire].first;
      lhit = WireHitRange[wire].second;
      for(iht = fhit; iht < lhit; ++iht) {
        if(fHits[iht]->Multiplicity() > 1) continue;
        if(WireHitRange[wire].first < 0) continue;
        if(fHits[iht]->GoodnessOfFit() < 0) continue;
        if(fHits[iht]->GoodnessOfFit() > 100) continue;
        avechg += fHits[iht]->Integral();
//        mf::LogVerbatim("TC")<<"hit "<<fHits[iht]->Integral()<<" "<<fHits[iht]->RMS()<<" "<<fHits[iht]->GoodnessOfFit();
        ++fcnt;
        if(fcnt == 50) break;
      } // iht
      if(fcnt == 50) break;
    } // wire
    
    avechg /= fcnt;
    fDefaultHitChgRMS = 0.3 * avechg;
    
  } // FindDefaults

  //////////////////////////////////
  void TrajClusterAlg::GetHitRange()
  {
    // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
    // Hits must have been sorted by increasing wire number
    fFirstHit = 0;
    geo::PlaneID planeID = DecodeCTP(fCTP);
    fNumWires = geom->Nwires(planeID.Plane, planeID.TPC, planeID.Cryostat);
    fMaxTime = (float)detprop->NumberTimeSamples();
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
    lariov::IChannelStatusProvider const& channelStatus
    = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();
    flag.first = -1; flag.second = -1;
    for(wire = 0; wire < fNumWires; ++wire) {
      raw::ChannelID_t chan = geom->PlaneWireToChannel
      ((int)planeID.Plane,(int)wire,(int)planeID.TPC,(int)planeID.Cryostat);
      if(channelStatus.IsBad(chan)) WireHitRange[wire] = flag;
    }
    
    // define the MergeAvailable vector and check for errors
//    if(mergeAvailable.size() < fHits.size()) throw art::Exception(art::errors::LogicError)
//      <<"GetHitRange: Invalid mergeAvailable vector size "<<mergeAvailable.size()<<fHits.size();
    unsigned int firstHit, lastHit;
    unsigned int cnt;
    cnt = 0;
//    float maxRMS, chiSep, peakCut;
    for(wire = 0; wire < fNumWires; ++wire) {
      // ignore dead wires and wires with no hits
      if(WireHitRange[wire].first < 0) continue;
      firstHit = WireHitRange[wire].first;
      lastHit = WireHitRange[wire].second;
      for(iht = firstHit; iht < lastHit; ++iht) {
        if(fHits[iht]->WireID().Wire != wire)
          throw art::Exception(art::errors::LogicError)<<"Bad WireHitRange on wire "<<wire<<"\n";
        ++cnt;
      } // iht
    } // wire
    if(cnt != nHitInPlane) mf::LogWarning("TC")<<"Bad WireHitRange count "<<cnt<<" "<<nHitInPlane<<"\n";
    
  } // GetHitRange()
  
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

    unsigned short itr, itj, jtr, jtj, nSameHits = 0;
    unsigned short nHitsInTraj;
    for(itr = 0; itr < trial.size(); ++itr) {
      nHitsInTraj = 0;
      for(itj = 0; itj < trial[itr].size(); ++itj) nHitsInTraj += trial[itr][itj].NHits;
      std::cout<<"Trial "<<itr<<" NHits in all trajectories "<<nHitsInTraj<<"\n";
    } // is

    std::vector<unsigned int> iHitVec, jHitVec;
    tjphs.clear();
    TjPairHitShare tmp;
    for(itr = 0; itr < trial.size() - 1; ++itr) {
      for(itj = 0; itj < trial[itr].size(); ++itj) {
        // ignore obsolete trajectories
        if(trial[itr][itj].ProcCode == USHRT_MAX) continue;
        // look at long trajectories for testing
        if(trial[itr][itj].Pts.size() < 5) continue;
        PutTrajHitsInVector(trial[itr][itj], iHitVec);
//        std::cout<<"itr "<<itr<<" itj "<<itj<<" hit size "<<iHitVec.size()<<" nSameHits ";
        for(jtr = itr + 1; jtr < trial.size(); ++jtr) {
          for(jtj = 0; jtj < trial[jtr].size(); ++jtj) {
            if(trial[jtr][jtj].ProcCode == USHRT_MAX) continue;
            PutTrajHitsInVector(trial[jtr][jtj], jHitVec);
            CountSameHits(iHitVec, jHitVec, nSameHits);
//            std::cout<<" "<<nSameHits;
            if(nSameHits == 0) continue;
            tmp.iTrial = itr; tmp.iTj = itj;
            tmp.jTrial = jtr; tmp.jTj = jtj;
            tmp.nSameHits = nSameHits;
            tjphs.push_back(tmp);
            // keep track of the
          } // jtj
        } // jtr
//        std::cout<<"\n";
      } // itj
    } // itr
    
    std::cout<<"  itr  itj iNHits  jtr  jtj jNHits nSameHits\n";
    for(auto& tmp : tjphs) {
      itj = tmp.iTj; jtj = tmp.jTj;
      std::cout<<std::setw(5)<<tmp.iTrial<<std::setw(5)<<itj<<std::setw(7)<<trial[tmp.iTrial][itj].NHits;
      std::cout<<std::setw(5)<<tmp.jTrial<<std::setw(5)<<jtj<<std::setw(7)<<trial[tmp.jTrial][jtj].NHits;
      std::cout<<std::setw(6)<<tmp.nSameHits<<"\n";
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
    nSameHits = 0;
    for(unsigned short ii = 0; ii < iHitVec.size(); ++ii) {
      if(std::find(jHitVec.begin(), jHitVec.end(), iHitVec[ii]) != jHitVec.end()) ++nSameHits;
    }
  } // CountSameHits

  ////////////////////////////////////////////////
  void TrajClusterAlg::PutTrajHitsInVector(Trajectory const& tj, std::vector<unsigned int>& hitVec)
  {
    hitVec.resize(tj.NHits);
    unsigned short ipt, iht, cnt;
    cnt = 0;
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      for(iht = 0; iht < tj.Pts[ipt].Hits.size(); ++iht) {
        if(cnt > tj.NHits-1) {
          mf::LogError("TC")<<"PutTrajHitsInVector: traj NHits not sized properly "<<tj.NHits<<" cnt "<<cnt;
          hitVec.clear();
          return;
        }
        hitVec[cnt] = tj.Pts[ipt].Hits[iht];
        ++cnt;
      }
    }
  } // PutTrajHitsInVector
   
  ////////////////////////////////////////////////
  void TrajClusterAlg::ReconstructAllTraj()
  {
    // Reconstruct clusters in fPlane and put them in allTraj
    
    unsigned int ii, iwire, jwire, iht, jht, oht;
    
    unsigned int nwires = fLastWire - fFirstWire - 1;
    unsigned int ifirsthit, ilasthit, jfirsthit, jlasthit;
    float fromWire, fromTick, toWire, toTick, deltaRms, qtot;
    bool success, SignalPresent, updateOK;
    bool didPrt = false;

    for(fPass = 0; fPass < fMinNPtsFit.size(); ++fPass) {
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
          prt = (fDebugPlane == (int)fPlane && (int)iwire == fDebugWire && std::abs((int)fHits[iht]->PeakTime() - fDebugHit) < 10);
          if(prt) mf::LogVerbatim("TC")<<"+++++++ Pass "<<fPass<<" Found debug hit "<<fPlane<<":"<<fHits[iht]->WireID().Wire<<":"<<(int)fHits[iht]->PeakTime()<<" inTraj "<<inTraj[iht]<<" RMS "<<fHits[iht]->RMS()<<" Multiplicity "<<fHits[iht]->Multiplicity()<<" LocalIndex "<<fHits[iht]->LocalIndex();
          if(prt) didPrt = true;
          if(inTraj[iht] != 0) continue;
          fromWire = fHits[iht]->WireID().Wire;
          fromTick = fHits[iht]->PeakTime();
          // decide whether to "merge" two hits in a multiplet
          if(fHits[iht]->Multiplicity() == 2) {
            // get the index of the other hit in the multiplet
            oht = 1 - fHits[iht]->LocalIndex();
            // Make sure it isn't used and check the separation
            if(inTraj[oht] == 0 && std::abs(fHits[oht]->PeakTime() - fromTick) < fMultHitSep * fHits[oht]->RMS()) {
              HitMultipletPosition(iht, fromTick, deltaRms, qtot);
              if(prt) mf::LogVerbatim("TC")<<" Hit doublet: back from HitMultipletPosition ";
            }
          }
          for(jht = jfirsthit; jht < jlasthit; ++jht) {
            if(inTraj[iht] != 0) continue;
            if(inTraj[jht] == -3) mf::LogVerbatim("TC")<<"ReconstructAllTraj bad jht flag "<<fPlane<<":"<<fHits[jht]->WireID().Wire<<":"<<(int)fHits[jht]->PeakTime()<<" using iht "<<fPlane<<":"<<fHits[iht]->WireID().Wire<<":"<<(int)fHits[iht]->PeakTime();
            if(inTraj[jht] != 0) continue;
            if(prt) mf::LogVerbatim("TC")<<"+++++++ checking ClusterHitsOK with jht "<<fPlane<<":"<<fHits[jht]->WireID().Wire<<":"<<(int)fHits[jht]->PeakTime()<<" RMS "<<fHits[jht]->RMS()<<" Multiplicity "<<fHits[jht]->Multiplicity()<<" LocalIndex "<<fHits[jht]->LocalIndex();
            // Ensure that the hits StartTick and EndTick have the proper overlap
            if(!TrajHitsOK(iht, jht)) continue;
/*
            // Ensure that we pick up all the hits in a iht multiplet.
            // We can do this by skipping this iht if we would be stepping away from
            // the other hits.
            if(fHits[iht]->Multiplicity() > 1) {
              // Will move in -time direction, away from the hits in the +time direction from iht
              if(fHits[iht]->Multiplicity() > 2 && fHits[iht]->LocalIndex() < (fHits[iht]->Multiplicity()-1) && fHits[iht]->PeakTime() > fHits[jht]->PeakTime()) continue;
              // The reverse case **shouldn't** be needed because the hit loop iterates from -time to +time
              // Also check to ensure that we don't start going in a non-physical direction when
              // the other hits in the multiplet are used. Do the simple case of a doublet
              if(fHits[iht]->Multiplicity() == 2) {
                if(fHits[iht]->LocalIndex() == 0 && inTraj[iht+1] > 0 && fHits[jht]->PeakTime() > fHits[iht]->PeakTime()) continue;
                if(fHits[iht]->LocalIndex() == 1 && inTraj[iht-1] > 0 && fHits[jht]->PeakTime() < fHits[iht]->PeakTime()) continue;
              }
              // Large multiplicity hits TODO deal with this when GausHit is fixed
            } // fHits[iht]->Multiplicity() > 1
*/
            if(prt) mf::LogVerbatim("TC")<<"  Starting trajectory";
            // start a trajectory in the direction from iht -> jht
            toWire = jwire;
            if(fHits[jht]->Multiplicity() == 1) {
              toTick = fHits[jht]->PeakTime();
            } else {
              HitMultipletPosition(jht, toTick, deltaRms, qtot);
            }
            StartTraj(fromWire, fromTick, toWire, toTick);
            // check for a failure
            if(work.Pts.size() == 0) {
              if(prt) mf::LogVerbatim("TC")<<"ReconstructAllTraj: StartTraj failed";
              continue;
            }
            // Adjust the TP delta RMS
            if(fHits[iht]->Multiplicity() == 2) work.Pts[0].DeltaRMS = deltaRms;
            // try to add close hits
            AddTrajHits(work.Pts[0], SignalPresent);
            if(!SignalPresent || work.Pts[0].Hits.size() == 0) {
              if(prt) mf::LogVerbatim("TC")<<" No hits at initial trajectory point ";
              FlagWorkHits(0);
              continue;
            }
            UpdateWork(updateOK);
            if(prt) PrintTrajectory(work, USHRT_MAX);
            // now try stepping away
            StepCrawl(success);
            FlagWorkHits(0);
            if(!success) {
              if(prt) mf::LogVerbatim("TC")<<" xxxxxxx StepCrawl failed ";
              continue;
            }
            if(prt) mf::LogVerbatim("TC")<<"StepCrawl done: work.Pts size "<<work.Pts.size();
            if(work.Pts.size() < fMinPts) continue;
            StoreWork();
            break;
          } // jht
        } // iht
      } // iwire
    } // fPass
    
    
    prt = false;

    FillTrajTruth();
    
    if(didPrt) {
      mf::LogVerbatim("TC")<<"Done in ReconstructAllTraj";
      PrintAllTraj(USHRT_MAX, 0);
    }
    
    work.Pts.clear();
    
    // need to fake out FindVertices
//    if(fFindTrajVertices) FindVertices();
    if(fTagAllTraj) TagAllTraj();
    
    
  } // ReconstructAllTraj

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
      if(allTraj[itj].ProcCode == USHRT_MAX) continue;
      sIDEs.clear();
      for(ii = 0; ii < tidcnt.size(); ++ii) tidcnt[ii] = 0;
      for(ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
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
      maxCnt = allTraj[itj].Pts.size() / 2;
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
      if(tj1.ProcCode == USHRT_MAX) continue;
      if(tj1.CTP != fPlane) continue;
      for(i2 = i1 + 1; i2 < allTraj.size(); ++i2) {
        Trajectory& tj2 = allTraj[i2];
        if(tj2.ProcCode == USHRT_MAX) continue;
        if(tj2.CTP != tj1.CTP) continue;
        // find the closest approach
        minSep2 = maxSep2;
        TrajTrajDOCA(tj1, tj2, ipt1, ipt2, minSep2);
        if(minSep2 == maxSep2) continue;
        // Count the number at each end and in the middle
        bin1 = (unsigned short)(3 * (float)ipt1 / (float)tj1.Pts.size());
        if(bin1 > 2) bin1 = 2;
        // only count if this end doesn't have a vertex
        if(tj1.Vtx[bin1] < 0) ++nCloseEnd[i1][bin1];
        bin2 = (unsigned short)(3 * (float)ipt2 / (float)tj2.Pts.size());
        if(bin2 > 2) bin2 = 2;
        if(tj2.Vtx[bin2] < 0) ++nCloseEnd[i2][bin2];
        // find the angle between the TPs at the intersection
        dang = std::abs(tj1.Pts[ipt1].Ang - tj2.Pts[ipt2].Ang);
        if(bin1 != 1  && bin2 != 1) {
          // the DOCA point is at the ends of the two TJs.
          // Find the intersection using the appropriate end points
          endPt1 = 0;
          if(bin1 == 2) endPt1 = allTraj[i1].Pts.size() - 1;
          endPt2 = 0;
          if(bin2 == 2) endPt2 = allTraj[i2].Pts.size() - 1;
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
/* This creates too many vertices. Better to do this after a shower is identified
        // make a vertex instead if the vertex position and the DOCA position
        // are close && the angle difference is sufficiently large and the
        // DOCA positions are not in the middle of either trajectory
        if(dvtx2 < 5 && dangSig > 3 && bin1 != 1 && bin2 != 1) {
          MakeTrajVertex(aTrjInt, bin1, bin2, madeVtx);
          if(madeVtx) continue;;
        }
*/
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
  void TrajClusterAlg::MakeTrajVertex(TrjInt const& aTrjInt, unsigned short bin1, unsigned short bin2, bool& madeVtx)
  {
    // Make 2D vertex
    madeVtx = false;
    
    // bin1 and bin2 are either 0, 1, or 2 which denote whether the
    // DOCA point is in the first third, second third or last third of
    // the trajectory points. Convert this into a index for referencing
    // the vertex end for trajectory 1 and trajectory 2
    unsigned short end1 = bin1 / 2;
    unsigned short end2 = bin2 / 2;
   
    unsigned short itj1 = aTrjInt.itj1;
    unsigned short itj2 = aTrjInt.itj2;
    
    VtxStore aVtx;
    aVtx.Wire = aTrjInt.vw;
    aVtx.WireErr = 1;
    aVtx.Time = aTrjInt.vt / fScaleF;
    aVtx.TimeErr = 1;
    aVtx.NClusters = 2;
    aVtx.ChiDOF = 0;
    aVtx.CTP = fCTP;
    aVtx.Topo = 9;
    vtx.push_back(aVtx);
    unsigned short ivx = vtx.size() - 1;
    // make the vertex - traj association
    allTraj[itj1].Vtx[end1] = ivx;
    allTraj[itj2].Vtx[end2] = ivx;
    
    // try to attach other TJs to it
    unsigned short end, ipt, oEndPt;
    float dw, dt, dr;
    for(unsigned short itj = 0; itj < allTraj.size(); ++itj) {
      if(allTraj[itj].ProcCode == USHRT_MAX) continue;
      if(allTraj[itj].CTP != fCTP) continue;
      for(end = 0; end < 2; ++end) {
        if(allTraj[itj].Vtx[end] >= 0) continue;
        if(end == 0) { ipt = 0; } else { ipt = allTraj[itj].Pts.size() - 1; }
        dw = aTrjInt.vw - allTraj[itj].Pts[ipt].Pos[0];
        // make some rough cuts
        if(std::abs(dw) > 3) continue;
        dt = aTrjInt.vt - allTraj[itj].Pts[ipt].Pos[1];
        if(std::abs(dt) > 3) continue;
        // distance between this end and the vertex
        dr = dw * dw + dt * dt;
        // TODO convert this to a chisq cut
        if(dr > 5) continue;
        // make sure that the other end isn't closer - short trajectories
        if(allTraj[itj].Pts.size() < 5) {
          if(end == 0) { oEndPt = allTraj[itj].Pts.size() - 1; } else { oEndPt = 0; }
          dw = aTrjInt.vw - allTraj[itj].Pts[oEndPt].Pos[0];
          dt = aTrjInt.vw - allTraj[itj].Pts[oEndPt].Pos[1];
          if((dw * dw + dt * dt) < dr) continue;
        } // short trajectory
        // attach it
        allTraj[itj].Vtx[end] = ivx;
        ++vtx[ivx].NClusters;
      } // end
    } // itj
    
    mf::LogVerbatim("TC")<<" Made vertex "<<ivx<<" with "<<vtx[ivx].NClusters;
   
    madeVtx = true;

    
  } // MakeTrajVertex

  //////////////////////////////////////////
  void TrajClusterAlg::DefineShowerTraj(unsigned short icot, std::vector<std::vector<unsigned short>> trjintIndices)
  {
    // This routine is called if there area significant (> 3) number
    // of trjints enabling an estimate of the shower direction, etc.
    
    unsigned short ii, iTji, i1, ipt1, i2, ipt2, nEndInside;
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
      arg = allTraj[i1].Pts.size() + allTraj[i2].Pts.size();
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
      ipt1 = allTraj[i1].Pts.size()-1;
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
      arg = 10 * (1 + allTraj[i1].Pass) * dang * allTraj[i1].AveTP3Chi * allTraj[i1].Pts.size();
//      mf::LogVerbatim("TC")<<"Candidate "<<i1<<" nCloseEnd "<<nCloseEnd[i1][0]<<" "<<nCloseEnd[i1][1]<<" "<<nCloseEnd[i1][2]<<" pass "<<allTraj[i1].Pass<<" dang "<<dang<<" arg "<<arg;
      if(arg < fom) {
        fom = arg;
        primTraj = i1;
      }
    } // ii
    // determine which end is closest to the shower end
    // distance between TP0 and the shower min position
    ipt1 = allTraj[primTraj].Pts.size()-1;
    pos0 = cs * allTraj[primTraj].Pts[0].Pos[0] - sn * allTraj[primTraj].Pts[0].Pos[1];
    pos1 = cs * allTraj[primTraj].Pts[ipt1].Pos[0] - sn * allTraj[primTraj].Pts[ipt1].Pos[1];
    minsep = std::abs(pos0 - minshlong);
    primTrajEnd = 0;
    //        mf::LogVerbatim("TC")<<"0min "<<pos0<<" minsep "<<minsep<<" end "<<primTrajEnd;
    arg = std::abs(pos1 - minshlong);
    if(arg < minsep) { minsep = arg;  primTrajEnd = ipt1; }
    //        mf::LogVerbatim("TC")<<"1min "<<pos1<<" arg "<<arg<<" end "<<primTrajEnd;
    arg = std::abs(pos0 - maxshlong);
    if(arg < minsep) { minsep = arg;  primTrajEnd = 0; }
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
  
  //////////////////////////////////////////
  float TrajClusterAlg::TrajLength(Trajectory& tj)
  {
    float len = 0, dx, dy;
    for(unsigned short ipt = 1; ipt < tj.Pts.size(); ++ipt) {
      dx = tj.Pts[ipt].Pos[0] - tj.Pts[ipt-1].Pos[0];
      dy = tj.Pts[ipt].Pos[1] - tj.Pts[ipt-1].Pos[1];
      len += sqrt(dx * dx + dy * dy);
    }
    return len;
  } // TrajLength
  
  //////////////////////////////////////////
  float TrajClusterAlg::TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
  {
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation
/*
  //////////////////////////////////////////
  void TrajClusterAlg::FindTrajVertices()
  {
    // Looks for vertices that ClusterCrawler can't find.
    // All trajectories **should_be** in ClusterCrawler order

    if(allTraj.size() < 2) return;

    vtxprt = (fDebugPlane == (int)fPlane && fDebugHit < 0);

    VtxStore aVtx;

    unsigned short itj, iend, itjPt, jtj, jend, jtjPt, ivx, clsIndex;
    unsigned short ktj, kend;
    float xint, yint, dwi, dti, dwj, dtj;
//    float dr2cut = 20;
    bool madeVtx;
    // first look for vertices at the ends of trajectories
    for(itj = 0; itj < allTraj.size() - 1; ++itj) {
      for(iend = 0; iend < 2; ++iend) {
        if(allTraj[itj].Vtx[iend] >= 0) continue;
        if(iend == 0) { itjPt = 0; } else { itjPt = allTraj[itj].Pts.size()-1; };
        madeVtx = false;
        for(jtj = itj + 1; jtj < allTraj.size(); ++jtj) {
          for(jend = 0; jend < 2; ++jend) {
            if(allTraj[jtj].Vtx[iend] >= 0) continue;
            if(jend == 0) { jtjPt = 0; } else { jtjPt = allTraj[jtj].Pts.size()-1; };
            TrajIntersection(allTraj[itj].Pts[itjPt], allTraj[jtj].Pts[jtjPt], xint, yint);
            dwi = xint - allTraj[itj].Pts[itjPt].Pos[0];
            if(std::abs(dwi) > 5) continue;
            dwj = xint - allTraj[jtj].Pts[jtjPt].Pos[0];
            if(std::abs(dwj) > 5) continue;
            dti = std::abs(allTraj[itj].Pts[itjPt].Pos[1] - yint);
            if(dti > 5) continue;
            dtj = std::abs(allTraj[jtj].Pts[jtjPt].Pos[1] - yint);
            if(dtj > 5) continue;
            aVtx.Wire = xint;
            aVtx.WireErr = 1;
            aVtx.Time = yint / fScaleF;
            aVtx.TimeErr = 1;
            aVtx.NClusters = 2; aVtx.CTP = fCTP; aVtx.Fixed = false;
            aVtx.ChiDOF = 0;    aVtx.Topo = 9;
            vtx.push_back(aVtx);
            ivx = vtx.size() - 1;
            // make the vertex - traj association
            allTraj[itj].Vtx[iend] = ivx;
            allTraj[jtj].Vtx[jend] = ivx;
            // make the vertex - cluster association
            clsIndex = allTraj[itj].ClusterIndex;
            if(VtxClusterEnd(ivx, clsIndex) == 0) { tcl[clsIndex].BeginVtx = ivx; } else {tcl[clsIndex].EndVtx = ivx;}
            clsIndex = allTraj[jtj].ClusterIndex;
            if(VtxClusterEnd(ivx, clsIndex) == 0) { tcl[clsIndex].BeginVtx = ivx; } else {tcl[clsIndex].EndVtx = ivx;}
            if(vtxprt) {
              mf::LogVerbatim("TC")<<"Vtx itj "<<itj<<" clsID "<<allTraj[itj].ClusterIndex+1<<" itjPt "<<itjPt<<" jtj "<<jtj<<" clsID "<<allTraj[jtj].ClusterIndex+1<<" jtjPt "<<jtjPt<<" Vtx "<<aVtx.Wire<<" "<<aVtx.Time;
              PrintAllTraj(itj, itjPt);
              PrintAllTraj(jtj, jtjPt);
            }
            madeVtx = true;
            // Look for another
            for(ktj = jtj + 1; ktj < allTraj.size(); ++ktj) {
              for(kend = 0; kend < 2; ++kend) {
                if(allTraj[ktj].Vtx[kend] >= 0) continue;
                clsIndex = allTraj[ktj].ClusterIndex;
                if(ClusterVertexChi(clsIndex, kend, ivx) < fVertex2DCut) {
                  allTraj[ktj].Vtx[kend] = ivx;
                  if(VtxClusterEnd(ivx, clsIndex) == 0) { tcl[clsIndex].BeginVtx = ivx; } else { tcl[clsIndex].EndVtx = ivx; }
                  FitVtx(ivx);
                  if(vtxprt) mf::LogVerbatim("TC")<<" add ktj"<<ktj<<" clsID "<<clsIndex+1<<" Vtx "<<vtx[ivx].Wire<<" "<<vtx[ivx].Time;
                  break;
                } // good chisq
              } // kend
            } // ktj
          } // jend
          if(madeVtx) break;
        } // jtj
      } // iend
    } // itj

  } // FindTrajVertices
*/

  //////////////////////////////////////////
  void TrajClusterAlg::TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& ipt, float& Distance)
  {
    // find the closest approach between a trajectory tj and a point (x,y). Returns
    // the index of the closest trajectory point and the distance
    
    float dx, dy, dist, best = 1E6;
    ipt = 0;
    Distance = best;
    
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      dx = tj.Pts[ii].Pos[0] - x;
      dy = tj.Pts[ii].Pos[1] - y;
      dist = dx * dx + dy * dy;
      if(dist < best) {
        best = dist;
        ipt = ii;
      }
      // TODO is this wise?
//      if(dist > best) break;
    } // ii
    
    Distance = sqrt(best);
    
  } // TrajClosestApproach

  
  ////////////////////////////////////////////////
  void TrajClusterAlg::AddTrajHits(TrajPoint& tp, bool& SignalPresent)
  {
    std::vector<unsigned int> closeHits;
    std::vector<float> closeDeltas;
    unsigned int wire, loWire, hiWire, iht, firstHit, lastHit;
    unsigned int lastPt = work.Pts.size() - 1;
    
    // figure out which wires to consider
    // On the first entry only consider the wire the TP is on
    if(work.Pts.size() == 1) {
      loWire = std::nearbyint(tp.Pos[0]);
      hiWire = loWire + 1;
    } else  {
      if(IsLargeAngle(tp)) {
        // look at adjacent wires for larger angle trajectories
        loWire = std::nearbyint(tp.Pos[0] - 1);
        hiWire = loWire + 3;
//        hiWire = std::nearbyint(tp.Pos[0] + 2);
      } else {
        loWire = std::nearbyint(tp.Pos[0]);
        hiWire = loWire + 1;
      }
    } // work.Pts.size > 1
    
    float fwire, ftime, delta;
    
    // find the projection error to this point
    float dw = tp.Pos[0] - work.Pts[lastPt].Pos[0];
    float dt = tp.Pos[1] - work.Pts[lastPt].Pos[1];
    float dpos = sqrt(dw * dw + dt * dt);
    float projErr = dpos * work.Pts[lastPt].AngErr;
    // Add this to the Delta RMS factor and construct a cut
    float deltaCut = 3 * (projErr + tp.DeltaRMS);
    if(IsLargeAngle(tp)) {
      if(deltaCut < 1.1) deltaCut = 1.1;
    } else {
      if(deltaCut < 0.1) deltaCut = 0.1;
    }
    deltaCut *= fProjectionErrFactor;
    // cut for nearby hits
    //    float bigDeltaCut = 3 * deltaCut;
    // make it big for study mode
    if(fStudyMode) deltaCut = 1.5;
    
    // projected time in ticks for testing the existence of a hit signal
    raw::TDCtick_t rawProjTick = (float)(tp.Pos[1] / fScaleF);
    
    // assume failure
    SignalPresent = false;
    if(prt) mf::LogVerbatim("TC")<<" AddTrajHits: loWire "<<loWire<<" tp.Pos[0] "<<tp.Pos[0]<<" hiWire "<<hiWire<<" projTick "<<rawProjTick<<" fHitChgRMS "<<fHitChgRMS<<" deltaRMS "<<tp.DeltaRMS<<" projErr "<<projErr<<" tp.Dir[0] "<<tp.Dir[0]<<" IsLA "<<IsLargeAngle(tp);
    
    for(wire = loWire; wire < hiWire; ++wire) {
      // Assume a signal exists on a dead wire
      if(WireHitRange[wire].first == -1) SignalPresent = true;
      if(WireHitRange[wire].first < 0) continue;
      firstHit = (unsigned int)WireHitRange[wire].first;
      lastHit = (unsigned int)WireHitRange[wire].second;
      fwire = wire;
      for(iht = firstHit; iht < lastHit; ++iht) {
        if(rawProjTick > fHits[iht]->StartTick() && rawProjTick < fHits[iht]->EndTick()) SignalPresent = true;
        ftime = fScaleF * fHits[iht]->PeakTime();
        delta = PointTrajDOCA(fwire, ftime, tp);
        float dt = std::abs(ftime - tp.Pos[1]);
        if(prt) {
          bool close = delta < fChgDiffCut;
          mf::LogVerbatim myprt("TC");
          myprt<<"  chk "<<fHits[iht]->WireID().Plane<<":"<<fHits[iht]->WireID().Wire<<":"<<(int)fHits[iht]->PeakTime();
          myprt<<" delta "<<std::fixed<<std::setprecision(2)<<delta<<" deltaCut "<<deltaCut<<" dt "<<dt;
          myprt<<" Mult "<<fHits[iht]->Multiplicity()<<" localIndex "<<fHits[iht]->LocalIndex()<<" RMS "<<std::setprecision(1)<<fHits[iht]->RMS();
          myprt<<" Chi "<<std::setprecision(1)<<fHits[iht]->GoodnessOfFit();
          myprt<<" inTraj "<<inTraj[iht];
          myprt<<" Chg "<<(int)fHits[iht]->Integral();
          myprt<<" Signal? "<<SignalPresent<<" close? "<<close;
          if(inTraj[iht] == 2) PrintAllTraj(1,USHRT_MAX);
        }
        if(inTraj[iht] != 0) continue;
        if(delta > deltaCut) continue;
        // The impact parameter delta may be good but we may be projecting
        // the trajectory too far away (in time) from the current position
        if(dt > 3) continue;
        SignalPresent = true;
        inTraj[iht] = -3;
        closeHits.push_back(iht);
        closeDeltas.push_back(delta);
      } // iht
    } // wire
    if(!SignalPresent) {
      if(prt) mf::LogVerbatim("TC")<<" no signal on any wire at tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" tick "<<(int)tp.Pos[1]/fScaleF<<" closeHits size "<<closeHits.size();
      return;
    }
    if(closeHits.size() == 0) return;
    // decide whether to use some or all of the hits using a figure of merit.
    // This only works after we get enough information
    std::vector<unsigned int> useHits;
    unsigned int jj, jht, kk, kht;
    float bestFOM = 50;
    // Consider single hits first
    float ddif, qdif, fom, qtot, dRms;
    for(jj = 0; jj < closeHits.size(); ++jj) {
      jht = closeHits[jj];
      ddif = closeDeltas[jj] / tp.DeltaRMS;
      if(IsLargeAngle(tp)) {
        // nothing fancy for large angles
        fom = closeDeltas[jj];
        qdif = 0;
      } else {
        // de-weight the charge contribution for small angle TPs
        qdif = std::abs(fHits[jht]->Integral() - tp.AveChg) / (2 * fHitChgRMS);
        fom = (ddif + qdif) / 2;
      }
      if(fom < bestFOM) {
        bestFOM = fom;
        useHits.resize(1);
        useHits[0] = jht;
      }
      if(prt) mf::LogVerbatim("TC")<<"Singles "<<fHits[jht]->WireID().Plane<<":"<<fHits[jht]->WireID().Wire<<":"<<(int)fHits[jht]->PeakTime()<<" ddif "<<ddif<<" qdif "<<qdif<<" fom "<<fom<<" bestFOM "<<bestFOM<<" deltaRMS "<<tp.DeltaRMS;
    } // jj
    // Next consider pairs of hits whose charge may sum to a value
    // that is consistent with a trajectory hit. The delta of this
    // pair is calculated using a charge weighted average.
    for(jj = 0; jj < closeHits.size()-1; ++jj) {
      jht = closeHits[jj];
      for(kk = jj + 1; kk < closeHits.size(); ++kk) {
        kht = closeHits[kk];
        // check the hit separation
        if(LargeHitSep(jht, kht)) continue;
        // only consider hits on the same wire in the same multiplet
        if(fHits[jht]->WireID().Wire != fHits[kht]->WireID().Wire) continue;
        // TODO: this needs to be done better...
        if(jht > kht && jht != kht + 1) continue;
        if(kht > jht && kht != jht + 1) continue;
        HitMultipletPosition(jht, ftime, dRms, qtot);
        ftime *= fScaleF;
        fwire = fHits[jht]->WireID().Wire;
        ddif = PointTrajDOCA(fwire, ftime, tp) / dRms;
        // Charge difference, de-weighted
        qdif = std::abs(qtot- tp.AveChg) / (2 * fHitChgRMS);
        fom = (ddif + qdif) / 2;
        if(prt) mf::LogVerbatim("TC")<<"Pairs "<<fHits[jht]->WireID().Plane<<":"<<fHits[jht]->WireID().Wire<<":"<<(int)fHits[jht]->PeakTime()<<" + "<<fHits[kht]->WireID().Plane<<":"<<fHits[kht]->WireID().Wire<<":"<<(int)fHits[kht]->PeakTime()<<" ddif "<<ddif<<" qdif "<<qdif<<" fom "<<fom<<" bestFOM "<<bestFOM;
        if(fom < bestFOM) {
          bestFOM = fom;
          useHits.resize(2);
          useHits[0] = jht;
          useHits[1] = kht;
        }
      } // kk
    } // jj
    tp.Hits.insert(tp.Hits.end(), useHits.begin(), useHits.end());
    // Find the combined hit time error
    tp.HitsTimeErr2 = HitsTimeErr2(tp.Hits);
    tp.NCloseNotUsed = closeHits.size() - useHits.size();
    // release all hits flagged in closeHit
    for(jj = 0; jj < closeHits.size(); ++jj) inTraj[closeHits[jj]] = 0;
    // re-flag the hits in useHits
    for(jj = 0; jj < useHits.size(); ++jj) inTraj[useHits[jj]] = -3;
    // specify the hits position
    tp.HitPos[0] = 0;
    tp.HitPos[1] = 0;
    tp.Chg = 0;
    if(tp.Hits.size() > 1) {
      // sort the new hits by distance from the previous traj point
      std::vector<SortEntry> sortVec;
      SortEntry sortEntry;
      unsigned short ii;
      float dw, dt;
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        sortEntry.index = ii;
        iht = tp.Hits[ii];
        dw = fHits[iht]->WireID().Wire - work.Pts[lastPt].Pos[0];
        dt = fHits[iht]->PeakTime() * fScaleF - work.Pts[lastPt].Pos[1];
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
    // Calculate the charge weighted position of the hits associated with this traj point
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      iht = tp.Hits[ii];
      tp.Chg += fHits[iht]->Integral();
      tp.HitPos[0] += fHits[iht]->Integral() * fHits[iht]->WireID().Wire;
      // convert to WSE units
      tp.HitPos[1] += fHits[iht]->Integral() * fHits[iht]->PeakTime() * fScaleF;
    }
    tp.HitPos[0] /= tp.Chg;
    tp.HitPos[1] /= tp.Chg;
    //    if(prt) mf::LogVerbatim("TC")<<" HitPos "<<tp.HitPos[0]<<" tick "<<(int)tp.HitPos[1]/fScaleF<<" nhits "<<closeHits.size();
    
  } // AddTrajHits
  
  //////////////////////////////////////////
  float TrajClusterAlg::HitsTimeErr2(std::vector<unsigned int> const& hitVec)
  {
    // Estimates the error^2 in the time component using all hits in hitVec
    
    if(hitVec.size() == 0) return 0;
    // single hit
    float err = fHits[hitVec[0]]->RMS() * fScaleF * fHitErrFac;
    // testing
//    float err = fHits[hitVec[0]]->SigmaPeakTime() * fScaleF;
    if(hitVec.size() == 1) return err * err;
    
    // This approximation works for two hits of roughly similar RMS
    // and with roughly similar (within 2X) amplitude. TODO deal with the
    // case of more than 2 hits if the need arises
    float averms = 0.5 * (fHits[hitVec[0]]->RMS() + fHits[hitVec[1]]->RMS());
    float hitsep = (fHits[hitVec[0]]->PeakTime() - fHits[hitVec[1]]->PeakTime()) / averms;
    // This will estimate the RMS of two hits separated by 1 sigma by 1.5 * RMS of one hit
    err = averms * (1 + 0.14 * hitsep * hitsep) * fScaleF * fHitErrFac;
    return err * err;
    
  } // HitsTimeErr2

  /*
  //////////////////////////////////////////
  void TrajClusterAlg::FindTrajVertices()
  {
    // Looks for vertices that ClusterCrawler can't find.
    // All trajectories **should_be** in ClusterCrawler order
    
    if(allTraj.size() < 2) return;
    
    vtxprt = (fDebugPlane == (int)plane && fDebugHit < 0);
    vtxprt = true;
    
    struct IntPoint {
      unsigned short iTraj;
      unsigned short iTrajPt;
      unsigned short jTraj;
      unsigned short jTrajPt;
      float DOCA;
    };
    std::vector<IntPoint> intPts;

    unsigned short itj, iend, itjEnd, jtj, jtjEnd, jtjPt, jint, ii;
    float best, dw, dt, dr;
    bool sameij;
    
    // first look for vertices at the ends of trajectories
    for(itj = 0; itj < allTraj.size() - 1; ++itj) {
      for(iend = 0; iend < 2; ++iend) {
        if()
        for(jtj = itj + 1; jtj < allTraj.size(); ++jtj) {
          for(jend = 0; jend < 2; ++jend) {
            
          } // jend
        } // jtj
      } // iend
    } // itj

    for(itj = 0; itj < allTraj.size(); ++itj) {
      if(allTraj[itj].ProcCode == USHRT_MAX) continue;
      for(iend = 0; iend < 2; ++iend) {
        if(iend == 0) { itjEnd = 0; } else { itjEnd = allTraj[itj].Pts.size() - 1; }
        for(jtj = 0; jtj < allTraj.size(); ++jtj) {
          if(allTraj[jtj].ProcCode == USHRT_MAX) continue;
          if(itj == jtj) continue;
          jtjEnd = allTraj[jtj].Pts.size() - 1;
          // check that the end position (wire) of the i trajectory that we are considering is within the range
          // of positions (wires) that is spanned by the j trajectory.
          if(allTraj[itj].Pts[itjEnd].Pos[0] > allTraj[jtj].Pts[0].Pos[0]) continue;
          if(allTraj[itj].Pts[itjEnd].Pos[0] < allTraj[jtj].Pts[jtjEnd].Pos[1]) continue;
          // check DOCA between the end of itj and the trajectory points in jtj
          // best = 10 is equivalent to a DOCA of sqrt(10) in WSE
          best = 10;
          jint = USHRT_MAX;
          for(jtjPt = 0; jtjPt <= jtjEnd; ++jtjPt) {
            dw = allTraj[jtj].Pts[jtjPt].Pos[0] - allTraj[itj].Pts[itjEnd].Pos[0];
            dt = allTraj[jtj].Pts[jtjPt].Pos[1] - allTraj[itj].Pts[itjEnd].Pos[1];
            dr = dw * dw + dt * dt;
            if(dr < best) {
              best = dr;
              jint = jtjPt;
            }
          } // jtjPt
          if(jint == USHRT_MAX) continue;
          // We have an intersection point
          IntPoint ip;
          ip.iTraj = itj; ip.iTrajEnd = end; ip.jTraj = jtj; ip.jTrajPt = jint; ip.DOCA = best;
          // see if there exists another intersection with the same two trajectories but at the other end.
          sameij = false;
          for(ii = 0; ii < intPts.size(); ++ii) {
            if(itj == intPts[ii].iTraj && jtj == intPts[ii].jTraj && best < intPts[ii].DOCA) {
              // This one is better. Overwrite it
              intPts[ii] = ip;
              sameij = true;
              break;
            }
          } // ii
          // Add a new intersection point
          if(!sameij) intPts.push_back(ip);
        } // jtj
      } // end
    } // itj

    if(vtxprt) mf::LogVerbatim("TC")<<"FindTrajVertices: intPts size "<<intPts.size();
    
    if(intPts.size() == 0) return;
    
    if(vtxprt) {
      for(ii = 0; ii < intPts.size(); ++ii) {
        itj = intPts[ii].iTraj;
        jtj = intPts[ii].iTraj;
        mf::LogVerbatim("TC")<<ii<<" itj "<<itj<<" cls "<<allTraj[itj].ClusterIndex<<" itjPt "<<intPts[ii].iTrajPt<<" jtj "<<jtj<<" cls "<<allTraj[jtj].ClusterIndex<<" jPt "<<intPts[ii].jTrajPt<<" DOCA "<<intPts[ii].DOCA;
        
      }
    } // vtxprt

    // We have intersection points to turn into 2D vertices. Decide whether to create
    // a vertex with one cluster (e.g. delta ray) TODO: need a mechanism for associating
    // a delta ray cluster with a muon cluster w/o assigning it to a Begin or End. Check the
    // angle of the j trajectory before and after the intersection point to determine
    // whether to split the j cluster into two clusters
    VtxStore aVtx;
    unsigned short icl, jcl, ivx;
    for(ii = 0; ii < intPts.size(); ++ii) {
      // make a vtx at the intersection point.
      // Index of the j trajectory
      jtj = intPts[ii].jTraj;
      // index of the split point on the j trajectory
      jtjPt = intPts[ii].jTrajPt;
      // the i trajectory
      itj = intPts[ii].iTraj;
      // the end of the i trajectory that may result in a vertex
      end = intPts[ii].iTrajEnd;
      // get the trajectory at the end of the itj trajectory
      if(end == 0) { itjEnd = 0; } else { itjEnd = allTraj[itj].Pts.size() - 1; }
      // make a simple estimate of the vertex position and make it. This will be
      // refined later
      aVtx.Wire = (allTraj[itj].Pts[itjEnd].Pos[0] + allTraj[jtj].Pts[jtjPt].Pos[0]) / 2;
      aVtx.WireErr = std::abs(allTraj[itj].Pts[itjEnd].Pos[0] - aVtx.Wire);
      aVtx.Time = (allTraj[itj].Pts[itjEnd].Pos[1] + allTraj[jtj].Pts[jtjPt].Pos[1]) / 2;
      aVtx.TimeErr = std::abs(allTraj[itj].Pts[itjEnd].Pos[1] - aVtx.Time);
       // convert to ticks
      aVtx.Time /= fScaleF; aVtx.TimeErr /= fScaleF;
      std::cout<<"vtx wire "<<aVtx.Wire<<" err "<<aVtx.WireErr<<" time "<<aVtx.Time<<" err "<<aVtx.TimeErr<<"\n";
      aVtx.NClusters = 2;
      aVtx.ChiDOF = 0;
      aVtx.Topo = 9;
      aVtx.CTP = clCTP;
      aVtx.Fixed = true;
      vtx.push_back(aVtx);
      // associate the i trajectory with the vertex
      icl = allTraj[itj].ClusterIndex;
//      icl = trajClusterIndex[itj];
      ivx = vtx.size() - 1;
      if(end == 0) { tcl[icl].EndVtx = ivx; } else { tcl[icl].BeginVtx = ivx; }
      // decide whether the j trajectory (and cluster) should be split or if the
      // vertex should be attached to an existing end
//      jcl = trajClusterIndex[jtj];
      jcl = allTraj[jtj].ClusterIndex;
      if(tcl[jcl].ID < 0) {
        mf::LogWarning("TC")<<"Cluster "<<tcl[jcl].ID<<" has already been split ";
        continue;
      }
      if((allTraj[jtj].Pts.size() - jtjPt) < 3) {
        // attach the vertex to the j cluster End end
        allTraj[jtj].Vtx[1] = ivx;
      } else if(jtjPt < 3) {
        // attach the vertex to the j cluster Begin end
        allTraj[jtj].Vtx[0] = ivx;
      }
      else {
        // Split point is not near an end. Check the angle to see if it
        // should be split
        unsigned short pt1, pt2;
        if(jtjPt > 2) { pt1 = jtjPt - 3; } else { pt1 = 0; }
        if(jtjPt < allTraj[jtj].Pts.size() - 4) { pt2 = jtjPt + 3; } else { pt2 = allTraj[jtj].Pts.size()-1; }
        float dang = std::abs(allTraj[jtj].Pts[pt1].Ang - allTraj[jtj].Pts[pt2].Ang);
        if(vtxprt) mf::LogVerbatim("TC")<<" kink angle at split point "<<dang;
        if(dang > fStepCrawlKinkAngCut) {
          prt = vtxprt;
          if(!SplitAllTraj(jtj, jtjPt, ivx)) {
            mf::LogWarning("TC")<<"FindTrajVertices: SplitAllTraj failed";
            return;
          }
          // fit the vertex position
          vtx[ivx].Fixed = false;
          FitVtx(ivx);
        } // dang > fStepCrawlKinkAngCut
      } // split point in the middle of jtj
    } // ii
    
  } // FindTrajVertices

  //////////////////////////////////////////
  bool TrajClusterAlg::SplitAllTraj(unsigned short itj, unsigned short pos, unsigned short ivx)
  {
    // Splits the trajectory itj in the allTraj vector into two trajectories at position pos. Splits
    // the trajectory cluster and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
    
    if(itj > allTraj.size()-1) return false;
    if(pos < 1 || pos > allTraj[itj].Pts.size()-2) return false;
    if(allTraj[itj].ClusterIndex > tcl.size()-1) return false;
    
    prt = true;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"*************************** Inside SplitAllTraj ***************************";
      PrintAllTraj(itj, USHRT_MAX);
    }
    
    unsigned short icl = allTraj[itj].ClusterIndex;
    MakeClusterObsolete(icl);
    // double check hit releasing
    for(unsigned short ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
      unsigned int hit = allTraj[itj].Pts[ipt].Hits[0];
      if(inTraj[hit] != 0) {
        mf::LogWarning("TC")<<"SplitAllTraj: MakeClusterObsolete problem on icl "<<icl<<" itj "<<itj<<" inTraj "<<inTraj[hit]<<" hit "<<plane<<":"<<fHits[hit]->WireID().Wire<<":"<<(int)fHits[hit]->PeakTime();
        return false;
      }
    } // check
    
    // flag the old traj as obsolete
    allTraj[itj].ProcCode = USHRT_MAX;
    
    work.Pts.clear();
    work.Pts.insert(work.Pts.end(), allTraj[itj].Pts.begin(), allTraj[itj].Pts.begin()+pos);
    work.ProcCode = 911;
    StoreWork();
    if(prt) {
      mf::LogVerbatim("TC")<<"Check DS end";
      PrintTrajectory(work, USHRT_MAX);
    }
    // associate the vertex if the one supplied is valid
    if(ivx < vtx.size()) {
      tcl[tcl.size()-1].BeginVtx = ivx;
      work.Vtx[1] = ivx;
    }
    // Do the same for the US end
    work.Pts.clear();
    work.Pts.insert(work.Pts.end(), allTraj[itj].Pts.begin()+pos, allTraj[itj].Pts.end());
    work.ProcCode = 911;
    if(prt) {
      mf::LogVerbatim("TC")<<"Check US end";
      PrintTrajectory(work, USHRT_MAX);
    }
    StoreWork();
    if(ivx < vtx.size()) {
      tcl[tcl.size()-1].EndVtx = ivx;
      work.Vtx[0] = ivx;
    }
    return true;
  } // SplitAllTraj
*/
  
  //////////////////////////////////////////
  void TrajClusterAlg::TrajTrajDOCA(Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
  {
    float best = minSep * minSep;
    float dw, dt, dp2;
    unsigned short i1, i2;
    for(i1 = 0; i1 < tp1.Pts.size(); ++i1) {
      for(i2 = 0; i2 < tp2.Pts.size(); ++i2) {
        dw = tp1.Pts[i1].Pos[0] - tp2.Pts[i2].Pos[0];
        dt = tp1.Pts[i1].Pos[1] - tp2.Pts[i2].Pos[1];
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
  void TrajClusterAlg::StepCrawl(bool& success)
  {
    // Crawl along the direction specified in the traj vector in steps of size step
    // (wire spacing equivalents). Find hits between the last trajectory point and
    // the last trajectory point + step. A new trajectory point is added if hits are
    // found. Crawling continues until no signal is found for two consecutive steps
    // or until a wire or time boundary is reached.
    
    success = false;
    if(work.Pts.size() == 0) return;
    // ensure that the direction is defined
    unsigned short lastPt = work.Pts.size() - 1;
    TrajPoint tp = work.Pts[lastPt];

    unsigned short iwt;
    
    unsigned int step;
    float stepSize;
    float bigChgDiffCut = 3 * fChgDiffCut;
    float maxPos0 = fNumWires;
    float maxPos1 = fMaxTime * fScaleF;
    
    // count number of steps taken with no trajectory point added
    unsigned short nMissedStep = 0;
    unsigned short missedStepCut = SetMissedStepCut(tp);
    
    if(prt) mf::LogVerbatim("TC")<<"Start StepCrawl with missedStepCut "<<missedStepCut<<" fChgDiffCut "<<fChgDiffCut;

    bool hitsPresent, updateSuccess, keepGoing, stopOnKink;
    unsigned short killPts;
    unsigned short noSignal = 0;
    for(step = 0; step < 1000; ++step) {
      if(IsLargeAngle(tp)) { stepSize = 1; } else { stepSize = 1 / tp.Dir[0]; }
      // move the position by one step in the right direction
      for(iwt = 0; iwt < 2; ++iwt) tp.Pos[iwt] += tp.Dir[iwt] * stepSize;
      // hit the boundary of the TPC?
      if(tp.Pos[0] < 0 || tp.Pos[0] > maxPos0) break;
      if(tp.Pos[1] < 0 || tp.Pos[1] > maxPos1) break;
      // remove the old hits
      tp.Hits.clear();
      // look for new hits
      AddTrajHits(tp, hitsPresent);
      if(!hitsPresent) {
        // check the number of missed wires
        lastPt = work.Pts.size() - 1;
        if(prt)  mf::LogVerbatim("TC")<<" No signal. Missed wires "<<std::abs(tp.Pos[0] - work.Pts[lastPt].Pos[0])<<" user cut "<<fMaxWireSkip;
        if(std::abs(tp.Pos[0] - work.Pts[lastPt].Pos[0]) < fMaxWireSkip) continue;
        break;
      }
      if(tp.Hits.size() == 0) {
        ++nMissedStep;
        // check for a signal at this TP
        if(!SignalAtTp(tp)) ++noSignal;
        if(prt) mf::LogVerbatim("TC")<<" nMissedStep "<<nMissedStep<<" missedStepCut "<<missedStepCut<<" noSignal count "<<noSignal;
        // Break if we took too many steps without seeing a signal
        if(IsLargeAngle(tp)) {
          if(noSignal > 100) break;
        } else {
          if(noSignal > 2) break;
        }
        // Break if we took too many steps without adding a traj point
        if(nMissedStep > missedStepCut) {
          if(prt) mf::LogVerbatim("TC")<<" Too many steps w/o traj point ";
          break;
        }
        // Otherwise keep stepping
        continue;
      } // tp.Hits.size() == 0
      // Have new traj hits. Add the trajectory point and update
      tp.Step = step;
      nMissedStep = 0;
      noSignal = 0;
      // this cut for the next step
      missedStepCut = SetMissedStepCut(tp);
      work.Pts.push_back(tp);
      UpdateWork(updateSuccess);
      if(!updateSuccess) break;
      lastPt = work.Pts.size()-1;
      // Quit if we are starting out poorly. This can happen on the 2nd trajectory point
      // when a hit is picked up in the wrong direction
      if(work.Pts.size() == 2 && std::signbit(work.Pts[1].Dir[0]) != std::signbit(work.Pts[0].Dir[0])) {
        if(prt) {
          mf::LogVerbatim("TC")<<" Picked wrong hit on 2nd traj point. Drop trajectory.";
          PrintTrajectory(work, USHRT_MAX);
        }
        return;
      }
      if(work.Pts.size() == 3) {
        // ensure that the last hit added is in the same direction as the first two.
        // This is a simple way of doing it
        if(prt) mf::LogVerbatim("TC")<<"StepCrawl: Third TP. 0-2 sep "<<TrajPointHitSep2(work.Pts[0], work.Pts[2])<<" 0-1 sep "<<TrajPointHitSep2(work.Pts[0], work.Pts[1]);
        if(TrajPointHitSep2(work.Pts[0], work.Pts[2]) < TrajPointHitSep2(work.Pts[0], work.Pts[1])) return;
      } // work.Pts.size() == 3
      // See if the Chisq/DOF exceeds the maximum.
      // UpdateWork should have reduced the number of points fit
      // as much as possible for this pass, so this trajectory is in trouble
      if(work.Pts[lastPt].FitChi > fMaxChi) {
        if(prt) mf::LogVerbatim("TC")<<"   bad FitChi "<<work.Pts[lastPt].FitChi<<" cut "<<fMaxChi<<" Dropping trajectory";
        success = false;
        return;
      }
      // assume that we aren't going to kill the point we just added, or any
      // of the previous points...
      killPts = 0;
      // assume that we should keep going after killing points
      keepGoing = true;
      // Bad 3 TP fit chisq
      if(work.Pts[lastPt].TP3Chi > fTP3ChiCut) {
        if(prt) mf::LogVerbatim("TC")<<"   bad TP3Chi "<<work.Pts[lastPt].TP3Chi<<" cut "<<fTP3ChiCut;
        killPts = 1;
      }
      // Look for two consecutive large Deltas
      if(killPts == 0 && lastPt > 3 && tp.Delta > fMaxDeltaJump * work.Pts[lastPt-2].Delta) {
        // check the previous Delta ratio
        if(work.Pts[lastPt-1].Delta > fMaxDeltaJump * work.Pts[lastPt-2].Delta) killPts = 2;
        // check for not nearby hits
        if(prt) mf::LogVerbatim("TC")<<"   NCloseNotUsed "<<tp.NCloseNotUsed<<" previous "<<work.Pts[lastPt-1].NCloseNotUsed;
        if(tp.NCloseNotUsed < 2 && work.Pts[lastPt-1].NCloseNotUsed < 2) killPts = 0;
        if(killPts > 0) keepGoing = false;
        if(prt) mf::LogVerbatim("TC")<<"   large Deltas? killPts "<<killPts<<" keepGoing? "<<keepGoing;
      }
      // check the charge ratio
      if(killPts == 0 && !fStudyMode) {
        // Kill this point if it has a high charge and the previous did as well
        if(tp.ChgDiff > fChgDiffCut && work.Pts[lastPt].ChgDiff > fChgDiffCut) killPts = 1;
        // or if it has extraordinarily high charge
        if(tp.ChgDiff > bigChgDiffCut) killPts = 1;
        // or if the charge difference and TP delta are large
        if(tp.ChgDiff > fChgDiffCut && tp.TP3Chi > 2) killPts = 1;
        if(prt) mf::LogVerbatim("TC")<<"   large ChgDiff? "<<tp.ChgDiff<<" killPts "<<killPts;
      } // fHitChgRMS > 0
      // check for a kink. Stop crawling if one is found
      if(GottaKink()) {
        stopOnKink = true;
        break;
      }
      // update the local tp unless we have killing to do
      if(killPts == 0) {
        tp = work.Pts[lastPt];
        if(prt) PrintTrajectory(work, lastPt);
      } else {
        // shorten the trajectory by the desired number of points
        // and correct the projected position.
        // release the hits
        FlagWorkHits(0);
        // truncate the work trajectory
        work.Pts.resize(lastPt+1-killPts);
        lastPt = work.Pts.size() - 1;
        float nSteps = (float)(step - work.Pts[lastPt].Step);
        if(prt) mf::LogVerbatim("TC")<<"TRP   killing "<<killPts<<" after "<<nSteps<<" steps from prev TP.  Current tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1];
        // re-assign the hits
        FlagWorkHits(-3);
        // copy to the local trajectory point
        tp = work.Pts[lastPt];
        // move the position
        tp.Pos[0] += nSteps * tp.Dir[0];
        tp.Pos[1] += nSteps * tp.Dir[1];
        if(!IsLargeAngle(tp)) {
          // put the TP at the last wire position + a bit
          float wire = fHits[tp.Hits[0]]->WireID().Wire + 0.5;
          float dw = wire - tp.Pos[0];
          tp.Pos[0] = wire;
          tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
        }
        if(prt) mf::LogVerbatim("TC")<<"  New tp.Pos     "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" ticks "<<(int)tp.Pos[1]/fScaleF;
        if(!keepGoing) break;
      }
    } // step
    
    if(prt) mf::LogVerbatim("TC")<<"End StepCrawl with "<<step<<" steps. work size "<<work.Pts.size();
    
    // release all the hits since we have stopped crawling. Doing this
    // here simplifies the task of truncating the trajectory
    FlagWorkHits(0);
    // check the cluster quality and possibly truncate it
    if(!stopOnKink) CheckWork();
    

//    if(prt) PrintTrajectory(work, USHRT_MAX);
    
    if(!CheckAllHitsFlag()) return;
    
    success = true;
    return;
    
  } // StepCrawl
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::CheckWork()
  {
    // look for a gap in the trajectory with only one or two
    // TPs at the end
    
    // Ignore LA trajectories
    if(IsLargeAngle(work.Pts[0])) return;
    
    // Compare the number of steps taken per TP near the beginning and
    // at the end
    
    if(work.Pts.size() < 4) return;
    short nStepBegin = work.Pts[2].Step - work.Pts[1].Step;
    short nStepEnd;
    unsigned short lastPt = work.Pts.size() - 1;
    unsigned short newSize = work.Pts.size();
    for(unsigned short ipt = lastPt; ipt > lastPt - 2; --ipt) {
      nStepEnd = work.Pts[ipt].Step - work.Pts[ipt - 1].Step;
      if(nStepEnd > 3 * nStepBegin) newSize = ipt;
      if(prt) mf::LogVerbatim("TC")<<"CheckWork "<<ipt<<" Begin "<<nStepBegin<<" End "<<nStepEnd<<" newSize "<<newSize;
    }
    if(newSize < work.Pts.size()) work.Pts.resize(newSize);
    
  } // CheckWorkEnd

  ////////////////////////////////////////////////
  void TrajClusterAlg::IsHitFlagSet(std::string someText, unsigned int iht)
  {
    if(inTraj[iht] == -3) {
      mf::LogVerbatim myprt("TC");
      myprt<<someText<<" "<<iht<<" flag set. ";
      if(work.Pts.size() > 0) {
        unsigned int iht = work.Pts[0].Hits[0];
        myprt<<"Seed hit "<<fPlane<<":"<<fHits[iht]->WireID().Wire<<":"<<(int)fHits[iht]->PeakTime();
      } else {
        myprt<<"No Pts in work";
      }
    } // inTraj -3
  } // IsHitFlagSet
  
  ////////////////////////////////////////////////
  void TrajClusterAlg::FindHit(std::string someText, unsigned int iht)
  {
    // finds a hit in work
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      TrajPoint& tp = work.Pts[ipt];
      if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) {
        unsigned int sht = work.Pts[0].Hits[0];
        mf::LogVerbatim("TC")<<someText<<" work hit "<<iht<<" inTraj "<<inTraj[iht]<<" Seed hit "<<fPlane<<":"<<fHits[sht]->WireID().Wire<<":"<<(int)fHits[sht]->PeakTime();
        return;
      }
    } // ipt
    // look in allTraj
    for(unsigned short itj = 0; itj < allTraj.size(); ++itj) {
      for(unsigned short ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
        TrajPoint& tp = allTraj[itj].Pts[ipt];
        if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) {
          unsigned int sht = work.Pts[0].Hits[0];
          mf::LogVerbatim("TC")<<someText<<" allTraj hit "<<itj<<" "<<iht<<" inTraj "<<inTraj[iht]<<" Seed hit "<<fPlane<<":"<<fHits[sht]->WireID().Wire<<":"<<(int)fHits[sht]->PeakTime();
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
      if(inTraj[iht] == -3) {
        itsBad = true;
        break;
      }
    } // iht
    
    if(!itsBad) return true;
    
    // print out the information
    mf::LogVerbatim myprt("TC");
    myprt<<"Not released hits ";
    for(unsigned int iht = firsthit; iht < lasthit; ++iht) {
      if(inTraj[iht] == -3) {
        myprt<<" "<<iht<<"_"<<fPlane<<":"<<fHits[iht]->WireID().Wire<<":"<<(int)fHits[iht]->PeakTime();
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
  bool TrajClusterAlg::GottaKink()
  {
    // Does a local fit of just-added traj points to identify a kink while step crawling.
    // Truncates the traj vector and returns true if one is found.
    
    if(work.Pts.size() < 10) return false;
    
    // don't look for a kink if the number of points fit is very low
    unsigned short toIndex = work.Pts.size() - 1;
    if(work.Pts[toIndex].NTPsFit < 6) return false;
    
    // The kink angle cut will be scaled by the normalized average hit resolution so this
    // must be known before proceeding
//    if(work.Pts[work.Pts.size()-1].AveHitRes < 0) return false;
    
    unsigned short nTPF = 4;
    
    unsigned short fromIndex = work.Pts.size() - nTPF;
     TrajPoint tp1;
    FitTrajMid(fromIndex, toIndex, tp1);
    if(prt) mf::LogVerbatim("TC")<<"Fit1 "<<fromIndex<<" "<<toIndex<<" Ang "<<tp1.Ang<<" chi "<<tp1.FitChi<<" NTPsFit "<<work.Pts[toIndex].NTPsFit;
    if(tp1.FitChi > 900) {
      unsigned int hit = tp1.Hits[0];
      mf::LogWarning("TC")<<"GottaKink: Bad traj fit1. Seed hit "<<fPlane<<":"<<fHits[hit]->WireID().Wire<<":"<<(int)fHits[hit]->PeakTime()<<" Pts size "<<work.Pts.size()<<" fromIndex "<<fromIndex<<" to "<<toIndex;
      PrintTrajectory(work, USHRT_MAX);
      return false;
    }
    fromIndex -= nTPF;
    toIndex -= nTPF;
    TrajPoint tp2;
    FitTrajMid(fromIndex, toIndex, tp2);
    if(tp2.FitChi > 900) {
      unsigned int hit = tp2.Hits[0];
      mf::LogWarning("TC")<<"GottaKink: Bad traj fit2 "<<fromIndex<<" "<<toIndex<<" work.Pts size "<<work.Pts.size()<<" Seed hit "<<fPlane<<":"<<fHits[hit]->WireID().Wire<<":"<<(int)fHits[hit]->PeakTime();
      return false;
    }
    float dang = std::abs(tp1.Ang - tp2.Ang);
    if(prt) mf::LogVerbatim("TC")<<"GottaKink: Ang1 "<<tp1.Ang<<" Ang2 "<<tp2.Ang<<" dang "<<dang<<" fKinkAngCut "<<fKinkAngCut;
    if(dang < fKinkAngCut) return false;
    
    unsigned short kinkIndex = work.Pts.size()-nTPF+1;
    // found a kink so truncate
    if(prt) mf::LogVerbatim("TC")<<" >>> resize work to size "<<kinkIndex;
    FlagWorkHits(0);
    work.Pts.resize(kinkIndex);
    return true;
    
  } // GottaKink

 
  //////////////////////////////////////////
  unsigned short TrajClusterAlg::SetMissedStepCut(TrajPoint const& tp)
  {
    // maximum of 100 WSE's
    unsigned short itmp = 100;
    if(std::abs(tp.Dir[0]) > 0.01) itmp = 1 + 4 * std::abs(1/tp.Dir[0]);
//    if(itmp > 10) itmp = 10;
    if(itmp < 2) itmp = 2;
    return itmp;
  }
  
  //////////////////////////////////////////
  void TrajClusterAlg::UpdateWork( bool& success)
  {
    // A new trajectory point was put on the traj vector with hits by StepCrawl.
    // This routine updates the last trajectory point fit, average hit rms, etc.
    // StepCrawl will decide what to do next with this information.
    
    success = false;
    
    if(work.Pts.size() == 0) return;
    
    // always update the charge difference
    UpdateWorkChgDiff();
    if(work.Pts.size() < 2) return;
    
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];

    UpdateWorkDelta();

    if(work.Pts.size() == 2) {
      // Handle the second trajectory point. No error calculation. Just update
      // the position and direction
      lastTP.NTPsFit = 2;
      FitTraj();
      lastTP.FitChi = 0.01;
      lastTP.AngErr = work.Pts[0].AngErr;
      if(prt) mf::LogVerbatim("TC")<<"UpdateWork: Second traj point pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1];
    } else if(work.Pts.size() == 3) {
       // Third trajectory point. Keep it simple
      lastTP.NTPsFit = 2;
      FitTraj();
    } else {
      lastTP.NTPsFit += 1;
      FitTraj();
      if(prt) mf::LogVerbatim("TC")<<"UpdateWork: First fit "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" FitChi "<<lastTP.FitChi;
      // check for a failure
      if(lastTP.FitChi > 900) return;
      // reduce the number of points fit to keep Chisq/DOF < 2 adhering to the pass constraint
      unsigned short newNTPSFit = lastTP.NTPsFit;
      while(lastTP.FitChi > 2 && lastTP.NTPsFit > 2) {
        if(lastTP.NTPsFit > 3) newNTPSFit -= 2;
        else if(lastTP.NTPsFit == 3) newNTPSFit = 2;
        lastTP.NTPsFit = newNTPSFit;
        if(prt) mf::LogVerbatim("TC")<<"  Bad FitChi "<<lastTP.FitChi<<" Reduced NTPsFit to "<<lastTP.NTPsFit;
        FitTraj();
        if(lastTP.NTPsFit <= fMinNPtsFit[fPass]) break;
      }
      if(prt) mf::LogVerbatim("TC")<<"  Fit done. Chi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;
    } // work.Pts.size
    
/*
    // Move the TP to the wire position for small angle clusters
    if(!IsLargeAngle(lastTP)) {
      float wire = fHits[lastTP.Hits[0]]->WireID().Wire;
      float dw = wire - lastTP.Pos[0];
      lastTP.Pos[0] = wire;
      lastTP.Pos[1] += dw * lastTP.Dir[1] / lastTP.Dir[0];
      if(prt) mf::LogVerbatim("TC")<<"  Moved tp to wire position by dw = "<<dw<<" Pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1];
    } // !LargeAngle
*/
    // update the first traj point if this point is in the fit
    if(lastTP.NTPsFit == work.Pts.size()) {
      work.Pts[0].AveChg = lastTP.AveChg;
      work.Pts[0].Ang = lastTP.Ang;
      work.Pts[0].AngErr = lastTP.AngErr;
    }
    // Calculate TP3Chi
    if(work.Pts.size() > 2) {
      lastTP.TP3Chi = 0;
      TrajPoint tp;
      FitTrajMid(lastPt, lastPt - 2, tp);
      // put this in the current TP. It will be over-written if a new TP is found
      lastTP.TP3Chi = tp.FitChi;
      // put it in the previous TP
      work.Pts[lastPt - 1].FitChi = tp.FitChi;
      // put it in the first 2 TPs if we are on the third
      if(work.Pts.size() == 3) {
        work.Pts[0].FitChi = tp.FitChi;
        work.Pts[1].FitChi = tp.FitChi;
        work.Pts[0].TP3Chi = tp.FitChi;
        work.Pts[1].TP3Chi = tp.FitChi;
      }
    } // work.Pts.size() > 2
    
    // Calculate the average TP3Chi
    work.AveTP3Chi = 0;
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) work.AveTP3Chi += work.Pts[ipt].TP3Chi;
    float fcnt = work.Pts.size();
    work.AveTP3Chi /= fcnt;
    
    // Estimate the multiple scattering angle
    float sum = 0, sum2 = 0;
    work.ThetaMCS = 1;
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      sum += work.Pts[ipt].Ang;
      sum2 += work.Pts[ipt].Ang * work.Pts[ipt].Ang;
    }
    sum /= fcnt;
    if(fcnt > 2) work.ThetaMCS = sqrt((sum2 - fcnt * sum * sum) / (fcnt-1));
    if(prt) mf::LogVerbatim("TC")<<"    AveTP3Chi "<<work.AveTP3Chi<<" ThetaMCS "<<work.ThetaMCS;

    success = true;
    return;
    
  } // UpdateWork

  //////////////////////////////////////////
  void TrajClusterAlg::UpdateWorkDelta()
  {
    // Find Delta for the last trajectory point.
    
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];

    lastTP.Delta = PointTrajDOCA(lastTP.HitPos[0], lastTP.HitPos[1] ,lastTP);
    
    if(work.Pts.size() < 3) return;
    
    float ave = 0, sum2 = 0, fcnt = 0;
    unsigned short ii, ipt;
    unsigned short npts = work.Pts.size() - 1;
    if(npts > fNPtsAve) npts = fNPtsAve;
    for(ii = 0; ii < npts; ++ii) {
      ipt = work.Pts.size() -1 - ii;
      ave += work.Pts[ipt].Delta;
      sum2 += work.Pts[ipt].Delta * work.Pts[ipt].Delta;
      ++fcnt;
    }
    ave /= fcnt;
    if(fcnt > 2) {
      lastTP.DeltaRMS = sqrt((sum2 - fcnt * ave * ave) / (fcnt - 1));
      // this is not well defined for trajectories with direction at
      // +/- pi/2. Create a lower limit which is sqrt(12). TODO check this
      if(IsLargeAngle(lastTP) && lastTP.DeltaRMS < 0.083) lastTP.DeltaRMS = 0.083;
    }

  } // UpdateWorkDelta

  
  //////////////////////////////////////////
  void TrajClusterAlg::FitTrajMid(unsigned short fromIndex, unsigned short toIndex, TrajPoint& tp)
  {
    // Fit the work trajectory points from fromIndex to (and including) toIndex. The origin
    // of the fit position is set to work.Pts[fromIndex].Pos. Fit results are stashed in
    // TrajPoint tp. This copies the trajectory point from work into tp. Note that the
    // AngleErr is the error on the fitted points (not the angle error for the last point)
    
    tp.FitChi = 9999;
    
    unsigned short lastPt = work.Pts.size() - 1;

    if(fromIndex > lastPt || toIndex > lastPt) return;
    
    tp = work.Pts[fromIndex];
    
    unsigned short lo, hi;
    if(fromIndex < toIndex) {
      lo = fromIndex; hi = toIndex;
    } else {
      lo = toIndex; hi = fromIndex;
    }
    unsigned short nTPF = hi - lo + 1;
    tp.NTPsFit = nTPF;
    if(nTPF < 1) return;
    ++hi;
    
    std::vector<float> x(nTPF), y(nTPF), yerr2(nTPF);
    
    // Rotate the traj hit position into the coordinate system defined by the
    // fromIndex traj point, where x = along the trajectory, y = transverse
    float rotAngle = work.Pts[fromIndex].Ang;
    float cs = cos(-rotAngle);
    float sn = sin(-rotAngle);
//    if(prt) mf::LogVerbatim("TC")<<"FTM: Ang "<<rotAngle<<" cs "<<cs<<" sn "<<sn<<" work size "<<work.Pts.size()<<" StepDir "<<work.StepDir;
    
    // fit origin is the position of the fromIndex trajectory point
    std::array<float, 2> dir, origin = work.Pts[fromIndex].Pos;
    float xx, yy;
    
    // protect against cases where we have a number of points with the same XX value
    unsigned short indx, nsame = 0;
    // Find the average charge
    tp.Chg = 0;
    float sum2 = 0;
//    terr = std::abs(fScaleF * fHitErrFac * fAveHitRMS);
//    if(terr < 0.01) terr = 0.01;
    for(unsigned short ipt = lo; ipt < hi; ++ipt) {
      indx = ipt - lo;
      if(ipt > work.Pts.size()-1 || indx > nTPF-1) {
        mf::LogError("TC")<<"bad ipt "<<ipt<<" work size "<<work.Pts.size()<<" or bad indx "<<indx<<" nTPF "<<nTPF;
        return;
      }
      xx = work.Pts[ipt].HitPos[0] - origin[0];
      yy = work.Pts[ipt].HitPos[1] - origin[1];
      x[indx] = cs * xx - sn * yy;
      y[indx] = sn * xx + cs * yy;
      // Time error
      // Wire error^2 is a constant [1/sqrt(12)]^2 in WSE units = 0.0833
      yerr2[indx] = std::abs(sn) * 0.0833 + std::abs(cs) * work.Pts[ipt].HitsTimeErr2;
      
//      if(prt) mf::LogVerbatim("TC")<<"FTM: ipt "<<ipt<<" indx "<<indx<<" xx "<<xx<<" yy "<<yy<<" x "<<x[indx]<<" y "<<y[indx]<<" err "<<yerr2[indx]<<" Chg "<<work.Pts[ipt].Chg;
      if(indx > 0 && std::abs(xx - x[indx-1]) < 1.E-3) ++nsame;
      tp.Chg += work.Pts[ipt].Chg;
      sum2 += work.Pts[ipt].Chg * work.Pts[ipt].Chg;
    } // ii
    float fcnt = nTPF;
    tp.Chg /= fcnt;
    float arg = sum2 - fcnt * tp.Chg * tp.Chg;
    if(arg > 1 && nTPF > 1) {
      tp.ChgDiff = sqrt(arg / (fcnt - 1));
    } else {
      tp.ChgDiff = 0.3;
    }
    
    if(nTPF < 4 && nsame > 0) return;
    
    float intcpt, slope, intcpterr, slopeerr, chidof;
    fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
//    if(prt) mf::LogVerbatim("TC")<<" LinFit chidof "<<chidof;
    if(chidof < 0.01) chidof = 0.01;
    tp.FitChi = chidof;
    if(chidof > 900) return;
    
    // calculate the new direction vector in the (xx, yy) coordinate system
    float newang = atan(slope);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
//    if(prt) mf::LogVerbatim("TC")<<"newdir XX,YY system "<<dir[0]<<" "<<dir[1]<<" slope "<<slope;
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAngle);
    sn = sin(rotAngle);
//    if(prt) mf::LogVerbatim("TC")<<" lastPT ang "<<rotAngle<<" dir "<<cs<<" "<<sn;
    tp.Dir[0] = cs * dir[0] - sn * dir[1];
    tp.Dir[1] = sn * dir[0] + cs * dir[1];
    // decide whether to reverse the direction
    if(work.StepDir < 0) {
      tp.Dir[0] = -tp.Dir[0];
      tp.Dir[1] = -tp.Dir[1];
    }
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    // a good enough approximation since the trajectory can't change
    // too much from point to point
    tp.AngErr = std::abs(atan(slopeerr));
//    if(prt) mf::LogVerbatim("TC")<<"  W,T dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" ang "<<tp.Ang;
    // rotate (0, intcpt) into W,T
    tp.Pos[0] = -sn * intcpt + origin[0];
    tp.Pos[1] =  cs * intcpt + origin[1];
    
//    if(prt) mf::LogVerbatim("TC")<<"FTM: fromIndex "<<fromIndex<<" Pos1 "<<tp.Pos[1]<<" Ang "<<tp.Ang<<" chi "<<chidof<<" old Ang "<<work.Pts[fromIndex].Ang;

  } // FitTrajMid
  
  //////////////////////////////////////////
  void TrajClusterAlg::FitTraj()
  {
    // Fit the trajectory to a line. The hit positions are rotated into a coordinate
    // system with X along the direction of the last trajectory point and Y transverse.
    // The number of trajectory points fitted (nTPF) is determined by the calling routine.
    
    if(work.Pts.size() < 2) return;
    
    unsigned short lastPt = work.Pts.size() - 1;
    unsigned short ipt, ii;
    TrajPoint& lastTP = work.Pts[lastPt];
    
    if(lastTP.NTPsFit == 2) {
      lastTP.Pos[0] = lastTP.HitPos[0];
      lastTP.Pos[1] = lastTP.HitPos[1];
      float dw = lastTP.HitPos[0] - work.Pts[0].Pos[0];
      float dt = lastTP.HitPos[1] - work.Pts[0].Pos[1];
      float nrm = sqrt(dw * dw + dt * dt);
      lastTP.Dir[0] = dw / nrm; lastTP.Dir[1] = dt / nrm;
      lastTP.Ang = atan2(lastTP.Dir[1], lastTP.Dir[0]);
      unsigned int iht = work.Pts[0].Hits[0];
      float minrms = fHits[iht]->RMS();
      iht = work.Pts[1].Hits[0];
      if(fHits[iht]->RMS() < minrms) minrms = fHits[iht]->RMS();
      lastTP.AngErr = atan(2 * fHitErrFac * fScaleF * minrms);
      lastTP.FitChi = 0.01;
//      std::cout<<"chk2 from "<<work.Pts[lastPt-1].Pos[0]<<" "<<work.Pts[lastPt-1].Pos[1]<<" to "<<work.Pts[lastPt].Pos[0]<<" "<<work.Pts[lastPt].Pos[1];
//      std::cout<<" prev Ang "<<work.Pts[lastPt-1].Ang<<" new "<<work.Pts[lastPt].Ang<<" dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<"\n";
      return;
    }
    
    // try to fit with the number of points specified
    unsigned short nTPF = lastTP.NTPsFit;
    // truncate if this is too many points
    if(nTPF > work.Pts.size()) {
      mf::LogWarning("TC")<<"FitTraj: NTPsFit is too large "<<nTPF<<" work.Pts size "<<work.Pts.size()<<" Truncating...";
      nTPF = work.Pts.size();
    }
    
    std::vector<float> x(nTPF), y(nTPF), yerr2(nTPF);
    
    // Rotate the traj hit position into the coordinate system defined by the
    // current traj point, where X = along the trajectory, Y = transverse
    float rotAng = work.Pts[lastPt].Ang;
    float cs = cos(-rotAng);
    float sn = sin(-rotAng);
    
    // fit origin is the hit position of the current trajectory point
    std::array<float, 2> dir, origin = work.Pts[lastPt].HitPos;
    float xx, yy;
    
//    if(prt) mf::LogVerbatim("TC")<<"Enter with dir "<<work.Pts[lastPt].Dir[0]<<" "<<work.Pts[lastPt].Dir[1];

    // protect against cases where we have a number of points with the same X value
    unsigned short nsame = 0;
    // weight by the average charge
    float qAve = 0, wght;
    for(ii = 0; ii < nTPF; ++ii) qAve += work.Pts[lastPt-ii].Chg;
    qAve /= (float)(nTPF - 1);
    for(ii = 0; ii < nTPF; ++ii) {
      ipt = lastPt - ii;
      xx = work.Pts[ipt].HitPos[0] - origin[0];
      yy = work.Pts[ipt].HitPos[1] - origin[1];
      x[ii] = cs * xx - sn * yy;
      y[ii] = sn * xx + cs * yy;
      // Time error
      // Wire error^2 is a constant [1/sqrt(12)]^2 in WSE units = 0.0833
      // TODO These should be added in quadrature...
      yerr2[ii] = std::abs(sn) * 2 * 0.0833 + std::abs(cs) * work.Pts[ipt].HitsTimeErr2;
      wght = std::abs((work.Pts[ipt].Chg / qAve) - 1);
      if(wght < 0.01) wght = 0.01;
      yerr2[ii] *= wght;
//      if(yerr2[ii] < 0.001) yerr2[ii] = 0.001;
//      if(prt) mf::LogVerbatim("TC")<<"pt "<<ipt<<" xx "<<xx<<" yy "<<yy<<" x "<<x[ipt]<<" y "<<y[ipt]<<" err "<<yerr2[ipt];
      if(ii > 0 && std::abs(xx - x[ii-1]) < 1.E-3) ++nsame;
    } // ii
    if(nTPF < 4 && nsame > 0) return;
    
    float intcpt, slope, intcpterr, slopeerr, chidof;
    fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    if(chidof < 0.01) chidof = 0.01;
    lastTP.FitChi = chidof;
    if(chidof > 900) return;
    
    // calculate the new direction vector in the (xx, yy) coordinate system
    // TODO: Store this with the trajectory so it doesn't need to be calculated on the next step?
    float newang = atan(slope);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAng);
    sn = sin(rotAng);
    lastTP.Dir[0] = cs * dir[0] - sn * dir[1];
    lastTP.Dir[1] = sn * dir[0] + cs * dir[1];
    lastTP.Ang = atan2(lastTP.Dir[1], lastTP.Dir[0]);
    lastTP.AngErr = std::abs(atan(slopeerr));
//    if(prt) mf::LogVerbatim("TC")<<"  W,T dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" ang "<<lastTP.Ang;
    // rotate (0, intcpt) into W,T
    lastTP.Pos[0] = -sn * intcpt + origin[0];
    lastTP.Pos[1] =  cs * intcpt + origin[1];
//    if(prt) mf::LogVerbatim("TC")<<"  intcpt "<<intcpt<<" new pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1];
    
    if(prt) mf::LogVerbatim("TC")<<"Fit: "<<nTPF<<" pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<" dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" chi "<<chidof;
    
  } // FitTraj
  
  
  //////////////////////////////////////////
  void TrajClusterAlg::UpdateWorkChgDiff()
  {
    // calculate trajectory charge using nTrajPoints at the leading edge of the
    // trajectory
    
    if(work.Pts.size() == 0) return;
    
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];
    lastTP.ChgDiff = 99;
    
    if(fHitChgRMS < 1) {
      fHitChgRMS = 100;
      mf::LogWarning("TC")<<"fHitChgRMS not defined. Setting it to "<<fHitChgRMS;
    } // error checking
    
    lastTP.ChgDiff = 0.1;
    
    float sum, sum2, fcnt;
    unsigned short ii;
    if(work.Pts.size() < 5) {
      // just average for short trajectories
      sum = 0;
      fcnt = work.Pts.size() - 1;
      for(ii = 1; ii < work.Pts.size(); ++ii) sum += work.Pts[ii].Chg;
      if(fcnt > 0) {
        lastTP.AveChg = sum / fcnt;
        lastTP.ChgDiff = (lastTP.Chg - lastTP.AveChg) / fHitChgRMS;
        if(std::abs(lastTP.ChgDiff) < 0.01) lastTP.ChgDiff = 0.01;
      }
      return;
    }
    
    unsigned short ipt, cnt;
    
    if(lastTP.NTPsFit > work.Pts.size()) {
      mf::LogError("TC")<<"UpdateWorkChgDiff: Invalid lastTP.NTPsFit = "<<lastTP.NTPsFit<<" work.Pts size "<<work.Pts.size();
      return;
    }
    
    // sort by charge and find the average using the lowest (<60%) charge hits at the
    // leading edge. Limit the number of trajectory points to allow for energy loss
    std::vector<float> chg;
    // Use the full length of the trajectory except for the first point
    unsigned short npts = lastPt;
    // no more than fNPtsAve points however
    if(npts > fNPtsAve) npts = fNPtsAve;
    for(ii = 0; ii < npts; ++ii) {
      ipt = lastPt - ii;
      chg.push_back(work.Pts[ipt].Chg);
    } // ii
    // default sort is lowest to highest
    std::sort(chg.begin(), chg.end());
    unsigned short maxcnt = (unsigned short)(0.8 * chg.size());
    if(maxcnt < 2) maxcnt = 2;
    sum = 0; sum2 = 0; cnt = 0;
    for(ipt = 0; ipt < maxcnt; ++ipt) {
      sum += chg[ipt];
      sum2 += chg[ipt] * chg[ipt];
      ++cnt;
    }
    fcnt = cnt;
    lastTP.AveChg = sum / fcnt;
    fHitChgRMS = sqrt((sum2 - fcnt * lastTP.AveChg * lastTP.AveChg) / (fcnt - 1));
    if(fHitChgRMS < 10) fHitChgRMS = 10;
    lastTP.ChgDiff = (lastTP.Chg - lastTP.AveChg) / fHitChgRMS;
    if(std::abs(lastTP.ChgDiff) < 0.01) lastTP.ChgDiff = 0.01;
    
  } // UpdateWorkChgDiff

  ////////////////////////////////////////////////
  void TrajClusterAlg::StartTraj(float fromWire, float fromTick, float toWire, float toTick)
  {
    // Start a simple (seed) trajectory going from a hit to a position (toWire, toTick).
    // The traj vector is cleared if an error occurs
    
    work.Pts.clear();
    work.CTP = fCTP;
    work.Pass = fPass;
    work.ProcCode = 901;
    work.Vtx[0] = -1;  work.Vtx[1] = -1;
    
    
    // create a trajectory point
    TrajPoint tp;
    MakeBareTrajPoint(fromWire, fromTick, toWire, toTick, tp);

    tp.AngErr = 0.1;
    if(prt) mf::LogVerbatim("TC")<<"StartTraj "<<(int)fromWire<<":"<<(int)fromTick<<" -> "<<(int)toWire<<":"<<(int)toTick<<" dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" ang "<<tp.Ang<<" angErr "<<tp.AngErr;
    work.Pts.push_back(tp);
    work.StepDir = fStepDir;
    
    // initialize tracking averages
    fHitChgRMS = fDefaultHitChgRMS;
    
  } // StartTraj
  
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
      if(tj.Pts[ipt].Hits.size() > 1) std::reverse(tj.Pts[ipt].Hits.begin(), tj.Pts[ipt].Hits.end());
    } // ipt
  }

  ////////////////////////////////////////////////
  void TrajClusterAlg::StoreWork()
  {
    
    // put trajectories in order of US -> DS
    if(work.StepDir < 0) ReverseTraj(work);
    work.NHits = 0;
    for(auto& tjPt : work.Pts) work.NHits += tjPt.Hits.size();
    allTraj.push_back(work);
    short trID = allTraj.size();
    FlagWorkHits(trID);
    if(prt) mf::LogVerbatim("TC")<<"StoreWork trID "<<trID<<" CTP "<<work.CTP<<" NHits "<<work.NHits;
    
  } // StoreWork
  

  ////////////////////////////////////////////////
  void TrajClusterAlg::MakeAllTrajClusters()
  {
    // Make clusters from all trajectories in allTraj
    
    ClusterStore cls;
    tcl.clear();
    inClus.resize(fHits.size());
    unsigned int iht;
    for(iht = 0; iht < inClus.size(); ++iht) inClus[iht] = 0;
    
    mf::LogVerbatim("TC")<<"MakeAllTrajClusters: allTraj size "<<allTraj.size();
    
    unsigned short itj, ncl, lastPt, ipt;
    ncl = 0;
    
    // Make one cluster for each trajectory. The indexing of trajectory parents
    // should map directly to cluster parents
    for(itj = 0; itj < allTraj.size(); ++itj) {
      if(allTraj[itj].ProcCode == USHRT_MAX) continue;
      if(allTraj[itj].StepDir > 0) ReverseTraj(allTraj[itj]);
      ++ncl;
      cls.ID = ncl;
      cls.CTP = allTraj[itj].CTP;
      cls.PDG = allTraj[itj].PDG;
      cls.ParentCluster = allTraj[itj].ParentTraj;
      cls.BeginWir = allTraj[itj].Pts[0].Pos[0];
      cls.BeginTim = allTraj[itj].Pts[0].Pos[1] / fScaleF;
      cls.BeginAng = allTraj[itj].Pts[0].Ang;
      cls.BeginChg = allTraj[itj].Pts[0].Chg;
      cls.BeginVtx = allTraj[itj].Vtx[0];
      lastPt = allTraj[itj].Pts.size() - 1;
      cls.EndWir = allTraj[itj].Pts[lastPt].Pos[0];
      cls.EndTim = allTraj[itj].Pts[lastPt].Pos[1] / fScaleF;
      cls.EndAng = allTraj[itj].Pts[lastPt].Ang;
      cls.EndChg = allTraj[itj].Pts[lastPt].Chg;
      cls.EndVtx = allTraj[itj].Vtx[1];
      cls.tclhits.clear();
      for(ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt)
        cls.tclhits.insert(cls.tclhits.end(), allTraj[itj].Pts[ipt].Hits.begin(), allTraj[itj].Pts[ipt].Hits.end());
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
        return;
      } //iht
      for(iht = 0; iht < cls.tclhits.size(); ++iht) inClus[cls.tclhits[iht]] = ncl;
    } // itj

    PrintAllTraj(USHRT_MAX, 0);
    PrintClusters();
    
  } // MakeAllTrajClusters

  //////////////////////////////////////////
  void TrajClusterAlg::PrintAllTraj(unsigned short itj, unsigned short ipt)
  {
    
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      mf::LogVerbatim myprt("TC");
      std::vector<unsigned int> tmp;
      myprt<<"TRJ CTP Pass Ind nPts    W:Tick     Ang   AveQ     W:T        Ang   AveQ TP3Chi  thMCS Hits/TP __Vtx__ PDG Parent TRuPDG  Prnt   KE  \n";
      for(unsigned short ii = 0; ii < allTraj.size(); ++ii) {
        if(fDebugPlane >=0 && (unsigned short)fDebugPlane != allTraj[ii].CTP) continue;
        myprt<<"TRJ"<<std::fixed;
        myprt<<std::setw(3)<<allTraj[ii].CTP;
        myprt<<std::setw(5)<<allTraj[ii].Pass;
        myprt<<std::setw(4)<<ii<<std::setw(5)<<allTraj[ii].Pts.size();
        TrajPoint tp = allTraj[ii].Pts[0];
        unsigned short itick = tp.Pos[1]/fScaleF;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 1000) myprt<<" ";
        myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.AveChg;
        unsigned short lastPt = allTraj[ii].Pts.size() - 1;
        tp = allTraj[ii].Pts[lastPt];
        itick = tp.Pos[1]/fScaleF;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 1000) myprt<<" ";
        myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<allTraj[ii].AveTP3Chi;
        myprt<<std::setw(7)<<std::setprecision(2)<<allTraj[ii].ThetaMCS;
        // find average number of hits / TP
        PutTrajHitsInVector(allTraj[ii], tmp);
        float ave = (float)tmp.size() / (float)allTraj[ii].Pts.size();
        myprt<<std::setw(8)<<std::setprecision(2)<<ave;
        myprt<<std::setw(4)<<allTraj[ii].Vtx[0];
        myprt<<std::setw(4)<<allTraj[ii].Vtx[1];
        myprt<<std::setw(6)<<allTraj[ii].PDG;
        myprt<<std::setw(6)<<allTraj[ii].ParentTraj;
        myprt<<std::setw(6)<<allTraj[ii].TruPDG;
        myprt<<std::setw(6)<<allTraj[ii].IsPrimary;
        myprt<<std::setw(7)<<(int)allTraj[ii].TruKE;
        myprt<<"\n";
      } // ii
      return;
    } // itj > allTraj.size()-1
    
    if(itj > allTraj.size()-1) return;
    
    Trajectory tj = allTraj[itj];
    
    mf::LogVerbatim("TC")<<"Print allTraj["<<itj<<"]: ClusterIndex "<<tj.ClusterIndex<<" ProcCode "<<tj.ProcCode<<" Vtx[0] "<<tj.Vtx[0]<<" Vtx[1] "<<tj.Vtx[1];
    
    PrintHeader();
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) PrintTrajPoint(ii, tj.StepDir, tj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(ipt, tj.StepDir, tj.Pts[ipt]);
    }
  } // PrintAllTraj

  
  //////////////////////////////////////////
  void TrajClusterAlg::PrintTrajectory(Trajectory& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    unsigned short first = 0;
    unsigned short last = tj.Pts.size();
    if(tPoint == USHRT_MAX) {
      mf::LogVerbatim("TC")<<"Trajectory: CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDG<<" TruPDG "<<tj.TruPDG<<" vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size();
      PrintHeader();
      for(unsigned short ipt = first; ipt < last; ++ipt) PrintTrajPoint(ipt, tj.StepDir, tj.Pts[ipt]);
    } else {
      // just print one traj point
      if(tPoint > tj.Pts.size() -1) {
        mf::LogVerbatim("TC")<<"Cant print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(tPoint, tj.StepDir, tj.Pts[tPoint]);
    }
  } // PrintTrajectory
  
  //////////////////////////////////////////
  void TrajClusterAlg::PrintHeader()
  {
    mf::LogVerbatim("TC")<<"TRP   Ind  Stp     W:T     Delta   RMS dTick   Ang   Err  Dir0  Dir1      Q  QDiff    AveQ TP3Chi FitChi NTPF NotUsd  Hits ";
  } // PrintHeader

  ////////////////////////////////////////////////
  void TrajClusterAlg::PrintTrajPoint(unsigned short ipt, short dir, TrajPoint tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<"TRP"<<std::fixed;
    myprt<<fPass;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(4)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]; // W:T
    if(tp.Pos[1] < 10) myprt<<"  "; if(tp.Pos[1] < 100) myprt<<" ";
//    myprt<<std::setw(5)<<std::setprecision(2)<<sqrt(tp.HitsTimeErr2);
//    myprt<<std::setw(6)<<(int)(tp.Pos[1]/fScaleF); // Tick
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.DeltaRMS;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.Delta/fScaleF; // dTick
    int itick = (int)(tp.Delta/fScaleF);
    if(itick < 100) myprt<<" ";
    myprt<<std::setw(5)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[0];
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[1];
    myprt<<std::setw(7)<<(int)tp.Chg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgDiff;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(7)<<std::setprecision(2)<<tp.TP3Chi;
     myprt<<std::setw(7)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NTPsFit;
    myprt<<std::setw(6)<<tp.NCloseNotUsed;
    // print the hits associated with this traj point
    for(unsigned short iht = 0; iht < tp.Hits.size(); ++iht)
      myprt<<" "<<fHits[tp.Hits[iht]]->WireID().Wire<<":"<<(int)fHits[tp.Hits[iht]]->PeakTime();
  } // PrintTrajPoint

  ////////////////////////////////////////////////
  void TrajClusterAlg::HitMultipletPosition(unsigned int iht, float& hitTick, float& deltaRms, float& qtot)
  {
    // returns the charge weighted wire, time position of hits in the multiplet which are within
    // fMultHitSep of iht
    
    unsigned int loHit = iht - fHits[iht]->LocalIndex();
    unsigned int hiHit = loHit + fHits[iht]->Multiplicity();
    qtot = 0;
    hitTick = 0;
    float hitSep = fMultHitSep * fHits[iht]->RMS();
    std::vector<unsigned int> closeHits;
    for(unsigned int jht = loHit; jht < hiHit; ++jht) {
      if(prt) mf::LogVerbatim("TC")<<"HitMultipletPosition: iht "<<iht<<" jht "<<jht<<" separation "<<std::abs(fHits[jht]->PeakTime() - fHits[iht]->PeakTime())<<" Ticks. Cut = "<<hitSep;
      if(std::abs(fHits[jht]->PeakTime() - fHits[iht]->PeakTime()) > hitSep) continue;
      qtot += fHits[jht]->Integral();
      hitTick += fHits[jht]->Integral() * fHits[jht]->PeakTime();
      closeHits.push_back(jht);
    } // jht
    hitTick /= qtot;
    if(prt) mf::LogVerbatim("TC")<<" hitTick "<<hitTick;
    
    deltaRms = sqrt(HitsTimeErr2(closeHits));

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
  void TrajClusterAlg::FlagWorkHits(short flag)
  {
    unsigned short ipt, iht;
    for(ipt = 0; ipt < work.Pts.size(); ++ipt) {
      for(iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inTraj[work.Pts[ipt].Hits[iht]] = flag;
    } // ipt
  } // FlagWorkHits

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
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<vtx[iv].CTP;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Wire<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].WireErr;
        myprt<<std::right<<std::setw(8)<<vtx[iv].Time<<" +/- ";
        myprt<<std::right<<std::setw(4)<<vtx[iv].TimeErr;
        myprt<<std::right<<std::setw(8)<<vtx[iv].ChiDOF;
        myprt<<std::right<<std::setw(5)<<vtx[iv].NClusters;
        myprt<<std::right<<std::setw(6)<<vtx[iv].Topo;
        myprt<<"    ";
        // display the cluster IDs
        for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
          if(fDebugPlane < 3 && fDebugPlane != (int)tcl[ii].CTP) continue;
          if(tcl[ii].ID < 0) continue;
          if(tcl[ii].BeginVtx == (short)iv) myprt<<std::right<<std::setw(4)<<tcl[ii].ID;
          if(tcl[ii].EndVtx == (short)iv) myprt<<std::right<<std::setw(4)<<tcl[ii].ID;
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
    MakeBareTrajPoint((float)fHits[fromHit]->WireID().Wire, fHits[fromHit]->PeakTime(),
                      (float)fHits[toHit]->WireID().Wire,   fHits[toHit]->PeakTime(), tp);
    
  } // MakeBareTrajPoint

  /////////////////////////////////////////
  void TrajClusterAlg::MakeBareTrajPoint(float fromWire, float fromTick, float toWire, float toTick, TrajPoint& tp)
  {
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
  bool TrajClusterAlg::LargeHitSep(unsigned int iht, unsigned int jht)
  {
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
    return std::abs(tp.Dir[0] < fLargeAngle);
  } // IsLargeAngle
  
  /////////////////////////////////////////
  bool TrajClusterAlg::SignalAtTp(TrajPoint const& tp)
  {
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

  
} // namespace cluster