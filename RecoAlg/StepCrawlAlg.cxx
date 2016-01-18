//////////////////////////////////////////////////////////////////////
///
/// Step crawling code used by ClusterCrawlerAlg
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
#include "RecoAlg/ClusterCrawlerAlg.h"

// TEMP for TagClusters
#include "MCCheater/BackTracker.h"


struct SortEntry{
  int index;
  int length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}

namespace cluster {
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::ClusterLoop2()
  {
    // Version 2 of ClusterLoop.
    unsigned int ii, iwire, jwire, iht, jht;
    
    unsigned int nwires = fLastWire - fFirstWire - 1;
    unsigned int ifirsthit, ilasthit, jfirsthit, jlasthit;
    float toWire, toTick;
    bool success;

    // turn off the hit width check in ClusterHitsOK
    fAveHitWidth = -1;
    for(ii = 0; ii < nwires; ++ii) {
      // decide which way to step given the sign of fStepCrawlDir
      if(fStepCrawlDir > 0) {
        // step DS
        iwire = fFirstWire + ii;
        jwire = iwire + 1;
      } else {
        // step US
        iwire = fLastWire - ii - 1;
        jwire = iwire - 1;
      }
/*
      if(iwire < fFirstWire || iwire > fLastWire - 1) {
        std::cout<<"Bad iwire indexing "<<fFirstWire<<" "<<iwire<<" "<<fLastWire<<"\n";
        exit(1);
      }
      if(jwire < fFirstWire || jwire > fLastWire - 1) {
        std::cout<<"Bad jwire indexing "<<fFirstWire<<" "<<jwire<<" "<<fLastWire<<"\n";
        exit(1);
      }
*/
      // skip bad wires or no hits on the wire
      if(WireHitRange[iwire].first < 0) continue;
      if(WireHitRange[jwire].first < 0) continue;
      ifirsthit = (unsigned int)WireHitRange[iwire].first;
      ilasthit = (unsigned int)WireHitRange[iwire].second;
      jfirsthit = (unsigned int)WireHitRange[jwire].first;
      jlasthit = (unsigned int)WireHitRange[jwire].second;
      for(iht = ifirsthit; iht < ilasthit; ++iht) {
        prt = (fDebugPlane == (int)plane && (int)iwire == fDebugWire && std::abs((int)fHits[iht].PeakTime() - fDebugHit) < 20);
        if(prt) mf::LogVerbatim("CC")<<"Found debug hit "<<plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" inClus "<<inClus[iht]<<" RMS "<<fHits[iht].RMS()<<" fAveHitWidth "<<fAveHitWidth;
        if(inClus[iht] != 0) continue;
        // BIG TEMP
//        if(!prt) continue;
        for(jht = jfirsthit; jht < jlasthit; ++jht) {
          if(inClus[iht] != 0) continue;
          if(inClus[jht] != 0) continue;
          if(prt) mf::LogVerbatim("CC")<<" +++++++ checking ClusterHitsOK with jht "<<plane<<":"<<fHits[jht].WireID().Wire<<":"<<(int)fHits[jht].PeakTime()<<" RMS "<<fHits[jht].RMS()<<" Multiplicity "<<fHits[jht].Multiplicity();
          // Ensure that the hits StartTick and EndTick have the proper overlap
          // stuff these into fcl2hits for ClusterHitsOK
          fcl2hits.resize(2);
          fcl2hits[0] = iht; fcl2hits[1] = jht;
          if(!ClusterHitsOK(-1)) continue;
          // Ensure that we pick up all the hits in a iht multiplet.
          // We can do this by skipping this iht if we would be stepping away from
          // the other hits.
          if(fHits[iht].Multiplicity() > 1) {
            // Will move in -time direction, away from the hits in the +time direction from iht
            if(fHits[iht].Multiplicity() > 2 && fHits[iht].LocalIndex() < (fHits[iht].Multiplicity()-1) && fHits[iht].PeakTime() > fHits[jht].PeakTime()) continue;
            // The reverse case **shouldn't** be needed because the hit loop iterates from -time to +time
            // Also check to ensure that we don't start going in a non-physical direction when
            // the other hits in the multiplet are used. Do the simple case of a doublet
            if(fHits[iht].Multiplicity() == 2) {
              if(fHits[iht].LocalIndex() == 0 && fHits[jht].PeakTime() > fHits[iht].PeakTime()) continue;
              if(fHits[iht].LocalIndex() == 1 && fHits[jht].PeakTime() < fHits[iht].PeakTime()) continue;
            }
          } // fHits[iht].Multiplicity() > 1
          if(prt) mf::LogVerbatim("CC")<<"  Starting trajectory";
          // start a trajectory in the direction from iht -> jht
          toWire = jwire;
          if(fHits[jht].Multiplicity() == 1) {
            toTick = fHits[jht].PeakTime();
          } else {
            HitMultipletPosition(jht, toTick);
          }
          StartTraj(iht, toWire, toTick);
          // check for a failure
          if(work.Pts.size() == 0) {
            if(prt) mf::LogVerbatim("CC")<<"ClusterLoop2: StartTraj failed";
            continue;
          }
          
          if(prt) PrintWork(USHRT_MAX);
          StepCrawl(success);
          if(!success) {
            if(prt) mf::LogVerbatim("CC")<<" xxxxxxx StepCrawl failed ";
            ReleaseAllWorkHits();
            continue;
          }
/*
          // decide if we should reverse and keep going on large angle trajectories
          if(std::abs(work.Pts[work.Pts.size()-1].Ang) > 1.2) {
            // re-tag the traj hits
            if(prt) mf::LogVerbatim("CC")<<"TRP  Large angle traj reverse trajectory and keep stepping";
            FlagAllWorkHits();
            ReverseTraj();
            if(prt) PrintWork(USHRT_MAX);
            StepCrawl(success);
          }
*/
          if(prt) mf::LogVerbatim("CC")<<"StepCrawl done: work.Pts size "<<work.Pts.size();
//          if(work.Pts.size() < 2) continue;
          if(work.Pts.size() < 3) continue;
          // temp
          if(fStepCrawlStudyMode && work.Pts.size() > 6 && work.CrawlDir > 0) {
            unsigned short lastPt = work.Pts.size() - 1;
            float ang = std::abs(work.Pts[lastPt].Ang);
            if(ang > M_PI/2) ang = M_PI - ang;
            mf::LogVerbatim("CC")<<"ANG "<<plane<<" "<<ang<<" "<<work.Pts[lastPt].AveHitRes<<" "<<TrajLength(work);
          }
          StoreTraj();
          break;
        } // jht
      } // iht
    } // iwire
    
    prt = false;

    FillTrajTruth();
    if((int)plane == fDebugPlane) PrintAllTraj(USHRT_MAX, 0);
    
    work.Pts.clear();
    
    // need to fake out FindVertices
    pass = 1;
    if(fFindTrajVertices) FindVertices();
    if(fTagClusters) TagClusters();
    
    allTraj.clear();

    CheckHitClusterAssociations();
    
  } // ClusterLoop2
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::FillTrajTruth()
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
//      mf::LogVerbatim("CC")<<"TrackID "<<trackID<<" pdg "<<part->PdgCode()<<" E "<<part->E()<<" mass "<<part->Mass()<<" KE "<<KE<<" Mother "<<part->Mother()<<" Proc "<<part->Process();
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
//      mf::LogVerbatim("CC")<<"Traj "<<itj<<" match with trackID "<<trackID;
      allTraj[itj].TruPDG = plist2[trackID]->PdgCode();
      allTraj[itj].TruKE = 1000 * (plist2[trackID]->E() - plist2[trackID]->Mass());
      if(plist2[trackID]->Process() == "primary") allTraj[itj].IsPrimary = true;
      if(plist2[trackID]->Mother() == 1) allTraj[itj].IsSecondary = true;
    } // itj

    
  } // FillTrajTruth
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::TagClusters()
  {
    // try to tag as shower-like or track-like
    unsigned short ii, itjPt, end, jj, jtjPt, iClosePt;
    
    if(allTraj.size() < 2) return;
    // make sure trajectories are in CC order
    if(allTraj[0].Pts[0].Pos[0] < allTraj[0].Pts[1].Pos[0]) return;
    
//    if(plane != 1) return;
    
    mf::LogVerbatim("CC")<<"Inside TagClusters: plane "<<plane;
    
    PrintAllTraj(USHRT_MAX, 0);
    
    struct PhotonConv {
      unsigned short LeadTraj;    // index of the lead trajectory (itj)
      unsigned short LeadClosePt; // index of the DOCA on the LeadTraj trajectory
      unsigned short Traj2;       // index of the second trajectory (jtj)
      float MinDist;              // Transverse separation distance at LeadClosePt
      float SepDist;              // Distance along itj where the trajectories separate
      float LeadDAng;             // Angle between the end of itj and the end of jtj
    };
    std::vector<PhotonConv> pcCand;
    
    float dx, dy, dang, dist, sep;
    
    for(ii = 0; ii < allTraj.size(); ++ii) {
      // Try to characterize the itj trajectory as a gamma conversion.
      Trajectory& itj = allTraj[ii];
      if(itj.ProcCode == USHRT_MAX) continue;
      if(itj.Pts.size() < 10) continue;
      itjPt = itj.Pts.size() - 1;
      // temp
//      if(ii > 0) continue;
      // temp
      if(!itj.IsSecondary) continue;
      // ensure that there are no hits near the End of itj
      // as should be expected for a photon converion
      for(jj = 0; jj < allTraj.size(); ++jj) {
        if(jj == ii) continue;
        Trajectory& jtj = allTraj[jj];
        if(jtj.ProcCode == USHRT_MAX) continue;
        if(jtj.Pts.size() < 5) continue;
        for(end = 0; end < 2; ++end) {
          // temp
          if(end == 0) continue;
          if(end == 0) {
            // require the jtj End to be US of the itj End
            if(jtj.Pts[0].Pos[0] > itj.Pts[0].Pos[0]) continue;
            itjPt = 0;
            jtjPt = 0;
          }
          else {
            itjPt = itj.Pts.size() - 1;
            jtjPt = jtj.Pts.size() - 1;
            // require the jtj End to be DS of the itj End
            if(jtj.Pts[jtjPt].Pos[0] < itj.Pts[itjPt].Pos[0]) continue;
          }
          // project jtj to the End of itj
          // index of the other end
          TrajClosestApproach(itj, jtj.Pts[jtjPt].Pos[0], jtj.Pts[jtjPt].Pos[1], iClosePt, dist);
          dang = std::abs(itj.Pts[itjPt].Ang - jtj.Pts[jtjPt].Ang);
          dx = itj.Pts[itjPt].Pos[0] - jtj.Pts[jtjPt].Pos[0];
          dy = itj.Pts[itjPt].Pos[1] - jtj.Pts[jtjPt].Pos[1];
          sep = sqrt(dx * dx + dy * dy);
          mf::LogVerbatim("CC")<<"PCC "<<ii<<" "<<jj<<" "<<dist<<" "<<dang<<" "<<sep<<" "<<(int)TrajLength(itj)<<" "<<(int)TrajLength(jtj);
        } // end
      } // jj
    } // ii
    
    mf::LogVerbatim("CC")<<"Photon conversion candidates";
    for(unsigned short ipc = 0; ipc < pcCand.size(); ++ipc) {
      mf::LogVerbatim("CC")<<"ipc "<<ipc<<" Lead "<<pcCand[ipc].LeadTraj<<" Second traj "<<pcCand[ipc].Traj2;
    } // ipc
  } // TagClusters
  
  //////////////////////////////////////////
  float ClusterCrawlerAlg::TrajLength(Trajectory& tj)
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
  float ClusterCrawlerAlg::TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
  {
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation
 
  //////////////////////////////////////////
  void ClusterCrawlerAlg::FindTrajVertices()
  {
    // Looks for vertices that ClusterCrawler can't find.
    // All trajectories **should_be** in ClusterCrawler order

    if(allTraj.size() < 2) return;

    vtxprt = (fDebugPlane == (int)plane && fDebugHit < 0);

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
            aVtx.NClusters = 2; aVtx.CTP = clCTP; aVtx.Fixed = false;
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
              mf::LogVerbatim("CC")<<"Vtx itj "<<itj<<" clsID "<<allTraj[itj].ClusterIndex+1<<" itjPt "<<itjPt<<" jtj "<<jtj<<" clsID "<<allTraj[jtj].ClusterIndex+1<<" jtjPt "<<jtjPt<<" Vtx "<<aVtx.Wire<<" "<<aVtx.Time;
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
                  if(vtxprt) mf::LogVerbatim("CC")<<" add ktj"<<ktj<<" clsID "<<clsIndex+1<<" Vtx "<<vtx[ivx].Wire<<" "<<vtx[ivx].Time;
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

  //////////////////////////////////////////
  unsigned short ClusterCrawlerAlg::VtxClusterEnd(unsigned short ivx, unsigned short icl)
  {
    // returns 0 if the vertex is closest to the cluster Begin and 1 if the vertex is closer
    // to the End
    
    if(ivx > vtx.size() - 1) return 0;
    if(icl > tcl.size() - 1) return 0;
    
    float dw, dt, dr0, dr1;
    
    dw = vtx[ivx].Wire - tcl[icl].BeginWir;
    dt = vtx[ivx].Time - tcl[icl].BeginTim;
    dr0 = dw * dw + dt * dt;
    
    dw = vtx[ivx].Wire - tcl[icl].EndWir;
    dt = vtx[ivx].Time - tcl[icl].EndTim;
    dr1 = dw * dw + dt * dt;
    
    if(dr0 < dr1) { return 0; } else { return 1; }
    
  } // VtxClusterEnd

  //////////////////////////////////////////
  void ClusterCrawlerAlg::TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& ipt, float& Distance)
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
  
  /*
  //////////////////////////////////////////
  void ClusterCrawlerAlg::FindTrajVertices()
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

    if(vtxprt) mf::LogVerbatim("CC")<<"FindTrajVertices: intPts size "<<intPts.size();
    
    if(intPts.size() == 0) return;
    
    if(vtxprt) {
      for(ii = 0; ii < intPts.size(); ++ii) {
        itj = intPts[ii].iTraj;
        jtj = intPts[ii].iTraj;
        mf::LogVerbatim("CC")<<ii<<" itj "<<itj<<" cls "<<allTraj[itj].ClusterIndex<<" itjPt "<<intPts[ii].iTrajPt<<" jtj "<<jtj<<" cls "<<allTraj[jtj].ClusterIndex<<" jPt "<<intPts[ii].jTrajPt<<" DOCA "<<intPts[ii].DOCA;
        
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
        mf::LogWarning("CC")<<"Cluster "<<tcl[jcl].ID<<" has already been split ";
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
        if(vtxprt) mf::LogVerbatim("CC")<<" kink angle at split point "<<dang;
        if(dang > fStepCrawlKinkAngCut) {
          prt = vtxprt;
          if(!SplitAllTraj(jtj, jtjPt, ivx)) {
            mf::LogWarning("CC")<<"FindTrajVertices: SplitAllTraj failed";
            return;
          }
          // fit the vertex position
          vtx[ivx].Fixed = false;
          FitVtx(ivx);
        } // dang > fStepCrawlKinkAngCut
      } // split point in the middle of jtj
    } // ii
    
  } // FindTrajVertices
*/
  //////////////////////////////////////////
  bool ClusterCrawlerAlg::SplitAllTraj(unsigned short itj, unsigned short pos, unsigned short ivx)
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
      mf::LogVerbatim("CC")<<"*************************** Inside SplitAllTraj ***************************";
      PrintAllTraj(itj, USHRT_MAX);
    }
    
    unsigned short icl = allTraj[itj].ClusterIndex;
    MakeClusterObsolete(icl);
    // double check hit releasing
    for(unsigned short ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
      unsigned int hit = allTraj[itj].Pts[ipt].Hits[0];
      if(inClus[hit] != 0) {
        mf::LogWarning("CC")<<"SplitAllTraj: MakeClusterObsolete problem on icl "<<icl<<" itj "<<itj<<" inClus "<<inClus[hit]<<" hit "<<plane<<":"<<fHits[hit].WireID().Wire<<":"<<(int)fHits[hit].PeakTime();
        return false;
      }
    } // check
    
    // flag the old traj as obsolete
    allTraj[itj].ProcCode = USHRT_MAX;
    
    work.Pts.clear();
    work.Pts.insert(work.Pts.end(), allTraj[itj].Pts.begin(), allTraj[itj].Pts.begin()+pos);
    work.ProcCode = 911;
    StoreTraj();
    if(prt) {
      mf::LogVerbatim("CC")<<"Check DS end";
      PrintWork(USHRT_MAX);
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
      mf::LogVerbatim("CC")<<"Check US end";
      PrintWork(USHRT_MAX);
    }
    StoreTraj();
    if(ivx < vtx.size()) {
      tcl[tcl.size()-1].EndVtx = ivx;
      work.Vtx[0] = ivx;
    }
    return true;
  } // SplitAllTraj
  
  //////////////////////////////////////////
  float ClusterCrawlerAlg::PointTrajDOCA(float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach between a (wire, tick) point
    // and a trajectory point
    
    // find the angle of a line between the point and the trajectory point
    float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
    float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    float dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return sqrt(dw * dw + dt * dt);

  } // PointTrajDOCA
  
  //////////////////////////////////////////
  bool ClusterCrawlerAlg::TrajIntersection(TrajPoint& tp1, TrajPoint tp2, float& x, float& y)
  {
    // returns the intersection position, intPos, of two trajectory point and
    // returns true if successfull
    
    x = -9999; y = -9999;
    
    double arg1 = tp1.Pos[0] * tp1.Dir[1] - tp1.Pos[1] * tp1.Dir[0];
    double arg2 = tp2.Pos[0] * tp1.Dir[1] - tp2.Pos[1] * tp1.Dir[0];
    double arg3 = tp2.Dir[0] * tp1.Dir[1] - tp2.Dir[1] * tp1.Dir[0];
    if(arg3 == 0) return false;
    double s = (arg1 - arg2) / arg3;
    
    x = (float)(tp1.Pos[0] + s * tp1.Dir[0]);
    y = (float)(tp1.Pos[1] + s * tp1.Dir[1]);
    return true;
    
  } // TrajIntersection
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::StepCrawl(bool& success)
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
    
    unsigned int nit;
    float bigChgDiffCut = 3 * fStepCrawlChgDiffCut;
    float maxPos0 = fNumWires;
    float maxPos1 = fMaxTime * fScaleF;
    
    // count number of steps taken with no trajectory point added
    unsigned short nMissedStep = 0;
    unsigned short missedStepCut = SetMissedStepCut(tp);

    bool SignalPresent, stopOnKink = false;
    unsigned short killPts;
    for(nit = 0; nit < 10000; ++nit) {
      // move the position by one step in the right direction
      for(iwt = 0; iwt < 2; ++iwt) tp.Pos[iwt] += tp.Dir[iwt];
      // hit the boundary of the TPC?
      if(tp.Pos[0] < 0 || tp.Pos[0] > maxPos0) break;
      if(tp.Pos[1] < 0 || tp.Pos[1] > maxPos1) break;
      // remove the old hits
      tp.Hits.clear();
      // look for new hits
      AddTrajHits(tp, SignalPresent);
      if(!SignalPresent) {
        // cut on the number of missed wires
        lastPt = work.Pts.size() - 1;
        if(prt)  mf::LogVerbatim("CC")<<" No signal. Missed wires "<<std::abs(tp.Pos[0] - work.Pts[lastPt].Pos[0])<<" user cut "<<fStepCrawlMaxWireSkip;
        if(std::abs(tp.Pos[0] - work.Pts[lastPt].Pos[0]) < fStepCrawlMaxWireSkip) continue;
        break;
      }
      if(tp.Hits.size() == 0) {
        ++nMissedStep;
        // Break if we took too many steps without adding a traj point
        if(prt) mf::LogVerbatim("CC")<<" nMissedStep "<<nMissedStep<<" missedStepCut "<<missedStepCut<<" Dir "<<tp.Dir[0]<<" "<<tp.Dir[1];
        if(nMissedStep > missedStepCut) {
          if(prt) mf::LogVerbatim("CC")<<" Too many steps w/o traj point ";
          break;
        }
        // Otherwise keep stepping
        continue;
      } // tp.Hits.size() == 0
      // Have new traj hits. Add the trajectory point and update
      tp.Step = nit;
      nMissedStep = 0;
      missedStepCut = SetMissedStepCut(tp);
      work.Pts.push_back(tp);
      UpdateTraj(success);
      if(!success) break;
      lastPt = work.Pts.size()-1;
      // IMPORTANT NOTE: Local tp hasn't been updated yet using results from UpdateTraj
      // Quit if we are starting out poorly
      if(work.Pts.size() == 3 && work.Pts[2].AveHitRes > 3) {
        success = false;
        if(prt) mf::LogVerbatim("CC")<<" Bad AveHitRes. Drop trajectory.";
        return;
      }
      // assume that we aren't going to kill the point we just added, or any
      // of the previous points...
      killPts = 0;
      // Look for two consecutive large Deltas
      if(lastPt > 3 && tp.Delta > fStepCrawlMaxDeltaJump * work.Pts[lastPt-2].Delta) {
        // check the previous Delta ratio
        if(work.Pts[lastPt-1].Delta > fStepCrawlMaxDeltaJump * work.Pts[lastPt-2].Delta) killPts = 2;
        // check for not nearby hits
        if(prt) mf::LogVerbatim("CC")<<"   NumNotNear "<<tp.NumNotNear<<" previous "<<work.Pts[lastPt-1].NumNotNear;
        if(tp.NumNotNear < 2 && work.Pts[lastPt-1].NumNotNear < 2) killPts = 0;
        if(prt) mf::LogVerbatim("CC")<<"   large Deltas? killPts "<<killPts;
      }
      // check the charge ratio
      if(killPts == 0 && work.Pts.size() > 10 && fChgRMS > 0) {
        // Kill this point if it has a high charge and the previous did as well
        if(tp.ChgDiff > fStepCrawlChgDiffCut && work.Pts[lastPt].ChgDiff > fStepCrawlChgDiffCut) killPts = 1;
        // or if it has extraordinarily high charge
        if(tp.ChgDiff > bigChgDiffCut) killPts = 1;
        if(prt) mf::LogVerbatim("CC")<<"   large ChgDiff? "<<tp.ChgDiff<<" killPts "<<killPts;
      } // fChgRMS > 0
      // check for a kink. Stop crawling if one is found
      if(GottaKink()) {
        stopOnKink = true;
        break;
      }
      // update the local tp unless we have killing to do
      if(killPts == 0) {
        tp = work.Pts[lastPt];
        if(prt) PrintWork(lastPt);
      } else {
        // shorten the trajectory by the desired number of points
        // and correct the projected position
        // release the hits
        ReleaseAllWorkHits();
        float nSteps = (float)(work.Pts[lastPt].Step - work.Pts[lastPt-killPts].Step);
        if(prt) {
          mf::LogVerbatim("CC")<<"   killing "<<killPts<<" points with "<<nSteps<<" steps";
          if(prt) mf::LogVerbatim("CC")<<"  Current tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1];
        }
        // truncate the work trajectory
        work.Pts.resize(lastPt+1-killPts);
        // re-assign the hits
        FlagAllWorkHits();
/*
        for(ipt = 0; ipt < work.Pts.size(); ++ipt) {
          for(iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = -3;
        } // ipt
*/
        // copy to the local trajectory point
        lastPt = work.Pts.size() - 1;
        tp = work.Pts[lastPt];
        // move the position
        tp.Pos[0] += nSteps * tp.Dir[0];
        tp.Pos[1] += nSteps * tp.Dir[1];
        nSteps = 0;
        if(prt) mf::LogVerbatim("CC")<<"  New tp.Pos     "<<tp.Pos[0]<<" "<<tp.Pos[1];
      }
    } // nit

    // check the cluster quality and possibly truncate it
    if(!stopOnKink) StepCrawlClusterCheck();

//    if(prt) PrintWork(USHRT_MAX);

    ReleaseAllWorkHits();
    
    // temp checking
//    if(prt) {
      unsigned int firsthit = (unsigned int)WireHitRange[fFirstWire].first;
      unsigned int lasthit = (unsigned int)WireHitRange[fLastWire-1].second;
      bool itsBad = false;
      bool first = true;
      mf::LogVerbatim myprt("CC");
      for(unsigned int iht = firsthit; iht < lasthit; ++iht) {
        if(inClus[iht] == -3) {
          if(first) {
            myprt<<"Not released hits ";
            first = false;
          }
          myprt<<" "<<plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime();
          inClus[iht] = 0;
          itsBad = true;
        }
      } // iht
    if(itsBad) {
      unsigned short iht = work.Pts[0].Hits[0];
      myprt<<" Seed hit "<<plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime();
    }
//    } // prt
    
    success = true;
    return;
    
  } // StepCrawl

  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::AddTrajHits(TrajPoint& tp, bool& SignalPresent)
  {
    std::vector<unsigned int> closeHits, nearbyHits;
    unsigned int wire, loWire, hiWire, iht, firstHit, lastHit;
    unsigned int lastPt = work.Pts.size() - 1;
    
    // figure out which wires to consider
    if(std::abs(tp.Dir[0]) > 0.1) {
      // smallish angle trajectory
      loWire = tp.Pos[0];
      hiWire = loWire + 1;
    } else {
      // look at adjacent wires for larger angle trajectories
      loWire = tp.Pos[0] - 1;
      hiWire = loWire + 2;
    }

    float fwire, ftime, delta;
    // find the projection error to this point
    float dw = tp.Pos[0] - work.Pts[lastPt].Pos[0];
    float dt = tp.Pos[1] - work.Pts[lastPt].Pos[1];
    float dpos = sqrt(dw * dw + dt * dt);
    float projErr = dpos * work.Pts[lastPt].AngErr;
    // Add this to the Delta RMS factor and construct a cut
    float deltaCut = 3 * (projErr + tp.DeltaRMS);
    // don't let it be too small (< 10% of a wire spacing)
    if(deltaCut < 0.1) deltaCut = 0.1;
    // cut for nearby hits
    float bigDeltaCut = 2 * deltaCut;
    // make it big for study mode
    if(fStepCrawlStudyMode) deltaCut = 1.5;
    // projected time in ticks for testing the existence of a hit signal
    raw::TDCtick_t rawProjTick = (float)(tp.Pos[1] / fScaleF);
    
    // assume failure
    SignalPresent = false;
    if(prt) mf::LogVerbatim("CC")<<" AddTrajHits: loWire "<<loWire<<" tp.Pos[0] "<<tp.Pos[0]<<" hiWire "<<hiWire<<" projTick "<<rawProjTick;
    
    for(wire = loWire; wire < hiWire; ++wire) {
      // Assume a signal exists on a dead wire
      if(WireHitRange[wire].first == -1) SignalPresent = true;
      if(WireHitRange[wire].first < 0) continue;
      firstHit = (unsigned int)WireHitRange[wire].first;
      lastHit = (unsigned int)WireHitRange[wire].second;
      fwire = wire;
      for(iht = firstHit; iht < lastHit; ++iht) {
        if(rawProjTick > fHits[iht].StartTick() && rawProjTick < fHits[iht].EndTick()) SignalPresent = true;
        ftime = fScaleF * fHits[iht].PeakTime();
        delta = PointTrajDOCA(fwire, ftime, tp);
        if(prt) {
          bool close = delta < deltaCut;
          bool nearby = delta < bigDeltaCut;
          mf::LogVerbatim("CC")<<"  chk "<<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" delta "<<delta<<" deltaCut "<<deltaCut<<" inClus "<<inClus[iht]<<" Chg "<<(int)fHits[iht].Integral()<<" SignalPresent "<<SignalPresent<<" close "<<close<<" nearby "<<nearby;
        }
        if(inClus[iht] != 0) continue;
        if(std::abs(delta) > deltaCut) {
          // missed the small window but inside the larger one
          if(std::abs(delta) < bigDeltaCut) {
            SignalPresent = true;
            nearbyHits.push_back(iht);
            inClus[iht] = -3;
          }
          continue;
        }
        closeHits.push_back(iht);
        inClus[iht] = -3;
        SignalPresent = true;
      } // iht
    } // wire
    if(!SignalPresent) {
      if(prt) mf::LogVerbatim("CC")<<" no signal on any wire at tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" tick "<<(int)tp.Pos[1]/fScaleF<<" closeHits size "<<closeHits.size();
      return;
    }
    if(closeHits.size() == 0 && nearbyHits.size() == 0) return;
    // don't use nearbyHits when the trajectory is short and not yet well defined
    if(closeHits.size() > 0) {
      // Add the hits to the traj
      tp.Hits = closeHits;
      if(prt) mf::LogVerbatim("CC")<<" Added "<<closeHits.size()<<" closeHits ";
    } else {
      // no close hits found but other hits are nearby.
      // Don't use nearbyHits when the trajectory is short and not yet well defined
      if(work.Pts.size() < 4) return;
      tp.Hits = nearbyHits;
      tp.UsedNotCloseHit = 0;
      if(prt) mf::LogVerbatim("CC")<<" Added "<<closeHits.size()<<" nearbyHits ";
    }
    tp.NumNotNear = nearbyHits.size();
    // specify the hits position
    tp.HitPos[0] = 0;
    tp.HitPos[1] = 0;
    tp.Chg = 0;
    if(tp.Hits.size() > 1) {
      // sort the new hits by distance from the previous traj point
      std::vector<SortEntry> sortVec;
      SortEntry sortEntry;
      unsigned short ii;
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        sortEntry.index = ii;
        iht = tp.Hits[ii];
        dw = fHits[iht].WireID().Wire - work.Pts[lastPt].Pos[0];
        dt = fHits[iht].PeakTime() * fScaleF - work.Pts[lastPt].Pos[1];
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
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      iht = tp.Hits[ii];
      tp.Chg += fHits[iht].Integral();
      tp.HitPos[0] += fHits[iht].Integral() * fHits[iht].WireID().Wire;
      tp.HitPos[1] += fHits[iht].Integral() * fHits[iht].PeakTime() * fScaleF;
    }
    tp.HitPos[0] /= tp.Chg;
    tp.HitPos[1] /= tp.Chg;
//    if(prt) mf::LogVerbatim("CC")<<" HitPos "<<tp.HitPos[0]<<" tick "<<(int)tp.HitPos[1]/fScaleF<<" nhits "<<closeHits.size();

  } // AddTrajHits

/*
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::AddTrajHits(TrajPoint& tp, bool& stopCrawl)
  {
    unsigned int loWire, hiWire, wire, firstHit, lastHit, iht;
    unsigned int loloWire, hihiWire;
    // loTime and hiTime are in Wire Space Equivalent (WSE) units
    float loTime, hiTime, loloTime, hihiTime;
    // Start at the end of the trajectory
    // pos is in wire spacing equivalent units
    std::array<float, 2> hitsPos;
    float sum;
    std::vector<unsigned int> closeHits;
    unsigned short nMissedSignal, nSignal;
    // index of a (single) nearby hit
    unsigned int nearbyHitIndex;
    float projTick, hitDir;
    // version for testing StartTick and EndTick
    raw::TDCtick_t rawProjTick;
    nMissedSignal = 0;
    // count number of steps taken with no trajectory point added
    bool largeAngle;
    // add hits to the trajectory point. Assume success
    stopCrawl = false;
    // assume that we will find a hit in the smaller window
    tp.UsedNotCloseHit = false;
    tp.NumNotNear = 0;
    largeAngle = (std::abs(tp.Dir[1]) > 0.7);
    GetStepCrawlWindow(tp, loWire, hiWire, loTime, hiTime, loloWire, hihiWire, loloTime, hihiTime);
    hitsPos[0] = 0; hitsPos[1] = 0; sum = 0;
    nSignal = 0;
    // projected time in tick units
    projTick = tp.Pos[1] / fScaleF;
    rawProjTick = projTick;
    if(prt) mf::LogVerbatim("CC")<<" loWire "<<loWire<<" tp.Pos[0] "<<tp.Pos[0]<<" hiWire "<<hiWire<<" loTime "<<(int)loTime<<" projTick "<<(int)projTick<<" hiTime "<<(int)hiTime;
    for(wire = loWire; wire < hiWire; ++wire) {
      // Assume a signal exists on a dead wire
      if(WireHitRange[wire].first == -1) ++nSignal;
      if(WireHitRange[wire].first < 0) continue;
      firstHit = (unsigned int)WireHitRange[wire].first;
      lastHit = (unsigned int)WireHitRange[wire].second;
      for(iht = firstHit; iht < lastHit; ++iht) {
        // count number of nearby hit signals
        if(rawProjTick > fHits[iht].StartTick() && rawProjTick < fHits[iht].EndTick()) ++nSignal;
        if(prt) mf::LogVerbatim("CC")<<"  chk "<<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" StartTick "<<fHits[iht].StartTick()<<" EndTick "<<fHits[iht].EndTick()<<" inClus "<<inClus[iht]<<" Chg "<<(int)fHits[iht].Integral()<<" nSignal "<<nSignal;
        if(fHits[iht].PeakTime() < loTime) continue;
        if(fHits[iht].PeakTime() > hiTime) break;
        if(inClus[iht] != 0) continue;
        // ensure that the hit is in the appropriate time direction
        if(largeAngle) {
          hitDir = fHits[iht].PeakTime() - work.Pts[work.Pts.size()-1].Pos[1];
          if(prt) mf::LogVerbatim("CC")<<" check hitDir "<<hitDir;
          if(hitDir > 0 && work.Pts[work.Pts.size()-1].Dir[1] < 0) continue;
          if(hitDir < 0 && work.Pts[work.Pts.size()-1].Dir[1] > 0) continue;
        }
        if(prt) mf::LogVerbatim("CC")<<"   ADD "<<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime();
        closeHits.push_back(iht);
        inClus[iht] = -3;
        // prepare to find the new position in WSE units
        sum += fHits[iht].Integral();
        hitsPos[0] += fHits[iht].Integral() * fHits[iht].WireID().Wire;
        hitsPos[1] += fHits[iht].Integral() * fHits[iht].PeakTime() * fScaleF;
      } // iht
    } // wire
    if(prt) mf::LogVerbatim("CC")<<" Done looking for hits in small window. closeHits size "<<closeHits.size();
    // didn't find a hit within the small window. Look for a hit in the larger window
    if(closeHits.size() == 0) {
      // Check for a hit in a wider window
      // Require nearby hit to be within (2 WSE)^2
      float best = 1000, dpos;
      nearbyHitIndex = UINT_MAX;
      tp.NumNotNear = 0;
      for(wire = loloWire; wire < hihiWire; ++wire) {
        if(WireHitRange[wire].first == -1) ++nSignal;
        if(WireHitRange[wire].first < 0) continue;
        firstHit = (unsigned int)WireHitRange[wire].first;
        lastHit = (unsigned int)WireHitRange[wire].second;
        for(iht = firstHit; iht < lastHit; ++iht) {
//            if(prt) mf::LogVerbatim("CC")<<"  chk nearby "<<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" StartTick "<<fHits[iht].StartTick()<<" EndTick "<<fHits[iht].EndTick()<<" inClus "<<inClus[iht]<<" Chg "<<(int)fHits[iht].Integral();
          if(inClus[iht] != 0) continue;
          if(fHits[iht].PeakTime() < loloTime) continue;
          if(fHits[iht].PeakTime() > hihiTime) continue;
          ++tp.NumNotNear;
          // hijack hitsPos to find the separation
          hitsPos[0] = fHits[iht].WireID().Wire - tp.Pos[0];
          hitsPos[1] = fHits[iht].PeakTime() * fScaleF - tp.Pos[1];
          dpos = hitsPos[0] * hitsPos[0] + hitsPos[1] * hitsPos[1];
          if(dpos < best) {
            best = dpos;
            nearbyHitIndex = iht;
            if(prt) mf::LogVerbatim("CC")<<"  Found nearby hit "
              <<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" StartTick "<<fHits[iht].StartTick()<<" EndTick "<<fHits[iht].EndTick()<<" inClus "<<inClus[iht]<<" Integral "<<(int)fHits[iht].Integral()<<" dpos "<<dpos;
          }
        } // iht
      } // wire
      if(prt) mf::LogVerbatim("CC")<<" Done looking for hits in large window. nearbyHitIndex "<<nearbyHitIndex;
      if(nearbyHitIndex < fHits.size()) {
        closeHits.push_back(nearbyHitIndex);
        inClus[nearbyHitIndex] = -3;
        // set the flag indicating this is a hit outside the smaller window
        tp.UsedNotCloseHit = true;
        // prepare to find the new position in WSE units
        sum += fHits[nearbyHitIndex].Integral();
        hitsPos[0] += fHits[nearbyHitIndex].Integral() * fHits[nearbyHitIndex].WireID().Wire;
        hitsPos[1] += fHits[nearbyHitIndex].Integral() * fHits[nearbyHitIndex].PeakTime() * fScaleF;
        nMissedSignal = 0;
      } // nearbyHitIndex < fHits.size()
      else {
        ++nMissedSignal;
        if(nMissedSignal > 4) {
          if(prt) mf::LogVerbatim("CC")<<" no signal on any wire at tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" tick "<<(int)tp.Pos[1]/fScaleF;
          // release the hits
          for(iht = 0; iht < closeHits.size(); ++iht) inClus[closeHits[iht]] = 0;
          stopCrawl = true;
          return;
        }
      } // nNearby == 0
    } // closeHits.size() == 0
    // No hit found in the larger window
    if(closeHits.size() == 0) return;
    // stuff the charge weighted position into tp and try to append it to the trajectory
    for(unsigned short iwt = 0; iwt < 2; ++iwt) tp.HitPos[iwt] = hitsPos[iwt] / sum;
    tp.Chg = sum;
    if(closeHits.size() > 1) {
      // sort the new hits by distance from the previous traj point
      std::vector<SortEntry> sortVec;
      SortEntry sortEntry;
      unsigned short ii;
      bool sortByTime = (std::abs(tp.Dir[1]) > std::abs(tp.Dir[0]));
      for(ii = 0; ii < closeHits.size(); ++ii) {
        sortEntry.index = ii;
        // xxx
        if(sortByTime) {
          // sort by time
          sortEntry.length = fHits[closeHits[ii]].PeakTime();
        } else {
          // sort by wire
          sortEntry.length = (float)fHits[closeHits[ii]].WireID().Wire;
        } // !SortByTime
        sortVec.push_back(sortEntry);
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), lessThan);
      // make a temp vector
      std::vector<unsigned int> tmp = closeHits;
      // overwrite with the sorted values
      for(ii = 0; ii < sortVec.size(); ++ii) tmp[ii] = closeHits[sortVec[ii].index];
      // replace
      closeHits = tmp;
      // swap the order? TODO: check
      if(sortByTime && tp.Dir[1] < 0) std::reverse(closeHits.begin(), closeHits.end());
      if(!sortByTime && tp.Dir[0] < 0) std::reverse(closeHits.begin(), closeHits.end());
    }
    // Now have a good trajectory point. Add the hits to the traj
    tp.Hits = closeHits;
    stopCrawl = false;
  } // AddTrajHits
*/
  ////////////////////////////////////////////////
  bool ClusterCrawlerAlg::GottaKink()
  {
    // Does a local fit of just-added traj points to identify a kink while step crawling.
    // Truncates the traj vector and returns true if one is found.
    
    if(work.Pts.size() < 8) return false;
    
    // The kink angle cut will be scaled by the normalized average hit resolution so this
    // must be known before proceeding
    if(work.Pts[work.Pts.size()-1].AveHitRes < 0) return false;
    
    unsigned short nTPF = 4;
    
    unsigned short fromIndex = work.Pts.size() - nTPF;
    unsigned short toIndex = work.Pts.size() - 1;
    TrajPoint tp1;
    FitTrajMid(fromIndex, toIndex, tp1);
    if(prt) mf::LogVerbatim("CC")<<"Fit1 "<<fromIndex<<" "<<toIndex<<" Ang "<<tp1.Ang<<" chi "<<tp1.FitChi;
    if(tp1.FitChi > 900) {
      unsigned int hit = tp1.Hits[0];
      mf::LogWarning("CC")<<"GottaKink: Bad traj fit1. Seed hit "<<plane<<":"<<fHits[hit].WireID().Wire<<":"<<(int)fHits[hit].PeakTime()<<" Pts size "<<work.Pts.size()<<" fromIndex "<<fromIndex<<" to "<<toIndex;
      PrintWork(USHRT_MAX);
      return false;
    }
    fromIndex -= nTPF;
    toIndex -= nTPF;
    TrajPoint tp2;
    FitTrajMid(fromIndex, toIndex, tp2);
    if(tp2.FitChi > 900) {
      unsigned int hit = tp2.Hits[0];
      mf::LogWarning("CC")<<"GottaKink: Bad traj fit2 "<<fromIndex<<" "<<toIndex<<" work.Pts size "<<work.Pts.size()<<" Seed hit "<<plane<<":"<<fHits[hit].WireID().Wire<<":"<<(int)fHits[hit].PeakTime();
      return false;
    }
    float dang = std::abs(tp1.Ang - tp2.Ang);
    float cut, cutBig;
    cut = fStepCrawlKinkAngCut;
    cutBig = 1.2 * fStepCrawlKinkAngCut;
    if(prt) mf::LogVerbatim("CC")<<"GottaKink: Ang1 "<<tp1.Ang<<" Ang2 "<<tp2.Ang<<" dang "<<dang<<" cut "<<cut;
    if(dang < cutBig) return false;
    
    unsigned short kinkIndex = work.Pts.size()-nTPF-1;
    // Look at bit closer if we are near the cut
    // Probably not a kink if there is no big increase in Delta
    if(work.Pts[kinkIndex].Delta < 2 * work.Pts[kinkIndex].Delta) return false;
    
    // found a kink so truncate
    // TODO Actually we should start a new trajectory using these points
    // and create a vertex between the two of them
    if(prt) mf::LogVerbatim("CC")<<"TRP kink angle "<<dang<<" cut "<<cut<<" work.Pts old size "<<work.Pts.size()<<" new size "<<kinkIndex;
    ReleaseAllWorkHits();
    work.Pts.resize(kinkIndex);
    return true;
    
  } // GottaKink
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::StepCrawlClusterCheck()
  {
    // analyze the end of a cluster after crawling and trim any
    // vagrant trajectory points and hits if they have a kink
    
    // ignore short clusters
    if(work.Pts.size() < 10) return;
    
    unsigned short lastPt = work.Pts.size()-1;
    
    if(prt) mf::LogVerbatim("CC")<<"StepCrawlClusterCheck: AveHitRes "<<work.Pts[work.Pts.size()-1].AveHitRes;
    // ignore lower quality clusters
    if(work.Pts[lastPt].AveHitRes > 1.5) return;
    
    float delrat;
    unsigned short ipt;
    bool lopped = false;
    for(ipt = work.Pts.size() - 3; ipt < work.Pts.size(); ++ipt) {
      // check for a largish change in Delta
      delrat = work.Pts[ipt].Delta / work.Pts[ipt-1].Delta;
//      if(prt) mf::LogVerbatim("CC")<<" ipt "<<ipt<<" delrat "<<delrat;
      if(delrat > 2) {
        ReleaseAllWorkHits();
        work.Pts.resize(ipt);
        if(prt) mf::LogVerbatim("CC")<<" new work.Pts size "<<work.Pts.size();
        lopped = true;
        break;
      }
    } // ipt
    
    if(!lopped) {
      if(prt) mf::LogVerbatim("CC")<<" SCCC: cluster OK ";
      return;
    }
    
    // Check for a missing wire at the end - small angle cluster
    ipt = work.Pts.size() - 1;
    if(prt) mf::LogVerbatim("CC")<<" SCCC: check last point. "<<ipt<<" Dir "<<work.Pts[ipt].Dir[1]<<" dPos[0] "<<std::abs(work.Pts[ipt].HitPos[0] - work.Pts[ipt-1].HitPos[0]);
    if(std::abs(work.Pts[ipt].Dir[1]) > 0.95) return;
    if(std::abs(work.Pts[ipt].HitPos[0] - work.Pts[ipt-1].HitPos[0]) > 1.5) work.Pts.pop_back();
    if(prt) {
      mf::LogVerbatim("CC")<<" SCCC: lopped one more work.Pts point new work.Pts size "<<work.Pts.size();
//      PrintWork(USHRT_MAX);
    }
    
  } // StepCrawlClusterCheck
  
  //////////////////////////////////////////
  unsigned short ClusterCrawlerAlg::SetMissedStepCut(TrajPoint const& tp)
  {
    unsigned short itmp = fNumWires;
    if(tp.Dir[1] != 0) itmp = 1 + 4 * std::abs(1/tp.Dir[1]);
    if(itmp > fNumWires) itmp = fNumWires;
    if(itmp < 2) itmp = 2;
    return itmp;
  }
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::UpdateTraj(bool& success)
  {
    // A new trajectory point was put on the traj vector with hits by StepCrawl.
    // This routine updates the last trajectory point fit, average hit rms, etc.
    // StepCrawl will decide what to do next with this information.
    
    success = false;
    if(work.Pts.size() < 2) return;
    
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];
    
    unsigned short ii, cnt, totCnt;
    
    // always update the charge difference
    UpdateTrajChgDiff();
    UpdateTrajDelta();

    bool smallAngle = (std::abs(lastTP.Ang) < 0.7);
    unsigned short nTPF = work.Pts[lastPt].NumTPsFit;
    if(work.Pts.size() == 2) {
      // Handle the second trajectory point. No error calculation. Just update
      // the position and direction
      lastTP.NumTPsFit = 2;
      FitTraj();
      lastTP.FitChi = 0.01;
      lastTP.AngErr = work.Pts[0].AngErr;
      if(prt) mf::LogVerbatim("CC")<<"UpdateTraj: Second traj point pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1];
      success = true;
      return;
    } else if(work.Pts.size() == 3) {
       // Third trajectory point. Keep it simple
      nTPF = 3;
      if(prt) mf::LogVerbatim("CC")<<"UpdateTraj: Third traj point "<<nTPF;
    } else {
      // many trajectory points available to fit. Determine how many
      // trajectory points should be fit to a line. Start with the number used in
      // the previous trajectory point
      if(prt) mf::LogVerbatim("CC")<<"UpdateTraj: Starting nTPF "<<nTPF;
      // decide if this should be shorter or longer
      // Count the number of previous traj points that used "nearby" hits
      cnt = 0;
      totCnt = 0;
      if(smallAngle && work.Pts.size() > 3) {
        for(ii = 0; ii < nTPF; ++ii) {
          if(!work.Pts[lastPt - ii].UsedNotCloseHit) break;
          ++cnt;
        } // ii
        // count of all nearby hits used in the trajectory fit
        for(ii = 0; ii < nTPF; ++ii) if(work.Pts[lastPt - ii].UsedNotCloseHit) ++totCnt;
      } // not a short trajectory
      // Handle the case where reducing the number of points fit doesn't work
      if(cnt == 0 && work.Pts[lastPt-2].FitChi > 0) {
        // Hit(s) just added were within the small window and there were no
        // other nearby hits added. Before adding another point to the fit
        // ensure that the fit chisq isn't increasing too much
        float chirat = 0;
        if(work.Pts.size() > 4) chirat = work.Pts[lastPt-1].FitChi / work.Pts[lastPt-2].FitChi;
        if(totCnt == 0 && chirat < 1.2) ++nTPF;
        if(prt) mf::LogVerbatim("CC")<<"  check chirat "<<chirat<<" totCnt "<<totCnt<<" new nTPF "<<nTPF;
      } else if(cnt < 3) {
        // Had a few consecutive hits outside the small window. Reduce the number of points fit
        nTPF = (unsigned short)(0.7 * (float)nTPF);
      } else {
        // Had many points outside the small window. Reduce significantly from the previous traj point
        unsigned short prevnTPF = work.Pts[lastPt - 1].NumTPsFit;
        if(prt) mf::LogVerbatim("CC")<<"  In trouble "<<prevnTPF<<" "<<lastTP.NumTPsFit<<" work.Pts size "<<work.Pts.size()<<" new "<<nTPF<<" cnt "<<cnt;
      }
    }

    // check for craziness
    if(nTPF < 2) nTPF = 2;
    lastTP.NumTPsFit = nTPF;
    FitTraj();
    // check for a failure
    if(lastTP.FitChi > 900) return;
    std::vector<unsigned int> tHits;
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt)
      tHits.insert(tHits.end(), work.Pts[ipt].Hits.begin(), work.Pts[ipt].Hits.end());
    // Calculate the average RMS width of hits on the cluster
    if(tHits.size() > 1) {
      fAveHitRMS = 0;
      for(unsigned short ii = 0; ii < tHits.size(); ++ii) fAveHitRMS += fHits[tHits[ii]].RMS();
      fAveHitRMS /= (float)tHits.size();
    }
    
    // Calculate the rms scatter of hit times
    float aveRes = 0, arg;
    float hitErr;
    // Angle correction factor found by fitting AveHitRes vs angle for protons
    // and muons. See ANG print statement in StepCrawl to write a text file.
    float ang = std::abs(work.Pts[lastPt].Ang);
    if(ang > M_PI/2) ang = M_PI - ang;
    float angCorr = 1 + 3.4 * ang * ang;
    if(fStepCrawlStudyMode) angCorr = 1;
    unsigned int hit0, hit1, hit2;
    cnt = 0;
    for(ii = 1; ii < tHits.size() - 1; ++ii) {
      hit0 = tHits[ii-1];
      hit1 = tHits[ii];
      hit2 = tHits[ii+1];
      // require hits on adjacent wires
      // Note std::abs doesn't work in the following two lines
      // TODO This doesn't work in both directions....
//      if(abs((int)fHits[hit1].WireID().Wire - (int)fHits[hit0].WireID().Wire) != 1) continue;
//      if(abs((int)fHits[hit1].WireID().Wire - (int)fHits[hit2].WireID().Wire) != 1) continue;
      arg = (fHits[hit0].PeakTime() + fHits[hit2].PeakTime())/2 - fHits[hit1].PeakTime();
      hitErr = fHitErrFac * fHits[hit1].RMS();
      arg /= hitErr;
      arg /= angCorr;
      aveRes += arg * arg;
      ++cnt;
    }
    if(cnt > 0) {
      aveRes /= (float)cnt;
      aveRes = sqrt(aveRes);
      // convert to a quality factor 1 = perfect, > 1
//      aveRes /= (fAveHitRMS * fHitErrFac);
      lastTP.AveHitRes = aveRes;
      //      if(prt) mf::LogVerbatim("CC")<<std::right<<std::setw(6)<<std::fixed<<std::setprecision(1)<<" aveRes "<<aveRes<<" fAveHitRMS "<<fAveHitRMS;
    } else {
      if(prt) mf::LogVerbatim("CC")<<"    NA tHits size "<<tHits.size()<<" cnt "<<cnt;
    }

    success = true;
    return;
    
  } // UpdateTraj

  //////////////////////////////////////////
  void ClusterCrawlerAlg::UpdateTrajDelta()
  {
    // Find Delta for the last trajectory point.
    
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];

    lastTP.Delta = PointTrajDOCA(lastTP.HitPos[0], lastTP.HitPos[1] ,lastTP);
    if(work.Pts.size() < 4) return;
    float ave = 0, sum2 = 0, fcnt = 0;
    for(unsigned short ipt = 2; ipt < work.Pts.size(); ++ipt) {
      ave += work.Pts[ipt].Delta;
      sum2 += work.Pts[ipt].Delta * work.Pts[ipt].Delta;
      ++fcnt;
    }
    ave /= fcnt;
    lastTP.DeltaRMS = sqrt((sum2 - fcnt * ave * ave) / (fcnt - 1));

  } // UpdateTrajDelta

  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::FitTrajMid(unsigned short fromIndex, unsigned short toIndex, TrajPoint& tp)
  {
    // Fit the work trajectory points from fromIndex to (and including) toIndex. The origin
    // of the fit position is set to work.Pts[fromIndex].Pos. Fit results are stashed in
    // TrajPoint tp. This copies the trajectory point from work into tp. Note that the
    // AngleErr is the error on the fitted points (not the angle error for the last point)
    
    tp.FitChi = 9999;
    
    unsigned short lastPt = work.Pts.size() - 1;
    clChisq = 999;
    if(fromIndex > lastPt || toIndex > lastPt) return;
    
    tp = work.Pts[fromIndex];
    
    unsigned short lo, hi;
    if(fromIndex < toIndex) {
      lo = fromIndex; hi = toIndex;
    } else {
      lo = toIndex; hi = fromIndex;
    }
    unsigned short nTPF = hi - lo + 1;
    tp.NumTPsFit = nTPF;
    if(nTPF < 1) return;
    ++hi;
    
    std::vector<float> x(nTPF), y(nTPF), yerr2(nTPF);
    
    // Rotate the traj hit position into the coordinate system defined by the
    // fromIndex traj point, where x = along the trajectory, y = transverse
    float rotAngle = work.Pts[fromIndex].Ang;
    float cs = cos(-rotAngle);
    float sn = sin(-rotAngle);
    if(prt) mf::LogVerbatim("CC")<<"FTM: Ang "<<rotAngle<<" cs "<<cs<<" sn "<<sn<<" work size "<<work.Pts.size()<<" CrawlDir "<<work.CrawlDir;
    
    // fit origin is the position of the fromIndex trajectory point
    std::array<float, 2> dir, origin = work.Pts[fromIndex].Pos;
    float xx, yy, terr;
    
    // protect against cases where we have a number of points with the same XX value
    unsigned short indx, nsame = 0;
    // Find the average charge
    tp.Chg = 0;
    float sum2 = 0;
    terr = std::abs(fScaleF * fHitErrFac * fAveHitRMS);
    if(terr < 0.01) terr = 0.01;
    for(unsigned short ipt = lo; ipt < hi; ++ipt) {
      indx = ipt - lo;
      if(ipt > work.Pts.size()-1 || indx > nTPF-1) {
        std::cout<<"bad ipt "<<ipt<<" work size "<<work.Pts.size()<<" or bad indx "<<indx<<" nTPF "<<nTPF<<" fAveHitRMS "<<fAveHitRMS<<"\n";
        return;
      }
      xx = work.Pts[ipt].HitPos[0] - origin[0];
      yy = work.Pts[ipt].HitPos[1] - origin[1];
      x[indx] = cs * xx - sn * yy;
      y[indx] = sn * xx + cs * yy;
      // Time error
//      terr = fScaleF * fHitErrFac * fAveHitRMS * work.Pts[ipt].ChgDiff;
      // Wire error^2 is a constant [1/sqrt(12)]^2 in WSE units = 0.0833
      yerr2[indx] = 0.0833 + terr * terr;
      if(prt) mf::LogVerbatim("CC")<<"FTM: ipt "<<ipt<<" indx "<<indx<<" xx "<<xx<<" yy "<<yy<<" x "<<x[indx]<<" y "<<y[indx]<<" err "<<yerr2[indx]<<" Chg "<<work.Pts[ipt].Chg;
      if(indx > 0 && std::abs(xx - x[indx-1]) < 1.E-3) ++nsame;
      tp.Chg += work.Pts[ipt].Chg;
      sum2 += work.Pts[ipt].Chg * work.Pts[ipt].Chg;
    } // ii
    float fcnt = nTPF;
    tp.Chg /= fcnt;
    float arg = sum2 - fcnt * tp.Chg * tp.Chg;
    if(arg > 0 && nTPF > 1) {
      tp.ChgDiff = sqrt(arg / (fcnt - 1));
    } else {
      tp.ChgDiff = 0.3;
    }
    
    if(nTPF < 4 && nsame > 0) return;
    
    float intcpt, slope, intcpterr, slopeerr, chidof;
    fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    if(chidof < 0.01) chidof = 0.01;
    tp.FitChi = chidof;
    if(chidof > 900) return;
    
    // calculate the new direction vector in the (xx, yy) coordinate system
    float newang = atan(slope);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    if(prt) mf::LogVerbatim("CC")<<"newdir XX,YY system "<<dir[0]<<" "<<dir[1]<<" slope "<<slope;
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAngle);
    sn = sin(rotAngle);
    if(prt) mf::LogVerbatim("CC")<<" lastPT ang "<<rotAngle<<" dir "<<cs<<" "<<sn;
    tp.Dir[0] = cs * dir[0] - sn * dir[1];
    tp.Dir[1] = sn * dir[0] + cs * dir[1];
    // decide whether to reverse the direction
    if(work.CrawlDir < 0) {
      tp.Dir[0] = -tp.Dir[0];
      tp.Dir[1] = -tp.Dir[1];
    }
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    // a good enough approximation since the trajectory can't change
    // too much from point to point
    tp.AngErr = std::abs(atan(slopeerr));
    if(prt) mf::LogVerbatim("CC")<<"  W,T dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" ang "<<tp.Ang;
    // rotate (0, intcpt) into W,T
    tp.Pos[0] = -sn * intcpt + origin[0];
    tp.Pos[1] =  cs * intcpt + origin[1];
    
    if(prt) {
      mf::LogVerbatim("CC")<<"FTM: fromIndex "<<fromIndex<<" Pos1 "<<tp.Pos[1]<<" Ang "<<tp.Ang<<" chi "<<chidof<<" old Ang "<<work.Pts[fromIndex].Ang;
//      PrintHeader();
//      PrintWork(fromIndex);
    }

  } // FitTrajMid
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::FitTraj()
  {
    // Fit the trajectory to a line. The hit positions are rotated into a coordinate
    // system with X along the direction of the last trajectory point and Y transverse.
    // The number of trajectory points fitted (nTPF) is determined by the calling routine.
    
    if(work.Pts.size() < 2) return;
    
    unsigned short lastPt = work.Pts.size() - 1;
    unsigned short ipt, ii;
    TrajPoint& lastTP = work.Pts[lastPt];
    
    if(lastTP.NumTPsFit == 2) {
      lastTP.Pos[0] = lastTP.HitPos[0];
      lastTP.Pos[1] = lastTP.HitPos[1];
      float dw = lastTP.HitPos[0] - work.Pts[0].Pos[0];
      float dt = lastTP.HitPos[1] - work.Pts[0].Pos[1];
      float nrm = sqrt(dw * dw + dt * dt);
      lastTP.Dir[0] = dw / nrm; lastTP.Dir[1] = dt / nrm;
      lastTP.Ang = atan2(lastTP.Dir[1], lastTP.Dir[0]);
      unsigned int iht = work.Pts[0].Hits[0];
      float minrms = fHits[iht].RMS();
      iht = work.Pts[1].Hits[0];
      if(fHits[iht].RMS() < minrms) minrms = fHits[iht].RMS();
      lastTP.AngErr = atan(2 * fHitErrFac * fScaleF * minrms);
      lastTP.FitChi = 0.01;
//      std::cout<<"chk2 from "<<work.Pts[lastPt-1].Pos[0]<<" "<<work.Pts[lastPt-1].Pos[1]<<" to "<<work.Pts[lastPt].Pos[0]<<" "<<work.Pts[lastPt].Pos[1];
//      std::cout<<" prev Ang "<<work.Pts[lastPt-1].Ang<<" new "<<work.Pts[lastPt].Ang<<" dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<"\n";
      return;
    }
    
    // try to fit with the number of points specified
    unsigned short nTPF = lastTP.NumTPsFit;
    // truncate if this is too many points
    if(nTPF > work.Pts.size()) {
      mf::LogWarning("CC")<<"FitTraj: NumTPsFit is too large "<<nTPF<<" work.Pts size "<<work.Pts.size()<<" Truncating...";
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
    float xx, yy, terr;
    
//    if(prt) mf::LogVerbatim("CC")<<"Enter with dir "<<work.Pts[lastPt].Dir[0]<<" "<<work.Pts[lastPt].Dir[1];

    // protect against cases where we have a number of points with the same X value
    unsigned short nsame = 0;
    terr = std::abs(fScaleF * fHitErrFac * fAveHitRMS);
    terr *= terr;
    if(terr < 0.01) terr = 0.01;
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
//      terr = fScaleF * fHitErrFac * fAveHitRMS * work.Pts[ipt].ChgDiff;
      // Wire error^2 is a constant [1/sqrt(12)]^2 in WSE units = 0.0833
      // TODO These should be added in quadrature...
      yerr2[ii] = std::abs(sn) * 0.0833 + std::abs(cs) * terr;
      wght = std::abs((work.Pts[ipt].Chg / qAve) - 1);
      if(wght < 0.01) wght = 0.01;
      yerr2[ii] *= wght;
//      if(yerr2[ii] < 0.001) yerr2[ii] = 0.001;
//      if(prt) mf::LogVerbatim("CC")<<"pt "<<ipt<<" xx "<<xx<<" yy "<<yy<<" x "<<x[ipt]<<" y "<<y[ipt]<<" err "<<yerr2[ipt];
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
//    if(prt) mf::LogVerbatim("CC")<<"  W,T dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" ang "<<lastTP.Ang;
    // rotate (0, intcpt) into W,T
    lastTP.Pos[0] = -sn * intcpt + origin[0];
    lastTP.Pos[1] =  cs * intcpt + origin[1];
//    if(prt) mf::LogVerbatim("CC")<<"  intcpt "<<intcpt<<" new pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1];
    
    if(prt) mf::LogVerbatim("CC")<<"Fit: "<<nTPF<<" pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<" dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" chi "<<chidof<<" AveHitRMS "<<fAveHitRMS;
    
  } // FitTraj
  
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::UpdateTrajChgDiff()
  {
    // calculate trajectory charge using nTrajPoints at the leading edge of the
    // trajectory
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];
    
    lastTP.ChgDiff = 99;
    unsigned short ii;
    
//    if(prt) std::cout<<"UpdateTrajChgDiff "<<work.Pts.size()<<"\n";
    
    if(work.Pts.size() < 2) return;
    
    float sum, sum2, arg, fcnt;
    if(work.Pts.size() < 4) {
      // just average for short trajectories
      sum = 0;
      sum2 = 0;
      fcnt = work.Pts.size();
      for(ii = 0; ii < work.Pts.size(); ++ii) {
        sum  += work.Pts[ii].Chg;
        sum2 += work.Pts[ii].Chg * work.Pts[ii].Chg;
      }
      lastTP.AveChg = sum / fcnt;
      arg = sum2 - fcnt * sum * sum;
      if(arg > 0 && work.Pts.size() > 2) {
        fChgRMS = sqrt(arg / (fcnt - 1));
      } else {
        fChgRMS = 0.3 * lastTP.AveChg;
      }
      if(fChgRMS < 10) fChgRMS = 10;
      lastTP.ChgDiff = (lastTP.Chg - lastTP.AveChg) / fChgRMS;
      if(std::abs(lastTP.ChgDiff) < 0.01) lastTP.ChgDiff = 0.01;
      return;
    }
    
    unsigned short ipt, cnt;
    
    if(lastTP.NumTPsFit > work.Pts.size()) {
      mf::LogError("CC")<<"UpdateTrajChgDiff: Invalid lastTP.NumTPsFit = "<<lastTP.NumTPsFit<<" work.Pts size "<<work.Pts.size();
      return;
    }
    
    // sort by charge and find the average using the lowest (<60%) charge hits at the
    // leading edge. Limit the number of trajectory points to allow for energy loss
    std::vector<float> chg;
    // Use the full length of the trajectory except for the first point
    unsigned short npts = lastPt;
    // no more than 20 points however
    if(npts > 20) npts = 20;
    for(ii = 0; ii < npts; ++ii) {
      ipt = lastPt - ii;
      chg.push_back(work.Pts[ipt].Chg);
    } // ii
    // default sort is lowest to highest
    std::sort(chg.begin(), chg.end());
    unsigned short maxcnt = (unsigned short)(0.6 * chg.size());
    if(maxcnt < 2) maxcnt = 2;
    sum = 0; sum2 = 0; cnt = 0;
    for(ipt = 0; ipt < maxcnt; ++ipt) {
      sum += chg[ipt];
      sum2 += chg[ipt] * chg[ipt];
      ++cnt;
    }
    fcnt = cnt;
    lastTP.AveChg = sum / fcnt;
    fChgRMS = sqrt((sum2 - fcnt * lastTP.AveChg * lastTP.AveChg) / (fcnt - 1));
    if(fChgRMS < 10) fChgRMS = 10;
    lastTP.ChgDiff = (lastTP.Chg - lastTP.AveChg) / fChgRMS;
    if(std::abs(lastTP.ChgDiff) < 0.01) lastTP.ChgDiff = 0.01;
//    if(prt) std::cout<<"ChgDiff "<<work.Pts.size()<<" "<<npts<<" "<<cnt<<" ave "<<lastTP.AveChg<<" rms "<<fChgRMS<<" diff "<<lastTP.ChgDiff<<"\n";
    
  } // UpdateTrajChgDiff
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::StartTraj(unsigned int fromHit, float toWire, float toTick)
  {
    // Start a simple (seed) trajectory going from a hit to a position (toWire, toTick).
    // The traj vector is cleared if an error occurs
    
    work.Pts.clear();
    work.ProcCode = 901; work.StopCode = 0;
    work.Vtx[0] = -1;  work.Vtx[1] = -1;
    if(fromHit > fHits.size() - 1) return;
    
    std::array<float, 2> pos, dir;
    
    // create a trajectory point
    TrajPoint tp;

    // position
    pos[0] = (float)fHits[fromHit].WireID().Wire;
    pos[1] = fHits[fromHit].PeakTime() * fScaleF;
    tp.HitPos = pos; tp.Pos = pos;
    // direction
    dir[0] = toWire - (float)fHits[fromHit].WireID().Wire;
    dir[1] = (toTick - fHits[fromHit].PeakTime()) * fScaleF;
    float ur = sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
    dir[0] /= ur; dir[1] /= ur;
    tp.Dir = dir;
    tp.Ang = atan2(dir[1], dir[0]);
    tp.AngErr = atan(fHitErrFac * fScaleF * fHits[fromHit].RMS());
    tp.Chg = fHits[fromHit].Integral();
    if(prt) mf::LogVerbatim("CC")<<"StartTraj "<<fHits[fromHit].WireID().Wire<<":"<<(int)fHits[fromHit].PeakTime()<<" -> "<<(int)toWire<<":"<<(int)toTick<<" dir "<<dir[0]<<" "<<dir[1]<<" ang "<<tp.Ang<<" angErr "<<tp.AngErr;
    
    // initialize other variables used by StepCrawl
    fAveHitRMS = fHits[fromHit].RMS();
    tp.AveChg = fHits[fromHit].Integral();
    fChgRMS = 0.3 * tp.AveChg;
    if(fChgRMS < 10) fChgRMS = 10;
    tp.ChgDiff = (tp.Chg - tp.AveChg) / fChgRMS;
    
    // Add the first trajectory hit
    tp.Hits.push_back(fromHit);
    inClus[fromHit] = -3;
    work.Pts.push_back(tp);
    work.CrawlDir = fStepCrawlDir;
    
  } // StartTraj
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::ReverseTraj()
  {
    // reverse the trajectory
    if(work.Pts.size() == 0) return;
    // reverse the crawling direction flag
    work.CrawlDir = -work.CrawlDir;
    std::reverse(work.Pts.begin(), work.Pts.end());
    // reverse the direction vector
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      if(work.Pts[ipt].Dir[0] != 0) work.Pts[ipt].Dir[0] = -work.Pts[ipt].Dir[0];
      if(work.Pts[ipt].Dir[1] != 0) work.Pts[ipt].Dir[1] = -work.Pts[ipt].Dir[1];
      work.Pts[ipt].Ang = atan2(work.Pts[ipt].Dir[1], work.Pts[ipt].Dir[0]);
      if(work.Pts[ipt].Hits.size() > 1) std::reverse(work.Pts[ipt].Hits.begin(), work.Pts[ipt].Hits.end());
    } // ipt
  }
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::StoreTraj()
  {
    // Creates a cluster using trajectory point hits
//    std::cout<<"StoreTraj "<<work.Pts.size()<<" Pos "<<work.Pts[0].Pos[0]<<" "<<work.Pts[0].Pos[1]<<"\n";
    // TODO: Handle the 2 point trajectory
    if(work.Pts.size() < 2) return;
    if(work.ProcCode == USHRT_MAX) return;
    if(work.Pts.size() == 0) {
      mf::LogWarning("CC")<<"StoreTraj: No points on work traj ";
      return;
    }
    
    // reverse the trajectory it's not in the proper ClusterCrawler order
    if(work.CrawlDir > 0) ReverseTraj();
    
    fcl2hits.clear();
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      if(work.Pts[ipt].Hits.size() == 0) {
        mf::LogWarning("CC")<<"StoreTraj: No hits on work traj point "<<ipt<<" Pts size "<<work.Pts.size();
        return;
      }
      fcl2hits.insert(fcl2hits.end(), work.Pts[ipt].Hits.begin(), work.Pts[ipt].Hits.end());
    } // ipt
    
    for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
      if(inClus[fcl2hits[ii]] > 0) {
        unsigned int hit = fcl2hits[ii];
        mf::LogWarning("CC")<<"StoreTraj trying to use an invalid hit "<<plane<<":"<<fHits[hit].WireID().Wire<<":"<<(int)fHits[hit].PeakTime()<<" with inClus "<<inClus[fcl2hits[ii]];
        return;
      }
      inClus[fcl2hits[ii]] = 0;
    } // ii
    
    if(work.Pts.size() > 2) {
      // Fit both ends of the work trajectory to get the slope errors right
      // Try to fit a length of 10 WSEs = 3 cm for uB
      float wantLength = 10;
      float Length = 0, dw, dt;
      unsigned short ii, toIndex;
      unsigned short fromIndex = work.Pts.size() - 1;
      unsigned short nPtsFit = 0;
      for(ii = 1; ii < work.Pts.size(); ++ii) {
        toIndex = work.Pts.size() - 1 - ii;
        dw = work.Pts[toIndex].Pos[0] - work.Pts[toIndex + 1].Pos[0];
        dt = work.Pts[toIndex].Pos[1] - work.Pts[toIndex + 1].Pos[1];
        Length += sqrt(dw * dw + dt * dt);
        ++nPtsFit;
        if(Length > wantLength) break;
      }
      TrajPoint tp;
/*
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - before";
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
*/
      FitTrajMid(fromIndex, toIndex, tp);
      if(tp.FitChi > 900) return;
      // Only use selected quantities from the fit
      work.Pts[fromIndex].Pos = tp.Pos;
      work.Pts[fromIndex].Chg = tp.Chg;
      work.Pts[fromIndex].ChgDiff = tp.ChgDiff;
      work.Pts[fromIndex].FitChi = tp.FitChi;
      work.Pts[fromIndex].NumTPsFit = tp.NumTPsFit;
      work.Pts[fromIndex].AngErr = tp.AngErr;
/*
      // Handle a sign change in the direction
      if(std::signbit(work.Pts[fromIndex].Ang) != std::signbit(tp.Ang)) {
        work.Pts[fromIndex].Dir[0] = -tp.Dir[0];
        work.Pts[fromIndex].Dir[1] = -tp.Dir[1];
        work.Pts[fromIndex].Ang = -tp.Ang;
      } else {
        work.Pts[fromIndex].Dir = tp.Dir;
        work.Pts[fromIndex].Ang = tp.Ang;
      }
*/
      
/*
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - after. fromIndex "<<fromIndex<<" size "<<work.Pts.size();
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
*/
      // fit the other end
      fromIndex = 0;
      toIndex = nPtsFit - 1;
/*
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - before";
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
*/
      FitTrajMid(fromIndex, toIndex, tp);
      if(tp.FitChi > 900) return;
      // Only use selected quantities from the fit
      work.Pts[fromIndex].Pos = tp.Pos;
      work.Pts[fromIndex].Chg = tp.Chg;
      work.Pts[fromIndex].ChgDiff = tp.ChgDiff;
      work.Pts[fromIndex].FitChi = tp.FitChi;
      work.Pts[fromIndex].NumTPsFit = tp.NumTPsFit;
      work.Pts[fromIndex].AngErr = tp.AngErr;
/*
      // Handle a sign change in the direction
      if(std::signbit(work.Pts[fromIndex].Dir[0]) != std::signbit(tp.Dir[0])) {
        work.Pts[fromIndex].Dir[0] = -tp.Dir[0];
        work.Pts[fromIndex].Dir[1] = -tp.Dir[1];
        work.Pts[fromIndex].Ang = -tp.Ang;
      } else {
        work.Pts[fromIndex].Dir = tp.Dir;
        work.Pts[fromIndex].Ang = tp.Ang;
      }
*/
 
/*
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - after. fromIndex "<<fromIndex<<" size "<<work.Pts.size();
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
*/
    }

    // Cluster Begin
    if(work.Pts[0].Dir[0] != 0) {
      clBeginSlp = work.Pts[0].Dir[1] / work.Pts[0].Dir[0];
    } else {
      if(work.Pts[0].Dir[1] > 0) {
        clBeginSlp = 999;
      } else {
        clBeginSlp = -999;
      }
    }
    clBeginSlpErr = std::abs(tan(work.Pts[0].AngErr));
    clBeginSlp /= fScaleF; clBeginSlpErr /= fScaleF;
    clBeginChg = work.Pts[0].AveChg;
    clBeginChgNear = 0;
    
    // Cluster End
    unsigned short lastPt = work.Pts.size() - 1;
    if(work.Pts[lastPt].Dir[0] != 0) {
      clEndSlp = work.Pts[lastPt].Dir[1] / work.Pts[lastPt].Dir[0];
    } else {
      if(work.Pts[lastPt].Dir[1] > 0) {
        clEndSlp = 999;
      } else {
        clEndSlp = -999;
      }
    }
    clEndSlpErr = std::abs(tan(work.Pts[lastPt].AngErr));
    clEndSlp /= fScaleF; clEndSlpErr /= fScaleF;
    clEndChg = work.Pts[lastPt].AveChg;
    clEndChgNear = 0;
    
    if(prt) mf::LogVerbatim("CC")<<"Store cluster of size "<<fcl2hits.size();
    clProcCode = work.ProcCode;
    // special stop code to indicate special crawling
    clStopCode = work.StopCode;
    // Store the cluster
    if(!TmpStore()) return;
    
//    std::cout<<"StoreTraj success \n";

    // store the trajectory and associated cluster index
    work.ClusterIndex = NClusters - 1;
    allTraj.push_back(work);
    
    // double check hit assignment
    unsigned short itj = allTraj.size() - 1;
    for(unsigned short ipt = 0; ipt < allTraj[itj].Pts.size(); ++ipt) {
      unsigned int hit = allTraj[itj].Pts[ipt].Hits[0];
      if(inClus[hit] != allTraj[itj].ClusterIndex+1) {
        mf::LogWarning("CC")<<"StoreTraj: Problem with inClus "<<inClus[hit]<<" != "<<allTraj[itj].ClusterIndex<<" itj "<<itj;
      }
    } // check

  } // StoreTraj
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::PrintAllTraj(unsigned short itj, unsigned short ipt)
  {
    
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      mf::LogVerbatim myprt("CC");
      myprt<<"TRJ Ind nPts     W:Tick    Ang    AveQ  Qual     W:T(WSE)  Ang    AveQ   Qual PDG     KE  Prim? 2ndry? \n";
      for(unsigned short ii = 0; ii < allTraj.size(); ++ii) {
        myprt<<"TRJ"<<std::fixed;
        myprt<<std::setw(4)<<ii<<std::setw(5)<<allTraj[ii].Pts.size();
        TrajPoint tp = allTraj[ii].Pts[0];
        unsigned short itick = tp.Pos[1]/fScaleF;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 1000) myprt<<" ";
        myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<tp.AveHitRes;
        unsigned short lastPt = allTraj[ii].Pts.size() - 1;
        tp = allTraj[ii].Pts[lastPt];
        itick = tp.Pos[1]/fScaleF;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 1000) myprt<<" ";
        myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<tp.AveHitRes;
        myprt<<std::setw(6)<<allTraj[ii].TruPDG;
        myprt<<std::setw(7)<<(int)allTraj[ii].TruKE;
        myprt<<std::setw(7)<<allTraj[ii].IsPrimary;
        myprt<<std::setw(7)<<allTraj[ii].IsSecondary;
        myprt<<"\n";
      } // ii
      return;
    } // itj > allTraj.size()-1
    
    if(itj > allTraj.size()-1) return;
    
    Trajectory tj = allTraj[itj];
    
    mf::LogVerbatim("CC")<<"Print allTraj["<<itj<<"]: ClusterIndex "<<tj.ClusterIndex<<" ProcCode "<<tj.ProcCode<<" StopCode "<<tj.StopCode<<" Vtx[0] "<<tj.Vtx[0]<<" Vtx[1] "<<tj.Vtx[1];
    
    PrintHeader();
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) PrintTrajPoint(ii, tj.CrawlDir, tj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(ipt, tj.CrawlDir, tj.Pts[ipt]);
    }
  } // PrintAllTraj

  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::PrintWork(unsigned short tPoint)
  {
    // prints one or all
    
    unsigned short first = 0;
    unsigned short last = work.Pts.size();
    if(tPoint == USHRT_MAX) {
      PrintHeader();
      for(unsigned short ipt = first; ipt < last; ++ipt) PrintTrajPoint(ipt, work.CrawlDir, work.Pts[ipt]);
    } else {
      // just print one traj point
      if(tPoint > work.Pts.size() -1) {
        mf::LogVerbatim("CC")<<"Cant print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(tPoint, work.CrawlDir, work.Pts[tPoint]);
    }
  } // PrintWork
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::PrintHeader()
  {
    mf::LogVerbatim("CC")<<"TRP  Ind  Stp     W:T(WSE) Tick Delta   RMS dTick    Ang   Err  Dir0       Q  QDiff    AveQ  FitChi NTPF NotNr UseNr?   Qual   Vtx Hits ";
  } // PrintHeader

  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::PrintTrajPoint(unsigned short ipt, short dir, TrajPoint tp)
  {
    mf::LogVerbatim myprt("CC");
    myprt<<"TRP"<<std::fixed;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(4)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]; // W:T
    if(tp.Pos[1] < 10) myprt<<"  "; if(tp.Pos[1] < 100) myprt<<" ";
    myprt<<std::setw(6)<<(int)(tp.Pos[1]/fScaleF); // Tick
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.DeltaRMS;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.Delta/fScaleF; // dTick
    int itick = (int)(tp.Delta/fScaleF);
    if(itick < 100) myprt<<" ";
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[0];
    myprt<<std::setw(8)<<(int)tp.Chg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgDiff;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(8)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NumTPsFit;
    myprt<<std::setw(6)<<tp.NumNotNear;
    myprt<<std::setw(6)<<tp.UsedNotCloseHit;
    myprt<<std::setw(8)<<std::setprecision(2)<<tp.AveHitRes;
    myprt<<std::setw(4)<<work.Vtx[0];
    myprt<<std::setw(4)<<work.Vtx[1];
    // print the hits associated with this traj point
    unsigned int hit;
    for(unsigned short iht = 0; iht < tp.Hits.size(); ++iht) {
      hit = tp.Hits[iht];
      if(hit > fHits.size()-1) {
        myprt<<" crazy hit ";
        continue;
      } // bad hit
      myprt<<" "<<fHits[hit].WireID().Wire<<":"<<(int)fHits[hit].PeakTime();
    }
  } // PrintTrajPoint

  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::HitMultipletPosition(unsigned int hit, float& hitTick)
  {
    // returns the charge weighted wire, time position of all hits in the multiplet
    // of which hit is a member
    
    std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(hit);
    float qtot = 0;
    hitTick = 0;
    for(size_t jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
      qtot += fHits[jht].Integral();
      hitTick += fHits[jht].Integral() * fHits[jht].PeakTime();
    }
    hitTick /= qtot;
    
  } // HitMultipletPosition
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::ReleaseTrajHits(unsigned short ipt)
  {
    // release the hits attached to trajectory point itj in the work vector
    if(ipt > work.Pts.size()-1) return;
    for(unsigned short iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = 0;
  } // ReleaseTrajHits
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::ReleaseAllWorkHits()
  {
    unsigned short ipt, iht;
    for(ipt = 0; ipt < work.Pts.size(); ++ipt) {
      for(iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = 0;
    } // ipt
  } // ReleaseAllWorkHits
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::FlagAllWorkHits()
  {
    unsigned short ipt, iht;
    for(ipt = 0; ipt < work.Pts.size(); ++ipt) {
      for(iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = -3;
    } // ipt
  } // ReleaseAllWorkHits

  
} // namespace cluster