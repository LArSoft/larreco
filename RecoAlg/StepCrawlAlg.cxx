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
    bool success;
    std::vector<std::pair<unsigned short, unsigned short>> kinkIndex;
    
    // Set this negative so it doesn't get used in ClusterHitsOK
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
      if(iwire < fFirstWire || iwire > fLastWire - 1) {
        std::cout<<"Bad iwire indexing "<<fFirstWire<<" "<<iwire<<" "<<fLastWire<<"\n";
        exit(1);
      }
      if(jwire < fFirstWire || jwire > fLastWire - 1) {
        std::cout<<"Bad jwire indexing "<<fFirstWire<<" "<<jwire<<" "<<fLastWire<<"\n";
        exit(1);
      }
      // skip bad wires or no hits on the wire
      if(WireHitRange[iwire].first < 0) continue;
      if(WireHitRange[jwire].first < 0) continue;
      unsigned int ifirsthit = (unsigned int)WireHitRange[iwire].first;
      unsigned int ilasthit = (unsigned int)WireHitRange[iwire].second;
      unsigned int jfirsthit = (unsigned int)WireHitRange[jwire].first;
      unsigned int jlasthit = (unsigned int)WireHitRange[jwire].second;
      for(iht = ifirsthit; iht < ilasthit; ++iht) {
        prt = (fDebugPlane == (int)plane && (int)iwire == fDebugWire && std::abs((int)fHits[iht].PeakTime() - fDebugHit) < 20);
        if(prt) mf::LogVerbatim("CC")<<"Found debug hit "<<plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime()<<" inClus "<<inClus[iht];
        if(inClus[iht] != 0) continue;
        // BIG TEMP
//        if(!prt) continue;
        for(jht = jfirsthit; jht < jlasthit; ++jht) {
          if(inClus[iht] != 0) continue;
          if(inClus[jht] != 0) continue;
          if(prt) mf::LogVerbatim("CC")<<" checking ClusterHitsOK with jht "<<plane<<":"<<fHits[jht].WireID().Wire<<":"<<(int)fHits[jht].PeakTime();
          // Ensure that the hits StartTick and EndTick have the proper overlap
          // stuff these into fcl2hits for ClusterHitsOK
          fcl2hits.resize(2);
          fcl2hits[0] = iht; fcl2hits[1] = jht;
          if(!ClusterHitsOK(-1)) continue;
          // start a trajectory in the direction from iht -> jht
          StartTraj(iht, jht);
/*
          if(prt) {
            mf::LogVerbatim("CC")<<"ClusterLoop2: First traj point";
            PrintWork(USHRT_MAX);
          }
*/
          // check for a failure
          if(work.Pts.size() == 0) {
            if(prt) mf::LogVerbatim("CC")<<"ClusterLoop2: StartTraj failed";
            continue;
          }
          StepCrawl(success);
          if(!success) {
            mf::LogVerbatim("CC")<<"ClusterLoop2: Bad return from StepCrawl";
//            mf::LogVerbatim("CC")<<"ClusterLoop2: Bad return from StepCrawl. Seed hit "<<plane<<":"<<fHits[work.Pts[0].Hits[0]].WireID().Wire<<":"<<(int)fHits[work.Pts[0].Hits[0]].PeakTime();
            ReleaseAllTrajHits();
            continue;
          }
          if(prt) mf::LogVerbatim("CC")<<"StepCrawl done: work.Pts size "<<work.Pts.size();
          if(work.Pts.size() < 2) {
            ReleaseAllTrajHits();
            continue;
          }
          StoreTraj();
          break;
        } // jht
      } // iht
    } // iwire
    
    prt = false;
    work.Pts.clear();
    
    // need to fake out FindVertices
    pass = 1;
    if(fFindTrajVertices) FindVertices();

    allTraj.clear();

    CheckHitClusterAssociations();
    
  } // ClusterLoop2
  
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
    
    if(prt) PrintWork(USHRT_MAX);
    
    unsigned int nit;
    // count number of steps taken with no trajectory point added
    unsigned short nMissedStep = 0, missedStepCut;
    float bigChgDiffCut = 3 * fStepCrawlChgDiffCut;
    float arg;
    bool stopCrawl, stopOnKink = false;
    unsigned short killPts, ipt, iht;
    for(nit = 0; nit < 10000; ++nit) {
      // move the position by one step in the right direction
      for(iwt = 0; iwt < 2; ++iwt) tp.Pos[iwt] += tp.Dir[iwt];
      // remove the old hits
      tp.Hits.clear();
      // look for new hits
      AddTrajHits(tp, stopCrawl);
      if(stopCrawl) break;
      if(tp.Hits.size() == 0) {
        ++nMissedStep;
        // Break if we took too many steps without adding a traj point
        arg = std::abs(tp.Dir[0]);
        if(arg < 0.1) arg = 0.1;
        missedStepCut = (unsigned short)(0.3 * (float)work.Pts.size());
        if(missedStepCut < 2) missedStepCut = 2;
        missedStepCut = (unsigned short)((float)missedStepCut / arg);
        if(prt) mf::LogVerbatim("CC")<<" nMissedStep "<<nMissedStep<<" cut "<<missedStepCut;
        if(nMissedStep > missedStepCut) {
          if(prt) mf::LogVerbatim("CC")<<" Too many steps w/o traj point ";
          break;
        }
        // See if we have missed more wires than the user has specified
        // get the last hit on the last traj point
        lastPt = work.Pts.size() - 1;
        if(std::abs(tp.Pos[0] - work.Pts[lastPt].Pos[0]) > fStepCrawlMaxWireSkip) {
          if(prt) mf::LogVerbatim("CC")<<" Too many wires missed. Last traj point wire "<<work.Pts[lastPt].Pos[0]<<" Pos[0] "<<tp.Pos[0]<<" fStepCrawlMaxWireSkip "<<fStepCrawlMaxWireSkip;
          break;
        }
        // Otherwise keep stepping
        continue;
      } // tp.Hits.size() == 0
      // Have new traj hits. Add the trajectory point and update
      tp.Step = nit;
      work.Pts.push_back(tp);
      UpdateTraj(success);
      if(!success) break;
      lastPt = work.Pts.size()-1;
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
/*
      // check for a kink. Stop crawling if one is found
      if(GottaKink()) {
        stopOnKink = true;
        break;
      }
*/
      // update the local tp unless we have killing to do
      if(killPts == 0) {
        tp = work.Pts[lastPt];
        if(prt) PrintWork(lastPt);
      } else {
        // shorten the trajectory by the desired number of points
        // and correct the projected position
        // release the hits
        ReleaseAllTrajHits();
        float nSteps = (float)(work.Pts[lastPt].Step - work.Pts[lastPt-killPts].Step);
        if(prt) {
          mf::LogVerbatim("CC")<<"   killing "<<killPts<<" points with "<<nSteps<<" steps";
          if(prt) mf::LogVerbatim("CC")<<"  Current tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1];
        }
        // truncate the work trajectory
        work.Pts.resize(lastPt+1-killPts);
        // re-assign the hits
        for(ipt = 0; ipt < work.Pts.size(); ++ipt) {
          for(iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = -3;
        } // ipt
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

    ReleaseAllTrajHits();
    
    // temp checking
    if(prt) {
      unsigned int firsthit = (unsigned int)WireHitRange[fFirstWire].first;
      unsigned int lasthit = (unsigned int)WireHitRange[fLastWire-1].second;
      mf::LogVerbatim myprt("CC");
      myprt<<"Not released hits ";
      for(unsigned int iht = firsthit; iht < lasthit; ++iht) {
        if(inClus[iht] == -3) {
          myprt<<" "<<plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime();
          inClus[iht] = 0;
        }
      }
    }
    
    success = true;
    return;
    
  } // StepCrawl
  
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
    std::vector<unsigned int> newHits;
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
    tp.UsedNotNearHit = false;
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
          if(hitDir > 0 && work.Pts[work.Pts.size()-1].Pos[1] < 0) continue;
          if(hitDir < 0 && work.Pts[work.Pts.size()-1].Pos[1] > 0) continue;
        }
        if(prt) mf::LogVerbatim("CC")<<"   ADD "<<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire<<":"<<(int)fHits[iht].PeakTime();
        newHits.push_back(iht);
        inClus[iht] = -3;
        // prepare to find the new position in WSE units
        sum += fHits[iht].Integral();
        hitsPos[0] += fHits[iht].Integral() * fHits[iht].WireID().Wire;
        hitsPos[1] += fHits[iht].Integral() * fHits[iht].PeakTime() * fScaleF;
      } // iht
    } // wire
    if(prt) mf::LogVerbatim("CC")<<" Done looking for hits in small window. newHits size "<<newHits.size();
    // didn't find a hit within the small window. Look for a hit in the larger window
    if(newHits.size() == 0) {
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
        newHits.push_back(nearbyHitIndex);
        inClus[nearbyHitIndex] = -3;
        // set the flag indicating this is a hit outside the smaller window
        tp.UsedNotNearHit = true;
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
          for(iht = 0; iht < newHits.size(); ++iht) inClus[newHits[iht]] = 0;
          stopCrawl = true;
          return;
        }
      } // nNearby == 0
    } // newHits.size() == 0
    // No hit found in the larger window
    if(newHits.size() == 0) return;
    // stuff the charge weighted position into tp and try to append it to the trajectory
    for(unsigned short iwt = 0; iwt < 2; ++iwt) tp.HitPos[iwt] = hitsPos[iwt] / sum;
    tp.Chg = sum;
    if(newHits.size() > 1) {
      // sort the new hits
      std::vector<SortEntry> sortVec;
      SortEntry sortEntry;
      unsigned short ii;
      bool sortByTime = (std::abs(tp.Dir[1]) > std::abs(tp.Dir[0]));
      for(ii = 0; ii < newHits.size(); ++ii) {
        sortEntry.index = ii;
        if(sortByTime) {
          // sort by time
          sortEntry.length = fHits[newHits[ii]].PeakTime();
        } else {
          // sort by wire
          sortEntry.length = (float)fHits[newHits[ii]].WireID().Wire;
        } // !SortByTime
        sortVec.push_back(sortEntry);
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), lessThan);
      // make a temp vector
      std::vector<unsigned int> tmp = newHits;
      // overwrite with the sorted values
      for(ii = 0; ii < sortVec.size(); ++ii) tmp[ii] = newHits[sortVec[ii].index];
      // replace
      newHits = tmp;
      // swap the order? TODO: check
      if(sortByTime && tp.Dir[1] < 0) std::reverse(newHits.begin(), newHits.end());
      if(!sortByTime && tp.Dir[0] < 0) std::reverse(newHits.begin(), newHits.end());
    }
    // Now have a good trajectory point. Add the hits to the traj
    tp.Hits = newHits;
/*
    work.Pts.push_back(tp);
    UpdateTraj(success);
    if(!success) {
      mf::LogVerbatim("CC")<<" UpdateTraj failed. Stop crawling ";
      ReleaseAllTrajHits();
      return;
    }
*/
    stopCrawl = false;
  } // AddTrajHits
  
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
    if(tp1.FitChi > 900) {
      unsigned int hit = tp1.Hits[0];
      mf::LogWarning("CC")<<"GottaKink: Bad traj fit1. Seed hit"<<plane<<":"<<fHits[hit].WireID().Wire<<":"<<(int)fHits[hit].PeakTime();
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
    if(prt) mf::LogVerbatim("CC")<<"TRJ kink angle "<<dang<<" cut "<<cut<<" work.Pts old size "<<work.Pts.size()<<" new size "<<kinkIndex;
    ReleaseAllTrajHits();
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
    
    if(prt) mf::LogVerbatim("CC")<<"StepCrawlClusterCheck: AveHitRes "<<work.Pts[work.Pts.size()-1].AveHitRes;
    // ignore lower quality clusters
    if(work.Pts[work.Pts.size()-1].AveHitRes > 1.2) return;
    
    float delrat;
    unsigned short ipt;
    bool lopped = false;
    for(ipt = work.Pts.size() - 3; ipt < work.Pts.size(); ++ipt) {
      // check for a largish change in Delta
      delrat = work.Pts[ipt].Delta / work.Pts[ipt-1].Delta;
//      if(prt) mf::LogVerbatim("CC")<<" ipt "<<ipt<<" delrat "<<delrat;
      if(delrat > 2) {
        ReleaseAllTrajHits();
        work.Pts.resize(ipt);
        if(prt) mf::LogVerbatim("CC")<<" new work.Pts size "<<work.Pts.size();
        lopped = true;
        break;
      }
    } // ipt
    
    if(!lopped) return;
    
    // Check for a missing wire at the end - small angle cluster
    ipt = work.Pts.size() - 1;
    if(prt) mf::LogVerbatim("CC")<<" check last point. "<<ipt<<" Dir "<<work.Pts[ipt].Dir[1]<<" dPos[0] "<<std::abs(work.Pts[ipt].HitPos[0] - work.Pts[ipt-1].HitPos[0]);
    if(std::abs(work.Pts[ipt].Dir[1]) > 0.95) return;
    if(std::abs(work.Pts[ipt].HitPos[0] - work.Pts[ipt-1].HitPos[0]) > 1.5) work.Pts.pop_back();
    if(prt) {
      mf::LogVerbatim("CC")<<" lopped one more work.Pts point new work.Pts size "<<work.Pts.size();
//      PrintWork(USHRT_MAX);
    }
    
  } // StepCrawlClusterCheck
  
  
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::GetStepCrawlWindow(TrajPoint& tp, unsigned int& loWire, unsigned int& hiWire, float& loTime, float& hiTime,
                                             unsigned int& loloWire, unsigned int& hihiWire, float& loloTime, float& hihiTime)
  {
    // Define two windows for accepting hits within +/- 0.5 wire spacing units
    loWire = 0; hiWire = 0; loTime = 0; hiTime = 0;
    loloWire = 0; hihiWire = 0; loloTime = 0; hihiTime = 0;
    
    if(tp.Dir[0] == 0 && tp.Dir[1] == 0) return;
    
    // Find the positions of the error ellipse
    std::array<float, 2> pos;
    float loPos0 = 99999, hiPos0 = 0;
    loTime = 99999, hiTime = 0;
    // move 1/2 step in each direction and +/- 3 sigma angle error
    float dang = 3 * tp.AngErr;
    //    if(prt) std::cout<<"tp.pos "<<tp.Pos[0]<<" "<<tp.Pos[0]<<" tp.Ang "<<tp.Ang<<" dang "<<dang<<"\n";
    pos[0] = tp.Pos[0] - 0.5 * cos(tp.Ang - dang);
    if(pos[0] < loPos0) loPos0 = pos[0]; if(pos[0] > hiPos0) hiPos0 = pos[0];
    pos[1] = tp.Pos[1] - 0.5 * sin(tp.Ang - dang);
    if(pos[1] < loTime) loTime = pos[1]; if(pos[1] > hiTime) hiTime = pos[1];
    //    if(prt) std::cout<<"pos "<<pos[0]<<" "<<pos[1]<<"\n";
    pos[0] = tp.Pos[0] - 0.5 * cos(tp.Ang + dang);
    if(pos[0] < loPos0) loPos0 = pos[0]; if(pos[0] > hiPos0) hiPos0 = pos[0];
    pos[1] = tp.Pos[1] - 0.5 * sin(tp.Ang + dang);
    if(pos[1] < loTime) loTime = pos[1]; if(pos[1] > hiTime) hiTime = pos[1];
    //    if(prt) std::cout<<"pos "<<pos[0]<<" "<<pos[1]<<"\n";
    
    pos[0] = tp.Pos[0] + 0.5 * cos(tp.Ang - dang);
    if(pos[0] < loPos0) loPos0 = pos[0]; if(pos[0] > hiPos0) hiPos0 = pos[0];
    pos[1] = tp.Pos[1] + 0.5 * sin(tp.Ang - dang);
    if(pos[1] < loTime) loTime = pos[1]; if(pos[1] > hiTime) hiTime = pos[1];
    //    if(prt) std::cout<<"pos "<<pos[0]<<" "<<pos[1]<<"\n";
    pos[0] = tp.Pos[0] + 0.5 * cos(tp.Ang + dang);
    if(pos[0] < loPos0) loPos0 = pos[0]; if(pos[0] > hiPos0) hiPos0 = pos[0];
    pos[1] = tp.Pos[1] + 0.5 * sin(tp.Ang + dang);
    if(pos[1] < loTime) loTime = pos[1]; if(pos[1] > hiTime) hiTime = pos[1];
    //    if(prt) std::cout<<"pos "<<pos[0]<<" "<<pos[1]<<"\n";
    
    
    // convert to ticks
    loTime /= fScaleF;
    hiTime /= fScaleF;
    //    std::cout<<"ticks loTime "<<loTime<<" "<<hiTime;
    // Add a hit time uncertainty in ticks
    float hitTimeErr = 3 * fHitErrFac * fAveHitRMS;
    loTime -= hitTimeErr;
    hiTime += hitTimeErr;
    
    // convert pos[0] to wire number
    loWire = (unsigned int)(loPos0 + 0.5);
    hiWire = (unsigned int)(hiPos0 + 1.5);
//    if(prt) std::cout<<"lo,hi Wire "<<loWire<<" "<<hiWire<<"\n";
    
    // prevent looking on an already considered wire for not-too-large angle clusters
    if(std::abs(tp.Dir[1]) < 0.7 && work.Pts.size() > 0 && work.Pts[work.Pts.size()-1].Hits.size() > 0) {
      unsigned short lastPt = work.Pts.size() - 1;
      unsigned short lastPtLastHit = work.Pts[lastPt].Hits.size() - 1;
      unsigned int lastHit = work.Pts[lastPt].Hits[lastPtLastHit];
      unsigned int lastHitWire = fHits[lastHit].WireID().Wire;
      if(tp.Dir[0] > 0) {
        // moving downstream
        if(loWire < lastHitWire + 1) loWire = lastHitWire + 1;
        hiWire = loWire + 1;
      } else {
        // moving upstream
        if(hiWire > lastHitWire - 1) hiWire = lastHitWire;
        loWire = hiWire - 1;
      }
      // define larger window in Time only
      loloWire = loWire; loloTime = loTime - 2 * hitTimeErr;
      hihiWire = hiWire; hihiTime = hiTime + 2 * hitTimeErr;
    } // !large angle
    else {
      // add an extra wire on both sides
      loloWire = loWire - 1;
      hihiWire = hiWire + 1;
      // increase the size of the time window. Ignore the small hit time uncertainty
      float time = (hiTime + loTime) /2;
      float win = 3 * (hiTime - time);
      loloTime = time - win;
      hihiTime = time + win;
    } // large angle
    
    if(prt) {
      mf::LogVerbatim("CC")<<"loWire   "<<loWire<<" hiWire   "<<hiWire<<" loTime   "<<(int)loTime<<" hiTime   "<<(int)hiTime<<" hitTimeErr "<<(int)hitTimeErr<<" ticks";
      mf::LogVerbatim("CC")<<"loloWire "<<loloWire<<" hihiWire "<<hihiWire<<" loloTime "<<(int)loloTime<<" hihiTime "<<(int)hihiTime;
      
    }
  } // GetStepCrawlWindow
  
  
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
    
    unsigned short ii, cnt;
    std::array<float, 2> dpos;
    float path;
    
    // Handle the second trajectory point. No error calculation. Just copy
    // the previous
    if(work.Pts.size() == 2) {
      // second trajectory point is the hit position
      lastTP.Pos = lastTP.HitPos;
      for(ii = 0; ii < 2; ++ii) dpos[ii] = lastTP.Pos[ii] - work.Pts[0].Pos[ii];
      path = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1]);
      if(path == 0) return;
      // calculate the direction and angle
      for(ii = 0; ii < 2; ++ii) lastTP.Dir[ii] = dpos[ii] / path;
      lastTP.Ang = atan2(lastTP.Dir[1], lastTP.Dir[0]);
      lastTP.AngErr = work.Pts[0].AngErr;
      lastTP.NumTPsFit = 2;
      UpdateTrajChgDiff();
      lastTP.FitChi = 0;
      if(prt) mf::LogVerbatim("CC")<<"UpdateTraj: Second traj point pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1];
      success = true;
      return;
    }

    
    // Find Delta for the last trajectory point before we fit it.
    // Ensure that we aren't moving exactly in the +/- X direction
    if(lastTP.Dir[0] != 0) {
      // Find the Distance Of Closest Approach
      float tslp = lastTP.Dir[1] / lastTP.Dir[0];
      float docaW = (lastTP.HitPos[0] + tslp * (lastTP.HitPos[1]-lastTP.Pos[1]) + lastTP.Pos[0] * tslp * tslp) /
      (1 + tslp * tslp);
      float docaT = lastTP.Pos[1] + (lastTP.HitPos[0] -  lastTP.Pos[0]) * tslp;
      dpos[0] = docaW - lastTP.HitPos[0];
      dpos[1] = docaT - lastTP.HitPos[1];
      lastTP.Delta = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1]);
      // don't let Delta be smaller than the user specified fHitErrFace * RMS
      float minDelta = fHitErrFac * fAveHitRMS * fScaleF;
      if(lastTP.Delta < minDelta) lastTP.Delta = minDelta;
    } else {
      // trajectory going in the drift direction
      lastTP.Delta = std::abs(lastTP.HitPos[0] - lastTP.Pos[0]);
    }

    // Re-fit the trajectory with the newly added point. Determine how many
    // trajectory points should be fit to a line. Start with the number used in
    // the previous trajectory point
    unsigned short nTPF = work.Pts[lastPt].NumTPsFit;
    if(prt) mf::LogVerbatim("CC")<<"UpdateTraj: Starting nTPF "<<nTPF;
    // decide if this should be shorter or longer
    // Count the number of previous traj points that used "nearby" hits
    cnt = 0;
    for(ii = 0; ii < nTPF; ++ii) {
      if(!work.Pts[lastPt - ii].UsedNotNearHit) break;
      ++cnt;
    } // ii
    // count of all nearby hits used in the trajectory fit
    unsigned short totCnt = 0;
    for(ii = 0; ii < nTPF; ++ii) {
      if(work.Pts[lastPt - ii].UsedNotNearHit) ++totCnt;
    } // ii
    // Handle the case where reducing the number of points fit doesn't work
    bool inTrouble = false;
    if(cnt == 0) {
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
      inTrouble = true;
    }
    // check for craziness
    if(nTPF < 2) nTPF = 2;
    lastTP.NumTPsFit = nTPF;
    
    // fill the charge information so we can use it to weight the fit
    UpdateTrajChgDiff();
    
    FitTraj();
    // check for a failure
    if(lastTP.FitChi > 900) return;
    
    if(inTrouble) {
      // set the angle error large
      lastTP.AngErr = 0.5;
    } else {
      // calculate new angle error
      float fcnt;
      if(work.Pts.size() > 4) {
        // Update the angle error the hits that were fit after we have a few points.
        float sum = 0, delta;
        unsigned short ipt;
        cnt = 0;
        // use small angle approximation
        for(ii = 0; ii < nTPF - 3; ++ii) {
          ipt = lastPt - ii;
          // The first two points don't have a real Delta
          if(ipt < 3) break;
          delta = work.Pts[ipt].Delta;
          if(delta == 0) delta = 0.1;
          dpos[0] = work.Pts[ipt].Pos[0] - work.Pts[ipt-1].Pos[0];
          dpos[1] = work.Pts[ipt].Pos[1] - work.Pts[ipt-1].Pos[1];
          path = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1]);
          sum += delta / path;
          ++cnt;
        }
        if(cnt > 1) {
          fcnt = cnt;
          lastTP.AngErr = sum / fcnt;
        }
      } // work.Pts.size() > 3
    }

    std::vector<unsigned int> tHits;
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt)
      tHits.insert(tHits.end(), work.Pts[ipt].Hits.begin(), work.Pts[ipt].Hits.end());
    // Calculate the average RMS width of hits on the cluster
    fAveHitRMS = 0;
    for(unsigned short ii = 0; ii < tHits.size(); ++ii) fAveHitRMS += fHits[tHits[ii]].RMS();
    fAveHitRMS /= (float)tHits.size();
    
    // Calculate the rms scatter of hit times
    float aveRes = 0, arg;
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
      aveRes += arg * arg;
      ++cnt;
    }
    if(cnt > 0) {
      aveRes /= (float)cnt;
      aveRes = sqrt(aveRes);
      // convert to a quality factor 1 = perfect, > 1
      aveRes /= (fAveHitRMS * fHitErrFac);
      lastTP.AveHitRes = aveRes;
      //      if(prt) mf::LogVerbatim("CC")<<std::right<<std::setw(6)<<std::fixed<<std::setprecision(1)<<" aveRes "<<aveRes<<" fAveHitRMS "<<fAveHitRMS;
    } else {
      if(prt) mf::LogVerbatim("CC")<<"    NA tHits size "<<tHits.size()<<" cnt "<<cnt;
    }
    
    success = true;
    return;
    
  } // UpdateTraj

  //////////////////////////////////////////
  void ClusterCrawlerAlg::FitTrajMid(unsigned short fromIndex, unsigned short toIndex, TrajPoint& tp)
  {
    // Fit the trajectory points from fromIndex to (and including) toIndex. The origin
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
    float cs = cos(-work.Pts[fromIndex].Ang);
    float sn = sin(-work.Pts[fromIndex].Ang);
    
    // fit origin is the position of the current trajectory point
    std::array<float, 2> dir, origin = work.Pts[fromIndex].HitPos;
    float xx, yy, terr;
    
    // protect against cases where we have a number of points with the same XX value
    unsigned short indx, nsame = 0;
    // Find the average charge
    tp.Chg = 0;
    for(unsigned short ipt = lo; ipt < hi; ++ipt) {
      indx = ipt - lo;
      xx = work.Pts[ipt].HitPos[0] - origin[0];
      yy = work.Pts[ipt].HitPos[1] - origin[1];
      x[indx] = cs * xx - sn * yy;
      y[indx] = sn * xx + cs * yy;
      // Time error
      terr = fScaleF * fHitErrFac * fAveHitRMS * work.Pts[ipt].ChgDiff;
      // Wire error^2 is a constant [1/sqrt(12)]^2 in WSE units = 0.0833
      yerr2[indx] = 0.0833 + terr * terr;
//      if(prt) mf::LogVerbatim("CC")<<"pt "<<ipt<<" xx "<<xx<<" yy "<<yy<<" x "<<x[ipt]<<" y "<<y[ipt]<<" err "<<yerr2[ipt];
      if(indx > 0 && std::abs(xx - x[indx-1]) < 1.E-3) ++nsame;
      tp.Chg += work.Pts[ipt].Chg;
    } // ii
    tp.Chg /= (float)nTPF;
    
    if(nTPF < 4 && nsame > 0) return;
    
    float intcpt, slope, intcpterr, slopeerr, chidof;
    fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    
    tp.FitChi = chidof;
    if(chidof > 900) return;
    
    // calculate the new direction vector in the (xx, yy) coordinate system
    float newang = atan(slope);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
//    if(prt) mf::LogVerbatim("CC")<<"newdir XX,YY system "<<dir[0]<<" "<<dir[1]<<" slope "<<slope;
    // rotate back into the (w,t) coordinate system
    cs = cos(work.Pts[fromIndex].Ang);
    sn = sin(work.Pts[fromIndex].Ang);
//    if(prt) mf::LogVerbatim("CC")<<" lastPT ang "<<work.Pts[fromIndex].Ang<<" dir "<<cs<<" "<<sn;
    tp.Dir[0] = cs * dir[0] - sn * dir[1];
    tp.Dir[1] = sn * dir[0] + cs * dir[1];
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    // a good enough approximation since the trajectory can't change
    // too much from point to point
    tp.AngErr = std::abs(atan(slopeerr));
//    if(prt) mf::LogVerbatim("CC")<<"  W,T dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" ang "<<lastTP.Ang;
    // rotate (0, intcpt) into W,T
    tp.Pos[0] = -sn * intcpt + origin[0];
    tp.Pos[1] =  cs * intcpt + origin[1];
    
//    if(prt) mf::LogVerbatim("CC")<<"  intcpt "<<intcpt<<" new pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1];

  } // FitTrajMid
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::FitTraj()
  {
    // Fit the trajectory to a line. The hit positions are rotated into a coordinate
    // system with X along the direction of the last trajectory point and Y transverse.
    // The number of trajectory points fitted (nTPF) is determined by the calling routine.
    
    unsigned short lastPt = work.Pts.size() - 1;
    unsigned short ipt;
    TrajPoint& lastTP = work.Pts[lastPt];
    
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
    float cs = cos(-work.Pts[lastPt].Ang);
    float sn = sin(-work.Pts[lastPt].Ang);
    
    // fit origin is the position of the current trajectory point
    // TODO Check this. HitPos or Pos?
    std::array<float, 2> dir, origin = work.Pts[lastPt].HitPos;
    float xx, yy, terr;
    
//    if(prt) mf::LogVerbatim("CC")<<"Enter with dir "<<work.Pts[lastPt].Dir[0]<<" "<<work.Pts[lastPt].Dir[1];

    // protect against cases where we have a number of points with the same X value
    unsigned short nsame = 0;
    for(unsigned short ii = 0; ii < nTPF; ++ii) {
      ipt = lastPt - ii;
      xx = work.Pts[ipt].HitPos[0] - origin[0];
      yy = work.Pts[ipt].HitPos[1] - origin[1];
      x[ii] = cs * xx - sn * yy;
      y[ii] = sn * xx + cs * yy;
      // Time error
      terr = fScaleF * fHitErrFac * fAveHitRMS * work.Pts[ipt].ChgDiff;
      // Wire error^2 is a constant [1/sqrt(12)]^2 in WSE units = 0.0833
      // TODO: This is the **right way** but check it
//      yerr2[ii] = cs * 0.0833 + sn * (terr * terr);
      yerr2[ii] = 0.0833 + terr * terr;
//      if(prt) mf::LogVerbatim("CC")<<"pt "<<ipt<<" xx "<<xx<<" yy "<<yy<<" x "<<x[ipt]<<" y "<<y[ipt]<<" err "<<yerr2[ipt];
      if(ii > 0 && std::abs(xx - x[ii-1]) < 1.E-3) ++nsame;
    } // ii
    if(nTPF < 4 && nsame > 0) return;
    
    float intcpt, slope, intcpterr, slopeerr, chidof;
    fLinFitAlg.LinFit(x, y, yerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    
    lastTP.FitChi = chidof;
    if(chidof > 900) return;
    
    // calculate the new direction vector in the (xx, yy) coordinate system
    // TODO: Store this with the trajectory so it doesn't need to be calculated on the next step?
    float newang = atan(slope);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    // rotate back into the (w,t) coordinate system
    cs = cos(work.Pts[lastPt].Ang);
    sn = sin(work.Pts[lastPt].Ang);
    lastTP.Dir[0] = cs * dir[0] - sn * dir[1];
    lastTP.Dir[1] = sn * dir[0] + cs * dir[1];
    lastTP.Ang = atan2(lastTP.Dir[1], lastTP.Dir[0]);
//    if(prt) mf::LogVerbatim("CC")<<"  W,T dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" ang "<<lastTP.Ang;
    // rotate (0, intcpt) into W,T
    lastTP.Pos[0] = -sn * intcpt + origin[0];
    lastTP.Pos[1] =  cs * intcpt + origin[1];
//    if(prt) mf::LogVerbatim("CC")<<"  intcpt "<<intcpt<<" new pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1];
    
    if(prt) mf::LogVerbatim("CC")<<"Fit: "<<nTPF<<" pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<" dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" chi "<<chidof;
    
  } // FitTraj
  
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::UpdateTrajChgDiff()
  {
    // calculate trajectory charge using nTrajPoints at the leading edge of the
    // trajectory
    unsigned int lastPt = work.Pts.size()-1;
    TrajPoint& lastTP = work.Pts[lastPt];
    
    lastTP.ChgDiff = -99;
    unsigned short ii;
    
//    if(prt) std::cout<<"UpdateTrajChgDiff "<<work.Pts.size()<<"\n";
    
    if(work.Pts.size() < 4) {
      // just average for short trajectories
      lastTP.AveChg = 0;
      for(ii = 0; ii < work.Pts.size(); ++ii) lastTP.AveChg += work.Pts[ii].Chg;
      lastTP.AveChg /= (float)work.Pts.size();
      fChgRMS = 0.3 * lastTP.AveChg;
      lastTP.ChgDiff = (lastTP.Chg - lastTP.AveChg) / fChgRMS;
      return;
    }
    
    float sum, sum2, fcnt;
    unsigned short ipt, cnt;
    
    if(lastTP.NumTPsFit > work.Pts.size()) {
      mf::LogError("CC")<<"UpdateTrajChgDiff: Invalid lastTP.NumTPsFit = "<<lastTP.NumTPsFit<<" work.Pts size "<<work.Pts.size();
      return;
    }
    
    // sort by charge and find the average using the lowest (<60%) charge hits at the
    // leading edge. Limit the number of trajectory points to allow for energy loss
    std::vector<float> chg;
    unsigned short npts = lastTP.NumTPsFit;
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
    lastTP.ChgDiff = (lastTP.Chg - lastTP.AveChg) / fChgRMS;
//    if(prt) std::cout<<"ChgDiff "<<work.Pts.size()<<" "<<npts<<" "<<cnt<<" ave "<<lastTP.AveChg<<" rms "<<fChgRMS<<" diff "<<lastTP.ChgDiff<<"\n";
    
  } // UpdateTrajChgDiff
  
  //////////////////////////////////////////
/*
  void ClusterCrawlerAlg::FindVLAClusters()
  {
    // Find Very Large Angle clusters after clustering is done in a plane.
    // Signatures are long hit multiplets that are not included in a cluster.
    // and possibly nearby large angle clusters
    
    unsigned int wire, iht, jht;
    unsigned short itj, end;
    unsigned int iFirstHit, iLastHit;
    //    std::vector<std::array<float, 2>> traj;
    TrajPoint trajPoint;
    std::array<float, 2> pos, dir;
    float slp;
    unsigned short nAvail;
    bool success;
    for(wire = fFirstWire; wire < fLastWire; ++wire) {
      if(WireHitRange[wire].first < 0) continue;
      iFirstHit = (unsigned int)WireHitRange[wire].first;
      iLastHit = (unsigned int)WireHitRange[wire].second;
      fcl2hits.clear();
      prt = (fDebugPlane == (int)plane && (unsigned int)fDebugWire == wire && fDebugHit == 7777);
      if(prt) mf::LogVerbatim("CC")<<"FindVLACluster: looking for multiplet";
      nAvail = 0;
      for(iht = iFirstHit; iht < iLastHit; ++iht) {
        // ignore obsolete hits
        if(inClus[iht] < 0) continue;
        if(fHits[iht].Multiplicity() == 1) continue;
        // Only look at the first hit in a multiplet
        if(fHits[iht].Multiplicity() > 1 && fHits[iht].LocalIndex() > 0) continue;
        std::pair<size_t, size_t> MultipletRange = FindHitMultiplet(iht);
        for(jht = MultipletRange.first; jht < MultipletRange.second; ++jht) {
          // look for 2 or more adjacent un-used hits
          if(prt) mf::LogVerbatim("CC")<<"Hit"<<fHits[jht].WireID().Plane<<":"<<fHits[jht].WireID().Wire
            <<":"<<(int)fHits[jht].PeakTime()<<" inClus "<<inClus[jht]<<" nAvail "<<nAvail;
          if(nAvail == 1 && inClus[jht] > 0) nAvail = 0;
          if(inClus[jht] != 0) continue;
          if(mergeAvailable[jht]) ++nAvail;
          // only grab the first 2 hits
          if(fcl2hits.size() < 2) fcl2hits.push_back(jht);
        }
      } // iht
      if(prt) mf::LogVerbatim("CC")<<"Found multiplet "<<fHits[iht].WireID().Plane<<":"<<fHits[iht].WireID().Wire
        <<":"<<(int)fHits[iht].PeakTime()<<" nAvail "<<nAvail;
      if(nAvail < 2) continue;
      iht = fcl2hits[0];
      // step crawl in the +time direction
      work.Pts.clear();
      // first traj point at the first hit
      iht = fcl2hits[0];
      pos[0] = (float)fHits[iht].WireID().Wire;
      pos[1] = fHits[iht].PeakTime() * fScaleF;
      trajPoint.HitPos = pos;
      trajPoint.Pos = pos;
      dir[0] = 0; dir[1] = 1;
      trajPoint.Dir = dir;
      trajPoint.Ang = atan2(dir[1], dir[0]);
      // Calculate the angular error using dir for dirErr
      // wire error
      dir[0] = 2;
      // time error in wire units
      dir[1] = (fHits[fcl2hits[1]].PeakTime() - fHits[fcl2hits[0]].PeakTime()) * fScaleF;
      trajPoint.AngErr = std::abs(atan(dir[0] / dir[1]));
      trajPoint.Delta = 0;
      trajPoint.Chg = fHits[iht].Integral();
      work.Pts.push_back(trajPoint);
      // clobber all the hits except for the first. Let StepCrawl pick up
      // hits on adjacent wires in the correct order
      fcl2hits.resize(1);
      StepCrawl(success);
      if(!success) continue;
      // reverse the direction and step in the -time direction
      //      ReverseTraj();
      //      StepCrawl(step);
      if(prt) mf::LogVerbatim("CC")<<"StepCrawl done: fcl2hits "<<fcl2hits.size()<<" traj size "<<traj.size();
      if(fcl2hits.size() < 3) continue;
      if(traj.size() < 2) continue;
      // Define the needed variables for storing
      clBeginChg = -1; clEndChg = -1; // let TmpStore do this
      // See if the trajectory is in proper order such that end 0 = Begin = larger wire number
      itj = work.Pts.size() - 1;
      bool properOrder = (work.Pts[itj].Pos[0] < work.Pts[0].Pos[0]);
      if(!properOrder) ReverseTraj();
      // Define the clBegin (end = 0) and clEnd (end = 1) cluster parameters
      for(end = 0; end < 2; ++end) {
        if(end == 0) {
          dir = work.Pts[0].Dir;
        } else {
          dir = work.Pts[work.Pts.size()-1].Dir;
        }
        if(std::abs(dir[0]) < 1e-3) {
          // Slope way too large - fake it
          if(dir[0] > 0) { slp = 999; } else { slp = 999; }
        } else {
          slp = dir[1] / dir[0];
        } // std::abs(dir[0]) < 1e-3
        // Define clBegin and clEnd parameters
        // TODO: More careful treatment of slope errors
        if(end == 0) {
          clBeginSlp = slp; clBeginSlpErr = 0.1 * std::abs(slp); clBeginChgNear = 0;
        } else {
          clEndSlp = slp; clEndSlpErr = 0.1 * std::abs(slp); clEndChgNear = 0;
        }
      } // end
      if(prt) mf::LogVerbatim("CC")<<"Store cluster with size "<<fcl2hits.size();
      clProcCode = 999;
      // special stop code to indicate special crawling
      clStopCode = 8;
      TmpStore();
    } // wire
    
    
    for(iht = 0; iht < fHits.size(); ++iht) if(inClus[iht] == -3) inClus[iht] = 0;
    fcl2hits.clear();
    
  } // FindVLAClusters
*/
  
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::StartTraj(unsigned int fromHit, unsigned int toHit)
  {
    // Start a simple (seed) trajectory using two hits.
    // The direction of the trajectory is from fromHit to toHit and consists of
    // one trajectory points composed of hits in the vicinity of fromHit. The traj
    // vector is cleared if an error occurs
    
    work.Pts.clear();
    work.ProcCode = 901; work.StopCode = 0;
    work.Vtx[0] = -1;  work.Vtx[1] = -1;
    if(fromHit > fHits.size() || toHit > fHits.size()) return;
    
    std::array<float, 2> pos, dir;
    
    // create a trajectory point
    TrajPoint tp;

    // position
    pos[0] = (float)fHits[fromHit].WireID().Wire;
    pos[1] = fHits[fromHit].PeakTime() * fScaleF;
    tp.HitPos = pos; tp.Pos = pos;
    // direction
    dir[0] = (float)fHits[toHit].WireID().Wire - (float)fHits[fromHit].WireID().Wire;
    dir[1] = (fHits[toHit].PeakTime() - fHits[fromHit].PeakTime()) * fScaleF;
    float ur = sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
    dir[0] /= ur; dir[1] /= ur;
    tp.Dir = dir;
    tp.Ang = atan2(dir[1], dir[0]);
    if(prt) mf::LogVerbatim("CC")<<"StartTraj "<<fHits[fromHit].WireID().Wire<<":"<<(int)fHits[fromHit].PeakTime()<<" -> "<<fHits[toHit].WireID().Wire<<":"<<(int)fHits[toHit].PeakTime()<<" dir "<<dir[0]<<" "<<dir[1]<<" ang "<<tp.Ang;
    tp.AngErr = 0.5;
    tp.Chg = fHits[fromHit].Integral();
    
    // initialize other variables used by StepCrawl
    fAveHitRMS = (fHits[fromHit].RMS() + fHits[toHit].RMS()) / 2;
    tp.AveChg = (fHits[fromHit].Integral() + fHits[toHit].Integral()) / 2;
    fChgRMS = 0.3 * tp.AveChg;
    tp.ChgDiff = (tp.Chg - tp.AveChg) / fChgRMS;
    
    // Add the first seed hit. The second one will be picked up
    // by StepCrawl
    tp.Hits.push_back(fromHit);
    inClus[fromHit] = -3;
    work.Pts.push_back(tp);
    
  } // StartTraj
  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::ReverseTraj()
  {
    // reverse the trajectory
    if(work.Pts.size() == 0) return;
    std::reverse(work.Pts.begin(), work.Pts.end());
    // reverse the direction vector
    for(unsigned short ipt = 0; ipt < work.Pts.size(); ++ipt) {
      if(work.Pts[ipt].Dir[0] != 0) work.Pts[ipt].Dir[0] = -work.Pts[ipt].Dir[0];
      if(work.Pts[ipt].Dir[1] != 0) work.Pts[ipt].Dir[1] = -work.Pts[ipt].Dir[1];
      if(work.Pts[ipt].Hits.size() > 1) std::reverse(work.Pts[ipt].Hits.begin(), work.Pts[ipt].Hits.end());
    } // ipt
  }
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::StoreTraj()
  {
    // Creates a cluster using trajectory point hits
    
    // TODO: Handle the 2 point trajectory
    if(work.Pts.size() < 2) return;
    if(work.ProcCode == USHRT_MAX) return;
    if(work.Pts.size() == 0) {
      mf::LogWarning("CC")<<"StoreTraj: No points on work traj ";
      return;
    }
    
    // reverse the trajectory it's not in the proper ClusterCrawler order
    if(fStepCrawlDir > 0) ReverseTraj();
    
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
      // Fit both ends of the work trajectory to get the slope
      // errors right
      unsigned short nPts = work.Pts.size();
      unsigned short fromIndex = nPts - 1;
      unsigned short toIndex, nPtsFit;
      // fit 5 points if the trajectory has 8 or more points
      if(nPts > 7) {
        nPtsFit = 5;
      } else {
        nPtsFit = nPts / 2;
        if(nPtsFit < 3) nPtsFit = 3;
      }
      toIndex = fromIndex - nPtsFit + 1;
      // fit fewer points for poor trajectories
      if(work.Pts[toIndex].AveHitRes > 1) toIndex = fromIndex - 2;
      TrajPoint tp;
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - before";
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
      FitTrajMid(fromIndex, toIndex, tp);
      if(tp.FitChi > 900) return;
//      work.Pts[fromIndex] = tp;
      // Only use selected quantities from the fit
      work.Pts[fromIndex].Pos = tp.Pos;
      work.Pts[fromIndex].Chg = tp.Chg;
      work.Pts[fromIndex].ChgDiff = tp.ChgDiff;
      work.Pts[fromIndex].FitChi = tp.FitChi;
      work.Pts[fromIndex].NumTPsFit = tp.NumTPsFit;
      work.Pts[fromIndex].AngErr = tp.AngErr;
      // Handle a sign change in the direction
      if(std::signbit(work.Pts[fromIndex].Ang) != std::signbit(tp.Ang)) {
        work.Pts[fromIndex].Dir[0] = -tp.Dir[0];
        work.Pts[fromIndex].Dir[1] = -tp.Dir[1];
        work.Pts[fromIndex].Ang = -tp.Ang;
      } else {
        work.Pts[fromIndex].Dir = tp.Dir;
        work.Pts[fromIndex].Ang = tp.Ang;
      }
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - after";
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
      
      // fit the other end
      fromIndex = 0;
      toIndex = nPtsFit - 1;
      if(work.Pts[toIndex].AveHitRes > 1) toIndex = 2;
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - before";
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }
      FitTrajMid(fromIndex, toIndex, tp);
      if(tp.FitChi > 900) return;
//      work.Pts[fromIndex] = tp;
      // Only use selected quantities from the fit
      work.Pts[fromIndex].Pos = tp.Pos;
      work.Pts[fromIndex].Chg = tp.Chg;
      work.Pts[fromIndex].ChgDiff = tp.ChgDiff;
      work.Pts[fromIndex].FitChi = tp.FitChi;
      work.Pts[fromIndex].NumTPsFit = tp.NumTPsFit;
      work.Pts[fromIndex].AngErr = tp.AngErr;
      // Handle a sign change in the direction
      if(std::signbit(work.Pts[fromIndex].Dir[0]) != std::signbit(tp.Dir[0])) {
        work.Pts[fromIndex].Dir[0] = -tp.Dir[0];
        work.Pts[fromIndex].Dir[1] = -tp.Dir[1];
        work.Pts[fromIndex].Ang = -tp.Ang;
      } else {
        work.Pts[fromIndex].Dir = tp.Dir;
        work.Pts[fromIndex].Ang = tp.Ang;
      }
      if(prt) {
        mf::LogVerbatim("CC")<<"Storing work. Refit trajpoint - after";
        PrintTrajPoint(fromIndex, work.Pts[fromIndex]);
      }

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
    
    if(prt) mf::LogVerbatim("CC")<<"Store cluster with size "<<fcl2hits.size();
    clProcCode = work.ProcCode;
    // special stop code to indicate special crawling
    clStopCode = work.StopCode;
    // Store the cluster
    if(!TmpStore()) return;
    

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
    
    if(itj > allTraj.size()-1) return;
    
    Trajectory& tj = allTraj[itj];
    
    mf::LogVerbatim("CC")<<"Print allTraj["<<itj<<"]: ClusterIndex "<<tj.ClusterIndex<<" ProcCode "<<tj.ProcCode<<" StopCode "<<tj.StopCode<<" Vtx[0] "<<tj.Vtx[0]<<" Vtx[1] "<<tj.Vtx[1];
    
    mf::LogVerbatim("CC")<<"TRJ Ind    W:T(WSE) Tick Delta  dTick   Ang  AngErr     Q  QDiff    AveQ  FitChi Nrm  NTPF NotNr UseNr?   Qual Vtx Hits";
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) PrintTrajPoint(ii, tj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(ipt, tj.Pts[ipt]);
    }
  } // PrintAllTraj

  
  //////////////////////////////////////////
  void ClusterCrawlerAlg::PrintWork(unsigned short tPoint)
  {
    // prints one or all
    
    mf::LogVerbatim myprt("CC");
    
    unsigned short first = 0;
    unsigned short last = work.Pts.size();
    if(tPoint == USHRT_MAX) {
      myprt<<"TRJ Ind    W:T(WSE) Tick Delta  dTick   Ang  AngErr     Q  QDiff    AveQ  FitChi Nrm  NTPF NotNr UseNr?   Qual Vtx Hits \n";
      for(unsigned short ipt = first; ipt < last; ++ipt) PrintTrajPoint(ipt, work.Pts[ipt]);
    } else {
      // just print one traj point
      if(tPoint > work.Pts.size() -1) {
        mf::LogVerbatim("CC")<<"Cant print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(tPoint, work.Pts[tPoint]);
    }
  } // PrintWork

  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::PrintTrajPoint(unsigned short ipt, TrajPoint& tp)
  {
    mf::LogVerbatim myprt("CC");
    myprt<<"TRJ"<<std::fixed;
    myprt<<std::setw(4)<<ipt;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]; // W:T
    if(tp.Pos[1] < 10) myprt<<" ";
    myprt<<std::setw(6)<<(int)(tp.Pos[1]/fScaleF); // Tick
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta; // Delta
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.Delta/fScaleF; // dTick
    myprt<<std::setw(8)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(8)<<(int)tp.Chg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgDiff;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(8)<<tp.FitChi;
    bool normalFit = (std::abs(tp.Dir[1]) < 0.95);
    myprt<<std::setw(3)<<normalFit;
    myprt<<std::setw(6)<<tp.NumTPsFit;
    myprt<<std::setw(6)<<tp.NumNotNear;
    myprt<<std::setw(6)<<tp.UsedNotNearHit;
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
  void ClusterCrawlerAlg::ReleaseTrajHits(unsigned short ipt)
  {
    // release the hits attached to trajectory point itj in the work vector
    if(ipt > work.Pts.size()-1) return;
    for(unsigned short iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = 0;
  } // ReleaseTrajHits
  
  ////////////////////////////////////////////////
  void ClusterCrawlerAlg::ReleaseAllTrajHits()
  {
    unsigned short ipt, iht;
    for(ipt = 0; ipt < work.Pts.size(); ++ipt) {
      for(iht = 0; iht < work.Pts[ipt].Hits.size(); ++iht) inClus[work.Pts[ipt].Hits[iht]] = 0;
    } // ipt
  } // ReleaseAllTrajHits

  
} // namespace cluster