#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {
  
  //////////////////////////////////////////
  bool SplitAllTraj(TjStuff& tjs, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt)
  {
    // Splits the trajectory itj in the tjs.allTraj vector into two trajectories at position pos. Splits
    // the trajectory and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
    
    if(itj > tjs.allTraj.size()-1) return false;
    if(pos < tjs.allTraj[itj].EndPt[0] + 1 || pos > tjs.allTraj[itj].EndPt[1] - 1) return false;
    
    Trajectory& tj = tjs.allTraj[itj];
    
    // ensure that there will be at least 3 TPs on each trajectory
    unsigned short ipt, ii, ntp = 0;
    for(ipt = 0; ipt < pos; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) return false;
    ntp = 0;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) return false;
    
    // make a copy
    Trajectory newTj = tjs.allTraj[itj];
    newTj.ID = tjs.allTraj.size() + 1;
    
    // Leave the first section of tj in place. Re-assign the hits
    // to the new trajectory
    unsigned int iht;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(!tj.Pts[ipt].UseHit[ii]) continue;
        iht = tj.Pts[ipt].Hits[ii];
        // This shouldn't happen but check anyway
        if(tjs.inTraj[iht] != tj.ID) continue;
        tjs.inTraj[iht] = newTj.ID;
        tj.Pts[ipt].UseHit[ii] = false;
      } // ii
    } // ipt
    SetEndPoints(tjs, tj);
    if(ivx != USHRT_MAX) tj.Vtx[1] = ivx;
    tj.AlgMod[kSplitTraj] = true;
    if(prt) {
      mf::LogVerbatim("TC")<<"Splittjs.allTraj: itj "<<tj.ID<<" EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
      PrintTrajectory(tjs, tjs.allTraj[itj], USHRT_MAX);
    }
    
    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if(eraseSize > newTj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"Splittjs.allTraj: Bad erase size ";
      return false;
    }
    
    // erase the TPs at the beginning of the new trajectory
    //    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + pos + 1);
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for(ipt = 0; ipt < 3; ++ipt) {
      for(ii = 0; ii < newTj.Pts[ipt].UseHit.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
    } // ipt
    SetEndPoints(tjs, newTj);
    if(ivx != USHRT_MAX) newTj.Vtx[0] = ivx;
    newTj.AlgMod[kSplitTraj] = true;
    tjs.allTraj.push_back(newTj);
    if(prt) {
      mf::LogVerbatim("TC")<<"Splittjs.allTraj: NewTj "<<newTj.ID<<" EndPts "<<newTj.EndPt[0]<<" to "<<newTj.EndPt[1];
      PrintTrajectory(tjs, newTj, USHRT_MAX);
    }
    return true;
    
  } // SplitAllTraj
  
  //////////////////////////////////////////
  void TrajPointTrajDOCA(TjStuff& tjs, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep)
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
  void TrajTrajDOCA(Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
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
  float TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2)
  {
    // returns the separation^2 between the hit position of two trajectory points
    float dw = tp1.HitPos[0] - tp2.HitPos[0];
    float dt = tp1.HitPos[1] - tp2.HitPos[1];
    return dw * dw + dt * dt;
  } // TrajPointHitSep2
  
  //////////////////////////////////////////
  float HitSep2(TjStuff& tjs, unsigned int iht, unsigned int jht)
  {
    // returns the separation^2 between two hits in WSE units
    if(iht > tjs.fHits.size()-1 || jht > tjs.fHits.size()-1) return 1E6;
    float dw = (float)tjs.fHits[iht]->WireID().Wire - (float)tjs.fHits[jht]->WireID().Wire;
    float dt = (tjs.fHits[iht]->PeakTime() - tjs.fHits[jht]->PeakTime()) * tjs.UnitsPerTick;
    return dw * dw + dt * dt;
  } // TrajPointHitSep2
  
  //////////////////////////////////////////
  float PointTrajSep2(float wire, float time, TrajPoint const& tp)
  {
    float dw = wire - tp.Pos[0];
    float dt = time - tp.Pos[1];
    return dw * dw + dt * dt;
  }
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff& tjs, unsigned int iht, TrajPoint const& tp)
  {
    float wire = tjs.fHits[iht]->WireID().Wire;
    float time = tjs.fHits[iht]->PeakTime() * tjs.UnitsPerTick;
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff& tjs, float wire, float time, TrajPoint const& tp)
  {
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA2(TjStuff& tjs, float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach squared between a (wire, time(WSE)) point
    // and a trajectory point
    
    float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
    float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    float dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return (dw * dw + dt * dt);
    
  } // PointTrajDOCA2
  
  //////////////////////////////////////////
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y)
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
  float TrajLength(Trajectory& tj)
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
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Returns the separation distance between two trajectory points
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation
  
  //////////////////////////////////////////
  void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& Distance)
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
  
  /////////////////////////////////////////
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Calculates the angle of a line between two TPs
    float dw = tp2.Pos[0] - tp1.Pos[0];
    float dt = tp2.Pos[1] - tp1.Pos[1];
    return atan2(dw, dt);
  } // TwoTPAngle
  
  ////////////////////////////////////////////////
  void PutTrajHitsInVector(Trajectory const& tj, bool onlyUsedHits, std::vector<unsigned int>& hitVec)
  {
    // Put hits in each trajectory point into a flat vector. Only hits with UseHit if onlyUsedHits == true
    hitVec.clear();
    unsigned short ipt, ii;
    unsigned int iht;
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Hits.size() != tj.Pts[ipt].UseHit.size()) {
        mf::LogVerbatim("TC")<<"Ooops: sizes different "<<ipt<<" "<<tj.Pts[ipt].Hits.size()<<" "<<tj.Pts[ipt].UseHit.size();
        hitVec.clear();
        return;
      }
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        iht = tj.Pts[ipt].Hits[ii];
        if(onlyUsedHits) {
          if(tj.Pts[ipt].UseHit[ii]) hitVec.push_back(iht);
        } else {
          hitVec.push_back(iht);
        }
      } // iht
    } // ipt
  } // PutTrajHitsInVector
  
  //////////////////////////////////////////
  void MoveTPToWire(TrajPoint& tp, float wire)
  {
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    if(tp.Dir[0] == 0) return;
    float dw = wire - tp.Pos[0];
    if(dw == 0) return;
    tp.Pos[0] = wire;
    tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
  } // MoveTPToWire
  
  ////////////////////////////////////////////////
  void MaskTrajEndPoints(TjStuff& tjs, Trajectory& tj, unsigned short nPts)
  {
    // Mask off some trajectory points at the leading edge. Keep track of the wire we are on at the current TP.
    // The trajectory has been corrupted by these hits, so we need to first remove them
    // and then recalculate the trajectory position of the current step.
    // don't remove points but simply set UseHit false
    unsigned short ii, ipt, iht;
    for(unsigned short kill = 0; kill < nPts; ++kill) {
      ipt = tj.Pts.size() - 1 - kill;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        iht = tj.Pts[ipt].Hits[ii];
        tj.Pts[ipt].UseHit[ii] = false;
        if(tjs.inTraj[iht] == tj.ID) tjs.inTraj[iht] = 0;
      } // ii
    } // kill
    
  } // MaskTrajEndPoints
  
  //////////////////////////////////////////
  float DeltaAngle(float Ang1, float Ang2) {
    return std::abs(std::remainder(Ang1 - Ang2, M_PI));
  }
  
  ////////////////////////////////////////////////
  void SetEndPoints(TjStuff& tjs, Trajectory& tj)
  {
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
    
    tj.EndPt[0] = 0; tj.EndPt[1] = 0;
    if(tj.Pts.size() == 0) return;
    
    // check the end point pointers
    unsigned short ipt, ii;
    // make sure that the Chg is set correctly
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].UseHit.size(); ++ii)
        if(tj.Pts[ipt].UseHit[ii]) tj.Pts[ipt].Chg += tjs.fHits[tj.Pts[ipt].Hits[ii]]->Integral();
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
  
  // ****************************** Printing  ******************************
  
  void PrintAllTraj(TjStuff& tjs, DebugStuff& Debug, unsigned short itj, unsigned short ipt)
  {
    
    mf::LogVerbatim myprt("TC");
    
    if(!tjs.vtx3.empty()) {
      // print out 3D vertices
      myprt<<"****** 3D vertices ******************************************__2DVtx_Indx__*******\n";
      myprt<<"Vtx  Cstat  TPC Proc     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for(unsigned short iv = 0; iv < tjs.vtx3.size(); ++iv) {
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(7)<<tjs.vtx3[iv].CStat;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].TPC;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].ProcCode;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].X;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Y;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Z;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].XErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].YErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].ZErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[0];
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[1];
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[2];
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Wire;
        if(tjs.vtx3[iv].Wire < 0) {
          myprt<<"    Matched in all planes";
        } else {
          myprt<<"    Incomplete";
        }
        myprt<<"\n";
      }
    } // tjs.vtx3.size
    
    if(!tjs.vtx.empty()) {
      // print out 2D vertices
      myprt<<"************ 2D vertices ************\n";
      myprt<<"Vtx   CTP    wire     error   tick     error  ChiDOF  NCl  topo  traj IDs\n";
      for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
        auto const& aVtx = tjs.vtx[iv];
        if(Debug.Plane < 3 && Debug.Plane != (int)aVtx.CTP) continue;
        if(aVtx.NTraj == 0) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<aVtx.CTP;
        myprt<<std::right<<std::setw(8)<<aVtx.Wire<<" +/- ";
        myprt<<std::right<<std::setw(4)<<aVtx.WireErr;
        myprt<<std::right<<std::setw(8)<<aVtx.Time/tjs.UnitsPerTick<<" +/- ";
        myprt<<std::right<<std::setw(4)<<aVtx.TimeErr/tjs.UnitsPerTick;
        myprt<<std::right<<std::setw(8)<<aVtx.ChiDOF;
        myprt<<std::right<<std::setw(5)<<aVtx.NTraj;
        myprt<<std::right<<std::setw(6)<<aVtx.Topo;
        myprt<<"    ";
        // display the traj indices
        for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
          auto const& aTj = tjs.allTraj[ii];
          if(Debug.Plane < 3 && Debug.Plane != (int)aTj.CTP) continue;
          if(aTj.AlgMod[kKilled]) continue;
          for(unsigned short end = 0; end < 2; ++end)
            if(aTj.Vtx[end] == (short)iv) myprt<<std::right<<std::setw(4)<<aTj.ID<<"_"<<end;
        }
        myprt<<"\n";
      } // iv
    } // tjs.vtx.size
    
    // Print all trajectories in tjs.allTraj if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    unsigned short endPt;
    unsigned int iht;
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      std::vector<unsigned int> tmp;
      myprt<<"TRJ  ID CTP Pass Pts frm  to     W:Tick   Ang   AveQ     W:T      Ang   AveQ  ChgRMS Hits/TP __Vtx__ PDG Parent TRuPDG   EP   KE  \n";
      for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
        auto const& aTj = tjs.allTraj[ii];
        if(Debug.Plane >=0 && Debug.Plane < 3 && (unsigned short)Debug.Plane != aTj.CTP) continue;
        //        if(tj.AlgMod[kKilled]) continue;
        if(aTj.AlgMod[kKilled]) { myprt<<"xxx"; } else { myprt<<"TRJ"; }
        myprt<<std::fixed<<std::setw(4)<<aTj.ID;
        myprt<<std::setw(3)<<aTj.CTP;
        myprt<<std::setw(5)<<aTj.Pass;
        myprt<<std::setw(5)<<aTj.Pts.size();
        myprt<<std::setw(4)<<aTj.EndPt[0];
        myprt<<std::setw(4)<<aTj.EndPt[1];
        endPt = aTj.EndPt[0];
        TrajPoint tp = aTj.Pts[endPt];
        unsigned short itick = tp.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.AveChg;
        endPt = aTj.EndPt[1];
        tp = aTj.Pts[endPt];
        itick = tp.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(7)<<(int)tp.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<aTj.ChgRMS;
        // find average number of used hits / TP
        PutTrajHitsInVector(aTj, true, tmp);
        float ave = (float)tmp.size() / (float)aTj.Pts.size();
        myprt<<std::setw(8)<<std::setprecision(2)<<ave;
        myprt<<std::setw(4)<<aTj.Vtx[0];
        myprt<<std::setw(4)<<aTj.Vtx[1];
        myprt<<std::setw(5)<<aTj.PDG;
        myprt<<std::setw(6)<<aTj.ParentTraj;
        myprt<<std::setw(6)<<aTj.TruPDG;
        myprt<<std::setw(6)<<std::setprecision(2)<<aTj.EffPur;
        myprt<<std::setw(5)<<(int)aTj.TruKE;
        if(!aTj.Pts[0].Hits.empty()) {
          iht = aTj.Pts[0].Hits[0];
          myprt<<" "<<PrintHit(tjs.fHits[iht]);
        }
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(aTj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
      } // ii
      return;
    } // itj > tjs.allTraj.size()-1
    
    if(itj > tjs.allTraj.size()-1) return;
    
    auto const& aTj = tjs.allTraj[itj];
    
    mf::LogVerbatim("TC")<<"Print tjs.allTraj["<<itj<<"]: ClusterIndex "<<aTj.ClusterIndex<<" Vtx[0] "<<aTj.Vtx[0]<<" Vtx[1] "<<aTj.Vtx[1];
    
    PrintHeader();
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < aTj.Pts.size(); ++ii) PrintTrajPoint(tjs, ii, aTj.StepDir, aTj.Pass, aTj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(tjs, ipt, aTj.StepDir, aTj.Pass, aTj.Pts[ipt]);
    }
  } // Printtjs.allTraj
  
  
  //////////////////////////////////////////
  void PrintTrajectory(TjStuff& tjs, Trajectory const& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    unsigned short first = 0;
    unsigned short last = tj.Pts.size();
    if(tPoint == USHRT_MAX) {
      if(tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Work:    ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDG<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      } else {
        mf::LogVerbatim myprt("TC");
        myprt<<"tjs.allTraj: ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDG<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[0];
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      }
      PrintHeader();
      for(unsigned short ipt = first; ipt < last; ++ipt) PrintTrajPoint(tjs, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
    } else {
      // just print one traj point
      if(tPoint > tj.Pts.size() -1) {
        mf::LogVerbatim("TC")<<"Can't print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(tjs, tPoint, tj.StepDir, tj.Pass, tj.Pts[tPoint]);
    }
  } // PrintTrajectory
  
  //////////////////////////////////////////
  void PrintHeader()
  {
    mf::LogVerbatim("TC")<<"TRP  CTP Ind  Stp      W:Tick   Delta  RMS    Ang   Err   Kink  Dir0  Dir1      Q    AveQ Pull FitChi NTPF  Hits ";
  } // PrintHeader
  
  ////////////////////////////////////////////////
  void PrintTrajPoint(TjStuff& tjs, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<"TRP"<<std::fixed;
    myprt<<pass;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(3)<<tp.CTP;
    myprt<<std::setw(4)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(7)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]/tjs.UnitsPerTick; // W:T
    if(tp.Pos[1] < 10) myprt<<"  "; if(tp.Pos[1] < 100) myprt<<" "; if(tp.Pos[1] < 1000) myprt<<" ";
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta;
    /*
     int itick = (int)(tp.Delta/tjs.UnitsPerTick);
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
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgPull;
    myprt<<std::setw(7)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NTPsFit;
    // print the hits associated with this traj point
    for(unsigned short iht = 0; iht < tp.Hits.size(); ++iht) {
      myprt<<" "<<tjs.fHits[tp.Hits[iht]]->WireID().Plane;
      myprt<<":"<<tjs.fHits[tp.Hits[iht]]->WireID().Wire;
      if(tp.UseHit[iht]) {
        // Distinguish used hits from nearby hits
        myprt<<":";
      } else {
        myprt<<"x";
      }
      myprt<<(int)tjs.fHits[tp.Hits[iht]]->PeakTime();
      myprt<<"_"<<tjs.inTraj[tp.Hits[iht]];
    } // iht
  } // PrintTrajPoint
  
  /////////////////////////////////////////
  std::string PrintHit(const art::Ptr<recob::Hit>& hit)
  {
    return std::to_string(hit->WireID().Wire) + ":" + std::to_string((int)hit->PeakTime());
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintPos(TjStuff& tjs, TrajPoint const& tp)
  {
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    unsigned int time = std::nearbyint(tp.Pos[1]/tjs.UnitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

  
} // namespace tca

