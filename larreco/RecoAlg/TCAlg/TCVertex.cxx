#include "larreco/RecoAlg/TCAlg/TCVertex.h"

namespace tca {
  
  struct SortEntry{
    unsigned int index;
    float val;
  };
  
  bool valDecreasing (SortEntry c1, SortEntry c2) { return (c1.val > c2.val);}
  bool valIncreasing (SortEntry c1, SortEntry c2) { return (c1.val < c2.val);}
  
  //////////////////////////////////////////
  void Find2DVertices(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Find 2D vertices between pairs of tjs that have a same-end topology. Using an example
    // where StepDir = 1 (end 0 is at small wire number) vertices will be found with Topo = 0
    // with a vertex US of the ends (<) or Topo = 2 with a vertex DS of the ends (>). This is reversed
    // if StepDir = -1. Vertices with Topo = 1 (/\) and (\/) are found in EndMerge.
    
    // tjs.Vertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction
    // 7 min Score
    // 8 Min charge fraction near a merge point (not a vertex)
    // 9 max MCSmom asymmetry for a merge
    
    if(tjs.Vertex2DCuts[0] <= 0) return;
    
    if(tjs.allTraj.size() < 2) return;
    
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    bool prt = (debug.Plane == (int)planeID.Plane && debug.Tick < 0);
    if(prt) {
      mf::LogVerbatim("TC")<<"prt set for plane "<<planeID.Plane<<" in Find2DVertices";
      PrintAllTraj("F2DVi", tjs, debug, USHRT_MAX, tjs.allTraj.size());
    }
    
    // Ensure that all tjs are in the same order
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != inCTP) continue;
      if(tj.StepDir != tjs.StepDir && !tj.AlgMod[kSetDir]) ReverseTraj(tjs, tj);
    } // tj
    
    unsigned short maxShortTjLen = tjs.Vertex2DCuts[0];
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size() - 1; ++it1) {
      auto& tj1 = tjs.allTraj[it1];
      if(tj1.AlgMod[kKilled]) continue;
      if(tj1.AlgMod[kInShower]) continue;
      if(tj1.CTP != inCTP) continue;
      bool tj1Short = (TrajLength(tj1) < maxShortTjLen);
      for(unsigned short end = 0; end < 2; ++end) {
        // vertex assignment exists?
        if(tj1.VtxID[end] > 0) continue;
        // default condition is to use the end point to define the trajectory and direction
        // at the end
        short endPt1 = tj1.EndPt[end];
        float wire1 = tj1.Pts[endPt1].Pos[0];
        // unless there are few points fitted, indicating that the trajectory fit
        // may have been biased by the presence of another trajectory at the vertex or by
        // other close unresolved tracks
        if(tj1.Pts.size() > 6 && tj1.Pts[endPt1].NTPsFit < 4) {
          if(end == 0 && endPt1 < int(tj1.Pts.size()) - 3) {
            endPt1 += 3;
          } else if (end == 1 && endPt1 >=3 ) {
            endPt1 -= 3;
          }
          if(tj1.Pts[endPt1].Chg == 0) endPt1 = NearestPtWithChg(tjs, tj1, endPt1);
        } // few points fit at end1
        TrajPoint tp1 = tj1.Pts[endPt1];
        MoveTPToWire(tp1, wire1);
        // re-purpose endPt1 to reference the end point. This will be used the find the point on
        // tj1 that is closest to the vertex position
        endPt1 = tj1.EndPt[end];
        short oendPt1 = tj1.EndPt[1-end];
        for(unsigned short it2 = it1 + 1; it2 < tjs.allTraj.size(); ++it2) {
          auto& tj2 = tjs.allTraj[it2];
          if(tj2.AlgMod[kKilled]) continue;
          if(tj2.AlgMod[kInShower]) continue;
          if(tj2.CTP != inCTP) continue;
          if(tj1.VtxID[end] > 0) continue;
          if(tj2.VtxID[end] > 0) continue;
          if(tj1.MCSMom < tjs.Vertex2DCuts[5] && tj2.MCSMom < tjs.Vertex2DCuts[5]) continue;
          bool tj2Short = (TrajLength(tj2) < maxShortTjLen);
          short endPt2 = tj2.EndPt[end];
          float wire2 = tj2.Pts[endPt2].Pos[0];
          if(tj2.Pts.size() > 6 && tj2.Pts[endPt2].NTPsFit < 4) {
            if(end == 0 && endPt2 < int(tj2.Pts.size()) - 3) {
              endPt2 += 3;
            } else if (end == 1 && endPt2 >= 3){
              endPt2 -= 3;
            }
            if(tj2.Pts[endPt2].Chg == 0) endPt2 = NearestPtWithChg(tjs, tj2, endPt2);
          } // few points fit at end1
          TrajPoint tp2 = tj2.Pts[endPt2];
          MoveTPToWire(tp2, wire2);
          // re-purpose endPt2
          endPt2 = tj2.EndPt[end];
          unsigned short oendPt2 = tj2.EndPt[1-end];
          // Rough first cut on the separation between the end points of the
          // two trajectories
          float sepCut = 100;
          if(std::abs(tp1.Pos[0] - tp2.Pos[0]) > sepCut) continue;
          if(std::abs(tp1.Pos[1] - tp2.Pos[1]) > sepCut) continue;
          float wint, tint;
          TrajIntersection(tp1, tp2, wint, tint);
          // make sure this is inside the TPC
          if(wint < 0 || wint > tjs.MaxPos0[planeID.Plane]) continue;
          if(tint < 0 || tint > tjs.MaxPos1[planeID.Plane]) continue;
          // Next cut on separation between the TPs and the intersection point
          if(tj1Short || tj2Short) { sepCut = tjs.Vertex2DCuts[1]; } else { sepCut = tjs.Vertex2DCuts[2]; }
          std::array<float, 2> vPos {wint, tint};
          float vt1Sep = PosSep(vPos, tp1.Pos);
          float vt2Sep = PosSep(vPos, tp2.Pos);
          float dwc1 = DeadWireCount(tjs, wint, tp1.Pos[0], tp1.CTP);
          float dwc2 = DeadWireCount(tjs, wint, tp2.Pos[0], tp1.CTP);
          vt1Sep -= dwc1;
          vt2Sep -= dwc2;
          bool vtxOnDeadWire = (DeadWireCount(tjs, wint, wint, tp1.CTP) == 1);            
          if(prt && vt1Sep < 200 && vt2Sep < 200) {
            mf::LogVerbatim myprt("TC");
            myprt<<"F2DV candidate tj1-tj2 "<<tj1.ID<<"_"<<end<<"-"<<tj2.ID<<"_"<<end;
            myprt<<" vtx pos "<<(int)wint<<":"<<(int)(tint/tjs.UnitsPerTick)<<" tp1 "<<PrintPos(tjs, tp1)<<" tp2 "<<PrintPos(tjs, tp2);
            myprt<<" dwc1 "<<dwc1<<" dwc2 "<<dwc2<<" on dead wire? "<<vtxOnDeadWire;
            myprt<<" vt1Sep "<<vt1Sep<<" vt2Sep "<<vt2Sep<<" sepCut "<<sepCut;
          }
          if(vt1Sep > sepCut || vt2Sep > sepCut) continue;
          // make sure that the other end isn't closer
          if(PosSep(vPos, tjs.allTraj[it1].Pts[oendPt1].Pos) < vt1Sep) {
            if(prt) mf::LogVerbatim("TC")<<" tj1 other end "<<PrintPos(tjs, tj1.Pts[oendPt1])<<" is closer to the vertex";
            continue;
          }
          if(PosSep(vPos, tjs.allTraj[it2].Pts[oendPt2].Pos) < vt2Sep) {
            if(prt) mf::LogVerbatim("TC")<<" tj2 other end "<<PrintPos(tjs, tj2.Pts[oendPt2])<<" is closer to the vertex";
            continue;
          }
          // Ensure that the vertex position is close to the end of each Tj
          unsigned short closePt1;
          float doca1;
          TrajClosestApproach(tjs.allTraj[it1], wint, tint, closePt1, doca1);
          // dpt1 (and dpt2) will be 0 if the vertex 
          short dpt1 = tjs.StepDir * (closePt1 - endPt1);
          if(prt) mf::LogVerbatim("TC")<<" endPt1 "<<endPt1<<" closePt1 "<<closePt1<<" dpt1 "<<dpt1<<" doca1 "<<doca1;
          if(tjs.allTraj[it1].EndPt[1] > 4) {
            if(dpt1 > 3) continue;
          } else {
            // tighter cut for short trajectories
            if(dpt1 > 2) continue;
          }
          unsigned short closePt2;
          float doca2;
          TrajClosestApproach(tjs.allTraj[it2], wint, tint, closePt2, doca2);
          short dpt2 = tjs.StepDir * (closePt2 - endPt2);
          if(prt) mf::LogVerbatim("TC")<<" endPt2 "<<endPt2<<" closePt2 "<<closePt2<<" dpt2 "<<dpt2<<" doca2 "<<doca2;
          if(tjs.allTraj[it2].EndPt[1] > 4) {
            if(dpt2 > 3) continue;
          } else {
            // tighter cut for short trajectories
            if(dpt2 > 2) continue;
          }
          if(prt) mf::LogVerbatim("TC")<<" wint:tint "<<(int)wint<<":"<<(int)(tint/tjs.UnitsPerTick);
          // ensure that there is a signal between these TPs and the vertex on most of the wires
          short dpt = abs(wint - tp1.Pos[0]);
          if(dpt > 2 && !SignalBetween(tjs, tp1, wint, tjs.Vertex2DCuts[6], prt)) {
            if(prt) mf::LogVerbatim("TC")<<" Fails SignalBetween for tp1 "<<dpt;
            continue;
          }
          dpt = abs(wint - tp2.Pos[0]);
          if(dpt > 2 && !SignalBetween(tjs, tp2, wint, tjs.Vertex2DCuts[6], prt)) {
            if(prt) mf::LogVerbatim("TC")<<" Fails SignalBetween for tp2 "<<dpt;
            continue;
          }
          // TODO: move the vertex a bit here?
          // make a new temporary vertex
          VtxStore aVtx;
          aVtx.Pos[0] = wint;
          aVtx.Pos[1] = tint;
          aVtx.NTraj = 0;
          aVtx.Pass = tj1.Pass;
          // Topo 0 has this topology (<) and Topo 2 has this (>)
          aVtx.Topo = 2 * end;
          aVtx.ChiDOF = 0;
          aVtx.CTP = inCTP;
          // fix the vertex position if we needed to move it significantly, or if it is on a dead wire
//          if(close2 > 1) aVtx.Stat[kFixed] = true;
          // try to fit it. We need to give it an ID to do that. Take the next
          // available ID
          unsigned short newVtxID = tjs.vtx.size() + 1;
          aVtx.ID = newVtxID;
          tj1.VtxID[end] = newVtxID;
          tj2.VtxID[end] = newVtxID;
          if(!FitVertex(tjs, aVtx, prt)) {
            tj1.VtxID[end] = 0;
            tj2.VtxID[end] = 0;
            continue;
          }
          // check proximity to nearby vertices
          unsigned short mergeMeWithVx = IsCloseToVertex(tjs, aVtx);
          if(mergeMeWithVx > 0 && MergeWithVertex(tjs, aVtx, mergeMeWithVx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<" Merged with close vertex "<<mergeMeWithVx;
            continue;
          }
          // Save it
          if(!StoreVertex(tjs, aVtx)) continue;
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<" New vtx Vx_"<<aVtx.ID;
            myprt<<" Tjs "<<tj1.ID<<"_"<<end<<"-"<<tj2.ID<<"_"<<end;
            myprt<<" at "<<std::fixed<<std::setprecision(1)<<aVtx.Pos[0]<<":"<<aVtx.Pos[1]/tjs.UnitsPerTick;
          }
          SetVx2Score(tjs, prt);
        } // it2
      } // end1
    } // it1
    
    // check the consistency of the Tjs for the newly added vertices
    ChkVxTjs(tjs, inCTP, prt);
    
    // Split trajectories that cross a vertex
    SplitTrajCrossingVertices(tjs, inCTP);
    FindHammerVertices(tjs, inCTP);
    FindHammerVertices2(tjs, inCTP);
    
    if(prt) PrintAllTraj("F2DVo", tjs, debug, USHRT_MAX, USHRT_MAX);
    
  } // Find2DVertices

  
  //////////////////////////////////////////
  void MakeJunkTjVertices(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Make 2D vertices between the end of Tjs that have a junk Tj near an end
    
    if(!tjs.UseAlg[kBlobVx]) return;
    
    if(tjs.allTraj.empty()) return;
    
    unsigned short plane = DecodeCTP(inCTP).Plane;
    bool prt = (debug.Plane >= 0 && debug.Tick == 99999);
    if(prt) {
      mf::LogVerbatim("TC")<<"MakeJunkTjVertices: prt set for plane "<<plane;
    }

    float window = 3;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kShowerTj]) continue;
      if(tj.Pts.size() < 10) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] > 0) continue;
        auto tp = tj.Pts[tj.EndPt[end]];
        auto closeTjs = FindCloseTjs(tjs, tp, tp, window);
        // make a subset of junk tj candidates
        std::vector<int> junkTjs;
        for(auto tjid : closeTjs) {
          auto& ctj = tjs.allTraj[tjid - 1];
          if(ctj.ID == tj.ID) continue;
          if(!tj.AlgMod[kJunkTj]) continue;
          junkTjs.push_back(tjid);
        } // tjid
        if(junkTjs.empty()) continue;
        std::cout<<"MakeJunkTjVertices:";
        for(auto tjid : junkTjs) std::cout<<" "<<tjid;
        std::cout<<"\n";
      } // end
    } //  tj
  } // MakeJunkTjVertices
  
  //////////////////////////////////////////
  bool MergeWithVertex(TjStuff& tjs, VtxStore& vx, unsigned short oVxID, bool prt)
  {
    // Attempts to merge the trajectories attached to vx with an existing 2D vertex
    // referenced by existingVxID. This function doesn't use the existing end0/end1 vertex association.
    // It returns true if the merging was successful in which case the calling function should
    // not store vx. The calling function needs to have set VtxID to vx.ID for tjs that are currently attached
    // to vx. It assumed that vx hasn't yet been pushed onto tjs.vtx
    
    if(!tjs.UseAlg[kVxMerge]) return false;
    
    if(oVxID > tjs.vtx.size()) return false;
    auto& oVx = tjs.vtx[oVxID - 1];
    if(vx.CTP != oVx.CTP) return false;
    
    // get a list of tjs attached to both vertices
    std::vector<int> tjlist = GetVtxTjIDs(tjs, vx);
    if(tjlist.empty()) return false;
    std::vector<int> tmp = GetVtxTjIDs(tjs, oVx);
    if(tmp.empty()) return false;
    for(auto tjid : tmp) {
      if(std::find(tjlist.begin(), tjlist.end(), tjid) == tjlist.end()) tjlist.push_back(tjid);
    } // tjid
    if(tjlist.size() < 2) return false;
    // handle the simple case
    if(tjlist.size() == 2) {
      // Unset the fixed bit
      vx.Stat[kFixed] = false;
      oVx.Stat[kFixed] = false;
      // assign the vx tjs to oVx
      for(auto tjid : tjlist) {
        auto& tj = tjs.allTraj[tjid - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == vx.ID) tj.VtxID[end] = oVx.ID;
        } // end
      } // tjid
      if(!FitVertex(tjs, oVx, prt)) {
        if(prt) mf::LogVerbatim("TC")<<"MWV: merge failed "<<vx.ID<<" and existing "<<oVx.ID;
        return false;
      }
      return true;
    } // size = 2
    
    // sort by decreasing length
    std::vector<SortEntry> sortVec(tjlist.size());
    for(unsigned int indx = 0; indx < sortVec.size(); ++indx) {
      sortVec[indx].index = indx;
      auto& tj = tjs.allTraj[tjlist[indx] - 1];
      sortVec[indx].val = tj.Pts.size();
    } // indx
    std::sort(sortVec.begin(), sortVec.end(), valDecreasing);
    // re-order the list of Tjs
    auto ttl = tjlist;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) tjlist[ii] = ttl[sortVec[ii].index];
    // Create a local vertex using the two longest tjs, then add the shorter ones
    // until the pull reaches the cut
    VtxStore aVx;
    aVx.CTP = vx.CTP;
    std::vector<TrajPoint> tjpts(tjlist.size());
    // determine which point on each Tj that will be used in the vertex fit and stash it in
    // the traj point Step variable. This requires knowing the real position of the merged vertex
    // which we estimate by averaging
    std::array<float, 2> vpos;
    vpos[0] = 0.5 * (vx.Pos[0] + oVx.Pos[0]);
    vpos[1] = 0.5 * (vx.Pos[1] + oVx.Pos[1]);
    for(unsigned short ii = 0; ii < tjpts.size(); ++ii) {
      auto& tj = tjs.allTraj[tjlist[ii] - 1];
      unsigned short end = CloseEnd(tjs, tj, vpos);
      // assume that we will use the end point of the tj
      unsigned short endPt = tj.EndPt[end];
      if(tj.Pts.size() > 6 && tj.Pts[endPt].NTPsFit < 4) {
        if(end == 0) {
          endPt += 3;
        } else {
          endPt -= 3;
        }
      } // few points fit at the end
      // define tjpts
      tjpts[ii].CTP = tj.CTP;
      tjpts[ii].Pos = tj.Pts[endPt].Pos;
      tjpts[ii].Dir = tj.Pts[endPt].Dir;
      tjpts[ii].Ang = tj.Pts[endPt].Ang;
      tjpts[ii].AngErr = tj.Pts[endPt].AngErr;
      // stash the point in Step
      tjpts[ii].Step = endPt;
      // and the end in AngleCode
      tjpts[ii].AngleCode = end;
    } // tjid
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MWV: "<<oVxID;
      myprt<<" Fit TPs";
      for(unsigned short ii = 0; ii < tjpts.size(); ++ii) {
        auto& tjpt = tjpts[ii];
        myprt<<" "<<tjlist[ii]<<"_"<<tjpt.Step<<"_"<<PrintPos(tjs, tjpt.Pos);
      }
    } // prt
    // create a subset of the first two for the first fit
    auto fitpts = tjpts;
    fitpts.resize(2);
    std::vector<int> fittjs(2);
    fittjs[0] = tjlist[0];
    fittjs[1] = tjlist[1];
    if(!FitVertex(tjs, aVx, fitpts, prt)) {
      if(prt) mf::LogVerbatim("TC")<<"MWV: first fit failed ";
      return false;
    }
    // Fit and add tjs to the vertex
    bool needsUpdate = false;
    for(unsigned short ii = 2; ii < tjlist.size(); ++ii) {
      fitpts.push_back(tjpts[ii]);
      fittjs.push_back(tjlist[ii]);
      if(FitVertex(tjs, aVx, fitpts, prt)) {
        needsUpdate = false;
      } else {
        // remove the last Tj point and keep going
        fitpts.pop_back();
        fittjs.pop_back();
        needsUpdate = true;
      }
    } // ii
    
    if(needsUpdate) FitVertex(tjs, aVx, fitpts, prt);
    if(prt) mf::LogVerbatim("TC")<<"MWV: done "<<vx.ID<<" and existing "<<oVx.ID;
    
    // update. Remove old associations
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx.CTP) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == vx.ID) tj.VtxID[end] = 0;
        if(tj.VtxID[end] == oVxID) tj.VtxID[end] = 0;
      }
    } // tj
    // set the new associations
    for(unsigned short ii = 0; ii < fitpts.size(); ++ii) {
      auto& tjpt = fitpts[ii];
      unsigned short end = tjpt.AngleCode;
      auto& tj = tjs.allTraj[fittjs[ii] - 1];
      if(tj.VtxID[end] != 0) {
        std::cout<<"MWV: coding error. tj "<<tj.ID<<" end "<<end<<" VtxID "<<tj.VtxID[end]<<" != 0\n";
        return false;
      }
      tj.VtxID[end] = oVxID;
    } // ii
    
    // Update oVx 
    oVx.Pos = aVx.Pos;
    oVx.PosErr = aVx.PosErr;
    oVx.ChiDOF = aVx.ChiDOF;
    oVx.NTraj = fitpts.size();
    // Update the score and the charge fraction
    SetVx2Score(tjs, oVx, prt);
    oVx.Stat[kVtxMerged] = true;
    oVx.Stat[kFixed] = false;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MWV: "<<oVxID;
      myprt<<" Done TPs";
      for(unsigned short ii = 0; ii < fitpts.size(); ++ii) {
        auto& tjpt = fitpts[ii];
        myprt<<" "<<fittjs[ii]<<"_"<<tjpt.AngleCode<<"_"<<PrintPos(tjs, tjpt.Pos);
      }
    } // prt

    return true;
  } // MergeWithVertex
    
  //////////////////////////////////////////
  void ChkVxTjs(TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // 
    
    if(!tjs.UseAlg[kChkVxTj]) return;
    
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      auto& vx2 = tjs.vtx[ivx];
      if(vx2.ID == 0) continue;
      if(vx2.CTP != inCTP) continue;
      auto vxtjs = GetVtxTjIDs(tjs, vx2);
      if(vxtjs.size() < 2) continue;
      for(unsigned short it1 = 0; it1 < vxtjs.size() - 1; ++it1) {
        auto& tj1 = tjs.allTraj[vxtjs[it1] - 1];
        if(tj1.AlgMod[kKilled]) continue;
        unsigned short end1 = 0;
        if(tj1.VtxID[1] == vx2.ID) end1 = 1;
        auto& vtp1 = tj1.Pts[tj1.EndPt[end1]];
        auto& otp1 = tj1.Pts[tj1.EndPt[1 - end1]];
        float tj1sep = PosSep(vtp1.Pos, vx2.Pos);
        for(unsigned short it2 = it1 + 1; it2 < vxtjs.size(); ++it2) {
          auto& tj2 = tjs.allTraj[vxtjs[it2] - 1];
          if(tj2.AlgMod[kKilled]) continue;
          unsigned short end2 = 0;
          if(tj2.VtxID[2] == vx2.ID) end2 = 1;
          auto& vtp2 = tj2.Pts[tj2.EndPt[end2]];
          auto& otp2 = tj2.Pts[tj2.EndPt[1 - end2]];
          float tj2sep = PosSep(vtp2.Pos, vx2.Pos);
          float otj1tj2 = PosSep(otp1.Pos, vtp2.Pos);
          float delta12 = PointTrajDOCA(tjs, otp1.Pos[0], otp1.Pos[1], vtp2);
          float dang12 = DeltaAngle(otp1.Ang, vtp2.Ang);
          if(otj1tj2 < tj2sep && delta12 < 1 && otj1tj2 < 4) {
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<"CVTjs: "<<vx2.ID<<" tj1 "<<tj1.ID<<" tj2 "<<tj2.ID;
              myprt<<" otj1tj2 "<<otj1tj2;
              myprt<<" delta12 "<<delta12;
              myprt<<" dang12 "<<dang12;
              myprt<<" Try to merge";
            }
            // End 1 of tj1 is closer to end0 of tj2 than tj2 is to the vertex
            tj2.VtxID[end2] = 0;
            if(CompatibleMerge(tjs, tj1, tj2, prt) && MergeAndStore(tjs, vxtjs[it1], vxtjs[it2], prt)) {
              auto& newTj = tjs.allTraj[tjs.allTraj.size()-1];
              newTj.AlgMod[kChkVxTj] = true;
              if(prt) mf::LogVerbatim("TC")<<"CVTjs: Merged tjs "<<tj1.ID<<" and "<<tj2.ID<<" -> "<<newTj.ID;
            } else {
              if(prt) mf::LogVerbatim("TC")<<"CVTjs: Merge failed";
            }
            continue;
          } // other end is closer
          // now check the other end of tj2
          float tj1otj2 = PosSep(vtp1.Pos, otp2.Pos);
          if(tj1otj2 < tj1sep && delta12 < 1 && tj1otj2 < 4) {
            // End 1 of tj1 is closer to end0 of tj2 than tj2 is to the vertex
            tj1.VtxID[end1] = 0;
            if(CompatibleMerge(tjs, tj2, tj1, prt) && MergeAndStore(tjs, vxtjs[it2], vxtjs[it1], prt)) {
              auto& newTj = tjs.allTraj[tjs.allTraj.size()-1];
              newTj.AlgMod[kChkVxTj] = true;
              if(prt) mf::LogVerbatim("TC")<<"CVTjs: Merged tjs "<<tj1.ID<<" and "<<tj2.ID<<" -> "<<newTj.ID;
            } else {
              if(prt) mf::LogVerbatim("TC")<<"CVTjs: Merge failed";
            }
          } // the other end is closer
        } // it2
      } // it1
      // Check for delta-rays that have a vertex when they should have been merged
      if(vx2.Topo == 1 && vxtjs.size() == 2) {
        auto& tj1 = tjs.allTraj[vxtjs[0] - 1];
        auto& tj2 = tjs.allTraj[vxtjs[1] - 1];
        // ensure that these weren't killed above
        if(tj1.AlgMod[kKilled] || tj2.AlgMod[kKilled]) continue;
        if(tj1.AlgMod[kDeltaRay] || tj2.AlgMod[kDeltaRay]) {
          if(prt) mf::LogVerbatim("TC")<<"CVTjs: Merge delta rays "<<tj1.ID<<" and "<<tj2.ID<<" CompatibleMerge? "<<CompatibleMerge(tjs, tj1, tj2, prt);
          MakeVertexObsolete(tjs, vx2, true);
          MergeAndStore(tjs, vxtjs[0], vxtjs[1], prt);
        } // one is a tagged delta-ray
      } // delta-ray check
    } // ivx
  } // ChkVxTjs

  //////////////////////////////////////////
  void FindHammerVertices2(TjStuff& tjs, const CTP_t& inCTP)
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
    
    if(!tjs.UseAlg[kHamVx2]) return;
    
    bool prt = (debug.Plane >= 0 && debug.Tick == 66666);
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      if(tjs.allTraj[it1].CTP != inCTP) continue;
      if(tjs.allTraj[it1].AlgMod[kKilled]) continue;
      if(tjs.allTraj[it1].AlgMod[kHamVx]) continue;
      if(tjs.allTraj[it1].AlgMod[kHamVx2]) continue;
      unsigned short numPtsWithCharge1 = NumPtsWithCharge(tjs, tjs.allTraj[it1], false);
      if(numPtsWithCharge1 < 6) continue;
      if(prt) mf::LogVerbatim("TC")<<"FindHammerVertices2: tj1 "<<tjs.allTraj[it1].ID<<" "<<tjs.allTraj[it1].MCSMom;
      // Check each end of tj1
      bool didaSplit = false;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tjs.allTraj[it1].VtxID[end1] > 0) continue;
        unsigned short endPt1 = tjs.allTraj[it1].EndPt[end1];
        for(unsigned short it2 = 0; it2 < tjs.allTraj.size(); ++it2) {
          if(it1 == it2) continue;
          if(tjs.allTraj[it2].AlgMod[kKilled]) continue;
          if(tjs.allTraj[it2].AlgMod[kHamVx]) continue;
          if(tjs.allTraj[it2].AlgMod[kHamVx2]) continue;
          // require that both be in the same CTP
          if(tjs.allTraj[it2].CTP != inCTP) continue;
          // Don't mess with muons
          if(tjs.allTraj[it2].PDGCode == 13) continue;
          unsigned short numPtsWithCharge2 = NumPtsWithCharge(tjs, tjs.allTraj[it2], true);
          if(numPtsWithCharge2 < 6) continue;
          // ignore if tj1 is a lot shorter than tj2
          // ignore if ChgRMS isn't known
          if(tjs.allTraj[it2].ChgRMS == 0) continue;
          if(numPtsWithCharge1 < 0.2 * numPtsWithCharge2) continue;
          // Find the minimum separation between tj1 and tj2
          float minDOCA = 5;
          float doca = minDOCA;
          unsigned short closePt2 = 0;
          TrajPointTrajDOCA(tjs, tjs.allTraj[it1].Pts[endPt1], tjs.allTraj[it2], closePt2, doca);
          if(doca == minDOCA) continue;
          if(prt) mf::LogVerbatim("TC")<<" FHV2 CTP "<<tjs.allTraj[it1].CTP<<" tj1 ID "<<tjs.allTraj[it1].ID<<"_"<<end1<<" tj2 ID "<<tjs.allTraj[it2].ID<<" doca "<<doca<<" tj2.EndPt[0] "<<tjs.allTraj[it2].EndPt[0]<<" closePt2 "<<closePt2<<" tj2.EndPt[1] "<<tjs.allTraj[it2].EndPt[1];
          // ensure that the closest point is not near an end
          if(closePt2 < tjs.allTraj[it2].EndPt[0] + 3) continue;
          if(closePt2 > tjs.allTraj[it2].EndPt[1] - 3) continue;
          // Find the intersection point between the tj1 end and tj2 closest Pt
          float wint, tint;
          TrajIntersection(tjs.allTraj[it1].Pts[endPt1], tjs.allTraj[it2].Pts[closePt2], wint, tint);
          // make an angle cut
          float dang = DeltaAngle(tjs.allTraj[it1].Pts[endPt1].Ang, tjs.allTraj[it2].Pts[closePt2].Ang);
          if(dang < 0.2) continue;
          // ensure that tj1 doesn't cross tj2 but ensuring that the closest point on tj1 is at closePt1
          doca = 5;
          unsigned short closePt1 = 0;
          TrajPointTrajDOCA(tjs, tjs.allTraj[it2].Pts[closePt2], tjs.allTraj[it1], closePt1, doca);
          if(closePt1 != endPt1) continue;
          if(prt) mf::LogVerbatim("TC")<<" intersection W:T "<<(int)wint<<":"<<(int)(tint/tjs.UnitsPerTick)<<" dang "<<dang;
          // Find the point on tj2 that is closest to this point, overwriting closePt
          doca = minDOCA;
          // the point on tj2 that is closest to the intersection point
          unsigned short intPt2;
          TrajClosestApproach(tjs.allTraj[it2], wint, tint, intPt2, doca);
          if(prt) mf::LogVerbatim("TC")<<" intPt2 on tj2 "<<intPt2<<" at Pos "<<PrintPos(tjs, tjs.allTraj[it2].Pts[intPt2])<<" doca "<<doca;
          if(doca == minDOCA) continue;
          // start scanning for the point on tj2 that has the best IP with the end of tj1 in the direction
          // from closePt2 -> endPt2
          short dir = 1;
          if(intPt2 < closePt2) dir = -1;
          unsigned short nit = 0;
          unsigned short ipt = intPt2;
          float mostChg = tjs.allTraj[it2].Pts[ipt].Chg;
          if(prt) mf::LogVerbatim("TC")<<" ipt "<<ipt<<" at Pos "<<PrintPos(tjs, tjs.allTraj[it2].Pts[ipt])<<" chg "<<mostChg;
          while(nit < 20) {
            ipt += dir;
            if(ipt < 3 || ipt > tjs.allTraj[it2].Pts.size() - 4) break;
            float delta = PointTrajDOCA(tjs, tjs.allTraj[it2].Pts[ipt].Pos[0], tjs.allTraj[it2].Pts[ipt].Pos[1], tjs.allTraj[it1].Pts[endPt1]);
            float sep = PosSep(tjs.allTraj[it2].Pts[ipt].Pos, tjs.allTraj[it1].Pts[endPt1].Pos);
            float dang = delta / sep;
            float pull = dang / tjs.allTraj[it1].Pts[endPt1].DeltaRMS;
            //            if(prt) mf::LogVerbatim("TC")<<" intPt2 "<<intPt2<<" at Pos "<<PrintPos(tjs, tjs.allTraj[it2].Pts[ipt])<<" delta "<<delta<<" chg "<<tjs.allTraj[it2].Pts[ipt].Chg<<" pull "<<pull;
            if(pull < 2 && tjs.allTraj[it2].Pts[ipt].Chg > mostChg) {
              mostChg = tjs.allTraj[it2].Pts[ipt].Chg;
              intPt2 = ipt;
            }
          }
          // require a lot of charge in tj2 in this vicinity compared with the average.
          float chgPull = (mostChg - tjs.allTraj[it2].AveChg) / tjs.allTraj[it2].ChgRMS;
          if(prt) mf::LogVerbatim("TC")<<" chgPull at intersection point "<<chgPull;
          if(chgPull < 10) continue;
          // require a signal on every wire between it1 and intPt2
          TrajPoint ltp = tjs.allTraj[it1].Pts[endPt1];
          unsigned int fromWire = std::nearbyint(tjs.allTraj[it1].Pts[endPt1].Pos[0]);
          unsigned int toWire =   std::nearbyint(tjs.allTraj[it2].Pts[intPt2].Pos[0]);
          if(fromWire > toWire) {
            unsigned int tmp = fromWire;
            fromWire = toWire;
            toWire = tmp;
          }
          bool skipIt = false;
          for(unsigned int wire = fromWire + 1; wire < toWire; ++wire) {
            MoveTPToWire(ltp, (float)wire);
            if(!SignalAtTp(tjs, ltp)) {
              skipIt = true;
              break;
            }
          } // wire
          if(skipIt) continue;
          // we have a winner
          // create a new vertex
          VtxStore aVtx;
          aVtx.Pos = tjs.allTraj[it2].Pts[intPt2].Pos;
          aVtx.NTraj = 3;
          aVtx.Pass = tjs.allTraj[it2].Pass;
          aVtx.Topo = 6;
          aVtx.ChiDOF = 0;
          aVtx.CTP = inCTP;
          aVtx.ID = tjs.vtx.size() + 1;
          unsigned short ivx = tjs.vtx.size();
          if(!StoreVertex(tjs, aVtx)) continue;
          if(!SplitAllTraj(tjs, it2, intPt2, ivx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<"FHV2: Failed to split trajectory";
            // we can just remove the vertex since no Tj VtxID association has been made yet
            tjs.vtx.pop_back();
            continue;
          }
          tjs.allTraj[it1].VtxID[end1] = tjs.vtx[ivx].ID;
          tjs.allTraj[it1].AlgMod[kHamVx2] = true;
          tjs.allTraj[it2].AlgMod[kHamVx2] = true;
          unsigned short newTjIndex = tjs.allTraj.size() - 1;
          tjs.allTraj[newTjIndex].AlgMod[kHamVx2] = true;
          SetVx2Score(tjs, prt);
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(tjs, it2);
          // and for the new trajectory
          SetPDGCode(tjs, newTjIndex);
          if(prt) mf::LogVerbatim("TC")<<" New vtx Vx_"<<tjs.vtx[ivx].ID;
          didaSplit = true;
          break;
        } // it2
        if(didaSplit) break;
      } // end1
    } // it1
  } // FindHammerVertices2
  
  //////////////////////////////////////////
  void FindHammerVertices(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Look for a trajectory that intersects another. Split
    // the trajectory and make a vertex. The convention used
    // is shown pictorially here. Trajectory tj1 must be longer
    // than tj2
    // tj2       ------
    // tj1         /
    // tj1        /
    // tj1       /
    
    if(!tjs.UseAlg[kHamVx]) return;
    
    bool prt = (debug.Plane >= 0 && debug.Tick == 55555);

    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      if(tjs.allTraj[it1].CTP != inCTP) continue;
      if(tjs.allTraj[it1].AlgMod[kKilled]) continue;
      // minimum length requirements
      unsigned short tj1len = tjs.allTraj[it1].EndPt[1] - tjs.allTraj[it1].EndPt[0];
      if(tj1len < 6) continue;
      // Check each end of tj1
      bool didaSplit = false;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tjs.allTraj[it1].VtxID[end1] > 0) continue;
        unsigned short endPt1 = tjs.allTraj[it1].EndPt[end1];
        for(unsigned short it2 = 0; it2 < tjs.allTraj.size(); ++it2) {
          if(tjs.allTraj[it2].CTP != inCTP) continue;
          if(it1 == it2) continue;
          if(tjs.allTraj[it2].AlgMod[kKilled]) continue;
          // length of tj2 cut
          unsigned short tj2len = tjs.allTraj[it2].EndPt[1] - tjs.allTraj[it2].EndPt[0];
          if(tj2len < 6) continue;
          // ignore if tj1 is a lot shorter than tj2
          if(tj1len < 0.5 * tj2len) continue;
          // ignore very long straight trajectories (probably a cosmic muon)
          unsigned short end20 = tjs.allTraj[it2].EndPt[0];
          unsigned short end21 = tjs.allTraj[it2].EndPt[1];
          if(tj2len > 100 && DeltaAngle(tjs.allTraj[it2].Pts[end20].Ang, tjs.allTraj[it2].Pts[end21].Ang) < 0.2) continue;
          // Require no vertex associated with itj2
          if(tjs.allTraj[it2].VtxID[0] > 0 || tjs.allTraj[it2].VtxID[1] > 0) continue;
          float minDOCA = 3;
          float doca = minDOCA;
          unsigned short closePt2 = 0;
          TrajPointTrajDOCA(tjs, tjs.allTraj[it1].Pts[endPt1], tjs.allTraj[it2], closePt2, doca);
          if(doca == minDOCA) continue;
          // ensure that the closest point is not near an end
          if(prt) mf::LogVerbatim("TC")<<"FindHammerVertices: Candidate "<<tjs.allTraj[it1].ID<<"  "<<tjs.allTraj[it2].ID<<" doca "<<doca<<" tj2.EndPt[0] "<<tjs.allTraj[it2].EndPt[0]<<" closePt2 "<<closePt2<<" tj2.EndPt[1] "<<tjs.allTraj[it2].EndPt[1];
          if(closePt2 < tjs.allTraj[it2].EndPt[0] + 3) continue;
          if(closePt2 > tjs.allTraj[it2].EndPt[1] - 3) continue;
          // make an angle cut
          float dang = DeltaAngle(tjs.allTraj[it1].Pts[endPt1].Ang, tjs.allTraj[it2].Pts[closePt2].Ang);
          if(prt) mf::LogVerbatim("TC")<<" dang "<<dang<<" imposing a hard cut of 0.4 for now ";
          if(dang < 0.4) continue;
          // we have a winner
          // create a new vertex
          VtxStore aVtx;
          aVtx.Pos = tjs.allTraj[it2].Pts[closePt2].Pos;
          aVtx.NTraj = 3;
          aVtx.Pass = tjs.allTraj[it2].Pass;
          aVtx.Topo = 5;
          aVtx.ChiDOF = 0;
          aVtx.CTP = inCTP;
          aVtx.ID = tjs.vtx.size() + 1;
          unsigned short ivx = tjs.vtx.size();
          if(!StoreVertex(tjs, aVtx)) continue;
          if(!SplitAllTraj(tjs, it2, closePt2, ivx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<"FindHammerVertices: Failed to split trajectory";
            tjs.vtx.pop_back();
            continue;
          }
          tjs.allTraj[it1].VtxID[end1] = tjs.vtx[ivx].ID;
          tjs.allTraj[it1].AlgMod[kHamVx] = true;
          tjs.allTraj[it2].AlgMod[kHamVx] = true;
          unsigned short newTjIndex = tjs.allTraj.size() - 1;
          tjs.allTraj[newTjIndex].AlgMod[kHamVx] = true;
          SetVx2Score(tjs, prt);
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(tjs, it2);
          // and for the new trajectory
          SetPDGCode(tjs, newTjIndex);
          didaSplit = true;
          break;
        } // tj2
        if(didaSplit) break;
      } // end1
    } // tj1
    
  } // FindHammerVertices

  //////////////////////////////////////////
  void SplitTrajCrossingVertices(TjStuff& tjs, const CTP_t& inCTP)
  {
    // This is kind of self-explanatory...

    if(!tjs.UseAlg[kSplitTjCVx]) return;

    if(tjs.vtx.empty()) return;
    if(tjs.allTraj.empty()) return;

    bool prt = (debug.Plane >= 0 && debug.Tick == 77777);

    geo::PlaneID planeID = DecodeCTP(inCTP);        

    // Splits trajectories in tjs.allTraj that cross a vertex
    unsigned short itj, iv, nTraj = tjs.allTraj.size();
    unsigned short tPass, closePt;
    float doca;
    for(itj = 0; itj < nTraj; ++itj) {
      if(tjs.allTraj[itj].CTP != inCTP) continue;
      // obsolete trajectory
      if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
      if(tjs.allTraj[itj].AlgMod[kSplitTjCVx]) continue;
      // too short
      if(tjs.allTraj[itj].EndPt[1] < 6) continue;
      tPass = tjs.allTraj[itj].Pass;
      for(iv = 0; iv < tjs.vtx.size(); ++iv) {
        // obsolete vertex
        if(tjs.vtx[iv].NTraj == 0) continue;
        // trajectory already associated with vertex
        if(tjs.allTraj[itj].VtxID[0] == tjs.vtx[iv].ID ||
           tjs.allTraj[itj].VtxID[1] == tjs.vtx[iv].ID) continue;
        // not in the cryostat/tpc/plane
        if(tjs.allTraj[itj].CTP != tjs.vtx[iv].CTP) continue;
        TrajClosestApproach(tjs.allTraj[itj], tjs.vtx[iv].Pos[0], tjs.vtx[iv].Pos[1], closePt, doca);
        if(prt)  mf::LogVerbatim("TC")<<" doca "<<doca<<" btw traj "<<tjs.allTraj[itj].ID<<" and tjs.vtx "<<tjs.vtx[iv].ID<<" closePt "<<closePt<<" in plane "<<planeID.Plane<<" CTP "<<tjs.vtx[iv].CTP;
        //if(doca > fMaxVertexTrajSep[tPass]) continue;
        if(doca > tjs.Vertex2DCuts[1]) continue;
        // compare the length of the Tjs used to make the vertex with the length of the
        // Tj that we want to split. Don't allow a vertex using very short Tjs to split a long
        // Tj in the 3rd plane
        auto vxtjs = GetVtxTjIDs(tjs, tjs.vtx[iv]);
        unsigned short maxPts = 0;
        for(auto tjid : vxtjs) {
          auto& tj = tjs.allTraj[tjid - 1];
          unsigned short npwc = NumPtsWithCharge(tjs, tj, false);
          if(npwc > maxPts) maxPts = npwc;
        } // tjid
        // skip this operation if any of the Tjs in the split list are > 3 * maxPts
        maxPts *= 3;
        bool skipit = false;
        if(NumPtsWithCharge(tjs,tjs.allTraj[itj] , false) > maxPts && maxPts < 100) skipit = true;
        if(prt) mf::LogVerbatim("TC")<<"  maxPts "<<maxPts<<" vxtjs[0] "<<vxtjs[0]<<" skipit? "<<skipit;
        if(skipit) continue;

        // improve closePt based on vertex position
        // check if closePt and EndPt[1] are the two sides of vertex
        // take dot product of closePt-vtx and EndPt[1]-vtx
        if ((tjs.allTraj[itj].Pts[closePt].Pos[0]-tjs.vtx[iv].Pos[0])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[1]].Pos[0]-tjs.vtx[iv].Pos[0]) + (tjs.allTraj[itj].Pts[closePt].Pos[1]-tjs.vtx[iv].Pos[1])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[1]].Pos[1]-tjs.vtx[iv].Pos[1]) <0 && closePt < tjs.allTraj[itj].EndPt[1] - 1) ++closePt;
        else if ((tjs.allTraj[itj].Pts[closePt].Pos[0]-tjs.vtx[iv].Pos[0])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[0]].Pos[0]-tjs.vtx[iv].Pos[0]) + (tjs.allTraj[itj].Pts[closePt].Pos[1]-tjs.vtx[iv].Pos[1])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[0]].Pos[1]-tjs.vtx[iv].Pos[1]) <0 && closePt > tjs.allTraj[itj].EndPt[0] + 1) --closePt;
        
        if(prt)  {
          mf::LogVerbatim("TC")<<"Good doca "<<doca<<" btw traj "<<itj<<" and tjs.vtx "<<iv<<" closePt "<<closePt<<" in plane "<<planeID.Plane<<" CTP "<<tjs.vtx[iv].CTP;
          PrintTrajPoint("STCV", tjs, closePt, 1, tPass, tjs.allTraj[itj].Pts[closePt]);
        }
        // ensure that the closest point is not near an end
        if(closePt < tjs.allTraj[itj].EndPt[0] + 3) continue;
        if(closePt > tjs.allTraj[itj].EndPt[1] - 3) continue;
        if(!SplitAllTraj(tjs, itj, closePt, iv, prt)) {
          if(prt) mf::LogVerbatim("TC")<<"SplitTrajCrossingVertices: Failed to split trajectory";
          continue;
        }
        tjs.allTraj[itj].AlgMod[kSplitTjCVx] = true;
        unsigned short newTjIndex = tjs.allTraj.size() - 1;
        tjs.allTraj[newTjIndex].AlgMod[kSplitTjCVx] = true;
      } // iv
    } // itj
    
  } // SplitTrajCrossingVertices

  //////////////////////////////////////
  void Find3DVertices(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // Create 3D vertices from 2D vertices. 3D vertices that are matched
    // in all three planes have Vtx2ID > 0 for all planes. This function re-scores all
    // 2D and 3D vertices and flags Tjs that have high-score 3D vertices    
    
    if(tjs.Vertex3DCuts[0] < 0) return;
    if(tjs.vtx.size() < 2) return;
    
    const unsigned int cstat = tpcid.Cryostat;
    const unsigned int tpc = tpcid.TPC;
    
    // create a array/vector of 2D vertex indices in each plane
    std::vector<std::vector<unsigned short>> vIndex(3);
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      // obsolete vertex
      if(tjs.vtx[ivx].ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(tjs.vtx[ivx].CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      unsigned short plane = planeID.Plane;
      if(plane > 2) continue;
      vIndex[plane].push_back(ivx);
    }
    
    unsigned short vtxInPln = 0;
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) if(vIndex[plane].size() > 0) ++vtxInPln;
    if(vtxInPln < 2) return;
    
    bool prt = (debug.Plane >= 0) && (debug.Tick == 2222);
    
    if(prt) {
      mf::LogVerbatim("TC")<<"Inside Find3DVertices";
      PrintAllTraj("F3DV", tjs, debug, USHRT_MAX, tjs.allTraj.size());
    }
    
    // wire spacing in cm
    float wirePitch = tjs.geom->WirePitch(0, 1, 0, tpcid.TPC, tpcid.Cryostat);
    
    size_t vsize = tjs.vtx.size();
    // vector of 2D vertices -> 3D vertices.
    std::vector<short> vPtr(vsize, -1);
    // fill temp vectors of 2D vertex X and X errors
    std::vector<float> vX(vsize, -100);
//    std::vector<float> vXsigma(vsize);
    
    for(unsigned short ivx = 0; ivx < vsize; ++ivx) {
      if(tjs.vtx[ivx].ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(tjs.vtx[ivx].CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      int plane = planeID.Plane;
      unsigned int wire = std::nearbyint(tjs.vtx[ivx].Pos[0]);
      if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, plane, wire))) continue;
      // Convert 2D vertex time error to X error
      double ticks = tjs.vtx[ivx].Pos[1] / tjs.UnitsPerTick;
      vX[ivx]  = tjs.detprop->ConvertTicksToX(ticks, plane, (int)tpc, (int)cstat);
/*
      ticks = (tjs.vtx[ivx].Pos[1] + tjs.vtx[ivx].PosErr[1]) / tjs.UnitsPerTick;
      vXp = tjs.detprop->ConvertTicksToX(ticks, (int)iplID.Plane, (int)tpc, (int)cstat);
      vXsigma[ivx] = std::abs(vXp - vX[ivx]);
*/
    } // ivx
    
    // temp vector of all 2D vertex matches
    std::vector<Vtx3Store> v3temp;
    
    TrajPoint tp;
    float maxScore = 0;
    // i, j, k indicates 3 different wire planes
    // compare vertices in each view
    for(unsigned short ipl = 0; ipl < 2; ++ipl) {
      for(unsigned short ii = 0; ii < vIndex[ipl].size(); ++ii) {
        unsigned short ivx = vIndex[ipl][ii];
        if(vX[ivx] < 0) continue;
        unsigned int iWire = std::nearbyint(tjs.vtx[ivx].Pos[0]);
        for(unsigned short jpl = ipl + 1; jpl < 3; ++jpl) {
          for(unsigned short jj = 0; jj < vIndex[jpl].size(); ++jj) {
            unsigned short jvx = vIndex[jpl][jj];
            if(vX[jvx] < 0) continue;
            unsigned int jWire = std::nearbyint(tjs.vtx[jvx].Pos[0]);
            float dX = std::abs(vX[ivx] - vX[jvx]);
            if(dX > tjs.Vertex3DCuts[0]) continue;
            
            if(prt) mf::LogVerbatim("TC")<<"F3DV: ipl "<<ipl<<" ivxID "<<tjs.vtx[ivx].ID<<" iX "<<vX[ivx]
              <<" jpl "<<jpl<<" jvxID "<<tjs.vtx[jvx].ID<<" jvX "<<vX[jvx]<<" W:T "<<(int)tjs.vtx[jvx].Pos[0]<<":"<<(int)tjs.vtx[jvx].Pos[1]<<" dX "<<dX;
            
            double y = -1000, z = -1000;
            tjs.geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
            if(y < tjs.YLo || y > tjs.YHi || z < tjs.ZLo || z > tjs.ZHi) continue;
            unsigned short kpl = 3 - ipl - jpl;
            float kX = 0.5 * (vX[ivx] + vX[jvx]);
            float kWire = -1;
            if(tjs.NumPlanes > 2) {
              kWire = tjs.geom->WireCoordinate(y, z, kpl, tpc, cstat) + 0.5;
              if(kWire < 0 || (unsigned int)kWire > tjs.NumWires[kpl]) continue;
              if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, kpl, kWire))) continue;
              tp.Pos[0] = kWire;
              // See if there is a wire signal nearby in kpl
              tp.Pos[1] = tjs.detprop->ConvertXToTicks(kX, kpl, tpc, cstat) * tjs.UnitsPerTick;
              tp.CTP = EncodeCTP(cstat, tpc, kpl);
              bool sigOK = SignalAtTp(tjs, tp);
              if(!sigOK) continue;
              if(prt) mf::LogVerbatim("TC")<<" signal exists at "<<kpl<<":"<<PrintPos(tjs, tp);
            }
            kpl = 3 - ipl - jpl;
            // save this incomplete 3D vertex
            Vtx3Store v3d;
            v3d.Vx2ID[ipl] = ivx + 1;
            v3d.Vx2ID[jpl] = jvx + 1;
            v3d.Vx2ID[kpl] = 0;
            // see if this is already in the list
            bool gotit = false;
            for(unsigned short i3t = 0; i3t < v3temp.size(); ++i3t) {
              if(v3temp[i3t].Vx2ID[0] == v3d.Vx2ID[0] && v3temp[i3t].Vx2ID[1] == v3d.Vx2ID[1] && v3temp[i3t].Vx2ID[2] == v3d.Vx2ID[2]) {
                gotit = true;
                break;
              }
            } // i3t
            if(gotit) continue;
            v3d.X = kX;
            // Use XErr to store dX
            v3d.XErr = dX;
            v3d.Y = y;
            v3d.Z = z;
            v3d.Wire = kWire;
            v3d.Score = dX / tjs.Vertex3DCuts[0];
            v3d.TPCID = tpcid;
            // push the incomplete vertex onto the list
            v3temp.push_back(v3d);
            
            if(prt) mf::LogVerbatim("TC")<<"F3DV: 2 Plane match ivxID "<<tjs.vtx[ivx].ID<<" P:W:T "<<ipl<<":"<<(int)tjs.vtx[ivx].Pos[0]<<":"<<(int)tjs.vtx[ivx].Pos[1]<<" jvxID "<<tjs.vtx[jvx].ID<<" P:W:T "<<jpl<<":"<<(int)tjs.vtx[jvx].Pos[0]<<":"<<(int)tjs.vtx[jvx].Pos[1]<<" dX "<<dX;
            
            if(tjs.NumPlanes == 2) continue;
            
            // look for a 3 plane match
            for(unsigned short kk = 0; kk < vIndex[kpl].size(); ++kk) {
              unsigned short kvx = vIndex[kpl][kk];
              if(vX[kvx] < 0) continue;
//              if(vPtr[kvx] >= 0) continue;
              float dX = std::abs(vX[kvx] - vX[ivx]);
              if(dX > tjs.Vertex3DCuts[0]) continue;
              // Wire difference error
              float dW = wirePitch * std::abs(tjs.vtx[kvx].Pos[0] - kWire);
              if(dW > tjs.Vertex3DCuts[1]) continue;
              // put the Y,Z difference in YErr and ZErr
              double y = -1000, z = -1000;
              tjs.geom->IntersectionPoint(iWire, kWire, ipl, kpl, cstat, tpc, y, z);
              v3d.YErr = y - v3d.Y;
              v3d.ZErr = z - v3d.Z;
              v3d.Vx2ID[kpl] = kvx + 1;
              v3d.Wire = -1;
              // hijack the Score variable to hold the separation^2, weighted by the
              // vertex3DCuts
              dX = (vX[kvx] - v3d.X) / tjs.Vertex3DCuts[0];
              float dY = v3d.YErr / tjs.Vertex3DCuts[1];
              float dZ = v3d.ZErr / tjs.Vertex3DCuts[1];
              v3d.Score = dX * dX + dY * dY + dZ * dZ;
              if(v3d.Score > maxScore) maxScore = v3d.Score;
              v3temp.push_back(v3d);
/*
              // X difference error
              dX = (vX[kvx] - kX) / dXSigma;
              float kChi = 0.5 * sqrt(dW * dW + dX * dX);
              if(kChi < tjs.Vertex3DCuts[1]) {
                // push all complete vertices onto the list
                v3d.X = (vX[kvx] + 2 * kX) / 3;
                v3d.XErr = kChi;
                v3d.Vx2ID[kpl] = kvx + 1;
                // see if this is already in the list
                gotit = false;
                for(unsigned short i3t = 0; i3t < v3temp.size(); ++i3t) {
                  if(v3temp[i3t].Vx2ID[0] == v3d.Vx2ID[0] && v3temp[i3t].Vx2ID[1] == v3d.Vx2ID[1] && v3temp[i3t].Vx2ID[2] == v3d.Vx2ID[2]) {
                    gotit = true;
                    break;
                  }
                } // i3t
                if(gotit) continue;
                v3d.Wire = -1;
                v3temp.push_back(v3d);
                if(prt) mf::LogVerbatim("TC")<<" kvx "<<kvx<<" kpl "<<kpl
                  <<" wire "<<(int)tjs.vtx[kvx].Pos[0]<<" kTime "<<(int)tjs.vtx[kvx].Pos[1]<<" kChi "<<kChi<<" dW "<<tjs.vtx[kvx].Pos[0] - kWire;
              } // kChi < best
*/
            } // kk
          } // jj
        } // jpl
      } // ii
    } // ipl
    
    if(v3temp.empty()) return;
    
    // We will sort this list by increasing score. First add the maxScore for 2-plane matches so that
    // they are considered after the 3-plane matches
    maxScore += 1;
    for(auto& v3 : v3temp) if(v3.Wire >= 0) v3.Score += maxScore;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"v3temp list";
      for(auto& v3 : v3temp) {
        mf::LogVerbatim("TC")<<v3.Vx2ID[0]<<" "<<v3.Vx2ID[1]<<" "<<v3.Vx2ID[2]<<" wire "<<v3.Wire<<" "<<v3.Score;
      } // v3
    }
    SortEntry sEntry;
    std::vector<SortEntry> sortVec(v3temp.size());
    for(unsigned short ivx = 0; ivx < v3temp.size(); ++ivx) {
      sEntry.index = ivx;
      sEntry.val = v3temp[ivx].Score;
      sortVec[ivx] = sEntry;
    } // ivx
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valIncreasing);
    // create a new vector of selected 3D vertices
    std::vector<Vtx3Store> v3sel;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short ivx = sortVec[ii].index;
      // ensure that all 2D vertices are unique
      bool skipit = false;
      for(auto& v3 : v3sel) {
        for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
          if(v3temp[ivx].Vx2ID[ipl] == 0) continue;
          if(v3temp[ivx].Vx2ID[ipl] == v3.Vx2ID[ipl]) {
            skipit = true;
            break;
          }
        } // ipl
        if(skipit) break;
      } // v3
      if(skipit) continue;
      v3sel.push_back(v3temp[ivx]);
    } // ii
    v3temp.clear();
    
    if(prt) {
      mf::LogVerbatim("TC")<<"v3sel list";
      for(auto& v3d : v3sel) {
        mf::LogVerbatim("TC")<<v3d.Vx2ID[0]<<" "<<v3d.Vx2ID[1]<<" "<<v3d.Vx2ID[2]<<" wire "<<v3d.Wire<<" "<<v3d.Score;
      } // v3d
    }
    
    // Count the number of incomplete vertices and store
    unsigned short ninc = 0;
    for(auto& vx3 : v3sel) {
      if(tjs.NumPlanes == 2) {
        vx3.Vx2ID[2] = 666;
      } else {
        if(vx3.Wire >= 0) ++ninc;
      }
      if(prt) mf::LogVerbatim("TC")<<"3D vtx "<<tjs.vtx3.size()<<" Vtx2ID "<<vx3.Vx2ID[0]<<" "<<vx3.Vx2ID[1]<<" "<<vx3.Vx2ID[2]
        <<" wire "<<vx3.Wire;
      vx3.ID = tjs.vtx3.size() + 1;
      tjs.vtx3.push_back(vx3);
      // make the 2D -> 3D associations
      for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
        if(vx3.Vx2ID[ipl] == 0) continue;
        VtxStore& vx2 = tjs.vtx[vx3.Vx2ID[ipl]-1];
        vx2.Vx3ID = vx3.ID;
      } // ipl
    } // ivx
    
    // Try to complete incomplete vertices
    if(ninc > 0) {
      CompleteIncomplete3DVerticesInGaps(tjs, tpcid);
      CompleteIncomplete3DVertices(tjs, tpcid);
    }
    
    // Score and flag Tjs that are attached to high-score vertices
    // First remove Tj vertex flags
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      tj.AlgMod[kTjHiVx3Score] = false;
    } // tj
    // Score the 2D vertices
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(vx2.CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      SetVx2Score(tjs, vx2, prt);
    } // vx2
    // and the 3D vertices
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      SetVx3Score(tjs, vx3, prt);
    } // vx3

  } // Find3DVertices

  ////////////////////////////////////////////////
  void Match3DVtxTjs(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // Matches Tjs that are attached to 2D vertices that are matched in 3D. This function does not attempt
    // to determine the appropriate ends of matched Tjs when there is a 3D vertex at both ends. This is done
    // in Find3DEndPoints
    
    if(tjs.vtx3.empty()) return;
    if(tjs.matchVec.empty()) return;

    // sort the vertices by decreasing score
    std::vector<SortEntry> sortVec;
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      // put the TjIDs into a vector for matching
      float score = 0;
      auto v3TjIDs = GetVtxTjIDs(tjs, vx3, score);
      if(v3TjIDs.empty()) continue;
      if(score < tjs.Vertex2DCuts[7]) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"M3DVTj vx3 "<<vx3.ID<<" score "<<score<<" TjIDs";
        for(auto& tjID : v3TjIDs) myprt<<" "<<tjID;
      }
      SortEntry se;
      se.index = vx3.ID - 1;
      se.val = score;
      sortVec.push_back(se);
    } // vx3
    if(sortVec.empty()) return;
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasing);
    
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      auto& vx3 = tjs.vtx3[sortVec[ii].index];
      float score = 0;
      auto v3TjIDs = GetVtxTjIDs(tjs, vx3, score);
      // flag Tjs that have a large separation from the 2D vertex
      unsigned short nfar = 0;
      for(unsigned short itj = 0; itj < v3TjIDs.size(); ++itj) {
        auto& tj = tjs.allTraj[v3TjIDs[itj] - 1];
        float sep = 0;
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == 0) continue;
          auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
          if(vx2.Vx3ID == vx3.ID) sep = PosSep(vx2.Pos, tj.Pts[tj.EndPt[end]].Pos);
        } // end
        if(sep > 5) {
          v3TjIDs[itj] *= -1;
          ++nfar;
        }
      } // itj
      std::sort(v3TjIDs.begin(), v3TjIDs.end());
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<ii<<" vx3 "<<vx3.ID<<" "<<score<<" TjIDs tagged";
        for(auto& tjID : v3TjIDs) myprt<<" "<<tjID;
        myprt<<" nfar "<<nfar;
      }
      // look for these in matchVec with two iterations. Find matching Tjs close to the
      // vertex on the first pass and those farther away on the second
      for(unsigned short nit = 0; nit < 2; ++nit) {
        for(unsigned int ims = 0; ims < tjs.matchVec.size(); ++ims) {
          // count the number of shower Tjs
          unsigned short nstj = 0;
          auto& ms = tjs.matchVec[ims];
          bool skipit = false;
          for(unsigned short ipl = 0; ipl < ms.TjIDs.size(); ++ipl) {
            unsigned short itj = ms.TjIDs[ipl] - 1;
            if(tjs.allTraj[itj].AlgMod[kMat3D]) skipit = true;
            if(tjs.allTraj[itj].AlgMod[kShowerTj]) ++nstj;
          }
          if(skipit) continue;
          // Don't consider shower Tjs
          if(nstj != 0) continue;
          // make a copy of the TjIDs so they can be sorted in increasing order so
          // that std::set_intersection works properly
          auto mstjids = ms.TjIDs;
          std::sort(mstjids.begin(), mstjids.end());
          std::vector<int> shared;
          std::set_intersection(v3TjIDs.begin(), v3TjIDs.end(), 
                                mstjids.begin(), mstjids.end(), std::back_inserter(shared));
          if(shared.size() != mstjids.size()) continue;
          // perfect match. Ensure that the points near the vertex are consistent
          PFPStruct pfp = CreatePFPStruct(tjs, tpcid);
          pfp.TjIDs = shared;
//          TagBragg(tjs, pfp, prt);
          // declare a start or end vertex and set the end points
          if(pfp.Vx3ID[0] == 0) {
            pfp.Vx3ID[0] = vx3.ID;
            if(!SetPFPEndPoints(tjs, pfp, 0, prt)) continue;
          } else {
            pfp.Vx3ID[1] = vx3.ID;
            if(!SetPFPEndPoints(tjs, pfp, 1, prt)) continue;
          }
          tjs.pfps.push_back(pfp);
          std::vector<int> leftover(v3TjIDs.size());
          auto it = std::set_difference(v3TjIDs.begin(), v3TjIDs.end(), shared.begin(), shared.end(), leftover.begin());
          leftover.resize(it - leftover.begin());
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<"nit "<<nit<<" perfect match with ims "<<ims<<" TjIDs";
            for(auto& tjID : mstjids) myprt<<" "<<tjID;
            myprt<<" leftover";
            for(auto& tjID : leftover) myprt<<" "<<tjID;
          }
          // flag these Tjs as matched
          for(auto id : mstjids) tjs.allTraj[id - 1].AlgMod[kMat3D] = true;
          if(leftover.empty()) break;
          // keep looking using the leftovers
          v3TjIDs = leftover;
        } // ims
        if(nfar == 0) break;
        for(auto& id : v3TjIDs) id = abs(id);
        std::sort(v3TjIDs.begin(), v3TjIDs.end());
      } // nit
      if(v3TjIDs.empty()) continue;
    } // ii (vx3 sorted)

  } // Match3DVtxTjs

  //////////////////////////////////////////
  unsigned short TPNearVertex(TjStuff& tjs, const TrajPoint& tp)
  {
    // Returns the index of a vertex if tp is nearby
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      if(tjs.vtx[ivx].ID == 0) continue;
      if(tjs.vtx[ivx].CTP != tp.CTP) continue;
      if(std::abs(tjs.vtx[ivx].Pos[0] - tp.Pos[0]) > 1.2) continue;
      if(std::abs(tjs.vtx[ivx].Pos[1] - tp.Pos[1]) > 1.2) continue;
      return ivx;
    } // ivx
    return USHRT_MAX;
  } // TPNearVertex
  
  //////////////////////////////////////////
  bool AttachPFPToVertex(TjStuff& tjs, PFPStruct& pfp, unsigned short end, unsigned short vx3ID, bool prt)
  {
    if(vx3ID > tjs.vtx3.size()) {
      mf::LogVerbatim("TC")<<"AttachPFPToVertex: vx3 "<<vx3ID<<" doesn't exist in tjs.vtx3";
      return false;
    }
    if(pfp.ID > tjs.pfps.size()) {
      mf::LogVerbatim("TC")<<"AttachPFPToVertex: pfp "<<pfp.ID<<" doesn't exist in tjs.pfps";
      return false;
    }
    if(end > 1) return false;
    
    auto& vx3 = tjs.vtx3[vx3ID - 1];
    
    pfp.Vx3ID[end] = vx3.ID;
    
    // We are done if this a PFP-only vertex
    if(vx3.Wire == -2) return true;
    
    if(prt) {
      std::cout<<"APTV: pfp.ID "<<pfp.ID<<" end "<<end<<" vx3.ID "<<vx3.ID<<" vx3.Vx2ID";
      for(unsigned short plane = 0; plane < 3; ++plane) std::cout<<" "<<vx3.Vx2ID[plane];
      std::cout<<"\n";
    }
    
    // Update the 2D and 3D vertex and tj associations
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      // TODO: Check to see if the Tjs have been ordered correctly? 
      if(tj.VtxID[end] == 0) {
        // tj is available to be attached to a 2D vertex. See if the 3D vertex is matched to 
        // an existing 2D vertex in this plane
        if(vx3.Vx2ID[plane] == 0) {
          // not matched. Look for one
          std::array<float, 2> pos;
          PosInPlane(tjs, vx3, plane, pos);
          if(prt) std::cout<<" tj "<<tj.ID<<" has no 2D vertex. Look for one vertex near "<<tj.CTP<<":"<<PrintPos(tjs, pos)<<" Events processed "<<tjs.EventsProcessed<<"\n";
        } else {
          // Existing 2D vertex matched to the 3D vertex
          if(prt) std::cout<<" tj "<<tj.ID<<" has no 2D vertex in CTP "<<tj.CTP<<" but vx3 is matched to 2D vertex"<<vx3.Vx2ID[plane]<<". Attach it? Events processed "<<tjs.EventsProcessed<<"\n";
        }
      }
    } // tjid
    
    return true;
  } // AttachPFPToVertex

  //////////////////////////////////////////
  bool AttachAnyTrajToVertex(TjStuff& tjs, unsigned short ivx, bool prt)
  {
    
    if(ivx > tjs.vtx.size() - 1) return false;
    if(tjs.vtx[ivx].ID == 0) return false;
    if(tjs.Vertex2DCuts[0] < 0) return false;
    
    VtxStore& vx = tjs.vtx[ivx];
    
    unsigned short nadd = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx.CTP) continue;
      if(tj.VtxID[0] == vx.ID || tj.VtxID[1] == vx.ID) continue;
      if(AttachTrajToVertex(tjs, tj, vx, prt)) ++nadd;
    } // tj
    if(prt) mf::LogVerbatim("TC")<<" AttachAnyTrajToVertex: nadd "<<nadd;
    if(nadd == 0) return false;
    return true;
    
  } // AttachAnyTrajToVertex

  //////////////////////////////////////////
  bool AttachTrajToVertex(TjStuff& tjs, Trajectory& tj, VtxStore& vx, bool prt)
  {
    // Note that this function does not require a signal between the end of the Tj and the vertex
    
    // tjs.Vertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction
    
    if(tj.AlgMod[kKilled]) return false;
    if(tj.CTP != vx.CTP) return false;
    // already attached?
    if(tj.VtxID[0] == vx.ID || tj.VtxID[1] == vx.ID) return false;
    
    unsigned short maxShortTjLen = tjs.Vertex2DCuts[0];
    // square the separation cut to simplify testing in the loop
    float maxSepCutShort2 = tjs.Vertex2DCuts[1] * tjs.Vertex2DCuts[1];
    float maxSepCutLong2 = tjs.Vertex2DCuts[2] * tjs.Vertex2DCuts[2];
    
    // assume that end 0 is closest to the vertex
    unsigned short end = 0;
    float vtxTjSep2 = PosSep2(vx.Pos, tj.Pts[tj.EndPt[0]].Pos);
    float sep1 = PosSep2(vx.Pos, tj.Pts[tj.EndPt[1]].Pos);
    if(sep1 < vtxTjSep2) {
      // End 1 is closer
      end = 1;
      vtxTjSep2 = sep1;
    }
    // is there a vertex already assigned to this end?
    if(tj.VtxID[end] > 0) return false;
    
    // is the trajectory short?
    bool tjShort = (tj.EndPt[1] - tj.EndPt[0] < maxShortTjLen);
    // ignore bad separation between the closest tj end and the vertex
    if(tjShort) {
      if(vtxTjSep2 > maxSepCutShort2) return false;
    } else {
      if(vtxTjSep2 > maxSepCutLong2) return false;
    }
    
    // Calculate the pull on the vertex
    TrajPoint& tp = tj.Pts[tj.EndPt[end]];
    float tpVxPull = TrajPointVertexPull(tjs, tp, vx);
    
    // See if the vertex position is close to an end
    unsigned short closePt;
    float closestApproach;
    TrajClosestApproach(tj, vx.Pos[0], vx.Pos[1], closePt, closestApproach);
    // count the number of points between the end of the trajectory and the vertex.
    // tj     -------------   tj ------------
    // vx  *   >> dpt = 0     vx   *  >> dpt = 2
    short dpt;
    if(end == 0) {
      dpt = closePt - tj.EndPt[end];
    } else {
      dpt = tj.EndPt[end] - closePt;
    }
    
    float length = TrajLength(tj);
    // don't attach it if the tj length is shorter than the separation distance
    if(length > 2 && length < closestApproach) return false;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"ATTV: vx.ID "<<vx.ID;
      myprt<<" NTraj "<<vx.NTraj;
      myprt<<" oldTJs";
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.AlgMod[kKilled]) continue;
        if(tj.CTP != vx.CTP) continue;
        if(tj.VtxID[0] == vx.ID) myprt<<" "<<tj.ID<<"_0";
        if(tj.VtxID[1] == vx.ID) myprt<<" "<<tj.ID<<"_1";
      }
      myprt<<" +tjID "<<tj.ID<<"_"<<end<<" vtxTjSep "<<sqrt(vtxTjSep2)<<" tpVxPull "<<tpVxPull<<" tjs.Vertex2DCuts[3] "<<tjs.Vertex2DCuts[3];
    }
    if(tpVxPull > tjs.Vertex2DCuts[3]) return false;
    if(dpt > 2) return false;
    
    // remove the fixed position flag if there are more than 2 tjs
    bool fixedBit = vx.Stat[kFixed];
    if(fixedBit && vx.NTraj < 2) vx.Stat[kFixed] = false;

    // Passed all the cuts. Attach it to the vertex and try a fit
    tj.VtxID[end] = vx.ID;
    if(FitVertex(tjs, vx, prt)) {
      SetVx2Score(tjs, vx, prt);
      if(prt) mf::LogVerbatim("TC")<<" success";
      return true;
    }
    
    // fit failed so remove the tj -> vx assignment
    tj.VtxID[end] = 0;
    // restore the fixed flag
    vx.Stat[kFixed] = fixedBit;
    // and refit
    if(prt) mf::LogVerbatim("TC")<<" failed. Re-fit w/o this tj ";
    FitVertex(tjs, vx, prt);
    return false;
    
  } // AttachTrajToVertex
  
  /////////////////////////////////////////
  float TrajPointVertexPull(TjStuff& tjs, const TrajPoint& tp, const VtxStore& vx)
  {
    // Calculates the position pull between a trajectory point and a vertex
    
    // impact parameter between them
    double ip = PointTrajDOCA(tjs, vx.Pos[0], vx.Pos[1], tp);
    // separation^2
    double sep2 = PosSep2(vx.Pos, tp.Pos);
    
    // Find the projection of the vertex error ellipse in a coordinate system
    // along the TP direction
    double vxErrW = vx.PosErr[0] * tp.Dir[1];
    double vxErrT = vx.PosErr[1] * tp.Dir[0];
    double vxErr2 = vxErrW * vxErrW + vxErrT * vxErrT;
    
    // close together so ignore the TP projection error and return
    // the pull using only the vertex error
    if(sep2 < 1) return (float)(ip/sqrt(vxErr2));
    
    double dang = ip / sqrt(sep2);
    
    // calculate the angle error.
    // Start with the vertex error^2
    double angErr = vxErr2 / sep2;
    // Add the TP angle error^2
    angErr += tp.AngErr * tp.AngErr;
    if(angErr == 0) return 999;
    angErr = sqrt(angErr);
    return (float)(dang / angErr);
    
  } // TrajPointVertexPull
  
  /////////////////////////////////////////
  float VertexVertexPull(TjStuff& tjs, const Vtx3Store& vx1, const Vtx3Store& vx2)
  {
    // Calculates the position pull between two vertices
    double dx = vx1.X - vx2.X;
    double dy = vx1.Y - vx2.Y;
    double dz = vx1.Z - vx2.Z;
    double dxErr2 = (vx1.XErr * vx1.XErr + vx2.XErr * vx2.XErr) / 2;
    double dyErr2 = (vx1.YErr * vx1.YErr + vx2.YErr * vx2.YErr) / 2;
    double dzErr2 = (vx1.ZErr * vx1.ZErr + vx2.ZErr * vx2.ZErr) / 2;
    dx = dx * dx / dxErr2;
    dy = dy * dy / dyErr2;
    dz = dz * dz / dzErr2;
    return (float)(sqrt(dx + dy + dz)/3);
  }
  
  /////////////////////////////////////////
  float VertexVertexPull(TjStuff& tjs, const VtxStore& vx1, const VtxStore& vx2)
  {
    // Calculates the position pull between two vertices
    double dw = vx1.Pos[0] - vx2.Pos[0];
    double dt = vx1.Pos[1] - vx2.Pos[1];
    double dwErr2 = (vx1.PosErr[0] * vx1.PosErr[0] + vx2.PosErr[0] * vx2.PosErr[0]) / 2;
    double dtErr2 = (vx1.PosErr[1] * vx1.PosErr[1] + vx2.PosErr[1] * vx2.PosErr[1]) / 2;
    dw = dw * dw / dwErr2;
    dt = dt * dt / dtErr2;
    return (float)sqrt(dw + dt);
  }
  
  ////////////////////////////////////////////////
  bool StoreVertex(TjStuff& tjs, VtxStore& vx)
  {
    // jacket around the push to ensure that the Tj and vtx CTP is consistent.
    // The calling function should score the vertex after the trajectories are attached
    
    if(vx.ID != tjs.vtx.size() + 1) {
      mf::LogVerbatim("TC")<<"StoreVertex: Invalid ID "<<vx.ID<<" It should be "<<tjs.vtx.size() + 1;
      return false;
    }
    
    unsigned short nvxtj = 0;
    unsigned short nok = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(vx.ID == tj.VtxID[0] || vx.ID == tj.VtxID[1]) ++nvxtj;
      if(vx.CTP != tj.CTP) continue;
      if(vx.ID == tj.VtxID[0] || vx.ID == tj.VtxID[1]) ++nok;
    } // tj
    
    if(nok != nvxtj) {
      mf::LogVerbatim("TC")<<"StoreVertex: vertex "<<vx.ID<<" has inconsistent CTP code "<<vx.CTP<<" with one or more Tjs\n";
      for(auto& tj : tjs.allTraj) {
        if(tj.AlgMod[kKilled]) continue;
        mf::LogVerbatim("TC")<<"Tj ID "<<tj.ID<<" CTP "<<tj.CTP<<"\n";
        if(tj.VtxID[0] == vx.ID) tj.VtxID[0] = 0;
        if(tj.VtxID[1] == vx.ID) tj.VtxID[1] = 0;
      }
      return false;
    }
    vx.NTraj = nok;
    tjs.vtx.push_back(vx);
    return true;
    
  } // StoreVertex
  
  /////////////////////////////////////////
  bool FitVertex(TjStuff& tjs, VtxStore& vx, bool prt)
  {
    // A poor-mans fitting scheme. If the fitted vertex position error is within the supplied
    // value, the position and errors are updated and we return true, otherwise the vertex is
    // left unchanged and we return false
    
    // tjs.Vertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction
    
    if(vx.Stat[kFixed]) {
      if(prt) mf::LogVerbatim("TC")<<" vertex position fixed. No fit allowed";
      return true;
    }
    
    // Create a vector of trajectory points that will be used to fit the vertex position
    std::vector<TrajPoint> vxTp;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx.CTP) continue;
      if(tj.VtxID[0] == vx.ID) vxTp.push_back(tj.Pts[tj.EndPt[0]]);
      if(tj.VtxID[1] == vx.ID) vxTp.push_back(tj.Pts[tj.EndPt[1]]);
    } // tj
    
    bool success = FitVertex(tjs, vx, vxTp, prt);
    
    if(!success) return false;
    return true;
    
  } // FitVertex
  
  /////////////////////////////////////////
  bool FitVertex(TjStuff& tjs, VtxStore& vx, std::vector<TrajPoint> vxTp, bool prt)
  {
    // Variant of FitVertex that fits the passed trajectory points to a vertex position but doesn't
    // require using information in TJStuff

    // tjs.Vertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction
    // 7 min Score
    // 8 Min charge fraction near a merge point (not a vertex)
    // 9 max MCSmom asymmetry for a merge

    
    vx.NTraj = vxTp.size();
    
    if(vxTp.size() < 2) return false;
    
    if(prt) {
      PrintHeader("FV");
      for(auto& tp : vxTp) PrintTrajPoint("FV", tjs, 0, 1, 1, tp);
    }
    
    // Find trajectory intersections pair-wise tweaking the angle and position(?) within
    // +/- 1 sigma
    double sum0 = 0, sum02 = 0;
    double sum1 = 0, sum12 = 0;
    double sumw = 0;
    double wgt;
    // a temporary TP for tweaking the angle
    TrajPoint tmp;
    // another point to check for a signal at each intersection
    TrajPoint intTp;
    intTp.CTP = vxTp[0].CTP;
    for(unsigned short itj = 0; itj < vxTp.size() - 1; ++itj) {
      for(unsigned short jtj = itj + 1; jtj < vxTp.size(); ++jtj) {
        float p0, p1;
        TrajIntersection(vxTp[itj], vxTp[jtj], p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        wgt = 1;
        // accumulate
        // Requiring SignalAtTp results in 1% to 3% drop in EP 
//        if(SignalAtTp(tjs, intTp)) {
          sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
//        }
        // tweak the itj angle +
        tmp = vxTp[itj];
        tmp.Ang += tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(tmp, vxTp[jtj], p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // adjust the weight for 4 points at +/1 1 sigma = 0.607 / 4
        wgt = 0.152;
        // accumulate
//        if(SignalAtTp(tjs, intTp)) {
          sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
//        }
        // tweak the itj angle -
        tmp = vxTp[itj];
        tmp.Ang -= 2 * tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(tmp, vxTp[jtj], p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // accumulate
//        if(SignalAtTp(tjs, intTp)) {
          sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
//        }
        // Repeat this process with jtj
        // tweak the jtj angle +
        tmp = vxTp[jtj];
        tmp.Ang += tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(vxTp[itj], tmp, p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // accumulate
//        if(SignalAtTp(tjs, intTp)) {
          sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
//        }
        // tweak the itj angle -
        tmp = vxTp[itj];
        tmp.Ang -= 2 * tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(vxTp[itj], tmp, p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // accumulate
//        if(SignalAtTp(tjs, intTp)) {
          sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
//        }
      } // jtj
    } // itj
    
    if(sumw == 0) return false;
    
    double vxP0 = sum0 / sumw;
    double vxP1 = sum1 / sumw;
    double vxP0rms = sqrt((sum02 - sumw * vxP0 * vxP0) / sumw);
    double vxP1rms = sqrt((sum12 - sumw * vxP1 * vxP1) / sumw);
    // don't let the errors get too small
    if(vxP0rms < 0.5) vxP0rms = 0.5;
    if(vxP1rms < 0.5) vxP1rms = 0.5;
    
    if(prt) mf::LogVerbatim("TC")<<"FitVertex Vx_"<<vx.ID<<" CTP "<<vx.CTP<<" NTraj "<<vx.NTraj<<" in "<<std::fixed<<std::setprecision(1)<<vx.Pos[0]<<":"<<vx.Pos[1]/tjs.UnitsPerTick<<" out wire "<<vxP0<<" +/- "<<vxP0rms<<" ticks "<<vxP1/tjs.UnitsPerTick<<"+/-"<<vxP1rms/tjs.UnitsPerTick;
    
    if(vxP0rms > tjs.Vertex2DCuts[4] || vxP1rms > tjs.Vertex2DCuts[4]) {
      if(prt) mf::LogVerbatim("TC")<<" fit failed. tjs.Vertex2DCuts[4] "<<tjs.Vertex2DCuts[4];
      return false;
    }
    
    vx.Pos[0] = vxP0;
    vx.PosErr[0] = vxP0rms;
    vx.Pos[1] = vxP1;
    vx.PosErr[1] = vxP1rms;
    
    // Calculate chisq
    vx.ChiDOF = 0;
    for(unsigned short itj = 0; itj < vxTp.size(); ++itj) {
      vx.ChiDOF += TrajPointVertexPull(tjs, vxTp[itj], vx);
    } // itj
    vx.ChiDOF /= (float)vxTp.size();
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Pull";
      for(unsigned short itj = 0; itj < vxTp.size(); ++itj) {
        float pull = TrajPointVertexPull(tjs, vxTp[itj], vx);
        myprt<<" "<<PrintPos(tjs, vxTp[itj])<<" - "<<std::fixed<<std::setprecision(2)<<pull;
      } // itj
      myprt<<" ChiDOF "<<vx.ChiDOF;
    }
    return true;

  } // FitVertex

  //////////////////////////////////////////
  bool ChkVtxAssociations(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Check the associations

    // check the 2D -> 3D associations
    geo::PlaneID planeID = DecodeCTP(inCTP);
    unsigned short plane = planeID.Plane;
    for(auto& vx2 : tjs.vtx) {
      if(vx2.CTP != inCTP) continue;
      if(vx2.ID == 0) continue;
      if(vx2.Vx3ID == 0) continue;
      if(vx2.Vx3ID > tjs.vtx3.size()) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: Invalid vx2.Vx3ID "<<vx2.Vx3ID<<" in 2D vtx "<<vx2.ID;
        return false;
      }
      auto& vx3 = tjs.vtx3[vx2.Vx3ID-1];
      if(vx3.ID == 0) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: vx2 "<<vx2.ID<<" thinks it is matched to vx3 "<<vx3.ID<<" but vx3 is obsolete";
        return false;
      }
      if(vx3.Vx2ID[plane] != vx2.ID) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: vx2 "<<vx2.ID<<" thinks it is matched to vx3 "<<vx3.ID<<" but vx3 says no!";
        return false;
      }
    } // vx2
    // check the 3D -> 2D associations
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.Vx2ID[plane] == 0) continue;
      if(vx3.Vx2ID[plane] > tjs.vtx.size()) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: Invalid vx3.Vx2ID "<<vx3.Vx2ID[plane]<<" in CTP "<<inCTP;
        return false;
      }
      auto& vx2 = tjs.vtx[vx3.Vx2ID[plane]-1];
      if(vx2.Vx3ID != vx3.ID) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: vx3 "<<vx3.ID<<" thinks it is matched to vx2 "<<vx2.ID<<" but vx2 says no!";
        return false;
      }
    } // vx3
    
    // check the Tj -> 2D associations
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == 0) continue;
        if(tj.VtxID[end] > tjs.vtx.size()) {
          mf::LogVerbatim("TC")<<"ChkVtxAssociations: tj "<<tj.ID<<" thinks it is matched to vx2 "<<tj.VtxID[end]<<" on end "<<end<<" but no vertex exists. Recovering";
          tj.VtxID[end] = 0;
          return false;
        }
        unsigned short ivx = tj.VtxID[end] - 1;
        auto& vx2 = tjs.vtx[ivx];
        if(vx2.ID == 0) {
          mf::LogVerbatim("TC")<<"ChkVtxAssociations: tj "<<tj.ID<<" thinks it is matched to vx2 "<<tj.VtxID[end]<<" on end "<<end<<" but the vertex is killed. Fixing the Tj";
          tj.VtxID[end] = 0;
          return false;
        }
      } // end
    } // tj
    
    return true;
    
  } // ChkVtxAssociations
  
  //////////////////////////////////////////
  void ScoreVertices(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // reset all 3D vertex, 2D vertex and Tj high-score vertex bits in tpcid 
    
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    
    // reset the 2D vertex status bits
    for(auto& vx : tjs.vtx) {
      if(vx.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(vx.CTP);
      if(planeID.Cryostat != cstat) continue;
      if(planeID.TPC != tpc) continue;
      vx.Stat[kHiVx3Score] = false;
    } // vx
    // and the tj bits
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.Cryostat != cstat) continue;
      if(planeID.TPC != tpc) continue;
      tj.AlgMod[kTjHiVx3Score] = false;
    } // tj
    // Score the 2D vertices
    for(auto& vx : tjs.vtx) {
      if(vx.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(vx.CTP);
      if(planeID.Cryostat != cstat) continue;
      if(planeID.TPC != tpc) continue;
      SetVx2Score(tjs, vx, prt);
      if(vx.Score < tjs.Vertex2DCuts[7]) MakeVertexObsolete(tjs, vx, false);
    } // vx
    // Score the 3D vertices
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
       SetVx3Score(tjs, vx3, prt);
    } // vx3
    
  } // ScoreVertices
  
  //////////////////////////////////////////
  void SetHighScoreBits(TjStuff& tjs, Vtx3Store& vx3)
  {
    // Sets the tj and 2D vertex score bits to true 
    
    if(vx3.ID == 0) return;
    
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      if(vx3.Vx2ID[ipl] <= 0) continue;
      VtxStore& vx2 = tjs.vtx[vx3.Vx2ID[ipl] - 1];
      vx2.Stat[kHiVx3Score] = false;
      // transfer this to all attached tjs and vertices attached to those tjs
      std::vector<int> tjlist = GetVtxTjIDs(tjs, vx2);
      std::vector<int> vxlist;
//      std::cout<<"SHSB: vx3 "<<vx3.ID<<" vx2 "<<vx2.ID<<"\n";
      while(true) {
        // tag Tjs and make a list of attached vertices whose high-score
        // bit needs to be set
        vxlist.clear();
//        std::cout<<" tjlist";
        for(auto tjid : tjlist) {
          auto& tj = tjs.allTraj[tjid - 1];
          tj.AlgMod[kTjHiVx3Score] = true;
//          std::cout<<" "<<tj.ID;
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] == 0) continue;
            auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
            if(vx2.Stat[kHiVx3Score]) continue;
            vx2.Stat[kHiVx3Score] = true;
            vxlist.push_back(vx2.ID);
          } // end
        } // tjid
/*
        std::cout<<" vxlist ";
        for(auto vxid : vxlist) std::cout<<" "<<vxid;
        std::cout<<"\n";
*/
        if(vxlist.empty()) break;
        // re-build tjlist using vxlist
        std::vector<int> newtjlist;
        for(auto vxid : vxlist) {
          auto& vx2 = tjs.vtx[vxid - 1];
          auto tmp = GetVtxTjIDs(tjs, vx2);
          for(auto tjid : tmp) {
            if(std::find(tjlist.begin(), tjlist.end(), tjid) == tjlist.end()) newtjlist.push_back(tjid);
          } // tjid
        } // vxid
        if(newtjlist.empty()) break;
        tjlist = newtjlist;
      } // true
    } // ipl

  } // SetHighScoreBits
  
  //////////////////////////////////////////
  void SetVx3Score(TjStuff& tjs, Vtx3Store& vx3, bool prt)
  {
    // Calculate the 3D vertex score and flag Tjs that are attached to high score vertices as defined 
    // by Vertex2DCuts 
    
    if(vx3.ID == 0) return;
    
    vx3.Score = 0;
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      if(vx3.Vx2ID[ipl] <= 0) continue;
      VtxStore& vx2 = tjs.vtx[vx3.Vx2ID[ipl] - 1];
      vx3.Score += vx2.Score;
    } // ipl
    vx3.Score /= (float)tjs.NumPlanes;
    if(vx3.Score > tjs.Vertex2DCuts[7]) SetHighScoreBits(tjs, vx3);
    
  } // SetVx3Score
  
  //////////////////////////////////////////
  void SetVx2Score(TjStuff& tjs, bool prt)
  {
    // A version that sets the score of the last added vertex
    if(tjs.vtx.empty()) return;
    auto& vx2 = tjs.vtx[tjs.vtx.size() - 1];
    SetVx2Score(tjs, vx2, prt);
  } // SetVx2Score
  
  //////////////////////////////////////////
  void SetVx2Score(TjStuff& tjs, VtxStore& vx2, bool prt)
  {
    // Calculate the 2D vertex score
    
    // Don't score vertices from CheckTrajBeginChg. Set to the minimum
    if(vx2.Topo == 8) {
      vx2.Score = tjs.Vertex2DCuts[7] + 0.1;
      auto vtxTjID = GetVtxTjIDs(tjs, vx2);
      vx2.TjChgFrac = ChgFracNearPos(tjs, vx2.Pos, vtxTjID);
      return;
    }
    
    // Cuts on Tjs attached to vertices
    constexpr float maxChgRMS = 0.25;
    constexpr float momBin = 50;

    vx2.Score = -1000;
    vx2.TjChgFrac = 0;
    if(vx2.ID == 0) return;    
    if(tjs.VertexScoreWeights.size() < 4) return;
    
    auto vtxTjID = GetVtxTjIDs(tjs, vx2);
    if(vtxTjID.empty()) return;

    // Vertex position error
    float vpeScore = -tjs.VertexScoreWeights[0] * (vx2.PosErr[0] + vx2.PosErr[1]);
    
    unsigned short m3Dcnt = 0;
    if(vx2.Vx3ID > 0) {
      m3Dcnt = 1;
      // Add another if the 3D vertex is complete
      unsigned short ivx3 = vx2.Vx3ID - 1;
      if(tjs.vtx3[ivx3].Wire < 0) m3Dcnt = 2;
    }
    float m3DScore = tjs.VertexScoreWeights[1] * m3Dcnt;
    
    vx2.TjChgFrac = ChgFracNearPos(tjs, vx2.Pos, vtxTjID);
    float cfScore = tjs.VertexScoreWeights[2] * vx2.TjChgFrac;
    
    // Define a weight for each Tj
    std::vector<float> tjwts(vtxTjID.size());
    for(unsigned short it1 = 0; it1 < vtxTjID.size(); ++it1) {
      unsigned short itj1 = vtxTjID[it1] - 1;
      Trajectory& tj1 = tjs.allTraj[itj1];
      float wght1 = (float)tj1.MCSMom / momBin;
      if(wght1 > 10) wght1 = 10;
      // weight by tagged muon
      if(tj1.PDGCode == 13) wght1 *= 2;
      // weight by charge rms
      if(tj1.ChgRMS < maxChgRMS) ++wght1;
      // Shower Tj
      if(tj1.AlgMod[kShowerTj]) ++wght1;
      tjwts[it1] = wght1;
    } // tjid
    
    float tjScore = 0;
    float sum = 0;
    float cnt = 0;
    for(unsigned short it1 = 0; it1 < vtxTjID.size() - 1; ++it1) {
      unsigned short itj1 = vtxTjID[it1] - 1;
      Trajectory& tj1 = tjs.allTraj[itj1];
      float wght1 = tjwts[it1];
      // the end that has a vertex
      unsigned short end1 = 0;
      if(tj1.VtxID[1] == vx2.ID) end1 = 1;
      unsigned short endPt1 = tj1.EndPt[end1];
      // bump up the weight if there is a Bragg peak at the other end
      unsigned short oend1 = 1 - end1;
      if(tj1.StopFlag[oend1][kBragg]) ++wght1;
      float ang1 = tj1.Pts[endPt1].Ang;
      float ang1Err2 = tj1.Pts[endPt1].AngErr * tj1.Pts[endPt1].AngErr;
      for(unsigned short it2 = it1 + 1; it2 < vtxTjID.size(); ++it2) {
        unsigned short itj2 = vtxTjID[it2]  - 1;
        Trajectory& tj2 = tjs.allTraj[itj2];
        float wght2 = tjwts[it2];
        unsigned end2 = 0;
        if(tj2.VtxID[1] == vx2.ID) end2 = 1;
        // bump up the weight if there is a Bragg peak at the other end
        unsigned short oend2 = 1 - end2;
        if(tj2.StopFlag[oend2][kBragg]) ++wght2;
        unsigned short endPt2 = tj2.EndPt[end2];
        float ang2 = tj2.Pts[endPt2].Ang;
        float ang2Err2 = tj2.Pts[endPt2].AngErr * tj2.Pts[endPt2].AngErr;
        float dang = DeltaAngle(ang1, ang2);
        float dangErr = 0.5 * sqrt(ang1Err2 + ang2Err2);
        if((dang / dangErr) > 3 && wght1 > 0 && wght2 > 0) {
          sum += wght1 + wght2;
          ++cnt;
        }
      } // it2
    } // it1
    if(cnt > 0) {
      sum /= cnt;
      tjScore = tjs.VertexScoreWeights[3] * sum;
    }
//    std::cout<<" vpe "<<vpeScore<<" m3d "<<m3DScore<<" cf "<<cfScore<<" tj "<<tjScore<<"\n";
    vx2.Score = vpeScore + m3DScore + cfScore + tjScore;
    if(prt) {
      // last call after vertices have been matched to the truth. Use to optimize VertexScoreWeights using
      // an ntuple
      mf::LogVerbatim myprt("TC");
      myprt<<"VSW "<<vx2.ID;
      myprt<<" m3Dcnt"<<m3Dcnt;
      myprt<<" "<<std::fixed<<std::setprecision(2)<<(vx2.PosErr[0] + vx2.PosErr[1]);
      myprt<<" "<<std::fixed<<std::setprecision(3)<<vx2.TjChgFrac;
      myprt<<" "<<std::fixed<<std::setprecision(1)<<sum;
      myprt<<" "<<(int)cnt;
    }
  } // SetVx2Score
  
  //////////////////////////////////////////
  unsigned short Vx3Topo(TjStuff& tjs, Vtx3Store& vx3)
  {
    // Returns the most common value of Topo for the 2D vertices that are matched
    // to this 3D vertex. This **might** be a useful measure to identify neutrino interaction
    // vertices
    
    if(vx3.ID == 0) return USHRT_MAX;
    // Consider Topo values between 0 and 9
    std::array<short, 10> cnts = {0};
    for(auto vx2id : vx3.Vx2ID) {
      if(vx2id == 0) continue;
      auto& vx2 = tjs.vtx[vx2id - 1];
      if(vx2.Topo < 0 || vx2.Topo > 9) continue;
      ++cnts[vx2.Topo];
    } // vx2id
    short most = 0;
    unsigned short theMost = USHRT_MAX;
    for(unsigned short itp = 0; itp < 10; ++itp) {
      if(cnts[itp] > most) {
        most = cnts[itp];
        theMost = itp;
      }
    } // itp
    return theMost;
  } // Vx3Topo

  //////////////////////////////////////////
  void CompleteIncomplete3DVerticesInGaps(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    
    if(!tjs.UseAlg[kComp3DVxIG]) return;

    bool prt = (debug.Plane >= 0 && debug.Tick == 44444);
    if(prt) mf::LogVerbatim("TC")<<"Inside CI3DVIG:";

    for(unsigned short iv3 = 0; iv3 < tjs.vtx3.size(); ++iv3) {
      Vtx3Store& vx3 = tjs.vtx3[iv3];
      // ignore obsolete vertices
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      // check for a completed 3D vertex
      if(vx3.Wire < 0) continue;
      unsigned short mPlane = USHRT_MAX;
      for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
        if(vx3.Vx2ID[ipl] > 0) continue;
        mPlane = ipl;
        break;
      } // ipl
      if(mPlane == USHRT_MAX) continue;
//      CTP_t mCTP = EncodeCTP(vx3.CStat, vx3.TPC, mPlane);
      CTP_t mCTP = EncodeCTP(vx3.TPCID.Cryostat, vx3.TPCID.TPC, mPlane);
      // require that the missing vertex be in a large block of dead wires
      float dwc = DeadWireCount(tjs, vx3.Wire - 4, vx3.Wire + 4, mCTP);
      if(dwc < 5) continue;
      // X position of the purported missing vertex
      // A TP for the missing 2D vertex
      VtxStore aVtx;
      aVtx.ID = tjs.vtx.size() + 1;
      aVtx.Pos[0] = vx3.Wire;
      aVtx.Pos[1] = tjs.detprop->ConvertXToTicks(vx3.X, mPlane, vx3.TPCID.TPC, vx3.TPCID.Cryostat) * tjs.UnitsPerTick;
      aVtx.CTP = mCTP;
      aVtx.Topo = 4;
      aVtx.NTraj = 0;
      // Give it a bogus pass to indicate it wasn't created while stepping
      aVtx.Pass = 9;
      if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: Incomplete vertex "<<iv3<<" in plane "<<mPlane<<" wire "<<vx3.Wire<<" Made 2D vertex ";
      std::vector<int> tjIDs;
      std::vector<unsigned short> tjEnds;
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        if(tjs.allTraj[itj].CTP != mCTP) continue;
        if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
        for(unsigned short end = 0; end < 2; ++end) {
          unsigned short ept = tjs.allTraj[itj].EndPt[end];
          TrajPoint& tp = tjs.allTraj[itj].Pts[ept];
          unsigned short oept = tjs.allTraj[itj].EndPt[1 - end];
          TrajPoint& otp = tjs.allTraj[itj].Pts[oept];
          // ensure that this is the end closest to the vertex
          if(std::abs(tp.Pos[0] - aVtx.Pos[0]) > std::abs(otp.Pos[0] - aVtx.Pos[0])) continue;
          float doca = PointTrajDOCA(tjs, aVtx.Pos[0], aVtx.Pos[1], tp);
          if(doca > 2) continue;
          float dwc = DeadWireCount(tjs, aVtx.Pos[0], tp.Pos[0], tp.CTP);
          float ptSep;
          if(aVtx.Pos[0] > tp.Pos[0]) {
            ptSep = aVtx.Pos[0] - tp.Pos[0] - dwc;
          } else {
            ptSep = tp.Pos[0] - aVtx.Pos[0] - dwc;
          }
          if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: tj ID "<<tjs.allTraj[itj].ID<<" doca "<<doca<<" ptSep "<<ptSep;
          if(ptSep < -2 || ptSep > 2) continue;
          // don't clobber an existing association
          if(tjs.allTraj[itj].VtxID[end] > 0) continue;
          tjIDs.push_back(tjs.allTraj[itj].ID);
          tjEnds.push_back(end);
        } // end
      } // itj
      if(!tjIDs.empty()) {
        // Determine how messy this region is
        aVtx.TjChgFrac = ChgFracNearPos(tjs, aVtx.Pos, tjIDs);
        if(aVtx.TjChgFrac < 0.7) continue;
        aVtx.Vx3ID = vx3.ID;
        // Save the 2D vertex
        if(!StoreVertex(tjs, aVtx)) continue;
        for(unsigned short ii = 0; ii < tjIDs.size(); ++ii) {
          unsigned short itj = tjIDs[ii] - 1;
          tjs.allTraj[itj].VtxID[tjEnds[ii]] = aVtx.ID;
          tjs.allTraj[itj].AlgMod[kComp3DVxIG] = true;
        } // ii
        SetVx2Score(tjs, prt);
        vx3.Vx2ID[mPlane] = aVtx.ID;
        vx3.Wire = -1;
        if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: new vtx Vx_"<<aVtx.ID<<" points to 3D tjs.vtx ";
      }
    } // vx3
    
  } // CompleteIncomplete3DVerticesInGaps
  
  //////////////////////////////////////////
  void CompleteIncomplete3DVertices(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // Look for trajectories in a plane that lack a 2D vertex as listed in
    // 2DVtxID that are near the projected wire. This may trigger splitting trajectories,
    // assigning them to a new 2D vertex and completing 3D vertices
    
    if(!tjs.UseAlg[kComp3DVx]) return;
    if(tjs.NumPlanes != 3) return;
    
    bool prt = (debug.Plane >= 0 && debug.Tick == 33333);
    if(prt) mf::LogVerbatim("TC")<<"Inside CI3DV";
    
    float maxdoca = 6;
    unsigned short ivx3 = 0;
    for(auto& vx3 : tjs.vtx3) {
      // ignore obsolete vertices
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      // check for a completed 3D vertex
      if(vx3.Wire < 0) continue;
      unsigned short mPlane = USHRT_MAX;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(vx3.Vx2ID[plane] > 0) continue;
        mPlane = plane;
      } // ipl
      if(mPlane == USHRT_MAX) continue;
      CTP_t mCTP = EncodeCTP(vx3.TPCID.Cryostat, vx3.TPCID.TPC, mPlane);
      // X position of the purported missing vertex
      // A TP for the missing 2D vertex
      TrajPoint vtp;
      vtp.Pos[0] = vx3.Wire;
      vtp.Pos[1] = tjs.detprop->ConvertXToTicks(vx3.X, mPlane, vx3.TPCID.TPC, vx3.TPCID.Cryostat) * tjs.UnitsPerTick;
      if(prt) mf::LogVerbatim("TC")<<"CI3DV vx3.ID "<<vx3.ID<<" Pos "<<mPlane<<":"<<PrintPos(tjs, vtp.Pos);
      std::vector<int> tjIDs;
      std::vector<unsigned short> tjPts;
      for(auto& tj : tjs.allTraj) {
        if(tj.CTP != mCTP) continue;
        if(tj.AlgMod[kKilled]) continue;
        if(tj.Pts.size() < 6) continue;
        if(tj.AlgMod[kComp3DVx]) continue;
        float doca = maxdoca;
        // find the closest distance between the vertex and the trajectory
        unsigned short closePt = 0;
        TrajPointTrajDOCA(tjs, vtp, tj, closePt, doca);
        if(closePt > tj.EndPt[1]) continue;
        // try to improve the location of the vertex by looking for a distinctive feature on the
        // trajectory, e.g. high multiplicity hits or larger than normal charge
        if(RefineVtxPosition(tjs, tj, closePt, 3, false)) vtp.Pos = tj.Pts[closePt].Pos;
        if(prt) mf::LogVerbatim("TC")<<"CI3DV vx3.ID "<<vx3.ID<<" candidate itj ID "<<tj.ID<<" vtx pos "<<PrintPos(tjs, vtp.Pos)<<" doca "<<doca;
        tjIDs.push_back(tj.ID);
        tjPts.push_back(closePt);
      } // itj
      if(tjIDs.empty()) continue;
      // compare the length of the Tjs used to make the vertex with the length of the
      // Tj that we want to split. Don't allow a vertex using very short Tjs to split a long
      // Tj in the 3rd plane
      float score;
      auto vxtjs = GetVtxTjIDs(tjs, vx3, score);
      unsigned short maxPts = 0;
      for(auto tjid : vxtjs) {
        auto& tj = tjs.allTraj[tjid - 1];
        unsigned short npwc = NumPtsWithCharge(tjs, tj, false);
        if(npwc > maxPts) maxPts = npwc;
      } // tjid
      // skip this operation if any of the Tjs in the split list are > 3 * maxPts
      maxPts *= 3;
      bool skipit = false;
      for(auto tjid : tjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(NumPtsWithCharge(tjs, tj, false) > maxPts) skipit = true;
      } // tjid
      if(prt) mf::LogVerbatim("TC")<<"  maxPts "<<maxPts<<" vxtjs[0] "<<vxtjs[0]<<" skipit? "<<skipit;
      if(skipit) continue;
      // 2D vertex
      VtxStore aVtx;
      unsigned short newVtxIndx = tjs.vtx.size();
      aVtx.ID = newVtxIndx + 1;
      aVtx.CTP = mCTP;
      aVtx.Topo = 3;
      aVtx.NTraj = 0;
      // Give it a bogus pass to indicate it wasn't created while stepping
      aVtx.Pass = 9;
      aVtx.Pos = vtp.Pos;
      // ensure this isn't in a messy region
      aVtx.TjChgFrac = ChgFracNearPos(tjs, aVtx.Pos, tjIDs);
      if(prt) mf::LogVerbatim("TC")<<" charge fraction near position "<<aVtx.TjChgFrac;
      if(aVtx.TjChgFrac < 0.6) continue;
      if(!StoreVertex(tjs, aVtx)) continue;
      // make a reference to the new vertex
      VtxStore& newVtx = tjs.vtx[tjs.vtx.size()-1];
      if(prt) mf::LogVerbatim("TC")<<" Stored 2D vertex "<<newVtx.ID;
      // make a temporary copy so we can nudge it a bit if there is only one Tj
      std::array<float, 2> vpos = aVtx.Pos;
      for(unsigned short ii = 0; ii < tjIDs.size(); ++ii) {
        unsigned short itj = tjIDs[ii] - 1;
        auto& tj = tjs.allTraj[itj];
        unsigned short closePt = tjPts[ii];
        // determine which end is the closest
        unsigned short end = 1;
        // closest to the beginning?
        if(fabs(closePt - tjs.allTraj[itj].EndPt[0]) < fabs(closePt - tj.EndPt[1])) end = 0;
        short dpt = fabs(closePt - tj.EndPt[end]);
        if(dpt < 3) {
          // close to an end
          if(tj.VtxID[end] > 0) {
            if(prt) mf::LogVerbatim("TC")<<" Tj has a vertex "<<tj.VtxID[end]<<" at this end "<<end;
            continue;
          }
          tj.VtxID[end] = tjs.vtx[newVtxIndx].ID;
          ++newVtx.NTraj;
          if(prt) mf::LogVerbatim("TC")<<" attach Traj ID "<<tj.ID<<" to end "<<end;
          tj.AlgMod[kComp3DVx] = true;
          vpos = tj.Pts[tj.EndPt[end]].Pos;
        } else {
          // closePt is not near an end, so split the trajectory
          if(SplitAllTraj(tjs, itj, closePt, newVtxIndx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<" SplitAllTraj success "<<tjs.vtx[newVtxIndx].ID<<" at closePt "<<closePt;
            // successfully split the Tj
            newVtx.NTraj += 2;
          } else {
            // split failed. Give up
            if(prt) mf::LogVerbatim("TC")<<" SplitAllTraj failed";
            newVtx.ID = 0;
            break;
          }
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(tjs, itj);
          // and for the new trajectory
          SetPDGCode(tjs, tjs.allTraj.size()-1);
        } // closePt is not near an end, so split the trajectory
        tj.AlgMod[kComp3DVx] = true;
        itj = tjs.allTraj.size() - 1;
        tjs.allTraj[itj].AlgMod[kComp3DVx] = true;
      } // ii
      if(newVtx.NTraj == 0) {
        // A failure occurred. Recover
        if(prt) mf::LogVerbatim("TC")<<"  Failed. Recover and delete vertex "<<newVtx.ID;
        MakeVertexObsolete(tjs, newVtx, true);
      } else {
        // success
        vx3.Vx2ID[mPlane] = newVtx.ID;
        newVtx.Vx3ID = vx3.ID;
        vx3.Wire = -1;
        // set the vertex position to the start of the Tj if there is only one and fix it
        if(newVtx.NTraj == 1) {
          newVtx.Pos = vpos;
          newVtx.Stat[kFixed] = true;
        }
        AttachAnyTrajToVertex(tjs, newVtx.ID - 1, prt);
        SetVx2Score(tjs, prt);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<" Success: new 2D vtx ID "<<newVtx.ID<<" at "<<(int)newVtx.Pos[0]<<":"<<(int)newVtx.Pos[1]/tjs.UnitsPerTick;
          myprt<<" points to 3D vtx "<<vx3.ID;
          myprt<<" TjIDs:";
          for(auto& tjID : tjIDs) myprt<<" "<<tjID;
        } // prt
      } // success
      ++ivx3;
    } // vx3
    
  } // CompleteIncomplete3DVertices
  
  ////////////////////////////////////////////////
  bool RefineVtxPosition(TjStuff& tjs, const Trajectory& tj, unsigned short& nearPt, short nPtsToChk, bool prt)
  {
    // The tj has been slated to be split somewhere near point nearPt. This function will move
    // the near point a bit to the most likely point of a vertex
    
    float maxChg = tj.Pts[nearPt].Chg;
    short maxChgPt = nearPt;
    
    for(short ipt = nearPt - nPtsToChk; ipt < nearPt + nPtsToChk; ++ipt) {
      if(ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) continue;
      auto& tp = tj.Pts[ipt];
      if(tp.Chg > maxChg) {
        maxChg = tp.Chg;
        maxChgPt = ipt;
      }
      if(prt) mf::LogVerbatim("TC")<<"RVP: ipt "<<ipt<<" Pos "<<tp.CTP<<":"<<PrintPos(tjs, tp.Pos)<<" chg "<<(int)tp.Chg<<" nhits "<<tp.Hits.size();
    } // ipt
    if(nearPt == maxChgPt) return false;
    nearPt = maxChgPt;
    return true;
  } //RefineVtxPosition
  
  /////////////////////TY:///////////////////////////
  void VtxHitsSwap(TjStuff& tjs, const CTP_t inCTP){
    
    if(!tjs.UseAlg[kVtxHitsSwap]) return;
    
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    bool vtxPrt = ((debug.Plane == (int)planeID.Plane || debug.Plane == 3) && debug.Tick < 0);
    if(vtxPrt) mf::LogVerbatim("TC")<<"inside VtxHitsSwap on plane "<<planeID.Plane;
    for (unsigned short iv = 0; iv < tjs.vtx.size(); ++iv){
      VtxStore& rvx = tjs.vtx[iv];
      if(rvx.CTP != inCTP) continue;
      // Only consider vertex with two trajectories
      if(rvx.NTraj != 2) continue;
      if (vtxPrt) mf::LogVerbatim("TC")<<"Vertex "<<iv<<" Pos[0]: "<<rvx.Pos[0]<<" "<<rvx.Pos[1];
      std::array<unsigned short, 2> tjlist{{0,0}};
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
        if(tjs.allTraj[itj].CTP != rvx.CTP) continue;
        Trajectory& tj = tjs.allTraj[itj];
        // ensure that the ID is OK so the code below doesn't choke
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == rvx.ID) {
            tjlist[end] = itj;
          }
        }
      }//all trajectories
      
      //Ignore short trajectories
      if (tjs.allTraj[tjlist[0]].EndPt[1]<5||
          tjs.allTraj[tjlist[1]].EndPt[1]<5) continue;
      
      for (unsigned short i = 0; i<2; ++i){
        
        //First check if first hit should be swapped
        Trajectory& tj0 = tjs.allTraj[tjlist[i]];
        Trajectory& tj1 = tjs.allTraj[tjlist[1-i]];
        unsigned short endPt0 = tjs.allTraj[tjlist[i]].EndPt[i];
        unsigned short endPt1 = tjs.allTraj[tjlist[1-i]].EndPt[1-i];
        //first TP on first trajectory
        float chg0 = TpSumHitChg(tjs, tj0.Pts[endPt0]);
        float w0 = tj0.Pts[endPt0].Pos[0];
        //if (vtxPrt) mf::LogVerbatim("TC")<<PrintPos(tjs, tj0.Pts[endPt0]);
        //second TP on first trajectory
        float chg1 = 0;
        float w1 = 0;
        unsigned short j = endPt0;
        while (j!=tj0.EndPt[1-i]){
          if (i==0) ++j;
          else --j;
          if (tj0.Pts[j].Chg){
            chg1 = TpSumHitChg(tjs, tj0.Pts[j]);
            w1 = tj0.Pts[j].Pos[0];
            //if (vtxPrt) mf::LogVerbatim("TC")<<PrintPos(tjs, tj0.Pts[j]);
            break;
          }
        }
        //first TP on second trajectory
        float chg2 = TpSumHitChg(tjs,tj1.Pts[endPt1]);
        float w2 = tj1.Pts[endPt1].Pos[0];
        //DOCA between first TP on first TJ and first TP on second TJ
        float delta = 1000;
        for (size_t k = 0; k<tj0.Pts[endPt0].Hits.size(); ++k){
          if (!tj0.Pts[endPt0].UseHit[k]) continue;
          float this_delta = PointTrajDOCA(tjs, tj0.Pts[endPt0].Hits[k], tj1.Pts[endPt1]);
          if (this_delta<delta) delta = this_delta;
        }
        //        if (vtxPrt) mf::LogVerbatim("TC")<<PrintPos(tjs, tj1.Pts[endPt1]);
        //if (vtxPrt) mf::LogVerbatim("TC")<<chg0<<" "<<chg1<<" "<<chg2<<" "<<delta;
        if (chg0>0&&std::abs((chg0-chg1)/chg0)-std::abs((chg0-chg2)/chg0)>0.2&&delta<1.5&&std::abs(w2-w1)<1.5){
          if (vtxPrt) {
            mf::LogVerbatim("TC")<<"chg0 = "<<chg0<<" chg1 = "<<chg1<<" chg2 = "<<chg2<<" delta = "<<delta<<" w0 = "<<w0<<" w1 = "<<w1<<" w2 = "<<w2;
            mf::LogVerbatim("TC")<<"VHS Moving TP "<<PrintPos(tjs, tj0.Pts[endPt0])<<" from TJ "<<tj0.ID<<" to TJ "<<tj1.ID;
          }
          //Add first TP of first trajectory to the second trajectory
          TrajPoint tp = tj0.Pts[endPt0];
          for (size_t j = 0; j<tp.Hits.size(); ++j){
            if (!tp.UseHit[j]) continue;
            tjs.fHits[tp.Hits[j]].InTraj = tj1.ID;
          }
          if (i==0){
            //append to the end
            tj1.Pts.push_back(tp);
          }
          else{
            //insert at the beginning
            tj1.Pts.insert(tj1.Pts.begin(), tp);
          }
          SetEndPoints(tjs, tj1); 
          
          //Remove first TP from first trajectory
          tj0.Pts[endPt0].Chg = 0;
          for (size_t j = 0; j<tj0.Pts[endPt0].Hits.size(); ++j){
            tj0.Pts[endPt0].UseHit[j] = false;
          }
          SetEndPoints(tjs, tj0);
          tj0.AlgMod[kVtxHitsSwap] = true;
          tj1.AlgMod[kVtxHitsSwap] = true;
          break;
        }
        
        //Now Check if the beginning of the first trajectory should be moved to the second trajectory.
        j = endPt0;
        std::vector<unsigned short> tplist;
        while (j!=tj0.EndPt[1-i]){
          if (tj0.Pts[j].Chg){
            float delta = 1000;
            for (size_t k = 0; k<tj0.Pts[j].Hits.size(); ++k){
              if (!tj0.Pts[j].UseHit[k]) continue;
              float this_delta = PointTrajDOCA(tjs, tj0.Pts[j].Hits[k], tj1.Pts[endPt1]);
              if (this_delta<delta) delta = this_delta;
              //if (vtxPrt) mf::LogVerbatim("TC")<<j<<" "<<k<<" "<<PrintPos(tjs, tj0.Pts[j])<<" "<<PointTrajDOCA(tjs, tj0.Pts[j].Hits[k], tj1.Pts[endPt1]);
            }
            if (delta < 0.3 && tj0.Pts[j].Delta > 1.0 && (j==endPt0 || !tplist.empty())){
              tplist.push_back(j);
            }
            else break;
          }
          if (i==0) ++j;
          else --j;
        }
        //if (vtxPrt) mf::LogVerbatim("TC")<<tplist.size();
        //Need at least two TPs to make this change
        if (tplist.size()>1){
          if (vtxPrt) mf::LogVerbatim("TC")<<"VHS Moving "<<tplist.size()<<" TPs from TJ "<<tj0.ID<<" to TJ "<<tj1.ID;
          //Append TPs to the second TJ
          for (unsigned short j = 0; j<tplist.size(); ++j){
            TrajPoint tp = tj0.Pts[tplist[j]];
            for (size_t k = 0; k<tp.Hits.size(); ++k){
              if (!tp.UseHit[k]) continue;
              tjs.fHits[tp.Hits[k]].InTraj = tj1.ID;
            }
            if (i==0){
              //append to the end
              tj1.Pts.push_back(tp);
            }
            else{
              //insert at the beginning
              tj1.Pts.insert(tj1.Pts.begin(), tp);
            }
          }
          SetEndPoints(tjs, tj1); 
          
          //Remove TPs from first trajectory
          for (unsigned short j = 0; j<tplist.size(); ++j){
            tj0.Pts[tplist[j]].Chg = 0;
            for (size_t k = 0; k<tj0.Pts[tplist[j]].Hits.size(); ++k){
              tj0.Pts[tplist[j]].UseHit[k] = false;
            }
          }
          SetEndPoints(tjs, tj0);
          tj0.AlgMod[kVtxHitsSwap] = true;
          tj1.AlgMod[kVtxHitsSwap] = true;
          break;
        }
      }//loop over two trajectories
    }//loop over vertices
  }
  
  ////////////////////////////////////////////////
  bool MakeVertexObsolete(TjStuff& tjs, VtxStore& vx2, bool forceKill)
  {
    // Makes a 2D vertex obsolete
    
    // check for a high-score 3D vertex
    bool hasHighScoreVx3 = (vx2.Vx3ID > 0);
    if(hasHighScoreVx3 && !forceKill && tjs.vtx3[vx2.Vx3ID - 1].Score >= tjs.Vertex2DCuts[7]) return false;
    
    // Kill it
    unsigned short vx2id = vx2.ID;
    vx2.ID = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] != vx2id) continue;
        tj.VtxID[end] = 0;
        if(tj.AlgMod[kTjHiVx3Score]) {
          // see if the vertex at the other end is high-score and if so, preserve the state
          unsigned short oend = 1 - end;
          if(tj.VtxID[oend] > 0) {
            auto& ovx2 = tjs.vtx[tj.VtxID[oend] - 1];
            if(!ovx2.Stat[kHiVx3Score]) tj.AlgMod[kTjHiVx3Score] = false;
          } // vertex at the other end
        } // tj.AlgMod[kTjHiVx3Score]
      } // end
    } // tj
    
    if(!hasHighScoreVx3) return true;
    
    // update the affected 3D vertex
    Vtx3Store& vx3 = tjs.vtx3[vx2.Vx3ID - 1];
    // make the 3D vertex incomplete
    geo::PlaneID planeID = DecodeCTP(vx2.CTP);
    unsigned short plane = planeID.Plane;
    if(vx3.Vx2ID[plane] != vx2id) return true;
    vx3.Vx2ID[plane] = 0;
    vx3.Wire = vx2.Pos[0] / tjs.UnitsPerTick;
    // Ensure that there are at least two 2D vertices left
    unsigned short n2D = 0;
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
      if(vx3.Vx2ID[plane] > 0) ++n2D;
    } // plane
    
    if(n2D > 1) {
      // 3D vertex is incomplete
      // correct the score
      SetVx3Score(tjs, vx3, false);
      return true;
    }
    
    // 3D vertex is obsolete
    // Detach the all remaining 2D vertices from the 3D vertex
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      if(vx2.Vx3ID == vx3.ID) vx2.Vx3ID = 0;
    } // vx2
    vx3.ID = 0;
    return true;
    
  } // MakeVertexObsolete
/*
  ////////////////////////////////////////////////
  bool MakeVertexObsolete(TjStuff& tjs, unsigned short vx2ID, bool forceKill)
  {
    // Deletes a 2D vertex and possibly a 3D vertex and 2D vertices in other planes
    // The 2D and 3D vertices are NOT killed if forceKill is false and the 3D vertex
    // has a high score
    if(vx2ID == 0) return true;
    if(vx2ID > tjs.vtx.size()) return false;
    VtxStore& vx2 = tjs.vtx[vx2ID - 1];
    
    return MakeVertexObsolete(tjs, vx2, forceKill);
   
  } // MakeVertexObsolete
*/
  //////////////////////////////////////////
  std::vector<int> GetVtxTjIDs(const TjStuff& tjs, const VtxStore& vx2)
  {
    // returns a list of trajectory IDs that are attached to vx2
    std::vector<int> tmp;
    if(vx2.ID == 0) return tmp;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx2.CTP) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == vx2.ID) tmp.push_back(tj.ID);
      } // end
    } // tj
    return tmp;
  } // GetVtxTjIDs
  
  
  //////////////////////////////////////////
  std::vector<int> GetVtxTjIDs(const TjStuff& tjs, const Vtx3Store& vx3, float& score)
  {
    // returns a list of Tjs in all planes that are attached to vx3
    std::vector<int> tmp;
    if(vx3.ID == 0) return tmp;
    float nvx2 = 0;
    score = 0;
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      if(vx2.Vx3ID != vx3.ID) continue;
      auto vtxTjID2 = GetVtxTjIDs(tjs, vx2);
      tmp.insert(tmp.end(), vtxTjID2.begin(), vtxTjID2.end());
      score += vx2.Score;
      ++nvx2;
    } // vx2
    if(nvx2 < 1) return tmp;
    // find the average score
    score /= nvx2;
    // sort by increasing ID
    std::sort(tmp.begin(), tmp.end());
    return tmp;
  } // GetVtxTjIDs

  //////////////////////////////////////////
  void PosInPlane(const TjStuff& tjs, const Vtx3Store& vx3, unsigned short plane, std::array<float, 2>& pos)
  {
    // returns the 2D position of the vertex in the plane
    pos[0] = tjs.geom->WireCoordinate(vx3.Y, vx3.Z, plane, vx3.TPCID.TPC, vx3.TPCID.Cryostat);
    pos[1] = tjs.detprop->ConvertXToTicks(vx3.X, plane, vx3.TPCID.TPC, vx3.TPCID.Cryostat) * tjs.UnitsPerTick;
    
  } // PosInPlane

  
  /////////////////////////////////////////
  unsigned short IsCloseToVertex(TjStuff& tjs, VtxStore& inVx2)
  {
    // Returns the ID of a 2D vertex having the minimum pull < user-specified cut
    
    float minPull = tjs.Vertex2DCuts[3];
    unsigned short imBest = 0;
    for(auto& vx2 : tjs.vtx) {
      float pull = VertexVertexPull(tjs, inVx2, vx2);
      if(pull < minPull) {
        minPull = pull;
        imBest = vx2.ID;
      }
    } // vx2
    return imBest;
  } // IsCloseToVertex
  
  /////////////////////////////////////////
  unsigned short IsCloseToVertex(TjStuff& tjs, Vtx3Store& vx3)
  {
    // Returns the ID of a 3D vertex having the minimum pull < user-specified cut
    
    float minPull = tjs.Vertex3DCuts[1];
    unsigned short imBest = 0;
    for(auto& oldvx3 : tjs.vtx3) {
      if(oldvx3.ID == 0) continue;
      if(std::abs(oldvx3.X - vx3.X) > tjs.Vertex3DCuts[0]) continue;
      float pull = VertexVertexPull(tjs, vx3, oldvx3);
      if(pull < minPull) {
        minPull = pull;
        imBest = oldvx3.ID;
      }
    } // oldvx3
    return imBest;
    
  } // IsCloseToVertex
  
/*
   //////////////////////////////////////////
   void Refine2DVertices()
   {
   // Improve trajectories near vertices in the current plane
   if(tjs.vtx.empty()) return;
   
   if(!tjs.UseAlg[kRefineVtx]) return;
   
   geo::PlaneID planeID = DecodeCTP(inCTP);
   unsigned short ipl = planeID.Plane;
   if(ipl != 0) return;
   
   // store a copy of everything so that we can recover gracefully if there is a major failure
   auto lHits = tjs.fHits;
   auto lWireHitRange = tjs.WireHitRange;
   auto lAllTraj = tjs.allTraj;
   bool majorFailure = false;
   
   if(prt) PrintHeader("R2D");
   
   for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
   VtxStore& rvx = tjs.vtx[ivx];
   if(rvx.CTP != inCTP) continue;
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
   if(prt) PrintTrajectory("R2D", tjs, tj, endPt);
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
   std::array<int, 2> iwROC {(int)loWire, (int)hiWire};
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
   */

  
} // namespace
