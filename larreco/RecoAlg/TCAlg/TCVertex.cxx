#include "larreco/RecoAlg/TCAlg/TCVertex.h"

namespace tca {
  
  struct SortEntry{
    unsigned int index;
    float val;
  };
  
  bool valDecreasing (SortEntry c1, SortEntry c2) { return (c1.val > c2.val);}
  bool valIncreasing (SortEntry c1, SortEntry c2) { return (c1.val < c2.val);}
  
  //////////////////////////////////////////
  void MakeJunkVertices(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Vertices between poorly reconstructed tjs (especially junk tjs) and normal
    // tjs can fail because the junk tj trajectory parameters are inaccurate. This function
    // uses proximity and not pointing to make junk vertices
    // Don't use this if standard vertex reconstruction is disabled
    if(tjs.Vertex2DCuts[0] <= 0) return;
    if(!tjs.UseAlg[kJunkVx]) return;
    if(tjs.allTraj.size() < 2) return;
    
    // Look for tjs that are within maxSep of the end of a Tj
    constexpr float maxSep = 4;
    
    geo::PlaneID planeID = DecodeCTP(inCTP);
    bool prt = (debug.Plane == (int)planeID.Plane && debug.Tick == 99999);
    if(prt) {
      mf::LogVerbatim("TC")<<"MakeJunkVertices: prt set for plane "<<planeID.Plane<<" maxSep btw tjs "<<maxSep;
//      PrintAllTraj("MJTi", tjs, debug, USHRT_MAX, tjs.allTraj.size());
    }
    
    // make a template vertex
    VtxStore junkVx;
    junkVx.CTP = inCTP;
    junkVx.Topo = 9;
    junkVx.Stat[kJunkVx] = true;
    junkVx.Stat[kFixed] = true;
    // set an invalid ID
    junkVx.ID = USHRT_MAX;
    // put in generous errors
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    junkVx.PosErr = {{2.0, 2.0}};
    // define a minimal score so it won't get clobbered
    junkVx.Score = tjs.Vertex2DCuts[7] + 0.1;

    // look at both ends of long tjs
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size() - 1; ++it1) {
      auto& tj1 = tjs.allTraj[it1];
      if(tj1.AlgMod[kKilled]) continue;
      if(tj1.AlgMod[kInShower]) continue;
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kJunkTj]) continue;
      if(TrajLength(tj1) < 10) continue;
      if(tj1.MCSMom < 100) continue;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // existing vertex?
        if(tj1.VtxID[end1] > 0) continue;
        auto& tp1 = tj1.Pts[tj1.EndPt[end1]];
        // get a list of tjs in this vicinity
        auto tjlist = FindCloseTjs(tjs, tp1, tp1, maxSep);
        if(tjlist.empty()) continue;
        // set to an invalid ID
        junkVx.ID = USHRT_MAX;
        for(auto tj2id : tjlist) {
          auto& tj2 = tjs.allTraj[tj2id - 1];
          if(tj2.CTP != inCTP) continue;
          if(tj2id == tj1.ID) continue;
//          if(tj2.MCSMom > 50) continue;
          if(tj2.AlgMod[kInShower]) continue;
//          if(tj2.Pts.size() > 20) continue;
          float close = maxSep;
          unsigned short closeEnd = USHRT_MAX;
          for(unsigned short end2 = 0; end2 < 2; ++end2) {
            auto& tp2 = tj2.Pts[tj2.EndPt[end2]];
            float sep = PosSep(tp1.Pos, tp2.Pos);
            if(sep < close) {
              close = sep;
              closeEnd = end2;
            } // sep
          } // end2
          if(closeEnd > 1) continue;
          auto& tp2 = tj2.Pts[tj2.EndPt[closeEnd]];
          bool signalBetween = SignalBetween(tjs, tp1, tp2, 0.8, prt);
          if(!signalBetween) continue;
          if(junkVx.ID == USHRT_MAX) {
            // define the new vertex
            junkVx.ID = tjs.vtx.size() + 1;
            junkVx.Pos = tp1.Pos;
          } // new vertex
          tj2.VtxID[closeEnd] = junkVx.ID;
          tj1.VtxID[end1] = junkVx.ID;
        } // tjid
        if(junkVx.ID == USHRT_MAX) continue;
        if(!StoreVertex(tjs, junkVx)) {
          mf::LogVerbatim("TC")<<"MJV: StoreVertex failed";
          for(auto& tj : tjs.allTraj) {
            if(tj.AlgMod[kKilled]) continue;
            if(tj.VtxID[0] == junkVx.ID) tj.VtxID[0] = 0;
            if(tj.VtxID[1] == junkVx.ID) tj.VtxID[1] = 0;
          } // tj
          continue;
        } // StoreVertex failed
        if(prt) {
          mf::LogVerbatim("TC")<<" New junk 2V"<<junkVx.ID<<" at "<<std::fixed<<std::setprecision(1)<<junkVx.Pos[0]<<":"<<junkVx.Pos[1]/tjs.UnitsPerTick;
        } // prt
        junkVx.ID = USHRT_MAX;
      } // end1
    } // it1
    
  } // MakeJunkVertices

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
    
    unsigned short maxShortTjLen = tjs.Vertex2DCuts[0];
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size() - 1; ++it1) {
      auto& tj1 = tjs.allTraj[it1];
      if(tj1.AlgMod[kKilled]) continue;
      if(tj1.AlgMod[kInShower]) continue;
      if(tj1.CTP != inCTP) continue;
      bool tj1Short = (TrajLength(tj1) < maxShortTjLen);
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tj1.VtxID[end1] > 0) continue;
        // default condition is to use the end point to define the trajectory and direction
        // at the end
        short endPt1 = tj1.EndPt[end1];
        float wire1 = tj1.Pts[endPt1].Pos[0];
        // unless there are few points fitted, indicating that the trajectory fit
        // may have been biased by the presence of another trajectory at the vertex or by
        // other close unresolved tracks
        if(tj1.Pts.size() > 6 && tj1.Pts[endPt1].NTPsFit < 4) {
          if(end1 == 0 && endPt1 < int(tj1.Pts.size()) - 3) {
            endPt1 += 3;
          } else if (end1 == 1 && endPt1 >=3 ) {
            endPt1 -= 3;
          }
          if(tj1.Pts[endPt1].Chg == 0) endPt1 = NearestPtWithChg(tjs, tj1, endPt1);
        } // few points fit at end1
        TrajPoint tp1 = tj1.Pts[endPt1];
        MoveTPToWire(tp1, wire1);
        // re-purpose endPt1 to reference the end point. This will be used the find the point on
        // tj1 that is closest to the vertex position
        endPt1 = tj1.EndPt[end1];
        short oendPt1 = tj1.EndPt[1-end1];
        // reference to the other end of tj1
        auto& otp1 = tj1.Pts[oendPt1];
        for(unsigned short it2 = it1 + 1; it2 < tjs.allTraj.size(); ++it2) {
          auto& tj2 = tjs.allTraj[it2];
          if(tj2.AlgMod[kKilled]) continue;
          if(tj2.AlgMod[kInShower]) continue;
          if(tj2.CTP != inCTP) continue;
          if(tj1.VtxID[end1] > 0) continue;
          if(tj1.MCSMom < tjs.Vertex2DCuts[5] && tj2.MCSMom < tjs.Vertex2DCuts[5]) continue;
          bool tj2Short = (TrajLength(tj2) < maxShortTjLen);
          // find the end that is closer to tp1
          unsigned short end2 = 0;
          if(PosSep2(tj2.Pts[tj2.EndPt[1]].Pos, tp1.Pos) < PosSep2(tj2.Pts[tj2.EndPt[0]].Pos, tp1.Pos)) end2 = 1;
          if(tj2.VtxID[end2] > 0) continue;
          // check for a vertex between these tjs at the other ends
          if(tj1.VtxID[1 - end1] > 0 && tj1.VtxID[1 - end1] == tj2.VtxID[1 - end2]) continue;
          // see if the other ends are closer 
          unsigned short oendPt2 = tj2.EndPt[1-end2];
          auto& otp2 = tj2.Pts[oendPt2];
          if(PosSep2(otp1.Pos, otp2.Pos) < PosSep2(tp1.Pos, tj2.Pts[tj2.EndPt[end2]].Pos)) continue;
          short endPt2 = tj2.EndPt[end2];
          float wire2 = tj2.Pts[endPt2].Pos[0];
          if(tj2.Pts.size() > 6 && tj2.Pts[endPt2].NTPsFit < 4) {
            if(end2 == 0 && endPt2 < int(tj2.Pts.size()) - 3) {
              endPt2 += 3;
            } else if (end2 == 1 && endPt2 >= 3){
              endPt2 -= 3;
            }
            if(tj2.Pts[endPt2].Chg == 0) endPt2 = NearestPtWithChg(tjs, tj2, endPt2);
          } // few points fit at end1
          TrajPoint tp2 = tj2.Pts[endPt2];
          MoveTPToWire(tp2, wire2);
          // re-purpose endPt2
          endPt2 = tj2.EndPt[end2];
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
          // BUG the double brace syntax is required to work around clang bug 21629
          // (https://bugs.llvm.org/show_bug.cgi?id=21629)
          Point2_t vPos {{wint, tint}};
          float vt1Sep = PosSep(vPos, tp1.Pos);
          float vt2Sep = PosSep(vPos, tp2.Pos);
          float dwc1 = DeadWireCount(tjs, wint, tp1.Pos[0], tp1.CTP);
          float dwc2 = DeadWireCount(tjs, wint, tp2.Pos[0], tp1.CTP);
          vt1Sep -= dwc1;
          vt2Sep -= dwc2;
          bool vtxOnDeadWire = (DeadWireCount(tjs, wint, wint, tp1.CTP) == 1);            
          if(prt && vt1Sep < 200 && vt2Sep < 200) {
            mf::LogVerbatim myprt("TC");
            myprt<<"F2DV candidate T"<<tj1.ID<<"_"<<end1<<"-T"<<tj2.ID<<"_"<<end2;
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
          float doca1 = sepCut;
          if(!TrajClosestApproach(tj1, wint, tint, closePt1, doca1)) continue;
          // dpt1 (and dpt2) will be 0 if the vertex is at the end
          short dpt1 = tjs.StepDir * (closePt1 - endPt1);
          if(prt) mf::LogVerbatim("TC")<<" endPt1 "<<endPt1<<" closePt1 "<<closePt1<<" dpt1 "<<dpt1<<" doca1 "<<doca1;
          // BB April 19, 2018: require vertex to be near the end
          if(dpt1 < -1) continue;
          if(tjs.allTraj[it1].EndPt[1] > 4) {
            if(dpt1 > 3) continue;
          } else {
            // tighter cut for short trajectories
            if(dpt1 > 2) continue;
          }
          unsigned short closePt2;
          float doca2 = sepCut;
          if(!TrajClosestApproach(tj2, wint, tint, closePt2, doca2)) continue;
          short dpt2 = tjs.StepDir * (closePt2 - endPt2);
          if(prt) mf::LogVerbatim("TC")<<" endPt2 "<<endPt2<<" closePt2 "<<closePt2<<" dpt2 "<<dpt2<<" doca2 "<<doca2;
          // BB April 19, 2018: require vertex to be near the end
          if(dpt2 < -1) continue;
          if(tjs.allTraj[it2].EndPt[1] > 4) {
            if(dpt2 > 3) continue;
          } else {
            // tighter cut for short trajectories
            if(dpt2 > 2) continue;
          }
          if(prt) mf::LogVerbatim("TC")<<" wint:tint "<<(int)wint<<":"<<(int)(tint/tjs.UnitsPerTick);
          // ensure that there is a signal between these TPs and the vertex on most of the wires
          bool signalBetween = true;
          bool fixVxPos = false;
          short dpt = abs(wint - tp1.Pos[0]);
          if(dpt > 2 && !SignalBetween(tjs, tp1, wint, tjs.Vertex2DCuts[6], prt)) {
            if(prt) mf::LogVerbatim("TC")<<" Fails SignalBetween for tp1 "<<dpt;
            signalBetween = false;
          }
          dpt = abs(wint - tp2.Pos[0]);
          if(dpt > 2 && !SignalBetween(tjs, tp2, wint, tjs.Vertex2DCuts[6], prt)) {
            if(prt) mf::LogVerbatim("TC")<<" Fails SignalBetween for tp2 "<<dpt;
            signalBetween = false;
          }
          // consider the case where the intersection point is wrong because the
          // end TP angles are screwed up but the Tjs are close to each other near the end
          if(!signalBetween) {
            unsigned short ipt1, ipt2;
            float maxSep = 3;
            bool isClose = TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, maxSep, false);
            if(prt) mf::LogVerbatim("TC")<<" TrajTrajDOCA close? "<<isClose<<" minSep "<<maxSep;
            if(isClose) {
              // put the vertex at the TP that is closest to the intersection point
              Point2_t vpos = {{wint, tint}};
              if(PosSep2(tp1.Pos, vpos) < PosSep2(tp2.Pos, vpos)) {
                wint = tp1.Pos[0];
                tint = tp1.Pos[1];
              } else {
                wint = tp2.Pos[0];
                tint = tp2.Pos[1];
              }
              fixVxPos = true;
              if(prt) mf::LogVerbatim("TC")<<" new wint:tint "<<(int)wint<<":"<<(int)(tint/tjs.UnitsPerTick);
            } else {
              // closest approach > 3
              continue;
            }
          } // no signal between
          // make a new temporary vertex
          VtxStore aVtx;
          aVtx.Pos[0] = wint;
          aVtx.Pos[1] = tint;
          aVtx.NTraj = 0;
          aVtx.Pass = tj1.Pass;
          // Topo 0 has this topology (<) and Topo 2 has this (>)
          aVtx.Topo = 2 * end1;
          aVtx.ChiDOF = 0;
          aVtx.CTP = inCTP;
          aVtx.Stat[kOnDeadWire] = vtxOnDeadWire;
          // fix the vertex position if we needed to move it significantly, or if it is on a dead wire
          aVtx.Stat[kFixed] = fixVxPos;
//          if(close2 > 1) aVtx.Stat[kFixed] = true;
          // try to fit it. We need to give it an ID to do that. Take the next
          // available ID
          unsigned short newVtxID = tjs.vtx.size() + 1;
          aVtx.ID = newVtxID;
          tj1.VtxID[end1] = newVtxID;
          tj2.VtxID[end2] = newVtxID;
          if(!FitVertex(tjs, aVtx, prt)) {
            tj1.VtxID[end1] = 0;
            tj2.VtxID[end2] = 0;
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
            myprt<<" New vtx 2V"<<aVtx.ID;
            myprt<<" Tjs "<<tj1.ID<<"_"<<end1<<"-"<<tj2.ID<<"_"<<end2;
            myprt<<" at "<<std::fixed<<std::setprecision(1)<<aVtx.Pos[0]<<":"<<aVtx.Pos[1]/tjs.UnitsPerTick;
          }
          AttachAnyTrajToVertex(tjs, tjs.vtx.size() - 1, prt);
          SetVx2Score(tjs, prt);
        } // it2
      } // end1
    } // it1
    
    // check the consistency of the Tjs for the newly added vertices
    ChkVxTjs(tjs, inCTP, prt);
    
    // Split trajectories that cross a vertex
    // BB April 20, 2018: move to ReconstructAllTraj
//    SplitTrajCrossingVertices(tjs, inCTP);
    FindHammerVertices(tjs, inCTP);
    FindHammerVertices2(tjs, inCTP);
    
    if(prt) PrintAllTraj("F2DVo", tjs, debug, USHRT_MAX, USHRT_MAX);
    
  } // Find2DVertices

  //////////////////////////////////////////
  void FindNeutralVertices(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // Look for 2D neutral vertices between Tjs 
    if(!tjs.UseAlg[kVxNeutral]) return;
    if(tjs.NeutralVxCuts.size() < 4) return;
    if(tjs.NumPlanes < 3) return;
    if(tjs.pfps.size() < 2) return;
    
    bool prt = (debug.Plane >= 0 && debug.Tick == 88888);
    
    struct CandVx {
      unsigned short ip1;
      unsigned short end1;
      unsigned short ip2;
      unsigned short end2;
      Point3_t intersect;
      float sepSum;
      bool isValid;
    };
    std::vector<CandVx> candVxs;

    for(unsigned short ip1 = 0; ip1 < tjs.pfps.size() - 1; ++ip1) {
      auto& p1 = tjs.pfps[ip1];
      if(p1.ID == 0) continue;
      if(p1.Tp3s.empty()) continue;
      if(p1.TPCID != tpcid) continue;
      float len1 = PosSep(p1.XYZ[0], p1.XYZ[1]);
      if(len1 < tjs.NeutralVxCuts[3]) continue;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        float cfne1 = ChgFracNearEnd(tjs, p1, end1);
        if(cfne1 < tjs.NeutralVxCuts[0]) continue;
        for(unsigned short ip2 = ip1 + 1; ip2 < tjs.pfps.size(); ++ip2) {
          auto& p2 = tjs.pfps[ip2];
          if(p2.ID == 0) continue;
          if(p2.Tp3s.empty()) continue;
          if(p2.TPCID != tpcid) continue;
          float len2 = PosSep(p2.XYZ[0], p2.XYZ[1]);
          if(len2 < tjs.NeutralVxCuts[3]) continue;
          for(unsigned short end2 = 0; end2 < 2; ++end2) {
            float cfne2 = ChgFracNearEnd(tjs, p2, end2);
            if(cfne2 < tjs.NeutralVxCuts[0]) continue;
            float vxDOCA = 1E6;
            Point3_t intersect;
            if(!PointDirIntersect(p1.XYZ[end1], p1.Dir[end1], p2.XYZ[end2], p2.Dir[end2], intersect, vxDOCA)) continue;
            if(intersect[0] < tjs.XLo || intersect[0] > tjs.XHi) continue;
            if(intersect[1] < tjs.YLo || intersect[1] > tjs.YHi) continue;
            if(intersect[2] < tjs.ZLo || intersect[2] > tjs.ZHi) continue;
            // ensure that the pfp end and the vertex are consistent
            float sep1 = PosSep(intersect, p1.XYZ[end1]);
            if(PosSep(intersect, p1.XYZ[1-end1]) < sep1) continue;
            float sep2 = PosSep(intersect, p2.XYZ[end2]);
            if(PosSep(intersect, p2.XYZ[1-end2]) < sep2) continue;
            if(vxDOCA > tjs.NeutralVxCuts[1]) continue;
            // find the DOCA between these two
            unsigned short closePt1, closePt2;
            float pfpDOCA = PFPDOCA(p1, p2, closePt1, closePt2);
            if(closePt1 == USHRT_MAX) continue;
            // ensure that there isn't a lot of charge between the end of each pfp and the intersection point
            float cfb1 = ChgFracBetween(tjs, p1.XYZ[end1], intersect, p1.TPCID);
            if(cfb1 > tjs.NeutralVxCuts[2]) continue;
            float cfb2 = ChgFracBetween(tjs, p2.XYZ[end2], intersect, p2.TPCID);
            if(cfb2 > tjs.NeutralVxCuts[2]) continue;
            // check existing candidates
            float sepSum = sep1 + sep2;
            bool skipit = false;
            for(auto& candVx : candVxs) {
              if(!candVx.isValid) continue;
              if(candVx.ip1 != ip1 && candVx.ip2 != ip2) continue;
              // see if the separation sum is smaller
              if(sepSum < candVx.sepSum) {
                // flag the saved one as not valid
                candVx.isValid = false;
              } else {
                skipit = true;
                break;
              }
            } // candVx
            if(skipit) continue;
            CandVx candVx;
            candVx.ip1 = ip1;
            candVx.end1 = end1;
            candVx.ip2 = ip2;
            candVx.end2 = end2;
            candVx.intersect = intersect;
            candVx.sepSum = sep1 + sep2;
            candVx.isValid = true;
            candVxs.push_back(candVx);
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<"FNV: P"<<p1.ID<<"_"<<end1<<" sep1 "<<std::fixed<<std::setprecision(2)<<sep1;
              myprt<<" cfne1 "<<cfne1<<" cfb1 "<<cfb1;
              myprt<<" P"<<p2.ID<<"_"<<end2<<" sep2 "<<sep2<<" cfne2 "<<cfne2<<" cfb2 "<<cfb2;
              myprt<<" intersect "<<intersect[0]<<" "<<intersect[1]<<" "<<intersect[2];
              myprt<<" vxDOCA "<<vxDOCA<<" pfpDOCA "<<pfpDOCA;
            } // prt
          } // end2
        } // ip2
      } // end1
    } // ip1
    
    if(candVxs.empty()) return;
    
    // Make vertices with the valid candidates
    for(auto& candVx : candVxs) {
      if(!candVx.isValid) continue;
      Vtx3Store vx3;
      auto& p1 = tjs.pfps[candVx.ip1];
      auto& p2 = tjs.pfps[candVx.ip2];
      vx3.TPCID = p1.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = candVx.intersect[0];
      vx3.Y = candVx.intersect[1];
      vx3.Z = candVx.intersect[2];
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = true;
      tjs.vtx3.push_back(vx3);
      p1.Vx3ID[candVx.end1] = vx3.ID;
      p2.Vx3ID[candVx.end2] = vx3.ID;
      if(prt) mf::LogVerbatim("TC")<<"FNV: P"<<p1.ID<<"_"<<candVx.end1<<" P"<<p2.ID<<"_"<<candVx.end2<<" 3V"<<vx3.ID;
    } // candVx
    
  } // FindNeutralVertices

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
      unsigned short npwc = NumPtsWithCharge(tjs, tj, false);
      unsigned short end = CloseEnd(tjs, tj, vpos);
      // assume that we will use the end point of the tj
      unsigned short endPt = tj.EndPt[end];
      if(npwc > 6 && tj.Pts[endPt].NTPsFit < 4) {
        if(end == 0) {
          endPt += 3;
        } else {
          endPt -= 3;
        }
        endPt = NearestPtWithChg(tjs, tj, endPt);
      } // few points fit at the end
      if(endPt < tj.EndPt[0]) endPt = tj.EndPt[0];
      if(endPt > tj.EndPt[1]) endPt = tj.EndPt[1];
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
//        std::cout<<"MWV: coding error. tj "<<tj.ID<<" end "<<end<<" VtxID "<<tj.VtxID[end]<<" != 0\n";
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
          MergeAndStore(tjs, vxtjs[0] - 1, vxtjs[1] - 1, prt);
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
    if(prt) mf::LogVerbatim("TC")<<"Inside FindHammerVertices2";
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      if(tjs.allTraj[it1].CTP != inCTP) continue;
      if(tjs.allTraj[it1].AlgMod[kKilled]) continue;
      if(tjs.allTraj[it1].AlgMod[kHamVx]) continue;
      if(tjs.allTraj[it1].AlgMod[kHamVx2]) continue;
      // Jan 22 Let tj1 be InShower but not tj2
//      if(tjs.allTraj[it1].AlgMod[kInShower]) continue;
      if(tjs.allTraj[it1].AlgMod[kJunkTj]) continue;
      unsigned short numPtsWithCharge1 = NumPtsWithCharge(tjs, tjs.allTraj[it1], false);
      if(numPtsWithCharge1 < 6) continue;
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
          if(tjs.allTraj[it2].AlgMod[kInShower]) continue;
          if(tjs.allTraj[it2].AlgMod[kJunkTj]) continue;
          unsigned short numPtsWithCharge2 = NumPtsWithCharge(tjs, tjs.allTraj[it2], true);
          if(numPtsWithCharge2 < 6) continue;
          // ignore if tj1 is a lot shorter than tj2
          // ignore if ChgRMS isn't known
          // Jan 22. Try this
//          if(tjs.allTraj[it2].ChgRMS == 0) continue;
//          if(numPtsWithCharge1 < 0.2 * numPtsWithCharge2) continue;
          // Find the minimum separation between tj1 and tj2
          float minDOCA = 5;
          float doca = minDOCA;
          unsigned short closePt2 = 0;
          TrajPointTrajDOCA(tjs, tjs.allTraj[it1].Pts[endPt1], tjs.allTraj[it2], closePt2, doca);
          if(doca == minDOCA) continue;
          if(prt) {
            mf::LogVerbatim myprt("TC");
            auto& tj1 = tjs.allTraj[it1];
            auto& tj2 = tjs.allTraj[it2];
            myprt<<" FHV2 CTP"<<tj1.CTP;
            myprt<<" T"<<tj1.ID<<"_"<<end1<<" MCSMom "<<tj1.MCSMom<<" ChgRMS "<<tj1.ChgRMS;
            myprt<<" split T"<<tj2.ID<<"? MCSMom "<<tj2.MCSMom<<" ChgRMS "<<tj2.ChgRMS;
            myprt<<" doca "<<doca<<" tj2.EndPt[0] "<<tj2.EndPt[0]<<" closePt2 "<<closePt2;
            myprt<<" tj2.EndPt[1] "<<tj2.EndPt[1];
          } // prt
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
          // find the MCSMom for both sections of tj2
          float mcsmom = tjs.allTraj[it2].MCSMom;
          float mcsmom1 = MCSMom(tjs, tjs.allTraj[it2], tjs.allTraj[it2].EndPt[0], intPt2);
          float mcsmom2 = MCSMom(tjs, tjs.allTraj[it2], intPt2, tjs.allTraj[it2].EndPt[1]);
          // require that the both MCSMoms be greater than 
          if(prt) mf::LogVerbatim("TC")<<" Check MCSMom after split: mcsmom1 "<<mcsmom1<<" mcsmom2 "<<mcsmom2;
          if(mcsmom1 < mcsmom || mcsmom2 < mcsmom) continue;
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
          if(!SplitTraj(tjs, it2, intPt2, ivx, prt)) {
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
          AttachAnyTrajToVertex(tjs, ivx, prt);
          SetVx2Score(tjs, prt);
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(tjs, it2);
          // and for the new trajectory
          SetPDGCode(tjs, newTjIndex);
          if(prt) mf::LogVerbatim("TC")<<" FHV2: New vtx 2V"<<tjs.vtx[ivx].ID<<" Score "<<tjs.vtx[ivx].Score;
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
      if(tjs.allTraj[it1].AlgMod[kInShower]) continue;
      if(tjs.allTraj[it1].AlgMod[kJunkTj]) continue;
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
          // Let tj1 be InShower but not tj2
//          if(tjs.allTraj[it2].AlgMod[kInShower]) continue;
          if(tjs.allTraj[it2].AlgMod[kJunkTj]) continue;
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
          if(!SplitTraj(tjs, it2, closePt2, ivx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<"FindHammerVertices: Failed to split trajectory";
            tjs.vtx.pop_back();
            continue;
          }
          tjs.allTraj[it1].VtxID[end1] = tjs.vtx[ivx].ID;
          tjs.allTraj[it1].AlgMod[kHamVx] = true;
          tjs.allTraj[it2].AlgMod[kHamVx] = true;
          unsigned short newTjIndex = tjs.allTraj.size() - 1;
          tjs.allTraj[newTjIndex].AlgMod[kHamVx] = true;
          AttachAnyTrajToVertex(tjs, ivx, prt);
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
  void SplitTrajCrossingVertices(TjStuff& tjs, CTP_t inCTP)
  {
    // This is kind of self-explanatory...

    if(!tjs.UseAlg[kSplitTjCVx]) return;

    if(tjs.vtx.empty()) return;
    if(tjs.allTraj.empty()) return;
    
    constexpr float docaCut = 4;

    bool prt = (debug.Plane >= 0 && debug.Tick == 77777);
    if(prt) mf::LogVerbatim("TC")<<"Inside SplitTrajCrossingVertices inCTP "<<inCTP;

    geo::PlaneID planeID = DecodeCTP(inCTP);        

    unsigned short nTraj = tjs.allTraj.size();
    for(unsigned short itj = 0; itj < nTraj; ++itj) {
      // NOTE: Don't use a reference variable because it may get lost if the tj is split
//      auto& tj = tjs.allTraj[itj];
      if(tjs.allTraj[itj].CTP != inCTP) continue;
      // obsolete trajectory
      if(tjs.allTraj[itj].AlgMod[kKilled]) continue;
      if(tjs.allTraj[itj].AlgMod[kSplitTjCVx]) continue;
      // too short
      if(tjs.allTraj[itj].EndPt[1] < 6) continue;
      for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
        auto& vx2 = tjs.vtx[iv];
        // obsolete vertex
        if(vx2.NTraj == 0) continue;
        // trajectory already associated with vertex
        if(tjs.allTraj[itj].VtxID[0] == vx2.ID ||
           tjs.allTraj[itj].VtxID[1] == vx2.ID) continue;
        // not in the cryostat/tpc/plane
        if(tjs.allTraj[itj].CTP != vx2.CTP) continue;
        // poor quality
        if(vx2.Score < tjs.Vertex2DCuts[7]) continue;
        float doca = docaCut;
        // make the cut significantly larger if the vertex is in a dead
        // wire gap to get the first TP that is just outside the gap.
        if(vx2.Stat[kOnDeadWire]) doca = 100;
        unsigned short closePt = 0;
        if(!TrajClosestApproach(tjs.allTraj[itj], vx2.Pos[0], vx2.Pos[1], closePt, doca)) continue;
        if(vx2.Stat[kOnDeadWire]) {
          // special handling for vertices in dead wire regions. Find the IP between
          // the closest point on the Tj and the vertex
          doca = PointTrajDOCA(tjs, vx2.Pos[0], vx2.Pos[1], tjs.allTraj[itj].Pts[closePt]);
        }
        if(doca > docaCut) continue;
        if(prt)  mf::LogVerbatim("TC")<<" doca "<<doca<<" btw T"<<tjs.allTraj[itj].ID<<" and 2V"<<tjs.vtx[iv].ID<<" closePt "<<closePt<<" in plane "<<planeID.Plane;
        // compare the length of the Tjs used to make the vertex with the length of the
        // Tj that we want to split. Don't allow a vertex using very short Tjs to split a long
        // Tj in the 3rd plane
        auto vxtjs = GetVtxTjIDs(tjs, vx2);
        if(vxtjs.empty()) continue;
        unsigned short maxPts = 0;
        // ensure that there is a large angle between a Tj already attached to the vertex and the
        // tj that we want to split. We might be considering a delta-ray here
        float maxdang = 0;
        float tjAng = tjs.allTraj[itj].Pts[closePt].Ang;
        for(auto tjid : vxtjs) {
          auto& vtj = tjs.allTraj[tjid - 1];
          if(vtj.AlgMod[kDeltaRay]) continue;
          unsigned short npwc = NumPtsWithCharge(tjs, vtj, false);
          if(npwc > maxPts) maxPts = npwc;
          unsigned short end = 0;
          if(vtj.VtxID[1] == tjs.vtx[iv].ID) end = 1;
          auto& vtp = vtj.Pts[vtj.EndPt[end]];
          float dang = DeltaAngle(vtp.Ang, tjAng);
          if(dang > maxdang) maxdang = dang; 
        } // tjid
        // skip this operation if any of the Tjs in the split list are > 3 * maxPts
        maxPts *= 3;
        bool skipit = false;
        if(NumPtsWithCharge(tjs, tjs.allTraj[itj], false) > maxPts && maxPts < 100) skipit = true;
        if(!skipit && maxdang < tjs.KinkCuts[0]) skipit = true;
        if(prt) mf::LogVerbatim("TC")<<"  maxPts "<<maxPts<<" vxtjs[0] "<<vxtjs[0]<<" maxdang "<<maxdang<<" skipit? "<<skipit;
        if(skipit) {
          // kill the vertex?
          if(doca < 1) MakeVertexObsolete(tjs, vx2, true);
          continue;
        }
        
        // make some adjustments to closePt
        if(vx2.Stat[kOnDeadWire]) {
          // ensure that the tj will be split at the gap. The closePt point may be
          // on the wrong side of it
          auto& closeTP = tjs.allTraj[itj].Pts[closePt];
          if(tjs.allTraj[itj].StepDir > 0 && closePt > tjs.allTraj[itj].EndPt[0]) {
            if(closeTP.Pos[0] > vx2.Pos[0]) --closePt;
          } else if(tjs.allTraj[itj].StepDir < 0 && closePt < tjs.allTraj[itj].EndPt[1]) {
            if(closeTP.Pos[0] < vx2.Pos[0]) ++closePt;
          }
        } else {
          // improve closePt based on vertex position
          // check if closePt and EndPt[1] are the two sides of vertex
          // take dot product of closePt-vtx and EndPt[1]-vtx
          if ((tjs.allTraj[itj].Pts[closePt].Pos[0]-vx2.Pos[0])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[1]].Pos[0]-vx2.Pos[0]) + (tjs.allTraj[itj].Pts[closePt].Pos[1]-vx2.Pos[1])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[1]].Pos[1]-vx2.Pos[1]) <0 && closePt < tjs.allTraj[itj].EndPt[1] - 1) ++closePt;
          else if ((tjs.allTraj[itj].Pts[closePt].Pos[0]-vx2.Pos[0])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[0]].Pos[0]-vx2.Pos[0]) + (tjs.allTraj[itj].Pts[closePt].Pos[1]-vx2.Pos[1])*(tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[0]].Pos[1]-tjs.vtx[iv].Pos[1]) <0 && closePt > tjs.allTraj[itj].EndPt[0] + 1) --closePt;
        }

        
        if(prt)  {
          mf::LogVerbatim("TC")<<"Good doca "<<doca<<" btw T"<<tjs.allTraj[itj].ID<<" and 2V"<<vx2.ID<<" closePt "<<closePt<<" in plane "<<planeID.Plane<<" CTP "<<tjs.vtx[iv].CTP;
          PrintTrajPoint("STCV", tjs, closePt, 1, tjs.allTraj[itj].Pass, tjs.allTraj[itj].Pts[closePt]);
        }
        // ensure that the closest point is not near an end
        if(closePt < tjs.allTraj[itj].EndPt[0] + 3) continue;
        if(closePt > tjs.allTraj[itj].EndPt[1] - 3) continue;
        if(!SplitTraj(tjs, itj, closePt, iv, prt)) {
          if(prt) mf::LogVerbatim("TC")<<"SplitTrajCrossingVertices: Failed to split trajectory";
          continue;
        }
        tjs.allTraj[itj].AlgMod[kSplitTjCVx] = true;
        unsigned short newTjIndex = tjs.allTraj.size() - 1;
        tjs.allTraj[newTjIndex].AlgMod[kSplitTjCVx] = true;
        // re-fit the vertex position
        FitVertex(tjs, vx2, prt);
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
    
    float thirdPlanedXCut = 2 * tjs.Vertex3DCuts[0];
    
    if(prt) {
      mf::LogVerbatim("TC")<<"Inside Find3DVertices. dX cut "<<tjs.Vertex3DCuts[0]<<" thirdPlanedXCut "<<thirdPlanedXCut;
      PrintAllTraj("F3DV", tjs, debug, USHRT_MAX, tjs.allTraj.size());
    }
    
    // wire spacing in cm
    float wirePitch = tjs.geom->WirePitch(0, 1, 0, tpcid.TPC, tpcid.Cryostat);
    
    size_t vsize = tjs.vtx.size();
    // vector of 2D vertices -> 3D vertices.
    std::vector<short> vPtr(vsize, -1);
    // fill temp vectors of 2D vertex X and X errors
    std::vector<float> vX(vsize, -100);
    
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
    } // ivx
    
    // temp vector of all 2D vertex matches
    std::vector<Vtx3Store> v3temp;
    
    TrajPoint tp;
    float maxScore = 0;
    constexpr float maxSep = 4;
    // i, j, k indicates 3 different wire planes
    // compare vertices in each view
    for(unsigned short ipl = 0; ipl < 2; ++ipl) {
      for(unsigned short ii = 0; ii < vIndex[ipl].size(); ++ii) {
        unsigned short ivx = vIndex[ipl][ii];
        if(vX[ivx] < 0) continue;
        auto& ivx2 = tjs.vtx[ivx];
        unsigned int iWire = std::nearbyint(ivx2.Pos[0]);
        for(unsigned short jpl = ipl + 1; jpl < 3; ++jpl) {
          for(unsigned short jj = 0; jj < vIndex[jpl].size(); ++jj) {
            unsigned short jvx = vIndex[jpl][jj];
            if(vX[jvx] < 0) continue;
            auto& jvx2 = tjs.vtx[jvx];
            unsigned int jWire = std::nearbyint(jvx2.Pos[0]);
            float dX = std::abs(vX[ivx] - vX[jvx]);
            if(dX > tjs.Vertex3DCuts[0]) continue;
            if(prt) {
              mf::LogVerbatim("TC")<<"F3DV: ipl "<<ipl<<" i2V"<<ivx2.ID<<" iX "<<vX[ivx]
              <<" jpl "<<jpl<<" j2V"<<jvx2.ID<<" jvX "<<vX[jvx]<<" W:T "<<(int)jvx2.Pos[0]<<":"<<(int)jvx2.Pos[1]<<" dX "<<dX;
            }
            double y = -1000, z = -1000;
            tjs.geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
            if(y < tjs.YLo || y > tjs.YHi || z < tjs.ZLo || z > tjs.ZHi) continue;
            unsigned short kpl = 3 - ipl - jpl;
            float kX = 0.5 * (vX[ivx] + vX[jvx]);
            int kWire = -1;
            if(tjs.NumPlanes > 2) {
              kWire = (int)(tjs.geom->WireCoordinate(y, z, kpl, tpc, cstat) + 0.5);
              if(kWire < 0 || (unsigned int)kWire > tjs.NumWires[kpl]) continue;
              if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, kpl, kWire))) continue;
              std::array<int, 2> wireWindow;
              std::array<float, 2> timeWindow;
              wireWindow[0] = kWire - maxSep;
              wireWindow[1] = kWire + maxSep;
              float time = tjs.detprop->ConvertXToTicks(kX, kpl, (int)tpc, (int)cstat) * tjs.UnitsPerTick;
              timeWindow[0] = time - maxSep;
              timeWindow[1] = time + maxSep;
              bool hitsNear;
              std::vector<unsigned int> closeHits = FindCloseHits(tjs, wireWindow, timeWindow, kpl, kAllHits, true, hitsNear);
              if(prt) {
                mf::LogVerbatim myprt("TC");
                myprt<<" Hits near "<<kpl<<":"<<kWire<<":"<<(int)(time/tjs.UnitsPerTick)<<" = ";
                for(auto iht : closeHits) myprt<<" "<<PrintHit(tjs.fHits[iht]);
              }
              if(!hitsNear) continue;
            } // 3-plane TPC
/* Don't require a Tj near this position. Look for hits instead
            if(tjs.NumPlanes > 2) {
              kWire = (int)(tjs.geom->WireCoordinate(y, z, kpl, tpc, cstat) + 0.5);
              if(kWire < 0 || (unsigned int)kWire > tjs.NumWires[kpl]) continue;
              if(!tjs.geom->HasWire(geo::WireID(cstat, tpc, kpl, kWire))) continue;
              tp.Pos[0] = kWire;
              // See if there is a wire signal nearby in kpl
              tp.Pos[1] = tjs.detprop->ConvertXToTicks(kX, kpl, tpc, cstat) * tjs.UnitsPerTick;
              tp.CTP = EncodeCTP(cstat, tpc, kpl);
              auto tjlist = FindCloseTjs(tjs, tp, tp, maxSep);
              if(prt) {
                mf::LogVerbatim myprt("TC");
                myprt<<" Tjs near "<<kpl<<":"<<PrintPos(tjs, tp);
                for(auto tjid : tjlist) myprt<<" T"<<tjid;
              }
              if(tjlist.empty()) continue;
            }
*/
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
            
            if(prt) mf::LogVerbatim("TC")<<"F3DV: 2 Plane match i2V"<<tjs.vtx[ivx].ID<<" P:W:T "<<ipl<<":"<<(int)tjs.vtx[ivx].Pos[0]<<":"<<(int)tjs.vtx[ivx].Pos[1]<<" j2V"<<tjs.vtx[jvx].ID<<" P:W:T "<<jpl<<":"<<(int)tjs.vtx[jvx].Pos[0]<<":"<<(int)tjs.vtx[jvx].Pos[1]<<" dX "<<dX;
            
            if(tjs.NumPlanes == 2) continue;
            
            // look for a 3 plane match
            for(unsigned short kk = 0; kk < vIndex[kpl].size(); ++kk) {
              unsigned short kvx = vIndex[kpl][kk];
              if(vX[kvx] < 0) continue;
              float dX = std::abs(vX[kvx] - v3d.X);
              // Wire difference error
              float dW = wirePitch * std::abs(tjs.vtx[kvx].Pos[0] - kWire);
              if(prt) mf::LogVerbatim("TC")<<" k2V"<<kvx+1<<" dX "<<dX<<" dW "<<dW;
              if(dX > thirdPlanedXCut) continue;
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
      vx3.ID = tjs.vtx3.size() + 1;
      if(prt) mf::LogVerbatim("TC")<<" 3V"<<vx3.ID<<"  2V"<<vx3.Vx2ID[0]<<" 2V"<<vx3.Vx2ID[1]<<" 2V"<<vx3.Vx2ID[2]
        <<" wire "<<vx3.Wire;
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
    for(unsigned short iv3 = 0; iv3 < tjs.vtx3.size(); ++iv3) {
      auto& vx3 = tjs.vtx3[iv3];
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      // put the TjIDs into a vector for matching
      float score = 0;
      auto v3TjIDs = GetVtxTjIDs(tjs, vx3, score);
      if(v3TjIDs.empty()) continue;
      if(score < tjs.Vertex2DCuts[7]) continue;
/*
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"M3DVTj 3V"<<vx3.ID<<" score "<<score<<" TjIDs";
        for(auto& tjID : v3TjIDs) myprt<<" T"<<std::to_string(tjID);
      }
*/
      SortEntry se;
      se.index = iv3;
      se.val = score;
      sortVec.push_back(se);
    } // vx3
    if(sortVec.empty()) return;
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasing);
    
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      auto& vx3 = tjs.vtx3[sortVec[ii].index];
      float score = 0;
      std::vector<int> v3TjIDs = GetVtxTjIDs(tjs, vx3, score);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"M3DVTj 3V"<<vx3.ID<<" score "<<score<<" Tjs";
        for(auto& tjID : v3TjIDs) myprt<<" T"<<tjID;
      }
      for(unsigned int ims = 0; ims < tjs.matchVec.size(); ++ims) {
        auto& ms = tjs.matchVec[ims];
        if(ms.Count == 0) continue;
        bool skipit = false;
        for(auto tjid : ms.TjIDs) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kMat3D]) skipit = true;
          if(tj.AlgMod[kKilled]) skipit = true;
        }
        if(skipit) continue;
        std::vector<int> shared = SetIntersection(ms.TjIDs, v3TjIDs);
        if(shared.size() < 2) continue;
        if(prt) mf::LogVerbatim("TC")<<" ims "<<ims<<" shared size "<<shared.size();
        PFPStruct pfp = CreatePFP(tjs, tpcid);
        pfp.TjIDs = ms.TjIDs;
        pfp.Vx3ID[0] = vx3.ID;
        // note that the ms position is the average position of all 3D matched Tp3s at this point.
        // It is not the start position. This will be determined in DefinePFP.
        pfp.XYZ[0] = ms.Pos;
        pfp.Dir[0] = ms.Dir;
        pfp.MatchVecIndex = ims;
        if(prt) mf::LogVerbatim("TC")<<"M3DVTj: pfp P"<<pfp.ID<<" 3V"<<vx3.ID<<" ims "<<ims;
        // Set the PDGCode so DefinePFP can ignore incompatible matches
        pfp.PDGCode = PDGCodeVote(tjs, pfp.TjIDs, prt);
        // make a copy of the TjIDs to test for changes
        auto saveTjIDs = pfp.TjIDs;
        if(!DefinePFP("M3DVTj1", tjs, pfp, prt)) continue;
        // separation distance (cm) for kink detection.
        double sep = 1;
        bool didSplit = Split3DKink(tjs, pfp, sep, prt);
        if(prt) PrintPFP("M3D", tjs, pfp, true);
        if(!didSplit && shared.size() != ms.TjIDs.size()) {
          // Try to repair the PFParticle by merging the Tj that was in the match list but
          // wasn't attached to the vertex. Hopefully there aren't more than one...
          auto tjNotInVx = SetDifference(ms.TjIDs, shared);
        }
        AnalyzePFP(tjs, pfp, prt);
        if(!StorePFP(tjs, pfp)) continue;
        if(pfp.TjIDs != saveTjIDs) {
          // v3TjIDs is wrong if Tjs were merged so re-make it.
          auto tmp = GetVtxTjIDs(tjs, vx3, score);
          v3TjIDs.clear();
          // then re-build it
          for(auto tjid : tmp) {
            auto& tj = tjs.allTraj[tjid - 1];
            if(!tj.AlgMod[kMat3D]) v3TjIDs.push_back(tjid);
          } // tjid
          if(v3TjIDs.size() < 2) break;
        } // pfp.TjIDs != saveTjIDs
        ms.Count = 0;
        // clobber MatchStructs that use the Tjs in this pfp
        for(auto& allms : tjs.matchVec) {
          auto mfpfp = SetIntersection(allms.TjIDs, pfp.TjIDs);
          if(!mfpfp.empty()) allms.Count = 0;
        } // allms
        auto leftover = SetDifference(v3TjIDs, pfp.TjIDs);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<" leftover";
          for(auto tjid : leftover) myprt<<" T"<<tjid;
        }
        if(leftover.size() < 2) break;
        // keep looking using the leftovers
        v3TjIDs = leftover;
      } // ims
      if(v3TjIDs.size() > 1 && v3TjIDs.size() <= tjs.NumPlanes) {
        // a last-ditch attempt
        PFPStruct pfp = CreatePFP(tjs, tpcid);
        // put the available leftovers in
        for(auto& tjid : v3TjIDs) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kMat3D]) continue;
          if(tj.AlgMod[kKilled]) continue;
          pfp.TjIDs.push_back(tjid);
        } // tjid
        if(pfp.TjIDs.size() < 2) continue;
        pfp.Vx3ID[0] = vx3.ID;
        if(!DefinePFP("M3DVTj2", tjs, pfp, prt)) continue;
        Split3DKink(tjs, pfp, 1, prt);
        AnalyzePFP(tjs, pfp, prt);
        if(prt) PrintPFP("left", tjs, pfp, true);
        if(!StorePFP(tjs, pfp)) continue;
      }
    } // ims
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
    if(vx3ID > tjs.vtx3.size()) return false;
    if(pfp.ID > tjs.pfps.size()) return false;
    if(pfp.PDGCode == 22) return false;
    if(end > 1) return false;
    
    auto& vx3 = tjs.vtx3[vx3ID - 1];
    
    pfp.Vx3ID[end] = vx3.ID;
    
    // We are done if this a PFP-only vertex
    if(vx3.Wire == -2) return true;
    
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
//          if(prt) std::cout<<" tj "<<tj.ID<<" has no 2D vertex. Look for one vertex near "<<tj.CTP<<":"<<PrintPos(tjs, pos)<<" Events processed "<<tjs.EventsProcessed<<"\n";
        } else {
          // Existing 2D vertex matched to the 3D vertex
//          if(prt) std::cout<<" tj "<<tj.ID<<" has no 2D vertex in CTP "<<tj.CTP<<" but vx3 is matched to 2D vertex"<<vx3.Vx2ID[plane]<<". Attach it? Events processed "<<tjs.EventsProcessed<<"\n";
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
    float closestApproach;
    // ignore bad separation between the closest tj end and the vertex
    if(tjShort) {
      if(vtxTjSep2 > maxSepCutShort2) return false;
      closestApproach = tjs.Vertex2DCuts[1];
    } else {
      closestApproach = tjs.Vertex2DCuts[2];
      if(vtxTjSep2 > maxSepCutLong2) return false;
    }
    
    // Calculate the pull on the vertex
    TrajPoint& tp = tj.Pts[tj.EndPt[end]];
    float tpVxPull = TrajPointVertexPull(tjs, tp, vx);
    bool signalBetween = SignalBetween(tjs, tp, vx.Pos[0], 0.8, prt);
    
    // See if the vertex position is close to an end
    unsigned short closePt;
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

    float pullCut = tjs.Vertex2DCuts[3];
    // Dec 21, 2017 Loosen up the pull cut for short close Tjs. These are likely to
    // be poorly reconstructed. It is better to have them associated with the vertex
    // than not.
    if(tjShort) pullCut = 10;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"ATTV: 2V"<<vx.ID;
      myprt<<" NTraj "<<vx.NTraj;
      myprt<<" oldTJs";
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.AlgMod[kKilled]) continue;
        if(tj.CTP != vx.CTP) continue;
        if(tj.VtxID[0] == vx.ID) myprt<<" "<<tj.ID<<"_0";
        if(tj.VtxID[1] == vx.ID) myprt<<" "<<tj.ID<<"_1";
      }
      myprt<<" +tjID "<<tj.ID<<"_"<<end<<" vtxTjSep "<<sqrt(vtxTjSep2)<<" tpVxPull "<<tpVxPull<<" pullCut "<<pullCut<<" dpt "<<dpt;
    }
//    if(tpVxPull > tjs.Vertex2DCuts[3]) return false;
    if(tpVxPull > pullCut) return false;
    if(dpt > 2) return true;
    
    // remove the fixed position flag if there are more than 2 tjs
    bool fixedBit = vx.Stat[kFixed];
    if(fixedBit && vx.NTraj < 2) vx.Stat[kFixed] = false;
    
    // don't allow a short Tj with a large pull to bias the fit
    if(tjShort && tpVxPull > tjs.Vertex2DCuts[3]) tj.AlgMod[kNoFitToVx] = true;

    // Passed all the cuts. Attach it to the vertex and try a fit
    tj.VtxID[end] = vx.ID;
    // flag as a photon Tj so it isn't included in the fit
    tj.AlgMod[kPhoton] = !signalBetween;
    // make a copy of the vertex and fit it
    auto vxTmp = vx;
    if(FitVertex(tjs, vxTmp, prt)) {
      SetVx2Score(tjs, vxTmp, prt);
      if(prt) mf::LogVerbatim("TC")<<" Success";
      vx = vxTmp;
      return true;
    }
    
    // test this again
    tj.AlgMod[kNoFitToVx] = true;
    if(prt) mf::LogVerbatim("TC")<<" Poor fit. Keep Tj "<<tj.ID<<" with kNoFitToVx";
    return true;
/*
    // fit failed so remove the tj -> vx assignment if it is long and
    // set noFitToVtx if it is short
    if(tjShort) {
      tj.AlgMod[kNoFitToVx] = true;
      if(prt) mf::LogVerbatim("TC")<<" Poor fit. Keep short Tj "<<tj.ID<<" with kNoFitToVx";
      return true;
    } else {
      tj.VtxID[end] = 0;
      // restore the fixed flag
      vx.Stat[kFixed] = fixedBit;
      if(prt) mf::LogVerbatim("TC")<<" Poor fit. Removed Tj "<<tj.ID;
      return false;
    }
*/
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
    // add the TP position error^2
    vxErr2 += tp.HitPosErr2;
    
    // close together so ignore the TP projection error and return
    // the pull using the vertex error and TP position error
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
      mf::LogVerbatim("TC")<<"StoreVertex: vertex "<<vx.ID<<" Topo "<<vx.Topo<<" has inconsistent CTP code "<<vx.CTP<<" with one or more Tjs\n";
      for(auto& tj : tjs.allTraj) {
        if(tj.AlgMod[kKilled]) continue;
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
      if(tj.AlgMod[kPhoton]) continue;
      if(tj.AlgMod[kNoFitToVx]) continue;
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
    double cnt = 0;
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
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; 
        sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt; ++cnt;
        // tweak the itj angle +
        tmp = vxTp[itj];
        tmp.Ang += tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(tmp, vxTp[jtj], p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // adjust the weight for 4 points at +/1 1 sigma = 0.607 / 4
        wgt = 0.152;
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; 
        sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt; ++cnt;
        // tweak the itj angle -
        tmp = vxTp[itj];
        tmp.Ang -= 2 * tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(tmp, vxTp[jtj], p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; 
        sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt; ++cnt;
        // Repeat this process with jtj
        // tweak the jtj angle +
        tmp = vxTp[jtj];
        tmp.Ang += tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(vxTp[itj], tmp, p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; 
        sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt; ++cnt;
        // tweak the itj angle -
        tmp = vxTp[itj];
        tmp.Ang -= 2 * tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(vxTp[itj], tmp, p0, p1);
        intTp.Pos[0] = p0; intTp.Pos[1] = p1;
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; 
        sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt; ++cnt;
      } // jtj
    } // itj
    
    if(sumw == 0) return false;
    
    double vxP0 = sum0 / sumw;
    double vxP1 = sum1 / sumw;
    double vxP0rms = sqrt((sum02 - sumw * vxP0 * vxP0) / sumw);
    double vxP1rms = sqrt((sum12 - sumw * vxP1 * vxP1) / sumw);
    double rootN = sqrt(cnt);
    vxP0rms /= rootN;
    vxP1rms /= rootN;
    // don't let the errors get too small
    if(vxP0rms < 0.5) vxP0rms = 0.5;
    if(vxP1rms < 0.5) vxP1rms = 0.5;
    
    if(prt) mf::LogVerbatim("TC")<<"FitVertex 2V"<<vx.ID<<" CTP "<<vx.CTP<<" NTraj "<<vx.NTraj<<" in "<<std::fixed<<std::setprecision(1)<<vx.Pos[0]<<":"<<vx.Pos[1]/tjs.UnitsPerTick<<" out wire "<<vxP0<<" +/- "<<vxP0rms<<" ticks "<<vxP1/tjs.UnitsPerTick<<"+/-"<<vxP1rms/tjs.UnitsPerTick;
    
    // apply Vertex2DCuts if this isn't a neutral vertex (which is expected to have very large
    // errors)
    if(vx.Topo != 11) {
      float inflate = 1;
      if(vx.Stat[kOnDeadWire]) inflate = 1.5;
      if(vxP0rms > inflate * tjs.Vertex2DCuts[4] || vxP1rms > inflate * tjs.Vertex2DCuts[4]) {
        if(prt) mf::LogVerbatim("TC")<<" fit failed. Max allowed position error "<<inflate * tjs.Vertex2DCuts[4];
        return false;
      }
    } // not a neutral vertex
    
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
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: 2V"<<vx2.ID<<" thinks it is matched to 3V"<<vx3.ID<<" but vx3 is obsolete";
        return false;
      }
      if(vx3.Vx2ID[plane] != vx2.ID) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: 2V"<<vx2.ID<<" thinks it is matched to 3V"<<vx3.ID<<" but vx3 says no!";
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
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: 3V"<<vx3.ID<<" thinks it is matched to 2V"<<vx2.ID<<" but vx2 says no!";
        return false;
      }
    } // vx3
    
    // check the Tj -> 2D associations
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == 0) continue;
        if(tj.VtxID[end] > tjs.vtx.size()) {
          mf::LogVerbatim("TC")<<"ChkVtxAssociations: T"<<tj.ID<<" thinks it is matched to 2V"<<tj.VtxID[end]<<" on end "<<end<<" but no vertex exists. Recovering";
          tj.VtxID[end] = 0;
          return false;
        }
        unsigned short ivx = tj.VtxID[end] - 1;
        auto& vx2 = tjs.vtx[ivx];
        if(vx2.ID == 0) {
          mf::LogVerbatim("TC")<<"ChkVtxAssociations: T"<<tj.ID<<" thinks it is matched to 2V"<<tj.VtxID[end]<<" on end "<<end<<" but the vertex is killed. Fixing the Tj";
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
    } // vx
    // Score the 3D vertices
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
       SetVx3Score(tjs, vx3, prt);
    } // vx3
  } // ScoreVertices
  
  //////////////////////////////////////////
  void KillPoorVertices(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // kill 2D vertices that have low score and are not attached to a high-score 3D vertex
    if(tjs.vtx.empty()) return;
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    for(auto& vx : tjs.vtx) {
      if(vx.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(vx.CTP);
      if(planeID.Cryostat != cstat) continue;
      if(planeID.TPC != tpc) continue;
      if(vx.Score > tjs.Vertex2DCuts[7]) continue;
      if(vx.Vx3ID > 0) {
        auto& vx3 = tjs.vtx3[vx.Vx3ID - 1];
        if(vx3.Primary) continue;
        if(tjs.vtx3[vx.Vx3ID - 1].Score >= tjs.Vertex2DCuts[7]) continue;
      }
      MakeVertexObsolete(tjs, vx, false);
    } // vx
    
  } // KillPoorVertices
  
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
      while(true) {
        // tag Tjs and make a list of attached vertices whose high-score
        // bit needs to be set
        vxlist.clear();
        for(auto tjid : tjlist) {
          auto& tj = tjs.allTraj[tjid - 1];
          tj.AlgMod[kTjHiVx3Score] = true;
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] == 0) continue;
            auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
            if(vx2.Stat[kHiVx3Score]) continue;
            vx2.Stat[kHiVx3Score] = true;
            vxlist.push_back(vx2.ID);
          } // end
        } // tjid

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
    if(vx2.ID == 0) return;
    
    // Don't score vertices from CheckTrajBeginChg, MakeJunkVertices or Neutral vertices. Set to the minimum
    if(vx2.Topo == 8 || vx2.Topo == 9 || vx2.Topo == 11) {
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
    
    auto vtxTjIDs = GetVtxTjIDs(tjs, vx2);
    if(vtxTjIDs.empty()) return;

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
    
    vx2.TjChgFrac = ChgFracNearPos(tjs, vx2.Pos, vtxTjIDs);
    float cfScore = tjs.VertexScoreWeights[2] * vx2.TjChgFrac;
    
    // Define a weight for each Tj
    std::vector<int> tjids;
    std::vector<float> tjwts;
    for(auto tjid : vtxTjIDs) {
      Trajectory& tj = tjs.allTraj[tjid - 1];
      // Feb 22 Ignore short Tjs and junk tjs
      if(tj.AlgMod[kJunkTj]) continue;
      unsigned short lenth = tj.EndPt[1] - tj.EndPt[0] + 1;
      if(lenth < 3) continue;
      float wght = (float)tj.MCSMom / momBin;
      if(wght > 10) wght = 10;
      // weight by tagged muon
      if(tj.PDGCode == 13) wght *= 2;
      // weight by charge rms
      if(tj.ChgRMS < maxChgRMS) ++wght;
      // Shower Tj
      if(tj.AlgMod[kShowerTj]) ++wght;
      // InShower
      if(tj.AlgMod[kInShower]) --wght;
      tjids.push_back(tjid);
      tjwts.push_back(wght);
    } // tjid
    
    if(tjids.empty()) return;
    
    float tjScore = 0;
    float sum = 0;
    float cnt = 0;
    for(unsigned short it1 = 0; it1 < tjids.size() - 1; ++it1) {
      Trajectory& tj1 = tjs.allTraj[tjids[it1] - 1];
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
      for(unsigned short it2 = it1 + 1; it2 < tjids.size(); ++it2) {
        Trajectory& tj2 = tjs.allTraj[tjids[it2] - 1];
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
    vx2.Score = vpeScore + m3DScore + cfScore + tjScore;
    if(prt) {
      // last call after vertices have been matched to the truth. Use to optimize VertexScoreWeights using
      // an ntuple
      mf::LogVerbatim myprt("TC");
      myprt<<" SVx2W 2V"<<vx2.ID;
      myprt<<" m3Dcnt "<<m3Dcnt;
      myprt<<" PosErr "<<std::fixed<<std::setprecision(2)<<(vx2.PosErr[0] + vx2.PosErr[1]);
      myprt<<" TjChgFrac "<<std::fixed<<std::setprecision(3)<<vx2.TjChgFrac;
      myprt<<" sum "<<std::fixed<<std::setprecision(1)<<sum;
      myprt<<" cnt "<<(int)cnt;
      myprt<<" Score "<<vx2.Score;
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
    std::array<short, 10> cnts;
    cnts.fill(0);
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
        if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: new vtx 2V"<<aVtx.ID<<" points to 3V"<<vx3.ID;
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
      if(prt) mf::LogVerbatim("TC")<<"CI3DV 3V"<<vx3.ID<<" Pos "<<mPlane<<":"<<PrintPos(tjs, vtp.Pos);
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
        if(prt) mf::LogVerbatim("TC")<<"CI3DV 3V"<<vx3.ID<<" candidate itj ID "<<tj.ID<<" vtx pos "<<PrintPos(tjs, vtp.Pos)<<" doca "<<doca;
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
        // Don't use a reference variable since it may get scrambled after SplitTraj
//        auto& tj = tjs.allTraj[itj];
        unsigned short closePt = tjPts[ii];
        // determine which end is the closest
        unsigned short end = 1;
        // closest to the beginning?
        if(fabs(closePt - tjs.allTraj[itj].EndPt[0]) < fabs(closePt - tjs.allTraj[itj].EndPt[1])) end = 0;
        short dpt = fabs(closePt - tjs.allTraj[itj].EndPt[end]);
        if(dpt < 3) {
          // close to an end
          if(tjs.allTraj[itj].VtxID[end] > 0) {
            if(prt) mf::LogVerbatim("TC")<<" T"<<tjs.allTraj[itj].ID<<" has a vertex "<<tjs.allTraj[itj].VtxID[end]<<" at end "<<end<<". Skip it";
            continue;
          }
          tjs.allTraj[itj].VtxID[end] = tjs.vtx[newVtxIndx].ID;
          ++newVtx.NTraj;
          if(prt) mf::LogVerbatim("TC")<<" attach Traj T"<<tjs.allTraj[itj].ID<<" at end "<<end;
          tjs.allTraj[itj].AlgMod[kComp3DVx] = true;
          vpos = tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[end]].Pos;
        } else {
          // closePt is not near an end, so split the trajectory
          if(SplitTraj(tjs, itj, closePt, newVtxIndx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<" SplitTraj success 2V"<<tjs.vtx[newVtxIndx].ID<<" at closePt "<<closePt;
            // successfully split the Tj
            newVtx.NTraj += 2;
          } else {
            // split failed. Give up
            if(prt) mf::LogVerbatim("TC")<<" SplitTraj failed";
            newVtx.NTraj = 0;
            break;
          }
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(tjs, itj);
          // and for the new trajectory
          SetPDGCode(tjs, tjs.allTraj.size()-1);
        } // closePt is not near an end, so split the trajectory
        tjs.allTraj[itj].AlgMod[kComp3DVx] = true;
        unsigned short newtj = tjs.allTraj.size() - 1;
        tjs.allTraj[newtj].AlgMod[kComp3DVx] = true;
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
          myprt<<" Success: new 2V"<<newVtx.ID<<" at "<<(int)newVtx.Pos[0]<<":"<<(int)newVtx.Pos[1]/tjs.UnitsPerTick;
          myprt<<" points to 3V"<<vx3.ID;
          myprt<<" TjIDs:";
          for(auto& tjID : tjIDs) myprt<<" T"<<std::to_string(tjID);
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
        tj.AlgMod[kPhoton] = false;
        // clear the kEnvOverlap bits on the TPs
        for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
          if(end == 0) { 
            unsigned short ipt = tj.EndPt[0] + ii;
            auto& tp = tj.Pts[ipt];
            if(!tp.Environment[kEnvOverlap]) break;
            if(ipt == tj.EndPt[1]) break;
          } else {
            unsigned short ipt = tj.EndPt[1] - ii;
            auto& tp = tj.Pts[ipt];
            if(!tp.Environment[kEnvOverlap]) break;
            if(ipt == tj.EndPt[0]) break;
          }
        } // ii
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
    vx3.Wire = vx2.Pos[0];
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
    for(auto& pfp : tjs.pfps) {
      for(unsigned short end = 0; end < 2; ++end) if(pfp.Vx3ID[end] == vx3.ID) pfp.Vx3ID[end] = 0;
    } // pfp
    vx3.ID = 0;
    return true;
    
  } // MakeVertexObsolete

  ////////////////////////////////////////////////
  bool MakeVertexObsolete(TjStuff& tjs, Vtx3Store& vx3)
  {
    // Deletes a 3D vertex and 2D vertices in all planes
    // The 2D and 3D vertices are NOT killed if forceKill is false and the 3D vertex
    // has a high score
    if(vx3.ID == 0) return true;
    if(vx3.ID > tjs.vtx3.size()) return false;
    
    // set the score to 0
    vx3.Score = 0;
    
    for(auto vx2id : vx3.Vx2ID) {
      if(vx2id == 0 || vx2id > tjs.vtx.size()) continue;
      auto& vx2 = tjs.vtx[vx2id - 1];
      MakeVertexObsolete(tjs, vx2, true);
    }
    return true;
  } // MakeVertexObsolete
  
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
  std::vector<unsigned short> GetPFPVertices(const TjStuff& tjs, const PFPStruct& pfp)
  {
    // returns a list of 3D vertices that are attached to Tjs in this pfp. No check is
    // made of the actual vertex attachment of the pfp.
    std::vector<unsigned short> tmp;
    if(pfp.TjIDs.empty()) return tmp;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == 0) continue;
        auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
        if(vx2.Vx3ID == 0) continue;
        if(std::find(tmp.begin(), tmp.end(), vx2.Vx3ID) != tmp.end()) continue;
        tmp.push_back(vx2.Vx3ID);
      } // end
    } // tjid
    return tmp;
  } // GetPFPVertices

  //////////////////////////////////////////
  void PosInPlane(const TjStuff& tjs, const Vtx3Store& vx3, unsigned short plane, Point2_t& pos)
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
  
} // namespace
