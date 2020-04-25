#include "larreco/RecoAlg/TCAlg/TCVertex.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <stdlib.h>
#include <string>
#include <vector>

namespace tca {

  struct SortEntry{
    unsigned int index;
    float val;
  };

  bool valDecreasing (SortEntry c1, SortEntry c2) { return (c1.val > c2.val);}
  bool valIncreasing (SortEntry c1, SortEntry c2) { return (c1.val < c2.val);}

  //////////////////////////////////////////
  void MakeJunkVertices(TCSlice& slc, const CTP_t& inCTP)
  {
    // Vertices between poorly reconstructed tjs (especially junk slc) and normal
    // tjs can fail because the junk tj trajectory parameters are inaccurate. This function
    // uses proximity and not pointing to make junk vertices
    // Don't use this if standard vertex reconstruction is disabled
    if(tcc.vtx2DCuts[0] <= 0) return;
    if(!tcc.useAlg[kJunkVx]) return;
    if(slc.tjs.size() < 2) return;

    // Look for tjs that are within maxSep of the end of a Tj
    constexpr float maxSep = 4;

    geo::PlaneID planeID = DecodeCTP(inCTP);
    bool prt = (tcc.dbgVxJunk && tcc.dbgSlc);
    if(prt) {
      mf::LogVerbatim("TC")<<"MakeJunkVertices: prt set for plane "<<planeID.Plane<<" maxSep btw tjs "<<maxSep;
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
    junkVx.PosErr = {{2.0, 2.0}};
    // define a minimal score so it won't get clobbered
    junkVx.Score = tcc.vtx2DCuts[7] + 0.1;

    // look at both ends of long tjs
    for(unsigned short it1 = 0; it1 < slc.tjs.size() - 1; ++it1) {
      auto& tj1 = slc.tjs[it1];
      if(tj1.AlgMod[kKilled] || tj1.AlgMod[kHaloTj]) continue;
      if(tj1.SSID > 0 || tj1.AlgMod[kShowerLike]) continue;
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kJunkTj]) continue;
      if(TrajLength(tj1) < 10) continue;
      if(tj1.MCSMom < 100) continue;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // existing vertex?
        if(tj1.VtxID[end1] > 0) continue;
        auto& tp1 = tj1.Pts[tj1.EndPt[end1]];
        // get a list of tjs in this vicinity
        auto tjlist = FindCloseTjs(slc, tp1, tp1, maxSep);
        if(tjlist.empty()) continue;
        // set to an invalid ID
        junkVx.ID = USHRT_MAX;
        for(auto tj2id : tjlist) {
          auto& tj2 = slc.tjs[tj2id - 1];
          if(tj2.CTP != inCTP) continue;
          if(tj2id == tj1.ID) continue;
          if(tj2.SSID > 0 || tj2.AlgMod[kShowerLike]) continue;
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
          bool signalBetween = SignalBetween(slc, tp1, tp2, 0.8);
          if(!signalBetween) continue;
          if(junkVx.ID == USHRT_MAX) {
            // define the new vertex
            junkVx.ID = slc.vtxs.size() + 1;
            junkVx.Pos = tp1.Pos;
          } // new vertex
          tj2.VtxID[closeEnd] = junkVx.ID;
          tj1.VtxID[end1] = junkVx.ID;
        } // tjid
        if(junkVx.ID == USHRT_MAX) continue;
        if(!StoreVertex(slc, junkVx)) {
          mf::LogVerbatim("TC")<<"MJV: StoreVertex failed";
          for(auto& tj : slc.tjs) {
            if(tj.AlgMod[kKilled]) continue;
            if(tj.VtxID[0] == junkVx.ID) tj.VtxID[0] = 0;
            if(tj.VtxID[1] == junkVx.ID) tj.VtxID[1] = 0;
          } // tj
          continue;
        } // StoreVertex failed
        if(prt) {
          mf::LogVerbatim("TC")<<" New junk 2V"<<junkVx.ID<<" at "<<std::fixed<<std::setprecision(1)<<junkVx.Pos[0]<<":"<<junkVx.Pos[1]/tcc.unitsPerTick;
        } // prt
        junkVx.ID = USHRT_MAX;
      } // end1
    } // it1

  } // MakeJunkVertices

  //////////////////////////////////////////
  void Find2DVertices(TCSlice& slc, const CTP_t& inCTP, unsigned short pass)
  {
    // Find 2D vertices between pairs of tjs that have a same-end topology. Using an example
    // where StepDir = 1 (end 0 is at small wire number) vertices will be found with Topo = 0
    // with a vertex US of the ends (<) or Topo = 2 with a vertex DS of the ends (>). This is reversed
    // if StepDir = -1. Vertices with Topo = 1 (/\) and (\/) are found in EndMerge.

    // tcc.vtx2DCuts fcl input usage
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
    // 10 Require charge on wires between a vtx and the start of the tjs in induction planes? (1 = yes)

    if(tcc.vtx2DCuts[0] <= 0) return;
    if(slc.tjs.size() < 2) return;

    bool firstPassCuts = (pass == 0);

    geo::PlaneID planeID = DecodeCTP(inCTP);

    // require charge between the vertex and the tj start points?
    bool requireVtxTjChg = true;
    if(tcc.vtx2DCuts[10] == 0 && int(planeID.Plane) < slc.nPlanes - 1) requireVtxTjChg = false;

    bool prt = (tcc.dbg2V && tcc.dbgSlc && inCTP == debug.CTP);
    if(prt) {
      mf::LogVerbatim("TC")<<"prt set for CTP "<<inCTP<<" in Find2DVertices. firstPassCuts? "<<firstPassCuts<<" requireVtxTjChg "<<requireVtxTjChg;
      PrintAllTraj("F2DVi", slc, USHRT_MAX, slc.tjs.size());
    }

    unsigned short maxShortTjLen = tcc.vtx2DCuts[0];
    for(unsigned short it1 = 0; it1 < slc.tjs.size() - 1; ++it1) {
      auto& tj1 = slc.tjs[it1];
      if(tj1.AlgMod[kKilled] || tj1.AlgMod[kHaloTj]) continue;
      if(tj1.SSID > 0 || tj1.AlgMod[kShowerLike]) continue;
      if(tj1.CTP != inCTP) continue;
      bool tj1Short = (TrajLength(tj1) < maxShortTjLen);
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(tj1.VtxID[end1] > 0) continue;
        // wrong end of a high energy electron?
        if(tj1.PDGCode == 111 && end1 != tj1.StartEnd) continue;
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
          if(tj1.Pts[endPt1].Chg == 0) endPt1 = NearestPtWithChg(slc, tj1, endPt1);
        } // few points fit at end1
        TrajPoint tp1 = tj1.Pts[endPt1];
        MoveTPToWire(tp1, wire1);
        // re-purpose endPt1 to reference the end point. This will be used the find the point on
        // tj1 that is closest to the vertex position
        endPt1 = tj1.EndPt[end1];
        short oendPt1 = tj1.EndPt[1-end1];
        // reference to the other end of tj1
        auto& otp1 = tj1.Pts[oendPt1];
        for(unsigned short it2 = it1 + 1; it2 < slc.tjs.size(); ++it2) {
          auto& tj2 = slc.tjs[it2];
          if(tj2.AlgMod[kKilled] || tj2.AlgMod[kHaloTj]) continue;
          if(tj2.SSID > 0 || tj2.AlgMod[kShowerLike]) continue;
          if(tj2.CTP != inCTP) continue;
          if(tj1.VtxID[end1] > 0) continue;
          if(tj1.MCSMom < tcc.vtx2DCuts[5] && tj2.MCSMom < tcc.vtx2DCuts[5]) continue;
          bool tj2Short = (TrajLength(tj2) < maxShortTjLen);
          // find the end that is closer to tp1
          unsigned short end2 = 0;
          if(PosSep2(tj2.Pts[tj2.EndPt[1]].Pos, tp1.Pos) < PosSep2(tj2.Pts[tj2.EndPt[0]].Pos, tp1.Pos)) end2 = 1;
          if(tj2.VtxID[end2] > 0) continue;
          // wrong end of a high energy electron?
          if(tj2.PDGCode == 111 && end2 != tj2.StartEnd) continue;
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
            if(tj2.Pts[endPt2].Chg == 0) endPt2 = NearestPtWithChg(slc, tj2, endPt2);
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
          // make sure this is inside the TPC.
          if(wint < 0 || wint > tcc.maxPos0[planeID.Plane] - 3) continue;
          if(tint < 0 || tint > tcc.maxPos1[planeID.Plane]) continue;
          // Next cut on separation between the TPs and the intersection point
          if(tj1Short || tj2Short) { sepCut = tcc.vtx2DCuts[1]; } else { sepCut = tcc.vtx2DCuts[2]; }
          // NewVtxCuts: require close separation on the first pass
          if(firstPassCuts) sepCut = tcc.vtx2DCuts[1];
          Point2_t vPos {{wint, tint}};
          float vt1Sep = PosSep(vPos, tp1.Pos);
          float vt2Sep = PosSep(vPos, tp2.Pos);
          float dwc1 = DeadWireCount(slc, wint, tp1.Pos[0], tp1.CTP);
          float dwc2 = DeadWireCount(slc, wint, tp2.Pos[0], tp1.CTP);
          vt1Sep -= dwc1;
          vt2Sep -= dwc2;
          bool vtxOnDeadWire = (DeadWireCount(slc, wint, wint, tp1.CTP) == 1);
          if(prt && vt1Sep < 200 && vt2Sep < 200) {
            mf::LogVerbatim myprt("TC");
            myprt<<"F2DV candidate T"<<tj1.ID<<"_"<<end1<<"-T"<<tj2.ID<<"_"<<end2;
            myprt<<" vtx pos "<<(int)wint<<":"<<(int)(tint/tcc.unitsPerTick)<<" tp1 "<<PrintPos(slc, tp1)<<" tp2 "<<PrintPos(slc, tp2);
            myprt<<" dwc1 "<<dwc1<<" dwc2 "<<dwc2<<" on dead wire? "<<vtxOnDeadWire;
            myprt<<" vt1Sep "<<vt1Sep<<" vt2Sep "<<vt2Sep<<" sepCut "<<sepCut;
          }
          if(vt1Sep > sepCut || vt2Sep > sepCut) continue;
          // make sure that the other end isn't closer
          if(PosSep(vPos, slc.tjs[it1].Pts[oendPt1].Pos) < vt1Sep) {
            if(prt) mf::LogVerbatim("TC")<<" tj1 other end "<<PrintPos(slc, tj1.Pts[oendPt1])<<" is closer to the vertex";
            continue;
          }
          if(PosSep(vPos, slc.tjs[it2].Pts[oendPt2].Pos) < vt2Sep) {
            if(prt) mf::LogVerbatim("TC")<<" tj2 other end "<<PrintPos(slc, tj2.Pts[oendPt2])<<" is closer to the vertex";
            continue;
          }
          // Ensure that the vertex position is close to the end of each Tj
          unsigned short closePt1;
          float doca1 = sepCut;
          if(!TrajClosestApproach(tj1, wint, tint, closePt1, doca1)) continue;
          // dpt1 (and dpt2) will be 0 if the vertex is at the end
          short stepDir = -1;
          if(tcc.modes[kStepDir]) stepDir = 1;
          short dpt1 = stepDir * (closePt1 - endPt1);
          if(prt) mf::LogVerbatim("TC")<<" endPt1 "<<endPt1<<" closePt1 "<<closePt1<<" dpt1 "<<dpt1<<" doca1 "<<doca1;
          if(dpt1 < -1) continue;
          if(slc.tjs[it1].EndPt[1] > 4) {
            if(dpt1 > 3) continue;
          } else {
            // tighter cut for short trajectories
            if(dpt1 > 2) continue;
          }
          unsigned short closePt2;
          float doca2 = sepCut;
          if(!TrajClosestApproach(tj2, wint, tint, closePt2, doca2)) continue;
          short dpt2 = stepDir * (closePt2 - endPt2);
          if(prt) mf::LogVerbatim("TC")<<" endPt2 "<<endPt2<<" closePt2 "<<closePt2<<" dpt2 "<<dpt2<<" doca2 "<<doca2;
          if(dpt2 < -1) continue;
          if(slc.tjs[it2].EndPt[1] > 4) {
            if(dpt2 > 3) continue;
          } else {
            // tighter cut for short trajectories
            if(dpt2 > 2) continue;
          }
          bool fixVxPos = false;
          // fix the vertex position if there is a charge kink here
          if(tj1.EndFlag[end1][kAtKink]) fixVxPos = true;
          if(prt) mf::LogVerbatim("TC")<<" wint:tint "<<(int)wint<<":"<<(int)(tint/tcc.unitsPerTick)<<" fixVxPos? "<<fixVxPos;
          if(requireVtxTjChg) {
            // ensure that there is a signal between these TPs and the vertex on most of the wires
            bool signalBetween = true;
            short dpt = abs(wint - tp1.Pos[0]);
            if(dpt > 2 && !SignalBetween(slc, tp1, wint, tcc.vtx2DCuts[6])) {
              if(prt) mf::LogVerbatim("TC")<<" Fails SignalBetween for tp1 "<<dpt;
              signalBetween = false;
            }
            dpt = abs(wint - tp2.Pos[0]);
            if(dpt > 2 && !SignalBetween(slc, tp2, wint, tcc.vtx2DCuts[6])) {
              if(prt) mf::LogVerbatim("TC")<<" Fails SignalBetween for tp2 "<<dpt;
              signalBetween = false;
            }
            // consider the case where the intersection point is wrong because the
            // end TP angles are screwed up but the Tjs are close to each other near the end
            if(!signalBetween) {
              unsigned short ipt1, ipt2;
              float maxSep = 3;
              bool isClose = TrajTrajDOCA(slc, tj1, tj2, ipt1, ipt2, maxSep, false);
              // require that they are close at the correct end
              if(isClose) isClose = (abs(ipt1 - endPt1) < 4 && abs(ipt2 - endPt2) < 4);
              if(isClose) {
                if(prt) mf::LogVerbatim("TC")<<" TrajTrajDOCA are close with minSep "<<maxSep<<" near "<<PrintPos(slc, tj1.Pts[ipt1].Pos)<<" "<<PrintPos(slc, tj2.Pts[ipt2].Pos);
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
                if(prt) mf::LogVerbatim("TC")<<" new wint:tint "<<(int)wint<<":"<<(int)(tint/tcc.unitsPerTick);
              } else {
                // closest approach > 3
                continue;
              }
            } // no signal between
          } // requireVtxTjChg
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
          aVtx.Stat[kVxIndPlnNoChg] = !requireVtxTjChg;
          // try to fit it. We need to give it an ID to do that. Take the next
          // available ID
          unsigned short newVtxID = slc.vtxs.size() + 1;
          aVtx.ID = newVtxID;
          tj1.VtxID[end1] = newVtxID;
          tj2.VtxID[end2] = newVtxID;
          if(!FitVertex(slc, aVtx, prt)) {
            tj1.VtxID[end1] = 0;
            tj2.VtxID[end2] = 0;
            continue;
          }
          // check proximity to nearby vertices
          unsigned short mergeMeWithVx = IsCloseToVertex(slc, aVtx);
          if(mergeMeWithVx > 0 && MergeWithVertex(slc, aVtx, mergeMeWithVx)) {
            if(prt) mf::LogVerbatim("TC")<<" Merged with close vertex "<<mergeMeWithVx;
            continue;
          }
          // Save it
          if(!StoreVertex(slc, aVtx)) continue;
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<" New vtx 2V"<<aVtx.ID;
            myprt<<" Tjs "<<tj1.ID<<"_"<<end1<<"-"<<tj2.ID<<"_"<<end2;
            myprt<<" at "<<std::fixed<<std::setprecision(1)<<aVtx.Pos[0]<<":"<<aVtx.Pos[1]/tcc.unitsPerTick;
          }
          AttachAnyTrajToVertex(slc, slc.vtxs.size() - 1, prt);
          SetVx2Score(slc);
        } // it2
      } // end1
    } // it1

    // only call these on the last pass
    if(pass == USHRT_MAX) {
      FindHammerVertices(slc, inCTP);
      FindHammerVertices2(slc, inCTP);
    }

    if(prt) PrintAllTraj("F2DVo", slc, USHRT_MAX, USHRT_MAX);

  } // Find2DVertices

  //////////////////////////////////////////
  bool MergeWithVertex(TCSlice& slc, VtxStore& vx, unsigned short oVxID)
  {
    // Attempts to merge the trajectories attached to vx with an existing 2D vertex
    // referenced by existingVxID. This function doesn't use the existing end0/end1 vertex association.
    // It returns true if the merging was successful in which case the calling function should
    // not store vx. The calling function needs to have set VtxID to vx.ID for tjs that are currently attached
    // to vx. It assumed that vx hasn't yet been pushed onto slc.vtxs

    if(!tcc.useAlg[kVxMerge]) return false;

    bool prt = tcc.dbgVxMerge && tcc.dbgSlc;

    if(oVxID > slc.vtxs.size()) return false;
    auto& oVx = slc.vtxs[oVxID - 1];
    if(vx.CTP != oVx.CTP) return false;

    // get a list of tjs attached to both vertices
    std::vector<int> tjlist = GetVtxTjIDs(slc, vx);
    if(tjlist.empty()) return false;
    std::vector<int> tmp = GetVtxTjIDs(slc, oVx);
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
        auto& tj = slc.tjs[tjid - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == vx.ID) tj.VtxID[end] = oVx.ID;
        } // end
      } // tjid
      if(!FitVertex(slc, oVx, prt)) {
        if(prt) mf::LogVerbatim("TC")<<"MWV: merge failed "<<vx.ID<<" and existing "<<oVx.ID;
        return false;
      }
      return true;
    } // size = 2

    // sort by decreasing length
    std::vector<SortEntry> sortVec(tjlist.size());
    for(unsigned int indx = 0; indx < sortVec.size(); ++indx) {
      sortVec[indx].index = indx;
      auto& tj = slc.tjs[tjlist[indx] - 1];
      sortVec[indx].val = tj.Pts.size();
    } // indx
    std::sort(sortVec.begin(), sortVec.end(), valDecreasing);
    // re-order the list of Tjs
    auto ttl = tjlist;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) tjlist[ii] = ttl[sortVec[ii].index];
    // Create a local vertex using the two longest slc, then add the shorter ones
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
      auto& tj = slc.tjs[tjlist[ii] - 1];
      unsigned short npwc = NumPtsWithCharge(slc, tj, false);
      unsigned short end = CloseEnd(slc, tj, vpos);
      // assume that we will use the end point of the tj
      unsigned short endPt = tj.EndPt[end];
      if(npwc > 6 && tj.Pts[endPt].NTPsFit < 4) {
        if(end == 0) {
          endPt += 3;
        } else {
          endPt -= 3;
        }
        endPt = NearestPtWithChg(slc, tj, endPt);
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
      // stash the ID in Hits
      tjpts[ii].Hits.resize(1, tj.ID);
    } // tjid
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MWV: "<<oVxID;
      myprt<<" Fit TPs";
      for(unsigned short ii = 0; ii < tjpts.size(); ++ii) {
        auto& tjpt = tjpts[ii];
        myprt<<" "<<tjlist[ii]<<"_"<<tjpt.Step<<"_"<<PrintPos(slc, tjpt.Pos);
      }
    } // prt
    // create a subset of the first two for the first fit
    auto fitpts = tjpts;
    fitpts.resize(2);
    if(!FitVertex(slc, aVx, fitpts, prt)) {
      if(prt) mf::LogVerbatim("TC")<<"MWV: first fit failed ";
      return false;
    }
    // Fit and add tjs to the vertex
    bool needsUpdate = false;
    for(unsigned short ii = 2; ii < tjlist.size(); ++ii) {
      fitpts.push_back(tjpts[ii]);
      if(FitVertex(slc, aVx, fitpts, prt)) {
        needsUpdate = false;
      } else {
        // remove the last Tj point and keep going
        fitpts.pop_back();
        needsUpdate = true;
      }
    } // ii

    if(needsUpdate) FitVertex(slc, aVx, fitpts, prt);
    if(prt) mf::LogVerbatim("TC")<<"MWV: done "<<vx.ID<<" and existing "<<oVx.ID;

    // update. Remove old associations
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
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
      auto& tj = slc.tjs[tjpt.Hits[0] - 1];
      if(tj.VtxID[end] != 0) return false;
      tj.VtxID[end] = oVxID;
    } // ii

    // Update oVx
    oVx.Pos = aVx.Pos;
    oVx.PosErr = aVx.PosErr;
    oVx.ChiDOF = aVx.ChiDOF;
    oVx.NTraj = fitpts.size();
    // Update the score and the charge fraction
    SetVx2Score(slc, oVx);
    oVx.Stat[kVxMerged] = true;
    oVx.Stat[kFixed] = false;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MWV: "<<oVxID;
      myprt<<" Done TPs";
      for(unsigned short ii = 0; ii < fitpts.size(); ++ii) {
        auto& tjpt = fitpts[ii];
        myprt<<" "<<tjpt.Hits[0]<<"_"<<tjpt.AngleCode<<"_"<<PrintPos(slc, tjpt.Pos);
      }
    } // prt

    return true;
  } // MergeWithVertex

  //////////////////////////////////////////
  void FindHammerVertices2(TCSlice& slc, const CTP_t& inCTP)
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

    if(!tcc.useAlg[kHamVx2]) return;
    // This wasn't written for test beam events
    if(tcc.modes[kTestBeam]) return;

    bool prt = (tcc.modes[kDebug] && tcc.dbgSlc && tcc.dbgAlg[kHamVx2]);
    if(prt) mf::LogVerbatim("TC")<<"Inside HamVx2";

    for(unsigned short it1 = 0; it1 < slc.tjs.size(); ++it1) {
      if(slc.tjs[it1].CTP != inCTP) continue;
      if(slc.tjs[it1].AlgMod[kKilled] || slc.tjs[it1].AlgMod[kHaloTj]) continue;
      if(slc.tjs[it1].AlgMod[kHamVx]) continue;
      if(slc.tjs[it1].AlgMod[kHamVx2]) continue;
      if(slc.tjs[it1].AlgMod[kJunkTj]) continue;
      if(slc.tjs[it1].PDGCode == 111) continue;
      unsigned short numPtsWithCharge1 = NumPtsWithCharge(slc, slc.tjs[it1], false);
      if(numPtsWithCharge1 < 6) continue;
      // require high MCSMom
      if(slc.tjs[it1].MCSMom < 200) continue;
      // Check each end of tj1
      bool didaSplit = false;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(slc.tjs[it1].VtxID[end1] > 0) continue;
        unsigned short endPt1 = slc.tjs[it1].EndPt[end1];
        for(unsigned short it2 = 0; it2 < slc.tjs.size(); ++it2) {
          if(it1 == it2) continue;
          if(slc.tjs[it2].AlgMod[kKilled] || slc.tjs[it2].AlgMod[kHaloTj]) continue;
          if(slc.tjs[it2].AlgMod[kHamVx]) continue;
          if(slc.tjs[it2].AlgMod[kHamVx2]) continue;
          // require that both be in the same CTP
          if(slc.tjs[it2].CTP != inCTP) continue;
          if(slc.tjs[it2].AlgMod[kShowerLike]) continue;
          if(slc.tjs[it2].AlgMod[kJunkTj]) continue;
          if(slc.tjs[it2].PDGCode == 111) continue;
          unsigned short numPtsWithCharge2 = NumPtsWithCharge(slc, slc.tjs[it2], true);
          if(numPtsWithCharge2 < 6) continue;
          // ignore muon-like tjs
          if(numPtsWithCharge2 > 100 && slc.tjs[it2].MCSMom > 500) continue;
          // Find the minimum separation between tj1 and tj2
          float minDOCA = 5;
          float doca = minDOCA;
          unsigned short closePt2 = 0;
          TrajPointTrajDOCA(slc, slc.tjs[it1].Pts[endPt1], slc.tjs[it2], closePt2, doca);
          if(doca == minDOCA) continue;
          if(prt) {
            mf::LogVerbatim myprt("TC");
            auto& tj1 = slc.tjs[it1];
            auto& tj2 = slc.tjs[it2];
            myprt<<" FHV2 CTP"<<tj1.CTP;
            myprt<<" t"<<tj1.ID<<"_"<<end1<<" MCSMom "<<tj1.MCSMom<<" ChgRMS "<<tj1.ChgRMS;
            myprt<<" split t"<<tj2.ID<<"? MCSMom "<<tj2.MCSMom<<" ChgRMS "<<tj2.ChgRMS;
            myprt<<" doca "<<doca<<" tj2.EndPt[0] "<<tj2.EndPt[0]<<" closePt2 "<<closePt2;
            myprt<<" tj2.EndPt[1] "<<tj2.EndPt[1];
          } // prt
          // ensure that the closest point is not near an end
          if(closePt2 < slc.tjs[it2].EndPt[0] + 3) continue;
          if(closePt2 > slc.tjs[it2].EndPt[1] - 3) continue;
          // Find the intersection point between the tj1 end and tj2 closest Pt
          float wint, tint;
          TrajIntersection(slc.tjs[it1].Pts[endPt1], slc.tjs[it2].Pts[closePt2], wint, tint);
          // make an angle cut
          float dang = DeltaAngle(slc.tjs[it1].Pts[endPt1].Ang, slc.tjs[it2].Pts[closePt2].Ang);
          if(dang < 0.2) continue;
          // ensure that tj1 doesn't cross tj2 but ensuring that the closest point on tj1 is at closePt1
          doca = 5;
          unsigned short closePt1 = 0;
          TrajPointTrajDOCA(slc, slc.tjs[it2].Pts[closePt2], slc.tjs[it1], closePt1, doca);
          if(closePt1 != endPt1) continue;
          if(prt) mf::LogVerbatim("TC")<<" intersection W:T "<<(int)wint<<":"<<(int)(tint/tcc.unitsPerTick)<<" dang "<<dang;
          // Find the point on tj2 that is closest to this point, overwriting closePt
          doca = minDOCA;
          // the point on tj2 that is closest to the intersection point
          unsigned short intPt2;
          TrajClosestApproach(slc.tjs[it2], wint, tint, intPt2, doca);
          if(prt) mf::LogVerbatim("TC")<<" intPt2 on tj2 "<<intPt2<<" at Pos "<<PrintPos(slc, slc.tjs[it2].Pts[intPt2])<<" doca "<<doca;
          if(doca == minDOCA) continue;
          // find the MCSMom for both sections of tj2
          float mcsmom = slc.tjs[it2].MCSMom;
          float mcsmom1 = MCSMom(slc, slc.tjs[it2], slc.tjs[it2].EndPt[0], intPt2);
          float mcsmom2 = MCSMom(slc, slc.tjs[it2], intPt2, slc.tjs[it2].EndPt[1]);
          // require that the both MCSMoms be greater than
          if(prt) mf::LogVerbatim("TC")<<" Check MCSMom after split: mcsmom1 "<<mcsmom1<<" mcsmom2 "<<mcsmom2;
          if(mcsmom1 < mcsmom || mcsmom2 < mcsmom) continue;
          // start scanning for the point on tj2 that has the best IP with the end of tj1 in the direction
          // from closePt2 -> endPt2
          short dir = 1;
          if(intPt2 < closePt2) dir = -1;
          unsigned short nit = 0;
          unsigned short ipt = intPt2;
          float mostChg = slc.tjs[it2].Pts[ipt].Chg;
          if(prt) mf::LogVerbatim("TC")<<" ipt "<<ipt<<" at Pos "<<PrintPos(slc, slc.tjs[it2].Pts[ipt])<<" chg "<<mostChg;
          while(nit < 20) {
            ipt += dir;
            if(ipt < 3 || ipt > slc.tjs[it2].Pts.size() - 4) break;
            float delta = PointTrajDOCA(slc, slc.tjs[it2].Pts[ipt].Pos[0], slc.tjs[it2].Pts[ipt].Pos[1], slc.tjs[it1].Pts[endPt1]);
            float sep = PosSep(slc.tjs[it2].Pts[ipt].Pos, slc.tjs[it1].Pts[endPt1].Pos);
            float dang = delta / sep;
            float pull = dang / slc.tjs[it1].Pts[endPt1].DeltaRMS;
            if(pull < 2 && slc.tjs[it2].Pts[ipt].Chg > mostChg) {
              mostChg = slc.tjs[it2].Pts[ipt].Chg;
              intPt2 = ipt;
            }
          }
          // require a lot of charge in tj2 in this vicinity compared with the average.
          float chgPull = (mostChg - slc.tjs[it2].AveChg) / slc.tjs[it2].ChgRMS;
          if(prt) mf::LogVerbatim("TC")<<" chgPull at intersection point "<<chgPull;
          if(chgPull < 10) continue;
          // require a signal on every wire between it1 and intPt2
          TrajPoint ltp = slc.tjs[it1].Pts[endPt1];
          if(slc.tjs[it1].Pts[endPt1].Pos[0] < -0.4) continue;
          unsigned int fromWire = std::nearbyint(slc.tjs[it1].Pts[endPt1].Pos[0]);
          if(slc.tjs[it2].Pts[intPt2].Pos[0] < -0.4) continue;
          unsigned int toWire =   std::nearbyint(slc.tjs[it2].Pts[intPt2].Pos[0]);
          if(fromWire > toWire) {
            unsigned int tmp = fromWire;
            fromWire = toWire;
            toWire = tmp;
          }
          bool skipIt = false;
          for(unsigned int wire = fromWire + 1; wire < toWire; ++wire) {
            MoveTPToWire(ltp, (float)wire);
            if(!SignalAtTp(ltp)) {
              skipIt = true;
              break;
            }
          } // wire
          if(skipIt) continue;
          // we have a winner
          // create a new vertex
          VtxStore aVtx;
          aVtx.Pos = slc.tjs[it2].Pts[intPt2].Pos;
          aVtx.NTraj = 3;
          aVtx.Pass = slc.tjs[it2].Pass;
          aVtx.Topo = 6;
          aVtx.ChiDOF = 0;
          aVtx.CTP = inCTP;
          aVtx.ID = slc.vtxs.size() + 1;
          unsigned short ivx = slc.vtxs.size();
          if(!StoreVertex(slc, aVtx)) continue;
          if(!SplitTraj(slc, it2, intPt2, ivx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<"FHV2: Failed to split trajectory";
            MakeVertexObsolete("HamVx2", slc, slc.vtxs[ivx], true);
            continue;
          }
          slc.tjs[it1].VtxID[end1] = slc.vtxs[ivx].ID;
          slc.tjs[it1].AlgMod[kHamVx2] = true;
          slc.tjs[it2].AlgMod[kHamVx2] = true;
          unsigned short newTjIndex = slc.tjs.size() - 1;
          slc.tjs[newTjIndex].AlgMod[kHamVx2] = true;
          AttachAnyTrajToVertex(slc, ivx, prt);
          SetVx2Score(slc);
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(slc, it2);
          // and for the new trajectory
          SetPDGCode(slc, newTjIndex);
          if(prt) mf::LogVerbatim("TC")<<" FHV2: New vtx 2V"<<slc.vtxs[ivx].ID<<" Score "<<slc.vtxs[ivx].Score;
          didaSplit = true;
          break;
        } // it2
        if(didaSplit) break;
      } // end1
    } // it1
  } // FindHammerVertices2

  //////////////////////////////////////////
  void FindHammerVertices(TCSlice& slc, const CTP_t& inCTP)
  {
    // Look for a trajectory that intersects another. Split
    // the trajectory and make a vertex. The convention used
    // is shown pictorially here. Trajectory tj1 must be longer
    // than tj2
    // tj2       ------
    // tj1         /
    // tj1        /
    // tj1       /

    if(!tcc.useAlg[kHamVx]) return;

    bool prt = (tcc.modes[kDebug] && tcc.dbgSlc && tcc.dbgAlg[kHamVx]);
    if(prt) {
      mf::LogVerbatim("TC")<<"Inside HamVx inCTP "<<inCTP;
    }

    for(unsigned short it1 = 0; it1 < slc.tjs.size(); ++it1) {
      if(slc.tjs[it1].CTP != inCTP) continue;
      if(slc.tjs[it1].AlgMod[kKilled] || slc.tjs[it1].AlgMod[kHaloTj]) continue;
      if(slc.tjs[it1].AlgMod[kShowerLike]) continue;
      if(slc.tjs[it1].AlgMod[kJunkTj]) continue;
      if(slc.tjs[it1].PDGCode == 111) continue;
      // minimum length requirements
      unsigned short tj1len = slc.tjs[it1].EndPt[1] - slc.tjs[it1].EndPt[0] + 1;
      if(tj1len < 5) continue;
      // require high MCSMom
      if(slc.tjs[it1].MCSMom < 50) continue;
      if(prt) mf::LogVerbatim("TC")<<"FHV: intersection T"<<slc.tjs[it1].ID<<" candidate";
      // Check each end of tj1
      bool didaSplit = false;
      for(unsigned short end1 = 0; end1 < 2; ++end1) {
        // vertex assignment exists?
        if(slc.tjs[it1].VtxID[end1] > 0) continue;
        unsigned short endPt1 = slc.tjs[it1].EndPt[end1];
        for(unsigned short it2 = 0; it2 < slc.tjs.size(); ++it2) {
          if(slc.tjs[it2].CTP != inCTP) continue;
          if(it1 == it2) continue;
          if(slc.tjs[it2].AlgMod[kKilled] || slc.tjs[it2].AlgMod[kHaloTj]) continue;
          if(slc.tjs[it2].AlgMod[kJunkTj]) continue;
          if(slc.tjs[it2].PDGCode == 111) continue;
          // length of tj2 cut
          unsigned short tj2len = slc.tjs[it2].EndPt[1] - slc.tjs[it2].EndPt[0] + 1;
          if(tj2len < 6) continue;
          // ignore very long straight trajectories (probably a cosmic muon)
          if(tj2len > 200 && slc.tjs[it2].PDGCode == 13 && slc.tjs[it2].MCSMom > 600) continue;
          float minDOCA = 5;
          minDOCA /= std::abs(slc.tjs[it1].Pts[endPt1].Dir[0]);
          float doca = minDOCA;
          unsigned short closePt2 = 0;
          TrajPointTrajDOCA(slc, slc.tjs[it1].Pts[endPt1], slc.tjs[it2], closePt2, doca);
          if(doca == minDOCA) continue;
          // ensure that the closest point is not near an end
          if(prt) mf::LogVerbatim("TC")<<"FHV: Candidate T"<<slc.tjs[it1].ID<<" T"<<slc.tjs[it2].ID<<" doca "<<doca<<" tj2.EndPt[0] "<<slc.tjs[it2].EndPt[0]<<" closePt2 "<<closePt2<<" tj2.EndPt[1] "<<slc.tjs[it2].EndPt[1];
          if(closePt2 < slc.tjs[it2].EndPt[0] + 3) continue;
          if(closePt2 > slc.tjs[it2].EndPt[1] - 3) continue;
          // make an angle cut
          float dang = DeltaAngle(slc.tjs[it1].Pts[endPt1].Ang, slc.tjs[it2].Pts[closePt2].Ang);
          if(prt) mf::LogVerbatim("TC")<<" dang "<<dang<<" imposing a hard cut of 0.4 for now ";
          if(dang < 0.4) continue;
          // check the cleanliness in this area
          std::vector<int> tjids(2);
          tjids[0] = slc.tjs[it1].ID;
          tjids[1] = slc.tjs[it2].ID;
          float chgFrac = ChgFracNearPos(slc, slc.tjs[it2].Pts[closePt2].Pos, tjids);
          if(prt) mf::LogVerbatim("TC")<<" chgFrac "<<chgFrac;
          if(chgFrac < 0.9) continue;
          Point2_t vxpos = slc.tjs[it2].Pts[closePt2].Pos;
          // get a better estimate of the vertex position
          TrajIntersection(slc.tjs[it1].Pts[endPt1], slc.tjs[it2].Pts[closePt2], vxpos[0], vxpos[1]);
          // and a better estimate of the point on tj2 where the split should be made
          doca = minDOCA;
          TrajClosestApproach(slc.tjs[it2], vxpos[0], vxpos[1], closePt2, doca);
          if(prt) mf::LogVerbatim("TC")<<" better pos "<<PrintPos(slc, vxpos)<<" new closePt2 "<<closePt2;
          // create a new vertex
          VtxStore aVtx;
          aVtx.Pos = vxpos;
          aVtx.NTraj = 3;
          aVtx.Pass = slc.tjs[it2].Pass;
          aVtx.Topo = 5;
          aVtx.ChiDOF = 0;
          aVtx.CTP = inCTP;
          aVtx.ID = slc.vtxs.size() + 1;
          if(!StoreVertex(slc, aVtx)) continue;
          unsigned short ivx = slc.vtxs.size() - 1;
          if(!SplitTraj(slc, it2, closePt2, ivx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<"FHV: Failed to split trajectory";
            MakeVertexObsolete("HamVx", slc, slc.vtxs[slc.vtxs.size() - 1], true);
            continue;
          }
          slc.tjs[it1].VtxID[end1] = slc.vtxs[ivx].ID;
          slc.tjs[it1].AlgMod[kHamVx] = true;
          slc.tjs[it2].AlgMod[kHamVx] = true;
          unsigned short newTjIndex = slc.tjs.size() - 1;
          slc.tjs[newTjIndex].AlgMod[kHamVx] = true;
          SetVx2Score(slc);
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(slc, it2);
          // and for the new trajectory
          SetPDGCode(slc, newTjIndex);
          if(prt) {
            auto& vx2 = slc.vtxs[ivx];
            mf::LogVerbatim myprt("TC");
            myprt<<" new 2V"<<vx2.ID<<" Score "<<vx2.Score<<" Tjs";
            auto tjlist = GetAssns(slc, "2V", vx2.ID, "T");
            for(auto tid : tjlist) myprt<<" t"<<tid;
          }
          didaSplit = true;
          break;
        } // tj2
        if(didaSplit) break;
      } // end1
    } // tj1

  } // FindHammerVertices

  //////////////////////////////////////////
  void SplitTrajCrossingVertices(TCSlice& slc, CTP_t inCTP)
  {
    // This is kind of self-explanatory...

    if(!tcc.useAlg[kSplitTjCVx]) return;

    if(slc.vtxs.empty()) return;
    if(slc.tjs.empty()) return;

    constexpr float docaCut = 4;

    bool prt = (tcc.modes[kDebug] && tcc.dbgSlc && tcc.dbgAlg[kSplitTjCVx]);
    if(prt) mf::LogVerbatim("TC")<<"Inside SplitTrajCrossingVertices inCTP "<<inCTP;

    geo::PlaneID planeID = DecodeCTP(inCTP);

    unsigned short nTraj = slc.tjs.size();
    for(unsigned short itj = 0; itj < nTraj; ++itj) {
      // NOTE: Don't use a reference variable because it may get lost if the tj is split
      if(slc.tjs[itj].CTP != inCTP) continue;
      // obsolete trajectory
      if(slc.tjs[itj].AlgMod[kKilled] || slc.tjs[itj].AlgMod[kHaloTj]) continue;
      if(slc.tjs[itj].AlgMod[kJunkTj]) continue;
      if(slc.tjs[itj].AlgMod[kSplitTjCVx]) continue;
      // too short
      if(slc.tjs[itj].EndPt[1] < 6) continue;
      for(unsigned short iv = 0; iv < slc.vtxs.size(); ++iv) {
        auto& vx2 = slc.vtxs[iv];
        // obsolete vertex
        if(vx2.NTraj == 0) continue;
        // trajectory already associated with vertex
        if(slc.tjs[itj].VtxID[0] == vx2.ID ||
           slc.tjs[itj].VtxID[1] == vx2.ID) continue;
        // not in the cryostat/tpc/plane
        if(slc.tjs[itj].CTP != vx2.CTP) continue;
        // poor quality
        if(vx2.Score < tcc.vtx2DCuts[7]) continue;
        float doca = docaCut;
        // make the cut significantly larger if the vertex is in a dead
        // wire gap to get the first TP that is just outside the gap.
        if(vx2.Stat[kOnDeadWire]) doca = 100;
        unsigned short closePt = 0;
        if(!TrajClosestApproach(slc.tjs[itj], vx2.Pos[0], vx2.Pos[1], closePt, doca)) continue;
        if(vx2.Stat[kOnDeadWire]) {
          // special handling for vertices in dead wire regions. Find the IP between
          // the closest point on the Tj and the vertex
          doca = PointTrajDOCA(slc, vx2.Pos[0], vx2.Pos[1], slc.tjs[itj].Pts[closePt]);
        }
        if(doca > docaCut) continue;
        if(prt)  mf::LogVerbatim("TC")<<" doca "<<doca<<" btw T"<<slc.tjs[itj].ID<<" and 2V"<<slc.vtxs[iv].ID<<" closePt "<<closePt<<" in plane "<<planeID.Plane;
        // compare the length of the Tjs used to make the vertex with the length of the
        // Tj that we want to split. Don't allow a vertex using very short Tjs to split a long
        // Tj in the 3rd plane
        auto vxtjs = GetVtxTjIDs(slc, vx2);
        if(vxtjs.empty()) continue;
        unsigned short maxPts = 0;
        // ensure that there is a large angle between a Tj already attached to the vertex and the
        // tj that we want to split. We might be considering a delta-ray here
        float maxdang = 0.3;
        float tjAng = slc.tjs[itj].Pts[closePt].Ang;
        for(auto tjid : vxtjs) {
          auto& vtj = slc.tjs[tjid - 1];
          if(vtj.AlgMod[kDeltaRay]) continue;
          unsigned short npwc = NumPtsWithCharge(slc, vtj, false);
          if(npwc > maxPts) maxPts = npwc;
          unsigned short end = 0;
          if(vtj.VtxID[1] == slc.vtxs[iv].ID) end = 1;
          auto& vtp = vtj.Pts[vtj.EndPt[end]];
          float dang = DeltaAngle(vtp.Ang, tjAng);
          if(dang > maxdang) maxdang = dang;
        } // tjid
        // skip this operation if any of the Tjs in the split list are > 3 * maxPts
        maxPts *= 3;
        bool skipit = false;
        if(NumPtsWithCharge(slc, slc.tjs[itj], false) > maxPts && maxPts < 100) skipit = true;
        if(!skipit && maxdang < 0.3) skipit = true;
        if(prt) mf::LogVerbatim("TC")<<"  maxPts "<<maxPts<<" vxtjs[0] "<<vxtjs[0]<<" maxdang "<<maxdang<<" skipit? "<<skipit;
        if(skipit) {
          // kill the vertex?
          if(doca < 1) MakeVertexObsolete("STCV", slc, vx2, true);
          continue;
        }

        // make some adjustments to closePt
        if(vx2.Stat[kOnDeadWire]) {
          // ensure that the tj will be split at the gap. The closePt point may be
          // on the wrong side of it
          auto& closeTP = slc.tjs[itj].Pts[closePt];
          if(slc.tjs[itj].StepDir > 0 && closePt > slc.tjs[itj].EndPt[0]) {
            if(closeTP.Pos[0] > vx2.Pos[0]) --closePt;
          } else if(slc.tjs[itj].StepDir < 0 && closePt < slc.tjs[itj].EndPt[1]) {
            if(closeTP.Pos[0] < vx2.Pos[0]) ++closePt;
          }
        } else {
          // improve closePt based on vertex position
          // check if closePt and EndPt[1] are the two sides of vertex
          // take dot product of closePt-vtx and EndPt[1]-vtx
          if ((slc.tjs[itj].Pts[closePt].Pos[0]-vx2.Pos[0])*(slc.tjs[itj].Pts[slc.tjs[itj].EndPt[1]].Pos[0]-vx2.Pos[0]) + (slc.tjs[itj].Pts[closePt].Pos[1]-vx2.Pos[1])*(slc.tjs[itj].Pts[slc.tjs[itj].EndPt[1]].Pos[1]-vx2.Pos[1]) <0 && closePt < slc.tjs[itj].EndPt[1] - 1) ++closePt;
          else if ((slc.tjs[itj].Pts[closePt].Pos[0]-vx2.Pos[0])*(slc.tjs[itj].Pts[slc.tjs[itj].EndPt[0]].Pos[0]-vx2.Pos[0]) + (slc.tjs[itj].Pts[closePt].Pos[1]-vx2.Pos[1])*(slc.tjs[itj].Pts[slc.tjs[itj].EndPt[0]].Pos[1]-slc.vtxs[iv].Pos[1]) <0 && closePt > slc.tjs[itj].EndPt[0] + 1) --closePt;
        }


        if(prt)  {
          mf::LogVerbatim("TC")<<"Good doca "<<doca<<" btw T"<<slc.tjs[itj].ID<<" and 2V"<<vx2.ID<<" closePt "<<closePt<<" in plane "<<planeID.Plane<<" CTP "<<slc.vtxs[iv].CTP;
          PrintTP("STCV", slc, closePt, 1, slc.tjs[itj].Pass, slc.tjs[itj].Pts[closePt]);
        }
        // ensure that the closest point is not near an end
        if(closePt < slc.tjs[itj].EndPt[0] + 3) continue;
        if(closePt > slc.tjs[itj].EndPt[1] - 3) continue;
        if(!SplitTraj(slc, itj, closePt, iv, prt)) {
          if(prt) mf::LogVerbatim("TC")<<"SplitTrajCrossingVertices: Failed to split trajectory";
          continue;
        }
        slc.tjs[itj].AlgMod[kSplitTjCVx] = true;
        unsigned short newTjIndex = slc.tjs.size() - 1;
        slc.tjs[newTjIndex].AlgMod[kSplitTjCVx] = true;
        // re-fit the vertex position
        FitVertex(slc, vx2, prt);
      } // iv
    } // itj

  } // SplitTrajCrossingVertices

  //////////////////////////////////////
  void Reconcile2Vs(TCSlice& slc)
  {
    // This function is called before Find3DVertices to identify (and possibly reconcile)
    // Tj and 2V inconsistencies using 2D and 3D(?) information
    if(!tcc.useAlg[kReconcile2Vs]) return;
    if(slc.vtxs.empty()) return;

    bool prt = (tcc.dbg2V && tcc.dbgSlc);

    // clusters of 2Vs
    std::vector<std::vector<int>> vx2Cls;

    // iterate over planes
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      // look for 2D vertices that are close to each other
      vx2Cls.clear();
      for(unsigned short ii = 0; ii < slc.vtxs.size() - 1; ++ii) {
        auto& i2v = slc.vtxs[ii];
        if(i2v.ID <= 0) continue;
        if(DecodeCTP(i2v.CTP).Plane != plane) continue;
        for(unsigned short jj = ii + 1; jj < slc.vtxs.size(); ++jj) {
          auto& j2v = slc.vtxs[jj];
          if(j2v.ID <= 0) continue;
          if(DecodeCTP(j2v.CTP).Plane != plane) continue;
          // make rough separation cuts
          float dp0 = std::abs(i2v.Pos[0] - j2v.Pos[0]);
          if(dp0 > 10) continue;
          float dp1 = std::abs(i2v.Pos[1] - j2v.Pos[1]);
          if(dp1 > 10) continue;
          // do a more careful look
          float err = i2v.PosErr[0];
          if(j2v.PosErr[0] > err) err = j2v.PosErr[0];
          float dp0Sig = dp0 / err;
          if(dp0Sig > 4) continue;
          err = i2v.PosErr[1];
          if(j2v.PosErr[1] > err) err = j2v.PosErr[1];
          float dp1Sig = dp1 / err;
          if(dp1Sig > 4) continue;
          // Look for one of the 2V IDs in a cluster
          bool gotit = false;
          for(auto& vx2cls : vx2Cls) {
            bool goti = (std::find(vx2cls.begin(), vx2cls.end(), i2v.ID) != vx2cls.end());
            bool gotj = (std::find(vx2cls.begin(), vx2cls.end(), j2v.ID) != vx2cls.end());
            if(goti && gotj) {
              gotit = true;
              break;
            } else if(goti) {
              vx2cls.push_back(j2v.ID);
              gotit = true;
              break;
            } else if(gotj) {
              gotit = true;
              vx2cls.push_back(i2v.ID);
              break;
            }
          } // vx2cls
          if(!gotit) {
            // start a new cluster with this pair
            std::vector<int> cls(2);
            cls[0] = i2v.ID;
            cls[1] = j2v.ID;
            vx2Cls.push_back(cls);
          } // !gotit
        } // jj
      } // ii
      if(vx2Cls.empty()) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"2V clusters in plane "<<plane;
        for(auto& vx2cls : vx2Cls) {
          myprt<<"\n";
          for(auto vx2id : vx2cls) myprt<<" 2V"<<vx2id;
        } // vx2cls
      } // prt
      for(auto& vx2cls : vx2Cls) {
        Reconcile2VTs(slc, vx2cls, prt);
      } // vx2cls
    } // plane

    // See if any of the vertices have been altered. If so the environment near them,
    // specifically tagging overlapping trajectories, should be re-done
    bool VxEnvironmentNeedsUpdate = false;
    for(auto& vx : slc.vtxs) {
      if(vx.ID <= 0) continue;
      if(!vx.Stat[kVxEnvOK]) VxEnvironmentNeedsUpdate = true;
    } // vx

    if(VxEnvironmentNeedsUpdate) UpdateVxEnvironment(slc);

  } // Reconcile2Vs

  //////////////////////////////////////
  bool Reconcile2VTs(TCSlice& slc, std::vector<int>& vx2cls, bool prt)
  {
    // The 2D vertices IDs in vx2cls were clustered by the calling function. This function
    // checks the T -> 2V assns and possibly changes it. It returns true if an assn is changed
    // or if a vertex in vx2cls is made obsolete, necessitating a change to the list of 2V
    // clusters
    if(vx2cls.size() < 2) return false;

    // Form a list of all Tjs associated with this 2V cluster
    std::vector<int> t2vList;

    CTP_t inCTP;
    for(auto vx2id : vx2cls) {
      auto& vx2 = slc.vtxs[vx2id - 1];
      inCTP = vx2.CTP;
      // vertex clobbered? If so, vertex clustering needs to be re-done
      if(vx2.ID <= 0) return true;
      auto tlist = GetAssns(slc, "2V", vx2.ID, "T");
      for(auto tid : tlist) if(std::find(t2vList.begin(), t2vList.end(), tid) == t2vList.end()) t2vList.push_back(tid);
    } // vx2id
    if(t2vList.size() < 3) return false;

    // Sum the T -> 2V pulls
    float sumPulls = 0;
    float cnt = 0;
    for(auto tid : t2vList) {
      auto& tj = slc.tjs[tid - 1];
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] <= 0) continue;
        if(std::find(vx2cls.begin(), vx2cls.end(), tj.VtxID[end]) == vx2cls.end()) continue;
        auto& vx = slc.vtxs[tj.VtxID[end] - 1];
        unsigned short nearEnd = 1 - FarEnd(slc, tj, vx.Pos);
        unsigned short fitPt = NearbyCleanPt(slc, tj, nearEnd);
        if(fitPt == USHRT_MAX) return false;
        auto& tp = tj.Pts[fitPt];
        sumPulls += TrajPointVertexPull(slc, tp, vx);
        ++cnt;
      } // end
    } // tid

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"R2VTs: cluster:";
      for(auto vid : vx2cls) myprt<<" 2V"<<vid;
      myprt<<" ->";
      for(auto tid : t2vList) myprt<<" T"<<tid;
      myprt<<" sumPulls "<<std::setprecision(2)<<sumPulls<<" cnt "<<cnt;
    } // prt

    // try to fit all Tjs to one vertex. Find the average position of all of the
    // vertices. This will be used to find the end of the Tjs that are closest to the
    // presumed single vertex
    VtxStore oneVx;
    oneVx.CTP = inCTP;
    for(auto vid : vx2cls) {
      auto& vx = slc.vtxs[vid - 1];
      oneVx.Pos[0] += vx.Pos[0];
      oneVx.Pos[1] += vx.Pos[2];
    } // vid
    oneVx.Pos[0] /= vx2cls.size();
    oneVx.Pos[1] /= vx2cls.size();
    std::vector<TrajPoint> oneVxTPs(t2vList.size());
    for(unsigned short itj = 0; itj < t2vList.size(); ++itj) {
      auto& tj = slc.tjs[t2vList[itj] - 1];
      unsigned short nearEnd = 1 - FarEnd(slc, tj, oneVx.Pos);
      unsigned short fitPt = NearbyCleanPt(slc, tj, nearEnd);
      if(fitPt == USHRT_MAX) return false;
      oneVxTPs[itj] = tj.Pts[fitPt];
      // inflate the TP angle angle error if a TP without an overlap wasn't found
      if(oneVxTPs[itj].Environment[kEnvOverlap]) oneVxTPs[itj].AngErr *= 4;
      oneVxTPs[itj].Step = tj.ID;
    } // ii
    if(!FitVertex(slc, oneVx, oneVxTPs, prt)) return false;

    if(oneVx.ChiDOF < 3) {
      // Update the position of the first 2V in the list
      auto& vx = slc.vtxs[vx2cls[0] - 1];
      vx.Pos = oneVx.Pos;
      vx.PosErr = oneVx.PosErr;
      vx.NTraj = t2vList.size();
      vx.ChiDOF = oneVx.ChiDOF;
      vx.Topo = 14;
      // Set a flag that the environment near this vertex (and all other vertices in this slice)
      // should be revisited
      vx.Stat[kVxEnvOK] = false;
      for(unsigned short ivx = 1; ivx < vx2cls.size(); ++ivx) {
        auto& vx = slc.vtxs[vx2cls[ivx] - 1];
        MakeVertexObsolete("R2VTPs", slc, vx, true);
      } // ivx
      // now attach the trajectories
      for(auto tid : t2vList) {
        auto& tj = slc.tjs[tid - 1];
        unsigned short nearEnd = 1 - FarEnd(slc, tj, vx.Pos);
        tj.VtxID[nearEnd] = vx.ID;
      } // tid
      return true;
    } // oneVx.ChiDOF < 3
    return false;
  } // Reconcile2VTs

  //////////////////////////////////////
  void Find3DVertices(TCSlice& slc)
  {
    // Create 3D vertices from 2D vertices. 3D vertices that are matched
    // in all three planes have Vtx2ID > 0 for all planes. This function re-scores all
    // 2D and 3D vertices and flags Tjs that have high-score 3D vertices

    if(tcc.vtx3DCuts[0] < 0) return;
    if(slc.vtxs.size() < 2) return;

    // create a array/vector of 2D vertex indices in each plane
    std::vector<std::vector<unsigned short>> vIndex(3);
    for(unsigned short ivx = 0; ivx < slc.vtxs.size(); ++ivx) {
      // obsolete vertex
      if(slc.vtxs[ivx].ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(slc.vtxs[ivx].CTP);
      if(planeID.TPC != slc.TPCID.TPC) continue;
      unsigned short plane = planeID.Plane;
      if(plane > 2) continue;
      vIndex[plane].push_back(ivx);
    }

    unsigned short vtxInPln = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(vIndex[plane].size() > 0) ++vtxInPln;
    if(vtxInPln < 2) return;

    float thirdPlanedXCut = 2 * tcc.vtx3DCuts[0];
    bool prt = (tcc.dbg3V && tcc.dbgSlc);
    if(prt) {
      mf::LogVerbatim("TC")<<"Inside Find3DVertices. dX cut "<<tcc.vtx3DCuts[0]<<" thirdPlanedXCut "<<thirdPlanedXCut;
    }

    size_t vsize = slc.vtxs.size();
    // vector of 2D vertices -> 3D vertices.
    std::vector<short> vPtr(vsize, -1);
    // fill temp vectors of 2D vertex X and X errors
    std::vector<float> vX(vsize, FLT_MAX);

    for(unsigned short ivx = 0; ivx < vsize; ++ivx) {
      if(slc.vtxs[ivx].ID <= 0) continue;
      if(slc.vtxs[ivx].Score < tcc.vtx2DCuts[7]) continue;
      if(slc.vtxs[ivx].Pos[0] < -0.4) continue;
      geo::PlaneID planeID = DecodeCTP(slc.vtxs[ivx].CTP);
      // Convert 2D vertex time error to X error
      double ticks = slc.vtxs[ivx].Pos[1] / tcc.unitsPerTick;
      vX[ivx]  = tcc.detprop->ConvertTicksToX(ticks, planeID);
    } // ivx

    // temp vector of all 2D vertex matches
    std::vector<Vtx3Store> v3temp;

    unsigned int cstat = slc.TPCID.Cryostat;
    unsigned int tpc = slc.TPCID.TPC;

    TrajPoint tp;
    constexpr float maxSep = 4;
    float maxScore = 0;
    // i, j, k indicates 3 different wire planes
    // compare vertices in each view
    for(unsigned short ipl = 0; ipl < 2; ++ipl) {
      for(unsigned short ii = 0; ii < vIndex[ipl].size(); ++ii) {
        unsigned short ivx = vIndex[ipl][ii];
        if(vX[ivx] == FLT_MAX) continue;
        auto& ivx2 = slc.vtxs[ivx];
        if(ivx2.Pos[0] < -0.4) continue;
        unsigned int iWire = std::nearbyint(ivx2.Pos[0]);
        for(unsigned short jpl = ipl + 1; jpl < 3; ++jpl) {
          for(unsigned short jj = 0; jj < vIndex[jpl].size(); ++jj) {
            unsigned short jvx = vIndex[jpl][jj];
            if(vX[jvx] == FLT_MAX) continue;
            auto& jvx2 = slc.vtxs[jvx];
            if(jvx2.Pos[0] < -0.4) continue;
            unsigned int jWire = std::nearbyint(jvx2.Pos[0]);
            float dX = std::abs(vX[ivx] - vX[jvx]);
            if(dX > tcc.vtx3DCuts[0]) continue;
            if(prt) {
              mf::LogVerbatim("TC")<<"F3DV: ipl "<<ipl<<" i2V"<<ivx2.ID<<" iX "<<vX[ivx]
              <<" jpl "<<jpl<<" j2V"<<jvx2.ID<<" jvX "<<vX[jvx]<<" W:T "<<(int)jvx2.Pos[0]<<":"<<(int)jvx2.Pos[1]<<" dX "<<dX;
            }
            double y = -1000, z = -1000;
            tcc.geom->IntersectionPoint(iWire, jWire, ipl, jpl, cstat, tpc, y, z);
            if(y < slc.yLo || y > slc.yHi || z < slc.zLo || z > slc.zHi) continue;
            unsigned short kpl = 3 - ipl - jpl;
            float kX = 0.5 * (vX[ivx] + vX[jvx]);
            int kWire = -1;
            if(slc.nPlanes > 2) {
              kWire = tcc.geom->WireCoordinate(y, z, kpl, tpc, cstat);
              std::array<int, 2> wireWindow;
              std::array<float, 2> timeWindow;
              wireWindow[0] = kWire - maxSep;
              wireWindow[1] = kWire + maxSep;
              float time = tcc.detprop->ConvertXToTicks(kX, kpl, (int)tpc, (int)cstat) * tcc.unitsPerTick;
              timeWindow[0] = time - maxSep;
              timeWindow[1] = time + maxSep;
              bool hitsNear;
              std::vector<unsigned int> closeHits = FindCloseHits(slc, wireWindow, timeWindow, kpl, kAllHits, true, hitsNear);
              if(prt) {
                mf::LogVerbatim myprt("TC");
                myprt<<" Hits near "<<kpl<<":"<<kWire<<":"<<(int)(time/tcc.unitsPerTick)<<" = ";
                for(auto iht : closeHits) myprt<<" "<<PrintHit(slc.slHits[iht]);
              }
              if(!hitsNear) continue;
            } // 3-plane TPC
            // save this incomplete 3D vertex
            Vtx3Store v3d;
            // give it a non-zero ID so that SetVx3Score returns a valid score
            v3d.ID = 666;
            v3d.Vx2ID.resize(slc.nPlanes);
            v3d.Vx2ID[ipl] = ivx2.ID;
            v3d.Vx2ID[jpl] = jvx2.ID;
            if(slc.nPlanes == 2) v3d.Vx2ID[2] = -1;
            v3d.X = kX;
            // Use XErr to store dX
            v3d.XErr = dX;
            v3d.Y = y;
            v3d.Z = z;
            v3d.Wire = kWire;
            float posError = dX / tcc.vtx3DCuts[0];
            float vxScoreWght = 0;
            SetVx3Score(slc, v3d);
            vxScoreWght = tcc.vtx3DCuts[2] / v3d.Score;
            if(posError < 0.5) posError = 0;
            v3d.Score = posError + vxScoreWght;
            v3d.TPCID = slc.TPCID;
            // push the incomplete vertex onto the list
            v3temp.push_back(v3d);

            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<"  2 Plane match i2V";
              myprt<<slc.vtxs[ivx].ID<<" P:W:T "<<ipl<<":"<<(int)slc.vtxs[ivx].Pos[0]<<":"<<(int)slc.vtxs[ivx].Pos[1];
              myprt<<" j2V"<<slc.vtxs[jvx].ID<<" P:W:T "<<jpl<<":"<<(int)slc.vtxs[jvx].Pos[0]<<":"<<(int)slc.vtxs[jvx].Pos[1];
              myprt<<std::fixed<<std::setprecision(3);
              myprt<<" dX "<<dX<<" posError "<<posError<<" vxScoreWght "<<vxScoreWght<<" Score "<<v3d.Score;
            }

            if(slc.nPlanes == 2) continue;

            // look for a 3 plane match
            for(unsigned short kk = 0; kk < vIndex[kpl].size(); ++kk) {
              unsigned short kvx = vIndex[kpl][kk];
              float dX = std::abs(vX[kvx] - v3d.X);
              // Wire difference error
              float dW = tcc.wirePitch * std::abs(slc.vtxs[kvx].Pos[0] - kWire);
              if(dX > thirdPlanedXCut) continue;
              if(dW > tcc.vtx3DCuts[1]) continue;
              // put the Y,Z difference in YErr and ZErr
              double y = -1000, z = -1000;
              tcc.geom->IntersectionPoint(iWire, kWire, ipl, kpl, cstat, tpc, y, z);
              v3d.YErr = y - v3d.Y;
              v3d.ZErr = z - v3d.Z;
              v3d.Vx2ID[kpl] = slc.vtxs[kvx].ID;
              v3d.Wire = -1;
              // hijack the Score variable to hold the separation^2, weighted by the
              // vertex3DCuts
              dX = (vX[kvx] - v3d.X) / tcc.vtx3DCuts[0];
              float dY = v3d.YErr / tcc.vtx3DCuts[1];
              float dZ = v3d.ZErr / tcc.vtx3DCuts[1];
              posError = dX * dX + dY * dY + dZ * dZ;
              vxScoreWght = 0;
              SetVx3Score(slc, v3d);
              vxScoreWght = tcc.vtx3DCuts[2] / v3d.Score;
              if(posError < 0.5) posError = 0;
              v3d.Score = posError + vxScoreWght;
              if(v3d.Score > maxScore) maxScore = v3d.Score;
              if(prt) mf::LogVerbatim("TC")<<"    k2V"<<kvx+1<<" dX "<<dX<<" dW "<<dW<<" 3D score "<<v3d.Score;
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
        if(slc.nPlanes == 2) {
          mf::LogVerbatim("TC")<<"2V"<<v3.Vx2ID[0]<<" 2V"<<v3.Vx2ID[1]<<" wire "<<v3.Wire<<" Score "<<v3.Score;
        } else {
          mf::LogVerbatim("TC")<<"2V"<<v3.Vx2ID[0]<<" 2V"<<v3.Vx2ID[1]<<" 2V"<<v3.Vx2ID[2]<<" wire "<<v3.Wire<<" Score "<<v3.Score;
        }
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
        for(unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
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
      mf::LogVerbatim myprt("TC");
      myprt<<"v3sel list\n";
      for(auto& v3d : v3sel) {
        for(auto vx2id : v3d.Vx2ID) if(vx2id > 0) myprt<<" 2V"<<vx2id;
        myprt<<" wire "<<v3d.Wire<<" Score "<<v3d.Score;
        myprt<<"\n";
      } // v3d
    } // prt

    // Count the number of incomplete vertices and store
    unsigned short ninc = 0;
    for(auto& vx3 : v3sel) {
      if(slc.nPlanes == 3 && vx3.Wire >= 0) ++ninc;
      vx3.ID = slc.vtx3s.size() + 1;
      ++evt.global3V_UID;
      vx3.UID = evt.global3V_UID;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" 3V"<<vx3.ID;
        for(auto v2id : vx3.Vx2ID) myprt<<" 2V"<<v2id;
        myprt<<" wire "<<vx3.Wire;
      } // prt
      slc.vtx3s.push_back(vx3);
      // make the 2D -> 3D associations
      for(unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
        if(vx3.Vx2ID[ipl] == 0) continue;
        VtxStore& vx2 = slc.vtxs[vx3.Vx2ID[ipl]-1];
        vx2.Vx3ID = vx3.ID;
      } // ipl
    } // ivx

    // Try to complete incomplete vertices
    if(ninc > 0) {
      CompleteIncomplete3DVerticesInGaps(slc);
      CompleteIncomplete3DVertices(slc);
    }

    // Score and flag Tjs that are attached to high-score vertices
    // First remove Tj vertex flags
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      tj.AlgMod[kTjHiVx3Score] = false;
    } // tj
    // Score the 2D vertices
    for(auto& vx2 : slc.vtxs) {
      if(vx2.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(vx2.CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      SetVx2Score(slc, vx2);
    } // vx2
    // and the 3D vertices
    for(auto& vx3 : slc.vtx3s) {
      if(vx3.ID == 0) continue;
      SetVx3Score(slc, vx3);
    } // vx3

  } // Find3DVertices

  //////////////////////////////////////////
  unsigned short TPNearVertex(const TCSlice& slc, const TrajPoint& tp)
  {
    // Returns the index of a vertex if tp is nearby
    for(unsigned short ivx = 0; ivx < slc.vtxs.size(); ++ivx) {
      if(slc.vtxs[ivx].ID == 0) continue;
      if(slc.vtxs[ivx].CTP != tp.CTP) continue;
      if(std::abs(slc.vtxs[ivx].Pos[0] - tp.Pos[0]) > 1.2) continue;
      if(std::abs(slc.vtxs[ivx].Pos[1] - tp.Pos[1]) > 1.2) continue;
      return ivx;
    } // ivx
    return USHRT_MAX;
  } // TPNearVertex

  //////////////////////////////////////////
  bool AttachToAnyVertex(TCSlice& slc, PFPStruct& pfp, float maxSep, bool prt)
  {
    // Attaches to any 3D vertex but doesn't require consistency with
    // PFP -> Tj -> 2V -> 3V assns
    if(pfp.ID <= 0) return false;

    float pLen = Length(pfp);
    if(pLen == 0) return false;

    // save the old assignents and clear them
    //    auto oldVx3ID = pfp.Vx3ID;
    for(unsigned short end = 0; end < 2; ++end) pfp.Vx3ID[end] = 0;
    std::array<Point3_t, 2> endPos;
    endPos[0] = PosAtEnd(pfp, 0);
    endPos[1] = PosAtEnd(pfp, 1);

    std::array<float, 2> foms {{100., 100.}};
    std::array<int, 2> vtxs {{0}};
    for(auto& vx3 : slc.vtx3s) {
      if(vx3.ID <= 0) continue;
      if(vx3.TPCID != pfp.TPCID) continue;
      std::array<float, 2> sep;
      Point3_t vpos = {{vx3.X, vx3.Y, vx3.Z}};
      sep[0] = PosSep(vpos, endPos[0]);
      sep[1] = PosSep(vpos, endPos[1]);
      unsigned short end = 0;
      if(sep[1] < sep[0]) end = 1;
      // ignore if separation is too large
      if(sep[end] > 100) continue;
      // find the direction vector between these points
      auto vpDir = PointDirection(vpos, endPos[end]);
      auto dir = DirAtEnd(pfp, end);
      double dotp = std::abs(DotProd(vpDir, dir));
      float fom = dotp * sep[end];
      if(prt) mf::LogVerbatim("TC")<<"ATAV: P"<<pfp.ID<<" end "<<end<<" 3V"<<vx3.ID<<" sep "<<sep[end]<<" fom "<<fom<<" maxSep "<<maxSep;
      // ignore if separation is too large
      if(sep[end] > maxSep) continue;
      if(fom < foms[end]) {
        foms[end] = fom;
        vtxs[end] = vx3.ID;
      }
    } // vx3
    bool bingo = false;
    for(unsigned short end = 0; end < 2; ++end) {
      if(vtxs[end] == 0) continue;
      if(prt) mf::LogVerbatim("TC")<<"ATAV: set P"<<pfp.ID<<" end "<<end<<" -> 3V"<<vtxs[end];
      pfp.Vx3ID[end] = vtxs[end];
      bingo = true;
    } // end
    return bingo;
  } // AttachToAnyVertex

  //////////////////////////////////////////
  bool AttachAnyVertexToTraj(TCSlice& slc, int tjID, bool prt)
  {
    // Try to attach a tj that is stored in slc.tjs with any vertex
    if(tjID <= 0 || tjID > (int)slc.tjs.size()) return false;
    if(slc.vtxs.empty()) return false;
    auto& tj = slc.tjs[tjID - 1];
    if(tj.AlgMod[kKilled]) return false;
    if(tcc.vtx2DCuts[0] <= 0) return false;

    unsigned short bestVx = USHRT_MAX;
    // Construct a FOM = (TP-Vtx pull) * (TP-Vtx sep + 1) * (Vtx Score).
    // The +1 keeps FOM from being 0
    float bestFOM = 2 * tcc.vtx2DCuts[3] * (tcc.vtx2DCuts[0] + 1) * tcc.vtx2DCuts[7];
    for(unsigned short ivx = 0; ivx < slc.vtxs.size(); ++ivx) {
      auto& vx = slc.vtxs[ivx];
      if(vx.ID == 0) continue;
      if(vx.CTP != tj.CTP) continue;
      // make some rough cuts
      std::array<float, 2> sep;
      sep[0] = PosSep(vx.Pos, tj.Pts[tj.EndPt[0]].Pos);
      sep[1] = PosSep(vx.Pos, tj.Pts[tj.EndPt[1]].Pos);
      unsigned short end = 0;
      if(sep[1] < sep[0]) end = 1;
      if(sep[end] > 100) continue;
      if(tj.VtxID[end] > 0) continue;
      auto& tp = tj.Pts[tj.EndPt[end]];
      // Pad the separation a bit so we don't get zero
      float fom = TrajPointVertexPull(slc, tp, vx) * (sep[end] + 1) * vx.Score;
      if(fom > bestFOM) continue;
      if(prt) mf::LogVerbatim("TC")<<"AAVTT: T"<<tjID<<" 2V"<<vx.ID<<" FOM "<<fom<<" cut "<<bestFOM;
      bestVx = ivx;
      bestFOM = fom;
    } // vx
    if(bestVx > slc.vtxs.size() - 1) return false;
    auto& vx = slc.vtxs[bestVx];
    return AttachTrajToVertex(slc, tj, vx, prt);
  } // AttachAnyVertexToTraj

  //////////////////////////////////////////
  bool AttachAnyTrajToVertex(TCSlice& slc, unsigned short ivx, bool prt)
  {

    if(ivx > slc.vtxs.size() - 1) return false;
    if(slc.vtxs[ivx].ID == 0) return false;
    if(tcc.vtx2DCuts[0] < 0) return false;

    VtxStore& vx = slc.vtxs[ivx];
    // Hammer vertices should be isolated and clean
    if(vx.Topo == 5 || vx.Topo == 6) return false;

    unsigned short bestTj = USHRT_MAX;
    // Construct a FOM = (TP-Vtx pull) * (TP-Vtx sep + 1).
    // The +1 keeps FOM from being 0
    float bestFOM = 2 * tcc.vtx2DCuts[3] * (tcc.vtx2DCuts[0] + 1);
    for(unsigned int itj = 0; itj < slc.tjs.size(); ++itj) {
      auto& tj = slc.tjs[itj];
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      if(tj.CTP != vx.CTP) continue;
      // make some rough cuts
      std::array<float, 2> sep;
      sep[0] = PosSep(vx.Pos, tj.Pts[tj.EndPt[0]].Pos);
      sep[1] = PosSep(vx.Pos, tj.Pts[tj.EndPt[1]].Pos);
      unsigned short end = 0;
      if(sep[1] < sep[0]) end = 1;
      if(sep[end] > 100) continue;
      if(tj.VtxID[end] > 0) continue;
      auto& tp = tj.Pts[tj.EndPt[end]];
      // Pad the separation a bit so we don't get zero
      float fom = TrajPointVertexPull(slc, tp, vx) * (sep[end] + 1);
      if(fom > bestFOM) continue;
      if(prt) {
        mf::LogVerbatim("TC")<<"AATTV: T"<<tj.ID<<" 2V"<<vx.ID<<" Topo "<<vx.Topo<<" FOM "<<fom<<" cut "<<bestFOM;
      }
      bestTj = itj;
      bestFOM = fom;
    } // tj
    if(bestTj > slc.tjs.size() - 1) return false;
    auto& tj = slc.tjs[bestTj];
    return AttachTrajToVertex(slc, tj, vx, prt);
  } // AttachAnyTrajToVertex

  //////////////////////////////////////////
  bool AttachTrajToVertex(TCSlice& slc, Trajectory& tj, VtxStore& vx, bool prt)
  {
    // Note that this function does not require a signal between the end of the Tj and the vertex

    // tcc.vtx2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction

    if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) return false;
    if(tj.CTP != vx.CTP) return false;
    // already attached?
    if(tj.VtxID[0] == vx.ID || tj.VtxID[1] == vx.ID) return false;

    unsigned short maxShortTjLen = tcc.vtx2DCuts[0];
    // square the separation cut to simplify testing in the loop
    float maxSepCutShort2 = tcc.vtx2DCuts[1] * tcc.vtx2DCuts[1];
    float maxSepCutLong2 = tcc.vtx2DCuts[2] * tcc.vtx2DCuts[2];

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
    // use the short Tj cut if the trajectory looks like an electron
    if(!tjShort && tj.ChgRMS > 0.5) tjShort = true;
    float closestApproach;
    // ignore bad separation between the closest tj end and the vertex
    if(tjShort) {
      if(vtxTjSep2 > maxSepCutShort2) return false;
      closestApproach = tcc.vtx2DCuts[1];
    } else {
      closestApproach = tcc.vtx2DCuts[2];
      if(vtxTjSep2 > maxSepCutLong2) return false;
    }

    // Calculate the pull on the vertex
    TrajPoint& tp = tj.Pts[tj.EndPt[end]];
    float tpVxPull = TrajPointVertexPull(slc, tp, vx);
    bool signalBetween = SignalBetween(slc, tp, vx.Pos[0], 0.8);

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
    if(length > 4 && length < closestApproach) return false;

    float pullCut = tcc.vtx2DCuts[3];
    // Dec 21, 2017 Loosen up the pull cut for short close slc. These are likely to
    // be poorly reconstructed. It is better to have them associated with the vertex
    // than not.
    if(tjShort) pullCut *= 2;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"ATTV: 2V"<<vx.ID;
      myprt<<" NTraj "<<vx.NTraj;
      myprt<<" oldTJs";
      for(unsigned short itj = 0; itj < slc.tjs.size(); ++itj) {
        Trajectory& tj = slc.tjs[itj];
        if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
        if(tj.CTP != vx.CTP) continue;
        if(tj.VtxID[0] == vx.ID) myprt<<" T"<<tj.ID<<"_0";
        if(tj.VtxID[1] == vx.ID) myprt<<" T"<<tj.ID<<"_1";
      }
      myprt<<" +  T"<<tj.ID<<"_"<<end<<" vtxTjSep "<<sqrt(vtxTjSep2)<<" tpVxPull "<<tpVxPull<<" pullCut "<<pullCut<<" dpt "<<dpt;
    }
    if(tpVxPull > pullCut) return false;
    if(dpt > 2) return true;

    // remove the fixed position flag if there are more than 2 tjs
    bool fixedBit = vx.Stat[kFixed];
    if(fixedBit && vx.NTraj < 2) vx.Stat[kFixed] = false;

    // Passed all the cuts. Attach it to the vertex and try a fit
    tj.VtxID[end] = vx.ID;
    // flag as a photon Tj so it isn't included in the fit
    tj.AlgMod[kPhoton] = !signalBetween;
    // make a copy of the vertex and fit it
    auto vxTmp = vx;
    if(FitVertex(slc, vxTmp, prt)) {
      SetVx2Score(slc, vxTmp);
      if(prt) mf::LogVerbatim("TC")<<" Success";
      vx = vxTmp;
      return true;
    }
    // Keep the Tj -> Vx assn since we got this far, but don't include this end of the Tj in the fit
    tj.EndFlag[end][kNoFitVx] = true;
    if(prt) mf::LogVerbatim("TC")<<" Poor fit. Keep T"<<tj.ID<<"-2V"<<vx.ID<<" assn with kNoFitVx";
    // restore the fixed flag
    vx.Stat[kFixed] = fixedBit;
    return true;
  } // AttachTrajToVertex

  /////////////////////////////////////////
  float TrajPointVertexPull(const TCSlice& slc, const TrajPoint& tp, const VtxStore& vx)
  {
    // Calculates the position pull between a trajectory point and a vertex

    // impact parameter between them
    double ip = PointTrajDOCA(slc, vx.Pos[0], vx.Pos[1], tp);
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
  float VertexVertexPull(const TCSlice& slc, const Vtx3Store& vx1, const Vtx3Store& vx2)
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
  float VertexVertexPull(const TCSlice& slc, const VtxStore& vx1, const VtxStore& vx2)
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
  bool StoreVertex(TCSlice& slc, VtxStore& vx)
  {
    // jacket around the push to ensure that the Tj and vtx CTP is consistent.
    // The calling function should score the vertex after the trajectories are attached

    if(vx.ID != int(slc.vtxs.size() + 1)) return false;

    ++evt.global2V_UID;
    vx.UID = evt.global2V_UID;

    unsigned short nvxtj = 0;
    unsigned short nok = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      if(vx.ID == tj.VtxID[0] || vx.ID == tj.VtxID[1]) ++nvxtj;
      if(vx.CTP != tj.CTP) continue;
      if(vx.ID == tj.VtxID[0] || vx.ID == tj.VtxID[1]) ++nok;
    } // tj

    if(nok != nvxtj) {
      mf::LogVerbatim("TC")<<"StoreVertex: vertex "<<vx.ID<<" Topo "<<vx.Topo<<" has inconsistent CTP code "<<vx.CTP<<" with one or more Tjs\n";
      for(auto& tj : slc.tjs) {
        if(tj.AlgMod[kKilled]) continue;
        if(tj.VtxID[0] == vx.ID) tj.VtxID[0] = 0;
        if(tj.VtxID[1] == vx.ID) tj.VtxID[1] = 0;
      }
      return false;
    }
    vx.NTraj = nok;
    slc.vtxs.push_back(vx);
    return true;
  } // StoreVertex

  /////////////////////////////////////////
  bool FitVertex(TCSlice& slc, VtxStore& vx, bool prt)
  {
    // Fit the vertex using T -> 2V assns

    // tcc.vtx2DCuts fcl input usage
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
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      if(tj.CTP != vx.CTP) continue;
      if(tj.AlgMod[kPhoton]) continue;
      bool added = false;
      if(tj.VtxID[0] == vx.ID && !tj.EndFlag[0][kNoFitVx]) {
        vxTp.push_back(tj.Pts[tj.EndPt[0]]);
        added = true;
      }
      if(tj.VtxID[1] == vx.ID && !tj.EndFlag[1][kNoFitVx]) {
        vxTp.push_back(tj.Pts[tj.EndPt[1]]);
        added = true;
      }
      // stash the ID in Step to help debugging
      if(added) {
        auto& tp = vxTp[vxTp.size()-1];
        if(tj.ID > 0) tp.Step = (int)tj.ID;
        // inflate the angle errors for Tjs with few fitted points
        if(tp.NTPsFit < 4) tp.AngErr *= 4;
      }
    } // tj

    bool success = FitVertex(slc, vx, vxTp, prt);

    if(!success) return false;
    return true;
  } // FitVertex

  /////////////////////////////////////////
  bool FitVertex(TCSlice& slc, VtxStore& vx, std::vector<TrajPoint>& vxTPs, bool prt)
  {
    // Version with LSQ fit. Each TP position (P0,P1) and slope S are fit to a vertex
    // at position (V0, V1), using the equation P1 = V1 + (P0 - V0) * S. This is put
    // in the form A * V = b. V is found using V = (A^T * A)^-1 * A^T * b. This is
    // usually done using the TDecompSVD Solve method but here we do it step-wise to
    // get access to the covariance matrix (A^T * A)^-1. The pull between the TP position
    // and the vertex position is stored in tp.Delta

    if(vxTPs.size() < 2) return false;
    if(vxTPs.size() == 2) {
      vx.ChiDOF = 0.;
      return true;
    }

    unsigned short npts = vxTPs.size();
    TMatrixD A(npts, 2);
    TVectorD b(npts);
    for(unsigned short itj = 0; itj < vxTPs.size(); ++itj) {
      auto& tp = vxTPs[itj];
      double dtdw = tp.Dir[1] / tp.Dir[0];
      double wt = 1 / (tp.AngErr * tp.AngErr);
      A(itj, 0) = -dtdw * wt;
      A(itj, 1) = 1. * wt;
      b(itj) = (tp.Pos[1] - tp.Pos[0] * dtdw) * wt;
    } // itj

    TMatrixD AT(2, npts);
    AT.Transpose(A);
    TMatrixD ATA = AT * A;
    double *det = 0;
    try{ ATA.Invert(det); }
    catch(...) { return false; }
    if(det == NULL) return false;
    TVectorD vxPos = ATA * AT * b;
    vx.PosErr[0] = sqrt(ATA[0][0]);
    vx.PosErr[1] = sqrt(ATA[1][1]);
    vx.Pos[0] = vxPos[0];
    vx.Pos[1] = vxPos[1];

    // Calculate Chisq
    vx.ChiDOF = 0;
    if(vxTPs.size() > 2) {
      for(auto& tp : vxTPs) {
        // highjack TP Delta for the vertex pull
        tp.Delta = TrajPointVertexPull(slc, tp, vx);
        vx.ChiDOF += tp.Delta;
      } // itj
      vx.ChiDOF /= (float)(vxTPs.size() - 2);
    } // vxTPs.size() > 2

    if(prt) {
      mf::LogVerbatim("TC")<<"Note: TP - 2V pull is stored in TP.Delta";
      PrintTPHeader("FV");
      for(auto& tp : vxTPs) PrintTP("FV", slc, 0, 1, 1, tp);
    }

    return true;
  } // FitVertex

//////////////////////////////////////////
  bool ChkVtxAssociations(TCSlice& slc, const CTP_t& inCTP)
  {
    // Check the associations

    // check the 2D -> 3D associations
    geo::PlaneID planeID = DecodeCTP(inCTP);
    unsigned short plane = planeID.Plane;
    for(auto& vx2 : slc.vtxs) {
      if(vx2.CTP != inCTP) continue;
      if(vx2.ID == 0) continue;
      if(vx2.Vx3ID == 0) continue;
      if(vx2.Vx3ID > int(slc.vtx3s.size())) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: Invalid vx2.Vx3ID "<<vx2.Vx3ID<<" in 2D vtx "<<vx2.ID;
        return false;
      }
      auto& vx3 = slc.vtx3s[vx2.Vx3ID-1];
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
    for(auto& vx3 : slc.vtx3s) {
      if(vx3.ID == 0) continue;
      if(vx3.Vx2ID[plane] == 0) continue;
      if(vx3.Vx2ID[plane] > (int)slc.vtxs.size()) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: Invalid vx3.Vx2ID "<<vx3.Vx2ID[plane]<<" in CTP "<<inCTP;
        return false;
      }
      auto& vx2 = slc.vtxs[vx3.Vx2ID[plane]-1];
      if(vx2.Vx3ID != vx3.ID) {
        mf::LogVerbatim("TC")<<"ChkVtxAssociations: 3V"<<vx3.ID<<" thinks it is matched to 2V"<<vx2.ID<<" but vx2 says no!";
        return false;
      }
    } // vx3

    // check the Tj -> 2D associations
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == 0) continue;
        if(tj.VtxID[end] > slc.vtxs.size()) {
          mf::LogVerbatim("TC")<<"ChkVtxAssociations: T"<<tj.ID<<" thinks it is matched to 2V"<<tj.VtxID[end]<<" on end "<<end<<" but no vertex exists. Recovering";
          tj.VtxID[end] = 0;
          return false;
        }
        unsigned short ivx = tj.VtxID[end] - 1;
        auto& vx2 = slc.vtxs[ivx];
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
  void ScoreVertices(TCSlice& slc)
  {
    // reset all 3D vertex, 2D vertex and Tj high-score vertex bits in tpcid

    // reset the 2D vertex status bits
    for(auto& vx : slc.vtxs) {
      if(vx.ID == 0) continue;
      vx.Stat[kHiVx3Score] = false;
    } // vx
    // and the tj bits
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      tj.AlgMod[kTjHiVx3Score] = false;
    } // tj
    // Score the 2D vertices
    for(auto& vx : slc.vtxs) {
      if(vx.ID == 0) continue;
      SetVx2Score(slc, vx);
    } // vx
    // Score the 3D vertices
    for(auto& vx3 : slc.vtx3s) {
      if(vx3.ID == 0) continue;
      SetVx3Score(slc, vx3);
    } // vx3
  } // ScoreVertices

  //////////////////////////////////////////
  void KillPoorVertices(TCSlice& slc)
  {
    // kill 2D vertices that have low score and are not attached to a high-score 3D vertex
    if(slc.vtxs.empty()) return;
    for(auto& vx : slc.vtxs) {
      if(vx.ID == 0) continue;
      if(vx.Score > tcc.vtx2DCuts[7]) continue;
      if(vx.Vx3ID > 0) {
        auto& vx3 = slc.vtx3s[vx.Vx3ID - 1];
        if(vx3.Primary) continue;
        if(slc.vtx3s[vx.Vx3ID - 1].Score >= tcc.vtx2DCuts[7]) continue;
      }
      MakeVertexObsolete("KPV", slc, vx, false);
    } // vx

  } // KillPoorVertices

  //////////////////////////////////////////
  void SetHighScoreBits(TCSlice& slc, Vtx3Store& vx3)
  {
    // Sets the tj and 2D vertex score bits to true

    if(vx3.ID == 0) return;

    for(unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
      if(vx3.Vx2ID[ipl] <= 0) continue;
      VtxStore& vx2 = slc.vtxs[vx3.Vx2ID[ipl] - 1];
      vx2.Stat[kHiVx3Score] = false;
      // transfer this to all attached tjs and vertices attached to those tjs
      std::vector<int> tjlist = GetVtxTjIDs(slc, vx2);
      std::vector<int> vxlist;
      while(true) {
        // tag Tjs and make a list of attached vertices whose high-score
        // bit needs to be set
        vxlist.clear();
        for(auto tjid : tjlist) {
          auto& tj = slc.tjs[tjid - 1];
          tj.AlgMod[kTjHiVx3Score] = true;
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] == 0) continue;
            auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
            if(vx2.Stat[kHiVx3Score]) continue;
            vx2.Stat[kHiVx3Score] = true;
            vxlist.push_back(vx2.ID);
          } // end
        } // tjid

        if(vxlist.empty()) break;
        // re-build tjlist using vxlist
        std::vector<int> newtjlist;
        for(auto vxid : vxlist) {
          auto& vx2 = slc.vtxs[vxid - 1];
          auto tmp = GetVtxTjIDs(slc, vx2);
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
  void SetVx3Score(TCSlice& slc, Vtx3Store& vx3)
  {
    // Calculate the 3D vertex score and flag Tjs that are attached to high score vertices as defined
    // by vtx2DCuts

    if(vx3.ID == 0) return;

    vx3.Score = 0;
    for(unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
      if(vx3.Vx2ID[ipl] <= 0) continue;
      VtxStore& vx2 = slc.vtxs[vx3.Vx2ID[ipl] - 1];
      vx3.Score += vx2.Score;
    } // ipl
    vx3.Score /= (float)slc.nPlanes;
    // don't allow it to get too small or negative
    if(vx3.Score < 0.001) vx3.Score = 0.001;
    if(vx3.Score > tcc.vtx2DCuts[7]) SetHighScoreBits(slc, vx3);

  } // SetVx3Score

  //////////////////////////////////////////
  void SetVx2Score(TCSlice& slc)
  {
    // A version that sets the score of the last added vertex
    if(slc.vtxs.empty()) return;
    auto& vx2 = slc.vtxs[slc.vtxs.size() - 1];
    SetVx2Score(slc, vx2);
  } // SetVx2Score

  //////////////////////////////////////////
  void SetVx2Score(TCSlice& slc, VtxStore& vx2)
  {
    // Calculate the 2D vertex score
    if(vx2.ID == 0) return;

    // Don't score vertices from CheckTrajBeginChg, MakeJunkVertices or Neutral vertices. Set to the minimum
    if(vx2.Topo == 8 || vx2.Topo == 9 || vx2.Topo == 11 || vx2.Topo == 12) {
      vx2.Score = tcc.vtx2DCuts[7] + 0.1;
      auto vtxTjID = GetVtxTjIDs(slc, vx2);
      vx2.TjChgFrac = ChgFracNearPos(slc, vx2.Pos, vtxTjID);
      return;
    }

    // Cuts on Tjs attached to vertices
    constexpr float maxChgRMS = 0.25;
    constexpr float momBin = 50;

    vx2.Score = -1000;
    vx2.TjChgFrac = 0;
    if(vx2.ID == 0) return;
    if(tcc.vtxScoreWeights.size() < 4) return;

    auto vtxTjIDs = GetVtxTjIDs(slc, vx2);
    if(vtxTjIDs.empty()) return;

    // Vertex position error
    float vpeScore = -tcc.vtxScoreWeights[0] * (vx2.PosErr[0] + vx2.PosErr[1]);

    unsigned short m3Dcnt = 0;
    if(vx2.Vx3ID > 0) {
      m3Dcnt = 1;
      // Add another if the 3D vertex is complete
      unsigned short ivx3 = vx2.Vx3ID - 1;
      if(slc.vtx3s[ivx3].Wire < 0) m3Dcnt = 2;
    }
    float m3DScore = tcc.vtxScoreWeights[1] * m3Dcnt;

    vx2.TjChgFrac = ChgFracNearPos(slc, vx2.Pos, vtxTjIDs);
    float cfScore = tcc.vtxScoreWeights[2] * vx2.TjChgFrac;

    // Define a weight for each Tj
    std::vector<int> tjids;
    std::vector<float> tjwts;
    unsigned short cnt13 = 0;
    for(auto tjid : vtxTjIDs) {
      Trajectory& tj = slc.tjs[tjid - 1];
      // Feb 22 Ignore short Tjs and junk tjs
      if(tj.AlgMod[kJunkTj]) continue;
      unsigned short lenth = tj.EndPt[1] - tj.EndPt[0] + 1;
      if(lenth < 3) continue;
      float wght = (float)tj.MCSMom / momBin;
      // weight by the first tagged muon
      if(tj.PDGCode == 13) {
        ++cnt13;
        if(cnt13 == 1) wght *= 2;
      }
      // weight by charge rms
      if(tj.ChgRMS < maxChgRMS) ++wght;
      // Shower Tj
      if(tj.AlgMod[kShowerTj]) ++wght;
      // ShowerLike
      if(tj.AlgMod[kShowerLike]) --wght;
      tjids.push_back(tjid);
      tjwts.push_back(wght);
    } // tjid

    if(tjids.empty()) return;

    float tjScore = 0;
    float sum = 0;
    float cnt = 0;
    for(unsigned short it1 = 0; it1 < tjids.size() - 1; ++it1) {
      Trajectory& tj1 = slc.tjs[tjids[it1] - 1];
      float wght1 = tjwts[it1];
      // the end that has a vertex
      unsigned short end1 = 0;
      if(tj1.VtxID[1] == vx2.ID) end1 = 1;
      unsigned short endPt1 = tj1.EndPt[end1];
      // bump up the weight if there is a Bragg peak at the other end
      unsigned short oend1 = 1 - end1;
      if(tj1.EndFlag[oend1][kBragg]) ++wght1;
      float ang1 = tj1.Pts[endPt1].Ang;
      float ang1Err2 = tj1.Pts[endPt1].AngErr * tj1.Pts[endPt1].AngErr;
      for(unsigned short it2 = it1 + 1; it2 < tjids.size(); ++it2) {
        Trajectory& tj2 = slc.tjs[tjids[it2] - 1];
        float wght2 = tjwts[it2];
        unsigned end2 = 0;
        if(tj2.VtxID[1] == vx2.ID) end2 = 1;
        // bump up the weight if there is a Bragg peak at the other end
        unsigned short oend2 = 1 - end2;
        if(tj2.EndFlag[oend2][kBragg]) ++wght2;
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
      tjScore = tcc.vtxScoreWeights[3] * sum;
    }
    vx2.Score = vpeScore + m3DScore + cfScore + tjScore;
    if(tcc.dbg2V && tcc.dbgSlc && vx2.CTP == debug.CTP) {
      // last call after vertices have been matched to the truth. Use to optimize vtxScoreWeights using
      // an ntuple
      mf::LogVerbatim myprt("TC");
      bool printHeader = true;
      Print2V("SVx2S", myprt, vx2, printHeader);
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<" vpeScore "<<vpeScore<<" m3DScore "<<m3DScore;
      myprt<<" cfScore "<<cfScore<<" tjScore "<<tjScore;
      myprt<<" Score "<<vx2.Score;
    }
  } // SetVx2Score

  //////////////////////////////////////////
  void CompleteIncomplete3DVerticesInGaps(TCSlice& slc)
  {

    if(!tcc.useAlg[kComp3DVxIG]) return;
    if(slc.nPlanes != 3) return;

    bool prt = (tcc.modes[kDebug] && tcc.dbgSlc && tcc.dbgAlg[kComp3DVxIG]);
    if(prt) mf::LogVerbatim("TC")<<"Inside CI3DVIG:";

    for(unsigned short iv3 = 0; iv3 < slc.vtx3s.size(); ++iv3) {
      Vtx3Store& vx3 = slc.vtx3s[iv3];
      // ignore obsolete vertices
      if(vx3.ID == 0) continue;
      // check for a completed 3D vertex
      if(vx3.Wire < 0) continue;
      unsigned short mPlane = USHRT_MAX;
      for(unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
        if(vx3.Vx2ID[ipl] > 0) continue;
        mPlane = ipl;
        break;
      } // ipl
      if(mPlane == USHRT_MAX) continue;
      CTP_t mCTP = EncodeCTP(vx3.TPCID.Cryostat, vx3.TPCID.TPC, mPlane);
      // require that the missing vertex be in a large block of dead wires
      float dwc = DeadWireCount(slc, vx3.Wire - 4, vx3.Wire + 4, mCTP);
      if(dwc < 5) continue;
      // X position of the purported missing vertex
      VtxStore aVtx;
      aVtx.ID = slc.vtxs.size() + 1;
      aVtx.Pos[0] = vx3.Wire;
      aVtx.Pos[1] = tcc.detprop->ConvertXToTicks(vx3.X, mPlane, vx3.TPCID.TPC, vx3.TPCID.Cryostat) * tcc.unitsPerTick;
      aVtx.CTP = mCTP;
      aVtx.Topo = 4;
      aVtx.NTraj = 0;
      // Give it a bogus pass to indicate it wasn't created while stepping
      aVtx.Pass = 9;
      if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: Incomplete vertex "<<iv3<<" in plane "<<mPlane<<" wire "<<vx3.Wire<<" Made 2D vertex ";
      std::vector<int> tjIDs;
      std::vector<unsigned short> tjEnds;
      for(unsigned short itj = 0; itj < slc.tjs.size(); ++itj) {
        if(slc.tjs[itj].CTP != mCTP) continue;
        if(slc.tjs[itj].AlgMod[kKilled] || slc.tjs[itj].AlgMod[kHaloTj]) continue;
        for(unsigned short end = 0; end < 2; ++end) {
          unsigned short ept = slc.tjs[itj].EndPt[end];
          TrajPoint& tp = slc.tjs[itj].Pts[ept];
          unsigned short oept = slc.tjs[itj].EndPt[1 - end];
          TrajPoint& otp = slc.tjs[itj].Pts[oept];
          // ensure that this is the end closest to the vertex
          if(std::abs(tp.Pos[0] - aVtx.Pos[0]) > std::abs(otp.Pos[0] - aVtx.Pos[0])) continue;
          float doca = PointTrajDOCA(slc, aVtx.Pos[0], aVtx.Pos[1], tp);
          if(doca > 2) continue;
          float dwc = DeadWireCount(slc, aVtx.Pos[0], tp.Pos[0], tp.CTP);
          float ptSep;
          if(aVtx.Pos[0] > tp.Pos[0]) {
            ptSep = aVtx.Pos[0] - tp.Pos[0] - dwc;
          } else {
            ptSep = tp.Pos[0] - aVtx.Pos[0] - dwc;
          }
          if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: tj ID "<<slc.tjs[itj].ID<<" doca "<<doca<<" ptSep "<<ptSep;
          if(ptSep < -2 || ptSep > 2) continue;
          // don't clobber an existing association
          if(slc.tjs[itj].VtxID[end] > 0) continue;
          tjIDs.push_back(slc.tjs[itj].ID);
          tjEnds.push_back(end);
        } // end
      } // itj
      if(!tjIDs.empty()) {
        // Determine how messy this region is
        aVtx.TjChgFrac = ChgFracNearPos(slc, aVtx.Pos, tjIDs);
        if(aVtx.TjChgFrac < 0.7) continue;
        aVtx.Vx3ID = vx3.ID;
        // Save the 2D vertex
        if(!StoreVertex(slc, aVtx)) continue;
        for(unsigned short ii = 0; ii < tjIDs.size(); ++ii) {
          unsigned short itj = tjIDs[ii] - 1;
          slc.tjs[itj].VtxID[tjEnds[ii]] = aVtx.ID;
          slc.tjs[itj].AlgMod[kComp3DVxIG] = true;
        } // ii
        SetVx2Score(slc);
        vx3.Vx2ID[mPlane] = aVtx.ID;
        vx3.Wire = -1;
        if(prt) mf::LogVerbatim("TC")<<"CI3DVIG: new vtx 2V"<<aVtx.ID<<" points to 3V"<<vx3.ID;
      }
    } // vx3

  } // CompleteIncomplete3DVerticesInGaps

  //////////////////////////////////////////
  void CompleteIncomplete3DVertices(TCSlice& slc)
  {
    // Look for trajectories in a plane that lack a 2D vertex as listed in
    // 2DVtxID that are near the projected wire. This may trigger splitting trajectories,
    // assigning them to a new 2D vertex and completing 3D vertices

    if(!tcc.useAlg[kComp3DVx]) return;
    if(slc.nPlanes != 3) return;

    bool prt = (tcc.modes[kDebug] && tcc.dbgSlc && tcc.dbgAlg[kComp3DVx]);

    float maxdoca = 3;
    if(prt) mf::LogVerbatim("TC")<<"Inside CI3DV with maxdoca set to "<<maxdoca;
    unsigned short ivx3 = 0;
    for(auto& vx3 : slc.vtx3s) {
      // ignore obsolete vertices
      if(vx3.ID == 0) continue;
      // check for a completed 3D vertex
      if(vx3.Wire < 0) continue;
      unsigned short mPlane = USHRT_MAX;
      // look for vertices in the induction plane in which the charge requirement wasn't imposed
      bool indPlnNoChgVtx = false;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        if(vx3.Vx2ID[plane] > 0) {
          auto& vx2 = slc.vtxs[vx3.Vx2ID[plane] - 1];
          if(vx2.Stat[kVxIndPlnNoChg]) indPlnNoChgVtx = true;
          continue;
        }
        mPlane = plane;
      } // ipl
      if(mPlane == USHRT_MAX) continue;
      if(indPlnNoChgVtx) continue;
      CTP_t mCTP = EncodeCTP(vx3.TPCID.Cryostat, vx3.TPCID.TPC, mPlane);
      // X position of the purported missing vertex
      // A TP for the missing 2D vertex
      TrajPoint vtp;
      vtp.Pos[0] = vx3.Wire;
      vtp.Pos[1] = tcc.detprop->ConvertXToTicks(vx3.X, mPlane, vx3.TPCID.TPC, vx3.TPCID.Cryostat) * tcc.unitsPerTick;
      if(prt) mf::LogVerbatim("TC")<<"CI3DV 3V"<<vx3.ID<<" Pos "<<mPlane<<":"<<PrintPos(slc, vtp.Pos);
      std::vector<int> tjIDs;
      std::vector<unsigned short> tjPts;
      for(auto& tj : slc.tjs) {
        if(tj.CTP != mCTP) continue;
        if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
        if(tj.Pts.size() < 6) continue;
        if(tj.AlgMod[kComp3DVx]) continue;
        float doca = maxdoca;
        // find the closest distance between the vertex and the trajectory
        unsigned short closePt = 0;
        TrajPointTrajDOCA(slc, vtp, tj, closePt, doca);
        if(closePt > tj.EndPt[1]) continue;
        // try to improve the location of the vertex by looking for a distinctive feature on the
        // trajectory, e.g. high multiplicity hits or larger than normal charge
        if(RefineVtxPosition(slc, tj, closePt, 3, false)) vtp.Pos = tj.Pts[closePt].Pos;
        if(prt) mf::LogVerbatim("TC")<<"CI3DV 3V"<<vx3.ID<<" candidate  T"<<tj.ID<<" vtx pos "<<PrintPos(slc, vtp.Pos)<<" doca "<<doca<<" closePt "<<closePt;
        tjIDs.push_back(tj.ID);
        tjPts.push_back(closePt);
      } // itj
      if(tjIDs.empty()) continue;
      // compare the length of the Tjs used to make the vertex with the length of the
      // Tj that we want to split. Don't allow a vertex using very short Tjs to split a long
      // Tj in the 3rd plane
      auto vxtjs = GetAssns(slc, "3V", vx3.ID, "T");
      unsigned short maxPts = 0;
      unsigned short minPts = USHRT_MAX;
      for(auto tjid : vxtjs) {
        auto& tj = slc.tjs[tjid - 1];
        unsigned short npwc = NumPtsWithCharge(slc, tj, false);
        if(npwc > maxPts) maxPts = npwc;
        if(npwc < minPts) minPts = npwc;
      } // tjid
      // skip this operation if any of the Tjs in the split list are > 3 * maxPts
      maxPts *= 3;
      bool skipit = false;
      for(auto tjid : tjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(NumPtsWithCharge(slc, tj, false) > maxPts) skipit = true;
      } // tjid
      if(prt) mf::LogVerbatim("TC")<<"  maxPts "<<maxPts<<" skipit? "<<skipit<<"  minPts "<<minPts;
      if(skipit) continue;
      // 2D vertex
      VtxStore aVtx;
      unsigned short newVtxIndx = slc.vtxs.size();
      aVtx.ID = newVtxIndx + 1;
      aVtx.CTP = mCTP;
      aVtx.Topo = 3;
      aVtx.NTraj = 0;
      // Give it a bogus pass to indicate it wasn't created while stepping
      aVtx.Pass = 9;
      aVtx.Pos = vtp.Pos;
      // ensure this isn't in a messy region
      aVtx.TjChgFrac = ChgFracNearPos(slc, aVtx.Pos, tjIDs);
      if(prt) mf::LogVerbatim("TC")<<" charge fraction near position "<<aVtx.TjChgFrac<<" cut if < 0.6";
      if(aVtx.TjChgFrac < 0.6) continue;
      if(!StoreVertex(slc, aVtx)) continue;
      // make a reference to the new vertex
      VtxStore& newVtx = slc.vtxs[slc.vtxs.size()-1];
      if(prt) mf::LogVerbatim("TC")<<" Stored new 2V"<<newVtx.ID;
      // make a temporary copy so we can nudge it a bit if there is only one Tj
      std::array<float, 2> vpos = aVtx.Pos;
      for(unsigned short ii = 0; ii < tjIDs.size(); ++ii) {
        unsigned short itj = tjIDs[ii] - 1;
        unsigned short closePt = tjPts[ii];
        // determine which end is the closest
        unsigned short end = 1;
        // closest to the beginning?
        if(fabs(closePt - slc.tjs[itj].EndPt[0]) < fabs(closePt - slc.tjs[itj].EndPt[1])) end = 0;
        short dpt = fabs(closePt - slc.tjs[itj].EndPt[end]);
        if(prt) mf::LogVerbatim("TC")<<" dpt "<<dpt<<" to end "<<end;
        if(dpt < 2) {
          // close to an end
          if(slc.tjs[itj].VtxID[end] > 0) {
            // find the distance btw the existing vertex and the end of this tj
            auto& oldTj = slc.tjs[itj];
            auto& oldVx = slc.vtxs[oldTj.VtxID[end] - 1];
            short oldSep = fabs(oldVx.Pos[0] - oldTj.Pts[oldTj.EndPt[end]].Pos[0]);
            if(prt) mf::LogVerbatim("TC")<<" T"<<slc.tjs[itj].ID<<" has vertex 2V"<<slc.tjs[itj].VtxID[end]<<" at end "<<end<<". oldSep "<<oldSep;
            if(dpt < oldSep) {
              MakeVertexObsolete("CI3DV", slc, oldVx, true);
            } else {
              continue;
            }
          } // slc.tjs[itj].VtxID[end] > 0
          slc.tjs[itj].VtxID[end] = slc.vtxs[newVtxIndx].ID;
          ++newVtx.NTraj;
          if(prt) mf::LogVerbatim("TC")<<" attach Traj T"<<slc.tjs[itj].ID<<" at end "<<end;
          slc.tjs[itj].AlgMod[kComp3DVx] = true;
          vpos = slc.tjs[itj].Pts[slc.tjs[itj].EndPt[end]].Pos;
        } else {
          // closePt is not near an end, so split the trajectory
          if(SplitTraj(slc, itj, closePt, newVtxIndx, prt)) {
            if(prt) mf::LogVerbatim("TC")<<" SplitTraj success 2V"<<slc.vtxs[newVtxIndx].ID<<" at closePt "<<closePt;
            // successfully split the Tj
            newVtx.NTraj += 2;
          } else {
            // split failed. Give up
            if(prt) mf::LogVerbatim("TC")<<" SplitTraj failed";
            newVtx.NTraj = 0;
            break;
          }
          // Update the PDGCode for the chopped trajectory
          SetPDGCode(slc, itj);
          // and for the new trajectory
          SetPDGCode(slc, slc.tjs.size()-1);
        } // closePt is not near an end, so split the trajectory
        slc.tjs[itj].AlgMod[kComp3DVx] = true;
        unsigned short newtj = slc.tjs.size() - 1;
        slc.tjs[newtj].AlgMod[kComp3DVx] = true;
      } // ii
      if(newVtx.NTraj == 0) {
        // A failure occurred. Recover
        if(prt) mf::LogVerbatim("TC")<<"  Failed. Recover and delete vertex "<<newVtx.ID;
        MakeVertexObsolete("CI3DV", slc, newVtx, true);
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
        AttachAnyTrajToVertex(slc, newVtx.ID - 1, prt);
        SetVx2Score(slc);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<" Success: new 2V"<<newVtx.ID<<" at "<<(int)newVtx.Pos[0]<<":"<<(int)newVtx.Pos[1]/tcc.unitsPerTick;
          myprt<<" points to 3V"<<vx3.ID;
          myprt<<" TjIDs:";
          for(auto& tjID : tjIDs) myprt<<" T"<<std::to_string(tjID);
        } // prt
      } // success
      ++ivx3;
    } // vx3

  } // CompleteIncomplete3DVertices

  ////////////////////////////////////////////////
  bool RefineVtxPosition(TCSlice& slc, const Trajectory& tj, unsigned short& nearPt, short nPtsToChk, bool prt)
  {
    // The tj has been slated to be split somewhere near point nearPt. This function will move
    // the near point a bit to the most likely point of a vertex

    float maxChg = tj.Pts[nearPt].Chg;
    short maxChgPt = nearPt;
    unsigned short fromPt = tj.EndPt[0];
    short spt = (short)nearPt - (short)nPtsToChk;
    if(spt > (short)fromPt) fromPt = nearPt - nPtsToChk;
    unsigned short toPt = nearPt + nPtsToChk;
    if(toPt > tj.EndPt[1]) toPt = tj.EndPt[1];

    for(short ipt = fromPt; ipt <= toPt; ++ipt) {
      if(ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) continue;
      auto& tp = tj.Pts[ipt];
      if(tp.Chg > maxChg) {
        maxChg = tp.Chg;
        maxChgPt = ipt;
      }
      if(prt) mf::LogVerbatim("TC")<<"RVP: ipt "<<ipt<<" Pos "<<tp.CTP<<":"<<PrintPos(slc, tp.Pos)<<" chg "<<(int)tp.Chg<<" nhits "<<tp.Hits.size();
    } // ipt
    if(nearPt == maxChgPt) return false;
    nearPt = maxChgPt;
    return true;
  } //RefineVtxPosition

  ////////////////////////////////////////////////
  bool MakeVertexObsolete(std::string fcnLabel, TCSlice& slc, VtxStore& vx2, bool forceKill)
  {
    // Makes a 2D vertex obsolete

    // check for a high-score 3D vertex
    bool hasHighScoreVx3 = (vx2.Vx3ID > 0);
    if(hasHighScoreVx3 && !forceKill && slc.vtx3s[vx2.Vx3ID - 1].Score >= tcc.vtx2DCuts[7]) return false;

    if(tcc.dbg2V || tcc.dbg3V) {
      mf::LogVerbatim("TC")<<fcnLabel<<" MVO: killing 2V"<<vx2.ID;
    }

    // Kill it
    int vx2id = vx2.ID;
    if(vx2.Vx3ID > 0) {
      auto& vx3 = slc.vtx3s[vx2.Vx3ID - 1];
      for(auto& v3v2id : vx3.Vx2ID) if(v3v2id == vx2.ID) v3v2id = 0;
    }
    vx2.ID = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
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
            auto& ovx2 = slc.vtxs[tj.VtxID[oend] - 1];
            if(!ovx2.Stat[kHiVx3Score]) tj.AlgMod[kTjHiVx3Score] = false;
          } // vertex at the other end
        } // tj.AlgMod[kTjHiVx3Score]
      } // end
    } // tj

    if(!hasHighScoreVx3) return true;

    // update the affected 3D vertex
    Vtx3Store& vx3 = slc.vtx3s[vx2.Vx3ID - 1];
    // make the 3D vertex incomplete
    geo::PlaneID planeID = DecodeCTP(vx2.CTP);
    unsigned short plane = planeID.Plane;
    if(vx3.Vx2ID[plane] != vx2id) return true;
    vx3.Vx2ID[plane] = 0;
    vx3.Wire = vx2.Pos[0];
    // Ensure that there are at least two 2D vertices left
    unsigned short n2D = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(vx3.Vx2ID[plane] > 0) ++n2D;

    if(n2D > 1) {
      // 3D vertex is incomplete
      // correct the score
      SetVx3Score(slc, vx3);
      return true;
    }

    // 3D vertex is obsolete
    // Detach all remaining 2D vertices from the 3D vertex
    for(auto& vx2 : slc.vtxs) {
      if(vx2.ID == 0) continue;
      if(vx2.Vx3ID == vx3.ID) vx2.Vx3ID = 0;
    } // vx2
    for(auto& pfp : slc.pfps) {
      for(unsigned short end = 0; end < 2; ++end) if(pfp.Vx3ID[end] == vx3.ID) pfp.Vx3ID[end] = 0;
    } // pfp
    vx3.ID = 0;
    return true;

  } // MakeVertexObsolete

  ////////////////////////////////////////////////
  bool MakeVertexObsolete(TCSlice& slc, Vtx3Store& vx3)
  {
    // Deletes a 3D vertex and 2D vertices in all planes
    // The 2D and 3D vertices are NOT killed if forceKill is false and the 3D vertex
    // has a high score
    if(vx3.ID <= 0) return true;
    if(vx3.ID > int(slc.vtx3s.size())) return false;

    for(auto vx2id : vx3.Vx2ID) {
      if(vx2id == 0 || vx2id > (int)slc.vtxs.size()) continue;
      auto& vx2 = slc.vtxs[vx2id - 1];
      MakeVertexObsolete("MVO3", slc, vx2, true);
    }
    vx3.ID = 0;
    return true;
  } // MakeVertexObsolete

  //////////////////////////////////////////
  std::vector<int> GetVtxTjIDs(const TCSlice& slc, const VtxStore& vx2)
  {
    // returns a list of trajectory IDs that are attached to vx2
    std::vector<int> tmp;
    if(vx2.ID == 0) return tmp;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx2.CTP) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == vx2.ID) tmp.push_back(tj.ID);
      } // end
    } // tj
    return tmp;
  } // GetVtxTjIDs


  //////////////////////////////////////////
  std::vector<int> GetVtxTjIDs(const TCSlice& slc, const Vtx3Store& vx3, float& score)
  {
    // returns a list of Tjs in all planes that are attached to vx3
    std::vector<int> tmp;
    if(vx3.ID == 0) return tmp;
    float nvx2 = 0;
    score = 0;
    for(auto& vx2 : slc.vtxs) {
      if(vx2.ID == 0) continue;
      if(vx2.Vx3ID != vx3.ID) continue;
      auto vtxTjID2 = GetVtxTjIDs(slc, vx2);
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
  void PosInPlane(const TCSlice& slc, const Vtx3Store& vx3, unsigned short plane, Point2_t& pos)
  {
    // returns the 2D position of the vertex in the plane
    pos[0] = tcc.geom->WireCoordinate(vx3.Y, vx3.Z, plane, vx3.TPCID.TPC, vx3.TPCID.Cryostat);
    pos[1] = tcc.detprop->ConvertXToTicks(vx3.X, plane, vx3.TPCID.TPC, vx3.TPCID.Cryostat) * tcc.unitsPerTick;

  } // PosInPlane

  /////////////////////////////////////////
  unsigned short IsCloseToVertex(const TCSlice& slc, const VtxStore& inVx2)
  {
    // Returns the ID of a 2D vertex having the minimum pull < user-specified cut

    float minPull = tcc.vtx2DCuts[3];
    unsigned short imBest = 0;
    for(auto& vx2 : slc.vtxs) {
      if(vx2.CTP != inVx2.CTP) continue;
      if(vx2.ID <= 0) continue;
      float pull = VertexVertexPull(slc, inVx2, vx2);
      if(pull < minPull) {
        minPull = pull;
        imBest = vx2.ID;
      }
    } // vx2
    return imBest;
  } // IsCloseToVertex

  /////////////////////////////////////////
  unsigned short IsCloseToVertex(const TCSlice& slc, const Vtx3Store& vx3)
  {
    // Returns the ID of a 3D vertex having the minimum pull < user-specified cut

    float minPull = tcc.vtx3DCuts[1];
    unsigned short imBest = 0;
    for(auto& oldvx3 : slc.vtx3s) {
      if(oldvx3.ID == 0) continue;
      if(std::abs(oldvx3.X - vx3.X) > tcc.vtx3DCuts[0]) continue;
      float pull = VertexVertexPull(slc, vx3, oldvx3);
      if(pull < minPull) {
        minPull = pull;
        imBest = oldvx3.ID;
      }
    } // oldvx3
    return imBest;

  } // IsCloseToVertex

} // namespace
