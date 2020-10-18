#include "larreco/RecoAlg/TCAlg/PostStepUtils.h"

#include <limits.h>                                        // for USHRT_MAX
#include <stdlib.h>                                        // for abs, size_t
#include <cmath>                                           // for sqrt, atan
#include <algorithm>                                       // for find, max
#include <array>                                           // for array, arr...
#include <bitset>                                          // for bitset<>::...
#include <iomanip>                                         // for operator<<
#include <iostream>                                        // for cout
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"   // for TDCtick_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"  // for PlaneID
#include "lardataobj/RecoBase/Hit.h"                       // for Hit
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"             // for DebugStuff
#include "larreco/RecoAlg/TCAlg/TCVertex.h"                // for tcc, evt
#include "larreco/RecoAlg/TCAlg/Utils.h"                   // for SetEndPoints
#include <math.h>                                          // for abs, nearb...
#include <numeric>                                         // for iota
#include <string>                                          // for basic_string
#include <utility>                                         // for pair
#include <vector>                                          // for vector

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace tca {

  ////////////////////////////////////////////////
  void CheckTraj(TCSlice& slc, Trajectory& tj)
  {
    // Check the quality of the trajectory and possibly trim it or flag it for deletion

    if(!tj.IsGood) return;

    // ensure that the end points are defined
    SetEndPoints(tj);
    if(tj.EndPt[0] == tj.EndPt[1]) return;

    if(tj.Strategy[kStiffEl]) {
      CheckStiffEl(slc, tj);
      return;
    }

    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"inside CheckTraj with NumPtsWithCharge = "<<NumPtsWithCharge(slc, tj, false);
    }

    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    if(npwc < tcc.minPts[tj.Pass]) {
      tj.IsGood = false;
      return;
    }

    // Find the average number of used hits/TP and try to use the
    // unused hits
//    UseUnusedHits(slc, tj);

    // correct the start TP of short Tjs
    if(npwc < 5) {
      auto& tp0 = tj.Pts[tj.EndPt[0]];
      auto& tp1 = tj.Pts[tj.EndPt[1]];
      tp0.Dir = tp1.Dir;
      float dw = tp0.Pos[0] - tp1.Pos[0];
      if(tp0.Dir[0] != 0) tp0.Pos[1] = tp1.Pos[1] + dw * tp0.Dir[1] / tp0.Dir[0];
      tp0.Delta = PointTrajDOCA(slc, tp0.HitPos[0], tp0.HitPos[1], tp0);
      tp0.Ang = tp1.Ang;
      tp0.AngErr = tp1.AngErr;
      tp0.AveChg = tp1.AveChg;
    } // npwc < 4

    // Look for a charge asymmetry between points on both sides of a high-
    // charge point and trim points in the vicinity
    ChkChgAsymmetry(slc, tj, tcc.dbgStp);

    // flag this tj as a junk Tj (even though it wasn't created in FindJunkTraj).
    // Drop it and let FindJunkTraj do it's job
    TagJunkTj(slc, tj, tcc.dbgStp);
    if(tj.AlgMod[kJunkTj]) {
      tj.IsGood = false;
      return;
    }

    tj.MCSMom = MCSMom(slc, tj);

    // See if the points at the stopping end can be included in the Tj
    ChkStopEnd1(slc, tj, tcc.dbgStp);

    // remove any points at the end that don't have charge
    tj.Pts.resize(tj.EndPt[1] + 1);

    // Ensure that a hit only appears once in the TJ
    if(HasDuplicateHits(slc, tj, tcc.dbgStp)) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" HasDuplicateHits ";
      tj.IsGood = false;
      return;
    }

    // See if this is a ghost trajectory
    if(IsGhost(slc, tj)) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" CT: Ghost trajectory - trimmed hits ";
      if(!tj.IsGood) return;
    }

    if(tj.AlgMod[kJunkTj]) return;

    // checks are different for Very Large Angle trajectories
    bool isVLA = (tj.Pts[tj.EndPt[1]].AngleCode == 2);

    tj.Pts.resize(tj.EndPt[1] + 1);

    // Fill in any gaps with hits that were skipped, most likely delta rays on muon tracks
    if(!isVLA) FillGaps(slc, tj);

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" CheckTraj MCSMom "<<tj.MCSMom<<" isVLA? "<<isVLA<<" NumPtsWithCharge "<<NumPtsWithCharge(slc, tj, false)<<" Min Req'd "<<tcc.minPts[tj.Pass];

    // Trim the end points until the TJ meets the quality cuts
    TrimEndPts("CT", slc, tj, tcc.qualityCuts, tcc.dbgStp);
    if(tj.AlgMod[kKilled]) {
      tj.IsGood = false;
      return;
    }

    // Check for a Bragg peak at both ends. This may be used by FixBegin.
    ChkStop(slc, tj);

    // Update the trajectory parameters at the beginning of the trajectory
    ChkBegin(slc, tj);

    // Look for a high charge plateau at the end (kLEPhys)
    TrimHiChgEndPts(slc, tj);
    // Look for high multiplicity hits at the end (kLEPhys)
    TrimHiMultEndPts(slc, tj);
    // Old algorithm to be removed
    ChkHiMultEndHits(slc, tj);
    // Check for a Bragg peak at both ends one last time
    ChkStop(slc, tj);
    // Check the ChiDOF of the last points and re-fit if necessary
    ChkEndPtFit(slc, tj);

    // check for garbage Tjs and ignore the standard quality cuts
    if(tcc.useAlg[kLEPhys] && npwc < 5 && !isVLA) {
      float tothits = 0;
      for(auto& tp : tj.Pts) tothits += tp.Hits.size();
      float hitsPerTP = tothits / (float)tj.Pts.size();
      if(tj.Pts[tj.EndPt[1]].AngleCode == 0) {
        if(hitsPerTP > 1.5) tj.IsGood = false;
      } else {
        if(hitsPerTP > 2.0) tj.IsGood = false;
      }
      if(tcc.dbgStp) mf::LogVerbatim("TC") << "CheckTraj: hitsPerTP " << hitsPerTP << " Good? " << tj.IsGood;
      return;
    } // npwc < 6

    // final quality check
    float npts = std::abs(tj.EndPt[1] - tj.EndPt[0]);
    npwc = NumPtsWithCharge(slc, tj, false);
    float frac = (float)npwc / npts;
    tj.IsGood = (frac >= tcc.qualityCuts[0]);
    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"CheckTraj: fraction of points with charge "
        << std::setprecision(2) << frac <<" good traj? "<<tj.IsGood
        << " chk npwc " << npwc << " npts " << npts;
    }
    if(!tj.IsGood || !slc.isValid) return;

  } // CheckTraj

  ////////////////////////////////////////////////
  void ChkStopEnd1(TCSlice& slc, Trajectory& tj, bool prt)
  {
    // Analyze the end of the Tj after crawling has stopped to see if any of the points
    // should be used

    // don't use this old algorithm
    if(tcc.useAlg[kLEPhys]) return;
    if(tj.EndFlag[1][kEndKink]) return;
    if(!tcc.useAlg[kChkStopEP]) return;
    if(tj.AlgMod[kJunkTj]) return;
    if(tj.Strategy[kStiffEl]) return;

    unsigned short endPt = tj.EndPt[1];
    // ignore VLA Tjs
    if(tj.Pts[endPt].AngleCode > 1) return;
    // don't get too carried away with this
    if(tj.Pts.size() - endPt > 10) return;

    // Get a list of hits a few wires beyond the last point on the Tj
    geo::PlaneID planeID = DecodeCTP(tj.CTP);
    unsigned short plane = planeID.Plane;

    unsigned short lastPt = tj.EndPt[1];
    auto& lastTP = tj.Pts[lastPt];
    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"CSEP: checking "<<tj.ID<<" end at "<<PrintPos(slc, lastTP.Pos);
    }
    // make a copy that we can move
    auto ltp = lastTP;

    double stepSize = std::abs(1/ltp.Dir[0]);
    std::array<int, 2> wireWindow;
    std::array<float, 2> timeWindow;
    std::vector<int> closeHits;
    bool isClean = true;
    for(unsigned short step = 0; step < 10; ++step) {
      for(unsigned short iwt = 0; iwt < 2; ++iwt) ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;
      int wire = std::nearbyint(ltp.Pos[0]);
      wireWindow[0] = wire;
      wireWindow[1] = wire;
      timeWindow[0] = ltp.Pos[1] - 5;
      timeWindow[1] = ltp.Pos[1] + 5;
      bool hitsNear;
      auto tmp = FindCloseHits(slc, wireWindow, timeWindow, plane, kAllHits, true, hitsNear);
      // add close hits that are not associated with this tj
      for(auto iht : tmp) if(slc.slHits[iht].InTraj != tj.ID) closeHits.push_back(iht);
      float nWiresPast = 0;
      // Check beyond the end of the trajectory to see if there are hits there
      if(ltp.Dir[0] > 0) {
        // stepping +
        nWiresPast = ltp.Pos[0] - lastTP.Pos[0];
      }  else {
        // stepping -
        nWiresPast = lastTP.Pos[0] - ltp.Pos[0];
      }
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Found "<<tmp.size()<<" hits near pos "<<PrintPos(slc, ltp.Pos)<<" nWiresPast "<<nWiresPast;
      if(nWiresPast > 0.5) {
        if(!tmp.empty()) isClean = false;
        if(nWiresPast > 1.5) break;
      } // nWiresPast > 0.5
    } // step

    // count the number of available hits
    unsigned short nAvailable = 0;
    for(auto iht : closeHits) if(slc.slHits[iht].InTraj == 0) ++nAvailable;

    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<"closeHits";
      for(auto iht : closeHits) myprt<<" "<<PrintHit(slc.slHits[iht]);
      myprt<<" nAvailable "<<nAvailable;
      myprt<<" isClean "<<isClean;
    } // prt

    if(!isClean || nAvailable != closeHits.size()) return;

    unsigned short originalEndPt = tj.EndPt[1] + 1;
    // looks clean so use all the hits
    for(unsigned short ipt = originalEndPt; ipt <= lastPt; ++ipt) {
      auto& tp = tj.Pts[ipt];
      bool hitsAdded = false;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        // This shouldn't happen but check anyway
        if(slc.slHits[tp.Hits[ii]].InTraj != 0) continue;
        tp.UseHit[ii] = true;
        slc.slHits[tp.Hits[ii]].InTraj = tj.ID;
        hitsAdded = true;
      } // ii
      if(hitsAdded) DefineHitPos(slc, tp);
    } // ipt
    tj.AlgMod[kChkStopEP] = true;
    SetEndPoints(tj);
    // Re-fitting the end might be a good idea but it's probably not necessary. The
    // values of Delta should have already been filled

    if(!tcc.useAlg[kLEPhys]) {
      // require a Bragg peak
      ChkStop(slc, tj);
      if(!tj.EndFlag[1][kEndBragg]) {
        // restore the original
        for(unsigned short ipt = originalEndPt; ipt <= lastPt; ++ipt) UnsetUsedHits(slc, tj.Pts[ipt]);
        SetEndPoints(tj);
      } // no Bragg Peak
    } // !tcc.useAlg[kLEPhys]

    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<"CSEP: Leaving with T"<<tj.ID<<" extent "<<PrintPos(slc, tj.Pts[tj.EndPt[0]]);
      myprt<<" to "<<PrintPos(slc, tj.Pts[tj.EndPt[1]]);
    }

    UpdateTjChgProperties("CSEP", slc, tj, prt);

  } // ChkStopEnd1

  ////////////////////////////////////////////////
  void ChkBegin(TCSlice& slc, Trajectory& tj)
  {
    // Check the parameters at the start of the trajectory. The first
    // points may not belong to this trajectory since they were added when there was
    // little information. This information may be updated later if ReversePropagate is used

    if(!tcc.useAlg[kFixBegin]) return;
    if(tj.AlgMod[kJunkTj]) return;

    // don't do anything if this tj has been modified by ReversePropagate
    if(tj.AlgMod[kRvPrp]) return;

    // don't bother with really short tjs
    if(tj.Pts.size() < 3) return;

    unsigned short atPt = tj.EndPt[1];
    unsigned short maxPtsFit = 0;
    unsigned short firstGoodChgPullPt = USHRT_MAX;
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      if(tp.AveChg > 0 && firstGoodChgPullPt == USHRT_MAX) {
        if(std::abs(tp.ChgPull) < tcc.chargeCuts[0]) firstGoodChgPullPt = ipt;
      } // find the first good charge pull point
      if(tp.NTPsFit > maxPtsFit) {
        maxPtsFit = tp.NTPsFit;
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
    unsigned short nPtsLeft = NumPtsWithCharge(slc, tj, false) - firstPtFit;
    if(needsRevProp) {
      needsRevProp = (nPtsLeft > 5);
    }
    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<"CB: firstPtFit "<<firstPtFit<<" at "<<PrintPos(slc, tj.Pts[firstPtFit].Pos);
      myprt<<" atPt "<<PrintPos(slc, tj.Pts[atPt].Pos);
      myprt<<" nPts with charge "<<nPtsLeft;
      myprt<<" firstGoodChgPullPt "<<firstGoodChgPullPt;
      if(firstGoodChgPullPt != USHRT_MAX) myprt<<" at "<<PrintPos(slc,tj.Pts[firstGoodChgPullPt]);
      myprt<<" needsRevProp? "<<needsRevProp;
    }

    if(!needsRevProp && firstGoodChgPullPt == tj.EndPt[0]) {
      // check one wire on the other side of EndPt[0] to see if there are hits that are available which could
      // be picked up by reverse propagation
      TrajPoint tp = tj.Pts[0];
      tp.Hits.clear();
      tp.UseHit.reset();
      // Move the TP "backwards"
      double stepSize = tcc.VLAStepSize;
      if(tp.AngleCode < 2) stepSize = std::abs(1/tp.Dir[0]);
      tp.Pos[0] -= tp.Dir[0] * stepSize * tj.StepDir;
      tp.Pos[1] -= tp.Dir[1] * stepSize * tj.StepDir;
      // launch RevProp if this wire is dead
      unsigned int wire = std::nearbyint(tp.Pos[0]);
      unsigned short plane = DecodeCTP(tp.CTP).Plane;
      needsRevProp = (wire < slc.nWires[plane] && !evt.goodWire[plane][wire]);
      if(tcc.dbgStp && needsRevProp) mf::LogVerbatim("TC")<<"CB: Previous wire "<<wire<<" is dead. Call ReversePropagate";
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"CB: Look for  hits near "<<PrintPos(slc, tp);
      if(!needsRevProp && firstGoodChgPullPt != USHRT_MAX) {
        // check for hits on a not-dead wire
        // BB Do this more carefully
        float maxDelta = 2 * tp.DeltaRMS;
        if(FindCloseHits(slc, tp, maxDelta, kAllHits) && !tp.Hits.empty()) {
          // count used and unused hits
          unsigned short nused = 0;
          for(auto iht : tp.Hits) if(slc.slHits[iht].InTraj > 0) ++nused;
          if(nused == 0) {
            needsRevProp = true;
            if(tcc.dbgStp) {
              mf::LogVerbatim("TC")<<"CB: Found "<<tp.Hits.size()-nused<<" close unused hits found near EndPt[0] "<<PrintPos(slc, tp)<<". Call ReversePropagate";
            } // tcc.dbgStp
          } // nused = 0
        } // Close hits exist
      } // !needsRevProp
    } // !needsRevProp

    if(tcc.dbgStp) {
      mf::LogVerbatim("TC") << "CB: maxPtsFit "<<maxPtsFit<<" at "<<PrintPos(slc, tj.Pts[atPt])
          <<" firstPtFit "<<firstPtFit<<" Needs ReversePropagate? "<<needsRevProp;
    }

    if(tcc.useAlg[kFTBRvProp] && needsRevProp) {
      // lop off the points before firstPtFit and reverse propagate
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"  clobber TPs "<<PrintPos(slc, tj.Pts[0])<<" to "<<PrintPos(slc, tj.Pts[firstPtFit])<<". Call TrimEndPts then ReversePropagate ";
      // first save the first TP on this trajectory. We will try to re-use it if
      // it isn't used during reverse propagation
      seeds.push_back(tj.Pts[0]);
      for(unsigned short ipt = 0; ipt <= firstPtFit; ++ipt) UnsetUsedHits(slc, tj.Pts[ipt]);
      SetEndPoints(tj);
      tj.AlgMod[kFTBRvProp] = true;
      // Check for quality and trim if necessary before reverse propagation
      TrimEndPts("RPi", slc, tj, tcc.qualityCuts, tcc.dbgStp);
      if(tj.AlgMod[kKilled]) {
        tj.IsGood = false;
        return;
      }
      ReversePropagate(slc, tj);
      ChkStopEnd1(slc, tj, tcc.dbgStp);
    }
    if(firstGoodChgPullPt != USHRT_MAX && firstGoodChgPullPt > atPt) atPt = firstGoodChgPullPt;
    // Update the fit parameters of the first points if no reverse propagation was done
    if(!tj.AlgMod[kRvPrp]) FixBegin(slc, tj, atPt);
  } // ChkBegin

  ////////////////////////////////////////////////
  void FixBegin(TCSlice& slc, Trajectory& tj, unsigned short atPt)
  {
    // Update the parameters at the beginning of the trajectory starting at point atPt

    if(!tcc.useAlg[kFixBegin]) return;
    // ignore short trajectories
    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    if(npwc < 6) return;
    // ignore shower-like trajectories
    if(tj.PDGCode == 11) return;
    // ignore junk trajectories
    if(tj.AlgMod[kJunkTj]) return;
    unsigned short firstPt = tj.EndPt[0];

    if(atPt == tj.EndPt[0]) return;

    // Default is to use DeltaRMS of the last point on the Tj
    float maxDelta = 4 * tj.Pts[tj.EndPt[1]].DeltaRMS;

    // Find the max DeltaRMS of points from atPt to EndPt[1]
    float maxDeltaRMS = 0;
    for(unsigned short ipt = atPt; ipt <= tj.EndPt[1]; ++ipt) {
      if(tj.Pts[ipt].DeltaRMS > maxDeltaRMS) maxDeltaRMS = tj.Pts[ipt].DeltaRMS;
    } // ipt
    maxDelta = 3 * maxDeltaRMS;

    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"FB: atPt "<<atPt<<" firstPt "<<firstPt<<" Stops at end 0? "<<PrintEndFlag(tj, 0)<<" start vertex "<<tj.VtxID[0]<<" maxDelta "<<maxDelta;
    }

    // update the trajectory for all the points up to atPt
    // assume that we will use all of these points
    bool maskedPts = false;
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
      tp.Delta = PointTrajDOCA(slc, tp.HitPos[0], tp.HitPos[1], tp);
      tp.DeltaRMS = tj.Pts[atPt].DeltaRMS;
      tp.NTPsFit = tj.Pts[atPt].NTPsFit;
      tp.FitChi = tj.Pts[atPt].FitChi;
      tp.AveChg = tj.Pts[atPt].AveChg;
      tp.ChgPull = (tj.Pts[ipt].Chg / tj.AveChg - 1) / tj.ChgRMS;
      bool badChg = (std::abs(tp.ChgPull) > tcc.chargeCuts[0]);
      bool maskThisPt = (tp.Delta > maxDelta || badChg);
      if(maskThisPt) {
        UnsetUsedHits(slc, tp);
        maskedPts = true;
      }
      if(tcc.dbgStp) {
        mf::LogVerbatim myprt("TC");
        myprt<<" Point "<<PrintPos(slc, tj.Pts[ipt].Pos)<<" Delta "<<tj.Pts[ipt].Delta<<" ChgPull "<<tj.Pts[ipt].ChgPull<<" maskThisPt? "<<maskThisPt;
      }
      if(ipt == 0) break;
    } // ii
    if(maskedPts) SetEndPoints(tj);
    tj.AlgMod[kFixBegin] = true;
  } // FixBegin

  ////////////////////////////////////////////////
  bool IsGhost(TCSlice& slc, Trajectory& tj)
  {
    // Sees if trajectory tj shares many hits with another trajectory and if so merges them.

    if(!tcc.useAlg[kUseGhostHits]) return false;
    // ensure that tj is not a saved trajectory
    if(tj.ID > 0) return true;
    // or an already killed trajectory
    if(tj.AlgMod[kKilled]) return true;
    if(tj.Pts.size() < 3) return false;
    if(tj.Strategy[kStiffEl]) return false;

    // vectors of traj IDs, and the occurrence count
    std::vector<int> tID;
    std::vector<unsigned short> tCnt;

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
        if(slc.slHits[iht].InTraj > 0 && (unsigned int)slc.slHits[iht].InTraj <= slc.tjs.size()) {
          int tjid = slc.slHits[iht].InTraj;
          unsigned short indx;
          for(indx = 0; indx < tID.size(); ++indx) if(tID[indx] == tjid) break;
          if(indx == tID.size()) {
            tID.push_back(tjid);
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

    if(tcc.dbgStp) {
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
    Trajectory& oTj = slc.tjs[oldTjIndex];
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

    int nwires = wire1 - wire0 + 1;
    std::vector<float> oTjPos1(nwires, -1);
    unsigned short nMissedWires = 0;
    for(unsigned short ipt = oTj.EndPt[0]; ipt <= oTj.EndPt[1]; ++ipt) {
      if(oTj.Pts[ipt].Chg == 0) continue;
      int wire = std::nearbyint(oTj.Pts[ipt].Pos[0]);
      int indx = wire - wire0;
      if(indx < 0 || indx > nwires - 1) continue;
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
          if(slc.slHits[iht].InTraj ==  oldTjID) HitInoTj = true;
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

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Number of missed wires in oTj gaps "<<nMissedWires<<" Number of ghost hits in these gaps "<<ngh<<" nghPlus "<<nghPlus<<" cut "<<0.2 * nMissedWires;

    if(ngh < 0.2 * nMissedWires) return false;
    if(firstPtInoTj > lastPtInoTj) return false;

    // require all of the tj TPs to be on either the + or - side of the oTj trajectory
    if(!(nghPlus > 0.8 * ngh || nghPlus < 0.2 * ngh) ) return false;

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Trajectory is a ghost of "<<oldTjID<<" first point in oTj "<<firstPtInoTj<<" last point "<<lastPtInoTj;

    // unset all of the shared hits
    for(unsigned short ipt = firstPtInoTj; ipt <= lastPtInoTj; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      UnsetUsedHits(slc, tj.Pts[ipt]);
      if(tcc.dbgStp) PrintTrajectory("IG", slc, tj, ipt);
    }
    // see how many points are left at the end
    ngh = 0;
    for(unsigned short ipt = lastPtInoTj; ipt <= tj.EndPt[1]; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ngh;
    } // ipt
    // clobber those too?
    if(ngh > 0 && ngh < tcc.minPts[tj.Pass]) {
      for(unsigned short ipt = lastPtInoTj; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg > 0) UnsetUsedHits(slc, tj.Pts[ipt]);
      } // ipt
    }
    SetEndPoints(tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    slc.tjs[oldTjIndex].AlgMod[kUseGhostHits] = true;
    TrimEndPts("IG", slc, tj, tcc.qualityCuts, tcc.dbgStp);
    if(tj.AlgMod[kKilled]) {
      tj.IsGood = false;
      if(tcc.dbgStp)  mf::LogVerbatim("TC")<<" Failed quality cuts";
      return true;
    }
    tj.MCSMom = MCSMom(slc, tj);
    if(tcc.dbgStp)  mf::LogVerbatim("TC")<<" New tj size "<<tj.Pts.size();
    return true;
  } // IsGhost

  ////////////////////////////////////////////////
  bool IsGhost(TCSlice& slc, std::vector<unsigned int>& tHits)
  {
    // Called by FindJunkTraj to see if the passed hits are close to an existing
    // trajectory and if so, they will be used in that other trajectory

    if(!tcc.useAlg[kUseGhostHits]) return false;

    if(tHits.size() < 2) return false;

    bool prt = (tcc.dbgStp || tcc.dbgAlg[kUseGhostHits]);

    // find all nearby hits
    std::vector<unsigned int> hitsInMuliplet, nearbyHits;
    for(auto iht : tHits) {
      GetHitMultiplet(slc, iht, hitsInMuliplet, false);
      // prevent double counting
      for(auto mht : hitsInMuliplet) {
        if(std::find(nearbyHits.begin(), nearbyHits.end(), mht) == nearbyHits.end()) {
          nearbyHits.push_back(mht);
        }
      } // mht
    } // iht

    // vectors of traj IDs, and the occurrence count
    std::vector<unsigned int> tID, tCnt;
    for(auto iht : nearbyHits) {
      if(slc.slHits[iht].InTraj <= 0) continue;
      unsigned int tid = slc.slHits[iht].InTraj;
      unsigned short indx = 0;
      for(indx = 0; indx < tID.size(); ++indx) if(tID[indx] == tid) break;
      if(indx == tID.size()) {
        tID.push_back(tid);
        tCnt.push_back(1);
      }  else {
        ++tCnt[indx];
      }
    } // iht
    if(tCnt.empty()) return false;

    // Call it a ghost if > 50% of the hits are used by another trajectory
    unsigned short tCut = 0.5 * tHits.size();
    int tid = INT_MAX;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"IsGhost tHits size "<<tHits.size()<<" cut fraction "<<tCut<<" tID_tCnt";
      for(unsigned short ii = 0; ii < tCnt.size(); ++ii) myprt<<" "<<tID[ii]<<"_"<<tCnt[ii];
    } // prt

    for(unsigned short ii = 0; ii < tCnt.size(); ++ii) {
      if(tCnt[ii] > tCut) {
        tid = tID[ii];
        break;
      }
    } // ii
    if(tid > (int)slc.tjs.size()) return false;

    if(prt) mf::LogVerbatim("TC")<<" is ghost of trajectory "<<tid;

    // Use all hits in tHits that are found in itj
    for(auto& tp : slc.tjs[tid - 1].Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(slc.slHits[iht].InTraj != 0) continue;
        for(unsigned short jj = 0; jj < tHits.size(); ++jj) {
          unsigned int tht = tHits[jj];
          if(tht != iht) continue;
          tp.UseHit[ii] = true;
          slc.slHits[iht].InTraj = tid;
          break;
        } // jj
      } // ii
    } // tp
    slc.tjs[tid - 1].AlgMod[kUseGhostHits] = true;
    return true;
  } // IsGhost

  ////////////////////////////////////////////////
  void FillGaps(TCSlice& slc, Trajectory& tj)
  {
    // Fill in any gaps in the trajectory with close hits regardless of charge (well maybe not quite that)

    if(!tcc.useAlg[kFillGaps]) return;
    if(tj.AlgMod[kJunkTj]) return;
    if(tj.ChgRMS <= 0) return;

    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    if(npwc < 8) return;

    // don't consider the last few points since they would be trimmed
    unsigned short toPt = tj.EndPt[1] - 2;
    if(!tj.EndFlag[1][kEndBragg]) {
      // Don't fill gaps (with high-charge hits) near the end. Find the last point near the
      // end that would have normal charge if all the hit were added
      unsigned short cnt = 0;
      for(unsigned short ipt = tj.EndPt[1] - 2; ipt > tj.EndPt[0]; --ipt) {
        auto& tp = tj.Pts[ipt];
        float chg = tp.Chg;
        if(chg == 0) {
          for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
            unsigned int iht = tp.Hits[ii];
            auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
            chg += hit.Integral();
          }
        } // chg == 0
        float chgPull = (chg / tj.AveChg - 1) / tj.ChgRMS;
        if(chgPull < 2) {
          toPt = ipt;
          break;
        }
        ++cnt;
        if(cnt > 20) break;
      } // ipt
    } // !tj.EndFlag[1][kEndBragg]

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"FG: Check Tj "<<tj.ID<<" from "<<PrintPos(slc, tj.Pts[tj.EndPt[0]])<<" to "<<PrintPos(slc, tj.Pts[toPt]);

    // start with the first point that has charge
    short firstPtWithChg = tj.EndPt[0];
    bool first = true;
    float maxDelta = 1;
    // don't let MCSMom suffer too much while filling gaps
    short minMCSMom = 0.7 * tj.MCSMom;
    while(firstPtWithChg < toPt) {
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
      // Found a gap. Require at least two consecutive points with charge after the gap
      if(nextPtWithChg < (tj.EndPt[1] - 1) && tj.Pts[nextPtWithChg + 1].Chg == 0) {
        firstPtWithChg = nextPtWithChg;
        continue;
      }
      // Make a bare trajectory point at firstPtWithChg that points to nextPtWithChg
      TrajPoint tp;
      if(!MakeBareTrajPoint(slc, tj.Pts[firstPtWithChg], tj.Pts[nextPtWithChg], tp)) {
        tj.IsGood = false;
        return;
      }
      // Find the maximum delta between hits and the trajectory Pos for all
      // hits on this trajectory
      if(first) {
        maxDelta = 2.5 * MaxHitDelta(slc, tj);
        first = false;
      } // first
      // define a loose charge cut using the average charge at the first point with charge
      float maxChg = tj.Pts[firstPtWithChg].AveChg * (1 + 2 * tcc.chargeCuts[0] * tj.ChgRMS);
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
          mf::LogVerbatim("TC")<<"FillGaps coding error: firstPtWithChg "<<firstPtWithChg<<" mpt "<<mpt<<" nextPtWithChg "<<nextPtWithChg;
          slc.isValid = false;
          return;
        }
        bool filled = false;
        float chg = 0;
        for(unsigned short ii = 0; ii < tj.Pts[mpt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[mpt].Hits[ii];
          if(slc.slHits[iht].InTraj > 0) continue;
          auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
          float delta = PointTrajDOCA(slc, iht, tp);
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<" FG: "<<PrintPos(slc,tj.Pts[mpt])<<" hit "<<PrintHit(slc.slHits[iht])<<" delta "<<delta<<" maxDelta "<<maxDelta<<" Chg "<<hit.Integral()<<" maxChg "<<maxChg;
          if(delta > maxDelta) continue;
          tj.Pts[mpt].UseHit[ii] = true;
          slc.slHits[iht].InTraj = tj.ID;
          chg += hit.Integral();
          filled = true;
        } // ii
        if(chg > maxChg || MCSMom(slc, tj) < minMCSMom) {
          // don't use these hits after all
          UnsetUsedHits(slc, tj.Pts[mpt]);
          filled = false;
        }
        if(filled) {
          DefineHitPos(slc, tj.Pts[mpt]);
          tj.AlgMod[kFillGaps] = true;
          if(tcc.dbgStp) {
            PrintTP("FG", slc, mpt, tj.StepDir, tj.Pass, tj.Pts[mpt]);
            mf::LogVerbatim("TC")<<"Check MCSMom "<<MCSMom(slc, tj);
          }
        } // filled
      } // mpt
      firstPtWithChg = nextPtWithChg;
    } // firstPtWithChg
    if(tj.AlgMod[kFillGaps]) tj.MCSMom = MCSMom(slc, tj);
  } // FillGaps

  //////////////////////////////////////////
  void
  TrimHiChgEndPts(TCSlice& slc, Trajectory& tj)
  {
    // look for a significant step increase at the end, indicating that two tracks are
    // overlapping in this view. 
    if(!tcc.useAlg[kLEPhys]) return;
    if(!tcc.useAlg[kTHiQEP]) return;
    if(tj.Strategy[kSlowing]) return;
    if(tj.EndFlag[1][kEndBragg]) return;
    if(tj.EndPt[1] - tj.EndPt[0] < 10) return;

    bool prt = tcc.dbgStp;

    unsigned short cnt = 0;
    float maxAsym = 0.5;
    short atPt = -1;
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      short ipt = tj.EndPt[1] - ii - 1;
      if(ipt < tj.EndPt[0] + 3) break;
      auto& tpp0 = tj.Pts[ipt];
      if(tpp0.Chg <= 0) return;
      auto& tpp1 = tj.Pts[ipt+1];
      if(tpp1.Chg <= 0) return;
      auto& tpm0 = tj.Pts[ipt-1];
      if(tpm0.Chg <= 0) return;
      auto& tpm1 = tj.Pts[ipt-2];
      if(tpm1.Chg <= 0) return;
      float chgp = tpp0.Chg + tpp1.Chg;
      float chgm = tpm0.Chg + tpm1.Chg;
      float asym = (chgp - chgm) / (chgp + chgm);
      if(asym > maxAsym) {
        maxAsym = asym;
        atPt = ipt;
      }
      if(prt) mf::LogVerbatim("TC")<<"THiQEP ipt " << ipt << " " << PrintPos(slc, tpp0)
              << " Chg asym " << asym;
      ++cnt;
      if(cnt > 10) break;
    } // ii
    if(atPt < 0) return;
    if(prt) mf::LogVerbatim("TC")<<" THiQEP trim points after "<<atPt;
    for(short ipt = atPt; ipt <= tj.EndPt[1]; ++ipt) UnsetUsedHits(slc, tj.Pts[ipt]);
    SetEndPoints(tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kTHiQEP] = true;
  } // TrimHiChgEndPts

  //////////////////////////////////////////
  void
  TrimHiMultEndPts(TCSlice& slc, Trajectory& tj)
  {
    // Trim high multiplicity points at the end of the trajectory
    if(!tcc.useAlg[kLEPhys]) return;
    if(!tcc.useAlg[kTHMEP]) return;
    if(tj.Strategy[kSlowing]) return;
    if(tj.EndFlag[1][kEndBragg]) return;
    if(tj.EndPt[1] - tj.EndPt[0] < 10) return;

    bool prt = tcc.dbgStp;

    unsigned short cnt = 0;
    float maxAsym = 0.4;
    short atPt = -1;
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      short ipt = tj.EndPt[1] - ii - 1;
      if(ipt < tj.EndPt[0] + 3) break;
      auto& tpp0 = tj.Pts[ipt];
      if(tpp0.Chg <= 0) return;
      auto& tpp1 = tj.Pts[ipt+1];
      if(tpp1.Chg <= 0) return;
      auto& tpm0 = tj.Pts[ipt-1];
      if(tpm0.Chg <= 0) return;
      auto& tpm1 = tj.Pts[ipt-2];
      if(tpm1.Chg <= 0) return;
      float nhtp = tpp0.Hits.size() + tpp1.Hits.size();
      float nhtm = tpm0.Hits.size() + tpm1.Hits.size();
      float asym = (nhtp - nhtm) / (nhtp + nhtm);
      if(asym > maxAsym) {
        maxAsym = asym;
        atPt = ipt;
      }
      ++cnt;
      if(cnt > 10) break;
    } // ii
    if(atPt < 0) return;
    auto& atTP = tj.Pts[atPt];
    if(atTP.Hits.size() == 1) ++atPt;
    if(prt) mf::LogVerbatim("TC") << "THMEP: Hit size asymmetry of " << maxAsym 
              << " exceeds 0.4 at " << PrintPos(slc, atTP) << ". Trim to the end";
    for(short ipt = atPt; ipt <= tj.EndPt[1]; ++ipt) UnsetUsedHits(slc, tj.Pts[ipt]);
    SetEndPoints(tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kTHMEP] = true;
  } // TrimHiMultEndPts

  ////////////////////////////////////////////////
  void ChkHiMultEndHits(TCSlice& slc, Trajectory& tj)
  {
    // mask off high multiplicity TPs at the end
    if(tcc.useAlg[kLEPhys]) return;
//    if(!tcc.useAlg[kCHMEH]) return;
    if(tj.EndFlag[1][kEndBragg]) return;
    if(tj.Pts.size() < 10) return;
//    if(tj.Pts[tj.EndPt[1]].AngleCode == 0) return;
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
        UnsetUsedHits(slc, tj.Pts[ipt]);
        ++cnt;
        continue;
      }
      break;
    } // ipt
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"CHMEH multiplicity cut "<<aveMult<<" number of TPs masked off "<<cnt;
    if(cnt > 0) {
//      tj.AlgMod[kCHMEH] = true;
      SetEndPoints(tj);
    }
  } // ChkHiMultEndHits

  ////////////////////////////////////////////////
  void LastEndMerge(TCSlice& slc, CTP_t inCTP)
  {
    // last ditch attempt to merge long straight broken trajectories by averaging
    // all points in the trajectory and applying tight angle and separation cuts.
    if(slc.tjs.size() < 2) return;
    if(!tcc.useAlg[kLastEndMerge]) return;

    bool prt = tcc.dbgAlg[kLastEndMerge];

    // create an averaged TP for each long Trajectory
    std::vector<TrajPoint> tjTP;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != inCTP) continue;
      if(tj.Pts.size() < 10) continue;
      if(tj.MCSMom < 100) continue;
      auto tjtp = CreateTPFromTj(slc, tj);
      if(tjtp.Chg < 0) continue;
      tjTP.push_back(tjtp);
    } // tj
    if(tjTP.size() < 2) return;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"inside LastEndMerge slice "<<slices.size()-1<<" inCTP "<<inCTP<<" tjTPs";
      for(auto& tjtp : tjTP) myprt<<" T"<<tjtp.Step;
    }

    for(unsigned short pt1 = 0; pt1 < tjTP.size() - 1; ++pt1) {
      auto& tp1 = tjTP[pt1];
      auto& tj1 = slc.tjs[tp1.Step - 1];
      if(tj1.AlgMod[kKilled]) continue;
      for(unsigned short pt2 = pt1 + 1; pt2 < tjTP.size(); ++pt2) {
        auto& tp2 = tjTP[pt2];
        auto& tj2 = slc.tjs[tp2.Step - 1];
        if(tj2.AlgMod[kKilled]) continue;
        float dang = DeltaAngle(tp1.Ang, tp2.Ang);
        // make an angle cut
        if(prt && dang < 0.5) mf::LogVerbatim("TC")<<" T"<<tj1.ID<<" T"<<tj2.ID<<" dang "<<dang;
        if(dang > 0.2) continue;
        // and an impact parameter cut
        unsigned short ipt1, ipt2;
        float ip12 = PointTrajDOCA(slc, tp1.Pos[0], tp1.Pos[1], tp2);
        float ip21 = PointTrajDOCA(slc, tp2.Pos[0], tp2.Pos[1], tp1);
        if(prt) mf::LogVerbatim("TC")<<" ip12 "<<ip12<<" ip21 "<<ip21;
        if(ip12 > 15 && ip21 > 15) continue;
        float minSep = 100;
        // find the separation considering dead wires
        TrajTrajDOCA(slc, tj1, tj2, ipt1, ipt2, minSep, false);
        if(minSep == 100) continue;
        if(ipt1 >= tj1.Pts.size() || ipt2 >= tj2.Pts.size()) continue;
        float dwc = DeadWireCount(slc, tj1.Pts[ipt1], tj2.Pts[ipt2]);
        if(prt) mf::LogVerbatim("TC")<<" minSep "<<minSep<<" dwc "<<dwc;
        minSep -= dwc;
        if(minSep > 5) continue;
        // finally require that the proximate points are close to the ends
        float sep10 = PosSep(tj1.Pts[ipt1].Pos, tj1.Pts[tj1.EndPt[0]].Pos);
        float sep11 = PosSep(tj1.Pts[ipt1].Pos, tj1.Pts[tj1.EndPt[1]].Pos);
        if(sep10 > 5 && sep11 > 5) continue;
        unsigned short end1 = 0;
        if(sep11 < sep10) end1 = 1;
        float sep20 = PosSep(tj2.Pts[ipt2].Pos, tj2.Pts[tj2.EndPt[0]].Pos);
        float sep21 = PosSep(tj2.Pts[ipt2].Pos, tj2.Pts[tj2.EndPt[1]].Pos);
        if(sep20 > 5 && sep21 > 5) continue;
        unsigned short end2 = 0;
        if(sep21 < sep20) end2 = 1;
        // don't merge if there is a kink
        if(tj1.EndFlag[end1][kEndKink] || tj2.EndFlag[end2][kEndKink]) continue;
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<"LEM: T"<<tj1.ID<<"_"<<PrintPos(slc, tp1);
          if(tj1.VtxID[end1] > 0) myprt<<"->2V"<<tj1.VtxID[end1];
          myprt<<" T"<<tj2.ID<<"_"<<PrintPos(slc, tp2);
          if(tj2.VtxID[end2] > 0) myprt<<"->2V"<<tj2.VtxID[end2];
          myprt<<" dang "<<std::setprecision(2)<<dang<<" ip12 "<<ip12;
          myprt<<" ip21 "<<ip21;
          myprt<<" minSep "<<minSep;
          myprt<<" end sep1 "<<sep10<<" "<<sep11;
          myprt<<" end sep2 "<<sep20<<" "<<sep21;
        } // prt
        if(tj1.VtxID[end1] > 0) {
          auto& vx2 = slc.vtxs[tj1.VtxID[end1] - 1];
          MakeVertexObsolete("LEM", slc, vx2, true);
        }
        if(tj2.VtxID[end2] > 0 && tj2.VtxID[end2] != tj1.VtxID[end1]) {
          auto& vx2 = slc.vtxs[tj2.VtxID[end2] - 1];
          MakeVertexObsolete("LEM", slc, vx2, true);
        }
        // remove Bragg flags
        tj1.EndFlag[end1][kEndBragg] = false;
        tj2.EndFlag[end2][kEndBragg] = false;
        unsigned int it1 = tj1.ID - 1;
        unsigned int it2 = tj2.ID - 1;
        if(!MergeAndStore(slc, it1, it2, tcc.dbgMrg)) continue;
        // set the AlgMod bit
        auto& ntj = slc.tjs[slc.tjs.size() - 1];
        ntj.AlgMod[kLastEndMerge] = true;
        // create a tp for this tj and add it to the list
        auto tjtp = CreateTPFromTj(slc, ntj);
        if(tjtp.Chg < 0) continue;
        if(prt) mf::LogVerbatim("TC")<<" added T"<<ntj.ID<<" to the merge list";
        tjTP.push_back(tjtp);
        break;
      } // pt1
    } // pt1
  } // LastEndMerge

  ////////////////////////////////////////////////
  void EndMerge(TCSlice& slc, CTP_t inCTP, bool lastPass)
  {
    // Merges trajectories end-to-end or makes vertices. Does a more careful check on the last pass

    if(slc.tjs.size() < 2) return;
    if(!tcc.useAlg[kMerge]) return;

    bool prt = (tcc.dbgMrg && tcc.dbgSlc && inCTP == debug.CTP);
    if(prt) mf::LogVerbatim("TC")<<"inside EndMerge slice "<<slices.size()-1<<" inCTP "<<inCTP<<" nTjs "<<slc.tjs.size()<<" lastPass? "<<lastPass;

    // Ensure that all tjs are in the same order
    short tccStepDir = 1;
    if(!tcc.modes[kStepPos]) tccStepDir = -1;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != inCTP) continue;
      if(tj.StepDir != tccStepDir) ReverseTraj(slc, tj);
    } // tj

    unsigned short maxShortTjLen = tcc.vtx2DCuts[0];

    // temp vector for checking the fraction of hits near a merge point
    std::vector<int> tjlist(2);

    float minChgRMS = 0.5 * (tcc.chargeCuts[1] + tcc.chargeCuts[2]);

    // iterate whenever a merge occurs since tjs will change. This is not necessary
    // when a vertex is created however.
    bool iterate = true;
    while(iterate) {
      iterate = false;
      for(unsigned int it1 = 0; it1 < slc.tjs.size(); ++it1) {
        auto& tj1 = slc.tjs[it1];
        if(tj1.AlgMod[kKilled]) continue;
        if(tj1.CTP != inCTP) continue;
        // don't try to merge high energy electrons
        if(tj1.PDGCode == 111) continue;
        for(unsigned short end1 = 0; end1 < 2; ++end1) {
          // no merge if there is a vertex at the end
          if(tj1.VtxID[end1] > 0) continue;
          // make a copy of tp1 so we can mess with it
          TrajPoint tp1 = tj1.Pts[tj1.EndPt[end1]];
          // do a local fit on the lastpass only using the last 3 points
          if(lastPass && tp1.NTPsFit > 3) {
            // make a local copy of the tj
            auto ttj = slc.tjs[it1];
            auto& lastTP = ttj.Pts[ttj.EndPt[end1]];
            // fit the last 3 points
            lastTP.NTPsFit = 3;
            FitTraj(slc, ttj);
            tp1 = ttj.Pts[ttj.EndPt[end1]];
          } // last pass
          bool isVLA = (tp1.AngleCode == 2);
          float bestFOM = 5;
          if(isVLA) bestFOM = 20;
          float bestDOCA;
          unsigned int imbest = UINT_MAX;
          for(unsigned int it2 = 0; it2 < slc.tjs.size(); ++it2) {
            if(it1 == it2) continue;
            auto& tj2 = slc.tjs[it2];
            // check for consistent direction
            if(tj1.StepDir != tj2.StepDir) continue;
            if(tj2.AlgMod[kKilled]) continue;
            if(tj2.CTP != inCTP) continue;
            // don't try to merge high energy electrons
            if(tj2.PDGCode == 111) continue;
            float olf = OverlapFraction(slc, tj1, tj2);
            if(olf > 0.25) continue;
            unsigned short end2 = 1 - end1;
            // check for a vertex at this end
            if(tj2.VtxID[end2] > 0) continue;
            TrajPoint& tp2 = tj2.Pts[tj2.EndPt[end2]];
            TrajPoint& tp2OtherEnd = tj2.Pts[tj2.EndPt[end1]];
            // ensure that the other end isn't closer
            if(std::abs(tp2OtherEnd.Pos[0] - tp1.Pos[0]) < std::abs(tp2.Pos[0] - tp1.Pos[0])) continue;
            // ensure that the order is correct
            if(tj1.StepDir > 0) {
              if(tp2.Pos[0] < tp1.Pos[0] - 2) continue;
            } else {
              if(tp2.Pos[0] > tp1.Pos[0] + 2) continue;
            }
            // ensure that there is a signal on most of the wires between these points
            if(!SignalBetween(slc, tp1, tp2, 0.8)) {
              continue;
            }
            // Find the distance of closest approach for small angle merging
            // Inflate the doca cut if we are bridging a block of dead wires
            float dang = DeltaAngle(tp1.Ang, tp2.Ang);
            float doca = 15;
            if(isVLA) {
              // compare the minimum separation between Large Angle trajectories using a generous cut
              unsigned short ipt1, ipt2;
              TrajTrajDOCA(slc, tj1, tj2, ipt1, ipt2, doca);
              //              if(prt) mf::LogVerbatim("TC")<<" isVLA check ipt1 "<<ipt1<<" ipt2 "<<ipt2<<" doca "<<doca;
            } else {
              // small angle
              doca = PointTrajDOCA(slc, tp1.Pos[0], tp1.Pos[1], tp2);
            }
            float fom = dang * doca;
            if(fom < bestFOM) {
              bestFOM = fom;
              bestDOCA = doca;
              imbest = it2;
            }
          } // it2
          // No merge/vertex candidates
          if(imbest == UINT_MAX) continue;

          // Make angle adjustments to tp1.
          unsigned int it2 = imbest;
          auto& tj2 = slc.tjs[imbest];
          unsigned short end2 = 1 - end1;
          bool loMCSMom = (tj1.MCSMom + tj2.MCSMom) < 150;
          // Don't use the angle at the end Pt for high momentum long trajectories in case there is a little kink at the end
          if(tj1.Pts.size() > 50 && tj1.MCSMom > 100) {
            if(end1 == 0) {
              tp1.Ang = tj1.Pts[tj1.EndPt[0] + 2].Ang;
            } else {
              tp1.Ang = tj1.Pts[tj1.EndPt[1] - 2].Ang;
            }
          } else if(loMCSMom) {
            // Low momentum - calculate the angle using the two Pts at the end
            unsigned short pt1, pt2;
            if(end1 == 0) {
              pt1 = tj1.EndPt[0];
              pt2 = pt1 + 1;
            } else {
              pt2 = tj1.EndPt[1];
              pt1 = pt2 - 1;
            }
            TrajPoint tpdir;
            if(MakeBareTrajPoint(slc, tj1.Pts[pt1], tj1.Pts[pt2], tpdir)) tp1.Ang = tpdir.Ang;
          } // low MCSMom
          // Now do the same for tj2
          TrajPoint tp2 = tj2.Pts[tj2.EndPt[end2]];
          if(tj2.Pts.size() > 50 && tj2.MCSMom > 100) {
            if(end1 == 0) {
              tp2.Ang = tj2.Pts[tj2.EndPt[0] + 2].Ang;
            } else {
              tp2.Ang = tj2.Pts[tj2.EndPt[1] - 2].Ang;
            }
          } else if(loMCSMom) {
            // Low momentum - calculate the angle using the two Pts at the end
            unsigned short pt1, pt2;
            if(end2 == 0) {
              pt1 = tj2.EndPt[0];
              pt2 = pt1 + 1;
            } else {
              pt2 = tj2.EndPt[1];
              pt1 = pt2 - 1;
            }
            TrajPoint tpdir;
            if(MakeBareTrajPoint(slc, tj2.Pts[pt1], tj2.Pts[pt2], tpdir)) tp2.Ang = tpdir.Ang;
          } // low MCSMom

          if(!isVLA && !SignalBetween(slc, tp1, tp2, 0.99)) continue;

          // decide whether to merge or make a vertex
          // protect against angles > pi/2
          float dang = acos(DotProd(tp1.Dir, tp2.Dir));
          float sep = PosSep(tp1.Pos, tp2.Pos);
          // ignore this pair if the gap between them is much longer than the length of the shortest Tj
          float len1 = TrajLength(slc.tjs[it1]);
          float len2 = TrajLength(slc.tjs[it2]);
          if(len1 < len2 && sep > 3 * len1) continue;
          if(len2 < len1 && sep > 3 * len2) continue;

          // default cuts for locMCSMom condition
          float dangCut = 1;
          float docaCut = 2;
          float kinkSig = -1;
          if(!loMCSMom) {
            unsigned short nPtsFit = tcc.kinkCuts[0];
            bool useChg = (tcc.kinkCuts[2] > 0);
            kinkSig = KinkSignificance(slc, tj1, end1, tj2, end2, nPtsFit, useChg, prt);
          }
          docaCut = 1.5;
          if(isVLA) docaCut = 15;
          float chgPull = 0;
          if(tp1.AveChg > tp2.AveChg) {
            chgPull = (tp1.AveChg / tp2.AveChg - 1) / minChgRMS;
          } else {
            chgPull = (tp2.AveChg / tp1.AveChg - 1) / minChgRMS;
          }
          // open up the cuts on the last pass
          float chgFracCut = tcc.vtx2DCuts[8];
          float chgPullCut = tcc.chargeCuts[0];
          if(lastPass) {
            docaCut *= 2;
            chgFracCut *= 0.5;
            chgPullCut *= 1.5;
          }

          // check the merge cuts. Start with doca and dang requirements
          bool doMerge = bestDOCA < docaCut && dang < dangCut;
          bool showerTjs = tj1.PDGCode == 11 || tj2.PDGCode == 11;
          bool hiMCSMom = tj1.MCSMom > 200 || tj2.MCSMom > 200;
          // add a charge similarity requirement if not shower-like or low momentum or not LA
          if(doMerge && !showerTjs && hiMCSMom && chgPull > tcc.chargeCuts[0] && !isVLA) doMerge = false;
          // ignore the charge pull cut if both are high momentum and dang is really small
          if(!doMerge && tj1.MCSMom > 900 && tj2.MCSMom > 900 && dang < 0.1 && bestDOCA < docaCut) doMerge = true;

          // do not merge if chgPull is really high
          if(doMerge && chgPull > 2 * chgPullCut) doMerge = false;
          float dwc = DeadWireCount(slc, tp1, tp2);

          if(doMerge) {
            if(lastPass) {
              // last pass cuts are looser but ensure that the tj after merging meets the quality cut
              float npwc = NumPtsWithCharge(slc, tj1, true) + NumPtsWithCharge(slc, tj2, true);
              auto& tp1OtherEnd = tj1.Pts[tj1.EndPt[1 - end1]];
              auto& tp2OtherEnd = tj2.Pts[tj2.EndPt[1 - end2]];
              float nwires = std::abs(tp1OtherEnd.Pos[0] - tp2OtherEnd.Pos[0]);
              if(nwires == 0) nwires = 1;
              float hitFrac = npwc / nwires;
              doMerge = (hitFrac > tcc.qualityCuts[0]);
            } else {
              // don't merge if the gap between them is longer than the length of the shortest Tj
              if(len1 < len2) {
                if(sep > len1) doMerge = false;
              } else {
                if(sep > len2) doMerge = false;
              }
              if(prt) mf::LogVerbatim("TC")<<" merge check sep "<<sep<<" len1 "<<len1<<" len2 "<<len2<<" dead wire count "<<dwc<<" Merge? "<<doMerge;
            } // not lastPass
          } // doMerge

          // Require a large charge fraction near a merge point
          tjlist[0] = slc.tjs[it1].ID;
          tjlist[1] = slc.tjs[it2].ID;
          float chgFrac = ChgFracNearPos(slc, tp1.Pos, tjlist);
          if(doMerge && bestDOCA > 1 && chgFrac < chgFracCut) doMerge = false;

          // Check the MCSMom asymmetry and don't merge if it is higher than the user-specified cut
          float momAsym = std::abs(tj1.MCSMom - tj2.MCSMom) / (float)(tj1.MCSMom + tj2.MCSMom);
          if(doMerge && momAsym > tcc.vtx2DCuts[9]) doMerge = false;

          // be more lenient with short, slowing trajectories
          if(!doMerge ) {
            bool isSlowing = (tj1.Strategy[kSlowing] || tj2.Strategy[kSlowing]);
            bool isShort = (len1 < 40 && len2 < 40);
            if(isSlowing && isShort && bestDOCA < 0.1 && dang < 0.3) doMerge = true;
            // Remove the Bragg Peak Flag?
            if(tj1.EndFlag[end1][kEndBragg]) tj1.EndFlag[end1][kEndBragg] = false;
            if(tj2.EndFlag[end2][kEndBragg]) tj2.EndFlag[end2][kEndBragg] = false;
          } // !doMerge

          if(doMerge && (tj1.EndFlag[end1][kEndKink] || tj2.EndFlag[end2][kEndKink])) {
            // don't merge if a kink exists and the tjs are not too long
            if(len1 < 40 && len2 < 40) doMerge = false;
            // Kink on one + Bragg at other end of the other
            if(tj1.EndFlag[end1][kEndKink] && tj2.EndFlag[1-end2][kEndBragg]) doMerge = false;
            if(tj1.EndFlag[1-end1][kEndBragg] && tj2.EndFlag[end2][kEndKink]) doMerge = false;
          }

          // decide if we should make a vertex instead
          bool doVtx = false;
          if(!doMerge) {
            // check for a significant kink
            doVtx = (kinkSig > tcc.kinkCuts[1]);
            // and a less significant kink but very close separation
            doVtx = (kinkSig > 0.5 * tcc.kinkCuts[1] && sep < 2);
          } // !doMerge

          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<"  EM: T"<<slc.tjs[it1].ID<<"_"<<end1<<" - T"<<slc.tjs[it2].ID<<"_"<<end2<<" tp1-tp2 "<<PrintPos(slc, tp1)<<"-"<<PrintPos(slc, tp2);
            myprt<<" FOM "<<std::fixed<<std::setprecision(2)<<bestFOM;
            myprt<<" DOCA "<<std::setprecision(1)<<bestDOCA;
            myprt<<" cut "<<docaCut<<" isVLA? "<<isVLA;
            myprt<<" dang "<<std::setprecision(2)<<dang<<" dangCut "<<dangCut;
            myprt<<" chgPull "<<std::setprecision(1)<<chgPull<<" Cut "<<chgPullCut;
            myprt<<" chgFrac "<<std::setprecision(2)<<chgFrac;
            myprt<<" momAsym "<<momAsym;
            myprt<<" kinkSig "<<std::setprecision(1)<<kinkSig;
            myprt<<" doMerge? "<<doMerge;
            myprt<<" doVtx? "<<doVtx;
          }

          if(bestDOCA > docaCut) continue;

          if(doMerge) {
            if(prt) mf::LogVerbatim("TC")<<"  Merge ";
            bool didMerge = false;
            if(end1 == 1) {
              didMerge = MergeAndStore(slc, it1, it2, tcc.dbgMrg);
            } else {
              didMerge = MergeAndStore(slc, it2, it1, tcc.dbgMrg);
            }
            if(didMerge) {
              // Set the end merge flag for the killed trajectories to aid tracing merges
              tj1.AlgMod[kMerge] = true;
              tj2.AlgMod[kMerge] = true;
              iterate = true;
            } // Merge and store successfull
            else {
              if(prt) mf::LogVerbatim("TC")<<"  MergeAndStore failed ";
            }
          } else if(doVtx) {
            // create a vertex instead if it passes the vertex cuts
            VtxStore aVtx;
            aVtx.CTP = slc.tjs[it1].CTP;
            aVtx.ID = slc.vtxs.size() + 1;
            // keep it simple if tp1 and tp2 are very close or if the angle between them
            // is small
            if(prt) {
              mf::LogVerbatim("TC")<<"  candidate 2V"<<aVtx.ID<<" dang "<<dang<<" sep "<<PosSep(tp1.Pos, tp2.Pos);
            }
            bool fix2V = (PosSep(tp1.Pos, tp2.Pos) < 3 || dang < 0.1);
            if(fix2V) {
              aVtx.Pos[0] = 0.5 * (tp1.Pos[0] + tp2.Pos[0]);
              aVtx.Pos[1] = 0.5 * (tp1.Pos[1] + tp2.Pos[1]);
              aVtx.Stat[kFixed] = true;
              aVtx.PosErr[0] = std::abs(tp1.Pos[0] - tp2.Pos[0]);
              aVtx.PosErr[1] = std::abs(tp1.Pos[1] - tp2.Pos[1]);
            } else {
              float sepCut = tcc.vtx2DCuts[1];
              bool tj1Short = (slc.tjs[it1].EndPt[1] - slc.tjs[it1].EndPt[0] < maxShortTjLen);
              bool tj2Short = (slc.tjs[it2].EndPt[1] - slc.tjs[it2].EndPt[0] < maxShortTjLen);
              if(tj1Short || tj2Short) sepCut = tcc.vtx2DCuts[1];
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
              // but not too short
              if(tj1Short && len1 > 4) {
                TrajPoint otp1 = slc.tjs[it1].Pts[slc.tjs[it1].EndPt[1-end1]];
                if(PosSep2(otp1.Pos, aVtx.Pos) < PosSep2(tp1.Pos, aVtx.Pos)) continue;
              }
              if(tj2Short && len2 > 4) {
                TrajPoint otp2 = slc.tjs[it2].Pts[slc.tjs[it2].EndPt[1-end2]];
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
            slc.tjs[it1].VtxID[end1] = aVtx.ID;
            slc.tjs[it2].VtxID[end2] = aVtx.ID;
            if(!aVtx.Stat[kFixed] && !FitVertex(slc, aVtx, prt)) {
              // back out
              slc.tjs[it1].VtxID[end1] = 0;
              slc.tjs[it2].VtxID[end2] = 0;
              if(prt) mf::LogVerbatim("TC")<<" Vertex fit failed ";
              continue;
            }
            aVtx.NTraj = 2;
            aVtx.Pass = slc.tjs[it1].Pass;
            aVtx.Topo = end1 + end2;
            tj1.AlgMod[kMerge] = true;
            tj2.AlgMod[kMerge] = true;
            if(!StoreVertex(slc, aVtx)) continue;
            SetVx2Score(slc);
            if(prt) {
              auto& newVx = slc.vtxs[slc.vtxs.size() - 1];
              mf::LogVerbatim("TC")<<"  New 2V"<<newVx.ID<<" at "<<(int)newVx.Pos[0]<<":"<<(int)(newVx.Pos[1]/tcc.unitsPerTick)<<" Score "<<newVx.Score;
            }
            // check the score and kill it if it is below the cut
            // Don't kill the vertex in this function since it is
            // called before short trajectories are reconstructed
            auto& newVx2 = slc.vtxs[slc.vtxs.size() - 1];
            if(newVx2.Score < tcc.vtx2DCuts[7] && CompatibleMerge(slc, tj1, tj2, prt)) {
              if(prt) {
                mf::LogVerbatim myprt("TC");
                myprt<<"  Bad vertex: Bad score? "<<(newVx2.Score < tcc.vtx2DCuts[7]);
                myprt<<" cut "<<tcc.vtx2DCuts[7];
                myprt<<" CompatibleMerge? "<<CompatibleMerge(slc, tj1, tj2, prt);
              }
              slc.tjs[it1].VtxID[end1] = 0;
              slc.tjs[it2].VtxID[end2] = 0;
              MakeVertexObsolete("EM", slc, newVx2, true);
              bool didMerge = false;
              if(end1 == 1) {
                didMerge = MergeAndStore(slc, it1, it2, tcc.dbgMrg);
              } else {
                didMerge = MergeAndStore(slc, it2, it1, tcc.dbgMrg);
              }
              if(didMerge) {
                // Set the end merge flag for the killed trajectories to aid tracing merges
                tj1.AlgMod[kMerge] = true;
                tj1.AlgMod[kMerge] = true;
                iterate = true;
              } // Merge and store successfull
              else {
                if(prt) mf::LogVerbatim("TC")<<"  MergeAndStore failed ";
              }
            } // OK score
          } // create a vertex
          if(tj1.AlgMod[kKilled]) break;
        } // end1
      } // it1
    } // iterate
  } // EndMerge

  ////////////////////////////////////////////////
  void ChkStop(TCSlice& slc, Trajectory& tj)
  {
    // Sets the EndFlag[kEndBragg] bits on the trajectory by identifying the Bragg peak
    // at each end. This function checks both ends, finding the point with the highest charge nearest the
    // end and considering the first (when end = 0) 4 points or last 4 points (when end = 1). The next
    // 5 - 10 points (fChkStop[0]) are fitted to a line, Q(x - x0) = Qo + (x - x0) * slope where x0 is the
    // wire position of the highest charge point. A large negative slope indicates that there is a Bragg
    // peak at the end.

    tj.EndFlag[0][kEndBragg] = false;
    tj.EndFlag[1][kEndBragg] = false;
    if(!tcc.useAlg[kChkStop]) return;
    if(tcc.chkStopCuts[0] < 0) return;

    if(tj.Strategy[kStiffEl]) return;

    // ignore trajectories that are very large angle at both ends
    if(tj.Pts[tj.EndPt[0]].AngleCode == 2 || tj.Pts[tj.EndPt[1]].AngleCode == 2) return;

    unsigned short nPtsToCheck = tcc.chkStopCuts[1];
    if(tj.Pts.size() < 6) return;

    bool prt = (tcc.dbgStp || tcc.dbgAlg[kChkStop]);

    // find the highest charge hit in the first 3 points at each end
    for(unsigned short end = 0; end < 2; ++end) {
      tj.EndFlag[end][kEndBragg] = false;
      short dir = 1 - 2 * end;
      // find the point with the highest charge considering the first 3 points
      float big = 0;
      unsigned short hiPt = 0;
      float wire0 = 0;
      for(unsigned short ii = 0; ii < 5; ++ii) {
        short ipt = tj.EndPt[end] + ii * dir;
        if(ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) break;
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg > big) {
          big = tp.Chg;
          wire0 = tp.Pos[0];
          hiPt = ipt;
        }
      } // ii
      // ensure that the high point is closer to the end that is being
      // considered than to the other end
      short dpt0 = hiPt - tj.EndPt[0];
      short dpt1 = tj.EndPt[1] - hiPt;
      if(end == 0 && dpt1 <= dpt0) continue;
      if(end == 1 && dpt0 <= dpt1) continue;
      float prevChg = big;
      // prepare to do the fit
      Point2_t inPt;
      Vector2_t outVec, outVecErr;
      float chgErr, chiDOF;
      // Initialize
      Fit2D(0, inPt, chgErr, outVec, outVecErr, chiDOF);
      unsigned short cnt = 0;
      for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        short ipt = hiPt + ii * dir;
        if(ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) break;
        TrajPoint& tp = tj.Pts[ipt];
        float tpChg = TpSumHitChg(slc, tp);
        if(tpChg == 0) continue;
        // quit if the charge is much larger than the previous charge
        if(tpChg > 1.5 * prevChg) continue;
        prevChg = tpChg;
        // Accumulate and save points
        inPt[0] = std::abs(tp.Pos[0] - wire0);
        inPt[1] = tpChg;
        // Assume 20% point-to-point charge fluctuations
        chgErr = 0.2 * tpChg;
        if(!Fit2D(2, inPt, chgErr, outVec, outVecErr, chiDOF)) break;
        ++cnt;
        if(cnt == nPtsToCheck) break;
      } // ii
      if(cnt < nPtsToCheck) continue;
      // do the fit and get the results
      if(!Fit2D(-1, inPt, chgErr, outVec, outVecErr, chiDOF)) continue;
      // check for really bad chidof indicating a major failure
      if(chiDOF > 500) continue;
      // The charge slope is negative for a stopping track in the way that the fit was constructed.
      // Flip the sign so we can make a cut against tcc.chkStopCuts[0] which is positive.
      outVec[1] = -outVec[1];
      // BB 9/16/2020 Lower confidence level from 99% to ~97%
      bool itStops = (outVec[1] > tcc.chkStopCuts[0] && chiDOF < tcc.chkStopCuts[2] && outVec[1] > 1.8 * outVecErr[1]);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"ChkStop: end "<<end;
        myprt<<" at "<<PrintPos(slc, tj.Pts[hiPt]);
        myprt<<" slope "<<outVec[1];
        if(outVec[1] > tcc.chkStopCuts[0]) myprt<<" OK";
        myprt<<" +/- "<<outVecErr[1];
        if(outVec[1] > 2.5 * outVecErr[1]) myprt<<" OK";
        myprt<<" fit chidof "<<chiDOF;
        if(chiDOF < tcc.chkStopCuts[2]) myprt<<" OK";
        myprt<<" AveChg "<<std::nearbyint(outVec[0]);
        if(itStops) myprt<<" -> Stops";
      } // prt
      if(itStops) {
        tj.EndFlag[end][kEndBragg] = true;
        tj.AlgMod[kChkStop] = true;
        if(tj.PDGCode == 11) tj.PDGCode = 0;
        // Put the charge at the end into tp.AveChg
        tj.Pts[tj.EndPt[end]].AveChg = outVec[0];
      } // itStops
    } // end
  } // ChkStop

  //////////////////////TY://////////////////////////
  bool ChkMichel(TCSlice& slc, Trajectory& tj, unsigned short& lastGoodPt){

    if(!tcc.useAlg[kMichel]) return false;
    if(tj.PDGCode == 11 || tj.PDGCode == 111) return false;

    bool prt = (tcc.dbgStp || tcc.dbgAlg[kMichel]);

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
      tj.AlgMod[kMichel] = true;
      return true;
    }
    else{
      return false;
    }
  } // ChkMichel

  //////////////////////////////////////////
  bool MakeJunkTraj(TCSlice& slc, std::vector<unsigned int> tHits)
  {
    if(!tcc.useAlg[kJunkTj]) return false;
    // Make a crummy trajectory using the provided hits

    if(tHits.size() < 2) return false;

    bool prt = false;
    if(tcc.dbgAlg[kJunkTj]) {
      for(unsigned short ii = 0; ii < tHits.size(); ++ii) {
        if(slc.slHits[tHits[ii]].allHitsIndex == debug.Hit) {
          prt = true;
          break;
        }
      } // ii
      if(prt) std::cout<<"MakeJunkTraj found debug hit\n";
    } // tcc.dbgAlg[kJunkTj]

    // Start the trajectory using the first and last hits to
    // define a starting direction. Use the last pass settings
    Trajectory work;
    unsigned short pass = tcc.minPts.size() - 1;
    if(!StartTraj(slc, work, tHits[0], tHits[tHits.size()-1], pass)) return false;
    // make a TP for every hit
    work.Pts.resize(tHits.size());
    // and put a hit into each one
    for(unsigned short ii = 0; ii < tHits.size(); ++ii) {
      auto& tp = work.Pts[ii];
      unsigned int iht = tHits[ii];
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      tp.CTP = EncodeCTP(hit.WireID());
      if(tp.CTP != work.CTP) return false;
      tp.Hits.push_back(iht);
      tp.UseHit[0] = true;
      // don't use DefineHitPos here because the angle isn't really known yet. Just
      // define enough information to do a fit
      tp.HitPos[0] = hit.WireID().Wire;
      tp.HitPos[1] = hit.PeakTime() * tcc.unitsPerTick;
      tp.HitPosErr2 = 100;
      tp.Chg = hit.Integral();
      tp.Step = ii;
      tp.NTPsFit = tHits.size();
      // flag long-pulse hits
      if(LongPulseHit(hit)) tp.Environment[kEnvUnusedHits] = true;
    } // ii
    work.EndPt[0] = 0;
    work.EndPt[1] = tHits.size() - 1;
    // do an initial fit. The fit results are put in the last TP.
    FitTraj(slc, work);
    auto& lastTP = work.Pts.back();
    // Prepare to sort along the general direction. First find the
    // along and transverse position (Delta) of each TP and calculate DeltaRMS
    double sum = 0.;
    double sum2 = 0.;
    for(auto& tp : work.Pts) {
      Point2_t at;
      FindAlongTrans(lastTP.Pos, lastTP.Dir, tp.HitPos, at);
      sum += at[1];
      sum2 += at[1] * at[1];
      // store the along distance in AveChg for now
      tp.AveChg = at[0];
      tp.Delta = at[1];
      if(tp.Step != lastTP.Step) {
        tp.FitChi = lastTP.FitChi;
        tp.Dir = lastTP.Dir;
        tp.Ang = lastTP.Ang;
        tp.Pos[0] = lastTP.Pos[0] + at[0] * lastTP.Dir[0];
        tp.Pos[1] = lastTP.Pos[1] + at[0] * lastTP.Dir[1];
      }
    } // tp
    double npts = tHits.size();
    sum /= npts;
    double arg = sum2 - npts * sum * sum;
    if(arg <= 0) return false;
    float rms = sqrt(arg) / (npts - 1);
    // apply a loose Delta cut
    float transCut = sum + 3 * rms;
    std::vector<SortEntry> sortVec;
    SortEntry se;
    work.TotChg = 0;
    work.NeedsUpdate = false;
    for(auto& tp : work.Pts) {
      if(tp.Delta > transCut && !tp.Environment[kEnvUnusedHits]) {
        work.NeedsUpdate = true;
        continue;
      }
      se.index = tp.Step;
      se.val = tp.AveChg;
      sortVec.push_back(se);
      tp.DeltaRMS = rms;
      work.TotChg += tp.Chg;
    } // tp
    if(sortVec.size() < 3) return false;
    std::sort(sortVec.begin(), sortVec.end(), valsDecreasing);
    std::vector<TrajPoint> ntps(sortVec.size());
    for(unsigned short ipt = 0; ipt < sortVec.size(); ++ipt) ntps[ipt] = work.Pts[sortVec[ipt].index];
    work.Pts = ntps;
    sum = work.TotChg / (double)ntps.size();
    if(work.NeedsUpdate) {
      work.EndPt[1] = work.Pts.size() - 1;
      UpdateTraj(slc, work);
    } // needs update
    work.AlgMod[kJunkTj] = true;
    work.IsGood = true;
    if(prt) {
      PrintTrajectory("MJT", slc, work, USHRT_MAX);
    }
    // Store it
    return StoreTraj(slc, work);
  } // MakeJunkTraj

  ////////////////////////////////////////////////
  void
  MergeShortWithJunk(TCSlice& slc, CTP_t inCTP)
  {
    // Merge short Tjs that are close to a Junk Tj into the Junk Tj
    if(!tcc.useAlg[kMrgShortJunk]) return;

    bool prt = tcc.dbgAlg[kMrgShortJunk];

    std::vector<int> junkList;
    std::vector<int> shortList;
    for(auto& tj : slc.tjs) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kJunkTj]) {
        junkList.push_back(tj.ID);
      } else if(NumPtsWithCharge(slc, tj, false) < 4) {
        shortList.push_back(tj.ID);
      }
    } // tj
    if(junkList.empty() || shortList.empty()) return;
    if(prt) mf::LogVerbatim("TC")<<"MSWJ: found "<<junkList.size()<<" junk and "<<shortList.size()<<" short tjs"
            " in CTP " << inCTP << " Separation cut " << tcc.showerTag[1];
    for(auto jtid : junkList) {
      auto& jtj = slc.tjs[jtid - 1];
      auto& jtp0 = jtj.Pts[jtj.EndPt[0]];
      for(auto stid : shortList) {
        auto& stj = slc.tjs[stid - 1];
        auto& stp0 = stj.Pts[0];
        // make a rough separation cut
        if(PosSep(jtp0.HitPos, stp0.Pos) > 50) continue;
        float doca = tcc.showerTag[1];
        unsigned short jpt, spt;
        TrajTrajDOCA(slc, jtj, stj, jpt, spt, doca, false);
        if (doca == tcc.showerTag[1]) continue;
        bool sb = SignalBetween(slc, jtj.Pts[jpt], stj.Pts[spt], 0.5);
        if(prt) mf::LogVerbatim("TC") << " Junk T" << jtid << " short T" << stid << " doca " << doca
                << " SignalBetween? " << sb;
        if (!sb) continue;
        // merge them
        auto jhits = PutTrajHitsInVector(jtj, kUsedHits);
        auto shits = PutTrajHitsInVector(stj, kUsedHits);
        MakeTrajectoryObsolete(slc, (unsigned int)(jtid - 1));
        MakeTrajectoryObsolete(slc, (unsigned int)(stid - 1));
        jhits.insert(jhits.end(), shits.begin(), shits.end());
        // make local TPs that will span the range of wires for both Tjs
        auto loPos = jtj.Pts[0].Pos;
        auto hiPos = loPos;
        for(auto& stp : stj.Pts) {
          if(stp.Pos[0] < loPos[0]) loPos = stp.Pos;
          if(stp.Pos[0] > hiPos[0]) hiPos = stp.Pos;
        }
        for(auto& jtp : jtj.Pts) {
          if(jtp.Pos[0] < loPos[0]) loPos = jtp.Pos;
          if(jtp.Pos[0] > hiPos[0]) hiPos = jtp.Pos;
        }
        unsigned int loWire = std::nearbyint(loPos[0]);
        unsigned int hiWire = std::nearbyint(hiPos[0]);
        float delta = std::abs(hiPos[1] - loPos[1]);
        TrajPoint ltp;
        ltp.CTP = jtj.CTP;
        if(!MakeBareTrajPoint(loPos, hiPos, ltp)) continue;
        for (unsigned int wire = loWire; wire <= hiWire; ++wire) {
          MoveTPToWire(ltp, (float)wire);
          if(!FindCloseHits(slc, ltp, delta, kUnusedHits)) continue;
          for(auto iht : ltp.Hits) {
            if(slc.slHits[iht].InTraj != 0) continue;
            if (std::find(jhits.begin(), jhits.end(), iht) == jhits.end()) jhits.push_back(iht);
          } // iht
        } // wire
        if (!MakeJunkTraj(slc, jhits)) {
          if (prt) mf::LogVerbatim("TC") << " MakeJunkTraj failed";
          break;
        }
        if (prt) mf::LogVerbatim("TC") << "Created junk T" << slc.tjs.back().ID;
        slc.tjs.back().AlgMod[kMrgShortJunk] = true;
        break;
      } // stid
    } // jtid
  } // MergeShortWithJunk

  //////////////////////////////////////////
  bool
  BraggSplit(TCSlice& slc, unsigned short itj)
  {
    // Searches the stored trajectory for a Bragg Peak and kink and splits it
    if (!tcc.useAlg[kBraggSplit]) return false;
    if (itj > slc.tjs.size() - 1) return false;
    if (tcc.chkStopCuts.size() < 4) return false;
    if (tcc.chkStopCuts[3] <= 0) return false;
    unsigned short nPtsToCheck = tcc.chkStopCuts[1];
    auto& tj = slc.tjs[itj];
    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    if (npwc < 4) return false;
    if (npwc < nPtsToCheck) nPtsToCheck = npwc;
    // do a rough ChgPull check first
    float maxPull = 2;
    unsigned short maxPullPt = USHRT_MAX;
    for (unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.ChgPull < maxPull) continue;
      maxPull = tp.ChgPull;
      maxPullPt = ipt;
    } // ipt
    if (maxPullPt == USHRT_MAX) return false;
    short dpt;
    if (maxPullPt < 0.5 * (tj.EndPt[0] + tj.EndPt[1])) { dpt = maxPullPt - tj.EndPt[0]; }
    else {
      dpt = tj.EndPt[1] - maxPullPt;
    }
    if (dpt < 3) return false;
    bool prt = (tcc.dbgSlc && (tcc.dbgStp || tcc.dbgAlg[kBraggSplit]));
    if (prt)
      mf::LogVerbatim("TC") << "BS: T" << tj.ID << " maxPull " << maxPull << " at "
                            << PrintPos(slc, tj.Pts[maxPullPt]) << " dpt " << dpt;
    unsigned short breakPt = USHRT_MAX;
    float bestFOM = tcc.chkStopCuts[3];
    unsigned short bestBragg = 0;
    unsigned short nPtsFit = tcc.kinkCuts[0];
    TrajPoint tp1, tp2;
    ParFit chgFit1, chgFit2;
    for (unsigned short ipt = maxPullPt - 2; ipt <= maxPullPt + 2; ++ipt) {
      FitTraj(slc, tj, ipt - 1, nPtsFit, -1, tp1);
      if (tp1.FitChi > 10) continue;
      FitTraj(slc, tj, ipt + 1, nPtsFit, 1, tp2);
      if (tp2.FitChi > 10) continue;
      float dang = std::abs(tp1.Ang - tp2.Ang);
      FitPar(slc, tj, ipt - 1, nPtsToCheck, -1, chgFit1, 1);
      if (chgFit1.ChiDOF > 100) continue;
      chgFit1.ParSlp = -chgFit1.ParSlp;
      FitPar(slc, tj, ipt + 1, nPtsToCheck, 1, chgFit2, 1);
      if (chgFit2.ChiDOF > 100) continue;
      chgFit2.ParSlp = -chgFit2.ParSlp;
      // require a large positive slope on at least one side
      if (chgFit1.ParSlp < tcc.chkStopCuts[0] && chgFit2.ParSlp < tcc.chkStopCuts[0]) continue;
      // assume it is on side 1
      unsigned short bragg = 1;
      float bchi = chgFit1.ChiDOF;
      if (chgFit2.ParSlp > chgFit1.ParSlp) {
        bragg = 2;
        bchi = chgFit2.ChiDOF;
      }
      float chgAsym = std::abs(chgFit1.Par0 - chgFit2.Par0) / (chgFit1.Par0 + chgFit2.Par0);
      float slpAsym = std::abs(chgFit1.ParSlp - chgFit2.ParSlp) / (chgFit1.ParSlp + chgFit2.ParSlp);
      if (bchi < 1) bchi = 1;
      float fom = 10 * dang * chgAsym * slpAsym / bchi;
      if (prt) {
        mf::LogVerbatim myprt("TC");
        myprt << "pt " << PrintPos(slc, tj.Pts[ipt]) << " " << std::setprecision(2) << dang;
        myprt << " chg1 " << (int)chgFit1.Par0 << " slp " << chgFit1.ParSlp << " chi "
              << chgFit1.ChiDOF;
        myprt << " chg2 " << (int)chgFit2.Par0 << " slp " << chgFit2.ParSlp << " chi "
              << chgFit2.ChiDOF;
        myprt << " chgAsym " << chgAsym;
        myprt << " slpAsym " << slpAsym;
        myprt << " fom " << fom;
        myprt << " bragg " << bragg;
      }
      if (fom < bestFOM) continue;
      bestFOM = fom;
      breakPt = ipt;
      bestBragg = bragg;
    } // ipt
    if (breakPt == USHRT_MAX) return false;
    if (prt)
      mf::LogVerbatim("TC") << " breakPt " << PrintPos(slc, tj.Pts[breakPt]) << " bragg "
                            << bestBragg;
    // Create a vertex at the break point
    VtxStore aVtx;
    aVtx.Pos = tj.Pts[breakPt].Pos;
    aVtx.NTraj = 2;
    aVtx.Pass = tj.Pass;
    aVtx.Topo = 12;
    aVtx.ChiDOF = 0;
    aVtx.CTP = tj.CTP;
    aVtx.ID = slc.vtxs.size() + 1;
    aVtx.Stat[kFixed] = true;
    unsigned short ivx = slc.vtxs.size();
    if (!StoreVertex(slc, aVtx)) return false;
    if (!SplitTraj(slc, itj, breakPt, ivx, prt)) {
      if (prt) mf::LogVerbatim("TC") << "BS: Failed to split trajectory";
      MakeVertexObsolete("BS", slc, slc.vtxs[ivx], false);
      return false;
    }
    SetVx2Score(slc);
    slc.tjs[itj].AlgMod[kBraggSplit] = true;
    unsigned short otj = slc.tjs.size() - 1;
    if (bestBragg == 2) std::swap(itj, otj);
    slc.tjs[itj].PDGCode = 211;
    slc.tjs[itj].EndFlag[1][kEndBragg] = true;
    slc.tjs[otj].PDGCode = 13;
    return true;
  } // BraggSplit

  /////////////////////////////////////////
  void
  ChkChgAsymmetry(TCSlice& slc, Trajectory& tj, bool prt)
  {
    // looks for a high-charge point in the trajectory which may be due to the
    // trajectory crossing an interaction vertex. The properties of points on the opposite
    // sides of the high-charge point are analyzed. If significant differences are found, all points
    // near the high-charge point are removed as well as those from that point to the end
    if (!tcc.useAlg[kChkChgAsym]) return;
    if (tj.PDGCode == 111) return;
    unsigned short npts = tj.EndPt[1] - tj.EndPt[0];
    if (prt) mf::LogVerbatim("TC") << " Inside ChkChgAsymmetry T" << tj.ID;
    // ignore long tjs
    if (npts > 50) return;
    // ignore short tjs
    if (npts < 8) return;
    // require the charge pull > 5
    float bigPull = 5;
    unsigned short atPt = 0;
    // Don't consider the first/last few points in case there is a Bragg peak
    for (unsigned short ipt = tj.EndPt[0] + 2; ipt <= tj.EndPt[1] - 2; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.ChgPull > bigPull) {
        bigPull = tp.ChgPull;
        atPt = ipt;
      }
    } // ipt
    if (atPt == 0) return;
    // require that this point be near the DS end
    if ((atPt - tj.EndPt[0]) < 0.5 * npts) return;
    if (prt)
      mf::LogVerbatim("TC") << "CCA: T" << tj.ID << " Large Chg point at " << atPt
                            << ". Check charge asymmetry around it.";
    unsigned short nchk = 0;
    unsigned short npos = 0;
    unsigned short nneg = 0;
    for (short ii = 1; ii < 5; ++ii) {
      short iplu = atPt + ii;
      if (iplu > tj.EndPt[1]) break;
      short ineg = atPt - ii;
      if (ineg < tj.EndPt[0]) break;
      if (tj.Pts[iplu].Chg == 0) continue;
      if (tj.Pts[ineg].Chg == 0) continue;
      float asym = (tj.Pts[iplu].Chg - tj.Pts[ineg].Chg) / (tj.Pts[iplu].Chg + tj.Pts[ineg].Chg);
      ++nchk;
      if (asym > 0.5) ++npos;
      if (asym < -0.5) ++nneg;
      if (prt)
        mf::LogVerbatim("TC") << " ineg " << ineg << " iplu " << iplu << " asym " << asym
                              << " nchk " << nchk;
    } // ii
    if (nchk < 3) return;
    // require most of the points be very positive or very negative
    nchk -= 2;
    bool doTrim = (nneg > nchk) || (npos > nchk);
    if (!doTrim) return;
    // remove all the points at the end starting at the one just before the peak if the pull is not so good
    auto& prevTP = tj.Pts[atPt - 1];
    if (std::abs(prevTP.ChgPull) > 2) --atPt;
    for (unsigned short ipt = atPt; ipt <= tj.EndPt[1]; ++ipt)
      UnsetUsedHits(slc, tj.Pts[ipt]);
    SetEndPoints(tj);
    tj.AlgMod[kChkChgAsym] = true;
    if (prt) PrintTrajectory("CCA", slc, tj, USHRT_MAX);
  } // ChkChgAsymmetry

  ////////////////////////////////////////////////
  void
  TagShowerLike(TCSlice& slc, const CTP_t& inCTP)
  {
    // Tag Tjs as PDGCode = 11 if they have MCSMom < ShowerTag[1] and there are more than
    // ShowerTag[6] other Tjs with a separation < ShowerTag[2].

    if (tcc.showerTag[0] <= 0) return;
    if (slc.tjs.size() > 20000) return;
    float typicalChgRMS = 0.5 * (tcc.chargeCuts[1] + tcc.chargeCuts[2]);

    bool prt = (tcc.modes[kShowerTag] && inCTP == debug.CTP);
    unsigned int minClusterSize = tcc.showerTag[2];

    // clear out old tags and make a list of Tjs to consider
    std::vector<std::vector<int>> tjLists;
    std::vector<int> tjids;
    for (auto& tj : slc.tjs) {
      if (tj.CTP != inCTP) continue;
      if (tj.AlgMod[kKilled]) continue;
      if (tj.AlgMod[kHaloTj]) continue;
      if (tj.PDGCode == 11) tj.PDGCode = 0;
      // ignore Tjs with Bragg peaks
      bool skipit = false;
      for (unsigned short end = 0; end < 2; ++end)
        if (tj.EndFlag[end][kEndBragg]) skipit = true;
      if (skipit) continue;
      short npwc = NumPtsWithCharge(slc, tj, false);
      // Don't expect any (primary) electron to be reconstructed as a single trajectory for
      // more than ~2 radiation lengths ~ 30 cm for uB ~ 100 wires
      if (npwc > 100) continue;
      // allow short Tjs.
      if (npwc > 5) {
        // Increase the MCSMom cut if the Tj is long and the charge RMS is high to reduce sensitivity
        // to the fcl configuration. A primary electron may be reconstructed as one long Tj with large
        // charge rms and possibly high MCSMom or as several nearby shorter Tjs with lower charge rms
        float momCut = tcc.showerTag[0];
        if (tj.ChgRMS > typicalChgRMS) momCut *= tj.ChgRMS / typicalChgRMS;
        if (tj.MCSMom > momCut) continue;
      }
      tjids.push_back(tj.ID);
    } // tj
    if (tjids.size() < minClusterSize) return;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "TagShowerLike candidates";
      for (auto tid : tjids) myprt << " T" << tid;
    } // prt

    for (unsigned short it1 = 0; it1 < tjids.size() - 1; ++it1) {
      Trajectory& tj1 = slc.tjs[tjids[it1] - 1];
      for (unsigned short it2 = it1 + 1; it2 < tjids.size(); ++it2) {
        Trajectory& tj2 = slc.tjs[tjids[it2] - 1];
        unsigned short ipt1, ipt2;
        float doca = tcc.showerTag[1];
        // Find the separation between Tjs without considering dead wires
        TrajTrajDOCA(slc, tj1, tj2, ipt1, ipt2, doca, false);
        if (doca == tcc.showerTag[1]) continue;
        // make tighter cuts for user-defined short Tjs
        // found a close pair. See if one of these is in an existing cluster of Tjs
        bool inlist = false;
        for (unsigned short it = 0; it < tjLists.size(); ++it) {
          bool tj1InList =
            (std::find(tjLists[it].begin(), tjLists[it].end(), tj1.ID) != tjLists[it].end());
          bool tj2InList =
            (std::find(tjLists[it].begin(), tjLists[it].end(), tj2.ID) != tjLists[it].end());
          if (tj1InList || tj2InList) {
            // add the one that is not in the list
            if (!tj1InList) tjLists[it].push_back(tj1.ID);
            if (!tj2InList) tjLists[it].push_back(tj2.ID);
            inlist = true;
            break;
          }
          if (inlist) break;
        } // it
        // start a new list with this pair?
        if (!inlist) {
          std::vector<int> newlist(2);
          newlist[0] = tj1.ID;
          newlist[1] = tj2.ID;
          tjLists.push_back(newlist);
        }
      } // it2
    }   // it1
    if (tjLists.empty()) return;

    // mark them all as ShowerLike Tjs
    for (auto& tjl : tjLists) {
      // ignore small clusters
      if (tjl.size() < minClusterSize) continue;
      for (auto& tjID : tjl) {
        auto& tj = slc.tjs[tjID - 1];
        tj.PDGCode = 11;
      } // tjid
    }   // tjl

    if (prt) {
      unsigned short nsh = 0;
      for (auto& tjl : tjLists) {
        for (auto& tjID : tjl) {
          auto& tj = slc.tjs[tjID - 1];
          if (tj.PDGCode == 11) ++nsh;
        } // tjid
      }   // tjl
      mf::LogVerbatim("TC") << "TagShowerLike tagged " << nsh << " Tjs vertices in CTP " << inCTP;
    } // prt
  }   // TagShowerLike

} // namespace