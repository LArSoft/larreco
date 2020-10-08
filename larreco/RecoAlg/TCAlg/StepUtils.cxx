#include "larreco/RecoAlg/TCAlg/StepUtils.h"
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

  //////////////////////////////////////////
  void StepAway(TCSlice& slc, Trajectory& tj)
  {
    // Step along the direction specified in the traj vector in steps of size step
    // (wire spacing equivalents). Find hits between the last trajectory point and
    // the last trajectory point + step. A new trajectory point is added if hits are
    // found. Stepping continues until no signal is found for two consecutive steps
    // or until a wire or time boundary is reached.

    tj.IsGood = false;
    if(tj.Pts.empty()) return;

    unsigned short plane = DecodeCTP(tj.CTP).Plane;

    unsigned short lastPtWithUsedHits = tj.EndPt[1];

    unsigned short lastPt = lastPtWithUsedHits;
    // Construct a local TP from the last TP that will be moved on each step.
    // Only the Pos and Dir variables will be used
    TrajPoint ltp;
    ltp.CTP = tj.CTP;
    ltp.Pos = tj.Pts[lastPt].Pos;
    ltp.Dir = tj.Pts[lastPt].Dir;
    // A second TP is cloned from the leading TP of tj, updated with hits, fit
    // parameters,etc and possibly pushed onto tj as the next TP
    TrajPoint tp;

    // assume it is good from here on
    tj.IsGood = true;

    unsigned short nMissedSteps = 0;
    // Use MaxChi chisq cut for stiff trajectories
    bool useMaxChiCut = (tj.PDGCode == 13 || !tj.Strategy[kSlowing]);

    // Get the first forecast when there are 6 points with charge
    tjfs.resize(1);
    tjfs[0].nextForecastUpdate = 6;

    for(unsigned short step = 1; step < 10000; ++step) {
      unsigned short npwc = NumPtsWithCharge(slc, tj, false);
      // analyze the Tj when there are 6 points to see if we should stop
      if(npwc == 6 && StopShort(slc, tj, tcc.dbgStp)) break;
      // Get a forecast of what is ahead.
      bool getForecast = (tcc.doForecast && !tj.AlgMod[kRvPrp] && npwc == tjfs[tjfs.size() - 1].nextForecastUpdate);
      if(!getForecast && tj.Pts[tj.EndPt[1]].FitChi > 2) getForecast = true;
      if(getForecast) {
        Forecast(slc, tj);
        SetStrategy(slc, tj);
        SetPDGCode(slc, tj);
      }
      // make a copy of the previous TP
      lastPt = tj.Pts.size() - 1;
      tp = tj.Pts[lastPt];
      ++tp.Step;
      double stepSize = tcc.VLAStepSize;
      if(tp.AngleCode < 2) stepSize = std::abs(1/ltp.Dir[0]);
      // move the local TP position by one step in the right direction
      for(unsigned short iwt = 0; iwt < 2; ++iwt) ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;
      // copy this position into tp
      tp.Pos = ltp.Pos;
      tp.Dir = ltp.Dir;
      if(tcc.dbgStp) {
        mf::LogVerbatim myprt("TC");
        myprt<<"StepAway "<<step<<" Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" Dir "<<tp.Dir[0]<<" "<<tp.Dir[1]<<" stepSize "<<stepSize<<" AngCode "<<tp.AngleCode<<" Strategy";
        for(unsigned short ibt = 0; ibt < StrategyBitNames.size(); ++ibt) {
          if(tj.Strategy[ibt]) myprt<<" "<<StrategyBitNames[ibt];
        } // ib
      } // tcc.dbgStp
      // hit the boundary of the TPC?
      if(tp.Pos[0] < 0 || tp.Pos[0] > tcc.maxPos0[plane] ||
         tp.Pos[1] < 0 || tp.Pos[1] > tcc.maxPos1[plane]) break;
      // remove the old hits and other stuff
      tp.Hits.clear();
      tp.UseHit.reset();
      tp.FitChi = 0; tp.Chg = 0;
      tp.Environment.reset();
      unsigned int wire = std::nearbyint(tp.Pos[0]);
      if(!evt.goodWire[plane][wire]) tp.Environment[kEnvNotGoodWire] = true;
      // append to the trajectory
      tj.Pts.push_back(tp);
      // update the index of the last TP
      lastPt = tj.Pts.size() - 1;
      // look for hits
      bool sigOK = false;
      AddHits(slc, tj, lastPt, sigOK);
      // If successfull, AddHits has defined UseHit for this TP,
      // set the trajectory endpoints, and define HitPos.
      if(tj.Pts[lastPt].Hits.empty() && !tj.Pts[lastPt].Environment[kEnvNearSrcHit]) {
        // Require three points with charge on adjacent wires for small angle
        // stepping.
        if(tj.Pts[lastPt].AngleCode == 0 && lastPt == 2) return;
        // No close hits added.
        ++nMissedSteps;
        // First check for no signal in the vicinity. AddHits checks the hit collection for
        // the current slice. This version of SignalAtTp checks the allHits collection.
        double dir1 = ltp.Dir[1];
        ltp.Dir[1] = 0;
        sigOK = SignalAtTp(ltp);
        ltp.Dir[1] = dir1;
        if(lastPt > 0) {
          // break if this is a reverse propagate activity and there was no signal (not on a dead wire)
          if(!sigOK && tj.AlgMod[kRvPrp]) break;
          // Ensure that there is a signal here after missing a number of steps on a LA trajectory
          if(tj.Pts[lastPt].AngleCode > 0 && nMissedSteps > 4 && !sigOK) break;
          // the last point with hits (used or not) is the previous point
          unsigned short lastPtWithHits = lastPt - 1;
          float tps = TrajPointSeparation(tj.Pts[lastPtWithHits], ltp);
          float dwc = DeadWireCount(slc, ltp, tj.Pts[lastPtWithHits]);
          float nMissedWires = tps * std::abs(ltp.Dir[0]) - dwc;
          float maxWireSkip = tcc.maxWireSkipNoSignal;
          if(sigOK) maxWireSkip = tcc.maxWireSkipWithSignal;
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<" StepAway: no hits found at ltp "<<PrintPos(slc, ltp)<<" nMissedWires "<<std::fixed<<std::setprecision(1)<<nMissedWires<<" dead wire count "<<dwc<<" maxWireSkip "<<maxWireSkip<<" tj.PDGCode "<<tj.PDGCode;
          if(nMissedWires > maxWireSkip) {
            // We passed a number of wires without adding hits and are ready to quit.
            // First see if there is one good unused hit on the end TP and if so use it
            // lastPtWithHits + 1 == lastPt && tj.Pts[lastPtWithHits].Chg == 0 && tj.Pts[lastPtWithHits].Hits.size() == 1
            if(tj.EndPt[1] < tj.Pts.size() - 1 && tj.Pts[tj.EndPt[1]+1].Hits.size() == 1) {
              unsigned short lastLonelyPoint = tj.EndPt[1] + 1;
              unsigned int iht = tj.Pts[lastLonelyPoint].Hits[0];
              if(slc.slHits[iht].InTraj == 0 && tj.Pts[lastLonelyPoint].Delta < 3 * tj.Pts[lastLonelyPoint].DeltaRMS) {
                slc.slHits[iht].InTraj = tj.ID;
                tj.Pts[lastLonelyPoint].UseHit[0] = true;
                DefineHitPos(slc, tj.Pts[lastLonelyPoint]);
                SetEndPoints(tj);
                if(tcc.dbgStp) {
                  mf::LogVerbatim("TC")<<" Added a Last Lonely Hit before breaking ";
                  PrintTP("LLH", slc, lastPt, tj.StepDir, tj.Pass, tj.Pts[lastLonelyPoint]);
                }
              }
            }
            break;
          }
        } // lastPt > 0
        // no sense keeping this TP on tj if no hits were added
        tj.Pts.pop_back();
        continue;
      } // tj.Pts[lastPt].Hits.empty()
      // ensure that we actually moved
      if(lastPt > 0 && PosSep2(tj.Pts[lastPt].Pos, tj.Pts[lastPt-1].Pos) < 0.1) return;
      // Found hits at this location so reset the missed steps counter
      nMissedSteps = 0;
      // Update the last point fit, etc using the just added hit(s)
      UpdateTraj(slc, tj);
      // a failure occurred
      if(tj.NeedsUpdate) return;
      if(tj.Pts[lastPt].Chg == 0) {
        // There are points on the trajectory by none used in the last step. See
        // how long this has been going on
        float tps = TrajPointSeparation(tj.Pts[tj.EndPt[1]], ltp);
        float dwc = DeadWireCount(slc, ltp, tj.Pts[tj.EndPt[1]]);
        float nMissedWires = tps * std::abs(ltp.Dir[0]) - dwc;
        if(tcc.dbgStp)  mf::LogVerbatim("TC")<<" Hits exist on the trajectory but are not used. Missed wires "<<std::nearbyint(nMissedWires)<<" dead wire count "<<(int)dwc;
        // break if this is a reverse propagate activity with no dead wires
        if(tj.AlgMod[kRvPrp] && dwc == 0) break;
        if(nMissedWires > tcc.maxWireSkipWithSignal) break;
        // try this out
        if(!MaskedHitsOK(slc, tj)) {
          return;
        }
        // check for a series of bad fits and stop stepping
        if(tcc.useAlg[kStopBadFits] && nMissedWires > 4 && StopIfBadFits(slc, tj)) break;
        // Keep stepping
        if(tcc.dbgStp) {
          if(tj.AlgMod[kRvPrp]) {
            PrintTrajectory("RP", slc, tj, lastPt);
          } else {
            PrintTrajectory("SC", slc, tj, lastPt);
          }
        }
        continue;
      } // tp.Hits.empty()
      if(tj.Pts.size() == 3) {
        // ensure that the last hit added is in the same direction as the first two.
        // This is a simple way of doing it
        bool badTj = (PosSep2(tj.Pts[0].HitPos, tj.Pts[2].HitPos) < PosSep2(tj.Pts[0].HitPos, tj.Pts[1].HitPos));
        // ensure that this didn't start as a small angle trajectory and immediately turn
        // into a large angle one
        if(!badTj && tj.Pts[lastPt].AngleCode > tcc.maxAngleCode[tj.Pass]) badTj = true;
        // check for a large change in angle
        if(!badTj) {
          float dang = DeltaAngle(tj.Pts[0].Ang, tj.Pts[2].Ang);
          if(dang > 0.5) badTj = false;
        }
        //check for a wacky delta
        if(!badTj && tj.Pts[2].Delta > 2) badTj = true;
        if(badTj) {
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Bad Tj found on the third point. Quit stepping.";
          tj.IsGood = false;
          return;
        }
      } // tj.Pts.size() == 3
      // Update the local TP with the updated position and direction
      ltp.Pos = tj.Pts[lastPt].Pos;
      ltp.Dir = tj.Pts[lastPt].Dir;
      if(tj.MaskedLastTP) {
        // see if TPs have been masked off many times and if the
        // environment is clean. If so, return and try with next pass
        // cuts
        if(!MaskedHitsOK(slc, tj)) {
          if(tcc.dbgStp) {
            if(tj.AlgMod[kRvPrp]) {
              PrintTrajectory("RP", slc, tj, lastPt);
            } else {
              PrintTrajectory("SC", slc, tj, lastPt);
            }
          }
          return;
        }
        if(tcc.dbgStp) {
          if(tj.AlgMod[kRvPrp]) {
            PrintTrajectory("RP", slc, tj, lastPt);
          } else {
            PrintTrajectory("SC", slc, tj, lastPt);
          }
        }
        continue;
      }
      // We have added a TP with hits
      // check for a kink. Stop crawling if one is found
      GottaKink(slc, tj, true);
      if(tj.EndFlag[1][kEndKink]) {
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<"   stop at kink";
        break;
      }
      // See if the Chisq/DOF exceeds the maximum.
      // UpdateTraj should have reduced the number of points fit
      // as much as possible for this pass, so this trajectory is in trouble.
      if(tj.Pts[lastPt].FitChi > tcc.maxChi && useMaxChiCut) {
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<"   bad FitChi "<<tj.Pts[lastPt].FitChi<<" cut "<<tcc.maxChi;
        // remove the last point before quitting
        UnsetUsedHits(slc, tj.Pts[lastPt]);
        SetEndPoints(tj);
        tj.IsGood = (NumPtsWithCharge(slc, tj, true) > tcc.minPtsFit[tj.Pass]);
        break;
      }
      if(tcc.dbgStp) {
        if(tj.AlgMod[kRvPrp]) {
          PrintTrajectory("RP", slc, tj, lastPt);
        } else {
          PrintTrajectory("SC", slc, tj, lastPt);
        }
      } // tcc.dbgStp
    } // step

    SetPDGCode(slc, tj);

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"End StepAway with tj size "<<tj.Pts.size();

  } // StepAway

//////////////////////////////////////////
  bool StopShort(TCSlice& slc, Trajectory& tj, bool prt)
  {
    // Analyze the trajectory when it is short (~6 points) to look for a pattern like
    // this QQQqqq, where Q is a large charge hit and q is a low charge hit. If this
    // pattern is found, this function removes the last 3 points and returns true.
    // The calling function (e.g. StepAway) should decide what to do (e.g. stop stepping)

    // don't use this function during reverse propagation
    if(tj.AlgMod[kRvPrp]) return false;
    if(!tcc.useAlg[kStopShort]) return false;

    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    if(npwc > 10) return false;
    ParFit chgFit;
    FitPar(slc, tj, tj.EndPt[0], npwc, 1, chgFit, 1);
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"StopShort: chgFit at "<<PrintPos(slc, tj.Pts[tj.EndPt[0]]);
      myprt<<" ChiDOF "<<chgFit.ChiDOF;
      myprt<<" chg0 "<<chgFit.Par0<<" +/- "<<chgFit.ParErr;
      myprt<<" slp "<<chgFit.ParSlp<<" +/- "<<chgFit.ParSlpErr;
    } // prt
    // Look for a poor charge fit ChiDOF and a significant negative slope. These cuts
    // **should** prevent removing a genuine Bragg peak at the start, which should have
    // a good ChiDOF
    if(chgFit.ChiDOF < 2) return false;
    if(chgFit.ParSlp > -20) return false;
    if(prt) PrintTrajectory("SS", slc, tj, USHRT_MAX);
    // Find the average charge using the first 3 points
    float cnt = 0;
    float aveChg = 0;
    unsigned short lastHiPt = 0;
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      auto& tp = tj.Pts[ipt];
      if(tp.Chg <= 0) continue;
      aveChg += tp.Chg;
      ++cnt;
      lastHiPt = ipt;
      if(cnt == 3) break;
    } // tp
    if(cnt < 3) return false;
    aveChg /= cnt;
    if(prt) mf::LogVerbatim("TC")<<"  aveChg "<<(int)aveChg<<" last TP AveChg "<<(int)tj.Pts[tj.EndPt[1]].AveChg;
    // Look for a sudden drop in the charge (< 1/2)
    unsigned short firstLoPt = lastHiPt + 1;
    for(unsigned short ipt = lastHiPt + 1; ipt < tj.Pts.size(); ++ipt) {
      auto& tp = tj.Pts[ipt];
      if(tp.Chg <= 0 || tp.Chg > 0.5 * aveChg) continue;
      firstLoPt = ipt;
      break;
    } // ipt
    if(prt) mf::LogVerbatim("TC")<<"    stop tracking at "<<PrintPos(slc, tj.Pts[firstLoPt]);
    // Remove everything from the firstLoPt to the end of the trajectory
    for(unsigned short ipt = firstLoPt; ipt <= tj.EndPt[1]; ++ipt) UnsetUsedHits(slc, tj.Pts[ipt]);
    SetEndPoints(tj);
    UpdateTjChgProperties("SS", slc, tj, prt);
    tj.AlgMod[kStopShort] = true;
    return true;
  } // StopShort

//////////////////////////////////////////
  void SetStrategy(TCSlice& slc, Trajectory& tj)
  {
    // Determine if the tracking strategy is appropriate and make some tweaks if it isn't
    if(tjfs.empty()) return;
    // analyze the last forecast
    auto& tjf = tjfs[tjfs.size() - 1];

    auto& lastTP = tj.Pts[tj.EndPt[1]];

    if(!tcc.useAlg[kNewCuts]) {
      // Stay in Slowing strategy if we are in it and keep the number of points fit constant
      if(tj.Strategy[kSlowing]) {
        lastTP.NTPsFit = 5;
        return;
      }
    } // !tcc.useAlg[kNewCuts]

    float npwc = NumPtsWithCharge(slc, tj, false);
    // Keep using the StiffMu strategy if the tj is long and MCSMom is high
    if(tj.Strategy[kStiffMu] && tj.MCSMom > 800 && npwc > 200) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: Keep using the StiffMu strategy";
      return;
    }
    bool tkLike = (tjf.outlook < 1.5);
    // A showering-electron-like trajectory
    bool chgIncreasing = (tjf.chgSlope > 0);
    // A showering-electron-like trajectory
    bool shLike = (tjf.outlook > 3 && chgIncreasing && !tkLike);
    if(!shLike) shLike = tjf.showerLikeFraction > 0.5;
    float momRat = 0;
    if(tj.MCSMom > 0) momRat = (float)tjf.MCSMom / (float)tj.MCSMom;
    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<"SetStrategy: npwc "<<npwc<<" outlook "<<tjf.outlook;
      myprt<<" tj MCSMom "<<tj.MCSMom<<" forecast MCSMom "<<tjf.MCSMom;
      myprt<<" momRat "<<std::fixed<<std::setprecision(2)<<momRat;
      myprt<<" tkLike? "<<tkLike<<" shLike? "<<shLike;
      myprt<<" chgIncreasing? "<<chgIncreasing;
      myprt<<" leavesBeforeEnd? "<<tjf.leavesBeforeEnd<<" endBraggPeak? "<<tjf.endBraggPeak;
      myprt<<" nextForecastUpdate "<<tjf.nextForecastUpdate;
    }
    if(tjf.outlook < 0) return;
    // Look for a long clean muon in the forecast
    bool stiffMu = (tkLike && tjf.MCSMom > 600 && tjf.nextForecastUpdate > 100 && !tjf.endBraggPeak);
    if(stiffMu) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: High MCSMom, long forecast. Use the StiffMu strategy";
      tj.Strategy.reset();
      tj.Strategy[kStiffMu] = true;
      return;
    } // StiffMu
    bool notStiff = (!tj.Strategy[kStiffEl] && !tj.Strategy[kStiffMu]);
    if(notStiff && !shLike && tj.MCSMom < 100 && tjf.MCSMom < 100 && tjf.leavesBeforeEnd && npwc < 100) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: Low MCSMom. Use the Slowing strategy";
      tj.Strategy.reset();
      tj.Strategy[kSlowing] = true;
      lastTP.NTPsFit = 5;
      return;
    } // Low MCSMom
    if(notStiff && !shLike && tj.MCSMom < 200 && momRat < 0.7 && chgIncreasing) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: Low MCSMom & low momRat. Use the Slowing strategy";
      tj.Strategy.reset();
      tj.Strategy[kSlowing] = true;
      lastTP.NTPsFit = 5;
      return;
    } // low MCSMom
    if(tjf.endBraggPeak) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: Found a Bragg peak. Use the Slowing strategy";
      tj.Strategy.reset();
      tj.Strategy[kSlowing] = true;
      lastTP.NTPsFit = 5;
      return;
    } // tracklike with Bragg peak
    // A short forecast window. No Bragg-peak (which would have been found above) but increasing charge
    if(!tjf.leavesBeforeEnd && chgIncreasing && !shLike && tjf.nextForecastUpdate < 30) {
      tj.Strategy.reset();
      tj.Strategy[kSlowing] = true;
      lastTP.NTPsFit = 5;
      return;
    }
    if(tkLike && tjf.nextForecastUpdate > 100 && tjf.leavesBeforeEnd && tjf.MCSMom < 500) {
      // A long track-like trajectory that has many points fit and the outlook is track-like and
      // it leaves the forecast polygon. Don't change the strategy but decrease the number of points fit
      lastTP.NTPsFit /= 2;
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: Long track-like wandered out of forecast envelope. Reduce NTPsFit to "<<lastTP.NTPsFit;
      return;
    } // fairly long and leaves the side
    // a track-like trajectory that has high MCSMom in the forecast and hits a shower
    if(tkLike && tjf.MCSMom > 600 && (tjf.foundShower || tjf.chgFitChiDOF > 20)) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: high MCSMom "<<tjf.MCSMom<<" and a shower ahead. Use the StiffEl strategy";
      tj.Strategy.reset();
      tj.Strategy[kStiffEl] = true;
      // we think we know the direction (towards the shower) so  startEnd is 0
      tj.StartEnd = 0;
      return;
    } // Stiff electron
    if(shLike && !tjf.leavesBeforeEnd) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"SetStrategy: Inside a shower. Use the StiffEl strategy";
      tj.Strategy.reset();
      tj.Strategy[kStiffEl] = true;
      // we think we know the direction (towards the shower) so  startEnd is 0
      tj.StartEnd = 0;
      return;
    } // StiffEl
    if(!tcc.useAlg[kNewCuts]) {
      tj.Strategy.reset();
      tj.Strategy[kNormal] = true;
    } // !tcc.useAlg[kNewCuts]
    // set to normal

  } // SetStrategy

  //////////////////////////////////////////
  void Forecast(TCSlice& slc, const Trajectory& tj)
  {
    // Extrapolate the last TP of tj by many steps and return a forecast of what is ahead
    // -1       error or not sure
    // ~1       track-like with a slight chance of showers
    // >2       shower-like
    // nextForecastUpdate is set to the number of points on the trajectory when this function should
    // be called again

    if(tj.Pts[tj.EndPt[1]].AngleCode == 2) return;

    // add a new forecast
    tjfs.resize(tjfs.size() + 1);
    // assume there is insufficient info to make a decision
    auto& tjf = tjfs[tjfs.size() - 1];
    tjf.outlook = -1;
    tjf.nextForecastUpdate = USHRT_MAX;

    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    unsigned short istp = 0;
    unsigned short nMissed = 0;

    bool doPrt = tcc.dbgStp;
    // turn off annoying output from DefineHitPos
    if(doPrt) tcc.dbgStp = false;
    // find the minimum average TP charge. This will be used to calculate the
    // 'effective number of hits' on a wire = total charge on the wire within the
    // window / (minimum average TP charge). This is intended to reduce the sensitivity
    // of this metric to the hit finder configuration
    float minAveChg = 1E6;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      if(tj.Pts[ipt].AveChg <= 0) continue;
      if(tj.Pts[ipt].AveChg > minAveChg) continue;
      minAveChg = tj.Pts[ipt].AveChg;
    } // ipt
    if(minAveChg <= 0 || minAveChg == 1E6) return;
    // start a forecast Tj comprised of the points in the forecast envelope
    Trajectory fctj;
    fctj.CTP = tj.CTP;
    fctj.ID = evt.WorkID;
    // make a local copy of the last point
    auto ltp = tj.Pts[tj.EndPt[1]];
    // Use the hits position instead of the fitted position so that a bad
    // fit doesn't screw up the forecast.
    float forecastWin0 = std::abs(ltp.Pos[1] - ltp.HitPos[1]);
    if(forecastWin0 < 1) forecastWin0 = 1;
    ltp.Pos = ltp.HitPos;
    double stepSize = std::abs(1/ltp.Dir[0]);
    float window = tcc.showerTag[7] * stepSize;
    if(doPrt) {
      mf::LogVerbatim("TC")<<"Forecast: T"<<tj.ID<<" PDGCode "<<tj.PDGCode<<" npwc "<<npwc
          <<" minAveChg "<<(int)minAveChg<<" stepSize "<<std::setprecision(2)<<stepSize<<" window "<<window;
      mf::LogVerbatim("TC")<<" stp ___Pos____  nTPH  Chg ChgPull  Delta    DRMS  chgWid nTkLk nShLk";
    }
    unsigned short plane = DecodeCTP(ltp.CTP).Plane;
    float totHits = 0;
    fctj.TotChg = 0;
    float maxChg = 0;
    unsigned short maxChgPt = 0;
    unsigned short leavesNear = USHRT_MAX;
    bool leavesBeforeEnd = false;
    unsigned short showerStartNear = USHRT_MAX;
    unsigned short showerEndNear = USHRT_MAX;
    float nShLike = 0;
    float nTkLike = 0;
    unsigned short trimPts = 0;
    for(istp = 0; istp < 1000; ++istp) {
      // move the local TP position by one step in the right direction
      for(unsigned short iwt = 0; iwt < 2; ++iwt) ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;
      unsigned int wire = std::nearbyint(ltp.Pos[0]);
      if(wire < slc.firstWire[plane]) break;
      if(wire > slc.lastWire[plane]-1) break;
      MoveTPToWire(ltp, (float)wire);
      ++ltp.Step;
      ltp.Environment.reset();
      if(FindCloseHits(slc, ltp, window, kAllHits)) {
        // Found hits or the wire is dead
        // set all hits used so that we can use DefineHitPos. Note that
        // the hit InTraj is not used or tested in DefineHitPos so this doesn't
        // screw up any assns
        if(!ltp.Environment[kEnvNotGoodWire]) {
          nMissed = 0;
          ltp.UseHit.set();
          DefineHitPos(slc, ltp);
          fctj.TotChg += ltp.Chg;
          ltp.Delta = PointTrajDOCA(slc, ltp.HitPos[0], ltp.HitPos[1], ltp);
          ltp.DeltaRMS = ltp.Delta / window;
          ltp.Environment.reset();
          totHits += ltp.Hits.size();
          if(ltp.Chg > maxChg) {
            maxChg = ltp.Chg;
            maxChgPt = fctj.Pts.size();
          }
          // Note that ChgPull uses the average charge and charge RMS of the
          // trajectory before it entered the forecast envelope
          ltp.ChgPull = (ltp.Chg / minAveChg - 1) / tj.ChgRMS;
          if(ltp.ChgPull > 3 && ltp.Hits.size() > 2) {
            ++nShLike;
            // break if it approaches the side of the envelope
            ltp.Environment[kEnvNearShower] = true;
            // flag a showerlike TP so it isn't used in the MCSMom calculation
            ltp.HitPosErr2 = 100;
          } else {
            ++nTkLike;
            ltp.Environment[kEnvNearShower] = false;
          }
          if(fctj.Pts.size() > 10) {
            float shFrac = nShLike / (nShLike + nTkLike);
            if(shFrac > 0.5) {
              if(doPrt) mf::LogVerbatim("TC")<<"Getting showerlike - break";
              break;
            }
          } // fctj.Pts.size() > 6
          // break if it approaches the side of the envelope
          if(ltp.DeltaRMS > 0.8) {
            leavesNear = npwc + fctj.Pts.size();
            if(doPrt) mf::LogVerbatim("TC")<<"leaves before end - break";
            leavesBeforeEnd = true;
            break;
          }
          fctj.Pts.push_back(ltp);
          if(doPrt) {
            mf::LogVerbatim myprt("TC");
            myprt<<std::setw(4)<<npwc + fctj.Pts.size()<<" "<<PrintPos(slc, ltp);
            myprt<<std::setw(5)<<ltp.Hits.size();
            myprt<<std::setw(5)<<(int)ltp.Chg;
            myprt<<std::fixed<<std::setprecision(1);
            myprt<<std::setw(8)<<ltp.ChgPull;
            myprt<<std::setw(8)<<ltp.Delta;
            myprt<<std::setw(8)<<std::setprecision(2)<<ltp.DeltaRMS;
            myprt<<std::setw(8)<<sqrt(ltp.HitPosErr2);
            myprt<<std::setw(6)<<(int)nTkLike;
            myprt<<std::setw(6)<<(int)nShLike;
          } // doPrt
        }
      } else {
        // no hits found
        ++nMissed;
        if(nMissed == 10) {
          if(doPrt) mf::LogVerbatim("TC")<<"No hits found after 10 steps - break";
          break;
        }
      } // no hits found
    } // istp
    // not enuf info to make a forecast
    tcc.dbgStp = doPrt;
    if(fctj.Pts.size() < 3) return;
    // truncate and re-calculate totChg?
    if(trimPts > 0) {
      // truncate the forecast trajectory
      fctj.Pts.resize(fctj.Pts.size() - trimPts);
      // recalculate the total charge
      fctj.TotChg = 0;
      for(auto& tp : fctj.Pts) fctj.TotChg += tp.Chg;
    } // showerEndNear != USHRT_MAX
    SetEndPoints(fctj);
    fctj.MCSMom = MCSMom(slc, fctj);
    tjf.MCSMom = fctj.MCSMom;
    ParFit chgFit;
    if(maxChgPt > 0.3 * fctj.Pts.size() && maxChg > 3 * tj.AveChg) {
      // find the charge slope from the beginning to the maxChgPt if it is well away
      // from the start and it is very large
      FitPar(slc, fctj, 0, maxChgPt, 1, chgFit, 1);
    } else {
      // BB 6/13/2020 Don't use the last point
      FitPar(slc, fctj, 0, fctj.Pts.size()-1, 1, chgFit, 1);
    }
    tjf.chgSlope = chgFit.ParSlp;
    tjf.chgSlopeErr = chgFit.ParSlpErr;
    tjf.chgFitChiDOF = chgFit.ChiDOF;
    ChkStop(slc, fctj);
    UpdateTjChgProperties("Fc", slc, fctj, false);
    tjf.chgRMS = fctj.ChgRMS;
    tjf.endBraggPeak = fctj.EndFlag[1][kEndBragg];
    // Set outlook = Estimate of the number of hits per wire
    tjf.outlook = fctj.TotChg / (fctj.Pts.size() * tj.AveChg);
    // assume we got to the end
    tjf.nextForecastUpdate = npwc + fctj.Pts.size();
    tjf.leavesBeforeEnd = leavesBeforeEnd;
    tjf.foundShower = false;
    if(leavesNear < tjf.nextForecastUpdate) {
      // left the side
      tjf.nextForecastUpdate = leavesNear;
    } else if(showerStartNear < tjf.nextForecastUpdate) {
      // found a shower start
      tjf.nextForecastUpdate = showerStartNear;
      tjf.foundShower = true;
    } else if(showerEndNear < tjf.nextForecastUpdate) {
      // found a shower end
      tjf.nextForecastUpdate = showerEndNear;
    }
    nShLike = 0;
    for(auto& tp : fctj.Pts) if(tp.Environment[kEnvNearShower]) ++nShLike;
    tjf.showerLikeFraction = (float)nShLike / (float)fctj.Pts.size();

    if(doPrt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Forecast T"<<tj.ID<<" tj.AveChg "<<(int)tj.AveChg;
      myprt<<" start "<<PrintPos(slc, tj.Pts[tj.EndPt[1]])<<" cnt "<<fctj.Pts.size()<<" totChg "<<(int)fctj.TotChg;
      myprt<<" last pos "<<PrintPos(slc, ltp);
      myprt<<" MCSMom "<<tjf.MCSMom;
      myprt<<" outlook "<<std::fixed<<std::setprecision(2)<<tjf.outlook;
      myprt<<" chgSlope "<<std::setprecision(1)<<tjf.chgSlope<<" +/- "<<tjf.chgSlopeErr;
      myprt<<" chgRMS "<<std::setprecision(1)<<tjf.chgRMS;
      myprt<<" endBraggPeak "<<tjf.endBraggPeak;
      myprt<<" chiDOF "<<tjf.chgFitChiDOF;
      myprt<<" showerLikeFraction "<<tjf.showerLikeFraction;
      myprt<<" nextForecastUpdate "<<tjf.nextForecastUpdate;
      myprt<<" leavesBeforeEnd? "<<tjf.leavesBeforeEnd;
      myprt<<" foundShower? "<<tjf.foundShower;
    }

  } // Forecast

  //////////////////////////////////////////
  void UpdateStiffEl(TCSlice& slc, Trajectory& tj)
  {
    // A different stategy for updating a high energy electron trajectories
    if(!tj.Strategy[kStiffEl]) return;
    TrajPoint& lastTP = tj.Pts[tj.EndPt[1]];
    // Set the lastPT delta before doing the fit
    lastTP.Delta = PointTrajDOCA(slc, lastTP.HitPos[0], lastTP.HitPos[1], lastTP);
    if(tj.Pts.size() < 30) lastTP.NTPsFit += 1;
    FitTraj(slc, tj);
    UpdateTjChgProperties("UET", slc, tj, tcc.dbgStp);
    UpdateDeltaRMS(slc, tj);
    tj.MCSMom = MCSMom(slc, tj);
    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"UpdateStiffEl: lastPt "<<tj.EndPt[1]<<" Delta "<<lastTP.Delta<<" AngleCode "<<lastTP.AngleCode<<" FitChi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit<<" MCSMom "<<tj.MCSMom;
    }
    tj.NeedsUpdate = false;
    tj.PDGCode = 111;
    return;
  } // UpdateStiffTj

  //////////////////////////////////////////
  void UpdateTraj(TCSlice& slc, Trajectory& tj)
  {
    // Updates the last added trajectory point fit, average hit rms, etc.

    tj.NeedsUpdate = true;
    tj.MaskedLastTP = false;

    if(tj.EndPt[1] < 1) return;

    if(tj.Strategy[kStiffEl]) {
      UpdateStiffEl(slc, tj);
      return;
    }
    unsigned int lastPt = tj.EndPt[1];
    TrajPoint& lastTP = tj.Pts[lastPt];
    // nothing needs to be done if the last point has no hits but is near a hit in the
    // srcHit collection
    if(lastTP.Hits.empty() && lastTP.Environment[kEnvNearSrcHit]) {
      tj.NeedsUpdate = false;
      return;
    }
    unsigned short npwc = NumPtsWithCharge(slc, tj, false);

    // find the previous TP that has hits (and was therefore in the fit)
    unsigned short prevPtWithHits = USHRT_MAX;
    unsigned short firstFitPt = tj.EndPt[0];
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg > 0) {
        prevPtWithHits = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(prevPtWithHits == USHRT_MAX) return;

    // define the FitChi threshold above which something will be done
    float maxChi = 2;
    unsigned short minPtsFit = tcc.minPtsFit[tj.Pass];
    // just starting out?
    if(lastPt < 4) minPtsFit = 2;
    bool cleanMuon = (tj.PDGCode == 13 && TrajIsClean(slc, tj, tcc.dbgStp) && !tj.Strategy[kSlowing]);
    // was !TrajIsClean...
    if(cleanMuon) {
      // Fitting a clean muon
      maxChi = tcc.maxChi;
      minPtsFit = lastPt / 3;
    }

    // Set the lastPT delta before doing the fit
    lastTP.Delta = PointTrajDOCA(slc, lastTP.HitPos[0], lastTP.HitPos[1], lastTP);

    // update MCSMom. First ensure that nothing bad has happened
    if(npwc > 3 && tj.Pts[lastPt].Chg > 0 && !tj.Strategy[kSlowing]) {
      short newMCSMom = MCSMom(slc, tj);
      short minMCSMom = 0.5 * tj.MCSMom;
      if(lastPt > 10 && newMCSMom < minMCSMom) {
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<"UT: MCSMom took a nose-dive "<<newMCSMom;
        UnsetUsedHits(slc, lastTP);
        DefineHitPos(slc, lastTP);
        SetEndPoints(tj);
        tj.NeedsUpdate = false;
        return;
      }
      tj.MCSMom = newMCSMom;
    } // npwc > 3

    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"UT: lastPt "<<lastPt<<" lastTP.Delta "<<lastTP.Delta<<" previous point with hits "<<prevPtWithHits<<" tj.Pts size "<<tj.Pts.size()<<" AngleCode "<<lastTP.AngleCode<<" PDGCode "<<tj.PDGCode<<" maxChi "<<maxChi<<" minPtsFit "<<minPtsFit<<" MCSMom "<<tj.MCSMom;
    }

    UpdateTjChgProperties("UT", slc, tj, tcc.dbgStp);

    if(lastPt == 1) {
      // Handle the second trajectory point. No error calculation. Just update
      // the position and direction
      lastTP.NTPsFit = 2;
      FitTraj(slc, tj);
      lastTP.FitChi = 0.01;
      lastTP.AngErr = tj.Pts[0].AngErr;
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"UT: Second traj point pos "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1];
      tj.NeedsUpdate = false;
      SetAngleCode(lastTP);
      return;
    }

    if(lastPt == 2) {
      // Third trajectory point. Keep it simple
      lastTP.NTPsFit = 3;
      FitTraj(slc, tj);
      tj.NeedsUpdate = false;
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"UT: Third traj point fit "<<lastTP.FitChi;
      SetAngleCode(lastTP);
      return;
    }

    // Fit with > 2 TPs
    // Keep adding hits until Chi/DOF exceeds 1
    if(tj.Pts[prevPtWithHits].FitChi < 1 && !tj.Strategy[kSlowing]) lastTP.NTPsFit += 1;
    // Reduce the number of points fit if the trajectory is long and chisq is getting a bit larger
    if(lastPt > 20 && tj.Pts[prevPtWithHits].FitChi > 1.5 && lastTP.NTPsFit > minPtsFit) lastTP.NTPsFit -= 2;
    // don't let long muon fits get too long
    if(cleanMuon && lastPt > 200 && tj.Pts[prevPtWithHits].FitChi > 1.0) lastTP.NTPsFit -= 2;
    FitTraj(slc, tj);

    // don't get too fancy when we are starting out
    if(lastPt < 6) {
      tj.NeedsUpdate = false;
      UpdateDeltaRMS(slc, tj);
      SetAngleCode(lastTP);
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Return with lastTP.FitChi "<<lastTP.FitChi<<" Chg "<<lastTP.Chg;
      return;
    }

    // find the first point that was fit.
    unsigned short cnt = 0;
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tj.Pts[ipt].Chg > 0) {
        firstFitPt = ipt;
        ++cnt;
      }
      if(cnt == lastTP.NTPsFit) break;
      if(ipt == 0) break;
    }

    unsigned short ndead = DeadWireCount(slc, lastTP.HitPos[0], tj.Pts[firstFitPt].HitPos[0], tj.CTP);
    if(lastTP.FitChi > 1.5 && tj.Pts.size() > 6) {
      // A large chisq jump can occur if we just jumped a large block of dead wires. In
      // this case we don't want to mask off the last TP but reduce the number of fitted points
      // This count will be off if there a lot of dead or missing wires...
      // reduce the number of points significantly
      if(ndead > 5 && !cleanMuon) {
        if(lastTP.NTPsFit > 5) lastTP.NTPsFit = 5;
      } else {
        // Have a longish trajectory and chisq was a bit large.
        // Was this a sudden occurrence and the fraction of TPs are included
        // in the fit? If so, we should mask off this
        // TP and keep going. If these conditions aren't met, we
        // should reduce the number of fitted points
        float chirat = 0;
        if(prevPtWithHits != USHRT_MAX && tj.Pts[prevPtWithHits].FitChi > 0) chirat = lastTP.FitChi / tj.Pts[prevPtWithHits].FitChi;
        // Don't mask hits when doing RevProp. Reduce NTPSFit instead
        tj.MaskedLastTP = (chirat > 1.5 && lastTP.NTPsFit > 0.3 * NumPtsWithCharge(slc, tj, false) && !tj.AlgMod[kRvPrp]);
        // BB Don't mask TPs on low MCSMom Tjs
        if(tj.MaskedLastTP && tj.MCSMom < 30) tj.MaskedLastTP = false;
        if(tcc.dbgStp) {
          mf::LogVerbatim("TC")<<" First fit chisq too large "<<lastTP.FitChi<<" prevPtWithHits chisq "<<tj.Pts[prevPtWithHits].FitChi<<" chirat "<<chirat<<" NumPtsWithCharge "<<NumPtsWithCharge(slc, tj, false)<<" tj.MaskedLastTP "<<tj.MaskedLastTP;
        }
        // we should also mask off the last TP if there aren't enough hits
        // to satisfy the minPtsFit constraint
        if(!tj.MaskedLastTP && NumPtsWithCharge(slc, tj, true) < minPtsFit) tj.MaskedLastTP = true;
      } // few dead wires
    } // lastTP.FitChi > 2 ...

    // Deal with a really long trajectory that is in trouble (uB cosmic).
    if(tj.PDGCode == 13 && lastTP.FitChi > tcc.maxChi) {
      if(lastTP.NTPsFit > 1.3 * tcc.muonTag[0]) {
        lastTP.NTPsFit *= 0.8;
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Muon - Reduce NTPsFit "<<lastPt;
      } else {
        tj.MaskedLastTP = true;
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Muon - mask last point "<<lastPt;
      }
    }

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"UT: First fit "<<lastTP.Pos[0]<<" "<<lastTP.Pos[1]<<"  dir "<<lastTP.Dir[0]<<" "<<lastTP.Dir[1]<<" FitChi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit<<" ndead wires "<<ndead<<" tj.MaskedLastTP "<<tj.MaskedLastTP;
      if(tj.MaskedLastTP) {
        UnsetUsedHits(slc, lastTP);
        DefineHitPos(slc, lastTP);
        SetEndPoints(tj);
        lastPt = tj.EndPt[1];
        lastTP.NTPsFit -= 1;
        FitTraj(slc, tj);
        UpdateTjChgProperties("UT", slc, tj, tcc.dbgStp);
        SetAngleCode(lastTP);
        return;
      }  else {
        // a more gradual change in chisq. Maybe reduce the number of points
        unsigned short newNTPSFit = lastTP.NTPsFit;
        // reduce the number of points fit to keep Chisq/DOF < 2 adhering to the pass constraint
        // and also a minimum number of points fit requirement for long muons
        float prevChi = lastTP.FitChi;
        unsigned short ntry = 0;
        float chiCut = 1.5;
        if(tj.Strategy[kStiffMu]) chiCut = 5;
        while(lastTP.FitChi > chiCut && lastTP.NTPsFit > minPtsFit) {
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
          // BB Try to add a last lonely hit on a low MCSMom tj on the last try
          if(newNTPSFit == minPtsFit && tj.MCSMom < 30) chiCut = 2;
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<"  Bad FitChi "<<lastTP.FitChi<<" Reduced NTPsFit to "<<lastTP.NTPsFit<<" Pass "<<tj.Pass<<" chiCut "<<chiCut;
          FitTraj(slc, tj);
          tj.NeedsUpdate = true;
          if(lastTP.FitChi > prevChi) {
            if(tcc.dbgStp) mf::LogVerbatim("TC")<<"  Chisq is increasing "<<lastTP.FitChi<<"  Try to remove an earlier bad hit";
            MaskBadTPs(slc, tj, chiCut);
            ++ntry;
            if(ntry == 2) break;
          }
          prevChi = lastTP.FitChi;
          if(lastTP.NTPsFit == minPtsFit) break;
        } // lastTP.FitChi > 2 && lastTP.NTPsFit > 2
      }
      // last ditch attempt if things look bad. Drop the last hit
      if(tj.Pts.size() > tcc.minPtsFit[tj.Pass] && lastTP.FitChi > maxChi) {
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<"  Last try. Drop last TP "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;
        UnsetUsedHits(slc, lastTP);
        DefineHitPos(slc, lastTP);
        SetEndPoints(tj);
        lastPt = tj.EndPt[1];
        FitTraj(slc, tj);
        tj.MaskedLastTP = true;
      }

    if(tj.NeedsUpdate) UpdateTjChgProperties("UT", slc, tj, tcc.dbgStp);

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"  Fit done. Chi "<<lastTP.FitChi<<" NTPsFit "<<lastTP.NTPsFit;

    if(tj.EndPt[0] == tj.EndPt[1]) return;

    // Don't let the angle error get too small too soon. Stepping would stop if the first
    // few hits on a low momentum wandering track happen to have a very good fit to a straight line.
    // We will do this by averaging the default starting value of AngErr of the first TP with the current
    // value from FitTraj.
    if(lastPt < 14) {
      float defFrac = 1 / (float)(tj.EndPt[1]);
      lastTP.AngErr = defFrac * tj.Pts[0].AngErr + (1 - defFrac) * lastTP.AngErr;
    }

    UpdateDeltaRMS(slc, tj);
    SetAngleCode(lastTP);

    tj.NeedsUpdate = false;
    return;

  } // UpdateTraj

  ////////////////////////////////////////////////
  void CheckStiffEl(TCSlice& slc, Trajectory& tj)
  {
    if(!tj.Strategy[kStiffEl]) return;
    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"inside CheckStiffTj with NumPtsWithCharge = "<<NumPtsWithCharge(slc, tj, false);
    }
    // Fill in any gaps with hits that were skipped, most likely delta rays on muon tracks
    FillGaps(slc, tj);
    // Update the trajectory parameters at the beginning of the trajectory
    ChkBegin(slc, tj);
  } // CheckStiffTj

  ////////////////////////////////////////////////
  void AddHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, bool& sigOK)
  {
    // Try to add hits to the trajectory point ipt on the supplied
    // trajectory

    // assume failure
    sigOK = false;

    if(tj.Pts.empty()) return;
    if(ipt > tj.Pts.size() - 1) return;

    // Call large angle hit finding if the last tp is large angle
    if(tj.Pts[ipt].AngleCode == 2) {
      AddLAHits(slc, tj, ipt, sigOK);
      return;
    }

    TrajPoint& tp = tj.Pts[ipt];
    std::vector<unsigned int> closeHits;
    unsigned int lastPtWithUsedHits = tj.EndPt[1];

    unsigned short plane = DecodeCTP(tj.CTP).Plane;
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if(wire < slc.firstWire[plane] || wire > slc.lastWire[plane]-1) return;
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

    // The delta cut shouldn't be less than the delta of hits added on the previous step
    float minDeltaCut = 1.1 * tj.Pts[lastPtWithUsedHits].Delta;
    if(deltaCut < minDeltaCut) deltaCut = minDeltaCut;

    deltaCut *= tcc.projectionErrFactor;
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" AddHits: calculated deltaCut "<<deltaCut<<" dw "<<dw<<" dpos "<<dpos;

    if(deltaCut < 0.5) deltaCut = 0.5;
    if(deltaCut > 3) deltaCut = 3;

    // TY: open it up for RevProp, since we might be following a stopping track
    if(tj.AlgMod[kRvPrp]) deltaCut *= 2;

    // loosen up a bit if we just passed a block of dead wires
    bool passedDeadWires = (abs(dw) > 20 && DeadWireCount(slc, tp.Pos[0], tj.Pts[lastPtWithUsedHits].Pos[0], tj.CTP) > 10);
    if(passedDeadWires) deltaCut *= 2;
    // open it up for StiffEl and Slowing strategies
    if(tj.Strategy[kStiffEl] || tj.Strategy[kSlowing]) deltaCut = 3;

    // Create a larger cut to use in case there is nothing close
    float bigDelta = 2 * deltaCut;
    unsigned int imBig = UINT_MAX;
    tp.Delta = deltaCut;
    // ignore all hits with delta larger than maxDeltaCut
    float maxDeltaCut = 2 * bigDelta;
    // apply some limits
    if(!passedDeadWires && maxDeltaCut > 3) {
      maxDeltaCut = 3;
      bigDelta = 1.5;
    }

    // projected time in ticks for testing the existence of a hit signal
    raw::TDCtick_t rawProjTick = (float)(tp.Pos[1] / tcc.unitsPerTick);
    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<" AddHits: wire "<<wire<<" tp.Pos[0] "<<tp.Pos[0]<<" projTick "<<rawProjTick<<" deltaRMS "<<tp.DeltaRMS<<" tp.Dir[0] "<<tp.Dir[0]<<" deltaCut "<<deltaCut<<" dpos "<<dpos<<" projErr "<<projErr<<" ExpectedHitsRMS "<<ExpectedHitsRMS(slc, tp);
    }

    std::vector<unsigned int> hitsInMultiplet;

    geo::PlaneID planeID = DecodeCTP(tj.CTP);
    unsigned int ipl = planeID.Plane;
    if(wire > slc.lastWire[ipl]) return;
    // Assume a signal exists on a dead wire
    if(!evt.goodWire[ipl][wire]) sigOK = true;
    if(slc.wireHitRange[ipl][wire].first == UINT_MAX) return;
    unsigned int firstHit = slc.wireHitRange[ipl][wire].first;
    unsigned int lastHit = slc.wireHitRange[ipl][wire].second;
    float fwire = wire;
    for(unsigned int iht = firstHit; iht <= lastHit; ++iht) {
      if(slc.slHits[iht].InTraj == tj.ID) continue;
      if(slc.slHits[iht].InTraj == SHRT_MAX) continue;
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      if(rawProjTick > hit.StartTick() && rawProjTick < hit.EndTick()) sigOK = true;
      float ftime = tcc.unitsPerTick * hit.PeakTime();
      float delta = PointTrajDOCA(slc, fwire, ftime, tp);
      // increase the delta cut if this is a long pulse hit
      bool longPulseHit = LongPulseHit(hit);
      if(longPulseHit) {
        if(delta > 3) continue;
      } else {
        if(delta > maxDeltaCut) continue;
      }
      float dt = std::abs(ftime - tp.Pos[1]);
      GetHitMultiplet(slc, iht, hitsInMultiplet, false);
      if(tcc.dbgStp && delta < 100 && dt < 100) {
        mf::LogVerbatim myprt("TC");
        myprt<<"  iht "<<iht;
        myprt<<" "<<PrintHit(slc.slHits[iht]);
        myprt<<" delta "<<std::fixed<<std::setprecision(2)<<delta<<" deltaCut "<<deltaCut<<" dt "<<dt;
        myprt<<" BB Mult "<<hitsInMultiplet.size()<<" RMS "<<std::setprecision(1)<<hit.RMS();
        myprt<<" Chi "<<std::setprecision(1)<<hit.GoodnessOfFit();
        myprt<<" InTraj "<<slc.slHits[iht].InTraj;
        myprt<<" Chg "<<(int)hit.Integral();
        myprt<<" Signal? "<<sigOK;
      }
      if(slc.slHits[iht].InTraj == 0 && delta < bigDelta && hitsInMultiplet.size() < 3 && !tj.AlgMod[kRvPrp]) {
        // An available hit that is just outside the window that is not part of a large multiplet
        bigDelta = delta;
        imBig = iht;
      }
      if(longPulseHit) {
        if(delta > 3) continue;
      } else {
        if(delta > deltaCut) continue;
      }

      if(std::find(closeHits.begin(), closeHits.end(), iht) != closeHits.end()) continue;
      closeHits.push_back(iht);
      if(hitsInMultiplet.size() > 1) {
        // include all the hits in a multiplet
        for(auto& jht : hitsInMultiplet) {
          if(slc.slHits[jht].InTraj == tj.ID) continue;
          if(std::find(closeHits.begin(), closeHits.end(), jht) != closeHits.end()) continue;
          closeHits.push_back(jht);
        } // jht
      } // multiplicity > 1
    } // iht

    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<"closeHits ";
      for(auto iht : closeHits) myprt<<" "<<PrintHit(slc.slHits[iht]);
      if(imBig < slc.slHits.size()) {
        myprt<<" imBig "<<PrintHit(slc.slHits[imBig]);
      } else {
        myprt<<" imBig "<<imBig;
      }
    }
    // check the srcHit collection if it is defined. Add the TP to the trajectory if
    // there is NO hit in the allHits collection but there is a hit in srcHit collection. We
    // can't use it for fitting, etc however
    bool nearSrcHit = false;
    if(!sigOK) nearSrcHit = NearbySrcHit(planeID, wire, (float)rawProjTick, (float)rawProjTick);
    sigOK = sigOK || !closeHits.empty() || nearSrcHit;

    if(!sigOK) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" no signal on any wire at tp.Pos "<<tp.Pos[0]<<" "<<tp.Pos[1]<<" tick "<<(int)tp.Pos[1]/tcc.unitsPerTick<<" closeHits size "<<closeHits.size();
      return;
    }
    if(imBig < slc.slHits.size() && closeHits.empty()) {
      closeHits.push_back(imBig);
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Added bigDelta hit "<<PrintHit(slc.slHits[imBig])<<" w delta = "<<bigDelta;
    }
    if(closeHits.size() > 16) closeHits.resize(16);
    if(nearSrcHit) tp.Environment[kEnvNearSrcHit] = true;
    tp.Hits.insert(tp.Hits.end(), closeHits.begin(), closeHits.end());

    // reset UseHit and assume that none of these hits will be used (yet)
    tp.UseHit.reset();
   // decide which of these hits should be used in the fit. Use a generous maximum delta
    // and require a charge check if we're not just starting out
    bool useChg = true;
    if(tj.Strategy[kStiffEl] || tj.Strategy[kSlowing]) useChg = false;
    FindUseHits(slc, tj, ipt, 10, useChg);
    DefineHitPos(slc, tp);
    SetEndPoints(tj);
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" number of close hits "<<closeHits.size()<<" used hits "<<NumHitsInTP(tp, kUsedHits);
  } // AddHits


  ////////////////////////////////////////////////
  void AddLAHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, bool& sigOK)
  {
    // Very Large Angle version of AddHits to be called for the last angle range

    if(ipt > tj.Pts.size() - 1) return;
    TrajPoint& tp = tj.Pts[ipt];
    tp.Hits.clear();
    tp.UseHit.reset();
    sigOK = false;

    unsigned short plane = DecodeCTP(tj.CTP).Plane;

    // look at adjacent wires for larger angle trajectories
    // We will check the most likely wire first
    std::vector<int> wires(1);
    wires[0] = std::nearbyint(tp.Pos[0]);
    if(wires[0] < 0 || wires[0] > (int)slc.lastWire[plane]-1) return;

    if(tp.AngleCode != 2) {
      mf::LogVerbatim("TC")<<"AddLAHits called with a bad angle code. "<<tp.AngleCode<<" Don't do this";
      return;
    }
    // and the adjacent wires next in the most likely order only
    // after the first two points have been defined
    if(ipt > 1) {
      if(tp.Dir[0] > 0) {
        if(wires[0] < (int)slc.lastWire[plane]-1) wires.push_back(wires[0] + 1);
        if(wires[0] > 0) wires.push_back(wires[0] - 1);
      } else {
        if(wires[0] > 0) wires.push_back(wires[0] - 1);
        if(wires[0] < (int)slc.lastWire[plane]-1) wires.push_back(wires[0] + 1);
      }
    } // ipt > 0 ...

    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<" AddLAHits: Pos "<<PrintPos(slc, tp)<<" tp.AngleCode "<<tp.AngleCode<<" Wires under consideration";
      for(auto& wire : wires) myprt<<" "<<wire;
    }

    // a temporary tp that we can move around
    TrajPoint ltp = tp;
    // do this while testing
    sigOK = false;

    tp.Hits.clear();
    std::array<int, 2> wireWindow;
    std::array<float, 2> timeWindow;
    float pos1Window = tcc.VLAStepSize/2;
    timeWindow[0] = ltp.Pos[1] - pos1Window;
    timeWindow[1] = ltp.Pos[1] + pos1Window;
    // Put the existing hits in to a vector so we can ensure that they aren't added again
    std::vector<unsigned int> oldHits = PutTrajHitsInVector(tj, kAllHits);

    for(unsigned short ii = 0; ii < wires.size(); ++ii) {
      int wire = wires[ii];
      if(wire < 0 || wire > (int)slc.lastWire[plane]) continue;
      // Assume a signal exists on a dead wire
      if(slc.wireHitRange[plane][wire].first == UINT_MAX) sigOK = true;
      if(slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
      wireWindow[0] = wire;
      wireWindow[1] = wire;
      bool hitsNear;
      // Look for hits using the requirement that the timeWindow overlaps with the hit StartTick and EndTick
      std::vector<unsigned int> closeHits = FindCloseHits(slc, wireWindow, timeWindow, plane, kAllHits, true, hitsNear);
      if(hitsNear) sigOK = true;
      for(auto& iht : closeHits) {
        // Ensure that none of these hits are already used by this trajectory
        if(slc.slHits[iht].InTraj == tj.ID) continue;
        // or in another trajectory in any previously added point
        if(std::find(oldHits.begin(), oldHits.end(), iht) != oldHits.end()) continue;
        tp.Hits.push_back(iht);
      }
    } // ii

    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<" LAPos "<<PrintPos(slc, ltp)<<" Tick window "<<(int)(timeWindow[0]/tcc.unitsPerTick)<<" to "<<(int)(timeWindow[1]/tcc.unitsPerTick);
      for(auto& iht : tp.Hits) myprt<<" "<<PrintHit(slc.slHits[iht]);
    } // prt

    // no hits found
    if(tp.Hits.empty()) return;

    if(tp.Hits.size() > 16) tp.Hits.resize(16);

    tp.UseHit.reset();
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      unsigned int iht = tp.Hits[ii];
      if(slc.slHits[iht].InTraj != 0) continue;
      tp.UseHit[ii] = true;
      slc.slHits[iht].InTraj = tj.ID;
    } // ii
    DefineHitPos(slc, tp);
    SetEndPoints(tj);
    UpdateTjChgProperties("ALAH", slc, tj, tcc.dbgStp);

  } // AddLAHits

  //////////////////////////////////////////
  void ReversePropagate(TCSlice& slc, Trajectory& tj)
  {
    // Reverse the trajectory and step in the opposite direction. The
    // updated trajectory is returned if this process is successful

    if(!tcc.useAlg[kRvPrp]) return;

    if(tj.Pts.size() < 6) return;
    // only do this once
    if(tj.AlgMod[kRvPrp]) return;

    // This code can't handle VLA trajectories
    if(tj.Pts[tj.EndPt[0]].AngleCode == 2) return;

    bool prt = (tcc.dbgStp || tcc.dbgAlg[kRvPrp]);

    // this function requires the first TP be included in the trajectory.
    if(tj.EndPt[0] > 0) {
      tj.Pts.erase(tj.Pts.begin(), tj.Pts.begin() + tj.EndPt[0]);
      SetEndPoints(tj);
    }

    if(prt) mf::LogVerbatim("TC")<<"ReversePropagate: Prepping Tj "<<tj.ID<<" incoming StepDir "<<tj.StepDir;

    short stepDir = tj.StepDir;

    // find the wire on which the first TP resides
    unsigned int wire0 = std::nearbyint(tj.Pts[0].Pos[0]);
    unsigned int nextWire = wire0 - tj.StepDir;

    // check for dead wires
    geo::PlaneID planeID = DecodeCTP(tj.CTP);
    unsigned short ipl = planeID.Plane;
    while(nextWire > slc.firstWire[ipl] && nextWire < slc.lastWire[ipl]) {
      if(evt.goodWire[ipl][nextWire]) break;
      nextWire -= tj.StepDir;
    }
    if(nextWire == slc.lastWire[ipl] - 1) return;
    // clone the first point
    TrajPoint tp = tj.Pts[0];
    // strip off the hits
    tp.Hits.clear(); tp.UseHit.reset();
    // move it to the next wire (in the opposite direction of the step direction)
    MoveTPToWire(tp, (float)nextWire);
    // find close unused hits near this position
    float maxDelta = 10 * tj.Pts[tj.EndPt[1]].DeltaRMS;
    if(!FindCloseHits(slc, tp, maxDelta, kUnusedHits)) return;
    if(prt) mf::LogVerbatim("TC")<<" nUnused hits "<<tp.Hits.size()<<" at Pos "<<PrintPos(slc, tp);
    if(tp.Hits.empty()) return;
    // There are hits on the next wire. Make a working copy of the trajectory, reverse it and try
    // to extend it with StepAway
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" tp.Hits ";
      for(auto& iht : tp.Hits) myprt<<" "<<PrintHit(slc.slHits[iht])<<"_"<<slc.slHits[iht].InTraj;
    } // tcc.dbgStp
    //
    // Make a working copy of tj
    Trajectory tjWork = tj;
    // So the first shall be last and the last shall be first
    ReverseTraj(slc, tjWork);
    // Flag it to use special cuts in StepAway
    tjWork.AlgMod[kRvPrp] = true;
    // save the strategy word and set it to normal
    auto saveStrategy = tjWork.Strategy;
    tjWork.Strategy.reset();
    tjWork.Strategy[kNormal] = true;
    // Reduce the number of fitted points to a small number
    unsigned short lastPt = tjWork.Pts.size() - 1;
    if(lastPt < 4) return;
    // update the charge
    float chg = 0;
    float cnt = 0;
    for(unsigned short ii = 0; ii < 4; ++ii) {
      unsigned short ipt = lastPt - ii;
      if(tjWork.Pts[ipt].Chg == 0) continue;
      chg += tjWork.Pts[ipt].Chg;
      ++cnt;
    } // ii
    if(cnt == 0) return;
    if(cnt > 1) tjWork.Pts[lastPt].AveChg = chg / cnt;
    StepAway(slc, tjWork);
    if(!tj.IsGood) {
      if(prt) mf::LogVerbatim("TC")<<" ReversePropagate StepAway failed";
      return;
    }
    tjWork.Strategy = saveStrategy;
    // check the new stopping point
    ChkStopEnd1(slc, tjWork, tcc.dbgStp);
    // restore the original direction
    if(tjWork.StepDir != stepDir) ReverseTraj(slc, tjWork);
    tj = tjWork;
    // TODO: Maybe UpdateTjChgProperties should be called here
    // re-check the ends
    ChkStop(slc, tj);

  } // ReversePropagate

  ////////////////////////////////////////////////
  void GetHitMultiplet(const TCSlice& slc, unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet, bool useLongPulseHits)
  {
    // This function attempts to return a list of hits in the current slice that are close to the
    // hit specified by theHit and that are similar to it. If theHit is a high-pulseheight hit (aka imTall)
    // and has an RMS similar to a hit on a small angle trajectory (aka Narrow) and is embedded in a series of
    // nearby low-pulseheight wide hits, the hit multiplet will consist of the single Tall and Narrow hit. On the
    // other hand, if theHit references a short and not-narrow hit, all of the hits in the series of nearby
    // hits will be returned. The localIndex is the index of theHit in hitsInMultiplet and shouldn't be
    // confused with the recob::Hit LocalIndex
    hitsInMultiplet.clear();
    // check for flagrant errors
    if(theHit >= slc.slHits.size()) return;
    if(slc.slHits[theHit].InTraj == INT_MAX) return;
    if(slc.slHits[theHit].allHitsIndex >= (*evt.allHits).size()) return;

    auto& hit = (*evt.allHits)[slc.slHits[theHit].allHitsIndex];
    // handle long-pulse hits
    if(useLongPulseHits && LongPulseHit(hit)) {
      // return everything in the multiplet as defined by the hit finder, but check for errors
      short int hitMult = hit.Multiplicity();
      unsigned int lIndex = hit.LocalIndex();
      unsigned int firstHit = 0;
      if(lIndex < theHit) firstHit = theHit - lIndex;
      for(unsigned int ii = firstHit; ii < firstHit + hitMult; ++ii) {
        if(ii >= slc.slHits.size()) break;
        auto& tmp = (*evt.allHits)[slc.slHits[ii].allHitsIndex];
        if(tmp.Multiplicity() == hitMult) hitsInMultiplet.push_back(ii);
      } // ii
      return;
    } // LongPulseHit

    hitsInMultiplet.resize(1);
    hitsInMultiplet[0] = theHit;
    unsigned int theWire = hit.WireID().Wire;
    unsigned short ipl = hit.WireID().Plane;

    float theTime = hit.PeakTime();
    float theRMS = hit.RMS();
    float narrowHitCut = 1.5 * evt.aveHitRMS[ipl];
    bool theHitIsNarrow = (theRMS < narrowHitCut);
    float maxPeak = hit.PeakAmplitude();
    unsigned int imTall = theHit;
    unsigned short nNarrow = 0;
    if(theHitIsNarrow) nNarrow = 1;
    // look for hits < theTime but within hitSep
    if(theHit > 0) {
      for(unsigned int iht = theHit - 1; iht != 0; --iht) {
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        if(hit.WireID().Wire != theWire) break;
        if(hit.WireID().Plane != ipl) break;
        float hitSep = tcc.multHitSep * theRMS;
        float rms = hit.RMS();
        if(rms > theRMS) {
          hitSep = tcc.multHitSep * rms;
          theRMS = rms;
        }
        float dTick = std::abs(hit.PeakTime() - theTime);
        if(dTick > hitSep) break;
        hitsInMultiplet.push_back(iht);
        if(rms < narrowHitCut) ++nNarrow;
        float peakAmp = hit.PeakAmplitude();
        if(peakAmp > maxPeak) {
          maxPeak = peakAmp;
          imTall = iht;
        }
        theTime = hit.PeakTime();
        if(iht == 0) break;
      } // iht
    } // iht > 0
    // reverse the order so that hitsInMuliplet will be
    // returned in increasing time order
    if(hitsInMultiplet.size() > 1) std::reverse(hitsInMultiplet.begin(), hitsInMultiplet.end());
    // look for hits > theTime but within hitSep
    theTime = hit.PeakTime();
    theRMS = hit.RMS();
    for(unsigned int iht = theHit + 1; iht < slc.slHits.size(); ++iht) {
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      if(hit.WireID().Wire != theWire) break;
      if(hit.WireID().Plane != ipl) break;
      if(slc.slHits[iht].InTraj == INT_MAX) continue;
      float hitSep = tcc.multHitSep * theRMS;
      float rms = hit.RMS();
      if(rms > theRMS) {
        hitSep = tcc.multHitSep * rms;
        theRMS = rms;
      }
      float dTick = std::abs(hit.PeakTime() - theTime);
      if(dTick > hitSep) break;
      hitsInMultiplet.push_back(iht);
      if(rms < narrowHitCut) ++nNarrow;
      float peakAmp = hit.PeakAmplitude();
      if(peakAmp > maxPeak) {
        maxPeak = peakAmp;
        imTall = iht;
      }
      theTime = hit.PeakTime();
    } // iht
    if(hitsInMultiplet.size() == 1) return;

    // Don't make a multiplet that includes a tall narrow hit with short fat hits
    if(nNarrow == hitsInMultiplet.size()) return;
    if(nNarrow == 0) return;

    if(theHitIsNarrow && theHit == imTall) {
      // theHit is narrow and it is the highest amplitude hit in the multiplet. Ignore any
      // others that are short and fat
      auto tmp = hitsInMultiplet;
      tmp.resize(1);
      tmp[0] = theHit;
      hitsInMultiplet = tmp;
    } else {
      // theHit is not narrow and it is not the tallest. Ignore a single hit if it is
      // the tallest and narrow
      auto& hit = (*evt.allHits)[slc.slHits[imTall].allHitsIndex];
      if(hit.RMS() < narrowHitCut) {
        unsigned short killMe = 0;
        for(unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
          if(hitsInMultiplet[ii] == imTall) {
            killMe = ii;
            break;
          }
        } // ii
        hitsInMultiplet.erase(hitsInMultiplet.begin() + killMe);
      } // slc.slHits[imTall].RMS < narrowHitCut
    } // narrow / tall test

  } // GetHitMultiplet

  //////////////////////////////////////////
  float HitTimeErr(const TCSlice& slc, unsigned int iht)
  {
    if(iht > slc.slHits.size() - 1) return 0;
    auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
    return hit.RMS() * tcc.unitsPerTick * tcc.hitErrFac * hit.Multiplicity();
  } // HitTimeErr

  //////////////////////////////////////////
  float HitsTimeErr2(const TCSlice& slc, const std::vector<unsigned int>& hitVec)
  {
    // Estimates the error^2 of the time using all hits in hitVec
    if(hitVec.empty()) return 0;
    float err = tcc.hitErrFac * HitsRMSTime(slc, hitVec, kUnusedHits);
    return err * err;
  } // HitsTimeErr2

  //////////////////////////////////////////
  void DefineHitPos(TCSlice& slc, TrajPoint& tp)
  {
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point

    tp.Chg = 0;
    if(tp.Hits.empty()) return;

    unsigned short nused = 0;
    unsigned int iht = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      ++nused;
      iht = tp.Hits[ii];
      if(iht >= slc.slHits.size()) return;
      if(slc.slHits[iht].allHitsIndex >= (*evt.allHits).size()) return;
    }
    if(nused == 0) return;

    // don't bother with rest of this if there is only one hit since it can
    // only reside on one wire
    if(nused == 1) {
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      tp.Chg = hit.Integral();
      tp.HitPos[0] = hit.WireID().Wire;
      tp.HitPos[1] = hit.PeakTime() * tcc.unitsPerTick;
      if(LongPulseHit(hit)) {
        // give it a huge error^2 since the position is not well defined
        tp.HitPosErr2 = 100;
      } else {
        // Normalize to 1 WSE path length
        float pathInv = std::abs(tp.Dir[0]);
        if(pathInv < 0.05) pathInv = 0.05;
        tp.Chg *= pathInv;
        float wireErr = tp.Dir[1] * 0.289;
        float timeErr = tp.Dir[0] * HitTimeErr(slc, iht);
        tp.HitPosErr2 = wireErr * wireErr + timeErr * timeErr;
      }
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"DefineHitPos: singlet "<<std::fixed<<std::setprecision(1)<<tp.HitPos[0]<<":"<<(int)(tp.HitPos[1]/tcc.unitsPerTick)<<" ticks. HitPosErr "<<sqrt(tp.HitPosErr2);
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
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      chg = hit.Integral();
      unsigned int wire = hit.WireID().Wire;
      if(wire < loWire) loWire = wire;
      if(wire > hiWire) hiWire = wire;
      newpos[0] += chg * wire;
      newpos[1] += chg * hit.PeakTime();
      tp.Chg += chg;
      hitVec.push_back(iht);
    } // ii

    if(tp.Chg == 0) return;

    tp.HitPos[0] = newpos[0] / tp.Chg;
    tp.HitPos[1] = newpos[1] * tcc.unitsPerTick / tp.Chg;
    // Normalize to 1 WSE path length
    float pathInv = std::abs(tp.Dir[0]);
    if(pathInv < 0.05) pathInv = 0.05;
    tp.Chg *= pathInv;
    // Error is the wire error (1/sqrt(12))^2 if all hits are on one wire.
    // Scale it by the wire range
    float dWire = 1 + hiWire - loWire;
    float wireErr = tp.Dir[1] * dWire * 0.289;
    float timeErr2 = tp.Dir[0] * tp.Dir[0] * HitsTimeErr2(slc, hitVec);
    tp.HitPosErr2 = wireErr * wireErr + timeErr2;
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"DefineHitPos: multiplet "<<std::fixed<<std::setprecision(1)<<tp.HitPos[0]<<":"<<(int)(tp.HitPos[1]/tcc.unitsPerTick)<<" ticks. HitPosErr "<<sqrt(tp.HitPosErr2);

  } // DefineHitPos


  //////////////////////////////////////////
  void FindUseHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg)
  {
    // Hits have been associated with trajectory point ipt but none are used. Here we will
    // decide which hits to use.

    if(ipt > tj.Pts.size() - 1) return;
    TrajPoint& tp = tj.Pts[ipt];

    if(tp.Hits.empty()) return;
    // don't check charge when starting out
    if(ipt < 5) useChg = false;
    float chgPullCut = 1000;
    if(useChg) chgPullCut = tcc.chargeCuts[0];
    // open it up for RevProp, since we might be following a stopping track
    if(tj.AlgMod[kRvPrp]) chgPullCut *= 2;
    if(tj.MCSMom < 30) chgPullCut *= 2;

    bool ignoreLongPulseHits = false;
    unsigned short npts = tj.EndPt[1] - tj.EndPt[0] + 1;
    if(npts < 10 || tj.AlgMod[kRvPrp]) ignoreLongPulseHits = true;
    float expectedHitsRMS = ExpectedHitsRMS(slc, tp);
    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"FUH:  maxDelta "<<maxDelta<<" useChg requested "<<useChg<<" Norm AveChg "<<(int)tp.AveChg<<" tj.ChgRMS "<<std::setprecision(2)<<tj.ChgRMS<<" chgPullCut "<<chgPullCut<<" TPHitsRMS "<<(int)TPHitsRMSTick(slc, tp, kUnusedHits)<<" ExpectedHitsRMS "<<(int)expectedHitsRMS<<" AngCode "<<tp.AngleCode;
    }

    // inverse of the path length for normalizing hit charge to 1 WSE unit
    float pathInv = std::abs(tp.Dir[0]);
    if(pathInv < 0.05) pathInv = 0.05;

    // Find the hit that has the smallest delta and the number of available hits
    tp.Delta = maxDelta;
    float delta;
    unsigned short imbest = USHRT_MAX;
    std::vector<float> deltas(tp.Hits.size(), 100);
    // keep track of the best delta - even if it is used
    float bestDelta = maxDelta;
    unsigned short nAvailable = 0;
    unsigned short firstAvailable = USHRT_MAX;
    unsigned short lastAvailable = USHRT_MAX;
    unsigned short firstUsed = USHRT_MAX;
    unsigned short imBadRecoHit = USHRT_MAX;
    bool bestHitLongPulse = false;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      tp.UseHit[ii] = false;
      unsigned int iht = tp.Hits[ii];
      if(iht >= slc.slHits.size()) continue;
      if(slc.slHits[iht].allHitsIndex >= (*evt.allHits).size()) continue;
      delta = PointTrajDOCA(slc, iht, tp);
      if(delta < bestDelta) bestDelta = delta;
      if(slc.slHits[iht].InTraj > 0) {
        if(firstUsed == USHRT_MAX) firstUsed = ii;
        continue;
      }
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      if(ignoreLongPulseHits && LongPulseHit(hit)) continue;
      if(hit.GoodnessOfFit() < 0 || hit.GoodnessOfFit() > 100) imBadRecoHit = ii;
      if(firstAvailable == USHRT_MAX) firstAvailable = ii;
      lastAvailable = ii;
      ++nAvailable;
      if(tcc.dbgStp) {
        if(useChg) {
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<" "<<ii<<"  "<<PrintHit(slc.slHits[iht])<<" delta "<<delta<<" Norm Chg "<<(int)(hit.Integral() * pathInv);
        } else {
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<" "<<ii<<"  "<<PrintHit(slc.slHits[iht])<<" delta "<<delta;
        }
      } // tcc.dbgStp
      deltas[ii] = delta;
      if(delta < tp.Delta) {
        tp.Delta = delta;
        imbest = ii;
        bestHitLongPulse = LongPulseHit(hit);
      }
    } // ii

    float chgWght = 0.5;

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" firstAvailable "<<firstAvailable<<" lastAvailable "<<lastAvailable<<" firstUsed "<<firstUsed<<" imbest "<<imbest<<" single hit. tp.Delta "<<std::setprecision(2)<<tp.Delta<<" bestDelta "<<bestDelta<<" path length "<<1 / pathInv<<" imBadRecoHit "<<imBadRecoHit;
    if(imbest == USHRT_MAX || nAvailable == 0) return;
    unsigned int bestDeltaHit = tp.Hits[imbest];

    // just use the best hit if we are tracking a high energy electron and the best hit is a long pulse hit
    if(tj.Strategy[kStiffEl] && bestHitLongPulse) {
      tp.UseHit[imbest] = true;
      slc.slHits[bestDeltaHit].InTraj = tj.ID;
      return;
    }

    // Don't try to use a multiplet if a hit in the middle is in a different trajectory
    if(tp.Hits.size() > 2 && nAvailable > 1 && firstUsed != USHRT_MAX && firstAvailable < firstUsed && lastAvailable > firstUsed) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" A hit in the middle of the multiplet is used. Use only the best hit";
      tp.UseHit[imbest] = true;
      slc.slHits[bestDeltaHit].InTraj = tj.ID;
      return;
    } // Used hit inside multiplet

    if(tp.AngleCode == 1) {
      // Get the hits that are in the same multiplet as bestDeltaHit
      std::vector<unsigned int> hitsInMultiplet;
      GetHitMultiplet(slc, bestDeltaHit, hitsInMultiplet, false);
      if(tcc.dbgStp) {
        mf::LogVerbatim myprt("TC");
        myprt<<" bestDeltaHit "<<PrintHit(slc.slHits[bestDeltaHit]);
        myprt<<" in multiplet:";
        for(auto& iht : hitsInMultiplet) myprt<<" "<<PrintHit(slc.slHits[iht]);
      }
      // Consider the case where a bad reco hit might be better. It is probably wider and
      // has more charge
      if(imBadRecoHit != USHRT_MAX) {
        unsigned int iht = tp.Hits[imBadRecoHit];
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        if(hit.RMS() > HitsRMSTick(slc, hitsInMultiplet, kUnusedHits)) {
          if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Using imBadRecoHit "<<PrintHit(slc.slHits[iht]);
          tp.UseHit[imBadRecoHit] = true;
          slc.slHits[iht].InTraj = tj.ID;
          return;
        }
      } // bad reco hit
      // Use the hits in the multiplet instead
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(slc.slHits[iht].InTraj > 0) continue;
        if(std::find(hitsInMultiplet.begin(), hitsInMultiplet.end(), iht) == hitsInMultiplet.end()) continue;
        tp.UseHit[ii] = true;
        slc.slHits[iht].InTraj = tj.ID;
      } // ii
      return;
    } // isLA

    // don't use the best UNUSED hit if the best delta is for a USED hit and it is much better
    // TY: ignore for RevProp
    if(bestDelta < 0.5 * tp.Delta && !tj.AlgMod[kRvPrp]) return;

    if(!useChg || (useChg && (tj.AveChg <= 0 || tj.ChgRMS <= 0))) {
      // necessary quantities aren't available for more careful checking
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" tj.AveChg "<<tj.AveChg<<" or tj.ChgRMS "<<tj.ChgRMS<<". Use the best hit";
      tp.UseHit[imbest] = true;
      slc.slHits[bestDeltaHit].InTraj = tj.ID;
      return;
    }

    // Don't try to get fancy if we are tracking a stiff tj
    if(tj.PDGCode == 13 && bestDelta < 0.5) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Tracking muon. Use the best hit";
      tp.UseHit[imbest] = true;
      slc.slHits[bestDeltaHit].InTraj = tj.ID;
      return;
    }

    // The best hit is the only one available or this is a small angle trajectory
    if(nAvailable == 1 || tp.AngleCode == 0) {
      auto& hit = (*evt.allHits)[slc.slHits[bestDeltaHit].allHitsIndex];
      float aveChg = tp.AveChg;
      if(aveChg <= 0) aveChg = tj.AveChg;
      if(aveChg <= 0) aveChg = hit.Integral();
      float chgRMS = tj.ChgRMS;
      if(chgRMS < 0.2) chgRMS = 0.2;
      float bestDeltaHitChgPull = std::abs(hit.Integral() * pathInv / aveChg - 1) / chgRMS;
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" bestDeltaHitChgPull "<<bestDeltaHitChgPull<<" chgPullCut "<<chgPullCut;
      if(bestDeltaHitChgPull < chgPullCut || tp.Delta < 0.1) {
        tp.UseHit[imbest] = true;
        slc.slHits[bestDeltaHit].InTraj = tj.ID;
      } // good charge or very good delta
      return;
    } // bestDeltaHitMultiplicity == 1

    // Find the expected width for the angle of this TP (ticks)
    float expectedWidth = ExpectedHitsRMS(slc, tp);

    // Handle two available hits
    if(nAvailable == 2) {
      // See if these two are in the same multiplet and both are available
      std::vector<unsigned int> tHits;
      GetHitMultiplet(slc, bestDeltaHit, tHits, false);
      // ombest is the index of the other hit in tp.Hits that is in the same multiplet as bestDeltaHit
      // if we find it
      unsigned short ombest = USHRT_MAX;
      unsigned int otherHit = INT_MAX;
      if(tHits.size() == 2) {
        unsigned short localIndex = 0;
        if(tHits[0] == bestDeltaHit) localIndex = 1;
        otherHit = tHits[1 - localIndex];
        // get the index of this hit in tp.Hits
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(slc.slHits[tp.Hits[ii]].InTraj > 0) continue;
          if(tp.Hits[ii] == otherHit) {
            ombest = ii;
            break;
          }
        } // ii
      } // tHits.size() == 2
      if(tcc.dbgStp) {
        mf::LogVerbatim("TC")<<" Doublet: imbest "<<imbest<<" ombest "<<ombest;
      }
      // The other hit exists in the tp and it is available
      if(ombest < tp.Hits.size()) {
        // compare the best delta hit and the other hit separately and the doublet as a merged pair
        float bestHitDeltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(slc, bestDeltaHit);
        // Construct a FOM starting with the delta pull
        float bestDeltaHitFOM = deltas[imbest] /  bestHitDeltaErr;
        if(bestDeltaHitFOM < 0.5) bestDeltaHitFOM = 0.5;
        // multiply by the charge pull if it is significant
        auto& hit = (*evt.allHits)[slc.slHits[bestDeltaHit].allHitsIndex];
        float bestDeltaHitChgPull = std::abs(hit.Integral() * pathInv / tp.AveChg - 1) / tj.ChgRMS;
        if(bestDeltaHitChgPull > 1) bestDeltaHitFOM *= chgWght * bestDeltaHitChgPull;
        // scale by the ratio
        float rmsRat = hit.RMS() / expectedWidth;
        if(rmsRat < 1) rmsRat = 1 / rmsRat;
        bestDeltaHitFOM *= rmsRat;
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<" bestDeltaHit FOM "<<deltas[imbest]/bestHitDeltaErr<<" bestDeltaHitChgPull "<<bestDeltaHitChgPull<<" rmsRat "<<rmsRat<<" bestDeltaHitFOM "<<bestDeltaHitFOM;
        // Now do the same for the other hit
        float otherHitDeltaErr = std::abs(tp.Dir[1]) * 0.17 + std::abs(tp.Dir[0]) * HitTimeErr(slc, otherHit);
        float otherHitFOM = deltas[ombest] /  otherHitDeltaErr;
        if(otherHitFOM < 0.5) otherHitFOM = 0.5;
        auto& ohit = (*evt.allHits)[slc.slHits[otherHit].allHitsIndex];
        float otherHitChgPull = std::abs(ohit.Integral() * pathInv / tp.AveChg - 1) / tj.ChgRMS;
        if(otherHitChgPull > 1) otherHitFOM *= chgWght * otherHitChgPull;
        rmsRat = ohit.RMS() / expectedWidth;
        if(rmsRat < 1) rmsRat = 1 / rmsRat;
        otherHitFOM *= rmsRat;
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<" otherHit FOM "<<deltas[ombest]/otherHitDeltaErr<<" otherHitChgPull "<<otherHitChgPull<<" rmsRat "<<rmsRat<<" otherHitFOM "<<otherHitFOM;
        // And for the doublet
        float doubletChg = hit.Integral() + ohit.Integral();
        float doubletTime = (hit.Integral() * hit.PeakTime() + ohit.Integral() * ohit.PeakTime()) / doubletChg;
        doubletChg *= pathInv;
        doubletTime *= tcc.unitsPerTick;
        float doubletWidthTick = TPHitsRMSTick(slc, tp, kUnusedHits);
        float doubletRMSTimeErr = doubletWidthTick * tcc.unitsPerTick;
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<" doublet Chg "<<doubletChg<<" doubletTime "<<doubletTime<<" doubletRMSTimeErr "<<doubletRMSTimeErr;
        float doubletFOM = PointTrajDOCA(slc, tp.Pos[0], doubletTime, tp) / doubletRMSTimeErr;
        if(doubletFOM < 0.5) doubletFOM = 0.5;
        float doubletChgPull = std::abs(doubletChg * pathInv / tp.AveChg - 1) / tj.ChgRMS;
        if(doubletChgPull > 1) doubletFOM *= chgWght * doubletChgPull;
        rmsRat = doubletWidthTick / expectedWidth;
        if(rmsRat < 1) rmsRat = 1 / rmsRat;
        doubletFOM *= rmsRat;
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<" doublet FOM "<<PointTrajDOCA(slc, tp.Pos[0], doubletTime, tp)/doubletRMSTimeErr<<" doubletChgPull "<<doubletChgPull<<" rmsRat "<<rmsRat<<" doubletFOM "<<doubletFOM;
        if(doubletFOM < bestDeltaHitFOM && doubletFOM < otherHitFOM) {
          tp.UseHit[imbest] = true;
          slc.slHits[bestDeltaHit].InTraj = tj.ID;
          tp.UseHit[ombest] = true;
          slc.slHits[otherHit].InTraj = tj.ID;
        } else {
          // the doublet is not the best
          if(bestDeltaHitFOM < otherHitFOM) {
            tp.UseHit[imbest] = true;
            slc.slHits[bestDeltaHit].InTraj = tj.ID;
          } else {
            tp.UseHit[ombest] = true;
            slc.slHits[otherHit].InTraj = tj.ID;
          } // otherHit is the best
        } // doublet is not the best
      } else {
        // the other hit isn't available. Just use the singlet
        tp.UseHit[imbest] = true;
        slc.slHits[bestDeltaHit].InTraj = tj.ID;
      }
      return;
    } // nAvailable == 2
    float hitsWidth = TPHitsRMSTick(slc, tp, kUnusedHits);
    float maxTick = tp.Pos[1] / tcc.unitsPerTick + 0.6 * expectedWidth;
    float minTick = tp.Pos[1] / tcc.unitsPerTick - 0.6 * expectedWidth;
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<" Multiplet: hitsWidth "<<hitsWidth<<" expectedWidth "<<expectedWidth<<" tick range "<<(int)minTick<<" "<<(int)maxTick;
    // use all of the hits in the tick window
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      unsigned int iht = tp.Hits[ii];
      if(slc.slHits[iht].InTraj > 0) continue;
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      if(hit.PeakTime() < minTick) continue;
      if(hit.PeakTime() > maxTick) continue;
      tp.UseHit[ii] = true;
      slc.slHits[iht].InTraj = tj.ID;
    }

  } // FindUseHits

  //////////////////////////////////////////
  void UpdateDeltaRMS(TCSlice& slc, Trajectory& tj)
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
      sum += PointTrajDOCA(slc, tj.Pts[ipt].Pos[0], tj.Pts[ipt].Pos[1], lastTP);
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
  void MaskBadTPs(TCSlice& slc, Trajectory& tj, float const& maxChi)
  {
    // Remove TPs that have the worst values of delta until the fit chisq < maxChi

    if(!tcc.useAlg[kMaskBadTPs]) return;
    //don't use this function for reverse propagation
    if(!tcc.useAlg[kRvPrp]) return;

    bool prt = (tcc.dbgStp || tcc.dbgAlg[kMaskBadTPs]);

    if(tj.Pts.size() < 3) {
      //      mf::LogError("TC")<<"MaskBadTPs: Trajectory ID "<<tj.ID<<" too short to mask hits ";
      tj.IsGood = false;
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
      UnsetUsedHits(slc, tj.Pts[imBad]);
      FitTraj(slc, tj);
      if(prt) mf::LogVerbatim("TC")<<"  after FitTraj "<<lastTP.FitChi;
      tj.AlgMod[kMaskBadTPs] = true;
      ++nit;
    } // lastTP.FItChi > maxChi && nit < 3

  } // MaskBadTPs

  ////////////////////////////////////////////////
  bool MaskedHitsOK(TCSlice& slc, Trajectory& tj)
  {
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration (true) or whether to stop and possibly try with the next pass settings (false)

    if(!tcc.useAlg[kMaskHits]) return true;

    unsigned short lastPt = tj.Pts.size() - 1;
    if(tj.Pts[lastPt].Chg > 0) return true;
    unsigned short endPt = tj.EndPt[1];

    // count the number of points w/o used hits and the number with one unused hit
    unsigned short nMasked = 0;
    unsigned short nOneHit = 0;
    unsigned short nOKChg = 0;
    unsigned short nOKDelta = 0;
    // number of points with Pos > HitPos
    unsigned short nPosDelta = 0;
    // number of points with Delta increasing vs ipt
    unsigned short nDeltaIncreasing = 0;
    // Fake this a bit to simplify comparing the counts
    float prevDelta = tj.Pts[endPt].Delta;
    float maxOKDelta = 10 * tj.Pts[endPt].DeltaRMS;
    float maxOKChg = 0;
    // find the maximum charge point on the trajectory
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) if(tj.Pts[ipt].Chg > maxOKChg) maxOKChg = tj.Pts[ipt].Chg;
    for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - ii;
      auto& tp = tj.Pts[ipt];
      if(tp.Chg > 0) break;
      unsigned short nUnusedHits = 0;
      float chg = 0;
      for(unsigned short jj = 0; jj < tp.Hits.size(); ++jj) {
        unsigned int iht = tp.Hits[jj];
        if(slc.slHits[iht].InTraj != 0) continue;
        ++nUnusedHits;
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        chg += hit.Integral();
      } // jj
      if(chg < maxOKChg) ++nOKChg;
      if(nUnusedHits == 1) ++nOneHit;
      if(tp.Delta < maxOKDelta) ++nOKDelta;
      // count the number of points with Pos time > HitPos time
      if(tp.Pos[1] > tp.HitPos[1]) ++nPosDelta;
      // The number of increasing delta points: Note implied absolute value
      if(tp.Delta < prevDelta) ++nDeltaIncreasing;
      prevDelta = tp.Delta;
      ++nMasked;
    } // ii

    // determine if the hits are wandering away from the trajectory direction. This will result in
    // nPosDelta either being ~0 or ~equal to the number of masked points. nPosDelta should have something
    // in between these two extremes if we are stepping through a messy region
    bool driftingAway = nMasked > 2 && (nPosDelta == 0 || nPosDelta == nMasked);
    // Note that nDeltaIncreasing is always positive
    if(driftingAway && nDeltaIncreasing < nMasked - 1) driftingAway = false;

    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"MHOK:  nMasked "<<nMasked<<" nOneHit "<<nOneHit<<" nOKChg "<<nOKChg<<" nOKDelta "<<nOKDelta<<" nPosDelta "<<nPosDelta<<" nDeltaIncreasing "<<nDeltaIncreasing<<" driftingAway? "<<driftingAway;
    }

    if(!driftingAway) {
      if(nMasked < 8 || nOneHit < 8) return true;
      if(nOKDelta != nMasked) return true;
      if(nOKChg != nMasked) return true;
    }

    // we would like to reduce the number of fitted points to a minimum and include
    // the masked hits, but we can only do that if there are enough points
    if(tj.Pts[endPt].NTPsFit <= tcc.minPtsFit[tj.Pass]) {
      // stop stepping if we have masked off more points than are in the fit
      if(nMasked > tj.Pts[endPt].NTPsFit) return false;
      return true;
    }
    // Reduce the number of points fit and try to include the points
    unsigned short newNTPSFit;
    if(tj.Pts[endPt].NTPsFit > 2 * tcc.minPtsFit[tj.Pass]) {
      newNTPSFit = tj.Pts[endPt].NTPsFit / 2;
    } else {
      newNTPSFit = tcc.minPtsFit[tj.Pass];
    }
    for(unsigned ipt = endPt + 1; ipt < tj.Pts.size(); ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(slc.slHits[iht].InTraj == 0) {
          tp.UseHit[ii] = true;
          slc.slHits[iht].InTraj = tj.ID;
          break;
        }
      } // ii
      DefineHitPos(slc, tp);
      SetEndPoints(tj);
      tp.NTPsFit = newNTPSFit;
      FitTraj(slc, tj);
      if(tcc.dbgStp) PrintTrajectory("MHOK", slc, tj, ipt);
    } // ipt

    tj.AlgMod[kMaskHits] = true;
    UpdateTjChgProperties("MHOK", slc, tj, tcc.dbgStp);
    return true;

  } // MaskedHitsOK

  ////////////////////////////////////////////////
  bool StopIfBadFits(TCSlice& slc, Trajectory& tj)
  {
    // Returns true if there are a number of Tps that were not used in the trajectory because the fit was poor and the
    // charge pull is not really high. This

    // don't consider muons
    if(tj.PDGCode == 13) return false;
    // or long straight Tjs
    if(tj.Pts.size() > 40 && tj.MCSMom > 200) return false;

    unsigned short nBadFit = 0;
    unsigned short nHiChg = 0;
    unsigned short cnt = 0;
    for(unsigned short ipt = tj.Pts.size() - 1; ipt > tj.EndPt[1]; --ipt ) {
      if(tj.Pts[ipt].FitChi > 2) ++nBadFit;
      if(tj.Pts[ipt].ChgPull > 3) ++nHiChg;
      ++cnt;
      if(cnt == 5) break;
    } // ipt

    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"StopIfBadFits: nBadFit "<<nBadFit<<" nHiChg "<<nHiChg;
    if(nBadFit > 3 && nHiChg == 0) return true;
    return false;

  } // StopIfBadFits

////////////////////////////////////////////////
  bool GottaKink(TCSlice& slc, Trajectory& tj, bool doTrim)
  {
    // This function returns true if it detects a kink in the trajectory
    // This function trims the points after a kink if one is found if doTrim is true.

    // tcc.kinkCuts[] fcl configuration
    // 0 = Number of TPs to fit at the end
    // 1 = Min kink significance
    // 2 = Use charge in significance calculation if > 0
    // 3 = 3D kink fit length (cm) - used in PFPUtils/SplitAtKinks

    // don't look for kinks if this looks a high energy electron
    // BB Return true if a kink was found but don't set the
    // stop-at-kink end flag
    if(tj.Strategy[kStiffEl]) return false;
    // Need at least 2 * kinkCuts[2] points with charge to find a kink
    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    unsigned short nPtsFit = tcc.kinkCuts[0];
    // Set nPtsFit for slowing tjs to the last TP NTPsFit
    if(tj.Strategy[kSlowing]) nPtsFit = tj.Pts[tj.EndPt[1]].NTPsFit;
    if(npwc < 2 * nPtsFit) return false;

    bool useCharge = (tcc.kinkCuts[2] > 0);

    // find the point where a kink is expected and fit the points after that point
    unsigned short fitPt = USHRT_MAX;
    unsigned short cnt = 0;
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.EndPt[1] - ii - 1;
      // stay away from the starting points which may be skewed if this is a
      // stopping track
      if(ipt <= tj.EndPt[0] + 2) break;
      if(tj.Pts[ipt].Chg <= 0) continue;
      ++cnt;
      // Note that the fitPt is not included in the fits in the kink significance so we need
      // one more point
      if(cnt > nPtsFit) {
        fitPt = ipt;
        break;
      }
    } // ii
    if(fitPt == USHRT_MAX) {
      if(tcc.dbgStp) {
        mf::LogVerbatim myprt("TC");
        myprt<<"GKv2 fitPt not valid. Counted "<<cnt<<" points. Need "<<nPtsFit;
      } // tcc.dbgStp
      return false;
    }

    tj.Pts[fitPt].KinkSig = KinkSignificance(slc, tj, fitPt, nPtsFit, useCharge, tcc.dbgStp);

    bool thisPtHasKink = (tj.Pts[fitPt].KinkSig > tcc.kinkCuts[1]);
    bool prevPtHasKink = (tj.Pts[fitPt - 1].KinkSig > tcc.kinkCuts[1]);
    if(tcc.dbgStp) {
      mf::LogVerbatim myprt("TC");
      myprt<<"GKv2 fitPt "<<fitPt<<" "<<PrintPos(slc, tj.Pts[fitPt]);
      myprt<<std::fixed<<std::setprecision(5);
      myprt<<" KinkSig "<<std::setprecision(5)<<tj.Pts[fitPt].KinkSig;
      myprt<<" prevPt significance "<<tj.Pts[fitPt - 1].KinkSig;
      if(!thisPtHasKink && !prevPtHasKink) myprt<<" no kink";
      if(thisPtHasKink && !prevPtHasKink) myprt<<" -> Start kink region";
      if(thisPtHasKink && prevPtHasKink) myprt<<" -> Inside kink region";
      if(!thisPtHasKink && prevPtHasKink) myprt<<" -> End kink region";
    } // dbgStp
    // See if we just passed a series of points having a high kink significance. If so,
    // then find the point with the maximum value and call that the kink point
    // Don't declare a kink (yet)
    if(thisPtHasKink) return false;
    // neither points are kink-like
    if(!prevPtHasKink) return false;

    // We have left a kink region. Find the point with the max likelihood and call
    // that the kink point
    float maxSig = tcc.kinkCuts[1];
    unsigned short kinkRegionLength = 0;
    unsigned short maxKinkPt = USHRT_MAX;
    for(unsigned short ipt = fitPt - 1; ipt > tj.EndPt[0]; --ipt) {
      auto& tp = tj.Pts[ipt];
      if(tp.KinkSig < 0) continue;
      if(tp.KinkSig > maxSig) {
        // track the max significance
        maxSig = tp.KinkSig;
        maxKinkPt = ipt;
      } // tp.KinkSig > maxSig
      // find the start of the kink region
      if(tp.KinkSig < tcc.kinkCuts[1]) break;
      ++kinkRegionLength;
    } // ipt
    if(maxKinkPt == USHRT_MAX) return false;
    // Require that the candidate kink be above the cut threshold for more than one point.
    // Scale the requirement by the number of points in the fit
    unsigned short kinkRegionLengthMin = 1 + nPtsFit / 5;
    if(tj.Strategy[kStiffMu]) kinkRegionLengthMin = 1 + nPtsFit / 3;
    if(kinkRegionLength < kinkRegionLengthMin) {
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<"GKv2: kink region too short "<<kinkRegionLength<<" Min "<<kinkRegionLengthMin;
      return false;
    }
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"GKv2:   kink at "<<PrintPos(slc, tj.Pts[maxKinkPt])<<std::setprecision(3)<<" maxSig "<<maxSig<<" kinkRegionLength "<<kinkRegionLength<<" Min "<<kinkRegionLengthMin;
    // don't alter the tj unless doTrim is true
    if(!doTrim) return true;
    // trim the points
    for(unsigned short ipt = maxKinkPt + 1; ipt <= tj.EndPt[1]; ++ipt) UnsetUsedHits(slc, tj.Pts[ipt]);
    SetEndPoints(tj);
    // trim another point if the charge of the last two points is wildly dissimilar
    float lastChg = tj.Pts[tj.EndPt[1]].Chg;
    float prevChg = tj.Pts[tj.EndPt[1] - 1].Chg;
    float chgAsym = std::abs(lastChg - prevChg) / (lastChg + prevChg);
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"GKv2: last point after trim "<<PrintPos(slc, tj.Pts[tj.EndPt[1]])<<" chgAsym "<<chgAsym;
    if(chgAsym > 0.1) {
      UnsetUsedHits(slc, tj.Pts[tj.EndPt[1]]);
      SetEndPoints(tj);
    }
    tj.EndFlag[1][kEndKink] = true;
    return true;

  } // GottaKink

  ////////////////////////////////////////////////
  TrajPoint CreateTPFromTj(TCSlice& slc, const Trajectory& tj)
  {
    // Create a trajectory point by averaging the position and direction of all
    // TPs in the trajectory. This is used in LastEndMerge
    TrajPoint tjtp;
    // set the charge invalid
    tjtp.Chg = -1;
    if(tj.AlgMod[kKilled]) return tjtp;
    // stash the ID in the Step
    tjtp.Step = tj.ID;
    tjtp.CTP = tj.CTP;
    tjtp.Pos[0] = 0;
    tjtp.Pos[1] = 0;
    tjtp.Dir[0] = 0;
    tjtp.Dir[1] = 0;
    float cnt = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if(tp.Chg <= 0) continue;
      tjtp.Pos[0] += tp.Pos[0];
      tjtp.Pos[1] += tp.Pos[1];
      tjtp.Dir[1] += tp.Dir[1];
      ++cnt;
    } // ipt
    tjtp.Pos[0] /= cnt;
    tjtp.Pos[1] /= cnt;
    tjtp.Dir[1] /= cnt;
    double arg = 1 - tjtp.Dir[1] * tjtp.Dir[1];
    if(arg < 0) arg = 0;
    tjtp.Dir[0] = sqrt(arg);
    tjtp.Ang = atan2(tjtp.Dir[1], tjtp.Dir[0]);
    tjtp.Chg = 1;
    return tjtp;
  } // CreateTjTP

  //////////////////////////////////////////
  void MaskTrajEndPoints(TCSlice& slc, Trajectory& tj, unsigned short nPts)
  {

    // Masks off (sets all hits not-Used) nPts trajectory points at the leading edge of the
    // trajectory, presumably because the fit including this points is poor. The position, direction
    // and Delta of the last nPts points is updated as well

    if(tj.EndFlag[1][kEndKink]) return;
    if(tj.Pts.size() < 3) {
      mf::LogError("TC")<<"MaskTrajEndPoints: Trajectory ID "<<tj.ID<<" too short to mask hits ";
      return;
    }
    if(nPts > tj.Pts.size() - 2) {
      mf::LogError("TC")<<"MaskTrajEndPoints: Trying to mask too many points "<<nPts<<" Pts.size "<<tj.Pts.size();
      return;
    }

    // find the last good point (with charge)
    unsigned short lastGoodPt = USHRT_MAX ;

    if(tcc.dbgStp) {
      mf::LogVerbatim("TC")<<"MTEP: lastGoodPt "<<lastGoodPt<<" Pts size "<<tj.Pts.size()<<" tj.IsGood "<<tj.IsGood;
    }
    if(lastGoodPt == USHRT_MAX) return;
    tj.EndPt[1] = lastGoodPt;

    //for(unsigned short ii = 0; ii < nPts; ++ii) {
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      if (ipt==lastGoodPt) break;
      UnsetUsedHits(slc, tj.Pts[ipt]);
      // Reset the position and direction of the masked off points
      tj.Pts[ipt].Dir = tj.Pts[lastGoodPt].Dir;
      if(tj.Pts[lastGoodPt].AngleCode == 2) {
        // Very large angle: Move by path length
        float path = TrajPointSeparation(tj.Pts[lastGoodPt], tj.Pts[ipt]);
        tj.Pts[ipt].Pos[0] = tj.Pts[lastGoodPt].Pos[0] + path * tj.Pts[ipt].Dir[0];
        tj.Pts[ipt].Pos[1] = tj.Pts[lastGoodPt].Pos[1] + path * tj.Pts[ipt].Dir[1];
      } else {
        // Not large angle: Move by wire
        float dw = tj.Pts[ipt].Pos[0] - tj.Pts[lastGoodPt].Pos[0];
        // Correct the projected time to the wire
        float newpos = tj.Pts[lastGoodPt].Pos[1] + dw * tj.Pts[ipt].Dir[1] / tj.Pts[ipt].Dir[0];
        if(tcc.dbgStp) mf::LogVerbatim("TC")<<"MTEP: ipt "<<ipt<<" Pos[0] "<<tj.Pts[ipt].Pos[0]<<". Move Pos[1] from "<<tj.Pts[ipt].Pos[1]<<" to "<<newpos;
        tj.Pts[ipt].Pos[1] = tj.Pts[lastGoodPt].Pos[1] + dw * tj.Pts[ipt].Dir[1] / tj.Pts[ipt].Dir[0];
      }
      tj.Pts[ipt].Delta = PointTrajDOCA(slc, tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tj.Pts[ipt]);
      if(tcc.dbgStp) mf::LogVerbatim("TC")<<" masked ipt "<<ipt<<" Pos "<<PrintPos(slc, tj.Pts[ipt])<<" Chg "<<tj.Pts[ipt].Chg;
    } // ii
    SetEndPoints(tj);

  } // MaskTrajEndPoints

} // namespace tca
