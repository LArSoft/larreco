#include "larreco/RecoAlg/TCAlg/Utils.h"

#include <boost/algorithm/string/classification.hpp>   // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp>            // Include for boost::split

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h" // for Hit
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
#include "larreco/RecoAlg/TCAlg/StepUtils.h"
#include "larreco/RecoAlg/TCAlg/PostStepUtils.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h" // for tcc

#include <array>
#include <stdlib.h>
#include <string>
#include <vector>

bool valsDecreasing (const SortEntry& c1, const SortEntry& c2) { return c1.val > c2.val;}
bool valsIncreasing (const SortEntry& c1, const SortEntry& c2) { return c1.val < c2.val;}

namespace tca {

  // dressed muons
  void
  MakeHaloTj(TCSlice& slc, Trajectory& muTj, bool prt)
  {
    // Creates a "halo trajectory" around a muon tj consisting of hits and trajectories
    // that are within MuonTag[4] distance. The halo tj is a virtual clone of muTj in the
    // sense that it has the same number of points and the same start and end points.

    if (tcc.muonTag.size() < 5) return;
    if (tcc.muonTag[4] <= 0) return;
    if (!tcc.useAlg[kHaloTj]) return;

    if (muTj.PDGCode != 13) return;

    // check for daughter delta-rays
    std::vector<int> dtrs;
    for (auto& dtj : slc.tjs) {
      if (dtj.AlgMod[kKilled]) continue;
      if (dtj.ParentID != muTj.ID) continue;
      dtrs.push_back(dtj.ID);
      if (!dtj.AlgMod[kDeltaRay]) continue;
      if (prt) mf::LogVerbatim("TC") << "MakeHaloTj: Killing delta-ray T" << dtj.ID;
      // Kill a delta-ray PFParticle?
      if (dtj.AlgMod[kMat3D]) {
        unsigned short pfpIndex = GetPFPIndex(slc, dtj.ID);
        if (pfpIndex == USHRT_MAX) {
          if (prt) mf::LogVerbatim("TC") << " No PFP found for 3D-matched delta-ray";
        }
        else {
          auto& pfp = slc.pfps[pfpIndex];
          if (prt) mf::LogVerbatim("TC") << " Killing delta-ray PFParticle P" << pfp.UID;
          pfp.ID = 0;
          // correct the parent -> daughter assn
          if (pfp.ParentUID > 0) {
            auto parentIndx = GetSliceIndex("P", pfp.ParentUID);
            if (parentIndx.first != USHRT_MAX) {
              auto& parent = slices[parentIndx.first].pfps[parentIndx.second];
              std::vector<int> newDtrUIDs;
              for (auto uid : parent.DtrUIDs)
                if (uid != dtj.UID) newDtrUIDs.push_back(uid);
              parent.DtrUIDs = newDtrUIDs;
            } // parent found
          }   // correct the parent
        }     // kill PFParticle
      }       // kill
      MakeTrajectoryObsolete(slc, (unsigned int)(dtj.ID - 1));
    } // dtj

    // make a copy
    Trajectory tj;
    tj.CTP = muTj.CTP;
    // We can't use StoreTraj so variables need to be defined here
    tj.ID = slc.tjs.size() + 1;
    tj.WorkID = muTj.WorkID;
    // increment the global ID
    ++evt.globalT_UID;
    tj.UID = evt.globalT_UID;
    tj.PDGCode = 11;
    tj.Pass = muTj.Pass;
    tj.StepDir = muTj.StepDir;
    tj.StartEnd = muTj.StartEnd;
    tj.TotChg = 0;
    tj.ChgRMS = 0;
    tj.EndPt[0] = 0;
    tj.ParentID = muTj.ID;
    tj.AlgMod.reset();
    tj.AlgMod[kHaloTj] = true;
    // start a list of tjs that have points near the muon
    std::vector<int> closeTjs;
    for (unsigned short ipt = muTj.EndPt[0]; ipt <= muTj.EndPt[1]; ++ipt) {
      auto tp = muTj.Pts[ipt];
      tp.Hits.resize(0);
      tp.UseHit.reset();
      tp.Chg = 0;
      tp.AveChg = 0;
      tp.ChgPull = 0;
      tp.Delta = 0;
      tp.DeltaRMS = 0;
      tp.FitChi = 0;
      tp.NTPsFit = 0;
      float window = tcc.muonTag[4];
      if (tp.Dir[0] != 0) window *= std::abs(1 / tp.Dir[0]);
      if (!FindCloseHits(slc, tp, window, kAllHits)) continue;
      // add unused hits to the point and look for close tjs
      bool hitsAdded = false;
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        auto inTraj = slc.slHits[iht].InTraj;
        if (inTraj < 0) continue;
        if (inTraj == 0) {
          tp.UseHit[ii] = true;
          slc.slHits[iht].InTraj = tj.ID;
          hitsAdded = true;
        }
        else {
          // add to the closeTjs list
          if (inTraj != muTj.ID &&
              std::find(closeTjs.begin(), closeTjs.end(), inTraj) == closeTjs.end())
            closeTjs.push_back(inTraj);
        }
      } // ii
      if (hitsAdded) {
        DefineHitPos(slc, tp);
        tp.Delta = PointTrajDOCA(slc, tp.HitPos[0], tp.HitPos[1], tp);
        tj.TotChg += tp.Chg;
        tj.Pts.push_back(tp);
      } // hitsAdded
    }   // ipt
    if (tj.Pts.empty()) return;
    tj.EndPt[1] = tj.Pts.size() - 1;
    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "MHTj: T" << muTj.ID << " npts " << tj.Pts.size() << " close";
      for (auto tid : closeTjs)
        myprt << " T" << tid;
      myprt << "\n";
      PrintTrajectory("DM", slc, tj, USHRT_MAX);
    }
    slc.tjs.push_back(tj);
  } // MakeHaloTj

  /////////////////////////////////////////
  float
  MaxChargeAsymmetry(TCSlice& slc, std::vector<int>& tjIDs)
  {
    // calculates the maximum charge asymmetry in all planes using the supplied list of Tjs
    if (tjIDs.size() < 2) return 1;
    std::vector<float> plnchg(slc.nPlanes);
    for (auto tjid : tjIDs) {
      if (tjid <= 0 || tjid > (int)slc.tjs.size()) return 1;
      auto& tj = slc.tjs[tjid - 1];
      if (tj.TotChg == 0) UpdateTjChgProperties("MCA", slc, tj, false);
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      plnchg[plane] += tj.TotChg;
    } // tjid
    float aveChg = 0;
    float cnt = 0;
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      if (plnchg[plane] == 0) continue;
      aveChg += plnchg[plane];
      ++cnt;
    } // plane
    if (cnt < 2) return 1;
    aveChg /= cnt;
    float maxAsym = 0;
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      // ignore zeros
      if (plnchg[plane] == 0) continue;
      float asym = std::abs(plnchg[plane] - aveChg) / (plnchg[plane] + aveChg);
      if (asym > maxAsym) maxAsym = asym;
    } // plane
    return maxAsym;
  } // MaxChargeAsymmetry

  /////////////////////////////////////////
  int
  PDGCodeVote(const TCSlice& slc, const std::vector<int>& tjIDs)
  {
    // Returns the most likely PDGCode for the set of Tjs provided
    // The PDG codes are:
    // 0 = your basic track-like trajectory
    // 11 = Tagged delta-ray
    // 13 = Tagged muon
    // 211 = pion-like. There exists a Bragg peak at an end with a vertex
    // 2212 = proton-like. There exists a Bragg peak at an end without a vertex
    std::array<int, 5> codeList = {{0, 11, 13, 111, 211}};
    unsigned short codeIndex = 0;
    if (tjIDs.empty()) return codeList[codeIndex];

    std::array<unsigned short, 5> cnts;
    cnts.fill(0);
    float maxLen = 0;
    for (auto tjid : tjIDs) {
      if (tjid <= 0 || tjid > (int)slc.tjs.size()) continue;
      auto& tj = slc.tjs[tjid - 1];
      for (unsigned short ii = 0; ii < 5; ++ii)
        if (tj.PDGCode == codeList[ii]) ++cnts[ii];
      float len = TrajLength(tj);
      if (len > maxLen) maxLen = len;
    } // tjid
    unsigned maxCnt = 0;
    // ignore the first PDG code in the list (the default)
    for (unsigned short ii = 1; ii < 5; ++ii) {
      if (cnts[ii] > maxCnt) {
        maxCnt = cnts[ii];
        codeIndex = ii;
      }
    } // ii
    return codeList[codeIndex];
  } // PDGCodeVote

  /////////////////////////////////////////
  int
  NeutrinoPrimaryTjID(const TCSlice& slc, const Trajectory& tj)
  {
    // Returns the ID of the grandparent of this tj that is a primary tj that is attached
    // to the neutrino vertex. 0 is returned if this condition is not met.
    if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) return -1;
    if (tj.ParentID <= 0) return -1;
    int primID = PrimaryID(slc, tj);
    if (primID <= 0 || primID > (int)slc.tjs.size()) return -1;

    // We have the ID of the primary tj. Now see if it is attached to the neutrino vertex
    auto& ptj = slc.tjs[primID - 1];
    for (unsigned short end = 0; end < 2; ++end) {
      if (ptj.VtxID[end] == 0) continue;
      auto& vx2 = slc.vtxs[ptj.VtxID[end] - 1];
      if (vx2.Vx3ID == 0) continue;
      auto& vx3 = slc.vtx3s[vx2.Vx3ID - 1];
      if (vx3.Neutrino) return primID;
    } // end
    return -1;
  } // NeutrinoPrimaryTjUID

  /////////////////////////////////////////
  int
  PrimaryID(const TCSlice& slc, const Trajectory& tj)
  {
    // Returns the ID of the grandparent trajectory of this trajectory that is a primary
    // trajectory (i.e. whose ParentID = 0).
    if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) return -1;
    if (tj.ParentID < 0 || tj.ParentID > (int)slc.tjs.size()) return -1;
    if (tj.ParentID == 0) return tj.ID;
    int parid = tj.ParentID;
    for (unsigned short nit = 0; nit < 10; ++nit) {
      if (parid < 1 || parid > (int)slc.tjs.size()) break;
      auto& tj = slc.tjs[parid - 1];
      if (tj.ParentID < 0 || tj.ParentID > (int)slc.tjs.size()) return -1;
      if (tj.ParentID == 0) return tj.ID;
      parid = tj.ParentID;
    } // nit
    return -1;
  } // PrimaryID

  /////////////////////////////////////////
  int
  PrimaryUID(const TCSlice& slc, const PFPStruct& pfp)
  {
    // returns the UID of the most upstream PFParticle (that is not a neutrino)

    if (int(pfp.ParentUID) == pfp.UID || pfp.ParentUID <= 0) return pfp.ID;
    int paruid = pfp.ParentUID;
    int dtruid = pfp.UID;
    unsigned short nit = 0;
    while (true) {
      auto slcIndx = GetSliceIndex("P", paruid);
      auto& parent = slices[slcIndx.first].pfps[slcIndx.second];
      // found a neutrino
      if (parent.PDGCode == 14 || parent.PDGCode == 12) return dtruid;
      // found a primary PFParticle?
      if (parent.ParentUID == 0) return parent.UID;
      if (int(parent.ParentUID) == parent.UID) return parent.UID;
      dtruid = parent.UID;
      paruid = parent.ParentUID;
      if (paruid < 0) return 0;
      ++nit;
      if (nit == 10) return 0;
    }
  } // PrimaryUID

  /////////////////////////////////////////
  float
  PointPull(TCSlice& slc, Point2_t pos, float chg, const Trajectory& tj)
  {
    // returns the combined position and charge pull for the charge at pos
    // relative to the Tj closest to that point using a loose requirement on position separation.
    if (tj.AlgMod[kKilled]) return 100;
    if (tj.AveChg <= 0) return 100;
    // find the closest point on the tj to pos
    unsigned short closePt = USHRT_MAX;
    float close = 1000;
    for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      float sep2 = PosSep2(pos, tp.Pos);
      if (sep2 > close) continue;
      close = sep2;
      closePt = ipt;
    } // ipt
    if (closePt == USHRT_MAX) return 100;
    // find the delta between the projection of the Tj close TP to inTP
    auto& tp = tj.Pts[closePt];
    float delta = PointTrajDOCA(slc, pos[0], pos[1], tp);
    // estimate the proejcted position error (roughly)
    float posErr = tp.DeltaRMS;
    if (tp.AngErr > 0 && close > 10) posErr += sqrt(tp.AngErr * sqrt(close));
    if (posErr < 0.1) posErr = 0.1;
    float posPull = delta / posErr;
    float chgErr = tj.ChgRMS;
    if (chgErr < 0.15) chgErr = 0.15;
    float chgPull = std::abs(chg / tj.AveChg - 1) / chgErr;
    // return a simple average
    return 0.5 * (posPull + chgPull);
  } // PointPull

  /////////////////////////////////////////
  bool
  CompatibleMerge(const TCSlice& slc, std::vector<int>& tjIDs, bool prt)
  {
    // Returns true if the last Tj in tjIDs has a topology consistent with it being
    // merged with other Tjs in the same plane in the list. This is done by requiring that
    // the closest TP between the last Tj and any other Tj is EndPt[0] or EndPt[1]. This is
    // shown graphically here where the numbers represent the ID of a Tj that has a TP on a wire.
    // Assume that TjIDs = {1, 2, 3, 4, 7} where T1 and T3 are in plane 0, T2 is in plane 1 and
    // T4 is in plane 2. T7, in plane 0, was added to TjIDs with the intent of merging it with
    // T1 and T3 into a single trajectory. This is a compatible merge if Tj7 has the following
    // topology:
    //  111111 333333 7777777
    // This is an incompatible topology
    //  111111 333333
    //      7777777777
    if (tjIDs.size() < 2) return false;
    unsigned short lasttj = tjIDs[tjIDs.size() - 1] - 1;
    auto& mtj = slc.tjs[lasttj];
    bool mtjIsShort = (mtj.Pts.size() < 5);
    // minimum separation from each end of mtj
    std::array<float, 2> minsep2{{1000, 1000}};
    // ID of the Tj with the minimum separation
    std::array<int, 2> minsepTj{{0, 0}};
    // and the index of the point on that Tj
    std::array<unsigned short, 2> minsepPt;
    // determine the end of the closest Tj point. Start by assuming
    // the closest Tj point is not near an end (end = 0);
    std::array<unsigned short, 2> minsepEnd;
    for (auto tjid : tjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      if (tj.CTP != mtj.CTP) continue;
      if (tj.ID == mtj.ID) continue;
      for (unsigned short mend = 0; mend < 2; ++mend) {
        Point2_t mendPos = mtj.Pts[mtj.EndPt[mend]].Pos;
        float sep2 = minsep2[mend];
        unsigned short closePt = 0;
        if (!TrajClosestApproach(tj, mendPos[0], mendPos[1], closePt, sep2)) continue;
        minsep2[mend] = sep2;
        minsepTj[mend] = tjid;
        minsepPt[mend] = closePt;
        // set the end to a bogus value (not near an end)
        minsepEnd[mend] = 2;
        short dend0 = abs((short)closePt - tj.EndPt[0]);
        short dend1 = abs((short)closePt - tj.EndPt[1]);
        if (dend0 < dend1 && dend0 < 3) minsepEnd[mend] = 0;
        if (dend1 < dend0 && dend1 < 3) minsepEnd[mend] = 1;
      } // mend
    }   // tjid
    // don't require that the minsepTjs be the same. This would reject this topology
    //  111111 333333 7777777
    // if mtj.ID = 3
    bool isCompatible = (minsepEnd[0] != 2 && minsepEnd[1] != 2);
    // check for large separation between the closest points for short Tjs
    if (isCompatible && mtjIsShort) {
      float minminsep = minsep2[0];
      if (minsep2[1] < minminsep) minminsep = minsep2[1];
      // require that the separation be less than sqrt(5)
      isCompatible = minminsep < 5;
    }
    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "CompatibleMerge: T" << mtj.ID << " end";
      for (unsigned short end = 0; end < 2; ++end)
        myprt << " T" << minsepTj[end] << "_I" << minsepPt[end] << "_E" << minsepEnd[end]
              << " minsep " << sqrt(minsep2[end]);
      myprt << " Compatible? " << isCompatible;
    } // prt
    return isCompatible;

  } // CompatibleMerge

  /////////////////////////////////////////
  bool
  CompatibleMerge(const TCSlice& slc, const Trajectory& tj1, const Trajectory& tj2, bool prt)
  {
    // returns true if the two Tjs are compatible with and end0-end1 merge. This function has many aspects of the
    // compatibility checks done in EndMerge but with looser cuts.
    if (tj1.AlgMod[kKilled] || tj2.AlgMod[kKilled]) return false;
    if (tj1.AlgMod[kHaloTj] || tj2.AlgMod[kHaloTj]) return false;
    if (tj1.CTP != tj2.CTP) return false;
    unsigned short end1 = -1, end2 = 0;
    float minLen = PosSep(tj1.Pts[tj1.EndPt[0]].Pos, tj1.Pts[tj1.EndPt[1]].Pos);
    float len2 = PosSep(tj2.Pts[tj2.EndPt[0]].Pos, tj2.Pts[tj2.EndPt[1]].Pos);
    if (len2 < minLen) minLen = len2;
    minLen *= 1.2;
    if (minLen > 10) minLen = 10;
    for (unsigned short e1 = 0; e1 < 2; ++e1) {
      auto& tp1 = tj1.Pts[tj1.EndPt[e1]];
      for (unsigned short e2 = 0; e2 < 2; ++e2) {
        auto& tp2 = tj2.Pts[tj2.EndPt[e2]];
        float sep = PosSep(tp1.Pos, tp2.Pos);
        if (sep < minLen) {
          minLen = sep;
          end1 = e1;
          end2 = e2;
        }
      } // e2
    }   // e1
    if (end1 < 0) return false;
    // require end to end
    if (end2 != 1 - end1) return false;

    float overlapFraction = OverlapFraction(slc, tj1, tj2);
    if (overlapFraction > 0.25) {
      if (prt)
        mf::LogVerbatim("TC") << "CM: " << tj1.ID << " " << tj2.ID << " overlapFraction "
                              << overlapFraction << " > 0.25 ";
      return false;
    }

    auto& tp1 = tj1.Pts[tj1.EndPt[end1]];
    auto& tp2 = tj2.Pts[tj2.EndPt[end2]];
    float doca1 = PointTrajDOCA(slc, tp1.Pos[0], tp1.Pos[1], tp2);
    float doca2 = PointTrajDOCA(slc, tp2.Pos[0], tp2.Pos[1], tp1);
    if (doca1 > 2 && doca2 > 2) {
      if (prt)
        mf::LogVerbatim("TC") << "CM: " << tj1.ID << " " << tj2.ID << " Both docas > 2 " << doca1
                              << " " << doca2;
      return false;
    }

    float dang = DeltaAngle(tp1.Ang, tp2.Ang);
    if (dang > 2 * tcc.kinkCuts[0]) {
      if (prt)
        mf::LogVerbatim("TC") << "CM: " << tj1.ID << " " << tj2.ID << " dang " << dang << " > "
                              << 2 * tcc.kinkCuts[0];
      return false;
    }

    return true;
  } // CompatibleMerge

  /////////////////////////////////////////
  float
  OverlapFraction(const TCSlice& slc, const Trajectory& tj1, const Trajectory& tj2)
  {
    // returns the fraction of wires spanned by two trajectories
    float minWire = 1E6;
    float maxWire = -1E6;

    float cnt1 = 0;
    for (auto& tp : tj1.Pts) {
      if (tp.Chg == 0) continue;
      if (tp.Pos[0] < 0) continue;
      if (tp.Pos[0] < minWire) minWire = tp.Pos[0];
      if (tp.Pos[0] > maxWire) maxWire = tp.Pos[0];
      ++cnt1;
    }
    if (cnt1 == 0) return 0;
    float cnt2 = 0;
    for (auto& tp : tj2.Pts) {
      if (tp.Chg == 0) continue;
      if (tp.Pos[0] < 0) continue;
      if (tp.Pos[0] < minWire) minWire = tp.Pos[0];
      if (tp.Pos[0] > maxWire) maxWire = tp.Pos[0];
      ++cnt2;
    }
    if (cnt2 == 0) return 0;
    int span = maxWire - minWire;
    if (span <= 0) return 0;
    std::vector<unsigned short> wcnt(span);
    for (auto& tp : tj1.Pts) {
      if (tp.Chg == 0) continue;
      if (tp.Pos[0] < -0.4) continue;
      int indx = std::nearbyint(tp.Pos[0] - minWire);
      if (indx < 0 || indx > span - 1) continue;
      ++wcnt[indx];
    }
    for (auto& tp : tj2.Pts) {
      if (tp.Chg == 0) continue;
      if (tp.Pos[0] < -0.4) continue;
      int indx = std::nearbyint(tp.Pos[0] - minWire);
      if (indx < 0 || indx > span - 1) continue;
      ++wcnt[indx];
    }
    float cntOverlap = 0;
    for (auto cnt : wcnt)
      if (cnt > 1) ++cntOverlap;
    if (cnt1 < cnt2) { return cntOverlap / cnt1; }
    else {
      return cntOverlap / cnt2;
    }

  } // OverlapFraction

  /////////////////////////////////////////
  unsigned short
  AngleRange(TrajPoint const& tp)
  {
    return AngleRange(tp.Ang);
  }

  /////////////////////////////////////////
  void
  SetAngleCode(TrajPoint& tp)
  {
    unsigned short ar = AngleRange(tp.Ang);
    if (ar == tcc.angleRanges.size() - 1) {
      // Very large angle
      tp.AngleCode = 2;
    }
    else if (tcc.angleRanges.size() > 2 && ar == tcc.angleRanges.size() - 2) {
      // Large angle
      tp.AngleCode = 1;
    }
    else {
      // Small angle
      tp.AngleCode = 0;
    }

  } // SetAngleCode

  /////////////////////////////////////////
  unsigned short
  AngleRange(float angle)
  {
    // returns the index of the angle range
    if (angle > M_PI) angle = M_PI;
    if (angle < -M_PI) angle = M_PI;
    if (angle < 0) angle = -angle;
    if (angle > M_PI / 2) angle = M_PI - angle;
    for (unsigned short ir = 0; ir < tcc.angleRanges.size(); ++ir) {
      if (angle < tcc.angleRanges[ir]) return ir;
    }
    return tcc.angleRanges.size() - 1;
  } // AngleRange

  //////////////////////////////////////////
  void
  FitTraj(TCSlice& slc, Trajectory& tj)
  {
    // Jacket around FitTraj to fit the leading edge of the supplied trajectory
    unsigned short originPt = tj.EndPt[1];
    unsigned short npts = tj.Pts[originPt].NTPsFit;
    TrajPoint tpFit;
    unsigned short fitDir = -1;
    FitTraj(slc, tj, originPt, npts, fitDir, tpFit);
    tj.Pts[originPt] = tpFit;

  } // FitTraj

  //////////////////////////////////////////
  void
  FitTraj(TCSlice& slc,
          Trajectory& tj,
          unsigned short originPt,
          unsigned short npts,
          short fitDir,
          TrajPoint& tpFit)
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

    if (originPt > tj.Pts.size() - 1) return;

    // copy the origin TP into the fit TP
    tpFit = tj.Pts[originPt];
    // Assume that the fit will fail
    tpFit.FitChi = 999;
    if (fitDir < -1 || fitDir > 1) return;

    std::vector<double> x, y;
    Point2_t origin = tj.Pts[originPt].HitPos;
    // Use TP position if there aren't any hits on it
    if (tj.Pts[originPt].Chg == 0) origin = tj.Pts[originPt].Pos;

    // simple two point case
    if (NumPtsWithCharge(slc, tj, false) == 2) {
      for (unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
        if (tj.Pts[ipt].Chg <= 0) continue;
        double xx = tj.Pts[ipt].HitPos[0] - origin[0];
        double yy = tj.Pts[ipt].HitPos[1] - origin[1];
        x.push_back(xx);
        y.push_back(yy);
      } // ii
      if (x.size() != 2) return;
      if (x[0] == x[1]) {
        // Either + or - pi/2
        tpFit.Ang = M_PI / 2;
        if (y[1] < y[0]) tpFit.Ang = -tpFit.Ang;
      }
      else {
        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        tpFit.Ang = atan2(dy, dx);
      }
      tpFit.Dir[0] = cos(tpFit.Ang);
      tpFit.Dir[1] = sin(tpFit.Ang);
      tpFit.Pos[0] += origin[0];
      tpFit.Pos[1] += origin[1];
      tpFit.AngErr = 0.01;
      tpFit.FitChi = 0.01;
      SetAngleCode(tpFit);
      return;
    } // two points

    std::vector<double> w, q;
    std::array<double, 2> dir;
    double xx, yy, xr, yr;
    double chgWt;

    // Rotate the traj hit position into the coordinate system defined by the
    // originPt traj point, where x = along the trajectory, y = transverse
    double rotAngle = tj.Pts[originPt].Ang;
    double cs = cos(-rotAngle);
    double sn = sin(-rotAngle);

    // enter the originPT hit info if it exists
    if (tj.Pts[originPt].Chg > 0) {
      xx = tj.Pts[originPt].HitPos[0] - origin[0];
      yy = tj.Pts[originPt].HitPos[1] - origin[1];
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      x.push_back(xr);
      y.push_back(yr);
      chgWt = tj.Pts[originPt].ChgPull;
      if (chgWt < 1) chgWt = 1;
      chgWt *= chgWt;
      w.push_back(chgWt * tj.Pts[originPt].HitPosErr2);
    }

    // correct npts to account for the origin point
    if (fitDir != 0) --npts;

    // step in the + direction first
    if (fitDir != -1) {
      unsigned short cnt = 0;
      for (unsigned short ipt = originPt + 1; ipt < tj.Pts.size(); ++ipt) {
        if (tj.Pts[ipt].Chg <= 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if (chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if (cnt == npts) break;
      } // ipt
    }   // fitDir != -1

    // step in the - direction next
    if (fitDir != 1 && originPt > 0) {
      unsigned short cnt = 0;
      for (unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = originPt - ii;
        if (ipt > tj.Pts.size() - 1) continue;
        if (tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if (chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if (cnt == npts) break;
        if (ipt == 0) break;
      } // ipt
    }   // fitDir != -1

    // Not enough points to define a line?
    if (x.size() < 2) return;

    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;

    // weight by the charge ratio and accumulate sums
    double wght;
    for (unsigned short ipt = 0; ipt < x.size(); ++ipt) {
      if (w[ipt] < 0.00001) w[ipt] = 0.00001;
      wght = 1 / w[ipt];
      sum += wght;
      sumx += wght * x[ipt];
      sumy += wght * y[ipt];
      sumx2 += wght * x[ipt] * x[ipt];
      sumy2 += wght * y[ipt] * y[ipt];
      sumxy += wght * x[ipt] * y[ipt];
    }
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if (delta == 0) return;
    // A is the intercept
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // B is the slope
    double B = (sumxy * sum - sumx * sumy) / delta;

    // The chisq will be set below if there are enough points. Don't allow it to be 0
    // so we can take Chisq ratios later
    tpFit.FitChi = 0.01;
    double newang = atan(B);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAngle);
    sn = sin(rotAngle);
    tpFit.Dir[0] = cs * dir[0] - sn * dir[1];
    tpFit.Dir[1] = sn * dir[0] + cs * dir[1];
    // ensure that the direction is consistent with the originPt direction
    bool flipDir = false;
    if (AngleRange(tj.Pts[originPt]) > 0) {
      flipDir = std::signbit(tpFit.Dir[1]) != std::signbit(tj.Pts[originPt].Dir[1]);
    }
    else {
      flipDir = std::signbit(tpFit.Dir[0]) != std::signbit(tj.Pts[originPt].Dir[0]);
    }
    if (flipDir) {
      tpFit.Dir[0] = -tpFit.Dir[0];
      tpFit.Dir[1] = -tpFit.Dir[1];
    }
    tpFit.Ang = atan2(tpFit.Dir[1], tpFit.Dir[0]);
    SetAngleCode(tpFit);

    // rotate (0, intcpt) into (W,T) coordinates
    tpFit.Pos[0] = -sn * A + origin[0];
    tpFit.Pos[1] = cs * A + origin[1];
    // force the origin to be at origin[0]
    if (tpFit.AngleCode < 2) MoveTPToWire(tpFit, origin[0]);

    if (x.size() < 3) return;

    // Calculate chisq/DOF
    double ndof = x.size() - 2;
    double varnce =
      (sumy2 + A * A * sum + B * B * sumx2 - 2 * (A * sumy + B * sumxy - A * B * sumx)) / ndof;
    if (varnce > 0.) {
      // Intercept error is not used
      //      InterceptError = sqrt(varnce * sumx2 / delta);
      double slopeError = sqrt(varnce * sum / delta);
      tpFit.AngErr = std::abs(atan(slopeError));
    }
    else {
      tpFit.AngErr = 0.01;
    }
    sum = 0;
    // calculate chisq
    double arg;
    for (unsigned short ii = 0; ii < y.size(); ++ii) {
      arg = y[ii] - A - B * x[ii];
      sum += arg * arg / w[ii];
    }
    tpFit.FitChi = sum / ndof;
  } // FitTraj

  ////////////////////////////////////////////////
  unsigned short
  GetPFPIndex(const TCSlice& slc, int tjID)
  {
    if (slc.pfps.empty()) return USHRT_MAX;
    for (unsigned int ipfp = 0; ipfp < slc.pfps.size(); ++ipfp) {
      const auto& pfp = slc.pfps[ipfp];
      if (std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjID) != pfp.TjIDs.end()) return ipfp;
    } // indx
    return USHRT_MAX;
  } // GetPFPIndex

  ////////////////////////////////////////////////
  void
  ReleaseHits(TCSlice& slc, Trajectory& tj)
  {
    // Sets InTraj[] = 0 for all TPs in work. Called when abandoning work
    for (auto& tp : tj.Pts) {
      for (auto iht : tp.Hits) {
        if (slc.slHits[iht].InTraj == tj.ID) slc.slHits[iht].InTraj = 0;
      }
    } // tp

  } // ReleaseWorkHits

  //////////////////////////////////////////
  void
  UnsetUsedHits(TCSlice& slc, TrajPoint& tp)
  {
    // Sets InTraj = 0 and UseHit false for all used hits in tp
    for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if (tp.UseHit[ii]) {
        slc.slHits[tp.Hits[ii]].InTraj = 0;
        tp.UseHit[ii] = false;
      } // UseHit
    }   // ii
    tp.Chg = 0;
  } // UnsetUsedHits

  ////////////////////////////////////////////////
  bool
  StoreTraj(TCSlice& slc, Trajectory& tj)
  {

    // check for errors
    if(!tj.IsGood) return false;
    if(tj.Pts.size() < 2) return false;
    for (auto& tp : tj.Pts) {
      if (tp.Hits.size() > 16) return false;
    } // tp

    // This shouldn't be necessary but do it anyway
    SetEndPoints(tj);

    if (tj.NeedsUpdate) UpdateTjChgProperties("ST", slc, tj, false);
    // update the flags at both ends
    SetEndFlags(slc, tj, 2);
    if (slc.tjs.size() >= USHRT_MAX || tj.EndPt[1] <= tj.EndPt[0] || tj.EndPt[1] > tj.Pts.size()) {
      ReleaseHits(slc, tj);
      return false;
    }

    unsigned short npts = tj.EndPt[1] - tj.EndPt[0] + 1;
    if (npts < 2) return false;

    auto& endTp0 = tj.Pts[tj.EndPt[0]];
    auto& endTp1 = tj.Pts[tj.EndPt[1]];

    // ensure that angle errors are defined at both ends, ignoring junk Tjs
    if (!tj.AlgMod[kJunkTj]) {
      if (endTp0.AngErr == 0.1 && endTp1.AngErr != 0.1) { endTp0.AngErr = endTp1.AngErr; }
      else if (endTp0.AngErr != 0.1 && endTp1.AngErr == 0.1) {
        endTp1.AngErr = endTp0.AngErr;
      }
    } // not a junk Tj

    // Calculate the charge near the end and beginning if necessary. This must be a short
    // trajectory. Find the average using 4 points
    if (endTp0.AveChg <= 0) {
      unsigned short cnt = 0;
      float sum = 0;
      for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        if (tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if (cnt == 4) break;
      }
      tj.Pts[tj.EndPt[0]].AveChg = sum / (float)cnt;
    }
    if (endTp1.AveChg <= 0 && npts < 5) endTp1.AveChg = endTp0.AveChg;
    if (endTp1.AveChg <= 0) {
      float sum = 0;
      unsigned short cnt = 0;
      for (unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        short ipt = tj.EndPt[1] - ii;
        if (ipt < 0) break;
        if (tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if (cnt == 4) break;
        if (ipt == 0) break;
      } // ii
      tj.Pts[tj.EndPt[1]].AveChg = sum / (float)cnt;
    } // begin charge == end charge

    // update the kink significance
    if (!tj.AlgMod[kJunkTj]) {
      unsigned short nPtsFit = tcc.kinkCuts[0];
      bool useChg = (tcc.kinkCuts[2] > 0);
      if (npts > 2 * nPtsFit) {
        for (unsigned short ipt = tj.EndPt[0] + nPtsFit; ipt < tj.EndPt[1] - nPtsFit; ++ipt) {
          auto& tp = tj.Pts[ipt];
          if (tp.KinkSig < 0) tp.KinkSig = KinkSignificance(slc, tj, ipt, nPtsFit, useChg, false);
        }
      } // long trajectory
    }   // not JunkTj

    UpdateTjChgProperties("ST", slc, tj, false);

    int trID = slc.tjs.size() + 1;

    // Define the Tj StartEnd. StartEnd = 0 means that the trajectory seems to be traveling
    // in the same order order as the points on the trajectory judging from the pattern of charge
    // at the beginning and the end
    float chg0 = 0, cnt0 = 0;
    float chg1 = 0, cnt1 = 0;
    unsigned short halfWay = (tj.EndPt[0] + tj.EndPt[1]) / 2;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      if(tj.Pts[ipt].Chg <= 0) continue;
      if(ipt < halfWay) {
        chg0 += tj.Pts[ipt].Chg;
        ++cnt0;
      } else {
        chg1 += tj.Pts[ipt].Chg;
        ++cnt1;
      }
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(tj.Pts[ipt].UseHit[ii]) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if (iht > slc.slHits.size() - 1) {
            ReleaseHits(slc, tj);
            return false;
          }
          if (slc.slHits[iht].InTraj > 0) {
            ReleaseHits(slc, tj);
            return false;
          } // error
          slc.slHits[iht].InTraj = trID;
        }
      } // ii
    } // ipt
    if(cnt0 > 0 && cnt1 > 0) {
      // assume that the start end is the end with the lower average charge
      chg0 /= cnt0;
      chg1 /= cnt1;
      if(chg1 > 1.2 * chg0) {
        tj.StartEnd = 0;
      } else if(chg0 > 1.2 * chg1) { 
        tj.StartEnd = 1;
      }
    } // valid cnt0 and cnt1

    // ensure that inTraj is clean for the ID
    for (unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      if (slc.slHits[iht].InTraj == tj.ID) slc.slHits[iht].InTraj = 0;
    } // iht

    tj.WorkID = tj.ID;
    tj.ID = trID;
    // increment the global ID
    ++evt.globalT_UID;
    tj.UID = evt.globalT_UID;
    // Don't clobber the ParentID if it was defined by the calling function
    if (tj.ParentID == 0) tj.ParentID = trID;
    slc.tjs.push_back(tj);
    if (tcc.modes[kModeDebug] && tcc.dbgSlc && debug.Hit != UINT_MAX) {
      // print some debug info
      for (unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
        for (unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if (slc.slHits[iht].allHitsIndex == debug.Hit) {
            std::cout << "Debug hit appears in trajectory w WorkID " << tj.WorkID << " UseHit "
                      << tj.Pts[ipt].UseHit[ii] << "\n";
          }
        } // ii
      }   // ipt
    }     // debug.Hit ...

    return true;
  } // StoreTraj

  ////////////////////////////////////////////////
  void SetEndFlags(TCSlice& slc, Trajectory& tj, unsigned short atEnd)
  {
    // Set the EndFlag for the trajectory at end = 0,1 or both ends if end > 1. Note
    // that some flags, e.g. kEndKink, kEndBragg are set elsewhere

    std::vector<unsigned short> ends;
    if(atEnd > 1) {
      ends.push_back(0);
      ends.push_back(1);
    } else {
      ends.push_back(atEnd);
    }
    unsigned short plane = DecodeCTP(tj.CTP).Plane;
    for(unsigned short itr = 0; itr < ends.size(); ++itr) {
      auto end = ends[itr];
      // set the no-nearby-hits-past-the-end flag to false
      tj.EndFlag[end][kHitsAfterEnd] = false;
      auto& tp = tj.Pts[tj.EndPt[end]];
      // ignore VLA TPs
      if(tp.AngleCode > 1) continue;
      // make a copy that we can move
      auto ntp = tj.Pts[tj.EndPt[end]];
      float wireIter = tj.StepDir;
      if(end == 0) wireIter *= -1;
      float nextWire = std::nearbyint(ntp.Pos[0]) + wireIter;
      if(nextWire <= 0 || (unsigned int)nextWire > slc.nWires[plane]) continue;
      MoveTPToWire(ntp, nextWire);
      // move another wire if this one is dead
      if(!evt.goodWire[plane][nextWire]) {
        nextWire += wireIter;
        if(nextWire <= 0 || (unsigned int)nextWire > slc.nWires[plane]) continue;
      } // not goodWire
      // Set the window to 2 WSE 
      float window = 2 / std::abs(ntp.Dir[0]);
      if (FindCloseHits(slc, ntp, window, kAllHits)) tj.EndFlag[end][kHitsAfterEnd] = true;
    } // itr

  } // SetEndFlags

  //////////////////////////////////////////
  void
  FitPar(const TCSlice& slc,
         const Trajectory& tj,
         unsigned short originPt,
         unsigned short npts,
         short fitDir,
         ParFit& pFit,
         unsigned short usePar)
  {
    // Fit a TP parameter, like Chg or Delta, to a line using the points starting at originPT.
    // Currently supported values of usePar are Chg (1) and Delta (2)

    pFit.ChiDOF = 999;
    pFit.AvePar = 0.;
    if (originPt > tj.Pts.size() - 1) return;
    if (fitDir != 1 && fitDir != -1) return;
    Point2_t inPt;
    Vector2_t outVec, outVecErr;
    float pErr, chiDOF;
    Fit2D(0, inPt, pErr, outVec, outVecErr, chiDOF);
    unsigned short cnt = 0;
    for (unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = originPt + ii * fitDir;
      if (ipt < tj.EndPt[0] || ipt > tj.EndPt[1]) break;
      auto& tp = tj.Pts[ipt];
      if (tp.Chg <= 0) continue;
      // Accumulate and save points
      inPt[0] = std::abs(tp.Pos[0] - tj.Pts[originPt].Pos[0]);
      float parVal = tp.Chg;
      // Assume errors are 10% for a charge fit
      pErr = 0.1 * parVal;
      if (usePar > 1) {
        parVal = tp.Delta;
        // use the TP hit position error for a Delta Fit
        pErr = sqrt(tp.HitPosErr2);
      }
      inPt[1] = parVal;
      pFit.AvePar += parVal;
      if (!Fit2D(2, inPt, pErr, outVec, outVecErr, chiDOF)) break;
      ++cnt;
      if (cnt == npts) break;
    } // ii
    if (cnt < npts) return;
    // do the fit and get the results
    if (!Fit2D(-1, inPt, pErr, outVec, outVecErr, chiDOF)) return;
    pFit.Pos = tj.Pts[originPt].Pos;
    pFit.Par0 = outVec[0];
    pFit.AvePar /= (float)cnt;
    pFit.ParErr = outVecErr[0];
    pFit.Pos = tj.Pts[originPt].Pos;
    pFit.ParSlp = outVec[1];
    pFit.ParSlpErr = outVecErr[1];
    pFit.ChiDOF = chiDOF;
    pFit.nPtsFit = cnt;
  } // FitPar
  
  //////////////////////////////////////////
  void
  ChkEndPtFit(TCSlice& slc, Trajectory& tj)
  {
    // Re-fit the end of the trajectory if it is a long track (made with loose cuts)
    // and the ChiDOF is high
    if(!tcc.useAlg[kNewCuts]) return;
    if(!tcc.useAlg[kEndPtFit]) return;
    if(tj.EndPt[1] - tj.EndPt[0] < 20) return;

    constexpr float maxChi = 1.5;
    if(tj.Pts[tj.EndPt[1]].FitChi < maxChi) return;

    unsigned short nPtsFit = 0;
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      short ipt = tj.EndPt[1] - ii;
      if(ipt < tj.EndPt[0] + 3) break;
      auto& tp = tj.Pts[ipt];
      if(tp.Chg <= 0) continue;
      ++nPtsFit;
      if(tp.FitChi < maxChi) break;
      // Don't re-fit with more than 20 points
      if(nPtsFit > 20) break;
    } // ii
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"ChkEndPtFit: Found "<<nPtsFit<<" points at end to re-fit";
    if(nPtsFit < 3) return;
    tj.Pts[tj.EndPt[1]].NTPsFit = nPtsFit;
    FitTraj(slc, tj);
    auto& lastTP = tj.Pts[tj.EndPt[1]];
    tj.AlgMod[kEndPtFit] = true;
    if(tcc.dbgStp) mf::LogVerbatim("TC")<<"  lastTP FitChi "<<tj.Pts[tj.EndPt[1]].FitChi;
    // update the last set of points
    for(unsigned short ipt = tj.EndPt[1] - nPtsFit; ipt < tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      tp.Dir = lastTP.Dir;
      tp.Ang = lastTP.Ang;
      tp.AngErr = lastTP.AngErr;
      tp.AngleCode = lastTP.AngleCode;
      float dw = tp.Pos[0] - lastTP.Pos[0];
      if(tp.Dir[0] != 0) tp.Pos[1] = lastTP.Pos[1] + dw * tp.Dir[1] / tp.Dir[0];
      tp.Delta = PointTrajDOCA(slc, tp.HitPos[0], tp.HitPos[1], tp);
      tp.DeltaRMS = lastTP.DeltaRMS;
      tp.NTPsFit = lastTP.NTPsFit;
      tp.FitChi = lastTP.FitChi;
      tp.AveChg = lastTP.AveChg;
      tp.ChgPull = (tp.Chg / tj.AveChg - 1) / tj.ChgRMS;
    } // ipt

  } // ChkEndPtFit

  //////////////////////////////////////////
  void
  TrimEndPts(std::string fcnLabel,
             TCSlice& slc,
             Trajectory& tj,
             const std::vector<float>& fQualityCuts,
             bool prt)
  {
    // Trim the hits off the end until there are at least MinPts consecutive hits at the end
    // and the fraction of hits on the trajectory exceeds fQualityCuts[0]
    // Minimum length requirement accounting for dead wires where - denotes a wire with a point
    // and D is a dead wire. Here is an example with minPts = 3
    //  ---DDDDD--- is OK
    //  ----DD-DD-- is OK
    //  ----DDD-D-- is OK
    //  ----DDDDD-- is not OK

    if (!tcc.useAlg[kTEP]) return;
    if (tj.PDGCode == 111) return;
    if (tj.EndFlag[1][kEndKink]) return;

    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    short minPts = fQualityCuts[1];
    if (minPts < 1) return;
    if (npwc < minPts) return;
    // don't consider short Tjs
    if (npwc < 8) return;

    // handle short tjs
    if (npwc == minPts + 1) {
      unsigned short endPt1 = tj.EndPt[1];
      auto& tp = tj.Pts[endPt1];
      auto& ptp = tj.Pts[endPt1 - 1];
      // remove the last point if the previous point has no charge or if
      // it isn't on the next wire
      float dwire = std::abs(ptp.Pos[0] - tp.Pos[0]);
      if (ptp.Chg == 0 || dwire > 1.1) {
        UnsetUsedHits(slc, tp);
        SetEndPoints(tj);
        tj.AlgMod[kTEP] = true;
      }
      return;
    } // short tj

    // find the separation between adjacent points, starting at the end
    short lastPt = 0;
    for (lastPt = tj.EndPt[1]; lastPt >= minPts; --lastPt) {
      // check for an error
      if (lastPt == 1) break;
      if (tj.Pts[lastPt].Chg == 0) continue;
      // number of points on adjacent wires
      unsigned short nadj = 0;
      unsigned short npwc = 0;
      for (short ipt = lastPt - minPts; ipt < lastPt; ++ipt) {
        if (ipt < 2) break;
        // the current point
        auto& tp = tj.Pts[ipt];
        // the previous point
        auto& ptp = tj.Pts[ipt - 1];
        if (tp.Chg > 0 && ptp.Chg > 0) {
          ++npwc;
          if (std::abs(tp.Pos[0] - ptp.Pos[0]) < 1.5) ++nadj;
        }
      } // ipt
      float ntpwc = NumPtsWithCharge(slc, tj, true, tj.EndPt[0], lastPt);
      float nwires = std::abs(tj.Pts[tj.EndPt[0]].Pos[0] - tj.Pts[lastPt].Pos[0]) + 1;
      float hitFrac = ntpwc / nwires;
      if (prt)
        mf::LogVerbatim("TC") << fcnLabel << "-TEP: T" << tj.ID << " lastPt " << lastPt << " npwc "
                              << npwc << " ntpwc " << ntpwc << " nadj " << nadj << " hitFrac "
                              << hitFrac;
      if (hitFrac > fQualityCuts[0] && npwc == minPts && nadj >= minPts - 1) break;
    } // lastPt

    if (prt) mf::LogVerbatim("TC") << " lastPt " << lastPt << " " << tj.EndPt[1] << "\n";
    // trim the last point if it just after a dead wire.
    if (tj.Pts[lastPt].Pos[0] > -0.4) {
      unsigned int prevWire = std::nearbyint(tj.Pts[lastPt].Pos[0]);
      if (tj.StepDir > 0) { --prevWire; }
      else {
        ++prevWire;
      }
      if (prt) {
        mf::LogVerbatim("TC") << fcnLabel << "-TEP: is prevWire " << prevWire << " dead? ";
      }
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      if (prevWire < slc.nWires[plane] && !evt.goodWire[plane][prevWire]) --lastPt;
    } // valid Pos[0]

    // Nothing needs to be done
    if (lastPt == tj.EndPt[1]) {
      if (prt) mf::LogVerbatim("TC") << fcnLabel << "-TEPo: Tj is OK";
      return;
    }

    // clear the points after lastPt
    for (unsigned short ipt = lastPt + 1; ipt <= tj.EndPt[1]; ++ipt)
      UnsetUsedHits(slc, tj.Pts[ipt]);
    SetEndPoints(tj);
    tj.AlgMod[kTEP] = true;
    tj.EndFlag[1][kEndBraggChkd] = false;
    if (prt) {
      fcnLabel += "-TEPo";
      PrintTrajectory(fcnLabel, slc, tj, USHRT_MAX);
    }

  } // TrimEndPts

  /////////////////////////////////////////
  bool
  SignalBetween(const TCSlice& slc,
                const TrajPoint& tp1,
                const TrajPoint& tp2,
                const float& MinWireSignalFraction)
  {
    // Returns true if there is a signal on > MinWireSignalFraction of the wires between tp1 and tp2.
    if (MinWireSignalFraction == 0) return true;

    if (tp1.Pos[0] < -0.4 || tp2.Pos[0] < -0.4) return false;
    int fromWire = std::nearbyint(tp1.Pos[0]);
    int toWire = std::nearbyint(tp2.Pos[0]);

    if (fromWire == toWire) {
      TrajPoint tp = tp1;
      // check for a signal midway between
      tp.Pos[1] = 0.5 * (tp1.Pos[1] + tp2.Pos[1]);
      return SignalAtTp(tp);
    }
    // define a trajectory point located at tp1 that has a direction towards tp2
    TrajPoint tp;
    if (!MakeBareTrajPoint(slc, tp1, tp2, tp)) return true;
    return SignalBetween(slc, tp, toWire, MinWireSignalFraction);
  } // SignalBetween

  /////////////////////////////////////////
  bool
  SignalBetween(const TCSlice& slc, TrajPoint tp, float toPos0, const float& MinWireSignalFraction)
  {
    // Returns true if there is a signal on > MinWireSignalFraction of the wires between tp and toPos0.
    return ChgFracBetween(slc, tp, toPos0) >= MinWireSignalFraction;
  } // SignalBetween

  /////////////////////////////////////////
  float
  ChgFracBetween(const TCSlice& slc, TrajPoint tp, float toPos0)
  {
    // Returns the fraction of wires between tp.Pos[0] and toPos0 that have a hit
    // on the line defined by tp.Pos and tp.Dir

    if (tp.Pos[0] < -0.4 || toPos0 < -0.4) return 0;
    int fromWire = std::nearbyint(tp.Pos[0]);
    int toWire = std::nearbyint(toPos0);

    if (fromWire == toWire) return SignalAtTp(tp);

    int nWires = abs(toWire - fromWire) + 1;

    if (std::abs(tp.Dir[0]) < 0.001) tp.Dir[0] = 0.001;
    float stepSize = std::abs(1 / tp.Dir[0]);
    // ensure that we step in the right direction
    if (toWire > fromWire && tp.Dir[0] < 0) stepSize = -stepSize;
    if (toWire < fromWire && tp.Dir[0] > 0) stepSize = -stepSize;
    float nsig = 0;
    float num = 0;
    for (unsigned short cnt = 0; cnt < nWires; ++cnt) {
      ++num;
      if (SignalAtTp(tp)) ++nsig;
      tp.Pos[0] += tp.Dir[0] * stepSize;
      tp.Pos[1] += tp.Dir[1] * stepSize;
    } // cnt
    float sigFrac = nsig / num;
    return sigFrac;
  } // ChgFracBetween

  ////////////////////////////////////////////////
  bool
  JTHitsOK(TCSlice& slc,
           const std::vector<unsigned int>& iHitsInMultiplet,
           const std::vector<unsigned int>& jHitsInMultiplet)
  {
    // Version of TrajHitsOK that uses hit separation instead of tick overlap
    // on adjacent wires
    for (auto& iht : iHitsInMultiplet) {
      auto const& ihit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      Point2_t ipos;
      ipos[0] = ihit.WireID().Wire;
      ipos[1] = ihit.PeakTime() * tcc.unitsPerTick;
      for (auto& jht : jHitsInMultiplet) {
        auto const& jhit = (*evt.allHits)[slc.slHits[jht].allHitsIndex];
        Point2_t jpos;
        jpos[0] = jhit.WireID().Wire;
        jpos[1] = jhit.PeakTime() * tcc.unitsPerTick;
        if(PosSep2(ipos, jpos) < tcc.JTMaxHitSep2) return true;
      } // jht
    } // iht
    return false;
  } // JTHitsOK

  ////////////////////////////////////////////////
  bool
  TrajHitsOK(TCSlice& slc,
             const std::vector<unsigned int>& iHitsInMultiplet,
             const std::vector<unsigned int>& jHitsInMultiplet)
  {
    // Hits (assume to be on adjacent wires have an acceptable signal overlap

    if (iHitsInMultiplet.empty() || jHitsInMultiplet.empty()) return false;

    float sum;
    float cvI = HitsPosTick(slc, iHitsInMultiplet, sum, kAllHits);
    if (cvI < 0) return false;
    float minI = 1E6;
    float maxI = 0;
    for (auto& iht : iHitsInMultiplet) {
      auto const& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      float cv = hit.PeakTime();
      float rms = hit.RMS();
      float arg = cv - 3.1 * rms;
      if (arg < minI) minI = arg;
      arg = cv + 3.1 * rms;
      if (arg > maxI) maxI = arg;
    }

    float cvJ = HitsPosTick(slc, jHitsInMultiplet, sum, kAllHits);
    if (cvJ < 0) return false;
    float minJ = 1E6;
    float maxJ = 0;
    for (auto& jht : jHitsInMultiplet) {
      auto& hit = (*evt.allHits)[slc.slHits[jht].allHitsIndex];
      float cv = hit.PeakTime();
      float rms = hit.RMS();
      float arg = cv - 3.1 * rms;
      if (arg < minJ) minJ = arg;
      arg = cv + 3.1 * rms;
      if (arg > maxJ) maxJ = arg;
    }

    if (cvI < cvJ) {
      if (maxI > minJ) return true;
    }
    else {
      if (minI < maxJ) return true;
    }
    return false;
  } // TrajHitsOK

  /////////////////////////////////////////
  bool
  TrajHitsOK(TCSlice& slc, const unsigned int iht, const unsigned int jht)
  {
    // ensure that two adjacent hits have an acceptable overlap
    if (iht > slc.slHits.size() - 1) return false;
    if (jht > slc.slHits.size() - 1) return false;
    // require that they be on adjacent wires
    auto& ihit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
    auto& jhit = (*evt.allHits)[slc.slHits[jht].allHitsIndex];
    int iwire = ihit.WireID().Wire;
    int jwire = jhit.WireID().Wire;
    if (std::abs(iwire - jwire) > 1) return false;
    if (ihit.PeakTime() > jhit.PeakTime()) {
      float minISignal = ihit.PeakTime() - 3 * ihit.RMS();
      float maxJSignal = jhit.PeakTime() + 3 * ihit.RMS();
      if (maxJSignal > minISignal) return true;
    }
    else {
      float maxISignal = ihit.PeakTime() + 3 * ihit.RMS();
      float minJSignal = jhit.PeakTime() - 3 * ihit.RMS();
      if (minJSignal > maxISignal) return true;
    }
    return false;
  } // TrajHitsOK

  ////////////////////////////////////////////////
  float
  ExpectedHitsRMS(TCSlice& slc, const TrajPoint& tp)
  {
    // returns the expected RMS of hits for the trajectory point in ticks
    if (std::abs(tp.Dir[0]) > 0.001) {
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      return 1.5 * evt.aveHitRMS[planeID.Plane] +
             2 * std::abs(tp.Dir[1] / tp.Dir[0]) / tcc.unitsPerTick;
    }
    else {
      return 500;
    }
  } // ExpectedHitsRMS

  ////////////////////////////////////////////////
  bool
  SignalAtTpInSlc(const TCSlice& slc, const TrajPoint& tp)
  {
    // Version of SignalAtTP that only checks the hit collection in the current slice

    if (tp.Pos[0] < -0.4) return false;
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short pln = planeID.Plane;
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if (wire > evt.goodWire[pln].size() - 1) return false;
    // assume there is a signal on a dead wire
    if (!evt.goodWire[pln][wire]) return true;
    // no signal here if there are no hits on this wire
    if (slc.wireHitRange[pln][wire].first == UINT_MAX) return false;
    // check the proximity of all of the hits in the range
    float projTick = (float)(tp.Pos[1] / tcc.unitsPerTick);
    float tickRange = 0;
    if (std::abs(tp.Dir[1]) != 0) {
      tickRange = std::abs(0.5 / tp.Dir[1]) / tcc.unitsPerTick;
      // don't let it get too large
      if (tickRange > 40) tickRange = 40;
    }
    float loTpTick = projTick - tickRange;
    float hiTpTick = projTick + tickRange;
    for (unsigned int iht = slc.wireHitRange[pln][wire].first;
         iht <= slc.wireHitRange[pln][wire].second;
         ++iht) {
      unsigned int ahi = slc.slHits[iht].allHitsIndex;
      auto& hit = (*evt.allHits)[ahi];
      if (projTick < hit.PeakTime()) {
        float loHitTick = hit.PeakTime() - 3 * hit.RMS();
        if (hiTpTick > loHitTick) return true;
      }
      else {
        float hiHitTick = hit.PeakTime() + 3 * hit.RMS();
        if (loTpTick < hiHitTick) return true;
      }
    } // iht
    return false;
  } // SignalAtTpInSlc

  /////////////////////////////////////////
  bool
  SignalAtTp(TrajPoint& tp)
  {
    // returns true if there is a hit near tp.Pos by searching through the full hit collection (if there
    // are multiple slices) or through the last slice (if there is only one slice)

    tp.Environment[kEnvNearSrcHit] = false;

    // just check the hits in the last slice
    if (evt.wireHitRange.empty()) {
      const auto& slc = slices[slices.size() - 1];
      return SignalAtTpInSlc(slc, tp);
    }

    if (tp.Pos[0] < -0.4) return false;
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short pln = planeID.Plane;
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if (wire > evt.goodWire[pln].size() - 1) return false;
    // assume there is a signal on a dead wire
    if (!evt.goodWire[pln][wire]) return true;

    // check the proximity of all of the hits in the range
    float projTick = (float)(tp.Pos[1] / tcc.unitsPerTick);
    float tickRange = 0;
    if (std::abs(tp.Dir[1]) != 0) {
      tickRange = std::abs(0.5 / tp.Dir[1]) / tcc.unitsPerTick;
      // don't let it get too large
      if (tickRange > 40) tickRange = 40;
    }
    float loTpTick = projTick - tickRange;
    float hiTpTick = projTick + tickRange;

    // no signal here if there are no hits on this wire
    if (evt.wireHitRange[pln][wire].first == UINT_MAX) return false;

    for (unsigned int iht = evt.wireHitRange[pln][wire].first;
         iht <= evt.wireHitRange[pln][wire].second;
         ++iht) {
      auto& hit = (*evt.allHits)[iht];
      // We wouldn't need to make this check if hits were sorted
      const auto& wid = hit.WireID();
      if (wid.Cryostat != planeID.Cryostat) continue;
      if (wid.TPC != planeID.TPC) continue;
      if (wid.Plane != planeID.Plane) continue;
      if (projTick < hit.PeakTime()) {
        float loHitTick = hit.PeakTime() - 3 * hit.RMS();
        if (hiTpTick > loHitTick) return true;
      }
      else {
        float hiHitTick = hit.PeakTime() + 3 * hit.RMS();
        if (loTpTick < hiHitTick) return true;
      }
    } // iht
    // No hit was found near projTick. Search through the source hits collection
    // (if it is defined) for a hit that may have been removed by disambiguation
    // Use the srcHit collection if it is available
    if (evt.srcHits != NULL) {
      if (NearbySrcHit(planeID, wire, loTpTick, hiTpTick)) {
        tp.Environment[kEnvNearSrcHit] = true;
        return true;
      } // NearbySrcHit
    }   // evt.srcHits != NULL
    return false;
  } // SignalAtTp

  //////////////////////////////////////////
  bool
  NearbySrcHit(geo::PlaneID plnID, unsigned int wire, float loTick, float hiTick)
  {
    // Look for a hit on wid in the srcHits collection that has a tick in the range. This
    // is a DUNE-specific function in which hit disambiguation is done in the U and V planes
    if (evt.srcHits == NULL) return false;
    unsigned int pln = plnID.Plane;
    if (pln == 2) return false;

    unsigned int cstat = plnID.Cryostat;
    unsigned int tpc = plnID.TPC;
    // get a valid range of hits to search
    if (evt.tpcSrcHitRange[tpc].first >= (*evt.srcHits).size()) return false;
    if (evt.tpcSrcHitRange[tpc].second >= (*evt.srcHits).size()) return false;
    raw::ChannelID_t chan = tcc.geom->PlaneWireToChannel((int)pln, (int)wire, (int)tpc, (int)cstat);
    float atTick = 0.5 * (loTick + hiTick);
    for (unsigned int iht = evt.tpcSrcHitRange[tpc].first; iht <= evt.tpcSrcHitRange[tpc].second;
         ++iht) {
      auto& hit = (*evt.srcHits)[iht];
      if (hit.Channel() != chan) continue;
      if (atTick < hit.PeakTime()) {
        float loHitTick = hit.PeakTime() - 3 * hit.RMS();
        if (hiTick > loHitTick) return true;
      }
      else {
        float hiHitTick = hit.PeakTime() + 3 * hit.RMS();
        if (loTick < hiHitTick) return true;
      }
    } // iht
    return false;
  } // NearbySrcHit

  //////////////////////////////////////////
  float
  TpSumHitChg(const TCSlice& slc, TrajPoint const& tp)
  {
    float totchg = 0;
    for (size_t i = 0; i < tp.Hits.size(); ++i) {
      if (!tp.UseHit[i]) continue;
      totchg += (*evt.allHits)[slc.slHits[tp.Hits[i]].allHitsIndex].Integral();
    }
    return totchg;
  } // TpSumHitChg

  //////////////////////////////////////////
  unsigned short
  NumPtsWithCharge(const TCSlice& slc, const Trajectory& tj, bool includeDeadWires)
  {
    unsigned short firstPt = tj.EndPt[0];
    unsigned short lastPt = tj.EndPt[1];
    return NumPtsWithCharge(slc, tj, includeDeadWires, firstPt, lastPt);
  }

  //////////////////////////////////////////
  unsigned short
  NumPtsWithCharge(const TCSlice& slc,
                   const Trajectory& tj,
                   bool includeDeadWires,
                   unsigned short firstPt,
                   unsigned short lastPt)
  {
    unsigned short ntp = 0;
    for (unsigned short ipt = firstPt; ipt <= lastPt; ++ipt)
      if (tj.Pts[ipt].Chg > 0) ++ntp;
    // Add the count of deadwires
    if (includeDeadWires) ntp += DeadWireCount(slc, tj.Pts[firstPt], tj.Pts[lastPt]);
    return ntp;
  } // NumPtsWithCharge

  //////////////////////////////////////////
  float
  DeadWireCount(const TCSlice& slc, const TrajPoint& tp1, const TrajPoint& tp2)
  {
    return DeadWireCount(slc, tp1.Pos[0], tp2.Pos[0], tp1.CTP);
  } // DeadWireCount

  //////////////////////////////////////////
  float
  DeadWireCount(const TCSlice& slc, const float& inWirePos1, const float& inWirePos2, CTP_t tCTP)
  {
    if (inWirePos1 < -0.4 || inWirePos2 < -0.4) return 0;
    unsigned int inWire1 = std::nearbyint(inWirePos1);
    unsigned int inWire2 = std::nearbyint(inWirePos2);
    geo::PlaneID planeID = DecodeCTP(tCTP);
    unsigned short plane = planeID.Plane;
    if (inWire1 > slc.nWires[plane] || inWire2 > slc.nWires[plane]) return 0;
    if (inWire1 > inWire2) {
      // put in increasing order
      unsigned int tmp = inWire1;
      inWire1 = inWire2;
      inWire2 = tmp;
    } // inWire1 > inWire2
    ++inWire2;
    unsigned int wire, ndead = 0;
    for (wire = inWire1; wire < inWire2; ++wire)
      if (!evt.goodWire[plane][wire]) ++ndead;
    return ndead;
  } // DeadWireCount

  ////////////////////////////////////////////////
  unsigned short
  PDGCodeIndex(int PDGCode)
  {
    unsigned short pdg = abs(PDGCode);
    if (pdg == 11) return 0;   // electron
    if (pdg == 13) return 1;   // muon
    if (pdg == 211) return 2;  // pion
    if (pdg == 321) return 3;  // kaon
    if (pdg == 2212) return 4; // proton
    return USHRT_MAX;
  } // PDGCodeIndex

  ////////////////////////////////////////////////
  void
  MakeTrajectoryObsolete(TCSlice& slc, unsigned int itj)
  {
    // Note that this does not change the state of UseHit to allow
    // resurrecting the trajectory later (RestoreObsoleteTrajectory)
    if (itj > slc.tjs.size() - 1) return;
    int killTjID = slc.tjs[itj].ID;
    for (auto& hit : slc.slHits)
      if (hit.InTraj == killTjID) hit.InTraj = 0;
    slc.tjs[itj].AlgMod[kKilled] = true;
  } // MakeTrajectoryObsolete

  ////////////////////////////////////////////////
  void
  RestoreObsoleteTrajectory(TCSlice& slc, unsigned int itj)
  {
    if (itj > slc.tjs.size() - 1) return;
    if (!slc.tjs[itj].AlgMod[kKilled]) {
      mf::LogWarning("TC")
        << "RestoreObsoleteTrajectory: Trying to restore not-obsolete trajectory "
        << slc.tjs[itj].ID;
      return;
    }
    unsigned int iht;
    for (auto& tp : slc.tjs[itj].Pts) {
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if (tp.UseHit[ii]) {
          iht = tp.Hits[ii];
          if (slc.slHits[iht].InTraj == 0) { slc.slHits[iht].InTraj = slc.tjs[itj].ID; }
        }
      } // ii
    }   // tp
    slc.tjs[itj].AlgMod[kKilled] = false;
  } // RestoreObsoleteTrajectory

  //////////////////////////////////////////
  bool
  SplitTraj(detinfo::DetectorPropertiesData const& detProp,
            TCSlice& slc,
            unsigned short itj,
            float XPos,
            bool makeVx2,
            bool prt)
  {
    // Splits the trajectory at an X position and optionally creates a 2D vertex
    // at the split point
    if (itj > slc.tjs.size() - 1) return false;

    auto& tj = slc.tjs[itj];
    geo::PlaneID planeID = DecodeCTP(tj.CTP);
    float atPos1 = detProp.ConvertXToTicks(XPos, planeID) * tcc.unitsPerTick;
    unsigned short atPt = USHRT_MAX;
    for (unsigned short ipt = tj.EndPt[0] + 1; ipt <= tj.EndPt[1]; ++ipt) {
      if (tj.Pts[ipt].Pos[1] > tj.Pts[ipt - 1].Pos[1]) {
        // positive slope
        if (tj.Pts[ipt - 1].Pos[1] < atPos1 && tj.Pts[ipt].Pos[1] >= atPos1) {
          atPt = ipt;
          break;
        }
      }
      else {
        // negative slope
        if (tj.Pts[ipt - 1].Pos[1] >= atPos1 && tj.Pts[ipt].Pos[1] < atPos1) {
          atPt = ipt;
          break;
        }
      } // negative slope
    }   // ipt
    if (atPt == USHRT_MAX) return false;
    unsigned short vx2Index = USHRT_MAX;
    if (makeVx2) {
      VtxStore newVx2;
      newVx2.CTP = tj.CTP;
      newVx2.Pos[0] = 0.5 * (tj.Pts[atPt - 1].Pos[0] + tj.Pts[atPt].Pos[0]);
      newVx2.Pos[1] = 0.5 * (tj.Pts[atPt - 1].Pos[1] + tj.Pts[atPt].Pos[1]);
      newVx2.Topo = 10;
      newVx2.NTraj = 2;
      if (StoreVertex(slc, newVx2)) vx2Index = slc.vtxs.size() - 1;
    } // makeVx2
    return SplitTraj(slc, itj, atPt, vx2Index, prt);
  } // SplitTraj

  //////////////////////////////////////////
  bool
  SplitTraj(TCSlice& slc, unsigned short itj, unsigned short atPt, bool prt)
  {
    // another version that splits the trajectory and makes a new 2D vertex
    if (itj > slc.tjs.size() - 1) return false;
    Trajectory& tj = slc.tjs[itj];
    if (atPt <= tj.EndPt[0] || atPt >= tj.EndPt[1]) return false;

    VtxStore new2V;
    new2V.CTP = tj.CTP;
    new2V.Pos[0] = 0.5 * (tj.Pts[atPt - 1].Pos[0] + tj.Pts[atPt].Pos[0]);
    new2V.Pos[1] = 0.5 * (tj.Pts[atPt - 1].Pos[1] + tj.Pts[atPt].Pos[1]);
    new2V.Topo = 10;
    new2V.NTraj = 2;
    if(!StoreVertex(slc, new2V)) return false;
    unsigned short ivx = slc.vtxs.size() - 1;
    return SplitTraj(slc, itj, atPt, ivx, prt);

  } // SplitTraj

  //////////////////////////////////////////
  bool
  SplitTraj(TCSlice& slc, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt)
  {
    // Splits the trajectory itj in the slc.tjs vector into two trajectories at position pos. Splits
    // the trajectory and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)

    if (itj > slc.tjs.size() - 1) return false;
    if (pos < slc.tjs[itj].EndPt[0] + 1 || pos > slc.tjs[itj].EndPt[1] - 1) return false;
    if (ivx != USHRT_MAX && ivx > slc.vtxs.size() - 1) return false;

    Trajectory& tj = slc.tjs[itj];

    // Reset the PDG Code if we are splitting a tagged muon
    bool splittingMuon = (tj.PDGCode == 13);
    if (splittingMuon) tj.PDGCode = 0;

    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "SplitTraj: Split T" << tj.ID << " at point " << PrintPos(slc, tj.Pts[pos]);
      if (ivx < slc.vtxs.size()) myprt << " with Vtx 2V" << slc.vtxs[ivx].ID;
    }

    // ensure that there will be at least 3 TPs on each trajectory
    unsigned short ntp = 0;
    for (unsigned short ipt = 0; ipt <= pos; ++ipt) {
      if (tj.Pts[ipt].Chg > 0) ++ntp;
      if (ntp > 2) break;
    } // ipt
    if (ntp < 3) {
      if (prt) mf::LogVerbatim("TC") << " Split point to small at begin " << ntp << " pos " << pos;
      return false;
    }
    ntp = 0;
    for (unsigned short ipt = pos + 1; ipt <= tj.EndPt[1]; ++ipt) {
      if (tj.Pts[ipt].Chg > 0) ++ntp;
      if (ntp > 2) break;
    } // ipt
    if (ntp < 3) {
      if (prt)
        mf::LogVerbatim("TC") << " Split point too small at end " << ntp << " pos " << pos
                              << " EndPt " << tj.EndPt[1];
      return false;
    }

    // make a copy that will become the Tj after the split point
    Trajectory newTj = tj;
    newTj.ID = slc.tjs.size() + 1;
    ++evt.globalT_UID;
    newTj.UID = evt.globalT_UID;
    // make another copy in case something goes wrong
    Trajectory oldTj = tj;

    // Leave the first section of tj in place. Re-assign the hits
    // to the new trajectory
    unsigned int iht;
    for (unsigned short ipt = pos + 1; ipt <= tj.EndPt[1]; ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for (unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if (!tj.Pts[ipt].UseHit[ii]) continue;
        iht = tj.Pts[ipt].Hits[ii];
        // This shouldn't happen but check anyway
        if (slc.slHits[iht].InTraj != tj.ID) continue;
        slc.slHits[iht].InTraj = newTj.ID;
        tj.Pts[ipt].UseHit[ii] = false;
      } // ii
    }   // ipt
    SetEndPoints(tj);
    // Update MCSMom and charge properties
    tj.MCSMom = MCSMom(slc, tj);
    UpdateTjChgProperties("ST", slc, tj, prt);
    if (splittingMuon) SetPDGCode(slc, tj);

    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if (eraseSize > newTj.Pts.size() - 1) {
      tj = oldTj;
      return false;
    }

    if (ivx < slc.vtxs.size()) tj.VtxID[1] = slc.vtxs[ivx].ID;
    tj.AlgMod[kSplit] = true;
    if (prt) {
      mf::LogVerbatim("TC") << " Splitting T" << tj.ID << " new EndPts " << tj.EndPt[0] << " to "
                            << tj.EndPt[1];
    }

    // erase the TPs at the beginning of the new trajectory
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for (unsigned short ipt = 0; ipt < 3; ++ipt) {
      for (unsigned short ii = 0; ii < newTj.Pts[ipt].Hits.size(); ++ii)
        newTj.Pts[ipt].UseHit[ii] = false;
      newTj.Pts[ipt].Chg = 0;
    } // ipt
    SetEndPoints(newTj);
    newTj.MCSMom = MCSMom(slc, newTj);
    UpdateTjChgProperties("ST", slc, newTj, prt);
    if (splittingMuon) SetPDGCode(slc, newTj);
    if (ivx < slc.vtxs.size()) newTj.VtxID[0] = slc.vtxs[ivx].ID;
    newTj.AlgMod[kSplit] = true;
    newTj.ParentID = 0;
    slc.tjs.push_back(newTj);

    if (prt) {
      mf::LogVerbatim("TC") << "  newTj T" << newTj.ID << " EndPts " << newTj.EndPt[0] << " to "
                            << newTj.EndPt[1];
    }
    return true;

  } // SplitTraj

  //////////////////////////////////////////
  void
  TrajPointTrajDOCA(const TCSlice& slc,
                    TrajPoint const& tp,
                    Trajectory const& tj,
                    unsigned short& closePt,
                    float& minSep)
  {
    // Finds the point, ipt, on trajectory tj that is closest to trajpoint tp
    float best = minSep * minSep;
    closePt = USHRT_MAX;
    float dw, dt, dp2;
    unsigned short ipt;
    for (ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      dw = tj.Pts[ipt].Pos[0] - tp.Pos[0];
      dt = tj.Pts[ipt].Pos[1] - tp.Pos[1];
      dp2 = dw * dw + dt * dt;
      if (dp2 < best) {
        best = dp2;
        closePt = ipt;
      }
    } // ipt
    minSep = sqrt(best);
  } // TrajPointTrajDOCA

  //////////////////////////////////////////
  bool
  TrajTrajDOCA(const TCSlice& slc,
               const Trajectory& tj1,
               const Trajectory& tj2,
               unsigned short& ipt1,
               unsigned short& ipt2,
               float& minSep)
  {
    return TrajTrajDOCA(slc, tj1, tj2, ipt1, ipt2, minSep, false);
  } // TrajTrajDOCA

  //////////////////////////////////////////
  bool
  TrajTrajDOCA(const TCSlice& slc,
               const Trajectory& tj1,
               const Trajectory& tj2,
               unsigned short& ipt1,
               unsigned short& ipt2,
               float& minSep,
               bool considerDeadWires)
  {
    // Find the Distance Of Closest Approach between two trajectories less than minSep
    // start with some rough cuts to minimize the use of the more expensive checking. This
    // function returns true if the DOCA is less than minSep
    for (unsigned short iwt = 0; iwt < 2; ++iwt) {
      // Apply box cuts on the ends of the trajectories
      // The Lo/Hi wire(time) at each end of tj1
      float wt0 = tj1.Pts[tj1.EndPt[0]].HitPos[iwt];
      float wt1 = tj1.Pts[tj1.EndPt[1]].HitPos[iwt];
      float lowt1 = wt0;
      float hiwt1 = wt1;
      if (wt1 < lowt1) {
        lowt1 = wt1;
        hiwt1 = wt0;
      }
      // The Lo/Hi wire(time) at each end of tj2
      wt0 = tj2.Pts[tj2.EndPt[0]].HitPos[iwt];
      wt1 = tj2.Pts[tj2.EndPt[1]].HitPos[iwt];
      float lowt2 = wt0;
      float hiwt2 = wt1;
      if (wt1 < lowt2) {
        lowt2 = wt1;
        hiwt2 = wt0;
      }
      // Check for this configuration
      //  loWire1.......hiWire1   minSep  loWire2....hiWire2
      //  loTime1.......hiTime1   minSep  loTime2....hiTime2
      if (lowt2 > hiwt1 + minSep) return false;
      // and the other
      if (lowt1 > hiwt2 + minSep) return false;
    } // iwt

    float best = minSep * minSep;
    ipt1 = 0;
    ipt2 = 0;
    float dwc = 0;
    bool isClose = false;
    for (unsigned short i1 = tj1.EndPt[0]; i1 < tj1.EndPt[1] + 1; ++i1) {
      for (unsigned short i2 = tj2.EndPt[0]; i2 < tj2.EndPt[1] + 1; ++i2) {
        if (considerDeadWires) dwc = DeadWireCount(slc, tj1.Pts[i1], tj2.Pts[i2]);
        float dw = tj1.Pts[i1].HitPos[0] - tj2.Pts[i2].HitPos[0] - dwc;
        if (std::abs(dw) > minSep) continue;
        float dt = tj1.Pts[i1].HitPos[1] - tj2.Pts[i2].HitPos[1];
        if (std::abs(dt) > minSep) continue;
        float dp2 = dw * dw + dt * dt;
        if (dp2 < best) {
          best = dp2;
          ipt1 = i1;
          ipt2 = i2;
          isClose = true;
        }
      } // i2
    }   // i1
    minSep = sqrt(best);
    return isClose;
  } // TrajTrajDOCA

  //////////////////////////////////////////
  float
  HitSep2(const TCSlice& slc, unsigned int iht, unsigned int jht)
  {
    // returns the separation^2 between two hits in WSE units
    if (iht > slc.slHits.size() - 1 || jht > slc.slHits.size() - 1) return 1E6;
    auto& ihit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
    auto& jhit = (*evt.allHits)[slc.slHits[jht].allHitsIndex];
    float dw = (float)ihit.WireID().Wire - (float)jhit.WireID().Wire;
    float dt = (ihit.PeakTime() - jhit.PeakTime()) * tcc.unitsPerTick;
    return dw * dw + dt * dt;
  } // HitSep2

  //////////////////////////////////////////
  unsigned short
  CloseEnd(const TCSlice& slc, const Trajectory& tj, const Point2_t& pos)
  {
    unsigned short endPt = tj.EndPt[0];
    auto& tp0 = tj.Pts[endPt];
    endPt = tj.EndPt[1];
    auto& tp1 = tj.Pts[endPt];
    if (PosSep2(tp0.Pos, pos) < PosSep2(tp1.Pos, pos)) return 0;
    return 1;
  } // CloseEnd

  //////////////////////////////////////////
  float
  PointTrajSep2(float wire, float time, TrajPoint const& tp)
  {
    float dw = wire - tp.Pos[0];
    float dt = time - tp.Pos[1];
    return dw * dw + dt * dt;
  }

  //////////////////////////////////////////
  float
  PointTrajDOCA(const TCSlice& slc, unsigned int iht, TrajPoint const& tp)
  {
    if (iht > slc.slHits.size() - 1) return 1E6;
    auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
    float wire = hit.WireID().Wire;
    float time = hit.PeakTime() * tcc.unitsPerTick;
    return sqrt(PointTrajDOCA2(slc, wire, time, tp));
  } // PointTrajDOCA

  //////////////////////////////////////////
  float
  PointTrajDOCA(const TCSlice& slc, float wire, float time, TrajPoint const& tp)
  {
    return sqrt(PointTrajDOCA2(slc, wire, time, tp));
  } // PointTrajDOCA

  //////////////////////////////////////////
  float
  PointTrajDOCA2(const TCSlice& slc, float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach squared between a (wire, time(WSE)) point
    // and a trajectory point

    double t = (double)(wire - tp.Pos[0]) * tp.Dir[0] + (double)(time - tp.Pos[1]) * tp.Dir[1];
    double dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    double dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return (float)(dw * dw + dt * dt);

  } // PointTrajDOCA2

  //////////////////////////////////////////
  void
  TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, Point2_t& pos)
  {
    TrajIntersection(tp1, tp2, pos[0], pos[1]);
  } // TrajIntersection
  //////////////////////////////////////////
  void
  TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y)
  {
    // returns the intersection position, (x,y), of two trajectory points

    x = -9999;
    y = -9999;

    double arg1 = tp1.Pos[0] * tp1.Dir[1] - tp1.Pos[1] * tp1.Dir[0];
    double arg2 = tp2.Pos[0] * tp1.Dir[1] - tp2.Pos[1] * tp1.Dir[0];
    double arg3 = tp2.Dir[0] * tp1.Dir[1] - tp2.Dir[1] * tp1.Dir[0];
    if (arg3 == 0) return;
    double s = (arg1 - arg2) / arg3;

    x = (float)(tp2.Pos[0] + s * tp2.Dir[0]);
    y = (float)(tp2.Pos[1] + s * tp2.Dir[1]);

  } // TrajIntersection

  //////////////////////////////////////////
  float
  MaxTjLen(const TCSlice& slc, std::vector<int>& tjIDs)
  {
    // returns the length of the longest Tj in the supplied list
    if (tjIDs.empty()) return 0;
    float maxLen = 0;
    for (auto tjid : tjIDs) {
      if (tjid < 1 || tjid > (int)slc.tjs.size()) continue;
      auto& tj = slc.tjs[tjid - 1];
      float sep2 = PosSep2(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
      if (sep2 > maxLen) maxLen = sep2;
    } // tj
    return sqrt(maxLen);
  } // MaxTjLen

  //////////////////////////////////////////
  float
  TrajLength(const Trajectory& tj)
  {
    float len = 0, dx, dy;
    unsigned short ipt;
    unsigned short prevPt = tj.EndPt[0];
    for (ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1] + 1; ++ipt) {
      if (tj.Pts[ipt].Chg == 0) continue;
      dx = tj.Pts[ipt].Pos[0] - tj.Pts[prevPt].Pos[0];
      dy = tj.Pts[ipt].Pos[1] - tj.Pts[prevPt].Pos[1];
      len += sqrt(dx * dx + dy * dy);
      prevPt = ipt;
    }
    return len;
  } // TrajLength

  //////////////////////////////////////////
  float
  PosSep(const Point2_t& pos1, const Point2_t& pos2)
  {
    return sqrt(PosSep2(pos1, pos2));
  } // PosSep

  //////////////////////////////////////////
  float
  PosSep2(const Point2_t& pos1, const Point2_t& pos2)
  {
    // returns the separation distance^2 between two positions
    float d0 = pos1[0] - pos2[0];
    float d1 = pos1[1] - pos2[1];
    return d0 * d0 + d1 * d1;
  } // PosSep2

  //////////////////////////////////////////
  float
  TrajPointSeparation(const TrajPoint& tp1, const TrajPoint& tp2)
  {
    // Returns the separation distance between two trajectory points
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation

  //////////////////////////////////////////
  bool
  TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& DOCA)
  {
    // find the closest approach between a trajectory tj and a point (x,y). Returns
    // the index of the closest trajectory point and the distance. Returns false if none
    // of the points on the tj are within DOCA

    float close2 = DOCA * DOCA;
    closePt = 0;
    bool foundClose = false;

    for (unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      if (tj.Pts[ipt].Chg == 0) continue;
      float dx = tj.Pts[ipt].Pos[0] - x;
      if (std::abs(dx) > DOCA) continue;
      float dy = tj.Pts[ipt].Pos[1] - y;
      if (std::abs(dy) > DOCA) continue;
      float sep2 = dx * dx + dy * dy;
      if (sep2 < close2) {
        close2 = sep2;
        closePt = ipt;
        foundClose = true;
      }
    } // ipt

    DOCA = sqrt(close2);
    return foundClose;

  } // TrajClosestApproach

  /////////////////////////////////////////
  float
  TwoTPAngle(const TrajPoint& tp1, const TrajPoint& tp2)
  {
    // Calculates the angle of a line between two TPs
    float dw = tp2.Pos[0] - tp1.Pos[0];
    float dt = tp2.Pos[1] - tp1.Pos[1];
    return atan2(dw, dt);
  } // TwoTPAngle

  ////////////////////////////////////////////////
  std::vector<unsigned int>
  PutHitsInVector(const TCSlice& slc, PFPStruct const& pfp, HitStatus_t hitRequest)
  {
    // Put hits with the assn P -> TP3D -> TP -> Hit into a vector
    std::vector<unsigned int> hitVec;
    if (pfp.TP3Ds.empty()) return hitVec;

    for(auto& tp3d : pfp.TP3Ds) {
      // ignore end points
      if(tp3d.TPIndex == USHRT_MAX) continue;
      if(tp3d.TjID <= 0) continue;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        bool useit = (hitRequest == kAllHits);
        if (tp.UseHit[ii] && hitRequest == kUsedHits) useit = true;
        if (!tp.UseHit[ii] && hitRequest == kUnusedHits) useit = true;
        if (useit) hitVec.push_back(iht);
      }
    } // tp3d
    return hitVec;
  } // PutHitsInVector

  ////////////////////////////////////////////////
  std::vector<unsigned int>
  PutTrajHitsInVector(const Trajectory& tj, HitStatus_t hitRequest)
  {
    // Put hits (which are indexed into slHits) in each trajectory point into a flat vector
    std::vector<unsigned int> hitVec;

    // reserve under the assumption that there will be one hit per point
    hitVec.reserve(tj.Pts.size());
    for (unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      for (unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        unsigned int iht = tj.Pts[ipt].Hits[ii];
        bool useit = (hitRequest == kAllHits);
        if (tj.Pts[ipt].UseHit[ii] && hitRequest == kUsedHits) useit = true;
        if (!tj.Pts[ipt].UseHit[ii] && hitRequest == kUnusedHits) useit = true;
        if (useit) hitVec.push_back(iht);
      } // iht
    }   // ipt
    return hitVec;
  } // PutTrajHitsInVector

  //////////////////////////////////////////
  void
  TagJunkTj(TCSlice& slc, Trajectory& tj, bool prt)
  {
    // Characterizes the trajectory as a junk tj even though it may not
    // have been reconstructed in FindJunkTraj. The distinguishing feature is
    // that it is short and has many used hits in each trajectory point.

    // Don't bother if it is too long
    if (tj.Pts.size() > 10) return;
    if (tj.PDGCode == 111) return;
    // count the number of points that have many used hits
    unsigned short nhm = 0;
    unsigned short npwc = 0;
    for (auto& tp : tj.Pts) {
      if (tp.Chg == 0) continue;
      ++npwc;
      unsigned short nused = 0;
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if (tp.UseHit[ii]) ++nused;
      } // ii
      if (nused > 3) ++nhm;
    } // tp
    // Set the junkTj bit if most of the hits are used in most of the tps
    if (nhm > 0.5 * npwc) tj.AlgMod[kJunkTj] = true;
    if (prt)
      mf::LogVerbatim("TC") << "TGT: T" << tj.ID << " npwc " << npwc << " nhm " << nhm << " junk? "
                            << tj.AlgMod[kJunkTj];
  } // TagJunkTj

  //////////////////////////////////////////
  bool
  HasDuplicateHits(const TCSlice& slc, Trajectory const& tj, bool prt)
  {
    // returns true if a hit is associated with more than one TP
    auto tjHits = PutTrajHitsInVector(tj, kAllHits);
    for (unsigned short ii = 0; ii < tjHits.size() - 1; ++ii) {
      for (unsigned short jj = ii + 1; jj < tjHits.size(); ++jj) {
        if (tjHits[ii] == tjHits[jj]) {
          if (prt)
            mf::LogVerbatim("TC") << "HDH: Hit " << PrintHit(slc.slHits[ii]) << " is a duplicate "
                                  << ii << " " << jj;
          return true;
        }
      } // jj
    }   // ii
    return false;
  } // HasDuplicateHits

  //////////////////////////////////////////
  void
  MoveTPToWire(TrajPoint& tp, float wire)
  {
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    if (tp.Dir[0] == 0) return;
    float dw = wire - tp.Pos[0];
    if (std::abs(dw) < 0.01) return;
    tp.Pos[0] = wire;
    tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
  } // MoveTPToWire

  //////////////////////////////////////////
  std::vector<unsigned int>
  FindCloseHits(const TCSlice& slc,
                std::array<int, 2> const& wireWindow,
                Point2_t const& timeWindow,
                const unsigned short plane,
                HitStatus_t hitRequest,
                bool usePeakTime,
                bool& hitsNear)
  {
    // returns a vector of hits that are within the Window[Pos0][Pos1] in plane.
    // Note that hits on wire wireWindow[1] are returned as well. The definition of close
    // depends on setting of usePeakTime. If UsePeakTime is true, a hit is considered nearby if
    // the PeakTime is within the window. This is shown schematically here where
    // the time is on the horizontal axis and a "-" denotes a valid entry
    // timeWindow     -----------------
    // hit PeakTime             +         close
    // hit PeakTime  +                    not close
    // If usePeakTime is false, a hit is considered nearby if the hit StartTick and EndTick overlap with the timeWindow
    // Time window                  ---------
    // Hit StartTick-EndTick      --------        close
    // Hit StartTick - EndTick                  --------  not close

    hitsNear = false;
    std::vector<unsigned int> closeHits;
    if (plane > slc.firstWire.size() - 1) return closeHits;
    // window in the wire coordinate
    int loWire = wireWindow[0];
    if (loWire < (int)slc.firstWire[plane]) loWire = slc.firstWire[plane];
    int hiWire = wireWindow[1];
    if (hiWire > (int)slc.lastWire[plane] - 1) hiWire = slc.lastWire[plane] - 1;
    // window in the time coordinate
    float minTick = timeWindow[0] / tcc.unitsPerTick;
    float maxTick = timeWindow[1] / tcc.unitsPerTick;
    for (int wire = loWire; wire <= hiWire; ++wire) {
      // Set hitsNear if the wire is dead
      if (!evt.goodWire[plane][wire]) hitsNear = true;
      if (slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
      unsigned int firstHit = slc.wireHitRange[plane][wire].first;
      unsigned int lastHit = slc.wireHitRange[plane][wire].second;
      for (unsigned int iht = firstHit; iht <= lastHit; ++iht) {
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        if (usePeakTime) {
          if (hit.PeakTime() < minTick) continue;
          if (hit.PeakTime() > maxTick) break;
        }
        else {
          int hiLo = minTick;
          if (hit.StartTick() > hiLo) hiLo = hit.StartTick();
          int loHi = maxTick;
          if (hit.EndTick() < loHi) loHi = hit.EndTick();
          if (loHi < hiLo) continue;
          if (hiLo > loHi) break;
        }
        hitsNear = true;
        bool takeit = (hitRequest == kAllHits);
        if (hitRequest == kUsedHits && slc.slHits[iht].InTraj > 0) takeit = true;
        if (hitRequest == kUnusedHits && slc.slHits[iht].InTraj == 0) takeit = true;
        if (takeit) closeHits.push_back(iht);
      } // iht
    }   // wire
    return closeHits;
  } // FindCloseHits

  //////////////////////////////////////////
  bool
  FindCloseHits(TCSlice& slc, TrajPoint& tp, float const& maxDelta, HitStatus_t hitRequest)
  {
    // Fills tp.Hits sets tp.UseHit true for hits that are close to tp.Pos. Returns true if there are
    // close hits OR if the wire at this position is dead

    tp.Hits.clear();
    tp.UseHit.reset();
    tp.Environment.reset();
    if (!WireHitRangeOK(slc, tp.CTP)) { return false; }

    if (tp.Pos[0] < -0.4) return false;
    unsigned short plane = DecodeCTP(tp.CTP).Plane;
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if (wire < slc.firstWire[plane]) return false;
    if (wire > slc.lastWire[plane] - 1) return false;

    // dead wire
    if (!evt.goodWire[plane][wire]) {
      tp.Environment[kEnvNotGoodWire] = true;
      return true;
    }
    tp.Environment[kEnvNotGoodWire] = false;
    // live wire with no hits
    if (slc.wireHitRange[plane][wire].first == UINT_MAX) return false;

    unsigned int firstHit = slc.wireHitRange[plane][wire].first;
    unsigned int lastHit = slc.wireHitRange[plane][wire].second;

    float fwire = wire;
    for (unsigned int iht = firstHit; iht <= lastHit; ++iht) {
      if ((unsigned int)slc.slHits[iht].InTraj > slc.tjs.size()) continue;
      bool useit = (hitRequest == kAllHits);
      if (hitRequest == kUsedHits && slc.slHits[iht].InTraj > 0) useit = true;
      if (hitRequest == kUnusedHits && slc.slHits[iht].InTraj == 0) useit = true;
      if (!useit) continue;
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      float ftime = tcc.unitsPerTick * hit.PeakTime();
      float delta = PointTrajDOCA(slc, fwire, ftime, tp);
      if (delta < maxDelta) tp.Hits.push_back(iht);
    } // iht
    if (tp.Hits.size() > 16) { tp.Hits.resize(16); }
    // Set UseHit false. The calling routine should decide if these hits should be used
    tp.UseHit.reset();
    return (!tp.Hits.empty());

  } // FindCloseHits

  //////////////////////////////////////////
  unsigned short
  NearbyCleanPt(const TCSlice& slc, const Trajectory& tj, unsigned short end)
  {
    // Searches for a TP near the end (or beginnin) that doesn't have the kEnvOverlap bit set
    // with the intent that a fit of a vertex position using this tj will be minimally
    // biased if there are no nearby hits from other tjs. A search is done from the
    // supplied nearPt moving in the + direction if nearPt == tj.EndPt[0] and moving in
    // the - direction if nearPt == tj.EndPt[1]
    if (end > 1) return USHRT_MAX;
    short dir = 1;
    if (end == 1) dir = -1;
    for (short ii = 0; ii < (short)tj.Pts.size(); ++ii) {
      short ipt = tj.EndPt[end] + dir * ii;
      if (ipt < 0 || ipt >= (short)tj.Pts.size()) return USHRT_MAX;
      auto& tp = tj.Pts[ipt];
      if (!tp.Environment[kEnvOverlap]) return ipt;
    } // ii
    return tj.EndPt[end];
  } // FindCleanPt

  //////////////////////////////////////////
  std::vector<int>
  FindCloseTjs(const TCSlice& slc,
               const TrajPoint& fromTp,
               const TrajPoint& toTp,
               const float& maxDelta)
  {
    // Returns a list of Tj IDs that have hits within distance maxDelta on a line drawn between the two Tps as shown
    // graphically here, where a "*" is a Tp and "|" and "-" are the boundaries of the region that is checked
    //
    //    ---------------
    //    |             |
    //    *             *
    //    |             |
    //    ---------------
    // If the wire positions of fromTp and toTp are the same, a different region is checked as shown here
    //
    //     -----------
    //     |         |
    //     |    *    |
    //     |         |
    //     -----------

    std::vector<int> tmp;
    if (fromTp.Pos[0] < -0.4 || toTp.Pos[0] < -0.4) return tmp;

    TrajPoint tp;
    // Make the tp so that stepping is positive
    unsigned int firstWire, lastWire;
    if (toTp.Pos[0] > fromTp.Pos[0]) {
      if (!MakeBareTrajPoint(slc, fromTp, toTp, tp)) return tmp;
      firstWire = std::nearbyint(fromTp.Pos[0]);
      lastWire = std::nearbyint(toTp.Pos[0]);
    }
    else if (toTp.Pos[0] < fromTp.Pos[0]) {
      if (!MakeBareTrajPoint(slc, toTp, fromTp, tp)) return tmp;
      firstWire = std::nearbyint(toTp.Pos[0]);
      lastWire = std::nearbyint(fromTp.Pos[0]);
    }
    else {
      tp.Pos = fromTp.Pos;
      float tmp = fromTp.Pos[0] - maxDelta;
      if (tmp < 0) tmp = 0;
      firstWire = std::nearbyint(tmp);
      tmp = fromTp.Pos[0] + maxDelta;
      lastWire = std::nearbyint(tmp);
    }

    unsigned short plane = DecodeCTP(tp.CTP).Plane;

    if (firstWire < slc.firstWire[plane]) firstWire = slc.firstWire[plane];
    if (firstWire > slc.lastWire[plane] - 1) return tmp;
    if (lastWire < slc.firstWire[plane]) return tmp;
    if (lastWire > slc.lastWire[plane] - 1) lastWire = slc.lastWire[plane] - 1;

    for (unsigned int wire = firstWire; wire <= lastWire; ++wire) {
      if (slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
      MoveTPToWire(tp, (float)wire);
      // Find the tick range at this position
      float minTick = (tp.Pos[1] - maxDelta) / tcc.unitsPerTick;
      float maxTick = (tp.Pos[1] + maxDelta) / tcc.unitsPerTick;
      unsigned int firstHit = slc.wireHitRange[plane][wire].first;
      unsigned int lastHit = slc.wireHitRange[plane][wire].second;
      for (unsigned int iht = firstHit; iht <= lastHit; ++iht) {
        if (slc.slHits[iht].InTraj <= 0) continue;
        if ((unsigned int)slc.slHits[iht].InTraj > slc.tjs.size()) continue;
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        if (hit.PeakTime() < minTick) continue;
        // Hits are sorted by increasing time so we can break when maxTick is reached
        if (hit.PeakTime() > maxTick) break;
        if (std::find(tmp.begin(), tmp.end(), slc.slHits[iht].InTraj) != tmp.end()) continue;
        tmp.push_back(slc.slHits[iht].InTraj);
      } // iht
    }   // wire

    return tmp;

  } // FindCloseTjs

  ////////////////////////////////////////////////
  float
  KinkSignificance(TCSlice& slc,
                   Trajectory& tj1,
                   unsigned short end1,
                   Trajectory& tj2,
                   unsigned short end2,
                   unsigned short nPtsFit,
                   bool useChg,
                   bool prt)
  {
    // returns the significance of a potential kink between the ends of two trajectories. This
    // is used when deciding to either merge trajectories or make a vertex between them

    if (tj1.CTP != tj2.CTP) return -1;
    if (end1 > 1 || end2 > 1) return -1;

    // construct a temporary trajectory to allow using the standard KinkSignificance function.
    // The first nPtsFit points are comprised of TPs from tj1 and the last nPtsFits points are from tj2
    Trajectory tj;
    tj.ID = 666;
    tj.CTP = tj1.CTP;
    short dir = 1;
    if (end1 == 1) dir = -1;
    unsigned short cnt = 0;
    // add tj1 points to the trajectory
    for (short ii = 0; ii < (short)tj1.Pts.size(); ++ii) {
      short ipt = tj1.EndPt[end1] + dir * ii;
      if (ipt < 0) break;
      if (ipt >= (short)tj1.Pts.size()) break;
      auto& tp = tj1.Pts[ipt];
      if (tp.Chg <= 0) continue;
      tj.Pts.push_back(tp);
      ++cnt;
      if (cnt == nPtsFit + 1) break;
    } // ipt
    if (cnt < nPtsFit) return -1;
    // add tj2 points to the trajectory
    dir = 1;
    if (end2 == 1) dir = -1;
    cnt = 0;
    for (short ii = 0; ii < (short)tj2.Pts.size(); ++ii) {
      short ipt = tj2.EndPt[end2] + dir * ii;
      if (ipt < 0) break;
      if (ipt >= (short)tj2.Pts.size()) break;
      auto& tp = tj2.Pts[ipt];
      if (tp.Chg <= 0) continue;
      tj.Pts.push_back(tp);
      ++cnt;
      if (cnt == nPtsFit + 1) break;
    } // ipt
    tj.EndPt[0] = 0;
    tj.EndPt[1] = tj.Pts.size() - 1;
    return KinkSignificance(slc, tj, nPtsFit, nPtsFit, useChg, prt);
  } // KinkSignificance

  ////////////////////////////////////////////////
  float
  KinkSignificance(TCSlice& slc,
                   Trajectory& tj,
                   unsigned short kinkPt,
                   unsigned short nPtsFit,
                   bool useChg,
                   bool prt)
  {
    // returns a kink significance in the trajectory at the presumed kink point kinkPt
    // using angle and (optional) charge asymmetry. The returned value is negative if there is insufficient
    // information.
    //
    // Check the limits
    if (kinkPt < tj.EndPt[0] + 2) return -1;
    if (kinkPt > tj.EndPt[1] - 2) return -1;

    // This function requires knowledge of the DOF of the line fit
    if (nPtsFit < 3) return -1;
    unsigned short npwc = NumPtsWithCharge(slc, tj, false);
    // need enough points to do a fit on each side of the presumed kink point
    if (npwc < 2 * nPtsFit + 1) return -1;

    // The hit charge uncertainty is 0.12 - 0.15 (neglecting 2ndry interactions) for hadrons.
    // This translates into an error on the charge
    // asymmetry of about 0.07, or about 0.6 * the charge uncertainty
    double chgRMS = 0.07;
    // An additional contribution to the rms is the dependence on the DOF of the fit.
    // Apply a factor to the significance similar to (and simpler than) the Students t-distribution
    // This will increase the angle and charge error rms by 1.3 (1.05) when nPtsFit = 3 (8)
    double tFactor = 1 + 0.3 / double(nPtsFit - 2);
    chgRMS *= tFactor;

    // Set a lower bound on the error expected from the fit, using some assumptions
    double minAngErr = 0.;
    if(tcc.useAlg[kNewCuts]) minAngErr = tcc.hitErrFac * evt.aveHitRMS[0] * tcc.unitsPerTick / nPtsFit;

    // Fit the trajectory direction on the + side
    short fitDir = 1;
    TrajPoint tpPos;
    FitTraj(slc, tj, kinkPt, nPtsFit, fitDir, tpPos);
    if (tpPos.FitChi > 900) return -1;
    // repeat the trajectory fit on the - side
    fitDir = -1;
    TrajPoint tpNeg;
    FitTraj(slc, tj, kinkPt, nPtsFit, fitDir, tpNeg);
    if (tpNeg.FitChi > 900) return -1;
    double angErr = tpNeg.AngErr;
    if (tpPos.AngErr > angErr) angErr = tpPos.AngErr;
    if(angErr < minAngErr) angErr = minAngErr;
    angErr *= tFactor;
    double dang = DeltaAngle(tpPos.Ang, tpNeg.Ang);
    double dangSig = dang / angErr;

    double chgAsym = 0;
    double chgSig = 0;
    if (useChg) {
      // Sum the charge Neg and Pos, excluding the kinkPt
      double chgNeg = 0;
      unsigned short cntNeg = 0;
      for (short ipt = kinkPt - 1; ipt >= tj.EndPt[0]; --ipt) {
        if(ipt < tj.EndPt[0]) break;
        auto& tp = tj.Pts[ipt];
        if (tp.Chg <= 0) continue;
        chgNeg += tp.Chg;
        ++cntNeg;
        if (cntNeg == nPtsFit) break;
        if (ipt == 0) break;
      } // ipt
      if (cntNeg != nPtsFit) return -1;
      // now Pos
      double chgPos = 0;
      unsigned short cntPos = 0;
      for (unsigned short ipt = kinkPt + 1; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if (tp.Chg <= 0) continue;
        chgPos += tp.Chg;
        ++cntPos;
        if (cntPos == nPtsFit) break;
      } // ipt
      if (cntPos != nPtsFit) return -1;
      chgNeg /= (float)nPtsFit;
      chgPos /= (float)nPtsFit;
      // The charge asymmetry varies between 0 and 1;
      chgAsym = std::abs(chgPos - chgNeg) / (chgPos + chgNeg);
      // calculate the charge asymmetry significance
      chgSig = chgAsym / chgRMS;
    } // useChg
    double kinkSig = sqrt(dangSig * dangSig + chgSig * chgSig);

    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "KL: T" << tj.ID << " kinkPt " << PrintPos(slc, tj.Pts[kinkPt]);
      myprt << " nPtsFit " << nPtsFit;
      myprt << " dang " << std::fixed << std::setprecision(3) << dang;
      myprt << std::fixed << std::setprecision(3) << " angErr " << angErr;
      myprt << std::setprecision(2) << " sig " << dangSig;
      myprt << " chgAsym " << chgAsym;
      myprt << " chgSig " << chgSig;
      myprt << " kinkSig " << kinkSig;
      myprt << " minAngErr " << std::setprecision(3) << minAngErr;
    }
    return (float)kinkSig;
  } // KinkSignificance

  ////////////////////////////////////////////////
  bool
  IsShowerLike(TCSlice& slc, const std::vector<int> TjIDs)
  {
    // Vote for the list of Tjs (assumed associated with a PFParticle) being shower-like
    if (TjIDs.empty()) return false;
    unsigned short cnt = 0;
    for (auto tid : TjIDs) {
      if (tid <= 0 || tid > (int)slc.tjs.size()) continue;
      if (slc.tjs[tid - 1].PDGCode == 11) ++cnt;
    } // tjid
    return (cnt > 1);
  } // IsInShower

  ////////////////////////////////////////////////
  float
  ElectronLikelihood(const TCSlice& slc, const Trajectory& tj)
  {
    // returns a number between 0 (not electron-like) and 1 (electron-like)
    if (NumPtsWithCharge(slc, tj, false) < 8) return -1;
    if (tj.EndFlag[0][kEndBragg] || tj.EndFlag[1][kEndBragg]) return 0;

    unsigned short midPt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
    double rms0 = 0, rms1 = 0;
    unsigned short cnt;
    TjDeltaRMS(slc, tj, tj.EndPt[0], midPt, rms0, cnt);
    TjDeltaRMS(slc, tj, midPt, tj.EndPt[1], rms1, cnt);
    float asym = std::abs(rms0 - rms1) / (rms0 + rms1);
    float chgFact = (tj.ChgRMS - 0.1) * 5;
    float elh = 5 * asym * chgFact;
    if (elh > 1) elh = 1;
    return elh;
  } // ElectronLikelihood

  ////////////////////////////////////////////////
  float
  ChgFracNearPos(const TCSlice& slc, const Point2_t& pos, const std::vector<int>& tjIDs)
  {
    // returns the fraction of the charge in the region around pos that is associated with
    // the list of Tj IDs
    if (tjIDs.empty()) return 0;
    std::array<int, 2> wireWindow;
    Point2_t timeWindow;
    // 1/2 size of the region
    constexpr float NNDelta = 5;
    wireWindow[0] = pos[0] - NNDelta;
    wireWindow[1] = pos[0] + NNDelta;
    timeWindow[0] = pos[1] - NNDelta;
    timeWindow[1] = pos[1] + NNDelta;
    // do some checking
    for (auto& tjID : tjIDs)
      if (tjID <= 0 || tjID > (int)slc.tjs.size()) return 0;
    // Determine which plane we are in
    geo::PlaneID planeID = DecodeCTP(slc.tjs[tjIDs[0] - 1].CTP);
    // get a list of all hits in this region
    bool hitsNear;
    std::vector<unsigned int> closeHits =
      FindCloseHits(slc, wireWindow, timeWindow, planeID.Plane, kAllHits, true, hitsNear);
    if (closeHits.empty()) return 0;
    float chg = 0;
    float tchg = 0;
    // Add the hit charge in the box
    // All hits in the box, and all hits associated with the Tjs
    for (auto& iht : closeHits) {
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      chg += hit.Integral();
      if (slc.slHits[iht].InTraj == 0) continue;
      if (std::find(tjIDs.begin(), tjIDs.end(), slc.slHits[iht].InTraj) != tjIDs.end())
        tchg += hit.Integral();
    } // iht
    if (chg == 0) return 0;
    return tchg / chg;
  } // ChgFracNearPos

  ////////////////////////////////////////////////
  float
  MaxHitDelta(TCSlice& slc, Trajectory& tj)
  {
    float delta, md = 0;
    unsigned short ii;
    unsigned int iht;
    for (auto& tp : tj.Pts) {
      for (ii = 0; ii < tp.Hits.size(); ++ii) {
        if (!tp.UseHit[ii]) continue;
        iht = tp.Hits[ii];
        delta = PointTrajDOCA(slc, iht, tp);
        if (delta > md) md = delta;
      } // ii
    }   // pts
    return md;
  } // MaxHitDelta

  //////////////////////////////////////////
  void
  ReverseTraj(TCSlice& slc, Trajectory& tj)
  {
    // reverse the trajectory
    if (tj.Pts.empty()) return;
    // reverse the crawling direction flag
    tj.StepDir = -tj.StepDir;
    // Vertices
    std::swap(tj.VtxID[0], tj.VtxID[1]);
    // trajectory points
    std::reverse(tj.Pts.begin(), tj.Pts.end());
    // reverse the stop flag
    std::swap(tj.EndFlag[0], tj.EndFlag[1]);
    std::swap(tj.dEdx[0], tj.dEdx[1]);
    // reverse the direction vector on all points
    for (unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if (tj.Pts[ipt].Dir[0] != 0) tj.Pts[ipt].Dir[0] = -tj.Pts[ipt].Dir[0];
      if (tj.Pts[ipt].Dir[1] != 0) tj.Pts[ipt].Dir[1] = -tj.Pts[ipt].Dir[1];
      if (tj.Pts[ipt].Ang > 0) { tj.Pts[ipt].Ang -= M_PI; }
      else {
        tj.Pts[ipt].Ang += M_PI;
      }
    } // ipt
    if (tj.StartEnd == 0 || tj.StartEnd == 1) tj.StartEnd = 1 - tj.StartEnd;
    SetEndPoints(tj);
    //    UpdateMatchStructs(slc, tj.ID, tj.ID);
  } // ReverseTraj

  //////////////////////////////////////////
  bool
  PointInsideEnvelope(const Point2_t& Point, const std::vector<Point2_t>& Envelope)
  {
    // returns true if the Point is within the Envelope polygon. Entries in Envelope are the
    // Pos[0], Pos[1] locations of the polygon vertices. This is based on the algorithm that the
    // sum of the angles of a vector between a point and the vertices will be 2 * pi for an interior
    // point and 0 for an exterior point

    Point2_t p1, p2;
    unsigned short nvx = Envelope.size();
    double angleSum = 0;
    for (unsigned short ii = 0; ii < Envelope.size(); ++ii) {
      p1[0] = Envelope[ii][0] - Point[0];
      p1[1] = Envelope[ii][1] - Point[1];
      p2[0] = Envelope[(ii + 1) % nvx][0] - Point[0];
      p2[1] = Envelope[(ii + 1) % nvx][1] - Point[1];
      angleSum += DeltaAngle(p1, p2);
    }
    if (abs(angleSum) < M_PI) return false;
    return true;

  } // InsideEnvelope

  //////////////////////////////////////////
  bool
  SetMag(Vector2_t& v1, double mag)
  {
    double den = v1[0] * v1[0] + v1[1] * v1[1];
    if (den == 0) return false;
    den = sqrt(den);

    v1[0] *= mag / den;
    v1[1] *= mag / den;
    return true;
  } // SetMag

  ////////////////////////////////////////////////
  void
  FindAlongTrans(Point2_t pos1, Vector2_t dir1, Point2_t pos2, Point2_t& alongTrans)
  {
    // Calculate the distance along and transverse to the direction vector dir1 from pos1 to pos2
    alongTrans[0] = 0;
    alongTrans[1] = 0;
    if (pos1[0] == pos2[0] && pos1[1] == pos2[1]) return;
    pos1[0] = pos2[0] - pos1[0];
    pos1[1] = pos2[1] - pos1[1];
    double sep = sqrt(pos1[0] * pos1[0] + pos1[1] * pos1[1]);
    if (sep < 1E-6) return;
    Vector2_t ptDir;
    ptDir[0] = pos1[0] / sep;
    ptDir[1] = pos1[1] / sep;
    SetMag(dir1, 1.0);
    double costh = DotProd(dir1, ptDir);
    if(costh < -1.) costh = -1.;
    if(costh > 1.) costh =  1.;
    alongTrans[0] = costh * sep;
    double sinth = sqrt(1 - costh * costh);
    alongTrans[1] = sinth * sep;
  } // FindAlongTrans

  //////////////////////////////////////////
  double
  DeltaAngle(const Point2_t& p1, const Point2_t& p2)
  {
    // angle between two points
    double ang1 = atan2(p1[1], p1[0]);
    double ang2 = atan2(p2[1], p2[0]);
    return DeltaAngle2(ang1, ang2);
  } // DeltaAngle

  //////////////////////////////////////////
  double
  DeltaAngle2(double Ang1, double Ang2)
  {
    constexpr double twopi = 2 * M_PI;
    double dang = Ang1 - Ang2;
    while (dang > M_PI)
      dang -= twopi;
    while (dang < -M_PI)
      dang += twopi;
    return dang;
  }

  //////////////////////////////////////////
  double
  DeltaAngle(double Ang1, double Ang2)
  {
    return std::abs(std::remainder(Ang1 - Ang2, M_PI));
  }

  ////////////////////////////////////////////////
  void
  SetEndPoints(Trajectory& tj)
  {
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge

    tj.EndPt[0] = 0;
    tj.EndPt[1] = 0;
    if (tj.Pts.size() == 0) return;

    // check the end point pointers
    for (unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if (tj.Pts[ipt].Chg != 0) {
        tj.EndPt[0] = ipt;
        break;
      }
    }
    for (unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      if (tj.Pts[ipt].Chg != 0) {
        tj.EndPt[1] = ipt;
        break;
      }
    }
  } // SetEndPoints

  ////////////////////////////////////////////////
  bool
  TrajIsClean(TCSlice& slc, Trajectory& tj, bool prt)
  {
    // Returns true if the trajectory has low hit multiplicity and is in a
    // clean environment
    unsigned short nUsed = 0;
    unsigned short nTotHits = 0;
    for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      nTotHits += tp.Hits.size();
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if (tp.UseHit[ii]) ++nUsed;
      } // ii
    }   // ipt
    if (nTotHits == 0) return false;
    float fracUsed = (float)nUsed / (float)nTotHits;
    if (prt)
      mf::LogVerbatim("TC") << "TrajIsClean: nTotHits " << nTotHits << " nUsed " << nUsed
                            << " fracUsed " << fracUsed;

    if (fracUsed > 0.9) return true;
    return false;

  } // TrajIsClean

  ////////////////////////////////////////////////
  short
  MCSMom(const TCSlice& slc, const std::vector<int>& tjIDs)
  {
    // Find the average MCSMom of the trajectories
    if (tjIDs.empty()) return 0;
    float summ = 0;
    float suml = 0;
    for (auto tjid : tjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.MCSMom<= 0) continue;
      float npts = tj.EndPt[1] - tj.EndPt[0] + 1;
      summ += npts * tj.MCSMom;
      suml += npts;
    } // tjid
    if(suml == 0) return 0;
    return (short)(summ / suml);
  } // MCSMom

  ////////////////////////////////////////////////
  short
  MCSMom(const TCSlice& slc, const Trajectory& tj)
  {
    // ignore the points near an end if there is a Bragg peak.
    // Muons and pions wander significantly near the stopping point
    // but are otherwise fairly straight
    unsigned short fromPt = tj.EndPt[0];
    unsigned short toPt = tj.EndPt[1];
    if(tcc.useAlg[kNewCuts] && tj.Pts.size() > 30) {
      if(tj.EndFlag[0][kEndBragg]) {
        fromPt += 10;
      } else if(tj.EndFlag[1][kEndBragg]) {
        toPt -= 10;
      }
    } // NewCuts
    return MCSMom(slc, tj, fromPt, toPt);
  } // MCSMom

  ////////////////////////////////////////////////
  short
  MCSMom(const TCSlice& slc, const Trajectory& tj, unsigned short firstPt, unsigned short lastPt)
  {
    // Estimate the trajectory momentum using Multiple Coulomb Scattering ala PDG RPP

    if (firstPt == lastPt) return 0;
    if (firstPt > lastPt) std::swap(firstPt, lastPt);

    firstPt = NearestPtWithChg(slc, tj, firstPt);
    lastPt = NearestPtWithChg(slc, tj, lastPt);
    if (firstPt >= lastPt) return 0;

    if (firstPt < tj.EndPt[0]) return 0;
    if (lastPt > tj.EndPt[1]) return 0;
    // Can't do this with only 2 points
    if (NumPtsWithCharge(slc, tj, false, firstPt, lastPt) < 3) return 0;
    // Ignore junk Tjs
    if (tj.AlgMod[kJunkTj]) return 0;

    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if (tjLen < 1) return 0;
    // mom calculated in MeV
    double thetaRMS = MCSThetaRMS(slc, tj, firstPt, lastPt);
    if (thetaRMS < 0.001) return 999;
    double mom = 13.8 * sqrt(tjLen / 14) / thetaRMS;
    if (mom > 999) mom = 999;
    return (short)mom;
  } // MCSMom

  ////////////////////////////////////////////////
  unsigned short
  NearestPtWithChg(const TCSlice& slc, const Trajectory& tj, unsigned short thePt)
  {
    // returns a point near thePt which has charge
    if (thePt > tj.EndPt[1]) return thePt;
    if (tj.Pts[thePt].Chg > 0) return thePt;

    short endPt0 = tj.EndPt[0];
    short endPt1 = tj.EndPt[1];
    for (short off = 1; off < 10; ++off) {
      short ipt = thePt + off;
      if (ipt <= endPt1 && tj.Pts[ipt].Chg > 0) return (unsigned short)ipt;
      ipt = thePt - off;
      if (ipt >= endPt0 && tj.Pts[ipt].Chg > 0) return (unsigned short)ipt;
    } // off
    return thePt;
  } // NearestPtWithChg

  /////////////////////////////////////////
  float
  MCSThetaRMS(const TCSlice& slc, const Trajectory& tj)
  {
    // This returns the MCS scattering angle expected for one WSE unit of travel along the trajectory.
    // It is used to define kink and vertex cuts. This should probably be named something different to
    // prevent confusion

    float tps = TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
    if (tps < 1) return 1;

    return MCSThetaRMS(slc, tj, tj.EndPt[0], tj.EndPt[1]) / sqrt(tps);

  } // MCSThetaRMS

  /////////////////////////////////////////
  double
  MCSThetaRMS(const TCSlice& slc,
              const Trajectory& tj,
              unsigned short firstPt,
              unsigned short lastPt)
  {
    // This returns the MCS scattering angle expected for the length of the trajectory
    // spanned by firstPt to lastPt. It is used primarily to calculate MCSMom

    if (firstPt < tj.EndPt[0]) return 1;
    if (lastPt > tj.EndPt[1]) return 1;

    firstPt = NearestPtWithChg(slc, tj, firstPt);
    lastPt = NearestPtWithChg(slc, tj, lastPt);
    if (firstPt >= lastPt) return 1;

    double sigmaS;
    unsigned short cnt;
    TjDeltaRMS(slc, tj, firstPt, lastPt, sigmaS, cnt);
    if (sigmaS < 0) return 1;
    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if (tjLen < 1) return 1;
    // Theta_o =  4 * sqrt(3) * sigmaS / path
    return (6.8 * sigmaS / tjLen);

  } // MCSThetaRMS

  /////////////////////////////////////////
  void
  TjDeltaRMS(const TCSlice& slc,
             const Trajectory& tj,
             unsigned short firstPt,
             unsigned short lastPt,
             double& rms,
             unsigned short& cnt)
  {
    // returns the rms scatter of points around a line formed by the firstPt and lastPt of the trajectory

    rms = -1;
    if (firstPt < tj.EndPt[0]) return;
    if (lastPt > tj.EndPt[1]) return;

    firstPt = NearestPtWithChg(slc, tj, firstPt);
    lastPt = NearestPtWithChg(slc, tj, lastPt);
    if (firstPt >= lastPt) return;

    TrajPoint tmp;
    // make a bare trajectory point to define a line between firstPt and lastPt.
    // Use the position of the hits at these points
    TrajPoint firstTP = tj.Pts[firstPt];
    firstTP.Pos = firstTP.HitPos;
    TrajPoint lastTP = tj.Pts[lastPt];
    lastTP.Pos = lastTP.HitPos;
    if (!MakeBareTrajPoint(slc, firstTP, lastTP, tmp)) return;
    // sum up the deviations^2
    double dsum = 0;
    cnt = 0;
    for (unsigned short ipt = firstPt + 1; ipt < lastPt; ++ipt) {
      if (tj.Pts[ipt].Chg == 0) continue;
      // ignore points with large error
      if (tj.Pts[ipt].HitPosErr2 > 4) continue;
      dsum += PointTrajDOCA2(slc, tj.Pts[ipt].HitPos[0], tj.Pts[ipt].HitPos[1], tmp);
      ++cnt;
    } // ipt
    if (cnt < 2) return;
    rms = sqrt(dsum / (double)cnt);

  } // TjDeltaRMS

  /////////////////////////////////////////
  void
  SetTPEnvironment(TCSlice& slc, CTP_t inCTP)
  {
    // This function is called after tj reconstruction is completed to set TP Environment
    // bits that are dependent on reconstruction, just kEnvNearMuon for now. This bit is
    // set for all TPs that are within 5 wire-equivalents of a muon

    std::array<int, 2> wireWindow;
    Point2_t timeWindow;
    unsigned short plane = DecodeCTP(inCTP).Plane;
    float delta = 5;

    for (auto& mutj : slc.tjs) {
      if (mutj.AlgMod[kKilled]) continue;
      if (mutj.CTP != inCTP) continue;
      if (mutj.PDGCode != 13) continue;
      for (auto& mutp : mutj.Pts) {
        if(mutp.Chg <= 0) continue;
        wireWindow[0] = mutp.Pos[0];
        wireWindow[1] = mutp.Pos[0];
        timeWindow[0] = mutp.Pos[1] - delta;
        timeWindow[1] = mutp.Pos[1] + delta;
        // get a list of all hits in this region
        bool hitsNear;
        auto closeHits =
          FindCloseHits(slc, wireWindow, timeWindow, plane, kAllHits, true, hitsNear);
        if (closeHits.empty()) continue;
        for (auto iht : closeHits) {
          auto inTraj = slc.slHits[iht].InTraj;
          if (inTraj <= 0 || inTraj > (int)slc.tjs.size()) continue;
          if (inTraj == mutj.ID) continue;
          auto& dtj = slc.tjs[inTraj - 1];
          if (dtj.AlgMod[kKilled]) continue;
          if (dtj.PDGCode == 13) continue;
          for (auto& dtp : dtj.Pts) {
            if(dtp.Chg <= 0) continue;
            if (std::find(dtp.Hits.begin(), dtp.Hits.end(), iht) == dtp.Hits.end()) continue;
            dtp.Environment[kEnvNearMuon] = true;
          } // jpt
        }   // iht
      }     // mutp
    }       // mutj
  }         // SetTPEnvironment

  /////////////////////////////////////////
  void
  UpdateTjChgProperties(std::string inFcnLabel, TCSlice& slc, Trajectory& tj, bool prt)
  {
    // Updates properties of the tj that are affected when the TP environment
    // is changed. The most likely reason for a change is when the tj is attached to a
    // vertex in which case the Environment kEnvOverlap bit may be set by the UpdateVxEnvironment
    // function in which case this function is called.
    if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) return;

    // first (un)set some bits
    for (auto& tp : tj.Pts) {
      if (tp.Chg <= 0) continue;
      tp.Environment[kEnvUnusedHits] = false;
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if (tp.UseHit[ii]) continue;
        if (slc.slHits[tp.Hits[ii]].InTraj == 0) tp.Environment[kEnvUnusedHits] = true;
      } // ii
    }   // tp

    // Update the tj charge variables. The concept is explained by this graphic where
    // each column is a wire, Q = a TP with charge, q = a TP with charge that is an
    // EnvOverlap region, x = a wire that has a TP with Chg = 0 or a wire that has no TP
    // because the wire is dead, o = an EnvOverlap region, V = vertex attached to end. You should
    // imagine that all 3 tjs come from the same vertex
    //   01234567890123456789   npwc  cnt range
    //   VooooQQQQxxxQQQ          7    7   0 - 14
    //   VqqqqQQQQxxxQQQQQQQQ    16   12   0 - 19
    //   VooQQQ                   3    3   0 - 5
    // The average is first calculated using Ave = sum(Q) / npwc
    tj.TotChg = 0;
    tj.AveChg = 0;
    tj.ChgRMS = 0.5;

    // Reject a single large charge TP if the Tj isn't short
    float bigChg = 0;
    unsigned short skipPt = USHRT_MAX;
    // count the number of overlap TPs
    short nOverlap = 0;
    // variables used to calculate backup averages in case the Tj is short
    // or many points are overlapping with a different trajectory
    short cnt = 0;
    double sum = 0;
    double sum2 = 0;
    for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.Chg <= 0) continue;
      if (tj.Pts[ipt].Environment[kEnvOverlap]) ++nOverlap;
      if (tp.Chg > bigChg) {
        bigChg = tp.Chg;
        skipPt = ipt;
      } // tp.Chg > bigChg
      ++cnt;
      sum += tp.Chg;
      sum2 += tp.Chg * tp.Chg;
    } // ipt
    // there must be at least 3 TPs with charge in a trajectory
    if(cnt < 3) return;
    // calculate the average and RMS into the Tj
    tj.AveChg = sum / (double)cnt;
    double arg = sum2 - cnt * tj.AveChg * tj.AveChg;
    double rms = 0.5;
    if (arg > 0) rms = sqrt(arg / (cnt - 1));
    // normalize with the average
    rms /= tj.AveChg;
    // don't let it be an unrealistically low value. It could be crazy large however.
    if (rms < 0.1) rms = 0.1;
    // Don't let the calculated charge RMS dominate until it is well known - after there are 10 TPs.
    if (cnt < 10) {
      double defFrac = 1 / cnt;
      rms = defFrac * 0.5 + (1 - defFrac) * rms;
    }
    tj.ChgRMS = rms;
    tj.TotChg = cnt * tj.AveChg;
    if(tj.AlgMod[kJunkTj]) {
      tj.Pts[tj.EndPt[0]].AveChg = tj.AveChg;
      tj.Pts[tj.EndPt[1]].AveChg = tj.AveChg;
      return;
    }

    // See if there are enough points to find a better average and RMS by
    // ignoring the end points (which may have a Bragg peak or may be low) 
    // and ignoring the overlapping TPs
    short nTPsToUse = cnt - 2 - nOverlap;
    // and maybe ignore an anomalously high-charge point if the point is
    // not at an end
    if(skipPt != tj.EndPt[0] && skipPt != tj.EndPt[1]) {
      float bigPull = (bigChg / tj.AveChg - 1) / tj.ChgRMS;
      if(bigPull > 2.5 && nTPsToUse > 3) {
        --nTPsToUse;
      } else {
        skipPt = USHRT_MAX;
      }
    } // skipPt != ...
    // Can't do any better
    if(nTPsToUse < 3) return;

    cnt = 0;
    sum = 0;
    sum2 = 0;
    // these are used as a backup sum in case there are an insufficient number of points
    // after the cuts.
    // don't include the end points in the average
    for (unsigned short ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.Chg <= 0) continue;
      // ignore the single large charge TP
      if (ipt == skipPt) continue;
      // Skip TPs that overlap with TPs on other Tjs
      if (tj.Pts[ipt].Environment[kEnvOverlap]) continue;
      double tpchg = 0;
      for (unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if (!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        tpchg += (*evt.allHits)[slc.slHits[iht].allHitsIndex].Integral();
      } // ii
      ++cnt;
      sum += tpchg;
      sum2 += tpchg * tpchg;
    } // ipt
    tj.AveChg = sum / cnt;
    // calculate the total charge using the tj wire range
    double nWires = std::abs(tj.Pts[tj.EndPt[1]].Pos[0] - tj.Pts[tj.EndPt[0]].Pos[0]);
    tj.TotChg = nWires * tj.AveChg;
    // calculate the rms
    arg = sum2 - cnt * tj.AveChg * tj.AveChg;
    rms = 0.5;
    if (arg > 0) rms = sqrt(arg / (cnt - 1));
    rms /= tj.AveChg;
    // don't let it be an unrealistically low value. It could be crazy large however.
    if (rms < 0.1) rms = 0.1;
    // Don't let the calculated charge RMS dominate until it is well known; after there are 5 - 10 valid TPs.
    // Set the starting charge rms = 0.5
    if (cnt < 10) {
      double defFrac = 1 / cnt;
      rms = defFrac * 0.5 + (1 - defFrac) * rms;
    }
    tj.ChgRMS = rms;
    if (prt)
      mf::LogVerbatim("TC") << inFcnLabel << ".UpdateTjChgProperties: Set tj.AveChg "
                            << (int)tj.AveChg << " ChgRMS " << tj.ChgRMS;

    // Update the TP charge pulls.
    // Don't let the calculated charge RMS dominate the default
    // RMS until it is well known. Start with 50% error on the
    // charge RMS
    for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.Chg <= 0) continue;
      tp.ChgPull = (tp.Chg / tj.AveChg - 1) / tj.ChgRMS;
    } // ipt

    // update the local charge average using NPtsAve of the preceding points.
    // Handle short Tjs first.
    if (cnt < tcc.nPtsAve) {
      for (auto& tp : tj.Pts)
        tp.AveChg = tj.AveChg;
      return;
    }

    // Set the local average to 0 first
    for (auto& tp : tj.Pts)
      tp.AveChg = 0;
    // Enter the local average on the points where an average can be calculated
    unsigned short nptsave = tcc.nPtsAve;
    unsigned short minPt = tj.EndPt[0] + nptsave;
    float lastAve = 0;
    for (unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.EndPt[1] - ii;
      if (ipt < minPt) break;
      float cnt = 0;
      float sum = 0;
      for (unsigned short iii = 0; iii < nptsave; ++iii) {
        unsigned short iipt = ipt - iii;
        // Don't include the charge of the first point
        if (iipt == tj.EndPt[0]) break;
        auto& tp = tj.Pts[iipt];
        if (tp.Chg <= 0) continue;
        sum += tp.Chg;
        ++cnt;
      } // iii
      if (cnt > 2) {
        tj.Pts[ipt].AveChg = sum / cnt;
        lastAve = tj.Pts[ipt].AveChg;
      }
    } // ii
    // Fill in the points where no average was calculated
    for (unsigned short ii = tj.EndPt[0]; ii <= tj.EndPt[1]; ++ii) {
      unsigned short ipt = tj.EndPt[1] - ii;
      auto& tp = tj.Pts[ipt];
      if (tp.AveChg == 0) { tp.AveChg = lastAve; }
      else {
        lastAve = tp.AveChg;
      }
    } // ii

    tj.NeedsUpdate = false;

  } // UpdateTjChgProperties

  /////////////////////////////////////////
  void
  UpdateVxEnvironment(TCSlice& slc)
  {
    // Set the kEnvOverlap bit true for all TPs that are close to other
    // trajectories that are close to vertices. The positions of TPs that
    // overlap are biased and shouldn't be used in a vertex fit. Also, these
    // TPs shouldn't be used to calculate dE/dx. The kEnvOverlap bit is first cleared
    // for ALL TPs and then set for ALL 2D vertices

    for (auto& tj : slc.tjs) {
      if (tj.AlgMod[kKilled]) continue;
      for (auto& tp : tj.Pts)
        tp.Environment[kEnvOverlap] = false;
    } // tj

    for (auto& vx : slc.vtxs) {
      if (vx.ID <= 0) continue;
      UpdateVxEnvironment(slc, vx, false);
    } // vx

  } // UpdateVxEnvironment

  /////////////////////////////////////////
  void
  UpdateVxEnvironment(TCSlice& slc, VtxStore& vx2, bool prt)
  {
    // Update the Environment each TP on trajectories near the vertex

    if (vx2.ID == 0) return;
    if (vx2.Stat[kOnDeadWire]) return;


    std::vector<int> tjlist;
    std::vector<unsigned short> tjends;
    if (vx2.Pos[0] < -0.4) return;
    unsigned int vxWire = std::nearbyint(vx2.Pos[0]);
    unsigned int loWire = vxWire;
    unsigned int hiWire = vxWire;
    for (auto& tj : slc.tjs) {
      if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      if (tj.CTP != vx2.CTP) continue;
      // ignore photon Tjs
      if (tj.AlgMod[kPhoton]) continue;
      for (unsigned short end = 0; end < 2; ++end) {
        if (tj.VtxID[end] != vx2.ID) continue;
        tjlist.push_back(tj.ID);
        tjends.push_back(end);
        if (tj.Pts[tj.EndPt[end]].Pos[0] < -0.4) return;
        unsigned int endWire = std::nearbyint(tj.Pts[tj.EndPt[end]].Pos[0]);
        if (endWire < loWire) loWire = endWire;
        if (endWire > hiWire) hiWire = endWire;
      } // end
    }   // tj
    if (tjlist.size() < 2) return;
    if (hiWire < loWire + 1) return;
    if (prt)
      mf::LogVerbatim("TC") << " check Tjs on wires in the range " << loWire << " to " << hiWire;

    // create a vector of TPs between loWire and hiWire for every tj in the list
    //   wire       TP
    std::vector<std::vector<TrajPoint>> wire_tjpt;
    // companion vector of IDs
    std::vector<int> tjids;
    // populate this vector with TPs on Tjs that are in this range
    unsigned short nwires = hiWire - loWire + 1;
    for (unsigned short itj = 0; itj < tjlist.size(); ++itj) {
      auto& tj = slc.tjs[tjlist[itj] - 1];
      unsigned short end = tjends[itj];
      std::vector<TrajPoint> tjpt(nwires);
      // first enter valid TPs in the range
      for (unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt;
        if (end == 0) { ipt = tj.EndPt[0] + ii; }
        else {
          ipt = tj.EndPt[1] - ii;
        }
        if (ipt > tj.Pts.size() - 1) break;
        // Make a copy of the TP so we can alter it
        auto tp = tj.Pts[ipt];
        if (tp.Chg <= 0) continue;
        tp.Chg = 1;
        tp.Hits.clear();
        if (tp.Pos[0] < -0.4) continue;
        unsigned int wire = std::nearbyint(tp.Pos[0]);
        unsigned short indx = wire - loWire;
        if (indx > nwires - 1) break;
        tp.Step = ipt;
        // We will use NTPsFit to count the number of neighboring TPs
        tp.NTPsFit = 0;
        tjpt[indx] = tp;
      } // ii
      // next make TPs on the wires that don't have real TPs
      TrajPoint ltp;
      // put ltp at the vertex position with direction towards the end point
      MakeBareTrajPoint(vx2.Pos, tj.Pts[tj.EndPt[end]].Pos, ltp);
      if (ltp.Dir[0] == 0) continue;
      if (ltp.Pos[0] < -0.4) continue;
      unsigned int wire = std::nearbyint(ltp.Pos[0]);
      ltp.Chg = 0;
      unsigned short indx = wire - loWire;
      // Break if we found a real TP
      if (tjpt[indx].Chg == 0) tjpt[indx] = ltp;
      double stepSize = std::abs(1 / ltp.Dir[0]);
      for (unsigned short ii = 0; ii < nwires; ++ii) {
        // move the local TP position by one step in the right direction
        for (unsigned short iwt = 0; iwt < 2; ++iwt)
          ltp.Pos[iwt] += ltp.Dir[iwt] * stepSize;
        if (ltp.Pos[0] < -0.4) break;
        wire = std::nearbyint(ltp.Pos[0]);
        if (wire < loWire || wire > hiWire) break;
        indx = wire - loWire;
        if (tjpt[indx].Chg > 0) continue;
        tjpt[indx] = ltp;
      } // ii
      if (prt) {
        mf::LogVerbatim myprt("TC");
        myprt << " T" << tj.ID;
        for (auto& tp : tjpt)
          myprt << " " << PrintPos(slc, tp.Pos) << "_" << tp.Step << "_" << (int)tp.Chg;
      }
      wire_tjpt.push_back(tjpt);
      tjids.push_back(tj.ID);
    } // itj

    // iterate over the wires in the range
    for (unsigned short indx = 0; indx < nwires; ++indx) {
      // count the number of valid points on this wire
      unsigned short npts = 0;
      // count the number of points on this wire that have charge
      unsigned short npwc = 0;
      for (unsigned short itj = 0; itj < wire_tjpt.size(); ++itj) {
        if (wire_tjpt[itj][indx].Pos[0] == 0) continue;
        // found a valid point
        ++npts;
        if (wire_tjpt[itj][indx].Chg > 0) ++npwc;
      } // itj
      // no valid points
      if (npts == 0) continue;
      // all valid points have charge
      if (npwc == npts) continue;
      // re-find the valid points with charge and set the kEnvOverlap bit
      for (unsigned short itj = 0; itj < wire_tjpt.size(); ++itj) {
        if (wire_tjpt[itj][indx].Pos[0] == 0) continue;
        if (wire_tjpt[itj][indx].Chg == 0) continue;
        auto& tj = slc.tjs[tjids[itj] - 1];
        unsigned short ipt = wire_tjpt[itj][indx].Step;
        tj.Pts[ipt].Environment[kEnvOverlap] = true;
        tj.NeedsUpdate = true;
        if (prt) mf::LogVerbatim("TC") << " Set kEnvOverlap bit on T" << tj.ID 
                << " TP " <<PrintPos(slc, tj.Pts[ipt]);
      } // itj
    }   // indx

    // update the charge rms for those tjs whose environment was changed above
    // (or elsewhere)
    for (auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if (!tj.NeedsUpdate) continue;
      if (tj.CTP != vx2.CTP) continue;
      UpdateTjChgProperties("UVxE", slc, tj, prt);
    } // tjid

  } // UpdateVxEnvironment

  /////////////////////////////////////////
  TrajPoint
  MakeBareTP(detinfo::DetectorPropertiesData const& detProp,
             const TCSlice& slc,
             const Point3_t& pos,
             CTP_t inCTP)
  {
    // A version to use when the 2D direction isn't required
    TrajPoint tp;
    tp.Pos = {{0, 0}};
    tp.Dir = {{0, 1}};
    tp.CTP = inCTP;
    geo::PlaneID planeID = DecodeCTP(inCTP);

    tp.Pos[0] = tcc.geom->WireCoordinate(pos[1], pos[2], planeID);
    tp.Pos[1] = detProp.ConvertXToTicks(pos[0], planeID) * tcc.unitsPerTick;
    return tp;
  } // MakeBareTP

  /////////////////////////////////////////
  TrajPoint
  MakeBareTP(detinfo::DetectorPropertiesData const& detProp,
             const TCSlice& slc,
             const Point3_t& pos,
             const Vector3_t& dir,
             CTP_t inCTP)
  {
    // Projects the space point defined by pos and dir into the CTP and returns
    // it in the form of a trajectory point. The TP Pos[0] is set to a negative
    // number if the point has an invalid wire position but doesn't return an
    // error if the position is on a dead wire. The projection of the direction
    // vector in CTP is stored in tp.Delta.
    TrajPoint tp;
    tp.Pos = {{-1, 0}};
    tp.Dir = {{0, 1}};
    tp.CTP = inCTP;
    geo::PlaneID planeID = DecodeCTP(inCTP);

    tp.Pos[0] = tcc.geom->WireCoordinate(pos[1], pos[2], planeID);
    tp.Pos[1] = detProp.ConvertXToTicks(pos[0], planeID) * tcc.unitsPerTick;

    // now find the direction if dir is defined
    if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) return tp;

    // Make a point at the origin and one 100 units away
    Point3_t ori3 = {{0.0, 0.0, 0.0}};
    Point3_t pos3 = {{100 * dir[0], 100 * dir[1], 100 * dir[2]}};
    // 2D position of ori3 and the pos3 projection
    std::array<double, 2> ori2;
    std::array<double, 2> pos2;
    std::array<double, 2> dir2;
    // the wire coordinates
    ori2[0] = tcc.geom->WireCoordinate(ori3[1], ori3[2], planeID);
    pos2[0] = tcc.geom->WireCoordinate(pos3[1], pos3[2], planeID);
    // the time coordinates
    ori2[1] = detProp.ConvertXToTicks(ori3[0], planeID) * tcc.unitsPerTick;
    pos2[1] = detProp.ConvertXToTicks(pos3[0], planeID) * tcc.unitsPerTick;

    dir2[0] = pos2[0] - ori2[0];
    dir2[1] = pos2[1] - ori2[1];

    double norm = sqrt(dir2[0] * dir2[0] + dir2[1] * dir2[1]);
    tp.Dir[0] = dir2[0] / norm;
    tp.Dir[1] = dir2[1] / norm;
    tp.Ang = atan2(dir2[1], dir2[0]);
    tp.Delta = norm / 100;

    // The Orth vectors are not unit normalized so we need to correct for this
    double w0 = tcc.geom->WireCoordinate(0, 0, planeID);
    // cosine-like component
    double cs = tcc.geom->WireCoordinate(1, 0, planeID) - w0;
    // sine-like component
    double sn = tcc.geom->WireCoordinate(0, 1, planeID) - w0;
    norm = sqrt(cs * cs + sn * sn);
    tp.Delta /= norm;

    // Stasb dt/dWire in DeltaRMS. This is used in PFPUtils/FitSection to find the
    // distance along a 3D line given the wire number in a plane
    tp.DeltaRMS = 100 / (pos2[0] - ori2[0]);
    return tp;

  } // MakeBareTP

  /////////////////////////////////////////
  bool
  MakeBareTrajPoint(const TCSlice& slc, unsigned int fromHit, unsigned int toHit, TrajPoint& tp)
  {
    if (fromHit > slc.slHits.size() - 1) return false;
    if (toHit > slc.slHits.size() - 1) return false;
    auto& fhit = (*evt.allHits)[slc.slHits[fromHit].allHitsIndex];
    auto& thit = (*evt.allHits)[slc.slHits[toHit].allHitsIndex];
    CTP_t tCTP = EncodeCTP(fhit.WireID());
    return MakeBareTrajPoint(slc,
                             (float)fhit.WireID().Wire,
                             fhit.PeakTime(),
                             (float)thit.WireID().Wire,
                             thit.PeakTime(),
                             tCTP,
                             tp);

  } // MakeBareTrajPoint

  /////////////////////////////////////////
  bool
  MakeBareTrajPoint(const TCSlice& slc,
                    float fromWire,
                    float fromTick,
                    float toWire,
                    float toTick,
                    CTP_t tCTP,
                    TrajPoint& tp)
  {
    tp.CTP = tCTP;
    tp.Pos[0] = fromWire;
    tp.Pos[1] = tcc.unitsPerTick * fromTick;
    tp.Dir[0] = toWire - fromWire;
    tp.Dir[1] = tcc.unitsPerTick * (toTick - fromTick);
    double norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    if (norm == 0) return false;
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    return true;
  } // MakeBareTrajPoint

  /////////////////////////////////////////
  bool
  MakeBareTrajPoint(const Point2_t& fromPos, const Point2_t& toPos, TrajPoint& tpOut)
  {
    tpOut.Pos = fromPos;
    tpOut.Dir = PointDirection(fromPos, toPos);
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
    return true;

  } // MakeBareTrajPoint

  /////////////////////////////////////////
  bool
  MakeBareTrajPoint(const TCSlice& slc,
                    const TrajPoint& tpIn1,
                    const TrajPoint& tpIn2,
                    TrajPoint& tpOut)
  {
    tpOut.CTP = tpIn1.CTP;
    tpOut.Pos = tpIn1.Pos;
    tpOut.Dir = PointDirection(tpIn1.Pos, tpIn2.Pos);
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
    return true;
  } // MakeBareTrajPoint

  ////////////////////////////////////////////////
  unsigned short
  FarEnd(TCSlice& slc, const Trajectory& tj, const Point2_t& pos)
  {
    // Returns the end (0 or 1) of the Tj that is furthest away from the position pos
    if (tj.ID == 0) return 0;
    if (PosSep2(tj.Pts[tj.EndPt[1]].Pos, pos) > PosSep2(tj.Pts[tj.EndPt[0]].Pos, pos)) return 1;
    return 0;
  } // FarEnd

  ////////////////////////////////////////////////
  Vector2_t
  PointDirection(const Point2_t p1, const Point2_t p2)
  {
    // Finds the direction vector between the two points from p1 to p2
    Vector2_t dir;
    for (unsigned short xyz = 0; xyz < 2; ++xyz)
      dir[xyz] = p2[xyz] - p1[xyz];
    if (dir[0] == 0 && dir[1] == 0) return dir;
    double norm = sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
    dir[0] /= norm;
    dir[1] /= norm;
    return dir;
  } // PointDirection

  ////////////////////////////////////////////////
  float
  TPHitsRMSTime(const TCSlice& slc, const TrajPoint& tp, HitStatus_t hitRequest)
  {
    return tcc.unitsPerTick * TPHitsRMSTick(slc, tp, hitRequest);
  } // TPHitsRMSTime

  ////////////////////////////////////////////////
  float
  TPHitsRMSTick(const TCSlice& slc, const TrajPoint& tp, HitStatus_t hitRequest)
  {
    // Estimate the RMS of all hits associated with a trajectory point
    // without a lot of calculation. Note that this returns a value that is
    // closer to a FWHM, not the RMS
    if (tp.Hits.empty()) return 0;
    float minVal = 9999;
    float maxVal = 0;
    for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      bool useit = (hitRequest == kAllHits);
      if (hitRequest == kUsedHits && tp.UseHit[ii]) useit = true;
      if (hitRequest == kUnusedHits && !tp.UseHit[ii]) useit = true;
      if (!useit) continue;
      unsigned int iht = tp.Hits[ii];
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      float cv = hit.PeakTime();
      float rms = hit.RMS();
      float arg = cv - rms;
      if (arg < minVal) minVal = arg;
      arg = cv + rms;
      if (arg > maxVal) maxVal = arg;
    } // ii
    if (maxVal == 0) return 0;
    return (maxVal - minVal) / 2;
  } // TPHitsRMSTick

  ////////////////////////////////////////////////
  float
  HitsRMSTime(const TCSlice& slc,
              const std::vector<unsigned int>& hitsInMultiplet,
              HitStatus_t hitRequest)
  {
    return tcc.unitsPerTick * HitsRMSTick(slc, hitsInMultiplet, hitRequest);
  } // HitsRMSTick

  ////////////////////////////////////////////////
  float
  HitsRMSTick(const TCSlice& slc,
              const std::vector<unsigned int>& hitsInMultiplet,
              HitStatus_t hitRequest)
  {
    if (hitsInMultiplet.empty()) return 0;

    if (hitsInMultiplet.size() == 1) {
      auto& hit = (*evt.allHits)[slc.slHits[hitsInMultiplet[0]].allHitsIndex];
      return hit.RMS();
    }

    float minVal = 9999;
    float maxVal = 0;
    for (unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
      unsigned int iht = hitsInMultiplet[ii];
      bool useit = (hitRequest == kAllHits);
      if (hitRequest == kUsedHits && slc.slHits[iht].InTraj > 0) useit = true;
      if (hitRequest == kUnusedHits && slc.slHits[iht].InTraj == 0) useit = true;
      if (!useit) continue;
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      float cv = hit.PeakTime();
      float rms = hit.RMS();
      float arg = cv - rms;
      if (arg < minVal) minVal = arg;
      arg = cv + rms;
      if (arg > maxVal) maxVal = arg;
    } // ii
    if (maxVal == 0) return 0;
    return (maxVal - minVal) / 2;
  } // HitsRMSTick

  ////////////////////////////////////////////////
  float
  HitsPosTime(const TCSlice& slc,
              const std::vector<unsigned int>& hitsInMultiplet,
              float& sum,
              HitStatus_t hitRequest)
  {
    return tcc.unitsPerTick * HitsPosTick(slc, hitsInMultiplet, sum, hitRequest);
  } // HitsPosTime

  ////////////////////////////////////////////////
  float
  HitsPosTick(const TCSlice& slc,
              const std::vector<unsigned int>& hitsInMultiplet,
              float& sum,
              HitStatus_t hitRequest)
  {
    // returns the position and the charge
    float pos = 0;
    sum = 0;
    for (unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
      unsigned int iht = hitsInMultiplet[ii];
      bool useit = (hitRequest == kAllHits);
      if (hitRequest == kUsedHits && slc.slHits[iht].InTraj > 0) useit = true;
      if (hitRequest == kUnusedHits && slc.slHits[iht].InTraj == 0) useit = true;
      if (!useit) continue;
      auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
      float chg = hit.Integral();
      pos += chg * hit.PeakTime();
      sum += chg;
    } // ii
    if (sum <= 0) return -1;
    return pos / sum;
  } // HitsPosTick

  //////////////////////////////////////////
  unsigned short
  NumNotUsedHitsInTP(const TCSlice& slc, const TrajPoint& tp)
  {
    // counts the number of hits not used by the TJ in this TP
    unsigned short notUsed = 0;
    for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if (!tp.UseHit[ii]) ++notUsed;
    return notUsed;
  } // NumNotHitsInTP

  //////////////////////////////////////////
  unsigned short
  NumHitsInTP(const TrajPoint& tp, HitStatus_t hitRequest)
  {
    // Counts the number of hits of the specified type in tp
    if (tp.Hits.empty()) return 0;

    if (hitRequest == kAllHits) return tp.Hits.size();

    unsigned short nhits = 0;
    for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if (hitRequest == kUsedHits) {
        if (tp.UseHit[ii]) ++nhits;
      }
      else {
        // looking for unused hits
        if (!tp.UseHit[ii]) ++nhits;
      }
    } // ii
    return nhits;
  } // NumHitsInTP

  ////////////////////////////////////////////////
  void
  SetPDGCode(TCSlice& slc, unsigned short itj)
  {
    if (itj > slc.tjs.size() - 1) return;
    SetPDGCode(slc, slc.tjs[itj]);
  }

  ////////////////////////////////////////////////
  void
  SetPDGCode(TCSlice& slc, Trajectory& tj)
  {
    // Sets the PDG code for the supplied trajectory. Note that the existing
    // PDG code is left unchanged if these cuts are not met

    // this PDGCode is set in TagShowerLike and shouldn't be altereed here
    if(tj.PDGCode == 11) return;

    short npwc = NumPtsWithCharge(slc, tj, false);
    if (npwc < 6) {
      tj.PDGCode = 0;
      return;
    }
    if (tj.Strategy[kStiffMu]) {
      tj.PDGCode = 13;
      return;
    }
    if (tj.Strategy[kStiffEl]) {
      tj.PDGCode = 111;
      return;
    }
    if (tcc.muonTag[0] <= 0) return;
    // Special handling of very long straight trajectories, e.g. uB cosmic rays
    bool isAMuon = (npwc > (unsigned short)tcc.muonTag[0] && tj.MCSMom > tcc.muonTag[1]);
    // anything really really long must be a muon
    if (npwc > 500) isAMuon = true;
    if (isAMuon) tj.PDGCode = 13;

  } // SetPDGCode

  ////////////////////////////////////////////////
  bool
  AnalyzeHits()
  {
    // Find the average hit rms by analyzing the full hit collection. This
    // only needs to be done once per job.

    if ((*evt.allHits).empty()) return true;
    // no sense re-calculating it if it's been done
    if (evt.aveHitRMSValid) return true;

    unsigned short cstat = (*evt.allHits)[0].WireID().Cryostat;
    unsigned short tpc = (*evt.allHits)[0].WireID().TPC;

    unsigned short nplanes = tcc.geom->Nplanes(tpc, cstat);
    evt.aveHitRMS.resize(nplanes);
    std::vector<float> cnt(nplanes, 0);
    for (unsigned short iht = 0; iht < (*evt.allHits).size(); ++iht) {
      auto& hit = (*evt.allHits)[iht];
      unsigned short plane = hit.WireID().Plane;
      if (plane > nplanes - 1) return false;
      if (cnt[plane] > 200) continue;
      // require multiplicity one
      if (hit.Multiplicity() != 1) continue;
      // not-crazy Chisq/DOF
      if (hit.GoodnessOfFit() < 0 || hit.GoodnessOfFit() > 500) continue;
      // don't let a lot of runt hits screw up the calculation
      if (hit.PeakAmplitude() < 1) continue;
      evt.aveHitRMS[plane] += hit.RMS();
      ++cnt[plane];
      // quit if enough hits are found
      bool allDone = true;
      for (unsigned short plane = 0; plane < nplanes; ++plane)
        if (cnt[plane] < 200) allDone = false;
      if (allDone) break;
    } // iht

    // assume there are enough hits in each plane
    evt.aveHitRMSValid = true;
    for (unsigned short plane = 0; plane < nplanes; ++plane) {
      if (cnt[plane] > 4) { evt.aveHitRMS[plane] /= cnt[plane]; }
      else {
        evt.aveHitRMS[plane] = 10;
        evt.aveHitRMSValid = false;
      } // cnt too low
    }   // plane

    if (tcc.modes[kModeDebug]) {
      std::cout << "Analyze hits aveHitRMS";
      std::cout << std::fixed << std::setprecision(1);
      for (auto rms : evt.aveHitRMS)
        std::cout << " " << rms;
      std::cout << " aveHitRMSValid? " << evt.aveHitRMSValid << "\n";
    }

    return true;
  } // Analyze hits

  ////////////////////////////////////////////////
  bool
  LongPulseHit(const recob::Hit& hit)
  {
    // return true if the hit is in a long pulse indicating that it's position
    // and charge are not well known
    return ((hit.GoodnessOfFit() < 0 || hit.GoodnessOfFit() > 50) && hit.Multiplicity() > 5);
  }

  ////////////////////////////////////////////////
  void
  FillWireHitRange(geo::TPCID inTPCID)
  {
    // Defines the local vector of dead wires and the low-high range of hits in each wire in
    // the TPCID in TCEvent. Note that there is no requirement that the allHits collection is sorted. Care should
    // be taken when looping over hits using this range - see SignalAtTp

    // see if this function was called in the current TPCID. There is nothing that needs to
    // be done if that is the case
    if (inTPCID == evt.TPCID) return;

    evt.TPCID = inTPCID;
    unsigned short nplanes = tcc.geom->Nplanes(inTPCID);
    unsigned int cstat = inTPCID.Cryostat;
    unsigned int tpc = inTPCID.TPC;
    if (tcc.useChannelStatus) {
      lariov::ChannelStatusProvider const& channelStatus =
        art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
      evt.goodWire.resize(nplanes);
      for (unsigned short pln = 0; pln < nplanes; ++pln) {
        unsigned int nwires = tcc.geom->Nwires(pln, tpc, cstat);
        // set all wires dead
        evt.goodWire[pln].resize(nwires, false);
        for (unsigned int wire = 0; wire < nwires; ++wire) {
          raw::ChannelID_t chan =
            tcc.geom->PlaneWireToChannel((int)pln, (int)wire, (int)tpc, (int)cstat);
          evt.goodWire[pln][wire] = channelStatus.IsGood(chan);
        } // wire
      }   // pln
    }
    else {
      // resize and set every channel good
      evt.goodWire.resize(nplanes);
      for (unsigned short pln = 0; pln < nplanes; ++pln) {
        unsigned int nwires = tcc.geom->Nwires(pln, tpc, cstat);
        evt.goodWire[pln].resize(nwires, true);
      } // pln
    }   // don't use channelStatus

    // there is no need to define evt.wireHitRange if the hit collection is not sliced. The function
    // SignalAtTP will then use the (smaller) slc.WireHitRange instead of evt.wireHitRange
    if (!evt.expectSlicedHits) return;

    // define the size of evt.wireHitRange
    evt.wireHitRange.resize(nplanes);
    for (unsigned short pln = 0; pln < nplanes; ++pln) {
      unsigned int nwires = tcc.geom->Nwires(pln, tpc, cstat);
      evt.wireHitRange[pln].resize(nwires);
      for (unsigned int wire = 0; wire < nwires; ++wire)
        evt.wireHitRange[pln][wire] = {UINT_MAX, UINT_MAX};
    } // pln

    // next define the wireHitRange values. Make one loop through the allHits collection
    unsigned int nBadWireFix = 0;
    for (unsigned int iht = 0; iht < (*evt.allHits).size(); ++iht) {
      auto& hit = (*evt.allHits)[iht];
      auto wid = hit.WireID();
      if (wid.Cryostat != cstat) continue;
      if (wid.TPC != tpc) continue;
      unsigned short pln = wid.Plane;
      unsigned int wire = wid.Wire;
      // Check the goodWire status and correct it if it's wrong
      if (!evt.goodWire[pln][wire]) {
        evt.goodWire[pln][wire] = true;
        ++nBadWireFix;
      } // not goodWire
      if (evt.wireHitRange[pln][wire].first == UINT_MAX) evt.wireHitRange[pln][wire].first = iht;
      evt.wireHitRange[pln][wire].second = iht;
    } // iht
    if (nBadWireFix > 0 && tcc.modes[kModeDebug]) {
      std::cout << "FillWireHitRange found hits on " << nBadWireFix
                << " wires that were declared not-good by the ChannelStatus service. Fixed it...\n";
    }
  } // FillWireHitRange

  ////////////////////////////////////////////////
  bool
  FillWireHitRange(detinfo::DetectorClocksData const& clockData,
                   detinfo::DetectorPropertiesData const& detProp,
                   TCSlice& slc)
  {
    // fills the WireHitRange vector. Slightly modified version of the one in ClusterCrawlerAlg.
    // Returns false if there was a serious error

    // determine the number of planes
    unsigned int cstat = slc.TPCID.Cryostat;
    unsigned int tpc = slc.TPCID.TPC;
    unsigned short nplanes = tcc.geom->Nplanes(tpc, cstat);
    slc.nPlanes = nplanes;
    if (nplanes > 3) return false;

    // Y,Z limits of the detector
    double local[3] = {0., 0., 0.};
    double world[3] = {0., 0., 0.};
    const geo::TPCGeo& thetpc = tcc.geom->TPC(tpc, cstat);
    thetpc.LocalToWorld(local, world);
    // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
    slc.xLo = world[0] - tcc.geom->DetHalfWidth(tpc, cstat) + 1;
    slc.xHi = world[0] + tcc.geom->DetHalfWidth(tpc, cstat) - 1;
    slc.yLo = world[1] - tcc.geom->DetHalfHeight(tpc, cstat) + 1;
    slc.yHi = world[1] + tcc.geom->DetHalfHeight(tpc, cstat) - 1;
    slc.zLo = world[2] - tcc.geom->DetLength(tpc, cstat) / 2 + 1;
    slc.zHi = world[2] + tcc.geom->DetLength(tpc, cstat) / 2 - 1;

    // initialize everything
    slc.wireHitRange.resize(nplanes);
    slc.firstWire.resize(nplanes);
    slc.lastWire.resize(nplanes);
    slc.nWires.resize(nplanes);
    tcc.maxPos0.resize(nplanes);
    tcc.maxPos1.resize(nplanes);
    evt.aveHitRMS.resize(nplanes, nplanes);

    std::pair<unsigned int, unsigned int> flag;
    flag.first = UINT_MAX;
    flag.second = UINT_MAX;

    // Calculate tcc.unitsPerTick, the scale factor to convert a tick into
    // Wire Spacing Equivalent (WSE) units where the wire spacing in this plane = 1.
    // Strictly speaking this factor should be calculated for each plane to handle the
    // case where the wire spacing is different in each plane. Deal with this later if
    // the approximation used here fails.

    raw::ChannelID_t channel = tcc.geom->PlaneWireToChannel(0, 0, (int)tpc, (int)cstat);
    tcc.wirePitch = tcc.geom->WirePitch(tcc.geom->View(channel));
    float tickToDist = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
    tickToDist *= 1.e-3 * sampling_rate(clockData); // 1e-3 is conversion of 1/us to 1/ns
    tcc.unitsPerTick = tickToDist / tcc.wirePitch;
    for (unsigned short plane = 0; plane < nplanes; ++plane) {
      slc.firstWire[plane] = UINT_MAX;
      slc.lastWire[plane] = 0;
      slc.nWires[plane] = tcc.geom->Nwires(plane, tpc, cstat);
      slc.wireHitRange[plane].resize(slc.nWires[plane], flag);
      tcc.maxPos0[plane] = (float)slc.nWires[plane] - 0.5;
      tcc.maxPos1[plane] = (float)detProp.NumberTimeSamples() * tcc.unitsPerTick;
    }

    unsigned int lastWire = 0, lastPlane = 0;
    for (unsigned int iht = 0; iht < slc.slHits.size(); ++iht) {
      unsigned int ahi = slc.slHits[iht].allHitsIndex;
      if (ahi > (*evt.allHits).size() - 1) return false;
      auto& hit = (*evt.allHits)[ahi];
      if (hit.WireID().Cryostat != cstat) continue;
      if (hit.WireID().TPC != tpc) continue;
      unsigned short plane = hit.WireID().Plane;
      unsigned int wire = hit.WireID().Wire;
      if (wire > slc.nWires[plane] - 1) {
        mf::LogWarning("TC") << "FillWireHitRange: Invalid wire number " << wire << " > "
                             << slc.nWires[plane] - 1 << " in plane " << plane << " Quitting";
        return false;
      } // too large wire number
      if (plane == lastPlane && wire < lastWire) {
        mf::LogWarning("TC")
          << "FillWireHitRange: Hits are not in increasing wire order. Quitting ";
        return false;
      } // hits out of order
      lastWire = wire;
      lastPlane = plane;
      if (slc.firstWire[plane] == UINT_MAX) slc.firstWire[plane] = wire;
      if (slc.wireHitRange[plane][wire].first == UINT_MAX)
        slc.wireHitRange[plane][wire].first = iht;
      slc.wireHitRange[plane][wire].second = iht;
      slc.lastWire[plane] = wire + 1;
    } // iht
    // check
    unsigned int slhitsSize = slc.slHits.size();
    for (unsigned short plane = 0; plane < nplanes; ++plane) {
      for (unsigned int wire = slc.firstWire[plane]; wire < slc.lastWire[plane]; ++wire) {
        if (slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
        if (slc.wireHitRange[plane][wire].first > slhitsSize - 1 &&
            slc.wireHitRange[plane][wire].second > slhitsSize)
          return false;
      } // wire
    }   // plane

    // Find the average multiplicity 1 hit RMS and calculate the expected max RMS for each range
    if (tcc.modes[kModeDebug] && (int)tpc == debug.TPC) {
      // Note that this function is called before the slice is pushed into slices so the index
      // isn't decremented by 1
      std::cout << "Slice ID/Index " << slc.ID << "/" << slices.size() << " tpc " << tpc
                << " tcc.unitsPerTick " << std::setprecision(3) << tcc.unitsPerTick;
      std::cout << " Active volume (";
      std::cout << std::fixed << std::setprecision(1) << slc.xLo << " < X < " << slc.xHi << ") (";
      std::cout << std::fixed << std::setprecision(1) << slc.yLo << " < Y < " << slc.yHi << ") (";
      std::cout << std::fixed << std::setprecision(1) << slc.zLo << " < Z < " << slc.zHi << ")\n";
    }

    return true;

  } // FillWireHitRange

  ////////////////////////////////////////////////
  bool
  WireHitRangeOK(TCSlice& slc, const CTP_t& inCTP)
  {
    // returns true if the passed CTP code is consistent with the CT code of the WireHitRangeVector
    geo::PlaneID planeID = DecodeCTP(inCTP);
    if (planeID.Cryostat != slc.TPCID.Cryostat) return false;
    if (planeID.TPC != slc.TPCID.TPC) return false;
    return true;
  }

  ////////////////////////////////////////////////
  bool
  MergeAndStore(TCSlice& slc, std::vector<int> tidList)
  {
    // merge a list of trajectories (assumed to be mostly junk Tjs) and
    // make a new junk Tj
    std::vector<unsigned int> jtHits;
    for(auto tid : tidList) {
      if(tid <= 0) return false;
      unsigned int itj = tid - 1;
      auto& tj = slc.tjs[itj];
      auto tHits = PutTrajHitsInVector(tj, kUsedHits);
      MakeTrajectoryObsolete(slc, itj);
      jtHits.insert(jtHits.end(), tHits.begin(), tHits.end());
    } // tid
    return MakeJunkTraj(slc, jtHits);
  } // MergeAndStore

  ////////////////////////////////////////////////
  bool
  MergeAndStore(TCSlice& slc, unsigned int itj1, unsigned int itj2, bool doPrt)
  {
    // Merge the two trajectories in allTraj and store them. Returns true if it was successfull.
    // Merging is done between the end (end = 1) of tj1 and the beginning (end = 0) of tj2. This function preserves the
    // AlgMod state of itj1.
    // The itj1 -> itj2 merge order is reversed if end1 of itj2 is closer to end0 of itj1

    if (itj1 > slc.tjs.size() - 1) return false;
    if (itj2 > slc.tjs.size() - 1) return false;
    if (slc.tjs[itj1].AlgMod[kKilled] || slc.tjs[itj2].AlgMod[kKilled]) return false;
    if (slc.tjs[itj1].AlgMod[kHaloTj] || slc.tjs[itj2].AlgMod[kHaloTj]) return false;

    // Ensure that the order of 3D-matched Tjs is consistent with the convention that
    unsigned short pfp1 = GetPFPIndex(slc, slc.tjs[itj1].ID);
    unsigned short pfp2 = GetPFPIndex(slc, slc.tjs[itj2].ID);
    if (pfp1 != USHRT_MAX || pfp2 != USHRT_MAX) {
      if (pfp1 != USHRT_MAX && pfp2 != USHRT_MAX) return false;
      // Swap so that the order of tj1 is preserved. Tj2 may be reversed to be consistent
      if (pfp1 == USHRT_MAX) std::swap(itj1, itj2);
    } // one or both used in a PFParticle

    // make copies so they can be trimmed as needed
    Trajectory tj1 = slc.tjs[itj1];
    Trajectory tj2 = slc.tjs[itj2];

    // ensure that these are in the same step order
    if (tj2.StepDir != tj1.StepDir) ReverseTraj(slc, tj2);

    Point2_t tp1e0 = tj1.Pts[tj1.EndPt[0]].Pos;
    Point2_t tp1e1 = tj1.Pts[tj1.EndPt[1]].Pos;
    Point2_t tp2e0 = tj2.Pts[tj2.EndPt[0]].Pos;
    Point2_t tp2e1 = tj2.Pts[tj2.EndPt[1]].Pos;

    if (doPrt) {
      mf::LogVerbatim("TC") << "MergeAndStore: T" << tj1.ID << " and T" << tj2.ID
                            << " at merge points " << PrintPos(slc, tp1e1) << " "
                            << PrintPos(slc, tp2e0);
    }

    // swap the order so that abs(tj1end1 - tj2end0) is less than abs(tj2end1 - tj1end0)
    if (PosSep2(tp1e1, tp2e0) > PosSep2(tp2e1, tp1e0)) {
      std::swap(tj1, tj2);
      std::swap(tp1e0, tp2e0);
      std::swap(tp1e1, tp2e1);
      if (doPrt)
        mf::LogVerbatim("TC") << " swapped the order. Merge points " << PrintPos(slc, tp1e1) << " "
                              << PrintPos(slc, tp2e0);
    }

    // Here is what we are looking for, where - indicates a TP with charge.
    // Note that this graphic is in the stepping direction (+1 = +wire direction)
    // tj1:  0------------1
    // tj2:                  0-----------1
    // Another possibility with overlap
    // tj1:  0-------------1
    // tj2:               0--------------1

    if (tj1.StepDir > 1) {
      // Not allowed
      // tj1:  0---------------------------1
      // tj2:                  0------1
      if (tp2e0[0] > tp1e0[0] && tp2e1[0] < tp1e1[0]) return false;
      /// Not allowed
      // tj1:                  0------1
      // tj2:  0---------------------------1
      if (tp1e0[0] > tp2e0[0] && tp1e1[0] < tp2e1[0]) return false;
    }
    else {
      // same as above but with ends reversed
      if (tp2e1[0] > tp1e1[0] && tp2e0[0] < tp1e0[0]) return false;
      if (tp1e1[0] > tp2e1[0] && tp1e0[0] < tp2e0[0]) return false;
    }

    if (tj1.VtxID[1] > 0 && tj2.VtxID[0] == tj1.VtxID[1]) {
      auto& vx = slc.vtxs[tj1.VtxID[1] - 1];
      if (!MakeVertexObsolete("MAS", slc, vx, false)) {
        if (doPrt)
          mf::LogVerbatim("TC") << "MergeAndStore: Found a good vertex between Tjs " << tj1.VtxID[1]
                                << " No merging";
        return false;
      }
    }

    // remove any points at the end of tj1 that don't have used hits
    tj1.Pts.resize(tj1.EndPt[1] + 1);

    // determine if they overlap by finding the point on tj2 that is closest
    // to the end point of tj1.
    TrajPoint& endtj1TP = tj1.Pts[tj1.EndPt[1]];
    // Set minSep large so that dead wire regions are accounted for
    float minSep = 1000;
    unsigned short tj2ClosePt = 0;
    // Note that TrajPointTrajDOCA only considers TPs that have charge
    TrajPointTrajDOCA(slc, endtj1TP, tj2, tj2ClosePt, minSep);
    if (doPrt)
      mf::LogVerbatim("TC") << " Merge point tj1 " << PrintPos(slc, endtj1TP) << " tj2ClosePt "
                            << tj2ClosePt << " Pos " << PrintPos(slc, tj2.Pts[tj2ClosePt]);
    // check for full overlap
    if (tj2ClosePt > tj2.EndPt[1]) return false;

    // The approach is to append tj2 to tj1, store tj1 as a new trajectory,
    // and re-assign all hits to the new trajectory

    // First ensure that any hit will appear only once in the merged trajectory in the overlap region
    // whether it is used or unused. The point on tj2 where the merge will begin, tj2ClosePt, will be
    // increased until this condition is met.
    // Make a temporary vector of tj1 hits in the end points for simpler searching
    std::vector<unsigned int> tj1Hits;
    for (unsigned short ii = 0; ii < tj1.Pts.size(); ++ii) {
      // only go back a few points in tj1
      if (ii > 10) break;
      unsigned short ipt = tj1.Pts.size() - 1 - ii;
      tj1Hits.insert(tj1Hits.end(), tj1.Pts[ipt].Hits.begin(), tj1.Pts[ipt].Hits.end());
      if (ipt == 0) break;
    } // ii

    bool bumpedPt = true;
    while (bumpedPt) {
      bumpedPt = false;
      for (unsigned short ii = 0; ii < tj2.Pts[tj2ClosePt].Hits.size(); ++ii) {
        unsigned int iht = tj2.Pts[tj2ClosePt].Hits[ii];
        if (std::find(tj1Hits.begin(), tj1Hits.end(), iht) != tj1Hits.end()) bumpedPt = true;
      } // ii
      if (bumpedPt && tj2ClosePt < tj2.EndPt[1]) { ++tj2ClosePt; }
      else {
        break;
      }
    } // bumpedPt
    if (doPrt) mf::LogVerbatim("TC") << " revised tj2ClosePt " << tj2ClosePt;
    // append tj2 hits to tj1

    tj1.Pts.insert(tj1.Pts.end(), tj2.Pts.begin() + tj2ClosePt, tj2.Pts.end());
    // re-define the end points
    SetEndPoints(tj1);
    tj1.EndFlag[1] = tj2.EndFlag[1];

    // A more exhaustive check that hits only appear once
    if (HasDuplicateHits(slc, tj1, doPrt)) return false;
    if (tj2.VtxID[1] > 0) {
      // move the end vertex of tj2 to the end of tj1
      tj1.VtxID[1] = tj2.VtxID[1];
    }
    // Transfer some of the AlgMod bits
    if (tj2.AlgMod[kMichel]) tj1.AlgMod[kMichel] = true;
    if (tj2.AlgMod[kDeltaRay]) {
      tj1.AlgMod[kDeltaRay] = true;
      tj1.ParentID = tj2.ParentID;
    }
    // keep track of the IDs before they are clobbered
    int tj1ID = tj1.ID;
    int tj2ID = tj2.ID;
    // kill the original trajectories
    MakeTrajectoryObsolete(slc, itj1);
    MakeTrajectoryObsolete(slc, itj2);
    // Do this so that StoreTraj keeps the correct WorkID (of itj1)
    tj1.ID = tj1.WorkID;
    SetPDGCode(slc, tj1);
    tj1.NeedsUpdate = true;
    UpdateTjChgProperties("MAS", slc, tj1, doPrt);
    ChkStop(slc, tj1);
    if (!StoreTraj(slc, tj1)) return false;
    int newTjID = slc.tjs.size();
    // Use the ParentID to trace which new Tj is superseding the merged ones
    tj1.ParentID = newTjID;
    tj2.ParentID = newTjID;
    if (doPrt) mf::LogVerbatim("TC") << " MAS success. Created T" << newTjID;
    // Transfer the ParentIDs of any other Tjs that refer to Tj1 and Tj2 to the new Tj
    for (auto& tj : slc.tjs)
      if (tj.ParentID == tj1ID || tj.ParentID == tj2ID) tj.ParentID = newTjID;
    // try to attach it to a vertex
    AttachAnyVertexToTraj(slc, newTjID, doPrt);
    return true;
  } // MergeAndStore

  ////////////////////////////////////////////////
  std::vector<int>
  GetAssns(TCSlice& slc, std::string type1Name, int id, std::string type2Name)
  {
    // returns a list of IDs of objects (slc, vertices, pfps, etc) with type1Name that are in slc with
    // type2Name. This is intended to be a general purpose replacement for specific functions like GetVtxTjIDs, etc

    std::vector<int> tmp;
    if (id <= 0) return tmp;
    unsigned int uid = id;

    if (type1Name == "T" && uid <= slc.tjs.size() && type2Name == "P") {
      // return a list of PFPs that have the tj in TjIDs, P -> T<ID>
      for (auto& pfp : slc.pfps) {
        if (pfp.ID <= 0) continue;
        if (std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), id) != pfp.TjIDs.end())
          tmp.push_back(pfp.ID);
      } // pf
      return tmp;
    } // P -> T

    if (type1Name == "2V" && uid <= slc.vtxs.size() && type2Name == "T") {
      // 2V -> T
      for (auto& tj : slc.tjs) {
        if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
        for (unsigned short end = 0; end < 2; ++end) {
          if (tj.VtxID[end] != id) continue;
          if (std::find(tmp.begin(), tmp.end(), tj.ID) == tmp.end()) tmp.push_back(tj.ID);
        } // end
      }   // tj
      return tmp;
    } // 2V -> T

    if (type1Name == "3V" && uid <= slc.vtx3s.size() && type2Name == "P") {
      for (auto& pfp : slc.pfps) {
        if (pfp.ID == 0) continue;
        for (unsigned short end = 0; end < 2; ++end) {
          if (pfp.Vx3ID[end] != id) continue;
          // encode the end with the ID
          if (std::find(tmp.begin(), tmp.end(), pfp.ID) == tmp.end()) tmp.push_back(pfp.ID);
        } // end
      }   // pfp
      return tmp;
    } // 3V -> P

    if (type1Name == "3V" && uid <= slc.vtx3s.size() && type2Name == "T") {
      // 3V -> T
      for (auto& tj : slc.tjs) {
        if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
        for (unsigned short end = 0; end < 2; ++end) {
          if (tj.VtxID[end] > 0 && tj.VtxID[end] <= slc.vtxs.size()) {
            auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
            if (vx2.Vx3ID != id) continue;
            if (std::find(tmp.begin(), tmp.end(), tj.ID) == tmp.end()) tmp.push_back(tj.ID);
          }
        } // end
      }   // tj
      return tmp;
    } // 3V -> T

    if (type1Name == "3V" && uid <= slc.vtx3s.size() && type2Name == "2V") {
      // 3V -> 2V
      for (auto& vx2 : slc.vtxs) {
        if (vx2.ID == 0) continue;
        if (vx2.Vx3ID == id) tmp.push_back(vx2.ID);
      } // vx2
      return tmp;
    } // 3V -> 2V

    return tmp;
  } // GetAssns

  ////////////////////////////////////////////////
  bool
  StartTraj(TCSlice& slc,
            Trajectory& tj,
            unsigned int fromhit,
            unsigned int tohit,
            unsigned short pass)
  {
    // Start a trajectory located at fromHit with direction pointing to toHit

    auto& fromHit = (*evt.allHits)[slc.slHits[fromhit].allHitsIndex];
    auto& toHit = (*evt.allHits)[slc.slHits[tohit].allHitsIndex];
    float fromWire = fromHit.WireID().Wire;
    float fromTick = fromHit.PeakTime();
    float toWire = toHit.WireID().Wire;
    float toTick = toHit.PeakTime();
    CTP_t tCTP = EncodeCTP(fromHit.WireID());
    bool success = StartTraj(slc, tj, fromWire, fromTick, toWire, toTick, tCTP, pass);
    if (!success) return false;
    // turn on debugging using the WorkID?
    if (tcc.modes[kModeDebug] && !tcc.dbgStp && !tcc.dbgDump && tcc.dbgSlc && tj.ID == debug.WorkID)
      tcc.dbgStp = true;
    if (tcc.dbgStp) {
      auto& tp = tj.Pts[0];
      mf::LogVerbatim("TC") << "StartTraj T" << tj.ID << " from " << (int)fromWire << ":"
                            << (int)fromTick << " -> " << (int)toWire << ":" << (int)toTick
                            << " StepDir " << tj.StepDir << " dir " << tp.Dir[0] << " " << tp.Dir[1]
                            << " ang " << tp.Ang << " AngleCode " << tp.AngleCode << " angErr "
                            << tp.AngErr << " ExpectedHitsRMS " << ExpectedHitsRMS(slc, tp);
    } // tcc.dbgStp
    return true;
  } // StartTraj

  ////////////////////////////////////////////////
  bool
  StartTraj(TCSlice& slc,
            Trajectory& tj,
            float fromWire,
            float fromTick,
            float toWire,
            float toTick,
            CTP_t& tCTP,
            unsigned short pass)
  {
    // Start a simple (seed) trajectory going from (fromWire, toTick) to (toWire, toTick).

    // decrement the work ID so we can use it for debugging problems
    --evt.WorkID;
    if (evt.WorkID == INT_MIN) evt.WorkID = -1;
    tj.ID = evt.WorkID;
    tj.Pass = pass;
    // Assume we are stepping in the positive WSE units direction
    short stepdir = 1;
    int fWire = std::nearbyint(fromWire);
    int tWire = std::nearbyint(toWire);
    if (tWire < fWire) { stepdir = -1; }
    else if (tWire == fWire) {
      // on the same wire
      if (toTick < fromTick) stepdir = -1;
    }
    tj.StepDir = stepdir;
    tj.CTP = tCTP;
    tj.ParentID = -1;
    tj.Strategy.reset();
    tj.Strategy[kNormal] = true;

    // create a trajectory point
    TrajPoint tp;
    if (!MakeBareTrajPoint(slc, fromWire, fromTick, toWire, toTick, tCTP, tp)) return false;
    SetAngleCode(tp);
    tp.AngErr = 0.1;
    tj.Pts.push_back(tp);
    // turn on debugging using the WorkID?
    if (tcc.modes[kModeDebug] && !tcc.dbgStp && !tcc.dbgDump && tcc.dbgSlc && tj.ID == debug.WorkID)
      tcc.dbgStp = true;
    if (tcc.dbgStp) {
      auto& tp = tj.Pts[0];
      mf::LogVerbatim("TC") << "StartTraj T" << tj.ID << " from " << (int)fromWire << ":"
                            << (int)fromTick << " -> " << (int)toWire << ":" << (int)toTick
                            << " StepDir " << tj.StepDir << " dir " << tp.Dir[0] << " " << tp.Dir[1]
                            << " ang " << tp.Ang << " AngleCode " << tp.AngleCode << " angErr "
                            << tp.AngErr << " ExpectedHitsRMS " << ExpectedHitsRMS(slc, tp);
    } // tcc.dbgStp
    return true;

  } // StartTraj

  ////////////////////////////////////////////////
  std::pair<unsigned short, unsigned short>
  GetSliceIndex(std::string typeName, int uID)
  {
    // returns the slice index and product index of a data product having typeName and unique ID uID
    for (unsigned short isl = 0; isl < slices.size(); ++isl) {
      auto& slc = slices[isl];
      if (typeName == "T") {
        for (unsigned short indx = 0; indx < slc.tjs.size(); ++indx) {
          if (slc.tjs[indx].UID == uID) { return std::make_pair(isl, indx); }
        }
      } // T
      if (typeName == "P") {
        for (unsigned short indx = 0; indx < slc.pfps.size(); ++indx) {
          if (slc.pfps[indx].UID == uID) { return std::make_pair(isl, indx); }
        }
      } // P
      if (typeName == "2V") {
        for (unsigned short indx = 0; indx < slc.vtxs.size(); ++indx) {
          if (slc.vtxs[indx].UID == uID) { return std::make_pair(isl, indx); }
        }
      } // 2V
      if (typeName == "3V") {
        for (unsigned short indx = 0; indx < slc.vtx3s.size(); ++indx) {
          if (slc.vtx3s[indx].UID == uID) { return std::make_pair(isl, indx); }
        }
      } // 3V
    }   // isl
    return std::make_pair(USHRT_MAX, USHRT_MAX);
  } // GetSliceIndex

  ////////////////////////////////////////////////
  bool
  Fit2D(short mode,
        Point2_t inPt,
        float& inPtErr,
        Vector2_t& outVec,
        Vector2_t& outVecErr,
        float& chiDOF)
  {
    // Fit points to a 2D line.
    // Mode = 0: Initialize
    // Mode = 1: Accumulate
    // Mode = 2: Accumulate and store to calculate chiDOF
    // Mode = -1: Fit and put results in outVec and chiDOF

    static double sum, sumx, sumy, sumx2, sumy2, sumxy;
    static unsigned short cnt;
    static std::vector<Point2_t> fitPts;
    static std::vector<double> fitWghts;

    if (mode == 0) {
      // initialize
      cnt = 0;
      sum = 0.;
      sumx = 0.;
      sumy = 0.;
      sumx2 = 0.;
      sumy2 = 0.;
      sumxy = 0;
      fitPts.resize(0);
      fitWghts.resize(0);
      return true;
    } // mode == 0

    if (mode > 0) {
      if (inPtErr <= 0.) return false;
      ++cnt;
      double wght = 1 / (inPtErr * inPtErr);
      sum += wght;
      sumx += wght * inPt[0];
      sumx2 += wght * inPt[0] * inPt[0];
      sumy += wght * inPt[1];
      sumy2 += wght * inPt[1] * inPt[1];
      sumxy += wght * inPt[0] * inPt[1];
      if (mode == 1) return true;
      fitPts.push_back(inPt);
      fitWghts.push_back(wght);
      return true;
    } // Accumulate

    if (cnt < 2) return false;
    // do the fit
    double delta = sum * sumx2 - sumx * sumx;
    if (delta == 0.) return false;
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    double B = (sumxy * sum - sumx * sumy) / delta;
    outVec[0] = A;
    outVec[1] = B;
    chiDOF = 0;
    if (cnt == 2 || fitPts.empty()) return true;

    // calculate errors and chiDOF
    if (fitPts.size() != cnt) return false;
    double ndof = cnt - 2;
    double varnce =
      (sumy2 + A * A * sum + B * B * sumx2 - 2 * (A * sumy + B * sumxy - A * B * sumx)) / ndof;
    if (varnce > 0.) {
      outVecErr[0] = sqrt(varnce * sumx2 / delta);
      outVecErr[1] = sqrt(varnce * sum / delta);
    }
    else {
      outVecErr[0] = 0.;
      outVecErr[1] = 0.;
    }
    sum = 0.;
    // calculate chisq
    for (unsigned short ii = 0; ii < fitPts.size(); ++ii) {
      double arg = fitPts[ii][1] - A - B * fitPts[ii][0];
      sum += fitWghts[ii] * arg * arg;
    }
    chiDOF = sum / ndof;
    fitPts.resize(0);
    fitWghts.resize(0);
    return true;

  } // Fit2D

} // namespace tca
