#include "larreco/RecoAlg/TCAlg/Utils.h"

#include "larsim/MCCheater/ParticleInventoryService.h" // for printing
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
#include "larreco/RecoAlg/TCAlg/TCShower.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h" // for tcc
#include "nusimdata/SimulationBase/MCParticle.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

namespace {
  struct SortEntry {
    unsigned int index;
    float val;
  };

  bool
  valDecreasing(SortEntry c1, SortEntry c2)
  {
    return (c1.val > c2.val);
  }
}

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
  void
  DefineTjParents(TCSlice& slc, bool prt)
  {
    /*
    This function sets the ParentUID of Tjs in this tpcid to create a hierarchy. The highest Score
    3D vertex in a chain of Tjs and vertices is declared the primary vertex; vx3.Primary = true. Tjs directly attached
    to that vertex are declared Primary trajectories with ParentUID = 0. All other Tjs in the chain have ParentUID
    set to the next upstream Tj to which it is attached by a vertex. In the graphical description below, V1 and V4 are
    2D vertices that are matched to a high-score 3D vertex. The V1 Score is greater than the V2 Score and V3 Score.
    V1 and V4 are declared to be primary vertices. T1, T2, T6 and T7 are declared to be primary Tjs

      V1 - T1 - V2 - T3          V4 - T6         / T8
         \                          \           /
           T2 - V3 - T4               T7
                   \
                     T5

    This is represented as follows. The NeutrinoPrimaryTjID is defined by a function.
     Tj   ParentUID   NeutrinoPrimaryTjID
     -----------------------------------
     T1      0          T1
     T2      0          T2
     T3     T1          T2
     T4     T2          T2
     T5     T2          T2
     T6      0          -1
     T7      0          -1
     T8     -1          -1
*/

    // don't do anything if this is test beam data
    if (tcc.modes[kTestBeam]) return;

    // clear old information
    for (auto& tj : slc.tjs) {
      if (tj.AlgMod[kKilled]) continue;
      // ignore delta rays
      if (tj.AlgMod[kDeltaRay] || tj.AlgMod[kHaloTj]) continue;
      tj.ParentID = 0;
    } // tj

    // sort vertice by decreasing score
    std::vector<int> temp;
    for (auto& vx3 : slc.vtx3s) {
      if (vx3.ID == 0) continue;
      // clear the Primary flag while we are here
      vx3.Primary = false;
      temp.push_back(vx3.ID);
    } // vx3
    if (temp.empty()) return;

    // Make a master list of all Tjs that are attached to these vertices
    std::vector<int> masterlist;
    for (auto vx3id : temp) {
      auto& vx3 = slc.vtx3s[vx3id - 1];
      float score;
      auto tjlist = GetVtxTjIDs(slc, vx3, score);
      for (auto tjid : tjlist) {
        auto& tj = slc.tjs[tjid - 1];
        if (tj.ParentID != 0) tj.ParentID = 0;
        if (std::find(masterlist.begin(), masterlist.end(), tjid) == masterlist.end())
          masterlist.push_back(tjid);
      } // tjid
    }   // vxid
    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "DTP: masterlist Tjs";
      for (auto tjid : masterlist)
        myprt << " " << tjid;
    }

    // Do the sort
    std::vector<SortEntry> sortVec(temp.size());
    for (unsigned short indx = 0; indx < temp.size(); ++indx) {
      auto& vx3 = slc.vtx3s[temp[indx] - 1];
      sortVec[indx].index = indx;
      sortVec[indx].val = vx3.Score;
    } // indx
    if (sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasing);
    // put them into order
    auto vlist = temp;
    for (unsigned short indx = 0; indx < temp.size(); ++indx)
      vlist[indx] = temp[sortVec[indx].index];

    // make a neutrino PFParticle to associate with the highest score vertex if it is high enough
    if (tcc.match3DCuts[0] > 0) {
      auto& vx3 = slc.vtx3s[vlist[0] - 1];
      if (vx3.Score > tcc.vtx2DCuts[7]) {
        auto neutrinoPFP = CreatePFP(slc);
        // call it the neutrino vertex
        vx3.Neutrino = true;
        // put the vertex at the end of the neutrino
        auto& sf = neutrinoPFP.SectionFits[0];
        sf.Pos[0] = vx3.X;
        sf.Pos[1] = vx3.Y;
        sf.Pos[2] = vx3.Z;
        sf.Dir[2] = 1;
        // This may be set to 12 later on if a primary shower is reconstructed
        neutrinoPFP.PDGCode = 14;
        neutrinoPFP.Vx3ID[1] = vx3.ID;
        neutrinoPFP.Vx3ID[0] = vx3.ID;
        neutrinoPFP.Flags[kNeedsUpdate] = false;
        // the rest of this will be defined later
        if (!StorePFP(slc, neutrinoPFP)) return;
      }
    } // User wants to make PFParticles
    // a temp vector to ensure that we only consider a vertex once
    std::vector<bool> lookedAt3(slc.vtx3s.size() + 1, false);
    std::vector<bool> lookedAt2(slc.vtxs.size() + 1, false);
    // vector of parent-daughter pairs
    std::vector<std::pair<int, int>> pardtr;
    // Start with the highest score vertex
    for (unsigned short indx = 0; indx < vlist.size(); ++indx) {
      auto& vx3 = slc.vtx3s[vlist[indx] - 1];
      if (lookedAt3[vx3.ID]) continue;
      vx3.Primary = true;
      lookedAt3[vx3.ID] = true;
      // make a list of Tjs attached to this vertex
      float score;
      auto primTjList = GetVtxTjIDs(slc, vx3, score);
      if (primTjList.empty()) continue;
      pardtr.clear();
      for (auto primTjID : primTjList) {
        auto& primTj = slc.tjs[primTjID - 1];
        // This isn't a primary tj if the parent ID isn't -1
        if (primTj.ParentID != -1) continue;
        if (prt) mf::LogVerbatim("TC") << "Vx3 " << vx3.ID << " Primary tj " << primTj.ID;
        // declare this a primary tj
        primTj.ParentID = 0;
        // look for daughter tjs = those that are attached to a 2D vertex
        // at the other end
        for (unsigned short end = 0; end < 2; ++end) {
          if (primTj.VtxID[end] == 0) continue;
          auto& vx2 = slc.vtxs[primTj.VtxID[end] - 1];
          if (vx2.Vx3ID == vx3.ID) continue;
          // found a 2D vertex. Check for daughters
          auto dtrList = GetVtxTjIDs(slc, vx2);
          for (auto dtrID : dtrList) {
            // ignore the primary tj
            if (dtrID == primTjID) continue;
            auto& dtj = slc.tjs[dtrID - 1];
            if (dtj.ParentID != -1) continue;
            pardtr.push_back(std::make_pair(primTjID, dtrID));
            if (prt) mf::LogVerbatim("TC") << "  primTj " << primTjID << " dtrID " << dtrID;
          } // tjid
        }   // end
        // Ensure that end 0 of the trajectory is attached to the primary vertex
        for (unsigned short end = 0; end < 2; ++end) {
          if (primTj.VtxID[end] == 0) continue;
          auto& vx2 = slc.vtxs[primTj.VtxID[end] - 1];
          if (vx2.Vx3ID == vx3.ID && end != 0) ReverseTraj(slc, primTj);
        } // end
      }   // tjid
      if (pardtr.empty()) continue;
      if (prt) {
        mf::LogVerbatim myprt("TC");
        myprt << " par_dtr";
        for (auto pdtr : pardtr)
          myprt << " " << pdtr.first << "_" << pdtr.second;
      }
      // iterate through the parent - daughter stack, removing the last pair when a
      // ParentID is updated and adding pairs for new daughters
      for (unsigned short nit = 0; nit < 100; ++nit) {
        auto lastPair = pardtr[pardtr.size() - 1];
        auto& dtj = slc.tjs[lastPair.second - 1];
        dtj.ParentID = lastPair.first;
        // reverse the daughter trajectory if necessary so that end 0 is closest to the parent
        float doca = 100;
        unsigned short dpt = 0, ppt = 0;
        auto& ptj = slc.tjs[lastPair.first - 1];
        // find the point on the daughter tj that is closest to the parent
        TrajTrajDOCA(slc, dtj, ptj, dpt, ppt, doca);
        // reverse the daughter if the closest point is near end 1 of the daughter
        if (prt) mf::LogVerbatim("TC") << "Set parent " << ptj.ID << " dtr " << dtj.ID;
        // remove that entry
        pardtr.pop_back();
        // Add entries for new daughters
        for (unsigned short end = 0; end < 2; ++end) {
          if (dtj.VtxID[end] == 0) continue;
          auto& vx2 = slc.vtxs[dtj.VtxID[end] - 1];
          if (lookedAt2[vx2.ID]) continue;
          lookedAt2[vx2.ID] = true;
          auto tjlist = GetVtxTjIDs(slc, vx2);
          for (auto tjid : tjlist) {
            if (tjid == dtj.ID || tjid == ptj.ID) continue;
            pardtr.push_back(std::make_pair(dtj.ID, tjid));
            if (prt) {
              mf::LogVerbatim myprt("TC");
              myprt << " add par_dtr";
              for (auto pdtr : pardtr)
                myprt << " " << pdtr.first << "_" << pdtr.second;
            }
          }
        } // end
        if (pardtr.empty()) break;
      } // nit
    }   // indx
    // check the master list
    for (auto tjid : masterlist) {
      auto& tj = slc.tjs[tjid - 1];
      if (tj.ParentID < 0) tj.ParentID = tj.ID;
    } // tjid

  } // DefineTjParents

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
  bool
  MergeTjIntoPFP(TCSlice& slc, int mtjid, PFPStruct& pfp, bool prt)
  {
    // Tries to merge Tj with ID tjid into PFParticle pfp
    if (mtjid > (int)slc.tjs.size()) return false;
    auto& mtj = slc.tjs[mtjid - 1];
    // find the Tj in pfp.TjIDs which it should be merged with
    int otjid = 0;
    for (auto tjid : pfp.TjIDs) {
      auto& otj = slc.tjs[tjid - 1];
      if (otj.CTP == mtj.CTP) {
        otjid = tjid;
        break;
      }
    } // tjid
    if (otjid == 0) return false;
    if (MergeAndStore(slc, otjid - 1, mtjid - 1, prt)) {
      int newtjid = slc.tjs.size();
      if (prt)
        mf::LogVerbatim("TC") << "MergeTjIntoPFP: merged T" << otjid << " with T" << mtjid
                              << " -> T" << newtjid;
      std::replace(pfp.TjIDs.begin(), pfp.TjIDs.begin(), otjid, newtjid);
      return true;
    }
    else {
      if (prt)
        mf::LogVerbatim("TC") << "MergeTjIntoPFP: merge T" << otjid << " with T" << mtjid
                              << " failed ";
      return false;
    }
  } // MergeTjIntoPFP

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

    //    static const float twoPi = 2 * M_PI;

    if (originPt > tj.Pts.size() - 1) {
      mf::LogWarning("TC") << "FitTraj: Requesting fit of invalid TP " << originPt;
      return;
    }

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
    for (auto& tp : tj.Pts) {
      if (tp.Hits.size() > 16) return false;
    } // tp

    if (tj.NeedsUpdate) UpdateTjChgProperties("ST", slc, tj, false);

    // This shouldn't be necessary but do it anyway
    SetEndPoints(tj);

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
      if (slc.slHits[iht].InTraj == tj.ID) {
//        mf::LogWarning("TC") << "StoreTraj: Hit " << PrintHit(slc.slHits[iht])
//                             << " thinks it belongs to T" << tj.ID 
//                             << " but it isn't in the Tj... Fixed";
        slc.slHits[iht].InTraj = 0;
      }
    } // iht

    tj.WorkID = tj.ID;
    tj.ID = trID;
    // increment the global ID
    ++evt.globalT_UID;
    tj.UID = evt.globalT_UID;
    // Don't clobber the ParentID if it was defined by the calling function
    if (tj.ParentID == 0) tj.ParentID = trID;
    slc.tjs.push_back(tj);
    if (tcc.modes[kDebug] && tcc.dbgSlc && debug.Hit != UINT_MAX) {
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

  ////////////////////////////////////////////////
  bool
  InTrajOK(TCSlice& slc, std::string someText)
  {
    // Check slc.tjs -> InTraj associations

    unsigned short tID;
    unsigned int iht;
    unsigned short itj = 0;
    std::vector<unsigned int> tHits;
    std::vector<unsigned int> atHits;
    for (auto& tj : slc.tjs) {
      // ignore abandoned trajectories
      if (tj.AlgMod[kKilled]) continue;
      tID = tj.ID;
      tHits = PutTrajHitsInVector(tj, kUsedHits);
      if (tHits.size() < 2) continue;
      std::sort(tHits.begin(), tHits.end());
      atHits.clear();
      for (iht = 0; iht < slc.slHits.size(); ++iht) {
        if (slc.slHits[iht].InTraj == tID) atHits.push_back(iht);
      } // iht
      if (atHits.size() < 2) continue;
      if (!std::equal(tHits.begin(), tHits.end(), atHits.begin())) {
        mf::LogVerbatim myprt("TC");
        myprt << someText << " ChkInTraj failed: inTraj - UseHit mis-match for T" << tID
              << " tj.WorkID " << tj.WorkID << " atHits size " << atHits.size() << " tHits size "
              << tHits.size() << " in CTP " << tj.CTP << "\n";
        myprt << "AlgMods: ";
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
        myprt << "\n";
        myprt << "index     inTraj     UseHit \n";
        for (iht = 0; iht < atHits.size(); ++iht) {
          myprt << "iht " << iht << " " << PrintHit(slc.slHits[atHits[iht]]);
          if (iht < tHits.size()) myprt << " " << PrintHit(slc.slHits[tHits[iht]]);
          if (atHits[iht] != tHits[iht]) myprt << " <<< " << atHits[iht] << " != " << tHits[iht];
          myprt << "\n";
        } // iht
        if (tHits.size() > atHits.size()) {
          for (iht = atHits.size(); iht < atHits.size(); ++iht) {
            myprt << "atHits " << iht << " " << PrintHit(slc.slHits[atHits[iht]]) << "\n";
          } // iht
          PrintTrajectory("CIT", slc, tj, USHRT_MAX);
        } // tHit.size > atHits.size()
        return false;
      }
      // check the VtxID
      for (unsigned short end = 0; end < 2; ++end) {
        if (tj.VtxID[end] > slc.vtxs.size()) {
          mf::LogVerbatim("TC") << someText << " ChkInTraj: Bad VtxID " << tj.ID;
          tj.AlgMod[kKilled] = true;
          return false;
        }
      } // end
      ++itj;
    } // tj
    return true;

  } // InTrajOK

  //////////////////////////////////////////
  void
  CheckTrajBeginChg(TCSlice& slc, unsigned short itj)
  {
    // This function is called after the beginning of the tj has been inspected to see if
    // reverse propagation was warranted. Trajectory points at the beginning were removed by
    // this process.
    // A search has been made for a Bragg peak with nothing
    // found. Here we look for a charge pattern like the following, where C means large charge
    // and c means lower charge:
    // CCCCCCccccccc
    // The charge in the two regions should be fairly uniform.

    // This function may split the trajectory so it needs to have been stored
    if (itj > slc.tjs.size() - 1) return;
    auto& tj = slc.tjs[itj];

    if (!tcc.useAlg[kBeginChg]) return;
    if (tj.EndFlag[0][kBragg]) return;
    if (tj.AlgMod[kFTBRvProp]) return;
    if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) return;
    if (tj.Pts.size() < 20) return;

    bool prt = (tcc.dbgSlc && (tcc.dbgStp || tcc.dbgAlg[kBeginChg]));

    // look for a large drop between the average charge near the beginning
    float chg2 = tj.Pts[tj.EndPt[0] + 2].AveChg;
    // and the average charge 15 points away
    float chg15 = tj.Pts[tj.EndPt[0] + 15].AveChg;
    if (chg2 < 3 * chg15) return;

    // find the point where the charge falls below the mid-point
    float midChg = 0.5 * (chg2 + chg15);

    unsigned short breakPt = USHRT_MAX;
    for (unsigned short ipt = tj.EndPt[0] + 3; ipt < 15; ++ipt) {
      float chgm2 = tj.Pts[ipt - 2].Chg;
      if (chgm2 == 0) continue;
      float chgm1 = tj.Pts[ipt - 1].Chg;
      if (chgm1 == 0) continue;
      float chgp1 = tj.Pts[ipt + 1].Chg;
      if (chgp1 == 0) continue;
      float chgp2 = tj.Pts[ipt + 2].Chg;
      if (chgp2 == 0) continue;
      if (chgm2 > midChg && chgm1 > midChg && chgp1 < midChg && chgp2 < midChg) {
        breakPt = ipt;
        break;
      }
    } // breakPt
    if (breakPt == USHRT_MAX) return;
    // check the charge and rms before and after the split
    std::array<double, 2> cnt, sum, sum2;
    for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.Chg <= 0) continue;
      unsigned short end = 0;
      if (ipt > breakPt) end = 1;
      ++cnt[end];
      sum[end] += tp.Chg;
      sum2[end] += tp.Chg * tp.Chg;
    } // ipt
    for (unsigned short end = 0; end < 2; ++end) {
      if (cnt[end] < 3) return;
      double ave = sum[end] / cnt[end];
      double arg = sum2[end] - cnt[end] * ave * ave;
      if (arg <= 0) return;
      sum2[end] = sqrt(arg / (cnt[end] - 1));
      sum2[end] /= ave;
      sum[end] = ave;
    } // region
    bool doSplit = true;
    // don't split if this looks like an electron - no significant improvement
    // in the charge rms before and after
    if (tj.ChgRMS > 0.5 && sum2[0] > 0.3 && sum2[1] > 0.3) doSplit = false;
    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "CTBC: T" << tj.ID << " chgRMS " << tj.ChgRMS;
      myprt << " AveChg before split point " << (int)sum[0] << " rms " << sum2[0];
      myprt << " after " << (int)sum[1] << " rms " << sum2[1] << " doSplit? " << doSplit;
    } // prt
    if (!doSplit) return;
    // Create a vertex at the break point
    VtxStore aVtx;
    aVtx.Pos = tj.Pts[breakPt].Pos;
    aVtx.NTraj = 2;
    aVtx.Pass = tj.Pass;
    aVtx.Topo = 8;
    aVtx.ChiDOF = 0;
    aVtx.CTP = tj.CTP;
    aVtx.ID = slc.vtxs.size() + 1;
    aVtx.Stat[kFixed] = true;
    unsigned short ivx = slc.vtxs.size();
    if (!StoreVertex(slc, aVtx)) return;
    if (!SplitTraj(slc, itj, breakPt, ivx, prt)) {
      if (prt) mf::LogVerbatim("TC") << "CTBC: Failed to split trajectory";
      MakeVertexObsolete("CTBC", slc, slc.vtxs[ivx], false);
      return;
    }
    SetVx2Score(slc);
    slc.tjs[itj].AlgMod[kBeginChg] = true;

    if (prt)
      mf::LogVerbatim("TC") << "CTBC: Split T" << tj.ID << " at "
                            << PrintPos(slc, tj.Pts[breakPt].Pos) << "\n";

  } // CheckTrajBeginChg

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
    slc.tjs[itj].EndFlag[1][kBragg] = true;
    slc.tjs[otj].PDGCode = 13;
    return true;
  } // BraggSplit

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
    if (tj.EndFlag[1][kAtKink]) return;

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
    if (prt) {
      fcnLabel += "-TEPo";
      PrintTrajectory(fcnLabel, slc, tj, USHRT_MAX);
    }

  } // TrimEndPts

  /////////////////////////////////////////
  void
  ChkEndKink(TCSlice& slc, Trajectory& tj, bool prt)
  {
    // look for large-angle kink near the end
    if (!tcc.useAlg[kEndKink]) return;
    if (tj.PDGCode == 111) return;
    if (tj.EndPt[1] - tj.EndPt[0] < 6) return;

    if (prt) mf::LogVerbatim("TC") << "CEK: Inside ChkEndKinks T" << tj.ID << " ";

    float maxSig = tcc.kinkCuts[1];
    unsigned short withNptsFit = 0;
    unsigned short nPtsFit = tcc.kinkCuts[0];
    bool useChg = (tcc.kinkCuts[2] > 0);
    for (unsigned short nptsf = 3; nptsf < nPtsFit; ++nptsf) {
      unsigned short ipt = tj.EndPt[1] - nptsf;
      float ks = KinkSignificance(slc, tj, ipt, nptsf, useChg, prt);
      if (ks > maxSig) {
        maxSig = ks;
        withNptsFit = nptsf;
      }
    } // nptsf
    if (withNptsFit > 0) {
      unsigned short ipt = tj.EndPt[1] - withNptsFit;
      std::cout << "CEK: T" << tj.ID << " ipt " << ipt;
      float ks = KinkSignificance(slc, tj, ipt, withNptsFit, false, prt);
      auto& tp = tj.Pts[ipt];
      std::cout << " " << PrintPos(slc, tp) << " withNptsFit " << withNptsFit << " ks " << ks
                << "\n";
    }

  } // ChkEndKink

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
  void
  MergeGhostTjs(TCSlice& slc, CTP_t inCTP)
  {
    // Merges short Tjs that share many hits with a longer Tj
    if (!tcc.useAlg[kMrgGhost]) return;

    for (auto& shortTj : slc.tjs) {
      if (shortTj.AlgMod[kKilled] || shortTj.AlgMod[kHaloTj]) continue;
      if (shortTj.CTP != inCTP) continue;
      unsigned short spts = shortTj.EndPt[1] - shortTj.EndPt[0];
      if (spts > 20) continue;
      // ignore delta rays
      if (shortTj.PDGCode == 11) continue;
      // ignore InShower Tjs
      if (shortTj.SSID > 0) continue;
      auto tjhits = PutTrajHitsInVector(shortTj, kAllHits);
      if (tjhits.empty()) continue;
      std::vector<int> tids;
      std::vector<unsigned short> tcnt;
      for (auto iht : tjhits) {
        auto& hit = slc.slHits[iht];
        if (hit.InTraj <= 0) continue;
        if ((unsigned int)hit.InTraj > slc.tjs.size()) continue;
        if (hit.InTraj == shortTj.ID) continue;
        unsigned short indx = 0;
        for (indx = 0; indx < tids.size(); ++indx)
          if (hit.InTraj == tids[indx]) break;
        if (indx == tids.size()) {
          tids.push_back(hit.InTraj);
          tcnt.push_back(1);
        }
        else {
          ++tcnt[indx];
        }
      } // iht
      if (tids.empty()) continue;
      // find the max count for Tjs that are longer than this one
      unsigned short maxcnt = 0;
      for (unsigned short indx = 0; indx < tids.size(); ++indx) {
        if (tcnt[indx] > maxcnt) {
          auto& ltj = slc.tjs[tids[indx] - 1];
          unsigned short lpts = ltj.EndPt[1] - ltj.EndPt[0];
          if (lpts < spts) continue;
          maxcnt = tcnt[indx];
        }
      } // indx
      float hitFrac = (float)maxcnt / (float)tjhits.size();
      if (hitFrac < 0.1) continue;
    } // shortTj
  }   // MergeGhostTjs

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
      float wt0 = tj1.Pts[tj1.EndPt[0]].Pos[iwt];
      float wt1 = tj1.Pts[tj1.EndPt[1]].Pos[iwt];
      float lowt1 = wt0;
      float hiwt1 = wt1;
      if (wt1 < lowt1) {
        lowt1 = wt1;
        hiwt1 = wt0;
      }
      // The Lo/Hi wire(time) at each end of tj2
      wt0 = tj2.Pts[tj2.EndPt[0]].Pos[iwt];
      wt1 = tj2.Pts[tj2.EndPt[1]].Pos[iwt];
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
        float dw = tj1.Pts[i1].Pos[0] - tj2.Pts[i2].Pos[0] - dwc;
        if (std::abs(dw) > minSep) continue;
        float dt = tj1.Pts[i1].Pos[1] - tj2.Pts[i2].Pos[1];
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

    // special handling for shower trajectories. UseHit isn't valid
    if (tj.AlgMod[kShowerTj]) {
      for (auto& tp : tj.Pts)
        hitVec.insert(hitVec.end(), tp.Hits.begin(), tp.Hits.end());
      return hitVec;
    } // shower Tj

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
    // need enough points to do a fit on each sideof the presumed kink point
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
    angErr *= tFactor;
    double dang = DeltaAngle(tpPos.Ang, tpNeg.Ang);
    double dangSig = dang / angErr;

    double chgAsym = 0;
    double chgSig = 0;
    if (useChg) {
      // Sum the charge Neg and Pos, excluding the kinkPt
      double chgNeg = 0;
      unsigned short cntNeg = 0;
      for (unsigned short ipt = kinkPt - 1; ipt >= tj.EndPt[0]; --ipt) {
        auto& tp = tj.Pts[ipt];
        if (tp.Chg <= 0) continue;
        chgNeg += tp.Chg;
        ++cntNeg;
        if (cntNeg == nPtsFit) break;
        if (ipt == 0) break;
      } // ipt
      if (cntNeg != nPtsFit) {
        if (prt) mf::LogVerbatim("TC") << " KL: Bad cntNeg " << cntNeg << " != " << nPtsFit;
        return -1;
      }
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
      if (cntPos != nPtsFit) {
        if (prt) mf::LogVerbatim("TC") << " KL: Bad cntPos " << cntPos << " != " << nPtsFit;
        return -1;
      }
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
    }
    return (float)kinkSig;
  } // KinkSignificance

  ////////////////////////////////////////////////
  float
  ElectronLikelihood(const TCSlice& slc, const Trajectory& tj)
  {
    // returns a number between 0 (not electron-like) and 1 (electron-like)
    if (NumPtsWithCharge(slc, tj, false) < 8) return -1;
    if (tj.EndFlag[0][kBragg] || tj.EndFlag[1][kBragg]) return 0;

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
    std::reverse(tj.EndFlag.begin(), tj.EndFlag.end());
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
    if (costh > 1.0 || costh < -1.0) return;
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

    // don't mess with showerTjs or halo tjs
    if (tj.AlgMod[kShowerTj] || tj.AlgMod[kHaloTj]) return;

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
    return MCSMom(slc, tj, tj.EndPt[0], tj.EndPt[1]);
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
    //
    float delta = 5;

    for (auto& mutj : slc.tjs) {
      if (mutj.AlgMod[kKilled]) continue;
      if (mutj.CTP != inCTP) continue;
      if (mutj.PDGCode != 13) continue;
      unsigned short nnear = 0;
      for (unsigned short ipt = mutj.EndPt[0]; ipt <= mutj.EndPt[1]; ++ipt) {
        auto& tp = mutj.Pts[ipt];
        wireWindow[0] = tp.Pos[0];
        wireWindow[1] = tp.Pos[0];
        timeWindow[0] = tp.Pos[1] - delta;
        timeWindow[1] = tp.Pos[1] + delta;
        // get a list of all hits in this region
        bool hitsNear;
        auto closeHits =
          FindCloseHits(slc, wireWindow, timeWindow, plane, kAllHits, true, hitsNear);
        if (closeHits.empty()) continue;
        for (auto iht : closeHits) {
          auto inTraj = slc.slHits[iht].InTraj;
          if (inTraj <= 0) continue;
          if (inTraj == mutj.ID) continue;
          auto& dtj = slc.tjs[inTraj - 1];
          if (dtj.PDGCode == 13) continue;
          for (unsigned short jpt = dtj.EndPt[0]; jpt <= dtj.EndPt[1]; ++jpt) {
            auto& dtp = dtj.Pts[jpt];
            if (std::find(dtp.Hits.begin(), dtp.Hits.end(), iht) == dtp.Hits.end()) continue;
            dtp.Environment[kEnvNearMuon] = true;
            ++nnear;
          } // jpt
        }   // iht
      }     // ipt
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
        unsigned int iht = tp.Hits[ii];
        if (slc.slHits[iht].InTraj == 0) tp.Environment[kEnvUnusedHits] = true;
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
    // TotChg is calculated using
    tj.TotChg = 0;
    tj.AveChg = 0;
    tj.ChgRMS = 0.5;

    // These variables are used to calculate the average and rms using valid points with charge
    double vcnt = 0;
    double vsum = 0;
    double vsum2 = 0;
    // Reject a single large charge TP
    float bigChg = 0;
    for (unsigned short ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.Chg > bigChg) bigChg = tp.Chg;
    } // ipt
    //  variables for calculating the backup quanties. These are only used if npwc < 3
    double bcnt = 0;
    double bsum = 0;
    double bsum2 = 0;
    // don't include the end points
    for (unsigned short ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if (tp.Chg <= 0) continue;
      // ignore the single large charge TP
      if (tp.Chg == bigChg) continue;
      // accumulate a backup sum in case most of the points are overlapped. Note that
      // tp.Chg has an angle correction, which is why the hit integral is summed
      // below. We don't care about this detail for the backup sum
      bsum += tp.Chg;
      bsum2 += tp.Chg * tp.Chg;
      if (tp.Chg > bigChg) bigChg = tp.Chg;
      ++bcnt;
      // Skip TPs that overlap with TPs on other Tjs. A correction will be made below
      if (tj.Pts[ipt].Environment[kEnvOverlap]) continue;
      ++vcnt;
      double tpchg = 0;
      for (unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if (!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        tpchg += (*evt.allHits)[slc.slHits[iht].allHitsIndex].Integral();
      } // ii
      vsum += tpchg;
      vsum2 += tpchg * tpchg;
    } // ipt

    if (bcnt == 0) return;

    if (vcnt < 3) {
      // use the backup sum
      tj.TotChg = bsum;
      tj.AveChg = bsum / bcnt;
      if (vcnt > 2) {
        double arg = bsum2 - bcnt * tj.AveChg * tj.AveChg;
        if (arg > 0) tj.ChgRMS = sqrt(arg / (bcnt - 1));
      }
      for (auto& tp : tj.Pts)
        tp.AveChg = tj.AveChg;
      if (prt)
        mf::LogVerbatim("TC") << inFcnLabel << ".UpdateTjChgProperties: backup sum Set tj.AveChg "
                              << (int)tj.AveChg << " ChgRMS " << tj.ChgRMS;
      return;
    } // low npwc

    double nWires = tj.EndPt[1] - tj.EndPt[0] + 1;
    if (nWires < 2) return;
    // correct for wires missing near vertices.
    // Count the number of wires between vertices at the ends and the first wire
    // that has charge. This code assumes that there should be one TP on each wire
    if (!tj.AlgMod[kPhoton]) {
      for (unsigned short end = 0; end < 2; ++end) {
        if (tj.VtxID[end] == 0) continue;
        auto& tp = tj.Pts[tj.EndPt[end]];
        auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
        int dw = std::abs(tp.Pos[0] - vx2.Pos[0]);
        // This assumes that the vertex is not inside the wire boundaries of the tj
        nWires += dw;
      } // end
    }   // not a photon Tj

    tj.AveChg = vsum / vcnt;
    // calculate the total charge using the tj wire range
    tj.TotChg = nWires * tj.AveChg;
    // calculate the rms
    double arg = vsum2 - vcnt * tj.AveChg * tj.AveChg;
    double rms = 0.5;
    if (arg > 0) rms = sqrt(arg / (vcnt - 1));
    rms /= tj.AveChg;
    // don't let it be an unrealistically low value. It could be crazy large however.
    if (rms < 0.1) rms = 0.1;
    // Don't let the calculated charge RMS dominate until it is well known; after there are 5 - 10 valid TPs.
    // Set the starting charge rms = 0.5
    if (vcnt < 10) {
      double defFrac = 1 / vcnt;
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
    if (vcnt < tcc.nPtsAve) {
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
  NumUsedHitsInTj(const TCSlice& slc, const Trajectory& tj)
  {
    if (tj.AlgMod[kKilled]) return 0;
    if (tj.Pts.empty()) return 0;
    unsigned short nhits = 0;
    for (auto& tp : tj.Pts) {
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii)
        if (tp.UseHit[ii]) ++nhits;
    } // tp
    return nhits;
  } // NumHitsInTj

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

    short npwc = NumPtsWithCharge(slc, tj, false);
    if (npwc < 6) {
      tj.PDGCode = 0;
      return;
    }

    if (tj.Strategy[kStiffEl] && ElectronLikelihood(slc, tj) > tcc.showerTag[6]) {
      tj.PDGCode = 111;
      return;
    }
    if (tj.Strategy[kStiffMu]) {
      tj.PDGCode = 13;
      return;
    }

    if (tcc.showerTag[6] > 0 && ElectronLikelihood(slc, tj) > tcc.showerTag[6]) {
      tj.PDGCode = 11;
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

    if (tcc.modes[kDebug]) {
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
    if (nBadWireFix > 0 && tcc.modes[kDebug]) {
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
    if (tcc.modes[kDebug] && (int)tpc == debug.TPC) {
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

    // Merging shower Tjs requires merging the showers as well.
    if (slc.tjs[itj1].AlgMod[kShowerTj] || slc.tjs[itj2].AlgMod[kShowerTj])
      return MergeShowerTjsAndStore(slc, itj1, itj2, doPrt);

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

    if (tj1.EndFlag[1][kBragg]) {
      if (doPrt)
        mf::LogVerbatim("TC") << "MergeAndStore: You are merging the end of trajectory T" << tj1.ID
                              << " with a Bragg peak. Not merging\n";
      return false;
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

    if (type1Name == "P" && uid <= slc.pfps.size() && (type2Name == "2S" || type2Name == "3S")) {
      // return a list of 3D or 2D showers with the assn 3S -> 2S -> T -> P<ID> or 2S -> T -> P.
      auto& pfp = slc.pfps[uid - 1];
      // First form a list of 2S -> T -> P<ID>
      std::vector<int> ssid;
      for (auto& ss : slc.cots) {
        if (ss.ID <= 0) continue;
        auto shared = SetIntersection(ss.TjIDs, pfp.TjIDs);
        if (!shared.empty() && std::find(ssid.begin(), ssid.end(), ss.ID) == ssid.end())
          ssid.push_back(ss.ID);
      } // ss
      if (type2Name == "2S") return ssid;
      for (auto& ss3 : slc.showers) {
        if (ss3.ID <= 0) continue;
        auto shared = SetIntersection(ss3.CotIDs, ssid);
        if (!shared.empty() && std::find(tmp.begin(), tmp.end(), ss3.ID) == tmp.end())
          tmp.push_back(ss3.ID);
      } // ss3
      return tmp;
    } // 3S -> 2S -> T -> P

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

    if (type1Name == "3S" && uid <= slc.showers.size() && type2Name == "T") {
      // 3S -> T
      auto& ss3 = slc.showers[uid - 1];
      if (ss3.ID == 0) return tmp;
      for (auto cid : ss3.CotIDs) {
        auto& ss = slc.cots[cid - 1];
        if (ss.ID == 0) continue;
        tmp.insert(tmp.end(), ss.TjIDs.begin(), ss.TjIDs.end());
      } // cid
      return tmp;
    } // 3S -> T

    // This isn't strictly necessary but do it for consistency
    if (type1Name == "2S" && uid <= slc.cots.size() && type2Name == "T") {
      // 2S -> T
      auto& ss = slc.cots[uid - 1];
      return ss.TjIDs;
    } // 2S -> T

    if (type1Name == "3S" && uid <= slc.showers.size() && type2Name == "P") {
      // 3S -> P
      auto& ss3 = slc.showers[uid - 1];
      if (ss3.ID == 0) return tmp;
      for (auto cid : ss3.CotIDs) {
        auto& ss = slc.cots[cid - 1];
        if (ss.ID == 0) continue;
        for (auto tid : ss.TjIDs) {
          auto& tj = slc.tjs[tid - 1];
          if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
          if (!tj.AlgMod[kMat3D]) continue;
          for (auto& pfp : slc.pfps) {
            if (pfp.ID <= 0) continue;
            if (std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tj.ID) == pfp.TjIDs.end()) continue;
            if (std::find(tmp.begin(), tmp.end(), pfp.ID) == tmp.end()) tmp.push_back(pfp.ID);
          } // pf
        }   // tid
      }     // cid
      return tmp;
    } // 3S -> P

    if (type1Name == "T" && uid <= slc.tjs.size() && type2Name == "2S") {
      // T -> 2S
      for (auto& ss : slc.cots) {
        if (ss.ID == 0) continue;
        if (std::find(ss.TjIDs.begin(), ss.TjIDs.end(), id) != ss.TjIDs.end()) tmp.push_back(ss.ID);
      } // ss
      return tmp;
    } // T -> 2S

    if (type1Name == "T" && uid <= slc.tjs.size() && type2Name == "3S") {
      // T -> 3S
      for (auto& ss : slc.cots) {
        if (ss.ID == 0) continue;
        if (std::find(ss.TjIDs.begin(), ss.TjIDs.end(), id) == ss.TjIDs.end()) continue;
        if (ss.SS3ID > 0) tmp.push_back(ss.SS3ID);
      } // ss
      return tmp;
    } // T -> 3S

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
    if (tcc.modes[kDebug] && !tcc.dbgStp && !tcc.dbgDump && tcc.dbgSlc && tj.ID == debug.WorkID)
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
    if (tcc.modes[kDebug] && !tcc.dbgStp && !tcc.dbgDump && tcc.dbgSlc && tj.ID == debug.WorkID)
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
      if (typeName == "2S") {
        for (unsigned short indx = 0; indx < slc.cots.size(); ++indx) {
          if (slc.cots[indx].UID == uID) { return std::make_pair(isl, indx); }
        }
      } // 2S
      if (typeName == "3S") {
        for (unsigned short indx = 0; indx < slc.showers.size(); ++indx) {
          if (slc.showers[indx].UID == uID) { return std::make_pair(isl, indx); }
        }
      } // T
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

  ////////////////////////////////////////////////
  bool
  DecodeDebugString(std::string strng)
  {
    // try to unpack the string as Cryostat:TPC:Plane:Wire:Tick or something
    // like Slice:<slice index>

    if (strng == "instruct") {
      std::cout << "****** Unrecognized DebugConfig. Here are your options\n";
      std::cout << " 'C:T:P:W:Tick' where C = cryostat, T = TPC, W = wire, Tick (+/-5) to debug "
                   "stepping (DUNE)\n";
      std::cout << " 'P:W:Tick' for single cryostat/TPC detectors (uB, LArIAT, etc)\n";
      std::cout << " 'WorkID <id> <slice index>' where <id> is a tj work ID (< 0) in slice <slice "
                   "index> (default = 0)\n";
      std::cout << " 'Merge <CTP>' to debug trajectory merging\n";
      std::cout << " '2V <CTP>' to debug 2D vertex finding\n";
      std::cout << " '3V' to debug 3D vertex finding\n";
      std::cout << " 'VxMerge' to debug 2D vertex merging\n";
      std::cout << " 'JunkVx' to debug 2D junk vertex finder\n";
      std::cout << " 'PFP' to debug 3D matching and PFParticles\n";
      std::cout << " 'MVI <MVI> <MVI Iteration>' for detailed debugging of one PFP MatchVecIndex\n";
      std::cout << " 'DeltaRay' to debug delta ray tagging\n";
      std::cout << " 'Muon' to debug muon tagging\n";
      std::cout << " '2S <CTP>' to debug a 2D shower in CTP\n";
      std::cout << " 'Reco TPC <TPC>' to only reconstruct hits in the specified TPC\n";
      std::cout << " 'Reco Slice <ID>' to reconstruct all sub-slices in the recob::Slice with the "
                   "specified ID\n";
      std::cout << " 'SubSlice <sub-slice index>' where <slice index> restricts output to the "
                   "specified sub-slice index\n";
      std::cout << " 'Stitch' to debug PFParticle stitching between TPCs\n";
      std::cout << " 'Sum' or 'Summary' to print a debug summary report\n";
      std::cout << " 'Dump <WorkID>' or 'Dump <UID>' to print all TPs in the trajectory to "
                   "tcdump<UID>.csv\n";
      std::cout << " Note: Algs with debug printing include HamVx, HamVx2, SplitTjCVx, Comp3DVx, "
                   "Comp3DVxIG, VtxHitsSwap\n";
      std::cout << " Set SkipAlgs: [\"bogusText\"] to print a list of algorithm names\n";
      return false;
    } // instruct

    // handle the simple cases that don't need decoding
    if (strng.find("3V") != std::string::npos) {
      tcc.dbg3V = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("3S") != std::string::npos) {
      tcc.dbg3S = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("VxMerge") != std::string::npos) {
      tcc.dbgVxMerge = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("JunkVx") != std::string::npos) {
      tcc.dbgVxJunk = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("DeltaRay") != std::string::npos) {
      tcc.dbgDeltaRayTag = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("Muon") != std::string::npos) {
      tcc.dbgMuonTag = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("Stitch") != std::string::npos) {
      tcc.dbgStitch = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("HamVx") != std::string::npos) {
      tcc.dbgAlg[kHamVx] = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("HamVx2") != std::string::npos) {
      tcc.dbgAlg[kHamVx2] = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (strng.find("Sum") != std::string::npos) {
      tcc.dbgSummary = true;
      tcc.modes[kDebug] = true;
      return true;
    }

    std::vector<std::string> words;
    boost::split(words, strng, boost::is_any_of(" :"), boost::token_compress_on);
    if (words.size() == 5) {
      // configure for DUNE
      debug.Cryostat = std::stoi(words[0]);
      debug.TPC = std::stoi(words[1]);
      debug.Plane = std::stoi(words[2]);
      debug.Wire = std::stoi(words[3]);
      debug.Tick = std::stoi(words[4]);
      tcc.modes[kDebug] = true;
      tcc.dbgStp = true;
      // also dump this tj
      tcc.dbgDump = true;
      return true;
    } // nums.size() == 5
    if (words[0] == "PFP" || words[0] == "MVI") {
      tcc.dbgPFP = true;
      tcc.modes[kDebug] = true;
      // Use debug.Hit to identify the matchVec index
      if (words.size() > 2) {
        debug.MVI = std::stoi(words[1]);
        if (words.size() == 3) debug.MVI_Iter = std::stoi(words[2]);
      }
      return true;
    } // PFP
    if (words.size() == 2 && words[0] == "Dump") {
      debug.WorkID = std::stoi(words[1]);
      debug.Slice = 0;
      tcc.modes[kDebug] = true;
      tcc.dbgDump = true;
      return true;
    }
    if (words.size() > 1 && words[0] == "WorkID") {
      debug.WorkID = std::stoi(words[1]);
      if (debug.WorkID >= 0) return false;
      // default to sub-slice index 0
      debug.Slice = 0;
      if (words.size() > 2) debug.Slice = std::stoi(words[2]);
      tcc.modes[kDebug] = true;
      // dbgStp is set true after debug.WorkID is found
      tcc.dbgStp = false;
      return true;
    } // words.size() == 3 && words[0] == "WorkID"
    if (words.size() == 3 && words[0] == "Reco" && words[1] == "TPC") {
      tcc.recoTPC = std::stoi(words[2]);
      tcc.modes[kDebug] = true;
      std::cout << "Reconstructing only in TPC " << tcc.recoTPC << "\n";
      return true;
    }
    if(words.size() == 3 && words[0] == "Reco" && words[1] == "Slice") {
      tcc.recoSlice = std::stoi(words[2]);
      std::cout<<"Reconstructing Slice "<<tcc.recoSlice<<"\n";
      return true;
    }
    if (words.size() == 3) {
      // configure for uB, LArIAT, etc
      debug.Cryostat = 0;
      debug.TPC = 0;
      debug.Plane = std::stoi(words[0]);
      debug.Wire = std::stoi(words[1]);
      debug.Tick = std::stoi(words[2]);
      debug.CTP = EncodeCTP(debug.Cryostat, debug.TPC, debug.Plane);
      tcc.modes[kDebug] = true;
      tcc.dbgStp = true;
      return true;
    }
    if (words.size() == 2 && words[0] == "Merge") {
      debug.CTP = std::stoi(words[1]);
      tcc.dbgMrg = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (words.size() == 2 && words[0] == "2V") {
      debug.CTP = std::stoi(words[1]);
      tcc.dbg2V = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    if (words.size() == 2 && words[0] == "2S") {
      debug.CTP = std::stoi(words[1]);
      tcc.dbg2S = true;
      tcc.modes[kDebug] = true;
      return true;
    }
    // Slice could apply to several debug options.
    if (words.size() == 2 && words[0] == "SubSlice") {
      debug.Slice = std::stoi(words[1]);
      return true;
    }
    return false;
  } // DecodeDebugString

  // ****************************** Printing  ******************************

  void
  DumpTj()
  {
    // Dump all of the points in a trajectory to the output in a form that can
    // be imported by another application, e.g. Excel
    // Search for the trajectory with the specified WorkID or Unique ID

    for (auto& slc : slices) {
      for (auto& tj : slc.tjs) {
        if (tj.WorkID != debug.WorkID && tj.UID != debug.WorkID) continue;
        // print a header
        std::ofstream outfile;
        std::string fname = "tcdump" + std::to_string(tj.UID) + ".csv";
        outfile.open(fname, std::ios::out | std::ios::trunc);
        outfile << "Dump trajectory T" << tj.UID << " WorkID " << tj.WorkID;
        outfile << " ChgRMS " << std::setprecision(2) << tj.ChgRMS;
        outfile << "\n";
        outfile << "Wire, Chg T" << tj.UID
                << ", totChg, Tick, Delta, NTPsFit, Ang, ChiDOF, KinkSig, HitPosErr\n";
        for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          auto& tp = tj.Pts[ipt];
          outfile << std::fixed;
          outfile << std::setprecision(0) << std::nearbyint(tp.Pos[0]);
          outfile << "," << (int)tp.Chg;
          // total charge near the TP
          float totChg = 0;
          for (auto iht : tp.Hits) {
            auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
            totChg += hit.Integral();
          }
          outfile << "," << (int)totChg;
          outfile << "," << std::setprecision(0) << std::nearbyint(tp.Pos[1] / tcc.unitsPerTick);
          outfile << "," << std::setprecision(2) << tp.Delta;
          outfile << "," << tp.NTPsFit;
          outfile << "," << std::setprecision(3) << tp.Ang;
          outfile << "," << std::setprecision(2) << tp.FitChi;
          outfile << "," << std::setprecision(2) << tp.KinkSig;
          outfile << "," << std::setprecision(2) << sqrt(tp.HitPosErr2);
          outfile << "\n";
        } // ipt
        outfile.close();
        std::cout<<"Points on T"<<tj.UID<<" dumped to "<<fname<<"\n";
//        tcc.dbgDump = false;
//        return;
      } // tj
    }   // slc

  } // DumpTj

  ////////////////////////////////////////////////
  void
  PrintDebugMode()
  {
    // print the debug mode configuration to the screen
    std::cout << "*** TrajCluster debug mode configuration in";
    std::cout << " CTP=";
    if (debug.CTP == UINT_MAX) { std::cout << "NA"; }
    else {
      std::cout << debug.CTP;
    }
    std::cout << " Cryostat=" << debug.Cryostat;
    std::cout << " TPC=" << debug.TPC;
    std::cout << " Plane=" << debug.Plane;
    std::cout << " Wire=" << debug.Wire;
    std::cout << " Tick=" << debug.Tick;
    std::cout << " Hit=";
    if (debug.Hit == UINT_MAX) { std::cout << "NA"; }
    else {
      std::cout << debug.Hit;
    }
    std::cout << " WorkID=";
    if (debug.WorkID == 0) { std::cout << "NA"; }
    else {
      std::cout << debug.WorkID;
    }
    std::cout << " Slice=";
    if (debug.Slice == -1) { std::cout << "All"; }
    else {
      std::cout << debug.Slice;
    }
    std::cout << "\n";
    std::cout << "*** tcc.dbg modes:";
    if (tcc.dbgSlc) std::cout << " dbgSlc";
    if (tcc.dbgStp) std::cout << " dbgStp";
    if (tcc.dbgMrg) std::cout << " dbgMrg";
    if (tcc.dbg2V) std::cout << " dbg2V";
    if (tcc.dbg2S) std::cout << " dbg2S";
    if (tcc.dbgVxNeutral) std::cout << " dbgVxNeutral";
    if (tcc.dbgVxMerge) std::cout << " dbgVxMerge";
    if (tcc.dbgVxJunk) std::cout << " dbgVxJunk";
    if (tcc.dbg3V) std::cout << " dbg3V";
    if (tcc.dbgPFP) std::cout << " dbgPFP";
    if (tcc.dbgDeltaRayTag) std::cout << " dbgDeltaRayTag";
    if (tcc.dbgMuonTag) std::cout << " dbgMuonTag";
    if (tcc.dbgStitch) std::cout << " dbgStitch";
    if (tcc.dbgSummary) std::cout << " dbgSummary";
    if (tcc.dbgDump) std::cout << " dbgDump";
    std::cout << "\n";
    std::cout << "*** Using algs:";
    unsigned short cnt = 0;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
      if (tcc.useAlg[ib] && ib != kKilled) {
        ++cnt;
        if (cnt % 10 == 0) std::cout << "\n   ";
        std::cout << " " << AlgBitNames[ib];
      }
    }
    std::cout << "\n";
    std::cout << "*** Skipping algs:";
    cnt = 0;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) {
      if (!tcc.useAlg[ib] && ib != kKilled) {
        ++cnt;
        if (cnt % 10 == 0) std::cout << "\n   ";
        std::cout << " " << AlgBitNames[ib];
      }
    }
    std::cout << "\n";
  } // PrintDebugMode

  ////////////////////////////////////////////////
  void PrintAll(detinfo::DetectorPropertiesData const& detProp, std::string someText)
  {
    // print everything in all slices
    bool prt3V = false;
    bool prt2V = false;
    bool prtT = false;
    bool prtP = false;
    bool prtS3 = false;
    for (size_t isl = 0; isl < slices.size(); ++isl) {
      if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
      auto& slc = slices[isl];
      if (!slc.vtx3s.empty()) prt3V = true;
      if (!slc.vtxs.empty()) prt2V = true;
      if (!slc.tjs.empty()) prtT = true;
      if (!slc.pfps.empty()) prtP = true;
      if (!slc.showers.empty()) prtS3 = true;
    } // slc
    mf::LogVerbatim myprt("TC");
    myprt << "Debug report from caller " << someText << "\n";
    myprt << " 'prodID' = <sliceID>:<subSliceIndex>:<productID>/<productUID>\n";
    if (prtS3) {
      myprt << "************ Showers ************\n";
      myprt << "     prodID      Vtx  parUID  ___ChgPos____ ______Dir_____ ____posInPln____ "
               "___projInPln____ 2D shower UIDs\n";
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.showers.empty()) continue;
        for (auto& ss3 : slc.showers)
          Print3S(detProp, someText, myprt, ss3);
      } // slc
    }   // prtS3
    if (prtP) {
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.pfps.empty()) continue;
        for (auto& pfp : slc.pfps)
          PrintP(someText, myprt, pfp, printHeader);
      } // slc
    }   // prtS3
    if (prt3V) {
      bool printHeader = true;
      myprt << "****** 3D vertices "
               "******************************************__2DVtx_UID__*******\n";
      myprt << "     prodID    Cstat TPC     X       Y       Z    XEr  YEr  "
               "ZEr pln0 pln1 pln2 Wire score Prim? Nu? nTru";
      myprt << " ___________2D_Pos____________ _____Tj UIDs________\n";
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.vtx3s.empty()) continue;
        for (auto& vx3 : slc.vtx3s)
          Print3V(detProp, someText, myprt, vx3, printHeader);
      } // slc
    }   // prt3V
    if (prt2V) {
      bool printHeader = true;
      myprt << "************ 2D vertices ************\n";
      myprt << "     prodID      CTP  wire  err   tick   err  ChiDOF  NTj Pass "
               " Topo ChgFrac Score  v3D Tj UIDs\n";
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.vtxs.empty()) continue;
        for (auto& vx2 : slc.vtxs)
          Print2V(someText, myprt, vx2, printHeader);
      } // slc
    }   // prt2V
    if (prtT) {
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.tjs.empty()) continue;
        for (auto& tj : slc.tjs)
          PrintT(someText, myprt, tj, printHeader);
      } // slc
    }   // prtT
  }     // PrintAll

  ////////////////////////////////////////////////
  void
  PrintP(std::string someText, mf::LogVerbatim& myprt, PFPStruct& pfp, bool& printHeader)
  {
    if (pfp.ID <= 0) return;
    if (printHeader) {
      myprt << "************ PFParticles ************\n";
      myprt << "     prodID    sVx  _____sPos____ CS _______sDir______ ____sdEdx_____    eVx  "
               "_____ePos____ CS ____edEdx_____  MVI MCSMom  Len nTP3 nSec SLk? PDG  Par \n";
      printHeader = false;
    } // printHeader
    auto sIndx = GetSliceIndex("P", pfp.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    std::string str =
      std::to_string(slc.ID) + ":" + std::to_string(sIndx.first) + ":" + std::to_string(pfp.ID);
    str += "/" + std::to_string(pfp.UID);
    myprt << std::setw(12) << str;
    // start and end stuff
    for (unsigned short end = 0; end < 2; ++end) {
      str = "--";
      if (pfp.Vx3ID[end] > 0) str = "3V" + std::to_string(slc.vtx3s[pfp.Vx3ID[end] - 1].UID);
      myprt << std::setw(6) << str;
      myprt << std::fixed << std::right << std::setprecision(0);
      auto pos = PosAtEnd(pfp, end);
      myprt << std::setw(5) << pos[0];
      myprt << std::setw(5) << pos[1];
      myprt << std::setw(5) << pos[2];
      // print character for Outside or Inside the FV
      if (InsideFV(slc, pfp, end)) { myprt << "  I"; }
      else {
        myprt << "  O";
      }
      // only print the starting direction
      if (end == 0) {
        myprt << std::fixed << std::right << std::setprecision(2);
        auto dir = DirAtEnd(pfp, end);
        myprt << std::setw(6) << dir[0];
        myprt << std::setw(6) << dir[1];
        myprt << std::setw(6) << dir[2];
      } // end == 0
      for (auto& dedx : pfp.dEdx[end]) {
        if (dedx < 50) { myprt << std::setw(5) << std::setprecision(1) << dedx; }
        else {
          myprt << std::setw(5) << std::setprecision(0) << dedx;
        }
      } // dedx
      if (pfp.dEdx[end].size() < 3) {
        for (size_t i = 0; i < 3 - pfp.dEdx[end].size(); ++i) {
          myprt << std::setw(6) << ' ';
        }
      }
    } // startend
    myprt << std::setw(6) << pfp.MVI;
    // global stuff
    myprt << std::setw(7) << MCSMom(slc, pfp.TjIDs);
    float length = Length(pfp);
    if (length < 100) { myprt << std::setw(5) << std::setprecision(1) << length; }
    else {
      myprt << std::setw(5) << std::setprecision(0) << length;
    }
    myprt << std::setw(5) << pfp.TP3Ds.size();
    myprt << std::setw(5) << pfp.SectionFits.size();
    myprt << std::setw(5) << IsShowerLike(slc, pfp.TjIDs);
    myprt << std::setw(5) << pfp.PDGCode;
    myprt << std::setw(4) << pfp.ParentUID;
    if (!pfp.TjIDs.empty()) {
      if (pfp.TjUIDs.empty()) {
        // print Tjs in one TPC
        for (auto tjid : pfp.TjIDs)
          myprt << " TU" << slc.tjs[tjid - 1].UID;
      }
      else {
        // print Tjs in all TPCs (if this is called after FinishEvent)
        for (auto tjuid : pfp.TjUIDs)
          myprt << " TU" << tjuid;
      }
    } // TjIDs exist
    if (!pfp.DtrUIDs.empty()) {
      myprt << " dtrs";
      for (auto dtruid : pfp.DtrUIDs)
        myprt << " PU" << dtruid;
    } // dtr ids exist
    myprt << "\n";
  } // PrintP

  ////////////////////////////////////////////////
  void
  Print3V(detinfo::DetectorPropertiesData const& detProp,
          std::string someText,
          mf::LogVerbatim& myprt,
          Vtx3Store& vx3,
          bool& printHeader)
  {
    // print a 3D vertex on one line
    if (vx3.ID <= 0) return;
    auto sIndx = GetSliceIndex("3V", vx3.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    if (printHeader) {
      myprt
        << "****** 3D vertices ******************************************__2DVtx_UID__*******\n";
      myprt << "     prodID    Cstat TPC     X       Y       Z    pln0   pln1   pln2 Wire score "
               "Prim? Nu? nTru";
      myprt << " ___________2D_Pos____________ _____Tj UIDs________\n";
      printHeader = false;
    }
    std::string str = "3V" + std::to_string(vx3.ID) + "/3VU" + std::to_string(vx3.UID);
    myprt << std::right << std::setw(12) << std::fixed << str;
    myprt << std::setprecision(0);
    myprt << std::right << std::setw(7) << vx3.TPCID.Cryostat;
    myprt << std::right << std::setw(5) << vx3.TPCID.TPC;
    myprt << std::right << std::setw(8) << vx3.X;
    myprt << std::right << std::setw(8) << vx3.Y;
    myprt << std::right << std::setw(8) << vx3.Z;
    for (auto vx2id : vx3.Vx2ID) {
      if (vx2id > 0) {
        str = "2VU" + std::to_string(slc.vtxs[vx2id - 1].UID);
        myprt << std::right << std::setw(7) << str;
      }
      else {
        myprt << "   --";
      }
    } // vx2id
    myprt << std::right << std::setw(5) << vx3.Wire;
    unsigned short nTruMatch = 0;
    for (unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
      if (vx3.Vx2ID[ipl] == 0) continue;
      unsigned short iv2 = vx3.Vx2ID[ipl] - 1;
      if (slc.vtxs[iv2].Stat[kVxTruMatch]) ++nTruMatch;
    } // ipl
    myprt << std::right << std::setw(6) << std::setprecision(1) << vx3.Score;
    myprt << std::setw(6) << vx3.Primary;
    myprt << std::setw(4) << vx3.Neutrino;
    myprt << std::right << std::setw(5) << nTruMatch;
    Point2_t pos;
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      PosInPlane(detProp, slc, vx3, plane, pos);
      myprt << " " << PrintPos(slc, pos);
    } // plane
    if (vx3.Wire == -2) {
      // find the Tjs that are attached to it
      for (unsigned short end = 0; end < 2; ++end) {
        for (auto& pfp : slc.pfps) {
          if (pfp.Vx3ID[end] == vx3.ID) {
            for (auto tjID : pfp.TjIDs) {
              auto& tj = slc.tjs[tjID - 1];
              myprt << " T" << tj.UID;
            } // tjID
          }   // pfp.Vx3ID[0] == vx3.ID
        }     // pfp
      }       // end
    }
    else {
      auto vxtjs = GetAssns(slc, "3V", vx3.ID, "T");
      for (auto tjid : vxtjs) {
        auto& tj = slc.tjs[tjid - 1];
        myprt << " TU" << tj.UID;
      }
    } // vx3.Wire != -2
    myprt << "\n";
  } // Print3V

  ////////////////////////////////////////////////
  void
  Print2V(std::string someText, mf::LogVerbatim& myprt, VtxStore& vx2, bool& printHeader)
  {
    // print a 2D vertex on one line
    if (vx2.ID <= 0) return;
    if (debug.CTP != UINT_MAX && vx2.CTP != debug.CTP) return;
    auto sIndx = GetSliceIndex("2V", vx2.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    if (printHeader) {
      myprt << "************ 2D vertices ************\n";
      myprt << "     prodID    CTP    wire  err   tick   err  ChiDOF  NTj Pass  Topo ChgFrac Score "
               " v3D Tj UIDs\n";
      printHeader = false;
    }
    std::string str = "2V" + std::to_string(vx2.ID) + "/2VU" + std::to_string(vx2.UID);
    myprt << std::right << std::setw(12) << std::fixed << str;
    myprt << std::right << std::setw(6) << vx2.CTP;
    myprt << std::right << std::setw(8) << std::setprecision(0) << std::nearbyint(vx2.Pos[0]);
    myprt << std::right << std::setw(5) << std::setprecision(1) << vx2.PosErr[0];
    myprt << std::right << std::setw(8) << std::setprecision(0)
          << std::nearbyint(vx2.Pos[1] / tcc.unitsPerTick);
    myprt << std::right << std::setw(5) << std::setprecision(1) << vx2.PosErr[1] / tcc.unitsPerTick;
    myprt << std::right << std::setw(7) << vx2.ChiDOF;
    myprt << std::right << std::setw(5) << vx2.NTraj;
    myprt << std::right << std::setw(5) << vx2.Pass;
    myprt << std::right << std::setw(6) << vx2.Topo;
    myprt << std::right << std::setw(9) << std::setprecision(2) << vx2.TjChgFrac;
    myprt << std::right << std::setw(6) << std::setprecision(1) << vx2.Score;
    int v3id = 0;
    if (vx2.Vx3ID > 0) v3id = slc.vtx3s[vx2.Vx3ID - 1].UID;
    myprt << std::right << std::setw(5) << v3id;
    myprt << "    ";
    // display the traj IDs
    for (unsigned short ii = 0; ii < slc.tjs.size(); ++ii) {
      auto const& tj = slc.tjs[ii];
      if (tj.AlgMod[kKilled]) continue;
      for (unsigned short end = 0; end < 2; ++end) {
        if (tj.VtxID[end] != (short)vx2.ID) continue;
        std::string tid = " TU" + std::to_string(tj.UID) + "_" + std::to_string(end);
        myprt << std::right << std::setw(6) << tid;
      } // end
    }   // ii
    myprt << " Stat:";
    // Special flags. Ignore the first flag bit (0 = kVxTrjTried) which is done for every vertex
    for (unsigned short ib = 1; ib < VtxBitNames.size(); ++ib)
      if (vx2.Stat[ib]) myprt << " " << VtxBitNames[ib];
    myprt << "\n";
  } // Print2V

  ////////////////////////////////////////////////
  void
  Print3S(detinfo::DetectorPropertiesData const& detProp,
          std::string someText,
          mf::LogVerbatim& myprt,
          ShowerStruct3D& ss3)
  {
    if (ss3.ID <= 0) return;
    auto sIndx = GetSliceIndex("3S", ss3.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    std::string str =
      std::to_string(slc.ID) + ":" + std::to_string(sIndx.first) + ":" + std::to_string(ss3.ID);
    str += "/" + std::to_string(ss3.UID);
    myprt << std::fixed << std::setw(12) << str;
    str = "--";
    if (ss3.Vx3ID > 0) str = "3V" + std::to_string(slc.vtx3s[ss3.Vx3ID - 1].UID);
    myprt << std::setw(6) << str;
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      myprt << std::setprecision(0) << std::setw(5) << ss3.ChgPos[xyz];
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      myprt << std::setprecision(2) << std::setw(5) << ss3.Dir[xyz];
    std::vector<float> projInPlane(slc.nPlanes);
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(ss3.TPCID.Cryostat, ss3.TPCID.TPC, plane);
      auto tp = MakeBareTP(detProp, slc, ss3.ChgPos, ss3.Dir, inCTP);
      myprt << " " << PrintPos(slc, tp.Pos);
      projInPlane[plane] = tp.Delta;
    } // plane
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      myprt << std::setprecision(2) << std::setw(5) << projInPlane[plane];
    } // plane
    for (auto cid : ss3.CotIDs) {
      auto& ss = slc.cots[cid - 1];
      str = "2SU" + std::to_string(ss.UID);
      myprt << std::setw(5) << str;
    } // ci
    if (ss3.NeedsUpdate) myprt << " *** Needs update";
    myprt << "\n";
  } // Print3S

  ////////////////////////////////////////////////
  void
  PrintT(std::string someText, mf::LogVerbatim& myprt, Trajectory& tj, bool& printHeader)
  {
    // print a 2D vertex on one line
    if(tj.ID <= 0) return;
    if(debug.CTP != UINT_MAX && tj.CTP != debug.CTP) return;
    if(printHeader) {
      myprt<<"************ Trajectories ************\n";
      myprt<<"Tj AngleCode-EndFlag decoder (EF): <AngleCode> + <end flag>";
      myprt<<" (B=Bragg Peak, V=Vertex, A=AngleKink, C=ChargeKink, T=Trajectory, S=StartEnd)\n";
      myprt<<"     prodID    CTP Pass  Pts     W:T      Ang EF AveQ     W:T      Ang EF AveQ Chg(k) chgRMS  Mom __Vtx__  PDG eLike  Par Pri NuPar   WorkID \n";
      printHeader = false;
    }
    auto sIndx = GetSliceIndex("T", tj.UID);
    if (sIndx.first == USHRT_MAX) return;
    auto& slc = slices[sIndx.first];
    std::string str = "T" + std::to_string(tj.ID) + "/TU" + std::to_string(tj.UID);
    myprt << std::fixed << std::setw(12) << str;
    myprt << std::setw(6) << tj.CTP;
    myprt << std::setw(5) << tj.Pass;
    myprt << std::setw(5) << tj.EndPt[1] - tj.EndPt[0] + 1;
    unsigned short endPt0 = tj.EndPt[0];
    auto& tp0 = tj.Pts[endPt0];
    int itick = tp0.Pos[1]/tcc.unitsPerTick;
    if(itick < 0) itick = 0;
    myprt<<std::setw(6)<<(int)(tp0.Pos[0]+0.5)<<":"<<itick; // W:T
    if(itick < 10) { myprt<<" "; }
    if(itick < 100) { myprt<<" "; }
    if(itick < 1000) { myprt<<" "; }
    myprt<<std::setw(6)<<std::setprecision(2)<<tp0.Ang;
    myprt<<std::setw(2)<<tp0.AngleCode;
    if(tj.EndFlag[0][kBragg]) {
      myprt<<"B";
    } else if(tj.EndFlag[0][kAtVtx]) {
      myprt<<"V";
    } else if(tj.EndFlag[0][kAtKink]) {
      myprt<<"K";
    } else if(tj.EndFlag[0][kAtTj]) {
      myprt<<"T";
    } else {
      myprt<<" ";
    }
    if(tj.StartEnd == 0) { myprt<<"S"; } else { myprt<<" "; }
    myprt<<std::setw(5)<<(int)tp0.AveChg;
    unsigned short endPt1 = tj.EndPt[1];
    auto& tp1 = tj.Pts[endPt1];
    itick = tp1.Pos[1]/tcc.unitsPerTick;
    myprt<<std::setw(6)<<(int)(tp1.Pos[0]+0.5)<<":"<<itick; // W:T
    if(itick < 10) { myprt<<" "; }
    if(itick < 100) { myprt<<" "; }
    if(itick < 1000) { myprt<<" "; }
    myprt<<std::setw(6)<<std::setprecision(2)<<tp1.Ang;
    myprt<<std::setw(2)<<tp1.AngleCode;
    if(tj.EndFlag[1][kBragg]) {
      myprt<<"B";
    } else if(tj.EndFlag[1][kAtVtx]) {
      myprt<<"V";
    } else if(tj.EndFlag[1][kAtKink]) {
      myprt<<"K";
    } else if(tj.EndFlag[1][kAtTj]) {
      myprt<<"T";
    } else {
      myprt<<" ";
    }
    if(tj.StartEnd == 1) myprt<<"S";
    myprt<<std::setw(5)<<(int)tp1.AveChg;
    myprt<<std::setw(7)<<std::setprecision(1)<<tj.TotChg/1000;
    myprt<<std::setw(7)<<std::setprecision(2)<<tj.ChgRMS;
    myprt<<std::setw(5)<<tj.MCSMom;
    int vxid = 0;
    if (tj.VtxID[0] > 0) vxid = slc.vtxs[tj.VtxID[0] - 1].UID;
    myprt << std::setw(4) << vxid;
    vxid = 0;
    if (tj.VtxID[1] > 0) vxid = slc.vtxs[tj.VtxID[1] - 1].UID;
    myprt << std::setw(4) << vxid;
    myprt << std::setw(5) << tj.PDGCode;
    myprt << std::setw(7) << std::setprecision(2) << ElectronLikelihood(slc, tj);
    myprt << std::setw(5) << tj.ParentID;
    myprt << std::setw(5) << PrimaryID(slc, tj);
    myprt << std::setw(6) << NeutrinoPrimaryTjID(slc, tj);
    myprt << std::setw(7) << tj.WorkID;
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
      if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
    for (unsigned short ib = 0; ib < StrategyBitNames.size(); ++ib)
      if (tj.Strategy[ib]) myprt << " " << StrategyBitNames[ib];
    myprt << "\n";
  } // PrintT

  ////////////////////////////////////////////////
  void
  PrintAllTraj(detinfo::DetectorPropertiesData const& detProp,
               std::string someText,
               TCSlice& slc,
               unsigned short itj,
               unsigned short ipt,
               bool prtVtx)
  {

    mf::LogVerbatim myprt("TC");

    if (prtVtx) {
      if (!slc.vtx3s.empty()) {
        // print out 3D vertices
        myprt
          << someText
          << "****** 3D vertices ******************************************__2DVtx_ID__*******\n";
        myprt << someText
              << "  Vtx  Cstat  TPC     X       Y       Z    XEr  YEr  ZEr pln0 pln1 pln2 Wire "
                 "score Prim? Nu? nTru";
        myprt << " ___________2D_Pos____________ _____Tjs________\n";
        for (unsigned short iv = 0; iv < slc.vtx3s.size(); ++iv) {
          if (slc.vtx3s[iv].ID == 0) continue;
          const Vtx3Store& vx3 = slc.vtx3s[iv];
          myprt << someText;
          std::string vid = "3v" + std::to_string(vx3.ID);
          myprt << std::right << std::setw(5) << std::fixed << vid;
          myprt << std::setprecision(1);
          myprt << std::right << std::setw(7) << vx3.TPCID.Cryostat;
          myprt << std::right << std::setw(5) << vx3.TPCID.TPC;
          myprt << std::right << std::setw(8) << vx3.X;
          myprt << std::right << std::setw(8) << vx3.Y;
          myprt << std::right << std::setw(8) << vx3.Z;
          myprt << std::right << std::setw(5) << vx3.XErr;
          myprt << std::right << std::setw(5) << vx3.YErr;
          myprt << std::right << std::setw(5) << vx3.ZErr;
          myprt << std::right << std::setw(5) << vx3.Vx2ID[0];
          myprt << std::right << std::setw(5) << vx3.Vx2ID[1];
          myprt << std::right << std::setw(5) << vx3.Vx2ID[2];
          myprt << std::right << std::setw(5) << vx3.Wire;
          unsigned short nTruMatch = 0;
          for (unsigned short ipl = 0; ipl < slc.nPlanes; ++ipl) {
            if (vx3.Vx2ID[ipl] == 0) continue;
            unsigned short iv2 = vx3.Vx2ID[ipl] - 1;
            if (slc.vtxs[iv2].Stat[kVxTruMatch]) ++nTruMatch;
          } // ipl
          myprt << std::right << std::setw(6) << std::setprecision(1) << vx3.Score;
          myprt << std::setw(6) << vx3.Primary;
          myprt << std::setw(4) << vx3.Neutrino;
          myprt << std::right << std::setw(5) << nTruMatch;
          Point2_t pos;
          for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
            PosInPlane(detProp, slc, vx3, plane, pos);
            myprt << " " << PrintPos(slc, pos);
          } // plane
          if (vx3.Wire == -2) {
            // find the Tjs that are attached to it
            for (auto& pfp : slc.pfps) {
              if (pfp.Vx3ID[0] == slc.vtx3s[iv].ID) {
                for (auto& tjID : pfp.TjIDs)
                  myprt << " t" << tjID;
              }
              if (pfp.Vx3ID[1] == slc.vtx3s[iv].ID) {
                for (auto& tjID : pfp.TjIDs)
                  myprt << " t" << tjID;
              }
            } // ipfp
          }
          else {
            auto vxtjs = GetAssns(slc, "3V", vx3.ID, "T");
            for (auto tjid : vxtjs)
              myprt << " t" << tjid;
          }
          myprt << "\n";
        }
      } // slc.vtx3s.size
      if (!slc.vtxs.empty()) {
        bool foundOne = false;
        for (unsigned short iv = 0; iv < slc.vtxs.size(); ++iv) {
          auto& vx2 = slc.vtxs[iv];
          if (debug.Plane < 3 && debug.Plane != (int)DecodeCTP(vx2.CTP).Plane) continue;
          if (vx2.NTraj == 0) continue;
          foundOne = true;
        } // iv
        if (foundOne) {
          // print out 2D vertices
          myprt << someText << "************ 2D vertices ************\n";
          myprt << someText
                << " ID   CTP    wire  err   tick   err  ChiDOF  NTj Pass  Topo ChgFrac Score  v3D "
                   "TjIDs\n";
          for (auto& vx2 : slc.vtxs) {
            if (vx2.ID == 0) continue;
            if (debug.Plane < 3 && debug.Plane != (int)DecodeCTP(vx2.CTP).Plane) continue;
            myprt << someText;
            std::string vid = "2v" + std::to_string(vx2.ID);
            myprt << std::right << std::setw(5) << std::fixed << vid;
            myprt << std::right << std::setw(6) << vx2.CTP;
            myprt << std::right << std::setw(8) << std::setprecision(0)
                  << std::nearbyint(vx2.Pos[0]);
            myprt << std::right << std::setw(5) << std::setprecision(1) << vx2.PosErr[0];
            myprt << std::right << std::setw(8) << std::setprecision(0)
                  << std::nearbyint(vx2.Pos[1] / tcc.unitsPerTick);
            myprt << std::right << std::setw(5) << std::setprecision(1)
                  << vx2.PosErr[1] / tcc.unitsPerTick;
            myprt << std::right << std::setw(7) << vx2.ChiDOF;
            myprt << std::right << std::setw(5) << vx2.NTraj;
            myprt << std::right << std::setw(5) << vx2.Pass;
            myprt << std::right << std::setw(6) << vx2.Topo;
            myprt << std::right << std::setw(9) << std::setprecision(2) << vx2.TjChgFrac;
            myprt << std::right << std::setw(6) << std::setprecision(1) << vx2.Score;
            myprt << std::right << std::setw(5) << vx2.Vx3ID;
            myprt << "    ";
            // display the traj IDs
            for (unsigned short ii = 0; ii < slc.tjs.size(); ++ii) {
              auto const& aTj = slc.tjs[ii];
              if (debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
              if (aTj.AlgMod[kKilled]) continue;
              for (unsigned short end = 0; end < 2; ++end) {
                if (aTj.VtxID[end] != (short)vx2.ID) continue;
                std::string tid = " t" + std::to_string(aTj.ID) + "_" + std::to_string(end);
                myprt << std::right << std::setw(6) << tid;
              } // end
            }   // ii
            // Special flags. Ignore the first flag bit (0 = kVxTrjTried) which is done for every vertex
            for (unsigned short ib = 1; ib < VtxBitNames.size(); ++ib)
              if (vx2.Stat[ib]) myprt << " " << VtxBitNames[ib];
            myprt << "\n";
          } // iv
        }
      } // slc.vtxs.size
    }

    if (slc.tjs.empty()) {
      mf::LogVerbatim("TC") << someText << " No allTraj trajectories to print";
      return;
    }

    // Print all trajectories in slc.tjs if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    if (itj == USHRT_MAX) {
      // Print summary trajectory information
      myprt << "Tj AngleCode-EndFlag (EF) decoder: <AngleCode> + <reason for stopping>";
      myprt << " (B=Bragg Peak, V=Vertex, A=AngleKink, C=ChargeKink, T=Trajectory)\n";
      std::vector<unsigned int> tmp;
      myprt << someText
            << "   UID   CTP Pass  Pts     W:T      Ang EF AveQ     W:T      Ang EF AveQ Chg(k) "
               "chgRMS  Mom SDr __Vtx__  PDG  Par Pri NuPar   WorkID \n";
      for (unsigned short ii = 0; ii < slc.tjs.size(); ++ii) {
        auto& aTj = slc.tjs[ii];
        if (debug.CTP != UINT_MAX && aTj.CTP != debug.CTP) continue;
        myprt << someText << " ";
        std::string tid;
        if (aTj.AlgMod[kKilled]) { tid = "k" + std::to_string(aTj.UID); }
        else {
          tid = "t" + std::to_string(aTj.UID);
        }
        myprt << std::fixed << std::setw(5) << tid;
        myprt << std::setw(6) << aTj.CTP;
        myprt << std::setw(5) << aTj.Pass;
        myprt << std::setw(5) << aTj.EndPt[1] - aTj.EndPt[0] + 1;
        unsigned short endPt0 = aTj.EndPt[0];
        auto& tp0 = aTj.Pts[endPt0];
        int itick = tp0.Pos[1] / tcc.unitsPerTick;
        if (itick < 0) itick = 0;
        myprt << std::setw(6) << (int)(tp0.Pos[0] + 0.5) << ":" << itick; // W:T
        if (itick < 10) { myprt << " "; }
        if (itick < 100) { myprt << " "; }
        if (itick < 1000) { myprt << " "; }
        myprt << std::setw(6) << std::setprecision(2) << tp0.Ang;
        myprt << std::setw(2) << tp0.AngleCode;
        if (aTj.EndFlag[0][kBragg]) { myprt << "B"; }
        else if (aTj.EndFlag[0][kAtVtx]) {
          myprt << "V";
        }
        else if (aTj.EndFlag[0][kAtKink]) {
          myprt << "K";
        }
        else if (aTj.EndFlag[0][kAtTj]) {
          myprt << "T";
        }
        else {
          myprt << " ";
        }
        myprt << std::setw(5) << (int)tp0.AveChg;
        unsigned short endPt1 = aTj.EndPt[1];
        auto& tp1 = aTj.Pts[endPt1];
        itick = tp1.Pos[1] / tcc.unitsPerTick;
        myprt << std::setw(6) << (int)(tp1.Pos[0] + 0.5) << ":" << itick; // W:T
        if (itick < 10) { myprt << " "; }
        if (itick < 100) { myprt << " "; }
        if (itick < 1000) { myprt << " "; }
        myprt << std::setw(6) << std::setprecision(2) << tp1.Ang;
        myprt << std::setw(2) << tp1.AngleCode;
        if (aTj.EndFlag[1][kBragg]) { myprt << "B"; }
        else if (aTj.EndFlag[1][kAtVtx]) {
          myprt << "V";
        }
        else {
          myprt << " ";
        }
        myprt << std::setw(5) << (int)tp1.AveChg;
        myprt << std::setw(7) << std::setprecision(1) << aTj.TotChg / 1000;
        myprt << std::setw(7) << std::setprecision(2) << aTj.ChgRMS;
        myprt << std::setw(5) << aTj.MCSMom;
        myprt << std::setw(4) << aTj.StepDir;
        myprt << std::setw(4) << aTj.VtxID[0];
        myprt << std::setw(4) << aTj.VtxID[1];
        myprt << std::setw(5) << aTj.PDGCode;
        myprt << std::setw(5) << aTj.ParentID;
        myprt << std::setw(5) << PrimaryID(slc, aTj);
        myprt << std::setw(6) << NeutrinoPrimaryTjID(slc, aTj);
        myprt << std::setw(7) << aTj.WorkID;
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (aTj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
        myprt << "\n";
      } // ii
      return;
    } // itj > slc.tjs.size()-1

    if (itj > slc.tjs.size() - 1) return;

    auto const& aTj = slc.tjs[itj];

    mf::LogVerbatim("TC") << "Print slc.tjs[" << itj << "] Vtx[0] " << aTj.VtxID[0] << " Vtx[1] "
                          << aTj.VtxID[1];
    myprt << "AlgBits";
    for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
      if (aTj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
    myprt << "\n";

    PrintTPHeader(someText);
    if (ipt == USHRT_MAX) {
      // print all points
      for (unsigned short ii = 0; ii < aTj.Pts.size(); ++ii)
        PrintTP(someText, slc, ii, aTj.StepDir, aTj.Pass, aTj.Pts[ii]);
    }
    else {
      // print just one
      PrintTP(someText, slc, ipt, aTj.StepDir, aTj.Pass, aTj.Pts[ipt]);
    }
  } // PrintAllTraj

  //////////////////////////////////////////
  void
  PrintTrajectory(std::string someText,
                  const TCSlice& slc,
                  const Trajectory& tj,
                  unsigned short tPoint)
  {
    // prints one or all trajectory points on tj

    if (tPoint == USHRT_MAX) {
      if (tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt << someText << " ";
        myprt << "Work:   UID " << tj.UID << "    CTP " << tj.CTP << " StepDir " << tj.StepDir
              << " PDG " << tj.PDGCode << " slc.vtxs " << tj.VtxID[0] << " " << tj.VtxID[1]
              << " nPts " << tj.Pts.size() << " EndPts " << tj.EndPt[0] << " " << tj.EndPt[1];
        myprt << " MCSMom " << tj.MCSMom;
        myprt << " EndFlags " << PrintEndFlag(tj, 0) << " " << PrintEndFlag(tj, 1);
        myprt << " AlgMods:";
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
      }
      else {
        mf::LogVerbatim myprt("TC");
        myprt << someText << " ";
        myprt << "slcID " << slc.ID << " T" << tj.ID << " uT" << tj.UID << " WorkID " << tj.WorkID
              << " StepDir " << tj.StepDir << " PDG " << tj.PDGCode << " VtxID " << tj.VtxID[0]
              << " " << tj.VtxID[1] << " nPts " << tj.Pts.size() << " EndPts " << tj.EndPt[0] << " "
              << tj.EndPt[1];
        myprt << " MCSMom " << tj.MCSMom;
        myprt << " EndFlags " << PrintEndFlag(tj, 0) << " " << PrintEndFlag(tj, 1);
        myprt << " AlgMods:";
        for (unsigned short ib = 0; ib < AlgBitNames.size(); ++ib)
          if (tj.AlgMod[ib]) myprt << " " << AlgBitNames[ib];
      }
      PrintTPHeader(someText);
      for (unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt)
        PrintTP(someText, slc, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
      // See if this trajectory is a shower Tj
      if (tj.AlgMod[kShowerTj]) {
        for (unsigned short ic = 0; ic < slc.cots.size(); ++ic) {
          if (slc.cots[ic].TjIDs.empty()) continue;
          // only print out the info for the correct Tj
          if (slc.cots[ic].ShowerTjID != tj.ID) continue;
          const ShowerStruct& ss = slc.cots[ic];
          mf::LogVerbatim myprt("TC");
          myprt << "cots index " << ic << " ";
          myprt << someText << " Envelope";
          if (ss.Envelope.empty()) { myprt << " NA"; }
          else {
            for (auto& vtx : ss.Envelope)
              myprt << " " << (int)vtx[0] << ":" << (int)(vtx[1] / tcc.unitsPerTick);
          }
          myprt << " Energy " << (int)ss.Energy;
          myprt << " Area " << std::fixed << std::setprecision(1) << (int)ss.EnvelopeArea
                << " ChgDensity " << ss.ChgDensity;
          myprt << "\nInShower TjIDs";
          for (auto& tjID : ss.TjIDs) {
            myprt << " " << tjID;
          } // tjID

          myprt << "\n";
          myprt << "NearTjIDs";
          for (auto& tjID : ss.NearTjIDs) {
            myprt << " " << tjID;
          } // tjID
          myprt << "\n";
          myprt << "\n";
          myprt << "Angle " << std::fixed << std::setprecision(2) << ss.Angle << " +/- "
                << ss.AngleErr;
          myprt << " AspectRatio " << std::fixed << std::setprecision(2) << ss.AspectRatio;
          myprt << " DirectionFOM " << std::fixed << std::setprecision(2) << ss.DirectionFOM;
          if (ss.ParentID > 0) { myprt << " Parent Tj " << ss.ParentID << " FOM " << ss.ParentFOM; }
          else {
            myprt << " No parent";
          }
          myprt << " TruParentID " << ss.TruParentID << " SS3ID " << ss.SS3ID << "\n";
          if (ss.NeedsUpdate) myprt << "*********** This shower needs to be updated ***********";
          myprt << "................................................";
        } // ic
      }   // Shower Tj
    }
    else {
      // just print one traj point
      if (tPoint > tj.Pts.size() - 1) {
        mf::LogVerbatim("TC") << "Can't print non-existent traj point " << tPoint;
        return;
      }
      PrintTP(someText, slc, tPoint, tj.StepDir, tj.Pass, tj.Pts[tPoint]);
    }
  } // PrintTrajectory

  //////////////////////////////////////////
  void
  PrintTPHeader(std::string someText)
  {
    mf::LogVerbatim("TC") << someText
                          << " TRP     CTP  Ind  Stp Delta  RMS    Ang C   Err  Dir0  Dir1      Q  "
                             "  AveQ  Pull FitChi  NTPF KinkSig  Hits ";
  } // PrintTPHeader

  ////////////////////////////////////////////////
  void
  PrintTP(std::string someText,
          const TCSlice& slc,
          unsigned short ipt,
          short dir,
          unsigned short pass,
          const TrajPoint& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt << someText << " TRP" << std::fixed;
    myprt << pass;
    if (dir > 0) { myprt << "+"; }
    else {
      myprt << "-";
    }
    myprt << std::setw(6) << tp.CTP;
    myprt << std::setw(5) << ipt;
    myprt << std::setw(5) << tp.Step;
    myprt << std::setw(6) << std::setprecision(2) << tp.Delta;
    myprt << std::setw(6) << std::setprecision(2) << tp.DeltaRMS;
    myprt << std::setw(6) << std::setprecision(2) << tp.Ang;
    myprt << std::setw(2) << tp.AngleCode;
    myprt << std::setw(6) << std::setprecision(2) << tp.AngErr;
    myprt << std::setw(6) << std::setprecision(2) << tp.Dir[0];
    myprt << std::setw(6) << std::setprecision(2) << tp.Dir[1];
    myprt << std::setw(7) << (int)tp.Chg;
    myprt << std::setw(8) << (int)tp.AveChg;
    myprt << std::setw(6) << std::setprecision(1) << tp.ChgPull;
    myprt << std::setw(7) << tp.FitChi;
    myprt << std::setw(6) << tp.NTPsFit;
    myprt << std::setw(7) << std::setprecision(3) << tp.KinkSig;
    // print the hits associated with this traj point
    if (tp.Hits.size() > 16) {
      // don't print too many hits (e.g. from a shower Tj)
      myprt << " " << tp.Hits.size() << " shower hits";
    }
    else {
      for (unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        auto& hit = (*evt.allHits)[slc.slHits[iht].allHitsIndex];
        myprt << " " << hit.WireID().Wire << ":" << (int)hit.PeakTime();
        if (tp.UseHit[ii]) {
          // Distinguish used hits from nearby hits
          myprt << "_";
        }
        else {
          myprt << "x";
        }
        myprt << "T" << slc.slHits[iht].InTraj;
      } // iht
      if (tp.InPFP > 0) myprt << " inP" << tp.InPFP;
    }
    // print Environment
    if (tp.Environment.any()) myprt << " Env: " << TPEnvString(tp);
  } // PrintTP

  /////////////////////////////////////////
  std::string
  TPEnvString(const TrajPoint& tp)
  {
    // Print environment bits in human-readable format
    std::string str = "";
    for (unsigned short ib = 0; ib < 8; ++ib) {
      // There aren't any bit names for Environment_t
      if (!tp.Environment[ib]) continue;
      if (ib == kEnvNotGoodWire) str += " NoGdwire";
      if (ib == kEnvNearMuon) str += " NearMuon";
      if (ib == kEnvNearShower) str += " NearShower";
      if (ib == kEnvOverlap) str += " Overlap";
      if (ib == kEnvUnusedHits) str += " UnusedHits";
      if (ib == kEnvNearSrcHit) str += " NearSrcHit";
      if (ib == kEnvFlag) str += " Flag";
    } // ib
    return str;
  } // TPEnvironment

  /////////////////////////////////////////
  void
  PrintPFP(std::string someText, TCSlice& slc, const PFPStruct& pfp, bool printHeader)
  {
    mf::LogVerbatim myprt("TC");
    if (printHeader) {
      myprt << someText;
      myprt << "  PFP sVx  ________sPos_______ EF _______sDir______ ____sdEdx_____ eVx  "
               "________ePos_______ EF _______eDir______ ____edEdx____   Len nTp3 MCSMom ShLike? "
               "PDG Par Prim\n";
    }
    myprt << someText;
    std::string pid = "P" + std::to_string(pfp.ID);
    myprt << std::setw(5) << pid;
    // start and end stuff
    for (unsigned short end = 0; end < 2; ++end) {
      myprt << std::setw(4) << pfp.Vx3ID[end];
      myprt << std::fixed << std::right << std::setprecision(1);
      auto pos = PosAtEnd(pfp, end);
      myprt << std::setw(7) << pos[0];
      myprt << std::setw(7) << pos[1];
      myprt << std::setw(7) << pos[2];
      // print characters that encode the EndFlag
      std::string ef;
      if (pfp.EndFlag[end][kOutFV]) { ef = "O"; }
      else {
        ef = "I";
      }
      if (pfp.EndFlag[end][kBragg]) ef += "B";
      myprt << std::setw(6) << ef;
      myprt << std::fixed << std::right << std::setprecision(2);
      auto dir = DirAtEnd(pfp, end);
      myprt << std::setw(6) << dir[0];
      myprt << std::setw(6) << dir[1];
      myprt << std::setw(6) << dir[2];
      for (auto& dedx : pfp.dEdx[end]) {
        if (dedx < 50) { myprt << std::setw(5) << std::setprecision(1) << dedx; }
        else {
          myprt << std::setw(5) << std::setprecision(0) << dedx;
        }
      } // dedx
      if (pfp.dEdx[end].size() < 3) {
        for (size_t i = 0; i < 3 - pfp.dEdx[end].size(); ++i) {
          myprt << std::setw(6) << ' ';
        }
      }
    } // startend
    // global stuff
    float length = Length(pfp);
    if (length < 100) { myprt << std::setw(5) << std::setprecision(1) << length; }
    else {
      myprt << std::setw(5) << std::setprecision(0) << length;
    }
    myprt << std::setw(5) << std::setprecision(2) << pfp.TP3Ds.size();
    myprt << std::setw(7) << MCSMom(slc, pfp.TjIDs);
    myprt << std::setw(5) << IsShowerLike(slc, pfp.TjIDs);
    myprt << std::setw(5) << pfp.PDGCode;
    myprt << "      NA";
    myprt << std::setw(4) << pfp.ParentUID;
    myprt << std::setw(5) << PrimaryUID(slc, pfp);
    if (!pfp.TjIDs.empty()) {
      for (auto& tjID : pfp.TjIDs)
        myprt << " T" << tjID;
    }
    if (!pfp.DtrUIDs.empty()) {
      myprt << " dtrs";
      for (auto& dtrUID : pfp.DtrUIDs)
        myprt << " P" << dtrUID;
    }
  } // PrintPFP

  /////////////////////////////////////////
  void
  PrintPFPs(std::string someText, TCSlice& slc)
  {
    if (slc.pfps.empty()) return;

    mf::LogVerbatim myprt("TC");
    myprt << someText;
    myprt
      << "  PFP sVx  ________sPos_______  ______sDir______  ______sdEdx_____ eVx  "
         "________ePos_______  ______eDir______  ______edEdx_____ BstPln PDG TruPDG Par Prim E*P\n";
    bool printHeader = true;
    for (auto& pfp : slc.pfps) {
      PrintPFP(someText, slc, pfp, printHeader);
      printHeader = false;
    } // im

  } // PrintPFPs

  /////////////////////////////////////////
  std::string
  PrintEndFlag(const PFPStruct& pfp, unsigned short end)
  {
    if (end > 1) return "Invalid end";
    std::string tmp;
    bool first = true;
    for (unsigned short ib = 0; ib < EndFlagNames.size(); ++ib) {
      if (pfp.EndFlag[end][ib]) {
        if (first) {
          tmp = std::to_string(end) + ":" + EndFlagNames[ib];
          first = false;
        }
        else {
          tmp += "," + EndFlagNames[ib];
        }
      }
    } // ib
    if (first) tmp = " none";
    return tmp;
  } // PrintEndFlag

  /////////////////////////////////////////
  std::string
  PrintEndFlag(const Trajectory& tj, unsigned short end)
  {
    if (end > 1) return "Invalid end";
    std::string tmp;
    bool first = true;
    for (unsigned short ib = 0; ib < EndFlagNames.size(); ++ib) {
      if (tj.EndFlag[end][ib]) {
        if (first) {
          tmp = std::to_string(end) + ":" + EndFlagNames[ib];
          first = false;
        }
        else {
          tmp += "," + EndFlagNames[ib];
        }
      }
    } // ib
    return tmp;
  } // PrintEndFlag

  /////////////////////////////////////////
  std::string
  PrintHitShort(const TCHit& tch)
  {
    if (tch.allHitsIndex > (*evt.allHits).size() - 1) return "NA";
    auto& hit = (*evt.allHits)[tch.allHitsIndex];
    return std::to_string(hit.WireID().Plane) + ":" + std::to_string(hit.WireID().Wire) + ":" +
           std::to_string((int)hit.PeakTime());
  } // PrintHit

  /////////////////////////////////////////
  std::string
  PrintHit(const TCHit& tch)
  {
    if (tch.allHitsIndex > (*evt.allHits).size() - 1) return "NA";
    auto& hit = (*evt.allHits)[tch.allHitsIndex];
    return std::to_string(hit.WireID().Plane) + ":" + std::to_string(hit.WireID().Wire) + ":" +
           std::to_string((int)hit.PeakTime()) + "_" + std::to_string(tch.InTraj);
  } // PrintHit

  /////////////////////////////////////////
  std::string
  PrintPos(const TCSlice& slc, const TrajPoint& tp)
  {
    return std::to_string(DecodeCTP(tp.CTP).Plane) + ":" + PrintPos(slc, tp.Pos);
  } // PrintPos

  /////////////////////////////////////////
  std::string
  PrintPos(const TCSlice& slc, const Point2_t& pos)
  {
    unsigned int wire = 0;
    if (pos[0] > -0.4) wire = std::nearbyint(pos[0]);
    int time = std::nearbyint(pos[1] / tcc.unitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

} // namespace tca
