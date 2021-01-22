#include "larreco/RecoAlg/TCAlg/PFPUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TDecompSVD.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <limits>
#include <stdlib.h>
#include <utility>
#include <vector>

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/StepUtils.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace tca {
  /////////////////////////////////////////
  void
  StitchPFPs()
  {
    // Stitch PFParticles in different TPCs. This does serious damage to PFPStruct and should
    // only be called from TrajCluster module just before making PFParticles to put in the event
    if (slices.size() < 2) return;
    if (tcc.geom->NTPC() == 1) return;
    if (tcc.pfpStitchCuts.size() < 2) return;
    if (tcc.pfpStitchCuts[0] <= 0) return;

    bool prt = tcc.dbgStitch;

    if (prt) {
      mf::LogVerbatim myprt("TC");
      std::string fcnLabel = "SP";
      myprt << fcnLabel << " cuts " << sqrt(tcc.pfpStitchCuts[0]) << " " << tcc.pfpStitchCuts[1]
            << "\n";
      bool printHeader = true;
      for (size_t isl = 0; isl < slices.size(); ++isl) {
        if (debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if (slc.pfps.empty()) continue;
        for (auto& pfp : slc.pfps)
          PrintP(fcnLabel, myprt, pfp, printHeader);
      } // slc
    }   // prt

    // lists of pfp UIDs to stitch
    std::vector<std::vector<int>> stLists;
    for (std::size_t sl1 = 0; sl1 < slices.size() - 1; ++sl1) {
      auto& slc1 = slices[sl1];
      for (std::size_t sl2 = sl1 + 1; sl2 < slices.size(); ++sl2) {
        auto& slc2 = slices[sl2];
        // look for PFParticles in the same recob::Slice
        if (slc1.ID != slc2.ID) continue;
        for (auto& p1 : slc1.pfps) {
          if (p1.ID <= 0) continue;
          for (auto& p2 : slc2.pfps) {
            if (p2.ID <= 0) continue;
            float maxSep2 = tcc.pfpStitchCuts[0];
            float maxCth = tcc.pfpStitchCuts[1];
            bool gotit = false;
            for (unsigned short e1 = 0; e1 < 2; ++e1) {
              auto tp3d1 = EndTP3D(p1, e1);
              if(tp3d1.Flags[kTP3DBad]) continue;
              auto pos1 = tp3d1.Pos;
              // require the end to be close to a TPC boundary
              if (InsideFV(slc1, p1, e1)) continue;
              auto dir1 = tp3d1.Dir;
              for (unsigned short e2 = 0; e2 < 2; ++e2) {
                auto tp3d2 = EndTP3D(p2, e2);
                if(tp3d2.Flags[kTP3DBad]) continue;
                auto pos2 = tp3d2.Pos;
                // require the end to be close to a TPC boundary
                if (InsideFV(slc2, p2, e2)) continue;
                auto dir2 = tp3d2.Dir;
                float sep = PosSep2(pos1, pos2);
                if (sep > maxSep2) continue;
                float cth = std::abs(DotProd(dir1, dir2));
                if (cth < maxCth) continue;
                maxSep2 = sep;
                maxCth = cth;
                gotit = true;
              } // e2
            }   // e1
            if (!gotit) continue;
            if (prt) {
              mf::LogVerbatim myprt("TC");
              myprt << "Stitch slice " << slc1.ID << " P" << p1.UID << " TPC " << p1.TPCID.TPC;
              myprt << " and P" << p2.UID << " TPC " << p2.TPCID.TPC;
              myprt << " sep " << sqrt(maxSep2) << " maxCth " << maxCth;
            }
            // see if either of these are in a list
            bool added = false;
            for (auto& pm : stLists) {
              bool p1InList = (std::find(pm.begin(), pm.end(), p1.UID) != pm.end());
              bool p2InList = (std::find(pm.begin(), pm.end(), p2.UID) != pm.end());
              if (p1InList || p2InList) {
                if (p1InList) pm.push_back(p2.UID);
                if (p2InList) pm.push_back(p1.UID);
                added = true;
              }
            } // pm
            if (added) continue;
            // start a new list
            std::vector<int> tmp(2);
            tmp[0] = p1.UID;
            tmp[1] = p2.UID;
            stLists.push_back(tmp);
            break;
          } // p2
        }   // p1
      }     // sl2
    }       // sl1
    if (stLists.empty()) return;

    for (auto& stl : stLists) {
      // Find the endpoints of the stitched pfp
      float minZ = 1E6;
      std::pair<unsigned short, unsigned short> minZIndx;
      unsigned short minZEnd = 2;
      for (auto puid : stl) {
        auto slcIndex = GetSliceIndex("P", puid);
        if (slcIndex.first == USHRT_MAX) continue;
        auto& pfp = slices[slcIndex.first].pfps[slcIndex.second];
        for (unsigned short end = 0; end < 2; ++end) {
          auto pos = EndTP3D(pfp, end).Pos;
          if (pos[2] < minZ) {
            minZ = pos[2];
            minZIndx = slcIndex;
            minZEnd = end;
          }
        } // end
      }   // puid
      if (minZEnd > 1) continue;
      // preserve the pfp with the min Z position
      auto& pfp = slices[minZIndx.first].pfps[minZIndx.second];
      if (prt) mf::LogVerbatim("TC") << "SP: P" << pfp.UID;
      // add the Tjs in the other slices to it
      for (auto puid : stl) {
        if (puid == pfp.UID) continue;
        auto sIndx = GetSliceIndex("P", puid);
        if (sIndx.first == USHRT_MAX) continue;
        auto& opfp = slices[sIndx.first].pfps[sIndx.second];
        if (prt) mf::LogVerbatim("TC") << " +P" << opfp.UID;
        pfp.TjUIDs.insert(pfp.TjUIDs.end(), opfp.TjUIDs.begin(), opfp.TjUIDs.end());
        if (prt) mf::LogVerbatim();
        // Check for parents and daughters
        if (opfp.ParentUID > 0) {
          auto pSlcIndx = GetSliceIndex("P", opfp.ParentUID);
          if (pSlcIndx.first < slices.size()) {
            auto& parpfp = slices[pSlcIndx.first].pfps[pSlcIndx.second];
            std::replace(parpfp.DtrUIDs.begin(), parpfp.DtrUIDs.begin(), opfp.UID, pfp.UID);
          } // valid pSlcIndx
        }   // has a parent
        for (auto dtruid : opfp.DtrUIDs) {
          auto dSlcIndx = GetSliceIndex("P", dtruid);
          if (dSlcIndx.first < slices.size()) {
            auto& dtrpfp = slices[dSlcIndx.first].pfps[dSlcIndx.second];
            dtrpfp.ParentUID = pfp.UID;
          } // valid dSlcIndx
        }   // dtruid
        // declare it obsolete
        opfp.ID = 0;
      } // puid
    }   // stl

  } // StitchPFPs

  void
  FindPFParticles(detinfo::DetectorClocksData const& clockData,
                  detinfo::DetectorPropertiesData const& detProp,
                  TCSlice& slc)
  {
    // Match Tjs in 3D and create PFParticles

    if (tcc.match3DCuts[0] <= 0) return;
    if (!tcc.useAlg[kMakePFPTjs]) return;

    FillWireIntersections(slc);

    // clear the TP -> P assn Tjs so that all are considered
    for (auto& tj : slc.tjs) {
      for (auto& tp : tj.Pts)
        tp.InPFP = 0;
    } // tj

    bool prt = (tcc.dbgPFP && tcc.dbgSlc);

    // Match these points in 3D
    std::vector<MatchStruct> matVec;

    // iterate twice (at most), looking for 3-plane matches in long tjs in 3-plane TPCs on the
    // first iteration, 3-plane matches in short tjs on the second iteration.
    // and 2-plane matches + dead regions in 3-plane TPCs on the last iteration
    unsigned short maxNit = 2;
    if (slc.nPlanes == 2) maxNit = 1;
    if (std::nearbyint(tcc.match3DCuts[2]) == 0) maxNit = 1;
    // fill the mAllTraj vector with TPs if we aren't using SpacePoints
    for (unsigned short nit = 0; nit < maxNit; ++nit) {
      if (evt.sptHits.empty()) FillmAllTraj(detProp, slc);
      matVec.clear();
      if (slc.nPlanes == 3 && nit == 0) {
        // look for match triplets
        Match3Planes(slc, matVec);
      }
      else {
        // look for match doublets
        Match2Planes(slc, matVec);
      }
      if (matVec.empty()) continue;
      if (prt) {
        mf::LogVerbatim myprt("TC");
        myprt << "nit " << nit << " MVI  Count  Tjs in TPC " << slc.TPCID.TPC <<"\n";
        for (unsigned int indx = 0; indx < matVec.size(); ++indx) {
          auto& ms = matVec[indx];
          myprt << std::setw(5) << indx << std::setw(6) << (int)ms.Count;
          for (auto tid : ms.TjIDs)
            myprt << " T" << tid;
          // count the number of TPs in all Tjs
          float tpCnt = 0;
          for (auto tid : ms.TjIDs) {
            auto& tj = slc.tjs[tid - 1];
            tpCnt += NumPtsWithCharge(slc, tj, false);
          } // tid
          float frac = ms.Count / tpCnt;
          myprt << " matFrac " << std::fixed << std::setprecision(3) << frac;
          myprt << "\n";
        } // indx
      }   // prt
      MakePFParticles(clockData, detProp, slc, matVec, nit);
    } // nit

    // a last debug print
    if (tcc.dbgPFP && debug.MVI != UINT_MAX) {
      for (auto& pfp : slc.pfps)
        if (tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds(clockData, detProp, "FPFP", slc, pfp, -1);
    } // last debug print

    slc.mallTraj.resize(0);

    if(tcc.dbgPFP && tcc.dbgSummary) {
      // print a list of Tjs that weren't used to make PFPs
      unsigned short nprt = 0;
      for(auto& tj : slc.tjs) {
        if(tj.AlgMod[kKilled]) continue;
        float inpfp = 0;
        float npwc = 0;
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          auto& tp = tj.Pts[ipt];
          if(tp.Chg <= 0) continue;
          ++npwc;
          if(tp.InPFP > 0) ++inpfp;
        } // tp
        float fracUsed = inpfp / npwc;
        if(fracUsed > 0.9) continue;
        ++nprt;
        if(npwc > 10 || nprt < 50) {
          mf::LogVerbatim("TC")<<"FPFP: T"<<tj.ID<<" npwc "<<npwc<<" low InPFP/npwc "<<fracUsed;
        }
      } // tj
    } // tcc.dbgPFP

  } // FindPFParticles

  ////////////////////////////////////////////////
  void
  MakePFParticles(detinfo::DetectorClocksData const& clockData,
                  detinfo::DetectorPropertiesData const& detProp,
                  TCSlice& slc,
                  std::vector<MatchStruct> matVec,
                  unsigned short matVec_Iter)
  {
    // Makes PFParticles using Tjs listed in matVec
    if (matVec.empty()) return;

    bool prt = (tcc.dbgPFP && tcc.dbgSlc);

    // create a PFParticle for each valid match combination
    for (std::size_t indx = 0; indx < matVec.size(); ++indx) {
      // tone down the level of printing in ReSection
      bool foundMVI = (tcc.dbgPFP && indx == debug.MVI && matVec_Iter == debug.MVI_Iter);
      if (foundMVI) prt = true;
      auto& ms = matVec[indx];
      if (foundMVI) {
        prt = foundMVI;
        std::cout << "found MVI " << indx << " in MakePFParticles ms.Count = " << ms.Count << "\n";
      }
      // ignore dead matches
      if (ms.Count == 0) continue;
      // count the number of TPs that are available (not already 3D-matched) and used in a pfp
      float npts = 0;
      // count the number of Bragg peaks
      unsigned short nBragg = 0;
      for(std::size_t itj = 0; itj < ms.TjIDs.size(); ++itj) {
        auto& tj = slc.tjs[ms.TjIDs[itj] - 1];
        if(tj.EndFlag[0][kEndBragg] || tj.EndFlag[1][kEndBragg]) ++nBragg;
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) if(tj.Pts[ipt].InPFP == 0) ++npts;
      } // tjID
      // Create a vector of PFPs for this match so that we can split it later on if a kink is found
      std::vector<PFPStruct> pfpVec(1);
      pfpVec[0] = CreatePFP(slc);
      if(nBragg > 1) pfpVec[0].Flags[kStops] = true;
      // Define the starting set of tjs that were matched. TPs from other tjs may be added later
      pfpVec[0].TjIDs = ms.TjIDs;
      pfpVec[0].MVI = indx;
      // fill the TP3D points using the 2D trajectory points for Tjs in TjIDs. All
      // points are put in one section
      if (!MakeTP3Ds(detProp, slc, pfpVec[0], foundMVI)) {
        if (foundMVI) mf::LogVerbatim("TC") << " MakeTP3Ds failed. ";
        continue;
      }
      // fit all the points to get the general direction
      if (!FitSection(clockData, detProp, slc, pfpVec[0], 0)) {
        if (foundMVI) {
          mf::LogVerbatim("TC") << " FitSection failed\n";
          PrintTP3Ds(clockData, detProp, "Fail", slc, pfpVec[0], -1);
        }
        continue;
      }
      if (pfpVec[0].SectionFits[0].ChiDOF > 100) {
        if (foundMVI)
          mf::LogVerbatim("TC") << " crazy high ChiDOF P" << pfpVec[0].ID << " "
                                << pfpVec[0].SectionFits[0].ChiDOF;
        Recover(clockData, detProp, slc, pfpVec[0], foundMVI);
      }
      // sort the points by the distance along the general direction vector
      for(unsigned short sfi = 0; sfi < pfpVec[0].SectionFits.size(); ++sfi) {
        if (!SortSection(pfpVec[0], sfi)) {
          if (foundMVI) mf::LogVerbatim("TC") << " SortSection "<<sfi<<" failed. ";
          continue;
        }
      } // sfi
      npts = pfpVec[0].TP3Ds.size();
      if (prt) {
        auto& pfp = pfpVec[0];
        mf::LogVerbatim myprt("TC");
        myprt << " indx " << matVec_Iter << "/" << indx << " Count " << std::setw(5)
              << (int)ms.Count;
        myprt << " P" << pfpVec[0].ID;
        myprt << " ->";
        for (auto& tjid : pfp.TjIDs)
          myprt << " T" << tjid;
        myprt << " projInPlane";
        for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(detProp, slc, pfp.SectionFits[0].Pos, pfp.SectionFits[0].Dir, inCTP);
          myprt << " " << std::setprecision(2) << tp.Delta;
        } // plane
        myprt << " maxTjLen " << (int)MaxTjLen(slc, pfp.TjIDs);
        myprt << " MCSMom " << MCSMom(slc, pfp.TjIDs);
        myprt << " PDGCodeVote " << PDGCodeVote(clockData, detProp, slc, pfp);
        myprt << " nTP3Ds " << pfp.TP3Ds.size();
        myprt << " Reco3DRange "
              << Find3DRecoRange(slc, pfp, 0, (unsigned short)tcc.match3DCuts[3], 1);
      } // prt
      if (foundMVI && !pfpVec[0].AlgMod[kSmallAng3D]) { PrintTP3Ds(clockData, detProp, "FF", slc, pfpVec[0], -1); }
      for (unsigned short ip = 0; ip < pfpVec.size(); ++ip) {
        auto& pfp = pfpVec[ip];
        if(pfp.TP3Ds.empty()) continue;
        // set the end flag bits
        geo::TPCID tpcid;
        for (unsigned short end = 0; end < 2; ++end) {
          // first set them all to 0
          pfp.EndFlag[end].reset();
          auto pos = EndTP3D(pfp, end).Pos;
          if (!InsideTPC(pos, tpcid)) pfp.EndFlag[end][kEndOutFV] = true;
        } // end
        // Set kink flag and create a vertex between this pfp and the previous one that was stored
        if (ip > 0) {
          pfp.EndFlag[0][kEndKink] = true;
          Vtx3Store vx3;
          vx3.TPCID = pfp.TPCID;
          vx3.X = pfp.TP3Ds[0].Pos[0];
          vx3.Y = pfp.TP3Ds[0].Pos[1];
          vx3.Z = pfp.TP3Ds[0].Pos[2];
          // TODO: Errors, Score?
          vx3.Score = 100;
          vx3.Vx2ID.resize(slc.nPlanes);
          vx3.Wire = -2;
          vx3.ID = slc.vtx3s.size() + 1;
          vx3.Primary = false;
          ++evt.global3V_UID;
          vx3.UID = evt.global3V_UID;
          slc.vtx3s.push_back(vx3);
          pfp.Vx3ID[0] = vx3.ID;
          auto& prevPFP = slc.pfps[slc.pfps.size() - 1];
          prevPFP.Vx3ID[1] = vx3.ID;
        } // ip > 0
        // check for a valid two-plane match with a Tj in the third plane for long pfps.
        // For short pfps, it is possible that a Tj would be too short to be reconstructed
        // in the third plane.
        if (pfp.TjIDs.size() == 2 && slc.nPlanes == 3 && pfp.TP3Ds.size() > 20 &&
            !ValidTwoPlaneMatch(detProp, slc, pfp)) {
          continue;
        }
        // Skip this combination if it isn't reconstructable in 3D
        if (Find3DRecoRange(slc, pfp, 0, (unsigned short)tcc.match3DCuts[3], 1) == USHRT_MAX)
          continue;
        // See if it possible to reconstruct in more than one section
        pfp.Flags[kCanSection] = CanSection(slc, pfp);
        // Do a fit in multiple sections if the initial fit is poor
        if (pfp.SectionFits[0].ChiDOF < tcc.match3DCuts[5]) {
          // Good fit with one section
          pfp.Flags[kNeedsUpdate] = false;
        }
        else if (pfp.Flags[kCanSection]) {
          if (!ReSection(clockData, detProp, slc, pfp, foundMVI)) continue;
        } // CanSection
        if (foundMVI && !pfp.AlgMod[kSmallAng3D]) { PrintTP3Ds(clockData, detProp, "RS", slc, pfp, -1); }
        // FillGaps3D looks for gaps in the TP3Ds vector caused by broken trajectories and
        // inserts new TP3Ds if there are hits in the gaps. This search is only done in a
        // plane if the projection of the pfp results in a large angle where 2D reconstruction
        // is likely to be poor
        FillGaps3D(clockData, detProp, slc, pfp, foundMVI);
        // Check the TP3D -> TP assn, resolve conflicts and set TP -> InPFP
        if (!ReconcileTPs(slc, pfp, foundMVI)) continue;
        // Look for mis-placed 2D and 3D vertices
        ReconcileVertices(slc, pfp, foundMVI);
        // Set isGood
        for(auto& tp3d : pfp.TP3Ds) {
          if(tp3d.Flags[kTP3DBad] || tp3d.TPIndex == USHRT_MAX) continue;
          auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
          if(tp.Environment[kEnvOverlap] || tp.Environment[kEnvUnusedHits]) tp3d.Flags[kTP3DGood] = false;
        } // tp3d
        SetPFPdEdx(clockData, detProp, slc, pfp);
        // Set the direction using dE/dx
        SetDirection(clockData, detProp, slc, pfp);
        pfp.PDGCode = PDGCodeVote(clockData, detProp, slc, pfp);
        if(pfp.PDGCode == 14) std::cout<<"MPFP: store neutrino";
        if(!Store(slc, pfp)) {
          if(tcc.dbgPFP) mf::LogVerbatim("TC")<<" Store failed P"<<pfp.ID;
          break;
        }
        if(tcc.dbgPFP && pfp.MVI == debug.MVI && !pfp.AlgMod[kSmallAng3D]) PrintTP3Ds(clockData, detProp, "STORE", slc, pfp, -1);
      } // ip (iterate over split pfps)
    }   // indx (iterate over matchVec entries)
    slc.mallTraj.resize(0);
  } // MakePFParticles

  ////////////////////////////////////////////////
  bool
  ReconcileTPs(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Reconcile TP -> P assns before the pfp is stored. The TP3D -> TP is defined but
    // the TP -> P assn may not have been done. This function overwrites the TjIDs
    // vector to be the list of Tjs that contribute > 80% of their TPs to this pfp.
    // This function returns true if the assns are consistent.

    if(!tcc.useAlg[kRTPs3D]) return true;
    if(pfp.TjIDs.empty()) return false;
    if(pfp.TP3Ds.empty()) return false;
    if(pfp.ID <= 0) return false;

    //                  Tj ID, TP count
    std::vector<std::pair<int, float>> tjTPCnt;
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.Flags[kTP3DBad] || tp3d.TPIndex == USHRT_MAX) continue;
      if(tp3d.TjID <= 0) return false;
      // compare the TP3D -> TP -> P assn with the P -> TP assn
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if (tp.InPFP > 0 && tp.InPFP != pfp.ID) return false;
      // find the (Tj ID, TP count) pair in the list
      unsigned short indx = 0;
      for (indx = 0; indx < tjTPCnt.size(); ++indx)
        if (tjTPCnt[indx].first == tp3d.TjID) break;
      if (indx == tjTPCnt.size()) tjTPCnt.push_back(std::make_pair(tp3d.TjID, 0));
      ++tjTPCnt[indx].second;
      // make the TP -> P assn
      tp.InPFP = pfp.ID;
    } // tp3d

    if(prt) mf::LogVerbatim("TC")<<"RTPs3D: P" << pfp.ID;
    std::vector<int> nTjIDs;
    for (auto& tjtpcnt : tjTPCnt) {
      auto& tj = slc.tjs[tjtpcnt.first - 1];
      float npwc = NumPtsWithCharge(slc, tj, false);
      if(prt) mf::LogVerbatim("TC")<<" T"<<tj.ID<<" npwc "<<npwc<<" in PFP "<<tjtpcnt.second;
      if (tjtpcnt.second > 0.8 * npwc) nTjIDs.push_back(tjtpcnt.first);
    } // tjtpcnt
    if (prt) { mf::LogVerbatim("TC") << "RTPs3D: P" << pfp.ID << " nTjIDs " << nTjIDs.size(); }
    // TODO: is this really a failure?
//    if (nTjIDs.size() < 2) { return false; }
//    pfp.TjIDs = nTjIDs;

    return true;
  } // ReconcileTPs

  ////////////////////////////////////////////////
  void
  ReconcileTPs(TCSlice& slc)
  {
    // Reconciles TP ownership conflicts between PFParticles
    // Make a one-to-one TP -> P assn and look for one-to-many assns.
    // Note: Comparing the pulls for a TP to two different PFParticles generally results
    // in selecting the first PFParticle that was made which is not too surprising considering
    // the order in which they were created. This comparison has been commented out in favor
    // of simply keeping the old assn and removing the new one by setting IsBad true.

    if (!tcc.useAlg[kRTPs3D]) return;

    // make a list of T -> P assns
    std::vector<int> TinP;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(std::size_t ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(tp3d.Flags[kTP3DBad] || tp3d.TPIndex == USHRT_MAX) continue;
        if(tp3d.TjID <= 0) continue;
        if(std::find(TinP.begin(), TinP.end(), tp3d.TjID) == TinP.end()) TinP.push_back(tp3d.TjID);
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if (tp.InPFP > 0) {
          // an assn exists. Set the overlap bit and check consistency
          tp.Environment[kEnvOverlap] = true;
          // keep the previous assn (since it was created earlier and is more credible) and remove the new one
          tp3d.Flags[kTP3DBad] = true;
          tp3d.Flags[kTP3DGood] = false;
          tp.InPFP = 0;
        }
        else {
          // no assn exists
          tp.InPFP = pfp.ID;
        } // tp.InPFP > 0
      }   // ipt
    }     // pfp
  }       // ReconcileTPs

  /////////////////////////////////////////
  void
  MakePFPTjs(TCSlice& slc)
  {
    // This function clobbers all of the tjs that are used in TP3Ds in the pfp and replaces
    // them with new tjs that have a consistent set of TPs to prepare for putting them
    // into the event. Note that none of the Tjs are attached to 2D vertices.
    if (!tcc.useAlg[kMakePFPTjs]) return;

    // kill trajectories
    std::vector<int> killme;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.TjID <= 0) continue;
        if(std::find(killme.begin(), killme.end(), tp3d.TjID) == killme.end()) killme.push_back(tp3d.TjID);
      } // tp3d
    } // pfp
    bool prt = (tcc.dbgPFP);
    for (auto tid : killme)
      MakeTrajectoryObsolete(slc, (unsigned int)(tid - 1));

    // Make template trajectories in each plane. These will be re-used by
    // each PFParticle
    std::vector<Trajectory> ptjs(slc.nPlanes);
    // define the basic tj variables
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      ptjs[plane].Pass = 0;
      ptjs[plane].CTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      // This Tj wasn't created by stepping
      ptjs[plane].StepDir = 0;
      // It was created by this function however
      ptjs[plane].AlgMod[kMakePFPTjs] = true;
      // and is 3D matched
      ptjs[plane].AlgMod[kMat3D] = true;
    } // plane

    // now make the new Tjs
    for (auto& pfp : slc.pfps) {
      if (pfp.ID <= 0) continue;
      // initialize the tjs
      for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        ptjs[plane].Pts.clear();
        --evt.WorkID;
        if (evt.WorkID == INT_MIN) evt.WorkID = -1;
        ptjs[plane].ID = evt.WorkID;
      } // plane
      pfp.TjIDs.clear();
      // iterate through all of the TP3Ds, adding TPs to the TJ in the appropriate plane.
      // The assumption here is that TP order reflects the TP3D order
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.TjID <= 0) continue;
        if(tp3d.Flags[kTP3DBad] || tp3d.TPIndex == USHRT_MAX) continue;
        // make a copy of the 2D TP
        auto tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if (tp.InPFP > 0 && tp.InPFP != pfp.ID) continue;
        tp.InPFP = pfp.ID;
        // the TP Step isn't useful anymore, so stash the original TJ ID into it
        tp.Step = tp3d.TjID;
        unsigned short plane = DecodeCTP(tp.CTP).Plane;
        // append it to Pts
        ptjs[plane].Pts.push_back(tp);
      } // tp3d
      // finish defining each of the Tjs and store them
      // new tj ID indexed by plane
      std::vector<int> tids(ptjs.size(), 0);
      for(unsigned short plane = 0; plane < ptjs.size(); ++plane) {
        auto& tj = ptjs[plane];
        if(tj.Pts.size() < 2) continue;
        SetEndPoints(tj);
        tj.PDGCode = pfp.PDGCode;
        tj.MCSMom = MCSMom(slc, tj);
        if (!StoreTraj(slc, tj)) continue;
        // associate it with the pfp
        auto& newTj = slc.tjs.back();
        pfp.TjIDs.push_back(newTj.ID);
        tids[plane] = newTj.ID;
      } // tj
      // preserve the PFP -> 3V -> 2V -> T assns
      for(unsigned short end = 0; end < 2; ++end) {
        if(pfp.Vx3ID[end] <= 0) continue;
        auto& vx3 = slc.vtx3s[pfp.Vx3ID[end] - 1];
        for(unsigned short plane = 0; plane < ptjs.size(); ++plane) {
          if(tids[plane] == 0) continue;
          if(vx3.Vx2ID[plane] <= 0) continue;
          auto& vx2 = slc.vtxs[vx3.Vx2ID[plane] - 1];
          auto& tj = slc.tjs[tids[plane] - 1];
          auto tend = CloseEnd(slc, tj, vx2.Pos);
          tj.VtxID[tend] = vx2.ID;
          if(prt) mf::LogVerbatim("TC") << "MPFPTjs: 3V" << vx3.ID << " -> 2V" << vx2.ID
                   << " -> T" << tj.ID << "_" << tend << " in plane " << plane;
        } // plane
      } // end
    }   // pfp
  }     // MakePFPTjs

  /////////////////////////////////////////
  void
  FillWireIntersections(TCSlice& slc)
  {
    // Find wire intersections and put them in evt.wireIntersections

    // see if anything needs to be done
    if (!evt.wireIntersections.empty() && evt.wireIntersections[0].tpc == slc.TPCID.TPC) return;

    evt.wireIntersections.clear();

    unsigned int cstat = slc.TPCID.Cryostat;
    unsigned int tpc = slc.TPCID.TPC;
    // find the minMax number of wires in each plane of the TPC
    unsigned int maxWire = slc.nWires[0];
    for (auto nw : slc.nWires)
      if (nw < maxWire) maxWire = nw;
    // Start looking for intersections in the middle
    unsigned int firstWire = maxWire / 2;

    // find a valid wire intersection in all plane combinations
    std::vector<std::pair<unsigned short, unsigned short>> pln1pln2;
    for (unsigned short pln1 = 0; pln1 < slc.nPlanes - 1; ++pln1) {
      for (unsigned short pln2 = pln1 + 1; pln2 < slc.nPlanes; ++pln2) {
        auto p1p2 = std::make_pair(pln1, pln2);
        if (std::find(pln1pln2.begin(), pln1pln2.end(), p1p2) != pln1pln2.end()) continue;
        // find two wires that have a valid intersection
        for (unsigned int wire = firstWire; wire < maxWire; ++wire) {
          double y00, z00;
          if (!tcc.geom->IntersectionPoint(wire, wire, pln1, pln2, cstat, tpc, y00, z00)) continue;
          // increment by one wire in pln1 and find another valid intersection
          double y10, z10;
          if (!tcc.geom->IntersectionPoint(wire + 10, wire, pln1, pln2, cstat, tpc, y10, z10))
            continue;
          // increment by one wire in pln2 and find another valid intersection
          double y01, z01;
          if (!tcc.geom->IntersectionPoint(wire, wire + 10, pln1, pln2, cstat, tpc, y01, z01))
            continue;
          TCWireIntersection tcwi;
          tcwi.tpc = tpc;
          tcwi.pln1 = pln1;
          tcwi.pln2 = pln2;
          tcwi.wir1 = wire;
          tcwi.wir2 = wire;
          tcwi.y = y00;
          tcwi.z = z00;
          tcwi.dydw1 = (y10 - y00) / 10;
          tcwi.dzdw1 = (z10 - z00) / 10;
          tcwi.dydw2 = (y01 - y00) / 10;
          tcwi.dzdw2 = (z01 - z00) / 10;
          evt.wireIntersections.push_back(tcwi);
          break;
        } // wire
      }   // pln2
    }     // pln1
  }       // FillWireIntersections

  /////////////////////////////////////////
  bool
  TCIntersectionPoint(unsigned int wir1,
                      unsigned int wir2,
                      unsigned int pln1,
                      unsigned int pln2,
                      float& y,
                      float& z)
  {
    // A TrajCluster analog of geometry IntersectionPoint that uses local wireIntersections with
    // float precision. The (y,z) position is only used to match TPs between planes - not for 3D fitting
    if (evt.wireIntersections.empty()) return false;
    if (pln1 == pln2) return false;

    if (pln1 > pln2) {
      std::swap(pln1, pln2);
      std::swap(wir1, wir2);
    }

    for (auto& wi : evt.wireIntersections) {
      if (wi.pln1 != pln1) continue;
      if (wi.pln2 != pln2) continue;
      // estimate the position using the wire differences
      double dw1 = wir1 - wi.wir1;
      double dw2 = wir2 - wi.wir2;
      y = (float)(wi.y + dw1 * wi.dydw1 + dw2 * wi.dydw2);
      z = (float)(wi.z + dw1 * wi.dzdw1 + dw2 * wi.dzdw2);
      return true;
    } // wi
    return false;
  } // TCIntersectionPoint

  /////////////////////////////////////////
  void
  Match3PlanesSpt(TCSlice& slc, std::vector<MatchStruct>& matVec)
  {
    // fill matVec using SpacePoint -> Hit -> TP -> tj assns
    if (evt.sptHits.empty()) return;

    // create a local vector of allHit -> Tj assns and populate it
    std::vector<int> inTraj((*evt.allHits).size(), 0);
    for (auto& tch : slc.slHits)
      inTraj[tch.allHitsIndex] = tch.InTraj;

    // the TJ IDs for one match
    std::array<int, 3> tIDs;
    // vector for matched Tjs
    std::vector<std::array<int, 3>> mtIDs;
    // and a matching vector for the count
    std::vector<unsigned short> mCnt;
    // ignore Tj matches after hitting a user-defined limit
    unsigned short maxCnt = USHRT_MAX;
    if (tcc.match3DCuts[1] < (float)USHRT_MAX) maxCnt = (unsigned short)tcc.match3DCuts[1];
    // a list of those Tjs
    std::vector<unsigned short> tMaxed;

    unsigned int tpc = slc.TPCID.TPC;

    for (auto& sptHits : evt.sptHits) {
      if (sptHits.size() != 3) continue;
      // ensure that the SpacePoint is in the requested TPC
      if (!SptInTPC(sptHits, tpc)) continue;
      unsigned short cnt = 0;
      for (unsigned short plane = 0; plane < 3; ++plane) {
        unsigned int iht = sptHits[plane];
        if (iht == UINT_MAX) continue;
        if (inTraj[iht] <= 0) continue;
        tIDs[plane] = inTraj[iht];
        ++cnt;
      } // iht
      if (cnt != 3) continue;
      // look for it in the list of tj combinations
      unsigned short indx = 0;
      for (indx = 0; indx < mtIDs.size(); ++indx)
        if (tIDs == mtIDs[indx]) break;
      if (indx == mtIDs.size()) {
        // not found so add it to mtIDs and add another element to mCnt
        mtIDs.push_back(tIDs);
        mCnt.push_back(0);
      }
      ++mCnt[indx];
      if (mCnt[indx] == maxCnt) {
        // add the Tjs to the list
        tMaxed.insert(tMaxed.end(), tIDs[0]);
        tMaxed.insert(tMaxed.end(), tIDs[1]);
        tMaxed.insert(tMaxed.end(), tIDs[2]);
        break;
      } // hit maxCnt
      ++cnt;
    } // sptHit

    std::vector<SortEntry> sortVec;
    for (unsigned short indx = 0; indx < mCnt.size(); ++indx) {
      auto& tIDs = mtIDs[indx];
      // find the fraction of TPs on the shortest tj that are matched
      float minTPCnt = USHRT_MAX;
      for (auto tid : tIDs) {
        auto& tj = slc.tjs[tid - 1];
        float tpcnt = NumPtsWithCharge(slc, tj, false);
        if (tpcnt < minTPCnt) minTPCnt = tpcnt;
      } // tid
      float frac = (float)mCnt[indx] / minTPCnt;
      // ignore matches with a very low match fraction
      if (frac < 0.05) continue;
      SortEntry se;
      se.index = indx;
      se.val = mCnt[indx];
      sortVec.push_back(se);
    } // ii
    if (sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valsDecreasing);

    matVec.resize(sortVec.size());

    for (unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      auto& ms = matVec[ii];
      ms.Count = mCnt[indx];
      ms.TjIDs.resize(3);
      for (unsigned short plane = 0; plane < 3; ++plane)
        ms.TjIDs[plane] = mtIDs[indx][plane];
    } // indx

  } // Match3PlanesSpt

  /////////////////////////////////////////
  bool
  SptInTPC(const std::array<unsigned int, 3>& sptHits, unsigned int tpc)
  {
    // returns true if a hit referenced in sptHits resides in the requested tpc. We assume
    // that if one does, then all of them do

    unsigned int ahi = UINT_MAX;
    for (auto ii : sptHits)
      if (ii != UINT_MAX) {
        ahi = ii;
        break;
      }
    if (ahi >= (*evt.allHits).size()) return false;
    // get a reference to the hit and see if it is in the desired tpc
    auto& hit = (*evt.allHits)[ahi];
    if (hit.WireID().TPC == tpc) return true;
    return false;

  } // SptInTPC

  /////////////////////////////////////////
  void
  Match3Planes(TCSlice& slc, std::vector<MatchStruct>& matVec)
  {
    // A simpler and faster version of MatchPlanes that only creates three plane matches

    if (slc.nPlanes != 3) return;

    // use SpacePoint -> Hit -> TP assns?
    if (!evt.sptHits.empty()) {
      Match3PlanesSpt(slc, matVec);
      return;
    }

    if (slc.mallTraj.empty()) return;
    float xcut = tcc.match3DCuts[0];
    double yzcut = 1.5 * tcc.wirePitch;

    // the TJ IDs for one match
    std::array<unsigned short, 3> tIDs;
    // vector for matched Tjs
    std::vector<std::array<unsigned short, 3>> mtIDs;
    // and a matching vector for the count
    std::vector<unsigned short> mCnt;
    // ignore Tj matches after hitting a user-defined limit
    unsigned short maxCnt = USHRT_MAX;
    if (tcc.match3DCuts[1] < (float)USHRT_MAX) maxCnt = (unsigned short)tcc.match3DCuts[1];
    // a list of those Tjs
    std::vector<unsigned short> tMaxed;

    for (std::size_t ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // see if we hit the maxCnt limit
      if (std::find(tMaxed.begin(), tMaxed.end(), iTjPt.id) != tMaxed.end()) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      if(tcc.useAlg[kNewCuts] && itp.InPFP > 0) continue;
      unsigned int iPlane = iTjPt.plane;
      unsigned int iWire = std::nearbyint(itp.Pos[0]);
      tIDs[iPlane] = iTjPt.id;
      bool hitMaxCnt = false;
      for (std::size_t jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if (jTjPt.plane == iTjPt.plane) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if (jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large
        if (jTjPt.xlo > iTjPt.xhi + xcut) break;
        // see if we hit the maxCnt limit
        if (std::find(tMaxed.begin(), tMaxed.end(), jTjPt.id) != tMaxed.end()) continue;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        if(tcc.useAlg[kNewCuts] && jtp.InPFP > 0) continue;
        unsigned short jPlane = jTjPt.plane;
        unsigned int jWire = jtp.Pos[0];
        Point2_t ijPos;
        if (!TCIntersectionPoint(iWire, jWire, iPlane, jPlane, ijPos[0], ijPos[1])) continue;
        tIDs[jPlane] = jTjPt.id;
        for (std::size_t kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
          auto& kTjPt = slc.mallTraj[kpt];
          // ensure that the planes are different
          if (kTjPt.plane == iTjPt.plane || kTjPt.plane == jTjPt.plane) continue;
          if (kTjPt.xlo > iTjPt.xhi) continue;
          // break out if the x range difference becomes large
          if (kTjPt.xlo > iTjPt.xhi + xcut) break;
          // see if we hit the maxCnt limit
          if (std::find(tMaxed.begin(), tMaxed.end(), kTjPt.id) != tMaxed.end()) continue;
          auto& ktp = slc.tjs[kTjPt.id - 1].Pts[kTjPt.ipt];
          if(tcc.useAlg[kNewCuts] && ktp.InPFP > 0) continue;
          unsigned short kPlane = kTjPt.plane;
          unsigned int kWire = ktp.Pos[0];
          Point2_t ikPos;
          if (!TCIntersectionPoint(iWire, kWire, iPlane, kPlane, ikPos[0], ikPos[1])) continue;
          if (std::abs(ijPos[0] - ikPos[0]) > yzcut) continue;
          if (std::abs(ijPos[1] - ikPos[1]) > yzcut) continue;
          // we have a match
          tIDs[kPlane] = kTjPt.id;
          // look for it in the list
          unsigned int indx = 0;
          for (indx = 0; indx < mtIDs.size(); ++indx)
            if (tIDs == mtIDs[indx]) break;
          if (indx == mtIDs.size()) {
            // not found so add it to mtIDs and add another element to mCnt
            mtIDs.push_back(tIDs);
            mCnt.push_back(0);
          }
          ++mCnt[indx];
          if (mCnt[indx] == maxCnt) {
            // add the Tjs to the list
            tMaxed.insert(tMaxed.end(), tIDs[0]);
            tMaxed.insert(tMaxed.end(), tIDs[1]);
            tMaxed.insert(tMaxed.end(), tIDs[2]);
            hitMaxCnt = true;
            break;
          } // hit maxCnt
        }   // kpt
        if (hitMaxCnt) break;
      } // jpt
    }   // ipt

    if (mCnt.empty()) return;

    std::vector<SortEntry> sortVec;
    for (std::size_t indx = 0; indx < mCnt.size(); ++indx) {
      auto& tIDs = mtIDs[indx];
      // count the number of TPs in all Tjs
      float tpCnt = 0;
      for (auto tid : tIDs) {
        auto& tj = slc.tjs[tid - 1];
        tpCnt += NumPtsWithCharge(slc, tj, false);
      } // tid
      float frac = mCnt[indx] / tpCnt;
      frac /= 3;
      // ignore matches with a very low match fraction
      if (frac < 0.05) continue;
      SortEntry se;
      se.index = indx;
      se.val = mCnt[indx];
      sortVec.push_back(se);
    } // ii
    if (sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valsDecreasing);

    matVec.resize(sortVec.size());

    for (std::size_t ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      auto& ms = matVec[ii];
      ms.Count = mCnt[indx];
      ms.TjIDs.resize(3);
      for (unsigned short plane = 0; plane < 3; ++plane)
        ms.TjIDs[plane] = (int)mtIDs[indx][plane];
    } // indx

  } // Match3Planes

  /////////////////////////////////////////
  void
  Match2Planes(TCSlice& slc, std::vector<MatchStruct>& matVec)
  {
    // A simpler faster version of MatchPlanes that only creates two plane matches

    matVec.clear();
    if (slc.mallTraj.empty()) return;

    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;

    float xcut = tcc.match3DCuts[0];

    // the TJ IDs for one match
    std::array<unsigned short, 2> tIDs;
    // vector for matched Tjs
    std::vector<std::array<unsigned short, 2>> mtIDs;
    // and a matching vector for the count
    std::vector<unsigned short> mCnt;
    // ignore Tj matches after hitting a user-defined limit
    unsigned short maxCnt = USHRT_MAX;
    if (tcc.match3DCuts[1] < (float)USHRT_MAX) maxCnt = (unsigned short)tcc.match3DCuts[1];
    // a list of those Tjs
    std::vector<unsigned short> tMaxed;

    for (std::size_t ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // see if we hit the maxCnt limit
      if (std::find(tMaxed.begin(), tMaxed.end(), iTjPt.id) != tMaxed.end()) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      if(tcc.useAlg[kNewCuts] && itp.InPFP > 0) continue;
      unsigned short iPlane = iTjPt.plane;
      unsigned int iWire = itp.Pos[0];
      bool hitMaxCnt = false;
      for (std::size_t jpt = ipt + 1; jpt < slc.mallTraj.size(); ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if (jTjPt.plane == iTjPt.plane) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if (jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large
        if (jTjPt.xlo > iTjPt.xhi + xcut) break;
        // see if we hit the maxCnt limit
        if (std::find(tMaxed.begin(), tMaxed.end(), jTjPt.id) != tMaxed.end()) continue;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        if(tcc.useAlg[kNewCuts] && jtp.InPFP > 0) continue;
        unsigned short jPlane = jTjPt.plane;
        unsigned int jWire = jtp.Pos[0];
        Point3_t ijPos;
        ijPos[0] = itp.Pos[0];
        if (!tcc.geom->IntersectionPoint(
              iWire, jWire, iPlane, jPlane, cstat, tpc, ijPos[1], ijPos[2]))
          continue;
        tIDs[0] = iTjPt.id;
        tIDs[1] = jTjPt.id;
        // swap the order so that the == operator works correctly
        if (tIDs[0] > tIDs[1]) std::swap(tIDs[0], tIDs[1]);
        // look for it in the list
        std::size_t indx = 0;
        for (indx = 0; indx < mtIDs.size(); ++indx)
          if (tIDs == mtIDs[indx]) break;
        if (indx == mtIDs.size()) {
          // not found so add it to mtIDs and add another element to mCnt
          mtIDs.push_back(tIDs);
          mCnt.push_back(0);
        }
        ++mCnt[indx];
        if (mCnt[indx] == maxCnt) {
          // add the Tjs to the list
          tMaxed.insert(tMaxed.end(), tIDs[0]);
          tMaxed.insert(tMaxed.end(), tIDs[1]);
          hitMaxCnt = true;
          break;
        } // hit maxCnt
        if (hitMaxCnt) break;
      } // jpt
    }   // ipt

    if (mCnt.empty()) return;

    std::vector<SortEntry> sortVec;
    for (std::size_t indx = 0; indx < mCnt.size(); ++indx) {
      auto& tIDs = mtIDs[indx];
      // count the number of TPs in all Tjs
      float tpCnt = 0;
      for (auto tid : tIDs) {
        auto& tj = slc.tjs[tid - 1];
        tpCnt += NumPtsWithCharge(slc, tj, false);
      } // tid
      float frac = mCnt[indx] / tpCnt;
      frac /= 2;
      // ignore matches with a very low match fraction
      if (frac < 0.05) continue;
      SortEntry se;
      se.index = indx;
      se.val = mCnt[indx];
      sortVec.push_back(se);
    } // ii
    if (sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valsDecreasing);

    matVec.resize(sortVec.size());

    for (std::size_t ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      auto& ms = matVec[ii];
      ms.Count = mCnt[indx];
      ms.TjIDs.resize(2);
      for (unsigned short plane = 0; plane < 2; ++plane)
        ms.TjIDs[plane] = (int)mtIDs[indx][plane];
    } // indx

  } // Match2Planes

  /////////////////////////////////////////
  bool
  Update(detinfo::DetectorClocksData const& clockData,
         detinfo::DetectorPropertiesData const& detProp,
         const TCSlice& slc,
         PFPStruct& pfp,
         bool prt)
  {
    // This function only updates SectionFits that need to be re-sorted or re-fit. It returns
    // false if there was a serious error indicating that the pfp should be abandoned
    if (pfp.TP3Ds.empty() || pfp.SectionFits.empty()) return false;

    if(!pfp.Flags[kdEdxDefined]) SetTP3DdEdx(clockData, detProp, slc, pfp);

    // special handling for small angle tracks
    if(pfp.AlgMod[kSmallAng3D]) {
      for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        auto& sf = pfp.SectionFits[sfi];
        if (!sf.NeedsUpdate) continue;
        if (!SortSection(pfp, sfi)) return false;
        sf.NPts = 0;
        sf.ChiDOF = 0;
        for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
          auto& tp3d = pfp.TP3Ds[ipt];
          if(tp3d.SFIndex < sfi) continue;
          if(tp3d.SFIndex > sfi) break;
          if(tp3d.Flags[kTP3DBad]) continue;
          ++sf.NPts;
          double weight = 0.2;
          // de-weight points that don't have the expected dE/dx
          if(tp3d.dEdx <= 0) tp3d.dEdx = dEdx(clockData, detProp, slc, tp3d);
          if(tp3d.dEdx > 0) {
            weight = 1;
            if(tp3d.dEdx < 1.5) weight = 0.2;
            if(tp3d.Flags[kTP3DHiDEdx]) weight = 0.2;
          } // valid dEdx
          weight /= tp3d.TPXErr2;
          double delta = tp3d.Pos[0] - tp3d.TPX;
          sf.ChiDOF += weight * delta * delta;
        } // ipt
        if(sf.NPts < 5) {
          sf.ChiDOF = 0;
        } else {
          sf.ChiDOF /= (float)(sf.NPts - 4);
        }
        sf.NeedsUpdate = false;
      } // sfi
      pfp.Flags[kNeedsUpdate] = false;
      return true;
    } // kSmallAng3D

    for (unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      if (!sf.NeedsUpdate) continue;
      if (!FitSection(clockData, detProp, slc, pfp, sfi)) return false;
      if (!SortSection(pfp, sfi)) return false;
      sf.NeedsUpdate = false;
    } // sfi

    // ensure that all points (good or not) have a valid SFIndex
    for (auto& tp3d : pfp.TP3Ds) {
      if (tp3d.SFIndex >= pfp.SectionFits.size()) SetSection(detProp, slc, pfp, tp3d);
    } // tp3d
    pfp.Flags[kNeedsUpdate] = false;
    return true;
  } // Update

  /////////////////////////////////////////
  bool
  ReSection(detinfo::DetectorClocksData const& clockData,
            detinfo::DetectorPropertiesData const& detProp,
            const TCSlice& slc,
            PFPStruct& pfp,
            bool prt)
  {
    // Re-fit the TP3Ds in sections and add/remove sections to keep ChiDOF of each section close to 1.
    // This function only fails when there is a serious error, otherwise if reasonable fits cannot be
    // achieved, the CanSection flag is set false.
    if (pfp.SectionFits.empty()) return false;
    // This function shouldn't be called if this is the case but it isn't a major failure if it is
    if(!pfp.Flags[kCanSection]) return true;
    if(pfp.AlgMod[kSmallAng3D]) return true;
    // Likewise this shouldn't be attempted if there aren't at least 3 points in 2 planes in 2 sections
    // but it isn't a failure
    if (pfp.TP3Ds.size() < 12) {
      pfp.Flags[kCanSection] = false;
      return true;
    }

    prt = (pfp.MVI == debug.MVI);

    // try to keep ChiDOF between chiLo and chiHi
    float chiLo = 0.5 * tcc.match3DCuts[5];
    float chiHi = 1.5 * tcc.match3DCuts[5];
    if(prt) mf::LogVerbatim("TC")<<"RS: Section chisq target range "<<chiLo<<" - "<<chiHi;

    // clobber the old sections if more than one exists
    if (pfp.SectionFits.size() > 1) {
      // make one section
      pfp.SectionFits.resize(1);
      // put all of the points in it and fit
      for (auto& tp3d : pfp.TP3Ds) {
        tp3d.SFIndex = 0;
        tp3d.Flags[kTP3DGood] = true;
      }
      auto& sf = pfp.SectionFits[0];
      if (!FitSection(clockData, detProp, slc, pfp, 0)) { return false; }
      if (sf.ChiDOF < tcc.match3DCuts[5]) return true;
    } // > 1 SectionFit
    // sort by distance from the start
    if (!SortSection(pfp, 0)) return false;
    // require a minimum of 3 points in 2 planes
    unsigned short min2DPts = 3;
    unsigned short fromPt = 0;
    // set the section index to invalid for all points
    for (auto& tp3d : pfp.TP3Ds)
      tp3d.SFIndex = USHRT_MAX;
    // Guess how many points should be added in each iteration
    unsigned short nPtsToAdd = pfp.TP3Ds.size() / 4;
    // the actual number of points that will be fit in the section
    unsigned short nPts = nPtsToAdd;
    // the minimum number of points
    unsigned short nPtsMin = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1) - fromPt + 1;
    if (nPtsMin >= pfp.TP3Ds.size()) {
      pfp.Flags[kCanSection] = false;
      return true;
    }
    float chiDOF = 0;
    if (nPts < nPtsMin) nPts = nPtsMin;
    // Try to reduce the number of iterations for long pfps
    if (pfp.TP3Ds.size() > 100) {
      unsigned short nhalf = pfp.TP3Ds.size() / 2;
      FitPFP(clockData, detProp, slc, pfp, fromPt, nhalf, USHRT_MAX, chiDOF);
      if (chiDOF < tcc.match3DCuts[5]) nPts = nhalf;
    }
    bool lastSection = false;
    for (unsigned short sfIndex = 0; sfIndex < 20; ++sfIndex) {
      // Try to add/remove points in each section no more than 20 times
      float chiDOFPrev = FLT_MAX;
      short nHiChi = 0;
      unsigned short maxNit = 10;
      if(pfp.Flags[kStops]) maxNit = 1;
      for(unsigned short nit = 0; nit < maxNit; ++nit) {
        // Decide how many points to add or subtract after doing the fit
        unsigned short nPtsNext = nPts;
        if (!FitPFP(clockData, detProp, slc, pfp, fromPt, nPts, USHRT_MAX, chiDOF)) {
          nPtsNext += 1.5 * nPtsToAdd;
        }
        else if (chiDOF < chiLo) {
          // low chiDOF
          if (nHiChi > 2) {
            // declare it close enough if several attempts were made
            nPtsNext = 0;
          }
          else {
            nPtsNext += nPtsToAdd;
          } // nHiChi < 2
          nHiChi = 0;
        }
        else if (chiDOF > chiHi) {
          // high chiDOF
          ++nHiChi;
          if (nHiChi == 1 && chiDOFPrev > tcc.match3DCuts[5]) {
            // reduce the number of points by 1/2 on the first attempt
            nPtsNext /= 2;
          }
          else {
            // that didn't work so start subtracting groups of points
            short npnext = (short)nPts - nHiChi * 5;
            // assume this won't work
            nPtsNext = 0;
            if (npnext > nPtsMin) nPtsNext = npnext;
          }
        }
        else {
          // just right
          nPtsNext = 0;
        }
        // check for passing the end
        if (fromPt + nPtsNext >= pfp.TP3Ds.size()) {
          nPtsNext = pfp.TP3Ds.size() - fromPt;
          lastSection = true;
        }
        if (prt) {
          mf::LogVerbatim myprt("TC");
          myprt << " RS: P" << pfp.ID << " sfi/nit/npts " << sfIndex << "/" << nit << "/" << nPts;
          myprt << std::fixed << std::setprecision(2) << " chiDOF " << chiDOF;
          myprt << " fromPt " << fromPt;
          myprt << " nPtsNext " << nPtsNext;
          myprt << " nHiChi " << nHiChi;
          myprt << " lastSection? " << lastSection;
        }
        if (nPtsNext == 0) break;
        // see if this is the last section
        if (lastSection) break;
        if (chiDOF == chiDOFPrev) {
          if (prt) mf::LogVerbatim("TC") << " MVI " << pfp.MVI << " chiDOF "<<chiDOF<<" not changing";
          break;
        }
        nPts = nPtsNext;
        chiDOFPrev = chiDOF;
      } // nit
      // finished this section. Assign the points to it
      unsigned short toPt = fromPt + nPts;
      if (toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
      for (unsigned short ipt = fromPt; ipt < toPt; ++ipt)
        pfp.TP3Ds[ipt].SFIndex = sfIndex;
      // See if there are enough points remaining to reconstruct another section if this isn't known
      // to be the last section
      if (!lastSection) {
        // this will be the first point in the next section
        unsigned short nextFromPt = fromPt + nPts;
        // See if it will have enough points to be reconstructed
        unsigned short nextToPtMin = Find3DRecoRange(slc, pfp, nextFromPt, min2DPts, 1);
        if (nextToPtMin == USHRT_MAX) {
          // not enough points so this is the last section
          lastSection = true;
          // assign the remaining points to the last section
          for (std::size_t ipt = nextFromPt; ipt < pfp.TP3Ds.size(); ++ipt)
            pfp.TP3Ds[ipt].SFIndex = sfIndex;
        }
      } // !lastSection
      // Do a final fit and update the points. Don't worry about a poor ChiDOF
      FitSection(clockData, detProp, slc, pfp, sfIndex);
      if (!SortSection(pfp, 0)) { return false; }
      if (lastSection) break;
      // Prepare for the next section.
      fromPt = fromPt + nPts;
      nPts = nPtsToAdd;
      nPtsMin = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1) - fromPt + 1;
      if (nPtsMin >= pfp.TP3Ds.size()) break;
      // add a new section
      pfp.SectionFits.resize(pfp.SectionFits.size() + 1);
    } // snit

    // see if the last sf is valid
    if (pfp.SectionFits.size() > 1 && pfp.SectionFits.back().ChiDOF < 0) {
      unsigned short badSFI = pfp.SectionFits.size() - 1;
      // remove it
      pfp.SectionFits.pop_back();
      for (std::size_t ipt = pfp.TP3Ds.size() - 1; ipt > 0; --ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if (tp3d.SFIndex < badSFI) break;
        --tp3d.SFIndex;
      }
      pfp.SectionFits.back().NeedsUpdate = true;
    } // bad last SF

    // Ensure that the points at the end are in the last section
    for (std::size_t ipt = pfp.TP3Ds.size() - 1; ipt > 0; --ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if (tp3d.SFIndex < pfp.SectionFits.size()) break;
      tp3d.SFIndex = pfp.SectionFits.size() - 1;
      pfp.Flags[kNeedsUpdate] = true;
      pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
    } // tp3d

    Update(clockData, detProp, slc, pfp, prt);

    // set CanSection false if the chisq is poor in any section
    for (auto& sf : pfp.SectionFits) {
      if (sf.ChiDOF > tcc.match3DCuts[5]) pfp.Flags[kCanSection] = false;
    }

    return true;
  } // resection

  /////////////////////////////////////////
  bool
  CanSection(const TCSlice& slc, const PFPStruct& pfp)
  {
    // analyze the TP3D vector to determine if it can be reconstructed in 3D in more than one section with
    // the requirement that there are at least 3 points in two planes
    if(pfp.AlgMod[kJunk3D]) return false;
    if(pfp.AlgMod[kSmallAng3D]) return false;
    if(pfp.TP3Ds.size() < 12) return false;
    unsigned short toPt = Find3DRecoRange(slc, pfp, 0, 3, 1);
    if (toPt > pfp.TP3Ds.size()) return false;
    unsigned short nextToPt = Find3DRecoRange(slc, pfp, toPt, 3, 1);
    if (nextToPt > pfp.TP3Ds.size()) return false;
    return true;
  } // CanSection

  /////////////////////////////////////////
  unsigned short
  Find3DRecoRange(const TCSlice& slc,
                  const PFPStruct& pfp,
                  unsigned short fromPt,
                  unsigned short min2DPts,
                  short dir)
  {
    // Scans the TP3Ds vector starting at fromPt until it finds min2DPts in two planes. It returns
    // with the index of that point (+1) in the TP3Ds vector. The dir variable defines the scan direction in
    // the TP3Ds vector
    if (fromPt > pfp.TP3Ds.size() - 1) return USHRT_MAX;
    if (pfp.TP3Ds.size() < 2 * min2DPts) return USHRT_MAX;
    if (dir == 0) return USHRT_MAX;

    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for (std::size_t ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      unsigned short ipt = fromPt + ii;
      if (dir < 0) ipt = fromPt - ii;
      if (ipt >= pfp.TP3Ds.size()) break;
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.Flags[kTP3DGood]) continue;
      unsigned short plane = DecodeCTP(slc.tjs[tp3d.TjID - 1].CTP).Plane;
      ++cntInPln[plane];
      unsigned short enufInPlane = 0;
      for (unsigned short plane = 0; plane < slc.nPlanes; ++plane)
        if (cntInPln[plane] >= min2DPts) ++enufInPlane;
      if (enufInPlane > 1) return ipt + 1;
      if (dir < 0 && ipt == 0) break;
    } // ipt
    return USHRT_MAX;
  } // Find3DRecoRange

  /////////////////////////////////////////
  void
  GetRange(const PFPStruct& pfp,
           unsigned short sfIndex,
           unsigned short& fromPt,
           unsigned short& npts)
  {
    fromPt = USHRT_MAX;
    if (sfIndex >= pfp.SectionFits.size()) return;
    if (pfp.TP3Ds.empty()) return;
    fromPt = USHRT_MAX;
    npts = 0;
    // Note that no test is made for not-good TP3Ds here since that would give a wrong npts count
    for (std::size_t ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if (tp3d.SFIndex < sfIndex) continue;
      if (tp3d.SFIndex > sfIndex) break;
      if (fromPt == USHRT_MAX) fromPt = ipt;
      ++npts;
    } // ipt
  }   // GetRange

  /////////////////////////////////////////
  bool
  FitSection(detinfo::DetectorClocksData const& clockData,
             detinfo::DetectorPropertiesData const& detProp,
             const TCSlice& slc,
             PFPStruct& pfp,
             unsigned short sfIndex)
  {
    // Fits the TP3D points in the selected section to a 3D line with the origin at the center of
    // the section
    if(pfp.TP3Ds.size() < 4) return false;
    if(sfIndex >= pfp.SectionFits.size()) return false;
    // don't fit a small angle PFP
    if(pfp.AlgMod[kSmallAng3D]) return true;

    unsigned short fromPt = USHRT_MAX;
    unsigned short npts = 0;
    GetRange(pfp, sfIndex, fromPt, npts);
    if (fromPt == USHRT_MAX) return false;
    if (npts < 4) return false;

    // check for errors
    for (unsigned short ipt = fromPt; ipt < fromPt + npts; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if (tp3d.SFIndex != sfIndex) return false;
    } // ipt

    // fit these points and update
    float chiDOF = 999;
    return FitPFP(clockData, detProp, slc, pfp, fromPt, npts, sfIndex, chiDOF);

  } // FitSection

  /////////////////////////////////////////
  SectionFit
  FitTP3Ds(detinfo::DetectorClocksData const& clockData,
           detinfo::DetectorPropertiesData const& detProp,
           const TCSlice& slc,
           const std::vector<TP3D>& tp3ds,
           unsigned short fromPt,
           short fitDir,
           unsigned short nPtsFit, bool prt)
  {
    // fits the points and returns the fit results in a SectionFit struct. This function assumes that the
    // vector of TP3Ds exists in the slc.TPCID

    SectionFit sf;
    sf.ChiDOF = FLT_MAX;
    // put the offset, cosine-like and sine-like components in a vector
    std::vector<std::array<double, 3>> ocs(slc.nPlanes);
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      auto planeID = geo::PlaneID(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      // plane offset
      ocs[plane][0] = tcc.geom->WireCoordinate(0, 0, planeID);
      // get the "cosine-like" component
      ocs[plane][1] = tcc.geom->WireCoordinate(1, 0, planeID) - ocs[plane][0];
      // the "sine-like" component
      ocs[plane][2] = tcc.geom->WireCoordinate(0, 1, planeID) - ocs[plane][0];
    } // plane

    // count the number of TPs in each plane
    std::vector<unsigned short> cntInPln(slc.nPlanes, 0);
    // and define the X position for the fit origin
    double x0 = 0.;
    // and define a vector of points to use
    std::vector<TP3D> useTP3Ds;
    // and their fit weights
    std::vector<double> weights;
    for(unsigned short ii = 0; ii < tp3ds.size(); ++ii) {
      short ipt = fromPt + fitDir * (short)ii;
      if (ipt < 0 || ipt >= (short)tp3ds.size()) break;
      auto& tp3d = tp3ds[ipt];
      if(!tp3d.Flags[kTP3DGood] || tp3d.TPIndex == USHRT_MAX) continue;
      if(tp3d.TPXErr2 < 0.0001) continue;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvUnusedHits]) continue;
      x0 += tp3d.TPX;
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      ++cntInPln[plane];
      useTP3Ds.push_back(tp3d);
      double weight = 1;
      // de-weight points that don't have the expected dE/dx
      if(tp3d.dEdx > 0) {
        if(tp3d.dEdx < 1.5) weight = 0.2;
        if(tp3d.Flags[kTP3DHiDEdx]) weight = 0.2;
      } // valid dEdx
      weight /= tp3d.TPXErr2;
      weights.push_back(weight);
      if(useTP3Ds.size() == nPtsFit) break;
    } // ipt
    if(useTP3Ds.size() < 4) return sf;
    nPtsFit = useTP3Ds.size();
    // ensure there are at least three points in at least two planes
    unsigned short enufInPlane = 0;
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane)
      if (cntInPln[plane] > 2) ++enufInPlane;
    if (enufInPlane < 2) return sf;
    // The X origin is the average X of all the points
    x0 /= (double)nPtsFit;

    TMatrixD A(nPtsFit, 4);
    TVectorD w(nPtsFit);
    for(unsigned short cnt = 0; cnt < useTP3Ds.size(); ++cnt) {
      auto& tp3d = useTP3Ds[cnt];
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double x = tp3d.TPX - x0;
      A[cnt][0] = weights[cnt] * ocs[plane][1];
      A[cnt][1] = weights[cnt] * ocs[plane][2];
      A[cnt][2] = weights[cnt] * ocs[plane][1] * x;
      A[cnt][3] = weights[cnt] * ocs[plane][2] * x;
      w[cnt] = weights[cnt] * (tp3d.Wire - ocs[plane][0]);
    } // ipt

    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    if(!ok) return sf;

    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    norm *= -1;
    sf.Dir[0] = 1 / norm;
    sf.Dir[1] = tVec[2] / norm;
    sf.Dir[2] = tVec[3] / norm;
    sf.Pos[0] = x0;
    sf.Pos[1] = tVec[0];
    sf.Pos[2] = tVec[1];
    sf.NPts = nPtsFit;

    // calculate ChiDOF
    sf.ChiDOF = 0;
    // project this 3D vector into a TP in every plane
    std::vector<TrajPoint> plnTP(slc.nPlanes);
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      plnTP[plane] = MakeBareTP(detProp, slc, sf.Pos, sf.Dir, inCTP);
    } // plane
    // a local position
    Point3_t pos;
    sf.DirErr[0] = 0.;
    for(unsigned short cnt = 0; cnt < useTP3Ds.size(); ++cnt) {
      auto& tp3d = useTP3Ds[cnt];
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double dw = tp3d.Wire - plnTP[plane].Pos[0];
      // dt/dW was stored in DeltaRMS by MakeBareTP
      double t = dw * plnTP[plane].DeltaRMS;
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
      // Note that the tp3d position is directly above the wire position and not the
      // point at the distance of closest approach. Delta is the difference in the
      // drift direction in cm
      double delta = pos[0] - tp3d.TPX;
//      sf.ChiDOF += fitWghts[cnt] * delta * delta;
      // estimate the X slope error ~ X direction vector with an overly simple average
      double dangErr = delta / t;
      sf.DirErr[0] += dangErr * dangErr;
    } // indx
    sf.ChiDOF = FitChiDOF(useTP3Ds, weights);
    sf.DirErr[0] = sqrt(sf.DirErr[0]) / (double)nPtsFit;
    sf.DirErr[1] = 3 * sf.DirErr[0];
    sf.DirErr[2] = 3 * sf.DirErr[0];
    sf.ChiDOF /= (float)(nPtsFit - 4);
    return sf;

  } // FitTP3Ds

  /////////////////////////////////////////
  double
  FitChiDOF(const std::vector<TP3D>& tp3ds, 
            const std::vector<double> weights)
  {
    if(tp3ds.size() < 5 ) return FLT_MAX;
    if(tp3ds.size() != weights.size()) return FLT_MAX;
    double chidof = 0.;
    double cnt = 0;
    for(unsigned short ipt = 0; ipt < tp3ds.size(); ++ipt) {
      if(weights[ipt] <= 0) break;
      auto& tp3d = tp3ds[ipt];
      double delta = tp3d.Pos[0] - tp3d.TPX;
      chidof += weights[ipt] * delta * delta;
      ++cnt;
    }
    if(cnt < 5) return FLT_MAX;
    return chidof / cnt;
  } // FitChiDOF
  /////////////////////////////////////////
  bool
  FitPFP(detinfo::DetectorClocksData const& clockData,
           detinfo::DetectorPropertiesData const& detProp,
           const TCSlice& slc,
           PFPStruct& pfp,
           unsigned short fromPt,
           unsigned short nPtsFit,
           unsigned short sfIndex,
           float& chiDOF)
  {
    // Fit points in the pfp.TP3Ds vector fromPt. This function
    // doesn't update the TP3Ds unless sfIndex refers to a valid SectionFit in the pfp.
    // No check is made to ensure that the TP3D SFIndex variable is compatible with sfIndex

    chiDOF = FLT_MAX;
    if (nPtsFit < 5) return false;
    if (fromPt + nPtsFit > pfp.TP3Ds.size()) return false;

    bool prt = (tcc.dbgPFP && pfp.MVI == debug.MVI);
    auto sf = FitTP3Ds(clockData, detProp, slc, pfp.TP3Ds, fromPt, 1, nPtsFit, prt);
    if (sf.ChiDOF == FLT_MAX) return false;

    // don't update the pfp?
    if (sfIndex >= pfp.SectionFits.size()) return true;

    // update the pfp Sectionfit
    pfp.SectionFits[sfIndex] = sf;
    // update the TP3Ds next
    // project this 3D vector into a TP in every plane
    std::vector<TrajPoint> plnTP(slc.nPlanes);
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      plnTP[plane] = MakeBareTP(detProp, slc, sf.Pos, sf.Dir, inCTP);
    } // plane

    Point3_t pos;
    bool needsSort = false;
    double prevAlong = 0;
    float nGood = 0;
    for(unsigned short ipt = fromPt; ipt < fromPt + nPtsFit; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.TPIndex == USHRT_MAX) continue;
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double dw = tp3d.Wire - plnTP[plane].Pos[0];
      // dt/dW was stored in DeltaRMS by MakeBareTP
      double t = dw * plnTP[plane].DeltaRMS;
      if (ipt == fromPt) { 
        prevAlong = t;
      }
      else {
        if (t < prevAlong) needsSort = true;
        prevAlong = t;
      }
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
      // Note that the tp3d position is directly above the wire position and not the
      // distance of closest approach. The Delta variable is the difference in the
      // drift direction in cm
      double delta = pos[0] - tp3d.TPX;
      tp3d.Pos = pos;
      tp3d.Dir = sf.Dir;
      tp3d.along = t;
      if(!tp3d.Flags[kTP3DGood]) continue;
      sf.ChiDOF += delta * delta / tp3d.TPXErr2;
      ++nGood;
    } // ipt
    if (needsSort) SortSection(pfp, sfIndex);
    pfp.SectionFits[sfIndex].NeedsUpdate = false;
    return true;

  } // FitPFP

  /////////////////////////////////////////
  void
  ReconcileVertices(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Checks for mis-placed 2D and 3D vertices and either attaches them
    // to a vertex or deletes(?) the vertex while attempting to preserve or
    // correct the P -> T -> 2V -> 3V assn. After this is done, the function
    // TCVertex/AttachToAnyVertex is called.
    // This function returns true if something was done to the pfp that requires
    // a re-definition of the pfp, e.g. adding or removing TP3Ds. Note that this
    // never occurs as the function is currently written

    if(tcc.vtx3DCuts.size() < 3) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.AlgMod[kJunk3D]) return;
//    if(pfp.AlgMod[kSmallAng3D]) return;

    // first make a list of all Tjs
    std::vector<int> tjList;
    for(auto& tp3d : pfp.TP3Ds) {
      if(!tp3d.Flags[kTP3DGood]) continue;
      // ignore single hits
      if (tp3d.TjID <= 0) continue;
      if (std::find(tjList.begin(), tjList.end(), tp3d.TjID) == tjList.end())
        tjList.push_back(tp3d.TjID);
    } // tp3d
    // look for 3D vertices associated with these Tjs and list of
    // orphan 2D vertices - those that are not matched to 3D vertices
    std::vector<int> vx2List, vx3List;
    for (auto tid : tjList) {
      auto& tj = slc.tjs[tid - 1];
      for (unsigned short end = 0; end < 2; ++end) {
        if (tj.VtxID[end] <= 0) continue;
        auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
        if (vx2.Vx3ID > 0) {
          if (std::find(vx3List.begin(), vx3List.end(), vx2.Vx3ID) == vx3List.end())
            vx3List.push_back(vx2.Vx3ID);
          // 3D vertex exists
        }
        else {
          // no 3D vertex
          if (std::find(vx2List.begin(), vx2List.end(), tj.VtxID[end]) == vx2List.end())
            vx2List.push_back(tj.VtxID[end]);
        } // no 3D vertex
      }   // end
    }     // tid
    // no vertex reconciliation is necessary
    if (vx2List.empty() && vx3List.empty()) return;
    if (prt) {
      mf::LogVerbatim myprt("TC");
      myprt << "RV: P" << pfp.ID << " ->";
      for (auto tid : tjList)
        myprt << " T" << tid;
      myprt << " ->";
      for (auto vid : vx3List)
        myprt << " 3V" << vid;
      if (!vx2List.empty()) {
        myprt << " orphan";
        for (auto vid : vx2List)
          myprt << " 2V" << vid;
      }
    } // prt
    // Just kill the orphan 2D vertices regardless of their score.
    // This is an indicator that the vertex was created between two tjs
    // that maybe should have been reconstructed as one or alternatively
    // as two Tjs. This decision presumes the existence of a 3D kink
    // algorithm that doesn't yet exist...
    for (auto vid : vx2List) {
      auto& vx2 = slc.vtxs[vid - 1];
      MakeVertexObsolete("RV", slc, vx2, true);
    } // vx2List
    // ignore the T -> 2V -> 3V assns (if any exist) and try to directly
    // attach to 3D vertices at both ends
    AttachToAnyVertex(slc, pfp, tcc.vtx3DCuts[2], prt);
    // check for differences and while we are here, see if the pfp was attached
    // to a neutrino vertex and the direction is wrong
    int neutrinoVx = 0;
    if (!slc.pfps.empty()) {
      auto& npfp = slc.pfps[0];
      bool neutrinoPFP = (npfp.PDGCode == 12 || npfp.PDGCode == 14);
      if (neutrinoPFP) neutrinoVx = npfp.Vx3ID[0];
    } // pfps exist
    unsigned short neutrinoVxEnd = 2;
    for (unsigned short end = 0; end < 2; ++end) {
      // see if a vertex got attached
      if (pfp.Vx3ID[end] <= 0) continue;
      if (pfp.Vx3ID[end] == neutrinoVx) neutrinoVxEnd = end;
      // see if this is a vertex in the list using the T -> 2V -> 3V assns
      if (std::find(vx3List.begin(), vx3List.end(), pfp.Vx3ID[end]) != vx3List.end()) continue;
    } // end
    if (neutrinoVxEnd < 2 && neutrinoVxEnd != 0) Reverse(slc, pfp);

    return;
  } // ReconcileVertices

  /////////////////////////////////////////
  void FillGaps3D(detinfo::DetectorClocksData const& clockData,
                  detinfo::DetectorPropertiesData const& detProp,
                  TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Look for gaps in each plane in the TP3Ds vector. Hits
    // reconstructed at large angles are poorly reconstructed which results
    // in poorly reconstructed 2D trajectories

    if(pfp.ID <= 0) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.SectionFits.empty()) return;
    if(!tcc.useAlg[kFillGaps3D]) return;
    if(pfp.AlgMod[kJunk3D]) return;
    if(tcc.useAlg[kNewCuts]) return;

    // make a copy in case something goes wrong
    auto pSave = pfp;

    if(!pfp.Flags[kdEdxDefined]) SetTP3DdEdx(clockData, detProp, slc, pfp);

    // Find the average dE/dx so we can apply a generous min/max dE/dx cut
    float dEdXAve = 0;
    float dEdXRms = 0;
    MPV_dEdX(clockData, detProp, slc, pfp, dEdXAve, dEdXRms);
    float dEdxMin = 0.5, dEdxMax = 50.;
    if (dEdXAve > 0.5) {
      dEdxMin = dEdXAve - 5 * dEdXRms;
      if (dEdxMin < 0.5) dEdxMin = 0.5;
      dEdxMax = dEdXAve + 5 * dEdXRms;
      if (dEdxMax > 50.) dEdxMax = 50.;
    } // dEdXAve > 0.5

    // iterate over all SectionFits
    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      // iterate over all planes
      unsigned short nadd = 0;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
        // find the start and end TP3D in this plane in this section
        unsigned short startPt = USHRT_MAX;
        unsigned short endPt = 0;
        // a list of (TjID, ipt) pairs that are used in this section and plane
        std::vector<std::pair<int, unsigned short>> tpUsed;
        for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
          auto& tp3d = pfp.TP3Ds[ipt];
          if(tp3d.SFIndex < sfi) continue;
          if(tp3d.SFIndex > sfi) break;
          if(startPt == USHRT_MAX) startPt = ipt;
          endPt = ipt;
          // add the (Tj, ipt) pair to the list 
          if(tp3d.CTP == inCTP) tpUsed.push_back(std::make_pair(tp3d.TjID, tp3d.TPIndex));
        } //ipt
        if(startPt == USHRT_MAX) continue;
        // Project the 3D position at both ends into the plane 
        auto startTP = MakeBareTP(detProp, slc, pfp.TP3Ds[startPt].Pos, inCTP);
        auto endTP = MakeBareTP(detProp, slc, pfp.TP3Ds[endPt].Pos, inCTP);
        if(startTP.Pos[0] < -0.5 || endTP.Pos[0] < -0.5) {
          if(prt) std::cout<<"FillGaps3D: Invalid TPs in P"<<pfp.ID
            <<" "<<startTP.Pos[0]<<" "<<endTP.Pos[0]<<"\n";
          continue;
        }
        float startWire = startTP.Pos[0];
        if(startWire < -0.5 || startWire > slc.nWires[plane]) {
          if(prt) std::cout<<"FillGaps3D: invalid start wire in P"<<pfp.ID<<" "<<startWire<<"\n";
          continue;
        }
        float endWire = endTP.Pos[0];
        if(endWire < -0.5 || endWire > slc.nWires[plane]) {
          if(prt) std::cout<<"FillGaps3D: invalid start wire in P"<<pfp.ID<<" "<<endWire<<"\n";
          continue;
        }
        float nChkWires = std::abs(endWire - startWire);
        double along = pfp.TP3Ds[startPt].along;
        double dtdw = PosSep(pfp.TP3Ds[startPt].Pos, pfp.TP3Ds[endPt].Pos) / (nChkWires - 1);
        if(startWire > endWire) {
          std::swap(startWire, endWire);
          std::swap(startTP, endTP);
          along = pfp.TP3Ds[endPt].along;
          dtdw *= -1;
        }
        if(nChkWires < 3) continue;
        if((float)tpUsed.size() > 0.7 * nChkWires) continue;
        unsigned int iStartWire = std::nearbyint(startWire);
        unsigned int iEndWire = std::nearbyint(endWire);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<"FG3D P"<<pfp.ID<<" sfi "<<sfi<<" plane "<<plane;
          myprt<<" nChkWires "<<nChkWires;
          myprt<<" wires "<<iStartWire<<" - "<<iEndWire;
          myprt<<" dT/dw "<<dtdw;
          myprt<<" along start "<<along;
          myprt<<" projInPln "<<startTP.Delta;
        } // prt
        // look for missing hits in the wire range
        for(unsigned int wire = iStartWire; wire <= iEndWire; ++wire) {
          // move startTP to this wire position
          MoveTPToWire(startTP, (float)wire);
          if(wire > iStartWire) along += dtdw;
          // set a generous search window in WSE units
          float window = 30;
          FindCloseHits(slc, startTP, window, kUsedHits);
          if(startTP.Hits.empty()) continue;
          float bestPull = tcc.match3DCuts[4];
          TP3D bestTP3D;
          for (auto iht : startTP.Hits) {
            if (slc.slHits[iht].InTraj <= 0) continue;
            // this hit is used in a TP so find the tpIndex
            auto& utj = slc.tjs[slc.slHits[iht].InTraj - 1];
            unsigned short tpIndex = 0;
            for (tpIndex = utj.EndPt[0]; tpIndex <= utj.EndPt[1]; ++tpIndex) {
              auto& utp = utj.Pts[tpIndex];
              if (utp.Chg <= 0) continue;
              // This doesn't check for UseHit true but that is probably ok here
              if (std::find(utp.Hits.begin(), utp.Hits.end(), iht) != utp.Hits.end()) break;
            } // ipt
            if (tpIndex > utj.EndPt[1]) continue;
            // see if it is already used in this pfp
            std::pair<int, unsigned short> tppr = std::make_pair(utj.ID, tpIndex);
            if (std::find(tpUsed.begin(), tpUsed.end(), tppr) != tpUsed.end()) continue;
            tpUsed.push_back(tppr);
            auto& utp = utj.Pts[tpIndex];
            // see if it is used in a different PFP
            if (utp.InPFP > 0) continue;
            // or if it overlaps another trajectory near a 2D vertex
            if (utp.Environment[kEnvOverlap]) continue;
            auto newTP3D = CreateTP3D(detProp, slc, utj.ID, tpIndex);
            if(newTP3D.Flags[kTP3DBad]) continue;
            // inflate the X error if the projection in this plane is small
            if(startTP.Delta < 0.3) newTP3D.TPXErr2 *= 4;
            newTP3D.SFIndex = sfi;
            newTP3D.along = along;
            // set the direction to the direction of the SectionFit it is in so we can calculate dE/dx
            auto& sf = pfp.SectionFits[sfi];
            newTP3D.Dir = sf.Dir;
            for(unsigned short xyz = 0; xyz < 3; ++xyz) newTP3D.Pos[xyz] = sf.Pos[xyz] + along * sf.Dir[xyz];
            float pull = PointPull(pfp, newTP3D);
            if(pull > bestPull) continue;
            newTP3D.dEdx = dEdx(clockData, detProp, slc, newTP3D);
            if(newTP3D.dEdx < dEdxMin || newTP3D.dEdx > dEdxMax) continue;
            newTP3D.Flags[kTP3DGood] = true;
            bestTP3D = newTP3D;
            bestPull = pull;
          } // iht
          if(bestTP3D.TjID <= 0) continue;
          if(bestPull > tcc.match3DCuts[4]) continue;
/*
          if (prt && bestPull < 10) {
            mf::LogVerbatim myprt("TC");
            auto& tp = slc.tjs[bestTP3D.TjID - 1].Pts[bestTP3D.TPIndex];
            myprt << "FG3D: P" << pfp.ID << " added T"<<bestTP3D.TjID;
            myprt << " TP " << PrintPos(slc, tp);
            myprt << " pull " << std::fixed << std::setprecision(2) << bestPull;
            myprt << " dx " << bestTP3D.TPX - bestTP3D.Pos[0] << " in section " << bestTP3D.SFIndex;
          }
*/
          if (InsertTP3D(pfp, bestTP3D) == USHRT_MAX) continue;
          ++nadd;
          pfp.SectionFits[sfi].NeedsUpdate = true;
          pfp.Flags[kNeedsUpdate] = true;
        } // wire
      } // plane
      if(prt) mf::LogVerbatim("TC")<<" added "<<nadd<<" points in section "<<sfi;
    } // sf

    if(!pfp.Flags[kNeedsUpdate]) return;
    pfp.AlgMod[kFillGaps3D] = true;
    if(!Update(clockData, detProp, slc, pfp, prt))  {
      if(prt) mf::LogVerbatim("TC")<<" Update failed. Recovering...";
      pfp = pSave;
    }

  } // FillGaps3D

  /////////////////////////////////////////
  bool
  ValidTwoPlaneMatch(detinfo::DetectorPropertiesData const& detProp,
                     const TCSlice& slc,
                     const PFPStruct& pfp)
  {
    // This function checks the third plane in the PFP when only two Tjs are 3D-matched to
    // ensure that the reason for the lack of a 3rd plane match is that it is in a dead region.
    // This function should be used after an initial fit is done and the TP3Ds are sorted
    if (pfp.TjIDs.size() != 2) return false;
    if (slc.nPlanes != 3) return false;
    if (pfp.TP3Ds.empty()) return false;

    // find the third plane
    std::vector<unsigned short> planes;
    for (auto tid : pfp.TjIDs)
      planes.push_back(DecodeCTP(slc.tjs[tid - 1].CTP).Plane);
    unsigned short thirdPlane = 3 - planes[0] - planes[1];
    CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, thirdPlane);
    // Project the 3D position at the start into the third plane
    auto tp = MakeBareTP(detProp, slc, pfp.TP3Ds[0].Pos, inCTP);
    unsigned int wire0 = 0;
    if (tp.Pos[0] > 0) wire0 = std::nearbyint(tp.Pos[0]);
    if (wire0 > slc.nWires[thirdPlane]) wire0 = slc.nWires[thirdPlane];
    // Do the same for the end
    unsigned short lastPt = pfp.TP3Ds.size() - 1;
    tp = MakeBareTP(detProp, slc, pfp.TP3Ds[lastPt].Pos, inCTP);
    unsigned int wire1 = 0;
    if (tp.Pos[0] > 0) wire1 = std::nearbyint(tp.Pos[0]);
    if (wire1 > slc.nWires[thirdPlane]) wire1 = slc.nWires[thirdPlane];
    if (wire0 == wire1) return !evt.goodWire[thirdPlane][wire0];
    if (wire1 < wire0) std::swap(wire0, wire1);
    // count the number of good wires
    int dead = 0;
    int wires = wire1 - wire0;
    for (unsigned int wire = wire0; wire < wire1; ++wire)
      if (!evt.goodWire[thirdPlane][wire]) ++dead;
    // require that most of the wires are dead
    return (dead > 0.8 * wires);
  } // ValidTwoPlaneMatch

  /////////////////////////////////////////
  unsigned short
  InsertTP3D(PFPStruct& pfp, TP3D& tp3d)
  {
    // inserts the tp3d into the section defined by tp3d.SFIndex
    if (tp3d.SFIndex >= pfp.SectionFits.size()) return USHRT_MAX;
    // Find the first occurrence of this SFIndex
    size_t ipt = 0;
    for (ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt)
      if (tp3d.SFIndex == pfp.TP3Ds[ipt].SFIndex) break;
    if (ipt == pfp.TP3Ds.size()) {
      // Must be adding the first TP3D in a SectionFit
      if(tp3d.SFIndex != pfp.SectionFits.size()-1) return USHRT_MAX;
      pfp.TP3Ds.push_back(tp3d);
      pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
      pfp.Flags[kNeedsUpdate] = true;
      return pfp.TP3Ds.size() - 1;      
    }
    // next see if we can insert it so that re-sorting of this section isn't required
    auto lastTP3D = pfp.TP3Ds.back();
    if (ipt == 0 && tp3d.along <= pfp.TP3Ds[0].along) {
      // insert at the beginning. No search needs to be done
    }
    else if (tp3d.SFIndex == lastTP3D.SFIndex && tp3d.along > lastTP3D.along) {
      // insert at the end. Use push_back and return
      pfp.TP3Ds.push_back(tp3d);
      pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
      pfp.Flags[kNeedsUpdate] = true;
      return pfp.TP3Ds.size() - 1;
    }
    else {
      for (std::size_t iipt = ipt; iipt < pfp.TP3Ds.size() - 1; ++iipt) {
        // break out if the next point is in a different section
        if (pfp.TP3Ds[iipt + 1].SFIndex != tp3d.SFIndex) break;
        if (tp3d.along > pfp.TP3Ds[iipt].along && tp3d.along < pfp.TP3Ds[iipt + 1].along) {
          ipt = iipt + 1;
          break;
        }
      } // iipt
    }   // insert in the middle
    pfp.TP3Ds.insert(pfp.TP3Ds.begin() + ipt, tp3d);
    pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
    pfp.Flags[kNeedsUpdate] = true;
    return ipt;
  } // InsertTP3D

  /////////////////////////////////////////
  bool
  SortSection(PFPStruct& pfp, unsigned short sfIndex)
  {
    // sorts the TP3Ds by the distance from the start of a fit section

    if(sfIndex > pfp.SectionFits.size() - 1) return false;
    auto& sf = pfp.SectionFits[sfIndex];
    if (sf.Pos[0] == 0.0 && sf.Pos[1] == 0.0 && sf.Pos[2] == 0.0) return false;

    // a temp vector of points in this section
    std::vector<TP3D> temp;
    // and the index into TP3Ds
    std::vector<unsigned short> indx;
    // See if the along variable is monotonically increasing
    float prevAlong = 0;
    bool first = true;
    bool needsSort = false;
    for (std::size_t ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      auto& tp3d = pfp.TP3Ds[ii];
      if(tp3d.TPIndex == USHRT_MAX) {
        // don't move the end TPs
        temp.push_back(tp3d);
        indx.push_back(ii);
        continue;
      }
      if(tp3d.SFIndex != sfIndex) continue;
      if(first) {
        first = false;
        prevAlong = tp3d.along;
      }
      else {
        if (tp3d.along < prevAlong) needsSort = true;
        prevAlong = tp3d.along;
      }
      temp.push_back(tp3d);
      indx.push_back(ii);
    } // tp3d
    if (temp.empty()) return false;
    // no sort needed?
    if (temp.size() == 1) return true;
    if (!needsSort) {
      sf.NeedsUpdate = false;
      return true;
    }
    // see if the points are not-contiguous
    bool contiguous = true;
    for (std::size_t ipt = 1; ipt < indx.size(); ++ipt) {
      if (indx[ipt] != indx[ipt - 1] + 1) contiguous = false;
    } // ipt
    if (!contiguous) { return false; }

    std::vector<SortEntry> sortVec(temp.size());
    for (std::size_t ii = 0; ii < temp.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].val = temp[ii].along;
    } // ipt
    std::sort(sortVec.begin(), sortVec.end(), valsIncreasing);
    for (std::size_t ii = 0; ii < temp.size(); ++ii) {
      // overwrite the tp3d
      auto& tp3d = pfp.TP3Ds[indx[ii]];
      tp3d = temp[sortVec[ii].index];
    } // ii
    sf.NeedsUpdate = false;
    return true;
  } // SortSection

  /////////////////////////////////////////
  void
  SortByX(std::vector<TP3D>& tp3ds)
  {
    // Sorts TP3Ds by increasing X hit position
    if(tp3ds.empty()) return;
    std::vector<SortEntry> sortVec(tp3ds.size());
    for(std::size_t ii = 0; ii < tp3ds.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].val = tp3ds[ii].TPX;
    } // ii
    std::sort(sortVec.begin(), sortVec.end(), valsIncreasing);
    auto temp = tp3ds;
    for(std::size_t ii = 0; ii < tp3ds.size(); ++ii) tp3ds[ii] = temp[sortVec[ii].index];
  } // SortByX

  /////////////////////////////////////////
  void
  Recover(detinfo::DetectorClocksData const& clockData,
          detinfo::DetectorPropertiesData const& detProp,
          TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // try to recover from a poor initial fit
    if(pfp.AlgMod[kSmallAng3D]) return;
    if(pfp.SectionFits.size() != 1) return;
    if(pfp.TP3Ds.size() < 20) return;
    if(!CanSection(slc, pfp)) return;

    // make a copy
    auto p2 = pfp;
    // try two sections
    p2.SectionFits.resize(2);
    unsigned short halfPt = p2.TP3Ds.size() / 2;
    for(unsigned short ipt = halfPt; ipt < p2.TP3Ds.size(); ++ipt) p2.TP3Ds[ipt].SFIndex = 1;
    // Confirm that both sections can be reconstructed
    unsigned short toPt = Find3DRecoRange(slc, p2, 0, 3, 1);
    if(toPt > p2.TP3Ds.size()) return;
    toPt = Find3DRecoRange(slc, p2, halfPt, 3, 1);
    if(toPt > p2.TP3Ds.size()) return;
    if(!FitSection(clockData, detProp, slc, p2, 0) || !FitSection(clockData, detProp, slc, p2, 1)) {
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt << "Recover failed MVI " << p2.MVI << " in TPC " << p2.TPCID.TPC;
        for(auto tid : p2.TjIDs) myprt << " T" << tid; 
      } // prt
      return;
    }
    if(prt) mf::LogVerbatim("TC")<<"Recover: P" << pfp.ID << " success";
    pfp = p2;

  } // Recover

  /////////////////////////////////////////
  bool
  MakeTP3Ds(detinfo::DetectorPropertiesData const& detProp, TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Create and populate the TP3Ds vector. This function is called before the first
    // fit is done so the TP3D along variable can't be determined. It returns false
    // if a majority of the tj points in TjIDs are already assigned to a different pfp
    pfp.TP3Ds.clear();
    if (!pfp.TP3Ds.empty() || pfp.SectionFits.size() != 1) return false;

    // Look for TPs that are ~parallel to the wire plane and trajectory points
    // where the min/max Pos[1] value is not near an end.
    // number of Small Angle Tjs
    // stash the inflection point index in the TjUIDs vector
    pfp.TjUIDs.resize(pfp.TjIDs.size(), -1);
    // count TPs with charge in all of the Tjs
    float cnt = 0;
    // and the number of TPs available for use
    float avail = 0;
    // and the number of junk Tjs
    unsigned short nJunk = 0;
    unsigned short nSA = 0;
    for(unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
      auto& tj = slc.tjs[pfp.TjIDs[itj] - 1];
      if(tj.AlgMod[kJunkTj]) ++nJunk;
      float posMin = 1E6;
      unsigned short iptMin = USHRT_MAX;
      float posMax = -1E6;
      unsigned short iptMax = USHRT_MAX;
      float aveAng = 0;
      float npwc = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        ++cnt;
        if (tp.InPFP > 0) continue;
        ++avail;
        if(tp.Pos[1] > posMax) { posMax = tp.Pos[1]; iptMax = ipt; }
        if(tp.Pos[1] < posMin) { posMin = tp.Pos[1]; iptMin = ipt; }
        aveAng += tp.Ang;
        ++npwc;
      } // ipt
      if(npwc == 0) continue;
      aveAng /= npwc;
      if(std::abs(aveAng) < 0.05) ++nSA;
      // No problem if the min/max points are near the ends
      if(iptMin > tj.EndPt[0] + 4 && iptMin < tj.EndPt[1] - 4) pfp.TjUIDs[itj] = iptMin;
      if(iptMax > tj.EndPt[0] + 4 && iptMax < tj.EndPt[1] - 4) pfp.TjUIDs[itj] = iptMax;
    } // tid
    if(avail < 0.8 * cnt) return false;
    // small angle trajectory?
    if(nSA > 1) pfp.AlgMod[kSmallAng3D] = true;
    if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" MVI "<<pfp.MVI<<" nJunkTj "<<nJunk<<" SmallAngle? "<<pfp.AlgMod[kSmallAng3D];

    if(pfp.AlgMod[kSmallAng3D]) return MakeSmallAnglePFP(detProp, slc, pfp, prt);

    // Add the points associated with the Tjs that were used to create the PFP
    for (auto tid : pfp.TjIDs) {
      auto& tj = slc.tjs[tid - 1];
      // There is one TP for every hit in a junk Tj so we can skip one, if there is only one
      if(nJunk == 1 && tj.AlgMod[kJunkTj]) continue;
      // All of the Tj's may be junk, especially for those at very high angle, so the
      // X position of the TP's isn't high quality. Inflate the errors below.
      bool isJunk = tj.AlgMod[kJunkTj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if (tp.Chg <= 0) continue;
        if (tp.InPFP > 0) continue;
        ++avail;
        auto tp3d = CreateTP3D(detProp, slc, tid, ipt);
        if(tp3d.Flags[kTP3DBad]) continue;
        tp3d.SFIndex = 0;
        if(isJunk) tp3d.TPXErr2 *= 4;
        // We need to assume that all points are good or the first fit will fail
        tp3d.Flags[kTP3DGood] = true;
        pfp.TP3Ds.push_back(tp3d);
      } // ipt
    } // tid
    if(prt) mf::LogVerbatim("TC")<<" has "<<pfp.TP3Ds.size()<<" TP3Ds";
    return true;
  } // MakeTP3Ds

  /////////////////////////////////////////
  bool MakeSmallAnglePFP(detinfo::DetectorPropertiesData const& detProp,
                         TCSlice& slc, PFPStruct& pfp,
                         bool prt)
  {
    // Create and populate the TP3Ds vector for a small-angle track. The standard track fit
    // will fail for these tracks. The kSmallAng3D AlgMod bit
    // is set true. Assume that the calling function, MakeTP3Ds, has decided that this is a
    // small-angle track. 

    if(!tcc.useAlg[kSmallAng3D]) return false;
    if(pfp.TjIDs.size() < 2) return false;

    std::vector<SortEntry> sortVec(pfp.TjIDs.size());
    unsigned short sbCnt = 0;
    for (unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
      sortVec[itj].index = itj;
      auto& tj = slc.tjs[pfp.TjIDs[itj] - 1];
      sortVec[itj].val = NumPtsWithCharge(slc, tj, false);
      if(pfp.TjUIDs[itj] > 0) ++sbCnt;
    } // ipt
    std::sort(sortVec.begin(), sortVec.end(), valsDecreasing);

    // Decide whether to use the inflection points to add another section. Inflection
    // points must exist in the two longest Tjs
    unsigned short tlIndex = sortVec[0].index;
    unsigned short nlIndex = sortVec[1].index;
    auto& tlong = slc.tjs[pfp.TjIDs[tlIndex] - 1];
    auto& nlong = slc.tjs[pfp.TjIDs[nlIndex] - 1];
    bool twoSections = (sbCnt > 1 && pfp.TjUIDs[tlIndex] > 0 && pfp.TjUIDs[nlIndex] > 0);
    unsigned short tStartPt = tlong.EndPt[0];
    unsigned short tEndPt = tlong.EndPt[1];
    unsigned short nStartPt = nlong.EndPt[0];
    unsigned short nEndPt = nlong.EndPt[1];
    if(twoSections) {
      pfp.SectionFits.resize(2);
      tEndPt = pfp.TjUIDs[tlIndex];
      nEndPt = pfp.TjUIDs[nlIndex];
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"MakeSmallAnglePFP: creating two sections using points";
        myprt<<" T"<<tlong.ID<<"_"<<tEndPt;
        myprt<<" T"<<nlong.ID<<"_"<<nEndPt;
      } // prt
    } // two Sections
    std::vector<Point3_t> sfEndPos;
    for(unsigned short isf = 0; isf < pfp.SectionFits.size(); ++isf) {
      // get the start and end TPs in this section
      auto& ltp0 = tlong.Pts[tStartPt];
      auto& ltp1 = tlong.Pts[tEndPt];
      auto& ntp0 = nlong.Pts[nStartPt];
      auto& ntp1 = nlong.Pts[nEndPt];
      // Get the 3D end points
      auto start = MakeTP3D(detProp, slc, ltp0, ntp0);
      auto end = MakeTP3D(detProp, slc, ltp1, ntp1);
      if(!start.Flags[kTP3DGood] || !end.Flags[kTP3DGood]) {
//        std::cout<<" Start/end fail in section "<<isf<<". Add recovery code\n";
        return false;
      } // failure
      if(isf == 0) sfEndPos.push_back(start.Pos);
      sfEndPos.push_back(end.Pos);
      auto& sf = pfp.SectionFits[isf];
      // Find the start and end positions
      for(unsigned short xyz = 0; xyz < 3; ++xyz) {
        sf.Dir[xyz] = end.Pos[xyz] - start.Pos[xyz];
        sf.Pos[xyz] = (end.Pos[xyz] + start.Pos[xyz]) / 2.;
      }
      SetMag(sf.Dir, 1.);
      sf.ChiDOF = 0.;
      sf.NPts = 0;
      // move the start/end point indices
      tStartPt = tEndPt + 1; tEndPt = tlong.EndPt[1];
      nStartPt = nEndPt + 1; nEndPt = nlong.EndPt[1];
    } // isf
    // Create TP3Ds
    // a temporary vector to hold TP3Ds for the second SectionFit
    std::vector<TP3D> sf2pts;
    for(unsigned short itj = 0; itj < sortVec.size(); ++itj) {
      int tid = pfp.TjIDs[sortVec[itj].index];
      // don't add points for the Tj that doesn't have an inflection point. It is
      // probably broken and would probably be put in the wrong section
      if(twoSections && pfp.TjUIDs[sortVec[itj].index] < 0) continue;
      auto& tj = slc.tjs[tid - 1];
      unsigned short sb = tj.EndPt[1];
      if(twoSections && pfp.TjUIDs[sortVec[itj].index] > 0) sb = pfp.TjUIDs[sortVec[itj].index];
      // count the number of good TPs in each section
      std::vector<double> npwc(pfp.SectionFits.size(), 0);
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(ipt > sb) { ++npwc[1]; } else { ++npwc[0]; }
      } // ipt
      double length = PosSep(sfEndPos[0], sfEndPos[1]);
      double step = length / npwc[0];
      double along = -length / 2;
      unsigned short sfi = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        auto tp3d = CreateTP3D(detProp, slc, tid, ipt);
        if(tp3d.Flags[kTP3DBad]) continue;
        if(ipt == sb + 1) {
          sfi = 1;
          length = PosSep(sfEndPos[1], sfEndPos[2]);
          step = length / npwc[1];
          along = -length / 2;
        }
        tp3d.SFIndex = sfi;
        auto& sf = pfp.SectionFits[sfi];
        ++sf.NPts;
        tp3d.along = along;
        for(unsigned short xyz = 0; xyz < 3; ++xyz) tp3d.Pos[xyz] = sf.Pos[xyz] + along * sf.Dir[xyz];
        tp3d.Dir = sf.Dir;
        along += step;
        double delta = tp3d.Pos[0] - tp3d.TPX;
        sf.ChiDOF += delta * delta / tp3d.TPXErr2;
        // Assume that all points are good
        tp3d.Flags[kTP3DGood] = true;
        if(sfi == 0) {
          pfp.TP3Ds.push_back(tp3d);
        } else {
          sf2pts.push_back(tp3d);
        }
      } // ipt
    } // tid
    if(pfp.TP3Ds.size() < 4) return false;
    for(auto& sf : pfp.SectionFits) {
      if(sf.NPts < 5) return false;
      sf.ChiDOF /= (float)(sf.NPts - 4);
    } // sf
    if(!SortSection(pfp, 0)) return false;
    if(!sf2pts.empty()) {
      // append the points and sort
      pfp.TP3Ds.insert(pfp.TP3Ds.end(), sf2pts.begin(), sf2pts.end());
      if(!SortSection(pfp, 1)) return false;
    } // two sections
    pfp.Flags[kCanSection] = false;
    pfp.AlgMod[kSmallAng3D] = true;
    pfp.Flags[kdEdxDefined] = false;
    if(prt) {
      mf::LogVerbatim("TC")<<"Created SmallAngle P"<<pfp.ID
          <<" with "<<pfp.TP3Ds.size()
          <<" points in "<<pfp.SectionFits.size()<<" sections\n";
    }
    return true;
  } // MakeSmallAnglePFP


  /////////////////////////////////////////
  void
  Reverse(TCSlice& slc, PFPStruct& pfp)
  {
    // reverse the PFParticle
    std::reverse(pfp.TP3Ds.begin(), pfp.TP3Ds.end());
    if(pfp.SectionFits.size() > 1) std::reverse(pfp.SectionFits.begin(), pfp.SectionFits.end());
    for (std::size_t sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      // flip the direction vector
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        sf.Dir[xyz] *= -1;
    } // sf
    // correct the along variable
    for(auto& tp3d : pfp.TP3Ds) {
      tp3d.along *= -1;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) tp3d.Dir[xyz] *= -1;
      tp3d.SFIndex = pfp.SectionFits.size() - 1 - tp3d.SFIndex;
    } 
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    std::swap(pfp.EndFlag[0], pfp.EndFlag[1]);
  } // Reverse

  /////////////////////////////////////////
  void
  FillmAllTraj(detinfo::DetectorPropertiesData const& detProp, TCSlice& slc)
  {
    // Fills the mallTraj vector with trajectory points in the tpc and sorts
    // them by increasing X

    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;

    // define mallTraj
    slc.mallTraj.clear();
    Tj2Pt tj2pt;
    // Count the number of shower-like trajectories
    float nShLike = 0;
    for (auto& tj : slc.tjs) {
      if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      // ignore already matched
      if (tj.AlgMod[kMat3D]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if ((int)planeID.Cryostat != cstat) continue;
      if ((int)planeID.TPC != tpc) continue;
      if (tj.ID <= 0) continue;
      if (tj.PDGCode == 11) nShLike += tj.Pts.size();
    } // tj
    bool ignoreShLike = (nShLike > tcc.match3DCuts[6]);

    float rms = tcc.match3DCuts[0];
    for (auto& tj : slc.tjs) {
      if (tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      // ignore already matched
      if (tj.AlgMod[kMat3D]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if ((int)planeID.Cryostat != cstat) continue;
      if ((int)planeID.TPC != tpc) continue;
      int plane = planeID.Plane;
      if (tj.ID <= 0) continue;
      if(ignoreShLike && tj.PDGCode == 11) continue;
      unsigned short tjID = tj.ID;
      for (unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if (tp.Chg <= 0) continue;
        if (tp.Pos[0] < -0.4) continue;
        // ignore already matched
        if (tp.InPFP > 0) continue;
        tj2pt.wire = std::nearbyint(tp.Pos[0]);
        // don't try matching if the wire doesn't exist
        if (!tcc.geom->HasWire(geo::WireID(cstat, tpc, plane, tj2pt.wire))) continue;
        float xpos = detProp.ConvertTicksToX(tp.Pos[1] / tcc.unitsPerTick, plane, tpc, cstat);
        tj2pt.xlo = xpos - rms;
        tj2pt.xhi = xpos + rms;
        tj2pt.plane = plane;
        tj2pt.id = tjID;
        tj2pt.ipt = ipt;
        tj2pt.npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        slc.mallTraj.push_back(tj2pt);
      } // tp
    }   // tj

    // sort by increasing x
    std::vector<SortEntry> sortVec(slc.mallTraj.size());
    for (std::size_t ipt = 0; ipt < slc.mallTraj.size(); ++ipt) {
      // populate the sort vector
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = slc.mallTraj[ipt].xlo;
    } // ipt
    // sort by increasing xlo
    std::sort(sortVec.begin(), sortVec.end(), valsIncreasing);
    // put slc.mallTraj into sorted order
    auto tallTraj = slc.mallTraj;
    for (std::size_t ii = 0; ii < sortVec.size(); ++ii)
      slc.mallTraj[ii] = tallTraj[sortVec[ii].index];

  } // FillmAllTraj

  /////////////////////////////////////////
  TP3D MakeTP3D(detinfo::DetectorPropertiesData const& detProp, 
                TCSlice& slc, const TrajPoint& itp, const TrajPoint& jtp)
  {
    // Make a 3D trajectory point using two 2D trajectory points. The TP3D Pos and Wire
    // variables are defined using itp. The SectionFit variables are un-defined
    TP3D tp3d;
    tp3d.TPIndex = 0;
    tp3d.TjID = 0;
    tp3d.CTP = itp.CTP;
    // assume failure
    tp3d.Flags[kTP3DGood] = false;
    tp3d.Dir = {{0.0, 0.0, 1.0}};
    tp3d.Pos = {{999.0, 999.0, 999.0}};
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);
    if(iPlnID == jPlnID) return tp3d;
    double upt = tcc.unitsPerTick;
    double ix = detProp.ConvertTicksToX(itp.Pos[1] / upt, iPlnID);
    double jx = detProp.ConvertTicksToX(jtp.Pos[1] / upt, jPlnID);
    
    // don't continue if the points are wildly far apart in X
    double dx = std::abs(ix - jx);
    if(dx > 20) return tp3d;
    tp3d.Pos[0] = (ix + jx) / 2;
    tp3d.TPX = ix;
    // Fake the error
    tp3d.TPXErr2 = dx;
    // determine the wire orientation and offsets using WireCoordinate
    // wire = yp * OrthY + zp * OrthZ - Wire0 = cs * yp + sn * zp - wire0
    // wire offset
    double iw0 = tcc.geom->WireCoordinate(0, 0, iPlnID);
    // cosine-like component
    double ics = tcc.geom->WireCoordinate(1, 0, iPlnID) - iw0;
    // sine-like component
    double isn = tcc.geom->WireCoordinate(0, 1, iPlnID) - iw0;
    double jw0 = tcc.geom->WireCoordinate(0, 0, jPlnID);
    double jcs = tcc.geom->WireCoordinate(1, 0, jPlnID) - jw0;
    double jsn = tcc.geom->WireCoordinate(0, 1, jPlnID) - jw0;
    double den = isn * jcs - ics * jsn;
    if(den == 0) return tp3d;
    double iPos0 = itp.Pos[0];
    double jPos0 = jtp.Pos[0];
    // Find the Z position of the intersection
    tp3d.Pos[2] = (jcs * (iPos0 - iw0) - ics * (jPos0 - jw0)) / den;
    // and the Y position
    bool useI = std::abs(ics) > std::abs(jcs);
    if(useI) {
      tp3d.Pos[1] = (iPos0 - iw0 - isn * tp3d.Pos[2]) / ics;
    } else {
      tp3d.Pos[1] = (jPos0 - jw0 - jsn * tp3d.Pos[2]) / jcs;
    }
    
    // Now find the direction. Protect against large angles first
    if(jtp.Dir[1] == 0) {
      // Going either in the +X direction or -X direction
      if(jtp.Dir[0] > 0) { tp3d.Dir[0] = 1; } else { tp3d.Dir[0] = -1; }
      tp3d.Dir[1] = 0;
      tp3d.Dir[2] = 0;
      return tp3d;
    } // jtp.Dir[1] == 0
    
    tp3d.Wire = iPos0;
    
    // make a copy of itp and shift it by many wires to avoid precision problems
    double itp2_0 = itp.Pos[0] + 100;
    double itp2_1 = itp.Pos[1];
    if(std::abs(itp.Dir[0]) > 0.01) itp2_1 += 100 * itp.Dir[1] / itp.Dir[0];
    // Create a second Point3 for the shifted point
    Point3_t pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = detProp.ConvertTicksToX(itp2_1 / upt, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = detProp.ConvertXToTicks(pos2[0], jPlnID) * upt;
    // Find the wire position (Pos0) in the j plane that this corresponds to
    double jtp2Pos0 = (jtp2Pos1 - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    // Find the Y,Z position using itp2 and jtp2Pos0
    pos2[2] = (jcs * (itp2_0 - iw0) - ics * (jtp2Pos0 - jw0)) / den;
    if(useI) {
      pos2[1] = (itp2_0 - iw0 - isn * pos2[2]) / ics;
    } else {
      pos2[1] = (jtp2Pos0 - jw0 - jsn * pos2[2]) / jcs;
    }
    double sep = PosSep(tp3d.Pos, pos2);
    if(sep == 0) return tp3d;
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) tp3d.Dir[ixyz] = (pos2[ixyz] - tp3d.Pos[ixyz]) /sep;
    tp3d.Flags[kTP3DGood] = true;
    return tp3d;
    
  } // MakeTP3D

  ////////////////////////////////////////////////
  double
  DeltaAngle(const Vector3_t v1, const Vector3_t v2)
  {
    if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) return 0;
    return acos(DotProd(v1, v2));
  }

  ////////////////////////////////////////////////
  Vector3_t
  PointDirection(const Point3_t p1, const Point3_t p2)
  {
    // Finds the direction vector between the two points from p1 to p2
    Vector3_t dir;
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      dir[xyz] = p2[xyz] - p1[xyz];
    if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) return dir;
    if (!SetMag(dir, 1)) {
      dir[0] = 0;
      dir[1] = 0;
      dir[3] = 0;
    }
    return dir;
  } // PointDirection

  //////////////////////////////////////////
  double
  PosSep(const Point3_t& pos1, const Point3_t& pos2)
  {
    return sqrt(PosSep2(pos1, pos2));
  } // PosSep

  //////////////////////////////////////////
  double
  PosSep2(const Point3_t& pos1, const Point3_t& pos2)
  {
    // returns the separation distance^2 between two positions in 3D
    double d0 = pos1[0] - pos2[0];
    double d1 = pos1[1] - pos2[1];
    double d2 = pos1[2] - pos2[2];
    return d0 * d0 + d1 * d1 + d2 * d2;
  } // PosSep2

  //////////////////////////////////////////
  bool
  SetMag(Vector3_t& v1, double mag)
  {
    double den = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
    if (den == 0) return false;
    den = sqrt(den);

    v1[0] *= mag / den;
    v1[1] *= mag / den;
    v1[2] *= mag / den;
    return true;
  } // SetMag

  /////////////////////////////////////////
  void SetDirection(detinfo::DetectorClocksData const& clockData,
                detinfo::DetectorPropertiesData const& detProp, TCSlice& slc, PFPStruct& pfp)
  {
    // Reverse the PFP if needed so that dE/dx increases from start to end
    // Do a simple average of points in the first (last) half

    if(!pfp.Flags[kdEdxDefined]) SetTP3DdEdx(clockData, detProp, slc, pfp);

    float sumNeg = 0, cntNeg = 0, sumPos = 0, cntPos = 0;
    // ignore the end points
    unsigned short halfwayPt = pfp.TP3Ds.size() / 2;
    for(unsigned short ipt = 1; ipt < pfp.TP3Ds.size() - 1; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.Flags[kTP3DGood] || tp3d.TPIndex == USHRT_MAX) continue;
      if(tp3d.dEdx < 0.5) continue;
      if(ipt < halfwayPt) {
        sumNeg += tp3d.dEdx;
        ++cntNeg;
      } else {
        sumPos += tp3d.dEdx;
        ++cntPos;
      }
    } // ipt
    if(cntNeg > 0 && cntPos > 0) {
      sumNeg /= cntNeg;
      sumPos /= cntPos;
      bool prt = (tcc.dbgPFP && pfp.MVI == debug.MVI);
      if(sumNeg > 1.2 * sumPos) {
        if(prt) mf::LogVerbatim("TC")<<"SetDirection found decreasing dE/dx and reversed P"<<pfp.ID;
        Reverse(slc, pfp);
      }
    } // cntNeg > 0 && cntPos > 0
 } // SetDirection
 
  /////////////////////////////////////////
  void SetPFPdEdx(detinfo::DetectorClocksData const& clockData,
                detinfo::DetectorPropertiesData const& detProp,
                const TCSlice& slc,
                PFPStruct& pfp)
  {
    if(pfp.ID <= 0 || pfp.TP3Ds.empty()) return;
    // Fills dE/dx variables in the pfp struct
    if(!pfp.Flags[kdEdxDefined]) SetTP3DdEdx(clockData, detProp, slc, pfp); 

    for (unsigned short end = 0; end < 2; ++end) {
      for (unsigned short plane = 0; plane < slc.nPlanes; ++plane)
        pfp.dEdx[end][plane] = 0;
    } // end

    // square of the maximum length that is used for finding the average dE/dx
    float maxSep2 = 5 * tcc.wirePitch;
    maxSep2 *= maxSep2;

    for (unsigned short end = 0; end < 2; ++end) {
      std::vector<float> cnt(slc.nPlanes);
      short dir = 1 - 2 * end;
      auto endPos = EndTP3D(pfp, end).Pos;
      for (std::size_t ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
        unsigned short ipt;
        if (dir > 0) { ipt = ii; }
        else {
          ipt = pfp.TP3Ds.size() - ii - 1;
        }
        if (ipt >= pfp.TP3Ds.size()) break;
        auto& tp3d = pfp.TP3Ds[ipt];
        if(PosSep2(tp3d.Pos, endPos) > maxSep2) break;
        // require good points
        if(!tp3d.Flags[kTP3DGood]) continue;
        if(tp3d.dEdx < 0.5) continue;
        unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
        pfp.dEdx[end][plane] += tp3d.dEdx;
        ++cnt[plane];
      } // ii
      for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        if (cnt[plane] > 0) {
          pfp.dEdx[end][plane] /= cnt[plane];
        } else {
          pfp.dEdx[end][plane] = 0;
        }
      } // plane
    }   // end

  } // SetPFPdEdx

  /////////////////////////////////////////
  void
  MPV_dEdX(detinfo::DetectorClocksData const& clockData,
               detinfo::DetectorPropertiesData const& detProp,
               const TCSlice& slc,
               PFPStruct& pfp,
               float& dEdXMPV,
               float& dEdXRms)
  {
    // Return a simple average of dE/dx and rms using all points in all planes, not
    // just those at the ends
    dEdXMPV = -1.;
    dEdXRms = -1.;
    if(!pfp.Flags[kdEdxDefined]) SetTP3DdEdx(clockData, detProp, slc, pfp);

    // Poor man's histogram. 50 bins * 1 MeV/cm per bin
    std::vector<int> hist(50, 0);
    unsigned short cnt = 0;
    for(auto& tp3d : pfp.TP3Ds) {
      if(!tp3d.Flags[kTP3DGood]) continue;
      if(tp3d.dEdx > 80.) {
        tp3d.Flags[kTP3DHiDEdx] = true;
        continue;
      } else {
        tp3d.Flags[kTP3DHiDEdx] = false;
      }
      if(tp3d.dEdx < 0.5) continue;
      int bin = std::nearbyint(tp3d.dEdx + 0.5);
      if(bin > 49) bin = 49;
      ++hist[bin];
      ++cnt;
    } // tp3d
    if(cnt < 3) return; 

    int maxBin = 0;
    int maxCnt = 0;

    for(unsigned short bin = 0; bin < hist.size(); ++bin) {
      if(hist[bin] > maxCnt) {
        maxCnt = hist[bin];
        maxBin = bin;
      }
    } // bin
    // construct a maximum dE/dx value for a truncated mean. Set it to
    // the central value of the max bin + 4 * the expected rms (0.15 * maxBin)
    float dedxMax = 4 * maxBin;

    double sum = 0;
    double sum2 = 0;
    double scnt = 0;
    for(auto& tp3d : pfp.TP3Ds) {
      if(!tp3d.Flags[kTP3DGood]) continue;
      if(tp3d.dEdx < 0.5 || tp3d.dEdx > 80.) continue;
      if(tp3d.dEdx < dedxMax) {
        sum += tp3d.dEdx;
        sum2 += tp3d.dEdx * tp3d.dEdx;
        ++scnt;
      } else {
        tp3d.Flags[kTP3DHiDEdx] = true;
      } 
    } // dedx
    if(scnt == 0) return;
    dEdXMPV = sum / cnt;
    // Use a default rms of 30% of the average
    dEdXRms = 0.3 * dEdXMPV;
    double arg = sum2 - scnt * dEdXMPV * dEdXMPV;
    if (arg < 0) return;
    dEdXRms = sqrt(arg) / (scnt - 1);
    // don't return a too-small rms
    double minRms = 0.05 * dEdXMPV;
    if (dEdXRms < minRms) dEdXRms = minRms;
  } // MPV_dEdX

  /////////////////////////////////////////
  void SetTP3DdEdx(detinfo::DetectorClocksData const& clockData,
               detinfo::DetectorPropertiesData const& detProp,
               const TCSlice& slc,
               PFPStruct& pfp) {
    // Calculate dE/dx for each TP3D, put it into the PFPStruct and set the kdEdxDefined flag
    if(pfp.TP3Ds.empty()) return;
    if(pfp.Flags[kdEdxDefined]) return;

    for(auto& tp3d : pfp.TP3Ds) tp3d.dEdx = dEdx(clockData, detProp, slc, tp3d);
    pfp.Flags[kdEdxDefined] = true;
  } // SetTP3DdEdx

  /////////////////////////////////////////
  float
  dEdx(detinfo::DetectorClocksData const& clockData,
       detinfo::DetectorPropertiesData const& detProp,
       const TCSlice& slc,
       const TP3D& tp3d)
  {
    if(tp3d.TPIndex == USHRT_MAX) return 0;
    if(tp3d.TjID > (int)slc.slHits.size()) return 0;
    if(tp3d.TjID <= 0) return 0;

    auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
    if (tp.Environment[kEnvOverlap]) return 0;

    double dQ = TpSumHitChg(slc, tp);
    if(dQ == 0) return 0;
    double time = tp.Pos[1] / tcc.unitsPerTick;
    geo::PlaneID plnID = DecodeCTP(tp.CTP);
    double dx = GetPitch(tp3d);
    if(dx <= 0.) return 0;
    double dQdx = dQ / dx;
    double t0 = 0;
    float dedx = tcc.caloAlg->dEdx_AREA(clockData, detProp, dQdx, time, plnID.Plane, t0);
    if (std::isinf(dedx)) dedx = 0;
    return dedx;
  } // dEdx

  ////////////////////////////////////////////////
  double GetPitch(const TP3D& tp3d)
  {
    geo::PlaneID plnID = DecodeCTP(tp3d.CTP);
    double angleToVert = tcc.geom->Plane(plnID).ThetaZ() - 0.5 * ::util::pi<>();
    double cosgamma = std::abs(std::sin(angleToVert) * tp3d.Dir[1] + std::cos(angleToVert) * tp3d.Dir[2]);
    if(cosgamma < 1.E-5) return -1;
    return tcc.geom->WirePitch(plnID) / cosgamma;
  } // GetPitch

  ////////////////////////////////////////////////
  TP3D CreateTP3D(detinfo::DetectorPropertiesData const& detProp, 
                  const TCSlice& slc, int tjID, unsigned short tpIndex)
  {
    // create a TP3D with a single TP

    TP3D tp3d;
    tp3d.Flags.reset();
    tp3d.Flags[kTP3DBad] = true;
    if(tjID <= 0 || tjID > (int)slc.tjs.size()) return tp3d;
    tp3d.TjID = tjID;
    tp3d.TPIndex = tpIndex;

    auto& tj = slc.tjs[tjID - 1];
    if(tpIndex < tj.EndPt[0] || tpIndex > tj.EndPt[1]) return tp3d;
    auto& tp2 = tj.Pts[tp3d.TPIndex];
    tp3d.CTP = tp2.CTP;
    auto plnID = DecodeCTP(tp2.CTP);
    tp3d.TPX = detProp.ConvertTicksToX(tp2.HitPos[1]/tcc.unitsPerTick, plnID);
    // Get the RMS of the TP in WSE units and convert to cm
    float rms = TPHitsRMSTime(slc, tp2, kAllHits) * tcc.wirePitch;
    // inflate the error for large angle TPs
    if (tp2.AngleCode == 1) rms *= 2;
    // a more careful treatment for long-pulse hits
    if (tp2.AngleCode > 1) {
      std::vector<unsigned int> hitMultiplet;
      for (std::size_t ii = 0; ii < tp2.Hits.size(); ++ii) {
        if (!tp2.UseHit[ii]) continue;
        GetHitMultiplet(slc, tp2.Hits[ii], hitMultiplet, true);
        if (hitMultiplet.size() > 1) break;
      } // ii
      rms = HitsRMSTime(slc, hitMultiplet, kAllHits) * tcc.wirePitch;
      // the returned RMS is closer to the FWHM, so divide by 2
      rms /= 2;
    } // tp2.AngleCode > 1
    tp3d.TPXErr2 = rms * rms;
    tp3d.Wire = tp2.Pos[0];
    // Can't declare it good since Pos and SFIndex aren't defined
    tp3d.Flags[kTP3DGood] = false;
    tp3d.Flags[kTP3DBad] = false;
    return tp3d;
  } // CreateTP3D

  /////////////////////////////////////////
  bool
  SetSection(detinfo::DetectorPropertiesData const& detProp,
             const TCSlice& slc,
             PFPStruct& pfp,
             TP3D& tp3d)
  {
    // Determine which SectionFit this tp3d should reside in, then calculate
    // the 3D position and the distance from the center of the SectionFit

    if(tp3d.Wire < 0) return false;
    // don't change the end points
    if(tp3d.TPIndex == USHRT_MAX) return true;
    if(pfp.SectionFits.empty()) return false;
    if(pfp.SectionFits[0].Pos[0] == -10.0) return false;
    if(pfp.AlgMod[kSmallAng3D]) return true;

    auto plnID = DecodeCTP(tp3d.CTP);

    if (pfp.SectionFits.size() == 1) { tp3d.SFIndex = 0; }
    else {
      // Find the section center that is closest to this point in the wire coordinate
      float best = 1E6;
      for (std::size_t sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        auto& sf = pfp.SectionFits[sfi];
        float sfWire = tcc.geom->WireCoordinate(sf.Pos[1], sf.Pos[2], plnID);
        float sep = std::abs(sfWire - tp3d.Wire);
        if (sep < best) {
          best = sep;
          tp3d.SFIndex = sfi;
        }
      } // sfi
    }   // pfp.SectionFits.size() > 1
    auto& sf = pfp.SectionFits[tp3d.SFIndex];
    auto plnTP = MakeBareTP(detProp, slc, sf.Pos, sf.Dir, tp3d.CTP);
    // the number of wires relative to the SectionFit center
    double dw = tp3d.Wire - plnTP.Pos[0];
    // dt/dW was stored in DeltaRMS
    double t = dw * plnTP.DeltaRMS;
    // define the 3D position
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      tp3d.Pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
    tp3d.along = t;
    return true;
  } // SetSection

  ////////////////////////////////////////////////
  float
  PointPull(const PFPStruct& pfp, const TP3D& tp3d)
  {
    // returns the pull that the tp3d will cause in the pfp section fit. This
    // currently only uses position but eventually will include charge
    return std::abs(tp3d.Pos[0] - tp3d.TPX) / sqrt(tp3d.TPXErr2);
  } // PointPull

  ////////////////////////////////////////////////
  PFPStruct
  CreatePFP(const TCSlice& slc)
  {
    // The calling function should define the size of pfp.TjIDs
    PFPStruct pfp;
    pfp.ID = slc.pfps.size() + 1;
    pfp.ParentUID = 0;
    pfp.TPCID = slc.TPCID;
    // initialize arrays for both ends
    if (slc.nPlanes < 4) {
      pfp.dEdx[0].resize(slc.nPlanes, -1);
      pfp.dEdx[1].resize(slc.nPlanes, -1);
      pfp.dEdxErr[0].resize(slc.nPlanes, -1);
      pfp.dEdxErr[1].resize(slc.nPlanes, -1);
    }
    // create a single section fit to hold the start/end positions and direction
    pfp.SectionFits.resize(1);
    return pfp;
  } // CreatePFP

  /////////////////////////////////////////
  void
  PFPVertexCheck(TCSlice& slc)
  {
    // Ensure that all PFParticles have a start vertex. It is possible for
    // PFParticles to be attached to a 3D vertex that is later killed.
    if (!slc.isValid) return;
    if (slc.pfps.empty()) return;

    for (auto& pfp : slc.pfps) {
      if (pfp.ID == 0) continue;
      if (pfp.Vx3ID[0] > 0) continue;
      if (pfp.SectionFits.empty()) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      vx3.Vx2ID.resize(slc.nPlanes);
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      Point3_t startPos;
      if (pfp.TP3Ds.empty()) {
        // must be a neutrino pfp
        startPos = pfp.SectionFits[0].Pos;
      }
      else if (!pfp.TP3Ds.empty()) {
        // normal pfp
        startPos = pfp.TP3Ds[0].Pos;
      }
      vx3.X = startPos[0];
      vx3.Y = startPos[1];
      vx3.Z = startPos[2];
      vx3.ID = slc.vtx3s.size() + 1;
      vx3.Primary = false;
      ++evt.global3V_UID;
      vx3.UID = evt.global3V_UID;
      slc.vtx3s.push_back(vx3);
      pfp.Vx3ID[0] = vx3.ID;
    } // pfp
  }   // PFPVertexCheck

  /////////////////////////////////////////
  bool
  Store(TCSlice& slc, PFPStruct& pfp)
  {
    // stores the PFParticle in the slice
    bool neutrinoPFP = (pfp.PDGCode == 12 || pfp.PDGCode == 14);
    if (!neutrinoPFP) {
      if (pfp.TjIDs.empty()) return false;
      if (pfp.TP3Ds.size() < 2) return false;
    }
    if(pfp.AlgMod[kSmallAng3D]) {
      // Make the PFP -> TP assn
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.TPIndex != USHRT_MAX) slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex].InPFP = pfp.ID;
      }
    }

    // make some quality checks
    unsigned short nNotSet = 0;
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.Flags[kTP3DBad]) return false;
      if(tp3d.TPIndex == USHRT_MAX) continue;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if (tp.InPFP != pfp.ID) ++nNotSet;
    } // tp3d
    if (nNotSet > 0) return false;
    // check the ID and correct it if it is wrong
    if (pfp.ID != (int)slc.pfps.size() + 1) pfp.ID = slc.pfps.size() + 1;
    ++evt.globalP_UID;
    pfp.UID = evt.globalP_UID;

    // set the 3D match flag
    for (auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      tj.AlgMod[kMat3D] = true;
    } // tjid

    slc.pfps.push_back(pfp);
    return true;
  } // Store

  ////////////////////////////////////////////////
  bool InsideFV(const TCSlice& slc, const PFPStruct& pfp, unsigned short end)
  {
    // returns true if the end of the pfp is inside the fiducial volume of the TPC
    if (pfp.ID <= 0) return false;
    if (end > 1) return false;
    if (pfp.SectionFits.empty()) return false;
    // require that the points are sorted which ensures that the start and end points
    // are the first and last points in the TP3Ds vector
    if (pfp.Flags[kNeedsUpdate]) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;

    float abit = 5;
    Point3_t pos;
    if (neutrinoPFP) { pos = pfp.SectionFits[0].Pos; }
    else if (end == 0) {
      pos = pfp.TP3Ds[0].Pos;
    }
    else {
      pos = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    }
    return (pos[0] > slc.xLo + abit && pos[0] < slc.xHi - abit && pos[1] > slc.yLo + abit &&
            pos[1] < slc.yHi - abit && pos[2] > slc.zLo + abit && pos[2] < slc.zHi - abit);

  } // InsideFV

  ////////////////////////////////////////////////
  bool
  InsideTPC(const Point3_t& pos, geo::TPCID& inTPCID)
  {
    // determine which TPC this point is in. This function returns false
    // if the point is not inside any TPC
    float abit = 5;
    for (const geo::TPCID& tpcid : tcc.geom->IterateTPCIDs()) {
      const geo::TPCGeo& TPC = tcc.geom->TPC(tpcid);
      double local[3] = {0., 0., 0.};
      double world[3] = {0., 0., 0.};
      TPC.LocalToWorld(local, world);
      // reduce the active area of the TPC by a bit to be consistent with FillWireHitRange
      if (pos[0] < world[0] - tcc.geom->DetHalfWidth(tpcid) + abit) continue;
      if (pos[0] > world[0] + tcc.geom->DetHalfWidth(tpcid) - abit) continue;
      if (pos[1] < world[1] - tcc.geom->DetHalfHeight(tpcid) + abit) continue;
      if (pos[1] > world[1] + tcc.geom->DetHalfHeight(tpcid) - abit) continue;
      if (pos[2] < world[2] - tcc.geom->DetLength(tpcid) / 2 + abit) continue;
      if (pos[2] > world[2] + tcc.geom->DetLength(tpcid) / 2 - abit) continue;
      inTPCID = tpcid;
      return true;
    } // tpcid
    return false;
  } // InsideTPC

  ////////////////////////////////////////////////
  void
  FindAlongTrans(Point3_t pos1, Vector3_t dir1, Point3_t pos2, Point2_t& alongTrans)
  {
    // Calculate the distance along and transvers to the direction vector from pos1 to pos2
    alongTrans[0] = 0;
    alongTrans[1] = 0;
    if (pos1[0] == pos2[0] && pos1[1] == pos2[1] && pos1[2] == pos2[2]) return;
    auto ptDir = PointDirection(pos1, pos2);
    SetMag(dir1, 1.0);
    double costh = DotProd(dir1, ptDir);
    if (costh > 1) costh = 1;
    double sep = PosSep(pos1, pos2);
    alongTrans[0] = costh * sep;
    double sinth = sqrt(1 - costh * costh);
    alongTrans[1] = sinth * sep;
  } // FindAlongTrans

  ////////////////////////////////////////////////
  bool
  PointDirIntersect(Point3_t p1,
                    Vector3_t p1Dir,
                    Point3_t p2,
                    Vector3_t p2Dir,
                    Point3_t& intersect,
                    float& doca)
  {
    // Point - vector version
    Point3_t p1End, p2End;
    for (unsigned short xyz = 0; xyz < 3; ++xyz) {
      p1End[xyz] = p1[xyz] + 10 * p1Dir[xyz];
      p2End[xyz] = p2[xyz] + 10 * p2Dir[xyz];
    }
    return LineLineIntersect(p1, p1End, p2, p2End, intersect, doca);
  } // PointDirIntersect

  ////////////////////////////////////////////////
  bool
  LineLineIntersect(Point3_t p1,
                    Point3_t p2,
                    Point3_t p3,
                    Point3_t p4,
                    Point3_t& intersect,
                    float& doca)
  {
    /*
     Calculate the line segment PaPb that is the shortest route between
     two lines P1P2 and P3P4. Calculate also the values of mua and mub where
     Pa = P1 + mua (P2 - P1)
     Pb = P3 + mub (P4 - P3)
     Return FALSE if no solution exists.
     http://paulbourke.net/geometry/pointlineplane/
     */

    Point3_t p13, p43, p21;
    double d1343, d4321, d1321, d4343, d2121;
    double numer, denom;
    constexpr double EPS = std::numeric_limits<double>::min();

    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    if (std::abs(p43[0]) < EPS && std::abs(p43[1]) < EPS && std::abs(p43[2]) < EPS) return (false);
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    if (std::abs(p21[0]) < EPS && std::abs(p21[1]) < EPS && std::abs(p21[2]) < EPS) return (false);

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < EPS) return (false);
    numer = d1343 * d4321 - d1321 * d4343;

    double mua = numer / denom;
    double mub = (d1343 + d4321 * mua) / d4343;

    intersect[0] = p1[0] + mua * p21[0];
    intersect[1] = p1[1] + mua * p21[1];
    intersect[2] = p1[2] + mua * p21[2];
    Point3_t pb;
    pb[0] = p3[0] + mub * p43[0];
    pb[1] = p3[1] + mub * p43[1];
    pb[2] = p3[2] + mub * p43[2];
    doca = PosSep(intersect, pb);
    // average the closest points
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      intersect[xyz] += pb[xyz];
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      intersect[xyz] /= 2;
    return true;
  } // LineLineIntersect

  ////////////////////////////////////////////////
  float
  ChgFracBetween(detinfo::DetectorPropertiesData const& detProp,
                 const TCSlice& slc,
                 Point3_t pos1,
                 Point3_t pos2)
  {
    // Step between pos1 and pos2 and find the fraction of the points that have nearby hits
    // in each plane. This function returns -1 if something is fishy, but this doesn't mean
    // that there is no charge. Note that there is no check for charge precisely at the pos1 and pos2
    // positions
    float sep = PosSep(pos1, pos2);
    if (sep == 0) return -1;
    unsigned short nstep = sep / tcc.wirePitch;
    auto dir = PointDirection(pos1, pos2);
    float sum = 0;
    float cnt = 0;
    TrajPoint tp;
    for (unsigned short step = 0; step < nstep; ++step) {
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        pos1[xyz] += tcc.wirePitch * dir[xyz];
      for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        tp.CTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
        tp.Pos[0] =
          tcc.geom->WireCoordinate(pos1[1], pos1[2], plane, slc.TPCID.TPC, slc.TPCID.Cryostat);
        tp.Pos[1] = detProp.ConvertXToTicks(pos1[0], plane, slc.TPCID.TPC, slc.TPCID.Cryostat) *
                    tcc.unitsPerTick;
        ++cnt;
        if (SignalAtTp(tp)) ++sum;
      } // plane
    }   // step
    if (cnt == 0) return -1;
    return sum / cnt;
  } // ChgFracBetween

  ////////////////////////////////////////////////
  float
  ChgFracNearEnd(detinfo::DetectorPropertiesData const& detProp,
                 const TCSlice& slc,
                 const PFPStruct& pfp,
                 unsigned short end)
  {
    // returns the charge fraction near the end of the pfp. Note that this function
    // assumes that there is only one Tj in a plane.
    if (pfp.ID == 0) return 0;
    if (pfp.TjIDs.empty()) return 0;
    if (end < 0 || end > 1) return 0;
    if (pfp.TPCID != slc.TPCID) return 0;
    if (pfp.SectionFits.empty()) return 0;

    float sum = 0;
    float cnt = 0;
    // keep track of the lowest value and maybe reject it
    float lo = 1;
    float hi = 0;
    auto pos3 = EndTP3D(pfp, end).Pos;
    for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      std::vector<int> tjids(1);
      for (auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if (tj.CTP != inCTP) continue;
        tjids[0] = tjid;
        Point2_t pos2;
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        pos2[0] = tcc.geom->WireCoordinate(pos3[1], pos3[2], planeID);
        if (pos2[0] < -0.4) continue;
        // check for dead wires
        unsigned int wire = std::nearbyint(pos2[0]);
        if (wire > slc.nWires[plane]) continue;
        if (slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
        pos2[1] = detProp.ConvertXToTicks(pos3[0], planeID) * tcc.unitsPerTick;
        float cf = ChgFracNearPos(slc, pos2, tjids);
        if (cf < lo) lo = cf;
        if (cf > hi) hi = cf;
        sum += cf;
        ++cnt;
      } // tjid
    }   // plane
    if (cnt == 0) return 0;
    if (cnt > 1 && lo < 0.3 && hi > 0.8) {
      sum -= lo;
      --cnt;
    }
    return sum / cnt;
  } // ChgFracNearEnd

  ////////////////////////////////////////////////
  TP3D
  EndTP3D(const PFPStruct& pfp, unsigned short end)
  {
    if(end > 1 || pfp.TP3Ds.empty()) {
      TP3D tmp;
      tmp.Flags[kTP3DBad] = true;
      return tmp;
    }
    if(end == 0) return pfp.TP3Ds[0];
    return pfp.TP3Ds[pfp.TP3Ds.size()-1];
  } // TP3DAtEnd

  ////////////////////////////////////////////////
  float
  Length(const PFPStruct& pfp)
  {
    if(pfp.TP3Ds.empty()) return 0;
    // find the first good point
    unsigned short firstGood = USHRT_MAX;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      if(pfp.TP3Ds[ipt].Flags[kTP3DGood]) {
        firstGood = ipt;
        break;
      }
    } // ipt
    if(firstGood == USHRT_MAX) return 0;
    unsigned short lastGood = USHRT_MAX;
    for(unsigned short ipt = pfp.TP3Ds.size() - 1; ipt > 1; --ipt) {
      if(pfp.TP3Ds[ipt].Flags[kTP3DGood]) {
        lastGood = ipt;
        break;
      }
    } // ipt
    if(lastGood == USHRT_MAX) return 0;
    return PosSep(pfp.TP3Ds[firstGood].Pos, pfp.TP3Ds[lastGood].Pos);
  } // Length

  ////////////////////////////////////////////////
  bool
  SectionStartEnd(const PFPStruct& pfp,
                  unsigned short sfIndex,
                  unsigned short& startPt,
                  unsigned short& endPt)
  {
    // this assumes that the TP3Ds vector is sorted
    startPt = USHRT_MAX;
    endPt = USHRT_MAX;
    if (sfIndex >= pfp.SectionFits.size()) return false;

    bool first = true;
    for (std::size_t ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if (tp3d.SFIndex < sfIndex) continue;
      if (first) {
        first = false;
        startPt = ipt;
      } // first
      if (tp3d.SFIndex > sfIndex) break;
      endPt = ipt;
    } // ipt
    return true;

  } // SectionStartEnd

  ////////////////////////////////////////////////
  unsigned short
  FarEnd(const TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos)
  {
    // Returns the end (0 or 1) of the pfp that is furthest away from the position pos
    if (pfp.ID == 0) return 0;
    if (pfp.TP3Ds.empty()) return 0;
    auto& pos0 = pfp.TP3Ds[0].Pos;
    auto& pos1 = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    if (PosSep2(pos1, pos) > PosSep2(pos0, pos)) return 1;
    return 0;
  } // FarEnd

  /////////////////////////////////////////
  int
  PDGCodeVote(detinfo::DetectorClocksData const& clockData,
              detinfo::DetectorPropertiesData const& detProp,
              const TCSlice& slc,
              PFPStruct& pfp)
  {
    // returns a vote using PDG code assignments from dE/dx. A PDGCode of -1 is
    // returned if there was a failure and returns 0 if no decision can be made
    if (pfp.TP3Ds.empty()) return -1;

    // try to do better using dE/dx
    float dEdXAve = 0;
    float dEdXRms = 0;
    MPV_dEdX(clockData, detProp, slc, pfp, dEdXAve, dEdXRms);
    if (dEdXAve < 0) return 0;
    // looks like a proton if dE/dx is high and the rms is low
    dEdXRms /= dEdXAve;
    float length = Length(pfp);
    float mcsmom = 0;
    float chgrms = 0;
    float cnt = 0;
    for (auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      float el = ElectronLikelihood(slc, tj);
      if (el <= 0) continue;
      mcsmom += MCSMom(slc, tj);
      chgrms += tj.ChgRMS;
      ++cnt;
    } // tjid
    if (cnt < 2) return 0;
    mcsmom /= cnt;
    chgrms /= cnt;
    int vote = 0;
    // call anything longer than 150 cm a muon
    if (length > 150) vote = 13;
    // or shorter with low dE/dx and really straight
    if (vote == 0 && length > 50 && dEdXAve < 2.5 && mcsmom > 500) vote = 13;
    // protons have high dE/dx, high MCSMom and low charge rms
    if (vote == 0 && dEdXAve > 3.0 && mcsmom > 200 && chgrms < 0.4) vote = 2212;
    // electrons have low MCSMom and large charge RMS
    if (vote == 0 && mcsmom < 50 && chgrms > 0.4) vote = 11;
    return vote;
  } // PDGCodeVote

  ////////////////////////////////////////////////
  void
  PrintTP3Ds(detinfo::DetectorClocksData const& clockData,
             detinfo::DetectorPropertiesData const& detProp,
             std::string someText,
             const TCSlice& slc,
             const PFPStruct& pfp,
             short printPts)
  {
    if (pfp.TP3Ds.empty()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" pfp P"<<pfp.ID<<" MVI "<<pfp.MVI;
    for(auto tid : pfp.TjIDs) myprt<<" T"<<tid;
    myprt<<" Flags:";
    if(pfp.Flags[kCanSection]) myprt<<" CanSection";
    if(pfp.Flags[kNeedsUpdate]) myprt<<" NeedsUpdate";
    if(pfp.Flags[kStops]) myprt<<" Stops";
    myprt<<" Algs:";
    for(unsigned short ib = 0; ib < pAlgModSize; ++ib) {
      if(pfp.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
    } // ib
    myprt << "\n";
    if (!pfp.SectionFits.empty()) {
      myprt << someText
            << "  SFI ________Pos________   ________Dir_______ _____EndPos________ ChiDOF  NPts "
               "NeedsUpdate?\n";
      for (std::size_t sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        myprt << someText << std::setw(4) << sfi;
        auto& sf = pfp.SectionFits[sfi];
        myprt << std::fixed << std::setprecision(1);
        unsigned short startPt = 0, endPt = 0;
        if (SectionStartEnd(pfp, sfi, startPt, endPt)) {
          auto& start = pfp.TP3Ds[startPt].Pos;
          myprt << std::setw(7) << start[0] << std::setw(7) << start[1] << std::setw(7) << start[2];
        }
        else {
          myprt << " Invalid";
        }
        myprt << std::fixed << std::setprecision(2);
        myprt << std::setw(7) << sf.Dir[0] << std::setw(7) << sf.Dir[1] << std::setw(7)
              << sf.Dir[2];
        myprt << std::fixed << std::setprecision(1);
        if (endPt < pfp.TP3Ds.size()) {
          auto& end = pfp.TP3Ds[endPt].Pos;
          myprt << std::setw(7) << end[0] << std::setw(7) << end[1] << std::setw(7) << end[2];
        }
        else {
          myprt << " Invalid";
        }
        myprt << std::setprecision(1) << std::setw(6) << sf.ChiDOF;
        myprt << std::setw(6) << sf.NPts;
        myprt << std::setw(6) << sf.NeedsUpdate;
        myprt << "\n";
      } // sec
    }   // SectionFits
    if (printPts < 0) {
      // print the head if we print all points
      myprt<<someText<<" Note: GBH = TP3D Flags. G = Good, B = Bad, H = High dE/dx \n";
      myprt<<someText<<"  ipt SFI ________Pos________  Delta Pull  GBH   Path  along dE/dx S?    T_ipt_P:W:T\n";
    }
    unsigned short fromPt = 0;
    unsigned short toPt = pfp.TP3Ds.size() - 1;
    if (printPts >= 0) fromPt = toPt;
    // temp kink angle for each point
    std::vector<float> dang(pfp.TP3Ds.size(), -1);
    for (unsigned short ipt = fromPt; ipt <= toPt; ++ipt) {
      auto tp3d = pfp.TP3Ds[ipt];
      myprt << someText << std::setw(4) << ipt;
      myprt << std::setw(4) << tp3d.SFIndex;
      myprt << std::fixed << std::setprecision(1);
      myprt << std::setw(7) << tp3d.Pos[0] << std::setw(7) << tp3d.Pos[1] << std::setw(7)
            << tp3d.Pos[2];
      myprt << std::setprecision(1) << std::setw(6) << (tp3d.Pos[0] - tp3d.TPX);
      float pull = PointPull(pfp, tp3d);
      myprt<<std::setprecision(1)<<std::setw(6)<<pull;
      myprt<<std::setw(3)<<tp3d.Flags[kTP3DGood]<<tp3d.Flags[kTP3DBad]<<tp3d.Flags[kTP3DHiDEdx];
      myprt<<std::setw(7)<<std::setprecision(2)<<PosSep(tp3d.Pos, pfp.TP3Ds[0].Pos);
      myprt<<std::setw(7)<<std::setprecision(1)<<tp3d.along;
      myprt<<std::setw(6)<<std::setprecision(2)<<tp3d.dEdx;
      // print SignalAtTP in each plane
      myprt << " ";
      for (unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        auto tp = MakeBareTP(detProp, slc, tp3d.Pos, inCTP);
        myprt<<SignalAtTp(tp);
      } // plane
      if(tp3d.TPIndex != USHRT_MAX) {
        if(tp3d.TjID > 0) {
          auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
          myprt<<" T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"_"<<PrintPos(slc, tp)<<" "<<TPEnvString(tp);
        } else {
          myprt<<" UNDEFINED";
        }
      } // tp3d.TPIndex != USHRT_MAX
      if(tp3d.Flags[kTP3DHiDEdx]) myprt<<" Hi_dEdx";
      myprt<<"\n";
    } // ipt
  }   // PrintTP3Ds
} // namespace
