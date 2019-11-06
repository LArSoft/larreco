#include "larreco/RecoAlg/TCAlg/PFPUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TDecompSVD.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <limits.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <array>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/StepUtils.h"
#include "larreco/RecoAlg/TCAlg/TCShower.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {

  struct SortEntry{
    unsigned int index;
    float val;
  };
  // TODO: Fix the sorting mess
  bool valDecreasings (SortEntry c1, SortEntry c2) { return (c1.val > c2.val);}
  bool valIncreasings (SortEntry c1, SortEntry c2) { return (c1.val < c2.val);}

  /////////////////////////////////////////
  void StitchPFPs()
  {
    // Stitch PFParticles in different TPCs. This does serious damage to PFPStruct and should
    // only be called from TrajCluster module just before making PFParticles to put in the event
    if(slices.size() < 2) return;
    if(tcc.geom->NTPC() == 1) return;
    if(tcc.pfpStitchCuts.size() < 2) return;
    if(tcc.pfpStitchCuts[0] <= 0) return;

    bool prt = tcc.dbgStitch;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      std::string fcnLabel = "SP";
      myprt<<fcnLabel<<" cuts "<<sqrt(tcc.pfpStitchCuts[0])<<" "<<tcc.pfpStitchCuts[1]<<"\n";
      bool printHeader = true;
      for(size_t isl = 0; isl < slices.size(); ++isl) {
        if(debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if(slc.pfps.empty()) continue;
        for(auto& pfp : slc.pfps) PrintP(fcnLabel, myprt, pfp, printHeader);
      } // slc
    } // prt

    // lists of pfp UIDs to stitch
    std::vector<std::vector<int>> stLists;
    for(unsigned short sl1 = 0; sl1 < slices.size() - 1; ++sl1) {
      auto& slc1 = slices[sl1];
      for(unsigned short sl2 = sl1 + 1; sl2 < slices.size(); ++sl2) {
        auto& slc2 = slices[sl2];
        // look for PFParticles in the same recob::Slice
        if(slc1.ID != slc2.ID) continue;
        for(auto& p1 : slc1.pfps) {
          if(p1.ID <= 0) continue;
          // Can't stitch shower PFPs
          if(p1.PDGCode == 1111) continue;
          for(auto& p2 : slc2.pfps) {
            if(p2.ID <= 0) continue;
            // Can't stitch shower PFPs
            if(p2.PDGCode == 1111) continue;
            float maxSep2 = tcc.pfpStitchCuts[0];
            float maxCth = tcc.pfpStitchCuts[1];
            bool gotit = false;
            for(unsigned short e1 = 0; e1 < 2; ++e1) {
              auto pos1 = PosAtEnd(p1, e1);
              // require the end to be close to a TPC boundary
              if(InsideFV(slc1, p1, e1)) continue;
              auto dir1 = DirAtEnd(p1, e1);
              for(unsigned short e2 = 0; e2 < 2; ++e2) {
                auto pos2 = PosAtEnd(p2, e2);
                // require the end to be close to a TPC boundary
                if(InsideFV(slc2, p2, e2)) continue;
                auto dir2 = DirAtEnd(p2, e2);
                float sep = PosSep2(pos1, pos2);
                if(sep > maxSep2) continue;
                float cth = std::abs(DotProd(dir1, dir2));
                if(cth < maxCth) continue;
                maxSep2 = sep;
                maxCth = cth;
                gotit = true;
              } // e2
            } // e1
            if(!gotit) continue;
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<"Stitch slice "<<slc1.ID<<" P"<<p1.UID<<" TPC "<<p1.TPCID.TPC;
              myprt<<" and P"<<p2.UID<<" TPC "<<p2.TPCID.TPC;
              myprt<<" sep "<<sqrt(maxSep2)<<" maxCth "<<maxCth;
            }
            // see if either of these are in a list
            bool added = false;
            for(auto& pm : stLists) {
              bool p1InList = (std::find(pm.begin(), pm.end(), p1.UID) != pm.end());
              bool p2InList = (std::find(pm.begin(), pm.end(), p2.UID) != pm.end());
              if(p1InList || p2InList) {
                if(p1InList) pm.push_back(p2.UID);
                if(p2InList) pm.push_back(p1.UID);
                added = true;
              }
            } // pm
            if(added) continue;
            // start a new list
            std::vector<int> tmp(2);
            tmp[0] = p1.UID;
            tmp[1] = p2.UID;
            stLists.push_back(tmp);
            break;
          } // p2
        } // p1
      } // sl2
    } // sl1
    if(stLists.empty()) return;

    for(auto& stl : stLists) {
      // Find the endpoints of the stitched pfp
      float minZ = 1E6;
      std::pair<unsigned short, unsigned short> minZIndx;
      unsigned short minZEnd = 2;
      for(auto puid : stl) {
        auto slcIndex = GetSliceIndex("P", puid);
        if(slcIndex.first == USHRT_MAX) continue;
        auto& pfp = slices[slcIndex.first].pfps[slcIndex.second];
        for(unsigned short end = 0; end < 2; ++end) {
          auto pos = PosAtEnd(pfp, end);
          if(pos[2] < minZ) { minZ = pos[2]; minZIndx = slcIndex;  minZEnd = end; }
        } // end
      } // puid
      if(minZEnd > 1) continue;
      // preserve the pfp with the min Z position
      auto& pfp = slices[minZIndx.first].pfps[minZIndx.second];
      if(prt) mf::LogVerbatim("TC")<<"SP: P"<<pfp.UID;
      // reverse it if necessary
//      if(minZEnd != 0) ReversePFP(slices[minZIndx.first], pfp);
      // add the Tjs in the other slices to it
      for(auto puid : stl) {
        if(puid == pfp.UID) continue;
        auto sIndx = GetSliceIndex("P", puid);
        if(sIndx.first == USHRT_MAX) continue;
        auto& opfp = slices[sIndx.first].pfps[sIndx.second];
        if(prt) mf::LogVerbatim("TC")<<" +P"<<opfp.UID;
        pfp.TjUIDs.insert(pfp.TjUIDs.end(), opfp.TjUIDs.begin(), opfp.TjUIDs.end());
        if(prt) mf::LogVerbatim();
        // Check for parents and daughters
        if(opfp.ParentUID > 0) {
          auto pSlcIndx = GetSliceIndex("P", opfp.ParentUID);
          if(pSlcIndx.first < slices.size()) {
            auto& parpfp = slices[pSlcIndx.first].pfps[pSlcIndx.second];
            std::replace(parpfp.DtrUIDs.begin(), parpfp.DtrUIDs.begin(), opfp.UID, pfp.UID);
          } // valid pSlcIndx
        } // has a parent
        for(auto dtruid : opfp.DtrUIDs) {
          auto dSlcIndx = GetSliceIndex("P", dtruid);
          if(dSlcIndx.first < slices.size()) {
            auto& dtrpfp = slices[dSlcIndx.first].pfps[dSlcIndx.second];
            dtrpfp.ParentUID = pfp.UID;
          } // valid dSlcIndx
        } // dtruid
        // declare it obsolete
        opfp.ID = 0;
      } // puid
    } // stl

  } // StitchPFPs

  /////////////////////////////////////////
  void FindSptPFParticles(TCSlice& slc)
  {
    // Find 3D Tj matches using SpacePoints
    if(!evt.sptHandle) return;
    // ensure that allHitsSptIndex was sized correctly
    if(evt.allHitsSptIndex.size() != (*evt.allHits).size()) return;
    // This code will choke if there are too many spts
    if(evt.allHitsSptIndex.size() > INT_MAX) return;
    
    // Create a list of spt -> Tjs in three planes
    std::vector<std::vector<int>> sptAssns;
    sptAssns.resize((*evt.sptHandle).size(), std::vector<int>(3, INT_MAX));

    unsigned int nspts = (*evt.sptHandle).size() - 1;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
          if(evt.allHitsSptIndex[ahi] > nspts) continue;
          sptAssns[evt.allHitsSptIndex[ahi]][plane] = tj.ID;
        } // ii
      } // ipt
    } // tj
    
    bool prt = (tcc.dbgPFP && tcc.dbgSlc);
    
    std::vector<MatchStruct> matVec;
    MatchStruct ms;
    for(auto& sptAssn : sptAssns) {
      // require a triple match
      if(sptAssn[0] == INT_MAX || sptAssn[1] == INT_MAX || sptAssn[2] == INT_MAX) continue;
      // look for this triplet in matVec
      unsigned short indx = 0;
      for(indx = 0; indx < matVec.size(); ++indx) if(matVec[indx].TjIDs == sptAssn) break;
      if(indx == matVec.size()) {
        ms.TjIDs = sptAssn;
        matVec.push_back(ms);
      }
      ++matVec[indx].Count;
    } // isp
    
    // sort by decreasing count
    std::vector<SortEntry> sortVec;
    for(unsigned short indx = 0; indx < matVec.size(); ++indx) {
      auto& ms = matVec[indx];
      // count the number of TPs in all Tjs
      float tpCnt = 0;
      for(auto tid : ms.TjIDs) {
        auto& tj = slc.tjs[tid - 1];
        tpCnt += NumPtsWithCharge(slc, tj, false);
      } // tid
      float frac = ms.Count / tpCnt;
      // ignore matches with a very low match fraction
      if(frac < 0.01) continue;
      SortEntry se;
      se.index = indx;
      se.val = matVec[indx].Count;
      sortVec.push_back(se);
    } // ii
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
    std::vector<MatchStruct> tmp(sortVec.size());
    for(unsigned short indx = 0; indx < sortVec.size(); ++indx) tmp[indx] = matVec[sortVec[indx].index];
    matVec = tmp;
    tmp.resize(0);
    sortVec.resize(0);
    
    if(matVec.empty()) return;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MVI  Count  Tjs\n";
      for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
        auto& ms = matVec[indx];
        myprt<<std::setw(5)<<indx<<std::setw(6)<<(int)ms.Count;
        for(auto tid : ms.TjIDs) myprt<<" T"<<tid;
        // count the number of TPs in all Tjs
        float tpCnt = 0;
        for(auto tid : ms.TjIDs) {
          auto& tj = slc.tjs[tid - 1];
          tpCnt += NumPtsWithCharge(slc, tj, false);
        } // tid
        float frac = ms.Count / tpCnt;
        myprt<<" matFrac "<<std::fixed<<std::setprecision(3)<<frac;
        myprt<<"\n";
      } // indx
    } // prt
    
    MakePFParticles(slc, matVec, 0);
    
    // a last debug print
    if(tcc.dbgPFP && debug.MVI != UINT_MAX) {
      for(auto& pfp : slc.pfps) if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("FSPFP", slc, pfp, -1);
    } // last debug print

    
  } // FindSptPFParticles

  /////////////////////////////////////////
  void FindPFParticles(TCSlice& slc)
  {
    // Match Tjs in 3D and create PFParticles
    
    if(tcc.match3DCuts[0] <= 0) return;
   
    // clear the TP -> P assn Tjs so that all are considered
    for(auto& tj : slc.tjs) {
      for(auto& tp : tj.Pts) tp.InPFP = 0;
    } // tj
    
    bool prt = (tcc.dbgPFP && tcc.dbgSlc);

    // Match these points in 3D
    std::vector<MatchStruct> matVec;
    
    // iterate twice (at most), looking for 3-plane matches in 3-plane TPCs on the
    // first iteration and 2-plane matches + dead regions in 3-plane TPCs on the second
    
    unsigned short maxNit = 2;
    if(slc.nPlanes == 2) maxNit = 1;
    if(std::nearbyint(tcc.match3DCuts[2]) == 0) maxNit = 1;
    for(unsigned short nit = 0; nit < maxNit; ++nit) {
      // fill the mAllTraj vector with TPs that aren't matched in 3D
      FillmAllTraj(slc);
      if(slc.mallTraj.empty()) break;
      matVec.clear();
      if(slc.nPlanes == 3 && nit == 0) {
        // look for match triplets
        Match3Planes(slc, matVec);
      } else {
        // look for match doublets requiring a dead region in the 3rd plane for 3-plane TPCs
        Match2Planes(slc, matVec);
      }
      if(matVec.empty()) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"nit "<<nit<<" MVI  Count  Tjs\n";
        for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
          auto& ms = matVec[indx];
          myprt<<std::setw(5)<<indx<<std::setw(6)<<(int)ms.Count;
          for(auto tid : ms.TjIDs) myprt<<" T"<<tid;
          // count the number of TPs in all Tjs
          float tpCnt = 0;
          for(auto tid : ms.TjIDs) {
            auto& tj = slc.tjs[tid - 1];
            tpCnt += NumPtsWithCharge(slc, tj, false);
          } // tid
          float frac = ms.Count / tpCnt;
          myprt<<" matFrac "<<std::fixed<<std::setprecision(3)<<frac;
          myprt<<"\n";
        } // indx
      } // prt
      MakePFParticles(slc, matVec, nit);
    } // nit
    
    // reconcile TP -> P assns in all pfps in this slice
//    ReconcileTPs(slc);
    
    // a last debug print
    if(tcc.dbgPFP && debug.MVI != UINT_MAX) {
      for(auto& pfp : slc.pfps) if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("FPFP", slc, pfp, -1);
    } // last debug print
    
    slc.mallTraj.resize(0);
    
  } // FindPFParticles
  
  ////////////////////////////////////////////////
  void MakePFParticles(TCSlice& slc, std::vector<MatchStruct> matVec, unsigned short matVec_Iter)
  {
    // Makes PFParticles using Tjs listed in matVec
    if(matVec.empty()) return;

    bool prt = (tcc.dbgPFP && tcc.dbgSlc);

    // create a PFParticle for each valid match combination
    // define a glabal matVec index for debugging
    for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
      // tone down the level of printing in ReSection
      bool foundMVI = (tcc.dbgPFP && indx == debug.MVI && matVec_Iter == debug.MVI_Iter);
      if(foundMVI) prt = true;
      auto& ms = matVec[indx];
      if(foundMVI) {
        std::cout<<"found MVI "<<indx<<" in MakePFParticles ms.Count = "<<ms.Count<<"\n";
      }
      // ignore dead matches
      if(ms.Count == 0) continue;
      // count the number of TPs that are available (not already 3D-matched) and used in a pfp
      float npts = 0;
      for(unsigned short itj = 0; itj < ms.TjIDs.size(); ++itj) {
        auto& tj = slc.tjs[ms.TjIDs[itj] - 1];
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) if(tj.Pts[ipt].InPFP == 0) ++npts;
      } // tjID
      int pdgCode = PDGCodeVote(slc, ms.TjIDs, prt);
      // Create a vector of PFPs for this match so that we can split it later on if a kink is found
      std::vector<PFPStruct> pfpVec(1);
      pfpVec[0] = CreatePFP(slc);
      // Define the starting set of tjs that were matched. TPs from other tjs may be added later
      pfpVec[0].TjIDs = ms.TjIDs;
      pfpVec[0].MVI = indx;
      // fill the TP3D points using the 2D trajectory points for Tjs in TjIDs. All
      // points are put in one section
      MakeTP3Ds(slc, pfpVec[0]);
      // fit all the points to get the general direction
      if(!FitSection(slc, pfpVec[0], 0)) continue;
      if(pfpVec[0].SectionFits[0].ChiDOF > 200) {
        std::cout<<"Found a crazy high ChiDOF P"<<pfpVec[0].ID<<" "<<pfpVec[0].SectionFits[0].ChiDOF<<"\n";
        continue;
      }
      // sort the points by the distance along the general direction vector
      if(!SortSection(pfpVec[0], 0)) continue;
      // define a junk pfp to be short with low MCSMom. These are likely to be shower-like
      // pfps. A simple 3D line fit will be done. No attempt will be made to reconstruct it
      // in sections or to look for kinks
      npts = pfpVec[0].TP3Ds.size();
      pfpVec[0].AlgMod[kJunk3D] = (npts < 20 && MCSMom(slc, pfpVec[0].TjIDs) < 50) || (npts < 10);
      if(prt) {
        auto& pfp = pfpVec[0];
        mf::LogVerbatim myprt("TC");
        myprt<<" indx "<<matVec_Iter<<"/"<<indx<<" Count "<<std::setw(5)<<(int)ms.Count;
        myprt<<" P"<<pfpVec[0].ID;
        myprt<<" ->";
        for(auto& tjid : pfp.TjIDs) myprt<<" T"<<tjid;
        myprt<<" projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, pfp.SectionFits[0].Pos, pfp.SectionFits[0].Dir, inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, pfp.TjIDs);
        myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(slc, pfp.TjIDs, false);
        myprt<<" nTP3Ds "<<pfp.TP3Ds.size();
        myprt<<" Reco3DRange "<<Find3DRecoRange(slc, pfp, 0, (unsigned short)tcc.match3DCuts[3], 1);
      } // prt
      if(foundMVI) {
        PrintTP3Ds("FF", slc, pfpVec[0], -1);
      }
      // Pass the pfp vector to SplitAtKinks which may produce pfpVec[1]
      SplitAtKinks(slc, pfpVec, foundMVI);
      bool success = true;
      for(unsigned short ip = 0; ip < pfpVec.size(); ++ip) {
        auto& pfp = pfpVec[ip];
        // set the end flag bits
        geo::TPCID tpcid;
        for(unsigned short end = 0; end < 2; ++end) {
          // first set them all to 0
          pfp.EndFlag[end].reset();
          auto pos = PosAtEnd(pfp, end);
          if(!InsideTPC(pos, tpcid)) pfp.EndFlag[end][kOutFV] = true;
        } // end
        // Set kink flag and create a vertex between this pfp and the previous one that was stored
        if(ip > 0) {
          pfp.EndFlag[0][kAtKink] = true;
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
          ++evt.globalP_UID;
          vx3.UID = evt.globalP_UID;
          slc.vtx3s.push_back(vx3);
          pfp.Vx3ID[0] = vx3.ID;
          auto& prevPFP = slc.pfps[slc.pfps.size() - 1];
          prevPFP.Vx3ID[1] = vx3.ID;
        } // ip > 0
        // remove really bad TP3Ds
        KillBadPoints(slc, pfp, 50., foundMVI);
        // check for a valid two-plane match with a Tj in the third plane for long pfps.
        // For short pfps, it is possible that a Tj would be too short to be reconstructed
        // in the third plane.
        if(pfp.TjIDs.size() == 2 && slc.nPlanes == 3 && pfp.TP3Ds.size() > 20 && !ValidTwoPlaneMatch(slc, pfp)) {
//          std::cout<<"  Invalid two-plane match\n";
          continue;
        }
        // Skip this combination if it isn't reconstructable in 3D
        if(Find3DRecoRange(slc, pfp, 0, (unsigned short)tcc.match3DCuts[3], 1) == USHRT_MAX) continue;
        // See if it possible to reconstruct in more than one section
        pfp.Flags[kCanSection] = CanSection(slc, pfp);
           // Do a fit in multiple sections if the initial fit is poor
           if(pfpVec[0].SectionFits[0].ChiDOF < tcc.match3DCuts[5]) {
          // Good fit with one section
          pfpVec[0].Flags[kNeedsUpdate] = false;
        } else if(pfpVec[0].Flags[kCanSection]) {
          if(!ReSection(slc, pfpVec[0], foundMVI)) continue;
          KillBadPoints(slc, pfpVec[0], tcc.match3DCuts[4], foundMVI);
          // Try to remove bad points if there was a ReSection problem
/*
          if(pfpVec[0].SectionFits.size() == 1 && pfpVec[0].SectionFits[0].ChiDOF > tcc.match3DCuts[5] && !pfpVec[0].AlgMod[kJunk3D]) {
            KillBadPoints(slc, pfpVec[0], tcc.match3DCuts[4], foundMVI);
            // try again
            if(pfpVec[0].Flags[kNeedsUpdate]) {
              pfpVec[0].Flags[kCanSection] = true;
              ReSection(slc, pfpVec[0], foundMVI);
            }
          } // ReSection problem
*/
        } // CanSection
        if(foundMVI) {
          PrintTP3Ds("RS", slc, pfpVec[0], -1);
        }
        pfpVec[0].PDGCode = pdgCode;
        // FillGaps3D looks for gaps in the TP3Ds vector caused by broken trajectories and
        // inserts new TP3Ds if there are hits in the gaps. This search is only done in a
        // plane if the projection of the pfp results in a large angle where 2D reconstruction
        // is likely to be poor - not true for TCWork2
        FillGaps3D(slc, pfp, foundMVI);
        // Trim points from the ends until there is a 3D point where there is a signal in at least two planes
        TrimEndPts(slc, pfpVec[0], foundMVI);
        // Check the TP3D -> TP assn, resolve conflicts and set TP -> InPFP
        if(!ReconcileTPs(slc, pfp, foundMVI)) continue;
        // Look for mis-placed 2D and 3D vertices
        ReconcileVertices(slc, pfpVec[0], foundMVI);
        FilldEdx(slc, pfp);
        if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("STORE", slc, pfp, -1);
        if(!StorePFP(slc, pfp)) {
          std::cout<<"MPFP: StorePFP failed\n";
          success = false;
          break;
        }
        // Oct30, 2019. It seems like a bad idea to reject a tj match combination if only
        // a few TPs were matched prior to this point. 
        if(!tcc.useAlg[kTCWork2]) {
          // clobber later entries that have these Tjs
          for(unsigned short jndx = indx + 1; jndx < matVec.size(); ++jndx) {
            for(auto tid : pfp.TjIDs) {
              auto& jms = matVec[jndx];
              if(jms.Count <= 0) continue;
              if(std::find(jms.TjIDs.begin(), jms.TjIDs.end(), tid) != jms.TjIDs.end()) jms.Count = 0;
            } // tid
          } // jndx        
        } // !tcc.useAlg[kTCWork2]
        if(!success) {
          std::cout<<"MPFP: Failed to split and store PFPs. Write some recovery code\n";
        }
      } // ip (iterate over split pfps)
    } // indx (iterate over matchVec entries)
    slc.mallTraj.resize(0);
    
  } // MakePFParticles

  ////////////////////////////////////////////////
  void ChkPFPMC(TCSlice& slc, PFPStruct& pfp)
  {
    // This function is used to decide what ChiDOF cut should be made to reject
    // invalid 3D matches 
    if(evt.allHitsMCPIndex.empty()) return;
    if(pfp.SectionFits.size() != 1 || pfp.TP3Ds.empty()) {
      std::cout<<"ChkPFPMC: Something wrong with P"<<pfp.ID<<"\n";
      return;
    }
    
    // mcpIndex and count
    std::vector<std::pair<unsigned int, unsigned short>> mcpi_cnt;
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.TjID <= 0) {
        std::cout<<"oops\n";
        exit(1);
      }
      unsigned int mcpIndex = UINT_MAX;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
        mcpIndex = evt.allHitsMCPIndex[ahi];
        break;
      } // ii
      if(mcpIndex == UINT_MAX) continue;
      // look for it in the list
      unsigned short indx = 0;
      for(indx = 0; indx < mcpi_cnt.size(); ++indx) if(mcpi_cnt[indx].first == mcpIndex) break;
      // not found so add it
      if(indx == mcpi_cnt.size()) mcpi_cnt.push_back(std::make_pair(mcpIndex, 0));
      ++mcpi_cnt[indx].second;
    } // tp3d
    auto& sf = pfp.SectionFits[0];
    std::cout<<"ChkPFPMC: P"<<pfp.ID;
    std::cout<<std::setprecision(1)<<" ChiDOF "<<sf.ChiDOF;
    std::cout<<" MCP_cnt";
    for(auto mc : mcpi_cnt) std::cout<<" "<<mc.first<<"_"<<mc.second;
    std::cout<<"\n";

  } // ChkPFPMC

  ////////////////////////////////////////////////
  bool ReconcileTPs(TCSlice& slc, PFPStruct& pfp, bool prt)
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
      if(tp3d.IsBad) continue;
      if(tp3d.TjID <= 0) {
        std::cout<<"RTPs3D P"<<pfp.ID<<" found an invalid TP3D -> Tj assn: tp3d.TjID = "<<tp3d.TjID<<"\n";
        return false;
      }
      // compare the TP3D -> TP -> P assn with the P -> TP assn
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.InPFP > 0 && tp.InPFP != pfp.ID) {
        std::cout<<"RTPs3D P"<<pfp.ID<<" existing TP -> P assn "<<tp.InPFP<<" isn't this pfp\n";
        return false;
      }
      // find the (Tj ID, TP count) pair in the list
      unsigned short indx = 0;
      for(indx = 0; indx < tjTPCnt.size(); ++indx) if(tjTPCnt[indx].first == tp3d.TjID) break;
      if(indx == tjTPCnt.size()) tjTPCnt.push_back(std::make_pair(tp3d.TjID, 0));
      ++tjTPCnt[indx].second;
      // make the TP -> P assn
      tp.InPFP = pfp.ID;
    } // tp3d

    std::vector<int> nTjIDs;
    for(auto& tjtpcnt : tjTPCnt) {
      auto& tj = slc.tjs[tjtpcnt.first - 1];
      float npwc = NumPtsWithCharge(slc, tj, false);
      if(prt) std::cout<<" T"<<tjtpcnt.first<<" npwc "<<npwc<<" cnt "<<tjtpcnt.second<<"\n";
      if(tjtpcnt.second > 0.8 * npwc) nTjIDs.push_back(tjtpcnt.first);
    } // tjtpcnt
    // TODO: is this really a failure?
    if(prt) {
      mf::LogVerbatim("TC")<<"RTPs3D: P"<<pfp.ID<<" nTjIDs "<<nTjIDs.size();
    }
    if(nTjIDs.size() < 2) {
      return false;
    }
    pfp.TjIDs = nTjIDs;
    
    return true;
  } // ReconcileTPs

  ////////////////////////////////////////////////
  void ReconcileTPs(TCSlice& slc)
  {
    // Reconciles TP ownership conflicts between PFParticles
    // Make a one-to-one TP -> P assn and look for one-to-many assns.
    // Note: Comparing the pulls for a TP to two different PFParticles generally results
    // in selecting the first PFParticle that was made which is not too surprising considering
    // the order in which they were created. This comparison has been commented out in favor
    // of simply keeping the old assn and removing the new one by setting IsBad true.
    
//    bool prt = false;
    if(!tcc.useAlg[kRTPs3D]) return;

    // make a list of T -> P assns
    std::vector<int> TinP;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(tp3d.TjID <= 0) continue;
        if(std::find(TinP.begin(), TinP.end(), tp3d.TjID) == TinP.end()) TinP.push_back(tp3d.TjID);
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if(tp.InPFP > 0) {
          // an assn exists. Set the overlap bit and check consistency
          tp.Environment[kEnvOverlap] = true;
/*
          auto& oldp = slc.pfps[tp.InPFP - 1];
          // find the TP3D index
          unsigned short otp = 0;
          for(otp = 0; otp < oldp.TP3Ds.size(); ++otp) {
            auto& otp3d = oldp.TP3Ds[otp];
            if(otp3d.TjID == tp3d.TjID && otp3d.TPIndex == tp3d.TPIndex) break;
          }
          auto& otp3d = oldp.TP3Ds[otp];
*/
          // keep the old assn and remove the new one
          tp3d.IsBad = true;
          tp3d.IsGood = false;
          tp.InPFP = 0;
/*
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<"RTPs: TP "<<PrintPos(slc, tp)<<" one-to-many -> P"<<tp.InPFP;
            myprt<<" pull "<<PointPull(oldp, otp3d);
            myprt<<" and P"<<pfp.ID;
            myprt<<" pull "<<PointPull(pfp, tp3d)<<". Keeping the first";
          } // prt
*/
        } else {
          // no assn exists
          tp.InPFP = pfp.ID;
        } // tp.InPFP > 0
      } // ipt
    } // pfp
    
  } // ReconcileTPs

  /////////////////////////////////////////
  void MakePFPTjs(TCSlice& slc)
  {
    // This function clobbers all of the tjs that are used in TP3Ds in the pfp and replaces
    // them with new tjs that have a consistent set of TPs to prepare for putting them
    // into the event. Note that none of the Tjs are attached to 2D vertices.
    if(!tcc.useAlg[kMakePFPTjs]) return;
    
    // kill trajectories
    std::vector<int> killme;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.TjID <= 0) continue;
        if(tp3d.IsBad) continue;
        if(std::find(killme.begin(), killme.end(), tp3d.TjID) == killme.end()) killme.push_back(tp3d.TjID);
      } // tp3d
    } // pfp
    
    for(auto tid : killme) MakeTrajectoryObsolete(slc, (unsigned int)(tid - 1));

    // Make template trajectories in each plane. These will be re-used by
    // each PFParticle
    std::vector<Trajectory> ptjs(slc.nPlanes);
    // define the basic tj variables
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
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
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      // initialize the tjs
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        ptjs[plane].Pts.clear();
        --evt.WorkID;
        if(evt.WorkID == INT_MIN) evt.WorkID = -1;
        ptjs[plane].ID = evt.WorkID;
      } // plane
      pfp.TjIDs.clear();
      // iterate through all of the TP3Ds, adding TPs to the TJ in the appropriate plane.
      // The assumption here is that TP order reflects the TP3D order
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.TjID <= 0) continue;
        if(tp3d.IsBad) continue;
        // Get a reference to the 2D TP
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if(tp.Chg <= 0) continue;
        if(tp.InPFP > 0 && (int)tp.InPFP != pfp.ID) continue;
        tp.InPFP = pfp.ID;
        unsigned short plane = DecodeCTP(tp.CTP).Plane;
        // append it to Pts
        ptjs[plane].Pts.push_back(tp);
      } // tp3d
      // finish defining each of the Tjs and store them
      for(auto& tj : ptjs) {
        if(tj.Pts.size() < 2) continue;
        tj.PDGCode = pfp.PDGCode;
        tj.MCSMom = MCSMom(slc, tj);
        if(!StoreTraj(slc, tj)) continue;
        // associate it with the pfp
        auto& newTj = slc.tjs[slc.tjs.size() - 1];
        pfp.TjIDs.push_back(newTj.ID);
      } // tj
    } // pfp
  } // MakePFPTjs
  
  /////////////////////////////////////////
  void Match3Planes(TCSlice& slc, std::vector<MatchStruct>& matVec)
  {
    // A simpler and faster version of MatchPlanes that only creates three plane matches
    
    if(slc.mallTraj.empty()) return;
    if(slc.nPlanes != 3) return;
    
    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;

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
    if(tcc.match3DCuts[1] < (float)USHRT_MAX) maxCnt = (unsigned short)tcc.match3DCuts[1];
    // a list of those Tjs
    std::vector<unsigned short> tMaxed;
    
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // see if we hit the maxCnt limit
      if(std::find(tMaxed.begin(), tMaxed.end(), iTjPt.id) != tMaxed.end()) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      unsigned short iPlane = iTjPt.plane;
      unsigned int iWire = itp.Pos[0];
      tIDs[iPlane] = iTjPt.id;
      bool hitMaxCnt = false;
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.plane == iTjPt.plane) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large
        if(jTjPt.xlo > iTjPt.xhi + xcut) break;
        // see if we hit the maxCnt limit
        if(std::find(tMaxed.begin(), tMaxed.end(), jTjPt.id) != tMaxed.end()) continue;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jPlane = jTjPt.plane;
        unsigned int jWire = jtp.Pos[0];
        Point3_t ijPos;
        ijPos[0] = itp.Pos[0];
        if(!tcc.geom->IntersectionPoint(iWire, jWire, iPlane, jPlane, cstat, tpc, ijPos[1], ijPos[2])) continue;
        tIDs[jPlane] = jTjPt.id;
        for(unsigned int kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
          auto& kTjPt = slc.mallTraj[kpt];
          // ensure that the planes are different
          if(kTjPt.plane == iTjPt.plane || kTjPt.plane == jTjPt.plane) continue;
          if(kTjPt.xlo > iTjPt.xhi) continue;
          // break out if the x range difference becomes large
          if(kTjPt.xlo > iTjPt.xhi + xcut) break;
          // see if we hit the maxCnt limit
          if(std::find(tMaxed.begin(), tMaxed.end(), kTjPt.id) != tMaxed.end()) continue;
          auto& ktp = slc.tjs[kTjPt.id - 1].Pts[kTjPt.ipt];
          unsigned short kPlane = kTjPt.plane;
          unsigned int kWire = ktp.Pos[0];
          Point3_t ikPos;
          ikPos[0] = ktp.Pos[0];
          if(!tcc.geom->IntersectionPoint(iWire, kWire, iPlane, kPlane, cstat, tpc, ikPos[1], ikPos[2])) continue;
          if(!tcc.geom->IntersectionPoint(iWire, kWire, iPlane, kPlane, cstat, tpc, ikPos[1], ikPos[2])) continue;
          if(std::abs(ijPos[1] - ikPos[1]) > yzcut) continue;
          if(std::abs(ijPos[2] - ikPos[2]) > yzcut) continue;
          // we have a match
          tIDs[kPlane] = kTjPt.id;
          // look for it in the list
          unsigned short indx = 0;
          for(indx = 0; indx < mtIDs.size(); ++indx) if(tIDs == mtIDs[indx]) break;
          if(indx == mtIDs.size()) {
            // not found so add it to mtIDs and add another element to mCnt
            mtIDs.push_back(tIDs);
            mCnt.push_back(0);
          }
          ++mCnt[indx];
          if(mCnt[indx] == maxCnt) {
            // add the Tjs to the list
            tMaxed.insert(tMaxed.end(), tIDs[0]);
            tMaxed.insert(tMaxed.end(), tIDs[1]);
            tMaxed.insert(tMaxed.end(), tIDs[2]);
            hitMaxCnt = true;
            break;
          } // hit maxCnt
        } // kpt
        if(hitMaxCnt) break;
      } // jpt
    } // ipt
    
    if(mCnt.empty()) return;
    
    std::vector<SortEntry> sortVec;
    for(unsigned short indx = 0; indx < mCnt.size(); ++indx) {
      auto& tIDs = mtIDs[indx];
      // count the number of TPs in all Tjs
      float tpCnt = 0;
      for(auto tid : tIDs) {
        auto& tj = slc.tjs[tid - 1];
        tpCnt += NumPtsWithCharge(slc, tj, false);
      } // tid
      float frac = mCnt[indx] / tpCnt;
      frac /= 3;
      // ignore matches with a very low match fraction
      if(frac < 0.05) continue;
      SortEntry se;
      se.index = indx;
      se.val = mCnt[indx];
      sortVec.push_back(se);
    } // ii
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasings);

    matVec.resize(sortVec.size());
    
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      auto& ms = matVec[ii];
      ms.Count = mCnt[indx];
      ms.TjIDs.resize(3);
      for(unsigned short plane = 0; plane < 3; ++plane) ms.TjIDs[plane] = (int)mtIDs[indx][plane];
    } // indx
    
  } // Match3Planes

  /////////////////////////////////////////
  void Match2Planes(TCSlice& slc, std::vector<MatchStruct>& matVec)
  {
    // A simpler faster version of MatchPlanes that only creates two plane matches
    
    std::cout<<"Inside Match2Planes mallTraj size "<<slc.mallTraj.size()<<"\n";
    matVec.clear();
    if(slc.mallTraj.empty()) return;
    
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
    if(tcc.match3DCuts[1] < (float)USHRT_MAX) maxCnt = (unsigned short)tcc.match3DCuts[1];
    // a list of those Tjs
    std::vector<unsigned short> tMaxed;
    
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // see if we hit the maxCnt limit
      if(std::find(tMaxed.begin(), tMaxed.end(), iTjPt.id) != tMaxed.end()) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      unsigned short iPlane = iTjPt.plane;
      unsigned int iWire = itp.Pos[0];
      bool hitMaxCnt = false;
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.plane == iTjPt.plane) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large
        if(jTjPt.xlo > iTjPt.xhi + xcut) break;
        // see if we hit the maxCnt limit
        if(std::find(tMaxed.begin(), tMaxed.end(), jTjPt.id) != tMaxed.end()) continue;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jPlane = jTjPt.plane;
        unsigned int jWire = jtp.Pos[0];
        Point3_t ijPos;
        ijPos[0] = itp.Pos[0];
        if(!tcc.geom->IntersectionPoint(iWire, jWire, iPlane, jPlane, cstat, tpc, ijPos[1], ijPos[2])) continue;
/*
        // require that this be in a dead region in the 3rd plane
        if(slc.nPlanes == 3) {
          unsigned short kPlane = 3 - iPlane - jPlane;
          float fkwire = tcc.geom->WireCoordinate(ijPos[1], ijPos[2], kPlane, tpc, cstat);
          if(fkwire < 0 || fkwire > tcc.maxPos0[kPlane]) continue;
          unsigned int kWire = std::nearbyint(fkwire);
          if(evt.goodWire[kPlane][kWire]) continue;
        } // slc.nPlanes == 3
*/
        tIDs[0] = iTjPt.id;
        tIDs[1] = jTjPt.id;
        // swap the order so that the == operator works correctly
        if(tIDs[0] > tIDs[1]) std::swap(tIDs[0], tIDs[1]);
        // look for it in the list
        unsigned short indx = 0;
        for(indx = 0; indx < mtIDs.size(); ++indx) if(tIDs == mtIDs[indx]) break;
        if(indx == mtIDs.size()) {
          // not found so add it to mtIDs and add another element to mCnt
          mtIDs.push_back(tIDs);
          mCnt.push_back(0);
        }
        ++mCnt[indx];
        if(mCnt[indx] == maxCnt) {
          // add the Tjs to the list
          tMaxed.insert(tMaxed.end(), tIDs[0]);
          tMaxed.insert(tMaxed.end(), tIDs[1]);
          hitMaxCnt = true;
          break;
        } // hit maxCnt
        if(hitMaxCnt) break;
      } // jpt
    } // ipt
    
    if(mCnt.empty()) return;
    
    std::vector<SortEntry> sortVec;
    for(unsigned short indx = 0; indx < mCnt.size(); ++indx) {
      auto& tIDs = mtIDs[indx];
      // count the number of TPs in all Tjs
      float tpCnt = 0;
      for(auto tid : tIDs) {
        auto& tj = slc.tjs[tid - 1];
        tpCnt += NumPtsWithCharge(slc, tj, false);
      } // tid
      float frac = mCnt[indx] / tpCnt;
      frac /= 2;
      // ignore matches with a very low match fraction
      std::cout<<" chk";
      for(auto tid : tIDs) std::cout<<" T"<<tid;
      std::cout<<" frac "<<frac<<"\n";
      if(frac < 0.05) continue;
      SortEntry se;
      se.index = indx;
      se.val = mCnt[indx];
      sortVec.push_back(se);
    } // ii
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
    
    matVec.resize(sortVec.size());
    
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      auto& ms = matVec[ii];
      ms.Count = mCnt[indx];
      ms.TjIDs.resize(2);
      for(unsigned short plane = 0; plane < 2; ++plane) ms.TjIDs[plane] = (int)mtIDs[indx][plane];
    } // indx
    
  } // Match2Planes

  /////////////////////////////////////////
  bool Update(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function only updates SectionFits that need to be re-sorted or re-fit. It returns
    // false if there was a serious error indicating that the pfp should be abandoned
    if(pfp.TP3Ds.empty() || pfp.SectionFits.empty()) return false;
    
    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      if(!sf.NeedsUpdate) continue;
      if(!FitSection(slc, pfp, sfi)) return false;
      if(!SortSection(pfp, sfi)) return false;
      sf.NeedsUpdate = false;
    } // sfi
    
    // ensure that all points (good or not) have a valid SFIndex
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.SFIndex >= pfp.SectionFits.size()) {
//        std::cout<<"Update: P"<<pfp.ID<<" MVI "<<pfp.MVI<<" invalid SFIndex. Fixing it...\n";
        SetSection(slc, pfp, tp3d);
      } // bad SFIndex
    } // tp3d
    
    pfp.Flags[kNeedsUpdate] = false;
    return true;
  } // Update
  
  /////////////////////////////////////////
  bool ReSection(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Re-fit the TP3Ds in sections and add/remove sections to keep ChiDOF of each section close to 1.
    // This function only fails when there is a serious error, otherwise if reasonable fits cannot be
    // achieved, the CanSection flag is set false.
    if(pfp.SectionFits.empty()) return false;
    // This function shouldn't be called if this is the case but it isn't a major failure if it is
    if(!pfp.Flags[kCanSection]) return true;
    // Likewise this shouldn't be attempted if there aren't at least 3 points in 2 planes in 2 sections
    // but it isn't a failure
    if(pfp.TP3Ds.size() < 12) {
      pfp.Flags[kCanSection] = false;
      return true;
    }
    
    prt = (pfp.MVI == debug.MVI);
    
    constexpr float chiLow = 0.5;
    
    // clobber the old sections if more than one exists
    if(pfp.SectionFits.size() > 1) {
      // make one section
      pfp.SectionFits.resize(1);
      // put all of the points in it and fit
      for(auto& tp3d : pfp.TP3Ds) {
        tp3d.SFIndex = 0;
        tp3d.IsGood = true;
      }
      auto& sf = pfp.SectionFits[0];
      if(!FitSection(slc, pfp, 0)) {
//        std::cout<<"ReSection: First fit failed\n";
        return false;
      } // fit failed
      if(sf.ChiDOF < tcc.match3DCuts[5]) return true;
    } // > 1 SectionFit
    // sort by distance from the start
    if(!SortSection(pfp, 0)) return false;
    // require a minimum of 3 points in 2 planes
    unsigned short min2DPts = 3;
    unsigned short fromPt = 0;
    // set the section index to invalid for all points
    for(auto& tp3d : pfp.TP3Ds) tp3d.SFIndex = USHRT_MAX;
    // Guess how many points should be added in each iteration
    unsigned short nPtsToAdd = pfp.TP3Ds.size() / 4;
    // the actual number of points that will be fit in the section
    unsigned short nPts = nPtsToAdd;
    // the minimum number of points
    unsigned short nPtsMin = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1) - fromPt + 1;
    if(nPtsMin >= pfp.TP3Ds.size()) {
      pfp.Flags[kCanSection] = false;
      return true;
    }
    float chiDOF = 0;
    if(nPts < nPtsMin) nPts = nPtsMin;
    // Try to reduce the number of iterations for long pfps
    if(pfp.TP3Ds.size() > 100) {
      unsigned short nhalf = pfp.TP3Ds.size() / 2;
      FitTP3Ds(slc, pfp, fromPt, nhalf, USHRT_MAX, chiDOF);
      if(chiDOF < tcc.match3DCuts[5]) nPts = nhalf;
    }
    bool lastSection = false;
    for(unsigned short sfIndex = 0; sfIndex < 20; ++sfIndex) {
      // Try to add/remove points in each section no more than 20 times
      float chiDOFPrev = 0;
      short nHiChi = 0;
      for(unsigned short nit = 0; nit < 10; ++nit) {
        // Decide how many points to add or subtract after doing the fit
        unsigned short nPtsNext = nPts;
        if(!FitTP3Ds(slc, pfp, fromPt, nPts, USHRT_MAX, chiDOF)) {
//          std::cout<<"RS: MVI "<<pfp.MVI<<" sfi/nit/npts "<<sfIndex<<"/"<<nit<<"/"<<nPts<<" fit failed\n";
          nPtsNext += 1.5 * nPtsToAdd;
        } else if(chiDOF < chiLow) {
          // low chiDOF
          nPtsNext += nPtsToAdd;
          nHiChi = 0;
        } else if(chiDOF > tcc.match3DCuts[5]) {
          // high chiDOF
          ++nHiChi;
          if(nHiChi == 1 && chiDOFPrev > tcc.match3DCuts[5]) {
            // reduce the number of points by 1/2 on the first attempt
            nPtsNext /= 2;
          } else {
            // that didn't work so start subtracting groups of points
            short npnext = (short)nPts - nHiChi * 5;
            // assume this won't work
            nPtsNext = 0;
            if(npnext > nPtsMin) nPtsNext = npnext;
          }
        } else {
          // just right
          nPtsNext = 0;
        }
        // check for passing the end
        if(fromPt + nPtsNext >= pfp.TP3Ds.size()) {
          nPtsNext = pfp.TP3Ds.size() - fromPt;
          lastSection = true;
        }
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<" RS: P"<<pfp.ID<<" sfi/nit/npts "<<sfIndex<<"/"<<nit<<"/"<<nPts;
          myprt<<std::fixed<<std::setprecision(1)<<" chiDOF "<<chiDOF;
          myprt<<" fromPt "<<fromPt;
          myprt<<" nPtsNext "<<nPtsNext<<" lastSection? "<<lastSection;
        }
        if(nPtsNext == 0) break;
        // see if this is the last section
        if(lastSection) break;
        if(chiDOF == chiDOFPrev) {
          if(prt) mf::LogVerbatim("TC")<<" MVI "<<pfp.MVI<<" chiDOF not changing\n";
          break;
        }
        nPts = nPtsNext;
        chiDOFPrev = chiDOF;
      } // nit
      // finished this section. Assign the points to it
      unsigned short toPt = fromPt + nPts;
      if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
      for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) pfp.TP3Ds[ipt].SFIndex = sfIndex;
      // See if there are enough points remaining to reconstruct another section if this isn't known
      // to be the last section
      if(!lastSection) {
        // this will be the first point in the next section
        unsigned short nextFromPt = fromPt + nPts;
        // See if it will have enough points to be reconstructed
        unsigned short nextToPtMin = Find3DRecoRange(slc, pfp, nextFromPt, min2DPts, 1);
//        std::cout<<" sfi "<<sfIndex<<" nextFromPt "<<nextFromPt<<" nextToPtMin "<<nextToPtMin<<"\n";
        if(nextToPtMin == USHRT_MAX) {
          // not enough points so this is the last section
          lastSection = true;
          // assign the remaining points to the last section
          for(unsigned short ipt = nextFromPt; ipt < pfp.TP3Ds.size(); ++ipt) pfp.TP3Ds[ipt].SFIndex = sfIndex;
        } 
      } // !lastSection
      // Do a final fit and update the points. Don't worry about a poor ChiDOF
      FitSection(slc, pfp, sfIndex);
      if(!SortSection(pfp, 0)) {
//        std::cout<<"RS: SortSection failed\n";
        return false;
      }
      if(lastSection) break;
      // Prepare for the next section.
      fromPt = fromPt + nPts;
      nPts = nPtsToAdd;
      nPtsMin = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1) - fromPt + 1;
      if(nPtsMin >= pfp.TP3Ds.size()) break;
      // add a new section
      pfp.SectionFits.resize(pfp.SectionFits.size() + 1);
    } // snit
    
    // see if the last sf is valid
    if(pfp.SectionFits.size() > 1 && pfp.SectionFits.back().ChiDOF < 0) {
      unsigned short badSFI = pfp.SectionFits.size() - 1;
      // remove it
      pfp.SectionFits.pop_back();
      for(unsigned short ipt = pfp.TP3Ds.size() - 1; ipt > 0; --ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(tp3d.SFIndex < badSFI) break;
        --tp3d.SFIndex;
      }
      pfp.SectionFits.back().NeedsUpdate = true;
    } // bad last SF
    
    // Ensure that the points at the end are in the last section
    for(unsigned short ipt = pfp.TP3Ds.size() - 1; ipt > 0; --ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.SFIndex < pfp.SectionFits.size()) break;
      tp3d.SFIndex = pfp.SectionFits.size() - 1;
      pfp.Flags[kNeedsUpdate] = true;
      pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
    } // tp3d
    
    Update(slc, pfp, prt);
    
    // set CanSection false if the chisq is poor in any section
    for(auto& sf : pfp.SectionFits) {
      if(sf.ChiDOF > tcc.match3DCuts[5]) pfp.Flags[kCanSection] = false;
    }
    
    return true;
  } // resection
  
  /////////////////////////////////////////
  void CountBadPoints(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, unsigned short& nBadPts, unsigned short& firstBadPt)
  {
    // Count the number of points whose pull exceeds tcc.match3DCuts[4]
    firstBadPt = USHRT_MAX;
    nBadPts = 0;
    if(fromPt > pfp.TP3Ds.size() - 1) {
      nBadPts = USHRT_MAX;
      return;
    }
    if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
    bool first = true;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      // don't clobber a point if it is on a TP that is overlapping another Tj. This will
      // happen for points close to a vertex and when trajectories cross
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvOverlap]) continue;
      if(PointPull(pfp, tp3d) < tcc.match3DCuts[4]) continue;
      ++nBadPts;
      if(first) {
        first = false;
        firstBadPt = ipt;
      }
    } // ipt
  } // CountBadPoints
  
  /////////////////////////////////////////
  void KillBadPoints(TCSlice& slc, PFPStruct& pfp, float pullCut, bool prt)
  {
    // Find bad TP3Ds points and remove them in all sections
    unsigned short nbad = 0;
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.IsBad) {
        ++nbad;
        pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
        continue;
      } // IsBad
      if(PointPull(pfp, tp3d) > pullCut) {
        tp3d.IsBad = true;
        pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
        ++nbad;
      }
    } // tp3d
    if(nbad == 0) return;
    pfp.Flags[kNeedsUpdate] = true;
    std::vector<TP3D> ntp3ds;
    for(auto& tp3d : pfp.TP3Ds) if(!tp3d.IsBad) ntp3ds.push_back(tp3d);
    pfp.TP3Ds = ntp3ds;
    Update(slc, pfp, prt);

  } // KillBadPoints
  
  /////////////////////////////////////////
  bool CanSection(TCSlice& slc, PFPStruct& pfp)
  {
    // analyze the TP3D vector to determine if it can be reconstructed in 3D in more than one section with
    // the requirement that there are at least 3 points in two planes
    if(pfp.AlgMod[kJunk3D]) return false;
    if(pfp.TP3Ds.size() < 12) return false;
    unsigned short toPt = Find3DRecoRange(slc, pfp, 0, 3, 1);
    if(toPt > pfp.TP3Ds.size()) return false;
    unsigned short nextToPt = Find3DRecoRange(slc, pfp, toPt, 3, 1);
    if(nextToPt > pfp.TP3Ds.size()) return false;
    return true;
  } // CanSection

  /////////////////////////////////////////
  unsigned short Find3DRecoRange(TCSlice& slc, const PFPStruct& pfp, unsigned short fromPt, unsigned short min2DPts, short dir)
  {
    // Scans the TP3Ds vector starting at fromPt until it finds min2DPts in two planes. It returns
    // with the index of that point (+1) in the TP3Ds vector. The dir variable defines the scan direction in
    // the TP3Ds vector
    if(fromPt > pfp.TP3Ds.size() - 1) return USHRT_MAX;
    if(pfp.TP3Ds.size() < 2 * min2DPts) return USHRT_MAX;
    if(dir == 0) return USHRT_MAX;
    
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      unsigned short ipt = fromPt + ii;
      if(dir < 0) ipt = fromPt - ii;
      if(ipt >= pfp.TP3Ds.size()) break;
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      unsigned short plane = DecodeCTP(slc.tjs[tp3d.TjID - 1].CTP).Plane;
      ++cntInPln[plane];
      unsigned short enufInPlane = 0;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] >= min2DPts) ++enufInPlane;
      if(enufInPlane > 1) return ipt + 1;
      if(dir < 0 && ipt == 0) break;
    } // ipt
    return USHRT_MAX;
  } // Find3DRecoRange
  
  /////////////////////////////////////////
  void GetRange(PFPStruct& pfp, unsigned short sfIndex, unsigned short& fromPt, unsigned short& npts)
  {
    fromPt = USHRT_MAX;
    if(sfIndex >= pfp.SectionFits.size()) return;
    if(pfp.TP3Ds.empty()) return;
    fromPt = USHRT_MAX;
    npts = 0;
    // Note that no test is made for not-good TP3Ds here since that would give a wrong npts count
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.SFIndex < sfIndex) continue;
      if(tp3d.SFIndex > sfIndex) break;
      if(fromPt == USHRT_MAX) fromPt = ipt;
      ++npts;
    } // ipt
  } // GetRange
  
  /////////////////////////////////////////
  bool FitSection(TCSlice& slc, PFPStruct& pfp, unsigned short sfIndex)
  {
    // Fits the TP3D points in the selected section to a 3D line with the origin at the center of
    // the section
    if(pfp.TP3Ds.size() < 4) return false;
    if(sfIndex >= pfp.SectionFits.size()) return false;
//    if(pfp.Flags[kJunk3D]) return true;

    unsigned short fromPt = USHRT_MAX;
    unsigned short npts = 0;
    GetRange(pfp, sfIndex, fromPt, npts);
    if(fromPt == USHRT_MAX) return false;
    if(npts < 4) return false;
    
    // check for errors
    for(unsigned short ipt = fromPt; ipt < fromPt + npts; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.SFIndex != sfIndex) {
//        std::cout<<"FitSection: MVI "<<pfp.MVI<<" sfIndex "<<sfIndex<<" points aren't contiguous\n";
        return false;
      }
    } // ipt
    
    // fit these points and update
    float chiDOF = 999;
    return FitTP3Ds(slc, pfp, fromPt, npts, sfIndex, chiDOF);

  } // FitSection

/////////////////////////////////////////
  SectionFit FitTP3Ds(TCSlice& slc, const std::vector<TP3D>& tp3ds, unsigned short fromPt, short fitDir, unsigned short nPtsFit)
  {
    // fits the points and returns the fit results in a SectionFit struct. This function assumes that the
    // vector of TP3Ds exists in the slc.TPCID
    
    SectionFit sf;
    sf.ChiDOF = 999;
    if(nPtsFit < 5) return sf;
    if(!(fitDir == -1 || fitDir == 1)) return sf;
    if(fitDir ==  1 && fromPt + nPtsFit >tp3ds.size()) return sf;
    if(fitDir == -1 && fromPt < 3) return sf;
    
    // put the offset, cosine-like and sine-like components in a vector
    std::vector<std::array<double, 3>> ocs(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      auto planeID = geo::PlaneID(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      // plane offset
      ocs[plane][0] = tcc.geom->WireCoordinate(0, 0, planeID);
      // get the "cosine-like" component
      ocs[plane][1] = tcc.geom->WireCoordinate(1, 0, planeID) - ocs[plane][0];
      // the "sine-like" component
      ocs[plane][2] = tcc.geom->WireCoordinate(0, 1, planeID) - ocs[plane][0];
    } // plane

    const unsigned int nvars = 4;
    unsigned int npts = 0;

    // count the number of TPs in each plane
    std::vector<unsigned short> cntInPln(slc.nPlanes, 0);
    // and define the X position for the fit origin
    double x0 = 0.;
    for(short ii = 0; ii < nPtsFit; ++ii) {
      short ipt = fromPt + fitDir * ii;
      if(ipt < 0 || ipt >= tp3ds.size()) break;
      auto& tp3d = tp3ds[ipt];
      if(!tp3d.IsGood) continue;
      if(tp3d.TPXErr2 < 0.0001) return sf;
      x0 += tp3d.TPX;
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      ++cntInPln[plane];
      ++npts;
    } // ipt
    if(npts < 6) return sf;
    // ensure there are at least three points in at least two planes
    unsigned short enufInPlane = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] > 2) ++enufInPlane;
    if(enufInPlane < 2) return sf;

    x0 /= (double)npts;
        
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
        
    unsigned short cnt = 0;
    double weight = 1;
    for(short ii = 0; ii < nPtsFit; ++ii) {
      short ipt = fromPt + fitDir * ii;
      auto& tp3d = tp3ds[ipt];
      if(!tp3d.IsGood) continue;
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double x = tp3d.TPX - x0;
      A[cnt][0] = weight * ocs[plane][1];
      A[cnt][1] = weight * ocs[plane][2];
      A[cnt][2] = weight * ocs[plane][1] * x;
      A[cnt][3] = weight * ocs[plane][2] * x;
      w[cnt] = weight * (tp3d.Wire - ocs[plane][0]);
      ++cnt;
    } // ipt
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    if(!ok) return sf;
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    
    // TODO: The direction is reversed for some reason
    norm *= -1;

    sf.Dir[0] = 1 / norm;
    sf.Dir[1] = tVec[2] / norm;
    sf.Dir[2] = tVec[3] / norm;
    sf.Pos[0] = x0;
    sf.Pos[1] = tVec[0];
    sf.Pos[2] = tVec[1];
    sf.NPts = npts;
    
    // Calculate errors from sigma * (A^T * A)^(-1) where sigma is the
    // error on the wire number (= 1)
    // TODO: is this an expensive calculation?
    TMatrixD AT(nvars, npts); 
    AT.Transpose(A);
    TMatrixD ATA = AT * A;
    double *det = 0;
    ATA.Invert(det);
    sf.DirErr[1] = -sqrt(ATA[2][2]) / norm;
    sf.DirErr[2] = -sqrt(ATA[3][3]) / norm;
    
    // calculate ChiDOF
    sf.ChiDOF = 0;
    // project this 3D vector into a TP in every plane
    std::vector<TrajPoint> plnTP(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      plnTP[plane] = MakeBareTP(slc, sf.Pos, sf.Dir, inCTP);
    } // plane
    // a local position
    Point3_t pos;
    sf.DirErr[0] = 0.;
    for(short ii = 0; ii < nPtsFit; ++ii) {
      short ipt = fromPt + fitDir * ii;
      auto& tp3d = tp3ds[ipt];
      if(!tp3d.IsGood) continue;
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double dw = tp3d.Wire - plnTP[plane].Pos[0];
      // dt/dW was stored in DeltaRMS by MakeBareTP
      double t = dw * plnTP[plane].DeltaRMS;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
      // Note that the tp3d position is directly above the wire position and not the
      // point at the distance of closest approach. Delta is the difference in the
      // drift direction in cm
      double delta = pos[0] - tp3d.TPX;
      sf.ChiDOF += delta * delta / tp3d.TPXErr2;
      // estimate the X slope error ~ X direction vector with an overly simple average
      double dangErr = delta / dw;
      sf.DirErr[0] += dangErr * dangErr;
    } // indx
    sf.DirErr[0] = sqrt(sf.DirErr[0]) / (double)nPtsFit;
    sf.ChiDOF /= (float)(npts - 4);
    return sf;

  } // FitTP3Ds

  /////////////////////////////////////////
  bool FitTP3Ds(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short nPtsFit, unsigned short sfIndex, float& chiDOF)
  {
    // Fit points in the pfp.TP3Ds vector fromPt. This function returns chiDOF but
    // doesn't update the TP3Ds unless sfIndex refers to a valid SectionFit in the pfp. 
    // No check is made to ensure that the TP3D SFIndex variable is compatible with sfIndex
    
    chiDOF = 999;
    if(nPtsFit < 5) return false;
    if(fromPt + nPtsFit > pfp.TP3Ds.size()) return false;
    
    auto sf = FitTP3Ds(slc, pfp.TP3Ds, fromPt, 1, nPtsFit);
    chiDOF = sf.ChiDOF;
    if(chiDOF > 900) return false;
    
    // don't update the pfp?
    if(sfIndex >= pfp.SectionFits.size()) return true;

    // update the pfp Sectionfit
    pfp.SectionFits[sfIndex] = sf;
    // update the TP3Ds
    // project this 3D vector into a TP in every plane
    std::vector<TrajPoint> plnTP(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      plnTP[plane] = MakeBareTP(slc, sf.Pos, sf.Dir, inCTP);
    } // plane
    Point3_t pos;
    bool needsSort = false;
    double prevAlong;
    for(unsigned short ipt = fromPt; ipt < fromPt + nPtsFit; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double dw = tp3d.Wire - plnTP[plane].Pos[0];
      // dt/dW was stored in DeltaRMS by MakeBareTP
      double t = dw * plnTP[plane].DeltaRMS;
      if(ipt == fromPt) {
        prevAlong = t;
      } else {
        if(t < prevAlong) needsSort = true;
        prevAlong = t;
      }
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
      // Note that the tp3d position is directly above the wire position and not the
      // distance of closest approach. The Delta variable is the difference in the
      // drift direction in cm
      double delta = pos[0] - tp3d.TPX;
      tp3d.Pos = pos;
      tp3d.Dir = sf.Dir;
      tp3d.along = t;
      if(tp3d.IsGood) sf.ChiDOF += delta * delta / tp3d.TPXErr2;
    } // ipt
    if(needsSort) {
//      std::cout<<"FitTP3Ds: P"<<pfp.ID<<" sf "<<sfIndex<<" needs to be sorted\n";
      SortSection(pfp, sfIndex);
    }
    pfp.SectionFits[sfIndex].NeedsUpdate = false;
    return true;

  } // FitTP3Ds

  /////////////////////////////////////////
  void KinkFit(TCSlice& slc, const PFPStruct& pfp, unsigned short atPt, double fitLen, double& dang, double& dangSig)
  {
    // calculates a kink angle at the point atPt in the pfp TP3Ds vector by doing a fit
    // of points within fitLen distance on both sides of that point
        
    dang = -1;
    dangSig = 0;
    if(atPt < 4) return;
    if(atPt > pfp.TP3Ds.size() - 5) return;
    
    auto& atPos = pfp.TP3Ds[atPt].Pos;
    double fLen2 = fitLen * fitLen;
    unsigned short nPtsMinus = 0;
    bool gotenuf = false;
    for(unsigned short ipt = atPt - 1; ipt > 0; --ipt) {
      ++nPtsMinus;
      if(PosSep2(pfp.TP3Ds[ipt].Pos, atPos) > fLen2) {
        gotenuf = true;
        break;
      }
    } // ipt
    if(!gotenuf) return;
    unsigned short nPtsPlus = 0;
    gotenuf = false;
    for(unsigned short ipt = atPt + 1; ipt < pfp.TP3Ds.size(); ++ipt) {
      ++nPtsPlus;
      if(PosSep2(pfp.TP3Ds[ipt].Pos, atPos) > fLen2) {
        gotenuf = true;
        break;
      }
    } // ipt
    if(!gotenuf) return;
    
    auto fitMinus = FitTP3Ds(slc, pfp.TP3Ds, atPt - 1, -1, nPtsMinus);
    if(fitMinus.ChiDOF > 10) return;
    auto fitPlus = FitTP3Ds(slc, pfp.TP3Ds, atPt + 1, 1, nPtsMinus);
    if(fitPlus.ChiDOF > 10) return;
    dang = DeltaAngle(fitMinus.Dir, fitPlus.Dir);
    double cnt = 0;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) {
      double maxErr = fitMinus.DirErr[xyz];
      if(fitPlus.DirErr[xyz] > maxErr) maxErr = fitPlus.DirErr[xyz];
      if(maxErr <= 0) continue;
      double sig = std::abs(fitMinus.Dir[xyz] - fitPlus.Dir[xyz]) / maxErr;
      dangSig += sig * sig;
      ++cnt;
    } // xyz
    dangSig = sqrt(dangSig) / cnt;
/* Debug output to port into Excel
    mf::LogVerbatim myprt("TC");
    myprt<<"KF, "<<atPt;
    myprt<<","<<PrintPos(slc, slc.tjs[pfp.TP3Ds[atPt].TjID - 1].Pts[pfp.TP3Ds[atPt].TPIndex]);
    myprt<<std::fixed<<std::setprecision(3);
    for(unsigned short xyz = 0; xyz < 3; ++xyz) myprt<<","<<fitMinus.Dir[xyz]<<","<<fitMinus.DirErr[xyz];
    myprt<<", "<<dang;
    myprt<<", "<<dangSig;
*/
  } // KinkFit

  /////////////////////////////////////////
  void SplitAtKinks(TCSlice& slc, std::vector<PFPStruct>& pfpVec, bool prt)
  {
    // Looks for a kink in pfpVec[0]. If one is found, the pfp is split and the
    // second one is put in pfpVec[1]
    
    if(tcc.useAlg[kKink3D]) return;
    if(tcc.kinkCuts.size() < 5) return;
    if(pfpVec.size() != 1) return;
    if(pfpVec[0].AlgMod[kJunk3D]) return;

    // Cut on a FOM = (tp3d 3D kink angle) * (tp 2D kink significance)
    double big = tcc.kinkCuts[0] * tcc.kinkCuts[3];
    unsigned short atPt = USHRT_MAX;
    double dang, dangSig;
    auto& p1 = pfpVec[0];
    if(prt) mf::LogVerbatim("TC")<<"SAK: P"<<p1.ID<<" min (3D angle) * (2D significance) cut "<<big<<" fit Len "<<tcc.kinkCuts[4];
    for(unsigned short ipt = 0; ipt < p1.TP3Ds.size(); ++ipt) {
      KinkFit(slc, p1, ipt, tcc.kinkCuts[4], dang, dangSig);
      if(dang < 0) continue;
      // multiply the kink angle by the 2D kink significance
      auto& tp3d = p1.TP3Ds[ipt];
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      // TODO: This should be replaced by a 3D kink significance after the SVD fit errors
      // are done correctly
      double fom = dang * tp.KinkSig;
      if(fom < big) continue;
      big = fom;
      atPt = ipt;
    } // ipt
    if(atPt == USHRT_MAX) return;
    auto& tp3d = p1.TP3Ds[atPt];
    auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
    if(prt) mf::LogVerbatim("TC")<<"SAK:  found kink big "<<big<<" TP3Ds["<<atPt<<"] 2D pos "<<PrintPos(slc, tp);
    auto p2 = CreatePFP(slc);
    p2.MVI = p1.MVI;
    p2.ID = pfpVec[0].ID + 1;
    if(!Split(slc, p1, atPt, p2, prt)) return;
    pfpVec.push_back(p2);

  } // SplitAtKinks

  /////////////////////////////////////////
  bool Split(TCSlice& slc, PFPStruct& p1, unsigned short atPt, PFPStruct& p2, bool prt)
  {
    // Splits pfp p1 at the specified point atPt of the TP3Ds vector.
    // Points < atPt are kept in p1. Points >= atPt are moved to p2. This
    // function returns true if the split is successful
    if(p1.TP3Ds.size() < 4) return false;
    if(atPt > p1.TP3Ds.size() - 3) return false;
    // ensure that p2 was created properly
    if(p2.SectionFits.size() != 1 || !p2.TP3Ds.empty()) return false;    
    // ensure that p1 is updated
    if(p1.Flags[kNeedsUpdate] && !Update(slc, p1, prt)) return false;
    // determine if the TP -> P assn has been made
    bool doTPAssns = false;
    auto& tp3d0 = p1.TP3Ds[0];
    if(tp3d0.TjID > 0 && slc.tjs[tp3d0.TjID - 1].Pts[tp3d0.TPIndex].InPFP > 0) doTPAssns = true;
    for(unsigned short ipt = atPt; ipt < p1.TP3Ds.size(); ++ipt) p2.TP3Ds.push_back(p1.TP3Ds[ipt]);
    p2.TjIDs.clear();
    for(auto& tp3d : p2.TP3Ds) {
      tp3d.SFIndex = 0;
      if(std::find(p2.TjIDs.begin(), p2.TjIDs.end(), tp3d.TjID) == p2.TjIDs.end()) p2.TjIDs.push_back(tp3d.TjID);
    } // tp3d
    // ReSection and Update if necessary
    p1.TP3Ds.resize(atPt);
    p1.TjIDs.clear();
    for(auto& tp3d : p1.TP3Ds) {
      if(std::find(p1.TjIDs.begin(), p1.TjIDs.end(), tp3d.TjID) == p1.TjIDs.end()) p1.TjIDs.push_back(tp3d.TjID);
    }
    if(!ReSection(slc, p1, prt)) return false;
    if(!ReSection(slc, p2, prt)) return false;
    // Correct the TP -> P2 assn 
    if(doTPAssns) {
      for(auto& tp3d : p2.TP3Ds) {
        if(tp3d.TjID <= 0) continue;
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if(tp.InPFP != p1.ID) {
          std::cout<<"Split: oops P"<<p2.ID<<" TP "<<PrintPos(slc, tp)<<" InPFP = "<<tp.InPFP<<" != "<<p1.ID<<"\n";
        }
        tp.InPFP = p2.ID;
      } // tp3d
    }
    mf::LogVerbatim("TC")<<"Split: P"<<p1.ID<<" -> P"<<p2.ID;
    return true;
    
  } // Split

  /////////////////////////////////////////
  void ReconcileVertices(TCSlice& slc, PFPStruct& pfp, bool prt)
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
    if(pfp.Flags[kJunk3D]) return;
    
    // first make a list of all Tjs 
    std::vector<int> tjList;
    for(auto& tp3d : pfp.TP3Ds) {
      if(!tp3d.IsGood) continue;
      // ignore single hits
      if(tp3d.TjID <= 0) continue;
      if(std::find(tjList.begin(), tjList.end(), tp3d.TjID) == tjList.end()) tjList.push_back(tp3d.TjID);
    } // tp3d
    // look for 3D vertices associated with these Tjs and list of
    // orphan 2D vertices - those that are not matched to 3D vertices
    std::vector<int> vx2List, vx3List;
    for(auto tid : tjList) {
      auto& tj = slc.tjs[tid - 1];
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] <= 0) continue;
        auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
        if(vx2.Vx3ID > 0) {
          if(std::find(vx3List.begin(), vx3List.end(), vx2.Vx3ID) == vx3List.end()) vx3List.push_back(vx2.Vx3ID);
          // 3D vertex exists
        } else {
          // no 3D vertex 
          if(std::find(vx2List.begin(), vx2List.end(), tj.VtxID[end]) == vx2List.end()) vx2List.push_back(tj.VtxID[end]);
        } // no 3D vertex
      } // end
    } // tid
    // no vertex reconciliation is necessary
    if(vx2List.empty() && vx3List.empty()) return;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"RV: P"<<pfp.ID<<" ->";
      for(auto tid : tjList) myprt<<" T"<<tid;
      myprt<<" ->";
      for(auto vid : vx3List) myprt<<" 3V"<<vid;
      if(!vx2List.empty()) {
        myprt<<" orphan";
        for(auto vid : vx2List) myprt<<" 2V"<<vid;
      }
    } // prt
    // Just kill the orphan 2D vertices regardless of their score. 
    // This is an indicator that the vertex was created between two tjs 
    // that maybe should have been reconstructed as one or alternatively 
    // as two Tjs. This decision presumes the existence of a 3D kink 
    // algorithm that doesn't yet exist...
    for(auto vid : vx2List) {
      auto& vx2 = slc.vtxs[vid - 1];
      MakeVertexObsolete("RV", slc, vx2, true);
    } // vx2List
    // ignore the T -> 2V -> 3V assns (if any exist) and try to directly
    // attach to 3D vertices at both ends
    AttachToAnyVertex(slc, pfp, tcc.vtx3DCuts[2], prt);
    // check for differences and while we are here, see if the pfp was attached
    // to a neutrino vertex and the direction is wrong
    int neutrinoVx = 0;
    if(!slc.pfps.empty()) {
      auto& npfp = slc.pfps[0];
      bool neutrinoPFP = (npfp.PDGCode == 12 || npfp.PDGCode == 14);
      if(neutrinoPFP) neutrinoVx = npfp.Vx3ID[0];
    } // pfps exist
    unsigned short neutrinoVxEnd = 2;
    for(unsigned short end = 0; end < 2; ++end) {
      // see if a vertex got attached
      if(pfp.Vx3ID[end] <= 0) continue;
      if(pfp.Vx3ID[end] == neutrinoVx) neutrinoVxEnd = end;
      // see if this is a vertex in the list using the T -> 2V -> 3V assns
      if(std::find(vx3List.begin(), vx3List.end(), pfp.Vx3ID[end]) != vx3List.end()) continue;
//      std::cout<<"RV: P"<<pfp.ID<<" was attached to 3V"<<pfp.Vx3ID[end]<<" but a P -> T -> 2V -> 3V assn exists. Write some code to clobber this assn or deal with it somehow.\n";
    } // end
    if(neutrinoVxEnd < 2 && neutrinoVxEnd != 0) Reverse(slc, pfp);
    
    return;
  } // ReconcileVertices

  /////////////////////////////////////////
  void TrimEndPts(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Check for a wire signal and trim points starting
    // at the ends until there at least 2 planes that have a wire signal. The
    // sections that have points removed are re-fit without those points and
    // another check is made in a second iteration

    if(!tcc.useAlg[kTEP3D]) return;
    if(pfp.ID <= 0) return;
    // Trimming short tracks that are barely reconstructable isn't a good idea
    if(pfp.TP3Ds.size() < 10) return;
    // don't trim shower-like pfps
//    if(IsShowerLike(slc, pfp.TjIDs)) return;
    if(pfp.Flags[kJunk3D]) return;
    // don't trim if the pfp failed in ReSection
//    if(!pfp.Flags[kCanSection]) return;
    
    auto pWork = pfp;
    
    pWork.Flags[kNeedsUpdate] = false;
    for(unsigned short nit = 0; nit < 2; ++nit) {
      // create a vector of valid points
      std::vector<bool> validPt(pWork.TP3Ds.size(), true);
      // inspect end 0
      for(unsigned short ipt = 0; ipt < pWork.TP3Ds.size(); ++ipt) {
        auto& tp3d = pWork.TP3Ds[ipt];
        if(!tp3d.IsGood) continue;
        unsigned short cnt = 0;
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pWork.TPCID.Cryostat, pWork.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, tp3d.Pos, inCTP);
          if(SignalAtTp(tp)) ++cnt;
        } // plane
        if(cnt >= 2) break;
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<"TEP: P"<<pWork.ID<<" nit "<<nit<<" trim TP3D "<<tp3d.TjID<<"_"<<tp3d.TPIndex;
          unsigned int mcp = FindMCPIndex(slc, tp3d);
          if(mcp != UINT_MAX) myprt<<" mcp "<<mcp;
        } // prt
        validPt[ipt] = false;
        pWork.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
        pWork.Flags[kNeedsUpdate] = true;
      } // ipt
      // inspect the other end
      for(unsigned short ipt = pWork.TP3Ds.size() - 1; ipt > 0; --ipt) {
        auto& tp3d = pWork.TP3Ds[ipt];
        if(!tp3d.IsGood) continue;
        unsigned short cnt = 0;
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pWork.TPCID.Cryostat, pWork.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, tp3d.Pos, inCTP);
          if(SignalAtTp(tp)) ++cnt;
        } // plane
        if(cnt >= 2) break;
        if(prt) mf::LogVerbatim("TC")<<"TEP: P"<<pWork.ID<<" nit "<<nit<<" trim TP3D "<<tp3d.TjID<<"_"<<tp3d.TPIndex;
        validPt[ipt] = false;
        pWork.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
        pWork.Flags[kNeedsUpdate] = true;
      } // ipt
      if(pWork.Flags[kNeedsUpdate]) {
        // trim the points
        // find the first good point
        unsigned short firstGood = 0;
        for(firstGood = 0; firstGood < pWork.TP3Ds.size(); ++firstGood) if(validPt[firstGood]) break;
        if(firstGood == pWork.TP3Ds.size()) break;
        // and the last good point
        unsigned short lastGood = pWork.TP3Ds.size();
        for(lastGood = pWork.TP3Ds.size() - 1; lastGood > 0; --lastGood) if(validPt[lastGood]) break;
        ++lastGood;
        if(firstGood == 0 && lastGood == pWork.TP3Ds.size()) break;
        std::vector<TP3D> temp(pWork.TP3Ds.begin() + firstGood, pWork.TP3Ds.begin() + lastGood);
        pWork.TP3Ds = temp;
        if(!Update(slc, pWork, prt)) {
//          std::cout<<"Update P"<<pWork.ID<<" failed in TrimEndPts. Recovering...\n";
          return;
        }
      }
    } // nit

    if(pWork.Flags[kNeedsUpdate]) Update(slc, pWork, prt);

    pfp = pWork;
    
  } // TrimEndPts

  /////////////////////////////////////////
  void FillGaps3D(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Look for gaps in each plane in the TP3Ds vector in planes in which 
    // the projection of the pfp angle is large (~> 60 degrees). Hits
    // reconstructed at large angles are poorly reconstructed which results
    // in poorly reconstructed 2D trajectories
    
    if(pfp.ID <= 0) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.SectionFits.empty()) return;
    if(!tcc.useAlg[kFillGaps3D]) return;
    if(pfp.Flags[kJunk3D]) return;
    
    // Only print APIR details if MVI is set
    bool foundMVI = (tcc.dbgPFP && pfp.MVI == debug.MVI);
    
    // make a copy in case something goes wrong
    auto pWork = pfp;
    
    if(tcc.useAlg[kTCWork2]) {
      unsigned short nPtsAdded = 0;
      unsigned short fromPt = 0;
      unsigned short toPt = pWork.TP3Ds.size();
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(pWork.TPCID.Cryostat, pWork.TPCID.TPC, plane);
        unsigned short nWires, nAdd;
        AddPointsInRange(slc, pWork, fromPt, toPt, inCTP, tcc.match3DCuts[4], nWires, nAdd, foundMVI);
        if(pWork.Flags[kNeedsUpdate]) Update(slc, pWork, prt);
        nPtsAdded += nAdd;
      } // plane
/*
      // Start looking for this pattern near the beginning of the TP3Ds vector
      unsigned short atPt = 0;
      for(atPt = 1; atPt < pWork.TP3Ds.size(); ++atPt) {
        if(PosSep2(pWork.TP3Ds[atPt].Pos, pWork.TP3Ds[0].Pos) > span2) break;
      } // ipt
      if(foundMVI) std::cout<<"FG3D P"<<pWork.ID<<" atPt "<<atPt<<"\n";
      if(atPt >= pWork.TP3Ds.size() - 1) return;
      // count of TPs in each plane within +/- span
      std::vector<unsigned short> plnCnt(slc.nPlanes);
      for(unsigned short ipt = atPt; ipt < pWork.TP3Ds.size(); ++ipt) {
        std::fill(plnCnt.begin(), plnCnt.end(), 0);
        for(unsigned short jpt = ipt; jpt > 0; --jpt) {
          unsigned short plane = DecodeCTP(pWork.TP3Ds[jpt].CTP).Plane;
          ++plnCnt[plane];
          if(PosSep2(pWork.TP3Ds[jpt].Pos, pWork.TP3Ds[ipt].Pos) > span2) break;
        } // jpt
        for(unsigned short jpt = ipt; jpt < pWork.TP3Ds.size(); ++jpt) {
          unsigned short plane = DecodeCTP(pWork.TP3Ds[jpt].CTP).Plane;
          ++plnCnt[plane];
          if(PosSep2(pWork.TP3Ds[jpt].Pos, pWork.TP3Ds[ipt].Pos) > span2) break;
        } // jpt
        // look for plnCnt == 0
        for(unsigned short plane = 0; plane < plnCnt.size(); ++plane) {
          if(plnCnt[plane] > 0) continue;
          // Look for a TP near this position
          CTP_t inCTP = EncodeCTP(pWork.TPCID.Cryostat, pWork.TPCID.TPC, plane);
          unsigned short nWires, nAdd;
          AddPointsInRange(slc, pWork, ipt, ipt, inCTP, tcc.match3DCuts[4], nWires, nAdd, foundMVI);
          nPtsAdded += nAdd;
        } // plane
      } // ipt
*/
      if(prt) mf::LogVerbatim("TC")<<"FG3D P"<<pWork.ID<<" added "<<nPtsAdded<<" points";
      if(pWork.Flags[kNeedsUpdate] && !Update(slc, pWork, prt)) {
        mf::LogVerbatim("TC")<<"  Update failed after adding points. Restored P"<<pfp.ID;
        return;
      }
      pfp = pWork;
      return;
    } // TCWork2
    
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pWork.TPCID.Cryostat, pWork.TPCID.TPC, plane);
      // check the start
      unsigned short fromPt = 0;
      unsigned short toPt;
      unsigned short nWires = 0, nAdd = 0;
      for(toPt = 0; toPt < pWork.TP3Ds.size(); ++toPt) if(pWork.TP3Ds[toPt].CTP == inCTP) break;
      if(toPt > 5) {
        AddPointsInRange(slc, pWork, fromPt, toPt, inCTP, tcc.match3DCuts[4], nWires, nAdd, foundMVI);
        if(prt) {
          mf::LogVerbatim("TC")<<"FillGaps3D: P"<<pWork.ID<<" Search for gaps in plane "<<plane<<" nWires "<<nWires<<" nAdd "<<nAdd;
        }
      }
      // now check the end
      toPt = pWork.TP3Ds.size() - 1;
      for(fromPt = toPt; fromPt > 0; --fromPt) if(pWork.TP3Ds[fromPt].CTP == inCTP) break;
      if(fromPt < toPt - 5) {
        if(prt) {
          mf::LogVerbatim("TC")<<"FillGaps3D: P"<<pWork.ID<<" Search for gaps in plane "<<plane<<" nWires "<<nWires<<" nAdd "<<nAdd;
        }
        AddPointsInRange(slc, pWork, fromPt, toPt, inCTP, tcc.match3DCuts[4], nWires, nAdd, foundMVI);
      }
    } // plane
        
    if(pWork.Flags[kNeedsUpdate] && !Update(slc, pWork, prt)) {
      mf::LogVerbatim("TC")<<"  Update failed after adding points. Restored P"<<pfp.ID;
      return;
    }
    pfp = pWork;
    
  } // FillGaps3D
  
  /////////////////////////////////////////
  bool ValidTwoPlaneMatch(TCSlice& slc, PFPStruct& pfp)
  {
    // This function checks the third plane in the PFP when only two Tjs are 3D-matched to
    // ensure that the reason for the lack of a 3rd plane match is that it is in a dead region.
    // This function should be used after an initial fit is done and the TP3Ds are sorted
    if(pfp.TjIDs.size() != 2) return false;
    if(slc.nPlanes != 3) return false;
    if(pfp.TP3Ds.empty()) return false;
    
    // find the third plane
    std::vector<unsigned short> planes;
    for(auto tid : pfp.TjIDs) planes.push_back(DecodeCTP(slc.tjs[tid - 1].CTP).Plane);
    unsigned short thirdPlane = 3 - planes[0] - planes[1];
    CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, thirdPlane);
    // Project the 3D position at the start into the third plane
    auto tp = MakeBareTP(slc, pfp.TP3Ds[0].Pos, inCTP);
    unsigned int wire0 = 0;
    if(tp.Pos[0] > 0) wire0 = std::nearbyint(tp.Pos[0]);
    if(wire0 > slc.nWires[thirdPlane]) wire0 = slc.nWires[thirdPlane];
    // Do the same for the end
    unsigned short lastPt = pfp.TP3Ds.size() - 1;
    tp = MakeBareTP(slc, pfp.TP3Ds[lastPt].Pos, inCTP);
    unsigned int wire1 = 0;
    if(tp.Pos[0] > 0) wire1 = std::nearbyint(tp.Pos[0]);
    if(wire1 > slc.nWires[thirdPlane]) wire1 = slc.nWires[thirdPlane];
    if(wire0 == wire1) return !evt.goodWire[thirdPlane][wire0];
    if(wire1 < wire0) std::swap(wire0, wire1);
    // count the number of good wires
    int dead = 0;
    int wires = wire1 - wire0;
    for(unsigned int wire = wire0; wire < wire1; ++wire) if(!evt.goodWire[thirdPlane][wire]) ++dead;
    // require that most of the wires are dead
    return (dead > 0.8 * wires);
  } // ValidTwoPlaneMatch

  /////////////////////////////////////////
  void AddPointsInRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, 
                        CTP_t inCTP, float maxPull, unsigned short& nWires, unsigned short& nAdd, bool prt)
  {
    // Try to insert 2D trajectory points into the 3D trajectory point vector pfp.TP3Ds. 
    // This function inserts new TP3Ds and sets the NeedsUpdate flags true. 
    // The calling function should call Update
    nWires = 0;
    nAdd = 0;
    if(fromPt > toPt) return;
    if(toPt >= pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size() - 1;
    
    // Make a TP in this plane using the fromPt 3D position. This code assumes that the
    // fromPt and toPt 3D positions are in the same section or alternatively that the
    // fits aren't too dissimilar if they are in different sections. We will move this
    // tp along the trajectory direction in this CTP
    auto fromTP = MakeBareTP(slc, pfp.TP3Ds[fromPt].Pos, pfp.TP3Ds[fromPt].Dir, inCTP);
//    if(prt) std::cout<<"APIR: inCTP "<<inCTP<<" fromPt "<<fromPt<<" "<<PrintPos(slc, fromTP);
    unsigned short plane = DecodeCTP(inCTP).Plane;
    if(fromTP.Pos[0] < 0.) MoveTPToWire(fromTP, 0.);
    if(fromTP.Pos[0] > slc.nWires[plane]) MoveTPToWire(fromTP, (float)slc.nWires[plane]);
    // and another using the toPt 3D position
    auto toTP = MakeBareTP(slc, pfp.TP3Ds[toPt].Pos, pfp.TP3Ds[toPt].Dir, inCTP);
//    if(prt) std::cout<<" toPt "<<toPt<<" "<<PrintPos(slc, toTP)<<"\n";
    if(toTP.Pos[0] < 0.) MoveTPToWire(toTP, 0.);
    if(toTP.Pos[0] > slc.nWires[plane]) MoveTPToWire(toTP, (float)slc.nWires[plane]);
    // We now have two 2D points that may have been found using 3D -> 2D positions, perhaps in
    // different Sections. Use fromTP which has direction information to decide which wires
    // to consider
    if(toTP.Pos[0] < fromTP.Pos[0]) {
      std::swap(fromTP, toTP);
    }
    fromTP.Dir = PointDirection(fromTP.Pos, toTP.Pos);
    fromTP.Ang = atan2(fromTP.Dir[1], fromTP.Dir[0]);
    SetAngleCode(fromTP);
    unsigned int fromWire = std::nearbyint(fromTP.Pos[0]);
    unsigned int toWire = std::nearbyint(toTP.Pos[0]);
//    if(prt) std::cout<<"APIR P"<<pfp.ID<<" CTP "<<inCTP<<" fromWire "<<fromWire<<" toWire "<<toWire<<"\n";
    if(fromWire > toWire) return;
    nWires = 1 + toWire - fromWire;
    // vector of already used Tj,TPIndex pairs
    std::vector<std::pair<int, unsigned short>> tpsUsed;
    // and wires that have TPs
    std::vector<unsigned int> wiresUsed;
    // populate the vectors
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.CTP != inCTP) continue;
      tpsUsed.push_back(std::make_pair(tp3d.TjID, tp3d.TPIndex));
      wiresUsed.push_back(std::nearbyint(tp3d.Wire));
    } // tp3d
    // set a generous search window in WSE units
    float window = 5;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"APIR: P"<<pfp.ID<<" plane "<<plane;
      myprt<<" fromPt "<<fromPt<<" toPt "<<toPt;
      myprt<<" TP range "<<PrintPos(slc, fromTP)<<" to "<<PrintPos(slc, toTP);
      myprt<<" AngleCode "<<fromTP.AngleCode;
      myprt<<" Hit search window "<<window/tcc.unitsPerTick<<" ticks";
      myprt<<" pull cut "<<maxPull;
    }

    // Move ltp to each wire and check for hits that are used in a trajectory
    for(unsigned int wire = fromWire; wire <= toWire; ++wire) {
      if(std::find(wiresUsed.begin(), wiresUsed.end(), wire) != wiresUsed.end()) continue;
      MoveTPToWire(fromTP, (float)wire);
//      if(prt) std::cout<<"APIR stp "<<PrintPos(slc, tp)<<"\n";
      if(!FindCloseHits(slc, fromTP, window, kUsedHits)) continue;
      if(fromTP.Environment[kEnvDeadWire]) continue;
      // create a template TP3D. This will be over-written with the best
      TP3D bestTP3D;
      float bestPull = maxPull;
      for(auto iht : fromTP.Hits) {
        if(slc.slHits[iht].InTraj <= 0) continue;
        // this hit is used in a TP so find the tpIndex
        auto& utj = slc.tjs[slc.slHits[iht].InTraj - 1];
        unsigned short tpIndex = 0;
        for(tpIndex = utj.EndPt[0]; tpIndex <= utj.EndPt[1]; ++tpIndex) {
          auto& utp = utj.Pts[tpIndex];
          if(utp.Chg <= 0) continue;
          // This doesn't check for UseHit true but that is probably ok here
          if(std::find(utp.Hits.begin(), utp.Hits.end(), iht) != utp.Hits.end()) break;
        } // ipt
        if(tpIndex > utj.EndPt[1]) continue;
        auto npr = std::make_pair(utj.ID, tpIndex);
        if(std::find(tpsUsed.begin(), tpsUsed.end(), npr) != tpsUsed.end()) continue;
//        tpsUsed.push_back(npr);
        auto& utp = utj.Pts[tpIndex];
        // see if it is used in a different PFP
        if(utp.InPFP > 0) continue;
        auto newTP3D = CreateTP3D(slc, utj.ID, tpIndex);
        if(!SetSection(slc, pfp, newTP3D)) continue;
        float pull = PointPull(pfp, newTP3D);
        if(prt && pull < 10) {
          mf::LogVerbatim myprt("TC");
          myprt<<"APIR: P"<<pfp.ID<<" TP "<<PrintHit(slc.slHits[iht])<<" pull "<<pull<<" dx "<<newTP3D.TPX - newTP3D.Pos[0]<<" in section "<<newTP3D.SFIndex;
          auto mcpi = FindMCPIndex(slc, newTP3D);
          if(mcpi != UINT_MAX) myprt<<" mcpIndex "<<mcpi;
        }
        if(pull > bestPull) continue;
        bestTP3D = newTP3D;
        bestPull = pull;
      } // iht
      if(bestPull < maxPull) {
        if(InsertTP3D(pfp, bestTP3D) == USHRT_MAX) {
          std::cout<<"APIR: InsertTP3D failed\n";
          continue;
        }
        ++nAdd;
      } // bestPull < maxPull
    } // wire

  } // AddPointsInRange

  /////////////////////////////////////////
  unsigned short InsertTP3D(PFPStruct& pfp, TP3D& tp3d)
  {
    // inserts the tp3d into the section defined by tp3d.SFIndex
    if(tp3d.SFIndex >= pfp.SectionFits.size()) return USHRT_MAX;
    // Find the first occurrence of this SFIndex
    unsigned short ipt = 0;
    for(ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) if(tp3d.SFIndex == pfp.TP3Ds[ipt].SFIndex) break;
    if(ipt == pfp.TP3Ds.size()) return USHRT_MAX;
    // next see if we can insert it so that re-sorting of this section isn't required
    auto lastTP3D = pfp.TP3Ds.back();
    if(ipt == 0 && tp3d.along < pfp.TP3Ds[0].along) {
      // insert at the beginning. No search needs to be done
    } else if(tp3d.SFIndex == lastTP3D.SFIndex && tp3d.along > lastTP3D.along) {
      // insert at the end. Use push_back and return
      pfp.TP3Ds.push_back(tp3d);
      pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
      pfp.Flags[kNeedsUpdate] = true;
      return pfp.TP3Ds.size() - 1;
    } else {
      for(unsigned short iipt = ipt; iipt < pfp.TP3Ds.size() - 1; ++iipt) {
        // break out if the next point is in a different section
        if(pfp.TP3Ds[iipt + 1].SFIndex != tp3d.SFIndex) break;
        if(tp3d.along > pfp.TP3Ds[iipt].along && tp3d.along < pfp.TP3Ds[iipt + 1].along) {
          ipt = iipt + 1;
          break;
        }
      } // iipt
    } // insert in the middle
    pfp.TP3Ds.insert(pfp.TP3Ds.begin() + ipt, tp3d);
    pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
    pfp.Flags[kNeedsUpdate] = true;
    return ipt;
  } // InsertTP3D

  /////////////////////////////////////////
  bool SortSection(PFPStruct& pfp, unsigned short sfIndex)
  {
    // sorts the TP3Ds by the distance from the start of a fit section
    
    if(sfIndex > pfp.SectionFits.size() - 1) return false;
    auto& sf = pfp.SectionFits[sfIndex];
    if(sf.Pos[0] == 0.0 && sf.Pos[1] == 0.0 && sf.Pos[2] == 0.0) {
//      std::cout<<"P"<<pfp.ID<<" section fit position not defined\n";
      return false;
    }
    
    // a temp vector of points in this section
    std::vector<TP3D> temp;
    // and the index into TP3Ds
    std::vector<unsigned short> indx;
    // See if the along variable is monotonically increasing
    float prevAlong = 0;
    bool first = true;
    bool needsSort = false;
    for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      auto& tp3d = pfp.TP3Ds[ii];
      if(tp3d.SFIndex != sfIndex) continue;
      if(first) {
        first = false;
        prevAlong = tp3d.along;
      } else {
        if(tp3d.along < prevAlong) needsSort = true;
        prevAlong = tp3d.along;
      }
      temp.push_back(tp3d);
      indx.push_back(ii);
    } // tp3d
    if(temp.empty()) return false;
    // no sort needed?
    if(temp.size() == 1) return true;
    if(!needsSort) {
//      std::cout<<"SortSection: P"<<pfp.ID<<" section "<<sfIndex<<" doesn't need sorting\n";
      sf.NeedsUpdate = false;
      return true;
    }    
    // see if the points are not-contiguous
    bool contiguous = true;
    for(unsigned short ipt = 1; ipt < indx.size(); ++ipt) {
      if(indx[ipt] != indx[ipt - 1] + 1) {
        contiguous = false;
        std::cout<<"SortSection: MVI "<<pfp.MVI<<" Points aren't contiguous in sfi "<<sfIndex<< " ipt "<<ipt<<". print and quit\n";
        for(unsigned short ipt = 1; ipt < pfp.TP3Ds.size(); ++ipt) {
          auto& tp3d = pfp.TP3Ds[ipt];
          std::cout<<ipt<<" sfi "<<tp3d.SFIndex<<" along "<<tp3d.along<<" good? "<<tp3d.IsGood<<"\n";
        } // tp3d
      } // not contiguous
    } // ipt
    if(!contiguous) {
      return false;
    }
    
    std::vector<SortEntry> sortVec(temp.size());
    for(unsigned short ii = 0; ii < temp.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].val = temp[ii].along;
    } // ipt
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    for(unsigned short ii = 0; ii < temp.size(); ++ii) {
      // overwrite the tp3d
      auto& tp3d = pfp.TP3Ds[indx[ii]];
      tp3d = temp[sortVec[ii].index];
    } // ii
    sf.NeedsUpdate = false;
    return true;
  } // SortSection
  
  /////////////////////////////////////////
  void MakeTP3Ds(TCSlice& slc, PFPStruct& pfp)
  {
    // Create and populate the TP3Ds vector. This function is called before the first
    // fit is done so the TP3D along variable can't be determined
    if(!pfp.TP3Ds.empty() || pfp.SectionFits.size() != 1) return;

    // Add the points associated with the Tjs that were used to create the PFP
    for(auto tid : pfp.TjIDs) {
      auto& tj = slc.tjs[tid - 1];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(tp.InPFP > 0) continue;
        auto tp3d = CreateTP3D(slc, tid, ipt);
        tp3d.SFIndex = 0;
        // We need to assume that all points are good or the first fit will fail
        tp3d.IsGood = true;
        // unless the tp is on a very large angle trajectory
        if(tp.AngleCode == 2) tp3d.IsGood = false;
        pfp.TP3Ds.push_back(tp3d);
      } // ipt
    } // tid
  } // MakeTP3Ds

  /////////////////////////////////////////
  void Reverse(TCSlice& slc, PFPStruct& pfp)
  {
    // reverse the PFParticle
    std::reverse(pfp.TP3Ds.begin(), pfp.TP3Ds.end());
    std::reverse(pfp.SectionFits.begin(), pfp.SectionFits.end());
    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      // flip the direction vector
      for(unsigned short xyz = 0; xyz < 3; ++xyz) sf.Dir[xyz] *= -1;
    } // sf
    // correct the along variable
    for(auto& tp3d : pfp.TP3Ds) tp3d.along *= -1;
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    std::swap(pfp.EndFlag[0], pfp.EndFlag[1]);
  } // Reverse

  /////////////////////////////////////////
  void FillmAllTraj(TCSlice& slc) 
  {
    // Fills the mallTraj vector with trajectory points in the tpc and sorts
    // them by increasing X
    
    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;

    // define mallTraj
    slc.mallTraj.clear();
    Tj2Pt tj2pt;
    unsigned short cnt = 0;

    float rms = tcc.match3DCuts[0];
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      // ignore already matched
      if(tj.AlgMod[kMat3D]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      int plane = planeID.Plane;
      int tjID = tj.ID;
      if(tjID <= 0) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(tp.Pos[0] < -0.4) continue;
        // ignore already matched
        if(tp.InPFP > 0) continue;
        tj2pt.wire = std::nearbyint(tp.Pos[0]);
        ++cnt;
        // don't try matching if the wire doesn't exist
        if(!tcc.geom->HasWire(geo::WireID(cstat, tpc, plane, tj2pt.wire))) continue;
        float xpos = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, plane, tpc, cstat);
        tj2pt.xlo = xpos - rms;
        tj2pt.xhi = xpos + rms;
        tj2pt.plane = plane;
        tj2pt.id = tjID;
        tj2pt.ipt = ipt;
        tj2pt.npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        slc.mallTraj.push_back(tj2pt);
      } // tp
    } // tj
    
    // sort by increasing x
    std::vector<SortEntry> sortVec(slc.mallTraj.size());
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size(); ++ipt) {
      // populate the sort vector
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = slc.mallTraj[ipt].xlo;
    } // ipt
    // sort by increasing xlo
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put slc.mallTraj into sorted order
    auto tallTraj = slc.mallTraj;
    for(unsigned int ii = 0; ii < sortVec.size(); ++ii) slc.mallTraj[ii] = tallTraj[sortVec[ii].index];
    
  } // FillmAllTraj

  /////////////////////////////////////////
  bool SharesHighScoreVx(TCSlice& slc, const PFPStruct& pfp, const Trajectory& tj)
  {
    // returns true if tj with tjID shares a high-score 3D vertex with any
    // tj in pfp.TjIDs
    for(unsigned short end = 0; end < 2; ++end) {
      if(tj.VtxID[end] == 0) continue;
      auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
      if(!vx2.Stat[kHiVx3Score]) continue;
      std::vector<int> vtjlist = GetVtxTjIDs(slc, vx2);
      auto shared = SetIntersection(vtjlist, pfp.TjIDs);
      if(!shared.empty()) return true;
    } // end
    return false;
  } // SharesHighScoreVx
  

  ////////////////////////////////////////////////
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2)
  {
    if(v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) return 0;
    return acos(DotProd(v1, v2));
  } 
  
  ////////////////////////////////////////////////
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2)
  {
    // Finds the direction vector between the two points from p1 to p2
    Vector3_t dir;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) dir[xyz] = p2[xyz] - p1[xyz];
    if(dir[0] == 0 && dir[1] == 0 && dir[2] == 0) return dir;
    if(!SetMag(dir, 1)) { dir[0] = 0; dir[1] = 0; dir[3] = 0; }
    return dir;
  } // PointDirection

  //////////////////////////////////////////
  double PosSep(const Point3_t& pos1, const Point3_t& pos2)
  {
    return sqrt(PosSep2(pos1, pos2));
  } // PosSep
  
  //////////////////////////////////////////
  double PosSep2(const Point3_t& pos1, const Point3_t& pos2)
  {
    // returns the separation distance^2 between two positions in 3D
    double d0 = pos1[0] - pos2[0];
    double d1 = pos1[1] - pos2[1];
    double d2 = pos1[2] - pos2[2];
    return d0*d0 + d1*d1 + d2*d2;
  } // PosSep2
  
  //////////////////////////////////////////
  bool SetMag(Vector3_t& v1, double mag)
  {
    double den = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
    if(den == 0) return false;
    den = sqrt(den);
    
    v1[0] *= mag / den;
    v1[1] *= mag / den;
    v1[2] *= mag / den;
    return true;
  } // SetMag

  /////////////////////////////////////////
  void FilldEdx(TCSlice& slc, PFPStruct& pfp)
  {
    // Fills dE/dx variables in the pfp struct
    // TODO: Do this correctly instead of just quickly
    
    // don't attempt to find dE/dx at the end of a shower
    unsigned short numEnds = 2;
    if(pfp.PDGCode == 1111) numEnds = 1;
    
    // set dE/dx to 0 to indicate that a valid dE/dx is expected
    for(unsigned short end = 0; end < numEnds; ++end) {
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) pfp.dEdx[end][plane] = 0;
    } // end
    
    // square of the maximum length that is used for finding the average dE/dx
    float maxSep2 = 5 * tcc.wirePitch;
    maxSep2 *= maxSep2;
    
    for(unsigned short end = 0; end < numEnds; ++end) {
      std::vector<float> cnt(slc.nPlanes);
      short dir = 1 - 2 * end;
      auto endPos = PosAtEnd(pfp, end);
      for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
        unsigned short ipt;
        if(dir > 0) {
          ipt = ii;
        } else {
          ipt = pfp.TP3Ds.size() - ii - 1;
        }
        if(ipt >= pfp.TP3Ds.size()) break;
        auto& tp3d = pfp.TP3Ds[ipt];
        if(tp3d.IsBad) continue;
        if(PosSep2(tp3d.Pos, endPos) > maxSep2) break;
        // require good points
        if(!tp3d.IsGood) continue;
        float dedx = dEdx(slc, tp3d);
        if(dedx < 0.5) continue;
        unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
        pfp.dEdx[end][plane] += dedx;
        ++cnt[plane];
      } // ii
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        if(cnt[plane] == 0) continue;
        pfp.dEdx[end][plane] /= cnt[plane];
      } // plane
    } // end

  } // FilldEdx

  /////////////////////////////////////////
  float dEdx(TCSlice& slc, TP3D& tp3d)
  {
    if(!tp3d.IsGood) return 0;
    if(tp3d.TjID > (int)slc.slHits.size()) return 0;
    if(tp3d.TjID <= 0) return 0;

    double dQ = 0.;
    double time = 0;
    geo::PlaneID plnID;
    auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      auto& hit = (*evt.allHits)[slc.slHits[tp.Hits[ii]].allHitsIndex];
      dQ += hit.Integral();
    } // ii
    time = tp.Pos[1] / tcc.unitsPerTick;
    plnID = DecodeCTP(tp.CTP);
    if(dQ == 0) return 0;
    double angleToVert = tcc.geom->Plane(plnID).ThetaZ() - 0.5 * ::util::pi<>();
    double cosgamma = std::abs(std::sin(angleToVert) * tp3d.Dir[1] + std::cos(angleToVert) * tp3d.Dir[2]);
    if(cosgamma < 1.E-5) return 0;
    double dx = tcc.geom->WirePitch(plnID) / cosgamma;
    double dQdx = dQ / dx;
    double t0 = 0;
    float dedx = tcc.caloAlg->dEdx_AREA(dQdx, time, plnID.Plane, t0);
    if(std::isinf(dedx)) dedx = 0;
    return dedx;
  } // dEdx

  ////////////////////////////////////////////////
  TP3D CreateTP3D(TCSlice& slc, int tjID, unsigned short tpIndex)
  {
    // create a TP3D with a single TP. Note that the SectionFit in which it
    // should be placed and the 3D position can't be determined until the the TP3D is 
    // associated with a pfp. See SetSection()
    
    if(tjID <= 0 || tjID > (int)slc.tjs.size()) {
      std::cout<<"bad tjID\n";
      exit(1);
    }
    auto& tj = slc.tjs[tjID - 1];
    if(tpIndex < tj.EndPt[0] || tpIndex > tj.EndPt[1]) {
      std::cout<<"bad tpIndex\n";
      exit(1);
    }
    
    TP3D tp3d;
    tp3d.TjID = tjID;
    tp3d.TPIndex = tpIndex;
    auto& tp2 = tj.Pts[tp3d.TPIndex];
    auto plnID = DecodeCTP(tp2.CTP);
    tp3d.CTP = tp2.CTP;
    double tick = tp2.HitPos[1]/tcc.unitsPerTick;
    tp3d.TPX = tcc.detprop->ConvertTicksToX(tick, plnID);
    // Get the RMS of the TP in WSE units and convert to cm
    float rms = TPHitsRMSTime(slc, tp2, kAllHits) * tcc.wirePitch;
    // inflate the error for large angle TPs
    if(tp2.AngleCode == 1) rms *= 2;
    // a more careful treatment for long-pulse hits
    if(tp2.AngleCode > 1) {
      std::vector<unsigned int> hitMultiplet;
      for(unsigned short ii = 0; ii < tp2.Hits.size(); ++ii) {
        if(!tp2.UseHit[ii]) continue;
        GetHitMultiplet(slc, tp2.Hits[ii], hitMultiplet);
        if(hitMultiplet.size() > 1) break;
      } // ii
      rms = HitsRMSTime(slc, hitMultiplet, kAllHits) * tcc.wirePitch;
      // the returned RMS is closer to the FWHM, so divide by 2
      rms /= 2;
    } // tp2.AngleCode > 1
    tp3d.TPXErr2 = rms * rms;
    tp3d.Wire = tp2.Pos[0];
    // Can't declare it good since Pos and SFIndex aren't defined
    tp3d.IsGood = false;
    return tp3d;
  } // CreateTP3D

  /////////////////////////////////////////
  bool SetSection(TCSlice& slc, PFPStruct& pfp, TP3D& tp3d)
  {
    // Determine which SectionFit this tp3d should reside, then calculate
    // the 3D position and the distance from the center of the SectionFit

    if(tp3d.Wire < 0) return false;
    if(pfp.SectionFits.empty()) return false;
    if(pfp.SectionFits[0].Pos[0] == -10.0) return false;
    
    auto plnID = DecodeCTP(tp3d.CTP);
    
    if(pfp.SectionFits.size() == 1) {
      tp3d.SFIndex = 0;
    } else {
      // Find the section center that is closest to this point in the wire coordinate
      float best = 1E6;
      for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        auto& sf = pfp.SectionFits[sfi];
        float sfWire = tcc.geom->WireCoordinate(sf.Pos[1], sf.Pos[2], plnID);
        float sep = std::abs(sfWire - tp3d.Wire);
        if(sep < best) {
          best = sep;
          tp3d.SFIndex = sfi;
        }
      } // sfi
    } // pfp.SectionFits.size() > 1
    auto& sf = pfp.SectionFits[tp3d.SFIndex];
    auto plnTP = MakeBareTP(slc, sf.Pos, sf.Dir, tp3d.CTP);
    // the number of wires relative to the SectionFit center
    double dw = tp3d.Wire - plnTP.Pos[0];
    // dt/dW was stored in DeltaRMS
    double t = dw * plnTP.DeltaRMS;
    // define the 3D position 
    for(unsigned short xyz = 0; xyz < 3; ++xyz) tp3d.Pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
    tp3d.along = t;
    tp3d.IsGood = true;
    return true;
  } // SetSection
  
  ////////////////////////////////////////////////
  float PointPull(const PFPStruct& pfp, TP3D& tp3d)
  {
    // returns the pull that the tp3d will cause in the pfp section fit. This
    // currently only uses position but eventually will include charge
    return std::abs(tp3d.Pos[0] - tp3d.TPX) / sqrt(tp3d.TPXErr2);
  } // PointPull

  ////////////////////////////////////////////////
  PFPStruct CreatePFP(TCSlice& slc)
  {
    // The calling function should define the size of pfp.TjIDs
    PFPStruct pfp;
    pfp.ID = slc.pfps.size() + 1;
    pfp.ParentUID = 0;
    pfp.TPCID = slc.TPCID;
    // initialize arrays for both ends
    if(slc.nPlanes < 4) {
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
  void PFPVertexCheck(TCSlice& slc)
  {
    // Ensure that all PFParticles have a start vertex. It is possible for
    // PFParticles to be attached to a 3D vertex that is later killed.
    if(!slc.isValid) return;

    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.Vx3ID[0] > 0) continue;
      if(pfp.SectionFits.empty()) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      vx3.Vx2ID.resize(slc.nPlanes);
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      Point3_t startPos;
      if(pfp.TP3Ds.empty()) {
        // must be a neutrino pfp
        startPos = pfp.SectionFits[0].Pos;
      } else if(!pfp.TP3Ds.empty()) {
        // normal pfp
        startPos = pfp.TP3Ds[0].Pos;
      }
      vx3.X = startPos[0];
      vx3.Y = startPos[1];
      vx3.Z = startPos[2];
      vx3.ID = slc.vtx3s.size() + 1;
      vx3.Primary = false;
      ++evt.globalP_UID;
      vx3.UID = evt.globalP_UID;
      slc.vtx3s.push_back(vx3);
//      std::cout<<"PFPVertexCheck: P"<<pfp.ID<<" create 3V"<<vx3.ID<<"\n";
      pfp.Vx3ID[0] = vx3.ID;
    } // pfp
  } // PFPVertexCheck

  /////////////////////////////////////////
  void DefinePFPParents(TCSlice& slc, bool prt)
  {
    /*
     This function reconciles vertices, PFParticles and slc, then
     defines the parent (j) - daughter (i) relationship and PDGCode. Here is a
     description of the conventions:

     V1 is the highest score 3D vertex in this tpcid so a neutrino PFParticle P1 is defined.
     V4 is a high-score vertex that has lower score than V1. It is declared to be a
     primary vertex because its score is higher than V5 and it is not associated with the
     neutrino interaction
     V6 was created to adhere to the convention that all PFParticles, in this case P9,
     be associated with a start vertex. There is no score for V6. P9 is it's own parent
     but is not a primary PFParticle.

     P1 - V1 - P2 - V2 - P4 - V3 - P5        V4 - P6                  V6 - P9
     \                                  \
     P3                                 P7 - V5 - P8

     The PrimaryID in this table is the ID of the PFParticle that is attached to the
     primary vertex, which may or may not be a neutrino interaction vertex.
     The PrimaryID is returned by the PrimaryID function
     PFP  parentID  DtrIDs     PrimaryID
     -----------------------------------
     P1     P1     P2, P3        P1
     P2     P1     P4            P2
     P3     P1     none          P3
     P4     P2     P5            P2
     P5     P4     none          P2

     P6     P6     none          P6
     P7     P7     P8            P7

     P9     P9     none          0

     */
    if(slc.pfps.empty()) return;
    if(tcc.modes[kTestBeam]) return;

    int neutrinoPFPID = 0;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      if(!tcc.modes[kTestBeam] && neutrinoPFPID == 0 && (pfp.PDGCode == 12 || pfp.PDGCode == 14)) {
        neutrinoPFPID = pfp.ID;
        break;
      }
    } // pfp

    // define the end vertex if the Tjs have end vertices
    constexpr unsigned short end1 = 1;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      // already done?
      if(pfp.Vx3ID[end1] > 0) continue;
      // ignore shower-like pfps
      if(IsShowerLike(slc, pfp.TjIDs)) continue;
      // count 2D -> 3D matched vertices
      unsigned short cnt3 = 0;
      unsigned short vx3id = 0;
      // list of unmatched 2D vertices that should be merged
      std::vector<int> vx2ids;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.VtxID[end1] == 0) continue;
        auto& vx2 = slc.vtxs[tj.VtxID[end1] - 1];
        if(vx2.Vx3ID == 0) {
          if(vx2.Topo == 1 && vx2.NTraj == 2) vx2ids.push_back(vx2.ID);
          continue;
        }
        if(vx3id == 0) vx3id = vx2.Vx3ID;
        if(vx2.Vx3ID == vx3id) ++cnt3;
      } // tjid
      if(cnt3 > 1) {
        // ensure it isn't attached at the other end
        if(pfp.Vx3ID[1 - end1] == vx3id) continue;
        pfp.Vx3ID[end1] = vx3id;
        if(cnt3 != slc.nPlanes && tcc.modes[kDebug]) {
//          std::cout<<"DPFPR: Missed an end vertex for PFP "<<pfp.ID<<" Write some code\n";
        }
      } // cnt3 > 1
    } // pfp

    // Assign a PDGCode to each PFParticle and look for a parent
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      // skip a neutrino PFParticle
      if(pfp.PDGCode == 12 || pfp.PDGCode == 14 || pfp.PDGCode == 22) continue;
      pfp.PDGCode = PDGCodeVote(slc, pfp.TjIDs, prt);
      // Define a PFP parent if there are two or more Tjs that are daughters of
      // Tjs that are used by the same PFParticle
      int pfpParentID = INT_MAX;
      unsigned short nParent = 0;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.ParentID <= 0) continue;
        auto parPFP = GetAssns(slc, "T", tj.ParentID, "P");
        if(parPFP.empty()) continue;
        if(pfpParentID == INT_MAX) pfpParentID = parPFP[0];
        if(parPFP[0] == pfpParentID) ++nParent;
      } // ii
      if(nParent > 1) {
        auto& ppfp = slc.pfps[pfpParentID - 1];
        // set the parent UID
        pfp.ParentUID = ppfp.UID;
        // add to the parent daughters list
        ppfp.DtrUIDs.push_back(pfp.UID);
      } // nParent > 1
    } // ipfp
    // associate primary PFParticles with a neutrino PFParticle
    if(neutrinoPFPID > 0) {
      auto& neutrinoPFP = slc.pfps[neutrinoPFPID - 1];
      int vx3id = neutrinoPFP.Vx3ID[1];
      for(auto& pfp : slc.pfps) {
        if(pfp.ID == 0 || pfp.ID == neutrinoPFPID) continue;
        if(pfp.Vx3ID[0] != vx3id) continue;
        pfp.ParentUID = (size_t)neutrinoPFPID;
//        pfp.Primary = true;
        neutrinoPFP.DtrUIDs.push_back(pfp.ID);
        if(pfp.PDGCode == 111) neutrinoPFP.PDGCode = 12;
      } // pfp
    } // neutrino PFP exists
  } // DefinePFPParents

  ////////////////////////////////////////////////
  bool StorePFP(TCSlice& slc, PFPStruct& pfp)
  {
    // stores the PFParticle in the slice
//    if(pfp.ID < int(slc.pfps.size())) return false;
    bool neutrinoPFP = (pfp.PDGCode == 12 || pfp.PDGCode == 14);
    if(!neutrinoPFP) {
      if(pfp.TjIDs.empty()) return false;
      if(pfp.PDGCode != 1111 && pfp.TP3Ds.size() < 2) return false;
    }
    // ensure that the InPFP flag is set
    unsigned short nNotSet = 0;
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.IsBad) continue;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.InPFP != pfp.ID) ++nNotSet;
    } // tp3d
    if(nNotSet > 0) {
      std::cout<<"StorePFP Warning: "<<nNotSet<<" TP3Ds have incorrect InPFP\n";
    } // nNotSet > 0
    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)slc.pfps.size() + 1) pfp.ID = slc.pfps.size() + 1;
    ++evt.globalP_UID;
    pfp.UID = evt.globalP_UID;
    
    // set the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      tj.AlgMod[kMat3D] = true;
    } // tjid
//    if(pfp.Flags[kNeedsUpdate]) std::cout<<"StorePFP: stored P"<<pfp.ID<<" but NeedsUpdate is true...\n";

    slc.pfps.push_back(pfp);
    return true;
  } // StorePFP
  
  ////////////////////////////////////////////////
  bool InsideFV(TCSlice& slc, PFPStruct& pfp, unsigned short end)
  {
    // returns true if the end of the pfp is inside the fiducial volume of the TPC
    if(pfp.ID <= 0) return false;
    if(end > 1) return false;
    if(pfp.SectionFits.empty()) return false;
    // require that the points are sorted which ensures that the start and end points
    // are the first and last points in the TP3Ds vector
    if(pfp.Flags[kNeedsUpdate]) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    
    float abit = 5;
    Point3_t pos;
    if(neutrinoPFP) {
      pos = pfp.SectionFits[0].Pos;
    } else if(end == 0) {
      pos = pfp.TP3Ds[0].Pos;
    } else {
      pos = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    }
    return (pos[0] > slc.xLo + abit && pos[0] < slc.xHi - abit && 
            pos[1] > slc.yLo + abit && pos[1] < slc.yHi - abit &&
            pos[2] > slc.zLo + abit && pos[2] < slc.zHi - abit);
    
  } // InsideFV

  ////////////////////////////////////////////////
  bool InsideTPC(const Point3_t& pos, geo::TPCID& inTPCID)
  {
    // determine which TPC this point is in. This function returns false
    // if the point is not inside any TPC
    float abit = 5;
    for (const geo::TPCID& tpcid: tcc.geom->IterateTPCIDs()) {
      const geo::TPCGeo& TPC = tcc.geom->TPC(tpcid);
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      TPC.LocalToWorld(local,world);
      // reduce the active area of the TPC by a bit to be consistent with FillWireHitRange
      if(pos[0] < world[0]-tcc.geom->DetHalfWidth(tpcid) + abit) continue;
      if(pos[0] > world[0]+tcc.geom->DetHalfWidth(tpcid) - abit) continue;
      if(pos[1] < world[1]-tcc.geom->DetHalfHeight(tpcid) + abit) continue;
      if(pos[1] > world[1]+tcc.geom->DetHalfHeight(tpcid) - abit) continue;
      if(pos[2] < world[2]-tcc.geom->DetLength(tpcid)/2 + abit) continue;
      if(pos[2] > world[2]+tcc.geom->DetLength(tpcid)/2 - abit) continue;
      inTPCID = tpcid;
      return true;
    } // tpcid
    return false;
  } // InsideTPC

  ////////////////////////////////////////////////
  void FindAlongTrans(Point3_t pos1, Vector3_t dir1, Point3_t pos2, Point2_t& alongTrans)
  {
    // Calculate the distance along and transvers to the direction vector from pos1 to pos2
    alongTrans[0] = 0;
    alongTrans[1] = 0;
    if(pos1[0] == pos2[0] && pos1[1] == pos2[1] && pos1[2] == pos2[2]) return;
    auto ptDir = PointDirection(pos1, pos2);
    SetMag(dir1, 1.0);
    double costh = DotProd(dir1, ptDir);
    if(costh > 1) costh = 1;
    double sep = PosSep(pos1, pos2);
    alongTrans[0] = costh * sep;
    double sinth = sqrt(1 - costh * costh);
    alongTrans[1] = sinth * sep;
  } // FindAlongTrans

  ////////////////////////////////////////////////
  bool PointDirIntersect(Point3_t p1, Vector3_t p1Dir, Point3_t p2, Vector3_t p2Dir, Point3_t& intersect, float& doca)
  {
    // Point - vector version
    Point3_t p1End, p2End;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) {
      p1End[xyz] = p1[xyz] + 10 * p1Dir[xyz];
      p2End[xyz] = p2[xyz] + 10 * p2Dir[xyz];
    }
    return LineLineIntersect(p1, p1End, p2, p2End, intersect, doca);
  } // PointDirIntersect

  ////////////////////////////////////////////////
  bool LineLineIntersect(Point3_t p1, Point3_t p2, Point3_t p3, Point3_t p4, Point3_t& intersect, float& doca)
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
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;
    constexpr double EPS = std::numeric_limits<double>::min();

    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    if (std::abs(p43[0]) < EPS && std::abs(p43[1]) < EPS && std::abs(p43[2]) < EPS) return(false);
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    if (std::abs(p21[0]) < EPS && std::abs(p21[1]) < EPS && std::abs(p21[2]) < EPS) return(false);

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < EPS) return(false);
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
    for(unsigned short xyz = 0; xyz < 3; ++xyz) intersect[xyz] += pb[xyz];
    for(unsigned short xyz = 0; xyz < 3; ++xyz) intersect[xyz] /= 2;
    return true;
  } // LineLineIntersect

  ////////////////////////////////////////////////
  float ChgFracBetween(TCSlice& slc, Point3_t pos1, Point3_t pos2)
  {
    // Step between pos1 and pos2 and find the fraction of the points that have nearby hits
    // in each plane. This function returns -1 if something is fishy, but this doesn't mean
    // that there is no charge. Note that there is no check for charge precisely at the pos1 and pos2
    // positions
    float sep = PosSep(pos1, pos2);
    if(sep == 0) return -1;
    unsigned short nstep = sep / tcc.wirePitch;
    auto dir = PointDirection(pos1, pos2);
    float sum = 0;
    float cnt = 0;
    TrajPoint tp;
    for(unsigned short step = 0; step < nstep; ++step) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pos1[xyz] += tcc.wirePitch * dir[xyz];
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        tp.CTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
        tp.Pos[0] = tcc.geom->WireCoordinate(pos1[1], pos1[2], plane, slc.TPCID.TPC, slc.TPCID.Cryostat);
        tp.Pos[1] = tcc.detprop->ConvertXToTicks(pos1[0], plane, slc.TPCID.TPC, slc.TPCID.Cryostat) * tcc.unitsPerTick;
        ++cnt;
        if(SignalAtTp(tp)) ++sum;
      } // plane
    } // step
    if(cnt == 0) return -1;
    return sum / cnt;

  } // ChgFracBetween

  ////////////////////////////////////////////////
  float ChgFracNearEnd(TCSlice& slc, PFPStruct& pfp, unsigned short end)
  {
    // returns the charge fraction near the end of the pfp. Note that this function
    // assumes that there is only one Tj in a plane.
    if(pfp.ID == 0) return 0;
    if(pfp.TjIDs.empty()) return 0;
    if(end < 0 || end > 1) return 0;
    if(pfp.TPCID != slc.TPCID) return 0;
    if(pfp.SectionFits.empty()) return 0;

    float sum = 0;
    float cnt = 0;
    // keep track of the lowest value and maybe reject it
    float lo = 1;
    float hi = 0;
    auto pos3 = PosAtEnd(pfp, end);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      std::vector<int> tjids(1);
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.CTP != inCTP) continue;
        tjids[0] = tjid;
        Point2_t pos2;
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        pos2[0] = tcc.geom->WireCoordinate(pos3[1], pos3[2], planeID);
        if(pos2[0] < -0.4) continue;
        // check for dead wires
        unsigned int wire = std::nearbyint(pos2[0]);
        if(wire > slc.nWires[plane]) continue;
        if(slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
        pos2[1] = tcc.detprop->ConvertXToTicks(pos3[0], planeID) * tcc.unitsPerTick;
        float cf = ChgFracNearPos(slc, pos2, tjids);
        if(cf < lo) lo = cf;
        if(cf > hi) hi = cf;
        sum += cf;
        ++cnt;
      } // tjid
    } // plane
    if(cnt == 0) return 0;
    if(cnt > 1 && lo < 0.3 && hi > 0.8) {
      sum -= lo;
      --cnt;
    }
    return sum / cnt;
  } // ChgFracNearEnd

  ////////////////////////////////////////////////
  Vector3_t DirAtEnd(const PFPStruct& pfp, unsigned short end)
  {
    if(end > 1 || pfp.SectionFits.empty()) return {{0., 0., 0.}};
    if(end == 0) return pfp.SectionFits[0].Dir;
    return pfp.SectionFits[pfp.SectionFits.size() - 1].Dir;
  } // PosAtEnd
  
  ////////////////////////////////////////////////
  Point3_t PosAtEnd(const PFPStruct& pfp, unsigned short end)
  {
    if(end > 1 || pfp.SectionFits.empty()) return {{0., 0., 0.}};
    // handle a neutrino pfp that doesn't have any TP3Ds
    if(pfp.TP3Ds.empty()) return pfp.SectionFits[0].Pos;
    if(end == 0) return pfp.TP3Ds[0].Pos;
    return pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
  } // PosAtEnd
  
  ////////////////////////////////////////////////
  float Length(const PFPStruct& pfp)
  {
    if(pfp.TP3Ds.empty()) return 0;
    return PosSep(pfp.TP3Ds[0].Pos, pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos);
  } // Length
  
  ////////////////////////////////////////////////
  bool SectionStartEnd(const PFPStruct& pfp, unsigned short sfIndex, unsigned short& startPt, unsigned short& endPt)
  {
    // this assumes that the TP3Ds vector is sorted
    startPt = USHRT_MAX;
    endPt = USHRT_MAX;
    if(sfIndex >= pfp.SectionFits.size()) return false;
    
    bool first = true;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.SFIndex < sfIndex) continue;
      if(first) {
        first = false;
        startPt = ipt;
      } // first
      if(tp3d.SFIndex > sfIndex) break;
      endPt = ipt;
    } // ipt
    return true;
    
  } // SectionStartEnd

  ////////////////////////////////////////////////
  unsigned short FarEnd(TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos)
  {
    // Returns the end (0 or 1) of the pfp that is furthest away from the position pos
    if(pfp.ID == 0) return 0;
    if(pfp.TP3Ds.empty()) return 0;
    auto& pos0 = pfp.TP3Ds[0].Pos;
    auto& pos1 = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    if(PosSep2(pos1, pos) > PosSep2(pos0, pos)) return 1;
    return 0;
  } // FarEnd
  
  ////////////////////////////////////////////////
  unsigned int FindMCPIndex(TCSlice& slc, TP3D tp3d)
  {
    // look for a mcp match
    if(evt.allHitsMCPIndex.empty()) return UINT_MAX;
    auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      unsigned ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
      return evt.allHitsMCPIndex[ahi];
    } // ii
    return UINT_MAX;
  } // FindMCPIndex

  ////////////////////////////////////////////////
  void PrintTP3Ds(std::string someText, TCSlice& slc, const PFPStruct& pfp, short printPts)
  {
    if(pfp.TP3Ds.empty()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" pfp P"<<pfp.ID<<" MVI "<<pfp.MVI;
    for(auto tid : pfp.TjIDs) myprt<<" T"<<tid;
    myprt<<" Flags: CanSection? "<<pfp.Flags[kCanSection];
    myprt<<" NeedsUpdate? "<<pfp.Flags[kNeedsUpdate];
    myprt<<" Algs:";
    for(unsigned short ib = 0; ib < pAlgModSize; ++ib) {
      if(pfp.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
    } // ib
    myprt<<"\n";
    if(!pfp.SectionFits.empty()) {
      myprt<<someText<<"  SFI ________Pos________   ________Dir_______ _____EndPos________ ChiDOF  NPts NeedsUpdate?\n";
      for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        myprt<<someText<<std::setw(4)<<sfi;
        auto& sf = pfp.SectionFits[sfi];
        myprt<<std::fixed<<std::setprecision(1);
        unsigned short startPt = 0, endPt = 0;
        if(SectionStartEnd(pfp, sfi, startPt, endPt)) {
          auto& start = pfp.TP3Ds[startPt].Pos;
          myprt<<std::setw(7)<<start[0]<<std::setw(7)<<start[1]<<std::setw(7)<<start[2];
        } else {
          myprt<<" Invalid";
        }
        myprt<<std::fixed<<std::setprecision(2);
        myprt<<std::setw(7)<<sf.Dir[0]<<std::setw(7)<<sf.Dir[1]<<std::setw(7)<<sf.Dir[2];
        myprt<<std::fixed<<std::setprecision(1);
        if(endPt < pfp.TP3Ds.size()) {
          auto& end = pfp.TP3Ds[endPt].Pos;
          myprt<<std::setw(7)<<end[0]<<std::setw(7)<<end[1]<<std::setw(7)<<end[2];
        } else {
          myprt<<" Invalid";
        }
        myprt<<std::setprecision(1)<<std::setw(6)<<sf.ChiDOF;
        myprt<<std::setw(6)<<sf.NPts;
        myprt<<std::setw(6)<<sf.NeedsUpdate;
        myprt<<"\n";
      } // sec
    } // SectionFits
    if(printPts < 0) {
      // print the head if we print all points
      myprt<<someText<<"  ipt SFI ________Pos________  Delta Pull GB?  Path along dE/dx Signal? KnkSig T_ipt_P:W:T     MCPIndex\n";
    }
    unsigned short fromPt = 0;
    unsigned short toPt = pfp.TP3Ds.size() - 1;
    if(printPts >= 0) fromPt = toPt;
    // temp kink angle for each point
    std::vector<float> dang(pfp.TP3Ds.size(), -1);
    for(unsigned short ipt = fromPt; ipt <= toPt; ++ipt) {
      auto tp3d = pfp.TP3Ds[ipt];
      myprt<<someText<<std::setw(4)<<ipt;
      myprt<<std::setw(4)<<tp3d.SFIndex;
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<std::setw(7)<<tp3d.Pos[0]<<std::setw(7)<<tp3d.Pos[1]<<std::setw(7)<<tp3d.Pos[2];
      myprt<<std::setprecision(1)<<std::setw(6)<<(tp3d.Pos[0] - tp3d.TPX);
      float pull = PointPull(pfp, tp3d);
      myprt<<std::setprecision(1)<<std::setw(6)<<pull;
      myprt<<std::setw(3)<<tp3d.IsGood<<tp3d.IsBad;
      myprt<<std::setw(7)<<std::setprecision(1)<<PosSep(tp3d.Pos, pfp.TP3Ds[0].Pos);
      myprt<<std::setw(7)<<std::setprecision(1)<<tp3d.along;
      myprt<<std::setw(6)<<std::setprecision(2)<<dEdx(slc, tp3d);
      // print SignalAtTP in each plane
      myprt<<" ";
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        auto tp = MakeBareTP(slc, tp3d.Pos, inCTP);
        myprt<<" "<<SignalAtTp(tp);
      } // plane
      
      double dang = -1, dangSig = -1;
      KinkFit(slc, pfp, ipt, tcc.kinkCuts[4], dang, dangSig);
      myprt<<std::setw(7)<<std::setprecision(3)<<dang;
      // print the 2D TP kink significance
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      myprt<<std::setw(7)<<std::setprecision(2)<<tp.KinkSig * dang;
      if(tp3d.TjID > 0) {
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        myprt<<" T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"_"<<PrintPos(slc, tp);
      } else {
        myprt<<" UNDEFINED";
      }
      unsigned int mcpIndex = FindMCPIndex(slc, tp3d);
      if(mcpIndex != UINT_MAX) myprt<<" "<<mcpIndex;
      myprt<<"\n";
    } // ipt
  } // PrintTP3Ds
} // namespace
