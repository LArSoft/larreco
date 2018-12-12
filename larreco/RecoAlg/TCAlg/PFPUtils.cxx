#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
// temp include for study
#include "larreco/RecoAlg/TCAlg/TCTruth.h"

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
              auto& pos1 = p1.XYZ[e1];
              // require the end to be close to a TPC boundary
              if(InsideFV(slc1, p1, e1)) continue;
              auto& dir1 = p1.Dir[e1];
              for(unsigned short e2 = 0; e2 < 2; ++e2) {
                auto& pos2 = p2.XYZ[e2];
                // require the end to be close to a TPC boundary
                if(InsideFV(slc2, p2, e2)) continue;
                auto& dir2 = p2.Dir[e2];
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
          if(pfp.XYZ[end][2] < minZ) { minZ = pfp.XYZ[end][2]; minZIndx = slcIndex;  minZEnd = end; }
        } // end
      } // puid
      if(minZEnd > 1) continue;
      // preserve the pfp with the min Z position
      auto& pfp = slices[minZIndx.first].pfps[minZIndx.second];
      if(prt) mf::LogVerbatim("TC")<<"SP: P"<<pfp.UID;
      // reverse it if necessary
      if(minZEnd != 0) ReversePFP(slices[minZIndx.first], pfp);
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
  void UpdateMatchStructs(TCSlice& slc, int oldTj, int newTj)
  {
    // Replaces tjid and ipt references in slc.matchVec and slc.pfps from
    // oldTj to newTj. This function is called when Tjs are split or merged
    // or if a Tj is reversed (in which case oldTj = newTj).
    // The method used is to match the trajectory point positions
    if(oldTj <= 0 || oldTj > (int)slc.tjs.size()) return;
    if(newTj <= 0 || newTj > (int)slc.tjs.size()) return;
    if(slc.mallTraj.empty() && slc.pfps.empty()) return;
    
    // convert from int to unsigned short
    int oldtjid = oldTj;
    int newtjid = newTj;
    auto& ntj = slc.tjs[newTj - 1];
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;
    // put the X positions of the new Tj into a vector for matching
    std::vector<float> xpos(ntj.Pts.size());
    geo::PlaneID planeID = DecodeCTP(ntj.CTP);
    for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
      auto& ntp = ntj.Pts[npt];
      if(ntp.Chg <= 0) continue;
      xpos[npt] = tcc.detprop->ConvertTicksToX(ntp.Pos[1]/tcc.unitsPerTick, planeID);
    } // npt
    
    if(!slc.mallTraj.empty()) {
      for(unsigned int ipt = 0; ipt < slc.mallTraj.size(); ++ipt) {
        auto& tj2pt = slc.mallTraj[ipt];
        if(tj2pt.id > slc.tjs.size()) continue;
        if(tj2pt.id != oldtjid) continue;
        // Found the old Tj. Now find the point
        for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
          auto& ntp = ntj.Pts[npt];
          if(ntp.Chg <= 0) continue;
          if(std::nearbyint(ntp.Pos[0]) == tj2pt.wire && xpos[npt] > tj2pt.xlo && xpos[npt] < tj2pt.xhi) {
            tj2pt.id = newtjid;
            tj2pt.ipt = npt;
            tj2pt.npts = npts;
            break;
          } // points match
        } // npt
      } // ipt
    } // !slc.mallTraj.empty()
    
    // Update pfp space points
    if(!slc.pfps.empty()) {
      for(auto& pfp : slc.pfps) {
        for(auto& tp3 : pfp.Tp3s) {
          // check each of the Tj2Pts associated with this space point
          for(auto& tj2pt : tp3.Tj2Pts) {
            if(tj2pt.id > slc.tjs.size()) continue;
            if(tj2pt.id != oldtjid) continue;
            // look for the corresponding point (wire) on the new Tj
            for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
              auto& ntp = ntj.Pts[npt];
              if(std::nearbyint(ntp.Pos[0]) == tj2pt.wire && xpos[npt] > tj2pt.xlo && xpos[npt] < tj2pt.xhi) {
                tj2pt.id = newtjid;
                tj2pt.ipt = npt;
                tj2pt.npts = npts;
                break;
              }
            } // npt
          } // tj2pt
        } // tp3
      } // pfp
    } // pfps exists

  } // UpdateMatchStructs
  
  /////////////////////////////////////////
  void UpdateTp3s(TCSlice& slc, PFPStruct& pfp, int oldTj, int newTj)
  {
    // Replaces occurrences of oldTj with newTj in the pfp vector of Tp3s
    if(oldTj <= 0 || oldTj > (int)slc.tjs.size()) return;
    if(newTj <= 0 || newTj > (int)slc.tjs.size()) return;
    if(slc.mallTraj.empty() && pfp.Tp3s.empty()) return;
    
    // convert from int to unsigned short
    unsigned short oldtjid = oldTj;
    unsigned short newtjid = newTj;
    auto& ntj = slc.tjs[newTj - 1];
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;
    // put the X positions of the new Tj into a vector for matching
    std::vector<float> xpos(ntj.Pts.size());
    geo::PlaneID planeID = DecodeCTP(ntj.CTP);
    for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
      auto& ntp = ntj.Pts[npt];
      if(ntp.Chg <= 0) continue;
      xpos[npt] = tcc.detprop->ConvertTicksToX(ntp.Pos[1]/tcc.unitsPerTick, planeID);
    } // npt

    for(auto& tp3 : pfp.Tp3s) {
      // check each of the Tj2Pts associated with this space point
      for(auto& tj2pt : tp3.Tj2Pts) {
        if(tj2pt.id > slc.tjs.size()) continue;
        if(tj2pt.id != oldtjid) continue;
        // look for the corresponding point (wire) on the new Tj
        for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
          auto& ntp = ntj.Pts[npt];
          if(std::nearbyint(ntp.Pos[0]) == tj2pt.wire && xpos[npt] > tj2pt.xlo && xpos[npt] < tj2pt.xhi) {
            tj2pt.id = newtjid;
            tj2pt.ipt = npt;
            tj2pt.npts = npts;
            break;
          }
        } // npt
      } // tj2pt
    } // tp3
    
  } // UpdateTp3s
  
  /////////////////////////////////////////
  void FillmAllTraj(TCSlice& slc) 
  {
    // Fills the mallTraj vector with trajectory points in the tpc and sorts
    // them by increasing X
    slc.matchVec.clear();
    
    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;
    
    // count the number of TPs and clear out any old 3D match flags
    unsigned int ntp = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      if(tj.ID <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      ntp += NumPtsWithCharge(slc, tj, false);
      tj.AlgMod[kMat3D] = false;
    } // tj
    if(ntp < 2) return;
    
    slc.mallTraj.resize(ntp);
    
    // define mallTraj
    unsigned int icnt = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      int plane = planeID.Plane;
      int tjID = tj.ID;
      if(tjID <= 0) continue;
      short score = 1;
      if(tj.AlgMod[kTjHiVx3Score]) score = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        if(icnt > slc.mallTraj.size() - 1) break;
        if(tp.Pos[0] < -0.4) continue;
        slc.mallTraj[icnt].wire = std::nearbyint(tp.Pos[0]);
        bool hasWire = tcc.geom->HasWire(geo::WireID(cstat, tpc, plane, slc.mallTraj[icnt].wire));
        // don't try matching if the wire doesn't exist
        if(!hasWire) continue;
        float xpos = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, plane, tpc, cstat);
        float posPlusRMS = tp.Pos[1] + TPHitsRMSTime(slc, tp, kUsedHits);
        float rms = tcc.detprop->ConvertTicksToX(posPlusRMS/tcc.unitsPerTick, plane, tpc, cstat) - xpos;
        if(rms < tcc.match3DCuts[0]) rms = tcc.match3DCuts[0];
        if(icnt == slc.mallTraj.size()) break;
        slc.mallTraj[icnt].xlo = xpos - rms;
        slc.mallTraj[icnt].xhi = xpos + rms;
        slc.mallTraj[icnt].dir = tp.Dir;
        slc.mallTraj[icnt].ctp = tp.CTP;
        slc.mallTraj[icnt].id = tjID;
        slc.mallTraj[icnt].ipt = ipt;
        slc.mallTraj[icnt].npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        slc.mallTraj[icnt].score = score;
        ++icnt;
      } // tp
    } // tj
    
    if(icnt < slc.mallTraj.size()) slc.mallTraj.resize(icnt);
    
    // This is pretty self-explanatory
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
  bool SetStart(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Analyzes the space point collection and the Tjs in the pfp to find a new start
    // position. The Tp3s found in FindCompleteness are ordered by increasing X. The general direction 
    // pfp.Dir[0] and the average position of all points in Tp3s was stored in pfp.XYZ[0]. This function
    // rotates each tp3 into this coordinate system to determine (along, trans) for each point. The min (max)
    // value of along defines the start (end) of the trajectory.
    
    if(pfp.ID <= 0 || pfp.TjIDs.empty()) return false;
    if(pfp.Tp3s.size() < 2) return false;

    // The projection along the general direction relative to the average position was found
    // in FillCompleteness. Now 
    float minAlong = 1E6;
    unsigned short minPt = 0;
    float maxAlong = -1E6;
    unsigned short maxPt = 0;
    std::vector<SortEntry> sortVec(pfp.Tp3s.size());
    for(unsigned short ipt = 0; ipt < pfp.Tp3s.size(); ++ipt) {
      auto& tp3 = pfp.Tp3s[ipt];
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = tp3.AlongTrans[0];
      // find the min (max) 
      if(tp3.AlongTrans[0] < minAlong) {
        minAlong = tp3.AlongTrans[0];
        minPt = ipt;
      }
      if(tp3.AlongTrans[0] > maxAlong) {
        maxAlong = tp3.AlongTrans[0];
        maxPt = ipt;
      }
    } // tp3
    
    pfp.XYZ[0] = pfp.Tp3s[minPt].Pos;
    pfp.XYZ[1] = pfp.Tp3s[maxPt].Pos;

    if(prt) {
      mf::LogVerbatim("TC")<<"SNS: P"<<pfp.ID<<" minPt "<<minPt<<" maxPt "<<maxPt<<" dir "<<std::fixed<<std::setprecision(2)<<pfp.Dir[0][0]<<" "<<pfp.Dir[0][1]<<" "<<pfp.Dir[0][2];
      PrintTp3("minPt", slc, pfp.Tp3s[minPt]);
      PrintTp3("maxPt", slc, pfp.Tp3s[maxPt]);
    }
    
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put them into order
    std::vector<TrajPoint3> temp;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) temp.push_back(pfp.Tp3s[sortVec[ii].index]);
    pfp.Tp3s = temp;
//    PrintTp3s("SNS", slc, pfp, -1);

    return true;
    
  } // SetStart

  /////////////////////////////////////////
  void FollowTp3s(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Step through the set of Tp3s on this pfp to create a trajectory. The start and end points
    // are assumed to be Tp3s[0] and Tp3s[Tp3s.size()-1] respectively. 
    
    if(pfp.Tp3s.size() < 2) return;
    
    unsigned short startPt = 0;
    unsigned short endPt = pfp.Tp3s.size() - 1;
    // divide the trajectory in 5 cm long sections. The average position of the Tp3s in
    // each section will be found. Tp3s 
    constexpr float sectionLen = 5;
    float endAlong = pfp.Tp3s[0].AlongTrans[0] + sectionLen;
    std::vector<Vector3_t> sectionPos;
    sectionPos.push_back(pfp.Tp3s[0].Pos);
    std::vector<unsigned short> sectionPt;
    sectionPt.push_back(0);
    for(unsigned short section = 0; section < 100; ++section) {
      // a point to find the average position in this section
      Point3_t avePos {{0,0,0}};
      Vector3_t aveDir {{0,0,0}};
      unsigned short cnt = 0;
      for(unsigned short ipt = startPt; ipt < endPt; ++ipt) {
        auto& tp3 = pfp.Tp3s[ipt];
        // The path length along the direction vector from the start point to the end
        // point was stashed in dEdxErr 
        // remove outliers
        if(tp3.AlongTrans[1] > 2) continue;
        if(tp3.AlongTrans[0] < endAlong) {
          // still in the same section - sum and continue
          for(unsigned short xyz = 0; xyz < 3; ++xyz) {
            avePos[xyz] += tp3.Pos[xyz];
            aveDir[xyz] += tp3.Dir[xyz];
          }
          ++cnt;
          continue;
        }
        // entered the next section. Check for a failure
        if(cnt == 0) continue;
        // calculate the average position
        for(unsigned short xyz = 0; xyz < 3; ++xyz) avePos[xyz] /= cnt;
        SetMag(aveDir, 1);
/*
        std::cout<<"Section "<<section<<" cnt "<<cnt<<" avePos"<<std::fixed<<std::setprecision(1);
        std::cout<<" "<<avePos[0]<<" "<<avePos[1]<<" "<<avePos[2];
        std::cout<<" aveDir "<<std::setprecision(2)<<aveDir[0]<<" "<<aveDir[1]<<" "<<aveDir[2]<<"\n";
*/
        sectionPos.push_back(avePos);
        sectionPt.push_back(ipt);
        startPt = ipt;
        endAlong += sectionLen;
        break;
      } // ipt
    } // section
    sectionPos.push_back(pfp.Tp3s[endPt].Pos);
    sectionPt.push_back(pfp.Tp3s.size() - 1);
/*
    for(unsigned short ipt = 0; ipt < sectionPos.size(); ++ipt) {
      std::cout<<ipt<<" sectionPt "<<sectionPt[ipt]<<" sectionPos "<<" "<<sectionPos[ipt][0]<<" "<<sectionPos[ipt][1]<<" "<<sectionPos[ipt][2]<<"\n";
    } // ipt
*/
    // set the general purpose flag bit false (unused) for all Tj Pts. This will be set true
    // when a Tp is used in a Tp3
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      for(auto& tp : tj.Pts) tp.Environment[kEnvFlag] = false;
    } // tjid
    // set the bits true for the first point
    for(auto tj2pt : pfp.Tp3s[0].Tj2Pts) {
      slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt].Environment[kEnvFlag] = true;
    }
    // create a vector of new Tp3s that will replace pfp.Tp3s
    std::vector<TrajPoint3> ntp3;
    ntp3.push_back(pfp.Tp3s[0]);
    // 2D position (WSE units) of the TPs at the start of this section. We will require that all 2D TPs are
    // less than sectionLen (in WSE units) from this point.
    std::vector<Point2_t> startPos2D(slc.nPlanes);
    // collect Tp3s in each section
    unsigned short lastPtAdded = 0;
    for(unsigned short section = 1; section < sectionPt.size(); ++section) {
      Point3_t startPos = sectionPos[section - 1];
      Point3_t endPos = sectionPos[section];
      auto startDir = PointDirection(startPos, endPos);
      // define the pfp start direction
      if(section == 1) pfp.Dir[0] = startDir;
      // and the end direction
      pfp.Dir[1] = startDir;
//      std::cout<<"Section "<<section<<" startDir "<<std::fixed<<std::setprecision(2)<<startDir[0]<<" "<<startDir[1]<<" "<<startDir[2]<<"\n";
      // define the 2D positions for this point in each plane
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        startPos2D[plane][0] = tcc.geom->WireCoordinate(sectionPos[section - 1][1], sectionPos[section - 1][2], planeID);
        startPos2D[plane][1] = tcc.detprop->ConvertXToTicks(sectionPos[section - 1][0], planeID) * tcc.unitsPerTick;
      } // plane
      for(unsigned short ipt = sectionPt[section - 1]; ipt < sectionPt[section]; ++ipt) {
        auto& tp3 = pfp.Tp3s[ipt];
        // count the number of Tps in this Tp3 that are already used in the trajectory
        unsigned short nused = 0;
        bool big2DSep = false;
        for(auto tj2pt : pfp.Tp3s[ipt].Tj2Pts) {
          auto& tp = slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt];
          if(tp.Environment[kEnvFlag]) ++nused;
          unsigned short plane = DecodeCTP(tp.CTP).Plane;
          float sep2D = PosSep(startPos2D[plane], tp.Pos) * tcc.wirePitch;
          if(sep2D > sectionLen) big2DSep = true;
        } // tj2pt
        if(big2DSep || nused > 1) continue;
        Point2_t alongTrans;
        FindAlongTrans(startPos, startDir, tp3.Pos, alongTrans);
        if(alongTrans[1] > 0.5) continue;
/*
        std::cout<<section<<" ipt "<<ipt<<" trans "<<alongTrans[1]<<" tj_ipt";
        for(auto tj2pt : tp3.Tj2Pts) std::cout<<" "<<tj2pt.id - 1<<"_"<<tj2pt.ipt;
        std::cout<<"\n";
*/
        tp3.AlongTrans = alongTrans;
        // don't clobber the original direction
//        tp3.Dir = dir;
        ntp3.push_back(tp3);
        // set the flag
        for(auto tj2pt : tp3.Tj2Pts) slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt].Environment[kEnvFlag] = true;
        lastPtAdded = ipt;
      } // ipt
    } // section
    
    if(lastPtAdded != endPt) ntp3.push_back(pfp.Tp3s[endPt]);
    
    if(prt) {
      float len = PosSep(ntp3[0].Pos, ntp3[ntp3.size()-1].Pos);
      mf::LogVerbatim("TC")<<"FollowTp3s: Tp3s size in "<<pfp.Tp3s.size()<<" size out "<<ntp3.size()<<" len "<<std::fixed<<std::setprecision(2)<<len;
    }
    
    // reverse if necessary to be consistent with a vertex
    if(pfp.Vx3ID[0] > 0) {
      auto& vx3 = slc.vtx3s[pfp.Vx3ID[0] - 1];
      Point3_t vpos = {{vx3.X, vx3.Y, vx3.Z}};
      auto& firstTp3Pos = ntp3[0].Pos;
      auto& lastTp3Pos = ntp3[ntp3.size() - 1].Pos;
      if(PosSep2(lastTp3Pos, vpos) < PosSep2(firstTp3Pos, vpos)) std::reverse(ntp3.begin(), ntp3.end());
    } // pfp.Vx3ID[0] > 0
    
    pfp.Tp3s = ntp3;
    // The directions were set above. Set the start and end positions. Note that the start position
    // may have been previously determined by a vertex but that is now superseded by the actual start
    // of the pfp
    pfp.XYZ[0] = ntp3[0].Pos;
    pfp.XYZ[1] = ntp3[ntp3.size()-1].Pos;
//    if(prt) PrintTp3s("FTp3o", slc, pfp, -1);

  } // FollowTp3s
  /////////////////////////////////////////
  bool FitTp3s(TCSlice& slc, const std::vector<TrajPoint3>& tp3s, Point3_t& pos, Vector3_t& dir, float& rCorr)
  {
    return FitTp3s(slc, tp3s, 0, tp3s.size(), pos, dir, rCorr);
  } // FitTp3s
  
  /////////////////////////////////////////
  bool FitTp3s(TCSlice& slc, const std::vector<TrajPoint3>& tp3s, unsigned short fromPt, unsigned short toPt, Point3_t& pos, Vector3_t& dir, float& rCorr)
  {
    // Fits the Tj2Pts points in Tp3s to a line
    if(tp3s.size() < 3) return false;
    if(fromPt >= toPt) return false;
    if(toPt > tp3s.size()) return false;

    // temp vectors to ensure that a TP is only used once
    std::vector<unsigned short> useID;
    std::vector<unsigned short> useIpt;
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3 = tp3s[ipt];
      for(auto& tj2pt : tp3.Tj2Pts) {
        bool isUsed = false;
        for(unsigned short ii = 0; ii < useID.size(); ++ii) {
          if(tj2pt.id == useID[ii] && tj2pt.ipt == useIpt[ii]) isUsed = true;
        } // ii
        if(isUsed) continue;
        // add it to the list
        useID.push_back(tj2pt.id);
        useIpt.push_back(tj2pt.ipt);
        auto& tj = slc.tjs[tj2pt.id - 1];
        ++cntInPln[DecodeCTP(tj.CTP).Plane];
      } // tj2pt
    } // ipt
    // ensure there are at least two points in at least two planes
    unsigned short enufInPlane = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] > 1) ++enufInPlane;
    if(enufInPlane < 2) return false;
    
    const unsigned int nvars = 4;
    unsigned int npts = useID.size();
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    
    // X origin
    double x0 = 0;
    for(unsigned short ipt = 0; ipt < useID.size(); ++ipt) {
      auto& tp = slc.tjs[useID[ipt] - 1].Pts[useIpt[ipt]];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      x0 += tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, planeID);
    }
    x0 /= (double)useID.size();

    double wght = 1;
    for(unsigned short ipt = 0; ipt < useID.size(); ++ipt) {
      auto& tp = slc.tjs[useID[ipt] - 1].Pts[useIpt[ipt]];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      unsigned int cstat = planeID.Cryostat;
      unsigned int tpc = planeID.TPC;
      unsigned int plane = planeID.Plane;
      // get the wire plane offset
      double off = tcc.geom->WireCoordinate(0, 0, plane, tpc, cstat);
      // get the "cosine-like" component
      double cw = tcc.geom->WireCoordinate(1, 0, plane, tpc, cstat) - off;
      // the "sine-like" component
      double sw = tcc.geom->WireCoordinate(0, 1, plane, tpc, cstat) - off;
      double x = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, planeID) - x0;
      A[ipt][0] = wght * cw;
      A[ipt][1] = wght * sw;
      A[ipt][2] = wght * cw * x;
      A[ipt][3] = wght * sw * x;
      w[ipt] = wght * (tp.Pos[0] - off);
    } // ipt
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    dir[0] = 1 / norm;
    dir[1] = tVec[2] / norm;
    dir[2] = tVec[3] / norm;
    pos[0] = x0;
    pos[1] = tVec[0];
    pos[2] = tVec[1];
    rCorr = 1;
//    std::cout<<"FTP3s: "<<useID.size()<<" cntInPln "<<cntInPln[0]<<" "<<cntInPln[1]<<" "<<cntInPln[2]<<"\n";
    return true;

  } // FitTp3s
  
  /////////////////////////////////////////
  bool FitTp3(TCSlice& slc, TrajPoint3& tp3, const std::vector<Tj2Pt>& tj2pts)
  {
    // Fits the vector of Tj2Pts points and puts the results into tp3. This code is adapted
    // from TrackLineFitAlg: SVD fit adapted from $ROOTSYS/tutorials/matrix/solveLinear.C
    // Fit equation is w = A(X)v, where w is a vector of hit wires, A is
    // a matrix to calculate a track projected to a point at X, and v is
    // a vector (Yo, Zo, dY/dX, dZ/dX).
    if(tj2pts.size() < 4) return false;

    const unsigned int nvars = 4;
    unsigned int npts = tj2pts.size();
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    
    double x0 = 0;
    for(auto& tj2pt : tj2pts) {
      auto& tp = slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      x0 += tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, planeID);
    }
    x0 /= (double)tj2pts.size();
    
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    double wght = 1;
    for(unsigned short ipt = 0; ipt < tj2pts.size(); ++ipt) {
      auto& tj2pt = tj2pts[ipt];
      auto& tp = slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      unsigned int cstat = planeID.Cryostat;
      unsigned int tpc = planeID.TPC;
      unsigned int plane = planeID.Plane;
      // get the wire plane offset
      double off = tcc.geom->WireCoordinate(0, 0, plane, tpc, cstat);
      // get the "cosine-like" component
      double cw = tcc.geom->WireCoordinate(1, 0, plane, tpc, cstat) - off;
      // the "sine-like" component
      double sw = tcc.geom->WireCoordinate(0, 1, plane, tpc, cstat) - off;
      double x = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, planeID) - x0;
      A[ipt][0] = wght * cw;
      A[ipt][1] = wght * sw;
      A[ipt][2] = wght * cw * x;
      A[ipt][3] = wght * sw * x;
      w[ipt] = wght * (tp.Pos[0] - off);
      ++ninpl[plane];
      // need at least two points in a plane
      if(ninpl[plane] == 2) ++nok;
    } // ipt

    // need at least 2 planes with at least two points
    if(nok < 2) return false;
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    
    // Calculate Chi/DOF here
//    tp3.ChiDOF = 1;
    
    Vector3_t fitDir;
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    fitDir[0] = 1 / norm;
    fitDir[1] = tVec[2] / norm;
    fitDir[2] = tVec[3] / norm;
    
    Point3_t fitPos;
    fitPos[0] = x0;
    fitPos[1] = tVec[0];
    fitPos[2] = tVec[1];
    // move it to the same Z position as tp3.Pos
    if(tp3.Pos[2] != 0) {
      double dz = tp3.Pos[2] - fitPos[2];
      fitPos[0] += dz * fitDir[0] / fitDir[2];
      fitPos[1] += dz * fitDir[1] / fitDir[2];
      fitPos[2] += dz;
    }
    
    if(PosSep2(fitPos, tp3.Pos) > 5) {
      std::cout<<"Crazy fitPos "<<PosSep(fitPos, tp3.Pos)<<"\n";
//      tp3.ChiDOF = 10;
      return false;
    }
    
    tp3.Pos = fitPos;
    tp3.Dir = fitDir;

    return true;
  } // FitTp3
  
  /////////////////////////////////////////
  void FindCompleteness(TCSlice& slc, PFPStruct& pfp, bool doFit, bool fillTp3s, bool prt)
  {
    // Calculate the 3D-matching completeness of the set of Tjs in pfp.TjIDs and store in pfp.EffPur.
    // The completeness for each Tj is put in pfp.TjCompleteness. The TP-weighted average completeness
    // is put in pfp.EffPur. This function also fits the matching points to a 3D line and puts the
    // position and direction in pfp.XYZ[0] and pfp.Dir[0]. The pfp.TP3s vector is optionally filled.
    
    if(pfp.TjIDs.size() < 2) return;
    if(tcc.match3DCuts[0] <= 0) return;
    // This function uses mallTraj but it isn't necessarily a failure if it doesn't exist
    if(slc.mallTraj.size() < 6) return;

    pfp.TjCompleteness.resize(pfp.TjIDs.size());
    std::fill(pfp.TjCompleteness.begin(), pfp.TjCompleteness.end(), 0);
    if(fillTp3s) pfp.Tp3s.clear();

    bool twoPlanes = (slc.nPlanes == 2);
    bool twoTjs = (pfp.TjIDs.size() == 2);
    // decide if special handling of small angle pfps is required when filling Tp3s
    bool smallAngle = false;
    if(fillTp3s) {
      smallAngle = (pfp.Dir[0][0] != 0 && std::abs(pfp.Dir[0][0]) < 0.1);
//      if(pfp.Dir[0][0] == 0 && slc.DebugMode) std::cout<<"P"<<pfp.ID<<" Dir[0] isn't defined\n";
    }
    double yzcut = 1.5 * tcc.match3DCuts[0];
    
    // initialize the fit sums
    Point3_t point;
    Vector3_t dir;
    if(doFit) Fit3D(0, point, dir, point, dir);
    
    // create a vector of bools for each tj for points that are matched in 3D 
    // cast the IDs into an unsigned short for faster comparing
    std::vector<unsigned short> tjids(pfp.TjIDs.size());
    // This vector is for matches in 3 planes
    std::vector<std::vector<bool>> tjptMat3;
    // This vector is for matches in 2 planes
    std::vector<std::vector<bool>> tjptMat2;
    // and the plane index
    std::vector<unsigned short> tjplane;
    // Set a maximum size for the TP3s vector
    unsigned int maxTp3Size = 10000;
    // Initialize the vectors
    for(unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
      if(pfp.TjIDs[itj] <= 0) {
        std::cout<<"FindCompleteness: Bad tjid "<<pfp.TjIDs[itj]<<"\n";
        return;
      }
      tjids[itj] = pfp.TjIDs[itj];
      auto& tj = slc.tjs[pfp.TjIDs[itj] - 1];
      std::vector<bool> tmp(tj.Pts.size(), false);
      tjptMat2.push_back(tmp);
      if(!twoPlanes) tjptMat3.push_back(tmp);
      tjplane.push_back(DecodeCTP(tj.CTP).Plane);
    } // tjid
    
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      unsigned short indx = 0;
      for(indx = 0; indx < tjids.size(); ++indx) if(iTjPt.id == tjids[indx]) break;
      // require that the Tj ID of this point be in the list
      if(indx == tjids.size()) continue;
      auto& itj = slc.tjs[iTjPt.id - 1];
//      if(itj.AlgMod[kMat3D]) continue;
      auto& itp = itj.Pts[iTjPt.ipt];
      unsigned short iplane = DecodeCTP(itp.CTP).Plane;
      unsigned short tpc = DecodeCTP(itp.CTP).TPC;
      unsigned short cstat = DecodeCTP(itp.CTP).Cryostat;
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.ctp == iTjPt.ctp) continue;
        unsigned short jndx = 0;
        for(jndx = 0; jndx < tjids.size(); ++jndx) if(jTjPt.id == tjids[jndx]) break;
        // require that the Tj ID of this point be in the list
        if(jndx == tjids.size()) continue;
        // check for x range overlap. We know that jTjPt.xlo is > iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtj = slc.tjs[jTjPt.id - 1];
//        if(jtj.AlgMod[kMat3D]) continue;
        auto& jtp = jtj.Pts[jTjPt.ipt];
        TrajPoint3 ijtp3;
        if(!MakeTp3(slc, itp, jtp, ijtp3, true)) continue;
        ijtp3.Tj2Pts.resize(2);
        ijtp3.Tj2Pts[0] = iTjPt;
        ijtp3.Tj2Pts[1] = jTjPt;
        // Set the 2-plane match bits
        tjptMat2[indx][iTjPt.ipt] = true;
        tjptMat2[jndx][jTjPt.ipt] = true;
        if(twoPlanes) {
          if(fillTp3s && pfp.Tp3s.size() < maxTp3Size) {
            bool saveIt = true;
            FindAlongTrans(pfp.XYZ[0], pfp.Dir[0], ijtp3.Pos, ijtp3.AlongTrans);
            // cut on transverse distance
            if(smallAngle) saveIt = ijtp3.AlongTrans[1] < 1;
            if(saveIt) pfp.Tp3s.push_back(ijtp3);
          }
          if(doFit) Fit3D(1, ijtp3.Pos, ijtp3.Dir, point, dir);
          continue;
        }
        // count it as a triple if this point is in a dead region
        unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
        unsigned short kplane = 3 - iplane - jplane;
        float fwire = tcc.geom->WireCoordinate(ijtp3.Pos[1], ijtp3.Pos[2], kplane, tpc, cstat);
        if(fwire < -0.4) continue;
        unsigned int kwire = std::nearbyint(fwire);
        if(kwire < slc.wireHitRange[kplane].size() && slc.wireHitRange[kplane][kwire].first == -1) {
          // accumulate the fit sums
          if(doFit) Fit3D(1, ijtp3.Pos, ijtp3.Dir, point, dir);
          // fill Tp3s?
          if(fillTp3s && pfp.Tp3s.size() < maxTp3Size) {
            bool saveIt = true;
            FindAlongTrans(pfp.XYZ[0], pfp.Dir[0], ijtp3.Pos, ijtp3.AlongTrans);
            if(smallAngle) saveIt = ijtp3.AlongTrans[1] < 1;
            if(saveIt) pfp.Tp3s.push_back(ijtp3);
          }
          continue;
        } // dead wire in kplane
        for(unsigned int kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
          auto& kTjPt = slc.mallTraj[kpt];
          // ensure that the planes are different
          if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
          // Look for this tj point in tjids
          unsigned short kndx = 0;
          for(kndx = 0; kndx < tjids.size(); ++kndx) if(kTjPt.id == tjids[kndx]) break;
          // require that the Tj ID of this point be in the list if we aren't filling the Tp3s
          if(!fillTp3s && kndx == tjids.size()) continue;
          if(kTjPt.xlo > iTjPt.xhi) continue;
          // break out if the x range difference becomes large
          if(kTjPt.xlo > iTjPt.xhi + 5) break;
          auto& ktj = slc.tjs[kTjPt.id - 1];
//          if(ktj.AlgMod[kMat3D]) continue;
          auto& ktp = ktj.Pts[kTjPt.ipt];
          TrajPoint3 iktp3;
          if(!MakeTp3(slc, itp, ktp, iktp3, true)) continue;
          if(std::abs(ijtp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
          if(std::abs(ijtp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
          // make a copy of ijtp3 -> ijktp3
          auto ijktp3 = ijtp3;
          // add the Tj2Pt to it
          ijktp3.Tj2Pts.push_back(kTjPt);
          // accumulate the fit sums
          if(doFit) Fit3D(1, iktp3.Pos, iktp3.Dir, point, dir);
          // fill Tp3s?
          if(fillTp3s && pfp.Tp3s.size() < maxTp3Size) {
            // update the charge
            ijktp3.dEdx = (2 * ijktp3.dEdx + ktp.Chg) / 3;
            bool saveIt = true;
            FindAlongTrans(pfp.XYZ[0], pfp.Dir[0], ijktp3.Pos, ijktp3.AlongTrans);
            if(smallAngle) saveIt = ijktp3.AlongTrans[1] < 1;
            if(saveIt) pfp.Tp3s.push_back(ijktp3);
          }
          // Set the 3-plane match bits
          if(kndx == tjids.size()) continue;
          tjptMat3[indx][iTjPt.ipt] = true;
          tjptMat3[jndx][jTjPt.ipt] = true;
          tjptMat3[kndx][kTjPt.ipt] = true;
        } // kpt
      } // jpt
    } // ipt
    // do the fit and put the results into the pfp
    Fit3D(2, point, dir, pfp.XYZ[0], pfp.Dir[0]);
    if(prt && doFit) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FC: P"<<pfp.ID<<" fit pos "<<std::fixed<<std::setprecision(1)<<pfp.XYZ[0][0]<<" "<<pfp.XYZ[0][1]<<" "<<pfp.XYZ[0][2];
      myprt<<" fit dir "<<std::setprecision(2)<<pfp.Dir[0][0]<<" "<<pfp.Dir[0][1]<<" "<<pfp.Dir[0][2];
      myprt<<" Note: fit pos is the average position of all Tp3s - not the start or end.";
    }
    // now count the number of tj points were matched
    // total number of points with charge in all Tjs
    float tnpwc = 0;
    // total number that are matched in 3D in 3 planes
    float tcnt3 = 0;
    // total number that are matched in 3D in 2 planes
    float tcnt2 = 0;
    for(unsigned short itj = 0; itj < tjids.size(); ++itj) {
      auto& tj = slc.tjs[tjids[itj] - 1];
      // counts for each tj
      float npwc = 0;
      float cnt2 = 0;
      float cnt3 = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg <= 0) continue;
        ++npwc;
        if(tjptMat2[itj][ipt]) ++cnt2;
        if(!twoPlanes && tjptMat3[itj][ipt]) ++cnt3;
      } // ipt
      if(twoTjs) {
        pfp.TjCompleteness[itj] = cnt2 / npwc;
      } else {
        pfp.TjCompleteness[itj] = cnt3 / npwc;
      }
      tnpwc += npwc;
      tcnt3 += cnt3;
      tcnt2 += cnt2;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"FC: P"<<pfp.ID<<" T"<<tj.ID<<" npwc "<<npwc<<" cnt2 "<<cnt2<<" cnt3 "<<cnt3<<" PDGCode "<<tj.PDGCode;
        myprt<<" MCSMom "<<tj.MCSMom<<" InShower? "<<tj.SSID;
        myprt<<" TjCompleteness "<<std::setprecision(2)<<pfp.TjCompleteness[itj];
      } // prt
    } // itj
    if(twoTjs) {
      pfp.EffPur = tcnt2 / tnpwc;
    } else {
      pfp.EffPur = tcnt3 / tnpwc;
    }

  } // FindCompleteness
  
  /////////////////////////////////////////
  void FindMissedTjsInTp3s(TCSlice& slc, PFPStruct& pfp, std::vector<int>& missTjs, std::vector<float>& missFrac)
  {
    // compare the Tjs in pfp.TjIDs with the Tjs in Tp3s and return a list of Tjs
    // in Tp3s that aren't in pfp.TjIDs
    missTjs.clear();
    missFrac.clear();
    if(pfp.TjIDs.empty() || pfp.Tp3s.empty()) return;
    
    // determine the projection of the pfp direction vector in each plane.
    // Don't try to merge if the projection is small
    std::vector<float> projInPlane(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      auto tp = MakeBareTP(slc, pfp.XYZ[0], pfp.Dir[0], inCTP);
      projInPlane[plane] = tp.Delta;
    } // plane
    
    std::vector<unsigned short> pfpTjs;
    std::vector<unsigned short> usMissTjs;
    std::vector<std::vector<bool>> misTjPtMat;
    for(auto tjid : pfp.TjIDs) pfpTjs.push_back((unsigned short)tjid);
    for(auto& tp3 : pfp.Tp3s) {
      for(auto& tj2pt : tp3.Tj2Pts) {
        if(std::find(pfpTjs.begin(), pfpTjs.end(), tj2pt.id) != pfpTjs.end()) continue;
        // Tj isn't in pfp.TjIDs. See if we have it in the missed list
        unsigned short mtjIndx = 0;
        for(mtjIndx = 0; mtjIndx < usMissTjs.size(); ++mtjIndx) if(tj2pt.id == usMissTjs[mtjIndx]) break;
        if(mtjIndx == usMissTjs.size()) {
          // not in the misTjs list. Ensure that it isn't matched
          auto& mtj = slc.tjs[tj2pt.id - 1];
          if(mtj.AlgMod[kKilled] || mtj.AlgMod[kMat3D]) continue;
          // add it to the list
          usMissTjs.push_back(tj2pt.id);
          // create the point match vector
          std::vector<bool> ptMat(mtj.Pts.size(), false);
          ptMat[tj2pt.ipt] = true;
          misTjPtMat.push_back(ptMat);
        } else {
          if(tj2pt.ipt < misTjPtMat[mtjIndx].size()) misTjPtMat[mtjIndx][tj2pt.ipt] = true;
        }
      } // tj2pt
    } // tp3
    for(unsigned short im = 0; im < usMissTjs.size(); ++im) {
      int mtjid = usMissTjs[im];
      // calculate the fraction of points that are in Tp3s
      float cnt = 0;
      float mat = 0;
      auto& mtj = slc.tjs[mtjid - 1];
      // ignore if there is a high-score vertex between the missed tj and those in the pfp list
      if(SharesHighScoreVx(slc, pfp, mtj)) continue;
      for(unsigned short ipt = mtj.EndPt[0]; ipt <= mtj.EndPt[1]; ++ipt) {
        auto& mtp = mtj.Pts[ipt];
        if(mtp.Chg <= 0) continue;
        ++cnt;
        if(misTjPtMat[im][ipt]) ++mat;
      } // ipt
      float frac = mat / cnt;
      // ignore if low fraction matched
      if(frac < 0.1) continue;
      // ignore if this would only extend the tj in this plane by a small amount
      float lenInPlane = 0;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.CTP != mtj.CTP) continue;
        float len = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
        if(len > lenInPlane) lenInPlane = len;
      } // tjid
      if(cnt < 0.05 * lenInPlane) continue;
      // check the direction vector projection in this plane
      if(projInPlane[DecodeCTP(mtj.CTP).Plane] < 0.1) continue;
      missTjs.push_back(mtjid);
      missFrac.push_back(frac);
    } // im
  } // FindMissedTjsInTp3s
  
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
  
  /////////////////////////////////////////
  void Fit3D(unsigned short mode, Point3_t point, Vector3_t dir, Point3_t& fitPos, Vector3_t& fitDir)
  {
    // initialize, accumulate and fit the points. The code to fit the direction using the positions
    // of the points is commented out and replaced with a simple average of the directions of the points
    
    // 3D fit sum variables
    static double fSum, fSumx, fSumy, fSumz;
//    static double fSumx2, fSumy2, fSumz2, fSumxz, fSumyz;
    // average the direction vectors
    static double fSumDir0, fSumDir1, fSumDir2;

    if(mode == 0) {
      fSum = 0; fSumx = 0; fSumy = 0; fSumz = 0; 
//      fSumx2 = 0; fSumy2 = 0; fSumz2 = 0; fSumxz = 0; fSumyz = 0;
      fSumDir0 = 0; fSumDir1 = 0; fSumDir2 = 0;
      return;
    }
    // accumulate
    if(mode == 1) {
      fSum += 1;
      fSumx += point[0];
      fSumy += point[1];
      fSumz += point[2];
/*
      fSumx2 += point[0] * point[0];
      fSumy2 += point[1] * point[1];
      fSumz2 += point[2] * point[2];
      fSumxz += point[0] * point[2];
      fSumyz += point[1] * point[2];
*/
      fSumDir0 += dir[0];
      fSumDir1 += dir[1];
      fSumDir2 += dir[2];
      return;
    }
    
    if(fSum < 2) return;
    // just use the average for the position
    fitPos[0] = fSumx / fSum;
    fitPos[1] = fSumy / fSum;
    fitPos[2] = fSumz / fSum;
    // and for the direction
    fitDir = {{fSumDir0, fSumDir1, fSumDir2}};
    SetMag(fitDir, 1);
  } // Fit3D

  /////////////////////////////////////////
  unsigned short WiresSkippedInCTP(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP)
  {
    // counts the number of wires between the end points of all Tjs in the list of tjids
    // in inCTP where there is no TP with charge
    if(tjids.empty()) return 0;
    
    // find the min and max Pos[0] positions
    float fLoWire = 1E6;
    float fHiWire = -1E6;
    for(auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.CTP != inCTP) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        float endWire = tj.Pts[tj.EndPt[end]].Pos[0];
        if(endWire < -0.4) continue;
        if(endWire < fLoWire) fLoWire = endWire;
        if(endWire > fHiWire) fHiWire = endWire;
      } // end
    } // tjid
    if(fLoWire >= fHiWire) return 0;
    unsigned int loWire = std::nearbyint(fLoWire);
    unsigned short nWires = std::nearbyint(fHiWire) - loWire + 1;
    std::vector<bool> ptOnWire(nWires, false);
    
    // count the number of points with charge on all Tjs
    float npwc = 0;
    for(auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.CTP != inCTP) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(tp.Pos[0] < -0.4) continue;
        ++npwc;
        unsigned short indx = std::nearbyint(tp.Pos[0]) - loWire;
        if(indx < nWires) ptOnWire[indx] = true;
      } // ipt
    } // tjid
    if(npwc == 0) return 0;
    float nskip = 0;
    for(unsigned short indx = 0; indx < nWires; ++indx) if(!ptOnWire[indx]) ++nskip;
    return (nskip / npwc);
    
  } // WiresSkippedInCTP
  
  /////////////////////////////////////////
  float LengthInCTP(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP)
  {
    // Calculates the maximum length between the end points of Tjs in the list of tjids in inCTP
    if(tjids.empty()) return 0;
    // put the end point positions into a vector
    std::vector<Point2_t> endPos;
    for(auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.CTP != inCTP) continue;
      endPos.push_back(tj.Pts[tj.EndPt[0]].Pos);
      endPos.push_back(tj.Pts[tj.EndPt[1]].Pos);
    } // tjid
    if(endPos.size() < 2) return 0;
    float extent = 0;
    for(unsigned short pt1 = 0; pt1 < endPos.size() - 1; ++pt1) {
      for(unsigned short pt2 = pt1 + 1; pt2 < endPos.size(); ++pt2) {
        float sep = PosSep2(endPos[pt1], endPos[pt2]);
        if(sep > extent) extent = sep;
      } // pt2
    } // pt1
    return sqrt(extent);
  } // LengthInCTP

  /////////////////////////////////////////
  bool AddMissedTj(TCSlice& slc, PFPStruct& pfp, unsigned short itj, bool looseCuts, bool prt)
  {
    // The Tj pfp.TjIDs[itj] has poor completeness. Search matchVec for
    // the occurrence of this tj with a large completeness AND the occurrence 
    // of another tj in pfp.TjIDs.
    if(itj > pfp.TjIDs.size() - 1) return false;
    if(slc.matchVec.empty()) return false;
    
    int theTj = pfp.TjIDs[itj];
//    bool pfpInShower = (pfp.PDGCode == 11);
    
    unsigned short oldSize = pfp.TjIDs.size();
    
    for(unsigned int ims = 0; ims < slc.matchVec.size(); ++ims) {
      auto& ms = slc.matchVec[ims];
      // look for theTj in the match struct
      unsigned short tjIndex = 0;
      for(tjIndex = 0; tjIndex < ms.TjIDs.size(); ++tjIndex) if(ms.TjIDs[tjIndex] == theTj) break;
      if(tjIndex == ms.TjIDs.size()) continue;
      auto shared = SetIntersection(pfp.TjIDs, ms.TjIDs);
      if(shared.empty()) continue;
      if(looseCuts) {
        // Look for shared size at least 2 (theTj and another tj) or size 1 and higher TjCompleteness
        bool isWorse = (ms.TjCompleteness[tjIndex] < pfp.TjCompleteness[itj]);
        if(shared.size() < 2 && isWorse) continue;
      } else {
        // Look for shared size at least 2 (theTj and another tj in pfp.TjIDs)
        if(shared.size() < 2) continue;
      }
      // Add the tjs that are not in pfp.TjIDs
      for(auto tjid : ms.TjIDs) {
        if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjid) != pfp.TjIDs.end()) continue;
        pfp.TjIDs.push_back(tjid);
        // check vertex - tj consistency
        if(PFPVxTjOK(slc, pfp, prt)) continue;
        pfp.TjCompleteness.push_back(0);
        if(prt) mf::LogVerbatim("TC")<<"AMT: P"<<pfp.ID<<" T"<<theTj<<" Add T"<<tjid;
      } // mtjid
    } // ims
    if(pfp.TjIDs.size() > oldSize) return true;
    return false;
  } // AddMissedTj

  /////////////////////////////////////////
  bool MergePFPTjs(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Checks the list of Tjs in pfp.TjIDs and merges those that are in 
    // the same plane. This function uses the ordering of Tps which should
    // have been sorted
    if(pfp.TjIDs.empty()) return false;
    if(pfp.Tp3s.empty()) return false;
    
    geo::TPCGeo const& TPC = tcc.geom->TPC(pfp.TPCID);
    unsigned short nplanes = TPC.Nplanes();
    
    // see if anything needs to be done
    std::vector<unsigned short> cntInPln(nplanes);
    bool itsOK = true;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      ++cntInPln[plane];
      if(cntInPln[plane] > 1) itsOK = false;
    }
    if(itsOK) return true;
    
    // vector of tj IDs that will replace pfp.TjIDs
    std::vector<int> newTjIDs;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MergePFPTjs: P"<<pfp.ID<<" in";
      for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
    }
    
    for(unsigned short plane = 0; plane < nplanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      // save the TjIDs as unsigned short to match with Tj2Pts
      std::vector<unsigned short> tjids;
      for(auto tjid : pfp.TjIDs) if(slc.tjs[tjid - 1].CTP == inCTP) tjids.push_back((unsigned short)tjid);
      // Only one tj in this plane. No need to merge
      if(tjids.size() == 1) {
        newTjIDs.push_back((int)tjids[0]);
        continue;
      }
      // no tjs in this plane
      if(tjids.size() == 0) continue;
      // find the first ID and ipt of Tjs in this plane as they are
      // encountered while iterating through Tp3s. This scheme assumes that the Tp3s have
      // been sorted by distance from the start and the Tjs are broken end-to-end. This 
      // information will be used to determine if Tjs need to be reversed before inserting
      // the points in to the merged trajectory
      //                     Tj ID   first ipt
      std::vector<std::array<unsigned short, 2>> firstPts;
      for(unsigned short itp3 = 0; itp3 < pfp.Tp3s.size(); ++itp3) {
        auto& tp3 = pfp.Tp3s[itp3];
        for(auto& tj2pt : tp3.Tj2Pts) {
          unsigned short tjIndx = 0;
          for(tjIndx = 0; tjIndx < tjids.size(); ++tjIndx) if(tj2pt.id == tjids[tjIndx]) break;
          if(tjIndx == tjids.size()) continue;
          // look for this tj in firstPts
          unsigned short firstPtsIndx = 0;
          for(firstPtsIndx = 0; firstPtsIndx < firstPts.size(); ++firstPtsIndx) if(tj2pt.id == firstPts[firstPtsIndx][0]) break;
          if(firstPtsIndx == firstPts.size()) {
            // not found so add it
            std::array<unsigned short, 2> firstPt {{tj2pt.id, tj2pt.ipt}};
            firstPts.push_back(firstPt);
          }
        } // tj2pt
      } // itp3
      if(firstPts.empty()) continue;
      // create a new merged trajectory
      Trajectory mtj;
      // give it a bogus ID
      mtj.ID = -6666;
      mtj.CTP = inCTP;
      mtj.StepDir = 1;
      bool first = true;
      for(auto firstPt : firstPts) {
        // make a copy so we can reverse it and drop it if the merge fails
        auto tj = slc.tjs[firstPt[0] - 1];
        unsigned short midPt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
        if(firstPt[1] > midPt) ReverseTraj(slc, tj);
        // Transfer vertices to mtj
        if(first) {
          first = false;
          mtj.VtxID[0] = tj.VtxID[0];
        }
        mtj.VtxID[1] = tj.VtxID[1];
        // insert the points at the end
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          auto& tp = tj.Pts[ipt];
          if(tp.Chg <= 0) continue;
          mtj.Pts.push_back(tp);
        }
      } // firstPt
      mtj.AlgMod[kMat3DMerge] = true;
      SetEndPoints(mtj);
      mtj.MCSMom = MCSMom(slc, mtj);
      SetPDGCode(slc, mtj);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" P"<<pfp.ID<<" try to merge";
        for(auto tjid : tjids) {
          auto& tj = slc.tjs[tjid - 1];
          myprt<<" T"<<tjid<<" MCSMom "<<tj.MCSMom;
        }
        myprt<<" -> T"<<slc.tjs.size() + 1;
        myprt<<" MCSMom "<<mtj.MCSMom;
      }
      // kill the broken tjs and update the pfp TP3s
      for(auto tjid : tjids) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.SSID > 0) mtj.SSID = tj.SSID;
        MakeTrajectoryObsolete(slc, tjid - 1);
      }
      // save the new one
      if(!StoreTraj(slc, mtj)) {
        std::cout<<"MergePFPTjs: StoreTraj failed P"<<pfp.ID<<"\n";
        return false;
      }
      int newTjID = slc.tjs.size();
      newTjIDs.push_back(newTjID);
      // prepare to clobber vertices
      std::vector<unsigned short> vxlist;
      for(auto tjid : tjids) {
        // update the stored match struct and Tp3s
        UpdateMatchStructs(slc, tjid, newTjID);
        // Update the Tp3s of this pfp
        UpdateTp3s(slc, pfp, tjid, newTjID);
        auto& tj = slc.tjs[tjid - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == 0) continue;
          if(std::find(vxlist.begin(), vxlist.end(), tj.VtxID[end]) != vxlist.end()) {
            auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
//            std::cout<<"P"<<pfp.ID<<" Clobber 2V"<<vx2.ID<<"\n";
            MakeVertexObsolete("MPTJ", slc, vx2, true);
          } else {
            vxlist.push_back(tj.VtxID[end]);
          }
        } // end
      } // tjid
    } // plane
    
    pfp.TjIDs = newTjIDs;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MergePFPTjs: P"<<pfp.ID<<" out";
      for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
      PrintPFP("MPTJ", slc, pfp, true);
    }
    return true;
  } // MergePFPTjs

  /////////////////////////////////////////
  void FindXMatches(TCSlice& slc, unsigned short numPlanes, short maxScore, std::vector<MatchStruct>& matVec, bool prt)
  {
    // This function matches trajectory points in slc.mallTraj using the X position. These points should
    // have already been sorted by increasing X by the function that created mallTraj. 
    
    if(slc.mallTraj.empty()) return;
    if(tcc.match3DCuts[0] <= 0) return;
    if(numPlanes < 2) return;
    
    int cstat = DecodeCTP(slc.mallTraj[0].ctp).Cryostat;
    int tpc = DecodeCTP(slc.mallTraj[0].ctp).TPC;
    constexpr float twopi = 2 * M_PI;
    constexpr float piOver2 = M_PI / 2;
    
    // create a temp vector to check for duplicates
    auto inMatVec = matVec;
    std::vector<MatchStruct> temp;
    
    // the minimum number of points for matching
    unsigned short minPts = 2;
    // override this with the user minimum for 2-plane matches
    if(numPlanes == 2) minPts = tcc.match3DCuts[2];
    
    // max number of match combos left
    unsigned int nAvailable = 0;
    if(matVec.size() < tcc.match3DCuts[4]) nAvailable = tcc.match3DCuts[4] - matVec.size();
    if(nAvailable == 0 || nAvailable > tcc.match3DCuts[4]) return;
    
    // these cuts presume that the resolution in X is better than it is in Y and Z
    float xcut = tcc.match3DCuts[0];
    double yzcut = 1.5 * xcut;
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // length cut
      if(iTjPt.npts < minPts) continue;
      // look for matches using Tjs that have the correct score
      if(iTjPt.score < 0 || iTjPt.score > maxScore) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      unsigned short iplane = DecodeCTP(itp.CTP).Plane;
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.ctp == iTjPt.ctp) continue;
        // length cut
        if(jTjPt.npts < minPts) continue;
        // score cut
        if(jTjPt.score < 0 || jTjPt.score > maxScore) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
        TrajPoint3 tp3;
        if(!MakeTp3(slc, itp, jtp, tp3, false)) continue;
        // count weight is one for a two-plane match
        float cntWght = 1;
        if(numPlanes == 3) {
          // numPlanes == 3
          for(unsigned int kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
            auto& kTjPt = slc.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
            if(kTjPt.score < 0 || kTjPt.score > maxScore) continue;
            if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 5) break;
            auto& ktp = slc.tjs[kTjPt.id - 1].Pts[kTjPt.ipt];
            unsigned short kplane = DecodeCTP(ktp.CTP).Plane;
            TrajPoint3 iktp3;
            if(!MakeTp3(slc, itp, ktp, iktp3, false)) continue;
            if(std::abs(tp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
            if(std::abs(tp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
            float dang = 0;
            if(tcc.match3DCuts[1] > 0) {
              dang = std::abs(DeltaAngle(tp3.Dir, iktp3.Dir));
              while(dang >  M_PI) dang -= twopi;
              if(dang >  piOver2) dang = M_PI - dang;
              float mcsmom = slc.tjs[iTjPt.id - 1].MCSMom + slc.tjs[jTjPt.id - 1].MCSMom + slc.tjs[kTjPt.id - 1].MCSMom;
              mcsmom /= 3;
              if(mcsmom > 150 && dang > tcc.match3DCuts[1]) continue;
            }
            // we have a match.
            // Just fill temp. See if the Tj IDs are in the match list.
            // first check the input matVec
            bool gotit = false;
            for(auto& ms : inMatVec) {
              if(ms.TjIDs.size() != 3) continue;
              if(iTjPt.id == ms.TjIDs[iplane] && jTjPt.id == ms.TjIDs[jplane] && kTjPt.id == ms.TjIDs[kplane]) {
                gotit = true;
                break;
              } 
            } // ms
            if(gotit) continue;
            // Triple match count = 2 de-weighted by delta angle
            cntWght = 2 - dang;
            if(cntWght <= 0) continue;
            // next check the temp vector
            unsigned short indx = 0;
            for(indx = 0; indx < temp.size(); ++indx) {
              auto& ms = temp[indx];
              if(iTjPt.id != ms.TjIDs[iplane]) continue;
              if(jTjPt.id != ms.TjIDs[jplane]) continue;
              if(kTjPt.id != ms.TjIDs[kplane]) continue;
              ms.Count += cntWght;
              break;
            } // indx
            if(indx == temp.size()) {
              // not found in the match vector so add it
              MatchStruct ms;
              ms.TjIDs.resize(3);
              // Note that we can put the Tj IDs in plane-order since there are 3 of them
              // This is not the case when there are 2 planes
              ms.TjIDs[iplane] = iTjPt.id;
              ms.TjIDs[jplane] = jTjPt.id;
              ms.TjIDs[kplane] = kTjPt.id;
              ms.Count = cntWght;
              temp.push_back(ms);
              // give up if there are too many
              if(temp.size() > nAvailable) break;
            } // not found in the list
          } // kpt
          // numPlanes == 3
        } else {
          // 2-plane TPC or 2-plane match in a 3-plane TPC
          if(slc.nPlanes == 3) {
            // See if there is a signal at this point.
            unsigned short kplane = 3 - iplane - jplane;
            float fkwire = tcc.geom->WireCoordinate(tp3.Pos[1], tp3.Pos[2], kplane, tpc, cstat);
            if(fkwire < 0 || fkwire > tcc.maxPos0[kplane]) continue;
            TrajPoint tpk;
            tpk.CTP = EncodeCTP(cstat, tpc, kplane);
            tpk.Pos[0] = fkwire;
            float xp = 0.5 * (iTjPt.xlo + iTjPt.xhi);
            tpk.Pos[1] = tcc.detprop->ConvertXToTicks(xp, kplane, tpc, cstat) * tcc.unitsPerTick;
            // Note that SignalAtTp assumes that a signal exists if the wire is dead
            if(!SignalAtTp(slc, tpk)) continue;
          }
          // Just fill temp. See if the Tj IDs are in the match list
          bool gotit = false;
          for(auto& ms : inMatVec) {
            if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), iTjPt.id) != ms.TjIDs.end() && 
               std::find(ms.TjIDs.begin(), ms.TjIDs.end(), jTjPt.id) != ms.TjIDs.end()) {
              gotit = true;
              break;
            } 
          } // ms
          if(gotit) continue;
          unsigned short indx = 0;
          for(indx = 0; indx < temp.size(); ++indx) {
            auto& ms = temp[indx];
            if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), iTjPt.id) != ms.TjIDs.end() &&
               std::find(ms.TjIDs.begin(), ms.TjIDs.end(), jTjPt.id) != ms.TjIDs.end()) break;
          } // indx
          if(indx == temp.size()) {
            MatchStruct ms;
            ms.TjIDs.resize(2);
            // Here we put the Tj IDs in no particular order
            ms.TjIDs[0] = iTjPt.id;
            ms.TjIDs[1] = jTjPt.id;
            ms.Count = 1;
            temp.push_back(ms);
          } // not found in the list
          else {
            ++temp[indx].Count;
          }
        } // 2-plane TPC
        // give up if there are too many
        if(temp.size() > nAvailable) break;
      } // jpt
      // give up if there are too many
      if(temp.size() > nAvailable) break;
    } // ipt
    
    // temp
    
    if(temp.empty()) return;

    if(temp.size() == 1) {
      matVec.push_back(temp[0]);
    } else {
      // multiple entries - need to sort by decreasing match count
      std::vector<SortEntry> sortVec(temp.size());
      for(unsigned int indx = 0; indx < sortVec.size(); ++indx) {
        sortVec[indx].index = indx;
        sortVec[indx].val = temp[indx].Count;
      } // indx
      std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
      // Re-order temp
      auto tmpVec = temp;
      for(unsigned int ii = 0; ii < sortVec.size(); ++ii) temp[ii] = tmpVec[sortVec[ii].index];
      // insert it after the triple matches
      matVec.insert(matVec.end(), temp.begin(), temp.end());
    } // temp size > 1
    
    if(prt) mf::LogVerbatim("TC")<<"FindXMatches: Found "<<temp.size()<<" matches";
    
  } // FindXMatches

  /////////////////////////////////////////
  bool MakeTp3(TCSlice& slc, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3, bool findDirection)
  {
    // Make a 3D trajectory point using two 2D trajectory points
    tp3.Dir = {{999.0, 999.0, 999.0}};
    tp3.Pos = {{999.0, 999.0, 999.0}};
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);
    double upt = tcc.unitsPerTick;
    double ix = tcc.detprop->ConvertTicksToX(itp.Pos[1] / upt, iPlnID);
    double jx = tcc.detprop->ConvertTicksToX(jtp.Pos[1] / upt, jPlnID);
    
    // don't continue if the points are wildly far apart in X
    if(std::abs(ix - jx) > 20) return false;
    tp3.Pos[0] = 0.5 * (ix + jx);
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
    if(den == 0) return false;
    double iPos0 = itp.Pos[0];
    double jPos0 = jtp.Pos[0];
    // Find the Z position of the intersection
    tp3.Pos[2] = (jcs * (iPos0 - iw0) - ics * (jPos0 - jw0)) / den;
    // and the Y position
    bool useI = std::abs(ics) > std::abs(jcs);
    if(useI) {
      tp3.Pos[1] = (iPos0 - iw0 - isn * tp3.Pos[2]) / ics;
    } else {
      tp3.Pos[1] = (jPos0 - jw0 - jsn * tp3.Pos[2]) / jcs;
    }
    
    // Now find the direction. Protect against large angles first
    if(jtp.Dir[1] == 0) {
      // Going either in the +X direction or -X direction
      if(jtp.Dir[0] > 0) { tp3.Dir[0] = 1; } else { tp3.Dir[0] = -1; }
      tp3.Dir[1] = 0;
      tp3.Dir[2] = 0;
      return true;
    } // jtp.Dir[1] == 0
    
    // stuff the average TP charge into dEdx
    tp3.dEdx = 0.5 * (itp.Chg + jtp.Chg);
    
    if(!findDirection) return true;

    // make a copy of itp and shift it by many wires to avoid precision problems
    double itp2_0 = itp.Pos[0] + 100;
    double itp2_1 = itp.Pos[1];
    if(std::abs(itp.Dir[0]) > 0.01) itp2_1 += 100 * itp.Dir[1] / itp.Dir[0];
    // Create a second Point3 for the shifted point
    Point3_t pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = tcc.detprop->ConvertTicksToX(itp2_1 / upt, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = tcc.detprop->ConvertXToTicks(pos2[0], jPlnID) * upt;
    // Find the wire position (Pos0) in the j plane that this corresponds to
    double jtp2Pos0 = (jtp2Pos1 - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    // Find the Y,Z position using itp2 and jtp2Pos0
    pos2[2] = (jcs * (itp2_0 - iw0) - ics * (jtp2Pos0 - jw0)) / den;
    if(useI) {
      pos2[1] = (itp2_0 - iw0 - isn * pos2[2]) / ics;
    } else {
      pos2[1] = (jtp2Pos0 - jw0 - jsn * pos2[2]) / jcs;
    }
    double sep = PosSep(tp3.Pos, pos2);
    if(sep == 0) return false;
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) tp3.Dir[ixyz] = (pos2[ixyz] - tp3.Pos[ixyz]) /sep;

    return true;
    
  } // MakeTP3
  
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
    // Fills the dEdX vector in the match struct. This function should be called after the
    // matched trajectory points are ordered so that dE/dx is calculated at the start of the PFParticle
    if(pfp.ID == 0) return;
    // error check
    bool notgood = false;
    for(unsigned short end = 0; end < 2; ++end) {
      if(pfp.dEdx[end].size() != slc.nPlanes) notgood = true;
      if(pfp.dEdxErr[end].size() != slc.nPlanes) notgood = true;
    }
    if(notgood) {
      //      if(prt) mf::LogVerbatim("TC")<<"FilldEdx found inconsistent sizes for dEdx\n";
      return;
    }
    
    double t0 = 0;
    
    // don't attempt to find dE/dx at the end of a shower
    unsigned short numEnds = 2;
    if(pfp.PDGCode == 1111) numEnds = 1;
    
    unsigned short maxlen = 0;
    for(auto tjID : pfp.TjIDs) {
      
      Trajectory& tj = slc.tjs[tjID - 1];
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      double angleToVert = tcc.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
      for(unsigned short end = 0; end < numEnds; ++end) {
        pfp.dEdx[end][planeID.Plane] = 0;
        tj.dEdx[end] = 0;
        double cosgamma = std::abs(std::sin(angleToVert) * pfp.Dir[end][1] + std::cos(angleToVert) * pfp.Dir[end][2]);
        if(cosgamma == 0) continue;
        double dx = tcc.geom->WirePitch(planeID) / cosgamma;
        if(dx == 0) continue;
        double dQ = tj.Pts[tj.EndPt[end]].AveChg;
        if(dQ == 0) continue;
        // convert to dQ/dx
        dQ /= dx;
        double time = tj.Pts[tj.EndPt[end]].Pos[1] / tcc.unitsPerTick;
        float dedx = tcc.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
        if(dedx > 999) dedx = -1;
        pfp.dEdx[end][planeID.Plane] = dedx;
        tj.dEdx[end] = dedx;
        // ChgRMS is the fractional error
        pfp.dEdxErr[end][planeID.Plane] = dedx * tj.ChgRMS;
        
      } // end
      // Grab the best plane iusing the start f 1 < dE/dx < 50 MeV/cm
      if(pfp.dEdx[0][planeID.Plane] > 1 && pfp.dEdx[0][planeID.Plane] < 50) {
        if(tj.Pts.size() > maxlen) {
          maxlen = tj.Pts.size();
          pfp.BestPlane = planeID.Plane;
        }
      } // valid dE/dx
      
    } // tj
  } // FilldEdX

  ////////////////////////////////////////////////
  void FilldEdx(TCSlice& slc, TrajPoint3& tp3)
  {
    // fills the dEdx variables in tp3 (MeV/cm)
    tp3.dEdx = 0;
    if(tp3.Tj2Pts.empty()) return;
    double t0 = 0;
    float sum = 0;
    float sum2 = 0;
    float cnt = 0;
    for(auto& tj2pt : tp3.Tj2Pts) {
      auto& tp = slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt];
      if(tp.Chg <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      double angleToVert = tcc.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
      double cosgamma = std::abs(std::sin(angleToVert) * tp3.Dir[1] + std::cos(angleToVert) * tp3.Dir[2]);
      if(cosgamma == 0) continue;
      double dx = tcc.geom->WirePitch(planeID) / cosgamma;
      // sum the hit charge
      double chg = 0;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        chg += (*evt.allHits)[slc.slHits[iht].allHitsIndex].Integral();
      }
      double dQ = chg / dx;
      double time = tp.Pos[1] / tcc.unitsPerTick;
      float dedx = tcc.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
      if(dedx > 200) continue;
      sum += dedx;
      sum2 += dedx * dedx;
      ++cnt;
    } // tj2pt
    if(cnt == 0) return;
    tp3.dEdx = sum / cnt;
/*
    if(cnt > 1) {
      float arg = sum2 - cnt * tp3.dEdx * tp3.dEdx;
      if(arg > 0) {
        tp3.dEdxErr = sqrt(arg / (cnt - 1));
        // convert to error on the mean
        tp3.dEdxErr /= sqrt(cnt);
      }
    } // cnt > 1
*/
  } // FilldEdx

  ////////////////////////////////////////////////
  float PFPDOCA(const PFPStruct& pfp1,  const PFPStruct& pfp2, unsigned short& close1, unsigned short& close2)
  {
    // returns the Distance of Closest Approach between two PFParticles. 
    close1 = USHRT_MAX;
    close2 = USHRT_MAX;
    float minSep2 = 1E8;
    for(unsigned short ipt1 = 0; ipt1 < pfp1.Tp3s.size(); ++ipt1) {
      auto& tp1 = pfp1.Tp3s[ipt1];
      for(unsigned short ipt2 = 0; ipt2 < pfp2.Tp3s.size(); ++ipt2) {
        auto& tp2 = pfp2.Tp3s[ipt2];
        float sep2 = PosSep2(tp1.Pos, tp2.Pos);
        if(sep2 > minSep2) continue;
        minSep2 = sep2;
        close1 = ipt1;
        close2 = ipt2;
      } // tp2
    } // tp1
    return sqrt(minSep2);
 } // PFPDOCA

  ////////////////////////////////////////////////
  bool Split3DKink(TCSlice& slc, PFPStruct& pfp, double sep, bool prt)
  {
    // Finds kinks in the PFParticle, splits slc, creates 2D vertices and forces a rebuild if any are found
    if(pfp.Tp3s.empty()) return false;
    if(!tcc.useAlg[kSplit3DKink]) return false;
    
    auto kinkPts = FindKinks(slc, pfp, sep, prt);
    if(kinkPts.empty()) return false;
    if(prt) mf::LogVerbatim("TC")<<"Split3DKink found a kink at Tp3s point "<<kinkPts[0];
    
    // Only split the biggest angle kink
    double big = 0;
    unsigned short kpt = 0;
    for(auto ipt : kinkPts) {
      double dang = KinkAngle(slc, pfp.Tp3s, ipt, sep);
      if(dang > big) {
        big = dang;
        kpt = ipt;
      }
    } // ipt
    if(kpt < 1 || kpt > pfp.Tp3s.size() - 1) return false;
    // determine which tjs need to be split
    std::vector<unsigned short> tjids;
    std::vector<unsigned short> vx2ids;
    // inspect a few Tp3s near the kink point to get a list of Tjs
    for(unsigned short ipt = kpt; ipt < kpt + 2; ++ipt) {
      auto& tp3 = pfp.Tp3s[ipt];
      for(auto& tp2 : tp3.Tj2Pts) {
        // see if this Tj id is in the list
        if(std::find(tjids.begin(), tjids.end(), tp2.id) != tjids.end()) continue;
        // ensure that it is pfp.TjIDs
        if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tp2.id) == pfp.TjIDs.end()) continue;
        tjids.push_back(tp2.id);
        auto& tj = slc.tjs[tp2.id - 1];
        auto& tp = tj.Pts[tp2.ipt];
        unsigned short closeEnd = USHRT_MAX;
        if(tp2.ipt < tj.EndPt[0] + 2) closeEnd = 0;
        if(tp2.ipt > tj.EndPt[1] - 2) closeEnd = 1;
        if(closeEnd < 2) {
          // No split is needed and there should be a vertex at this end of the Tj that
          // should be associated with a 3D vertex that we will construct
          if(tj.VtxID[closeEnd] == 0) {
//            std::cout<<Split3DKink: TODO Tj "<<tj.ID<<" has no vertex attached on end "<<closeEnd<<". Write some code.\n";
            return false;
          }
          vx2ids.push_back(tj.VtxID[closeEnd]);
          if(prt) mf::LogVerbatim("TC")<<" tj "<<tj.ID<<" use existing 2V"<<tj.VtxID[closeEnd];
        } else {
          // make a 2D vertex at this point
          VtxStore vx2;
          vx2.ID = slc.vtxs.size() + 1;
          vx2.CTP = tj.CTP;
          vx2.Topo = 10;
          vx2.Pos = tp.Pos;
          if(!StoreVertex(slc, vx2)) return false;
          if(!SplitTraj(slc, tp2.id - 1, tp2.ipt, slc.vtxs.size() - 1, prt)) return false;
          vx2ids.push_back(vx2.ID);
          AttachAnyTrajToVertex(slc, slc.vtxs.size() - 1, prt);
          if(prt) mf::LogVerbatim("TC")<<" tj "<<tj.ID<<" new 2V"<<vx2.ID;
        }
      } // tp2
    } // ipt
    
    if(vx2ids.size() != slc.nPlanes) {
//      std::cout<<"Split3DKink: TODO pfp "<<pfp.ID<<" only has "<<vx2ids.size()<<" 2D vertices. \n";
      return false;
    }
    Vtx3Store vx3;
    vx3.TPCID = pfp.TPCID;
    vx3.Vx2ID.resize(slc.nPlanes);
    vx3.ID = slc.vtx3s.size() + 1;
    vx3.X = pfp.Tp3s[kpt].Pos[0];
    vx3.Y = pfp.Tp3s[kpt].Pos[1];
    vx3.Z = pfp.Tp3s[kpt].Pos[2];
    for(auto vx2id : vx2ids) {
      if(vx2id == 0) continue;
      auto& vx2 = slc.vtxs[vx2id - 1];
      unsigned short plane = DecodeCTP(vx2.CTP).Plane;
      vx3.Vx2ID[plane] = vx2id;
      vx2.Vx3ID = vx3.ID;
    } // vx2id
    ++evt.globalS3ID;
    vx3.UID = evt.globalS3ID;
    slc.vtx3s.push_back(vx3);
    // mark this as needing an update
    pfp.NeedsUpdate = true;
    return true;
  } // Split3DKink
  
  ////////////////////////////////////////////////
  std::vector<unsigned short> FindKinks(TCSlice& slc, PFPStruct& pfp, double sep, bool prt)
  {
    // returns a vector of indices in pfp.Tp3s where kinks exist. The kink angle is calculated using
    // Tp3s separated by +/- sep (cm)
    std::vector<unsigned short> kinkPts;
    // look for a kink angle greater than angCut
    double angCut = 0.3;
    double kang = 0;
    unsigned short kStart = USHRT_MAX;
    // foundKink is set true after a kink is found to skip past some number of points after the kink
    bool foundKink = false;
    double kinkSep2 = 2 * sep * sep;
    for(unsigned short ipt = 1; ipt < pfp.Tp3s.size(); ++ipt) {
      // skip ahead after a kink?
      if(foundKink) {
        // location of the previously found kink
        unsigned short kpt = kinkPts[kinkPts.size() - 1];
        if(foundKink && PosSep2(pfp.Tp3s[ipt].Pos, pfp.Tp3s[kpt].Pos) < kinkSep2) continue;
      }
      foundKink = false;
      double dang = KinkAngle(slc, pfp.Tp3s, ipt, sep);
      if(dang < angCut && kStart == USHRT_MAX) continue;
      // found a kink larger than the cut. See if this is the onset of a kink
      if(kStart == USHRT_MAX) {
        // onset of a kink
        kStart = ipt;
      } else {
        // a kink was found. Keep scanning until delta angle (dang) is less than the maximum (kang)
        if(dang < kang) {
          unsigned short klen = ipt - kStart;
          if(prt) mf::LogVerbatim("TC")<<" findKinks: kink angle "<<kang<<" at point "<<ipt<<" klen "<<klen;
          kinkPts.push_back(ipt - 1);
          foundKink = true;
          kStart = USHRT_MAX;
        } // dang < kang
      } // kink found
      kang = dang;
    } // ipt
    return kinkPts;
  } // findKinks

  ////////////////////////////////////////////////
  double KinkAngle(TCSlice& slc, const std::vector<TrajPoint3>& tp3s, unsigned short atPt, double sep)
  {
    // calculate a kink angle at the TjPt 
    if(tp3s.empty()) return -1;
    if(atPt < 1 || atPt > tp3s.size() - 2) return -1;
    double sep2 = sep * sep;
    unsigned short pt1 = USHRT_MAX;
    for(unsigned short ii = 1; ii < tp3s.size(); ++ii) {
      unsigned short ipt = atPt - ii;
      if(PosSep2(tp3s[atPt].Pos, tp3s[ipt].Pos) > sep2) {
        pt1 = ipt;
        break;
      }
      if(ipt == 0) break;
    } // ii
    if(pt1 == USHRT_MAX) return -1;
    unsigned short pt2 = USHRT_MAX;
    for(unsigned short ii = 1; ii < tp3s.size(); ++ii) {
      unsigned short ipt = atPt + ii;
      if(ipt == tp3s.size()) break;
      if(PosSep2(tp3s[atPt].Pos, tp3s[ipt].Pos) > sep2) {
        pt2 = ipt;
        break;
      }
    } // ii
    if(pt2 == USHRT_MAX) return -1;
    return DeltaAngle(tp3s[pt1].Dir, tp3s[pt2].Dir);
  } // KinkAngle

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
      pfp.dEdx[0].resize(slc.nPlanes, 0);
      pfp.dEdx[1].resize(slc.nPlanes, 0);
      pfp.dEdxErr[0].resize(slc.nPlanes, 0);
      pfp.dEdxErr[1].resize(slc.nPlanes, 0);
    }
    for(unsigned short end = 0; end < 2; ++end) {
      // BUG the double brace syntax is required to work around clang bug 21629
      // (https://bugs.llvm.org/show_bug.cgi?id=21629)
      pfp.Dir[end] = {{0.0, 0.0, 0.0}};
      pfp.DirErr[end] = {{0.0, 0.0, 0.0}};
      pfp.XYZ[end] = {{0.0, 0.0, 0.0}};
    }
    return pfp;
  } // CreatePFP
  
  //////////////////////////////////////////
  void FindPFParticles(TCSlice& slc)
  {
    // Match Tjs in 3D and create PFParticles
    
    if(tcc.match3DCuts[0] <= 0) return;
    
    bool prt = (tcc.dbgPFP && tcc.dbgSlc);
    
    if(prt) mf::LogVerbatim("TC")<<" inside FindPFParticles";
    // clear matchVec
    slc.matchVec.clear();
    // clear the kEnvFlag bits on all Tjs. The bit will be set true when a TP is
    // used in a PFParticle Tp3
    for(auto& tj : slc.tjs) {
      for(auto& tp : tj.Pts) tp.Environment[kEnvFlag] = false;
    } // tj
    
    // Match these points in 3D and put the results in slc.matchVec
    std::vector<MatchStruct> matVec;
    // first look for 3-plane matches in a 3-plane TPC
    if(slc.nPlanes == 3) {
      // Match Tjs with high quality vertices first and the leftovers next
      for(short maxScore = 0; maxScore < 2; ++maxScore) FindXMatches(slc, 3, maxScore, matVec, prt);
    } // 3-plane TPC
    // Make 2-plane matches if we haven't hit the user-defined size limit
    if(matVec.size() < tcc.match3DCuts[4]) {
      // 2-plane TPC or 2-plane matches in a 3-plane TPC
      if(slc.nPlanes == 2) {
        for(short maxScore = 0; maxScore < 2; ++maxScore) FindXMatches(slc, 2, maxScore, matVec, prt);
      } else {
        // Make one attempt at 2-plane matches in a 3-plane TPC, setting maxScore large
        FindXMatches(slc, 2, 3, matVec, prt);
      }
    } // can add more combinations
    
    if(matVec.empty()) return;
    
    // sort by decreasing number of matched points
    if(matVec.size() > 1) {
      PFPStruct pfp = CreatePFP(slc);
      std::vector<int> dum1;
      std::vector<float> dum2;
      std::vector<SortEntry> sortVec(matVec.size());
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        sortVec[ii].index = ii;
        auto& ms = matVec[ii];
        sortVec[ii].val = ms.Count;
        // de-weight two-plane matches
        if(ms.TjIDs.size() == 2) sortVec[ii].val /= 2;
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
      std::vector<MatchStruct> tmpVec;
      tmpVec.reserve(matVec.size());
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        tmpVec.push_back(matVec[sortVec[ii].index]);
      } // ii
      matVec = tmpVec;
    } // sort matVec
    
    // put the maybe OK matches into tjs
    PFPStruct pfp = CreatePFP(slc);
    for(auto& ms : matVec) {
      if(ms.Count < 2) continue;
      // check for duplicates
      bool skipit = false;
      for(auto& oms : slc.matchVec) {
        if(ms.TjIDs == oms.TjIDs) {
          skipit = true;
          break;
        }
      } // oms
      if(skipit) continue;
      // Find completeness, do the fit, don't fill Tp3s, don't print
      pfp.TjIDs = ms.TjIDs;
      FindCompleteness(slc, pfp, true, false, false);
      // save the info in matchStruct
      ms.TjCompleteness = pfp.TjCompleteness;
      ms.Pos = pfp.XYZ[0];
      ms.Dir = pfp.Dir[0];
      slc.matchVec.push_back(ms);
    }
    if(slc.matchVec.empty()) return;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FPFP: slc.matchVec\n";
      unsigned short cnt = 0;
      PFPStruct pfp = CreatePFP(slc);
      std::vector<int> dum1;
      std::vector<float> dum2;
      for(unsigned int ii = 0; ii < slc.matchVec.size(); ++ii) {
        auto& ms = slc.matchVec[ii];
        if(ms.Count == 0) continue;
        myprt<<std::setw(4)<<ii<<" Count "<<std::setw(5)<<(int)ms.Count;
        myprt<<" Tj ID-UID:";
        for(auto& tjid : ms.TjIDs) {
          myprt<<" t"<<tjid<<"-T"<<slc.tjs[tjid-1].UID;
        }
        myprt<<" Comp ";
        for(unsigned short itj = 0; itj < ms.TjCompleteness.size(); ++itj) {
          myprt<<std::setprecision(2)<<std::setw(6)<<ms.TjCompleteness[itj];
        }
        myprt<<" Pos ("<<std::setprecision(0)<<std::fixed;
        myprt<<ms.Pos[0]<<", "<<ms.Pos[1]<<", "<<ms.Pos[2];
        myprt<<") Dir "<<std::setprecision(2)<<std::setw(6)<<ms.Dir[0]<<std::setw(6)<<ms.Dir[1]<<std::setw(6)<<ms.Dir[2];
        myprt<<" projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, ms.Pos, ms.Dir, inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, ms.TjIDs);
        myprt<<" MCSMom "<<MCSMom(slc, ms.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(slc, ms.TjIDs, false);
        myprt<<"\n";
        ++cnt;
        if(cnt == 500 || ms.Count < 2) {
          myprt<<"...stopped printing after 500 entries or Count < 2";
          break;
        }
      } // ii
    } // prt

    // create the list of associations to matches that will be converted to PFParticles
    // Start with large count tj matches that have a consistent PDGCode and no vertex attachments
    // and high completeness
    if(slc.matchVec.size() > 1 && slc.matchVec[0].Count > 2 * slc.matchVec[1].Count) {
      auto& ms = slc.matchVec[0];
      int pdgCode = PDGCodeVote(slc, ms.TjIDs, prt);
      bool hasVx = false;
      for(auto tid : ms.TjIDs) {
        auto& tj = slc.tjs[tid - 1];
        if(tj.VtxID[0] > 0 || tj.VtxID[1] > 0) hasVx = true;
      } // tid
      if(pdgCode != 0 && !hasVx) {
        float minCompleteness = 1;
        for(unsigned short itj = 0; itj < ms.TjCompleteness.size(); ++itj) {
          if(ms.TjCompleteness[itj] < minCompleteness) minCompleteness = ms.TjCompleteness[itj];
        } // itj
        if(minCompleteness > 0.5) {
          PFPStruct pfp = CreatePFP(slc);
          pfp.TjIDs = ms.TjIDs;
          // note that the ms position is the average position of all 3D matched Tp3s at this point.
          // It is not the start position. This will be determined in DefinePFP.
          pfp.XYZ[0] = ms.Pos;
          pfp.Dir[0] = ms.Dir;
          pfp.MatchVecIndex = 0;
          // Set the PDGCode so DefinePFP can ignore incompatible matches
          pfp.PDGCode = pdgCode;
          if(DefinePFP("FPFP", slc, pfp, prt) && AnalyzePFP(slc, pfp, prt) && StorePFP(slc, pfp)) {
            ms.Count = 0;
            // clobber MatchStructs that use the Tjs in this pfp
            for(auto& allms : slc.matchVec) {
              auto shared = SetIntersection(allms.TjIDs, pfp.TjIDs);
              if(!shared.empty()) allms.Count = 0;
            } // allms
          } // define/analyze/store PFP
          else {
            if(prt) mf::LogVerbatim("TC")<<" Define/Analyze/Store PFP failed";
          } // define/analyze/store PFP
        } // minCompleteness > 0.5
      } // pdgCode != 0
    } // large count on the first matchVec
    
    // Next consider Tjs attached to 3D vertices. This is only done when reconstructing neutrino events
    if(!tcc.modes[kTestBeam]) {
      Match3DVtxTjs(slc, prt);
    }

    // define the PFParticleList
    for(unsigned int indx = 0; indx < slc.matchVec.size(); ++indx) {
      auto& ms = slc.matchVec[indx];
      // ignore dead matches
      if(ms.Count == 0) continue;
      // skip this match if any of the trajectories is already matched or merged and killed
      bool skipit = false;
      // check for muons and delta rays or InShower Tjs
      bool has13 = false;
      bool has11 = false;
      for(unsigned short itj = 0; itj < ms.TjIDs.size(); ++itj) {
        if(ms.TjIDs[itj] <= 0 || ms.TjIDs[itj] > (int)slc.tjs.size()) {
          std::cout<<"FindPFParticles: bogus ms TjID "<<ms.TjIDs[itj]<<"\n";
          return;
        }
        auto& tj = slc.tjs[ms.TjIDs[itj] - 1];
        if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) skipit = true;
        // skip low TjCompleteness
        if(ms.TjCompleteness.empty()) {
          std::cout<<"FindPFParticles: ms.TjCompleteness is empty\n";
          return;
        }
        if(ms.TjCompleteness[itj] < 0.1) skipit = true;
        if(tj.PDGCode == 13) has13 = true;
        if(tj.PDGCode == 11) has11 = true;
      } // tjID
      if(skipit) continue;
      if(has13 && has11) continue;
      int pdgCode = PDGCodeVote(slc, ms.TjIDs, prt);
      PFPStruct pfp = CreatePFP(slc);
      pfp.TjIDs = ms.TjIDs;
      // note that the ms position is the average position of all 3D matched Tp3s at this point.
      // It is not the start position. This will be determined in DefinePFP.
      pfp.XYZ[0] = ms.Pos;
      pfp.Dir[0] = ms.Dir;
      pfp.MatchVecIndex = indx;
      // Set the PDGCode so DefinePFP can ignore incompatible matches
      pfp.PDGCode = pdgCode;
      // set the PDGCode to ensure that delta rays aren't merged with muons. PDGCodeVote
      // returns 0 if the vote is mixed
      if(has13) pfp.PDGCode = 13;
      if(has11) pfp.PDGCode = 11;
      if(!DefinePFP("FPFP", slc, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" DefinePFP failed";
        continue;
      }
      double sep = 1;
      Split3DKink(slc, pfp, sep, prt);
      if(!AnalyzePFP(slc, pfp, prt)) continue;
      if(!StorePFP(slc, pfp)) {
        if(prt) mf::LogVerbatim("TC")<<" StorePFP failed P"<<pfp.ID;
        continue;
      }
      ms.Count = 0;
      // clobber MatchStructs that use the Tjs in this pfp
      for(auto& allms : slc.matchVec) {
        auto shared = SetIntersection(allms.TjIDs, pfp.TjIDs);
        if(!shared.empty()) allms.Count = 0;
      } // allms
    } // indx
    
//    MatchMissedTjs(slc);

  } // FindPFParticles
  
  /////////////////////////////////////////
  bool DefinePFP(std::string inFcnLabel, TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function is called after the 3D matched TjIDs have been specified and optionally
    // a start or end vertex ID. It defines the PFParticle but doesn't store it
    
    if(pfp.PDGCode == 1111) return false;
    if(pfp.TjIDs.size() < 2) return false;
    
    std::string fcnLabel = inFcnLabel + ".DPFP";
    
    // require at least one tj in at least two planes
    std::vector<unsigned short> nInPln(slc.nPlanes);
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      ++nInPln[DecodeCTP(tj.CTP).Plane];
      if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) {
        std::cout<<fcnLabel<<" P"<<pfp.ID<<" uses T"<<tj.ID<<" but kMat3D is set true\n";
        return false;
      }
    } // tjid
    unsigned short npl = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(nInPln[plane] > 0) ++npl;
    if(npl < 2) return false;
    
    for(unsigned short end = 0; end < 2; ++end) {
      if(pfp.Vx3ID[end] < 0 || pfp.Vx3ID[end] > (int)slc.vtxs.size()) {
        std::cout<<fcnLabel<<" P"<<pfp.ID<<" end "<<end<<" is invalid\n";
        return false;
      } 
    } // end

    if(pfp.Vx3ID[0] == 0 && pfp.Vx3ID[1] > 0) {
      std::cout<<fcnLabel<<" P"<<pfp.ID<<" end 1 has a vertex but end 0 doesn't. No endpoints defined\n";
      return false;
    }
    
    // check for vertex consistency. There should only be one tj in each plane
    // that is attached to a vertex. Remove the shorter tj from the list
    // if that is not the case 
    if(pfp.Vx3ID[0] > 0) PFPVxTjOK(slc, pfp, prt);
    
    bool pfpTrackLike = (MaxTjLen(slc, pfp.TjIDs) > tcc.match3DCuts[5] && MCSMom(slc, pfp.TjIDs) > tcc.match3DCuts[3]);
    // don't look for broken Tjs for primary electron PFPs
    if(pfp.PDGCode == 111) pfpTrackLike = false;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" Vx3ID 3V"<<pfp.Vx3ID[0]<<" 3V"<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
      myprt<<" matchVec index "<<pfp.MatchVecIndex;
      myprt<<" max Tj len "<<MaxTjLen(slc, pfp.TjIDs);
      myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
      myprt<<" pfpTrackLike? "<<pfpTrackLike;
    } // prt
    
    if(tcc.useAlg[kMat3DMerge] && pfpTrackLike && pfp.MatchVecIndex < slc.matchVec.size()) {
      // The index of slc.matchVec has been specified for this pfp so we can look for evidence of
      // broken Tjs starting at the beginning
      for(unsigned short ims = 0; ims < pfp.MatchVecIndex + 10; ++ims) {
        if(ims >= slc.matchVec.size()) break;
        auto& ms = slc.matchVec[ims];
        if(ms.Count == 0) continue;
        std::vector<int> shared = SetIntersection(ms.TjIDs, pfp.TjIDs);
        if(shared.size() < 2) continue;
        // check the max length Tj and cut on MCSMom
        if(MaxTjLen(slc, ms.TjIDs) < tcc.match3DCuts[5]) continue;
        if(MCSMom(slc, ms.TjIDs) < tcc.match3DCuts[3]) continue;
        for(auto tjid : ms.TjIDs) {
          if(std::find(shared.begin(), shared.end(), tjid) != shared.end()) continue;
          auto& tj = slc.tjs[tjid - 1];
          if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
          if(tj.AlgMod[kMat3D]) continue;
          // check for PDGCode compatibility - muons and delta rays
          if(pfp.PDGCode == 13 && tj.PDGCode == 11) continue;
          if(pfp.PDGCode == 11 && tj.PDGCode == 13) continue;
          if(SharesHighScoreVx(slc, pfp, tj)) continue;
          if(tj.MCSMom < tcc.match3DCuts[3]) continue;
          float len = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
          if(len < tcc.match3DCuts[5]) continue;
          // check for a compatible merge
          bool skipit = false;
          for(auto tjid : pfp.TjIDs) {
            auto& ptj = slc.tjs[tjid - 1];
            if(ptj.CTP != tj.CTP) continue;
            if(!CompatibleMerge(slc, tj, ptj, prt)) {
              skipit = true;
              break;
            }
          } // tjid
          if(skipit) continue;
          if(prt) mf::LogVerbatim("TC")<<" add T"<<tjid<<" MCSMom "<<tj.MCSMom<<" length "<<len;
          pfp.TjIDs.push_back(tjid);
          PFPVxTjOK(slc, pfp, prt);
        } // tjid
      } // ims
    } // matchVec index defined
    
/* BB April 25. CheckAndMerge needs work. Shut it off for now.
    // check the completeness of matching points in this set of Tjs and possibly
    // merge Tjs
    if(tcc.useAlg[kMat3DMerge] && pfpTrackLike) {
      if(!CheckAndMerge(slc, pfp, prt)) return false;
    } else {
      // not track-like. Find completeness and fill the TP3s
      FindCompleteness(slc, pfp, true, true, prt);
    }
*/
    // Find completeness and fill the TP3s
    FindCompleteness(slc, pfp, true, true, prt);
    if(pfp.Tp3s.empty()) return false;

    // Set the starting position in 3D if it isn't already defined by a 3D vertex
    SetStart(slc, pfp, prt);
    FollowTp3s(slc, pfp, prt);
    // Check the list of Tjs and merge those that are in the same plane
    MergePFPTjs(slc, pfp, prt);

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" Vx3ID 3V"<<pfp.Vx3ID[0]<<" 3V"<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
      myprt<<" Tp3s size "<<pfp.Tp3s.size();
    }
    
    pfp.NeedsUpdate = false;
    return true;
  } // DefinePFP
  
  /////////////////////////////////////////
  bool PFPVxTjOK(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Checks the PFP Vx3 -> Vx2 -> Tj assignment to see if there is more
    // than one tj in a plane in the pfp.TjIDs list that is attached to the same 2D vertex. 
    // This problem is fixed by removing the shorter tj from the TjIDs list. This function
    // return true if nothing was done to TjIDs
    if(pfp.ID == 0) return true;
    if(pfp.TjIDs.empty()) return true;
    if(pfp.Vx3ID[0] <= 0 || pfp.Vx3ID[0] > (int)slc.vtx3s.size()) return true;
    
    auto& vx3 = slc.vtx3s[pfp.Vx3ID[0] - 1];
    std::vector<int> killMe;
    for(auto vx2id : vx3.Vx2ID) {
      if(vx2id == 0) continue;
      auto& vx2 = slc.vtxs[vx2id - 1];
      auto tjlist = GetVtxTjIDs(slc, vx2);
      auto setInt = SetIntersection(pfp.TjIDs, tjlist);
/*
      std::cout<<"PVTC: P"<<pfp.ID<<" Tjs";
      for(auto tid : pfp.TjIDs) std::cout<<" T"<<tid;
      std::cout<<" set Intersection";
      for(auto tid : setInt) std::cout<<" T"<<tid;
      std::cout<<"\n";
*/
      if(setInt.size() < 2) continue;
      // find the longest one
      int imLong = 0;
      unsigned short lenth = 0;
      for(auto tid : setInt) {
        auto& tj = slc.tjs[tid - 1];
        unsigned short npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        if(npts < lenth) continue;
        lenth = npts;
        imLong = tj.ID;
      } // tid
      if(imLong == 0) continue;
      // add the others to the killMe list
      for(auto tid : setInt) if(tid != imLong) killMe.push_back(tid);
    } // vx2id
    if(killMe.empty()) return true;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"PVTC: P"<<pfp.ID<<" removing short tjs attached to a vertex:";
      for(auto tid : killMe) myprt<<" T"<<tid;
    }
    // re-create the TjIDs vector
    std::vector<int> tmp;
    for(auto tid : pfp.TjIDs) {
      if(std::find(killMe.begin(), killMe.end(), tid) == killMe.end()) tmp.push_back(tid);
    } // tid
    pfp.TjIDs = tmp;
    return false;
  } // PFPVxTjOK
  
  /////////////////////////////////////////
  bool AnalyzePFP(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Analyzes the PFP for oddities and tries to fix them
    if(pfp.ID == 0) return false;
    if(pfp.TjIDs.empty()) return false;
    if(pfp.Tp3s.empty()) return false;
    
    // don't bother analyzing this pfp has been altered
    if(pfp.NeedsUpdate) {
      if(prt) mf::LogVerbatim("TC")<<"AnalyzePFP: P"<<pfp.ID<<" needs to be updated. Skip analysis ";
      return true;
    }
    
    // check the tj completeness
    float minCompleteness = 0.95;
    for(auto tjc : pfp.TjCompleteness) if(tjc < minCompleteness) minCompleteness = tjc;
    if(prt) mf::LogVerbatim("TC")<<"inside AnalyzePFP P"<<pfp.ID<<" minCompleteness "<<minCompleteness;
    if(minCompleteness == 0.95) return true;
    // don't analyze electrons
    if(pfp.PDGCode == 111) return true;
    
    // compare the Tjs in Tp3s with those in TjIDs
    std::vector<int> tjIDs;
    std::vector<unsigned short> tjCnt;
    for(auto& tp3 : pfp.Tp3s) {
      for(auto& tp2 : tp3.Tj2Pts) {
        // convert to int for std::find
        int itjID = tp2.id;
        unsigned short indx = 0;
        for(indx = 0; indx < tjIDs.size(); ++indx) if(tjIDs[indx] == tp2.id) break;
        if(indx == tjIDs.size()) {
          tjIDs.push_back(itjID);
          tjCnt.push_back(1);
        } else {
          ++tjCnt[indx];
        }
      } // tp2
    } // tp3
    
    // look for differences
    for(unsigned short ii = 0; ii < tjIDs.size(); ++ii) {
      if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjIDs[ii]) != pfp.TjIDs.end()) continue;
      auto& missTj = slc.tjs[tjIDs[ii] - 1];
      if(missTj.AlgMod[kMat3D]) continue;
      unsigned short npwc = NumPtsWithCharge(slc, missTj, false);
      if(prt) mf::LogVerbatim("TC")<<" missed T"<<missTj.ID<<" npwc "<<npwc<<" tjCnt "<<tjCnt[ii];
      if(tjCnt[ii] < 0.5 * npwc) continue;
      // June 4, 2018. 3D merging is turned off so require the missed tj to be in
      // a missing plane
      bool skipit = false;
      for(auto tid : pfp.TjIDs) {
        auto& tj = slc.tjs[tid - 1];
        if(tj.CTP == missTj.CTP) skipit = true;
      } // tid
      if(skipit) continue;
      // add the missed Tj to the pfp and flag it as needing an update
      pfp.TjIDs.push_back(missTj.ID);
      if(PFPVxTjOK(slc, pfp, prt)) pfp.NeedsUpdate = true;
    } // ii
    
    if(pfp.NeedsUpdate) DefinePFP("APFP", slc, pfp, prt);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"APFP: Tjs in pfp\n";
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        myprt<<"T"<<tj.ID<<" npwc "<<NumPtsWithCharge(slc, tj, false);
        unsigned short indx = 0;
        for(indx = 0; indx < tjIDs.size(); ++indx) if(tjIDs[indx] == tjid) break;
        if(indx == tjIDs.size()) {
          myprt<<" not found in P"<<pfp.ID<<"\n";
          continue;
        }
        myprt<<" nTp3 "<<tjCnt[indx]<<"\n";
      } // tjid
    } // prt
    return true;
  } // AnalyzePFP
  
  /////////////////////////////////////////
  void PFPVertexCheck(TCSlice& slc)
  {
    // Ensure that all PFParticles have a start vertex. It is possible for
    // PFParticles to be attached to a 3D vertex that is later killed.
    if(!slc.isValid) return;
    
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.Vx3ID[0] > 0) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      vx3.Vx2ID.resize(slc.nPlanes);
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      vx3.ID = slc.vtx3s.size() + 1;
      vx3.Primary = false;
      ++evt.globalS3ID;
      vx3.UID = evt.globalS3ID;
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
          std::cout<<"DPFPR: Missed an end vertex for PFP "<<pfp.ID<<" Write some code\n";
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

/* TODO: This needs work
    if(tcc.modes[kTestBeam]) {
      DefinePFPParentsTestBeam(slc, prt);
      return;
    }
*/
    // associate primary PFParticles with a neutrino PFParticle
    if(neutrinoPFPID > 0) {
      auto& neutrinoPFP = slc.pfps[neutrinoPFPID - 1];
      int vx3id = neutrinoPFP.Vx3ID[1];
      for(auto& pfp : slc.pfps) {
        if(pfp.ID == 0 || pfp.ID == neutrinoPFPID) continue;
        if(pfp.Vx3ID[0] != vx3id) continue;
        pfp.ParentUID = (size_t)neutrinoPFPID;
        pfp.Primary = true;
        neutrinoPFP.DtrUIDs.push_back(pfp.ID);
        if(pfp.PDGCode == 111) neutrinoPFP.PDGCode = 12;
      } // pfp
    } // neutrino PFP exists    
  } // DefinePFPParents
  
  /////////////////////////////////////////
  void DefinePFPParentsTestBeam(TCSlice& slc, bool prt)
  {
    // analog of the one above that was written for neutrino interactions. This differs in that
    // the Tj parent - daughter relationship isn't known yet. If one exists, it is ignored...
    // The assumption here is that all PFParticles that enter (end0) from upstream Z are parents and 
    // any PFParticles attached to them at end1 are daughters. 

    // create a list (stack) of parent ID <-> daughter IDs. The idea is similar to that
    // used in DefineTjParents. A parent-daughter association is made for each entry. After
    // it is made, 1) that entry is removed from the stack, 2) the daughter is checked to see
    // if it a parent of a grand-daughter and if so that pair is added to the stack. 
    std::vector<std::pair<unsigned short, unsigned short>> pardtr;

    // Fill the stack with parents that enter the TPC and have daughters attached to
    // 3D vertices at the other end
    double fidZCut = slc.zLo + 2;
    for(auto& parPFP : slc.pfps) {
      if(parPFP.ID == 0) continue;
      parPFP.Primary = false;
      if(parPFP.XYZ[0][2] > fidZCut) continue;
      parPFP.Primary = true;
      // we found a pfp that entered the TPC. Call it the parent and look for a daughter
      if(prt) mf::LogVerbatim("TC")<<"DPFPTestBeam: parent "<<parPFP.ID<<" end1 vtx "<<parPFP.Vx3ID[1];
      if(parPFP.Vx3ID[1] == 0) continue;
      // There must be other Tjs attached to this vertex which are the daughters. Find them
      // and add them to the pardtr stack
      float score = 0;
      auto& vx3 = slc.vtx3s[parPFP.Vx3ID[1] - 1];
      // ensure that it is valid
      if(vx3.ID == 0) continue;
      // get a list of Tjs attached to this vertex. This will include the Tjs in the parent.
      auto vx3TjList = GetVtxTjIDs(slc, vx3, score);
      if(vx3TjList.empty()) continue;
      // filter out the parent Tjs
      auto dtrTjlist = SetDifference(vx3TjList, parPFP.TjIDs);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" Dtrs:";
        for(auto dtjID : dtrTjlist) myprt<<" "<<dtjID<<"_"<<GetPFPIndex(slc, dtjID);
      }
      // Add to the stack
      for(auto dtjID : dtrTjlist) {
        unsigned short pfpIndex = GetPFPIndex(slc, dtjID);
        if(pfpIndex > slc.pfps.size() - 1) continue;
        unsigned short dtrID = pfpIndex + 1;
        // See if this is a duplicate
        bool duplicate = false;
        for(auto& pd : pardtr) if(parPFP.ID == pd.first && dtrID == pd.second) duplicate = true;
        if(!duplicate) pardtr.push_back(std::make_pair(parPFP.ID, dtrID));
      } // dtjID
    } // parPFP
    
    // iterate through the parent - daughter stack, removing the last pair when a 
    // ParentID is updated and adding pairs for new daughters
    for(unsigned short nit = 0; nit < 100; ++nit) {
      if(pardtr.empty()) break;
      auto lastPair = pardtr[pardtr.size() - 1];
      auto& dtr = slc.pfps[lastPair.second - 1];
      auto& par = slc.pfps[lastPair.first - 1];
      dtr.ParentUID = par.UID;
      par.DtrUIDs.push_back(dtr.UID);
      // remove the last pair
      pardtr.pop_back();
      // Now see if the daughter is a parent. First check for a vertex at the other end.
      // To do that we need to know which end has the vertex between the parent and daughter
      unsigned short dtrEnd = USHRT_MAX;
      for(unsigned short ep = 0; ep < 2; ++ep) {
        if(par.Vx3ID[ep] == 0) continue;
        for(unsigned short ed = 0; ed < 2; ++ed) if(dtr.Vx3ID[ed] == par.Vx3ID[ep]) dtrEnd = ed;
      } // ep
      if(dtrEnd > 1) continue;
      // look at the other end of the daughter
      dtrEnd = 1 - dtrEnd;
      // check for a vertex
      if(dtr.Vx3ID[dtrEnd] == 0) continue;
      // get the list of Tjs attached to it
      auto& vx3 = slc.vtx3s[dtr.Vx3ID[dtrEnd] - 1];
      float score = 0;
      auto vx3TjList = GetVtxTjIDs(slc, vx3, score);
      if(vx3TjList.empty()) continue;
      // filter out the new parent
      auto dtrTjlist = SetDifference(vx3TjList, dtr.TjIDs);
      // put these onto the stack
      for(auto tjid : dtrTjlist) pardtr.push_back(std::make_pair(dtr.ID, tjid));
    } // nit
  } // DefinePFPParentsTestBeam

  ////////////////////////////////////////////////
  bool StorePFP(TCSlice& slc, PFPStruct& pfp)
  {
    // stores the PFParticle in TJStuff
    if(pfp.ID < int(slc.pfps.size())) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    if(!neutrinoPFP) {
      if(pfp.TjIDs.empty()) return false;
      if(pfp.PDGCode != 1111 && pfp.Tp3s.size() < 2) return false;
    }
    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)slc.pfps.size() + 1) pfp.ID = slc.pfps.size() + 1;
    ++evt.globalPFPID;
    pfp.UID = evt.globalPFPID;
    // check the Tjs and set the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.AlgMod[kMat3D]) return false;
      tj.AlgMod[kMat3D] = true;
    } // tjid
    
    if(!pfp.Tp3s.empty()) {
      // ensure that the Tj points are in increasing order and reverse them if they aren't. This
      // presumes that the space points have been ordered from pfp start to pfp end
      std::vector<int> tjids;
      // list of tj points to check for increasing (or decreasing) order
      std::vector<short> firstIpt;
      std::vector<short> lastIpt;
      for(auto& Tp3 : pfp.Tp3s) {
        for(auto& tj2pt : Tp3.Tj2Pts) {
          int tjid = tj2pt.id;
          // check for the first occurrence
          unsigned short ii = 0;
          for(ii = 0; ii < tjids.size(); ++ii) if(tjid == tjids[ii]) break;
          if(ii < tjids.size()) {
            // exists in the list. Keep track of the last point
            lastIpt[ii] = tj2pt.ipt;
            continue;
          }
          tjids.push_back(tjid);
          firstIpt.push_back((short)tj2pt.ipt);
          lastIpt.push_back((short)tj2pt.ipt);
        } // tjpt
      } // spt
      // reverse Tjs if necessary so that end0 is at the start of the pfp
      for(unsigned short ii = 0; ii < tjids.size(); ++ii) {
        // ignore Tjs that aren't associated with this pfp
        if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjids[ii]) == pfp.TjIDs.end()) continue;
        auto& tj = slc.tjs[tjids[ii] - 1];
        if(lastIpt[ii] < firstIpt[ii]) ReverseTraj(slc, tj);
      } // ii
    } // Tp3s exist    
    
    if(pfp.BestPlane < 0) FilldEdx(slc, pfp);
    
    if(pfp.NeedsUpdate) std::cout<<"StorePFP: stored P"<<pfp.ID<<" but NeedsUpdate is true...\n";
    
    slc.pfps.push_back(pfp);
    return true;
  } // StorePFP

  ////////////////////////////////////////////////
  bool InsideFV(TCSlice& slc, PFPStruct& pfp, unsigned short end)
  {
    // returns true if the end of the pfp is inside the fiducial volume of the TPC
    if(pfp.ID <= 0) return false;
    if(end > 1) return false;
    
    float abit = 5;
    auto& pos1 = pfp.XYZ[end];
    return (pos1[0] > slc.xLo + abit && pos1[0] < slc.xHi - abit && 
            pos1[1] > slc.yLo + abit && pos1[1] < slc.yHi - abit &&
            pos1[2] > slc.zLo + abit && pos1[2] < slc.zHi - abit);
    
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
    if(costh > 1) return;
    double sep = PosSep(pos1, pos2);
    if(sep < 1E-6) return;
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
  void ReversePFP(TCSlice& slc, PFPStruct& pfp)
  {
    std::swap(pfp.XYZ[0], pfp.XYZ[1]);
    std::swap(pfp.Dir[0], pfp.Dir[1]);
    for(unsigned short xyz = 0; xyz < 3; ++xyz) {
      pfp.Dir[0][xyz] *= -1;
      pfp.Dir[1][xyz] *= -1;
    }
    std::swap(pfp.DirErr[0], pfp.DirErr[1]);
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    std::reverse(pfp.Tp3s.begin(), pfp.Tp3s.end());
    for(auto& tp3 : pfp.Tp3s) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) {
        tp3.Dir[xyz] *= -1;
        tp3.Dir[xyz] *= -1;
      } // xyz
    } // tp3
  } // ReversePFP
  
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
        if(SignalAtTp(slc, tp)) ++sum;
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

    float sum = 0;
    float cnt = 0;
    // keep track of the lowest value and maybe reject it
    float lo = 1;
    float hi = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      std::vector<int> tjids(1);
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.CTP != inCTP) continue;
        tjids[0] = tjid;
        Point2_t pos;
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        pos[0] = tcc.geom->WireCoordinate(pfp.XYZ[end][1], pfp.XYZ[end][2], planeID);
        if(pos[0] < -0.4) continue;
        // check for dead wires
        unsigned int wire = std::nearbyint(pos[0]);
        if(wire > slc.nWires[plane]) continue;
        if(slc.wireHitRange[plane][wire].first == -1) continue;
        pos[1] = tcc.detprop->ConvertXToTicks(pfp.XYZ[end][0], planeID) * tcc.unitsPerTick;
        float cf = ChgFracNearPos(slc, pos, tjids);
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
  unsigned short FarEnd(TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos)
  {
    // Returns the end (0 or 1) of the pfp that is furthest away from the position pos
    if(pfp.ID == 0) return 0;
    if(PosSep2(pfp.XYZ[1], pos) > PosSep2(pfp.XYZ[0], pos)) return 1;
    return 0;
  } // FarEnd

  ////////////////////////////////////////////////
  void PrintTp3(std::string someText, TCSlice& slc, const TrajPoint3& tp3)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" Pos"<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(6)<<tp3.Pos[0]<<std::setw(6)<<tp3.Pos[1]<<std::setw(6)<<tp3.Pos[2];
    myprt<<" Dir"<<std::fixed<<std::setprecision(3);
    myprt<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
    myprt<<" tj_ipt";
    for(auto tj2pt : tp3.Tj2Pts) {
      auto& tj = slc.tjs[tj2pt.id - 1];
      auto& tp = tj.Pts[tj2pt.ipt];
      myprt<<" "<<tj.ID<<"_"<<PrintPos(slc, tp);
    } // tj2pt
  } // PrintTp3
  
  ////////////////////////////////////////////////
  void PrintTp3s(std::string someText, TCSlice& slc, const PFPStruct& pfp, short printPts)
  {
    if(pfp.Tp3s.empty()) return;
    mf::LogVerbatim myprt("TC");
    if(printPts < 0) {
      // print the head if we are print all points
      myprt<<someText<<" pfp P"<<pfp.ID<<"\n";
      myprt<<someText<<"  ipt ________Pos________ Path   ________Dir_______ along trans  dang  Kink  Tj_ipt \n";
    }
    // print the start
    myprt<<someText<<"    ";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(7)<<pfp.XYZ[0][0]<<std::setw(7)<<pfp.XYZ[0][1]<<std::setw(7)<<pfp.XYZ[0][2];
    myprt<<"     ";
    myprt<<std::fixed<<std::setprecision(2);
    myprt<<std::setw(7)<<pfp.Dir[0][0]<<std::setw(7)<<pfp.Dir[0][1]<<std::setw(7)<<pfp.Dir[0][2];
    myprt<<" <--- pfp.XYZ[0] \n";
    
    unsigned short fromPt = 0;
    unsigned short toPt = pfp.Tp3s.size() - 1;
    if(printPts >= 0) fromPt = toPt;
    Vector3_t prevDir = pfp.Dir[0];
    for(unsigned short ipt = fromPt; ipt <= toPt; ++ipt) {
      auto tp3 = pfp.Tp3s[ipt];
      myprt<<someText<<std::setw(4)<<ipt;
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<std::setw(7)<<tp3.Pos[0]<<std::setw(7)<<tp3.Pos[1]<<std::setw(7)<<tp3.Pos[2];
      myprt<<std::setprecision(1)<<std::setw(5)<<PosSep(tp3.Pos, pfp.XYZ[0]);
      myprt<<std::setprecision(2)<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
      myprt<<std::setprecision(1)<<std::setw(6)<<tp3.AlongTrans[0];
      myprt<<std::setprecision(1)<<std::setw(6)<<tp3.AlongTrans[1];
      myprt<<std::setprecision(2)<<std::setw(7)<<DeltaAngle(prevDir, tp3.Dir);
      prevDir = tp3.Dir;
      // Calculate the kink angle at point ipt, using the two points that are
      // +/- 1 cm on either side of that point
      double sep = 1;
      myprt<<std::setprecision(2)<<std::setw(7)<<KinkAngle(slc, pfp.Tp3s, ipt, sep);
      for(auto tj2pt : tp3.Tj2Pts) {
        auto& tj = slc.tjs[tj2pt.id - 1];
        auto& tp = tj.Pts[tj2pt.ipt];
        myprt<<" "<<tj.ID<<"_"<<tj2pt.ipt<<"_"<<PrintPos(slc, tp);
      } // tj2pt
      myprt<<"\n";
    } // ipt
    // print the end
    myprt<<someText<<"    ";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(7)<<pfp.XYZ[1][0]<<std::setw(7)<<pfp.XYZ[1][1]<<std::setw(7)<<pfp.XYZ[1][2];
    myprt<<"     ";
    myprt<<std::fixed<<std::setprecision(2);
    myprt<<std::setw(7)<<pfp.Dir[1][0]<<std::setw(7)<<pfp.Dir[1][1]<<std::setw(7)<<pfp.Dir[1][2];
    myprt<<" <--- pfp.XYZ[1]. Length "<<PosSep(pfp.XYZ[0], pfp.XYZ[1])<<"\n";
  } // PrintTp3s

} // namespace
