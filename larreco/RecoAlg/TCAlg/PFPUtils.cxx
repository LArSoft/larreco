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
  void UpdateMatchStructs(TjStuff& tjs, int oldTj, int newTj)
  {
    // Replaces tjid and ipt references in tjs.matchVec and tjs.pfps from
    // oldTj to newTj. This function is called when Tjs are split or merged
    // or if a Tj is reversed (in which case oldTj = newTj).
    // The method used is to match the trajectory point positions
    if(oldTj <= 0 || oldTj > (int)tjs.allTraj.size()) return;
    if(newTj <= 0 || newTj > (int)tjs.allTraj.size()) return;
    if(tjs.mallTraj.empty() && tjs.pfps.empty()) return;
    
    // convert from int to unsigned short
    unsigned short oldtjid = oldTj;
    auto& ntj = tjs.allTraj[newTj - 1];
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;

    if(!tjs.mallTraj.empty()) {
      for(unsigned int ipt = 0; ipt < tjs.mallTraj.size(); ++ipt) {
        auto& tj2pt = tjs.mallTraj[ipt];
        if(tj2pt.id > tjs.allTraj.size()) continue;
        if(tj2pt.id != oldtjid) continue;
        // Found the old Tj. Now find the point
        for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
          auto& ntp = ntj.Pts[npt];
          if(std::nearbyint(ntp.Pos[0]) == tj2pt.wire && ntp.Pos[1] > tj2pt.xlo && ntp.Pos[1] < tj2pt.xhi) {
            tj2pt.id = newTj;
            tj2pt.ipt = npt;
            tj2pt.npts = npts;
            break;
          } // points match
        } // npt
      } // ipt
    } // !tjs.mallTraj.empty()
    
    // Update pfp space points
    if(!tjs.pfps.empty()) {
      for(auto& pfp : tjs.pfps) {
        for(auto& tp3 : pfp.Tp3s) {
          // check each of the Tj2Pts associated with this space point
          for(auto& tj2pt : tp3.Tj2Pts) {
            if(tj2pt.id > tjs.allTraj.size()) continue;
            if(tj2pt.id != oldtjid) continue;
            // look for the corresponding point (wire) on the new Tj
            for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
              auto& ntp = ntj.Pts[npt];
              if(std::nearbyint(ntj.Pts[npt].Pos[0]) == tj2pt.wire && ntp.Pos[1] > tj2pt.xlo && ntp.Pos[1] < tj2pt.xhi) {
                tj2pt.id = newTj;
                tj2pt.ipt = npt;
                tj2pt.npts = npts;
                break;
              }
            } // npt
          } // tj2pt
        } // spt
      } // pfp
    } // pfps exists

  } // UpdateMatchStructs
  
  /////////////////////////////////////////
  void FillmAllTraj(TjStuff& tjs, const geo::TPCID& tpcid) 
  {
    // Fills the tjs.mallTraj vector with trajectory points in the tpc and sorts
    // them by increasing X
    tjs.matchVec.clear();
    
    int cstat = tpcid.Cryostat;
    int tpc = tpcid.TPC;
    
    // count the number of TPs and clear out any old 3D match flags
    unsigned int ntp = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      // don't match InShower Tjs
//      if(tj.AlgMod[kInShower]) continue;
      // or Shower Tjs
//      if(tj.AlgMod[kShowerTj]) continue;
      if(tj.ID <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      ntp += NumPtsWithCharge(tjs, tj, false);
      tj.AlgMod[kMat3D] = false;
    } // tj
    if(ntp < 2) return;
    
    tjs.mallTraj.resize(ntp);
    
    // define mallTraj
    unsigned int icnt = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      // don't match shower-like Tjs
//      if(tj.AlgMod[kInShower]) continue;
      // or Shower Tjs
//      if(tj.AlgMod[kShowerTj]) continue;
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
        if(icnt > tjs.mallTraj.size() - 1) break;
        tjs.mallTraj[icnt].wire = std::nearbyint(tp.Pos[0]);
        bool hasWire = tjs.geom->HasWire(geo::WireID(cstat, tpc, plane, tjs.mallTraj[icnt].wire));
        // don't try matching if the wire doesn't exist
        if(!hasWire) continue;
        float xpos = tjs.detprop->ConvertTicksToX(tp.Pos[1]/tjs.UnitsPerTick, plane, tpc, cstat);
        float posPlusRMS = tp.Pos[1] + TPHitsRMSTime(tjs, tp, kUsedHits);
        float rms = tjs.detprop->ConvertTicksToX(posPlusRMS/tjs.UnitsPerTick, plane, tpc, cstat) - xpos;
        if(rms < tjs.Match3DCuts[0]) rms = tjs.Match3DCuts[0];
        if(icnt == tjs.mallTraj.size()) break;
        tjs.mallTraj[icnt].xlo = xpos - rms;
        tjs.mallTraj[icnt].xhi = xpos + rms;
        tjs.mallTraj[icnt].dir = tp.Dir;
        tjs.mallTraj[icnt].ctp = tp.CTP;
        tjs.mallTraj[icnt].id = tjID;
        tjs.mallTraj[icnt].ipt = ipt;
        tjs.mallTraj[icnt].npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        tjs.mallTraj[icnt].score = score;
        ++icnt;
      } // tp
    } // tj
    
    if(icnt < tjs.mallTraj.size()) tjs.mallTraj.resize(icnt);
    
    // This is pretty self-explanatory
    std::vector<SortEntry> sortVec(tjs.mallTraj.size());
    for(unsigned int ipt = 0; ipt < tjs.mallTraj.size(); ++ipt) {
      // populate the sort vector
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = tjs.mallTraj[ipt].xlo;
    } // ipt
    // sort by increasing xlo
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put tjs.mallTraj into sorted order
    auto tallTraj = tjs.mallTraj;
    for(unsigned int ii = 0; ii < sortVec.size(); ++ii) tjs.mallTraj[ii] = tallTraj[sortVec[ii].index];
    
  } // FillmAllTraj
  
  /////////////////////////////////////////
  void SetEndPoints(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Use the end points of the pfp to see if the track is very small angle
    // wrt to the wire planes. If that is the case, set the DirectionFixed flag true
    // and use the end points to find the start and end xyz positions
    
    if(pfp.TjIDs.empty()) return;
    // do an initial simple check using the Tj angles
    pfp.DirectionFixed = true;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      unsigned short endPt = tj.EndPt[0];
      auto& tp0 = tj.Pts[endPt];
      endPt = tj.EndPt[1];
      auto& tp1 = tj.Pts[endPt];
      float slp = (tp1.Pos[1] - tp0.Pos[1]) / (tp1.Pos[0] - tp0.Pos[0]);
      if(std::abs(slp) > 0.1) pfp.DirectionFixed = false;
    } // tj2pt
    
//    if(!pfp.DirectionFixed) return;
    
    if(prt) mf::LogVerbatim("TC")<<"SEP: pfp "<<pfp.ID<<" is small angle";
    
    std::vector<unsigned short> tjEnd0(pfp.TjIDs.size());
    
    double maxSep = 0;
    for(unsigned short ii = 0; ii < pfp.TjIDs.size() - 1; ++ii) {
      auto& itj = tjs.allTraj[pfp.TjIDs[ii] - 1];
      geo::PlaneID iPlaneID = DecodeCTP(itj.CTP);
      unsigned int cstat = iPlaneID.Cryostat;
      unsigned int tpc = iPlaneID.TPC;
      for(unsigned short jj = ii + 1; jj < pfp.TjIDs.size(); ++jj) {
        auto& jtj = tjs.allTraj[pfp.TjIDs[jj] - 1];
        geo::PlaneID jPlaneID = DecodeCTP(jtj.CTP);
        for(unsigned short iend = 0; iend < 2; ++iend) {
          unsigned short jend = iend;
          if(itj.StepDir != jtj.StepDir) jend = 1 - iend;
          unsigned int iwire = std::nearbyint(itj.Pts[itj.EndPt[iend]].Pos[0]);
          unsigned int jwire = std::nearbyint(jtj.Pts[jtj.EndPt[jend]].Pos[0]);
          Point3_t pt1, pt2;
          tjs.geom->IntersectionPoint(iwire, jwire, iPlaneID.Plane, jPlaneID.Plane, cstat, tpc, pt1[1], pt1[2]);
          // now find the intersection point using the other ends
          iwire = std::nearbyint(itj.Pts[itj.EndPt[1 - iend]].Pos[0]);
          jwire = std::nearbyint(jtj.Pts[jtj.EndPt[1 - jend]].Pos[0]);
          tjs.geom->IntersectionPoint(iwire, jwire, iPlaneID.Plane, jPlaneID.Plane, cstat, tpc, pt2[1], pt2[2]);
          double sep2 = PosSep2(pt1, pt2);
          if(sep2 > maxSep) {
            maxSep = sep2;
            double ticks = itj.Pts[itj.EndPt[iend]].Pos[1] / tjs.UnitsPerTick;
            pt1[0] = tjs.detprop->ConvertTicksToX(ticks, iPlaneID);
            ticks = jtj.Pts[jtj.EndPt[jend]].Pos[1] / tjs.UnitsPerTick;
            pt1[0] += tjs.detprop->ConvertTicksToX(ticks, jPlaneID);
            pt1[0] /= 2;
            // and the other ends
            ticks = itj.Pts[itj.EndPt[1 - iend]].Pos[1] / tjs.UnitsPerTick;
            pt2[0] = tjs.detprop->ConvertTicksToX(ticks, iPlaneID);
            ticks = jtj.Pts[jtj.EndPt[1 - jend]].Pos[1] / tjs.UnitsPerTick;
            pt2[0] += tjs.detprop->ConvertTicksToX(ticks, jPlaneID);
            pt2[0] /= 2;
            // Declare the start to be at the low Z end
            if(pt1[2] > pt2[2]) std::swap(pt1, pt2);
            pfp.XYZ[0] = pt1;
            pfp.XYZ[1] = pt2;
            Vector3_t dir;
            for(unsigned short xyz = 0; xyz < 3; ++xyz) dir[xyz] = pt2[xyz] - pt1[xyz];
            SetMag(dir, 1);
            pfp.Dir[0] = dir;
            pfp.Dir[1] = dir;
            // keep track of which ends these are
            tjEnd0[ii] = iend;
            tjEnd0[jj] = jend;
          } // sep2 > maxSep
        } // iend
      } // jj
    } // ii
    if(maxSep == 0) return;
    // try to attach vertices to the ends
    AttachVertices(tjs, pfp, prt);
    if(prt) {
      for(unsigned short ii = 0; ii < pfp.TjIDs.size(); ++ii) {
        unsigned short iend = tjEnd0[ii];
        auto& tj = tjs.allTraj[pfp.TjIDs[ii] - 1];
        unsigned short endPt0 = tj.EndPt[iend];
        unsigned short endPt1 = tj.EndPt[1 - iend];
        mf::LogVerbatim("TC")<<"Tj "<<tj.ID<<" endPt0 "<<PrintPos(tjs, tj.Pts[endPt0].Pos)<<" endPt1 "<<PrintPos(tjs, tj.Pts[endPt1].Pos);
      } // ii
      PrintPFP("SEP", tjs, pfp, true);
    } // prt
  } // SetEndPoints
  
  /////////////////////////////////////////
  void AttachVertices(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // try to attach a vertex on both ends of the pfp using the positions XYZ[end]. This function
    // doesn't need to have Tp3s associated with the pfp but probably doesn't work as well
    if(pfp.ID == 0) return;
    std::array<unsigned short, 2> imbest {0};
    // Ignore any separation larger than (10 cm)^2
    std::array<float, 2> best {100};
    for(unsigned short end = 0; end < 2; ++end) pfp.Vx3ID[end] = 0;
    auto vx3list = GetPFPVertices(tjs, pfp);
    for(unsigned short end = 0; end < 2; ++end) {
      for(auto vx3id : vx3list) {
        auto& vx3 = tjs.vtx3[vx3id - 1];
        Point3_t vxpos = {vx3.X, vx3.Y, vx3.Z};
        float sep2 = PosSep2(pfp.XYZ[end], vxpos);
        if(sep2 > best[end]) continue;
        best[end] = sep2;
        imbest[end] = vx3id;
      } // vx3id
    } // end
    if(imbest[0] == 0 && imbest[1] == 0) return;
    if(imbest[0] == imbest[1]) {
      // The same vertex meets the cuts for both ends. Take the best one.
      if(best[0] < best[1]) {
        pfp.Vx3ID[0] = imbest[0];
        if(prt) mf::LogVerbatim("TC")<<"AV0: attach vx "<<imbest[0]<<" to pfp "<<pfp.ID<<" to end 0";
      } else {
        pfp.Vx3ID[1] = imbest[0];
        if(prt) mf::LogVerbatim("TC")<<"AV1: attach vx "<<imbest[0]<<" to pfp "<<pfp.ID<<" to end 1";
      }
      return;
    } else {
      // have different vertices on each end
      for(unsigned short end = 0; end < 2; ++end) {
        if(imbest[end] == 0) continue;
        pfp.Vx3ID[end] = imbest[end];
        if(prt) mf::LogVerbatim("TC")<<"AV: attach vx "<<imbest[end]<<" to pfp "<<pfp.ID<<" end "<<end;
      } // end
    } // have different vertices (or no vertices) on each end
  } // AttachVertices

/*
  /////////////////////////////////////////
  void MakePFPTp3s(TjStuff& tjs, PFPStruct& pfp, bool anyTj)
  {
    // Creates a vector of 3D trajectory points
    // TODO This code needs serious re-thinking...
    
    pfp.Tp3s.clear();
    if(pfp.ID == 0) return;
    if(pfp.TjIDs.size() < 2) return;
    if(tjs.mallTraj.empty()) return;
    
    // Make a Tp3 for every point on every trajectory
    unsigned short tp3Size = 0;
    // We will populate pfp.Tp3s starting with the longest Tj so sort by
    // decreasing length
    std::vector<SortEntry> sortVec(pfp.TjIDs.size());
    unsigned short cnt = 0;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      unsigned short tjlen = tj.EndPt[1] - tj.EndPt[0] + 1;
      sortVec[cnt].index = cnt;
      sortVec[cnt].val = tjlen;
      ++cnt;
      tp3Size += tjlen;
    } // tjid
    // make a temp vector that is somewhat oversized to reduce the probability of
    // getting multiple Tj2Pts in a Tp3
    unsigned short vecSize = 1.2 * tp3Size;
    std::vector<TrajPoint3> tp3Vec(vecSize);
    std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
    // This code assumes that the Tjs are in the correct order away from the
    // start point which is at end0 of each Tj. The only thing we need to populate
    // in the Tp3 Tj2Pt is tj ID and ipt
    Tj2Pt tj2pt;
    for(unsigned short ii = 0; ii < pfp.TjIDs.size(); ++ii) {
      int tjid = pfp.TjIDs[sortVec[ii].index];
      auto& tj = tjs.allTraj[tjid - 1];
      // calculate an index scale factor
      float tjlen = sortVec[ii].val;
      // This is structured so that the first (longest) tj has it's first entry in tp3Vec[0],
      // the second longest has it's first entry in tp3Vec[1], etc. The last entry in tp3Vec
      // should be the Tj2Pt of the longest Tj, etc
      float scaleF = (float)(vecSize - 2 * ii) / (tjlen - 1);
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        tj2pt.ctp = tj.CTP;
        tj2pt.id = tjid;
        tj2pt.ipt = ipt;
        geo::PlaneID planeID = DecodeCTP(tp.CTP);
        tj2pt.xlo = tjs.detprop->ConvertTicksToX(tp.Pos[1]/tjs.UnitsPerTick, planeID);
        // calculate an index into tp3Vec so that the points are distributed uniformly (not
        // necessarily uniformly in path length)
        unsigned short indx = std::nearbyint(ii + scaleF * (ipt - tj.EndPt[0]));
        if(indx > vecSize) indx = vecSize;
        // After the first tj points have been entered, refine the index for the second and third
        // tjs so that the X positions are consistent
        if(ii > 0) std::cout<<pfp.ID<<" indx in "<<indx;
        if(ii > 0) NudgeIndex(tjs, tp3Vec, tj2pt, indx);
        if(ii > 0) std::cout<<pfp.ID<<" out "<<indx<<"\n";
        if(indx == USHRT_MAX) continue;
        if(indx > tp3Vec.size() - 1) tp3Vec.resize(indx+1);
        tp3Vec[indx].Tj2Pts.push_back(tj2pt);
        tp3Vec[indx].Pos[0] = tj2pt.xlo;
        ++cnt;
      } // ipt
    } // tj
    // Transfer the valid points into pfp.Tp3s
    for(auto& tp3 : tp3Vec) {
      if(!tp3.Tj2Pts.empty()) pfp.Tp3s.push_back(tp3);
    } //
    // Find the positions and directions at each point using adjacent points
    TrajPoint3 tmp;
    // make local copies of the points so they can be over-written
    TrajPoint itp2, jtp2;
    for(unsigned short indx = 0; indx < pfp.Tp3s.size(); ++indx) {
      auto& itp3 = pfp.Tp3s[indx];
      if(itp3.Tj2Pts.empty()) continue;
      // already defined?
      if(itp3.Pos[0] != 0 && itp3.Pos[1] != 0) continue;
      // get a reference to the Tj point
      auto& itj2pt = itp3.Tj2Pts[0];
      itp2 = tjs.allTraj[itj2pt.id - 1].Pts[itj2pt.ipt];
      // have a singlet. Check the next Tp3 for a Tj point
      unsigned short jndx = indx + 1;
      if(jndx == pfp.Tp3s.size() - 1) break;
      // get a reference to the next Tp3
      auto& jtp3 = pfp.Tp3s[jndx];
      if(jtp3.Tj2Pts.empty()) continue;
      // get a reference to the Tj point
      auto& jtj2pt = jtp3.Tj2Pts[0];
      // try an earlier point if the next one is in the same plane
      if(jtj2pt.ctp == itj2pt.ctp) {
        if(indx == 0) continue;
        auto& jtp3 = pfp.Tp3s[indx - 1];
        auto& jtj2pt = jtp3.Tj2Pts[0];
        jtp2 = tjs.allTraj[jtj2pt.id - 1].Pts[jtj2pt.ipt];
      } else {
        jtp2 = tjs.allTraj[jtj2pt.id - 1].Pts[jtj2pt.ipt];
      }
      // Don't attempt if they are too far apart
      if(std::abs(jtp3.Pos[0] - itp3.Pos[0]) > tjs.Match3DCuts[0]) continue;
      if(MakeTp3(tjs, itp2, jtp2, tmp)) {
        itp3.Pos = tmp.Pos;
        itp3.Dir = tmp.Dir;
      }
    } // indx
    FixDirection(tjs, pfp);
    PrintTp3s("MTp", tjs, pfp, -1);
  } // MakePFPTp3s
  
  /////////////////////////////////////////
  void NudgeIndex(TjStuff& tjs, const std::vector<TrajPoint3>& tp3Vec, const Tj2Pt& tj2pt, unsigned short& indx)
  {
    // We want to put the tj2pt in tp3Vec at an appropriate point. The index into the vector was calculated
    // using the assumption that all Tj points on all trajectories should be uniformly distributed in this
    // vector but that is only an approximation. This function revises the indx such that 1) there is not
    // already a point at that position and 2) that the X positions of nearby points are consistent
    
    // find the x range
    float x0 = -1E6, x1 = 0;
    for(unsigned short ipt = 0; ipt < tp3Vec.size(); ++ipt) {
      auto& tp3 = tp3Vec[ipt];
      if(tp3.Tj2Pts.empty()) continue;
      if(x0 == -1) x0 = tp3.Pos[0];
      x1 = tp3.Pos[0];
    } // tp3
    if(x0 == -1E6) return;
    bool smallXRange = (std::abs(x1 - x0) < 2 * tjs.Match3DCuts[0]);
    
    if(smallXRange) {
      // Nothing to be done if small X range and the slot is available
      if(tp3Vec[indx].Tj2Pts.empty()) return;
      // find an available slot
      for(unsigned short ii = 0; ii < 5; ++ii) {
        unsigned short next = indx + ii;
        if(next < tp3Vec.size() && tp3Vec[next].Tj2Pts.empty()) {
          indx = next;
          return;
        }
      } // ii
    } // small X range
    
    unsigned short closeXPt = USHRT_MAX;
    for(unsigned short ii = 0; ii < tp3Vec.size(); ++ii) {
      unsigned short ipt = indx + ii;
      if(ipt < tp3Vec.size() && !tp3Vec[ipt].Tj2Pts.empty()) {
        float dx = std::abs(tj2pt.xlo - tp3Vec[ipt].Pos[0]);
        if(dx < tjs.Match3DCuts[0]) {
          closeXPt = ipt;
          break;
        } // good dx
      } // valid ipt
      ipt = indx - ii;
      if(ipt < tp3Vec.size() && !tp3Vec[ipt].Tj2Pts.empty()) {
        float dx = std::abs(tj2pt.xlo - tp3Vec[ipt].Pos[0]);
        if(dx < tjs.Match3DCuts[0]) {
          closeXPt = ipt;
          break;
        } // good dx
      } // valid ipt
    } // ii
    
    if(closeXPt == USHRT_MAX) {
      indx = USHRT_MAX;
      return;
    }
    // now find an available slot near this point
    for(unsigned short ii = 1; ii < tp3Vec.size(); ++ii) {
      unsigned short ipt = closeXPt + ii;
      if(ipt == tp3Vec.size()) break;
      // look for a slot with a valid Tp3 and ensure that we haven't moved too far away
      if(!tp3Vec[ipt].Tj2Pts.empty() && std::abs(tj2pt.xlo - tp3Vec[ipt].Pos[0]) > tjs.Match3DCuts[0]) break;
      if(ipt < tp3Vec.size() && tp3Vec[ipt].Tj2Pts.empty()) {
        indx = ipt;
        return;
      }
    } // ii
    // Try going in the other direction
    for(unsigned short ii = 1; ii < tp3Vec.size(); ++ii) {
      unsigned short ipt = closeXPt - ii;
      if(ipt > tp3Vec.size() - 1) break;
      // look for a slot with a valid Tp3 and ensure that we haven't moved too far away
      if(!tp3Vec[ipt].Tj2Pts.empty() && std::abs(tj2pt.xlo - tp3Vec[ipt].Pos[0]) > tjs.Match3DCuts[0]) break;
      if(ipt < tp3Vec.size() && tp3Vec[ipt].Tj2Pts.empty()) {
        indx = ipt;
        return;
      }
    } // ii
        
  } // NudgeIndex
*/
  /////////////////////////////////////////
  bool SetNewStart(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Analyzes the space point collection and the Tjs in the pfp to create a start
    // position. This function returns false if a failure occurs
    if(pfp.ID == 0 || pfp.TjIDs.empty()) return false;
    if(pfp.Tp3s.size() < 2) return false;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"SNS: pfp "<<pfp.ID<<" Vx3ID[0] "<<pfp.Vx3ID[0]<<" DirectionFixed? "<<pfp.DirectionFixed;
      PrintTp3("First", tjs, pfp.Tp3s[0]);
      PrintTp3("Last", tjs, pfp.Tp3s[pfp.Tp3s.size() - 1]);
    }

    unsigned short pt1 = 0;
    unsigned short pt2 = pfp.Tp3s.size() - 1;
    
    // swap them to put Z_pt2 > Z_pt1
    if(pfp.Tp3s[pt1].Pos[2] > pfp.Tp3s[pt2].Pos[2]) std::swap(pt1, pt2);

    float fidCut = 2;
    bool entersUS = pfp.Tp3s[pt1].Pos[2] < tjs.ZLo + fidCut;
    bool leavesDS = pfp.Tp3s[pt2].Pos[2] > tjs.ZHi - fidCut;
    bool insideZ = !entersUS && !leavesDS;
    
    bool entersLoX = pfp.Tp3s[pt1].Pos[0] < tjs.XLo + fidCut;
    bool leavesHiX = pfp.Tp3s[pt2].Pos[0] > tjs.XHi - fidCut;
    bool insideX = !entersLoX && !leavesHiX;
    
    bool entersLoY = pfp.Tp3s[pt1].Pos[1] < tjs.YLo + fidCut;
    bool leavesHiY = pfp.Tp3s[pt2].Pos[1] > tjs.YHi - fidCut;
    bool insideY = !entersLoY && !leavesHiY;
    
    unsigned short startPt = 0, endPt = 0;
    if(entersUS && leavesDS) {
      // Through-going - probably a beam particle
      startPt = pt1;
      endPt = pt2;
      pfp.StopFlag[0][kOutFV] = true;
      pfp.StopFlag[1][kOutFV] = true;
    } else if(entersUS && insideX && insideY) {
      // enters US and doesn't leave
      startPt = pt1;
      endPt = pt2;
      pfp.StopFlag[0][kOutFV] = true;
      pfp.StopFlag[1][kOutFV] = false;
    } else if(insideX && insideY && insideZ) {
      // fully contained - assume it is beam related for now
      startPt = pt1;
      endPt = pt2;
      pfp.StopFlag[0][kOutFV] = false;
      pfp.StopFlag[1][kOutFV] = false;
    } else {
      // assume it is a cosmic ray entering from above
      if(pfp.Tp3s[pt1].Pos[1] > pfp.Tp3s[pt2].Pos[1]) {
        startPt = pt1;
        endPt = pt2;
        pfp.StopFlag[0][kOutFV] = (entersUS || entersLoX || entersLoX);
        pfp.StopFlag[1][kOutFV] = (leavesDS || leavesHiX || leavesHiX);
      } else {
        startPt = pt2;
        endPt = pt1;
        pfp.StopFlag[0][kOutFV] = (leavesDS || leavesHiX || leavesHiX);
        pfp.StopFlag[1][kOutFV] = (entersUS || entersLoX || entersLoX);
      }
    } // Not through-going or fully contained
    
    // TODO: Make additional tests for stopping tracks and muons with delta-rays
    
    // Define the start and end positions
    pfp.XYZ[0] = pfp.Tp3s[startPt].Pos;
    pfp.XYZ[1] = pfp.Tp3s[endPt].Pos;
    // and the general direction
    for(unsigned short xyz = 0; xyz < 3; ++xyz) pfp.Dir[0][xyz] = pfp.XYZ[1][xyz] - pfp.XYZ[0][xyz];
    SetMag(pfp.Dir[0], 1);
    pfp.Dir[1] = pfp.Dir[0];

    if(prt) {
      PrintTp3("SNSs", tjs, pfp.Tp3s[startPt]);
      PrintTp3("SNSe", tjs, pfp.Tp3s[endPt]);
      mf::LogVerbatim myprt("TC");
      myprt<<" SNS: XYZ[0] "<<std::fixed<<std::setprecision(1);
      myprt<<std::setw(7)<<pfp.XYZ[0][0];
      myprt<<std::setw(7)<<pfp.XYZ[0][1];
      myprt<<std::setw(7)<<pfp.XYZ[0][2];
      myprt<<" XYZ[0] outFV? "<<pfp.StopFlag[0][kOutFV];
      myprt<<" XYZ[1] ";
      myprt<<std::setw(7)<<pfp.XYZ[1][0];
      myprt<<std::setw(7)<<pfp.XYZ[1][1];
      myprt<<std::setw(7)<<pfp.XYZ[1][2];
      myprt<<" XYZ[1] outFV? "<<pfp.StopFlag[1][kOutFV];
      myprt<<" Dir[0] "<<std::fixed<<std::setprecision(3);
      myprt<<std::setw(7)<<pfp.Dir[0][0];
      myprt<<std::setw(7)<<pfp.Dir[0][1];
      myprt<<std::setw(7)<<pfp.Dir[0][2];
    }

    return true;
    
  } // SetNewStart

  /////////////////////////////////////////
  void SetEndVx(TjStuff& tjs, PFPStruct& pfp, unsigned short atEnd, bool prt)
  {
    // Analyzes the requested end of Tp3s to see if Tjs are attached to the same 3D vertex
    // and if so, attach the PFParticle that vertex
    if(atEnd > 1) return;
    if(pfp.Tp3s.empty()) return;
    //  already attached to a valid?
    if(pfp.Vx3ID[atEnd] > 0 && tjs.vtx3[pfp.Vx3ID[atEnd] - 1].Wire != -2) {
      if(prt) mf::LogVerbatim("TC")<<"SEV: pfp "<<pfp.ID<<" Vx3ID["<<atEnd<<"] = "<<pfp.Vx3ID[atEnd];
      return;
    }
    
    // clobber the PFP-only vertex that was used to define the start position
    if(pfp.Vx3ID[atEnd] > 0 && tjs.vtx3[pfp.Vx3ID[atEnd] - 1].Wire == -2) {
      auto& vx3 = tjs.vtx3[pfp.Vx3ID[atEnd] - 1];
      for(auto& opfp : tjs.pfps) {
        for(unsigned short end = 0; end < 2; ++end) {
          if(opfp.Vx3ID[end] == vx3.ID) opfp.Vx3ID[end] = 0;
        } // end
      } // opfp
      vx3.ID = 0;
    } // PFP-only vertex exists
    
    // make a list of 3D vertices at the end of Tp3s
    std::vector<unsigned short> endVxList;
    // and the Tjs to which they are attached
    std::vector<int> endTjList;
    // check a number of points near the end
    for(unsigned short ii = 0; ii < 5; ++ii) {
      short ipt = 0;
      if(atEnd == 1) {
        ipt = pfp.Tp3s.size() - ii - 1;
        // don't get too close to the start if this is a short pfp
        if(ipt < 3) break;
      } else {
        ipt = ii;
        if(ipt > (short)pfp.Tp3s.size() - 3) break;
      }
      auto& tp3 = pfp.Tp3s[ipt];
      for(auto& tj2pt : tp3.Tj2Pts) {
        auto& tj = tjs.allTraj[tj2pt.id - 1];
        unsigned short end = 0;
        unsigned short midPt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
        if(ipt > midPt) end = 1;
        // see if there is a 2D vertex at this end
        if(tj.VtxID[end] == 0 || tj.VtxID[end] > tjs.vtx.size()) continue;
        auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
        // see if this is matched to a 3D vertex
        if(vx2.Vx3ID == 0 || vx2.Vx3ID > tjs.vtx3.size()) continue;
        // Ignore it if it is already in the Tj list
        if(std::find(endTjList.begin(), endTjList.end(), tj.ID) != endTjList.end()) continue;
        // Ignore it if it is attached to the other end
        if(vx2.Vx3ID == pfp.Vx3ID[1 - atEnd]) continue;
        // Ignore it if it is closer to the other end
        unsigned short opt = 0;
        if(atEnd == 1) opt = pfp.Tp3s.size() - 1;
        auto& otp3 = pfp.Tp3s[opt];
        // put the position into a Point3_t
        auto& vx3 = tjs.vtx3[vx2.Vx3ID - 1];
        Point3_t vx3pos = {vx3.X, vx3.Y, vx3.Z};
        if(PosSep2(tp3.Pos, vx3pos) > PosSep2(otp3.Pos, vx3pos)) continue;
        // add it to the list
//        std::cout<<"pfp_atEnd "<<pfp.ID<<"_"<<atEnd<<" tj_end "<<tj.ID<<"_"<<end<<" vx3 "<<vx2.Vx3ID<<"\n";
        endVxList.push_back(vx2.Vx3ID);
        endTjList.push_back(tj.ID);
      } // tj2pt
      // break out if we have them all. 
      // TODO: Compare the endTjList with pfp.TjIDs here? And do what?...
      if(endTjList.size() == pfp.TjIDs.size()) break;
    } // ii
    if(endVxList.empty()) return;
    // Just take the first one
    // TODO: This may need to be done more carefully if there is vertex confusion at this end
    pfp.Vx3ID[atEnd] = endVxList[0];
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"SEV: pfp "<<pfp.ID<<" checking atEnd "<<atEnd<<" End Vx_Tj:";
      for(unsigned short iev = 0; iev < endVxList.size(); ++iev) {
        myprt<<" "<<endVxList[iev]<<"_"<<endTjList[iev];
      }
      myprt<<" Setting Vx3ID "<<pfp.Vx3ID[atEnd];
    } // prt
  } // SetEndVx

  /////////////////////////////////////////
  void SortByDistanceFromStart(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // sorts the TC spacepoints by distance from the start vertex pfp.XYZ[0]. This will result
    // in sorting by the distance from (0, 0, 0) if no start position has been specified
    
    if(pfp.Tp3s.size() < 2) return;
    
    if(pfp.Dir[0][0] == 0 && pfp.Dir[0][2] == 0) {
      std::cout<<"SBDFS: direction not set... \n";
      return;
    }

    std::vector<SortEntry> sortVec(pfp.Tp3s.size());
    
    for(unsigned short ii = 0; ii < pfp.Tp3s.size(); ++ii) {
      sortVec[ii].index = ii;
      // Use the distance along the direction vector from start to end
      Vector3_t sep;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) sep[xyz] = pfp.Tp3s[ii].Pos[xyz] - pfp.XYZ[0][xyz];
      sortVec[ii].val = DotProd(pfp.Dir[0], sep);
    } // ii

    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put them into order
    std::vector<TrajPoint3> temp;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) temp.push_back(pfp.Tp3s[sortVec[ii].index]);
    pfp.Tp3s = temp;

    // Don't attempt to set the direction since it has already been done
    if(pfp.DirectionFixed) return;
    
    // We know the start position. Take an average of the direction using points
    // at the start to get an estimate of the direction
    // Put this in pfp.Dir[0]
    auto& startPos = pfp.Tp3s[0].Pos;
    auto& endPos = pfp.Tp3s[pfp.Tp3s.size() - 1].Pos;
    pfp.Dir[0] = PointDirection(startPos, endPos);

    // Set the direction vectors of all points to be consistent with the general
    // start direction
    FixDirection(tjs, pfp);
    Vector3_t startDir = {0};
    // Find the average direction using the points in the first 10 cm
    for(auto& tp3 : pfp.Tp3s) {
      if(PosSep2(tp3.Pos, startPos) > 100) break;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) startDir[xyz] += tp3.Dir[xyz];
    } //  tp3
    SetMag(startDir, 1);
    pfp.Dir[0] = startDir;
    if(prt) mf::LogVerbatim("TC")<<"SBDFS: start direction "<<std::fixed<<std::setprecision(2)<<pfp.Dir[0][0]<<" "<<pfp.Dir[0][1]<<" "<<pfp.Dir[0][2];
    // do the same at the other end
    Vector3_t endDir = {0};
    for(unsigned short ii = 0; ii < pfp.Tp3s.size(); ++ii) {
      auto& tp3 = pfp.Tp3s[pfp.Tp3s.size() - 1 - ii];
      if(PosSep2(tp3.Pos, endPos) > 100) break;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) endDir[xyz] += tp3.Dir[xyz];
    } // ii
    SetMag(endDir, 1);
    pfp.Dir[1] = endDir;
    if(prt) mf::LogVerbatim("TC")<<"SBDFS: end direction "<<std::fixed<<std::setprecision(2)<<pfp.Dir[1][0]<<" "<<pfp.Dir[1][1]<<" "<<pfp.Dir[1][2];
    
  } // SortByDistanceFromStart
  
  /////////////////////////////////////////
  bool FitTp3(TjStuff& tjs, TrajPoint3& tp3, const std::vector<Tj2Pt>& tj2pts)
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
      auto& tp = tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      x0 += tjs.detprop->ConvertTicksToX(tp.Pos[1]/tjs.UnitsPerTick, planeID);
    }
    x0 /= (double)tj2pts.size();
    
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    double wght = 1;
    for(unsigned short ipt = 0; ipt < tj2pts.size(); ++ipt) {
      auto& tj2pt = tj2pts[ipt];
      auto& tp = tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      unsigned int cstat = planeID.Cryostat;
      unsigned int tpc = planeID.TPC;
      unsigned int plane = planeID.Plane;
      // get the wire plane offset
      double off = tjs.geom->WireCoordinate(0, 0, plane, tpc, cstat);
      // get the "cosine-like" component
      double cw = tjs.geom->WireCoordinate(1, 0, plane, tpc, cstat) - off;
      // the "sine-like" component
      double sw = tjs.geom->WireCoordinate(0, 1, plane, tpc, cstat) - off;
      double x = tjs.detprop->ConvertTicksToX(tp.Pos[1]/tjs.UnitsPerTick, planeID) - x0;
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
    tp3.ChiDOF = 1;
    
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
      tp3.ChiDOF = 10;
      return false;
    }
    
    tp3.Pos = fitPos;
    tp3.Dir = fitDir;

    return true;
  } // FitTp3
  
  /////////////////////////////////////////
  void MoveTp3ToZ(TjStuff& tjs, TrajPoint3& tp3, double z)
  {
    double dz = z - tp3.Pos[2];
    if(dz == 0) return;
    if(tp3.Dir[2] == 0) return;
    tp3.Pos[0] += dz * tp3.Dir[0] / tp3.Dir[2];
    tp3.Pos[1] += dz * tp3.Dir[1] / tp3.Dir[2];
    tp3.Pos[2] += dz;
  } // MoveTp3ToZ

  /////////////////////////////////////////
  bool CheckAndMerge(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Calculate the 3D-matching completeness of the set of Tjs in pfp.TjIDs.
    // Look for missing Tjs that are matched in 3D to the known set and merge them
    // if appropriate
    
    if(pfp.TjIDs.size() < 2) return false;
    if(tjs.Match3DCuts[0] <= 0) return false;
    // This function uses mallTraj but it isn't necessarily a failure if it doesn't exist
    if(tjs.mallTraj.size() < 6) return true;

    bool twoPlanes = (tjs.NumPlanes == 2);
    double yzcut = 1.5 * tjs.Match3DCuts[0];
    
    // list of Tjs that have a poor 3D match which should be ignored
    std::vector<unsigned short> ignore;
    
    for(unsigned short nit = 0; nit < 6; ++nit) {
      // put the pfp TjIDs into a vector of unsigned shorts for faster searching yyy
      std::vector<unsigned short> tjids;
      for(auto tjid : pfp.TjIDs) tjids.push_back((unsigned short)tjid);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" nit "<<nit;
        for(auto tjid : tjids) myprt<<" T"<<tjid;
      }
      // create a vector of bools for each tj for points that are matched in 3D 
      // This vector is for matches in 3 planes
      std::vector<std::vector<bool>> tjptMat3;
      // This vector is for matches in 2 planes
      std::vector<std::vector<bool>> tjptMat2;
      // and the plane index
      std::vector<unsigned short> tjplane;
      // Initialize the vectors
      for(unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
        auto& tj = tjs.allTraj[pfp.TjIDs[itj] - 1];
        std::vector<bool> tmp(tj.Pts.size(), false);
        tjptMat2.push_back(tmp);
        if(!twoPlanes) tjptMat3.push_back(tmp);
        tjplane.push_back(DecodeCTP(tj.CTP).Plane);
      } // tjid
      // count of triple matches on dead wires in the 3rd plane
      std::vector<unsigned short> deadCnt(pfp.TjIDs.size());
      // list of Tjs not in pfp.TjIDs that match in the 3rd plane
      std::vector<unsigned short> missTjIDs;
      std::vector<unsigned short> missTjCnt;
      pfp.Tp3s.clear();
      
      for(unsigned int ipt = 0; ipt < tjs.mallTraj.size() - 1; ++ipt) {
        auto& iTjPt = tjs.mallTraj[ipt];
        // ignore a tj
        if(std::find(ignore.begin(), ignore.end(), iTjPt.id) != ignore.end()) continue;
        unsigned short indx = 0;
        for(indx = 0; indx < tjids.size(); ++indx) if(iTjPt.id == tjids[indx]) break;
        // require that the Tj ID of this point be in the list
        if(indx == tjids.size()) continue;
        auto& itj = tjs.allTraj[iTjPt.id - 1];
        if(itj.AlgMod[kMat3D]) continue;
        auto& itp = itj.Pts[iTjPt.ipt];
        unsigned short iplane = DecodeCTP(itp.CTP).Plane;
        unsigned short tpc = DecodeCTP(itp.CTP).TPC;
        unsigned short cstat = DecodeCTP(itp.CTP).Cryostat;
        for(unsigned int jpt = ipt + 1; jpt < tjs.mallTraj.size() - 1; ++jpt) {
          auto& jTjPt = tjs.mallTraj[jpt];
          // ensure that the planes are different
          if(jTjPt.ctp == iTjPt.ctp) continue;
          // ignore a tj
          if(std::find(ignore.begin(), ignore.end(), jTjPt.id) != ignore.end()) continue;
          unsigned short jndx = 0;
          for(jndx = 0; jndx < tjids.size(); ++jndx) if(jTjPt.id == tjids[jndx]) break;
          // require that the Tj ID of this point be in the list
          if(jndx == tjids.size()) continue;
          // check for x range overlap. We know that jTjPt.xlo is > iTjPt.xlo because of the sort
          if(jTjPt.xlo > iTjPt.xhi) continue;
          // break out if the x range difference becomes large (5 cm)
          if(jTjPt.xlo > iTjPt.xhi + 5) break;
          auto& jtj = tjs.allTraj[jTjPt.id - 1];
          if(jtj.AlgMod[kMat3D]) continue;
          auto& jtp = jtj.Pts[jTjPt.ipt];
          TrajPoint3 ijtp3;
          if(!MakeTp3(tjs, itp, jtp, ijtp3)) continue;
          ijtp3.Tj2Pts.resize(2);
          ijtp3.Tj2Pts[0] = iTjPt;
          ijtp3.Tj2Pts[1] = jTjPt;
          // Set the 2-plane match bits
          tjptMat2[indx][iTjPt.ipt] = true;
          tjptMat2[jndx][jTjPt.ipt] = true;
          if(twoPlanes) continue;
          // count it as a triple if this point is in a dead region
          unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
          unsigned short kplane = 3 - iplane - jplane;
          unsigned int kwire = std::nearbyint(tjs.geom->WireCoordinate(ijtp3.Pos[1], ijtp3.Pos[2], kplane, tpc, cstat));
          if(kwire < tjs.WireHitRange[kplane].size() && tjs.WireHitRange[kplane][kwire].first == -1) {
            unsigned short kndx = 0;
            for(kndx = 0; kndx < tjplane.size(); ++kndx) if(kplane == tjplane[kndx]) break;
            if(kndx == deadCnt.size()) continue;
            ++deadCnt[kndx];
            // store the Tp3
            pfp.Tp3s.push_back(ijtp3);
            continue;
          } // dead wire in kplane
          for(unsigned int kpt = jpt + 1; kpt < tjs.mallTraj.size(); ++kpt) {
            auto& kTjPt = tjs.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
            // ignore a tj
            if(std::find(ignore.begin(), ignore.end(), kTjPt.id) != ignore.end()) continue;
            if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 5) break;
            auto& ktj = tjs.allTraj[kTjPt.id - 1];
            if(ktj.AlgMod[kMat3D]) continue;
            auto& ktp = ktj.Pts[kTjPt.ipt];
            TrajPoint3 iktp3;
            if(!MakeTp3(tjs, itp, ktp, iktp3)) continue;
            if(std::abs(ijtp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
            if(std::abs(ijtp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
            // make a copy of ijtp3 -> ijktp3
            auto ijktp3 = ijtp3;
            // add the Tj2Pt to it
            ijktp3.Tj2Pts.push_back(kTjPt);
            // Look for this point in pfp.TjIDs
            unsigned short kndx = 0;
            for(kndx = 0; kndx < tjids.size(); ++kndx) if(kTjPt.id == tjids[kndx]) break;
            if(kndx == tjids.size()) {
              // look for this tj in the missed tjs list
              unsigned short mndx = 0;
              for(mndx = 0; mndx < missTjIDs.size(); ++mndx) if(kTjPt.id == missTjIDs[mndx]) break;
              if(mndx == missTjIDs.size()) {
                // didn't find it so add it
                missTjIDs.push_back(kTjPt.id);
                missTjCnt.push_back(1);
              } else {
                ++missTjCnt[mndx];
              }
              // save the tp3 before continuing so we can check for missed Tjs
              pfp.Tp3s.push_back(ijktp3);
              continue;
            } // kndx check
            // Set the 3-plane match bits
            tjptMat3[indx][iTjPt.ipt] = true;
            tjptMat3[jndx][jTjPt.ipt] = true;
            tjptMat3[kndx][kTjPt.ipt] = true;
            // No sense refining the position of ijtp3 since it may be clobbered later
            pfp.Tp3s.push_back(ijktp3);
          } // kpt
        } // jpt
      } // ipt
      if(twoPlanes) {
        std::cout<<"CheckAndMerge: write some 2-plane code\n";
        return false;
      }
      // now count the number of tj points were matched      
      int pdgcode = PDGCodeVote(tjs, pfp.TjIDs, true);
      if(prt) mf::LogVerbatim("TC")<<"CheckAndMerge: pfp P"<<pfp.ID<<" pdgcode "<<pdgcode<<" matchVec index "<<pfp.BestPlane;
      // total number of points with charge in all Tjs
      float tnpwc = 0;
      // total number that are matched in 3D in 3 planes
      float tcnt3 = 0;
      // total number that are matched in 3D in 2 planes
      float tcnt2 = 0;
      unsigned short ignoreMe = USHRT_MAX;
      unsigned short nInShower = 0;
      // list of Tjs to add to pfp.TjIDs
      std::vector<int> addTjs;
      for(unsigned short itj = 0; itj < tjids.size(); ++itj) {
        auto& tj = tjs.allTraj[tjids[itj] - 1];
        // counts for each tj
        float npwc = 0;
        float cnt2 = 0;
        float cnt3 = 0;
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          if(tj.Pts[ipt].Chg <= 0) continue;
          ++npwc;
          if(tjptMat2[itj][ipt]) ++cnt2;
          if(tjptMat3[itj][ipt]) ++cnt3;
        } // ipt
        float completeness = cnt3 / npwc;
        tnpwc += npwc;
        tcnt3 += cnt3;
        tcnt2 += cnt2;
        if(prt) mf::LogVerbatim("TC")<<"T"<<tj.ID<<" npwc "<<npwc<<" cnt2 "<<cnt2<<" cnt3 "<<cnt3<<" deadCnt "<<deadCnt[itj]<<" completeness "<<completeness;
        if(npwc > 50 && completeness < 0.8 && pfp.BestPlane > 0) {
          // long Tj with poor completeness. This is probably a broken tj that isn't in the
          // pfp Tj list. Look for a match in the vicinity in tjs.matchVec
          unsigned int imsLo = pfp.BestPlane - 10;
          if(imsLo > tjs.matchVec.size() - 1) imsLo = 0;
          unsigned int imsHi = pfp.BestPlane + 10;
          if(imsHi > tjs.matchVec.size()) imsHi = tjs.matchVec.size();
          for(unsigned int ims = imsLo; ims < imsHi; ++ims) {
            auto& ms = tjs.matchVec[ims];
            // look for the ID in this matchvec entry
            if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tj.ID) == ms.TjIDs.end()) continue;
            // add the other tjs to pfp.TjIDs if they are available using the assumption
            // that these are broken Tjs
            for(auto btjid : ms.TjIDs) {
              if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), btjid) != pfp.TjIDs.end()) continue;
              auto& btj = tjs.allTraj[btjid - 1];
              if(btj.AlgMod[kMat3D]) continue;
              if(btj.AlgMod[kKilled]) continue;
              mf::LogVerbatim("TC")<<" Add T"<<btjid;
              addTjs.push_back(btj.ID);
            } // btjid
            break;
          } // ims
        } // long broken
        // Found a Tj that doesn't belong in this pfp. Prepare to ignore it
//        if(nit > 0 && completeness < 0.1) ignoreMe = itj;
        // count InShower tjs
        if(tj.AlgMod[kInShower] && completeness > 0.2) ++nInShower;
      } // itj
      float completeness = tcnt3 / tnpwc;
      if(prt) mf::LogVerbatim("TC")<<" completeness "<<completeness;
      // inspect the list of missed Tjs
      for(unsigned short mndx = 0; mndx < missTjIDs.size(); ++mndx) {
        auto& mtj = tjs.allTraj[missTjIDs[mndx] - 1];
        unsigned short mtjid = mtj.ID;
        // count the number of points matched to this pfp
        std::vector<bool> mtjMat3(mtj.Pts.size(), false);
        for(auto& tp3 : pfp.Tp3s) {
          for(auto& tj2pt : tp3.Tj2Pts) {
            if(tj2pt.id == mtjid) mtjMat3[tj2pt.ipt] = true;
          } // tj2pt
        } // tp3
        // now do the count
        float npwc = 0;
        float cnt3 = 0;
        for(unsigned short ipt = mtj.EndPt[0]; ipt <= mtj.EndPt[1]; ++ipt) {
          auto& mtp = mtj.Pts[ipt];
          if(mtp.Chg <= 0) continue;
          ++npwc;
          if(mtjMat3[ipt]) ++cnt3;
        } // ipt
        if(npwc == 0) continue;
        float missFrac = cnt3 / npwc;
        // ensure that we don't merge tjs that are incompatible with each other
        bool doMerge = true;
        // don't merge delta rays with non-delta rays
        if(pdgcode != 11 && mtj.AlgMod[kDeltaRay]) doMerge = false;
        // don't merge InShower Tjs with not-InShower
        if(nInShower == 0 && mtj.AlgMod[kInShower]) doMerge = false;
        if(nInShower > 1 && !mtj.AlgMod[kInShower]) doMerge = false;
        if(prt) mf::LogVerbatim("TC")<<" miss T"<<missTjIDs[mndx]<<" count "<<missTjCnt[mndx]<<" missFrac "<<missFrac<<" pdgcode "<<mtj.PDGCode<<" inShower? "<<mtj.AlgMod[kInShower]<<" doMerge? "<<doMerge;
        if(missFrac < 0.1 || !doMerge) {
          // Add this tj to the ignore list
          ignore.push_back(mtjid);
          continue;
        } 
        // found a missing tj. Add it to the list
        addTjs.push_back(mtj.ID);
      } // mndx
      // put the completeness into EffPur
      pfp.EffPur = completeness;
      // drop a Tj and ignore it?
      if(ignoreMe < pfp.TjIDs.size()) {
        ignore.push_back(pfp.TjIDs[ignoreMe]);
        pfp.TjIDs.erase(pfp.TjIDs.begin() + ignoreMe);
      } else {
        // break out of the loop?
        if(completeness > 0.95) break;
        if(addTjs.empty()) break;
      }
      if(!addTjs.empty()) pfp.TjIDs.insert(pfp.TjIDs.end(), addTjs.begin(), addTjs.end());
    } // nit
    
    // re-build the list to ensure that all Tjs are indeed valid
    std::vector<int> tmp;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kMat3D]) continue;
      tmp.push_back(tjid);
    } // tjid
    pfp.TjIDs = tmp;
    
    // something bad happened.
    if(pfp.TjIDs.size() < 2) return false;
    
    // Set the starting position in 3D if it isn't already defined by a 3D vertex
    SetNewStart(tjs, pfp, prt);
    // Sort Tp3s by distance from the start
    SortByDistanceFromStart(tjs, pfp, prt);
    // Check the list of Tjs and merge those that are in the same plane
    return MergePFPTjs(tjs, pfp, prt);

  } // CheckAndMerge
  
  /////////////////////////////////////////
  bool MergePFPTjs(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Checks the list of Tjs in pfp.TjIDs and merges those that are in 
    // the same plane. This function uses the ordering of Tps which should
    // have been sorted
    if(pfp.TjIDs.empty()) return false;
    if(pfp.Tp3s.empty()) return false;
    
    geo::TPCGeo const& TPC = tjs.geom->TPC(pfp.TPCID);
    unsigned short nplanes = TPC.Nplanes();
    // vector of tj IDs that will replace pfp.TjIDs
    std::vector<int> newTjIDs;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MergePFPTjs: in";
      for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
    }
    
    for(unsigned short plane = 0; plane < nplanes; ++plane) {
      std::vector<unsigned short> tjids;
      for(auto tjid : pfp.TjIDs) if(DecodeCTP(tjs.allTraj[tjid - 1].CTP).Plane == plane) tjids.push_back((unsigned short)tjid);
      // No Tjs in this plane or just one
      if(tjids.size() == 1) newTjIDs.push_back((int)tjids[0]);
      if(tjids.size() < 2) continue;
      // Need to merge. create a new Tj
      Trajectory mtj;
      // give it a bogus ID
      mtj.ID = -6666;
      mtj.CTP = tjs.allTraj[tjids[0] - 1].CTP;
      mtj.StepDir = 1;
      // create vectors to ensure that TPs are only used once
      std::vector<std::vector<bool>> ptused;
      for(unsigned short ii = 0; ii < tjids.size(); ++ii) {
        auto& tj = tjs.allTraj[tjids[ii] - 1];
        std::vector<bool> tmp(tj.Pts.size(), false);
        ptused.push_back(tmp);
      } // ii
      // iterate over the Tp3s and add TPs to it that are in the list of tjids
      // Check the ends of the Tjs for a 2D vertex while we are doing this
      for(unsigned short itp3 = 0; itp3 < pfp.Tp3s.size(); ++itp3) {
        auto& tp3 = pfp.Tp3s[itp3];
        for(auto& tj2pt : tp3.Tj2Pts) {
          unsigned short tjIndx = 0;
          for(tjIndx = 0; tjIndx < tjids.size(); ++tjIndx) if(tj2pt.id == tjids[tjIndx]) break;
          if(tjIndx == tjids.size()) continue;
          // see if the point is used
          if(ptused[tjIndx][tj2pt.ipt]) continue;
          // reference to the broken Tj
          auto& btj = tjs.allTraj[tjids[tjIndx] - 1];
          auto mtp = btj.Pts[tj2pt.ipt];
          // check for a vertex at the end closest to this point on the broken tj
          // and attach it to the merged tj
          unsigned short btjEnd = 0;
          if(tj2pt.ipt > 0.5 * (btj.EndPt[0] + btj.EndPt[1])) btjEnd = 1;
          unsigned short mtjEnd = 0;
          if(itp3 > 0.5 * pfp.Tp3s.size()) mtjEnd = 1;
          mtj.VtxID[mtjEnd] = btj.VtxID[btjEnd];;
          mtj.Pts.push_back(mtp);
          ptused[tjIndx][tj2pt.ipt] = true;
        } // tj2pt
      } // tp3
      // check for a serious error
      if(mtj.Pts.size() < 3) return false;
      mtj.AlgMod[kMat3DMerge] = true;
      // kill the broken tjs
      for(auto tjid : tjids) MakeTrajectoryObsolete(tjs, tjid - 1);
      // save the new one
      if(!StoreTraj(tjs, mtj)) {
        std::cout<<"MergePFPTjs: StoreTraj failed\n";
        return false;
      }
      int newTjID = tjs.allTraj.size();
      newTjIDs.push_back(newTjID);
    } // plane
    pfp.TjIDs = newTjIDs;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MergePFPTjs: out";
      for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
    }
    return true;
  } // MergePFPTjs
  
  /////////////////////////////////////////
  void FindXMatches(TjStuff& tjs, unsigned short numPlanes, short maxScore, std::vector<MatchStruct>& matVec, bool prt)
  {
    // This function matches trajectory points in tjs.mallTraj using the X position. These points should
    // have already been sorted by increasing X by the function that created mallTraj. 
    
    if(tjs.mallTraj.empty()) return;
    if(tjs.Match3DCuts[0] <= 0) return;
    if(numPlanes < 2) return;
    
    int cstat = DecodeCTP(tjs.mallTraj[0].ctp).Cryostat;
    int tpc = DecodeCTP(tjs.mallTraj[0].ctp).TPC;
    constexpr float twopi = 2 * M_PI;
    constexpr float piOver2 = M_PI / 2;
    
    // create a temp vector to check for duplicates
    auto inMatVec = matVec;
    std::vector<MatchStruct> temp;
    
    // the minimum number of points for matching
    unsigned short minPts = 2;
    // override this with the user minimum for 2-plane matches
    if(numPlanes == 2) minPts = tjs.Match3DCuts[2];
    
    // max number of match combos left
    unsigned int nAvailable = 0;
    if(matVec.size() < tjs.Match3DCuts[4]) nAvailable = tjs.Match3DCuts[4] - matVec.size();
    if(nAvailable == 0 || nAvailable > tjs.Match3DCuts[4]) return;
    
    // these cuts presume that the resolution in X is better than it is in Y and Z
    float xcut = tjs.Match3DCuts[0];
    double yzcut = 1.5 * xcut;
    for(unsigned int ipt = 0; ipt < tjs.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = tjs.mallTraj[ipt];
      // length cut
      if(iTjPt.npts < minPts) continue;
      // look for matches using Tjs that have the correct score
      if(iTjPt.score < 0 || iTjPt.score > maxScore) continue;
      auto& itp = tjs.allTraj[iTjPt.id - 1].Pts[iTjPt.ipt];
      unsigned short iplane = DecodeCTP(itp.CTP).Plane;
      for(unsigned int jpt = ipt + 1; jpt < tjs.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = tjs.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.ctp == iTjPt.ctp) continue;
        // length cut
        if(jTjPt.npts < minPts) continue;
        if(jTjPt.score < 0 || jTjPt.score > maxScore) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtp = tjs.allTraj[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
        TrajPoint3 tp3;
        if(!MakeTp3(tjs, itp, jtp, tp3)) continue;
//        bool dijOK = (useAngle && iTjPt.npts > 5 && jTjPt.npts > 5);
        // count weight is one for a two-plane match
        float cntWght = 1;
        if(numPlanes == 3) {
          // numPlanes == 3
          for(unsigned int kpt = jpt + 1; kpt < tjs.mallTraj.size(); ++kpt) {
            auto& kTjPt = tjs.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
            if(kTjPt.score < 0 || kTjPt.score > maxScore) continue;
            if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 5) break;
            auto& ktp = tjs.allTraj[kTjPt.id - 1].Pts[kTjPt.ipt];
            unsigned short kplane = DecodeCTP(ktp.CTP).Plane;
            TrajPoint3 iktp3;
            if(!MakeTp3(tjs, itp, ktp, iktp3)) continue;
            if(std::abs(tp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
            if(std::abs(tp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
            float dang = 0;
            if(tjs.Match3DCuts[1] > 0) {
              dang = std::abs(DeltaAngle(tp3.Dir, iktp3.Dir));
              while(dang >  M_PI) dang -= twopi;
              if(dang >  piOver2) dang = M_PI - dang;
              float mcsmom = tjs.allTraj[iTjPt.id - 1].MCSMom + tjs.allTraj[jTjPt.id - 1].MCSMom + tjs.allTraj[kTjPt.id - 1].MCSMom;
              mcsmom /= 3;
              //              mf::LogVerbatim("TC")<<"dang "<<std::fixed<<std::setprecision(3)<<dang<<std::setprecision(2)<<" "<<(iktp3.Pos[1] - tp3.Pos[1])<<" "<<(iktp3.Pos[2] - tp3.Pos[2])<<" "<<(int)mcsmom<<" "<<ijkMatch<<" "<<TMeV;
              if(mcsmom > 150 && dang > tjs.Match3DCuts[1]) continue;
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
//              ++ms.Count;
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
//              ms.Count = 1;
              temp.push_back(ms);
              // give up if there are too many
              if(temp.size() > nAvailable) break;
            } // not found in the list
          } // kpt
          // numPlanes == 3
        } else {
          // 2-plane TPC or 2-plane match in a 3-plane TPC
          if(tjs.NumPlanes == 3) {
            // See if there is a signal at this point.
            unsigned short kplane = 3 - iplane - jplane;
            float fkwire = tjs.geom->WireCoordinate(tp3.Pos[1], tp3.Pos[2], kplane, tpc, cstat);
//            float fkwire = tjs.geom->WireCoordinate(jyp, jzp, kpl, tpc, cstat);
            if(fkwire < 0 || fkwire > tjs.MaxPos0[kplane]) continue;
            TrajPoint tpk;
            tpk.CTP = EncodeCTP(cstat, tpc, kplane);
            tpk.Pos[0] = fkwire;
            float xp = 0.5 * (iTjPt.xlo + iTjPt.xhi);
            tpk.Pos[1] = tjs.detprop->ConvertXToTicks(xp, kplane, tpc, cstat) * tjs.UnitsPerTick;
            // Note that SignalAtTp assumes that a signal exists if the wire is dead
            if(!SignalAtTp(tjs, tpk)) continue;
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
    
    // Fill MatchFrac = number of 3D matched points divided by the length of the
    // longest trajectory. Note that this can be larger than 1.
    for(auto& ms : temp) {
      if(ms.Count == 0) continue;
      float maxlen = 1;
      for(auto& tjID : ms.TjIDs) {
        float len = NumUsedHitsInTj(tjs, tjs.allTraj[tjID-1]);
        if(len > maxlen) maxlen = len;
      }
      ms.MatchFrac = ms.Count / maxlen;
    } // ms
    
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
  bool MakeTp3(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3)
  {
    // Make a 3D trajectory point using two 2D trajectory points
    tp3.Dir = {999};
    tp3.Pos = {999};
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);
    double upt = tjs.UnitsPerTick;
    double ix = tjs.detprop->ConvertTicksToX(itp.Pos[1] / upt, iPlnID);
    double jx = tjs.detprop->ConvertTicksToX(jtp.Pos[1] / upt, jPlnID);
    
    // don't continue if the points are wildly far apart in X
    if(std::abs(ix - jx) > 10) return false;
    tp3.Pos[0] = 0.5 * (ix + jx);
    // determine the wire orientation and offsets using WireCoordinate
    // wire = yp * OrthY + zp * OrthZ - Wire0 = cs * yp + sn * zp - wire0
    // wire offset
    double iw0 = tjs.geom->WireCoordinate(0, 0, iPlnID);
    // cosine-like component
    double ics = tjs.geom->WireCoordinate(1, 0, iPlnID) - iw0;
    // sine-like component
    double isn = tjs.geom->WireCoordinate(0, 1, iPlnID) - iw0;
    double jw0 = tjs.geom->WireCoordinate(0, 0, jPlnID);
    double jcs = tjs.geom->WireCoordinate(1, 0, jPlnID) - jw0;
    double jsn = tjs.geom->WireCoordinate(0, 1, jPlnID) - jw0;
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
    
    // make a copy of itp and shift it by many wires to avoid precision problems
    double itp2_0 = itp.Pos[0] + 100;
    double itp2_1 = itp.Pos[1];
    if(std::abs(itp.Dir[0]) > 0.01) itp2_1 += 100 * itp.Dir[1] / itp.Dir[0];
    // Create a second Point3 for the shifted point
    Point3_t pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = tjs.detprop->ConvertTicksToX(itp2_1 / upt, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = tjs.detprop->ConvertXToTicks(pos2[0], jPlnID) * upt;
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
  
  ////////////////////////////////////////////////
  void FilldEdx(TjStuff& tjs, TrajPoint3& tp3)
  {
    // fills the dEdx variables in tp3 (MeV/cm)
    tp3.dEdx = 0;
    tp3.dEdxErr = 0.5;
    if(!tp3.IsValid) return;
    if(tp3.Tj2Pts.empty()) return;
    double t0 = 0;
    float sum = 0;
    float sum2 = 0;
    float cnt = 0;
    for(auto& tj2pt : tp3.Tj2Pts) {
      auto& tp = tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt];
      if(tp.Chg <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      double angleToVert = tjs.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
      double cosgamma = std::abs(std::sin(angleToVert) * tp3.Dir[1] + std::cos(angleToVert) * tp3.Dir[2]);
      if(cosgamma == 0) continue;
      double dx = tjs.geom->WirePitch(planeID) / cosgamma;
      // sum the hit charge
      double chg = 0;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        chg += tjs.fHits[iht].Integral;
      }
      double dQ = chg / dx;
      double time = tp.Pos[1] / tjs.UnitsPerTick;
      float dedx = tjs.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
      if(dedx > 200) continue;
      sum += dedx;
      sum2 += dedx * dedx;
      ++cnt;
    } // tj2pt
    if(cnt == 0) return;
    tp3.dEdx = sum / cnt;
    if(cnt > 1) {
      float arg = sum2 - cnt * tp3.dEdx * tp3.dEdx;
      if(arg > 0) {
        tp3.dEdxErr = sqrt(arg / (cnt - 1));
        // convert to error on the mean
        tp3.dEdxErr /= sqrt(cnt);
      }
    } // cnt > 1
  } // FilldEdx
  
  ////////////////////////////////////////////////
  bool Split3DKink(TjStuff& tjs, PFPStruct& pfp, double sep, bool prt)
  {
    // Finds kinks in the PFParticle, splits Tjs, creates 2D vertices and forces a rebuild if any are found
    if(pfp.Tp3s.empty()) return false;
    if(pfp.DirectionFixed) return false;
    if(!tjs.UseAlg[kSplit3DKink]) return false;
    
    auto kinkPts = FindKinks(tjs, pfp, sep, prt);
    if(kinkPts.empty()) return false;
    if(prt) mf::LogVerbatim("TC")<<"Split3DKink found a kink at Tp3s point "<<kinkPts[0];
    
    // Only split the biggest angle kink
    double big = 0;
    unsigned short kpt = 0;
    for(auto ipt : kinkPts) {
      double dang = KinkAngle(tjs, pfp.Tp3s, ipt, sep);
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
        auto& tj = tjs.allTraj[tp2.id - 1];
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
          vx2.ID = tjs.vtx.size() + 1;
          vx2.CTP = tj.CTP;
          vx2.Topo = 10;
          vx2.Pos = tp.Pos;
          if(!StoreVertex(tjs, vx2)) return false;
          if(!SplitTraj(tjs, tp2.id - 1, tp2.ipt, tjs.vtx.size() - 1, prt)) return false;
          vx2ids.push_back(vx2.ID);
          AttachAnyTrajToVertex(tjs, tjs.vtx.size() - 1, prt);
          if(prt) mf::LogVerbatim("TC")<<" tj "<<tj.ID<<" new 2V"<<vx2.ID;
          tjs.NeedsRebuild = true;
        }
      } // tp2
    } // ipt
    
    if(vx2ids.size() != tjs.NumPlanes) {
//      std::cout<<"Split3DKink: TODO pfp "<<pfp.ID<<" only has "<<vx2ids.size()<<" 2D vertices. \n";
      return false;
    }
    Vtx3Store vx3;
    vx3.TPCID = pfp.TPCID;
    vx3.ID = tjs.vtx3.size() + 1;
    vx3.X = pfp.Tp3s[kpt].Pos[0];
    vx3.Y = pfp.Tp3s[kpt].Pos[1];
    vx3.Z = pfp.Tp3s[kpt].Pos[2];
    for(auto vx2id : vx2ids) {
      if(vx2id == 0) continue;
      auto& vx2 = tjs.vtx[vx2id - 1];
      unsigned short plane = DecodeCTP(vx2.CTP).Plane;
      vx3.Vx2ID[plane] = vx2id;
      vx2.Vx3ID = vx3.ID;
    } // vx2id
    std::cout<<"Split3DKink add 3V"<<vx3.ID<<"\n";
    tjs.vtx3.push_back(vx3);
    // mark this as needing an update
    pfp.NeedsUpdate = true;
    return true;
  } // Split3DKink
  
  ////////////////////////////////////////////////
  std::vector<unsigned short> FindKinks(const TjStuff& tjs, PFPStruct& pfp, double sep, bool prt)
  {
    // returns a vector of indices in pfp.Tp3s where kinks exist. The kink angle is calculated using
    // Tp3s separated by +/- sep (cm)
    std::vector<unsigned short> kinkPts;
//    double kang = tjs.KinkCuts[0];
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
      double dang = KinkAngle(tjs, pfp.Tp3s, ipt, sep);
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
  double KinkAngle(const TjStuff& tjs, const std::vector<TrajPoint3>& tp3s, unsigned short atPt, double sep)
  {
    // calculate a kink angle at the TjPt 
    if(tp3s.empty()) return -1;
    if(atPt < 1 || atPt > tp3s.size() - 2) return -1;
    double sep2 = sep * sep;
    unsigned short pt1 = USHRT_MAX;
    for(unsigned short ii = 1; ii < tp3s.size(); ++ii) {
      unsigned short ipt = atPt - ii;
      if(tp3s[ipt].IsValid && PosSep2(tp3s[atPt].Pos, tp3s[ipt].Pos) > sep2) {
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
      if(tp3s[ipt].IsValid && PosSep2(tp3s[atPt].Pos, tp3s[ipt].Pos) > sep2) {
        pt2 = ipt;
        break;
      }
    } // ii
    if(pt2 == USHRT_MAX) return -1;
    return DeltaAngle(tp3s[pt1].Dir, tp3s[pt2].Dir);
  } // KinkAngle

  ////////////////////////////////////////////////
  PFPStruct CreatePFP(const TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // The calling function should define the size of pfp.TjIDs
    PFPStruct pfp;
    pfp.ID = tjs.pfps.size() + 1;
    // assume it is it's own parent
    pfp.ParentID = pfp.ID;
    pfp.TPCID = tpcid;
    // initialize arrays for both ends
    pfp.dEdx[0].resize(tjs.NumPlanes, 0);
    pfp.dEdx[1].resize(tjs.NumPlanes, 0);
    pfp.dEdxErr[0].resize(tjs.NumPlanes, 0);
    pfp.dEdxErr[1].resize(tjs.NumPlanes, 0);
    for(unsigned short startend = 0; startend < 2; ++startend) {
      pfp.Dir[startend] = {0, 0, 0};
      pfp.DirErr[startend] = {0, 0, 0};
      pfp.XYZ[startend] = {0, 0, 0};
    }
    return pfp;
  } // CreatePFP
  
  //////////////////////////////////////////
  void FindPFParticles(std::string fcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // Match Tjs in 3D and create PFParticles
    
    if(tjs.Match3DCuts[0] <= 0) return;
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" called FindPFParticles";
    // clear matchVec
    tjs.matchVec.clear();
    // assume that this will not have to be re-done
    tjs.NeedsRebuild = false;
    
    // Match these points in 3D and put the results in tjs.matchVec
    std::vector<MatchStruct> matVec;
    // first look for 3-plane matches in a 3-plane TPC
    if(tjs.NumPlanes == 3) {
      // Match Tjs with high quality vertices first and the leftovers next
      for(short maxScore = 0; maxScore < 2; ++maxScore) FindXMatches(tjs, 3, maxScore, matVec, prt);
    } // 3-plane TPC
    // Make 2-plane matches if we haven't hit the user-defined size limit
    if(matVec.size() < tjs.Match3DCuts[4]) {
      // 2-plane TPC or 2-plane matches in a 3-plane TPC
      if(tjs.NumPlanes == 2) {
        for(short maxScore = 0; maxScore < 2; ++maxScore) FindXMatches(tjs, 2, maxScore, matVec, prt);
      } else {
        // Make one attempt at 2-plane matches in a 3-plane TPC, setting maxScore large
        FindXMatches(tjs, 2, 3, matVec, prt);
      }
    } // can add more combinations
//    if(matVec.size() >= tjs.Match3DCuts[4]) std::cout<<"FMV: Hit the max combo limit "<<matVec.size()<<" events processed "<<tjs.EventsProcessed<<"\n";
    
    // sort by decreasing match count
    if(matVec.size() > 1) {
      std::vector<SortEntry> sortVec(matVec.size());
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        sortVec[ii].index = ii;
        sortVec[ii].val = matVec[ii].Count;
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
      std::vector<MatchStruct> tmpVec;
      tmpVec.reserve(matVec.size());
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        tmpVec.push_back(matVec[sortVec[ii].index]);
      } // ii
      matVec = tmpVec;
    } // sort matVec
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FMV: matVec\n";
      unsigned short cnt = 0;
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        auto& ms = matVec[ii];
        if(ms.Count == 0) continue;
        myprt<<ii<<" Count "<<ms.Count<<" TjIDs:";
        for(auto& tjID : ms.TjIDs) myprt<<" T"<<std::to_string(tjID);
        myprt<<" NumUsedHitsInTj ";
        for(auto& tjID : ms.TjIDs) myprt<<" "<<NumUsedHitsInTj(tjs, tjs.allTraj[tjID-1]);
        myprt<<" MatchFrac "<<std::fixed<<std::setprecision(2)<<ms.MatchFrac;
        myprt<<" TjChgAsymmetry "<<std::fixed<<std::setprecision(2)<<MaxChargeAsymmetry(tjs, ms.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(tjs, ms.TjIDs, false);
        myprt<<"\n";
        ++cnt;
        if(cnt == 1000) {
          myprt<<"...stopped printing after 500 entries.";
          break;
        }
      } // ii
    } // prt
    
    // put the maybe OK matches into tjs
    for(auto& ms : matVec) {
      if(ms.Count < 2) continue;
      // require that at least 20% of the hits are matched in the longest Tj. Note that MatchFrac may be > 1
      // in particular for small angle trajectories
      if(ms.MatchFrac < 0.2) continue;
      if(MaxChargeAsymmetry(tjs, ms.TjIDs) > tjs.Match3DCuts[5]) continue;
      // check for duplicates
      bool skipit = false;
      for(auto& oms : tjs.matchVec) {
        if(ms.TjIDs == oms.TjIDs) {
          skipit = true;
          break;
        }
      } // oms
      if(skipit) continue;
      tjs.matchVec.push_back(ms);
    }
    if(tjs.matchVec.empty()) return;
    
    // create the list of associations to matches that will be converted to PFParticles
    // Start with Tjs attached to 3D vertices. This is only done when reconstructing neutrino events
    if(!tjs.TestBeam) {
      Match3DVtxTjs(tjs, tpcid, prt);
      // Re-do the Tj hierarchy and re-find PFParticles if Match3DVtxTjs merged/split Tjs or
      // added/removed vertices
      if(tjs.NeedsRebuild) return;
    }
    // Re-check matchVec with the user cut matchfrac cut to eliminate poor combinations
    for(unsigned int indx = 0; indx < tjs.matchVec.size(); ++indx) {
      auto& ms = tjs.matchVec[indx];
      if(ms.Count == 0) continue;
      // check for a minimum user-defined match fraction
      if(ms.MatchFrac < tjs.Match3DCuts[3]) ms.Count = 0;
    } // ms
    // define the PFParticleList
    for(unsigned int indx = 0; indx < tjs.matchVec.size(); ++indx) {
      auto& ms = tjs.matchVec[indx];
      // ignore dead matches
      if(ms.Count == 0) continue;
      // skip this match if any of the trajectories is already matched or merged and killed
      bool skipit = false;
      for(auto tjID : ms.TjIDs) {
        auto& tj = tjs.allTraj[tjID - 1];
        if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) skipit = true;
      } // tjID
      if(skipit) continue;
      // count the number of shower Tjs
      unsigned short nstj = 0;
      for(unsigned short ipl = 0; ipl < ms.TjIDs.size(); ++ipl) {
        unsigned short itj = ms.TjIDs[ipl] - 1;
        if(tjs.allTraj[itj].AlgMod[kMat3D]) skipit = true;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) ++nstj;
      }
      if(skipit) continue;
      // Require 0 or a matched shower Tj in all planes
      if(nstj != 0 && nstj != ms.TjIDs.size()) continue;
      PFPStruct pfp = CreatePFP(tjs, tpcid);
      pfp.TjIDs = ms.TjIDs;
      // re-purpose BestPlane so we can search nearby entries of tjs.matchVec
      pfp.BestPlane = indx;
      // do a first look for broken tjs in the lower-rank entries
      for(unsigned int jj = indx + 1; jj < indx + 11; ++indx) {
        if(jj == tjs.matchVec.size()) break;
        auto& jms = tjs.matchVec[jj];
        auto shared = SetIntersection(jms.TjIDs, ms.TjIDs);
        if(shared.size() < 2) continue;
        for(auto tjid : jms.TjIDs) {
          if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjid) != pfp.TjIDs.end()) continue;
          pfp.TjIDs.push_back(tjid);
        }
      } // jj
      if(!DefinePFP("FPFP", tjs, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" DefinePFP failed";
        continue;
      }
      double sep = 1;
      Split3DKink(tjs, pfp, sep, prt);
      AnalyzePFP(tjs, pfp, prt);
      pfp.PDGCode = PDGCodeVote(tjs, pfp.TjIDs, prt);
      if(!StorePFP(tjs, pfp)) {
        if(prt) mf::LogVerbatim("TC")<<" StorePFP failed "<<pfp.ID;
      }
    } // indx
    //    CheckNoMatchTjs(tjs, tpcid, prt);
    
  } // FindPFParticles
  
  /////////////////////////////////////////
  bool DefinePFP(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // This function is called after the 3D matched TjIDs have been specified and optionally
    // a start or end vertex ID. It defines the PFParticle but doesn't store it
    
    if(pfp.PDGCode == 1111) return false;
    if(pfp.TjIDs.size() < 2) return false;
    
    std::string fcnLabel = inFcnLabel + ".DPFP";

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" Vx3ID "<<pfp.Vx3ID[0]<<" "<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
    }

    if(pfp.Vx3ID[0] == 0 && pfp.Vx3ID[1] > 0) {
      std::cout<<fcnLabel<<" pfp P"<<pfp.ID<<" end 1 has a vertex but end 0 doesn't. No endpoints defined\n";
      return false;
    }
    
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kMat3D]) {
        std::cout<<fcnLabel<<" pfp "<<pfp.ID<<" uses tj "<<tj.ID<<" but kMat3D is set true\n";
        return false;
      }
    } // tjid
    
    // check the completeness of matching points in this set of Tjs and possibly
    // merge Tjs
    if(!CheckAndMerge(tjs, pfp, prt)) return false;

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" Vx3ID "<<pfp.Vx3ID[0]<<" "<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
    }
    
    // Use the ends of the pfp Tjs to find small angle tracks which have a
    // small X range in which case the DirectionFixed variable is set true. The
    // variables pfp.Dir[end] are defined.
    // The standard X-matching scheme doesn't work for this situation.
    // SetEndPoints will determine the start and end points (pfp.XYZ[end]) in this case.
    SetEndPoints(tjs, pfp, prt);
//    MakePFPTp3s(tjs, pfp, false);
    pfp.NeedsUpdate = false;
    return true;
  } // DefinePFP
  
  /////////////////////////////////////////
  void AnalyzePFP(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Analyzes the PFP for oddities and tries to fix them
    if(pfp.ID == 0) return;
    if(pfp.TjIDs.empty()) return;
    
    // don't bother analyzing this pfp has been altered
    if(pfp.NeedsUpdate) {
      if(prt) mf::LogVerbatim("TC")<<"AnalyzePFP: "<<pfp.ID<<" needs to be updated. Skip analysis ";
      return;
    }
    if(prt) mf::LogVerbatim("TC")<<"inside AnalyzePFP "<<pfp.ID<<" NeedsUpdate? "<<pfp.NeedsUpdate<<" NeedsRebuild? "<<tjs.NeedsRebuild;
    
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
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"APFP: Tjs in pfp\n";
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        myprt<<tj.ID<<" npwc "<<NumPtsWithCharge(tjs, tj, false);
        unsigned short indx = 0;
        for(indx = 0; indx < tjIDs.size(); ++indx) if(tjIDs[indx] == tjid) break;
        if(indx == tjIDs.size()) {
          myprt<<" not found in pfp "<<pfp.ID<<"\n";
          continue;
        }
        myprt<<" nTp3 "<<tjCnt[indx]<<"\n";
      } // tjid
      myprt<<" Tjs not in pfp:";
      auto different = SetDifference(pfp.TjIDs, tjIDs);
      if(different.empty()) {
        myprt<<" none\n";
      } else {
        for(auto tjid : different) myprt<<" "<<tjid;
      }
    } // prt
    
  } // AnalyzePFP
  
  /////////////////////////////////////////
  void PFPVertexCheck(TjStuff& tjs)
  {
    // Ensure that all PFParticles have a start vertex. It is possible for
    // PFParticles to be attached to a 3D vertex that is later killed.
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.Vx3ID[0] > 0) continue;
      if(pfp.Vx3ID[1] == 0 && !pfp.Tp3s.empty()) {
        // See if the direction needs to be changed
        SetNewStart(tjs, pfp, false);
        SortByDistanceFromStart(tjs, pfp, false);      }
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = true;
      tjs.vtx3.push_back(vx3);
      std::cout<<"PFPVertexCheck: add 3V"<<vx3.ID<<"\n";
      pfp.Vx3ID[0] = vx3.ID;
    } // pfp
  } // PFPVertexCheck
  
  /////////////////////////////////////////
  void DefinePFPParents(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    /*
     This function reconciles vertices, PFParticles and Tjs, then
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
    if(tjs.pfps.empty()) return;
    
    int neutrinoPFPID = 0;
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != tpcid) continue;
      if(!tjs.TestBeam && neutrinoPFPID == 0 && (pfp.PDGCode == 12 || pfp.PDGCode == 14)) neutrinoPFPID = pfp.ID;
      if(pfp.Vx3ID[0] > 0) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = true;
      // TODO: we need to have PFP track position errors defined 
      unsigned short mergeToVx3ID = IsCloseToVertex(tjs, vx3);
      if(mergeToVx3ID > 0) {
        if(prt) mf::LogVerbatim("TC")<<"Merge PFP vertex "<<vx3.ID<<" with existing 3V"<<mergeToVx3ID;
        if(!AttachPFPToVertex(tjs, pfp, 0, mergeToVx3ID, prt)) {
          if(prt) mf::LogVerbatim("TC")<<" Failed to attach pfp "<<pfp.ID<<". Make new vertex \n";
          mergeToVx3ID = 0;
        }
      } // mergeMe > 0
      if(mergeToVx3ID == 0) {
        // Add the new vertex and attach the PFP to it
        tjs.vtx3.push_back(vx3);
        if(!AttachPFPToVertex(tjs, pfp, 0, vx3.ID, prt)) {
          if(prt) mf::LogVerbatim("TC")<<"Merge PFP vertex 3V"<<vx3.ID<<" with new vtx 3V"<<mergeToVx3ID;
        }
      } // merge to new vertex
    } // pfp
    
    // define the end vertex if the Tjs have end vertices
    constexpr unsigned short end1 = 1;
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != tpcid) continue;
      // already done?
      if(pfp.Vx3ID[end1] > 0) continue;
      // count 2D -> 3D matched vertices
      unsigned short cnt3 = 0;
      unsigned short vx3id = 0;
      // list of unmatched 2D vertices that should be merged
      std::vector<unsigned short> vx2ids;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.VtxID[end1] == 0) continue;
        auto& vx2 = tjs.vtx[tj.VtxID[end1] - 1];
        if(vx2.Vx3ID == 0) {
          if(vx2.Topo == 1 && vx2.NTraj == 2) vx2ids.push_back(vx2.ID);
          continue;
        }
        if(vx3id == 0) vx3id = vx2.Vx3ID;
        if(vx2.Vx3ID == vx3id) ++cnt3;
      } // tjid
      if(cnt3 > 1) {
        pfp.Vx3ID[end1] = vx3id;
        if(cnt3 != tjs.NumPlanes) mf::LogVerbatim("TC")<<"DPFPR: Missed an end vertex for PFP "<<pfp.ID<<" Write some code";
      }
    } // pfp
    
    // Assign a PDGCode to each PFParticle and look for a parent
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != tpcid) continue;
      // skip a neutrino PFParticle
      if(pfp.PDGCode == 12 || pfp.PDGCode == 14) continue;
      pfp.PDGCode = PDGCodeVote(tjs, pfp.TjIDs, prt);
      // next look for a parent
      int pfpParentID = INT_MAX;
      unsigned short nParent = 0;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.ParentID == tj.ID) continue;
        unsigned short ppindex = GetPFPIndex(tjs, tj.ParentID);
        if(ppindex == USHRT_MAX) continue;
        int ppid = ppindex + 1;
        if(pfpParentID == INT_MAX) pfpParentID = ppid;
        if(ppid == pfpParentID) ++nParent;
      } // ii
      // look for a parent
      if(nParent > 1) {
        pfp.ParentID = (size_t)pfpParentID;
        auto& parpfp = tjs.pfps[pfpParentID - 1];
        parpfp.DtrIDs.push_back(pfp.ID);
      } // nParent > 1
    } // ipfp

    
    if(tjs.TestBeam) {
      DefinePFPParentsTestBeam(tjs, tpcid, prt);
      return;
    }

    // associate primary PFParticles with a neutrino PFParticle
    if(neutrinoPFPID > 0) {
      auto& neutrinoPFP = tjs.pfps[neutrinoPFPID - 1];
      int vx3id = neutrinoPFP.Vx3ID[1];
      for(auto& pfp : tjs.pfps) {
        if(pfp.ID == 0 || pfp.ID == neutrinoPFPID) continue;
        if(pfp.TPCID != tpcid) continue;
        if(pfp.Vx3ID[0] != vx3id) continue;
        pfp.ParentID = (size_t)neutrinoPFPID;
        pfp.Primary = true;
        neutrinoPFP.DtrIDs.push_back(pfp.ID);
      } // pfp
    } // neutrino PFP exists    
  } // DefinePFPParents
  
  /////////////////////////////////////////
  void DefinePFPParentsTestBeam(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
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
    double fidZCut = tjs.ZLo + 2;
    for(auto& parPFP : tjs.pfps) {
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
      auto& vx3 = tjs.vtx3[parPFP.Vx3ID[1] - 1];
      // ensure that it is valid
      if(vx3.ID == 0) continue;
      // get a list of Tjs attached to this vertex. This will include the Tjs in the parent.
      auto vx3TjList = GetVtxTjIDs(tjs, vx3, score);
      if(vx3TjList.empty()) continue;
      // filter out the parent Tjs
      auto dtrTjlist = SetDifference(vx3TjList, parPFP.TjIDs);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" Dtrs:";
        for(auto dtjID : dtrTjlist) myprt<<" "<<dtjID<<"_"<<GetPFPIndex(tjs, dtjID);
      }
      // Add to the stack
      for(auto dtjID : dtrTjlist) {
        unsigned short pfpIndex = GetPFPIndex(tjs, dtjID);
        if(pfpIndex > tjs.pfps.size() - 1) continue;
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
      auto& dtr = tjs.pfps[lastPair.second - 1];
      auto& par = tjs.pfps[lastPair.first - 1];
      dtr.ParentID = par.ID;
      par.DtrIDs.push_back(dtr.ID);
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
      auto& vx3 = tjs.vtx3[dtr.Vx3ID[dtrEnd] - 1];
      float score = 0;
      auto vx3TjList = GetVtxTjIDs(tjs, vx3, score);
      if(vx3TjList.empty()) continue;
      // filter out the new parent
      auto dtrTjlist = SetDifference(vx3TjList, dtr.TjIDs);
      // put these onto the stack
      for(auto tjid : dtrTjlist) pardtr.push_back(std::make_pair(dtr.ID, tjid));
    } // nit
/*
    // deal with shower-like PFParticles, e.g. delta-rays
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != tpcid) continue;
      // look for delta-rays
      if(pfp.PDGCode != 11) continue;
      // ignore already assigned
      if(pfp.ParentID != pfp.ID) continue;
      // next look for a parent
      int pfpParentID = INT_MAX;
      unsigned short nParent = 0;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.ParentID == tj.ID) continue;
        unsigned short ppindex = GetPFPIndex(tjs, tj.ParentID);
        if(ppindex == USHRT_MAX) continue;
        int ppid = ppindex + 1;
        if(pfpParentID == INT_MAX) pfpParentID = ppid;
        if(ppid == pfpParentID) ++nParent;
      } // ii
      // look for a parent
      if(nParent > 1) {
        pfp.ParentID = (size_t)pfpParentID;
        auto& parpfp = tjs.pfps[pfpParentID - 1];
        parpfp.DtrIDs.push_back(pfp.ID);
      } // nParent > 1
    } // ipfp
*/
  } // DefinePFPParentsTestBeam

  ////////////////////////////////////////////////
  bool StorePFP(TjStuff& tjs, PFPStruct& pfp)
  {
    // stores the PFParticle in TJStuff
    if(pfp.ID < tjs.pfps.size()) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    if(!neutrinoPFP && pfp.TjIDs.empty()) return false;
    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)tjs.pfps.size() + 1) {
//      std::cout<<"StorePFP ID is wrong. Fixing it\n";
      pfp.ID = tjs.pfps.size() + 1;
    }
    // check the Tjs and set the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
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
        auto& tj = tjs.allTraj[tjids[ii] - 1];
        if(lastIpt[ii] < firstIpt[ii]) {
          if(tj.AlgMod[kSetDir]) {
//            std::cout<<"StorePFP "<<pfp.ID<<" Violating the SetDir flag for Tj "<<tj.ID<<"\n";
            tj.AlgMod[kSetDir] = false;
          }
          ReverseTraj(tjs, tj);
        } // lastIpt[ii] > firstIpt[ii]
      } // ii
    } // Tp3s exist    
    
    tjs.pfps.push_back(pfp);
    return true;
  } // StorePFP
  
  ////////////////////////////////////////////////
  bool InsideTPC(const TjStuff& tjs, Point3_t& pos, geo::TPCID& inTPCID)
  {
    // determine which TPC this point is in. This function returns false
    // if the point is not inside any TPC
    for (const geo::TPCID& tpcid: tjs.geom->IterateTPCIDs()) {
      const geo::TPCGeo& TPC = tjs.geom->TPC(tpcid);
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      TPC.LocalToWorld(local,world);
      unsigned int cstat = tpcid.Cryostat;
      unsigned int tpc = tpcid.TPC;
      // reduce the active area of the TPC by 1 cm to be consistent with FillWireHitRange
      if(pos[0] < world[0]-tjs.geom->DetHalfWidth(tpc,cstat) + 1) continue;
      if(pos[0] > world[0]+tjs.geom->DetHalfWidth(tpc,cstat) - 1) continue;
      if(pos[1] < world[1]-tjs.geom->DetHalfHeight(tpc,cstat) + 1) continue;
      if(pos[1] > world[1]+tjs.geom->DetHalfHeight(tpc,cstat) - 1) continue;
      if(pos[2] < world[2]-tjs.geom->DetLength(tpc,cstat)/2 + 1) continue;
      if(pos[2] > world[2]+tjs.geom->DetLength(tpc,cstat)/2 - 1) continue;
      inTPCID = tpcid;
      return true;
    } // tpcid
    return false;
  } // InsideTPC
  
  ////////////////////////////////////////////////
  void ReversePFP(TjStuff& tjs, PFPStruct& pfp)
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
    std::swap(pfp.StopFlag[0], pfp.StopFlag[1]);
    std::reverse(pfp.Tp3s.begin(), pfp.Tp3s.end());
    for(auto& tp3 : pfp.Tp3s) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) {
        tp3.Dir[xyz] *= -1;
        tp3.Dir[xyz] *= -1;
      } // xyz
    } // tp3
  } // ReversePFP
  
  ////////////////////////////////////////////////
  void FixDirection(TjStuff& tjs, PFPStruct& pfp)
  {
    // Reverse the direction vector of all Tp3s to be consistent with the start direction
    if(pfp.Tp3s.empty()) return;
    // ensure that the start direction is defined by finding a large component
    unsigned short useComp = 0;
    for(useComp = 0; useComp < 3; ++useComp) if(std::abs(pfp.Dir[0][useComp]) > 0.5) break;
    if(useComp == 3) {
      std::cout<<"FixDirection: pfp "<<pfp.ID<<" start direction not defined\n";
      return;
    }
    for(auto& tp3 : pfp.Tp3s) {
      FixDirection(tp3.Dir, pfp.Dir[0]);
      FixDirection(tp3.Dir, pfp.Dir[0]);
    } // tp3    
  } // FixDirection
  
  ////////////////////////////////////////////////
  void FixDirection(Vector3_t& ofDir, const Vector3_t& usingDir)
  {
    // Reverse the direction vector ofDir to be within pi/2 of the direction
    // vector usingDir.
    // Find a major component so we can compare the sign
    unsigned short useComp = 0;
    for(useComp = 0; useComp < 3; ++useComp) if(std::abs(usingDir[useComp]) > 0.5) break;
    if(useComp == 3) return;
    if(std::signbit(ofDir[useComp]) == std::signbit(usingDir[useComp])) return;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) ofDir[xyz] *= -1;
  } // FixDirection
  
  ////////////////////////////////////////////////
  void PrintTp3(std::string fcnLabel, const TjStuff& tjs, const TrajPoint3& tp3)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<fcnLabel<<" Pos";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(6)<<tp3.Pos[0]<<std::setw(6)<<tp3.Pos[1]<<std::setw(6)<<tp3.Pos[2];
    myprt<<" Pos";
    myprt<<std::setw(6)<<tp3.Pos[0]<<std::setw(6)<<tp3.Pos[1]<<std::setw(6)<<tp3.Pos[2];
    myprt<<std::fixed<<std::setprecision(3);
    myprt<<" Dir";
    myprt<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
    myprt<<" ChiDOF "<<std::setw(4)<<std::setprecision(1)<<tp3.ChiDOF;
    myprt<<" IsValid? "<<tp3.IsValid;
    myprt<<" nPtsFit "<<std::setw(4)<<tp3.nPtsFit;
    myprt<<" tj_ipt";
    for(auto tj2pt : tp3.Tj2Pts) {
      auto& tj = tjs.allTraj[tj2pt.id - 1];
      auto& tp = tj.Pts[tj2pt.ipt];
      myprt<<" "<<tj.ID<<"_"<<PrintPos(tjs, tp);
    } // tj2pt
  } // PrintTp3
  
  ////////////////////////////////////////////////
  void PrintTp3s(std::string someText, const TjStuff& tjs, const PFPStruct& pfp, short printPts)
  {
    if(pfp.Tp3s.empty()) return;
    mf::LogVerbatim myprt("TC");
    if(printPts < 0) {
      // print the head if we are print all points
      myprt<<someText<<" pfp "<<pfp.ID<<" DirectionFixed? "<<pfp.DirectionFixed<<"\n";
      myprt<<someText<<"  ipt ________Pos________ Path   ________Dir________ ChiDOF Val? nFit dang  Kink  Tj_ipt \n";
    }
    // print the start
    myprt<<someText<<"    ";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(7)<<pfp.XYZ[0][0]<<std::setw(7)<<pfp.XYZ[0][1]<<std::setw(7)<<pfp.XYZ[0][2];
    myprt<<"            ";
    myprt<<std::fixed<<std::setprecision(3);
    myprt<<std::setw(7)<<pfp.Dir[0][0]<<std::setw(7)<<pfp.Dir[0][1]<<std::setw(7)<<pfp.Dir[0][2];
    myprt<<" <--- pfp.XYZ[0] \n";
    
    unsigned short fromPt = 0;
    unsigned short toPt = pfp.Tp3s.size() - 1;
    if(printPts >= 0) fromPt = toPt;
    for(unsigned short ipt = fromPt; ipt <= toPt; ++ipt) {
      auto tp3 = pfp.Tp3s[ipt];
      myprt<<someText<<std::setw(4)<<ipt;
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<std::setw(7)<<tp3.Pos[0]<<std::setw(7)<<tp3.Pos[1]<<std::setw(7)<<tp3.Pos[2];
      myprt<<std::setprecision(1)<<std::setw(5)<<PosSep(tp3.Pos, pfp.XYZ[0]);
      myprt<<std::setprecision(3)<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
      myprt<<std::setprecision(1)<<std::setw(5)<<tp3.ChiDOF;
      myprt<<std::setw(5)<<tp3.IsValid;
      myprt<<std::setw(5)<<tp3.nPtsFit;
      myprt<<std::setprecision(3)<<std::setw(7)<<DeltaAngle(pfp.Dir[0], tp3.Dir);
      // Calculate the kink angle at point ipt, using the two points that are
      // +/- 1 cm on either side of that point
      double sep = 1;
      myprt<<std::setprecision(3)<<std::setw(7)<<KinkAngle(tjs, pfp.Tp3s, ipt, sep);
      for(auto tj2pt : tp3.Tj2Pts) {
        auto& tj = tjs.allTraj[tj2pt.id - 1];
        auto& tp = tj.Pts[tj2pt.ipt];
        myprt<<" "<<tj.ID<<"_"<<PrintPos(tjs, tp);
      } // tj2pt
      myprt<<"\n";
    } // ipt
    // print the end
    myprt<<someText<<"    ";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(7)<<pfp.XYZ[1][0]<<std::setw(7)<<pfp.XYZ[1][1]<<std::setw(7)<<pfp.XYZ[1][2];
    myprt<<"            ";
    myprt<<std::fixed<<std::setprecision(3);
    myprt<<std::setw(7)<<pfp.Dir[1][0]<<std::setw(7)<<pfp.Dir[1][1]<<std::setw(7)<<pfp.Dir[1][2];
    myprt<<" <--- pfp.XYZ[1] \n";
  } // PrintTp3s

} // namespace
