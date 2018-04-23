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
    unsigned short newtjid = newTj;
    auto& ntj = tjs.allTraj[newTj - 1];
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;
    // put the X positions of the new Tj into a vector for matching
    std::vector<float> xpos(ntj.Pts.size());
    geo::PlaneID planeID = DecodeCTP(ntj.CTP);
    for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
      auto& ntp = ntj.Pts[npt];
      if(ntp.Chg <= 0) continue;
      xpos[npt] = tjs.detprop->ConvertTicksToX(ntp.Pos[1]/tjs.UnitsPerTick, planeID);
    } // npt
    
    if(!tjs.mallTraj.empty()) {
      for(unsigned int ipt = 0; ipt < tjs.mallTraj.size(); ++ipt) {
        auto& tj2pt = tjs.mallTraj[ipt];
        if(tj2pt.id > tjs.allTraj.size()) continue;
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
  void UpdateTp3s(TjStuff& tjs, PFPStruct& pfp, int oldTj, int newTj)
  {
    // Replaces occurrences of oldTj with newTj in the pfp vector of Tp3s
    if(oldTj <= 0 || oldTj > (int)tjs.allTraj.size()) return;
    if(newTj <= 0 || newTj > (int)tjs.allTraj.size()) return;
    if(tjs.mallTraj.empty() && pfp.Tp3s.empty()) return;
    
    // convert from int to unsigned short
    unsigned short oldtjid = oldTj;
    unsigned short newtjid = newTj;
    auto& ntj = tjs.allTraj[newTj - 1];
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;
    // put the X positions of the new Tj into a vector for matching
    std::vector<float> xpos(ntj.Pts.size());
    geo::PlaneID planeID = DecodeCTP(ntj.CTP);
    for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
      auto& ntp = ntj.Pts[npt];
      if(ntp.Chg <= 0) continue;
      xpos[npt] = tjs.detprop->ConvertTicksToX(ntp.Pos[1]/tjs.UnitsPerTick, planeID);
    } // npt

    for(auto& tp3 : pfp.Tp3s) {
      // check each of the Tj2Pts associated with this space point
      for(auto& tj2pt : tp3.Tj2Pts) {
        if(tj2pt.id > tjs.allTraj.size()) continue;
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
        tjs.mallTraj[icnt].inShower = tj.AlgMod[kInShower];
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
  void AttachVertices(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // try to attach a vertex on both ends of the pfp using the positions XYZ[end]. This function
    // doesn't need to have Tp3s associated with the pfp but probably doesn't work as well
    if(pfp.ID == 0) return;
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<unsigned short, 2> imbest {{ 0, 0 }};
    // Ignore any separation larger than (10 cm)^2
    std::array<float, 2> best {{ 100.0f, 100.0f }};
    for(unsigned short end = 0; end < 2; ++end) pfp.Vx3ID[end] = 0;
    auto vx3list = GetPFPVertices(tjs, pfp);
    for(unsigned short end = 0; end < 2; ++end) {
      for(auto vx3id : vx3list) {
        auto& vx3 = tjs.vtx3[vx3id - 1];
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        Point3_t vxpos = {{ vx3.X, vx3.Y, vx3.Z}};
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
        if(prt) mf::LogVerbatim("TC")<<"AV0: attach 3V"<<imbest[0]<<" to P"<<pfp.ID<<" to end 0";
      } else {
        pfp.Vx3ID[1] = imbest[0];
        if(prt) mf::LogVerbatim("TC")<<"AV1: attach 3V"<<imbest[0]<<" to P"<<pfp.ID<<" to end 1";
      }
      return;
    } else {
      // have different vertices on each end
      for(unsigned short end = 0; end < 2; ++end) {
        if(imbest[end] == 0) continue;
        pfp.Vx3ID[end] = imbest[end];
        if(prt) mf::LogVerbatim("TC")<<"AV: attach 3V"<<imbest[end]<<" to P"<<pfp.ID<<" end "<<end;
      } // end
    } // have different vertices (or no vertices) on each end
  } // AttachVertices

  /////////////////////////////////////////
  bool SetNewStart(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Analyzes the space point collection and the Tjs in the pfp to find a new start
    // position. The Tp3s found in FindCompleteness are ordered by increasing X. The general direction 
    // pfp.Dir[0] and the average position of all points in Tp3s was stored in pfp.XYZ[0]. This function
    // rotates each tp3 into this coordinate system to determine (along, trans) for each point. The min (max)
    // value of along defines the start (end) of the trajectory. The along coordinate is stashed in the
    // pfp dEdxErr variable for use later on.
    if(pfp.ID == 0 || pfp.TjIDs.empty()) return false;
    if(pfp.Tp3s.size() < 2) return false;

    // find the projection along the general direction relative to the average position
    float minAlong = 1E6;
    unsigned short minPt = 0;
    float maxAlong = -1E6;
    unsigned short maxPt = 0;
    std::vector<SortEntry> sortVec(pfp.Tp3s.size());
    for(unsigned short ipt = 0; ipt < pfp.Tp3s.size(); ++ipt) {
      auto& tp3 = pfp.Tp3s[ipt];
      auto ptDir = PointDirection(pfp.XYZ[0], tp3.Pos);
      double along = DotProd(pfp.Dir[0], ptDir) * PosSep(pfp.XYZ[0], tp3.Pos);
      // stash along in dEdxErr for later use
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = along;
      tp3.dEdxErr = along;
      // find the min (max) 
      if(tp3.dEdxErr < minAlong) {
        minAlong = tp3.dEdxErr;
        minPt = ipt;
      }
      if(tp3.dEdxErr > maxAlong) {
        maxAlong = tp3.dEdxErr;
        maxPt = ipt;
      }
    } // tp3
    
    pfp.XYZ[0] = pfp.Tp3s[minPt].Pos;
    pfp.XYZ[1] = pfp.Tp3s[maxPt].Pos;

    if(prt) {
      mf::LogVerbatim("TC")<<"SNS: P"<<pfp.ID<<" minPt "<<minPt<<" maxPt "<<maxPt;
      PrintTp3("minPt", tjs, pfp.Tp3s[minPt]);
      PrintTp3("maxPt", tjs, pfp.Tp3s[maxPt]);
    }
    
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put them into order
    std::vector<TrajPoint3> temp;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) temp.push_back(pfp.Tp3s[sortVec[ii].index]);
    pfp.Tp3s = temp;

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
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        Point3_t vx3pos = {{vx3.X, vx3.Y, vx3.Z}};
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
  void FollowTp3s(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Step through the set of Tp3s on this pfp to create a trajectory. The start and end points
    // are assumed to be Tp3s[0] and Tp3s[Tp3s.size()-1] respectively. 
    
    if(pfp.Tp3s.size() < 2) return;
    
    unsigned short startPt = 0;
    unsigned short endPt = pfp.Tp3s.size() - 1;
    // divide the trajectory in 5 cm long sections. The average position of the Tp3s in
    // each section will be found. Tp3s 
    constexpr float sectionLen = 5;
    float endPathLen = sectionLen;
    std::vector<Vector3_t> sectionPos;
    sectionPos.push_back(pfp.Tp3s[0].Pos);
    std::vector<unsigned short> sectionPt;
    sectionPt.push_back(0);
    for(unsigned short section = 0; section < 100; ++section) {
      // a point to find the average position in this section
      Point3_t avePos {{0,0,0}};
      unsigned short cnt = 0;
      for(unsigned short ipt = startPt; ipt < endPt; ++ipt) {
        auto& tp3 = pfp.Tp3s[ipt];
        // The path length along the direction vector from the start point to the end
        // point was stashed in dEdxErr 
        if(tp3.dEdxErr < endPathLen) {
          // still in the same section - sum and continue
          for(unsigned short xyz = 0; xyz < 3; ++xyz) avePos[xyz] += tp3.Pos[xyz];
          ++cnt;
          continue;
        }
        // entered the next section. Check for a failure
        if(cnt == 0) continue;
        // calculate the average position
        for(unsigned short xyz = 0; xyz < 3; ++xyz) avePos[xyz] /= cnt;
//        std::cout<<"Section "<<section<<" cnt "<<cnt<<" pos"<<std::fixed<<std::setprecision(1);
//        std::cout<<" "<<avePos[0]<<" "<<avePos[1]<<" "<<avePos[2]<<"\n";
        sectionPos.push_back(avePos);
        sectionPt.push_back(ipt);
        startPt = ipt;
        endPathLen += sectionLen;
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
      auto& tj = tjs.allTraj[tjid - 1];
      for(auto& tp : tj.Pts) tp.Environment[kEnvFlag] = false;
    } // tjid
    // set the bits true for the first point
    for(auto tj2pt : pfp.Tp3s[0].Tj2Pts) {
      tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt].Environment[kEnvFlag] = true;
    }
    // create a vector of new Tp3s that will replace pfp.Tp3s
    std::vector<TrajPoint3> ntp3;
    ntp3.push_back(pfp.Tp3s[0]);
    // 2D position (WSE units) of the TPs at the start of this section. We will require that all 2D TPs are
    // less than sectionLen (in WSE units) from this point.
    std::vector<Point2_t> startPos2D(tjs.NumPlanes);
    // collect Tp3s in each section
    unsigned short lastPtAdded = 0;
    for(unsigned short section = 1; section < sectionPt.size(); ++section) {
      Point3_t startPos = sectionPos[section - 1];
      Point3_t endPos = sectionPos[section];
      auto dir = PointDirection(startPos, endPos);
      // define the pfp start direction
      if(section == 1) pfp.Dir[0] = dir;
      // and the end direction
      pfp.Dir[1] = dir;
      // define the 2D positions for this point in each plane
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        startPos2D[plane][0] = tjs.geom->WireCoordinate(sectionPos[section - 1][1], sectionPos[section - 1][2], planeID);
        startPos2D[plane][1] = tjs.detprop->ConvertXToTicks(sectionPos[section - 1][0], planeID) * tjs.UnitsPerTick;
      } // plane
      for(unsigned short ipt = sectionPt[section - 1]; ipt < sectionPt[section]; ++ipt) {
        auto& tp3 = pfp.Tp3s[ipt];
        // count the number of Tps in this Tp3 that are already used in the trajectory
        unsigned short nused = 0;
        bool big2DSep = false;
        for(auto tj2pt : pfp.Tp3s[ipt].Tj2Pts) {
          auto& tp = tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt];
          if(tp.Environment[kEnvFlag]) ++nused;
          unsigned short plane = DecodeCTP(tp.CTP).Plane;
          float sep2D = PosSep(startPos2D[plane], tp.Pos) * tjs.WirePitch;
          if(sep2D > sectionLen) big2DSep = true;
        } // tj2pt
        if(big2DSep || nused > 1) continue;
        double sep = PosSep(startPos, tp3.Pos);
        auto ptDir = PointDirection(startPos, tp3.Pos);
        double costh = DotProd(dir, ptDir);
        if(costh < 0) continue;
        double sinth = sqrt(1 - costh * costh);
        double trans = sinth * sep;
        // cut on the transverse separation (cm)
        if(trans > 0.5) continue;
/*
        std::cout<<section<<" ipt "<<ipt<<" trans "<<trans<<" tj_ipt";
        for(auto tj2pt : tp3.Tj2Pts) std::cout<<" "<<tj2pt.id - 1<<"_"<<tj2pt.ipt;
        std::cout<<"\n";
*/
        tp3.Trans = trans;
        tp3.Dir = dir;
        ntp3.push_back(tp3);
        // set the flag
        for(auto tj2pt : tp3.Tj2Pts) tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt].Environment[kEnvFlag] = true;
        lastPtAdded = ipt;
      } // ipt
    } // section
    
    if(lastPtAdded != endPt) ntp3.push_back(pfp.Tp3s[endPt]);
    
    if(prt) {
      float len = PosSep(ntp3[0].Pos, ntp3[ntp3.size()-1].Pos);
      mf::LogVerbatim("TC")<<"FollowTp3s: Tp3s size in "<<pfp.Tp3s.size()<<" size out "<<ntp3.size()<<" len "<<std::fixed<<std::setprecision(2)<<len;
    }
    
    pfp.Tp3s = ntp3;
    // The directions wer set above. Set the start and end positions. Note that the start position
    // may have been previously determined by a vertex but that is now superseded by the actual start
    // of the pfp
    pfp.XYZ[0] = ntp3[0].Pos;
    pfp.XYZ[1] = ntp3[ntp3.size()-1].Pos;
//    if(prt) PrintTp3s("FTp3o", tjs, pfp, -1);

  } // FollowTp3s
  /////////////////////////////////////////
  bool FitTp3s(TjStuff& tjs, const std::vector<TrajPoint3>& tp3s, Point3_t& pos, Vector3_t& dir, float& rCorr)
  {
    return FitTp3s(tjs, tp3s, 0, tp3s.size(), pos, dir, rCorr);
  } // FitTp3s
  
  /////////////////////////////////////////
  bool FitTp3s(TjStuff& tjs, const std::vector<TrajPoint3>& tp3s, unsigned short fromPt, unsigned short toPt, Point3_t& pos, Vector3_t& dir, float& rCorr)
  {
    // Fits the Tj2Pts points in Tp3s to a line
    if(tp3s.size() < 3) return false;
    if(fromPt >= toPt) return false;
    if(toPt > tp3s.size()) return false;

    // temp vectors to ensure that a TP is only used once
    std::vector<unsigned short> useID;
    std::vector<unsigned short> useIpt;
    std::vector<unsigned short> cntInPln(tjs.NumPlanes);
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
        auto& tj = tjs.allTraj[tj2pt.id - 1];
        ++cntInPln[DecodeCTP(tj.CTP).Plane];
      } // tj2pt
    } // ipt
    // ensure there are at least two points in at least two planes
    unsigned short enufInPlane = 0;
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) if(cntInPln[plane] > 1) ++enufInPlane;
    if(enufInPlane < 2) return false;
    
    const unsigned int nvars = 4;
    unsigned int npts = useID.size();
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    
    // X origin
    double x0 = 0;
    for(unsigned short ipt = 0; ipt < useID.size(); ++ipt) {
      auto& tp = tjs.allTraj[useID[ipt] - 1].Pts[useIpt[ipt]];
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      x0 += tjs.detprop->ConvertTicksToX(tp.Pos[1]/tjs.UnitsPerTick, planeID);
    }
    x0 /= (double)useID.size();

//<<<<<<< HEAD
    double wght = 1;
    for(unsigned short ipt = 0; ipt < useID.size(); ++ipt) {
      auto& tp = tjs.allTraj[useID[ipt] - 1].Pts[useIpt[ipt]];
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
    std::cout<<"FTP3s: "<<useID.size()<<" cntInPln "<<cntInPln[0]<<" "<<cntInPln[1]<<" "<<cntInPln[2]<<"\n";
/*
    =======
    // Set the direction vectors of all points to be consistent with the general
    // start direction
    FixDirection(tjs, pfp);
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    Vector3_t startDir = {{0.0, 0.0, 0.0}};
    // Find the average direction using the points in the first 10 cm
    for(auto& tp3 : pfp.Tp3s) {
      if(PosSep2(tp3.Pos, startPos) > 100) break;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) startDir[xyz] += tp3.Dir[xyz];
    } //  tp3
    SetMag(startDir, 1);
    pfp.Dir[0] = startDir;
    if(prt) mf::LogVerbatim("TC")<<"SBDFS: start direction "<<std::fixed<<std::setprecision(2)<<pfp.Dir[0][0]<<" "<<pfp.Dir[0][1]<<" "<<pfp.Dir[0][2];
    // do the same at the other end
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    Vector3_t endDir = {{0.0, 0.0, 0.0}};
    for(unsigned short ii = 0; ii < pfp.Tp3s.size(); ++ii) {
      auto& tp3 = pfp.Tp3s[pfp.Tp3s.size() - 1 - ii];
      if(PosSep2(tp3.Pos, endPos) > 100) break;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) endDir[xyz] += tp3.Dir[xyz];
    } // ii
    SetMag(endDir, 1);
    pfp.Dir[1] = endDir;
    if(prt) mf::LogVerbatim("TC")<<"SBDFS: end direction "<<std::fixed<<std::setprecision(2)<<pfp.Dir[1][0]<<" "<<pfp.Dir[1][1]<<" "<<pfp.Dir[1][2];
>>>>>>> develop
*/    
    return true;

  } // FitTp3s
  
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
  void FindCompleteness(TjStuff& tjs, PFPStruct& pfp, bool doFit, bool fillTp3s, bool prt)
  {
    // Calculate the 3D-matching completeness of the set of Tjs in pfp.TjIDs and store in pfp.EffPur.
    // The completeness for each Tj is put in pfp.TjCompleteness. The TP-weighted average completeness
    // is put in pfp.EffPur. This function also fits the matching points to a 3D line and puts the
    // position and direction in pfp.XYZ[0] and pfp.Dir[0]. The absolute value of the linear correlation
    // coefficients are averaged and stored in pfp.AspecRatio. The pfp.TP3s vector is optionally filled.
    
    if(pfp.TjIDs.size() < 2) return;
    if(tjs.Match3DCuts[0] <= 0) return;
    // This function uses mallTraj but it isn't necessarily a failure if it doesn't exist
    if(tjs.mallTraj.size() < 6) return;
    if(tjs.NumPlanes < 3) return;

    pfp.TjCompleteness.resize(pfp.TjIDs.size());
    std::fill(pfp.TjCompleteness.begin(), pfp.TjCompleteness.end(), 0);
    if(fillTp3s) pfp.Tp3s.clear();

    bool twoPlanes = (tjs.NumPlanes == 2);
    bool twoTjs = (pfp.TjIDs.size() == 2);
    double yzcut = 1.5 * tjs.Match3DCuts[0];
    
    // initialize the fit sums
    Point3_t point;
    if(doFit) Fit3D(0, point, pfp.XYZ[0], pfp.Dir[0]);
    
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
    unsigned int maxTp3Size = 0;
    // Initialize the vectors
    for(unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
      if(pfp.TjIDs[itj] <= 0) {
        std::cout<<"FindCompleteness: Bad tjid "<<pfp.TjIDs[itj]<<"\n";
        return;
      }
      tjids[itj] = pfp.TjIDs[itj];
      auto& tj = tjs.allTraj[pfp.TjIDs[itj] - 1];
      // allow each point to be matched no more than N times
      maxTp3Size += 10 * (tj.EndPt[1] - tj.EndPt[0] + 1);
      std::vector<bool> tmp(tj.Pts.size(), false);
      tjptMat2.push_back(tmp);
      if(!twoPlanes) tjptMat3.push_back(tmp);
      tjplane.push_back(DecodeCTP(tj.CTP).Plane);
    } // tjid
    if(maxTp3Size < 100) maxTp3Size = 100;
    // count of triple matches on dead wires in the 3rd plane
//    std::vector<unsigned short> deadCnt(pfp.TjIDs.size());
    
    for(unsigned int ipt = 0; ipt < tjs.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = tjs.mallTraj[ipt];
      unsigned short indx = 0;
      for(indx = 0; indx < tjids.size(); ++indx) if(iTjPt.id == tjids[indx]) break;
      // require that the Tj ID of this point be in the list
      if(indx == tjids.size()) continue;
      auto& itj = tjs.allTraj[iTjPt.id - 1];
//      if(itj.AlgMod[kMat3D]) continue;
      auto& itp = itj.Pts[iTjPt.ipt];
      unsigned short iplane = DecodeCTP(itp.CTP).Plane;
      unsigned short tpc = DecodeCTP(itp.CTP).TPC;
      unsigned short cstat = DecodeCTP(itp.CTP).Cryostat;
      for(unsigned int jpt = ipt + 1; jpt < tjs.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = tjs.mallTraj[jpt];
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
        auto& jtj = tjs.allTraj[jTjPt.id - 1];
//        if(jtj.AlgMod[kMat3D]) continue;
        auto& jtp = jtj.Pts[jTjPt.ipt];
        TrajPoint3 ijtp3;
        if(!MakeTp3(tjs, itp, jtp, ijtp3, true)) continue;
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
/*
          unsigned short kndx = 0;
          for(kndx = 0; kndx < tjplane.size(); ++kndx) if(kplane == tjplane[kndx]) break;
          if(kndx == deadCnt.size()) continue;
          ++deadCnt[kndx];
*/
          // accumulate the fit sums?
          if(doFit) Fit3D(1, ijtp3.Pos, pfp.XYZ[0], pfp.Dir[0]);
          // fill Tp3s?
          if(fillTp3s && pfp.Tp3s.size() < maxTp3Size) pfp.Tp3s.push_back(ijtp3);
          continue;
        } // dead wire in kplane
        for(unsigned int kpt = jpt + 1; kpt < tjs.mallTraj.size(); ++kpt) {
          auto& kTjPt = tjs.mallTraj[kpt];
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
          auto& ktj = tjs.allTraj[kTjPt.id - 1];
//          if(ktj.AlgMod[kMat3D]) continue;
          auto& ktp = ktj.Pts[kTjPt.ipt];
          TrajPoint3 iktp3;
          if(!MakeTp3(tjs, itp, ktp, iktp3, true)) continue;
          if(std::abs(ijtp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
          if(std::abs(ijtp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
          // make a copy of ijtp3 -> ijktp3
          auto ijktp3 = ijtp3;
          // add the Tj2Pt to it
          ijktp3.Tj2Pts.push_back(kTjPt);
          // accumulate the fit sums
          if(doFit) Fit3D(1, iktp3.Pos, pfp.XYZ[0], pfp.Dir[0]);
          // fill Tp3s?
          if(fillTp3s && pfp.Tp3s.size() < maxTp3Size) {
            // update the charge
            ijktp3.dEdx = (2 * ijktp3.dEdx + ktp.Chg) / 3;
            pfp.Tp3s.push_back(ijktp3);
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
    Fit3D(2, point, pfp.XYZ[0], pfp.Dir[0]);
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
        myprt<<" MCSMom "<<tj.MCSMom<<" InShower? "<<tj.AlgMod[kInShower];
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
  void FindMissedTjsInTp3s(TjStuff& tjs, PFPStruct& pfp, std::vector<int>& missTjs, std::vector<float>& missFrac)
  {
    // compare the Tjs in pfp.TjIDs with the Tjs in Tp3s and return a list of Tjs
    // in Tp3s that aren't in pfp.TjIDs
    missTjs.clear();
    missFrac.clear();
    if(pfp.TjIDs.empty() || pfp.Tp3s.empty()) return;
    
    // determine the projection of the pfp direction vector in each plane.
    // Don't try to merge if the projection is small
    std::vector<float> projInPlane(tjs.NumPlanes);
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      auto tp = MakeBareTP(tjs, pfp.XYZ[0], pfp.Dir[0], inCTP);
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
          auto& mtj = tjs.allTraj[tj2pt.id - 1];
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
      auto& mtj = tjs.allTraj[mtjid - 1];
      // ignore if there is a high-score vertex between the missed tj and those in the pfp list
      if(SharesHighScoreVx(tjs, pfp, mtj)) continue;
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
        auto& tj = tjs.allTraj[tjid - 1];
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
  bool SharesHighScoreVx(TjStuff& tjs, const PFPStruct& pfp, const Trajectory& tj)
  {
    // returns true if tj with tjID shares a high-score 3D vertex with any
    // tj in pfp.TjIDs
    for(unsigned short end = 0; end < 2; ++end) {
      if(tj.VtxID[end] == 0) continue;
      auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
      if(!vx2.Stat[kHiVx3Score]) continue;
      std::vector<int> vtjlist = GetVtxTjIDs(tjs, vx2);
      auto shared = SetIntersection(vtjlist, pfp.TjIDs);
      if(!shared.empty()) return true;
    } // end
    return false;
  } // SharesHighScoreVx
  
  /////////////////////////////////////////
  void Fit3D(unsigned short mode, Point3_t point, Point3_t& fitPos, Vector3_t& fitDir)
  {
    // initialize, accumulate and fit the points
    
    // 3D fit sum variables
    static double fSum, fSumx, fSumy, fSumz, fSumx2, fSumy2, fSumz2, fSumxz, fSumyz;

    if(mode == 0) {
      fSum = 0; fSumx = 0; fSumy = 0; fSumz = 0; fSumx2 = 0; fSumy2 = 0; fSumz2 = 0; fSumxz = 0; fSumyz = 0;
      return;
    }
    // accumulate
    if(mode == 1) {
      fSum += 1;
      fSumx += point[0];
      fSumy += point[1];
      fSumz += point[2];
      fSumx2 += point[0] * point[0];
      fSumy2 += point[1] * point[1];
      fSumz2 += point[2] * point[2];
      fSumxz += point[0] * point[2];
      fSumyz += point[1] * point[2];
      return;
    }
    
    if(fSum < 2) return;
    // just use the average for the position
    fitPos[0] = fSumx / fSum;
    fitPos[1] = fSumy / fSum;
    fitPos[2] = fSumz / fSum;
    // calculate the direction
    double delta = fSum * fSumz2 - fSumz * fSumz;
    if(delta == 0) return;
    double Bx = (fSumxz * fSum - fSumz * fSumx) / delta;
    double By = (fSumyz * fSum - fSumz * fSumy) / delta;
    double norm = sqrt(1 + Bx * Bx + By * By);
    fitDir[0] = Bx / norm;
    fitDir[1] = By / norm;
    fitDir[2] = 1 / norm;
/* 
    double den = delta*(fSum*fSumx2 - fSumx*fSumx);
    double rx = 1;
    if(den > 0) rx = std::abs(fSum*fSumxz - fSumx*fSumz) / sqrt(den);
    den = delta*(fSum*fSumy2 - fSumy*fSumy);
    double ry = 1;
    if(den > 0) ry = std::abs(fSum*fSumyz - fSumy*fSumz) / sqrt(den);
    aspectRatio = 0.5 * (rx + ry);
*/    
  } // Fit3D

  /////////////////////////////////////////
  bool CheckAndMerge(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Check the completeness of this pfp and try to add Tjs to improve it
    
    if(!tjs.UseAlg[kMat3DMerge]) return true;
    if(pfp.TjIDs.size() < 2) return false;
    if(tjs.Match3DCuts[0] <= 0) return false;
    std::vector<int> missTjs;
    std::vector<float> missFrac;
    bool tryThis = true;
    
    for(unsigned short nit = 0; nit < 4; ++nit) {
      unsigned short oldSize = pfp.TjIDs.size();
      // re-find the completeness, do the fit, fill Tp3s so that
      // FindMissedTjsInTp3s can check for matched points in Tjs that
      // aren't in pfp.TjIDs
      FindCompleteness(tjs, pfp, true, true, prt);
      FindMissedTjsInTp3s(tjs, pfp, missTjs, missFrac);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"CAM: nit "<<nit<<" Tjs";
        for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
        myprt<<" Completeness "<<std::fixed<<std::setprecision(2)<<pfp.EffPur;
        myprt<<" Tp3s size "<<pfp.Tp3s.size();
        if(!missTjs.empty()) {
          myprt<<"\n missTj_Frac :";
          for(unsigned short ii = 0; ii < missTjs.size(); ++ii) {
            myprt<<" T"<<missTjs[ii]<<"_"<<std::fixed<<std::setprecision(2)<<missFrac[ii];
          }
        } // missTjs exist
      } // prt
      if(pfp.Tp3s.empty()) return false;
/* This doesn't seem to do anything
      if(!missTjs.empty() && pfp.MatchVecIndex < tjs.matchVec.size()) {
        // look for the set intersection of the missed Tjs with the matchVec list
        for(unsigned short ims = 0; ims < pfp.MatchVecIndex + 10; ++ims) {
          if(ims >= tjs.matchVec.size()) break;
          auto& ms = tjs.matchVec[ims];
          if(ms.Count == 0) continue;
          std::vector<int> shared = SetIntersection(ms.TjIDs, missTjs);
          if(shared.size() < 2) continue;
          // check the max length Tj and cut on the minimum aspect ratio
          float mtjl = MaxTjLen(tjs, ms.TjIDs);
          float mcsmom = MCSMom(tjs, ms.TjIDs);
          if(prt) mf::LogVerbatim("TC")<<" chk ims "<<ims<<" mtjl "<<mtjl<<" MCSMom "<<mcsmom;
          bool tryMerge = (MaxTjLen(tjs, ms.TjIDs) > 10 && MCSMom(tjs, ms.TjIDs) > tjs.Match3DCuts[3]);
          if(!tryMerge) continue;
          for(auto tjid : ms.TjIDs) {
            if(std::find(shared.begin(), shared.end(), tjid) != shared.end()) continue;
            auto& tj = tjs.allTraj[tjid - 1];
            if(tj.AlgMod[kKilled]) continue;
            if(tj.AlgMod[kMat3D]) continue;
            // check for PDGCode compatibility - muons and delta rays
            if(pfp.PDGCode == 13 && tj.PDGCode == 11) continue;
            if(pfp.PDGCode == 11 && tj.PDGCode == 13) continue;
            float dotProd = DotProd(ms.Dir, pfp.Dir[0]);
            if(prt) mf::LogVerbatim("TC")<<" add T"<<tjid<<" DotProd "<<std::setprecision(3)<<dotProd;
            if(dotProd < tjs.Match3DCuts[6]) continue;
            // make a trial pfp with this tj added
            auto trial = pfp;
            trial.Tp3s.clear();
            trial.TjIDs.push_back(tjid);
            // find the completeness, do the fit, don't fill Tp3s
            FindCompleteness(tjs, trial, true, false, true);
            if(prt) mf::LogVerbatim("TC")<<" do something with Trial "<<trial.EffPur;
            std::cout<<" do something with Trial "<<trial.EffPur<<"\n";
            break;
          } // tjid
        } // ims
      } // try to add missed tjs
*/
      // Check TjCompleteness if nothing was added
      if(tryThis && pfp.TjIDs.size() == oldSize) {
        // look for the situation where the Tj with the worst completeness is also
        // the longest one - a bit
        float maxTjLen = 0.9 * MaxTjLen(tjs, pfp.TjIDs);
        int tjWithWorstCompleteness = 0;
        float worstCompleteness = 0.5;
        for(unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
          if(pfp.TjCompleteness[itj] > worstCompleteness) continue;
          auto& tj = tjs.allTraj[pfp.TjIDs[itj] - 1];
          float tjLen = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
          if(tjLen < maxTjLen) continue;
          worstCompleteness = pfp.TjCompleteness[itj];
          tjWithWorstCompleteness = pfp.TjIDs[itj];
        } // itj
        if(tjWithWorstCompleteness > 0 && pfp.MatchVecIndex < tjs.matchVec.size()) {
          for(unsigned short ims = 0; ims < pfp.MatchVecIndex + 10; ++ims) {
            if(ims >= tjs.matchVec.size()) break;
            // ignore the self matchVecIndex
            if(ims == pfp.MatchVecIndex) continue; 
            auto& ms = tjs.matchVec[ims];
            if(ms.Count == 0) continue;
            if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tjWithWorstCompleteness) == ms.TjIDs.end()) continue;
            // check the max length Tj and cut on the minimum aspect ratio
            bool tryMerge = (MaxTjLen(tjs, ms.TjIDs) > tjs.Match3DCuts[5] && MCSMom(tjs, ms.TjIDs) > tjs.Match3DCuts[3]);
            if(!tryMerge) continue;
            float dotProd = DotProd(ms.Dir, pfp.Dir[0]);
            if(prt) mf::LogVerbatim("TC")<<" ms "<<ims<<" looks interesting. dotProd "<<std::setprecision(3)<<dotProd;
            if(dotProd < tjs.Match3DCuts[6]) continue;
            for(auto& tjid : ms.TjIDs) {
              if(tjid == tjWithWorstCompleteness) continue;
              auto& tj = tjs.allTraj[tjid - 1];
              if(tj.AlgMod[kKilled]) continue;
              if(tj.AlgMod[kMat3D]) continue;
              // check for PDGCode compatibility - muons and delta rays
              if(pfp.PDGCode == 13 && tj.PDGCode == 11) continue;
              if(pfp.PDGCode == 11 && tj.PDGCode == 13) continue;
              if(prt) mf::LogVerbatim("TC")<<" add T"<<tjid;
              pfp.TjIDs.push_back(tjid);
              PFPVxTjOK(tjs, pfp, prt);
            } // tjid
          } // ims
          tryThis = false;
        } // Tj has poor completeness
      } // nothing was added
      if(pfp.TjIDs.size() == oldSize) break;
    } // nit
    // At this point

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"CAM done P"<<pfp.ID<<" Tjs";
      for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
      myprt<<" Completeness "<<std::fixed<<std::setprecision(2)<<pfp.EffPur;
      myprt<<" Tp3s size "<<pfp.Tp3s.size();
      myprt<<" MCSMom "<<MCSMom(tjs, pfp.TjIDs);
    } // prt

    // something bad happened.
    if(pfp.TjIDs.size() < 2) return false;
    return true;

  } // CheckAndMerge

  /////////////////////////////////////////
  unsigned short WiresSkippedInCTP(TjStuff& tjs, std::vector<int>& tjids, CTP_t inCTP)
  {
    // counts the number of wires between the end points of all Tjs in the list of tjids
    // in inCTP where there is no TP with charge
    if(tjids.empty()) return 0;
    
    // find the min and max Pos[0] positions
    float fLoWire = 1E6;
    float fHiWire = -1E6;
    for(auto tjid : tjids) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.CTP != inCTP) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        float endWire = tj.Pts[tj.EndPt[end]].Pos[0];
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
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.CTP != inCTP) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
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
  float LengthInCTP(TjStuff& tjs, std::vector<int>& tjids, CTP_t inCTP)
  {
    // Calculates the maximum length between the end points of Tjs in the list of tjids in inCTP
    if(tjids.empty()) return 0;
    // put the end point positions into a vector
    std::vector<Point2_t> endPos;
    for(auto tjid : tjids) {
      auto& tj = tjs.allTraj[tjid - 1];
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
  bool AddMissedTj(TjStuff& tjs, PFPStruct& pfp, unsigned short itj, bool looseCuts, bool prt)
  {
    // The Tj pfp.TjIDs[itj] has poor completeness. Search tjs.matchVec for
    // the occurrence of this tj with a large completeness AND the occurrence 
    // of another tj in pfp.TjIDs.
    if(itj > pfp.TjIDs.size() - 1) return false;
    if(tjs.matchVec.empty()) return false;
    
    int theTj = pfp.TjIDs[itj];
//    bool pfpInShower = (pfp.PDGCode == 11);
    
    unsigned short oldSize = pfp.TjIDs.size();
    
    for(unsigned int ims = 0; ims < tjs.matchVec.size(); ++ims) {
      auto& ms = tjs.matchVec[ims];
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
        if(PFPVxTjOK(tjs, pfp, prt)) continue;
        pfp.TjCompleteness.push_back(0);
        if(prt) mf::LogVerbatim("TC")<<"AMT: P"<<pfp.ID<<" T"<<theTj<<" Add T"<<tjid;
      } // mtjid
    } // ims
    if(pfp.TjIDs.size() > oldSize) return true;
    return false;
  } // AddMissedTj

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
    
    // see if anything needs to be done
    std::vector<unsigned short> cntInPln(nplanes);
    bool itsOK = true;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
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
      for(auto tjid : pfp.TjIDs) if(tjs.allTraj[tjid - 1].CTP == inCTP) tjids.push_back((unsigned short)tjid);
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
        auto tj = tjs.allTraj[firstPt[0] - 1];
        unsigned short midPt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
        if(firstPt[1] > midPt) ReverseTraj(tjs, tj);
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
      SetEndPoints(tjs, mtj);
      mtj.MCSMom = MCSMom(tjs, mtj);
      SetPDGCode(tjs, mtj);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" P"<<pfp.ID<<" try to merge";
        for(auto tjid : tjids) {
          auto& tj = tjs.allTraj[tjid - 1];
          myprt<<" T"<<tjid<<" MCSMom "<<tj.MCSMom;
        }
        myprt<<" -> T"<<tjs.allTraj.size() + 1;
        myprt<<" MCSMom "<<mtj.MCSMom;
      }
      // kill the broken tjs and update the pfp TP3s
      for(auto tjid : tjids) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.AlgMod[kInShower]) mtj.AlgMod[kInShower] = true;
        MakeTrajectoryObsolete(tjs, tjid - 1);
      }
      // save the new one
      if(!StoreTraj(tjs, mtj)) {
        std::cout<<"MergePFPTjs: StoreTraj failed P"<<pfp.ID<<" EventsProcessed "<<tjs.EventsProcessed<<"\n";
        return false;
      }
      int newTjID = tjs.allTraj.size();
      newTjIDs.push_back(newTjID);
      // prepare to clobber vertices
      std::vector<unsigned short> vxlist;
      for(auto tjid : tjids) {
        // update the stored match struct and Tp3s
        UpdateMatchStructs(tjs, tjid, newTjID);
        // Update the Tp3s of this pfp
        UpdateTp3s(tjs, pfp, tjid, newTjID);
        auto& tj = tjs.allTraj[tjid - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == 0) continue;
          if(std::find(vxlist.begin(), vxlist.end(), tj.VtxID[end]) != vxlist.end()) {
            auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
//            std::cout<<"P"<<pfp.ID<<" Clobber 2V"<<vx2.ID<<"\n";
            MakeVertexObsolete(tjs, vx2, true);
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
      PrintPFP("MPTJ", tjs, pfp, true);
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
        // InShower compatibility
        if(jTjPt.inShower != iTjPt.inShower) continue;
        // length cut
        if(jTjPt.npts < minPts) continue;
        // score cut
        if(jTjPt.score < 0 || jTjPt.score > maxScore) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtp = tjs.allTraj[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
        TrajPoint3 tp3;
        if(!MakeTp3(tjs, itp, jtp, tp3, false)) continue;
        // count weight is one for a two-plane match
        float cntWght = 1;
        if(numPlanes == 3) {
          // numPlanes == 3
          for(unsigned int kpt = jpt + 1; kpt < tjs.mallTraj.size(); ++kpt) {
            auto& kTjPt = tjs.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
            // InShower compatibility - not in the third plane
//            if(kTjPt.inShower != iTjPt.inShower) continue;
            if(kTjPt.score < 0 || kTjPt.score > maxScore) continue;
            if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 5) break;
            auto& ktp = tjs.allTraj[kTjPt.id - 1].Pts[kTjPt.ipt];
            unsigned short kplane = DecodeCTP(ktp.CTP).Plane;
            TrajPoint3 iktp3;
            if(!MakeTp3(tjs, itp, ktp, iktp3, false)) continue;
            if(std::abs(tp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
            if(std::abs(tp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
            float dang = 0;
            if(tjs.Match3DCuts[1] > 0) {
              dang = std::abs(DeltaAngle(tp3.Dir, iktp3.Dir));
              while(dang >  M_PI) dang -= twopi;
              if(dang >  piOver2) dang = M_PI - dang;
              float mcsmom = tjs.allTraj[iTjPt.id - 1].MCSMom + tjs.allTraj[jTjPt.id - 1].MCSMom + tjs.allTraj[kTjPt.id - 1].MCSMom;
              mcsmom /= 3;
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
/*
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
*/
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
  bool MakeTp3(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3, bool findDirection)
  {
    // Make a 3D trajectory point using two 2D trajectory points
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    tp3.Dir = {{999.0, 999.0, 999.0}};
    tp3.Pos = {{999.0, 999.0, 999.0}};
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
  float PFPDOCA(const PFPStruct& pfp1,  const PFPStruct& pfp2, unsigned short& close1, unsigned short& close2)
  {
    // returns the Distance of Closest Approach between two PFParticles. 
    
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
  bool Split3DKink(TjStuff& tjs, PFPStruct& pfp, double sep, bool prt)
  {
    // Finds kinks in the PFParticle, splits Tjs, creates 2D vertices and forces a rebuild if any are found
    if(pfp.Tp3s.empty()) return false;
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
  PFPStruct CreatePFP(const TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // The calling function should define the size of pfp.TjIDs
    PFPStruct pfp;
    pfp.ID = tjs.pfps.size() + 1;
    // assume it is it's own parent
    pfp.ParentID = pfp.ID;
    pfp.TPCID = tpcid;
    // initialize arrays for both ends
    if(tjs.NumPlanes < 4) {
      pfp.dEdx[0].resize(tjs.NumPlanes, 0);
      pfp.dEdx[1].resize(tjs.NumPlanes, 0);
      pfp.dEdxErr[0].resize(tjs.NumPlanes, 0);
      pfp.dEdxErr[1].resize(tjs.NumPlanes, 0);
    }
    for(unsigned short startend = 0; startend < 2; ++startend) {
      // BUG the double brace syntax is required to work around clang bug 21629
      // (https://bugs.llvm.org/show_bug.cgi?id=21629)
      pfp.Dir[startend] = {{0.0, 0.0, 0.0}};
      pfp.DirErr[startend] = {{0.0, 0.0, 0.0}};
      pfp.XYZ[startend] = {{0.0, 0.0, 0.0}};
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
    // clear the kEnvFlag bits on all Tjs. The bit will be set true when a TP is
    // used in a PFParticle Tp3
    for(auto& tj : tjs.allTraj) {
      for(auto& tp : tj.Pts) tp.Environment[kEnvFlag] = false;
    } // tj
    
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
    
    // sort by decreasing number of matched points
    if(matVec.size() > 1) {
      PFPStruct pfp = CreatePFP(tjs, tpcid);
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
    PFPStruct pfp = CreatePFP(tjs, tpcid);
    for(auto& ms : matVec) {
      if(ms.Count < 2) continue;
      // check for duplicates
      bool skipit = false;
      for(auto& oms : tjs.matchVec) {
        if(ms.TjIDs == oms.TjIDs) {
          skipit = true;
          break;
        }
      } // oms
      if(skipit) continue;
      // Find completeness, do the fit, don't fill Tp3s, don't print
      pfp.TjIDs = ms.TjIDs;
      FindCompleteness(tjs, pfp, true, false, false);
      // save the info in matchStruct
      ms.TjCompleteness = pfp.TjCompleteness;
      ms.Pos = pfp.XYZ[0];
      ms.Dir = pfp.Dir[0];
      tjs.matchVec.push_back(ms);
    }
    if(tjs.matchVec.empty()) return;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FMV: tjs.matchVec\n";
      unsigned short cnt = 0;
      PFPStruct pfp = CreatePFP(tjs, tpcid);
      std::vector<int> dum1;
      std::vector<float> dum2;
      for(unsigned int ii = 0; ii < tjs.matchVec.size(); ++ii) {
        auto& ms = tjs.matchVec[ii];
        if(ms.Count == 0) continue;
        myprt<<std::setw(4)<<ii<<" Count "<<std::setw(5)<<(int)ms.Count;
        myprt<<" TjIDs:";
        for(auto& tjid : ms.TjIDs) myprt<<" T"<<tjid;
        myprt<<" Comp ";
        for(unsigned short itj = 0; itj < ms.TjCompleteness.size(); ++itj) {
          myprt<<std::setprecision(2)<<std::setw(6)<<ms.TjCompleteness[itj];
        }
        myprt<<" Pos ("<<std::setprecision(0)<<std::fixed;
        myprt<<ms.Pos[0]<<", "<<ms.Pos[1]<<", "<<ms.Pos[2];
        myprt<<") Dir "<<std::setprecision(2)<<std::setw(6)<<ms.Dir[0]<<std::setw(6)<<ms.Dir[1]<<std::setw(6)<<ms.Dir[2];
        myprt<<" projInPlane";
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(tjs, ms.Pos, ms.Dir, inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(tjs, ms.TjIDs);
        myprt<<" MCSMom "<<MCSMom(tjs, ms.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(tjs, ms.TjIDs, false);
        myprt<<"\n";
        ++cnt;
        if(cnt == 1000 || ms.Count < 2) {
          myprt<<"...stopped printing after 500 entries or Count < 2";
          break;
        }
      } // ii
    } // prt

    // create the list of associations to matches that will be converted to PFParticles
    // Start with Tjs attached to 3D vertices. This is only done when reconstructing neutrino events
    if(!tjs.TestBeam) {
      Match3DVtxTjs(tjs, tpcid, prt);
    }

    // define the PFParticleList
    for(unsigned int indx = 0; indx < tjs.matchVec.size(); ++indx) {
      auto& ms = tjs.matchVec[indx];
      // ignore dead matches
      if(ms.Count == 0) continue;
      // skip this match if any of the trajectories is already matched or merged and killed
      bool skipit = false;
      // check for muons and delta rays or InShower Tjs
      bool has13 = false;
      bool has11 = false;
      for(unsigned short itj = 0; itj < ms.TjIDs.size(); ++itj) {
        auto& tj = tjs.allTraj[ms.TjIDs[itj] - 1];
        if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) skipit = true;
        // skip low TjCompleteness
        if(ms.TjCompleteness[itj] < 0.1) skipit = true;
        if(tj.PDGCode == 13) has13 = true;
        if(tj.PDGCode == 11) has11 = true;
      } // tjID
      if(skipit) continue;
      if(has13 && has11) continue;
      int pdgCode = PDGCodeVote(tjs, ms.TjIDs, prt);
//      if(pdgCode != 11 && maxTjLen > 10 && ms.AspectRatio < tjs.Match3DCuts[3]) continue;
      PFPStruct pfp = CreatePFP(tjs, tpcid);
      pfp.TjIDs = ms.TjIDs;
      pfp.MatchVecIndex = indx;
      // Set the PDGCode so DefinePFP can ignore incompatible matches
      pfp.PDGCode = pdgCode;
      // set the PDGCode to ensure that delta rays aren't merged with muons. PDGCodeVote
      // returns 0 if the vote is mixed
      if(has13) pfp.PDGCode = 13;
      if(has11) pfp.PDGCode = 11;
      if(!DefinePFP("FPFP", tjs, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" DefinePFP failed";
        if(prt) std::cout<<" DefinePFP P"<<pfp.ID<<" failed\n";
        continue;
      }
      double sep = 1;
      Split3DKink(tjs, pfp, sep, prt);
      AnalyzePFP(tjs, pfp, prt);
      if(!StorePFP(tjs, pfp)) {
        if(prt) mf::LogVerbatim("TC")<<" StorePFP failed "<<pfp.ID;
        if(prt) std::cout<<" StorePFP failed "<<pfp.ID<<"\n";
        continue;
      }
      ms.Count = 0;
      // clobber MatchStructs that use the Tjs in this pfp
      for(auto& allms : tjs.matchVec) {
        auto shared = SetIntersection(allms.TjIDs, pfp.TjIDs);
        if(!shared.empty()) allms.Count = 0;
      } // allms
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
    
    // require at least one tj in at least two planes
    std::vector<unsigned short> nInPln(tjs.NumPlanes);
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      ++nInPln[DecodeCTP(tj.CTP).Plane];
      if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) {
        std::cout<<fcnLabel<<" pfp "<<pfp.ID<<" uses tj T"<<tj.ID<<" but kMat3D is set true\n";
        return false;
      }
    } // tjid
    unsigned short npl = 0;
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) if(nInPln[plane] > 0) ++npl;
    if(npl < 2) return false;

    if(pfp.Vx3ID[0] == 0 && pfp.Vx3ID[1] > 0) {
      std::cout<<fcnLabel<<" pfp P"<<pfp.ID<<" end 1 has a vertex but end 0 doesn't. No endpoints defined\n";
      return false;
    }
    
    // check for vertex consistency. There should only be one tj in each plane
    // that is attached to a vertex. Remove the shorter tj from the list
    // if that is not the case 
    if(pfp.Vx3ID[0] > 0) PFPVxTjOK(tjs, pfp, prt);
    
    bool pfpTrackLike = (MaxTjLen(tjs, pfp.TjIDs) > tjs.Match3DCuts[5] && MCSMom(tjs, pfp.TjIDs) > tjs.Match3DCuts[3]);

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" Vx3ID "<<pfp.Vx3ID[0]<<" "<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
      myprt<<" matchVec index "<<pfp.MatchVecIndex;
      myprt<<" max Tj len "<<MaxTjLen(tjs, pfp.TjIDs);
      myprt<<" MCSMom "<<MCSMom(tjs, pfp.TjIDs);
      myprt<<" pfpTrackLike? "<<pfpTrackLike;
    } // prt
    
    if(tjs.UseAlg[kMat3DMerge] && pfpTrackLike && pfp.MatchVecIndex < tjs.matchVec.size()) {
      // The index of tjs.matchVec has been specified for this pfp so we can look for evidence of
      // broken Tjs starting at the beginning
      for(unsigned short ims = 0; ims < pfp.MatchVecIndex + 10; ++ims) {
        if(ims >= tjs.matchVec.size()) break;
        auto& ms = tjs.matchVec[ims];
        if(ms.Count == 0) continue;
        std::vector<int> shared = SetIntersection(ms.TjIDs, pfp.TjIDs);
        if(shared.size() < 2) continue;
        // check the max length Tj and cut on MCSMom
        if(MaxTjLen(tjs, ms.TjIDs) < tjs.Match3DCuts[5]) continue;
        if(MCSMom(tjs, ms.TjIDs) < tjs.Match3DCuts[3]) continue;
        for(auto tjid : ms.TjIDs) {
          if(std::find(shared.begin(), shared.end(), tjid) != shared.end()) continue;
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kKilled]) continue;
          if(tj.AlgMod[kMat3D]) continue;
          // check for PDGCode compatibility - muons and delta rays
          if(pfp.PDGCode == 13 && tj.PDGCode == 11) continue;
          if(pfp.PDGCode == 11 && tj.PDGCode == 13) continue;
          if(SharesHighScoreVx(tjs, pfp, tj)) continue;
          if(tj.MCSMom < tjs.Match3DCuts[3]) continue;
          float len = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
          if(len < tjs.Match3DCuts[5]) continue;
          // check for a compatible merge
          bool skipit = false;
          for(auto tjid : pfp.TjIDs) {
            auto& ptj = tjs.allTraj[tjid - 1];
            if(ptj.CTP != tj.CTP) continue;
            if(!CompatibleMerge(tjs, tj, ptj, prt)) {
              skipit = true;
              break;
            }
          } // tjid
          if(skipit) continue;
          if(prt) mf::LogVerbatim("TC")<<" add T"<<tjid<<" MCSMom "<<tj.MCSMom<<" length "<<len;
          pfp.TjIDs.push_back(tjid);
          PFPVxTjOK(tjs, pfp, prt);
        } // tjid
      } // ims
    } // matchVec index defined
    
    // check the completeness of matching points in this set of Tjs and possibly
    // merge Tjs
    if(tjs.UseAlg[kMat3DMerge] && pfpTrackLike) {
      if(!CheckAndMerge(tjs, pfp, prt)) return false;
    } else {
      // not track-like. Find completeness and fill the TP3s
      FindCompleteness(tjs, pfp, true, true, prt);
    }

    // Set the starting position in 3D if it isn't already defined by a 3D vertex
    SetNewStart(tjs, pfp, prt);
    FollowTp3s(tjs, pfp, prt);
    // Check the list of Tjs and merge those that are in the same plane
    MergePFPTjs(tjs, pfp, prt);

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" Vx3ID "<<pfp.Vx3ID[0]<<" "<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
      myprt<<" Tp3s size "<<pfp.Tp3s.size();
    }
    
    pfp.NeedsUpdate = false;
    return true;
  } // DefinePFP
  
  /////////////////////////////////////////
  bool PFPVxTjOK(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Checks the PFP Vx3 -> Vx2 -> Tj assignment to see if there is more
    // than one tj in a plane in the pfp.TjIDs list that is attached to the same 2D vertex. 
    // This problem is fixed by removing the shorter tj from the TjIDs list. This function
    // return true if nothing was done to TjIDs
    if(pfp.ID == 0) return true;
    if(pfp.TjIDs.empty()) return true;
    if(pfp.Vx3ID[0] == 0) return true;
    
    auto& vx3 = tjs.vtx3[pfp.Vx3ID[0] - 1];
    std::vector<int> killMe;
    for(auto vx2id : vx3.Vx2ID) {
      if(vx2id == 0) continue;
      auto& vx2 = tjs.vtx[vx2id - 1];
      auto tjlist = GetVtxTjIDs(tjs, vx2);
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
        auto& tj = tjs.allTraj[tid - 1];
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
      myprt<<"PVTC: P"<<pfp.ID<<" removing short tjs attached to a vertex with a longer tj:";
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
  void AnalyzePFP(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Analyzes the PFP for oddities and tries to fix them
    if(pfp.ID == 0) return;
    if(pfp.TjIDs.empty()) return;
    
    // don't bother analyzing this pfp has been altered
    if(pfp.NeedsUpdate) {
      if(prt) mf::LogVerbatim("TC")<<"AnalyzePFP: P"<<pfp.ID<<" needs to be updated. Skip analysis ";
      return;
    }
    if(prt) mf::LogVerbatim("TC")<<"inside AnalyzePFP P"<<pfp.ID<<" NeedsUpdate? "<<pfp.NeedsUpdate;
    
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
      auto& missTj = tjs.allTraj[tjIDs[ii] - 1];
      if(missTj.AlgMod[kMat3D]) continue;
      unsigned short npwc = NumPtsWithCharge(tjs, missTj, false);
      if(prt) mf::LogVerbatim("TC")<<" missed T"<<missTj.ID<<" npwc "<<npwc<<" tjCnt "<<tjCnt[ii];
      if(tjCnt[ii] < 0.5 * npwc) continue;
      // add the missed Tj to the pfp and flag it as needing an update
      pfp.TjIDs.push_back(missTj.ID);
      if(PFPVxTjOK(tjs, pfp, prt)) pfp.NeedsUpdate = true;
    } // ii
    
    if(pfp.NeedsUpdate) DefinePFP("APFP", tjs, pfp, prt);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"APFP: Tjs in pfp\n";
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        myprt<<"T"<<tj.ID<<" npwc "<<NumPtsWithCharge(tjs, tj, false);
        unsigned short indx = 0;
        for(indx = 0; indx < tjIDs.size(); ++indx) if(tjIDs[indx] == tjid) break;
        if(indx == tjIDs.size()) {
          myprt<<" not found in P"<<pfp.ID<<"\n";
          continue;
        }
        myprt<<" nTp3 "<<tjCnt[indx]<<"\n";
      } // tjid
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
      // ignore truth photons
      if(pfp.PDGCode == 22) continue;
      if(pfp.Vx3ID[1] == 0 && !pfp.Tp3s.empty()) {
        // See if the direction needs to be changed
        std::cout<<"PFPVertexCheck needs revision\n";
//        SetNewStart(tjs, pfp, false);
      }
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = false;
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
      // ignore truth photon
      if(pfp.PDGCode == 22) continue;
      if(pfp.Vx3ID[0] > 0) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = false;
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
        if(cnt3 != tjs.NumPlanes) {
          mf::LogVerbatim("TC")<<"DPFPR: Missed an end vertex for PFP "<<pfp.ID<<" Write some code";
        }
      } // cnt3 > 1
    } // pfp
    
    // Assign a PDGCode to each PFParticle and look for a parent
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != tpcid) continue;
      // skip a neutrino PFParticle
      if(pfp.PDGCode == 12 || pfp.PDGCode == 14 || pfp.PDGCode == 22) continue;
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
  } // DefinePFPParentsTestBeam

  ////////////////////////////////////////////////
  bool StorePFP(TjStuff& tjs, PFPStruct& pfp)
  {
    // stores the PFParticle in TJStuff
    if(pfp.ID < tjs.pfps.size()) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    if(!neutrinoPFP) {
      if(pfp.TjIDs.empty()) return false;
      if(pfp.Tp3s.size() < 2) return false;
    }
    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)tjs.pfps.size() + 1) pfp.ID = tjs.pfps.size() + 1;
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
  void PrintTp3(std::string someText, const TjStuff& tjs, const TrajPoint3& tp3)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" Pos"<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(6)<<tp3.Pos[0]<<std::setw(6)<<tp3.Pos[1]<<std::setw(6)<<tp3.Pos[2];
    myprt<<" Dir"<<std::fixed<<std::setprecision(3);
    myprt<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
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
      myprt<<someText<<" pfp P"<<pfp.ID<<"\n";
      myprt<<someText<<"  ipt ________Pos________ Path   _______Dir______ dE/dx dang  Kink  Tj_ipt \n";
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
      myprt<<std::setprecision(1)<<std::setw(6)<<tp3.dEdx;
      myprt<<std::setprecision(3)<<std::setw(7)<<DeltaAngle(prevDir, tp3.Dir);
      prevDir = tp3.Dir;
      // Calculate the kink angle at point ipt, using the two points that are
      // +/- 1 cm on either side of that point
      double sep = 1;
      myprt<<std::setprecision(2)<<std::setw(7)<<KinkAngle(tjs, pfp.Tp3s, ipt, sep);
      for(auto tj2pt : tp3.Tj2Pts) {
        auto& tj = tjs.allTraj[tj2pt.id - 1];
        auto& tp = tj.Pts[tj2pt.ipt];
        myprt<<" "<<tj.ID<<"_"<<tj2pt.ipt<<"_"<<PrintPos(tjs, tp);
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
