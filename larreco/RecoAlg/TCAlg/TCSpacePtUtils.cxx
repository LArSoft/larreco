#include "larreco/RecoAlg/TCAlg/TCSpacePtUtils.h"

namespace tca {

  struct SortEntry{
    unsigned int index;
    float val;
  };
  // TODO: Fix the sorting mess
  bool valDecreasings (SortEntry c1, SortEntry c2) { return (c1.val > c2.val);}
  bool valIncreasings (SortEntry c1, SortEntry c2) { return (c1.val < c2.val);}
  
  //////////////////////////////////////////
  void Match3DSpts(TjStuff& tjs, const art::Event& evt, const geo::TPCID& tpcid, 
                   const art::InputTag& fSpacePointModuleLabel)
  {
    // Match Tjs in 3D using SpacePoint <-> Hit associations.
    if(!tjs.UseAlg[kHitsOrdered]) {
//      std::cout<<"Match3DSpts: Write some code to deal with un-ordered hits\n";
      return;
    }
    
    bool prt = (debug.Plane >= 0) && (debug.Tick == 3333);
    if(prt) mf::LogVerbatim("TC")<<"Inside Match3DSpts";
    
    // sort the Tjs in this tpc by decreasing length
    std::vector<SortEntry> sortVec;
    SortEntry se;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(DecodeCTP(tj.CTP).Cryostat != tpcid.Cryostat) continue;
      if(DecodeCTP(tj.CTP).TPC != tpcid.TPC) continue;
      if(tj.AlgMod[kInShower]) continue;
      se.index = tj.ID;
      se.val = tj.EndPt[1] - tj.EndPt[0];
      sortVec.push_back(se);
    } // tj
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
    
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      int tjid = sortVec[ii].index;
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kMat3D]) continue;
      unsigned short npts = tj.EndPt[1] - tj.EndPt[0] + 1;
      // try to find a space point at each Tj point
      std::vector<std::vector<unsigned int>> sptLists(npts);
      // keep track of tj IDs at each point
      std::vector<std::vector<int>> sptTjLists(npts);
      // list of Tjs and the occurrence count
      std::vector<std::pair<int, int>> tjcnts;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        unsigned short indx = ipt - tj.EndPt[0];
        auto& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        // Find space points using hits associated with this Tp. 
        sptLists[indx] = SpacePtsAssociatedWith(tjs, tp);
        if(sptLists[indx].empty()) continue;
        // next get a list of Tjs that use hits that are in those space points
        auto tmp = TjsNearSpacePts(tjs, sptLists[indx]);
        for(auto tjid : tmp) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(!tj.AlgMod[kMat3D]) sptTjLists[indx].push_back(tjid);
        } // tjid
        for(auto tjid : sptTjLists[indx]) {
          unsigned short jndx = 0;
          for(jndx = 0; jndx < tjcnts.size(); ++jndx) if(tjid == tjcnts[jndx].first) break;
          if(jndx == tjcnts.size()) {
            // not in the list so add it
            tjcnts.push_back(std::make_pair(tjid, 1));
          } else {
            ++tjcnts[jndx].second;
          }
        } // tjid
      } // ipt
      if(tjcnts.empty()) continue;
      // put the list of matching tjs into plane ordered lists
      std::array<std::vector<int>, 3> tjInPln;
      // require that a broken Tj be not too short
      unsigned short minpts = 0.1 * (tj.EndPt[1] - tj.EndPt[0]);
      if(minpts < 2) minpts = 2;
      for(auto& tjcnt : tjcnts) {
        auto& mtj = tjs.allTraj[tjcnt.first - 1];
        if(mtj.EndPt[1] - mtj.EndPt[0] < minpts) continue;
        // ignore if killed or already matched
        if(mtj.AlgMod[kKilled]) continue;
        if(mtj.AlgMod[kMat3D]) continue;
        unsigned short plane = DecodeCTP(mtj.CTP).Plane;
        tjInPln[plane].push_back(tjcnt.first);
      }
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<tj.ID<<" tjInPln\n";
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          myprt<<" pln "<<plane<<" Tjs";
          for(auto tjid : tjInPln[plane]) myprt<<" "<<tjid;
        } // plane
      } // prt
      // Look for broken Tjs, that is the existence of more than one
      // tj that is matched to the space points in a plane
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(tjInPln[plane].size() < 2) continue;
        // TODO: Check for the existence of a vertex here and decide if
        // Tjs should split instead of merging. 
        if(!MergeBrokenTjs(tjs, tjInPln[plane], prt)) {
          if(prt) mf::LogVerbatim("TC")<<"Merge failed. Skipping this combination";
          continue;
        }
      } // plane
      // Ensure that there is an un-matched Tj in at least 2 planes
      unsigned short nPlnWithTj = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(tjInPln[plane].empty()) continue;
        bool tjAvailable = false;
        for(auto tjid : tjInPln[plane]) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(!tj.AlgMod[kMat3D]) tjAvailable = true;
        } // tjid
        if(tjAvailable) ++nPlnWithTj;
      } // plane
      if(nPlnWithTj < 2) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<tj.ID<<" After merging\n";
        for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
          myprt<<" pln "<<plane<<" Tjs";
          for(auto tjid : tjInPln[plane]) myprt<<" "<<tjid;
        } // plane
      } // prt
      PFPStruct pfp = CreatePFPStruct(tjs, tpcid);
      pfp.TjIDs.clear();
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(tjInPln[plane].empty()) continue;
        if(tjInPln[plane][0] > 0) {
          pfp.TjIDs.push_back(tjInPln[plane][0]);
          auto& tj = tjs.allTraj[tjInPln[plane][0] - 1];
          tj.AlgMod[kMat3D] = true;
        }
      } // plane
      // set the end points using the local version that uses SpacePoints instead
      // of tjs.malltraj
/*
      if(!SetPFPEndPoints(tjs, pfp, sptLists, tj.ID, prt)) {
        std::cout<<"SetPFPEndPoints failed";
      }
*/
      tjs.pfps.push_back(pfp);
      if(ii > 10) break;
    } // ii (tj)
    
  } // Match3DSpts

  /////////////////////////////////////////
  bool SetPFPEndPoints(TjStuff& tjs, PFPStruct& pfp, std::vector<std::vector<unsigned int>>& sptLists, int tjID, bool prt)
  {
    // Find PFParticle end points using the lists of spacepoints that are associated with each TP
    // on Trajectory tjID. This is the longest trajectory of the PFParticle
    
    // this code doesn't handle the ends of showers
    if(pfp.PDGCode == 1111) return false;
    if(pfp.TjIDs.size() < 2) return false;
    if(sptLists.empty()) return false;
    // this function needs spacepoints
    if(tjs.pfps.empty()) return false;
    unsigned short tjIndex = 0;
    for(tjIndex = 0; tjIndex < pfp.TjIDs.size(); ++tjIndex) if(tjID == pfp.TjIDs[tjIndex]) break;
    if(tjIndex == pfp.TjIDs.size()) return false;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"SPEP: PFP "<<pfp.ID;
//      myprt<<" Vx3ID "<<pfp.Vx3ID[end];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" "<<id;
      if(pfp.PDGCode == 1111) myprt<<" This is a shower PFP ";
    }
    
    // find the first point that has all Tjs in a space point. We will use this to decide
    // which ends of the Tjs match
    std::array<unsigned short, 2> tj2PtWithSpt {USHRT_MAX};
    for(unsigned short ipt = 0; ipt < sptLists.size(); ++ipt) {
      // sptList is a vector of all spacepoints that are associated with the trajectory point ipt
      auto tmp = TjsNearSpacePts(tjs, sptLists[ipt]);
      unsigned short cnt = 0;
      for(auto tjn : tmp) {
        if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjn) != pfp.TjIDs.end()) ++cnt;
      } // tjn
/*
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<ipt<<" tjn ";
        for(auto tjid : tmp) myprt<<" "<<tjid;
        myprt<<" cnt "<<cnt;
      }
*/
      if(cnt != pfp.TjIDs.size()) continue;
      if(tj2PtWithSpt[0] == USHRT_MAX) tj2PtWithSpt[0] = ipt;
      tj2PtWithSpt[1] = ipt;
    } // ipt
    // Do something quick and dirty
    for(unsigned short end = 0; end < 2; ++end) {
      unsigned short ipt = tj2PtWithSpt[end];
      if(ipt > sptLists.size() - 1) continue;
      unsigned int isp = sptLists[ipt][0];
      std::cout<<"ipt "<<ipt<<" "<<isp<<"\n";
      auto& spt = pfp.Tp3s[isp];
      pfp.XYZ[end] = spt.Pos;
    }
    if(prt) mf::LogVerbatim("TC")<<" firstPt "<<tj2PtWithSpt[0]<<" lastPt "<<tj2PtWithSpt[1];

    return true;
  } // SetPFPEndPoints
  
  /////////////////////////////////////////
  bool MergeBrokenTjs(TjStuff& tjs, std::vector<int>& tjInPln, bool prt)
  {
    // The vector tjInPln contains a list of Tjs that are in SpacePoints that are matched to
    // trajectories in the other planes. This function will merge them and correct the tjInPln vector
    
    if(tjInPln.size() < 2) return false;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MergeBrokenTjs working on";
      for(auto tjid : tjInPln) myprt<<" "<<tjid;
    }
    // put these into the same wire order
    for(auto tjid : tjInPln) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.StepDir != tjs.StepDir) ReverseTraj(tjs, tj);
    } // tjid
    // Put them into a consistent order - the step direction
    std::vector<SortEntry> sortVec(tjInPln.size());
    for(unsigned short ii = 0; ii < tjInPln.size(); ++ii) {
      auto& tj = tjs.allTraj[tjInPln[ii] - 1];
      if(tj.AlgMod[kKilled]) return false;
      if(tj.AlgMod[kMat3D]) return false;
      auto& tp0 = tj.Pts[tj.EndPt[0]];
      sortVec[ii].index = ii;
      sortVec[ii].val = tp0.Pos[0];
    } // ii
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put them into the right order
    auto temp = tjInPln;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) tjInPln[ii] = temp[sortVec[ii].index];
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"After sort";
      for(auto tjid : tjInPln) myprt<<" "<<tjid;
    }
    // merge them
    int mergedTjID = tjInPln[0];
    for(unsigned short nit = 1; nit < tjInPln.size(); ++nit) {
      auto& tj1 = tjs.allTraj[mergedTjID - 1];
      auto& tj2 = tjs.allTraj[tjInPln[nit] - 1];
      // check for a compatible merge
      if(!CompatibleMerge(tjs, tj1, tj2, prt)) continue;
      // check for a vertex between them and clobber it. The decision to break this into
      // two tracks should be done later using 3D angles
      if(tj1.VtxID[1] > 0 && tj1.VtxID[1] == tj2.VtxID[0]) {
        auto& vx2 = tjs.vtx[tj1.VtxID[1] - 1];
        MakeVertexObsolete(tjs, vx2, true);
      }
      // check for Bragg peaks and clobber them
      if(tj1.StopFlag[1][kBragg]) {
        if(prt) mf::LogVerbatim("TC")<<" clobbered Bragg Stop flag on Tj "<<tj1.ID;
        tj1.StopFlag[1][kBragg] = false;
      }
      if(tj2.StopFlag[0][kBragg]) {
        if(prt) mf::LogVerbatim("TC")<<" clobbered Bragg Stop flag on Tj "<<tj2.ID;
        tj2.StopFlag[0][kBragg] = false;
      }
      // TODO: The above actions should be done more carefully
      if(!MergeAndStore(tjs, tj1.ID - 1, tj2.ID - 1, prt)) return false;
      // Successfull merge
      mergedTjID = tjs.allTraj.size();
    } // nit
    // Revise tjs.mallTraj
    UpdateMatchStructs(tjs, tjInPln, mergedTjID);
    tjInPln.resize(1);
    tjInPln[0] = mergedTjID;
    if(prt) mf::LogVerbatim("TC")<<"MergeBrokenTjs returns "<<mergedTjID;
    return true;
  } // MergeBrokenTjs
  
  /////////////////////////////////////////
  std::vector<int> TjsNearSpacePts(TjStuff& tjs, std::vector<unsigned int> sptlist)
  {
    // Returns a list of Tjs that use hits that are associated with space
    // points in the list
    std::vector<int> tjlist;
    if(sptlist.empty()) return tjlist;
    for(auto isp : sptlist) {
      if(isp > tjs.spts.size() - 1) continue;
      auto& spt = tjs.spts[isp];
      for(auto iht : spt.Hits) {
        if(iht > tjs.fHits.size() - 1) continue;
        auto& hit = tjs.fHits[iht];
        if(hit.InTraj <= 0) continue;
        if(std::find(tjlist.begin(), tjlist.end(), hit.InTraj) == tjlist.end()) tjlist.push_back(hit.InTraj);
      } // iht
    } // isp
    if(tjlist.size() > 1) std::sort(tjlist.begin(), tjlist.end());
    return tjlist;
  } // TjsNearSpacePts

  /////////////////////////////////////////
  std::vector<unsigned int> SpacePtsAssociatedWith(TjStuff& tjs, const Trajectory& tj)
  {
    // returns a list of space point indices associated with the trajectory
    std::vector<unsigned int> allSptList;
    if(tj.AlgMod[kKilled]) return allSptList;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      auto& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      auto sptlist = SpacePtsAssociatedWith(tjs, tp);
      for(auto isp : sptlist) if(std::find(allSptList.begin(), allSptList.end(), isp) == allSptList.end()) allSptList.push_back(isp);
    } // ipt
    return allSptList;
  } // SpacePtsAssociatedWith

  /////////////////////////////////////////
  std::vector<unsigned int> SpacePtsAssociatedWith(TjStuff& tjs, const TrajPoint& tp)
  {
    // returns a list of space point indices associated with the trajectory point
    std::vector<unsigned int> allSptList;
    if(tp.Chg <= 0) return allSptList;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      unsigned int iht = tp.Hits[ii];
      auto sptlist = SpacePtsAtHit(tjs, iht);
      if(sptlist.empty()) continue;
      for(auto isp : sptlist) if(std::find(allSptList.begin(), allSptList.end(), isp) == allSptList.end()) allSptList.push_back(isp);
    } // ii
    return allSptList;
  } // SpacePtsAssociatedWith
  
  /////////////////////////////////////////
  std::vector<unsigned int> SpacePtsAtHit(TjStuff& tjs, unsigned int iht)
  {
    // returns a vector of indices into the tjs.Tp3s vector of SpacePoints that
    // contain hit with index iht
    std::vector<unsigned int> temp;
    if(tjs.spts.empty()) return temp;
    for(unsigned int isp = 0; isp < tjs.spts.size(); ++isp) {
      auto& spt = tjs.spts[isp];
      for(auto hit : spt.Hits) {
        if(hit == iht) {
          temp.push_back(isp);
          break;
        }
      } // hit
    } // isp
    return temp;
  } // SpacePtsAtHit
  
  /////////////////////////////////////////
  void UpdateMatchStructs(TjStuff& tjs, int tjid)
  {
    // Update the match structs when a single Tj has been reversed
    if(tjs.mallTraj.empty() && tjs.pfps.empty()) return;
    
    std::vector<int> tjids(1);
    tjids[0] = tjid;
    UpdateMatchStructs(tjs, tjids, tjid);
  } // UpdateMatchStructs
  
  /////////////////////////////////////////
  void UpdateMatchStructs(TjStuff& tjs, std::vector<int> oldTjs, int newTj)
  {
    // Updates the points in the match structs that reference the list of oldTjs with the newTj.
    // This is called after merging the oldTjs into the newTj.
    // It is assumed that the oldTjs have been or will be killed elsewhere
    
    if(oldTjs.empty()) return;
    if(newTj <= 0 || newTj > (int)tjs.allTraj.size()) return;
    if(tjs.mallTraj.empty() && tjs.pfps.empty()) return;
    
    auto& ntj = tjs.allTraj[newTj - 1];
    
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;
    for(auto tjid : oldTjs) {
      if(tjid <= 0 || tjid > (int)tjs.allTraj.size()) continue;
      if(!tjs.mallTraj.empty()) {
        for(unsigned int ipt = 0; ipt < tjs.mallTraj.size(); ++ipt) {
          auto& tj2pt = tjs.mallTraj[ipt];
          if(tj2pt.id != tjid) continue;
          // look for the corresponding point (wire) on the new Tj
          for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
            // TODO: Also check Pos[1]?
            if(std::nearbyint(ntj.Pts[npt].Pos[0]) == tj2pt.wire) {
              tj2pt.id = newTj;
              tj2pt.ipt = npt;
              tj2pt.npts = npts;
              break;
            }
          } // npt
        } // ipt
      } // mallTraj exists
      // Update pfp space points
      if(!tjs.pfps.empty()) {
        for(auto& pfp : tjs.pfps) {
          for(auto& tp3 : pfp.Tp3s) {
            // check each of the Tj2Pts associated with this space point
            for(auto& tj2pt : tp3.Tj2Pts) {
              if(tj2pt.id != tjid) continue;
              // look for the corresponding point (wire) on the new Tj
              for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
                if(std::nearbyint(ntj.Pts[npt].Pos[0]) == tj2pt.wire) {
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
    } // tjid (otj)
    
  } // UpdateMatchStructs
  
  /////////////////////////////////////////
  void FillmAllTraj(TjStuff& tjs, const geo::TPCID& tpcid) 
  {
    // Fills the tjs.mallTraj vector with trajectory points in the tpc
    tjs.matchVec.clear();
    
    int cstat = tpcid.Cryostat;
    int tpc = tpcid.TPC;
    
    // count the number of TPs and clear out any old 3D match flags
    unsigned int ntp = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      // don't match InShower Tjs
      if(tj.AlgMod[kInShower]) continue;
      // or Shower Tjs
      if(tj.AlgMod[kShowerTj]) continue;
      if(tj.ID <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      ntp += NumPtsWithCharge(tjs, tj, false);
      tj.AlgMod[kMat3D] = false;
    } // tj
    if(ntp < 2) return;
    
    tjs.mallTraj.resize(ntp);
    std::vector<SortEntry> sortVec(ntp);
    
    // define mallTraj
    unsigned int icnt = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      // don't match shower-like Tjs
      if(tj.AlgMod[kInShower]) continue;
      // or Shower Tjs
      if(tj.AlgMod[kShowerTj]) continue;
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
        if(icnt == tjs.mallTraj.size()) {
          std::cout<<"Match3D: indexing error\n";
          break;
        }
        tjs.mallTraj[icnt].xlo = xpos - rms;
        tjs.mallTraj[icnt].xhi = xpos + rms;
        tjs.mallTraj[icnt].dir = tp.Dir;
        tjs.mallTraj[icnt].ctp = tp.CTP;
        tjs.mallTraj[icnt].id = tjID;
        tjs.mallTraj[icnt].ipt = ipt;
        tjs.mallTraj[icnt].npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        tjs.mallTraj[icnt].score = score;
        if(tj.PDGCode == 11) {
          tjs.mallTraj[icnt].showerlike = true;
        } else {
          tjs.mallTraj[icnt].showerlike = false;
        }
        // populate the sort vector
        sortVec[icnt].index = icnt;
        sortVec[icnt].val = tjs.mallTraj[icnt].xlo;
        ++icnt;
      } // tp
    } // tj
    
    if(icnt < tjs.mallTraj.size()) {
      tjs.mallTraj.resize(icnt);
      sortVec.resize(icnt);
    }
    
    // sort by increasing xlo
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put tjs.mallTraj into sorted order
    auto tallTraj = tjs.mallTraj;
    for(unsigned int ii = 0; ii < sortVec.size(); ++ii) tjs.mallTraj[ii] = tallTraj[sortVec[ii].index];
    
  } // FillmAllTraj

  /////////////////////////////////////////
  void MakePFPSpacePts(TjStuff& tjs, PFPStruct& pfp, bool anyTj)
  {
    // Creates a vector of TrajCluster space points using 3D matches in mallTraj. If anyTj is set true
    // any Trajectory point in the 3rd plane that is consistent will be added as well, otherwise
    // only those Tjs that are in the pfp.TjIDs list are considered
    
    pfp.Tp3s.clear();
    if(pfp.ID == 0) return;
    if(pfp.TjIDs.empty()) return;
    if(tjs.mallTraj.empty()) return;
    double xcut = tjs.Match3DCuts[0];
    double yzcut = 1.5 * xcut;

    for(unsigned int ipt = 0; ipt < tjs.mallTraj.size() - 1; ++ipt) {
      auto& iTj2Pt = tjs.mallTraj[ipt];
      if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), iTj2Pt.id) == pfp.TjIDs.end()) continue;
      auto& itp = tjs.allTraj[iTj2Pt.id - 1].Pts[iTj2Pt.ipt];
      for(unsigned int jpt = ipt + 1; jpt < tjs.mallTraj.size(); ++jpt) {
        auto& jTj2Pt = tjs.mallTraj[jpt];
        // break out if we are well past the X matching region
        if(jTj2Pt.xlo > iTj2Pt.xhi + 5) break;
        if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), jTj2Pt.id) == pfp.TjIDs.end()) continue;
        // ensure that the planes are different
        if(jTj2Pt.ctp == iTj2Pt.ctp) continue;
        bool got2 = ((iTj2Pt.id == 2 && iTj2Pt.ipt < 5) || (jTj2Pt.id == 2 && jTj2Pt.ipt < 5));
        bool got6 = ((iTj2Pt.id == 6 && iTj2Pt.ipt < 5) || (jTj2Pt.id == 6 && jTj2Pt.ipt < 5));
        if(got2 && got6) {
          std::cout<<"Got2+6 ";
          std::cout<<ipt<<"_"<<iTj2Pt.id<<"_"<<iTj2Pt.ipt<<" "<<iTj2Pt.xlo<<" "<<iTj2Pt.xhi<<" ";
          std::cout<<jpt<<"_"<<jTj2Pt.id<<"_"<<jTj2Pt.ipt<<" "<<jTj2Pt.xlo<<" "<<jTj2Pt.xhi;
          std::cout<<"\n";
        }
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTj2Pt.xlo > iTj2Pt.xhi) continue;
        auto& jtp = tjs.allTraj[jTj2Pt.id - 1].Pts[jTj2Pt.ipt];
        Tp3Struct tcspt;
        if(!MakeSpt(tjs, itp, jtp, tcspt.Pos, tcspt.Dir)) continue;
        tcspt.Tj2Pts.resize(2);
        tcspt.Tj2Pts[0] = iTj2Pt;
        tcspt.Tj2Pts[1] = jTj2Pt;
        // try for a triple w/o the requirement for a Tj in the pfp
        if(tjs.NumPlanes == 3) {
          for(unsigned int kpt = jpt + 1; kpt < tjs.mallTraj.size(); ++kpt) {
            auto& kTj2Pt = tjs.mallTraj[kpt];
            if(kTj2Pt.xlo > iTj2Pt.xhi + 5) break;
            if(kTj2Pt.ctp == iTj2Pt.ctp) continue;
            if(kTj2Pt.ctp == jTj2Pt.ctp) continue;
            // Require the third point is on a Tj that is in the PFP list?
            if(!anyTj && std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), kTj2Pt.id) == pfp.TjIDs.end()) continue;
            if(kTj2Pt.xlo > iTj2Pt.xhi) continue;
            auto& ktp = tjs.allTraj[kTj2Pt.id - 1].Pts[kTj2Pt.ipt];
            Point3_t pos;
            Vector3_t dir;
            if(!MakeSpt(tjs, itp, ktp, pos, dir)) continue;
            if(std::abs(pos[1] - tcspt.Pos[1]) > yzcut) continue;
            if(std::abs(pos[2] - tcspt.Pos[2]) > yzcut) continue;
            // make a rough angle cut
            if(DotProd(tcspt.Dir, dir) < 0.8) continue;
            tcspt.Tj2Pts.push_back(kTj2Pt);
            // update the position and direction
            for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
              tcspt.Pos[ixyz] += pos[ixyz]; tcspt.Pos[ixyz] /= 2;
              tcspt.Dir[ixyz] += dir[ixyz]; tcspt.Dir[ixyz] /= 2;
            }
            // only allow one addition
            break;
          } // kpt
        } // 3 planes
        // sort the Tj2Pts by increasing Tj ID
        std::vector<SortEntry> sortVec(tcspt.Tj2Pts.size());
        for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
          sortVec[ii].index = ii;
          sortVec[ii].val = tcspt.Tj2Pts[ii].id;
        } // ii
        std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
        auto temp = tcspt.Tj2Pts;
        for(unsigned short ii = 0; ii < temp.size(); ++ii) tcspt.Tj2Pts[ii] = temp[sortVec[ii].index];
        tcspt.Chg = 0;
        for(auto tj2pt : tcspt.Tj2Pts) {
          auto& tp = tjs.allTraj[tj2pt.id - 1].Pts[tj2pt.ipt];
          tcspt.Chg += tp.Chg;
        } // mi
        tcspt.Chg /= (float)tcspt.Tj2Pts.size();
        pfp.Tp3s.push_back(tcspt);
        break;
      } // jpt
    } // ipt
    
  } // MakePFPSpacePts

  
  /////////////////////////////////////////
  bool SetNewStart(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Analyzes the space point collection and the Tjs in the pfp to create a start
    // vertex if one doesn't exist. This function returns true if the space points need to be sorted
    if(pfp.ID == 0 || pfp.TjIDs.empty()) return false;
    
    if(prt) mf::LogVerbatim("TC")<<"SNS: pfp "<<pfp.ID<<" Vx3ID[0] "<<pfp.Vx3ID[0];
    
    Point3_t minXYZ = {1E6, 1E6, 1E6};
    unsigned int minXYZPt[3] = {0, 0, 0};
    Point3_t maxXYZ = {-1E6, -1E6, -1E6};
    unsigned int maxXYZPt[3] = {0, 0, 0};
    for(unsigned short ipt = 0; ipt < pfp.Tp3s.size(); ++ipt) {
      auto& tcspt = pfp.Tp3s[ipt];
      auto& pos = tcspt.Pos;
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
        if(pos[ixyz] < minXYZ[ixyz]) {
          minXYZ[ixyz] = pos[ixyz];
          minXYZPt[ixyz] = ipt;
        }
        if(pos[ixyz] > maxXYZ[ixyz]) {
          maxXYZ[ixyz] = pos[ixyz];
          maxXYZPt[ixyz] = ipt;
        }
      } // ixyz
    } // tcspt
    
    // check for a fully contained track
    bool insideX = minXYZ[0] > tjs.XLo && maxXYZ[0] < tjs.XHi;
    bool insideY = minXYZ[1] > tjs.YLo && maxXYZ[1] < tjs.YHi;
    bool insideZ = minXYZ[2] > tjs.ZLo && maxXYZ[2] < tjs.ZHi;
    bool fullyContained = insideX && insideY && insideZ;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<" Xlo "<<tjs.XLo<<" minX "<<minXYZ[0]<<" maxX "<<maxXYZ[0]<<" XHi "<<tjs.XHi<<" insideX? "<<insideX<<"\n";
      myprt<<" Ylo "<<tjs.YLo<<" minY "<<minXYZ[1]<<" maxY "<<maxXYZ[1]<<" YHi "<<tjs.YHi<<" insideY? "<<insideY<<"\n";
      myprt<<" Zlo "<<tjs.ZLo<<" minZ "<<minXYZ[2]<<" maxZ "<<maxXYZ[2]<<" ZHi "<<tjs.ZHi<<" insideZ? "<<insideZ<<"\n";
    }
    
    // default is to do nothing. Start at the low X end
    unsigned short startSptIndex = 0;
    bool newStart = false;

    if(fullyContained) {
      // default is that this is beam-related. Start at minZ
      startSptIndex = minXYZPt[2];
      newStart = true;
    }  else {
      // not fully contained. Set the start to be the maximum Y position
      startSptIndex = maxXYZPt[1];
      newStart = true;
    } // not fully contained
    
    // TODO: insert code here to decide on the direction using delta-rays, etc
    
    if(prt) mf::LogVerbatim("TC")<<" fullyContained? "<<fullyContained<<" newStart? "<<newStart;
    
    // create a start vertex
    if(pfp.Vx3ID[0] == 0) {
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.Tp3s[startSptIndex].Pos[0];
      vx3.Y = pfp.Tp3s[startSptIndex].Pos[1];
      vx3.Z = pfp.Tp3s[startSptIndex].Pos[2];
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = true;
      tjs.vtx3.push_back(vx3);
      pfp.Vx3ID[0] = vx3.ID;
      pfp.XYZ[0] = pfp.Tp3s[startSptIndex].Pos;
      if(prt) mf::LogVerbatim("TC")<<" Created Vx_"<<vx3.ID;
    }
    return newStart;
    
  } // SetNewStart

  /////////////////////////////////////////
  void SortByDistanceFromStart(TjStuff& tjs, PFPStruct& pfp)
  {
    // sorts the TC spacepoints by distance from the start vertex pfp.XYZ[0]. This will result
    // in sorting by the distance from (0, 0, 0) if no start position has been specified
    
    if(pfp.Tp3s.size() < 2) return;
    
    Point3_t startPos = {0};
    if(pfp.Vx3ID[0] > 0) {
      auto& vx3 = tjs.vtx3[pfp.Vx3ID[0] - 1];
      startPos[0] = vx3.X;
      startPos[1] = vx3.Y;
      startPos[2] = vx3.Z;
    }
    mf::LogVerbatim("TC")<<"SBDFS: StartPos "<<std::fixed<<std::setprecision(1)<<startPos[0]<<" "<<startPos[1]<<" "<<startPos[2]<<"\n";

    std::vector<SortEntry> sortVec(pfp.Tp3s.size());
    for(unsigned short ii = 0; ii < pfp.Tp3s.size(); ++ii) {
      sortVec[ii].index = ii;
      float sep = 0;
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
        float arg = pfp.Tp3s[ii].Pos[ixyz] - startPos[ixyz];
        sep += arg * arg;
      } // ixyz
      sortVec[ii].val = sqrt(sep);
    } // ii

    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put them into order
    std::vector<Tp3Struct> temp;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) temp.push_back(pfp.Tp3s[sortVec[ii].index]);
    pfp.Tp3s = temp;
    
  } // SortByDistanceFromStart
  
  /////////////////////////////////////////
  void CheckSptValidity(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Checks and corrects the ordering of Tj points in the pfp space points. Checks for
    // consistency between the tjs listed in pfp.TjIDs and those appearing in the space points.
    // Sets the IsValid flag false for space points in the vector that are inconsistent
    // with their neighbors
    
    // ensure that the Tj points are in increasing order and reverse them if they aren't. This
    // presumes that the space points have been ordered from pfp start to pfp end
    std::vector<int> tjids;
    for(auto& spt : pfp.Tp3s) {
      for(auto& tj2pt : spt.Tj2Pts) {
        int tjid = tj2pt.id;
        // check for the first occurrence
        if(std::find(tjids.begin(), tjids.end(), tjid) != tjids.end()) continue;
        tjids.push_back(tjid);
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj2pt.ipt > tj.EndPt[0] + 3) {
          std::cout<<"pfp "<<pfp.ID<<" tj "<<tjid<<" first point not at end0 "<<tj2pt.ipt<<" tj SetDir "<<tj.AlgMod[kSetDir]<<"\n";
          if(tj.AlgMod[kSetDir]) {
            tj.AlgMod[kSetDir] = false;
            ReverseTraj(tjs, tj);
          }
        } // first tj point not at the beginning
        tj.AlgMod[kSetDir] = true;
      } // tjpt
    } // spt
    if(tjids.size() > 1) std::sort(tjids.begin(), tjids.end());
    if(tjids != pfp.TjIDs) {
      std::cout<<"CSV: Tjs are inconsistent";
      std::cout<<" pfp.TjIDs";
      for(auto tjid : pfp.TjIDs) std::cout<<" "<<tjid;
      std::cout<<" tjs in space points";
      for(auto tjid : tjids) std::cout<<" "<<tjid;
      std::cout<<"\n";
    }
    
    if(pfp.Tp3s.size() < 5) {
      return;
    } // short vector
    std::cout<<"write some code here\n";
    
  } // CheckSptValidity

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
    
    // create a temp vector to check for duplicates
    auto inMatVec = matVec;
    std::vector<MatchStruct> temp;
    
    // the minimum number of points for matching
    unsigned short minPts = 3;
    // override this with the user minimum for 2-plane matches
    if(numPlanes == 2) minPts = tjs.Match3DCuts[2];
    
    // max number of match combos left
    unsigned int nAvailable = 0;
    if(matVec.size() < tjs.Match3DCuts[4]) nAvailable = tjs.Match3DCuts[4] - matVec.size();
    if(nAvailable == 0 || nAvailable > tjs.Match3DCuts[4]) return;
    
    // temp TPs used to find 3D directions. The positions are not used
    TrajPoint tpi; tpi.Pos = {0, 0};
    TrajPoint tpj; tpj.Pos = {0, 0};
    TrajPoint tpk; tpk.Pos = {0, 0};
    // Direction vectors found using the i,j and i,k TPs
    TVector3 dij, dik, pos3;
    float piOver2 = M_PI / 2;
    // these cuts presume that the resolution in X is better than it is in Y and Z
    float xcut = tjs.Match3DCuts[0];
    float yzcut = 1.5 * xcut;
    bool useAngle = tjs.Match3DCuts[1] > 0;
    for(unsigned int ipt = 0; ipt < tjs.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = tjs.mallTraj[ipt];
      // length cut
      if(iTjPt.npts < minPts) continue;
      // look for matches using Tjs that have the correct score
      if(iTjPt.score < 0 || iTjPt.score > maxScore) continue;
      unsigned short iplane = DecodeCTP(iTjPt.ctp).Plane;
      // load the CTP and direction so we can find matching angles
      tpi.CTP = iTjPt.ctp;
      tpi.Dir = iTjPt.dir;
      for(unsigned int jpt = ipt + 1; jpt < tjs.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = tjs.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.ctp == iTjPt.ctp) continue;
        // length cut
        if(jTjPt.npts < minPts) continue;
        // ensure they are both showerlike or both not showerlike
        if(jTjPt.showerlike != iTjPt.showerlike) continue;
        if(jTjPt.score < 0 || jTjPt.score > maxScore) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (10 cm)
        if(jTjPt.xlo > iTjPt.xhi + 10) break;
        // ensure the intersection is inside the TPC
        unsigned short jplane = DecodeCTP(jTjPt.ctp).Plane;
        tpj.CTP = jTjPt.ctp;
        tpj.Dir = jTjPt.dir;
        double jyp = 0, jzp = 0;
        tjs.geom->IntersectionPoint(iTjPt.wire, jTjPt.wire, iplane, jplane, (unsigned int)cstat, (unsigned int)tpc, jyp, jzp);
        // get the direction. If this works, IntersectionPoint could be totally replaced with TrajPoint3D
        bool dijOK = false;
        if(useAngle && iTjPt.npts > 5 && jTjPt.npts > 5) dijOK = TrajPoint3D(tjs, tpi, tpj, pos3, dij, prt);
        if(numPlanes == 3) {
          // numPlanes == 3
          for(unsigned int kpt = jpt + 1; kpt < tjs.mallTraj.size(); ++kpt) {
            auto& kTjPt = tjs.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
            // ensure they are all showerlike or all not showerlike
            if(kTjPt.showerlike != iTjPt.showerlike) continue;
            if(kTjPt.score < 0 || kTjPt.score > maxScore) continue;
            if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 10) break;
            unsigned short kplane = DecodeCTP(kTjPt.ctp).Plane;
            tpk.CTP = kTjPt.ctp;
            tpk.Dir = kTjPt.dir;
            double kyp, kzp;
            tjs.geom->IntersectionPoint(iTjPt.wire, kTjPt.wire, iplane, kplane, (unsigned int)cstat, (unsigned int)tpc, kyp, kzp);
            if(std::abs(kyp - jyp) > yzcut || std::abs(kzp - jzp) > yzcut) continue;
            if(useAngle && dijOK && kTjPt.npts > 5 && TrajPoint3D(tjs, tpi, tpk, pos3, dik, prt)) {
              // compare the angles between the
              float dang = dij.Angle(dik);
              if(dang > piOver2) dang = piOver2 - dang;
              if(dang > tjs.Match3DCuts[1]) continue;
            } // useAngle etc
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
            // next check the temp vector
            unsigned short indx = 0;
            for(indx = 0; indx < temp.size(); ++indx) {
              auto& ms = temp[indx];
              if(iTjPt.id != ms.TjIDs[iplane]) continue;
              if(jTjPt.id != ms.TjIDs[jplane]) continue;
              if(kTjPt.id != ms.TjIDs[kplane]) continue;
              ++ms.Count;
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
              ms.Count = 1;
              temp.push_back(ms);
              // give up if there are too many
              if(temp.size() > nAvailable) break;
            } // not found in the list
          } // kpt
          // numPlanes == 3
        } else {
          // 2-plane TPC or 2-plane match in a 3-plane TPC
          if(tjs.NumPlanes == 3) {
            // See if there is a signal at this point. Use tpk but first ensure that the intersection
            // is reasonably OK
            unsigned short kpl = 3 - iplane - jplane;
            float fkwire = tjs.geom->WireCoordinate(jyp, jzp, kpl, tpc, cstat);
            if(fkwire < 0 || fkwire > tjs.MaxPos0[kpl]) continue;
            tpk.CTP = EncodeCTP(cstat, tpc, kpl);
            geo::PlaneID planeID = DecodeCTP(tpi.CTP);
            float xp = 0.5 * (iTjPt.xlo + iTjPt.xhi);
            tpk.Pos[0] = fkwire;
            tpk.Pos[1] = tjs.detprop->ConvertXToTicks(xp, planeID) * tjs.UnitsPerTick;
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
  bool MakeSpt(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, Point3_t& pos, Vector3_t& dir)
  {
    // A variant (replacement) for TrajPoint3D that doesn't use TVector3
    dir = {999};
    pos = {999};
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);
    double ix = tjs.detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iPlnID);
    double jx = tjs.detprop->ConvertTicksToX(jtp.Pos[1] / tjs.UnitsPerTick, jPlnID);
    //    std::cout<<"MSPT: "<<PrintPos(tjs, itp.Pos)<<" X "<<ix<<" "<<PrintPos(tjs, jtp.Pos)<<" "<<jx<<"\n";
    
    // don't continue if the points are wildly far apart in X
    if(std::abs(ix - jx) > 10) return false;
    pos[0] = 0.5 * (ix + jx);
    
    unsigned int iWire = std::nearbyint(itp.Pos[0]);
    if(!tjs.geom->HasWire(geo::WireID(iPlnID.Cryostat, iPlnID.TPC, iPlnID.Plane, iWire))) return false;
    unsigned int jWire = std::nearbyint(jtp.Pos[0]);
    if(!tjs.geom->HasWire(geo::WireID(jPlnID.Cryostat, jPlnID.TPC, jPlnID.Plane, jWire))) return false;
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
    pos[2] = (jcs * (iPos0 - iw0) - ics * (jPos0 - jw0)) / den;
    // and the Y position
    bool useI = std::abs(ics) > std::abs(jcs);
    if(useI) {
      pos[1] = (iPos0 - iw0 - isn * pos[2]) / ics;
    } else {
      pos[1] = (jPos0 - jw0 - jsn * pos[2]) / jcs;
    }
    
    // Now find the direction. Protect against large angles first
    if(jtp.Dir[1] == 0) {
      // Going either in the +X direction or -X direction
      if(jtp.Dir[0] > 0) { dir[0] = 1; } else { dir[0] = -1; }
      dir[1] = 0;
      dir[2] = 0;
      return true;
    } // jtp.Dir[1] == 0
    
    // make a copy of itp and shift it by many wires to avoid precision problems
    TrajPoint itp2 = itp;
    MoveTPToWire(itp2, itp2.Pos[0] + 100);
    // Create a second TVector3 for the shifted point
    Point3_t pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = tjs.detprop->ConvertTicksToX(itp2.Pos[1] / tjs.UnitsPerTick, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = tjs.detprop->ConvertXToTicks(pos2[0], jPlnID) * tjs.UnitsPerTick;
    // Find the wire position (Pos0) in the j plane that this corresponds to
    double jtp2Pos0 = (jtp2Pos1 - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    // Find the Y,Z position using itp2 and jtp2Pos0
    pos2[2] = (jcs * (itp2.Pos[0] - iw0) - ics * (jtp2Pos0 - jw0)) / den;
    if(useI) {
      pos2[1] = (itp2.Pos[0] - iw0 - isn * pos2[2]) / ics;
    } else {
      pos2[1] = (jtp2Pos0 - jw0 - jsn * pos2[2]) / jcs;
    }
    double sep = PosSep(pos, pos2);
    if(sep == 0) return false;
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) dir[ixyz] = (pos2[ixyz] - pos[ixyz]) /sep;
    
    return true;
  } // MakeSpt
  
  ////////////////////////////////////////////////
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2)
  {
    return acos(DotProd(v1, v2));
  } 
  
  ////////////////////////////////////////////////
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2)
  {
    // Finds the direction vector between the two points
    Vector3_t dir;
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
      dir[ixyz] = p2[ixyz] - p1[ixyz];
    }
/*
    std::cout<<"chk "<<std::setprecision(2)<<p1[0]<<" "<<p1[1]<<" "<<p1[2];
    std::cout<<" -> "<<std::setprecision(2)<<p2[0]<<" "<<p2[1]<<" "<<p2[2];
    std::cout<<" dir "<<std::setprecision(3)<<dir[0]<<" "<<dir[1]<<" "<<dir[2];
*/
    if(!SetMag(dir, 1)) { dir[0] = 0; dir[1] = 0; dir[3] = 999; }
/*
    std::cout<<" norm "<<std::setprecision(3)<<dir[0]<<" "<<dir[1]<<" "<<dir[2];
    std::cout<<"\n";
*/
    return dir;
  } // PointDirection

  /////////////////////////////////////////
  bool TrajPoint3D(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, TVector3& pos, TVector3& dir, bool prt)
  {
    // Calculate a 3D position and direction from two trajectory points
    
    dir.SetX(999);
    pos = {0, 0, 0};
    
    if(itp.CTP == jtp.CTP) {
      if(prt) mf::LogVerbatim("TC")<<"TP3D: points "<<itp.CTP<<":"<<PrintPos(tjs, itp.Pos)<<" "<<jtp.CTP<<":"<<PrintPos(tjs, jtp.Pos)<<" are in the same CTP";
      return false;
    }
    
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);
    
    double ix = tjs.detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iPlnID);
    double jx = tjs.detprop->ConvertTicksToX(jtp.Pos[1] / tjs.UnitsPerTick, jPlnID);
    //    std::cout<<"TP3D: "<<PrintPos(tjs, itp.Pos)<<" X "<<ix<<" "<<PrintPos(tjs, jtp.Pos)<<" "<<jx<<"\n";
    
    // don't continue if the points are too far apart in X
    if(std::abs(ix - jx) > 10) {
      if(prt) mf::LogVerbatim("TC")<<"TP3D: points "<<iPlnID.Plane<<":"<<PrintPos(tjs, itp.Pos)<<" "<<jPlnID.Plane<<":"<<PrintPos(tjs, jtp.Pos)<<" too far apart";
      return false;
    }
    pos[0] = 0.5 * (ix + jx);
    
    unsigned int iWire = std::nearbyint(itp.Pos[0]);
    if(!tjs.geom->HasWire(geo::WireID(iPlnID.Cryostat, iPlnID.TPC, iPlnID.Plane, iWire))) {
      if(prt) mf::LogVerbatim("TC")<<"TP3D: No wire at this position "<<iPlnID.Plane<<":"<<PrintPos(tjs, itp.Pos);
      return false;
    }
    unsigned int jWire = std::nearbyint(jtp.Pos[0]);
    if(!tjs.geom->HasWire(geo::WireID(jPlnID.Cryostat, jPlnID.TPC, jPlnID.Plane, jWire))) {
      if(prt) mf::LogVerbatim("TC")<<"TP3D: No wire at this position "<<jPlnID.Plane<<":"<<PrintPos(tjs, jtp.Pos);
      return false;
    }
    
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
    // Find the Z position of the intersection
    pos[2] = (jcs * (itp.Pos[0] - iw0) - ics * (jtp.Pos[0] - jw0)) / den;
    // and the Y position
    if(ics != 0) {
      pos[1] = (itp.Pos[0] - iw0 - isn * pos[2]) / ics;
    } else {
      pos[1] = (jtp.Pos[0] - jw0 - jsn * pos[2]) / jcs;
    }
    
    // Now find the direction. Protect against large angles first
    if(jtp.Dir[1] == 0) {
      // Going either in the +X direction or -X direction
      if(jtp.Dir[0] > 0) {
        dir.SetX(1);
      } else {
        dir.SetX(-1);
      }
      dir.SetY(0);
      dir.SetZ(0);
      return true;
    } // jtp.Dir[1] == 0
    
    // make a copy of itp and shift it by many wires to avoid precision problems
    TrajPoint itp2 = itp;
    MoveTPToWire(itp2, itp2.Pos[0] + 100);
    // Create a second TVector3 for the shifted point
    TVector3 pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = tjs.detprop->ConvertTicksToX(itp2.Pos[1] / tjs.UnitsPerTick, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = tjs.detprop->ConvertXToTicks(pos2[0], jPlnID) * tjs.UnitsPerTick;
    // Find the wire position (Pos0) in the j plane that this corresponds to
    double jtp2Pos0 = (jtp2Pos1 - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    // Find the Y,Z position using itp2 and jtp2Pos0
    pos2[2] = (jcs * (itp2.Pos[0] - iw0) - ics * (jtp2Pos0 - jw0)) / den;
    if(ics != 0) {
      pos2[1] = (itp2.Pos[0] - iw0 - isn * pos2[2]) / ics;
    } else {
      pos2[1] = (jtp2Pos0 - jw0 - jsn * pos2[2]) / jcs;
    }
    dir = pos2 - pos;
    if(dir.Mag() == 0) {
      if(prt) mf::LogVerbatim("TC")<<"TP3D: points "<<iPlnID.Plane<<":"<<PrintPos(tjs, itp.Pos)<<" "<<jPlnID.Plane<<":"<<PrintPos(tjs, jtp.Pos)<<" Magnitude is 0";
      return false;
    }
    dir.SetMag(1);
    
    return true;
    
  } // TrajPoint3D

} // namespace
