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
      // set the end points
/*
      auto& spt = tjs.spts[firstPt];
      pfp.XYZ[0][0] = spt.Pos.X();
      pfp.XYZ[0][1] = spt.Pos.Y();
      pfp.XYZ[0][2] = spt.Pos.Z();
*/
      tjs.pfps.push_back(pfp);
      if(ii > 10) break;
    } // ii (tj)
    
  } // Match3DSpts

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
    // returns a vector of indices into the tjs.spts vector of SpacePoints that
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

} // namespace
