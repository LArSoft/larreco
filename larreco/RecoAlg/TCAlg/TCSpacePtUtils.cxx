#include "larreco/RecoAlg/TCAlg/TCSpacePtUtils.h"
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
  
  //////////////////////////////////////////
  bool Repair(TjStuff& tjs, PFPStruct& pfp, int tjNotInVx, bool prt)
  {
    // The tj tjNotInVx was matched in 3D to other Tjs in the pfp but it wasn't 
    // attached to the vertex that the pfp is attached to. Try to repair it. This
    // function returns false if there is a failure. The NeedsUpdate pfp flag is set
    // true if an update is needed
    if(tjNotInVx <= 0 || tjNotInVx > (int)tjs.allTraj.size()) return false;
    if(pfp.ID == 0) return false;
    if(pfp.TjIDs.empty()) return false;
    auto& missTj = tjs.allTraj[tjNotInVx - 1];
    int mtjid = 0;
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.CTP == missTj.CTP) mtjid = tjid;
    }
    if(mtjid == 0) {
      // This is a strange situation but is easy to repair
      if(prt) mf::LogVerbatim("TC")<<" No other pfp Tj found in the same plane. Simple update needed";
      pfp.TjIDs.push_back(tjNotInVx);
      pfp.NeedsUpdate = true;
      return true;
    }
    // reference to the tj that is 3D-matched
    auto& m3Tj = tjs.allTraj[mtjid - 1];
    // look for a 2D vertex between them
    unsigned short vx2id = 0;
    for(unsigned short e1 = 0; e1 < 2; ++e1) {
      if(m3Tj.VtxID[e1] == 0) continue;
      for(unsigned short e2 = 0; e2 < 2; ++e2) {
        if(missTj.VtxID[e2] == m3Tj.VtxID[e1]) vx2id = missTj.VtxID[e2];
      } // e2
    } // e1
    if(vx2id == 0) {
      std::cout<<"Repair: No common vertex between "<<missTj.ID<<" and "<<m3Tj.ID<<". Deal with this case. \n";
      return false;
    }
    // separation distance (cm) for kink detection
    double sep = 1;
    auto kinkPts = FindKinks(tjs, pfp, sep, prt);
    if(prt) mf::LogVerbatim("TC")<<" number of kinks found "<<kinkPts.size();
/*
    // The vertex exists because there is a significant kink in one plane
    // but not a significant kink in the other planes. We need to split the
    // Tjs in the other planes and re-build
    if(prt) mf::LogVerbatim("TC")<<" Two Tjs needs to be split at the X position of vx2id "<<vx2id;
    auto& vx2 = tjs.vtx[vx2id - 1];
    geo::PlaneID planeID = DecodeCTP(vx2.CTP);
    float splitX = tjs.detprop->ConvertTicksToX(vx2.Pos[1] / tjs.UnitsPerTick, planeID);
    // create a 3D vertex
    Vtx3Store vx3;
    vx3.TPCID = pfp.TPCID;
    vx3.Vx2ID[planeID.Plane] = vx2id;
    for(auto tjid : pfp.TjIDs) {
      if(tjid == mtjid) continue;
      if(!SplitTraj(tjs, tjid - 1, splitX, true, prt)) continue;
      // If successful, SplitTraj will have created a new 2D vertex and put it at the end of tjs.vtx
      geo::PlaneID planeID = DecodeCTP(tjs.allTraj[tjid - 1].CTP);
      vx3.Vx2ID[planeID.Plane] = tjs.vtx.size();
    } // tjid
    std::cout<<"Repair: Find 3D vertex position and re-build everything\n";
    tjs.Needs3DRebuild
*/
    return true;
  } // Repair
/*
  //////////////////////////////////////////
  void Match3DSpts(TjStuff& tjs, const art::Event& evt, const geo::TPCID& tpcid, 
                   const art::InputTag& fSpacePointModuleLabel)
  {
    // Match Tjs in 3D using SpacePoint <-> Hit associations.
    if(!tjs.UseAlg[kHitsOrdered]) {
      std::cout<<"Match3DSpts: Write some code to deal with un-ordered hits\n";
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
      PFPStruct pfp = CreatePFP(tjs, tpcid);
      pfp.TjIDs.clear();
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        if(tjInPln[plane].empty()) continue;
        if(tjInPln[plane][0] > 0) pfp.TjIDs.push_back(tjInPln[plane][0]);
      } // plane
      // set the end points using the local version that uses SpacePoints instead
      // of tjs.malltraj
      if(!SetPFPEndPoints(tjs, pfp, sptLists, tj.ID, prt)) {
        std::cout<<"SetPFPEndPoints failed";
      }
      if(!StorePFP(tjs, pfp)) continue;
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
*/
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
  void MakePFPTp3s(TjStuff& tjs, PFPStruct& pfp, bool anyTj)
  {
    // Creates a vector of TrajCluster space points using 3D matches in mallTraj. If anyTj is set true
    // any Trajectory point in the 3rd plane that is consistent will be added as well, otherwise
    // only those Tjs that are in the pfp.TjIDs list are considered. Note that the Tp3 directions may
    // be 
    
    pfp.Tp3s.clear();
    if(pfp.ID == 0) return;
    if(pfp.TjIDs.empty()) return;
    if(tjs.mallTraj.empty()) return;
    double xcut = tjs.Match3DCuts[0];
    double yzcut = 1.5 * xcut;
    constexpr double twopi = 2 * M_PI;
    
//    MCParticleListUtils mcplu{tjs};
    
    for(unsigned int ipt = 0; ipt < tjs.mallTraj.size() - 1; ++ipt) {
      auto& iTj2Pt = tjs.mallTraj[ipt];
      if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), iTj2Pt.id) == pfp.TjIDs.end()) continue;
      auto& itp = tjs.allTraj[iTj2Pt.id - 1].Pts[iTj2Pt.ipt];
      // temp 
/*
      unsigned int imcpi = mcplu.GetMCPartListIndex(itp);
      int TMeV = 0;
      if(imcpi != UINT_MAX) {
        auto& mcp = tjs.MCPartList[imcpi];
        TMeV = 1000 * (mcp->E() - mcp->Mass());
      }
*/
      for(unsigned int jpt = ipt + 1; jpt < tjs.mallTraj.size(); ++jpt) {
        auto& jTj2Pt = tjs.mallTraj[jpt];
        // break out if we are well past the X matching region
        if(jTj2Pt.xlo > iTj2Pt.xhi + 5) break;
        if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), jTj2Pt.id) == pfp.TjIDs.end()) continue;
        // ensure that the planes are different
        if(jTj2Pt.ctp == iTj2Pt.ctp) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTj2Pt.xlo > iTj2Pt.xhi) continue;
        auto& jtp = tjs.allTraj[jTj2Pt.id - 1].Pts[jTj2Pt.ipt];
//        unsigned int jmcpi = mcplu.GetMCPartListIndex(jtp);
        TrajPoint3 tp3;
        if(!MakeTp3(tjs, itp, jtp, tp3)) continue;
        tp3.Tj2Pts.resize(2);
        tp3.Tj2Pts[0] = iTj2Pt;
        tp3.Tj2Pts[1] = jTj2Pt;
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
//            unsigned int kmcpi = mcplu.GetMCPartListIndex(ktp);
//            bool ijkMatch = (imcpi != UINT_MAX && jmcpi == imcpi && kmcpi == imcpi);
            TrajPoint3 iktp3;
            if(!MakeTp3(tjs, itp, ktp, iktp3)) continue;
            if(std::abs(iktp3.Pos[1] - tp3.Pos[1]) > yzcut) continue;
            if(std::abs(iktp3.Pos[2] - tp3.Pos[2]) > yzcut) continue;
            if(tjs.Match3DCuts[1] > 0) {
              double dang = DeltaAngle(tp3.Dir, iktp3.Dir);
              while(dang >  M_PI) dang -= twopi;
              while(dang < -M_PI) dang += twopi;
              if(dang >  M_PI/2) dang = M_PI - dang;
              if(dang < -M_PI/2) dang = M_PI + dang;
              float mcsmom = tjs.allTraj[iTj2Pt.id - 1].MCSMom + tjs.allTraj[jTj2Pt.id - 1].MCSMom + tjs.allTraj[kTj2Pt.id - 1].MCSMom;
              mcsmom /= 3;
//              mf::LogVerbatim("TC")<<"dang "<<std::fixed<<std::setprecision(3)<<dang<<std::setprecision(2)<<" "<<(iktp3.Pos[1] - tp3.Pos[1])<<" "<<(iktp3.Pos[2] - tp3.Pos[2])<<" "<<(int)mcsmom<<" "<<ijkMatch<<" "<<TMeV;
              if(mcsmom > 150 && dang > tjs.Match3DCuts[1]) continue;
            }
            // update the position and direction
            for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
              tp3.Pos[ixyz] += iktp3.Pos[ixyz]; tp3.Pos[ixyz] /= 2;
              tp3.Dir[ixyz] += iktp3.Dir[ixyz]; tp3.Dir[ixyz] /= 2;
            }
            tp3.Tj2Pts.push_back(kTj2Pt);
            // only allow one addition
            break;
          } // kpt
        } // 3 planes
        FilldEdx(tjs, tp3);
        pfp.Tp3s.push_back(tp3);
        break;
      } // jpt
    } // ipt
    
  } // MakePFPTp3s

  
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
      auto& tp3 = pfp.Tp3s[ipt];
      auto& pos = tp3.Pos;
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
    } // tp3
    
    // check for a fully contained track
    bool insideX = minXYZ[0] > tjs.XLo && maxXYZ[0] < tjs.XHi;
    bool insideY = minXYZ[1] > tjs.YLo && maxXYZ[1] < tjs.YHi;
    bool insideZ = minXYZ[2] > tjs.ZLo && maxXYZ[2] < tjs.ZHi;
    bool fullyContained = insideX && insideY && insideZ;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<" Xlo "<<std::setw(7)<<tjs.XLo<<" minX "<<std::setw(7)<<minXYZ[0]<<" maxX "<<std::setw(7)<<maxXYZ[0]<<" XHi "<<std::setw(7)<<tjs.XHi<<" insideX? "<<insideX<<"\n";
      myprt<<" Ylo "<<std::setw(7)<<tjs.YLo<<" minY "<<std::setw(7)<<minXYZ[1]<<" maxY "<<std::setw(7)<<maxXYZ[1]<<" YHi "<<std::setw(7)<<tjs.YHi<<" insideY? "<<insideY<<"\n";
      myprt<<" Zlo "<<std::setw(7)<<tjs.ZLo<<" minZ "<<std::setw(7)<<minXYZ[2]<<" maxZ "<<std::setw(7)<<maxXYZ[2]<<" ZHi "<<std::setw(7)<<tjs.ZHi<<" insideZ? "<<insideZ;
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
//    mf::LogVerbatim("TC")<<"SBDFS: StartPos "<<std::fixed<<std::setprecision(1)<<startPos[0]<<" "<<startPos[1]<<" "<<startPos[2];

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
    std::vector<TrajPoint3> temp;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) temp.push_back(pfp.Tp3s[sortVec[ii].index]);
    pfp.Tp3s = temp;
    
  } // SortByDistanceFromStart
  
  /////////////////////////////////////////
  bool CheckTp3Validity(TjStuff& tjs, PFPStruct& pfp, Vector3_t generalDirection, bool prt)
  {
    // Checks and corrects the ordering of Tj points in the pfp space points. Checks for
    // consistency between the tjs listed in pfp.TjIDs and those appearing in the space points.
    // Sets the IsValid flag false for space points in the vector that are inconsistent
    // with their neighbors and with the general direction
    
    if(prt) mf::LogVerbatim("TC")<<"CheckTp3Validity: pfp "<<pfp.ID;

    // find the average separation between adjacent points
    float sum = 0;
    float sum2 = 0;
    float cnt = 0;
    std::vector<float> seps(pfp.Tp3s.size());
    for(unsigned short ipt = 1; ipt < pfp.Tp3s.size(); ++ipt) {
      float sep = PosSep(pfp.Tp3s[ipt].Pos, pfp.Tp3s[ipt - 1].Pos);
      seps[ipt] = sep;
      // ignore large separations (> ~1.5 cm)
      if(sep > 1.5) continue;
      sum += sep;
      sum2 += sep * sep;
      ++cnt;
    } // ii
    // return if something is seriously wrong
    if(cnt < 0.5 * pfp.Tp3s.size()) {
      if(prt) mf::LogVerbatim("TC")<<" Tp3 separations too large. cnt "<<cnt<<" minimum required "<<0.5 * pfp.Tp3s.size();
      return false;
    }
    float aveSep = sum / cnt;
    float arg = sum2 - cnt * aveSep * aveSep;
    float maxSep = 2 * aveSep;
    if(arg > 0) {
      float rms = sqrt(arg / (cnt - 1));
      maxSep = aveSep + 3 * rms;
    }
    // decide whether to use angle in the decision
    float mcsmom = 0;
    for(auto tjid : pfp.TjIDs) mcsmom += tjs.allTraj[tjid - 1].MCSMom;
    mcsmom /= (float)pfp.TjIDs.size();
    bool useAngle = tjs.Match3DCuts[1] > 0 && mcsmom > 150;
    if(prt) mf::LogVerbatim("TC")<<" ave Tj mcsmom "<<mcsmom<<" use Angle (< 0.3) to check validity? "<<useAngle;
    unsigned short nValid = 0;
    for(unsigned short ipt = 0; ipt < pfp.Tp3s.size(); ++ipt) {
      // Set large sep points not valid
      if(seps[ipt] > maxSep) pfp.Tp3s[ipt].IsValid = false;
      // Set not valid if there is a large change in angle
      float dang = 0;
      if(useAngle) {
        dang = DeltaAngle(pfp.Tp3s[ipt].Dir, generalDirection);
        if(dang > 0.3) pfp.Tp3s[ipt].IsValid = false;
      }
//      if(prt) mf::LogVerbatim("TC")<<ipt<<" sep "<<seps[ipt]<<" dang "<<dang<<" isValid? "<<pfp.Tp3s[ipt].IsValid;
      // update the general direction vector
      if(pfp.Tp3s[ipt].IsValid) generalDirection = pfp.Tp3s[ipt].Dir;
      if(pfp.Tp3s[ipt].IsValid) ++nValid;
    } // ipt
    // require most of the point valid
    if(nValid < 0.5 * pfp.Tp3s.size()) {
      if(prt) mf::LogVerbatim("TC")<<" Not enough valid points "<<nValid<<" minimum required "<<0.5 * pfp.Tp3s.size();
      return false;
    }
    if(prt) mf::LogVerbatim("TC")<<" Found "<<nValid<<" Tp3s. Looks good. ";
    return true;
  } // CheckTp3Validity
  
  /////////////////////////////////////////
  bool FitTp3(TjStuff& tjs, std::vector<TrajPoint3> tp3s, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint3& outTp3)
  {
    // fits a section of the vector of Tp3s from point fromPt to toPt and returns the
    // result in outTp3. This function returns false if there
    // was a failure
    outTp3.IsValid = false;
    if(originPt > tp3s.size() - 1) return false;
    if(!(fitDir == 1 || fitDir == -1)) return false;
    unsigned short toPt = originPt + fitDir * npts;
    if(toPt > tp3s.size() - 1) return false;
    // This doesn't really do a fit (for now) but does a simple average of the directions defined by pairs.
    // get a loop lower and upper limit
    unsigned short fromPt = originPt;
    if(toPt < fromPt) std::swap(toPt, fromPt);
    outTp3.Dir = {0, 0, 0};
    double wgtsum = 0;
    std::cout<<"Start summing\n";
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      if(!tp3s[ipt].IsValid) continue;
      auto& itpt0 = tp3s[ipt].Tj2Pts[0];
      float x0 = 0.5 * (itpt0.xlo + itpt0.xhi);
      auto& itpt1 = tp3s[ipt].Tj2Pts[1];
      float x1 = 0.5 * (itpt1.xlo + itpt1.xhi);
      float idx = std::abs(x1 - x0);
      if(idx < 0.01) idx = 0.01;
      for(unsigned short jpt = ipt + 2; jpt <= toPt; ++jpt) {
        if(!tp3s[jpt].IsValid) continue;
        auto& jtpt0 = tp3s[jpt].Tj2Pts[0];
        float x0 = 0.5 * (jtpt0.xlo + jtpt0.xhi);
        auto& jtpt1 = tp3s[jpt].Tj2Pts[1];
        float x1 = 0.5 * (jtpt1.xlo + jtpt1.xhi);
        float jdx = std::abs(x1 - x0);
        if(jdx < 0.01) jdx = 0.01;
        // Find the direction vector between these points
        double sep2 = PosSep2(tp3s[ipt].Pos, tp3s[jpt].Pos);
        // ignore a separation less than 0.1 cm
        if(sep2 < 0.01) continue;
        auto dir = PointDirection(tp3s[ipt].Pos, tp3s[jpt].Pos);
        std::cout<<"ipt "<<ipt<<" "<<jpt<<" dir "<<std::fixed<<std::setprecision(3)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<" sep "<<sqrt(sep2);
        std::cout<<" idx "<<idx<<" "<<jdx;
        std::cout<<"\n";
        // weight by the separation^2
        double wgt = sep2 / (idx * jdx);
        for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) outTp3.Dir[ixyz] += dir[ixyz] * wgt;
        wgtsum += wgt;
      } // jpt
    } // ipt
    if(wgtsum == 0) return false;
    outTp3.Pos = tp3s[originPt].Pos;
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) outTp3.Dir[ixyz] /= wgtsum;
    std::cout<<"result "<<std::fixed<<std::setprecision(3)<<outTp3.Dir[0]<<" "<<outTp3.Dir[1]<<" "<<outTp3.Dir[2]<<"\n";
    outTp3.IsValid = SetMag(outTp3.Dir, 1);
    return outTp3.IsValid;

  } // FitTP3

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
    unsigned short minPts = 3;
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
        // ensure they are both showerlike or both not showerlike
        if(jTjPt.showerlike != iTjPt.showerlike) continue;
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
            // ensure they are all showerlike or all not showerlike
            if(kTjPt.showerlike != iTjPt.showerlike) continue;
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
    //    std::cout<<"MSPT: "<<PrintPos(tjs, itp.Pos)<<" X "<<ix<<" "<<PrintPos(tjs, jtp.Pos)<<" "<<jx<<"\n";
    
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
    if(!SetMag(dir, 1)) { dir[0] = 0; dir[1] = 0; dir[3] = 999; }
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
//      std::cout<<cnt<<" "<<PrintPos(tjs, tp)<<" dedx "<<dedx<<" dx "<<dx<<"\n";
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
          if(prt) mf::LogVerbatim("TC")<<" findKinks: kink angle "<<kang<<" at point "<<ipt;
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
  void PrintTp3(std::string fcnLabel, const TjStuff& tjs, const TrajPoint3& tp3)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<fcnLabel<<" Tp3";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(6)<<tp3.Pos[0]<<std::setw(6)<<tp3.Pos[1]<<std::setw(6)<<tp3.Pos[2];
    myprt<<std::fixed<<std::setprecision(3);
    myprt<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
    myprt<<" dEdx "<<std::setw(4)<<std::setprecision(1)<<tp3.dEdx;
    myprt<<" Err "<<std::setw(4)<<std::setprecision(1)<<tp3.dEdxErr;
    myprt<<" IsValid? "<<tp3.IsValid;
    myprt<<" tj_ipt";
    for(auto tj2pt : tp3.Tj2Pts) {
      auto& tj = tjs.allTraj[tj2pt.id - 1];
      auto& tp = tj.Pts[tj2pt.ipt];
      myprt<<" "<<tj.ID<<"_"<<PrintPos(tjs, tp);
    } // tj2pt
  } // PrintTp3

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
  } // MatchVecIndex
  
  ////////////////////////////////////////////////
  bool StorePFP(TjStuff& tjs, PFPStruct& pfp)
  {
    // stores the PFParticle in TJStuff
    if(pfp.ID < tjs.pfps.size()) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    if(!neutrinoPFP && pfp.TjIDs.empty()) return false;
    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)tjs.pfps.size() + 1) {
      std::cout<<"StorePFP ID is wrong. Fixing it\n";
      pfp.ID = tjs.pfps.size() + 1;
    }
    // check the Tjs and set the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kMat3D]) {
        int matchedTo = 0;
        for(auto& opfp : tjs.pfps) {
          if(std::find(opfp.TjIDs.begin(), opfp.TjIDs.end(), tjid) != opfp.TjIDs.end()) {
            matchedTo = opfp.ID;
            break;
          }
        } // old pfp
        std::cout<<"StorePFP "<<pfp.ID<<" Tj "<<tj.ID<<" is already 3D matched to pfp "<<matchedTo<<"\n";
        return false;
      }
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
            std::cout<<"StorePFP "<<pfp.ID<<" Violating the SetDir flag for Tj "<<tj.ID<<"\n";
            tj.AlgMod[kSetDir] = false;
          }
          ReverseTraj(tjs, tj);
        } // lastIpt[ii] > firstIpt[ii]
        tj.AlgMod[kSetDir] = true;
      } // ii
    } // Tp3s exist    
    
    tjs.pfps.push_back(pfp);
    return true;
  } // StorePFP

} // namespace
