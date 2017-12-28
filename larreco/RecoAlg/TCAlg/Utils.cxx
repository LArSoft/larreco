#include "larreco/RecoAlg/TCAlg/Utils.h"

struct SortEntry{
  unsigned int index;
  float val;
};

bool valDecreasing (SortEntry c1, SortEntry c2) { return (c1.val > c2.val);}
bool valIncreasing (SortEntry c1, SortEntry c2) { return (c1.val < c2.val);}

namespace tca {
  
  /////////////////////////////////////////
  void DefineTjParents(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
/*
    This function sets the ParentID of Tjs in this tpcid to create a hierarchy. The highest Score
    3D vertex in a chain of Tjs and vertices is declared the primary vertex; vx3.Primary = true. Tjs directly attached
    to that vertex are declared Primary trajectories with ParentID = 0. All other Tjs in the chain have ParentID
    set to the next upstream Tj to which it is attached by a vertex. In the graphical description below, V1 and V4 are 
    2D vertices that are matched to a high-score 3D vertex. The V1 Score is greater than the V2 Score and V3 Score. 
    V1 and V4 are declared to be primary vertices. T1, T2, T6 and T7 are declared to be primary Tjs.

      V1 - T1 - V2 - T3          V4 - T6         / T8
         \                          \           /
           T2 - V3 - T4               T7
                   \
                     T5
 
    This is represented as follows. The NeutrinoPrimaryTjID is defined by a function.
     Tj   ParentID   NeutrinoPrimaryTjID
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
    // sort vertice by decreasing score
    std::vector<int> temp;
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      // clear the Primary flag while we are here
      vx3.Primary = false;
      temp.push_back(vx3.ID);
    } // vx3
    if(temp.empty()) return;
    
    // Make a master list of all Tjs that are attached to these vertices
    std::vector<int> masterlist;
    for(auto vx3id : temp) {
      auto& vx3 = tjs.vtx3[vx3id - 1];
      float score;
      auto tjlist = GetVtxTjIDs(tjs, vx3, score);
      for(auto tjid : tjlist) {
        // Temp? Check for an existing parentID
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.ParentID != -1) {
          std::cout<<"**** Tj "<<tj.ID<<" Existing parent "<<tj.ParentID<<" PDGCode "<<tj.PDGCode<<". with a vertex... \n";
          tj.ParentID = -1;
        }
        if(std::find(masterlist.begin(), masterlist.end(), tjid) == masterlist.end()) masterlist.push_back(tjid);
      } // tjid
    } // vxid
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"DTP: masterlist Tjs";
      for(auto tjid : masterlist) myprt<<" "<<tjid;
    }
    
    // Do the sort
    std::vector<SortEntry> sortVec(temp.size());
    for(unsigned short indx = 0; indx < temp.size(); ++indx) {
      auto& vx3 = tjs.vtx3[temp[indx] - 1];
      sortVec[indx].index = indx;
      sortVec[indx].val = vx3.Score;
    } // indx
    if(sortVec.size() > 1) std::sort(sortVec.begin(), sortVec.end(), valDecreasing);
    // put them into order
    auto vlist = temp;
    for(unsigned short indx = 0; indx < temp.size(); ++indx) vlist[indx] = temp[sortVec[indx].index];
    
    // make a neutrino PFParticle to associate with the highest score vertex
    auto neutrinoPFP = CreatePFPStruct(tjs, tpcid);
    auto& vx3 = tjs.vtx3[vlist[0] - 1];
    // call it the neutrino vertex
    vx3.Neutrino = true;
    // put the vertex at the end of the neutrino
    neutrinoPFP.XYZ[1][0] = vx3.X;
    neutrinoPFP.XYZ[1][1] = vx3.Y;
    neutrinoPFP.XYZ[1][2] = vx3.Z;
    neutrinoPFP.XYZ[0] = neutrinoPFP.XYZ[1];
    neutrinoPFP.Dir[1][2] = 1;
    neutrinoPFP.Dir[0][2] = 1;
    // This may be set to 12 later on if a primary shower is reconstructed 
    neutrinoPFP.PDGCode = 14;
    neutrinoPFP.Vx3ID[1] = vx3.ID;
    neutrinoPFP.Vx3ID[0] = vx3.ID;
    // the rest of this will be defined later
    tjs.pfps.push_back(neutrinoPFP);
    
    // a temp vector to ensure that we only consider a vertex once
    std::vector<bool> lookedAt3(tjs.vtx3.size() + 1, false);
    std::vector<bool> lookedAt2(tjs.vtx.size() + 1, false);
    // vector of parent-daughter pairs
    std::vector<std::pair<int, int>> pardtr;
    // Start with the highest score vertex 
    for(unsigned short indx = 0; indx < vlist.size(); ++indx) {
      auto& vx3 = tjs.vtx3[vlist[indx] - 1];
      if(lookedAt3[vx3.ID]) continue;
      vx3.Primary = true;
      lookedAt3[vx3.ID] = true;
      // make a list of Tjs attached to this vertex
      float score;
      auto primTjList = GetVtxTjIDs(tjs, vx3, score);
      if(primTjList.empty()) continue;
      pardtr.clear();
      for(auto primTjID : primTjList) {
        auto& primTj = tjs.allTraj[primTjID - 1];
        // This isn't a primary tj if the parent ID isn't -1
        if(primTj.ParentID != -1) continue;
        if(prt) mf::LogVerbatim("TC")<<"Vx3 "<<vx3.ID<<" Primary tj "<<primTj.ID;
        // declare this a primary tj
        primTj.ParentID = 0;
        // look for daughter tjs = those that are attached to a 2D vertex
        // at the other end
        for(unsigned short end = 0; end < 2; ++end) {
          if(primTj.VtxID[end] == 0) continue;
          auto& vx2 = tjs.vtx[primTj.VtxID[end] - 1];
          if(vx2.Vx3ID == vx3.ID) continue;
          // found a 2D vertex. Check for daughters
          auto dtrList = GetVtxTjIDs(tjs, vx2);
          for(auto dtrID : dtrList) {
            // ignore the primary tj
            if(dtrID == primTjID) continue;
            auto& dtj = tjs.allTraj[dtrID - 1];
            if(dtj.ParentID != -1) {
              std::cout<<"DTP Error: dtr "<<dtrID<<" already has a parent "<<dtj.ParentID<<". Can't make it daughter of "<<primTjID<<"\n";
              continue;
            }
            pardtr.push_back(std::make_pair(primTjID, dtrID));
            if(prt) mf::LogVerbatim("TC")<<"  primTj "<<primTjID<<" dtrID "<<dtrID;
          } // tjid
        } // end
        // Ensure that end 0 of the trajectory is attached to the primary vertex
        for(unsigned short end = 0; end < 2; ++end) {
          if(primTj.VtxID[end] == 0) continue;
          auto& vx2 = tjs.vtx[primTj.VtxID[end] - 1];
          if(vx2.Vx3ID == vx3.ID && end != 0) ReverseTraj(tjs, primTj);
        } // end
        primTj.AlgMod[kSetDir] = true;
      } // tjid
      if(pardtr.empty()) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" par_dtr";
        for(auto pdtr : pardtr) myprt<<" "<<pdtr.first<<"_"<<pdtr.second;
      }
      // iterate through the parent - daughter stack, removing the last pair when a 
      // ParentID is updated and adding pairs for new daughters
      for(unsigned short nit = 0; nit < 100; ++nit) {
        auto lastPair = pardtr[pardtr.size() - 1];
        auto& dtj = tjs.allTraj[lastPair.second - 1];
        dtj.ParentID = lastPair.first;
        // reverse the daughter trajectory if necessary so that end 0 is closest to the parent
        float doca = 100;
        unsigned short dpt = 0, ppt = 0;
        auto& ptj = tjs.allTraj[lastPair.first - 1];
        // find the point on the daughter tj that is closest to the parent
        TrajTrajDOCA(tjs, dtj, ptj, dpt, ppt, doca);
        // reverse the daughter if the closest point is near end 1 of the daughter
        unsigned short midPt = (dtj.EndPt[0] + dtj.EndPt[1]) / 2;
        if(dpt > midPt && !dtj.AlgMod[kSetDir]) ReverseTraj(tjs, dtj);
        if(prt) mf::LogVerbatim("TC")<<"Set parent "<<ptj.ID<<" dtr "<<dtj.ID;
        dtj.AlgMod[kSetDir] = true;
        // remove that entry
        pardtr.pop_back();
        // Add entries for new daughters
        for(unsigned short end = 0; end < 2; ++end) {
          if(dtj.VtxID[end] == 0) continue;
          auto& vx2 = tjs.vtx[dtj.VtxID[end] - 1];
          if(lookedAt2[vx2.ID]) continue;
          lookedAt2[vx2.ID] = true;
          auto tjlist = GetVtxTjIDs(tjs, vx2);
          for(auto tjid : tjlist) {
            if(tjid == dtj.ID || tjid == ptj.ID) continue;
            pardtr.push_back(std::make_pair(dtj.ID, tjid));
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<" add par_dtr";
              for(auto pdtr : pardtr) myprt<<" "<<pdtr.first<<"_"<<pdtr.second;
            }
          }
        } // end
        if(pardtr.empty()) break;
      } // nit
    } // indx
    if(!pardtr.empty()) {
      std::cout<<"DefineTjParents: pardtr isn't empty...\n";
    }
    // check the master list
    for(auto tjid : masterlist) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.ParentID < 0) {
        std::cout<<"Tj "<<tj.ID<<" is in the master list but doesn't have a Parent\n";
        tj.ParentID = tj.ID;
      }
    } // tjid

  } // DefineTjParents
  
  /////////////////////////////////////////
  void DefinePFParticleRelationships(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
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
      if(neutrinoPFPID == 0 && (pfp.PDGCode == 12 || pfp.PDGCode == 14)) neutrinoPFPID = pfp.ID;
      if(pfp.Vx3ID[0] > 0) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      std::cout<<"DPFPR: Making a bogus PFP vertex\n";
      vx3.ID = tjs.vtx3.size() + 1;
      vx3.Primary = true;
      // TODO: we need to have PFP track position errors defined 
      unsigned short mergeToVx3ID = IsCloseToVertex(tjs, vx3);
      if(mergeToVx3ID > 0) {
        if(prt) mf::LogVerbatim("TC")<<"Merge PFP vertex "<<vx3.ID<<" with existing vtx "<<mergeToVx3ID;
        if(!AttachPFPToVertex(tjs, pfp, 0, mergeToVx3ID, prt)) {
          if(prt) mf::LogVerbatim("TC")<<" Failed to attach pfp "<<pfp.ID<<". Make new vertex \n";
          mergeToVx3ID = 0;
        }
      } // mergeMe > 0
      if(mergeToVx3ID == 0) {
        // Add the new vertex and attach the PFP to it
        tjs.vtx3.push_back(vx3);
        if(!AttachPFPToVertex(tjs, pfp, 0, vx3.ID, prt)) {
          if(prt) mf::LogVerbatim("TC")<<"Merge PFP vertex "<<vx3.ID<<" with new vtx "<<mergeToVx3ID;
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
  } // DefinePFParticleRelationships
  
  /////////////////////////////////////////
  float MaxChargeAsymmetry(TjStuff& tjs, std::vector<int>& tjIDs)
  {
    // calculates the maximum charge asymmetry of the supplied list of Tjs
    if(tjIDs.size() < 2) return 1;
    std::vector<float> tjchg(tjIDs.size(), 0);
    // the average charge has already been calculated. Multiply this by the
    // total number of points to get an estimate that includes dead wires
    for(unsigned short indx = 0; indx < tjIDs.size(); ++indx) {
      int tjid = tjIDs[indx];
      if(tjid <= 0 || tjid > (int)tjs.allTraj.size()) return 1;
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AveChg <= 0) UpdateAveChg(tjs, tj);
      tjchg[indx] = (tj.EndPt[1] - tj.EndPt[0]) * tj.AveChg;
    } // indx
    float maxAsym = 0;
    for(unsigned short indx = 0; indx < tjIDs.size() - 1; ++indx) {
      for(unsigned short jndx = indx + 1; jndx < tjIDs.size(); ++jndx) {
        float asym = std::abs(tjchg[indx] - tjchg[jndx]) / (tjchg[indx] + tjchg[jndx]);
        if(asym > maxAsym) maxAsym = asym;
      } // jtj
    } // itj
    return maxAsym;
  } // MaxChargeAsymmetry
  
  /////////////////////////////////////////
  int PDGCodeVote(TjStuff& tjs, std::vector<int>& tjIDs, bool prt)
  {
    // Returns the most likely PDGCode for the set of Tjs provided
    // The PDG codes are:
    // 0 = your basic track-like trajectory
    // 11 = Tagged delta-ray
    // 13 = Tagged muon
    // 211 = pion-like. There exists a Bragg peak at an end with a vertex
    // 2212 = proton-like. There exists a Bragg peak at an end without a vertex
    std::array<int, 5> codeList = {0, 11, 13, 211, 2212};
    unsigned short codeIndex = 0;
    if(tjIDs.empty()) return codeList[codeIndex];
    
    std::array<unsigned short, 5> cnts = {0};
    // Count Bragg peaks. This assumes that the Tjs are in order...
    std::array<unsigned short, 2> stopCnt {0};
    float maxLen = 0;
    for(auto tjid : tjIDs) {
      if(tjid <= 0 || tjid > (int)tjs.allTraj.size()) continue;
      auto& tj = tjs.allTraj[tjid - 1];
      for(unsigned short ii = 0; ii < 5; ++ii) if(tj.PDGCode == codeList[ii]) ++cnts[ii];
      for(unsigned short end = 0; end < 2; ++end) if(tj.StopFlag[kBragg][end]) ++stopCnt[end];
      float len = TrajLength(tj);
      if(len > maxLen) maxLen = len;
    } // tjid
    unsigned maxCnt = 0;
    // ignore the first PDG code in the list which is the default
    for(unsigned short ii = 1; ii < 5; ++ii) {
      if(cnts[ii] > maxCnt) {
        maxCnt = cnts[ii];
        codeIndex = ii;
      }
    } // ii
    // check for an inconsistent code
    bool confused = false;
    for(unsigned short ii = 1; ii < 5; ++ii) {
      if(ii == codeIndex) continue;
      if(cnts[ii] == 0) continue;
      // make it pion-like if one of the codes is pion-like
      if(ii == 3) {
        codeIndex = 3;
      } else {
        confused = true;
      }
    } // ii
    if(confused) {
      // Check for a muon called it a proton
      if(cnts[4] > 0 && stopCnt[2] > 0 && NumDeltaRays(tjs, tjIDs) == 0) {
        codeIndex = 4;
        confused = false;
      }
    } // confused
    if(confused) {
      codeIndex = 0;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"PDGCodeVote: mixed vote on the PDGCode: Tj_PDGCode";
        for(auto tjid : tjIDs) {
          if(tjid <= 0 || tjid > (int)tjs.allTraj.size()) continue;
          auto& tj = tjs.allTraj[tjid - 1];
          myprt<<" "<<tj.ID<<"_"<<tj.PDGCode<<"_"<<tj.StopFlag[kBragg][1];
        } // tjid
      }
    } // confused
    return codeList[codeIndex];
  } // PDGCodeVote
  
  /////////////////////////////////////////
  unsigned short NumDeltaRays(const TjStuff& tjs, const Trajectory& tj)
  {
    // returns the number of delta rays that have this tj as a parent
    unsigned short cnt = 0;
    for(auto& dtj : tjs.allTraj) {
      if(dtj.AlgMod[kKilled]) continue;
      if(!dtj.AlgMod[kDeltaRay]) continue;
      if(dtj.ParentID == tj.ID) ++cnt;
    } // tj
    return cnt;
  } // NumDeltaRays
  
  /////////////////////////////////////////
  unsigned short NumDeltaRays(const TjStuff& tjs, std::vector<int>& tjIDs)
  {
    // Count the number of delta-rays that have a Tj in the list of TjIDs as a parent.
    if(tjIDs.empty()) return 0;
    if(tjIDs[0] <= 0 || tjIDs[0] > (int)tjs.allTraj.size()) return 0;
    unsigned short cnt = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(!tj.AlgMod[kDeltaRay]) continue;
      if(std::find(tjIDs.begin(), tjIDs.end(), tj.ParentID) != tjIDs.end()) ++cnt;
    } // tj
    return cnt;
  } // NumDeltaRays

  /////////////////////////////////////////
  int NeutrinoPrimaryTjID(const TjStuff& tjs, const Trajectory& tj)
  {
    // Returns the ID of the grandparent of this tj that is a primary tj that is attached
    // to the neutrino vertex. 0 is returned if this condition is not met.
    if(tj.AlgMod[kKilled]) return -1;
    if(tj.ParentID <= 0) return -1;
    int primID = PrimaryID(tjs, tj);
    if(primID <= 0 || primID > (int)tjs.allTraj.size()) return -1;

    // We have the ID of the primary tj. Now see if it is attached to the neutrino vertex
    auto& ptj = tjs.allTraj[primID - 1];
    for(unsigned short end = 0; end < 2; ++end) {
      if(ptj.VtxID[end] == 0) continue;
      auto& vx2 = tjs.vtx[ptj.VtxID[end] - 1];
      if(vx2.Vx3ID == 0) continue;
      auto& vx3 = tjs.vtx3[vx2.Vx3ID - 1];
      if(vx3.Neutrino) return primID;
    } // end
    return -1;
  } // NeutrinoPrimaryTjID
  
  /////////////////////////////////////////
  int PrimaryID(const TjStuff& tjs, const Trajectory& tj)
  {
    // Returns the ID of the grandparent trajectory of this trajectory that is a primary
    // trajectory (i.e. whose ParentID = 0). 
    if(tj.AlgMod[kKilled]) return -1;
    if(tj.ParentID < 0 || tj.ParentID > (int)tjs.allTraj.size()) return -1;
    if(tj.ParentID == 0) return tj.ID;
    int parid = tj.ParentID;
    for(unsigned short nit = 0; nit < 10; ++nit) {
      if(parid < 1 || parid > (int)tjs.allTraj.size()) break;
      auto& tj = tjs.allTraj[parid - 1];
      if(tj.ParentID < 0 || tj.ParentID > (int)tjs.allTraj.size()) return -1;
      if(tj.ParentID == 0) return tj.ID;
      parid = tj.ParentID;
    } // nit
    return -1;
  } // PrimaryID
  
  /////////////////////////////////////////
  int PrimaryID(const TjStuff& tjs, const PFPStruct& pfp)
  {
    // returns the ID of the most upstream PFParticle (that is not a neutrino)
    
    if(pfp.ParentID == pfp.ID || pfp.ParentID <= 0) return pfp.ID;
    int parid = pfp.ParentID;
    int dtrid = pfp.ID;
    unsigned short nit = 0;
    while(true) {
      auto& parent = tjs.pfps[parid - 1];
      // found a neutrino
      if(parent.PDGCode == 14 || parent.PDGCode == 12) return dtrid;
      // found a primary PFParticle?
      if(parent.ParentID == 0) return parent.ID;
      if(parent.ParentID == parent.ID) return parent.ID;
      dtrid = parent.ID;
      parid = parent.ParentID;
      if(parid < 0) return 0;
      ++nit;
      if(nit == 10) return 0;
    }
  } // PrimaryID

  /////////////////////////////////////////
  std::vector<int> MergeChain(TjStuff& tjs, std::vector<int> mergeList, bool prt)
  {
    // Merges a chain of Tjs that are referenced in mergeList using a more careful approach
    // than that used in EndMerge. The returned vector contains the list of Tjs that were made
    // after merging. An empty vector is returned if there was an error.
    std::vector<int> mergedList;
    if(mergeList.size() < 2) return mergedList;
    if(!tjs.UseAlg[kMergeChain]) return mergedList;
    
    // start a list of vertices attached to these Tjs
    std::vector<unsigned short> vxids;
    // find the extent in Pos[0] = wire number
    float loPos0 = 1E6;
    float hiPos0 = 0;
    CTP_t inCTP = 0;
    for(auto tjid : mergeList) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kKilled] || tj.AlgMod[kMat3D]) return mergedList;
      inCTP = tj.CTP;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] > 0 && std::find(vxids.begin(), vxids.end(), tj.VtxID[end]) == vxids.end()) vxids.push_back(tj.VtxID[end]);
        auto& tp = tj.Pts[tj.EndPt[end]];
        if(tp.Pos[0] < loPos0) loPos0 = tp.Pos[0];
        if(tp.Pos[0] > hiPos0) hiPos0 = tp.Pos[0];
      } // end
    } // tjid
    // ensure that any vertices in the list are near the ends of this range. These
    // need to be preserved and re-attached to the new merged Tj. The existence of one
    // in the middle of the range needs to be dealt with by either deleting the vertex
    // or abandoning the merge
    std::array<unsigned short, 2> vtxID {0,0};
    if(!vxids.empty()) {
      for(auto vxid : vxids) {
        if(vxid == 0 || vxid > tjs.vtx.size()) continue;
        auto& vx2 = tjs.vtx[vxid - 1];
        if(std::abs(vx2.Pos[0] - loPos0) < 4) {
          vtxID[0] = vxid;
        } else if(std::abs(vx2.Pos[0] < hiPos0) < 4) {
          vtxID[1] = vxid;
        } else {
          MakeVertexObsolete(tjs, vx2, true);
        }
      } // vxid
    } // vertices exist
    unsigned int loWire = std::nearbyint(loPos0);
    unsigned int hiWire = std::nearbyint(hiPos0);
    unsigned int mergedSize = hiWire - loWire + 1;
    // Put the position of each TP in the vector if it is unambiguous
    std::vector<Point2_t> pos(mergedSize);
    std::vector<unsigned short> cnt(mergedSize);
    std::vector<int> mtjid(mergedSize);
    std::vector<unsigned short> mtjpt(mergedSize);
    for(auto tjid : mergeList) {
      auto& tj = tjs.allTraj[tjid - 1];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        unsigned short indx = std::nearbyint(tp.Pos[0]) - loWire;
        if(indx > mergedSize - 1) continue;
        pos[indx] = tp.Pos;
        mtjid[indx] = tj.ID;
        mtjpt[indx] = ipt;
        // Preserve Tps that have hits on a wire but none of them are used.
        if(tp.Chg > 0) ++cnt[indx];
      } // ipt
    } // tjid
    // look for overlapping regions; where the count is > 1
    Point2_t goodPos;
    for(unsigned short startIndx = 1; startIndx < mergedSize; ++startIndx) {
      if(cnt[startIndx] == 0) continue;
      if(cnt[startIndx] == 1) {
        goodPos = pos[startIndx];
        continue;
      } // cnt == 1
      // find the end
      unsigned short endIndx = 0;
      for(endIndx = startIndx + 1; endIndx < mergedSize; ++endIndx) if(cnt[endIndx] == 1) break;
      if(endIndx == mergedSize) endIndx = mergedSize - 1;
      // make a bare TP between the previous good position and this good position
      if(prt) mf::LogVerbatim("TC")<<" MergeChain inCTP "<<inCTP<<" overlap from "<<PrintPos(tjs, pos[startIndx])<<"to "<<PrintPos(tjs, pos[endIndx]);
      TrajPoint tmp;
      MakeBareTrajPoint(tjs, goodPos[0], goodPos[1], pos[endIndx][0], pos[endIndx][1], 0, tmp);
      // revise the intervening points with the one that has the best delta
      for(unsigned short jndx = startIndx; jndx < endIndx; ++jndx) {
        float bestDelta = 5;
        for(auto tjid : mergeList) {
          auto& ktj = tjs.allTraj[tjid - 1];
          for(unsigned short kpt = ktj.EndPt[0]; kpt <= ktj.EndPt[1]; ++kpt) {
            auto& ktp = ktj.Pts[kpt];
            unsigned short kndx = std::nearbyint(ktp.Pos[0]) - loWire;
            if(kndx == jndx) {
              float delta = PointTrajDOCA(tjs, pos[kndx][0], pos[kndx][1], tmp);
              if(delta < bestDelta) {
                // overwrite the contents
                cnt[kndx] = 1;
                pos[kndx] = ktp.Pos;
                mtjid[kndx] = ktj.ID;
                mtjpt[kndx] = kpt;
                bestDelta = delta;
              } // best delta
              break;
            } // kndx = jndx
          } // kpt
        } // tjid
      } // jndx
      startIndx = endIndx;
    } // startIndx
/*
    for(unsigned short indx = 0; indx < mergedSize; ++indx) {
      mf::LogVerbatim("TC")<<"wire "<<loWire+indx<<" mtjid "<<mtjid[indx]<<" ipt "<<mtjpt[indx]<<" pos "<<(int)pos[indx][0]<<":"<<(int)pos[indx][1]<<" cnt "<<cnt[indx];
    } // indx
*/
    // Create a new merged Tj
    Trajectory ntj;
    // StoreTraj expects to see a Tj with a negative ID
    ntj.ID = -666;
    ntj.StepDir = 1;
    ntj.CTP = inCTP;
    ntj.Pass = 9;
    ntj.ParentID = -1;
    ntj.VtxID = vtxID;
    ntj.AlgMod[kMergeChain] = true;
    int newID = tjs.allTraj.size() + 1;
    // TODO: Transfer stop flags
    for(auto tjid : mergeList) {
      auto& mtj = tjs.allTraj[tjid - 1];
      auto tHits = PutTrajHitsInVector(mtj, kUsedHits);
      MakeTrajectoryObsolete(tjs, tjid - 1);
      mtj.ParentID = newID;
    } // tjid
    for(unsigned short indx = 0; indx < mergedSize; ++indx) {
      if(mtjid[indx] <= 0) continue;
      auto& mtj = tjs.allTraj[mtjid[indx] - 1];
      auto& mtp = mtj.Pts[mtjpt[indx]];
      ntj.Pts.push_back(mtp);
    } // indx
    // re-assign the hits to ntj
    auto tHits = PutTrajHitsInVector(ntj, kUsedHits);
    for(auto iht : tHits) tjs.fHits[iht].InTraj = ntj.ID;
    SetEndPoints(tjs, ntj);
    if(!StoreTraj(tjs, ntj)) {
      std::cout<<"StoreTraj failed\n";
      return mergedList;
    }
    mergedList.push_back(ntj.ID);
    return mergedList;
  } // MergeChain
  
  /////////////////////////////////////////
  void CheckNoMatchTjs(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // Finds long-ish Tjs that are not 3D-matched and does something about it

    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    
    if(!tjs.UseAlg[kMat3DMerge]) return;
    
    if(prt) mf::LogVerbatim("TC")<<"Inside CheckNoMatchTjs";
    
    std::vector<unsigned short> pfpsToFix;
    for(auto& tj : tjs.allTraj) {
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.Cryostat != cstat) continue;
      if(planeID.TPC != tpc) continue;
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kMat3D]) continue;
      if(tj.AlgMod[kInShower]) continue;
      unsigned short npts = tj.EndPt[1] - tj.EndPt[0];
      if(npts < 10) continue;
      if(prt) mf::LogVerbatim("TC")<<"CNMT: Tj "<<tj.ID<<" npts "<<npts<<" is not matched in 3D. Look for it in matchVec ";
      // look for this Tj in matchvec
      unsigned short firstMS = 0;
      for(firstMS = 0; firstMS < tjs.matchVec.size(); ++firstMS) {
        auto& ms = tjs.matchVec[firstMS];
        if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tj.ID) != ms.TjIDs.end()) break;
      } // ms
      // not found for some reason. Deal with this later
      if(firstMS == tjs.matchVec.size()) continue;
      auto& ms = tjs.matchVec[firstMS];
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" First entry has tjs:";
        for(auto tjid : ms.TjIDs) myprt<<" "<<tjid;
      }
      // skip 2-plane matches for now
      if(ms.TjIDs.size() != tjs.NumPlanes) continue;
      // make a list of the Tjs that were matched
      std::vector<int> matched;
      for(auto tjid : ms.TjIDs) if(tjid != tj.ID) matched.push_back(tjid);
      int btjID = INT_MAX;
      // look for the broken tj in an earlier entry. 
      for(unsigned short ims = 0; ims < firstMS; ++ims) {
        auto& ms = tjs.matchVec[ims];
        if(ms.Count < 3) break;
        if(ms.TjIDs.size() < tjs.NumPlanes) break;
        std::vector<int> leftover(ms.TjIDs.size());
        auto it = std::set_difference(ms.TjIDs.begin(), ms.TjIDs.end(), matched.begin(), matched.end(), leftover.begin());
        leftover.resize(it - leftover.begin());
        if(leftover.size() != 1) continue;
        btjID = leftover[0];
        break;
      } // ims
      if(btjID == INT_MAX) continue;
      unsigned short pfpIndex = GetPFPIndex(tjs, btjID);
      if(pfpIndex > tjs.pfps.size() - 1) continue;
      if(prt) mf::LogVerbatim("TC")<<"  try to merge with broken Tj "<<btjID<<" count "<<ms.Count<<" pfpIndex "<<pfpIndex;
      if(MergeAndStore(tjs, tj.ID - 1, btjID - 1, prt)) {
        auto& newTj = tjs.allTraj[tjs.allTraj.size() - 1];
        newTj.AlgMod[kMat3DMerge] = true;
        // Update the PFParticle TjIDs
        if(pfpIndex < tjs.pfps.size()) {
          auto& pfp = tjs.pfps[pfpIndex];
          std::replace(pfp.TjIDs.begin(), pfp.TjIDs.end(), btjID, newTj.ID);
        }
        // add the pfp index to the list of those that need fixing
        if(std::find(pfpsToFix.begin(), pfpsToFix.end(), pfpIndex) == pfpsToFix.end()) pfpsToFix.push_back(pfpIndex);
        // update matchVec
        for(auto& ms : tjs.matchVec) {
          std::replace(ms.TjIDs.begin(), ms.TjIDs.end(), tj.ID, newTj.ID);
          std::replace(ms.TjIDs.begin(), ms.TjIDs.end(), btjID, newTj.ID);
        } // ms
        // update mallTraj
        std::vector<int> oldTjs(2);
        oldTjs[0] = tj.ID;
        oldTjs[1] = btjID;
        UpdateMatchStructs(tjs, oldTjs, newTj.ID);
        if(prt) mf::LogVerbatim("TC")<<"  success "<<tj.ID<<" merged with "<<btjID<<" -> "<<newTj.ID;
      }
      if(prt && !tj.AlgMod[kKilled]) {
        mf::LogVerbatim myprt("TC");
        myprt<<" CheckNoMatchTjs: No 3D match "<<tj.ID;
        myprt<<" in matchVec with other Tjs";
        for(auto tjid : matched) myprt<<" "<<tjid;
        myprt<<" btjID "<<btjID;
      } // prt
    } // tj
    
    // update the PFParticle end points
    std::cout<<"CNMTj is disabled\n";
/*
    for(auto pfpIndex : pfpsToFix) {
      auto& pfp = tjs.pfps[pfpIndex];
      for(unsigned short end = 0; end < 2; ++end) {
        if(!SetPFPEndPoints(tjs, pfp, end, prt)) {
          if(prt) mf::LogVerbatim("TC")<<" SetPFPEndPoints failed";
          pfp.ID = 0;
          continue;
        }
      } // end
    } // pfpIndex
*/
  } // CheckNoMatchTjs

  /////////////////////////////////////////
  bool DefinePFP(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // This function is called after the 3D matched TjIDs have been specified and optionally
    // a start or end vertex ID. It defins the PFParticle but doesn't store it
    
    if(pfp.PDGCode == 1111) return false;
    if(pfp.TjIDs.size() < 2) return false;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"DPFP: PFP "<<pfp.ID;
      myprt<<" Vx3ID "<<pfp.Vx3ID[0]<<" "<<pfp.Vx3ID[1];
      myprt<<" Tjs";
      for(auto id : pfp.TjIDs) myprt<<" "<<id;
      if(pfp.PDGCode == 1111) myprt<<" This is a shower PFP ";
    }
    
    if(pfp.Vx3ID[0] == 0 && pfp.Vx3ID[1] > 0) {
      std::cout<<"DPFP: pfp "<<pfp.ID<<" end 1 has a vertex but end 0 doesn't. No endpoints defined\n";
      return false;
    }
    
    // set the start and end positions on any call
    for(unsigned short startend = 0; startend < 2; ++startend) {
      if(pfp.Vx3ID[startend] > 0) {
        // 3D vertex exists
        auto& vx3 = tjs.vtx3[pfp.Vx3ID[startend] - 1];
        pfp.XYZ[startend][0] = vx3.X;
        pfp.XYZ[startend][1] = vx3.Y;
        pfp.XYZ[startend][2] = vx3.Z;
      }
    } // startend
    
    if(pfp.Vx3ID[0] > 0) {
      // reverse Tjs so that end 0 is at the start vertex if one exists
      auto& vx3 = tjs.vtx3[pfp.Vx3ID[0] - 1];
      for(auto id : pfp.TjIDs) {
        auto& tj = tjs.allTraj[id - 1];
        unsigned short vtxEnd = USHRT_MAX;
        for(unsigned short tjEnd = 0; tjEnd < 2; ++tjEnd) {
          if(tj.VtxID[tjEnd] == 0) continue;
          auto& vx2 = tjs.vtx[tj.VtxID[tjEnd] - 1];
          if(vx2.Vx3ID == vx3.ID) vtxEnd = tjEnd;
        } // tjEnd
        if(vtxEnd == USHRT_MAX) continue;
        if(vtxEnd != 0) ReverseTraj(tjs, tj); 
//        if(prt) mf::LogVerbatim("TC")<<" Reverse Tj "<<id<<" end0 Pos "<<PrintPos(tjs, tj.Pts[tj.EndPt[0]].Pos);
      } // id
    } // end = 0
 
    // Make a list of TC space points allowing for Tjs that are not in the pfp.TjIDs list in
    // the 3rd plane
    MakePFPTp3s(tjs, pfp, true);
    if(prt) mf::LogVerbatim("TC")<<"DPFP: found "<<pfp.Tp3s.size()<<" space points";
    if(pfp.Tp3s.size() < 2) return false;
    pfp.PDGCode = PDGCodeVote(tjs, pfp.TjIDs, prt);
    // The space points are naturally ordered by increasing X. We want to select the first
    // space point as the start if no start vertex exists so pick a point using other criteria,
    // such as muon direction, TPC boundaries etc. This function will set pfp.XYZ[0].
    // SortByDistanceFrom Start will sort the space points by distance from this point
    if(SetNewStart(tjs, pfp, prt)) SortByDistanceFromStart(tjs, pfp);
    CheckTp3Validity(tjs, pfp, prt);
    pfp.Dir[0] = pfp.Tp3s[0].Dir;
    if(pfp.Vx3ID[1] == 0) pfp.XYZ[1] = pfp.Tp3s[pfp.Tp3s.size() - 1].Pos;
    pfp.Dir[1] = pfp.Tp3s[pfp.Tp3s.size() - 1].Dir;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" TC space points\n";
      for(unsigned short ipt = 0; ipt < pfp.Tp3s.size(); ++ipt) {
        auto tp3 = pfp.Tp3s[ipt];
        myprt<<std::setw(3)<<ipt<<" pos ";
        myprt<<" "<<std::fixed<<std::setprecision(1);
        myprt<<std::setw(6)<<tp3.Pos[0]<<std::setw(6)<<tp3.Pos[1]<<std::setw(6)<<tp3.Pos[2];
        myprt<<" sep "<<std::setprecision(1)<<std::setw(5)<<PosSep(tp3.Pos, pfp.Tp3s[0].Pos);
        float sep2 = -1;
        if(ipt > 0) sep2 = PosSep(tp3.Pos, pfp.Tp3s[ipt - 1].Pos);
        myprt<<" sep2 "<<std::setprecision(2)<<std::setw(5)<<sep2;
        myprt<<" dir ";
        myprt<<" "<<std::setprecision(3)<<std::setw(7)<<tp3.Dir[0]<<std::setw(7)<<tp3.Dir[1]<<std::setw(7)<<tp3.Dir[2];
        myprt<<" dEdx "<<std::setw(4)<<std::setprecision(1)<<tp3.dEdx;
        myprt<<" Err "<<std::setw(4)<<std::setprecision(1)<<tp3.dEdxErr;
        myprt<<" IsValid? "<<tp3.IsValid;
        // Calculate the kink angle at point ipt, using the two points that are
        // +/- 1 cm on either side of that point
        double sep = 1;
        myprt<<" kinkAngle "<<std::setprecision(3)<<std::setw(7)<<KinkAngle(tjs, pfp.Tp3s, ipt, sep);
        myprt<<" tj_ipt";
        for(auto tj2pt : tp3.Tj2Pts) {
          auto& tj = tjs.allTraj[tj2pt.id - 1];
          auto& tp = tj.Pts[tj2pt.ipt];
          myprt<<" "<<tj.ID<<"_"<<PrintPos(tjs, tp);
        } // tj2pt
        myprt<<"\n";
      } // ipt
    } // prt
    FilldEdx(tjs, pfp);
    return true;
  } // DefinePFP
  
  /////////////////////////////////////////
  bool CompatibleMerge(TjStuff& tjs, const Trajectory& tj1, const Trajectory& tj2, bool prt)
  {
    // returns true if the two Tjs are compatible with and end0-end1 merge. This function has many aspects of the
    // compatibility checks done in EndMerge but with looser cuts.
    if(tj1.AlgMod[kKilled] || tj2.AlgMod[kKilled]) return false;
    if(tj1.CTP != tj2.CTP) return false;
    unsigned short end1 = -1, end2 = 0;
    float minLen = PosSep(tj1.Pts[tj1.EndPt[0]].Pos, tj1.Pts[tj1.EndPt[1]].Pos);
    float len2 = PosSep(tj2.Pts[tj2.EndPt[0]].Pos, tj2.Pts[tj2.EndPt[1]].Pos);
    if(len2 < minLen) minLen = len2;
    minLen *= 1.2;
    if(minLen > 10) minLen = 10;
    for(unsigned short e1 = 0; e1 < 2; ++e1) {
      auto& tp1 = tj1.Pts[tj1.EndPt[e1]];
      for(unsigned short e2 = 0; e2 < 2; ++e2) {
        auto& tp2 = tj2.Pts[tj2.EndPt[e2]];
        float sep = PosSep(tp1.Pos, tp2.Pos);
        if(sep < minLen) {
          minLen = sep;
          end1 = e1; end2 = e2;
        }
      } // e2
    } // e1
    if(end1 < 0) return false;
    // require end to end
    if(end2 != 1 - end1) return false;
    
    float overlapFraction = OverlapFraction(tjs, tj1, tj2);
    if(overlapFraction > 0.25) {
      if(prt) mf::LogVerbatim("TC")<<"CM: "<<tj1.ID<<" "<<tj2.ID<<" overlapFraction "<<overlapFraction<<" > 0.25 ";
      return false;
    }
    
    auto& tp1 = tj1.Pts[tj1.EndPt[end1]];
    auto& tp2 = tj2.Pts[tj2.EndPt[end2]];
/* This causes problems with hit collections that have cosmics removed
    if(!SignalBetween(tjs, tp1, tp2, 0.8, false)) {
      if(prt) mf::LogVerbatim("TC")<<"CM: "<<tj1.ID<<" "<<tj2.ID<<" no signal between these points "<<PrintPos(tjs, tp1.Pos)<<" "<<PrintPos(tjs, tp2.Pos);
      return false;
    }
*/
    float doca1 = PointTrajDOCA(tjs, tp1.Pos[0], tp1.Pos[1], tp2);
    float doca2 = PointTrajDOCA(tjs, tp2.Pos[0], tp2.Pos[1], tp1);
    if(doca1 > 2 && doca2 > 2) {
      if(prt) mf::LogVerbatim("TC")<<"CM: "<<tj1.ID<<" "<<tj2.ID<<" Both docas > 2 "<<doca1<<" "<<doca2;
      return false;
    }
    
    float dang = DeltaAngle(tp1.Ang, tp2.Ang);
    if(dang > 2 * tjs.KinkCuts[0]) {
      if(prt) mf::LogVerbatim("TC")<<"CM: "<<tj1.ID<<" "<<tj2.ID<<" dang "<<dang<<" > "<<2 * tjs.KinkCuts[0];
      return false;
    }
    
    return true;
  } // CompatibleMerge

  /////////////////////////////////////////
  float OverlapFraction(TjStuff& tjs, const Trajectory& tj1, const Trajectory& tj2)
  {
    // returns the fraction of wires spanned by two trajectories
    float minWire = 1E6;
    float maxWire = -1E6;
    
    float cnt1 = 0;
    for(auto& tp : tj1.Pts) {
      if(tp.Chg == 0) continue;
      if(tp.Pos[0] < minWire) minWire = tp.Pos[0];
      if(tp.Pos[0] > maxWire) maxWire = tp.Pos[0];
      ++cnt1;
    }
    if(cnt1 == 0) return 0;
    float cnt2 = 0;
    for(auto& tp : tj2.Pts) {
      if(tp.Chg == 0) continue;
      if(tp.Pos[0] < minWire) minWire = tp.Pos[0];
      if(tp.Pos[0] > maxWire) maxWire = tp.Pos[0];
      ++cnt2;
    }
    if(cnt2 == 0) return 0;
    int span = maxWire - minWire;
    if(span <= 0) return 0;
    std::vector<unsigned short> wcnt(span);
    for(auto& tp : tj1.Pts) {
      if(tp.Chg == 0) continue;
      int indx = std::nearbyint(tp.Pos[0] - minWire);
      if(indx < 0 || indx > span - 1) continue;
      ++wcnt[indx];
    }
    for(auto& tp : tj2.Pts) {
      if(tp.Chg == 0) continue;
      int indx = std::nearbyint(tp.Pos[0] - minWire);
      if(indx < 0 || indx > span - 1) continue;
      ++wcnt[indx];
    }
    float cntOverlap = 0;
    for(auto cnt : wcnt) if(cnt > 1) ++cntOverlap;
    if(cnt1 < cnt2) {
      return cntOverlap / cnt1;
    } else {
      return cntOverlap / cnt2;
    }
    
  } // OverlapFraction
  
  /////////////////////////////////////////
  void FilldEdx(TjStuff& tjs, PFPStruct& pfp)
  {
    // Fills the dEdX vector in the match struct. This function should be called after the
    // matched trajectory points are ordered so that dE/dx is calculated at the start of the PFParticle
    if(pfp.ID == 0) return;
    // error check
    bool notgood = false;
    for(unsigned short startend = 0; startend < 2; ++startend) {
      if(pfp.dEdx[startend].size() != tjs.NumPlanes) notgood = true;
      if(pfp.dEdxErr[startend].size() != tjs.NumPlanes) notgood = true;
    }
    if(notgood) {
//      if(prt) mf::LogVerbatim("TC")<<"FilldEdx found inconsistent sizes for dEdx\n";
      return;
    }

    double t0 = 0;
    
    unsigned short numEnds = 2;
    // don't attempt to find dE/dx at the end of a shower
    if(pfp.PDGCode == 1111) numEnds = 1;
    
    unsigned short maxlen = 0;
    for(auto tjID : pfp.TjIDs) {

      Trajectory& tj = tjs.allTraj[tjID - 1];
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      double angleToVert = tjs.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
      for(unsigned short startend = 0; startend < numEnds; ++startend) {
        pfp.dEdx[startend][planeID.Plane] = 0;
        tj.dEdx[startend] = 0;
        double cosgamma = std::abs(std::sin(angleToVert) * pfp.Dir[startend][1] + std::cos(angleToVert) * pfp.Dir[startend][2]);
        if(cosgamma == 0) continue;
        double dx = tjs.geom->WirePitch(planeID) / cosgamma;
        if(dx == 0) continue;
        double dQ = tj.Pts[tj.EndPt[startend]].AveChg;
        if(dQ == 0) continue;
        // convert to dQ/dx
        dQ /= dx;
        double time = tj.Pts[tj.EndPt[startend]].Pos[1] / tjs.UnitsPerTick;
        float dedx = tjs.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
        if(dedx > 999) dedx = 999;
        pfp.dEdx[startend][planeID.Plane] = dedx;
        tj.dEdx[startend] = dedx;
        // ChgRMS is the fractional error
        pfp.dEdxErr[startend][planeID.Plane] = dedx * tj.ChgRMS;
	
      } // startend
      // Grab the best plane iusing the start f 1 < dE/dx < 50 MeV/cm
      if(pfp.dEdx[0][planeID.Plane] > 1 && pfp.dEdx[0][planeID.Plane] < 50) {
        if(tj.Pts.size() > maxlen) {
          maxlen = tj.Pts.size();
          pfp.BestPlane = planeID.Plane;
        }
      } // valid dE/dx

    } // tj
  } // FilldEdX

  /////////////////////////////////////////
  unsigned short AngleRange(TjStuff& tjs, TrajPoint const& tp)
  {
    return AngleRange(tjs, tp.Ang);
  }
  
  /////////////////////////////////////////
  void SetAngleCode(TjStuff& tjs, TrajPoint& tp)
  {
    unsigned short ar = AngleRange(tjs, tp.Ang);
    if(ar == tjs.AngleRanges.size() - 1) {
      // Very large angle
      tp.AngleCode = 2;
    } else if(tjs.AngleRanges.size() > 2 && ar == tjs.AngleRanges.size() - 2) {
      // Large angle
      tp.AngleCode = 1;
    } else {
      // Small angle
      tp.AngleCode = 0;
    }
    
  } // SetAngleCode
  
  /////////////////////////////////////////
  unsigned short AngleRange(TjStuff& tjs, float angle)
  {
    // returns the index of the angle range
    if(angle > M_PI) angle = M_PI;
    if(angle < -M_PI) angle = M_PI;
    if(angle < 0) angle = -angle;
    if(angle > M_PI/2) angle = M_PI - angle;
    for(unsigned short ir = 0; ir < tjs.AngleRanges.size(); ++ir) {
      if(angle < tjs.AngleRanges[ir]) return ir;
    }
    return tjs.AngleRanges.size() - 1;
  } // AngleRange
  
  //////////////////////////////////////////
  void FitTraj(TjStuff& tjs, Trajectory& tj)
  {
    // Jacket around FitTraj to fit the leading edge of the supplied trajectory
    unsigned short originPt = tj.EndPt[1];
    unsigned short npts = tj.Pts[originPt].NTPsFit;
    TrajPoint tpFit;
    unsigned short fitDir = -1;
    FitTraj(tjs, tj, originPt, npts, fitDir, tpFit);
    tj.Pts[originPt] = tpFit;
    
  } // FitTraj
  
  //////////////////////////////////////////
  void FitTraj(TjStuff& tjs, Trajectory& tj, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint& tpFit)
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
    
    if(originPt > tj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"FitTraj: Requesting fit of invalid TP "<<originPt;
      return;
    }
    
    // copy the origin TP into the fit TP
    tpFit = tj.Pts[originPt];
    // Assume that the fit will fail
    tpFit.FitChi = 999;
    if(fitDir < -1 || fitDir > 1) return;
    
    std::vector<double> x, y;
    Point2_t origin = tj.Pts[originPt].HitPos;
    // Use TP position if there aren't any hits on it
    if(tj.Pts[originPt].Chg == 0) origin = tj.Pts[originPt].Pos;
    
    // simple two point case
    if(NumPtsWithCharge(tjs, tj, false) == 2) {
      for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        double xx = tj.Pts[ipt].HitPos[0] - origin[0];
        double yy = tj.Pts[ipt].HitPos[1] - origin[1];
        x.push_back(xx);
        y.push_back(yy);
      } // ii
      if(x.size() != 2) return;
      if(x[0] == x[1]) {
        // Either + or - pi/2
        tpFit.Ang = M_PI/2;
        if(y[1] < y[0]) tpFit.Ang = -tpFit.Ang;
      } else {
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
      SetAngleCode(tjs, tpFit);
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
    if(tj.Pts[originPt].Chg > 0) {
      xx = tj.Pts[originPt].HitPos[0] - origin[0];
      yy = tj.Pts[originPt].HitPos[1] - origin[1];
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      x.push_back(xr);
      y.push_back(yr);
      chgWt = tj.Pts[originPt].ChgPull;
      if(chgWt < 1) chgWt = 1;
      chgWt *= chgWt;
      w.push_back(chgWt * tj.Pts[originPt].HitPosErr2);
    }
    
    // correct npts to account for the origin point
    if(fitDir != 0) --npts;
    
    // step in the + direction first
    if(fitDir != -1) {
      unsigned short cnt = 0;
      for(unsigned short ipt = originPt + 1; ipt < tj.Pts.size(); ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if(chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if(cnt == npts) break;
      } // ipt
    } // fitDir != -1
    
    // step in the - direction next
    if(fitDir != 1 && originPt > 0) {
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = originPt - ii;
        if(ipt > tj.Pts.size() - 1) continue;
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if(chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if(cnt == npts) break;
        if(ipt == 0) break;
      } // ipt
    } // fitDir != -1
    
    // Not enough points to define a line?
    if(x.size() < 2) return;
    
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    
    // weight by the charge ratio and accumulate sums
    double wght;
    for(unsigned short ipt = 0; ipt < x.size(); ++ipt) {
      if(w[ipt] < 0.00001) w[ipt] = 0.00001;
      wght = 1 / w[ipt];
      sum   += wght;
      sumx  += wght * x[ipt];
      sumy  += wght * y[ipt];
      sumx2 += wght * x[ipt] * x[ipt];
      sumy2 += wght * y[ipt] * y[ipt];
      sumxy += wght * x[ipt] * y[ipt];
    }
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0) return;
    // A is the intercept
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // B is the slope
    double B = (sumxy * sum  - sumx * sumy) / delta;
    
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
    if(AngleRange(tjs, tj.Pts[originPt]) > 0) {
      flipDir = std::signbit(tpFit.Dir[1]) != std::signbit(tj.Pts[originPt].Dir[1]);
    } else {
      flipDir = std::signbit(tpFit.Dir[0]) != std::signbit(tj.Pts[originPt].Dir[0]);
    }
    if(flipDir) {
      tpFit.Dir[0] = -tpFit.Dir[0];
      tpFit.Dir[1] = -tpFit.Dir[1];
    }
    tpFit.Ang = atan2(tpFit.Dir[1], tpFit.Dir[0]);
    SetAngleCode(tjs, tpFit);
    //    if(prt) mf::LogVerbatim("TC")<<"FitTraj "<<originPt<<" originPt Dir "<<tj.Pts[originPt].Dir[0]<<" "<<tj.Pts[originPt].Dir[1]<<" rotAngle "<<rotAngle<<" tpFit.Dir "<<tpFit.Dir[0]<<" "<<tpFit.Dir[1]<<" Ang "<<tpFit.Ang<<" flipDir "<<flipDir<<" fit vector size "<<x.size();
    
    // rotate (0, intcpt) into (W,T) coordinates
    tpFit.Pos[0] = -sn * A + origin[0];
    tpFit.Pos[1] =  cs * A + origin[1];
    // force the origin to be at origin[0]
    if(tpFit.AngleCode < 2) MoveTPToWire(tpFit, origin[0]);
    
    if(x.size() < 3) return;
    
    // Calculate chisq/DOF
    double ndof = x.size() - 2;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    if(varnce > 0.) {
      // Intercept error is not used
      //      InterceptError = sqrt(varnce * sumx2 / delta);
      double slopeError = sqrt(varnce * sum / delta);
      tpFit.AngErr = std::abs(atan(slopeError));
    } else {
      tpFit.AngErr = 0.01;
    }
    sum = 0;
    // calculate chisq
    double arg;
    for(unsigned short ii = 0; ii < y.size(); ++ii) {
      arg = y[ii] - A - B * x[ii];
      sum += arg * arg / w[ii];
    }
    tpFit.FitChi = sum / ndof;
    
  } // FitTraj

   ////////////////////////////////////////////////
  void WatchHit(std::string someText, TjStuff& tjs, const unsigned int& wHit, short& wInTraj, const unsigned short& tjID)
  {
    // a temp routine to watch when inTraj changes for the supplied hit index, watchHit
    if(wHit > tjs.fHits.size() - 1) return;
    
    if(tjs.fHits[wHit].InTraj != wInTraj) {
      std::cout<<someText<<" Hit "<<PrintHitShort(tjs.fHits[wHit])<<" was InTraj "<<wInTraj<<" now InTraj "<<tjs.fHits[wHit].InTraj<<" tjID = "<<tjID<<"\n";
      wInTraj = tjs.fHits[wHit].InTraj;
    }
  } // WatchHit
  
  ////////////////////////////////////////////////
  void TagProtons(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    const unsigned int cstat = tpcid.Cryostat;
    const unsigned int tpc = tpcid.TPC;
    std::vector<int> tjlist(1);
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if(planeID.TPC != tpc || planeID.Cryostat != cstat) continue;
      // ignore tagged muons
      if(tj.PDGCode == 13) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] != 0) continue;
        if(!tj.StopFlag[end][kBragg]) continue;
        // check the environment near this end
        tjlist[0] = tj.ID;
        float chgFrac = ChgFracNearPos(tjs, tj.Pts[tj.EndPt[end]].Pos, tjlist);
        if(prt) mf::LogVerbatim("TC")<<"TagProtons: Tj "<<tj.ID<<" Charge fraction near end "<<end<<" "<<chgFrac;
        if(chgFrac > 0.9) tj.PDGCode = 2212;
      } // end
    } // tj
  } // TagProtons
/*
  ////////////////////////////////////////////////
  void TagBragg(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // sets the PDG code to 2212 if there are Bragg peaks on the Tjs
    if(pfp.PDGCode == 11 || pfp.PDGCode == 1111) return;
    
    unsigned short braggCnt0 = 0;
    unsigned short braggCnt1 = 0;
    for(auto& tjID : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjID - 1];
      if(tj.StopFlag[0][kBragg]) ++braggCnt0;
      if(tj.StopFlag[1][kBragg]) ++braggCnt1;
    }
    if(braggCnt0 > 1 || braggCnt1 > 1) pfp.PDGCode = 2212;
    
  } // TagBragg
*/
  ////////////////////////////////////////////////
  void Reverse3DMatchTjs(TjStuff& tjs, PFPStruct& pfp, bool prt)
  {
    // Return true if the 3D matched hits in the trajectories in tjs.pfps are in the wrong order in terms of the
    // physics standpoint, e.g. dQ/dx, muon delta-ray tag, cosmic rays entering the detector, etc. 
    
    // Don't reverse showers
    if(pfp.PDGCode == 1111) return;
    
    bool reverseMe = false;

    // look for stopping Tjs for contained PFParticles
    if(!reverseMe) {
      unsigned short braggCnt0 = 0;
      unsigned short braggCnt1 = 0;
      for(auto& tjID : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjID - 1];
        if(tj.StopFlag[0][kBragg]) ++braggCnt0;
        if(tj.StopFlag[1][kBragg]) ++braggCnt1;
      }
      if(braggCnt0 > 0 || braggCnt1 > 0) {
        pfp.PDGCode = 2212;
        // Vote for a Bragg peak at the beginning. It should be at the end
        if(braggCnt0 > braggCnt1) reverseMe = true;
      } // found a Bragg Peak 
    } // look for stopping Tjs 
    
    if(!reverseMe) return;
    
    // All of the trajectories should be reversed
    for(auto& tjID : pfp.TjIDs) {
      unsigned short itj = tjID - 1;
      Trajectory& tj = tjs.allTraj[itj];
      tj.AlgMod[kMat3D] = false;
      ReverseTraj(tjs, tj);
      tj.AlgMod[kMat3D] = true;
    } // tjID
    // swap the matchVec end info also
    std::swap(pfp.XYZ[0], pfp.XYZ[1]);
    std::swap(pfp.Dir[0], pfp.Dir[1]);
    std::swap(pfp.DirErr[0], pfp.DirErr[1]);
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    
    return;
    
  } // Reverse3DMatchTjs
  
  ////////////////////////////////////////////////
  unsigned short GetPFPIndex(const TjStuff& tjs, int tjID)
  {
    // returns the index into the tjs.matchVec vector of the first 3D match that
    // includes tjID
    if(tjs.pfps.empty()) return USHRT_MAX;
    for(unsigned int ipfp = 0; ipfp < tjs.pfps.size(); ++ipfp) {
      const auto& pfp = tjs.pfps[ipfp];
      if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tjID) != pfp.TjIDs.end()) return ipfp;
    } // indx
    return USHRT_MAX;
  } // GetPFPIndex

  ////////////////////////////////////////////////
  unsigned short MatchVecIndex(const TjStuff& tjs, int tjID)
  {
    // returns the index into the tjs.matchVec vector of the first 3D match that
    // includes tjID
    for(unsigned int ims = 0; ims < tjs.matchVec.size(); ++ims) {
      const auto& ms = tjs.matchVec[ims];
      if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tjID) != ms.TjIDs.end()) return ims;
    } // indx
    return USHRT_MAX;
  } // MatchVecIndex

  ////////////////////////////////////////////////
  PFPStruct CreatePFPStruct(const TjStuff& tjs, const geo::TPCID& tpcid)
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
  void ReleaseHits(TjStuff& tjs, Trajectory& tj)
  {
    // Sets InTraj[] = 0 for all TPs in work. Called when abandoning work
    for(auto& tp : tj.Pts) {
      for(auto iht : tp.Hits) {
        if(tjs.fHits[iht].InTraj == tj.ID) tjs.fHits[iht].InTraj = 0;
      }
    } // tp
    
  } // ReleaseWorkHits
  
  //////////////////////////////////////////
  void UnsetUsedHits(TjStuff& tjs, TrajPoint& tp)
  {
    // Sets InTraj = 0 and UseHit false for all used hits in tp
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(tp.UseHit[ii]) {
        tjs.fHits[tp.Hits[ii]].InTraj = 0;
        tp.UseHit[ii] = false;
      } // UseHit
    } // ii
    tp.Chg = 0;
  } // UnsetUsedHits

  ////////////////////////////////////////////////
  bool StoreTraj(TjStuff& tjs, Trajectory& tj)
  {
    
    if(tj.EndPt[1] <= tj.EndPt[0]) return false;
    
    if(!(tj.StepDir == 1 || tj.StepDir == -1)) {
      mf::LogError("TC")<<"StoreTraj: Invalid StepDir "<<tj.StepDir;
      return false;
    }
    
    if(tjs.allTraj.size() >= USHRT_MAX) {
      mf::LogError("TC")<<"StoreTraj: Too many trajectories "<<tjs.allTraj.size();
      return false;
    }
    
    // This shouldn't be necessary but do it anyway
    SetEndPoints(tjs, tj);
    UpdateAveChg(tjs, tj);
    
    auto& endTp0 = tj.Pts[tj.EndPt[0]];
    auto& endTp1 = tj.Pts[tj.EndPt[1]];
    
    // Calculate the charge near the end and beginning if necessary. This must be a short
    // trajectory. Find the average using 4 points
    if(endTp0.AveChg <= 0) {
      unsigned short cnt = 0;
      float sum = 0;
      for(unsigned short ipt = tj.EndPt[0] + 1; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == 4) break;
      }
      tj.Pts[tj.EndPt[0]].AveChg = sum / (float)cnt;
    }
    if(endTp1.AveChg <= 0) {
      float sum = 0;
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.EndPt[1] - ii;
        if(tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == 4) break;
        if(ipt == 0) break;
      } // ii
      tj.Pts[tj.EndPt[1]].AveChg = sum / (float)cnt;
    } // begin charge == end charge
    
    int trID = tjs.allTraj.size() + 1;

    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(tj.Pts[ipt].UseHit[ii]) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(tjs.fHits[iht].InTraj > 0) {
            mf::LogWarning("TC")<<"StoreTraj: Failed trying to store hit "<<PrintHit(tjs.fHits[iht])<<" in new tjs.allTraj "<<trID<<" but it is used in traj ID = "<<tjs.fHits[iht].InTraj<<" with WorkID "<<tjs.allTraj[tjs.fHits[iht].InTraj-1].WorkID<<" Print and quit";
            PrintTrajectory("SW", tjs, tj, USHRT_MAX);
            ReleaseHits(tjs, tj);
            return false;
          } // error
          tjs.fHits[iht].InTraj = trID;
        }
      } // ii
    } // ipt
    
    // ensure that inTraj is clean for the ID
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].InTraj == tj.ID) {
        mf::LogWarning("TC")<<"StoreTraj: Hit "<<PrintHit(tjs.fHits[iht])<<" thinks it belongs to traj ID "<<tj.ID<<" but it wasn't stored\n";
        PrintTrajectory("SW", tjs, tj, USHRT_MAX);
        return false;
      }
    } // iht
    
    tj.WorkID = tj.ID;
    tj.ID = trID;
    // Don't clobber the ParentID if it was defined by the calling function
    if(tj.ParentID == 0) tj.ParentID = trID;
    // Calculate the overall charge RMS relative to a linear
    UpdateChgRMS(tjs, tj);
    tjs.allTraj.push_back(tj);
//    if(prt) mf::LogVerbatim("TC")<<"StoreTraj trID "<<trID<<" CTP "<<tj.CTP<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
    if(debug.Hit != UINT_MAX) {
      // print out some debug info
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
        for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(iht == debug.Hit) std::cout<<"Debug hit appears in trajectory w WorkID "<<tj.WorkID<<" UseHit "<<tj.Pts[ipt].UseHit[ii]<<"\n";
        } // ii
      } // ipt
    } // debug.Hit ...
    
    return true;
    
  } // StoreTraj
  
  //////////////////////////////////////////
  void UpdateAveChg(TjStuff& tjs, Trajectory& tj)
  {
    
    if(tj.EndPt[1] == 0) return;
    unsigned short lastPt = tj.EndPt[1];
    tj.AveChg = 0;
    tj.Pts[lastPt].AveChg = 0;
    
    // calculate ave charge and charge RMS using hits in the trajectory
    unsigned short ii, ipt, cnt = 0;
    float fcnt, sum = 0;
    float sum2 = 0;
    // Don't include the first point in the average. It will be too
    // low if this is a stopping/starting particle
    for(ii = 0; ii < tj.Pts.size(); ++ii) {
      ipt = tj.EndPt[1] - ii;
      if(ipt == 0) break;
      if(tj.Pts[ipt].Chg == 0) continue;
      ++cnt;
      sum += tj.Pts[ipt].Chg;
      sum2 += tj.Pts[ipt].Chg * tj.Pts[ipt].Chg;
      if(cnt == tjs.NPtsAve) break;
    } // iii
    if(cnt == 0) return;
    fcnt = cnt;
    sum /= fcnt;
    tj.AveChg = sum;
    tj.Pts[lastPt].AveChg = sum;
    // define the first point average charge if necessary
    if(tj.Pts[tj.EndPt[0]].AveChg <= 0) tj.Pts[tj.EndPt[0]].AveChg = sum;
    if(cnt > 3) {
      float arg = sum2 - fcnt * sum * sum;
      if(arg < 0) arg = 0;
      float rms = sqrt(arg / (fcnt - 1));
      // convert this to a normalized RMS
      rms /= sum;
      // don't let the calculated charge RMS dominate the default
      // RMS until it is well known. Start with 50% error on the
      // charge RMS
      float defFrac = 1 / (float)(tj.EndPt[1]);
      tj.ChgRMS = defFrac * 0.5 + (1 - defFrac) * rms;
      if(tj.EndPt[1] > 10) {
        // don't let it get crazy small
        if(tj.ChgRMS < tjs.ChargeCuts[1]) tj.ChgRMS = tjs.ChargeCuts[1];
        // or crazy large
        if(tj.ChgRMS > tjs.ChargeCuts[2]) tj.ChgRMS = tjs.ChargeCuts[2];
      }
      tj.Pts[lastPt].ChgPull = (tj.Pts[lastPt].Chg / tj.AveChg - 1) / tj.ChgRMS;
    } // cnt > 3
  } // UpdateAveChg
  
  ////////////////////////////////////////////////
  void UpdateChgRMS(TjStuff& tjs, Trajectory& tj)
  {
    // Calculates the ChgRMS variable using all points on the trajectory except a few at the end
    double ave = 0;
    double sum2 = 0;
    double cnt = 0;
    for(short ipt = tj.EndPt[0] + 5; ipt < tj.EndPt[1] - 5; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      ave += tp.Chg;
      sum2 += tp.Chg * tp.Chg;
      ++cnt;
    } // tp
    if(cnt < 5) return;
    ave /= cnt;
    sum2 = sum2 - cnt * ave * ave;
    if(sum2 < 0) return;
    tj.ChgRMS = sqrt(sum2 / (cnt - 1));
    tj.ChgRMS /= ave;
  } // UpdateChgRMS
  
  ////////////////////////////////////////////////
  bool InTrajOK(TjStuff& tjs, std::string someText)
  {
    // Check tjs.allTraj -> InTraj associations
    
    unsigned short tID;
    unsigned int iht;
    unsigned short itj = 0;
    std::vector<unsigned int> tHits;
    std::vector<unsigned int> atHits;
    for(auto& tj : tjs.allTraj) {
      // ignore abandoned trajectories
      if(tj.AlgMod[kKilled]) continue;
      tID = tj.ID;
      if(tj.AlgMod[kKilled]) {
        std::cout<<someText<<" ChkInTraj hit size mis-match in tj ID "<<tj.ID<<" AlgBitNames";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) std::cout<<" "<<AlgBitNames[ib];
        std::cout<<"\n";
        continue;
      }
      tHits = PutTrajHitsInVector(tj, kUsedHits);
      if(tHits.size() < 2) {
        std::cout<<someText<<" ChkInTraj: Insufficient hits in traj "<<tj.ID<<"\n";
        PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
        continue;
      }
      std::sort(tHits.begin(), tHits.end());
      atHits.clear();
      for(iht = 0; iht < tjs.fHits.size(); ++iht) {
        if(tjs.fHits[iht].InTraj == tID) atHits.push_back(iht);
      } // iht
      if(atHits.size() < 2) {
        std::cout<<someText<<" ChkInTraj: Insufficient hits in atHits in traj "<<tj.ID<<" Killing it\n";
        tj.AlgMod[kKilled] = true;
        continue;
      }
      if(!std::equal(tHits.begin(), tHits.end(), atHits.begin())) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ChkInTraj failed: inTraj - UseHit mis-match for tj ID "<<tID<<" tj.WorkID "<<tj.WorkID<<" atHits size "<<atHits.size()<<" tHits size "<<tHits.size()<<" in CTP "<<tj.CTP<<"\n";
        myprt<<"AlgMods: ";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
        myprt<<"index     inTraj     UseHit \n";
        for(iht = 0; iht < atHits.size(); ++iht) {
          myprt<<"iht "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]]);
          if(iht < tHits.size()) myprt<<" "<<PrintHit(tjs.fHits[tHits[iht]]);
          if(atHits[iht] != tHits[iht]) myprt<<" <<< "<<atHits[iht]<<" != "<<tHits[iht];
          myprt<<"\n";
        } // iht
        if(tHits.size() > atHits.size()) {
          for(iht = atHits.size(); iht < atHits.size(); ++iht) {
            myprt<<"atHits "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]])<<"\n";
          } // iht
          PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
        } // tHit.size > atHits.size()
        return false;
      }
      // check the VtxID
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] > tjs.vtx.size()) {
          mf::LogVerbatim("TC")<<someText<<" ChkInTraj: Bad VtxID "<<tj.ID;
          std::cout<<someText<<" ChkInTraj: Bad VtxID "<<tj.ID<<" vtx size "<<tjs.vtx.size()<<"\n";
          tj.AlgMod[kKilled] = true;
          PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
          return false;
        }
      } // end
      ++itj;
    } // tj
    return true;
    
  } // InTrajOK
  
  //////////////////////////////////////////
  void CheckTrajBeginChg(TjStuff& tjs, unsigned short itj, bool prt)
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
    if(itj > tjs.allTraj.size() - 1) return;
    auto& tj = tjs.allTraj[itj];
    
    if(!tjs.UseAlg[kBeginChg]) return;
    if(tj.StopFlag[0][kBragg]) return;
    if(tj.AlgMod[kFTBRvProp]) return;
    if(tj.AlgMod[kKilled]) return;
    if(tj.Pts.size() < 20) return;
    
    // look for a large drop between the average charge near the beginning
    float chg2 = tj.Pts[tj.EndPt[0] + 2].AveChg;
    // and the average charge 15 points away
    float chg15 = tj.Pts[tj.EndPt[0] + 15].AveChg;
    if(chg2 < 3 * chg15) return;
    
    // find the point where the charge falls below the mid-point
    float midChg = 0.5 * (chg2 + chg15);
    
    unsigned short breakPt = USHRT_MAX;
    for(unsigned short ipt = tj.EndPt[0] + 3; ipt < 15; ++ipt) {
      float chgm2 = tj.Pts[ipt - 2].Chg;
      if(chgm2 == 0) continue;
      float chgm1 = tj.Pts[ipt - 1].Chg;
      if(chgm1 == 0) continue;
      float chgp1 = tj.Pts[ipt + 1].Chg;
      if(chgp1 == 0) continue;
      float chgp2 = tj.Pts[ipt + 2].Chg;
      if(chgp2 == 0) continue;
      if(chgm2 > midChg && chgm1 > midChg && chgp1 < midChg && chgp2 < midChg) {
        breakPt = ipt;
        break;
      }
    } // breakPt
    if(breakPt == USHRT_MAX) return;
    // Create a vertex at the break point
    VtxStore aVtx;
    aVtx.Pos = tj.Pts[breakPt].Pos;
    aVtx.NTraj = 2;
    aVtx.Pass = tj.Pass;
    aVtx.Topo = 8;
    aVtx.ChiDOF = 0;
    aVtx.CTP = tj.CTP;
    aVtx.ID = tjs.vtx.size() + 1;
    aVtx.Stat[kFixed] = true;
    unsigned short ivx = tjs.vtx.size();
    if(!StoreVertex(tjs, aVtx)) return;
    if(!SplitAllTraj(tjs, itj, breakPt, ivx, prt)) {
      if(prt) mf::LogVerbatim("TC")<<"CTBC: Failed to split trajectory";
      MakeVertexObsolete(tjs, tjs.vtx[ivx], false);
      return;
    }
    SetVx2Score(tjs, prt);
    
    if(prt) mf::LogVerbatim("TC")<<"CTBC: Split Tj "<<tj.ID<<" at "<<PrintPos(tjs, tj.Pts[breakPt].Pos)<<"\n";
    
  } // CheckTrajBeginChg
  
  //////////////////////////////////////////
  void TrimEndPts(TjStuff& tjs, Trajectory& tj, const std::vector<float>& fQualityCuts, bool prt)
  {
    // Trim the hits off the end until there are at least fMinPts consecutive hits at the end
    // and the fraction of hits on the trajectory exceeds fQualityCuts[0]
    // Minimum length requirement accounting for dead wires where - denotes a wire with a point
    // and D is a dead wire. Here is an example with minPts = 3
    //  ---DDDDD--- is OK
    //  ----DD-DD-- is OK
    //  ----DDD-D-- is OK
    //  ----DDDDD-- is not OK
    
    if(!tjs.UseAlg[kTEP]) return;
    
    unsigned short npwc = NumPtsWithCharge(tjs, tj, false);
    unsigned short minPts = fQualityCuts[1];
    if(minPts < 1) return;
    if(npwc < minPts) return;
    
    // handle short tjs
    if(npwc == minPts + 1) {
      unsigned short endPt1 = tj.EndPt[1];
      auto& tp = tj.Pts[endPt1];
      auto& ptp = tj.Pts[endPt1 - 1];
      // remove the last point if the previous point has no charge or if
      // it isn't on the next wire
      float dwire = std::abs(ptp.Pos[0] - tp.Pos[0]);
      if(ptp.Chg == 0 || dwire > 1.1) {
        UnsetUsedHits(tjs, tp);
        SetEndPoints(tjs, tj);
        tj.AlgMod[kTEP] = true;
      }
      return;
    } // short tj
    
    // find the separation between adjacent points, starting at the end
    unsigned short lastPt = 0;
    for(lastPt = tj.EndPt[1]; lastPt > minPts; --lastPt) {
      // check for an error
      if(lastPt == 1) break;
      if(tj.Pts[lastPt].Chg == 0) continue;
      // number of adjacent points on adjacent wires
      unsigned short nadj = 0;
      unsigned short npwc = 0;
      for(unsigned short ipt = lastPt - minPts; ipt < lastPt; ++ipt) {
        if(ipt == 1) break;
        // the current point
        auto& tp = tj.Pts[ipt];
        // the previous point
        auto& ptp = tj.Pts[ipt - 1];
        if(tp.Chg > 0 && ptp.Chg > 0) {
          ++npwc;
          if(std::abs(tp.Pos[0] - ptp.Pos[0]) < 1.5) ++nadj;
        }
//        std::cout<<" "<<PrintPos(tjs, ptp.Pos)<<"_"<<(int)ptp.Chg<<" "<<PrintPos(tjs, tp.Pos)<<"_"<<(int)tp.Chg<<"\n";
      } // ipt
      float ntpwc = NumPtsWithCharge(tjs, tj, true, tj.EndPt[0], lastPt);
      float nwires = std::abs(tj.Pts[tj.EndPt[0]].Pos[0] - tj.Pts[lastPt].Pos[0]) + 1;
      float hitFrac = ntpwc / nwires;
      if(prt) mf::LogVerbatim("TC")<<"TEP: ID "<<tj.ID<<" lastPt "<<lastPt<<" npwc "<<npwc<<" nadj "<<nadj<<" hitFrac "<<hitFrac;
      if(hitFrac > fQualityCuts[0] && npwc == minPts && nadj == minPts) break;
    } // lastPt
    
    // Nothing needs to be done
    if(lastPt == tj.EndPt[1]) return;
    
    // clear the points after lastPt
    for(unsigned short ipt = lastPt + 1; ipt <= tj.EndPt[1]; ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
    SetEndPoints(tjs, tj);
//    tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kTEP] = true;
    if(prt) PrintTrajectory("TEPo", tjs, tj, USHRT_MAX);
    
  } // TrimEndPts
/*  
  //////////////////////////////////////////
  void TrimEndPts(TjStuff& tjs, Trajectory& tj, const std::vector<float>& fQualityCuts, bool prt)
  {
    // Trim the hits off the end until there are at least fMinPts consecutive hits at the end
    // and the fraction of hits on the trajectory exceeds fQualityCuts[0]
    // Minimum length requirement accounting for dead wires where - denotes a wire with a point
    // and D is a dead wire. Here is an example with minPts = 3
    //  ---DDDDD--- is OK
    //  ----DD-DD-- is OK
    //  ----DDD-D-- is OK
    //  ----DDDDD-- is not OK
    
    if(!tjs.UseAlg[kTEP]) return;
    
    float npwc = NumPtsWithCharge(tjs, tj, false);
    if(npwc < fQualityCuts[1] + 1) return;
    
    // consider short Tjs that the code below doesn't handle
    if(npwc == fQualityCuts[1] + 1) {
      float sep = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
      std::cout<<"TEP: "<<tj.ID<<" npwc "<<npwc<<" sep "<<sep<<"\n";
    } // short tj
    
    unsigned short minPts = fQualityCuts[1];
    float maxPtSep = minPts + 2;

    if(prt) {
      mf::LogVerbatim("TC")<<"TrimEndPts: minPts "<<minPts<<" required. maxPtSep "<<maxPtSep<<" Minimum hit fraction "<<fQualityCuts[0];
      if(tj.Pts.size() < 50) PrintTrajectory("TEPi", tjs, tj, USHRT_MAX);
    }

    unsigned short newEndPt = tj.EndPt[1];
    unsigned short nPtsWithCharge;
    float hitFrac = 0;
    while(newEndPt > minPts) {
      nPtsWithCharge = 0;
      if(tj.Pts[newEndPt].Chg == 0) {
        --newEndPt;
        continue;
      }
      for(unsigned short jj = 0; jj < minPts && jj<= newEndPt; ++jj) {
        unsigned short jpt = newEndPt - jj;
        if(tj.Pts[jpt].Chg > 0) ++nPtsWithCharge; 
        if(jpt < minPts) break; //TY: so trajectory with 4 points won't be killed
      } // jj
      
      float ptSep = std::abs(tj.Pts[newEndPt - minPts].Pos[0] - tj.Pts[newEndPt].Pos[0]);
      if(prt) mf::LogVerbatim("TC")<<" newEndPt "<<newEndPt<<" ptSep "<<ptSep<<" nPtsWithCharge "<<nPtsWithCharge;
      // allow only one dead wire at the end
      if(nPtsWithCharge == minPts && ptSep < maxPtSep) {
        // minPts consecutive points have charge. Check the TP Chg fraction
        float npwc = NumPtsWithCharge(tjs, tj, true, tj.EndPt[0], newEndPt);
        float nwires = std::abs(tj.Pts[tj.EndPt[0]].Pos[0] - tj.Pts[newEndPt].Pos[0]) + 1;
        hitFrac = npwc / nwires;
        if(prt) mf::LogVerbatim("TC")<<" check hitFrac "<<newEndPt<<" nwires "<<(int)nwires<<" npwc "<<(int)npwc<<" hitFrac "<<hitFrac;
        if(hitFrac > fQualityCuts[0]) break;
        newEndPt -= minPts;
      }
      --newEndPt;
      // check for a serious failure
      if(newEndPt > tj.EndPt[1]) return;
    } // newEndPt

    // passed the cuts with no modifications
    if(newEndPt == tj.EndPt[1]) return;

    // newEndPt is now the last point that satisfies these conditions
    // dead wire check
    nPtsWithCharge = 0;
    unsigned short nConsecutivePts = 0;
    for(unsigned short jj = 0; jj < minPts; ++jj) {
      unsigned short jpt = newEndPt - jj;
      if(tj.Pts[jpt].Chg > 0) ++nPtsWithCharge;
      if(jj > 0 && std::abs(tj.Pts[jpt+1].Pos[0] - tj.Pts[jpt].Pos[0]) < 1.5) ++nConsecutivePts;
      if(jpt == 0) break;
    } // jj
    
    if(prt) mf::LogVerbatim("TC")<<" newEndPt "<<newEndPt<<" nConsecutivePts "<<nConsecutivePts<<" Required "<<minPts - 1;
    
    // lop off the last point if the consecutive point condition isn't met and re-calculate
    if(nConsecutivePts < minPts - 1 && newEndPt > minPts) {
      --newEndPt;
      nPtsWithCharge = 0;
      unsigned short nConsecutivePts = 0;
      for(unsigned short jj = 0; jj < minPts; ++jj) {
        unsigned short jpt = newEndPt - jj;
        if(tj.Pts[jpt].Chg > 0) ++nPtsWithCharge;
        if(jj > 0 && std::abs(tj.Pts[jpt+1].Pos[0] - tj.Pts[jpt].Pos[0]) < 1.5) ++nConsecutivePts;
        if(jpt == 0) break;
      } // jj
      if(prt) mf::LogVerbatim("TC")<<"   newEndPt "<<newEndPt<<" nConsecutivePts "<<nConsecutivePts<<" Required "<<minPts - 1;
    }
    
    if(newEndPt < minPts) {
      tj.AlgMod[kKilled] = true;
      return;
    }
    
    float nwires = std::abs(tj.Pts[tj.EndPt[0]].Pos[0] - tj.Pts[newEndPt].Pos[0]) + 1;
    npwc = NumPtsWithCharge(tjs, tj, true, tj.EndPt[0], newEndPt);
    hitFrac = npwc / nwires;
    
    if(hitFrac < fQualityCuts[0]) tj.AlgMod[kKilled] = true;
    if(prt) mf::LogVerbatim("TC")<<" Old endpoint "<<tj.EndPt[1]<<" newEndPt "<<newEndPt<<" nwires "<<nwires<<" npwc "<<npwc<<" nConsecutivePts "<<nConsecutivePts<<" hitFrac "<<hitFrac<<" Killed? "<<tj.AlgMod[kKilled];
    
    // failed the cuts
    if(tj.AlgMod[kKilled]) return;
    
    // modifications required
    tj.EndPt[1] = newEndPt;    
    for(unsigned short ipt = newEndPt + 1; ipt < tj.Pts.size(); ++ipt) UnsetUsedHits(tjs, tj.Pts[ipt]);
    SetEndPoints(tjs, tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kTEP] = true;
    if(prt) PrintTrajectory("TEPo", tjs, tj, USHRT_MAX);
    
  } // TrimEndPts
*/
  /////////////////////////////////////////
  bool SignalBetween(TjStuff& tjs, const TrajPoint& tp1, const TrajPoint& tp2, const float& MinWireSignalFraction, bool prt)
  {
    // Returns true if there is a signal on > MinWireSignalFraction of the wires between tp1 and tp2.
    if(MinWireSignalFraction == 0) return true;
    
    int fromWire = std::nearbyint(tp1.Pos[0]);
    int toWire = std::nearbyint(tp2.Pos[0]);
    
    if(fromWire == toWire) {
      TrajPoint tp = tp1;
      // check for a signal midway between
      tp.Pos[1] = 0.5 * (tp1.Pos[1] + tp2.Pos[1]);
      if(prt) mf::LogVerbatim("TC")<<" SignalBetween fromWire = toWire = "<<fromWire<<" SignalAtTp? "<<SignalAtTp(tjs, tp);
      return SignalAtTp(tjs, tp);
    }
    // define a trajectory point located at tp1 that has a direction towards tp2
    TrajPoint tp;
    if(!MakeBareTrajPoint(tjs, tp1, tp2, tp)) return true;
    return SignalBetween(tjs, tp, toWire, MinWireSignalFraction, prt);
  } // SignalBetween

  /////////////////////////////////////////
  bool SignalBetween(TjStuff& tjs, TrajPoint tp, float toPos0, const float& MinWireSignalFraction, bool prt)
  {
    // Returns true if there is a signal on > MinWireSignalFraction of the wires between tp and toPos0.
    // Note that this uses the direction vector of the tp
    
    if(MinWireSignalFraction == 0) return true;
    
    int fromWire = std::nearbyint(tp.Pos[0]);
    int toWire = std::nearbyint(toPos0);
    
    if(fromWire == toWire) {
      if(prt) mf::LogVerbatim("TC")<<" SignalBetween fromWire = toWire = "<<fromWire<<" SignalAtTp? "<<SignalAtTp(tjs, tp);
      return SignalAtTp(tjs, tp);
    }
    
    int nWires = abs(toWire - fromWire) + 1;
    
    unsigned short maxWiresNoSignal = (1 - MinWireSignalFraction) * nWires;
    if(std::abs(tp.Dir[0]) < 0.001) tp.Dir[0] = 0.001;
    float stepSize = std::abs(1/tp.Dir[0]);
    // ensure that we step in the right direction
    if(toWire > fromWire && tp.Dir[0] < 0) stepSize = -stepSize;
    if(toWire < fromWire && tp.Dir[0] > 0) stepSize = -stepSize;
    unsigned short nsig = 0;
    unsigned short num = 0;
    unsigned short nmissed = 0;
    for(unsigned short cnt = 0; cnt < nWires; ++cnt) {
      ++num;
      if(SignalAtTp(tjs, tp)) {
        ++nsig;
      } else {
        ++nmissed;
        if(nmissed == maxWiresNoSignal) return false;
      }
      tp.Pos[0] += tp.Dir[0] * stepSize;
      tp.Pos[1] += tp.Dir[1] * stepSize;
    } // cnt
    float sigFrac = (float)nsig / (float)nWires;
    if(prt) mf::LogVerbatim("TC")<<"  SignalBetween fromWire "<<fromWire<<" toWire "<<toWire<<" nWires "<<nWires<<" nsig "<<nsig<<" "<<sigFrac;
    return (sigFrac >= MinWireSignalFraction);
    
  } // SignalBetween
  
  ////////////////////////////////////////////////
  bool TrajHitsOK(TjStuff& tjs, const std::vector<unsigned int>& iHitsInMultiplet, const std::vector<unsigned int>& jHitsInMultiplet)
  {
    // Hits (assume to be on adjacent wires have an acceptable signal overlap
    
    if(iHitsInMultiplet.empty() || jHitsInMultiplet.empty()) return false;
    
    float sum;
    float cvI = HitsPosTick(tjs, iHitsInMultiplet, sum, kAllHits);
    float minI = 1E6;
    float maxI = 0;
    for(auto& iht : iHitsInMultiplet) {
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - 3 * rms;
      if(arg < minI) minI = arg;
      arg = cv + 3 * rms;
      if(arg > maxI) maxI = arg;
    }
    
    float cvJ = HitsPosTick(tjs, jHitsInMultiplet, sum, kAllHits);
    float minJ = 1E6;
    float maxJ = 0;
    for(auto& jht : jHitsInMultiplet) {
      float cv = tjs.fHits[jht].PeakTime;
      float rms = tjs.fHits[jht].RMS;
      float arg = cv - 3 * rms;
      if(arg < minJ) minJ = arg;
      arg = cv + 3 * rms;
      if(arg > maxJ) maxJ = arg;
    }
    
    if(cvI < cvJ) {
      if(maxI > minJ) return true;
    } else {
      if(minI < maxJ) return true;
    }
    return false;
  } // TrajHitsOK
  
  /////////////////////////////////////////
  bool TrajHitsOK(TjStuff& tjs, const unsigned int iht, const unsigned int jht)
  {
    // ensure that two adjacent hits have an acceptable overlap
    if(iht > tjs.fHits.size() - 1) return false;
    if(jht > tjs.fHits.size() - 1) return false;
    // require that they be on adjacent wires
    TCHit& ihit = tjs.fHits[iht];
    TCHit& jhit = tjs.fHits[jht];
    unsigned int iwire = ihit.ArtPtr->WireID().Wire;
    unsigned int jwire = jhit.ArtPtr->WireID().Wire;
    if(abs(iwire - jwire) > 1) return false;
    if(ihit.PeakTime > jhit.PeakTime) {
      float minISignal = ihit.PeakTime - 3 * ihit.RMS;
      float maxJSignal = jhit.PeakTime + 3 * ihit.RMS;
      if(maxJSignal > minISignal) return true;
    } else {
      float maxISignal = ihit.PeakTime + 3 * ihit.RMS;
      float minJSignal = jhit.PeakTime - 3 * ihit.RMS;
      if(minJSignal > maxISignal) return true;
    }
    return false;
  } // TrajHitsOK

  ////////////////////////////////////////////////
  float ExpectedHitsRMS(TjStuff& tjs, const TrajPoint& tp)
  {
    // returns the expected RMS of hits for the trajectory point in ticks
    if(std::abs(tp.Dir[0]) > 0.001) {
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      return 1.5 * tjs.AveHitRMS[planeID.Plane] + 2 * std::abs(tp.Dir[1]/tp.Dir[0])/tjs.UnitsPerTick;
    } else {
      return 500;
    }
  } // ExpectedHitsRMS

  /////////////////////////////////////////
  bool SignalAtTp(TjStuff& tjs, const TrajPoint& tp)
  {
    // returns true if there is a hit near tp.Pos
    
    if(tp.Pos[0] < 0) return false;
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned int ipl = planeID.Plane;
    if(wire >= tjs.NumWires[ipl]) return false;
    if(tp.Pos[1] > tjs.MaxPos1[ipl]) return false;
    // Assume dead wires have a signal
    if(tjs.WireHitRange[ipl][wire].first == -1) return true;
    float projTick = (float)(tp.Pos[1] / tjs.UnitsPerTick);
    // estimate the tick range for non-zero angle
//    float tickRange = ExpectedHitsRMS(tjs, tp);
    float tickRange = 0;
    if(std::abs(tp.Dir[1]) != 0) {
      tickRange = std::abs(0.5 / tp.Dir[1]) / tjs.UnitsPerTick;
      // don't let it get too large
      if(tickRange > 40) tickRange = 40;
    }
    float loTpTick = projTick - tickRange;
    float hiTpTick = projTick + tickRange;
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
    
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      TCHit& hit = tjs.fHits[iht];
      if(projTick < hit.PeakTime) {
        float loHitTick = hit.PeakTime - 3 * hit.RMS;
        if(hiTpTick > loHitTick) return true;
      } else {
        float hiHitTick = hit.PeakTime + 3 * hit.RMS;
        if(loTpTick < hiHitTick) return true;
      }
    } // iht
    return false;
    
  } // SignalAtTp

  //////////////////////////////////////////
  float TpSumHitChg(TjStuff& tjs, TrajPoint const& tp){
    float totchg = 0;
    for (size_t i = 0; i<tp.Hits.size(); ++i){
      if (!tp.UseHit[i]) continue;
      totchg += tjs.fHits[tp.Hits[i]].Integral;
    }
    return totchg;
  } // TpSumHitChg
/*
  //////////////////////////////////////////
  bool CheckHitClusterAssociations(TjStuff& tjs)
  {
    // check hit - cluster associations
    
    if(tjs.fHits.size() != tjs.inClus.size()) {
      mf::LogWarning("TC")<<"CHCA: Sizes wrong "<<tjs.fHits.size()<<" "<<tjs.inClus.size();
      return false;
    }
    
    // check cluster -> hit association
    for(unsigned short icl = 0; icl < tjs.tcl.size(); ++icl) {
      if(tjs.tcl[icl].ID <= 0) continue;
      int clID = tjs.tcl[icl].ID;
      for(unsigned short ii = 0; ii < tjs.tcl[icl].tclhits.size(); ++ii) {
        unsigned int iht = tjs.tcl[icl].tclhits[ii];
        if(iht > tjs.fHits.size() - 1) {
          mf::LogWarning("CC")<<"CHCA: Bad tclhits index "<<iht<<" tjs.fHits size "<<tjs.fHits.size();
          return false;
        } // iht > tjs.fHits.size() - 1
        if(tjs.inClus[iht] != clID) {
          mf::LogError("TC")<<"CHCA: Bad cluster -> hit association. clID "<<clID<<" hit "<<PrintHit(tjs.fHits[iht])<<" tjs.inClus "<<tjs.inClus[iht]<<" CTP "<<tjs.tcl[icl].CTP;
          return false;
        }
      } // ii
    } // icl
    
    // check hit -> cluster association
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.inClus[iht] <= 0) continue;
      unsigned int icl = tjs.inClus[iht] - 1;
      // see if the cluster is obsolete
      if(tjs.tcl[icl].ID == 0) {
        mf::LogError("TC")<<"CHCA: Hit "<<PrintHit(tjs.fHits[iht])<<" associated with an obsolete cluster tjs.tcl[icl].ID "<<tjs.tcl[icl].ID;
        return false;
      }
      if (std::find(tjs.tcl[icl].tclhits.begin(), tjs.tcl[icl].tclhits.end(), iht) == tjs.tcl[icl].tclhits.end()) {
        mf::LogError("TC")<<"CHCA: Hit "<<":"<<PrintHit(tjs.fHits[iht])<<" -> tjs.inClus "<<tjs.inClus[iht]<<" but isn't in tjs.tcl[icl].ID "<<tjs.tcl[icl].ID<<" list of hits. icl "<<icl<<" iht "<<iht;
        for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
          if(tjs.allTraj[itj].ClusterIndex == icl) mf::LogError("TC")<<"CHCA: Cluster index "<<icl<<" found in traj ID "<<tjs.allTraj[itj].ID;
        } // itj
        PrintAllTraj("CHCA", tjs, debug, USHRT_MAX, USHRT_MAX);
        return false;
      }
    } // iht
    
    return true;
    
  } // CheckHitClusterAssociations()
*/
  //////////////////////////////////////////
  unsigned short NumPtsWithCharge(TjStuff& tjs, const Trajectory& tj, bool includeDeadWires)
  {
    unsigned short firstPt = tj.EndPt[0];
    unsigned short lastPt = tj.EndPt[1];
    return NumPtsWithCharge(tjs, tj, includeDeadWires, firstPt, lastPt);
  }
  
  //////////////////////////////////////////
  unsigned short NumPtsWithCharge(TjStuff& tjs, const Trajectory& tj, bool includeDeadWires, unsigned short firstPt, unsigned short lastPt)
  {
    unsigned short ntp = 0;
    for(unsigned short ipt = firstPt; ipt <= lastPt; ++ipt) if(tj.Pts[ipt].Chg > 0) ++ntp;
    // Add the count of deadwires
    if(includeDeadWires) ntp += DeadWireCount(tjs, tj.Pts[firstPt], tj.Pts[lastPt]);
    return ntp;
  } // NumPtsWithCharge
  
  //////////////////////////////////////////
  float DeadWireCount(TjStuff& tjs, const TrajPoint& tp1, const TrajPoint& tp2)
  {
    return DeadWireCount(tjs, tp1.Pos[0], tp2.Pos[0], tp1.CTP);
  } // DeadWireCount
  
  //////////////////////////////////////////
  float DeadWireCount(TjStuff& tjs, const float& inWirePos1, const float& inWirePos2, CTP_t tCTP)
  {
    if(inWirePos1 < -0.4 || inWirePos2 < -0.4) return 0;
    unsigned int inWire1 = std::nearbyint(inWirePos1);
    unsigned int inWire2 = std::nearbyint(inWirePos2);
    geo::PlaneID planeID = DecodeCTP(tCTP);
    unsigned short plane = planeID.Plane;
    if(inWire1 > tjs.NumWires[plane] || inWire2 > tjs.NumWires[plane]) return 0;
    if(inWire1 > inWire2) {
      // put in increasing order
      unsigned int tmp = inWire1;
      inWire1 = inWire2;
      inWire2 = tmp;
    } // inWire1 > inWire2
    ++inWire2;
    unsigned int wire, ndead = 0;
    for(wire = inWire1; wire < inWire2; ++wire) if(tjs.WireHitRange[plane][wire].first == -1) ++ndead;
    return ndead;
  } // DeadWireCount

  ////////////////////////////////////////////////
  unsigned short PDGCodeIndex(TjStuff& tjs, int PDGCode)
  {
    unsigned short pdg = abs(PDGCode);
    if(pdg == 11) return 0; // electron
    if(pdg == 22) return 0; // call photons electrons
    if(pdg == 13) return 1; // muon
    if(pdg == 211) return 2; // pion
    if(pdg == 321) return 3; // kaon
    if(pdg == 2212) return 4; // proton
    
//    std::cout<<"PDGCodeIndex: unknown code "<<PDGCode<<"\n";
    
    return USHRT_MAX;
    
  } // PDGCodeIndex

  ////////////////////////////////////////////////
  void MakeTrajectoryObsolete(TjStuff& tjs, unsigned int itj)
  {
    // Note that this does not change the state of UseHit to allow
    // resurrecting the trajectory later (RestoreObsoleteTrajectory)
    if(itj > tjs.allTraj.size() - 1) return;
    unsigned int iht;
    for(auto& tp : tjs.allTraj[itj].Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        iht = tp.Hits[ii];
        if(tjs.fHits[iht].InTraj == tjs.allTraj[itj].ID) tjs.fHits[iht].InTraj = 0;
      } // ii
    } // tp
    tjs.allTraj[itj].AlgMod[kKilled] = true;
  } // MakeTrajectoryObsolete
  
  ////////////////////////////////////////////////
  void RestoreObsoleteTrajectory(TjStuff& tjs, unsigned int itj)
  {
    if(itj > tjs.allTraj.size() - 1) return;
    if(!tjs.allTraj[itj].AlgMod[kKilled]) {
      mf::LogWarning("TC")<<"RestoreObsoleteTrajectory: Trying to restore not-obsolete trajectory "<<tjs.allTraj[itj].ID;
      return;
    }
    unsigned int iht;
    for(auto& tp : tjs.allTraj[itj].Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) {
          iht = tp.Hits[ii];
          if(tjs.fHits[iht].InTraj == 0) {
            tjs.fHits[iht].InTraj = tjs.allTraj[itj].ID;
          }
        }
      } // ii
    } // tp
    tjs.allTraj[itj].AlgMod[kKilled] = false;
  } // RestoreObsoleteTrajectory
  
  //////////////////////////////////////////
  void MergeGhostTjs(TjStuff& tjs, CTP_t inCTP)
  {
    // Merges short Tjs that share many hits with a longer Tj
    if(!tjs.UseAlg[kMrgGhost]) return;
    
    for(auto& shortTj : tjs.allTraj) {
      if(shortTj.AlgMod[kKilled]) continue;
      if(shortTj.CTP != inCTP) continue;
      unsigned short spts = shortTj.EndPt[1] - shortTj.EndPt[0];
      if(spts > 20) continue;
      // ignore delta rays
      if(shortTj.PDGCode == 11) continue;
      // ignore InShower Tjs
      if(shortTj.AlgMod[kInShower]) continue;
      auto tjhits = PutTrajHitsInVector(shortTj, kAllHits);
      if(tjhits.empty()) continue;
      std::vector<int> tids;
      std::vector<unsigned short> tcnt;
      for(auto iht : tjhits) {
        auto& hit = tjs.fHits[iht];
        if(hit.InTraj <= 0) continue;
        if(hit.InTraj == shortTj.ID) continue;
        unsigned short indx = 0;
        for(indx = 0; indx < tids.size(); ++indx) if(hit.InTraj == tids[indx]) break;
        if(indx == tids.size()) {
          tids.push_back(hit.InTraj);
          tcnt.push_back(1);
        } else {
          ++tcnt[indx];
        }
      } // iht
      if(tids.empty()) continue;
      // find the max count for Tjs that are longer than this one
      unsigned short maxcnt = 0;
      unsigned short ltjID = 0;
      for(unsigned short indx = 0; indx < tids.size(); ++indx) {
        if(tcnt[indx] > maxcnt) {
          auto& ltj = tjs.allTraj[tids[indx] - 1];
          unsigned short lpts = ltj.EndPt[1] - ltj.EndPt[0];
          if(lpts < spts) continue;
          maxcnt = tcnt[indx];
          ltjID = tids[indx];
        }
      } // indx
      float hitFrac = (float)maxcnt / (float)tjhits.size();
      if(hitFrac < 0.1) continue;
      std::cout<<"MergeGhostTjs: tj "<<shortTj.ID<<" ghost of "<<ltjID<<"?  cnt "<<maxcnt<<" hitFrac "<<hitFrac<<"\n";
    } // shortTj
  } // MergeGhostTjs

  //////////////////////////////////////////
  bool SplitAllTraj(TjStuff& tjs, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt)
  {
    // Splits the trajectory itj in the tjs.allTraj vector into two trajectories at position pos. Splits
    // the trajectory and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
    
    if(itj > tjs.allTraj.size()-1) return false;
    if(pos < tjs.allTraj[itj].EndPt[0] + 1 || pos > tjs.allTraj[itj].EndPt[1] - 1) return false;
    if(ivx != USHRT_MAX && ivx > tjs.vtx.size() - 1) return false;
    
    Trajectory& tj = tjs.allTraj[itj];
    
    // Reset the PDG Code if we are splitting a tagged muon
    bool splittingMuon = (tj.PDGCode == 13);
    if(splittingMuon) tj.PDGCode = 0;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"SplitAllTraj: Split Tj ID "<<tj.ID<<" at point "<<pos;
      if(ivx < tjs.vtx.size()) myprt<<" with Vtx ID "<<tjs.vtx[ivx].ID;
    }

    // ensure that there will be at least 3 TPs on each trajectory
    unsigned short ipt, ii, ntp = 0;
    for(ipt = 0; ipt < pos; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) {
      if(prt) mf::LogVerbatim("TC")<<" Split point to small at begin "<<ntp<<" pos "<<pos<<" ID ";
      return false;
    }
    ntp = 0;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) {
      if(prt) mf::LogVerbatim("TC")<<" Split point too small at end "<<ntp<<" pos "<<pos<<" EndPt "<<tj.EndPt[1];
      return false;
    }
    
    // make a copy that will become the Tj after the split point
    Trajectory newTj = tj;
    newTj.ID = tjs.allTraj.size() + 1;
    // make another copy in case something goes wrong
    Trajectory oldTj = tj;
    
    // Leave the first section of tj in place. Re-assign the hits
    // to the new trajectory
    unsigned int iht;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(!tj.Pts[ipt].UseHit[ii]) continue;
        iht = tj.Pts[ipt].Hits[ii];
        // This shouldn't happen but check anyway
        if(tjs.fHits[iht].InTraj != tj.ID) continue;
        tjs.fHits[iht].InTraj = newTj.ID;
        tj.Pts[ipt].UseHit[ii] = false;
      } // ii
    } // ipt
    SetEndPoints(tjs, tj);
    if(splittingMuon) SetPDGCode(tjs, tj);
    
    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if(eraseSize > newTj.Pts.size() - 1) {
      tj = oldTj;
      return false;
    }
    
    if(ivx < tjs.vtx.size()) tj.VtxID[1] = tjs.vtx[ivx].ID;
    tj.AlgMod[kSplit] = true;
    if(prt) {
      mf::LogVerbatim("TC")<<" Splitting trajectory ID "<<tj.ID<<" new EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
    }
    
    // erase the TPs at the beginning of the new trajectory
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for(ipt = 0; ipt < 3; ++ipt) {
      for(ii = 0; ii < newTj.Pts[ipt].Hits.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
      newTj.Pts[ipt].Chg = 0;
    } // ipt
    SetEndPoints(tjs, newTj);
    if(splittingMuon) SetPDGCode(tjs, newTj);
    if(ivx < tjs.vtx.size()) newTj.VtxID[0] = tjs.vtx[ivx].ID;
    newTj.AlgMod[kSplit] = true;
    newTj.ParentID = -1;
    tjs.allTraj.push_back(newTj);

    if(prt) {
      mf::LogVerbatim("TC")<<"  newTj ID "<<newTj.ID<<" EndPts "<<newTj.EndPt[0]<<" to "<<newTj.EndPt[1];
    }
    return true;
    
  } // SplitAllTraj
  
  //////////////////////////////////////////
  void TrajPointTrajDOCA(TjStuff& tjs, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep)
  {
    // Finds the point, ipt, on trajectory tj that is closest to trajpoint tp
    float best = minSep * minSep;
    closePt = USHRT_MAX;
    float dw, dt, dp2;
    unsigned short ipt;
    for(ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      dw = tj.Pts[ipt].Pos[0] - tp.Pos[0];
      dt = tj.Pts[ipt].Pos[1] - tp.Pos[1];
      dp2 = dw * dw + dt * dt;
      if(dp2 < best) {
        best = dp2;
        closePt = ipt;
      }
    } // ipt
    minSep = sqrt(best);
  } // TrajPointTrajDOCA
  
  //////////////////////////////////////////
  bool TrajTrajDOCA(TjStuff& tjs, Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
  {
    return TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, minSep, false);
  } // TrajTrajDOCA
  
  //////////////////////////////////////////
  bool TrajTrajDOCA(TjStuff& tjs, Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep, bool considerDeadWires)
  {
    // Find the Distance Of Closest Approach between two trajectories less than minSep
    // start with some rough cuts to minimize the use of the more expensive checking. This
    // function returns true if the DOCA is less than minSep
    for(unsigned short iwt = 0; iwt < 2; ++iwt) {
      // Apply box cuts on the ends of the trajectories
      // The Lo/Hi wire(time) at each end of tj1
      float wt0 = tj1.Pts[tj1.EndPt[0]].Pos[iwt];
      float wt1 = tj1.Pts[tj1.EndPt[1]].Pos[iwt];
      float lowt1 = wt0;
      float hiwt1 = wt1;
      if(wt1 < lowt1) {
        lowt1 = wt1;
        hiwt1 = wt0;
      }
      // The Lo/Hi wire(time) at each end of tj2
      wt0 = tj2.Pts[tj2.EndPt[0]].Pos[iwt];
      wt1 = tj2.Pts[tj2.EndPt[1]].Pos[iwt];
      float lowt2 = wt0;
      float hiwt2 = wt1;
      if(wt1 < lowt2) {
        lowt2 = wt1;
        hiwt2 = wt0;
      }
      // Check for this configuration
      //  loWire1.......hiWire1   minSep  loWire2....hiWire2
      //  loTime1.......hiTime1   minSep  loTime2....hiTime2
      if(lowt2 > hiwt1 + minSep) return false;
      // and the other
      if(lowt1 > hiwt2 + minSep) return false;
    } // iwt

    float best = minSep * minSep;
    ipt1 = 0; ipt2 = 0;
    float dwc = 0;
    bool isClose = false;
    for(unsigned short i1 = tj1.EndPt[0]; i1 < tj1.EndPt[1] + 1; ++i1) {
      for(unsigned short i2 = tj2.EndPt[0]; i2 < tj2.EndPt[1] + 1; ++i2) {
        if(considerDeadWires) dwc = DeadWireCount(tjs, tj1.Pts[i1], tj2.Pts[i2]);
        float dw = tj1.Pts[i1].Pos[0] - tj2.Pts[i2].Pos[0] - dwc;
        if(std::abs(dw) > minSep) continue;
        float dt = tj1.Pts[i1].Pos[1] - tj2.Pts[i2].Pos[1];
        if(std::abs(dt) > minSep) continue;
        float dp2 = dw * dw + dt * dt;
        if(dp2 < best) {
          best = dp2;
          ipt1 = i1;
          ipt2 = i2;
          isClose = true;
        }
      } // i2
    } // i1
    minSep = sqrt(best);
    return isClose;
  } // TrajTrajDOCA

  //////////////////////////////////////////
  float HitSep2(TjStuff& tjs, unsigned int iht, unsigned int jht)
  {
    // returns the separation^2 between two hits in WSE units
    if(iht > tjs.fHits.size()-1 || jht > tjs.fHits.size()-1) return 1E6;
    float dw = (float)tjs.fHits[iht].ArtPtr->WireID().Wire - (float)tjs.fHits[jht].ArtPtr->WireID().Wire;
    float dt = (tjs.fHits[iht].PeakTime - tjs.fHits[jht].PeakTime) * tjs.UnitsPerTick;
    return dw * dw + dt * dt;
  } // HitSep2
  
  //////////////////////////////////////////
  unsigned short CloseEnd(TjStuff& tjs, const Trajectory& tj, const Point2_t& pos)
  {
    unsigned short endPt = tj.EndPt[0];
    auto& tp0 = tj.Pts[endPt];
    endPt = tj.EndPt[1];
    auto& tp1 = tj.Pts[endPt];
    if(PosSep2(tp0.Pos, pos) < PosSep2(tp1.Pos, pos)) return 0;
    return 1;
  } // CloseEnd

  //////////////////////////////////////////
  float PointTrajSep2(float wire, float time, TrajPoint const& tp)
  {
    float dw = wire - tp.Pos[0];
    float dt = time - tp.Pos[1];
    return dw * dw + dt * dt;
  }
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff const& tjs, unsigned int iht, TrajPoint const& tp)
  {
    float wire = tjs.fHits[iht].ArtPtr->WireID().Wire;
    float time = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff const& tjs, float wire, float time, TrajPoint const& tp)
  {
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA2(TjStuff const& tjs, float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach squared between a (wire, time(WSE)) point
    // and a trajectory point
    
    float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
    float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    float dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return (dw * dw + dt * dt);
    
  } // PointTrajDOCA2
  
  //////////////////////////////////////////
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, Point2_t& pos)
  {
    TrajIntersection(tp1, tp2, pos[0], pos[1]);
  } // TrajIntersection
  //////////////////////////////////////////
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y)
  {
    // returns the intersection position, (x,y), of two trajectory points
    
    x = -9999; y = -9999;
    
    double arg1 = tp1.Pos[0] * tp1.Dir[1] - tp1.Pos[1] * tp1.Dir[0];
    double arg2 = tp2.Pos[0] * tp1.Dir[1] - tp2.Pos[1] * tp1.Dir[0];
    double arg3 = tp2.Dir[0] * tp1.Dir[1] - tp2.Dir[1] * tp1.Dir[0];
    if(arg3 == 0) return;
    double s = (arg1 - arg2) / arg3;
    
    x = (float)(tp2.Pos[0] + s * tp2.Dir[0]);
    y = (float)(tp2.Pos[1] + s * tp2.Dir[1]);
    
  } // TrajIntersection
  
  //////////////////////////////////////////
  float TrajLength(Trajectory& tj)
  {
    float len = 0, dx, dy;
    unsigned short ipt;
    unsigned short prevPt = tj.EndPt[0];
    for(ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1] + 1; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dx = tj.Pts[ipt].Pos[0] - tj.Pts[prevPt].Pos[0];
      dy = tj.Pts[ipt].Pos[1] - tj.Pts[prevPt].Pos[1];
      len += sqrt(dx * dx + dy * dy);
      prevPt = ipt;
    }
    return len;
  } // TrajLength

  //////////////////////////////////////////
  float PosSep(const Point2_t& pos1, const Point2_t& pos2)
  {
    return sqrt(PosSep2(pos1, pos2));
  } // PosSep
  
  //////////////////////////////////////////
  float PosSep2(const Point2_t& pos1, const Point2_t& pos2)
  {
    // returns the separation distance^2 between two positions
    float d0 = pos1[0] - pos2[0];
    float d1 = pos1[1] - pos2[1];
    return d0*d0+d1*d1;
  } // PosSep2
  
  //////////////////////////////////////////
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Returns the separation distance between two trajectory points
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation
  
  //////////////////////////////////////////
  bool TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& DOCA)
  {
    // find the closest approach between a trajectory tj and a point (x,y). Returns
    // the index of the closest trajectory point and the distance. Returns false if none
    // of the points on the tj are within DOCA
    
    float close2 = DOCA * DOCA;
    closePt = 0;
    bool foundClose = false;
    
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      float dx = tj.Pts[ipt].Pos[0] - x;
      if(std::abs(dx) > DOCA) continue;
      float dy = tj.Pts[ipt].Pos[1] - y;
      if(std::abs(dy) > DOCA) continue;
      float sep2 = dx * dx + dy * dy;
      if(sep2 < close2) {
        close2 = sep2;
        closePt = ipt;
        foundClose = true;
      }
    } // ipt
    
    DOCA = sqrt(close2);
    return foundClose;
    
  } // TrajClosestApproach
  
  /////////////////////////////////////////
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Calculates the angle of a line between two TPs
    float dw = tp2.Pos[0] - tp1.Pos[0];
    float dt = tp2.Pos[1] - tp1.Pos[1];
    return atan2(dw, dt);
  } // TwoTPAngle
  
  ////////////////////////////////////////////////
  std::vector<unsigned int> PutTrajHitsInVector(Trajectory const& tj, HitStatus_t hitRequest)
  {
    // Put hits in each trajectory point into a flat vector
    std::vector<unsigned int> hitVec;
    
    // special handling for shower trajectories. UseHit isn't valid
    if(tj.AlgMod[kShowerTj]) {
      for(auto& tp : tj.Pts) hitVec.insert(hitVec.end(), tp.Hits.begin(), tp.Hits.end());
      return hitVec;
    } // shower Tj
    
    // reserve under the assumption that there will be one hit per point
    hitVec.reserve(tj.Pts.size());
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        unsigned int iht = tj.Pts[ipt].Hits[ii];
        bool useit = (hitRequest == kAllHits);
        if(tj.Pts[ipt].UseHit[ii] && hitRequest == kUsedHits) useit = true;
        if(!tj.Pts[ipt].UseHit[ii] && hitRequest == kUnusedHits) useit = true;
        if(useit) hitVec.push_back(iht);
      } // iht
    } // ipt
    return hitVec;
  } // PutTrajHitsInVector
  
  //////////////////////////////////////////
  void TagJunkTj(TjStuff const& tjs, Trajectory& tj, bool prt)
  {
    // Characterizes the trajectory as a junk tj even though it may not
    // have been reconstructed in FindJunkTraj. The distinguishing feature is
    // that it is short and has many used hits in each trajectory point.
    
    // Don't bother if it is too long
    if(tj.Pts.size() > 10) return;
    // count the number of points that have many used hits
    unsigned short nhm = 0;
    unsigned short npwc = 0;
    for(auto& tp : tj.Pts) {
      if(tp.Chg == 0) continue;
      ++npwc;
      unsigned short nused = 0;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) ++nused;
      } // ii
      if(nused > 3) ++nhm;
    } // tp
    // Set the junkTj bit if most of the hits are used in most of the tps
    if(nhm > 0.5 * npwc) tj.AlgMod[kJunkTj] = true;
    if(prt) mf::LogVerbatim("TC")<<"TGT: "<<tj.ID<<" npwc "<<npwc<<" nhm "<<nhm<<" junk? "<<tj.AlgMod[kJunkTj];
  } // TagJunkTj

  //////////////////////////////////////////
  bool HasDuplicateHits(TjStuff const& tjs, Trajectory const& tj, bool prt)
  {
    // returns true if a hit is associated with more than one TP
    auto tjHits = PutTrajHitsInVector(tj, kAllHits);
    for(unsigned short ii = 0; ii < tjHits.size() - 1; ++ii) {
      for(unsigned short jj = ii + 1; jj < tjHits.size(); ++jj) {
        if(tjHits[ii] == tjHits[jj]) {
          if(prt) mf::LogVerbatim("TC")<<"HDH: Hit "<<PrintHit(tjs.fHits[ii])<<" is a duplicate "<<ii<<" "<<jj;
          return true;
        }
      } // jj
    } // ii
    return false;
  } // HasDuplicateHits
  
  //////////////////////////////////////////
  void MoveTPToWire(TrajPoint& tp, float wire)
  {
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    if(tp.Dir[0] == 0) return;
    float dw = wire - tp.Pos[0];
    if(std::abs(dw) < 0.01) return;
    tp.Pos[0] = wire;
    tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
  } // MoveTPToWire
  
  //////////////////////////////////////////
  std::vector<unsigned int> FindCloseHits(TjStuff const& tjs, std::array<int, 2> const& wireWindow, Point2_t const& timeWindow, const unsigned short plane, HitStatus_t hitRequest, bool usePeakTime, bool& hitsNear)
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
    if(plane > tjs.FirstWire.size() - 1) return closeHits;
    // window in the wire coordinate
    int loWire = wireWindow[0];
    if(loWire < (int)tjs.FirstWire[plane]) loWire = tjs.FirstWire[plane];
    int hiWire = wireWindow[1];
    if(hiWire > (int)tjs.LastWire[plane]-1) hiWire = tjs.LastWire[plane]-1;
    // window in the time coordinate
    float minTick = timeWindow[0] / tjs.UnitsPerTick;
    float maxTick = timeWindow[1] / tjs.UnitsPerTick;
    for(int wire = loWire; wire <= hiWire; ++wire) {
      // Set hitsNear if the wire is dead
      if(tjs.WireHitRange[plane][wire].first == -2) hitsNear = true;
      if(tjs.WireHitRange[plane][wire].first < 0) continue;
      unsigned int firstHit = (unsigned int)tjs.WireHitRange[plane][wire].first;
      unsigned int lastHit = (unsigned int)tjs.WireHitRange[plane][wire].second;
      for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
        if(usePeakTime) {
          if(tjs.fHits[iht].PeakTime < minTick) continue;
          if(tjs.fHits[iht].PeakTime > maxTick) break;
        } else {
          int hiLo = minTick;
          if(tjs.fHits[iht].StartTick > hiLo) hiLo = tjs.fHits[iht].StartTick;
          int loHi = maxTick;
          if(tjs.fHits[iht].EndTick < loHi) loHi = tjs.fHits[iht].EndTick;
          if(loHi < hiLo) continue;
          if(hiLo > loHi) break;
        }
        hitsNear = true;
        bool takeit = (hitRequest == kAllHits);
        if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) takeit = true;
        if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) takeit = true;
        if(takeit) closeHits.push_back(iht);
      } // iht
    } // wire
    return closeHits;
  } // FindCloseHits
  
  //////////////////////////////////////////
  bool FindCloseHits(TjStuff const& tjs, TrajPoint& tp, float const& maxDelta, HitStatus_t hitRequest)
  {
    // Fills tp.Hits sets tp.UseHit true for hits that are close to tp.Pos. Returns true if there are
    // close hits OR if the wire at this position is dead
    
    tp.Hits.clear();
    tp.UseHit.reset();
    if(!WireHitRangeOK(tjs, tp.CTP)) {
      std::cout<<"FindCloseHits: WireHitRange not valid for CTP "<<tp.CTP<<". tjs.WireHitRange Cstat "<<tjs.TPCID.Cryostat<<" TPC "<<tjs.TPCID.TPC<<"\n";
      return false;
    }
    
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short ipl = planeID.Plane;
    
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if(wire < tjs.FirstWire[ipl]) return false;
    if(wire > tjs.LastWire[ipl]-1) return false;
    
    // dead wire
    if(tjs.WireHitRange[ipl][wire].first == -1) return true;
    // live wire with no hits
    if(tjs.WireHitRange[ipl][wire].first == -2) return false;
    
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;

    float fwire = wire;
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) useit = true;
      if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) useit = true;
      if(!useit) continue;
      float ftime = tjs.UnitsPerTick * tjs.fHits[iht].PeakTime;
      float delta = PointTrajDOCA(tjs, fwire, ftime, tp);
      if(delta < maxDelta) tp.Hits.push_back(iht);
    } // iht
    if(tp.Hits.size() > 16) {
      tp.Hits.resize(16);
    }
    // Set UseHit false. The calling routine should decide if these hits should be used
    tp.UseHit.reset();
    return true;
    
  } // FindCloseHits
  
  //////////////////////////////////////////
  std::vector<int> FindCloseTjs(const TjStuff& tjs, const TrajPoint& fromTp, const TrajPoint& toTp, const float& maxDelta)
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
    
    TrajPoint tp;
    // Make the tp so that stepping is positive
    unsigned int firstWire, lastWire;
    if(toTp.Pos[0] > fromTp.Pos[0]) {
      if(!MakeBareTrajPoint(tjs, fromTp, toTp, tp)) return tmp;
      firstWire = std::nearbyint(fromTp.Pos[0]);
      lastWire = std::nearbyint(toTp.Pos[0]);
    } else if(toTp.Pos[0] < fromTp.Pos[0]) {
      if(!MakeBareTrajPoint(tjs, toTp, fromTp, tp)) return tmp;
      firstWire = std::nearbyint(toTp.Pos[0]);
      lastWire = std::nearbyint(fromTp.Pos[0]);
    } else {
      tp.Pos = fromTp.Pos;
      float tmp = fromTp.Pos[0] - maxDelta;
      if(tmp < 0) tmp = 0;
      firstWire = std::nearbyint(tmp);
      tmp = fromTp.Pos[0] + maxDelta;
      lastWire = std::nearbyint(tmp);
    }
    
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short ipl = planeID.Plane;
    
    if(firstWire < tjs.FirstWire[ipl]) firstWire = tjs.FirstWire[ipl];
    if(firstWire > tjs.LastWire[ipl]-1) return tmp;
    if(lastWire < tjs.FirstWire[ipl]) return tmp;
    if(lastWire > tjs.LastWire[ipl]-1) lastWire = tjs.LastWire[ipl]-1;
    
    for(unsigned int wire = firstWire; wire <= lastWire; ++wire) {
      if(tjs.WireHitRange[ipl][wire].first == -1) continue;
      if(tjs.WireHitRange[ipl][wire].first == -2) continue;
      MoveTPToWire(tp, (float)wire);
      // Find the tick range at this position
      float minTick = (tp.Pos[1] - maxDelta) / tjs.UnitsPerTick;
      float maxTick = (tp.Pos[1] + maxDelta) / tjs.UnitsPerTick;
      unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
      unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
      for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
        if(tjs.fHits[iht].InTraj <= 0) continue;
        if(tjs.fHits[iht].PeakTime < minTick) continue;
        // Hits are sorted by increasing time so we can break when maxTick is reached
        if(tjs.fHits[iht].PeakTime > maxTick) break;
        if(std::find(tmp.begin(), tmp.end(), tjs.fHits[iht].InTraj) != tmp.end()) continue;
        tmp.push_back(tjs.fHits[iht].InTraj);
      } // iht
    } // wire
    
    return tmp;
    
  } // FindCloseTjs
  
  ////////////////////////////////////////////////
  float ChgFracNearPos(TjStuff& tjs, const Point2_t& pos, const std::vector<int>& tjIDs)
  {
    // returns the fraction of the charge in the region around pos that is associated with
    // the list of Tj IDs
    if(tjIDs.empty()) return 0;
    std::array<int, 2> wireWindow;
    Point2_t timeWindow;
    // 1/2 size of the region
    constexpr float NNDelta = 5; 
    wireWindow[0] = pos[0] - NNDelta;
    wireWindow[1] = pos[0] + NNDelta;
    timeWindow[0] = pos[1] - NNDelta;
    timeWindow[1] = pos[1] + NNDelta;
    // do some checking
    for(auto& tjID : tjIDs) if(tjID <= 0 || tjID > (int)tjs.allTraj.size()) return 0;
    // Determine which plane we are in
    geo::PlaneID planeID = DecodeCTP(tjs.allTraj[tjIDs[0]-1].CTP);
    // get a list of all hits in this region
    bool hitsNear;
    std::vector<unsigned int> closeHits = FindCloseHits(tjs, wireWindow, timeWindow, planeID.Plane, kAllHits, true, hitsNear);
    if(closeHits.empty()) return 0;
    float chg = 0;
    float tchg = 0;
    // Add the hit charge in the box
    // All hits in the box, and all hits associated with the Tjs
    for(auto& iht : closeHits) {
      chg += tjs.fHits[iht].Integral;
      if(tjs.fHits[iht].InTraj <= 0) continue;
      if(std::find(tjIDs.begin(), tjIDs.end(), tjs.fHits[iht].InTraj) != tjIDs.end()) tchg += tjs.fHits[iht].Integral;
    } // iht
    if(chg == 0) return 0;
    return tchg / chg;
  } // ChgFracNearPos
  
  ////////////////////////////////////////////////
  float MaxHitDelta(TjStuff& tjs, Trajectory& tj)
  {
    float delta, md = 0;
    unsigned short ii;
    unsigned int iht;
    for(auto& tp : tj.Pts) {
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        iht = tp.Hits[ii];
        delta = PointTrajDOCA(tjs, iht, tp);
        if(delta > md) md = delta;
      } // ii
    } // pts
    return md;
  } // MaxHitDelta

  //////////////////////////////////////////
  void ReverseTraj(TjStuff& tjs, Trajectory& tj)
  {
    // reverse the trajectory
    if(tj.Pts.empty()) return;
    if(tj.AlgMod[kMat3D]) {
      std::cout<<"Trying to reverse a 3D matched Tj. Need to modify other Tjs and the MatchStruct\n";
      return;
    }
    if(tj.AlgMod[kSetDir]) {
      std::cout<<"Trying to reverse Tj "<<tj.ID<<" whose direction has been set. Not doing it.\n";
      return;
    }
    // reverse the crawling direction flag
    tj.StepDir = -tj.StepDir;
    // Vertices
    std::swap(tj.VtxID[0], tj.VtxID[1]);
    // trajectory points
    std::reverse(tj.Pts.begin(), tj.Pts.end());
    // reverse the stop flag
    std::reverse(tj.StopFlag.begin(), tj.StopFlag.end());
    std::swap(tj.dEdx[0], tj.dEdx[1]);
    // reverse the direction vector on all points
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Dir[0] != 0) tj.Pts[ipt].Dir[0] = -tj.Pts[ipt].Dir[0];
      if(tj.Pts[ipt].Dir[1] != 0) tj.Pts[ipt].Dir[1] = -tj.Pts[ipt].Dir[1];
      if(tj.Pts[ipt].Ang > 0) {
        tj.Pts[ipt].Ang -= M_PI;
      } else {
        tj.Pts[ipt].Ang += M_PI;
      }
//      tj.Pts[ipt].Ang = std::atan2(tj.Pts[ipt].Dir[1], tj.Pts[ipt].Dir[0]);
    } // ipt
    SetEndPoints(tjs, tj);
    // correct mallTraj if it exists
    if(tjs.mallTraj.empty()) return;
    for(auto& tj2pt : tjs.mallTraj) {
      if(tj2pt.id != tj.ID) continue;
      tj2pt.ipt = tj.Pts.size() - tj2pt.ipt - 1;
    } // tj2pt
  } // ReverseTraj
  
  //////////////////////////////////////////
  bool PointInsideEnvelope(const Point2_t& Point, const std::vector<Point2_t>& Envelope)
  {
    // returns true if the Point is within the Envelope polygon. Entries in Envelope are the
    // Pos[0], Pos[1] locations of the polygon vertices. This is based on the algorithm that the
    // sum of the angles of a vector between a point and the vertices will be 2 * pi for an interior
    // point and 0 for an exterior point
    
    Point2_t p1, p2;
    unsigned short nvx = Envelope.size();
    double angleSum = 0;
    for(unsigned short ii = 0; ii < Envelope.size(); ++ii) {
      p1[0] = Envelope[ii][0] - Point[0];
      p1[1] = Envelope[ii][1] - Point[1];
      p2[0] = Envelope[(ii+1)%nvx][0] - Point[0];
      p2[1] = Envelope[(ii+1)%nvx][1] - Point[1];
      angleSum += DeltaAngle(p1, p2);
    }
    if(abs(angleSum) < M_PI) return false;
    return true;
      
  } // InsideEnvelope

  //////////////////////////////////////////
  double DeltaAngle(const Point2_t& p1, const Point2_t& p2)
  {
    // angle between two points
    double ang1 = atan2(p1[1], p1[0]);
    double ang2 = atan2(p2[1], p2[0]);
    return DeltaAngle2(ang1, ang2);
  } // DeltaAngle
  
  //////////////////////////////////////////
  double DeltaAngle2(double Ang1, double Ang2)
  {
    constexpr double twopi = 2 * M_PI;
    double dang = Ang1 - Ang2;
    while(dang >  M_PI) dang -= twopi;
    while(dang < -M_PI) dang += twopi;
    return dang;
  }

  //////////////////////////////////////////
  double DeltaAngle(double Ang1, double Ang2) 
  {
    return std::abs(std::remainder(Ang1 - Ang2, M_PI));
  }
  
  ////////////////////////////////////////////////
  void SetEndPoints(TjStuff& tjs, Trajectory& tj)
  {
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
    
    // don't mess with showerTjs
    if(tj.AlgMod[kShowerTj]) return;
    tj.EndPt[0] = 0; tj.EndPt[1] = 0;
    if(tj.Pts.size() == 0) return;
    
    // check the end point pointers
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[0] = ipt;
        break;
      }
    }
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[1] = ipt;
        break;
      }
    }
  } // SetEndPoints
  
  ////////////////////////////////////////////////
  bool TrajIsClean(TjStuff& tjs, Trajectory& tj, bool prt)
  {
    // Returns true if the trajectory has low hit multiplicity and is in a
    // clean environment
    unsigned short nUsed = 0;
    unsigned short nTotHits = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      nTotHits += tp.Hits.size();
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) ++nUsed;
      } // ii
    } // ipt
    if(nTotHits == 0) return false;
    float fracUsed = (float)nUsed / (float)nTotHits;
    if(prt) mf::LogVerbatim("TC")<<"TrajIsClean: nTotHits "<<nTotHits<<" nUsed "<<nUsed<<" fracUsed "<<fracUsed;
    
    if(fracUsed > 0.9) return true;
    return false;
    
  } // TrajIsClean
  
  ////////////////////////////////////////////////
  short MCSMom(TjStuff& tjs, Trajectory& tj)
  {
    return MCSMom(tjs, tj, tj.EndPt[0], tj.EndPt[1]);
  } // MCSMom
  
  
  ////////////////////////////////////////////////
  short MCSMom(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt)
  {
    // Estimate the trajectory momentum using Multiple Coulomb Scattering ala PDG RPP
    
    firstPt = NearestPtWithChg(tjs, tj, firstPt);
    lastPt = NearestPtWithChg(tjs, tj, lastPt);
    if(firstPt >= lastPt) return 0;
    
    if(firstPt < tj.EndPt[0]) return 0;
    if(lastPt > tj.EndPt[1]) return 0;
    // Can't do this with only 2 points
    if(NumPtsWithCharge(tjs, tj, false, firstPt, lastPt) < 3) return 0;
    // Ignore junk Tjs
    if(tj.AlgMod[kJunkTj]) return 0;
        
    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if(tjLen == 0) return 0;
    // mom calculated in MeV
    double mom = 13.8 * sqrt(tjLen / 14) / MCSThetaRMS(tjs, tj, firstPt, lastPt);
    if(mom > 999) mom = 999;
    return (short)mom;
  } // MCSMom
    
  
  ////////////////////////////////////////////////
  unsigned short NearestPtWithChg(TjStuff& tjs, Trajectory& tj, unsigned short thePt)
  {
    // returns a point near thePt which has charge
    if(thePt > tj.EndPt[1]) return thePt;
    if(tj.Pts[thePt].Chg > 0) return thePt;
    
    short endPt0 = tj.EndPt[0];
    short endPt1 = tj.EndPt[1];
    for(short off = 1; off < 10; ++off) {
      short ipt = thePt + off;
      if(ipt <= endPt1 && tj.Pts[ipt].Chg > 0) return (unsigned short)ipt;
      ipt = thePt - off;
      if(ipt >= endPt0 && tj.Pts[ipt].Chg > 0) return (unsigned short)ipt;
    } // off
    return thePt;
  } // NearestPtWithChg
  
  /////////////////////////////////////////
  float MCSThetaRMS(TjStuff& tjs, Trajectory& tj)
  {
    // This returns the MCS scattering angle expected for one WSE unit of travel along the trajectory.
    // It is used to define kink and vertex cuts. This should probably be named something different to
    // prevent confusion
    
    float tps = TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
    if(tps == 0) return 1;
    
    return MCSThetaRMS(tjs, tj, tj.EndPt[0], tj.EndPt[1]) / sqrt(tps);
    
  } // MCSThetaRMS
  
  /////////////////////////////////////////
  double MCSThetaRMS(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt)
  {
    // This returns the MCS scattering angle expected for the length of the trajectory
    // spanned by firstPt to lastPt. It is used primarily to calculate MCSMom
    
    if(firstPt < tj.EndPt[0]) return 1;
    if(lastPt > tj.EndPt[1]) return 1;
    
    firstPt = NearestPtWithChg(tjs, tj, firstPt);
    lastPt = NearestPtWithChg(tjs, tj, lastPt);
    if(firstPt >= lastPt) return 1;
    
    double sigmaS;
    unsigned short cnt;
    TjDeltaRMS(tjs, tj, firstPt, lastPt, sigmaS, cnt);
    if(sigmaS < 0) return 1;
    // require that cnt is a significant fraction of the total number of charged points
    // so that we don't get erroneously high MCSMom when there are large gaps.
    // This is the number of points expected in the count if there are no gaps
    unsigned short numPts = lastPt - firstPt - 1;
    // return the previously calculated value of MCSMom
    if(numPts > 5 && cnt < 0.7 * numPts) return tj.MCSMom;
    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if(tjLen == 0) return 1;
    // Theta_o =  4 * sqrt(3) * sigmaS / path
    return (6.8 * sigmaS / tjLen);
    
  } // MCSThetaRMS
  
  /////////////////////////////////////////
  void TjDeltaRMS(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt, double& rms, unsigned short& cnt)
  {
    // returns the rms scatter of points around a line formed by the firstPt and lastPt of the trajectory
    
    rms = -1;
    if(firstPt < tj.EndPt[0]) return;
    if(lastPt > tj.EndPt[1]) return;
    
    firstPt = NearestPtWithChg(tjs, tj, firstPt);
    lastPt = NearestPtWithChg(tjs, tj, lastPt);
    if(firstPt >= lastPt) return;
    
    TrajPoint tmp;
    // make a bare trajectory point to define a line between firstPt and lastPt.
    // Use the position of the hits at these points
    TrajPoint firstTP = tj.Pts[firstPt];
    firstTP.Pos = firstTP.HitPos;
    TrajPoint lastTP = tj.Pts[lastPt];
    lastTP.Pos = lastTP.HitPos;
    if(!MakeBareTrajPoint(tjs, firstTP, lastTP, tmp)) return;
    // sum up the deviations^2
    double dsum = 0;
    cnt = 0;
    for(unsigned short ipt = firstPt + 1; ipt < lastPt; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dsum += PointTrajDOCA2(tjs, tj.Pts[ipt].HitPos[0],  tj.Pts[ipt].HitPos[1], tmp);
      ++cnt;
    } // ipt
    if(cnt < 2) return;
    rms = sqrt(dsum / (double)cnt);

  } // TjDeltaRMS

  /////////////////////////////////////////
  void TagDeltaRays(TjStuff& tjs, const CTP_t& inCTP)
  {
    // DeltaRayTag vector elements
    // [0] = max separation of both endpoints from a muon
    // [1] = minimum MCSMom
    // [2] = maximum MCSMom
    
    if(!tjs.UseAlg[kDeltaRay]) return;
    if(tjs.DeltaRayTag[0] < 0) return;
    if(tjs.DeltaRayTag.size() < 3) return;
    
    bool prt = (debug.CTP == inCTP) && (debug.Tick == 31313);

    // double the user-defined separation cut. We will require that at least one of the ends of 
    // a delta ray be within the user-defined cut and allow
    float maxSep = 2 * tjs.DeltaRayTag[0];
    float maxMinSep = tjs.DeltaRayTag[0];
    unsigned short minMom = tjs.DeltaRayTag[1];
    unsigned short maxMom = tjs.DeltaRayTag[2];
    unsigned short minMuonLength = 2 * tjs.Vertex2DCuts[2];
    unsigned short minpts = 4;
    if(prt) mf::LogVerbatim("TC")<<"TagDeltaRays: maxSep "<<maxSep<<" maxMinSep "<<maxMinSep<<" Mom range "<<minMom<<" to "<<maxMom<<" minpts "<<minpts;
    
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& muTj = tjs.allTraj[itj];
      if(muTj.CTP != inCTP) continue;
      if(muTj.AlgMod[kKilled]) continue;
      if(muTj.PDGCode != 13) continue;
      if(prt) mf::LogVerbatim("TC")<<"TagDeltaRays: Muon "<<muTj.ID<<" EndPts "<<PrintPos(tjs, muTj.Pts[muTj.EndPt[0]])<<"-"<<PrintPos(tjs, muTj.Pts[muTj.EndPt[1]]);
      // min length
      if(muTj.EndPt[1] - muTj.EndPt[0] < minMuonLength) continue;
      auto& mtp0 = muTj.Pts[muTj.EndPt[0]];
      auto& mtp1 = muTj.Pts[muTj.EndPt[1]];
      // Found a muon, now look for delta rays
      for(unsigned short jtj = 0; jtj < tjs.allTraj.size(); ++jtj) {
        if(jtj == itj) continue;
        Trajectory& dtj = tjs.allTraj[jtj];
        if(dtj.AlgMod[kKilled]) continue;
        if(dtj.CTP != inCTP) continue;
        if(dtj.PDGCode == 13) continue;
        // MCSMom cut
        if(dtj.MCSMom < minMom) continue;
        if(dtj.MCSMom > maxMom) continue;
        if(dtj.EndPt[1] - dtj.EndPt[0] < minpts) continue;
        // some rough cuts to require that the delta ray is within the
        // ends of the muon
        auto& dtp0 = dtj.Pts[dtj.EndPt[0]];
        auto& dtp1 = dtj.Pts[dtj.EndPt[1]];
        if(muTj.StepDir > 0) {
          if(dtp0.Pos[0] < mtp0.Pos[0]) continue;
          if(dtp1.Pos[0] > mtp1.Pos[0]) continue;
        } else {
          if(dtp0.Pos[0] > mtp0.Pos[0]) continue;
          if(dtp1.Pos[0] < mtp1.Pos[0]) continue;
        }
        // find the minimum separation
        float doca = maxMinSep;
        unsigned short mpt = 0;
        unsigned short dpt = 0;
        TrajTrajDOCA(tjs, muTj, dtj, mpt, dpt, doca, false);
        if(doca == maxMinSep) continue;
        auto& dTp = dtj.Pts[dpt];
        // cut on the distance from the muon ends
        if(PosSep(dTp.Pos, mtp0.Pos) < tjs.Vertex2DCuts[2]) continue;
        if(PosSep(dTp.Pos, mtp1.Pos) < tjs.Vertex2DCuts[2]) continue;
       // make an angle cut at this point. A delta-ray should have a small angle
        float dang = DeltaAngle(muTj.Pts[mpt].Ang, dtj.Pts[dpt].Ang);
        if(prt) mf::LogVerbatim("TC")<<" dRay? "<<dtj.ID<<" at "<<PrintPos(tjs, dtj.Pts[dpt].Pos)<<" dang "<<dang;
        if(dang > tjs.KinkCuts[0]) continue;
        unsigned short oend = 0;
        // check the delta at the end of the delta-ray that is farthest away from the
        // closest point
        if(dpt > 0.5 * (dtj.EndPt[0] + dtj.EndPt[1])) oend = 1;
        auto& farEndTP = dtj.Pts[dtj.EndPt[oend]];
        float farEndDelta = PointTrajDOCA(tjs, farEndTP.Pos[0], farEndTP.Pos[1], muTj.Pts[mpt]);
        if(prt) mf::LogVerbatim("TC")<<"  farEnd "<<PrintPos(tjs, farEndTP.Pos)<<" farEndDelta "<<farEndDelta;
        if(farEndDelta > maxSep) continue;
        if(prt) mf::LogVerbatim("TC")<<"   delta ray "<<dtj.ID<<" parent -> "<<muTj.ID;
        dtj.ParentID = muTj.ID;
        dtj.PDGCode = 11;
        dtj.AlgMod[kDeltaRay] = true;
        // Set the start of the delta-ray to be end 0
        if(oend != 1) ReverseTraj(tjs, dtj);
        dtj.AlgMod[kSetDir] = true;
      } // jtj
    } // itj
    
  } // TagDeltaRays
  
  /////////////////////////////////////////
  void TagMuonDirections(TjStuff& tjs, short debugWorkID)
  {
    // Determine muon directions delta-ray proximity to muon trajectories
    
    if(tjs.MuonTag[0] < 0) return;
    
    unsigned short minLen = tjs.MuonTag[3];
    
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& muTj = tjs.allTraj[itj];
      if(muTj.AlgMod[kKilled]) continue;
      bool prt = (debugWorkID < 0 && muTj.WorkID == debugWorkID);
      if(prt) {
        mf::LogVerbatim("TC")<<"TagMuonDirection: Muon "<<muTj.CTP<<" "<<PrintPos(tjs, muTj.Pts[muTj.EndPt[0]])<<"-"<<PrintPos(tjs, muTj.Pts[muTj.EndPt[1]]);
      }
      if(muTj.PDGCode != 13) continue;
      // look for delta ray trajectories and count the number of times that
      unsigned short nPos = 0;
      unsigned short nNeg = 0;
      for(auto& dtj : tjs.allTraj) {
        if(dtj.AlgMod[kKilled]) continue;
        if(dtj.ParentID != muTj.ID) continue;
        if(dtj.EndPt[1] - dtj.EndPt[0] > minLen) continue;
        if(dtj.StepDir > 0) {
          ++nPos;
        } else {
          ++nNeg;
        }
      } // dtj
      if(nPos == nNeg) continue;
      if(nPos > nNeg) {
        if(muTj.StepDir < 0) ReverseTraj(tjs, muTj);
      } else {
        if(muTj.StepDir > 0) ReverseTraj(tjs, muTj);
      }
      muTj.AlgMod[kSetDir] = true;
    } // itj
  } // TagMuonDirections
  
  /////////////////////////////////////////
  TrajPoint MakeBareTrajPoint(TjStuff& tjs, Point3_t& pos, Vector3_t& dir, CTP_t inCTP)
  {
    // Projects the space point defined by pos and dir into the CTP and returns it in the form of a trajectory point.
    // The TP Pos[0] is set to a negative number if the point has an invalid wire position but doesn't return an
    // error if the position is on a dead wire
    TrajPoint tp;
    tp.CTP = inCTP;
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    tp.Pos[0] = tjs.geom->WireCoordinate(pos[1], pos[2], planeID);
    if(tp.Pos[0] < 0 || tp.Pos[0] > tjs.MaxPos0[planeID.Plane]) {
      tp.Pos[0] = -1;
      return tp;
    }
    tp.Pos[1] = tjs.detprop->ConvertXToTicks(pos[0], planeID) * tjs.UnitsPerTick;
    // now find the direction
    // Make a point at the origin and one 100 units away
    Point3_t ori3 = {0, 0, 0};
    Point3_t pos3 = {100 * dir[0], 100 * dir[1], 100 * dir[2]};
    // 2D position of ori3 and the pos3 projection
    std::array<double, 2> ori2;
    std::array<double, 2> pos2;
    std::array<double, 2> dir2;
    // the wire coordinates
    ori2[0] = tjs.geom->WireCoordinate(ori3[1], ori3[2], planeID);
    pos2[0] = tjs.geom->WireCoordinate(pos3[1], pos3[2], planeID);
    // the time coordinates
    ori2[1] = tjs.detprop->ConvertXToTicks(ori3[0], planeID) * tjs.UnitsPerTick;
    pos2[1] = tjs.detprop->ConvertXToTicks(pos3[0], planeID) * tjs.UnitsPerTick;
    
    dir2[0] = pos2[0] - ori2[0];
    dir2[1] = pos2[1] - ori2[1];
    
    double norm = sqrt(dir2[0] * dir2[0] + dir2[1] * dir2[1]);
    tp.Dir[0] = dir2[0] / norm;
    tp.Dir[1] = dir2[1] / norm;
    tp.Ang = atan2(dir2[1], dir2[0]);
    return tp;
  } // MakeBareTP

  /////////////////////////////////////////
  bool MakeBareTrajPoint(const TjStuff& tjs, unsigned int fromHit, unsigned int toHit, TrajPoint& tp)
  {
    CTP_t tCTP = EncodeCTP(tjs.fHits[fromHit].ArtPtr->WireID());
    return MakeBareTrajPoint(tjs, (float)tjs.fHits[fromHit].ArtPtr->WireID().Wire, tjs.fHits[fromHit].PeakTime,
                                  (float)tjs.fHits[toHit].ArtPtr->WireID().Wire,   tjs.fHits[toHit].PeakTime, tCTP, tp);
    
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(const TjStuff& tjs, float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp)
  {
    tp.CTP = tCTP;
    tp.Pos[0] = fromWire;
    tp.Pos[1] = tjs.UnitsPerTick * fromTick;
    tp.Dir[0] = toWire - fromWire;
    tp.Dir[1] = tjs.UnitsPerTick * (toTick - fromTick);
    double norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    if(norm == 0) return false;
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    return true;
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(const Point2_t& fromPos, const Point2_t& toPos, TrajPoint& tpOut)
  {
    tpOut.Pos = fromPos;
    tpOut.Dir[0] = toPos[0] - fromPos[0];
    tpOut.Dir[1] = toPos[1] - fromPos[1];
    double norm = sqrt(tpOut.Dir[0] * tpOut.Dir[0] + tpOut.Dir[1] * tpOut.Dir[1]);
    if(norm == 0) return false;
    tpOut.Dir[0] /= norm;
    tpOut.Dir[1] /= norm;
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
    return true;
    
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(const TjStuff& tjs, const TrajPoint& tpIn1, const TrajPoint& tpIn2, TrajPoint& tpOut)
  {
    tpOut.CTP = tpIn1.CTP;
    tpOut.Pos = tpIn1.Pos;
    tpOut.Dir[0] = tpIn2.Pos[0] - tpIn1.Pos[0];
    tpOut.Dir[1] = tpIn2.Pos[1] - tpIn1.Pos[1];
    double norm = sqrt(tpOut.Dir[0] * tpOut.Dir[0] + tpOut.Dir[1] * tpOut.Dir[1]);
    if(norm == 0) return false;
    tpOut.Dir[0] /= norm;
    tpOut.Dir[1] /= norm;
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
    return true;
  } // MakeBareTrajPoint
  
  ////////////////////////////////////////////////
  float TPHitsRMSTime(TjStuff& tjs, TrajPoint& tp, HitStatus_t hitRequest)
  {
    return tjs.UnitsPerTick * TPHitsRMSTick(tjs, tp, hitRequest);
  } // TPHitsRMSTime

  ////////////////////////////////////////////////
  float TPHitsRMSTick(TjStuff& tjs, TrajPoint& tp, HitStatus_t hitRequest)
  {
    // Estimate the RMS of all hits associated with a trajectory point
    // without a lot of calculation
    if(tp.Hits.empty()) return 0;
    float minVal = 9999;
    float maxVal = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tp.UseHit[ii]) useit = true;
      if(hitRequest == kUnusedHits && !tp.UseHit[ii]) useit = true;
      if(!useit) continue;
      unsigned int iht = tp.Hits[ii];
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - rms;
      if(arg < minVal) minVal = arg;
      arg = cv + rms;
      if(arg > maxVal) maxVal = arg;
    } // ii
    if(maxVal == 0) return 0;
    return (maxVal - minVal) / 2;
  } // TPHitsRMSTick
  
  ////////////////////////////////////////////////
  float HitsRMSTime(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest)
  {
    return tjs.UnitsPerTick * HitsRMSTick(tjs, hitsInMultiplet, hitRequest);
  } // HitsRMSTick

  ////////////////////////////////////////////////
  float HitsRMSTick(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest)
  {
    if(hitsInMultiplet.empty()) return 0;
    
    if(hitsInMultiplet.size() == 1) return tjs.fHits[hitsInMultiplet[0]].RMS;
 
    float minVal = 9999;
    float maxVal = 0;
    for(unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
      unsigned int iht = hitsInMultiplet[ii];
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) useit = true;
      if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) useit = true;
      if(!useit) continue;
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - rms;
      if(arg < minVal) minVal = arg;
      arg = cv + rms;
      if(arg > maxVal) maxVal = arg;
    } // ii
    if(maxVal == 0) return 0;
    return (maxVal - minVal) / 2;
  } // HitsRMSTick
  
  ////////////////////////////////////////////////
  float HitsPosTime(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, float& sum, HitStatus_t hitRequest)
  {
    return tjs.UnitsPerTick * HitsPosTick(tjs, hitsInMultiplet, sum, hitRequest);
  } // HitsPosTime
  
  ////////////////////////////////////////////////
  float HitsPosTick(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, float& sum, HitStatus_t hitRequest)
  {
    // returns the position and the charge
    float pos = 0;
    sum = 0;
    for(unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
      unsigned int iht = hitsInMultiplet[ii];
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) useit = true;
      if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) useit = true;
      if(!useit) continue;
      float chg = tjs.fHits[iht].Integral;
      pos += chg * tjs.fHits[iht].PeakTime;
      sum += chg;
    } // ii
    if(sum == 0) return 0;
    return pos / sum;
  } // HitsPosTick
  
  //////////////////////////////////////////
  unsigned short NumUsedHitsInTj(const TjStuff& tjs, const Trajectory& tj)
  {
    if(tj.AlgMod[kKilled]) return 0;
    if(tj.Pts.empty()) return 0;
    unsigned short nhits = 0;
    for(auto& tp : tj.Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if(tp.UseHit[ii]) ++nhits;
    } // tp
    return nhits;
  } // NumHitsInTj
  
  //////////////////////////////////////////
  unsigned short NumHitsInTP(const TrajPoint& tp, HitStatus_t hitRequest)
  {
    // Counts the number of hits of the specified type in tp
    if(tp.Hits.empty()) return 0;
    
    if(hitRequest == kAllHits) return tp.Hits.size();
    
    unsigned short nhits = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(hitRequest == kUsedHits) {
        if(tp.UseHit[ii]) ++nhits;
      } else {
        // looking for unused hits
        if(!tp.UseHit[ii]) ++nhits;
      }
    } // ii
    return nhits;
  } // NumHitsInTP
  
  ////////////////////////////////////////////////
  void SetPDGCode(TjStuff& tjs, unsigned short itj)
  {
    if(itj > tjs.allTraj.size() - 1) return;
    SetPDGCode(tjs, tjs.allTraj[itj]);
  }
  
  ////////////////////////////////////////////////
  void SetPDGCode(TjStuff& tjs, Trajectory& tj)
  {
    // Sets the PDG code for the supplied trajectory. Note that the existing
    // PDG code is left unchanged if these cuts are not met
    
    tj.MCSMom = MCSMom(tjs, tj);
    if(tjs.MuonTag[0] <= 0) return;
    // Special handling of very long straight trajectories, e.g. uB cosmic rays
    unsigned short npts = tj.EndPt[1] - tj.EndPt[0];
    bool isAMuon = (npts > (unsigned short)tjs.MuonTag[0] && tj.MCSMom > tjs.MuonTag[1]);
    // anything really really long must be a muon
    if(npts > 500) isAMuon = true;
    if(isAMuon) tj.PDGCode = 13;
    
  } // SetPDGCode
  
  
  ////////////////////////////////////////////////
  bool FillWireHitRange(TjStuff& tjs, const geo::TPCID& tpcid, bool debugMode)
  {
    // fills the WireHitRange vector. Slightly modified version of the one in ClusterCrawlerAlg.
    // Returns false if there was a serious error
    
    // determine the number of planes
    geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    unsigned short nplanes = TPC.Nplanes();
    tjs.NumPlanes = nplanes;
    tjs.TPCID = tpcid;
    
    // Y,Z limits of the detector
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &thetpc = tjs.geom->TPC(tpc, cstat);
    thetpc.LocalToWorld(local,world);
    // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
    tjs.XLo = world[0]-tjs.geom->DetHalfWidth(tpc,cstat) + 1;
    tjs.XHi = world[0]+tjs.geom->DetHalfWidth(tpc,cstat) - 1;
    tjs.YLo = world[1]-tjs.geom->DetHalfHeight(tpc,cstat) + 1;
    tjs.YHi = world[1]+tjs.geom->DetHalfHeight(tpc,cstat) - 1;
    tjs.ZLo = world[2]-tjs.geom->DetLength(tpc,cstat)/2 + 1;
    tjs.ZHi = world[2]+tjs.geom->DetLength(tpc,cstat)/2 - 1;
    
    lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    if(!tjs.WireHitRange.empty()) tjs.WireHitRange.clear();
    
    // initialize everything
    tjs.WireHitRange.resize(nplanes);
    tjs.FirstWire.resize(nplanes);
    tjs.LastWire.resize(nplanes);
    tjs.NumWires.resize(nplanes);
    tjs.MaxPos0.resize(nplanes);
    tjs.MaxPos1.resize(nplanes);
    tjs.AveHitRMS.resize(nplanes, nplanes);
    
    std::pair<int, int> flag;
    flag.first = -2; flag.second = -2;
    
    // Calculate tjs.UnitsPerTick, the scale factor to convert a tick into
    // Wire Spacing Equivalent (WSE) units where the wire spacing in this plane = 1.
    // Strictly speaking this factor should be calculated for each plane to handle the
    // case where the wire spacing is different in each plane. Deal with this later if
    // the approximation used here fails.
    
    raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel(0, 0, (int)tpc, (int)cstat);
    float wirePitch = tjs.geom->WirePitch(tjs.geom->View(channel));
    float tickToDist = tjs.detprop->DriftVelocity(tjs.detprop->Efield(),tjs.detprop->Temperature());
    tickToDist *= 1.e-3 * tjs.detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    tjs.UnitsPerTick = tickToDist / wirePitch;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      tjs.FirstWire[ipl] = INT_MAX;
      tjs.LastWire[ipl] = 0;
      tjs.NumWires[ipl] = tjs.geom->Nwires(ipl, tpc, cstat);
      tjs.WireHitRange[ipl].resize(tjs.NumWires[ipl], flag);
      tjs.MaxPos0[ipl] = (float)tjs.NumWires[ipl] - 0.5;
      tjs.MaxPos1[ipl] = (float)tjs.detprop->NumberTimeSamples() * tjs.UnitsPerTick;
    }
    
    // overwrite with the "dead wires" condition
    flag.first = -1; flag.second = -1;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        raw::ChannelID_t chan = tjs.geom->PlaneWireToChannel((int)ipl, (int)wire, (int)tpc, (int)cstat);
        if(!channelStatus.IsGood(chan)) tjs.WireHitRange[ipl][wire] = flag;
      } // wire
    } // ipl
    
    unsigned int lastwire = 0, lastipl = 0;
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].ArtPtr->WireID().Cryostat != cstat) continue;
      if(tjs.fHits[iht].ArtPtr->WireID().TPC != tpc) continue;
      unsigned short ipl = tjs.fHits[iht].ArtPtr->WireID().Plane;
      unsigned int wire = tjs.fHits[iht].ArtPtr->WireID().Wire;
      if(wire > tjs.NumWires[ipl] - 1) {
        mf::LogWarning("TC")<<"FillWireHitRange: Invalid wire number "<<wire<<" > "<<tjs.NumWires[ipl] - 1<<" in plane "<<ipl<<" Quitting";
        return false;
      } // too large wire number
      if(ipl == lastipl && wire < lastwire) {
        mf::LogWarning("TC")<<"FillWireHitRange: Hits are not in increasing wire order. Quitting ";
        return false;
      } // hits out of order
      lastwire = wire;
      lastipl = ipl;
      if(tjs.FirstWire[ipl] == INT_MAX) tjs.FirstWire[ipl] = wire;
      if(tjs.WireHitRange[ipl][wire].first < 0) tjs.WireHitRange[ipl][wire].first = iht;
      tjs.WireHitRange[ipl][wire].second = iht + 1;
      tjs.LastWire[ipl] = wire + 1;
    } // iht
    
    if(!CheckWireHitRange(tjs)) return false;
    
    // Find the average multiplicity 1 hit RMS and calculate the expected max RMS for each range
    if(debugMode && (int)tpc == debug.TPC) {
      std::cout<<"tpc "<<tpc<<" tjs.UnitsPerTick "<<std::setprecision(3)<<tjs.UnitsPerTick<<"\n";
      std::cout<<"Fiducial volume (";
      std::cout<<std::fixed<<std::setprecision(1)<<tjs.XLo<<" < X < "<<tjs.XHi<<") (";
      std::cout<<std::fixed<<std::setprecision(1)<<tjs.YLo<<" < Y < "<<tjs.YHi<<") (";
      std::cout<<std::fixed<<std::setprecision(1)<<tjs.ZLo<<" < Z < "<<tjs.ZHi<<")\n";
    }
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      float sumRMS = 0;
      float sumAmp = 0;
      unsigned int cnt = 0;
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        // don't let noisy wires screw up the calculation
        if(lastHit - firstHit > 100) continue;
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].Multiplicity != 1) continue;
          if(tjs.fHits[iht].GoodnessOfFit < 0 || tjs.fHits[iht].GoodnessOfFit > 100) continue;
          // don't let a lot of runt hits screw up the calculation
          if(tjs.fHits[iht].PeakAmplitude < 1) continue;
          ++cnt;
          sumRMS += tjs.fHits[iht].RMS;
          sumAmp += tjs.fHits[iht].PeakAmplitude;
        } // iht
      } // wire
      if(cnt < 4) continue;
      tjs.AveHitRMS[ipl] = sumRMS/(float)cnt;
      sumAmp  /= (float)cnt;
      if(debugMode) std::cout<<"Pln "<<ipl<<" tjs.AveHitRMS "<<tjs.AveHitRMS[ipl]<<" Ave PeakAmplitude "<<sumAmp<<"\n";
    } // ipl
    return true;
    
  } // FillWireHitRange
  
  ////////////////////////////////////////////////
  bool CheckWireHitRange(const TjStuff& tjs)
  {
    // do a QC check
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        // No hits or dead wire
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        if(lastHit > tjs.fHits.size()) {
          mf::LogWarning("TC")<<"CheckWireHitRange: Invalid lastHit "<<lastHit<<" > fHits.size "<<tjs.fHits.size()<<" in plane "<<ipl;
          std::cout<<"CheckWireHitRange: Invalid lastHit "<<lastHit<<" > fHits.size "<<tjs.fHits.size()<<" in plane "<<ipl<<"\n";
          return false;
        }
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].ArtPtr->WireID().Plane != ipl) {
            mf::LogWarning("TC")<<"CheckWireHitRange: Invalid plane "<<tjs.fHits[iht].ArtPtr->WireID().Plane<<" != "<<ipl;
            std::cout<<"CheckWireHitRange: Invalid plane "<<tjs.fHits[iht].ArtPtr->WireID().Plane<<" != "<<ipl<<"\n";
            return false;
          }
          if(tjs.fHits[iht].ArtPtr->WireID().Wire != wire) {
            mf::LogWarning("TC")<<"CheckWireHitRange: Invalid wire "<<tjs.fHits[iht].ArtPtr->WireID().Wire<<" != "<<wire<<" in plane "<<ipl;
            std::cout<<"CheckWireHitRange: Invalid wire "<<tjs.fHits[iht].ArtPtr->WireID().Wire<<" != "<<wire<<" in plane "<<ipl<<"\n";
            return false;
          }
        } // iht
      } // wire
    } // ipl
    
    return true;
    
  } // CheckWireHitRange

  ////////////////////////////////////////////////
  bool WireHitRangeOK(const TjStuff& tjs, const CTP_t& inCTP)
  {
    // returns true if the passed CTP code is consistent with the CT code of the WireHitRangeVector
    geo::PlaneID planeID = DecodeCTP(inCTP);
    if(planeID.Cryostat != tjs.TPCID.Cryostat) return false;
    if(planeID.TPC != tjs.TPCID.TPC) return false;
    return true;
  }
  
  ////////////////////////////////////////////////
  bool MergeAndStore(TjStuff& tjs, unsigned int itj1, unsigned int itj2, bool doPrt)
  {
    // Merge the two trajectories in allTraj and store them. Returns true if it was successfull.
    // Merging is done between the end (end = 1) of tj1 and the beginning (end = 0) of tj2. This function preserves the
    // AlgMod state of itj1.
    // The itj1 -> itj2 merge order is reversed if end1 of itj2 is closer to end0 of itj1
    
    if(itj1 > tjs.allTraj.size() - 1) return false;
    if(itj2 > tjs.allTraj.size() - 1) return false;
    if(tjs.allTraj[itj1].AlgMod[kKilled] || tjs.allTraj[itj2].AlgMod[kKilled]) return false;
    
    // Merging shower Tjs requires merging the showers as well.
    if(tjs.allTraj[itj1].AlgMod[kShowerTj] || tjs.allTraj[itj2].AlgMod[kShowerTj]) return MergeShowerTjsAndStore(tjs, itj1, itj2, doPrt);
    
    // Ensure that the order of 3D-matched Tjs is consistent with the convention that 
    unsigned short pfp1 = GetPFPIndex(tjs, tjs.allTraj[itj1].ID);
    unsigned short pfp2 = GetPFPIndex(tjs, tjs.allTraj[itj2].ID);
    if(pfp1 == USHRT_MAX || pfp2 == USHRT_MAX) {
      if(pfp1 != USHRT_MAX && pfp2 != USHRT_MAX) {
        std::cout<<"MAS: Both tjs are used in a PFParticle. Need PFParticle merging code to do this. pfps size "<<tjs.pfps.size()<<"\n";
        return false;
      }
      // Swap so that the order of tj1 is preserved. Tj2 may be reversed to be consistent
      if(pfp1 == USHRT_MAX) std::swap(itj1, itj2);
    } // one or both used in a PFParticle
        
    // make copies so they can be trimmed as needed
    Trajectory tj1 = tjs.allTraj[itj1];
    Trajectory tj2 = tjs.allTraj[itj2];
    
    // ensure that these are in the same step order
    if(tj1.StepDir != tj2.StepDir) {
      // See if the direction has been set elsewhere
      if(tj1.AlgMod[kSetDir] || tj2.AlgMod[kSetDir]) {
        if(tj1.AlgMod[kSetDir]) {
          ReverseTraj(tjs, tj2);
        } else {
          ReverseTraj(tjs, tj1);
        }
      } else {
        ReverseTraj(tjs, tj2);
      }
    } // inconsistent step direction
    
    Point2_t tp1e0 = tj1.Pts[tj1.EndPt[0]].Pos;
    Point2_t tp1e1 = tj1.Pts[tj1.EndPt[1]].Pos;
    Point2_t tp2e0 = tj2.Pts[tj2.EndPt[0]].Pos;
    Point2_t tp2e1 = tj2.Pts[tj2.EndPt[1]].Pos;
    
    if(doPrt) {
      mf::LogVerbatim("TC")<<"MergeAndStore: tj1.ID "<<tj1.ID<<" tj2.ID "<<tj2.ID<<" at merge points "<<PrintPos(tjs, tp1e1)<<" "<<PrintPos(tjs, tp2e0);
    }
    
    // swap the order so that abs(tj1end1 - tj2end0) is less than abs(tj2end1 - tj1end0)
    if(PosSep2(tp1e1, tp2e0) > PosSep2(tp2e1, tp1e0)) {
      std::swap(tj1, tj2);
      std::swap(tp1e0, tp2e0);
      std::swap(tp1e1, tp2e1);
      if(doPrt) mf::LogVerbatim("TC")<<" swapped the order. Merge points "<<PrintPos(tjs, tp1e1)<<" "<<PrintPos(tjs, tp2e0);
    }
    
    // Here is what we are looking for, where - indicates a TP with charge.
    // Note that this graphic is in the stepping direction (+1 = +wire direction)
    // tj1:  0------------1
    // tj2:                  0-----------1
    // Another possibility with overlap
    // tj1:  0-------------1
    // tj2:               0--------------1
    
    if(tj1.StepDir > 1) {
      // Not allowed
      // tj1:  0---------------------------1
      // tj2:                  0------1
      if(tp2e0[0] > tp1e0[0] && tp2e1[0] < tp1e1[0]) return false;
      /// Not allowed
      // tj1:                  0------1
      // tj2:  0---------------------------1
      if(tp1e0[0] > tp2e0[0] && tp1e1[0] < tp2e1[0]) return false;
    } else {
      // same as above but with ends reversed
      if(tp2e1[0] > tp1e1[0] && tp2e0[0] < tp1e0[0]) return false;
      if(tp1e1[0] > tp2e1[0] && tp1e0[0] < tp2e0[0]) return false;
    }
    
    if(tj1.VtxID[1] > 0 && tj2.VtxID[0] == tj1.VtxID[1]) {
      auto& vx = tjs.vtx[tj1.VtxID[1] - 1];
      if(!MakeVertexObsolete(tjs, vx, false)) {
        if(doPrt) mf::LogVerbatim("TC")<<"MergeAndStore: Found a good vertex between Tjs "<<tj1.VtxID[1]<<" No merging";
        return false;
      }
    }
    
    if(tj1.StopFlag[1][kBragg]) {
      if(doPrt) mf::LogVerbatim("TC")<<"MergeAndStore: You are merging the end of a trajectory "<<tj1.ID<<" with a Bragg peak. Not merging\n";
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
    TrajPointTrajDOCA(tjs, endtj1TP, tj2, tj2ClosePt, minSep);
    if(doPrt) mf::LogVerbatim("TC")<<" Merge point tj1 "<<PrintPos(tjs, endtj1TP)<<" tj2ClosePt "<<tj2ClosePt<<" Pos "<<PrintPos(tjs, tj2.Pts[tj2ClosePt]);
    // check for full overlap
    if(tj2ClosePt > tj2.EndPt[1]) return false;
    
    // The approach is to append tj2 to tj1, store tj1 as a new trajectory,
    // and re-assign all hits to the new trajectory
    
    // First ensure that any hit will appear only once in the merged trajectory in the overlap region
    // whether it is used or unused. The point on tj2 where the merge will begin, tj2ClosePt, will be
    // increased until this condition is met.
    // Make a temporary vector of tj1 hits in the end points for simpler searching
    std::vector<unsigned int> tj1Hits;
    for(unsigned short ii = 0; ii < tj1.Pts.size(); ++ii) {
      // only go back a few points in tj1
      if(ii > 10) break;
      unsigned short ipt = tj1.Pts.size() - 1 - ii;
      tj1Hits.insert(tj1Hits.end(), tj1.Pts[ipt].Hits.begin(), tj1.Pts[ipt].Hits.end());
      if(ipt == 0) break;
    } // ii
    
    bool bumpedPt = true;
    while(bumpedPt) {
      bumpedPt = false;
      for(unsigned short ii = 0; ii < tj2.Pts[tj2ClosePt].Hits.size(); ++ii) {
        unsigned int iht = tj2.Pts[tj2ClosePt].Hits[ii];
        if(std::find(tj1Hits.begin(), tj1Hits.end(), iht) != tj1Hits.end()) bumpedPt = true;
      } // ii
      if(bumpedPt && tj2ClosePt < tj2.EndPt[1]) {
        ++tj2ClosePt;
      } else {
        break;
      }
    } // bumpedPt
    if(doPrt) mf::LogVerbatim("TC")<<" revised tj2ClosePt "<<tj2ClosePt;
    // append tj2 hits to tj1
    
    tj1.Pts.insert(tj1.Pts.end(), tj2.Pts.begin() + tj2ClosePt, tj2.Pts.end());
    // re-define the end points
    SetEndPoints(tjs, tj1);
    tj1.StopFlag[1] = tj2.StopFlag[1];
    
    // A more exhaustive check that hits only appear once
    if(HasDuplicateHits(tjs, tj1, doPrt)) {
      if(doPrt) {
        mf::LogVerbatim("TC")<<"MergeAndStore found duplicate hits. Coding error";
        PrintTrajectory("MAS", tjs, tj1, USHRT_MAX);
        PrintTrajectory("tj1", tjs, tjs.allTraj[itj1], USHRT_MAX);
        PrintTrajectory("tj2", tjs, tjs.allTraj[itj2], USHRT_MAX);
      }
      return false;
    }
    if(tj2.VtxID[1] > 0) {
//      std::cout<<"MAS: Preserve vertex "<<tj2.VtxID[1]<<" merging "<<tj1.ID<<" and "<<tj2.ID<<" \n";
      // move the end vertex of tj2 to the end of tj1
      tj1.VtxID[1] = tj2.VtxID[1];
    }   
    // Transfer some of the AlgMod bits
    if(tj2.AlgMod[kMichel]) tj1.AlgMod[kMichel] = true;
    if(tj2.AlgMod[kDeltaRay]) {
      tj1.AlgMod[kDeltaRay] = true;
      tj1.ParentID = tj2.ParentID;
    }
    // keep track of the IDs before they are clobbered
    int tj1ID = tj1.ID;
    int tj2ID = tj2.ID;
    // kill the original trajectories
    MakeTrajectoryObsolete(tjs, itj1);
    MakeTrajectoryObsolete(tjs, itj2);
    // Do this so that StoreTraj keeps the correct WorkID (of itj1)
    tj1.ID = tj1.WorkID;
    SetPDGCode(tjs, tj1);
    if(!StoreTraj(tjs, tj1)) return false;
    int newTjID = tjs.allTraj.size();
    // Use the ParentID to trace which new Tj is superseding the merged ones
    tj1.ParentID = newTjID;
    tj2.ParentID = newTjID;
    // update match structs if they exist
    std::vector<int> oldTjs(2);
    oldTjs[0] = tj1.ID;
    oldTjs[1] = tj2.ID;
    UpdateMatchStructs(tjs, oldTjs, newTjID);
    if(doPrt) mf::LogVerbatim("TC")<<" MAS success. New TjID "<<newTjID;
    // Transfer the ParentIDs of any other Tjs that refer to Tj1 and Tj2 to the new Tj
    for(auto& tj : tjs.allTraj) if(tj.ParentID == tj1ID || tj.ParentID == tj2ID) tj.ParentID = newTjID;

    return true;
  } // MergeAndStore
  
  // ****************************** Printing  ******************************
  
  ////////////////////////////////////////////////
  void PrintAllTraj(std::string someText, const TjStuff& tjs, const DebugStuff& debug, unsigned short itj, unsigned short ipt, bool prtVtx)
  {
    
    mf::LogVerbatim myprt("TC");
    
    if(prtVtx) {
      if(!tjs.vtx3.empty()) {
        // print out 3D vertices
        myprt<<someText<<"****** 3D vertices ******************************************__2DVtx_ID__*******\n";
        myprt<<someText<<"Vtx  Cstat  TPC     X       Y       Z    XEr  YEr  ZEr pln0 pln1 pln2 Wire score Prim? Nu? nTru";
        myprt<<" ___________2D_Pos____________ _____Tjs________\n";
        for(unsigned short iv = 0; iv < tjs.vtx3.size(); ++iv) {
          if(tjs.vtx3[iv].ID == 0) continue;
          const Vtx3Store& vx3 = tjs.vtx3[iv];
          myprt<<someText;
          myprt<<std::right<<std::setw(3)<<std::fixed<<vx3.ID<<std::setprecision(1);
          myprt<<std::right<<std::setw(7)<<vx3.TPCID.Cryostat;
          myprt<<std::right<<std::setw(5)<<vx3.TPCID.TPC;
          myprt<<std::right<<std::setw(8)<<vx3.X;
          myprt<<std::right<<std::setw(8)<<vx3.Y;
          myprt<<std::right<<std::setw(8)<<vx3.Z;
          myprt<<std::right<<std::setw(5)<<vx3.XErr;
          myprt<<std::right<<std::setw(5)<<vx3.YErr;
          myprt<<std::right<<std::setw(5)<<vx3.ZErr;
          myprt<<std::right<<std::setw(5)<<vx3.Vx2ID[0];
          myprt<<std::right<<std::setw(5)<<vx3.Vx2ID[1];
          myprt<<std::right<<std::setw(5)<<vx3.Vx2ID[2];
          myprt<<std::right<<std::setw(5)<<vx3.Wire;
          unsigned short nTruMatch = 0;
          for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
            if(vx3.Vx2ID[ipl] == 0) continue;
            unsigned short iv2 = vx3.Vx2ID[ipl] - 1;
            if(tjs.vtx[iv2].Stat[kVtxTruMatch]) ++nTruMatch;
          } // ipl
          myprt<<std::right<<std::setw(6)<<std::setprecision(1)<<vx3.Score;
          myprt<<std::setw(6)<<vx3.Primary;
          myprt<<std::setw(4)<<vx3.Neutrino;
          myprt<<std::right<<std::setw(5)<<nTruMatch;
          Point2_t pos;
          for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
            PosInPlane(tjs, vx3, plane, pos);
            myprt<<" "<<PrintPos(tjs, pos);
          } // plane
          if(vx3.Wire == -2) {
            // find the Tjs that are attached to it
            myprt<<" Tjs";
            for(auto& pfp : tjs.pfps) {
              if(pfp.Vx3ID[0] == tjs.vtx3[iv].ID) {
                for(auto& tjID : pfp.TjIDs) myprt<<" "<<tjID;
              }
              if(pfp.Vx3ID[1] == tjs.vtx3[iv].ID) {
                for(auto& tjID : pfp.TjIDs) myprt<<" "<<tjID;
              }
            } // ipfp
          } else {
            float score;
            auto vxtjs = GetVtxTjIDs(tjs, vx3, score);
            myprt<<" Tjs";
            for(auto tjid : vxtjs) myprt<<" "<<tjid;
          }
          myprt<<"\n";
        }
      } // tjs.vtx3.size
      if(!tjs.vtx.empty()) {
        bool foundOne = false;
        for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
          auto& vx2 = tjs.vtx[iv];
          if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(vx2.CTP).Plane) continue;
          if(vx2.NTraj == 0) continue;
          foundOne = true;
        } // iv
        if(foundOne) {
          // print out 2D vertices
          myprt<<someText<<"************ 2D vertices ************\n";
          myprt<<someText<<"VtxID  CTP   wire  err   tick   err  ChiDOF  NTj Pass  Topo ChgFrac Score  v3D TjIDs\n";
          for(auto& vx2 : tjs.vtx) {
            if(vx2.ID == 0) continue;
            if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(vx2.CTP).Plane) continue;
            myprt<<someText;
            myprt<<std::right<<std::setw(3)<<std::fixed<<vx2.ID;
            myprt<<std::right<<std::setw(6)<<vx2.CTP;
            myprt<<std::right<<std::setw(8)<<std::setprecision(0)<<std::nearbyint(vx2.Pos[0]);
            myprt<<std::right<<std::setw(5)<<std::setprecision(1)<<vx2.PosErr[0];
            myprt<<std::right<<std::setw(8)<<std::setprecision(0)<<std::nearbyint(vx2.Pos[1]/tjs.UnitsPerTick);
            myprt<<std::right<<std::setw(5)<<std::setprecision(1)<<vx2.PosErr[1]/tjs.UnitsPerTick;
            myprt<<std::right<<std::setw(7)<<vx2.ChiDOF;
            myprt<<std::right<<std::setw(5)<<vx2.NTraj;
            myprt<<std::right<<std::setw(5)<<vx2.Pass;
            myprt<<std::right<<std::setw(6)<<vx2.Topo;
            myprt<<std::right<<std::setw(9)<<std::setprecision(2)<<vx2.TjChgFrac;
            myprt<<std::right<<std::setw(6)<<std::setprecision(1)<<vx2.Score;
            myprt<<std::right<<std::setw(5)<<vx2.Vx3ID;
            myprt<<"    ";
            // display the traj IDs
            for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
              auto const& aTj = tjs.allTraj[ii];
              if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
              if(aTj.AlgMod[kKilled]) continue;
              for(unsigned short end = 0; end < 2; ++end)
                if(aTj.VtxID[end] == (short)vx2.ID) myprt<<std::right<<std::setw(4)<<aTj.ID<<"_"<<end;
            }
            // Special flags. Ignore the first flag bit (0 = kVtxTrjTried) which is done for every vertex
            for(unsigned short ib = 1; ib < VtxBitNames.size(); ++ib) if(vx2.Stat[ib]) myprt<<" "<<VtxBitNames[ib];
            myprt<<"\n";
          } // iv
        }
      } // tjs.vtx.size
    }
     
    if(tjs.allTraj.empty()) {
      mf::LogVerbatim("TC")<<someText<<" No allTraj trajectories to print";
      return;
    }
    
    // Print all trajectories in tjs.allTraj if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      std::vector<unsigned int> tmp;
      myprt<<someText<<" TRJ  ID   CTP Pass  Pts     W:T      Ang CS AveQ dEdx     W:T      Ang CS AveQ dEdx chgRMS Mom SDr NN __Vtx__  PDG  Par Pri NuPar TRuPDG  E*P TruKE  WorkID \n";
      for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
        auto& aTj = tjs.allTraj[ii];
        if(debug.Plane >=0 && debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
        myprt<<someText<<" ";
        if(aTj.AlgMod[kKilled]) { myprt<<"xxx"; } else { myprt<<"TRJ"; }
        myprt<<std::fixed<<std::setw(4)<<aTj.ID;
        myprt<<std::setw(6)<<aTj.CTP;
        myprt<<std::setw(5)<<aTj.Pass;
//        myprt<<std::setw(5)<<aTj.Pts.size();
        myprt<<std::setw(5)<<aTj.EndPt[1] - aTj.EndPt[0] + 1;
        unsigned short endPt0 = aTj.EndPt[0];
        auto& tp0 = aTj.Pts[endPt0];
        int itick = tp0.Pos[1]/tjs.UnitsPerTick;
        if(itick < 0) itick = 0;
        myprt<<std::setw(6)<<(int)(tp0.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) { myprt<<" "; }
        if(itick < 100) { myprt<<" "; }
        if(itick < 1000) { myprt<<" "; }
        myprt<<std::setw(6)<<std::setprecision(2)<<tp0.Ang;
        myprt<<std::setw(2)<<tp0.AngleCode;
        if(aTj.StopFlag[0][kBragg]) {
          myprt<<"B";
        } else if(aTj.StopFlag[0][kAtVtx]) {
          myprt<<"V";
        } else {
          myprt<<" ";
        }
        myprt<<std::setw(5)<<(int)tp0.AveChg;
        unsigned short prec = 1;
        if(aTj.dEdx[0] > 99) prec = 0;
        myprt<<std::setw(5)<<std::setprecision(prec)<<aTj.dEdx[0];
        unsigned short endPt1 = aTj.EndPt[1];
        auto& tp1 = aTj.Pts[endPt1];
        itick = tp1.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp1.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) { myprt<<" "; }
        if(itick < 100) { myprt<<" "; }
        if(itick < 1000) { myprt<<" "; }
        myprt<<std::setw(6)<<std::setprecision(2)<<tp1.Ang;
        myprt<<std::setw(2)<<tp1.AngleCode;
        if(aTj.StopFlag[1][kBragg]) {
          myprt<<"B";
        } else if(aTj.StopFlag[1][kAtVtx]) {
          myprt<<"V";
        } else {
          myprt<<" ";
        }
        myprt<<std::setw(5)<<(int)tp1.AveChg;
        prec = 1;
        if(aTj.dEdx[1] > 99) prec = 0;
        myprt<<std::setw(5)<<std::setprecision(prec)<<aTj.dEdx[1];
        myprt<<std::setw(7)<<std::setprecision(2)<<aTj.ChgRMS;
        myprt<<std::setw(5)<<aTj.MCSMom;
        myprt<<std::setw(4)<<aTj.StepDir;
        myprt<<std::setw(3)<<aTj.NNeighbors;
        myprt<<std::setw(4)<<aTj.VtxID[0];
        myprt<<std::setw(4)<<aTj.VtxID[1];
        myprt<<std::setw(5)<<aTj.PDGCode;
        myprt<<std::setw(5)<<aTj.ParentID;
        myprt<<std::setw(5)<<PrimaryID(tjs, aTj);
        myprt<<std::setw(6)<<NeutrinoPrimaryTjID(tjs, aTj);
        int truKE = 0;
        int pdg = 0;
        if(aTj.MCPartListIndex < tjs.MCPartList.size()) {
          auto& mcp = tjs.MCPartList[aTj.MCPartListIndex];
          truKE = 1000 * (mcp->E() - mcp->Mass());
          pdg = mcp->PdgCode();
        }
        myprt<<std::setw(6)<<pdg;
        myprt<<std::setw(6)<<std::setprecision(2)<<aTj.EffPur;
        myprt<<std::setw(5)<<truKE;
        myprt<<std::setw(7)<<aTj.WorkID;
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(aTj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
      } // ii
      return;
    } // itj > tjs.allTraj.size()-1
    
    if(itj > tjs.allTraj.size()-1) return;
    
    auto const& aTj = tjs.allTraj[itj];
    
    mf::LogVerbatim("TC")<<"Print tjs.allTraj["<<itj<<"]: ClusterIndex "<<aTj.ClusterIndex<<" Vtx[0] "<<aTj.VtxID[0]<<" Vtx[1] "<<aTj.VtxID[1];
    myprt<<"AlgBits";
    for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(aTj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
    myprt<<"\n";
    
    PrintHeader(someText);
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < aTj.Pts.size(); ++ii) PrintTrajPoint(someText, tjs, ii, aTj.StepDir, aTj.Pass, aTj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(someText, tjs, ipt, aTj.StepDir, aTj.Pass, aTj.Pts[ipt]);
    }
  } // PrintAllTraj
  
  
  //////////////////////////////////////////
  void PrintTrajectory(std::string someText, const TjStuff& tjs, const Trajectory& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    int trupdg = 0;
    if(tj.MCPartListIndex < tjs.MCPartList.size()) trupdg = tjs.MCPartList[tj.MCPartListIndex]->PdgCode();
    
    if(tPoint == USHRT_MAX) {
      if(tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ";
        myprt<<"Work:    ID "<<tj.ID<<"    CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<trupdg<<" tjs.vtx "<<tj.VtxID[0]<<" "<<tj.VtxID[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        myprt<<" MCSMom "<<tj.MCSMom;
        myprt<<" StopFlags "<<PrintStopFlag(tj, 0)<<" "<<PrintStopFlag(tj, 1);
        myprt<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      } else {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ";
        myprt<<"tjs.allTraj: ID "<<tj.ID<<" WorkID "<<tj.WorkID<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<trupdg<<" tjs.vtx "<<tj.VtxID[0]<<" "<<tj.VtxID[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        myprt<<" MCSMom "<<tj.MCSMom;
        myprt<<" StopFlags "<<PrintStopFlag(tj, 0)<<" "<<PrintStopFlag(tj, 1);
        myprt<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      }
      PrintHeader(someText);
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) PrintTrajPoint(someText, tjs, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
      // See if this trajectory is a shower Tj
      if(tj.AlgMod[kShowerTj]) {
        for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
          if(tjs.cots[ic].TjIDs.empty()) continue;
          // only print out the info for the correct Tj
          if(tjs.cots[ic].ShowerTjID != tj.ID) continue;
          const ShowerStruct& ss = tjs.cots[ic];
          mf::LogVerbatim myprt("TC");
          myprt<<"cots index "<<ic<<" ";
          myprt<<someText<<" Envelope";
          if(ss.Envelope.empty()) {
            myprt<<" NA";
          } else {
            for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
          }
          myprt<<" Energy "<<(int)ss.Energy;
          myprt<<" Area "<<std::fixed<<std::setprecision(1)<<(int)ss.EnvelopeArea<<" ChgDensity "<<ss.ChgDensity;
          myprt<<"\nInShower TjIDs";
          for(auto& tjID : ss.TjIDs) {
            myprt<<" "<<tjID;
          } // tjIDA_Klystron_4U
          
          myprt<<"\n";
          myprt<<"NearTjIDs";
          for(auto& tjID : ss.NearTjIDs) {
            myprt<<" "<<tjID;
          } // tjID
          myprt<<"\n";
          myprt<<"\n";
          myprt<<"Angle "<<std::fixed<<std::setprecision(2)<<ss.Angle<<" +/- "<<ss.AngleErr;
          myprt<<" AspectRatio "<<std::fixed<<std::setprecision(2)<<ss.AspectRatio;
          myprt<<" DirectionFOM "<<std::fixed<<std::setprecision(2)<<ss.DirectionFOM;
          if(ss.ParentID > 0) {
            myprt<<" Parent Tj "<<ss.ParentID<<" FOM "<<ss.ParentFOM;
          } else {
            myprt<<" No parent";
          }
          myprt<<" TruParentID "<<ss.TruParentID<<" SS3ID "<<ss.SS3ID<<"\n";
          if(ss.NeedsUpdate) myprt<<"*********** This shower needs to be updated ***********";
          myprt<<"................................................";
        } // ic
      } // Shower Tj
    } else {
      // just print one traj point
      if(tPoint > tj.Pts.size() -1) {
        mf::LogVerbatim("TC")<<"Can't print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(someText, tjs, tPoint, tj.StepDir, tj.Pass, tj.Pts[tPoint]);
    }
  } // PrintTrajectory
  
  //////////////////////////////////////////
  void PrintHeader(std::string someText)
  {
    mf::LogVerbatim("TC")<<someText<<" TRP     CTP  Ind  Stp      W:Tick    Delta  RMS    Ang C   Err  Dir0  Dir1      Q    AveQ  Pull FitChi  NTPF  Hits ";
  } // PrintHeader
  
  ////////////////////////////////////////////////
  void PrintTrajPoint(std::string someText, const TjStuff& tjs, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" TRP"<<std::fixed;
    myprt<<pass;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(6)<<tp.CTP;
    myprt<<std::setw(5)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(7)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]/tjs.UnitsPerTick; // W:T
    if(tp.Pos[1] < 10) { myprt<<"  "; }
    if(tp.Pos[1] < 100) { myprt<<" "; }
    if(tp.Pos[1] < 1000) { myprt<<" "; }
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.DeltaRMS;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(2)<<tp.AngleCode;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[0];
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[1];
    myprt<<std::setw(7)<<(int)tp.Chg;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgPull;
    myprt<<std::setw(7)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NTPsFit;
    // print the hits associated with this traj point
    if(tp.Hits.size() > 16) {
      // don't print too many hits (e.g. from a shower Tj)
      myprt<<" "<<tp.Hits.size()<<" shower hits";
    } else {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        if(iht > tjs.fHits.size() - 1) {
          std::cout<<"crazy hit "<<iht<<" CTP "<<tp.CTP<<"\n";
          continue;
        }
        myprt<<" "<<tjs.fHits[iht].ArtPtr->WireID().Wire<<":"<<(int)tjs.fHits[iht].PeakTime;
        if(tp.UseHit[ii]) {
          // Distinguish used hits from nearby hits
          myprt<<"_";
        } else {
          myprt<<"x";
        }
        myprt<<tjs.fHits[iht].InTraj;
      } // iht
    }
  } // PrintTrajPoint
  
  /////////////////////////////////////////
  void PrintPFParticles(std::string someText, const TjStuff& tjs)
  {
    if(tjs.pfps.empty()) return;
    
    mf::LogVerbatim myprt("TC");
    myprt<<someText;
    myprt<<"  PFP sVx  ________sPos_______  ______sDir______  ______sdEdx_____ eVx  ________ePos_______  ______eDir______  ______edEdx_____ BstPln PDG TruPDG Par Prim E*P\n";
    unsigned short indx = 0;
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      myprt<<someText;
      myprt<<std::setw(5)<<pfp.ID;
      // start and end stuff
      for(unsigned short startend = 0; startend < 2; ++startend) {
        myprt<<std::setw(4)<<pfp.Vx3ID[startend];
        myprt<<std::fixed<<std::right<<std::setprecision(1);
        myprt<<std::setw(7)<<pfp.XYZ[startend][0];
        myprt<<std::setw(7)<<pfp.XYZ[startend][1];
        myprt<<std::setw(7)<<pfp.XYZ[startend][2];
        myprt<<std::fixed<<std::right<<std::setprecision(2);
        myprt<<std::setw(6)<<pfp.Dir[startend][0];
        myprt<<std::setw(6)<<pfp.Dir[startend][1];
        myprt<<std::setw(6)<<pfp.Dir[startend][2];
        for(auto& dedx : pfp.dEdx[startend]) {
          if(dedx < 50) {
            myprt<<std::setw(6)<<std::setprecision(1)<<dedx;
          } else {
            myprt<<std::setw(6)<<std::setprecision(0)<<dedx;
          }
        } // dedx
        if (pfp.dEdx[startend].size()<3){
          for(size_t i = 0; i<3-pfp.dEdx[startend].size(); ++i){
            myprt<<std::setw(6)<<' ';
          }
        }
      }
      // global stuff
      myprt<<std::setw(5)<<pfp.BestPlane;
      myprt<<std::setw(6)<<pfp.PDGCode;
      if(pfp.MCPartListIndex < tjs.MCPartList.size()) {
        myprt<<std::setw(6)<<tjs.MCPartList[pfp.MCPartListIndex]->PdgCode();
      } else {
        myprt<<"    NA";
      }
      myprt<<std::setw(4)<<pfp.ParentID;
      myprt<<std::setw(5)<<PrimaryID(tjs, pfp);
      myprt<<std::setw(5)<<std::setprecision(2)<<pfp.EffPur;
      if(!pfp.TjIDs.empty()) {
        myprt<<" tjs";
        for(auto& tjID : pfp.TjIDs) myprt<<" "<<tjID;
      }
      if(!pfp.DtrIDs.empty()) {
        myprt<<" dtrs";
        for(auto& dtrID : pfp.DtrIDs) myprt<<" "<<dtrID;
      }
      myprt<<"\n";
      ++indx;
    } // im
    
  } // PrintPFParticles
  
  /////////////////////////////////////////
  std::string PrintStopFlag(const Trajectory& tj, unsigned short end)
  {
    if(end > 1) return "Invalid end";
    std::string tmp;
    bool first = true;
    for(unsigned short ib = 0; ib < StopFlagNames.size(); ++ib) {
      if(tj.StopFlag[end][ib]) {
        if(first) {
          tmp = std::to_string(end) + ":" + StopFlagNames[ib];
          first = false;
        } else {
          tmp += "," + StopFlagNames[ib];
        }
      }
    } // ib
    return tmp;
  } // PrintStopFlag
  
  /////////////////////////////////////////
  std::string PrintHitShort(const TCHit& hit)
  {
    return std::to_string(hit.ArtPtr->WireID().Plane) + ":" + std::to_string(hit.ArtPtr->WireID().Wire) + ":" + std::to_string((int)hit.PeakTime);
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintHit(const TCHit& hit)
  {
    return std::to_string(hit.ArtPtr->WireID().Plane) + ":" + std::to_string(hit.ArtPtr->WireID().Wire) + ":" + std::to_string((int)hit.PeakTime) + "_" + std::to_string(hit.InTraj);
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintPos(const TjStuff& tjs, const TrajPoint& tp)
  {
    return std::to_string(tp.CTP) + ":" + PrintPos(tjs, tp.Pos);
  } // PrintPos
  
  /////////////////////////////////////////
  std::string PrintPos(const TjStuff& tjs, const Point2_t& pos)
  {
    unsigned int wire = std::nearbyint(pos[0]);
    int time = std::nearbyint(pos[1]/tjs.UnitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

  
} // namespace tca

