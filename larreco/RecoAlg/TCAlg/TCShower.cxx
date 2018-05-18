#include "larreco/RecoAlg/TCAlg/TCShower.h"

struct SortEntry{
  unsigned int index;
  float length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}

namespace tca {

  ////////////////////////////////////////////////
  bool FindShowerStart(TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // The shower ChgPos and Dir were found by the calling function but Dir 
    // may be inconsistent with the 2D shower directions
    if(ss3.ID == 0) return false;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Inside FSS: 3S"<<ss3.ID<<" ->";
      for(auto cid : ss3.CotIDs) myprt<<" 2S"<<cid;
      myprt<<" Vx 3V"<<ss3.Vx3ID;
    } // prt
    
    // find a parent Tj in the list of 2D showers that is the farthest away from the
    // shower center
    unsigned short useParentCID = 0;
    float maxParentSep = 0;
    unsigned short usePtSepCID = 0;
    float maxPtSep = 0;
    // assume that the 2D shower direction is consistent with the 3D shower direction. This
    // variable is only used when no parent exists
    bool dirOK = true;
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      // Find the position, direction and projection in this plane
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      auto chgCtrTP = MakeBareTP(tjs, ss3.ChgPos, ss3.Dir, stj.CTP);
      // projection too small in this view?
      if(chgCtrTP.Delta < 0.5) continue;
      auto& startTP = stj.Pts[0];
      float sep = PosSep(startTP.Pos, chgCtrTP.Pos);
      if(ss.ParentID > 0) {
        if(sep > maxParentSep) {
          maxParentSep = sep;
          useParentCID = cid;
        }
      } else if(sep > maxPtSep) {
        // no parent exists.
        maxPtSep = sep;
        usePtSepCID = cid;
        float costh = DotProd(chgCtrTP.Dir, startTP.Dir);
        if(costh < 0) dirOK = false;
      }
    } // ci
    if(useParentCID == 0 && usePtSepCID == 0) return false;
    
    unsigned short useCID = useParentCID;
    if(useCID == 0) {
      useCID = usePtSepCID;
      if(!dirOK) ReverseShower("FSS", tjs, useCID, prt);
    }

    // now define the start and length
    auto& ss = tjs.cots[useCID - 1];
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];

    auto chgCtrTP = MakeBareTP(tjs, ss3.ChgPos, ss3.Dir, stj.CTP);
    if(ss3.Vx3ID > 0) {
      auto& vx3 = tjs.vtx3[ss3.Vx3ID - 1];
      ss3.Start[0] = vx3.X;
      ss3.Start[1] = vx3.Y;
      ss3.Start[2] = vx3.Z;
    } else {
      // no start vertex
      auto& startTP = stj.Pts[0];
      // The 2D separation btw the shower start and the shower center, converted
      // to the 3D separation
//      std::cout<<"useCI "<<useCID<<" sep "<<PosSep(startTP.Pos, stj.Pts[1].Pos)<<" projInPln "<<chgCtrTP.Delta<<"\n";
      float sep = tjs.WirePitch * PosSep(startTP.Pos, stj.Pts[1].Pos) / chgCtrTP.Delta;
      // set the start position
      for(unsigned short xyz = 0; xyz < 3; ++xyz) ss3.Start[xyz] = ss3.ChgPos[xyz] - sep * ss3.Dir[xyz];
    }
    // now do the end position
    auto& endTP = stj.Pts[2];
    float sep = tjs.WirePitch * PosSep(endTP.Pos, chgCtrTP.Pos) / chgCtrTP.Delta;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) ss3.End[xyz] = ss3.ChgPos[xyz] + sep * ss3.Dir[xyz];
    ss3.Len = PosSep(ss3.Start, ss3.End);
    auto& startTP = stj.Pts[0];
    sep = PosSep(startTP.Pos, endTP.Pos);
    ss3.OpenAngle = (endTP.DeltaRMS - startTP.DeltaRMS) / sep;
    ss3.OpenAngle /= chgCtrTP.Delta;
    return true;
    
  } // FindShowerStart

  ////////////////////////////////////////////////
  void Finish3DShowers(TjStuff& tjs)
  {
    // Finish defining the showers, create a companion PFParticle for each one and define the mother
    // daughter relationships
    // Note to the reader: This code doesn't use MakeVertexObsolete to kill vertices using the assumption
    // that Finish3DShowers is being called after reconstruction is complete, in which case there is no
    // need to re-calculate the 2D and 3D vertex score which could potentially screw up the decisions that have
    // already been made.
    
    // See if any need to be finished
    bool noShowers = true;
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      noShowers = false;
    }
    if(noShowers) return;
    
    // Use special care if a neutrino PFParticle exists so that we don't clobber the neutrino vertex
    bool foundNeutrino = !tjs.pfps.empty() && (tjs.pfps[0].PDGCode == 12 || tjs.pfps[0].PDGCode == 14);
    
    // create a pfp and define the mother-daughter relationship. At this point, the shower parent PFP (if
    // one exists) is a track-like pfp that might be the daughter of another pfp, e.g. the neutrino. This
    // association is changed from shower ParentID -> parent pfp, to shower PFP -> parent pfp
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      if(ss3.PFPIndex != USHRT_MAX) {
        std::cout<<"Finish3DShowers 3S"<<ss3.ID<<" already has a pfp associated with it...\n";
        continue;
      }
      auto showerPFP = CreatePFP(tjs, ss3.TPCID);
      showerPFP.TjIDs.resize(ss3.CotIDs.size());
      for(unsigned short ii = 0; ii < ss3.CotIDs.size(); ++ii) {
        unsigned short cid = ss3.CotIDs[ii];
        if(cid == 0 || cid > tjs.cots.size()) {
          std::cout<<"Finish3DShowers 3S"<<ss3.ID<<" has an invalid cots ID"<<cid<<"\n";
          return;
        }
        auto& ss = tjs.cots[cid - 1];
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        showerPFP.TjIDs[ii] = stj.ID;
      } // ci
      showerPFP.PDGCode = 1111;
      showerPFP.XYZ[0] = ss3.Start;
      showerPFP.Dir[0] = ss3.Dir;
      showerPFP.DirErr[0] = ss3.DirErr;
      showerPFP.Vx3ID[0] = ss3.Vx3ID;
      showerPFP.XYZ[1] = ss3.End;
      showerPFP.Dir[1] = ss3.Dir;
      // dEdx is indexed by plane for pfps and by 2D shower index for 3D showers
      for(auto cid : ss3.CotIDs) {
        auto& ss = tjs.cots[cid - 1];
        unsigned short plane = DecodeCTP(ss.CTP).Plane;
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        showerPFP.dEdx[0][plane] = stj.dEdx[0];
        showerPFP.dEdxErr[0][plane] = 0.3 * stj.dEdx[0];
      } // ci
      ss3.PFPIndex = tjs.pfps.size();
      if(ss3.ParentID > 0) {
        // see if this is a daughter
        auto& dtrPFP = tjs.pfps[ss3.ParentID - 1];
        if(dtrPFP.ParentID > 0) {
          // Transfer the daughter <-> parent assn
          auto& parPFP = tjs.pfps[dtrPFP.ParentID - 1];
          showerPFP.ParentID = parPFP.ID;
          std::replace(parPFP.DtrIDs.begin(), parPFP.DtrIDs.end(), dtrPFP.ID, showerPFP.ID);
          dtrPFP.ParentID = 0;
        } // dtrPFP.ParentID > 0
      } // ss3.ParentID > 0
      tjs.pfps.push_back(showerPFP);
    } // ss3

    // Transfer Tj hits from InShower Tjs to the shower Tj. This kills the InShower Tjs but doesn't consider how
    // this action affects vertices
    if(!TransferTjHits(tjs, false)) return;
    
    // Associate shower Tj hits with 3D showers
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      for(unsigned short ii = 0; ii < ss3.CotIDs.size(); ++ii) {
        unsigned short cid = ss3.CotIDs[ii];
        auto& ss = tjs.cots[cid - 1];
        for(auto tjid : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjid - 1];
          auto tjHits = PutTrajHitsInVector(tj, kUsedHits);
          ss3.Hits.insert(ss3.Hits.end(), tjHits.begin(), tjHits.end());
          // kill vertices associated with the Tj unless it is the neutrino primary vtx
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] == 0) continue;
            auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
            if(vx2.Vx3ID <= 0) {
              // This is a 2D vertex that is not matched in 3D. Kill it. Shouldn't need to
              // use MakeVertexObsolete here...
              vx2.ID = 0;
              continue;
            }
            // vx2 is matched in 3D. Kill it if it is NOT the neutrino primary
            auto& vx3 = tjs.vtx3[vx2.Vx3ID - 1];
            if(vx3.Neutrino) continue;
            vx3.ID = 0;
          } // end
        } // tjid
      } // ii
    } // ss3
    
    // kill PFParticles
    if(!tjs.pfps.empty()) {
      for(auto& pfp : tjs.pfps) {
        if(pfp.ID == 0) continue;
        unsigned short ndead = 0;
        for(auto tjid : pfp.TjIDs) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kKilled]) ++ndead;
        } // tjid
        if(ndead == 0) continue;
        if(foundNeutrino && pfp.ParentID == 1) {
          std::cout<<"Finish3DShowers wants to kill neutrino primary P"<<pfp.ID<<". Not doing it.\n";
          continue;
        }
        pfp.ID = 0;
      } // pfp
    } // pfps not empty
    
    // kill orphan 2D vertices
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      auto vxtjs = GetVtxTjIDs(tjs, vx2);
      if(vxtjs.empty()) vx2.ID = 0;
    } // vx2
    
    // kill orphan vertices
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      float score;
      auto vxtjs = GetVtxTjIDs(tjs, vx3, score);
      if(vxtjs.empty()) {
        vx3.ID = 0;
        continue;
      }
    } // vx3
    
    // check
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.PDGCode != 1111) continue;
      if(pfp.Vx3ID[0] > 0 && tjs.vtx3[pfp.Vx3ID[0] - 1].ID == 0) {
        std::cout<<"Finish3DShowers shower P"<<pfp.ID<<" has Vx3ID[0] = "<<pfp.Vx3ID[0]<<" but the vertex is killed\n";
      } // Vx3 check
    } // pfp

  } // Finish3DShowers
  
  ////////////////////////////////////////////////
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // Find 2D showers using 3D-matched trajectories. This returns true if showers were found
    // which requires re-doing the 3D trajectory match
    
    bool reconstruct = (tjs.ShowerTag[0] == 2) || (tjs.ShowerTag[0] == 4);
    if(!reconstruct) return false;
    
    bool prt = false;
    // Add 10 for more detailed printing 
    short dbgPlane = ((int)tjs.ShowerTag[12] % 10);
    CTP_t dbgCTP = UINT_MAX;
    if(dbgPlane >= 0 && dbgPlane <= tjs.NumPlanes) dbgCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, dbgPlane);
    
    std::string fcnLabel = "FS";

    geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
    // check for already-existing showers
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      for(auto& ss : tjs.cots) if(ss.CTP == inCTP) return false;
    }
    
    // rebuild the hit range references if necessary
    if(tpcid != tjs.TPCID && !FillWireHitRange(tjs, tpcid, false)) return false;

    if(prt) {
      PrintPFPs("FSi", tjs);
      PrintAllTraj("FSi", tjs, debug, USHRT_MAX, 0);
    }

    // lists of Tj IDs in plane, (list1, list2, list3, ...)
    std::vector<std::vector<std::vector<int>>> bigList(tjs.NumPlanes);
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      std::vector<std::vector<int>> tjList;
      TagShowerLike(fcnLabel, tjs, inCTP, tjList, false);
      SaveTjInfo(tjs, inCTP, tjList, "TJL");
      if(tjList.empty()) continue;
      bigList[plane] = tjList;
    } // plane
    unsigned short nPlanesWithShowers = 0;
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) if(!bigList.empty()) ++nPlanesWithShowers;
    if(nPlanesWithShowers < 2) return false;
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      // print detailed debug info for one plane
      prt = (inCTP == dbgCTP);
      // Create a shower for each one
      for(auto& tjl : bigList[plane]) {
        if(tjl.empty()) continue;
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<"Plane "<<plane<<" tjList";
          for(auto& tjID : tjl) myprt<<" "<<tjID;
        }
        int cotID = Create2DShower(tjs, tjl);
        if(cotID == 0) continue;
        if(!DefineShower(fcnLabel, tjs, cotID, prt)) {
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failure from DefineShower "<<cotID<<" failed ";
          MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
          continue;
        }
        SaveTjInfo(tjs, inCTP, cotID, "DS");
        // Find nearby Tjs that were not included because they had too-high
        // MCSMom, etc. This will be used to decide if showers should be merged
        AddTjsInsideEnvelope(fcnLabel, tjs, cotID, prt);
        FindNearbyTjs(fcnLabel, tjs, cotID, prt);
      } // tjl
      // try to merge showers in this plane using the lists of nearby Tjs
      if(inCTP == UINT_MAX) continue;
      if(tjs.cots.empty()) continue;
      prt = (inCTP == dbgCTP || dbgPlane == 3);
      if(prt) Print2DShowers("tjl", tjs, inCTP, false);
      MergeShowerChain(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "MSC");
      if(prt) Print2DShowers("MSCo", tjs, inCTP, false);
      MergeOverlap(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "MO");
      if(prt) Print2DShowers("MO", tjs, inCTP, false);
      SaveAllCots(tjs, inCTP, "MSS");
      MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);
      MergeSubShowers(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "MNrby");
      if(prt) Print2DShowers("Nrby", tjs, inCTP, false);
      for(unsigned short ii = 0; ii < tjs.cots.size(); ++ii) {
        auto& ss = tjs.cots[ii];
        if(ss.ID == 0) continue;
        if(ss.CTP != inCTP) continue;
        if(AddTjsInsideEnvelope(fcnLabel, tjs, ss.ID, prt)) MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);
        if (tjs.SaveShowerTree) SaveAllCots(tjs, inCTP, "Merge");
      }
      SaveAllCots(tjs, inCTP, "ATj2");
      if(prt) Print2DShowers("ATIE", tjs, inCTP, false);
    } // plane
    if(tjs.cots.empty()) return false;
    prt = (dbgPlane > 2);
    if(prt) Print2DShowers("B4", tjs, USHRT_MAX, false);
    // Match in 3D, make 3D showers and define them
    Match2DShowers(fcnLabel, tjs, tpcid, prt);
    if(prt) Print2DShowers("M2DS", tjs, USHRT_MAX, false);
    // Kill vertices in 2D showers that weren't matched in 3D
//    KillVerticesInShowers(fcnLabel, tjs, tpcid, prt);
    SaveAllCots(tjs, "RS");
     
    unsigned short nNewShowers = 0;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      SaveTjInfo(tjs, ss.CTP, ss.ID, "Done");
     ++nNewShowers;
    } // ss
    
//    if(prt) Print2DShowers("FSo", tjs, USHRT_MAX, false);
    
    return (nNewShowers > 0);
    
  } // FindShowers3D

  ////////////////////////////////////////////////
  void CheckInShowerProb(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // Checks the InShower likelihood for all PFParticles in the shower
    if(ss3.PFPIDs.empty()) return;

    std::string fcnLabel = inFcnLabel + ".CISP";

    for(auto id : ss3.PFPIDs) {
      auto& pfp = tjs.pfps[id - 1];
      double inShowerProb = InShowerProb(fcnLabel, tjs, ss3, pfp);
      if(prt) mf::LogVerbatim("TC")<<" 3S"<<ss3.ID<<" P"<<pfp.ID<<" prob "<<inShowerProb;
    } // id
  } // CheckInShowerProb
  
  ////////////////////////////////////////////////
  void FindInShowerPFPs(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, std::vector<std::vector<unsigned short>>& plists)
  {
    // 
    
    plists.clear();
    if(tjs.pfps.size() < 2) return;
    
    std::string fcnLabel = inFcnLabel + ".FISP";
    
    struct ClosePair {
      float doca;
      unsigned short id1;
      unsigned short close1;
      unsigned short id2;
      unsigned short close2;
      bool used;
    };

    // min separation between pfps
    std::vector<ClosePair> cps;
    // list of pfps that are close
    std::vector<unsigned short> pfpList;
    for(unsigned short ip = 0; ip < tjs.pfps.size() - 1; ++ip) {
      auto& p1 = tjs.pfps[ip];
      if(p1.ID == 0) continue;
      if(p1.TPCID != tpcid) continue;
      // ignore neutrinos
      if(p1.PDGCode == 14 || p1.PDGCode == 12) continue;
      if(p1.TjIDs.empty()) continue;
      bool p1TrackLike = (MCSMom(tjs, p1.TjIDs) > tjs.ShowerTag[1]);
      for(unsigned short jp = ip + 1; jp < tjs.pfps.size(); ++jp) {
        auto& p2 = tjs.pfps[jp];
        if(p2.ID == 0) continue;
        if(p2.TPCID != tpcid) continue;
        // ignore neutrinos
        if(p2.PDGCode == 14 || p2.PDGCode == 12) continue;
        if(p2.TjIDs.empty()) continue;
        bool p2TrackLike = (MCSMom(tjs, p2.TjIDs) > tjs.ShowerTag[1]);
        // require one of them to be not tracklike
        if(p1TrackLike && p2TrackLike) continue;
        unsigned short close1, close2;
        float doca = PFPDOCA(p1, p2, close1, close2);
        std::cout<<fcnLabel<<" P"<<p1.ID<<" TrackLike? "<<p1TrackLike<<" P"<<p2.ID<<" TrackLike? "<<p2TrackLike;
        std::cout<<" doca "<<std::fixed<<std::setprecision(1)<<doca;
        std::cout<<" dang "<<DeltaAngle(p1.Dir[0], p2.Dir[0]);
        std::cout<<"\n";
        if(doca > tjs.ShowerTag[2]) continue;
        // add a new one
        ClosePair cp;
        cp.doca = doca;
        cp.id1 = p1.ID;
        cp.close1 = close1;
        cp.id2 = p2.ID;
        cp.close2 = close2;
        cp.used = false;
        cps.push_back(cp);
        if(std::find(pfpList.begin(), pfpList.end(), p1.ID) == pfpList.end()) pfpList.push_back(p1.ID);
        if(std::find(pfpList.begin(), pfpList.end(), p2.ID) == pfpList.end()) pfpList.push_back(p2.ID);
      } // jp
    } // ip

    if(cps.empty()) return;
    // sort pfpList by decreasing length
    std::vector<SortEntry> sortVec(pfpList.size());
    for(unsigned short ii = 0; ii < pfpList.size(); ++ii) {
      sortVec[ii].index = ii;
      auto& pfp = tjs.pfps[pfpList[ii] - 1];
      sortVec[ii].length = PosSep(pfp.XYZ[0], pfp.XYZ[1]);
    } // ii
    std::sort(sortVec.begin(), sortVec.end(), greaterThan);
    auto tmp = pfpList;
    for(unsigned short ii = 0; ii < pfpList.size(); ++ii) pfpList[ii] = tmp[sortVec[ii].index];
  
    // cluster pfparticles starting with the longest
    for(unsigned short ii = 0; ii < pfpList.size(); ++ii) {
      int pid1 = pfpList[ii];
      auto& p1 = tjs.pfps[pid1 - 1];
      float p1Len = PosSep(p1.XYZ[0], p1.XYZ[1]);
      std::cout<<ii<<" P"<<pid1<<" len "<<p1Len<<"\n";
      std::vector<unsigned short> plist;
      plist.push_back(pid1);
      // try to add ids to the list
      bool added = true;
      while(added) {
        added = false;
        for(auto& cp : cps) {
          if(cp.used) continue;
          bool isID1 = (std::find(plist.begin(), plist.end(), cp.id1) != plist.end());
          bool isID2 = (std::find(plist.begin(), plist.end(), cp.id2) != plist.end());
          if(!(isID1 || isID2)) continue;
          unsigned short pid2 = cp.id1;
          if(isID1) pid2 = cp.id2;
          auto& p2 = tjs.pfps[pid2 - 1];
          float p2Len = PosSep(p2.XYZ[0], p2.XYZ[1]);
          float dang = DeltaAngle(p1.Dir[0], p2.Dir[0]);
          // don't cluster two long pfparticles with a large angle difference
          if(p1Len > 5 && p2Len > 5 && cp.doca < 2 && dang > 0.3) continue;
          // don't cluster if this pfparticle is closer to a different long pfparticle
          bool isCloser = false;
          for(auto& pcp : cps) {
            if(pid1 == pcp.id1 || pid1 == pcp.id2) continue;
            if(!(pid2 == pcp.id1 || pid2 == pcp.id2)) continue;
            unsigned short oid = pcp.id1;
            if(oid == pid2) oid = pcp.id2;
            auto opfp = tjs.pfps[oid - 1];
            float pcpLen = PosSep(opfp.XYZ[0], opfp.XYZ[1]);
            std::cout<<"pid1 P"<<pid1<<" pid2 P"<<pid2<<" oid P"<<oid<<"\n";
            if(pcp.doca < cp.doca && pcpLen > 5) isCloser = true;
          } // kk
          if(isCloser) continue;
          if(std::find(plist.begin(), plist.end(), pid2) != plist.end()) continue;
          std::cout<<"  add P"<<pid2<<" doca "<<cp.doca<<" dang "<<dang<<"\n";
          plist.push_back(pid2);
          // call it used
          cp.used = true;
          added = true;
        } // cp
      } // added
      if(plist.size() > 1) plists.push_back(plist);
    } // ii
    
    // Check for leftover pfps and add them if they are shower-like (using the Tj InShower tag)
    std::vector<unsigned short> flat;
    for(auto& plist : plists) flat.insert(flat.end(), plist.begin(), plist.end());
    auto notClustered = SetDifference(pfpList, flat);
    for(auto nc : notClustered) {
      auto& pfp = tjs.pfps[nc - 1];
      if(!IsInShower(tjs, pfp.TjIDs)) continue;
      std::vector<unsigned short> plist(1, nc);
      plists.push_back(plist);
    }

    mf::LogVerbatim myprt("TC");
    myprt<<fcnLabel<<" plists size "<<plists.size()<<"\n";
    for(unsigned short ip = 0; ip < plists.size(); ++ip) {
      auto& plist = plists[ip];
      myprt<<"ip "<<ip;
      for(auto pid : plist) {
        myprt<<" P"<<pid;
/*
        auto& pfp = tjs.pfps[pid - 1];
        float plen = PosSep(pfp.XYZ[0], pfp.XYZ[1]);
        myprt<<"_"<<std::fixed<<std::setprecision(1)<<plen;
*/
      }
      myprt<<"\n";
    } // plist

  } // FindInShowerPFPs

  ////////////////////////////////////////////////
  void KillVerticesInShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // Kill vertices inside showers 
    
    std::string fcnLabel = inFcnLabel + ".KVIS";
    
    // first kill the vertices
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      unsigned short parentTjVxID = stj.VtxID[0];
      // look for vertices inside the envelope
      for(auto& vx2 : tjs.vtx) {
        if(vx2.CTP != ss.CTP) continue;
        if(vx2.ID == 0) continue;
        if(vx2.ID == parentTjVxID) continue;
        // make sure it isn't associated with the neutrino vertex
        if(vx2.Vx3ID > 0) {
          auto& vx3 = tjs.vtx3[vx2.Vx3ID - 1];
          if(vx3.Neutrino) continue;
        }
        bool insideEnvelope = PointInsideEnvelope(vx2.Pos, ss.Envelope);
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" ss "<<ss.ID<<" vx2 "<<vx2.ID<<" Score "<<std::fixed<<std::setprecision(1)<<vx2.Score<<" insideEnvelope? "<<insideEnvelope;
        if(!insideEnvelope) continue;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Clobber vtx "<<vx2.ID;
        MakeVertexObsolete(tjs, vx2, true);
        ss.NeedsUpdate = true;
      } // vx2
/*
      // check InShower tjs next
      for(auto tjid : ss.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == 0) continue;
          if(tj.VtxID[end] == parentTjVxID) continue;
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Clobber vx "<<tj.VtxID[end]<<" on tj "<<tj.ID;
          auto& vx = tjs.vtx[tj.VtxID[end] - 1];
          MakeVertexObsolete(tjs, vx, true);
          ss.NeedsUpdate = true;
        } // end
      } // tj
*/
    } // ss
/*
    // reset the scores
    ScoreVertices(tjs, tpcid, prt);
    
    // update the showers
    for(auto& ss : tjs.cots) {
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      if(!ss.NeedsUpdate) continue;
      DefineShower(fcnLabel, tjs, ss.ID - 1, prt);
      // Add Tjs with high-score vertices inside the shower and kill those vertices
      AddTjsInsideEnvelope(fcnLabel, tjs, ss.ID - 1, prt);
      if(ss.NeedsUpdate) DefineShower(fcnLabel, tjs, ss.ID - 1, prt);
    } // ci
*/
  } // KillVerticesInShowers
  
  ////////////////////////////////////////////////
  void Match2DShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    if(prt) mf::LogVerbatim("TC")<<"Inside M2DS";
    // Match 2D showers using position and direction to create 3D showers
    
    std::string fcnLabel = inFcnLabel + ".M2DS";
    
    float fomCut = 2;
    
    // sort the showers by decreasing energy and increasing AspectRatio so that the 3D direction is defined
    // by the first matching pair
    std::vector<SortEntry> sortVec;
    for(unsigned short indx = 0; indx < tjs.cots.size(); ++indx) {
      auto& ss = tjs.cots[indx];
      if(ss.ID == 0) continue;
      // already matched?
      if(ss.SS3ID > 0) continue;
      if(ss.TjIDs.empty()) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      SortEntry se;
      se.index = indx;
      se.length = ss.Energy / ss.AspectRatio;
      sortVec.push_back(se);
    } // indx
    if(sortVec.size() < 2) return;
    std::sort(sortVec.begin(), sortVec.end(), greaterThan);
    
    // Look for a 3D match using the 2D shower charge centers
    for(unsigned short ii = 0; ii < sortVec.size() - 1; ++ii) {
      unsigned short iIndx = sortVec[ii].index;
      auto& iss = tjs.cots[iIndx];
      // already matched?
      if(iss.SS3ID > 0) continue;
      Trajectory& istj = tjs.allTraj[iss.ShowerTjID - 1];
      geo::PlaneID iplaneID = DecodeCTP(iss.CTP);
      for(unsigned short jj = 0; jj < sortVec.size(); ++jj) {
        unsigned short jIndx = sortVec[jj].index;
        ShowerStruct& jss = tjs.cots[jIndx];
        // already matched?
        if(jss.SS3ID > 0) continue;
        if(jss.CTP == iss.CTP) continue;
        Trajectory& jstj = tjs.allTraj[jss.ShowerTjID - 1];
        TrajPoint3 tp3;
        if(!MakeTp3(tjs, istj.Pts[1], jstj.Pts[1], tp3, true)) {
//          if(prt) mf::LogVerbatim("TC")<<" MakeTp3s failed";
          continue;
        }
        float fomij = Match3DFOM(fcnLabel, tjs, iss.ID, jss.ID, prt);
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" i2S"<<iss.ID<<" j2S"<<jss.ID<<" fomij "<<fomij<<" fomCut "<<fomCut;
        if(fomij > fomCut) continue;
        geo::PlaneID jplaneID = DecodeCTP(jss.CTP);
        if(tjs.NumPlanes == 2) {
          ShowerStruct3D ss3 = CreateSS3(tjs, tpcid);
          ss3.ChgPos = tp3.Pos;
          ss3.Dir = tp3.Dir;
          ss3.CotIDs.resize(2);
          ss3.CotIDs[0] = iss.ID;
          ss3.CotIDs[1] = jss.ID;
          ss3.Energy.resize(2);
          ss3.Energy[0] = iss.Energy;
          ss3.Energy[1] = jss.Energy;
          ss3.FOM = fomij;
          ss3.PFPIndex = USHRT_MAX;
          // don't fill or use the rest of the variables
          tjs.showers.push_back(ss3);
          iss.SS3ID = ss3.ID;
          jss.SS3ID = ss3.ID;
          if(prt) mf::LogVerbatim("TC")<<" new 3S"<<ss3.ID<<" with fomij "<<fomij;
          continue;
        } // 2-plane TPC
        float bestFOM = fomCut;
        unsigned short bestck = USHRT_MAX;
        for(unsigned short ck = 0; ck < tjs.cots.size(); ++ck) {
          ShowerStruct& kss = tjs.cots[ck];
          if(kss.ID == iss.ID || kss.ID == jss.ID) continue;
          if(kss.CTP == iss.CTP || kss.CTP == jss.CTP) continue;
          if(kss.ID == 0) continue;
          if(kss.TjIDs.empty()) continue;
          if(kss.SS3ID > 0) continue;
          geo::PlaneID kplaneID = DecodeCTP(kss.CTP);
          if(kplaneID.Cryostat != tpcid.Cryostat) continue;
          if(kplaneID.TPC != tpcid.TPC) continue;
          Trajectory& kstj = tjs.allTraj[kss.ShowerTjID - 1];
          TrajPoint3 iktp3;
          MakeTp3(tjs, istj.Pts[1], kstj.Pts[1], iktp3, true);
          float fomik = Match3DFOM(fcnLabel, tjs, iss.ID, kss.ID, prt);
          if(fomik > bestFOM) continue;
          float sep = PosSep(tp3.Pos, iktp3.Pos);
          if(sep > 50) {
            if(prt) mf::LogVerbatim("TC")<<" Large stp[1] point separation "<<sep;
            continue;
          }
          bestFOM = fomik;
          bestck = ck;
        } // ck
        // 3-plane TPC below
        ShowerStruct3D ss3 = CreateSS3(tjs, tpcid);
        // Define ss3 using the tp3 found with the first pair
        ss3.ChgPos = tp3.Pos;
        ss3.Dir = tp3.Dir;
        iss.SS3ID = ss3.ID;
        jss.SS3ID = ss3.ID;
        ss3.FOM = bestFOM;
        if(bestck == USHRT_MAX) {
          // showers match in 2 planes
          ss3.CotIDs.resize(2);
          ss3.CotIDs[0] = iss.ID;
          ss3.CotIDs[1] = jss.ID;
          ss3.Energy[iplaneID.Plane] = iss.Energy;
          ss3.Energy[jplaneID.Plane] = jss.Energy;
          if(prt) mf::LogVerbatim("TC")<<" new 2-plane 3S"<<ss3.ID<<" using 2S"<<iss.ID<<" 2S"<<jss.ID<<" with FOM "<<ss3.FOM;
        } else {
          // showers match in 3 planes
          unsigned short ck = bestck;
          ShowerStruct& kss = tjs.cots[ck];
          ss3.CotIDs.resize(3);
          ss3.CotIDs[0] = iss.ID;
          ss3.CotIDs[1] = jss.ID;
          ss3.CotIDs[2] = kss.ID;
          geo::PlaneID kplaneID = DecodeCTP(kss.CTP);
          ss3.Energy[iplaneID.Plane] = iss.Energy;
          ss3.Energy[jplaneID.Plane] = jss.Energy;
          ss3.Energy[kplaneID.Plane] = kss.Energy;
          tjs.cots[ck].SS3ID = ss3.ID;
          if(prt) mf::LogVerbatim("TC")<<" new 3-plane 3S"<<ss3.ID<<" using 2S"<<iss.ID<<" 2S"<<jss.ID<<" 2S"<<tjs.cots[ck].ID<<" with FOM "<<ss3.FOM;
        }
        ss3.FOM = 0.5 * (fomij + bestFOM);
        // sort the IDs
        std::sort(ss3.CotIDs.begin(), ss3.CotIDs.end());
        // look for a 3D shower parent
        FindParent(fcnLabel, tjs, ss3, prt);
        // check for a failure
        if(ss3.ID == 0) continue;
        if(ss3.NeedsUpdate) UpdateShower(fcnLabel, tjs, ss3, prt);
        tjs.showers.push_back(ss3);
      } // cj
    } // ci
    if(tjs.showers.empty()) return;
    
    if(prt) PrintShowers("M2DS", tjs);

  } // Match2DShowers

  ////////////////////////////////////////////////
  void UpdateShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // Updates the 3D shower presumably because the 2D showers were changed or need to be updated. 
    
    if(ss3.ID == 0) return;
    if(ss3.CotIDs.size() < 2) return;
    
    std::string fcnLabel = inFcnLabel + ".U3S";
    
    // see if the 2D showers need an update
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      if(ss.NeedsUpdate) std::cout<<fcnLabel<<" ********* 3S"<<ss3.ID<<" 2S"<<ss.ID<<" needs an update...\n";
    } // ci
    
    // check consistency 
    if(ss3.ParentID > 0) {
      auto& pfp = tjs.pfps[ss3.ParentID - 1];
      unsigned short pend = FarEnd(tjs, pfp, ss3.ChgPos);
      if(pfp.Vx3ID[pend] != ss3.Vx3ID) {
        std::cout<<fcnLabel<<" ********* 3S"<<ss3.ID<<" has parent P"<<ss3.ParentID<<" with a vertex that is not attached the shower\n";
      }
    } // ss3.ParentID > 0
    
    // Find the average position and direction using pairs of 2D shower Tjs
    std::array<Point3_t, 3> pos;
    // the direction of all points in 2D showers is the same
    Vector3_t dir;
    std::array<double, 3> chg;
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      chg[ipt] = 0;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) {
        pos[ipt][xyz] = 0;
        dir[xyz] = 0;
      }
    } // ipt
    
    for(unsigned short ii = 0; ii < ss3.CotIDs.size() - 1; ++ii) {
      unsigned short ciid = ss3.CotIDs[ii];
      auto& iss = tjs.cots[ciid - 1];
      // make a local copy of the trajectory so we can replace Pos with HitPos = charge center
      auto istj = tjs.allTraj[iss.ShowerTjID - 1];
      for(auto& tp : istj.Pts) tp.Pos = tp.HitPos;
      for(unsigned short jj = ii + 1; jj < ss3.CotIDs.size(); ++jj) {
        unsigned short cjid = ss3.CotIDs[jj];
        auto& jss = tjs.cots[cjid - 1];
        auto jstj = tjs.allTraj[jss.ShowerTjID - 1];
        for(auto& tp : jstj.Pts) tp.Pos = tp.HitPos;
        TrajPoint3 tp3;
        for(unsigned short ipt = 0; ipt < 3; ++ipt) {
          if(!MakeTp3(tjs, istj.Pts[ipt], jstj.Pts[ipt], tp3, true)) {
            if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" MakeTp3s failure for 3S"<<ss3.ID<<". Killing it";
            ss3.ID = 0;
            return;
          }
          double ptchg = 0.5 * (istj.Pts[ipt].Chg + jstj.Pts[ipt].Chg);
          chg[ipt] += ptchg;
          for(unsigned short xyz = 0; xyz < 3; ++xyz) {
            pos[ipt][xyz] += ptchg * tp3.Pos[xyz];
            dir[xyz] += ptchg * tp3.Dir[xyz];
          } // xyz
        } // ipt
      } // jj
    } // ii
    
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      if(chg[ipt] == 0) continue;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pos[ipt][xyz] /= chg[ipt];
      SetMag(dir, 1);
    } // ipt
    ss3.ChgPos = pos[1];
    
    if(ss3.ParentID > 0) {
      // There is a 3D-matched pfp at the shower start. The end that is farthest away from the
      // shower center should be shower start
      auto& pfp = tjs.pfps[ss3.ParentID - 1];
      unsigned short pend = FarEnd(tjs, pfp, ss3.ChgPos);
      ss3.Start = pfp.XYZ[pend];
      // TODO: use charge weighted average of shower direction and pfp direction?
      ss3.Dir = dir;
    } else {
      ss3.Dir = dir;
      ss3.Start = pos[0];
    }
    // define the end
    ss3.End = pos[2];
    ss3.Len = PosSep(ss3.Start, ss3.End);
    
    // dE/dx, energy, etc
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid];
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      unsigned short plane = DecodeCTP(ss.CTP).Plane;
      ss3.Energy[plane] = ss.Energy;
      // TODO: calculate the errors in some better way
      ss3.EnergyErr[plane] = 0.3 * ss.Energy;
      // TODO: what does MIPEnergy mean anyway?
      ss3.MIPEnergy[plane] = ss3.EnergyErr[plane];
      ss3.MIPEnergyErr[plane] = ss.Energy;
      ss3.dEdx[plane] = stj.dEdx[0];
      ss3.dEdxErr[plane] = 0.3 * stj.dEdx[0];
    } // ci
             
    ss3.NeedsUpdate = false;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" updated";

  } // UpdateShower

  ////////////////////////////////////////////////
  void ReconcileParents(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // Ensures that a TJ is only used as a Parent in one shower
    
    std::string fcnLabel = inFcnLabel + ".RP";

    std::vector<int> parent;
    std::vector<std::vector<int>> inss;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.ParentID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      unsigned short indx = 0;
      for(indx = 0; indx < parent.size(); ++indx) if(ss.ParentID == parent[indx]) break;
      if(indx == parent.size()) {
        // not found add it
        parent.push_back(ss.ParentID);
        std::vector<int> tmp(1);
        tmp[0] = ss.ID;
        inss.push_back(tmp);
      } else {
        inss[indx].push_back(ss.ID);
      }
    } // ss
    if(parent.empty()) return;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      for(unsigned short pid = 0; pid < parent.size(); ++pid) {
        myprt<<fcnLabel<<" parent "<<parent[pid]<<" ssids";
        for(auto ssid : inss[pid]) myprt<<" "<<ssid;
        myprt<<"\n";
      }  // pid
    } // prt

    for(unsigned short pid = 0; pid < parent.size(); ++pid) {
      if(inss[pid].size() < 2) continue;
      std::cout<<fcnLabel<<" Tj "<<parent[pid]<<" is a parent in more than one shower";
      for(auto ssid : inss[pid]) std::cout<<" "<<ssid;
      std::cout<<". Do something about this.\n";
    } // pid

  } // ReconcileParents

  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    float fom = 0;
    float cnt = 0;
    for(unsigned short ii = 0; ii < ss3.CotIDs.size() - 1; ++ii) {
      unsigned short icid = ss3.CotIDs[ii];
      for(unsigned short jj = ii + 1; jj < ss3.CotIDs.size(); ++jj) {
        unsigned short jcid = ss3.CotIDs[jj];
        fom += Match3DFOM(inFcnLabel, tjs, icid, jcid, prt);
        ++cnt;
      } // cj
    } // ci
    if(cnt == 0) return 100;
    return fom / cnt;
  } // Match3DFOM

  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs,
                   int icid, int jcid, int kcid, bool prt)
  {
    if(icid == 0 || icid > tjs.cots.size()) return 100;
    if(jcid == 0 || jcid > tjs.cots.size()) return 100;
    if(kcid == 0 || kcid > tjs.cots.size()) return 100;
    
    float ijfom = Match3DFOM(inFcnLabel, tjs, icid, jcid, prt);
    float jkfom = Match3DFOM(inFcnLabel, tjs, jcid, kcid, prt);
    
    return 0.5 * (ijfom + jkfom);
    
  } // Match3DFOM
  
  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, int icid, int jcid, bool prt)
  {
    // returns a Figure of Merit for a 3D match of two showers
    if(icid == 0 || icid > tjs.cots.size()) return 100;
    if(jcid == 0 || jcid > tjs.cots.size()) return 100;
    
    auto& iss = tjs.cots[icid - 1];
    auto& istj = tjs.allTraj[iss.ShowerTjID - 1];    
    auto& jss = tjs.cots[jcid - 1];
    auto& jstj = tjs.allTraj[jss.ShowerTjID - 1];
    
    if(iss.CTP == jss.CTP) return 100;
    
    std::string fcnLabel = inFcnLabel + ".MFOM";
    
    float energyAsym = std::abs(iss.Energy - jss.Energy) / (iss.Energy + jss.Energy);
    
    if(energyAsym > 0.5) return 50;
    
    geo::PlaneID iPlnID = DecodeCTP(iss.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jss.CTP);

    // compare match at the charge center
    float ix = tjs.detprop->ConvertTicksToX(istj.Pts[1].Pos[1] / tjs.UnitsPerTick, iPlnID);
    float jx = tjs.detprop->ConvertTicksToX(jstj.Pts[1].Pos[1] / tjs.UnitsPerTick, jPlnID);
    float pos1fom = std::abs(ix - jx) / 10;

    float mfom = energyAsym * pos1fom;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" i2S"<<iss.ID<<" j2S"<<jss.ID;
      myprt<<std::fixed<<std::setprecision(2);
//      myprt<<" pos0fom "<<pos0fom;
      myprt<<" pos1fom "<<pos1fom;
//      myprt<<" pos2fom "<<pos2fom;
      myprt<<" energyAsym "<<energyAsym;
      myprt<<" mfom "<<mfom;
    }
    
    return mfom;
  } // Madtch3DFOM

  ////////////////////////////////////////////////
  void MergeTjList(std::vector<std::vector<int>>& tjList)
  {
    // Merge the lists of Tjs in the lists if they share a common Tj ID
    
    if(tjList.size() < 2) return;
    
    bool didMerge = true;
    while(didMerge) {
      didMerge = false;
      for(unsigned short itl = 0; itl < tjList.size() - 1; ++itl) {
        if(tjList[itl].empty()) continue;
        for(unsigned short jtl = itl + 1; jtl < tjList.size(); ++jtl) {
          if(tjList[itl].empty()) continue;
          auto& itList = tjList[itl];
          auto& jtList = tjList[jtl];
          // See if the j Tj is in the i tjList
          bool jtjInItjList = false;
          for(auto& jtj : jtList) {
            if(std::find(itList.begin(), itList.end(), jtj) != itList.end()) {
              jtjInItjList = true;
              break;
            }
            if(jtjInItjList) break;
          } // jtj
          if(jtjInItjList) {
            // append the jtList to itList
            itList.insert(itList.end(), jtList.begin(), jtList.end());
            // clear jtList
            jtList.clear();
            didMerge = true;
          }
        } // jtl
      } // itl
    } // didMerge
    
    // erase the deleted elements
    unsigned short imEmpty = 0;
    while(imEmpty < tjList.size()) {
      for(imEmpty = 0; imEmpty < tjList.size(); ++imEmpty) if(tjList[imEmpty].empty()) break;
      if(imEmpty < tjList.size()) tjList.erase(tjList.begin() + imEmpty);
    } // imEmpty < tjList.size()
    
    // sort the lists by increasing ID and remove duplicates
    for(auto& tjl : tjList) {
      std::sort(tjl.begin(), tjl.end());
      auto last = std::unique(tjl.begin(), tjl.end());
      tjl.erase(last, tjl.end());
    } // tjl
    
  } // MergeTjList

  ////////////////////////////////////////////////
  void FillPts(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    if(ss.ShowerTjID == 0) return;
    ss.ShPts.clear();
    
    std::string fcnLabel = inFcnLabel + ".FPts";
    
    float totChg = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      if(itj > tjs.allTraj.size() - 1) {
        mf::LogWarning("TC")<<"Bad TjID "<<ss.TjIDs[it];
        MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
        return;
      }
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.CTP != ss.CTP) {
        mf::LogWarning("TC")<<"Tj "<<tj.ID<<" is in the wrong CTP "<<tj.CTP<<" "<<ss.CTP;
        MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
        return;
      }
      if(tj.AlgMod[kShowerTj]) {
        mf::LogWarning("TC")<<fcnLabel<<" Tj "<<tj.ID<<" is in TjIDs in 2S"<<cotID<<" but is a ShowerTj! Killing it";
        MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
        return;
      }
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(!tp.UseHit[ii]) continue;
          unsigned int iht = tp.Hits[ii];
          if(tjs.fHits[iht].Integral <= 0) continue;
          ShowerPoint shpt;
          shpt.HitIndex = iht;
          shpt.TID = tj.ID;
          shpt.Chg = tjs.fHits[iht].Integral;
          shpt.Pos[0] = tjs.fHits[iht].ArtPtr->WireID().Wire;
          shpt.Pos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
          ss.ShPts.push_back(shpt);
          totChg += shpt.Chg;
        } // ii
      } // ipt
    } // it

    // Put the total charge into the shower Tj
    tjs.allTraj[ss.ShowerTjID - 1].AveChg = totChg;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 2S"<<ss.ID<<" filled "<<ss.ShPts.size()<<" points. Total charge "<<(int)totChg;
      for(auto tid : ss.TjIDs) myprt<<" T"<<tid;
    } // prt
  } // FillPts

  ////////////////////////////////////////////////
  bool DefinePFPShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // defines the 3D shower position, direction, length and opening angle. The PFParticles referenced
    // in plist were fit to a 3D line in PFPUtils/FindCompleteness. The origin of the fit is at the
    // center of matching points of the trajectories. The idea here is to first find the charge weighted position
    // of all these points to find the charge center and the direction. The next step is to find the longitudinal
    // extent and then the transverse width at each end. The end that has the smallest transverse width is declared
    // to be the start
    if(ss3.PFPIDs.empty()) return false;
    
    // Find the charge center of all pfps
    ss3.Energy.resize(tjs.NumPlanes);
    ss3.EnergyErr.resize(tjs.NumPlanes);
    ss3.MIPEnergy.resize(tjs.NumPlanes);
    ss3.MIPEnergyErr.resize(tjs.NumPlanes);
    ss3.dEdx.resize(tjs.NumPlanes);
    ss3.dEdxErr.resize(tjs.NumPlanes);
    
    std::string fcnLabel = inFcnLabel + ".DPS";
    
    // put all the tjs into a monster shower pfparticle so we can get the
    // direction and position. Sum the Tj charge and convert to energy while we are here
    PFPStruct shPFP = CreatePFP(tjs, ss3.TPCID);
    std::vector<float> chgSum(tjs.NumPlanes);
    for(auto pid : ss3.PFPIDs) {
      auto& pfp = tjs.pfps[pid - 1];
      shPFP.TjIDs.insert(shPFP.TjIDs.end(), pfp.TjIDs.begin(), pfp.TjIDs.end());
      shPFP.Tp3s.insert(shPFP.Tp3s.end(), pfp.Tp3s.begin(), pfp.Tp3s.end());
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        unsigned short plane = DecodeCTP(tj.CTP).Plane;
        chgSum[plane] += tj.TotChg;
      } // tjid
    } // pid
    if(shPFP.Tp3s.empty()) return false;
    // Initialize and pass dummy variables
    Point3_t dump;
    Vector3_t dumd;
    Fit3D(0, dump, dumd, dump, dumd);
    // Fill the fit sums
    for(auto& tp3 : shPFP.Tp3s) {
      Fit3D(1, tp3.Pos, tp3.Dir, dump, dumd);
    } // tp3
    Fit3D(2, dump, dumd, shPFP.XYZ[0], shPFP.Dir[0]);
    for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) ss3.Energy[plane] = ChgToMeV(chgSum[plane]);
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 3S"<<ss3.ID<<" shPFP fit nTp3s "<<shPFP.Tp3s.size();
      myprt<<" Energy[plane]";
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) myprt<<" "<<(int)ss3.Energy[plane];
      myprt<<" chgCtr "<<std::fixed<<std::setprecision(1)<<shPFP.XYZ[0][0]<<" "<<shPFP.XYZ[0][1]<<" "<<shPFP.XYZ[0][2];
      myprt<<" dir "<<std::setprecision(1)<<shPFP.Dir[0][0]<<" "<<shPFP.Dir[0][1]<<" "<<shPFP.Dir[0][2];
    } // prt
    ss3.ChgPos = shPFP.XYZ[0];
    ss3.Dir = shPFP.Dir[0];
    // find the points that are farthest away from the charge center on the shower axis
    double minAlong = 1E6;
    double maxAlong = -1E6;
    for(auto& tp3 : shPFP.Tp3s) {
      FindAlongTrans(ss3.ChgPos, ss3.Dir, tp3.Pos, tp3.AlongTrans);
      if(tp3.AlongTrans[0] < minAlong) minAlong = tp3.AlongTrans[0];
      if(tp3.AlongTrans[0] > maxAlong) maxAlong = tp3.AlongTrans[0];
    } // tp3
    // set the start and end positions
    for(unsigned short xyz = 0; xyz < 3; ++xyz) {
      ss3.Start[xyz] = ss3.ChgPos[xyz] + minAlong * ss3.Dir[xyz];
      ss3.End[xyz] = ss3.ChgPos[xyz] + maxAlong * ss3.Dir[xyz];
    }
    ss3.Len = maxAlong - minAlong;
    if(prt) mf::LogVerbatim("TC")<<" minAlong "<<minAlong<<" maxAlong "<<maxAlong<<" Len "<<ss3.Len;
    std::vector<double> chg(3);
    std::vector<double> sumt(3);
    std::vector<double> suml(3);
    for(auto& tp3 : shPFP.Tp3s) {
      unsigned short section = 3 * (tp3.AlongTrans[0] - minAlong) / ss3.Len;
      if(section > 2) section = 2;
      chg[section] += tp3.dEdx;
      suml[section] += tp3.dEdx * tp3.AlongTrans[0];
      sumt[section] += tp3.dEdx * tp3.AlongTrans[1];
    } // tp3
    // something is wrong if there is no charge in any section
    if(chg[0] == 0 || chg[1] == 0 || chg[2] == 0) return false;
    for(unsigned short section = 0; section < 3; ++section) {
      if(chg[section] <= 0) continue;
      sumt[section] /= chg[section];
      suml[section] /= chg[section];
      if(prt) mf::LogVerbatim("TC")<<" section "<<section<<" long "<<suml[section]<<" trans "<<sumt[section]<<" chg "<<(int)chg[section];
    } // section
    // reverse the direction?
    if(sumt[2] < sumt[0]) {
      std::swap(ss3.Start, ss3.End);
      for(unsigned short xyz = 0; xyz < 3; ++xyz) ss3.Dir[xyz] *= -1;
      // we don't need the Tp3s anymore so don't reverse them
      std::swap(sumt[0], sumt[2]);
      std::swap(suml[0], suml[2]);
      std::swap(chg[0], chg[2]);
    }
    // calculate the opening angle - Note that this is 1/2 the tangent of 
    // the opening angle...
    ss3.OpenAngle = std::abs((sumt[2] - sumt[0]) / (suml[2] - suml[0]));
    if(ss3.OpenAngle < 0.1) {
      std::cout<<"DefinePFPShower: Too small opening angle "<<ss3.OpenAngle<<". Setting it to 0.1\n";
      ss3.OpenAngle = 0.1;
    }
    if(ss3.OpenAngle > 0.3) {
      std::cout<<"DefinePFPShower: Too large opening angle "<<ss3.OpenAngle<<". Setting it to 0.3\n";
      ss3.OpenAngle = 0.3;
    }
    // check the direction vector using a largish component
    unsigned short useComp = 0;
    if(std::abs(ss3.Dir[2]) > std::abs(ss3.Dir[0])) useComp = 2;
    if(ss3.Dir[useComp] > 0 && ss3.Start[useComp] > ss3.End[useComp]) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) ss3.Dir[xyz] *= -1;
    } // reverse direction?
    // stash the charge in each section into PosErr. It will be transferred into
    // the shower Tj TP charge in MakePFPShowers. Put the "rms" into DirErr. Note that
    // the "rms" is really the average of the abs(transverse) positions. This is close
    // enough to the actual rms for this purpose.
    for(unsigned short section = 0; section < 3; ++section) {
      ss3.StartErr[section] = chg[section];
      ss3.DirErr[section] = sumt[section];
    }

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 3S"<<ss3.ID<<" OpenAngle "<<std::setprecision(2)<<ss3.OpenAngle;
      myprt<<" Dir "<<std::setprecision(2)<<ss3.Dir[0]<<" "<<ss3.Dir[1]<<" "<<ss3.Dir[2];
      float aveE = 0;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) aveE += ss3.Energy[plane];
      aveE /= tjs.NumPlanes;
      myprt<<" Length "<<PosSep(ss3.Start, ss3.End)<<" aveE "<<(int)aveE<<" MeV";
      double shMaxAlong, along95;
      ShowerParams(aveE, shMaxAlong, along95);
      myprt<<" shMaxAlong "<<std::fixed<<std::setprecision(1)<<shMaxAlong<<" along95 "<<along95<<"\n";
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        CTP_t inCTP = plane;
        auto tp = MakeBareTP(tjs, ss3.Start, ss3.Dir, inCTP);
        myprt<<" Start "<<plane<<":"<<PrintPos(tjs, tp.Pos);
        tp = MakeBareTP(tjs, ss3.ChgPos, ss3.Dir, inCTP);
        myprt<<" ChgPos "<<PrintPos(tjs, tp.Pos);
        tp = MakeBareTP(tjs, ss3.End, ss3.Dir, inCTP);
        myprt<<" End "<<PrintPos(tjs, tp.Pos);
        myprt<<" Ang "<<std::setprecision(2)<<tp.Ang<<" projInPlane "<<tp.Delta;
        myprt<<" E "<<(int)ss3.Energy[plane];
        myprt<<"\n";
      } // plane
    } // prt
    return true;
  } // DefinePFPShower
  
  ////////////////////////////////////////////////
  bool UpdatePFPShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // updates the properties of ss3 when the list of pfparticles in this shower is changed
    if(ss3.PFPIDs.empty()) return false;
    
    std::string fcnLabel = inFcnLabel + ".UPS";

    std::cout<<"UpdatePFPShower needs revision\n";
    return false;
    
    // put all the tjs into a monster shower pfparticle so that
    // we can determine the start and end. This pfp is not stored
    PFPStruct shPFP = CreatePFP(tjs, ss3.TPCID);
    for(auto pid : ss3.PFPIDs) {
      auto& pfp = tjs.pfps[pid - 1];
      shPFP.TjIDs.insert(shPFP.TjIDs.end(), pfp.TjIDs.begin(), pfp.TjIDs.end());
    } // pid
    // Fit all the points to a line and fill the TP3s vector
    FindCompleteness(tjs, shPFP, true, true, false);
    if(shPFP.Tp3s.empty()) return false;
    // Sort Tp3s by distance from the start. Don't let the sort function
    // change the Tp3 direction vectors since we aren't using them.
    shPFP.Dir[0] = ss3.Dir;
//    SortByDistanceFromStart(tjs, shPFP, prt);
    ss3.Start = shPFP.Tp3s[0].Pos;
    ss3.End = shPFP.Tp3s[shPFP.Tp3s.size() - 1].Pos;
    ss3.Len = PosSep(ss3.Start, ss3.End);
    // Find the charge and accumulate the transverse position in three sections
    std::vector<double> chg(3);
    std::vector<double> sumt(3);
    std::vector<double> suml(3);
    // determine the separations that define the section boundaries
    double sec1 = ss3.Len / 3;
    double sec2 = 2 * ss3.Len / 3;
    Point2_t alongTran;
    for(auto& tp3 : shPFP.Tp3s) {
      if(tp3.dEdx <= 0) {
        std::cout<<"Ooops. No charge in Tp3\n";
        return false;
      }
      // Find the longitudinal and transverse distance from the shower start to the tp3
      FindAlongTrans(ss3.Start, ss3.Dir, tp3.Pos, alongTran);
      // determine which section this is in
      unsigned short section = 0;
      if(alongTran[0] > sec2) {
        section = 2;
      } else if(alongTran[0] > sec1) {
        section = 1;
      }
      chg[section] += tp3.dEdx;
      suml[section] += tp3.dEdx * alongTran[0];
      sumt[section] += tp3.dEdx * alongTran[1];
    } // tp3
    // something is wrong if there is no charge in the end sections
    if(chg[0] == 0 || chg[2] == 0) return false;
    for(unsigned short section = 0; section < 3; ++section) {
      if(chg[section] <= 0) continue;
      sumt[section] /= chg[section];
      suml[section] /= chg[section];
      std::cout<<"section "<<section<<" long "<<suml[section]<<" trans "<<sumt[section]<<"\n";
    } // section
    // reverse the direction?
    if(sumt[2] < sumt[0]) {
      std::swap(ss3.Start, ss3.End);
      // we don't need the Tp3s anymore so don't reverse them
      std::swap(sumt[0], sumt[2]);
      std::swap(suml[0], suml[2]);
      std::swap(chg[0], chg[2]);
    }
    // calculate the opening angle - Note that this is 1/2 the tangent of 
    // the opening angle...
    ss3.OpenAngle = std::abs((sumt[2] - sumt[0]) / (suml[2] - suml[0]));
    if(ss3.OpenAngle < 0.05 || ss3.OpenAngle > 0.2) {
      std::cout<<"UpdatePFPShower: Crazy openining angle "<<ss3.OpenAngle<<". Setting it to 0.15\n";
      ss3.OpenAngle = 0.15;
    }
    // check the direction vector using a largish component
    unsigned short useComp = 0;
    if(std::abs(ss3.Dir[2]) > std::abs(ss3.Dir[0])) useComp = 2;
    if(ss3.Dir[useComp] > 0 && ss3.Start[useComp] > ss3.End[useComp]) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) ss3.Dir[xyz] *= -1;
    } // reverse direction?
    // stash the charge in each section into PosErr. It will be transferred into
    // the shower Tj TP charge in MakePFPShowers. Put the "rms" into DirErr. Note that
    // the "rms" is really the average of the abs(transverse) positions. This is close
    // enough to the actual rms for this purpose.
    for(unsigned short section = 0; section < 3; ++section) {
      ss3.StartErr[section] = chg[section];
      ss3.DirErr[section] = sumt[section];
    }

    return true;
  } // UpdatePFPShower

  ////////////////////////////////////////////////
  bool DefineShower(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Defines the properties of a 2D shower using the trajectory points within the trajectories listed
    // in TjIDs. This wipes out any other information that may exist
    
    if(cotID > tjs.cots.size()) return false;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    if(ss.ShowerTjID == 0) return false;

    std::string fcnLabel = inFcnLabel + ".DS";
    
    if(ss.ParentID > 0) {
      std::cout<<fcnLabel<<" shouldn't be used for showers with a parent tj\n";
      return false;
    }
    
    FillPts(fcnLabel, tjs, cotID, prt);
    if(!FindChargeCenter(fcnLabel, tjs, cotID, prt)) {
      if(prt) {
        mf::LogVerbatim("TC")<<"Failed to find shower charge center 2S"<<cotID;
        Print2DShowers("DS", tjs, ss.CTP, false);
      }
      MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
      return false;
    }
    FindAngle(fcnLabel, tjs, cotID, prt);
    FillRotPos(fcnLabel, tjs, cotID, prt);
    if(!DefineShowerTj(fcnLabel, tjs, cotID, prt)) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failed to define Shower Tj. Killed the shower";
      MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
      return false;
    }
    FindNearbyTjs(fcnLabel, tjs, cotID, prt);
    DefineEnvelope(fcnLabel, tjs, cotID, prt);
    
    if(!tjs.MCPartList.empty()) {
      // Find the MC particle that matches with these InShower Tjs
      MCParticleListUtils tm{tjs};
      unsigned short nTruHits;
      unsigned int mcpIndex = tm.GetMCPartListIndex(ss, nTruHits);
      // Find the Tj that is closest to the start of this MC Particle
      if(mcpIndex != UINT_MAX) ss.TruParentID = tm.MCParticleStartTjID(mcpIndex, ss.CTP);
    }

    return true;

  } // DefineShower
  
  ////////////////////////////////////////////////
  bool AddTj(std::string inFcnLabel, TjStuff& tjs, int tjID, int cotID, bool doUpdate, bool prt)
  {
    // Adds the Tj to the shower and optionally updates the shower variables
    
    if(tjID <= 0 || tjID > (int)tjs.allTraj.size()) return false;
    if(cotID > tjs.cots.size()) return false;
    
    std::string fcnLabel = inFcnLabel + ".ATj";
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    Trajectory& tj = tjs.allTraj[tjID - 1];
    
    if(tj.CTP != ss.CTP) {
      mf::LogVerbatim("TC")<<fcnLabel<<" T"<<tjID<<" is in the wrong CTP "<<ss.ID;
      return false;
    }
    
    if(tj.AlgMod[kInShower]) {
      // see if it is in this shower
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjID) != ss.TjIDs.end()) {
        mf::LogVerbatim("TC")<<fcnLabel<<" T"<<tjID<<" is already in 2S"<<ss.ID;
        return true;
      } else {
        mf::LogVerbatim("TC")<<fcnLabel<<" T"<<tjID<<" is already in a different shower ";
        return false;
      }
    }

    ss.TjIDs.push_back(tjID);
    std::sort(ss.TjIDs.begin(), ss.TjIDs.end());
    // remove this ID from the NearTjIDs list
    for(unsigned short ii = 0; ii < ss.NearTjIDs.size(); ++ii) {
      if(ss.NearTjIDs[ii] == tjID) ss.NearTjIDs.erase(ss.NearTjIDs.begin() + ii);
    }
    // count the hits in the TPs
    unsigned short cnt = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if(tp.UseHit[ii]) ++cnt;
    } // ipt
    unsigned short newSize = ss.ShPts.size() + cnt;
    cnt = ss.ShPts.size();
    ss.ShPts.resize(newSize);
    // now add them
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) {
          unsigned int iht = tp.Hits[ii];
          ss.ShPts[cnt].HitIndex = iht;
          ss.ShPts[cnt].TID = tj.ID;
          ss.ShPts[cnt].Chg = tjs.fHits[iht].Integral;
          ss.ShPts[cnt].Pos[0] = tjs.fHits[iht].ArtPtr->WireID().Wire;
          ss.ShPts[cnt].Pos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
          ++cnt;
        }
      }
    } // ipt
    tj.AlgMod[kInShower] = true;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added T"<<tj.ID<<" to 2S"<<ss.ID;
    
    if(doUpdate) return DefineShower(fcnLabel, tjs, cotID, prt);
    return true;
    
  } // AddTj
  
  ////////////////////////////////////////////////
  bool RemoveTj(std::string inFcnLabel, TjStuff& tjs, int TjID, int cotID, bool doUpdate, bool prt)
  {
    // Removes the Tj from a shower
    
    if(TjID > (int)tjs.allTraj.size()) return false;
    if(cotID > tjs.cots.size()) return false;
    
    std::string fcnLabel = inFcnLabel + ".RTj";
    
    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[TjID - 1];
    if(!tj.AlgMod[kInShower]) return false;
    tj.AlgMod[kInShower] = false;
    tj.AlgMod[kShwrParent] = false;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    bool gotit = false;
    for(unsigned short ii = 0; ii < ss.TjIDs.size(); ++ii) {
      if(TjID == ss.TjIDs[ii]) {
        ss.TjIDs.erase(ss.TjIDs.begin() + ii);
        gotit = true;
        break;
      }
    } // ii
    if(!gotit) return false;
    // Removing a parent Tj?
    if(TjID == ss.ParentID) ss.ParentID = 0;
    // re-build everything?
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Remove T"<<TjID<<" from 2S"<<ss.ID;
    if(doUpdate) return DefineShower(fcnLabel, tjs, cotID, prt);
    return true;
  } // RemoveTj
  
  ////////////////////////////////////////////////
  void FindParent(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // look for a parent pfp for the shower. 
    
    if(ss3.ID == 0) return;
    if(ss3.CotIDs.size() < 2) return;
    
    std::string fcnLabel = inFcnLabel + ".FPar";
    
    double energy = 0;
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      if(ss.Energy <= 0) std::cout<<"3S"<<ss3.ID<<" energy not defined\n";
      energy += ss.Energy;
    } // ci
    energy /= (double)ss3.CotIDs.size();
    // the energy is probably under-estimated since there isn't a parent yet.
    energy *= 1.2;
    double shMaxAlong, along95;
    ShowerParams(energy, shMaxAlong, along95);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" Estimated energy "<<(int)energy<<" MeV shMaxAlong "<<shMaxAlong<<" along95 "<<along95;
    
    // look for the pfp that has a reasonable probability of being in the shower but with the
    // minimum along distance from the shower center
    float minAlong = -0.5 * shMaxAlong;
    // min cos angle btw the shower direction and the pfp direction (0.1 radians)
    float maxCosth = 0.995;
    unsigned short bestPFP = 0;
    std::cout<<"3S"<<ss3.ID<<" ->";
    for(auto ssid : ss3.CotIDs) std::cout<<" 2S"<<ssid;
    std::cout<<"\n";
    // temp vector to flag pfps that are parents - indexed by ID
    std::vector<bool> isParent(tjs.pfps.size() + 1, false);
    for(auto& oldSS3 : tjs.showers) {
      if(oldSS3.ID == 0) continue;
      isParent[oldSS3.ParentID] = true;
    } // pfp
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != ss3.TPCID) continue;
      // ignore neutrinos
      if(pfp.PDGCode == 14 || pfp.PDGCode == 14) continue;
      // ignore shower pfps
      if(pfp.PDGCode == 1111) continue;
      // ignore existing parents
      if(isParent[pfp.ID]) continue;
      // get a list of 2D showers that are associated with this pfp
      auto pfp_ss = GetAssns(tjs, "P", pfp.ID, "2S");
      if(!pfp_ss.empty()) {
        // There are 2D showers associated with this pfp. See which ones are
        // associated with this shower and which are different
        auto shared = SetIntersection(pfp_ss, ss3.CotIDs);
        // this pfp is associated with 2D showers that aren't associated with this 3D shower
        if(shared.empty()) continue;
        auto setDiff = SetDifference(pfp_ss, ss3.CotIDs);
        if(!setDiff.empty()) {
          std::cout<<" P"<<pfp.ID<<" ->";
          for(auto tjid : pfp.TjIDs) std::cout<<" T"<<tjid;
          std::cout<<" -> sslist:";
          for(auto ssid : pfp_ss) std::cout<<" 2S"<<ssid;
          std::cout<<" shared:";
          for(auto id : shared) std::cout<<" 2S"<<id;
          std::cout<<" setDiff:";
          for(auto id : setDiff) std::cout<<" 2S"<<id;
          std::cout<<"\n";
        } // !setDiff.empty()
      } // !sslist.empty()
      // find the end that is farthest away
      unsigned short pend = FarEnd(tjs, pfp, ss3.ChgPos);
      auto pToS = PointDirection(pfp.XYZ[pend], ss3.ChgPos);
      double costh = DotProd(pToS, ss3.Dir);
      if(costh < 0.5) continue;
      Point2_t alongTrans;
      // find the longitudinal and transverse components of the pfp start point relative to the
      // shower center
      FindAlongTrans(ss3.ChgPos, ss3.Dir, pfp.XYZ[pend], alongTrans);
      // ignore pfps that are not near the beginning of the shower
      if(energy > 1000 && pfp.ID < 5) {
        std::cout<<fcnLabel<<" chk high energy 3S"<<ss3.ID<<" E "<<(int)energy<<" shMaxAlong "<<shMaxAlong<<" along[0] "<<alongTrans[0]<<"\n";
      }
      if(alongTrans[0] > minAlong) continue;
      // offset the longitudinal distance by the expected shower max distance
      alongTrans[0] += shMaxAlong;
      float prob = InShowerProbLong(energy, alongTrans[0]);
      if(prob < 0.1) continue;
      // use the overall pfp direction instead of the starting direction. It may not be so
      // good if the shower develops quickly
      auto pfpDir = PointDirection(pfp.XYZ[pend], pfp.XYZ[1 - pend]);
      costh = DotProd(pfpDir, ss3.Dir);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<fcnLabel;
        myprt<<" 3S"<<ss3.ID;
        myprt<<" P"<<pfp.ID<<"_"<<pend;
        myprt<<" along "<<std::fixed<<std::setprecision(1)<<alongTrans[0]<<" trans "<<alongTrans[1];
        myprt<<std::setprecision(2)<<" prob "<<prob;
        myprt<<" costh "<<costh;
      } // prt
      if(bestPFP == 0) {
        bestPFP = pfp.ID;
        minAlong = alongTrans[0];
        maxCosth = costh;
      }
      if(alongTrans[0] > minAlong) continue;
      if(costh < maxCosth) continue;
      minAlong = alongTrans[0];
      maxCosth = costh;
      bestPFP = pfp.ID;
    } // pfp
    
    if(bestPFP == 0) return;
    mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" best P"<<bestPFP<<" minAlong "<<minAlong;
    ss3.ParentID = bestPFP;
    auto& pfp = tjs.pfps[bestPFP - 1];
    unsigned short pend = FarEnd(tjs, pfp, ss3.ChgPos);
    ss3.Vx3ID = pfp.Vx3ID[pend];
    
    // remove old parent tjs in the 2D showers
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      if(ss.ParentID > 0) UpdateShowerWithParent(fcnLabel, tjs, cid, 0, 10, prt);
    } // ci
    
    // add tj parents to 2D showers
    for(auto tjid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      for(auto cid : ss3.CotIDs) {
        auto& ss = tjs.cots[cid - 1];
        if(ss.CTP != tj.CTP) continue;
        ss.NeedsUpdate = true;
        // Add with parentFOM = 0;
        UpdateShowerWithParent(fcnLabel, tjs, cid, tjid, 0, prt);
      } // ci
    } // tjid
    ss3.NeedsUpdate = true;
    
  } // FindParent

  ////////////////////////////////////////////////
  void FindExternalParent(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Look for a parent trajectory that starts outside the shower and ends inside.

    /*
     # 0 Mode (<= 0 OFF, 1 = find showers before 3D match, 2 = find showers after 3D match)
     # 1 Max Tj MCSMom for a shower tag
     # 2 Max separation
     # 3 Min energy (MeV)
     # 4 rms width factor
     # 5 Min shower 1/2 width (WSE units)
     # 6 Min total Tj Pts
     # 7 Min Tjs
     # 8 max parent FOM
     # 9 max direction FOM
     # 10 max aspect ratio
     # 11 Debug in CTP (>10 debug cotID + 10)
     */

    if(cotID > tjs.cots.size()) return;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    // Ensure that it is valid
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".FEP";
    
    // See if anything needs to be done
    if(ss.SS3ID > 0) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" already matched in 3D";
      return;
    }
    
    // References to shower Tj points
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    int oldParent = ss.ParentID;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      if(oldParent > 0) {
        myprt<<fcnLabel<<" 2S"<<ss.ID<<" ST"<<stj.ID<<" Existing parent ID T"<<oldParent<<" parent FOM "<<ss.ParentFOM;
        myprt<<" attached to vertex "<<stj.VtxID[0];
        myprt<<" end0 at "<<PrintPos(tjs, stj.Pts[0].Pos);
      } else {
        myprt<<fcnLabel<<" ss.ID "<<ss.ID<<" ST"<<stj.ID<<" has no existing parent";
      }
    } // prt
    
//    if(ss.AspectRatio > tjs.ShowerTag[10] || ss.DirectionFOM > tjs.ShowerTag[9]) {
    if(ss.AspectRatio > tjs.ShowerTag[10]) {
      if(prt) mf::LogVerbatim("TC")<<" poor AspectRatio "<<ss.AspectRatio<<" or ss.DirectionFOM "<<ss.DirectionFOM;
      return;
    }
    
    float bestFOM = tjs.ShowerTag[8];
    int imTheBest = 0;
    unsigned short imTheBestEnd = 0;
    float bestTp1Sep = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
//      if(tj.AlgMod[kKilled] && !tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kKilled]) continue;
      // ignore shower Tjs. Note that this also rejects parent Tjs of other showers
      if(tj.AlgMod[kShowerTj]) continue;
      bool isInThisShower = (std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end());
      // ignore in-shower Tjs that aren't in this shower
      if(tj.AlgMod[kInShower] && !isInThisShower) continue;
      // ignore existing shower parents
      if(tj.AlgMod[kShwrParent]) continue;
      // Ignore short Tjs
      if(tj.Pts.size() < 5) continue;
      // See if this tj has an end that is near and end of the shower if it is inside the shower
      bool goodDirectionFOM = ss.DirectionFOM < tjs.ShowerTag[9];
      unsigned short useEnd = 0;
      // Check to see if the end of the Tj that will be used to determine if it is a parent is anywhere
      // near the end of the shower. It can't be a parent if there is another Tj that is further away...
      if(isInThisShower && tjs.UseAlg[kChkShwrParEnd]) {
        // find the end of the tj that is farthest away from the shower center
        useEnd = FarEnd(tjs, tj, stj.Pts[1].Pos);
        // determine which end of the shower is closest to that point
        auto& tp = tj.Pts[tj.EndPt[useEnd]];
        unsigned short shEnd = 0;
        if(PosSep2(stj.Pts[2].Pos, tp.Pos) < PosSep2(stj.Pts[0].Pos, tp.Pos)) shEnd = 1;
        // Ignore candidate parent Tjs whose endpoint is closer to the wrong end of the shower
        // if the shower direction is well known
        if(goodDirectionFOM && shEnd == 1) continue;
        // inspect the list of Tjs that are in the shower near that end.
        // Check the first (last) 20% of the points
        bool nearShowerEnd = false;
        short start = 0;
        short end = 0.2 * ss.ShPts.size();
        if(shEnd == 1) {
          // near end 1
          start = 0.8 * ss.ShPts.size();
          end = ss.ShPts.size();
        } // shEnd == 1
        unsigned short usid = tj.ID;
        for(unsigned short ipt = start; ipt < end; ++ipt) {
          if(ss.ShPts[ipt].TID == usid) {
            nearShowerEnd = true;
            break;
          }
        } // ipt
//        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tj.ID<<" is inside shower "<<ss.ID<<". Is it near an end? "<<nearShowerEnd;
        if(!nearShowerEnd) continue;
      } // isInThisShower
      // Check trajectories that were split by 3D vertex matching
      if(WrongSplitTj(fcnLabel, tjs, tj, useEnd, ss, prt)) continue;
      float tp1Sep, vx3Score;
      float fom = ParentFOM(fcnLabel, tjs, tj, useEnd, ss, tp1Sep, vx3Score, prt);
      // ignore the really bad ones
      if(fom > 2 * bestFOM) continue;
      // See if this tj is matched to a pfp that is the daughter of a neutrino
      bool attachedToNeutrino = false;
      if(tj.AlgMod[kMat3D]) {
        unsigned short indx = GetPFPIndex(tjs, tj.ID);
        if(indx < tjs.pfps.size()) {
          auto& pfp = tjs.pfps[indx];
          if(pfp.ParentID > 0) {
            auto& parPFP = tjs.pfps[pfp.ParentID - 1];
            if(parPFP.PDGCode == 14) {
              // set the fom low but keep looking in case this is not the real primary
              fom = 0.1;
              attachedToNeutrino = true;
              if(prt) mf::LogVerbatim("TC")<<"  T"<<tj.ID<<" is attached to a neutrino pfp";
            } // parent is a neutrino
          } // pfp has a parent
        } // valid indx
      } // matched in 3D
      // Allow a lower FOM if the separation is larger
      if(!attachedToNeutrino && bestTp1Sep > 0 && tp1Sep < bestTp1Sep) continue;
      bestFOM = fom;
      imTheBest = tj.ID;
      imTheBestEnd = useEnd;
      bestTp1Sep = tp1Sep;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" current best T"<<imTheBest<<" useEnd "<<useEnd<<" bestTp1Sep "<<bestTp1Sep;
    } // tj

    if(imTheBest < 1 || imTheBest > (int)tjs.allTraj.size()) return;
    if(bestFOM > tjs.ShowerTag[8]) imTheBest = 0;
    UpdateShowerWithParent(fcnLabel, tjs, cotID, imTheBest, bestFOM, prt);
    
  } // FindExternalParent
  
  ////////////////////////////////////////////////
  bool UpdateShowerWithParent(std::string inFcnLabel, TjStuff& tjs, int cotID, unsigned short newParent, float newParentFOM, bool prt)
  {
    // This updates all shower and shower Tj parameters when a new shower parent Tj is specified. This function also
    // removes an existing parent if newParent = 0. This function returns NeedsUpdate false if it is successful 
    
    if(cotID > tjs.cots.size()) return false;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    // Ensure that everything is valid
    if(ss.ID == 0) return false;
    if(newParent > tjs.allTraj.size()) return false;
    if(ss.TjIDs.empty()) return false;
    if(ss.ShowerTjID == 0) return false;
    
    // no existing parent and no new parent. Dumb but not an error
    if(ss.ParentID == 0 && newParent == 0) return true;
    
    std::string fcnLabel = inFcnLabel + ".USWP";
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" old Parent T"<<ss.ParentID<<" new T"<<newParent;
       
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    // remove an existing parent?
    if(ss.ParentID > 0) {
      auto& oldParent = tjs.allTraj[ss.ParentID - 1];
      oldParent.AlgMod[kShwrParent] = false;
      ss.ParentFOM = 10;
      stj.VtxID[0] = 0;
      // Not adding a new parent?
      if(newParent == 0) return true;
    } // remove an existing parent
    
    // make copies to allow error recovery
    auto oldSS = ss;
    auto oldStj = stj;
    
    ss.ParentID = newParent;
    ss.ParentFOM = newParentFOM;
    
    if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), ss.ParentID) == ss.TjIDs.end()) {
      // add it to the Shower Tj but don't update. This will be done below
      if(!AddTj(fcnLabel, tjs, ss.ParentID, cotID, false, prt)) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failed to add parent T"<<newParent<<" for some reason. Removing it from the shower";
        // remove it and update the shower
        RemoveTj(fcnLabel, tjs, ss.ParentID, cotID, true, prt);
        ss.ParentID = 0;
        ss.ParentFOM = 10;
        return false;
      }
    } // not in the shower yet

    ss.NeedsUpdate = true;
    
    Trajectory& ptj = tjs.allTraj[ss.ParentID - 1];
    unsigned short pend = FarEnd(tjs, ptj, stj.Pts[1].Pos);
    // make a copy of the start point so we can reverse the direction if needed
    auto ptp = ptj.Pts[ptj.EndPt[pend]];

    if(!tjs.MCPartList.empty()) {
      // get the truth if it exists
      MCParticleListUtils tm{tjs};
      unsigned short nTruHits;
      unsigned int mcpIndex = tm.GetMCPartListIndex(ss, nTruHits);
      // Find the Tj that is closest to the start of this MC Particle
      if(mcpIndex != UINT_MAX) ss.TruParentID = tm.MCParticleStartTjID(mcpIndex, ss.CTP);
    }

    // set the start vertex
    stj.VtxID[0] = ptj.VtxID[pend];
    // and dE/dx of the shower Tj
    stj.dEdx[0] = ptj.dEdx[pend];
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<< "  ParentID T" << ss.ParentID << " dEdx " << stj.dEdx[0]<<" attached to 2V"<<stj.VtxID[0];

    // reference to the point on the parent Tj that is furthest away from the shower
    // update the angle if the parent Tj is high-quality 
    if(ptj.Pts.size() > 20 && ptj.MCSMom > 100) ss.Angle = ptp.Ang;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"   shower angle "<<ss.Angle;
    
    // determine how much space we need for the points (hits) ala FindChargeCenter
    unsigned int cnt = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if(tp.UseHit[ii]) ++cnt;
      } // ipt
    } // it
    ss.ShPts.resize(cnt);

    // Populate the Pts vector
    cnt = 0;
    float totChg = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      Trajectory& tj = tjs.allTraj[ss.TjIDs[it] - 1];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        // create a point for every hit
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tp.UseHit[ii]) {
            unsigned int iht = tp.Hits[ii];
            if(cnt > ss.ShPts.size() - 1) {
              std::cout<<fcnLabel<<" bad cnt "<<cnt<<" "<<ss.ShPts.size()<<"\n";
              return false;
            }
            ss.ShPts[cnt].HitIndex = iht;
            ss.ShPts[cnt].TID = tj.ID;
            ss.ShPts[cnt].Chg = tjs.fHits[iht].Integral;
            ss.ShPts[cnt].Pos[0] = tjs.fHits[iht].ArtPtr->WireID().Wire;
            ss.ShPts[cnt].Pos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
            totChg += ss.ShPts[cnt].Chg;
            ++cnt;
          }
        } // ii
      } // ipt
    } // it
    // Put the total charge into the shower Tj
    stj.AveChg = totChg;
    // Update the shower energy
    ss.Energy = ShowerEnergy(tjs, ss);
    
    // ensure that the shower direction is consistent with the parent start TP
    TrajPoint tp;
    MakeBareTrajPoint(tjs, ptp, stj.Pts[1], tp);
    float costh = DotProd(tp.Dir, stj.Pts[1].Dir);
    if(costh < 0) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) ptp.Dir[xyz] *= -1;
      ptp.Ang = atan2(ptp.Dir[1], ptp.Dir[0]);
      ss.Angle = ptp.Ang;
    }

    // define the angle of the shower Tj points
    for(auto& stp : stj.Pts) {
      stp.Ang = ptp.Ang;
      stp.Dir = ptp.Dir;
    } // stp

    // fill the rotated points
    FillRotPos(fcnLabel, tjs, cotID, prt);

    // AnalyzeRotPos calculates the charge (Chg), the charge center (HitPos) and
    // the transverse shower rms at each Tp in the shower Tj
    if(!AnalyzeRotPos(fcnLabel, tjs, cotID, prt)) {
      mf::LogVerbatim("TC")<<fcnLabel<< " Failure from AnalyzeRotPos. Killing this shower";
      MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
      return false;
    }
    
    // See if the shower with the new parent is lousy
    if(ss.AspectRatio > 0.8) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Bad AspectRatio "<<ss.AspectRatio<<". Removing the parent";
      ss = oldSS;
      stj = oldStj;
      ss.NeedsUpdate = false;
      return true;
    }
    
    // now define the positions of the shower Tj points.
    float minAlong = ss.ShPts[0].RotPos[0];
    float maxAlong = ss.ShPts[ss.ShPts.size()-1].RotPos[0];
    
    // put the shower Tj start position at the parent Tj start. Note that this is not
    // necessarily the same as HitPos[0]
    TrajPoint& stp0 = stj.Pts[0];
    stp0.Pos = ptj.Pts[pend].Pos;
    
    // rotate the other points
    TrajPoint& stp1 = stj.Pts[1];
    stp1.Pos[0] = stp0.Pos[0] + stp1.Dir[0] * (0 - minAlong);
    stp1.Pos[1] = stp0.Pos[1] + stp1.Dir[1] * (0 - minAlong);
    
    TrajPoint& stp2 = stj.Pts[2];
    stp2.Pos[0] = stp0.Pos[0] + stp1.Dir[0] * (maxAlong - minAlong);
    stp2.Pos[1] = stp0.Pos[1] + stp1.Dir[1] * (maxAlong - minAlong);

    FindNearbyTjs(fcnLabel, tjs, cotID, prt);
    // modify the shower tj points so that DefineEnvelope gives a reasonable result
    // TODO: Do this correctly, perhaps scaling by the aspect ratio
    if(stp0.DeltaRMS < 1) stp0.DeltaRMS = 1;
    float expectedRMS = 0.07 * PosSep(stp0.Pos, stp2.Pos);
//    std::cout<<"RMS "<<ss.ID<<" stp0 rms "<<stp0.DeltaRMS<<" Expected stp2 rms"<<expectedRMS<<" "<<stp2.DeltaRMS<<"\n";
    if(stp2.DeltaRMS < expectedRMS) stp2.DeltaRMS = expectedRMS;
    for(unsigned short nit = 0; nit < 2; ++nit) {
      DefineEnvelope(fcnLabel, tjs, cotID, prt);
      if(AddTjsInsideEnvelope(fcnLabel, tjs, cotID, prt)) ss.NeedsUpdate = true;
      if(!ss.NeedsUpdate) break;
    } // nit

    stj.AlgMod[kInShower] = true;
    stj.AlgMod[kShwrParent] = true;

    return true;
  } // UpdateShowerWithParent

  ////////////////////////////////////////////////
  unsigned short FarEnd(TjStuff& tjs, const Trajectory& tj, Point2_t& pos)
  {
    // Returns the end (0 or 1) of the Tj that is furthest away from the position pos
    if(tj.ID == 0) return 0;
    if(PosSep2(tj.Pts[tj.EndPt[1]].Pos, pos) > PosSep2(tj.Pts[tj.EndPt[0]].Pos, pos)) return 1;
    return 0;
  } // FarEnd

  ////////////////////////////////////////////////
  unsigned short FarEnd(TjStuff& tjs, const PFPStruct& pfp, Point3_t& pos)
  {
    // Returns the end (0 or 1) of the pfp that is furthest away from the position pos
    if(pfp.ID == 0) return 0;
    if(PosSep2(pfp.XYZ[1], pos) > PosSep2(pfp.XYZ[0], pos)) return 1;
    return 0;
  } // FarEnd
  
  ////////////////////////////////////////////////
  PFPStruct CreateFakePFP(const TjStuff& tjs, const ShowerStruct3D& ss3)
  {
    // Creates a fake PFParticle on the shower axis with a TP3S point every cm. This is
    // used to find the separation between the shower core and a PFParticle or trajectory
    auto pfp = CreatePFP(tjs, ss3.TPCID);
    pfp.XYZ[0] = ss3.Start;
    pfp.Dir[0] = ss3.Dir;
    pfp.XYZ[1] = ss3.End;
    pfp.Dir[1] = ss3.Dir;
    pfp.Dir[0] = PointDirection(ss3.Start, ss3.End);
    unsigned short npts = ss3.Len;
    if(npts < 2) return pfp;
    pfp.Tp3s.resize(npts);
    pfp.Tp3s[0].Pos = ss3.Start;
    for(unsigned short ipt = 1; ipt < npts; ++ipt) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pfp.Tp3s[ipt].Pos[xyz] = pfp.Tp3s[ipt-1].Pos[xyz] + pfp.Dir[0][xyz];
    } // ipt
    return pfp;
  } // CreateFakePFP

  ////////////////////////////////////////////////
  bool IsInShower(const TjStuff& tjs, const std::vector<int> TjIDs)
  {
    // Vote for the list of Tjs (assumed associated with a PFParticle) being in a shower
    if(TjIDs.empty()) return false;
    unsigned short cnt = 0;
    for(auto tjid : TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kInShower]) ++cnt;
    } // tjid
    return (cnt > 1);
  } // IsInShower

  ////////////////////////////////////////////////
  void ShowerParams(double showerEnergy, double& shMaxAlong, double& along95)
  {
    // Returns summary properties of photon showers parameterized in the energy range 50 MeV < E_gamma < 1 GeV:
    // shMaxAlong = the longitudinal distance (cm) between the start of the shower and the center of charge
    // along95 = the longitudinal distance (cm) between the start of the shower and 95% energy containment
    // all units are in cm
    if(showerEnergy < 50) {
      shMaxAlong = 0;
      along95 = 0;
      return;
    }
    shMaxAlong = 7.0 * log(showerEnergy / 15);
    // This needs to be increased. TODO: Study higher energy showers
    shMaxAlong *= 1.7;
    // The 95% containment is reduced a bit at higher energy
    double scale = 2.75 - 9.29E-4 * showerEnergy;
    if(scale < 2) scale = 2;
    along95 = scale * shMaxAlong;
  } // ShowerParams
  
  ////////////////////////////////////////////////
  double ShowerParamTransRMS(double showerEnergy, double along)
  {
    // returns the pareameterized width rms of a shower at along relative to the shower start
    double shMaxAlong, shE95Along;
    ShowerParams(showerEnergy, shMaxAlong, shE95Along);
    if(shMaxAlong <= 0) return 0;
    double tau = along / shMaxAlong;
    // The shower width is modeled as a simple cone that scales with tau
    double rms = -0.4 + 2.5 * tau;
    if(rms < 0.5) rms = 0.5;
    return rms;
  } // ShowerParamTransRMS

  ////////////////////////////////////////////////
  double InShowerProbLong(double showerEnergy, double along)
  {
    // Returns the likelihood that the point at position along (cm) is inside an EM shower
    // having showerEnergy (MeV). The variable along is relative to shower start.
    
    if(showerEnergy < 50) return 0;
    
    double shMaxAlong, shE95Along;
    ShowerParams(showerEnergy, shMaxAlong, shE95Along);
    // 50% of the shower energy is deposited between 0 < shMaxAlong < 1, which should be obvious considering
    // that is the definition of the shower max, so the probability should be ~1 at shMaxAlong = 1.
    // The Geant study shows that 95% of the energy is contained within 2.5 * shMax and has a small dependence 
    // on the shower energy, which is modeled in ShowerParams. This function uses a
    // sigmoid likelihood function is constructed with these constraints using the scaling variable tau 
    double tau = along / shMaxAlong;
    if(tau < -1 || tau > 4) return 0;
    
    double tauHalf, width;
    if(tau > 0) {
      tauHalf = 1.7;
      width = 0.35;
    } else {
      // Allow for some uncertainty in the shower start position
      tau = -tau;
      tauHalf = 0.2;
      width = 0.1;
    }

    double prob = 1 / (1 + exp((tau - tauHalf)/width));
    return prob;
    
  } // InShowrProbLong

  ////////////////////////////////////////////////
  double InShowerProbTrans(double showerEnergy, double along, double trans)
  {
    // Returns the likelihood that the point, (along, trans) (cm), is inside an EM shower having energy showerEnergy (MeV)
    // where along is relative to the shower start position and trans is the radial distance. 

    if(showerEnergy < 50) return 0;
    double rms = ShowerParamTransRMS(showerEnergy, along);
    trans = std::abs(trans);
    double prob = exp(-0.5 * trans / rms);
    return prob;
    
  } // InShowerProbTrans

  ////////////////////////////////////////////////
  double InShowerProbParam(double showerEnergy, double along, double trans)
  {
    return InShowerProbLong(showerEnergy, along) * InShowerProbTrans(showerEnergy, along, trans);
  } // InShowerProbParam

  ////////////////////////////////////////////////
  float InShowerProb(std::string inFcnLabel, const TjStuff& tjs, const ShowerStruct3D& ss3, const PFPStruct& pfp)
  {
    // returns a likelihood (0 - 1) that the pfp particle belongs in shower ss3
    
    std::string fcnLabel = inFcnLabel + ".ISP";
    
    // Find the average energy
    double energy = 0;
    double cnt = 0;
    for(auto nrg : ss3.Energy) {
      if(nrg <= 0) continue;
      energy += nrg;
      ++cnt;
    } // nrg
    if(cnt == 0) return 0;
    energy /= cnt;
    mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" energy "<<energy;
    Point2_t alongTrans;
    for(unsigned short end = 0; end < 2; ++end) {
      FindAlongTrans(ss3.Start, ss3.Dir, pfp.XYZ[end], alongTrans);
      mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" end "<<end<<" prob "<<InShowerProbParam(energy, alongTrans[0], alongTrans[1]);
    } // end
    return 0;
/*
    // pull on the angle difference
    double dangPull = 0;
    double pfpLen = PosSep(pfp.XYZ[0], pfp.XYZ[1]);
    // calculate the angle pull if it is long enough
    if(pfpLen > 5) {
      double costh = std::abs(DotProd(ss3.Dir, pfp.Dir[0]));
      // assume the shower opening angle should be 0.15 for 90% containment which is about
      // 2+ sigma
      double dangErr = 0.07;
      dangPull = acos(costh) / dangErr;
      if(dangPull) return 0;
    } // pfpLen > 5 cm
    // find the closest approach
    // Create a temp pfp with Tp3s populated along the shower axis
    auto shPFP = CreateFakePFP(tjs, ss3);
    double minSep2 = 1E6;
    unsigned short minSepPt = 0;
    double maxSep2 = 0;
    unsigned short maxSepPt = 0;
    for(unsigned short end = 0; end < 2; ++end) {
      for(unsigned short ipt = 0; ipt < shPFP.Tp3s.size(); ++ipt) {
        auto& tp3 = shPFP.Tp3s[ipt];
        // find separation^2 btw each point on the fake PFP and the end points of the pfp that
        // we are evaluating
        double sep2 = PosSep2(tp3.Pos, pfp.XYZ[end]); 
        if(sep2 < minSep2) {
          minSep2 = sep2;
          minSepPt = ipt;
        } // minSep2 test
        if(sep2 > maxSep2) {
          maxSep2 = sep2;
          maxSepPt = ipt;
        } // maxSep2 test
      } // tp3
    } // end
    std::cout<<"minSepPt "<<minSepPt<<" maxSepPt "<<maxSepPt<<"\n";
    return 0;
*/
  } // InShowerProb
  
  
  ////////////////////////////////////////////////
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, float& tp1Sep, float& vx2Score, bool prt)
  {
    // returns a FOM for the trajectory at the end point being the parent of ss and the end which
    // was matched.
    
    vx2Score = 0;
    tp1Sep = 0;
    
    if(tjEnd > 1) return 1000;
    if(ss.Energy == 0) return 1000;

    if(ss.ID == 0) return 1000;
    if(ss.TjIDs.empty()) return 1000;
    if(ss.ShowerTjID == 0) return 1000;
    
    std::string fcnLabel = inFcnLabel + ".PFOM";

    float fom = 0;
    float cnt = 0;
    // Shower charge center TP
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
    // get the end that is farthest away from the shower center
    tjEnd = FarEnd(tjs, tj, stp1.Pos);
    // prospective parent TP
    TrajPoint& ptp = tj.Pts[tj.EndPt[tjEnd]];
    // find the along and trans components in WSE units relative to the
    // shower center
    Point2_t alongTrans;
    FindAlongTrans(stp1.Pos, stp1.Dir, ptp.Pos, alongTrans);
    // TODO: we could return here if the shower direction is well defined and
    // alongTrans[0] is > 0
    tp1Sep = std::abs(alongTrans[0]);
    // Find the expected shower max along distance (cm)
    double shMaxAlong, shE95Along;
    ShowerParams(ss.Energy, shMaxAlong, shE95Along);
    double along = tjs.WirePitch * tp1Sep;
    float prob = InShowerProbLong(ss.Energy, along);
    if(prob < 0.05) return 100;
    // The transverse position must certainly be less than the longitudinal distance
    // to shower max
    if(alongTrans[1] > shMaxAlong) return 100;
    float sepFOM = std::abs(along - shMaxAlong) / 5;
    ++cnt;
    // make a tp between the supposed parent TP and the shower center
    TrajPoint tp;
    if(!MakeBareTrajPoint(tjs, ptp, stp1, tp)) return 100;
    // we have three angles to compare. The ptp angle, the shower angle and
    // the tp angle. 
    float dang1 = DeltaAngle(ptp.Ang, stp1.Ang);
    float dang1FOM = dang1 / 0.1;
    fom += dang1FOM;
    ++cnt;
    float dang2 = DeltaAngle(ptp.Ang, tp.Ang);
    float dang2FOM = dang1 / 0.1;
    fom += dang2FOM;
    ++cnt;
    // the environment near the parent start should be clean.
    std::vector<int> tjlist(1);
    tjlist[0] = tj.ID;
    // check for a vertex at this end and include the vertex tjs if the vertex is close
    // to the expected shower max position
    float vx2Sep = 0;
    // put in a largish FOM value for Tjs that don't have a vertex
    float vxFOM = 10;
    if(tj.VtxID[tjEnd] > 0) {
      VtxStore& vx2 = tjs.vtx[tj.VtxID[tjEnd] - 1];
      vx2Sep = PosSep(vx2.Pos, stp1.Pos);
      vx2Score = vx2.Score;
      tjlist = GetVtxTjIDs(tjs, vx2);
      vxFOM = std::abs(shMaxAlong - vx2Sep) / 20;
    } // 2D vertex exists
    fom += vxFOM;
    ++cnt;
    float chgFrac = ChgFracNearPos(tjs, ptp.Pos, tjlist);
    float chgFracFOM = (1 - chgFrac) / 0.1;
    fom += chgFracFOM;
    ++cnt;
    // Fraction of wires that have a signal between the parent start and the shower center
    float chgFracBtw = ChgFracBetween(tjs, ptp, stp1.Pos[0], false);
    float chgFrcBtwFOM = (1 - chgFrac) / 0.1;
    fom += chgFrcBtwFOM;
    ++cnt;

    // take the average
    fom /= cnt;
    // divide by the InShowerProbability
    fom /= prob;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel;
      myprt<<" 2S"<<ss.ID;
      myprt<<" T"<<tj.ID<<"_"<<tjEnd<<" Pos "<<PrintPos(tjs, ptp);
      myprt<<std::fixed<<std::setprecision(2);
      myprt<<" tp1Sep "<<std::fixed<<std::setprecision(1)<<tp1Sep<<" sepPull "<<sepFOM;
      myprt<<" InShowrProb "<<prob;
      myprt<<" dang1 "<<dang1<<" FOM "<<dang1FOM;
      myprt<<" dang2 "<<dang2<<" FOM "<<dang2FOM;
      myprt<<" vx2Score "<<vx2Score<<" vxFOM "<<vxFOM;
      myprt<<" chgFrac "<<chgFrac<<" chgFracFOM "<<chgFracFOM;
      myprt<<" chgFracBtw "<<chgFracBtw<<" chgFrcBtwFOM "<<chgFrcBtwFOM;
      myprt<<" FOM "<<fom;
    }
    return fom;
    
  } // ParentFOM

  ////////////////////////////////////////////////
  bool WrongSplitTj(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt)
  {
    // Returns true if the trajectory was split by a 3D vertex match and the end of this trajectory is further
    // away from the shower than the partner trajectory
    // Here is a cartoon showing what we are trying to prevent. The shower is represented by a box. The trajectory
    // that is shown as (---*---) was originally reconstructed as a single trajectory. It was later split at the * position
    // by matching in 3D into two trajectories with ID = 1 and 2. We don't want to consider Tj 1 using end 0 as a parent for
    // the shower. Tj is more likely to be the real parent
    //
    //  1111111111 2222222  TjID
    //  0        1 0     1  Tj end
    //               --------------
    //               |            |
    //  ----------*-------        |
    //               |            |
    //               --------------
    if(!tj.AlgMod[kComp3DVx]) return false;
    if(tjEnd > 1) return false;
    
    std::string fcnLabel = inFcnLabel + ".WSTj";
    
    // See if the other end is the end that was split. It should have a vertex with Topo = 8 or 11
    unsigned short otherEnd = 1 - tjEnd;
//    if(prt) mf::LogVerbatim("TC")<<"WSTj: otherEnd "<<otherEnd<<" vtxID "<<tj.VtxID[otherEnd];
    if(tj.VtxID[otherEnd] == 0) return false;
    unsigned short ivx = tj.VtxID[otherEnd] - 1;
    // A vertex exists but not a 3D split vertex
    if(tjs.vtx[ivx].Topo != 8 && tjs.vtx[ivx].Topo != 10) return false;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Primary candidate "<<tj.ID<<" was split by a 3D vertex";
    return true;
    
  } // WrongSplitTj

  ////////////////////////////////////////////////
  void MergeNearby2DShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    if(!tjs.UseAlg[kMergeNrShowers]) return;
    if(tjs.cots.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".MNS";
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" list";
      for(auto& ss : tjs.cots) {
        if(ss.CTP != inCTP) continue;
        if(ss.ID == 0) continue;
        myprt<<"  ss.ID "<<ss.ID<<" NearTjs";
        for(auto& id : ss.NearTjIDs) myprt<<" "<<id;
        myprt<<"\n";
      }
    } // prt
    
    bool keepMerging = true;
    while(keepMerging) {
      keepMerging = false;
      for(unsigned short ci1 = 0; ci1 < tjs.cots.size() - 1; ++ci1) {
        ShowerStruct& ss1 = tjs.cots[ci1];
        if(ss1.CTP != inCTP) continue;
        if(ss1.ID == 0) continue;
        if(ss1.TjIDs.empty()) continue;
        // put the inshower tjs and the nearby tjs into one list
        std::vector<int> ss1list = ss1.TjIDs;
        ss1list.insert(ss1list.end(), ss1.NearTjIDs.begin(), ss1.NearTjIDs.end());
        for(unsigned short ci2 = ci1 + 1; ci2 < tjs.cots.size(); ++ci2) {
          ShowerStruct& ss2 = tjs.cots[ci2];
          if(ss2.CTP != inCTP) continue;
          if(ss2.ID == 0) continue;
          if(ss2.TjIDs.empty()) continue;
          std::vector<int> ss2list = ss2.TjIDs;
          ss2list.insert(ss2list.end(), ss2.NearTjIDs.begin(), ss2.NearTjIDs.end());
          std::vector<int> shared = SetIntersection(ss1list, ss2list);
          if(shared.empty()) continue;
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<fcnLabel<<" Merge 2S"<<ss2.ID<<" into 2S"<<ss1.ID<<"? shared nearby:";
            for(auto tjid : shared) myprt<<" T"<<tjid;
          } // prt
          // add the shared Tjs to ss1 if they meet the requirements
          bool doMerge = false;
          for(auto& tjID : shared) {
            bool inSS1 = (std::find(ss1.TjIDs.begin(), ss1.TjIDs.end(), tjID) != ss1.TjIDs.end());
            bool inSS2 = (std::find(ss2.TjIDs.begin(), ss2.TjIDs.end(), tjID) != ss2.TjIDs.end());
            if(inSS1 && !inSS2) doMerge = true;
            if(!inSS1 && inSS2) doMerge = true;
            // need to add it?
            if(inSS1 || inSS2) continue;
            auto& tj = tjs.allTraj[tjID - 1];
            // ignore long muons
            if(tj.PDGCode == 13 && tj.Pts.size() > 100 && tj.ChgRMS < 0.5) {
              if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" T"<<tj.ID<<" looks like a muon. Don't add it";
              continue;
            }
            if(AddTj(fcnLabel, tjs, tjID, ss1.ID, false, prt)) doMerge = true;
          } // tjID
          if(!doMerge) continue;
          if(MergeShowersAndStore(fcnLabel, tjs, ss1.ID, ss2.ID, prt)) {
            Trajectory& stj = tjs.allTraj[ss1.ShowerTjID - 1];
            stj.AlgMod[kMergeNrShowers] = true;
            if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" success";
            keepMerging = true;
            break;
          }
        } // ci2
      } // ci1
    } // keepMerging
    return;
  } //MergeNearby2DShowers

  ////////////////////////////////////////////////
  void MergeOverlap(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Merge showers whose envelopes overlap each other
    
    /*
     # 0 Mode (<= 0 OFF, 1 = find showers before 3D match, 2 = find showers after 3D match)
     # 1 Max Tj MCSMom for a shower tag
     # 2 Max separation
     # 3 Min energy (MeV)
     # 4 rms width factor
     # 5 Min shower 1/2 width (WSE units)
     # 6 Min total Tj Pts
     # 7 Min Tjs
     # 8 max parent FOM
     # 9 max direction FOM
     # 10 max aspect ratio
     # 11 Debug in CTP (>10 debug cotID + 10)
     */
    
    if(tjs.ShowerTag[2] <= 0) return;
    if(!tjs.UseAlg[kMergeOverlap]) return;
    if(tjs.cots.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".MO";
    
    // Require that the maximum separation is about two radiation lengths
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" checking using separation cut "<<tjs.ShowerTag[2];
    
    float sepCut2 = tjs.ShowerTag[2] * tjs.ShowerTag[2];
    
    // Iterate if a merge is done
    bool didMerge = true;
    while(didMerge) {
      didMerge = false;
      // See if the envelopes overlap
      for(unsigned short ict = 0; ict < tjs.cots.size() - 1; ++ict) {
        auto& iss = tjs.cots[ict];
        if(iss.ID == 0) continue;
        if(iss.TjIDs.empty()) continue;
        if(iss.CTP != inCTP) continue;
        for(unsigned short jct = ict + 1; jct < tjs.cots.size(); ++jct) {
          auto& jss = tjs.cots[jct];
          if(jss.ID == 0) continue;
          if(jss.TjIDs.empty()) continue;
          if(jss.CTP != iss.CTP) continue;
          bool doMerge = false;
          for(auto& ivx : iss.Envelope) {
            doMerge = PointInsideEnvelope(ivx, jss.Envelope);
            if(doMerge) break;
          } // ivx
          if(!doMerge) {
            for(auto& jvx : jss.Envelope) {
              doMerge = PointInsideEnvelope(jvx, iss.Envelope);
              if(doMerge) break;
            } // ivx
          }
          if(!doMerge) {
            // check proximity between the envelopes
            for(auto& ivx : iss.Envelope) {
              for(auto& jvx : jss.Envelope) {
                if(PosSep2(ivx, jvx) < sepCut2) {
                  if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Envelopes "<<ict<<" "<<jct<<" are close "<<PosSep(ivx, jvx)<<" cut "<<tjs.ShowerTag[2];
                  doMerge = true;
                  break;
                }
              } // jvx
              if(doMerge) break;
            } // ivx
          } // !domerge
          if(!doMerge) continue;
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Merge "<<iss.ID<<" and "<<jss.ID;
          if(MergeShowersAndStore(fcnLabel, tjs, iss.ID, jss.ID, prt)) {
            Trajectory& stj = tjs.allTraj[iss.ShowerTjID - 1];
            stj.AlgMod[kMergeOverlap] = true;
            didMerge = true;
            break;
          } else {
            if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Merge failed";
          }
        } // jct
      } // ict
    } // didMerge

  } // MergeOverlap
  
  ////////////////////////////////////////////////
  void MergeShowerChain(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Merge chains of 3 or more showers that lie on a line
    
    if(!tjs.UseAlg[kMergeShChain]) return;
    
    std::string fcnLabel = inFcnLabel + ".MSC";
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<": MergeShowerChain inCTP "<<inCTP;
    
    std::vector<int> sids;
    std::vector<TrajPoint> tpList;
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      ShowerStruct& iss = tjs.cots[ict];
      if(iss.ID == 0) continue;
      if(iss.TjIDs.empty()) continue;
      if(iss.CTP != inCTP) continue;
      // save the shower ID
      sids.push_back(iss.ID);
      // and the shower center TP
      tpList.push_back(tjs.allTraj[iss.ShowerTjID - 1].Pts[1]);
    } // ict
    if(sids.size() < 3) return;
    
    // sort by wire so the chain order is reasonable
    std::vector<SortEntry> sortVec(sids.size());
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].length = tpList[ii].Pos[0];
    }
    std::sort(sortVec.begin(), sortVec.end(), lessThan);
    auto tsids = sids;
    auto ttpList = tpList;
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      sids[ii] = tsids[indx];
      tpList[ii] = ttpList[indx];
    }
    
    // TODO: These cuts should be generalized somehow
    float minSep = 150;
    float maxDelta = 30;
    for(unsigned short ii = 0; ii < sids.size() - 2; ++ii) {
      auto& iss = tjs.cots[sids[ii] - 1];
      if(iss.ID == 0) continue;
      unsigned short jj = ii + 1;
      auto& jss = tjs.cots[sids[jj] - 1];
      if(jss.ID == 0) continue;
      std::vector<int> chain;
      float sepij = PosSep(tpList[ii].Pos, tpList[jj].Pos);
      if(sepij > minSep) continue;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" i2S"<<iss.ID<<" "<<PrintPos(tjs, tpList[ii].Pos)<<" j2S"<<jss.ID<<" "<<PrintPos(tjs, tpList[jj].Pos)<<" sepij "<<sepij;
      // draw a line between these points
      TrajPoint tp;
      MakeBareTrajPoint(tjs, tpList[ii], tpList[jj], tp);
//      PrintTrajPoint("ij", tjs, 0, 1, 0, tp);
      for(unsigned short kk = jj + 1; kk < sids.size(); ++kk) {
        auto& kss = tjs.cots[sids[kk] - 1];
        if(kss.ID == 0) continue;
        float sepjk = PosSep(tpList[jj].Pos, tpList[kk].Pos);
        float delta = PointTrajDOCA(tjs, tpList[kk].Pos[0], tpList[kk].Pos[1], tp);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<"   k2S"<<kss.ID<<" "<<PrintPos(tjs, tpList[kk].Pos)<<" sepjk "<<sepjk<<" delta "<<delta;
          if(sepjk > minSep || delta > maxDelta) {
            myprt<<" failed separation "<<minSep<<" or delta cut "<<maxDelta;
          } else {
            myprt<<" add to the chain";
          }
        }
        if(sepjk > minSep || delta > maxDelta) {
          // clear a short chain?
          if(chain.size() > 2) {
            // merge this chain
            int newID = MergeShowers(fcnLabel, tjs, chain, prt);
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<fcnLabel<<" merged chain";
              for(auto ssID : chain) myprt<<" "<<ssID;
              myprt<<" -> new ssID "<<newID;
            } // prt
          } // long chain
          chain.clear();
          break;
        }  else {
          // add this shower to the chain
          if(chain.empty()) {
            chain.resize(3);
            chain[0] = sids[ii]; chain[1] = sids[jj]; chain[2] = sids[kk]; 
          } else {
            chain.push_back(sids[kk]);
          }
          // Refine the TP position and direction
          MakeBareTrajPoint(tjs, tpList[ii], tpList[kk], tp);
//          PrintTrajPoint("ik", tjs, 0, 0, chain.size(), tp);
        } // add to an existing chain
      } // kk
      // push the last one
      if(chain.size() > 2) {
        int newID = MergeShowers(fcnLabel, tjs, chain, prt);
        if(newID > 0 && AddTjsInsideEnvelope(fcnLabel, tjs, newID, prt)) DefineShower(fcnLabel, tjs, newID - 1, prt);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<" merged chain";
          for(auto ssID : chain) myprt<<" "<<ssID;
          myprt<<" -> new ssID "<<newID;
        } // prt
      } // long chain
    } // ii
    
  } // MergeShowerChain
  
  ////////////////////////////////////////////////
  void MergeSubShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Merge small showers that are downstream of larger showers
    
    if(!tjs.UseAlg[kMergeSubShowers]) return;
    
    std::string fcnLabel = inFcnLabel + ".MSS";
    

    bool newCuts = (tjs.ShowerTag[0] == 4);
    
    if(prt) {
      if(newCuts) {
        mf::LogVerbatim("TC")<<fcnLabel<<" MergeSubShowers checking using ShowerParams";
      } else {
        mf::LogVerbatim("TC")<<fcnLabel<<" MergeSubShowers checking using radiation length cut ";
      }
    } // prt

    constexpr float radLen = 14 / 0.3;
    
    bool keepMerging = true;
    while(keepMerging) {
      keepMerging = false;
      // sort by decreasing energy
      std::vector<SortEntry> sortVec;
      for(auto& ss : tjs.cots) {
        if(ss.ID == 0) continue;
        if(ss.CTP != inCTP) continue;
        SortEntry se;
        se.index = ss.ID - 1;
        se.length = ss.Energy;
        sortVec.push_back(se);
      } // ss
      if(sortVec.size() < 2) return;
      std::sort(sortVec.begin(), sortVec.end(), greaterThan);
      for(unsigned short ii = 0; ii < sortVec.size() - 1; ++ii) {
        if(inCTP == 0) std::cout<<"MSS ii "<<ii<<" "<<sortVec[ii].index<<" E "<<(int)sortVec[ii].length<<"\n";
        ShowerStruct& iss = tjs.cots[sortVec[ii].index];
        if(iss.ID == 0) continue;
        // this shouldn't be done to showers that are ~round
        if(iss.AspectRatio > 0.5) continue;
        TrajPoint& istp1 = tjs.allTraj[iss.ShowerTjID - 1].Pts[1];
        double shMaxAlong, along95;
        ShowerParams((double)iss.Energy, shMaxAlong, along95);
        // convert along95 to a separation btw shower max and along95
        along95 -= shMaxAlong;
        // convert to WSE
        along95 /= tjs.WirePitch;
        for(unsigned short jj = ii + 1; jj < sortVec.size(); ++jj) {
          if(inCTP == 0) std::cout<<"   jj "<<jj<<" "<<sortVec[jj].index<<" E "<<(int)sortVec[jj].length<<"\n";
          ShowerStruct& jss = tjs.cots[sortVec[jj].index];
          if(jss.ID == 0) continue;
          TrajPoint& jstp1 = tjs.allTraj[jss.ShowerTjID - 1].Pts[1];
          if(newCuts) {
            // find the longitudinal and transverse separation using the higher energy
            // shower which probably is better defined.
            Point2_t alongTrans;
            FindAlongTrans(istp1.Pos, istp1.Dir, jstp1.Pos, alongTrans);
            // the lower energy shower is at the wrong end of the higher energy shower if alongTrans[0] < 0
            if(alongTrans[0] < 0) continue;
            // increase the cut if the second shower is very low energy (< 10% of the first shower)
            float alongCut = along95;
            if(jss.Energy < 0.1 * iss.Energy) alongCut *= 1.5;
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<fcnLabel<<" Candidate i2S"<<iss.ID<<" j2S"<<jss.ID;
              myprt<<" along "<<std::fixed<<std::setprecision(1)<<alongTrans[0]<<" trans "<<alongTrans[1];
              myprt<<" alongCut "<<alongCut;
            } // prt
            if(alongTrans[0] > alongCut) continue;
          } else {
            // old cuts
            float sep = PosSep(istp1.Pos, jstp1.Pos);
            float trad = sep / radLen;
            // Find the IP between them using the projection of the one with the lowest aspect ratio
            float delta = 9999;
            if(iss.AspectRatio < jss.AspectRatio) {
              delta = PointTrajDOCA(tjs, jstp1.Pos[0], jstp1.Pos[1], istp1);
            } else {
              delta = PointTrajDOCA(tjs, istp1.Pos[0], istp1.Pos[1], jstp1);
            }
            // See if delta is consistent with the cone angle of the i shower
            float dang = delta / sep;
            if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Candidate i2S"<<iss.ID<<" j2S"<<jss.ID<<" separation "<<(int)sep<<" radiation lengths "<<trad<<" delta "<<delta<<" dang "<<dang;
            if(trad > 3) continue;
            // There must be a correlation between dang and the energy of these showers...
            if(dang > 0.3) continue;
          } // old cuts

          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Merge them. Re-find shower center, etc";
          if(MergeShowersAndStore(fcnLabel, tjs, iss.ID, jss.ID, prt)) {
            Trajectory& stj = tjs.allTraj[iss.ShowerTjID - 1];
            stj.AlgMod[kMergeSubShowers] = true;
            keepMerging = true;
            break;
          }
        } // jj
        if(keepMerging) break;
      } // ii
    } // keepMerging
    
  } // MergeSubShowers
  
  ////////////////////////////////////////////////
  int MergeShowers(std::string inFcnLabel, TjStuff& tjs, std::vector<int> showerIDs, bool prt)
  {
    // merge a list of showers and return the ID of the merged shower.
    // Returns 0 if there was a failure. 
    
    std::string fcnLabel = inFcnLabel + ".MS";
    if(showerIDs.size() < 2) return 0;
    // check for a valid ID
    for(auto ssID : showerIDs) if(ssID <= 0 || ssID > (int)tjs.cots.size()) return 0;
    // check for the same CTP
    auto& ss0 = tjs.cots[showerIDs[0] - 1];
    std::vector<int> tjl;
    for(auto ssID : showerIDs) {
      auto& ss = tjs.cots[ssID - 1];
      if(ss.CTP != ss0.CTP) return 0;
      tjl.insert(tjl.end(), ss.TjIDs.begin(), ss.TjIDs.end());
    }
    // ensure the InShower Tjs are valid
    for(auto tjID : tjl) {
      auto& tj = tjs.allTraj[tjID - 1];
      if(tj.CTP != ss0.CTP || tj.AlgMod[kKilled] || !tj.AlgMod[kInShower]) {
        std::cout<<fcnLabel<<" bad InShower Tj "<<tjID<<"\n";
        return 0;
      }
    } // tjID
    
    // mark the old showers killed
    for(auto ssID : showerIDs) {
      auto& ss = tjs.cots[ssID - 1];
      ss.ID = 0;
      // kill the shower Tj
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      stj.AlgMod[kKilled] = true;
    } // tjID

    // in with the new
    int cotID = Create2DShower(tjs, tjl);
    if(cotID == 0) return 0;
    
    // define the new shower
    if(!DefineShower(fcnLabel, tjs, cotID, prt)) {
      std::cout<<fcnLabel<<" DefineShower failed\n";
      MakeShowerObsolete(fcnLabel, tjs, cotID, prt);
      return 0;
    }
    
    return tjs.cots[cotID - 1].ID;
    
  } // MergeShowers
  
  ////////////////////////////////////////////////
  bool MergeShowersAndStore(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, bool prt)
  {
    // Merge showers using shower indices. The icotID shower is modified in-place.
    // The jcotID shower is declared obsolete. This function also re-defines the shower and
    // sets the Parent ID to 0.
    
    if(icotID <= 0 || icotID > tjs.cots.size()) return false;
    ShowerStruct& iss = tjs.cots[icotID - 1];
    if(iss.ID == 0) return false;
    if(iss.TjIDs.empty()) return false;
    if(iss.ShowerTjID <= 0) return false;
    
    if(jcotID <= 0 || jcotID > tjs.cots.size()) return false;
    ShowerStruct& jss = tjs.cots[jcotID - 1];
    if(jss.TjIDs.empty()) return false;
    if(jss.ID == 0) return false;
    if(jss.ShowerTjID <= 0) return false;

    if(iss.CTP != jss.CTP) return false;
    
    std::string fcnLabel = inFcnLabel + ".MSAS";
    
    Trajectory& itj = tjs.allTraj[iss.ShowerTjID - 1];
    Trajectory& jtj = tjs.allTraj[jss.ShowerTjID - 1];
    if(!itj.Pts[1].Hits.empty() || !jtj.Pts[1].Hits.empty()) {
      std::cout<<fcnLabel<<" Warning: These shower Tjs have hits! "<<itj.ID<<" "<<jtj.ID<<"\n";
      return false;
    }
    
    iss.TjIDs.insert(iss.TjIDs.end(), jss.TjIDs.begin(), jss.TjIDs.end());
    // make a new trajectory using itj as a template
    Trajectory ktj = itj;
    ktj.ID = tjs.allTraj.size() + 1;

    tjs.allTraj.push_back(ktj);
    // kill jtj
    MakeTrajectoryObsolete(tjs, iss.ShowerTjID - 1);
    MakeTrajectoryObsolete(tjs, jss.ShowerTjID - 1);
    tjs.allTraj[iss.ShowerTjID - 1].ParentID = ktj.ID;
    tjs.allTraj[jss.ShowerTjID - 1].ParentID = ktj.ID;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" killed stj T"<<iss.ShowerTjID<<" and T"<<jss.ShowerTjID<<" new T"<<ktj.ID;
    // revise the shower
    iss.ShowerTjID = ktj.ID;
    // transfer the list of Tj IDs
    iss.TjIDs.insert(iss.TjIDs.end(), jss.TjIDs.begin(), jss.TjIDs.begin());
    std::sort(iss.TjIDs.begin(), iss.TjIDs.end());
    // clear the list of nearby Tjs
    iss.NearTjIDs.clear();
    // append the list of matched Tjs
    iss.ParentID = 0;
    iss.NeedsUpdate = true;
    jss.ID = 0;
    bool success = DefineShower(fcnLabel, tjs, icotID, prt);
 
    return success;

  } // MergeShowersAndStore
  
  ////////////////////////////////////////////////
  bool MergeShowerTjsAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt)
  {
    // Merge showers using showerTj indices
    // This function is called from MergeAndStore whose function is to merge two line-like
    // trajectories and store them. This function was called because at least one of the
    // trajectories is a shower Tj. Assume that the decision to merge them has been made elsewhere.
    
    if(istj > tjs.allTraj.size() - 1) return false;
    if(jstj > tjs.allTraj.size() - 1) return false;
    
    Trajectory& itj = tjs.allTraj[istj];
    Trajectory& jtj = tjs.allTraj[jstj];
    
    std::string fcnLabel = "MSTJ";
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" MergeShowerTjsAndStore Tj IDs "<<itj.ID<<"  "<<jtj.ID;
    
    // First we check to make sure that both are shower Tjs.
    if(!itj.AlgMod[kShowerTj] && !jtj.AlgMod[kShowerTj]) {
      if(prt) mf::LogVerbatim("TC")<<" One of these isn't a shower Tj";
      return false;
    }
    
    // We need to keep the convention used in MergeAndStore to create a new merged trajectory
    // and kill the two fragments. This doesn't require making a new shower however. We can just
    // re-purpose one of the existing showers
    int icotID = GetCotID(tjs, itj.ID);
    if(icotID == 0) return false;
    ShowerStruct& iss = tjs.cots[icotID - 1];
    if(iss.ID == 0) return false;
    if(iss.TjIDs.empty()) return false;
    int jcotID = GetCotID(tjs, jtj.ID);
    if(jcotID == 0) return false;
    ShowerStruct& jss = tjs.cots[jcotID - 1];
    if(jss.ID == 0) return false;
    if(jss.TjIDs.empty()) return false;
    
    return MergeShowersAndStore(fcnLabel, tjs, icotID, jcotID, prt);
    
  } // MergeShowerTjsAndStore

  ////////////////////////////////////////////////
  bool FindChargeCenter(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Finds the charge center using all sub-structure trajectories in the cot. All of the shower
    // charge is assigned to the second TP and the charge weighted position is put in stp1.HitPos
    // and stp1.Pos
    // The charge will later be distributed between TP0 - TP2.
    // The total charge is stored in  shower Tj AveChg.
    
    if(cotID > tjs.cots.size()) return false;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    
    unsigned short stjIndex = ss.ShowerTjID - 1;
    if(stjIndex > tjs.allTraj.size() - 1) return false;
    if(tjs.allTraj[stjIndex].Pts.size() != 3) return false;
    
    std::string fcnLabel = inFcnLabel + ".FCC";
    
    // Showers with parent Tjs are dealt with elsewhere
    if(ss.ParentID > 0) return true;
    
    // initialize all of the points
    for(unsigned short ii = 0; ii < 3; ++ii) {
      TrajPoint& tp = tjs.allTraj[stjIndex].Pts[ii];
      tp.Chg = 0;
      tp.HitPos[0] = 0;
      tp.HitPos[1] = 0;
    }
    
    TrajPoint& stp1 = tjs.allTraj[stjIndex].Pts[1];
    
    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) {
      stp1.Chg += ss.ShPts[ii].Chg;
      stp1.HitPos[0] += ss.ShPts[ii].Chg * ss.ShPts[ii].Pos[0];
      stp1.HitPos[1] += ss.ShPts[ii].Chg * ss.ShPts[ii].Pos[1];
    } // ii
    
    stp1.HitPos[0] /= stp1.Chg;
    stp1.HitPos[1] /= stp1.Chg;
    // Define the position only if it hasn't been done already
    if(ss.Angle == 0) stp1.Pos = stp1.HitPos;

    // Use the trajectory AveChg variable to store the total charge
    // if it isn't identified as an InShower Tj
    tjs.allTraj[stjIndex].AveChg = stp1.Chg;
    ss.Energy = ShowerEnergy(tjs, ss);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<cotID<<" HitPos "<<(int)stp1.HitPos[0]<<":"<<(int)stp1.HitPos[1]/tjs.UnitsPerTick<<" stp1.Chg "<<(int)stp1.Chg<<" Energy "<<(int)ss.Energy<<" MeV";
    return true;
  } // FindChargeCenter

  ////////////////////////////////////////////////
  void FindAngle(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Find the angle of the shower using the position of all of the TPs
    
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    if(ss.ParentID > 0) return;
    
    unsigned short stjIndex = ss.ShowerTjID - 1;
    if(stjIndex > tjs.allTraj.size() - 1) return;
    if(tjs.allTraj[stjIndex].Pts.size() != 3) return;
    
    std::string fcnLabel = inFcnLabel + ".FA";

    TrajPoint& stp1 = tjs.allTraj[stjIndex].Pts[1];
    
    // Do a least squares fit using all the points
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    double xx, yy, xr, yr, wt;

    // rotate into a coordinate system along the current shower axis
    // TODO: Use and replace RotPos?
    double cs = cos(-ss.Angle);
    double sn = sin(-ss.Angle);
    
    unsigned short nptsFit = 0;
    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) {
      // Weight by charge
      wt = ss.ShPts[ii].Chg;
      sum  += wt;
      xx = wt * (ss.ShPts[ii].Pos[0] - stp1.Pos[0]);
      yy = wt * (ss.ShPts[ii].Pos[1] - stp1.Pos[1]);
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      sumx += wt * xr;
      sumy += wt * yr;
      sumx2 += wt * xr * xr;
      sumy2 += wt * yr * yr;
      sumxy += wt * xr * yr;
      ++nptsFit;
    } // ii
    
    if(nptsFit < 3) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Not enough points to fit";
      return;
    }
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0) {
      ss.Angle = 0;
      ss.AngleErr = 1.5;
      return;
    }
    // A is the intercept (This should be ~0 if the charge center is known)
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // B is the slope
    double B = (sumxy * sum  - sumx * sumy) / delta;
    float dang = atan(B);
    ss.Angle += dang;
/* TODO: This gives errors that are much too small
    double ndof = nptsFit - 1;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    ss.AngleErr = sqrt(varnce * sum / delta);
*/    
    // Fake a correction to the angle error instead
    ss.AngleErr = ss.AspectRatio * ss.AspectRatio * 1.57;
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<cotID<<" Pos "<<ss.CTP<<":"<<PrintPos(tjs, stp1)<<" Intercept "<<(int)A<<" dang "<<dang<<" Angle "<<ss.Angle<<" Err "<<ss.AngleErr<<" npts fit "<<nptsFit;
    
  } // FindAngle
  
  ////////////////////////////////////////////////
  void FillRotPos(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Fills the RotPos vector and sorts the points along the shower axis. Note that the rotation is
    // done around stp1.Pos but the charge center is at stp1.HitPos. Pos and HitPos will be exactly the
    // same if there is no parent. The Pos position may be shifted slightly in FindExternalParent so that
    // the parent trajectory lies on the central axis of the shower. This is done so that the charge at the
    // start of the shower is calculated correctly using the parent trajectory points
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return;
    
    std::string fcnLabel = inFcnLabel + ".FRP";

    // Determine the size of the shower along the axis and transverse to it. 
    // Rotate and translate each point into the coordinate system defined by tp[1]
    float cs = cos(-ss.Angle);
    float sn = sin(-ss.Angle);

    TrajPoint& stp1 = stj.Pts[1];

    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) {
      ss.ShPts[ii].RotPos[0] = ss.ShPts[ii].Pos[0] - stp1.Pos[0];
      ss.ShPts[ii].RotPos[1] = ss.ShPts[ii].Pos[1] - stp1.Pos[1];
      // Rotate into the stp1 direction
      float along = cs * ss.ShPts[ii].RotPos[0] - sn * ss.ShPts[ii].RotPos[1];
      float trans = sn * ss.ShPts[ii].RotPos[0] + cs * ss.ShPts[ii].RotPos[1];
      ss.ShPts[ii].RotPos[0] = along;
      ss.ShPts[ii].RotPos[1] = trans;
    } // ii
    
    std::vector<SortEntry> sortVec(ss.ShPts.size());
    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].length = ss.ShPts[ii].RotPos[0];
    }
    std::sort(sortVec.begin(), sortVec.end(), lessThan);
    
    // put the points vector into the sorted order
    auto tPts = ss.ShPts;
    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      ss.ShPts[ii] = tPts[indx];
    } // ii

    // Calculate the the aspect ratio
    float alongSum = 0;
    float transSum = 0;
    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) {
      alongSum += std::abs(ss.ShPts[ii].RotPos[0]);
      transSum += std::abs(ss.ShPts[ii].RotPos[1]);
    } // ii
    
    if(alongSum > 0) {
      ss.AspectRatio = transSum / alongSum;
    } else {
      ss.AspectRatio = 100;
    }
    // Fake a correction to the angle error
    ss.AngleErr = ss.AspectRatio * ss.AspectRatio * 1.57;
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  2S"<<ss.ID<<" Rotation origin "<<PrintPos(tjs, stp1.Pos)<<" Angle "<<std::setprecision(2)<<ss.Angle<<" AspectRatio "<<ss.AspectRatio<<" AngleErr "<<ss.AngleErr;
    
  } // FillRotPos
  
  ////////////////////////////////////////////////
  bool AnalyzeRotPos(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // The RotPos vector was filled and sorted by increasing distance along the shower axis in FillRotPos.
    // This function divides the RotPos points into 3 sections and puts the transverse rms width in the
    // three sections into the shower Tj TrajPoint DeltaRMS variable. It also calculates the charge and number of shower
    // points closest to each TrajPoint. It make some crude quality cuts and returns false if the shower
    // fails these cuts
    
    if(cotID > tjs.cots.size()) return false;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return false;
    
    std::string fcnLabel = inFcnLabel + ".ARP";
    
    for(auto& tp : stj.Pts) {
      tp.Chg = 0;
      tp.DeltaRMS = 0;
      tp.NTPsFit = 0;
      tp.HitPos = {{0.0, 0.0}};
    }
    
    float minAlong = ss.ShPts[0].RotPos[0];
    float maxAlong = ss.ShPts[ss.ShPts.size()-1].RotPos[0];
    float sectionLength = (maxAlong - minAlong) / 3;
    float sec0 = minAlong + sectionLength;
    float sec2 = maxAlong - sectionLength;
    // iterate over the shower points (aka hits)
    for(auto& spt : ss.ShPts) {
      // The point on the shower Tj to which the charge will assigned
      unsigned short ipt = 0;
      if(spt.RotPos[0] < sec0) {
        // closest to point 0
        ipt = 0;
      } else if(spt.RotPos[0] > sec2) {
        // closest to point 2
        ipt = 2;
      } else {
        // closest to point 1
        ipt = 1;
      }
      stj.Pts[ipt].Chg += spt.Chg;
      // Average the absolute value of the transverse position in lieu of
      // using the sum of the squares. The result is ~15% higher than the actual
      // rms which is OK since this is used to find the transverse size of the shower
      // which is not a precisely defined quantity anyway
      stj.Pts[ipt].DeltaRMS += spt.Chg * std::abs(spt.RotPos[1]);
      ++stj.Pts[ipt].NTPsFit;
      // Average the charge center at each point
      stj.Pts[ipt].HitPos[0] += spt.Chg * spt.Pos[0];
      stj.Pts[ipt].HitPos[1] += spt.Chg * spt.Pos[1];
    } // spt
    
    for(auto& tp : stj.Pts) {
      if(tp.Chg > 0) {
        tp.DeltaRMS /= tp.Chg;
        tp.HitPos[0] /= tp.Chg;
        tp.HitPos[1] /= tp.Chg;
      }
    } // tp
    
    // require that there is charge in point 0 and 2. Point 1 may not have charge if
    // we are constructing a sparse shower that is not yet well-defined
    if(stj.Pts[0].Chg == 0 || stj.Pts[2].Chg == 0) return false;
    
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) {
      // do a simple interpolation
      stj.Pts[1].HitPos[0] = 0.5 * (stj.Pts[0].HitPos[0] + stj.Pts[2].HitPos[0]);
      stj.Pts[1].HitPos[1] = 0.5 * (stj.Pts[0].HitPos[1] + stj.Pts[2].HitPos[1]);
    }
    if(stj.Pts[2].DeltaRMS > 0) {
      ss.DirectionFOM = stj.Pts[0].DeltaRMS / stj.Pts[2].DeltaRMS;
    } else {
      ss.DirectionFOM = 10;
    }
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 2S"<<ss.ID;
      myprt<<" HitPos "<<std::fixed<<std::setprecision(1);
      myprt<<stj.Pts[1].HitPos[0]<<" "<<stj.Pts[1].HitPos[1]<<" "<<stj.Pts[1].HitPos[2];
      myprt<<" DeltaRMS "<<std::setprecision(2);
      myprt<<stj.Pts[0].DeltaRMS<<" "<<stj.Pts[1].DeltaRMS<<" "<<stj.Pts[2].DeltaRMS;
      myprt<<" DirectionFOM "<<std::fixed<<std::setprecision(2)<<ss.DirectionFOM;
    }
    return true;

  } // AnalyzeRotPos
  
  ////////////////////////////////////////////////
  bool DefineShowerTj(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Defines the Shower Tj, calculates the shower aspect ratio, etc. This function
    // doesn't change the state of Parent
    
    if(cotID > tjs.cots.size()) return false;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    
    // don't alter showers with parents. Assume that the shower is defined
    if(ss.ParentID > 0) return true;

    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return false;
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return false;
    
    std::string fcnLabel = inFcnLabel + ".DSTj";
    
    // This function calculates the charge (Chg), the shower Tj points (HitPos) and
    // the transverse shower rms at each Tp in the shower Tj. This is used below to determine
    // the shower direction when a parent has not been identified. 
    if(!AnalyzeRotPos(fcnLabel, tjs, cotID, prt)) return false;

    // startsNeg is true if this assumption is correct
    bool startsNeg = (stj.Pts[0].DeltaRMS < stj.Pts[2].DeltaRMS);
    
    // reverse the points vector so that the narrow end of the shower is near Pts.begin()
    if(!startsNeg) ReverseShower(fcnLabel, tjs, cotID, prt);

    // define the shower start and end positions
    stj.Pts[0].Pos = ss.ShPts[0].Pos;
    unsigned short endPt = ss.ShPts.size()-1;
    stj.Pts[2].Pos = ss.ShPts[endPt].Pos;
    // and the charge center
    stj.Pts[1].Pos = stj.Pts[1].HitPos;

    // define the angle of all the shower Tps
    for(auto& stp : stj.Pts) {
      stp.Ang = ss.Angle;
      stp.AngErr = ss.AngleErr;
      stp.Dir[0] = cos(stp.Ang);
      stp.Dir[1] = sin(stp.Ang);
    } // stp
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" stj "<<stj.ID<<" Pos "<<PrintPos(tjs, stj.Pts[0].Pos)<<" "<<PrintPos(tjs, stj.Pts[1].Pos)<<" "<<PrintPos(tjs, stj.Pts[2].Pos);
    // Guard against too-narrow shower widths at the charge center
    if(stj.Pts[1].DeltaRMS < stj.Pts[0].DeltaRMS || stj.Pts[1].DeltaRMS < stj.Pts[2].DeltaRMS) {
      // take the average
      stj.Pts[1].DeltaRMS = 0.5 * (stj.Pts[0].DeltaRMS + stj.Pts[2].DeltaRMS);
    }

    // use Tj Pass to count the number of times this function is used
    ++stj.Pass;
    
    return true;

  } // DefineShowerTj
  
  ////////////////////////////////////////////////
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Reverses the shower and the shower tj
    
    if(cotID > tjs.cots.size()) return;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".RevSh";
    
    std::reverse(ss.ShPts.begin(), ss.ShPts.end());
    // change the sign of RotPos
    for(auto& sspt : ss.ShPts) {
      sspt.RotPos[0] = -sspt.RotPos[0];
      sspt.RotPos[1] = -sspt.RotPos[1];
    }
    // flip the shower angle
    if(ss.Angle > 0) {
      ss.Angle -= M_PI;
    } else {
      ss.Angle += M_PI;
    }
    ss.DirectionFOM = 1 / ss.DirectionFOM;
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    ReverseTraj(tjs, stj);
    DefineEnvelope(fcnLabel, tjs, cotID, prt);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Reversed shower. Shower angle = "<<ss.Angle;
  }
  ////////////////////////////////////////////////
  void RefineShowerTj(TjStuff& tjs, int cotID, bool prt)
  {
    // Checks the properties of Shower Tj and revises them if necessary. Returns true if the
    // shower needs to be updated
    
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    // no refinement necessary
    if(ss.ParentID > 0) return;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return;
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return;
    
    // Ignore fat showers
    if(ss.AspectRatio > tjs.ShowerTag[10] || ss.DirectionFOM > tjs.ShowerTag[9]) {
      if(prt) mf::LogVerbatim("TC")<<"RSTj: Not possible due to poor AspectRatio "<<ss.AspectRatio<<" or ss.DirectionFOM "<<ss.DirectionFOM;
      return;
    }

    if(ss.ParentID > 0) {
      if(prt) mf::LogVerbatim("TC")<<"RSTj: Use existing Parent Tj "<<ss.ParentID<<" info.";
      return;
    }
    
    // check the beginning of the shower to see if the points are close to the shower axis
    float sum = 0;
    float sum2 = 0;
    float chg = 0;
    unsigned short cnt = 0;
    // check the first 1/3 of the way along the shower
    unsigned short lastPt = ss.ShPts.size() - 1;
    float maxAlong = ss.ShPts[0].RotPos[0] + 0.3 * (ss.ShPts[lastPt].RotPos[0] - ss.ShPts[0].RotPos[0]);
    unsigned short firstTID = ss.ShPts[0].TID;
    unsigned short cntTID = 0;
    for(auto& sspt : ss.ShPts) {
      if(sspt.RotPos[0] > maxAlong) break;
      chg += sspt.Chg;
      sum += sspt.Chg * sspt.RotPos[1];
      sum2 += sspt.Chg * sspt.RotPos[1] * sspt.RotPos[1];
      if(prt && cnt < 10) mf::LogVerbatim("TC")<<"RSTj "<<sspt.Pos[0]<<":"<<(int)sspt.Pos[1]<<" RP "<<(int)sspt.RotPos[0]<<" "<<sspt.RotPos[1]<<" Chg "<<(int)sspt.Chg<<" TID "<< sspt.TID;
      ++cnt;
      if(sspt.TID == firstTID) ++cntTID;
    } // sspt
    if(chg == 0) return;
    float transAve = sum / chg;
    float arg = sum2 - chg * transAve * transAve;
    if(arg == 0) return;
    float transRMS = sqrt(arg / chg);
    transRMS /= sqrt((float)cnt);

    if(prt) mf::LogVerbatim("TC")<<"RSTj shower begin transAve "<<transAve<<" transRMS "<<transRMS;

  } // RefineShowerTj
  
  ////////////////////////////////////////////////
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Gracefully kills the shower and the associated shower Tj
    
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".MSO";
    
    if(ss.SS3ID > 0) {
      std::cout<<fcnLabel<<" Trying to kill shower "<<ss.ID<<" that is matched to a 3D shower...\n";
    }
    
    // Kill the shower Tj if it exists. This also releases the hits
    if(ss.ShowerTjID > 0) MakeTrajectoryObsolete(tjs, ss.ShowerTjID - 1);
    
    // Restore the original InShower Tjs
    // Unset the killed bit
    for(auto& tjID : ss.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjID - 1];
      tj.AlgMod[kKilled] = false;
      // clear all of the shower-related bits
      tj.AlgMod[kInShower] = false;
      tj.AlgMod[kShwrParent] = false;
      tj.AlgMod[kMergeOverlap] = false;
      tj.AlgMod[kMergeSubShowers] = false;
      tj.AlgMod[kMergeNrShowers] = false;
      tj.AlgMod[kMergeShChain] = false;
     // Restore the hit -> tj association. This is strictly only necessary if the
      // hits were re-assigned to the shower Tj but do it anyway just to be sure
      for(auto& tp : tj.Pts) {
        for(auto& iht : tp.Hits) tjs.fHits[iht].InTraj = tj.ID;
      } // tp
    } // tjID
    ss.ID = 0;
    ss.TjIDs.clear();
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Killed ShowerTj "<<ss.ShowerTjID<<" and restored InShower Tjs.";
    
  } // MakeShowerObsolete
  
  ////////////////////////////////////////////////
  void TagShowerLike(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList, bool applyMinTjCuts)
  {
    // Tag Tjs as InShower if they have MCSMom < ShowerTag[0] and there are more than
    // ShowerTag[6] other Tjs with a separation < ShowerTag[1]. Returns a list of Tjs that meet this criteria.
    // The shower cuts that are applicable are done here if applyMinTjCuts is true, for example when this function is called
    // from RunTrajClusterAlg before 3D vertex and tj matching is done. This reduces the number of spurious vertices and 3D matches.
    // When called from FindShowers3D, applyMinTjCuts is set false and the cuts are applied after 2D showers are reconstructed.
    
    tjList.clear();
    
    if(tjs.ShowerTag[0] <= 0) return;
    
    if(tjs.allTraj.size() > 20000) return;
    
    // evaluate different cuts
    bool newCuts = (tjs.ShowerTag[0] == 3);
    float typicalChgRMS = 0.5 * (tjs.ChargeCuts[1] + tjs.ChargeCuts[2]);
    
    // clear out old tags and make a list of Tjs to consider
    std::vector<int> tjids;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      tj.AlgMod[kShowerLike] = false;
      if(tj.AlgMod[kShowerTj]) continue;
      short npwc = NumPtsWithCharge(tjs, tj, false);
      if(newCuts) {
        // evaluate different cuts
        // Don't expect any (primary) electron to be reconstructed as a single trajectory for
        // more than ~2 radiation lengths ~ 30 cm for uB ~ 100 wires 
        if(npwc > 100) continue;
        // allow short Tjs.
        if(npwc > 5) {
          // Increase the MCSMom cut if the Tj is long and the charge RMS is high to reduce sensitivity 
          // to the fcl configuration. A primary electron may be reconstructed as one long Tj with large
          // charge rms and possibly high MCSMom or as several nearby shorter Tjs with lower charge rms
          float momCut = tjs.ShowerTag[1];
          if(tj.ChgRMS > typicalChgRMS) momCut *= tj.ChgRMS / typicalChgRMS;
          if(tj.MCSMom > momCut) continue;
        }
      } else {
        if(npwc < 3) continue;
        if(npwc > 4 && tj.MCSMom > tjs.ShowerTag[1]) continue;
      }
      tjids.push_back(tj.ID);
    } // tj
    
    if(tjids.size() < 2) return;
/*
    std::cout<<inCTP<<"TIST tjids";
    for(auto tjid : tjids) std::cout<<" T"<<tjid;
    std::cout<<"\n";
*/
    for(unsigned short it1 = 0; it1 < tjids.size() - 1; ++it1) {
      Trajectory& tj1 = tjs.allTraj[tjids[it1] - 1];
      float len1 = PosSep(tj1.Pts[tj1.EndPt[1]].Pos, tj1.Pts[tj1.EndPt[0]].Pos);
      for(unsigned short it2 = it1 + 1; it2 < tjids.size(); ++it2) {
        Trajectory& tj2 = tjs.allTraj[tjids[it2] - 1];
        unsigned short ipt1, ipt2;
        float doca = tjs.ShowerTag[2];
        // Find the separation between Tjs without considering dead wires
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, doca, false);
        if(doca == tjs.ShowerTag[2]) continue;
        // make tighter cuts for user-defined short Tjs
        float len2 = PosSep(tj2.Pts[tj2.EndPt[1]].Pos, tj2.Pts[tj2.EndPt[0]].Pos);
//        if(doca < 20) std::cout<<inCTP<<" T"<<tj1.ID<<" len "<<(int)len1<<" T"<<tj2.ID<<" len "<<(int)len2<<" doca "<<doca<<"\n"; 
        if(!newCuts) {
          if(len1 < len2 && len1 < doca) {
            if(len1 < doca) continue;
          } else {
            if(len2 < doca) continue;
          }
        } // !newCuts
        // found a close pair. See if one of these is in an existing cluster of Tjs
        bool inlist = false;
        for(unsigned short it = 0; it < tjList.size(); ++it) {
          bool tj1InList = (std::find(tjList[it].begin(), tjList[it].end(), tj1.ID) != tjList[it].end());
          bool tj2InList = (std::find(tjList[it].begin(), tjList[it].end(), tj2.ID) != tjList[it].end());
          if(tj1InList || tj2InList) {
            // add the one that is not in the list
            if(!tj1InList) tjList[it].push_back(tj1.ID);
            if(!tj2InList) tjList[it].push_back(tj2.ID);
            inlist = true;
            break;
          }
          if(inlist) break;
        } // it
        // start a new list with this pair?
        if(!inlist) {
          std::vector<int> newlist(2);
          newlist[0] = tj1.ID;
          newlist[1] = tj2.ID;
          tjList.push_back(newlist);
        }
      } // it2
    } // it1
    if(tjList.empty()) return;
/*
    std::cout<<"tjLists\n";
    for(auto& tjl : tjList) {
      for(auto tjid : tjl) std::cout<<" T"<<tjid;
      std::cout<<"\n";
    } // tjl
*/
    MergeTjList(tjList);
    
    // eliminate entries that fail ShowerTag[6]
    std::vector<std::vector<int>> newList;
    if(newCuts) {
      for(auto& tjl : tjList) {
        float npts = 0;
        for(auto tjid : tjl) {
          auto& tj = tjs.allTraj[tjid - 1];
          npts += NumPtsWithCharge(tjs, tj, false);
        } // tjid
        if(npts >= tjs.ShowerTag[6]) newList.push_back(tjl);
      } // tjl
    } else {
      // old cuts
      for(auto& tjl : tjList) if(tjl.size() >= tjs.ShowerTag[7]) newList.push_back(tjl);
    } // old cuts
    tjList = newList;
/*
    std::cout<<"tjLists final\n";
    for(auto& tjl : tjList) {
      for(auto tjid : tjl) std::cout<<" T"<<tjid;
      std::cout<<"\n";
    } // tjl
*/
    // mark them all as ShowerLike Tjs
    unsigned short nsh = 0;
    for(auto& tjl : tjList) {
      if(applyMinTjCuts && tjl.size() < tjs.ShowerTag[7]) continue;
      nsh += tjl.size();
      for(auto& tjID : tjl) {
        auto& tj = tjs.allTraj[tjID - 1];
        tj.AlgMod[kShowerLike] = true;
        // unset flags
        tj.AlgMod[kSetDir] = false;
        for(unsigned short end = 0; end < 2; ++end) tj.StopFlag[end][kBragg] = false;
      } // tjid
    } // tjl
    if(tjs.ShowerTag[12] >= 0) mf::LogVerbatim("TC")<<"TagShowerLike tagged "<<nsh<<" Tjs in CTP "<<inCTP;
    // Set the NearInShower bit on all Tjs in this TPC
    std::vector<unsigned int> closeHits;
    // Set the 
    std::array<int, 2> wireWindow;
    Point2_t timeWindow;
    // Define close to be 1/2 of the ShowerTag separation cut
    float deltaCut = 0.5 * tjs.ShowerTag[2];
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Pos[0] < -0.4) continue;
        // Find all hits near this time on this wire
        wireWindow[0] = std::nearbyint(tp.Pos[0]);
        wireWindow[1] = wireWindow[0];
        timeWindow[0] = tp.Pos[1] - deltaCut;
        timeWindow[1] = tp.Pos[1] + deltaCut;
        bool hitsNear = false;
        closeHits = FindCloseHits(tjs, wireWindow, timeWindow, plane, kAllHits, true, hitsNear);
        for(auto iht : closeHits) {
          auto& hit = tjs.fHits[iht];
          // ignore hits used in this TP
          if(hit.InTraj == tj.ID) continue;
          if(hit.InTraj > 0) {
            // Nearby hit used in a nearby Tj
            tp.Environment[kEnvNearTj] = true;
            auto& ntj = tjs.allTraj[hit.InTraj - 1];
            if(ntj.AlgMod[kInShower]) tp.Environment[kEnvNearShower] = true;
          }
        } // iht
      } // ipt
    } // tj
    
    if(tjs.UseAlg[kKillInShowerVx]) {
      // make a list of 2D vertices
      std::vector<unsigned short> vxids;
      for(auto& tjl : tjList) {
        for(auto& tjID : tjl) {
          auto& tj = tjs.allTraj[tjID - 1];
          for(unsigned short end = 0; end < 2; ++end) {
            if(!tj.AlgMod[kShowerLike]) continue;
            if(tj.VtxID[end] == 0) continue;
            if(std::find(vxids.begin(), vxids.end(), tj.VtxID[end]) == vxids.end()) vxids.push_back(tj.VtxID[end]);
          }
        } // tjid
      } // tjl
      if(vxids.empty()) return;
      for(auto vxid : vxids) {
        auto& vx2 = tjs.vtx[vxid - 1];
        // already killed?
        if(vx2.ID == 0) continue;
        // don't kill a primary vertex
        if(vx2.Vx3ID > 0) {
          auto& vx3 = tjs.vtx3[vx2.ID - 1];
          if(vx3.Primary) continue;
        }
        // get a list of Tjs attached to this vertex
        auto vxtjs = GetVtxTjIDs(tjs, vx2);
        // count the number that are InShower
        unsigned short ninsh = 0;
        for(auto tjid : vxtjs) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kShowerLike]) ++ninsh;
        } // tjid
        // Jan 22
//        if(ninsh > 1) MakeVertexObsolete(tjs, vx2, false);
        if(ninsh > 1) MakeVertexObsolete(tjs, vx2, true);
      } // vxid
    } // 
    
  } // TagShowerLike
  
  ////////////////////////////////////////////////
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Find Tjs that are near the shower but are not included in it
    if(cotID > tjs.cots.size()) return;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    ss.NearTjIDs.clear();
    
    // check for a valid envelope
    if(ss.Envelope.empty()) return;
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    
    std::string fcnLabel = inFcnLabel + ".FNTj";
    
    std::vector<int> ntj;
    
    // set max distance of closest approach ~ 5 radiation lengths ~200 WSE 
    constexpr float fiveRadLen = 200;
    
    // look for vertices inside the envelope
    for(auto vx : tjs.vtx) {
      if(vx.CTP != ss.CTP) continue;
      if(vx.ID == 0) continue;
      if(!PointInsideEnvelope(vx.Pos, ss.Envelope)) continue;
      auto vxTjIDs = GetVtxTjIDs(tjs, vx);
      for(auto tjID : vxTjIDs) {
        // ignore those that are in the shower
        if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjID) != ss.TjIDs.end()) continue;
        // or already in the list
        if(std::find(ntj.begin(), ntj.end(), tjID) != ntj.end()) continue;
        ntj.push_back(tjID);
      } // tjID
    } // vx

    // Check for tj points inside the envelope
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      // not a showerTj
      if(tj.AlgMod[kShowerTj]) continue;
      // make sure it's not in the shower
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
      // or already in the list
      if(std::find(ntj.begin(), ntj.end(), tj.ID) != ntj.end()) continue;
      // check proximity of long high MCSMom Tjs to the shower center
      if(tj.Pts.size() > 40 && tj.MCSMom > 200) {
        float delta = PointTrajDOCA(tjs, stj.Pts[1].Pos[0], stj.Pts[1].Pos[1], tj.Pts[tj.EndPt[0]]);
        // TODO: This could be done much better
        if(delta < 20) {
          float doca = fiveRadLen;
          unsigned short spt = 0, ipt = 0;
          TrajTrajDOCA(tjs, stj, tj, spt, ipt, doca);
          if(doca < fiveRadLen) {
            ntj.push_back(tj.ID);
            continue;
          }
        }
      } // long hi-MCSMom tj
      // don't need to check every point. Every third should be enough
      bool isInside = false;
      for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ipt += 3) {
        if(PointInsideEnvelope(tj.Pts[ipt].Pos, ss.Envelope)) {
          isInside = true;
          break;
        }
      } // ipt
      // check the last point which was probably missed above
      if(!isInside) {
        unsigned short ipt = tj.EndPt[1];
        isInside = PointInsideEnvelope(tj.Pts[ipt].Pos, ss.Envelope);
      }
      if(isInside) ntj.push_back(tj.ID);
    } // tj
    if(ntj.size() > 1) std::sort(ntj.begin(), ntj.end());
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Found "<<ntj.size()<<" Tjs near ss.ID "<<ss.ID;
    ss.NearTjIDs = ntj;
    
  } // FindNearbyTjs
  
  ////////////////////////////////////////////////
  void AddCloseTjsToList(TjStuff& tjs, unsigned short itj, std::vector<int> list)
  {
    // Searches the trajectory points for hits that are used in a different trajectory and add
    // them to the list if any are found, and the MCSMomentum is not too large
    if(itj > tjs.allTraj.size() - 1) return;
    
    //short maxMom = (short)(2 * tjs.ShowerTag[1]);
    short maxMom = tjs.ShowerTag[1];
    //XL: why is maxMom is twice of the shower tag [1]? 
    for(auto& tp : tjs.allTraj[itj].Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        // ignore hits that are used in this trajectory
        if(tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        // ignore if there is no hit -> Tj association
        if(tjs.fHits[iht].InTraj <= 0) continue;
        if((unsigned int)tjs.fHits[iht].InTraj > tjs.allTraj.size()) continue;
        // check the momentum
        Trajectory& tj = tjs.allTraj[tjs.fHits[iht].InTraj - 1];
        if(tj.MCSMom > maxMom) continue;
//        if(tj.AlgMod[kTjHiVx3Score]) continue;
        // see if it is already in the list
        if(std::find(list.begin(), list.end(), tjs.fHits[iht].InTraj) != list.end()) continue;
        list.push_back(tjs.fHits[iht].InTraj);
      } // ii
    } // tp
  } // AddCloseTjsToList

  ////////////////////////////////////////////////
  void DefineEnvelope(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".DE";
    
    ss.Envelope.resize(4);
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    TrajPoint& stp0 = stj.Pts[0];
    TrajPoint& stp1 = stj.Pts[1];
    TrajPoint& stp2 = stj.Pts[2];
    
    // construct the Envelope polygon. Start with a rectangle using the fixed 1/2 width fcl input
    // expanded by the rms width at each end to create a polygon. The polygon is constructed along
    // the Pos[0] direction and then rotated into the ShowerTj direction. Use sTp1 as the origin.
    // First vertex
    ss.Envelope[0][0] = -PosSep(stp0.Pos, stp1.Pos);
    ss.Envelope[0][1] = tjs.ShowerTag[5] + tjs.ShowerTag[4] * stp0.DeltaRMS;
    // second vertex
    ss.Envelope[1][0] = PosSep(stp1.Pos, stp2.Pos);
    ss.Envelope[1][1] = tjs.ShowerTag[5] + tjs.ShowerTag[4] * stp2.DeltaRMS;
    // third and fourth are reflections of the first and second
    ss.Envelope[2][0] =  ss.Envelope[1][0];
    ss.Envelope[2][1] = -ss.Envelope[1][1];
    ss.Envelope[3][0] =  ss.Envelope[0][0];
    ss.Envelope[3][1] = -ss.Envelope[0][1];
    
    float length = ss.Envelope[1][0] - ss.Envelope[0][0];
    float width = ss.Envelope[0][1] + ss.Envelope[1][1];
    ss.EnvelopeArea = length * width;

    // Rotate into the stp1 coordinate system
    float cs = cos(stp1.Ang);
    float sn = sin(stp1.Ang);
    for(auto& vtx : ss.Envelope) {
      // Rotate along the stj shower axis
      float pos0 = cs * vtx[0] - sn * vtx[1];
      float pos1 = sn * vtx[0] + cs * vtx[1];
      // translate
      vtx[0] = pos0 + stp1.Pos[0];
      vtx[1] = pos1 + stp1.Pos[1];
    } // vtx
    // Find the charge density inside the envelope
    ss.ChgDensity = (stp0.Chg + stp1.Chg + stp2.Chg) / ss.EnvelopeArea;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 2S"<<cotID<<" Envelope";
      for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      myprt<<" Area "<<(int)ss.EnvelopeArea;
      myprt<<" ChgDensity "<<ss.ChgDensity;
    }
    // This is the last function used to update a shower
    ss.NeedsUpdate = false;
  } // DefineEnvelope  
  
  ////////////////////////////////////////////////
  bool AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
   {
    // This function adds Tjs to the shower. It updates the shower parameters.
    
     if(cotID > tjs.cots.size()) return false;
    
     ShowerStruct& ss = tjs.cots[cotID - 1];
     if(ss.Envelope.empty()) return false;
     if(ss.ID == 0) return false;
     if(ss.TjIDs.empty()) return false;
     
     std::string fcnLabel = inFcnLabel + ".ATIE";

     if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Checking 2S"<<ss.ID;

     unsigned short nadd = 0;
     for(auto& tj : tjs.allTraj) {
       if(tj.CTP != ss.CTP) continue;
       if(tj.AlgMod[kKilled]) continue;
       if(tj.AlgMod[kInShower]) continue;
       if(tj.AlgMod[kShowerTj]) continue;
       // See if this Tjs is attached to a neutrino vertex. 
       if(tj.ParentID == 0) continue;
       int neutPrimTj = NeutrinoPrimaryTjID(tjs, tj);
       if(neutPrimTj > 0 && neutPrimTj != tj.ID) {
         // The Tj is connected to a primary Tj that is associated with a neutrino primary.
         // Don't allow tjs to be added to the shower that are not connected to this neutrino primary (if
         // one exists)
         if(ss.ParentID > 0 && neutPrimTj != ss.ParentID) continue;
       } // neutrino primary tj
       // This shouldn't be necessary but ensure that the Tj ID appears only once in ss.TjIDs
       if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
       // See if both ends are outside the envelope
       bool end0Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[0]].Pos, ss.Envelope);
       bool end1Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[1]].Pos, ss.Envelope);
       if(!end0Inside && !end1Inside) continue;
       if(end0Inside && end1Inside) {
         // TODO: See if the Tj direction is compatible with the shower?
         if(AddTj(fcnLabel, tjs, tj.ID, cotID, false, prt)) ++nadd;
         ++nadd;
         continue;
       } // both ends inside
       // Require high momentum Tjs be aligned with the shower axis
       // TODO also require high momentum Tjs close to the shower axis?

       if(tj.MCSMom > 200) {
         float tjAngle = tj.Pts[tj.EndPt[0]].Ang;
         float dangPull = std::abs(tjAngle - ss.AngleErr) / ss.AngleErr;
         if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" high MCSMom "<<tj.MCSMom<<" dangPull "<<dangPull;
         if(dangPull > 2) continue;
       } // high momentum
       if(AddTj(fcnLabel, tjs, tj.ID, cotID, false, prt)) {
         ++nadd;
       } else {
         if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" AddTj failed to add T"<<tj.ID;
       }
    } // tj
    
    if(nadd > 0) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added "<<nadd<<" trajectories ";
      ss.NeedsUpdate = true;
      if(ss.ParentID == 0) DefineShower(fcnLabel, tjs, cotID, prt);
      return true;
    } else {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" No new trajectories added to envelope ";
      ss.NeedsUpdate = false;
      return false;
    }
        
  } // AddTjsInsideEnvelope
  
  ////////////////////////////////////////////////
  bool AddLooseHits(TjStuff& tjs, int cotID, bool prt)
  {
    // Add hits that are inside the envelope to the shower if they are loose, i.e. not
    // used by any trajectory. This function returns true if the set of hits is different than
    // the current set. The calling function should update the shower if this is the case.
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.Envelope.empty()) return false;
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;

    geo::PlaneID planeID = DecodeCTP(ss.CTP);
    unsigned short ipl = planeID.Plane;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    TrajPoint& stp0 = stj.Pts[0];
    // start a list of new hits
    std::vector<unsigned int> newHits;
    
    // look for hits inside the envelope. Find the range of wires that spans the envelope
    float fLoWire = 1E6;
    float fHiWire = 0;
    // and the range of ticks
    float loTick = 1E6;
    float hiTick = 0;
    for(auto& vtx : ss.Envelope) {
      if(vtx[0] < fLoWire) fLoWire = vtx[0];
      if(vtx[0] > fHiWire) fHiWire = vtx[0];
      if(vtx[1] < loTick) loTick = vtx[1];
      if(vtx[1] > hiTick) hiTick = vtx[1];
    } // vtx
    loTick /= tjs.UnitsPerTick;
    hiTick /= tjs.UnitsPerTick;
    unsigned int loWire = std::nearbyint(fLoWire);
    unsigned int hiWire = std::nearbyint(fHiWire) + 1;
    if(hiWire > tjs.LastWire[ipl]-1) hiWire = tjs.LastWire[ipl]-1;
    std::array<float, 2> point;
    for(unsigned int wire = loWire; wire < hiWire; ++wire) {
      // skip bad wires or no hits on the wire
      if(tjs.WireHitRange[ipl][wire].first < 0) continue;
      unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
      unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
      for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
        // used in a trajectory?
        if(tjs.fHits[iht].InTraj != 0) continue;
        // inside the tick range?
        if(tjs.fHits[iht].PeakTime < loTick) continue;
        // Note that hits are sorted by increasing time so we can break here
        if(tjs.fHits[iht].PeakTime > hiTick) break;
        // see if this hit is inside the envelope
        point[0] = tjs.fHits[iht].ArtPtr->WireID().Wire;
        point[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
        if(!PointInsideEnvelope(point, ss.Envelope)) continue;
        newHits.push_back(iht);
      } // iht
    } // wire
    
    // no new hits and no old hits. Nothing to do
    if(newHits.empty()) {
      if(prt) mf::LogVerbatim("TC")<<"ALH: No new loose hits found";
      return false;
    }
    
    // Update
    stp0.Hits.insert(stp0.Hits.end(), newHits.begin(), newHits.end());
    for(auto& iht: newHits) tjs.fHits[iht].InTraj = stj.ID;
    
    if(prt) mf::LogVerbatim("TC")<<"ALH: Added "<<stp0.Hits.size()<<" hits to stj "<<stj.ID;
    return true;

  } // AddLooseHits
  
  ////////////////////////////////////////////////
  void FindStartChg(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Finds the charge at the start of a shower and puts it in AveChg of the first
    // point of the shower Tj. This is only done when there is no parent.
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    if(ss.ShowerTjID == 0) return;
    if(ss.ParentID > 0) return;
    auto& stp0 = tjs.allTraj[ss.ShowerTjID - 1].Pts[0];
    
    std::string fcnLabel = inFcnLabel + ".FSC";
    
    stp0.AveChg = 0;
    
    if(ss.AspectRatio > tjs.ShowerTag[10] || ss.DirectionFOM > tjs.ShowerTag[9]) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Not possible due to poor AspectRatio "<<ss.AspectRatio<<" or ss.DirectionFOM "<<ss.DirectionFOM;
      return;
    }

    // Create and fill a vector of the charge at the beginning of the shower in 1 WSE bins
    auto schg = StartChgVec(tjs, cotID, prt);
    if(schg.empty()) return;

    // Look for two consecutive charge entries. Use the second one
    // for the initial guess at the charge
    unsigned short startPt = USHRT_MAX;
    float chg = 0;
    for(unsigned short ii = 0; ii < schg.size() - 1; ++ii) {
      if(schg[ii] > 0 && schg[ii + 1] > 0) {
        startPt = ii + 1;
        chg = schg[ii + 1];
        break;
      }
    }
    if(startPt == USHRT_MAX) return;
    
    // get an initial average and rms using all the points
    float ave = 0;
    float rms = 0;
    float cnt = 0;
    for(unsigned short ii = startPt; ii < schg.size() - 1; ++ii) {
      ave += schg[ii];
      rms += schg[ii] * schg[ii];
      ++cnt;
    } // ii
    ave /= cnt;
    rms = rms - cnt * ave * ave;
    if(rms < 0) return;
    rms = sqrt(rms / (cnt - 1));
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" schg:";
      for(unsigned short ii = 0; ii < 20; ++ii) myprt<<" "<<(int)schg[ii];
      myprt<<"\n First guess at the charge "<<(int)chg<<" average charge of all bins "<<(int)ave<<" rms "<<(int)rms;
    }
    
    // initial guess at the charge rms
    rms = 0.8 * chg;
    
    // Correct for dead wires in this region - maybe later...
//    unsigned short nDeadWires = DeadWireCount();
    
    unsigned short nBinsAverage = 5;
    double maxChg = 2 * chg;
    for(unsigned short nit = 0; nit < 2; ++nit) {
      double sum = 0;
      double sum2 = 0;
      double cnt = 0;
      for(unsigned short ii = startPt; ii < schg.size() - 1; ++ii) {
        // break if we find 2 consecutive high charge points
        if(schg[ii] > maxChg && schg[ii + 1] > maxChg) break;
        // or two zeros
        if(schg[ii] == 0 && schg[ii + 1] == 0) break;
        if(schg[ii] > maxChg) continue;
        sum += schg[ii];
        sum2 += schg[ii] * schg[ii];
        ++cnt;
        if(cnt == nBinsAverage) break;
      } // ii
      // check for a failure
      if(cnt < 3) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" nit "<<nit<<" cnt "<<cnt<<" is too low. sum "<<(int)sum<<" maxChg "<<(int)maxChg;
        // try to recover. Pick the next point
        ++startPt;
        chg = schg[startPt];
        maxChg = 2 * chg;
        continue;
      }
      // convert sum to the average charge
      chg = sum / cnt;
      double arg = sum2 - cnt * chg * chg;
      if(arg < 0) break;
      rms = sqrt(arg / (cnt - 1));
      // don't allow the rms to get crazy
      double maxrms = 0.5 * sum;
      if(rms > maxrms) rms = maxrms;
      maxChg = chg + 2 * rms;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" nit "<<nit<<" cnt "<<cnt<<" chg "<<(int)chg<<" rms "<<(int)rms<<" maxChg "<<(int)maxChg<<" nBinsAverage "<<nBinsAverage;
      nBinsAverage = 20;
    } // nit
    
    stp0.AveChg = chg;
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<cotID<<" Starting charge "<<(int)stp0.AveChg<<" startPt  "<<startPt;
    
  } // FindStartChg
  
  ////////////////////////////////////////////////
  std::vector<float> StartChgVec(TjStuff& tjs, int cotID, bool prt)
  {
    // Returns a histogram vector of the charge in bins of 1 WSE unit at the start of the shower

    ShowerStruct& ss = tjs.cots[cotID - 1];
    constexpr unsigned short nbins = 20;
    std::vector<float> schg(nbins);
    if(ss.ID == 0) return schg; 
    if(ss.TjIDs.empty()) return schg; 
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID-1].Pts[1];
    
    // move the min along point back by 2 WSE so that most of the charge in the hits in the
    // first point is included in the histogram
    float minAlong = ss.ShPts[0].RotPos[0] - 2;
    
    float maxTrans = 4;
    // Tighten up on the maximum allowed transverse position if there is a parent
    if(ss.ParentID > 0) maxTrans = 1;
    float cs = cos(-ss.Angle);
    float sn = sin(-ss.Angle);
    std::array<float, 2> chgPos;
    float along, arg;

    for(auto& sspt : ss.ShPts) {
      unsigned short indx = (unsigned short)((sspt.RotPos[0] - minAlong));
      if(indx > nbins - 1) break;
      // Count the charge if it is within a few WSE transverse from the shower axis
      if(std::abs(sspt.RotPos[1]) > maxTrans) continue;
      unsigned int iht = sspt.HitIndex;
      float& peakTime = tjs.fHits[iht].PeakTime;
      float& amp = tjs.fHits[iht].PeakAmplitude;
      float& rms = tjs.fHits[iht].RMS;
      chgPos[0] = tjs.fHits[iht].ArtPtr->WireID().Wire - stp1.Pos[0];
      for(float time = peakTime - 2.5 * rms; time < peakTime + 2.5 * rms; ++time) {
        chgPos[1] = time * tjs.UnitsPerTick - stp1.Pos[1];
        along = cs * chgPos[0] - sn * chgPos[1];
        if(along < minAlong) continue;
        indx = (unsigned short)(along - minAlong);
        if(indx > nbins - 1) continue;
        arg = (time - peakTime) / rms;
        schg[indx] += amp * exp(-0.5 * arg * arg);
      } // time
    } // sspt

    return schg;
  } // StartChgVec
  
  ////////////////////////////////////////////////
  void DumpShowerPts(TjStuff& tjs, int cotID)
  {
    // Print the shower points to the screen. The user should probably pipe the output to a text file
    // then grep this file for the character string PTS which is piped to a text file which can then be
    // imported into Excel, etc
    // Finds the charge at the start of a shower
    if(cotID > tjs.cots.size()) return;
    
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    std::cout<<"PTS Pos0  Pos1   RPos0 RPos1 Chg TID\n";
    for(auto& pt : ss.ShPts) {
      std::cout<<"PTS "<<std::fixed<<std::setprecision(1)<<pt.Pos[0]<<" "<<pt.Pos[1]<<" "<<pt.RotPos[0]<<" "<<pt.RotPos[1];
      std::cout<<" "<<(int)pt.Chg<<" "<<pt.TID;
      std::cout<<"\n";
    }
    
  } // DumpShowerPts
  
  ////////////////////////////////////////////////
  void CheckQuality(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // drop those that don't meet the requirements. This should be called before 3D matching
    
    std::string fcnLabel = inFcnLabel + ".CQ";
    
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      // enough Tjs?
      unsigned short ntjs = ss.TjIDs.size();
      bool killit = (ntjs < tjs.ShowerTag[7]);
      // Kill runt showers
      if(!killit) killit = (ss.Energy < tjs.ShowerTag[3]);
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  2S"<<ss.ID<<" nTjs "<<ss.TjIDs.size()<<" nTjs "<<ss.TjIDs.size()<<" killit? "<<killit;
      if(!killit) {
        // count the number of Tj points
        unsigned short nTjPts = 0;
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          nTjPts += NumPtsWithCharge(tjs, tj, false);
        }  // tjID
        if(nTjPts < tjs.ShowerTag[6]) killit = true;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"    "<<" nTjPts "<<nTjPts<<" killit? "<<killit;
      } // !killit
      if(killit) MakeShowerObsolete(fcnLabel, tjs, ss.ID, prt);
      
    } // ic
    
    // check for duplicates
    std::vector<int> allIDs;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      allIDs.insert(allIDs.end(), ss.TjIDs.begin(), ss.TjIDs.end());
    } // ss
    if(allIDs.empty()) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  No showers left after quality check";
      return;
    }
    std::sort(allIDs.begin(), allIDs.end());
    for(unsigned short ii = 0; ii < allIDs.size() - 1; ++ii) {
      if(allIDs[ii] == allIDs[ii + 1]) {
        std::cout<<fcnLabel<<" Found duplicate Tjs "<<allIDs[ii]<<"\n";
      }
      auto& tj = tjs.allTraj[allIDs[ii] - 1];
      if(!tj.AlgMod[kInShower] || tj.AlgMod[kKilled]) {
        std::cout<<fcnLabel<<" Tj "<<tj.ID<<" isn't an inShower "<<tj.AlgMod[kInShower]<<" or is Killed "<<tj.AlgMod[kKilled]<<"\n";
      }
    } // ii
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Quality checks complete";
    
  } // CheckQuality

  ////////////////////////////////////////////////
  bool TransferTjHits(TjStuff& tjs, bool prt)
  {
    // Transfer InShower hits in all TPCs and planes to the shower Tjs
    
    bool newShowers = false;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.ShowerTjID == 0) continue;
      // Tp 1 of stj will get all of the shower hits
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      if(!stj.Pts[1].Hits.empty()) {
        std::cout<<"TTjH: ShowerTj T"<<stj.ID<<" already has "<<stj.Pts[1].Hits.size()<<" hits\n";
        continue;
      }
      // Note that UseHit is not used since the size is limited to 16
      for(auto& tjID : ss.TjIDs) {
        unsigned short itj = tjID - 1;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) {
          std::cout<<"TTjH: Coding error. T"<<tjID<<" is a ShowerTj but is in TjIDs\n";
          continue;
        }
        if(!tjs.allTraj[itj].AlgMod[kInShower]) {
          std::cout<<"TTjH: Coding error. Trying to transfer T"<<tjID<<" hits but it isn't an InShower Tj\n";
          continue;
        }
        auto thits = PutTrajHitsInVector(tjs.allTraj[itj], kUsedHits);
        // associate all hits with the charge center TP
        stj.Pts[1].Hits.insert(stj.Pts[1].Hits.end(), thits.begin(), thits.end());
        // kill Tjs that are in showers
        tjs.allTraj[itj].AlgMod[kKilled] = true;
      } //  tjID
      // re-assign the hit -> stj association
      for(auto& iht : stj.Pts[1].Hits) tjs.fHits[iht].InTraj = stj.ID;
      newShowers = true;
     } // ss

    if(prt) mf::LogVerbatim("TC")<<"TTJH: success? "<<newShowers;
    return newShowers;
  } // TransferTjHits

  ////////////////////////////////////////////////
  int GetCotID(TjStuff& tjs, int ShowerTjID)
  {
    for(unsigned short ii = 0; ii < tjs.cots.size(); ++ii) {
      if(ShowerTjID == tjs.cots[ii].ShowerTjID) return ii + 1;
    } // iii
    return 0;
    
  } // GetCotID

  ////////////////////////////////////////////////
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss)
  {
    if(ss.ID == 0) return 0;
    if(ss.TjIDs.empty()) return 0;
    if(ss.ShowerTjID == 0) return 0;
    
    const Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    return ChgToMeV(stj.AveChg);
    
  } // ShowerEnergy
  
  ////////////////////////////////////////////////
  float ChgToMeV(float chg)
  {
    // Conversion from shower charge to energy in MeV. The calibration factor
    // was found by processing 500 pizero events with StudyPiZeros using StudyMode
    return 0.012 * chg;
  }
  
  ////////////////////////////////////////////////
  ShowerStruct3D CreateSS3(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // create a 3D shower and size the vectors that are indexed by plane
    
    ShowerStruct3D ss3;
    ss3.TPCID = tpcid;
    ss3.ID = tjs.showers.size() + 1;
    ss3.Energy.resize(tjs.NumPlanes);
    ss3.EnergyErr.resize(tjs.NumPlanes);
    ss3.MIPEnergy.resize(tjs.NumPlanes);
    ss3.MIPEnergyErr.resize(tjs.NumPlanes);
    ss3.dEdx.resize(tjs.NumPlanes);
    ss3.dEdxErr.resize(tjs.NumPlanes);
    
    return ss3;
    
  } // CreateSS3
  
  ////////////////////////////////////////////////
  unsigned short Create2DShower(TjStuff& tjs, const std::vector<int>& tjl)
  {
    // Create a shower and shower Tj using Tjs in the list
    if(tjl.empty()) return USHRT_MAX;
    
    // Get the CTP from the first Tj
    Trajectory stj;
    stj.CTP = tjs.allTraj[tjl[0]-1].CTP;
    // with three points
    stj.Pts.resize(3);
    for(auto& stp : stj.Pts) {
      stp.CTP = stj.CTP;
      // set all UseHit bits true so we don't get confused
      stp.UseHit.set();
    }
    stj.EndPt[0] = 0;
    stj.EndPt[1] = 2;
    stj.ID = tjs.allTraj.size() + 1;
    // Declare that stj is a shower Tj
    stj.AlgMod[kShowerTj] = true;
    stj.PDGCode = 1111;
    tjs.allTraj.push_back(stj);
    // Create the shower struct
    ShowerStruct ss;
    ss.ID = tjs.cots.size() + 1;
    ss.CTP = stj.CTP;
    // assign all TJ IDs to this ShowerStruct
    ss.TjIDs = tjl;
    // declare them to be InShower
    for(auto tjid : tjl) {
      auto& tj = tjs.allTraj[tjid - 1];
      tj.AlgMod[kInShower] = true;
    } // tjid
    ss.ShowerTjID = stj.ID;
    ss.Envelope.resize(4);
    // put it in TJ stuff. The rest of the info will be added later
    tjs.cots.push_back(ss);
    return ss.ID;
    
  } // Create2DShower
  
  ////////////////////////////////////////////////
  void PrintShowers(std::string fcnLabel, TjStuff& tjs)
  {
    if(tjs.showers.empty()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<fcnLabel<<" 3D showers \n";
    for(auto& ss3 : tjs.showers) {
      myprt<<fcnLabel<<" 3S"<<ss3.ID<<" 3V"<<ss3.Vx3ID<<" ChgPos"<<std::fixed;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) myprt<<" "<<std::setprecision(0)<<ss3.ChgPos[xyz];
      myprt<<" Dir";
      for(unsigned short xyz = 0; xyz < 3; ++xyz) myprt<<" "<<std::setprecision(2)<<ss3.Dir[xyz];
      myprt<<" posInPlane";
      std::vector<float> projInPlane(tjs.NumPlanes);
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(ss3.TPCID.Cryostat, ss3.TPCID.TPC, plane);
        auto tp = MakeBareTP(tjs, ss3.ChgPos, ss3.Dir, inCTP);
        myprt<<" "<<PrintPos(tjs, tp.Pos);
        projInPlane[plane] = tp.Delta;
      } // plane
      myprt<<" projInPlane";
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        myprt<<" "<<std::fixed<<std::setprecision(2)<<projInPlane[plane];
      } // plane
      for(auto cid : ss3.CotIDs) {
        auto& ss = tjs.cots[cid - 1];
        myprt<<"\n  2S"<<ss.ID;
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        myprt<<" T"<<stj.ID;
        myprt<<" "<<PrintPos(tjs, stj.Pts[stj.EndPt[0]].Pos)<<" - "<<PrintPos(tjs, stj.Pts[stj.EndPt[1]].Pos);
      } // ci
      if(ss3.NeedsUpdate) myprt<<" *** Needs update";
      myprt<<"\n";
    } // sss3
  } // PrintShowers
  
  ////////////////////////////////////////////////
  void Print2DShowers(std::string someText, const TjStuff& tjs, CTP_t inCTP, bool printKilledShowers)
  {
    // Prints a one-line summary of 2D showers
    if(tjs.cots.empty()) return;

    mf::LogVerbatim myprt("TC");
    
    // see how many lines were are going to print
    bool printAllCTP = (inCTP == USHRT_MAX);
    if(!printAllCTP) {
      unsigned short nlines = 0;
      for(const auto& ss : tjs.cots) {
        if(!printAllCTP && ss.CTP != inCTP) continue;
        if(!printKilledShowers && ss.ID == 0) continue;
        ++nlines;
      } // ss
      if(nlines == 0) {
        myprt<<someText<<" Print2DShowers: Nothing to print";
        return;
      }
    } // !printAllCTP
    
    // print a header
//    myprt<<someText<<"  ID   CTP  ParID TruParID Energy nTjs  dFOM AspRat  stj  vx0 __Pos0___ nPts dRMS __Pos1___ nPts dRMS __Pos2___ nPts dRMS Angle SS3ID PFPID\n";
    myprt<<someText<<"  ID   CTP  ParID ParFOM TruParID Energy nTjs  dFOM AspRat  stj  vx0 __Pos0___ Chg(k) dRMS __Pos1___ Chg(k) dRMS __Pos2___ Chg(k) dRMS Angle SS3ID PFPID\n";

    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      myprt<<someText<<std::fixed;
      std::string sid = "2S" + std::to_string(ss.ID);
      myprt<<std::setw(4)<<sid;
      myprt<<std::setw(6)<<ss.CTP;
      sid = "NA";
      if(ss.ParentID > 0) sid = "T" + std::to_string(ss.ParentID);
      myprt<<std::setw(7)<<sid;
      myprt<<std::setw(7)<<std::setprecision(2)<<ss.ParentFOM;
      myprt<<std::setw(9)<<ss.TruParentID;
      myprt<<std::setw(7)<<(int)ss.Energy;
      myprt<<std::setw(5)<<ss.TjIDs.size();
      myprt<<std::setw(6)<<std::setprecision(2)<<ss.DirectionFOM;
      myprt<<std::setw(7)<<std::setprecision(2)<<ss.AspectRatio;
      const auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      std::string tid = "T" + std::to_string(stj.ID);
      myprt<<std::setw(5)<<tid;
      std::string vid = "2V" + std::to_string(stj.VtxID[0]);
      myprt<<std::setw(5)<<vid;
      for(auto& spt : stj.Pts) {
        myprt<<std::setw(10)<<PrintPos(tjs, spt.Pos);
        myprt<<std::setw(7)<<std::fixed<<std::setprecision(1)<<spt.Chg/1000;
//        myprt<<std::setw(5)<<spt.NTPsFit;
        myprt<<std::setw(5)<<std::setprecision(1)<<spt.DeltaRMS;
      } // spt
      myprt<<std::setw(6)<<std::setprecision(2)<<stj.Pts[1].Ang;
      std::string sss = "NA";
      if(ss.SS3ID > 0) sss = "3S" + std::to_string(ss.SS3ID);
      myprt<<std::setw(6)<<sss;
      if(ss.SS3ID > 0) {
        auto& ss3 = tjs.showers[ss.SS3ID - 1];
        if(ss3.PFPIndex < tjs.pfps.size()) {
          auto& pfp = tjs.pfps[ss3.PFPIndex];
          std::string pid = "P" + std::to_string(pfp.ID);
          myprt<<std::setw(6)<<pid;
        } else {
          myprt<<std::setw(6)<<"NA";
        }
      } else {
        myprt<<std::setw(6)<<"NA";
      }
      if(ss.NeedsUpdate) myprt<<" *** Needs update";
      myprt<<"\n";
    } // ss
    if(tjs.ShowerTag[12] > 9) {
      // List of Tjs
      for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
        const auto& ss = tjs.cots[ict];
        if(!printAllCTP && ss.CTP != inCTP) continue;
        if(!printKilledShowers && ss.ID == 0) continue;
        myprt<<someText<<std::fixed;
        std::string sid = "2S" + std::to_string(ss.ID);
        myprt<<std::setw(4)<<sid;
        myprt<<" Tjs";
        for(auto id : ss.TjIDs) myprt<<" T"<<id;
        myprt<<"\n";
      } // ict
      
    }
    // Print the envelopes
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      myprt<<someText<<std::fixed;
      std::string sid = "2S" + std::to_string(ss.ID);
      myprt<<std::setw(4)<<sid;
      myprt<<" Envelope";
      for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      myprt<<"\n";
    } // ict
    // List of nearby Tjs
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      myprt<<someText<<std::fixed;
      std::string sid = "2S" + std::to_string(ss.ID);
      myprt<<std::setw(4)<<sid;
      myprt<<" Nearby";
      for(auto id : ss.NearTjIDs) myprt<<" T"<<id;
      myprt<<"\n";
    } // ict
  } // Print2DShowers

} // namespace tca
