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
    // Finish defining the showers, create a companion PFParticle for each one.
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
        if(pfp.TjIDs.empty()) continue;
        unsigned short ndead = 0;
        for(auto tjid : pfp.TjIDs) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kKilled]) ++ndead;
        } // tjid
        if(ndead == 0) continue;
        if(ndead != pfp.TjIDs.size()) {
          std::cout<<"Finish3DShowers: Not all Tjs in P"<<pfp.ID<<" are killed";
          for(auto tid : pfp.TjIDs) {
            auto& tj = tjs.allTraj[tid - 1];
            std::cout<<" T"<<tid<<" dead? "<<tj.AlgMod[kKilled];
          }
          std::cout<<"\n";
        } // ndead
        pfp.ID = 0;
      } // pfp
    } // pfps not empty
    
    // kill orphan 2D vertices
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      auto vxtjs = GetAssns(tjs, "2V", vx2.ID, "T");
      if(vxtjs.empty()) vx2.ID = 0;
    } // vx2
    
    // kill orphan vertices
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      auto vxtjs = GetAssns(tjs, "3V", vx3.ID, "T");
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
    
    bool prt = tjs.ShowerTag[12] >= 0;
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
    if(tpcid != tjs.TPCID && !FillWireHitRange(tjs, tpcid)) return false;

    if(prt) {
      PrintPFPs("FSi", tjs);
      PrintAllTraj("FSi", tjs, debug, USHRT_MAX, 0);
    }
    
    // define a list of tjs that shouldn't be clustered in the same shower because
    // they are matched to pfp particles that are likely parents of 3D showers
    DefineDontCluster(tjs, tpcid, prt);

    // lists of Tj IDs in plane, (list1, list2, list3, ...)
    std::vector<std::vector<std::vector<int>>> bigList(tjs.NumPlanes);
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      std::vector<std::vector<int>> tjList;
      // Make lists of lists of ShowerLike tjs that will become showers
      FindCots(fcnLabel, tjs, inCTP, tjList, prt);
      SaveTjInfo(tjs, tjList, "TJL");
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
        auto ss = CreateSS(tjs, tjl);
        if(ss.ID == 0) continue;
        if(!UpdateShower(fcnLabel, tjs, ss, prt)) continue;
        SaveTjInfo(tjs, ss, "DS");
        FindNearbyTjs(fcnLabel, tjs, ss, prt);
        // don't try to do anything else here until all of the showers are defined
        if(!StoreShower(fcnLabel, tjs, ss)) {
          std::cout<<fcnLabel<<" StoreShower failed 2S"<<ss.ID<<"\n";
          MakeShowerObsolete(fcnLabel, tjs, ss, prt);
        }
      } // tjl
      ChkAssns(fcnLabel, tjs);
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
      MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "MSS");
      // merge sub-showers with other showers
      MergeSubShowers(fcnLabel, tjs, inCTP, prt);
      // merge sub-showers with shower-like tjs
      MergeSubShowersTj(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "MNrby");
      if(prt) Print2DShowers("Nrby", tjs, inCTP, false);
      bool tryMerge = false;
      for(unsigned short ii = 0; ii < tjs.cots.size(); ++ii) {
        auto& ss = tjs.cots[ii];
        if(ss.ID == 0) continue;
        if(ss.CTP != inCTP) continue;
        if(AddTjsInsideEnvelope(fcnLabel, tjs, ss, prt)) tryMerge = true;
        if (tjs.SaveShowerTree) SaveAllCots(tjs, inCTP, "Merge");
      }
      if(tryMerge) MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "ATj2");
      if(prt) Print2DShowers("ATIE", tjs, inCTP, false);
    } // plane
    if(tjs.cots.empty()) return false;
    prt = (dbgPlane > 2);
    if(prt) Print2DShowers("B4", tjs, USHRT_MAX, false);
    // Match in 3D, make 3D showers and define them
    Match2DShowers(fcnLabel, tjs, tpcid, prt);
    SaveAllCots(tjs, "M2DS");
    // Reconcile pfp and shower assns before the Parent search
    Reconcile3D(fcnLabel, tjs, tpcid, false, prt);
    SaveAllCots(tjs, "R3D");
    std::cout<<"Add function to check for neutrino vertex inside showers\n";
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      if(ss3.TPCID != tpcid) continue;
      FindParent(fcnLabel, tjs, ss3, prt);
    } // ss3
    // Reconcile pfp and shower assns again
    Reconcile3D(fcnLabel, tjs, tpcid, true, prt);
    if(prt) Print2DShowers("M2DS", tjs, USHRT_MAX, false);
    SaveAllCots(tjs, "FP");
    
    // kill single Tj showers that aren't matched in 3D
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.TjIDs.size() != 1) continue;
      if(ss.SS3ID > 0) continue;
      MakeShowerObsolete(fcnLabel, tjs, ss, prt);
    } // ss
     
    unsigned short nNewShowers = 0;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      SaveTjInfo(tjs, ss, "Done");
     ++nNewShowers;
    } // ss

    // temp for testing with MC
    if(!tjs.MCPartList.empty()) {
      // get some truth information
      MCParticleListUtils mcpu{tjs};
      for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
        CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
        int truElectronTID = mcpu.PrimaryElectronTjID(inCTP);
        std::cout<<"Plane "<<plane<<" trueElectron T"<<truElectronTID;
        for(auto& ss : tjs.cots) {
          if(ss.ID == 0) continue;
          if(ss.CTP != inCTP) continue;
          if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), truElectronTID) == ss.TjIDs.end()) continue;
          std::cout<<" in 2S"<<ss.ID<<" with Parent T"<<ss.ParentID;
          if(ss.SS3ID > 0) std::cout<<" -> 3S"<<ss.SS3ID;
        } // ss
        std::cout<<"\n";
      } // plane
      std::cout<<"Primary P"<<mcpu.PrimaryElectronPFPID(tpcid)<<"\n";
      Point3_t truStart;
      Vector3_t truDir;
      float truE;
      if(mcpu.PrimaryElectronStart(truStart, truDir, truE)) {
        // assume that the highest energy shower is from the primary electron
        float big = 0;
        int imBig = 0;
        for(auto& ss3 : tjs.showers) {
          float energy = ShowerEnergy(ss3);
          if(energy < big) continue;
          big = energy;
          imBig = ss3.ID;
        } // ss3
        if(imBig > 0) {
          auto& ss3 = tjs.showers[imBig - 1];
          float shMaxSep = PosSep(truStart, ss3.ChgPos);
          float dang = acos(DotProd(truDir, ss3.Dir));
          float energy = ShowerEnergy(ss3);
          std::cout<<"Big 3S"<<ss3.ID<<" E "<<(int)energy<<" truE "<<(int)truE<<" shMax "<<(int)shMaxSep<<" cm.";
          std::cout<<" dang "<<std::fixed<<std::setprecision(2)<<dang<<"\n";
        } // imBig > 0
      } // PrimaryElectronStart success
    } // !tjs.MCPartList.empty()

//    if(prt) Print2DShowers("FSo", tjs, USHRT_MAX, false);
    
    return (nNewShowers > 0);
    
  } // FindShowers3D
  
  ////////////////////////////////////////////////
  bool Reconcile3D(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool parentSearchDone, bool prt)
  {
    // Reconcile pfp and shower assns
    
    if(tjs.showers.size() < 2) return true;
    std::string fcnLabel = inFcnLabel + ".R3D2";
    
    if(prt) Print2DShowers("R3D2i", tjs, USHRT_MAX, false);
    
    for(unsigned short ii = 0; ii < tjs.showers.size() - 1; ++ii) {
      auto iss3 = tjs.showers[ii];
      if(iss3.ID == 0) continue;
      if(iss3.TPCID != tpcid) continue;
      auto iPInSS3 = GetAssns(tjs, "3S", iss3.ID, "P");
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<fcnLabel<<" 3S"<<iss3.ID<<" ->";
        for(auto pid : iPInSS3) myprt<<" P"<<pid;
      } // prt
      for(unsigned short jj = ii + 1; jj < tjs.showers.size(); ++jj) {
        auto jss3 = tjs.showers[jj];
        if(jss3.ID == 0) continue;
        auto jPInSS3 = GetAssns(tjs, "3S", jss3.ID, "P");
        auto shared = SetIntersection(iPInSS3, jPInSS3);
        if(shared.empty()) continue;
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<" Conflict i3S"<<iss3.ID<<" and j3S"<<jss3.ID<<" share";
          for(auto pid : shared) myprt<<" P"<<pid;
        } // prt
        // Compare the InShower likelihoods
        for(auto pid : shared) {
          auto& pfp = tjs.pfps[pid - 1];
          float iProb = InShowerProb(tjs, iss3, pfp);
          float jProb = InShowerProb(tjs, jss3, pfp);
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" i3S"<<iss3.ID<<" prob "<<std::setprecision(3)<<iProb<<"  j3S"<<jss3.ID<<" prob "<<jProb;
          if(iProb > jProb) {
            // remove the remnants of pfp from jss3
            RemovePFP(fcnLabel, tjs, pfp, jss3, true, prt);
            // and add them to iss3
            AddPFP(fcnLabel, tjs, pfp.ID, iss3, true, prt);
          } else {
            RemovePFP(fcnLabel, tjs, pfp, iss3, true, prt);
            AddPFP(fcnLabel, tjs, pfp.ID, jss3, true, prt);
          }
        } // pid
      } // jj
    } // ii
        
    // Look for an in-shower pfp that is not the shower parent that is attached to a vertex. 
    // Remove the attachment and any parent - daughter assn
    if(parentSearchDone) {
      for(auto& ss3 : tjs.showers) {
        if(ss3.ID == 0) continue;
        if(ss3.TPCID != tpcid) continue;
        auto PIn3S = GetAssns(tjs, "3S", ss3.ID, "P");
        for(auto pid : PIn3S) {
          if(pid == ss3.ParentID) continue;
          auto& pfp = tjs.pfps[pid - 1];
          for(unsigned short end = 0; end < 2; ++end) {
            if(pfp.Vx3ID[end] <= 0) continue;
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<fcnLabel<<" Detach 3S"<<ss3.ID<<" -> P"<<pfp.ID<<"_"<<end<<" -> 3V"<<pfp.Vx3ID[end];
              if(pfp.ParentID > 0) myprt<<" ->Parent P"<<pfp.ParentID;
            }
            // remove P -> P parent-daughter assn
            pfp.Vx3ID[end] = 0;
            if(pfp.ParentID > 0) {
              auto& parentPFP = tjs.pfps[pfp.ParentID - 1];
              std::vector<int> newDtrIDs;
              for(auto did : parentPFP.DtrIDs) if(did != pfp.ID) newDtrIDs.push_back(did);
              parentPFP.DtrIDs = newDtrIDs;
            } // pfp Parent exists
          } // end
        } // pid
      } // ss3
    } // parentSearchDone
    
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    // now look for 2D showers that not matched in 3D and have tjs
    // that are 3D-matched
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.SS3ID > 0) continue;
      if(DecodeCTP(ss.CTP).TPC != tpc) continue;
      if(DecodeCTP(ss.CTP).Cryostat != cstat) continue;
      std::vector<int> matchedTjs;
      for(auto tid : ss.TjIDs) if(tjs.allTraj[tid - 1].AlgMod[kMat3D]) matchedTjs.push_back(tid);
      if(matchedTjs.empty()) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<fcnLabel<<" 2S"<<ss.ID<<" is not 3D-matched but has 3D-matched Tjs:";
        for(auto tid : matchedTjs) {
          myprt<<" T"<<tid;
          auto TInP = GetAssns(tjs, "T", tid, "P");
          if(!TInP.empty()) myprt<<" -> P"<<TInP[0];
        } // tid
      } // prt
    } // ss
    
    
    if(prt) Print2DShowers("R3D2o", tjs, USHRT_MAX, false);
    
    ChkAssns(fcnLabel, tjs);
    
    return true;
  } // Reconcile3D
  
  ////////////////////////////////////////////////
  bool Reconcile3D(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // checks consistency between pfparticles, showers and tjs associated with ss3
    if(ss3.ID == 0) return false;
    // it isn't a failure if there is a 3D shower in two planes
    if(ss3.CotIDs.size() < 3) return true;
    std::string fcnLabel = inFcnLabel + ".R3D";

    if(prt) Print2DShowers("R3Di", tjs, USHRT_MAX, false);

    // make local copies so we can recover from a failure
    auto oldSS3 = ss3;
    std::vector<ShowerStruct> oldSS(ss3.CotIDs.size());
    for(unsigned short ii = 0; ii < ss3.CotIDs.size(); ++ii) {
      oldSS[ii] = tjs.cots[ss3.CotIDs[ii] - 1];
    }
    
    std::vector<std::vector<int>> plist(ss3.CotIDs.size());
    for(unsigned short ci = 0; ci < ss3.CotIDs.size(); ++ci) {
      auto& ss = tjs.cots[ss3.CotIDs[ci] - 1];
      for(auto tid : ss.TjIDs) {
        auto tToP = GetAssns(tjs, "T", tid, "P");
        if(tToP.empty()) continue;
        // there should only be one pfp for a tj
        int pid = tToP[0];
        if(std::find(plist[ci].begin(), plist[ci].end(), pid) == plist[ci].end()) plist[ci].push_back(pid);
      } // tid
    } // ci
    // count the occurrence of each pfp
    std::vector<std::array<int, 2>> p_cnt;
    for(auto& pl : plist) {
      for(auto pid : pl) {
        unsigned short indx = 0;
        for(indx = 0; indx < p_cnt.size(); ++indx) if(p_cnt[indx][0] == pid) break;
        if(indx == p_cnt.size()) {
          // not found so add it
          p_cnt.push_back(std::array<int,2> {{pid, 1}});
        } else {
          ++p_cnt[indx][1];
        }
      } // pid
    } // pl
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 3S"<<ss3.ID<<"\n";
      for(unsigned short ci = 0; ci < ss3.CotIDs.size(); ++ci) {
        myprt<<" -> 2S"<<ss3.CotIDs[ci]<<" ->";
        for(auto pid : plist[ci]) myprt<<" P"<<pid;
        myprt<<"\n";
      } // ci
      myprt<<" P<ID>_count:";
      for(auto& pc : p_cnt) myprt<<" P"<<pc[0]<<"_"<<pc[1];
    } // prt
    
    for(auto& pc : p_cnt) {
      // matched in all planes?
      if(pc[1] == (int)ss3.CotIDs.size()) continue;
      if(pc[1] == 2) {
        // missing a tj in a plane or is this a two-plane pfp?
        auto& pfp = tjs.pfps[pc[0] - 1];
        if(pfp.TjIDs.size() > 2) {
          // ensure that none of the tjs in this pfp are included in a different shower
          auto PIn2S = GetAssns(tjs, "P", pfp.ID, "2S");
          auto sDiff = SetDifference(PIn2S, ss3.CotIDs);
          // Not sure if this can happen
//          if(sDiff.size() > 1) {std::cout<<fcnLabel<<" sDiff size > 2. Is this possible?\n";}
          if(!sDiff.empty() && std::find(ss3.CotIDs.begin(), ss3.CotIDs.end(), sDiff[0]) == ss3.CotIDs.end()) continue;
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<fcnLabel<<" 3S"<<ss3.ID<<" P"<<pfp.ID<<" ->";
            for(auto sid : PIn2S) myprt<<" 2S"<<sid;
            myprt<<" sDiff";
            for(auto sid : sDiff) myprt<<" 2S"<<sid;
          } // prt
          // missed a tj in a 2D shower so "add the PFP to the shower" and update it
          if(AddPFP(fcnLabel, tjs, pfp.ID, ss3, true, prt)) {
            // Update the local copies
            oldSS3 = ss3;
            if(ss3.CotIDs.size() != oldSS.size()) {
              std::cout<<fcnLabel<<" Major failure...";
              return false;
            }
            for(unsigned short ii = 0; ii < ss3.CotIDs.size(); ++ii) oldSS[ii] = tjs.cots[ss3.CotIDs[ii] - 1];
          } else {
            // restore the previous state
            ss3 = oldSS3;
            for(unsigned short ii = 0; ii < oldSS.size(); ++ii) {
              auto& ss = oldSS[ii];
              tjs.cots[ss.ID - 1] = ss;
            } // ii
          } // AddPFP failed
        } // pfp.TjIDs.size() > 2
      } else {
        // only one occurrence. check proximity to ss3
        auto& pfp = tjs.pfps[pc[0] - 1];
        unsigned short nearEnd = 1 - FarEnd(tjs, pfp, ss3.ChgPos);
        float prob = InShowerProb(tjs, ss3, pfp);
        float sep = PosSep(pfp.XYZ[nearEnd], ss3.ChgPos);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<" one occurrence: P"<<pfp.ID<<"_"<<nearEnd<<" closest to ChgPos";
          myprt<<" ChgPos "<<std::fixed<<std::setprecision(1)<<ss3.ChgPos[0]<<" "<<ss3.ChgPos[1]<<" "<<ss3.ChgPos[2];
          myprt<<" sep "<<sep;
          myprt<<" InShowerProb "<<prob;
        } // prt 
        if(sep > 30 && !RemovePFP(fcnLabel, tjs, pfp, ss3, true, prt)) {
          std::cout<<"RemovePFP failed \n";
        }
      } // only one occurrence.
    } // pc
    
    return UpdateShower(fcnLabel, tjs, ss3, prt);
    
    ChkAssns(fcnLabel, tjs);
    
    if(prt) Print2DShowers("R3Do", tjs, USHRT_MAX, false);
    
  } // Reconcile3D
/*
  ////////////////////////////////////////////////
  void FindInShowerPFPs(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, std::vector<std::vector<int>>& plists)
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
    std::vector<int> pfpList;
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
      std::vector<int> plist;
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
          int pid2 = cp.id1;
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
            int oid = pcp.id1;
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
    std::vector<int> flat;
    for(auto& plist : plists) flat.insert(flat.end(), plist.begin(), plist.end());
    auto notClustered = SetDifference(pfpList, flat);
    for(auto nc : notClustered) {
      auto& pfp = tjs.pfps[nc - 1];
      if(!IsInShower(tjs, pfp.TjIDs)) continue;
      std::vector<int> plist(1, nc);
      plists.push_back(plist);
    }

    mf::LogVerbatim myprt("TC");
    myprt<<fcnLabel<<" plists size "<<plists.size()<<"\n";
    for(unsigned short ip = 0; ip < plists.size(); ++ip) {
      auto& plist = plists[ip];
      myprt<<"ip "<<ip;
      for(auto pid : plist) myprt<<" P"<<pid;
      myprt<<"\n";
    } // plist

  } // FindInShowerPFPs
*/
  ////////////////////////////////////////////////
  void KillVerticesInShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    // make the vertices inside the shower envelope obsolete and update dontCluster
    if(ss.ID == 0) return;
    if(!tjs.UseAlg[kKillInShowerVx]) return;
    std::string fcnLabel = inFcnLabel + ".KVIS";
    
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      if(vx2.CTP != ss.CTP) continue;
      // ensure it isn't associated with a neutrino vertex
      if(vx2.Vx3ID > 0 && tjs.vtx3[vx2.Vx3ID - 1].Neutrino) continue;
      if(!PointInsideEnvelope(vx2.Pos, ss.Envelope)) continue;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Clobber 2V"<<vx2.ID<<" -> 3V"<<vx2.Vx3ID<<" inside 2S"<<ss.ID;
      // update dontCluster
      for(auto& dc : tjs.dontCluster) {
        if(dc.TjIDs[0] == 0) continue;
        if(dc.Vx2ID != vx2.ID) continue;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Remove T"<<dc.TjIDs[0]<<"-T"<<dc.TjIDs[0]<<" in dontCluster";
        dc.TjIDs[0] = 0;
        dc.TjIDs[1] = 0;
      } // dc
      if(vx2.Vx3ID > 0) {
        auto TIn3V = GetAssns(tjs, "3V", vx2.Vx3ID, "T");
        for(auto tid : TIn3V) tjs.allTraj[tid - 1].AlgMod[kKillInShowerVx] = true;
        auto& vx3 = tjs.vtx3[vx2.Vx3ID - 1];
        MakeVertexObsolete(tjs, vx3);
      } else {
        auto TIn2V = GetAssns(tjs, "2V", vx2.ID, "T");
        for(auto tid : TIn2V) tjs.allTraj[tid - 1].AlgMod[kKillInShowerVx] = true;
        MakeVertexObsolete(tjs, vx2, true);
      }
    } // vx2
    
  } // KillVerticesInShower

  ////////////////////////////////////////////////
  void CompleteIncompleteShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // Find low-energy two-plane showers and try to complete it by making a 2D shower in the third
    // plane using 3D matched tjs
    
    if(tjs.NumPlanes != 3) return;
    if(ss3.CotIDs.size() != 2) return;
    
    if(!tjs.UseAlg[kCompleteShower]) return;
    
    std::string fcnLabel = inFcnLabel + ".CIS";
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID;
    
    auto& iss = tjs.cots[ss3.CotIDs[0] - 1];
    auto& jss = tjs.cots[ss3.CotIDs[1] - 1];
    // make a list of pfps for each SS
    std::vector<int> iplist;
    for(auto tid : iss.TjIDs) {
      auto plist = GetAssns(tjs, "T", tid, "P");
      if(!plist.empty()) iplist.insert(iplist.end(), plist.begin(), plist.end());
    } // tid
    std::vector<int> jplist;
    for(auto tid : jss.TjIDs) {
      auto plist = GetAssns(tjs, "T", tid, "P");
      if(!plist.empty()) jplist.insert(jplist.end(), plist.begin(), plist.end());
    } // tid
    // look for pfps that have tjs in both showers
    auto shared = SetIntersection(iplist, jplist);
    if(shared.empty()) return;
    // put the list of tjs for both SS into a flat vector to simplify searching
    std::vector<int> flat = iss.TjIDs;
    flat.insert(flat.end(), jss.TjIDs.begin(), jss.TjIDs.end());
    // make a list of tjs in the k plane that maybe should made into a shower if they
    // aren't already in a shower that failed the 3D match
    std::vector<int> ktlist;
    for(auto pid : shared) {
      auto& pfp = tjs.pfps[pid - 1];
      for(auto tid : pfp.TjIDs) {
        // ignore the tjs that are already in the shower in the other planes
        if(std::find(flat.begin(), flat.end(), tid) != flat.end()) continue;
        if(std::find(ktlist.begin(), ktlist.end(), tid) == ktlist.end()) ktlist.push_back(tid);
        // look for 2D vertices attached to this tj and add all attached tjs to ktlist
        auto& tj = tjs.allTraj[tid - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] <= 0) continue;
          auto& vx2 = tjs.vtx[tj.VtxID[end] - 1];
          auto TIn2V = GetAssns(tjs, "2V", vx2.ID, "T");
          for(auto vtid : TIn2V) {
            if(std::find(ktlist.begin(), ktlist.end(), vtid) == ktlist.end()) ktlist.push_back(vtid);
          }
        } // end
      } // tid
    } // pid
    if(ktlist.empty()) return;
    // list of 2D showers that include tjs in ktlist
    std::vector<int> ksslist;
    for(auto tid : ktlist) {
      auto& tj = tjs.allTraj[tid - 1];
      if(tj.SSID == 0) continue;
      // ignore showers that are 3D-matched. This case should be handled elsewhere by a merging function
      auto& ss = tjs.cots[tj.SSID - 1];
      if(ss.SS3ID > 0) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Found existing T"<<tid<<" -> 2S"<<ss.ID<<" -> 3S"<<ss.SS3ID<<" assn. Give up";
        return;
      }
      if(std::find(ksslist.begin(), ksslist.end(), ss.ID) == ksslist.end()) ksslist.push_back(ss.ID);
    } // tid
    // find the shower energy for this list
    float ktlistEnergy = ShowerEnergy(tjs, ktlist);
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 3S"<<ss3.ID<<"\n";
      myprt<<" -> i2S"<<iss.ID<<" ->";
      for(auto pid : iplist) myprt<<" P"<<pid;
      myprt<<"\n";
      myprt<<" -> j2S"<<jss.ID<<" ->";
      for(auto pid : jplist) myprt<<" P"<<pid;
      myprt<<"\n";
      geo::PlaneID iPlaneID = DecodeCTP(iss.CTP);
      geo::PlaneID jPlaneID = DecodeCTP(jss.CTP);
      unsigned short kplane = 3 - iPlaneID.Plane - jPlaneID.Plane;
      myprt<<" kplane "<<kplane<<" ktlist:";
      for(auto tid : ktlist) myprt<<" T"<<tid;
      myprt<<" ktlistEnergy "<<ktlistEnergy;
      if(ksslist.empty()) {
        myprt<<"\n No matching showers in kplane";
      }  else {
        myprt<<"\n";
        myprt<<" Candidate showers:";
        for(auto ssid : ksslist) {
          myprt<<" 2S"<<ssid;
          auto& sst = tjs.cots[ssid - 1];
          if(sst.SS3ID > 0) myprt<<"_3S"<<sst.SS3ID;
        } // ssid
      } // ssList not empty
    } // prt
    if(ksslist.size() > 1) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Found more than 1 shower. Need some better code here";
      return;
    }
    if(ktlistEnergy > 2 * ShowerEnergy(ss3)) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" ktlistEnergy exceeds 2 * ss3 energy. Need some better code here";
      return;
    } // ktlistEnergy too high
    
    if(ksslist.empty()) {
      // no 2D shower so make one using ktlist
      auto kss = CreateSS(tjs, ktlist);
      if(kss.ID == 0) return;
      kss.SS3ID = ss3.ID;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" create new 2S"<<kss.ID<<" from ktlist";
      if(!UpdateShower(fcnLabel, tjs, kss, prt)) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" UpdateShower failed 2S"<<kss.ID;
        MakeShowerObsolete(fcnLabel, tjs, kss, prt);
        return;
      } // UpdateShower failed
      if(!StoreShower(fcnLabel, tjs, kss)) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" StoreShower failed";
        MakeShowerObsolete(fcnLabel, tjs, kss, prt);
        return;
      } // StoreShower failed
      ss3.CotIDs.push_back(kss.ID);
      auto& stj = tjs.allTraj[kss.ShowerTjID - 1];
      stj.AlgMod[kCompleteShower] = true;
      ss3.NeedsUpdate = true;
      return;
    } // ksslist empty
    
    // associate ksslist[0] with 3S
    auto& ss = tjs.cots[ksslist[0] - 1];
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" found pfp-matched 2S"<<ss.ID;
    ss.SS3ID = ss3.ID;
    ss3.CotIDs.push_back(ss.ID);
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    stj.AlgMod[kCompleteShower] = true;
    ss3.NeedsUpdate = true;
    
    ChkAssns(fcnLabel, tjs);

 } // CompleteIncompleteShower
  
  ////////////////////////////////////////////////
  void Match2DShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // Match 2D showers using position and direction to create 3D showers
    
    std::string fcnLabel = inFcnLabel + ".M2DS";
    if(prt) mf::LogVerbatim("TC")<<fcnLabel;
    
    float fomCut = 2;
    
    ChkAssns(fcnLabel, tjs);
    
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
        if(iss.SS3ID > 0) break;
        unsigned short jIndx = sortVec[jj].index;
        ShowerStruct& jss = tjs.cots[jIndx];
        // already matched?
        if(iss.SS3ID > 0) break;
        if(jss.SS3ID > 0) continue;
        if(jss.CTP == iss.CTP) continue;
        Trajectory& jstj = tjs.allTraj[jss.ShowerTjID - 1];
        TrajPoint3 tp3;
        if(!MakeTp3(tjs, istj.Pts[1], jstj.Pts[1], tp3, true)) continue;
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
          ss3.MatchFOM = fomij;
          ss3.PFPIndex = USHRT_MAX;
          if(!StoreShower(fcnLabel, tjs, ss3)) continue;
          if(prt) mf::LogVerbatim("TC")<<" new 2-plane TPC 3S"<<ss3.ID<<" with fomij "<<fomij;
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
        ss3.MatchFOM = bestFOM;
        if(bestck == USHRT_MAX) {
          // showers match in 2 planes
          ss3.CotIDs.resize(2);
          ss3.CotIDs[0] = iss.ID;
          ss3.CotIDs[1] = jss.ID;
          ss3.Energy[iplaneID.Plane] = iss.Energy;
          ss3.Energy[jplaneID.Plane] = jss.Energy;
          if(prt) mf::LogVerbatim("TC")<<" new 2-plane 3S"<<ss3.ID<<" using 2S"<<iss.ID<<" 2S"<<jss.ID<<" with FOM "<<ss3.MatchFOM<<" try to complete it";
          CompleteIncompleteShower(fcnLabel, tjs, ss3, prt);
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
          if(prt) mf::LogVerbatim("TC")<<" new 3-plane 3S"<<ss3.ID<<" using 2S"<<iss.ID<<" 2S"<<jss.ID<<" 2S"<<tjs.cots[ck].ID<<" with FOM "<<ss3.MatchFOM;
        }
        ss3.MatchFOM = 0.5 * (fomij + bestFOM);
        // sort the IDs
        std::sort(ss3.CotIDs.begin(), ss3.CotIDs.end());
        // Set the 3S -> 2S assns and store it
        if(!StoreShower(fcnLabel, tjs, ss3)) {
          MakeShowerObsolete(fcnLabel, tjs, ss3, prt);
          continue;
        }
        // make a reference to the stored shower
        auto& nss3 = tjs.showers[tjs.showers.size() - 1];
        if(nss3.NeedsUpdate) UpdateShower(fcnLabel, tjs, nss3, prt);
        // reconcile tj -> 2S -> 3S and tj -> pfps
        if(!Reconcile3D(fcnLabel, tjs, nss3, prt)) {
          MakeShowerObsolete(fcnLabel, tjs, nss3, prt);
          continue;
        }
        if(nss3.NeedsUpdate) UpdateShower(fcnLabel, tjs, nss3, prt);
        if(prt) mf::LogVerbatim("TC")<<" 3S"<<nss3.ID<<" updated";
        break;
      } // cj
    } // ci
    
    ChkAssns(fcnLabel, tjs);
    
    if(prt) PrintShowers("M2DS", tjs);

  } // Match2DShowers

  ////////////////////////////////////////////////
  bool UpdateShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    // This is intended to be a single function replacement for FCC, FA, USWP, etc. The calling
    // function should have set NeedsUpdate true. A complete re-build is done if the ShPts vector
    // is empty. This is only required if a tj is removed from the shower. When adding a tj
    // to the shower the AddTj function appends the tj points to ShPts but doesn't fill
    // the ShPts RotPos values.
    // This function doesn't alter or check associations btw showers and tjs.
    
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    if(ss.ShowerTjID <= 0 || ss.ShowerTjID > (int)tjs.allTraj.size()) return false;
    if(ss.ParentID > 0 && ss.ParentID > (int)tjs.allTraj.size()) return false;
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return false;

    std::string fcnLabel = inFcnLabel + ".U2S";

    if(!ss.NeedsUpdate && !ss.ShPts.empty()) {
//      std::cout<<fcnLabel<<" 2S"<<ss.ID<<" doesn't need an update\n";
      return true;
    }

    // initialize the variables that will be defined in this function
    ss.Energy = 0; // This is just ShowerEnergy(stj.TotChg) and could be deleted
    ss.AspectRatio = 10;
    // Direction FOM (0 = good). This is a property of the shower shape and is not
    // defined by the presence or absence of a parent tj start point
    ss.DirectionFOM = 10;
    // Total charge of all hits in the shower
    stj.TotChg = 0;
    for(auto& stp : stj.Pts) {
      // Shower start, charge center, and shower end
      stp.Pos = {{0.0, 0.0}};
      // Charge weighted average of hits this section (TP) along the shower
      stp.HitPos = {{0.0, 0.0}};
      // Direction from the start to the charge center - same for all TPs
      stp.Dir = {{0.0, 0.0}};
      // Hit charge in each section
      stp.Chg = 0;
      // transverse rms of hit positions relative to HitPos in this section
      stp.DeltaRMS = 0;
      // number of hits in this section
      stp.NTPsFit = 0;
    } // stp
    
    ss.ShPts.clear();
    for(auto tjid : ss.TjIDs) {
      if(tjid <= 0 || tjid > (int)tjs.allTraj.size()) return false;
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.CTP != ss.CTP) return false;
      if(tj.AlgMod[kShowerTj]) return false;
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
        } // ii
      } // ipt
    } // tjid
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" nShPts "<<ss.ShPts.size();
    
    if(ss.ShPts.size() < 3) return false;
    
    // find the charge center and total charge
    auto& stp1 = stj.Pts[1];
    for(auto& shpt : ss.ShPts) {
      stp1.Pos[0] += shpt.Chg * shpt.Pos[0];
      stp1.Pos[1] += shpt.Chg * shpt.Pos[1];
      stj.TotChg += shpt.Chg;
    } // shpt
    if(stj.TotChg <= 0) return false;
    stp1.Pos[0] /= stj.TotChg;
    stp1.Pos[1] /= stj.TotChg;
    ss.Energy = ChgToMeV(stj.TotChg);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" Chg ctr "<<PrintPos(tjs, stp1.Pos)<<" Energy "<<(int)ss.Energy<<" MeV";
    
    // find the direction using the shower parent if one exists
    if(ss.ParentID > 0) {
      // Set the direction to be the start of the parent to the shower center
      auto& ptj = tjs.allTraj[ss.ParentID - 1];
      // find the parent end farthest away from the charge center
      unsigned short pend = FarEnd(tjs, ptj, stp1.Pos);
      auto& ptp = ptj.Pts[ptj.EndPt[pend]];
      stp1.Dir = PointDirection(ptp.Pos, stp1.Pos);
      stp1.Ang = atan2(stp1.Dir[1], stp1.Dir[0]);
    } else {
      // find the shower direction using the points
      double sum = 0.;
      double sumx = 0.;
      double sumy = 0.;
      double sumxy = 0.;
      double sumx2 = 0.;
      double sumy2 = 0.;
      for(auto& shpt : ss.ShPts) {
        sum += shpt.Chg;
        double xx = shpt.Pos[0] - stp1.Pos[0];
        double yy = shpt.Pos[1] - stp1.Pos[1];
        sumx += shpt.Chg * xx;
        sumy += shpt.Chg * yy;
        sumxy += shpt.Chg * xx * yy;
        sumx2 += shpt.Chg * xx * xx;
        sumy2 += shpt.Chg * yy * yy;
      } // shpt
      double delta = sum * sumx2 - sumx * sumx;
      if(delta == 0) return false;
      // A is the intercept (This should be ~0 )
//      double A = (sumx2 * sumy - sumx * sumxy) / delta;
      // B is the slope
      double B = (sumxy * sum  - sumx * sumy) / delta;
      stp1.Ang = atan(B);
      stp1.Dir[0] = cos(stp1.Ang);
      stp1.Dir[1] = sin(stp1.Ang);
    } // no shower parent
    
    // TODO: ss.Angle should be eliminated. The shower tj Ang should be used instead
    ss.Angle = stp1.Ang;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" dir "<<std::fixed<<std::setprecision(2)<<stp1.Dir[0]<<" "<<stp1.Dir[1]<<" Angle "<<stp1.Ang;
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      if(ipt == 1) continue;
      stj.Pts[ipt].Dir = stp1.Dir;
      stj.Pts[ipt].Ang = stp1.Ang;
    } // ipt
    
    // fill the RotPos vector and sort
    std::vector<SortEntry> sortVec(ss.ShPts.size());
    unsigned short indx = 0;
    double cs = cos(-stp1.Ang);
    double sn = sin(-stp1.Ang);
    for(auto& shpt : ss.ShPts) {
      double xx = shpt.Pos[0] - stp1.Pos[0];
      double yy = shpt.Pos[1] - stp1.Pos[1];
      shpt.RotPos[0] = cs * xx - sn * yy;
      shpt.RotPos[1] = sn * xx + cs * yy;
      sortVec[indx].index = indx;
      sortVec[indx].length = shpt.RotPos[0];
      ++indx;
    } // shpt
    std::sort(sortVec.begin(), sortVec.end(), lessThan);
    // put the points vector into the sorted order
    auto tPts = ss.ShPts;
    for(unsigned short ii = 0; ii < ss.ShPts.size(); ++ii) ss.ShPts[ii] = tPts[sortVec[ii].index];
    
    // Calculate the aspect ratio
    Point2_t alongTrans {{0.0, 0.0}};
    for(auto& shpt : ss.ShPts) {
      alongTrans[0] += shpt.Chg * std::abs(shpt.RotPos[0]);
      alongTrans[1] += shpt.Chg * std::abs(shpt.RotPos[1]);
    } // shpt
    alongTrans[0] /= stj.TotChg;
    alongTrans[1] /= stj.TotChg;
    if(alongTrans[1] == 0) return false;
    ss.AspectRatio = alongTrans[1] / alongTrans[0];
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" AspectRatio "<<ss.AspectRatio;
    
    // analyze the charge in three sections. Fill the stj HitPos and find DeltaRMS
    if(!AnalyzeRotPos(fcnLabel, tjs, ss, prt)) return false;
    
    // Reverse the shower direction if needed and define the start point
    if(ss.ParentID > 0) {
      // The direction was defined by the start of a parent to the charge center. Check the consistency
      // with ShPts and reverse if needed
      auto& ptj = tjs.allTraj[ss.ParentID - 1];
      // find the parent end farthest away from the charge center
      unsigned short pend = FarEnd(tjs, ptj, stp1.Pos);
      auto& ptp = ptj.Pts[ptj.EndPt[pend]];
      auto& firstShPt = ss.ShPts[0];
      auto& lastShPt = ss.ShPts[ss.ShPts.size() - 1];
      if(PosSep2(ptp.Pos, lastShPt.Pos) < PosSep2(ptp.Pos, firstShPt.Pos)) ReverseShower(fcnLabel, tjs, ss, prt);
      stj.Pts[0].Pos = ptp.Pos;
    } else {
      // no parent exists. Compare the DeltaRMS at the ends
      if(stj.Pts[2].DeltaRMS < stj.Pts[0].DeltaRMS) ReverseShower(fcnLabel, tjs, ss, prt);
      stj.Pts[0].Pos = ss.ShPts[0].Pos;
    } // no parent
    
    if(stj.Pts[2].DeltaRMS > 0) ss.DirectionFOM = stj.Pts[0].DeltaRMS / stj.Pts[2].DeltaRMS;
    // define the end point
    stj.Pts[2].Pos = ss.ShPts[ss.ShPts.size() - 1].Pos;
    
    DefineEnvelope(fcnLabel, tjs, ss, prt);
    ss.NeedsUpdate = false;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" updated";
    return true;
    
  } // UpdateShower
  
  ////////////////////////////////////////////////
  bool UpdateShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // Updates the 3D shower presumably because the 2D showers were changed or need to be updated. 
    // This function returns false if there was a failure.
    
    if(ss3.ID == 0) return false;
    if(ss3.CotIDs.size() < 2) return false;
    
    std::string fcnLabel = inFcnLabel + ".U3S";
    
    // see if any of the 2D showers need an update
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      if(ss.NeedsUpdate && prt) std::cout<<fcnLabel<<" ********* 3S"<<ss3.ID<<" 2S"<<ss.ID<<" needs an update...\n";
      UpdateShower(fcnLabel, tjs, ss, prt);
    } // ci
    
    // check consistency 
    if(ss3.ParentID > 0) {
      auto& pfp = tjs.pfps[ss3.ParentID - 1];
      unsigned short pend = FarEnd(tjs, pfp, ss3.ChgPos);
      if(pfp.Vx3ID[pend] != ss3.Vx3ID) {
        if(prt) std::cout<<fcnLabel<<" ********* 3S"<<ss3.ID<<" has parent P"<<ss3.ParentID<<" with a vertex that is not attached the shower\n";
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
          if(!MakeTp3(tjs, istj.Pts[ipt], jstj.Pts[ipt], tp3, true)) continue;
          double ptchg = 0.5 * (istj.Pts[ipt].Chg + jstj.Pts[ipt].Chg);
          chg[ipt] += ptchg;
          for(unsigned short xyz = 0; xyz < 3; ++xyz) {
            pos[ipt][xyz] += ptchg * tp3.Pos[xyz];
            dir[xyz] += ptchg * tp3.Dir[xyz];
          } // xyz
        } // ipt
      } // jj
    } // ii
    
    unsigned short nok = 0;
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      if(chg[ipt] == 0) continue;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pos[ipt][xyz] /= chg[ipt];
      SetMag(dir, 1);
      ++nok;
    } // ipt
    
    if(nok != 3) {
//      if(prt) std::cout<<fcnLabel<<" ********* 3S"<<ss3.ID<<" Can't find 3 points that match in 3D. \n";
      return false;
    }
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
      auto& ss = tjs.cots[cid - 1];
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
    return true;

  } // UpdateShower

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
    if(icid == 0 || icid > (int)tjs.cots.size()) return 100;
    if(jcid == 0 || jcid > (int)tjs.cots.size()) return 100;
    if(kcid == 0 || kcid > (int)tjs.cots.size()) return 100;
    
    float ijfom = Match3DFOM(inFcnLabel, tjs, icid, jcid, prt);
    float jkfom = Match3DFOM(inFcnLabel, tjs, jcid, kcid, prt);
    
    return 0.5 * (ijfom + jkfom);
    
  } // Match3DFOM
  
  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, int icid, int jcid, bool prt)
  {
    // returns a Figure of Merit for a 3D match of two showers
    if(icid == 0 || icid > (int)tjs.cots.size()) return 100;
    if(jcid == 0 || jcid > (int)tjs.cots.size()) return 100;
    
    auto& iss = tjs.cots[icid - 1];
    auto& istj = tjs.allTraj[iss.ShowerTjID - 1];    
    auto& jss = tjs.cots[jcid - 1];
    auto& jstj = tjs.allTraj[jss.ShowerTjID - 1];
    
    if(iss.CTP == jss.CTP) return 100;
    
    std::string fcnLabel = inFcnLabel + ".MFOM";
    
    float energyAsym = std::abs(iss.Energy - jss.Energy) / (iss.Energy + jss.Energy);
    
    // don't apply the asymmetry cut on low energy showers
    if((iss.Energy > 200 || jss.Energy > 200) && energyAsym > 0.5) return 50;
    
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
  bool RemovePFP(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, ShowerStruct3D& ss3, bool doUpdate, bool prt)
  {
    // removes the tjs in the pfp from the ss3 2D showers and optionally update. This function only returns
    // false if there was a failure. The absence of any pfp Tjs in ss3 is not considered a failure

    if(pfp.ID == 0 || ss3.ID == 0) return false;

    std::string fcnLabel = inFcnLabel + ".RemP";
    for(auto tid : pfp.TjIDs) {
      for(auto cid : ss3.CotIDs) {
        auto& ss = tjs.cots[cid - 1];
        if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tid) == ss.TjIDs.end()) continue;
        if(!RemoveTj(fcnLabel, tjs, tid, ss, doUpdate, prt)) return false;
        ss3.NeedsUpdate = true;
      } // cid
    } // ptid
    
    if(doUpdate && ss3.NeedsUpdate) UpdateShower(fcnLabel, tjs, ss3, prt);
    return true;
    
  } // Remove PFP
  
  ////////////////////////////////////////////////
  bool AddPFP(std::string inFcnLabel, TjStuff& tjs, int pID, ShowerStruct3D& ss3, bool doUpdate, bool prt)
  {
    // Add the tjs in the pfp with id = pID to the 2D showers in ss3 and optionally update everything. This
    // function returns true if the addition was successful or if the Tjs in the pfp are already in ss3.
    // This function returns false if there was a failure. There isn't any error recovery.
    
    std::string fcnLabel = inFcnLabel + ".AddP";

    if(pID <= 0 || pID > (int)tjs.pfps.size()) return false;
    auto& pfp = tjs.pfps[pID - 1];
    
    if(pfp.TPCID != ss3.TPCID) {
      mf::LogVerbatim("TC")<<fcnLabel<<" P"<<pID<<" is in the wrong TPC for 3S"<<ss3.ID;
      return false;
    }
    
    for(auto tid : pfp.TjIDs) {
      auto& tj = tjs.allTraj[tid - 1];
      // is this tj in a 2D shower that is in a 3D shower that is not this shower?
      if(tj.SSID > 0) {
        auto& ss = tjs.cots[tj.SSID - 1];
        if(ss.SS3ID > 0 && ss.SS3ID != ss3.ID) {
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Conflict: 3S"<<ss3.ID<<" adding P"<<pfp.ID<<" -> T"<<tid<<" is in 2S"<<tj.SSID<<" that is in 3S"<<ss.SS3ID<<" that is not this shower";
          return false;
        } // conflict
        // tj is in the correct 2D shower so nothing needs to be done
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" adding P"<<pfp.ID<<" -> T"<<tid<<" is in the correct shower 2S"<<tj.SSID;
        continue;
      } // pfp tj is in a shower
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<fcnLabel<<" 3S"<<ss3.ID<<" adding P"<<pfp.ID<<" -> T"<<tid;
        for(unsigned short ii = 0; ii < pfp.TjIDs.size(); ++ii) {
          if(pfp.TjIDs[ii] == tid) myprt<<" pfp TjCompleteness "<<std::fixed<<std::setprecision(2)<<pfp.TjCompleteness[ii];
        } // itj
        myprt<<" tj.SSID 2S"<<tj.SSID;
      } // prt
      // add it to the shower in the correct CTP
      for(auto& cid : ss3.CotIDs) {
        auto& ss = tjs.cots[cid - 1];
        if(ss.CTP != tj.CTP) continue;
        // Add it to the shower.
        AddTj(fcnLabel, tjs, tid, ss, doUpdate, prt);
        ss3.NeedsUpdate = true;
        break;
      } // cid
    } // tid
    
    if(doUpdate && ss3.NeedsUpdate) UpdateShower(fcnLabel, tjs, ss3, prt);
    return true;
    
  } // AddPFP

  ////////////////////////////////////////////////
  bool AddTj(std::string inFcnLabel, TjStuff& tjs, int tjID, ShowerStruct& ss, bool doUpdate, bool prt)
  {
    // Adds the Tj to the shower and optionally updates the shower variables
    
    if(tjID <= 0 || tjID > (int)tjs.allTraj.size()) return false;
    
    std::string fcnLabel = inFcnLabel + ".AddT";
    
    Trajectory& tj = tjs.allTraj[tjID - 1];
    
    if(tj.CTP != ss.CTP) {
      mf::LogVerbatim("TC")<<fcnLabel<<" T"<<tjID<<" is in the wrong CTP for 2S"<<ss.ID;
      return false;
    }
    
    if(tj.SSID > 0) {
      if(tj.SSID == ss.ID) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" T"<<tjID<<" is already in 2S"<<ss.ID;
        return true;
      } else {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Can't add T"<<tjID<<" to 2S"<<ss.ID<<". it is already used in 2S"<<tj.SSID;
        return false;
      }
    } // tj.SSID > 0

    ss.TjIDs.push_back(tjID);
    tj.SSID = ss.ID;
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
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added T"<<tj.ID<<" to 2S"<<ss.ID;
    ss.NeedsUpdate = true;
    
    if(doUpdate) return UpdateShower(fcnLabel, tjs, ss, prt);
    return true;
    
  } // AddTj
  
  ////////////////////////////////////////////////
  bool RemoveTj(std::string inFcnLabel, TjStuff& tjs, int TjID, ShowerStruct& ss, bool doUpdate, bool prt)
  {
    // Removes the Tj from a shower
    
    if(TjID > (int)tjs.allTraj.size()) return false;
    
    std::string fcnLabel = inFcnLabel + ".RTj";
    
    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[TjID - 1];

    if(tj.SSID != ss.ID) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Can't Remove T"<<TjID<<" from 2S"<<ss.ID<<" because it's not in this shower";
      // This isn't a failure
      return true;
    }
    tj.AlgMod[kShwrParent] = false;
    
    bool gotit = false;
    for(unsigned short ii = 0; ii < ss.TjIDs.size(); ++ii) {
      if(TjID == ss.TjIDs[ii]) {
        ss.TjIDs.erase(ss.TjIDs.begin() + ii);
        gotit = true;
        break;
      }
    } // ii
    if(!gotit) return false;
    tj.SSID = 0;
    // Removing a parent Tj?
    if(TjID == ss.ParentID) ss.ParentID = 0;
    // re-build everything?
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Remove T"<<TjID<<" from 2S"<<ss.ID;
    // removed the only tj
    if(ss.TjIDs.empty()) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Removed the last Tj. Killing 2S"<<ss.ID;
      MakeShowerObsolete(fcnLabel, tjs, ss, prt);
      return true;
    }
    // clear out the shower points to force a complete update when UpdateShower is next called
    ss.ShPts.clear();
    if(doUpdate) {
      ss.NeedsUpdate = true;
      return UpdateShower(fcnLabel, tjs, ss, prt);
    }
    return true;
  } // RemoveTj
  
  ////////////////////////////////////////////////
  bool FindNeutrinoParent(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // Look for a pfp in the shower that is the daughter of a neutrino and call it the parent
    
    std::string fcnLabel = inFcnLabel + ".FNuPar";
    
    int npid = 0;
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != ss3.TPCID) continue;
      if(pfp.PDGCode != 14) continue;
      npid = pfp.ID;
      break;
    } // pfp
    // no neutrino
    if(npid == 0) return false;
    auto& nuPFP = tjs.pfps[npid - 1];
    if(nuPFP.DtrIDs.empty()) {
      std::cout<<"Found neutrino with no daughters: "<<tjs.EventsProcessed<<" events processed\n";
      return false;
    }
    auto PIn3S = GetAssns(tjs, "3S", ss3.ID, "P");
    auto candDtrs = SetIntersection(nuPFP.DtrIDs, PIn3S);
    if(candDtrs.empty()) return false;
    int dtrID = 0;
    float bestFOM = 100;
    for(auto did : candDtrs) {
      auto& dtrPFP = tjs.pfps[did - 1];
      // ensure it isn't the parent of a shower already
      bool skipit = false;
      for(auto& oldSS3 : tjs.showers) {
        if(oldSS3.ID == 0) continue;
        if(did == oldSS3.ParentID) skipit = true;
      } // oldSS3
      if(skipit) continue;
      unsigned short pEnd = FarEnd(tjs, dtrPFP, ss3.ChgPos);
      float fom = ParentFOM(fcnLabel, tjs, dtrPFP, pEnd, ss3, prt);
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" -> nu P"<<npid<<" -> dtr P"<<did<<" fom "<<fom;
      if(fom < bestFOM) {
        bestFOM = fom;
        dtrID = did;
      }
    } // did
    if(dtrID == 0) return false;
    
    auto& dtrPFP = tjs.pfps[dtrID - 1];
    unsigned short pEnd = FarEnd(tjs, dtrPFP, ss3.ChgPos);
    float fom = ParentFOM(fcnLabel, tjs, dtrPFP, pEnd, ss3, prt);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" -> nu P"<<npid<<" -> dtr P"<<dtrID<<" fom "<<fom;
    
    
    // make local copies so we can recover from a failure
    auto oldSS3 = ss3;
    std::vector<ShowerStruct> oldSS(ss3.CotIDs.size());
    for(unsigned short ii = 0; ii < ss3.CotIDs.size(); ++ii) {
      oldSS[ii] = tjs.cots[ss3.CotIDs[ii] - 1];
    }
    
    if(SetParent(fcnLabel, tjs, dtrPFP, ss3, prt) && UpdateShower(fcnLabel, tjs, ss3, prt)) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" success";
      return true;
    }

    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" failed. Recovering";
    ss3 = oldSS3;
    for(unsigned short ii = 0; ii < oldSS.size(); ++ii) {
      auto& ss = oldSS[ii];
      tjs.cots[ss.ID - 1] = ss;
    } // ii
    std::cout<<fcnLabel<<" Failed to make P"<<dtrID<<" parent of 3S"<<ss3.ID<<". Need to remove it from the shower?\n";
    return false;
    
  } // FindNeutrinoParent
  
  ////////////////////////////////////////////////
  bool FindParent(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // look for a parent pfp for the shower.The 2D showers associated with it 
    // The parent should be at the start of the shower (shend = 0) if it is well-defined
    // (has small AspectRatio and small DirectionFOM). A search is also made for a parent at
    // the "wrong" end of the shower (shend = 1). The best one at the wrong end is used if
    // no parent is found at the shower start and the shower is poorly defined.
    //
    // This function returns false if there was a failure. Not finding a parent is not a failure
    
    if(ss3.ID == 0) return false;
    if(ss3.CotIDs.size() < 2) return false;
    
    std::string fcnLabel = inFcnLabel + ".FPar";
    MCParticleListUtils mcpu{tjs};
    int truPFP = mcpu.PrimaryElectronPFPID(ss3.TPCID);
    
    // look for a pfp in the shower that has a neutrino as a parent
    if(FindNeutrinoParent(fcnLabel, tjs, ss3, prt)) return true;
    
    double energy = ShowerEnergy(ss3);
    // the energy is probably under-estimated since there isn't a parent yet.
    energy *= 1.2;
    double shMaxAlong, along95;
    ShowerParams(energy, shMaxAlong, along95);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" Estimated energy "<<(int)energy<<" MeV shMaxAlong "<<shMaxAlong<<" along95 "<<along95<<" truPFP "<<truPFP;
    
    // look for the pfp that has a reasonable probability of being in the shower but with the
    // minimum along distance from the shower center. 
    // This image explains the concept. The *'s represents the points in 2D showers that define
    // the charge center in 3D, ChgPos. We are looking for a pfp parent denoted by ----. The end 
    // that is farthest from ChgPos is labeled P (pfp.XYZ[pend] in the code). The expected distance
    // from the shower start to shower Max, shMaxAlong, is found from ShowerParams. The longitudinal
    // and transverse distance of P relative to the shower center is alongTrans. The first cut on a
    // candidate parent is made requiring that D (alongTrans[0]) > 0.5 * shMaxAlong.
    //
    //       __________shend = 0________________      __________shend = 1________________
    //                                **********
    //                    ******  ******                   ******  ******
    //       P-----*****ChgPos******   **                *****ChgPos******-------P
    //                     ********     *******              *******
    //                     |----> along                          |---> along
    //       |<----D------>|                                     |<----D------>|
    // |<----shMaxAlong--->|                                     |<----shMaxAlong--->|
    //
    // Candidate parent ID for each end and the FOM 
    std::array<int, 2> parID {{0, 0}};
    std::array<float, 2> parFOM {{tjs.ShowerTag[8], tjs.ShowerTag[8]}};
    
    // temp vector to flag pfps that are already parents - indexed by ID
    std::vector<bool> isParent(tjs.pfps.size() + 1, false);
    for(auto& oldSS3 : tjs.showers) {
      if(oldSS3.ID == 0) continue;
      isParent[oldSS3.ParentID] = true;
    } // pfp
    
    // put the tjs associated with this shower in a flat vector
    auto TjsInSS3 = GetAssns(tjs, "3S", ss3.ID, "T");
    if(TjsInSS3.empty()) return false;
    
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      bool dprt = (pfp.ID == truPFP);
      if(pfp.TPCID != ss3.TPCID) continue;
      // ignore neutrinos
      if(pfp.PDGCode == 14 || pfp.PDGCode == 14) continue;
      // ignore shower pfps
      if(pfp.PDGCode == 1111) continue;
      // ignore existing parents
      if(isParent[pfp.ID]) continue;
      // check for inconsistent pfp - shower tjs 
      if(DontCluster(tjs, pfp.TjIDs, TjsInSS3)) continue;
      // ignore if the pfp energy is larger than the shower energy
      float pfpEnergy = 0;
      float minEnergy = 1E6;
      for(auto tid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tid - 1];
        float energy = ChgToMeV(tj.TotChg);
        pfpEnergy += energy;
        if(energy < minEnergy) minEnergy = energy;
      }
      pfpEnergy -= minEnergy;
      pfpEnergy /= (float)(pfp.TjIDs.size() - 1);
      if(dprt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" P"<<pfp.ID<<" E "<<pfpEnergy;
      if(pfpEnergy > energy) continue;
      // find the end that is farthest away
      unsigned short pEnd = FarEnd(tjs, pfp, ss3.ChgPos);
      auto pToS = PointDirection(pfp.XYZ[pEnd], ss3.ChgPos);
      double costh = DotProd(pToS, ss3.Dir);
      if(std::abs(costh) < 0.4) continue;
      // distance^2 between the pfp end and the shower start, charge center, and shower end
      float distToStart2 = PosSep2(pfp.XYZ[pEnd], ss3.Start);
      float distToChgPos2 = PosSep2(pfp.XYZ[pEnd], ss3.ChgPos);
      float distToEnd2 = PosSep2(pfp.XYZ[pEnd], ss3.End);
      if(dprt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" P"<<pfp.ID<<"_"<<pEnd<<" distToStart "<<sqrt(distToStart2)<<" distToChgPos "<<sqrt(distToChgPos2)<<" distToEnd "<<sqrt(distToEnd2);
      // find the end of the shower closest to the pfp
      unsigned short shEnd = 0;
      if(distToEnd2 < distToStart2) shEnd = 1;
      // This can't be a parent if the pfp end is closer to the shower center than the start or the end
      if(shEnd == 0 && distToChgPos2 < distToStart2) continue;
      if(shEnd == 1 && distToChgPos2 < distToEnd2) continue;
      if(dprt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<"_"<<shEnd<<" P"<<pfp.ID<<"_"<<pEnd<<" costh "<<costh;
      Point2_t alongTrans;
      // find the longitudinal and transverse components of the pfp start point relative to the
      // shower center
      FindAlongTrans(ss3.ChgPos, ss3.Dir, pfp.XYZ[pEnd], alongTrans);
      if(dprt) mf::LogVerbatim("TC")<<fcnLabel<<"   alongTrans "<<alongTrans[0]<<" "<<alongTrans[1];
      // find the probability this point is inside the shower. Offset by the expected
      // shower max distance. distToShowerMax will be > 0 if the pfp end is closer to
      // ChgPos than expected from the parameterization
      float distToShowerMax = shMaxAlong - std::abs(alongTrans[0]);
      float prob = InShowerProbLong(energy, distToShowerMax);
      if(dprt) mf::LogVerbatim("TC")<<fcnLabel<<"        prob "<<prob;
      if(prob < 0.1) continue;
/*
      // use the overall pfp direction instead of the starting direction. It may not be so
      // good if the shower develops quickly
      auto pfpDir = PointDirection(pfp.XYZ[pEnd], pfp.XYZ[1 - pEnd]);
      costh = DotProd(pfpDir, ss3.Dir);
      if(dprt) mf::LogVerbatim("TC")<<fcnLabel<<" pfpDir "<<pfpDir[0]<<" "<<pfpDir[1]<<" "<<pfpDir[2]<<" costh "<<costh;
      if(std::abs(costh) < 0.6) continue;
*/
      // find the parentFOM
      float candParFOM = ParentFOM(fcnLabel, tjs, pfp, pEnd, ss3, prt);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<fcnLabel;
        myprt<<" 3S"<<ss3.ID<<"_"<<shEnd;
        myprt<<" P"<<pfp.ID<<"_"<<pEnd;
        myprt<<" distToShowerMax "<<std::fixed<<std::setprecision(1)<<distToShowerMax<<" trans "<<alongTrans[1];
        myprt<<std::setprecision(2)<<" prob "<<prob;
        myprt<<" costh "<<costh;
        myprt<<" candParFOM "<<candParFOM;
      } // prt
      if(candParFOM < parFOM[shEnd]) {
        parFOM[shEnd] = candParFOM;
        parID[shEnd] = pfp.ID;
      } 
    } // pfp
    
    if(parID[0] == 0 && parID[1] == 0) return true;

    // decide which to use
    int bestPFP = 0;
    // find the average DirectionFOM to help decide
    float aveDirFOM = 0;
    float fom3D = 0;
    for(auto cid : ss3.CotIDs) aveDirFOM += tjs.cots[cid - 1].DirectionFOM;
    aveDirFOM /= (float)ss3.CotIDs.size();
    if(prt) {
      mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" parID[0] "<<parID[0]<<" fom "<<parFOM[0]<<" parID[1] "<<parID[1]<<" fom "<<parFOM[1]<<" aveDirFOM "<<aveDirFOM;
    }
    if(parID[0] > 0 && parID[1] > 0 && aveDirFOM > 0.3) {
      // candidates at both ends and the direction is not well known. Take
      // the one with the best FOM
      bestPFP = parID[0];
      fom3D = parFOM[0];
      if(parFOM[1] < parFOM[0]) {
        bestPFP = parID[1];
        fom3D = parFOM[1];
      }
    } else if(parID[0] > 0) {
      bestPFP = parID[0];
      fom3D = parFOM[0];
    } else {
      bestPFP = parID[1];
      fom3D = parFOM[1];
    }
    if(bestPFP == 0) return true;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" setting P"<<bestPFP<<" as the parent "<<fom3D;
    
    // make local copies so we can recover from a failure
    auto oldSS3 = ss3;
    std::vector<ShowerStruct> oldSS(ss3.CotIDs.size());
    for(unsigned short ii = 0; ii < ss3.CotIDs.size(); ++ii) {
      oldSS[ii] = tjs.cots[ss3.CotIDs[ii] - 1];
    }
    
    ss3.ParentID = bestPFP;
    auto& pfp = tjs.pfps[bestPFP - 1];
    unsigned short pend = FarEnd(tjs, pfp, ss3.ChgPos);
    ss3.Vx3ID = pfp.Vx3ID[pend];
    
    if(SetParent(fcnLabel, tjs, pfp, ss3, prt) && UpdateShower(fcnLabel, tjs, ss3, prt)) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" successful update";
      return true;
    }

    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" Failed. Recovering...";
    ss3 = oldSS3;
    for(unsigned short ii = 0; ii < oldSS.size(); ++ii) {
      auto& ss = oldSS[ii];
      tjs.cots[ss.ID - 1] = ss;
    } // ii
    return false;
    
  } // FindParent
  
  ////////////////////////////////////////////////
  bool SetParent(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, ShowerStruct3D& ss3, bool prt)
  {
    // set the pfp as the parent of ss3. The calling function should do the error recovery
    if(pfp.ID == 0 || ss3.ID == 0) return false;
    if(ss3.CotIDs.empty()) return false;
    
    std::string fcnLabel = inFcnLabel + ".SP";
    
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      stj.VtxID[0] = 0;
      if(ss.ParentID > 0) {
        auto& oldParent = tjs.allTraj[ss.ParentID - 1];
        oldParent.AlgMod[kShwrParent] = false;
        ss.ParentID = 0;
        ss.ParentFOM = 10;
      } // remove old parents
      // add new parents
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.CTP != ss.CTP) continue;
        if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjid) == ss.TjIDs.end()) {
          // Add the tj but don't update yet
          if(!AddTj(fcnLabel, tjs, tjid, ss, false, prt)) return false;
        } // parent not in ss
        // Don't define it to be the parent if the pfp projection in this plane is low
        unsigned short pEnd = FarEnd(tjs, pfp, ss3.ChgPos);
        auto tp = MakeBareTP(tjs, pfp.XYZ[0], pfp.Dir[pEnd], tj.CTP);
        if(tp.Delta > 0.5) {
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" parent P"<<pfp.ID<<" -> T"<<tjid<<" -> 2S"<<ss.ID<<" parent";
        } else {
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" parent P"<<pfp.ID<<" -> T"<<tjid<<" low projection in plane. Not a parent";
          continue;
        }
        ss.ParentID = tjid;
        ss.NeedsUpdate = true;
        // set the ss start vertex
        if(ss3.Vx3ID > 0) {
          auto& vx3 = tjs.vtx3[ss3.Vx3ID - 1];
          auto v2list = GetAssns(tjs, "3V", vx3.ID, "2V");
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] <= 0) continue;
            if(std::find(v2list.begin(), v2list.end(), tj.VtxID[end]) != v2list.end()) stj.VtxID[0] = tj.VtxID[end];
          } // end
        } // ss3.Vx3ID > 0
        // and update
        if(!UpdateShower(fcnLabel, tjs, ss, prt)) return false;
      } // tjid
    } // cid
    ss3.ParentID = pfp.ID;
    
    unsigned short pEnd = FarEnd(tjs, pfp, ss3.ChgPos);
    ss3.Vx3ID = pfp.Vx3ID[pEnd];
    float fom3D = ParentFOM(fcnLabel, tjs, pfp, pEnd, ss3, prt);
    for(auto cid : ss3.CotIDs) tjs.cots[cid - 1].ParentFOM = fom3D;

    return true;
  } // SetParent
  
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
  bool IsShowerLike(const TjStuff& tjs, const std::vector<int> TjIDs)
  {
    // Vote for the list of Tjs (assumed associated with a PFParticle) being shower-like
    if(TjIDs.empty()) return false;
    unsigned short cnt = 0;
    for(auto tid : TjIDs) {
      if(tid <= 0 || tid > (int)tjs.allTraj.size()) continue;
      if(tjs.allTraj[tid - 1].AlgMod[kShowerLike] > 0) ++cnt;
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
    if(showerEnergy < 10) {
      shMaxAlong = 0;
      along95 = 0;
      return;
    }
//    shMaxAlong = 7.0 * log(showerEnergy / 15);
    shMaxAlong = 16 * log(showerEnergy / 15);
    // The 95% containment is reduced a bit at higher energy
    double scale = 2.75 - 9.29E-4 * showerEnergy;
    if(scale < 2) scale = 2;
    along95 = scale * shMaxAlong;
  } // ShowerParams
  
  ////////////////////////////////////////////////
  double ShowerParamTransRMS(double showerEnergy, double along)
  {
    // returns the pareameterized width rms of a shower at along relative to the shower max
    double shMaxAlong, shE95Along;
    ShowerParams(showerEnergy, shMaxAlong, shE95Along);
    if(shMaxAlong <= 0) return 0;
    double tau = (along + shMaxAlong) / shMaxAlong;
    // The shower width is modeled as a simple cone that scales with tau
    double rms = -0.4 + 2.5 * tau;
    if(rms < 0.5) rms = 0.5;
    return rms;
  } // ShowerParamTransRMS

  ////////////////////////////////////////////////
  double InShowerProbLong(double showerEnergy, double along)
  {
    // Returns the likelihood that the point at position along (cm) is inside an EM shower
    // having showerEnergy (MeV). The variable along is relative to shower max.
    
    if(showerEnergy < 10) return 0;
    
    double shMaxAlong, shE95Along;
    ShowerParams(showerEnergy, shMaxAlong, shE95Along);
    // 50% of the shower energy is deposited between 0 < shMaxAlong < 1, which should be obvious considering
    // that is the definition of the shower max, so the probability should be ~1 at shMaxAlong = 1.
    // The Geant study shows that 95% of the energy is contained within 2.5 * shMax and has a small dependence 
    // on the shower energy, which is modeled in ShowerParams. This function uses a
    // sigmoid likelihood function is constructed with these constraints using the scaling variable tau 
    double tau = (along + shMaxAlong) / shMaxAlong;
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

    if(showerEnergy < 10) return 0;
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
  float InShowerProb(const TjStuff& tjs, const ShowerStruct3D& ss3, const PFPStruct& pfp)
  {
    // returns a likelihood (0 - 1) that the pfp particle belongs in shower ss3
    
    if(ss3.ID == 0 || pfp.ID == 0) return 0;
    float sum = 0;
    float cnt = 0;
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      if(ss.ID == 0) continue;
      for(auto tid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tid - 1];
        if(tj.CTP != ss.CTP) continue;
//        std::cout<<"3S"<<ss3.ID<<" P"<<pfp.ID<<" ";
        sum += InShowerProb(tjs, ss, tj);
        ++cnt;
      } // tid
    } //cid
    if(cnt == 0) return 0;
    return sum / cnt;

  } // InShowerProb
  
  ////////////////////////////////////////////////
  float InShowerProb(const TjStuff& tjs, const ShowerStruct& ss, const Trajectory& tj)
  {
    // returns a likelihood (0 - 1) that the tj particle belongs in shower ss
    // Keep it simple: construct a FOM, take the inverse and limit it to the range 0 - 1
    if(ss.ID == 0 || tj.ID == 0) return 0;
    if(ss.CTP != tj.CTP) return 0;
    
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return 0;
    unsigned short closePt1, closePt2;
    float doca = 1E6;
    TrajTrajDOCA(tjs, stj, tj, closePt1, closePt2, doca);
    if(doca == 1E6) return 0;
    float showerLen = PosSep(stj.Pts[0].Pos, stj.Pts[2].Pos);
    // make a rough separation cut. Return a small but non-zero value
    if(doca > 5 * showerLen) return 0.01;
    auto& stp = stj.Pts[closePt1];
    if(stp.DeltaRMS == 0) return 0;
    auto& ttp = tj.Pts[closePt2];
    Point2_t alongTrans;
    FindAlongTrans(stp.Pos, stp.Dir, ttp.Pos, alongTrans);
//    std::cout<<"ISP: 2S"<<ss.ID<<" T"<<tj.ID<<" showerLen "<<(int)showerLen<<" closePt "<<closePt1;
//    std::cout<<" closePt "<<closePt2<<" along "<<std::fixed<<std::setprecision(1)<<alongTrans[0]<<" "<<alongTrans[1];
    float rms = stp.DeltaRMS;
    if(rms < 1) rms = 1;
    float arg = alongTrans[1] / rms;
    float radProb = exp(-0.5 * arg * arg);
//    std::cout<<" rms "<<rms<<" radProb "<<std::setprecision(3)<<radProb;
    // This is a fake but may be OK if this function is called before the shower is well-defined
    rms = showerLen;
    arg = alongTrans[0] / rms;
    float longProb = exp(-0.5 * arg * arg);
//    std::cout<<" longProb "<<std::setprecision(3)<<longProb;
    float costh = std::abs(DotProd(stp.Dir, ttp.Dir));
//    std::cout<<" costh "<<std::setprecision(3)<<costh;
    float prob = radProb * longProb * costh;
//    std::cout<<" InShowerProb "<<std::setprecision(3)<<prob<<"\n";
    return prob;

  } // InShowerProb
  
  ////////////////////////////////////////////////
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, unsigned short pend, ShowerStruct3D& ss3, bool prt)
  {
    // Returns an average weighted parent FOM for all trajectories in the pfp being a parent of the 2D showers in ss3
    if(ss3.ID == 0) return 1000;
    float sum = 0;
    float wsum = 0;
    std::string fcnLabel = inFcnLabel + ".P3FOM";
    float dum1, dum2;
    for(auto cid : ss3.CotIDs) {
      auto& ss = tjs.cots[cid - 1];
      if(ss.ID == 0) continue;
      // look for the 3D matched tj in this CTP
      int tjid = 0;
      for(auto tid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tid - 1];
        if(tj.ID == 0) continue;
        if(tj.CTP == ss.CTP) tjid = tid;
      } // tid
      if(tjid == 0) continue;
      auto& ptj = tjs.allTraj[tjid - 1];
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      // determine which end is farthest away from the shower center
      unsigned short ptjEnd = FarEnd(tjs, ptj, stj.Pts[1].Pos);
      auto& farTP = ptj.Pts[ptj.EndPt[ptjEnd]];
      float chgCtrSep2 = PosSep2(farTP.Pos, stj.Pts[1].Pos);
      if(chgCtrSep2 < PosSep2(farTP.Pos, stj.Pts[0].Pos) && chgCtrSep2 < PosSep2(farTP.Pos, stj.Pts[2].Pos)) continue;
      float fom = ParentFOM(fcnLabel, tjs, ptj, ptjEnd, ss, dum1, dum2, prt);
      // ignore failures
      if(fom > 50) continue;
      // weight by the 1/aspect ratio
      float wt = 1;
      if(ss.AspectRatio > 0) wt = 1 / ss.AspectRatio;
      sum += wt * fom;
      wsum += wt;
    } // cid
    if(wsum == 0) return 100;
    float fom = sum / wsum;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 3S"<<ss3.ID<<" P"<<pfp.ID<<" fom "<<std::fixed<<std::setprecision(3)<<fom;
    return fom;
  } // ParentFOM
  
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
    
    if(ss.AspectRatio > 0.5) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" 2S"<<ss.ID<<" poor AspectRatio "<<ss.AspectRatio<<" FOM not calculated";
      return 100;
    }

    float fom = 0;
    float cnt = 0;
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    TrajPoint& stp0 = stj.Pts[0];
    // Shower charge center TP
    TrajPoint& stp1 = stj.Pts[1];
    // get the end that is farthest away from the shower center
    tjEnd = FarEnd(tjs, tj, stp1.Pos);
    // prospective parent TP
    TrajPoint& ptp = tj.Pts[tj.EndPt[tjEnd]];
    // find the along and trans components in WSE units relative to the
    // shower center
    Point2_t alongTrans;
    FindAlongTrans(stp1.Pos, stp1.Dir, ptp.Pos, alongTrans);
    // We can return here if the shower direction is well defined and
    // alongTrans[0] is > 0
    if(ss.AspectRatio < 0.2 && ss.DirectionFOM < 0.5 && alongTrans[0] > 0) return 100;
    tp1Sep = std::abs(alongTrans[0]);
    // Find the expected shower start relative to shower max (cm)
    double shMaxAlong, shE95Along;
    ShowerParams(ss.Energy, shMaxAlong, shE95Along);
    double alongcm = tjs.WirePitch * tp1Sep;
    // InShowerProbLong expects the longitudinal distance relative to shower max so it
    // should be < 0
    float prob = InShowerProbLong(ss.Energy, -alongcm);
    if(prob < 0.05) return 100;
    // The transverse position must certainly be less than the longitudinal distance
    // to shower max.
    if(alongTrans[1] > shMaxAlong) return 100;
    // longitudinal contribution to fom with 1 Xo error error (14 cm)
    float longFOM = std::abs(alongcm + shMaxAlong) / 14;
    fom += longFOM;
    ++cnt;
    // transverse contribution
    float transFOM = -1;
    if(stp0.DeltaRMS > 0) {
      transFOM = alongTrans[1] / stp0.DeltaRMS;
      fom += transFOM;
      ++cnt;
    }
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
      tjlist = GetAssns(tjs, "2V", vx2.ID, "T");
//      tjlist = GetVtxTjIDs(tjs, vx2);
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
    float chgFrcBtwFOM = (1 - chgFrac) / 0.05;
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
      myprt<<" along "<<std::fixed<<std::setprecision(1)<<alongTrans[0]<<" fom "<<longFOM;
      myprt<<" trans "<<alongTrans[1]<<" fom "<<transFOM;
      myprt<<" prob "<<prob;
      myprt<<" dang1 "<<dang1<<" fom "<<dang1FOM;
      myprt<<" dang2 "<<dang2<<" fom "<<dang2FOM;
      myprt<<" vx2Score "<<vx2Score<<" fom "<<vxFOM;
      myprt<<" chgFrac "<<chgFrac<<" fom "<<chgFracFOM;
      myprt<<" chgFracBtw "<<chgFracBtw<<" fom "<<chgFrcBtwFOM;
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
          if(DontCluster(tjs, ss1.TjIDs, ss2.TjIDs)) continue;
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
            // see if it looks like a muon in 3D
            if(tj.AlgMod[kMat3D]) {
              auto TInP = GetAssns(tjs, "T", tj.ID, "P");
              if(!TInP.empty()) {
                auto& pfp = tjs.pfps[TInP[0] - 1];
                if(pfp.PDGCode == 13 && MCSMom(tjs, pfp.TjIDs) > 500) continue;
              } // TInP not empty
            } // 3D matched
            if(AddTj(fcnLabel, tjs, tjID, ss1, false, prt)) doMerge = true;
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
    
    ChkAssns(fcnLabel, tjs);
    
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
          if(DontCluster(tjs, iss.TjIDs, jss.TjIDs)) continue;
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
                  if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Envelopes 2S"<<iss.ID<<" 2S"<<jss.ID<<" are close "<<PosSep(ivx, jvx)<<" cut "<<tjs.ShowerTag[2];
                  doMerge = true;
                  break;
                }
              } // jvx
              if(doMerge) break;
            } // ivx
          } // !domerge
          if(!doMerge) continue;
          // check the relative positions and angle differences. Determine which tps are the
          // closest. Don't merge if the closest points are at the shower start and the angle
          // difference is large
          unsigned short iClosePt = 0;
          unsigned short jClosePt = 0;
          float close = 1E6;
          auto& istj = tjs.allTraj[iss.ShowerTjID - 1];
          auto& jstj = tjs.allTraj[jss.ShowerTjID - 1];
          for(unsigned short ipt = 0; ipt < 3; ++ipt) {
            for(unsigned short jpt = 0; jpt < 3; ++jpt) {
              float sep = PosSep2(istj.Pts[ipt].Pos, jstj.Pts[jpt].Pos);
              if(sep < close) {
                close = sep;
                iClosePt = ipt;
                jClosePt = jpt;
              }
            } // jpt
          } // ipt
          float costh = DotProd(istj.Pts[0].Dir, jstj.Pts[0].Dir);
          if(iClosePt == 0 && jClosePt == 0 && costh < 0.955) {
            if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" showers are close at the start points with costh "<<costh<<". Don't merge";
            continue;
          }
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
    
    ChkAssns(fcnLabel, tjs);

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
      // ignore wimpy showers
      if(iss.Energy < 50) continue;
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
      bool skipit = DontCluster(tjs, iss.TjIDs, jss.TjIDs);
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" i2S"<<iss.ID<<" "<<PrintPos(tjs, tpList[ii].Pos)<<" j2S"<<jss.ID<<" "<<PrintPos(tjs, tpList[jj].Pos)<<" sepij "<<sepij<<" skipit? "<<skipit;
      if(skipit) continue;
      // draw a line between these points
      TrajPoint tp;
      MakeBareTrajPoint(tjs, tpList[ii], tpList[jj], tp);
//      PrintTrajPoint("ij", tjs, 0, 1, 0, tp);
      for(unsigned short kk = jj + 1; kk < sids.size(); ++kk) {
        auto& kss = tjs.cots[sids[kk] - 1];
        if(kss.ID == 0) continue;
        if(DontCluster(tjs, iss.TjIDs, kss.TjIDs)) continue;
        if(DontCluster(tjs, jss.TjIDs, kss.TjIDs)) continue;
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
        } // prt
        if(sepjk > minSep || delta > maxDelta) {
          // clear a short chain?
          if(chain.size() > 2) {
            // merge this chain
            int newID = MergeShowers(fcnLabel, tjs, chain, prt);
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<fcnLabel<<" merged chain";
              for(auto ssID : chain) myprt<<" 2S"<<ssID;
              myprt<<" -> 2S"<<newID;
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
/*
        if(newID > 0) {
          auto& newss = tjs.cots[newID - 1];
          if(AddTjsInsideEnvelope(fcnLabel, tjs, newss, prt)) UpdateShower(fcnLabel, tjs, newss, prt);
        }
*/
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<" merged chain";
          for(auto ssID : chain) myprt<<" "<<ssID;
          myprt<<" -> new ssID "<<newID;
        } // prt
      } // long chain
    } // ii
    
    ChkAssns(fcnLabel, tjs);
    
  } // MergeShowerChain
  
  ////////////////////////////////////////////////
  void MergeSubShowersTj(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // merge small showers that are downstream of shower-like tjs. This algorithm is written 
    // for low-energy showers with are likely to be sparse and poorly defined.
    
    if(!tjs.UseAlg[kMergeSubShowersTj]) return;
    
    std::string fcnLabel = inFcnLabel + ".MSSTj";
    
    struct TjSS {
      int ssID;
      int tjID;
      float dang;
    };
    std::vector<TjSS> tjss;
    
    // temp vector for DontCluster
    std::vector<int> tjid(1);
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.CTP != inCTP) continue;
      // TODO: Evaluate this cut
      if(ss.Energy > 300) continue;
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      auto stp0 = stj.Pts[0];
      float bestDang = 0.3;
      int bestTj = 0;
      // look for a Tj that has higher energy than the shower
      for(auto& tj : tjs.allTraj) {
        if(tj.AlgMod[kKilled]) continue;
        if(tj.CTP != ss.CTP) continue;
        // require that it isn't in any shower
        if(tj.SSID > 0) continue;
        // require it to be not short
        if(NumPtsWithCharge(tjs, tj, false) < 10) continue;
        // and satisfy the ShowerLike MCSMom cut. It is unlikely to be tagged shower-like
        if(tj.MCSMom > tjs.ShowerTag[1]) continue;
        // check consistency
        tjid[0] = tj.ID;
        if(DontCluster(tjs, tjid, ss.TjIDs)) continue;
        float tjEnergy = ChgToMeV(tj.TotChg);
        // find the end that is furthest away from the shower center
        unsigned short farEnd = FarEnd(tjs, tj, stj.Pts[1].Pos);
        // compare MCSMom at the far end and the near end
        unsigned short midpt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
        float mom1 = MCSMom(tjs, tj, tj.EndPt[farEnd], midpt);
        float mom2 = MCSMom(tjs, tj, tj.EndPt[1 - farEnd], midpt);
        float asym = (mom1 - mom2) / (mom1 + mom2);
        auto& farTP = tj.Pts[tj.EndPt[farEnd]];
        // IP btw the far end TP and the shower center
        float doca = PointTrajDOCA(tjs, stp0.Pos[0], stp0.Pos[1], farTP);
        float sep = PosSep(farTP.Pos, stp0.Pos);
        float dang = doca / sep;
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<" Candidate 2S"<<ss.ID<<" T"<<tj.ID<<"_"<<farEnd;
          myprt<<" ShEnergy "<<(int)ss.Energy<<" tjEnergy "<<(int)tjEnergy;
          myprt<<" doca "<<doca<<" sep "<<sep<<" dang "<<dang<<" asym "<<asym;
        }
        if(tjEnergy < ss.Energy) continue;
        if(asym < 0.5) continue;
        // TODO: This should be done more carefully
        // separation cut 100 WSE ~ 30 cm in uB
        if(sep > 100) continue;
        if(dang > bestDang) continue;
        bestDang = dang;
        bestTj = tj.ID;
      } // tj
      if(bestTj == 0) continue;
      TjSS match;
      match.ssID = ss.ID;
      match.tjID = bestTj;
      match.dang = bestDang;
      tjss.push_back(match);
    } // ss
    
    if(tjss.empty()) return;
    
    // ensure that a tj is only put in one shower
    bool keepGoing = true;
    while(keepGoing) {
      keepGoing = false;
      float bestDang = 0.3;
      int bestMatch = 0;
      for(unsigned short mat = 0; mat < tjss.size(); ++mat) {
        auto& match = tjss[mat];
        // already used
        if(match.dang < 0) continue;
        if(match.dang < bestDang) bestMatch = mat;
      } // mat
      if(bestMatch > 0) {
        auto& match = tjss[bestMatch];
        auto& ss = tjs.cots[match.ssID - 1];
        if(!AddTj(fcnLabel, tjs, match.tjID, ss, true, prt)) {
          if(prt) mf::LogVerbatim("TC")<<" Failed";
          continue;
        }
        match.dang = -1;
        // set the AlgMod bit
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        stj.AlgMod[kMergeSubShowersTj] = true;
        keepGoing = true;
      } // found bestMatch
    } // keepGoing
    
    ChkAssns(fcnLabel, tjs);
    
  } // MergeSubShowersTj
  
  ////////////////////////////////////////////////
  void MergeSubShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Merge small showers that are downstream of larger showers
    
    if(!tjs.UseAlg[kMergeSubShowers]) return;
    
    std::string fcnLabel = inFcnLabel + ".MSS";
    bool newCuts = (tjs.ShowerTag[0] == 4);
    constexpr float radLen = 14 / 0.3;
    
    if(prt) {
      if(newCuts) {
        mf::LogVerbatim("TC")<<fcnLabel<<" MergeSubShowers checking using ShowerParams";
      } else {
        mf::LogVerbatim("TC")<<fcnLabel<<" MergeSubShowers checking using radiation length cut ";
      }
    } // prt
    
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
          ShowerStruct& jss = tjs.cots[sortVec[jj].index];
          if(jss.ID == 0) continue;
          if(DontCluster(tjs, iss.TjIDs, jss.TjIDs)) continue;
          TrajPoint& jstp1 = tjs.allTraj[jss.ShowerTjID - 1].Pts[1];
          if(newCuts) {
            // find the longitudinal and transverse separation using the higher energy
            // shower which probably is better defined.
            Point2_t alongTrans;
            FindAlongTrans(istp1.Pos, istp1.Dir, jstp1.Pos, alongTrans);
            // the lower energy shower is at the wrong end of the higher energy shower if alongTrans[0] < 0
            if(alongTrans[0] < 0) continue;
            // increase the cut if the second shower is < 10% of the first shower
            float alongCut = along95;
            if(jss.Energy < 0.1 * iss.Energy) alongCut *= 1.5;
            float probLong = InShowerProbLong(iss.Energy, alongTrans[0]);
            float probTran = InShowerProbTrans(iss.Energy, alongTrans[0], alongTrans[1]);
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<fcnLabel<<" Candidate i2S"<<iss.ID<<" E = "<<(int)iss.Energy<<" j2S"<<jss.ID<<" E = "<<(int)jss.Energy;
              myprt<<" along "<<std::fixed<<std::setprecision(1)<<alongTrans[0]<<" trans "<<alongTrans[1];
              myprt<<" alongCut "<<alongCut<<" probLong "<<probLong<<" probTran "<<probTran;
            } // prt
            // TODO: fix ShowerParams so we can use the likelihood cut instead
            if(alongTrans[0] > alongCut) continue;
            if(alongTrans[1] > alongTrans[0]) continue;
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
    
    ChkAssns(fcnLabel, tjs);
    
  } // MergeSubShowers
  
  ////////////////////////////////////////////////
  int MergeShowers(std::string inFcnLabel, TjStuff& tjs, std::vector<int> ssIDs, bool prt)
  {
    // merge a list of showers and return the ID of the merged shower.
    // Returns 0 if there was a failure. 
    
    std::string fcnLabel = inFcnLabel + ".MS";
    if(ssIDs.size() < 2) return 0;
    // check for a valid ID
    for(auto ssID : ssIDs) if(ssID <= 0 || ssID > (int)tjs.cots.size()) return 0;
    // check for the same CTP and consistent assns
    int ss3Assn = 0;
    auto& ss0 = tjs.cots[ssIDs[0] - 1];
    std::vector<int> tjl;
    for(auto ssID : ssIDs) {
      auto& ss = tjs.cots[ssID - 1];
      if(ss.CTP != ss0.CTP) return 0;
      tjl.insert(tjl.end(), ss.TjIDs.begin(), ss.TjIDs.end());
      if(ss.SS3ID > 0 && ss3Assn == 0) ss3Assn = ss.SS3ID;
      if(ss.SS3ID > 0 && ss.SS3ID != ss3Assn) {
        std::cout<<fcnLabel<<" Assn conflict \n";
        return 0;
      }
    } // ssID
    // ensure the InShower Tjs are valid
    for(auto tjID : tjl) {
      auto& tj = tjs.allTraj[tjID - 1];
      if(tj.CTP != ss0.CTP || tj.AlgMod[kKilled]) {
        std::cout<<fcnLabel<<" bad InShower T"<<tjID<<"\n";
        return 0;
      }
    } // tjID
    
    // mark the old showers killed
    for(auto ssID : ssIDs) {
      auto& ss = tjs.cots[ssID - 1];
      ss.ID = 0;
      // kill the shower Tj
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      stj.AlgMod[kKilled] = true;
    } // tjID

    // in with the new
    auto newss = CreateSS(tjs, tjl);
    if(newss.ID == 0) return 0;
    
    for(auto tid : tjl) {
      auto& tj = tjs.allTraj[tid - 1];
      tj.SSID = newss.ID;
    } // tid
    newss.SS3ID = ss3Assn;
    
    // define the new shower
    if(!UpdateShower(fcnLabel, tjs, newss, prt)) {
      std::cout<<fcnLabel<<" UpdateShower failed\n";
      MakeShowerObsolete(fcnLabel, tjs, newss, prt);
      return 0;
    }
    // store it
    if(!StoreShower(fcnLabel, tjs, newss)) {
      std::cout<<fcnLabel<<" StoreShower failed\n";
      MakeShowerObsolete(fcnLabel, tjs, newss, prt);
      return 0;
    }
    return newss.ID;
    
  } // MergeShowers
  
  ////////////////////////////////////////////////
  bool MergeShowersAndStore(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, bool prt)
  {
    // Merge showers using shower indices. The icotID shower is modified in-place.
    // The jcotID shower is declared obsolete. This function also re-defines the shower and
    // sets the Parent ID to 0.
    
    if(icotID <= 0 || icotID > (int)tjs.cots.size()) return false;
    ShowerStruct& iss = tjs.cots[icotID - 1];
    if(iss.ID == 0) return false;
    if(iss.TjIDs.empty()) return false;
    if(iss.ShowerTjID <= 0) return false;
    
    if(jcotID <= 0 || jcotID > (int)tjs.cots.size()) return false;
    ShowerStruct& jss = tjs.cots[jcotID - 1];
    if(jss.TjIDs.empty()) return false;
    if(jss.ID == 0) return false;
    if(jss.ShowerTjID <= 0) return false;

    if(iss.CTP != jss.CTP) return false;
    
    std::string fcnLabel = inFcnLabel + ".MSAS";
    
    if(iss.SS3ID > 0 && jss.SS3ID > 0 && iss.SS3ID != jss.SS3ID) {
      std::cout<<fcnLabel<<" Error: 2S"<<iss.ID<<" and S"<<jss.ID<<" have different 2S -> 3S assns\n";
      return false;
    }
    
    Trajectory& itj = tjs.allTraj[iss.ShowerTjID - 1];
    Trajectory& jtj = tjs.allTraj[jss.ShowerTjID - 1];
    if(!itj.Pts[1].Hits.empty() || !jtj.Pts[1].Hits.empty()) {
      std::cout<<fcnLabel<<" Error: These shower Tjs have hits! T"<<itj.ID<<" T"<<jtj.ID<<"\n";
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
    // correct the assn
    for(auto tid : iss.TjIDs) {
      auto& tj = tjs.allTraj[tid - 1];
      tj.SSID = iss.ID;
    } // tid
    // transfer a 2S -> 3S assn
    if(iss.SS3ID == 0 && jss.SS3ID > 0) iss.SS3ID = jss.SS3ID;
    // merge the list of nearby Tjs
    iss.NearTjIDs.insert(iss.NearTjIDs.end(), jss.NearTjIDs.begin(), jss.NearTjIDs.end());
    // transfer the TruParentID if it is in jss
    if(jss.TruParentID > 0) iss.TruParentID = jss.TruParentID;
    // append the list of matched Tjs
    iss.ParentID = 0;
    iss.NeedsUpdate = true;
    // force a full update
    iss.ShPts.clear();
    jss.ID = 0;
    bool success = UpdateShower(fcnLabel, tjs, iss, prt);
    KillVerticesInShower(fcnLabel, tjs, iss, prt);
 
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
  bool AnalyzeRotPos(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    // The RotPos vector was filled and sorted by increasing distance along the shower axis.
    // This function divides the RotPos points into 3 sections and puts the transverse rms width in the
    // three sections into the shower Tj TrajPoint DeltaRMS variable. It also calculates the charge and number of shower
    // points closest to each TrajPoint. The 
    
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
      myprt<<" HitPos[0] "<<std::fixed<<std::setprecision(1);
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
    
    if(cotID > (int)tjs.cots.size()) return false;
    
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
    if(!AnalyzeRotPos(fcnLabel, tjs, ss, prt)) return false;

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
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    // Reverses the shower and the shower tj
    
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
    if(ss.DirectionFOM != 0) ss.DirectionFOM = 1 / ss.DirectionFOM;
    auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
    ReverseTraj(tjs, stj);
    DefineEnvelope(fcnLabel, tjs, ss, prt);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Reversed shower. Shower angle = "<<ss.Angle;
  } // ReverseShower
  
  ////////////////////////////////////////////////
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt)
  {
    // Reverses the shower and the shower tj
    
    if(cotID > (int)tjs.cots.size()) return;
    ShowerStruct& ss = tjs.cots[cotID - 1];
    if(ss.ID == 0) return;
    ReverseShower(inFcnLabel, tjs, ss, prt);

  }
  
  ////////////////////////////////////////////////
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // set the ss3 ID = 0 and remove 2D shower -> 3D shower associations. The 2D showers are not
    // declared obsolete
    for(auto cid : ss3.CotIDs) {
      if(cid == 0 || (unsigned short)cid > tjs.cots.size()) continue;
      auto& ss = tjs.cots[cid - 1];
      if(ss.SS3ID > 0 && ss.SS3ID != ss3.ID) {
        std::cout<<"MakeShowerObsolete:  3S"<<ss3.ID<<" -> 2S"<<ss.ID<<" SS3ID 3S"<<ss.SS3ID<<" != "<<ss3.ID<<"\n";
        continue;
      }
      ss.SS3ID = 0;
    } // cid
    if(prt) {
      std::string fcnLabel = inFcnLabel + ".MSO";
      mf::LogVerbatim("TC")<<fcnLabel<<" Killed  3S"<<ss3.ID;
    }
    if(ss3.PFPIndex < tjs.pfps.size()) {
      std::cout<<"MakeShowerObsolete:  3S"<<ss3.ID<<" -> P"<<ss3.PFPIndex+1<<" assn exists but maybe shouldn't...";
    }
    ss3.ID = 0;
  } // MakeShowerObsolete
  
  ////////////////////////////////////////////////
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    // Gracefully kills the shower and the associated shower Tj
    
    if(ss.ID == 0) return;
    
    std::string fcnLabel = inFcnLabel + ".MSO";
    
    auto& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
    if(!stp1.Hits.empty()) {
      std::cout<<fcnLabel<<" Trying to kill shower "<<ss.ID<<" that has hits associated with it. Don't do this...\n";
    }
    
    // clear a 3S -> 2S assn
    if(ss.SS3ID > 0 && ss.SS3ID <= (int)tjs.showers.size()) {
      auto& ss3 = tjs.showers[ss.SS3ID - 1];
      std::vector<int> newCIDs;
      for(auto cid : ss3.CotIDs) {
        if(cid != ss.ID) newCIDs.push_back(cid);
      } // cid
      ss3.CotIDs = newCIDs;
    } // ss3 assn exists
    
    // Kill the shower Tj if it exists. This also releases the hits
    if(ss.ShowerTjID > 0) MakeTrajectoryObsolete(tjs, ss.ShowerTjID - 1);
    
    // Restore the original InShower Tjs
    // Unset the killed bit
    for(auto& tjID : ss.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjID - 1];
      tj.AlgMod[kKilled] = false;
      // clear all of the shower-related bits
      tj.SSID = 0;
      tj.AlgMod[kShwrParent] = false;
      tj.AlgMod[kMergeOverlap] = false;
      tj.AlgMod[kMergeSubShowers] = false;
      tj.AlgMod[kMergeNrShowers] = false;
      tj.AlgMod[kMergeShChain] = false;
    } // tjID
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Killed 2S"<<ss.ID<<" and ST"<<ss.ShowerTjID;
    ss.ID = 0;
    // No reason to do this
//    ss.TjIDs.clear();
    
  } // MakeShowerObsolete
  
  ////////////////////////////////////////////////
  bool DontCluster(const TjStuff& tjs, const std::vector<int>& tjlist1, const std::vector<int>& tjlist2)
  {
    // returns true if a pair of tjs in the two lists are in the dontCluster vector
    if(tjlist1.empty() || tjlist2.empty()) return false;
    if(tjs.dontCluster.empty()) return false;
    for(auto tid1 : tjlist1) {
      for(auto tid2 : tjlist2) {
        int ttid1 = tid1;
        if(ttid1 > tid2) std::swap(ttid1, tid2);
        for(auto& dc : tjs.dontCluster) if(dc.TjIDs[0] == ttid1 && dc.TjIDs[1] == tid2) return true;
      } // dc
    } // tid1
    return false;
  } // DontCluster

  ////////////////////////////////////////////////
  void DefineDontCluster(TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // make a list of tj pairs that pass the cuts, but shouldn't be clustered for
    // different reasons, e.g. are likely parent tjs in different showers, are both attached
    // to the same high-score vertex, etc
    
    tjs.dontCluster.clear();

    DontClusterStruct dc;
    for(auto& vx3 : tjs.vtx3) {
      if(vx3.ID == 0) continue;
      if(vx3.TPCID != tpcid) continue;
      if(vx3.Score < tjs.Vertex2DCuts[7]) continue;
      auto PIn3V = GetAssns(tjs, "3V", vx3.ID, "P");
      if(PIn3V.size() < 2) continue;
//      Point3_t v3pos = {{vx3.X, vx3.Y, vx3.Z}};
      for(unsigned short ip1 = 0; ip1 < PIn3V.size() - 1; ++ip1) {
        auto& p1 = tjs.pfps[PIn3V[ip1]];
        // ignore the neutrino pfp
        if(p1.TjIDs.empty()) continue;
        unsigned short p1End = 0;
        if(p1.Vx3ID[1] == vx3.ID) p1End = 1;
        bool p1ShowerLike = IsShowerLike(tjs, p1.TjIDs);
//        float p1Sep = PosSep(p1.XYZ[p1End], v3pos);
        for(unsigned short ip2 = ip1 + 1; ip2 < PIn3V.size(); ++ip2) {
          auto& p2 = tjs.pfps[PIn3V[ip2]];
          if(p2.TjIDs.empty()) continue;
          unsigned short p2End = 0;
          if(p2.Vx3ID[1] == vx3.ID) p2End = 1;
          // Look for the case where an electron starts to shower close to the
          // vertex, creating a daughter that is also attached to the vertex. This
          // pair is OK to include in a shower. The signature is that the PFP doca between them is less
          // than the pfp - vx3 separation and both are shower like - something like this
          // where 3V is a 3D vertex
          //   \
          //    \ P3 (not shower-like) -> 3V
          //     3V ----------- P1 -> 3V (shower-like or NOT shower-like)
          //            ----------- P2 (shower-like doca closer to P1 than the vertex) -> 3V
          // The tjs in the P1 - P3 pair shouldn't be clustered
          // The tjs in the P2 - P3 pair shouldn't be clustered
          // The tjs in the P1 - P2 pair can be clustered
          bool p2ShowerLike = IsShowerLike(tjs, p2.TjIDs);
/*
          std::cout<<"DDC: P"<<p1.ID<<" p1ShowerLike "<<p1ShowerLike;
          std::cout<<" P"<<p2.ID<<" p2ShowerLike "<<p2ShowerLike<<"\n";
*/
          if(p1ShowerLike && p2ShowerLike) continue;
          // now enter the Tj pairs
          for(auto tid1 : p1.TjIDs) {
            auto& t1 = tjs.allTraj[tid1 - 1];
            for(auto tid2 : p2.TjIDs) {
              auto& t2 = tjs.allTraj[tid2 - 1];
              if(t1.CTP != t2.CTP) continue;
              dc.TjIDs[0] = tid1;
              dc.TjIDs[1] = tid2;
              if(dc.TjIDs[0] > dc.TjIDs[1]) std::swap(dc.TjIDs[0], dc.TjIDs[1]);
              dc.Vx2ID = vx3.Vx2ID[DecodeCTP(t1.CTP).Plane];
              dc.Vx3ID = vx3.ID;
              tjs.dontCluster.push_back(dc);
            } // tid2
          } // tid1
        } // ip2
      } // ip1
    } // vx3

  } // DefineDontCluster
  
  ////////////////////////////////////////////////
  void FindCots(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjLists, bool prt)
  {
    // Version 2 of TagShowerLike to try to improve the separation between close showers, e.g. from pi-zeros
    tjLists.clear();    
    if(tjs.ShowerTag[0] <= 0) return;

    // take the average of the low and high charge RMS range to use as a cut
    float typicalChgRMS = 0.5 * (tjs.ChargeCuts[1] + tjs.ChargeCuts[2]);

    // clear out old tags and make a list of Tjs to consider
    std::vector<int> tjids;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      tj.AlgMod[kShowerLike] = false;
      if(tj.AlgMod[kShowerTj]) continue;
      short npwc = NumPtsWithCharge(tjs, tj, false);
      // Don't expect any (primary) electron to be reconstructed as a single trajectory for
      // more than ~2 radiation lengths ~ 30 cm for uB ~ 100 wires 
      if(npwc > 100) continue;
      // check MCSMom for longish Tjs
      if(npwc > 5) {
        // Increase the MCSMom cut if the Tj is long and the charge RMS is high to reduce sensitivity 
        // to the fcl configuration. A primary electron may be reconstructed as one long Tj with large
        // charge rms and possibly high MCSMom or as several nearby shorter Tjs with lower charge rms
        float momCut = tjs.ShowerTag[1];
        if(tj.ChgRMS > typicalChgRMS) momCut *= tj.ChgRMS / typicalChgRMS;
        if(tj.MCSMom > momCut) continue;
      }
      // see if this tj is in a muon pfparticle that looks shower-like in this view
      if(tj.AlgMod[kMat3D]) {
        auto TInP = GetAssns(tjs, "T", tj.ID, "P");
        if(!TInP.empty()) {
          auto& pfp = tjs.pfps[TInP[0] - 1];
          if(pfp.PDGCode == 13 && MCSMom(tjs, pfp.TjIDs) > 500) continue;
        } // TInP not empty
      } // 3D-matched
      tjids.push_back(tj.ID);
    } // tj
    
    if(tjids.size() < 2) return;
    
    struct ClosePair {
      float doca;               // distance of closest approach between the tj pair
      int id1;       // id of the first tj
      int id2;       // id of the second tj
      unsigned short closePt1;  // index of the closest point on tj1
      unsigned short closePt2;  // index of the closest point on tj2
      bool used;                // set true when this pair is used
    };
    // min separation between tjs
    std::vector<ClosePair> cps;    
    // list of tjs that are close
    std::vector<int> closeTjs;
    for(unsigned short it1 = 0; it1 < tjids.size() - 1; ++it1) {
      Trajectory& t1 = tjs.allTraj[tjids[it1] - 1];
      bool t1TrackLike = (t1.MCSMom > tjs.ShowerTag[1]);
      auto T1InP = GetAssns(tjs, "T", t1.ID, "P");
      for(unsigned short it2 = it1 + 1; it2 < tjids.size(); ++it2) {
        Trajectory& t2 = tjs.allTraj[tjids[it2] - 1];
        // require one of them to be not tracklike
        bool t2TrackLike = (t2.MCSMom > tjs.ShowerTag[1]);
        if(t1TrackLike && t2TrackLike) continue;
        unsigned short ipt1, ipt2;
        float doca = tjs.ShowerTag[2];
        // Find the separation between Tjs without considering dead wires
        TrajTrajDOCA(tjs, t1, t2, ipt1, ipt2, doca, false);
        if(doca == tjs.ShowerTag[2]) continue;
        // see if they are close in 3D if that information exists
        if(!T1InP.empty()) {
          auto T2InP = GetAssns(tjs, "T", t2.ID, "P");
          if(!T2InP.empty()) {
            auto& p1 = tjs.pfps[T1InP[0] - 1];
            auto& p2 = tjs.pfps[T2InP[0] - 1];
            unsigned short closePt1, closePt2;
            float doca = PFPDOCA(p1, p2, closePt1, closePt2);
//            float costh = DotProd(p1.Dir[0], p2.Dir[0]);
//            std::cout<<"chk T"<<t1.ID<<" T"<<t2.ID<<" doca "<<doca<<" costh "<<costh<<"\n";
            if(doca > tjs.ShowerTag[2] * tjs.WirePitch) continue;
          } // !T2InP.empty()
        } // !T1InP.empty()

        // add a new one
        ClosePair cp;
        cp.doca = doca;
        cp.id1 = t1.ID;
        cp.closePt1 = ipt1;
        cp.id2 = t2.ID;
        cp.closePt2 = ipt2;
        cp.used = false;
        cps.push_back(cp);
        if(std::find(closeTjs.begin(), closeTjs.end(), t1.ID) == closeTjs.end()) closeTjs.push_back(t1.ID);
        if(std::find(closeTjs.begin(), closeTjs.end(), t2.ID) == closeTjs.end()) closeTjs.push_back(t2.ID);
      } // it2 (t2)
    } // it1 (t1)
    
    if(cps.empty()) return;

    // sort tjList by decreasing length
    std::vector<SortEntry> sortVec(closeTjs.size());
    for(unsigned short ii = 0; ii < closeTjs.size(); ++ii) {
      sortVec[ii].index = ii;
      auto& tj = tjs.allTraj[closeTjs[ii] - 1];
      sortVec[ii].length = PosSep(tj.Pts[tj.EndPt[0]].Pos, tj.Pts[tj.EndPt[1]].Pos);
    } // ii
    std::sort(sortVec.begin(), sortVec.end(), greaterThan);


    // cluster them starting with the longest
    // a temp vector for DontCluster
    std::vector<int> tmp(1);
    for(unsigned short ii = 0; ii < sortVec.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      auto& t1 = tjs.allTraj[closeTjs[indx] - 1];
      // already tagged?
      if(t1.AlgMod[kShowerLike]) continue;
//      float t1Len = sortVec[ii].length;
      // get the general direction
//      auto t1Dir = PointDirection(t1.Pts[t1.EndPt[0]].Pos, t1.Pts[t1.EndPt[1]].Pos);
      // start a list of clustered tjs
      std::vector<int> tlist;
      tlist.push_back(t1.ID);
      // try to add other close tjs
      bool added = true;
      while(added) {
        added = false;
        for(auto& cp : cps) {
          if(cp.used) continue;
          // is any ID in tlist equal to cp.id1 or cp.id2?
          bool isID1 = (std::find(tlist.begin(), tlist.end(), cp.id1) != tlist.end());
          bool isID2 = (std::find(tlist.begin(), tlist.end(), cp.id2) != tlist.end());
          if(!(isID1 || isID2)) continue;
          // determine which one is not in tlist and call it t2
          unsigned short t2id = cp.id1;
          if(isID1) t2id = cp.id2;
          auto& t2 = tjs.allTraj[t2id - 1];
          // already tagged?
          if(t2.AlgMod[kShowerLike]) continue;
          tmp[0] = t2.ID;
          if(DontCluster(tjs, tmp, tlist)) continue;
          // don't cluster if this tj is closer to another long tj
          bool isCloser = false;
          for(auto& pcp : cps) {
            if(t1.ID == pcp.id1 || t1.ID == pcp.id2) continue;
            if(!(t2.ID == pcp.id1 || t2.ID == pcp.id2)) continue;
            unsigned short oid = pcp.id1;
            if(oid == t2.ID) oid = pcp.id2;
            auto otj = tjs.allTraj[oid - 1];
            float otjLen = PosSep(otj.Pts[otj.EndPt[0]].Pos, otj.Pts[otj.EndPt[1]].Pos);
//            std::cout<<"tid1 T"<<t1.ID<<" tid2 T"<<t2.ID<<" oid T"<<oid<<"\n";
            if(pcp.doca < cp.doca && otjLen > 10) isCloser = true;
          } // pcp
          if(isCloser) continue;
          if(std::find(tlist.begin(), tlist.end(), t2.ID) != tlist.end()) continue;
//          std::cout<<"  add T"<<t2.ID<<" to "<<tjLists.size()<<"\n";
          tlist.push_back(t2.ID);
          // call it used
          cp.used = true;
          added = true;
        } // cp
      } // added
      if(tlist.size() > 1) {
        // tag them
        for(auto tjid : tlist) {
          auto& tj = tjs.allTraj[tjid - 1];
          tj.AlgMod[kShowerLike] = true;
        } // tjid
        // ignore wimpy cots (< 10 MeV)
        if(ShowerEnergy(tjs, tlist) < 10) continue;
        tjLists.push_back(tlist);
      } // tlist.size() > 1
    } // ii
/* This causes  problems later on
    // Check for leftover tjs and add them if they are shower-like
    for(auto tjid : closeTjs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.AlgMod[kShowerLike]) continue;
      std::vector<int> tlist(1, tjid);
      tj.AlgMod[kShowerLike] = true;
      // ignore wimpy cots (< 10 MeV)
      if(ShowerEnergy(tjs, tlist) < 10) continue;
      tjLists.push_back(tlist);
    } // tjid
*/
    // check consistency
    for(unsigned short ip = 0; ip < tjLists.size() - 1; ++ip) {
      auto& ilist = tjLists[ip];
      for(unsigned short jp = ip + 1; jp < tjLists.size(); ++jp) {
        auto& jlist = tjLists[jp];
        auto sij = SetIntersection(ilist, jlist);
        if(!sij.empty()) {
          std::cout<<"******** FindCots conflict:";
          for(auto tid : sij) std::cout<<" T"<<tid;
          std::cout<<" appears in multiple lists\n";
        }
      } // jp
    } // ip
/*
    std::cout<<"FindCots inCTP "<<inCTP<<" tjLists size "<<tjLists.size()<<"\n";
    for(unsigned short ip = 0; ip < tjLists.size(); ++ip) {
      auto& tlist = tjLists[ip];
      std::cout<<"ip "<<ip;
      for(auto tid : tlist) {
        std::cout<<" T"<<tid;
        auto plist = GetAssns(tjs, "T", tid, "P");
        if(!plist.empty()) std::cout<<"_P"<<plist[0];
      }
      std::cout<<"\n";
    } // ip
*/

  } // FindCots

  ////////////////////////////////////////////////
  void TagShowerLike(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP)
  {
    // Tag Tjs as InShower if they have MCSMom < ShowerTag[0] and there are more than
    // ShowerTag[6] other Tjs with a separation < ShowerTag[1].     
    
    if(tjs.ShowerTag[0] <= 0) return;
    
    if(tjs.allTraj.size() > 20000) return;
    
    // evaluate different cuts
    bool newCuts = (tjs.ShowerTag[0] > 2);
    float typicalChgRMS = 0.5 * (tjs.ChargeCuts[1] + tjs.ChargeCuts[2]);
    
    // clear out old tags and make a list of Tjs to consider
    std::vector<std::vector<int>> tjLists;
    std::vector<int> tjids;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      tj.AlgMod[kShowerLike] = false;
      if(tj.AlgMod[kShowerTj]) continue;
      // ignore Tjs with Bragg peaks
      bool skipit = false;
      for(unsigned short end = 0; end < 2; ++end) if(tj.StopFlag[end][kBragg]) skipit = true;
      if(skipit) continue;
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
        if(!newCuts) {
          if(len1 < len2 && len1 < doca) {
            if(len1 < doca) continue;
          } else {
            if(len2 < doca) continue;
          }
        } // !newCuts
        // found a close pair. See if one of these is in an existing cluster of Tjs
        bool inlist = false;
        for(unsigned short it = 0; it < tjLists.size(); ++it) {
          bool tj1InList = (std::find(tjLists[it].begin(), tjLists[it].end(), tj1.ID) != tjLists[it].end());
          bool tj2InList = (std::find(tjLists[it].begin(), tjLists[it].end(), tj2.ID) != tjLists[it].end());
          if(tj1InList || tj2InList) {
            // add the one that is not in the list
            if(!tj1InList) tjLists[it].push_back(tj1.ID);
            if(!tj2InList) tjLists[it].push_back(tj2.ID);
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
          tjLists.push_back(newlist);
        }
      } // it2
    } // it1
    if(tjLists.empty()) return;

    // mark them all as ShowerLike Tjs
    unsigned short nsh = 0;
    for(auto& tjl : tjLists) {
      nsh += tjl.size();
      for(auto& tjID : tjl) {
        auto& tj = tjs.allTraj[tjID - 1];
        tj.AlgMod[kShowerLike] = true;
        // unset flags
        tj.AlgMod[kSetDir] = false;
      } // tjid
    } // tjl
    
    // kill vertices with more than 1 shower-like tj that is close to the
    // vertex
    unsigned short nkill = 0;
    for(auto& vx2 : tjs.vtx) {
      if(vx2.ID == 0) continue;
      if(vx2.CTP != inCTP) continue;
      auto TInV2 = GetAssns(tjs, "2V", vx2.ID, "T");
      unsigned short nsl = 0;
      for(auto tid : TInV2) {
        auto& tj = tjs.allTraj[tid - 1];
        if(tj.AlgMod[kShowerLike]) {
          unsigned short nearEnd = 1 - FarEnd(tjs, tj, vx2.Pos);
          if(PosSep(tj.Pts[tj.EndPt[nearEnd]].Pos, vx2.Pos) < 6) ++nsl;
        }
      } // tid
      if(nsl < 2) continue;
      MakeVertexObsolete(tjs, vx2, true);
      ++nkill;
    } // vx2
    if(tjs.ShowerTag[12] >= 0) mf::LogVerbatim("TC")<<"TagShowerLike tagged "<<nsh<<" Tjs and killed "<<nkill<<" vertices in CTP "<<inCTP;

  } // TagShowerLike
  
  ////////////////////////////////////////////////
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    // Find Tjs that are near the shower but are not included in it
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
//      auto vxTjIDs = GetVtxTjIDs(tjs, vx);
      auto vxTjIDs = GetAssns(tjs, "2V", vx.ID, "T");
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
  void DefineEnvelope(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
  {
    
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    // shower Tj isn't fully defined yet
    if(stj.Pts[0].Pos[0] == 0) return;
    
    std::string fcnLabel = inFcnLabel + ".DE";
    
    ss.Envelope.resize(4);
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
      myprt<<fcnLabel<<" 2S"<<ss.ID<<" Envelope";
      for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      myprt<<" Area "<<(int)ss.EnvelopeArea;
      myprt<<" ChgDensity "<<ss.ChgDensity;
    }
    // This is the last function used to update a shower
    ss.NeedsUpdate = false;
  } // DefineEnvelope  
  
  ////////////////////////////////////////////////
  bool AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt)
   {
    // This function adds Tjs to the shower. It updates the shower parameters.
    
     if(ss.Envelope.empty()) return false;
     if(ss.ID == 0) return false;
     if(ss.TjIDs.empty()) return false;
     
     std::string fcnLabel = inFcnLabel + ".ATIE";

     if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Checking 2S"<<ss.ID;

     std::vector<int> tmp(1);
     unsigned short nadd = 0;
     for(auto& tj : tjs.allTraj) {
       if(tj.CTP != ss.CTP) continue;
       if(tj.AlgMod[kKilled]) continue;
       if(tj.SSID > 0) continue;
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
       // check consistency
       tmp[0] = tj.ID;
       if(DontCluster(tjs, tmp, ss.TjIDs)) continue;
       // See if both ends are outside the envelope
       bool end0Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[0]].Pos, ss.Envelope);
       bool end1Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[1]].Pos, ss.Envelope);
       if(!end0Inside && !end1Inside) continue;
       if(end0Inside && end1Inside) {
         // TODO: See if the Tj direction is compatible with the shower?
         if(AddTj(fcnLabel, tjs, tj.ID, ss, false, prt)) ++nadd;
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
       if(AddTj(fcnLabel, tjs, tj.ID, ss, false, prt)) {
         ++nadd;
       } else {
         if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" AddTj failed to add T"<<tj.ID;
       }
    } // tj
    
    if(nadd > 0) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added "<<nadd<<" trajectories ";
      ss.NeedsUpdate = true;
      UpdateShower(fcnLabel, tjs, ss, prt);
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
    if(cotID > (int)tjs.cots.size()) return;
    
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
    if(cotID > (int)tjs.cots.size()) return;
    
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
      if(killit) MakeShowerObsolete(fcnLabel, tjs, ss, prt);
      
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
      if(tj.SSID <= 0 || tj.AlgMod[kKilled]) {
        std::cout<<fcnLabel<<" Tj "<<tj.ID<<" isn't an inShower tj or is Killed "<<tj.AlgMod[kKilled]<<"\n";
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
        if(tjs.allTraj[itj].SSID <= 0) {
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
  double ShowerEnergy(const ShowerStruct3D& ss3)
  {
    if(ss3.ID == 0) return 0;
    if(ss3.Energy.empty()) return 0;
    double ave = 0;
    for(auto e : ss3.Energy) {
      ave += e;
    } // e
    ave /= ss3.Energy.size();
    return ave;
  } // ShowerEnergy

  ////////////////////////////////////////////////
  float ShowerEnergy(const TjStuff& tjs, const std::vector<int> tjIDs)
  {
    // Calculate energy using the total charge of all hits in each tj in the shower
    if(tjIDs.empty()) return 0;
    float sum = 0;
    for(auto tid : tjIDs) {
      auto& tj = tjs.allTraj[tid - 1];
      sum += tj.TotChg;
    } // tid
    return ChgToMeV(sum);
  } // ShowerEnergy
  
  ////////////////////////////////////////////////
  float ChgToMeV(float chg)
  {
    // Conversion from shower charge to energy in MeV. The calibration factor
    // was found by processing 500 pizero events with StudyPiZeros using StudyMode
    return 0.012 * chg;
  }
  
  ////////////////////////////////////////////////
  bool StoreShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3)
  {
    // Store a 3D shower. This function sets the 3S -> 2S assns using CotIDs and ensures
    // that the 2S -> 3S assns are OK. 
    
    std::string fcnLabel = inFcnLabel + ".S3S";
    if(ss3.ID <= 0) {
      std::cout<<fcnLabel<<" Invalid ID";
      return false;
    }
    if(ss3.CotIDs.size() < 2) {
      std::cout<<fcnLabel<<" not enough CotIDs";
      return false;
    }
    
    // check the 2S -> 3S assns
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.SS3ID == ss3.ID && std::find(ss3.CotIDs.begin(), ss3.CotIDs.end(), ss.ID) == ss3.CotIDs.end()) {
        std::cout<<fcnLabel<<" Bad assn: 2S"<<ss.ID<<" -> 3S"<<ss3.ID<<" but it's not inCotIDs.\n";
        return false;
      }
    } // ss
    
    // check the 3S -> 2S assns
    for(auto cid : ss3.CotIDs) {
      if(cid <= 0 || cid > (int)tjs.cots.size()) return false;
      auto& ss = tjs.cots[cid - 1];
      if(ss.SS3ID > 0 && ss.SS3ID != ss3.ID) {
        std::cout<<fcnLabel<<" Bad assn: 3S"<<ss3.ID<<" -> 2S"<<cid<<" but 2S -> 3S"<<ss.SS3ID<<"\n";
        return false;
      }
    } // cid
    
    // set the 2S -> 3S assns
    for(auto cid : ss3.CotIDs) tjs.cots[cid - 1].SS3ID = ss3.ID;
    
    tjs.showers.push_back(ss3);
    return true;
    
  } // StoreShower
  
  ////////////////////////////////////////////////
  bool StoreShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss)
  {
    // Store a ShowerStruct
    std::string fcnLabel = inFcnLabel + ".S2S";
    if(ss.ID <= 0) {
      std::cout<<fcnLabel<<" Invalid ID";
      return false;
    }
    if(ss.TjIDs.empty()) {
      std::cout<<fcnLabel<<" Fail: No TjIDs in 2S"<<ss.ID<<"\n";
      return false;
    }
    if(ss.ParentID > 0) {
      if(ss.ParentID > (int)tjs.allTraj.size()) {
        std::cout<<fcnLabel<<" Fail: 2S"<<ss.ID<<" has an invalid ParentID T"<<ss.ParentID<<"\n";
        return false;
      }
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), ss.ParentID) != ss.TjIDs.end()) {
        std::cout<<fcnLabel<<" Fail: 2S"<<ss.ID<<" ParentID is not in TjIDs.\n";
        return false;
      }
    } // ss.ParentID > 0

    // check the ID
    if(ss.ID != (int)tjs.cots.size() + 1) {
      std::cout<<fcnLabel<<" Correcting the ID 2S"<<ss.ID<<" -> 2S"<<tjs.cots.size() + 1;
      ss.ID = tjs.cots.size() + 1;
    }
    
    // set the tj shower bits
    for(auto& tjID : ss.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjID - 1];
      tj.SSID = ss.ID;
      tj.AlgMod[kShwrParent] = false;
      if(tj.ID == ss.ParentID) tj.AlgMod[kShwrParent] = true;
    } // tjID
    tjs.cots.push_back(ss);
    return true;
    
  } // StoreShower
  
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
  ShowerStruct CreateSS(TjStuff& tjs, const std::vector<int>& tjl)
  {
    // Create a shower and shower Tj using Tjs in the list
    ShowerStruct ss;
    
    if(tjl.empty()) {
      ss.ID = 0;
      return ss;
    }
    // Create the shower tj
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
    // define the ss
    ss.ID = tjs.cots.size() + 1;
    ss.CTP = stj.CTP;
    // assign all TJ IDs to this ShowerStruct
    ss.TjIDs = tjl;
    // declare them to be InShower
    for(auto tjid : tjl) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(tj.CTP != stj.CTP) {
        ss.ID = 0;
        return ss;
      }
      tj.SSID = ss.ID;
    } // tjid
    ss.ShowerTjID = stj.ID;
    ss.Envelope.resize(4);
    return ss;
    
  } // CreateSS
  
  ////////////////////////////////////////////////
  bool ChkAssns(std::string inFcnLabel, TjStuff& tjs)
  {
    // check tj - ss assns
    
    std::string fcnLabel = inFcnLabel + ".ChkAssns";
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      for(auto tid : ss.TjIDs) {
        auto& tj = tjs.allTraj[tid - 1];
        if(tj.SSID != ss.ID) {
          std::cout<<fcnLabel<<" ***** Error: 2S"<<ss.ID<<" -> TjIDs T"<<tid<<" != tj.SSID 2S"<<tj.SSID<<"\n";
          return false;
        }
      } // tid
      // check 2S -> 3S
      if(ss.SS3ID > 0 && ss.SS3ID <= (int)tjs.showers.size()) {
        auto& ss3 = tjs.showers[ss.SS3ID - 1];
        if(std::find(ss3.CotIDs.begin(), ss3.CotIDs.end(), ss.ID) == ss3.CotIDs.end()) {
          std::cout<<fcnLabel<<" ***** Error: 2S"<<ss.ID<<" -> 3S"<<ss.SS3ID<<" but the shower says no\n";
          return false;
        }
      } // ss.SS3ID > 0
    } // ss
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.SSID < 0) {
        std::cout<<fcnLabel<<" ***** Error: T"<<tj.ID<<" tj.SSID is fubar\n";
        tj.SSID = 0;
        return false;
      }
      if(tj.SSID == 0) continue;
      auto& ss = tjs.cots[tj.SSID - 1];
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
      std::cout<<fcnLabel<<" ***** Error: T"<<tj.ID<<" tj.SSID = 2S"<<tj.SSID<<" but the shower says no\n";
      return false;
    } // tj
    
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      for(auto cid : ss3.CotIDs) {
        auto& ss = tjs.cots[cid - 1];
        if(ss.SS3ID != ss3.ID) {
          std::cout<<fcnLabel<<" ***** Error: 3S"<<ss3.ID<<" -> 2S"<<cid<<" but it thinks it belongs to 3S"<<ss.SS3ID<<"\n";
          return false;
        }
      } // cid
    } // ss3
    return true;
  } // ChkAssns
  
  ////////////////////////////////////////////////
  void PrintShowers(std::string fcnLabel, TjStuff& tjs)
  {
    if(tjs.showers.empty()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<fcnLabel<<" 3D showers \n";
    for(auto& ss3 : tjs.showers) {
      myprt<<fcnLabel<<" 3S"<<ss3.ID<<" 3V"<<ss3.Vx3ID;
      myprt<<" parentID "<<ss3.ParentID;
      myprt<<" ChgPos"<<std::fixed;
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
        myprt<<" ST"<<stj.ID;
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
    
    bool printHeader = true;
    bool printExtras = false;
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      PrintShower(someText, tjs, ss, printHeader, printExtras);
      printHeader = false;
    } // ss
    if(tjs.ShowerTag[12] > 9) {
      // List of Tjs
      for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
        const auto& ss = tjs.cots[ict];
        if(!printAllCTP && ss.CTP != inCTP) continue;
        if(!printKilledShowers && ss.ID == 0) continue;
        myprt<<someText<<std::fixed;
        std::string sid = "2S" + std::to_string(ss.ID);
        myprt<<std::setw(5)<<sid;
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
      myprt<<std::setw(5)<<sid;
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
      myprt<<std::setw(5)<<sid;
      myprt<<" Nearby";
      for(auto id : ss.NearTjIDs) myprt<<" T"<<id;
      myprt<<"\n";
    } // ict
    // don't cluster list
    myprt<<"DontCluster";
    for(auto& dc : tjs.dontCluster) {
      if(dc.TjIDs[0] > 0) myprt<<" T"<<dc.TjIDs[0]<<"-T"<<dc.TjIDs[1];
    } // dc
    myprt<<"\nDontCluster";
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& iss = tjs.cots[ict];
      if(iss.ID == 0) continue;
      for(unsigned short jct = ict + 1; jct < tjs.cots.size(); ++jct) {
        const auto& jss = tjs.cots[jct];
        if(jss.ID == 0) continue;
        if(DontCluster(tjs, iss.TjIDs, jss.TjIDs)) myprt<<" 2S"<<iss.ID<<"-2S"<<jss.ID;
      } // jct
    } // ict
  } // Print2DShowers
  
  ////////////////////////////////////////////////
  void PrintShower(std::string someText, const TjStuff& tjs, const ShowerStruct& ss, bool printHeader, bool printExtras)
  {
    // print a single shower and a header (optional) and the extra variables like TjIDs, envelope, etc
    mf::LogVerbatim myprt("TC");
    
    if(printHeader) {
      myprt<<someText<<"   ID   CTP  ParID ParFOM TruParID Energy nTjs  dFOM AspRat  stj  vx0 __Pos0___ Chg(k) dRMS __Pos1___ Chg(k) dRMS __Pos2___ Chg(k) dRMS Angle SS3ID PFPID\n";
    } // printHeader

    myprt<<someText<<std::fixed;
    std::string sid = "2S" + std::to_string(ss.ID);
    myprt<<std::setw(5)<<sid;
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
    std::string vid = "NA";
    if(stj.VtxID[0] > 0) vid = "2V" + std::to_string(stj.VtxID[0]);
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
    if(ss.SS3ID > 0 && ss.SS3ID < (int)tjs.showers.size()) {
      auto& ss3 = tjs.showers[ss.SS3ID - 1];
      if(ss3.PFPIndex >= 0 && ss3.PFPIndex < tjs.pfps.size()) {
        std::string pid = "P" + std::to_string(ss3.PFPIndex + 1);
        myprt<<std::setw(6)<<pid;
      } else {
        myprt<<std::setw(6)<<"NA";
      }
    } else {
      myprt<<std::setw(6)<<"NA";
    }
    if(ss.NeedsUpdate) myprt<<" *** Needs update";
    
    if(!printExtras) return;
    myprt<<"\n";
    
    myprt<<someText<<std::fixed;
    sid = "2S" + std::to_string(ss.ID);
    myprt<<std::setw(5)<<sid;
    myprt<<" Tjs";
    for(auto id : ss.TjIDs) myprt<<" T"<<id;
    myprt<<"\n";
    myprt<<someText<<std::fixed;
    sid = "2S" + std::to_string(ss.ID);
    myprt<<std::setw(5)<<sid;
    myprt<<" Envelope";
    for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
    myprt<<"\n";
    myprt<<someText<<std::fixed;
    sid = "2S" + std::to_string(ss.ID);
    myprt<<std::setw(5)<<sid;
    myprt<<" Nearby";
    for(auto id : ss.NearTjIDs) myprt<<" T"<<id;

  } // PrintShower

} // namespace tca
