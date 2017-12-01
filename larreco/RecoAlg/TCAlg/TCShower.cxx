#include "larreco/RecoAlg/TCAlg/TCShower.h"

struct SortEntry{
  unsigned int index;
  float length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}


namespace tca {

  ////////////////////////////////////////////////
  bool Find3DShowerEndPoints(TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    
    if(ss3.ID == 0) return false;
    if(ss3.PFPIndex > tjs.pfps.size() - 1) return false;
    
    auto& pfp = tjs.pfps[ss3.PFPIndex];
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Inside F3DSEP: ss3 "<<ss3.ID<<" tjIDs";
      for(auto tjID : pfp.TjIDs) myprt<<" "<<tjID;
      myprt<<" start vtx ID "<<pfp.Vx3ID[0];
    }
    
    // See if the start end points are consistent
    std::vector<TrajPoint> spts;
    for(auto tjID : pfp.TjIDs) {
      unsigned short cotIndex = GetCotsIndex(tjs, tjID);
      if(cotIndex > tjs.cots.size() - 1) continue;
      auto& ss = tjs.cots[cotIndex];
      // Require good aspect ratio and knowledge of the direction
      if(ss.AspectRatio > tjs.ShowerTag[10] || ss.DirectionFOM > tjs.ShowerTag[9]) continue;
      auto& stj = tjs.allTraj[tjID - 1];
      // define dE/dx just in case it wasn't done before
      if(stj.dEdx[0] <= 0 && ss.ParentID > 0) {
        auto& ptj = tjs.allTraj[ss.ParentID - 1];
        unsigned short pend = FarEnd(tjs, ptj, ss);
        stj.dEdx[0] = ptj.dEdx[pend];
      }
      pfp.dEdx[0][DecodeCTP(stj.CTP).Plane] = stj.dEdx[0];
      TrajPoint stp = stj.Pts[0];
      if(prt) mf::LogVerbatim("TC")<<" start point "<<PrintPos(tjs, stp)<<" dir "<<std::fixed<<std::setprecision(2)<<stp.Dir[0]<<" "<<stp.Dir[1];
      spts.push_back(stp);
    } // tjID
    if(spts.size() < 2) {
      if(prt) mf::LogVerbatim("TC")<<" Couldn't find at least 2 2D showers with good AspectRatio and DirectionFOM";
      return false;
    }
//    if(prt) mf::LogVerbatim("TC")<<" stps size "<<spts.size();
    TVector3 pos, dir;
    if(!TrajPoint3D(tjs, spts[0], spts[1], pos, dir, prt)) {
      if(prt) mf::LogVerbatim("TC")<<"  TrajPoint3D failed. Maybe the shower direction is fubar";
      return false;
    }
    // set the start position using a 3D vertex if one exists
    if(pfp.Vx3ID[0] > 0) {
      auto& vx3 = tjs.vtx3[pfp.Vx3ID[0] - 1];
      pfp.XYZ[0][0] = vx3.X; pfp.XYZ[0][1] = vx3.Y; pfp.XYZ[0][2] = vx3.Z; 
    } else {
      // 
      pfp.XYZ[0][0] = pos[0]; pfp.XYZ[0][1] = pos[1]; pfp.XYZ[0][2] = pos[2]; 
    }
    pfp.Dir[0] = dir;
    
    // Now find the end point using the longest 2D shower
    double maxlen = 0;
    unsigned int maxID = 0;
    for(auto tjID : pfp.TjIDs) {
      auto& stj = tjs.allTraj[tjID - 1];
      float length = TrajLength(stj);
      if(length > maxlen) {
        maxlen = length;
        maxID = tjID;
      }
    } // tjID
    
    auto& longTj = tjs.allTraj[maxID - 1];
    geo::PlaneID planeID = DecodeCTP(longTj.CTP);
    pfp.BestPlane = planeID.Plane;
    ss3.BestPlane = planeID.Plane;
    double angleToVert = tjs.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
    double cosgamma = std::abs(std::sin(angleToVert) * pfp.Dir[0].Y() + std::cos(angleToVert) * pfp.Dir[0].Z());
    if(cosgamma == 0) return false;
    // convert maxlen from WSE units (1 wire spacing) to cm and find the 3D distance
    maxlen *= tjs.geom->WirePitch(planeID) / cosgamma;
    ss3.Len = maxlen;
    
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
      pfp.XYZ[1][ixyz] = pfp.XYZ[0][ixyz] + pfp.Dir[0][ixyz] * maxlen;
      ss3.Pos[ixyz] = pfp.XYZ[0][ixyz];
    }
    // Set the end direction to the start direction
    pfp.Dir[1] = pfp.Dir[0];
    ss3.Dir = pfp.Dir[0];
    ss3.OpenAngle = 0.05;
    if(prt) mf::LogVerbatim("TC")<<" ss3.Len "<<ss3.Len;
    
    return true;
  } // Find3DShowerEndPoints

  ////////////////////////////////////////////////
  void Finish3DShowers(TjStuff& tjs)
  {
    // Finish defining the shower properties that were not done previously.
    // Note to the reader: This code doesn't use MakeVertexObsolete to kill vertices using the assumption
    // that Finish3DShowers is being called after reconstruction is complete, in which case there is no
    // need to re-calculate the 2D and 3D vertex score which could potentially screw up the decisions that have
    // already been made.
    
    if(tjs.ShowerTag[0] < 2) return;

    // Transfer Tj hits from InShower Tjs to the shower Tj. This kills the InShower Tjs but doesn't consider how
    // this action affects vertices
    if(!TransferTjHits(tjs, false)) return;
    
    // Use special care if a neutrino PFParticle exists so that we don't clobber the neutrino vertex
    bool foundNeutrino = !tjs.pfps.empty() && (tjs.pfps[0].PDGCode == 12 || tjs.pfps[0].PDGCode == 14);
    // Associate shower Tj hits with 3D showers
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      for(unsigned short cotIndex = 0; cotIndex < ss3.CotIndices.size(); ++cotIndex) {
        auto& ss = tjs.cots[ss3.CotIndices[cotIndex]];
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
      } // cotIndex
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
          std::cout<<"Finish3DShowers wants to kill primary PFParticle "<<pfp.ID<<". Not doing it.\n";
          continue;
        }
        pfp.ID = 0;
        // kill the 3D vertex it is associated with
        for(unsigned short end = 0; end < 2; ++end) {
          if(pfp.Vx3ID[end] <=  0) continue;
          auto& vx3 = tjs.vtx3[pfp.Vx3ID[end] - 1];
          vx3.ID = 0;
        } // end
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
      for(auto tjid : vxtjs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.AlgMod[kKilled]) {
          vx3.ID = 0;
        }
      } // tjid
    } // vx3

  } // Finish3DShowers
  
  ////////////////////////////////////////////////
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // Find 2D showers using 3D-matched trajectories. This returns true if showers were found
    // which requires re-doing the 3D trajectory match
    
    if(tjs.ShowerTag[0] < 2) return false;
    
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
      PrintPFParticles("FSi", tjs);
      PrintAllTraj("FSi", tjs, debug, USHRT_MAX, 0);
    }

    // lists of Tj IDs in plane, (list1, list2, list3, ...)
    std::vector<std::vector<std::vector<int>>> bigList(tjs.NumPlanes);
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      std::vector<std::vector<int>> tjList;
      TagInShowerTjs(fcnLabel, tjs, inCTP, tjList, false);
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
        unsigned short cotIndex = Create2DShower(tjs, tjl);
        if(cotIndex == USHRT_MAX) continue;
        if(!DefineShower(fcnLabel, tjs, cotIndex, prt)) {
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failure from DefineShower "<<cotIndex<<" failed ";
          MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
          continue;
        }
        SaveTjInfo(tjs, inCTP, cotIndex, "DS");
        // Find nearby Tjs that were not included because they had too-high
        // MCSMom, etc. This will be used to decide if showers should be merged
        AddTjsInsideEnvelope(fcnLabel, tjs, cotIndex, prt);
        FindNearbyTjs(fcnLabel, tjs, cotIndex, prt);
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
      MergeSubShowers(fcnLabel, tjs, inCTP, prt);
      MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);
      SaveAllCots(tjs, inCTP, "MNrby");
      if(prt) Print2DShowers("Nrby", tjs, inCTP, false);
      for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
        auto& ss = tjs.cots[cotIndex];
        if(ss.ID == 0) continue;
        if(ss.CTP != inCTP) continue;
        if(AddTjsInsideEnvelope(fcnLabel, tjs, cotIndex, prt)) MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);
        if (tjs.SaveShowerTree) SaveAllCots(tjs, inCTP, "Merge");
      }
      SaveAllCots(tjs, inCTP, "ATj2");
      if(prt) Print2DShowers("ATIE", tjs, inCTP, false);
    } // plane
    if(tjs.cots.empty()) return false;
    
    prt = (dbgPlane > 2);
    // See if any of the primary PFParticles associated with a neutrino PFParticle
    // match with the 2D showers. 
    FindPrimaryShower("FS", tjs, prt);
    // Look for a 3D-matched parent in un-matched 2D showers
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      auto& ss = tjs.cots[cotIndex];
      if(ss.ID == 0) continue;
      FindExternalParent("FS", tjs, cotIndex, prt);
      SaveTjInfo(tjs, ss.CTP, cotIndex, "FEP");
      if(ss.ParentID == 0) FindStartChg(fcnLabel, tjs, cotIndex, prt);
    } // cotIndex
    CheckQuality(fcnLabel, tjs, tpcid, prt);
    
    SaveAllCots(tjs, "CQ");
    if(prt) Print2DShowers("FEP", tjs, USHRT_MAX, false);
    Match2DShowers(fcnLabel, tjs, tpcid, prt);
    if(prt) Print2DShowers("M2DS", tjs, USHRT_MAX, false);
    // Kill vertices in 2D showers that weren't matched in 3D
//    KillVerticesInShowers(fcnLabel, tjs, tpcid, prt);
    // Reconcile any differences by merging, etc and make the final set of
    // 3D showers and matching PFParticles. Note that the hits on InShower Tjs can't
    // be transferred to the shower until after hit merging has been done. This is done
    // in Finish3DShowers
    ReconcileShowers(fcnLabel, tjs, tpcid, prt);
    SaveAllCots(tjs, "RS");
     
    unsigned short nNewShowers = 0;
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      auto& ss = tjs.cots[cotIndex];
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      SaveTjInfo(tjs, ss.CTP, cotIndex, "Done");
     ++nNewShowers;
    } // cotIndex
    
    if(prt) Print2DShowers("FSo", tjs, USHRT_MAX, false);
    
    return (nNewShowers > 0);
    
  } // FindShowers3D
  
  ////////////////////////////////////////////////
  bool FindPrimaryShower(std::string inFcnLabel, TjStuff& tjs, bool prt)
  {
    // Try to match 2D showers to a PFParticle that is the daughter of a neutrino PFParticle
    if(tjs.cots.empty()) return false;
    if(tjs.pfps.empty()) return false;
    if(tjs.pfps[0].PDGCode != 14) return false;
    
    std::string fcnLabel = inFcnLabel + ".FPS";
    auto& neutrinoPFP = tjs.pfps[0];

    if(tjs.UseAlg[kKillShwrNuPFP]) {
      // ensure that what we think is the neutrino vertex is not inside a shower
      if(neutrinoPFP.Vx3ID[0] == 0 || neutrinoPFP.Vx3ID[0] > tjs.vtx3.size()) return false;
      auto& vx3 = tjs.vtx3[neutrinoPFP.Vx3ID[0] - 1];
      unsigned short ninsh = 0;
      std::array<float, 2> v2pos;
      for(unsigned short plane = 0; plane < tjs.NumPlanes; ++plane) {
        PosInPlane(tjs, vx3, plane, v2pos);
        for(auto& ss : tjs.cots) {
          if(ss.ID == 0) continue;
          CTP_t inCTP = EncodeCTP(vx3.TPCID.Cryostat, vx3.TPCID.TPC, plane);
          if(ss.CTP != inCTP) continue;
          // TODO: This maybe should be done more carefully
          bool insideEnvelope = PointInsideEnvelope(v2pos, ss.Envelope);
          if(!insideEnvelope) continue;
          // See if the vertex is close to the ends of the shower
          auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
          if(PosSep(stj.Pts[0].Pos, v2pos) < 5) continue;
          if(PosSep(stj.Pts[2].Pos, v2pos) < 5) continue;
//          std::cout<<fcnLabel<<" plane "<<plane<<" PFP vertex is inside shower ss "<<ss.ID<<"\n";
          ++ninsh;
        } // ss
      } // plane
      if(ninsh == tjs.NumPlanes) {
        // Vertex is inside a shower in all planes. Clobber the neutrino PFParticle
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" neutrino PFP is inside showers in all planes. Killing it ";
        neutrinoPFP.ID = 0;
        for(auto dtrid : neutrinoPFP.DtrIDs) {
          auto& dtr = tjs.pfps[dtrid - 1];
          if(dtr.ParentID != neutrinoPFP.ID) {
            std::cout<<fcnLabel<<" neutrino - daughter association error\n";
            continue;
          }
          dtr.ParentID = 0;
          dtr.Primary = false;
          for(auto tjid : dtr.TjIDs) {
            auto& dtj = tjs.allTraj[tjid - 1];
            dtj.AlgMod[kSetDir] = false;
          }
        }
        MakeVertexObsolete(tjs, vx3);
        return false;
      }
    } // kKillShwrNuPFP

    // Find the figure of merit for each trajectory in each neutrino daughter PFParticle to be the
    // parent of each 2D shower. Find the best FOM.
    float tp1sep, score;
    float bestFOM = tjs.ShowerTag[8];
    int shParentPFP = 0;
    int ssID = 0;
    for(auto pfpid : neutrinoPFP.DtrIDs) {
      // iterate through the daughter PFParticles
      auto& pfp = tjs.pfps[pfpid - 1];
      for(auto tjid : pfp.TjIDs) {
        // iterate through the daughter trajectories
        auto& tj = tjs.allTraj[tjid - 1];
        for(auto& ss : tjs.cots) {
          // iterate through the 2D showers
          if(ss.ID == 0) continue;
          if(ss.CTP != tj.CTP) continue;
          unsigned short end = 0;
          std::string pfpLabel = fcnLabel + "_" + std::to_string(pfp.ID);
          float fom = ParentFOM(pfpLabel, tjs, tj, end, ss, tp1sep, score, prt);
          if(fom < bestFOM) {
            bestFOM = fom;
            shParentPFP = pfpid;
            ssID = ss.ID;
          }
        } // cot
      } // tjid
    } // pfpid
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" bestFOM "<<bestFOM<<" pfp "<<shParentPFP<<" is a shower matched to ssid "<<ssID;
    if(shParentPFP == 0) return false;
    // We now know which PFParticle should be converted into a shower parent, hopefully one in each plane. Create
    // a 3D shower
    auto& showerParentPFP = tjs.pfps[shParentPFP - 1];
    ShowerStruct3D ss3;
    ss3.ID = tjs.showers.size() + 1;
    ss3.TPCID = showerParentPFP.TPCID;
    ss3.Pos[0] = showerParentPFP.XYZ[0][0];
    ss3.Pos[1] = showerParentPFP.XYZ[0][1];
    ss3.Pos[2] = showerParentPFP.XYZ[0][2];
    ss3.Dir = showerParentPFP.Dir[0];
    ss3.OpenAngle = 0.05;
    ss3.CotIndices.clear();
    ss3.Energy.resize(tjs.NumPlanes);
    ss3.EnergyErr.resize(tjs.NumPlanes);
    ss3.MIPEnergy.resize(tjs.NumPlanes);
    ss3.MIPEnergyErr.resize(tjs.NumPlanes);
    ss3.dEdx.resize(tjs.NumPlanes);
    ss3.dEdxErr.resize(tjs.NumPlanes);
    ss3.FOM = bestFOM;
    ss3.PFPIndex = shParentPFP - 1;
    ss3.Vx3ID = neutrinoPFP.Vx3ID[0];
    std::cout<<"FPS: ss3 "<<ss3.ID<<" ss3.PFPIndex "<<ss3.PFPIndex<<" ss3.Vx3ID "<<ss3.Vx3ID<<"\n";
    // Convert the muon neutrino to an electron neutrino
    neutrinoPFP.PDGCode = 12;
    // set the showerPFP pdgcode
    showerParentPFP.PDGCode = 1111;
    showerParentPFP.DtrIDs.clear();
    for(auto tjid : showerParentPFP.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      bestFOM = tjs.ShowerTag[8];
      ssID = 0;
      for(auto& ss : tjs.cots) {
        // iterate through the 2D showers
        if(ss.ID == 0) continue;
        if(ss.CTP != tj.CTP) continue;
        unsigned short end = 0;
        std::string pfpLabel = fcnLabel + "_" + std::to_string(shParentPFP);
        float fom = ParentFOM(pfpLabel, tjs, tj, end, ss, tp1sep, score, prt);
        if(fom < bestFOM) {
          bestFOM = fom;
          ssID = ss.ID;
        }
      } // ss
      if(ssID == 0) continue;
      auto& ss = tjs.cots[ssID - 1];
      ss.ParentID = tjid;
      // see if it is already in the shower
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjid) != ss.TjIDs.end()) {
        UpdateShowerWithParent(fcnLabel, tjs, ssID - 1, prt);
      } else if(!AddTj(fcnLabel, tjs, tjid, ssID - 1, true, prt)) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failed to add primary tj "<<tjid<<" to ss "<<ssID;
        ss.ParentID = 0;
        continue;
      }
      ss.SS3ID = ss3.ID;
      ss3.CotIndices.push_back(ssID - 1);
      unsigned short plane = DecodeCTP(tj.CTP).Plane;
      ss3.Energy[plane] = ss.Energy;
      ss3.EnergyErr[plane] = 0.3 * ss.Energy;
      ss3.dEdx[plane] = tj.dEdx[0];
      ss3.dEdxErr[plane] = 0.3 * tj.dEdx[0];
      // update the shower parent PFParticle and the neutrino PFParticle
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      std::replace(showerParentPFP.TjIDs.begin(), showerParentPFP.TjIDs.end(), tjid, stj.ID);
      stj.ParentID = tjid;
      stj.AlgMod[kMat3D] = true;
      stj.VtxID[0] = tj.VtxID[0];
      tj.PDGCode = 11;
    } // tjid
    
    if(!Find3DShowerEndPoints(tjs, ss3, prt)) {
      std::cout<<fcnLabel<<" Find3DShowerEndPoints failed\n";
      return false;
    }
    
    tjs.showers.push_back(ss3);
    std::cout<<"FPS: showers size "<<tjs.showers.size()<<" "<<ss3.ID;
    for(auto ci : ss3.CotIndices) std::cout<<" "<<ci;
    std::cout<<"\n";
    
    return true;
  } // FindPrimaryShower
  
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
    
    float fomCut = 1;
    
    // A simple matching scheme 
    for(unsigned short ci = 0; ci < tjs.cots.size() - 1; ++ci) {
      auto& iss = tjs.cots[ci];
      if(iss.ID == 0) continue;
      // already matched?
      if(iss.SS3ID > 0) continue;
      if(iss.TjIDs.empty()) continue;
      geo::PlaneID iplaneID = DecodeCTP(iss.CTP);
      if(iplaneID.Cryostat != tpcid.Cryostat) continue;
      if(iplaneID.TPC != tpcid.TPC) continue;
      Trajectory& istj = tjs.allTraj[iss.ShowerTjID - 1];
      for(unsigned short cj = ci + 1; cj < tjs.cots.size(); ++cj) {
        ShowerStruct& jss = tjs.cots[cj];
        if(jss.CTP == iss.CTP) continue;
        if(jss.ID == 0) continue;
        // already matched?
        if(iss.SS3ID > 0) continue;
        if(jss.SS3ID > 0) continue;
        if(jss.TjIDs.empty()) continue;
        geo::PlaneID jplaneID = DecodeCTP(jss.CTP);
        if(jplaneID.Cryostat != tpcid.Cryostat) continue;
        if(jplaneID.TPC != tpcid.TPC) continue;
        Trajectory& jstj = tjs.allTraj[jss.ShowerTjID - 1];
        TVector3 posij, dirij;
        // Use shower Tj point 0 which should yield the start point of the 3D shower
        if(!TrajPoint3D(tjs, istj.Pts[0], jstj.Pts[0], posij, dirij, prt)) continue;
        float fomij = Match3DFOM(fcnLabel, tjs, ci, cj, prt);
        if(fomij > fomCut) continue;
        if(tjs.NumPlanes == 2) {
          ShowerStruct3D ss3;
          ss3.ID = tjs.showers.size() + 1;
          ss3.TPCID = tpcid;
          ss3.Pos = posij;
          ss3.Dir = dirij;
          ss3.CotIndices.resize(2);
          ss3.CotIndices[0] = ci;
          ss3.CotIndices[1] = cj;
          ss3.Energy.resize(2);
          ss3.Energy[0] = tjs.cots[ci].Energy;
          ss3.Energy[1] = tjs.cots[cj].Energy;
          ss3.FOM = fomij;
          ss3.PFPIndex = USHRT_MAX;
          // don't fill or use the rest of the variables
          tjs.showers.push_back(ss3);
          iss.SS3ID = ss3.ID;
          jss.SS3ID = ss3.ID;
          if(prt) mf::LogVerbatim("TC")<<" new ss3 "<<ss3.ID<<" with fomij "<<fomij;
          continue;
        } // 2-plane TPC
        float bestFOM = fomCut;
        unsigned short bestck = USHRT_MAX;
        for(unsigned short ck = cj + 1; ck < tjs.cots.size(); ++ck) {
          ShowerStruct& kss = tjs.cots[ck];
          if(kss.CTP == iss.CTP || kss.CTP == jss.CTP) continue;
          if(kss.ID == 0) continue;
          if(kss.TjIDs.empty()) continue;
          if(kss.SS3ID > 0) continue;
          geo::PlaneID kplaneID = DecodeCTP(kss.CTP);
          if(kplaneID.Cryostat != tpcid.Cryostat) continue;
          if(kplaneID.TPC != tpcid.TPC) continue;
          Trajectory& kstj = tjs.allTraj[kss.ShowerTjID - 1];
          TVector3 posik, dirik;
          // Use shower Tj point 0 which should yield the start point of the 3D shower
          if(!TrajPoint3D(tjs, istj.Pts[0], kstj.Pts[0], posik, dirik, prt)) continue;
          float fomik = Match3DFOM(fcnLabel, tjs, ci, ck, prt);
          if(fomik > bestFOM) continue;
          TVector3 tmp = posik - posij;
          float sep = tmp.Mag();
          if(sep > 50) {
            if(prt) mf::LogVerbatim("TC")<<" Large stp[0] point separation "<<sep;
            continue;
          }
          bestFOM = fomik;
          bestck = ck;
        } // ck
        // 3-plane TPC below
        ShowerStruct3D ss3;
        ss3.ID = tjs.showers.size() + 1;
        ss3.TPCID = tpcid;
        // TODO: average posij and posik, etc here
        ss3.Pos = posij;
        ss3.Dir = dirij;
        iss.SS3ID = ss3.ID;
        jss.SS3ID = ss3.ID;
        // energy is entered by plane number
        ss3.Energy.resize(tjs.NumPlanes);
        ss3.FOM = bestFOM;
        if(bestck == USHRT_MAX) {
          // showers match in 2 planes
          ss3.CotIndices.resize(2);
          ss3.CotIndices[0] = ci;
          ss3.CotIndices[1] = cj;
          ss3.Energy[iplaneID.Plane] = tjs.cots[ci].Energy;
          ss3.Energy[jplaneID.Plane] = tjs.cots[cj].Energy;
          if(prt) mf::LogVerbatim("TC")<<" new 2-plane ss3 "<<ss3.ID<<" ssIDs "<<iss.ID<<" "<<jss.ID<<" with FOM "<<ss3.FOM;
        } else {
          // showers match in 3 planes
          unsigned short ck = bestck;
          ss3.CotIndices.resize(3);
          ss3.CotIndices[0] = ci;
          ss3.CotIndices[1] = cj;
          ss3.CotIndices[2] = bestck;
          ShowerStruct& kss = tjs.cots[ck];
          geo::PlaneID kplaneID = DecodeCTP(kss.CTP);
          ss3.Energy[iplaneID.Plane] = tjs.cots[ci].Energy;
          ss3.Energy[jplaneID.Plane] = tjs.cots[cj].Energy;
          ss3.Energy[kplaneID.Plane] = tjs.cots[ck].Energy;
          tjs.cots[ck].SS3ID = ss3.ID;
          if(prt) mf::LogVerbatim("TC")<<" new 3-plane ss3 "<<ss3.ID<<" ssIDs "<<iss.ID<<" "<<jss.ID<<" "<<tjs.cots[ck].ID<<" with FOM "<<ss3.FOM;
        }
        ss3.FOM = 0.5 * (fomij + bestFOM);
        if(bestck == USHRT_MAX && tjs.NumPlanes == 3 && FindMissingShowers1(fcnLabel, tjs, ss3, prt)) {
          // success
          if(prt) mf::LogVerbatim("TC")<<" Found missing showers";
        }
        // don't fill or use the rest of the variables yet
        tjs.showers.push_back(ss3);
      } // cj
    } // ci
    if(tjs.showers.empty()) return;
    
    // try to set Vx3ID
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      if(ss3.Vx3ID != 0) continue;
      unsigned short vid = USHRT_MAX;
      unsigned short cnt = 0;
      for(auto ci : ss3.CotIndices) {
        auto& ss = tjs.cots[ci];
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        if(stj.VtxID[0] > 0) {
          auto& vx2 = tjs.vtx[stj.VtxID[0] - 1];
          if(vx2.Vx3ID > 0) {
            if(vid == USHRT_MAX) {
              vid = vx2.Vx3ID;
              cnt = 1;
            } else if(vx2.Vx3ID == vid) {
              ++cnt;
            }
          }
        } // stj Vtx0 exists
      } // ci
      if(cnt == ss3.CotIndices.size()) ss3.Vx3ID = vid;
    } // ss3
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" 3D showers \n";
      for(auto& ss3 : tjs.showers) {
        myprt<<fcnLabel<<" "<<ss3.ID<<" Vx3ID "<<ss3.Vx3ID<<" cots\n";
        for(auto ci : ss3.CotIndices) {
          auto& ss = tjs.cots[ci];
          myprt<<fcnLabel<<"  "<<ss.ID;
          auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
          myprt<<" "<<stj.ID;
          myprt<<" "<<ss.ParentID;
          if(ss.ParentID > 0) {
            auto& ptj = tjs.allTraj[ss.ParentID - 1];
            myprt<<" ptjVtxID "<<ptj.VtxID[0];
          } else {
            myprt<<" Shower has no parent.";
          }
          myprt<<"\n";
        } // ci
      } // sss3
    } // prt

  } // Match2DShowers
  
  ////////////////////////////////////////////////
  bool FindMissingShowers2(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // The 2D shower referenced by cotIndex was not matched in 3D, presumably because it is low energy
    // and 2D showers were not found in the other two planes or alternatively because this blob of energy
    // was included in 2D showers in the other planes.
    // Look for a match using the previously 3D matched Tjs
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    auto& ss = tjs.cots[cotIndex];
    if(ss.ID == 0) return false;
    if(tjs.NumPlanes != 3) return false;
    
    std::string fcnLabel = inFcnLabel + ".FMS2";
    
    std::vector<unsigned short> pfpIndices;
    for(auto tjid : ss.TjIDs) {
      auto& tj = tjs.allTraj[tjid - 1];
      if(!tj.AlgMod[kMat3D]) continue;
      unsigned short pfpi = GetPFPIndex(tjs, tjid);
      if(pfpi == USHRT_MAX) continue;
      if(std::find(pfpIndices.begin(), pfpIndices.end(), pfpi) == pfpIndices.end()) pfpIndices.push_back(pfpi);
    } // tjid
    if(pfpIndices.empty()) return false;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" ss "<<ss.ID<<" pfpIndices";
      for(auto pfpi : pfpIndices) myprt<<" "<<pfpi;
    }
    
    // make a list of tjs in the other planes
    std::vector<std::vector<int>> tjlist(3);
    // and see if those Tjs are in showers
    std::vector<std::vector<int>> sslist(3);
    for(auto pfpi : pfpIndices) {
      auto& pfp = tjs.pfps[pfpi];
      if(pfp.ID == 0) continue;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        // ignore Tjs in this plane
        if(tj.CTP == ss.CTP) continue;
        unsigned short plane = DecodeCTP(tj.CTP).Plane;
        tjlist[plane].push_back(tjid);
        if(!tj.AlgMod[kInShower]) continue;
        for(auto ss : tjs.cots) {
          if(ss.ID == 0) continue;
          if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjid) != ss.TjIDs.end()) {
            sslist[plane].push_back(ss.ID);
            break;
          }
        } // ss
      } // tjid
    } // pfpi
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" Candidate Tjs\n";
      for(unsigned short plane = 0; plane < 3; ++plane) {
        if(tjlist[plane].empty()) continue;
        myprt<<" plane "<<plane<<" tjids:";
        for(auto tjid : tjlist[plane]) myprt<<" "<<tjid;
        myprt<<" ssids:";
        for(auto ssid : sslist[plane]) myprt<<" "<<ssid;
        myprt<<"\n";
      } // plane
    } // prt
    // See if these Tjs are in showers
    
    return false;
  } // FindMissingShowers2
  
  ////////////////////////////////////////////////
  bool FindMissingShowers1(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    // This function was called because only two showers were matched in a 3-plane TPC. Look for
    // missing showers in the third plane
    
    if(ss3.ID == 0) return false;
    if(ss3.CotIndices.size() == tjs.NumPlanes) return false;
    if(tjs.NumPlanes != 3) return false;
    
    std::string fcnLabel = inFcnLabel + ".FMS1";
    
    // determine which plane has the missed shower(s)
    std::array<bool, 3> hasShower {false};
    unsigned int cstat = 0;
    unsigned int tpc = 0;
    // check for 3D-matched parent Tjs
    unsigned short pfpIndex = USHRT_MAX;
    unsigned short cnt = 0;
    for(unsigned short ci = 0; ci < ss3.CotIndices.size(); ++ci) {
      auto& ss = tjs.cots[ss3.CotIndices[ci]];
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      unsigned short plane = planeID.Plane;
      cstat = planeID.Cryostat;
      tpc = planeID.TPC;
      hasShower[plane] = true;
      if(ss.ParentID == 0) continue;
      unsigned short pfpi = GetPFPIndex(tjs, ss.ParentID);
      if(pfpi == USHRT_MAX) continue;
      if(pfpIndex == USHRT_MAX) {
        pfpIndex = pfpi;
        cnt = 1;
      } else if(pfpi == pfpIndex) {
        ++cnt;
      }
     } // ci
    unsigned short missingPlane = 0;
    for(missingPlane = 0; missingPlane < 3; ++missingPlane) if(!hasShower[missingPlane]) break;
    if(missingPlane == tjs.NumPlanes) return false;
    CTP_t mCTP = EncodeCTP(cstat, tpc, missingPlane);
    // Find the start point in this plane
    TrajPoint mtp = MakeBareTrajPoint(tjs, ss3.Pos, ss3.Dir, mCTP);
    if(mtp.Pos[0] < 0) {
      mf::LogVerbatim("TC")<<fcnLabel<<" MakeBareTrajPoint failed";
      return false;
    }
    
    if(pfpIndex != USHRT_MAX && cnt < 2) {
      std::cout<<fcnLabel<<" 2D showers in this 3D shower have parent Tjs from different PFParticles\n";
      return false;
    }
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Look for missing showers in CTP "<<mCTP<<" pfpIndex "<<pfpIndex<<" cnt "<<cnt<<" starting at Pos "<<PrintPos(tjs, mtp.Pos)<<" dir "<<mtp.Dir[0]<<" "<<mtp.Dir[1];

    // make a list of showers that are close
    std::vector<int> cotIDs;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.CTP != mCTP) continue;
      // ignore already matched
      if(ss.SS3ID > 0) continue;
      // Find the IP between purported missing shower Tj and the shower center
      auto& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
      float delta = PointTrajDOCA(tjs, stp1.HitPos[0], stp1.HitPos[1], mtp);
      float sep = PosSep(stp1.HitPos, mtp.Pos);
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" chk ss "<<ss.ID<<" sep "<<sep<<" delta "<<delta;
      if(sep > 300 || delta > 20) continue;
      cotIDs.push_back(ss.ID);
    } // ss
    
    if(cotIDs.empty() && pfpIndex > 0 && cnt == 2) {
      std::cout<<fcnLabel<<" No matching showers found in CTP "<<mCTP<<" but 3D shower "<<ss3.ID<<" has 3D-matched parent Tjs\n";
      return false;
    }
    
    if(cotIDs.empty()) {
      std::cout<<fcnLabel<<" No matching showers found in CTP "<<mCTP<<".\n";
      return false;
    }
    
    // merge all of these showers
    int newID = MergeShowers(fcnLabel, tjs, cotIDs, prt);
    if(newID == 0) return false;
    
    auto& ss = tjs.cots[newID - 1];
    
    // see if a 3D-matched Tj in this plane should be declared the parent Tj of the new Tj
    ss.ParentID = 0;
    bool addedParent = false;
    if(pfpIndex != USHRT_MAX) {
      auto& pfp = tjs.pfps[pfpIndex];
      for(auto tjid : pfp.TjIDs) {
        auto& tj = tjs.allTraj[tjid - 1];
        if(tj.CTP == mCTP) ss.ParentID = tj.ID;
      } // tjid
      addedParent = UpdateShowerWithParent(fcnLabel, tjs, ss.ID - 1, prt);
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<ss.ParentID<<" is the parent of the new shower. Update success? "<<addedParent;
    }
    if(!addedParent) FindExternalParent(fcnLabel, tjs, ss.ID - 1, prt);
    if(ss.NeedsUpdate) DefineShower(fcnLabel, tjs, ss.ID - 1, prt);
    ss.SS3ID = ss3.ID;
    // add this to ss3
    ss3.CotIndices.push_back(ss.ID - 1);

    return true;
  } // FindMissingShowers1
  
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
  void ReconcileShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    // check the showers for consistency and reconcile differences
    
    if(tjs.showers.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".RS";

    if(prt) PrintPFParticles("RSi", tjs);

    // look for mis-matched shower parents
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      if(ss3.PFPIndex < tjs.pfps.size()) continue;
      std::vector<int> ptjs(ss3.CotIndices.size());
      std::vector<unsigned short> pfpis(ss3.CotIndices.size());
      unsigned short parentPFPIndex = USHRT_MAX;
      unsigned short cnt = 0;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" ss3.ID "<<ss3.ID;
      for(unsigned short ci = 0; ci < ss3.CotIndices.size(); ++ci) {
        auto& ss = tjs.cots[ss3.CotIndices[ci]];
        ptjs[ci] = ss.ParentID;
        auto& ptj = tjs.allTraj[ss.ParentID - 1];
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  ss.ID "<<ss.ID<<" ptj "<<ptj.ID<<" kMat3D? "<<ptj.AlgMod[kMat3D];
        if(!ptj.AlgMod[kMat3D]) continue;
        unsigned short pfpi = GetPFPIndex(tjs, ptj.ID);
        pfpis[ci] = pfpi;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"     PFPIndex "<<pfpi;
        if(pfpi == USHRT_MAX) continue;
        if(parentPFPIndex == USHRT_MAX) {
          parentPFPIndex = pfpi;
          cnt = 1;
        } else if(parentPFPIndex == pfpi) {
          ++cnt;
        }
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  ss.ID "<<ss.ID<<" ptj "<<ptjs[ci]<<" pfpi "<<pfpis[ci]<<" cnt "<<cnt;
      } // ci
      if(cnt < ss3.CotIndices.size()) {
        std::cout<<"Inconsistent PFPs for ss3 "<<ss3.ID<<" parent Tjs\n";
        mf::LogVerbatim("TC")<<"Inconsistent PFPs for ss3 "<<ss3.ID<<" parent Tjs";
      }
    } // ss3
    
    // look for missing 2D -> 3D matches
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      if(ss.SS3ID > 0) continue;
      if(!FindMissingShowers2(fcnLabel, tjs, ss.ID - 1, prt)) {
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        std::cout<<"Found incompletely matched 2D shower "<<ss.ID<<" at Pos "<<ss.CTP<<":"<<PrintPos(tjs, stj.Pts[1].Pos)<<"\n";
        mf::LogVerbatim("TC")<<"Found incompletely matched 2D shower "<<ss.ID<<" at Pos "<<ss.CTP<<":"<<PrintPos(tjs, stj.Pts[1].Pos);
      }
    } // ss

    // clobber PFParticles that have InShower tjs. Start by making a list of
    // all InShower Tjs that are in showers
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    std::vector<int> inShowerTjs;
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != cstat) continue;
      if(planeID.TPC != tpc) continue;
      inShowerTjs.insert(inShowerTjs.end(), ss.TjIDs.begin(), ss.TjIDs.end());
    } // ss
    // now do the clobbering
    for(auto& pfp : tjs.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.TPCID != tpcid) continue;
      for(auto tjid : pfp.TjIDs) {
        if(std::find(inShowerTjs.begin(), inShowerTjs.end(), tjid) != inShowerTjs.end()) {
          // remove the 3D match flag
          auto tj = tjs.allTraj[tjid - 1];
          tj.AlgMod[kMat3D] = false;
          pfp.ID = 0;
          break;
        }
      } // tjid
    } // pfp
    
    // Define the rest of the 3D shower struct (except for transferring the hits) 
    // and make a PFParticle for each one if none exists
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      if(ss3.TPCID != tpcid) continue;
      // Create a PFParticle?
      if(ss3.PFPIndex < tjs.pfps.size()) continue;
      ss3.Energy.resize(tjs.NumPlanes);
      ss3.EnergyErr.resize(tjs.NumPlanes);
      ss3.MIPEnergy.resize(tjs.NumPlanes);
      ss3.MIPEnergyErr.resize(tjs.NumPlanes);
      ss3.dEdx.resize(tjs.NumPlanes);
      ss3.dEdxErr.resize(tjs.NumPlanes);
      auto pfp = CreatePFPStruct(tjs, tpcid);
      // the list of of pfp TjIDs are for the Shower Tjs
      pfp.TjIDs.resize(ss3.CotIndices.size());
      float maxLen = 0;
      for(unsigned short cotIndex = 0; cotIndex < ss3.CotIndices.size(); ++cotIndex) {
        auto& ss = tjs.cots[ss3.CotIndices[cotIndex]];
        pfp.TjIDs[cotIndex] = ss.ShowerTjID;
        // set the match flag
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        stj.AlgMod[kMat3D] = true;
        unsigned short plane = DecodeCTP(ss.CTP).Plane;
        ss3.dEdx[plane] = stj.dEdx[0];
        pfp.dEdx[0][plane] = stj.dEdx[0];
        // dE/dx at the end is not defined
        // TODO: do this correctly
        ss3.dEdxErr[plane] = 0.3 * stj.dEdx[0];
        pfp.dEdxErr[0][plane] = ss3.dEdxErr[plane];
        ss3.Energy[plane] = ss.Energy;
        ss3.EnergyErr[plane] = 0.3 * ss.Energy;
        // just divide by 2.3 MeV/cm
        ss3.MIPEnergy[plane] = ss.Energy / 2.3;
        ss3.MIPEnergyErr[plane] = 0.3 * ss3.MIPEnergy[plane];
        // Set the best plane to be the longest one that has a reasonable start dE/dx
        float length = TrajLength(stj);
        if(length > maxLen && stj.dEdx[0] < 10) {
          maxLen = length;
          ss3.BestPlane = plane;
        }
      } // cotIndex
      pfp.Vx3ID[0] = ss3.Vx3ID;
      ss3.PFPIndex = tjs.pfps.size();
      pfp.PDGCode = 1111;
      tjs.pfps.push_back(pfp);
      if(!Find3DShowerEndPoints(tjs, ss3, prt)) {
        std::cout<<fcnLabel<<" Find3DShowerEndPoints failed\n";
        continue;
      }
      ss3.OpenAngle = 0.05;
      // change the neutrino PFParticle PDG code if this is a primary shower
      
    } // ss3
   
  } // ReconcileShowers
  
  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt)
  {
    float fom = 0;
    float cnt = 0;
    for(unsigned short ii = 0; ii < ss3.CotIndices.size() - 1; ++ii) {
      unsigned short ci = ss3.CotIndices[ii];
      for(unsigned short jj = ii + 1; jj < ss3.CotIndices.size(); ++jj) {
        unsigned short cj = ss3.CotIndices[jj];
        fom += Match3DFOM(inFcnLabel, tjs, ci, cj, prt);
        ++cnt;
      } // cj
    } // ci
    if(cnt == 0) return 100;
    return fom / cnt;
  } // Match3DFOM

  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs,
                   unsigned short icotIndex, unsigned short jcotIndex, unsigned short kcotIndex, bool prt)
  {
    if(icotIndex > tjs.cots.size() - 1) return 100;
    if(jcotIndex > tjs.cots.size() - 1) return 100;
    if(kcotIndex > tjs.cots.size() - 1) return 100;
    
    float ijfom = Match3DFOM(inFcnLabel, tjs, icotIndex, jcotIndex, prt);
    float jkfom = Match3DFOM(inFcnLabel, tjs, jcotIndex, kcotIndex, prt);
    
    return 0.5 * (ijfom + jkfom);
    
  } // Match3DFOM
  
  ////////////////////////////////////////////////
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, unsigned short icotIndex, unsigned short jcotIndex, bool prt)
  {
    // returns a Figure of Merit for a 3D match of two showers
    if(icotIndex > tjs.cots.size() - 1) return 100;
    if(jcotIndex > tjs.cots.size() - 1) return 100;
    
    auto& iss = tjs.cots[icotIndex];
    auto& istj = tjs.allTraj[iss.ShowerTjID - 1];    
    auto& jss = tjs.cots[jcotIndex];
    auto& jstj = tjs.allTraj[jss.ShowerTjID - 1];
    
    if(iss.CTP == jss.CTP) return 100;
    
    std::string fcnLabel = inFcnLabel + ".MFOM";
    
    float energyAsym = std::abs(iss.Energy - jss.Energy) / (iss.Energy + jss.Energy);
    
    geo::PlaneID iPlnID = DecodeCTP(iss.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jss.CTP);
    
    // compare match at the start of the showers
    float ix = tjs.detprop->ConvertTicksToX(istj.Pts[0].Pos[1] / tjs.UnitsPerTick, iPlnID);
    float jx = tjs.detprop->ConvertTicksToX(jstj.Pts[0].Pos[1] / tjs.UnitsPerTick, jPlnID);
    float pos0fom = std::abs(ix - jx) / 5;
    
    // compare match at the charge center
    ix = tjs.detprop->ConvertTicksToX(istj.Pts[1].Pos[1] / tjs.UnitsPerTick, iPlnID);
    jx = tjs.detprop->ConvertTicksToX(jstj.Pts[1].Pos[1] / tjs.UnitsPerTick, jPlnID);
    float pos1fom = std::abs(ix - jx) / 10;
    
    // and the end
    ix = tjs.detprop->ConvertTicksToX(istj.Pts[2].Pos[1] / tjs.UnitsPerTick, iPlnID);
    jx = tjs.detprop->ConvertTicksToX(jstj.Pts[2].Pos[1] / tjs.UnitsPerTick, jPlnID);
    float pos2fom = std::abs(ix - jx) / 20;
    
    float mfom = energyAsym * energyAsym;
    mfom += pos0fom * pos0fom;
    mfom += pos1fom * pos1fom;
    mfom += pos2fom * pos2fom;
    mfom = sqrt(mfom) / 4;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" iss.ID "<<iss.ID<<" jss.ID "<<jss.ID;
      myprt<<std::fixed<<std::setprecision(2);
      myprt<<" pos0fom "<<pos0fom;
      myprt<<" pos1fom "<<pos1fom;
      myprt<<" pos2fom "<<pos2fom;
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
  void FillPts(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    if(ss.ShowerTjID == 0) return;
    
    std::string fcnLabel = inFcnLabel + ".FP";
    
    unsigned int cnt = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      if(itj > tjs.allTraj.size() - 1) {
        mf::LogWarning("TC")<<"Bad TjID "<<ss.TjIDs[it];
        MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
        return;
      }
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.CTP != ss.CTP) {
        mf::LogWarning("TC")<<"Tj "<<tj.ID<<" is in the wrong CTP "<<tj.CTP<<" "<<ss.CTP;
        MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
        return;
      }
      if(tj.AlgMod[kShowerTj]) {
        mf::LogWarning("TC")<<fcnLabel<<" Tj "<<tj.ID<<" is in TjIDs in cotIndex "<<cotIndex<<" but is a ShowerTj! Killing it";
        MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
        return;
      }
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if(tp.UseHit[ii]) ++cnt;
      } // ipt
    } // it
    
    // Add any loose hits (those not in trajectory points) that are stashed in shower Tj Pt[0]
    TrajPoint& stp0 = tjs.allTraj[ss.ShowerTjID - 1].Pts[0];
    cnt += stp0.Hits.size();
    
    ss.ShPts.resize(cnt);
    
    // Now populate the vectors with the information we currently have 
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
    tjs.allTraj[ss.ShowerTjID - 1].AveChg = totChg;
    if(prt) mf::LogVerbatim("TC'")<<fcnLabel<<" cotIndex "<<cotIndex<<" filled "<<cnt<<" points including "<<stp0.Hits.size()<<" loose hits. Total charge "<<(int)totChg;
    
  } // FillPts

  ////////////////////////////////////////////////
  bool DefineShower(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Defines the properties of a shower using the trajectory points within the trajectories listed
    // in TjIDs. This wipes out any other information that may exist
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    if(ss.ShowerTjID == 0) return false;

    std::string fcnLabel = inFcnLabel + ".DS";
    
    if(ss.ParentID > 0) return UpdateShowerWithParent(fcnLabel, tjs, cotIndex, prt);
    
    FillPts(fcnLabel, tjs, cotIndex, prt);
    if(!FindChargeCenter(fcnLabel, tjs, cotIndex, prt)) {
      mf::LogWarning("TC")<<"Failed to find shower charge center";
      MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
      return false;
    }
    FindAngle(fcnLabel, tjs, cotIndex, prt);
    FillRotPos(fcnLabel, tjs, cotIndex, prt);
    if(!DefineShowerTj(fcnLabel, tjs, cotIndex, prt)) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failed to define Shower Tj. Killed the shower";
      MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
      return false;
    }
    FindNearbyTjs(fcnLabel, tjs, cotIndex, prt);
    DefineEnvelope(fcnLabel, tjs, cotIndex, prt);
    
    if(!tjs.MCPartList.empty()) {
      // get the truth if it exists
      auto& ss = tjs.cots[cotIndex];
      // Find the MC particle that matches with these InShower Tjs
      MCParticleListUtils tm{tjs};
      unsigned short nTruHits;
      unsigned short mcpIndex = tm.GetMCPartListIndex(ss, nTruHits);
      // Find the Tj that is closest to the start of this MC Particle
      if(mcpIndex != USHRT_MAX) ss.TruParentID = tm.MCParticleStartTjID(mcpIndex, ss.CTP);
    }

    return true;

  } // DefineShower
  
  ////////////////////////////////////////////////
  bool AddTj(std::string inFcnLabel, TjStuff& tjs, int tjID, unsigned short cotIndex, bool doUpdate, bool prt)
  {
    // Adds the Tj to the shower and optionally updates the shower variables
    
    if(tjID <= 0 || tjID > (int)tjs.allTraj.size()) return false;
    if(cotIndex > tjs.cots.size() - 1) return false;
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.ShowerTjID == 0) return false;
    
    std::string fcnLabel = inFcnLabel + ".ATj";

    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[tjID - 1];
    if(tj.AlgMod[kInShower]) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is already an InShower Tj";
      return false;
    }
    if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjID) != ss.TjIDs.end()) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is already in shower "<<ss.ID;
      return false;
    }
    if(tj.CTP != ss.CTP) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is in the wrong CTP "<<ss.ID;
      return false;
    }
    
    if(tj.ParentID == 0 && tj.ID != ss.ParentID) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is Primary and is not the shower parent "<<ss.ID;
      return false;
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
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Add Tj "<<tj.ID<<" parent "<<tj.ParentID<<" NeutrinoPrimaryTjID "<<NeutrinoPrimaryTjID(tjs, tj);
    
    if(doUpdate) return DefineShower(fcnLabel, tjs, cotIndex, prt);
    return true;
    
  } // AddTj
  
  ////////////////////////////////////////////////
  bool RemoveTj(std::string inFcnLabel, TjStuff& tjs, int TjID, unsigned short cotIndex, bool doUpdate, bool prt)
  {
    // Removes the Tj from a shower
    
    if(TjID > (int)tjs.allTraj.size()) return false;
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    std::string fcnLabel = inFcnLabel + ".RTj ";
    
    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[TjID - 1];
    if(!tj.AlgMod[kInShower]) return false;
    tj.AlgMod[kInShower] = false;
    tj.AlgMod[kShwrParent] = false;
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" RemoveTj "<<TjID<<" from cotIndex "<<cotIndex;
    if(doUpdate) return DefineShower(fcnLabel, tjs, cotIndex, prt);
    return true;
  } // RemoveTj

  ////////////////////////////////////////////////
  void FindExternalParent(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
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
     # 11 Debug in CTP (>10 debug cotIndex + 10)
     */

    if(cotIndex > tjs.cots.size() - 1) return;
    ShowerStruct& ss = tjs.cots[cotIndex];
    // Ensure that it is valid
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".FEP";
    
    // See if anything needs to be done
    if(ss.SS3ID > 0) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" ss "<<ss.ID<<" already matched in 3D";
      return;
    }
    
    // References to shower Tj points
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    int oldParent = ss.ParentID;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" ss.ID "<<ss.ID<<" stj.ID "<<stj.ID<<" Existing parent ID "<<oldParent<<" parent FOM "<<ss.ParentFOM;
      myprt<<" attached to vertex "<<stj.VtxID[0];
      myprt<<" end0 at "<<PrintPos(tjs, stj.Pts[0].Pos);
//      myprt<<" Tjs";
//      for(auto& tid : ss.TjIDs) myprt<<" "<<tid;
    }
    
    if(ss.AspectRatio > tjs.ShowerTag[10] || ss.DirectionFOM > tjs.ShowerTag[9]) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Don't search for a parent due to poor AspectRatio "<<ss.AspectRatio<<" or ss.DirectionFOM "<<ss.DirectionFOM;
      return;
    }
    
    float bestFOM = ss.ParentFOM;
    int imTheBest = 0;
    float bestVx3Score = 500;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled] && !tj.AlgMod[kInShower]) continue;
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
        useEnd = FarEnd(tjs, tj, ss);
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
      if(fom > bestFOM) continue;
      bestFOM = fom;
      imTheBest = tj.ID;
      bestVx3Score = vx3Score;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" current best "<<imTheBest<<" useEnd "<<useEnd<<" bestVx3Score "<<bestVx3Score;
    } // tj

    if(imTheBest < 1 || imTheBest > (int)tjs.allTraj.size()) return;
    
    if(bestFOM > tjs.ShowerTag[8]) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Best parent candidate FOM "<<bestFOM<<" exceeds the cut "<<tjs.ShowerTag[8];
      // Remove an old parent?
      if(oldParent > 0) {
        Trajectory& oldParentTj = tjs.allTraj[oldParent-1];
        // remove the parent flag
        oldParentTj.AlgMod[kShwrParent] = false;
        // remove it from the shower and update if it is attached to a high-score 3D vertex
//        if(oldParentTj.AlgMod[kTjHiVx3Score]) RemoveTj(fcnLabel, tjs, oldParent, cotIndex, true, prt);
      }
      ss.ParentID = 0;
      // detach the showerTj from a vertex if one exists
      if(prt && stj.VtxID[0] > 0) mf::LogVerbatim("TC")<<fcnLabel<<" Setting stj.VtxID[0] from "<<stj.VtxID[0]<<" to 0";
      stj.VtxID[0] = 0;
      return;
    }
    
    ss.ParentID = imTheBest;
    ss.ParentFOM = bestFOM;
    UpdateShowerWithParent(fcnLabel, tjs, cotIndex, prt);
    // do it again if Tjs were added during the update
    if(ss.NeedsUpdate) UpdateShowerWithParent(fcnLabel, tjs, cotIndex, prt);
    
  } // FindExternalParent
  
  ////////////////////////////////////////////////
  bool UpdateShowerWithParent(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // This updates all shower and shower Tj parameters when a new shower parent Tj is identified.
    if(cotIndex > tjs.cots.size() - 1) return false;
    ShowerStruct& ss = tjs.cots[cotIndex];
    // Ensure that everything is valid
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    if(ss.ParentID == 0) return false;
    if(ss.ShowerTjID == 0) return false;
    
    std::string fcnLabel = inFcnLabel + ".USWP";
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" ssID  "<<ss.ID<<" ParentID "<<ss.ParentID;
    
    if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), ss.ParentID) == ss.TjIDs.end()) {
      // add it to the Shower Tj but don't update. This will be done below
      if(!AddTj(fcnLabel, tjs, ss.ParentID, cotIndex, false, prt)) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failed to add the parent for some reason. Removing the Tj from the shower";
        // remove it and update the shower
        RemoveTj(fcnLabel, tjs, ss.ParentID, cotIndex, true, prt);
        ss.ParentID = 0;
        ss.ParentFOM = 99;
        return false;
      }
    } // not in the shower yet

    ss.NeedsUpdate = true;

    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    stj.AlgMod[kInShower] = true;
    stj.AlgMod[kShwrParent] = true;
    
    Trajectory& ptj = tjs.allTraj[ss.ParentID - 1];
    // determine the appropriate start point of the parent
    unsigned short pend = FarEnd(tjs, ptj, ss);
    auto& ptp = ptj.Pts[ptj.EndPt[pend]];
    // and the start point of the shower
/*
    unsigned short shend = 0;
    if(PosSep2(stj.Pts[2].Pos, ptp.Pos) < PosSep2(stj.Pts[0].Pos, ptp.Pos)) shend = 1;
    // ensure the shower and parent are in the same direction. Reverse the shower or the
    // parent Tj to be consistent with each other. This code assumes that the selection of
    // the parent Tj takes into account the shower direction FOM
    if(shend != pend) {
      // the ends are different
      if(ptj.AlgMod[kSetDir]) {
        // The primary tj direction has been set elsewhere, so reverse the shower
        ReverseShower(fcnLabel, tjs, cotIndex, prt);
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  Reversed the shower to be consistent with parent Tj";
      } else {
        // reverse the primary tj
        ReverseTraj(tjs, ptj);
        pend = shend;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"  Reversed the parent Tj to be consistent with shower";
      }
    } // shend != pend
*/
    if(!tjs.MCPartList.empty()) {
      // get the truth if it exists
      MCParticleListUtils tm{tjs};
      unsigned short nTruHits;
      unsigned short mcpIndex = tm.GetMCPartListIndex(ss, nTruHits);
      // Find the Tj that is closest to the start of this MC Particle
      if(mcpIndex != USHRT_MAX) ss.TruParentID = tm.MCParticleStartTjID(mcpIndex, ss.CTP);
    }

    // set the start vertex
    stj.VtxID[0] = ptj.VtxID[pend];
    // and dE/dx of the shower Tj
    stj.dEdx[0] = ptj.dEdx[pend];
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<< "  ParentID " << ss.ParentID << " dEdx " << stj.dEdx[0]<<" attached to vtx "<<stj.VtxID[0];

    // reference to the point on the parent Tj that is furthest away from the shower
    ss.Angle = ptp.Ang;
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
    
    // fill the rotated points
    FillRotPos(fcnLabel, tjs, cotIndex, prt);

    // AnalyzeRotPos calculates the charge (Chg), the charge center (HitPos) and
    // the transverse shower rms at each Tp in the shower Tj
    if(!AnalyzeRotPos(fcnLabel, tjs, cotIndex, prt)) {
      mf::LogVerbatim("TC")<<fcnLabel<< " Failure from AnalyzeRotPos. Killing this shower";
      MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
      return false;
    }
    
    if(ss.ID == 1) {
      std::cout<<fcnLabel<<" HitPos "<<PrintPos(tjs, stj.Pts[1].HitPos)<<" Pos "<<PrintPos(tjs, stj.Pts[1].Pos)<<"\n";
    }
    
    // define the angle of the shower Tj points
    for(auto& stp : stj.Pts) {
      stp.Ang = ptp.Ang;
      stp.Dir = ptp.Dir;
    } // stp
    
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

    
    if(ss.ID == 1) {
      std::cout<<fcnLabel<<" new "<<PrintPos(tjs, stj.Pts[1].HitPos)<<" Pos "<<PrintPos(tjs, stj.Pts[1].Pos)<<"\n";
    }

    FindNearbyTjs(fcnLabel, tjs, cotIndex, prt);
    // modify the shower tj points so that DefineEnvelope gives a reasonable result
    // TODO: Do this correctly, perhaps scaling by the aspect ratio
    if(stp0.DeltaRMS < 1) stp0.DeltaRMS = 1;
    float expectedRMS = 0.07 * PosSep(stp0.Pos, stp2.Pos);
//    std::cout<<"RMS "<<ss.ID<<" stp0 rms "<<stp0.DeltaRMS<<" Expected stp2 rms"<<expectedRMS<<" "<<stp2.DeltaRMS<<"\n";
    if(stp2.DeltaRMS < expectedRMS) stp2.DeltaRMS = expectedRMS;
    for(unsigned short nit = 0; nit < 2; ++nit) {
      DefineEnvelope(fcnLabel, tjs, cotIndex, prt);
      if(AddTjsInsideEnvelope(fcnLabel, tjs, cotIndex, prt)) ss.NeedsUpdate = true;
      if(!ss.NeedsUpdate) break;
    } // nit
    
    return true;
  } // UpdateShowerWithParent

  ////////////////////////////////////////////////
  unsigned short FarEnd(TjStuff& tjs, const Trajectory& tj, ShowerStruct& ss)
  {
    // Returns the end (0 or 1) of the Tj that is furthest away from the shower center
    if(ss.ShowerTjID == 0) return 0;
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID-1].Pts[1];
    unsigned short endPt0 = tj.EndPt[0];
    unsigned short endPt1 = tj.EndPt[1];
    if(PosSep2(tj.Pts[endPt1].Pos, stp1.Pos) > PosSep2(tj.Pts[endPt0].Pos, stp1.Pos)) return 1;
    return 0;
  } // FarEnd

  ////////////////////////////////////////////////
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, float& tp1Sep, float& vx3Score, bool prt)
  {
    // returns a FOM for the trajectory at the end point being the parent of ss and the end which
    // was matched.
    
    vx3Score = 0;
    tp1Sep = 0;
    
    if(tjEnd > 1) return 1000;
    if(ss.Energy == 0) return 1000;
    
    std::string fcnLabel = inFcnLabel + ".PFOM";
    
    // Radiation length converted to WSE units (for uB)
    constexpr float radLen = 14 / 0.3;
    constexpr float tenRadLen2 = 100 * radLen * radLen;

    if(ss.ID == 0) return 1000;
    if(ss.TjIDs.empty()) return 1000;
    if(ss.ShowerTjID == 0) return 1000;
    // get the end that is farthest away from the shower center
    tjEnd = FarEnd(tjs, tj, ss);
    // prospective parent TP
    unsigned short endPt = tj.EndPt[tjEnd];
    TrajPoint& ptp = tj.Pts[endPt];
    // Shower charge center TP
    unsigned short istj = ss.ShowerTjID - 1;
    TrajPoint& stp1 = tjs.allTraj[istj].Pts[1];
    tp1Sep = PosSep2(ptp.Pos, stp1.Pos);
    // Make a rough cut on radiation lengths
    if(tp1Sep > tenRadLen2) {
//      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" tjID "<<tj.ID<<" failed sep cut "<<(int)sqrt(tp1Sep)<<" 10 radiation lengths "<<(int)sqrt(tenRadLen2);
      return 100;
    }
    tp1Sep = sqrt(tp1Sep);
    
    // impact parameter between the projection of ptp and the charge center
    float delta = PointTrajDOCA(tjs, stp1.HitPos[0], stp1.HitPos[1], ptp);
    // make a rough cut
    if(delta > 100) return 50;
    
    // Estimate shower max. This parameterization comes from an Excel spreadsheet that uses the PDG shower max parameterization
    // from EGS4. Shower max, tmax, is calculated and used to generate a table of dE/dt vs t, which is then summed.
    // The value of t which yields 90% energy containment is used as a bound to find the shower center which is where the
    // shower center stp1 should be relative to the start of the parent. The parameterization is in units of the number of
    // radiation lengths with shEnergy in MeV.
    // Here is a summary table
    // E    t
    //----------
    //  70   1.90
    // 100   2.34
    // 150   2.68
    // 200   2.93
    // 250   3.10
    // 300   3.26
    // 350   3.34
    // 400   3.47
    // 600   3.79
    // 800   4.00
    //1000   4.11
    // Expected separation (cm) between the start of the parent trajectory and the shower charge center
    float expectedTPSep = 0.85 * log(3 * ss.Energy) - 2.65;
    // Convert it to the distance in WSE units 
    // We don't need great accuracy here because we don't know the projection of the shower in this view
    expectedTPSep *= radLen;
    // Assume that the projection of the shower in this view will be ~2/3
    expectedTPSep *= 0.6;
    // Guess that the RMS of the separation will be ~50% of the separation
    float expectedTPSepRMS = 0.8 * expectedTPSep;
    float sepPull = (tp1Sep - expectedTPSep) / expectedTPSepRMS;
    // The error on delta is probably dominated by not getting all of the shower Tjs and hits included
    // + missing energy (photons, etc). Just use the shower width at the shower center 
    float deltaErr = tjs.allTraj[istj].Pts[1].DeltaRMS;
    // protect against errors
    if(deltaErr < 0) deltaErr = 1;
    float deltaPull = delta / deltaErr;
    float dang = DeltaAngle(ptp.Ang, stp1.Ang);
    float dangErr = ss.AngleErr;
    if(dangErr < 0.1) dangErr = 0.1;
    // weight by the direction FOM?
    dangErr *= ss.DirectionFOM;
    float dangPull = dang / dangErr;
    // don't trust angle if shower is small
    if (ss.TjIDs.size() < 10) dangPull = 0;
    float mom = tj.MCSMom;
    if(mom > 500) mom = 500;
    float momPull = (mom - 500) / 100;
    // Pull due to the minimum separation between the other end of the parent Tj and the first shower Tj point.
    float tp0Sep2 = 0;
    float sep0Pull2 = 0;
    unsigned short otherEndPt = tj.EndPt[1-tjEnd];
    TrajPoint& optp = tj.Pts[otherEndPt];
    TrajPoint& stp0 = tjs.allTraj[istj].Pts[0];
    // don't bother taking the sqrt since it would be squared in a few lines anyway
    tp0Sep2 = PosSep2(optp.Pos, stp0.Pos);
    // expect this to be 1/2 of the separation between shower Tj point 0 and point 1
    float expectTp0Sep = 0.25 * PosSep2(stp0.Pos, stp1.Pos);
    sep0Pull2 = tp0Sep2 / expectTp0Sep;
    // primary Tj length
    float lenPull = (TrajLength(tj) - expectedTPSep) / expectedTPSepRMS;
    float fom = sqrt(sepPull * sepPull + deltaPull * deltaPull + dangPull * dangPull + momPull * momPull + sep0Pull2 + lenPull * lenPull);
    fom /= 6;
    if(tj.VtxID[tjEnd] > 0) {
      // check for a high-score 2D vertex that is outside the envelope with a high-score 3D vertex at this end
      VtxStore& vx2 = tjs.vtx[tj.VtxID[tjEnd] - 1];
      bool insideEnvelope = PointInsideEnvelope(vx2.Pos, ss.Envelope);
      if(!insideEnvelope && vx2.ID > 0 && vx2.Vx3ID > 0 && vx2.Vx3ID < tjs.vtx3.size() && vx2.Stat[kHiVx3Score]) {
        vx3Score = tjs.vtx3[vx2.Vx3ID - 1].Score;
        if(vx3Score > 0) fom /= sqrt(vx3Score);
      }
    }
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel;
      myprt<<" ssID "<<ss.ID;
      myprt<<" Tj "<<tj.ID<<" Pos "<<PrintPos(tjs, ptp);
/*
      myprt<<" VtxID "<<tj.VtxID[tjEnd];
      if(tj.VtxID[tjEnd] > 0) {
        VtxStore& vx2 = tjs.vtx[tj.VtxID[tjEnd] - 1];
        if(vx2.Vx3ID > 0) myprt<<" Vtx3ID "<<vx2.Vx3ID;
      }
*/
      myprt<<" end "<<tjEnd;
      myprt<<std::fixed<<std::setprecision(2);
      myprt<<" tp1Sep "<<std::fixed<<std::setprecision(1)<<tp1Sep;
      myprt<<" pull "<<sepPull;
      myprt<<" delta "<<delta<<" pull "<<deltaPull;
      myprt<<" dang "<<dang<<" pull "<<dangPull;
      myprt<<" mcsmom "<<(int)mom<<" pull "<<momPull;
      myprt<<" sep0 "<<sqrt(tp0Sep2)<<" pull "<<sqrt(sep0Pull2);
      myprt<<" length "<<(int)TrajLength(tj)<<" pull "<<lenPull;
      myprt<<" vx3Score "<<vx3Score;
      myprt<<" FOM "<<fom;
    }
/* Rarely used code
    if(tjs.ShowerTag[12] == -5) {
      // special output for creating an ntuple
      MCParticleListUtils tm{tjs};
      unsigned short nTruHits;
      unsigned short mcPtclIndex = tm.GetMCPartListIndex(ss, nTruHits);
      if(mcPtclIndex == 0 && sepPull < 4 && deltaPull < 10 && dangPull < 10) {
        // Print variables to create an ntuple
        float trueEnergy = tjs.MCPartList[mcPtclIndex]->E();
        mf::LogVerbatim myprt("TC");
        myprt<<"NTPL "<<ss.CTP<<" "<<std::fixed<<std::setprecision(2)<<trueEnergy;
        TrajPoint truTP;
        truTP.CTP = ss.CTP;
        tm.MakeTruTrajPoint(mcPtclIndex, truTP);
        myprt<<" "<<truTP.Ang<<" "<<truTP.Delta;
//        myprt<<" "<<tj.ID;
        // number of times this DefineShowerTj was called for this shower
        Trajectory& stj = tjs.allTraj[istj];
        myprt<<" "<<stj.Pass;
        // number of points (hits) in the shower
        myprt<<" "<<ss.ShPts.size();
        // shower charge
        myprt<<" "<<(int)stj.AveChg;
        myprt<<" "<<tp1Sep<<" "<<delta<<" "<<std::setprecision(3)<<dang<<" "<<(int)mom;
        myprt<<" "<<std::setprecision(2)<<sqrt(tp0Sep2)<<" "<<TrajLength(tj)<<" "<<fom;
        if(tj.ID == ss.TruParentID) {
          myprt<<" 1";
        } else {
          myprt<<" 0";
        }
      }
    } // tjs.ShowerTag[12] == -5
*/
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
    
    // merges showers if they share nearby Tjs
    for(unsigned short ci1 = 0; ci1 < tjs.cots.size() - 1; ++ci1) {
      ShowerStruct& ss1 = tjs.cots[ci1];
      if(ss1.CTP != inCTP) continue;
      if(ss1.ID == 0) continue;
      if(ss1.TjIDs.empty()) continue;
      for(unsigned short ci2 = ci1 + 1; ci2 < tjs.cots.size(); ++ci2) {
        ShowerStruct& ss2 = tjs.cots[ci2];
        if(ss2.CTP != inCTP) continue;
        if(ss2.ID == 0) continue;
        if(ss2.TjIDs.empty()) continue;
        std::vector<int> shared;
        std::set_intersection(ss1.NearTjIDs.begin(), ss1.NearTjIDs.end(), 
                              ss2.NearTjIDs.begin(), ss2.NearTjIDs.end(), std::back_inserter(shared));
        if(shared.empty()) continue;
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<" Merge ss2 "<<ss2.ID<<" into "<<ss1.ID<<"? shared nearby Tjs:";
          for(auto tjid : shared) myprt<<" "<<tjid;
        }
        // add the shared Tjs to ss1 if they meet the requirements
        // ensure that the shower isn't InShower already
        unsigned short nadd = 0;
        for(auto& tjID : shared) {
          auto& tj = tjs.allTraj[tjID - 1];
          if(tj.AlgMod[kInShower]) continue;
          // ignore long muons
          if(tj.PDGCode == 13 && tj.Pts.size() > 100) continue;
          if(AddTj(fcnLabel, tjs, tjID, ci1, false, prt)) ++nadd;
        } // tjID
        if(nadd == 0) continue;
        if(MergeShowersAndStore(fcnLabel, tjs, ci1, ci2, prt)) {
          Trajectory& stj = tjs.allTraj[ss1.ShowerTjID - 1];
          stj.AlgMod[kMergeNrShowers] = true;
          if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" success";
          break;
        }
      } // ci2
    } // ci1
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
     # 11 Debug in CTP (>10 debug cotIndex + 10)
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
          if(MergeShowersAndStore(fcnLabel, tjs, ict, jct, prt)) {
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
    
    std::vector<int> shList;
    std::vector<TrajPoint> tpList;
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      ShowerStruct& iss = tjs.cots[ict];
      if(iss.ID == 0) continue;
      if(iss.TjIDs.empty()) continue;
      if(iss.CTP != inCTP) continue;
      // Test an alternate cut
//      if(newCuts && iss.NearTjIDs.empty()) continue;
      // save the shower ID
      shList.push_back(iss.ID);
      // and the shower center TP
      tpList.push_back(tjs.allTraj[iss.ShowerTjID - 1].Pts[1]);
    } // ict
    if(shList.size() < 3) return;
    
    // sort by wire so the chain order is reasonable
    std::vector<SortEntry> sortVec(shList.size());
    for(unsigned short ii = 0; ii < shList.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].length = tpList[ii].Pos[0];
    }
    std::sort(sortVec.begin(), sortVec.end(), lessThan);
    auto tshList = shList;
    auto ttpList = tpList;
    for(unsigned short ii = 0; ii < shList.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      shList[ii] = tshList[indx];
      tpList[ii] = ttpList[indx];
    }
    
    // TODO: These cuts should be generalized somehow
    float minSep = 150;
    float maxDelta = 30;
    for(unsigned short ii = 0; ii < shList.size() - 2; ++ii) {
      auto& iss = tjs.cots[shList[ii] - 1];
      if(iss.ID == 0) continue;
      unsigned short jj = ii + 1;
      auto& jss = tjs.cots[shList[jj] - 1];
      if(jss.ID == 0) continue;
      std::vector<int> chain;
      float sepij = PosSep(tpList[ii].Pos, tpList[jj].Pos);
      if(sepij > minSep) continue;
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" sep ii "<<ii<<" "<<PrintPos(tjs, tpList[ii].Pos)<<" jj "<<jj<<" "<<PrintPos(tjs, tpList[jj].Pos)<<" sepij "<<sepij;
      // draw a line between these points
      TrajPoint tp;
      MakeBareTrajPoint(tjs, tpList[ii], tpList[jj], tp);
//      PrintTrajPoint("ij", tjs, 0, 1, 0, tp);
      for(unsigned short kk = jj + 1; kk < shList.size(); ++kk) {
        auto& kss = tjs.cots[shList[kk] - 1];
        if(kss.ID == 0) continue;
        float sepjk = PosSep(tpList[jj].Pos, tpList[kk].Pos);
        float delta = PointTrajDOCA(tjs, tpList[kk].Pos[0], tpList[kk].Pos[1], tp);
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<fcnLabel<<"   kk "<<kk<<" "<<PrintPos(tjs, tpList[kk].Pos)<<" sepjk "<<sepjk<<" delta "<<delta;
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
            chain[0] = shList[ii]; chain[1] = shList[jj]; chain[2] = shList[kk]; 
          } else {
            chain.push_back(shList[kk]);
          }
          // Refine the TP position and direction
          MakeBareTrajPoint(tjs, tpList[ii], tpList[kk], tp);
//          PrintTrajPoint("ik", tjs, 0, 0, chain.size(), tp);
        } // add to an existing chain
      } // kk
      // push the last one
      if(chain.size() > 2) {
        int newID = MergeShowers(fcnLabel, tjs, chain, prt);
        if(newID > 0 && AddTjsInsideEnvelope(fcnLabel, tjs, newID - 1, prt)) DefineShower(fcnLabel, tjs, newID - 1, prt);
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
    
    // Require that the maximum separation is about two radiation lengths
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" MergeSubShowers checking using radiation length cut ";
    
    constexpr float radLen = 14 / 0.3;
    
    for(unsigned short ict = 0; ict < tjs.cots.size() - 1; ++ict) {
      ShowerStruct& iss = tjs.cots[ict];
      if(iss.ID == 0) continue;
      if(iss.TjIDs.empty()) continue;
      if(iss.CTP != inCTP) continue;
      TrajPoint& istp1 = tjs.allTraj[iss.ShowerTjID - 1].Pts[1];
      for(unsigned short jct = ict + 1; jct < tjs.cots.size(); ++jct) {
        ShowerStruct& jss = tjs.cots[jct];
        if(jss.ID == 0) continue;
        if(jss.TjIDs.empty()) continue;
        if(jss.CTP != iss.CTP) continue;
        TrajPoint& jstp1 = tjs.allTraj[jss.ShowerTjID - 1].Pts[1];
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
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Candidate "<<iss.ID<<" "<<jss.ID<<" separation "<<sep<<" radiation lengths "<<trad<<" delta "<<delta<<" dang "<<dang;
        if(trad > 3) continue;
        // There must be a correlation between dang and the energy of these showers...
        if(dang > 0.3) continue;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Merge them. Re-find shower center, etc";
        if(MergeShowersAndStore(fcnLabel, tjs, ict, jct, prt)) {
          Trajectory& stj = tjs.allTraj[iss.ShowerTjID - 1];
          stj.AlgMod[kMergeSubShowers] = true;
          break;
        }
      } // jct
    } // ict
    
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
    unsigned short cotIndex = Create2DShower(tjs, tjl);
    if(cotIndex == USHRT_MAX) return 0;
    
    // define the new shower
    if(!DefineShower(fcnLabel, tjs, cotIndex, prt)) {
      std::cout<<fcnLabel<<" DefineShower failed\n";
      MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
      return 0;
    }
    
    return tjs.cots[cotIndex].ID;
    
  } // MergeShowers
  
  ////////////////////////////////////////////////
  bool MergeShowersAndStore(std::string inFcnLabel, TjStuff& tjs, unsigned short icotIndex, unsigned short jcotIndex, bool prt)
  {
    // Merge showers using shower indices. The icotIndex shower is modified in-place.
    // The jcotIndex shower is declared obsolete. This function also re-defines the shower and
    // sets the Parent ID to 0.
    
    if(icotIndex > tjs.cots.size() - 1) return false;
    ShowerStruct& iss = tjs.cots[icotIndex];
    if(iss.ID == 0) return false;
    if(iss.TjIDs.empty()) return false;
    if(iss.ShowerTjID <= 0) return false;
    
    if(jcotIndex > tjs.cots.size() - 1) return false;
    ShowerStruct& jss = tjs.cots[jcotIndex];
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
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" killed ShowerTjs "<<iss.ShowerTjID<<" and "<<jss.ShowerTjID<<" new Tj "<<ktj.ID;
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
    bool success = DefineShower(fcnLabel, tjs, icotIndex, prt);
 
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
    unsigned short icotIndex = GetCotsIndex(tjs, itj.ID);
    if(icotIndex == USHRT_MAX) return false;
    ShowerStruct& iss = tjs.cots[icotIndex];
    if(iss.ID == 0) return false;
    if(iss.TjIDs.empty()) return false;
    unsigned short jcotIndex = GetCotsIndex(tjs, jtj.ID);
    if(jcotIndex == USHRT_MAX) return false;
    ShowerStruct& jss = tjs.cots[jcotIndex];
    if(jss.ID == 0) return false;
    if(jss.TjIDs.empty()) return false;
    
    return MergeShowersAndStore(fcnLabel, tjs, icotIndex, jcotIndex, prt);
    
  } // MergeShowerTjsAndStore

  ////////////////////////////////////////////////
  bool FindChargeCenter(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Finds the charge center using all sub-structure trajectories in the cot. All of the shower
    // charge is assigned to the second TP and the charge weighted position is put in stp1.HitPos
    // and stp1.Pos
    // The charge will later be distributed between TP0 - TP2.
    // The total charge is stored in  shower Tj AveChg.
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
      if(ss.ShPts[ii].Chg <= 0) {
        std::cout<<fcnLabel<<" Found point with no charge. This shouldn't happen\n";
        ss.ID = 0;
        return false;
      }
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
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" cotIndex "<<cotIndex<<" HitPos "<<(int)stp1.HitPos[0]<<":"<<(int)stp1.HitPos[1]/tjs.UnitsPerTick<<" stp1.Chg "<<(int)stp1.Chg<<" Energy "<<(int)ss.Energy<<" MeV";
    return true;
  } // FindChargeCenter

  ////////////////////////////////////////////////
  void FindAngle(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Find the angle of the shower using the position of all of the TPs
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" "<<cotIndex<<" Pos "<<ss.CTP<<":"<<PrintPos(tjs, stp1)<<" Intercept "<<(int)A<<" dang "<<dang<<" Angle "<<ss.Angle<<" Err "<<ss.AngleErr<<" npts fit "<<nptsFit;
    
  } // FindAngle
  
  ////////////////////////////////////////////////
  void FillRotPos(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Fills the RotPos vector and sorts the points along the shower axis. Note that the rotation is
    // done around stp1.Pos but the charge center is at stp1.HitPos. Pos and HitPos will be exactly the
    // same if there is no parent. The Pos position may be shifted slightly in FindExternalParent so that
    // the parent trajectory lies on the central axis of the shower. This is done so that the charge at the
    // start of the shower is calculated correctly using the parent trajectory points
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" cotIndex "<<cotIndex<<" Rotation origin "<<PrintPos(tjs, stp1.Pos)<<" Angle "<<std::setprecision(2)<<ss.Angle<<" AspectRatio "<<ss.AspectRatio<<" AngleErr "<<ss.AngleErr;
    
  } // FillRotPos
  
  ////////////////////////////////////////////////
  bool AnalyzeRotPos(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // The RotPos vector was filled and sorted by increasing distance along the shower axis in FillRotPos.
    // This function divides the RotPos points into 3 sections and puts the transverse rms width in the
    // three sections into the shower Tj TrajPoint DeltaRMS variable. It also calculates the charge and number of shower
    // points closest to each TrajPoint. It make some crude quality cuts and returns false if the shower
    // fails these cuts
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.ID == 0) return false;
    if(ss.TjIDs.empty()) return false;
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return false;
    
    std::string fcnLabel = inFcnLabel + ".ARP";
    
    for(auto& tp : stj.Pts) {
      tp.Chg = 0;
      tp.DeltaRMS = 0;
      tp.NTPsFit = 0;
      tp.HitPos = {0, 0};
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
      myprt<<fcnLabel<<" ss "<<ss.ID;
      myprt<<" DeltaRMS "<<std::fixed<<std::setprecision(2)<<stj.Pts[0].DeltaRMS<<" "<<stj.Pts[1].DeltaRMS<<" "<<stj.Pts[2].DeltaRMS;
      myprt<<" DirectionFOM "<<std::fixed<<std::setprecision(2)<<ss.DirectionFOM;
    }
    return true;

  } // AnalyzeRotPos
  
  ////////////////////////////////////////////////
  bool DefineShowerTj(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Defines the Shower Tj, calculates the shower aspect ratio, etc. This function
    // doesn't change the state of Parent
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    if(!AnalyzeRotPos(fcnLabel, tjs, cotIndex, prt)) return false;

    // startsNeg is true if this assumption is correct
    bool startsNeg = (stj.Pts[0].DeltaRMS < stj.Pts[2].DeltaRMS);
    
    // reverse the points vector so that the narrow end of the shower is near Pts.begin()
    if(!startsNeg) ReverseShower(fcnLabel, tjs, cotIndex, prt);

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
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Reverses the shower and the shower tj
    
    if(cotIndex > tjs.cots.size() - 1) return;
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    DefineEnvelope(fcnLabel, tjs, cotIndex, prt);
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Reversed shower. Shower angle = "<<ss.Angle;
  }
  ////////////////////////////////////////////////
  void RefineShowerTj(TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Checks the properties of Shower Tj and revises them if necessary. Returns true if the
    // shower needs to be updated
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Gracefully kills the shower and the associated shower Tj
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
  void TagInShowerTjs(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList, bool applyMinTjCuts)
  {
    // Tag Tjs as InShower if they have MCSMom < ShowerTag[0] and there are more than
    // ShowerTag[6] other Tjs with a separation < ShowerTag[1]. Returns a list of Tjs that meet this criteria.
    // The shower cuts that are applicable are done here if applyMinTjCuts is true, for example when this function is called
    // from RunTrajClusterAlg before 3D vertex and tj matching is done. This reduces the number of spurious vertices and 3D matches.
    // When called from FindShowers3D, applyMinTjCuts is set false and the cuts are applied after 2D showers are reconstructed.
    
    tjList.clear();
    
    if(tjs.ShowerTag[0] <= 0) return;
    
    if(tjs.allTraj.size() > 20000) {
//      std::cout<<"TagInShowerTjs: Crazy number of Tjs "<<tjs.allTraj.size()<<". No shower tagging. Events processed "<<tjs.EventsProcessed<<" \n";
      return;
    }
    
    // clear out old tags and make a list of Tjs to consider
    std::vector<int> tjids;
    short maxMCSMom = tjs.ShowerTag[1];
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      tj.AlgMod[kInShower] = false;
      tj.NNeighbors = 0;
      if(tj.AlgMod[kShowerTj]) continue;
      if(tj.Pts.size() < 3) continue;
      if(tj.Pts.size() > 4 && tj.MCSMom > maxMCSMom) continue;
      tjids.push_back(tj.ID);
    } // tj
    
    if(tjids.size() < 2) return;
    
    for(unsigned short it1 = 0; it1 < tjids.size() - 1; ++it1) {
      Trajectory& tj1 = tjs.allTraj[tjids[it1] - 1];
      if(tj1.CTP != inCTP) continue;
      float len1 = TrajLength(tj1);
      for(unsigned short it2 = it1 + 1; it2 < tjids.size(); ++it2) {
        Trajectory& tj2 = tjs.allTraj[tjids[it2] - 1];
        if(tj2.CTP != inCTP) continue;
        unsigned short ipt1, ipt2;
        float doca = tjs.ShowerTag[2];
        // Find the separation between Tjs without considering dead wires
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, doca, false);
        if(doca == tjs.ShowerTag[2]) continue;
        // make tighter cuts for user-defined short Tjs
        float len2 = TrajLength(tj2);
        if(len1 < len2) {
          if(len1 < doca) continue;
        } else {
          if(len2 < doca) continue;
        }
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
    
    // eliminate entries that fail ShowerTag[7]
    std::vector<std::vector<int>> newList;
    for(auto& tjl : tjList) if(tjl.size() >= tjs.ShowerTag[7]) newList.push_back(tjl);
    tjList = newList;

    MergeTjList(tjList);

    // mark them all as InShower Tjs
    unsigned short nsh = 0;
    for(auto& tjl : tjList) {
      if(applyMinTjCuts && tjl.size() < tjs.ShowerTag[7]) continue;
      nsh += tjl.size();
      for(auto& tjID : tjl) {
        auto& tj = tjs.allTraj[tjID - 1];
        tj.AlgMod[kInShower] = true;
        // unset flags
        tj.AlgMod[kSetDir] = false;
        for(unsigned short end = 0; end < 2; ++end) tj.StopFlag[end][kBragg] = false;
      } // tjid
    } // tjl
    if(tjs.ShowerTag[12] >= 0) mf::LogVerbatim("TC")<<"TagInShowerTjs tagged "<<nsh<<" InShower Tjs in CTP "<<inCTP;
    
    if(tjs.UseAlg[kKillInShowerVx]) {
      // make a list of 2D vertices
      std::vector<unsigned short> vxids;
      for(auto& tjl : tjList) {
        for(auto& tjID : tjl) {
          auto& tj = tjs.allTraj[tjID - 1];
          for(unsigned short end = 0; end < 2; ++end) {
            if(!tj.AlgMod[kInShower]) continue;
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
        // get a list of Tjs attached to this vertex
        auto vxtjs = GetVtxTjIDs(tjs, vx2);
        // count the number that are InShower
        unsigned short ninsh = 0;
        for(auto tjid : vxtjs) {
          auto& tj = tjs.allTraj[tjid - 1];
          if(tj.AlgMod[kInShower]) ++ninsh;
        } // tjid
        if(ninsh > 1) MakeVertexObsolete(tjs, vx2, true);
      } // vxid
    } // 
    
  } // TagInShowerTjs
  
  ////////////////////////////////////////////////
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Find Tjs that are near the shower but are not included in it
    if(cotIndex > tjs.cots.size() - 1) return;
    ShowerStruct& ss = tjs.cots[cotIndex];
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
  void DefineEnvelope(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
      myprt<<fcnLabel<<" "<<cotIndex<<" Envelope";
      for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      myprt<<" Area "<<(int)ss.EnvelopeArea;
      myprt<<" ChgDensity "<<ss.ChgDensity;
    }
    // This is the last function used to update a shower
    ss.NeedsUpdate = false;
  } // DefineEnvelope  
  
  ////////////////////////////////////////////////
  bool AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
   {
    // This function adds Tjs to the shower. It updates the shower parameters.
    
     if(cotIndex > tjs.cots.size() - 1) return false;
    
     ShowerStruct& ss = tjs.cots[cotIndex];
     if(ss.Envelope.empty()) return false;
     if(ss.ID == 0) return false;
     if(ss.TjIDs.empty()) return false;
     
     std::string fcnLabel = inFcnLabel + ".ATIE";

     if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" ss.ID "<<ss.ID;

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
         if(AddTj(fcnLabel, tjs, tj.ID, cotIndex, false, prt)) ++nadd;
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
       if(AddTj(fcnLabel, tjs, tj.ID, cotIndex, false, prt)) {
         ++nadd;
       } else {
         if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" AddTj failed to add Tj "<<tj.ID;
       }
    } // tj
    
    if(nadd > 0) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added "<<nadd<<" trajectories ";
      ss.NeedsUpdate = true;
      if(ss.ParentID == 0) DefineShower(fcnLabel, tjs, cotIndex, prt);
      return true;
    } else {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" No new trajectories added to envelope ";
      ss.NeedsUpdate = false;
      return false;
    }
        
  } // AddTjsInsideEnvelope
  
  ////////////////////////////////////////////////
  bool AddLooseHits(TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Add hits that are inside the envelope to the shower if they are loose, i.e. not
    // used by any trajectory. This function returns true if the set of hits is different than
    // the current set. The calling function should update the shower if this is the case.
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
  void FindStartChg(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Finds the charge at the start of a shower and puts it in AveChg of the first
    // point of the shower Tj. This is only done when there is no parent.
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    auto schg = StartChgVec(tjs, cotIndex, prt);
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
    
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" cotIndex "<<cotIndex<<" Starting charge "<<(int)stp0.AveChg<<" startPt  "<<startPt;
    
  } // FindStartChg
  
  ////////////////////////////////////////////////
  std::vector<float> StartChgVec(TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Returns a histogram vector of the charge in bins of 1 WSE unit at the start of the shower

    ShowerStruct& ss = tjs.cots[cotIndex];
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
  void DumpShowerPts(TjStuff& tjs, unsigned short cotIndex)
  {
    // Print the shower points to the screen. The user should probably pipe the output to a text file
    // then grep this file for the character string PTS which is piped to a text file which can then be
    // imported into Excel, etc
    // Finds the charge at the start of a shower
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
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
    
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      ShowerStruct& ss = tjs.cots[cotIndex];
      if(ss.ID == 0) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
      // enough Tjs?
      unsigned short ntjs = ss.TjIDs.size();
      bool killit = (ntjs < tjs.ShowerTag[7]);
      // Kill runt showers
      if(!killit) killit = (ss.Energy < tjs.ShowerTag[3]);
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" cotIndex "<<cotIndex<<" nTjs "<<ss.TjIDs.size()<<" nTjs "<<ss.TjIDs.size()<<" killit? "<<killit;
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
      if(killit) MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
      
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
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      auto& ss = tjs.cots[cotIndex];
      if(ss.ID == 0) continue;
      if(ss.ShowerTjID == 0) continue;
      // Tp 1 of stj will get all of the shower hits
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      if(!stj.Pts[1].Hits.empty()) {
        std::cout<<"TTjH: ShowerTj "<<stj.ID<<" already has "<<stj.Pts[1].Hits.size()<<" hits\n";
        continue;
      }
      // Note that UseHit is not used since the size is limited to 16
      for(auto& tjID : ss.TjIDs) {
        unsigned short itj = tjID - 1;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) {
          std::cout<<"TTjH: Coding error. Tj "<<tjID<<" is a ShowerTj but is in TjIDs\n";
          continue;
        }
        if(!tjs.allTraj[itj].AlgMod[kInShower]) {
          std::cout<<"TTjH: Coding error. Trying to transfer Tj "<<tjID<<" hits but it isn't an InShower Tj\n";
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
     } // cotIndex

    if(prt) mf::LogVerbatim("TC")<<"TTJH: success? "<<newShowers;
    return newShowers;
  } // TransferTjHits

  ////////////////////////////////////////////////
  unsigned short GetCotsIndex(TjStuff& tjs, unsigned short ShowerTjID)
  {
    for(unsigned short ii = 0; ii < tjs.cots.size(); ++ii) {
      if(ShowerTjID == tjs.cots[ii].ShowerTjID) return ii;
    } // iii
    return USHRT_MAX;
    
  } // GetCotsIndex

  ////////////////////////////////////////////////
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss)
  {
    if(ss.ID == 0) return 0;
    if(ss.TjIDs.empty()) return 0;
    if(ss.ShowerTjID == 0) return 0;
    
    // Conversion from shower charge to energy in MeV. 0.0143 comes from an eye-bal fit.
    // Divide by the expected shower containment of 90%. This needs to be calculated directly
    constexpr float fShMeVPerChg = 0.0143 / 0.9;
    
    const Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    return fShMeVPerChg * stj.AveChg;
    
  } // ShowerEnergy
  
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
    ss.ShowerTjID = stj.ID;
    // put it in TJ stuff. The rest of the info will be added later
    tjs.cots.push_back(ss);
    return tjs.cots.size() - 1;
    
  } // Create2DShower
  
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
    myprt<<someText<<"  ID   CTP  ParID TruParID Energy nTjs  dFOM AspRat  stj  vx0 __Pos0___   Chg dRMS __Pos1___   Chg dRMS __Pos2___   Chg dRMS Angle SS3ID PFPID\n";

    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      myprt<<someText<<std::fixed;
//      myprt<<std::setw(4)<<ict;
      myprt<<std::setw(4)<<ss.ID;
      myprt<<std::setw(6)<<ss.CTP;
      myprt<<std::setw(7)<<ss.ParentID;
      myprt<<std::setw(9)<<ss.TruParentID;
      myprt<<std::setw(7)<<(int)ss.Energy;
      myprt<<std::setw(5)<<ss.TjIDs.size();
      myprt<<std::setw(6)<<std::setprecision(2)<<ss.DirectionFOM;
      myprt<<std::setw(7)<<std::setprecision(2)<<ss.AspectRatio;
      const auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      myprt<<std::setw(5)<<stj.ID;
      myprt<<std::setw(5)<<stj.VtxID[0];
      for(auto& spt : stj.Pts) {
        myprt<<std::setw(10)<<PrintPos(tjs, spt.Pos);
        myprt<<std::setw(6)<<(int)spt.Chg;
//        myprt<<std::setw(5)<<spt.NTPsFit;
        myprt<<std::setw(5)<<std::setprecision(1)<<spt.DeltaRMS;
      } // spt
      myprt<<std::setw(6)<<std::setprecision(2)<<stj.Pts[1].Ang;
      myprt<<std::setw(6)<<ss.SS3ID;
      if(ss.SS3ID < tjs.showers.size()) {
        auto& ss3 = tjs.showers[ss.SS3ID - 1];
        if(ss3.PFPIndex < tjs.pfps.size()) {
          auto& pfp = tjs.pfps[ss3.PFPIndex];
          myprt<<std::setw(6)<<pfp.ID;
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
        myprt<<std::setw(4)<<ss.ID;
        myprt<<" Tjs";
        for(auto id : ss.TjIDs) myprt<<" "<<id;
        myprt<<"\n";
      } // ict
      
    }
    // Print the envelopes
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      myprt<<someText<<std::fixed;
      myprt<<std::setw(4)<<ss.ID;
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
      myprt<<std::setw(4)<<ss.ID;
      myprt<<" Nearby";
      for(auto id : ss.NearTjIDs) myprt<<" "<<id;
      myprt<<"\n";
    } // ict
  } // Print2DShowers

} // namespace tca
