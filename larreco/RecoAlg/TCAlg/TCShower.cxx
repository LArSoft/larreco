#include "larreco/RecoAlg/TCAlg/TCShower.h"

struct SortEntry{
  unsigned int index;
  float length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}


namespace tca {

  ////////////////////////////////////////////////
  bool Find3DShowerEndPoints(TjStuff& tjs, MatchStruct& ms)
  {
    // The MatchStruct ms represents a 3D shower under construction. The 3D start position
    // and direction were found in Find3DEndpoints in TrajClusterAlg.cxx but that code
    // won't work well for showers since the ends in 2D are not well defined
    
    if(ms.PDGCode != 1111) return false;
    if(ms.Count == 0) return false;
    
    bool prt = (tjs.ShowerTag[12] == 3);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Inside F3DSEP: tjIDs";
      for(auto tjID : ms.TjIDs) myprt<<" "<<tjID;
      myprt<<" start vtx ID "<<ms.Vx3ID[0];
    }
    
    // See if the start end points are consistent
    std::vector<TrajPoint> spts;
    for(auto tjID : ms.TjIDs) {
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
      ms.dEdx[0][DecodeCTP(stj.CTP).Plane] = stj.dEdx[0];
      TrajPoint stp = stj.Pts[0];
      if(prt) mf::LogVerbatim("TC")<<" start point "<<PrintPos(tjs, stp)<<" dir "<<std::fixed<<std::setprecision(2)<<stp.Dir[0]<<" "<<stp.Dir[1];
      spts.push_back(stp);
    } // tjID
    if(spts.size() < 2) {
      if(prt) mf::LogVerbatim("TC")<<" Couldn't find at least 2 2D showers with good AspectRatio and DirectionFOM";
      return false;
    }
    if(prt) mf::LogVerbatim("TC")<<" stps size "<<spts.size();
    TVector3 pos, dir;
    if(!TrajPoint3D(tjs, spts[0], spts[1], pos, dir)) {
      if(prt) mf::LogVerbatim("TC")<<"  TrajPoint3D failed. Maybe the shower direction is fubar";
      return false;
    }
    // set the start position using a 3D vertex if one exists
    if(ms.Vx3ID[0] > 0) {
      auto& vx3 = tjs.vtx3[ms.Vx3ID[0] - 1];
      ms.XYZ[0][0] = vx3.X; ms.XYZ[0][1] = vx3.Y; ms.XYZ[0][2] = vx3.Z; 
    } else {
      // 
      ms.XYZ[0][0] = pos[0]; ms.XYZ[0][1] = pos[1]; ms.XYZ[0][2] = pos[2]; 
    }
    ms.Dir[0] = dir;
    
    // Now find the end point using the longest 2D shower
    double maxlen = 0;
    unsigned int maxID = 0;
    for(auto tjID : ms.TjIDs) {
      auto& stj = tjs.allTraj[tjID - 1];
      float length = TrajLength(stj);
      if(length > maxlen) {
        maxlen = length;
        maxID = tjID;
      }
    } // tjID
    
    auto& longTj = tjs.allTraj[maxID - 1];
    geo::PlaneID planeID = DecodeCTP(longTj.CTP);
    ms.BestPlane = planeID.Plane;
    double angleToVert = tjs.geom->WireAngleToVertical(tjs.geom->View(planeID), planeID.TPC, planeID.Cryostat) - 0.5 * ::util::pi<>();
    double cosgamma = std::abs(std::sin(angleToVert) * ms.Dir[0].Y() + std::cos(angleToVert) * ms.Dir[0].Z());
    if(cosgamma == 0) return false;
    // convert maxlen from WSE units (1 wire spacing) to cm and find the 3D distance
    maxlen *= tjs.geom->WirePitch(planeID) / cosgamma;
    
    for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
      ms.XYZ[1][ixyz] = ms.XYZ[0][ixyz] + ms.Dir[0][ixyz] * maxlen;
    }
    // Set the end direction to the start direction
    ms.Dir[1] = ms.Dir[0];
    
    return true;
  } // Find3DShowerEndPoints

  ////////////////////////////////////////////////
  void Finish3DShowers(TjStuff& tjs)
  {
    // Finish defining the shower properties that were not done previously
    
    if(tjs.ShowerTag[0] < 0) return;
    
    bool prt = (tjs.ShowerTag[12] >= 0);
    
    // Transfer the InShower Tj hits to the shower Tj
    for(const geo::TPCID& tpcid: tjs.geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
      for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
        CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
        if(!TransferTjHits(tjs, inCTP, prt)) return;
      } // plane
    } // tpcid
    
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      // Find a 2D shower that is matched to the 3D shower
      unsigned short iss = 0;
      for(unsigned short iss = 0; iss < tjs.cots.size(); ++iss) if(tjs.cots[iss].SS3ID == ss3.ID) break;
      if(iss == tjs.cots.size()) {
        std::cout<<"F3DS: Failed to find a 2D shower matched to ss3 "<<ss3.ID<<"\n";
        continue;
      }
      // now look for the associated shower Tj in the PFP list
      auto& ss = tjs.cots[iss];
      unsigned short mvpi = MatchVecPFPIndex(tjs, ss.ShowerTjID);
      if(mvpi == USHRT_MAX) {
        std::cout<<"F3DS: Failed to find a Shower Tj matchVecIndex for ss3 "<<ss3.ID<<"\n";
        continue;
      }
      unsigned short imv = tjs.matchVecPFPList[mvpi];
      auto& ms = tjs.matchVec[imv];
      if(ms.PDGCode != 1111) {
        std::cout<<"F3DS: The matchVecIndex for ss3 "<<ss3.ID<<" doesn't have a shower PDGCode\n";
        continue;
      }
      ms.PDGCode = 1111;
      ss3.Energy.resize(tjs.NumPlanes);
      ss3.EnergyErr.resize(tjs.NumPlanes);
      ss3.MIPEnergy.resize(tjs.NumPlanes);
      ss3.MIPEnergyErr.resize(tjs.NumPlanes);
      ss3.dEdx.resize(tjs.NumPlanes);
      ss3.dEdxErr.resize(tjs.NumPlanes);
      // fill the start position
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) ss3.Pos[ixyz] = ms.XYZ[0][ixyz];
      // and direction
      ss3.Dir = ms.Dir[0];
      if(prt) mf::LogVerbatim("TC")<<" Shower start "<<ss3.Pos.X()<<" "<<ss3.Pos.Y()<<" "<<ss3.Pos.Z()<<" dir "<<ss3.Dir.X()<<" "<<ss3.Dir.Y()<<" "<<ss3.Dir.Z();
      ss3.DirErr = ms.DirErr[0];
      // Find the shower length.
      ss3.Len = 0;
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
        double dpos = ms.XYZ[1][ixyz] - ms.XYZ[0][ixyz];
        ss3.Len += dpos * dpos;
      }
      ss3.Len = sqrt(ss3.Len);
      // TODO Calculate the opening angle here 
      ss3.OpenAngle = 0.1;
      // Fill shower energy and hits
      for(auto tjID : ms.TjIDs) {
        unsigned short cotIndex = GetCotsIndex(tjs, tjID);
        if(cotIndex > tjs.cots.size() - 1) continue;
        auto& ss = tjs.cots[cotIndex];
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        unsigned short plane = DecodeCTP(ss.CTP).Plane;
        ss3.dEdx[plane] = stj.dEdx[0];
        // TODO: do this correctly
        ss3.dEdxErr[plane] = 0.3 * stj.dEdx[0];
        ss3.Energy[plane] = ss.Energy;
        ss3.EnergyErr[plane] = 0.3 * ss.Energy;
        // just divide by 2.3 MeV/cm
        ss3.MIPEnergy[plane] = ss.Energy / 2.3;
        ss3.MIPEnergyErr[plane] = 0.3 * ss3.MIPEnergy[plane];
        // associate the shower Tj hits with the 3D shower
        auto tHits = PutTrajHitsInVector(stj, kUsedHits);
        ss3.Hits.insert(ss3.Hits.end(), tHits.begin(), tHits.end());
      } // tjID
      ss3.BestPlane = ms.BestPlane;
    } // ss3
    
  } // FinishShowers
  
  ////////////////////////////////////////////////
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    // Find 2D showers using 3D-matched trajectories. This returns true if showers were found
    // which requires re-doing the 3D trajectory match
    
    if(tjs.ShowerTag[0] != 2) return false;
    
    bool prt = false;
    short dbgPlane = tjs.ShowerTag[12];
    CTP_t dbgCTP = UINT_MAX;
    if(dbgPlane >= 0 && dbgPlane <= tjs.NumPlanes) dbgCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, dbgPlane);
    
    std::string fcnLabel = "FS";

    if(tjs.matchVec.empty()) {
      if(prt) mf::LogVerbatim("TC")<<"FindShowers3D: Give up because matchVec is empty";
      return false;
    }
    
    geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
    // check for already-existing showers
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      for(auto& ss : tjs.cots) if(ss.CTP == inCTP) return false;
    }
    
    // rebuild the hit range references if necessary
    if(tpcid != tjs.TPCID && !FillWireHitRange(tjs, tpcid, false)) return false;
    
    // lists of Tj IDs in plane, (list1, list2, list3, ...)
    std::vector<std::vector<std::vector<int>>> bigList(tjs.NumPlanes);
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      std::vector<std::vector<int>> tjList;
      TagShowerTjs(fcnLabel, tjs, inCTP, tjList);
      if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, tjList, fcnLabel);
      if(tjList.empty()) continue;
      MergeTjList(tjList);
      if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, tjList, "MTJL");
      bigList[plane] = tjList;
    } // plane
    unsigned short nPlanesWithShowers = 0;
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) if(!bigList.empty()) ++nPlanesWithShowers;
    if(nPlanesWithShowers < 2) return false;

    // mark them all as InShower Tjs
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      for(auto& tjl : bigList[plane]) {
        for(auto& tjID : tjl) tjs.allTraj[tjID - 1].AlgMod[kInShower] = true;
      } // tjl
      if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, bigList[plane], "MISTJ");
    } // plane

    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      prt = false;
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
        DefineShower(fcnLabel, tjs, cotIndex, prt);

        // Find nearby Tjs that were not included because they had too-high
        // MCSMom, etc. This will be used to decide if showers should be merged
        //AddTjsInsideEnvelope(tjs, cotIndex, prt, 1);

        if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, cotIndex, "DS");
        // Find nearby Tjs that were not included because they had too-high
        // MCSMom, etc. This will be used to decide if showers should be merged
        AddTjsInsideEnvelope(fcnLabel, tjs, cotIndex,prt);
        if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, cotIndex, "ATIE");
        FindNearbyTjs(fcnLabel, tjs, cotIndex, prt);
        FindMatchingTjs(fcnLabel, tjs, cotIndex, prt);
      } // tjl
      // try to merge showers in this plane using the lists of nearby Tjs
      if(inCTP == UINT_MAX) continue;
      if(tjs.cots.empty()) continue;
      prt = (inCTP == dbgCTP || dbgPlane == 3);
      if(prt) Print2DShowers("tjl", tjs, inCTP, true);
      MergeOverlap(fcnLabel, tjs, inCTP, prt);
      if(prt) Print2DShowers("MO", tjs, inCTP, true);
      MergeSubShowers(fcnLabel, tjs, inCTP, prt);
      if(prt) Print2DShowers("MSCi", tjs, inCTP, true);
      MergeShowerChain(fcnLabel, tjs, inCTP, prt);
      if(prt) Print2DShowers("MSCo", tjs, inCTP, true);
      MergeNearby2DShowers(fcnLabel, tjs, inCTP, prt);

    } // plane
    if(tjs.cots.empty()) return false;
    
    prt = false;
    
    // Look for a 3D-matched parent
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      auto& ss = tjs.cots[cotIndex];
      if(ss.ID == 0) continue;
      if (tjs.SaveShowerTree) SaveTjInfo(tjs, ss.CTP, cotIndex, "M2DS");
      prt = (ss.CTP == dbgCTP || dbgPlane == 3);
      FindExternalParent("FS", tjs, cotIndex, prt);
      if (tjs.SaveShowerTree) SaveTjInfo(tjs, ss.CTP, cotIndex, "FEP");
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      if(prt) std::cout<<cotIndex<<" Pos "<<ss.CTP<<":"<<PrintPos(tjs, stj.Pts[1].Pos)<<" ParID "<<ss.ParentID<<" TruParID "<<ss.TruParentID<<"\n";
      if(ss.ParentID == 0) FindStartChg(fcnLabel, tjs, cotIndex, prt);
    } // cotIndex
    
    prt = (dbgPlane > 3);
    if(prt) Print2DShowers("FEP", tjs, USHRT_MAX, false);
    Match2DShowers(fcnLabel, tjs, tpcid, prt);
    if(prt) Print2DShowers("M2DSo", tjs, USHRT_MAX, false);
    
    for(unsigned short plane = 0; plane < TPC.Nplanes(); ++plane) {
      CTP_t inCTP = EncodeCTP(tpcid.Cryostat, tpcid.TPC, plane);
      prt = (plane == dbgPlane || dbgPlane == 3);
      CheckQuality(fcnLabel, tjs, inCTP, prt);
    } // plane
    
    unsigned short nNewShowers = 0;
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      auto& ss = tjs.cots[cotIndex];
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      geo::PlaneID planeID = DecodeCTP(ss.CTP);
      if(planeID.Cryostat != tpcid.Cryostat) continue;
      if(planeID.TPC != tpcid.TPC) continue;
     ++nNewShowers;
    } // cotIndex
    
    if(prt) Print2DShowers("FSo", tjs, USHRT_MAX, false);
/*
    if(tjs.ShowerTag[12] >= 0) {
      for(auto& tj : tjs.allTraj) {
        if(tj.AlgMod[kKilled]) continue;
        if(!tj.AlgMod[kShowerTj]) continue;
        PrintTrajectory("FSO", tjs, tj, USHRT_MAX);
      }
    } // print trajectories
*/
    // clobber ss3 MatchVecPFPIndex since it will not be valid after re-matching in 3D
    for(auto& ss3 : tjs.showers) ss3.MatchVecPFPIndex = USHRT_MAX;
    
    return (nNewShowers > 0);
    
  } // FindShowers3D
  
  ////////////////////////////////////////////////
  void Match2DShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt)
  {
    if(prt) mf::LogVerbatim("TC")<<"Inside M2DS";
    // Use the tjs.showers vector to hold interim results. Note that some references in the
    // Shower3D struct, e.g. MatchVecPFPIndex will be invalid after 3D matching is redone
    
    // Clear out any old
    // stuff (which shouldn't exist...)
    for(unsigned short iss = 0; iss < tjs.showers.size(); ++iss) {
      // clear out any showers in this tpcid
      if(tjs.showers[iss].TPCID == tpcid) tjs.showers.erase(tjs.showers.begin() + iss); 
    } // ss3
    
    // convert the 2D separation cut into a 3D cut (later...)
    double sepCut = tjs.ShowerTag[2];
    
    std::string fcnLabel = inFcnLabel + ".M2DS";
    
/*
     // temp for testing
     if(inCTP == 0 && !tjs.MCPartList.empty()) {
     // Print some info on the first primary particle (electron)
     const simb::MCParticle* part = tjs.MCPartList[0];
     std::cout<<"Primary E = "<<std::setprecision(2)<<part->E();
     TVector3 dir;
     dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
     if(dir.Mag() != 0) dir.SetMag(1);
     std::cout<<" dir "<<std::setprecision(2)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
     for(CTP_t ctp = 0; ctp < 3; ++ctp) {
     std::cout<<" CTP "<<ctp;
     // print the primary trajectory ID
     std::cout<<" primary TjID "<<tm.MCParticleStartTjID(0, ctp);
     TrajPoint tp;
     tp.CTP = ctp;
     tm.MakeTruTrajPoint(0, tp);
     std::cout<<" Tru start pos "<<PrintPos(tjs, tp.Pos)<<" Ang "<<tp.Ang<<" Projection in plane "<<tp.Delta;
     std::cout<<"\n";
     } // ctp
     } // temp testing
*/
    
    for(unsigned short ci = 0; ci < tjs.cots.size() - 1; ++ci) {
      ShowerStruct& iss = tjs.cots[ci];
      if(iss.ID == 0) continue;
      if(iss.TjIDs.empty()) continue;
      geo::PlaneID iplaneID = DecodeCTP(iss.CTP);
      if(iplaneID.Cryostat != tpcid.Cryostat) continue;
      if(iplaneID.TPC != tpcid.TPC) continue;
      Trajectory& istj = tjs.allTraj[iss.ShowerTjID - 1];
//      if(prt) mf::LogVerbatim("TC")<<"ci "<<ci<<" istj "<<istj.ID<<" Energy "<<(int)iss.Energy;
      for(unsigned short cj = ci + 1; cj < tjs.cots.size(); ++cj) {
        ShowerStruct& jss = tjs.cots[cj];
        if(jss.CTP == iss.CTP) continue;
        if(jss.ID == 0) continue;
        if(jss.TjIDs.empty()) continue;
        geo::PlaneID jplaneID = DecodeCTP(jss.CTP);
        if(jplaneID.Cryostat != tpcid.Cryostat) continue;
        if(jplaneID.TPC != tpcid.TPC) continue;
        Trajectory& jstj = tjs.allTraj[jss.ShowerTjID - 1];
        TVector3 pos, dir;
        // Use shower Tj point 0 which should yield the start point of the 3D shower
        if(!TrajPoint3D(tjs, istj.Pts[0], jstj.Pts[0], pos, dir)) continue;
/*
        if(prt) {
          mf::LogVerbatim myprt("TC");
          myprt<<" cj "<<cj<<" jstj "<<jstj.ID<<" Energy "<<(int)jss.Energy;
          myprt<<" 3D Pos "<<std::fixed<<std::setprecision(1)<<pos[0]<<" "<<pos[1]<<" "<<pos[2];
          myprt<<" dir "<<std::setprecision(2)<<dir[0]<<" "<<dir[1]<<" "<<dir[2];
        }
*/
        // see if this position is consistent with an existing 3D shower
        unsigned short bestIndex = USHRT_MAX;
        float bestSep = sepCut;
        for(unsigned short ss3Index = 0; ss3Index < tjs.showers.size(); ++ ss3Index) {
          auto& ss3 = tjs.showers[ss3Index];
          if(ss3.ID == 0) continue;
          // Use a temporary TVector3 to get the separation
          TVector3 tmp = pos - ss3.Pos;
          float sep = tmp.Mag();
          if(sep < bestSep) {
            bestIndex = ss3Index;
            bestSep = sep;
          } // sep < sepCut
        } // ss3Index
        // matched to an existing 3D shower
        if(bestIndex < tjs.showers.size()) {
          auto& ss3 = tjs.showers[bestIndex];
          // add the cots indices to the list
          if(std::find(ss3.CotIndices.begin(), ss3.CotIndices.end(), ci) == ss3.CotIndices.end()) ss3.CotIndices.push_back(ci);
          if(std::find(ss3.CotIndices.begin(), ss3.CotIndices.end(), cj) == ss3.CotIndices.end()) ss3.CotIndices.push_back(cj);
//          if(prt) mf::LogVerbatim("TC")<<" bestSep "<<bestSep<<" (cm). Update ss3 "<<ss3.ID;
          continue;
        } else {
          // make a 3D shower struct to hold this information
          ShowerStruct3D ss3;
          ss3.ID = tjs.showers.size() + 1;
          ss3.TPCID = tpcid;
          ss3.Pos = pos;
          // TODO: would be nice to have a position and direction error here
          ss3.Dir = dir;
          ss3.CotIndices.resize(2);
          ss3.CotIndices[0] = ci;
          ss3.CotIndices[1] = cj;
          ss3.MatchVecPFPIndex = USHRT_MAX;
          // don't fill or use the rest of the variables
          tjs.showers.push_back(ss3);
//          if(prt) mf::LogVerbatim("TC")<<" new ss3 "<<ss3.ID;
        }
      } // ci2
    } // ci1
    if(tjs.showers.empty()) return;
    
    // Try to set MatchVecPFPIndex
    for(auto& ss3 : tjs.showers) {
      for(auto cotIndex : ss3.CotIndices) {
        auto& ss = tjs.cots[cotIndex];
        if(ss.ParentID > 0) {
          auto& ptj = tjs.allTraj[ss.ParentID - 1];
          // look for this parent ID in matchVec
          auto mvi = MatchVecIndex(tjs, ptj.ID);
          ss.SS3ID = ss3.ID;
          if(mvi < tjs.matchVecPFPList.size()) {
            // a matched parent exists
            if(ss3.MatchVecPFPIndex == USHRT_MAX) {
              // first occurrence
              ss3.MatchVecPFPIndex = mvi;
            } else if(ss3.MatchVecPFPIndex != mvi) {
              std::cout<<" Found inconsistent parent Tjs in a matched 3D shower "<<ss3.ID<<". Put some code here\n";
              // Need to decide which is better, the shower match or the parent match
            } 
          } // a matched parent exists
        } // a parent exists
      } // cotIndex
    } // ss3
    
    // look for incompletely defined 3D vertices
    for(auto& ss3 : tjs.showers) {
      if(ss3.ID == 0) continue;
      if(ss3.MatchVecPFPIndex < tjs.matchVecPFPList.size()) continue;
      std::cout<<"Found incomplete 3D shower "<<ss3.ID;
      std::cout<<" 2D shower info:";
      for(auto cotIndex : ss3.CotIndices) {
        auto& ss = tjs.cots[cotIndex];
        std::cout<<" ss.ID "<<ss.ID;
        auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
        std::cout<<" Pos "<<ss.CTP<<":"<<PrintPos(tjs, stj.Pts[1].Pos);
      }
      std::cout<<"\n";
    } // ss3
    
    // look for missing 2D -> 3D matches
    for(auto& ss : tjs.cots) {
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      if(ss.SS3ID > 0) continue;
      auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      std::cout<<"Found incompletely matched 2D shower "<<ss.ID<<" at Pos "<<ss.CTP<<":"<<PrintPos(tjs, stj.Pts[1].Pos)<<"\n";
    } // ss
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" Found "<<tjs.showers.size()<<" 3D showers";
      for(auto& ss3 : tjs.showers) {
        myprt<<"\nss3 "<<ss3.ID<<" MatchVecPFPIndex "<<ss3.MatchVecPFPIndex;
        myprt<<" 2D shower info:";
        for(auto cotIndex : ss3.CotIndices) {
          myprt<<"\n";
          auto& ss = tjs.cots[cotIndex];
          myprt<<" ss.ID "<<ss.ID;
          myprt<<" TruParent "<<ss.TruParentID;
          myprt<<" Parent "<<ss.ParentID;
          if(ss.ParentID > 0) {
            auto& ptj = tjs.allTraj[ss.ParentID - 1];
            if(ptj.AlgMod[kMat3D]) myprt<<" -> 3D match ";
            auto mvi = MatchVecIndex(tjs, ptj.ID);
            myprt<<" matchVecIndex "<<mvi;
          }
          myprt<<" Energy "<<(int)ss.Energy;
          auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
          myprt<<" stj "<<stj.ID;
          if(stj.AlgMod[kMat3D]) myprt<<" with 3D match";
        } // ss3cot
      } // ss3
    } // prt
    
  } // Match2DShowers

  ////////////////////////////////////////////////
  void FindMatchingTjs(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // fill the vector of Tjs in other planes that are matched to the Tjs in the shower
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    ss.MatchedTjIDs.clear();
    
    std::vector<int> mtj;
    for(auto& tjID : ss.TjIDs) {
      unsigned short mvIndex = MatchVecIndex(tjs, tjID);
      if(mvIndex > tjs.matchVec.size() - 1) continue;
      for(auto& tjID : tjs.matchVec[mvIndex].TjIDs) {
        Trajectory& tj = tjs.allTraj[tjID - 1];
        // not in the same plane
        if(tj.CTP == ss.CTP) continue;
        if(std::find(mtj.begin(), mtj.end(), tj.ID) == mtj.end()) mtj.push_back(tj.ID);
      } // otjID
    } // tjID
    if(mtj.empty()) return;
    if(mtj.size() > 1) std::sort(ss.MatchedTjIDs.begin(), ss.MatchedTjIDs.end());
    ss.MatchedTjIDs = mtj;

  } // FindMatchingTjs

  ////////////////////////////////////////////////
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Construct clusters of trajectories (cots) which will become shower PFParticles
    
    
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
     # 11 min Score to preserve a vertex and its Tjs
     # 12 Debug in CTP (>10 debug cotIndex + 10)
     */
    
    if(tjs.ShowerTag[0] != 1) return;
    
    MCParticleListUtils tm{tjs};
    
    bool prt = false;
    // print only one shower?
    unsigned short prtShower = USHRT_MAX;
    if(tjs.ShowerTag[12] >= 0) {
      geo::PlaneID planeID = DecodeCTP(inCTP);
      CTP_t printCTP = EncodeCTP(planeID.Cryostat, planeID.TPC, std::nearbyint(tjs.ShowerTag[12]));
      prt = (printCTP == inCTP);
      if(printCTP > 2) prt = true;
      if(printCTP > 9) prtShower = printCTP - 10;
    }
    // save the requested print state in case it gets changed
    bool saveprt = prt;
    /*
    // temp for testing
    if(inCTP == 0 && !tjs.MCPartList.empty()) {
      // Print some info on the first primary particle (electron)
      const simb::MCParticle* part = tjs.MCPartList[0];
      std::cout<<"Primary E = "<<std::setprecision(2)<<part->E();
      TVector3 dir;
      dir[0] = part->Px(); dir[1] = part->Py(); dir[2] = part->Pz();
      if(dir.Mag() != 0) dir.SetMag(1);
      std::cout<<" dir "<<std::setprecision(2)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<"\n";
      for(CTP_t ctp = 0; ctp < 3; ++ctp) {
        std::cout<<" CTP "<<ctp;
        // print the primary trajectory ID
        std::cout<<" primary TjID "<<tm.MCParticleStartTjID(0, ctp);
        TrajPoint tp;
        tp.CTP = ctp;
        tm.MakeTruTrajPoint(0, tp);
        std::cout<<" Tru start pos "<<PrintPos(tjs, tp.Pos)<<" Ang "<<tp.Ang<<" Projection in plane "<<tp.Delta;
        std::cout<<"\n";
      } // ctp
    } // temp testing
    */

    std::vector<std::vector<int>> tjList;
    TagShowerTjs(tjs, inCTP, tjList);
    
    //if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, tjList, "TSTJ");
    
    if(prt) std::cout<<"Inside FindShowers inCTP "<<inCTP<<" tjList size "<<tjList.size()<<"\n";
    if(tjList.empty()) return;
    MergeTjList(tjList);

    //if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, tjList, "MTJL");

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"tjlist\n";
      for(auto& tjl : tjList) {
        if(tjl.empty()) continue;
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
        myprt<<"\n";
      } // tjl
    } // prt


    //MergeTjList2(tjs, tjList, prt);
    
    // remove Tjs that don't have enough neighbors = ShowerTag[7] unless the shower
    // has few Tjs
    unsigned short minNeighbors = tjs.ShowerTag[7];
    if(minNeighbors > 0) --minNeighbors;
    for(auto& tjl : tjList) {
      bool didErase = true;
      while(didErase) {
        didErase = false;
        unsigned short indx = 0;
        for(indx = 0; indx < tjl.size(); ++indx) {
          unsigned short itj = tjl[indx] - 1;
          if(tjl.size() > 5 && tjs.allTraj[itj].NNeighbors < minNeighbors) break;
        } // indx
        if(indx < tjl.size()) {
          tjl.erase(tjl.begin() + indx);
          didErase = true;
        }
      } // didErase
    } // tjl
    
    //if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, tjList, "CNN");

    // mark all of these as InShower Tjs
    for(auto& tjl : tjList) {
      for(auto& tjID : tjl) {
        tjs.allTraj[tjID - 1].AlgMod[kInShower] = true;
      } // tjID
    } // tjl

    // Convert each one into a shower with a shower Tj
    for(auto& tjl : tjList) {
      if(tjl.empty()) continue;
      // Create the shower Tj
      unsigned short cotIndex = Create2DShower(tjs, tjl);
      if(cotIndex == USHRT_MAX) continue;
      ShowerStruct& ss = tjs.cots[cotIndex];
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID-1];
      if(prt && prtShower != USHRT_MAX) prt = (cotIndex == prtShower);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Make cots "<<cotIndex<<" in CTP "<<ss.CTP<<" TjID_NN";
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
      }
      if(!DefineShower(tjs, cotIndex, prt)) continue;
      if(tjs.cots[cotIndex].TjIDs.empty()) continue;

      //if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, cotIndex, "DS");
      // Fill the vector of Tjs that are close to this shower but were not included in it, most
      // likely because the MCSMom is too high. These will be used to merge showers
      // FindNearbyTjs(tjs, cotIndex, prt);
      // Try to add more Tjs to the shower

      //AddTjsInsideEnvelope(tjs, cotIndex, prt,1);
      AddTjsInsideEnvelope(tjs, cotIndex, prt);
      //if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, cotIndex, "ATJIE");

      FindExternalParent(tjs, cotIndex, prt);

      //if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, cotIndex, "FEP");
      // If no external parent was found, try to refine the direction and look for
      // an internal parent
      if(tjs.cots[cotIndex].ShowerTjID == 0) RefineShowerTj(tjs, cotIndex, prt);
      if(prt) PrintTrajectory("FS", tjs, tjs.allTraj[stj.ID-1], USHRT_MAX);



    } // tjl
    
    if(tjs.cots.empty()) return;

    prt = saveprt;
    // Merge showers whose envelopes overlap
    MergeOverlap(tjs, inCTP, prt);
    // Merge small showers with larger ones
    MergeSubShowers(tjs, inCTP, prt);
    
    // drop those that don't meet the requirements
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      ShowerStruct& ss = tjs.cots[cotIndex];
      if(ss.CTP != inCTP) continue;
      if(ss.TjIDs.empty()) continue;
      // enough Tjs?
      unsigned short ntjs = ss.TjIDs.size();
      bool killit = (ntjs < tjs.ShowerTag[7]);
      // Kill runt showers
      if(prt && prtShower != USHRT_MAX) prt = (cotIndex == prtShower);
      if(!killit) killit = (ss.Energy < tjs.ShowerTag[3]);
      if(prt) mf::LogVerbatim("TC")<<"cotIndex "<<cotIndex<<" nTjs "<<ss.TjIDs.size()<<" nTjs "<<ss.TjIDs.size()<<" killit? "<<killit;
      if(!killit) {
        // count the number of Tj points
        unsigned short nTjPts = 0;
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          nTjPts += NumPtsWithCharge(tjs, tj, false);
        }  // tjID
        if(nTjPts < tjs.ShowerTag[6]) killit = true;
        if(prt) mf::LogVerbatim("TC")<<"    "<<" nTjPts "<<nTjPts<<" killit? "<<killit;
      } // !killit
      if(killit) {
        MakeShowerObsolete(tjs, cotIndex, prt);
      } else {
        if(tjs.allTraj[ss.ShowerTjID - 1].AlgMod[kKilled]) {
          std::cout<<"FS logic error: ShowerTj "<<tjs.allTraj[ss.ShowerTjID - 1].ID<<" is killed\n";
          tjs.allTraj[ss.ShowerTjID - 1].AlgMod[kKilled] = false;
        }
        // A good shower. Set the pdgcode of InShower Tjs to 11
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          tj.PDGCode = 11;
          // Clobber 2D vertices that are inside the shower
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] > 0) {
              VtxStore& vx2 = tjs.vtx[tj.VtxID[end]-1];
              bool killMe = (vx2.Score < tjs.ShowerTag[11]);
              // don't kill the 2D vertex if it is attached to the far end of the parent Tj
              if(killMe && ss.ParentID == tjID) {
                unsigned short farEnd = FarEnd(tjs, tj, ss);
                if(farEnd == end) killMe = false;
              }
              if(killMe) {
                if(prt) mf::LogVerbatim("TC")<<"Clobber vtx "<<tj.VtxID[end]<<" Score "<<vx2.Score<<" Vtx3ID "<<tjs.vtx[tj.VtxID[end]-1].Vtx3ID;
                MakeVertexObsolete(tjs, tj.VtxID[end], true);
              }
            }
          } // end
        }
      } // don't killit

    } // ic

    CheckQuality(tjs, inCTP, prt);
    
    // find the start charge
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      ShowerStruct& ss = tjs.cots[cotIndex];
      if(ss.CTP != inCTP) continue;
      if(ss.TjIDs.empty()) continue;
      FindStartChg(tjs, cotIndex, prt);
    }
    
    // Finish up in this CTP. 
    // Re-assign hits from the InShower Tjs to the ShowerTj.
    TransferTjHits(tjs, inCTP, prt);
//    std::cout<<"Final calculation shower energy...\n";

    // check for consistency
    unsigned short cotIndex = 0; 
    for(auto& ss : tjs.cots) {
      if(ss.TjIDs.empty()) continue;
      if(ss.CTP != inCTP) continue;
      if(ss.ShowerTjID == 0) {
        std::cout<<"FindShowers: ShowerTjID not defined in CTP "<<ss.CTP<<"\n";
      }
      for(auto& tjID : ss.TjIDs) {
        if(tjID > (int)tjs.allTraj.size()) {
          std::cout<<"FindShowers: Bad tjID "<<tjID<<"\n";
        }
        Trajectory& tj = tjs.allTraj[tjID - 1];
        if(tj.CTP != ss.CTP) {
          std::cout<<"FindShowers: Bad CTP "<<ss.CTP<<" "<<tj.CTP<<" tjID "<<tjID<<"\n";
        }
        if(!tj.AlgMod[kKilled] || !tj.AlgMod[kInShower]) {
          std::cout<<"FindShowers: InShower TjID "<<tjID<<" invalid kKilled "<<tj.AlgMod[kKilled]<<" or kInShower "<<tj.AlgMod[kInShower]<<"\n";
        }
      } // tjID

      cotIndex++;
    } // ss
    
    if(tjs.ShowerTag[12] >= 0) {
      for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
        if(tjs.cots[ic].TjIDs.empty()) continue;
        unsigned short itj = tjs.cots[ic].ShowerTjID - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(prt || (tjs.ShowerTag[12] == 3 && tj.CTP == inCTP)) PrintTrajectory("FSO", tjs, tj, USHRT_MAX);
      } // ic
    } // print trajectories
    
    if(tjs.ShowerTag[12] >= 100) {
      unsigned short ic = tjs.ShowerTag[12] - 100;
      if(ic < tjs.cots.size() && tjs.cots[ic].CTP == inCTP) DumpShowerPts(tjs, ic);
    }
    
  } // FindShowers
>>>>>>> feature/XL_TJWork
  
  ////////////////////////////////////////////////
  void MergeTjList(std::vector<std::vector<int>>& tjList)
  {
    // Merge the lists of Tjs in the lists if they share a common Tj ID
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
            ss.ShPts[cnt].Pos[0] = tjs.fHits[iht].WireID.Wire;
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
    
    std::string fcnLabel = inFcnLabel + ".DS";
    
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
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    
    bool addParent = (ss.ParentID == tjID);
    
    std::string fcnLabel = inFcnLabel + ".ATj";

    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[tjID - 1];
    if(tj.AlgMod[kInShower]) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is already an InShower Tj";
      return false;
    }
    if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tjID) != ss.TjIDs.end()) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is already in this shower "<<cotIndex;
      return false;
    }
    if(tj.CTP != ss.CTP) {
      mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" is in the wrong CTP "<<cotIndex;
      return false;
    }
    // don't add a Tj to the shower if it has a nice vertex unless it is the parent
    if(!addParent && tj.AlgMod[kTjHiVx3Score]) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" in cotIndex "<<cotIndex<<" has a high score 3D vertex and is not the ParentID "<<ss.ParentID<<". Not adding it";
      return false;
    }
    // check for high score vertex
    if(!addParent && TjHasNiceVtx(tjs, tj, tjs.ShowerTag[11])) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Tj "<<tjID<<" in cotIndex "<<cotIndex<<" has a high score 2D vertex. Not adding it";
      return false;
    }
    // Ensure that we are not trying to add a Tj that is attached to the shower Tj vertex
    if(stj.VtxID[0] > 0) {
      // Get a list of TjIDs that are attached to the shower Tj vertex
      auto vxTjIDs = GetVtxTjIDs(tjs, tjs.vtx[stj.VtxID[0] - 1]);
      // make sure this Tj ID isn't in the list
      if(std::find(vxTjIDs.begin(), vxTjIDs.end(), tjID) != vxTjIDs.end()) return false;
    }
    ss.TjIDs.push_back(tjID);
    std::sort(ss.TjIDs.begin(), ss.TjIDs.end());
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
          ss.ShPts[cnt].Pos[0] = tjs.fHits[iht].WireID.Wire;
          ss.ShPts[cnt].Pos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
          ++cnt;
        }
      }
    } // ipt
    tj.AlgMod[kInShower] = true;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Add Tj "<<tj.ID;
    
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
    if(ss.Envelope.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".FEP";
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" FindExternalParent ss. ID "<<ss.ID;
    
    // require 3D matched parent?
    bool require3DMatch = (tjs.ShowerTag[0] == 2);
    
    // References to shower Tj points
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    TrajPoint& stp0 = stj.Pts[0];
    TrajPoint& stp1 = stj.Pts[1];
    int oldParent = ss.ParentID;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" cotIndex "<<cotIndex<<" stj.ID "<<stj.ID<<" Existing parent ID "<<oldParent<<" parent FOM "<<ss.ParentFOM;
      myprt<<" attached to vertex "<<stj.VtxID[0];
      myprt<<" Tjs";
      for(auto& tid : ss.TjIDs) myprt<<" "<<tid;
    }
    
    if(ss.AspectRatio > tjs.ShowerTag[10] || ss.DirectionFOM > tjs.ShowerTag[9]) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Don't search for a parent due to poor AspectRatio "<<ss.AspectRatio<<" or ss.DirectionFOM "<<ss.DirectionFOM;
      return;
    }
    
    float bestFOM = ss.ParentFOM;
    int imTheBest = 0;
    unsigned short imTheBestPt = 0;
    unsigned short imTheBestEnd = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled] && !tj.AlgMod[kInShower]) continue;
      // ignore shower Tjs. Note that this also rejects parent Tjs of other showers
      if(tj.AlgMod[kShowerTj]) continue;
      // ignore in-shower Tjs that aren't in this shower
      if(tj.AlgMod[kInShower] && std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) == ss.TjIDs.end()) continue;
      // ignore existing shower parents
      if(tj.AlgMod[kShwrParent]) continue;
      // Ignore short Tjs
      if(tj.Pts.size() < 5) continue;
      
      if(require3DMatch && !tj.AlgMod[kMat3D]) continue;

      unsigned short useEnd = 0;
      // Check trajectories that were split by 3D vertex matching
      if(WrongSplitTj(tjs, tj, useEnd, ss, prt)) continue;
      float fom = ParentFOM(fcnLabel, tjs, tj, useEnd, ss, prt);
      if(fom > bestFOM) continue;
      bestFOM = fom;
      imTheBest = tj.ID;
      imTheBestPt = tj.EndPt[useEnd];
      imTheBestEnd = useEnd;
    } // tj

    if(imTheBest < 1 || imTheBest > (int)tjs.allTraj.size()) return;
    
    if(bestFOM > tjs.ShowerTag[8]) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Best parent candidate FOM "<<bestFOM<<" exceeds the cut "<<tjs.ShowerTag[8];
      // Remove an old parent?
      if(oldParent > 0) {
        Trajectory& oldParentTj = tjs.allTraj[oldParent-1];
        // remove the parent flag
        oldParentTj.AlgMod[kShwrParent] = false;
        // remove it from the shower and update if it attached to a high-score 3D vertex
        if(oldParentTj.AlgMod[kTjHiVx3Score]) RemoveTj(fcnLabel, tjs, oldParent, cotIndex, true, prt);
      }
      ss.ParentID = 0;
      // detach the showerTj from a vertex if one exists
      stj.VtxID[0] = 0;
      return;
    }
    
    ss.ParentID = imTheBest;
    ss.ParentFOM = bestFOM;
    Trajectory& parentTj = tjs.allTraj[imTheBest - 1];
    parentTj.AlgMod[kShwrParent] = true;
    // Move the shower start point to the parent
    TrajPoint& ptp = parentTj.Pts[imTheBestPt];
    stp0.Pos = ptp.Pos;
    // Set the shower angle to the parent angle. Ensure that they are in the same direction. Assume that
    // 3D-matched Tjs are in the correct direction
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Compare parent angle "<<ptp.Ang<<" with shower angle "<<ss.Angle;
    if(!require3DMatch && std::abs(ptp.Ang - ss.Angle) > M_PI/2) {
      ReverseShower(fcnLabel, tjs, cotIndex, prt);
    } // angle consistency check
    ss.Angle = ptp.Ang;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<"   shower angle "<<ss.Angle;
    double cs = cos(ptp.Ang);
    double sn = sin(ptp.Ang);
    // move all the points onto the parent trajectory
    for(auto& tp : stj.Pts) {
      double dist = PosSep(ptp.Pos, tp.Pos);
      tp.Pos[0] = ptp.Pos[0] + cs * dist;
      tp.Pos[1] = ptp.Pos[1] + sn * dist;
    } // tp
    
    // we re-found the old parent
    if(ss.ParentID == oldParent) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Existing parent is good ";
      ss.NeedsUpdate = false;
      return;
    }
    
    // determine if the new parent is external or internal to the shower
    if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), parentTj.ID) == ss.TjIDs.end()) {
      // Parent is external
      // Remove the old parent from the shower if it has a nice vertex. Don't update since that will be done in a bit
      if(oldParent > 0) {
        Trajectory& oldParentTj = tjs.allTraj[oldParent-1];
        // remove the parent flag
        oldParentTj.AlgMod[kShwrParent] = false;
        // remove it from the shower if it attached to a high-score 3D vertex
        if(oldParentTj.AlgMod[kTjHiVx3Score]) RemoveTj(fcnLabel, tjs, oldParent, cotIndex, false, prt);
      } // oldParent exists
      // Add the Tj to the shower but don't bother defining it yet
      if(!AddTj(fcnLabel, tjs, imTheBest, cotIndex, false, prt)) {
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Failed to add the parent for some reason. Recovering...";
        ss.ParentID = 0;
        ss.ParentFOM = 999;
        parentTj.AlgMod[kShwrParent] = false;
        return;
      }
      // attach the shower Tj to a vertex (if one exists)
      stj.VtxID[0] = parentTj.VtxID[imTheBestEnd];
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added new parent "<<imTheBest<<" to the shower with FOM = "<<bestFOM<<" stp1 Pos "<<PrintPos(tjs, stp1)<<" vtxID "<<stj.VtxID[0];
    } else {
      // Parent is internal
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Existing InShower Tj "<<imTheBest<<" promoted to parent status with FOM = "<<bestFOM<<" stp1 Pos "<<PrintPos(tjs, stp1);
    } // end of external/internal decision
    ss.NeedsUpdate = true;
    // Update the shower with this new information
    UpdateShowerWithParent(fcnLabel, tjs, cotIndex, prt);
    
  } // FindExternalParent
  
  ////////////////////////////////////////////////
  void UpdateShowerWithParent(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // This updates all shower and shower Tj parameters when a new shower parent Tj is identified.
    if(cotIndex > tjs.cots.size() - 1) return;
    ShowerStruct& ss = tjs.cots[cotIndex];
    // Ensure that everything is valid
    if(ss.ID == 0) return;
    if(ss.TjIDs.empty()) return;
    if(ss.ParentID == 0) return;
    if(ss.ShowerTjID == 0) return;
    if(!ss.NeedsUpdate) return;
    
    std::string fcnLabel = inFcnLabel + ".USWP";
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" UpdateShowerWithParent ssID  "<<ss.ID;
    
    // The ParentID will be set to 0 if any failure occurs after this point
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    Trajectory& ptj = tjs.allTraj[ss.ParentID - 1];
    // determine the appropriate start end of the parent
    unsigned short pend = FarEnd(tjs, ptj, ss);
    if(pend != 0) {
      std::cout<<fcnLabel<<" Parent end is not 0... Is this bad?\n";
    }
    
    // set dE/dx of the shower Tj
    stj.dEdx[0] = ptj.dEdx[pend];
    std::cout << "USWP: ParentID " << ss.ParentID << " dEdx " << stj.dEdx[0] << std::endl;
    // Clear
    for(auto& stp : stj.Pts) {
      stp.Chg = 0;
      stp.Hits.clear();
    }
    
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
            ss.ShPts[cnt].HitIndex = iht;
            ss.ShPts[cnt].TID = tj.ID;
            ss.ShPts[cnt].Chg = tjs.fHits[iht].Integral;
            ss.ShPts[cnt].Pos[0] = tjs.fHits[iht].WireID.Wire;
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

    // This function calculates the charge (Chg), the charge center (HitPos) and
    // the transverse shower rms at each Tp in the shower Tj
    AnalyzeRotPos(fcnLabel, tjs, cotIndex, prt);
    
    float minAlong = ss.ShPts[0].RotPos[0];
    float maxAlong = ss.ShPts[ss.ShPts.size()-1].RotPos[0];
    
    // put the shower Tj start position at the parent Tj start
    TrajPoint& stp0 = stj.Pts[0];
    stp0.Pos = ptj.Pts[pend].Pos;
    
    // rotate the other points
    TrajPoint& stp1 = stj.Pts[1];
    stp1.Pos[0] = stp0.Pos[0] + stp1.Dir[0] * (0 - minAlong);
    stp1.Pos[1] = stp0.Pos[1] + stp1.Dir[1] * (0 - minAlong);
    
    TrajPoint& stp2 = stj.Pts[2];
    stp2.Pos[0] = stp0.Pos[0] + stp1.Dir[0] * (maxAlong - minAlong);
    stp2.Pos[1] = stp0.Pos[1] + stp1.Dir[1] * (maxAlong - minAlong);
    
    FindNearbyTjs(fcnLabel, tjs, cotIndex, prt);
    DefineEnvelope(fcnLabel, tjs, cotIndex, prt);
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
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, bool prt)
  {
    // returns a FOM for the trajectory at the end point being the parent of ss and the end which
    // was matched
    
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
    float tp1Sep = PosSep2(ptp.Pos, stp1.Pos);
    // Make a rough cut on radiation lengths
    if(tp1Sep > tenRadLen2) {
//      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" tjID "<<tj.ID<<" failed sep cut "<<(int)sqrt(tp1Sep)<<" 10 radiation lengths "<<(int)sqrt(tenRadLen2);
      return 100;
    }
    tp1Sep = sqrt(tp1Sep);
    
    // impact parameter between the projection of ptp and the charge center
    float delta = PointTrajDOCA(tjs, stp1.HitPos[0], stp1.HitPos[1], ptp);
    // make a rough cut
    if(delta > 100) {
//      if(prt) mf::LogVerbatim("TC")<<"PFOM "<<tj.ID<<" failed delta cut "<<delta<<" cut = 100";
      return 50;
    }
    
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
    float vx3Score = 0;
    if(tj.AlgMod[kMat3D] && tj.VtxID[tjEnd] > 0) {
      // check for a 3D vertex at this end
      VtxStore& vx2 = tjs.vtx[tj.VtxID[tjEnd] - 1];
      if(vx2.Vtx3ID > 0 && vx2.Vtx3ID < tjs.vtx3.size() && vx2.Stat[kHiVx3Score]) {
        vx3Score = tjs.vtx3[vx2.Vtx3ID - 1].Score;
        // prevent nuttiness
        fom /= 2;
      }
    }
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel;
      myprt<<" ssID "<<ss.ID;
      myprt<<" Tj "<<tj.ID<<" Pos "<<PrintPos(tjs, ptp);
      myprt<<" VtxID "<<tj.VtxID[tjEnd];
      if(tj.VtxID[tjEnd] > 0) {
        VtxStore& vx2 = tjs.vtx[tj.VtxID[tjEnd-1]];
        if(vx2.Vtx3ID > 0) myprt<<" Vtx3ID "<<vx2.Vtx3ID;
      }
      myprt<<std::fixed<<std::setprecision(2);
      myprt<<" tp1Sep "<<std::fixed<<std::setprecision(1)<<tp1Sep<<" pull "<<sepPull;
      myprt<<" delta "<<delta<<" pull "<<deltaPull;
      myprt<<" dang "<<dang<<" pull "<<dangPull;
      myprt<<" mcsmom "<<(int)mom<<" pull "<<momPull;
      myprt<<" sep0 "<<sqrt(tp0Sep2)<<" pull "<<sqrt(sep0Pull2);
      myprt<<" length "<<(int)TrajLength(tj)<<" pull "<<lenPull;
      myprt<<" vx3Score "<<vx3Score;
      myprt<<" FOM "<<fom;
    }

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

    return fom;
  } // ParentFOM
  
  ////////////////////////////////////////////////
  bool WrongSplitTj(TjStuff& tjs, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt)
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
    
    // See if the other end is the end that was split. It should have a vertex with Topo = 8 or 11
    unsigned short otherEnd = 1 - tjEnd;
//    if(prt) mf::LogVerbatim("TC")<<"WSTj: otherEnd "<<otherEnd<<" vtxID "<<tj.VtxID[otherEnd];
    if(tj.VtxID[otherEnd] == 0) return false;
    unsigned short ivx = tj.VtxID[otherEnd] - 1;
    // A vertex exists but not a 3D split vertex
    if(tjs.vtx[ivx].Topo != 8 && tjs.vtx[ivx].Topo != 10) return false;
    return true;
    
  } // WrongSplitTj

  ////////////////////////////////////////////////
  void MergeNearby2DShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    if(!tjs.UseAlg[kMergeNrShowers]) return;
    if(tjs.cots.empty()) return;
    
    std::string fcnLabel = inFcnLabel + ".MNS";
    
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
        // add the shared Tjs to ss1
        // ensure that the shower isn't InShower already
        unsigned short nadd = 0;
        for(auto& tjID : shared) {
          if(tjs.allTraj[tjID - 1].AlgMod[kInShower]) continue;
          // don't put it in the shower if it has a nice Tj
          if(tjs.allTraj[tjID - 1].AlgMod[kTjHiVx3Score]) continue;
          if(AddTj(fcnLabel, tjs, tjID, ci1, false, prt)) ++nadd;
        } // tjID
        if(nadd == 0) continue;
        if(MergeShowersAndStore(fcnLabel, tjs, ci1, ci2, prt)) {
          Trajectory& stj = tjs.allTraj[ss1.ShowerTjID - 1];
          stj.AlgMod[kMergeNrShowers] = true;
          if(prt) mf::LogVerbatim("TC")<<"MN2DS: Merged ss2 "<<ci2<<" into "<<ci1;
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
    
    // See if the envelopes overlap
//    unsigned short maxict = tjs.cots.size() - 1;
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
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Merge them. Re-find shower center, etc";
        if(MergeShowersAndStore(fcnLabel, tjs, ict, jct, prt)) {
          Trajectory& stj = tjs.allTraj[iss.ShowerTjID - 1];
          stj.AlgMod[kMergeOverlap] = true;
          break;
        }
      } // jct
    } // ict

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
    
    float minSep = 150;
    float maxDelta = 20;
    for(unsigned short ii = 0; ii < shList.size() - 2; ++ii) {
      auto& iss = tjs.cots[shList[ii] - 1];
      if(iss.ID == 0) continue;
      unsigned short jj = ii + 1;
      auto& jss = tjs.cots[shList[jj] - 1];
      if(jss.ID == 0) continue;
      std::vector<int> chain;
      float sepij = PosSep(tpList[ii].Pos, tpList[jj].Pos);
      if(sepij > minSep) continue;
//      std::cout<<" sep ii "<<ii<<" "<<PrintPos(tjs, tpList[ii].Pos)<<" jj "<<jj<<" "<<PrintPos(tjs, tpList[jj].Pos)<<" sepij "<<sepij<<"\n";
      // draw a line between these points
      TrajPoint tp;
      MakeBareTrajPoint(tjs, tpList[ii], tpList[jj], tp);
      for(unsigned short kk = jj + 1; kk < shList.size(); ++kk) {
        auto& kss = tjs.cots[shList[kk] - 1];
        if(kss.ID == 0) continue;
        float sepjk = PosSep(tpList[jj].Pos, tpList[jj].Pos);
        float delta = PointTrajDOCA(tjs, tpList[kk].Pos[0], tpList[kk].Pos[1], tp);
//        std::cout<<"   kk "<<kk<<" "<<PrintPos(tjs, tpList[kk].Pos)<<" sepjk "<<sepjk<<"\n";
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
        } // add to an existing chain
      } // kk
      // push the last one
      if(chain.size() > 2) {
        int newID = MergeShowers(fcnLabel, tjs, chain, prt);
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
    
    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      ShowerStruct& iss = tjs.cots[ict];
      if(iss.ID == 0) continue;
      if(iss.TjIDs.empty()) continue;
      if(iss.CTP != inCTP) continue;
      TrajPoint& istp0 = tjs.allTraj[iss.ShowerTjID - 1].Pts[0];
      TrajPoint& istp2 = tjs.allTraj[iss.ShowerTjID - 1].Pts[2];
      for(unsigned short jct = 0; jct < tjs.cots.size(); ++jct) {
        if(jct == ict) continue;
        ShowerStruct& jss = tjs.cots[jct];
        if(jss.ID == 0) continue;
        if(jss.TjIDs.empty()) continue;
        if(jss.CTP != iss.CTP) continue;
        // require that the j shower be lower energy than the i shower
        if(jss.Energy > iss.Energy) continue;
        // require that it be downstream of the i shower
        TrajPoint& jstp0 = tjs.allTraj[jss.ShowerTjID - 1].Pts[0];
        float sepj0i2 = PosSep2(jstp0.Pos, istp2.Pos);
        if(sepj0i2 > PosSep2(jstp0.Pos, istp0.Pos)) continue;
        sepj0i2 = sqrt(sepj0i2);
        float trad = sepj0i2 / radLen;
        // impact parameter between the projection of istj and the jstj charge center
        float delta = PointTrajDOCA(tjs, jstp0.Pos[0], jstp0.Pos[1], istp2);
        // See if delta is consistent with the cone angle of the i shower
        float dang = delta / sepj0i2;
        if(trad > 3) continue;
        if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Candidate "<<iss.ID<<" "<<jss.ID<<" separation "<<sepj0i2<<" radiation lengths "<<trad<<" delta "<<delta<<" dang "<<dang;
        // TODO This needs more work
        // Require that the j energy be much lower
        if(jss.Energy > 0.3 * iss.Energy) continue;
        // There must be a correlation between dang and the energy of the j shower...
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
    
    iss.TjIDs.insert(iss.TjIDs.end(), jss.TjIDs.begin(), jss.TjIDs.end());
    // make a new trajectory using itj as a template
    Trajectory ktj = itj;
    ktj.ID = tjs.allTraj.size() + 1;
    // re-assign the hits from itj to ktj
    for(auto& kpt : ktj.Pts) std::replace(kpt.Hits.begin(), kpt.Hits.end(), itj.ID, ktj.ID);
    // transfer the jtj hits to ktj
    ktj.Pts[0].Hits.insert(ktj.Pts[0].Hits.end(), jtj.Pts[0].Hits.begin(), jtj.Pts[0].Hits.end());
    ktj.Pts[1].Hits.insert(ktj.Pts[1].Hits.end(), jtj.Pts[1].Hits.begin(), jtj.Pts[1].Hits.end());
    ktj.Pts[2].Hits.insert(ktj.Pts[2].Hits.end(), jtj.Pts[2].Hits.begin(), jtj.Pts[2].Hits.end());
    // re-assign the hits from jtj to ktj
    for(auto& kpt : ktj.Pts) {
      std::replace(kpt.Hits.begin(), kpt.Hits.end(), jtj.ID, ktj.ID);
      // Fix InTraj
      for(auto& kht : kpt.Hits) {
        if(tjs.fHits[kht].InTraj == itj.ID || tjs.fHits[kht].InTraj == jtj.ID) tjs.fHits[kht].InTraj = ktj.ID;
      }
    } //  kpt
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
    iss.MatchedTjIDs.insert(iss.MatchedTjIDs.end(), jss.MatchedTjIDs.begin(), jss.MatchedTjIDs.begin());
    std::sort(iss.MatchedTjIDs.begin(), iss.MatchedTjIDs.end());
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
    
    if(!MergeShowersAndStore(fcnLabel, tjs, icotIndex, jcotIndex, prt)) return false;
    // merge was successful. re-find the external parent and start charge
    FindExternalParent(fcnLabel, tjs, icotIndex, prt);
    return true;
    
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
        exit(1);
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
    // points closest to each TrajPoint
    
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
    
    bool success = true;
    for(auto& tp : stj.Pts) {
      // require that every point have charge
      if(tp.Chg == 0) success = false;
      tp.DeltaRMS /= tp.Chg;
      tp.HitPos[0] /= tp.Chg;
      tp.HitPos[1] /= tp.Chg;
    }
    if(success) return true;
    if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" failed. Found no charge in a point on stj "<<stj.ID;
    for(auto& tp : stj.Pts) {
      tp.Chg = 0;
      tp.DeltaRMS = 0;
      tp.NTPsFit = 0;
      tp.HitPos = {0, 0};
    }
    return false;

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
    
    // Expect the shower to be narrow at the start (point 0) and wider at the end (point 2)
    ss.DirectionFOM = stj.Pts[0].DeltaRMS / stj.Pts[2].DeltaRMS;
    // startsNeg is true if this assumption is correct
    bool startsNeg = (stj.Pts[0].DeltaRMS < stj.Pts[2].DeltaRMS);
    
    // reverse the points vector so that the narrow end of the shower is near Pts.begin()
    if(!startsNeg) ReverseShower(fcnLabel, tjs, cotIndex, prt);

    // define the shower start and end positions
    stj.Pts[0].Pos = ss.ShPts[0].Pos;
    unsigned short endPt = ss.ShPts.size()-1;
    stj.Pts[2].Pos = ss.ShPts[endPt].Pos;
    
    // put the charge center where RotPos[0] changes sign
    for(unsigned short ipt = 0; ipt < endPt; ++ipt) {
      if(ss.ShPts[ipt].RotPos[0] * ss.ShPts[ipt + 1].RotPos[0] < 0) {
        stj.Pts[1].Pos = ss.ShPts[ipt].Pos;
        break;
      }
    } // spt
    
    // define the angle of all the shower Tps
    for(auto& stp : stj.Pts) {
      stp.Ang = ss.Angle;
      stp.AngErr = ss.AngleErr;
      stp.Dir[0] = cos(stp.Ang);
      stp.Dir[1] = sin(stp.Ang);
    } // stp
    
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
    
    std::string fcnLabel = inFcnLabel + ".RS";
    
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
  void TagShowerTjs(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList)
  {
    // Tag Tjs with PDGCode = 11 if they have MCSMom < ShowerTag[0] and there are more than
    // ShowerTag[6] other Tjs with a separation < ShowerTag[1]. Returns a list of Tjs that meet this criteria
    
    tjList.clear();
    
    short maxMCSMom = tjs.ShowerTag[1];
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size() - 1; ++it1) {
      Trajectory& tj1 = tjs.allTraj[it1];
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) continue;
      tj1.NNeighbors = 0;
      // identified as a parent
      // ignore shower Tjs
      if(tj1.AlgMod[kShowerTj]) continue;
      // ignore stubby Tjs
      if(tj1.Pts.size() < 3) continue;
      // Cut on length and MCSMom
      if(tj1.Pts.size() > 6 && tj1.MCSMom > maxMCSMom) continue;
      if(tj1.AlgMod[kTjHiVx3Score]) continue;
      if(TjHasNiceVtx(tjs, tj1, tjs.ShowerTag[11])) continue;
      for(unsigned short it2 = it1 + 1; it2 < tjs.allTraj.size(); ++it2) {
        Trajectory& tj2 = tjs.allTraj[it2];
        if(tj2.CTP != inCTP) continue;
        if(tj2.AlgMod[kKilled]) continue;
        // identified as a parent
        // ignore shower Tjs
        if(tj2.AlgMod[kShowerTj]) continue;
        // ignore stubby Tjs
        if(tj2.Pts.size() < 3) continue;
        if(tj2.AlgMod[kTjHiVx3Score]) continue;
        if(TjHasNiceVtx(tjs, tj2, tjs.ShowerTag[11])) continue;
        // Cut on length and MCSMom
        if(tj2.Pts.size() > 6 && tj2.MCSMom > maxMCSMom) continue;
        unsigned short ipt1, ipt2;
        float doca = tjs.ShowerTag[2];
        // Find the separation between Tjs without considering dead wires
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, doca, false);
        if(doca == tjs.ShowerTag[2]) continue;
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
    
    // Add Tjs that are attached to vertices that are attached to Tjs in tjList
    for(unsigned short it = 0; it < tjList.size(); ++it) {
      if(tjList[it].empty()) continue;
      std::vector<int> list;
      for(auto tjID : tjList[it]) {
        if(tjID < 1 || tjID > (int)tjs.allTraj.size()) continue;
        Trajectory& tj = tjs.allTraj[tjID - 1];
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == 0) continue;
          VtxStore& vx2 = tjs.vtx[tj.VtxID[end] - 1];
          // this shouldn't happen but check anyway
          if(vx2.ID == 0) continue;
          // ignore high score vertices
          if(vx2.Score > tjs.ShowerTag[11]) continue;
          // get a list of Tjs attached to this vertex
          auto vxTjs = GetVtxTjIDs(tjs, vx2);
          if(vxTjs.empty()) continue;
          for(auto vtjID : vxTjs) {
            if(std::find(tjList[it].begin(), tjList[it].end(), vtjID) != tjList[it].end()) continue;
            if(std::find(list.begin(), list.end(), vtjID) != list.end()) continue;
            list.push_back(vtjID);
          } // vtj
        } // end
      } // tjID
      if(!list.empty()) tjList[it].insert(tjList[it].end(), list.begin(), list.end());
    } // it
    
  } // TagShowerTjs
  
  ////////////////////////////////////////////////
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
  {
    // Find Tjs that are near the shower but are not included in it
    if(cotIndex > tjs.cots.size() - 1) return;
    ShowerStruct& ss = tjs.cots[cotIndex];
    ss.NearTjIDs.clear();
    
    std::string fcnLabel = inFcnLabel + ".FNTj";
    
    std::vector<int> ntj;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      // not a showerTj
      if(tj.AlgMod[kShowerTj]) continue;
      // ignore stubby Tjs
      if(tj.Pts.size() < 10) continue;
      // make sure it's not in the shower
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
      // methodical but probably  slow
      for(auto& inTjID : ss.TjIDs) {
        unsigned short ipt1, ipt2;
        float doca = tjs.ShowerTag[2];
        Trajectory& inTj = tjs.allTraj[inTjID - 1];
        TrajTrajDOCA(tjs, tj, inTj, ipt1, ipt2, doca, false);
        if(doca < tjs.ShowerTag[2]) {
          ntj.push_back(tj.ID);
          break;
        }
      } // inTjID
    } // tj
    if(!ntj.empty()) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" cotIndex "<<cotIndex<<" stj.ID "<<ss.ShowerTjID<<" found "<<ntj.size()<<" nearby Tjs";
      ss.NearTjIDs = ntj;
    }
    
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
        if(tj.AlgMod[kTjHiVx3Score]) continue;
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
  void AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt)
   {
    // This function adds Tjs to the shower. It updates the shower parameters.
    
     if(cotIndex > tjs.cots.size() - 1) return;
    
     ShowerStruct& ss = tjs.cots[cotIndex];
     if(ss.Envelope.empty()) return;
     if(ss.ID == 0) return;
     if(ss.TjIDs.empty()) return;
     
     std::string fcnLabel = inFcnLabel + ".ATIE";
    
    unsigned short nadd = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kShowerTj]) continue;
     
      if(tj.AlgMod[kTjHiVx3Score]) continue;
      if(TjHasNiceVtx(tjs, tj, tjs.ShowerTag[11])) continue;
      // This shouldn't be necessary but do it for now
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
      // See if both ends are outside the envelope
      bool end0Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[0]].Pos, ss.Envelope);
      bool end1Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[1]].Pos, ss.Envelope);
      if(!end0Inside && !end1Inside) continue;
      if(end0Inside && end1Inside) {
        // Fully contained
        // TODO: See if the Tj direction is compatible with the shower?
        if(AddTj(fcnLabel, tjs, tj.ID, cotIndex, false, prt)) ++nadd;
        ++nadd;
        continue;
      } // both ends inside
      // Require high momentum Tjs be aligned with the shower axis
      // TODO also require high momentum Tjs close to the shower axis?

      if(tj.MCSMom > 500) {
        float tjAngle = tj.Pts[tj.EndPt[0]].Ang;
        float dangPull = std::abs(tjAngle -ss.AngleErr) / ss.AngleErr;
        if(dangPull > 2) continue;
      } // high momentum
      if(AddTj(fcnLabel, tjs, tj.ID, cotIndex, false, prt)) ++nadd;
    } // tj
    
    if(nadd > 0) {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Added "<<nadd<<" trajectories. Calling DefineShower... ";
      DefineShower(fcnLabel, tjs, cotIndex, prt);
      return;
    } else {
      if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" No new trajectories added to envelope ";
      return;
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
        if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
        // inside the tick range?
        if(tjs.fHits[iht].PeakTime < loTick) continue;
        // Note that hits are sorted by increasing time so we can break here
        if(tjs.fHits[iht].PeakTime > hiTick) break;
        // see if this hit is inside the envelope
        point[0] = tjs.fHits[iht].WireID.Wire;
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
      chgPos[0] = tjs.fHits[iht].WireID.Wire - stp1.Pos[0];
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
    
  } // DumpShower
  
  ////////////////////////////////////////////////
  void CheckQuality(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // drop those that don't meet the requirements
    
    std::string fcnLabel = inFcnLabel + ".CQ";
    
    for(unsigned short cotIndex = 0; cotIndex < tjs.cots.size(); ++cotIndex) {
      ShowerStruct& ss = tjs.cots[cotIndex];
      if(ss.CTP != inCTP) continue;
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
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
      if(killit) {
        MakeShowerObsolete(fcnLabel, tjs, cotIndex, prt);
      } else {
        // A good shower. Set the pdgcode of InShower Tjs to 11
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          tj.PDGCode = 11;
          // Clobber 2D vertices that are inside the shower
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] > 0) {
              VtxStore& vx2 = tjs.vtx[tj.VtxID[end]-1];
              bool killMe = (vx2.Score < tjs.ShowerTag[11]);
              // don't kill the 2D vertex if it is attached to the far end of the parent Tj
              if(killMe && ss.ParentID == tjID) {
                unsigned short farEnd = FarEnd(tjs, tj, ss);
                if(farEnd == end) killMe = false;
              }
              if(killMe) {
                if(prt) mf::LogVerbatim("TC")<<fcnLabel<<" Clobber vtx "<<vx2.ID<<" Score "<<vx2.Score<<" Vtx3ID "<<tjs.vtx[tj.VtxID[end]-1].Vtx3ID;
                // force killing the vertex, possibly overriding the settings of Vertex2DCuts
                MakeVertexObsolete(tjs, tj.VtxID[end], true);
              }
            }
          } // end
        }
      } // don't killit
      if (tjs.SaveShowerTree) SaveTjInfo(tjs, inCTP, cotIndex, "CQ");
    } // ic
    
  } // CheckQuality

  ////////////////////////////////////////////////
  bool TransferTjHits(TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Transfer InShower hits to the shower Tj 
    
    bool newShowers = false;
    for(unsigned short ish = 0; ish < tjs.cots.size(); ++ish) {
      ShowerStruct& ss = tjs.cots[ish];
      // Ensure that this is the correct CTP
      if(ss.CTP != inCTP) continue;
      // Ensure that it is valid
      if(ss.ID == 0) continue;
      if(ss.TjIDs.empty()) continue;
      if(ss.ShowerTjID == 0) continue;
      // Tp 1 of stj will get all of the shower hits
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      if(!stj.Pts[1].Hits.empty()) {
        std::cout<<"TTjH: ShowerTj "<<stj.ID<<" already has hits. This can't be right\n";
        continue;
      }
      stj.PDGCode = 11;
      // Note that UseHit is not used since the size is limited to 16
      for(auto& tjID : ss.TjIDs) {
        unsigned short itj = tjID - 1;
        if(tjs.allTraj[itj].VtxID[0] != 0 && tjs.allTraj[itj].VtxID[0] != stj.VtxID[0] && tjs.allTraj[itj].ID != ss.ParentID) {
          std::cout<<"TTjH: Trying to transfer hits on Tj "<<tjID<<" attached to a vertex "<<tjs.allTraj[itj].VtxID[0]<<" at end0 that is not the showerTj vertex or the parent Tj "<<stj.VtxID[0]<<"\n";
          continue;
        }
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
    } // ish
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
    
    MCParticleListUtils tm{tjs};
    
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
    tjs.allTraj.push_back(stj);
    // Create the shower struct
    ShowerStruct ss;
    ss.ID = tjs.cots.size() + 1;
    ss.CTP = stj.CTP;
    // assign all TJ IDs to this ShowerStruct
    ss.TjIDs = tjl;
    ss.ShowerTjID = stj.ID;
    // try to define the true shower parent Tj
    unsigned short nTruHits;
    // Find the MC particle that matches with these InShower Tjs
    unsigned short mcpIndex = tm.GetMCPartListIndex(ss, nTruHits);
    // Find the Tj that is closest to the start of this MC Particle
    if(mcpIndex != USHRT_MAX) ss.TruParentID = tm.MCParticleStartTjID(mcpIndex, ss.CTP);
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
    }
    
    // print a header
    myprt<<someText<<" ict  ID  ParID Energy nTjs  dFOM AspRat   stj __Pos0___ nPts dRMS __Pos1___ nPts dRMS __Pos2___ nPts dRMS Angle SS3ID  > inCTP "<<inCTP<<"\n";

    for(unsigned short ict = 0; ict < tjs.cots.size(); ++ict) {
      const auto& ss = tjs.cots[ict];
      if(!printAllCTP && ss.CTP != inCTP) continue;
      if(!printKilledShowers && ss.ID == 0) continue;
      myprt<<someText<<std::fixed;
      myprt<<std::setw(4)<<ict;
      myprt<<std::setw(4)<<ss.ID;
      myprt<<std::setw(7)<<ss.ParentID;
      myprt<<std::setw(7)<<(int)ss.Energy;
      myprt<<std::setw(5)<<ss.TjIDs.size();
      const auto& stj = tjs.allTraj[ss.ShowerTjID - 1];
      myprt<<std::setw(6)<<std::setprecision(2)<<ss.DirectionFOM;
      myprt<<std::setw(7)<<std::setprecision(2)<<ss.AspectRatio;
      myprt<<std::setw(6)<<stj.ID;
      for(auto& spt : stj.Pts) {
        myprt<<std::setw(10)<<PrintPos(tjs, spt.Pos);
        myprt<<std::setw(5)<<spt.NTPsFit;
        myprt<<std::setw(5)<<std::setprecision(1)<<spt.DeltaRMS;
      } // spt
      myprt<<std::setw(6)<<std::setprecision(2)<<stj.Pts[1].Ang;
      myprt<<std::setw(6)<<ss.SS3ID;
      myprt<<"\n";
    } // ss
  } // Print2DShowers

} // namespace tca
