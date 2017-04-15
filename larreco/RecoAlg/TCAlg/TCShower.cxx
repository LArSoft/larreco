#include "larreco/RecoAlg/TCAlg/TCShower.h"


struct SortEntry{
  unsigned int index;
  float length;
};

bool greaterThan (SortEntry c1, SortEntry c2) { return (c1.length > c2.length);}
bool lessThan (SortEntry c1, SortEntry c2) { return (c1.length < c2.length);}


namespace tca {

  ////////////////////////////////////////////////
  void Find3DShowerEndPoints(TjStuff& tjs, const geo::TPCID& tpcid)
  {
    
    if(tjs.ShowerTag[0] < 0) return;
    
    bool prt = (tjs.ShowerTag[9] == 3);
    
    if(prt) mf::LogVerbatim("TC")<<"Inside FindShowerEndPoints";
    
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;

    for(auto& im : tjs.matchVecPFPList) {
      // a reference to a set of 3D matched trajectories
      auto& ms = tjs.matchVec[im];
      if(ms.TjIDs.empty()) continue;
      // Only consider shower Tjs
      if(ms.PDGCode != 1111) continue;
      // ensure we are in the correct tpcid using the first Tj CTP
      unsigned short it1 = ms.TjIDs[0] - 1;
      geo::PlaneID plane1ID = DecodeCTP(tjs.allTraj[it1].CTP);
      if(plane1ID.Cryostat != cstat) continue;
      if(plane1ID.TPC != tpc) continue;
      // sum of the matching weights
      float wsum = 0;
      ms.sXYZ = {0, 0, 0};
      ms.sDir = {0, 0, 0};
      // the number of end position entries
      unsigned short ne = 0;
      ms.eXYZ = {0, 0, 0};
      TVector3 prevDir = {0, 0, 0};
      for(unsigned short ii = 0; ii < ms.TjIDs.size() - 1; ++ii) {
        // find this ID in the shower struct vector
        unsigned short issIndex = ShowerTjCotsIndex(tjs, ms.TjIDs[ii]);
        if(issIndex == USHRT_MAX) continue;
        ShowerStruct& iss = tjs.cots[issIndex];
        Trajectory& itj = tjs.allTraj[iss.ShowerTjID - 1];
        TrajPoint& itp = itj.Pts[0];
        unsigned short iPln = DecodeCTP(itp.CTP).Plane;
        for(unsigned short jj = ii + 1; jj < ms.TjIDs.size(); ++jj) {
          unsigned short jssIndex = ShowerTjCotsIndex(tjs, ms.TjIDs[jj]);
          if(jssIndex == USHRT_MAX) continue;
          ShowerStruct& jss = tjs.cots[jssIndex];
          Trajectory& jtj = tjs.allTraj[jss.ShowerTjID - 1];
          TrajPoint& jtp = jtj.Pts[0];
          unsigned short jPln = DecodeCTP(jtp.CTP).Plane;
          TVector3 pos, dir;
          TrajPoint3D(tjs, itp, jtp, pos, dir);
          if(dir.X() > 1) {
            if(prt) mf::LogVerbatim("TC")<<" FSEP: TrajPoint3D failed using points "<<iPln<<":"<<PrintPos(tjs, itp.Pos)<<" "<<jPln<<":"<<PrintPos(tjs, jtp.Pos);
            continue;
          }
          // de-weight matches in which the start charge isn't known
          float wght = 1;
          if(iss.StartChg == 0 || jss.StartChg == 0) wght = 0.1;
          wsum += wght;
          ms.sXYZ[0] += wght * pos[0];
          ms.sXYZ[1] += wght * pos[1];
          ms.sXYZ[2] += wght * pos[2];
          // don't allow sign changes to the major direction cosine.
          if(prevDir[0] == 0) {
            prevDir = dir;
          } else {
            // a valid direction cosine exists
            for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
              if(std::abs(prevDir[ixyz]) > 0.5) {
                if(prevDir[ixyz] * dir[ixyz] < 0) dir *= -1;
                break;
              }
            }
          } //  check for direction sign change
          if(prt) {
            mf::LogVerbatim myprt("TC");
            myprt<<" itp "<<iPln<<":"<<PrintPos(tjs, itp.Pos)<<std::fixed<<std::setprecision(2)<<" Angle "<<itp.Ang;
            myprt<<" jtp "<<jPln<<":"<<PrintPos(tjs, jtp.Pos)<<" Angle "<<jtp.Ang;
            myprt<<" dir "<<std::fixed<<std::setprecision(2)<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<" wght "<<wght;
            myprt<<" prevDir "<<std::fixed<<std::setprecision(2)<<prevDir[0]<<" "<<prevDir[1]<<" "<<prevDir[2];
          }
          ms.sDir += wght * dir;
          // now find eXYZ
          TrajPoint& ietp = itj.Pts[2];
          TrajPoint& jetp = jtj.Pts[2];
          float ieX = tjs.detprop->ConvertTicksToX(ietp.Pos[1]/tjs.UnitsPerTick, iPln, tpc, cstat);
          float jeX = tjs.detprop->ConvertTicksToX(jetp.Pos[1]/tjs.UnitsPerTick, iPln, tpc, cstat);
          double yp, zp;
          tjs.geom->IntersectionPoint(std::nearbyint(ietp.Pos[0]), std::nearbyint(jetp.Pos[0]), iPln, jPln, cstat, tpc, yp, zp);
          if(yp < tjs.YLo || yp > tjs.YHi || zp < tjs.ZLo || zp > tjs.ZHi) continue;
          ++ne;
          ms.eXYZ[0] += 0.5 * (ieX + jeX);
          ms.eXYZ[1] += yp;
          ms.eXYZ[2] += zp;
        } // jj
      } // ii
      if(wsum == 0) {
        if(prt) mf::LogVerbatim("TC")<<" FSEP: Failed to find shower start. Skip this match";
        continue;
      }
      for(auto& xyz : ms.sXYZ) xyz /= wsum;
      for(auto& xyz : ms.eXYZ) xyz /= (float)ne;
      ms.sDir.SetMag(1);
      // correct the direction using a large direction cosine
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
        if(std::abs(ms.sDir[ixyz]) > 0.5) {
          if(ms.eXYZ[ixyz] > ms.sXYZ[ixyz] && ms.sDir[ixyz] < 0) ms.sDir *= -1;
          break;
        }
      }
//      if(ms.eXYZ[2] > ms.sXYZ[2] && ms.sDir.Z() < 0) ms.sDir *= -1;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" FSEP: TjIDs";
        for(auto& tjid : ms.TjIDs) myprt<<" "<<tjid;
        myprt<<" Start "<<std::fixed<<std::setprecision(1)<<ms.sXYZ[0]<<" "<<ms.sXYZ[1]<<" "<<ms.sXYZ[2];
        myprt<<" End "<<std::fixed<<std::setprecision(1)<<ms.eXYZ[0]<<" "<<ms.eXYZ[1]<<" "<<ms.eXYZ[2];
        myprt<<" 3D Dir "<<std::setprecision(2)<<ms.sDir.X()<<" "<<ms.sDir.Y()<<" "<<ms.sDir.Z();
      }
      std::cout<<"Shower Dir "<<std::setprecision(2)<<ms.sDir.X()<<" "<<ms.sDir.Y()<<" "<<ms.sDir.Z()<<"\n";
    } // im

  } // Find3DShowerEndPoints

  ////////////////////////////////////////////////
  void MakeShowers(TjStuff& tjs, const calo::CalorimetryAlg& fCaloAlg)
  {
    // Fill 3D shower variables. First look for matching shower Tjs, then use this
    // information to find matching parent Tjs
    
    if(tjs.ShowerTag[0] < 0) return;
    
    // Get the calibration constants
    
    bool prt = (tjs.ShowerTag[9] >= 0);
    
    int shID = 0;
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned short imv = tjs.matchVecPFPList[ipfp];
      auto& ms = tjs.matchVec[imv];
      if(ms.TjIDs.empty()) continue;
      if(ms.PDGCode == 13) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"MS: "<<ipfp<<" ParentMSIndex "<<ms.ParentMSIndex<<" PDGCode "<<ms.PDGCode<<" Dtr size "<<ms.DtrIndices.size();
        myprt<<" TjIDs:";
        for(auto& tjID : ms.TjIDs) myprt<<" "<<tjID;
      } // prt
      // look for matched shower Tjs
      if(ms.PDGCode != 1111) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" Shower Tj";
        for(auto& tjID : ms.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID-1];
          unsigned short endPt = tj.EndPt[0];
          myprt<<"\n "<<tj.ID<<" start "<<PrintPos(tjs, tj.Pts[endPt]);
          endPt = tj.EndPt[1];
          myprt<<" "<<tj.ID<<" end "<<PrintPos(tjs, tj.Pts[endPt]);
        } //  tjID
      } // prt
      ++shID;
      ShowerStruct3D ss3;
      ss3.Energy.resize(tjs.NumPlanes);
      ss3.EnergyErr.resize(tjs.NumPlanes);
      ss3.MIPEnergy.resize(tjs.NumPlanes);
      ss3.MIPEnergyErr.resize(tjs.NumPlanes);
      ss3.dEdx.resize(tjs.NumPlanes);
      ss3.dEdxErr.resize(tjs.NumPlanes);
      ss3.ID = shID;
      // fill the start position
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) ss3.Pos[ixyz] = ms.sXYZ[ixyz];
      // and direction
      ss3.Dir = ms.sDir;
      ss3.DirErr = ms.sDirErr;
      // Find the shower length.
      ss3.Len = 0;
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
        double dpos = ms.eXYZ[ixyz] - ms.sXYZ[ixyz];
        ss3.Len += dpos * dpos;
      }
      ss3.Len = sqrt(ss3.Len);
      // We need the shower structs to fill the variables in each plane
      for(unsigned short ii = 0; ii < ms.TjIDs.size(); ++ii) {
        unsigned short istjID = ms.TjIDs[ii];
        // find this ID in the shower struct vector
        unsigned short iss = ShowerTjCotsIndex(tjs, istjID);
        if(iss == USHRT_MAX) continue;
        geo::PlaneID planeID = DecodeCTP(tjs.cots[iss].CTP);
        unsigned short iPln = planeID.Plane;
        // TODO Calculate energy using fCaloAlg
        ss3.Energy[iPln] = tjs.cots[iss].Energy;
        // This is just a guess for now
        ss3.EnergyErr[iPln] = 0.3 * tjs.cots[iss].Energy;
        // This is probably wrong also...
        ss3.MIPEnergy[iPln] = ss3.Energy[iPln] / 2.3;
        ss3.MIPEnergyErr[iPln] = ss3.EnergyErr[iPln] / 2.3;
        ss3.dEdx[iPln] = 0;
        // Calculate dE/dx
        if(tjs.cots[iss].StartChg > 0) {
          double angleToVert = tjs.geom->WireAngleToVertical(tjs.geom->View(planeID), planeID.TPC, planeID.Cryostat) - 0.5 * ::util::pi<>();
          double cosgamma = std::abs(std::sin(angleToVert) * ss3.Dir.Y() + std::cos(angleToVert) * ss3.Dir.Z());
          if(cosgamma == 0) continue;
          double dx = tjs.geom->WirePitch(planeID) / cosgamma;
          double dQ = tjs.cots[iss].StartChg;
          // Get the time using the shower charge center position
          Trajectory& stj = tjs.allTraj[tjs.cots[iss].ShowerTjID - 1];
          double time = stj.Pts[1].Pos[1] / tjs.UnitsPerTick;
          ss3.dEdx[iPln] = fCaloAlg.dEdx_AREA(dQ, time, dx, iPln);
          if(prt) mf::LogVerbatim("TC")<<"Shower index "<<iss<<" plane "<<iPln<<" plane "<<iPln<<" dQ "<<(int)dQ<<" dx "<<dx<<" dE/dx "<<ss3.dEdx[iPln];
        }
      } // ii
      // Calculate the opening angle here - somehow
      ss3.OpenAngle = 0.1;
      // We shouldn't define the ss3 Hit vector until hit merging is done
      tjs.showers.push_back(ss3);
      
    } // ipfp
  } // MakeShowers

  ////////////////////////////////////////////////
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP)
  {
    // Construct clusters of trajectories (cots) which will become shower PFParticles
    
    // ShowerTag[] parameters
    // 0 Mode (<= 0 OFF, 1 = tag only, 2 = find showers)
    // 1 Max Tj MCSMom for a shower tag (< 0 = no shower-like Tj tagging or shower finding)
    // 2 Max separation
    // 3 Max delta angle
    // 4 rms width factor
    // 5 Min shower 1/2 width (WSE units)
    // 6 Min total Tj Pts
    // 7 Min Tjs
    // 8 Debug in CTP
    
    if(tjs.ShowerTag[0] <= 0) return;
    
    bool prt = false;
    // print only one shower?
    unsigned short prtShower = USHRT_MAX;
    if(tjs.ShowerTag[9] >= 0) {
      geo::PlaneID planeID = DecodeCTP(inCTP);
      CTP_t printCTP = EncodeCTP(planeID.Cryostat, planeID.TPC, std::nearbyint(tjs.ShowerTag[9]));
      prt = (printCTP == inCTP);
      if(printCTP > 2) prt = true;
      if(printCTP > 9) prtShower = printCTP - 10;
    }
    // save the requested print state in case it gets changed
    bool saveprt = prt;
    
    std::vector<std::vector<unsigned short>> tjList;
    TagShowerTjs(tjs, inCTP, tjList);
    if(prt) std::cout<<"Inside FindShowers inCTP "<<inCTP<<" tjList size "<<tjList.size()<<"\n";
    if(tjs.ShowerTag[0] == 1) return;
    if(tjList.empty()) return;

    // Merge the lists of Tjs in showers
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

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"tjlist after merging and removing duplicates\n";
      for(auto& tjl : tjList) {
        if(tjl.empty()) continue;
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
        myprt<<"\n";
      } // tjl
    } // prt
    
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
      Trajectory stj;
      stj.CTP = inCTP;
      // with three points
      stj.Pts.resize(3);
      for(auto& stp : stj.Pts) {
        stp.CTP = stj.CTP;
        // set all UseHit bits true so we don't get a confusing with few hits
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
      ss.CTP = stj.CTP;
      // assign all TJ IDs to this ShowerStruct
      ss.TjIDs = tjl;
      ss.ShowerTjID = stj.ID;
      // put it in TJ stuff. The rest of the info will be added later
      tjs.cots.push_back(ss);
      unsigned short cotIndex = tjs.cots.size() - 1;
      if(prt && prtShower != USHRT_MAX) prt = (cotIndex == prtShower);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Make cots "<<cotIndex<<" in CTP "<<ss.CTP<<" TjID_NN";
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
      }
      DefineShower(tjs, cotIndex, prt);
      FindExternalParent(tjs, cotIndex, prt);
      FindStartChg(tjs, cotIndex, prt);
      if(prt) PrintTrajectory("FS", tjs, tjs.allTraj[stj.ID-1], USHRT_MAX);
    } // tjl
    
    if(tjs.cots.empty()) return;

    // merge showers?
    prt = saveprt;
    if(tjs.ShowerTag[0] > 2) MergeShowers(tjs, inCTP, prt);
    
    // drop those that don't meet the requirements
    for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
      ShowerStruct& ss = tjs.cots[ic];
      if(ss.CTP != inCTP) continue;
      if(ss.TjIDs.empty()) continue;
      // enough Tjs?
      unsigned short ntjs = ss.TjIDs.size();
      bool killit = (ntjs < tjs.ShowerTag[7]);
      // Kill runt showers
      if(prt && prtShower != USHRT_MAX) prt = (ic == prtShower);
      if(!killit) killit = (ss.Energy < tjs.ShowerTag[3]);
      if(prt) mf::LogVerbatim("TC")<<"ic "<<ic<<" nTjs "<<ss.TjIDs.size()<<" nTjs "<<ss.TjIDs.size()<<" killit? "<<killit;
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
        // Unset shower and killed bits. Trajectories that are in showers haven't had their hits re-assigned to the
        // shower Tj yet so nothing needs to be done to them
        for(auto& tjID : ss.TjIDs) tjs.allTraj[tjID - 1].AlgMod[kKilled] = false;
        // kill the shower Tj
        ss.TjIDs.clear();
        unsigned short itj = ss.ShowerTjID - 1;
        if(prt) mf::LogVerbatim("TC")<<" killing ShowerTj "<<tjs.allTraj[itj].ID<<". Restored InShower Tjs.";
        MakeTrajectoryObsolete(tjs, itj);
      } else {
        if(tjs.allTraj[ss.ShowerTjID - 1].AlgMod[kKilled]) {
          std::cout<<"FS logic error: ShowerTj "<<tjs.allTraj[ss.ShowerTjID - 1].ID<<" is killed\n";
          tjs.allTraj[ss.ShowerTjID - 1].AlgMod[kKilled] = false;
        }
        // A good shower. Set the pdgcode of InShower Tjs to 11
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          tj.PDGCode = 11;
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] > 0) MakeVertexObsolete(tjs, tj.VtxID[end]);
          } // end
        }
      } // don't killit
    } // ic
    
    // Finish up in this CTP. 
    // Re-assign hits from the InShower Tjs to the ShowerTj.
    TransferTjHits(tjs, inCTP, prt);
    std::cout<<"Final calculation shower energy...\n";

    // check for consistency
    for(auto& ss : tjs.cots) {
      if(ss.TjIDs.empty()) continue;
      if(ss.CTP != inCTP) continue;
      if(ss.ShowerTjID == 0) {
        std::cout<<"FindShowers: ShowerTjID not defined in CTP "<<ss.CTP<<"\n";
      }
      for(auto& tjID : ss.TjIDs) {
        if(tjID > tjs.allTraj.size()) {
          std::cout<<"FindShowers: Bad tjID "<<tjID<<"\n";
        }
        Trajectory& tj = tjs.allTraj[tjID - 1];
        if(tj.CTP != ss.CTP) {
          std::cout<<"FindShowers: Bad CTP "<<ss.CTP<<" "<<tj.CTP<<" tjID "<<tjID<<"\n";
        }
        if(!tj.AlgMod[kKilled] || !tj.AlgMod[kInShower]) {
          std::cout<<"FindShowers: InShower TjID "<<tjID<<" invalid kKilled "<<tj.AlgMod[kKilled]<<" or kInShower "<<tj.AlgMod[kInShower]<<"\n";
          PrintTrajectory("FS", tjs, tj, USHRT_MAX);
        }
      } // tjID
    } // ss
    
    if(tjs.ShowerTag[9] >= 0) {
      for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
        if(tjs.cots[ic].TjIDs.empty()) continue;
        unsigned short itj = tjs.cots[ic].ShowerTjID - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(prt || (tjs.ShowerTag[9] == 3 && tj.CTP == inCTP)) PrintTrajectory("FSO", tjs, tj, USHRT_MAX);
      } // ic
    } // print trajectories

    
  } // FindShowers

  ////////////////////////////////////////////////
  void DefineShower(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Defines the properties of a shower using the trajectory points within the trajectories listed
    // in TjIDs. This wipes out any other information that may exist
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    unsigned short cnt = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      if(itj > tjs.allTraj.size() - 1) {
        mf::LogWarning("TC")<<"Bad TjID "<<ss.TjIDs[it];
        ss.TjIDs.clear();
        return;
      }
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.CTP != ss.CTP) {
        mf::LogWarning("TC")<<"Tj "<<tj.ID<<" is in the wrong CTP "<<tj.CTP<<" "<<ss.CTP;
        ss.TjIDs.clear();
        return;
      }
      if(tj.AlgMod[kShowerTj]) {
        mf::LogWarning("TC")<<"DSTj: Tj "<<tj.ID<<" is in TjIDs in cotIndex "<<cotIndex<<" but is a ShowerTj! Killing it";
        ss.TjIDs.clear();
        return;
      }
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg > 0) ++cnt;
      } // ipt
    } // it
    
    // Add any loose hits (those not in trajectory points) that are stashed in shower Tj Pt[0]
    TrajPoint& stp0 = tjs.allTraj[ss.ShowerTjID - 1].Pts[0];
    cnt += stp0.Hits.size();
    
    ss.Pts.resize(cnt);
    
    // Now populate the vectors with the information we currently have
    cnt = 0;
    float totChg = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      Trajectory& tj = tjs.allTraj[ss.TjIDs[it] - 1];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        ss.Pts[cnt].Pos = tp.HitPos;
        ss.Pts[cnt].Chg = tp.Chg;
        ss.Pts[cnt].TID = tj.ID;
        totChg += tp.Chg;
        ++cnt;
      } // ipt
    } // it
    
    // include the loose hits
    for(auto& iht : stp0.Hits) {
      ss.Pts[cnt].Pos[0] = tjs.fHits[iht].WireID.Wire;
      ss.Pts[cnt].Pos[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
      ss.Pts[cnt].Chg = tjs.fHits[iht].Integral;
      totChg += ss.Pts[cnt].Chg;
      ss.Pts[cnt].TID = 0;
      ++cnt;
    }
    
    // Put the total charge into the shower Tj
    tjs.allTraj[ss.ShowerTjID - 1].AveChg = totChg;
    
    if(prt) mf::LogVerbatim("TC'")<<"DS: cotIndex "<<cotIndex<<" filled "<<cnt<<" points including "<<stp0.Hits.size()<<" loose hits. Total charge "<<(int)totChg;

    UpdateShower(tjs, cotIndex, prt);

  } // DefineShower
  
  ////////////////////////////////////////////////
  bool AddTj(TjStuff& tjs, unsigned short TjID, const unsigned short& cotIndex, bool doUpdate, bool prt)
  {
    // Adds the Tj to the shower and optionally updates the shower variables
    
    if(TjID > tjs.allTraj.size()) return false;
    if(cotIndex > tjs.cots.size() - 1) return false;

    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[TjID - 1];
    if(tj.AlgMod[kInShower]) {
      mf::LogWarning("TC")<<"AddTj: Tj "<<TjID<<" is already an InShower Tj";
      return false;
    }
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), TjID) != ss.TjIDs.end()) {
      mf::LogWarning("TC")<<"AddTj: Tj "<<TjID<<" is already in this shower "<<cotIndex;
      return false;
    }
    if(tj.CTP != ss.CTP) {
      mf::LogWarning("TC")<<"AddTj: Tj "<<TjID<<" is in the wrong CTP "<<cotIndex;
      return false;
    }
    ss.TjIDs.push_back(TjID);
    // count the TPs
    unsigned short cnt = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      ++cnt;
    } // ipt
    unsigned short newSize = ss.Pts.size() + cnt;
    cnt = ss.Pts.size();
    ss.Pts.resize(newSize);
    // now add them
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      ss.Pts[cnt].Pos = tp.HitPos;
      ss.Pts[cnt].Chg = tp.Chg;
      ss.Pts[cnt].TID = tj.ID;
      ++cnt;
    } // ipt
    tj.AlgMod[kInShower] = true;
    
    if(doUpdate) UpdateShower(tjs, cotIndex, prt);
    
    return true;
  } // AddTj
  
  ////////////////////////////////////////////////
  bool RemoveTj(TjStuff& tjs, unsigned short TjID, const unsigned short& cotIndex, bool prt)
  {
    // Removes the Tj from a shower
    
    if(TjID > tjs.allTraj.size()) return false;
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    // make sure it isn't already in a shower
    Trajectory& tj = tjs.allTraj[TjID - 1];
    if(!tj.AlgMod[kInShower]) {
//      mf::LogWarning("TC")<<"RemoveTj: Tj "<<TjID<<" is not an InShower Tj";
      return false;
    }
    tj.AlgMod[kInShower] = false;
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
    // re-build everything
    DefineShower(tjs, cotIndex, prt);
    
    return true;
  } // RemoveTj
  
  ////////////////////////////////////////////////
  bool UpdateShower(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    if(!FindChargeCenter(tjs, cotIndex, prt)) {
      mf::LogWarning("TC")<<"Failed to find shower charge center";
      return false;
    }
    FindAngle(tjs, cotIndex, prt);
    FillRotPos(tjs, cotIndex, prt);
    if(!DefineShowerTj(tjs, cotIndex, prt)) {
      mf::LogWarning("TC")<<"Failed to define Shower Tj";
      return false;
    }
    // iterate a few times to nudge the angle into position
    unsigned short nit = 0;
    bool needsUpdate = true;
    RefineShowerTj(tjs, cotIndex, prt);
    while(needsUpdate && nit < 3) {
      FillRotPos(tjs, cotIndex, prt);
      DefineShowerTj(tjs, cotIndex, prt);
      needsUpdate = RefineShowerTj(tjs, cotIndex, prt);
      ++nit;
    } // needsUpdate
    DefineEnvelope(tjs, cotIndex, prt);
    
    // Try to add more Tjs to the shower
    if(AddTjsInsideEnvelope(tjs, cotIndex, prt)) {
      FindChargeCenter(tjs, cotIndex, prt);
      FindAngle(tjs, cotIndex, prt);
      FillRotPos(tjs, cotIndex, prt);
      DefineShowerTj(tjs, cotIndex, prt);
      unsigned short nit = 0;
      bool needsUpdate = RefineShowerTj(tjs, cotIndex, prt);
      while(needsUpdate && nit < 3) {
        FillRotPos(tjs, cotIndex, prt);
        DefineShowerTj(tjs, cotIndex, prt);
        needsUpdate = RefineShowerTj(tjs, cotIndex, prt);
        ++nit;
      } // needsUpdate
      DefineEnvelope(tjs, cotIndex, prt);
    } // Tjs added to the shower
    
    // Try to add loose hits to the shower
    if(AddLooseHits(tjs, cotIndex, prt)) {
      FindChargeCenter(tjs, cotIndex, prt);
      FindAngle(tjs, cotIndex, prt);
      FillRotPos(tjs, cotIndex, prt);
      DefineShowerTj(tjs, cotIndex, prt);
    }
   
    return true;
  } // UpdateShower

  
  ////////////////////////////////////////////////
  void FindExternalParent(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Look for a parent trajectory that starts outside the shower and ends inside
    
    // ShowerTag[] parameters
    // 0 Mode (<= 0 OFF, 1 = tag only, 2 = find showers)
    // 1 Max Tj MCSMom for a shower tag (< 0 = no shower-like Tj tagging or shower finding)
    // 2 Max separation
    // 3 Max delta angle
    // 4 rms width factor
    // 5 Min shower 1/2 width (WSE units)
    // 6 Min total Tj Pts
    // 7 Min Tjs
    // 8 Debug in CTP

    if(cotIndex > tjs.cots.size() - 1) return;
    ShowerStruct& ss = tjs.cots[cotIndex];
    // Ensure that it is valid
    if(ss.TjIDs.empty()) return;
    if(ss.Envelope.empty()) return;
    // An external parent was already added so don't try again
    if(ss.ParentID > 0) return;
    // References to shower Tj points
    TrajPoint& stp0 = tjs.allTraj[ss.ShowerTjID - 1].Pts[0];
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
    TrajPoint& stp2 = tjs.allTraj[ss.ShowerTjID - 1].Pts[2];
    
    float bestFOM = 20;
    unsigned short imTheBest = USHRT_MAX;
    unsigned short imTheBestPt = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      // ignore Tjs that are already in showers or are shower Tjs
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
      if(tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kShowerTj]) continue;
      // Ignore short Tjs
      if(tj.Pts.size() < 5) continue;
      // find the point that is farthest from stp1
      float sep0 = PosSep2(tj.Pts[tj.EndPt[0]].Pos, stp1.Pos);
      float sep1 = PosSep2(tj.Pts[tj.EndPt[1]].Pos, stp1.Pos);
      unsigned short useEnd = 0;
      if(sep1 > sep0) useEnd = 1;
      // ensure that this end is closer to point 0 of the shower Tj than point 2
      TrajPoint& ptp = tj.Pts[tj.EndPt[useEnd]];
      if(PosSep2(ptp.Pos, stp0.Pos) > PosSep2(ptp.Pos, stp2.Pos)) continue;
      float fom = ParentFOM(tjs, tj, useEnd, ss, prt);
      if(fom > bestFOM) continue;
      bestFOM = fom;
      imTheBest = tj.ID;
      imTheBestPt = tj.EndPt[useEnd];
    } // tj
    
    if(imTheBest == USHRT_MAX) return;
    if(bestFOM > tjs.ShowerTag[8]) {
      if(prt) mf::LogVerbatim("TC")<<"FEP: No external parent Tj found ";
      return;
    }
    
    if(prt) mf::LogVerbatim("TC")<<"FEP: Adding external parent "<<imTheBest<<" to the shower with FOM = "<<bestFOM;
    ss.ParentID = imTheBest;
    // Move the shower start point to the parent
    TrajPoint& ptp = tjs.allTraj[imTheBest - 1].Pts[imTheBestPt];
    stp0.Pos = ptp.Pos;
    // Set the shower angle to the parent angle
    ss.Angle = ptp.Ang;
    // Put stp1.Pos on the parent trajectory
    double dist = PosSep(ptp.HitPos, stp1.HitPos);
    stp1.Pos[0] = ptp.HitPos[0] + cos(ptp.Ang) * dist;
    stp1.Pos[1] = ptp.HitPos[1] + sin(ptp.Ang) * dist;
    if(prt) mf::LogVerbatim("TC")<<" stp1 HitPos "<<(int)stp1.HitPos[0]<<":"<<(int)(stp1.HitPos[1]/tjs.UnitsPerTick)<<" Pos "<<PrintPos(tjs, stp1);
    // Add it to the shower and update
    AddTj(tjs, imTheBest, cotIndex, true, prt);
    
  } // FindExternalParent

  ////////////////////////////////////////////////
  float ParentFOM(TjStuff& tjs, Trajectory& tj, const unsigned short& tjEnd, ShowerStruct& ss, bool prt)
  {
    // returns a FOM for the trajectory at the end point being the parent of ss
    
    if(tjEnd > 1) return 1000;
    if(ss.Energy == 0) return 1000;
    
    // Radiation length converted to WSE units (for uB)
    constexpr float radLen = 14 / 0.3;
    constexpr float tenRadLen2 = 100 * radLen * radLen;

    if(ss.TjIDs.empty()) return 1000;
    if(ss.ShowerTjID == 0) return 1000;
    
    // prospective parent TP
    unsigned short endPt = tj.EndPt[tjEnd];
    TrajPoint& ptp = tj.Pts[endPt];
    // Shower charge center TP
    unsigned short istj = ss.ShowerTjID - 1;
    TrajPoint& stp1 = tjs.allTraj[istj].Pts[1];
    float tp1Sep = PosSep2(ptp.Pos, stp1.Pos);
    // Make a rough cut on radiation lengths
    if(tp1Sep > tenRadLen2) {
      if(prt) mf::LogVerbatim("TC")<<"PFOM "<<tj.ID<<" failed sep cut "<<(int)sqrt(tp1Sep)<<" 10 radiation lengths "<<(int)sqrt(tenRadLen2);
      return 100;
    }
    tp1Sep = sqrt(tp1Sep);
    
    // impact parameter between the projection of stp1 and the ptp position
    float delta = PointTrajDOCA(tjs, ptp.Pos[0], ptp.Pos[1], stp1);
    // make a rough cut
    if(delta > 100) {
      if(prt) mf::LogVerbatim("TC")<<"PFOM "<<tj.ID<<" failed delta cut "<<delta<<" cut = 100";
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
    float expectedTPSepRMS = 0.5 * expectedTPSep;
    
    float sepPull = (tp1Sep - expectedTPSep) / expectedTPSepRMS;
    float deltaErr = tp1Sep * ss.AngleErr;
    // don't let the delta error exceed the shower width
    float showerWidth = tjs.allTraj[istj].Pts[0].DeltaRMS;
    if(deltaErr > showerWidth) deltaErr = showerWidth;
    float deltaPull = delta / deltaErr;
    float dang = DeltaAngle(ptp.Ang, stp1.Ang);
    float dangErr = sqrt(ss.AngleErr * ss.AngleErr + ptp.AngErr * ptp.AngErr);
    float dangPull = dang / dangErr;
    float mom = tj.MCSMom;
    if(mom > 300) mom = 300;
    float momPull = (mom - 300) / 100;
    float fom = sqrt(sepPull * sepPull + deltaPull * deltaPull + dangPull * dangPull + momPull * momPull);
    fom /= 4;
    
    // Add properties of the shower. We expect the shower to be narrow at the beginning and containing
    // a reasonable fraction of the charge. It should be wide at the end and have less charge.
    // Determine which shower Tj point is closest to the end point
/*
    unsigned short spt = 0;
    if(PosSep2(ptp.Pos, tjs.allTraj[istj].Pts[2].Pos) < PosSep2(ptp.Pos, tjs.allTraj[istj].Pts[0].Pos)) spt = 2;
    // index of the other shower point end
    unsigned short ospt = 2 - spt;
    TrajPoint& nearPt = tjs.allTraj[istj].Pts[spt];
    TrajPoint& farPt = tjs.allTraj[istj].Pts[ospt];
    float widRat = farPt.DeltaRMS / nearPt.DeltaRMS;
*/
    float widRat = tjs.allTraj[istj].Pts[2].DeltaRMS / tjs.allTraj[istj].Pts[0].DeltaRMS;
    // TODO this needs to be done better
    if(widRat > 1.5) widRat = 1.5;
    if(widRat < 0.66) widRat = 0.66;
//    float chgRat = nearPt.Chg / farPt.Chg;
    float chgRat = tjs.allTraj[istj].Pts[0].Chg / tjs.allTraj[istj].Pts[2].Chg;
    if(chgRat > 1.5) chgRat = 1.5;
    if(chgRat < 0.66) chgRat = 0.66;
    fom /= widRat;
    fom /= chgRat;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"PFOM: Tj "<<tj.ID<<" Pos "<<PrintPos(tjs, ptp)<<" Energy "<<(int)ss.Energy;
      myprt<<std::fixed<<std::setprecision(2);
      myprt<<" tp1Sep "<<tp1Sep<<" sepPull "<<sepPull;
      myprt<<" delta "<<delta<<" deltaPull "<<deltaPull;
      myprt<<" dang "<<dang<<" dangPull "<<dangPull;
      myprt<<" mcsmom "<<(int)mom<<" momPull "<<momPull;
      myprt<<" widRat "<<widRat<<" chgRat "<<chgRat;
      myprt<<" FOM "<<fom;
    }
    return fom;
  } // ParentFOM

  ////////////////////////////////////////////////
  void MergeShowers(TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Merge showers that point roughly in the same direction and aren't too far apart
    
    // ShowerTag[] parameters
    // 0 Mode (<= 0 OFF, 1 = tag only, 2 = find showers)
    // 1 Max Tj MCSMom for a shower tag (< 0 = no shower-like Tj tagging or shower finding)
    // 2 Max separation
    // 3 Max delta angle
    // 4 rms width factor
    // 5 Min shower 1/2 width (WSE units)
    // 6 Min total Tj Pts
    // 7 Min Tjs
    // 8 Debug in CTP
    
    // Require that the maximum separation is about two radiation lengths
    if(prt) mf::LogVerbatim("TC")<<"MergeShowers checking using separation cut "<<tjs.ShowerTag[2];
    
    float sepCut2 = tjs.ShowerTag[2] * tjs.ShowerTag[2];
    
    bool didMerge = true;
    while(didMerge) {
      didMerge = false;
      // See if the envelopes overlap
      for(unsigned short ict = 0; ict < tjs.cots.size() - 1; ++ict) {
        ShowerStruct& iss = tjs.cots[ict];
        if(iss.TjIDs.empty()) continue;
        if(iss.CTP != inCTP) continue;
        for(unsigned short jct = ict + 1; jct < tjs.cots.size(); ++jct) {
          ShowerStruct& jss = tjs.cots[jct];
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
          if(prt) mf::LogVerbatim("TC")<<" Envelopes "<<ict<<" "<<jct<<" overlap? "<<doMerge;
          if(!doMerge) {
            // check proximity between the envelopes
            for(auto& ivx : iss.Envelope) {
              for(auto& jvx : jss.Envelope) {
                if(PosSep2(ivx, jvx) < sepCut2) {
                  if(prt) mf::LogVerbatim("TC")<<" Envelopes "<<ict<<" "<<jct<<" are close "<<PosSep(ivx, jvx)<<" cut "<<tjs.ShowerTag[2];
                  doMerge = true;
                  break;
                }
              } // jvx
              if(doMerge) break;
            } // ivx
          } // !domerge
          if(!doMerge) continue;
          if(prt) mf::LogVerbatim("TC")<<" Merge them. Re-find shower center, etc. \n";
          // Move all of the Tjs from jct to ict
          iss.TjIDs.insert(iss.TjIDs.end(), jss.TjIDs.begin(), jss.TjIDs.end());
          DefineShower(tjs, ict, prt);
          // kill the shower Tj
          Trajectory& jtj = tjs.allTraj[jss.ShowerTjID - 1];
          jtj.AlgMod[kKilled] = true;
          // erase jct
          tjs.cots.erase(tjs.cots.begin() + jct);
          if(prt) {
            mf::LogVerbatim("TC")<<" ShowerTj after merge";
            Trajectory& itj = tjs.allTraj[iss.ShowerTjID - 1];
            PrintTrajectory("jToi", tjs, itj, USHRT_MAX);
          }
          didMerge = true;
        } // jct
        if(didMerge) break;
      } // ict
    } // didMerge
    
  } // MergeShowers
  
  ////////////////////////////////////////////////
  bool MergeShowersAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt)
  {
    // This function is called from MergeAndStore whose function is to merge two line-like
    // trajectories and store them. This function was called because at least one of the
    // trajectories is a shower Tj. Assume that the decision to merge them has been made elsewhere.
    
    if(prt) mf::LogVerbatim("TC")<<"MSAS: MergeShowerAndStore istj "<<istj<<" jstj "<<jstj;
    
    if(istj > tjs.allTraj.size() - 1) return false;
    if(jstj > tjs.allTraj.size() - 1) return false;
    
    Trajectory& itj = tjs.allTraj[istj];
    Trajectory& jtj = tjs.allTraj[jstj];
    
    // First we check to make sure that both are shower Tjs.
    if(!itj.AlgMod[kShowerTj] && !jtj.AlgMod[kShowerTj]) return false;
    
    // We need to keep the convention used in MergeAndStore to create a new merged trajectory
    // and kill the two fragments. This doesn't require making a new shower however. We can just
    // re-purpose one of the showers
    
    unsigned short icotIndex = ShowerTjCotsIndex(tjs, itj.ID);
    ShowerStruct& iss = tjs.cots[icotIndex];
    unsigned short jcotIndex = ShowerTjCotsIndex(tjs, jtj.ID);
    ShowerStruct& jss = tjs.cots[jcotIndex];
    
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
    MakeTrajectoryObsolete(tjs, istj);
    MakeTrajectoryObsolete(tjs, jstj);
    if(prt) mf::LogVerbatim("TC")<<" killed "<<istj+1<<" and "<<jstj+1<<" new Tj "<<ktj.ID;
    // revise the shower
    iss.ShowerTjID = ktj.ID;
    iss.TjIDs.insert(iss.TjIDs.end(), jss.TjIDs.begin(), jss.TjIDs.begin());
    iss.ParentID = 0;
    jss.TjIDs.clear();
    DefineShower(tjs, icotIndex, prt);
    FindExternalParent(tjs, icotIndex, prt);
    FindStartChg(tjs, icotIndex, prt);
    return true;
    
  } // MergeShowersAndStore

  ////////////////////////////////////////////////
  bool FindChargeCenter(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Finds the charge center using all sub-structure trajectories in the cot. All of the shower
    // charge is assigned to the second TP and the charge weighted position is put in stp1.HitPos
    // and stp1.Pos
    // The charge will later be distributed between TP0 - TP2.
    // The total charge is stored in  shower Tj AveChg.
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return false;
    
    unsigned short stjIndex = ss.ShowerTjID - 1;
    if(stjIndex > tjs.allTraj.size() - 1) return false;
    if(tjs.allTraj[stjIndex].Pts.size() != 3) return false;
    
    // initialize all of the points, except the first one if there is an external parent
    for(unsigned short ii = 0; ii < 3; ++ii) {
      if(ss.ParentID > 0 && ii == 0) continue;
      TrajPoint& tp = tjs.allTraj[stjIndex].Pts[ii];
      tp.Chg = 0;
      tp.HitPos[0] = 0;
      tp.HitPos[1] = 0;
    }
    
    TrajPoint& stp1 = tjs.allTraj[stjIndex].Pts[1];
    
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      if(ss.Pts[ii].Chg <= 0) {
        std::cout<<"FCC: Found point with no charge. This shouldn't happen\n";
        exit(1);
      }
      stp1.Chg += ss.Pts[ii].Chg;
      stp1.HitPos[0] += ss.Pts[ii].Chg * ss.Pts[ii].Pos[0];
      stp1.HitPos[1] += ss.Pts[ii].Chg * ss.Pts[ii].Pos[1];
    } // ii
    
    stp1.HitPos[0] /= stp1.Chg;
    stp1.HitPos[1] /= stp1.Chg;
    stp1.Pos = stp1.HitPos;
    // Use the trajectory AveChg variable to store the total charge including that of a primary Tj
    // if it isn't identified as an InShower Tj
    tjs.allTraj[stjIndex].AveChg = stp1.Chg;
    ss.Energy = ShowerEnergy(tjs, ss);
    if(prt) mf::LogVerbatim("TC")<<"FCC: "<<cotIndex<<" Pos "<<PrintPos(tjs, stp1.Pos)<<" stp1.Chg "<<(int)stp1.Chg<<" Energy "<<(int)ss.Energy<<" MeV";
    return true;
  } // FindChargeCenter

  ////////////////////////////////////////////////
  void FindAngle(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Find the angle of the shower using the position of all of the TPs
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    unsigned short stjIndex = ss.ShowerTjID - 1;
    if(stjIndex > tjs.allTraj.size() - 1) return;
    if(tjs.allTraj[stjIndex].Pts.size() != 3) return;
    
    if(ss.ParentID > 0) {
      if(prt) mf::LogVerbatim("TC")<<"FA: Using external parent Tj angle ";
      return;
    }
    
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
    
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      // Weight by distance^2 from the charge center
      wt = PosSep2(ss.Pts[ii].Pos, stp1.Pos);
      sum  += wt;
      xx = wt * (ss.Pts[ii].Pos[0] - stp1.Pos[0]);
      yy = wt * (ss.Pts[ii].Pos[1] - stp1.Pos[1]);
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      sumx += wt * xr;
      sumy += wt * yr;
      sumx2 += wt * xr * xr;
      sumy2 += wt * yr * yr;
      sumxy += wt * xr * yr;
    } // ii
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
/*
    // reverse by pi to be consistent with a previously found angle
    if(std::abs(angle - ss.Angle) > 1.57) {
      if(angle < 0) { angle += M_PI; } else { angle -= M_PI; }
    }
*/
    ss.Angle += dang;
    ss.AngleErr = 1.5;
    double ndof = ss.Pts.size() - 1;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    if(varnce > 0) ss.AngleErr = sqrt(varnce * sum / delta);
    // TODO need a scheme for estimating the angle error
    if(ss.AngleErr < 0.1) ss.AngleErr = 0.1;
    if(prt) mf::LogVerbatim("TC")<<"FA: Pos "<<PrintPos(tjs, stp1)<<" dang "<<dang<<" Angle "<<ss.Angle<<" Err "<<ss.AngleErr<<" npts "<<ss.Pts.size();
    
  } // FindAngle
  
  ////////////////////////////////////////////////
  void FillRotPos(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Fills the RotPos vector and sorts the points along the shower axis. Note that the rotation is
    // done around stp1.Pos but the charge center is at stp1.HitPos. Pos and HitPos will be exactly the
    // same if there is no parent. The Pos position may be shifted slightly in FindExternalParent so that
    // the parent trajectory lies on the central axis of the shower. This is done so that the charge at the
    // start of the shower is calculated correctly using the parent trajectory points
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return;

    // Determine the size of the shower along the axis and transverse to it. 
    // Rotate and translate each point into the coordinate system defined by tp[1]
    float cs = cos(-ss.Angle);
    float sn = sin(-ss.Angle);

    TrajPoint& stp1 = stj.Pts[1];
    
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      ss.Pts[ii].RotPos[0] = ss.Pts[ii].Pos[0] - stp1.Pos[0];
      ss.Pts[ii].RotPos[1] = ss.Pts[ii].Pos[1] - stp1.Pos[1];
      // Rotate into the stp1 direction
      float along = cs * ss.Pts[ii].RotPos[0] - sn * ss.Pts[ii].RotPos[1];
      float trans = sn * ss.Pts[ii].RotPos[0] + cs * ss.Pts[ii].RotPos[1];
      ss.Pts[ii].RotPos[0] = along;
      ss.Pts[ii].RotPos[1] = trans;
    } // ii
    
    std::vector<SortEntry> sortVec(ss.Pts.size());
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].length = ss.Pts[ii].RotPos[0];
    }
    std::sort(sortVec.begin(), sortVec.end(), lessThan);
    
    // put the points vector into the sorted order
    auto tPts = ss.Pts;
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      unsigned short indx = sortVec[ii].index;
      ss.Pts[ii] = tPts[indx];
    } // ii

  } // FillRotPos
  
  ////////////////////////////////////////////////
  bool DefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Defines the Shower Tj, calculates the shower aspect ratio, etc. This function
    // doesn't change the state of Parent
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return false;

    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return false;
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return false;
    
    bool hasParent = (ss.ParentID > 0);
    
    // Analyse RotPos to determine the shower direction. Calculate the aspect ratio while we are here
    float minAlong = ss.Pts[0].RotPos[0];
    float maxAlong = ss.Pts[ss.Pts.size()-1].RotPos[0];
    float sectionLength = (maxAlong - minAlong) / 3;
    float sec0 = minAlong + sectionLength;
    float sec2 = maxAlong - sectionLength;
    float chgNeg = 0;
    float transRMSNeg = 0;
    float chgPos = 0;
    float transRMSPos = 0;
    float alongSum = 0;
    float transSum = 0;
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      alongSum += std::abs(ss.Pts[ii].RotPos[0]);
      transSum += std::abs(ss.Pts[ii].RotPos[1]);
      if(ss.Pts[ii].RotPos[0] < sec0) {
        chgNeg += ss.Pts[ii].Chg;
        transRMSNeg += ss.Pts[ii].Chg * std::abs(ss.Pts[ii].RotPos[1]);
      } else if(ss.Pts[ii].RotPos[0] > sec2) {
        chgPos += ss.Pts[ii].Chg;
        transRMSPos += ss.Pts[ii].Chg * std::abs(ss.Pts[ii].RotPos[1]);
      }
    } // ii
    if(chgNeg == 0 || chgPos == 0) return false;
    transRMSNeg /= chgNeg;
    transRMSPos /= chgPos;
    ss.AspectRatio = transSum / alongSum;

    if(prt) mf::LogVerbatim("TC")<<"DSTj: "<<cotIndex<<" transRMSNeg "<<transRMSNeg<<" transRMSPos "<<transRMSPos<<" AspectRatio "<<ss.AspectRatio;
    
    if(hasParent && transRMSPos < transRMSNeg) {
      if(prt) mf::LogVerbatim("TC")<<" The external parent Tj and the shower are inconsistent for cotIndex "<<cotIndex<<". Dropping the parent constraint";
      hasParent = false;
      ss.ParentID = 0;
    }
    
    // reverse the points vector so that the narrow end of the shower is near Pts.begin()
    if(!hasParent && (transRMSPos < transRMSNeg)) {
      std::reverse(ss.Pts.begin(), ss.Pts.end());
      // change the sign of RotPos
      for(auto& sspt : ss.Pts) {
        sspt.RotPos[0] = -sspt.RotPos[0];
        sspt.RotPos[1] = -sspt.RotPos[1];
      }
      std::swap(transRMSNeg, transRMSPos);
      // flip the shower angle
      if(ss.Angle > 0) {
        ss.Angle -= M_PI;
      } else {
        ss.Angle += M_PI;
      }
      if(prt) mf::LogVerbatim("TC")<<" Reversed everything. Shower angle = "<<ss.Angle;
    } // reverse everything
    
    // Find the shower start by looking for the first point that has a small transverse distance from the shower spine. Grab
    // the longitudinal position of that point
    transRMSNeg /= 2;
    minAlong = 0;
    for(unsigned short ii = 0; ii < ss.Pts.size(); ++ii) {
      if(std::abs(ss.Pts[ii].RotPos[1]) < transRMSNeg) {
        minAlong = ss.Pts[ii].RotPos[0];
        break;
      }
    } // sspt
    if(minAlong == 0) return false;
    
    // define the angle of all the shower Tps
    for(auto& stp : stj.Pts) {
      stp.Ang = ss.Angle;
      stp.AngErr = ss.AngleErr;
      stp.Dir[0] = cos(stp.Ang);
      stp.Dir[1] = sin(stp.Ang);
      stp.Chg = 0;
      stp.DeltaRMS = 0;
      stp.NTPsFit = 0;
    } // stp
    
    TrajPoint& stp0 = stj.Pts[0];
    TrajPoint& stp1 = stj.Pts[1];

    // Put first shower Tj point on the shower axis using the minAlong position
    if(!hasParent) {
      stp0.Pos[0] = minAlong * stp1.Dir[0] + stp1.Pos[0];
      stp0.Pos[1] = minAlong * stp1.Dir[1] + stp1.Pos[1];
    }
    
    // Put the third shower Tj point on the shower axis at the maximum along position
    maxAlong = ss.Pts[ss.Pts.size()-1].RotPos[0];
    TrajPoint& stp2 = stj.Pts[2];
    stp2.Pos[0] = maxAlong * stp1.Dir[0] + stp1.Pos[0];
    stp2.Pos[1] = maxAlong * stp1.Dir[1] + stp1.Pos[1];
    
    // divide the longitudinal distance into 3 sections. Assign shower variables in these
    // sections to the 3 TPs
    sectionLength = (maxAlong - minAlong) / 3;
    sec0 = minAlong + sectionLength;
    sec2 = maxAlong - sectionLength;
    if(prt) mf::LogVerbatim("TC")<<" minAlong "<<(int)minAlong<<" maxAlong "<<(int)maxAlong<<" Section boundaries "<<(int)sec0<<" "<<(int)sec2<<" stp0 Pos "<<PrintPos(tjs, stp0.Pos)<<" stp2 Pos "<<PrintPos(tjs, stp2.Pos);
    
    // Calculate the charge and shower width in each section
    for(auto& sspt : ss.Pts) {
      unsigned short ipt = 1;
      if(sspt.RotPos[0] < sec0) ipt = 0;
      if(sspt.RotPos[0] > sec2) ipt = 2;
      stj.Pts[ipt].Chg += sspt.Chg;
      stj.Pts[ipt].DeltaRMS += sspt.Chg * sspt.RotPos[1] * sspt.RotPos[1];
      ++stj.Pts[ipt].NTPsFit;
    } // RotPt
    
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      TrajPoint& spt = stj.Pts[ipt];
      if(spt.NTPsFit < 2) continue;
      spt.DeltaRMS = sqrt(spt.DeltaRMS / spt.Chg);
      spt.AveChg = spt.Chg;
    } // ipt
    
    return true;

  } // DefineShowerTj
  
  ////////////////////////////////////////////////
  bool RefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Checks the properties of Shower Tj and revises them if necessary. Returns true if the
    // shower needs to be updated
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return false;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    if(stj.Pts.size() != 3) return false;
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return false;
    
    // Ignore fat showers
    if(ss.AspectRatio > 0.5) return false;
    // Ignore showers with an external parent
//    if(ss.ParentID > 0) return false;
    
    // check the beginning of the shower to see if the points are most on the shower axis
    float sum = 0;
    float sum2 = 0;
    float chg = 0;
    unsigned short cnt = 0;
    // Decide on the number of points to check = 1/2 of the points associated with the first shower Tj point
    unsigned short nchk = stj.Pts[0].NTPsFit / 2;
    if(nchk < 5) return false;
    for(auto& sspt : ss.Pts) {
      chg += sspt.Chg;
      sum += sspt.Chg * sspt.RotPos[1];
      sum2 += sspt.Chg * sspt.RotPos[1] * sspt.RotPos[1];
      ++cnt;
      if(cnt == nchk) break;
    } // sspt
    if(chg == 0) return false;
    float transAve = sum / chg;
    float arg = sum2 - chg * transAve * transAve;
    if(arg == 0) return false;
    float transRMS = sqrt(arg / chg);
    transRMS /= sqrt((float)cnt);
    
    if(std::abs(transAve)/transRMS < 2) return false;
    // tweak the angle
    ss.Angle += transAve / ss.Pts[0].RotPos[0];
    if(prt) mf::LogVerbatim("TC")<<"RSTj shower begin transAve "<<transAve<<" transRMS "<<transRMS<<" New shower angle "<<ss.Angle;
    
    return true;
  } // RefineShowerTj
      
  ////////////////////////////////////////////////
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<unsigned short>>& tjList)
  {
    // Tag Tjs with PDGCode = 11 if they have MCSMom < ShowerTag[0] and there are more than
    // ShowerTag[6] other Tjs with a separation < ShowerTag[1]. Returns a list of Tjs that meet this criteria
    
    tjList.clear();
    
    short maxMCSMom = tjs.ShowerTag[1];
    unsigned short minCnt = tjs.ShowerTag[7];
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      Trajectory& tj1 = tjs.allTraj[it1];
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) continue;
      tj1.NNeighbors = 0;
      // identified as a parent
      // ignore shower Tjs
      if(tj1.AlgMod[kShowerTj]) continue;
      // and Tjs that are already in showers
      if(tj1.AlgMod[kInShower]) continue;
      // ignore muons
      if(tj1.PDGCode == 13) continue;
      // ignore stubby Tjs
      if(tj1.Pts.size() < 3) continue;
      // Cut on length and MCSMom
      if(tj1.Pts.size() > 10 && tj1.MCSMom > maxMCSMom) continue;
      if(TjHasNiceVtx(tjs, tj1)) continue;
      tj1.PDGCode = 0;
      std::vector<unsigned short> list;
      for(unsigned short it2 = 0; it2 < tjs.allTraj.size(); ++it2) {
        if(it1 == it2) continue;
        Trajectory& tj2 = tjs.allTraj[it2];
        if(tj2.CTP != inCTP) continue;
        if(tj2.AlgMod[kKilled]) continue;
        // identified as a parent
        // ignore shower Tjs
        if(tj2.AlgMod[kShowerTj]) continue;
        // and Tjs that are already in showers
        if(tj2.AlgMod[kInShower]) continue;
        // ignore muons
        if(tj2.PDGCode == 13) continue;
        // ignore stubby Tjs
        if(tj2.Pts.size() < 3) continue;
        // Cut on length and MCSMom
        if(tj2.Pts.size() > 10 && tj2.MCSMom > maxMCSMom) continue;
        if(TjHasNiceVtx(tjs, tj2)) continue;
        unsigned short ipt1, ipt2;
        float doca = tjs.ShowerTag[2];
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, doca);
        if(doca < tjs.ShowerTag[2]) {
          // start the list with the ID of tj1
          if(list.empty()) list.push_back(tj1.ID);
          list.push_back(tj2.ID);
          ++tj1.NNeighbors;
        }
      } // it2
      if(list.size() > minCnt) tjList.push_back(list);
    } // it1
    
  } // TagShowerTjs
  
  ////////////////////////////////////////////////
  void AddMissedTjs(TjStuff& tjs, const CTP_t& inCTP, std::vector<unsigned short>& tjl)
  {
    // Just what it says
    if(tjl.empty()) return;
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      Trajectory& tj1 = tjs.allTraj[it1];
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) continue;
      // identified as a parent
      // ignore shower Tjs
      if(tj1.AlgMod[kShowerTj]) continue;
      // and Tjs that are already in showers
      if(tj1.AlgMod[kInShower]) continue;
      // Cut on length
      if(tj1.Pts.size() > 10) continue;
      if(TjHasNiceVtx(tjs, tj1)) continue;
      // already included in tjl?
      if(std::find(tjl.begin(), tjl.end(), tj1.ID) != tjl.end()) continue;
      // check proximity to Tjs in tjl
      unsigned short ipt1, ipt2;
      for(auto& tjID : tjl) {
        Trajectory& tjInShower = tjs.allTraj[tjID - 1];
        float doca = tjs.ShowerTag[2];
        TrajTrajDOCA(tjs, tj1, tjInShower, ipt1, ipt2, doca);
        if(doca < tjs.ShowerTag[2]) {
          tjl.push_back(tj1.ID);
          tj1.AlgMod[kInShower] = true;
          break;
        } // its a keeper
      } // tjid
    } // it1
    
  } // AddMissedTjs
  
  ////////////////////////////////////////////////
  void DefineEnvelope(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    ss.Envelope.resize(4);
    if(ss.TjIDs.empty()) return;
    
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
      myprt<<"DE: "<<cotIndex<<" Envelope";
      for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      myprt<<" Area "<<(int)ss.EnvelopeArea;
      myprt<<" ChgDensity "<<ss.ChgDensity;
    }
    
  } // DefineEnvelope  
  
  ////////////////////////////////////////////////
  bool AddTjsInsideEnvelope(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // This function adds Tjs to the shower. It updates the shower parameters.
    
    if(cotIndex > tjs.cots.size() - 1) return false;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.Envelope.empty()) return false;
    if(ss.TjIDs.empty()) return false;
    
     unsigned short nadd = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kShowerTj]) continue;
      // This shouldn't be necessary but do it for now
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) continue;
      // See if both ends are outside the envelope
      bool end0Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[0]].Pos, ss.Envelope);
      bool end1Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[1]].Pos, ss.Envelope);
      if(!end0Inside && !end1Inside) continue;
      // at least one end is inside. See if both are inside
      if(end0Inside && end1Inside) {
        // Fully contained
        // TODO: See if the Tj direction is compatible with the shower?
        if(AddTj(tjs, tj.ID, cotIndex, false, prt)) ++nadd;
        ++nadd;
//        if(prt) mf::LogVerbatim("TC")<<" Add contained Tj "<<tj.ID;
        continue;
      } // both ends inside
      // Require high momentum Tjs be aligned with the shower axis
      // TODO also require high momentum Tjs close to the shower axis?
      if(tj.MCSMom > 500) {
        float tjAngle = tj.Pts[tj.EndPt[0]].Ang;
        float dangPull = std::abs(tjAngle -ss.AngleErr) / ss.AngleErr;
        if(dangPull > 2) continue;
      } // high momentum
      if(AddTj(tjs, tj.ID, cotIndex, false, prt)) ++nadd;
    } // tj
    
    if(nadd > 0) {
      if(prt) mf::LogVerbatim("TC")<<"ATIE:  Added "<<nadd<<" trajectories. Calling UpdateShower... ";
      UpdateShower(tjs, cotIndex, prt);
      return true;
    } else {
      if(prt) mf::LogVerbatim("TC")<<"ATIE:  No new trajectories added to envelope ";
      return false;
    }
        
  } // AddTjsInsideEnvelope
  
  ////////////////////////////////////////////////
  bool AddLooseHits(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Add hits that are inside the envelope to the shower if they are loose, i.e. not
    // used by any trajectory. This function returns true if the set of hits is different than
    // the current set. The calling function should update the shower if this is the case.
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.Envelope.empty()) return false;
    if(ss.TjIDs.empty()) return false;

    if(prt) mf::LogVerbatim("TC")<<"ALH: AddLooseHits ";

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
        point[0] = tjs.fHits[iht].WireID.Wire;
        point[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
        if(!PointInsideEnvelope(point, ss.Envelope)) continue;
        newHits.push_back(iht);
      } // iht
    } // wire
    
    // no new hits and no old hits. Nothing to do
    if(newHits.empty()) {
      if(prt) mf::LogVerbatim("TC")<<" No new loose hits found";
      return false;
    }
    
    // Update
    stp0.Hits.insert(stp0.Hits.end(), newHits.begin(), newHits.end());
    for(auto& iht: newHits) tjs.fHits[iht].InTraj = stj.ID;
    
    if(prt) mf::LogVerbatim("TC")<<" Added "<<stp0.Hits.size()<<" hits to stj "<<stj.ID;
    return true;

  } // AddLooseHits
  
  
  ////////////////////////////////////////////////
  void FindStartChg(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Finds the charge at the start of a shower
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    ss.StartChg = 0;
    ss.StartChgErr = 1;
    float minAlong = ss.Pts[0].RotPos[0];
    
    // define a temp vector to sum the charge in bins of 4 WSE units ~ 1.2 cm in uB
    std::vector<float> schg(100);
    for(auto& sspt : ss.Pts) {
      unsigned short indx = (unsigned short)((sspt.RotPos[0] - minAlong)/4);
//      if(prt) mf::LogVerbatim("TC")<<indx<<"  "<<sspt.RotPos[1]<<" "<<(int)sspt.Chg;
      if(indx > 99) break;
      // Count the charge if it is within 4 WSE transverse
      if(std::abs(sspt.RotPos[1]) < 4) schg[indx] += sspt.Chg;
    }
    // find the first index with the first non-zero charge entry
    // TODO this needs to be done much better
    unsigned short ii = 0;
    for(ii = 0; ii < schg.size(); ++ii) if(schg[ii] > 0) break;
    if(ii > schg.size() / 2) return;
    ss.StartPt = ii;
    
    // average the first 3 bins
    float loChg = schg[ii];
    float hiChg = schg[ii];
    for(ii = ss.StartPt; ii < ss.StartPt + 3; ++ii) {
      ss.StartChg += schg[ii];
      if(schg[ii] < loChg) loChg = schg[ii];
      if(schg[ii] > hiChg) hiChg = schg[ii];
    }
    // Calculate the charge in 1 WSE unit of distance. We averaged over 3 bins each
    // of which is 4 WSE in size
    ss.StartChg /= 12;
    ss.StartChgErr = (hiChg - loChg) / 12;
    // Then correct for the path length in the wire cell
    TrajPoint& stp0 = tjs.allTraj[ss.ShowerTjID-1].Pts[0];
    float path = 50;
    if(stp0.Dir[0] != 0) path = 1 / std::abs(stp0.Dir[0]);
    if(path > 50) path = 50;
    ss.StartChg *= path;
    ss.StartChgErr *= path;
    
//    ss.Chg *= path;
//    ss.ChgErr *= path;
    
    if(prt) mf::LogVerbatim("TC")<<"FSC: cotIndex "<<cotIndex<<" path length "<<path<<" Starting charge "<<(int)ss.StartChg<<" +/- "<<(int)ss.StartChgErr<<" StartPt  "<<ss.StartPt<<" at pos "<<PrintPos(tjs, ss.Pts[ss.StartPt].Pos);
    
  } // FindStartChg
  
  ////////////////////////////////////////////////
  void TransferTjHits(TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Transfer InShower hits to the shower Tj 
    
    for(unsigned short ish = 0; ish < tjs.cots.size(); ++ish) {
      ShowerStruct& ss = tjs.cots[ish];
      // Ensure that this is the correct CTP
      if(ss.CTP != inCTP) continue;
      // Ensure that it is valid
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
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) {
          std::cout<<"TTjH: Coding error. Tj "<<tjID<<" is a ShowerTj but is in TjIDs\n";
          continue;
        }
        if(!tjs.allTraj[itj].AlgMod[kInShower]) {
          std::cout<<"TTjH: Coding error. Trying to transfer Tj "<<tjID<<" hits but it isn't an InShower Tj\n";
          continue;
        }
        auto thits = PutTrajHitsInVector(tjs.allTraj[itj], kUsedHits);
        stj.Pts[1].Hits.insert(stj.Pts[1].Hits.end(), thits.begin(), thits.end());
        // kill Tjs that are in showers
        tjs.allTraj[itj].AlgMod[kKilled] = true;
      } //  tjID
      // re-assign the hit -> stj association
      for(auto& iht : stj.Pts[1].Hits) tjs.fHits[iht].InTraj = stj.ID;
    } // ish
  } // TransferTjHits

  ////////////////////////////////////////////////
  void CollectLooseHits(TjStuff& tjs, const CTP_t& inCTP, bool prt)
  {
    // Collect hits in the vicinity of the shower
    
    geo::PlaneID planeID = DecodeCTP(inCTP);
    unsigned short ipl = planeID.Plane;
    
    for(unsigned short ish = 0; ish < tjs.cots.size(); ++ish) {
      ShowerStruct& ss = tjs.cots[ish];
      // Ensure that this is the correct CTP
      if(ss.CTP != inCTP) continue;
      // Ensure that it is valid
      if(ss.TjIDs.empty()) continue;
      if(ss.ShowerTjID == 0) continue;
      if(ss.ShowerTjID == 0) continue;
      // Tp 1 of stj has all of the shower hits
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      // look for other hits inside the envelope. Find the range of wires that spans the envelope
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
          point[0] = tjs.fHits[iht].WireID.Wire;
          point[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
          if(!PointInsideEnvelope(point, ss.Envelope)) continue;
          stj.Pts[0].Hits.push_back(iht);
          tjs.fHits[iht].InTraj = stj.ID;
        } // iht
      } // wire
    } // ish
  } // CollectLooseHits

  ////////////////////////////////////////////////
  unsigned short ShowerTjCotsIndex(TjStuff& tjs, const unsigned short& ShowerTjID)
  {
    for(unsigned short ii = 0; ii < tjs.cots.size(); ++ii) {
      if(ShowerTjID == tjs.cots[ii].ShowerTjID) return ii;
    } // iii
    return USHRT_MAX;
    
  } // ShowerTjCotsIndex

  ////////////////////////////////////////////////
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss)
  {
    if(ss.TjIDs.empty()) return 0;
    if(ss.ShowerTjID == 0) return 0;
    
    // Conversion from shower charge to energy in MeV. 0.0143 comes from an eye-bal fit.
    // Divide by the expected shower containment of 90%. This needs to be calculated directly
    constexpr float fShMeVPerChg = 0.0143 / 0.9;
    
    const Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    return fShMeVPerChg * stj.AveChg;
    
  } // ShowerEnergy

}