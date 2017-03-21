#include "larreco/RecoAlg/TCAlg/TCShower.h"

namespace tca {
  
  
  ////////////////////////////////////////////////
  void MakeShowers(TjStuff& tjs, const std::vector<float>& fShowerTag, const art::ServiceHandle<geo::Geometry>& geom, const detinfo::DetectorProperties*& detprop, const calo::CalorimetryAlg& fCaloAlg)
  {
    // Fill shower variables
    
    if(fShowerTag[0] < 0) return;
    
    // Get the calibration constants
    
    int shID = 0;
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned short imv = tjs.matchVecPFPList[ipfp];
      auto& msp = tjs.matchVec[imv];
      if(msp.TjIDs.empty()) continue;
      // look for a shower parent
      if(true) {
        mf::LogVerbatim myprt("TC");
        myprt<<"MS: "<<ipfp<<" Parent "<<msp.Parent<<" PDGCode "<<msp.PDGCode<<" Dtr size "<<msp.DtrIndices.size();
        myprt<<" TjIDs:";
        for(auto& tjID : msp.TjIDs) myprt<<" "<<tjID;
      }
      
      if(msp.PDGCode != 11) continue;
      if(msp.DtrIndices.empty()) continue;
      mf::LogVerbatim("TC")<<" Shower parent. Dtr matchVecPFPList index "<<msp.DtrIndices[0];
      ++shID;
      ShowerStruct3D ss3;
      ss3.Energy.resize(tjs.NumPlanes);
      ss3.EnergyErr.resize(tjs.NumPlanes);
      ss3.MIPEnergy.resize(tjs.NumPlanes);
      ss3.MIPEnergyErr.resize(tjs.NumPlanes);
      ss3.dEdx.resize(tjs.NumPlanes);
      ss3.dEdxErr.resize(tjs.NumPlanes);
      ss3.ID = shID;
      // Start the list of parent and daughter Tjs in all planes
      ss3.TjIDs = msp.TjIDs;
      // Find the daughter matchVec
      unsigned short dtrIndex = tjs.matchVecPFPList[msp.DtrIndices[0]];
      auto& msd = tjs.matchVec[dtrIndex];
      // Add the daughter trajectory IDs
      ss3.TjIDs.insert(ss3.TjIDs.end(), msd.TjIDs.begin(), msd.TjIDs.end());
      // fill the start position
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) ss3.Pos[ixyz] = msp.sXYZ[ixyz];
      // cycle through pairs of Tj IDs to get the direction
      unsigned short npair = 0;
      for(unsigned short ii = 0; ii < msp.TjIDs.size() - 1; ++ii) {
        unsigned short itj = msp.TjIDs[ii] - 1;
        for(unsigned short jj = ii + 1; jj < msp.TjIDs.size(); ++jj) {
          unsigned short jtj = msp.TjIDs[jj] - 1;
          if(tjs.allTraj[jtj].CTP == tjs.allTraj[itj].CTP) continue;
          TrajPoint& itp = tjs.allTraj[itj].Pts[tjs.allTraj[itj].EndPt[0]];
          TrajPoint& jtp = tjs.allTraj[jtj].Pts[tjs.allTraj[jtj].EndPt[0]];
          PrintHeader("itp");
          PrintTrajPoint("itp", tjs, tjs.allTraj[itj].EndPt[0], 0, 0, itp);
          PrintTrajPoint("jtp", tjs, tjs.allTraj[jtj].EndPt[0], 0, 0, jtp);
          TVector3 dir, dirErr;
          SpacePtDir(tjs, geom, detprop, itp, jtp, dir, dirErr);
          if(dir.X() < -10) continue;
          ss3.Dir += dir;
          ss3.DirErr += dirErr;
          ++npair;
          if(npair == 2) {
            ss3.Dir *= 0.5;
            ss3.DirErr *= 0.5;
            break;
          }
        } // jj
        if(npair == 2) break;
      } // ii
      // Calculate the shower length. Use the start XYZ position of the parent
      // and the XYZ position of the Shower Tj daughter that is furthest away. 
      // assume the far position is the end
      std::array<float, 3> farPos = msd.eXYZ;
      // and correct it if it's not
      if(std::abs(msd.sXYZ[2] - msp.sXYZ[2]) > std::abs(msd.eXYZ[2] - msp.sXYZ[2])) farPos = msd.sXYZ;
      ss3.Len = 0;
      for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) {
        double arg = farPos[ixyz] - msp.sXYZ[ixyz];
        ss3.Len += arg * arg;
      } // ixyz
      ss3.Len = sqrt(ss3.Len);
      mf::LogVerbatim("TC")<<"Pos "<<ss3.Pos.X()<<" "<<ss3.Pos.Y()<<" "<<ss3.Pos.Z()<<" Dir "<<ss3.Dir.X()<<" "<<ss3.Dir.Y()<<" "<<ss3.Dir.Z()<<" Len "<<ss3.Len;
      // Calculate the opening angle
      // Calculate dEdx, etc
      for(unsigned short ii = 0; ii < msp.TjIDs.size(); ++ii) {
        unsigned short itj = msp.TjIDs[ii] - 1;
        Trajectory& tj = tjs.allTraj[itj];
        geo::PlaneID planeID = DecodeCTP(tj.CTP);
        double wirePitch = geom->WirePitch(planeID);
        double angleToVert = geom->WireAngleToVertical(geom->View(planeID), planeID.TPC, planeID.Cryostat) - 0.5 * ::util::pi<>();
        double cosgamma = std::abs(std::sin(angleToVert) * ss3.Dir.Y() + std::cos(angleToVert) * ss3.Dir.Z());
        if(cosgamma == 0) continue;
        double dx = geom->WirePitch(planeID) / cosgamma;
        // get the average charge from the first 4 points
        double dQ = 0;
        double dQErr = 0;
        double cnt = 0;
        // Don't use the first point
        for(unsigned short ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1]; ++ipt) {
          if(tj.Pts[ipt].Chg == 0) continue;
          ++cnt;
          dQ += tj.Pts[ipt].Chg;
          dQErr += tj.Pts[ipt].Chg * tj.Pts[ipt].Chg;
          if(cnt == 5) break;
        }
        dQ /= cnt;
        double arg = dQErr - cnt * dQ * dQ;
        if(arg < 0) arg = 0;
        dQErr = sqrt(arg / (cnt - 1));
        double time = tj.Pts[tj.EndPt[0]].Pos[1] / tjs.UnitsPerTick;
        ss3.dEdx[planeID.Plane] = fCaloAlg.dEdx_AREA(dQ, time, dx, planeID.Plane);
        std::cout<<"ID "<<tj.ID<<" CTP "<<tj.CTP<<" pitch "<<wirePitch<<" angleToVert "<<angleToVert<<" dQ "<<dQ<<" dx "<<dx;
        std::cout<<" dEdx "<<ss3.dEdx[planeID.Plane];
        // Estimate the error
        dQ += dQErr;
        ss3.dEdxErr[planeID.Plane] = fCaloAlg.dEdx_AREA(dQ, time, dx, planeID.Plane);
        ss3.dEdxErr[planeID.Plane] -= ss3.dEdx[planeID.Plane];
        std::cout<<" dEdxErr "<<ss3.dEdxErr[planeID.Plane]<<"\n";
      } // ii
      // Note that we shouldn't define the ss3 Hit vector until hit merging is done
      tjs.showers.push_back(ss3);
    } // ipfp
  } // MakeShowers

  ////////////////////////////////////////////////
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag, std::vector<std::vector<unsigned short>>& tjList)
  {
    // Tag Tjs with PDGCode = 12 if they have MCSMom < fShowerTag[0] and there are more than
    // fShowerTag[6] other Tjs with a separation < fShowerTag[1]. Returns a list of Tjs that meet this criteria
    
    tjList.clear();
    
    short maxMCSMom = fShowerTag[1];
    unsigned short minCnt = fShowerTag[7];
    
    for(unsigned short it1 = 0; it1 < tjs.allTraj.size(); ++it1) {
      Trajectory& tj1 = tjs.allTraj[it1];
      if(tj1.CTP != inCTP) continue;
      if(tj1.AlgMod[kKilled]) continue;
      tj1.NNeighbors = 0;
      // identified as a parent
      if(tj1.AlgMod[kShowerParent]) continue;
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
        if(tj2.AlgMod[kShowerParent]) continue;
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
        float doca = fShowerTag[2];
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, doca);
        if(doca < fShowerTag[2]) {
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
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag)
  {
    // Construct clusters of trajectories (cots) which will become shower PFParticles
    
    // fShowerTag[] parameters
    // 0 Mode (<= 0 OFF, 1 = tag only, 2 = find showers)
    // 1 Max Tj MCSMom for a shower tag (< 0 = no shower-like Tj tagging or shower finding)
    // 2 Max separation
    // 3 Max delta angle
    // 4 rms width factor
    // 5 Min shower 1/2 width (WSE units)
    // 6 Min total Tj Pts
    // 7 Min Tjs
    // 8 Debug in CTP
    
    if(fShowerTag[0] <= 0) return;
    
    bool prt = false;
    if(fShowerTag[8] >= 0) {
      geo::PlaneID planeID = DecodeCTP(inCTP);
      CTP_t printCTP = EncodeCTP(planeID.Cryostat, planeID.TPC, std::nearbyint(fShowerTag[8]));
      prt = (printCTP == inCTP);
    }
    
    std::vector<std::vector<unsigned short>> tjList;
    TagShowerTjs(tjs, inCTP, fShowerTag, tjList);
    if(prt) std::cout<<"Inside FindShowers inCTP "<<inCTP<<" tjList size "<<tjList.size()<<"\n";
    if(fShowerTag[0] == 1) return;
    if(tjList.empty()) return;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"TagShower tjlist\n";
      for(auto& tjl : tjList) {
        if(tjl.empty()) continue;
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
        myprt<<"\n";
      } // tjl
    } // printCTP
    
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
    unsigned short minNeighbors = fShowerTag[7];
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
      // put it in TJ stuff. The rest of the info will be added in DefineShower
      tjs.cots.push_back(ss);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Make cots "<<tjs.cots.size()<<" using";
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
      }
      unsigned short cotIndex = tjs.cots.size() - 1;
      FindShowerCenter(tjs, cotIndex, prt);
      FindShowerAxis(tjs, cotIndex, prt);
      // Define the shower Tj using the shower axis
      DefineShowerTj(tjs, cotIndex, prt);
      // look for a parent
      FindShowerParent(tjs, cotIndex, fShowerTag, prt);
      // re-define the shower if a parent was found
//      if(tjs.cots[cotIndex].ParentTrajID > 0) DefineShowerTj(tjs, cotIndex, prt);
      DefineShowerEnvelope(tjs, cotIndex, fShowerTag, prt);
      if(prt) PrintTrajectory("FS", tjs, tjs.allTraj[stj.ID-1], USHRT_MAX);
    } // tjl
    
    if(tjs.cots.empty()) return;
    
    // merge showers?
    if(fShowerTag[0] > 2) MergeShowers(tjs, inCTP, fShowerTag, prt);
    
    // drop those that don't meet the requirements
    for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
      ShowerStruct& ss = tjs.cots[ic];
      if(ss.CTP != inCTP) continue;
      if(ss.TjIDs.empty()) continue;
      // enough Tjs?
      unsigned short ntjs = ss.TjIDs.size();
      if(!ss.Parent.empty()) ++ntjs;
      bool killit = (ntjs < fShowerTag[7]);
      if(prt) mf::LogVerbatim("TC")<<"ic "<<ic<<" nTjs "<<ss.TjIDs.size()<<" killit? "<<killit;
      unsigned short nTjWithVtx = 0;
      // require a parent Traj
      if(!killit) killit = (ss.Parent.empty());
      // with a good FOM
      if(!killit) killit = (ss.Parent[0].FOM > 1);
      // Kill runt showers
      if(!killit) killit = (ShowerEnergy(tjs, ss) < fShowerTag[3]);
      if(!killit) killit = (ss.ChgDensity < 0.5);
      if(!killit) killit = (ss.EnvelopeAspectRatio > 2);
      if(!killit) {
        // count the number of Tj points
        unsigned short nTjPts = 0;
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          nTjPts += NumPtsWithCharge(tjs, tj, false);
          if(tj.VtxID[0] > 0 || tj.VtxID[1] > 0) ++nTjWithVtx;
        }  // tjID
        // add the parent point count
        if(!ss.Parent.empty()) nTjPts += NumPtsWithCharge(tjs, tjs.allTraj[ss.Parent[0].ID - 1], false);
        if(nTjPts < fShowerTag[6]) killit = true;
        if(prt) mf::LogVerbatim("TC")<<"    "<<" nTjPts "<<nTjPts<<" killit? "<<killit;
      } // !killit
      if(killit) {
        // unset shower and killed bits
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          tj.AlgMod[kInShower] = false;
          tj.AlgMod[kKilled] = false;
        } // tjID
        // kill the shower Tj
        ss.TjIDs.clear();
        unsigned short itj = ss.ShowerTjID - 1;
        MakeTrajectoryObsolete(tjs, itj);
        // Trajectories that are in showers haven't had their hits re-assigned to the
        // shower Tj yet so nothing needs to be done to them
      }
      // Set the PDGCode of the parent to shower-like
      if(!killit && !ss.Parent.empty()) {
        unsigned short ptj = ss.Parent[0].ID - 1;
        tjs.allTraj[ptj].AlgMod[kShowerParent] = true;
        tjs.allTraj[ptj].PDGCode = 11;
        std::cout<<"FindShowers parent "<<ss.Parent[0].ID<<"\n";
      }
      // kill vertices in the showers that are left
      if(!killit && nTjWithVtx > 0) {
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] > 0) MakeVertexObsolete(tjs, tj.VtxID[end]);
          } // end
          tj.AlgMod[kInShower] = true;
        } // tjID
      } // !killit
    } // ic
    
    // Finish up in this CTP
    CollectHits(tjs, inCTP, prt);
    
    // check for consistency
    for(auto& ss : tjs.cots) {
      if(ss.TjIDs.empty()) continue;
      if(ss.CTP != inCTP) continue;
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
      if(ss.Parent.empty() == 0) {
        std::cout<<"No parent Tj found for shower. Just saying...\n";
        continue;
      }
      Trajectory& ptj = tjs.allTraj[ss.Parent[0].ID - 1];
      if(ptj.CTP != ss.CTP) {
        std::cout<<"FindShowers: Bad parent CTP "<<ss.CTP<<" "<<ptj.CTP<<" parent ID "<<ptj.ID<<"\n";
      }
    } // ss
    
    if(fShowerTag[8] >= 0) {
      for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
        if(tjs.cots[ic].TjIDs.empty()) continue;
        unsigned short itj = tjs.cots[ic].ShowerTjID - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(prt) PrintTrajectory("FS", tjs, tj, USHRT_MAX);
      } // ic
    }
    
  } // FindShowers
  
  ////////////////////////////////////////////////
  void FindShowerParent(TjStuff& tjs, const unsigned short& cotIndex, const std::vector<float>& fShowerTag, bool prt)
  {
    // look for a parent trajectory for the cluster of trajectories, cotIndex. This should have the
    // signature
    //                *   This represents the shower charge center = Point 1 of the shower trajectory
    //   -----------      The next 3 lines represent the parent trajectory that is well reconstructed
    //               \    before it enters the shower, but wanders after shower hits are added to it
    
    // fShowerTag[] parameters
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
    
    // Reference the Tp charge center of the shower Tj
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
    unsigned short onTj = USHRT_MAX;
    unsigned short atEnd = 0;
    float maxSep = 0;
    
    // Find the TP on an InShower Tj that is farthest away from the charge center
    TrajPoint farTP = stp1;
    for(auto& tjid : ss.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjid - 1];
      for(unsigned short end = 0; end < 2; ++end) {
        unsigned short endPt = tj.EndPt[end];
        float sep2 = PosSep2(tj.Pts[endPt].Pos, stp1.Pos);
        if(sep2 > maxSep) {
          maxSep = sep2;
          onTj = tjid;
          atEnd = end;
          farTP = tj.Pts[endPt];
        }
      } // end
    } // tjid
    mf::LogVerbatim("TC")<<"FSP maxSep point "<<PrintPos(tjs, farTP.Pos)<<" onTj "<<onTj<<" atEnd "<<atEnd;
    
  } //FindShowerParent
/*
  ////////////////////////////////////////////////
  void FindShowerParent(TjStuff& tjs, const unsigned short& cotIndex, const std::vector<float>& fShowerTag, bool prt)
  {
    // look for a parent trajectory for the cluster of trajectories, cotIndex. This should have the
    // signature
    //                *   This represents the shower charge center = Point 1 of the shower trajectory
    //   -----------      The next 3 lines represent the parent trajectory that is well reconstructed
    //               \    before it enters the shower, but wanders after shower hits are added to it
    
    // fShowerTag[] parameters
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
    
    unsigned short minParentLength = 3;
    float maxDelta = 2 * fShowerTag[2];
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    // clobber any previous Shower -> Parent assignment
    if(ss.ParentTrajID != 0) {
      unsigned short stj = ss.ShowerTjID - 1;
      tjs.allTraj[stj].ParentTrajID = 0;
    }
    ss.ParentTrajID = 0;
    
    // Ensure that it is valid
    if(ss.TjIDs.empty()) return;
    // Reference the Tp charge center of the shower Tj
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
    
    ss.ParentFOM = 5;    
    ss.FailedParentFOM = ss.ParentFOM;
    
    // A rough estimate of the radiation length in WSE units. A minimum separation between a trajectory and the
    // charge center should not be much larger than 5 radiation lengths. 
    constexpr float radLen = 14 / 0.3;
    // square it so we don't need to take a square root below to make the selection
    constexpr float radLen2Cut = 25 * radLen * radLen;
    
    constexpr float dangRMS = 0.3;
    // Expected rms if the impact parameter between the projected parent and the shower center
    constexpr float deltaRMS = 9;
    // Expect the parent Tj to have high MCSMom
//    constexpr float mcsmomAve = 300;
    // with generous errors
//    constexpr float mcsmomRMS = 200;
    // The shower max should be between 1 and 3 radiation lengths for electrons in the range of 0.1 - 1 GeV
    // Estimate the shower energy
    float shEnergy = ShowerEnergy(tjs, ss);
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
    float expectedParentTP1Sep = 0.85 * log(3 * shEnergy) - 2.65;
    // Convert it to the expected number of points
    expectedParentTP1Sep *= radLen;
    // Assume that the projection of the shower in this view will be ~2/3
    expectedParentTP1Sep *= 0.6;
    // Guess that the RMS of the separation will be ~30%
    float maxSepRMS = 0.5 * expectedParentTP1Sep;
//    constexpr float parentLengthAve = radLen;
//    constexpr float parentLengthRMS = radLen;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"FindShowerParent: Charge Center "<<PrintPos(tjs, stp1.Pos)<<" Shower energy estimate "<<(int)shEnergy;
    }
     // temp variables
//     float bdang = 0;
     float bdangp = 0;
     float bdelt = 0;
     float bdeltp = 0;
     float bmcsm = 0;
     float bmcsmp = 0;
     float blen = 0;
     float blenp = 0;
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != ss.CTP) continue;
      // it can't be a parent if it is a shower tj
      if(tj.AlgMod[kShowerTj]) continue;
      // or if it has low MCSMom TODO: Is this a good idea?
      //      if(tj.MCSMom < fShowerTag[1]) continue;
      // or if it is too short
      float npwc = NumPtsWithCharge(tjs, tj, true);
      if(npwc < minParentLength) continue;
      // determine which end is furthest from stp1
      unsigned short endPt0 = tj.EndPt[0];
      float sep0 = PosSep2(tj.Pts[endPt0].Pos, stp1.Pos);
      unsigned short endPt1 = tj.EndPt[1];
      float sep1 = PosSep2(tj.Pts[endPt1].Pos, stp1.Pos);
      // Use this end to check the DOCA with stp1. Use the other end to check the separation with stp1
      unsigned short farEnd = 0;
      float maxSep = sep0;
      float minSep = sep1;
      if(sep1 > sep0) {
        farEnd = 1;
        maxSep = sep1;
        minSep = sep0;
      }
//      if(prt) mf::LogVerbatim("TC")<<" Tj "<<tj.ID<<" minSep "<<minSep<<" radLen2Cut "<<radLen2Cut;
      if(minSep > radLen2Cut) continue;
      // re-purpose endPt0 to mean the end farthest away
      endPt0 = tj.EndPt[farEnd];
      TrajPoint& farTP = tj.Pts[endPt0];
      float delta = PointTrajDOCA(tjs, stp1.Pos[0], stp1.Pos[1], farTP);
      // Make a rough cut using the max separation cut
//      if(prt) mf::LogVerbatim("TC")<<"  delta "<<delta<<" maxDelta "<<maxDelta;
      if(delta > maxDelta) continue;
      float dang = DeltaAngle(farTP.Ang, stp1.Ang);
      float dangPull = dang / dangRMS;
      float deltaPull = delta / deltaRMS;
      maxSep = sqrt(maxSep);
      float maxSepPull = (maxSep - expectedParentTP1Sep) / maxSepRMS;
//      float farMom = tj.MCSMom;
//      float mcsmomPull = (farMom - mcsmomAve) / mcsmomRMS;
//      float lengthPull = (tj.Pts.size() - parentLengthAve) / parentLengthRMS;
      float fom = dangPull * dangPull + deltaPull * deltaPull + maxSepPull * maxSepPull;
//      fom += lengthPull * lengthPull;
      fom = 0.2 * sqrt(fom);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" tj.ID "<<tj.ID<<" farTP "<<PrintPos(tjs, farTP);
        myprt<<" expectedParentTP1Sep "<<(int)expectedParentTP1Sep;
        //        std::cout<<" minSep "<<minSep;
        myprt<<" dang  "<<std::fixed<<std::setprecision(2)<<dang<<" Pull  "<<dangPull;
        myprt<<" delta "<<delta<<" Pull "<<deltaPull;
        myprt<<" maxSep "<<maxSep<<" Pull "<<maxSepPull;
//        myprt<<" MCSMom "<<tj.MCSMom<<" farMom "<<(int)farMom<<" Pull "<<mcsmomPull;
        myprt<<" fom  "<<fom;
      }
      // check for a signal between these points
      TrajPoint& nearTP = tj.Pts[tj.EndPt[1 - farEnd]];
      if(!SignalBetween(tjs, nearTP, stp1, 0.2, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" No significant signal between "<<PrintPos(tjs, nearTP.Pos)<<" and "<<PrintPos(tjs, stp1.Pos);
        continue;
      }
      // keep track of the second best parent
      if(ss.ParentTrajID != 0) {
        if(fom < ss.ParentFOM) {
          // found a candidate parent that is worse than the previously found one
          ss.FailedParentTrajID = tj.ID;
          ss.FailedParentTrajEnd = farEnd;
          ss.FailedParentFOM = fom;
        } else {
          // found a candidate parent that is better than the previously found one
          ss.FailedParentTrajID = ss.ParentTrajID;
          ss.FailedParentTrajEnd = ss.ParentTrajEnd;
          ss.FailedParentFOM = ss.ParentFOM;
        }
      } // foundParent
      if(fom < ss.ParentFOM) {
        ss.ParentTrajID = tj.ID;
        ss.ParentTrajEnd = farEnd;
        ss.ParentFOM = fom;
        //        std::cout<<"PFOM1 "<<minSep<<" "<<dang<<" "<<delta<<" "<<maxSep<<" fom "<<fom<<"\n";

//         bdang = dang;
//         bdangp = dangPull;
//         bdelt = delta;
//         bdeltp = deltaPull;
//         bmcsm = farMom;
//         bmcsmp = mcsmomPull;
//         blen = tj.Pts.size();
//         blenp = lengthPull;

      }
    } // itj
    if(prt) mf::LogVerbatim("TC")<<"FSP: Set Parent ID "<<ss.ParentTrajID;
    
    if(ss.FailedParentTrajID == ss.ParentTrajID) ss.FailedParentTrajID = 0;
    
    // Look for the parent in the list of shower Tjs and remove it
    auto deleteMe = std::find(ss.TjIDs.begin(), ss.TjIDs.end(), ss.ParentTrajID);
    if(deleteMe != ss.TjIDs.end()) {
      ss.TjIDs.erase(deleteMe);
      // unset the InShower bit
      tjs.allTraj[ss.ParentTrajID - 1].AlgMod[kInShower] = false;
      FindShowerCenter(tjs, cotIndex, prt);
    }
    
//    if(bdang > 0) mf::LogVerbatim("TC")<<"FOM "<<ShowerEnergy(tjs, ss)<<" "<<bdang<<" "<<bdangp<<" "<<bdelt<<" "<<bdeltp<<" "<<bmcsm<<" "<<bmcsmp<<" "<<blen<<" "<<blenp<<" "<<ss.ParentFOM;
    
  } // FindShowerParent
 */
  
  ////////////////////////////////////////////////
  void MergeShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag, bool prt)
  {
    // Merge showers that point roughly in the same direction and aren't too far apart
    
    // fShowerTag[] parameters
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
    float maxSep = 2 * 14 / 0.3;
    
    bool didMerge = true;
    while(didMerge) {
      didMerge = false;
      for(unsigned short ict = 0; ict < tjs.cots.size() - 1; ++ict) {
        ShowerStruct& iss = tjs.cots[ict];
        if(iss.TjIDs.empty()) continue;
        if(iss.CTP != inCTP) continue;
        unsigned short itjIndex = iss.ShowerTjID - 1;
        Trajectory& itj = tjs.allTraj[itjIndex];
        float iChg = itj.Pts[0].Chg + itj.Pts[1].Chg + itj.Pts[2].Chg;
        for(unsigned short jct = ict + 1; jct < tjs.cots.size(); ++jct) {
          ShowerStruct& jss = tjs.cots[jct];
          if(jss.TjIDs.empty()) continue;
          if(jss.CTP != iss.CTP) continue;
          unsigned short jtjIndex = jss.ShowerTjID - 1;
          Trajectory& jtj = tjs.allTraj[jtjIndex];
          float jChg = jtj.Pts[0].Chg + jtj.Pts[1].Chg + jtj.Pts[2].Chg;
          //          if(prt) mf::LogVerbatim("TC")<<"MS "<<itj.ID<<" "<<jtj.ID;
          float sepi0j2 = PosSep(itj.Pts[0].Pos, jtj.Pts[2].Pos);
          float sepi2j0 = PosSep(itj.Pts[2].Pos, jtj.Pts[0].Pos);
          if(sepi0j2 > maxSep && sepi2j0 > maxSep) {
            //            if(prt) mf::LogVerbatim("TC")<<" Separation too large "<<sepi0j2<<" "<<sepi2j0<<" maxSep "<<maxSep;
            continue;
          }
          float delta;
          // find delta using the trajectory with the highest charge. The error on delta should be calculated
          // more carefully here. For now just assume that the error is the largest DeltaRMS of the higher charge shower
          float deltaErr = 0;
          // Also check to see if the direction of the higher
          // charge Tj is towards the lower charge Tj
          bool chgOK = false;
          if(iChg > jChg) {
            delta = PointTrajDOCA(tjs, jtj.Pts[1].Pos[0], jtj.Pts[1].Pos[1], itj.Pts[1]);
            TrajPoint itoj;
            MakeBareTrajPoint(tjs, itj.Pts[1], jtj.Pts[1], itoj);
            chgOK = (std::abs(itj.Pts[1].Ang - itoj.Ang) < M_PI / 2);
            for(auto& pt : itj.Pts) if(pt.DeltaRMS > deltaErr) deltaErr = pt.DeltaRMS;
          } else {
            delta = PointTrajDOCA(tjs, itj.Pts[1].Pos[0], itj.Pts[1].Pos[1], jtj.Pts[1]);
            TrajPoint jtoi;
            MakeBareTrajPoint(tjs, jtj.Pts[1], itj.Pts[1], jtoi);
            chgOK = (std::abs(jtj.Pts[1].Ang - jtoi.Ang) < M_PI / 2);
            for(auto& pt : jtj.Pts) if(pt.DeltaRMS > deltaErr) deltaErr = pt.DeltaRMS;
          }
          float dang = DeltaAngle(itj.Pts[1].Ang, jtj.Pts[1].Ang);
          // estimate the error on the relative angle between the showers
          float dangErr = sqrt(itj.Pts[1].Ang * itj.Pts[1].Ang + jtj.Pts[1].Ang * jtj.Pts[1].Ang);
          float dangSig = dang / dangErr;
          float deltaSig = delta / deltaErr;
          if(prt) {
            mf::LogVerbatim("TC")<<" merge candidates printed below:  dang "<<dang<<" dangSig "<<dangSig<<" delta "<<delta<<" deltaErr "<<deltaErr<<" deltaSig "<<deltaSig<<" chgOK? "<<chgOK;
            PrintTrajectory("itj", tjs, itj, USHRT_MAX);
            PrintTrajectory("jtj", tjs, jtj, USHRT_MAX);
          }
          // These cuts could use more care
          bool doMerge = dangSig < 2 && chgOK && deltaSig < 2;
          // merge if they have the same parent
          // deal with this later
//          if(iss.ParentTrajID == jss.ParentTrajID) doMerge = true;
          if(!doMerge) continue;
          if(prt) mf::LogVerbatim("TC")<<" Merge them";
          // Merge em. Put all the InShower Tjs in the higher charge shower Tj
          if(iChg > jChg) {
            iss.TjIDs.insert(iss.TjIDs.end(), jss.TjIDs.begin(), jss.TjIDs.end());
            FindShowerCenter(tjs, ict, prt);
            FindShowerParent(tjs, ict, fShowerTag, prt);
            DefineShowerTj(tjs, ict, prt);
            DefineShowerEnvelope(tjs, ict, fShowerTag, prt);
            jss.TjIDs.clear();
            // kill the shower Tj
            jtj.AlgMod[kKilled] = true;
            if(prt) {
              mf::LogVerbatim("TC")<<" merge done";
              PrintTrajectory("itj", tjs, itj, USHRT_MAX);
            }
          } else {
            jss.TjIDs.insert(jss.TjIDs.end(), iss.TjIDs.begin(), iss.TjIDs.end());
            FindShowerCenter(tjs, jct, prt);
            FindShowerParent(tjs, jct, fShowerTag, prt);
            DefineShowerTj(tjs, jct, prt);
            DefineShowerEnvelope(tjs, jct, fShowerTag, prt);
            iss.TjIDs.clear();
            itj.AlgMod[kKilled] = true;
            if(prt) {
              mf::LogVerbatim("TC")<<" merge done";
              PrintTrajectory("jtj", tjs, jtj, USHRT_MAX);
            }
          }
          didMerge = true;
        } // jct
        if(didMerge) break;
      } // ict
    } // didMerge
  } // MergeShowers
  
  ////////////////////////////////////////////////
  void FindShowerAxis(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Find the angle of the shower using the position of all of the TPs

    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    unsigned short stjIndex = ss.ShowerTjID - 1;
    if(stjIndex > tjs.allTraj.size() - 1) return;
    if(tjs.allTraj[stjIndex].Pts.size() != 3) return;

    // Do a least squares fit using all the points
    unsigned short cnt = 0;
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    double wt, xx, yy;
    
    TrajPoint& stp1 = tjs.allTraj[stjIndex].Pts[1];
    
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        wt = 1 / tp.Chg;
        cnt += 1;
        sum  += wt;
        xx = wt * tp.Pos[0] - stp1.Pos[0];
        yy = wt * tp.Pos[1] - stp1.Pos[1];
        sumx += wt * xx;
        sumy += wt * yy;
        sumx2 += wt * xx * xx;
        sumy2 += wt * yy * yy;
        sumxy += wt * xx * yy;
      } // ipt
    } // it
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0 || cnt < 3) {
      ss.Angle = 10;
      ss.AngleErr = 1.5;
      return;
    }
    // A is the intercept
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // B is the slope
    double B = (sumxy * sum  - sumx * sumy) / delta;
    double ndof = cnt - 2;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    double BErr = sqrt(varnce * sum / delta);
    ss.Angle = atan(B);
    ss.AngleErr = std::abs(atan(B + BErr) - ss.Angle);
    if(prt) mf::LogVerbatim("TC")<<"FindShowerAxis: Pos "<<PrintPos(tjs, stp1)<<" Angle "<<ss.Angle<<" Err "<<ss.AngleErr;

  } // FindShowerAxis

  ////////////////////////////////////////////////
  void FindShowerCenter(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Finds the charge center using all sub-structure trajectories in the cot. All of the shower
    // charge is assigned to the second TP. The charge will later be distributed between TP0 - TP2.
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    unsigned short stjIndex = ss.ShowerTjID - 1;
    if(stjIndex > tjs.allTraj.size() - 1) return;
    if(tjs.allTraj[stjIndex].Pts.size() != 3) return;
    
    TrajPoint& stp1 = tjs.allTraj[stjIndex].Pts[1];
    stp1.Chg = 0;
    stp1.Pos[0] = 0;
    stp1.Pos[1] = 0;
    
    unsigned short cnt = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        stp1.Chg += tp.Chg;
        stp1.Pos[0] += tp.Chg * tp.Pos[0];
        stp1.Pos[1] += tp.Chg * tp.Pos[1];
        ++cnt;
      } // ipt
    } // it
    
    if(stp1.Chg == 0) return;
    stp1.Pos[0] /= stp1.Chg;
    stp1.Pos[1] /= stp1.Chg;
    if(prt) mf::LogVerbatim("TC")<<"FindShowerCenter: Pos "<<PrintPos(tjs, stp1.Pos);
    
  } // FindShowerCenter
  
  ////////////////////////////////////////////////
  float ShowerEnergy(TjStuff& tjs, const ShowerStruct& ss)
  {
    if(ss.TjIDs.empty()) return 0;
    if(ss.ShowerTjID == 0) return 0;
    
    // Conversion from shower charge to energy in MeV. 0.0143 comes from an eye-bal fit.
    // Divide by the expected shower containment of 90%. This needs to be calculated directly
    constexpr float fShMeVPerChg = 0.0143 / 0.9;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    float nrg = fShMeVPerChg * (stj.Pts[0].Chg + stj.Pts[1].Chg + stj.Pts[2].Chg);
    return nrg;
    
  } // ShowerEnergy
  
  ////////////////////////////////////////////////
  void DefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // Finishes the definition of the shower. This also calculates the shower aspect ratio
     
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.TjIDs.empty()) return;
    
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    // ensure that the charge center is defined
    if(stj.Pts[1].Chg == 0) return;
    
    // decide which angle to use to define the shower Tj axis
    float showerAng = ss.Angle;
    float showerAngErr = ss.AngleErr;
    if(!ss.Parent.empty()) {
      // A parent is identified so use that angle instead
      // Use the first one for now
      unsigned short ptj = ss.Parent[0].ID - 1;
      unsigned short end = ss.Parent[0].End;
      showerAng = tjs.allTraj[ptj].Pts[end].Ang;
      showerAngErr = 0.2;
    }
    if(prt) mf::LogVerbatim("TC")<<"DefineShowerTj Parent Tj "<<ss.Parent[0].ID<<" showerAng "<<showerAng;
    float cs = cos(showerAng);
    float sn = sin(showerAng);
    // use this angle in all shower Tj TPs
    for(auto& stp : stj.Pts) {
      stp.Ang = showerAng;
      stp.AngErr = showerAngErr;
      stp.Dir[0] = cs;
      stp.Dir[1] = sn;
    } // stp
    
    // Determine the size of the shower along the axis and transverse to it. 
    // Rotate and translate each point into the coordinate system defined by tp[1]
    cs = cos(-showerAng);
    sn = sin(-showerAng);
    float minAlong = 0;
    float maxAlong = 0;
    std::array<float, 3> rotPos;
    // keep a copy of the rotated point and the charge of each to use later
    // rotPos (distance along, distance transverse, charge)
    std::vector<std::array<float, 3>> rotPts;
    
    TrajPoint& stp1 = stj.Pts[1];
    
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        // Position of this point relative to stp1
        rotPos[0] = tp.Pos[0] - stp1.Pos[0];
        rotPos[1] = tp.Pos[1] - stp1.Pos[1];
        // Rotated into the stp1 direction
        float along = cs * rotPos[0] - sn * rotPos[1];
        float trans = sn * rotPos[0] + cs * rotPos[1];
        rotPos[0] = along;
        rotPos[1] = std::abs(trans);
        rotPos[2] = tp.Chg;
        rotPts.push_back(rotPos);
        //        if(prt) std::cout<<PrintPos(tjs, tp)<<" along "<<along<<" trans "<<trans<<"\n";
        if(along < 0) {
          // along < 0 stj Pts[0]
          if(along < minAlong) {
            minAlong = along;
            stj.Pts[0].Pos = tp.Pos;
          }
        }  else {
          // along > 0 stj Pts[2]
          if(along > maxAlong) {
            maxAlong = along;
            stj.Pts[2].Pos = tp.Pos;
          }
        } // along > 0
      } // ipt
    } // it
    
    // Add the charge from the parent inside the shower
    if(!ss.Parent.empty()) {
      unsigned short itj = ss.Parent[0].ID - 1;
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        // Position of this point relative to stp1
        rotPos[0] = tp.Pos[0] - stp1.Pos[0];
        rotPos[1] = tp.Pos[1] - stp1.Pos[1];
        // Rotated into the stp1 direction
        float along = cs * rotPos[0] - sn * rotPos[1];
        if(along < minAlong || along > maxAlong) continue;
        float trans = sn * rotPos[0] + cs * rotPos[1];
        rotPos[0] = along;
        rotPos[1] = std::abs(trans);
        rotPos[2] = tp.Chg;
        rotPts.push_back(rotPos);
//        if(prt) std::cout<<PrintPos(tjs, tp)<<" along "<<along<<" trans "<<trans<<" Chg"<<(int)tp.Chg<<" \n";
        if(along < 0) {
          // along < 0 stj Pts[0]
          if(along < minAlong) {
            minAlong = along;
            stj.Pts[0].Pos = tp.Pos;
          }
        }  else {
          // along > 0 stj Pts[2]
          if(along > maxAlong) {
            maxAlong = along;
            stj.Pts[2].Pos = tp.Pos;
          }
        } // along > 0
      } // ipt
    } // Parent traj exists
    
    //    if(prt) std::cout<<"rotPts size "<<rotPts.size()<<"\n";
    
    
    // divide the longitudinal distance into 3 sections. Assign shower variables in these
    // sections to the 3 TPs
    float sectionLength = (maxAlong - minAlong) / 3;
    float sec0 = minAlong + sectionLength;
    float sec2 = maxAlong - sectionLength;
    // initialize the stj points
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      stj.Pts[ipt].Chg = 0;
      stj.Pts[ipt].Delta = 0;
      stj.Pts[ipt].DeltaRMS = 0;
      stj.Pts[ipt].NTPsFit = 0;
    } // ipt
    
    for(auto& rotPos : rotPts) {
      unsigned short ipt = 1;
      if(rotPos[0] < sec0) ipt = 0;
      if(rotPos[0] > sec2) ipt = 2;
      TrajPoint& spt = stj.Pts[ipt];
      spt.Chg += rotPos[2];
      if(rotPos[1] > spt.Delta) spt.Delta = rotPos[1];
      spt.DeltaRMS += rotPos[2] * rotPos[1] * rotPos[1];
      ++spt.NTPsFit;
    } // rotPos
    
    for(unsigned short ipt = 0; ipt < 3; ++ipt) {
      TrajPoint& spt = stj.Pts[ipt];
      if(spt.NTPsFit < 2) continue;
      spt.DeltaRMS = sqrt(spt.DeltaRMS / spt.Chg);
      spt.AveChg = spt.Chg;
    } // ipt
    
    // calculate the aspect ratio
    float along = 0;
    float trans = 0;
    for(auto& rotPos : rotPts) {
      along += rotPos[2] * std::abs(rotPos[0]);
      trans += rotPos[2] * std::abs(rotPos[1]);
    }
    ss.AspectRatio = trans / along;
    if(prt) mf::LogVerbatim("TC")<<"ShowerAspectRatio "<<ss.AspectRatio;

    // Here would be a good place to move Pts[0] and Pts[2] onto the shower axis
    
  } // DefineShowerTj
  
  ////////////////////////////////////////////////
  void DefineShowerEnvelope(TjStuff& tjs, const unsigned short& cotIndex, const std::vector<float>& fShowerTag, bool prt)
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
    ss.Envelope[0][0] = 1.1 * (stp0.Pos[0] - stp1.Pos[0]);
    ss.Envelope[0][1] = fShowerTag[5] + fShowerTag[4] * stp0.DeltaRMS;
    // second vertex
    ss.Envelope[1][0] = 1.1 * (stp2.Pos[0] - stp1.Pos[0]);
    ss.Envelope[1][1] = fShowerTag[5] + fShowerTag[4] * stp2.DeltaRMS;
    // third and fourth are reflections of the first and second
    ss.Envelope[2][0] =  ss.Envelope[1][0];
    ss.Envelope[2][1] = -ss.Envelope[1][1];
    ss.Envelope[3][0] =  ss.Envelope[0][0];
    ss.Envelope[3][1] = -ss.Envelope[0][1];
    
    // Find the aspect ratio
    
    ss.EnvelopeLength = ss.Envelope[1][0] - ss.Envelope[0][0];
    float width  = ss.Envelope[1][1] + ss.Envelope[0][1];
    ss.EnvelopeArea = ss.EnvelopeLength * width;
    ss.EnvelopeAspectRatio = width / ss.EnvelopeLength;
    
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
    
  } // DefineShowerEnvelope  
  
  ////////////////////////////////////////////////
  void CollectHits(TjStuff& tjs, const CTP_t& inCTP, bool prt)
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
      // Tp 1 of stj will get all of the shower hits
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      // this shouldn't be necessary but do it anyway
      ReleaseHits(tjs, stj);
      // make stj a daughter
/* deal with this
      if(ss.ParentTrajID != 0) {
        stj.ParentTrajID = ss.ParentTrajID;
        tjs.allTraj[ss.ParentTrajID - 1].PDGCode = 13;
      }
*/
      stj.PDGCode = 11;
      // Note that UseHit is not used since the size is limited.
      for(auto& tjID : ss.TjIDs) {
        unsigned short itj = tjID - 1;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) {
          std::cout<<"CollectHits: Coding error. Tj "<<tjID<<" is a ShowerTj but is in TjIDs\n";
          continue;
        }
        auto thits = PutTrajHitsInVector(tjs.allTraj[itj], kUsedHits);
        stj.Pts[1].Hits.insert(stj.Pts[1].Hits.end(), thits.begin(), thits.end());
        // kill Tjs that are in showers
        tjs.allTraj[itj].AlgMod[kInShower] = true;
        MakeTrajectoryObsolete(tjs, itj);
      } //  tjID
      // re-assign the hits to stj
      for(auto& iht : stj.Pts[1].Hits) tjs.fHits[iht].InTraj = stj.ID;
      // look for other hits inside the envelope
      float fLoWire = 1E6;
      float fHiWire = 0;
      for(auto& vtx : ss.Envelope) {
        if(vtx[0] < fLoWire) fLoWire = vtx[0];
        if(vtx[0] > fHiWire) fHiWire = vtx[0];
      } // vtx
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
          // already in the shower?
          if(std::find(stj.Pts[1].Hits.begin(), stj.Pts[1].Hits.end(), iht) != stj.Pts[1].Hits.end()) continue;
          // see if this hit is inside the envelope
          point[0] = tjs.fHits[iht].WireID.Wire;
          point[1] = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
          if(!PointInsideEnvelope(point, ss.Envelope)) continue;
          // Unused hit?
          if(tjs.fHits[iht].InTraj == 0) {
            // assign it to the shower tj
            stj.Pts[1].Hits.push_back(iht);
            tjs.fHits[iht].InTraj = stj.ID;
            continue;
          } // unused hit
          if(tjs.fHits[iht].InTraj < 0) continue;
          // Hit is inside the Envelope and is used in a Tj
          unsigned short itj = tjs.fHits[iht].InTraj - 1;
          Trajectory& tj = tjs.allTraj[itj];
          if(tj.AlgMod[kKilled]) continue;
          // correct CTP?
          if(tj.CTP != inCTP) continue;
          // Ignore muons
          if(tj.PDGCode == 13) continue;
          // ignore Tjs that are either a shower Tj or are in a shower
          if(tj.AlgMod[kShowerTj]) continue;
          if(tj.AlgMod[kInShower]) continue;
          // ignore Tjs with a nice vertex
          //          if(TjHasNiceVtx(tjs, tj)) continue;
          // See if the Tj is fully contained inside the envelope
          bool isContained = true;
          for(unsigned short end = 0; end < 2; ++end) {
            unsigned short endPt = tj.EndPt[end];
            if(!PointInsideEnvelope(tj.Pts[endPt].Pos, ss.Envelope)) {
              isContained = false;
              break;
            }
          } // end
          if(isContained) {
            // re-assign the Tj hits to the shower
            auto thits = PutTrajHitsInVector(tj, kUsedHits);
            MakeTrajectoryObsolete(tjs, itj);
            stj.Pts[1].Hits.insert(stj.Pts[1].Hits.end(), thits.begin(), thits.end());
            //            std::cout<<"Put contained Tj "<<tj.ID<<" hits in shower Tj "<<stj.ID<<"\n";
            for(auto& tht : thits) tjs.fHits[tht].InTraj = stj.ID;
            //            WatchHit("CH3", tjs, wHit, wInTraj, stj.ID);
          } // isContained
        } // iht
      } // wire
    } // ish
  } // CollectHits

  
  /////////////////////////////////////////
  void SpacePtDir(TjStuff& tjs, const art::ServiceHandle<geo::Geometry>& geom, const detinfo::DetectorProperties*& detprop, TrajPoint itp, TrajPoint jtp, TVector3& dir, TVector3& dirErr)
  {
    
    if(itp.CTP == jtp.CTP) {
      dir.SetX(-999);
      return;
    }
    TVector3 pt1, pt2;
    geo::PlaneID iplnID = DecodeCTP(itp.CTP);
    geo::PlaneID jplnID = DecodeCTP(jtp.CTP);
    //    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double y, z;
    double xi  = detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iplnID);
    //    std::cout<<"xi "<<xi;
    double xj = detprop->ConvertTicksToX(jtp.Pos[1] / tjs.UnitsPerTick, jplnID);
    //    std::cout<<" xj "<<xj<<"\n";
    // don't continue if the points are too far apart in X
    if(std::abs(xi - xj) > 5) {
      dir.SetX(-999);
      return;
    }
    xi = 0.5 * (xi + xj);
    
    unsigned int wire1 = (unsigned int)(itp.Pos[0] + 0.5);
    unsigned int wire2 = (unsigned int)(jtp.Pos[0] + 0.5);
    //    std::cout<<"wire1 "<<iplnID.Plane<<":"<<wire1<<" wire2 "<<jplnID.Plane<<":"<<wire2;
    geom->IntersectionPoint(wire1, wire2, iplnID.Plane, jplnID.Plane, iplnID.Cryostat, iplnID.TPC, y, z);
    pt1.SetX(xi);
    pt1.SetY(y);
    pt1.SetZ(z);
    //    std::cout<<" pt1.X "<<pt1.X()<<" pt1.Y "<<pt1.Y()<<" pt1.z "<<pt1.Z()<<"\n";
    
    // Move itp by 100 wires. It doesn't matter if we end up outside the TPC bounds since
    // bounds checking isn't done by these utility functions
    //    std::cout<<"itp "<<PrintPos(tjs, itp.Pos);
    MoveTPToWire(itp, itp.Pos[0] + 100);
    //    std::cout<<" -> "<<PrintPos(tjs, itp.Pos);
    wire1 = (unsigned int)(itp.Pos[0] + 0.5);
    xi  = detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iplnID);
    //    std::cout<<" new x "<<xi<<"\n";
    //    std::cout<<" jtp Dir "<<jtp.Dir[0]<<" "<<jtp.Dir[1]<<"\n";
    // Determine the number of wires to move jtp to get to the same X position
    std::array<float, 2> newPos;
    newPos[1] = detprop->ConvertXToTicks(xi, jplnID) * tjs.UnitsPerTick;
    //    std::cout<<" jtp.Pos[1] "<<jtp.Pos[1]<<" newPos[1] "<<newPos[1]<<" UPT "<<tjs.UnitsPerTick<<"\n";
    newPos[0] = (newPos[1] - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    //    std::cout<<" newPos "<<PrintPos(tjs, newPos)<<"\n";
    wire2 = (unsigned int)(newPos[0] + 0.5);
    geom->IntersectionPoint(wire1, wire2, iplnID.Plane, jplnID.Plane, iplnID.Cryostat, iplnID.TPC, y, z);
    pt2.SetX(xi);
    pt2.SetY(y);
    pt2.SetZ(z);
    dir = pt2 - pt1;
    dir.SetMag(1);
  } // MakeSpacePt

}