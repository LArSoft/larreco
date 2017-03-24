#include "larreco/RecoAlg/TCAlg/TCShower.h"

namespace tca {
  
  bool isBetter(const ShowerParentStruct& x, const ShowerParentStruct& y) { return x.FOM < y.FOM; }

  ////////////////////////////////////////////////
  void MakeShowers(TjStuff& tjs, const calo::CalorimetryAlg& fCaloAlg)
  {
    // Fill 3D shower variables
    
    if(tjs.ShowerTag[0] < 0) return;
    
    // Get the calibration constants
    
    bool prt = (tjs.ShowerTag[8] >= 0);
    
    int shID = 0;
    for(unsigned short ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned short imv = tjs.matchVecPFPList[ipfp];
      auto& msp = tjs.matchVec[imv];
      if(msp.TjIDs.empty()) continue;
      // look for a shower parent
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"MS: "<<ipfp<<" Parent "<<msp.Parent<<" PDGCode "<<msp.PDGCode<<" Dtr size "<<msp.DtrIndices.size();
        myprt<<" TjIDs:";
        for(auto& tjID : msp.TjIDs) myprt<<" "<<tjID;
      }
      
      if(msp.PDGCode != 11) continue;
      if(msp.DtrIndices.empty()) continue;
      if(prt) mf::LogVerbatim("TC")<<" Shower parent. Dtr matchVecPFPList index "<<msp.DtrIndices[0];
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
          SpacePtDir(tjs, itp, jtp, dir, dirErr);
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
//        double wirePitch = tjs.geom->WirePitch(planeID);
        double angleToVert = tjs.geom->WireAngleToVertical(tjs.geom->View(planeID), planeID.TPC, planeID.Cryostat) - 0.5 * ::util::pi<>();
        double cosgamma = std::abs(std::sin(angleToVert) * ss3.Dir.Y() + std::cos(angleToVert) * ss3.Dir.Z());
        if(cosgamma == 0) continue;
        double dx = tjs.geom->WirePitch(planeID) / cosgamma;
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
//        std::cout<<"ID "<<tj.ID<<" CTP "<<tj.CTP<<" pitch "<<wirePitch<<" angleToVert "<<angleToVert<<" dQ "<<dQ<<" dx "<<dx;
//        std::cout<<" dEdx "<<ss3.dEdx[planeID.Plane];
        // Estimate the error
        dQ += dQErr;
        ss3.dEdxErr[planeID.Plane] = fCaloAlg.dEdx_AREA(dQ, time, dx, planeID.Plane);
        ss3.dEdxErr[planeID.Plane] -= ss3.dEdx[planeID.Plane];
//        std::cout<<" dEdxErr "<<ss3.dEdxErr[planeID.Plane]<<"\n";
      } // ii
      // Note that we shouldn't define the ss3 Hit vector until hit merging is done
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
    if(tjs.ShowerTag[8] >= 0) {
      geo::PlaneID planeID = DecodeCTP(inCTP);
      CTP_t printCTP = EncodeCTP(planeID.Cryostat, planeID.TPC, std::nearbyint(tjs.ShowerTag[8]));
      prt = (printCTP == inCTP);
    }
    
    std::vector<std::vector<unsigned short>> tjList;
    TagShowerTjs(tjs, inCTP, tjList);
    if(prt) std::cout<<"Inside FindShowers inCTP "<<inCTP<<" tjList size "<<tjList.size()<<"\n";
    if(tjs.ShowerTag[0] == 1) return;
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
    
    // Add Tjs that were missed
    for(auto& tjl : tjList) {
      if(tjl.empty()) continue;
      AddMissedTjs(tjs, inCTP, tjl);
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
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"Make cots "<<tjs.cots.size()<<" in CTP "<<ss.CTP<<" TjID_NN";
        for(auto& tjID : tjl) myprt<<" "<<tjID<<"_"<<tjs.allTraj[tjID-1].NNeighbors;
      }
      unsigned short cotIndex = tjs.cots.size() - 1;
      FindShowerCenter(tjs, cotIndex, prt);
      // check for a failure
      if(tjs.cots[cotIndex].TjIDs.empty()) continue;
      FindShowerAxis(tjs, cotIndex, prt);
      // Define the shower Tj using the shower axis
      DefineShowerTj(tjs, cotIndex, prt);
      DefineEnvelope(tjs, cotIndex, prt);
      AddTjsInsideEnvelope(tjs, cotIndex, prt);
      // look for a parent that starts outside the envelope and ends inside the envelope
      FindParent(tjs, cotIndex, prt);
      if(prt) PrintTrajectory("FS", tjs, tjs.allTraj[stj.ID-1], USHRT_MAX);
    } // tjl
    
    if(tjs.cots.empty()) return;
    
    // merge showers?
    if(tjs.ShowerTag[0] > 2) MergeShowers(tjs, inCTP, prt);
    
    // drop those that don't meet the requirements
    for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
      ShowerStruct& ss = tjs.cots[ic];
      if(ss.CTP != inCTP) continue;
      if(ss.TjIDs.empty()) continue;
      // enough Tjs?
      unsigned short ntjs = ss.TjIDs.size();
      if(!ss.Parent.empty()) ++ntjs;
      bool killit = (ntjs < tjs.ShowerTag[7]);
      if(prt) mf::LogVerbatim("TC")<<"ic "<<ic<<" nTjs "<<ss.TjIDs.size()<<" killit? "<<killit;
      unsigned short nTjWithVtx = 0;
      // require a parent Traj
      if(!killit) killit = (ss.Parent.empty());
      // with a good FOM
      if(!killit) killit = (ss.Parent[0].FOM > 1);
      // Kill runt showers
      if(!killit) killit = (ss.Energy < tjs.ShowerTag[3]);
      if(!killit) killit = (ss.ChgDensity < 0.5);
//      if(!killit) killit = (ss.EnvelopeAspectRatio > 2);
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
        if(nTjPts < tjs.ShowerTag[6]) killit = true;
        if(prt) mf::LogVerbatim("TC")<<"    "<<" nTjPts "<<nTjPts<<" killit? "<<killit;
      } // !killit
      if(killit) {
        // Unset shower and killed bits. Trajectories that are in showers haven't had their hits re-assigned to the
        // shower Tj yet so nothing needs to be done to them
        for(auto& tjID : ss.TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          tj.AlgMod[kInShower] = false;
          tj.AlgMod[kKilled] = false;
        } // tjID
        // kill the shower Tj
        ss.TjIDs.clear();
        unsigned short itj = ss.ShowerTjID - 1;
        MakeTrajectoryObsolete(tjs, itj);
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
    
    // Finish up in this CTP. 
    // Re-assign hits from the InShower Tjs to the ShowerTj.
    TransferTjHits(tjs, inCTP, prt);
    // Assign unused hits inside the envelope to the ShowerTjs
    CollectLooseHits(tjs, inCTP, prt);
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
      if(ss.Parent.empty()) {
        std::cout<<"No parent Tj found for shower Tj "<<ss.ShowerTjID<<". Just saying...\n";
        unsigned short itj = ss.ShowerTjID - 1;
        PrintTrajectory("FS", tjs, tjs.allTraj[itj], USHRT_MAX);
        continue;
      }
      Trajectory& ptj = tjs.allTraj[ss.Parent[0].ID - 1];
      if(ptj.CTP != ss.CTP) {
        std::cout<<"FindShowers: Bad parent CTP "<<ss.CTP<<" "<<ptj.CTP<<" parent ID "<<ptj.ID<<"\n";
      }
    } // ss
    
    if(prt) {
      for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
        if(tjs.cots[ic].TjIDs.empty()) continue;
        unsigned short itj = tjs.cots[ic].ShowerTjID - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(prt) PrintTrajectory("FS", tjs, tj, USHRT_MAX);
      } // ic
    }
    
  } // FindShowers
  
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
      if(tj1.AlgMod[kShowerParent]) continue;
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
      for(auto& tjid : tjl) {
        Trajectory& tj2 = tjs.allTraj[tjid - 1];
        float doca = tjs.ShowerTag[2];
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, doca);
        if(doca < tjs.ShowerTag[2]) {
          tjl.push_back(tj2.ID);
          tj2.AlgMod[kInShower] = true;
        } // its a keeper
      } // tjid
    } // it1
    
  } // AddMissedTjs
  
  ////////////////////////////////////////////////
  void FindParent(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    // look for a parent trajectory for the cluster of trajectories, cotIndex. This should have the
    // signature
    //                *   This represents the shower charge center = Point 1 of the shower trajectory
    //   -----------      The next 3 lines represent the parent trajectory that is well reconstructed
    //               \    before it enters the shower, but wanders after shower hits are added to it
    
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
    
    // Reference the Tp charge center of the shower Tj
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
    
    if(prt) mf::LogVerbatim("TC")<<"FSP: cotIndex "<<cotIndex<<" Parent candidate size "<<ss.Parent.size();
    
    if(ss.Parent.empty()) {
      std::cout<<"FSP: No parent candidates for cotIndex "<<cotIndex<<" have been found. Assume this is bad\n";
      ss.TjIDs.clear();
      return;
    }
    
    float bestFOM = 20;
    unsigned short imTheBest = USHRT_MAX;
    unsigned short imTheBestEnd = 0;
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      // ignore Tjs that are already in showers or are shower Tjs
      // Candidate parents have already been identified
      if(tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kShowerTj]) continue;
      // find the point that is farthest from stp1
      float sep0 = PosSep2(tj.Pts[tj.EndPt[0]].Pos, stp1.Pos);
      float sep1 = PosSep2(tj.Pts[tj.EndPt[1]].Pos, stp1.Pos);
      unsigned short useEnd = 0;
      if(sep1 > sep0) useEnd = 1;
      unsigned short endPt = tj.EndPt[useEnd];
      float fom = ParentFOM(tjs, tj, endPt, ss, prt);
      if(fom > 30) continue;
      if(prt) mf::LogVerbatim("TC")<<"FSP: Tj_end "<<tj.ID<<"_"<<useEnd<<" fom "<<fom;
      if(fom > bestFOM) continue;
      bestFOM = fom;
      imTheBest = tj.ID;
      imTheBestEnd = useEnd;
    } // tj
    
    if(prt) mf::LogVerbatim("TC")<<"FSP: Best Tj_end "<<imTheBest<<"_"<<imTheBestEnd<<" fom "<<bestFOM;
    
    if(imTheBest == USHRT_MAX) return;
    
    // remove any old parent definition
    for(auto& parentCandidate : ss.Parent) {
      unsigned short itj = parentCandidate.ID - 1;
      tjs.allTraj[itj].AlgMod[kShowerParent] = false;
    }
    // Add it to the list of parent candidates
    ShowerParentStruct sps;
    sps.ID = imTheBest;
    sps.End = imTheBestEnd;
    sps.FOM = bestFOM;
    ss.Parent.push_back(sps);
/*
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Parent candidate list\n";
      for(auto& sps : ss.Parent) myprt<<" "<<sps.ID<<" "<<sps.FOM<<"\n";
    }
*/
    // sort by increasing FOM
    std::sort(ss.Parent.begin(), ss.Parent.end(), isBetter);
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Parent candidate list sorted\n";
      for(auto& sps : ss.Parent) myprt<<" "<<sps.ID<<" "<<sps.FOM<<"\n";
    }
    
    // See if the parent is an InShower Tj and if so change its status
    unsigned short iptj = ss.Parent[0].ID - 1;
    Trajectory& ptj = tjs.allTraj[iptj];
    if(ptj.AlgMod[kInShower]) {
      // ensure that it is associated with this shower
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), ss.Parent[0].ID) == ss.TjIDs.end()) {
        std::cout<<"FindParent: Tj "<<ss.Parent[0].ID<<" is an InShower Tj but it is not used in shower with cotIndex "<<cotIndex<<". Can't make it a parent\n";
        ss.Parent.clear();
        return;
      } // Parent is an InShower TJ but not in the correct shower
      ptj.AlgMod[kInShower] = false;
      ptj.AlgMod[kShowerParent] = true;
    } // Parent is an InShower Tj
    
  } //FindParent
  
  ////////////////////////////////////////////////
  float ParentFOM(TjStuff& tjs, Trajectory& tj, const unsigned short& tjPt, ShowerStruct& ss, bool prt)
  {
    // returns a FOM for the trajectory at point tjPt being the parent of ss
    
    if(tjPt > tj.Pts.size() - 1) return 1000;
    if(ss.Energy == 0) return 1000;
    
    // Radiation length converted to WSE units (for uB)
    constexpr float radLen = 14 / 0.3;
    constexpr float tenRadLen2 = 100 * radLen * radLen;

    // expected error on the angle of parent TP and the shower axis
    constexpr float dangRMS = 0.3;
    // Expected rms if the impact parameter between the projected parent and the shower center
    constexpr float deltaRMS = 9;

    if(ss.TjIDs.empty()) return 1000;
    if(ss.ShowerTjID == 0) return 1000;
    
    // prospective parent TP
    TrajPoint& ptp = tj.Pts[tjPt];
    // Shower charge center TP
    TrajPoint& stp1 = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
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
    float deltaPull = delta / deltaRMS;
    float dang = DeltaAngle(ptp.Ang, stp1.Ang);
    float dangPull = dang / dangRMS;
    float fom = 0.33 * sqrt(sepPull * sepPull + deltaPull * deltaPull + dangPull * dangPull);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"PFOM: Tj "<<tj.ID<<" Pos "<<PrintPos(tjs, ptp)<<" Energy "<<(int)ss.Energy;
      myprt<<std::fixed<<std::setprecision(2);
      myprt<<" tp1Sep "<<tp1Sep<<" sepPull "<<sepPull;
      myprt<<" delta "<<delta<<" deltaPull "<<deltaPull;
      myprt<<" dang "<<dang<<" dangPull "<<dangPull;
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
            FindParent(tjs, ict, prt);
            DefineShowerTj(tjs, ict, prt);
            DefineEnvelope(tjs, ict, prt);
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
            FindParent(tjs, jct, prt);
            DefineShowerTj(tjs, jct, prt);
            DefineEnvelope(tjs, jct, prt);
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
    if(prt) mf::LogVerbatim("TC")<<"FSA: Pos "<<PrintPos(tjs, stp1)<<" Angle "<<ss.Angle<<" +/- "<<ss.AngleErr;

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
    
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      if(itj > tjs.allTraj.size() - 1) {
        std::cout<<"Bad TjID "<<ss.TjIDs[it]<<"\n";
        ss.TjIDs.clear();
        return;
      }
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        stp1.Chg += tp.Chg;
        stp1.Pos[0] += tp.Chg * tp.Pos[0];
        stp1.Pos[1] += tp.Chg * tp.Pos[1];
      } // ipt
    } // it
    
    if(stp1.Chg <= 0) {
      std::cout<<"FSC: Crazy charge for cotIndex "<<cotIndex<<"\n";
      ss.TjIDs.clear();
      return;
    }
    stp1.Pos[0] /= stp1.Chg;
    stp1.Pos[1] /= stp1.Chg;
    ss.Energy = ShowerEnergy(tjs, ss);
    if(prt) mf::LogVerbatim("TC")<<"FSC: "<<cotIndex<<" Pos "<<PrintPos(tjs, stp1.Pos)<<" Chg "<<(int)stp1.Chg<<" Energy "<<(int)ss.Energy<<" MeV";
    
  } // FindShowerCenter
  
  ////////////////////////////////////////////////
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss)
  {
    if(ss.TjIDs.empty()) return 0;
    if(ss.ShowerTjID == 0) return 0;
    
    // Conversion from shower charge to energy in MeV. 0.0143 comes from an eye-bal fit.
    // Divide by the expected shower containment of 90%. This needs to be calculated directly
    constexpr float fShMeVPerChg = 0.0143 / 0.9;
    
    const Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    float totChg = (stj.Pts[0].Chg + stj.Pts[1].Chg + stj.Pts[2].Chg);
    // Add the energy of the parent
    if(!ss.Parent.empty()) {
      unsigned short iptj = ss.Parent[0].ID - 1;
      if(tjs.allTraj[iptj].AlgMod[kShowerParent]) {
        for(auto& tp : tjs.allTraj[iptj].Pts) totChg += tp.Chg;
      } // first in the list is defined to be the parent
    } // parent candidates exist
    return fShMeVPerChg * totChg;
    
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
      unsigned short endPt = tjs.allTraj[ptj].EndPt[end];
      showerAng = tjs.allTraj[ptj].Pts[endPt].Ang;
      showerAngErr = 0.2;
    }
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
    
    // Struct of all in-shower Tj points rotated into the shower coordinate system
    struct rotStruct {
      std::array<float, 2> Pos;
      float Chg;
      unsigned short TjID;
      unsigned short TjPt;
      unsigned short TjEnd;
    };
    
    std::vector<rotStruct> rotPos;
    
    TrajPoint& stp0 = stj.Pts[0];
    TrajPoint& stp1 = stj.Pts[1];
    TrajPoint& stp2 = stj.Pts[2];
    rotStruct rs;
    
    // Keep track of the index of the rotPos vector for the point which is lowest and highest
    // One of the Tjs with these points might be the parent
    unsigned short minAlongRPIndex = 0;
    unsigned short maxAlongRPIndex = 0;
    for(unsigned short it = 0; it < ss.TjIDs.size(); ++it) {
      unsigned short itj = ss.TjIDs[it] - 1;
      Trajectory& tj = tjs.allTraj[itj];
      unsigned short middlePt = 0.5 * (tj.EndPt[0] + tj.EndPt[1]);
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        TrajPoint& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        // Position of this point relative to stp1
        rs.Pos[0] = tp.Pos[0] - stp1.Pos[0];
        rs.Pos[1] = tp.Pos[1] - stp1.Pos[1];
        // Rotated into the stp1 direction
        float along = cs * rs.Pos[0] - sn * rs.Pos[1];
        float trans = sn * rs.Pos[0] + cs * rs.Pos[1];
        rs.Pos[0] = along;
        rs.Pos[1] = std::abs(trans);
        rs.Chg = tp.Chg;
        rs.TjID = tj.ID;
        rs.TjPt = ipt;
        rs.TjEnd = 0;
        if(ipt > middlePt) rs.TjEnd = 1;
        if(along < 0) {
          // along < 0 stj Pts[0]
          if(along < minAlong) {
            minAlong = along;
            minAlongRPIndex = rotPos.size();
            stp0.Pos = tp.Pos;
          }
        }  else {
          // along > 0 stj Pts[2]
          if(along > maxAlong) {
            maxAlong = along;
            maxAlongRPIndex = rotPos.size();
            stp2.Pos = tp.Pos;
          }
        } // along > 0
        rotPos.push_back(rs);
      } // ipt
    } // it
    
    if(minAlong == 0 || maxAlong == 0) return;
    
    // Place stp0 and stp2 on the shower axis
    stp0.Pos[0] = minAlong * stp1.Dir[0] + stp1.Pos[0];
    stp0.Pos[1] = minAlong * stp1.Dir[1] + stp1.Pos[1];
    stp2.Pos[0] = maxAlong * stp1.Dir[0] + stp1.Pos[0];
    stp2.Pos[1] = maxAlong * stp1.Dir[1] + stp1.Pos[1];
    
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
    
    for(auto& rp : rotPos) {
      unsigned short ipt = 1;
      if(rp.Pos[0] < sec0) ipt = 0;
      if(rp.Pos[0] > sec2) ipt = 2;
      TrajPoint& spt = stj.Pts[ipt];
      spt.Chg += rp.Pos[2];
      if(rp.Pos[1] > spt.Delta) spt.Delta = rp.Pos[1];
      spt.DeltaRMS += rp.Pos[2] * rp.Pos[1] * rp.Pos[1];
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
    for(auto& rp : rotPos) {
      along += std::abs(rp.Pos[0]);
      trans += std::abs(rp.Pos[1]);
    }
    ss.AspectRatio = trans / along;
    
    // push the min/max Tj IDs onto the Parent vector
    if(ss.Parent.empty()) {
      struct ShowerParentStruct sps;
      sps.ID = rotPos[minAlongRPIndex].TjID;
      sps.End = rotPos[minAlongRPIndex].TjEnd;
      unsigned short endPt = tjs.allTraj[sps.ID - 1].EndPt[sps.End];
      sps.FOM = ParentFOM(tjs, tjs.allTraj[sps.ID-1], endPt, ss, prt);
      ss.Parent.push_back(sps);
      sps.ID = rotPos[maxAlongRPIndex].TjID;
      sps.End = rotPos[maxAlongRPIndex].TjEnd;
      endPt = tjs.allTraj[sps.ID - 1].EndPt[sps.End];
      sps.FOM = ParentFOM(tjs, tjs.allTraj[sps.ID-1], endPt, ss, prt);
      ss.Parent.push_back(sps);
    } // parent empty

    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"DSTj cotIndex "<<cotIndex<<" showerAng "<<showerAng;
      myprt<<" stj Positions "<<PrintPos(tjs, stp0.Pos)<<" "<<PrintPos(tjs, stp1.Pos)<<" "<<PrintPos(tjs, stp2.Pos);
      myprt<<" AspectRatio "<<ss.AspectRatio;
    };
    
  } // DefineShowerTj
  
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
      myprt<<"DefineEnvelope "<<cotIndex;
      for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      myprt<<" ChgDensity "<<ss.ChgDensity;
    }
    
  } // DefineEnvelope  
  
  ////////////////////////////////////////////////
  void AddTjsInsideEnvelope(TjStuff& tjs, const unsigned short& cotIndex, bool prt)
  {
    
    if(cotIndex > tjs.cots.size() - 1) return;
    
    ShowerStruct& ss = tjs.cots[cotIndex];
    if(ss.Envelope.empty()) return;
    if(ss.TjIDs.empty()) return;
     
    Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
    
    // Get the charge before any more Tjs are added
    float showerChg = stj.Pts[0].Chg + stj.Pts[1].Chg + stj.Pts[2].Chg;
    float addedChg = 0;
    if(prt) mf::LogVerbatim("TC")<<"ATIE: cotIndex "<<cotIndex<<" shower charge "<<(int)showerChg;
    
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != ss.CTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      if(tj.AlgMod[kInShower]) continue;
      if(tj.AlgMod[kShowerTj]) continue;
      // This shouldn't be necessary but do it for now
      if(std::find(ss.TjIDs.begin(), ss.TjIDs.end(), tj.ID) != ss.TjIDs.end()) {
        std::cout<<"AddTjsInsideEnvelope: Tj "<<tj.ID<<" is already inside this envelope "<<cotIndex<<"\n";
        continue;
      }
      // See if both ends are outside the envelope
      bool end0Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[0]].Pos, ss.Envelope);
      bool end1Inside = PointInsideEnvelope(tj.Pts[tj.EndPt[1]].Pos, ss.Envelope);
      if(!end0Inside || !end1Inside) continue;
      // at least one end is inside. See if both are inside
      if(end0Inside || end1Inside) {
        // TODO: See if the Tj direction is compatible with the shower
        // both ends are inside. Add it to the shower
        ss.TjIDs.push_back(tj.ID);
        tjs.allTraj[tj.ID-1].AlgMod[kInShower] = true;
        // Count the charge added
        for(auto& tp : tj.Pts) addedChg += tp.Chg;
        if(prt) mf::LogVerbatim("TC")<<" Add contained Tj "<<tj.ID<<" addedChg sum "<<addedChg;
        continue;
      } // both ends inside
      // One end is outside. Count the charge inside and the total charge
      float tjChg = 0;
      float tjChgInside = 0;
      for(auto& tp : tj.Pts) {
        if(tp.Chg == 0) continue;
        tjChg += tp.Chg;
        if(PointInsideEnvelope(tp.Pos, ss.Envelope)) tjChgInside += tp.Chg;
      } // tp
      if(tjChg == 0) continue;
      float insideFrac = tjChgInside / tjChg;
      float tjFrac = tjChg / showerChg;
      // Ignore it if the total charge is small fraction of the shower charge and < 20% of the charge is inside.
      if(prt) mf::LogVerbatim("TC")<<" Partially contained Tj "<<tj.ID<<" tjChg "<<(int)tjChg<<" tjChgInside "<<(int)tjChgInside<<" insideFrac "<<insideFrac<<" tj charge fraction "<<tjFrac;
      if(tjFrac < 0.2 && insideFrac < 0.2) continue;
      ss.TjIDs.push_back(tj.ID);
      tjs.allTraj[tj.ID-1].AlgMod[kInShower] = true;
      addedChg += tjChg;
    } // tj
    
    if(prt) mf::LogVerbatim("TC")<<" Beginning charge "<<(int)showerChg<<" "<<(int)addedChg;
    // decide whether to update the 
    
  } // AddTjsInsideEnvelope
  
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
        std::cout<<"CollectHits: ShowerTj "<<stj.ID<<" already has hits. This can't be right\n";
        continue;
      }
      // make stj a daughter if there is a parent
      if(!ss.Parent.empty()) {
        unsigned short parentTjID = ss.Parent[0].ID;
        stj.ParentTrajID = parentTjID;
        tjs.allTraj[parentTjID - 1].PDGCode = 11;
      }
      stj.PDGCode = 11;
      // Note that UseHit is not used since the size is limited to 16
      for(auto& tjID : ss.TjIDs) {
        unsigned short itj = tjID - 1;
        if(tjs.allTraj[itj].AlgMod[kShowerTj]) {
          std::cout<<"CollectHits: Coding error. Tj "<<tjID<<" is a ShowerTj but is in TjIDs\n";
          continue;
        }
        auto thits = PutTrajHitsInVector(tjs.allTraj[itj], kUsedHits);
        stj.Pts[1].Hits.insert(stj.Pts[1].Hits.end(), thits.begin(), thits.end());
        // kill Tjs that are in showers
        // Set the InShower bit to indicate the reason that it was killed
        tjs.allTraj[itj].AlgMod[kInShower] = true;
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
          stj.Pts[1].Hits.push_back(iht);
          tjs.fHits[iht].InTraj = stj.ID;
        } // iht
      } // wire
    } // ish
  } // CollectLooseHits

  
  /////////////////////////////////////////
  void SpacePtDir(TjStuff& tjs, TrajPoint itp, TrajPoint jtp, TVector3& dir, TVector3& dirErr)
  {
    
    dir.SetX(-999);

    if(itp.CTP == jtp.CTP) return;

    TVector3 pt1, pt2;
    geo::PlaneID iplnID = DecodeCTP(itp.CTP);
    geo::PlaneID jplnID = DecodeCTP(jtp.CTP);
    
    double y, z;
    double xi  = tjs.detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iplnID);
    //    std::cout<<"xi "<<xi;
    double xj = tjs.detprop->ConvertTicksToX(jtp.Pos[1] / tjs.UnitsPerTick, jplnID);
    //    std::cout<<" xj "<<xj<<"\n";
    // don't continue if the points are too far apart in X
    if(std::abs(xi - xj) > 5) return;
    xi = 0.5 * (xi + xj);
    
    unsigned int wire1 = (unsigned int)(itp.Pos[0] + 0.5);
    unsigned int wire2 = (unsigned int)(jtp.Pos[0] + 0.5);
/*
    geo::WireID wireID1(iplnID, wire1), wireID2(jplnID, wire2);
    bool exists = tjs.geom->WireIDIntersection(wireID1, wireID2, pt1);
    if(!exists) return;
    geom->Plane(iplnID)->DriftPoint(-xi);
*/
    
    //    std::cout<<"wire1 "<<iplnID.Plane<<":"<<wire1<<" wire2 "<<jplnID.Plane<<":"<<wire2;
    tjs.geom->IntersectionPoint(wire1, wire2, iplnID.Plane, jplnID.Plane, iplnID.Cryostat, iplnID.TPC, y, z);
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
    xi  = tjs.detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iplnID);
    //    std::cout<<" new x "<<xi<<"\n";
    //    std::cout<<" jtp Dir "<<jtp.Dir[0]<<" "<<jtp.Dir[1]<<"\n";
    // Determine the number of wires to move jtp to get to the same X position
    std::array<float, 2> newPos;
    newPos[1] = tjs.detprop->ConvertXToTicks(xi, jplnID) * tjs.UnitsPerTick;
    //    std::cout<<" jtp.Pos[1] "<<jtp.Pos[1]<<" newPos[1] "<<newPos[1]<<" UPT "<<tjs.UnitsPerTick<<"\n";
    newPos[0] = (newPos[1] - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    //    std::cout<<" newPos "<<PrintPos(tjs, newPos)<<"\n";
    wire2 = (unsigned int)(newPos[0] + 0.5);
    tjs.geom->IntersectionPoint(wire1, wire2, iplnID.Plane, jplnID.Plane, iplnID.Cryostat, iplnID.TPC, y, z);
    pt2.SetX(xi);
    pt2.SetY(y);
    pt2.SetZ(z);
    dir = pt2 - pt1;
    dir.SetMag(1);
  } // MakeSpacePt

}