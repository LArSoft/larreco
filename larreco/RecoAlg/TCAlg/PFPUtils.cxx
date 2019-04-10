#include "larreco/RecoAlg/TCAlg/PFPUtils.h"
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

  /////////////////////////////////////////
  void StitchPFPs()
  {
    // Stitch PFParticles in different TPCs. This does serious damage to PFPStruct and should
    // only be called from TrajCluster module just before making PFParticles to put in the event
    if(slices.size() < 2) return;
    if(tcc.geom->NTPC() == 1) return;
    if(tcc.pfpStitchCuts.size() < 2) return;
    if(tcc.pfpStitchCuts[0] <= 0) return;
    
    bool prt = tcc.dbgStitch;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      std::string fcnLabel = "SP";
      myprt<<fcnLabel<<" cuts "<<sqrt(tcc.pfpStitchCuts[0])<<" "<<tcc.pfpStitchCuts[1]<<"\n";
      bool printHeader = true;
      for(size_t isl = 0; isl < slices.size(); ++isl) {
        if(debug.Slice >= 0 && int(isl) != debug.Slice) continue;
        auto& slc = slices[isl];
        if(slc.pfps.empty()) continue;
        for(auto& pfp : slc.pfps) PrintP(fcnLabel, myprt, pfp, printHeader);
      } // slc
    } // prt
    
    // lists of pfp UIDs to stitch
    std::vector<std::vector<int>> stLists;
    for(unsigned short sl1 = 0; sl1 < slices.size() - 1; ++sl1) {
      auto& slc1 = slices[sl1];
      for(unsigned short sl2 = sl1 + 1; sl2 < slices.size(); ++sl2) {
        auto& slc2 = slices[sl2];
        // look for PFParticles in the same recob::Slice
        if(slc1.ID != slc2.ID) continue;
        for(auto& p1 : slc1.pfps) {
          if(p1.ID <= 0) continue;
          // Can't stitch shower PFPs
          if(p1.PDGCode == 1111) continue;
          for(auto& p2 : slc2.pfps) {
            if(p2.ID <= 0) continue;
            // Can't stitch shower PFPs
            if(p2.PDGCode == 1111) continue;
            float maxSep2 = tcc.pfpStitchCuts[0];
            float maxCth = tcc.pfpStitchCuts[1];
            bool gotit = false;
            for(unsigned short e1 = 0; e1 < 2; ++e1) {
              auto pos1 = PosAtEnd(p1, e1);
              // require the end to be close to a TPC boundary
              if(InsideFV(slc1, p1, e1)) continue;
              auto dir1 = DirAtEnd(p1, e1);
              for(unsigned short e2 = 0; e2 < 2; ++e2) {
                auto pos2 = PosAtEnd(p2, e2);
                // require the end to be close to a TPC boundary
                if(InsideFV(slc2, p2, e2)) continue;
                auto dir2 = DirAtEnd(p2, e2);
                float sep = PosSep2(pos1, pos2);
                if(sep > maxSep2) continue;
                float cth = std::abs(DotProd(dir1, dir2));
                if(cth < maxCth) continue;
                maxSep2 = sep;
                maxCth = cth;
                gotit = true;
              } // e2
            } // e1
            if(!gotit) continue;
            if(prt) {
              mf::LogVerbatim myprt("TC");
              myprt<<"Stitch slice "<<slc1.ID<<" P"<<p1.UID<<" TPC "<<p1.TPCID.TPC;
              myprt<<" and P"<<p2.UID<<" TPC "<<p2.TPCID.TPC;
              myprt<<" sep "<<sqrt(maxSep2)<<" maxCth "<<maxCth;
            }
            // see if either of these are in a list
            bool added = false;
            for(auto& pm : stLists) {
              bool p1InList = (std::find(pm.begin(), pm.end(), p1.UID) != pm.end());
              bool p2InList = (std::find(pm.begin(), pm.end(), p2.UID) != pm.end());
              if(p1InList || p2InList) {
                if(p1InList) pm.push_back(p2.UID);
                if(p2InList) pm.push_back(p1.UID);
                added = true;
              }
            } // pm
            if(added) continue;
            // start a new list
            std::vector<int> tmp(2);
            tmp[0] = p1.UID;
            tmp[1] = p2.UID;
            stLists.push_back(tmp);
            break;
          } // p2
        } // p1
      } // sl2
    } // sl1
    if(stLists.empty()) return;
    
    for(auto& stl : stLists) {
      // Find the endpoints of the stitched pfp
      float minZ = 1E6;
      std::pair<unsigned short, unsigned short> minZIndx;
      unsigned short minZEnd = 2;
      for(auto puid : stl) {
        auto slcIndex = GetSliceIndex("P", puid);
        if(slcIndex.first == USHRT_MAX) continue;
        auto& pfp = slices[slcIndex.first].pfps[slcIndex.second];
        for(unsigned short end = 0; end < 2; ++end) {
          auto pos = PosAtEnd(pfp, end);
          if(pos[2] < minZ) { minZ = pos[2]; minZIndx = slcIndex;  minZEnd = end; }
        } // end
      } // puid
      if(minZEnd > 1) continue;
      // preserve the pfp with the min Z position
      auto& pfp = slices[minZIndx.first].pfps[minZIndx.second];
      if(prt) mf::LogVerbatim("TC")<<"SP: P"<<pfp.UID;
      // reverse it if necessary
//      if(minZEnd != 0) ReversePFP(slices[minZIndx.first], pfp);
      // add the Tjs in the other slices to it
      for(auto puid : stl) {
        if(puid == pfp.UID) continue;
        auto sIndx = GetSliceIndex("P", puid);
        if(sIndx.first == USHRT_MAX) continue;
        auto& opfp = slices[sIndx.first].pfps[sIndx.second];
        if(prt) mf::LogVerbatim("TC")<<" +P"<<opfp.UID;
        pfp.TjUIDs.insert(pfp.TjUIDs.end(), opfp.TjUIDs.begin(), opfp.TjUIDs.end());
        if(prt) mf::LogVerbatim();
        // Check for parents and daughters
        if(opfp.ParentUID > 0) {
          auto pSlcIndx = GetSliceIndex("P", opfp.ParentUID);
          if(pSlcIndx.first < slices.size()) {
            auto& parpfp = slices[pSlcIndx.first].pfps[pSlcIndx.second];
            std::replace(parpfp.DtrUIDs.begin(), parpfp.DtrUIDs.begin(), opfp.UID, pfp.UID);
          } // valid pSlcIndx
        } // has a parent
        for(auto dtruid : opfp.DtrUIDs) {
          auto dSlcIndx = GetSliceIndex("P", dtruid);
          if(dSlcIndx.first < slices.size()) {
            auto& dtrpfp = slices[dSlcIndx.first].pfps[dSlcIndx.second];
            dtrpfp.ParentUID = pfp.UID;
          } // valid dSlcIndx
        } // dtruid
        // declare it obsolete
        opfp.ID = 0;
      } // puid
    } // stl

  } // StitchPFPs
  
  /////////////////////////////////////////
  void FindPFParticles(TCSlice& slc)
  {
    // Match Tjs in 3D and create PFParticles
    
    if(tcc.match3DCuts[0] <= 0) return;
    
    bool prt = (tcc.dbgPFP && tcc.dbgSlc);
    if(prt) mf::LogVerbatim("TC")<<" inside FindPFParticles";
    // clear the kEnvInPFP bits on all Tjs. The bit will be set true when a TP is
    // used in a PFParticle
    for(auto& tj : slc.tjs) {
      for(auto& tp : tj.Pts) tp.Environment[kEnvInPFP] = false;
    } // tj
    
    // Match these points in 3D
    std::vector<MatchStruct> matVec;
    // Look for matches in all planes in the TPC
    MatchPlanes(slc, slc.nPlanes, matVec, prt);
    // Make 2-plane matches if we haven't hit the size limit
    if(matVec.size() < tcc.match3DCuts[1] && slc.nPlanes == 3) MatchPlanes(slc, 2, matVec, prt);
    
    if(matVec.empty()) return;
    
    // sort by decreasing number of matched points
    if(matVec.size() > 1) {
      PFPStruct pfp = CreatePFP(slc);
      std::vector<int> dum1;
      std::vector<float> dum2;
      std::vector<SortEntry> sortVec(matVec.size());
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        sortVec[ii].index = ii;
        auto& ms = matVec[ii];
        sortVec[ii].val = ms.Count;
        // de-weight two-plane matches
        if(ms.TjIDs.size() == 2) sortVec[ii].val /= 2;
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
      std::vector<MatchStruct> tmpVec;
      tmpVec.reserve(matVec.size());
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        tmpVec.push_back(matVec[sortVec[ii].index]);
      } // ii
      matVec = tmpVec;
    } // sort matVec
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FPFP: matVec\n";
      unsigned short cnt = 0;
      PFPStruct pfp = CreatePFP(slc);
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        auto& ms = matVec[ii];
        if(ms.Count == 0) continue;
        pfp.TjIDs = ms.TjIDs;
        FindCompleteness(slc, pfp, true, false);
        myprt<<std::setw(4)<<ii<<" Count "<<std::setw(5)<<(int)ms.Count;
        myprt<<" Tj ID-UID:";
        for(auto& tjid : ms.TjIDs) myprt<<" T"<<tjid;
        myprt<<" Comp ";
        for(unsigned short itj = 0; itj < pfp.TjCompleteness.size(); ++itj) {
          myprt<<std::setprecision(2)<<std::setw(6)<<pfp.TjCompleteness[itj];
        }
        myprt<<" Pos ("<<std::setprecision(0)<<std::fixed;
        auto pos = PosAtEnd(pfp, 0);
        myprt<<pos[0]<<", "<<pos[1]<<", "<<pos[2];
        auto dir = DirAtEnd(pfp, 0);
        myprt<<") Dir "<<std::setprecision(2)<<std::setw(6)<<dir[0]<<std::setw(6)<<dir[1]<<std::setw(6)<<dir[2];
        myprt<<" projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, pos, dir, inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
//        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, pfp.TjIDs);
        myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(slc, pfp.TjIDs, false);
        myprt<<"\n";
        ++cnt;
        if(cnt == 500 || ms.Count < 2) {
          myprt<<"...stopped printing after 500 entries or Count < 2";
          break;
        }
      } // ii
    } // prt
    // create a PFParticle for each valid match combination
    for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
      // skip making all PFParticles except for the selected one
      if(debug.MVI < matVec.size() && indx != debug.MVI) continue;
      auto& ms = matVec[indx];
      // ignore dead matches
      if(ms.Count == 0) continue;
      // skip this match if any of the trajectories is already matched or merged or killed
      bool skipit = false;
      // check for incompatible muons, delta rays and InShower Tjs
      bool has13 = false;
      bool has11 = false;
      for(unsigned short itj = 0; itj < ms.TjIDs.size(); ++itj) {
        auto& tj = slc.tjs[ms.TjIDs[itj] - 1];
        if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) skipit = true;
        if(tj.PDGCode == 13) has13 = true;
        if(tj.PDGCode == 11) has11 = true;
      } // tjID
      if(skipit) continue;
      if(has13 && has11) continue;
      int pdgCode = PDGCodeVote(slc, ms.TjIDs, prt);
      PFPStruct pfp = CreatePFP(slc);
      pfp.TjIDs = ms.TjIDs;
      pfp.MVI = indx;
      FindCompleteness(slc, pfp, true, false);
      // skip low TjCompleteness
      if(pfp.TjCompleteness.empty()) {
        std::cout<<"FindPFParticles: pfp.TjCompleteness is empty\n";
        continue;
      }
      skipit = false;
      for(auto tjcomp : pfp.TjCompleteness) if(tjcomp < 0.1) skipit = true;
      if(skipit) {
        if(prt) mf::LogVerbatim()<<" skip it - low completeness";
        continue;
      }
      pfp.PDGCode = pdgCode;
      // set the PDGCode to ensure that delta rays aren't merged with muons. PDGCodeVote
      // returns 0 if the vote is mixed
      if(has13) pfp.PDGCode = 13;
      if(has11) pfp.PDGCode = 11;
      // fill the TP3D points using the 2D trajectory points for Tjs in TjIDs. This function
      // also finds the 3D position of those points using the global fit results from FindCompleteness
      FillTP3Ds(slc, pfp);
      if(pfp.TP3Ds.empty()) continue;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<std::setw(4)<<"MVI "<<indx<<" Count "<<std::setw(5)<<(int)ms.Count;
        myprt<<" P"<<pfp.ID;
        myprt<<" ->";
        for(auto& tjid : pfp.TjIDs) myprt<<" T"<<tjid;
        myprt<<" Comp ";
        for(unsigned short itj = 0; itj < pfp.TjCompleteness.size(); ++itj) {
          myprt<<" "<<std::setprecision(2)<<pfp.TjCompleteness[itj];
        }
        myprt<<" Pos ("<<std::setprecision(0)<<std::fixed;
        auto pos = PosAtEnd(pfp, 0);
        myprt<<pos[0]<<", "<<pos[1]<<", "<<pos[2];
        auto dir = DirAtEnd(pfp, 0);
        myprt<<") Dir ("<<std::setprecision(2)<<dir[0]<<", "<<dir[1]<<", "<<dir[2];
        myprt<<") projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, pos, dir, inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, pfp.TjIDs);
        myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(slc, pfp.TjIDs, false);
        myprt<<" nTP3Ds "<<pfp.TP3Ds.size();
      } // prt
      // Define the PFP, the most important aspect being the 3D trajectory
      if(!Define(slc, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" DefinePFP failed";
        continue;
      }
//      if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("TMP", slc, pfp, -1);
      // FillGaps3D looks for gaps in the TP3Ds vector caused by broken trajectories and
      // inserts new TP3Ds if there are hits in the gaps
      FillGaps3D(slc, pfp, prt);
      // Trim points from the ends until there is one point in two planes
      TrimEndPts(slc, pfp, prt);
      // Look for evidence of broken tjs in one plane
      ChkEndPts(slc, pfp, prt);
      // Look for mis-placed 2D and 3D vertices
//      ReconcileVertices(slc, pfp, prt);
      // set the end flag bits
      geo::TPCID tpcid;
      for(unsigned short end = 0; end < 2; ++end) {
        // first set them all to 0
        pfp.EndFlag[end].reset();
        auto pos = PosAtEnd(pfp, end);
        if(!InsideTPC(pos, tpcid)) pfp.EndFlag[end][kOutFV] = true;
      } // end
      
      // check the ends of the Tjs for Bragg peaks
      unsigned short braggCnt = 0;
      for(auto tid : pfp.TjIDs) {
        auto& tj = slc.tjs[tid - 1];
        if(tj.EndFlag[0][kBragg]) ++braggCnt;
        if(tj.EndFlag[1][kBragg]) ++braggCnt;
      } // tid
      if(braggCnt > 1) {
        std::cout<<"P"<<pfp.ID<<" has Tjs with more than one Bragg peak. Write some code\n";
      }
      // debug?
      if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("FPFP", slc, pfp, -1);
      if(!StorePFP(slc, pfp)) {
        std::cout<<"StorePFP failed P"<<pfp.ID<<"\n";
        continue;
      }
    } // indx
    
    // count the number of TPs and the number that are in PFPs
    unsigned short ntp = 0;
    unsigned short ntpUsed = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        ++ntp;
        if(tp.Environment[kEnvInPFP]) ++ntpUsed;
      } // ipt
    } // tj
    std::cout<<"ntp "<<ntp<<" ntpUsed "<<ntpUsed;
    // do the same for hits
    unsigned short nhtUsed = 0;
    for(unsigned short iht = 0; iht < slc.slHits.size(); ++iht) if(slc.slHits[iht].InPFP > 0) ++nhtUsed;
    std::cout<<" nhits "<<slc.slHits.size()<<" nhtUsed "<<nhtUsed<<"\n";
    
    
    slc.mallTraj.resize(0);

  } // FindPFParticles
  
  /////////////////////////////////////////
  void MakePFPTjs(TCSlice& slc)
  {
    // This function clobbers all of the tjs that are used in TP3Ds in the pfp and replaces
    // them with new tjs that have a consistent set of TPs to prepare for putting them
    // into the event. Note that none of the Tjs are attached to 2D vertices.
    if(!tcc.useAlg[kMakePFPTjs]) return;
    
    bool prt = (tcc.dbgPFP && tcc.dbgSlc);
    
    // kill trajectories
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(auto& tp3d : pfp.TP3Ds) {
        if(!tp3d.IsGood) continue;
        if(tp3d.TjID <= 0) continue;
        unsigned int itj = tp3d.TjID - 1;
        auto& tj = slc.tjs[itj];
        if(tj.AlgMod[kKilled]) continue;
        MakeTrajectoryObsolete(slc, itj);
      } // tp3d
    } // pfp
    
    // Make template trajectories in each plane. These will be re-used by
    // each PFParticle
    std::vector<Trajectory> ptjs(slc.nPlanes);
    // define the basic tj variables
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      ptjs[plane].Pass = 0;
      ptjs[plane].CTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      // This Tj wasn't created by stepping
      ptjs[plane].StepDir = 0;
      // It was created by this function however
      ptjs[plane].AlgMod[kMakePFPTjs] = true;
      // and will be 3D matched
      ptjs[plane].AlgMod[kMat3D] = true;
    } // plane

    // now make the new Tjs
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      // initialize the tjs
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        ptjs[plane].Pts.clear();
        --evt.WorkID;
        if(evt.WorkID == INT_MIN) evt.WorkID = -1;
        ptjs[plane].ID = evt.WorkID;
      } // plane
      pfp.TjIDs.clear();
      // TjUIDs is resized in FinishEvent after this function is called
      pfp.TjCompleteness.clear();
      // iterate through all of the TP3Ds, adding TPs to the TJ in the appropriate plane.
      // The assumption here is that TP order reflects the TP3D order
      for(auto& tp3d : pfp.TP3Ds) {
        if(!tp3d.IsGood) continue;
        if(tp3d.TjID > 0) {
          // a 2D TP
          // Get a reference to the 2D TP
          auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
          unsigned short plane = DecodeCTP(tp.CTP).Plane;
          // append it to Pts
          ptjs[plane].Pts.push_back(tp);
        } else if(tp3d.slHitsIndex < slc.slHits.size()) {
          // a single hit - put it in a TP
          TrajPoint tp;
          unsigned int ahi = slc.slHits[tp3d.slHitsIndex].allHitsIndex;
          auto& hit = (*evt.allHits)[ahi];
          unsigned short plane = hit.WireID().Plane;
          tp.CTP = ptjs[plane].CTP;
          tp.Hits.push_back(tp3d.slHitsIndex);
          tp.UseHit[0] = true;
          DefineHitPos(slc, tp);
          // append it to Pts
          ptjs[plane].Pts.push_back(tp);
        } // a single hit
      } // tp3d
      // finish defining each of the Tjs and store it
      for(auto& tj : ptjs) {
        if(tj.Pts.size() < 2) continue;
        tj.PDGCode = pfp.PDGCode;
        tj.MCSMom = MCSMom(slc, tj);
        if(!StoreTraj(slc, tj)) {
          std::cout<<"MakePFPTjs: StoreTraj failed\n";
          continue;
        } // StoreTraj failed
        // associate it with the pfp
        auto& newTj = slc.tjs[slc.tjs.size() - 1];
        pfp.TjIDs.push_back(newTj.ID);
        pfp.TjCompleteness.push_back(1.);
      } // tj
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"MPFPTjs: P"<<pfp.ID<<" ->";
        for(auto tjid : pfp.TjIDs) myprt<<" T"<<tjid;
      } // prt
    } // pfp
  } // MakePFPTjs
  
  /////////////////////////////////////////
  void MatchPlanes(TCSlice& slc, unsigned short numPlanes, std::vector<MatchStruct>& matVec, bool prt)
  {
    // This function matches trajectory points in slc.mallTraj using the X position. These points should
    // have already been sorted by increasing X by the function that created mallTraj. 
    
    if(slc.mallTraj.empty()) return;
    if(tcc.match3DCuts[0] <= 0) return;
    if(numPlanes < 2) return;
    
    int cstat = DecodeCTP(slc.mallTraj[0].ctp).Cryostat;
    int tpc = DecodeCTP(slc.mallTraj[0].ctp).TPC;
    
    // create a temp vector to check for duplicates
    auto inMatVec = matVec;
    std::vector<MatchStruct> temp;
    
    // the minimum number of points for matching
    unsigned short minPts = 2;
    // override this with the user minimum for 2-plane matches
    if(numPlanes == 2) minPts = tcc.match3DCuts[2];
    
    // max number of match combos left
    unsigned int nAvailable = 0;
    if(matVec.size() < tcc.match3DCuts[1]) nAvailable = tcc.match3DCuts[1] - matVec.size();
    if(nAvailable == 0 || nAvailable > tcc.match3DCuts[1]) return;
    
    // these cuts presume that the resolution in X is better than it is in Y and Z
    float xcut = tcc.match3DCuts[0];
    double yzcut = 1.5 * xcut;
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // length cut
      if(iTjPt.npts < minPts) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      unsigned short iplane = DecodeCTP(itp.CTP).Plane;
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.ctp == iTjPt.ctp) continue;
        // length cut
        if(jTjPt.npts < minPts) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
        TrajPoint3 tp3;
        if(!MakeTp3(slc, itp, jtp, tp3, false)) continue;
        // count weight is one for a two-plane match
        float cntWght = 1;
        if(numPlanes == 3) {
          // numPlanes == 3
          for(unsigned int kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
            auto& kTjPt = slc.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
            if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 5) break;
            auto& ktp = slc.tjs[kTjPt.id - 1].Pts[kTjPt.ipt];
            unsigned short kplane = DecodeCTP(ktp.CTP).Plane;
            TrajPoint3 iktp3;
            if(!MakeTp3(slc, itp, ktp, iktp3, false)) continue;
            if(std::abs(tp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
            if(std::abs(tp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
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
            // Triple match count = 2 
            cntWght = 2;
            // next check the temp vector
            unsigned short indx = 0;
            for(indx = 0; indx < temp.size(); ++indx) {
              auto& ms = temp[indx];
              if(iTjPt.id != ms.TjIDs[iplane]) continue;
              if(jTjPt.id != ms.TjIDs[jplane]) continue;
              if(kTjPt.id != ms.TjIDs[kplane]) continue;
              ms.Count += cntWght;
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
              temp.push_back(ms);
              // give up if there are too many
              if(temp.size() > nAvailable) break;
            } // not found in the list
          } // kpt
          // numPlanes == 3
        } else {
          // 2-plane TPC or 2-plane match in a 3-plane TPC
          if(slc.nPlanes == 3) {
            // See if there is a signal at this point.
            unsigned short kplane = 3 - iplane - jplane;
            float fkwire = tcc.geom->WireCoordinate(tp3.Pos[1], tp3.Pos[2], kplane, tpc, cstat);
            if(fkwire < 0 || fkwire > tcc.maxPos0[kplane]) continue;
            TrajPoint tpk;
            tpk.CTP = EncodeCTP(cstat, tpc, kplane);
            tpk.Pos[0] = fkwire;
            float xp = 0.5 * (iTjPt.xlo + iTjPt.xhi);
            tpk.Pos[1] = tcc.detprop->ConvertXToTicks(xp, kplane, tpc, cstat) * tcc.unitsPerTick;
            // Note that SignalAtTp assumes that a signal exists if the wire is dead
            if(!SignalAtTp(tpk)) continue;
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
    
    if(prt) mf::LogVerbatim("TC")<<"MatchPlanes: Found "<<temp.size()<<" matches in "<<numPlanes<<" planes";
    
  } // MatchPlanes

  /////////////////////////////////////////
  void FindCompleteness(TCSlice& slc, PFPStruct& pfp, bool doFit, bool prt)
  {
    // Calculate the 3D-matching completeness of the set of Tjs in pfp.TjIDs and store in pfp.EffPur.
    // The completeness for each Tj is put in pfp.TjCompleteness. The TP-weighted average completeness
    // is put in pfp.EffPur. This function also fits the matching points to a 3D line and puts the
    // position and direction in pfp.XYZ[0] and pfp.Dir[0].
    
    if(pfp.TjIDs.size() < 2) return;
    if(tcc.match3DCuts[0] <= 0) return;
    if(pfp.SectionFits.empty()) return;
    // This function uses mallTraj but it isn't necessarily a failure if it doesn't exist
    if(slc.mallTraj.size() < 6) return;
    
    pfp.TjCompleteness.resize(pfp.TjIDs.size());
    std::fill(pfp.TjCompleteness.begin(), pfp.TjCompleteness.end(), 0);
    
    bool twoPlanes = (slc.nPlanes == 2);
    bool twoTjs = (pfp.TjIDs.size() == 2);
    double yzcut = 1.5 * tcc.match3DCuts[0];
    
    // initialize the fit sums
    Point3_t point;
    Vector3_t dir;
    if(doFit) Fit3D(0, point, dir, point, dir);
    
    // create a vector of bools for each tj for points that are matched in 3D 
    // cast the IDs into an unsigned short for faster comparing
    std::vector<unsigned short> tjids(pfp.TjIDs.size());
    // This vector is for matches in 3 planes
    std::vector<std::vector<bool>> tjptMat3;
    // This vector is for matches in 2 planes
    std::vector<std::vector<bool>> tjptMat2;
    // and the plane index
    std::vector<unsigned short> tjplane;
    // Initialize the vectors
    for(unsigned short itj = 0; itj < pfp.TjIDs.size(); ++itj) {
      if(pfp.TjIDs[itj] <= 0) {
        std::cout<<"FindCompleteness: Bad tjid "<<pfp.TjIDs[itj]<<"\n";
        return;
      }
      tjids[itj] = pfp.TjIDs[itj];
      auto& tj = slc.tjs[pfp.TjIDs[itj] - 1];
      std::vector<bool> tmp(tj.Pts.size(), false);
      tjptMat2.push_back(tmp);
      if(!twoPlanes) tjptMat3.push_back(tmp);
      tjplane.push_back(DecodeCTP(tj.CTP).Plane);
    } // tjid
    
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      unsigned short indx = 0;
      for(indx = 0; indx < tjids.size(); ++indx) if(iTjPt.id == tjids[indx]) break;
      // require that the Tj ID of this point be in the list
      if(indx == tjids.size()) continue;
      auto& itj = slc.tjs[iTjPt.id - 1];
      //      if(itj.AlgMod[kMat3D]) continue;
      auto& itp = itj.Pts[iTjPt.ipt];
      unsigned short iplane = DecodeCTP(itp.CTP).Plane;
      unsigned short tpc = DecodeCTP(itp.CTP).TPC;
      unsigned short cstat = DecodeCTP(itp.CTP).Cryostat;
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.ctp == iTjPt.ctp) continue;
        unsigned short jndx = 0;
        for(jndx = 0; jndx < tjids.size(); ++jndx) if(jTjPt.id == tjids[jndx]) break;
        // require that the Tj ID of this point be in the list
        if(jndx == tjids.size()) continue;
        // check for x range overlap. We know that jTjPt.xlo is > iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtj = slc.tjs[jTjPt.id - 1];
        //        if(jtj.AlgMod[kMat3D]) continue;
        auto& jtp = jtj.Pts[jTjPt.ipt];
        TrajPoint3 ijtp3;
        if(!MakeTp3(slc, itp, jtp, ijtp3, true)) continue;
        ijtp3.Tj2Pts.resize(2);
        ijtp3.Tj2Pts[0] = iTjPt;
        ijtp3.Tj2Pts[1] = jTjPt;
        // Set the 2-plane match bits
        tjptMat2[indx][iTjPt.ipt] = true;
        tjptMat2[jndx][jTjPt.ipt] = true;
        if(twoPlanes) {
          if(doFit) Fit3D(1, ijtp3.Pos, ijtp3.Dir, point, dir);
          continue;
        }
        // count it as a triple if this point is in a dead region
        unsigned short jplane = DecodeCTP(jtp.CTP).Plane;
        unsigned short kplane = 3 - iplane - jplane;
        float fwire = tcc.geom->WireCoordinate(ijtp3.Pos[1], ijtp3.Pos[2], kplane, tpc, cstat);
        if(fwire < -0.4) continue;
        unsigned int kwire = std::nearbyint(fwire);
        if(kwire < slc.wireHitRange[kplane].size() && !evt.goodWire[kplane][kwire]) {
          // accumulate the fit sums
          if(doFit) Fit3D(1, ijtp3.Pos, ijtp3.Dir, point, dir);
          continue;
        } // dead wire in kplane
        for(unsigned int kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
          auto& kTjPt = slc.mallTraj[kpt];
          // ensure that the planes are different
          if(kTjPt.ctp == iTjPt.ctp || kTjPt.ctp == jTjPt.ctp) continue;
          // Look for this tj point in tjids
          unsigned short kndx = 0;
          for(kndx = 0; kndx < tjids.size(); ++kndx) if(kTjPt.id == tjids[kndx]) break;
          if(kTjPt.xlo > iTjPt.xhi) continue;
          // break out if the x range difference becomes large
          if(kTjPt.xlo > iTjPt.xhi + 5) break;
          auto& ktj = slc.tjs[kTjPt.id - 1];
          auto& ktp = ktj.Pts[kTjPt.ipt];
          TrajPoint3 iktp3;
          if(!MakeTp3(slc, itp, ktp, iktp3, true)) continue;
          if(std::abs(ijtp3.Pos[1] - iktp3.Pos[1]) > yzcut) continue;
          if(std::abs(ijtp3.Pos[2] - iktp3.Pos[2]) > yzcut) continue;
          // make a copy of ijtp3 -> ijktp3
          auto ijktp3 = ijtp3;
          // add the Tj2Pt to it
          ijktp3.Tj2Pts.push_back(kTjPt);
          // accumulate the fit sums
          if(doFit) Fit3D(1, iktp3.Pos, iktp3.Dir, point, dir);
          // Set the 3-plane match bits
          if(kndx == tjids.size()) continue;
          tjptMat3[indx][iTjPt.ipt] = true;
          tjptMat3[jndx][jTjPt.ipt] = true;
          tjptMat3[kndx][kTjPt.ipt] = true;
        } // kpt
      } // jpt
    } // ipt
    // do the fit and put the results into the pfp
    auto& sf = pfp.SectionFits[0];
    Fit3D(2, point, dir, sf.Pos, sf.Dir);
    if(prt && doFit) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FC: P"<<pfp.ID<<" fit pos "<<std::fixed<<std::setprecision(1)<<sf.Pos[0]<<" "<<sf.Pos[1]<<" "<<sf.Pos[2];
      myprt<<" fit dir "<<std::setprecision(2)<<sf.Dir[0]<<" "<<sf.Dir[1]<<" "<<sf.Dir[2];
      myprt<<" Note: fit pos is the average position of all Tp3s - not the start or end.";
    }
    // now count the number of tj points were matched
    // total number of points with charge in all Tjs
    float tnpwc = 0;
    // total number that are matched in 3D in 3 planes
    float tcnt3 = 0;
    // total number that are matched in 3D in 2 planes
    float tcnt2 = 0;
    for(unsigned short itj = 0; itj < tjids.size(); ++itj) {
      auto& tj = slc.tjs[tjids[itj] - 1];
      // counts for each tj
      float npwc = 0;
      float cnt2 = 0;
      float cnt3 = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg <= 0) continue;
        ++npwc;
        if(tjptMat2[itj][ipt]) ++cnt2;
        if(!twoPlanes && tjptMat3[itj][ipt]) ++cnt3;
      } // ipt
      if(twoTjs) {
        pfp.TjCompleteness[itj] = cnt2 / npwc;
      } else {
        pfp.TjCompleteness[itj] = cnt3 / npwc;
      }
      tnpwc += npwc;
      tcnt3 += cnt3;
      tcnt2 += cnt2;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"FC: P"<<pfp.ID<<" T"<<tj.ID<<" npwc "<<npwc<<" cnt2 "<<cnt2<<" cnt3 "<<cnt3<<" PDGCode "<<tj.PDGCode;
        myprt<<" MCSMom "<<tj.MCSMom<<" InShower? "<<tj.SSID;
        myprt<<" TjCompleteness "<<std::setprecision(2)<<pfp.TjCompleteness[itj];
      } // prt
    } // itj
    if(twoTjs) {
      pfp.EffPur = tcnt2 / tnpwc;
    } else {
      pfp.EffPur = tcnt3 / tnpwc;
    }
    
  } // FindCompleteness

  /////////////////////////////////////////
  bool Define(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function is called after the 3D matched TjIDs have been specified and optionally
    // a start or end vertex ID. It defines the PFParticle but doesn't store it
    
    if(pfp.PDGCode == 1111) return false;
    if(pfp.TjIDs.size() < 2) return false;
    
    // require at least one tj in at least two planes
    std::vector<unsigned short> nInPln(slc.nPlanes);
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      ++nInPln[DecodeCTP(tj.CTP).Plane];
      if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) {
        std::cout<<"DPFP: P"<<pfp.ID<<" uses T"<<tj.ID<<" but kMat3D is set true\n";
        return false;
      }
    } // tjid
    unsigned short npl = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(nInPln[plane] > 0) ++npl;
    if(npl < 2) return false;
    
    // and sort by the distance along the trajectory (just a 3D line at this point) in
    // section 0 (the only one)
    if(!SortSection(pfp, 0)) return false;
    
    // do a fit of all points to find chiDOF and update the TP3Ds
    if(!FitSection(slc, pfp, 0)) {
      std::cout<<"DPFP: First fit failed\n";
      return false;
    }
    // sort again (TODO Is this necessary?)
    if(!SortSection(pfp, 0)) return false;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"DPFP: FitSection: P"<<pfp.ID<<" first fit";
      myprt<<std::fixed<<std::setprecision(1);
      auto pos = PosAtEnd(pfp, 0);
      myprt<<" Pos ("<<pos[0]<<" "<<pos[1]<<" "<<pos[2];
      myprt<<std::fixed<<std::setprecision(2);
      auto dir = DirAtEnd(pfp, 0);
      myprt<<") Dir ("<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<")";
    }
    
    // See if the TP3Ds need to be fit in sections. First test for a bad chiDOF
    auto& sf = pfp.SectionFits[0];
    if(sf.ChiDOF > 2 && !ReSection(slc, pfp, prt)) {
      std::cout<<"ReSection failed\n";
      return false;
      if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("RS", slc, pfp, -1);
    } // ReSection

    // Check the ordering of the TPs and possibly declare the out-of-order
    // TPs invalid. Update is called within FixOrder
    FixOrder(slc, pfp, prt);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"DPFP: P"<<pfp.ID;
      myprt<<" ->";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
      myprt<<" max Tj len "<<MaxTjLen(slc, pfp.TjIDs);
      myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
      myprt<<" TP3Ds "<<pfp.TP3Ds.size();
      myprt<<" EndFlags ";
      for(unsigned short end = 0; end < 2; ++end) myprt<<PrintEndFlag(pfp, end);
    } // prt
    
    return true;
  } // Define

  /////////////////////////////////////////
  void FixOrder(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Check the wire order of the 2D TPs and fix them
    
    std::vector<short> plnOrder(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      auto pos = PosAtEnd(pfp, 0);
      auto dir = DirAtEnd(pfp, 0);
      auto tp = MakeBareTP(slc, pos, dir, inCTP);
      plnOrder[plane] = 1;
      if(tp.Dir[0] < 0) plnOrder[plane] = -1;
      float prevPos0 = 0;
      // count of wire minus/plus
      unsigned short nmi = 0;
      unsigned short npl = 0;
      bool first = true;
      unsigned short cnt = 0;
      for(auto& tp3d : pfp.TP3Ds) {
        if(!tp3d.IsGood) continue;
        // ignore single hits
        if(tp3d.TjID <= 0) continue;
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if(tp.CTP != inCTP) continue;
        ++cnt;
        if(first) {
          first = false;
          prevPos0 = tp.Pos[0];
        } else {
          float dw = tp.Pos[0] - prevPos0;
          if(dw > 0) { ++npl; } else { ++nmi; }
          prevPos0 = tp.Pos[0];
        } // !first
      } // tp3d
      // Nothing to do if the order is all positive or all negative
      if(npl == 0 || nmi == 0) continue;
      short ordr = 1;
      if(nmi > npl) ordr = -1;
      std::cout<<"FixOrder found a problem: plane "<<plane<<" nmi "<<nmi<<" npl "<<npl<<" ordr "<<ordr;
      std::cout<<" plnOrder "<<plnOrder[plane]<<" dir "<<tp.Dir[0];
      std::cout<<" cnt "<<cnt<<"\n";
      // create a list of indices into TP3D of points in this plane to simplify checking
      // adjacent TPs
      std::vector<unsigned short> indxInPln;
      for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(!tp3d.IsGood) continue;
        if(tp3d.TjID <= 0) continue;
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if(tp.CTP != inCTP) continue;
        indxInPln.push_back(ipt);
      } // tp3d
      // check adjacent pairs of TPs to find the order
      for(unsigned short ii = 1; ii < indxInPln.size() - 1; ++ii) {
        auto& thisTP3 = pfp.TP3Ds[indxInPln[ii]];
        auto& thisTP2 = slc.tjs[thisTP3.TjID - 1].Pts[thisTP3.TPIndex];
        auto& prevTP3 = pfp.TP3Ds[indxInPln[ii-1]];
        auto& prevTP2 = slc.tjs[prevTP3.TjID - 1].Pts[prevTP3.TPIndex];
        auto& nextTP3 = pfp.TP3Ds[indxInPln[ii+1]];
        auto& nextTP2 = slc.tjs[nextTP3.TjID - 1].Pts[nextTP3.TPIndex];
        short prevdw = std::nearbyint(thisTP2.Pos[0] - prevTP2.Pos[0]);
        short nextdw = std::nearbyint(nextTP2.Pos[0] - thisTP2.Pos[0]);
        if(prevdw == ordr && nextdw == ordr) continue;
        if(prevdw != ordr && nextdw != ordr) {
          // mark both as not good
          prevTP3.IsGood = false;
          thisTP3.IsGood = false;
          auto& sf = pfp.SectionFits[thisTP3.SFIndex];
          sf.NeedsUpdate = true;
          pfp.NeedsUpdate = true;
        } else if(prevdw != ordr && nextdw == ordr) {
          // mark this not-good
          thisTP3.IsGood = false;
          auto& sf = pfp.SectionFits[thisTP3.SFIndex];
          sf.NeedsUpdate = true;
          pfp.NeedsUpdate = true;
        } else {
          std::cout<<"Weird indxInPln "<<indxInPln[ii]<<" "<<PrintPos(slc, thisTP2)<<" prevdw "<<prevdw<<" nextdw "<<nextdw<<" ordr "<<ordr<<"\n";
        }
      } // ii
    } // plane
    if(pfp.NeedsUpdate) Update(slc, pfp, prt);
  } // FixOrder

  /////////////////////////////////////////
  bool Update(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function only updates SectionFits that need to be re-sorted or re-fit. It returns
    // false if there was a serious error
    if(pfp.TP3Ds.empty()) return false;
    
    // Re-build the TP3Ds vector if the first or last points are not Good. This ensures
    // that references to the start and end position of the PFP are valid. Only the end
    // TP3Ds points are checked. Not good-Points interior to the vector may later be deemed to
    // be good so don't eliminate them.
    bool reBuild = !pfp.TP3Ds[0].IsGood || !pfp.TP3Ds[pfp.TP3Ds.size() - 1].IsGood;
    if(reBuild) {
      if(prt) mf::LogVerbatim("TC")<<"Update: P"<<pfp.ID<<" needs a re-build. First/last points are !IsGood";
      // find the first good point
      unsigned short firstGood = 0;
      for(firstGood = 0; firstGood < pfp.TP3Ds.size(); ++firstGood) if(pfp.TP3Ds[firstGood].IsGood) break;
      if(firstGood == pfp.TP3Ds.size()) return false;
      // and the last good point
      unsigned short lastGood = pfp.TP3Ds.size();
      for(unsigned short ii = 1; ii < pfp.TP3Ds.size(); ++ii) {
        lastGood = pfp.TP3Ds.size() - ii;
        if(pfp.TP3Ds[lastGood].IsGood) break;
      } // ii
      ++lastGood;
      std::vector<TP3D> temp(pfp.TP3Ds.begin() + firstGood, pfp.TP3Ds.begin() + lastGood);
      pfp.TP3Ds = temp;
    } // reBuild

    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      if(sf.NPts == 0) continue;
      if(sf.ChiDOF <= 0) {
        std::cout<<"Update: P"<<pfp.ID<<" NPts is > 0 but ChiDOF <= 0. Does this make sense?\n";
      }
      if(!sf.NeedsUpdate) continue;
      if(!FitSection(slc, pfp, sfi)) return false;
      if(!SortSection(pfp, sfi)) return false;
      sf.NeedsUpdate = false;
    } // sf
    pfp.NeedsUpdate = false;
    return true;
  } // Update
  
  /////////////////////////////////////////
  bool ReSection(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Re-fit the TP3Ds in sections and add/remove sections to keep ChiDOF of each section close to 1
    if(pfp.TP3Ds.size() < 5 || pfp.SectionFits.empty()) return false;
    
    prt = false;
    
    // clobber the old sections if more than one exists
    if(pfp.SectionFits.size() > 1) {
      // make one section
      pfp.SectionFits.resize(1);
      // put all of the points in it and fit
      for(auto& tp3d : pfp.TP3Ds) {
        tp3d.SFIndex = 0;
        tp3d.IsGood = true;
      }
      auto& sf = pfp.SectionFits[0];
      if(!FitSection(slc, pfp, 0)) {
        std::cout<<"ReSection: First fit failed\n";
        return false;
      } // fit failed
      if(sf.ChiDOF < 2) return true;
    } // > 1 SectionFit
    // sort by distance from the start
    if(!SortSection(pfp, 0)) return false;
    // temp vector for debugging
    std::vector<bool> used(pfp.TP3Ds.size(), false);
    // Next see if it is long enough to be reconstructed in more than one section
    // Scale by the MCSMom to have 20 points in a plane in each section for a PFParticle
    // with MCSMom = 1000;
    unsigned short min2DPts = MCSMom(slc, pfp.TjIDs) / 50;
    if(min2DPts < (unsigned short)tcc.match3DCuts[3]) min2DPts = tcc.match3DCuts[3];
    if(prt) mf::LogVerbatim("TC")<<"ReSection: min2DPts "<<min2DPts;
    unsigned short fromPt = 0;
    unsigned short toPt = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1);
    if(toPt > pfp.TP3Ds.size()) return true;
    // put all points into an "undefined" section
    for(auto& tp3d : pfp.TP3Ds) tp3d.SFIndex = USHRT_MAX;
    // iterate to create sections
    for(unsigned short snit = 0; snit < 40; ++snit) {
      unsigned short sfi = pfp.SectionFits.size() - 1;
      auto& sf = pfp.SectionFits[sfi];
      if(prt) mf::LogVerbatim("TC")<<" snit "<<snit<<" sfi "<<sfi<<" fromPt "<<fromPt<<" toPt "<<toPt;
      // put these points into this section between fromPt and toPt
      for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) pfp.TP3Ds[ipt].SFIndex = sfi;
      // debugging
      for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) used[ipt] = true;
      // Iterate to fit the section, then add or remove points to get a reasonable ChiDOF
      float prevChi = 0;
      for(unsigned short nit = 0; nit < 20; ++nit) {
        if(!FitSection(slc, pfp, sfi)) {
          std::cout<<"section fit "<<sfi<<" failed\n";
          sf.ChiDOF = -1;
          break;
        } // fit failed
        if(prt) mf::LogVerbatim("TC")<<"  nit "<<nit<<" sf.ChiDOF "<<sf.ChiDOF<<" sf.NPts "<<sf.NPts<<" fromPt "<<fromPt<<" toPt "<<toPt<<" prevChi "<<prevChi;
        if(sf.ChiDOF == prevChi) break;
        if(sf.ChiDOF > 2) {
          // largish chisq. Try to remove a bad point and break if one wasn't found
          if(!KillBadPoint(slc, pfp, fromPt, toPt, prt)) break;
        } else if(sf.ChiDOF < 0.5) {
          // add about 50% more points
          unsigned short nadd = 0.5 * (toPt - fromPt);
          toPt += nadd;
          if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
          for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) pfp.TP3Ds[ipt].SFIndex = sfi;
          // debugging
          for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) used[ipt] = true;
          // sf.ChiDOF < 0.5
        } else {
          // ChiDOF is reasonable so break
          break;
        }
        prevChi = sf.ChiDOF;
      } // nit
      if(toPt > pfp.TP3Ds.size() || sf.NPts == 0) break;
      // sort by distance from the start
      if(!SortSection(pfp, sfi)) {
        std::cout<<"SortSection failed in section "<<sfi<<"\n";
        continue;
      }
      // set the maxAlong variable for this section
//      sf.maxAlong = pfp.TP3Ds[toPt - 1].along;
      fromPt = toPt;
      toPt = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1);
      if(toPt > pfp.TP3Ds.size()) break;
      // add a new section
      pfp.SectionFits.resize(pfp.SectionFits.size() + 1);
    } // snit
    
    // It is likely that the last points weren't included in the last section
    unsigned short lastUsed = 0;
    for(lastUsed = 0; lastUsed < pfp.TP3Ds.size(); ++lastUsed) if(pfp.TP3Ds[lastUsed].IsGood && !used[lastUsed]) break;
    if(lastUsed < pfp.TP3Ds.size()) {
      // Keep it simple: just add all of the points to the last section and
      // don't worry about ChiDOF
      unsigned short sfi = pfp.SectionFits.size() - 1;
      for(unsigned short ipt = lastUsed; ipt < pfp.TP3Ds.size(); ++ipt) pfp.TP3Ds[ipt].SFIndex = sfi;
      // just fit, sort and find end points without checking
      FitSection(slc, pfp, sfi);
      SortSection(pfp, sfi);
    } // 
    
    return true;
  } // resection
  
  /////////////////////////////////////////
  bool KillBadPoint(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, bool prt)
  {
    // Sets IsGood false if pull < match3DCuts[4] for a single point in the range and
    // returns true if one is found
    if(fromPt > pfp.TP3Ds.size() - 1) return false;
    if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
    float maxPull = tcc.match3DCuts[4];
    unsigned short badPt = USHRT_MAX;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      // don't clobber a point if it is on a TP that is overlapping another Tj. This will
      // happen for points close to a vertex and when trajectories cross
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvOverlap]) continue;
      float pull = std::abs(tp3d.Pos[0] - tp3d.TPX) / tp3d.TPXErr2;
      if(pull < maxPull) continue;
      maxPull = pull;
      badPt = ipt;
    } // ipt
    if(badPt == USHRT_MAX) return false;
    // clobber that point
    pfp.TP3Ds[badPt].IsGood = false;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<" KBP: P"<<pfp.ID<<" clobber "<<pfp.TP3Ds[badPt].TjID<<"_"<<pfp.TP3Ds[badPt].TPIndex;
    } // prt
    return true;
  } // KillBadPoint
  
  /////////////////////////////////////////
  unsigned short Find3DRecoRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short min2DPts, short dir)
  {
    // Scans the TP3Ds vector starting at fromPt until it finds min2DPts in two planes. It returns
    // with the index of that point (+1) in the TP3Ds vector. The dir variable defines the scan direction in
    // the TP3Ds vector
    if(fromPt > pfp.TP3Ds.size() - 1) return USHRT_MAX;
    if(pfp.TP3Ds.size() < 2 * min2DPts) return USHRT_MAX;
    if(dir == 0) return USHRT_MAX;
    
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      unsigned short ipt = fromPt + ii;
      if(dir < 0) ipt = fromPt - ii;
      if(ipt >= pfp.TP3Ds.size()) break;
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      unsigned short plane = 0;
      if(tp3d.TjID > 0) {
        plane = DecodeCTP(slc.tjs[tp3d.TjID - 1].CTP).Plane;
      } else {
        auto& hit = (*evt.allHits)[slc.slHits[tp3d.slHitsIndex].allHitsIndex];
        plane = hit.WireID().Plane;
      }
      ++cntInPln[plane];
      unsigned short enufInPlane = 0;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] >= min2DPts) ++enufInPlane;
      if(enufInPlane > 1) return ipt + 1;
      if(dir < 0 && ipt == 0) break;
    } // ipt
    return USHRT_MAX;
  } // Find3DRecoRange

  /////////////////////////////////////////
  bool FitSection(TCSlice& slc, PFPStruct& pfp, unsigned short sfIndex)
  {
    // Fits the TP3D points in the selected section to a 3D line with the origin at the center of
    // the section
    if(pfp.TP3Ds.size() < 4) return false;
    if(sfIndex >= pfp.SectionFits.size()) return false;
    auto& sf = pfp.SectionFits[sfIndex];
    sf.ChiDOF = 999;

    // put the offset, cosine-like and sine-like components in a vector
    std::vector<std::array<double, 3>> ocs(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      auto planeID = geo::PlaneID(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
      // plane offset
      ocs[plane][0] = tcc.geom->WireCoordinate(0, 0, planeID);
      // get the "cosine-like" component
      ocs[plane][1] = tcc.geom->WireCoordinate(1, 0, planeID) - ocs[plane][0];
      // the "sine-like" component
      ocs[plane][2] = tcc.geom->WireCoordinate(0, 1, planeID) - ocs[plane][0];
    } // plane

    // count the number of points in this section
    sf.NPts = 0;
    // and define the X position for the fit origin
    double x0 = 0.;
    // vector of points in this section
    std::vector<unsigned short> tpIndex;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      if(tp3d.SFIndex != sfIndex) continue;
      tpIndex.push_back(ipt);
      x0 += tp3d.TPX;
    } //indx
    sf.NPts = tpIndex.size();
    // Need at least 4 + 1 points for a chisq calculation
    if(sf.NPts < 5) return false;
    x0 /= (double)sf.NPts;
    
    const unsigned int nvars = 4;
    TMatrixD A(sf.NPts, nvars);
    // vector holding the Wire number
    TVectorD w(sf.NPts);
    
    // count the number of TPs in each plane
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for(unsigned short indx = 0; indx < tpIndex.size(); ++indx) {
      auto& tp3d = pfp.TP3Ds[tpIndex[indx]];
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      ++cntInPln[plane];
      double x = tp3d.TPX - x0;
      double wght = 1 / tp3d.TPXErr2;
      A[indx][0] = wght * ocs[plane][1];
      A[indx][1] = wght * ocs[plane][2];
      A[indx][2] = wght * ocs[plane][1] * x;
      A[indx][3] = wght * ocs[plane][2] * x;
      w[indx] = wght * (tp3d.Wire - ocs[plane][0]);
    } // ipt
    // ensure there are at least two points in at least two planes
    unsigned short enufInPlane = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] > 1) ++enufInPlane;
    if(enufInPlane < 2) return false;
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    // TODO: The direction is reversed for some reason
    norm *= -1;
    sf.Dir[0] = 1 / norm;
    sf.Dir[1] = tVec[2] / norm;
    sf.Dir[2] = tVec[3] / norm;
    sf.Pos[0] = x0;
    sf.Pos[1] = tVec[0];
    sf.Pos[2] = tVec[1];
    // calculate ChiDOF
    sf.ChiDOF = 0;
    // project this 3D vector into a TP in every plane
    std::vector<TrajPoint> plnTP(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      plnTP[plane] = MakeBareTP(slc, sf.Pos, sf.Dir, inCTP);
    } // plane
    for(unsigned short indx = 0; indx < tpIndex.size(); ++indx) {
      auto& tp3d = pfp.TP3Ds[tpIndex[indx]];
      tp3d.Dir = sf.Dir;
      // interpoloate to find the position
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      double dw = tp3d.Wire - plnTP[plane].Pos[0];
      // dt/dW was stored in DeltaRMS
      double t = dw * plnTP[plane].DeltaRMS;
      for(unsigned short xyz = 0; xyz < 3; ++xyz) tp3d.Pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
      tp3d.along = t;
      // Note that the tp3d position is directly above the wire position and not the
      // distance of closest approach. The Delta variable is the difference in the
      // drift direction in cm
      double delta = tp3d.Pos[0] - tp3d.TPX;
      sf.ChiDOF += delta * delta / tp3d.TPXErr2;
    } // indx

    sf.ChiDOF /= (float)(sf.NPts - 4);
    
    return true;
    
  } // FitSection

  /////////////////////////////////////////
  void ReconcileVertices(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Checks for mis-placed 2D and 3D vertices and either attaches them
    // to a vertex or deletes(?) the vertex while attempting to preserve or
    // correct the P -> T -> 2V -> 3V assn. After this is done, the function
    // TCVertex/AttachToAnyVertex is called.
    // This function returns true if something was done to the pfp that requires
    // a re-definition of the pfp, e.g. adding or removing TP3Ds. Note that this
    // never occurs as the function is currently written
    
    if(pfp.NeedsUpdate) {
      std::cout<<"ReconcileVertices: P"<<pfp.ID<<" NeedsUpdate is true at entry. Fixing it...\n";
      if(!Update(slc, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" Define failed inside ReconcileVertices"; 
        return;
      }
    }
    pfp.NeedsUpdate = false;
    if(pfp.TP3Ds.empty()) return;
    
    // first make a list of all Tjs 
    std::vector<int> tjList;
    for(auto& tp3d : pfp.TP3Ds) {
      if(!tp3d.IsGood) continue;
      // ignore single hits
      if(tp3d.TjID <= 0) continue;
      if(std::find(tjList.begin(), tjList.end(), tp3d.TjID) == tjList.end()) tjList.push_back(tp3d.TjID);
    } // tp3d
    // look for 3D vertices associated with these Tjs and list of
    // orphan 2D vertices - those that are not matched to 3D vertices
    std::vector<int> vx2List, vx3List;
    for(auto tid : tjList) {
      auto& tj = slc.tjs[tid - 1];
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] <= 0) continue;
        auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
        if(vx2.Vx3ID > 0) {
          if(std::find(vx3List.begin(), vx3List.end(), vx2.Vx3ID) == vx3List.end()) vx3List.push_back(vx2.Vx3ID);
          // 3D vertex exists
        } else {
          // no 3D vertex 
          if(std::find(vx2List.begin(), vx2List.end(), tj.VtxID[end]) == vx2List.end()) vx2List.push_back(tj.VtxID[end]);
        } // no 3D vertex
      } // end
    } // tid
    // no vertex reconciliation is necessary
    if(vx2List.empty() && vx3List.empty()) return;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"RV: P"<<pfp.ID<<" ->";
      for(auto tid : tjList) myprt<<" T"<<tid;
      myprt<<" ->";
      for(auto vid : vx3List) myprt<<" 3V"<<vid;
      if(!vx2List.empty()) {
        myprt<<" orphan";
        for(auto vid : vx2List) myprt<<" 2V"<<vid;
      }
    } // prt
    // Just kill the orphan 2D vertices regardless of their score. 
    // This is an indicator that the vertex was created between two tjs 
    // that maybe should have been reconstructed as one or alternatively 
    // as two Tjs. This decision presumes the existence of a 3D kink 
    // algorithm that doesn't yet exist...
    for(auto vid : vx2List) {
      auto& vx2 = slc.vtxs[vid - 1];
      MakeVertexObsolete("RV", slc, vx2, true);
    } // vx2List
    // ignore the T -> 2V -> 3V assns (if any exist) and try to directly
    // attach to 3D vertices at both ends
    AttachToAnyVertex(slc, pfp, prt);
    // check for differences
    for(unsigned short end = 0; end < 2; ++end) {
      // see if a vertex got attached
      if(pfp.Vx3ID[end] <= 0) continue;
      // see if this is a vertex in the list using the T -> 2V -> 3V assns
      if(std::find(vx3List.begin(), vx3List.end(), pfp.Vx3ID[end]) != vx3List.end()) continue;
      std::cout<<"RV: P"<<pfp.ID<<" was attached to 3V"<<pfp.Vx3ID[end]<<" but a P -> T -> 2V -> 3V assn exists. Write some code to clobber this assn or deal with it somehow.\n";
    } // end
    
    return;
  } // ReconcileVertices
  
  /////////////////////////////////////////
  void ChkEndPts(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function looks for 
    if(pfp.ID <= 0) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.SectionFits.empty()) return;
    
    for(unsigned short end = 0; end < 2; ++end) {
      // find the first good point in each plane from the end
      std::vector<unsigned short> firstPt(slc.nPlanes, USHRT_MAX);
      unsigned short minPlnPt = USHRT_MAX;
      unsigned short maxPlnPt = 0;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        firstPt[plane] = FirstPointInPlane(pfp, plane, end);
        if(firstPt[plane] < minPlnPt) minPlnPt = firstPt[plane];
        if(firstPt[plane] > maxPlnPt) maxPlnPt = firstPt[plane];
      } // plane
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"CEP: P"<<pfp.ID<<" end "<<end<<" plane_firstPt";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) myprt<<" "<<plane<<"_"<<firstPt[plane];
      } // prt
      if(maxPlnPt > minPlnPt + 8) {
        unsigned short fromPt = 0;
        unsigned short toPt = pfp.TP3Ds.size();
        unsigned short inPlane = 0;
        if(end == 0) {
          toPt = maxPlnPt;
          for(inPlane = 0; inPlane < slc.nPlanes; ++inPlane) if(firstPt[inPlane] == maxPlnPt) break;
        } else {
          fromPt = minPlnPt;
          for(inPlane = 0; inPlane < slc.nPlanes; ++inPlane) if(firstPt[inPlane] == minPlnPt) break;
        }
        std::cout<<"CEP: P"<<pfp.ID<<" check points "<<fromPt<<" to "<<toPt<<" in plane "<<inPlane<<"\n";
        AddPointsInRange(slc, pfp, fromPt, toPt, inPlane, tcc.match3DCuts[4], prt);
      } // maxPlnPt > minPlnPt + 8
    } // end
    if(!Update(slc, pfp, prt)) {
      std::cout<<"Update failed in ChkEndPts. Debug this\n";
    }
  } // ChkEndPts
  
  /////////////////////////////////////////
  void TrimEndPts(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Check the end TP3Ds points for consistency and trim if necessary.
    // Ideally there will be a TP3D in each plane in any local region of ~few cm
    // along the 3D trajectory, unless there are dead wires or there are overlapping
    // trajectories in this view, mostly likely due to the presence of a nearby
    // vertex. This function only checks the ends.
    if(pfp.ID <= 0) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.SectionFits.empty()) return;
    
    for(unsigned short end = 0; end < 2; ++end) {
      // find the first good point in each plane from the end
      std::vector<unsigned short> firstPt(slc.nPlanes, USHRT_MAX);
      unsigned short minPlnPt = USHRT_MAX;
      unsigned short maxPlnPt = 0;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        firstPt[plane] = FirstPointInPlane(pfp, plane, end);
        if(firstPt[plane] < minPlnPt) minPlnPt = firstPt[plane];
        if(firstPt[plane] > maxPlnPt) maxPlnPt = firstPt[plane];
      } // plane
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<"TEP: P"<<pfp.ID<<" end "<<end<<" plane_firstPt";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) myprt<<" "<<plane<<"_"<<firstPt[plane];
      } // prt
      // expect that the min and max points are similar if the ends are well-reconstructed
      if(maxPlnPt > minPlnPt + 8) {
        unsigned short fromPt = 0;
        unsigned short toPt = pfp.TP3Ds.size();
        unsigned short trimPt = 0;
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          if(firstPt[plane] != minPlnPt && firstPt[plane] != maxPlnPt) trimPt = firstPt[plane];
        } // plane
        if(end == 0) {
          toPt = trimPt;
        } else {
          fromPt = trimPt;
        }
        for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
          auto& tp3d = pfp.TP3Ds[ipt];
          // already declared not good?
          if(!tp3d.IsGood) continue;
          if(tp3d.TjID > 0) {
            auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
            // don't trim the point if it is overlapped with a different tj (e.g. close to a vertex)
            if(tp.Environment[kEnvOverlap]) continue;
          } // tp3d.TjID > 0
          // check the pull
//          double pull = std::abs(tp3d.Pos[0] - tp3d.TPX) / sqrt(tp3d.TPXErr2);
//          std::cout<<" pull "<<pull<<"\n";
          tp3d.IsGood = false;
          auto& sfi = tp3d.SFIndex;
          pfp.SectionFits[sfi].NeedsUpdate = true;
          pfp.NeedsUpdate = true;
        } // ipt
      } // maxPlnPt > minPlnPt + 5
    } // end
    if(!Update(slc, pfp, prt)) {
      std::cout<<"Update failed in TrimEndPts. Debug this\n";
    }
  } // TrimEndPts
  
  /////////////////////////////////////////
  unsigned short FirstPointInPlane(PFPStruct& pfp, unsigned short plane, unsigned short end)
  {
    // returns the first (last) point in the TP3Ds that is in the plane when end = 0 (1)
    if(pfp.ID <= 0) return USHRT_MAX;
    if(pfp.TP3Ds.empty()) return USHRT_MAX;
    // assume end = 0
    short dir = 1;
    short endPt = 0;
    if(end == 1) {
      dir = -1;
      endPt = pfp.TP3Ds.size() - 1;
    }
    for(short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      short sipt = endPt + ii * dir;
      if(sipt < 0 || sipt >= pfp.TP3Ds.size()) break;
      unsigned short ipt = sipt;
      if(DecodeCTP(pfp.TP3Ds[ipt].CTP).Plane == plane) return ipt;
    } // ii
    return USHRT_MAX;
  } // FirstPointInPlane
  
  /////////////////////////////////////////
  void FillGaps3D(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Look for gaps in each plane in the TP3Ds vector. Note that this only looks for
    // interior gaps. 
    if(pfp.ID <= 0) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.SectionFits.empty()) return;
    if(!tcc.useAlg[kFillGaps3D]) return;
    
    // Ignore gaps that are less than (10 cm)^2 in length
    constexpr float minGapLength2 = 100;
    
    // project the first section fit into a TP in each plane
    std::vector<TrajPoint> plnTP(slc.nPlanes);
    unsigned short sfi = 0;
    auto& sf = pfp.SectionFits[sfi];
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      plnTP[plane] = MakeBareTP(slc, sf.Pos, sf.Dir, inCTP);
    } // plane
    
    // index of the previous point in each plane
    std::vector<unsigned short> prevPt(slc.nPlanes, USHRT_MAX);
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      // Don't test for a bad point since we are looking for missing points
      // update the plnTPs?
      if(tp3d.SFIndex != sfi) {
        sfi = tp3d.SFIndex;
        auto& sf = pfp.SectionFits[sfi];
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          plnTP[plane] = MakeBareTP(slc, sf.Pos, sf.Dir, inCTP);
        } // plane
      } // new section
      // check after the previous point is defined
      unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
      if(plane < ipt) {
        // check the separation^2
        auto& fromTP3D = pfp.TP3Ds[prevPt[plane]];
        float sep2 = PosSep2(fromTP3D.Pos, tp3d.Pos);
        if(sep2 > minGapLength2) {
          // do a more careful check. Find the start and end points in 2D in this plane
          auto& tp0 = slc.tjs[fromTP3D.TjID - 1].Pts[fromTP3D.TPIndex];
          auto& tp1 = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
          if(tp0.CTP != tp1.CTP) {
            std::cout<<"AMP: Coding error\n";
            exit(1);
          }
          // check the wire separation and correct for dead wires
          float dw = std::abs(tp0.Pos[0] - tp1.Pos[0]);
          float dwc = DeadWireCount(slc, tp0, tp1);
          float missed = dw - dwc;
          if(missed > 3) {
            std::cout<<"AMP: P"<<pfp.ID<<" check missed in range ";
            std::cout<<prevPt[plane]<<" "<<ipt<<" in CTP "<<tp3d.CTP;
            std::cout<<" "<<PrintPos(slc, tp0)<<" - "<<PrintPos(slc, tp1);
            std::cout<<" dw "<<dw<<" dwc "<<dwc;
            std::cout<<" >>>>>>>. write code to do something here";
            std::cout<<"\n";
          } // missed > 3
        } // sep2 > 25
      } // prevTP3D defined
      prevPt[plane] = ipt;
    } // ipt
    
  } // FillGaps3D

  /////////////////////////////////////////
  void AddPointsInRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, 
                        unsigned short inPlane, float maxPull, bool prt)
  {
    // Try to insert 2D trajectory points into the 3D trajectory point vector pfp.TP3Ds. This function inserts
    // new Tp3Ds and sets the NeedsUpdate flags true. The calling function should call Update
    if(fromPt > toPt) return;
    if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
    if(pfp.SectionFits.empty()) return;
    
    // put the new points into a temporary vector so we can decide if adding them
    // makes sense
    std::vector<TP3D> temp;
    
    bool added = false;
    // max deviation in 2D (WSE units) for considering a hit to be "close" in time. This
    // is a rough cut that is made before a TP3D is created and the pull cut applied
    float maxDelta = 5;
    
    // find the range of wires in each plane that need to be considered
    unsigned int fromWire, toWire;
    auto& fromTP3D = pfp.TP3Ds[fromPt];
    auto& toTP3D = pfp.TP3Ds[toPt - 1];

    // flag the TPs that are already in this PFP.
    // This requires clearing the flag on ALL TPs first
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      for(auto& tp : tj.Pts) tp.Environment[kEnvFlag] = false;
    } // tj
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.TjID <= 0) continue;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      tp.Environment[kEnvFlag] = true;
    } // tp3d
    
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      // restrict the search to one plane?
      if(inPlane < slc.nPlanes && plane != inPlane) continue;
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      // TODO: Note that no check is made for the SectionFit. The assumption is
      // that the trajectory doesn't change significantly between sections
      auto tp = MakeBareTP(slc, fromTP3D.Pos, fromTP3D.Dir, inCTP);
      fromWire = std::nearbyint(tp.Pos[0]);
      tp = MakeBareTP(slc, toTP3D.Pos, toTP3D.Dir, inCTP);
      toWire = std::nearbyint(tp.Pos[0]);
      if(prt) mf::LogVerbatim("TC")<<"APIR: checking wire range "<<fromWire<<" -> "<<toWire;
      if(toWire < fromWire) std::swap(fromWire, toWire);
      for(float wire = fromWire; wire <= toWire; ++wire) {
        // ignore this wire if there is a point on it in this PFP and CTP
        bool skipit = false;
        for(auto& tp3d : pfp.TP3Ds) {
          if(tp3d.Wire == wire && tp3d.CTP == inCTP) {
            skipit = true;
            break;
          }
        } // tp3d
        if(skipit) continue;
        MoveTPToWire(tp, wire);
        if(!FindCloseHits(slc, tp, maxDelta, kAllHits)) continue;
        // Either found one or more hits or we are on a dead wire
        if(tp.Hits.empty()) continue;
        // inspect the hits. Assume that there is only one
        unsigned int iht = tp.Hits[0];
        if(tp.Hits.size() > 1) {
          // Use the closest hit if there are several
          float best = maxDelta;
          for(auto iiht : tp.Hits) {
            // ignore hits that are used in a PFP
            if(slc.slHits[iiht].InPFP > 0) continue;
            float delta = PointTrajDOCA(slc, iiht, tp);
            if(delta > best) continue;
            best = delta;
            iht = iiht;
          } // iiht
          continue;
        } // > 1 hit
        // ignore hits that are used in a PFP
        if(slc.slHits[iht].InPFP > 0) continue;
        // Is used in a trajectory?
        int inTj = slc.slHits[iht].InTraj;
        if(inTj <= 0) {
          // an unused hit. Create a new TP3D
          auto newTP3D = CreateTP3D(slc, iht);
          SetSection(slc, pfp, newTP3D);
          float pull = (newTP3D.Pos[0] - newTP3D.TPX) / sqrt(newTP3D.TPXErr2);
          if(pull > maxPull) continue;
          temp.push_back(newTP3D);
        } else {
          // hit is used in a trajectory. Find the TP index in the trajectory
          unsigned short tjPt = USHRT_MAX;
          auto& tj = slc.tjs[inTj - 1];
          for(tjPt = tj.EndPt[0]; tjPt <= tj.EndPt[1]; ++tjPt) {
            if(std::find(tj.Pts[tjPt].Hits.begin(), tj.Pts[tjPt].Hits.end(), iht) != tj.Pts[tjPt].Hits.end()) break;
          } // tjPt
          if(tjPt == USHRT_MAX) continue;
          auto& tpm = tj.Pts[tjPt];
          // see if the TP is already used in this function
          if(tpm.Environment[kEnvFlag]) continue;
          // or if used in a different PFP
          if(tpm.Environment[kEnvInPFP]) continue;
          auto newTP3D = CreateTP3D(slc, inTj, tjPt);
          SetSection(slc, pfp, newTP3D);
          float pull = (newTP3D.Pos[0] - newTP3D.TPX) / sqrt(newTP3D.TPXErr2);
          if(pull > maxPull) continue;
          tpm.Environment[kEnvFlag] = true;
          temp.push_back(newTP3D);
        } // hit is used in a trajectory
      } // wire
    } // plane
    
    if(temp.empty()) return;
    
    // make a list of Tjs that are not known to belong to this pfp but have a significant fraction 
    // of points used in this PFP
    std::vector<int> tjList;
    for(auto& tp3d : temp) {
      // ignore single hits for this check
      if(tp3d.TjID <= 0) continue;
      // ignore the known Tjs
      if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tp3d.TjID) != pfp.TjIDs.end()) continue;
      if(std::find(tjList.begin(), tjList.end(), tp3d.TjID) != tjList.end()) continue;
      // We found a Tj point that maybe should be associated with this pfp. Determine
      // the fraction of points in this Tj that are used in this PFP
      auto& tj = slc.tjs[tp3d.TjID - 1];
      float cnt = 0;
      float use = 0;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        ++cnt;
        if(tp.Environment[kEnvFlag]) ++use;
      }
      if(cnt == 0) continue;
      float fracUsed = use/cnt;
      // ignore long Tjs with low fracUsed
      if(cnt > 10 && fracUsed < 0.1) continue;
      tjList.push_back(tp3d.TjID);
    } // tp3d
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"APIR: P"<<pfp.ID<<" Found "<<temp.size()<<" candidate points:";
      for(auto tid : tjList) {
        auto& tj = slc.tjs[tid - 1];
        float cnt = 0;
        float use = 0;
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          auto& tp = tj.Pts[ipt];
          if(tp.Chg <= 0) continue;
          ++cnt;
          if(tp.Environment[kEnvFlag]) ++use;
        }
        myprt<<" T"<<tid<<" cnt "<<(int)cnt<<" fracUsed "<<std::setprecision(2)<<use/cnt;
      } // tid
    } // prt
    // Add the TP3Ds that are in the list
    for(auto& tp3d : temp) {
      // don't add the ones that aren't in the list
      if(tp3d.TjID > 0 && std::find(tjList.begin(), tjList.end(), tp3d.TjID) == tjList.end()) continue;
//      float pull = PointPull(slc, pfp, tp3d);
/*
      // We need to find the tj that is in the same CTP as the tp3d
      int useTj = 0;
      for(auto tjid : pfp.TjIDs) {
        if(slc.tjs[tjid - 1].CTP != slc.tjs[tp3d.TjID - 1].CTP) continue;
        useTj = tjid;
        break;
      } // tjid
      if(useTj == 0) continue;
      if(tp3d.slHitsIndex < slc.slHits.size()) {
        if(slc.slHits[tp3d.slHitsIndex].InPFP > 0) continue;
        // single hit
        // check the pull
        float pull = -1;
        if(useTj > 0) {
          auto& tj = slc.tjs[useTj - 1];
          auto& hit = (*evt.allHits)[slc.slHits[tp3d.slHitsIndex].allHitsIndex];
          Point2_t pos;
          pos[0] = hit.WireID().Wire;
          pos[1] = hit.PeakTime() * tcc.unitsPerTick;
          pull = PointPull(slc, pos, hit.Integral(), tj);
        } // useTj > 0
        if(pull > maxPull) continue;
        // associate the hit with this PFP
        slc.slHits[tp3d.slHitsIndex].InPFP = pfp.ID;
        if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" add hit "<<tp3d.slHitsIndex<<" "<<PrintHit(slc.slHits[tp3d.slHitsIndex])<<" pull "<<pull<<" in section "<<tp3d.SFIndex;
        if(!InsertTP3D(pfp, tp3d)) {
          std::cout<<"APIR: insert failed\n";
          continue;
        }
        added = true;
      } else {
        // don't add the ones that aren't in the list
        if(std::find(tjList.begin(), tjList.end(), tp3d.TjID) == tjList.end()) continue;
        // TP
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        // check the pull
        float pull = 100;
        if(useTj > 0) {
          auto& tj = slc.tjs[useTj - 1];
          pull = PointPull(slc, tp.Pos, tp.Chg, tj);
        } // useTj > 0
        if(pull > maxPull) continue;
        if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" add TP T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<" "<<PrintPos(slc, tp)<<" in section "<<tp3d.SFIndex;
        if(!InsertTP3D(pfp, tp3d)) continue;
        added = true;
      }
*/
      if(!InsertTP3D(pfp, tp3d)) continue;
      added = true;
    } // tp3d
    if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" added TP3Ds? "<<added;
    if(!pfp.NeedsUpdate && added) pfp.NeedsUpdate = true;
  } // AddPointsInRange

  /////////////////////////////////////////
  bool InsertTP3D(PFPStruct& pfp, TP3D& tp3d)
  {
    // inserts the tp3d into the correct section
    if(tp3d.SFIndex >= pfp.SectionFits.size()) return false;
    unsigned short ipt = 0;
    for(ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) if(tp3d.SFIndex == pfp.TP3Ds[ipt].SFIndex) break;
    if(ipt == pfp.TP3Ds.size()) return false;
    pfp.TP3Ds.insert(pfp.TP3Ds.begin() + ipt, tp3d);
    pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
    return true;
  } // InsertTP3D

  /////////////////////////////////////////
  bool SortSection(PFPStruct& pfp, unsigned short sfIndex)
  {
    // sorts the TP3Ds by the distance from the start of a fit section
    
    if(sfIndex > pfp.SectionFits.size() - 1) return false;
    auto& sf = pfp.SectionFits[sfIndex];
    if(sf.Pos[0] == 0.0 && sf.Pos[1] == 0.0 && sf.Pos[2] == 0.0) {
      std::cout<<"P"<<pfp.ID<<" section fit start position not defined\n";
      return false;
    }
    
    // a temp vector of points in this section
    std::vector<TP3D> temp;
    // and the index into TP3Ds
    std::vector<unsigned short> indx;
    for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      auto& tp3d = pfp.TP3Ds[ii];
      if(tp3d.SFIndex != sfIndex) continue;
      temp.push_back(tp3d);
      indx.push_back(ii);
    } // tp3d
    if(temp.empty()) return false;
    // no sort needed
    if(temp.size() == 1) return true;
    
    // see if the points are not-contiguous
    bool contiguous = true;
    for(unsigned short ipt = 1; ipt < indx.size(); ++ipt) {
      if(indx[ipt] != indx[ipt - 1] + 1) {
        contiguous = false;
        std::cout<<"SortSection: Points aren't contiguous. Need to re-sort everything\n";
        for(unsigned short ii = 0; ii < indx.size(); ++ii) std::cout<<ii<<" "<<indx[ii]<<"\n";
      }
    }
    if(!contiguous) {
      return false;
    }
    
    std::vector<SortEntry> sortVec(temp.size());
    for(unsigned short ii = 0; ii < temp.size(); ++ii) {
      sortVec[ii].index = ii;
      sortVec[ii].val = temp[ii].along;
//      sortVec[ii].val = PosSep(sf.Pos, temp[ii].Pos);
    } // ipt
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    for(unsigned short ii = 0; ii < temp.size(); ++ii) {
      // overwrite the tp3d
      auto& tp3d = pfp.TP3Ds[indx[ii]];
      tp3d = temp[sortVec[ii].index];
    } // ii
    sf.NeedsUpdate = false;
    return true;
  } // SortSection
  
  /////////////////////////////////////////
  void CountOrder(TCSlice& slc, int tid, const std::vector<TP3D>& tp3ds, unsigned short& nNeg, unsigned short& nPos)
  {
    // returns a count of Tj points that are in negative and positive order
    nNeg = 0; nPos = 0;
    if(tid <= 0) return;
    if(tp3ds.empty()) return;
    unsigned short prevTPIndex = USHRT_MAX;
    for(auto& tp3d : tp3ds) {
      if(!tp3d.IsGood) continue;
      if(tp3d.TjID != tid) continue;
      if(prevTPIndex == USHRT_MAX) {
        prevTPIndex = tp3d.TPIndex;
      } else {
        if(tp3d.TPIndex > prevTPIndex) { ++nPos; } else { ++nNeg; }
        prevTPIndex = tp3d.TPIndex;
      }
    } // tp3d
  } // CountOrder
  
  /////////////////////////////////////////
  void FillTP3Ds(TCSlice& slc, PFPStruct& pfp)
  {
    // Fill the TP3Ds vector. This function is called after FindCompleteness has
    // done a fit of the 3D TP intersections and put the results in SectionFits[0] Pos and Dir.
    if(!pfp.TP3Ds.empty() || pfp.SectionFits.size() != 1) {
      std::cout<<"FillTP3Ds: invalid call P"<<pfp.ID<<". TP3Ds is not empty or SectionFits size != 1\n";
      return;
    } // error
    auto& sf = pfp.SectionFits[0];
    for(auto tid : pfp.TjIDs) {
      auto& tj = slc.tjs[tid - 1];
      // Project the 3D fit parameters found in FindCompleteness into this plane
      TrajPoint ptp = MakeBareTP(slc, sf.Pos, sf.Dir, tj.CTP);
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(tp.Environment[kEnvInPFP]) continue;
        auto tp3d = CreateTP3D(slc, tid, ipt);
        SetSection(slc, pfp, tp3d);
/*
        tp3d.IsGood = true;
        tp3d.SFIndex = 0;
        // find the point of closest approach in 2D. Note that t is
        // re-defined in FitSection to be the closest wire
        double t = (double)(tp.Pos[0] - ptp.Pos[0]) * ptp.Dir[0] + (double)(tp.Pos[1] - ptp.Pos[1]) * ptp.Dir[1];
        double dp0 = ptp.Pos[0] + t * ptp.Dir[0] - tp.Pos[0];
        double dp1 = ptp.Pos[1] + t * ptp.Dir[1] - tp.Pos[1];
        double delta = tcc.wirePitch * (float)sqrt(dp0 * dp0 + dp1 * dp1);
        // find the position in 3D by correcting t (constructed in 2D) for the direction vector projection in this plane.
        // This was stashed in ptp.Delta by MakeBareTP
        t /= delta;
        tp3d.along = t;
        for(unsigned short xyz = 0; xyz < 3; ++xyz) tp3d.Pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
*/
        pfp.TP3Ds.push_back(tp3d);
      } // ipt
    } // tid
    sf.NPts = pfp.TP3Ds.size();
  } // FillTP3Ds
  
  /////////////////////////////////////////
  void Reverse(TCSlice& slc, PFPStruct& pfp)
  {
    // reverse the PFParticle. Note that there is no need to reverse the TP3Ds vector.
    std::reverse(pfp.SectionFits.begin(), pfp.SectionFits.end());
    // swap the start and end points of each SectionFit
    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      // flip the direction vector
      for(unsigned short xyz = 0; xyz < 3; ++xyz) sf.Dir[xyz] *= -1;
      // correct the along variable
      float maxAlong = -1E6;
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.SFIndex != sfi) continue;
        if(tp3d.along > maxAlong) maxAlong = tp3d.along;
      } // tp3d
      for(auto& tp3d : pfp.TP3Ds) {
        if(tp3d.SFIndex != sfi) continue;
        tp3d.along = maxAlong - tp3d.along;
      } // tp3d
    } // sf
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    std::swap(pfp.EndFlag[0], pfp.EndFlag[1]);
  } // Reverse
  
  /////////////////////////////////////////
  void UpdateMatchStructs(TCSlice& slc, int oldTj, int newTj)
  {
    // Replaces tjid and ipt references in slc.pfps from
    // oldTj to newTj. This function is called when Tjs are split or merged
    // or if a Tj is reversed (in which case oldTj = newTj).
    // The method used is to match the trajectory point positions
    if(oldTj <= 0 || oldTj > (int)slc.tjs.size()) return;
    if(newTj <= 0 || newTj > (int)slc.tjs.size()) return;
    if(slc.mallTraj.empty() && slc.pfps.empty()) return;
    
    // convert from int to unsigned short
    int oldtjid = oldTj;
    int newtjid = newTj;
    auto& ntj = slc.tjs[newTj - 1];
    unsigned short npts = ntj.EndPt[1] - ntj.EndPt[0] + 1;
    // put the X positions of the new Tj into a vector for matching
    std::vector<float> xpos(ntj.Pts.size());
    geo::PlaneID planeID = DecodeCTP(ntj.CTP);
    for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
      auto& ntp = ntj.Pts[npt];
      if(ntp.Chg <= 0) continue;
      xpos[npt] = tcc.detprop->ConvertTicksToX(ntp.Pos[1]/tcc.unitsPerTick, planeID);
    } // npt
    
    if(!slc.mallTraj.empty()) {
      for(unsigned int ipt = 0; ipt < slc.mallTraj.size(); ++ipt) {
        auto& tj2pt = slc.mallTraj[ipt];
        if(tj2pt.id > slc.tjs.size()) continue;
        if(tj2pt.id != oldtjid) continue;
        // Found the old Tj. Now find the point
        for(unsigned short npt = ntj.EndPt[0]; npt <= ntj.EndPt[1]; ++npt) {
          auto& ntp = ntj.Pts[npt];
          if(ntp.Chg <= 0) continue;
          if(std::nearbyint(ntp.Pos[0]) == tj2pt.wire && xpos[npt] > tj2pt.xlo && xpos[npt] < tj2pt.xhi) {
            tj2pt.id = newtjid;
            tj2pt.ipt = npt;
            tj2pt.npts = npts;
            break;
          } // points match
        } // npt
      } // ipt
    } // !slc.mallTraj.empty()
  } // UpdateMatchStructs

  /////////////////////////////////////////
  void FillmAllTraj(TCSlice& slc) 
  {
    // Fills the mallTraj vector with trajectory points in the tpc and sorts
    // them by increasing X
    
    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;
    
    // count the number of TPs and clear out any old 3D match flags
    unsigned int ntp = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      if(tj.ID <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      ntp += NumPtsWithCharge(slc, tj, false);
      tj.AlgMod[kMat3D] = false;
    } // tj
    if(ntp < 2) return;
    
    slc.mallTraj.resize(ntp);
    
    // define mallTraj
    unsigned int icnt = 0;
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled] || tj.AlgMod[kHaloTj]) continue;
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      if((int)planeID.Cryostat != cstat) continue;
      if((int)planeID.TPC != tpc) continue;
      int plane = planeID.Plane;
      int tjID = tj.ID;
      if(tjID <= 0) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg == 0) continue;
        if(icnt > slc.mallTraj.size() - 1) break;
        if(tp.Pos[0] < -0.4) continue;
        slc.mallTraj[icnt].wire = std::nearbyint(tp.Pos[0]);
        bool hasWire = tcc.geom->HasWire(geo::WireID(cstat, tpc, plane, slc.mallTraj[icnt].wire));
        // don't try matching if the wire doesn't exist
        if(!hasWire) continue;
        float xpos = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, plane, tpc, cstat);
        float posPlusRMS = tp.Pos[1] + TPHitsRMSTime(slc, tp, kUsedHits);
        float rms = tcc.detprop->ConvertTicksToX(posPlusRMS/tcc.unitsPerTick, plane, tpc, cstat) - xpos;
        if(rms < tcc.match3DCuts[0]) rms = tcc.match3DCuts[0];
        if(icnt == slc.mallTraj.size()) break;
        slc.mallTraj[icnt].xlo = xpos - rms;
        slc.mallTraj[icnt].xhi = xpos + rms;
        slc.mallTraj[icnt].dir = tp.Dir;
        slc.mallTraj[icnt].ctp = tp.CTP;
        slc.mallTraj[icnt].id = tjID;
        slc.mallTraj[icnt].ipt = ipt;
        slc.mallTraj[icnt].npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        ++icnt;
        if(icnt == ntp) break;
      } // tp
    } // tj
    
    if(icnt < slc.mallTraj.size()) slc.mallTraj.resize(icnt);
    
    // sort by increasing x
    std::vector<SortEntry> sortVec(slc.mallTraj.size());
    for(unsigned int ipt = 0; ipt < slc.mallTraj.size(); ++ipt) {
      // populate the sort vector
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = slc.mallTraj[ipt].xlo;
    } // ipt
    // sort by increasing xlo
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    // put slc.mallTraj into sorted order
    auto tallTraj = slc.mallTraj;
    for(unsigned int ii = 0; ii < sortVec.size(); ++ii) slc.mallTraj[ii] = tallTraj[sortVec[ii].index];
    
  } // FillmAllTraj

  /////////////////////////////////////////
  bool SharesHighScoreVx(TCSlice& slc, const PFPStruct& pfp, const Trajectory& tj)
  {
    // returns true if tj with tjID shares a high-score 3D vertex with any
    // tj in pfp.TjIDs
    for(unsigned short end = 0; end < 2; ++end) {
      if(tj.VtxID[end] == 0) continue;
      auto& vx2 = slc.vtxs[tj.VtxID[end] - 1];
      if(!vx2.Stat[kHiVx3Score]) continue;
      std::vector<int> vtjlist = GetVtxTjIDs(slc, vx2);
      auto shared = SetIntersection(vtjlist, pfp.TjIDs);
      if(!shared.empty()) return true;
    } // end
    return false;
  } // SharesHighScoreVx

  /////////////////////////////////////////
  void Fit3D(unsigned short mode, Point3_t point, Vector3_t dir, Point3_t& fitPos, Vector3_t& fitDir)
  {
    // initialize, accumulate and fit the points. The code to fit the direction using the positions
    // of the points is commented out and replaced with a simple average of the directions of the points
    
    // 3D fit sum variables
    static double fSum, fSumx, fSumy, fSumz;
//    static double fSumx2, fSumy2, fSumz2, fSumxz, fSumyz;
    // average the direction vectors
    static double fSumDir0, fSumDir1, fSumDir2;

    if(mode == 0) {
      fSum = 0; fSumx = 0; fSumy = 0; fSumz = 0; 
//      fSumx2 = 0; fSumy2 = 0; fSumz2 = 0; fSumxz = 0; fSumyz = 0;
      fSumDir0 = 0; fSumDir1 = 0; fSumDir2 = 0;
      return;
    }
    // accumulate
    if(mode == 1) {
      fSum += 1;
      fSumx += point[0];
      fSumy += point[1];
      fSumz += point[2];
      fSumDir0 += dir[0];
      fSumDir1 += dir[1];
      fSumDir2 += dir[2];
      return;
    }
    
    if(fSum < 2) return;
    // just use the average for the position
    fitPos[0] = fSumx / fSum;
    fitPos[1] = fSumy / fSum;
    fitPos[2] = fSumz / fSum;
    // and for the direction
    fitDir = {{fSumDir0, fSumDir1, fSumDir2}};
    SetMag(fitDir, 1);
  } // Fit3D
/*
  /////////////////////////////////////////
  unsigned short WiresSkippedInCTP(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP)
  {
    // counts the number of wires between the end points of all Tjs in the list of tjids
    // in inCTP where there is no TP with charge
    if(tjids.empty()) return 0;
    
    // find the min and max Pos[0] positions
    float fLoWire = 1E6;
    float fHiWire = -1E6;
    for(auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.CTP != inCTP) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        float endWire = tj.Pts[tj.EndPt[end]].Pos[0];
        if(endWire < -0.4) continue;
        if(endWire < fLoWire) fLoWire = endWire;
        if(endWire > fHiWire) fHiWire = endWire;
      } // end
    } // tjid
    if(fLoWire >= fHiWire) return 0;
    unsigned int loWire = std::nearbyint(fLoWire);
    unsigned short nWires = std::nearbyint(fHiWire) - loWire + 1;
    std::vector<bool> ptOnWire(nWires, false);
    
    // count the number of points with charge on all Tjs
    float npwc = 0;
    for(auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.CTP != inCTP) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(tp.Pos[0] < -0.4) continue;
        ++npwc;
        unsigned short indx = std::nearbyint(tp.Pos[0]) - loWire;
        if(indx < nWires) ptOnWire[indx] = true;
      } // ipt
    } // tjid
    if(npwc == 0) return 0;
    float nskip = 0;
    for(unsigned short indx = 0; indx < nWires; ++indx) if(!ptOnWire[indx]) ++nskip;
    return (nskip / npwc);
    
  } // WiresSkippedInCTP
*/
/*
  /////////////////////////////////////////
  float LengthInCTP(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP)
  {
    // Calculates the maximum length between the end points of Tjs in the list of tjids in inCTP
    if(tjids.empty()) return 0;
    // put the end point positions into a vector
    std::vector<Point2_t> endPos;
    for(auto tjid : tjids) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.CTP != inCTP) continue;
      endPos.push_back(tj.Pts[tj.EndPt[0]].Pos);
      endPos.push_back(tj.Pts[tj.EndPt[1]].Pos);
    } // tjid
    if(endPos.size() < 2) return 0;
    float extent = 0;
    for(unsigned short pt1 = 0; pt1 < endPos.size() - 1; ++pt1) {
      for(unsigned short pt2 = pt1 + 1; pt2 < endPos.size(); ++pt2) {
        float sep = PosSep2(endPos[pt1], endPos[pt2]);
        if(sep > extent) extent = sep;
      } // pt2
    } // pt1
    return sqrt(extent);
  } // LengthInCTP
*/
  /////////////////////////////////////////
  bool MakeTp3(TCSlice& slc, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3, bool findDirection)
  {
    // Make a 3D trajectory point using two 2D trajectory points
    tp3.Dir = {{999.0, 999.0, 999.0}};
    tp3.Pos = {{999.0, 999.0, 999.0}};
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);
    double upt = tcc.unitsPerTick;
    double ix = tcc.detprop->ConvertTicksToX(itp.Pos[1] / upt, iPlnID);
    double jx = tcc.detprop->ConvertTicksToX(jtp.Pos[1] / upt, jPlnID);
    
    // don't continue if the points are wildly far apart in X
    if(std::abs(ix - jx) > 20) return false;
    tp3.Pos[0] = 0.5 * (ix + jx);
    // determine the wire orientation and offsets using WireCoordinate
    // wire = yp * OrthY + zp * OrthZ - Wire0 = cs * yp + sn * zp - wire0
    // wire offset
    double iw0 = tcc.geom->WireCoordinate(0, 0, iPlnID);
    // cosine-like component
    double ics = tcc.geom->WireCoordinate(1, 0, iPlnID) - iw0;
    // sine-like component
    double isn = tcc.geom->WireCoordinate(0, 1, iPlnID) - iw0;
    double jw0 = tcc.geom->WireCoordinate(0, 0, jPlnID);
    double jcs = tcc.geom->WireCoordinate(1, 0, jPlnID) - jw0;
    double jsn = tcc.geom->WireCoordinate(0, 1, jPlnID) - jw0;
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
    
    if(!findDirection) return true;

    // make a copy of itp and shift it by many wires to avoid precision problems
    double itp2_0 = itp.Pos[0] + 100;
    double itp2_1 = itp.Pos[1];
    if(std::abs(itp.Dir[0]) > 0.01) itp2_1 += 100 * itp.Dir[1] / itp.Dir[0];
    // Create a second Point3 for the shifted point
    Point3_t pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = tcc.detprop->ConvertTicksToX(itp2_1 / upt, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = tcc.detprop->ConvertXToTicks(pos2[0], jPlnID) * upt;
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
    if(v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) return 0;
    return acos(DotProd(v1, v2));
  } 
  
  ////////////////////////////////////////////////
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2)
  {
    // Finds the direction vector between the two points from p1 to p2
    Vector3_t dir;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) dir[xyz] = p2[xyz] - p1[xyz];
    if(dir[0] == 0 && dir[1] == 0 && dir[2] == 0) return dir;
    if(!SetMag(dir, 1)) { dir[0] = 0; dir[1] = 0; dir[3] = 0; }
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

  /////////////////////////////////////////
  void FilldEdx(TCSlice& slc, PFPStruct& pfp)
  {
    // Fills dE/dx variables in the pfp struct
    // TODO: Do this correctly instead of just quickly
    
    // don't attempt to find dE/dx at the end of a shower
    unsigned short numEnds = 2;
    if(pfp.PDGCode == 1111) numEnds = 1;
    
    // set dE/dx to 0 to indicate that a valid dE/dx is expected
    for(unsigned short end = 0; end < numEnds; ++end) {
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) pfp.dEdx[end][plane] = 0;
    } // end
    
//    Fit2D();
    float rngOff = tcc.wirePitch / 2;
    
    // find dE/dx for points in each plane at both ends
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      for(unsigned short end = 0; end < 2; ++end) {
        short dir = 1 - 2 * end;
        auto endPos = PosAtEnd(pfp, end);
        for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
          unsigned short ipt;
          if(dir > 0) {
            ipt = ii;
          } else {
            ipt = pfp.TP3Ds.size() - ii - 1;
          }
          if(ipt >= pfp.TP3Ds.size()) break;
          auto& tp3d = pfp.TP3Ds[ipt];
          if(!tp3d.IsGood) continue;
          // ignore single hits
          if(tp3d.TjID <= 0) continue;
          auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
          // look for a TP in this plane
          if(DecodeCTP(tp.CTP).Plane != plane) continue;
          // check for overlap condition
          if(tp.Environment[kEnvOverlap]) continue;
          float dedx = dEdx(slc, tp3d);
          if(dedx < 0.5) continue;
          // find the residual range and add the offset
          float rr = PosSep(tp3d.Pos, endPos) + rngOff;
          std::cout<<"FdEdx: P"<<pfp.ID<<" T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<" rr "<<rr<<" dEdx "<<dedx<<"\n";
          if(ipt == 0) break;
          // break after getting 30 cm away from the end
          if(rr > 30) break;
        } // ii
      } // end
    } // plane
  } // FilldEdx

  /////////////////////////////////////////
  float dEdx(TCSlice& slc, TP3D& tp3d)
  {
    if(!tp3d.IsGood) return 0;
    if(tp3d.TjID <= 0) {
      std::cout<<"Fix dEdx to handle single hits\n";
      return 0;
    }
    auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    double angleToVert = tcc.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
    double cosgamma = std::abs(std::sin(angleToVert) * tp3d.Dir[1] + std::cos(angleToVert) * tp3d.Dir[2]);
    if(cosgamma == 0) return 0;
    double dx = tcc.geom->WirePitch(planeID) / cosgamma;
    double dQ = tp.Chg;
    // Correct for the 1 WSE path length normalization
    if(tp.Dir[0] != 0) {
      double pathInv = std::abs(tp.Dir[0]);
      if(pathInv < 0.05) pathInv = 0.05;
      dQ /= pathInv;
    }
    double dQdx = dQ / dx;
    double time = tp.Pos[1] / tcc.unitsPerTick;
    double t0 = 0;
    return tcc.caloAlg->dEdx_AREA(dQdx, time, planeID.Plane, t0);
  } // dEdx
  
  ////////////////////////////////////////////////
  TP3D CreateTP3D(TCSlice& slc, unsigned short slHitsIndex)
  {
    // create a TP3D with a single hit. Note that the SectionFit in which it
    // should be placed and the 3D position can't be determined until the the TP3D is 
    // associated with a pfp. See SetSection()
    TP3D tp3d;
    if(slHitsIndex < slc.slHits.size()) {
      tp3d.slHitsIndex = slHitsIndex;
      auto& hit = (*evt.allHits)[slc.slHits[tp3d.slHitsIndex].allHitsIndex];
      auto plnID = geo::PlaneID(slc.TPCID.Cryostat, slc.TPCID.TPC, hit.WireID().Plane);
      tp3d.CTP = EncodeCTP(plnID);
      tp3d.TPX = tcc.detprop->ConvertTicksToX(hit.PeakTime(), plnID);
      double tickErr = hit.RMS() * tcc.hitErrFac * hit.Multiplicity();
      double xdx = tcc.detprop->ConvertTicksToX(hit.PeakTime() + tickErr, plnID);
      // inflate the errors for single hit. It must be poor quality since it isn't used in a TP
      double xErr = 2 * (xdx - tp3d.TPX);
      tp3d.TPXErr2 = xErr * xErr;
      tp3d.Wire = hit.WireID().Wire;
      // Can't declare it good since Pos and SFIndex aren't defined
      tp3d.IsGood = false;
    } // valid slHitsIndex
    return tp3d;
  } // CreateTP3D
  
  ////////////////////////////////////////////////
  TP3D CreateTP3D(TCSlice& slc, int tjID, unsigned short tpIndex)
  {
    // create a TP3D with a single TP. Note that the SectionFit in which it
    // should be placed and the 3D position can't be determined until the the TP3D is 
    // associated with a pfp. See SetSection()
    
    TP3D tp3d;    
    // add the 2D TP
    if(tjID > 0 && tjID <= slc.tjs.size()) {
      tp3d.TjID = tjID;
      auto& tj = slc.tjs[tp3d.TjID - 1];
      if(tpIndex < tj.Pts.size()) {
        tp3d.TPIndex = tpIndex;
        auto& tp2 = tj.Pts[tp3d.TPIndex];
        auto plnID = DecodeCTP(tp2.CTP);
        tp3d.CTP = tp2.CTP;
        tp3d.TPX = tcc.detprop->ConvertTicksToX(tp2.Pos[1]/tcc.unitsPerTick, plnID);
        // The TP2 hit error is in WSE units and includes the wire error. See
        // StepUtils/DefineHitPos. Ignore these details for now and fake it.
        double tickErr = tp2.HitPosErr2 / (tcc.unitsPerTick * tcc.unitsPerTick);
        double xdx = tcc.detprop->ConvertTicksToX(tp2.Pos[1]/tcc.unitsPerTick + tickErr, plnID);
        double xErr = xdx - tp3d.TPX;
        // increase the error if the TP overlaps with another trajectory
        if(tp2.Environment[kEnvOverlap]) xErr *= 2;
        // TODO: increase the error if this near an end that has a Bragg peak. How is this
        // known here?
        tp3d.TPXErr2 = xErr * xErr;
        tp3d.Wire = tp2.Pos[0];
        // Can't declare it good since Pos and SFIndex aren't defined
        tp3d.IsGood = false;
      } // valid tpIndex
    } // valid tjID
    return tp3d;
  } // CreateTP3D

  /////////////////////////////////////////
  void SetSection(TCSlice& slc, PFPStruct& pfp, TP3D& tp3d)
  {
    // Determine which SectionFit this tp3d should reside, then calculate
    // the 3D position and the distance fromt the center of the SectionFit
    tp3d.IsGood = false;
    if(tp3d.TPIndex == USHRT_MAX && tp3d.SFIndex == USHRT_MAX) return;
    if(tp3d.Wire < 0) return;
    if(pfp.SectionFits.empty()) return;
    
    auto plnID = DecodeCTP(tp3d.CTP);
    
    if(pfp.SectionFits.size() == 1) {
      tp3d.SFIndex = 0;
    } else {
      // Find the section center that is closest to this point in the wire coordinate
      float best = 1E6;
      for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        auto& sf = pfp.SectionFits[sfi];
        float sfWire = tcc.geom->WireCoordinate(sf.Pos[1], sf.Pos[2], plnID);
        float sep = std::abs(sfWire - tp3d.Wire);
        if(sep < best) {
          best = sep;
          tp3d.SFIndex = sfi;
        }
      } // sfi
    } // pfp.SectionFits.size() > 1
    auto& sf = pfp.SectionFits[tp3d.SFIndex];
    auto plnTP = MakeBareTP(slc, sf.Pos, sf.Dir, tp3d.CTP);
    // the number of wires relative to the SectionFit center
    double dw = tp3d.Wire - plnTP.Pos[0];
    // dt/dW was stored in DeltaRMS
    double t = dw * plnTP.DeltaRMS;
    // define the 3D position 
    for(unsigned short xyz = 0; xyz < 3; ++xyz) tp3d.Pos[xyz] = sf.Pos[xyz] + t * sf.Dir[xyz];
    tp3d.along = t;
    tp3d.IsGood = true;
  } // SetSection

  ////////////////////////////////////////////////
  PFPStruct CreatePFP(TCSlice& slc)
  {
    // The calling function should define the size of pfp.TjIDs
    PFPStruct pfp;
    pfp.ID = slc.pfps.size() + 1;
    pfp.ParentUID = 0;
    pfp.TPCID = slc.TPCID;
    // initialize arrays for both ends
    if(slc.nPlanes < 4) {
      pfp.dEdx[0].resize(slc.nPlanes, -1);
      pfp.dEdx[1].resize(slc.nPlanes, -1);
      pfp.dEdxErr[0].resize(slc.nPlanes, -1);
      pfp.dEdxErr[1].resize(slc.nPlanes, -1);
    }
    // create a single section fit to hold the start/end positions and direction
    pfp.SectionFits.resize(1);
    return pfp;
  } // CreatePFP
  /////////////////////////////////////////
  void PFPVertexCheck(TCSlice& slc)
  {
    // Ensure that all PFParticles have a start vertex. It is possible for
    // PFParticles to be attached to a 3D vertex that is later killed.
    if(!slc.isValid) return;
    
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      if(pfp.Vx3ID[0] > 0) continue;
      if(pfp.SectionFits.empty()) continue;
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      vx3.Vx2ID.resize(slc.nPlanes);
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      auto& startPos = pfp.SectionFits[0].Pos;
      vx3.X = startPos[0];
      vx3.Y = startPos[1];
      vx3.Z = startPos[2];
      vx3.ID = slc.vtx3s.size() + 1;
      vx3.Primary = false;
      ++evt.globalP_UID;
      vx3.UID = evt.globalP_UID;
      slc.vtx3s.push_back(vx3);
//      std::cout<<"PFPVertexCheck: P"<<pfp.ID<<" create 3V"<<vx3.ID<<"\n";
      pfp.Vx3ID[0] = vx3.ID;
    } // pfp
  } // PFPVertexCheck
  
  /////////////////////////////////////////
  void DefinePFPParents(TCSlice& slc, bool prt)
  {
    /*
     This function reconciles vertices, PFParticles and slc, then
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
    if(slc.pfps.empty()) return;
    if(tcc.modes[kTestBeam]) return;
    
    int neutrinoPFPID = 0;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      if(!tcc.modes[kTestBeam] && neutrinoPFPID == 0 && (pfp.PDGCode == 12 || pfp.PDGCode == 14)) {
        neutrinoPFPID = pfp.ID;
        break;
      }
    } // pfp
    
    // define the end vertex if the Tjs have end vertices
    constexpr unsigned short end1 = 1;
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      // already done?
      if(pfp.Vx3ID[end1] > 0) continue;
      // ignore shower-like pfps
      if(IsShowerLike(slc, pfp.TjIDs)) continue;
      // count 2D -> 3D matched vertices
      unsigned short cnt3 = 0;
      unsigned short vx3id = 0;
      // list of unmatched 2D vertices that should be merged
      std::vector<int> vx2ids;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.VtxID[end1] == 0) continue;
        auto& vx2 = slc.vtxs[tj.VtxID[end1] - 1];
        if(vx2.Vx3ID == 0) {
          if(vx2.Topo == 1 && vx2.NTraj == 2) vx2ids.push_back(vx2.ID);
          continue;
        }
        if(vx3id == 0) vx3id = vx2.Vx3ID;
        if(vx2.Vx3ID == vx3id) ++cnt3;
      } // tjid
      if(cnt3 > 1) {
        // ensure it isn't attached at the other end
        if(pfp.Vx3ID[1 - end1] == vx3id) continue;
        pfp.Vx3ID[end1] = vx3id;
        if(cnt3 != slc.nPlanes && tcc.modes[kDebug]) {
          std::cout<<"DPFPR: Missed an end vertex for PFP "<<pfp.ID<<" Write some code\n";
        }
      } // cnt3 > 1
    } // pfp
    
    // Assign a PDGCode to each PFParticle and look for a parent
    for(auto& pfp : slc.pfps) {
      if(pfp.ID == 0) continue;
      // skip a neutrino PFParticle
      if(pfp.PDGCode == 12 || pfp.PDGCode == 14 || pfp.PDGCode == 22) continue;
      pfp.PDGCode = PDGCodeVote(slc, pfp.TjIDs, prt);
      // Define a PFP parent if there are two or more Tjs that are daughters of
      // Tjs that are used by the same PFParticle
      int pfpParentID = INT_MAX;
      unsigned short nParent = 0;
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.ParentID <= 0) continue;
        auto parPFP = GetAssns(slc, "T", tj.ParentID, "P");
        if(parPFP.empty()) continue;
        if(pfpParentID == INT_MAX) pfpParentID = parPFP[0];
        if(parPFP[0] == pfpParentID) ++nParent;
      } // ii
      if(nParent > 1) {
        auto& ppfp = slc.pfps[pfpParentID - 1];
        // set the parent UID
        pfp.ParentUID = ppfp.UID;
        // add to the parent daughters list
        ppfp.DtrUIDs.push_back(pfp.UID);
      } // nParent > 1
    } // ipfp

/* TODO: This needs work
    if(tcc.modes[kTestBeam]) {
      DefinePFPParentsTestBeam(slc, prt);
      return;
    }
*/
    // associate primary PFParticles with a neutrino PFParticle
    if(neutrinoPFPID > 0) {
      auto& neutrinoPFP = slc.pfps[neutrinoPFPID - 1];
      int vx3id = neutrinoPFP.Vx3ID[1];
      for(auto& pfp : slc.pfps) {
        if(pfp.ID == 0 || pfp.ID == neutrinoPFPID) continue;
        if(pfp.Vx3ID[0] != vx3id) continue;
        pfp.ParentUID = (size_t)neutrinoPFPID;
        pfp.Primary = true;
        neutrinoPFP.DtrUIDs.push_back(pfp.ID);
        if(pfp.PDGCode == 111) neutrinoPFP.PDGCode = 12;
      } // pfp
    } // neutrino PFP exists    
  } // DefinePFPParents
/*
  /////////////////////////////////////////
  void DefinePFPParentsTestBeam(TCSlice& slc, bool prt)
  {
    // analog of the one above that was written for neutrino interactions. This differs in that
    // the Tj parent - daughter relationship isn't known yet. If one exists, it is ignored...
    // The assumption here is that all PFParticles that enter (end0) from upstream Z are parents and 
    // any PFParticles attached to them at end1 are daughters. 

    // create a list (stack) of parent ID <-> daughter IDs. The idea is similar to that
    // used in DefineTjParents. A parent-daughter association is made for each entry. After
    // it is made, 1) that entry is removed from the stack, 2) the daughter is checked to see
    // if it a parent of a grand-daughter and if so that pair is added to the stack. 
    std::vector<std::pair<unsigned short, unsigned short>> pardtr;

    // Fill the stack with parents that enter the TPC and have daughters attached to
    // 3D vertices at the other end
    double fidZCut = slc.zLo + 2;
    for(auto& parPFP : slc.pfps) {
      if(parPFP.ID == 0) continue;
      parPFP.Primary = false;
      if(parPFP.XYZ[0][2] > fidZCut) continue;
      parPFP.Primary = true;
      // we found a pfp that entered the TPC. Call it the parent and look for a daughter
      if(prt) mf::LogVerbatim("TC")<<"DPFPTestBeam: parent "<<parPFP.ID<<" end1 vtx "<<parPFP.Vx3ID[1];
      if(parPFP.Vx3ID[1] == 0) continue;
      // There must be other Tjs attached to this vertex which are the daughters. Find them
      // and add them to the pardtr stack
      float score = 0;
      auto& vx3 = slc.vtx3s[parPFP.Vx3ID[1] - 1];
      // ensure that it is valid
      if(vx3.ID == 0) continue;
      // get a list of Tjs attached to this vertex. This will include the Tjs in the parent.
      auto vx3TjList = GetVtxTjIDs(slc, vx3, score);
      if(vx3TjList.empty()) continue;
      // filter out the parent Tjs
      auto dtrTjlist = SetDifference(vx3TjList, parPFP.TjIDs);
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<" Dtrs:";
        for(auto dtjID : dtrTjlist) myprt<<" "<<dtjID<<"_"<<GetPFPIndex(slc, dtjID);
      }
      // Add to the stack
      for(auto dtjID : dtrTjlist) {
        unsigned short pfpIndex = GetPFPIndex(slc, dtjID);
        if(pfpIndex > slc.pfps.size() - 1) continue;
        unsigned short dtrID = pfpIndex + 1;
        // See if this is a duplicate
        bool duplicate = false;
        for(auto& pd : pardtr) if(parPFP.ID == pd.first && dtrID == pd.second) duplicate = true;
        if(!duplicate) pardtr.push_back(std::make_pair(parPFP.ID, dtrID));
      } // dtjID
    } // parPFP
    
    // iterate through the parent - daughter stack, removing the last pair when a 
    // ParentID is updated and adding pairs for new daughters
    for(unsigned short nit = 0; nit < 100; ++nit) {
      if(pardtr.empty()) break;
      auto lastPair = pardtr[pardtr.size() - 1];
      auto& dtr = slc.pfps[lastPair.second - 1];
      auto& par = slc.pfps[lastPair.first - 1];
      dtr.ParentUID = par.UID;
      par.DtrUIDs.push_back(dtr.UID);
      // remove the last pair
      pardtr.pop_back();
      // Now see if the daughter is a parent. First check for a vertex at the other end.
      // To do that we need to know which end has the vertex between the parent and daughter
      unsigned short dtrEnd = USHRT_MAX;
      for(unsigned short ep = 0; ep < 2; ++ep) {
        if(par.Vx3ID[ep] == 0) continue;
        for(unsigned short ed = 0; ed < 2; ++ed) if(dtr.Vx3ID[ed] == par.Vx3ID[ep]) dtrEnd = ed;
      } // ep
      if(dtrEnd > 1) continue;
      // look at the other end of the daughter
      dtrEnd = 1 - dtrEnd;
      // check for a vertex
      if(dtr.Vx3ID[dtrEnd] == 0) continue;
      // get the list of Tjs attached to it
      auto& vx3 = slc.vtx3s[dtr.Vx3ID[dtrEnd] - 1];
      float score = 0;
      auto vx3TjList = GetVtxTjIDs(slc, vx3, score);
      if(vx3TjList.empty()) continue;
      // filter out the new parent
      auto dtrTjlist = SetDifference(vx3TjList, dtr.TjIDs);
      // put these onto the stack
      for(auto tjid : dtrTjlist) pardtr.push_back(std::make_pair(dtr.ID, tjid));
    } // nit
  } // DefinePFPParentsTestBeam
*/
  ////////////////////////////////////////////////
  bool StorePFP(TCSlice& slc, PFPStruct& pfp)
  {
    // stores the PFParticle in the slice
    if(pfp.ID < int(slc.pfps.size())) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    if(!neutrinoPFP) {
      if(pfp.TjIDs.empty()) return false;
      if(pfp.PDGCode != 1111 && pfp.TP3Ds.size() < 2) return false;
    }
    // check the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.AlgMod[kMat3D]) return false;
    } // tjid

    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)slc.pfps.size() + 1) pfp.ID = slc.pfps.size() + 1;
    ++evt.globalP_UID;
    pfp.UID = evt.globalP_UID;
    
    // set the TP used in PFP bit
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.TjID <= 0) continue;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvInPFP]) {
        std::cout<<"StorePFP: Trying to use already-used TP T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"\n";
        return false;
      } // error
      tp.Environment[kEnvInPFP] = true;
      // associate all of the used hits in the TP with this pfp
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) slc.slHits[tp.Hits[ii]].InPFP = pfp.ID;
      } // ii
    } // tp3d
    
    // set the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      tj.AlgMod[kMat3D] = true;
    } // tjid
    
    if(pfp.NeedsUpdate) std::cout<<"StorePFP: stored P"<<pfp.ID<<" but NeedsUpdate is true...\n";
    
    slc.pfps.push_back(pfp);
    return true;
  } // StorePFP
  ////////////////////////////////////////////////
  bool InsideFV(TCSlice& slc, PFPStruct& pfp, unsigned short end)
  {
    // returns true if the end of the pfp is inside the fiducial volume of the TPC
    if(pfp.ID <= 0) return false;
    if(end > 1) return false;
    if(pfp.SectionFits.empty()) return false;
    // require that the points are sorted which ensures that the start and end points
    // are the first and last points in the TP3Ds vector
    if(pfp.NeedsUpdate) return false;
    
    float abit = 5;
    Point3_t pos;
    if(end == 0) {
      pos = pfp.TP3Ds[0].Pos;
    } else {
      pos = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    }
    return (pos[0] > slc.xLo + abit && pos[0] < slc.xHi - abit && 
            pos[1] > slc.yLo + abit && pos[1] < slc.yHi - abit &&
            pos[2] > slc.zLo + abit && pos[2] < slc.zHi - abit);
    
  } // InsideFV
  
  ////////////////////////////////////////////////
  bool InsideTPC(const Point3_t& pos, geo::TPCID& inTPCID)
  {
    // determine which TPC this point is in. This function returns false
    // if the point is not inside any TPC
    float abit = 5;
    for (const geo::TPCID& tpcid: tcc.geom->IterateTPCIDs()) {
      const geo::TPCGeo& TPC = tcc.geom->TPC(tpcid);
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      TPC.LocalToWorld(local,world);
      // reduce the active area of the TPC by a bit to be consistent with FillWireHitRange
      if(pos[0] < world[0]-tcc.geom->DetHalfWidth(tpcid) + abit) continue;
      if(pos[0] > world[0]+tcc.geom->DetHalfWidth(tpcid) - abit) continue;
      if(pos[1] < world[1]-tcc.geom->DetHalfHeight(tpcid) + abit) continue;
      if(pos[1] > world[1]+tcc.geom->DetHalfHeight(tpcid) - abit) continue;
      if(pos[2] < world[2]-tcc.geom->DetLength(tpcid)/2 + abit) continue;
      if(pos[2] > world[2]+tcc.geom->DetLength(tpcid)/2 - abit) continue;
      inTPCID = tpcid;
      return true;
    } // tpcid
    return false;
  } // InsideTPC
  
  ////////////////////////////////////////////////
  void FindAlongTrans(Point3_t pos1, Vector3_t dir1, Point3_t pos2, Point2_t& alongTrans)
  {
    // Calculate the distance along and transvers to the direction vector from pos1 to pos2
    alongTrans[0] = 0; 
    alongTrans[1] = 0;
    if(pos1[0] == pos2[0] && pos1[1] == pos2[1] && pos1[2] == pos2[2]) return;
    auto ptDir = PointDirection(pos1, pos2);
    SetMag(dir1, 1.0);
    double costh = DotProd(dir1, ptDir);
    if(costh > 1) costh = 1;
    double sep = PosSep(pos1, pos2);
    alongTrans[0] = costh * sep;
    double sinth = sqrt(1 - costh * costh);
    alongTrans[1] = sinth * sep;
  } // FindAlongTrans
  
  ////////////////////////////////////////////////
  bool PointDirIntersect(Point3_t p1, Vector3_t p1Dir, Point3_t p2, Vector3_t p2Dir, Point3_t& intersect, float& doca)
  {
    // Point - vector version
    Point3_t p1End, p2End;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) {
      p1End[xyz] = p1[xyz] + 10 * p1Dir[xyz];
      p2End[xyz] = p2[xyz] + 10 * p2Dir[xyz];
    }
    return LineLineIntersect(p1, p1End, p2, p2End, intersect, doca);
  } // PointDirIntersect
  
  ////////////////////////////////////////////////
  bool LineLineIntersect(Point3_t p1, Point3_t p2, Point3_t p3, Point3_t p4, Point3_t& intersect, float& doca)
  {
    /*
     Calculate the line segment PaPb that is the shortest route between
     two lines P1P2 and P3P4. Calculate also the values of mua and mub where
     Pa = P1 + mua (P2 - P1)
     Pb = P3 + mub (P4 - P3)
     Return FALSE if no solution exists.
     http://paulbourke.net/geometry/pointlineplane/
     */

    Point3_t p13, p43, p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;
    constexpr double EPS = std::numeric_limits<double>::min();
    
    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    if (std::abs(p43[0]) < EPS && std::abs(p43[1]) < EPS && std::abs(p43[2]) < EPS) return(false);
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    if (std::abs(p21[0]) < EPS && std::abs(p21[1]) < EPS && std::abs(p21[2]) < EPS) return(false);
    
    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];
    
    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < EPS) return(false);
    numer = d1343 * d4321 - d1321 * d4343;
    
    double mua = numer / denom;
    double mub = (d1343 + d4321 * mua) / d4343;
    
    intersect[0] = p1[0] + mua * p21[0];
    intersect[1] = p1[1] + mua * p21[1];
    intersect[2] = p1[2] + mua * p21[2];
    Point3_t pb;
    pb[0] = p3[0] + mub * p43[0];
    pb[1] = p3[1] + mub * p43[1];
    pb[2] = p3[2] + mub * p43[2];
    doca = PosSep(intersect, pb);
    // average the closest points
    for(unsigned short xyz = 0; xyz < 3; ++xyz) intersect[xyz] += pb[xyz];
    for(unsigned short xyz = 0; xyz < 3; ++xyz) intersect[xyz] /= 2;
    return true;
  } // LineLineIntersect

  ////////////////////////////////////////////////
  float ChgFracBetween(TCSlice& slc, Point3_t pos1, Point3_t pos2)
  {
    // Step between pos1 and pos2 and find the fraction of the points that have nearby hits
    // in each plane. This function returns -1 if something is fishy, but this doesn't mean
    // that there is no charge. Note that there is no check for charge precisely at the pos1 and pos2
    // positions 
    float sep = PosSep(pos1, pos2);
    if(sep == 0) return -1;
    unsigned short nstep = sep / tcc.wirePitch;
    auto dir = PointDirection(pos1, pos2);
    float sum = 0;
    float cnt = 0;
    TrajPoint tp;
    for(unsigned short step = 0; step < nstep; ++step) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) pos1[xyz] += tcc.wirePitch * dir[xyz];
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        tp.CTP = EncodeCTP(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
        tp.Pos[0] = tcc.geom->WireCoordinate(pos1[1], pos1[2], plane, slc.TPCID.TPC, slc.TPCID.Cryostat);
        tp.Pos[1] = tcc.detprop->ConvertXToTicks(pos1[0], plane, slc.TPCID.TPC, slc.TPCID.Cryostat) * tcc.unitsPerTick;
        ++cnt;
        if(SignalAtTp(tp)) ++sum;
      } // plane
    } // step
    if(cnt == 0) return -1;
    return sum / cnt;
    
  } // ChgFracBetween
  
  ////////////////////////////////////////////////
  float ChgFracNearEnd(TCSlice& slc, PFPStruct& pfp, unsigned short end)
  {
    // returns the charge fraction near the end of the pfp. Note that this function
    // assumes that there is only one Tj in a plane.
    if(pfp.ID == 0) return 0;
    if(pfp.TjIDs.empty()) return 0;
    if(end < 0 || end > 1) return 0;
    if(pfp.TPCID != slc.TPCID) return 0;
    if(pfp.SectionFits.empty()) return 0;

    float sum = 0;
    float cnt = 0;
    // keep track of the lowest value and maybe reject it
    float lo = 1;
    float hi = 0;
    auto pos3 = PosAtEnd(pfp, end);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      std::vector<int> tjids(1);
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.CTP != inCTP) continue;
        tjids[0] = tjid;
        Point2_t pos2;
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        pos2[0] = tcc.geom->WireCoordinate(pos3[1], pos3[2], planeID);
        if(pos2[0] < -0.4) continue;
        // check for dead wires
        unsigned int wire = std::nearbyint(pos2[0]);
        if(wire > slc.nWires[plane]) continue;
        if(slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
        pos2[1] = tcc.detprop->ConvertXToTicks(pos3[0], planeID) * tcc.unitsPerTick;
        float cf = ChgFracNearPos(slc, pos2, tjids);
        if(cf < lo) lo = cf;
        if(cf > hi) hi = cf;
        sum += cf;
        ++cnt;
      } // tjid
    } // plane
    if(cnt == 0) return 0;
    if(cnt > 1 && lo < 0.3 && hi > 0.8) {
      sum -= lo;
      --cnt;
    }
    return sum / cnt;
  } // ChgFracNearEnd
  
  ////////////////////////////////////////////////
  Vector3_t DirAtEnd(const PFPStruct& pfp, unsigned short end)
  {
    if(end > 1 || pfp.SectionFits.empty()) return {{0., 0., 0.}};
    if(end == 0) return pfp.SectionFits[0].Dir;
    return pfp.SectionFits[pfp.SectionFits.size() - 1].Dir;
  } // PosAtEnd
  
  ////////////////////////////////////////////////
  Point3_t PosAtEnd(const PFPStruct& pfp, unsigned short end)
  {
    if(end > 1 || pfp.SectionFits.empty()) return {{0., 0., 0.}};
    // handle a neutrino pfp that doesn't have any TP3Ds
    if(pfp.TP3Ds.empty()) return pfp.SectionFits[0].Pos;
    if(end == 0) return pfp.TP3Ds[0].Pos;
    return pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
  } // PosAtEnd
  
  ////////////////////////////////////////////////
  float Length(const PFPStruct& pfp)
  {
    if(pfp.TP3Ds.empty()) return 0;
    return PosSep(pfp.TP3Ds[0].Pos, pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos);
  } // Length
  
  ////////////////////////////////////////////////
  bool FindStartEnd(const PFPStruct& pfp, unsigned short sfIndex, unsigned short& startPt, unsigned short& endPt)
  {
    // this assumes that the TP3Ds vector is sorted
    startPt = USHRT_MAX;
    endPt = USHRT_MAX;
    if(sfIndex >= pfp.SectionFits.size()) return false;
    
    bool first = true;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      if(tp3d.SFIndex < sfIndex) continue;
      if(first) {
        first = false;
        startPt = ipt;
      } // first
      if(tp3d.SFIndex > sfIndex) break;
      endPt = ipt;
    } // ipt
    return true;
    
  } // FindStartEnd
  
  ////////////////////////////////////////////////
  unsigned short FarEnd(TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos)
  {
    // Returns the end (0 or 1) of the pfp that is furthest away from the position pos
    if(pfp.ID == 0) return 0;
    if(pfp.TP3Ds.empty()) return 0;
    auto& pos0 = pfp.TP3Ds[0].Pos;
    auto& pos1 = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    if(PosSep2(pos1, pos) > PosSep2(pos0, pos)) return 1;
    return 0;
  } // FarEnd

  
  ////////////////////////////////////////////////
  void PrintTP3Ds(std::string someText, TCSlice& slc, const PFPStruct& pfp, short printPts)
  {
    if(pfp.TP3Ds.empty()) return;
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" pfp P"<<pfp.ID<<" NeedsUpdate? "<<pfp.NeedsUpdate<<"\n";
    if(!pfp.SectionFits.empty()) {
      myprt<<someText<<"  SFI ________Pos________   ________Dir_______ _____EndPos________ ChiDOF  NPts NeedsUpdate?\n";
      for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
        myprt<<someText<<std::setw(4)<<sfi;
        auto& sf = pfp.SectionFits[sfi];
        myprt<<std::fixed<<std::setprecision(1);
        unsigned short startPt = 0, endPt = 0;
        if(FindStartEnd(pfp, sfi, startPt, endPt)) {
          auto& start = pfp.TP3Ds[startPt].Pos;
          myprt<<std::setw(7)<<start[0]<<std::setw(7)<<start[1]<<std::setw(7)<<start[2];
        } else {
          myprt<<" Invalid";
        }
        myprt<<std::fixed<<std::setprecision(2);
        myprt<<std::setw(7)<<sf.Dir[0]<<std::setw(7)<<sf.Dir[1]<<std::setw(7)<<sf.Dir[2];
        myprt<<std::fixed<<std::setprecision(1);
        if(endPt < pfp.TP3Ds.size()) {
          auto& end = pfp.TP3Ds[endPt].Pos;
          myprt<<std::setw(7)<<end[0]<<std::setw(7)<<end[1]<<std::setw(7)<<end[2];
        } else {
          myprt<<" Invalid";
        }
        myprt<<std::setprecision(1)<<std::setw(6)<<sf.ChiDOF;
        myprt<<std::setw(6)<<sf.NPts;
        myprt<<std::setw(6)<<sf.NeedsUpdate;
        myprt<<"\n";
      } // sec
    } // SectionFits
    if(printPts < 0) {
      // print the head if we print all points
      myprt<<someText<<"  ipt SFI ________Pos________  Delta   Pull IsGood  along dE/dx  T_ipt_P:W:T    3D->2D\n";
    }
    unsigned short fromPt = 0;
    unsigned short toPt = pfp.TP3Ds.size() - 1;
    if(printPts >= 0) fromPt = toPt;
    for(unsigned short ipt = fromPt; ipt <= toPt; ++ipt) {
      auto tp3d = pfp.TP3Ds[ipt];
      myprt<<someText<<std::setw(4)<<ipt;
      myprt<<std::setw(4)<<tp3d.SFIndex;
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<std::setw(7)<<tp3d.Pos[0]<<std::setw(7)<<tp3d.Pos[1]<<std::setw(7)<<tp3d.Pos[2];
      myprt<<std::setprecision(1)<<std::setw(6)<<(tp3d.Pos[0] - tp3d.TPX);
      float pull = (tp3d.Pos[0] - tp3d.TPX) / sqrt(tp3d.TPXErr2);
      myprt<<std::setprecision(1)<<std::setw(6)<<pull;
      myprt<<std::setw(7)<<tp3d.IsGood;
      myprt<<std::setw(7)<<std::setprecision(1)<<tp3d.along;
      myprt<<std::setw(6)<<std::setprecision(2)<<dEdx(slc, tp3d);
      if(tp3d.slHitsIndex < slc.slHits.size()) {
        myprt<<" "<<PrintHit(slc.slHits[tp3d.slHitsIndex]);
      } else if(tp3d.TjID > 0) {
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        myprt<<" T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"_"<<PrintPos(slc, tp);
        // print the 2D position
        auto tmp = MakeBareTP(slc, tp3d.Pos, tp3d.Dir, tp.CTP);
        myprt<<" "<<PrintPos(slc, tmp.Pos);
      } else {
        myprt<<" UNDEFINED";
      }

      myprt<<"\n";
    } // ipt
    // print the end
  } // PrintTP3Ds
} // namespace
