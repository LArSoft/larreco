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
              auto& pos1 = p1.XYZ[e1];
              // require the end to be close to a TPC boundary
              if(InsideFV(slc1, p1, e1)) continue;
              auto& dir1 = p1.Dir[e1];
              for(unsigned short e2 = 0; e2 < 2; ++e2) {
                auto& pos2 = p2.XYZ[e2];
                // require the end to be close to a TPC boundary
                if(InsideFV(slc2, p2, e2)) continue;
                auto& dir2 = p2.Dir[e2];
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
          if(pfp.XYZ[end][2] < minZ) { minZ = pfp.XYZ[end][2]; minZIndx = slcIndex;  minZEnd = end; }
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
    // clear matchVec
    slc.matchVec.clear();
    // clear the kEnvInPFP bits on all Tjs. The bit will be set true when a TP is
    // used in a PFParticle
    for(auto& tj : slc.tjs) {
      for(auto& tp : tj.Pts) tp.Environment[kEnvInPFP] = false;
    } // tj
    
    // Match these points in 3D and put the results in slc.matchVec
    std::vector<MatchStruct> matVec;
    // Look for matches in all planes in the TPC
    MatchPlanes(slc, slc.nPlanes, matVec, prt);
    // Make 2-plane matches if we haven't hit the user-defined size limit
    if(matVec.size() < tcc.match3DCuts[4] && slc.nPlanes == 3) MatchPlanes(slc, 2, matVec, prt);
    
    if(matVec.empty()) return;
    std::cout<<"FPFP matVec size "<<matVec.size()<<"\n";
    
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
/*
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
        myprt<<pfp.XYZ[0][0]<<", "<<pfp.XYZ[0][1]<<", "<<pfp.XYZ[0][2];
        myprt<<") Dir "<<std::setprecision(2)<<std::setw(6)<<pfp.Dir[0][0]<<std::setw(6)<<pfp.Dir[0][1]<<std::setw(6)<<pfp.Dir[0][2];
        myprt<<" projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, pfp.XYZ[0], pfp.Dir[0], inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, pfp.TjIDs);
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
*/
    // create a PFParticle for each valid match combination
    for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
      prt = (tcc.dbgPFP && debug.MVI < matVec.size() && indx == debug.MVI);
      if(!prt) continue;
      std::cout<<"Detailed output mvi "<<indx<<"\n";
      auto& ms = matVec[indx];
      // ignore dead matches
      if(ms.Count == 0) continue;
      // skip this match if any of the trajectories is already matched or merged and killed
      bool skipit = false;
      // check for muons and delta rays or InShower Tjs
      bool has13 = false;
      bool has11 = false;
      for(unsigned short itj = 0; itj < ms.TjIDs.size(); ++itj) {
        if(ms.TjIDs[itj] <= 0 || ms.TjIDs[itj] > (int)slc.tjs.size()) {
          std::cout<<"FindPFParticles: bogus ms TjID "<<ms.TjIDs[itj]<<"\n";
          return;
        }
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
      // Set the PDGCode so DefinePFP can ignore incompatible matches
      pfp.PDGCode = pdgCode;
      // set the PDGCode to ensure that delta rays aren't merged with muons. PDGCodeVote
      // returns 0 if the vote is mixed
      if(has13) pfp.PDGCode = 13;
      if(has11) pfp.PDGCode = 11;
      // fill the TP3D points using the 2D trajectory points for Tjs in TjIDs
      FillTP3Ds(slc, pfp);
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
        myprt<<pfp.XYZ[0][0]<<", "<<pfp.XYZ[0][1]<<", "<<pfp.XYZ[0][2];
        myprt<<") Dir ("<<std::setprecision(2)<<pfp.Dir[0][0]<<", "<<pfp.Dir[0][1]<<", "<<pfp.Dir[0][2];
        myprt<<") projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, pfp.XYZ[0], pfp.Dir[0], inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, pfp.TjIDs);
        myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(slc, pfp.TjIDs, false);
        myprt<<" nTP3Ds "<<pfp.TP3Ds.size();
      } // prt
      if(!DefinePFP("FPFP", slc, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" DefinePFP failed";
        continue;
      }
      // Look for problems and fix them
      AddMissedTP3Ds(slc, pfp, prt);
      if(pfp.NeedsUpdate && !DefinePFP("FPFP", slc, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" DefinePFP failed after AddMissedTP3Ds";
        continue;
      }
    } // indx

  } // FindPFParticles
  
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
    constexpr float twopi = 2 * M_PI;
    constexpr float piOver2 = M_PI / 2;
    
    // create a temp vector to check for duplicates
    auto inMatVec = matVec;
    std::vector<MatchStruct> temp;
    
    // the minimum number of points for matching
    unsigned short minPts = 2;
    // override this with the user minimum for 2-plane matches
    if(numPlanes == 2) minPts = tcc.match3DCuts[2];
    
    // max number of match combos left
    unsigned int nAvailable = 0;
    if(matVec.size() < tcc.match3DCuts[4]) nAvailable = tcc.match3DCuts[4] - matVec.size();
    if(nAvailable == 0 || nAvailable > tcc.match3DCuts[4]) return;
    
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
            float dang = 0;
            if(tcc.match3DCuts[1] > 0) {
              dang = std::abs(DeltaAngle(tp3.Dir, iktp3.Dir));
              while(dang >  M_PI) dang -= twopi;
              if(dang >  piOver2) dang = M_PI - dang;
              float mcsmom = slc.tjs[iTjPt.id - 1].MCSMom + slc.tjs[jTjPt.id - 1].MCSMom + slc.tjs[kTjPt.id - 1].MCSMom;
              mcsmom /= 3;
              if(mcsmom > 150 && dang > tcc.match3DCuts[1]) continue;
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
    
    if(prt) mf::LogVerbatim("TC")<<"FindXMatches: Found "<<temp.size()<<" matches in "<<numPlanes<<" planes";
    
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
          // require that the Tj ID of this point be in the list if we aren't filling the Tp3s
//          if(!fillTp3s && kndx == tjids.size()) continue;
          if(kTjPt.xlo > iTjPt.xhi) continue;
          // break out if the x range difference becomes large
          if(kTjPt.xlo > iTjPt.xhi + 5) break;
          auto& ktj = slc.tjs[kTjPt.id - 1];
          //          if(ktj.AlgMod[kMat3D]) continue;
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
    Fit3D(2, point, dir, pfp.XYZ[0], pfp.Dir[0]);
    if(prt && doFit) {
      mf::LogVerbatim myprt("TC");
      myprt<<"FC: P"<<pfp.ID<<" fit pos "<<std::fixed<<std::setprecision(1)<<pfp.XYZ[0][0]<<" "<<pfp.XYZ[0][1]<<" "<<pfp.XYZ[0][2];
      myprt<<" fit dir "<<std::setprecision(2)<<pfp.Dir[0][0]<<" "<<pfp.Dir[0][1]<<" "<<pfp.Dir[0][2];
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
  bool DefinePFP(std::string inFcnLabel, TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function is called after the 3D matched TjIDs have been specified and optionally
    // a start or end vertex ID. It defines the PFParticle but doesn't store it
    
    if(pfp.PDGCode == 1111) return false;
    if(pfp.TjIDs.size() < 2) return false;
    
    std::string fcnLabel = inFcnLabel + ".DPFP";
    
    // require at least one tj in at least two planes
    std::vector<unsigned short> nInPln(slc.nPlanes);
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      ++nInPln[DecodeCTP(tj.CTP).Plane];
      if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) {
        std::cout<<fcnLabel<<" P"<<pfp.ID<<" uses T"<<tj.ID<<" but kMat3D is set true\n";
        return false;
      }
    } // tjid
    unsigned short npl = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(nInPln[plane] > 0) ++npl;
    if(npl < 2) return false;
    
    // Fit all points to get the general direction and central position and stuff
    // that into XYZ and Dir at the "start end"
    float chiDOF;
    if(!FitTP3Ds(slc, pfp.TP3Ds, 0, pfp.TP3Ds.size(), pfp.XYZ[0], pfp.Dir[0], chiDOF, true)) return false;
    // Correct the start position
    if(!SetStartEnd(slc, pfp)) return false;
    SortTP3Ds(slc, pfp);
    
    // Do a series of fits in sections with the requirement that there are at least
    // 5 2D points in at least 2 planes
    unsigned short min2DPts = 5;
    unsigned short fromPt = 0;
    Point3_t dumPos;
    Vector3_t dumDir;
    // Set ChiDOF negative to indicate a fit failure
    for(auto& tp3d : pfp.TP3Ds) tp3d.ChiDOF = -1;
    dumPos = pfp.TP3Ds[fromPt].Pos;
    for(auto& tp3d : pfp.TP3Ds) tp3d.flag = false;
    while(fromPt < pfp.TP3Ds.size()) {
      // iterate several times to reject outliers
      unsigned short toPt = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1);
      for(unsigned short nit = 0; nit < 3; ++nit) {
//        std::cout<<nit<<" chk "<<fromPt<<" - "<<toPt<<"\n";
        if(toPt > pfp.TP3Ds.size()) break;
        if(!FitTP3Ds(slc, pfp.TP3Ds, fromPt, toPt, dumPos, dumDir, chiDOF, true)) {
          std::cout<<fromPt<<"-"<<toPt<<" fit failed. Do something about it\n";
          break;
        }
        // no sense iterating if the fit is good
        if(pfp.TP3Ds[fromPt].ChiDOF < 2) break;
        // find the point with the largest Delta
        float maxDelta = 0;
        unsigned short badPt = fromPt;
        for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
          auto& tp3d = pfp.TP3Ds[ipt];
          if(!tp3d.IsGood) continue;
          if(std::abs(tp3d.Delta) < maxDelta) continue;
          maxDelta = std::abs(tp3d.Delta);
          badPt = ipt;
        } // ipt
        // clobber that point
        pfp.TP3Ds[badPt].IsGood = false;
        // re-find the toPt
        toPt = Find3DRecoRange(slc, pfp, fromPt, min2DPts, 1);
      } // nit
      if(toPt > pfp.TP3Ds.size()) break;
      // update all the points in each section
      fromPt = toPt;
      if(fromPt > pfp.TP3Ds.size()) break;
    } // ipt
    // The last section was most likely not fit because it wasn't reconstructable in 3D
    // so reverse the fromPt and toPt points and do the fit
    if(fromPt < pfp.TP3Ds.size()) {
      // Find the reconstructable range working backwards from the end
      fromPt = pfp.TP3Ds.size() - 1;
      unsigned short toPt = Find3DRecoRange(slc, pfp, fromPt, min2DPts, -1);
      // swap to and from
      std::swap(fromPt, toPt);
      // and increment to
      ++toPt;
      if(!FitTP3Ds(slc, pfp.TP3Ds, fromPt, toPt, dumPos, dumDir, chiDOF, true)) {
        std::cout<<fromPt<<"-"<<toPt<<" fit failed. Do something about it\n";
      }
    }
    // debug check that all the points were fit in each section
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.IsGood && !tp3d.flag) std::cout<<"point "<<ipt<<" not fit\n";
    }
    PrintTP3Ds("DP", slc, pfp, -1);
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<fcnLabel<<" pfp P"<<pfp.ID;
      myprt<<" ->";
      for(auto id : pfp.TjIDs) myprt<<" T"<<id;
      myprt<<" max Tj len "<<MaxTjLen(slc, pfp.TjIDs);
      myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
      myprt<<" Tp3s TP3Ds "<<pfp.TP3Ds.size();
    } // prt
    
    pfp.NeedsUpdate = false;
    return true;
  } // DefinePFP
  
  /////////////////////////////////////////
  unsigned short Find3DRecoRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short min2DPts, short dir)
  {
    // Scans the TP3Ds vector starting at fromPt until it finds min2DPts in two planes. It returns
    // with the index of that point (+1) in the TP3Ds vector. The dir variable defines the direction in
    // the TP3Ds vector
    if(fromPt > pfp.TP3Ds.size() - 1) return USHRT_MAX;
    if(pfp.TP3Ds.size() < 2 * min2DPts) return USHRT_MAX;
    if(dir == 0) return USHRT_MAX;
    
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      unsigned short ipt = fromPt + ii;
      if(dir > 0 && ipt == pfp.TP3Ds.size()) break;
      if(dir < 0) ipt = fromPt - ii;
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      unsigned short plane = 0;
      if(tp3d.TjID > 0) {
        plane = DecodeCTP(slc.tjs[tp3d.TjID - 1].CTP).Plane;
      } else {
        unsigned int ahi = slc.slHits[tp3d.slHitsIndex].allHitsIndex;
        auto& hit = (*evt.allHits)[ahi];
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
  bool FitTP3Ds(TCSlice& slc, PFPStruct& pfp, unsigned short atPt, unsigned short requestedNumPts)
  {
    // Fit the trajectory at point atPt using nTPsfit points and define the PFPstruct variables at
    // that point
    if(pfp.TP3Ds.size() < 3) return false;
    if(atPt > pfp.TP3Ds.size()) return false;
    if(requestedNumPts < 2) return false;
    
    unsigned short fromPt = 0;
    unsigned short toPt = 0;
    // ensure a valid value of fromPt
    short tmp = (short)atPt - (short)requestedNumPts / 2;
    if(tmp > 0) fromPt = (unsigned short)tmp;
    // find the upper bound (toPt - 1)
    unsigned short cnt = 0;
    for(toPt = fromPt; toPt < pfp.TP3Ds.size(); ++toPt) {
      if(!pfp.TP3Ds[toPt].IsGood) continue;
      ++cnt;
      if(cnt == requestedNumPts) break;
    } // toPt
    // check for failures
    if(toPt <= fromPt) return false;
    if(cnt < 3) return false;
    Point3_t pos;
    // define the fit origin
    pos[0] = pfp.TP3Ds[atPt].Pos[0];
    Vector3_t dir;
    float chiDOF;
    // Don't update all of the points between fromPt to toPt.
    if(FitTP3Ds(slc, pfp.TP3Ds, fromPt, toPt, pos, dir, chiDOF, false)) {
      // successfull fit. Update the point
      auto& tp3d = pfp.TP3Ds[atPt];
      tp3d.ChiDOF = chiDOF;
      tp3d.nTPsFit = cnt;
      tp3d.Pos = pos;
      tp3d.Dir = dir;
      Point2_t hitPos;
      CTP_t inCTP = 0;
      if(tp3d.TjID > 0) {
        // a real TP exists
        inCTP = slc.tjs[tp3d.TjID - 1].CTP;
        // transfer the hit position
        hitPos = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex].HitPos;
      } else {
        unsigned int ahi = slc.slHits[tp3d.slHitsIndex].allHitsIndex;
        auto& hit = (*evt.allHits)[ahi];
        inCTP = EncodeCTP(hit.WireID().Cryostat, hit.WireID().TPC, hit.WireID().Plane);
        // define the hit position
        hitPos[0] = hit.WireID().Wire;
        hitPos[1] = hit.PeakTime() * tcc.unitsPerTick;
      }
      // get the 2D position for this 3D point
      auto tp = MakeBareTP(slc, pos, dir, inCTP);
      tp3d.Delta = PosSep(hitPos, tp.Pos);
    } // success
    return true;
  } // FitTP3Ds

  /////////////////////////////////////////
  bool FitTP3Ds(TCSlice& slc, std::vector<TP3D>& tp3s, unsigned short fromPt, unsigned short toPt, Point3_t& pos, Vector3_t& dir, float& chiDOF, bool doUpdate)
  {
    // Fits the TP3D points in the selected range to a 3D line with X origin at pos[0]. Note that the
    // toPt is not included in the fit
    chiDOF = 999;
    if(tp3s.size() < 4) return false;
    if(fromPt >= toPt) return false;
    if(toPt > tp3s.size()) return false;
    
    double x0 = pos[0];
    
    unsigned short npts = 0;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) if(tp3s[ipt].IsGood) ++npts;
    // Need at least 4 + 1 points for a chisq calculation
    if(npts < 5) return false;
    
    const unsigned int nvars = 4;
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    
    // count the number of TPs in each plane
    std::vector<unsigned short> cntInPln(slc.nPlanes);
    std::vector<geo::PlaneID> plnIDs(slc.nPlanes);
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      plnIDs[plane] = geo::PlaneID(slc.TPCID.Cryostat, slc.TPCID.TPC, plane);
    } // plane
    double wght = 1;
    unsigned short cnt = 0;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = tp3s[ipt];
      if(!tp3d.IsGood) continue;
      double wire, x;
      unsigned short plane = 0;
      if(tp3d.slHitsIndex < slc.slHits.size()) {
        // single hit
        unsigned int ahi = slc.slHits[tp3d.slHitsIndex].allHitsIndex;
        auto& hit = (*evt.allHits)[ahi];
        plane = hit.WireID().Plane;
        wire = hit.WireID().Wire;
        x = tcc.detprop->ConvertTicksToX(hit.PeakTime(), plnIDs[plane]) - x0;
        double tickErr = hit.RMS() * tcc.hitErrFac * hit.Multiplicity();
        double xdx = tcc.detprop->ConvertTicksToX(hit.PeakTime() + tickErr, plnIDs[plane]) - x0;
        double xErr = xdx - x;
        // TODO: This isn't correct. We need to know the error on the wire
        // position - not the error on X
        wght = 1 / (xErr * xErr);
      } else {
        // TP
        auto& tp2 = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        plane = DecodeCTP(tp2.CTP).Plane;
        wire = tp2.Pos[0];
        x = tcc.detprop->ConvertTicksToX(tp2.Pos[1]/tcc.unitsPerTick, plnIDs[plane]) - x0;
        // The TP2 hit error is in WSE units and includes the wire error. See
        // StepUtils/DefineHitPos. Ignore these details for now and fake it.
        double tickErr = tp2.HitPosErr2 / (tcc.unitsPerTick * tcc.unitsPerTick);
        double xdx = tcc.detprop->ConvertTicksToX(tp2.Pos[1]/tcc.unitsPerTick + tickErr, plnIDs[plane]) - x0;
        double xErr = xdx - x;
        wght = 1 / (xErr * xErr);
      }
      ++cntInPln[plane];
      // get the wire plane offset
      double off = tcc.geom->WireCoordinate(0, 0, plnIDs[plane]);
      // get the "cosine-like" component
      double cw = tcc.geom->WireCoordinate(1, 0, plnIDs[plane]) - off;
      // the "sine-like" component
      double sw = tcc.geom->WireCoordinate(0, 1, plnIDs[plane]) - off;
      A[cnt][0] = wght * cw;
      A[cnt][1] = wght * sw;
      A[cnt][2] = wght * cw * x;
      A[cnt][3] = wght * sw * x;
      w[cnt] = wght * (wire - off);
      ++cnt;
    } // ipt
    // ensure there are at least two points in at least two planes
    unsigned short enufInPlane = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] > 1) ++enufInPlane;
    if(enufInPlane < 2) return false;
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    dir[0] = 1 / norm;
    dir[1] = tVec[2] / norm;
    dir[2] = tVec[3] / norm;
    pos[0] = x0;
    pos[1] = tVec[0];
    pos[2] = tVec[1];
    // calculate ChiDOF
    chiDOF = 0;
    // a temporary position to calculate ChiDOF
    Point3_t tmpPos;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = tp3s[ipt];
      if(!tp3d.IsGood) continue;
      unsigned short plane = 0;
      double wire;
      if(tp3d.slHitsIndex < slc.slHits.size()) {
        unsigned int ahi = slc.slHits[tp3d.slHitsIndex].allHitsIndex;
        auto& hit = (*evt.allHits)[ahi];
        plane = hit.WireID().Plane;
        wire = hit.WireID().Wire;
        tmpPos[0] = tcc.detprop->ConvertTicksToX(hit.PeakTime(), plnIDs[plane]);
      } else {
        auto& tp2 = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        plane = DecodeCTP(tp2.CTP).Plane;
        wire = tp2.Pos[0];
        tmpPos[0] = tcc.detprop->ConvertTicksToX(tp2.Pos[1]/tcc.unitsPerTick, plnIDs[plane]);
      }
      double off = tcc.geom->WireCoordinate(0, 0, plnIDs[plane]);
      double cw = tcc.geom->WireCoordinate(1, 0, plnIDs[plane]) - off;
      double sw = tcc.geom->WireCoordinate(0, 1, plnIDs[plane]) - off;
      double x = tmpPos[0] - x0;
      if(doUpdate) {
        // update the vector
        tp3d.Pos[0] = tmpPos[0];
        tp3d.Pos[1] = tVec[0] + tVec[2] * x;
        tp3d.Pos[2] = tVec[1] + tVec[3] * x;
        tp3d.Delta = tp3d.Pos[1] * cw + tp3d.Pos[2] * sw - (wire - off);
        tp3d.Dir = dir;
        tp3d.flag = true;
        chiDOF += tp3d.Delta * tp3d.Delta;
      } else {
        // don't update
        tmpPos[1] = tVec[0] + tVec[2] * x;
        tmpPos[2] = tVec[1] + tVec[3] * x;
        double delta = tmpPos[1] * cw + tmpPos[2] * sw - (wire - off);
        chiDOF += delta * delta;
      }
   } // ipt
    chiDOF /= (float)(npts - 4);
    // Update chidof and nTPsFit
    if(doUpdate) {
      for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
        auto& tp3d = tp3s[ipt];
        if(!tp3d.IsGood) continue;
        tp3d.ChiDOF = chiDOF;
        tp3d.nTPsFit = npts;
      } // ipt
    } // doUpdate
/*
    std::cout<<"FitTP3Ds: Pos "<<std::fixed<<std::setprecision(1)<<pos[0]<<" "<<pos[1]<<" "<<pos[2];
    std::cout<<" Dir "<<std::setprecision(2)<<dir[0]<<" "<<dir[1]<<" "<<dir[2];
    std::cout<<" npts "<<npts;
    std::cout<<" chiDOF "<<chiDOF<<"\n";
*/
    return true;
    
  } // FitTP3Ds

  /////////////////////////////////////////
  void AddMissedTP3Ds(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Analyze the TP3Ds vector to look for missing Tj points. Add them if any
    // are found and set pfp.NeedsUpdate true
    if(pfp.TP3Ds.empty()) return;
    
    // find the first and last occurrences of TPs on the 3D trajectory
    std::vector<std::pair<unsigned short, unsigned short>> firstLast(pfp.TjIDs.size());
    for(auto& fl : firstLast) fl.first = USHRT_MAX;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      unsigned short indx = 0;
      for(indx = 0; indx < pfp.TjIDs.size(); ++indx) if(pfp.TjIDs[indx] == tp3d.TjID) break;
      if(indx == pfp.TjIDs.size()) {
        std::cout<<"AddMissedTP3Ds: T"<<tp3d.TjID<<" isn't in TjIDs in P"<<pfp.ID<<"\n";
        continue;
      }
      if(firstLast[indx].first == USHRT_MAX) firstLast[indx].first = ipt;
      firstLast[indx].second = ipt;
    } // ipt

    // look for evidence of missing Tj points
    for(unsigned short indx = 0; indx < pfp.TjIDs.size(); ++indx) {
      auto& tj = slc.tjs[pfp.TjIDs[indx] - 1];
      if(firstLast[indx].first > 5) {
//        std::cout<<"Look in CTP "<<tj.CTP<<" start 0 -> "<<firstLast[indx].first<<"\n";
        AddMissedTP3Ds(slc, pfp, 0, firstLast[indx].first, tj.CTP, prt);
      }
      if(firstLast[indx].second < pfp.TP3Ds.size() - 5) {
//        std::cout<<"Look in CTP "<<tj.CTP<<" start "<<firstLast[indx].second<<" -> "<<pfp.TP3Ds.size()<<"\n";
        AddMissedTP3Ds(slc, pfp, firstLast[indx].second, USHRT_MAX, tj.CTP, prt);
      }
    } // indx
    
  } // AddMissedTP3Ds
  
  /////////////////////////////////////////
  void AddMissedTP3Ds(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, CTP_t inCTP, bool prt)
  {
    // Try to insert 2D trajectory points into the 3D trajectory point vector pfp.TP3Ds. This function appends
    // new Tp3Ds and sets the NeedsUpdate flag true. The calling function should call DefinePFP to finish
    if(fromPt > toPt) return;
    if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
    geo::PlaneID planeID = DecodeCTP(inCTP);
    
    constexpr float maxPull = 3; 
    
    // reset the general purpose flag on all Tjs in this CTP so we can use it
    for(auto& tj : slc.tjs) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        tj.Pts[ipt].Environment[kEnvFlag] = false;
      } // ipt
    } // tj
    
    // put the new points into a temporary vector so we can decide if adding them
    // makes sense
    std::vector<TP3D> temp;
    
    bool added = false;
    
    TrajPoint tp2;
    float maxDelta = 2;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      tp2 = MakeBareTP(slc, tp3d.Pos, tp3d.Dir, inCTP);
      // check for a failure
      if(tp2.Pos[0] < 0) continue;
      tp2.Hits.clear();
      if(FindCloseHits(slc, tp2, maxDelta, kAllHits)) {
        for(auto iht : tp2.Hits) {
          // ignore hits that are used in a PFP
          if(slc.slHits[iht].InPFP > 0) continue;
          int inTj = slc.slHits[iht].InTraj;
          // create a new TP3D
          TP3D tp3d;
          if(inTj <= 0) {
            // using a single unused hit
            tp3d.slHitsIndex = iht;
            unsigned int ahi = slc.slHits[iht].allHitsIndex;
            auto& hit = (*evt.allHits)[ahi];
            tp3d.Pos[0] = tcc.detprop->ConvertTicksToX(hit.PeakTime(), planeID);
            temp.push_back(tp3d);
          } else {
            // using a TP
            tp3d.TjID = inTj;
            // find the TP index
            unsigned short tjPt = USHRT_MAX;
            auto& tj = slc.tjs[inTj - 1];
            for(tjPt = tj.EndPt[0]; tjPt <= tj.EndPt[1]; ++tjPt) {
              if(std::find(tj.Pts[tjPt].Hits.begin(), tj.Pts[tjPt].Hits.end(), iht) != tj.Pts[tjPt].Hits.end()) break;
            } // tjPt
            if(tjPt == USHRT_MAX) continue;
            auto& tpm = tj.Pts[tjPt];
            if(tpm.Environment[kEnvFlag]) continue;
            if(tpm.Environment[kEnvInPFP]) {
              std::cout<<"FindPFParticles: coding error. Trying to make a TP -> PFP assn when one already exists\n";
              return;
            } // kEnvInPFP
            // see if the TP is already used
            if(tpm.Environment[kEnvFlag]) continue;
            tp3d.TPIndex = tjPt;
            tp3d.Pos[0] = tcc.detprop->ConvertTicksToX(tpm.Pos[1]/tcc.unitsPerTick, planeID);
            tpm.Environment[kEnvFlag] = true;
            temp.push_back(tp3d);
          } // inTj > 0
        } // iht
      } // FindCloseHits
    } // ipt
    if(temp.empty()) return;
    
    // make a list of Tjs that are not in pfp.TjIDs and have a significant fraction of points used in this PFP
    std::vector<int> tjList;
    for(auto& tp3d : temp) {
      // ignore single hits for this check
      if(tp3d.TjID <= 0) continue;
      if(std::find(pfp.TjIDs.begin(), pfp.TjIDs.end(), tp3d.TjID) != pfp.TjIDs.end()) continue;
      if(std::find(tjList.begin(), tjList.end(), tp3d.TjID) != tjList.end()) continue;
      // determine the fraction of points used in this PFP
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
      myprt<<"AddMissedTP3Ds: P"<<pfp.ID<<" Found "<<temp.size()<<" TP3D candidates inCTP "<<inCTP;
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
        myprt<<" T"<<tid<<" cnt "<<(int)cnt<<" fracUsed "<<use/cnt<<"\n";
      } // tid
    } // prt
    // Add the TP3Ds that are in the list and the single-hit TP3Ds after making
    // a rough check on the pulls
    for(auto& tp3d : temp) {
      if(tp3d.TjID > 0 && tp3d.slHitsIndex < slc.slHits.size()) {
        std::cout<<"Invalid combination\n";
        exit(1);
      }
      // We need to find the tj that is in the same CTP as the tp3d
      int useTj = 0;
      for(auto tjid : pfp.TjIDs) {
        if(slc.tjs[tjid - 1].CTP != slc.tjs[tp3d.TjID - 1].CTP) continue;
        useTj = tjid;
        break;
      } // tjid
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
        if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" candidate hit "<<tp3d.slHitsIndex<<" "<<PrintHit(slc.slHits[tp3d.slHitsIndex])<<" pull "<<pull<<" maxPull "<<maxPull;
        if(pull > maxPull) continue;
        slc.slHits[tp3d.slHitsIndex].InPFP = pfp.ID;
        pfp.TP3Ds.push_back(tp3d);
        added = true;
      } else {
        // don't add the ones that aren't in the list
        if(std::find(tjList.begin(), tjList.end(), tp3d.TjID) == tjList.end()) continue;
        // TP
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        // check the pull
        float pull = -1;
        if(useTj > 0) {
          auto& tj = slc.tjs[useTj - 1];
          pull = PointPull(slc, tp.Pos, tp.Chg, tj);
        } // useTj > 0
        if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" candidate TP T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<" "<<PrintPos(slc, tp)<<" pull "<<pull<<" maxPull "<<maxPull;
        if(pull > maxPull) continue;
        // flag all of the used hits in this TP
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tp.UseHit[ii]) slc.slHits[tp.Hits[ii]].InPFP = pfp.ID;
        } // ii
        pfp.TP3Ds.push_back(tp3d);
        added = true;
      }
    } // tp3d
    if(prt) mf::LogVerbatim("TC")<<" P"<<pfp.ID<<" added TP3Ds? "<<added;
    if(!pfp.NeedsUpdate && added) pfp.NeedsUpdate = true;
  } // AddMissedTP3Ds
  
  /////////////////////////////////////////
  bool SetStartEnd(TCSlice& slc, PFPStruct& pfp)
  {
    // Find the minimum value of alongTrans[0] and set that as the
    // PFP start position
    
    // first ensure that the direction is defined
    bool itsok = false;
    for(unsigned short xyz = 0; xyz < 3; ++xyz) if(pfp.Dir[0][xyz] != 0) itsok = true;
    if(!itsok) return false;
    float minAlong = 1E6;
    unsigned short minAlongPt = 0;
    float maxAlong = -1E6;
    unsigned short maxAlongPt = 0;
    Point2_t alongTrans;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      FindAlongTrans(pfp.XYZ[0], pfp.Dir[0], tp3d.Pos, alongTrans);
      if(alongTrans[0] < minAlong) {
        minAlong = alongTrans[0];
        minAlongPt = ipt;
      }
      if(alongTrans[0] > maxAlong) {
        maxAlong = alongTrans[0];
        maxAlongPt = ipt;
      }
    } // tp3d
    pfp.XYZ[0] = pfp.TP3Ds[minAlongPt].Pos;
    pfp.Dir[0] = pfp.TP3Ds[minAlongPt].Dir;
    pfp.XYZ[1] = pfp.TP3Ds[maxAlongPt].Pos;
    pfp.Dir[1] = pfp.TP3Ds[maxAlongPt].Dir;
    // update the along variable
    for(auto& tp3d : pfp.TP3Ds) tp3d.along -= minAlong;
    return true;
  } // SetStartEnd

  /////////////////////////////////////////
  void SortTP3Ds(TCSlice& slc, PFPStruct& pfp)
  {
    // sorts the TP3Ds by the distance from the start. Note that the !IsGood TP3s are
    // sorted as well.
    Point2_t alongTrans;
    std::vector<SortEntry> sortVec(pfp.TP3Ds.size());
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      FindAlongTrans(pfp.XYZ[0], pfp.Dir[0], tp3d.Pos, alongTrans);
      sortVec[ipt].index = ipt;
      sortVec[ipt].val = alongTrans[0];
      tp3d.along = alongTrans[0];
    } // tp3d
    std::sort(sortVec.begin(), sortVec.end(), valIncreasings);
    std::vector<TP3D> temp(pfp.TP3Ds.size());
    for(unsigned int ii = 0; ii < sortVec.size(); ++ii) temp[ii] = pfp.TP3Ds[sortVec[ii].index];
    pfp.TP3Ds = temp;
    // set the start and end points
    pfp.XYZ[0] = pfp.TP3Ds[0].Pos;
    pfp.Dir[0] = pfp.TP3Ds[0].Dir;
    pfp.XYZ[1] = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Pos;
    pfp.Dir[1] = pfp.TP3Ds[pfp.TP3Ds.size() - 1].Dir;
  } // SortTP3Ds
  
  /////////////////////////////////////////
  void FillTP3Ds(TCSlice& slc, PFPStruct& pfp)
  {
    // Fill the TP3Ds vector and flag hits
    pfp.TP3Ds.clear();
    for(auto tid : pfp.TjIDs) {
      auto& tj = slc.tjs[tid - 1];
      TP3D tp3d;
      tp3d.TjID = tid;
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        geo::PlaneID planeID = DecodeCTP(tp.CTP);
        // Define the X position. Y and Z will be determined later in a fit
        tp3d.Pos[0] = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, planeID);
        tp3d.TPIndex = ipt;
        pfp.TP3Ds.push_back(tp3d);
/* This should be done later in StorePFP
        // flag the hits used
        for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
          if(tp.UseHit[ii]) slc.slHits[tp.Hits[ii]].InPFP = pfp.ID;
        } // ii
*/
      } // ipt
    } // tid
  } // FillTP3Ds
  
  /////////////////////////////////////////
  void UpdateMatchStructs(TCSlice& slc, int oldTj, int newTj)
  {
    // Replaces tjid and ipt references in slc.matchVec and slc.pfps from
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
    
    slc.matchVec.clear();
    
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
    
    if(tcc.match3DCuts.size() > 7 && tcc.match3DCuts[7] > 0 && ntp > tcc.match3DCuts[7]) ntp = tcc.match3DCuts[7];
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
    
    // This is pretty self-explanatory
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
/*
  /////////////////////////////////////////
  void FilldEdx(TCSlice& slc, PFPStruct& pfp)
  {
    // Fills the dEdX vector in the match struct. This function should be called after the
    // matched trajectory points are ordered so that dE/dx is calculated at the start of the PFParticle
    if(pfp.ID == 0) return;
    // error check
    bool notgood = false;
    for(unsigned short end = 0; end < 2; ++end) {
      if(pfp.dEdx[end].size() != slc.nPlanes) notgood = true;
      if(pfp.dEdxErr[end].size() != slc.nPlanes) notgood = true;
    }
    if(notgood) {
      //      if(prt) mf::LogVerbatim("TC")<<"FilldEdx found inconsistent sizes for dEdx\n";
      return;
    }
    
    double t0 = 0;
    
    // don't attempt to find dE/dx at the end of a shower
    unsigned short numEnds = 2;
    if(pfp.PDGCode == 1111) numEnds = 1;
    
    unsigned short maxlen = 0;
    for(auto tjID : pfp.TjIDs) {
      
      Trajectory& tj = slc.tjs[tjID - 1];
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      double angleToVert = tcc.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
      for(unsigned short end = 0; end < numEnds; ++end) {
        pfp.dEdx[end][planeID.Plane] = 0;
        tj.dEdx[end] = 0;
        double cosgamma = std::abs(std::sin(angleToVert) * pfp.Dir[end][1] + std::cos(angleToVert) * pfp.Dir[end][2]);
        if(cosgamma == 0) continue;
        double dx = tcc.geom->WirePitch(planeID) / cosgamma;
        if(dx == 0) continue;
        double dQ = tj.Pts[tj.EndPt[end]].AveChg;
        if(dQ == 0) continue;
        // convert to dQ/dx
        dQ /= dx;
        double time = tj.Pts[tj.EndPt[end]].Pos[1] / tcc.unitsPerTick;
        float dedx = tcc.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
        if(dedx > 999) dedx = -1;
        pfp.dEdx[end][planeID.Plane] = dedx;
        tj.dEdx[end] = dedx;
        // ChgRMS is the fractional error
        pfp.dEdxErr[end][planeID.Plane] = dedx * tj.ChgRMS;
        
      } // end
      // Grab the best plane iusing the start f 1 < dE/dx < 50 MeV/cm
      if(pfp.dEdx[0][planeID.Plane] > 1 && pfp.dEdx[0][planeID.Plane] < 50) {
        if(tj.Pts.size() > maxlen) {
          maxlen = tj.Pts.size();
          pfp.BestPlane = planeID.Plane;
        }
      } // valid dE/dx
      
    } // tj
  } // FilldEdX

  ////////////////////////////////////////////////
  void FilldEdx(TCSlice& slc, TrajPoint3& tp3)
  {
    // fills the dEdx variables in tp3 (MeV/cm)
    tp3.dEdx = 0;
    if(tp3.Tj2Pts.empty()) return;
    double t0 = 0;
    float sum = 0;
    float sum2 = 0;
    float cnt = 0;
    for(auto& tj2pt : tp3.Tj2Pts) {
      auto& tp = slc.tjs[tj2pt.id - 1].Pts[tj2pt.ipt];
      if(tp.Chg <= 0) continue;
      geo::PlaneID planeID = DecodeCTP(tp.CTP);
      double angleToVert = tcc.geom->Plane(planeID).ThetaZ() - 0.5 * ::util::pi<>();
      double cosgamma = std::abs(std::sin(angleToVert) * tp3.Dir[1] + std::cos(angleToVert) * tp3.Dir[2]);
      if(cosgamma == 0) continue;
      double dx = tcc.geom->WirePitch(planeID) / cosgamma;
      // sum the hit charge
      double chg = 0;
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned int iht = tp.Hits[ii];
        chg += (*evt.allHits)[slc.slHits[iht].allHitsIndex].Integral();
      }
      double dQ = chg / dx;
      double time = tp.Pos[1] / tcc.unitsPerTick;
      float dedx = tcc.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
      if(dedx > 200) continue;
      sum += dedx;
      sum2 += dedx * dedx;
      ++cnt;
    } // tj2pt
    if(cnt == 0) return;
    tp3.dEdx = sum / cnt;
  } // FilldEdx
*/
/*
  ////////////////////////////////////////////////
  float PFPDOCA(const PFPStruct& pfp1,  const PFPStruct& pfp2, unsigned short& close1, unsigned short& close2)
  {
    // returns the Distance of Closest Approach between two PFParticles. 
    close1 = USHRT_MAX;
    close2 = USHRT_MAX;
    float minSep2 = 1E8;
    for(unsigned short ipt1 = 0; ipt1 < pfp1.Tp3s.size(); ++ipt1) {
      auto& tp1 = pfp1.Tp3s[ipt1];
      for(unsigned short ipt2 = 0; ipt2 < pfp2.Tp3s.size(); ++ipt2) {
        auto& tp2 = pfp2.Tp3s[ipt2];
        float sep2 = PosSep2(tp1.Pos, tp2.Pos);
        if(sep2 > minSep2) continue;
        minSep2 = sep2;
        close1 = ipt1;
        close2 = ipt2;
      } // tp2
    } // tp1
    return sqrt(minSep2);
 } // PFPDOCA
*/

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
      pfp.dEdx[0].resize(slc.nPlanes, 0);
      pfp.dEdx[1].resize(slc.nPlanes, 0);
      pfp.dEdxErr[0].resize(slc.nPlanes, 0);
      pfp.dEdxErr[1].resize(slc.nPlanes, 0);
    }
    for(unsigned short end = 0; end < 2; ++end) {
      // BUG the double brace syntax is required to work around clang bug 21629
      // (https://bugs.llvm.org/show_bug.cgi?id=21629)
      pfp.Dir[end] = {{0.0, 0.0, 0.0}};
      pfp.DirErr[end] = {{0.0, 0.0, 0.0}};
      pfp.XYZ[end] = {{0.0, 0.0, 0.0}};
    }
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
      Vtx3Store vx3;
      vx3.TPCID = pfp.TPCID;
      vx3.Vx2ID.resize(slc.nPlanes);
      // Flag it as a PFP vertex that isn't required to have matched 2D vertices
      vx3.Wire = -2;
      vx3.X = pfp.XYZ[0][0];
      vx3.Y = pfp.XYZ[0][1];
      vx3.Z = pfp.XYZ[0][2];
      vx3.ID = slc.vtx3s.size() + 1;
      vx3.Primary = false;
      ++evt.globalS3ID;
      vx3.UID = evt.globalS3ID;
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
    // stores the PFParticle in TJStuff
    if(pfp.ID < int(slc.pfps.size())) return false;
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    if(!neutrinoPFP) {
      if(pfp.TjIDs.empty()) return false;
      if(pfp.PDGCode != 1111 && pfp.TP3Ds.size() < 2) return false;
    }
    // check the ID and correct it if it is wrong
    if(pfp.ID != (int)slc.pfps.size() + 1) pfp.ID = slc.pfps.size() + 1;
    ++evt.globalPFPID;
    pfp.UID = evt.globalPFPID;
    
    // set the TP used in PFP bit
    for(auto& tp3d : pfp.TP3Ds) {
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvInPFP]) {
        std::cout<<"StorePFP: Trying to use already-used TP T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"\n";
        return false;
      } // error
      tp.Environment[kEnvInPFP] = true;
    } // tp3d
    
    // check the Tjs and set the 3D match flag
    for(auto tjid : pfp.TjIDs) {
      auto& tj = slc.tjs[tjid - 1];
      if(tj.AlgMod[kMat3D]) return false;
      tj.AlgMod[kMat3D] = true;
    } // tjid
    
//    if(pfp.BestPlane < 0) FilldEdx(slc, pfp);
    
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
    
    float abit = 5;
    auto& pos1 = pfp.XYZ[end];
    return (pos1[0] > slc.xLo + abit && pos1[0] < slc.xHi - abit && 
            pos1[1] > slc.yLo + abit && pos1[1] < slc.yHi - abit &&
            pos1[2] > slc.zLo + abit && pos1[2] < slc.zHi - abit);
    
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
/*
  ////////////////////////////////////////////////
  void ReversePFP(TCSlice& slc, PFPStruct& pfp)
  {
    std::swap(pfp.XYZ[0], pfp.XYZ[1]);
    std::swap(pfp.Dir[0], pfp.Dir[1]);
    for(unsigned short xyz = 0; xyz < 3; ++xyz) {
      pfp.Dir[0][xyz] *= -1;
      pfp.Dir[1][xyz] *= -1;
    }
    std::swap(pfp.DirErr[0], pfp.DirErr[1]);
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    std::reverse(pfp.Tp3s.begin(), pfp.Tp3s.end());
    for(auto& tp3 : pfp.Tp3s) {
      for(unsigned short xyz = 0; xyz < 3; ++xyz) {
        tp3.Dir[xyz] *= -1;
        tp3.Dir[xyz] *= -1;
      } // xyz
    } // tp3
  } // ReversePFP
*/
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

    float sum = 0;
    float cnt = 0;
    // keep track of the lowest value and maybe reject it
    float lo = 1;
    float hi = 0;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      std::vector<int> tjids(1);
      for(auto tjid : pfp.TjIDs) {
        auto& tj = slc.tjs[tjid - 1];
        if(tj.CTP != inCTP) continue;
        tjids[0] = tjid;
        Point2_t pos;
        geo::PlaneID planeID = geo::PlaneID(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        pos[0] = tcc.geom->WireCoordinate(pfp.XYZ[end][1], pfp.XYZ[end][2], planeID);
        if(pos[0] < -0.4) continue;
        // check for dead wires
        unsigned int wire = std::nearbyint(pos[0]);
        if(wire > slc.nWires[plane]) continue;
        if(slc.wireHitRange[plane][wire].first == UINT_MAX) continue;
        pos[1] = tcc.detprop->ConvertXToTicks(pfp.XYZ[end][0], planeID) * tcc.unitsPerTick;
        float cf = ChgFracNearPos(slc, pos, tjids);
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
  unsigned short FarEnd(TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos)
  {
    // Returns the end (0 or 1) of the pfp that is furthest away from the position pos
    if(pfp.ID == 0) return 0;
    if(PosSep2(pfp.XYZ[1], pos) > PosSep2(pfp.XYZ[0], pos)) return 1;
    return 0;
  } // FarEnd

  
  ////////////////////////////////////////////////
  void PrintTP3Ds(std::string someText, TCSlice& slc, const PFPStruct& pfp, short printPts)
  {
    if(pfp.TP3Ds.empty()) return;
    mf::LogVerbatim myprt("TC");
    if(printPts < 0) {
      // print the head if we are print all points
      myprt<<someText<<" pfp P"<<pfp.ID<<"\n";
      myprt<<someText<<"  ipt ________Pos________  Path   ________Dir_______ Delta IsGood ChiDOF nTPsFit  T_ipt_P:W:T\n";
    }
    // print the start
    myprt<<someText<<"    ";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(7)<<pfp.XYZ[0][0]<<std::setw(7)<<pfp.XYZ[0][1]<<std::setw(7)<<pfp.XYZ[0][2];
    myprt<<"      ";
    myprt<<std::fixed<<std::setprecision(2);
    myprt<<std::setw(7)<<pfp.Dir[0][0]<<std::setw(7)<<pfp.Dir[0][1]<<std::setw(7)<<pfp.Dir[0][2];
    myprt<<" <--- pfp.XYZ[0] \n";
    
    unsigned short fromPt = 0;
    unsigned short toPt = pfp.TP3Ds.size() - 1;
    if(printPts >= 0) fromPt = toPt;
    //    Vector3_t prevDir = pfp.Dir[0];
    for(unsigned short ipt = fromPt; ipt <= toPt; ++ipt) {
      auto tp3d = pfp.TP3Ds[ipt];
      myprt<<someText<<std::setw(4)<<ipt;
      myprt<<std::fixed<<std::setprecision(1);
      myprt<<std::setw(7)<<tp3d.Pos[0]<<std::setw(7)<<tp3d.Pos[1]<<std::setw(7)<<tp3d.Pos[2];
      myprt<<std::setprecision(1)<<std::setw(6)<<PosSep(tp3d.Pos, pfp.XYZ[0]);
      myprt<<std::setprecision(2)<<std::setw(7)<<tp3d.Dir[0]<<std::setw(7)<<tp3d.Dir[1]<<std::setw(7)<<tp3d.Dir[2];
      myprt<<std::setprecision(1)<<std::setw(6)<<tp3d.Delta;
      myprt<<std::setw(7)<<tp3d.IsGood;
      myprt<<std::setprecision(1)<<std::setw(6)<<tp3d.ChiDOF;
      myprt<<std::setw(6)<<tp3d.nTPsFit;
      if(tp3d.slHitsIndex < slc.slHits.size()) {
        myprt<<" "<<PrintHit(slc.slHits[tp3d.slHitsIndex]);
      } else {
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        myprt<<" T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"_"<<PrintPos(slc, tp);
      }
      myprt<<" along "<<std::setprecision(3)<<tp3d.along;
      myprt<<"\n";
    } // ipt
    // print the end
    myprt<<someText<<"    ";
    myprt<<std::fixed<<std::setprecision(1);
    myprt<<std::setw(7)<<pfp.XYZ[1][0]<<std::setw(7)<<pfp.XYZ[1][1]<<std::setw(7)<<pfp.XYZ[1][2];
    myprt<<"      ";
    myprt<<std::fixed<<std::setprecision(2);
    myprt<<std::setw(7)<<pfp.Dir[1][0]<<std::setw(7)<<pfp.Dir[1][1]<<std::setw(7)<<pfp.Dir[1][2];
    myprt<<" <--- pfp.XYZ[1]. Length "<<PosSep(pfp.XYZ[0], pfp.XYZ[1])<<"\n";
    
  } // PrintTP3Ds
} // namespace
