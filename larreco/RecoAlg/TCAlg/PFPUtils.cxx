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
<<<<<<< HEAD
  void FindPFParticles(TCSlice& slc)
  {
    // Match Tjs in 3D and create PFParticles
    
    if(tcc.match3DCuts[0] <= 0) return;
    
    FillmAllTraj(slc);
    
    bool prt = (tcc.dbgPFP && tcc.dbgSlc);

    // clear the TP -> P assn Tjs
    for(auto& tj : slc.tjs) {
      for(auto& tp : tj.Pts) tp.InPFP = 0;
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
      std::vector<SortEntry> sortVec;
      for(unsigned int ii = 0; ii < matVec.size(); ++ii) {
        auto& ms = matVec[ii];
        // count the number of TPs in all Tjs
        float tpCnt = 0;
        for(auto tid : ms.TjIDs) {
          auto& tj = slc.tjs[tid - 1];
          tpCnt += NumPtsWithCharge(slc, tj, false);
        } // tid
        float frac = ms.Count / tpCnt;
        frac /= ms.TjIDs.size();
        // ignore matches with a very low match fraction
        if(frac < 0.05) continue;
        SortEntry se;
        se.index = ii;
        se.val = ms.Count;
        // de-weight two-plane matches
        if(ms.TjIDs.size() == 2) se.val /= 2;
        sortVec.push_back(se);
      } // ii
      std::sort(sortVec.begin(), sortVec.end(), valDecreasings);
      std::vector<MatchStruct> tmpVec(sortVec.size());
      for(unsigned int ii = 0; ii < sortVec.size(); ++ii) {
        tmpVec[ii] = matVec[sortVec[ii].index];
      } // ii
      matVec = tmpVec;
    } // sort matVec
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"MVI  Count  Tjs\n";
      for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
        auto& ms = matVec[indx];
        myprt<<std::setw(5)<<indx<<std::setw(6)<<(int)ms.Count;
        for(auto tid : ms.TjIDs) myprt<<" T"<<tid;
        myprt<<"\n";
      } // indx
    } // prt

    // create a PFParticle for each valid match combination
    for(unsigned int indx = 0; indx < matVec.size(); ++indx) {
      // tone down the level of printing in ReSection
      bool foundMVI = (tcc.dbgPFP && indx == debug.MVI);
      auto& ms = matVec[indx];
      if(foundMVI) {
        std::cout<<"found MVI "<<indx<<"\n";
      }
      // ignore dead matches
      if(ms.Count == 0) continue;
      // skip this match if any of the trajectories is already matched or merged or killed
      bool skipit = false;
      for(unsigned short itj = 0; itj < ms.TjIDs.size(); ++itj) {
        auto& tj = slc.tjs[ms.TjIDs[itj] - 1];
        if(tj.AlgMod[kMat3D] || tj.AlgMod[kKilled]) skipit = true;
      } // tjID
      if(skipit) continue;
      int pdgCode = PDGCodeVote(slc, ms.TjIDs, prt);
      PFPStruct pfp = CreatePFP(slc);
      pfp.TjIDs = ms.TjIDs;
      pfp.MVI = indx;
      // fill the TP3D points using the 2D trajectory points for Tjs in TjIDs. All
      // points are put in one section
      MakeTP3Ds(slc, pfp);
      // fit all the points to get the general direction
      if(!FitSection(slc, pfp, 0)) {
        std::cout<<"FPFP: First fit failed\n";
        continue;
      }
      // a temp function for determining the chisq cut
      //      ChkPFPMC(slc, pfp);
      if(pfp.SectionFits[0].ChiDOF > 100) continue;
      // sort the points by the distance along the general direction vector
      if(!SortSection(pfp, 0)) continue;
      // Skip this combination if it isn't reconstructable in 3D
      if(Find3DRecoRange(slc, pfp, 0, (unsigned short)tcc.match3DCuts[3], 1) == USHRT_MAX) continue;
      // Do a fit in multiple sections if the initial fit is poor
      auto& sf0 = pfp.SectionFits[0];
      if(foundMVI) {
        std::cout<<"stop here before ReSection\n";
      }
      if(foundMVI) {
        PrintTP3Ds("FF", slc, pfp, -1);
      }
      if(sf0.ChiDOF < 2) {
        // Good fit with one section
        pfp.NeedsUpdate = false;
      } else if(!ReSection(slc, pfp, foundMVI)) {
        std::cout<<"ReSection failed. This is bad...\n";
        continue;
      }
      pfp.PDGCode = pdgCode;
      if(prt) {
        mf::LogVerbatim myprt("TC");
        myprt<<std::setw(4)<<"MVI "<<indx<<" Count "<<std::setw(5)<<(int)ms.Count;
        myprt<<" P"<<pfp.ID;
        myprt<<" ->";
        for(auto& tjid : pfp.TjIDs) myprt<<" T"<<tjid;
        myprt<<" projInPlane";
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, pfp.SectionFits[0].Pos, pfp.SectionFits[0].Dir, inCTP);
          myprt<<" "<<std::setprecision(2)<<tp.Delta;
        } // plane
        myprt<<" maxTjLen "<<(int)MaxTjLen(slc, pfp.TjIDs);
        myprt<<" MCSMom "<<MCSMom(slc, pfp.TjIDs);
        myprt<<" PDGCodeVote "<<PDGCodeVote(slc, pfp.TjIDs, false);
        myprt<<" nTP3Ds "<<pfp.TP3Ds.size();
        myprt<<" Reco3DRange "<<Find3DRecoRange(slc, pfp, 0, (unsigned short)tcc.match3DCuts[3], 1);
      } // prt
      // Trim points from the ends until there is a 3D point where there is a signal in at least two planes
      TrimEndPts(slc, pfp, prt);
      // FillGaps3D looks for gaps in the TP3Ds vector caused by broken trajectories and
      // inserts new TP3Ds if there are hits in the gaps. This search is only done in a
      // plane if the projection of the pfp results in a large angle where 2D reconstruction
      // is likely to be poor
      // This needs work...
      FillGaps3D(slc, pfp, prt);
      // Look for mis-placed 2D and 3D vertices
      ReconcileVertices(slc, pfp, prt);
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
        std::cout<<"P"<<pfp.ID<<" has Tjs with more than one Bragg peak. Write some code MVI = "<<pfp.MVI<<"\n";
      }
      // debug?
      FilldEdx(slc, pfp);
      if(tcc.dbgPFP && pfp.MVI == debug.MVI) PrintTP3Ds("FPFP", slc, pfp, -1);
      if(!StorePFP(slc, pfp)) {
        std::cout<<"StorePFP failed P"<<pfp.ID<<"\n";
        continue;
      }
    } // indx
    slc.mallTraj.resize(0);
    
    ReconcileTPs(slc);

  } // FindPFParticles
/*
  ////////////////////////////////////////////////
  void ChkPFPMC(TCSlice& slc, PFPStruct& pfp)
  {
    // This function is used to decide what ChiDOF cut should be made to reject
    // invalid 3D matches 
    if(evt.allHitsMCPIndex.empty()) return;
    if(pfp.SectionFits.size() != 1 || pfp.TP3Ds.empty()) {
      std::cout<<"ChkPFPMC: Something wrong with P"<<pfp.ID<<"\n";
      return;
    }
    
    // mcpIndex and count
    std::vector<std::pair<unsigned int, unsigned short>> mcpi_cnt;
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.TjID <= 0) {
        std::cout<<"oops\n";
        exit(1);
      }
      unsigned int mcpIndex = UINT_MAX;
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        unsigned ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
        mcpIndex = evt.allHitsMCPIndex[ahi];
        break;
      } // ii
      if(mcpIndex == UINT_MAX) continue;
      // look for it in the list
      unsigned short indx = 0;
      for(indx = 0; indx < mcpi_cnt.size(); ++indx) if(mcpi_cnt[indx].first == mcpIndex) break;
      // not found so add it
      if(indx == mcpi_cnt.size()) mcpi_cnt.push_back(std::make_pair(mcpIndex, 0));
      ++mcpi_cnt[indx].second;
    } // tp3d
    auto& sf = pfp.SectionFits[0];
    std::cout<<"ChkPFPMC: P"<<pfp.ID;
    std::cout<<std::setprecision(1)<<" ChiDOF "<<sf.ChiDOF;
    std::cout<<" MCP_cnt";
    for(auto mc : mcpi_cnt) std::cout<<" "<<mc.first<<"_"<<mc.second;
    std::cout<<"\n";

  } // ChkPFPMC
*/
  ////////////////////////////////////////////////
  void ReconcileTPs(TCSlice& slc)
  {
    // Reconciles TP ownership conflicts between PFParticles
    
    // Make a one-to-one TP -> P assn and look for one-to-many assns
/*
    // TODO: This check may not be necessary
    for(auto& tj : slc.tjs) {
      if(tj.AlgMod[kKilled]) continue;
      for(auto& tp : tj.Pts) {
        if(tp.InPFP != 0) {
          std::cout<<"not 0\n";
          exit(1);
        }
      } // tp
    } // tj
*/
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(!tp3d.IsGood) continue;
        if(tp3d.TjID <= 0) {
          std::cout<<"RTP: oops\n";
          continue;
        } // oops
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        if(tp.InPFP > slc.pfps.size()) {
          std::cout<<"RTP: oops1\n";
          continue;
        }
        if(tp.InPFP > 0) {
          // an assn exists. Set the overlap bit and check consistency
          tp.Environment[kEnvOverlap] = true;
          auto& oldp = slc.pfps[tp.InPFP - 1];
          // find the TP3D index
          unsigned short otp = 0;
          for(otp = 0; otp < oldp.TP3Ds.size(); ++otp) {
            auto& otp3d = oldp.TP3Ds[otp];
            if(otp3d.TjID == tp3d.TjID && otp3d.TPIndex == tp3d.TPIndex) break;
          }
          if(otp == oldp.TP3Ds.size()) {
            std::cout<<"RTP: oops2\n";
            continue;
          } // oops
          auto& otp3d = oldp.TP3Ds[otp];
          // keep the old assn and remove the new one
          if(tp3d.VLAInduction && !otp3d.VLAInduction) tp3d.IsGood = false;
          if(!tp3d.VLAInduction && otp3d.VLAInduction) {
            // keep the new assn and remove the old one
            otp3d.IsGood = false;
            tp.InPFP = pfp.ID;
          }
          // both vlaInduction. Keep the first
          if(tp3d.VLAInduction && otp3d.VLAInduction) tp3d.IsGood = false;
          if(!tp3d.VLAInduction && !otp3d.VLAInduction) {
            std::cout<<"RTPs: TP "<<PrintPos(slc, tp)<<" one-to-many -> P"<<tp.InPFP<<" and P"<<tp3d.TjID<<"\n";
            tp3d.IsGood = false;
          }
        } else {
          // no assn exists
          tp.InPFP = pfp.ID;
        } // tp.InPFP > 0
      } // ipt
    } // pfp
  } // ReconcileTPs

  /////////////////////////////////////////
  void MakePFPTjs(TCSlice& slc)
  {
    // This function clobbers all of the tjs that are used in TP3Ds in the pfp and replaces
    // them with new tjs that have a consistent set of TPs to prepare for putting them
    // into the event. Note that none of the Tjs are attached to 2D vertices.
    if(!tcc.useAlg[kMakePFPTjs]) return;
    
    // kill trajectories
    for(auto& pfp : slc.pfps) {
      if(pfp.ID <= 0) continue;
      for(auto& tp3d : pfp.TP3Ds) {
        if(!tp3d.IsGood) continue;
        if(tp3d.TjID > 0) {
          unsigned int itj = tp3d.TjID - 1;
          auto& tj = slc.tjs[itj];
          if(tj.AlgMod[kKilled]) continue;
          MakeTrajectoryObsolete(slc, itj);
        }
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
      // and is 3D matched
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
        } //
      } // tp3d
      // finish defining each of the Tjs and store them
      for(auto& tj : ptjs) {
        if(tj.Pts.size() < 2) continue;
        tj.PDGCode = pfp.PDGCode;
        tj.MCSMom = MCSMom(slc, tj);
        if(!StoreTraj(slc, tj)) {
          std::cout<<"MakePFPTjs: StoreTraj failed P"<<pfp.ID<<" T"<<tj.ID<<"\n";
          continue;
        } // StoreTraj failed
        // associate it with the pfp
        auto& newTj = slc.tjs[slc.tjs.size() - 1];
        pfp.TjIDs.push_back(newTj.ID);
      } // tj
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
    
    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;
    
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

    double yzcut = 1.5 * tcc.wirePitch;

    for(unsigned int ipt = 0; ipt < slc.mallTraj.size() - 1; ++ipt) {
      auto& iTjPt = slc.mallTraj[ipt];
      // length cut
      if(iTjPt.npts < minPts) continue;
      auto& itp = slc.tjs[iTjPt.id - 1].Pts[iTjPt.ipt];
      unsigned short iPlane = iTjPt.plane;
      unsigned int iWire = itp.Pos[0];
      for(unsigned int jpt = ipt + 1; jpt < slc.mallTraj.size() - 1; ++jpt) {
        auto& jTjPt = slc.mallTraj[jpt];
        // ensure that the planes are different
        if(jTjPt.plane == iTjPt.plane) continue;
        // length cut
        if(jTjPt.npts < minPts) continue;
        // check for x range overlap. We know that jTjPt.xlo is >= iTjPt.xlo because of the sort
        if(jTjPt.xlo > iTjPt.xhi) continue;
        // break out if the x range difference becomes large (5 cm)
        if(jTjPt.xlo > iTjPt.xhi + 5) break;
        auto& jtp = slc.tjs[jTjPt.id - 1].Pts[jTjPt.ipt];
        unsigned short jPlane = jTjPt.plane;
        unsigned int jWire = jtp.Pos[0];
        Point3_t ijPos;
//        ijPos[0] = 0.5 * (itp.Pos[0] + jtp.Pos[0]);
        ijPos[0] = itp.Pos[0];
        if(!tcc.geom->IntersectionPoint(iWire, jWire, iPlane, jPlane, cstat, tpc, ijPos[1], ijPos[2])) continue;
        // count weight is one for a two-plane match
//        float cntWght = 1;
        if(numPlanes == 3) {
          // numPlanes == 3
          for(unsigned int kpt = jpt + 1; kpt < slc.mallTraj.size(); ++kpt) {
            auto& kTjPt = slc.mallTraj[kpt];
            // ensure that the planes are different
            if(kTjPt.plane == iTjPt.plane || kTjPt.plane == jTjPt.plane) continue;
           if(kTjPt.xlo > iTjPt.xhi) continue;
            // break out if the x range difference becomes large
            if(kTjPt.xlo > iTjPt.xhi + 5) break;
            auto& ktp = slc.tjs[kTjPt.id - 1].Pts[kTjPt.ipt];
            unsigned short kPlane = kTjPt.plane;
            unsigned int kWire = ktp.Pos[0];
            Point3_t ikPos;
            ikPos[0] = ktp.Pos[0];
            if(!tcc.geom->IntersectionPoint(iWire, kWire, iPlane, kPlane, cstat, tpc, ikPos[1], ikPos[2])) continue;
            if(std::abs(ijPos[1] - ikPos[1]) > yzcut) continue;
            if(std::abs(ijPos[2] - ikPos[2]) > yzcut) continue;
            // we have a match.
            // Just fill temp. See if the Tj IDs are in the match list.
            // first check the input matVec
            bool gotit = false;
            for(auto& ms : inMatVec) {
              if(ms.TjIDs.size() != 3) continue;
              if(iTjPt.id == ms.TjIDs[iPlane] && jTjPt.id == ms.TjIDs[jPlane] && kTjPt.id == ms.TjIDs[kPlane]) {
                gotit = true;
                break;
              }
            } // ms
            if(gotit) continue;
            // next check the temp vector
            unsigned short indx = 0;
            for(indx = 0; indx < temp.size(); ++indx) {
              auto& ms = temp[indx];
              if(iTjPt.id != ms.TjIDs[iPlane]) continue;
              if(jTjPt.id != ms.TjIDs[jPlane]) continue;
              if(kTjPt.id != ms.TjIDs[kPlane]) continue;
              ++ms.Count;
              // Stop checking these Tjs if the count is very high
              if(ms.Count > tcc.match3DCuts[1]) {
                // set npts = 0 for all tj2pts in these Tjs so that they
                // aren't considered anymore
                for(auto& tj2pt : slc.mallTraj) {
                  if(tj2pt.id == ms.TjIDs[iPlane]) tj2pt.npts = 0;
                  if(tj2pt.id == ms.TjIDs[jPlane]) tj2pt.npts = 0;
                  if(tj2pt.id == ms.TjIDs[kPlane]) tj2pt.npts = 0;
                } // lpt
              } // high count
              break;
            } // indx
            if(indx == temp.size()) {
              // not found in the match vector so add it
              MatchStruct ms;
              ms.TjIDs.resize(3);
              // Note that we can put the Tj IDs in plane-order since there are 3 of them
              // This is not the case when there are 2 planes
              ms.TjIDs[iPlane] = iTjPt.id;
              ms.TjIDs[jPlane] = jTjPt.id;
              ms.TjIDs[kPlane] = kTjPt.id;
              ms.Count = 1;
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
            unsigned short kPlane = 3 - iPlane - jPlane;
            float fkwire = tcc.geom->WireCoordinate(ijPos[1], ijPos[2], kPlane, tpc, cstat);
            if(fkwire < 0 || fkwire > tcc.maxPos0[kPlane]) continue;
            TrajPoint tpk;
            tpk.CTP = EncodeCTP(cstat, tpc, kPlane);
            tpk.Pos[0] = fkwire;
            float xp = 0.5 * (iTjPt.xlo + iTjPt.xhi);
            tpk.Pos[1] = tcc.detprop->ConvertXToTicks(xp, kPlane, tpc, cstat) * tcc.unitsPerTick;
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

/*
  /////////////////////////////////////////
  void FixOrder(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Check the wire order of the 2D TPs and fix them.
    // TODO: This is an early algorithm that might not be useful now
    
    if(IsShowerLike(slc, pfp.TjIDs)) return;
    
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
      std::cout<<"FixOrder found a problem: P"<<pfp.ID<<" plane "<<plane<<" nmi "<<nmi<<" npl "<<npl<<" ordr "<<ordr;
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
*/
  /////////////////////////////////////////
  bool Update(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // This function only updates SectionFits that need to be re-sorted or re-fit. It returns
    // false if there was a serious error
    if(pfp.TP3Ds.empty()) return false;
    
    // make some checks
    unsigned short nsu = 0;
    for(auto& sf : pfp.SectionFits) {
      if(sf.NeedsUpdate) ++nsu;
    } // sf
    if(!pfp.NeedsUpdate && nsu > 0) {
      std::cout<<"Update: P"<<pfp.ID<<" has "<<nsu<<" sections that need an update, but the PFP doesn't...\n";
    }

    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      if(sf.NPts == 0) continue;
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
    
    if(prt) PrintTP3Ds("RS", slc, pfp, -1);
    
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
    // Scale by the MCSMom to have 50 points in a plane in each section for a PFParticle
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
        unsigned short nBadPts = 0, firstBadPt;
        CountBadPoints(slc, pfp, fromPt, toPt, nBadPts, firstBadPt);
        if(prt) mf::LogVerbatim("TC")<<"  nit "<<nit<<" sf.ChiDOF "<<sf.ChiDOF<<" sf.NPts "<<sf.NPts<<" fromPt "<<fromPt<<" toPt "<<toPt<<" nbadPts "<<nBadPts<<" firstbadPt "<<firstBadPt<<" prevChi "<<prevChi;
        if(sf.ChiDOF == prevChi) break;
        if(nBadPts > 2) break;
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
    pfp.NeedsUpdate = false;
    
    return true;
  } // resection
  
  /////////////////////////////////////////
  void CountBadPoints(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, unsigned short& nBadPts, unsigned short& firstBadPt)
  {
    // Count the number of points whose pull exceeds tcc.match3DCuts[4]
    firstBadPt = USHRT_MAX;
    nBadPts = 0;
    if(fromPt > pfp.TP3Ds.size() - 1) {
      nBadPts = USHRT_MAX;
      return;
    }
    if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
    bool first = true;
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      if(tp3d.VLAInduction) continue;
      // don't clobber a point if it is on a TP that is overlapping another Tj. This will
      // happen for points close to a vertex and when trajectories cross
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvOverlap]) continue;
      if(PointPull(pfp, tp3d) < tcc.match3DCuts[4]) continue;
      ++nBadPts;
      if(first) {
        first = false;
        firstBadPt = ipt;
      }
    } // ipt
  } // CountBadPoints
  
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
      if(tp3d.VLAInduction) continue;
      // don't clobber a point if it is on a TP that is overlapping another Tj. This will
      // happen for points close to a vertex and when trajectories cross
      auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
      if(tp.Environment[kEnvOverlap]) continue;
      float pull = PointPull(pfp, tp3d);
      if(pull < maxPull) continue;
      maxPull = pull;
      badPt = ipt;
    } // ipt
    if(badPt == USHRT_MAX) return false;
    // clobber that point
    pfp.TP3Ds[badPt].IsGood = false;
    if(prt) {
      mf::LogVerbatim myprt("TC");
      auto& tp3d = pfp.TP3Ds[badPt];
      if(tp3d.TjID > 0) {
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        myprt<<" KBP: P"<<pfp.ID<<" clobber ipt "<<badPt<<" "<<PrintPos(slc, tp)<<" "<<maxPull;
      } else {
        myprt<<" KBP: P"<<pfp.ID<<" clobber ipt "<<badPt<<" "<<maxPull;
      }
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
      unsigned short plane = DecodeCTP(slc.tjs[tp3d.TjID - 1].CTP).Plane;
      ++cntInPln[plane];
      unsigned short enufInPlane = 0;
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) if(cntInPln[plane] >= min2DPts) ++enufInPlane;
      if(enufInPlane > 1) return ipt + 1;
      if(dir < 0 && ipt == 0) break;
    } // ipt
    return USHRT_MAX;
  } // Find3DRecoRange

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
      A[indx][0] = ocs[plane][1] / tp3d.TPXErr2;
      A[indx][1] = ocs[plane][2] / tp3d.TPXErr2;
      A[indx][2] = ocs[plane][1] * x / tp3d.TPXErr2;
      A[indx][3] = ocs[plane][2] * x / tp3d.TPXErr2;
      w[indx] = (tp3d.Wire - ocs[plane][0]) / tp3d.TPXErr2;
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
    
    if(tcc.vtx3DCuts.size() < 3) return;
    
    if(pfp.NeedsUpdate) {
      std::cout<<"ReconcileVertices: P"<<pfp.ID<<" NeedsUpdate is true at entry. Fixing it...\n";
      if(!Update(slc, pfp, prt)) {
        if(prt) mf::LogVerbatim("TC")<<" Update failed inside ReconcileVertices"; 
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
    AttachToAnyVertex(slc, pfp, tcc.vtx3DCuts[2], prt);
    // check for differences and while we are here, see if the pfp was attached
    // to a neutrino vertex and the direction is wrong
    int neutrinoVx = 0;
    if(!slc.pfps.empty()) {
      auto& npfp = slc.pfps[0];
      bool neutrinoPFP = (npfp.PDGCode == 12 || npfp.PDGCode == 14);
      if(neutrinoPFP) neutrinoVx = npfp.Vx3ID[0];
    } // pfps exist
    unsigned short neutrinoVxEnd = 2;
    for(unsigned short end = 0; end < 2; ++end) {
      // see if a vertex got attached
      if(pfp.Vx3ID[end] <= 0) continue;
      if(pfp.Vx3ID[end] == neutrinoVx) neutrinoVxEnd = end;
      // see if this is a vertex in the list using the T -> 2V -> 3V assns
      if(std::find(vx3List.begin(), vx3List.end(), pfp.Vx3ID[end]) != vx3List.end()) continue;
//      std::cout<<"RV: P"<<pfp.ID<<" was attached to 3V"<<pfp.Vx3ID[end]<<" but a P -> T -> 2V -> 3V assn exists. Write some code to clobber this assn or deal with it somehow.\n";
    } // end
    if(neutrinoVxEnd < 2 && neutrinoVxEnd != 0) {
//      std::cout<<"RV: reversing P"<<pfp.ID<<"\n";
      Reverse(slc, pfp);
    }
    
    return;
  } // ReconcileVertices

  /////////////////////////////////////////
  void TrimEndPts(TCSlice& slc, PFPStruct& pfp, bool prt)
  {
    // Check for a wire signal and trim points starting
    // at the ends until there at least 2 planes that have a wire signal. The
    // sections that have points removed are re-fit without those points and
    // another check is made in a second iteration

    if(pfp.ID <= 0) return;
    // Trimming short tracks that are barely reconstructable isn't a good idea
    if(pfp.TP3Ds.size() < 10) return;
    // don't trim shower-like pfps
    if(IsShowerLike(slc, pfp.TjIDs)) return;
    
    pfp.NeedsUpdate = false;
    for(unsigned short nit = 0; nit < 2; ++nit) {
      // create a vector valid points
      std::vector<bool> validPt(pfp.TP3Ds.size(), true);
      // inspect end 0
      for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(!tp3d.IsGood) continue;
        unsigned short cnt = 0;
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, tp3d.Pos, inCTP);
          if(SignalAtTp(tp)) ++cnt;
        } // plane
        if(cnt >= 2) break;
//        if(prt) mf::LogVerbatim("TC")<<"TEP: P"<<pfp.ID<<" nit "<<nit<<" trim TP3D "<<tp3d.TjID<<"_"<<tp3d.TPIndex;
        validPt[ipt] = false;
        pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
        pfp.NeedsUpdate = true;
      } // ipt
      // inspect the other end
      for(unsigned short ipt = pfp.TP3Ds.size() - 1; ipt > 0; --ipt) {
        auto& tp3d = pfp.TP3Ds[ipt];
        if(!tp3d.IsGood) continue;
        unsigned short cnt = 0;
        for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
          CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
          auto tp = MakeBareTP(slc, tp3d.Pos, inCTP);
          if(SignalAtTp(tp)) ++cnt;
        } // plane
        if(cnt >= 2) break;
//        if(prt) mf::LogVerbatim("TC")<<"TEP: P"<<pfp.ID<<" nit "<<nit<<" trim TP3D "<<tp3d.TjID<<"_"<<tp3d.TPIndex;
        validPt[ipt] = false;
        pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
        pfp.NeedsUpdate = true;
      } // ipt
      if(pfp.NeedsUpdate) {
        // trim the points
        // find the first good point
        unsigned short firstGood = 0;
        for(firstGood = 0; firstGood < pfp.TP3Ds.size(); ++firstGood) if(validPt[firstGood]) break;
        // and the last good point
        unsigned short lastGood = pfp.TP3Ds.size();
        for(lastGood = pfp.TP3Ds.size() - 1; lastGood > 0; --lastGood) if(validPt[lastGood]) break;
        ++lastGood;
        if(firstGood == 0 && lastGood == pfp.TP3Ds.size()) break;
        std::vector<TP3D> temp(pfp.TP3Ds.begin() + firstGood, pfp.TP3Ds.begin() + lastGood);
        pfp.TP3Ds = temp;
        if(!Update(slc, pfp, prt)) {
          std::cout<<"Update P"<<pfp.ID<<" failed in TrimEndPts\n";
        }
      }
    } // nit
    
    if(pfp.NeedsUpdate) Update(slc, pfp, prt);
    
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
    // Look for gaps in each plane in the TP3Ds vector in planes in which 
    // the projection of the pfp angle is large (~> 60 degrees). Hits
    // reconstructed at large angles are poorly reconstructed which results
    // in poorly reconstructed 2D trajectories
    
    if(pfp.ID <= 0) return;
    if(pfp.TP3Ds.empty()) return;
    if(pfp.SectionFits.empty()) return;
    if(!tcc.useAlg[kFillGaps3D]) return;
    
    // Only print if MVI is set
    prt = (tcc.dbgPFP && pfp.MVI == debug.MVI);
    
    // project the first 3D point into a TP in each plane and decide
    // whether to continue;
    auto& fromTP3 = pfp.TP3Ds[0];
    TrajPoint fromTP;
    // Use a delta cut instead of a pull cut
    float maxPull = -1;
    for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
      CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
      fromTP = MakeBareTP(slc, fromTP3.Pos, fromTP3.Dir, inCTP);
      double dang = DeltaAngle(0., fromTP.Ang);
      if(std::abs(dang) < 1) continue;
      float expectedRMS = ExpectedHitsRMS(slc, fromTP);
      if(prt) mf::LogVerbatim("TC")<<"FillGaps3D: P"<<pfp.ID<<" Search for gaps in plane "<<plane<<" angle "<<fromTP.Ang<<" expected rms (ticks) "<<expectedRMS;
      AddPointsInRange(slc, pfp, 0, pfp.TP3Ds.size(), inCTP, maxPull, expectedRMS, prt);
    } // plane
    
    if(!pfp.NeedsUpdate) return;
    
    if(!Update(slc, pfp, prt)) {
      std::cout<<"FillGaps3D: Update failed\n";
    }
    
  } // FillGaps3D

  /////////////////////////////////////////
  void AddPointsInRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, 
                        CTP_t inCTP, float maxPull, float maxDelta, bool prt)
  {
    // Try to insert 2D trajectory points into the 3D trajectory point vector pfp.TP3Ds. 
    // This function inserts new Tp3Ds and sets the NeedsUpdate flags true. 
    // The calling function should call Update
    if(fromPt > toPt) return;
    if(toPt > pfp.TP3Ds.size()) toPt = pfp.TP3Ds.size();
    
    TrajPoint ltp;
    // having multiple hits or TPs on a wire is OK but we need to ensure that
    // no wires were missed
    std::vector<unsigned int> wiresFound;
    std::vector<std::pair<int, unsigned short>> tpsFound;
    // populate the vectors
    for(auto& tp3d : pfp.TP3Ds) {
      if(tp3d.CTP != inCTP) continue;
      if(tp3d.TjID > 0) {
        auto npr = std::make_pair(tp3d.TjID, tp3d.TPIndex);
        tpsFound.push_back(npr);
      }
      unsigned int wire = std::nearbyint(tp3d.Wire);
      if(std::find(wiresFound.begin(), wiresFound.end(), wire) == wiresFound.end()) wiresFound.push_back(wire);
    } // tp3d
    if(wiresFound.size() < 2) return;
    // flag TP3Ds on very-large-angle induction plane hits since they are poorly
    // reconstructed and would bias the fit
    bool vlaInduction = false;
    if(DecodeCTP(inCTP).Plane < slc.nPlanes - 1) {
      // get the 2D angle from the first section
      auto& sf = pfp.SectionFits[0];
      ltp = MakeBareTP(slc, sf.Pos, sf.Dir, inCTP);
      vlaInduction = (std::abs(ltp.Ang) > 1);
    }
    // iterate through the vector and project each point into inPlane.
    for(unsigned short ipt = fromPt; ipt < toPt; ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(!tp3d.IsGood) continue;
      // check for already found TPs
      if(tp3d.TjID > 0) {
        auto npr = std::make_pair(tp3d.TjID, tp3d.TPIndex);
        if(std::find(tpsFound.begin(), tpsFound.end(), npr) != tpsFound.end()) continue;
        tpsFound.push_back(npr);
      }
      ltp = MakeBareTP(slc, tp3d.Pos, inCTP);
//      if(tcc.dbgPFP && pfp.MVI == debug.MVI) std::cout<<"TP "<<PrintPos(slc, ltp)<<"\n";
      unsigned int wire = std::nearbyint(ltp.Pos[0]);
      // check the list of found wires
      if(std::find(wiresFound.begin(), wiresFound.end(), wire) != wiresFound.end()) continue;
      wiresFound.push_back(wire);
      // 
      if(!FindCloseHits(slc, ltp, maxDelta, kUsedHits)) continue;
      for(auto iht : ltp.Hits) {
        if(slc.slHits[iht].InTraj <= 0) continue;
        // this hit is used in a TP so find the tpIndex
        auto& tj = slc.tjs[slc.slHits[iht].InTraj - 1];
        unsigned short tpIndex = 0;
        for(tpIndex = tj.EndPt[0]; tpIndex <= tj.EndPt[1]; ++tpIndex) {
          auto& tp = tj.Pts[tpIndex];
          if(tp.Chg <= 0) continue;
          // This doesn't check for UseHit true but that is probably ok here
          if(std::find(tp.Hits.begin(), tp.Hits.end(), iht) != tp.Hits.end()) break;
        } // ipt
        if(tpIndex > tj.EndPt[1]) continue;
        auto npr = std::make_pair(tj.ID, tpIndex);
        if(std::find(tpsFound.begin(), tpsFound.end(), npr) != tpsFound.end()) continue;
        tpsFound.push_back(npr);
        auto& tp = tj.Pts[tpIndex];
        // see if it is used in a different PFP
        if(tp.InPFP > 0) continue;
        auto newTP3D = CreateTP3D(slc, tj.ID, tpIndex, vlaInduction);
        if(!SetSection(slc, pfp, newTP3D)) {
          std::cout<<"APIR: SetSection failed\n";
          continue;
        }
        float pull = PointPull(pfp, newTP3D);
        if(!vlaInduction && pull > maxPull) continue;
        if(prt && pfp.MVI == debug.MVI) mf::LogVerbatim("TC")<<"APIR: P"<<pfp.ID<<" TP "<<PrintPos(slc, tp)<<" pull "<<pull<<" dx "<<newTP3D.TPX - newTP3D.Pos[0];
        if(!InsertTP3D(pfp, newTP3D)) continue;
      } // iht
      
    } // ipt
    
  } // AddPointsInRange

  /////////////////////////////////////////
  bool InsertTP3D(PFPStruct& pfp, TP3D& tp3d)
  {
    // inserts the tp3d into the section defined by tp3d.SFIndex
    if(tp3d.SFIndex >= pfp.SectionFits.size()) return false;
    // Find the first occurrence of this SFIndex
    unsigned short ipt = 0;
    for(ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) if(tp3d.SFIndex == pfp.TP3Ds[ipt].SFIndex) break;
    if(ipt == pfp.TP3Ds.size()) return false;
    // next see if we can insert it so that re-sorting of this section isn't required
    auto& along = tp3d.along;
    for(unsigned short iipt = ipt; iipt < pfp.TP3Ds.size() - 1; ++iipt) {
      // break out if the next point is in a different section
      if(pfp.TP3Ds[iipt + 1].SFIndex != tp3d.SFIndex) break;
      if(along > pfp.TP3Ds[iipt].along && along < pfp.TP3Ds[iipt + 1].along) {
        ipt = iipt + 1;
        break;
      }
    } // iipt
    pfp.TP3Ds.insert(pfp.TP3Ds.begin() + ipt, tp3d);
    pfp.SectionFits[tp3d.SFIndex].NeedsUpdate = true;
    pfp.NeedsUpdate = true;
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
    // See if the along variable is monotonically increasing
    float prevAlong = 0;
    bool first = true;
    bool needsSort = false;
    for(unsigned short ii = 0; ii < pfp.TP3Ds.size(); ++ii) {
      auto& tp3d = pfp.TP3Ds[ii];
      if(tp3d.SFIndex != sfIndex) continue;
      if(first) {
        first = false;
        prevAlong = tp3d.along;
      } else {
        if(tp3d.along < prevAlong) needsSort = true;
        prevAlong = tp3d.along;
      }
      temp.push_back(tp3d);
      indx.push_back(ii);
    } // tp3d
    if(temp.empty()) return false;
    // no sort needed?
    if(temp.size() == 1) return true;
    if(!needsSort) {
//      std::cout<<"SortSection: P"<<pfp.ID<<" section "<<sfIndex<<" doesn't need sorting\n";
      sf.NeedsUpdate = false;
      return true;
    }    
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
  void MakeTP3Ds(TCSlice& slc, PFPStruct& pfp)
  {
    // Create and populate the TP3Ds vector. This function is called before the first
    // fit is done so the TP3D along variable can't be determined
    if(!pfp.TP3Ds.empty() || pfp.SectionFits.size() != 1) {
      std::cout<<"MakeTP3Ds: invalid call P"<<pfp.ID<<". TP3Ds is not empty or SectionFits size != 1\n";
      return;
    } // error
//    auto& sf = pfp.SectionFits[0];
    // Ddd the points associated with the Tjs that were used to create the PFP
    for(auto tid : pfp.TjIDs) {
      auto& tj = slc.tjs[tid - 1];
      auto& tp0 = tj.Pts[tj.EndPt[0]];
      bool vlaInduction = (DecodeCTP(tj.CTP).Plane < slc.nPlanes - 1 && std::abs(tp0.Ang) > 1);
      for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
        auto& tp = tj.Pts[ipt];
        if(tp.Chg <= 0) continue;
        if(tp.InPFP > 0) continue;
        auto tp3d = CreateTP3D(slc, tid, ipt, vlaInduction);
        tp3d.SFIndex = 0;
        // We need to assume that all points are good or the first fit will fail
        tp3d.IsGood = true;
        pfp.TP3Ds.push_back(tp3d);
      } // ipt
    } // tid
  } // MakeTP3Ds

  /////////////////////////////////////////
  void Reverse(TCSlice& slc, PFPStruct& pfp)
  {
    // reverse the PFParticle
    std::reverse(pfp.TP3Ds.begin(), pfp.TP3Ds.end());
    std::reverse(pfp.SectionFits.begin(), pfp.SectionFits.end());
    for(unsigned short sfi = 0; sfi < pfp.SectionFits.size(); ++sfi) {
      auto& sf = pfp.SectionFits[sfi];
      // flip the direction vector
      for(unsigned short xyz = 0; xyz < 3; ++xyz) sf.Dir[xyz] *= -1;
    } // sf
    // correct the along variable
    for(auto& tp3d : pfp.TP3Ds) tp3d.along *= -1;
    std::swap(pfp.dEdx[0], pfp.dEdx[1]);
    std::swap(pfp.dEdxErr[0], pfp.dEdxErr[1]);
    std::swap(pfp.Vx3ID[0], pfp.Vx3ID[1]);
    std::swap(pfp.EndFlag[0], pfp.EndFlag[1]);
  } // Reverse

  /////////////////////////////////////////
  void FillmAllTraj(TCSlice& slc) 
  {
    // Fills the mallTraj vector with trajectory points in the tpc and sorts
    // them by increasing X
    
    int cstat = slc.TPCID.Cryostat;
    int tpc = slc.TPCID.TPC;

    // define mallTraj
    Tj2Pt tj2pt;
    unsigned short cnt = 0;
  } // FillmAllTraj

  //////////////////////////////////////////
  void FillmAllTraj(TCSlice& slc)
  {

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
        if(tp.Chg <= 0) continue;
        if(tp.Pos[0] < -0.4) continue;
        tj2pt.wire = std::nearbyint(tp.Pos[0]);
        ++cnt;
        bool hasWire = tcc.geom->HasWire(geo::WireID(cstat, tpc, plane, tj2pt.wire));
        // don't try matching if the wire doesn't exist
        if(!hasWire) continue;
        float xpos = tcc.detprop->ConvertTicksToX(tp.Pos[1]/tcc.unitsPerTick, plane, tpc, cstat);
        float posPlusRMS = tp.Pos[1] + TPHitsRMSTime(slc, tp, kUsedHits);
        float rms = tcc.detprop->ConvertTicksToX(posPlusRMS/tcc.unitsPerTick, plane, tpc, cstat) - xpos;
        if(rms < tcc.match3DCuts[0]) rms = tcc.match3DCuts[0];
        tj2pt.xlo = xpos - rms;
        tj2pt.xhi = xpos + rms;
        tj2pt.plane = plane;
        tj2pt.id = tjID;
        tj2pt.ipt = ipt;
        tj2pt.npts = tj.EndPt[1] - tj.EndPt[0] + 1;
        slc.mallTraj.push_back(tj2pt);
      } // tp
    } // tj
    
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
    
    // square of the maximum length that is used for finding the average dE/dx
    float maxSep2 = 5 * tcc.wirePitch;
    maxSep2 *= maxSep2;
    
    for(unsigned short end = 0; end < numEnds; ++end) {
      std::vector<float> cnt(slc.nPlanes);
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
        if(PosSep2(tp3d.Pos, endPos) > maxSep2) break;
        // require good points
        if(!tp3d.IsGood) continue;
        float dedx = dEdx(slc, tp3d);
        if(dedx < 0.5) continue;
        unsigned short plane = DecodeCTP(tp3d.CTP).Plane;
        pfp.dEdx[end][plane] += dedx;
        ++cnt[plane];
      } // ii
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        if(cnt[plane] == 0) continue;
        pfp.dEdx[end][plane] /= cnt[plane];
      } // plane
    } // end

  } // FilldEdx

  /////////////////////////////////////////
  float dEdx(TCSlice& slc, TP3D& tp3d)
  {
    if(!tp3d.IsGood) return 0;
    if(tp3d.TjID > slc.slHits.size()) return 0;
    if(tp3d.TjID <= 0) return 0;

    double dQ = 0.;
    double time = 0;
    geo::PlaneID plnID;
    auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(!tp.UseHit[ii]) continue;
      auto& hit = (*evt.allHits)[slc.slHits[tp.Hits[ii]].allHitsIndex];
      dQ += hit.Integral();
    } // ii
    time = tp.Pos[1] / tcc.unitsPerTick;
    plnID = DecodeCTP(tp.CTP);
    if(dQ == 0) return 0;
    double angleToVert = tcc.geom->Plane(plnID).ThetaZ() - 0.5 * ::util::pi<>();
    double cosgamma = std::abs(std::sin(angleToVert) * tp3d.Dir[1] + std::cos(angleToVert) * tp3d.Dir[2]);
    if(cosgamma == 0) return 0;
    double dx = tcc.geom->WirePitch(plnID) / cosgamma;
    double dQdx = dQ / dx;
    double t0 = 0;
    return tcc.caloAlg->dEdx_AREA(dQdx, time, plnID.Plane, t0);
  } // dEdx

  ////////////////////////////////////////////////
  TP3D CreateTP3D(TCSlice& slc, int tjID, unsigned short tpIndex, bool vlaInduction)
  {
    // create a TP3D with a single TP. Note that the SectionFit in which it
    // should be placed and the 3D position can't be determined until the the TP3D is 
    // associated with a pfp. See SetSection()
    
    if(tjID > slc.tjs.size()) {
      std::cout<<"bad tjID\n";
      exit(1);
    }
    auto& tj = slc.tjs[tjID - 1];
    if(tpIndex < tj.EndPt[0] || tpIndex > tj.EndPt[1]) {
      std::cout<<"bad tpIndex\n";
      exit(1);
    }
    
    TP3D tp3d;
    tp3d.TjID = tjID;
    tp3d.TPIndex = tpIndex;
    auto& tp2 = tj.Pts[tp3d.TPIndex];
    auto plnID = DecodeCTP(tp2.CTP);
    tp3d.CTP = tp2.CTP;
    double tick = tp2.HitPos[1]/tcc.unitsPerTick;
    tp3d.TPX = tcc.detprop->ConvertTicksToX(tick, plnID);
    if(vlaInduction) {
      tp3d.TPXErr2 = 1.;
      tp3d.VLAInduction = true;
    } else {
      // Find the rms in WSE units assuming that all hits in the TP are used in the TP
      float rms = TPHitsRMSTime(slc, tp2, kAllHits);
      // convert to X
      rms *= tcc.wirePitch;
      tp3d.TPXErr2 = rms * rms;
      tp3d.VLAInduction = false;
    }
    tp3d.Wire = tp2.Pos[0];
    // Can't declare it good since Pos and SFIndex aren't defined
    tp3d.IsGood = false;
    return tp3d;
  } // CreateTP3D

  /////////////////////////////////////////
  bool SetSection(TCSlice& slc, PFPStruct& pfp, TP3D& tp3d)
  {
<<<<<<< HEAD
    // Determine which SectionFit this tp3d should reside, then calculate
    // the 3D position and the distance fromt the center of the SectionFit

    if(tp3d.Wire < 0) return false;
    if(pfp.SectionFits.empty()) return false;
    if(pfp.SectionFits[0].Pos[0] == -10.0) {
      std::cout<<"SetSection: P"<<pfp.ID<<" Can't determine the section when no fit exists\n";
      return false;
    }
    
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
    return true;
  } // SetSection
  
  ////////////////////////////////////////////////
  float PointPull(const PFPStruct& pfp, TP3D& tp3d)
  {
    // returns the pull that the tp3d will cause in the pfp section fit. This
    // currently only uses position but eventually will include charge
    return std::abs(tp3d.Pos[0] - tp3d.TPX) / sqrt(tp3d.TPXErr2);
  } // PointPull

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
      Point3_t startPos;
      if(pfp.TP3Ds.empty()) {
        // must be a neutrino pfp
        startPos = pfp.SectionFits[0].Pos;
      } else if(!pfp.TP3Ds.empty()) {
        // normal pfp
        startPos = pfp.TP3Ds[0].Pos;
      }
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
    bool neutrinoPFP = pfp.PDGCode == 12 || pfp.PDGCode == 14;
    
    float abit = 5;
    Point3_t pos;
    if(neutrinoPFP) {
      pos = pfp.SectionFits[0].Pos;
    } else if(end == 0) {
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
  bool SectionStartEnd(const PFPStruct& pfp, unsigned short sfIndex, unsigned short& startPt, unsigned short& endPt)
  {
    // this assumes that the TP3Ds vector is sorted
    startPt = USHRT_MAX;
    endPt = USHRT_MAX;
    if(sfIndex >= pfp.SectionFits.size()) return false;
    
    bool first = true;
    for(unsigned short ipt = 0; ipt < pfp.TP3Ds.size(); ++ipt) {
      auto& tp3d = pfp.TP3Ds[ipt];
      if(tp3d.SFIndex < sfIndex) continue;
      if(first) {
        first = false;
        startPt = ipt;
      } // first
      if(tp3d.SFIndex > sfIndex) break;
      endPt = ipt;
    } // ipt
    return true;
    
  } // SectionStartEnd

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
        if(SectionStartEnd(pfp, sfi, startPt, endPt)) {
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
      myprt<<someText<<"  ipt SFI ________Pos________  Delta Pull GI  along dE/dx  T_ipt_P:W:T  Signal? MCPIndex\n";
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
      float pull = PointPull(pfp, tp3d);
      myprt<<std::setprecision(1)<<std::setw(6)<<pull;
      myprt<<std::setw(3)<<tp3d.IsGood<<tp3d.VLAInduction;
      myprt<<std::setw(7)<<std::setprecision(1)<<tp3d.along;
      myprt<<std::setw(6)<<std::setprecision(2)<<dEdx(slc, tp3d);
      if(tp3d.TjID > 0) {
        auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
        myprt<<" T"<<tp3d.TjID<<"_"<<tp3d.TPIndex<<"_"<<PrintPos(slc, tp);
      } else {
        myprt<<" UNDEFINED";
      }
      // print SignalAtTP in each plane
      myprt<<" ";
      for(unsigned short plane = 0; plane < slc.nPlanes; ++plane) {
        CTP_t inCTP = EncodeCTP(pfp.TPCID.Cryostat, pfp.TPCID.TPC, plane);
        auto tp = MakeBareTP(slc, tp3d.Pos, inCTP);
        myprt<<" "<<SignalAtTp(tp);
      } // plane
      
      // look for a mcp match
      if(!evt.allHitsMCPIndex.empty()) {
        unsigned int mcpIndex = UINT_MAX;
        if(tp3d.TjID > 0) {
          auto& tp = slc.tjs[tp3d.TjID - 1].Pts[tp3d.TPIndex];
          for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
            if(!tp.UseHit[ii]) continue;
            unsigned ahi = slc.slHits[tp.Hits[ii]].allHitsIndex;
            mcpIndex = evt.allHitsMCPIndex[ahi];
            break;
          } // ii
          if(mcpIndex != UINT_MAX) myprt<<" "<<mcpIndex;
        } // tp3d.TjID > 0
      } // !evt.allHitsMCPIndex.empty()

      myprt<<"\n";
    } // ipt
  } // PrintTP3Ds
} // namespace
