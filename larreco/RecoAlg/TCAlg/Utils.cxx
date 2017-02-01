#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {
  
   ////////////////////////////////////////////////
  void WatchHit(std::string someText, TjStuff& tjs, const unsigned int& wHit, short& wInTraj, const unsigned short& tjID)
  {
    // a temp routine to watch when inTraj changes for the supplied hit index, watchHit
    if(wHit > tjs.fHits.size() - 1) return;
    
    if(tjs.fHits[wHit].InTraj != wInTraj) {
      std::cout<<someText<<" Hit "<<PrintHitShort(tjs.fHits[wHit])<<" was InTraj "<<wInTraj<<" now InTraj "<<tjs.fHits[wHit].InTraj<<" tjID = "<<tjID<<"\n";
      wInTraj = tjs.fHits[wHit].InTraj;
    }
  } // WatchHit

  ////////////////////////////////////////////////
  bool Reverse3DMatchTjs(TjStuff& tjs, unsigned short im, bool prt)
  {
    // Return true if the 3D matched hits in the trajectories in tjs.matchVecPFPList are in the wrong order in terms of the
    // physics standpoint, e.g. dQ/dx, muon delta-ray tag, cosmic rays entering the detector, etc. 
    
    if(im > tjs.matchVecPFPList.size() - 1) return false;
    
    auto& mv = tjs.matchVec[im];

    // through-going track? Check for outside the Fiducial Volume at the start (s) and end (e).
    // These variables assume that the TPC is exposed to a beam that contains muons entering at the front and
    // a background of cosmic rays that enter from the top
    bool sAtSide = (mv.sXYZ[0] < tjs.XLo || mv.sXYZ[0] > tjs.XHi);
    bool sAtTop = (mv.sXYZ[1] > tjs.YHi);
    bool sAtBottom = (mv.sXYZ[1] < tjs.YLo);
    bool sAtFront = (mv.sXYZ[2] < tjs.ZLo);
    bool sAtBack = (mv.sXYZ[2] > tjs.ZHi);
    
    bool eAtSide = (mv.eXYZ[0] < tjs.XLo || mv.eXYZ[0] > tjs.XHi);
    bool eAtTop = (mv.sXYZ[1] > tjs.YHi);
    bool eAtBottom = (mv.sXYZ[1] < tjs.YLo);
    bool eAtFront = (mv.sXYZ[2] < tjs.ZLo);
    bool eAtBack = (mv.sXYZ[2] > tjs.ZHi);
    
    // the start (end) is outside the FV
    bool sOutsideFV = sAtBottom || sAtTop || sAtFront || sAtBack;
    bool eOutsideFV = eAtBottom || eAtTop || eAtFront || eAtBack;
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"sXYZ ("<<(int)mv.sXYZ[0]<<", "<<(int)mv.sXYZ[1]<<", "<<(int)mv.sXYZ[2]<<") sOutsideFV "<<sOutsideFV;
      myprt<<" eXYZ ("<<(int)mv.eXYZ[0]<<", "<<(int)mv.eXYZ[1]<<", "<<(int)mv.eXYZ[2]<<") eOutsideFV "<<eOutsideFV;
    } // prt
    
    if(sOutsideFV && eOutsideFV) {
      // both ends are outside the FV - probably a through-going muon. See if it enters at the top or the front.
      if(sAtFront && eAtBack) return false;
      if(eAtFront && sAtBack) return true;
      // Next consider cosmic rays entering the top
      if(sAtTop && eAtBottom) return false;
      if(eAtTop && sAtBottom) return true;
      // entering/leaving the sides
      if(sAtSide && eAtBottom) return false;
      if(eAtSide && sAtBottom) return true;
      return false;
    } // outside the FV

    // Use a simple voting scheme using charge and muon tag direction
    // ngt is the number of times that something (charge) has the correct behavior (for a stopping track)
    unsigned short ngt = 0;
    unsigned short nlt = 0;
    for(auto& tjID : mv.TjIDs) {
      unsigned short itj = tjID - 1;
      unsigned short endPt0 = tjs.allTraj[itj].EndPt[0];
      unsigned short endPt1 = tjs.allTraj[itj].EndPt[1];
      float chgrat = tjs.allTraj[itj].Pts[endPt1].AveChg / tjs.allTraj[itj].Pts[endPt0].AveChg;
      if(chgrat > 1.1) {
        ++ngt;
      } else if(chgrat < 0.9) {
        ++nlt;
      }
    } // itj
    
    // everything seems to be in proper order
    if(ngt >= nlt) return false;
    
    return true;
    
  } // Reverse3DMatchTjs
  
  ////////////////////////////////////////////////
  unsigned short Matched3DVtx(TjStuff& tjs, unsigned short im)
  {
    // Checks for a 3D vertex associated with trajectoris in the MatchStruct. If one or more are found,
    // define sVtx3DIndex and eVtx3DIndex (if there are a 2 vertices) and return true
    
    if(im > tjs.matchVecPFPList.size() - 1) return 0;
    
    auto& ms = tjs.matchVec[im];
    if(ms.TjIDs.empty()) return 0;
    
    // There should be at most 2 unless there is a problem
    std::vector<unsigned short> vIndex;
    
    for(unsigned short ii = 0; ii < ms.TjIDs.size(); ++ii) {
      unsigned short itj = ms.TjIDs[ii] - 1;
      Trajectory& tj = tjs.allTraj[itj];
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] <= 0) continue;
        // Has a 2D vertex
        unsigned short iv2 = tj.VtxID[end] - 1;
        if(tjs.vtx[iv2].Ptr3D == SHRT_MAX) continue;
        if(tjs.vtx[iv2].Ptr3D < 0) continue;
        // Has a 3D vertex
        unsigned short iv3 = tjs.vtx[iv2].Ptr3D;
        // already in the list?
        if(std::find(vIndex.begin(), vIndex.end(), iv3) != vIndex.end()) continue;
        vIndex.push_back(iv3);
      } // end
    } // ii
    
    if(vIndex.empty()) return 0;
    
    // Need to do something here when there are more than 2 vertices
    if(vIndex.size() > 2) {
      mf::LogVerbatim("TC")<<"MatchHas3DVtx found more than 2 3D vertices. Ignore for now - Write some code";
      vIndex.resize(2);
    }
    
    if(vIndex.size() == 2) {
      // Determine which should be the start vertex. Pick the one at larger X
      if(tjs.vtx3[vIndex[0]].X > tjs.vtx3[vIndex[1]].X) {
        ms.sVtx3DIndex = vIndex[0];
        ms.eVtx3DIndex = vIndex[1];
      } else {
        ms.sVtx3DIndex = vIndex[1];
        ms.eVtx3DIndex = vIndex[0];
      }
      // This shouldn't be necessary but do it anyway
      ms.sXYZ[0] = tjs.vtx3[ms.sVtx3DIndex].X;
      ms.sXYZ[1] = tjs.vtx3[ms.sVtx3DIndex].Y;
      ms.sXYZ[2] = tjs.vtx3[ms.sVtx3DIndex].X;
      ms.eXYZ[0] = tjs.vtx3[ms.eVtx3DIndex].X;
      ms.eXYZ[1] = tjs.vtx3[ms.eVtx3DIndex].Y;
      ms.eXYZ[2] = tjs.vtx3[ms.eVtx3DIndex].X;
      return 2;
    } // vIndex.size() == 2
    
    // Have 1 3D vertex. Make it the start vertex
    ms.sVtx3DIndex = vIndex[0];
    ms.sXYZ[0] = tjs.vtx3[ms.sVtx3DIndex].X;
    ms.sXYZ[1] = tjs.vtx3[ms.sVtx3DIndex].Y;
    ms.sXYZ[2] = tjs.vtx3[ms.sVtx3DIndex].X;
    return 1;
    
  } // Matched3DVtx

  
  ////////////////////////////////////////////////
  void ReleaseHits(TjStuff& tjs, Trajectory& tj)
  {
    // Sets InTraj[] = 0 and UseHit false for all TPs in work. Called when abandoning work
    for(auto& tp : tj.Pts) {
      for(auto& iht : tp.Hits) {
        if(tjs.fHits[iht].InTraj == tj.ID) tjs.fHits[iht].InTraj = 0;
      }
    } // tp
    
  } // ReleaseWorkHits
  
  //////////////////////////////////////////
  void UnsetUsedHits(TjStuff& tjs, TrajPoint& tp)
  {
    // Sets InTraj = 0 and UseHit false for all used hits in tp
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(tp.UseHit[ii]) {
        tjs.fHits[tp.Hits[ii]].InTraj = 0;
        tp.UseHit[ii] = false;
      } // UseHit
    } // ii
    tp.Chg = 0;
  } // UnsetUsedHits

  //////////////////////////////////////////
  void TrimEndPts(TjStuff& tjs, Trajectory& tj, const std::vector<float>& fQualityCuts, bool prt)
  {
    // Trim the hits off the end until there are at least fMinPts consecutive hits at the end
    // and the fraction of hits on the trajectory exceeds fQualityCuts[0]
    // Minimum length requirement accounting for dead wires where - denotes a wire with a point
    // and D is a dead wire. Here is an example with minPts = 3
    //  ---DDDDD--- is OK
    //  ----DD-DD-- is OK
    //  ----DDD-D-- is OK
    //  ----DDDDD-- is not OK
    
    unsigned short minPts = fQualityCuts[1];
    float maxPtSep = minPts + 2;
    if(NumPtsWithCharge(tjs, tj, false) < minPts) return;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"TrimEndPts: minPts "<<minPts<<" required. maxPtSep "<<maxPtSep<<" Minimum hit fraction "<<fQualityCuts[0];
      PrintTrajectory("TEPi", tjs, tj, USHRT_MAX);
    }

    unsigned short newEndPt = tj.EndPt[1];
    unsigned short nPtsWithCharge;
    float hitFrac = 0;
    while(newEndPt > minPts) {
      nPtsWithCharge = 0;
      if(tj.Pts[newEndPt].Chg == 0) {
        --newEndPt;
        continue;
      }
      for(unsigned short jj = 0; jj < minPts; ++jj) {
        unsigned short jpt = newEndPt - jj;
        if(tj.Pts[jpt].Chg > 0) ++nPtsWithCharge; 
        if(jpt < minPts) break;
      } // jj
      float ptSep = std::abs(tj.Pts[newEndPt - minPts].Pos[0] - tj.Pts[newEndPt].Pos[0]);
      if(prt) mf::LogVerbatim("TC")<<" newEndPt "<<newEndPt<<" ptSep "<<ptSep<<" nPtsWithCharge "<<nPtsWithCharge;
      // allow only one dead wire at the end
      if(nPtsWithCharge == minPts && ptSep < maxPtSep) {
        // minPts consecutive points have charge. Check the hit fraction
        float npwc = NumPtsWithCharge(tjs, tj, true, tj.EndPt[0], newEndPt);
        float nwires = std::abs(tj.Pts[tj.EndPt[0]].Pos[0] - tj.Pts[newEndPt].Pos[0]) + 1;
        hitFrac = npwc / nwires;
        if(prt) mf::LogVerbatim("TC")<<" check hitFrac "<<newEndPt<<" nwires "<<(int)nwires<<" npwc "<<(int)npwc<<" hitFrac "<<hitFrac;
        if(hitFrac > fQualityCuts[0]) break;
        newEndPt -= minPts;
      }
      --newEndPt;
    } // newEndPt

    // passed the cuts with no modifications
    if(newEndPt == tj.EndPt[1]) return;

    // newEndPt is now the last point that satisfies these conditions
    // dead wire check
    nPtsWithCharge = 0;
    unsigned short nConsecutivePts = 0;
    for(unsigned short jj = 0; jj < minPts; ++jj) {
      unsigned short jpt = newEndPt - jj;
      if(tj.Pts[jpt].Chg > 0) ++nPtsWithCharge;
      if(jj > 0 && std::abs(tj.Pts[jpt+1].Pos[0] - tj.Pts[jpt].Pos[0]) < 1.5) ++nConsecutivePts;
      if(jpt == 0) break;
    } // jj
    
    if(prt) mf::LogVerbatim("TC")<<" newEndPt "<<newEndPt<<" nConsecutivePts "<<nConsecutivePts<<" Required "<<minPts - 1;
    
    // lop off the last point if the consecutive point condition isn't met and re-calculate
    if(nConsecutivePts < minPts - 1) {
      --newEndPt;
      nPtsWithCharge = 0;
      unsigned short nConsecutivePts = 0;
      for(unsigned short jj = 0; jj < minPts; ++jj) {
        unsigned short jpt = newEndPt - jj;
        if(tj.Pts[jpt].Chg > 0) ++nPtsWithCharge;
        if(jj > 0 && std::abs(tj.Pts[jpt+1].Pos[0] - tj.Pts[jpt].Pos[0]) < 1.5) ++nConsecutivePts;
        if(jpt == 0) break;
      } // jj
      if(prt) mf::LogVerbatim("TC")<<"   newEndPt "<<newEndPt<<" nConsecutivePts "<<nConsecutivePts<<" Required "<<minPts - 1;
    }
    
    float nwires = std::abs(tj.Pts[tj.EndPt[0]].Pos[0] - tj.Pts[newEndPt].Pos[0]) + 1;
    float npwc = NumPtsWithCharge(tjs, tj, true, tj.EndPt[0], newEndPt);
    hitFrac = npwc / nwires;
    
    if(hitFrac < fQualityCuts[0]) tj.AlgMod[kKilled] = true;
    if(prt) mf::LogVerbatim("TC")<<" Old endpoint "<<tj.EndPt[1]<<" newEndPt "<<newEndPt<<" nwires "<<nwires<<" npwc "<<npwc<<" nConsecutivePts "<<nConsecutivePts<<" hitFrac "<<hitFrac<<" Killed? "<<tj.AlgMod[kKilled];
    
    // failed the cuts
    if(tj.AlgMod[kKilled]) return;
    
    // modifications required
    tj.EndPt[1] = newEndPt;    
    for(unsigned short ipt = newEndPt + 1; ipt < tj.Pts.size(); ++ipt) {
      if(prt) mf::LogVerbatim("TC")<<" unset "<<ipt;
      UnsetUsedHits(tjs, tj.Pts[ipt]);
    }
    SetEndPoints(tjs, tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kTrimEndPts] = true;
    if(prt) PrintTrajectory("TEPo", tjs, tj, USHRT_MAX);
    
  } // TrimEndPts
  
  /////////////////////////////////////////
  bool SignalBetween(TjStuff& tjs, const TrajPoint& tp1, const TrajPoint& tp2, const float& MinWireSignalFraction, bool prt)
  {
    // Returns true if there is a signal on > MinWireSignalFraction of the wires between tp1 and tp2.
    if(MinWireSignalFraction == 0) return true;
    
    int fromWire = std::nearbyint(tp1.Pos[0]);
    int toWire = std::nearbyint(tp2.Pos[0]);
    
    if(fromWire == toWire) {
      if(prt) mf::LogVerbatim("TC")<<" SignalBetween fromWire = toWire = "<<fromWire<<" SignalAtTp? "<<SignalAtTp(tjs, tp1);
      return SignalAtTp(tjs, tp1);
    }

    // define a trajectory point located at tp1 that has a direction towards tp2
    TrajPoint tp;
    if(!MakeBareTrajPoint(tjs, tp1, tp2, tp)) return true;
    
    return SignalBetween(tjs, tp, toWire, MinWireSignalFraction, prt);

  } // SignalBetween

  /////////////////////////////////////////
  bool SignalBetween(TjStuff& tjs, TrajPoint tp, float toPos0, const float& MinWireSignalFraction, bool prt)
  {
    // Returns true if there is a signal on > MinWireSignalFraction of the wires between tp and toPos0.
    // Note that this uses the direction vector of the tp
    
    if(MinWireSignalFraction == 0) return true;
    
    int fromWire = std::nearbyint(tp.Pos[0]);
    int toWire = std::nearbyint(toPos0);
    
    if(fromWire == toWire) {
      if(prt) mf::LogVerbatim("TC")<<" SignalBetween fromWire = toWire = "<<fromWire<<" SignalAtTp? "<<SignalAtTp(tjs, tp);
      return SignalAtTp(tjs, tp);
    }
    
    int nWires = abs(toWire - fromWire) + 1;
    
    unsigned short maxWiresNoSignal = (1 - MinWireSignalFraction) * nWires;
    if(std::abs(tp.Dir[0]) < 0.001) tp.Dir[0] = 0.001;
    float stepSize = std::abs(1/tp.Dir[0]);
    // ensure that we step in the right direction
    if(toWire > fromWire && tp.Dir[0] < 0) stepSize = -stepSize;
    if(toWire < fromWire && tp.Dir[0] > 0) stepSize = -stepSize;
    unsigned short nsig = 0;
    unsigned short num = 0;
    unsigned short nmissed = 0;
    for(unsigned short cnt = 0; cnt < nWires; ++cnt) {
      ++num;
      if(SignalAtTp(tjs, tp)) {
        ++nsig;
      } else {
        ++nmissed;
        if(nmissed == maxWiresNoSignal) return false;
      }
      tp.Pos[0] += tp.Dir[0] * stepSize;
      tp.Pos[1] += tp.Dir[1] * stepSize;
    } // cnt
    float sigFrac = (float)nsig / (float)nWires;
    if(prt) mf::LogVerbatim("TC")<<"  SignalBetween fromWire "<<fromWire<<" toWire "<<toWire<<" nWires "<<nWires<<" nsig "<<nsig<<" "<<sigFrac;
    return (sigFrac >= MinWireSignalFraction);
    
  } // SignalBetween
  
  /////////////////////////////////////////
  bool SignalAtTp(TjStuff& tjs, const TrajPoint& tp)
  {
    return SignalAtPos(tjs, tp.Pos[0], tp.Pos[1], tp.CTP);
  } // SignalAtTp
  
  /////////////////////////////////////////
  bool SignalAtPos(TjStuff& tjs, const float& pos0, const float& pos1, CTP_t tCTP)
  {
    // Returns true if the TP is near the position
    
    if(pos0 < 0) return false;
    if(pos1 < 0) return false;
    unsigned int wire = std::nearbyint(pos0);
    geo::PlaneID planeID = DecodeCTP(tCTP);
    unsigned int ipl = planeID.Plane;
    if(wire >= tjs.NumWires[ipl]) return false;
    if(pos1 > tjs.MaxPos1[ipl]) return false;
    // Assume dead wires have a signal
    if(tjs.WireHitRange[ipl][wire].first == -1) return true;
    raw::TDCtick_t rawProjTick = (float)(pos1 / tjs.UnitsPerTick);
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      if(rawProjTick > tjs.fHits[iht].StartTick && rawProjTick < tjs.fHits[iht].EndTick) return true;
    } // iht
    return false;
  } // SignalAtPos

  //////////////////////////////////////////
  bool TjHasNiceVtx(TjStuff& tjs, const Trajectory& tj)
  {
    // returns true if there is a high-quality vertex at either end
    for(unsigned short end = 0; end < 2; ++end) {
      if(tj.VtxID[end] > 0) {
        unsigned short ivx = tj.VtxID[end] - 1;
        if(tjs.vtx[ivx].Stat[kNiceVtx]) return true;
      }
    } // end
    return false;
    
  } // TjHasNiceVtx
  
  //////////////////////////////////////////
  bool CheckHitClusterAssociations(TjStuff& tjs)
  {
    // check hit - cluster associations
    
    if(tjs.fHits.size() != tjs.inClus.size()) {
      mf::LogWarning("TC")<<"CHCA: Sizes wrong "<<tjs.fHits.size()<<" "<<tjs.inClus.size();
      return false;
    }
    
    unsigned int iht;
    short clID;
    
    // check cluster -> hit association
    for(unsigned short icl = 0; icl < tjs.tcl.size(); ++icl) {
      if(tjs.tcl[icl].ID < 0) continue;
      clID = tjs.tcl[icl].ID;
      for(unsigned short ii = 0; ii < tjs.tcl[icl].tclhits.size(); ++ii) {
        iht = tjs.tcl[icl].tclhits[ii];
        if(iht > tjs.fHits.size() - 1) {
          mf::LogWarning("CC")<<"CHCA: Bad tclhits index "<<iht<<" tjs.fHits size "<<tjs.fHits.size();
          return false;
        } // iht > tjs.fHits.size() - 1
        if(tjs.inClus[iht] != clID) {
          mf::LogError("TC")<<"CHCA: Bad cluster -> hit association. clID "<<clID<<" hit "<<PrintHit(tjs.fHits[iht])<<" tjs.inClus "<<tjs.inClus[iht]<<" CTP "<<tjs.tcl[icl].CTP;
          return false;
        }
      } // ii
    } // icl
    
    // check hit -> cluster association
    unsigned short icl;
    for(iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.inClus[iht] <= 0) continue;
      icl = tjs.inClus[iht] - 1;
      // see if the cluster is obsolete
      if(tjs.tcl[icl].ID < 0) {
        mf::LogError("TC")<<"CHCA: Hit "<<PrintHit(tjs.fHits[iht])<<" associated with an obsolete cluster tjs.tcl[icl].ID "<<tjs.tcl[icl].ID;
        return false;
      }
      if (std::find(tjs.tcl[icl].tclhits.begin(), tjs.tcl[icl].tclhits.end(), iht) == tjs.tcl[icl].tclhits.end()) {
        mf::LogError("TC")<<"CHCA: Hit "<<":"<<PrintHit(tjs.fHits[iht])<<" -> tjs.inClus "<<tjs.inClus[iht]<<" but isn't in tjs.tcl[icl].ID "<<tjs.tcl[icl].ID<<" list of hits. icl "<<icl<<" iht "<<iht;
        for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
          if(tjs.allTraj[itj].ClusterIndex == icl) mf::LogError("TC")<<"CHCA: Cluster index "<<icl<<" found in traj ID "<<tjs.allTraj[itj].ID;
        } // itj
        PrintAllTraj("CHCA", tjs, debug, USHRT_MAX, USHRT_MAX);
        return false;
      }
    } // iht
    
    return true;
    
  } // CheckHitClusterAssociations()
  
  //////////////////////////////////////////
  unsigned short NumPtsWithCharge(TjStuff& tjs, const Trajectory& tj, bool includeDeadWires)
  {
    unsigned short firstPt = tj.EndPt[0];
    unsigned short lastPt = tj.EndPt[1];
    return NumPtsWithCharge(tjs, tj, includeDeadWires, firstPt, lastPt);
  }
  
  //////////////////////////////////////////
  unsigned short NumPtsWithCharge(TjStuff& tjs, const Trajectory& tj, bool includeDeadWires, unsigned short firstPt, unsigned short lastPt)
  {
    unsigned short ntp = 0;
    for(unsigned short ipt = firstPt; ipt <= lastPt; ++ipt) if(tj.Pts[ipt].Chg > 0) ++ntp;
    // Add the count of deadwires
    if(includeDeadWires) ntp += DeadWireCount(tjs, tj.Pts[firstPt], tj.Pts[lastPt]);
    return ntp;
  } // NumPtsWithCharge
  
  //////////////////////////////////////////
  float DeadWireCount(TjStuff& tjs, const TrajPoint& tp1, const TrajPoint& tp2)
  {
    return DeadWireCount(tjs, tp1.Pos[0], tp2.Pos[0], tp1.CTP);
  } // DeadWireCount
  
  //////////////////////////////////////////
  float DeadWireCount(TjStuff& tjs, const float& inWirePos1, const float& inWirePos2, CTP_t tCTP)
  {
    if(inWirePos1 < -0.4 || inWirePos2 < -0.4) return 0;
    unsigned int inWire1 = std::nearbyint(inWirePos1);
    unsigned int inWire2 = std::nearbyint(inWirePos2);
    geo::PlaneID planeID = DecodeCTP(tCTP);
    unsigned short plane = planeID.Plane;
    if(inWire1 > tjs.NumWires[plane] || inWire2 > tjs.NumWires[plane]) return 0;
    if(inWire1 > inWire2) {
      // put in increasing order
      unsigned int tmp = inWire1;
      inWire1 = inWire2;
      inWire2 = tmp;
    } // inWire1 > inWire2
    ++inWire2;
    unsigned int wire, ndead = 0;
    for(wire = inWire1; wire < inWire2; ++wire) if(tjs.WireHitRange[plane][wire].first == -1) ++ndead;
    return ndead;
  } // DeadWireCount

  ////////////////////////////////////////////////
  unsigned short PDGCodeIndex(TjStuff& tjs, int PDGCode)
  {
    unsigned short pdg = abs(PDGCode);
    if(pdg == 11) return 0; // electron
    if(pdg == 13) return 1; // muon
    if(pdg == 211) return 2; // pion
    if(pdg == 321) return 3; // kaon
    if(pdg == 2212) return 4; // proton
    
    return USHRT_MAX;
    
  } // PDGCodeIndex
  
  ////////////////////////////////////////////////
  bool WireHitRangeOK(const TjStuff& tjs, const CTP_t& inCTP)
  {
    // returns true if the passed CTP code is consistent with the CT code of the WireHitRangeVector
    geo::PlaneID planeID = DecodeCTP(inCTP);
    if(planeID.Cryostat != tjs.WireHitRangeCstat) return false;
    if(planeID.TPC != tjs.WireHitRangeTPC) return false;
    return true;
  }
  
  ////////////////////////////////////////////////
  void MakeVertexObsolete(TjStuff& tjs, unsigned short vtxID)
  {
    // deletes a 2D vertex and possibly a 3D vertex and 2D vertices in other planes
    if(vtxID > tjs.vtx.size()) return;
    unsigned short ivx = vtxID - 1;
    tjs.vtx[ivx].NTraj = 0;
    for(auto& tj : tjs.allTraj) {
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] == vtxID) tj.VtxID[end] = 0;
      } // end
    } // tj
    // check for a 3D vertex
    if(tjs.vtx[ivx].Ptr3D == SHRT_MAX) return;
    unsigned short ivx3 = tjs.vtx[ivx].Ptr3D;
    // make this obsolete
    tjs.vtx3[ivx3].Wire = SHRT_MAX;
    // look for the matched 2D vertices
    for(unsigned short ipl = 0; ipl < 3; ++ipl) {
      // no 2D vertex in this plane
      if(tjs.vtx3[ivx3].Ptr2D[ipl] < 0) continue;
      // already killed 2D vertex in this plane
      if(tjs.vtx3[ivx3].Ptr2D[ipl] > (short)(tjs.vtx.size() - 1)) continue;
      unsigned short ivx2 = tjs.vtx3[ivx3].Ptr2D[ipl];
      unsigned short vtx2ID = ivx2 + 1;
      tjs.vtx[ivx2].NTraj = 0;
      for(auto& tj : tjs.allTraj) {
        for(unsigned short end = 0; end < 2; ++end) {
          if(tj.VtxID[end] == vtx2ID) tj.VtxID[end] = 0;
        } // end
      } // tj
    } // ipl
    
  } // MakeVertexObsolete

  ////////////////////////////////////////////////
  void MakeTrajectoryObsolete(TjStuff& tjs, unsigned short itj)
  {
    // Note that this does not change the state of UseHit to allow
    // resurrecting the trajectory later (RestoreObsoleteTrajectory)
    if(itj > tjs.allTraj.size() - 1) return;
    unsigned int iht;
    for(auto& tp : tjs.allTraj[itj].Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        iht = tp.Hits[ii];
        if(tjs.fHits[iht].InTraj == tjs.allTraj[itj].ID) tjs.fHits[iht].InTraj = 0;
      } // ii
    } // tp
    tjs.allTraj[itj].AlgMod[kKilled] = true;
  } // MakeTrajectoryObsolete
  
  ////////////////////////////////////////////////
  void RestoreObsoleteTrajectory(TjStuff& tjs, unsigned short itj)
  {
    if(itj > tjs.allTraj.size() - 1) return;
    if(!tjs.allTraj[itj].AlgMod[kKilled]) {
      mf::LogWarning("TC")<<"RestoreObsoleteTrajectory: Trying to restore not-obsolete trajectory "<<tjs.allTraj[itj].ID;
      return;
    }
    unsigned int iht;
    for(auto& tp : tjs.allTraj[itj].Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) {
          iht = tp.Hits[ii];
          if(tjs.fHits[iht].InTraj == 0) {
            tjs.fHits[iht].InTraj = tjs.allTraj[itj].ID;
          }
        }
      } // ii
    } // tp
    tjs.allTraj[itj].AlgMod[kKilled] = false;
  } // RestoreObsoleteTrajectory
  
  //////////////////////////////////////////
  void CheckVtxAssociations(TjStuff& tjs, const CTP_t& inCTP)
  {
    // start by setting NTraj = 0 for all vertices
    for(auto& vtx : tjs.vtx) {
      if(vtx.CTP == inCTP) {
        vtx.NTraj = 0;
        vtx.Stat[kNiceVtx] = false;
      } // CTP check
    } // vtx
    for(auto& tj : tjs.allTraj) {
      if(tj.CTP != inCTP) continue;
      if(tj.AlgMod[kKilled]) continue;
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] > 0) {
          unsigned short ivx = tj.VtxID[end] - 1;
          ++tjs.vtx[ivx].NTraj;
          // set the NiceVtx bit?
          bool niceVtx = (tj.PDGCode == 13);
          if(NumPtsWithCharge(tjs, tj, false) > 20 && tj.MCSMom > 100) niceVtx = true;
          if(niceVtx) tjs.vtx[ivx].Stat[kNiceVtx] = true;
          if(tjs.vtx[ivx].CTP != inCTP) std::cout<<"CheckVertexAssociations: CTP tj-vtx mis-match\n";
        }
      } // end
    } // tj
  } // CheckVtxAssociations

  //////////////////////////////////////////
  bool SplitAllTraj(TjStuff& tjs, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt)
  {
    // Splits the trajectory itj in the tjs.allTraj vector into two trajectories at position pos. Splits
    // the trajectory and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
    
    if(itj > tjs.allTraj.size()-1) return false;
    if(pos < tjs.allTraj[itj].EndPt[0] + 1 || pos > tjs.allTraj[itj].EndPt[1] - 1) return false;
    
    Trajectory& tj = tjs.allTraj[itj];
    
    // ensure that there will be at least 3 TPs on each trajectory
    unsigned short ipt, ii, ntp = 0;
    for(ipt = 0; ipt < pos; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) {
      if(prt) mf::LogVerbatim("TC")<<"SplitAllTraj: Split point to small at begin "<<ntp<<" pos "<<pos<<" ID ";
      return false;
    }
    ntp = 0;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) {
      if(prt) mf::LogVerbatim("TC")<<"SplitAllTraj: Split point too small at end "<<ntp<<" pos "<<pos<<" EndPt "<<tj.EndPt[1];
      return false;
    }
    
    // make a copy
    Trajectory newTj = tjs.allTraj[itj];
    newTj.ID = tjs.allTraj.size() + 1;
    
    // Leave the first section of tj in place. Re-assign the hits
    // to the new trajectory
    unsigned int iht;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(!tj.Pts[ipt].UseHit[ii]) continue;
        iht = tj.Pts[ipt].Hits[ii];
        // This shouldn't happen but check anyway
        if(tjs.fHits[iht].InTraj != tj.ID) continue;
        tjs.fHits[iht].InTraj = newTj.ID;
        tj.Pts[ipt].UseHit[ii] = false;
      } // ii
    } // ipt
    SetEndPoints(tjs, tj);
    if(ivx != USHRT_MAX) tj.VtxID[1] = tjs.vtx[ivx].ID;
    tj.AlgMod[kSplitTraj] = true;
    if(prt) {
      mf::LogVerbatim("TC")<<"Splitting trajectory ID "<<tj.ID<<" new EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
    }
    
    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if(eraseSize > newTj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"SplitAllTraj: Bad erase size ";
      return false;
    }
    
    // erase the TPs at the beginning of the new trajectory
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for(ipt = 0; ipt < 3; ++ipt) {
      for(ii = 0; ii < newTj.Pts[ipt].Hits.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
      newTj.Pts[ipt].Chg = 0;
    } // ipt
    SetEndPoints(tjs, newTj);
    if(ivx != USHRT_MAX) newTj.VtxID[0] = tjs.vtx[ivx].ID;
    newTj.AlgMod[kSplitTraj] = true;
    tjs.allTraj.push_back(newTj);
    if(prt) {
      mf::LogVerbatim("TC")<<"  newTj ID "<<newTj.ID<<" EndPts "<<newTj.EndPt[0]<<" to "<<newTj.EndPt[1];
    }
    return true;
    
  } // SplitAllTraj
  
  //////////////////////////////////////////
  void TrajPointTrajDOCA(TjStuff& tjs, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep)
  {
    // Finds the point, ipt, on trajectory tj that is closest to trajpoint tp
    float best = minSep * minSep;
    closePt = USHRT_MAX;
    float dw, dt, dp2;
    unsigned short ipt;
    for(ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      dw = tj.Pts[ipt].Pos[0] - tp.Pos[0];
      dt = tj.Pts[ipt].Pos[1] - tp.Pos[1];
      dp2 = dw * dw + dt * dt;
      if(dp2 < best) {
        best = dp2;
        closePt = ipt;
      }
    } // ipt
    minSep = sqrt(best);
  } // TrajPointTrajDOCA
  
  //////////////////////////////////////////
  void TrajTrajDOCA(TjStuff& tjs, Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
  {
    TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, minSep, false);
  } // TrajTrajDOCA
  
  //////////////////////////////////////////
  void TrajTrajDOCA(TjStuff& tjs, Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep, bool considerDeadWires)
  {
    // Find the Distance Of Closest Approach between two trajectories less than minSep
    float best = minSep * minSep;
    ipt1 = 0; ipt2 = 0;
    float dwc = 0;
    for(unsigned short i1 = tj1.EndPt[0]; i1 < tj1.EndPt[1] + 1; ++i1) {
      for(unsigned short i2 = tj2.EndPt[0]; i2 < tj2.EndPt[1] + 1; ++i2) {
        if(considerDeadWires) dwc = DeadWireCount(tjs, tj1.Pts[i1], tj2.Pts[i2]);
        float dw = tj1.Pts[i1].Pos[0] - tj2.Pts[i2].Pos[0] - dwc;
        if(std::abs(dw) > minSep) continue;
        float dt = tj1.Pts[i1].Pos[1] - tj2.Pts[i2].Pos[1];
        if(std::abs(dt) > minSep) continue;
        float dp2 = dw * dw + dt * dt;
        if(dp2 < best) {
          best = dp2;
          ipt1 = i1;
          ipt2 = i2;
        }
      } // i2
    } // i1
    minSep = sqrt(best);
  } // TrajTrajDOCA

  //////////////////////////////////////////
  float HitSep2(TjStuff& tjs, unsigned int iht, unsigned int jht)
  {
    // returns the separation^2 between two hits in WSE units
    if(iht > tjs.fHits.size()-1 || jht > tjs.fHits.size()-1) return 1E6;
    float dw = (float)tjs.fHits[iht].WireID.Wire - (float)tjs.fHits[jht].WireID.Wire;
    float dt = (tjs.fHits[iht].PeakTime - tjs.fHits[jht].PeakTime) * tjs.UnitsPerTick;
    return dw * dw + dt * dt;
  } // HitSep2
  
  //////////////////////////////////////////
  float PointTrajSep2(float wire, float time, TrajPoint const& tp)
  {
    float dw = wire - tp.Pos[0];
    float dt = time - tp.Pos[1];
    return dw * dw + dt * dt;
  }
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff const& tjs, unsigned int iht, TrajPoint const& tp)
  {
    float wire = tjs.fHits[iht].WireID.Wire;
    float time = tjs.fHits[iht].PeakTime * tjs.UnitsPerTick;
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff const& tjs, float wire, float time, TrajPoint const& tp)
  {
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA2(TjStuff const& tjs, float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach squared between a (wire, time(WSE)) point
    // and a trajectory point
    
    float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
    float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    float dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return (dw * dw + dt * dt);
    
  } // PointTrajDOCA2
  
  //////////////////////////////////////////
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, std::array<float, 2>& pos)
  {
    TrajIntersection(tp1, tp2, pos[0], pos[1]);
  } // TrajIntersection
  //////////////////////////////////////////
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y)
  {
    // returns the intersection position, (x,y), of two trajectory points
    
    x = -9999; y = -9999;
    
    double arg1 = tp1.Pos[0] * tp1.Dir[1] - tp1.Pos[1] * tp1.Dir[0];
    double arg2 = tp2.Pos[0] * tp1.Dir[1] - tp2.Pos[1] * tp1.Dir[0];
    double arg3 = tp2.Dir[0] * tp1.Dir[1] - tp2.Dir[1] * tp1.Dir[0];
    if(arg3 == 0) return;
    double s = (arg1 - arg2) / arg3;
    
    x = (float)(tp2.Pos[0] + s * tp2.Dir[0]);
    y = (float)(tp2.Pos[1] + s * tp2.Dir[1]);
    
  } // TrajIntersection
  
  //////////////////////////////////////////
  float TrajLength(Trajectory& tj)
  {
    float len = 0, dx, dy;
    unsigned short ipt;
    unsigned short prevPt = tj.EndPt[0];
    for(ipt = tj.EndPt[0] + 1; ipt < tj.EndPt[1] + 1; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dx = tj.Pts[ipt].Pos[0] - tj.Pts[prevPt].Pos[0];
      dy = tj.Pts[ipt].Pos[1] - tj.Pts[prevPt].Pos[1];
      len += sqrt(dx * dx + dy * dy);
      prevPt = ipt;
    }
    return len;
  } // TrajLength

  //////////////////////////////////////////
  float PosSep(const std::array<float, 2>& pos1, const std::array<float, 2>& pos2)
  {
    return sqrt(PosSep2(pos1, pos2));
  } // PosSep
  
  //////////////////////////////////////////
  float PosSep2(const std::array<float, 2>& pos1, const std::array<float, 2>& pos2)
  {
    // returns the separation distance^2 between two positions
    float d0 = pos1[0] - pos2[0];
    float d1 = pos1[1] - pos2[1];
    return d0*d0+d1*d1;
  } // PosSep2
  
  //////////////////////////////////////////
  float PosSep2(const std::array<float, 3>& pos1, const std::array<float, 3>& pos2)
  {
    // returns the separation distance^2 between two positions in 3D
    float d0 = pos1[0] - pos2[0];
    float d1 = pos1[1] - pos2[1];
    float d2 = pos1[2] - pos2[2];
    return d0*d0 + d1*d1 + d2*d2;
  } // PosSep2
  
  //////////////////////////////////////////
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Returns the separation distance between two trajectory points
    float dx = tp1.Pos[0] - tp2.Pos[0];
    float dy = tp1.Pos[1] - tp2.Pos[1];
    return sqrt(dx * dx + dy * dy);
  } // TrajPointSeparation
  
  //////////////////////////////////////////
  void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& Distance)
  {
    // find the closest approach between a trajectory tj and a point (x,y). Returns
    // the index of the closest trajectory point and the distance
    
    float dx, dy, dist, best = 1E6;
    closePt = 0;
    Distance = best;
    
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dx = tj.Pts[ipt].Pos[0] - x;
      dy = tj.Pts[ipt].Pos[1] - y;
      dist = dx * dx + dy * dy;
      if(dist < best) {
        best = dist;
        closePt = ipt;
      }
      // TODO is this wise?
      //      if(dist > best) break;
    } // ipt
    
    Distance = sqrt(best);
    
  } // TrajClosestApproach
  
  /////////////////////////////////////////
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2)
  {
    // Calculates the angle of a line between two TPs
    float dw = tp2.Pos[0] - tp1.Pos[0];
    float dt = tp2.Pos[1] - tp1.Pos[1];
    return atan2(dw, dt);
  } // TwoTPAngle
  
  ////////////////////////////////////////////////
  std::vector<unsigned int> PutTrajHitsInVector(Trajectory const& tj, HitStatus_t hitRequest)
  {
    // Put hits in each trajectory point into a flat vector
    std::vector<unsigned int> hitVec;
    
    // special handling for shower trajectories
    if(tj.AlgMod[kShowerTj]) return tj.Pts[1].Hits;
    
    // reserve under the assumption that there will be one hit per point
    hitVec.reserve(tj.Pts.size());
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        unsigned int iht = tj.Pts[ipt].Hits[ii];
        bool useit = (hitRequest == kAllHits);
        if(tj.Pts[ipt].UseHit[ii] && hitRequest == kUsedHits) useit = true;
        if(!tj.Pts[ipt].UseHit[ii] && hitRequest == kUnusedHits) useit = true;
        if(useit) hitVec.push_back(iht);
      } // iht
    } // ipt
    return hitVec;
  } // PutTrajHitsInVector
  
  //////////////////////////////////////////
  bool HitIsInTj(Trajectory const& tj, const unsigned int& iht, short nPtsToCheck)
  {
    // returns true if hit iht is associated with trajectory tj. Checking starts at the
    // end of tj for nPtsToCheck points
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      if(std::find(tj.Pts[ipt].Hits.begin(), tj.Pts[ipt].Hits.end(), iht) != tj.Pts[ipt].Hits.end()) return true;
      // only go back a few points
      if(nPtsToCheck >= 0 && ii == nPtsToCheck) return false;
      if(ipt == 0) return false;
    } // ii
    return false;

  } // HitIsInTj
  
  //////////////////////////////////////////
  bool HasDuplicateHits(TjStuff const& tjs, Trajectory const& tj, bool prt)
  {
    // returns true if a hit is associated with more than one TP
    auto tjHits = PutTrajHitsInVector(tj, kAllHits);
    for(unsigned short ii = 0; ii < tjHits.size() - 1; ++ii) {
      for(unsigned short jj = ii + 1; jj < tjHits.size(); ++jj) {
        if(tjHits[ii] == tjHits[jj]) {
          if(prt) mf::LogVerbatim()<<"HDH: Hit "<<PrintHit(tjs.fHits[ii])<<" is a duplicate "<<ii<<" "<<jj;
          return true;
        }
      } // jj
    } // ii
    return false;
  } // HasDuplicateHits
  
  //////////////////////////////////////////
  void MoveTPToWire(TrajPoint& tp, float wire)
  {
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    if(tp.Dir[0] == 0) return;
    float dw = wire - tp.Pos[0];
    if(std::abs(dw) < 0.01) return;
    tp.Pos[0] = wire;
    tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
  } // MoveTPToWire
  
  //////////////////////////////////////////
  std::vector<unsigned int> FindCloseHits(TjStuff const& tjs, std::array<int, 2> const& wireWindow, std::array<float, 2> const& timeWindow, const unsigned short plane, HitStatus_t hitRequest, bool usePeakTime, bool& hitsNear)
  {
    // returns a vector of hits that are within the Window[Pos0][Pos1] in plane.
    // Note that hits on wire wireWindow[1] are returned as well. The definition of close
    // depends on setting of usePeakTime. If UsePeakTime is true, a hit is considered nearby if
    // the PeakTime is within the window. This is shown schematically here where
    // the time is on the horizontal axis and a "-" denotes a valid entry
    // timeWindow     -----------------
    // hit PeakTime             +         close
    // hit PeakTime  +                    not close
    // If usePeakTime is false, a hit is considered nearby if the hit StartTick and EndTick overlap with the timeWindow
    // Time window                  ---------
    // Hit StartTick-EndTick      --------        close
    // Hit StartTick - EndTick                  --------  not close
    
    hitsNear = false;
    std::vector<unsigned int> closeHits;
    if(plane > tjs.FirstWire.size() - 1) return closeHits;
    // window in the wire coordinate
    int loWire = wireWindow[0];
    if(loWire < (int)tjs.FirstWire[plane]) loWire = tjs.FirstWire[plane];
    int hiWire = wireWindow[1];
    if(hiWire > (int)tjs.LastWire[plane]-1) hiWire = tjs.LastWire[plane]-1;
    // window in the time coordinate
    float minTick = timeWindow[0] / tjs.UnitsPerTick;
    float maxTick = timeWindow[1] / tjs.UnitsPerTick;
    for(int wire = loWire; wire <= hiWire; ++wire) {
      // Set hitsNear if the wire is dead
      if(tjs.WireHitRange[plane][wire].first == -2) hitsNear = true;
      if(tjs.WireHitRange[plane][wire].first < 0) continue;
      unsigned int firstHit = (unsigned int)tjs.WireHitRange[plane][wire].first;
      unsigned int lastHit = (unsigned int)tjs.WireHitRange[plane][wire].second;
      for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
        if(usePeakTime) {
          if(tjs.fHits[iht].PeakTime < minTick) continue;
          if(tjs.fHits[iht].PeakTime > maxTick) break;
        } else {
          int hiLo = minTick;
          if(tjs.fHits[iht].StartTick > hiLo) hiLo = tjs.fHits[iht].StartTick;
          int loHi = maxTick;
          if(tjs.fHits[iht].EndTick < loHi) loHi = tjs.fHits[iht].EndTick;
          if(loHi < hiLo) continue;
          if(hiLo > loHi) break;
        }
        hitsNear = true;
        bool takeit = (hitRequest == kAllHits);
        if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) takeit = true;
        if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) takeit = true;
        if(takeit) closeHits.push_back(iht);
      } // iht
    } // wire
    return closeHits;
  } // FindCloseHits

  //////////////////////////////////////////
  bool FindCloseHits(TjStuff const& tjs, TrajPoint& tp, float const& maxDelta, HitStatus_t hitRequest)
  {
    // Fills tp.Hits sets tp.UseHit true for hits that are close to tp.Pos. Returns true if there are
    // close hits OR if the wire at this position is dead
    
    tp.Hits.clear();
    tp.UseHit.reset();
    if(!WireHitRangeOK(tjs, tp.CTP)) {
      std::cout<<"FindCloseHits: WireHitRange not valid for CTP "<<tp.CTP<<". tjs.WireHitRange Cstat "<<tjs.WireHitRangeCstat<<" TPC "<<tjs.WireHitRangeTPC<<"\n";
      return false;
    }
    
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short ipl = planeID.Plane;
    
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    if(wire < tjs.FirstWire[ipl]) return false;
    if(wire > tjs.LastWire[ipl]-1) return false;
    
    // dead wire
    if(tjs.WireHitRange[ipl][wire].first == -1) return true;
    // live wire with no hits
    if(tjs.WireHitRange[ipl][wire].first == -2) return false;
    
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;

    float fwire = wire;
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) useit = true;
      if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) useit = true;
      if(!useit) continue;
      float ftime = tjs.UnitsPerTick * tjs.fHits[iht].PeakTime;
      float delta = PointTrajDOCA(tjs, fwire, ftime, tp);
//      std::cout<<"chk "<<PrintHit(tjs.fHits[iht])<<" delta "<<delta<<" maxDelta "<<maxDelta<<"\n";
      if(delta < maxDelta) tp.Hits.push_back(iht);
    } // iht
    if(tp.Hits.size() > 16) {
//      mf::LogWarning("TC")<<"FindCloseHits: Found "<<tp.Hits.size()<<" hits. Truncating to 16";
      tp.Hits.resize(16);
    }
    // Set UseHit false. The calling routine should decide if these hits should be used
    tp.UseHit.reset();
    return true;
    
  } // FindCloseHits
  
  ////////////////////////////////////////////////
  float MaxHitDelta(TjStuff& tjs, Trajectory& tj)
  {
    float delta, md = 0;
    unsigned short ii;
    unsigned int iht;
    for(auto& tp : tj.Pts) {
      for(ii = 0; ii < tp.Hits.size(); ++ii) {
        if(!tp.UseHit[ii]) continue;
        iht = tp.Hits[ii];
        delta = PointTrajDOCA(tjs, iht, tp);
        if(delta > md) md = delta;
      } // ii
    } // pts
    return md;
  } // MaxHitDelta

  //////////////////////////////////////////
  void ReverseTraj(TjStuff& tjs, Trajectory& tj)
  {
    // reverse the trajectory
    if(tj.Pts.empty()) return;
    // reverse the crawling direction flag
    tj.StepDir = -tj.StepDir;
    // reverse the direction
    tj.TjDir = -tj.TjDir;
    // Vertices
    std::swap(tj.VtxID[0], tj.VtxID[1]);
    // trajectory points
    std::reverse(tj.Pts.begin(), tj.Pts.end());
    // reverse the stop flag
    std::reverse(tj.StopFlag.begin(), tj.StopFlag.end());
    // reverse the direction vector on all points
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Dir[0] != 0) tj.Pts[ipt].Dir[0] = -tj.Pts[ipt].Dir[0];
      if(tj.Pts[ipt].Dir[1] != 0) tj.Pts[ipt].Dir[1] = -tj.Pts[ipt].Dir[1];
      tj.Pts[ipt].Ang = std::atan2(tj.Pts[ipt].Dir[1], tj.Pts[ipt].Dir[0]);
    } // ipt
    SetEndPoints(tjs, tj);
  } // ReverseTraj
  
  //////////////////////////////////////////
  bool PointInsideEnvelope(const std::array<float, 2>& Point, const std::vector<std::array<float, 2>>& Envelope)
  {
    // returns true if the Point is within the Envelope polygon. Entries in Envelope are the
    // Pos[0], Pos[1] locations of the polygon vertices. This is based on the algorithm that the
    // sum of the angles of a vector between a point and the vertices will be 2 * pi for an interior
    // point and 0 for an exterior point
    
    std::array<float, 2> p1, p2;
    unsigned short nvx = Envelope.size();
    double angleSum = 0;
    for(unsigned short ii = 0; ii < Envelope.size(); ++ii) {
      p1[0] = Envelope[ii][0] - Point[0];
      p1[1] = Envelope[ii][1] - Point[1];
      p2[0] = Envelope[(ii+1)%nvx][0] - Point[0];
      p2[1] = Envelope[(ii+1)%nvx][1] - Point[1];
      angleSum += DeltaAngle(p1, p2);
    }
    if(abs(angleSum) < M_PI) return false;
    return true;
      
  } // InsideEnvelope

  //////////////////////////////////////////
  double DeltaAngle(const std::array<float,2>& p1, const std::array<float,2>& p2)
  {
    // angle between two points
    double ang1 = atan2(p1[1], p1[0]);
    double ang2 = atan2(p2[1], p2[0]);
    return DeltaAngle2(ang1, ang2);
  } // DeltaAngle
  
  //////////////////////////////////////////
  double DeltaAngle2(double Ang1, double Ang2)
  {
    constexpr double twopi = 2 * M_PI;
    double dang = Ang1 - Ang2;
    while(dang >  M_PI) dang -= twopi;
    while(dang < -M_PI) dang += twopi;
    return dang;
  }

  //////////////////////////////////////////
  double DeltaAngle(double Ang1, double Ang2) 
  {
    return std::abs(std::remainder(Ang1 - Ang2, M_PI));
  }
  
  ////////////////////////////////////////////////
  void SetEndPoints(TjStuff& tjs, Trajectory& tj)
  {
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
    
    tj.EndPt[0] = 0; tj.EndPt[1] = 0;
    if(tj.Pts.size() == 0) return;
    
    // check the end point pointers
    for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[0] = ipt;
        break;
      }
    }
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[1] = ipt;
        break;
      }
    }
  } // SetEndPoints
  
  ////////////////////////////////////////////////
  bool TrajIsClean(TjStuff& tjs, Trajectory& tj, bool prt)
  {
    // Returns true if the trajectory has low hit multiplicity and is in a
    // clean environment
    unsigned short nUsed = 0;
    unsigned short nTotHits = 0;
    for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      nTotHits += tp.Hits.size();
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        if(tp.UseHit[ii]) ++nUsed;
      } // ii
    } // ipt
    if(nTotHits == 0) return false;
    float fracUsed = (float)nUsed / (float)nTotHits;
    if(prt) mf::LogVerbatim("TC")<<"TrajIsClean: nTotHits "<<nTotHits<<" nUsed "<<nUsed<<" fracUsed "<<fracUsed;
    
    if(fracUsed > 0.9) return true;
    return false;
    
  } // TrajIsClean
  
  ////////////////////////////////////////////////
  short MCSMom(TjStuff& tjs, Trajectory& tj)
  {
    return MCSMom(tjs, tj, tj.EndPt[0], tj.EndPt[1]);
  } // MCSMom
  
  
  ////////////////////////////////////////////////
  short MCSMom(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt)
  {
    // Estimate the trajectory momentum using Multiple Coulomb Scattering ala PDG RPP
    
    firstPt = NearestPtWithChg(tjs, tj, firstPt);
    lastPt = NearestPtWithChg(tjs, tj, lastPt);
    if(firstPt >= lastPt) return 0;
    
    if(firstPt < tj.EndPt[0]) return 0;
    if(lastPt > tj.EndPt[1]) return 0;
    // Can't do this with only 2 points
    if(NumPtsWithCharge(tjs, tj, false, firstPt, lastPt) < 3) return 0;
        
    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if(tjLen == 0) return 0;
    // mom calculated in MeV
    double mom = 13.8 * sqrt(tjLen / 14) / MCSThetaRMS(tjs, tj, firstPt, lastPt);
    if(mom > 999) mom = 999;
    return (short)mom;
  } // MCSMom
    
  
  ////////////////////////////////////////////////
  unsigned short NearestPtWithChg(TjStuff& tjs, Trajectory& tj, unsigned short thePt)
  {
    // returns a point near thePt which has charge
    if(thePt > tj.EndPt[1]) return thePt;
    if(tj.Pts[thePt].Chg > 0) return thePt;
    
    short endPt0 = tj.EndPt[0];
    short endPt1 = tj.EndPt[1];
    for(short off = 1; off < 10; ++off) {
      short ipt = thePt + off;
      if(ipt <= endPt1 && tj.Pts[ipt].Chg > 0) return (unsigned short)ipt;
      ipt = thePt - off;
      if(ipt >= endPt0 && tj.Pts[ipt].Chg > 0) return (unsigned short)ipt;
    } // off
    return thePt;
  } // NearestPtWithChg
  
  /////////////////////////////////////////
  float MCSThetaRMS(TjStuff& tjs, Trajectory& tj)
  {
    // This returns the MCS scattering angle expected for one WSE unit of travel along the trajectory.
    // It is used to define kink and vertex cuts. This should probably be named something different to
    // prevent confusion
    
    return MCSThetaRMS(tjs, tj, tj.EndPt[0], tj.EndPt[1]) / sqrt(TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]));
    
  } // MCSThetaRMS
  
  /////////////////////////////////////////
  double MCSThetaRMS(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt)
  {
    // This returns the MCS scattering angle expected for the length of the trajectory
    // spanned by firstPt to lastPt. It is used primarily to calculate MCSMom
    
    if(firstPt < tj.EndPt[0]) return 1;
    if(lastPt > tj.EndPt[1]) return 1;
    
    firstPt = NearestPtWithChg(tjs, tj, firstPt);
    lastPt = NearestPtWithChg(tjs, tj, lastPt);
    if(firstPt >= lastPt) return 1;
    
    TrajPoint tmp;
    // make a bare trajectory point to define a line between firstPt and lastPt.
    // Use the position of the hits at these points
    TrajPoint firstTP = tj.Pts[firstPt];
    firstTP.Pos = firstTP.HitPos;
    TrajPoint lastTP = tj.Pts[lastPt];
    lastTP.Pos = lastTP.HitPos;
    if(!MakeBareTrajPoint(tjs, firstTP, lastTP, tmp)) return 1;
    // sum up the deviations^2
    double dsum = 0;
    unsigned short cnt = 0;
    for(unsigned short ipt = firstPt + 1; ipt < lastPt; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dsum += PointTrajDOCA2(tjs, tj.Pts[ipt].HitPos[0],  tj.Pts[ipt].HitPos[1], tmp);
      ++cnt;
    } // ipt
    if(cnt == 0) return 1;
    // require that cnt is a significant fraction of the total number of charged points
    // so that we don't get erroneously high MCSMom when there are large gaps.
    // This is the number of points expected in the count if there are no gaps
    unsigned short numPts = lastPt - firstPt - 1;
    // return the previously calculated value of MCSMom
    if(numPts > 5 && cnt < 0.7 * numPts) return tj.MCSMom;
    double sigmaS = sqrt(dsum / (double)cnt);
    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if(tjLen == 0) return 0;
    // Theta_o =  4 * sqrt(3) * sigmaS / path
    return (6.8 * sigmaS / tjLen);
    
  } // MCSThetaRMS

  /////////////////////////////////////////
  void TagDeltaRays(TjStuff& tjs, const CTP_t& inCTP, const std::vector<short>& fDeltaRayTag, short debugWorkID)
  {
    // fDeltaRayTag vector elements
    // [0] = max separation of both endpoints from a muon
    // [1] = minimum MCSMom
    // [2] = maximum MCSMom
    
    if(fDeltaRayTag[0] < 0) return;
    if(fDeltaRayTag.size() < 3) return;
    
    float sepCut = fDeltaRayTag[0];
    unsigned short minMom = fDeltaRayTag[1];
    unsigned short maxMom = fDeltaRayTag[2];
    
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& muTj = tjs.allTraj[itj];
      if(muTj.CTP != inCTP) continue;
      if(muTj.AlgMod[kKilled]) continue;
      bool prt = (muTj.WorkID == debugWorkID);
      if(prt) mf::LogVerbatim("TC")<<"TagDeltaRays: Muon "<<muTj.CTP<<" "<<PrintPos(tjs, muTj.Pts[muTj.EndPt[0]])<<"-"<<PrintPos(tjs, muTj.Pts[muTj.EndPt[1]]);
      if(muTj.PDGCode != 13) continue;
      // Found a muon, now look for delta rays
      for(unsigned short jtj = 0; jtj < tjs.allTraj.size(); ++jtj) {
        Trajectory& drTj = tjs.allTraj[jtj];
        if(drTj.AlgMod[kKilled]) continue;
        if(drTj.CTP != inCTP) continue;
        if(drTj.PDGCode == 13) continue;
        // already tagged
        if(drTj.PDGCode == 11) continue;
        // MCSMom cut
        if(drTj.MCSMom < minMom) continue;
        if(drTj.MCSMom > maxMom) continue;
        // some rough cuts to require that the delta ray is within the
        // ends of the muon
        if(muTj.StepDir > 0) {
          if(drTj.Pts[drTj.EndPt[0]].Pos[0] < muTj.Pts[muTj.EndPt[0]].Pos[0]) continue;
          if(drTj.Pts[drTj.EndPt[1]].Pos[0] > muTj.Pts[muTj.EndPt[1]].Pos[0]) continue;
        } else {
          if(drTj.Pts[drTj.EndPt[0]].Pos[0] > muTj.Pts[muTj.EndPt[0]].Pos[0]) continue;
          if(drTj.Pts[drTj.EndPt[1]].Pos[0] < muTj.Pts[muTj.EndPt[1]].Pos[0]) continue;
        }
        unsigned short muPt0, muPt1;
        float sep0 = sepCut;
        // check both ends of the prospective delta ray
        TrajPointTrajDOCA(tjs, drTj.Pts[drTj.EndPt[0]], muTj, muPt0, sep0);
        if(sep0 == sepCut) continue;
        if(prt) mf::LogVerbatim("TC")<<"  ID "<<drTj.ID<<" "<<PrintPos(tjs, drTj.Pts[drTj.EndPt[0]])<<" muPt0 "<<muPt0<<" sep0 "<<sep0;
        // stay away from the ends
        if(muPt0 < muTj.EndPt[0] + 5) continue;
        if(muPt0 > muTj.EndPt[1] - 5) continue;
        float sep1 = sepCut;
        TrajPointTrajDOCA(tjs, drTj.Pts[drTj.EndPt[1]], muTj, muPt1, sep1);
        if(prt) mf::LogVerbatim("TC")<<"      "<<PrintPos(tjs, drTj.Pts[drTj.EndPt[1]])<<" muPt1 "<<muPt1<<" sep1 "<<sep1;
        if(sep1 == sepCut) continue;
        // stay away from the ends
        if(muPt1 < muTj.EndPt[0] + 5) continue;
        if(muPt1 > muTj.EndPt[1] - 5) continue;
        if(prt) mf::LogVerbatim("TC")<<" delta ray "<<drTj.ID<<" near "<<PrintPos(tjs, muTj.Pts[muPt0]);
        drTj.ParentTrajID = muTj.ID;
        drTj.PDGCode = 11;
        // check for a vertex with another tj and if one is found, kill it
        for(unsigned short end = 0; end < 2; ++end) if(drTj.VtxID[end] > 0) MakeVertexObsolete(tjs, drTj.VtxID[end]);
      } // jtj
    } // itj
    
  } // TagDeltaRays
  
  /////////////////////////////////////////
  void TagMuonDirections(TjStuff& tjs, const short& minDeltaRayLength, short debugWorkID)
  {
    // Determine muon directions delta-ray proximity to muon trajectories
    
    if(minDeltaRayLength < 0) return;
    
    unsigned short minLen = minDeltaRayLength;
    
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& muTj = tjs.allTraj[itj];
      if(muTj.AlgMod[kKilled]) continue;
      bool prt = (debugWorkID < 0 && muTj.WorkID == debugWorkID);
      if(prt) {
        mf::LogVerbatim("TC")<<"TagMuonDirection: Muon "<<muTj.CTP<<" "<<PrintPos(tjs, muTj.Pts[muTj.EndPt[0]])<<"-"<<PrintPos(tjs, muTj.Pts[muTj.EndPt[1]]);
      }
      if(muTj.PDGCode != 13) continue;
      // look for delta ray trajectories and count the number of times that
      // one end is closer than the other to the muon
      unsigned short n0 = 0;
      unsigned short n1 = 0;
      for(unsigned short jtj = 0; jtj < tjs.allTraj.size(); ++jtj) {
        Trajectory& drTj = tjs.allTraj[jtj];
        if(drTj.AlgMod[kKilled]) continue;
        if(drTj.PDGCode != 11) continue;
        if(drTj.ParentTrajID != muTj.ID) continue;
        // ignore short delta rays
        if(drTj.Pts.size() < minLen) continue;
        float sep0 = 100;
        unsigned short muPt0;
        TrajPointTrajDOCA(tjs, drTj.Pts[drTj.EndPt[0]], muTj, muPt0, sep0);
        if(muPt0 > muTj.EndPt[1]) continue;
        float sep1 = 100;
        unsigned short muPt1;
        TrajPointTrajDOCA(tjs, drTj.Pts[drTj.EndPt[1]], muTj, muPt1, sep1);
        if(prt) mf::LogVerbatim("TC")<<" drTj.ID "<<drTj.ID<<" sep 0 "<<sep0<<" sep1 "<<sep1;
        if(muPt1 > muTj.EndPt[1]) continue;
        if(sep0 < sep1) { ++n0; } else { ++n1; }
      } // unsigned short jtj
      // Can't tell the direction using this method, so leave the current assignment unchanged
      if(prt) mf::LogVerbatim("TC")<<" n0 "<<n0<<" n1 "<<n1;
      if(n0 == n1) continue;
      if(n0 > n1) {
        // Delta-rays are closer to the beginning (0) end than the end (1) end
        muTj.TjDir = 1;
      } else {
        muTj.TjDir = -1;
      }
      if(muTj.StepDir < 0) muTj.TjDir = -muTj.TjDir;
    } // itj
  } // TagMuonDirections

  ////////////////////////////////////////////////
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag)
  {
    // Construct clusters of trajectories (cots) which will become shower PFParticles
    
    // fShowerTag[] parameters
    //  0 max Tj MCSMom for a shower tag
    //  1 max separation
    //  2 max angle
    //  3 rms width factor WF (width = WF * rms)
    //  4 Min shower 1/2 width (WSE units)
    //  5 Min total number of Tps
    //  6 Min number of Tjs
    //  7 print in plane
    
    if(fShowerTag[0] < 0) return;
    
    CTP_t printCTP = UINT_MAX;
    if(fShowerTag[7] >= 0) {
      geo::PlaneID planeID = DecodeCTP(inCTP);
      printCTP = EncodeCTP(planeID.Cryostat, planeID.TPC, std::nearbyint(fShowerTag[7]));
    }
    
    unsigned short minLenMCSMomCut = 10;

    short maxMCSMom = fShowerTag[0];
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
      // ignore nearby muons
      if(tj1.PDGCode == 13) continue;
      // Cut on length and MCSMom
      if(tj1.Pts.size() >= minLenMCSMomCut && tj1.MCSMom > maxMCSMom) continue;
      // check for a high-quality vertex at either end
      if(TjHasNiceVtx(tjs, tj1)) continue;
      // Check for proximity to Tjs in existing cots
      for(unsigned short ish = 0; ish < tjs.cots.size(); ++ish) {
        for(auto& tjID : tjs.cots[ish].TjIDs) {
          Trajectory& stj = tjs.allTraj[tjID - 1];
          float minSep = fShowerTag[1];
          unsigned short ipt1, ipt2;
          TrajTrajDOCA(tjs, tj1, stj, ipt1, ipt2, minSep, true);
          if(minSep == fShowerTag[1]) continue;
          tjs.cots[ish].TjIDs.push_back(tj1.ID);
          tj1.AlgMod[kInShower] = true;
          break;
        } // tjID
        if(tj1.AlgMod[kInShower]) break;
      } // ish
      if(tj1.AlgMod[kInShower]) continue;
      // Look for nearby tjs and possibly start a new shower
      for(unsigned short it2 = it1 + 1; it2 < tjs.allTraj.size(); ++it2) {
        Trajectory& tj2 = tjs.allTraj[it2];
        if(tj2.CTP != inCTP) continue;
        if(tj2.AlgMod[kKilled]) continue;
        // already tagged
        if(tj2.AlgMod[kInShower]) continue;
        if(tj2.PDGCode == 13) continue;
        if(tj2.StepDir != tj1.StepDir) {
          std::cout<<"FindShowers: Tj StepDirs are different\n";
          continue;
        }
        // Cut on MCSMom
        if(tj2.Pts.size() >= minLenMCSMomCut && tj2.MCSMom > maxMCSMom) continue;
        // check for a high-quality vertex at either end
        if(TjHasNiceVtx(tjs, tj2)) continue;
        float minSep = fShowerTag[1];
        unsigned short ipt1, ipt2;
        TrajTrajDOCA(tjs, tj1, tj2, ipt1, ipt2, minSep, true);
        if(minSep == fShowerTag[1]) continue;
        // Start a new shower. This is the index to a "shower" Tj that will have 3 Tps, the first for the
        // start position, the second for the charge center and the third for the end position. The number
        // of points may be increased later in DefineShowerTj
        ShowerStruct ss;
        // Create the shower trajectory
        Trajectory stj;
        stj.CTP = tj1.CTP;
        stj.ID = tjs.allTraj.size() + 1;
        stj.StepDir = tj1.StepDir;
        // Declare that it is a shower Tj
        stj.AlgMod[kShowerTj] = true;
        // define the points
        stj.Pts.resize(3);
        for(auto& stp : stj.Pts) stp.CTP = stj.CTP;
        stj.EndPt[0] = 0;
        stj.EndPt[1] = 2;
        tjs.allTraj.push_back(stj);
        ss.ShowerTjID = stj.ID;
        // now specify the IDs of the Tjs which are in the shower Tj
        ss.TjIDs.resize(2);
        ss.TjIDs[0] = tj1.ID;
        ss.TjIDs[1] = tj2.ID;
        // put it in TJ stuff
        tjs.cots.push_back(ss);
        tj1.AlgMod[kInShower] = true;
        tj2.AlgMod[kInShower] = true;
        break;
      } // it2 (tj2)
    } // it1 (tj1)
    
    if(tjs.cots.empty()) return;
    
    // Define the charge, angle, etc of the shower Tj. DefineShowerTj does not re-assign hits to the ShowerTj.
    DefineShowerTj(tjs, inCTP, USHRT_MAX, fShowerTag, printCTP);
//    MergeShowers(tjs, inCTP, fShowerTag, printCTP);
    FindShowerParent(tjs, inCTP, USHRT_MAX, fShowerTag, printCTP);
    CollectHits(tjs, inCTP, USHRT_MAX, printCTP);
    
    // drop those that don't meet the requirements
    for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
      if(tjs.cots[ic].TjIDs.empty()) continue;
      // enough Tjs?
      bool killit = (tjs.cots[ic].TjIDs.size() < fShowerTag[6]);
      unsigned short nTjWithVtx = 0;
      if(!killit) {
        // count the number of Tj points
        unsigned short nTjPts = 0;
        for(auto& tjID : tjs.cots[ic].TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          nTjPts += NumPtsWithCharge(tjs, tj, false);
          if(tj.VtxID[0] > 0 || tj.VtxID[1] > 0) ++nTjWithVtx;
        }  // tjID
        if(nTjPts < fShowerTag[5]) killit = true;
      } // !killit
      if(killit) {
        tjs.cots[ic].TjIDs.clear();
        unsigned short itj = tjs.cots[ic].ShowerTjID - 1;
        MakeTrajectoryObsolete(tjs, itj);
      }
      if(!killit && nTjWithVtx > 0) {
        // kill vertices in showers
        for(auto& tjID : tjs.cots[ic].TjIDs) {
          Trajectory& tj = tjs.allTraj[tjID - 1];
          for(unsigned short end = 0; end < 2; ++end) {
            if(tj.VtxID[end] > 0) MakeVertexObsolete(tjs, tj.VtxID[end]);
          } // end
        } // tjID
      } // !killit
    } // ic

    if(fShowerTag[7] > 0) {
      for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
        if(tjs.cots[ic].TjIDs.empty()) continue;
        unsigned short itj = tjs.cots[ic].ShowerTjID - 1;
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.CTP != printCTP) continue;
        ShowerStruct& ss = tjs.cots[ic];
        PrintTrajectory("FS", tjs, tj, USHRT_MAX);
        mf::LogVerbatim myprt("TC");
        myprt<<" Envelope";
        for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
        myprt<<"\n Tjs";
        for(auto& tjID : ss.TjIDs) {
          myprt<<" "<<tjID;
        } // tjID
      } // ic
    }

  } // FindShowers
  
  ////////////////////////////////////////////////
  void DefineShowerTj(TjStuff& tjs, const CTP_t& inCTP, const unsigned short& showerIndex, const std::vector<float>& fShowerTag, const CTP_t& printCTP)
  {
    // Defines the shower Tj variables for the showerIndex shower or all showers if
    // showerIndex is USHRT_MAX. 
    
    // fShowerTag[] parameters
    //  0 max Tj MCSMom for a shower tag
    //  1 max separation
    //  2 max angle
    //  3 rms width factor WF (width = WF * rms)
    //  4 Min shower 1/2 width (WSE units)
    //  5 Min total number of Tps
    //  6 Min number of Tjs
    //  7 print in plane
    
    if(tjs.cots.empty()) return;
    
    unsigned short first = showerIndex;
    unsigned short last = showerIndex + 1;
    if(showerIndex == USHRT_MAX) {
      first = 0;
      last = tjs.cots.size() - 1;
    }
    if(last > tjs.cots.size() - 1) return;
    
    for(unsigned short ish = first; ish < last; ++ish) {
      ShowerStruct& ss = tjs.cots[ish];
      // Ensure that it is valid
      if(ss.TjIDs.empty()) continue;
      // Ensure that this is the correct CTP
      if(tjs.allTraj[ss.TjIDs[0] - 1].CTP != inCTP) continue;
      // Define a reference to the shower trajectory
      if(ss.ShowerTjID > tjs.allTraj.size()) {
        std::cout<<"DefineShowerTj: ss.ShowerTjID not defined "<<ss.ShowerTjID<<" ish "<<ish<<"\n";
        exit(1);
      }
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      if(stj.Pts.size() != 3) {
        std::cout<<"DefineShowerTj: ss.ShowerTjID "<<ss.ShowerTjID<<" size incorrect "<<stj.Pts.size()<<"\n";
        exit(1);
      }
      // Find the trajectory end points with the largest separation and put these positions
      // into the Shower Tj end points 0 and 2. 
      TrajPoint& sTp0 = stj.Pts[0];
      TrajPoint& sTp1 = stj.Pts[1];
      TrajPoint& sTp2 = stj.Pts[2];
      // Start with the end points of the first trajectory
      unsigned short itj = ss.TjIDs[0] - 1;
      unsigned short endPt = tjs.allTraj[itj].EndPt[0];
      sTp0.Pos = tjs.allTraj[itj].Pts[endPt].Pos;
      endPt = tjs.allTraj[itj].EndPt[1];
      sTp2.Pos = tjs.allTraj[itj].Pts[endPt].Pos;
      float maxSep = PosSep2(sTp0.Pos, sTp2.Pos);
      for(unsigned short it1 = 0; it1 < ss.TjIDs.size() - 1; ++it1) {
        unsigned short itj1 = ss.TjIDs[it1] - 1;
        Trajectory& tj1 = tjs.allTraj[itj1];
        for(unsigned short end1 = 0; end1 < 2; ++end1) {
          unsigned short endPt1 = tj1.EndPt[end1];
          TrajPoint& tp1 = tj1.Pts[endPt1];
          for(unsigned short it2 = it1 + 1; it2 < ss.TjIDs.size(); ++it2) {
            Trajectory& tj2 = tjs.allTraj[ss.TjIDs[it2] - 1];
            for(unsigned short end2 = 0; end2 < 2; ++end2) {
              unsigned short endPt2 = tj2.EndPt[end2];
              TrajPoint& tp2 = tj2.Pts[endPt2];
              float sep = PosSep2(tp1.Pos, tp2.Pos);
              if(sep > maxSep) {
                maxSep = sep;
                sTp0.Pos = tp1.Pos;
                sTp2.Pos = tp2.Pos;
                // define the Tp charge even though there are no hits at these points
                sTp0.Chg = tp1.Chg;
                sTp2.Chg = tp2.Chg;
              }
            } // end2
          } // it2
        } // end1
      } // it1
      // Put sTp0 and sTp2 into a standard order
      if(sTp0.Pos[0] > sTp2.Pos[0]) std::swap(sTp0, sTp2);
      // Find the charge center
      float chgSum = 0;
      float chgPos0 = 0;
      float chgPos1 = 0;
      for(unsigned short it1 = 0; it1 < ss.TjIDs.size(); ++it1) {
        Trajectory& tj = tjs.allTraj[ss.TjIDs[it1] - 1];
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          if(tj.Pts[ipt].Chg == 0) continue;
          chgPos0 += tj.Pts[ipt].Chg * tj.Pts[ipt].Pos[0];
          chgPos1 += tj.Pts[ipt].Chg * tj.Pts[ipt].Pos[1];
          chgSum += tj.Pts[ipt].Chg;
        } // ipt
      } // it1
      TrajPoint tmp;
      MakeBareTrajPoint(tjs, sTp0, sTp2, tmp);
      sTp1.Ang = tmp.Ang;
      sTp1.AngErr = 0.2;
      sTp1.Dir = tmp.Dir;
      sTp1.Pos[0] = chgPos0 / chgSum;
      sTp1.Pos[1] = chgPos1 / chgSum;
      sTp1.Chg = chgSum;
      // Calculate the rms width of the hits wrt to the shower axis. Try to determine the
      // shower direction by finding the rms width for trajectory Tps that are closest to the start
      // and end points of the shower Tj
      sTp0.Delta = 0;
      float cnt0 = 0;
      sTp1.Delta = 0;
      sTp2.Delta = 0;
      float cnt2 = 0;
      for(unsigned short it1 = 0; it1 < ss.TjIDs.size(); ++it1) {
        Trajectory& tj = tjs.allTraj[ss.TjIDs[it1] - 1];
        for(unsigned short ipt = tj.EndPt[0]; ipt <= tj.EndPt[1]; ++ipt) {
          if(tj.Pts[ipt].Chg == 0) continue;
          float delta = PointTrajDOCA2(tjs, tj.Pts[ipt].Pos[0], tj.Pts[ipt].Pos[1], sTp1);
          sTp1.Delta += delta;
          // find rms at each end
          if(PosSep2(tj.Pts[ipt].Pos, sTp0.Pos) < PosSep2(tj.Pts[ipt].Pos, sTp2.Pos)) {
            // Tj point is closer to shower Tj point 0 
            sTp0.Delta += delta;
            ++cnt0;
          } else {
            // Tj point is closer to shower Tj point 2
            sTp2.Delta += delta;
            ++cnt2;
          }
        } // ipt
      } // it1
      // The rms width will be used to define the envelope
      if(cnt0 > 2) sTp0.Delta = fShowerTag[3] * sqrt(sTp0.Delta / cnt0);
      sTp1.Delta = fShowerTag[3] * sqrt(sTp1.Delta / (float)(cnt0 + cnt2));
      if(cnt2 > 2) sTp2.Delta = fShowerTag[3] * sqrt(sTp2.Delta / cnt2);
      // set the shower direction equal to the step direction
      stj.TjDir = 1;
      // Reverse it if the end0 width is wider than the end1 width
      if(sTp0.Delta > sTp2.Delta) {
        // Don't use ReverseTraj here or it will screw up the references to sTp0 and sTp2.
        // Simply interchange the points
        std::swap(sTp0, sTp2);
        stj.TjDir = -1;
      }
      // construct the Envelope polygon. Start with a rectangle using the fixed 1/2 width fcl input
      // expanded by the rms width at each end to create a polygon. The polygon is constructed along
      // the Pos[0] direction and then rotated into the ShowerTj direction. Use sTp1 as the origin.
      ss.Envelope.resize(4);
      // First vertex
      ss.Envelope[0][0] = sTp0.Pos[0] - sTp1.Pos[0];
      ss.Envelope[0][1] = fShowerTag[4] + 0.5 * sTp0.Delta;
      // second vertex
      ss.Envelope[1][0] = sTp2.Pos[0] - sTp1.Pos[0];
      ss.Envelope[1][1] = fShowerTag[4] + 0.5 * sTp2.Delta;
      // third and fourth are reflections of the first and second
      ss.Envelope[2][0] =  ss.Envelope[1][0];
      ss.Envelope[2][1] = -ss.Envelope[1][1];
      ss.Envelope[3][0] =  ss.Envelope[0][0];
      ss.Envelope[3][1] = -ss.Envelope[0][1];
      double rotAngle = -sTp1.Ang;
      double cs = cos(-rotAngle);
      double sn = sin(-rotAngle);
      for(auto& vtx : ss.Envelope) {
        // Rotate along the stj shower axis
        float pos0 = cs * vtx[0] - sn * vtx[1];
        float pos1 = sn * vtx[0] + cs * vtx[1];
        vtx[0] = pos0;
        vtx[1] = pos1;
        // translate
        vtx[0] += sTp1.Pos[0];
        vtx[1] += sTp1.Pos[1];
      } // vtx
      
      // Transfer the angle info to Tp 0 and 2
      sTp1.HitPos = sTp1.Pos;
      for(unsigned short ipt = 0; ipt < 3; ++ipt) {
        if(ipt == 1) continue;
        stj.Pts[ipt].Ang = sTp1.Ang;
        stj.Pts[ipt].AngErr = sTp1.AngErr;
        stj.Pts[ipt].Dir = sTp1.Dir;
        // This is only necessary so the debug printing makes sense
        stj.Pts[ipt].HitPos = stj.Pts[ipt].Pos;
      }
      if(printCTP == stj.CTP) {
        // This is done purely for printing
        stj.MCSMom = ss.TjIDs.size();
        PrintTrajectory("DSTj", tjs, stj, USHRT_MAX);
        mf::LogVerbatim myprt("TC");
        myprt<<" Envelope";
        for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
      } // print
    } // ish
  } // DefineShowerTj

  ////////////////////////////////////////////////
  void MergeShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag, const CTP_t& printCTP)
  {
    
    // fShowerTag[] parameters
    //  0 max Tj MCSMom for a shower tag
    //  1 max separation
    //  2 max angle
    //  3 rms width factor WF (width = WF * rms)
    //  4 Min shower 1/2 width (WSE units)
    //  5 Min total number of Tps
    //  6 Min number of Tjs
    //  7 print in plane
    
    if(printCTP == inCTP) PrintHeader("MS");
    
    for(unsigned short ish1 = 0; ish1 < tjs.cots.size() - 1; ++ish1) {
      // shower already merged?
      if(tjs.cots[ish1].TjIDs.empty()) continue;
      // Ensure that this is the correct CTP
      if(tjs.allTraj[tjs.cots[ish1].TjIDs[0] - 1].CTP != inCTP) continue;
      unsigned short itj1 = tjs.cots[ish1].ShowerTjID - 1;
//      if(printCTP == inCTP) PrintTrajPoint("MS", tjs, ish1, 0, 9, chgTp1);
      for(unsigned short ish2 = ish1 + 1; ish2 < tjs.cots.size(); ++ish2) {
        ShowerStruct& ss2 = tjs.cots[ish2];
        if(ss2.TjIDs.empty()) continue;
        // Ensure that this is the correct CTP
        if(tjs.allTraj[ss2.TjIDs[0] - 1].CTP != inCTP) continue;
        unsigned short itj2 = ss2.ShowerTjID - 1;
        ShowerStruct& ss1 = tjs.cots[ish1];
        if(printCTP == inCTP) {
          mf::LogVerbatim myprt("TC");
          myprt<<"Merge check ID1 "<<tjs.allTraj[itj1].ID;
          myprt<<" Envelope";
          for(auto& vtx : ss1.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
          myprt<<" ID2 "<<tjs.allTraj[itj2].ID;
          myprt<<" Envelope";
          for(auto& vtx : ss2.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
        }
        bool mergeEm = false;
        // See if any vertices in tj2 are inside tj1
        for(auto& vx2 : ss2.Envelope) {
          if(PointInsideEnvelope(vx2, ss1.Envelope)) {
            mergeEm = true;
            break;
          }
        } // vx2
        if(printCTP == inCTP) mf::LogVerbatim("TC")<<" mergeEM? "<<mergeEm;
        // Merge ss2 into ss1?
        if(!mergeEm) continue;
        ss1.TjIDs.insert(ss1.TjIDs.end(), ss2.TjIDs.begin(), ss2.TjIDs.end());
        ss2.TjIDs.clear();
        // kill the shower Tj in ss2
        unsigned short itj = ss2.ShowerTjID - 1;
        tjs.allTraj[itj].AlgMod[kKilled] = true;
        // update ss1
        DefineShowerTj(tjs, inCTP, ish1, fShowerTag, printCTP);
      } // ish2
    } // ish1
  } // MergeShowers

  ////////////////////////////////////////////////
  void CollectHits(TjStuff& tjs, const CTP_t& inCTP, const unsigned short& showerIndex, const CTP_t& printCTP)
  {
    // Collect hits in the vicinity of the shower
    
    unsigned short first = showerIndex;
    unsigned short last = showerIndex + 1;
    if(showerIndex == USHRT_MAX) {
      first = 0;
      last = tjs.cots.size() - 1;
    }
    if(last > tjs.cots.size() - 1) return;
    
    geo::PlaneID planeID = DecodeCTP(inCTP);
    unsigned short ipl = planeID.Plane;
    
    for(unsigned short ish = first; ish < last; ++ish) {
      ShowerStruct& ss = tjs.cots[ish];
      // Ensure that this is the correct CTP
      if(tjs.allTraj[ss.TjIDs[0] - 1].CTP != inCTP) continue;
      // Ensure that it is valid
      if(ss.TjIDs.empty()) continue;
      // Tp 1 of stj will get all of the shower hits
      Trajectory& stj = tjs.allTraj[ss.ShowerTjID - 1];
      // this shouldn't be necessary but do it anyway
      ReleaseHits(tjs, stj);
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
          if(TjHasNiceVtx(tjs, tj)) continue;
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

  ////////////////////////////////////////////////
  void FindShowerParent(TjStuff& tjs, const CTP_t& inCTP, const unsigned short& showerIndex, const std::vector<float>& fShowerTag, const CTP_t& printCTP)
  {
  // look for a parent trajectory, for instance long Tj with high MCSMom that enters the Tj and has
  // a small angle wrt to the shower direction
    
    // fShowerTag[] parameters
    //  0 max Tj MCSMom for a shower tag
    //  1 max separation
    //  2 max angle
    //  3 rms width factor WF (width = WF * rms)
    //  4 Min shower 1/2 width (WSE units)
    //  5 Min total number of Tps
    //  6 Min number of Tjs
    //  7 print in plane
    
    unsigned short first = showerIndex;
    unsigned short last = showerIndex + 1;
    if(showerIndex == USHRT_MAX) {
      first = 0;
      last = tjs.cots.size() - 1;
    }
    if(last > tjs.cots.size() - 1) return;
    
    unsigned short minParentLength = 4;
    
    for(unsigned short ish = first; ish < last; ++ish) {
      ShowerStruct& ss = tjs.cots[ish];
      // Ensure that this is the correct CTP
      if(tjs.allTraj[ss.TjIDs[0] - 1].CTP != inCTP) continue;
      // Ensure that it is valid
      if(ss.TjIDs.empty()) continue;
      // Reference the Tp charge center of the shower Tj
      TrajPoint& stp = tjs.allTraj[ss.ShowerTjID - 1].Pts[1];
      // Construct a Figure of Merit for finding the parent which is the DOCA * DeltaAngle / Parent Length
      float bestFOM = fShowerTag[1] * fShowerTag[2] / (float)minParentLength;
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.CTP != inCTP) continue;
        if(tj.AlgMod[kKilled]) continue;
        // it can't be a parent if it is a shower tj
        if(tj.AlgMod[kShowerTj]) continue;
        // or in a shower
        if(tj.AlgMod[kInShower]) continue;
        // or if it has low MCSMom TODO: Is this a good idea?
        if(tj.MCSMom < fShowerTag[0]) continue;
        // or if it is too short
        float npwc = NumPtsWithCharge(tjs, tj, true);
        if(npwc < minParentLength) continue;
        // make a rough DOCA and angle cut using the tj angles at both ends
        bool skipit = true;
        for(unsigned short end = 0; end < 2; ++end) {
          unsigned short endPt = tj.EndPt[end];
          if(DeltaAngle(stp.Ang, tj.Pts[endPt].Ang) < fShowerTag[2] && PointTrajDOCA(tjs, stp.Pos[0], stp.Pos[1], tj.Pts[endPt]) < 10) skipit = false;
        } // end
        if(skipit) continue;
        // find the doca between the tj and the shower charge center
        unsigned short closePt = USHRT_MAX;
        float minSep = fShowerTag[1];
        TrajPointTrajDOCA(tjs, stp, tj, closePt, minSep);
        if(closePt == USHRT_MAX) continue;
        // determine which end of the parent should be considered the start (furthest away from the shower center)
        unsigned short parentStartEnd = 0;
        if(closePt > tj.Pts.size() / 2) parentStartEnd = 1;
        // Check the angle between the start of the purported parent and the shower angle
        unsigned short startEndPt = tj.EndPt[parentStartEnd];
        float dang = DeltaAngle(tj.Pts[startEndPt].Ang, stp.Ang);
        float fom = minSep * dang / npwc;
        if(fom < bestFOM) {
          ss.ParentTjID = tj.ID;
          bestFOM = fom;
          if(inCTP == printCTP) mf::LogVerbatim("TC")<<"FSP: ish "<<ish<<" Parent tj.ID "<<tj.ID<<" minSep "<<minSep<<" dang  "<<dang<<" fom  "<<fom;
        }
      } // itj
      if(ss.ParentTjID != USHRT_MAX) {
        // We found a parent. 
/*
        // Re-find the DOCA to the shower center and split the parent so that it doesn't cross
        // the shower center Tp
        Trajectory& tj = tjs.allTraj[ss.ParentTjID - 1];
        unsigned short closePt = USHRT_MAX;
        float minSep = fShowerTag[1];
        TrajPointTrajDOCA(tjs, stp, tj, closePt, minSep);
        if(closePt > tj.EndPt[0] + 3 && closePt < tj.EndPt[1] - 3) {
          if(SplitAllTraj()) {
            
          }
        } // closePt inside the parent
*/
      } // found a parent
    } // ish

  } // FindShowerParent

 /////////////////////////////////////////
  bool MakeBareTrajPoint(TjStuff& tjs, unsigned int fromHit, unsigned int toHit, TrajPoint& tp)
  {
    CTP_t tCTP = EncodeCTP(tjs.fHits[fromHit].WireID);
    return MakeBareTrajPoint(tjs, (float)tjs.fHits[fromHit].WireID.Wire, tjs.fHits[fromHit].PeakTime,
                                  (float)tjs.fHits[toHit].WireID.Wire,   tjs.fHits[toHit].PeakTime, tCTP, tp);
    
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(TjStuff& tjs, float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp)
  {
    tp.CTP = tCTP;
    tp.Pos[0] = fromWire;
    tp.Pos[1] = tjs.UnitsPerTick * fromTick;
    tp.Dir[0] = toWire - fromWire;
    tp.Dir[1] = tjs.UnitsPerTick * (toTick - fromTick);
    float norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    if(norm == 0) return false;
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    return true;
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(TjStuff& tjs, const TrajPoint& tpIn1, const TrajPoint& tpIn2, TrajPoint& tpOut)
  {
    tpOut.CTP = tpIn1.CTP;
    tpOut.Pos = tpIn1.Pos;
    tpOut.Dir[0] = tpIn2.Pos[0] - tpIn1.Pos[0];
    tpOut.Dir[1] = tpIn2.Pos[1] - tpIn1.Pos[1];
    float norm = sqrt(tpOut.Dir[0] * tpOut.Dir[0] + tpOut.Dir[1] * tpOut.Dir[1]);
    if(norm == 0) {
      mf::LogError myprt("TC");
      myprt<<"Bad Dir in MakeBareTrajPoint ";
      myprt<<" tpIn1 Pos "<<tpIn1.Pos[0]<<" "<<tpIn1.Pos[1];
      myprt<<" tpIn2 Pos "<<tpIn2.Pos[0]<<" "<<tpIn2.Pos[1];
      return false;
    }
    tpOut.Dir[0] /= norm;
    tpOut.Dir[1] /= norm;
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
    return true;
  } // MakeBareTrajPoint
  
  ////////////////////////////////////////////////
  float TPHitsRMSTime(TjStuff& tjs, TrajPoint& tp, HitStatus_t hitRequest)
  {
    return tjs.UnitsPerTick * TPHitsRMSTick(tjs, tp, hitRequest);
  } // TPHitsRMSTime

  ////////////////////////////////////////////////
  float TPHitsRMSTick(TjStuff& tjs, TrajPoint& tp, HitStatus_t hitRequest)
  {
    // Estimate the RMS of all hits associated with a trajectory point
    // without a lot of calculation
    if(tp.Hits.empty()) return 0;
    float minVal = 9999;
    float maxVal = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tp.UseHit[ii]) useit = true;
      if(hitRequest == kUnusedHits && !tp.UseHit[ii]) useit = true;
      if(!useit) continue;
      unsigned int iht = tp.Hits[ii];
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - rms;
      if(arg < minVal) minVal = arg;
      arg = cv + rms;
      if(arg > maxVal) maxVal = arg;
    } // ii
    if(maxVal == 0) return 0;
    return (maxVal - minVal) / 2;
  } // TPHitsRMSTick
  
  ////////////////////////////////////////////////
  float HitsRMSTime(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest)
  {
    return tjs.UnitsPerTick * HitsRMSTick(tjs, hitsInMultiplet, hitRequest);
  } // HitsRMSTick

  ////////////////////////////////////////////////
  float HitsRMSTick(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest)
  {
    if(hitsInMultiplet.empty()) return 0;
    
    if(hitsInMultiplet.size() == 1) return tjs.fHits[hitsInMultiplet[0]].RMS;
 
    float minVal = 9999;
    float maxVal = 0;
    for(unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
      unsigned int iht = hitsInMultiplet[ii];
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) useit = true;
      if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) useit = true;
      if(!useit) continue;
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - rms;
      if(arg < minVal) minVal = arg;
      arg = cv + rms;
      if(arg > maxVal) maxVal = arg;
    } // ii
    if(maxVal == 0) return 0;
    return (maxVal - minVal) / 2;
  } // HitsRMSTick
  
  ////////////////////////////////////////////////
  float HitsPosTime(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, float& sum, HitStatus_t hitRequest)
  {
    return tjs.UnitsPerTick * HitsPosTick(tjs, hitsInMultiplet, sum, hitRequest);
  } // HitsPosTime
  
  ////////////////////////////////////////////////
  float HitsPosTick(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, float& sum, HitStatus_t hitRequest)
  {
    // returns the position and the charge
    float pos = 0;
    sum = 0;
    for(unsigned short ii = 0; ii < hitsInMultiplet.size(); ++ii) {
      unsigned int iht = hitsInMultiplet[ii];
      bool useit = (hitRequest == kAllHits);
      if(hitRequest == kUsedHits && tjs.fHits[iht].InTraj > 0) useit = true;
      if(hitRequest == kUnusedHits && tjs.fHits[iht].InTraj == 0) useit = true;
      if(!useit) continue;
      float chg = tjs.fHits[iht].Integral;
      pos += chg * tjs.fHits[iht].PeakTime;
      sum += chg;
    } // ii
    if(sum == 0) return 0;
    return pos / sum;
  } // HitsPosTick
  
  //////////////////////////////////////////
  unsigned short NumHitsInTP(const TrajPoint& tp, HitStatus_t hitRequest)
  {
    // Counts the number of hits of the specified type in tp
    if(tp.Hits.empty()) return 0;
    
    if(hitRequest == kAllHits) return tp.Hits.size();
    
    unsigned short nhits = 0;
    for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
      if(hitRequest == kUsedHits) {
        if(tp.UseHit[ii]) ++nhits;
      } else {
        // looking for unused hits
        if(!tp.UseHit[ii]) ++nhits;
      }
    } // ii
    return nhits;
  } // NumHitsInTP
  
  //////////////////////////////////////////
  unsigned short TPNearVertex(TjStuff& tjs, const TrajPoint& tp)
  {
    // Returns the index of a vertex if tp is nearby
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      if(tjs.vtx[ivx].NTraj == 0) continue;
      if(tjs.vtx[ivx].CTP != tp.CTP) continue;
      if(std::abs(tjs.vtx[ivx].Pos[0] - tp.Pos[0]) > 1.2) continue;
      if(std::abs(tjs.vtx[ivx].Pos[1] - tp.Pos[1]) > 1.2) continue;
      return ivx;
    } // ivx
    return USHRT_MAX;
  } // TPNearVertex
  
  //////////////////////////////////////////
  bool AttachAnyTrajToVertex(TjStuff& tjs, unsigned short ivx, const std::vector<float>& fVertex2DCuts, bool vtxPrt)
  {
    
    if(ivx > tjs.vtx.size() - 1) return false;
    if(tjs.vtx[ivx].NTraj == 0) return false;
    if(fVertex2DCuts[0] < 0) return false;
    
    VtxStore& vx = tjs.vtx[ivx];
    
    unsigned short nadd = 0;
    for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
      Trajectory& tj = tjs.allTraj[itj];
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx.CTP) continue;
      if(tj.VtxID[0] == vx.ID || tj.VtxID[1] == vx.ID) continue;
      if(AttachTrajToVertex(tjs, tj, vx, fVertex2DCuts, vtxPrt)) ++nadd;
    } // itj
    if(vtxPrt) mf::LogVerbatim("TC")<<" AttachAnyTrajToVertex: nadd "<<nadd;
    if(nadd == 0) return false;
    return true;
    
  } // AttachAnyTrajToVertex
  
  //////////////////////////////////////////
  bool AttachTrajToAnyVertex(TjStuff& tjs, unsigned short itj, const std::vector<float>& fVertex2DCuts, bool vtxPrt)
  {
    
    if(itj > tjs.allTraj.size() - 1) return false;
    if(fVertex2DCuts[0] < 0) return false;
    if(tjs.vtx.size() == 0) return false;
    
    Trajectory& tj = tjs.allTraj[itj];
    
    unsigned short nadd = 0;
    for(unsigned short ivx = 0; ivx < tjs.vtx.size(); ++ivx) {
      VtxStore& vx = tjs.vtx[ivx];
      if(vx.NTraj == 0) continue;
      if(vx.CTP != tj.CTP) continue;
      if(tj.VtxID[0] == vx.ID || tj.VtxID[1] == vx.ID) continue;
      if(AttachTrajToVertex(tjs, tj, vx, fVertex2DCuts, vtxPrt)) ++nadd;
    } // ivx
    if(nadd == 0) return false;
    return true;
    
  } // AttachAnyTrajToVertex

  //////////////////////////////////////////
  bool AttachTrajToVertex(TjStuff& tjs, Trajectory& tj, VtxStore& vx, const std::vector<float>& fVertex2DCuts, bool prt)
  {
    
    // fVertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction

    if(tj.AlgMod[kKilled]) return false;
    if(tj.CTP != vx.CTP) return false;
    // already attached?
    if(tj.VtxID[0] == vx.ID || tj.VtxID[1] == vx.ID) return false;
    
    unsigned short maxShortTjLen = fVertex2DCuts[0];
    // square the separation cut to simplify testing in the loop
    float maxSepCutShort2 = fVertex2DCuts[1] * fVertex2DCuts[1];
    float maxSepCutLong2 = fVertex2DCuts[2] * fVertex2DCuts[2];
    
    // assume that end 0 is closest to the vertex
    unsigned short end = 0;
    float vtxTjSep2 = PosSep2(vx.Pos, tj.Pts[tj.EndPt[0]].Pos);
    float sep1 = PosSep2(vx.Pos, tj.Pts[tj.EndPt[1]].Pos);
    if(sep1 < vtxTjSep2) {
      // End 1 is closer
      end = 1;
      vtxTjSep2 = sep1;
    }
    // is the trajectory short?
    bool tjShort = (tj.EndPt[1] - tj.EndPt[0] < maxShortTjLen);
    // ignore bad separation between the closest tj end and the vertex
    if(tjShort) {
      if(vtxTjSep2 > maxSepCutShort2) return false;
    } else {
      if(vtxTjSep2 > maxSepCutLong2) return false;
    }
    
    // Calculate the pull on the vertex
    TrajPoint& tp = tj.Pts[tj.EndPt[end]];
    float tpVxPull = TrajPointVertexPull(tjs, tp, vx);

    // See if the vertex position is close to an end
    unsigned short closePt;
    float closestApproach;
    TrajClosestApproach(tj, vx.Pos[0], vx.Pos[1], closePt, closestApproach);
    // count the number of points between the end of the trajectory and the vertex.
    // tj     -------------   tj ------------
    // vx  *   >> dpt = 0     vx   *  >> dpt = 2
    short dpt;
    if(end == 0) {
      dpt = closePt - tj.EndPt[end];
    } else {
      dpt = tj.EndPt[end] - closePt;
    }
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"ATTV: vx.ID "<<vx.ID;
      myprt<<" oldTJs";
      for(unsigned short itj = 0; itj < tjs.allTraj.size(); ++itj) {
        Trajectory& tj = tjs.allTraj[itj];
        if(tj.AlgMod[kKilled]) continue;
        if(tj.CTP != vx.CTP) continue;
        if(tj.VtxID[0] == vx.ID) myprt<<" "<<tj.ID<<"_0";
        if(tj.VtxID[1] == vx.ID) myprt<<" "<<tj.ID<<"_1";
      }
      myprt<<" +tjID "<<tj.ID<<"_"<<end<<" vtxTjSep "<<sqrt(vtxTjSep2)<<" tpVxPull "<<tpVxPull<<" fVertex2DCuts[3] "<<fVertex2DCuts[3];
    }
    if(tpVxPull > fVertex2DCuts[3]) return false;
    if(dpt > 2) return false;

    // Passed all the cuts. Attach it to the vertex and try a fit
    tj.VtxID[end] = vx.ID;
    if(FitVertex(tjs, vx, fVertex2DCuts, prt)) {
      if(prt) mf::LogVerbatim("TC")<<" success";
      return true;
    }
    
    // fit failed so remove the tj -> vx assignment
    tj.VtxID[end] = 0;
    // and refit
    if(prt) mf::LogVerbatim("TC")<<" failed. Re-fit w/o this tj ";
    FitVertex(tjs, vx, fVertex2DCuts, prt);
    return false;
    
  } // AttachTrajToVertex

  /////////////////////////////////////////
  float TrajPointVertexPull(TjStuff& tjs, const TrajPoint& tp, const VtxStore& vx)
  {
    // Calculates the position pull between a trajectory point and a vertex

    // impact parameter between them
    double ip = PointTrajDOCA(tjs, vx.Pos[0], vx.Pos[1], tp);
    // separation^2
    double sep2 = PosSep2(vx.Pos, tp.Pos);
    
    // Find the projection of the vertex error ellipse in a coordinate system
    // along the TP direction
    double vxErrW = vx.PosErr[0] * tp.Dir[1];
    double vxErrT = vx.PosErr[1] * tp.Dir[0];
    double vxErr2 = vxErrW * vxErrW + vxErrT * vxErrT;
    
    // close together so ignore the TP projection error and return
    // the pull using only the vertex error
    if(sep2 < 1) return (float)(ip/sqrt(vxErr2));
    
    double dang = ip / sqrt(sep2);

    // calculate the angle error.
    // Start with the vertex error^2
    double angErr = vxErr2 / sep2;
    // Add the TP angle error^2
    angErr += tp.AngErr * tp.AngErr;
    if(angErr == 0) return 999;
    angErr = sqrt(angErr);
    return (float)(dang / angErr);
    
  } // TrajPointVertexPull
  
  /////////////////////////////////////////
  float VertexVertexPull(TjStuff& tjs, const VtxStore& vx1, const VtxStore& vx2)
  {
    // Calculates the position pull between two vertices
    double dw = vx1.Pos[0] - vx2.Pos[0];
    double dt = vx1.Pos[1] - vx2.Pos[1];
    double dwErr2 = (vx1.PosErr[0] * vx1.PosErr[0] + vx2.PosErr[0] * vx2.PosErr[0]) / 2;
    double dtErr2 = (vx1.PosErr[1] * vx1.PosErr[1] + vx2.PosErr[1] * vx2.PosErr[1]) / 2;
    dw = dw * dw / dwErr2;
    dt = dt * dt / dtErr2;
    return (float)sqrt(dw + dt);
  }
  
  /////////////////////////////////////////
  bool FitVertex(TjStuff& tjs, VtxStore& vx, const std::vector<float>& fVertex2DCuts, bool prt)
  {
    // A poor-mans fitting scheme. If the fitted vertex position error is within the supplied
    // value, the position and errors are updated and we return true, otherwise the vertex is
    // left unchanged and we return false
    
    // fVertex2DCuts fcl input usage
    // 0 = maximum length of a short trajectory
    // 1 = max vertex - trajectory separation for short trajectories
    // 2 = max vertex - trajectory separation for long trajectories
    // 3 = max position pull for adding TJs to a vertex
    // 4 = max allowed vertex position error
    // 5 = min MCSMom
    // 6 = min Pts/Wire fraction
    
    // Create a vector of trajectory points that will be used to fit the vertex position
    std::vector<TrajPoint> vxTp;
    for(auto& tj : tjs.allTraj) {
      if(tj.AlgMod[kKilled]) continue;
      if(tj.CTP != vx.CTP) continue;
      if(tj.VtxID[0] == vx.ID) vxTp.push_back(tj.Pts[tj.EndPt[0]]);
      if(tj.VtxID[1] == vx.ID) vxTp.push_back(tj.Pts[tj.EndPt[1]]);
    } // tj
    
    vx.NTraj = vxTp.size();
    
    if(vxTp.size() < 2) return false;
    
    if(prt) {
      PrintHeader("FV");
      for(auto& tp : vxTp) PrintTrajPoint("FV", tjs, 0, 1, 1, tp);
    }
    if(vx.Stat[kFixed]) {
      if(prt) mf::LogVerbatim("TC")<<" vertex position fixed. No fit.";
      return true;
    }
 
    // Find trajectory intersections pair-wise tweaking the angle and position(?) within
    // +/- 1 sigma
    double sum0 = 0, sum02 = 0;
    double sum1 = 0, sum12 = 0;
    double sumw = 0;
    double wgt;
    // a temporary TP for tweaking the angle
    TrajPoint tmp;
    for(unsigned short itj = 0; itj < vxTp.size() - 1; ++itj) {
      for(unsigned short jtj = itj + 1; jtj < vxTp.size(); ++jtj) {
        float p0, p1;
        TrajIntersection(vxTp[itj], vxTp[jtj], p0, p1);
        // accumulate
        wgt = 1;
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
        // tweak the itj angle +
        tmp = vxTp[itj];
        tmp.Ang += tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(tmp, vxTp[jtj], p0, p1);
        // accumulate
        // adjust the weight for 4 points at +/1 1 sigma = 0.607 / 4
        wgt = 0.152;
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
        // tweak the itj angle -
        tmp = vxTp[itj];
        tmp.Ang -= 2 * tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(tmp, vxTp[jtj], p0, p1);
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
        // Repeat this process with jtj
        // tweak the jtj angle +
        tmp = vxTp[jtj];
        tmp.Ang += tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(vxTp[itj], tmp, p0, p1);
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
        // tweak the itj angle -
        tmp = vxTp[itj];
        tmp.Ang -= 2 * tmp.AngErr;
        tmp.Dir[0] = cos(tmp.Ang); tmp.Dir[1] = sin(tmp.Ang);
        TrajIntersection(vxTp[itj], tmp, p0, p1);
        // accumulate
        sum0 += wgt * p0; sum02 += wgt * p0 * p0; sum1 += wgt * p1; sum12 += wgt * p1 * p1; sumw += wgt;
      } // jtj
    } // itj
    
    double vxP0 = sum0 / sumw;
    double vxP1 = sum1 / sumw;
    double vxP0rms = sqrt((sum02 - sumw * vxP0 * vxP0) / sumw);
    double vxP1rms = sqrt((sum12 - sumw * vxP1 * vxP1) / sumw);
    // don't let the errors get too small
    if(vxP0rms < 0.5) vxP0rms = 0.5;
    if(vxP1rms < 0.5) vxP1rms = 0.5;
    
    if(prt) mf::LogVerbatim("TC")<<"FitVertex "<<vx.ID<<" CTP "<<vx.CTP<<" NTraj "<<vx.NTraj<<" in "<<std::fixed<<std::setprecision(1)<<vx.Pos[0]<<":"<<vx.Pos[1]/tjs.UnitsPerTick<<" out "<<vxP0<<"+/-"<<vxP0rms<<":"<<vxP1/tjs.UnitsPerTick<<"+/-"<<vxP1rms/tjs.UnitsPerTick;
    
    if(vxP0rms > fVertex2DCuts[4] || vxP1rms > fVertex2DCuts[4]) {
      if(prt) mf::LogVerbatim("TC")<<" fit failed. fVertex2DCuts[4] "<<fVertex2DCuts[4];
      return false;
    }
    
    vx.Pos[0] = vxP0;
    vx.PosErr[0] = vxP0rms;
    vx.Pos[1] = vxP1;
    vx.PosErr[1] = vxP1rms;
    
    // Calculate chisq
    vx.ChiDOF = 0;
    for(unsigned short itj = 0; itj < vxTp.size(); ++itj) {
      vx.ChiDOF += TrajPointVertexPull(tjs, vxTp[itj], vx);
    } // itj
    vx.ChiDOF /= (float)vxTp.size();
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"Pull";
      for(unsigned short itj = 0; itj < vxTp.size(); ++itj) {
        float pull = TrajPointVertexPull(tjs, vxTp[itj], vx);
        myprt<<" "<<PrintPos(tjs, vxTp[itj])<<"-"<<std::fixed<<std::setprecision(2)<<pull;
      } // itj
      myprt<<" ChiDOF "<<vx.ChiDOF;
    }
    return true;
    
  } // FitVertex
  
  // ****************************** Printing  ******************************
  
  void PrintAllTraj(std::string someText, TjStuff& tjs, DebugStuff& debug, unsigned short itj, unsigned short ipt)
  {
    
    mf::LogVerbatim myprt("TC");
    
    if(!tjs.vtx3.empty()) {
      // print out 3D vertices
      myprt<<"****** 3D vertices ******************************************__2DVtx_ID__*******\n";
      myprt<<"Vtx  Cstat  TPC     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for(unsigned short iv = 0; iv < tjs.vtx3.size(); ++iv) {
        if(tjs.vtx3[iv].Wire == SHRT_MAX) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(7)<<tjs.vtx3[iv].CStat;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].TPC;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].X;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Y;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Z;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].XErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].YErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].ZErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[0]+1;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[1]+1;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[2]+1;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Wire;
        if(tjs.vtx3[iv].Wire < 0) {
          myprt<<"    Matched in all planes";
        } else {
          myprt<<"    Incomplete";
        }
        myprt<<"\n";
      }
    } // tjs.vtx3.size
    
    if(!tjs.vtx.empty()) {
      // print out 2D vertices
      myprt<<"************ 2D vertices ************\n";
      myprt<<"Vtx   CTP    wire     error   tick     error  ChiDOF  NTj Pass  Topo Nice? traj IDs\n";
      for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
        auto const& aVtx = tjs.vtx[iv];
        if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aVtx.CTP).Plane) continue;
        if(aVtx.NTraj == 0) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<aVtx.ID<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<aVtx.CTP;
        myprt<<std::right<<std::setw(8)<<aVtx.Pos[0]<<" +/- ";
        myprt<<std::right<<std::setw(4)<<aVtx.PosErr[0];
        myprt<<std::right<<std::setw(8)<<aVtx.Pos[1]/tjs.UnitsPerTick<<" +/- ";
        myprt<<std::right<<std::setw(4)<<aVtx.PosErr[1]/tjs.UnitsPerTick;
        myprt<<std::right<<std::setw(8)<<aVtx.ChiDOF;
        myprt<<std::right<<std::setw(5)<<aVtx.NTraj;
        myprt<<std::right<<std::setw(5)<<aVtx.Pass;
        myprt<<std::right<<std::setw(6)<<aVtx.Topo;
        myprt<<std::right<<std::setw(6)<<aVtx.Stat[kNiceVtx];
        myprt<<"    ";
        // display the traj indices
        for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
          auto const& aTj = tjs.allTraj[ii];
          if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
          if(aTj.AlgMod[kKilled]) continue;
          for(unsigned short end = 0; end < 2; ++end)
            if(aTj.VtxID[end] == (short)aVtx.ID) myprt<<std::right<<std::setw(4)<<aTj.ID<<"_"<<end;
        }
        myprt<<"\n";
      } // iv
    } // tjs.vtx.size
    
    if(tjs.allTraj.empty()) {
      mf::LogVerbatim("TC")<<someText<<" No allTraj trajectories to print";
      return;
    }
    
    // Print all trajectories in tjs.allTraj if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      std::vector<unsigned int> tmp;
      myprt<<someText<<" TRJ  ID CTP Pass Pts frm   to     W:Tick   Ang C AveQ     W:T      Ang C AveQ ChgRMS  Mom sDr tDr __Vtx__ PDG  Par TRuPDG  E*P TruKE  WorkID \n";
      for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
        auto const& aTj = tjs.allTraj[ii];
        if(debug.Plane >=0 && debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
        myprt<<someText<<" ";
        if(aTj.AlgMod[kKilled]) { myprt<<"xxx"; } else { myprt<<"TRJ"; }
        myprt<<std::fixed<<std::setw(4)<<aTj.ID;
        myprt<<std::setw(3)<<aTj.CTP;
        myprt<<std::setw(5)<<aTj.Pass;
        myprt<<std::setw(5)<<aTj.Pts.size();
        myprt<<std::setw(4)<<aTj.EndPt[0];
        myprt<<std::setw(5)<<aTj.EndPt[1];
        int endPt = aTj.EndPt[0];
        TrajPoint tp = aTj.Pts[endPt];
        int itick = tp.Pos[1]/tjs.UnitsPerTick;
        if(itick < 0) itick = 0;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(2)<<tp.AngleCode;
        myprt<<std::setw(5)<<(int)tp.AveChg;
        endPt = aTj.EndPt[1];
        tp = aTj.Pts[endPt];
        itick = tp.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(2)<<tp.AngleCode;
        myprt<<std::setw(5)<<(int)tp.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<aTj.ChgRMS;
        myprt<<std::setw(5)<<aTj.MCSMom;
        myprt<<std::setw(4)<<aTj.StepDir;
        myprt<<std::setw(3)<<aTj.TjDir;
        myprt<<std::setw(4)<<aTj.VtxID[0];
        myprt<<std::setw(4)<<aTj.VtxID[1];
        myprt<<std::setw(5)<<aTj.PDGCode;
        myprt<<std::setw(5)<<aTj.ParentTrajID;
        myprt<<std::setw(6)<<aTj.TruPDG;
        myprt<<std::setw(6)<<std::setprecision(2)<<aTj.EffPur;
        myprt<<std::setw(5)<<(int)aTj.TruKE;
        myprt<<std::setw(7)<<aTj.WorkID;
//        myprt<<" "<<PrintStopFlag(tjs, aTj, 0)<<" "<<PrintStopFlag(tjs, aTj, 1)<<" ";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(aTj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
      } // ii
      return;
    } // itj > tjs.allTraj.size()-1
    
    if(itj > tjs.allTraj.size()-1) return;
    
    auto const& aTj = tjs.allTraj[itj];
    
    mf::LogVerbatim("TC")<<"Print tjs.allTraj["<<itj<<"]: ClusterIndex "<<aTj.ClusterIndex<<" Vtx[0] "<<aTj.VtxID[0]<<" Vtx[1] "<<aTj.VtxID[1];
    myprt<<"AlgBits";
    for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(aTj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
    myprt<<"\n";
    
    PrintHeader(someText);
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < aTj.Pts.size(); ++ii) PrintTrajPoint(someText, tjs, ii, aTj.StepDir, aTj.Pass, aTj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(someText, tjs, ipt, aTj.StepDir, aTj.Pass, aTj.Pts[ipt]);
    }
  } // PrintAllTraj
  
  
  //////////////////////////////////////////
  void PrintTrajectory(std::string someText, TjStuff& tjs, Trajectory const& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    if(tPoint == USHRT_MAX) {
      if(tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ";
        myprt<<"Work:    ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.VtxID[0]<<" "<<tj.VtxID[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        myprt<<" MCSMom "<<tj.MCSMom;
        myprt<<" StopFlags "<<PrintStopFlag(tjs, tj, 0)<<" "<<PrintStopFlag(tjs, tj, 1);
        myprt<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      } else {
        mf::LogVerbatim myprt("TC");
        myprt<<"tjs.allTraj: ID "<<tj.ID<<" WorkID "<<tj.WorkID<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.VtxID[0]<<" "<<tj.VtxID[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        myprt<<" MCSMom "<<tj.MCSMom;
        myprt<<" StopFlags "<<PrintStopFlag(tjs, tj, 0)<<" "<<PrintStopFlag(tjs, tj, 1);
        myprt<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      }
      PrintHeader(someText);
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) PrintTrajPoint(someText, tjs, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
    } else {
      // just print one traj point
      if(tPoint > tj.Pts.size() -1) {
        mf::LogVerbatim("TC")<<"Can't print non-existent traj point "<<tPoint;
        return;
      }
      PrintTrajPoint(someText, tjs, tPoint, tj.StepDir, tj.Pass, tj.Pts[tPoint]);
    }
  } // PrintTrajectory
  
  //////////////////////////////////////////
  void PrintHeader(std::string someText)
  {
    mf::LogVerbatim("TC")<<someText<<" TRP  CTP  Ind  Stp      W:Tick    Delta  RMS    Ang C   Err  Dir0  Dir1      Q    AveQ  Pull FitChi  NTPF  Hits ";
  } // PrintHeader
  
  ////////////////////////////////////////////////
  void PrintTrajPoint(std::string someText, TjStuff& tjs, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" TRP"<<std::fixed;
    myprt<<pass;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(3)<<tp.CTP;
    myprt<<std::setw(5)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(7)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]/tjs.UnitsPerTick; // W:T
    if(tp.Pos[1] < 10) myprt<<"  "; if(tp.Pos[1] < 100) myprt<<" "; if(tp.Pos[1] < 1000) myprt<<" ";
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Delta;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.DeltaRMS;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(2)<<tp.AngleCode;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[0];
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[1];
    myprt<<std::setw(7)<<(int)tp.Chg;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgPull;
    myprt<<std::setw(7)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NTPsFit;
    // print the hits associated with this traj point
    if(tp.Hits.size() > 16) {
      // don't print too many hits (e.g. from a shower Tj)
      myprt<<" "<<tp.Hits.size()<<" shower hits";
    } else {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) {
        unsigned int iht = tp.Hits[ii];
        myprt<<" "<<tjs.fHits[iht].WireID.Wire<<":"<<(int)tjs.fHits[iht].PeakTime;
        if(tp.UseHit[ii]) {
          // Distinguish used hits from nearby hits
          myprt<<"_";
        } else {
          myprt<<"x";
        }
        myprt<<tjs.fHits[iht].InTraj;
      } // iht
    }
  } // PrintTrajPoint
  
  /////////////////////////////////////////
  std::string PrintStopFlag(TjStuff& tjs, const Trajectory& tj, unsigned short end)
  {
    if(end > 1) return "Invalid end";
    std::string tmp;
    bool first = true;
    for(unsigned short ib = 0; ib < StopFlagNames.size(); ++ib) {
      if(tj.StopFlag[end][ib]) {
        if(first) {
          tmp = std::to_string(end) + ":" + StopFlagNames[ib];
          first = false;
        } else {
          tmp += "," + StopFlagNames[ib];
        }
      }
    } // ib
    return tmp;
  } // PrintStopFlag
  
  /////////////////////////////////////////
  std::string PrintHitShort(const TCHit& hit)
  {
    return std::to_string(hit.WireID.Plane) + ":" + std::to_string(hit.WireID.Wire) + ":" + std::to_string((int)hit.PeakTime);
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintHit(const TCHit& hit)
  {
    return std::to_string(hit.WireID.Plane) + ":" + std::to_string(hit.WireID.Wire) + ":" + std::to_string((int)hit.PeakTime) + "_" + std::to_string(hit.InTraj);
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintPos(TjStuff& tjs, const TrajPoint& tp)
  {
    return PrintPos(tjs, tp.Pos);
  } // PrintPos
  
  /////////////////////////////////////////
  std::string PrintPos(TjStuff& tjs, const std::array<float, 2>& pos)
  {
    unsigned int wire = std::nearbyint(pos[0]);
    int time = std::nearbyint(pos[1]/tjs.UnitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

  
} // namespace tca

