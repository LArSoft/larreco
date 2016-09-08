#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {
/*
  /////////////////////////////////////////
  bool SignalPresent(TjStuff& tjs, TrajPoint const& tp, float minAmp)
  {
    return SignalPresent(tjs, tp.Pos[1], tp.Pos[0], tp.Pos[1], tp.Pos[0], tp.CTP, minAmp);
  }
  
  /////////////////////////////////////////
  bool SignalPresent(TjStuff& tjs, float wire1, float time1, TrajPoint const& tp, float minAmp)
  {
    unsigned int w1 = std::nearbyint(wire1);
    unsigned int w2 = std::nearbyint(tp.Pos[0]);
    return SignalPresent(tjs, w1, time1, w2, tp.Pos[1], tp.CTP, minAmp);
  }
  
  /////////////////////////////////////////
  bool SignalPresent(TjStuff& tjs, float wire1, float time1, float wire2, float time2, CTP_t pCTP, float minAmp)
  {
    unsigned int w1 = std::nearbyint(wire1);
    unsigned int w2 = std::nearbyint(wire2);
    return SignalPresent(tjs, w1, time1, w2, time2, pCTP, minAmp);
  } // SignalPresent
  
  /////////////////////////////////////////
  bool SignalPresent(TjStuff& tjs, unsigned int wire1, float time1, unsigned int wire2, float time2, CTP_t pCTP, float minAmp)
  {
    // returns  true if there is a signal on the line between (wire1, time1) and (wire2, time2).
    
    // Gaussian amplitude in bins of size 0.15
    const float gausAmp[20] = {1, 0.99, 0.96, 0.90, 0.84, 0.75, 0.67, 0.58, 0.49, 0.40, 0.32, 0.26, 0.20, 0.15, 0.11, 0.08, 0.06, 0.04, 0.03, 0.02};
    
    // convert time to tick
    time1 /= tjs.UnitsPerTick;
    time2 /= tjs.UnitsPerTick;
    //    mf::LogVerbatim("TC")<<"SignalPresent: check "<<wire1<<":"<<(int)time1<<" to "<<wire2<<":"<<(int)time2;
    
    // get the begin and end right
    unsigned int wireb = wire1;
    float timeb = time1;
    unsigned int wiree = wire2;
    float timee = time2;
    // swap them?
    if(wiree > wireb) {
      wireb = wire2;
      timeb = time2;
      wiree = wire1;
      timee = time1;
    }
    
    geo::PlaneID planeID = DecodeCTP(pCTP);
    unsigned int ipl = planeID.Plane;
    
    if(wiree < tjs.FirstWire[ipl] || wiree > tjs.LastWire[ipl]) return false;
    if(wireb < tjs.FirstWire[ipl] || wireb > tjs.LastWire[ipl]) return false;
    
    int wire0 = wiree;
    // checking a single wire?
    float slope = 0;
    bool oneWire = false;
    float prTime, prTimeLo = 0, prTimeHi = 0;
    if(wireb == wiree) {
      oneWire = true;
      if(time1 < time2) {
        prTimeLo = time1;
        prTimeHi = time2;
      } else {
        prTimeLo = time2;
        prTimeHi = time1;
      }
    } else {
      slope = (timeb - timee) / (wireb - wiree);
    }

    
    int bin;
    for(unsigned int wire = wiree; wire < wireb + 1; ++wire) {
      if(oneWire) {
        prTime = (prTimeLo + prTimeHi) / 2;
      } else {
        prTime = timee + (wire - wire0) * slope;
      }
      if (wire >= tjs.NumWires[ipl]) continue;
      // skip dead wires
      if(tjs.WireHitRange[ipl][wire].first == -1) continue;
      // no hits on this wire
      if(tjs.WireHitRange[ipl][wire].first == -2) return false;
      unsigned int firsthit = tjs.WireHitRange[ipl][wire].first;
      unsigned int lasthit = tjs.WireHitRange[ipl][wire].second;
      //      mf::LogVerbatim("TC")<<" wire "<<wire<<" Hit range "<<firsthit<<" "<<lasthit<<" prTime "<<prTime;
      float amp = 0;
      for(unsigned int khit = firsthit; khit < lasthit; ++khit) {
        //        mf::LogVerbatim("TC")<<"  hit "<<PrintHit(tjs.fHits[khit])<<" rms "<<tjs.fHits[khit]->RMS()<<" amp "<<(int)tjs.fHits[khit]->PeakAmplitude()<<" StartTick "<<tjs.fHits[khit]->StartTick()<<" EndTick "<<tjs.fHits[khit]->EndTick();
        if(oneWire) {
          // TODO: This sometimes doesn't work with overlapping hits
          //            if(prTimeHi > tjs.fHits[khit].EndTick()) continue;
          //            if(prTimeLo < tjs.fHits[khit].StartTick()) continue;
          // A not totally satisfactory solution
          if(prTime < tjs.fHits[khit]->StartTick()) continue;
          if(prTime > tjs.fHits[khit]->EndTick()) continue;
          return true;
        } else {
          // skip checking if we are far away from prTime on the positive side
          if(tjs.fHits[khit]->PeakTime() - prTime > 500) continue;
          bin = std::abs(tjs.fHits[khit]->PeakTime() - prTime) / tjs.fHits[khit]->RMS();
          bin /= 0.15;
          if(bin > 19) continue;
          if(bin < 0) continue;
          //          mf::LogVerbatim("CC")<<"  bin "<<bin<<" add "<<tjs.fHits[khit]->PeakAmplitude() * gausAmp[bin]<<" to amp "<<amp;
          // add amplitude from all hits
          amp += tjs.fHits[khit]->PeakAmplitude() * gausAmp[bin];
        }
      } // khit
      //      mf::LogVerbatim("TC")<<"Amp "<<amp<<" fMinAmp "<<fMinAmp;
      if(amp < minAmp) return false;
    } // wire
    return true;
    
  } // SignalPresent
*/
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
        if(tjs.inTraj[iht] == tjs.allTraj[itj].ID) tjs.inTraj[iht] = 0;
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
          if(tjs.inTraj[iht] == 0) {
            tjs.inTraj[iht] = tjs.allTraj[itj].ID;
          }
        }
      } // ii
    } // tp
    tjs.allTraj[itj].AlgMod[kKilled] = false;
  } // RestoreObsoleteTrajectory

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
        if(tjs.inTraj[iht] != tj.ID) continue;
        tjs.inTraj[iht] = newTj.ID;
        tj.Pts[ipt].UseHit[ii] = false;
      } // ii
    } // ipt
    SetEndPoints(tjs, tj);
    if(ivx != USHRT_MAX) tj.Vtx[1] = ivx;
    tj.AlgMod[kSplitTraj] = true;
    if(prt) {
      mf::LogVerbatim("TC")<<"Splittjs.allTraj: itj "<<tj.ID<<" EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
//      PrintTrajectory(tjs, tjs.allTraj[itj], USHRT_MAX);
    }
    
    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if(eraseSize > newTj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"Splittjs.allTraj: Bad erase size ";
      return false;
    }
    
    // erase the TPs at the beginning of the new trajectory
    //    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + pos + 1);
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for(ipt = 0; ipt < 3; ++ipt) {
      for(ii = 0; ii < newTj.Pts[ipt].UseHit.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
    } // ipt
    SetEndPoints(tjs, newTj);
    if(ivx != USHRT_MAX) newTj.Vtx[0] = ivx;
    newTj.AlgMod[kSplitTraj] = true;
    tjs.allTraj.push_back(newTj);
    if(prt) {
      mf::LogVerbatim("TC")<<"Splittjs.allTraj: NewTj "<<newTj.ID<<" EndPts "<<newTj.EndPt[0]<<" to "<<newTj.EndPt[1];
//      PrintTrajectory(tjs, newTj, USHRT_MAX);
    }
    return true;
    
  } // SplitAllTraj
  
  //////////////////////////////////////////
  void TrajPointTrajDOCA(TjStuff& tjs, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep)
  {
    // Finds the point, ipt, on trajectory tj that is closest to trajpoint tp
    float best = minSep * minSep;
    closePt = 0;
    float dw, dt, dp2;
    unsigned short ipt;
    for(ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
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
  void TrajTrajDOCA(Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
  {
    // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
    float best = minSep * minSep;
    ipt1 = 0; ipt2 = 0;
    float dw, dt, dp2;
    unsigned short i1, i2;
    for(i1 = tj1.EndPt[0]; i1 < tj1.EndPt[1] + 1; ++i1) {
      for(i2 = tj2.EndPt[0]; i2 < tj2.EndPt[1] + 1; ++i2) {
        // TODO: What about TPs with UseHit = false?
        dw = tj1.Pts[i1].Pos[0] - tj2.Pts[i2].Pos[0];
        dt = tj1.Pts[i1].Pos[1] - tj2.Pts[i2].Pos[1];
        dp2 = dw * dw + dt * dt;
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
  float TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2)
  {
    // returns the separation^2 between the hit position of two trajectory points
    float dw = tp1.HitPos[0] - tp2.HitPos[0];
    float dt = tp1.HitPos[1] - tp2.HitPos[1];
    return dw * dw + dt * dt;
  } // TrajPointHitSep2
  
  //////////////////////////////////////////
  float HitSep2(TjStuff& tjs, unsigned int iht, unsigned int jht)
  {
    // returns the separation^2 between two hits in WSE units
    if(iht > tjs.fHits.size()-1 || jht > tjs.fHits.size()-1) return 1E6;
    float dw = (float)tjs.fHits[iht]->WireID().Wire - (float)tjs.fHits[jht]->WireID().Wire;
    float dt = (tjs.fHits[iht]->PeakTime() - tjs.fHits[jht]->PeakTime()) * tjs.UnitsPerTick;
    return dw * dw + dt * dt;
  } // TrajPointHitSep2
  
  //////////////////////////////////////////
  float PointTrajSep2(float wire, float time, TrajPoint const& tp)
  {
    float dw = wire - tp.Pos[0];
    float dt = time - tp.Pos[1];
    return dw * dw + dt * dt;
  }
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff& tjs, unsigned int iht, TrajPoint const& tp)
  {
    float wire = tjs.fHits[iht]->WireID().Wire;
    float time = tjs.fHits[iht]->PeakTime() * tjs.UnitsPerTick;
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA(TjStuff& tjs, float wire, float time, TrajPoint const& tp)
  {
    return sqrt(PointTrajDOCA2(tjs, wire, time, tp));
  } // PointTrajDOCA
  
  //////////////////////////////////////////
  float PointTrajDOCA2(TjStuff& tjs, float wire, float time, TrajPoint const& tp)
  {
    // returns the distance of closest approach squared between a (wire, time(WSE)) point
    // and a trajectory point
    
    float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
    float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
    float dt = tp.Pos[1] + t * tp.Dir[1] - time;
    return (dw * dw + dt * dt);
    
  } // PointTrajDOCA2
  
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
  void PutTrajHitsInVector(Trajectory const& tj, bool onlyUsedHits, std::vector<unsigned int>& hitVec)
  {
    // Put hits in each trajectory point into a flat vector. Only hits with UseHit if onlyUsedHits == true
    hitVec.clear();
    hitVec.reserve(tj.Pts.size());
    unsigned short ipt, ii;
    unsigned int iht;
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        iht = tj.Pts[ipt].Hits[ii];
        if(onlyUsedHits) {
          if(tj.Pts[ipt].UseHit[ii]) hitVec.push_back(iht);
        } else {
          hitVec.push_back(iht);
        }
      } // iht
    } // ipt
  } // PutTrajHitsInVector
  
  //////////////////////////////////////////
  bool HitIsInTj(Trajectory const& tj, const unsigned int& iht, short nPtsToCheck, bool prt)
  {
    // returns true if hit iht is associated with trajectory tj. Checking starts at the
    // end of tj for nPtsToCheck points
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      unsigned short ipt = tj.Pts.size() - 1 - ii;
/*
      if(prt) {
        std::cout<<"HitIsInTJ: iht "<<iht<<" ipt "<<ipt<<" size "<<tj.Pts.size()<<" jhts";
        for(auto jht : tj.Pts[ipt].Hits) std::cout<<" "<<jht;
        std::cout<<"\n";
      }
*/
      if(std::find(tj.Pts[ipt].Hits.begin(), tj.Pts[ipt].Hits.end(), iht) != tj.Pts[ipt].Hits.end()) return true;
      // only go back a few points
      if(nPtsToCheck >= 0 && ii == nPtsToCheck) return false;
      if(ipt == 0) return false;
    } // ii
    return false;

  } // HitIsInTj
  
  //////////////////////////////////////////
  void MoveTPToWire(TrajPoint& tp, float wire)
  {
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    if(tp.Dir[0] == 0) return;
    float dw = wire - tp.Pos[0];
    if(dw == 0) return;
    tp.Pos[0] = wire;
    tp.Pos[1] += dw * tp.Dir[1] / tp.Dir[0];
  } // MoveTPToWire

  //////////////////////////////////////////
  float DeltaAngle(float Ang1, float Ang2) {
    return std::abs(std::remainder(Ang1 - Ang2, M_PI));
  }
  
  ////////////////////////////////////////////////
  void SetEndPoints(TjStuff& tjs, Trajectory& tj)
  {
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
    
    tj.EndPt[0] = 0; tj.EndPt[1] = 0;
    if(tj.Pts.size() == 0) return;
    
    // check the end point pointers
    unsigned short ipt, ii;
    // make sure that the Chg is set correctly
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      tj.Pts[ipt].Chg = 0;
      for(ii = 0; ii < tj.Pts[ipt].UseHit.size(); ++ii)
        if(tj.Pts[ipt].UseHit[ii]) tj.Pts[ipt].Chg += tjs.fHits[tj.Pts[ipt].Hits[ii]]->Integral();
    } // ipt
    
    for(ipt = 0; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[0] = ipt;
        break;
      }
    }
    for(unsigned short ii = 0; ii < tj.Pts.size(); ++ii) {
      ipt = tj.Pts.size() - 1 - ii;
      if(tj.Pts[ipt].Chg != 0) {
        tj.EndPt[1] = ipt;
        break;
      }
    }
  } // SetEndPoints
  
  ////////////////////////////////////////////////
  short MCSMom(TjStuff& tjs, Trajectory& tj)
  {
    unsigned short firstPt = tj.EndPt[0];
    unsigned short lastPt = tj.EndPt[1];
    return MCSMom(tjs, tj, firstPt, lastPt);
  } // MCSMom
  
  ////////////////////////////////////////////////
  short MCSMom(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt)
  {
    // Estimate the trajectory momentum using Multiple Coulomb Scattering ala PDG RPP
    
    if(firstPt < tj.EndPt[0]) return 0;
    if(lastPt > tj.EndPt[1]) return 0;
    
    TrajPoint tmp;
    // make a bare trajectory point to define a line between the first
    // and last points on the trajectory that have charge
    MakeBareTrajPoint(tjs, tj.Pts[firstPt], tj.Pts[lastPt], tmp);
    // sum up the deviations^2
    double dsum = 0;
    unsigned short cnt = 0;
    for(unsigned short ipt = firstPt + 1; ipt < lastPt; ++ipt) {
      if(tj.Pts[ipt].Chg == 0) continue;
      dsum += PointTrajDOCA2(tjs, tj.Pts[ipt].HitPos[0],  tj.Pts[ipt].HitPos[1], tmp);
      ++cnt;
    } // ipt
    if(cnt == 0) return 0;
    double sigmaS = sqrt(dsum / (float)cnt);
    double tjLen = TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
    // Theta_o =  4 * sqrt(3) * sigmaS / path
    double thetaRMS = 6.8 * sigmaS / tjLen;
    double mom = 14 * sqrt(tjLen / 14) / thetaRMS;
    if(mom > 999) mom = 999;
    return (short)mom;
  } // MCSMom

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
        if(drTj.PDGCode == 12) continue;
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
        drTj.PDGCode = 12;
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
      bool prt = (muTj.WorkID == debugWorkID);
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
        if(drTj.PDGCode != 12) continue;
        if(drTj.ParentTrajID != muTj.ID) continue;
        // ignore short delta rays
        if(drTj.Pts.size() < minLen) continue;
        float sep0 = 100;
        unsigned short muPt0;
        TrajPointTrajDOCA(tjs, drTj.Pts[drTj.EndPt[0]], muTj, muPt0, sep0);
        float sep1 = 100;
        unsigned short muPt1;
        TrajPointTrajDOCA(tjs, drTj.Pts[drTj.EndPt[1]], muTj, muPt1, sep1);
        if(prt) mf::LogVerbatim("TC")<<" drTj.ID "<<drTj.ID<<" sep 0 "<<sep0<<" sep1 "<<sep1;
        if(sep0 < sep1) { ++n0; } else { ++n1; }
      } // unsigned short jtj
      // Can't tell the direction using this method, so leave the current assignment unchanged
      if(prt) mf::LogVerbatim("TC")<<" n0 "<<n0<<" n1 "<<n1;
      if(n0 == n1) continue;
      if(n0 > n1) {
        // Delta-rays are closer to the beginning (0) end than the end (1) end
        muTj.Dir = 1;
      } else {
        muTj.Dir = -1;
      }
      if(muTj.StepDir < 0) muTj.Dir = -muTj.Dir;
    } // itj
  } // TagMuonDirections
  
  /////////////////////////////////////////
  void MakeBareTrajPoint(TjStuff& tjs, unsigned int fromHit, unsigned int toHit, TrajPoint& tp)
  {
    CTP_t tCTP = EncodeCTP(tjs.fHits[fromHit]->WireID());
    MakeBareTrajPoint(tjs, (float)tjs.fHits[fromHit]->WireID().Wire, tjs.fHits[fromHit]->PeakTime(),
                           (float)tjs.fHits[toHit]->WireID().Wire,   tjs.fHits[toHit]->PeakTime(), tCTP, tp);
    
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  void MakeBareTrajPoint(TjStuff& tjs, float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp)
  {
    tp.CTP = tCTP;
    tp.Pos[0] = fromWire;
    tp.Pos[1] = tjs.UnitsPerTick * fromTick;
    tp.Dir[0] = toWire - fromWire;
    tp.Dir[1] = tjs.UnitsPerTick * (toTick - fromTick);
    float norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  void MakeBareTrajPoint(TjStuff& tjs, TrajPoint& tpIn1, TrajPoint& tpIn2, TrajPoint& tpOut)
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
      tpOut.Pos[0] = -99;
      return;
    }
    tpOut.Dir[0] /= norm;
    tpOut.Dir[1] /= norm;
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
  } // MakeBareTrajPoint
  // ****************************** Printing  ******************************
  
  void PrintAllTraj(std::string someText, TjStuff& tjs, DebugStuff& debug, unsigned short itj, unsigned short ipt)
  {
    
    mf::LogVerbatim myprt("TC");
    
    if(!tjs.vtx3.empty()) {
      // print out 3D vertices
      myprt<<"****** 3D vertices ******************************************__2DVtx_Indx__*******\n";
      myprt<<"Vtx  Cstat  TPC Proc     X       Y       Z    XEr  YEr  ZEr  pln0 pln1 pln2  Wire\n";
      for(unsigned short iv = 0; iv < tjs.vtx3.size(); ++iv) {
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(7)<<tjs.vtx3[iv].CStat;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].TPC;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].ProcCode;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].X;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Y;
        myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Z;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].XErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].YErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].ZErr;
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[0];
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[1];
        myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Ptr2D[2];
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
      myprt<<"Vtx   CTP    wire     error   tick     error  ChiDOF  NCl  topo  traj IDs\n";
      for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
        auto const& aVtx = tjs.vtx[iv];
        if(debug.Plane < 3 && debug.Plane != (int)aVtx.CTP) continue;
        if(aVtx.NTraj == 0) continue;
        myprt<<std::right<<std::setw(3)<<std::fixed<<iv<<std::setprecision(1);
        myprt<<std::right<<std::setw(6)<<aVtx.CTP;
        myprt<<std::right<<std::setw(8)<<aVtx.Wire<<" +/- ";
        myprt<<std::right<<std::setw(4)<<aVtx.WireErr;
        myprt<<std::right<<std::setw(8)<<aVtx.Time/tjs.UnitsPerTick<<" +/- ";
        myprt<<std::right<<std::setw(4)<<aVtx.TimeErr/tjs.UnitsPerTick;
        myprt<<std::right<<std::setw(8)<<aVtx.ChiDOF;
        myprt<<std::right<<std::setw(5)<<aVtx.NTraj;
        myprt<<std::right<<std::setw(6)<<aVtx.Topo;
        myprt<<"    ";
        // display the traj indices
        for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
          auto const& aTj = tjs.allTraj[ii];
          if(debug.Plane < 3 && debug.Plane != (int)aTj.CTP) continue;
          if(aTj.AlgMod[kKilled]) continue;
          for(unsigned short end = 0; end < 2; ++end)
            if(aTj.Vtx[end] == (short)iv) myprt<<std::right<<std::setw(4)<<aTj.ID<<"_"<<end;
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
    unsigned short endPt;
//    unsigned int iht;
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      std::vector<unsigned int> tmp;
      myprt<<someText<<" TRJ  ID CTP Pass Pts frm  to     W:Tick   Ang AveQ     W:T      Ang AveQ ChgRMS  Mom Dir __Vtx__ PDG  Par TRuPDG   EP   KE  WorkID\n";
      for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
        auto const& aTj = tjs.allTraj[ii];
        if(debug.Plane >=0 && debug.Plane < 3 && (unsigned short)debug.Plane != aTj.CTP) continue;
        myprt<<someText<<" ";
        if(aTj.AlgMod[kKilled]) { myprt<<"xxx"; } else { myprt<<"TRJ"; }
        myprt<<std::fixed<<std::setw(4)<<aTj.ID;
        myprt<<std::setw(3)<<aTj.CTP;
        myprt<<std::setw(5)<<aTj.Pass;
        myprt<<std::setw(5)<<aTj.Pts.size();
        myprt<<std::setw(4)<<aTj.EndPt[0];
        myprt<<std::setw(4)<<aTj.EndPt[1];
        endPt = aTj.EndPt[0];
        TrajPoint tp = aTj.Pts[endPt];
        unsigned short itick = tp.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(5)<<(int)tp.AveChg;
        endPt = aTj.EndPt[1];
        tp = aTj.Pts[endPt];
        itick = tp.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
        myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
        myprt<<std::setw(5)<<(int)tp.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<aTj.ChgRMS;
//        myprt<<std::setw(7)<<std::setprecision(2)<<aTj.Quality;
        myprt<<std::setw(5)<<aTj.MCSMom;
        myprt<<std::setw(4)<<aTj.Dir;
        myprt<<std::setw(4)<<aTj.Vtx[0];
        myprt<<std::setw(4)<<aTj.Vtx[1];
        myprt<<std::setw(5)<<aTj.PDGCode;
        myprt<<std::setw(4)<<aTj.ParentTrajID;
        myprt<<std::setw(6)<<aTj.TruPDG;
        myprt<<std::setw(6)<<std::setprecision(2)<<aTj.EffPur;
        myprt<<std::setw(5)<<(int)aTj.TruKE;
        myprt<<std::setw(7)<<aTj.WorkID;
        // print the seed hit that started this trajectory
        if(aTj.StepDir > 0) {
          if(aTj.Pts[0].Chg == 0) myprt<<" "<<PrintPos(tjs, aTj.Pts[0]);
        } else {
          endPt = aTj.EndPt[1];
          if(aTj.Pts[0].Chg == 0) myprt<<" "<<PrintPos(tjs, aTj.Pts[endPt]);
        }
/*
        if(!aTj.Pts[0].Hits.empty()) {
          iht = aTj.Pts[0].Hits[0];
          myprt<<" "<<PrintHit(tjs.fHits[iht]);
        }
*/
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(aTj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
      } // ii
      return;
    } // itj > tjs.allTraj.size()-1
    
    if(itj > tjs.allTraj.size()-1) return;
    
    auto const& aTj = tjs.allTraj[itj];
    
    mf::LogVerbatim("TC")<<"Print tjs.allTraj["<<itj<<"]: ClusterIndex "<<aTj.ClusterIndex<<" Vtx[0] "<<aTj.Vtx[0]<<" Vtx[1] "<<aTj.Vtx[1];
    
    PrintHeader(someText);
    if(ipt == USHRT_MAX) {
      // print all points
      for(unsigned short ii = 0; ii < aTj.Pts.size(); ++ii) PrintTrajPoint(someText, tjs, ii, aTj.StepDir, aTj.Pass, aTj.Pts[ii]);
    } else {
      // print just one
      PrintTrajPoint(someText, tjs, ipt, aTj.StepDir, aTj.Pass, aTj.Pts[ipt]);
    }
  } // Printtjs.allTraj
  
  
  //////////////////////////////////////////
  void PrintTrajectory(std::string someText, TjStuff& tjs, Trajectory const& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    unsigned short first = 0;
    unsigned short last = tj.Pts.size();
    if(tPoint == USHRT_MAX) {
      if(tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ";
        myprt<<"Work:    ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1]<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      } else {
        mf::LogVerbatim myprt("TC");
        myprt<<"tjs.allTraj: ID "<<tj.ID<<" CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.Vtx[0]<<" "<<tj.Vtx[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[0]<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      }
      PrintHeader(someText);
      for(unsigned short ipt = first; ipt < last; ++ipt) PrintTrajPoint(someText, tjs, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
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
    mf::LogVerbatim("TC")<<someText<<" TRP  CTP  Ind  Stp      W:Tick    Delta  RMS    Ang   Err   Kink  Dir0  Dir1      Q    AveQ  Pull FitChi  NTPF  Hits ";
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
    /*
     int itick = (int)(tp.Delta/tjs.UnitsPerTick);
     //    if(itick < 10) myprt<<" "; if(itick < 100) myprt<<" "; if(itick < 1000) myprt<<" ";
     myprt<<std::setw(6)<<std::setprecision(1)<<itick; // dTick
     */
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.DeltaRMS;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Ang;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.AngErr;
    myprt<<std::setw(7)<<std::setprecision(2)<<tp.KinkAng;
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[0];
    myprt<<std::setw(6)<<std::setprecision(2)<<tp.Dir[1];
    myprt<<std::setw(7)<<(int)tp.Chg;
    myprt<<std::setw(8)<<(int)tp.AveChg;
    myprt<<std::setw(6)<<std::setprecision(1)<<tp.ChgPull;
    myprt<<std::setw(7)<<tp.FitChi;
    myprt<<std::setw(6)<<tp.NTPsFit;
    // print the hits associated with this traj point
    for(unsigned short iht = 0; iht < tp.Hits.size(); ++iht) {
      if(tjs.newHits.empty()) {
        // print old hits
        myprt<<" "<<PrintHit(tjs.fHits[tp.Hits[iht]]);
      } else {
        // print new hits
        myprt<<" "<<PrintHit(tjs.newHits[tp.Hits[iht]]);
      }
//      myprt<<" "<<tjs.fHits[tp.Hits[iht]]->WireID().Wire<<":"<<(int)tjs.fHits[tp.Hits[iht]]->PeakTime();
      if(tp.UseHit[iht]) {
        // Distinguish used hits from nearby hits
        myprt<<"_";
      } else {
        myprt<<"x";
      }
      if(tjs.newHits.empty()) { myprt<<tjs.inTraj[tp.Hits[iht]]; } else { myprt<<"NA"; }
    } // iht
  } // PrintTrajPoint
  
  /////////////////////////////////////////
  std::string PrintHit(const art::Ptr<recob::Hit>& hit)
  {
    return std::to_string(hit->WireID().Wire) + ":" + std::to_string((int)hit->PeakTime());
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintHit(const recob::Hit& hit)
  {
    return std::to_string(hit.WireID().Wire) + ":" + std::to_string((int)hit.PeakTime());
  } // PrintHit
  
  /////////////////////////////////////////
  std::string PrintPos(TjStuff& tjs, TrajPoint const& tp)
  {
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    unsigned int time = std::nearbyint(tp.Pos[1]/tjs.UnitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

  
} // namespace tca

