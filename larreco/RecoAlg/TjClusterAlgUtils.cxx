//////////////////////////////////////////////////////////////////////
///
/// TrajClusterAlg utilities
///
/// Bruce Baller, baller@fnal.gov
///
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TjClusterAlgUtils.h"


//////////////////////////////////////////
void TrajClusterAlg::SplitAllTraj(unsigned short itj, unsigned short pos, unsigned short ivx)
{
  // Splits the trajectory itj in the allTraj vector into two trajectories at position pos. Splits
  // the trajectory and associates the ends to the supplied vertex.
  // Here is an example where itj has 9 points and we will split at pos = 4
  // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
  
  fSplitTrajOK = false;
  if(itj > allTraj.size()-1) return;
  if(pos < allTraj[itj].EndPt[0] + 1 || pos > allTraj[itj].EndPt[1] - 1) return;
  
  Trajectory& tj = allTraj[itj];
  
  // ensure that there will be at least 3 TPs on each trajectory
  unsigned short ipt, ii, ntp = 0;
  for(ipt = 0; ipt < pos; ++ipt) {
    if(tj.Pts[ipt].Chg > 0) ++ntp;
    if(ntp > 2) break;
  } // ipt
  if(ntp < 3) return;
  ntp = 0;
  for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
    if(tj.Pts[ipt].Chg > 0) ++ntp;
    if(ntp > 2) break;
  } // ipt
  if(ntp < 3) return;
  
  // make a copy
  Trajectory newTj = allTraj[itj];
  newTj.ID = allTraj.size() + 1;
  
  // Leave the first section of tj in place. Re-assign the hits
  // to the new trajectory
  unsigned int iht;
  for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
    tj.Pts[ipt].Chg = 0;
    for(ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
      if(!tj.Pts[ipt].UseHit[ii]) continue;
      iht = tj.Pts[ipt].Hits[ii];
      // This shouldn't happen but check anyway
      if(inTraj[iht] != tj.ID) continue;
      inTraj[iht] = newTj.ID;
      tj.Pts[ipt].UseHit[ii] = false;
    } // ii
  } // ipt
  SetEndPoints(tj);
  if(ivx != USHRT_MAX) tj.Vtx[1] = ivx;
  tj.AlgMod[kSplitTraj] = true;
  if(vtxPrt) {
    mf::LogVerbatim("TC")<<"SplitAllTraj: itj "<<tj.ID<<" EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
    PrintTrajectory(allTraj[itj], USHRT_MAX);
  }
  
  // Append 3 points from the end of tj onto the
  // beginning of newTj so that hits can be swapped between
  // them later
  unsigned short eraseSize = pos - 2;
  if(eraseSize > newTj.Pts.size() - 1) {
    mf::LogWarning("TC")<<"SplitAllTraj: Bad erase size ";
    fSplitTrajOK = false;
    return;
  }
  
  // erase the TPs at the beginning of the new trajectory
  //    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + pos + 1);
  newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
  // unset the first 3 TP hits
  for(ipt = 0; ipt < 3; ++ipt) {
    for(ii = 0; ii < newTj.Pts[ipt].UseHit.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
  } // ipt
  SetEndPoints(newTj);
  if(ivx != USHRT_MAX) newTj.Vtx[0] = ivx;
  newTj.AlgMod[kSplitTraj] = true;
  allTraj.push_back(newTj);
  if(vtxPrt) {
    mf::LogVerbatim("TC")<<"SplitAllTraj: NewTj "<<newTj.ID<<" EndPts "<<newTj.EndPt[0]<<" to "<<newTj.EndPt[1];
    PrintTrajectory(newTj, USHRT_MAX);
  }
  fSplitTrajOK = true;
  
} // SplitAllTraj

//////////////////////////////////////////
void TrajClusterAlg::TrajPointTrajDOCA(TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep)
{
  // Finds the point, ipt, on trajectory tj that is closest to trajpoint tp
  float best = minSep * minSep;
  closePt = 0;
  float dw, dt, dp2;
  unsigned short ipt;
  for(ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
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
void TrajClusterAlg::TrajTrajDOCA(Trajectory const& tj1, Trajectory const& tj2, unsigned short& ipt1, unsigned short& ipt2, float& minSep)
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
float TrajClusterAlg::TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2)
{
  // returns the separation^2 between the hit position of two trajectory points
  float dw = tp1.HitPos[0] - tp2.HitPos[0];
  float dt = tp1.HitPos[1] - tp2.HitPos[1];
  return dw * dw + dt * dt;
} // TrajPointHitSep2

//////////////////////////////////////////
float TrajClusterAlg::HitSep2(unsigned int iht, unsigned int jht)
{
  // returns the separation^2 between two hits in WSE units
  if(iht > fHits.size()-1 || jht > fHits.size()-1) return 1E6;
  float dw = (float)fHits[iht]->WireID().Wire - (float)fHits[jht]->WireID().Wire;
  float dt = (fHits[iht]->PeakTime() - fHits[jht]->PeakTime()) * fScaleF;
  return dw * dw + dt * dt;
} // TrajPointHitSep2

//////////////////////////////////////////
float TrajClusterAlg::PointTrajSep2(float wire, float time, TrajPoint const& tp)
{
  float dw = wire - tp.Pos[0];
  float dt = time - tp.Pos[1];
  return dw * dw + dt * dt;
}

//////////////////////////////////////////
float TrajClusterAlg::PointTrajDOCA(unsigned int iht, TrajPoint const& tp)
{
  float wire = fHits[iht]->WireID().Wire;
  float time = fHits[iht]->PeakTime() * fScaleF;
  return sqrt(PointTrajDOCA2(wire, time, tp));
} // PointTrajDOCA

//////////////////////////////////////////
float TrajClusterAlg::PointTrajDOCA(float wire, float time, TrajPoint const& tp)
{
  return sqrt(PointTrajDOCA2(wire, time, tp));
} // PointTrajDOCA

//////////////////////////////////////////
float TrajClusterAlg::PointTrajDOCA2(float wire, float time, TrajPoint const& tp)
{
  // returns the distance of closest approach squared between a (wire, time(WSE)) point
  // and a trajectory point
  
  float t = (wire  - tp.Pos[0]) * tp.Dir[0] + (time - tp.Pos[1]) * tp.Dir[1];
  float dw = tp.Pos[0] + t * tp.Dir[0] - wire;
  float dt = tp.Pos[1] + t * tp.Dir[1] - time;
  return (dw * dw + dt * dt);
  
} // PointTrajDOCA2

//////////////////////////////////////////
void TrajClusterAlg::TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y)
{
  // returns the intersection position, intPos, of two trajectory points
  
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
float TrajClusterAlg::TrajLength(Trajectory& tj)
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
float TrajClusterAlg::TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2)
{
  // Returns the separation distance between two trajectory points
  float dx = tp1.Pos[0] - tp2.Pos[0];
  float dy = tp1.Pos[1] - tp2.Pos[1];
  return sqrt(dx * dx + dy * dy);
} // TrajPointSeparation

//////////////////////////////////////////
void TrajClusterAlg::TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& Distance)
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

