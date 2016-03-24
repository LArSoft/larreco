////////////////////////////////////////////////////////////////////////
//
//
// TrajClusterAlg Utilities
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGUTILS_H
#define TRAJCLUSTERALGUTILS_H

namespace cluster {
  
  void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& iClosePt, float& Distance);
  // returns the DOCA between a hit and a trajectory
  float PointTrajDOCA(unsigned int iht, TrajPoint const& tp);
  // returns the DOCA between a (W,T) point and a trajectory
  float PointTrajDOCA(float wire, float time, TrajPoint const& tp);
  // returns the DOCA^2 between a point and a trajectory
  float PointTrajDOCA2(float wire, float time, TrajPoint const& tp);
  // returns the separation^2 between two TPs
  float TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2);
  // returns the separation^2 between a point and a TP
  float PointTrajSep2(float wire, float time, TrajPoint const& tp);
  // finds the point on trajectory tj that is closest to trajpoint tp
  void TrajPointTrajDOCA(TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep);
  // returns the intersection position, intPos, of two trajectory points
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
  // Returns the separation distance between two trajectory points
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
  // returns the separation^2 between two hits in WSE units
  float HitSep2(unsigned int iht, unsigned int jht);
  // Returns true if the time separation of two hits exceeds fMultHitSep
  bool LargeHitSep(unsigned int iht, unsigned int jht);
  // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
  void TrajTrajDOCA(Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep);
  // Calculates the angle between two TPs
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2);
  // Put hits in each trajectory point into a flat vector. Only hits with UseHit if onlyUsedHits == true
  void PutTrajHitsInVector(Trajectory const& tj, bool onlyUsedHits, std::vector<unsigned int>& hitVec);
  // Returns a vector with the size of the number of trajectory points on iTj
  // with the separation distance between the closest traj points
//    void TrajSeparation(Trajectory& iTj, Trajectory& jTj, std::vector<float>& tSep);
  // Project TP to a "wire position" Pos[0] and update Pos[1]
  void MoveTPToWire(TrajPoint& tp, float wire);

} // namespace cluster

#endif