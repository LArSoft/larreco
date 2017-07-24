#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {

  /////////////////////////////////////////
  bool TrajPoint3D(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, TVector3& pos, TVector3& dir)
  {
    // Calculate a 3D position and direction from two trajectory points
    
    dir.SetX(999);
    pos = {0, 0, 0};
    
    if(itp.CTP == jtp.CTP) return false;
    
    geo::PlaneID iPlnID = DecodeCTP(itp.CTP);
    geo::PlaneID jPlnID = DecodeCTP(jtp.CTP);

    double ix = tjs.detprop->ConvertTicksToX(itp.Pos[1] / tjs.UnitsPerTick, iPlnID);
    double jx = tjs.detprop->ConvertTicksToX(jtp.Pos[1] / tjs.UnitsPerTick, jPlnID);

    // don't continue if the points are too far apart in X
    if(std::abs(ix - jx) > 10) return false;
    pos[0] = 0.5 * (ix + jx);

    unsigned int iWire = (unsigned int)(itp.Pos[0] + 0.5);
    if(!tjs.geom->HasWire(geo::WireID(iPlnID.Cryostat, iPlnID.TPC, iPlnID.Plane, iWire))) return false;
    unsigned int jWire = (unsigned int)(jtp.Pos[0] + 0.5);
    if(!tjs.geom->HasWire(geo::WireID(jPlnID.Cryostat, jPlnID.TPC, jPlnID.Plane, jWire))) return false;

    // determine the wire orientation and offsets using WireCoordinate
    // wire = yp * OrthY + zp * OrthZ - Wire0 = cs * yp + sn * zp - wire0
    // wire offset
    double iw0 = tjs.geom->WireCoordinate(0, 0, iPlnID);
    // cosine-like component
    double ics = tjs.geom->WireCoordinate(1, 0, iPlnID) - iw0;
    // sine-like component
    double isn = tjs.geom->WireCoordinate(0, 1, iPlnID) - iw0;
    
    double jw0 = tjs.geom->WireCoordinate(0, 0, jPlnID);
    double jcs = tjs.geom->WireCoordinate(1, 0, jPlnID) - jw0;
    double jsn = tjs.geom->WireCoordinate(0, 1, jPlnID) - jw0;
    
    double den = isn * jcs - ics * jsn;
    if(den == 0) return false;
    // Find the Z position of the intersection
    pos[2] = (jcs * (itp.Pos[0] - iw0) - ics * (jtp.Pos[0] - jw0)) / den;
    // and the Y position
    if(ics != 0) {
      pos[1] = (itp.Pos[0] - iw0 - isn * pos[2]) / ics;
    } else {
      pos[1] = (jtp.Pos[0] - jw0 - jsn * pos[2]) / jcs;
    }
    
    // Now find the direction. Protect against large angles first
    if(jtp.Dir[1] == 0) {
      // Going either in the +X direction or -X direction
      if(jtp.Dir[0] > 0) {
        dir.SetX(1);
      } else {
        dir.SetX(-1);
      }
      dir.SetY(0);
      dir.SetZ(0);
      return true;
    } // jtp.Dir[1] == 0
    
    // make a copy of itp and shift it by many wires to avoid precision problems
    TrajPoint itp2 = itp;
    MoveTPToWire(itp2, itp2.Pos[0] + 100);
    // Create a second TVector3 for the shifted point
    TVector3 pos2;
    // Find the X position corresponding to the shifted point 
    pos2[0] = tjs.detprop->ConvertTicksToX(itp2.Pos[1] / tjs.UnitsPerTick, iPlnID);
    // Convert X to Ticks in the j plane and then to WSE units
    double jtp2Pos1 = tjs.detprop->ConvertXToTicks(pos2[0], jPlnID) * tjs.UnitsPerTick;
    // Find the wire position (Pos0) in the j plane that this corresponds to
    double jtp2Pos0 = (jtp2Pos1 - jtp.Pos[1]) * (jtp.Dir[0] / jtp.Dir[1]) + jtp.Pos[0];
    // Find the Y,Z position using itp2 and jtp2Pos0
    pos2[2] = (jcs * (itp2.Pos[0] - iw0) - ics * (jtp2Pos0 - jw0)) / den;
    if(ics != 0) {
      pos2[1] = (itp2.Pos[0] - iw0 - isn * pos2[2]) / ics;
    } else {
      pos2[1] = (jtp2Pos0 - jw0 - jsn * pos2[2]) / jcs;
    }
    dir = pos2 - pos;
    if(dir.Mag() != 0) dir.SetMag(1);
    // Reverse the direction?
    if(dir[0] * itp.Dir[0] < 0) dir *= -1;
    return true;

  } // TrajPoint3D
  
  /////////////////////////////////////////
  void FindMatchingPts(const TjStuff& tjs, const MatchStruct& ms, std::vector<TrajPoint>& stps)
  {
    // Find the points on the list of Tjs matched in MatchStruct ms that best represent a 3D position
    stps.clear();
    if(ms.Count == 0) return;
    
    // find the two Tjs that have the largest extent
    // find the longest Tj
    unsigned short itjLong = USHRT_MAX;
    unsigned short ilen = 0;
    for(unsigned short ii = 0; ii < ms.TjIDs.size(); ++ii) {
      unsigned short itj = ms.TjIDs[ii] - 1;
      const Trajectory& tj = tjs.allTraj[itj];
      if(tj.Pts.size() > ilen) {
        ilen = tj.Pts.size();
        itjLong = itj;
      }
    } // ii
    const Trajectory& itj = tjs.allTraj[itjLong];
    // find the 2nd longest tj
    unsigned short jtjLong = USHRT_MAX;
    unsigned short jlen = 0;
    for(unsigned short jj = 0; jj < ms.TjIDs.size(); ++jj) {
      unsigned short jtj = ms.TjIDs[jj] - 1;
      if(jtj == itjLong) continue;
      const Trajectory& tj = tjs.allTraj[jtj];
      // ensure that this Tj is in a different plane than the longest one
      if(tj.CTP == itj.CTP) continue;
      if(tj.Pts.size() > jlen) {
        jlen = tj.Pts.size();
        jtjLong = jtj;
      }
    } // ii
    // consider both ends of each Tj
    unsigned short ipln = DecodeCTP(itj.CTP).Plane;
    const Trajectory& jtj = tjs.allTraj[jtjLong];
    unsigned short jpln = DecodeCTP(jtj.CTP).Plane;
    unsigned short iend = 0, jend = 0;
    float bestDx = 1E6;
    float bestDxX = 0;
    for(unsigned short ie = 0; ie < 2; ++ie) {
      unsigned short iEndPt = itj.EndPt[ie];
      const TrajPoint& itp = itj.Pts[iEndPt];
      float iX = tjs.detprop->ConvertTicksToX(itp.Pos[1]/tjs.UnitsPerTick, ipln, ms.TPCID.TPC, ms.TPCID.Cryostat);
      for(unsigned short je = 0; je < 2; ++je) {
        unsigned short jEndPt = jtj.EndPt[je];
        const TrajPoint& jtp = jtj.Pts[jEndPt];
        float jX = tjs.detprop->ConvertTicksToX(jtp.Pos[1]/tjs.UnitsPerTick, jpln, ms.TPCID.TPC, ms.TPCID.Cryostat);
        float dx = std::abs(iX - jX);
        if(dx < bestDx) {
          bestDx = dx;
          iend = ie;
          jend = je;
          bestDxX = iX;
        }
      } // e2
    } // e1
    if(bestDx > 5) return;
    unsigned short iEndPt = itj.EndPt[iend];
    TrajPoint tp = itj.Pts[iEndPt];
    // Highjack some TP variables
    tp.AngleCode = iend;
    tp.Delta = itj.ID;
    tp.Hits.clear();
    stps.push_back(tp);
    unsigned short jEndPt = jtj.EndPt[jend];
    tp = jtj.Pts[jEndPt];
    tp.AngleCode = jend;
    tp.Delta = jtj.ID;
    tp.Hits.clear();
    stps.push_back(tp);
    
    if(ms.TjIDs.size() == 2) return;
    
    // now find the matching end of the other Tj
    // find the 2nd longest tj
    for(unsigned short kk = 0; kk < ms.TjIDs.size(); ++kk) {
      unsigned short ktj = ms.TjIDs[kk] - 1;
      if(ktj == itjLong || ktj == jtjLong) continue;
      const Trajectory& tj = tjs.allTraj[ktj];
       // ensure that this Tj is in a different plane
      if(tj.CTP == itj.CTP || tj.CTP == jtj.CTP) continue;
      unsigned short kpln = DecodeCTP(tj.CTP).Plane;
      TrajPoint tp0 = tj.Pts[tj.EndPt[0]];
      float x0 = tjs.detprop->ConvertTicksToX(tp0.Pos[1]/tjs.UnitsPerTick, kpln, ms.TPCID.TPC, ms.TPCID.Cryostat);
      TrajPoint tp1 = tj.Pts[tj.EndPt[1]];
      float x1 = tjs.detprop->ConvertTicksToX(tp1.Pos[1]/tjs.UnitsPerTick, kpln, ms.TPCID.TPC, ms.TPCID.Cryostat);
      if(std::abs(x0 - bestDxX) < std::abs(x1 - bestDxX)) {
        tp0.AngleCode = 0;
        tp0.Delta = tj.ID;
        tp0.Hits.clear();
        stps.push_back(tp0);
      } else {
        tp1.AngleCode = 1;
        tp1.Delta = tj.ID;
        tp1.Hits.clear();
        stps.push_back(tp1);
      }
      break;
    } // ii

  } // FindMatchingPts
  
  /////////////////////////////////////////
  bool CompatibleMerge(TjStuff& tjs, const Trajectory& tj1, const Trajectory& tj2)
  {
    // Returns true if the two trajectories meet some basic requirements for merging
    
    if(tj1.CTP != tj2.CTP) return false;
    if(tj1.AlgMod[kKilled] || tj2.AlgMod[kKilled]) return false;
    // assume that the merge will occur between end1 of tj1 and end0 of tj2
    unsigned short end1 = 1;
    unsigned short end2 = 0;
    const TrajPoint& tp10 = tj1.Pts[tj1.EndPt[0]];
    const TrajPoint& tp11 = tj1.Pts[tj1.EndPt[1]];
    const TrajPoint& tp20 = tj2.Pts[tj2.EndPt[0]];
    const TrajPoint& tp21 = tj2.Pts[tj2.EndPt[1]];
    if(PosSep2(tp10.Pos, tp21.Pos) < PosSep2(tp11.Pos, tp20.Pos)) { end1 = 0; end2 = 1; }
    
    const TrajPoint& tp1 = tj1.Pts[tj1.EndPt[end1]];
    const TrajPoint& tp2 = tj1.Pts[tj1.EndPt[end2]];
    
    if(tp1.AngleCode == 2 && tp2.AngleCode == 2) return true;
    
    if(!SignalBetween(tjs, tp1, tp2, 0.8, false)) return false;
    
    float sep = PosSep2(tp1.Pos, tp2.Pos);
    // require < 5^2 wires
    if(sep > 25) return false;
    // either one short?
    if(tj1.Pts.size() < 5 || tj2.Pts.size() < 5) return false;
    if(sep > 15) return false;
    if(DeltaAngle(tp1.Ang, tp2.Ang) > 0.8) return false;
    
    return true;
    
  } // CompatibleMerge
  
  /////////////////////////////////////////
  void FilldEdx(TjStuff& tjs, MatchStruct& ms)
  {
    // Fills the dEdX vector in the match struct. This function should be called after the
    // matched trajectory points are ordered so that dE/dx is calculated at the start of the PFParticle
    if(ms.Count == 0) return;
    ms.dEdx.resize(tjs.NumPlanes);
    ms.dEdxErr.resize(tjs.NumPlanes);

    double t0 = 0;
    
    unsigned short maxlen = 0;
    for(auto tjID : ms.TjIDs) {
      Trajectory& tj = tjs.allTraj[tjID - 1];
      geo::PlaneID planeID = DecodeCTP(tj.CTP);
      double angleToVert = tjs.geom->WireAngleToVertical(tjs.geom->View(planeID), planeID.TPC, planeID.Cryostat) - 0.5 * ::util::pi<>();
      double cosgamma = std::abs(std::sin(angleToVert) * ms.sDir.Y() + std::cos(angleToVert) * ms.sDir.Z());
      if(cosgamma == 0) continue;
      double dx = tjs.geom->WirePitch(planeID) / cosgamma;
      if(dx == 0) continue;
      double dQ = tj.Pts[tj.EndPt[0]].AveChg;
      if(dQ == 0) continue;
      // convert to dQ/dx
      dQ /= dx;
      double time = tj.Pts[tj.EndPt[0]].Pos[1] / tjs.UnitsPerTick;
      ms.dEdx[planeID.Plane] = tjs.caloAlg->dEdx_AREA(dQ, time, planeID.Plane, t0);
      // Grab the best plane if 1 < dE/dx < 50 MeV/cm
      if(ms.dEdx[planeID.Plane] > 1 && ms.dEdx[planeID.Plane] < 50) {
        if(tj.Pts.size() > maxlen) {
          maxlen = tj.Pts.size();
          ms.BestPlane = planeID.Plane;
        }
      } // valid dE/dx
    } // tj
  } // FilldEdX

  /////////////////////////////////////////
  unsigned short AngleRange(TjStuff& tjs, TrajPoint const& tp)
  {
    return AngleRange(tjs, tp.Ang);
  }
  
  /////////////////////////////////////////
  void SetAngleCode(TjStuff& tjs, TrajPoint& tp)
  {
    unsigned short ar = AngleRange(tjs, tp.Ang);
    if(ar == tjs.AngleRanges.size() - 1) {
      // Very large angle
      tp.AngleCode = 2;
    } else if(tjs.AngleRanges.size() > 2 && ar == tjs.AngleRanges.size() - 2) {
      // Large angle
      tp.AngleCode = 1;
    } else {
      // Small angle
      tp.AngleCode = 0;
    }
    
  } // SetAngleCode
  
  /////////////////////////////////////////
  unsigned short AngleRange(TjStuff& tjs, float angle)
  {
    // returns the index of the angle range
    if(angle > M_PI) angle = M_PI;
    if(angle < -M_PI) angle = M_PI;
    if(angle < 0) angle = -angle;
    if(angle > M_PI/2) angle = M_PI - angle;
    for(unsigned short ir = 0; ir < tjs.AngleRanges.size(); ++ir) {
      if(angle < tjs.AngleRanges[ir]) return ir;
    }
    return tjs.AngleRanges.size() - 1;
  } // AngleRange
  
  //////////////////////////////////////////
  void FitTraj(TjStuff& tjs, Trajectory& tj)
  {
    // Jacket around FitTraj to fit the leading edge of the supplied trajectory
    unsigned short originPt = tj.EndPt[1];
    unsigned short npts = tj.Pts[originPt].NTPsFit;
    TrajPoint tpFit;
    unsigned short fitDir = -1;
    FitTraj(tjs, tj, originPt, npts, fitDir, tpFit);
    tj.Pts[originPt] = tpFit;
    
  } // FitTraj
  
  //////////////////////////////////////////
  void FitTraj(TjStuff& tjs, Trajectory& tj, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint& tpFit)
  {
    // Fit the supplied trajectory using HitPos positions with the origin at originPt.
    // The npts is interpreted as the number of points on each side of the origin
    // The allowed modes are as follows, where i denotes a TP that is included, . denotes
    // a TP with no hits, and x denotes a TP that is not included
    //TP 012345678  fitDir  originPt npts
    //   Oiiixxxxx   1        0       4 << npts in the fit
    //   xi.iiOxxx  -1        5       4
    //   xiiiOiiix   0        4       4 << 2 * npts + 1 points in the fit
    //   xxxiO.ixx   0        4       1
    //   0iiixxxxx   0        0       4
    // This routine puts the results into tp if the fit is successfull. The
    // fit "direction" is in increasing order along the trajectory from 0 to tj.Pts.size() - 1.
    
    //    static const float twoPi = 2 * M_PI;
    
    if(originPt > tj.Pts.size() - 1) {
      mf::LogWarning("TC")<<"FitTraj: Requesting fit of invalid TP "<<originPt;
      return;
    }
    
    // copy the origin TP into the fit TP
    tpFit = tj.Pts[originPt];
    // Assume that the fit will fail
    tpFit.FitChi = 999;
    if(fitDir < -1 || fitDir > 1) return;
    
    std::vector<double> x, y;
    std::array<float, 2> origin = tj.Pts[originPt].HitPos;
    // Use TP position if there aren't any hits on it
    if(tj.Pts[originPt].Chg == 0) origin = tj.Pts[originPt].Pos;
    
    // simple two point case
    if(NumPtsWithCharge(tjs, tj, false) == 2) {
      for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        double xx = tj.Pts[ipt].HitPos[0] - origin[0];
        double yy = tj.Pts[ipt].HitPos[1] - origin[1];
        x.push_back(xx);
        y.push_back(yy);
      } // ii
      if(x.size() != 2) return;
      if(x[0] == x[1]) {
        // Either + or - pi/2
        tpFit.Ang = M_PI/2;
        if(y[1] < y[0]) tpFit.Ang = -tpFit.Ang;
      } else {
        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        tpFit.Ang = atan2(dy, dx);
      }
      tpFit.Dir[0] = cos(tpFit.Ang);
      tpFit.Dir[1] = sin(tpFit.Ang);
      tpFit.Pos[0] += origin[0];
      tpFit.Pos[1] += origin[1];
      tpFit.AngErr = 0.01;
      tpFit.FitChi = 0.01;
      SetAngleCode(tjs, tpFit);
      return;
    } // two points
    
    std::vector<double> w, q;
    std::array<double, 2> dir;
    double xx, yy, xr, yr;
    double chgWt;
    
    // Rotate the traj hit position into the coordinate system defined by the
    // originPt traj point, where x = along the trajectory, y = transverse
    double rotAngle = tj.Pts[originPt].Ang;
    double cs = cos(-rotAngle);
    double sn = sin(-rotAngle);
    
    // enter the originPT hit info if it exists
    if(tj.Pts[originPt].Chg > 0) {
      xx = tj.Pts[originPt].HitPos[0] - origin[0];
      yy = tj.Pts[originPt].HitPos[1] - origin[1];
      xr = cs * xx - sn * yy;
      yr = sn * xx + cs * yy;
      x.push_back(xr);
      y.push_back(yr);
      chgWt = tj.Pts[originPt].ChgPull;
      if(chgWt < 1) chgWt = 1;
      chgWt *= chgWt;
      w.push_back(chgWt * tj.Pts[originPt].HitPosErr2);
    }
    
    // correct npts to account for the origin point
    if(fitDir != 0) --npts;
    
    // step in the + direction first
    if(fitDir != -1) {
      unsigned short cnt = 0;
      for(unsigned short ipt = originPt + 1; ipt < tj.Pts.size(); ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if(chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if(cnt == npts) break;
      } // ipt
    } // fitDir != -1
    
    // step in the - direction next
    if(fitDir != 1 && originPt > 0) {
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = originPt - ii;
        if(ipt > tj.Pts.size() - 1) continue;
        if(tj.Pts[ipt].Chg == 0) continue;
        xx = tj.Pts[ipt].HitPos[0] - origin[0];
        yy = tj.Pts[ipt].HitPos[1] - origin[1];
        xr = cs * xx - sn * yy;
        yr = sn * xx + cs * yy;
        x.push_back(xr);
        y.push_back(yr);
        chgWt = tj.Pts[ipt].ChgPull;
        if(chgWt < 1) chgWt = 1;
        chgWt *= chgWt;
        w.push_back(chgWt * tj.Pts[ipt].HitPosErr2);
        ++cnt;
        if(cnt == npts) break;
        if(ipt == 0) break;
      } // ipt
    } // fitDir != -1
    
    // Not enough points to define a line?
    if(x.size() < 2) return;
    
    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    
    // weight by the charge ratio and accumulate sums
    double wght;
    for(unsigned short ipt = 0; ipt < x.size(); ++ipt) {
      if(w[ipt] < 0.00001) w[ipt] = 0.00001;
      wght = 1 / w[ipt];
      sum   += wght;
      sumx  += wght * x[ipt];
      sumy  += wght * y[ipt];
      sumx2 += wght * x[ipt] * x[ipt];
      sumy2 += wght * y[ipt] * y[ipt];
      sumxy += wght * x[ipt] * y[ipt];
    }
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if(delta == 0) return;
    // A is the intercept
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    // B is the slope
    double B = (sumxy * sum  - sumx * sumy) / delta;
    
    // The chisq will be set below if there are enough points. Don't allow it to be 0
    // so we can take Chisq ratios later
    tpFit.FitChi = 0.01;
    double newang = atan(B);
    dir[0] = cos(newang);
    dir[1] = sin(newang);
    // rotate back into the (w,t) coordinate system
    cs = cos(rotAngle);
    sn = sin(rotAngle);
    tpFit.Dir[0] = cs * dir[0] - sn * dir[1];
    tpFit.Dir[1] = sn * dir[0] + cs * dir[1];
    // ensure that the direction is consistent with the originPt direction
    bool flipDir = false;
    if(AngleRange(tjs, tj.Pts[originPt]) > 0) {
      flipDir = std::signbit(tpFit.Dir[1]) != std::signbit(tj.Pts[originPt].Dir[1]);
    } else {
      flipDir = std::signbit(tpFit.Dir[0]) != std::signbit(tj.Pts[originPt].Dir[0]);
    }
    if(flipDir) {
      tpFit.Dir[0] = -tpFit.Dir[0];
      tpFit.Dir[1] = -tpFit.Dir[1];
    }
    tpFit.Ang = atan2(tpFit.Dir[1], tpFit.Dir[0]);
    SetAngleCode(tjs, tpFit);
    //    if(prt) mf::LogVerbatim("TC")<<"FitTraj "<<originPt<<" originPt Dir "<<tj.Pts[originPt].Dir[0]<<" "<<tj.Pts[originPt].Dir[1]<<" rotAngle "<<rotAngle<<" tpFit.Dir "<<tpFit.Dir[0]<<" "<<tpFit.Dir[1]<<" Ang "<<tpFit.Ang<<" flipDir "<<flipDir<<" fit vector size "<<x.size();
    
    // rotate (0, intcpt) into (W,T) coordinates
    tpFit.Pos[0] = -sn * A + origin[0];
    tpFit.Pos[1] =  cs * A + origin[1];
    // force the origin to be at origin[0]
    if(tpFit.AngleCode < 2) MoveTPToWire(tpFit, origin[0]);
    
    if(x.size() < 3) return;
    
    // Calculate chisq/DOF
    double ndof = x.size() - 2;
    double varnce = (sumy2 + A*A*sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    if(varnce > 0.) {
      // Intercept error is not used
      //      InterceptError = sqrt(varnce * sumx2 / delta);
      double slopeError = sqrt(varnce * sum / delta);
      tpFit.AngErr = std::abs(atan(slopeError));
    } else {
      tpFit.AngErr = 0.01;
    }
    sum = 0;
    // calculate chisq
    double arg;
    for(unsigned short ii = 0; ii < y.size(); ++ii) {
      arg = y[ii] - A - B * x[ii];
      sum += arg * arg / w[ii];
    }
    tpFit.FitChi = sum / ndof;
    
  } // FitTraj

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
    
    // Don't reverse showers
    if(mv.PDGCode == 1111) return false;

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
  unsigned int MatchVecIndex(const TjStuff& tjs, int tjID)
  {
    // returns the index into the tjs.matchVec vector of the first 3D match that
    // includes tjID
    for(unsigned int ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned int ims = tjs.matchVecPFPList[ipfp];
      const MatchStruct& ms = tjs.matchVec[ims];
      if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tjID) != ms.TjIDs.end()) return ims;
    } // indx
    return INT_MAX;
  } // MatchedTjs
  
  ////////////////////////////////////////////////
  unsigned int MatchVecIndex(const TjStuff& tjs, int tjID1, int tjID2)
  {
    // returns the index into the tjs.matchVec vector of the first 3D match that
    // includes both tjID1 and tjID2
    for(unsigned int ipfp = 0; ipfp < tjs.matchVecPFPList.size(); ++ipfp) {
      unsigned int ims = tjs.matchVecPFPList[ipfp];
      const MatchStruct& ms = tjs.matchVec[ims];
      if(std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tjID1) != ms.TjIDs.end() &&
         std::find(ms.TjIDs.begin(), ms.TjIDs.end(), tjID2) != ms.TjIDs.end()) return ims;
    } // indx
    return INT_MAX;
  } // MatchedTjs
  
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

  ////////////////////////////////////////////////
  bool StoreTraj(TjStuff& tjs, Trajectory& tj)
  {
    
    if(tj.EndPt[1] <= tj.EndPt[0]) return false;
    
    if(!(tj.StepDir == 1 || tj.StepDir == -1)) {
      mf::LogError("TC")<<"StoreTraj: Invalid StepDir "<<tj.StepDir;
      return false;
    }
    
    if(tjs.allTraj.size() >= USHRT_MAX) {
      mf::LogError("TC")<<"StoreTraj: Too many trajectories "<<tjs.allTraj.size();
      return false;
    }
    // This shouldn't be necessary but do it anyway
    SetEndPoints(tjs, tj);
    
    // Calculate the charge near the end and beginning if necessary. This must be a short
    // trajectory. Find the average using 4 points
    if(tj.Pts[tj.EndPt[0]].AveChg <= 0) {
      unsigned short cnt = 0;
      float sum = 0;
      for(unsigned short ipt = tj.EndPt[0] + 1; ipt <= tj.EndPt[1]; ++ipt) {
        if(tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == 4) break;
      }
      tj.Pts[tj.EndPt[0]].AveChg = sum / (float)cnt;
    }
    if(tj.Pts[tj.EndPt[1]].AveChg <= 0) {
      float sum = 0;
      unsigned short cnt = 0;
      for(unsigned short ii = 1; ii < tj.Pts.size(); ++ii) {
        unsigned short ipt = tj.EndPt[1] - ii;
        if(tj.Pts[ipt].Chg == 0) continue;
        sum += tj.Pts[ipt].Chg;
        ++cnt;
        if(cnt == 4) break;
        if(ipt == 0) break;
      } // ii
      tj.Pts[tj.EndPt[1]].AveChg = sum / (float)cnt;
    } // begin charge == end charge
    
    short trID = tjs.allTraj.size() + 1;
    for(unsigned short ipt = tj.EndPt[0]; ipt < tj.EndPt[1] + 1; ++ipt) {
      for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
        if(tj.Pts[ipt].UseHit[ii]) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(tjs.fHits[iht].InTraj > 0) {
            mf::LogWarning("TC")<<"StoreTraj: Failed trying to store hit "<<PrintHit(tjs.fHits[iht])<<" in new tjs.allTraj "<<trID<<" but it is used in traj ID = "<<tjs.fHits[iht].InTraj<<" with WorkID "<<tjs.allTraj[tjs.fHits[iht].InTraj-1].WorkID<<" Print and quit";
            PrintTrajectory("SW", tjs, tj, USHRT_MAX);
            ReleaseHits(tjs, tj);
            return false;
          } // error
          tjs.fHits[iht].InTraj = trID;
        }
      } // ii
    } // ipt
    
    // ensure that inTraj is clean for the ID
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].InTraj == tj.ID) {
        mf::LogWarning("TC")<<"StoreTraj: Hit "<<PrintHit(tjs.fHits[iht])<<" thinks it belongs to traj ID "<<tj.ID<<" but it wasn't stored\n";
        PrintTrajectory("SW", tjs, tj, USHRT_MAX);
        return false;
      }
    } // iht
    
    tj.WorkID = tj.ID;
    tj.ID = trID;
    // Calculate the overall charge RMS relative to a linear
    UpdateChgRMS(tjs, tj);
    tjs.allTraj.push_back(tj);
//    if(prt) mf::LogVerbatim("TC")<<"StoreTraj trID "<<trID<<" CTP "<<tj.CTP<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
    if(debug.Hit != UINT_MAX) {
      // print out some debug info
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) {
        for(unsigned short ii = 0; ii < tj.Pts[ipt].Hits.size(); ++ii) {
          unsigned int iht = tj.Pts[ipt].Hits[ii];
          if(iht == debug.Hit) std::cout<<"Debug hit appears in trajectory w WorkID "<<tj.WorkID<<" UseHit "<<tj.Pts[ipt].UseHit[ii]<<"\n";
        } // ii
      } // ipt
    } // debug.Hit ...
    
    return true;
    
  } // StoreTraj
  
  ////////////////////////////////////////////////
  void UpdateChgRMS(TjStuff& tjs, Trajectory& tj)
  {
    // Calculates the ChgRMS variable using all points on the trajectory except a few at the end
    double ave = 0;
    double sum2 = 0;
    double cnt = 0;
    for(short ipt = tj.EndPt[0] + 5; ipt < tj.EndPt[1] - 5; ++ipt) {
      TrajPoint& tp = tj.Pts[ipt];
      if(tp.Chg == 0) continue;
      ave += tp.Chg;
      sum2 += tp.Chg * tp.Chg;
      ++cnt;
    } // tp
    if(cnt < 5) return;
    ave /= cnt;
    sum2 = sum2 - cnt * ave * ave;
    if(sum2 < 0) return;
    tj.ChgRMS = sqrt(sum2 / (cnt - 1));
    tj.ChgRMS /= ave;
  } // UpdateChgRMS
  
  ////////////////////////////////////////////////
  bool InTrajOK(TjStuff& tjs, std::string someText)
  {
    // Check tjs.allTraj -> InTraj associations
    
    unsigned short tID;
    unsigned int iht;
    unsigned short itj = 0;
    std::vector<unsigned int> tHits;
    std::vector<unsigned int> atHits;
    for(auto& tj : tjs.allTraj) {
      // ignore abandoned trajectories
      if(tj.AlgMod[kKilled]) continue;
      tID = tj.ID;
      if(tj.AlgMod[kKilled]) {
        std::cout<<someText<<" ChkInTraj hit size mis-match in tj ID "<<tj.ID<<" AlgBitNames";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) std::cout<<" "<<AlgBitNames[ib];
        std::cout<<"\n";
        continue;
      }
      tHits = PutTrajHitsInVector(tj, kUsedHits);
      if(tHits.size() < 2) {
        std::cout<<someText<<" ChkInTraj: Insufficient hits in traj "<<tj.ID<<"\n";
        PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
        continue;
      }
      std::sort(tHits.begin(), tHits.end());
      atHits.clear();
      for(iht = 0; iht < tjs.fHits.size(); ++iht) {
        if(tjs.fHits[iht].InTraj == tID) atHits.push_back(iht);
      } // iht
      if(atHits.size() < 2) {
        std::cout<<someText<<" ChkInTraj: Insufficient hits in atHits in traj "<<tj.ID<<" Killing it\n";
        tj.AlgMod[kKilled] = true;
        continue;
      }
      if(!std::equal(tHits.begin(), tHits.end(), atHits.begin())) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ChkInTraj failed: inTraj - UseHit mis-match for tj ID "<<tID<<" tj.WorkID "<<tj.WorkID<<" atHits size "<<atHits.size()<<" tHits size "<<tHits.size()<<" in CTP "<<tj.CTP<<"\n";
        myprt<<"AlgMods: ";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
        myprt<<"\n";
        myprt<<"index     inTraj     UseHit \n";
        for(iht = 0; iht < atHits.size(); ++iht) {
          myprt<<"iht "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]]);
          if(iht < tHits.size()) myprt<<" "<<PrintHit(tjs.fHits[tHits[iht]]);
          if(atHits[iht] != tHits[iht]) myprt<<" <<< "<<atHits[iht]<<" != "<<tHits[iht];
          myprt<<"\n";
        } // iht
        if(tHits.size() > atHits.size()) {
          for(iht = atHits.size(); iht < atHits.size(); ++iht) {
            myprt<<"atHits "<<iht<<" "<<PrintHit(tjs.fHits[atHits[iht]])<<"\n";
          } // iht
          PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
        } // tHit.size > atHits.size()
        return false;
      }
      // check the VtxID
      for(unsigned short end = 0; end < 2; ++end) {
        if(tj.VtxID[end] > tjs.vtx.size()) {
          mf::LogVerbatim("TC")<<someText<<" ChkInTraj: Bad VtxID "<<tj.ID;
          std::cout<<someText<<" ChkInTraj: Bad VtxID "<<tj.ID<<" vtx size "<<tjs.vtx.size()<<"\n";
          tj.AlgMod[kKilled] = true;
          PrintTrajectory("CIT", tjs, tj, USHRT_MAX);
          return false;
        }
      } // end
      ++itj;
    } // tj
    return true;
    
  } // InTrajOK

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
    
    if(!tjs.UseAlg[kTEP]) return;
    
    unsigned short minPts = fQualityCuts[1];
    float maxPtSep = minPts + 2;
    if(NumPtsWithCharge(tjs, tj, false) < minPts) return;
    
    if(prt) {
      mf::LogVerbatim("TC")<<"TrimEndPts: minPts "<<minPts<<" required. maxPtSep "<<maxPtSep<<" Minimum hit fraction "<<fQualityCuts[0];
      if(tj.Pts.size() < 50) PrintTrajectory("TEPi", tjs, tj, USHRT_MAX);
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
        // minPts consecutive points have charge. Check the TP Chg fraction
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
    if(nConsecutivePts < minPts - 1 && newEndPt > minPts) {
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
    
    if(newEndPt < minPts) {
      tj.AlgMod[kKilled] = true;
      return;
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
//      if(prt) mf::LogVerbatim("TC")<<" unset "<<ipt;
      UnsetUsedHits(tjs, tj.Pts[ipt]);
    }
    SetEndPoints(tjs, tj);
    tj.Pts.resize(tj.EndPt[1] + 1);
    tj.AlgMod[kTEP] = true;
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
      TrajPoint tp = tp1;
      // check for a signal midway between
      tp.Pos[1] = 0.5 * (tp1.Pos[1] + tp2.Pos[1]);
      if(prt) mf::LogVerbatim("TC")<<" SignalBetween fromWire = toWire = "<<fromWire<<" SignalAtTp? "<<SignalAtTp(tjs, tp);
      return SignalAtTp(tjs, tp);
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
  
  ////////////////////////////////////////////////
  bool TrajHitsOK(TjStuff& tjs, const std::vector<unsigned int>& iHitsInMultiplet, const std::vector<unsigned int>& jHitsInMultiplet)
  {
    // Hits (assume to be on adjacent wires have an acceptable signal overlap
    
    if(iHitsInMultiplet.empty() || jHitsInMultiplet.empty()) return false;
    
    float sum;
    float cvI = HitsPosTick(tjs, iHitsInMultiplet, sum, kAllHits);
    float minI = 1E6;
    float maxI = 0;
    for(auto& iht : iHitsInMultiplet) {
      float cv = tjs.fHits[iht].PeakTime;
      float rms = tjs.fHits[iht].RMS;
      float arg = cv - 3 * rms;
      if(arg < minI) minI = arg;
      arg = cv + 3 * rms;
      if(arg > maxI) maxI = arg;
    }
    
    float cvJ = HitsPosTick(tjs, jHitsInMultiplet, sum, kAllHits);
    float minJ = 1E6;
    float maxJ = 0;
    for(auto& jht : jHitsInMultiplet) {
      float cv = tjs.fHits[jht].PeakTime;
      float rms = tjs.fHits[jht].RMS;
      float arg = cv - 3 * rms;
      if(arg < minJ) minJ = arg;
      arg = cv + 3 * rms;
      if(arg > maxJ) maxJ = arg;
    }
    
    if(cvI < cvJ) {
      if(maxI > minJ) return true;
    } else {
      if(minI < maxJ) return true;
    }
    return false;
  } // TrajHitsOK
  
  /////////////////////////////////////////
  bool TrajHitsOK(TjStuff& tjs, const unsigned int iht, const unsigned int jht)
  {
    // ensure that two adjacent hits have an acceptable overlap
    if(iht > tjs.fHits.size() - 1) return false;
    if(jht > tjs.fHits.size() - 1) return false;
    // require that they be on adjacent wires
    TCHit& ihit = tjs.fHits[iht];
    TCHit& jhit = tjs.fHits[jht];
    unsigned int iwire = ihit.WireID.Wire;
    unsigned int jwire = jhit.WireID.Wire;
    if(abs(iwire - jwire) > 1) return false;
    if(ihit.PeakTime > jhit.PeakTime) {
      float minISignal = ihit.PeakTime - 3 * ihit.RMS;
      float maxJSignal = jhit.PeakTime + 3 * ihit.RMS;
      if(maxJSignal > minISignal) return true;
    } else {
      float maxISignal = ihit.PeakTime + 3 * ihit.RMS;
      float minJSignal = jhit.PeakTime - 3 * ihit.RMS;
      if(minJSignal > maxISignal) return true;
    }
    return false;
  } // TrajHitsOK

  
  /////////////////////////////////////////
  bool SignalAtTp(TjStuff& tjs, const TrajPoint& tp)
  {
    // returns true if there is a hit near tp.Pos
    
    if(tp.Pos[0] < 0) return false;
    unsigned int wire = std::nearbyint(tp.Pos[0]);
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned int ipl = planeID.Plane;
    if(wire >= tjs.NumWires[ipl]) return false;
    if(tp.Pos[1] > tjs.MaxPos1[ipl]) return false;
    // Assume dead wires have a signal
    if(tjs.WireHitRange[ipl][wire].first == -1) return true;
    float projTick = (float)(tp.Pos[1] / tjs.UnitsPerTick);
    // estimate the tick range for non-zero angle
    float tickRange = 0;
    if(std::abs(tp.Dir[1]) != 0) {
      tickRange = std::abs(0.5 / tp.Dir[1]) / tjs.UnitsPerTick;
      // don't let it get too large
      if(tickRange > 40) tickRange = 40;
    }
    float loTpTick = projTick - tickRange;
    float hiTpTick = projTick + tickRange;
    unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
    unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
    
    for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
      TCHit& hit = tjs.fHits[iht];
      if(projTick < hit.PeakTime) {
        float loHitTick = hit.PeakTime - 3 * hit.RMS;
        if(hiTpTick > loHitTick) return true;
      } else {
        float hiHitTick = hit.PeakTime + 3 * hit.RMS;
        if(loTpTick < hiHitTick) return true;
      }
    } // iht
    return false;
    
  } // SignalAtTp

  //////////////////////////////////////////
  float TpSumHitChg(TjStuff& tjs, TrajPoint const& tp){
    float totchg = 0;
    for (size_t i = 0; i<tp.Hits.size(); ++i){
      if (!tp.UseHit[i]) continue;
      totchg += tjs.fHits[tp.Hits[i]].Integral;
    }
    return totchg;
  } // TpSumHitChg

  //////////////////////////////////////////
  bool TjHasNiceVtx(const TjStuff& tjs, const Trajectory& tj, unsigned short minScore)
  {
    // returns true if there is a high-quality vertex at either end
    for(unsigned short end = 0; end < 2; ++end) {
      if(tj.VtxID[end] > 0) {
        unsigned short ivx = tj.VtxID[end] - 1;
        if(tjs.vtx[ivx].Score >= minScore) return true;
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
  bool SplitAllTraj(TjStuff& tjs, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt)
  {
    // Splits the trajectory itj in the tjs.allTraj vector into two trajectories at position pos. Splits
    // the trajectory and associates the ends to the supplied vertex.
    // Here is an example where itj has 9 points and we will split at pos = 4
    // itj (0 1 2 3 4 5 6 7 8) -> new traj (0 1 2 3) + new traj (4 5 6 7 8)
    
    if(itj > tjs.allTraj.size()-1) return false;
    if(pos < tjs.allTraj[itj].EndPt[0] + 1 || pos > tjs.allTraj[itj].EndPt[1] - 1) return false;
    if(ivx != USHRT_MAX && ivx > tjs.vtx.size() - 1) return false;
    
    Trajectory& tj = tjs.allTraj[itj];
    
    if(prt) {
      mf::LogVerbatim myprt("TC");
      myprt<<"SplitAllTraj: Split Tj ID "<<tj.ID<<" at point "<<pos;
      if(ivx < tjs.vtx.size()) myprt<<" with Vtx ID "<<tjs.vtx[ivx].ID;
    }

    // ensure that there will be at least 3 TPs on each trajectory
    unsigned short ipt, ii, ntp = 0;
    for(ipt = 0; ipt < pos; ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) {
      if(prt) mf::LogVerbatim("TC")<<" Split point to small at begin "<<ntp<<" pos "<<pos<<" ID ";
      return false;
    }
    ntp = 0;
    for(ipt = pos + 1; ipt < tj.Pts.size(); ++ipt) {
      if(tj.Pts[ipt].Chg > 0) ++ntp;
      if(ntp > 2) break;
    } // ipt
    if(ntp < 3) {
      if(prt) mf::LogVerbatim("TC")<<" Split point too small at end "<<ntp<<" pos "<<pos<<" EndPt "<<tj.EndPt[1];
      return false;
    }
    
    // make a copy that will become the Tj after the split point
    Trajectory newTj = tjs.allTraj[itj];
    newTj.ID = tjs.allTraj.size() + 1;
    // make another copy in case something goes wrong
    Trajectory oldTj = tjs.allTraj[itj];
    
    
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
    
    // Append 3 points from the end of tj onto the
    // beginning of newTj so that hits can be swapped between
    // them later
    unsigned short eraseSize = pos - 2;
    if(eraseSize > newTj.Pts.size() - 1) {
      tj = oldTj;
      return false;
    }
    
    if(ivx < tjs.vtx.size()) tj.VtxID[1] = tjs.vtx[ivx].ID;
    tj.AlgMod[kSplit] = true;
    if(prt) {
      mf::LogVerbatim("TC")<<" Splitting trajectory ID "<<tj.ID<<" new EndPts "<<tj.EndPt[0]<<" to "<<tj.EndPt[1];
    }
    
    // erase the TPs at the beginning of the new trajectory
    newTj.Pts.erase(newTj.Pts.begin(), newTj.Pts.begin() + eraseSize);
    // unset the first 3 TP hits
    for(ipt = 0; ipt < 3; ++ipt) {
      for(ii = 0; ii < newTj.Pts[ipt].Hits.size(); ++ii) newTj.Pts[ipt].UseHit[ii] = false;
      newTj.Pts[ipt].Chg = 0;
    } // ipt
    SetEndPoints(tjs, newTj);
    if(ivx < tjs.vtx.size()) newTj.VtxID[0] = tjs.vtx[ivx].ID;
    newTj.AlgMod[kSplit] = true;
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
  bool HasDuplicateHits(TjStuff const& tjs, Trajectory const& tj, bool prt)
  {
    // returns true if a hit is associated with more than one TP
    auto tjHits = PutTrajHitsInVector(tj, kAllHits);
    for(unsigned short ii = 0; ii < tjHits.size() - 1; ++ii) {
      for(unsigned short jj = ii + 1; jj < tjHits.size(); ++jj) {
        if(tjHits[ii] == tjHits[jj]) {
          if(prt) mf::LogVerbatim("TC")<<"HDH: Hit "<<PrintHit(tjs.fHits[ii])<<" is a duplicate "<<ii<<" "<<jj;
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
        if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
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
      std::cout<<"FindCloseHits: WireHitRange not valid for CTP "<<tp.CTP<<". tjs.WireHitRange Cstat "<<tjs.TPCID.Cryostat<<" TPC "<<tjs.TPCID.TPC<<"\n";
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
      if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
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
  
  //////////////////////////////////////////
  std::vector<int> FindCloseTjs(const TjStuff& tjs, const TrajPoint& fromTp, const TrajPoint& toTp, const float& maxDelta)
  {
    // Returns a list of Tj IDs that have hits within distance maxDelta on a line drawn between the two Tps as shown
    // graphically here, where a "*" is a Tp and "|" and "-" are the boundaries of the region that is checked 
    //
    //    ---------------
    //    |             |
    //    *             *
    //    |             |
    //    ---------------
    // If the wire positions of fromTp and toTp are the same, a different region is checked as shown here
    //
    //     -----------
    //     |         |
    //     |    *    |
    //     |         |
    //     -----------
    
    std::vector<int> tmp;
    
    TrajPoint tp;
    // Make the tp so that stepping is positive
    unsigned int firstWire, lastWire;
    if(toTp.Pos[0] > fromTp.Pos[0]) {
      if(!MakeBareTrajPoint(tjs, fromTp, toTp, tp)) return tmp;
      firstWire = std::nearbyint(fromTp.Pos[0]);
      lastWire = std::nearbyint(toTp.Pos[0]);
    } else if(toTp.Pos[0] < fromTp.Pos[0]) {
      if(!MakeBareTrajPoint(tjs, toTp, fromTp, tp)) return tmp;
      firstWire = std::nearbyint(toTp.Pos[0]);
      lastWire = std::nearbyint(fromTp.Pos[0]);
    } else {
      tp.Pos = fromTp.Pos;
      float tmp = fromTp.Pos[0] - maxDelta;
      if(tmp < 0) tmp = 0;
      firstWire = std::nearbyint(tmp);
      tmp = fromTp.Pos[0] + maxDelta;
      lastWire = std::nearbyint(tmp);
    }
    
    geo::PlaneID planeID = DecodeCTP(tp.CTP);
    unsigned short ipl = planeID.Plane;
    
    if(firstWire < tjs.FirstWire[ipl]) firstWire = tjs.FirstWire[ipl];
    if(firstWire > tjs.LastWire[ipl]-1) return tmp;
    if(lastWire < tjs.FirstWire[ipl]) return tmp;
    if(lastWire > tjs.LastWire[ipl]-1) lastWire = tjs.LastWire[ipl]-1;
    
    for(unsigned int wire = firstWire; wire <= lastWire; ++wire) {
      if(tjs.WireHitRange[ipl][wire].first == -1) continue;
      if(tjs.WireHitRange[ipl][wire].first == -2) continue;
      MoveTPToWire(tp, (float)wire);
      // Find the tick range at this position
      float minTick = (tp.Pos[1] - maxDelta) / tjs.UnitsPerTick;
      float maxTick = (tp.Pos[1] + maxDelta) / tjs.UnitsPerTick;
      unsigned int firstHit = (unsigned int)tjs.WireHitRange[ipl][wire].first;
      unsigned int lastHit = (unsigned int)tjs.WireHitRange[ipl][wire].second;
      for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
        if(tjs.IgnoreNegChiHits && tjs.fHits[iht].GoodnessOfFit < 0) continue;
        if(tjs.fHits[iht].InTraj <= 0) continue;
        if(tjs.fHits[iht].PeakTime < minTick) continue;
        // Hits are sorted by increasing time so we can break when maxTick is reached
        if(tjs.fHits[iht].PeakTime > maxTick) break;
        if(std::find(tmp.begin(), tmp.end(), tjs.fHits[iht].InTraj) != tmp.end()) continue;
        tmp.push_back(tjs.fHits[iht].InTraj);
      } // iht
    } // wire
    
    return tmp;
    
  } // FindCloseTjs
  
  ////////////////////////////////////////////////
  float ChgFracNearPos(TjStuff& tjs, const std::array<float, 2>& pos, const std::vector<int>& tjIDs)
  {
    // returns the fraction of the charge in the region around pos that is associated with
    // the list of Tj IDs
    if(tjIDs.empty()) return 0;
    std::array<int, 2> wireWindow;
    std::array<float, 2> timeWindow;
    // 1/2 size of the region
    constexpr float NNDelta = 5; 
    wireWindow[0] = pos[0] - NNDelta;
    wireWindow[1] = pos[0] + NNDelta;
    timeWindow[0] = pos[1] - NNDelta;
    timeWindow[1] = pos[1] + NNDelta;
    // do some checking
    for(auto& tjID : tjIDs) if(tjID <= 0 || tjID > (int)tjs.allTraj.size()) return 0;
    // Determine which plane we are in
    geo::PlaneID planeID = DecodeCTP(tjs.allTraj[tjIDs[0]-1].CTP);
    // get a list of all hits in this region
    bool hitsNear;
    std::vector<unsigned int> closeHits = FindCloseHits(tjs, wireWindow, timeWindow, planeID.Plane, kAllHits, true, hitsNear);
    if(closeHits.empty()) return 0;
    float chg = 0;
    float tchg = 0;
    // Add the hit charge in the box
    // All hits in the box, and all hits associated with the Tjs
    for(auto& iht : closeHits) {
      chg += tjs.fHits[iht].Integral;
      if(tjs.fHits[iht].InTraj <= 0) continue;
      if(std::find(tjIDs.begin(), tjIDs.end(), tjs.fHits[iht].InTraj) != tjIDs.end()) tchg += tjs.fHits[iht].Integral;
    } // iht
    if(chg == 0) return 0;
    return tchg / chg;
  } // ChgFracNearPos
  
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
    
    // don't mess with showerTjs
    if(tj.AlgMod[kShowerTj]) return;
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
    
    float tps = TrajPointSeparation(tj.Pts[tj.EndPt[0]], tj.Pts[tj.EndPt[1]]);
    if(tps == 0) return 1;
    
    return MCSThetaRMS(tjs, tj, tj.EndPt[0], tj.EndPt[1]) / sqrt(tps);
    
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
    if(cnt < 3) return 1;
    // require that cnt is a significant fraction of the total number of charged points
    // so that we don't get erroneously high MCSMom when there are large gaps.
    // This is the number of points expected in the count if there are no gaps
    unsigned short numPts = lastPt - firstPt - 1;
    // return the previously calculated value of MCSMom
    if(numPts > 5 && cnt < 0.7 * numPts) return tj.MCSMom;
    double sigmaS = sqrt(dsum / (double)cnt);
    double tjLen = TrajPointSeparation(tj.Pts[firstPt], tj.Pts[lastPt]);
    if(tjLen == 0) return 1;
    // Theta_o =  4 * sqrt(3) * sigmaS / path
    return (6.8 * sigmaS / tjLen);
    
  } // MCSThetaRMS

  /////////////////////////////////////////
  void TagDeltaRays(TjStuff& tjs, const CTP_t& inCTP, short debugWorkID)
  {
    // DeltaRayTag vector elements
    // [0] = max separation of both endpoints from a muon
    // [1] = minimum MCSMom
    // [2] = maximum MCSMom
    
    if(tjs.DeltaRayTag[0] < 0) return;
    if(tjs.DeltaRayTag.size() < 3) return;
    
    float sepCut = tjs.DeltaRayTag[0];
    unsigned short minMom = tjs.DeltaRayTag[1];
    unsigned short maxMom = tjs.DeltaRayTag[2];
    
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
  void TagMuonDirections(TjStuff& tjs, short debugWorkID)
  {
    // Determine muon directions delta-ray proximity to muon trajectories
    
    if(tjs.MuonTag[0] < 0) return;
    
    unsigned short minLen = tjs.MuonTag[3];
    
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

  /////////////////////////////////////////
  bool MakeBareTrajPoint(const TjStuff& tjs, unsigned int fromHit, unsigned int toHit, TrajPoint& tp)
  {
    CTP_t tCTP = EncodeCTP(tjs.fHits[fromHit].WireID);
    return MakeBareTrajPoint(tjs, (float)tjs.fHits[fromHit].WireID.Wire, tjs.fHits[fromHit].PeakTime,
                                  (float)tjs.fHits[toHit].WireID.Wire,   tjs.fHits[toHit].PeakTime, tCTP, tp);
    
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(const TjStuff& tjs, float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp)
  {
    tp.CTP = tCTP;
    tp.Pos[0] = fromWire;
    tp.Pos[1] = tjs.UnitsPerTick * fromTick;
    tp.Dir[0] = toWire - fromWire;
    tp.Dir[1] = tjs.UnitsPerTick * (toTick - fromTick);
    double norm = sqrt(tp.Dir[0] * tp.Dir[0] + tp.Dir[1] * tp.Dir[1]);
    if(norm == 0) return false;
    tp.Dir[0] /= norm;
    tp.Dir[1] /= norm;
    tp.Ang = atan2(tp.Dir[1], tp.Dir[0]);
    return true;
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(const std::array<float, 2>& fromPos, const std::array<float, 2>& toPos, TrajPoint& tpOut)
  {
    tpOut.Pos = fromPos;
    tpOut.Dir[0] = toPos[0] - fromPos[0];
    tpOut.Dir[1] = toPos[1] - fromPos[1];
    double norm = sqrt(tpOut.Dir[0] * tpOut.Dir[0] + tpOut.Dir[1] * tpOut.Dir[1]);
    if(norm == 0) return false;
    tpOut.Dir[0] /= norm;
    tpOut.Dir[1] /= norm;
    tpOut.Ang = atan2(tpOut.Dir[1], tpOut.Dir[0]);
    return true;
    
  } // MakeBareTrajPoint
  
  /////////////////////////////////////////
  bool MakeBareTrajPoint(const TjStuff& tjs, const TrajPoint& tpIn1, const TrajPoint& tpIn2, TrajPoint& tpOut)
  {
    tpOut.CTP = tpIn1.CTP;
    tpOut.Pos = tpIn1.Pos;
    tpOut.Dir[0] = tpIn2.Pos[0] - tpIn1.Pos[0];
    tpOut.Dir[1] = tpIn2.Pos[1] - tpIn1.Pos[1];
    double norm = sqrt(tpOut.Dir[0] * tpOut.Dir[0] + tpOut.Dir[1] * tpOut.Dir[1]);
    if(norm == 0) return false;
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
  unsigned short NumUsedHitsInTj(const TjStuff& tjs, const Trajectory& tj)
  {
    if(tj.AlgMod[kKilled]) return 0;
    unsigned short nhits = 0;
    for(auto& tp : tj.Pts) {
      for(unsigned short ii = 0; ii < tp.Hits.size(); ++ii) if(tp.UseHit[ii]) ++nhits;
    } // tp
    return nhits;
  } // NumHitsInTj
  
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
  
  ////////////////////////////////////////////////
  void SetPDGCode(TjStuff& tjs, unsigned short itj)
  {
    if(itj > tjs.allTraj.size() - 1) return;
    SetPDGCode(tjs, tjs.allTraj[itj]);
  }
  
  ////////////////////////////////////////////////
  void SetPDGCode(TjStuff& tjs, Trajectory& tj)
  {
    // Sets the PDG code for the supplied trajectory
    
    // assume it is unknown
    tj.PDGCode = 0;
    
    tj.MCSMom = MCSMom(tjs, tj);
    
    if(tjs.MuonTag[0] <= 0) return;
/*
    if(tj.AlgMod[kShowerParent]) {
      tj.PDGCode = 11;
      return;
    }
*/
    // Special handling of very long straight trajectories, e.g. uB cosmic rays
    bool isAMuon = (tj.Pts.size() > (unsigned short)tjs.MuonTag[0] && tj.MCSMom > tjs.MuonTag[1]);
    // anything really really long must be a muon
    if(tj.Pts.size() > 200) isAMuon = true;
    if(isAMuon) {
      tj.PDGCode = 13;
//      if(prt) mf::LogVerbatim("TC")<<" SetPDGCode: MCSMom "<<tj.MCSMom<<" PDGCode set to 13";
    }
    
  } // SetPDGCode
  
  
  ////////////////////////////////////////////////
  bool FillWireHitRange(TjStuff& tjs, const geo::TPCID& tpcid, bool debugMode)
  {
    // fills the WireHitRange vector. Slightly modified version of the one in ClusterCrawlerAlg.
    // Returns false if there was a serious error
    
    // determine the number of planes
    geo::TPCGeo const& TPC = tjs.geom->TPC(tpcid);
    unsigned int cstat = tpcid.Cryostat;
    unsigned int tpc = tpcid.TPC;
    unsigned short nplanes = TPC.Nplanes();
    tjs.NumPlanes = nplanes;
    tjs.TPCID = tpcid;
    
    // Y,Z limits of the detector
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &thetpc = tjs.geom->TPC(tpc, cstat);
    thetpc.LocalToWorld(local,world);
    // reduce the active area of the TPC by 1 cm to prevent wire boundary issues
    tjs.XLo = world[0]-tjs.geom->DetHalfWidth(tpc,cstat) + 5;
    tjs.XHi = world[0]+tjs.geom->DetHalfWidth(tpc,cstat) - 5;
    tjs.YLo = world[1]-tjs.geom->DetHalfHeight(tpc,cstat) + 5;
    tjs.YHi = world[1]+tjs.geom->DetHalfHeight(tpc,cstat) - 5;
    tjs.ZLo = world[2]-tjs.geom->DetLength(tpc,cstat)/2 + 5;
    tjs.ZHi = world[2]+tjs.geom->DetLength(tpc,cstat)/2 - 5;
    
    lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    if(!tjs.WireHitRange.empty()) tjs.WireHitRange.clear();
    
    // initialize everything
    tjs.WireHitRange.resize(nplanes);
    tjs.FirstWire.resize(nplanes);
    tjs.LastWire.resize(nplanes);
    tjs.NumWires.resize(nplanes);
    tjs.MaxPos0.resize(nplanes);
    tjs.MaxPos1.resize(nplanes);
    tjs.AveHitRMS.resize(nplanes, nplanes);
    
    std::pair<int, int> flag;
    flag.first = -2; flag.second = -2;
    
    // Calculate tjs.UnitsPerTick, the scale factor to convert a tick into
    // Wire Spacing Equivalent (WSE) units where the wire spacing in this plane = 1.
    // Strictly speaking this factor should be calculated for each plane to handle the
    // case where the wire spacing is different in each plane. Deal with this later if
    // the approximation used here fails.
    
    raw::ChannelID_t channel = tjs.geom->PlaneWireToChannel(0, 0, (int)tpc, (int)cstat);
    float wirePitch = tjs.geom->WirePitch(tjs.geom->View(channel));
    float tickToDist = tjs.detprop->DriftVelocity(tjs.detprop->Efield(),tjs.detprop->Temperature());
    tickToDist *= 1.e-3 * tjs.detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    tjs.UnitsPerTick = tickToDist / wirePitch;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      tjs.FirstWire[ipl] = INT_MAX;
      tjs.LastWire[ipl] = 0;
      tjs.NumWires[ipl] = tjs.geom->Nwires(ipl, tpc, cstat);
      tjs.WireHitRange[ipl].resize(tjs.NumWires[ipl], flag);
      tjs.MaxPos0[ipl] = (float)tjs.NumWires[ipl] - 0.5;
      tjs.MaxPos1[ipl] = (float)tjs.detprop->NumberTimeSamples() * tjs.UnitsPerTick;
    }
    
    // overwrite with the "dead wires" condition
    flag.first = -1; flag.second = -1;
    for(unsigned short ipl = 0; ipl < nplanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        raw::ChannelID_t chan = tjs.geom->PlaneWireToChannel((int)ipl, (int)wire, (int)tpc, (int)cstat);
        if(!channelStatus.IsGood(chan)) tjs.WireHitRange[ipl][wire] = flag;
      } // wire
    } // ipl
    
    unsigned int lastwire = 0, lastipl = 0;
    for(unsigned int iht = 0; iht < tjs.fHits.size(); ++iht) {
      if(tjs.fHits[iht].WireID.Cryostat != cstat) continue;
      if(tjs.fHits[iht].WireID.TPC != tpc) continue;
      unsigned short ipl = tjs.fHits[iht].WireID.Plane;
      unsigned int wire = tjs.fHits[iht].WireID.Wire;
      if(wire > tjs.NumWires[ipl] - 1) {
        mf::LogWarning("TC")<<"FillWireHitRange: Invalid wire number "<<wire<<" > "<<tjs.NumWires[ipl] - 1<<" in plane "<<ipl<<" Quitting";
        return false;
      } // too large wire number
      if(ipl == lastipl && wire < lastwire) {
        mf::LogWarning("TC")<<"FillWireHitRange: Hits are not in increasing wire order. Quitting ";
        return false;
      } // hits out of order
      lastwire = wire;
      lastipl = ipl;
      if(tjs.FirstWire[ipl] == INT_MAX) tjs.FirstWire[ipl] = wire;
      if(tjs.WireHitRange[ipl][wire].first < 0) tjs.WireHitRange[ipl][wire].first = iht;
      tjs.WireHitRange[ipl][wire].second = iht + 1;
      tjs.LastWire[ipl] = wire + 1;
    } // iht
    
    if(!CheckWireHitRange(tjs)) return false;
    
    // Find the average multiplicity 1 hit RMS and calculate the expected max RMS for each range
    if(debugMode && (int)tpc == debug.TPC) std::cout<<"tjs.UnitsPerTick "<<std::setprecision(3)<<tjs.UnitsPerTick<<"\n";
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      float sumRMS = 0;
      float sumAmp = 0;
      unsigned int cnt = 0;
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        // don't let noisy wires screw up the calculation
        if(lastHit - firstHit > 100) continue;
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].Multiplicity != 1) continue;
          if(tjs.fHits[iht].GoodnessOfFit < 0 || tjs.fHits[iht].GoodnessOfFit > 100) continue;
          // don't let a lot of runt hits screw up the calculation
          if(tjs.fHits[iht].PeakAmplitude < 1) continue;
          ++cnt;
          sumRMS += tjs.fHits[iht].RMS;
          sumAmp += tjs.fHits[iht].PeakAmplitude;
        } // iht
      } // wire
      if(cnt < 4) continue;
      tjs.AveHitRMS[ipl] = sumRMS/(float)cnt;
      sumAmp  /= (float)cnt;
      if(debugMode) std::cout<<"Pln "<<ipl<<" tjs.AveHitRMS "<<tjs.AveHitRMS[ipl]<<" Ave PeakAmplitude "<<sumAmp<<"\n";
    } // ipl
    return true;
    
  } // FillWireHitRange
  
  ////////////////////////////////////////////////
  bool CheckWireHitRange(const TjStuff& tjs)
  {
    // do a QC check
    for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
      for(unsigned int wire = 0; wire < tjs.NumWires[ipl]; ++wire) {
        // No hits or dead wire
        if(tjs.WireHitRange[ipl][wire].first < 0) continue;
        unsigned int firstHit = tjs.WireHitRange[ipl][wire].first;
        unsigned int lastHit = tjs.WireHitRange[ipl][wire].second;
        if(lastHit > tjs.fHits.size()) {
          mf::LogWarning("TC")<<"CheckWireHitRange: Invalid lastHit "<<lastHit<<" > fHits.size "<<tjs.fHits.size()<<" in plane "<<ipl;
          std::cout<<"CheckWireHitRange: Invalid lastHit "<<lastHit<<" > fHits.size "<<tjs.fHits.size()<<" in plane "<<ipl<<"\n";
          return false;
        }
        for(unsigned int iht = firstHit; iht < lastHit; ++iht) {
          if(tjs.fHits[iht].WireID.Plane != ipl) {
            mf::LogWarning("TC")<<"CheckWireHitRange: Invalid plane "<<tjs.fHits[iht].WireID.Plane<<" != "<<ipl;
            std::cout<<"CheckWireHitRange: Invalid plane "<<tjs.fHits[iht].WireID.Plane<<" != "<<ipl<<"\n";
            return false;
          }
          if(tjs.fHits[iht].WireID.Wire != wire) {
            mf::LogWarning("TC")<<"CheckWireHitRange: Invalid wire "<<tjs.fHits[iht].WireID.Wire<<" != "<<wire<<" in plane "<<ipl;
            std::cout<<"CheckWireHitRange: Invalid wire "<<tjs.fHits[iht].WireID.Wire<<" != "<<wire<<" in plane "<<ipl<<"\n";
            return false;
          }
        } // iht
      } // wire
    } // ipl
    
    return true;
    
  } // CheckWireHitRange
  
  ////////////////////////////////////////////////
  bool WireHitRangeOK(const TjStuff& tjs, const CTP_t& inCTP)
  {
    // returns true if the passed CTP code is consistent with the CT code of the WireHitRangeVector
    geo::PlaneID planeID = DecodeCTP(inCTP);
    if(planeID.Cryostat != tjs.TPCID.Cryostat) return false;
    if(planeID.TPC != tjs.TPCID.TPC) return false;
    return true;
  }
  
  // ****************************** Printing  ******************************
  
  void PrintAllTraj(std::string someText, const TjStuff& tjs, const DebugStuff& debug, unsigned short itj, unsigned short ipt, bool prtVtx)
  {
    
    mf::LogVerbatim myprt("TC");
    
    if(prtVtx) {
      if(!tjs.vtx3.empty()) {
        // print out 3D vertices
        myprt<<someText<<"****** 3D vertices ******************************************__2DVtx_ID__*******\n";
        myprt<<someText<<"Vtx  Cstat  TPC     X       Y       Z    XEr  YEr  ZEr pln0 pln1 pln2 Wire score nTru  2D_Vtx_Pos\n";
        for(unsigned short iv = 0; iv < tjs.vtx3.size(); ++iv) {
          if(tjs.vtx3[iv].ID == 0) continue;
          myprt<<someText;
          myprt<<std::right<<std::setw(3)<<std::fixed<<tjs.vtx3[iv].ID<<std::setprecision(1);
          myprt<<std::right<<std::setw(7)<<tjs.vtx3[iv].TPCID.Cryostat;
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].TPCID.TPC;
          myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].X;
          myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Y;
          myprt<<std::right<<std::setw(8)<<tjs.vtx3[iv].Z;
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].XErr;
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].YErr;
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].ZErr;
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Vtx2ID[0];
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Vtx2ID[1];
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Vtx2ID[2];
          myprt<<std::right<<std::setw(5)<<tjs.vtx3[iv].Wire;
          // calculate the total score for all planes
          short score = 0;
          unsigned short nTruMatch = 0;
          for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
            if(tjs.vtx3[iv].Vtx2ID[ipl] == 0) continue;
            unsigned short iv2 = tjs.vtx3[iv].Vtx2ID[ipl] - 1;
            score += tjs.vtx[iv2].Score;
            if(tjs.vtx[iv2].Stat[kVtxTruMatch]) ++nTruMatch;
          } // ipl
          myprt<<std::right<<std::setw(6)<<score;
          myprt<<std::right<<std::setw(5)<<nTruMatch;
          if(tjs.vtx3[iv].Wire == -2) {
            myprt<<" PFParticle vertex. No 2D vertices created";
          } else {
            for(unsigned short ipl = 0; ipl < tjs.NumPlanes; ++ipl) {
              if(tjs.vtx3[iv].Vtx2ID[ipl] == 0) {
                myprt<<" NA";
              } else {
                unsigned short ivx = tjs.vtx3[iv].Vtx2ID[ipl] - 1;
                myprt<<" "<<ipl<<":"<<PrintPos(tjs, tjs.vtx[ivx].Pos);
              }
            } // ipl
          }
          myprt<<"\n";
        }
      } // tjs.vtx3.size
      if(!tjs.vtx.empty()) {
        bool foundOne = false;
        for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
          auto& aVtx = tjs.vtx[iv];
          if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aVtx.CTP).Plane) continue;
          if(aVtx.NTraj == 0) continue;
          foundOne = true;
        } // iv
        if(foundOne) {
          // print out 2D vertices
          myprt<<someText<<"************ 2D vertices ************\n";
          myprt<<someText<<"VtxID  CTP   wire  err   tick   err  ChiDOF  NTj Pass  Topo ChgFrac Score  v3D traj_IDs\n";
          for(unsigned short iv = 0; iv < tjs.vtx.size(); ++iv) {
            auto& aVtx = tjs.vtx[iv];
            if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aVtx.CTP).Plane) continue;
            if(aVtx.Stat[kVtxKilled]) continue;
            myprt<<someText;
            myprt<<std::right<<std::setw(3)<<std::fixed<<aVtx.ID;
            myprt<<std::right<<std::setw(6)<<aVtx.CTP;
            myprt<<std::right<<std::setw(8)<<std::setprecision(0)<<std::nearbyint(aVtx.Pos[0]);
            myprt<<std::right<<std::setw(5)<<std::setprecision(1)<<aVtx.PosErr[0];
            myprt<<std::right<<std::setw(8)<<std::setprecision(0)<<std::nearbyint(aVtx.Pos[1]/tjs.UnitsPerTick);
            myprt<<std::right<<std::setw(5)<<std::setprecision(1)<<aVtx.PosErr[1]/tjs.UnitsPerTick;
            myprt<<std::right<<std::setw(7)<<aVtx.ChiDOF;
            myprt<<std::right<<std::setw(5)<<aVtx.NTraj;
            myprt<<std::right<<std::setw(5)<<aVtx.Pass;
            myprt<<std::right<<std::setw(6)<<aVtx.Topo;
            myprt<<std::right<<std::setw(9)<<std::setprecision(2)<<aVtx.TjChgFrac;
            myprt<<std::right<<std::setw(6)<<aVtx.Score;
            myprt<<std::right<<std::setw(5)<<aVtx.Vtx3ID;
            myprt<<"    ";
            // display the traj IDs
            for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
              auto const& aTj = tjs.allTraj[ii];
              if(debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
              if(aTj.AlgMod[kKilled]) continue;
              for(unsigned short end = 0; end < 2; ++end)
                if(aTj.VtxID[end] == (short)aVtx.ID) myprt<<std::right<<std::setw(4)<<aTj.ID<<"_"<<end;
            }
            // Special flags. Ignore the first flag bit (0 = kVtxTrjTried) which is done for every vertex
            for(unsigned short ib = 1; ib < VtxBitNames.size(); ++ib) if(aVtx.Stat[ib]) myprt<<" "<<VtxBitNames[ib];
            myprt<<"\n";
          } // iv
        }
      } // tjs.vtx.size
    }
     
    if(tjs.allTraj.empty()) {
      mf::LogVerbatim("TC")<<someText<<" No allTraj trajectories to print";
      return;
    }
    
    // Print all trajectories in tjs.allTraj if itj == USHRT_MAX
    // Print a single traj (itj) and a single TP (ipt) or all TPs (USHRT_MAX)
    if(itj == USHRT_MAX) {
      // Print summary trajectory information
      std::vector<unsigned int> tmp;
      myprt<<someText<<" TRJ  ID   CTP Pass  Pts frm   to     W:Tick   Ang CS AveQ     W:T      Ang CS AveQ ChgRMS  Mom SDr TDr NN __Vtx__  PDG  Par TRuPDG  E*P TruKE  WorkID \n";
      for(unsigned short ii = 0; ii < tjs.allTraj.size(); ++ii) {
        auto& aTj = tjs.allTraj[ii];
        if(debug.Plane >=0 && debug.Plane < 3 && debug.Plane != (int)DecodeCTP(aTj.CTP).Plane) continue;
        myprt<<someText<<" ";
        if(aTj.AlgMod[kKilled]) { myprt<<"xxx"; } else { myprt<<"TRJ"; }
        myprt<<std::fixed<<std::setw(4)<<aTj.ID;
        myprt<<std::setw(6)<<aTj.CTP;
        myprt<<std::setw(5)<<aTj.Pass;
        myprt<<std::setw(5)<<aTj.Pts.size();
        myprt<<std::setw(4)<<aTj.EndPt[0];
        myprt<<std::setw(5)<<aTj.EndPt[1];
        unsigned short endPt0 = aTj.EndPt[0];
        auto& tp0 = aTj.Pts[endPt0];
        int itick = tp0.Pos[1]/tjs.UnitsPerTick;
        if(itick < 0) itick = 0;
        myprt<<std::setw(6)<<(int)(tp0.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) { myprt<<" "; }
        if(itick < 100) { myprt<<" "; }
        if(itick < 1000) { myprt<<" "; }
        myprt<<std::setw(6)<<std::setprecision(2)<<tp0.Ang;
        myprt<<std::setw(2)<<tp0.AngleCode;
        if(aTj.StopFlag[0][kBragg]) {
          myprt<<"B";
        } else {
          myprt<<" ";
        }
        myprt<<std::setw(5)<<(int)tp0.AveChg;
        unsigned short endPt1 = aTj.EndPt[1];
        auto& tp1 = aTj.Pts[endPt1];
        itick = tp1.Pos[1]/tjs.UnitsPerTick;
        myprt<<std::setw(6)<<(int)(tp1.Pos[0]+0.5)<<":"<<itick; // W:T
        if(itick < 10) { myprt<<" "; }
        if(itick < 100) { myprt<<" "; }
        if(itick < 1000) { myprt<<" "; }
        myprt<<std::setw(6)<<std::setprecision(2)<<tp1.Ang;
        myprt<<std::setw(2)<<tp1.AngleCode;
        if(aTj.StopFlag[1][kBragg]) {
          myprt<<"B";
        } else {
          myprt<<" ";
        }
        myprt<<std::setw(5)<<(int)tp1.AveChg;
        myprt<<std::setw(7)<<std::setprecision(2)<<aTj.ChgRMS;
        myprt<<std::setw(5)<<aTj.MCSMom;
        myprt<<std::setw(4)<<aTj.StepDir;
        myprt<<std::setw(3)<<aTj.TjDir;
        myprt<<std::setw(3)<<aTj.NNeighbors;
        myprt<<std::setw(4)<<aTj.VtxID[0];
        myprt<<std::setw(4)<<aTj.VtxID[1];
        myprt<<std::setw(5)<<aTj.PDGCode;
        myprt<<std::setw(5)<<aTj.ParentTrajID;
        int truKE = 0;
        int pdg = 0;
        if(aTj.MCPartListIndex < tjs.MCPartList.size()) {
          auto& mcp = tjs.MCPartList[aTj.MCPartListIndex];
          truKE = 1000 * (mcp->E() - mcp->Mass());
          pdg = mcp->PdgCode();
        }
        myprt<<std::setw(6)<<pdg;
        myprt<<std::setw(6)<<std::setprecision(2)<<aTj.EffPur;
        myprt<<std::setw(5)<<truKE;
        myprt<<std::setw(7)<<aTj.WorkID;
//        myprt<<" "<<PrintStopFlag(aTj, 0)<<" "<<PrintStopFlag(aTj, 1)<<" ";
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
  void PrintTrajectory(std::string someText, const TjStuff& tjs, const Trajectory& tj, unsigned short tPoint)
  {
    // prints one or all trajectory points on tj
    
    if(tPoint == USHRT_MAX) {
      if(tj.ID < 0) {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ";
        myprt<<"Work:    ID "<<tj.ID<<"    CTP "<<tj.CTP<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.VtxID[0]<<" "<<tj.VtxID[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        myprt<<" MCSMom "<<tj.MCSMom;
        myprt<<" StopFlags "<<PrintStopFlag(tj, 0)<<" "<<PrintStopFlag(tj, 1);
        myprt<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      } else {
        mf::LogVerbatim myprt("TC");
        myprt<<someText<<" ";
        myprt<<"tjs.allTraj: ID "<<tj.ID<<" WorkID "<<tj.WorkID<<" StepDir "<<tj.StepDir<<" PDG "<<tj.PDGCode<<" TruPDG "<<tj.TruPDG<<" tjs.vtx "<<tj.VtxID[0]<<" "<<tj.VtxID[1]<<" nPts "<<tj.Pts.size()<<" EndPts "<<tj.EndPt[0]<<" "<<tj.EndPt[1];
        myprt<<" MCSMom "<<tj.MCSMom;
        myprt<<" StopFlags "<<PrintStopFlag(tj, 0)<<" "<<PrintStopFlag(tj, 1);
        myprt<<" AlgMod names:";
        for(unsigned short ib = 0; ib < AlgBitNames.size(); ++ib) if(tj.AlgMod[ib]) myprt<<" "<<AlgBitNames[ib];
      }
      PrintHeader(someText);
      for(unsigned short ipt = 0; ipt < tj.Pts.size(); ++ipt) PrintTrajPoint(someText, tjs, ipt, tj.StepDir, tj.Pass, tj.Pts[ipt]);
      // See if this trajectory is a shower Tj
      if(tj.AlgMod[kShowerTj]) {
        for(unsigned short ic = 0; ic < tjs.cots.size(); ++ic) {
          if(tjs.cots[ic].TjIDs.empty() == 0) continue;
          // only print out the info for the correct Tj
          if(tjs.cots[ic].ShowerTjID != tj.ID) continue;
          const ShowerStruct& ss = tjs.cots[ic];
          mf::LogVerbatim myprt("TC");
          myprt<<"cots index "<<ic<<" ";
          myprt<<someText<<" Envelope";
          if(ss.Envelope.empty()) {
            myprt<<" NA";
          } else {
            for(auto& vtx : ss.Envelope) myprt<<" "<<(int)vtx[0]<<":"<<(int)(vtx[1]/tjs.UnitsPerTick);
          }
          myprt<<" Energy "<<(int)ss.Energy;
          myprt<<" Area "<<std::fixed<<std::setprecision(1)<<(int)ss.EnvelopeArea<<" ChgDensity "<<ss.ChgDensity;
          myprt<<" StartChg "<<(int)ss.StartChg<<" +/- "<<(int)ss.StartChgErr;
          myprt<<"\nInShower TjIDs";
          for(auto& tjID : ss.TjIDs) {
            myprt<<" "<<tjID;
          } // tjID
          myprt<<"\n";
          myprt<<"NearTjIDs";
          for(auto& tjID : ss.NearTjIDs) {
            myprt<<" "<<tjID;
          } // tjID
          myprt<<"\n";
          myprt<<"MatchedTjIDs";
          for(auto& tjID : ss.MatchedTjIDs) {
            myprt<<" "<<tjID;
          } // tjID
          myprt<<"\n";
          myprt<<"Angle "<<std::fixed<<std::setprecision(2)<<ss.Angle<<" +/- "<<ss.AngleErr;
          myprt<<" AspectRatio "<<std::fixed<<std::setprecision(2)<<ss.AspectRatio;
          myprt<<" DirectionFOM "<<std::fixed<<std::setprecision(2)<<ss.DirectionFOM;
          if(ss.ParentID > 0) {
            myprt<<" Parent Tj "<<ss.ParentID<<" FOM "<<ss.ParentFOM;
          } else {
            myprt<<" No external parent defined";
          }
          myprt<<" TruParentID "<<ss.TruParentID<<"\n";
          if(ss.NeedsUpdate) myprt<<"*********** This shower needs to be updated ***********";
          myprt<<"................................................";
        } // ic
      } // Shower Tj
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
    mf::LogVerbatim("TC")<<someText<<" TRP     CTP  Ind  Stp      W:Tick    Delta  RMS    Ang C   Err  Dir0  Dir1      Q    AveQ  Pull FitChi  NTPF  Hits ";
  } // PrintHeader
  
  ////////////////////////////////////////////////
  void PrintTrajPoint(std::string someText, const TjStuff& tjs, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp)
  {
    mf::LogVerbatim myprt("TC");
    myprt<<someText<<" TRP"<<std::fixed;
    myprt<<pass;
    if(dir > 0) { myprt<<"+"; } else { myprt<<"-"; }
    myprt<<std::setw(6)<<tp.CTP;
    myprt<<std::setw(5)<<ipt;
    myprt<<std::setw(5)<<tp.Step;
    myprt<<std::setw(7)<<std::setprecision(1)<<tp.Pos[0]<<":"<<tp.Pos[1]/tjs.UnitsPerTick; // W:T
    if(tp.Pos[1] < 10) { myprt<<"  "; }
    if(tp.Pos[1] < 100) { myprt<<" "; }
    if(tp.Pos[1] < 1000) { myprt<<" "; }
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
        if(iht > tjs.fHits.size() - 1) {
          std::cout<<"crazy hit "<<iht<<" CTP "<<tp.CTP<<"\n";
          continue;
        }
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
  std::string PrintStopFlag(const Trajectory& tj, unsigned short end)
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
  std::string PrintPos(const TjStuff& tjs, const TrajPoint& tp)
  {
    return PrintPos(tjs, tp.Pos);
  } // PrintPos
  
  /////////////////////////////////////////
  std::string PrintPos(const TjStuff& tjs, const std::array<float, 2>& pos)
  {
    unsigned int wire = std::nearbyint(pos[0]);
    int time = std::nearbyint(pos[1]/tjs.UnitsPerTick);
    return std::to_string(wire) + ":" + std::to_string(time);
  } // PrintPos

  
} // namespace tca

