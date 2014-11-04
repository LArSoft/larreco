//////////////////////////////////////////////////////////////////////
///
/// TrackLineFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 3D line given a number of points in 3 wire planes
///
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>
#include <iostream>
#include <iomanip>

#include "RecoAlg/TrackLineFitAlg.h"


namespace trkf{

  TrackLineFitAlg::TrackLineFitAlg() { }

  TrackLineFitAlg::~TrackLineFitAlg() { }


//------------------------------------------------------------------------------
  void TrackLineFitAlg::TrkLineFit(
    std::vector<std::pair<double, geo::WireID>>& hits, double XOrigin,
    TVector3& Pos, TVector3& Dir, float& ChiDOF)
  {
    // Linear fit using X as the independent variable. Hits to be fitted
    // are passed in the hits vector in a pair form (X, WireID). The
    // fitted track position at XOrigin is returned in the Pos vector. 
    // The direction cosines are returned in the Dir vector.
    //
    // SVD fit adapted from $ROOTSYS/tutorials/matrix/solveLinear.C
    // Fit equation is w = A(X)v, where w is a vector of hit wires, A is
    // a matrix to calculate a track projected to a point at X, and v is
    // a vector (Yo, Zo, dY/dX, dZ/dX).
    //
    // Note: The covariance matrix should also be returned
    // B. Baller August 2014

    // assume failure
    ChiDOF = -1;
    
    if(hits.size() < 4) return;
    
    const unsigned int nvars = 4;
    unsigned int npts = hits.size();
  
    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    unsigned short iht, cstat, tpc, ipl;
    double x, cw, sw, off;
    for(iht = 0; iht < hits.size(); ++iht) {
      cstat = hits[iht].second.Cryostat;
      tpc = hits[iht].second.TPC;
      ipl = hits[iht].second.Plane;
      // get the wire plane offset
      off = geom->WireCoordinate(0, 0, ipl, tpc, cstat);
      // get the "cosine-like" component
      cw = geom->WireCoordinate(1, 0, ipl, tpc, cstat) - off;
      // the "sine-like" component
      sw = geom->WireCoordinate(0, 1, ipl, tpc, cstat) - off;
      x = hits[iht].first - XOrigin;
      A[iht][0] = cw;
      A[iht][1] = sw;
      A[iht][2] = cw * x;
      A[iht][3] = sw * x;
      w[iht] = (hits[iht].second.Wire - off);
      ++ninpl[ipl];
      // need at least two points in a plane
      if(ninpl[ipl] == 2) ++nok;
    }
    
    // need at least 2 planes with at least two points
    if(nok < 2) return;
    
    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);

    ChiDOF = 0;
    
    // not enough points to calculate Chisq
    if(hits.size() == 4) return;
    
    double ypr, zpr, diff;
    for(iht = 0; iht < hits.size(); ++iht) {
      cstat = hits[iht].second.Cryostat;
      tpc = hits[iht].second.TPC;
      ipl = hits[iht].second.Plane;
      off = geom->WireCoordinate(0, 0, ipl, tpc, cstat);
      cw = geom->WireCoordinate(1, 0, ipl, tpc, cstat) - off;
      sw = geom->WireCoordinate(0, 1, ipl, tpc, cstat) - off;
      x = hits[iht].first - XOrigin;
      ypr = tVec[0] + tVec[2] * x;
      zpr = tVec[1] + tVec[3] * x;
      diff = ypr * cw + zpr * sw - (hits[iht].second.Wire - off);
      ChiDOF += diff * diff;
    }

    float werr2 = geom->WirePitch() * geom->WirePitch();
    ChiDOF /= werr2;
    ChiDOF /= (float)(npts - 4);
    
    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    Dir[0] = 1 / norm;
    Dir[1] = tVec[2] / norm;
    Dir[2] = tVec[3] / norm;
    
    Pos[0] = XOrigin;
    Pos[1] = tVec[0];
    Pos[2] = tVec[1];
    
  } // TrkLineFit()

} // namespace trkf
