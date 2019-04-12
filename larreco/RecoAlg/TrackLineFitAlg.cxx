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
#include <iomanip>

#include "larreco/RecoAlg/TrackLineFitAlg.h"


namespace trkf{

  TrackLineFitAlg::TrackLineFitAlg() { }

  TrackLineFitAlg::~TrackLineFitAlg() { }


//------------------------------------------------------------------------------
  void TrackLineFitAlg::TrkLineFit(std::vector<geo::WireID>& hitWID, std::vector<double>& hitX, std::vector<double>& hitXErr,
                                   double XOrigin, TVector3& Pos, TVector3& Dir, float& ChiDOF)
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

    // assume failure
    ChiDOF = -1;

    if(hitX.size() < 4) return;
    if(hitX.size() != hitWID.size()) return;
    if(hitX.size() != hitXErr.size()) return;

    const unsigned int nvars = 4;
    unsigned int npts = hitX.size();

    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    unsigned short iht;
    unsigned int ipl, tpc, cstat;
    double x, cw, sw, off, wght;
    for(iht = 0; iht < hitX.size(); ++iht) {
      cstat = hitWID[iht].Cryostat;
      tpc = hitWID[iht].TPC;
      ipl = hitWID[iht].Plane;
      // get the wire plane offset
      off = geom->WireCoordinate(0, 0, ipl, tpc, cstat);
      // get the "cosine-like" component
      cw = geom->WireCoordinate(1, 0, ipl, tpc, cstat) - off;
      // the "sine-like" component
      sw = geom->WireCoordinate(0, 1, ipl, tpc, cstat) - off;
      x = hitX[iht] - XOrigin;
      if(hitXErr[iht] > 0) {
        wght = 1 / hitXErr[iht];
      } else {
        wght = 1;
      }
      A[iht][0] = wght * cw;
      A[iht][1] = wght * sw;
      A[iht][2] = wght * cw * x;
      A[iht][3] = wght * sw * x;
      w[iht] = wght * (hitWID[iht].Wire - off);
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
    if(hitX.size() == 4) return;

    double ypr, zpr, diff;
    for(iht = 0; iht < hitX.size(); ++iht) {
      cstat = hitWID[iht].Cryostat;
      tpc = hitWID[iht].TPC;
      ipl = hitWID[iht].Plane;
      off = geom->WireCoordinate(0, 0, ipl, tpc, cstat);
      cw = geom->WireCoordinate(1, 0, ipl, tpc, cstat) - off;
      sw = geom->WireCoordinate(0, 1, ipl, tpc, cstat) - off;
      x = hitX[iht] - XOrigin;
      ypr = tVec[0] + tVec[2] * x;
      zpr = tVec[1] + tVec[3] * x;
      if(hitXErr[iht] > 0) {
        wght = 1 / hitXErr[iht];
      } else {
        wght = 1;
      }
      if(wght <= 0) wght = 1;
      diff = (ypr * cw + zpr * sw - (hitWID[iht].Wire - off)) / wght;
      ChiDOF += diff * diff;
    }

    double werr2 = geom->WirePitch(0, tpc, cstat);
    werr2 *= werr2;
    ChiDOF /= werr2;
    ChiDOF /= (double)(npts - 4);

    double norm = sqrt(1 + tVec[2] * tVec[2] + tVec[3] * tVec[3]);
    Dir[0] = 1 / norm;
    Dir[1] = tVec[2] / norm;
    Dir[2] = tVec[3] / norm;

    Pos[0] = XOrigin;
    Pos[1] = tVec[0];
    Pos[2] = tVec[1];
/*
    // covariance matrix
    TMatrixD fV = svd.GetV();
    PosCov(1, 1) = fV(0, 0);
    PosCov(2, 1) = fV(1, 0);
    PosCov(1, 2) = fV(0, 1);
    PosCov(2, 2) = fV(1, 1);
    // A conservative fake for Pos[0]
    PosCov(0, 0) = PosCov(1,1);
*/

  } // TrkLineFit()

} // namespace trkf
