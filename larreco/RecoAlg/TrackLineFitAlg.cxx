//////////////////////////////////////////////////////////////////////
///
/// TrackLineFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 3D line given a number of points in 3 wire planes
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrackLineFitAlg.h"

#include <math.h>

#include "TDecompSVD.h"
#include "TMatrixDfwd.h"
#include "TMatrixT.h"
#include "TMatrixTUtils.h"
#include "TVector3.h"
#include "TVectorDfwd.h"
#include "TVectorT.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace trkf {

  //------------------------------------------------------------------------------
  void TrackLineFitAlg::TrkLineFit(std::vector<geo::WireID>& hitWID,
                                   std::vector<double>& hitX,
                                   std::vector<double>& hitXErr,
                                   double XOrigin,
                                   TVector3& Pos,
                                   TVector3& Dir,
                                   float& ChiDOF) const
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

    if (hitX.size() < 4) return;
    if (hitX.size() != hitWID.size()) return;
    if (hitX.size() != hitXErr.size()) return;

    const unsigned int nvars = 4;
    unsigned int npts = hitX.size();

    TMatrixD A(npts, nvars);
    // vector holding the Wire number
    TVectorD w(npts);
    unsigned short ninpl[3] = {0};
    unsigned short nok = 0;
    for (std::size_t iht = 0; iht < hitX.size(); ++iht) {
      auto const& wid = hitWID[iht];
      // get the wire plane offset
      double const off = geom->WireCoordinate(geo::Point_t{0, 0, 0}, wid);
      // get the "cosine-like" component
      double const cw = geom->WireCoordinate(geo::Point_t{0, 1, 0}, wid) - off;
      // the "sine-like" component
      double const sw = geom->WireCoordinate(geo::Point_t{0, 0, 1}, wid) - off;
      double const x = hitX[iht] - XOrigin;
      double wght{1.};
      if (hitXErr[iht] > 0) { wght = 1 / hitXErr[iht]; }
      A[iht][0] = wght * cw;
      A[iht][1] = wght * sw;
      A[iht][2] = wght * cw * x;
      A[iht][3] = wght * sw * x;
      w[iht] = wght * (hitWID[iht].Wire - off);
      unsigned int ipl = wid.Plane;
      ++ninpl[ipl];
      // need at least two points in a plane
      if (ninpl[ipl] == 2) ++nok;
    }

    // need at least 2 planes with at least two points
    if (nok < 2) return;

    TDecompSVD svd(A);
    bool ok;
    TVectorD tVec = svd.Solve(w, ok);

    ChiDOF = 0;

    // not enough points to calculate Chisq
    if (hitX.size() == 4) return;

    // FIXME: The 'tpc' and 'cstat' variables are suspect as they are
    //        updated for each iteration below, but only the last
    //        value is used in the geom->WirePitch(...) calculation
    //        below.
    unsigned int tpc{-1u}, cstat{-1u};
    for (std::size_t iht = 0; iht < hitX.size(); ++iht) {
      auto const& wid = hitWID[iht];
      tpc = wid.TPC;
      cstat = wid.Cryostat;
      double const off = geom->WireCoordinate(geo::Point_t{0, 0, 0}, wid);
      double const cw = geom->WireCoordinate(geo::Point_t{0, 1, 0}, wid) - off;
      double const sw = geom->WireCoordinate(geo::Point_t{0, 0, 1}, wid) - off;
      double const x = hitX[iht] - XOrigin;
      double const ypr = tVec[0] + tVec[2] * x;
      double const zpr = tVec[1] + tVec[3] * x;
      double wght{1.};
      if (hitXErr[iht] > 0) { wght = 1 / hitXErr[iht]; }
      assert(wght > 0.);
      double const diff = (ypr * cw + zpr * sw - (hitWID[iht].Wire - off)) / wght;
      ChiDOF += diff * diff;
    }

    double werr2 = geom->WirePitch(geo::PlaneID{cstat, tpc, 0});
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
  } // TrkLineFit()

} // namespace trkf
