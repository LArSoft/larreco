//////////////////////////////////////////////////////////////////////
///
/// LinFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 2D line
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/LinFitAlg.h"
#include <math.h>

namespace trkf {

  void LinFitAlg::LinFit(std::vector<float>& x,
                         std::vector<float>& y,
                         std::vector<float>& ey2,
                         float& Intercept,
                         float& Slope,
                         float& InterceptError,
                         float& SlopeError,
                         float& ChiDOF) const
  {
    // fit a line ala Bevington linfit.F. The number of points fit is defined by
    // the size of the y vector.

    ChiDOF = 999.;

    if (y.size() < 2) return;
    if (x.size() < y.size() || ey2.size() < y.size()) return;

    double sum = 0.;
    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumx2 = 0.;
    double sumy2 = 0.;
    double weight, arg;
    unsigned short ii;

    for (ii = 0; ii < y.size(); ++ii) {
      weight = 1. / ey2[ii];
      sum += weight;
      sumx += weight * x[ii];
      sumy += weight * y[ii];
      sumx2 += weight * x[ii] * x[ii];
      sumxy += weight * x[ii] * y[ii];
      sumy2 += weight * y[ii] * y[ii];
    }
    // calculate coefficients and std dev
    double delta = sum * sumx2 - sumx * sumx;
    if (delta == 0.) return;
    double A = (sumx2 * sumy - sumx * sumxy) / delta;
    double B = (sumxy * sum - sumx * sumy) / delta;
    Intercept = A;
    Slope = B;
    if (x.size() == 2) {
      ChiDOF = 0.;
      return;
    }
    double ndof = x.size() - 2;
    double varnce =
      (sumy2 + A * A * sum + B * B * sumx2 - 2 * (A * sumy + B * sumxy - A * B * sumx)) / ndof;
    if (varnce > 0.) {
      InterceptError = sqrt(varnce * sumx2 / delta);
      SlopeError = sqrt(varnce * sum / delta);
    }
    else {
      InterceptError = 0.;
      SlopeError = 0.;
    }
    sum = 0.;
    // calculate chisq
    for (ii = 0; ii < y.size(); ++ii) {
      arg = y[ii] - A - B * x[ii];
      sum += arg * arg / ey2[ii];
    }
    ChiDOF = sum / ndof;
  } // LinFit

} // namespace trkf
