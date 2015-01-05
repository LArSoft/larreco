//////////////////////////////////////////////////////////////////////
///
/// LinFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 2D line
///
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "RecoAlg/LinFitAlg.h"


namespace trkf{

  LinFitAlg::LinFitAlg() { }

  LinFitAlg::~LinFitAlg() { }


  void LinFitAlg::LinFit(std::vector<float>& x, std::vector<float>& y, 
    std::vector<float>& ey2, float& Intercept, float& Slope, 
    float& InterceptError, float& SlopeError, float& ChiDOF) 
  {
    // fit a line ala Bevington linfit.F. The number of points fit is defined by
    // the size of the y vector. 

    ChiDOF = 999.;

    if(y.size() < 2) return;
    if(x.size() < y.size() || ey2.size() < y.size()) return;
    
    float sum = 0.;
    float sumx = 0.;
    float sumy = 0.;
    float sumxy = 0.;
    float sumx2 = 0.;
    float sumy2 = 0.;
    float weight, arg;
    unsigned short ii;

    for(ii = 0; ii < y.size(); ++ii) {
      weight = 1. / ey2[ii];
      sum += weight;
      sumx += weight * x[ii];
      sumy += weight * y[ii];
      sumx2 += weight * x[ii] * x[ii];
      sumxy += weight * x[ii] * y[ii];
      sumy2 += weight * y[ii] * y[ii];
    }
    // calculate coefficients and std dev
    float delta = sum * sumx2 - sumx * sumx;
    if(delta == 0.) return;
    float A = (sumx2 * sumy - sumx * sumxy) / delta;
    float B = (sumxy * sum  - sumx * sumy) / delta;
    Intercept = A;
    Slope = B;
    if(x.size() == 2) {
      ChiDOF = 0.;
      return;
    }
    float ndof = x.size() - 2;
    float varnce = (sumy2 + A*A*sum + B*B*sumx2 - 
                    2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
    if(varnce > 0.) {
      InterceptError = sqrt(varnce * sumx2 / delta);
      SlopeError = sqrt(varnce * sum / delta);
    } else {
      InterceptError = 0.;
      SlopeError = 0.;
    }
    sum = 0.;
    // calculate chisq
    for(ii = 0; ii < y.size(); ++ii) {
      arg = y[ii] - A - B * x[ii];
      sum += arg * arg / ey2[ii];
    }
    ChiDOF = sum / ndof;
  } // LinFit
  
} // namespace trkf
