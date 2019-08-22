//////////////////////////////////////////////////////////////////////
///
/// LinFitAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for fitting a 2D line
///
////////////////////////////////////////////////////////////////////////
#ifndef LINFITALG_H
#define LINFITALG_H

#include <vector>


namespace trkf{

  class LinFitAlg {
    public:

    void LinFit(std::vector<float>& x, std::vector<float>& y,
      std::vector<float>& ey2, float& Intercept, float& Slope,
      float& InterceptError, float& SlopeError, float& ChiDOF);

  }; // class LinFitAlg

} // namespace trkf

#endif // ifndef LINFITALG_H
