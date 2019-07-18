////////////////////////////////////////////////////////////////////////
// Class:       MarqFitAlg
// Purpose:     Fit gaussians
//              
//              
// Original code by Mike Wang, converted to a larsoft algorithm by S. Berkman
////////////////////////////////////////////////////////////////////////

#ifndef MARQFITALG_H
#define MARQFITALG_H

#include <iostream>
#include <cmath>
#include <vector>

//#include "fhiclcpp/ParameterSet.h"
//#include "lardataobj/RecoBase/Hit.h"

namespace gshf{

  class MarqFitAlg {
    public:
      explicit MarqFitAlg();
      virtual ~MarqFitAlg() {}

      int cal_perr(float p[], float y[], const int nParam, const int nData, float perr[]);
      int mrqdtfit(float &lambda, float p[], float y[], const int nParam, const int nData, float &chiSqr, float &dchiSqr);


    private:
      //these functions are  called by the public functions
      void fgauss(const float yd[], const float p[], const int npar, const int ndat, std::vector<float> &res);
      void dgauss(const float p[], const int npar, const int ndat, std::vector<float> &dydp);
      float cal_xi2(const std::vector<float> &res, const int ndat);
      void setup_matrix(const std::vector<float> &res, const std::vector<float> &dydp, const int npar, const int ndat, std::vector<float> &beta, std::vector<float> &alpha);
      void solve_matrix(const std::vector<float> &beta, const std::vector<float> &alpha, const int npar, std::vector<float> &dp);
      float invrt_matrix(std::vector<float> &alphaf, const int npar);

  };

}//end namespace gshf
#endif
