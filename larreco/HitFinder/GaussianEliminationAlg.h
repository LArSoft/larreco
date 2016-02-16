#ifndef GAUSSIANELIMINATIONALG_H
#define GAUSSIANELIMINATIONALG_H

/*!
 * Title:   GaussianEliminationAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Class that solves system of linear equations via Gaussian Elimination.
 * Intended for use with RFFHitFitter
 *
*/

#include <vector>

namespace util{

  const float SQRT_TWO_PI = 2.506628;

  class GaussianEliminationAlg {

  public:
    GaussianEliminationAlg(float, float);

    double GetDistance(float) const;
    const std::vector<float>& SolveEquations(const std::vector<float>& meanVector,
					     const std::vector<float>& sigmaVector,
					     const std::vector<float>& heightVector);

    void FillAugmentedMatrix(const std::vector<float>& meanVector,
			     const std::vector<float>& sigmaVector,
			     const std::vector<float>& heightVector);			     
    void GaussianElimination();
    const std::vector<float>& GetSolutions() { return fSolutions; }
    void Print();
    
  private:

    float fDistanceStepSize;
    float fDistanceMax;
    std::vector<double> fDistanceLookupTable;

    void FillDistanceLookupTable();

    std::vector< std::vector<double> > fMatrix;
    std::vector<float>                 fSolutions;

  };

}

#endif
