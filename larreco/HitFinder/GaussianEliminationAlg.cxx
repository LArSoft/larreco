/*!
 * Title:   GaussianEliminationAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class that solves system of linear equations via Gaussian Elimination.
 * Intended for use with RFFHitFitter
 *
*/

#include "GaussianEliminationAlg.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

util::GaussianEliminationAlg::GaussianEliminationAlg(float step, float max)
{
  fDistanceStepSize = step;
  fDistanceMax = max;

  if (fDistanceStepSize < 0 || fDistanceMax < 0)
    throw std::runtime_error(
      "Error in GaussianEliminationAlg: Cannot construct with negative step or max.");

  FillDistanceLookupTable();
}

void util::GaussianEliminationAlg::FillDistanceLookupTable()
{

  fDistanceLookupTable.clear();
  fDistanceLookupTable.reserve(std::ceil(fDistanceMax / fDistanceStepSize) + 2);

  //let's be extra safe, and make sure we have more sampling points than we will need
  double x_val = 0.0;
  while (x_val <= fDistanceMax + fDistanceStepSize) {
    fDistanceLookupTable.push_back(std::exp(x_val * x_val * 0.5 * -1));
    x_val += fDistanceStepSize;
  }
  //do one more to be sure to push beyond...
  fDistanceLookupTable.push_back(std::exp(x_val * x_val * 0.5 * -1));
}

double util::GaussianEliminationAlg::GetDistance(float d) const
{
  double d_abs = std::abs(d);
  if (d_abs > fDistanceMax) return 0.0;

  size_t low_bin = std::floor(d_abs / fDistanceStepSize);
  return fDistanceLookupTable[low_bin] -
         (d_abs / fDistanceStepSize - (double)low_bin) *
           (fDistanceLookupTable[low_bin] - fDistanceLookupTable[low_bin + 1]);
}

const std::vector<float>& util::GaussianEliminationAlg::SolveEquations(
  const std::vector<float>& meanVector,
  const std::vector<float>& sigmaVector,
  const std::vector<float>& heightVector)
{

  if (meanVector.size() != sigmaVector.size() || meanVector.size() != heightVector.size())
    throw std::runtime_error("Error in GaussianEliminationAlg: Vector sizes not equal.");

  FillAugmentedMatrix(meanVector, sigmaVector, heightVector);
  GaussianElimination();

  return fSolutions;
}

void util::GaussianEliminationAlg::FillAugmentedMatrix(const std::vector<float>& meanVector,
                                                       const std::vector<float>& sigmaVector,
                                                       const std::vector<float>& heightVector)
{

  fMatrix.resize(meanVector.size());
  for (size_t i = 0; i < meanVector.size(); i++) {

    fMatrix[i].resize(meanVector.size() + 1);
    for (size_t j = 0; j < meanVector.size(); j++) {
      if (sigmaVector[j] < std::numeric_limits<float>::epsilon()) {
        if (i == j)
          fMatrix[i][j] = 1.0;
        else
          fMatrix[i][j] = 0.0;
      }
      else
        fMatrix[i][j] = GetDistance((meanVector[i] - meanVector[j]) / sigmaVector[j]);
    }
    fMatrix[i][meanVector.size()] = heightVector[i];
  }
}

void util::GaussianEliminationAlg::GaussianElimination()
{

  fSolutions.resize(fMatrix.size(), 0.0);

  for (size_t i = 0; i < fMatrix.size(); i++) {

    for (size_t j = i + 1; j < (fMatrix[i].size() - 1); j++) {
      float scale_value = fMatrix[j][i] / fMatrix[i][i];
      for (size_t k = i; k < fMatrix[i].size(); k++)
        fMatrix[j][k] -= fMatrix[i][k] * scale_value;
    } //end column loop

  } //end row loop

  for (int i = fMatrix.size() - 1; i >= 0; i--) {
    fSolutions[i] = fMatrix[i].back();

    for (size_t j = i + 1; j < fMatrix.size(); j++)
      fSolutions[i] -= fMatrix[i][j] * fSolutions[j];

    fSolutions[i] /= fMatrix[i][i];
  }
}

void util::GaussianEliminationAlg::Print()
{
  std::cout << "GaussianEliminationAlg." << std::endl;

  std::cout << "\tLookup table (step=" << fDistanceStepSize << ", max=" << fDistanceMax << ")"
            << std::endl;
  for (size_t i = 0; i < fDistanceLookupTable.size(); i++)
    std::cout << "\t\tGaussian(" << fDistanceStepSize * i << ") = " << fDistanceLookupTable[i]
              << std::endl;

  std::cout << "\tAugmented matrix " << std::endl;
  for (size_t i = 0; i < fMatrix.size(); i++) {
    std::cout << "\t\t | ";
    for (size_t j = 0; j < fMatrix[i].size() - 1; j++)
      std::cout << fMatrix[i][j] << " ";
    std::cout << " | " << fMatrix[i][fMatrix[i].size() - 1] << " |" << std::endl;
  }

  std::cout << "\tSolutions" << std::endl;
  for (size_t i = 0; i < fSolutions.size(); i++)
    std::cout << "\t\t" << i << " " << fSolutions[i] << std::endl;
}
