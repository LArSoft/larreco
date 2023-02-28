/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#ifndef GFDAF_H
#define GFDAF_H

#include <map>
#include <vector>

#include "TMatrixT.h"

/** @brief Determinstic Annealing Filter (DAF) implementation
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * The DAF is an iterative Kalman filter with annealing. It is capable of fitting tracks which are
 * contaminated with noise hits. The alogrithm is taken from the references
 * R. Fruehwirth & A. Strandlie, Computer Physics Communications 120 (199) 197-214
 * and CERN thesis: Dissertation by Matthias Winkler.

 */

namespace genf {

  class GFTrack;

  class GFDaf {
  public:
    /** @brief Standard CTOR. Sets default values for annealing scheme and probablity cut.
   */
    GFDaf();

    ~GFDaf();

    /** @brief Performs DAF fit on all track representations in a GFTrack.
   *
   */
    void processTrack(GFTrack*);

    /** @brief Set the blowup factor (see blowUpCovs() )
   */
    void setBlowUpFactor(double f) { fBlowUpFactor = f; }

    /** @brief Set the probabilty cut for the weight calculation for the hits. Currently
   * supported are the values 0.01 0.005, and 0.001. The corresponding chi2 cuts for
   * different hits dimensionalities are hardcoded in the implementation because I did
   * not yet figure out how to calculate them. Please feel very welcome to change the
   * implementtion if you know how to do it.
   */
    void setProbCut(double val);

    /** @brief Configure the annealing scheme.
   * In the current implementation you need to provide at least two temperatures. The maximum would ten
   * tempertatures.
   */
    void setBetas(double b1,
                  double b2,
                  double b3 = -1.,
                  double b4 = -1.,
                  double b5 = -1.,
                  double b6 = -1.,
                  double b7 = -1.,
                  double b8 = -1.,
                  double b9 = -1.,
                  double b10 = -1.);
    // Private Methods -----------------
  private:
    /** @brief Calculate Kalman Gain
   */
    TMatrixT<Double_t> calcGain(const TMatrixT<Double_t>& cov,
                                const TMatrixT<Double_t>& HitCov,
                                const TMatrixT<Double_t>& H,
                                const double& p);

    /** @brief This is needed to blow up the covariance matrix before a fitting pass.
   * The method drops off-diagonal elements and blows up diagonal by blowUpFactor.
   */
    void blowUpCovs(GFTrack* trk);

    /** @brief invert a matrix. First argument is matrix to be inverted, second is return by ref.
   */
    void invertMatrix(const TMatrixT<Double_t>&, TMatrixT<Double_t>&);

    double fBlowUpFactor;
    std::vector<double> fBeta;
    std::map<int, double> chi2Cuts;
  };

} // namespace genf
#endif

/** @} */
