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
/** @addtogroup genfit */
/* @{ */

#ifndef GFABSRECOHIT_H
#define GFABSRECOHIT_H

#include <iostream>
#include <stdexcept> // std::logic_error

#include "TMatrixT.h"
#include "TObject.h"

#include "larreco/Genfit/GFDetPlane.h"
#include "larreco/Genfit/GFException.h" // PrintROOTobject()
#include <cmath>

namespace genf {
  class GFAbsTrackRep;
}

/** @brief Base Class for representing a Hit in GENFIT
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * A hit is defined as a single measurement of a detector. Each detector can
 * define it's own hit representation (geometry) by inherting from GFAbsRecoHit.
 * We call such a child object a "RecoHit" to make clear that
 * it inherits from GFAbsRecoHit.
 * All detector specific information is handled inside the RecoHit objects.
 * The idea behind this is that the fitting algorithms can work on different
 * detectors simultanously.
 *
 * GFAbsRecoHit defines the basic interface that is used by all genfit algorithms
 * to access hit-measurement information. It provides:
 *  - the matrices to store hit coordinates and covariances
 *  - defines methods to access these data
 *  - the interface to access a hit's detector plane object
 *
 * All hits have to inherit from this base class.
 * Inheritance can be direct or through template class
 * RecoHitIfc<GeometryPolicy>.
 * These interface classes (defined with a specific policy)
 * provide additional functionality for specific hit geometries,
 * such as space points, wires, etc. For details look
 * at the RecoHitIfc documentation.
 *
 * A simple example is given in VirtSpacePointRecoHit
 *
 * Background information: The main feature here is
 * that coordinates and covariances are available as general
 * TMatrixT<Double_t> objects. By using these general algebraic
 * objects it is possible to abstract from the concrete measurement and
 * thus define a unified framework to handle the data. This is a key ingredient
 * for the separation of data and algorithms which in turn enables us
 * to elegantly combine information from different detectors.
 */
namespace genf {

  class GFAbsRecoHit : public TObject {
  protected:
    /// Vector of raw coordinates of hit
    TMatrixT<Double_t> fHitCoord;

    /// Covariance of raw hit coordinates
    TMatrixT<Double_t> fHitCov;

  private:
    int fNparHit;

  public:
    virtual ~GFAbsRecoHit();

    /** @brief Constructor specifying dimension of coordinate vector
   *
   * One hit is generally represented by a vector of coordinates.
   * @param NparHit specifies the dimension of this vector.
   * (e.g. 3 for a spacepoint, 2 for a pixel, 1 for a strip)
   */
    GFAbsRecoHit(int NparHit);

    /** @brief Default constructor needed for compatibility with ROOT
   */
    GFAbsRecoHit();

    /** @brief Get transformation matrix. Transformation between hit
   * coordinates and track representation coordinates.
   *
   * This is a virtual abstract method which has to be implemented in the child
   * classes.
   *
   * In general there is a linear transformation between the coordinate system
   * of the hit (which is defined by the detector plane) and the coordinates of
   * the track representation in that plane. In the most simple case the track
   * representation has 5 parameters (space + momentum)
   * while a hit usually has less (one to three space coordinates).
   *
   * The transformation matrix is then simply projecting out the
   * space-components of the track representation.
   *
   * Its dimensions are NxM. Where N is the number of dimensions of the
   * hit in the detector plane (usually 2 or 1) and M is the dimension of the
   * track representation.
   *
   * In this method a hit has to define with which track representations it can
   * work together. It should be the only point where this explicit
   * coordination is necessary.
   *
   * For example code see implementing classes below:
   */
    virtual TMatrixT<Double_t> getHMatrix(const GFAbsTrackRep* stateVector) = 0;
    virtual TMatrixT<Double_t> getHMatrix(const GFAbsTrackRep* stateVector,
                                          const Double_t&,
                                          const Double_t&) = 0;

    /** @brief Calculate residual with respect to a track representation.
   *
   * Returns the N-dimensional residual of this vector to a given
   * track representation.
   *
   * This method is not doing any extrapolation. But it
   * creates the necessary detector plane object. See GFAbsRecoHit::getGFDetPlane
   *
   * @param stateVector pointer to track representation - used to synchronize
   * with the track repesentation
   * @param state parameter vector of the track representation
   *
   * @sa setHMatrix
   * @sa getGFDetPlane
   */
    virtual TMatrixT<Double_t> residualVector(const GFAbsTrackRep* stateVector,
                                              const TMatrixT<Double_t>& state,
                                              const GFDetPlane& d)
    {
      std::cout << "GFAbsRecoHit::residualVector(3args) Not correctly Using theta -- multiple "
                   "scattering -- information !!! Fix this if you really want to use getChi2Hit"
                << std::endl;
      TMatrixT<Double_t> H = getHMatrix(stateVector);
      return (getHitCoord(d) - (H * state));
    }

    virtual TMatrixT<Double_t> residualVector(const GFAbsTrackRep* stateVector,
                                              const TMatrixT<Double_t>& state,
                                              const GFDetPlane& d,
                                              const GFDetPlane& dPrev,
                                              const double& mass)
    {
      Double_t dist = (d.getO() - dPrev.getO()).Mag();
      //    Double_t mass = 0.104; // close enough for muons, pions, I think.
      Double_t mom = fabs(1.0 / state[0][0]);
      Double_t beta = mom / sqrt(mass * mass + mom * mom);
      if (std::isnan(dist) || dist < 0.2) dist = 0.2; // don't allow 0s here.
      if (std::isnan(beta) || beta < 0.04) beta = 0.04;

      TMatrixT<Double_t> H = getHMatrix(stateVector, beta, dist);
      return (getHitCoord(d, dPrev) - (H * state));
    }

    /** @brief Get raw hit covariances.
   *
   */
    TMatrixT<Double_t> getRawHitCov() const { return fHitCov; }

    /** @brief Get raw hit coordinates.
   *
   */
    TMatrixT<Double_t> getRawHitCoord() const { return fHitCoord; }

    /** @brief Get hit covariances in a specific detector plane
   *
   * Virtual abstract method has to be implemented by inherting classes.
   * Implementation involves transformation from raw coordinates in detector
   * coordinate system to detector plane coordinate system.
   * @sa getGFDetPlane
   */
    virtual TMatrixT<Double_t> getHitCov(const GFDetPlane&) = 0;
    virtual TMatrixT<Double_t> getHitCov(const GFDetPlane&,
                                         const GFDetPlane&,
                                         const TMatrixT<Double_t>&,
                                         const Double_t&) = 0;

    /** @brief Get hit coordinates in a specific detector plane
   *
   * Virtual abstract method has to be implemented by inherting classes.
   * Implementation involves transformation from raw coordinates in detector
   * coordinate system to detector plane coordinate system.
   * @sa getDetPlane
   */
    virtual TMatrixT<Double_t> getHitCoord(const GFDetPlane&, const GFDetPlane&) = 0;
    virtual TMatrixT<Double_t> getHitCoord(const GFDetPlane&) = 0;

    /** @brief Get detector plane for a given track representation.
   *
   * Virtual abstract method has to be implemented by inherting classes.
   *
   * In general the detector plane can depend both on the detector/hit geometry
   * as well as the track geometry. E.g. for a space point one usually chooses
   * a plane that is perpendicular to the current track, since in that case no
   * other plane is predefined.
   *
   * There are several implementations for this method in the HitPolicies.
   * In the most simple case (a planar detector) the method just returns a
   * fixed (detector module specific) plane. This behaviour for example is
   * implemented in PlanarHitPolicy.
   */
    virtual const GFDetPlane& getDetPlane(GFAbsTrackRep*) = 0;

    /** @brief Get clone of this object.
   *
   * Virtual abstract method. Has to be implemented by inherting classes.
   * Creates a deep copy of this object.
   * Ownership is trandsferred to the caller!
   */
    virtual GFAbsRecoHit* clone() = 0;

    /** @brief Print raw hit coordinates.
   */
    virtual void Print(std::ostream& out = std::cout) const { PrintROOTobject(out, fHitCoord); }

    virtual const std::string& getPolicyName();

    int getNparHit() { return fNparHit; }

    // public:
    //ClassDef(GFAbsRecoHit,3)

  private:
    virtual void Print(Option_t*) const
    {
      throw std::logic_error(std::string(__func__) + "::Print(Option_t*) not available");
    }
  };
} // namespace

#endif //FITTER_ABSHIT_H

/** @} */
