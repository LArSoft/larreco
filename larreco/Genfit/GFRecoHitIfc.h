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

#ifndef GFRECOHITIFC_H
#define GFRECOHITIFC_H

#include "TMatrixT.h"

#include "larreco/Genfit/GFAbsRecoHit.h"
#include "larreco/Genfit/GFDetPlane.h"

/** @brief RecoHit interface template class. Provides comfortable
 * interface to create RecoHits
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This class defines a comfortable interface to create hit classes in genfit.
 * It is a template class. The template parameter is used to specify a certain
 * basic type of hit:
 *  - GFRecoHitIfc<PlanarHitPolicy> a basic planar hit
 *  - GFRecoHitIfc<SpacepointHitPolicy> a basic space point hit
 *  - GFRecoHitIfc<WirepointHitPolicy> a basic hit on a wire
 *
 * To create a hit for a detector simply inherit from one of the options
 * above and fill in your data. For details look at the respective
 * HitPolicy documentations. You can also directly inherit from
 * GFAbsRecoHit though this is not recommended. If a new hit geometry is needed
 * one should think about implementing a new HitPolicy for this type of hit.
 *
 * @sa PlanarHitPolicy
 * @sa SpacepointHitPolicy
 * @sa WirepointHitPolicy
 *
 * Implementation details: The actual implementations of the methods
 * declared here can be found in the HitPolicy objects.
 */
namespace genf {

  template <class HitPolicy>
  class GFRecoHitIfc : public GFAbsRecoHit {
  protected:
    HitPolicy fPolicy;

  public:
    /** @brief Constructor specifying dimension of hit coordinate vector
   */
    GFRecoHitIfc(int dim) : GFAbsRecoHit(dim) { ; }
    virtual ~GFRecoHitIfc() { ; }

    /** @brief Returns the detector plane object for this hit and a given track
   * representation.
   *
   * The actutal code for this method depends on the hit geometry and is
   * implemented in the HitPolicy
   * @sa PlanarHitPolicy
   * @sa SpacepointHitPolicy
   * @sa WirepointHitPolicy
   */
    virtual const GFDetPlane& getDetPlane(GFAbsTrackRep* rep)
    {
      return fPolicy.detPlane(this, rep);
    }

    /** @brief Get hit coordinates in a specific detector plane
   *
   * Implementation in the HitPolicy
   */
    virtual TMatrixT<Double_t> getHitCoord(const GFDetPlane& plane, const GFDetPlane& planePrev)
    {
      return fPolicy.hitCoord(this, plane, planePrev);
    }
    virtual TMatrixT<Double_t> getHitCoord(const GFDetPlane& plane)
    {
      return fPolicy.hitCoord(this, plane);
    }

    /** @brief Get hit covariances in a specific detector plane
   *
   * Implementation in the HitPolicy
   */
    virtual TMatrixT<Double_t> getHitCov(const GFDetPlane& plane)
    {
      return fPolicy.hitCov(this, plane);
    }
    virtual TMatrixT<Double_t> getHitCov(const GFDetPlane& plane,
                                         const GFDetPlane& planePrev,
                                         const TMatrixT<Double_t>& state,
                                         const Double_t& mass)
    {
      return fPolicy.hitCov(this, plane, planePrev, state, mass);
    }

    const std::string& getPolicyName() { return fPolicy.getName(); }

    //public:
    //  ClassDef(GFRecoHitIfc,1);
  };

}

#endif

/** @} */
