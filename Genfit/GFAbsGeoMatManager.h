/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

/** @addtogroup RKTrackRep
 * @{
 */
 
#ifndef GFABSGEOMATMANAGER_H
#define GFABSGEOMATMANAGER_H

#include"TObject.h"
#include<iostream>
#include <assert.h>


/** @brief Base class for material and geometry interface
 * 
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, author)
 * 
 *  
 *  
 *  
 *  
 */

namespace genf {

class GFAbsGeoMatManager : public TObject{
 public:
  virtual ~GFAbsGeoMatManager(){}
  //! Gets material parameters (density, Z, A, radiation length, mean excitation energy)
  /**  
  */
  virtual void getMaterialParameters(double& matDensity,
                                     double& matZ,
                                     double& matA,
                                     double& radiationLength,
                                     double& mEE) = 0;
 
  //! Initializes the track
  /**  
  */
  virtual void initTrack(const double& posx,
                         const double& posy,
                         const double& posz,
                         const double& dirx,
                         const double& diry,
                         const double& dirz) = 0;

  //! Makes a step, limited to next material boundary
  /**  Tries to make a step with length maxDist along the track. If there is a material boundary, 
    *  the step is made to that boundary and the distance to that boundary is returned. 
    *  Otherwise the step is made with maxDist.
    * 
  */    
  virtual double stepOrNextBoundary(const double& maxDist) = 0;
     
  // public:
  //ClassDef(GFAbsGeoMatManager,1)
 
};

} // namespace genf
#endif

/* @} **/


