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
 
#ifndef GFGEOMATMANAGER_H
#define GFGEOMATMANAGER_H

#include"GFAbsGeoMatManager.h"

#include"TObject.h"
#include<iostream>
#include <assert.h>

  
/** @brief Material and geometry interface via TGeoMaterial and gGeoManager
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

class GFGeoMatManager : public GFAbsGeoMatManager{
 public:
   virtual ~GFGeoMatManager(){}
  void getMaterialParameters(double& matDensity,
                             double& matZ,
                             double& matA,
                             double& radiationLength,
                             double& mEE);
 
  void initTrack(const double& posx,
                 const double& posy,
                 const double& posz,
                 const double& dirx,
                 const double& diry,
                 const double& dirz);
 
  double stepOrNextBoundary(const double& maxDist);
     
  // public:
  //ClassDef(GFGeoMatManager,1)
 
};

} // namespace genf
#endif

/* @} **/


