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
 
#ifndef GFENERGYLOSSBETHEBLOCH_H
#define GFENERGYLOSSBETHEBLOCH_H

#include <iostream>
#include <assert.h>

#include"GFAbsEnergyLoss.h"

  
/** @brief Energy loss for charged particles (Bethe Bloch), energy loss straggeling
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

class GFEnergyLossBetheBloch : public GFAbsEnergyLoss{  
 public:
  //! Returns energy loss, optional calculation of energy loss straggeling
  /**  Uses Bethe Bloch formula to calculate energy loss. 
    *  For the energy loss straggeling, different formulas are used for different regions:
    *  - Vavilov-Gaussian regime
    *  - Urban/Landau approximation
    *  - truncated Landau distribution 
    *  - Urban model 
    * 
  */
  double energyLoss(const double& step,
                    const double& mom,
                    const int&    pdg,              
                    const double& matDensity,
                    const double& matZ,
                    const double& matA,
                    const double& radiationLength,
                    const double& meanExcitationEnergy,
                    const bool&   doNoise = false,
                          TMatrixT<Double_t>* noise = NULL,
                    const TMatrixT<Double_t>* jacobian = NULL,
                    const TVector3* directionBefore = NULL,
                    const TVector3* directionAfter = NULL);
  virtual ~GFEnergyLossBetheBloch();
 
  // public:
  //ClassDef(GFEnergyLossBetheBloch,1); 
 
};

} // namespace genf

#endif

/* @} **/

