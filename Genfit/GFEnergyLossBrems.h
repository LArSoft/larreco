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
 
#ifndef GFENERGYLOSSBREMS_H
#define GFENERGYLOSSBREMS_H

#include <iostream>
#include <assert.h>

#include"GFAbsEnergyLoss.h"

  
/** @brief Energy loss for electrons/positrons due to bremsstrahlung, energy loss straggeling
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

class GFEnergyLossBrems : public GFAbsEnergyLoss{  
 public:
  //! Returns energy loss, optional calculation of energy loss straggeling
  /** Can be called with any pdg, but only calculates energy loss and straggeling for electrons and positrons (otherwise returns 0).
    * Uses a gaussian approximation (Bethe-Heitler formula with Migdal corrections). 
    * For positrons the energy loss is weighed with a correction factor.
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
  virtual ~GFEnergyLossBrems();

  // public:
  //ClassDef(GFEnergyLossBrems,1);
  
};
}

#endif

/* @} **/

