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
 
#ifndef GFABSENERGYLOSS_H
#define GFABSENERGYLOSS_H

#include <iostream>
#include <assert.h>

#include"TObject.h"
#include<vector>
#include"TVector3.h"
#include"TMatrixT.h"  

 
/** @brief Base class for energy loss and noise matrix calculation
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

class GFAbsEnergyLoss : public TObject{  
 public:
  //! Calculates energy loss in a given step, optional calculation of noise matrix
  /** 
  */
  virtual double energyLoss(const double& step,
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
                            const TVector3* directionAfter = NULL) = 0;
  virtual ~GFAbsEnergyLoss();
 
 protected:
  //! Gets particle charge and mass (in GeV)
  void getParticleParameters(const int&    pdg,
                             double& charge,
                             double& mass);
  //! Returns particle mass (in GeV)
  double getParticleMass(const int& pdg);

  // public:
  //ClassDef(GFAbsEnergyLoss,1);
  
};

} // namespace genf
#endif

/* @} **/

