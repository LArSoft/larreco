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

#include "Genfit/GFEnergyLossCoulomb.h"
#include "assert.h"
#include "math.h"
#include <iostream>

genf::GFEnergyLossCoulomb::~GFEnergyLossCoulomb(){
}

double genf::GFEnergyLossCoulomb::energyLoss(const double& step,
                                       const double& mom,
                                       const int&    pdg,              
                                       const double& matDensity,
                                       const double& matZ,
                                       const double& matA,
                                       const double& radiationLength,
                                       const double& meanExcitationEnergy,
                                       const bool&   doNoise,
                                             TMatrixT<Double_t>* noise,
                                       const TMatrixT<Double_t>* jacobian,
                                       const TVector3* directionBefore,
                                       const TVector3* directionAfter){

  double charge, mass;
  getParticleParameters(pdg, charge, mass);
  
  const double beta = mom/sqrt(mass*mass+mom*mom);
 
  if (doNoise) {
    assert (noise!=NULL); // assert that optional argument noise exists
    assert (jacobian!=NULL); // assert that optional argument jacobian exists
    assert (directionBefore!=NULL); // assert that optional argument directionBefore exists
    assert (directionAfter!=NULL); // assert that optional argument directionAfter exists
    
    
    // MULTIPLE SCATTERING; calculate sigma^2 
    // PANDA report PV/01-07 eq(43); linear in step length
    double sigma2 = 225.E-6/(beta*beta*mom*mom) * step/radiationLength * matZ/(matZ+1) * log(159.*pow(matZ,-1./3.))/log(287.*pow(matZ,-0.5)); // sigma^2 = 225E-6/mom^2 * XX0/beta^2 * Z/(Z+1) * ln(159*Z^(-1/3))/ln(287*Z^(-1/2)
        
    
    // noiseBefore
      TMatrixT<Double_t> noiseBefore(7,7);
      
      // calculate euler angles theta, psi (so that directionBefore' points in z' direction)
      double psi = 0;
      if (fabs((*directionBefore)[1]) < 1E-14) {  // numerical case: arctan(+-inf)=+-PI/2
        if ((*directionBefore)[0] >= 0.) psi = M_PI/2.; 
        else psi = 3.*M_PI/2.;
      }
      else {
        if ((*directionBefore)[1] > 0.) psi = M_PI - atan((*directionBefore)[0]/(*directionBefore)[1]);
        else psi = -atan((*directionBefore)[0]/(*directionBefore)[1]);
      }
      // cache sin and cos
      double sintheta = sqrt(1-(*directionBefore)[2]*(*directionBefore)[2]);  // theta = arccos(directionBefore[2])
      double costheta = (*directionBefore)[2];  
      double sinpsi = sin(psi);
      double cospsi = cos(psi);
   
      // calculate NoiseBefore Matrix R M R^T  
      double noiseBefore34 =  sigma2 * cospsi * sinpsi * sintheta*sintheta; // noiseBefore_ij = noiseBefore_ji
      double noiseBefore35 = -sigma2 * costheta * sinpsi * sintheta;
      double noiseBefore45 =  sigma2 * costheta * cospsi * sintheta;

      noiseBefore[3][3] = sigma2 * (cospsi*cospsi + costheta*costheta - costheta*costheta * cospsi*cospsi);
      noiseBefore[4][3] = noiseBefore34;
      noiseBefore[5][3] = noiseBefore35;

      noiseBefore[3][4] = noiseBefore34;
      noiseBefore[4][4] = sigma2 * (sinpsi*sinpsi + costheta*costheta * cospsi*cospsi);
      noiseBefore[5][4] = noiseBefore45;

      noiseBefore[3][5] = noiseBefore35;
      noiseBefore[4][5] = noiseBefore45;
      noiseBefore[5][5] = sigma2 * sintheta*sintheta;
      
      TMatrixT<Double_t> jacobianT(7,7);
      jacobianT = (*jacobian);
      jacobianT.T();
    
      noiseBefore = jacobianT*noiseBefore*(*jacobian); //propagate
    
    // noiseAfter
      TMatrixT<Double_t> noiseAfter(7,7);
      
      // calculate euler angles theta, psi (so that A' points in z' direction)
      psi = 0;
      if (fabs((*directionAfter)[1]) < 1E-14) {  // numerical case: arctan(+-inf)=+-PI/2
        if ((*directionAfter)[0] >= 0.) psi = M_PI/2.; 
        else psi = 3.*M_PI/2.;
      }
      else {
        if ((*directionAfter)[1] > 0.) psi = M_PI - atan((*directionAfter)[0]/(*directionAfter)[1]);
        else psi = -atan((*directionAfter)[0]/(*directionAfter)[1]);
      }
      // cache sin and cos
      sintheta = sqrt(1-(*directionAfter)[2]*(*directionAfter)[2]);  // theta = arccos(directionAfter[2])
      costheta = (*directionAfter)[2];  
      sinpsi = sin(psi);
      cospsi = cos(psi);
   
      // calculate NoiseAfter Matrix R M R^T  
      double noiseAfter34 =  sigma2 * cospsi * sinpsi * sintheta*sintheta; // noiseAfter_ij = noiseAfter_ji
      double noiseAfter35 = -sigma2 * costheta * sinpsi * sintheta;
      double noiseAfter45 =  sigma2 * costheta * cospsi * sintheta;

      noiseAfter[3][3] = sigma2 * (cospsi*cospsi + costheta*costheta - costheta*costheta * cospsi*cospsi);
      noiseAfter[4][3] = noiseAfter34;
      noiseAfter[5][3] = noiseAfter35;

      noiseAfter[3][4] = noiseAfter34;
      noiseAfter[4][4] = sigma2 * (sinpsi*sinpsi + costheta*costheta * cospsi*cospsi);
      noiseAfter[5][4] = noiseAfter45;

      noiseAfter[3][5] = noiseAfter35;
      noiseAfter[4][5] = noiseAfter45;
      noiseAfter[5][5] = sigma2 * sintheta*sintheta;
      
    //calculate mean of noiseBefore and noiseAfter and update noise
      (*noise) += 0.5*noiseBefore + 0.5*noiseAfter;
  }
  
  return 0.;
}

//ClassImp(GFEnergyLossCoulomb)
