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

#include "Genfit/GFEnergyLossBetheBloch.h"
#include "assert.h"
#include "math.h"
#include <iostream>

genf::GFEnergyLossBetheBloch::~GFEnergyLossBetheBloch(){
}

double genf::GFEnergyLossBetheBloch::energyLoss(const double& step,
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

  static const double me = getParticleMass(11); // electron mass (GeV)  0.000519...

  double charge, mass;
  getParticleParameters(pdg, charge, mass);
  
  const double beta = mom/sqrt(mass*mass+mom*mom);
  double dedx = 0.307075*matZ/matA*matDensity/(beta*beta)*charge*charge;

  //for numerical stability
  double gammaSquare = 1.-beta*beta;
  if(gammaSquare>1.E-10) gammaSquare = 1./gammaSquare;
  else gammaSquare = 1.E10;
  double gamma = sqrt(gammaSquare);

  double massRatio = me/mass;
  double argument = gammaSquare*beta*beta*me*1.E3*2./((1.E-6*meanExcitationEnergy) * sqrt(1+2*sqrt(gammaSquare)*massRatio + massRatio*massRatio));
  if (argument <= exp(beta*beta)) 
    dedx = 0.;
  else{
    dedx *= (log(argument)-beta*beta); // Bethe-Bloch [MeV/cm]
    dedx *= 1.E-3;  // in GeV/cm, hence 1.e-3
    assert(dedx>0);
  }
  
  double DE = step * dedx; //always positive
  double momLoss = sqrt(mom*mom+2.*sqrt(mom*mom+mass*mass)*DE+DE*DE) - mom; //always positive
  
  //in vacuum it can numerically happen that momLoss becomes a small negative number. A cut-off at 0.01 eV for momentum loss seems reasonable
  if(fabs(momLoss)<1.E-11)momLoss=1.E-11;

    
  if (doNoise) {
    assert (noise!=NULL); // assert that optional argument noise exists
    
    // ENERGY LOSS FLUCTUATIONS; calculate sigma^2(E); 
    double sigma2E = 0.;
    double zeta  = 153.4E3 * charge*charge/(beta*beta) * matZ/matA * matDensity * step; // eV
    double Emax  = 2.E9*me*beta*beta*gammaSquare / (1. + 2.*gamma*me/mass + (me/mass)*(me/mass) ); // eV
    double kappa = zeta/Emax;
    
    if (kappa > 0.01) { // Vavilov-Gaussian regime
      sigma2E += zeta*Emax*(1.-beta*beta/2.);  // eV^2
    }
    else { // Urban/Landau approximation
      double alpha = 0.996;
      double sigmaalpha = 15.76;
      // calculate number of collisions Nc
      double I = 16. * pow(matZ, 0.9); // eV
      double f2 = 0.;
      if (matZ > 2.) f2 = 2./matZ;
      double f1 = 1. - f2;
      double e2 = 10.*matZ*matZ; // eV
      double e1 = pow( (I/pow(e2,f2)), 1./f1);  // eV
      
      double mbbgg2 = 2.E9*mass*beta*beta*gammaSquare; // eV
      double Sigma1 = dedx*1.0E9 * f1/e1 * (log(mbbgg2 / e1) - beta*beta) / (log(mbbgg2 / I) - beta*beta) * 0.6; // 1/cm
      double Sigma2 = dedx*1.0E9 * f2/e2 * (log(mbbgg2 / e2) - beta*beta) / (log(mbbgg2 / I) - beta*beta) * 0.6; // 1/cm
      double Sigma3 = dedx*1.0E9 * Emax / ( I*(Emax+I)*log((Emax+I)/I) ) * 0.4; // 1/cm
        
      double Nc = (Sigma1 + Sigma2 + Sigma3)*step;
      
      if (Nc > 50.) { // truncated Landau distribution 
        // calculate sigmaalpha  (see GEANT3 manual W5013)
        double RLAMED = -0.422784 - beta*beta - log(zeta/Emax);
        double RLAMAX =  0.60715 + 1.1934*RLAMED +(0.67794 + 0.052382*RLAMED)*exp(0.94753+0.74442*RLAMED);
        // from lambda max to sigmaalpha=sigma (empirical polynomial)
        if(RLAMAX <= 1010.) {
           sigmaalpha =  1.975560  
                        +9.898841e-02 *RLAMAX  
                        -2.828670e-04 *RLAMAX*RLAMAX
                        +5.345406e-07 *pow(RLAMAX,3.)   
                        -4.942035e-10 *pow(RLAMAX,4.) 
                        +1.729807e-13 *pow(RLAMAX,5.); 
        }
        else { sigmaalpha = 1.871887E+01 + 1.296254E-02 *RLAMAX; }
        // alpha=54.6  corresponds to a 0.9996 maximum cut
        if(sigmaalpha > 54.6) sigmaalpha=54.6;
        sigma2E += sigmaalpha*sigmaalpha * zeta*zeta;  // eV^2
      }
      else { // Urban model
        double Ealpha  = I / (1.-(alpha*Emax/(Emax+I)));   // eV
        double meanE32 = I*(Emax+I)/Emax * (Ealpha-I);     // eV^2
        sigma2E += step * (Sigma1*e1*e1 + Sigma2*e2*e2 + Sigma3*meanE32); // eV^2
      }
    }
        
    sigma2E*=1.E-18; // eV -> GeV
    
    // update noise matrix
    (*noise)[6][6] += (mom*mom+mass*mass)/pow(mom,6.)*sigma2E; 
  }
  
  return momLoss;
}

//ClassImp(GFEnergyLossBetheBloch)
