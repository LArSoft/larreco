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

#include "Genfit/GFEnergyLossBrems.h"
#include "assert.h"
#include "math.h"
#include <iostream>

//#define BETHE  // don't use MIGDAL correction for electron bremsstrahlung (only Bethe Heitler)

genf::GFEnergyLossBrems::~GFEnergyLossBrems(){
}

double genf::GFEnergyLossBrems::energyLoss(const double& step,
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

  if (fabs(pdg==11)) { // only for electrons and positrons
    #if !defined(BETHE)
      static const double C[101]={ 0.0,-0.960613E-01, 0.631029E-01,-0.142819E-01, 0.150437E-02,-0.733286E-04, 0.131404E-05, 0.859343E-01,-0.529023E-01, 0.131899E-01,-0.159201E-02, 0.926958E-04,-0.208439E-05,-0.684096E+01, 0.370364E+01,-0.786752E+00, 0.822670E-01,-0.424710E-02, 0.867980E-04,-0.200856E+01, 0.129573E+01,-0.306533E+00, 0.343682E-01,-0.185931E-02, 0.392432E-04, 0.127538E+01,-0.515705E+00, 0.820644E-01,-0.641997E-02, 0.245913E-03,-0.365789E-05, 0.115792E+00,-0.463143E-01, 0.725442E-02,-0.556266E-03, 0.208049E-04,-0.300895E-06,-0.271082E-01, 0.173949E-01,-0.452531E-02, 0.569405E-03,-0.344856E-04, 0.803964E-06, 0.419855E-02,-0.277188E-02, 0.737658E-03,-0.939463E-04, 0.569748E-05,-0.131737E-06,-0.318752E-03, 0.215144E-03,-0.579787E-04, 0.737972E-05,-0.441485E-06, 0.994726E-08, 0.938233E-05,-0.651642E-05, 0.177303E-05,-0.224680E-06, 0.132080E-07,-0.288593E-09,-0.245667E-03, 0.833406E-04,-0.129217E-04, 0.915099E-06,-0.247179E-07, 0.147696E-03,-0.498793E-04, 0.402375E-05, 0.989281E-07,-0.133378E-07,-0.737702E-02, 0.333057E-02,-0.553141E-03, 0.402464E-04,-0.107977E-05,-0.641533E-02, 0.290113E-02,-0.477641E-03, 0.342008E-04,-0.900582E-06, 0.574303E-05, 0.908521E-04,-0.256900E-04, 0.239921E-05,-0.741271E-07,-0.341260E-04, 0.971711E-05,-0.172031E-06,-0.119455E-06, 0.704166E-08, 0.341740E-05,-0.775867E-06,-0.653231E-07, 0.225605E-07,-0.114860E-08,-0.119391E-06, 0.194885E-07, 0.588959E-08,-0.127589E-08, 0.608247E-10};
      static const double xi=2.51, beta=0.99, vl=0.00004;
    #endif
    #if defined(BETHE) // no MIGDAL corrections
      static const double C[101]={ 0.0, 0.834459E-02, 0.443979E-02,-0.101420E-02, 0.963240E-04,-0.409769E-05, 0.642589E-07, 0.464473E-02,-0.290378E-02, 0.547457E-03,-0.426949E-04, 0.137760E-05,-0.131050E-07,-0.547866E-02, 0.156218E-02,-0.167352E-03, 0.101026E-04,-0.427518E-06, 0.949555E-08,-0.406862E-02, 0.208317E-02,-0.374766E-03, 0.317610E-04,-0.130533E-05, 0.211051E-07, 0.158941E-02,-0.385362E-03, 0.315564E-04,-0.734968E-06,-0.230387E-07, 0.971174E-09, 0.467219E-03,-0.154047E-03, 0.202400E-04,-0.132438E-05, 0.431474E-07,-0.559750E-09,-0.220958E-02, 0.100698E-02,-0.596464E-04,-0.124653E-04, 0.142999E-05,-0.394378E-07, 0.477447E-03,-0.184952E-03,-0.152614E-04, 0.848418E-05,-0.736136E-06, 0.190192E-07,-0.552930E-04, 0.209858E-04, 0.290001E-05,-0.133254E-05, 0.116971E-06,-0.309716E-08, 0.212117E-05,-0.103884E-05,-0.110912E-06, 0.655143E-07,-0.613013E-08, 0.169207E-09, 0.301125E-04,-0.461920E-04, 0.871485E-05,-0.622331E-06, 0.151800E-07,-0.478023E-04, 0.247530E-04,-0.381763E-05, 0.232819E-06,-0.494487E-08,-0.336230E-04, 0.223822E-04,-0.384583E-05, 0.252867E-06,-0.572599E-08, 0.105335E-04,-0.567074E-06,-0.216564E-06, 0.237268E-07,-0.658131E-09, 0.282025E-05,-0.671965E-06, 0.565858E-07,-0.193843E-08, 0.211839E-10, 0.157544E-04,-0.304104E-05,-0.624410E-06, 0.120124E-06,-0.457445E-08,-0.188222E-05,-0.407118E-06, 0.375106E-06,-0.466881E-07, 0.158312E-08, 0.945037E-07, 0.564718E-07,-0.319231E-07, 0.371926E-08,-0.123111E-09};
      static const double xi=2.10, beta=1.00, vl=0.001;
    #endif                                        
                                           
    double BCUT=10000.; // energy up to which soft bremsstrahlung energy loss is calculated
      
    static const double me = getParticleMass(11); // electron mass (GeV)  0.000519...

    double charge, mass;
    getParticleParameters(pdg, charge, mass);
    
    double THIGH=100., CHIGH=50.;

    double dedxBrems=0.;

    if(BCUT>0.){

      double T, kc;
      
      if(BCUT>=mom) BCUT=mom; // confine BCUT to mom
      
      // T=mom,  confined to THIGH
      // kc=BCUT, confined to CHIGH ??
      if(mom>=THIGH) {
        T=THIGH;
        if(BCUT>=THIGH) kc=CHIGH;
        else kc=BCUT;
      }
      else {
        T=mom;
        kc=BCUT;
      }

      double E=T+me; // total electron energy
      if(BCUT>T) kc=T;
      
      double X=log(T/me);
      double Y=log(kc/(E*vl));

      double XX;
      int    K;
      double S=0., YY=1.;

      for (unsigned int I=1; I<=2; ++I) {
        XX=1.;
        for (unsigned int J=1; J<=6; ++J) {
          K=6*I+J-6;
          S=S+C[K]*XX*YY;
          XX=XX*X;
        }
        YY=YY*Y;
      }

      for (unsigned int I=3; I<=6; ++I) {
         XX=1.;
         for (unsigned int J=1; J<=6; ++J) {
            K=6*I+J-6;
            if(Y<=0.) S=S+C[K]*XX*YY;
            else      S=S+C[K+24]*XX*YY;
            XX=XX*X;
         }
         YY=YY*Y;
      }

      double SS=0.;
      YY=1.;

      for (unsigned int I=1; I<=2; ++I) {
        XX=1.;
        for (unsigned int J=1; J<=5; ++J) {
          K=5*I+J+55;
          SS=SS+C[K]*XX*YY;
          XX=XX*X;
        }
        YY=YY*Y;
      }
       
      for (unsigned int I=3; I<=5; ++I) {
        XX=1.;
        for (unsigned int J=1; J<=5; ++J) {
          K=5*I+J+55;
          if(Y<=0.) SS=SS+C[K]*XX*YY;
          else      SS=SS+C[K+15]*XX*YY;
          XX=XX*X;
        }
        YY=YY*Y;
      }

      S=S+matZ*SS; 
      
      if(S>0.){
        double CORR=1.;
        #if !defined(CERNLIB_BETHE)
          CORR=1./(1.+0.805485E-10*matDensity*matZ*E*E/(matA*kc*kc)); // MIGDAL correction factor  
        #endif

        double FAC=matZ*(matZ+xi)*E*E * pow((kc*CORR/T),beta) / (E+me);
        if(FAC<=0.) return 0.;
        dedxBrems=FAC*S;
        
        double RAT;

        if(mom>THIGH) {
          if(BCUT<THIGH) {
            RAT=BCUT/mom;
            S=(1.-0.5*RAT+2.*RAT*RAT/9.);
            RAT=BCUT/T;
            S=S/(1.-0.5*RAT+2.*RAT*RAT/9.);
          }
          else {
            RAT=BCUT/mom;
            S=BCUT*(1.-0.5*RAT+2.*RAT*RAT/9.);
            RAT=kc/T;
            S=S/(kc*(1.-0.5*RAT+2.*RAT*RAT/9.));
          }
          dedxBrems=dedxBrems*S; // GeV barn
        }

        dedxBrems = 0.60221367*matDensity*dedxBrems/matA; // energy loss dE/dx [GeV/cm]
      }
    }                                         

    assert(dedxBrems>=0.);                                         

    double factor=1.; // positron correction factor                                       

    if (pdg==-11){
        static const double AA=7522100., A1=0.415, A3=0.0021, A5=0.00054;

        double ETA=0.;
        if(matZ>0.) {
          double X=log(AA*mom/matZ*matZ);
          if(X>-8.) {
            if(X>=+9.) ETA=1.;
            else {
               double W=A1*X+A3*pow(X,3.)+A5*pow(X,5.);
               ETA=0.5+atan(W)/M_PI;
            }
          }
        }

        double E0;
        
        if(ETA<0.0001) factor=1.E-10;
        else if (ETA>0.9999) factor=1.;
        else {
           E0=BCUT/mom;
           if(E0>1.) E0=1.;
           if(E0<1.E-8) factor=1.;
           else factor = ETA * ( 1.-pow(1.-E0, 1./ETA) ) / E0;
        }
    }

    double DE = step * factor*dedxBrems; //always positive
    double momLoss = sqrt(mom*mom+2.*sqrt(mom*mom+mass*mass)*DE+DE*DE) - mom; //always positive                                           

    
    if (doNoise) {
      assert (noise!=NULL); // assert that optional argument noise exists
      
      double LX  = 1.442695*step/radiationLength;
      double S2B = mom*mom * ( 1./pow(3.,LX) - 1./pow(4.,LX) );
      double DEDXB  = pow(fabs(S2B),0.5);
      DEDXB = 1.2E9*DEDXB; //eV
      double sigma2E = DEDXB*DEDXB; //eV^2
      sigma2E*=1.E-18; // eV -> GeV
      
      (*noise)[6][6] += (mom*mom+mass*mass)/pow(mom,6.)*sigma2E; 
    }

    return momLoss;
  }
  else return 0.; // if not an electron/positron
}

//ClassImp(GFEnergyLossBrems)

