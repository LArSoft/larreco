// \file  TrackMomentumCalculator.cxx
// 
// \author sowjanyag@phys.ksu.edu  

#include "RecoAlg/TrackMomentumCalculator.h"
#include "TMath.h"

namespace trkf{ 
 
  double TrackMomentumCalculator::GetTrackMomentum(double trkrange, int pdg) 
  {
   
   /* Muon range-momentum tables from CSDA (Argon density = 1.4 g/cm^3)
   website: http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf

   CSDA table values:
   Float_t Range_grampercm[30] = {9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2, 2.385E2, 4.934E2,
   				   6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
                                   4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5};
   Float_t KE_MeV[30] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000, 4000, 8000, 10000, 14000, 20000, 30000, 
   			 40000, 80000, 100000, 140000, 200000, 300000, 400000};
			 
   Functions below are obtained by fitting polynomial fits to KE_MeV vs Range (cm) graph. A better fit was obtained
   by splitting the graph into two: Below range<=200cm,a polynomial of power 4 was a good fit; above 200cm, a
   polynomial of power 6 was a good fit
   
   Fit errors for future purposes:
   Below 200cm, Forpoly4 fit: p0 err=1.38533;p1 err=0.209626; p2 err=0.00650077; p3 err=6.42207E-5; p4 err=1.94893E-7;  			 
   Above 200cm, Forpoly6 fit: p0 err=5.24743;p1 err=0.0176229; p2 err=1.6263E-5; p3 err=5.9155E-9; p4 err=9.71709E-13;
   			      p5 err=7.22381E-17;p6 err=1.9709E-21;*/ 
			 
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			      
   //*********For muon, the calculations are valid up to 1.91E4 cm range corresponding to a Muon KE of 40 GeV**********//
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
   website: http://projects-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=1696;filename=Stopping%20Protons%20-%20July.pdf;version=2
   
   CSDA values:
   Double_t KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
                              1500, 2000, 2500, 3000, 4000, 5000};
   
   Double_t Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296, 2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1, 9.184E1, 1.138E2,
                                     1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2, 2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2, 
				     7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3}; 

   Functions below are obtained by fitting power and polynomial fits to KE_MeV vs Range (cm) graph. A better fit was obtained
   by splitting the graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
   polynomial of power 6 was a good fit
   
   Fit errors for future purposes:
   For power function fit: a=0.388873; and b=0.00347075		 
   Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3 err=4.45542E-7; p4 err=4.16428E-10;
   	         p5 err=1.81679E-13;p6 err=2.96958E-17;*/ 

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			      
   //*********For proton, the calculations are valid up to 3.022E3 cm range corresponding to a Muon KE of 5 GeV**********//
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			      
		 					      			       			 			
   double KE, Momentum, Muon_M = 105.7, Proton_M = 938.272, M;
   
   if (abs(pdg) == 13){
     M = Muon_M;		
     if (trkrange>0 && trkrange<=200)
     	  KE = 10.658+ (3.71181*trkrange) + (-0.0302898*trkrange*trkrange) +
	           (0.000228501*trkrange*trkrange*trkrange)+
	           (-5.86201E-7*trkrange*trkrange*trkrange*trkrange);
     else if (trkrange>200 && trkrange<=1.91E4) 	     		
   	  KE = 30.6157+ (2.11089*trkrange) + (0.000255267*trkrange*trkrange)+
	          (-5.19254E-8*trkrange*trkrange*trkrange)+
	          (6.39454E-12*trkrange*trkrange*trkrange*trkrange)+
	          (-3.99392E-16*trkrange*trkrange*trkrange*trkrange*trkrange)+
	          (9.73174E-21*trkrange*trkrange*trkrange*trkrange*trkrange*trkrange);		  
      else
          KE = -999;
      } 	  
   else if (abs(pdg) == 2212){
      M = Proton_M;
      if (trkrange>0 && trkrange<=80)
   	 KE = 29.9317*TMath::Power(trkrange,0.586304);
      else if (trkrange>80 && trkrange<=3.022E3)	
   	 KE = 149.904+ (3.34146*trkrange) + (-0.00318856*trkrange*trkrange)+
	     (4.34587E-6*trkrange*trkrange*trkrange)+
	     (-3.18146E-9*trkrange*trkrange*trkrange*trkrange)+
	     (1.17854E-12*trkrange*trkrange*trkrange*trkrange*trkrange)+
	     (-1.71763E-16*trkrange*trkrange*trkrange*trkrange*trkrange*trkrange);
      else
         KE = -999;
      }
    else
      KE = -999;   
      	  	     
    if (KE<0) 	     			  
       Momentum = -999;  
    else
       Momentum = TMath::Sqrt((KE*KE)+(2*M*KE));
       
    return Momentum;       	  		 	     			  
  }
  
  /*double TrackMomentumCalculator::GetProtonTrackMomentum(double trkrange) 
  {	 		
   double ProtonKE;		
   if (trkrange<=80)
   	ProtonKE = 29.9317*TMath::Power(trkrange,0.586304);
   else	
   	ProtonKE = 149.904+ (3.34146*trkrange) + (-0.00318856*trkrange*trkrange)+
	     (4.34587E-6*trkrange*trkrange*trkrange)+
	     (-3.18146E-9*trkrange*trkrange*trkrange*trkrange)+
	     (1.17854E-12*trkrange*trkrange*trkrange*trkrange*trkrange)+
	     (-1.71763E-16*trkrange*trkrange*trkrange*trkrange*trkrange*trkrange);
	     			  
    return ProtonKE;
  }*/
} // namespace track
