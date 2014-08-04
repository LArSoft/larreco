// \file  TrackMomentumCalculator.cxx
// 
// \author sowjanyag@phys.ksu.edu  

#include "RecoAlg/TrackMomentumCalculator.h"

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
    
  
  // MultiScatter business ...
  
  // Author: Leonidas N. Kalousis (August 2014)
  
  double TrackMomentumCalculator::MyMCSChi2( const double *x )
  {
    Double_t result = 0.0;
    
    Double_t p = x[0]; 
    
    Double_t theta0 = x[1]; 
    
    // xmeas[10000]; Double_t ymeas[10000]; Double_t eymeas[10000]; Int_t n_gr;
    
    for ( Int_t i=0; i<n_gr; i++ )
      {
	Double_t xx = xmeas[i]; 
	
	Double_t yy = ymeas[i]; 
	
	Double_t ey = eymeas[i]; 
	
	Double_t rad_length = 14.0;
	
	Double_t l0 = xx/rad_length;
	
	Double_t res1 = 0.0;
	
	if ( xx>0 && p>0 ) res1 = ( 13.6/p )*sqrt( l0 )*( 1.0+0.038*TMath::Log( l0 ) );
	
	res1 = sqrt( res1*res1+theta0*theta0 );
	
	Double_t diff = yy-res1;
	
	if ( ey==0 ) { cout << " Zero denominator in my_mcs_chi2 ! " << endl; getchar(); }
	
	Double_t incre = pow( diff/ey, 2.0 );
	
	result += incre;
	
      }
        
    return result;
    
  }
  
  Int_t TrackMomentumCalculator::GetSegTracks( std::vector<Float_t> *xxx, std::vector<Float_t> *yyy, std::vector<Float_t> *zzz )
  {
    Int_t a1 = xxx->size(); Int_t a2 = yyy->size(); Int_t a3 = zzz->size();
    
    if ( ( a1!=a2 ) || ( a1!=a3 ) || ( a2!=a3 ) ) { cout << " ( Digitize reco tacks ) Error ! " << endl; return -1; }
    
    Int_t stopper = stop / seg_size; 
    	  
    Int_t a4 = a1-1;
        
    segx->clear(); segy->clear(); segz->clear(); segnx->clear(); segny->clear(); segnz->clear(); 
        
    Double_t x0 = xxx->at( 0 ); Double_t y0 = yyy->at( 0 ); Double_t z0 = zzz->at( 0 );
    
    segx->push_back( x0 ); segy->push_back( y0 ); segz->push_back( z0 );
    
    n_seg = 0; x_seg[ n_seg ] = x0; y_seg[ n_seg ] = y0; z_seg[ n_seg ] = z0;
    
    n_seg++;
    
    Int_t ntot = 0; 
    
    Double_t Sz = 0; Double_t Szz = 0;
    
    Double_t Szx = 0; Double_t Sx = 0;
    
    Double_t Szy = 0; Double_t Sy = 0;
    
    Sz += z0; Szz += z0*z0;
    
    Szx += z0*x0; Sx += x0;
    
    Szy += z0*y0; Sy += y0;
    
    ntot++;
    
    for ( Int_t i=1; i<a4; i++ )
      {
	Double_t x1 = xxx->at( i ); Double_t y1 = yyy->at( i );	Double_t z1 = zzz->at( i );
	
	Double_t dr1 = sqrt( pow( x1-x0, 2 )+pow( y1-y0, 2)+pow( z1-z0, 2 ) ); 
	
	Double_t x2 = xxx->at( i+1 ); Double_t y2 = yyy->at( i+1 ); Double_t z2 = zzz->at( i+1 );
	
	Double_t dr2 = sqrt( pow( x2-x0, 2 )+pow( y2-y0, 2)+pow( z2-z0, 2 ) ); 
	
	if ( dr1<seg_size )
	  {
	    Sz += z1; Szz += z1*z1;
    
	    Szx += z1*x1; Sx += x1;
    
	    Szy += z1*y1; Sy += y1;
	    
	    ntot++;
    
	  }
	
	if ( dr1<seg_size && dr2>seg_size )
	  {
	    // ..
	    
	    // cout << " 1 " << endl;
	    
	    Double_t t = 0.0;
	    
	    Double_t dx = x2-x1; Double_t dy = y2-y1; Double_t dz = z2-z1;
	    
	    Double_t dr = sqrt( dx*dx+dy*dy+dz*dz );
	    
	    if ( dr==0 ) { cout << " ( Zero ) Error ! " << endl; getchar(); }
	    
	    Double_t beta = 2.0*( (x1-x0)*dx+(y1-y0)*dy+(z1-z0)*dz )/( dr*dr );
	    	    
	    Double_t gamma = ( dr1*dr1 - seg_size*seg_size )/( dr*dr );
	    
	    Double_t delta = beta*beta - 4.0*gamma;
		
	    if ( delta<0.0 ) { cout << " ( Discriminant ) Error ! " << endl; getchar(); }
		
	    Double_t lysi1 = ( -beta+sqrt( delta ) )/2.0; 
	    	    
	    t=lysi1;
	    	    
	    Double_t xp = x1+t*dx;
	    
	    Double_t yp = y1+t*dy;
	    
	    Double_t zp = z1+t*dz;
	    	    	    	    
	    segx->push_back( xp ); segy->push_back( yp ); segz->push_back( zp );
	    
	    x_seg[ n_seg ] = xp; y_seg[ n_seg ] = yp; z_seg[ n_seg ] = zp; n_seg++; 
	    	    
	    x0 = xp; y0 = yp; z0 = zp;
	    	    	     
	    Sz += z0; Szz += z0*z0;
	    
	    Szx += z0*x0; Sx += x0;
	    
	    Szy += z0*y0; Sy += y0;
	    
	    ntot++;
	    
	    Double_t ax = ( Szx*ntot-Sz*Sx )/( Szz*ntot-Sz*Sz );
	      
	    Double_t ay = ( Szy*ntot-Sz*Sy )/( Szz*ntot-Sz*Sz );
	    
	    segnx->push_back( ax ); segny->push_back( ay ); segnz->push_back( 1.0 );
	    
	    // Double_t angx = find_angle( 1.0, ax ); Double_t angy = find_angle( 1.0, ay );
	    
	    // cout << angx*0.001*180.0/3.14 << endl;
	    
	    ntot = 0; 
    
	    Sz = 0; Szz = 0; 
	    
	    Szx = 0; Sx = 0; 
	    
	    Szy = 0; Sy = 0;
	    
	    Sz += z0; Szz += z0*z0;
	     
	    Szx += z0*x0; Sx += x0;
	    
	    Szy += z0*y0; Sy += y0;
	    
	    ntot++;
	    	    
	  }
	
	else if ( dr1>=seg_size )
	  {
	    // ..
	    
	    // cout << " 2 " << endl;
	    
	    Double_t t = 0.0;
	    	    
	    Double_t dx = x1-x0; Double_t dy = y1-y0; Double_t dz = z1-z0;
	    	    
	    Double_t dr = sqrt( dx*dx+dy*dy+dz*dz );
	    
	    if ( dr==0 ) { cout << " ( Zero ) Error ! " << endl; getchar(); }
	    
	    if ( dr!=0 ) t = seg_size/dr;
	    	    	    
	    Double_t xp = x0+t*dx;
	    
	    Double_t yp = y0+t*dy;
	    
	    Double_t zp = z0+t*dz;
	    	    	    	    
	    segx->push_back( xp ); segy->push_back( yp ); segz->push_back( zp );
	    
	    x_seg[ n_seg ] = xp; y_seg[ n_seg ] = yp; z_seg[ n_seg ] = zp; n_seg++; 
	    	    
	    x0 = xp; y0 = yp; z0 = zp;
	    
	    i=( i-1 );
	    
	    // ..
	    
	    Sz += z0; Szz += z0*z0;
	    
	    Szx += z0*x0; Sx += x0;
	    
	    Szy += z0*y0; Sy += y0;
	    
	    ntot++;
	    
	    Double_t ax = ( Szx*ntot-Sz*Sx )/( Szz*ntot-Sz*Sz );
	      
	    Double_t ay = ( Szy*ntot-Sz*Sy )/( Szz*ntot-Sz*Sz );
	    
	    segnx->push_back( ax ); segny->push_back( ay ); segnz->push_back( 1.0 );
	    
	    // Double_t angx = find_angle( 1.0, ax ); Double_t angy = find_angle( 1.0, ay );
	    
	    // cout << angx*0.001*180.0/3.14 << endl;
	    
	    Sz = 0; Szz = 0; 
	    
	    Szx = 0; Sx = 0; 
	    
	    Szy = 0; Sy = 0;
	    
	    Sz += z0; Szz += z0*z0;
	     
	    Szx += z0*x0; Sx += x0;
	    
	    Szy += z0*y0; Sy += y0;
	    
	    ntot++;
	    	    
	  }
	
	if ( n_seg>=( stopper+1.0 ) && stop!=-1 ) break;
	
      }
        
    gr_seg_xyz = new TPolyLine3D( n_seg, z_seg, x_seg, y_seg );
    
    gr_seg_yz = new TGraph( n_seg, z_seg, y_seg ); gr_seg_xz = new TGraph( n_seg, z_seg, x_seg ); gr_seg_xy = new TGraph( n_seg, x_seg, y_seg );
    
    cout << " - Segmented tracks generated ! Segment size : " << seg_size << " cm " << endl;
    
    cout << "" << endl;
    
    return 0;
    
  }
  
  void TrackMomentumCalculator::GetDeltaThetaRMS( Double_t &mean, Double_t &rms, Double_t &rmse, Double_t thick )
  {
    mean = 666.0; rms = 666.0; rmse = 666.0;
    
    return;
    
  }
  
  Double_t TrackMomentumCalculator::GetMomentumMultiScatterChi2( art::Ptr<recob::Track> &trk )
  {
    cout << " Nothing will come of nothing, speak again ! " << endl;
    
    Double_t p = -1.0; 
    
    std::vector<Float_t> *recoX = new std::vector<Float_t>; recoX->clear();
    
    std::vector<Float_t> *recoY = new std::vector<Float_t>; recoY->clear();
    
    std::vector<Float_t> *recoZ = new std::vector<Float_t>; recoZ->clear();
        
    steps->clear();
    
    for ( Int_t i=1; i<=nsteps; i++ ) { steps->push_back( do_steps*i ); }
    
    Double_t L = trk->Length( 0 ); 
    
    if ( L==0.0 ) return -1.0;
    
    Int_t npoints = trk->NumberTrajectoryPoints();
        
    for ( Int_t i=0; i<npoints; i++ )
      {
	const TVector3 &muP = trk->LocationAtPoint( i );
	
	recoX->push_back( muP.X() ); recoY->push_back( muP.Y() ); recoZ->push_back( muP.Z() ); 
	
	// cout << " muX, Y, Z : " << muP.X() << ", " << muP.Y() << ", " << muP.Z() << endl;
	
      }
    
    seg_size = do_steps;
    
    Int_t ch = GetSegTracks( recoX, recoY, recoZ );
    
    if ( ch!=0 ) return -1.0;
    
    Int_t seg_steps = segx->size();
    
    Int_t seg_steps0 = seg_steps-1;
    
    Double_t recoL = seg_steps0*seg_size; cout << recoL << endl;
    
    Double_t mean = 666.0; Double_t rms = 666.0; Double_t rmse = 666.0;
    
    n_gr = 0; Double_t max1=-999.0; Double_t min1=+999.0;
	      
    for ( Int_t j=0; j<nsteps; j++ )
      {
	Double_t trial = steps->at( j );
	
	GetDeltaThetaRMS( mean, rms, rmse, trial );
	
	xmeas[ n_gr ] = trial;
		  
	ymeas[ n_gr ] = rms;
		
	eymeas[ n_gr ] = sqrt( pow( rmse, 2.0 ) + pow( 0.05*rms, 2.0 ) ); // <--- conservative syst. error to fix chi^{2} behaviour !!!
	
	// ymeas[ n_gr ] = 10.0; eymeas[ nxy ] = 1.0; // <--- for debugging !
	
	if ( min1>ymeas[ n_gr ] ) min1=ymeas[ n_gr ];
	
	if ( max1<ymeas[ n_gr ] ) max1=ymeas[ n_gr ];
	
	n_gr++;
      
      }
    
    gr_meas = new TGraphErrors( n_gr, xmeas, ymeas, 0, eymeas );
	      
    gr_meas->SetTitle( "(#Delta#theta)_{rms} versus material thickness; Material thickness in cm; (#Delta#theta)_{rms} in mrad" );
    
    gr_meas->SetLineColor( kBlack ); gr_meas->SetMarkerColor( kBlack ); gr_meas->SetMarkerStyle( 20 ); gr_meas->SetMarkerSize( 1.2 ); 
    
    gr_meas->GetXaxis()->SetLimits( ( steps->at( 0 )-steps->at( 0 ) ), ( steps->at( nsteps-1 )+steps->at( 0 ) ) );
    
    gr_meas->SetMinimum( 0.0 );
    
    gr_meas->SetMaximum( 1.80*max1 );
    
    // c1->cd();
    
    gr_meas->Draw( "APEZ" );
    
    // c1->Update();
    
    // c1->WaitPrimitive();
    
    ROOT::Minuit2::Minuit2Minimizer *mP = new ROOT::Minuit2::Minuit2Minimizer( );
    
    ROOT::Math::Functor FCA( &trkf::TrackMomentumCalculator::MyMCSChi2, 2 ); 
    
    mP->SetFunction( FCA );
    
    mP->SetLimitedVariable( 0, "p_{reco}", 1.0, 0.001, 0.001, 7.5 ); 
    	      
    mP->SetLimitedVariable( 1, "#delta#theta_{0}", 0.0, 1.0, 0, 500.0 );
    	      
    mP->SetMaxFunctionCalls( 1.E9 );
    
    mP->SetMaxIterations( 1.E9 );
    
    mP->SetTolerance( 0.01 );
    
    mP->SetStrategy( 2 );
    
    mP->SetErrorDef( 1.0 );
    /*
    mP->Minimize();
    
    mP->Hesse();
    
    // const double *pars = mP->X();  
    
    // const double *erpars = mP->Errors(); 
    
    delete mP;
    */
    delete recoX; delete recoY; delete recoZ; 
    
    cout << " Speak again ! " << endl;
        
    return p;
        
  }
  

  
} // namespace track
