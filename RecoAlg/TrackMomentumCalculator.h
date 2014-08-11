/// \file  TrackMomentumCalculator.h
/// \author  sowjanyag@phys.ksu.edu

#ifndef TrackMomentumCalculator_H
#define TrackMomentumCalculator_H

#include "iostream"
#include "vector"
#include "TMath.h"
#include "RecoBase/Track.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TPolyLine3D.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/FunctionMinimum.h" 
#include "Minuit2/MnMigrad.h" 
#include "Minuit2/MnUserParameters.h" 
#include "Minuit2/MnPrint.h" 
#include "Minuit2/FCNBase.h" 
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

Double_t xmeas[100]; Double_t ymeas[100]; Double_t eymeas[100]; Int_t n_gr;

namespace trkf{

   class TrackMomentumCalculator
   {
     Double_t do_steps; Int_t nsteps; std::vector<Float_t> steps;
     
     Double_t seg_size; Double_t stop; Int_t n_seg;
     
     Double_t x_seg[100000]; Double_t y_seg[100000]; Double_t z_seg[100000];
     
     Double_t find_angle( Double_t vz, Double_t vy );
     
     void normalizer( std::vector<Double_t>& v1 );
     
     Double_t dot_prod( const std::vector<Double_t>& v1, const std::vector<Double_t>& v2 );
     
     void cross_prod( const std::vector<Double_t>& v1, const std::vector<Double_t>& v2, std::vector<Double_t>& v3 );
     
     std::vector<Double_t> basex; std::vector<Double_t> basey; std::vector<Double_t> basez;
     
   public:
     
     // Constructor and destructor 
     
     TrackMomentumCalculator();
     
     virtual ~TrackMomentumCalculator() {}
     
     double GetTrackMomentum(double trkrange, int pdg);
          
     TPolyLine3D *gr_seg_xyz; TGraph *gr_seg_xy; TGraph *gr_seg_yz; TGraph *gr_seg_xz; 
     
     std::vector<Float_t> segx; std::vector<Float_t> segy; std::vector<Float_t> segz; 
     
     std::vector<Float_t> segnx; std::vector<Float_t> segny; std::vector<Float_t> segnz;
               
     TGraphErrors *gr_meas;
          
     Int_t GetSegTracks( const std::vector<Float_t>& xxx, const std::vector<Float_t>& yyy, const std::vector<Float_t>& zzz );
     
     void GetDeltaThetaRMS( Double_t &mean, Double_t &rms, Double_t &rmse, Double_t thick );
          
     Double_t GetMomentumMultiScatterChi2(const art::Ptr<recob::Track> &trk );
          
     Double_t p_reco; 
     
     Double_t p_reco_e; 
     
     Double_t chi2;
     
   };
   
} //namespace trkf

#endif // TrackMomentumCalculator_H         
	
