/// \file  TrackMomentumCalculator.h
//  \author  sowjanyag@phys.ksu.edu

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
#include "RecoAlg/RootMathFunctor.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include <math.h>
#include <cmath>

using namespace std;

// Global variables/input for the TMinuit2 chi^2 minimization ..

Double_t xmeas[30]; Double_t ymeas[30]; Double_t eymeas[30]; Int_t nmeas;

// ..

namespace trkf{
  
  class TrackMomentumCalculator
  {
    Int_t n;
  
    Double_t x[10000]; Double_t y[10000]; Double_t z[10000];
        
    Int_t n_reco;
  
    Float_t x_reco[1000]; Float_t y_reco[1000]; Float_t z_reco[1000];
        
    Float_t seg_size; Float_t seg_stop; Int_t n_seg;
    
    Float_t x_seg[1000]; Float_t y_seg[1000]; Float_t z_seg[1000];
            
    TVector3 basex; TVector3 basey; TVector3 basez; 
       
    std::vector<Float_t> segx; std::vector<Float_t> segy; std::vector<Float_t> segz; 
  
    std::vector<Float_t> segnx; std::vector<Float_t> segny; std::vector<Float_t> segnz;
    
    std::vector<Float_t> segL;
    
    Double_t find_angle( Double_t vz, Double_t vy );
    
    Float_t steps_size; Int_t n_steps; std::vector<Float_t> steps;
    
    
  public:
    
    // Constructor and destructor  // 
    
    TrackMomentumCalculator();
    
    virtual ~TrackMomentumCalculator() {}
    
    double GetTrackMomentum(double trkrange, int pdg);
    
    TPolyLine3D *gr_xyz; TGraph *gr_xy; TGraph *gr_yz; TGraph *gr_xz; 
    
    Int_t GetTracks( const std::vector<Float_t> &xxx, const std::vector<Float_t> &yyy, const std::vector<Float_t> &zzz );
        
    TPolyLine3D *gr_reco_xyz; TGraph *gr_reco_xy; TGraph *gr_reco_yz; TGraph *gr_reco_xz; 
    
    Int_t GetRecoTracks( const std::vector<Float_t> &xxx, const std::vector<Float_t> &yyy, const std::vector<Float_t> &zzz );
        
    TPolyLine3D *gr_seg_xyz; TGraph *gr_seg_xy; TGraph *gr_seg_yz; TGraph *gr_seg_xz; 

    Int_t GetSegTracks2( const std::vector<Float_t> &xxx, const std::vector<Float_t> &yyy, const std::vector<Float_t> &zzz );
    
    void GetDeltaThetaRMS( Double_t &mean, Double_t &rms, Double_t &rmse, Double_t thick );
    
    TGraphErrors *gr_meas;

    Double_t GetMomentumMultiScatterChi2( const art::Ptr<recob::Track> &trk );
    
    Double_t p_mcs; Double_t p_mcs_e; Double_t chi2;
    
  };
  
  
} //namespace trkf

#endif // TrackMomentumCalculator_H  
