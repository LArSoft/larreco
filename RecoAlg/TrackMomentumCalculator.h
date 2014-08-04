/// \file  TrackMomentumCalculator.h
/// \author  sowjanyag@phys.ksu.edu

#ifndef TrackMomentumCalculator_H
#define TrackMomentumCalculator_H

#include "iostream"
#include "vector"
#include "TMath.h"
#include "RecoBase/Track.h"
#include "TGraph.h"
#include "TPolyLine3D.h" 

using namespace std;

namespace trkf{

   class TrackMomentumCalculator
   {
   public:
     // Constructor and destructor 
     TrackMomentumCalculator() {}
     virtual ~TrackMomentumCalculator() {}
     
     double GetTrackMomentum(double trkrange, int pdg);
     
     
     // MultiScatter business ...
     
     // Author: Leonidas N. Kalousis (August 2014)

     Double_t do_steps; Int_t nsteps; std::vector<Float_t> *steps;

     Double_t seg_size; Double_t stop; Int_t n_seg;

     Double_t x_seg[100000]; Double_t y_seg[100000]; Double_t z_seg[100000];

     TPolyLine3D *gr_seg_xyz; TGraph *gr_seg_xy; TGraph *gr_seg_yz; TGraph *gr_seg_xz; 
  
     std::vector<Float_t> *segx; std::vector<Float_t> *segy; std::vector<Float_t> *segz; 
     
     std::vector<Float_t> *segnx; std::vector<Float_t> *segny; std::vector<Float_t> *segnz;

     Int_t GetSegTracks( std::vector<Float_t> *xxx, std::vector<Float_t> *yyy, std::vector<Float_t> *zzz );
     
     Double_t GetMultiScatterChi2( recob::Track *trk );
          
     
   };
   
} //namespace trkf

#endif // TrackMomentumCalculator_H         
	
