/// \file  TrackMomentumCalculator.h
/// \author  sowjanyag@phys.ksu.edu

#ifndef TrackMomentumCalculator_H
#define TrackMomentumCalculator_H

#include "iostream"
#include "vector"
#include "TMath.h"
#include "RecoBase/Track.h"

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

     Double_t seg_size; 

     Int_t get_seg_tracks( std::vector<Float_t> *xxx, std::vector<Float_t> *yyy, std::vector<Float_t> *zzz );

     Double_t GetMultiScatterChi2( recob::Track *trk );
          
     
   };
   
} //namespace trkf

#endif // TrackMomentumCalculator_H         
	
