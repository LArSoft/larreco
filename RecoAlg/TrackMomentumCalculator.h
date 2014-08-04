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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

namespace trkf{

   class TrackMomentumCalculator
   {
   public:
     
     // Constructor and destructor 
     
     TrackMomentumCalculator() 
       {
	 do_steps = 5.0; nsteps = 9; 
	 
	 steps = new std::vector<Float_t>; steps->clear();
	 
	 for ( Int_t i=1; i<=nsteps; i++ ) { steps->push_back( do_steps*i ); }
	 
	 stop = -1.0; 
	 
	 n_seg = 0; 
	 
	 gr_seg_xyz = new TPolyLine3D(); gr_seg_xy = new TGraph(); gr_seg_yz = new TGraph(); gr_seg_xz = new TGraph(); 
	 
	 segx = new std::vector<Float_t>; segy = new std::vector<Float_t>; segz = new std::vector<Float_t>;
	 
	 segnx = new std::vector<Float_t>; segny = new std::vector<Float_t>; segnz = new std::vector<Float_t>;
	 
	 do_steps2 = 10.0;
	 
       }
     
     virtual ~TrackMomentumCalculator() {}
     
     double GetTrackMomentum(double trkrange, int pdg);
          
     // MultiScatter business ...
     
     // Author: Leonidas N. Kalousis (August 2014)
     
     Double_t do_steps; 
     
     Int_t nsteps; std::vector<Float_t> *steps;
     
     Double_t seg_size; Double_t stop; Int_t n_seg;

     Double_t x_seg[100000]; Double_t y_seg[100000]; Double_t z_seg[100000];

     TPolyLine3D *gr_seg_xyz; TGraph *gr_seg_xy; TGraph *gr_seg_yz; TGraph *gr_seg_xz; 
  
     std::vector<Float_t> *segx; std::vector<Float_t> *segy; std::vector<Float_t> *segz; 
     
     std::vector<Float_t> *segnx; std::vector<Float_t> *segny; std::vector<Float_t> *segnz;

     Int_t GetSegTracks( std::vector<Float_t> *xxx, std::vector<Float_t> *yyy, std::vector<Float_t> *zzz );
     
     Double_t GetMomentumMultiScatterChi2( art::Ptr<recob::Track> &trk );
     
     Double_t do_steps2; 
     
   };
   
} //namespace trkf

#endif // TrackMomentumCalculator_H         
	
