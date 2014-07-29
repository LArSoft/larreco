/// \file  TrackMomentumCalculator.h
/// \author  sowjanyag@phys.ksu.edu

#ifndef TrackMomentumCalculator_H
#define TrackMomentumCalculator_H

namespace trkf{

   class TrackMomentumCalculator
   {
   public:
     // Constructor and destructor 
     TrackMomentumCalculator() {}
     virtual ~TrackMomentumCalculator() {}
     
     double GetTrackMomentum(double trkrange, int pdg);
   };
   
} //namespace trkf

#endif // TrackMomentumCalculator_H         
	
