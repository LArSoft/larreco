////////////////////////////////////////////////////////////////////////
///
/// \file   Track3DKalmanHitAlg.h
///
/// \brief  Track3DKalmanHit Algorithm
///
/// \author
///
////////////////////////////////////////////////////////////////////////

#ifndef TRACK3DKALMANHIT_H
#define TRACK3DKALMANHIT_H

#include <vector>
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"


#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Seed.h"


namespace trkf {
   struct KalmanInput
   {
      art::Ptr<recob::PFParticle> pfPartPtr;
      art::PtrVector<recob::Hit> hits;
      art::PtrVector<recob::Seed> seeds;
      std::vector<art::PtrVector<recob::Hit>> seedhits;
      
      KalmanInput() {};
      explicit KalmanInput(art::PtrVector<recob::Hit>&& h): hits(std::move(h)){};
      
   };
   struct KalmanOutput {
      std::deque<trkf::KGTrack> tracks;
   };
   
   typedef typename std::vector<KalmanInput> KalmanInputs;
}

#endif
