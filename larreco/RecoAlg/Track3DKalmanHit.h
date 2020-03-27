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

#include <deque>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardata/RecoObjects/KGTrack.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Track.h"

namespace trkf {
  struct KalmanInput {
    art::Ptr<recob::PFParticle> pfPartPtr;
    art::PtrVector<recob::Hit> hits;
    art::PtrVector<recob::Seed> seeds;
    std::vector<art::PtrVector<recob::Hit>> seedhits;

    KalmanInput(){};
    explicit KalmanInput(art::PtrVector<recob::Hit>&& h) : hits(std::move(h)){};
  };
  struct KalmanOutput {
    std::deque<trkf::KGTrack> tracks;
  };

  typedef std::vector<KalmanInput> KalmanInputs;
  typedef art::PtrVector<recob::Hit> Hits;
}

#endif
