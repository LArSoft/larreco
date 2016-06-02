////////////////////////////////////////////////////////////////////////
///
/// \file   Track3DKalmanHitAlg.h
///
/// \brief  Track3DKalmanHit Algorithm
///
/// \author
///
////////////////////////////////////////////////////////////////////////
// Configuration parameters:
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// DoDedx             - Global dE/dx enable flag.
// SelfSeed           - Self seed flag.
// LineSurface        - Hits on line surfaces (true) or plane surfaces (false).
// MinSeedHits        - Minimum number of hits per track seed.
// MinSeedChopHits:   - Potentially chop seeds that exceed this length.
// MaxChopHits        - Maximum number of hits to chop from each end of seed.
// MaxSeedChiDF       - Maximum seed track chisquare/dof.
// MinSeedSlope       - Minimum seed slope (dx/dz).
// InitialMomentum    - Initial momentum guess.
// KalmanFilterAlg    - Parameter set for KalmanFilterAlg.
// SeedFinderAlg      - Parameter set for seed finder algorithm object.
////////////////////////////////////////////////////////////////////////


#ifndef TRACK3DKALMANHITALG_H
#define TRACK3DKALMANHITALG_H

#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Seed.h"
#include "larreco/RecoAlg/KalmanFilterAlg.h"
#include "larreco/RecoAlg/SeedFinderAlgorithm.h"
#include "larreco/RecoAlg/Track3DKalmanHit.h"
#include "lardata/RecoObjects/KHitContainerWireLine.h"
#include "lardata/RecoObjects/KHitContainerWireX.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/PropAny.h"
#include "lardata/RecoObjects/KHit.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TH1F.h"

namespace trkf {
   class Propagator;
   class Track3DKalmanHitAlg {
   public:
      
      /// Constructor.
      explicit Track3DKalmanHitAlg(const fhicl::ParameterSet& pset);
      
      /// Reconfigure method.
      void reconfigure(const fhicl::ParameterSet& pset);
      
      // Private member functions that do not use art::event or Assns
      // but use art::PtrVector and art::Ptr.
      std::vector<trkf::KalmanOutput> makeTracks(KalmanInputs &kalman_inputs);
      bool fetchPFParticleSeeds(const art::PtrVector<recob::Seed> &pfseeds,
                                const std::vector<Hits> &pfseedhits,
                                std::vector<recob::Seed>& seeds,
                                std::vector<Hits>& hitsperseed);
      recob::Seed makeSeed(const Hits& hits) const;
      void growSeedsIntoTracks(const bool pfseed,
                               const std::vector<recob::Seed>& seeds,
                               const std::vector<Hits >& hitsperseed,
                               Hits& unusedhits,
                               Hits& hits,
                               std::deque<KGTrack>& kalman_tracks);
      void growSeedIntoTracks(const bool pfseed,
                              const recob::Seed& seed,
                              const Hits& hpsit,
                              Hits& unusedhits,
                              Hits& hits,
                              std::deque<KGTrack>& kgtracks);
      void chopHitsOffSeeds(Hits const & hpsit,
                            bool pfseed,
                            Hits &seedhits) const;
      bool testSeedSlope(const double *dir) const;
      std::shared_ptr<Surface> makeSurface(const recob::Seed &seed, double *dir) const;
      bool makeKalmanTracks(const std::shared_ptr<trkf::Surface> psurf,
                            const Surface::TrackDirection trkdir,
                            Hits& seedhits,
                            Hits& hits,
                            std::deque<KGTrack>& kalman_tracks);
      bool smoothandextendTrack(KGTrack &trg0,
                                const Hits hits,
                                unsigned int prefplane,
                                std::deque<KGTrack>& kalman_tracks);
      bool extendandsmoothLoop(
                               KGTrack &trg1,
                               unsigned int prefplane,
                               Hits &trackhits);
      void filterHitsOnKalmanTrack(const KGTrack& trg,
                                   Hits& hits,
                                   Hits& seederhits) const;
      std::unique_ptr<KHitContainer> fillHitContainer(const Hits &hits) const;
      
      bool qualityCutsOnSeedTrack(const KGTrack &trg0) const;
      
      void fitnupdateMomentum(KGTrack& trg1,
                              KGTrack& trg2);
      
   private:
      
      // Fcl parameters.
      bool fDoDedx;                       ///< Global dE/dx enable flag.
      bool fSelfSeed;                     ///< Self seed flag.
      double fMaxTcut;                    ///< Maximum delta ray energy in MeV for restricted dE/dx.
      bool fLineSurface;                  ///< Line surface flag.
      size_t fMinSeedHits;                   ///< Minimum number of hits per track seed.
      int fMinSeedChopHits;               ///< Potentially chop seeds that exceed this length.
      int fMaxChopHits;                   ///< Maximum number of hits to chop from each end of seed.
      double fMaxSeedChiDF;               ///< Maximum seed track chisquare/dof.
      double fMinSeedSlope;               ///< Minimum seed slope (dx/dz).
      double fInitialMomentum;            ///< Initial (or constant) momentum.
      
      // Algorithm objects.
      
      KalmanFilterAlg fKFAlg;             ///< Kalman filter algorithm.
      SeedFinderAlgorithm fSeedFinderAlg; ///< Seed finder.
      
      /// Propagator.
      std::unique_ptr<const Propagator> fProp;
      
      // Statistics.
      int fNumTrack;    ///< Number of tracks produced.
      
   };
}

#endif
