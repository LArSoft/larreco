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

#include <deque>
#include <memory>
#include <stddef.h>
#include <vector>

#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/fwd.h"

#include "lardata/RecoObjects/KGTrack.h"
#include "lardata/RecoObjects/KHitContainer.h"
#include "lardata/RecoObjects/PropAny.h"
#include "lardata/RecoObjects/Surface.h"
#include "lardataobj/RecoBase/Seed.h"
#include "larreco/RecoAlg/KalmanFilterAlg.h"
#include "larreco/RecoAlg/SeedFinderAlgorithm.h"
#include "larreco/RecoAlg/Track3DKalmanHit.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace trkf {
  class KHitContainer;
  class Propagator;
}

namespace trkf {
  class Track3DKalmanHitAlg {
  public:
    explicit Track3DKalmanHitAlg(const fhicl::ParameterSet& pset);

    std::vector<trkf::KalmanOutput> makeTracks(detinfo::DetectorClocksData const& clockData,
                                               detinfo::DetectorPropertiesData const& detProp,
                                               KalmanInputs& kalman_inputs);
    void fetchPFParticleSeeds(const art::PtrVector<recob::Seed>& pfseeds,
                              const std::vector<Hits>& pfseedhits,
                              std::vector<recob::Seed>& seeds,
                              std::vector<Hits>& hitsperseed) const;
    recob::Seed makeSeed(detinfo::DetectorPropertiesData const& detProp, const Hits& hits) const;
    void growSeedsIntoTracks(detinfo::DetectorPropertiesData const& detProp,
                             const bool pfseed,
                             const std::vector<recob::Seed>& seeds,
                             const std::vector<Hits>& hitsperseed,
                             Hits& unusedhits,
                             Hits& hits,
                             std::deque<KGTrack>& kalman_tracks);
    void growSeedIntoTracks(detinfo::DetectorPropertiesData const& detProp,
                            const bool pfseed,
                            const recob::Seed& seed,
                            const Hits& hpsit,
                            Hits& unusedhits,
                            Hits& hits,
                            std::deque<KGTrack>& kgtracks);
    void chopHitsOffSeeds(Hits const& hpsit, bool pfseed, Hits& seedhits) const;
    bool testSeedSlope(const double* dir) const;
    std::shared_ptr<Surface> makeSurface(const recob::Seed& seed, double* dir) const;
    bool makeKalmanTracks(detinfo::DetectorPropertiesData const& detProp,
                          const std::shared_ptr<trkf::Surface> psurf,
                          const Surface::TrackDirection trkdir,
                          Hits& seedhits,
                          Hits& hits,
                          std::deque<KGTrack>& kalman_tracks);
    bool smoothandextendTrack(detinfo::DetectorPropertiesData const& detProp,
                              Propagator const& propagator,
                              KGTrack& trg0,
                              const Hits hits,
                              unsigned int prefplane,
                              std::deque<KGTrack>& kalman_tracks);
    bool extendandsmoothLoop(detinfo::DetectorPropertiesData const& detProp,
                             Propagator const& propagator,
                             KGTrack& trg1,
                             unsigned int prefplane,
                             Hits& trackhits) const;
    void filterHitsOnKalmanTrack(const KGTrack& trg, Hits& hits, Hits& seederhits) const;
    std::unique_ptr<KHitContainer> fillHitContainer(detinfo::DetectorPropertiesData const& detProp,
                                                    const Hits& hits) const;

    bool qualityCutsOnSeedTrack(const KGTrack& trg0) const;

    void fitnupdateMomentum(Propagator const& propagator, KGTrack& trg1, KGTrack& trg2) const;

  private:
    // Fcl parameters.
    bool fDoDedx;            ///< Global dE/dx enable flag.
    bool fSelfSeed;          ///< Self seed flag.
    double fMaxTcut;         ///< Maximum delta ray energy in MeV for restricted dE/dx.
    bool fLineSurface;       ///< Line surface flag.
    size_t fMinSeedHits;     ///< Minimum number of hits per track seed.
    int fMinSeedChopHits;    ///< Potentially chop seeds that exceed this length.
    int fMaxChopHits;        ///< Maximum number of hits to chop from each end of seed.
    double fMaxSeedChiDF;    ///< Maximum seed track chisquare/dof.
    double fMinSeedSlope;    ///< Minimum seed slope (dx/dz).
    double fInitialMomentum; ///< Initial (or constant) momentum.

    // Algorithm objects.

    KalmanFilterAlg fKFAlg;             ///< Kalman filter algorithm.
    SeedFinderAlgorithm fSeedFinderAlg; ///< Seed finder.

    // Statistics.
    int fNumTrack; ///< Number of tracks produced.
  };
}

#endif
