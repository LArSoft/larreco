////////////////////////////////////////////////////////////////////////
///
/// \file   Track3DKalmanHitAlg.h
///
/// \brief  Track3DKalmanHit
///
/// \author H. Greenlee
///
/// Configuration parameters:
///
/// Trace         - Trace flag.
/// MaxPErr       - Maximum pointing error for free propagation.
/// GoodPErr      - Pointing error threshold for switching to free propagation.
/// MaxIncChisq   - Maximum incremental chisquare to accept a hit.
/// MaxSeedIncChisq - Maximum incremental chisquare to accept in seed phase.
/// MaxSmoothIncChisq - Maximum incremental chisquare to accept in smooth phase.
/// MaxEndChisq   - Maximum incremental chisquare for endpoint hit.
/// MinLHits      - Minimum number of hits to turn off linearized propagation.
/// MaxLDist      - Maximum distance for linearized propagation.
/// MaxPredDist   - Maximum prediciton distance to accept a hit.
/// MaxSeedPredDist - Maximum prediciton distance to accept a hit in seed phase.
/// MaxPropDist   - Maximum propagation distance to candidate surface.
/// MinSortDist   - Sort low distance threshold.
/// MaxSortDist   - Sort high distance threshold.
/// MaxSamePlane  - Maximum consecutive hits in same plane.
/// GapDist       - Minimum gap distance.
/// MaxNoiseHits  - Maximum number of hits in noise cluster.
/// MinSampleDist - Minimum sample distance (for momentum measurement).
/// FitMomRange   - Fit momentum using range.
/// FitMomMS      - Fit momentum using multiple scattering.
/// GTrace        - Graphical trace flag.
/// GTraceWW      - Graphical trace window width (pixels).
/// GTraceWH      - Graphical trace window height (pixels).
/// GTraceXMin    - Graphical trace minimum x (same for each view).
/// GTraceXMax    - Graphical trace maximum x (same for each view).
/// GTraceZMin    - Graphical trace minimum z (vector).
/// GTraceZMax    - Graphical trace maximum z (vector).
///
////////////////////////////////////////////////////////////////////////

#ifndef TRACK3DKALMANHITALG_H
#define TRACK3DKALMANHITALG_H

#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"
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
#include "lardata/RecoObjects/KHitContainerWireLine.h"
#include "lardata/RecoObjects/KHitContainerWireX.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/PropAny.h"
#include "lardata/RecoObjects/KHit.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TH1F.h"





namespace trkf {
   //SS: move this to its own header
   struct KalmanInput
   {
      art::Ptr<recob::PFParticle> pfPartPtr;
      art::PtrVector<recob::Hit> hits;
      art::PtrVector<recob::Seed> seeds;
      std::vector<art::PtrVector<recob::Hit> > seedhits;
      std::deque<trkf::KGTrack> tracks;
   };
   struct KalmanOutput {
      std::deque<trkf::KGTrack> tracks;
   };
   
   class Propagator;
   class Track3DKalmanHitAlg {
   public:
      
      /// Constructor.
      explicit Track3DKalmanHitAlg(const fhicl::ParameterSet& pset);
      
      /// Reconfigure method.
      void reconfigure(const fhicl::ParameterSet& pset);
      
      // Private member functions that do not use art::event or Assns
      // but use art::PtrVector and art::Ptr.
      void makeTracks(std::list<trkf::KalmanInput> &kalman_inputs);
      void fetchSeeds(bool first,
                      bool &pfseed,
                      KalmanInput& kalman_input,
                      art::PtrVector<recob::Hit>& unusedhits,
                      std::vector<recob::Seed>& seeds,
                      std::vector<art::PtrVector<recob::Hit>>& hitsperseed);
      recob::Seed makeSeed(const art::PtrVector<recob::Hit>& hits) const;
      void growSeedsIntoTracks(bool pfseed,
                               std::vector<recob::Seed>& seeds,
                               std::vector<art::PtrVector<recob::Hit> >& hitsperseed,
                               art::PtrVector<recob::Hit>& unusedhits,
                               art::PtrVector<recob::Hit>& hits,
                               std::deque<KGTrack>& kalman_tracks);
      void growSeedIntoTracks(bool pfseed,
                              const recob::Seed& seed,
                              const art::PtrVector<recob::Hit>& hpsit,
                              art::PtrVector<recob::Hit>& unusedhits,
                              art::PtrVector<recob::Hit>& hits,
                              std::deque<KGTrack>& kgtracks);
      void chopHitsOffSeeds(art::PtrVector<recob::Hit>const & hpsit,
                            bool pfseed,
                            art::PtrVector<recob::Hit> &seedhits) const;
      bool testSeedSlope(const double *dir) const;
      std::shared_ptr<Surface> makeSurface(const recob::Seed &seed, double *dir) const;
      bool makeKalmanTracks(const std::shared_ptr<trkf::Surface> psurf,
                            const Surface::TrackDirection trkdir,
                            art::PtrVector<recob::Hit>& seedhits,
                            art::PtrVector<recob::Hit>& hits,
                            std::deque<KGTrack>& kalman_tracks);
      
      bool smoothandextendTrack(KGTrack &trg0,
                                const art::PtrVector<recob::Hit> hits,
                                unsigned int prefplane,
                                std::deque<KGTrack>& kalman_tracks);
      void filterHitsOnKalmanTrack(const KGTrack& trg,
                                   art::PtrVector<recob::Hit>& hits,
                                   art::PtrVector<recob::Hit>& seederhits) const;
      std::unique_ptr<KHitContainer> fillHitContainer(const art::PtrVector<recob::Hit> &hits) const;
      
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
      //SS: Recommendation use std::unique_ptr<const Propagator> fProp; instead
      //this has a ripple affect, this is for later then
      const Propagator *fProp;
      // Statistics.
      
      int fNumTrack;    ///< Number of tracks produced.
      
   };
}

#endif
