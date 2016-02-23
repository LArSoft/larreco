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



struct LocalKalmanStruct
{
   art::Ptr<recob::PFParticle> pfPartPtr;
   art::PtrVector<recob::Hit> hits;
   art::PtrVector<recob::Seed> seeds;
   std::vector<art::PtrVector<recob::Hit> > seedhits;
   std::deque<trkf::KGTrack> tracks; //SS: Why deque is used?
};

namespace trkf {
   
   class Propagator;
   class Track3DKalmanHitAlg {
   public:
      
      /// Constructor.
      Track3DKalmanHitAlg(const fhicl::ParameterSet& pset);
      
      /// Destructor.
      ~Track3DKalmanHitAlg();
      
      /// Reconfigure method.
      void reconfigure(const fhicl::ParameterSet& pset);
      
      
      // Private member functions that do not use art::event or Assns
      // but use art::PtrVector and art::Ptr.
      
      recob::Seed makeSeed(const art::PtrVector<recob::Hit>& hits) const;
      // the functions that are not modifying any class members should be const
      
      void filterHitsOnKalmanTrack(const KGTrack& trg,
                                   art::PtrVector<recob::Hit>& hits,
                                   art::PtrVector<recob::Hit>& seederhits) const;
      std::unique_ptr<KHitContainer> fillHitContainer(const art::PtrVector<recob::Hit> &hits) const;
      
      void chopHitsOffSeeds(art::PtrVector<recob::Hit>const & hpsit,
                            bool pfseed,
                            art::PtrVector<recob::Hit> &seedhits) const;
      bool qualityCutsOnSeedTrack(const KGTrack &trg0) const;
      bool smoothandextendTrack(KGTrack &trg0,
                                const art::PtrVector<recob::Hit> hits,
                                unsigned int prefplane,
                                std::deque<KGTrack>& kalman_tracks);
      void fitnupdateMomentum(KGTrack& trg1,
                              KGTrack& trg2);
      bool test1(const double *dir);
      std::shared_ptr<Surface> makeSurface(const recob::Seed &seed, double *dir);
      std::vector<KTrack> makeInitialKtracks(const std::shared_ptr<trkf::Surface> psurf);
      bool processInitialtracks(const trkf::KTrack &trk,
                                art::PtrVector<recob::Hit>& seedhits,
                                art::PtrVector<recob::Hit>& hits,
                                std::deque<KGTrack>& kalman_tracks);
      void processSeeds(bool pfseed,
                        std::vector<recob::Seed>& seeds,
                        std::vector<art::PtrVector<recob::Hit> >& hitsperseed,
                        art::PtrVector<recob::Hit>& unusedhits,
                        art::PtrVector<recob::Hit>& hits,
                        std::deque<KGTrack>& kalman_tracks);
      void generateKalmantracks(std::list<LocalKalmanStruct> &LocalKalmanStructList);
      
      
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
      const Propagator* fProp;
      
      // Statistics.
      
      int fNumTrack;    ///< Number of tracks produced.
      
   };
}

#endif
