////////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterAlg.h
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee 
///
/// Configuration parameters:
///
/// Trace         - Trace flag.
/// MaxPErr       - Maximum pointing error for free propagation.
/// GoodPErr      - Pointing error threshold for switching to free propagation.
/// MaxIncChisq   - Maximum incremental chisquare to accept a hit.
/// MaxEndChisq   - Maximum incremental chisquare for endpoint hit.
/// MinLHits      - Minimum number of hits to turn off linearized propagation.
/// MaxLDist      - Maximum distance for linearized propagation.
/// MaxPredDist   - Maximum prediciton distance to accept a hit.
/// MaxPropDist   - Maximum propagation distance to candidate surface.
/// MinSortDist   - Sort low distance threshold.
/// MaxSortDist   - Sort high distance threshold.
/// MaxSamePlane  - Maximum consecutive hits in same plane.
/// GapDist       - Minimum gap distance.
/// MaxNoiseHits  - Maximum number of hits in noise cluster.
/// MinSampleDist - Minimum sample distance (for momentum measurement).
/// FitMomRange   - Fit momentum using range.
/// FitMomMS      - Fit momentum using multiple scattering.
///
////////////////////////////////////////////////////////////////////////

#ifndef KALMANFILTERALG_H
#define KALMANFILTERALG_H

#include "RecoObjects/KGTrack.h"
#include "RecoObjects/Propagator.h"
#include "RecoObjects/KHitContainer.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"

namespace trkf {

  class KalmanFilterAlg {
  public:

    /// Constructor.
    KalmanFilterAlg(const fhicl::ParameterSet& pset);

    /// Destructor.
    ~KalmanFilterAlg();

    /// Reconfigure method.
    void reconfigure(const fhicl::ParameterSet& pset);

    // Accessors.

    bool getTrace() const {return fTrace;}      ///< Trace config parameters.
    int getPlane() const {return fPlane;}       ///< Preferred view plane.

    // Modifiers.

    void setTrace(bool trace) {fTrace = trace;} ///< Set trace config parameter.
    void setPlane(int plane) {fPlane = plane;}  ///< Set preferred view plane.

    // Methods.

    /// Make a new track.
    bool buildTrack(const KTrack& trk,                     // Starting track.
		    KGTrack& trg,                          // Result global track.
		    const Propagator* prop,                // Propagator.
		    const Propagator::PropDirection dir,   // Direction.
		    KHitContainer& hits) const;            // Candidate measurements.

    /// Smooth track.
    bool smoothTrack(KGTrack& trg,                         // Global track to be smoothed.
		     KGTrack* trg1,                        // Result of unidirectional fit.
		     const Propagator* prop) const;        // Propagator.

    /// Add hits to existing track.
    bool extendTrack(KGTrack& trg,                         // Global track.
		     const Propagator* prop,               // Propagator.
		     KHitContainer& hits) const;           // Candidate measurements.

    /// Estimate track momentum using range.
    bool fitMomentumRange(const KGTrack& trg,              // Global track.
			  const Propagator* prop,          // Propagator.
			  KETrack& tremom) const;          // Track with updated momentum.

    /// Estimate track momentum using multiple scattering.
    bool fitMomentumMS(const KGTrack& trg,                 // Global track.
		       const Propagator* prop,             // Propagator.
		       KETrack& tremom) const;             // Track with updated momentum.

    /// Estimate track momentum using either range or multiple scattering.
    bool fitMomentum(const KGTrack& trg,                   // Global track.
		     const Propagator* prop,               // Propagator.
		     KETrack& tremom) const;               // Track with updated momentum.

    /// Set track momentum at each track surface.
    bool updateMomentum(const KETrack& tremom,             // Track with momentum estimate.
			const Propagator* prop,            // Propagator.
			KGTrack& trg) const;               // Global track to be updated.

    /// Iteratively smooth a track.
    bool smoothTrackIter(int niter,                        // Number of iterations.
			 KGTrack& trg,                     // Global track.
			 const Propagator* prop) const;    // Propagator.

    /// Clean track by removing noise hits near endpoints.
    void cleanTrack(KGTrack& trg) const;

  private:

    // Fcl parameters.

    bool fTrace;             ///< Trace flag.
    double fMaxPErr;         ///< Maximum pointing error for free propagation.
    double fGoodPErr;        ///< Pointing error threshold for switching to free propagation.
    double fMaxIncChisq;     ///< Maximum incremental chisquare to accept a hit.
    double fMaxEndChisq;     ///< Maximum incremental chisquare for endpoint hit.
    int fMinLHits;           ///< Minimum number of hits to turn off linearized propagation.
    double fMaxLDist;        ///< Maximum distance for linearized propagation.
    double fMaxPredDist;     ///< Maximum prediciton distance to accept a hit.
    double fMaxPropDist;     ///< Maximum propagation distance to candidate surface.
    double fMinSortDist;     ///< Sort low distance threshold.
    double fMaxSortDist;     ///< Sort high distance threshold.
    int fMaxSamePlane;       ///< Maximum consecutive hits in same plane.
    double fGapDist;         ///< Minimum gap distance.
    int fMaxNoiseHits;       ///< Maximum number of hits in noise cluster.
    double fMinSampleDist;   ///< Minimum sample distance (for momentum measurement).
    bool fFitMomRange;       ///< Fit momentum using range.
    bool fFitMomMS;          ///< Fit momentum using multiple scattering.

    // Other attributes.

    int fPlane;          ///< Preferred view plane.
  };
}

#endif
