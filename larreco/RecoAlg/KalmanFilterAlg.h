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

#ifndef KALMANFILTERALG_H
#define KALMANFILTERALG_H

#include <map>
#include <memory>
#include <vector>

#include "lardata/RecoObjects/KETrack.h"
#include "lardata/RecoObjects/KGTrack.h"
#include "lardata/RecoObjects/KTrack.h"
#include "lardata/RecoObjects/Propagator.h"

namespace fhicl {
  class ParameterSet;
}
namespace trkf {
  class KHitContainer;
}

#include "TCanvas.h"

class TMarker;
class TPaveText;
class TVirtualPad;

namespace trkf {

  class KalmanFilterAlg {
  public:
    KalmanFilterAlg(const fhicl::ParameterSet& pset);

    // Accessors.

    bool
    getTrace() const
    {
      return fTrace;
    } ///< Trace config parameters.
    int
    getPlane() const
    {
      return fPlane;
    } ///< Preferred view plane.

    // Modifiers.

    void
    setTrace(bool trace)
    {
      fTrace = trace;
    } ///< Set trace config parameter.
    void
    setPlane(int plane)
    {
      fPlane = plane;
    } ///< Set preferred view plane.

    // Methods.

    /// Make a new track.
    bool buildTrack(const KTrack& trk,                   // Starting track.
                    KGTrack& trg,                        // Result global track.
                    const Propagator& prop,              // Propagator.
                    const Propagator::PropDirection dir, // Direction.
                    KHitContainer& hits,                 // Candidate measurements.
                    bool linear) const;                  // Linear flag.

    /// Smooth track.
    bool smoothTrack(KGTrack& trg,                  // Global track to be smoothed.
                     KGTrack* trg1,                 // Result of unidirectional fit.
                     const Propagator& prop) const; // Propagator.

    /// Add hits to existing track.
    bool extendTrack(KGTrack& trg,               // Global track.
                     const Propagator& prop,     // Propagator.
                     KHitContainer& hits) const; // Candidate measurements.

    /// Estimate track momentum using range.
    bool fitMomentumRange(const KGTrack& trg,     // Global track.
                          KETrack& tremom) const; // Track with updated momentum.

    /// Estimate track momentum using multiple scattering.
    bool fitMomentumMS(const KGTrack& trg,     // Global track.
                       const Propagator& prop, // Propagator.
                       KETrack& tremom) const; // Track with updated momentum.

    /// Estimate track momentum using either range or multiple scattering.
    bool fitMomentum(const KGTrack& trg,     // Global track.
                     const Propagator& prop, // Propagator.
                     KETrack& tremom) const; // Track with updated momentum.

    /// Set track momentum at each track surface.
    bool updateMomentum(const KETrack& tremom,  // Track with momentum estimate.
                        const Propagator& prop, // Propagator.
                        KGTrack& trg) const;    // Global track to be updated.

    /// Iteratively smooth a track.
    bool smoothTrackIter(int niter,                     // Number of iterations.
                         KGTrack& trg,                  // Global track.
                         const Propagator& prop) const; // Propagator.

    /// Clean track by removing noise hits near endpoints.
    void cleanTrack(KGTrack& trg) const;

  private:
    // Fcl parameters.

    bool fTrace;               ///< Trace flag.
    double fMaxPErr;           ///< Maximum pointing error for free propagation.
    double fGoodPErr;          ///< Pointing error threshold for switching to free propagation.
    double fMaxIncChisq;       ///< Maximum incremental chisquare to accept a hit.
    double fMaxSeedIncChisq;   ///< Maximum incremental chisquare to accept a hit in seed phase.
    double fMaxSmoothIncChisq; ///< Maximum incremental chisquare to accept a hit in smooth phase.
    double fMaxEndChisq;       ///< Maximum incremental chisquare for endpoint hit.
    int fMinLHits;             ///< Minimum number of hits to turn off linearized propagation.
    double fMaxLDist;          ///< Maximum distance for linearized propagation.
    double fMaxPredDist;       ///< Maximum prediciton distance to accept a hit.
    double fMaxSeedPredDist;   ///< Maximum prediciton distance to accept a hit in seed phase.
    double fMaxPropDist;       ///< Maximum propagation distance to candidate surface.
    double fMinSortDist;       ///< Sort low distance threshold.
    double fMaxSortDist;       ///< Sort high distance threshold.
    int fMaxSamePlane;         ///< Maximum consecutive hits in same plane.
    double fGapDist;           ///< Minimum gap distance.
    int fMaxNoiseHits;         ///< Maximum number of hits in noise cluster.
    double fMinSampleDist;     ///< Minimum sample distance (for momentum measurement).
    bool fFitMomRange;         ///< Fit momentum using range.
    bool fFitMomMS;            ///< Fit momentum using multiple scattering.
    bool fGTrace;              ///< Graphical trace flag.
    double fGTraceWW;          ///< Window width.
    double fGTraceWH;          ///< Window height.
    double fGTraceXMin;        ///< Graphical trace minimum x.
    double fGTraceXMax;        ///< Graphical trace maximum x.
    std::vector<double> fGTraceZMin; ///< Graphical trace minimum z for each view.
    std::vector<double> fGTraceZMax; ///< Graphical trace maximum z for each view.

    // Other attributes.

    int fPlane;                                              ///< Preferred view plane.
    mutable std::vector<std::unique_ptr<TCanvas>> fCanvases; ///< Graphical trace canvases.
    mutable std::vector<TVirtualPad*> fPads;                 ///< View pads in current canvas.
    mutable TVirtualPad* fInfoPad;                           ///< Information pad.
    mutable TPaveText* fMessages;                            ///< Message box.
    mutable std::map<int, TMarker*> fMarkerMap;              ///< Markers in current canvas.
  };
}

#endif
