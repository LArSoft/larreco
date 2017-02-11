#ifndef TRACKKALMANFITTER_H
#define TRACKKALMANFITTER_H

#include "canvas/Persistency/Common/Ptr.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoObjects/PropagatorToPlane.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

namespace recob {
  class Track;
  class Hit;
}

class TVector3;

namespace trkf {

  /**
   * @brief Fit tracks using Kalman Filter fit+smooth.
   *
   * This algorithm fits tracks using a Kalman Filter forward fit followed by a backward smoothing.
   *
   * Inputs are: track object, associated hits, and momentum estimate.
   * Output are: resulting track and associated hits. The resulting track will feature covariance matrices at start and end positions.
   *
   * Note: the fit assumes the track direction to be towards increasing z.
   *
   */

  class TrackKalmanFitter {

  public:
    TrackKalmanFitter(const PropagatorToPlane* prop, bool useRMS, bool sortHitsByPlane, bool sortOutputHitsMinLength, bool skipNegProp, float hitErrScaleFact){
      propToPlane=prop;
      useRMS_=useRMS;
      sortHitsByPlane_=sortHitsByPlane;
      sortOutputHitsMinLength_=sortOutputHitsMinLength;
      skipNegProp_=skipNegProp;
      hitErrScaleFact_=hitErrScaleFact;
      detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    }

    bool fitTrack(const recob::Track& inputTrack, const std::vector<art::Ptr<recob::Hit> >& hits,
		  const double pval, const int pdgid, const bool flipDirection,
		  recob::Track& outputTrack, art::PtrVector<recob::Hit>& outputHits,
		  std::vector<recob::TrackFitHitInfo>& trackFitHitInfos);

    bool getSkipNegProp() const     { return skipNegProp_; }
    void setSkipNegProp(bool value) { skipNegProp_=value; }

  private:

    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::DetectorProperties* detprop;
    const PropagatorToPlane* propToPlane;
    bool useRMS_;
    bool sortHitsByPlane_;
    bool sortOutputHitsMinLength_;
    bool skipNegProp_;
    float hitErrScaleFact_;
  };

}

#endif
