#ifndef TRACKKALMANFITTER_H
#define TRACKKALMANFITTER_H

#include "canvas/Persistency/Common/Ptr.h"

namespace recob {
  class Track;
  class Hit;
}

class TVector3;

namespace trkf {

  class Propagator;
  class KTrack;

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
    TrackKalmanFitter(const trkf::Propagator* prop, bool useRMS, bool sortHitsByPlane, bool sortOutputHitsMinLength, bool skipNegProp, float hitErrScaleFact){
      prop_=prop;
      useRMS_=useRMS;
      sortHitsByPlane_=sortHitsByPlane;
      sortOutputHitsMinLength_=sortOutputHitsMinLength;
      skipNegProp_=skipNegProp;
      hitErrScaleFact_=hitErrScaleFact;
    }

    bool fitTrack(const recob::Track& inputTrack, const std::vector<art::Ptr<recob::Hit> >& hits,
		  const double pval, const int pdgid, const bool flipDirection,
		  recob::Track& outputTrack, art::PtrVector<recob::Hit>& outputHits);

    trkf::KTrack convertRecobTrackIntoKTrack(const TVector3& position, const TVector3&  direction,  const double pval, const int pdgid);

  private:

    const trkf::Propagator* prop_;
    bool useRMS_;
    bool sortHitsByPlane_;
    bool sortOutputHitsMinLength_;
    bool skipNegProp_;
    float hitErrScaleFact_;
  };

}

#endif
