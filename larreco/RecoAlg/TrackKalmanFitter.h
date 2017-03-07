#ifndef TRACKKALMANFITTER_H
#define TRACKKALMANFITTER_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

namespace recob {
  class Track;
  class Hit;
  class TrajectoryPointFlags;
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

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> useRMS {
	Name("useRMSError"),
	Comment("Flag to replace the default hit error SigmaPeakTime() with RMS().")
      };
      fhicl::Atom<bool> sortHitsByPlane {
        Name("sortHitsByPlane"),
        Comment("Flag to sort hits along the forward fit. The hit order in each plane is preserved, the next hit to process in 3D is chosen as the one with shorter 3D propagation distance among the next hit in all planes.")
      };
      fhicl::Atom<bool> sortOutputHitsMinLength {
        Name("sortOutputHitsMinLength"),
        Comment("Flag to decide whether the hits are sorted before creating the output track in order to avoid tracks with huge length.")
      };
      fhicl::Atom<bool> skipNegProp {
        Name("skipNegProp"),
        Comment("Flag to decide whether, during the forward fit, the hits corresponding to a negative propagation distance should be dropped.")
      };
      fhicl::Atom<float> hitErr2ScaleFact {
        Name("hitErr2ScaleFact"),
	Comment(".") //fixme
      };
      fhicl::Atom<int> dumpLevel {
        Name("dumpLevel"),
	Comment("0 for no debug printouts, 1 for moderate, 2 for maximum.")
      };
    };
    using Parameters = fhicl::Table<Config>;
    
    TrackKalmanFitter(const TrackStatePropagator* prop, bool useRMS, bool sortHitsByPlane, bool sortOutputHitsMinLength, bool skipNegProp, float hitErr2ScaleFact, int dumpLevel){
      propagator=prop;
      useRMS_=useRMS;
      sortHitsByPlane_=sortHitsByPlane;
      sortOutputHitsMinLength_=sortOutputHitsMinLength;
      skipNegProp_=skipNegProp;
      hitErr2ScaleFact_=hitErr2ScaleFact;
      dumpLevel_=dumpLevel;
      detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    }
    explicit TrackKalmanFitter(const TrackStatePropagator* prop, Parameters const & p)
      : TrackKalmanFitter(prop,p().useRMS(),p().sortHitsByPlane(),p().sortOutputHitsMinLength(),p().skipNegProp(),p().hitErr2ScaleFact(),p().dumpLevel()) {}
    
    bool fitTrack(const recob::Trajectory& track, int tkID,
		  const SMatrixSym55& covVtx, const SMatrixSym55& covEnd,
		  const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags,
		  const double pval, const int pdgid, const bool flipDirection,
		  recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits,
		  std::vector<recob::TrackFitHitInfo>& trackFitHitInfos);

    bool getSkipNegProp() const     { return skipNegProp_; }
    void setSkipNegProp(bool value) { skipNegProp_=value; }

  private:

    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::DetectorProperties* detprop;
    const TrackStatePropagator* propagator;
    bool useRMS_;
    bool sortHitsByPlane_;
    bool sortOutputHitsMinLength_;
    bool skipNegProp_;
    float hitErr2ScaleFact_;
    int dumpLevel_;
  };

}

#endif
