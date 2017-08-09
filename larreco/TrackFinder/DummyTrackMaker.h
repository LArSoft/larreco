#ifndef DUMMYTRACKMAKER_H
#define DUMMYTRACKMAKER_H

#include "fhiclcpp/ParameterSet.h"
#include "larreco/TrackFinder/TrackMaker.h"

namespace trkmkr {

  class DummyTrackMaker : public TrackMaker {

  public:
    explicit DummyTrackMaker(fhicl::ParameterSet const& p) {}

    bool makeTrack(const recob::Trajectory& traj, const std::vector<recob::TrajectoryPointFlags>& flags,
		   const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
		   recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals,
		   const art::Event& e) const override;
  };

}

#endif
