#include "larreco/TrackFinder/DummyTrackMaker.h"
#include "larreco/RecoAlg/TrackCreationBookKeeper.h"
#include "art/Utilities/ToolMacros.h"

bool trkmkr::DummyTrackMaker::makeTrack(const recob::TrackTrajectory& traj, const int tkID, const std::vector<art::Ptr<recob::Hit> >& inHits,
					recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, OptionalOutputs& optionals,
					const art::Event& e) const {
  //
  TrackCreationBookKeeper tcbk(outTrack, outHits, optionals, tkID, -1, false);
  for (unsigned int i=0;i<inHits.size();++i) {
    tcbk.addPoint(traj.LocationAtPoint(i), traj.DirectionAtPoint(i), inHits[i], traj.FlagsAtPoint(i), 1);
  }
  recob::tracking::SMatrixSym55 cV;
  recob::tracking::SMatrixSym55 cE;
  tcbk.finalizeTrack(cV, cE);
  return true;
  //
}

DEFINE_ART_CLASS_TOOL(trkmkr::DummyTrackMaker)
