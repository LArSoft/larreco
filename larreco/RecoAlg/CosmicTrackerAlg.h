////////////////////////////////////////////////////////////////////////
// CosmicTrackerAlg.h
//
// tjyang@fnal.gov
//
// Use Bruce's TrackTrajectoryAlg to reconstruct track
//
// Input: a vector of recob::hits
// output: a vector of 3D points and their assocations wit the hits
//
/////////////////////////////////////////////////////////////////////////
#ifndef COSMICTRACKERALG_H
#define COSMICTRACKERALG_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/fwd.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TrackTrajectoryAlg.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
  class LArProperties;
}

#include <vector>

#include "TVector3.h"

namespace trkf {
  class CosmicTrackerAlg {
  public:
    explicit CosmicTrackerAlg(fhicl::ParameterSet const& pset);

    void SPTReco(detinfo::DetectorClocksData const& clockData,
                 detinfo::DetectorPropertiesData const& detProp,
                 std::vector<art::Ptr<recob::Hit>>& fHits);

    //trajectory position and direction returned by TrackTrajectoryAlg
    std::vector<TVector3> trajPos;
    std::vector<TVector3> trajDir;
    std::vector<std::vector<art::Ptr<recob::Hit>>> trajHit;

    //position and direction of each point on a track trajectory
    std::vector<TVector3> trkPos;
    std::vector<TVector3> trkDir;

  private:
    int fSPTAlg; //0: Use TrackTrajectoryAlg
                 //1: Use only Track3DReco alg

    bool fTrajOnly; // if true, only return trajectory points, if false, return
                    // a 3D point for every hit

    // use TrackTrajectoryAlg to get trajectory points
    void TrackTrajectory(detinfo::DetectorPropertiesData const& detProp,
                         std::vector<art::Ptr<recob::Hit>>& fHits);

    // use algorithm in Track3DReco
    void Track3D(detinfo::DetectorClocksData const& clockData,
                 detinfo::DetectorPropertiesData const& detProp,
                 std::vector<art::Ptr<recob::Hit>>& fHits);
    double ftmatch; ///< tolerance for time matching (in ticks)
    double fsmatch; ///< tolerance for distance matching (in cm)

    // create one 3D point for each hit using trajectory points
    void MakeSPT(detinfo::DetectorClocksData const& clockData,
                 detinfo::DetectorPropertiesData const& detProp,
                 std::vector<art::Ptr<recob::Hit>>& fHits);

    // track trajectory for a track under construction
    TrackTrajectoryAlg fTrackTrajectoryAlg;

    //projection of trajectory points on wire planes
    std::vector<std::vector<std::vector<std::vector<double>>>> vw;
    std::vector<std::vector<std::vector<std::vector<double>>>> vt;
    std::vector<std::vector<std::vector<std::vector<unsigned int>>>> vtraj;

    art::ServiceHandle<geo::Geometry const> geom;
    const detinfo::LArProperties* larprop;

  }; //class CosmicTrackerAlg
} // namespace trkf

#endif //ifndef COSMICTRACKERALG_H
