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
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Ptr.h" 

#include "lardata/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/TrackTrajectoryAlg.h"

#include <vector>

namespace trkf
{
  class CosmicTrackerAlg {

  public:

    CosmicTrackerAlg(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset);

    void SPTReco(std::vector<art::Ptr<recob::Hit> >&fHits);

    //trajectory position and direction returned by TrackTrajectoryAlg
    std::vector<TVector3> trajPos;
    std::vector<TVector3> trajDir;
    std::vector<std::vector<art::Ptr<recob::Hit>>> trajHit;

    //position and direction of each point on a track trajectory
    std::vector<TVector3> trkPos;
    std::vector<TVector3> trkDir;

  private:

    int fSPTAlg;    //0: Use TrackTrajectoryAlg
                    //1: Use only Track3DReco alg

    bool fTrajOnly; //if true, only return trajectory points, if false, return a 3D point for every hit
    
    //use TrackTrajectoryAlg to get trajectory points
    void TrackTrajectory(std::vector<art::Ptr<recob::Hit> >&fHits);

    //use algorithm in Track3DReco
    void   Track3D(std::vector<art::Ptr<recob::Hit> >&fHits);
    double ftmatch;             ///< tolerance for time matching (in ticks) 
    double fsmatch;             ///< tolerance for distance matching (in cm)

    //create one 3D point for each hit using trajectory points
    void MakeSPT(std::vector<art::Ptr<recob::Hit> >&fHits);

    // track trajectory for a track under construction
    TrackTrajectoryAlg fTrackTrajectoryAlg;

    //projection of trajectory points on wire planes
    std::vector<std::vector<std::vector<std::vector<double>>>> vw;
    std::vector<std::vector<std::vector<std::vector<double>>>> vt;
    std::vector<std::vector<std::vector<std::vector<unsigned int>>>> vtraj;


    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::LArProperties* larprop;
    const detinfo::DetectorProperties* detprop;


  }; //class CosmicTrackerAlg
}// namespace trkf

#endif //ifndef COSMICTRACKERALG_H
