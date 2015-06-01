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

#include "RecoBase/Hit.h"
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/TrackTrajectoryAlg.h"

#include <vector>

namespace trkf
{
  class CosmicTrackerAlg {

  public:

    CosmicTrackerAlg(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset);

    void SPTReco(std::vector<art::Ptr<recob::Hit> >&fHits);

    std::vector<unsigned int> usehit;
    std::vector<double> spx;
    std::vector<double> spy;
    std::vector<double> spz;

    //position and direction of each point on a track trajectory
    std::vector<TVector3> trkPos;
    std::vector<TVector3> trkDir;

  private:
    
    void MakeSPT(std::vector<art::Ptr<recob::Hit> >&fHits);

    // track trajectory for a track under construction

    TrackTrajectoryAlg fTrackTrajectoryAlg;

    //trajectory position and direction returned by TrackTrajectoryAlg
    std::vector<TVector3> trajPos;
    std::vector<TVector3> trajDir;

    //projection of trajectory points on wire planes
    std::vector<std::vector<std::vector<std::vector<double>>>> vw;
    std::vector<std::vector<std::vector<std::vector<double>>>> vt;
    std::vector<std::vector<std::vector<std::vector<unsigned int>>>> vtraj;


    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;


  }; //class CosmicTrackerAlg
}// namespace trkf

#endif //ifndef COSMICTRACKERALG_H
