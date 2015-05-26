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

    std::vector<TVector3> trkPos;
    std::vector<TVector3> trkDir;

  private:
    
    void GetClusterSegment(std::vector<art::Ptr<recob::Hit> >&fHits);

    void GetSPs(std::vector<art::Ptr<recob::Hit> >&fHits, 
		unsigned int ipl0, unsigned int ipl1);

    unsigned int fMinNumCluHits;  ///< minimal number of hits in the cluster
    unsigned int fSegmentSize;    ///< number of hits in cluster segments

    //index of hits on each plane
    std::vector<std::vector<unsigned int> > hitindex;

    //wire and time ranges for cluster segments
    std::vector<std::vector<double> > w0;
    std::vector<std::vector<double> > w1;
    std::vector<std::vector<double> > t0;
    std::vector<std::vector<double> > t1;

    // track trajectory for a track under construction

    TrackTrajectoryAlg fTrackTrajectoryAlg;

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;


  }; //class CosmicTrackerAlg
}// namespace trkf

#endif //ifndef COSMICTRACKERALG_H
