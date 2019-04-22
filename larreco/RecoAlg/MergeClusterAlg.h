////////////////////////////////////////////////////////////////////
// Merge Cluster algorithm
//
// Runs on the output of previous clustering algorithms to merge
// clusters together which lie on a straight line and are within
// some separation threshold.
// Runs recursively over all clusters, including new ones formed
// in the algorithm.
//
// M Wallbank (m.wallbank@sheffield.ac.uk), July 2015
////////////////////////////////////////////////////////////////////

#ifndef MergeCluster_h
#define MergeCluster_h

// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcore/Geometry/Geometry.h"

//tmp
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
//The backtracker is not actually used, so it is commented out.
//#include "larsim/MCCheater/BackTrackerService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TTree.h"
#include "TPrincipal.h"
#include "TVector2.h"

#include <vector>
#include <map>

namespace cluster {
  class MergeClusterAlg;
}

class cluster::MergeClusterAlg {
public:

  MergeClusterAlg(fhicl::ParameterSet const& pset);

  void     FindClusterEndPoints(art::PtrVector<recob::Hit> const& cluster, TVector2 const& centre, TVector2 const& direction, TVector2& start, TVector2& end);
  double   FindClusterOverlap(TVector2 const& direction, TVector2 const& centre, TVector2 const& start1, TVector2 const& end1, TVector2 const& start2, TVector2 const& end2);
  double   FindCrossingDistance(TVector2 const &direction1, TVector2 const &centre1, TVector2 const&direction2, TVector2 const &centre2);
  double   FindMinSeparation(art::PtrVector<recob::Hit> const &cluster1, art::PtrVector<recob::Hit> const &cluster2);
  double   FindProjectedWidth(TVector2 const& centre1, TVector2 const& start1, TVector2 const& end1, TVector2 const& centre2, TVector2 const& start2, TVector2 const& end2);
  double   GlobalWire(geo::WireID const& wireID);
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);
  int      MergeClusters(std::vector<art::PtrVector<recob::Hit> > const &planeClusters, std::vector<art::PtrVector<recob::Hit> > &clusters);
  void     reconfigure(fhicl::ParameterSet const& p);
  bool     PassCuts(double const& angle, double const& crossingDistance, double const& projectedWidth, double const& separation, double const& overlap, double const& longLength);

private:

  // Merging parameters
  unsigned int fMinMergeClusterSize; // Minimum size of a cluster for it to be considered for merging
  double fMaxMergeSeparation;        // Maximum separation of clusters for merging
  double fProjWidthThreshold;        // Maximum projected width (width of a tube parallel to the line connecting centres of clusters which just encompasses the clusters) for merging

  // Create geometry and detector property handle
  art::ServiceHandle<geo::Geometry const> fGeom;
//  const detinfo::DetectorProperties* fDetProp;
  art::ServiceHandle<art::TFileService const> tfs;
//  art::ServiceHandle<cheat::BackTracker const> backtracker;

  std::map<int,int> trueClusterMap;

  // Tree
  TTree *fTree;
  double fAngle;
  double fEigenvalue;
  int fCluster1Size;
  int fCluster2Size;
  double fLength1;
  double fLength2;
  double fSeparation;
  double fCrossingDistance;
  double fProjectedWidth;
  double fOverlap;
  bool fTrueMerge;
//  bool fMerge;

};

#endif
