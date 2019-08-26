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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
namespace fhicl { class ParameterSet; }

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcore/Geometry/Geometry.h"
namespace geo { struct WireID; }

//tmp
#include "lardataobj/RecoBase/Hit.h"

#include "TVector2.h"
class TTree;

#include <vector>
#include <map>

namespace cluster {
  class MergeClusterAlg;
}

class cluster::MergeClusterAlg {
public:

  MergeClusterAlg(fhicl::ParameterSet const& pset);

  void     FindClusterEndPoints(art::PtrVector<recob::Hit> const& cluster, TVector2 const& centre, TVector2 const& direction, TVector2& start, TVector2& end) const;
  double   FindClusterOverlap(TVector2 const& direction, TVector2 const& centre, TVector2 const& start1, TVector2 const& end1, TVector2 const& start2, TVector2 const& end2) const;
  double   FindCrossingDistance(TVector2 const &direction1, TVector2 const &centre1, TVector2 const&direction2, TVector2 const &centre2) const;
  double   FindMinSeparation(art::PtrVector<recob::Hit> const &cluster1, art::PtrVector<recob::Hit> const &cluster2) const;
  double   FindProjectedWidth(TVector2 const& centre1, TVector2 const& start1, TVector2 const& end1, TVector2 const& centre2, TVector2 const& start2, TVector2 const& end2) const;
  double   GlobalWire(geo::WireID const& wireID) const;
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit) const;
  int      MergeClusters(std::vector<art::PtrVector<recob::Hit> > const &planeClusters, std::vector<art::PtrVector<recob::Hit> > &clusters) const;
  void     reconfigure(fhicl::ParameterSet const& p);
  bool     PassCuts(double const& angle, double const& crossingDistance, double const& projectedWidth, double const& separation, double const& overlap, double const& longLength) const;

private:

  // Merging parameters
  unsigned int fMinMergeClusterSize; // Minimum size of a cluster for it to be considered for merging
  double fMaxMergeSeparation;        // Maximum separation of clusters for merging
  double fProjWidthThreshold;        // Maximum projected width (width of a tube parallel to the line connecting centres of clusters which just encompasses the clusters) for merging

  // Create geometry and detector property handle
  art::ServiceHandle<geo::Geometry const> fGeom;
  art::ServiceHandle<art::TFileService const> tfs;

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
