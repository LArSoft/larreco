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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Hit.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/Geometry.h"

#include "TPrincipal.h"
#include "TVector2.h"

#include <vector>

namespace cluster {
  class MergeClusterAlg;
}

class cluster::MergeClusterAlg {
public:

  MergeClusterAlg(fhicl::ParameterSet const& pset);

  TVector2 ConvertWireDriftToCm(unsigned int wire, float drift, unsigned int plane, unsigned int tpc, unsigned int cryo);
  int MergeClusters(std::vector<art::PtrVector<recob::Hit> > *planeClusters, std::vector<art::PtrVector<recob::Hit> > &clusters, unsigned int plane, unsigned int tpc, unsigned int cryo);
  double MinSeparation(art::PtrVector<recob::Hit> &cluster1, art::PtrVector<recob::Hit> &cluster2, unsigned int plane, unsigned int tpc, unsigned int cryo);
  double MinSeparation(art::PtrVector<recob::Hit> &cluster1, art::PtrVector<recob::Hit> &cluster2);
  void reconfigure(fhicl::ParameterSet const& p);

private:

  unsigned int fMinMergeClusterSize;
  double fMaxMergeSeparation;
  double fMergingThreshold;

  // Create geometry and detector property handle
  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

};

#endif
