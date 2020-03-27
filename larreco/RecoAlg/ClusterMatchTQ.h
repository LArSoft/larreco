//////////////////////////////////////////////////////////////////////
///
/// ClusterMatchTQ class
///
/// tjyang@fnal.gov
///
/// Algorithm for matching clusters between different views
/// based on time and charge information
///
/// Input: a list of clusters and all hits associated with clusters
/// Output: a vector of index vectors. Each group of indices represent
///         a particle candidate
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTERMATCHTQ_H
#define CLUSTERMATCHTQ_H

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/fwd.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

#include <vector>

namespace cluster {
  class ClusterMatchTQ {
  public:
    ClusterMatchTQ(fhicl::ParameterSet const& pset);

    std::vector<std::vector<unsigned int>> MatchedClusters(
      const detinfo::DetectorClocksData& clockdata,
      const detinfo::DetectorPropertiesData& detProp,
      const std::vector<art::Ptr<recob::Cluster>>& clusterlist,
      const art::FindManyP<recob::Hit>& fm) const;

  private:
    double fKSCut;
    bool fEnableU;
    bool fEnableV;
    bool fEnableZ;

  }; // class ClusterMatchTQ
} // namespace cluster

#endif //ifndef CLUSTERMATCHTQ_H
