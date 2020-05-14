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
namespace fhicl {
  class ParameterSet;
}

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include <vector>

namespace cluster {
  class ClusterMatchTQ {
  public:
    ClusterMatchTQ(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& p);

    void ClusterMatch(const std::vector<art::Ptr<recob::Cluster>>& clusterlist,
                      const art::FindManyP<recob::Hit>& fm);

    std::vector<std::vector<unsigned int>> matchedclusters;

  private:
    double fKSCut;
    bool fEnableU;
    bool fEnableV;
    bool fEnableZ;

  }; // class ClusterMatchTQ
} // namespace cluster

#endif //ifndef CLUSTERMATCHTQ_H
