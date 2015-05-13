//////////////////////////////////////////////////////////////////////
///
/// ClusterMatchTQ class
///
/// tjyang@fnal.gov
///
/// Algorithm for matching clusters between different views
/// based on time and charge information
///
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTERMATCHTQ_H
#define CLUSTERMATCHTQ_H
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h" 
#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"

#include <vector>

namespace cluster
{
  class ClusterMatchTQ {
  public:
    
    ClusterMatchTQ(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& p);

    void ClusterMatch(std::vector<art::Ptr<recob::Cluster> > clusterlist,
		      art::FindManyP<recob::Hit> fm);

    std::vector<std::vector<unsigned int> > matchedclusters;

  private:

    double fKSCut;
    bool   fEnableU;
    bool   fEnableV;
    bool   fEnableZ;

  }; // class ClusterMatchTQ
} // namespace cluster

#endif //ifndef CLUSTERMATCHTQ_H 
