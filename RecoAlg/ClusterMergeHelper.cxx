////////////////////////////////////////////////////////////////////////
//
//  ClusterMergeHelper source
//
////////////////////////////////////////////////////////////////////////

#ifndef CLUSTERMERGEHELPER_CXX
#define CLUSTERMERGEHELPER_CXX

#include "ClusterMergeHelper.h"

namespace cluster{
  
  
  //####################################################################################################  
  void ClusterMergeHelper::SetClusters(const std::vector<std::vector<art::Ptr<recob::Hit> > > &clusters)
  //####################################################################################################  
  {
    fInputClusters.clear();
    fOutputClusters.clear();

    std::vector<std::vector<util::PxHit> > px_clusters(clusters.size(),std::vector<util::PxHit>());

    fInputClusters.resize(clusters.size(),std::vector<art::Ptr<recob::Hit> >());

    for(size_t cluster_index=0; cluster_index < clusters.size(); ++cluster_index) {

      px_clusters.at(cluster_index).resize(clusters.at(cluster_index).size(),util::PxHit());

      fInputClusters.at(cluster_index).resize(clusters.at(cluster_index).size());

      for(size_t hit_index=0; hit_index < clusters.at(cluster_index).size(); ++hit_index) {
	
	px_clusters.at(cluster_index).at(hit_index).plane  = clusters.at(cluster_index).at(hit_index)->WireID().Plane;
	px_clusters.at(cluster_index).at(hit_index).w      = clusters.at(cluster_index).at(hit_index)->WireID().Wire * fGeoU.WireToCm();
	px_clusters.at(cluster_index).at(hit_index).t      = clusters.at(cluster_index).at(hit_index)->PeakTime() * fGeoU.TimeToCm();
	px_clusters.at(cluster_index).at(hit_index).charge = clusters.at(cluster_index).at(hit_index)->Charge();

	fInputClusters.at(cluster_index).at(hit_index) = clusters.at(cluster_index).at(hit_index);

      }
      
    }

    SetClusters(px_clusters);

  }

  //##################################################################################################
  void ClusterMergeHelper::SetClusters(const art::Event &evt, const std::string &cluster_module_label)
  //##################################################################################################
  {

    art::Handle<std::vector<recob::Cluster> >  clusters_h;
    evt.getByLabel(cluster_module_label, clusters_h);

    if(!(clusters_h.isValid()))

      throw cet::exception(__FUNCTION__) 
	<< "\033[93m"
	<< " Failed to retrieve recob::Cluster with label: " 
	<< cluster_module_label.c_str()
	<< "\033[00m"
	<< std::endl;

    std::vector<std::vector<art::Ptr<recob::Hit> > > cluster_hits_v;

    cluster_hits_v.reserve(clusters_h->size());

    art::FindManyP<recob::Hit> hit_m(clusters_h, evt, cluster_module_label);

    for(size_t i=0; i<clusters_h->size(); ++i)

      cluster_hits_v.push_back(hit_m.at(i));

    SetClusters(cluster_hits_v);
    
  }

  //################################
  void ClusterMergeHelper::Process()
  //################################
  {

  }

  //#############################################################################################
  const std::vector<std::vector<art::Ptr<recob::Hit> > >& ClusterMergeHelper::GetClusters() const
  //#############################################################################################
  { 
    if(!fOutputClusters.size()) 

      throw cet::exception(__FUNCTION__) 
	<< "\033[93m"
	<< "You must call Process() before calling GetClusters() to retrieve result."
	<< "\033[00m"
	<< std::endl;

    return fOutputClusters;
  }

} // namespace cluster

#endif 
