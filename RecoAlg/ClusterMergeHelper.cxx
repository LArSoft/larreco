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
    if(fMgr.GetClusters().size()) 

      throw cet::exception(__PRETTY_FUNCTION__) << "\033[93m"
						<< "Merged cluster set not empty... Called Process() twice?"
						<< "\033[00m"
						<< std::endl;

    // Process
    fMgr.Process();

    // Now create output clusters
    CBookKeeper res = fMgr.GetBookKeeper();
    
    std::vector<std::vector<unsigned short> > out_clusters = res.GetResult();

    fOutputClusters.clear();
    
    fOutputClusters.reserve(out_clusters.size());

    for(auto const &cluster_index_v : out_clusters) {
      
      std::vector<art::Ptr<recob::Hit> > out_cluster;

      for(auto const& cluster_index : cluster_index_v) {

	out_cluster.reserve(out_cluster.size() + fInputClusters.at(cluster_index).size());

	for(auto const& hit_ptr : fInputClusters.at(cluster_index))

	  out_cluster.push_back(hit_ptr);

      }

      fOutputClusters.push_back(out_cluster);

    }

  }

  //######################################################################################################
  const std::vector<std::vector<art::Ptr<recob::Hit> > >& ClusterMergeHelper::GetMergedClusterHits() const
  //######################################################################################################
  { 
    if(!fOutputClusters.size()) 

      throw cet::exception(__FUNCTION__) 
	<< "\033[93m"
	<< "You must call Process() before calling " << __FUNCTION__ << " to retrieve result."
	<< "\033[00m"
	<< std::endl;

    return fOutputClusters;
  }
  
  //#####################################################################################
  const std::vector<cluster::ClusterParamsAlg>& ClusterMergeHelper::GetMergedCPAN() const
  //#####################################################################################
  {
    if(!fOutputClusters.size()) 

      throw cet::exception(__FUNCTION__) 
	<< "\033[93m"
	<< "You must call Process() before calling " << __FUNCTION__ << " to retrieve result."
	<< "\033[00m"
	<< std::endl;

    return fMgr.GetClusters();
  }

  void ClusterMergeHelper::AppendResult(art::EDProducer &ed,
					art::Event      &ev,
					std::vector<recob::Cluster>           &out_clusters,
					art::Assns<recob::Cluster,recob::Hit> &assns) const
  {

    if(!fOutputClusters.size()) 

      throw cet::exception(__FUNCTION__) 
	<< "\033[93m"
	<< "You must call Process() before calling " << __FUNCTION__ << " to retrieve result."
	<< "\033[00m"
	<< std::endl;

    art::ServiceHandle<geo::Geometry> geo;

    // Store output
    for(size_t out_index=0; out_index < GetMergedCPAN().size(); ++out_index) {
      
      // To save typing let's just retrieve const cluster_params instance
      const cluster_params &res = GetMergedCPAN().at(out_index).GetParams();
      
      // View_t needed but not a part of cluster_params, so retrieve it here
      geo::View_t view_id = geo->Plane(GetMergedCPAN().at(out_index).Plane()).View();
      
      // Push back a new cluster data product with parameters copied from cluster_params
      out_clusters.push_back( recob::Cluster( res.start_point.w / fGeoU.WireToCm(), 0,  // start wire & error
					      res.start_point.t / fGeoU.TimeToCm(), 0,  // start time & error
					      res.end_point.w   / fGeoU.WireToCm(), 0,  // end   wire & error
					      res.end_point.t   / fGeoU.TimeToCm(), 0,  // end   time & error
					      res.cluster_angle_2d,                 0,  // dT/dW (slope)
					      0,                                    0,  // dQ/dW (what is that?)
					      res.sum_charge,                           // charge sum
					      view_id,                                  // geo::View_t
					      out_clusters.size()                       // Cluster ID
					      )
			      );
      
      util::CreateAssn(ed,
		       ev, 
		       out_clusters,
		       GetMergedClusterHits().at(out_index), 
		       assns
		       );
    }
    
  }

} // namespace cluster

#endif 
