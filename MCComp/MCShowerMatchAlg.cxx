#ifndef RECOTOOL_MCSHOWERMATCHALG_CXX
#define RECOTOOL_MCSHOWERMATCHALG_CXX

#include "MCShowerMatchAlg.h"

namespace btutil {

  MCShowerMatchAlg::MCShowerMatchAlg() {

    _name = "MCShowerMatchAlg";
    _view_to_plane.clear();
    //fBTAlgo.SetMaxEnergyCut(2000);
    //fBTAlgo.SetMinEnergyCut(20);
  }

  bool MCShowerMatchAlg::BuildMap(const std::vector< unsigned int>      &g4_trackid_v,
				  const std::vector< ::sim::SimChannel> &simch_v,
				  const std::vector< ::recob::Hit>      &hit_v,
				  const std::vector<std::vector<size_t> > &cluster_hit_association
				  )
  {
    fBTAlgo.Reset(g4_trackid_v,simch_v);

    //auto geo = ::larutil::Geometry::GetME();
    art::ServiceHandle<geo::Geometry> geo;

    //
    // Perform back-tracking
    //
    // (1) Load cluster/hit data product
    // (2) Loop over all clusters and find charge fraction to each MCShower
    // (3) For each MCShower, find a cluster per plane with the highest charge
    
    // Loop over clusters & get charge fraction per MCShower
    _summed_mcq.clear();
    _cluster_mcq_v.clear();
    _cluster_plane_id.clear();

    _summed_mcq.resize(g4_trackid_v.size()+1,std::vector<double>(geo->Nplanes(),0));
    _cluster_mcq_v.reserve(ev_cluster->size());
    _cluster_plane_id.reserve(ev_cluster->size());

    for(auto const& h_index_v : hit_ass) {

      size_t plane = geo->Nplanes();

      // Create hit list
      std::vector<WireRange_t> wr_v;
      wr_v.reserve(h_index_v.size());

      for(auto const& index : h_index_v) {
	
	WireRange_t wr;
	auto const& h = hit_v.at(index);
	wr.ch    = h.Channel();
	wr.start = h.StartTime();
	wr.end   = h.EndTime();
	wr_v.push_back(wr);
	if(plane==geo->Nplanes()) plane = h.WireID().Plane;
      }

      _cluster_plane_id.push_back(plane);

      auto mcq_v = fBTAlgo.MCQ(hit_v);

      for(size_t i=0; i<mcq_v.size(); ++i)
	_summed_mcq[i][plane] += mcq_v[i];

      _cluster_mcq_v.push_back(mcq_v);

    }

    //
    // Find the best matched pair (and its charge) per MCShower
    //
    _bmatch_id.clear();
    _bmatch_id.resize(g4_trackid_v.size(),std::vector<int>(geo->Nplanes(),-1));

    std::vector<std::vector<double> > bmatch_mcq (g4_trackid_v.size(),std::vector<double>(geo->Nplanes(),0));

    for(size_t shower_index=0; shower_index < mcshower_id.size(); ++shower_index) {

      for(size_t cluster_index=0; cluster_index < ev_cluster->size(); ++cluster_index) {

	if((*(_cluster_mcq_v.at(cluster_index).rbegin())) < 0) continue;

	auto plane = _cluster_plane_id.at(cluster_index);
	
	double q = _cluster_mcq_v.at(cluster_index).at(shower_index);
	
	if( bmatch_mcq.at(shower_index).at(plane) < q ) {

	  bmatch_mcq.at(shower_index).at(plane) = q;
	  _bmatch_id.at(shower_index).at(plane) = cluster_index;

	}
      }
    }

    _mcshower_index_v.resize(ev_mcshower->size(),-1);
    for(size_t i=0; i<fBTAlgo.UniqueShowerID().size(); ++i)

      _mcshower_index_v.at(fBTAlgo.UniqueShowerID()[i]) = i;

    fClusterProducer = cluster_producer;
    return true;
  }

  size_t MCShowerMatchAlg::__IndexConversion__(const size_t mcshower_index) const
  {
    if(mcshower_index >= _mcshower_index_v.size())
      throw MCShowerBTException(Form("MCShower index %zu exceeds # of MCShower data in the file (%zu)!",
				     mcshower_index,
				     _mcshower_index_v.size()));

    if(_mcshower_index_v[mcshower_index] < 0)
      throw MCShowerBTException(Form("MCShower index %zu not among the relevant MCShower list!",
				     mcshower_index)
				);
    return _mcshower_index_v[mcshower_index];

  }
  
  double MCShowerMatchAlg::ClusterCorrectness(const size_t cluster_index,
					      const size_t mcshower_index) const
  {    

    if(!_bmatch_id.size()) 
      throw MCShowerBTException("Preparation not done yet!");

    if( cluster_index >= _cluster_mcq_v.size()) 
      throw MCShowerBTException(Form("Input cluster index (%zu) out of range (%zu)!",
				     cluster_index,
				     _cluster_mcq_v.size())
				);

    auto shower_index = __IndexConversion__(mcshower_index);

    auto plane = _cluster_plane_id.at(cluster_index);

    auto best_cluster_index = _bmatch_id.at(shower_index).at(plane);
    
    return _cluster_mcq_v.at(cluster_index).at(shower_index) / _cluster_mcq_v.at(best_cluster_index).at(shower_index);

  }

  std::pair<size_t,double> MCShowerMatchAlg::ShowerCorrectness(const std::vector<unsigned int> cluster_indices) const
  {    
    
    if(!cluster_indices.size()) throw MCShowerBTException("Input cluster indices empty!");

    auto& mcshower_id = fBTAlgo.UniqueShowerID();

    // Compute efficiency per MCShower
    std::vector<double> match_eff(mcshower_id.size(),1);
    
    for(auto const& cluster_index : cluster_indices) {

      for(size_t shower_index=0; shower_index < mcshower_id.size(); ++shower_index)

	match_eff.at(shower_index) *= ClusterCorrectness(cluster_index, mcshower_id[shower_index]);
    }

    std::pair<size_t,double> result(0,-1);
    
    // Find the best qratio
    for(size_t shower_index=0; shower_index < mcshower_id.size(); ++shower_index) {
      
      if(match_eff.at(shower_index) < result.second) continue;

      result.second = match_eff.at(shower_index);
      
      result.first = mcshower_id.at(shower_index);
      
    }
    return result;
  }

  std::pair<double,double>  MCShowerMatchAlg::ClusterEP(const size_t cluster_index,
							const size_t mcshower_index) const
  {

    if(cluster_index >= _cluster_plane_id.size())
      throw MCShowerBTException(Form("Cluster index %zu exceeds the # of loaded clusters (by \"%s\")!",
				     cluster_index,
				     fClusterProducer.c_str())
				);

    auto shower_index = __IndexConversion__(mcshower_index);

    // efficiency = this cluster's mcq for specified mcshower / total mcq by this mcshower in all clusters
    // purity = this cluster's mcq for this mcshower / total mcq from all mcshower in this cluster

    std::pair<double,double> result;

    auto plane = _cluster_plane_id[cluster_index];

    result.first  = _cluster_mcq_v[cluster_index][shower_index] / _summed_mcq[shower_index][plane];

    double cluster_mcq_total = 0;
    for(auto const& q : _cluster_mcq_v[cluster_index]) cluster_mcq_total += q;

    result.second =  _cluster_mcq_v[cluster_index][shower_index] / cluster_mcq_total;

    return result;
  }

  const std::vector<int>& MCShowerMatchAlg::BestClusters(const size_t mcshower_index) const
  {
    auto shower_index = __IndexConversion__(mcshower_index);

    return _bmatch_id[shower_index];

  }

  std::pair<double,double> MCShowerMatchAlg::BestClusterEP(const size_t mcshower_index,
							   const size_t plane_id) const
  {

    auto c_index_v = BestClusters(mcshower_index);

    if(c_index_v.size() <= plane_id)
      throw MCShowerBTException(Form("Plane ID %zu exceeds # of planes recorded in data (%zu)...",
				     plane_id,
				     c_index_v.size())
				);

    std::pair<double,double> result(0,0);

    if(c_index_v[plane_id]<0) return result;

    return ClusterEP(c_index_v[plane_id],mcshower_index);

  }

}
#endif
