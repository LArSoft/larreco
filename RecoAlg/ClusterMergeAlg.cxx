////////////////////////////////////////////////////////////////////////
//
//  ClusterMergeAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef CLUSTERMERGEALG_CC
#define CLUSTERMERGEALG_CC

#include "ClusterMergeAlg.h"

namespace cluster{

  //##################################################################
  ClusterMergeAlg::ClusterMergeAlg(fhicl::ParameterSet const& pset)
  //##################################################################
  {    
          
    _verbose = pset.get<bool> ("VerboseMode");
    
    SetAngleCut(pset.get<double>("MaxAngleDiff",180.));

    SetSquaredDistanceCut(pset.get<double>("MaxSquaredDist",1e9));

    _min_distance_unit = -1;

    _det_params_prepared = false;

    _merge_tree = 0;

    ClearEventInfo();
  }

  //##################################################################
  void ClusterMergeAlg::ReportConfig() const 
  //##################################################################
  {

    std::ostringstream msg;
    msg
      << std::endl
      << " ClusterMergeAlg Configuration:              " << std::endl
      << "---------------------------------------------" << std::endl;
    msg 
      << " Verbose Mode ... " << (_verbose ? "enabled!" : "disabled!") << std::endl
      << " Wire => Cm Conversion: " << _wire_2_cm << std::endl
      << " Time => Cm Conversion: " << _time_2_cm << std::endl
      << std::endl
      << " Squared-Distance Cut: " << _max_2D_dist2 << std::endl
      << " Angle Difference Cut: " << _max_allowed_2D_angle_diff << std::endl
      << std::endl
      << "---------------------------------------------" << std::endl;

    mf::LogWarning(__FUNCTION__)<<msg.str();

  }

  //############################################################################################  
  void ClusterMergeAlg::AppendClusterInfo(const recob::Cluster &cl,
					  const std::vector<art::Ptr<recob::Hit> > &in_hit_v) 
  //############################################################################################  
  {

    PrepareDetParams();
    cluster_merge_info ci;
    ci.cluster_index = cl.ID();
    ci.view = cl.View();
    ci.start_wire = cl.StartPos()[0];
    ci.start_time = cl.StartPos()[1];
    ci.end_wire   = cl.EndPos()[0];
    ci.end_time   = cl.EndPos()[1];
    ci.angle      = cl.dTdW();

    AppendHitInfo(ci,in_hit_v);
    
    if(ci.view == geo::kU) _u_clusters.push_back(ci);
    else if(ci.view == geo::kV) _v_clusters.push_back(ci);
    else if(ci.view == geo::kW) _w_clusters.push_back(ci);
    else mf::LogError(__FUNCTION__)<<Form("Invalid plane ID: %d",ci.view);

    if(_cluster_merged_index.size() <= (size_t)(ci.cluster_index)) {

      _cluster_merged_index.resize((size_t)(ci.cluster_index+1),-1);
      
      mf::LogInfo(__FUNCTION__)<<Form("Extending the input cluster index ... length = %zu",_cluster_merged_index.size());
    }

  }

  //############################################################################################    
  void ClusterMergeAlg::AppendClusterInfo(const art::Ptr<recob::Cluster> cl,
					  const std::vector<art::Ptr<recob::Hit> > &in_hit_v) 
  //############################################################################################  
  {

    PrepareDetParams();
    cluster_merge_info ci;
    ci.cluster_index = cl->ID();
    ci.view = cl->View();
    ci.start_wire = cl->StartPos()[0];
    ci.start_time = cl->StartPos()[1];
    ci.end_wire   = cl->EndPos()[0];
    ci.end_time   = cl->EndPos()[1];
    ci.angle      = cl->dTdW();

    AppendHitInfo(ci,in_hit_v);
    
    if(ci.view == geo::kU) _u_clusters.push_back(ci);
    else if(ci.view == geo::kV) _v_clusters.push_back(ci);
    else if(ci.view == geo::kW) _w_clusters.push_back(ci);
    else mf::LogError(__FUNCTION__)<<Form("Invalid plane ID: %d",ci.view);

    if(_cluster_merged_index.size() <= (size_t)(ci.cluster_index)) {

      _cluster_merged_index.resize((size_t)(ci.cluster_index+1),-1);
      
      mf::LogInfo(__FUNCTION__)<<Form("Extending the input cluster index ... length = %zu",_cluster_merged_index.size());
    }

  }

  //####################################################################################  
  void ClusterMergeAlg::AppendHitInfo(cluster_merge_info &ci, 
				      const std::vector<art::Ptr<recob::Hit> > &in_hit_v) 
  //####################################################################################  
  {}

  //##################################################################
  void ClusterMergeAlg::PrepareTTree()
  //##################################################################
  {
    if(!_merge_tree){
      art::ServiceHandle<art::TFileService> fileService;
      _merge_tree = fileService->make<TTree>("merge_tree","");
    }
  }

  //##################################################################
  void ClusterMergeAlg::PrepareDetParams() 
  //##################################################################
  {
    if(!_det_params_prepared){
      // Total number of planes
      art::ServiceHandle<geo::Geometry> geo_h;
      art::ServiceHandle<util::DetectorProperties> det_h;
      art::ServiceHandle<util::LArProperties> lar_h;

      _wire_2_cm = geo_h->WirePitch(0,1,0);
      _time_2_cm = (det_h->SamplingRate()/1.e3) * (lar_h->DriftVelocity(lar_h->Efield(), lar_h->Temperature()));

      if(_min_distance_unit < 0) {

	double smaller_factor = (_wire_2_cm < _time_2_cm) ? _wire_2_cm : _time_2_cm;
	
	_min_distance_unit = pow(smaller_factor,2);

      }

      ReportConfig();
      _det_params_prepared = true;
    }
  }

  //##################################################################
  void ClusterMergeAlg::ClearInputInfo()
  //##################################################################
  {

    _u_clusters.clear();
    _v_clusters.clear();
    _w_clusters.clear();
    _cluster_merged_index.clear();

  }

  //##################################################################
  void ClusterMergeAlg::ClearOutputInfo()
  //##################################################################
  {

    _cluster_sets_v.clear();

  }

  //##################################################################
  void ClusterMergeAlg::ClearTTreeInfo()
  //##################################################################
  {}

  //##################################################################
  void ClusterMergeAlg::ClearEventInfo() 
  //##################################################################
  {

    // Clear input event data holders
    ClearInputInfo();

    // Clear result data holders
    ClearOutputInfo();

    /// Clear TTree variables
    ClearTTreeInfo();

  }

  //##################################################################
  void ClusterMergeAlg::ProcessMergeAlg() 
  //##################################################################
  {

    // Clear all algorithm's output
    ClearOutputInfo();
    
    //Loop over all possible combinations of clusters in u-view
    for(int iclus = 0; iclus < (int)_u_clusters.size(); ++iclus){

      for(int jclus = iclus+1; jclus < (int)_u_clusters.size(); ++jclus){

	BuildClusterSets(_u_clusters.at(iclus),_u_clusters.at(jclus));

      }//end loop over jclus _u_clusters

    }//end loop over iclus _u_clusters

    //Loop over all possible combinations of clusters in v-view
    for(int iclus = 0; iclus < (int)_v_clusters.size(); ++iclus){

      for(int jclus = iclus+1; jclus < (int)_v_clusters.size(); ++jclus){

	BuildClusterSets(_v_clusters.at(iclus),_v_clusters.at(jclus));

      }//end loop over jclus _v_clusters

    }//end loop over iclus _v_clusters

    //Loop over all possible combinations of clusters in w-view
    for(int iclus = 0; iclus < (int)_w_clusters.size(); ++iclus){

      for(int jclus = iclus+1; jclus < (int)_w_clusters.size(); ++jclus){

	BuildClusterSets(_w_clusters.at(iclus),_w_clusters.at(jclus));

      }//end loop over jclus _w_clusters

    }//end loop over iclus _w_clusters
    
    //make sure all of the un-mergable clusters are in _cluster_sets_v, individually
    FinalizeClusterSets();

    if(_merge_tree) _merge_tree->Fill();

  }

  //######################################################################
  bool ClusterMergeAlg::CompareClusters(const cluster_merge_info &clusA,
					 const cluster_merge_info &clusB)
  //######################################################################
  {
    
    if(_verbose) {
      std::cout
	<< std::endl << "Printing out two input cluster information..." << std::endl;
      PrintClusterVars(clusA);
      PrintClusterVars(clusB);
    }

    bool merge_clusters = true;

    merge_clusters = merge_clusters && Angle2DCompatibility(clusA, clusB);

    merge_clusters = merge_clusters &&  ShortestDistanceCompatibility(clusA, clusB);

    return merge_clusters;

  }

  //######################################################################################
  bool ClusterMergeAlg::Angle2DCompatibility(const cluster_merge_info &cluster_a, 
					      const cluster_merge_info &cluster_b) const 
  //######################################################################################
  {
        
    double angle1 = cluster_a.angle * _time_2_cm / _wire_2_cm;
    double angle2 = cluster_b.angle * _time_2_cm / _wire_2_cm;

    bool compatible = ( abs(angle1-angle2)     < _max_allowed_2D_angle_diff ||
			abs(angle1-angle2-180) < _max_allowed_2D_angle_diff ||
			abs(angle1-angle2+180) < _max_allowed_2D_angle_diff   );

    if(_verbose) {

      if(compatible) std::cout << " Compatible in angle." << std::endl;
      else std::cout << " NOT compatible in angle." << std::endl;

    }

    return compatible;

  }//end Angle2DCompatibility

  //################################################################################################
  bool ClusterMergeAlg::ShortestDistanceCompatibility(const cluster_merge_info &clus_info_A,
						       const cluster_merge_info &clus_info_B) const 
  //################################################################################################
  {

    double w_start1 = clus_info_A.start_wire * _wire_2_cm;
    double t_start1 = clus_info_A.start_time * _time_2_cm;
    double w_end1   = clus_info_A.end_wire   * _wire_2_cm;
    double t_end1   = clus_info_A.end_time   * _time_2_cm;

    double w_start2 = clus_info_B.start_wire * _wire_2_cm;
    double t_start2 = clus_info_B.start_time * _time_2_cm;
    double w_end2   = clus_info_B.end_wire   * _wire_2_cm;
    double t_end2   = clus_info_B.end_time   * _time_2_cm;
    
    //First, pretend the first cluster is a 2D line segment, from its start point to end point
    //Find the shortest distance between start point of the second cluster to this line segment.
    //Repeat for end point of second cluster to this line segment.
    //Then, pretend second cluster is a 2D line segment, from its start point to end point.
    //Find the shortest distance between start point of the first cluster to this line segment.
    //Repeat for end point of first cluster to this line segment.
    //If the shortest of these four distances is less than the cutoff, 
    //return true... the clusters are merge-compatible. else, return false.
    
    // Step 1: inspect (w_start1, t_start1) vs. line (w_start2, t_start2) => (w_end2, t_end2)
    double shortest_distance2 = ShortestDistanceSquared(w_start1, t_start1,
							w_start2, t_start2, w_end2, t_end2);
    
    // Step 2: inspect (w_end1, t_end1) vs. line (w_start2, t_start2) => (w_end2, t_end2)
    double shortest_distance2_tmp = ShortestDistanceSquared(w_end1, t_end1,
							    w_start2, t_start2, w_end2, t_end2);

    shortest_distance2 = (shortest_distance2_tmp < shortest_distance2) ? shortest_distance2_tmp : shortest_distance2;

    // Step 3: inspect (w_start2, t_start2) vs. line (w_start1, t_start1) => (w_end1, t_end1)
    shortest_distance2_tmp = ShortestDistanceSquared(w_start2, t_start2,
						     w_start1, t_start1, w_end1, t_end1);

    shortest_distance2 = (shortest_distance2_tmp < shortest_distance2) ? shortest_distance2_tmp : shortest_distance2;

    // Step 4: inspect (w_end2, t_end2) vs. line (w_start1, t_start1) => (w_end1, t_end1)
    shortest_distance2_tmp = ShortestDistanceSquared(w_end2, t_end2,
						     w_start1, t_start1, w_end1, t_end1);
    
    shortest_distance2 = (shortest_distance2_tmp < shortest_distance2) ? shortest_distance2_tmp : shortest_distance2;

    bool compatible = shortest_distance2 < _max_2D_dist2;

    if(_verbose) {

      if(compatible) std::cout << Form(" Compatible in distance (%g)",shortest_distance2) << std::endl;
      else std::cout << Form(" NOT compatible in distance (%g).",shortest_distance2) << std::endl;

    }

    return compatible;

  }//end startend2dcompatibility

  //####################################################################################
  double ClusterMergeAlg::ShortestDistanceSquared(double point_x, double point_y, 
						   double start_x, double start_y,
						   double end_x,   double end_y  ) const 
  //####################################################################################
  {
    
    //This code finds the shortest distance between a point and a line segment.    
    //code based off sample from 
    //http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    //note to self: rewrite this with TVector2 and compare time differences... 
    //TVector2 code might be more understandable

    double distance_squared = -1;

    // Line segment: from ("V") = (start_x, start_y) to ("W")=(end_x, end_y)
    double length_squared = pow((end_x - start_x), 2) + pow((end_y - start_y), 2);

    // Treat the case start & end point overlaps
    if(length_squared < _min_distance_unit) {

      std::ostringstream msg;
      msg 
	<< std::endl
	<< Form(" Provided very short line segment: (%g,%g) => (%g,%g)",start_x,start_y,end_x,end_y) << std::endl
	<< " Likely this means one of two clusters have start & end point identical." << std::endl
	<< " Check the cluster output!" << std::endl
	<< std::endl
	<< Form(" At this time, the algorithm uses a point (%g,%g)",start_x,start_y) << std::endl
	<< " to represent this cluster's location." << std::endl
	<< std::endl;
      
      mf::LogWarning(__FUNCTION__)<<msg.str();

      return (pow((point_x - start_x),2) + pow((point_y - start_y),2));

    }

    //Find shortest distance between point ("P")=(point_x,point_y) to this line segment
    double t = ( (point_x - start_x)*(end_x - start_x) + (point_y - start_y)*(end_y - start_y) ) / length_squared;

    if(t<0.0) distance_squared = pow((point_x - start_x), 2) + pow((point_y - start_y), 2);

    else if (t>1.0) distance_squared = pow((point_x - end_x), 2) + pow(point_y - end_y, 2);

    else distance_squared = pow((point_x - (start_x + t*(end_x - start_x))), 2) + pow((point_y - (start_y + t*(end_y - start_y))),2);

    if(_verbose)
      std::cout<<Form("(%g,%g) : (%g,%g) => (%g,%g) ... l2=%g, t=%g, d=%g",
		      point_x,point_y,start_x,start_y,end_x,end_y,length_squared,t,distance_squared)<<std::endl;
       
    return distance_squared;
  }  

  //#########################################################################
  void ClusterMergeAlg::PrintClusterVars(cluster_merge_info clus_info) const
  //#########################################################################
  {

    std::ostringstream msg;

    msg
      << std::endl
      << "***********CLUSTER INFO PARAMETERS***********" << std::endl
      << " ID: " << clus_info.cluster_index << std::endl
      << " start_wire : " << clus_info.start_wire << std::endl
      << " start_time : " << clus_info.start_time << std::endl
      << " end wire   : " << clus_info.end_wire << std::endl
      << " end_time   : " << clus_info.end_time << std::endl
      << " angle      : " << clus_info.angle << std::endl
      << std::endl;
    
    std::cout << msg.str() << std::endl;

  }//end PrintClusterVars function

  //#####################################################################################
  int ClusterMergeAlg::AppendToClusterSets(unsigned int cluster_index, int merged_index) 
  //#####################################################################################
  {

    // This function append the provided cluster_index into the _cluster_sets_v.
    // If merged_index is not provided (default=-1), then we assume this meands to append a new merged cluster set.
    // Note we also log the merged index in _cluster_merged_index vector to support fast search in isInClusterSets() function.

    // First, check if this cluster_index is already in the set or not.
    if(!(isInClusterSets(cluster_index)<0)) {

      mf::LogError(__FUNCTION__) << Form(" Cluster ID = %d already in the set!",cluster_index);

      return isInClusterSets(cluster_index);
      
    }

    if(merged_index < 0) {

      std::vector<unsigned int> tmp(1,cluster_index);

      _cluster_merged_index[cluster_index] = (unsigned int)(_cluster_sets_v.size());

      _cluster_sets_v.push_back(tmp);

    }
    else if(merged_index < (int)(_cluster_sets_v.size())) {

      _cluster_merged_index[cluster_index] = merged_index;

      _cluster_sets_v[merged_index].push_back(cluster_index);

    }
    else

      mf::LogError(__FUNCTION__)
	<< Form(" Requested to merge the cluster ID = %ud into the set = %d which does not exist!",cluster_index,merged_index);

    return isInClusterSets(cluster_index);
  }

  //######################################################################
  void ClusterMergeAlg::BuildClusterSets(const cluster_merge_info &clusA, 
					  const cluster_merge_info &clusB)
  //######################################################################
  {
    //this function builds the _cluster_sets_v vector, according to the format defined
    //in ClusterAlgo. here's an example of what _cluster_sets_v format is:
    //imagine an event with clusters 0, 1, 2, 3, 4, 5
    //pretend 0, 1 are found compatible for merging.
    //pretend 1, 2 are found compatible for merging.
    //pretend 4, 5 are found compatible for merging.
    //_cluster_sets_v is a vector<vector<int>>, and it should look like
    //_cluster_sets_v = ( (0, 1, 2), (4, 5), (3) )

    //if the two clusters are compatible for merging
    if( CompareClusters(clusA,clusB) ){
     
      //if this is the first compatibility found for this entire event
      if(_cluster_sets_v.empty()){
	//push both IDs in to the _cluster_sets_v vector, grouped together as one element
	AppendToClusterSets(clusA.cluster_index);
	AppendToClusterSets(clusB.cluster_index,0);
      }
      else {


	// Inspect if these clusters are already in the set or not.
	// Return value < 0 indicates the set does not hold the subject cluster.
	// Otherwise the function returns the index of a set that includes the cluster.
	int a_index = isInClusterSets(clusA.cluster_index);
	int b_index = isInClusterSets(clusB.cluster_index);

	// if only A index already somewhere in the _cluster_sets_v vector, and B index is not
	if(a_index >= 0 && b_index < 0) 

	  AppendToClusterSets(clusB.cluster_index, a_index);

	// else if only B index already somewhere in the _cluster_sets_v vector, and A index is not
	else if(a_index < 0 && b_index >= 0)

	  AppendToClusterSets(clusA.cluster_index, b_index);

	// else if neither of the two are already in the _cluster_sets_v vector, anywhere
	else if(a_index < 0 && b_index < 0) {
	  
	  a_index = AppendToClusterSets(clusA.cluster_index);
	  AppendToClusterSets(clusB.cluster_index,a_index);
	    
	}  

	// else if both indexes are already in the _cluster_sets_v vector but the indexes are not same
	else if(a_index != b_index)

	  mf::LogError(__FUNCTION__) << "LOGIC ERROR: Found two compatible clusters in different sets!";

	//at the end of all possible cluster matching permutations, if some are not in 
	//_cluster_sets_v yet, push them back individually as individual elements...
	//this is what FinalizeClusterSets function is for

      }//end if _cluster_sets_v is non-empty

    }//end if the two clusters are compatible for merging
    
  }//end BuildClusterSets

  //##########################################
  void ClusterMergeAlg::FinalizeClusterSets()
  //##########################################
  {

    //loop over all cluster ID's in the event... guess I have to do it view-by-view
    for(int iclus = 0; iclus < (int)_u_clusters.size(); ++iclus){    
      //if the cluster_id is not in _cluster_sets_v at all, push it back individually
      if(isInClusterSets(_u_clusters.at(iclus).cluster_index) < 0)

	AppendToClusterSets(_u_clusters.at(iclus).cluster_index);

    }

    for(int iclus = 0; iclus < (int)_v_clusters.size(); ++iclus){    
      //if the cluster_id is not in _cluster_sets_v at all, push it back individually
      if(isInClusterSets(_v_clusters.at(iclus).cluster_index) < 0)

	AppendToClusterSets(_v_clusters.at(iclus).cluster_index);

    }

    for(int iclus = 0; iclus < (int)_w_clusters.size(); ++iclus){    
      //if the cluster_id is not in _cluster_sets_v at all, push it back individually
      if(isInClusterSets(_w_clusters.at(iclus).cluster_index) < 0)

	AppendToClusterSets(_w_clusters.at(iclus).cluster_index);

    }

    
  } // end FinalizeClusterSets

  //##############################################################  
  int ClusterMergeAlg::isInClusterSets(unsigned int index) const
  //##############################################################  
  {
    //function to check if an index is in _cluster_sets_v anywhere,
    //and return where it is (-1 if it's not in _cluster_sets_v)

    // Check if the provided index does not exceed the book-keeping vector
    if((size_t)index < _cluster_merged_index.size())

      return _cluster_merged_index.at(index);

    else {

      mf::LogError(__FUNCTION__) << Form(" Cluster index %d not found among the input clusters!",index);

      return -2;

    }
    
  }//end isInClusterSets function

} // namespace cluster

#endif 
