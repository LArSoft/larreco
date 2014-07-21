#ifndef CMALGOSHORTESTDISTNONENDPOINT_CXX
#define CMALGOSHORTESTDISTNONENDPOINT_CXX

#include "CMAlgoShortestDistNonEndpoint.h"

namespace cluster {


  CMAlgoShortestDistNonEndpoint::CMAlgoShortestDistNonEndpoint() {

    //this just sets default values    
    SetVerbose(true);
    SetDebug(false);
    SetMinHits(0);
    
    //1e9 is huge; everything will be merged
    SetSquaredDistanceCut(1e9);

    _min_distance_unit = 1.e-3;
  } //end constructor

  bool CMAlgoShortestDistNonEndpoint::Bool(const ClusterParamsAlg &cluster1,
				const ClusterParamsAlg &cluster2)
  {
    
    //if number of hits not large enough skip
    if ( (_minHits > 0) and ((cluster1.GetParams().N_Hits < _minHits) or (cluster2.GetParams().N_Hits < _minHits)) ) {
      if (_debug) { std::cout << "Num of Hits below threshold..." << std::endl; }
      return false;
    }

    double w_start1 = cluster1.GetParams().start_point.w;
    double t_start1 = cluster1.GetParams().start_point.t;
    double w_end1   = cluster1.GetParams().end_point.w;
    double t_end1   = cluster1.GetParams().end_point.t;

    double w_start2 = cluster2.GetParams().start_point.w;
    double t_start2 = cluster2.GetParams().start_point.t;
    double w_end2   = cluster2.GetParams().end_point.w;
    double t_end2   = cluster2.GetParams().end_point.t;

    if (_debug){
      std::cout << "Start point Cluster 1: (" << w_start1 << ", " << t_start1 << ")"  << std::endl;
      std::cout << "End point Cluster 2: (" << w_end1 << ", " << t_end1 << ")"  << std::endl;
      std::cout << "Start point Cluster 1: (" << w_start2 << ", " << t_start2 << ")"  << std::endl;
      std::cout << "End point Cluster 2: (" << w_end2 << ", " << t_end2 << ")"  << std::endl;
    }
    
    //First, pretend the first cluster is a 2D line segment, from its start point to end point
    //Find the shortest distance between start point of the second cluster to this line segment.
    //Repeat for end point of second cluster to this line segment.
    //Then, pretend second cluster is a 2D line segment, from its start point to end point.
    //Find the shortest distance between start point of the first cluster to this line segment.
    //Repeat for end point of first cluster to this line segment.
    //If the shortest of these four distances is less than the cutoff, 
    //return true (the clusters are merge-compatible). else, return false.
    
    // Step 1: inspect (w_start1, t_start1) vs. line (w_start2, t_start2) => (w_end2, t_end2)
    double shortest_distance2 = ShortestDistanceSquared(w_start1, t_start1,
							w_start2, t_start2, 
							w_end2, t_end2);
    
    // Step 2: inspect (w_end1, t_end1) vs. line (w_start2, t_start2) => (w_end2, t_end2)
    double shortest_distance2_tmp = ShortestDistanceSquared(w_end1, t_end1,
							    w_start2, t_start2, 
							    w_end2, t_end2);
    
    shortest_distance2 = (shortest_distance2_tmp < shortest_distance2) ?
      shortest_distance2_tmp : shortest_distance2;
    
    // Step 3: inspect (w_start2, t_start2) vs. line (w_start1, t_start1) => (w_end1, t_end1)
    shortest_distance2_tmp = ShortestDistanceSquared(w_start2, t_start2,
						     w_start1, t_start1, 
						     w_end1, t_end1);

    shortest_distance2 = (shortest_distance2_tmp < shortest_distance2) ? 
      shortest_distance2_tmp : shortest_distance2;

    // Step 4: inspect (w_end2, t_end2) vs. line (w_start1, t_start1) => (w_end1, t_end1)
    shortest_distance2_tmp = ShortestDistanceSquared(w_end2, t_end2,
						     w_start1, t_start1, 
						     w_end1, t_end1);
    
    shortest_distance2 = (shortest_distance2_tmp < shortest_distance2) ? 
      shortest_distance2_tmp : shortest_distance2;

    bool compatible = shortest_distance2 < _max_2D_dist2;

    if(_verbose) {

      if(compatible) std::cout<<Form(" Compatible in distance (%g).\n",shortest_distance2);
      else std::cout<<Form(" NOT compatible in distance (%g).\n",shortest_distance2);
    
    }

    return compatible;


  }//end Bool function

  double CMAlgoShortestDistNonEndpoint::ShortestDistanceSquared(double point_x, double point_y, 
						     double start_x, double start_y,
						     double end_x,   double end_y  ) const {
    
    //This code finds the shortest distance between a point and a line segment.    
    //code based off sample from 
    //http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    //note to self: rewrite this with TVector2 and compare time differences... 
    //TVector2 code might be more understandable
    
    double distance_squared = 999999;
    
    // Line segment: from ("V") = (start_x, start_y) to ("W")=(end_x, end_y)
    double length_squared = pow((end_x - start_x), 2) + pow((end_y - start_y), 2);
    
    // Treat the case start & end point overlaps
    if(length_squared < _min_distance_unit) {
      if(_verbose){
	std::cout << std::endl;
	std::cout << Form(" Provided very short line segment: (%g,%g) => (%g,%g)",
			  start_x,start_y,end_x,end_y) << std::endl;
	std::cout << " Likely this means one of two clusters have start & end point identical." << std::endl;
	std::cout << " Check the cluster output!" << std::endl;
	std::cout << std::endl;
	std::cout << Form(" At this time, the algorithm uses a point (%g,%g)",start_x,start_y) << std::endl;
	std::cout << " to represent this cluster's location." << std::endl;
	std::cout << std::endl;
      }
      
      return (pow((point_x - start_x),2) + pow((point_y - start_y),2));
    }
    
    //Find shortest distance between point ("P")=(point_x,point_y) to this line segment
    double t = ( (point_x - start_x)*(end_x - start_x) + (point_y - start_y)*(end_y - start_y) ) / length_squared;
    
    //if(t<0.0) distance_squared = pow((point_x - start_x), 2) + pow((point_y - start_y), 2);
    
    //else if (t>1.0) distance_squared = pow((point_x - end_x), 2) + pow(point_y - end_y, 2);
    
    //else distance_squared = pow((point_x - (start_x + t*(end_x - start_x))), 2) + pow((point_y - (start_y + t*(end_y - start_y))),2);


    if(t>0.1 && t < 0.9) distance_squared = pow((point_x - (start_x + t*(end_x - start_x))), 2) + pow((point_y - (start_y + t*(end_y - start_y))),2);


    return distance_squared;
    
  }//end ShortestDistanceSquared function
  
}
#endif
