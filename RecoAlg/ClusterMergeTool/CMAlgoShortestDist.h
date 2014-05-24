/**
 * \file CMAlgoShortestDist.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoShortestDist
 *
 * @author davidkaleko
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOSHORTESTDIST_H
#define CMALGOSHORTESTDIST_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoShortestDist
     User defined class CMAlgoShortestDist ... these comments are used to generate
     doxygen documentation!
  */
  class CMAlgoShortestDist : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoShortestDist();
    
    /// Default destructor
    virtual ~CMAlgoShortestDist(){};
    
    /// Overloaded (from CBoolAlgoBase) Bool function
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);
    
    
    /// Method to set cut value in cm^2 for distance compatibility test
    void SetSquaredDistanceCut(double d) { _max_2D_dist2 = d; }

    /// Method to set verbose mode
    void SetVerbose(bool on) { _verbose = on; }

    /// Method to set debug mode
    void SetDebug(bool on) { _debug = on; }

    /// Set Minimum Number of Hits to consider Cluster
    void SetMinHits(int n) { _minHits = n; }
   
    /**
       Function to compute a distance between a 2D point (point_x, point_y) to a 2D finite line segment
       (start_x, start_y) => (end_x, end_y).
    */
    double ShortestDistanceSquared(double point_x, double point_y, 
				   double start_x, double start_y,
				   double end_x,   double end_y  ) const;



  protected:
    
    bool _verbose;         /// bool to suppress lots of output if you want

    bool _debug;         /// bool to suppress lots of output if you want

    int _minHits;        /// Min Number of hits for cluster to be considered
  
    double _min_distance_unit; /// minimum distance b/t start and end point of cluster to use it

    double _max_2D_dist2;      /// max distance b/t clusters for comability, in cm^2 (the main param of this algo)

  };
  

} //end namespace cluster

#endif
/** @} */ // end of doxygen group 
  
