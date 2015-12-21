/**
 * \file CBAlgoShortestDistSmallCluster.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CBAlgoShortestDistSmallCluster
 *
 * @author davidkaleko
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CBALGOSHORTESTDISTSMALLCLUSTER_H
#define CBALGOSHORTESTDISTSMALLCLUSTER_H

#include <iostream>
#include "larreco/RecoAlg/CMTool/CMToolBase/CBoolAlgoBase.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace cmtool {
  /**
     \class CBAlgoShortestDistSmallCluster
     User defined class CBAlgoShortestDistSmallCluster ... these comments are used to generate
     doxygen documentation!
  */
  class CBAlgoShortestDistSmallCluster : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CBAlgoShortestDistSmallCluster();
    
    /// Default destructor
    virtual ~CBAlgoShortestDistSmallCluster(){};
    
    /// Overloaded (from CBoolAlgoBase) Bool function
    virtual bool Bool(const ::cluster::ClusterParamsAlg &cluster1,
		      const ::cluster::ClusterParamsAlg &cluster2);
    
    
    /// Method to set cut value in cm^2 for distance compatibility test
    void SetSquaredDistanceCut(double d) { _max_2D_dist2 = d; }

    /// Method to set debug mode
    void SetDebug(bool on) { _debug = on; }

    /// Set Minimum Number of Hits to consider Cluster
    void SetMinHits(size_t n) { _minHits = n; }
   
    /// Set Maximum Number of Hits to consider Cluster
    void SetMaxHits(size_t n) { _maxHits = n; }
   
    /**
       Function to compute a distance between a 2D point (point_x, point_y) to a 2D finite line segment
       (start_x, start_y) => (end_x, end_y).
    */
    double ShortestDistanceSquared(double point_x, double point_y, 
				   double start_x, double start_y,
				   double end_x,   double end_y  ) const;



  protected:
    
    bool _debug;         /// bool to suppress lots of output if you want

    size_t _minHits;        /// Min Number of hits for cluster to be considered

    size_t _maxHits;        /// Max Number of hits for cluster to be considered

    double _wire_2_cm, _time_2_cm; /// Conversion factors ogtten from GeometryUtilities

    double _min_distance_unit; /// minimum distance b/t start and end point of cluster to use it

    double _max_2D_dist2;      /// max distance b/t clusters for comability, in cm^2 (the main param of this algo)

  };
  

} //end namespace cluster

#endif
/** @} */ // end of doxygen group 
  
