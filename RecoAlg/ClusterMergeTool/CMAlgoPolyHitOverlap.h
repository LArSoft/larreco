/**
 * \file CMalgoPolyHitOverlap.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoPolyHitOverlap
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOPOLYHITOVERLAP_H
#define CMALGOPOLYHITOVERLAP_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMalgoPolyContain
     Merge Polygons if one is completely inside the other
  */
  class CMAlgoPolyHitOverlap : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoPolyHitOverlap();
    
    /// Default destructor
    virtual ~CMAlgoPolyHitOverlap(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Method to re-configure the instance
    void reconfigure();

  };
}

#endif
/** @} */ // end of doxygen group 

