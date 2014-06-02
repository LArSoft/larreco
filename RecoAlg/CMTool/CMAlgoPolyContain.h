/**
 * \file CMalgoPolyContain.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoPolyContain
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOPOLYCONTAIN_H
#define CMALGOPOLYCONTAIN_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMalgoPolyContain
     Merge Polygons if one is completely inside the other
  */
  class CMAlgoPolyContain : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoPolyContain();
    
    /// Default destructor
    virtual ~CMAlgoPolyContain(){};
 
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

