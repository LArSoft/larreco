/**
 * \file CMAlgoMergeAll.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoMergeAll
 *
 * @author david caratelli
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOMERGEALL_H
#define CMALGOMERGEALL_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoMergeAll
     Merges all clusters: maybe useful to test how well a cluster-separating
     algorithm has performed
  */
  class CMAlgoMergeAll: public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoMergeAll();
    
    /// Default destructor
    virtual ~CMAlgoMergeAll(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset(){}

    /// Function to report what's going on per merging
    virtual void Report();

  protected:

  };
}

#endif
/** @} */ // end of doxygen group 

