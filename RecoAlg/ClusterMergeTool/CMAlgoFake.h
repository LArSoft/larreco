/**
 * \file CMAlgoFake.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoFake
 *
 * @author kazuhiro
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOFAKE_H
#define CMALGOFAKE_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoFake
     An abstract fake class for merging algorithm. Having this fake class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CMAlgoFake : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoFake();
    
    /// Default destructor
    virtual ~CMAlgoFake(){};
 
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

    bool _flip;
    int _ctr;
  };
}

#endif
/** @} */ // end of doxygen group 

