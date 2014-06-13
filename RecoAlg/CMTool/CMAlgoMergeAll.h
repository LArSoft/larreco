/**
 * \file CMAlgoMergeAll.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMAlgoMergeAll
 *
 * @author david caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMALGOMERGEALL_HH
#define CMALGOMERGEALL_HH

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

  protected:

  };
}

#endif
/** @} */ // end of doxygen group 

