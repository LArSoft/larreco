/**
 * \file CMalgoPolyOverlap.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoPolyOverlap
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterCluster    

    @{*/
#ifndef CMALGOPOLYOVERLAP_H
#define CMALGOPOLYOVERLAP_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMalgoPolyContain
     Merge Polygons if the two overlap even partially
  */
  class CMAlgoPolyOverlap : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoPolyOverlap();
    
    /// Default destructor
    virtual ~CMAlgoPolyOverlap(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    void SetDebug(bool debug) { _debug = debug; }

    void SetVerbose(bool verbose) { _verbose = verbose; }

    /// Method to re-configure the instance
    void reconfigure();

  private:
    
    bool _verbose;
    bool _debug;

  };
}

#endif
/** @} */ // end of doxygen group 

