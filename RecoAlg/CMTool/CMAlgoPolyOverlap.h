/**
 * \file CMalgoPolyOverlap.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMAlgoPolyOverlap
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMALGOPOLYOVERLAP_HH
#define CMALGOPOLYOVERLAP_HH

#include <iostream>
#include "CBoolAlgoBase.h"
#include "Utilities/GeometryUtilities.h"


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

    //both clusters must have > this # of hits to be considered for merging
    void SetMinNumHits(size_t nhits) { _min_hits = nhits; }

    /// Method to re-configure the instance
    void reconfigure();

  private:
    
    bool _debug;
    size_t _min_hits;
  };
}

#endif
/** @} */ // end of doxygen group 

