/**
 * \file CMalgoPolyHitOverlap.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMAlgoPolyHitOverlap
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMALGOPOLYHITOVERLAP_HH
#define CMALGOPOLYHITOVERLAP_HH

#include <iostream>
#include "CBoolAlgoBase.h"
#include "Utilities/GeometryUtilities.h"

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

