/**
 * \file CMalgoPolyStartInPoly.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMAlgoStartInPoly
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMALGOSTARTINPOLY_HH
#define CMALGOSTARTINPOLY_HH

#include <iostream>
#include "CBoolAlgoBase.h"
#include "Utilities/GeometryUtilities.h"

namespace cluster {
  /**
     \class CMalgoStartInPoly
     If start point of one cluster inside other's polygon -> merge
  */
  class CMAlgoStartInPoly : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoStartInPoly();
    
    /// Default destructor
    virtual ~CMAlgoStartInPoly(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */

    /// Method to set cut value on minimum number of hits considered
    void SetMinHitsCut(size_t n) { _MinHits = n; }

    void SetDebug(bool debug) { _debug = debug; }

    /// Merging Algorithm is Here
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Method to re-configure the instance
    void reconfigure();

  protected:

    double _wire_2_cm, _time_2_cm; /// Conversion factors ogtten from GeometryUtilities
    size_t _MinHits; /// Minimum number of hits for cluster whose start point is being considered. We want it to be a good start point...
    bool _debug;
  };
}

#endif
/** @} */ // end of doxygen group 

