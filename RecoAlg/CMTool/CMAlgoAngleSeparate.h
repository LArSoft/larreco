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
#ifndef CMALGOANGLESEPARATE_HH
#define CMALGOANGLESEPARATE_HH

#include <math.h>
#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoAngleSeparate
     Track Prohibit algorithm: if the angle between the direction of a cluster (end-start) and the line connecting the cluster's start point and the start point of t a second cluster is too large, then probihit merging between the two clusters. The first cluster needs to be a "good" and "large" cluster
     algorithm has performed
  */
  class CMAlgoAngleSeparate: public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoAngleSeparate();
    
    /// Default destructor
    virtual ~CMAlgoAngleSeparate(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Set Debug Mode on or off
    void SetDebug(bool on) { _debug = on; }

    /// Set Max Angle Separation for separation
    void SetMaxAngleSep(float angle) { _MaxAngle = angle; }

    /// Set Minimum length for "big" cluster
    void SetMinLength(float len) { _MinLen = len; }

    /// SetMinimum number of hits for small cluster
    void SetMinHits(size_t n) { _minHits = n; }

    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset(){}

    /// Function to report what's going on per merging
    virtual void Report();

  protected:

    bool _debug;
    float _MaxAngle;
    float _MinLen;
    size_t _minHits;

  };
}

#endif
/** @} */ // end of doxygen group 

