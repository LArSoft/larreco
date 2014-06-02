/**
 * \file CMAlgoFake.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoTrackSeparate
 *
 * @author david caratelli
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOTRACKSEPARATE_H
#define CMALGOTRACKSEPARATE_H

#include <iostream>
#include "math.h"
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoTrackSeparate
  */
  class CMAlgoTrackSeparate : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoTrackSeparate();
    
    /// Default destructor
    virtual ~CMAlgoTrackSeparate(){};

    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    void SetVerbose(bool on) { _verbose = on; }

    void SetDebug(bool on) { _debug = on; }

    void SetMinNumHits(int n) { _MinNumHits = n; }

    void SetMinAngleDiff(float anglesep) { _MinAngleDiff = anglesep; }
    
    void SetMaxOpeningAngle(float maxangle) { _MaxOpeningAngle = maxangle; }

    void SetMinLength(float minlength) { _MinLength = minlength; }

    void SetMinPolyHitDensity(float mindensity) { _MinDensity = mindensity; }

    void SetMaxWidth(float maxwidth) { _MaxWidth = maxwidth; }

    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset(){}

    /// Function to report what's going on per merging
    virtual void Report();

  protected:

    bool _verbose;
    bool _debug;
    int _MinNumHits;
    float _MinAngleDiff;
    float _MaxOpeningAngle;
    float _MinLength;
    float _MinDensity;
    float _MaxWidth;
  };
}

#endif
/** @} */ // end of doxygen group 

