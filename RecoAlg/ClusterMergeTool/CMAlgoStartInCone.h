/**
 * \file CMAlgoStartInCone.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoStartInCone
 *
 * @author david
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOSTARTINCONE_H
#define CMALGOSTARTINCONE_H

#include <iostream>
#include "CBoolAlgoBase.h"
#include <math.h>

namespace cluster {
  
  /**
     \class CMAlgoStartInCone
     If start point of Cluster B in Cone region of Cluster A then merge
  */
  class CMAlgoStartInCone : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CMAlgoStartInCone();
    
    /// Default destructor
    virtual ~CMAlgoStartInCone(){};
  
    /// Merging Algorithm is Here
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);
    
    /// Method to re-configure the instance
    void reconfigure();

    /// Set Minimum number of hits for cone-cluster
    void SetMinHits(int n){ _NhitsMin = n; }

    /// Set Minimum number of hits for cone-cluster
    void SetMinLen(double l){ _lenMin = l; }

    /// Set Verbosity of messages
    void SetVerbose(bool verbosity){ _verbose = verbosity; }

    /// Set Debug for messages
    void SetDebug(bool debug){ _debug = debug; }

    /// Set Angle Compatibility betweeen the clusters
    void SetAngleCompat(double deg){ _angleCompat = deg; }

    /// Set Length Reach: How for out the cone extends as percent of cluster length
    void SetLengthReach(double frac){ _lengthReach = frac; }

  protected:

    int _NhitsMin;     /// Larger cluster which determines cone must have this many hits
    double _lenMin;    /// Larger cluster which determines cone must be at least this long
    bool _verbose;
    bool _debug;
    double _angleCompat; /// Two clusters must have direction within this value of each other
    double _lengthReach; ///How four out - as percent of cluster length - cone will extend from start point
    
  };

}  
#endif
/** @} */ // end of doxygen group 

