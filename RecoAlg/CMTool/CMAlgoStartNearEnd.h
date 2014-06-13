/**
 * \file CMAlgoStartNearEnd.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMAlgoStartNearEnd
 *
 * @author david caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMALGOSTARTNEAREND_HH
#define CMALGOSTARTNEAREND_HH

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoStartNearEnd
     Merge if start point of one is near end point of another
     and require angle compatibility
  */
  class CMAlgoStartNearEnd : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoStartNearEnd();
    
    /// Default destructor
    virtual ~CMAlgoStartNearEnd(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */

    void SetMaxStartEndSeparation(double d) { _separation=d; }

    void SetMaxAngle(double a) { _maxopeningangle=a; }

    void SetMinHits(size_t n) { _MinHits=n; }

    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset(){}

    /// Function to report what's going on per merging
    virtual void Report();

  protected:

    double _maxopeningangle;
    double _separation;
    size_t _MinHits;

  };
}

#endif
/** @} */ // end of doxygen group 

