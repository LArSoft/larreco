#ifndef CMALGOPOLYHITOVERLAP_CXX
#define CMALGOPOLYHITOVERLAP_CXX

#include "CMAlgoPolyHitOverlap.h"

namespace cluster {

  CMAlgoPolyHitOverlap::CMAlgoPolyHitOverlap()
  {
    // Nothing to be done in the base class
    this->reconfigure();
  }


  void CMAlgoPolyHitOverlap::reconfigure(){

    //not sure what needs to be reset/reconfigured for this algo
    
  }//end reconfigure function

  
  bool CMAlgoPolyHitOverlap::Bool(const ClusterParamsAlg &cluster1,
				  const ClusterParamsAlg &cluster2)
  {

    //Check and see if a certain fraction of hits of a cluster
    //lie within polygon boundary of other cluster
    

    return false;
  }


}

#endif
