#ifndef CMALGOPOLYCONTAIN_CXX
#define CMALGOPOLYCONTAIN_CXX

#include "CMAlgoPolyContain.h"

namespace cluster {

  CMAlgoPolyContain::CMAlgoPolyContain()
  {
    // Nothing to be done in the base class
    this->reconfigure();
  }


  void CMAlgoPolyContain::reconfigure(){

    //not sure what needs to be reset/reconfigured for this algo
    
  }//end reconfigure function

  
  bool CMAlgoPolyContain::Bool(const ClusterParamsAlg &cluster1,
			       const ClusterParamsAlg &cluster2)
  {

    if ( (cluster1.GetParams().PolyObject.Size() < 3) or (cluster2.GetParams().PolyObject.Size() < 3) )
      return false;

    //if either polygon is fully contained in other
    //then return true! --> MERGE!
    if ( (cluster1.GetParams().PolyObject.Contained(cluster2.GetParams().PolyObject)) or
	 (cluster2.GetParams().PolyObject.Contained(cluster1.GetParams().PolyObject)) )
      return true;
    else
      return false;
  }
  

}

#endif
