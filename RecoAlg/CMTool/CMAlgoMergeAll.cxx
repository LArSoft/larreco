#ifndef CMALGOMERGEALL_CC
#define CMALGOMERGEALL_CC

#include "CMAlgoMergeAll.h"

namespace cluster {

  //----------------------------------------
  CMAlgoMergeAll::CMAlgoMergeAll() : CBoolAlgoBase()
  //----------------------------------------
  {

  }

  //--------------------------------------------------------
  bool CMAlgoMergeAll::Bool(const ClusterParamsAlg &cluster1,
			    const ClusterParamsAlg &cluster2)
  //--------------------------------------------------------
  {
    if(cluster1.GetNHits() && cluster2.GetNHits()) return true;
    else return false;
  }
  
}

#endif
