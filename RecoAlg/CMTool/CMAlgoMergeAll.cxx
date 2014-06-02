#ifndef CMALGOMERGEALL_CXX
#define CMALGOMERGEALL_CXX

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
    return true;
  }

  //-----------------------
  void CMAlgoMergeAll::Report()
  //-----------------------
  {

  }

}

#endif
