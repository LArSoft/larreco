#ifndef CMALGOFAKE_CXX
#define CMALGOFAKE_CXX

#include "CMAlgoFake.h"

namespace cluster {

  //----------------------------------------
  CMAlgoFake::CMAlgoFake() : CBoolAlgoBase()
  //----------------------------------------
  {
    _flip = false;
    _ctr  = 0;
    // Nothing to be done in the base class
  }

  //--------------------------------------------------------
  bool CMAlgoFake::Bool(const ClusterParamsAlg &cluster1,
			const ClusterParamsAlg &cluster2)
  //--------------------------------------------------------
  {
    _ctr++;
    if( (_ctr%64) == 0)
      _flip = (!_flip);
    return _flip;
  }

  //-----------------------
  void CMAlgoFake::Report()
  //-----------------------
  {
    std::cout<< "  I am just flpping every 64 counts... " << std::endl;
  }

}

#endif
