#ifndef CMALGOSTARTTRACK_CXX
#define CMALGOSTARTTRACK_CXX

#include "CMAlgoStartTrack.h"

namespace cluster {

  //-------------------------------------------------------
  CMAlgoStartTrack::CMAlgoStartTrack() : CBoolAlgoBase()
  //-------------------------------------------------------
  {

  }

  //-----------------------------
  void CMAlgoStartTrack::Reset()
  //-----------------------------
  {

  }

  //------------------------------------------------------------------------------------------
  void CMAlgoStartTrack::EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //------------------------------------------------------------------------------------------
  {

  }

  //-------------------------------
  void CMAlgoStartTrack::EventEnd()
  //-------------------------------
  {

  }

  //-----------------------------------------------------------------------------------------------
  void CMAlgoStartTrack::IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //-----------------------------------------------------------------------------------------------
  {

  }

  //------------------------------------
  void CMAlgoStartTrack::IterationEnd()
  //------------------------------------
  {

  }
  
  //----------------------------------------------------------------
  bool CMAlgoStartTrack::Bool(const ClusterParamsAlg &cluster1,
			       const ClusterParamsAlg &cluster2)
  //----------------------------------------------------------------
  {
    return false;
  }

  //------------------------------
  void CMAlgoStartTrack::Report()
  //------------------------------
  {

  }
    
}
#endif
