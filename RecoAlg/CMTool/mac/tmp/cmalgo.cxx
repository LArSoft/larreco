#ifndef CMALGO_CLASS_NAME_CXX
#define CMALGO_CLASS_NAME_CXX

#include "CMAlgo_Class_Name.h"

namespace cluster {

  //-------------------------------------------------------
  CMAlgo_Class_Name::CMAlgo_Class_Name() : CBoolAlgoBase()
  //-------------------------------------------------------
  {

  }

  //-----------------------------
  void CMAlgo_Class_Name::Reset()
  //-----------------------------
  {

  }

  //------------------------------------------------------------------------------------------
  void CMAlgo_Class_Name::EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //------------------------------------------------------------------------------------------
  {

  }

  //-------------------------------
  void CMAlgo_Class_Name::EventEnd()
  //-------------------------------
  {

  }

  //-----------------------------------------------------------------------------------------------
  void CMAlgo_Class_Name::IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //-----------------------------------------------------------------------------------------------
  {

  }

  //------------------------------------
  void CMAlgo_Class_Name::IterationEnd()
  //------------------------------------
  {

  }
  
  //----------------------------------------------------------------
  bool CMAlgo_Class_Name::Bool(const ClusterParamsAlg &cluster1,
			       const ClusterParamsAlg &cluster2)
  //----------------------------------------------------------------
  {
    return false;
  }

  //------------------------------
  void CMAlgo_Class_Name::Report()
  //------------------------------
  {

  }
    
}
#endif
