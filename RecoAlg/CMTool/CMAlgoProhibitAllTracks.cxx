#ifndef CMALGOPROHIBITALLTRACKS_CXX
#define CMALGOPROHIBITALLTRACKS_CXX

#include "CMAlgoProhibitAllTracks.h"

namespace cluster {

  //-------------------------------------------------------
  CMAlgoProhibitAllTracks::CMAlgoProhibitAllTracks() : CBoolAlgoBase()
  //-------------------------------------------------------
  {
    SetMinEP(.990000);
  }

  //-----------------------------
  void CMAlgoProhibitAllTracks::Reset()
  //-----------------------------
  {

  }

  //------------------------------------------------------------------------------------------
  //void CMAlgoProhibitAllTracks::EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //------------------------------------------------------------------------------------------
  //{
  //
  //}

  //-------------------------------
  //void CMAlgoProhibitAllTracks::EventEnd()
  //-------------------------------
  //{
  //
  //}

  //-----------------------------------------------------------------------------------------------
  //void CMAlgoProhibitAllTracks::IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //-----------------------------------------------------------------------------------------------
  //{
  //
  //}

  //------------------------------------
  //void CMAlgoProhibitAllTracks::IterationEnd()
  //------------------------------------
  //{
  //
  //}
  
  //----------------------------------------------------------------
  bool CMAlgoProhibitAllTracks::Bool(const ClusterParamsAlg &cluster1,
			       const ClusterParamsAlg &cluster2)
  //----------------------------------------------------------------
  {
    //return true means don't prohibit these two clusters
    if(cluster1.GetParams().eigenvalue_principal > _min_EP ||
       cluster2.GetParams().eigenvalue_principal > _min_EP) 
      {
	if(_verbose) 
	  std::cout<<"Prohibiting clusters with EP's of "
		   <<cluster1.GetParams().eigenvalue_principal
		   <<" and "
		   <<cluster2.GetParams().eigenvalue_principal
		   <<std::endl;
	return true;
      }
    return false;
  }

  //------------------------------
  void CMAlgoProhibitAllTracks::Report()
  //------------------------------
  {

  }
    
}
#endif
