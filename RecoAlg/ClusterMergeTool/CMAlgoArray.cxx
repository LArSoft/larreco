#ifndef CMALGOARRAY_CXX
#define CMALGOARRAY_CXX

#include "CMAlgoArray.h"

namespace cluster {

  //------------------------------------------
  CMAlgoArray::CMAlgoArray() : CBoolAlgoBase()
  //------------------------------------------
  {
    _algo_array.clear();
    _ask_and.clear();
  }

  //-----------------------
  void CMAlgoArray::Reset()
  //-----------------------
  {
    for(auto &algo : _algo_array) algo->Reset();
  }

  //-------------------------------------------------------------------------------------
  void CMAlgoArray::EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //-------------------------------------------------------------------------------------
  {
    for(auto &algo : _algo_array) algo->EventBegin(clusters);
  }

  //--------------------------
  void CMAlgoArray::EventEnd()
  //--------------------------
  {
    for(auto &algo : _algo_array) algo->EventEnd();
  }

  //-------------------------------------------------------------------------------------
  void CMAlgoArray::IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
  //-------------------------------------------------------------------------------------
  {
    for(auto &algo : _algo_array) algo->IterationBegin(clusters);
  }

  //--------------------------
  void CMAlgoArray::IterationEnd()
  //--------------------------
  {
    for(auto &algo : _algo_array) algo->IterationEnd();
  }

  //---------------------------------------------------------
  bool CMAlgoArray::Bool(const ClusterParamsAlg &cluster1,
			 const ClusterParamsAlg &cluster2)
  //---------------------------------------------------------
  {
    bool status = true;
    
    for(size_t i=0; i<_algo_array.size(); ++i) {

      if(!i) status = _algo_array.at(i)->Bool(cluster1,cluster2);

      else {

	//
	// Break before executing algo if possible
	//

	// Case 1: if AND and status==false, then break
	if(  _ask_and.at(i) && !status ) break;

	// Case 2: if OR and status==true, then break
	if( !_ask_and.at(i) &&  status ) break;

	// Case 3: the remaining algorithms are all OR condition and stauts==true
	if( i > _last_and_algo_index && status) break;

	//
	// Execute algorithm
	//
	if( _ask_and.at(i) ) 

	  status = status && _algo_array.at(i)->Bool(cluster1,cluster2);

	else 
	  
	  status = status || _algo_array.at(i)->Bool(cluster1,cluster2);
	
      }
    }

    return status;
  }

  //------------------------
  void CMAlgoArray::Report()
  //------------------------
  {
    for(auto &algo : _algo_array) algo->Report();
  }
    
}
#endif
