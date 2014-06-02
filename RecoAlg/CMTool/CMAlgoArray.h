/**
 * \file CMAlgoArray.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CMAlgoArray
 *
 * @author kazuhiro_NAME
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOARRAY_H
#define CMALGOARRAY_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoArray
     User implementation for CMAlgoArray class
     doxygen documentation!
  */
  class CMAlgoArray : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoArray();
    
    /// Default destructor
    virtual ~CMAlgoArray(){};

    /// A simple method to add a one merging step
    void AddAlgo(CBoolAlgoBase* algo,
		 bool ask_and = true)
    { 
      if(ask_and) _last_and_algo_index = _algo_array.size();
      _algo_array.push_back(algo);
      _ask_and.push_back(ask_and);
    }

    /**
       Optional function: called at the beginning of 1st iteration. This is called per event.
     */
    virtual void EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of event ... after the last merging iteration is over.
     */
    virtual void EventEnd();

    /**
       Optional function: Called at the beginning of each iteration over possible pairs of clusters
     */
    virtual void IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: Called at the end of each iteration over possible pairs of clusters
     */
    virtual void IterationEnd();

    /**
       Core function: given the CPAN input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2);

    

    
    /**
       Optional function: called after each Merge() function call by CMergeManager IFF
       CMergeManager is run with verbosity level kPerMerging. Maybe useful for debugging.
    */
    virtual void Report();
    
    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset();

  protected:

    /**
       A list of algorithms to be run over. Algorithms are executed in consecutive order
       in this vector, which is the order of calling AddMergeAlgo function. For each
       algorithm, the Merge() function return bool is taken as either AND or OR condition
       with the result of previous algorithm (or simply with "true" for the 1st algo). 
       Whether using AND or OR is a user configuration stored in _ask_and attribute
     */
    std::vector<CBoolAlgoBase*> _algo_array;

    /**
       A boolean vector that holds user configuration of whether asking AND or OR condition
       for algorithms in _algo_array.
     */
    std::vector<bool> _ask_and;

    /// Index of last AND condition algorithm to speed execution
    size_t _last_and_algo_index;
    
  };
}
#endif
/** @} */ // end of doxygen group 

