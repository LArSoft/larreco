/**
 * \file CMAlgoStartTrack.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief This merge algo is looking for short tracks from the 
 *        start of a shower that are overlapping a blob that is
 *        a cluster belonging to the same shower.
 *
 * @author davidkaleko
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CMALGOSTARTTRACK_H
#define CMALGOSTARTTRACK_H

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoStartTrack
     User implementation for CBoolAlgoBase class
     doxygen documentation!
  */
  class CMAlgoStartTrack : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoStartTrack();
    
    /// Default destructor
    virtual ~CMAlgoStartTrack(){};

    /**
       Optional function: called at the beginning of 1st iteration. This is called per event.
     */
    virtual void EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of event ... after the last merging iteration is over.
     */
    virtual void EventEnd();
 
    /**
       Optional function: called at the beggining of each iteration over all pairs of clusters. 
       This provides all clusters' information in case the algorithm need them. Note this
       is called per iteration which may be more than once per event.
     */
    virtual void IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of each iteration over all pairs of clusters.
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
    
  };
}
#endif
/** @} */ // end of doxygen group 

