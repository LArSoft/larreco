/**
 * \file CMAlgoProhibitAllTracks.hh
 *
 * \ingroup ClusterStudy
 * 
 * \brief Class def header for a class CMAlgoProhibitAllTracks
 *
 * @author davidkaleko_NAME
 */

/** \addtogroup ClusterStudy

    @{*/
#ifndef CMALGOPROHIBITALLTRACKS_HH
#define CMALGOPROHIBITALLTRACKS_HH

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoProhibitAllTracks
     User implementation for CBoolAlgoBase class
     doxygen documentation!
  */
  class CMAlgoProhibitAllTracks : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoProhibitAllTracks();
    
    /// Default destructor
    virtual ~CMAlgoProhibitAllTracks(){};

    /**
       Optional function: called at the beginning of 1st iteration. This is called per event.
     */
    //virtual void EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of event ... after the last merging iteration is over.
     */
    //virtual void EventEnd();
 
    /**
       Optional function: called at the beggining of each iteration over all pairs of clusters. 
       This provides all clusters' information in case the algorithm need them. Note this
       is called per iteration which may be more than once per event.
     */
    //virtual void IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of each iteration over all pairs of clusters.
     */
    //virtual void IterationEnd();
    
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
    

    void SetMinEP(double value) { _min_EP = value; }

  protected:

    double _min_EP;

  };


}
#endif
/** @} */ // end of doxygen group 

