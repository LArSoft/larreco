/**
 * \file CMAlgoPolyShortestDist.hh
 *
 * \ingroup ClusterStudy
 * 
 * \brief Class def header for a class CMAlgoPolyShortestDist
 *
 * @author davidkaleko_NAME
 */

/** \addtogroup ClusterStudy

    @{*/
#ifndef CMALGOPOLYSHORTESTDIST_HH
#define CMALGOPOLYSHORTESTDIST_HH

#include <iostream>
#include "CBoolAlgoBase.h"

namespace cluster {
  /**
     \class CMAlgoPolyShortestDist
     User implementation for CBoolAlgoBase class
     doxygen documentation!
  */
  class CMAlgoPolyShortestDist : public CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    CMAlgoPolyShortestDist();
    
    /// Default destructor
    virtual ~CMAlgoPolyShortestDist(){};

    /**
       Optional function: called at the beginning of 1st iteration. This is called per event.
    */
    virtual void EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

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

    //both clusters must have > this # of hits to be considered for merging
    void SetMinNumHits(size_t nhits) { _min_hits = nhits; }

    void SetMinDistSquared(double dist) { _dist_sqrd_cut = dist; }

    void SetDebug(bool flag) { _debug = flag; }

  private:

    size_t _min_hits;

    double _dist_sqrd_cut;

    bool _debug;

    double tmp_min_dist;
  };
}
#endif
/** @} */ // end of doxygen group 

