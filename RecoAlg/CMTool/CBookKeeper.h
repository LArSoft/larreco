/**
 * \file CBookKeeper.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CBookKeeper
 *
 * @author kazuhiro
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CBOOKKEEPER_HH
#define CBOOKKEEPER_HH

#include <iostream>
#include <vector>
#include <TString.h>
#include "RecoAlg/ClusterRecoUtil/CRUException.h"

namespace cluster {

  /**
     \class CBookKeeper
     A utility class for CMergeManager which merge clusters using merging algorithms.
     One of major task for CMergeManager is to keep track of which clusters to be merged
     in the original input. CBookKeeper handles this part. It works with indexes. The user
     (primarily CMergeManager) provides number of clusters to CBookKeeper's constructor.
     Then it can ask CBookKeeper to merge two specific clusters by specifying index number
     of the cluster which has to be smaller than the previously specified number of clusters.
     CBookKeeper keeps track of which clusters are asked to merge together, and it can be
     asked to return a vector of merged cluster indexes.
  */
  class CBookKeeper : public std::vector<unsigned short>{
    
  public:
    
    /// Default constructor
    CBookKeeper(unsigned short nclusters=0);
    
    /// Default destructor
    virtual ~CBookKeeper(){};

    /// Reset method
    void Reset(unsigned short nclusters=0);

    /// Method to set a pair of clusters to prohibit from merging
    void ProhibitMerge(unsigned short index1, unsigned short index2);

    /// Method to inqury if a combination is prohibited to merge
    bool MergeAllowed(unsigned short index1, unsigned short index2);

    /// Method to merge 2 clusters via index numbers
    void Merge(unsigned short index1, unsigned short index2);

    /**
       Method to retrieve a vector of cluster indexes which 
       is merged with the input cluster index. All indexes here
       are regarding the original cluster index.
    */
    std::vector<unsigned short> GetMergedSet(unsigned short index1) const;

    /**
       Method to ask if a given 2 clusters are already merged.
       This method is expected to be much faster than obtaining
       a merged cluster sets from GetMergedIndexes and check if
       two clusters are merged.
    */
    bool IsMerged(unsigned short index1, unsigned short index2) const;

    /**
       A method to get the full result. The return is a vector
       of merged cluster indexes (which is a vector of original cluster indexes).
    */
    void PassResult(std::vector<std::vector<unsigned short> > &result) const;


    std::vector<std::vector<unsigned short> > GetResult() const
    { 
      std::vector<std::vector<unsigned short> > result; 
      PassResult(result);
      return result;
    }

    /**
       Method to combine with another CBookKeeper instance.
     */
    
    void Combine(const CBookKeeper &another);

    void Report() const;

  protected:

    /**
       A 2D-map vector that stores pair of clusters for which merging is prohibited
     */
    std::vector<std::vector<bool> > _prohibit_merge;

    /// Number of output clusters
    size_t _out_cluster_count;
    
  };
  
}

#endif
/** @} */ // end of doxygen group 

