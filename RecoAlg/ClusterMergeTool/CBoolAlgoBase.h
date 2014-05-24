/**
 * \file CBoolAlgoBase.h
 *
 * \ingroup ClusterCluster
 * 
 * \brief Class def header for a class CAlgoBase, CBoolAlgoBase and CSeparateAlgoBase
 *
 * @author kazuhiro
 */

/** \addtogroup ClusterCluster

    @{*/
#ifndef CBOOLALGOBASE_H
#define CBOOLALGOBASE_H

#include <iostream>
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"

namespace cluster {

  /**
     \class CBoolAlgoBase
     An abstract base class for merging/mathcing algorithm. Having this base class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CBoolAlgoBase {
    
  public:
    
    /// Default constructor
    //CBoolAlgoBase(){ _fout = 0; }
    CBoolAlgoBase(){}
    
    /// Default destructor
    virtual ~CBoolAlgoBase(){}

    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset(){}

    /**
       Optional function: called at the beginning of 1st iteration. This is called per event.
     */
    virtual void EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
    { if(clusters.size()) return; }

    /**
       Optional function: called at the end of event ... after the last merging iteration is over.
     */
    virtual void EventEnd()
    {return;}
 
    /**
       Optional function: called at the beggining of each iteration over all pairs of clusters. 
       This provides all clusters' information in case the algorithm need them. Note this
       is called per iteration which may be more than once per event.
     */
    virtual void IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters)
    { if(clusters.size()) return;}

    /**
       Optional function: called at the end of each iteration over all pairs of clusters.
     */
    virtual void IterationEnd()
    {return; }

    /**
       Core function: given the CPAN input, return whether a cluster should be
       merged or not.
    */

    virtual bool Bool(const ClusterParamsAlg &cluster1,
		      const ClusterParamsAlg &cluster2)
    { if(cluster1.Plane() != cluster2.Plane()) return false;
      else return true;
    }

    /**
       Optional function: called after each Merge() function call by CMergeManager IFF
       CMergeManager is run with verbosity level kPerMerging. Maybe useful for debugging.
     */
    virtual void Report()
    {return;}

    /// Setter function for an output plot TFile pointer
    //void SetAnaFile(TFile* fout) { _fout = fout; }

  protected:

    /// TFile pointer to an output file
    //TFile* _fout;

  };

}

#endif
/** @} */ // end of doxygen group 

