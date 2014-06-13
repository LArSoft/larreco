/**
 * \file CMergeManager.hh
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class CMergeManager
 *
 * @author kazuhiro
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef CMERGEMANAGER_HH
#define CMERGEMANAGER_HH

#include <iostream>

#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include "CBookKeeper.h"
#include "CBoolAlgoBase.h"

namespace cluster {

  /**
     \class CMergeManager
     A class that instantiates merging algorithm(s) and run.
     The book-keeping of merged cluster sets are done by CBookKeeper.
  */
  class CMergeManager{

  public:

    /// Enum to specify the priority for running a merging.
    enum CMergePriority_t {
      /// Merge clusters along the vector index (from low to high)
      kIndex,
      /// Merge clusters with large areas first 
      kPolyArea,
      /// Merge clusters with high charge first
      kChargeSum,
      /// Merge clusters with many hits first
      kNHits
    };

    /// Enum to specify message output level
    enum CMergeMSGLevel_t {
      /// Extremely verbose (cout per individual algorithm execution)
      kPerMerging,
      /// Very verbose (cout per set of algorithms execution)
      kPerAlgoSet,
      /// Somewhat verbose (cout per merging iteration)
      kPerIteration,
      /// Bit verbose (cout per event)
      kPerEvent,
      /// Normal
      kNone
    };
    
  public:
    
    /// Default constructor
    CMergeManager(CMergePriority_t priority = kNHits);
    
    /// Default destructor
    virtual ~CMergeManager(){};

    /// Method to enable debug mode (lots of couts)
    void DebugMode(CMergeMSGLevel_t level) {_debug_mode=level;}

    /// Switch to continue merging till converges
    void MergeTillConverge(bool doit=true) {_merge_till_converge = doit;}

    /// Choose ordering for clusters to be merged
    void SetMergePriority(CMergePriority_t level) { _priority=level; }

    /// Method to reset itself
    void Reset();

    /// A simple method to add an algorithm for merging
    void AddMergeAlgo(CBoolAlgoBase* algo) {_merge_algo = algo;}

    /// A simple method to add an algorithm for separation
    void AddSeparateAlgo(CBoolAlgoBase* algo) {_separate_algo = algo;}

    /// A simple method to add a cluster
    void SetClusters(const std::vector<std::vector<util::PxHit> > &clusters);

    /// A method to execute merging algorithms
    void Process();

    /// A method to obtain output clusters
    const std::vector<cluster::ClusterParamsAlg>& GetClusters() const { return _out_clusters; }

    /// A method to obtain book keeper
    const CBookKeeper& GetBookKeeper() const { return _book_keeper; }

    /// A setter for an analysis output file
    void SetAnaFile(TFile* fout) { _fout = fout; }

  protected:

    void RunMerge(const std::vector<cluster::ClusterParamsAlg > &in_clusters,
		  CBookKeeper &book_keeper) const;

    void RunMerge(const std::vector<cluster::ClusterParamsAlg > &in_clusters,
		  const std::vector<bool> &merge_flag,
		  CBookKeeper &book_keeper) const;

    void RunSeparate(const std::vector<cluster::ClusterParamsAlg > &in_clusters,
		     CBookKeeper &book_keeper) const;

  protected:

    /// Iterative approach for merging
    bool _merge_till_converge;

    /// Debug mode switch
    CMergeMSGLevel_t _debug_mode;

    /// Input clusters
    std::vector<cluster::ClusterParamsAlg> _in_clusters;

    /// Output clusters
    std::vector<cluster::ClusterParamsAlg> _out_clusters;

    /// Book keeper instance
    CBookKeeper _book_keeper;

    /// Merging algorithm
    ::cluster::CBoolAlgoBase* _merge_algo;

    /// Separation algorithm
    ::cluster::CBoolAlgoBase* _separate_algo;
    
    /// Merging priority type
    CMergePriority_t _priority;

    /// Output analysis plot TFile
    TFile* _fout;

  };
}

#endif
/** @} */ // end of doxygen group 

