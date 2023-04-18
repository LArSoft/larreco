/**
 * \file CMergeManager.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CMergeManager
 *
 * @author kazuhiro
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CMERGEMANAGER_H
#define RECOTOOL_CMERGEMANAGER_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CMManagerBase.h"
#include "larreco/RecoAlg/CMTool/CMToolBase/CMergeBookKeeper.h"

#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"
#include <vector>

namespace cmtool {

  class CBoolAlgoBase;

  /**
     \class CMergeManager
     A class that instantiates merging algorithm(s) and run.
     The book-keeping of merged cluster sets are done by CMergeBookKeeper.
  */
  class CMergeManager : public CMManagerBase {
  public:
    CMergeManager();

    /// Default destructor
    virtual ~CMergeManager() {}

    /// Method to reset itself
    virtual void Reset();

    /// A simple method to add an algorithm for merging
    void AddMergeAlgo(CBoolAlgoBase* algo) { _merge_algo = algo; }

    /// A simple method to add an algorithm for separation
    void AddSeparateAlgo(CBoolAlgoBase* algo) { _separate_algo = algo; }

    /// A method to obtain output clusters
    const std::vector<cluster::ClusterParamsAlg>& GetClusters() const { return _out_clusters; }

    /// A method to obtain book keeper
    const CMergeBookKeeper& GetBookKeeper() const { return _book_keeper; }

  protected:
    //
    // FMWK functions override
    //

    /// FMWK function called @ beginning of Process()
    virtual void EventBegin();

    /// FMWK function called @ beginning of iterative loop inside Process()
    virtual void IterationBegin();

    /// FMWK function called @ iterative loop inside Process()
    virtual bool IterationProcess(util::GeometryUtilities const& gser);

    /// FMWK function called @ end of iterative loop inside Process()
    virtual void IterationEnd();

    /// FMWK function called @ end of Process()
    virtual void EventEnd();

  protected:
    void RunMerge(const std::vector<cluster::ClusterParamsAlg>& in_clusters,
                  CMergeBookKeeper& book_keeper) const;

    void RunMerge(const std::vector<cluster::ClusterParamsAlg>& in_clusters,
                  const std::vector<bool>& merge_flag,
                  CMergeBookKeeper& book_keeper) const;

    void RunSeparate(const std::vector<cluster::ClusterParamsAlg>& in_clusters,
                     CMergeBookKeeper& book_keeper) const;

  protected:
    /// Output clusters
    std::vector<cluster::ClusterParamsAlg> _out_clusters;

    /// Book keeper instance
    CMergeBookKeeper _book_keeper;

    /// Merging algorithm
    ::cmtool::CBoolAlgoBase* _merge_algo;

    /// Separation algorithm
    ::cmtool::CBoolAlgoBase* _separate_algo;

    size_t _iter_ctr;

    std::vector<CMergeBookKeeper> _book_keeper_v;

    std::vector<std::vector<unsigned short>> _tmp_merged_indexes;

    std::vector<cluster::ClusterParamsAlg> _tmp_merged_clusters;
  };
}

#endif
/** @} */ // end of doxygen group
