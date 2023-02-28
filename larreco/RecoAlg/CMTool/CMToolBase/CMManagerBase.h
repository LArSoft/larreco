/**
 * \file CMManagerBase.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CMManagerBase
 *
 * @author kazuhiro
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CMMANAGERBASE_H
#define RECOTOOL_CMMANAGERBASE_H

#include <map>
#include <set>
#include <stddef.h>
#include <vector>

#include "RtypesCore.h"
class TFile;

#include "lardata/Utilities/PxUtils.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"

namespace util {
  class GeometryUtilities;
}

namespace cmtool {

  class CPriorityAlgoBase;

  /**
     \class CMManagerBase
     A class that instantiates merging algorithm(s) and run.
     The book-keeping of merged cluster sets are done by CMergeBookKeeper.
  */
  class CMManagerBase {
  public:
    /// Enum to specify message output level
    enum CMMSGLevel_t {
      /// Extremely verbose (cout per individual algorithm execution)
      kPerMerging,
      /// Somewhat verbose (cout per merging iteration)
      kPerIteration,
      /// Bit verbose (cout per event)
      kPerEvent,
      /// Normal
      kNone
    };

    /// Default constructor
    CMManagerBase();

    /// Default destructor
    virtual ~CMManagerBase() = default;

    /// Method to enable debug mode (lots of couts)
    void DebugMode(CMMSGLevel_t level) { _debug_mode = level; }

    /// Method to enable timing profile cout
    void ReportTimings(bool time_report = true) { _time_report = time_report; }

    /// Method to reset itself
    void Reset();

    /// Setter to add an algorithm for priority determination
    void AddPriorityAlgo(CPriorityAlgoBase* algo) { _priority_algo = algo; }

    /// Switch to continue merging till converges
    void MergeTillConverge(bool doit = true) { _merge_till_converge = doit; }

    /// A simple method to add a cluster
    void SetClusters(util::GeometryUtilities const& gser,
                     const std::vector<std::vector<util::PxHit>>& clusters);

    /// A simple method to add a cluster
    void SetClusters(const std::vector<cluster::ClusterParamsAlg>& clusters);

    /// A getter for input clusters
    const std::vector<cluster::ClusterParamsAlg>& GetInputClusters() const { return _in_clusters; }

    /// A setter for minimum # of hits ... passed onto ClusterParamsAlg
    void SetMinNHits(unsigned int n) { _min_nhits = n; }

    /// A method to execute the main action, to be called per event
    void Process(util::GeometryUtilities const& gser);

    /// A setter for an analysis output file
    void SetAnaFile(TFile* fout) { _fout = fout; }

  protected:
    /// Function to compute priority
    void ComputePriority(const std::vector<cluster::ClusterParamsAlg>& clusters);

    /// FMWK function called @ beginning of Process()
    virtual void EventBegin() {}

    /// FMWK function called @ beginning of iterative loop inside Process()
    virtual void IterationBegin() {}

    /// FMWK function called @ iterative loop inside Process()
    virtual bool IterationProcess(util::GeometryUtilities const& gser) = 0;

    /// FMWK function called @ end of iterative loop inside Process()
    virtual void IterationEnd() {}

    /// FMWK function called @ end of Process()
    virtual void EventEnd() {}

  protected:
    /// Timing verbosity flag
    bool _time_report;

    /// Minimum number of hits: the limit set for ClusterParamsAlg
    unsigned int _min_nhits;

    /// Debug mode switch
    CMMSGLevel_t _debug_mode;

    /// Input clusters
    std::vector<cluster::ClusterParamsAlg> _in_clusters;

    /// Priority algorithm
    ::cmtool::CPriorityAlgoBase* _priority_algo;

    /// Output analysis plot TFile
    TFile* _fout;

    /// Priority record
    std::multimap<float, size_t> _priority;

    /// Iteration loop switch
    bool _merge_till_converge;

    /// A holder for # of unique planes in the clusters, computed in ComputePriority() function
    std::set<UChar_t> _planes;
  };
}

#endif
/** @} */ // end of doxygen group
