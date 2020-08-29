/**
 * \file CFAlgoZOverlap.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CFAlgoZOverlap
 *
 * @author ah673
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CFALGOZOVERLAP_H
#define RECOTOOL_CFALGOZOVERLAP_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CFloatAlgoBase.h"

namespace cmtool {
  /**
     \class CFAlgoZOverlap
     User implementation for CFloatAlgoBase class
     doxygen documentation!
  */
  class CFAlgoZOverlap : public CFloatAlgoBase {

  public:
    /// Default constructor
    CFAlgoZOverlap();

    //
    // Author should be aware of 3 functions at least: Float, Report,
    // and Reset. More possibly-useful functions can be found in the later
    // part but commented out. All of these functions are virtual and defined
    // in the base class.

    /**
       Core function: given a set of CPANs, return a float which indicates
       the compatibility the cluster combination.
    */
    float Float(util::GeometryUtilities const&,
                const std::vector<const cluster::ClusterParamsAlg*>& clusters) override;

    /**
       Optional function: called after each iterative approach if a manager class is
       run with verbosity level <= kPerIteration. Maybe useful for debugging.
    */
    void Report() override;

    /// Function to reset the algorithm instance, called together with manager's Reset()
    void Reset() override;

    /**
       Optional function: called at the beginning of 1st iteration. This is called per event.
     */
    //void EventBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of event ... after the last merging iteration is over.
     */
    //void EventEnd();

    /**
       Optional function: called at the beggining of each iterative loop.
       This provides all clusters' information in case the algorithm need them. Note this
       is called per iteration which may be more than once per event.
     */
    //void IterationBegin(const std::vector<cluster::ClusterParamsAlg> &clusters);

    /**
       Optional function: called at the end of each iterative loop.
     */
    //void IterationEnd();

  protected:
    float _wire_ratio_cut;
  };
}
#endif
/** @} */ // end of doxygen group
