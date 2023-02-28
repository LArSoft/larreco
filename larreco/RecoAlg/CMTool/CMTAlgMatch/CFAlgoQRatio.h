/**
 * \file CFAlgoQRatio.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CFAlgoQRatio
 *
 * @author kazuhiro
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CFALGOQRATIO_H
#define RECOTOOL_CFALGOQRATIO_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CFloatAlgoBase.h"

namespace cmtool {
  /**
     \class CFAlgoQRatio
     User implementation for CFloatAlgoBase class
     This algorithm compares charge ratio of clusters to find a match
  */
  class CFAlgoQRatio : public CFloatAlgoBase {

  public:
    /// Default constructor
    CFAlgoQRatio();

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

    /// Setter for the minimum value for charge ratio (below this value Float() returns -1)
    void SetQRatioCut(float cut) { _qratio_cut = cut; }

  protected:
    float _qratio_cut;
  };
}
#endif
/** @} */ // end of doxygen group
