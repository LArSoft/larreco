/**
 * \file CPAlgoPolyArea.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CPAlgoPolyArea
 *
 * @author kazuhiro
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CPALGOPOLYAREA_H
#define RECOTOOL_CPALGOPOLYAREA_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CPriorityAlgoBase.h"

namespace cmtool {
  /**
     \class CPAlgoPolyArea
     Simple algorithm to determine priority based on area of 2D polygon.
     If area < set cut value by a user, returns -1.
  */
  class CPAlgoPolyArea : public CPriorityAlgoBase {

  public:
    /// Default constructor
    CPAlgoPolyArea();

    /**
       Core function: given the CPAN input, return a float which indicates
       the user-defined priority for analysis.
    */
    virtual float Priority(const ::cluster::ClusterParamsAlg& cluster);

    /// Setter for minimum area
    void SetMinArea(double area) { _area_cut = area; }

  private:
    double _area_cut;
  };
}
#endif
/** @} */ // end of doxygen group
