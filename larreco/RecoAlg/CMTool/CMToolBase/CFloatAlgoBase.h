/**
 * \file CFloatAlgoBase.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CFloatAlgoBase
 *
 * @author kazuhiro
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CFLOATALGOBASE_H
#define RECOTOOL_CFLOATALGOBASE_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CMAlgoBase.h"
namespace util {
  class GeometryUtilities;
}

#include <vector>

namespace cmtool {

  /**
     \class CFloatAlgoBase
     An abstract base class for mathcing algorithm. Having this base class helps
     to have a better overall design of various merging for iterative approach.
     The algorithms are run through CMergeManager.
  */
  class CFloatAlgoBase : public CMAlgoBase {
  public:
    /**
       Core function: given the CPAN input, return whether a cluster should be
       merged or not.
    */
    virtual float Float(util::GeometryUtilities const&,
                        const std::vector<const cluster::ClusterParamsAlg*>& clusters)
    {
      if (clusters.size())
        return 1;
      else
        return -1;
    }
  };

}

#endif
/** @} */ // end of doxygen group
