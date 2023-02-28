/**
 * \file CBAlgoMergeAll.h
 *
 * \ingroup CMTool
 *
 * \brief Class def header for a class CBAlgoMergeAll
 *
 * @author david caratelli
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CBALGOMERGEALL_H
#define RECOTOOL_CBALGOMERGEALL_H

#include "larreco/RecoAlg/CMTool/CMToolBase/CBoolAlgoBase.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"

namespace cmtool {
  /**
     \class CBAlgoMergeAll
     Merges all clusters: maybe useful to test how well a cluster-separating
     algorithm has performed
  */
  class CBAlgoMergeAll : public CBoolAlgoBase {

  public:
    /// Default constructor
    CBAlgoMergeAll();

    /// Default destructor
    virtual ~CBAlgoMergeAll(){};

    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ::cluster::ClusterParamsAlg& cluster1,
                      const ::cluster::ClusterParamsAlg& cluster2);

    /// Function to reset the algorithm instance ... maybe implemented via child class
    virtual void Reset() {}

  protected:
  };
}

#endif
/** @} */ // end of doxygen group
