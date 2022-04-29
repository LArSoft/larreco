#include "larreco/RecoAlg/CMTool/CMTAlgPriority/CPAlgoPolyArea.h"

namespace cmtool {

  //----------------------------------------------------------
  CPAlgoPolyArea::CPAlgoPolyArea() : CPriorityAlgoBase()
  //----------------------------------------------------------
  {
    _area_cut = 0;
  }

  //------------------------------------------------------------------------------
  float CPAlgoPolyArea::Priority(const ::cluster::ClusterParamsAlg &cluster)
  //------------------------------------------------------------------------------
  {

    auto area = cluster.GetParams().PolyObject.Area();

    return ( area < _area_cut ? -1 : area);
  }

}
