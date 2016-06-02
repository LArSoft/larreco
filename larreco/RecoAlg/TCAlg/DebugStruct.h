////////////////////////////////////////////////////////////////////////
//
//
// TCAlg debug struct
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
#define TRAJCLUSTERALGDEBUGSTRUCT_H


// C/C++ standard libraries
#include <array>
#include <vector>
#include <bitset>

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardata/RecoBase/Hit.h"

namespace tca {
  
  struct DebugStuff {
    int Plane {-1}; ///< Select plane
    int Wire {-1};  ///< Select hit Wire for debugging
    int Tick {-1};   ///< Select hit PeakTime for debugging (< 0 for vertex finding)
  };
  extern DebugStuff Debug;
} // namespace tca

#endif // ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
