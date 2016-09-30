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
#include "lardataobj/RecoBase/Hit.h"

namespace tca {
  
  struct DebugStuff {
    int Plane {-1}; ///< Select plane
    int Wire {-1};  ///< Select hit Wire for debugging
    int Tick {-1};   ///< Select hit PeakTime for debugging (< 0 for vertex finding)
    unsigned int Hit{UINT_MAX};    ///< set to the hit index in fHits if a Plane:Wire:Tick match is found
    short WorkID {0}; ///< Select the StartWorkID for debugging
  };
  extern DebugStuff debug;
} // namespace tca

#endif // ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
