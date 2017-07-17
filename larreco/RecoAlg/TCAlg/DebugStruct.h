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
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "lardataobj/RecoBase/Hit.h"

namespace tca {
  
  struct DebugStuff {
    int Cryostat {0};    ///< Select Cryostat
    int TPC {0};    ///< Select TPC
    int Plane {-1}; ///< Select plane
    CTP_t CTP {UINT_MAX}; ///< set to an invalid CTP
    int Wire {-1};  ///< Select hit Wire for debugging
    int Tick {-1};   ///< Select hit PeakTime for debugging (< 0 for vertex finding)
    unsigned int Hit {UINT_MAX};    ///< set to the hit index in fHits if a Plane:Wire:Tick match is found
    int WorkID {0}; ///< Select the StartWorkID for debugging
  };
  extern DebugStuff debug;
} // namespace tca

#endif // ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
