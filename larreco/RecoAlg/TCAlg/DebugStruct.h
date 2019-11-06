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

#include <limits.h>

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {

  struct DebugStuff {
    int Cryostat {0};    ///< Select Cryostat
    int TPC {0};    ///< Select TPC
    int Plane {-1}; ///< Select plane
    CTP_t CTP {UINT_MAX}; ///< set to an invalid CTP
    int Wire {-1};  ///< Select hit Wire for debugging
    int Tick {-1};   ///< Select hit PeakTime for debugging (< 0 for vertex finding)
    unsigned int Hit {UINT_MAX};    ///< set to the hit index in evt.allHits if a Plane:Wire:Tick match is found
    int WorkID {0}; ///< Select the StartWorkID for debugging
    unsigned int MVI {UINT_MAX}; ///< MatchVec Index for detailed 3D matching
    unsigned short MVI_Iter {USHRT_MAX}; ///< MVI iteration - see FindPFParticles
    int Slice {-1};
  };
  extern DebugStuff debug;
} // namespace tca

#endif // ifndef TRAJCLUSTERALGDEBUGSTRUCT_H
