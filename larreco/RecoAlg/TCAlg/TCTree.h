////////////////////////////////////////////////////////////////////////                            
//                                                                                                  
//                                                                                                  
// TCAlg shower tree code
//                                                                           
// Rory Fitzpatrick
//
///////////////////////////////////////////////////////////////////////   
#ifndef TRAJCLUSTERALGTREE_H
#define TRAJCLUSTERALGTREE_H

// C/C++ standard libraries                                                                         
#include <array>
#include <vector>
#include <bitset>
#include <utility> // std::pair<>                                                                   
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries                                                                                
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/TCShower.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace tca {

  void SaveTjInfo(TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList, unsigned int stageNum);
  void SaveTjInfo(TjStuff& tjs,  const CTP_t& inCTP, const unsigned short& cotIndex, unsigned int stageNum);
  void SaveTjInfoStuff(TjStuff& tjs,  const CTP_t& inCTP, Trajectory& tj,  unsigned int stageNum);

  void ClearShowerTree(ShowerTreeVars& stv);

} // namespace tca

#endif
