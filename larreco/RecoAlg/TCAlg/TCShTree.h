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

  void SaveTjInfo(TCSlice& slc, std::vector<std::vector<int>>& tjList, std::string stageName);
  void SaveTjInfo(TCSlice& slc, const ShowerStruct& ss, std::string stageName);
  void SaveTjInfoStuff(TCSlice& slc, Trajectory& tj, int stageNum, std::string stageName);
  void SaveAllCots(TCSlice& slc, const CTP_t& inCTP, std::string someText);
  void SaveAllCots(TCSlice& slc, std::string someText);
  int GetStageNum(ShowerTreeVars& stv, std::string stageName);
  void ClearShowerTree(ShowerTreeVars& stv);

} // namespace tca

#endif
