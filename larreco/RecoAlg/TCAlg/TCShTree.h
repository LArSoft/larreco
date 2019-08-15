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
#include <string>
#include <vector>

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"

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
