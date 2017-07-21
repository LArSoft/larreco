////////////////////////////////////////////////////////////////////////
//
//
// TCAlg shower code
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGSHOWERS_H
#define TRAJCLUSTERALGSHOWERS_H


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
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace tca {

  void Find3DShowerEndPoints(TjStuff& tjs, const geo::TPCID& tpcid);
  void MakeShowers(TjStuff& tjs);
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP);
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid);
  void FindMatchingTjs(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void Match2DShowers(TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void FillPts(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  bool DefineShower(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void MakeShowerObsolete(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  bool AddTj(TjStuff& tjs, unsigned short TjID, const unsigned short& cotIndex, bool doUpdate, bool prt);
  bool RemoveTj(TjStuff& tjs, unsigned short TjID, const unsigned short& cotIndex, bool doUpdate, bool prt);
  bool FindChargeCenter(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void FindAngle(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void FillRotPos(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  bool DefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void FindExternalParent(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  bool WrongSplitTj(TjStuff& tjs, Trajectory& tj, const unsigned short& tjEnd, ShowerStruct& ss, bool prt);
  float ParentFOM(TjStuff& tjs, Trajectory& tj, const unsigned short& tjEnd, ShowerStruct& ss, bool prt);
  void DefineEnvelope(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void AddTjsInsideEnvelope(TjStuff& tjs, const unsigned short& cotIndex, bool prt, int mode);
  void RefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  bool AddLooseHits(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void FindStartChg(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  std::vector<float> StartChgVec(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void DumpShowerPts(TjStuff& tjs, const unsigned short& cotIndex);
  void CheckQuality(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList);
  void AddCloseTjsToList(TjStuff& tjs, const unsigned short& itj, std::vector<int> list);
  void MergeTjList(std::vector<std::vector<int>>& tjList);
  void MergeTjList2(TjStuff& tjs, std::vector<std::vector<int>>& tjList, bool prt);
  void MergeNearby2DShowers(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeOverlap(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  bool MergeShowersAndStore(TjStuff& tjs, unsigned short icotIndex, unsigned short jcotIndex, bool prt);
  bool MergeShowerTjsAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt);
  bool TransferTjHits(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  unsigned short GetCotsIndex(TjStuff& tjs, const unsigned short& ShowerTjID);
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss);
  unsigned short FarEnd(TjStuff& tjs, const Trajectory& tj, ShowerStruct& ss);
  unsigned short Create2DShower(TjStuff& tjs, const std::vector<int>& tjl);
  void FindNearbyTjs(TjStuff& tjs, const unsigned short& cotIndex, bool prt);

  void SaveTjInfo(TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList, unsigned int stageNum);
  void SaveTjInfo(TjStuff& tjs,  const CTP_t& inCTP, const unsigned short& cotIndex,unsigned int stageNum);
  void SaveTjInfoStuff(TjStuff& tjs,  const CTP_t& inCTP, Trajectory& tj,  unsigned int stageNum);
}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
