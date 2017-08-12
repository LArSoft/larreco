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
#include "larreco/RecoAlg/TCAlg/TCTree.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace tca {

  bool Find3DShowerEndPoints(TjStuff& tjs, MatchStruct& ms);
  void Finish3DShowers(TjStuff& tjs);
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid);
  void FindMatchingTjs(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool MergeShowerTjsAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt);
  bool TransferTjHits(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  unsigned short GetCotsIndex(TjStuff& tjs, unsigned short ShowerTjID);
  
  void Match2DShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void FillPts(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool DefineShower(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool AddTj(std::string inFcnLabel, TjStuff& tjs, int TjID, unsigned short cotIndex, bool doUpdate, bool prt);
  bool RemoveTj(std::string inFcnLabel, TjStuff& tjs, int TjID, unsigned short cotIndex, bool doUpdate, bool prt);
  bool FindChargeCenter(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FindAngle(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FillRotPos(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool AnalyzeRotPos(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool DefineShowerTj(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FindExternalParent(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void UpdateShowerWithParent(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool WrongSplitTj(TjStuff& tjs, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt);
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, bool prt);
  void DefineEnvelope(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void RefineShowerTj(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool AddLooseHits(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FindStartChg(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  std::vector<float> StartChgVec(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void DumpShowerPts(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex);
  void CheckQuality(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  
  void TagShowerTjs(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList);
  void AddCloseTjsToList(std::string inFcnLabel, TjStuff& tjs, unsigned short itj, std::vector<int> list);
  void MergeTjList(std::vector<std::vector<int>>& tjList);
  void MergeTjList2(std::string inFcnLabel, TjStuff& tjs, std::vector<std::vector<int>>& tjList, bool prt);
  void MergeNearby2DShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeOverlap(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeShowerChain(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  int MergeShowers(std::string inFcnLabel, TjStuff& tjs, std::vector<int> showerIDs, bool prt);
  bool MergeShowersAndStore(std::string inFcnLabel, TjStuff& tjs, unsigned short icotIndex, unsigned short jcotIndex, bool prt);
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss);
  unsigned short FarEnd(TjStuff& tjs, const Trajectory& tj, ShowerStruct& ss);
  unsigned short Create2DShower(TjStuff& tjs, const std::vector<int>& tjl);
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, unsigned short cotIndex, bool prt);
  void Print2DShowers(std::string someText, const TjStuff& tjs, CTP_t inCTP, bool printKilledShowers);

}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
