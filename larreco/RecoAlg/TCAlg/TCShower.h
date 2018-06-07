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
#include "larreco/RecoAlg/TCAlg/TCShTree.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace tca {

  bool FindShowerStart(TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  void KillVerticesInShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  void Finish3DShowers(TjStuff& tjs);
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid);
  bool Reconcile3D(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool parentSearchDone, bool prt);
  bool Reconcile3D(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
//  void FindInShowerPFPs(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, std::vector<std::vector<int>>& plists);
  bool MergeShowerTjsAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt);
  bool TransferTjHits(TjStuff& tjs, bool prt);
  int GetCotID(TjStuff& tjs, int ShowerTjID);
  
  void CompleteIncompleteShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  void Match2DShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  bool UpdateShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  bool UpdateShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, bool prt);
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, int kcotID, bool prt);
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  bool DontCluster(const TjStuff& tjs, const std::vector<int>& tjlist1, const std::vector<int>& tjlist2);
  void DefineDontCluster(TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  bool RemovePFP(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, ShowerStruct3D& ss3, bool doUpdate, bool prt);
  bool AddPFP(std::string inFcnLabel, TjStuff& tjs, int pID, ShowerStruct3D& ss3, bool doUpdate, bool prt);
  bool RemovePFP(std::string inFcnLabel, TjStuff& tjs, int pID, ShowerStruct3D& ss3, bool doUpdate, bool prt);
  bool AddTj(std::string inFcnLabel, TjStuff& tjs, int TjID, ShowerStruct& ss, bool doUpdate, bool prt);
  bool RemoveTj(std::string inFcnLabel, TjStuff& tjs, int TjID, ShowerStruct& ss, bool doUpdate, bool prt);
  bool AnalyzeRotPos(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  bool DefineShowerTj(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool FindNeutrinoParent(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool FindParent(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool SetParent(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, ShowerStruct3D& ss3, bool prt);
  bool WrongSplitTj(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt);
  bool IsShowerLike(const TjStuff& tjs, const std::vector<int> TjIDs);
  float InShowerProb(const TjStuff& tjs, const ShowerStruct3D& ss3, const PFPStruct& pfp);
  float InShowerProb(const TjStuff& tjs, const ShowerStruct& ss, const Trajectory& tj);
  void ShowerParams(double showerEnergy, double& shMaxAlong, double& shE95Along);
  double ShowerParamTransRMS(double showerEnergy, double along);
  double InShowerProbLong(double showerEnergy, double along);
  double InShowerProbTrans(double showerEnergy, double along, double trans);
  double InShowerProb(double showerEnergy, double along, double trans);
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, PFPStruct& pfp, unsigned short pend, ShowerStruct3D& ss3, bool prt);
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, float& tp1Sep, float& vx3Score, bool prt);
  void DefineEnvelope(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  bool AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  bool AddLooseHits(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void FindStartChg(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  std::vector<float> StartChgVec(TjStuff& tjs, int cotID, bool prt);
  void DumpShowerPts(std::string inFcnLabel, TjStuff& tjs, int cotID);
  void CheckQuality(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  
  void FindCots(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjLists, bool prt);
  void TagShowerLike(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP);
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss, bool prt);
  void AddCloseTjsToList(std::string inFcnLabel, TjStuff& tjs, unsigned short itj, std::vector<int> list);
  void MergeTjList(std::vector<std::vector<int>>& tjList);
  void MergeTjList2(std::string inFcnLabel, TjStuff& tjs, std::vector<std::vector<int>>& tjList, bool prt);
  void MergeNearby2DShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeOverlap(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeShowerChain(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowersTj(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  int MergeShowers(std::string inFcnLabel, TjStuff& tjs, std::vector<int> showerIDs, bool prt);
  bool MergeShowersAndStore(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, bool prt);
  double ShowerEnergy(const ShowerStruct3D& ss3);
  float ShowerEnergy(const TjStuff& tjs, const std::vector<int> tjIDs);
  float ChgToMeV(float chg);
  PFPStruct CreateFakePFP(const TjStuff& tjs, const ShowerStruct3D& ss3);
  bool StoreShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3);
  bool StoreShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct& ss);
  ShowerStruct3D CreateSS3(TjStuff& tjs, const geo::TPCID& tpcid);
  ShowerStruct CreateSS(TjStuff& tjs, const std::vector<int>& tjl);
  bool ChkAssns(std::string inFcnLabel, TjStuff& tjs);
  void PrintShowers(std::string someText, TjStuff& tjs);
  void Print2DShowers(std::string someText, const TjStuff& tjs, CTP_t inCTP, bool printKilledShowers);
  void PrintShower(std::string someText, const TjStuff& tjs, const ShowerStruct& ss, bool printHeader, bool printExtras);

}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
