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
  void KillVerticesInShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void Finish3DShowers(TjStuff& tjs);
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid);
  void FindInShowerPFPs(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, std::vector<std::vector<unsigned short>>& plists);
  void CheckInShowerProb(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool MergeShowerTjsAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt);
  bool TransferTjHits(TjStuff& tjs, bool prt);
  int GetCotID(TjStuff& tjs, int ShowerTjID);
  
  void Match2DShowers(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void UpdateShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool FindMissingShowers1(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool FindMissingShowers2(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void ReconcileParents(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, bool prt);
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, int kcotID, bool prt);
  float Match3DFOM(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  void FillPts(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool DefinePFPShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool UpdatePFPShower(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  bool DefineShower(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void MakeShowerObsolete(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool AddTj(std::string inFcnLabel, TjStuff& tjs, int TjID, int cotID, bool doUpdate, bool prt);
  bool RemoveTj(std::string inFcnLabel, TjStuff& tjs, int TjID, int cotID, bool doUpdate, bool prt);
  bool FindChargeCenter(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void FindAngle(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void FillRotPos(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool AnalyzeRotPos(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool DefineShowerTj(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void ReverseShower(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void FindParent(std::string inFcnLabel, TjStuff& tjs, ShowerStruct3D& ss3, bool prt);
  void FindExternalParent(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool UpdateShowerWithParent(std::string inFcnLabel, TjStuff& tjs, int cotID, unsigned short newParent, float newParentFOM, bool prt);
  bool WrongSplitTj(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt);
  bool IsInShower(const TjStuff& tjs, const std::vector<int> TjIDs);
  float InShowerProb(std::string inFcnLabel, const TjStuff& tjs, const ShowerStruct3D& ss3, const PFPStruct& pfp);
  void ShowerParams(double showerEnergy, double& shMaxAlong, double& shE95Along);
  double ShowerParamTransRMS(double showerEnergy, double along);
  double InShowerProbLong(double showerEnergy, double along);
  double InShowerProbTrans(double showerEnergy, double along, double trans);
  double InShowerProb(double showerEnergy, double along, double trans);
  float ParentFOM(std::string inFcnLabel, TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, float& tp1Sep, float& vx3Score, bool prt);
  void DefineEnvelope(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool AddTjsInsideEnvelope(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void RefineShowerTj(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  bool AddLooseHits(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void FindStartChg(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  std::vector<float> StartChgVec(TjStuff& tjs, int cotID, bool prt);
  void DumpShowerPts(std::string inFcnLabel, TjStuff& tjs, int cotID);
  void CheckQuality(std::string inFcnLabel, TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  
  void TagShowerLike(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList, bool applyMinTjCuts);
  void FindNearbyTjs(std::string inFcnLabel, TjStuff& tjs, int cotID, bool prt);
  void AddCloseTjsToList(std::string inFcnLabel, TjStuff& tjs, unsigned short itj, std::vector<int> list);
  void MergeTjList(std::vector<std::vector<int>>& tjList);
  void MergeTjList2(std::string inFcnLabel, TjStuff& tjs, std::vector<std::vector<int>>& tjList, bool prt);
  void MergeNearby2DShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeOverlap(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeShowerChain(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(std::string inFcnLabel, TjStuff& tjs, const CTP_t& inCTP, bool prt);
  int MergeShowers(std::string inFcnLabel, TjStuff& tjs, std::vector<int> showerIDs, bool prt);
  bool MergeShowersAndStore(std::string inFcnLabel, TjStuff& tjs, int icotID, int jcotID, bool prt);
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss);
  float ChgToMeV(float chg);
  unsigned short FarEnd(TjStuff& tjs, const Trajectory& tj, Point2_t& pos);
  unsigned short FarEnd(TjStuff& tjs, const PFPStruct& pfp, Point3_t& pos);
  PFPStruct CreateFakePFP(const TjStuff& tjs, const ShowerStruct3D& ss3);
  ShowerStruct3D CreateSS3(TjStuff& tjs, const geo::TPCID& tpcid);
  unsigned short Create2DShower(TjStuff& tjs, const std::vector<int>& tjl);
  void PrintShowers(std::string someText, TjStuff& tjs);
  void Print2DShowers(std::string someText, const TjStuff& tjs, CTP_t inCTP, bool printKilledShowers);

}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
