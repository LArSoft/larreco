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
#include <vector>
#include <string>

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {

  void ConfigureMVA(TCConfig& tcc, std::string fMVAShowerParentWeights);
  bool FindShowerStart(TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  void KillVerticesInShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  void Finish3DShowers(TCSlice& slc);
  bool FindShowers3D(TCSlice& slc);
  bool Reconcile3D(std::string inFcnLabel, TCSlice& slc, bool parentSearchDone, bool prt);
  bool Reconcile3D(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  bool MergeShowerTjsAndStore(TCSlice& slc, unsigned short istj, unsigned short jstj, bool prt);
  bool TransferTjHits(TCSlice& slc, bool prt);
  int GetCotID(TCSlice& slc, int ShowerTjID);

  bool CompleteIncompleteShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  void Match2DShowers(std::string inFcnLabel, TCSlice& slc, bool prt);
  bool UpdateShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  bool UpdateShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  float Match3DFOM(std::string inFcnLabel, TCSlice& slc, int icotID, int jcotID, bool prt);
  float Match3DFOM(std::string inFcnLabel, TCSlice& slc, int icotID, int jcotID, int kcotID, bool prt);
  float Match3DFOM(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  void MakeShowerObsolete(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  void MakeShowerObsolete(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  bool DontCluster(TCSlice& slc, const std::vector<int>& tjlist1, const std::vector<int>& tjlist2);
  void DefineDontCluster(TCSlice& slc, bool prt);
  bool RemovePFP(std::string inFcnLabel, TCSlice& slc, PFPStruct& pfp, ShowerStruct3D& ss3, bool doUpdate, bool prt);
  bool AddPFP(std::string inFcnLabel, TCSlice& slc, int pID, ShowerStruct3D& ss3, bool doUpdate, bool prt);
  bool RemovePFP(std::string inFcnLabel, TCSlice& slc, int pID, ShowerStruct3D& ss3, bool doUpdate, bool prt);
  bool AddTj(std::string inFcnLabel, TCSlice& slc, int TjID, ShowerStruct& ss, bool doUpdate, bool prt);
  bool RemoveTj(std::string inFcnLabel, TCSlice& slc, int TjID, ShowerStruct& ss, bool doUpdate, bool prt);
  bool AnalyzeRotPos(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  void ReverseShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  void ReverseShower(std::string inFcnLabel, TCSlice& slc, int cotID, bool prt);
  bool FindParent(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3, bool prt);
  bool SetParent(std::string inFcnLabel, TCSlice& slc, PFPStruct& pfp, ShowerStruct3D& ss3, bool prt);
  bool WrongSplitTj(std::string inFcnLabel, TCSlice& slc, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt);
  bool IsShowerLike(TCSlice& slc, const std::vector<int> TjIDs);
  float InShowerProb(TCSlice& slc, const ShowerStruct3D& ss3, const PFPStruct& pfp);
  float InShowerProb(TCSlice& slc, const ShowerStruct& ss, const Trajectory& tj);
  void ShowerParams(double showerEnergy, double& shMaxAlong, double& shE95Along);
  double ShowerParamTransRMS(double showerEnergy, double along);
  double InShowerProbLong(double showerEnergy, double along);
  double InShowerProbTrans(double showerEnergy, double along, double trans);
  double InShowerProb(double showerEnergy, double along, double trans);
  float ParentFOM(std::string inFcnLabel, TCSlice& slc, PFPStruct& pfp, unsigned short pend, ShowerStruct3D& ss3, bool prt);
  float ParentFOM(std::string inFcnLabel, TCSlice& slc, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, float& tp1Sep, float& vx3Score, bool prt);
  void DefineEnvelope(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  bool AddTjsInsideEnvelope(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  bool AddLooseHits(std::string inFcnLabel, TCSlice& slc, int cotID, bool prt);
  void FindStartChg(std::string inFcnLabel, TCSlice& slc, int cotID, bool prt);
  std::vector<float> StartChgVec(TCSlice& slc, int cotID, bool prt);
  void DumpShowerPts(std::string inFcnLabel, TCSlice& slc, int cotID);

  void FindCots(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP, std::vector<std::vector<int>>& tjLists, bool prt);
  void TagShowerLike(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP);
  void FindNearbyTjs(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss, bool prt);
  void AddCloseTjsToList(std::string inFcnLabel, TCSlice& slc, unsigned short itj, std::vector<int> list);
  void MergeTjList(std::vector<std::vector<int>>& tjList);
  void MergeTjList2(std::string inFcnLabel, TCSlice& slc, std::vector<std::vector<int>>& tjList, bool prt);
  void MergeNearby2DShowers(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP, bool prt);
  void MergeOverlap(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP, bool prt);
  void MergeShowerChain(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP, bool prt);
  void MergeSubShowersTj(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(std::string inFcnLabel, TCSlice& slc, const CTP_t& inCTP, bool prt);
  int MergeShowers(std::string inFcnLabel, TCSlice& slc, std::vector<int> showerIDs, bool prt);
  bool MergeShowersAndStore(std::string inFcnLabel, TCSlice& slc, int icotID, int jcotID, bool prt);
  double ShowerEnergy(const ShowerStruct3D& ss3);
  float ShowerEnergy(TCSlice& slc, const std::vector<int> tjIDs);
  float ChgToMeV(float chg);
 bool StoreShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct3D& ss3);
  bool StoreShower(std::string inFcnLabel, TCSlice& slc, ShowerStruct& ss);
  ShowerStruct3D CreateSS3(TCSlice& slc);
  ShowerStruct CreateSS(TCSlice& slc, const std::vector<int>& tjl);
  bool ChkAssns(std::string inFcnLabel, TCSlice& slc);
  void PrintShowers(std::string someText, TCSlice& slc);
  void Print2DShowers(std::string someText, TCSlice& slc, CTP_t inCTP, bool printKilledShowers);
  void PrintShower(std::string someText, TCSlice& slc, const ShowerStruct& ss, bool printHeader, bool printExtras);

}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
