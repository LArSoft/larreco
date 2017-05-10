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
#include "larreco/Calorimetry/CalorimetryAlg.h"

// #include "TPrincipal.h"

namespace tca {

  void Find3DShowerEndPoints(TjStuff& tjs, const geo::TPCID& tpcid);
  void MakeShowers(TjStuff& tjs, const calo::CalorimetryAlg& fCaloAlg);
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP);
  void FillPts(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void DefineShower(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void RefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
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
  void AddTjsInsideEnvelope(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  bool AddLooseHits(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void FindStartChg(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  std::vector<float> StartChgVec(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void DumpShowerPts(TjStuff& tjs, const unsigned short& cotIndex);
  
  void AddMissedTjs(TjStuff& tjs, const CTP_t& inCTP, std::vector<unsigned short>& tjl);
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<unsigned short>>& tjList);
  void MergeOverlap(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  bool MergeShowersAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt);
  void TransferTjHits(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  unsigned short GetCotsIndex(TjStuff& tjs, const unsigned short& ShowerTjID);
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss);

}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
