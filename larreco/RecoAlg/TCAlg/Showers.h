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
//#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace tca {
  void MakeShowers(TjStuff& tjs, const std::vector<float>& fShowerTag, const art::ServiceHandle<geo::Geometry>& geom, const detinfo::DetectorProperties*& detprop, const calo::CalorimetryAlg& fCaloAlg);
//  void MakeShowers(TjStuff& tjs, const std::vector<float>& fShowerTag, const art::ServiceHandle<geo::Geometry>& geom, const calo::CalorimetryAlg& fCaloAlg);
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag);
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag, std::vector<std::vector<unsigned short>>& tjList);
  void FindShowerCenter(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void FindShowerParent(TjStuff& tjs, const unsigned short& showerIndex, const std::vector<float>& fShowerTag, bool prt);
  void FindFirstTPAng(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void DefineShowerTj(TjStuff& tjs, const unsigned short& cotIndex, bool prt);
  void DefineShowerEnvelope(TjStuff& tjs, const unsigned short& cotIndex, const std::vector<float>& fShowerTag, bool prt);
  void MergeShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag, bool prt);
  void CollectHits(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  float ShowerEnergy(TjStuff& tjs, const ShowerStruct& ss);
  void SpacePtDir(TjStuff& tjs, const art::ServiceHandle<geo::Geometry>& geom, const detinfo::DetectorProperties*& detprop, TrajPoint itp, TrajPoint jtp, TVector3& dir, TVector3& dirErr);
}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
