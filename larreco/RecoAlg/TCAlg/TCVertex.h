////////////////////////////////////////////////////////////////////////
//
//
// TCAlg vertex code
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGVERTEX_H
#define TRAJCLUSTERALGVERTEX_H


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

namespace tca {
  
  void Find2DVertices(TjStuff& tjs, const CTP_t& inCTP);
  bool MergeWithNearbyVertex(TjStuff& tjs, VtxStore& vx2, unsigned short it1, unsigned short end1, unsigned short it2, unsigned short end2);
  void FindHammerVertices(TjStuff& tjs, const CTP_t& inCTP);
  void FindHammerVertices2(TjStuff& tjs, const CTP_t& inCTP);
  void Find3DVertices(TjStuff& tjs, const geo::TPCID& tpcid);
   void Match3DVtxTjs(TjStuff& tjsconst, const geo::TPCID& tpcid, bool prt);
  void CompleteIncomplete3DVertices(TjStuff& tjs, const geo::TPCID& tpcid);
  void CompleteIncomplete3DVerticesInGaps(TjStuff& tjs, const geo::TPCID& tpcid);
  // Improve hit assignments near vertex 
  void VtxHitsSwap(TjStuff& tjs, const CTP_t inCTP);

  unsigned short TPNearVertex(TjStuff& tjs, const TrajPoint& tp);
  bool AttachAnyTrajToVertex(TjStuff& tjs, unsigned short iv, bool prt);
  bool AttachTrajToAnyVertex(TjStuff& tjs, unsigned short itj, bool prt);
  bool AttachTrajToVertex(TjStuff& tjs, Trajectory& tj, VtxStore& vx, bool prt);
  float TrajPointVertexPull(TjStuff& tjs, const TrajPoint& tp, const VtxStore& vx);
  float VertexVertexPull(TjStuff& tjs, const VtxStore& vx1, const VtxStore& vx2);
  bool FitVertex(TjStuff& tjs, VtxStore& vx, bool prt);
  bool StoreVertex(TjStuff& tjs, VtxStore& vx);
  void ChkVtxAssociations(TjStuff& tjs, const CTP_t& inCTP);
  bool TjHasNiceVtx(const TjStuff& tjs, const Trajectory& tj, float scoreCut);
  void SetVtxScore(TjStuff& tjs, VtxStore& vx2, bool prt);
  void SetVtxScore(TjStuff& tjs, Vtx3Store& vx3, bool prt);
  void KillPoorVertices(TjStuff& tjs);
  bool MakeVertexObsolete(TjStuff& tjs, unsigned short ivx, bool forceKill);
  std::vector<int> GetVtxTjIDs(const TjStuff& tjs, const VtxStore& vx2);
  std::vector<int> GetVtxTjIDs(const TjStuff& tjs, const Vtx3Store& vx3, float& score);
  //    void Refine2DVertices();
} // namespace

#endif // ifndef TRAJCLUSTERALGVERTEX_H
