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
  
  void MakeJunkVertices(TjStuff& tjs, const CTP_t& inCTP);
  void Find2DVertices(TjStuff& tjs, const CTP_t& inCTP);
  void FindNeutralVertices(TjStuff& tjs, const geo::TPCID& tpcid);
  void MakeJunkTjVertices(TjStuff& tjs, const CTP_t& inCTP);
  void ChkVxTjs(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  bool MergeWithVertex(TjStuff& tjs, VtxStore& vx2, unsigned short existingVxID, bool prt);
  void SplitTrajCrossingVertices(TjStuff& tjs, CTP_t inCTP);
  void FindHammerVertices(TjStuff& tjs, const CTP_t& inCTP);
  void FindHammerVertices2(TjStuff& tjs, const CTP_t& inCTP);
  void Find3DVertices(TjStuff& tjs, const geo::TPCID& tpcid);
  void Match3DVtxTjs(TjStuff& tjsconst, const geo::TPCID& tpcid, bool prt);
  void CompleteIncomplete3DVertices(TjStuff& tjs, const geo::TPCID& tpcid);
  bool RefineVtxPosition(TjStuff& tjs, const Trajectory& tj, unsigned short& nearPt, short nPtsToChk, bool prt);
  void CompleteIncomplete3DVerticesInGaps(TjStuff& tjs, const geo::TPCID& tpcid);
  // Improve hit assignments near vertex 
  void VtxHitsSwap(TjStuff& tjs, const CTP_t inCTP);

  unsigned short TPNearVertex(TjStuff& tjs, const TrajPoint& tp);
  bool AttachPFPToVertex(TjStuff& tjs, PFPStruct& pfp, unsigned short end, unsigned short vx3ID, bool prt);
  bool AttachAnyTrajToVertex(TjStuff& tjs, unsigned short iv, bool prt);
  bool AttachTrajToVertex(TjStuff& tjs, Trajectory& tj, VtxStore& vx, bool prt);
  float TrajPointVertexPull(TjStuff& tjs, const TrajPoint& tp, const VtxStore& vx);
  float VertexVertexPull(TjStuff& tjs, const Vtx3Store& vx1, const Vtx3Store& vx2);
  float VertexVertexPull(TjStuff& tjs, const VtxStore& vx1, const VtxStore& vx2);
  bool FitVertex(TjStuff& tjs, VtxStore& vx, bool prt);
  bool FitVertex(TjStuff& tjs, VtxStore& vx, std::vector<TrajPoint> vxTp, bool prt);
  bool StoreVertex(TjStuff& tjs, VtxStore& vx);
  bool ChkVtxAssociations(TjStuff& tjs, const CTP_t& inCTP);
  void ScoreVertices(TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void KillPoorVertices(TjStuff& tjs, const geo::TPCID& tpcid);
  void SetVx2Score(TjStuff& tjs, bool prt);
  void SetVx2Score(TjStuff& tjs, VtxStore& vx2, bool prt);
  void SetVx3Score(TjStuff& tjs, Vtx3Store& vx3, bool prt);
  unsigned short Vx3Topo(TjStuff& tjs, Vtx3Store& vx3);
  void SetHighScoreBits(TjStuff& tjs, Vtx3Store& vx3);
  bool MakeVertexObsolete(TjStuff& tjs, VtxStore& vx2, bool forceKill);
  bool MakeVertexObsolete(TjStuff& tjs, Vtx3Store& vx3);
  std::vector<int> GetVtxTjIDs(const TjStuff& tjs, const VtxStore& vx2);
  std::vector<int> GetVtxTjIDs(const TjStuff& tjs, const Vtx3Store& vx3, float& score);
  std::vector<unsigned short> GetPFPVertices(const TjStuff& tjs, const PFPStruct& pfp);
  void PosInPlane(const TjStuff& tjs, const Vtx3Store& vx3, unsigned short plane, Point2_t& pos);
  unsigned short IsCloseToVertex(TjStuff& tjs, VtxStore& vx);
  unsigned short IsCloseToVertex(TjStuff& tjs, Vtx3Store& vx3);
  //    void Refine2DVertices();
} // namespace

#endif // ifndef TRAJCLUSTERALGVERTEX_H
