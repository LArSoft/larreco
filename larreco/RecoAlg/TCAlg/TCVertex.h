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
#include <string>
#include <vector>

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
namespace detinfo {
  class DetectorPropertiesData;
}

namespace tca {

  extern TCEvent evt;
  extern TCConfig tcc;
  // vector of hits, tjs, etc in each slice
  extern std::vector<TCSlice> slices;

  void MakeJunkVertices(TCSlice& slc, const CTP_t& inCTP);
  void Find2DVertices(detinfo::DetectorPropertiesData const& detProp,
                      TCSlice& slc,
                      const CTP_t& inCTP,
                      unsigned short pass);
  void MakeJunkTjVertices(TCSlice& slc, const CTP_t& inCTP);
  bool MergeWithVertex(TCSlice& slc, VtxStore& vx2, unsigned short existingVxID);
  void FindHammerVertices(TCSlice& slc, const CTP_t& inCTP);
  void FindHammerVertices2(TCSlice& slc, const CTP_t& inCTP);
  void SplitTrajCrossingVertices(TCSlice& slc, CTP_t inCTP);
  void Reconcile2Vs(TCSlice& slc);
  bool Reconcile2VTs(TCSlice& slc, std::vector<int>& vx2cls, bool prt);
  void Find3DVertices(detinfo::DetectorPropertiesData const& detProp, TCSlice& slc);
  void CompleteIncomplete3DVertices(detinfo::DetectorPropertiesData const& detProp, TCSlice& slc);
  void CompleteIncomplete3DVerticesInGaps(detinfo::DetectorPropertiesData const& detProp,
                                          TCSlice& slc);
  bool RefineVtxPosition(const Trajectory& tj, unsigned short& nearPt, short nPtsToChk, bool prt);
  unsigned short TPNearVertex(const TCSlice& slc, const TrajPoint& tp);
  bool AttachToAnyVertex(TCSlice& slc, PFPStruct& pfp, float maxSep, bool prt);
  bool AttachAnyVertexToTraj(TCSlice& slc, int tjID, bool prt);
  bool AttachAnyTrajToVertex(TCSlice& slc, unsigned short iv, bool prt);
  bool AttachTrajToVertex(TCSlice& slc, Trajectory& tj, VtxStore& vx, bool prt);
  float TrajPointVertexPull(const TrajPoint& tp, const VtxStore& vx);
  float VertexVertexPull(const Vtx3Store& vx1, const Vtx3Store& vx2);
  float VertexVertexPull(const VtxStore& vx1, const VtxStore& vx2);
  bool FitVertex(TCSlice& slc, VtxStore& vx, bool prt);
  bool FitVertex(TCSlice& slc, VtxStore& vx, std::vector<TrajPoint>& vxTp, bool prt);
  bool StoreVertex(TCSlice& slc, VtxStore& vx);
  bool ChkVtxAssociations(TCSlice& slc, const CTP_t& inCTP);
  void ScoreVertices(TCSlice& slc);
  void KillPoorVertices(TCSlice& slc);
  void SetVx2Score(TCSlice& slc);
  void SetVx2Score(TCSlice& slc, VtxStore& vx2);
  void SetVx3Score(TCSlice& slc, Vtx3Store& vx3);
  void SetHighScoreBits(TCSlice& slc, Vtx3Store& vx3);
  bool MakeVertexObsolete(std::string fcnLabel, TCSlice& slc, VtxStore& vx2, bool forceKill);
  bool MakeVertexObsolete(TCSlice& slc, Vtx3Store& vx3);
  std::vector<int> GetVtxTjIDs(const TCSlice& slc, const VtxStore& vx2);
  std::vector<int> GetVtxTjIDs(const TCSlice& slc, const Vtx3Store& vx3, float& score);
  void PosInPlane(detinfo::DetectorPropertiesData const& detProp,
                  const Vtx3Store& vx3,
                  unsigned short plane,
                  Point2_t& pos);
  unsigned short IsCloseToVertex(const TCSlice& slc, const VtxStore& vx);
  unsigned short IsCloseToVertex(const TCSlice& slc, const Vtx3Store& vx3);
} // namespace

#endif // ifndef TRAJCLUSTERALGVERTEX_H
