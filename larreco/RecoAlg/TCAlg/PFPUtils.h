////////////////////////////////////////////////////////////////////////
//
//
// PFParticle  utilities
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGSPTUTILS_H
#define TRAJCLUSTERALGSPTUTILS_H

// C/C++ standard libraries
#include <string>

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}
namespace geo {
  struct TPCID;
}

namespace tca {

  void StitchPFPs();
  void FindPFParticles(detinfo::DetectorClocksData const& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       TCSlice& slc);
  void MakePFParticles(detinfo::DetectorClocksData const& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       TCSlice& slc,
                       std::vector<MatchStruct> matVec,
                       unsigned short matVec_Iter);
  bool ReconcileTPs(TCSlice& slc, PFPStruct& pfp, bool prt);
  void ReconcileTPs(TCSlice& slc);
  void MakePFPTjs(TCSlice& slc);
  void FillWireIntersections(TCSlice& slc);
  bool TCIntersectionPoint(unsigned int wir1,
                           unsigned int wir2,
                           unsigned int pln1,
                           unsigned int pln2,
                           float& y,
                           float& z);
  void Match3Planes(TCSlice& slc, std::vector<MatchStruct>& matVec);
  bool SptInTPC(const std::array<unsigned int, 3>& sptHits, unsigned int tpc);
  void Match2Planes(TCSlice& slc, std::vector<MatchStruct>& matVec);
  bool Update(detinfo::DetectorPropertiesData const& detProp,
              const TCSlice& slc,
              PFPStruct& pfp,
              bool prt);
  bool ReSection(detinfo::DetectorPropertiesData const& detProp,
                 const TCSlice& slc,
                 PFPStruct& pfp,
                 bool prt);
  void CountBadPoints(const TCSlice& slc,
                      const PFPStruct& pfp,
                      unsigned short fromPt,
                      unsigned short toPt,
                      unsigned short& nBadPts,
                      unsigned short& firstBadPt);
  void KillBadPoints(detinfo::DetectorPropertiesData const& detProp,
                     TCSlice& slc,
                     PFPStruct& pfp,
                     float pullCut,
                     bool prt);
  bool CanSection(const TCSlice& slc, const PFPStruct& pfp);
  unsigned short Find3DRecoRange(const TCSlice& slc,
                                 const PFPStruct& pfp,
                                 unsigned short fromPt,
                                 unsigned short min2DPts,
                                 short dir);
  void GetRange(const PFPStruct& pfp,
                unsigned short sfIndex,
                unsigned short& fromPt,
                unsigned short& npts);
  bool FitSection(detinfo::DetectorPropertiesData const& detProp,
                  const TCSlice& slc,
                  PFPStruct& pfp,
                  unsigned short sfIndex);
  SectionFit FitTP3Ds(detinfo::DetectorPropertiesData const& detProp,
                      const TCSlice& slc,
                      const std::vector<TP3D>& tp3ds,
                      unsigned short fromPt,
                      short fitDir,
                      unsigned short nPtsFit);
  bool FitTP3Ds(detinfo::DetectorPropertiesData const& detProp,
                const TCSlice& slc,
                PFPStruct& pfp,
                unsigned short fromPt,
                unsigned short npts,
                unsigned short sfIndex,
                float& chiDOF);
  void SplitAtKinks(detinfo::DetectorPropertiesData const& detProp,
                    const TCSlice& slc,
                    std::vector<PFPStruct>& pfpVec,
                    bool prt);
  void KinkFit(detinfo::DetectorPropertiesData const& detProp,
               const TCSlice& slc,
               const PFPStruct& pfp,
               unsigned short atPt,
               double fitLen,
               double& dang,
               double& dangErr);
  bool Split(detinfo::DetectorPropertiesData const& detProp,
             const TCSlice& slc,
             PFPStruct& p1,
             unsigned short atPt,
             PFPStruct& p2,
             bool prt);
  void ReconcileVertices(TCSlice& slc, PFPStruct& pfp, bool prt);
  void TrimEndPts(detinfo::DetectorPropertiesData const& detProp,
                  TCSlice& slc,
                  PFPStruct& pfp,
                  bool prt);
  void FillGaps3D(detinfo::DetectorClocksData const& clockData,
                  detinfo::DetectorPropertiesData const& detProp,
                  TCSlice& slc,
                  PFPStruct& pfp,
                  bool prt);
  bool ValidTwoPlaneMatch(detinfo::DetectorPropertiesData const& detProp,
                          const TCSlice& slc,
                          const PFPStruct& pfp);
  void AddPointsInRange(detinfo::DetectorClocksData const& clockData,
                        detinfo::DetectorPropertiesData const& detProp,
                        TCSlice& slc,
                        PFPStruct& pfp,
                        unsigned short fromPt,
                        unsigned short toPt,
                        CTP_t inCTP,
                        float maxPull,
                        unsigned short& nWires,
                        unsigned short& nAdd,
                        bool prt);
  unsigned short InsertTP3D(PFPStruct& pfp, TP3D& tp3d);
  bool SortSection(PFPStruct& pfp, unsigned short sectionFitIndex);
  bool MakeTP3Ds(detinfo::DetectorPropertiesData const& detProp, TCSlice& slc, PFPStruct& pfp);
  void Reverse(TCSlice& slc, PFPStruct& pfp);
  void FillmAllTraj(detinfo::DetectorPropertiesData const& detProp, TCSlice& slc);
  bool MakeTp3(TCSlice& slc,
               const TrajPoint& itp,
               const TrajPoint& jtp,
               TrajPoint3& tp3,
               bool findDirection);
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2);
  inline double
  DotProd(const Vector3_t& v1, const Vector3_t& v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2);
  double PosSep(const Point3_t& pos1, const Point3_t& pos2);
  double PosSep2(const Point3_t& pos1, const Point3_t& pos2);
  bool SetMag(Vector3_t& v1, double mag);
  void FilldEdx(detinfo::DetectorClocksData const& clockData,
                detinfo::DetectorPropertiesData const& detProp,
                const TCSlice& slc,
                PFPStruct& pfp);
  float dEdx(detinfo::DetectorClocksData const& clockData,
             detinfo::DetectorPropertiesData const& detProp,
             const TCSlice& slc,
             TP3D& tp3d);
  void Average_dEdX(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    const TCSlice& slc,
                    PFPStruct& pfp,
                    float& dEdXAve,
                    float& dEdXRms);
  TP3D CreateTP3D(detinfo::DetectorPropertiesData const& detProp,
                  const TCSlice& slc,
                  int tjID,
                  unsigned short tjPt);
  bool SetSection(detinfo::DetectorPropertiesData const& detProp,
                  const TCSlice& slc,
                  PFPStruct& pfp,
                  TP3D& tp3d);
  float PointPull(const PFPStruct& pfp, const TP3D& tp3d);
  PFPStruct CreatePFP(const TCSlice& slc);
  void PFPVertexCheck(TCSlice& tcs);
  void DefinePFPParents(TCSlice& slc, bool prt);
  bool StorePFP(TCSlice& slc, PFPStruct& pfp);
  bool InsideFV(const TCSlice& slc, const PFPStruct& pfp, unsigned short end);
  bool InsideTPC(const Point3_t& pos, geo::TPCID& inTPCID);
  void FindAlongTrans(Point3_t pos1, Vector3_t dir1, Point3_t pos2, Point2_t& alongTrans);
  bool PointDirIntersect(Point3_t p1,
                         Vector3_t p1Dir,
                         Point3_t p2,
                         Vector3_t p2Dir,
                         Point3_t& intersect,
                         float& doca);
  bool LineLineIntersect(Point3_t p1,
                         Point3_t p2,
                         Point3_t p3,
                         Point3_t p4,
                         Point3_t& intersect,
                         float& doca);
  float ChgFracBetween(detinfo::DetectorPropertiesData const& detProp,
                       const TCSlice& slc,
                       Point3_t pos1,
                       Point3_t pos2);
  float ChgFracNearEnd(detinfo::DetectorPropertiesData const& detProp,
                       const TCSlice& slc,
                       const PFPStruct& pfp,
                       unsigned short end);
  Point3_t PosAtEnd(const PFPStruct& pfp, unsigned short end);
  Vector3_t DirAtEnd(const PFPStruct& pfp, unsigned short end);
  float Length(const PFPStruct& pfp);
  bool SectionStartEnd(const PFPStruct& pfp,
                       unsigned short sfIndex,
                       unsigned short& startPt,
                       unsigned short& endPt);
  unsigned short FarEnd(const TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos);
  int PDGCodeVote(detinfo::DetectorClocksData const& clockData,
                  detinfo::DetectorPropertiesData const& detProp,
                  const TCSlice& slc,
                  PFPStruct& pfp);
  void PrintTP3Ds(detinfo::DetectorClocksData const& clockData,
                  detinfo::DetectorPropertiesData const& detProp,
                  std::string someText,
                  const TCSlice& slc,
                  const PFPStruct& pfp,
                  short printPts);
} // namespace tca

#endif // ifndef TRAJCLUSTERALGSPTUTILS_H
