/**
 *  @file   Utilities.h
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Some geometrical functions and sorting helpers.
 *          See PmaTrack3D.h file for details.
 */

#ifndef Utilities_h
#define Utilities_h

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
namespace detinfo {
  class DetectorPropertiesData;
}

#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/DisplacementVector2D.h"

#include "TVector2.h"
#include "TVector3.h"

#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace pma {
  typedef ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>> Vector2D;
  typedef recob::tracking::Vector_t Vector3D;

  typedef std::map<size_t, std::vector<double>> dedx_map;

  class Hit3D;
  class TrkCandidate;
  class bSegmentProjLess;
  class bDistCenterLess2D;
  class bDistCenterLess3D;
  struct bTrajectory3DOrderLess;
  struct bTrajectory3DDistLess;
  struct bTrack3DLonger;

  double Dist2(const TVector2& v1, const TVector2& v2);
  double Dist2(const Vector2D& v1, const Vector2D& v2);

  template <typename T, typename U>
  double
  Dist2(const T& v1, const U& v2)
  {
    double dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y(), dz = v1.Z() - v2.Z();
    return dx * dx + dy * dy + dz * dz;
  }

  size_t GetHitsCount(const std::vector<pma::Hit3D*>& hits, unsigned int view);
  double GetSummedADC(const std::vector<pma::Hit3D*>& hits, unsigned int view = geo::kUnknown);
  double GetSummedAmpl(const std::vector<pma::Hit3D*>& hits, unsigned int view = geo::kUnknown);

  double GetHitsRadius3D(const std::vector<pma::Hit3D*>& hits, bool exact = false);
  double GetHitsRadius2D(const std::vector<pma::Hit3D*>& hits, bool exact = false);

  double GetSegmentProjVector(const TVector2& p, const TVector2& p0, const TVector2& p1);
  double GetSegmentProjVector(const Vector2D& p, const Vector2D& p0, const Vector2D& p1);

  double GetSegmentProjVector(const TVector3& p, const TVector3& p0, const TVector3& p1);
  double GetSegmentProjVector(const Vector3D& p, const Vector3D& p0, const Vector3D& p1);

  TVector2 GetProjectionToSegment(const TVector2& p, const TVector2& p0, const TVector2& p1);
  TVector3 GetProjectionToSegment(const TVector3& p, const TVector3& p0, const TVector3& p1);

  double SolveLeastSquares3D(const std::vector<std::pair<TVector3, TVector3>>& lines,
                             TVector3& result);

  TVector2 GetProjectionToPlane(const TVector3& p,
                                unsigned int plane,
                                unsigned int tpc,
                                unsigned int cryo);
  TVector2 GetVectorProjectionToPlane(const TVector3& v,
                                      unsigned int plane,
                                      unsigned int tpc,
                                      unsigned int cryo);
  TVector2 WireDriftToCm(detinfo::DetectorPropertiesData const& detProp,
                         unsigned int wire,
                         float drift,
                         unsigned int plane,
                         unsigned int tpc,
                         unsigned int cryo);
  TVector2 CmToWireDrift(detinfo::DetectorPropertiesData const& detProp,
                         float xw,
                         float yd,
                         unsigned int plane,
                         unsigned int tpc,
                         unsigned int cryo);
}

struct pma::bTrajectory3DOrderLess : public std::binary_function<pma::Hit3D*, pma::Hit3D*, bool> {
  bool operator()(pma::Hit3D* h1, pma::Hit3D* h2);
};

struct pma::bTrajectory3DDistLess : public std::binary_function<pma::Hit3D*, pma::Hit3D*, bool> {
  bool operator()(pma::Hit3D* h1, pma::Hit3D* h2);
};

struct pma::bTrack3DLonger
  : public std::binary_function<const pma::TrkCandidate&, const pma::TrkCandidate&, bool> {
  bool operator()(const pma::TrkCandidate& t1, const pma::TrkCandidate& t2);
};

class pma::bSegmentProjLess : public std::binary_function<TVector3*, TVector3*, bool> {
public:
  bSegmentProjLess(const TVector3& s0, const TVector3& s1);

  bool
  operator()(TVector3* p1, TVector3* p2)
  {
    if (p1 && p2) {
      double b1 = pma::GetSegmentProjVector(*p1, segStart, segStop);
      double b2 = pma::GetSegmentProjVector(*p1, segStart, segStop);
      return b1 < b2;
    }
    else
      return false;
  }

private:
  TVector3 segStart, segStop;
};

class pma::bDistCenterLess2D : public std::binary_function<TVector2, TVector2, bool> {
public:
  bDistCenterLess2D(const TVector2& c) : center(c) {}

  bool
  operator()(TVector2 p1, TVector2 p2)
  {
    double b1 = pma::Dist2(p1, center);
    double b2 = pma::Dist2(p2, center);
    return b1 < b2;
  }

private:
  TVector2 center;
};

class pma::bDistCenterLess3D : public std::binary_function<TVector3, TVector3, bool> {
public:
  bDistCenterLess3D(const TVector3& c) : center(c) {}

  bool
  operator()(TVector3 p1, TVector3 p2)
  {
    double b1 = pma::Dist2(p1, center);
    double b2 = pma::Dist2(p2, center);
    return b1 < b2;
  }

private:
  TVector3 center;
};

#endif
