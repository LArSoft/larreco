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

#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMatrixT.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorT.h"

#include <cmath>
#include <utility>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larreco/RecoAlg/PMAlg/PmaHit3D.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"

#include "range/v3/algorithm.hpp"
#include "range/v3/numeric.hpp"
#include "range/v3/view.hpp"

double
pma::Dist2(const TVector2& v1, const TVector2& v2)
{
  double const dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y();
  return cet::sum_of_squares(dx, dy);
}
double
pma::Dist2(const Vector2D& v1, const Vector2D& v2)
{
  double const dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y();
  return cet::sum_of_squares(dx, dy);
}

size_t
pma::GetHitsCount(const std::vector<pma::Hit3D*>& hits, unsigned int view)
{
  if (view == geo::kUnknown) { return hits.size(); }
  return ranges::count_if(hits, [view](auto hit) { return view == hit->View2D(); });
}

double
pma::GetSummedADC(const std::vector<pma::Hit3D*>& hits, unsigned int view)
{
  using namespace ranges;
  auto to_summed_adc = [](auto hit) { return hit->SummedADC(); };
  if (view == geo::kUnknown) { return accumulate(hits | views::transform(to_summed_adc), 0.); }
  return accumulate(hits | views::filter([view](auto hit) { return view == hit->View2D(); }) |
                      views::transform(to_summed_adc),
                    0.);
}

double
pma::GetSummedAmpl(const std::vector<pma::Hit3D*>& hits, unsigned int view)
{
  using namespace ranges;
  auto to_amplitude = [](auto hit) { return hit->GetAmplitude(); };
  if (view == geo::kUnknown) { return accumulate(hits | views::transform(to_amplitude), 0.); }
  return accumulate(hits | views::filter([view](auto hit) { return view == hit->View2D(); }) |
                      views::transform(to_amplitude),
                    0.);
}

double
pma::GetHitsRadius3D(const std::vector<pma::Hit3D*>& hits, bool exact)
{
  if (hits.empty()) return 0.0;

  if (!exact && (hits.size() < 5)) return 0.0;

  using namespace ranges;
  auto to_3d_point = [](auto hit) -> decltype(auto) { return hit->Point3D(); };
  auto const mean_point =
    accumulate(hits | views::transform(to_3d_point), TVector3{}) * (1. / hits.size());

  auto to_dist2_from_mean = [&mean_point](auto hit) {
    return pma::Dist2(hit->Point3D(), mean_point);
  };
  auto const max_r2 = max(hits | views::transform(to_dist2_from_mean));
  return sqrt(max_r2);
}

double
pma::GetHitsRadius2D(const std::vector<pma::Hit3D*>& hits, bool exact)
{
  if (hits.empty()) return 0.0;

  if (!exact && (hits.size() < 5)) return 0.0;

  using namespace ranges;
  auto to_2d_point = [](auto hit) -> decltype(auto) { return hit->Point2D(); };
  auto const mean_point =
    accumulate(hits | views::transform(to_2d_point), TVector2{}) * (1. / hits.size());

  auto to_dist2_from_mean = [&mean_point](auto hit) {
    return pma::Dist2(hit->Point2D(), mean_point);
  };
  auto const max_r2 = max(hits | views::transform(to_dist2_from_mean));
  return sqrt(max_r2);
}

double
pma::GetSegmentProjVector(const TVector2& p, const TVector2& p0, const TVector2& p1)
{
  TVector2 const v0(p - p0);
  TVector2 const v1(p1 - p0);
  return v0 * v1 / v1.Mod2();
}

double
pma::GetSegmentProjVector(const pma::Vector2D& p, const pma::Vector2D& p0, const pma::Vector2D& p1)
{
  pma::Vector2D const v0(p - p0);
  pma::Vector2D const v1(p1 - p0);
  return v0.Dot(v1) / v1.Mag2();
}

double
pma::GetSegmentProjVector(const TVector3& p, const TVector3& p0, const TVector3& p1)
{
  TVector3 const v0(p - p0);
  TVector3 const v1(p1 - p0);
  return v0.Dot(v1) / v1.Mag2();
}

double
pma::GetSegmentProjVector(const pma::Vector3D& p, const pma::Vector3D& p0, const pma::Vector3D& p1)
{
  pma::Vector3D const v0(p - p0);
  pma::Vector3D const v1(p1 - p0);
  return v0.Dot(v1) / v1.Mag2();
}

TVector2
pma::GetProjectionToSegment(const TVector2& p, const TVector2& p0, const TVector2& p1)
{
  TVector2 const v1(p1 - p0);
  double const b = GetSegmentProjVector(p, p0, p1);
  return p0 + v1 * b;
}

TVector3
pma::GetProjectionToSegment(const TVector3& p, const TVector3& p0, const TVector3& p1)
{
  TVector3 const v1(p1 - p0);
  double const b = GetSegmentProjVector(p, p0, p1);
  return p0 + v1 * b;
}

double
pma::SolveLeastSquares3D(const std::vector<std::pair<TVector3, TVector3>>& lines, TVector3& result)
{
  // RS: please, ask me if you need examples/explanation of formulas as they
  // are not easy to derive from the code solely; I have Mathcad sources that
  // were used to test the solving method, weighting, etc.

  result.SetXYZ(0., 0., 0.);
  if (lines.size() < 2) {
    mf::LogError("pma::SolveLeastSquares3D") << "Need min. two lines.";
    return -1.0;
  }

  double m;
  std::vector<TVectorT<double>> U, P;
  for (size_t v = 0; v < lines.size(); v++) {
    TVector3 point = lines[v].first;
    TVector3 dir = lines[v].second;
    dir -= point;
    m = dir.Mag();
    if (m > 0.0) {
      dir *= 1.0 / m;

      P.push_back(TVectorT<double>(3));
      P.back()[0] = point.X();
      P.back()[1] = point.Y();
      P.back()[2] = point.Z();

      U.push_back(TVectorT<double>(3));
      U.back()[0] = dir.X();
      U.back()[1] = dir.Y();
      U.back()[2] = dir.Z();
    }
    else
      mf::LogWarning("pma::SolveLeastSquares3D") << "Line undefined.";
  }
  if (P.size() < 2) {
    mf::LogError("pma::SolveLeastSquares3D") << "Need min. two lines.";
    return -1.0;
  }

  TVectorT<double> x(3), y(3), w(3);
  TMatrixT<double> A(3, 3);
  double ur, uc, pc;
  double s_uc2[3], s_ur_uc[3];
  double s_p_uc2[3], s_p_ur_uc[3];

  w[0] = 1.0;
  w[1] = 1.0;
  w[2] = 1.0;
  for (size_t r = 0; r < 3; r++) {
    y[r] = 0.0;
    for (size_t c = 0; c < 3; c++) {
      s_uc2[c] = 0.0;
      s_ur_uc[c] = 0.0;
      s_p_uc2[c] = 0.0;
      s_p_ur_uc[c] = 0.0;

      for (size_t v = 0; v < P.size(); v++) {
        //w[1] = fWeights[v]; // to remember that individual coordinates can be supressed...
        //w[2] = fWeights[v];

        ur = U[v][r];
        uc = U[v][c];
        pc = P[v][c];

        s_uc2[c] += w[r] * w[c] * (1 - uc * uc);
        s_p_uc2[c] += w[r] * w[r] * pc * (1 - uc * uc);

        s_ur_uc[c] += w[r] * w[c] * ur * uc;
        s_p_ur_uc[c] += w[r] * w[r] * pc * ur * uc;
      }

      if (r == c) {
        y[r] += s_p_uc2[c];
        A(r, c) = s_uc2[c];
      }
      else {
        y[r] -= s_p_ur_uc[c];
        A(r, c) = -s_ur_uc[c];
      }
    }
  }
  try {
    x = A.InvertFast() * y;
  }
  catch (...) {
    result.SetXYZ(0., 0., 0.);
    return 1.0e12;
  }

  result.SetXYZ(x[0], x[1], x[2]);

  double mse = 0.0;
  for (size_t v = 0; v < lines.size(); v++) {
    TVector3 const pproj = pma::GetProjectionToSegment(result, lines[v].first, lines[v].second);

    double const dx = result.X() - pproj.X(); // dx, dy, dz and the result point can be weighted
    double const dy = result.Y() - pproj.Y(); // here (linearly) by each line uncertainty
    double const dz = result.Z() - pproj.Z();
    mse += cet::sum_of_squares(dx, dy, dz);
  }
  return mse / lines.size();
}

TVector2
pma::GetProjectionToPlane(const TVector3& p,
                          unsigned int plane,
                          unsigned int tpc,
                          unsigned int cryo)
{
  art::ServiceHandle<geo::Geometry const> geom;

  return TVector2(geom->TPC(tpc, cryo).Plane(plane).PlaneCoordinate(p), p.X());
}

TVector2
pma::GetVectorProjectionToPlane(const TVector3& v,
                                unsigned int plane,
                                unsigned int tpc,
                                unsigned int cryo)
{
  TVector3 v0_3d(0., 0., 0.);
  TVector2 v0_2d = GetProjectionToPlane(v0_3d, plane, tpc, cryo);
  TVector2 v1_2d = GetProjectionToPlane(v, plane, tpc, cryo);

  return v1_2d - v0_2d;
}

TVector2
pma::WireDriftToCm(detinfo::DetectorPropertiesData const& detProp,
                   unsigned int wire,
                   float drift,
                   unsigned int plane,
                   unsigned int tpc,
                   unsigned int cryo)
{
  art::ServiceHandle<geo::Geometry const> geom;
  return TVector2(geom->TPC(tpc, cryo).Plane(plane).WirePitch() * wire,
                  detProp.ConvertTicksToX(drift, plane, tpc, cryo));
}

TVector2
pma::CmToWireDrift(detinfo::DetectorPropertiesData const& detProp,
                   float xw,
                   float yd,
                   unsigned int plane,
                   unsigned int tpc,
                   unsigned int cryo)
{
  art::ServiceHandle<geo::Geometry const> geom;
  return TVector2(xw / geom->TPC(tpc, cryo).Plane(plane).WirePitch(),
                  detProp.ConvertXToTicks(yd, plane, tpc, cryo));
}

bool
pma::bTrajectory3DOrderLess::operator()(pma::Hit3D* h1, pma::Hit3D* h2)
{
  if (h1 && h2)
    return h1->fSegFraction < h2->fSegFraction;
  else
    return false;
}

bool
pma::bTrajectory3DDistLess::operator()(pma::Hit3D* h1, pma::Hit3D* h2)
{
  if (h1 && h2)
    return h1->GetDist2ToProj() < h2->GetDist2ToProj();
  else
    return false;
}

bool
pma::bTrack3DLonger::operator()(const pma::TrkCandidate& t1, const pma::TrkCandidate& t2)
{
  pma::Track3D* trk1 = t1.Track();
  pma::Track3D* trk2 = t2.Track();
  if (trk1 && trk2) {
    double l1 = pma::Dist2(trk1->front()->Point3D(), trk1->back()->Point3D());
    double l2 = pma::Dist2(trk2->front()->Point3D(), trk2->back()->Point3D());
    return l1 > l2;
  }
  else
    return false;
}

pma::bSegmentProjLess::bSegmentProjLess(const TVector3& s0, const TVector3& s1)
  : segStart(s0), segStop(s1)
{
  if (s0 == s1) mf::LogError("pma::bSegmentProjLess") << "Vectors equal!";
}
