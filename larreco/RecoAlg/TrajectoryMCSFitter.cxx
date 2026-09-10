#include "TrajectoryMCSFitter.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

using namespace std;
using namespace trkf;
using namespace recob::tracking;

namespace {
  inline Vector_t crossProduct(const Vector_t& a, const Vector_t& b)
  {
    return Vector_t(a.Y() * b.Z() - a.Z() * b.Y(),
                    a.Z() * b.X() - a.X() * b.Z(),
                    a.X() * b.Y() - a.Y() * b.X());
  }

  inline double dotProduct(const Vector_t& a, const Vector_t& b)
  {
    return a.X() * b.X() + a.Y() * b.Y() + a.Z() * b.Z();
  }

  inline double magnitude(const Vector_t& v) { return std::sqrt(dotProduct(v, v)); }
}

recob::MCSFitResult TrajectoryMCSFitter::fitMcs(const recob::TrackTrajectory& traj, int pid) const
{
  // Step 1: split the reconstructed trajectory into fixed-length segments.
  // Each segment contributes a material thickness in radiation lengths.
  vector<size_t> breakpoints;
  vector<float> segradlengths;
  vector<float> cumseglens;
  breakTrajInSegments(traj, breakpoints, segradlengths, cumseglens);

  // Step 2: fit one straight direction per segment and measure the signed
  // projected angle between each neighboring pair of segments.
  if (segradlengths.size() < 2) return recob::MCSFitResult();

  vector<float> dtheta;
  vector<float> dthetaX;
  vector<float> dthetaY;
  vector<float> angleDriftFrac;

  Vector_t pcdir0;
  Vector_t pcdir1;
  const Vector_t driftAxis(1.0, 0.0, 0.0);

  for (unsigned int p = 0; p < segradlengths.size(); ++p) {
    linearRegression(traj, breakpoints[p], breakpoints[p + 1], pcdir1);

    if (p > 0) {
      if (segradlengths[p] < -100. || segradlengths[p - 1] < -100.) {
        dtheta.push_back(-999.);
        dthetaX.push_back(-999.);
        dthetaY.push_back(-999.);
        angleDriftFrac.push_back(-999.);
      }
      else {
        const double cosval = dotProduct(pcdir0, pcdir1);
        const double dt = 1000. * std::acos(std::clamp(cosval, -1.0, 1.0));
        dtheta.push_back(dt);

        // Local frame for this segment pair:
        //   zhat = previous segment direction (pcdir0)
        //   yhat = zhat x detector-drift axis
        //   xhat = zhat x yhat
        // The next segment direction is projected onto xhat/yhat to form
        // theta_x'z' and theta_y'z'.  angleDriftFrac stores |v_x|, the absolute
        // drift-axis component of the average pair direction, for optional
        // direction-dependent y'-z' calibration.
        Vector_t yhat = crossProduct(pcdir0, driftAxis);
        const double ymag = magnitude(yhat);

        if (ymag < 1e-6) {
          // Singular when previous segment is parallel to drift direction.
          dthetaX.push_back(-999.);
          dthetaY.push_back(-999.);
          angleDriftFrac.push_back(-999.);
        }
        else {
          yhat = Vector_t(yhat.X() / ymag, yhat.Y() / ymag, yhat.Z() / ymag);

          Vector_t xhat = crossProduct(pcdir0, yhat);
          const double xmag = magnitude(xhat);

          if (xmag < 1e-6) {
            dthetaX.push_back(-999.);
            dthetaY.push_back(-999.);
            angleDriftFrac.push_back(-999.);
          }
          else {
            xhat = Vector_t(xhat.X() / xmag, xhat.Y() / xmag, xhat.Z() / xmag);

            const double vx = dotProduct(pcdir1, xhat);
            const double vy = dotProduct(pcdir1, yhat);
            const double vz = dotProduct(pcdir1, pcdir0);

            dthetaX.push_back(1000. * std::atan2(vx, vz)); // mrad, drift-sensitive
            dthetaY.push_back(1000. * std::atan2(vy, vz)); // mrad, drift-insensitive

            const Vector_t avgdir(pcdir0.X() + pcdir1.X(),
                                  pcdir0.Y() + pcdir1.Y(),
                                  pcdir0.Z() + pcdir1.Z());
            const double avgmag = magnitude(avgdir);
            const double dirVx =
              (avgmag > 1e-6 ? dotProduct(avgdir, driftAxis) / avgmag :
                                dotProduct(pcdir0, driftAxis));
            angleDriftFrac.push_back(static_cast<float>(std::abs(dirVx)));
          }
        }
      }
    }

    pcdir0 = pcdir1;
  }

  // Step 3: run the likelihood scan twice.  The forward scan assumes the track
  // starts at the trajectory start; the backward scan assumes it starts at the
  // trajectory end.  Energy loss is propagated from the assumed start point.
  vector<float> cumLenFwd;
  vector<float> cumLenBwd;
  for (unsigned int i = 0; i < cumseglens.size() - 2; ++i) {
    cumLenFwd.push_back(cumseglens[i]);
    cumLenBwd.push_back(cumseglens.back() - cumseglens[i + 2]);
  }

  const ScanResult fwdResult =
    doLikelihoodScan(dthetaX, dthetaY, angleDriftFrac, segradlengths, cumLenFwd, true, pid);
  const ScanResult bwdResult =
    doLikelihoodScan(dthetaX, dthetaY, angleDriftFrac, segradlengths, cumLenBwd, false, pid);

  return recob::MCSFitResult(pid,
                             fwdResult.p,
                             fwdResult.pUnc,
                             fwdResult.logL,
                             bwdResult.p,
                             bwdResult.pUnc,
                             bwdResult.logL,
                             segradlengths,
                             dtheta);
}

void TrajectoryMCSFitter::breakTrajInSegments(const recob::TrackTrajectory& traj,
                                              vector<size_t>& breakpoints,
                                              vector<float>& segradlengths,
                                              vector<float>& cumseglens) const
{
  // Split the input trajectory into approximately equal path-length segments.
  // breakpoints holds the trajectory-point index at each segment boundary.
  // segradlengths holds each segment length divided by the LAr radiation length.
  // cumseglens holds the cumulative path length from the start of the track.
  art::ServiceHandle<geo::Geometry const> geom;
  auto const* _SCE = (applySCEcorr_ ? lar::providerFrom<spacecharge::SpaceChargeService>() : NULL);

  const double trajlen = traj.Length();
  const double thisSegLen =
    (trajlen > (segLen_ * minNSegs_) ? segLen_ : trajlen / double(minNSegs_));

  // Liquid argon radiation length is approximated as 14 cm here, matching the
  // historical fitter convention used by this algorithm.
  constexpr double lar_radl_inv = 1. / 14.0;
  cumseglens.push_back(0.);
  double thislen = 0.;
  double totlen = 0.;
  auto nextValid = traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  auto pos0 = traj.LocationAtPoint(nextValid);
  if (applySCEcorr_) {
    geo::TPCID tpcid = geom->FindTPCAtPosition(pos0);
    geo::Vector_t pos0_offset(0., 0., 0.);
    if (tpcid.isValid) { pos0_offset = _SCE->GetCalPosOffsets(pos0, tpcid.TPC); }
    pos0.SetX(pos0.X() - pos0_offset.X());
    pos0.SetY(pos0.Y() + pos0_offset.Y());
    pos0.SetZ(pos0.Z() + pos0_offset.Z());
  }
  auto dir0 = traj.DirectionAtPoint(nextValid);
  nextValid = traj.NextValidPoint(nextValid + 1);
  int npoints = 0;
  while (nextValid != recob::TrackTrajectory::InvalidIndex) {
    if (npoints == 0) dir0 = traj.DirectionAtPoint(nextValid);
    auto pos1 = traj.LocationAtPoint(nextValid);
    if (applySCEcorr_) {
      geo::TPCID tpcid = geom->FindTPCAtPosition(pos1);
      geo::Vector_t pos1_offset(0., 0., 0.);
      if (tpcid.isValid) { pos1_offset = _SCE->GetCalPosOffsets(pos1, tpcid.TPC); }
      pos1.SetX(pos1.X() - pos1_offset.X());
      pos1.SetY(pos1.Y() + pos1_offset.Y());
      pos1.SetZ(pos1.Z() + pos1_offset.Z());
    }
    auto step = (pos1 - pos0).R();
    thislen += dir0.Dot(pos1 - pos0);
    totlen += step;
    pos0 = pos1;
    npoints++;
    if (thislen >= (thisSegLen - segLenTolerance_)) {
      breakpoints.push_back(nextValid);
      if (npoints >= minHitsPerSegment_)
        segradlengths.push_back(thislen * lar_radl_inv);
      else
        segradlengths.push_back(-999.);
      cumseglens.push_back(totlen);
      thislen = 0.;
      npoints = 0;
    }
    nextValid = traj.NextValidPoint(nextValid + 1);
  }
  // Add the final partial segment if any path length remains.
  if (thislen > 0.) {
    breakpoints.push_back(traj.LastValidPoint() + 1);
    segradlengths.push_back(thislen * lar_radl_inv);
    cumseglens.push_back(cumseglens.back() + thislen);
  }
  return;
}

const TrajectoryMCSFitter::ScanResult TrajectoryMCSFitter::doLikelihoodScan(
  std::vector<float>& dthetaX,
  std::vector<float>& dthetaY,
  std::vector<float>& angleDriftFrac,
  std::vector<float>& seg_nradlengths,
  std::vector<float>& cumLen,
  bool fwdFit,
  int pid,
  float pmin,
  float pmax,
  float pstep) const
{
  int best_idx = -1;
  float best_logL = std::numeric_limits<float>::max();
  float best_p = -1.0;
  std::vector<float> vlogL;
  for (float p_test = pmin; p_test <= pmax; p_test += pstep) {
    const float logL =
      mcsLikelihood(p_test, dthetaX, dthetaY, angleDriftFrac, seg_nradlengths, cumLen, fwdFit, pid);
    if (logL < best_logL) {
      best_p = p_test;
      best_logL = logL;
      best_idx = vlogL.size();
    }
    vlogL.push_back(logL);
  }
  //
  // uncertainty from left side scan
  float lunc = -1.;
  if (best_idx > 0) {
    for (int j = best_idx - 1; j >= 0; --j) {
      const float dLL = vlogL[j] - vlogL[best_idx];
      if (dLL >= 0.5) {
        lunc = (best_idx - j) * pstep;
        break;
      }
    }
  }
  // uncertainty from right side scan
  float runc = -1.;
  if (best_idx < int(vlogL.size() - 1)) {
    for (unsigned int j = best_idx + 1; j < vlogL.size(); ++j) {
      const float dLL = vlogL[j] - vlogL[best_idx];
      if (dLL >= 0.5) {
        runc = (j - best_idx) * pstep;
        break;
      }
    }
  }
  return ScanResult(best_p, std::max(lunc, runc), best_logL);
}

const TrajectoryMCSFitter::ScanResult TrajectoryMCSFitter::doLikelihoodScan(
  std::vector<float>& dthetaX,
  std::vector<float>& dthetaY,
  std::vector<float>& angleDriftFrac,
  std::vector<float>& seg_nradlengths,
  std::vector<float>& cumLen,
  bool fwdFit,
  int pid) const
{
  // First pass: coarse momentum scan over the configured full range.
  const ScanResult coarseRes = doLikelihoodScan(
    dthetaX, dthetaY, angleDriftFrac, seg_nradlengths, cumLen, fwdFit, pid,
    pMin_, pMax_, pStepCoarse_);

  float pmax = std::min(coarseRes.p + fineScanWindow_, pMax_);
  float pmin = std::max(coarseRes.p - fineScanWindow_, pMin_);
  if (coarseRes.pUnc < (std::numeric_limits<float>::max() - 1.)) {
    pmax = std::min(coarseRes.p + 2 * coarseRes.pUnc, pMax_);
    pmin = std::max(coarseRes.p - 2 * coarseRes.pUnc, pMin_);
  }

  // Second pass: fine scan around the coarse minimum.
  const ScanResult refineRes =
    doLikelihoodScan(dthetaX, dthetaY, angleDriftFrac, seg_nradlengths, cumLen, fwdFit, pid,
                     pmin, pmax, pStep_);

  return refineRes;
}

void TrajectoryMCSFitter::linearRegression(const recob::TrackTrajectory& traj,
                                           const size_t firstPoint,
                                           const size_t lastPoint,
                                           Vector_t& pcdir) const
{
  // Fit a straight direction to one segment by principal-component analysis.
  // The largest eigenvector of the point covariance matrix is taken as the
  // segment direction, then flipped to follow the reconstructed track order.
  art::ServiceHandle<geo::Geometry const> geom;
  auto const* _SCE = (applySCEcorr_ ? lar::providerFrom<spacecharge::SpaceChargeService>() : NULL);

  int npoints = 0;
  geo::vect::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  while (nextValid < lastPoint) {
    auto tempP = traj.LocationAtPoint(nextValid);
    if (applySCEcorr_) {
      geo::TPCID tpcid = geom->FindTPCAtPosition(tempP);
      geo::Vector_t tempP_offset(0., 0., 0.);
      if (tpcid.isValid) { tempP_offset = _SCE->GetCalPosOffsets(tempP, tpcid.TPC); }
      tempP.SetX(tempP.X() - tempP_offset.X());
      tempP.SetY(tempP.Y() + tempP_offset.Y());
      tempP.SetZ(tempP.Z() + tempP_offset.Z());
    }
    middlePointCalc.add(tempP);
    nextValid = traj.NextValidPoint(nextValid + 1);
    npoints++;
  }
  const auto avgpos = middlePointCalc.middlePoint();
  const double norm = 1. / double(npoints);
  //
  TMatrixDSym m(3);
  nextValid = firstPoint;
  while (nextValid < lastPoint) {
    auto p = traj.LocationAtPoint(nextValid);
    if (applySCEcorr_) {
      geo::TPCID tpcid = geom->FindTPCAtPosition(p);
      geo::Vector_t p_offset(0., 0., 0.);
      if (tpcid.isValid) { p_offset = _SCE->GetCalPosOffsets(p, tpcid.TPC); }
      p.SetX(p.X() - p_offset.X());
      p.SetY(p.Y() + p_offset.Y());
      p.SetZ(p.Z() + p_offset.Z());
    }
    const double xxw0 = p.X() - avgpos.X();
    const double yyw0 = p.Y() - avgpos.Y();
    const double zzw0 = p.Z() - avgpos.Z();
    m(0, 0) += xxw0 * xxw0 * norm;
    m(0, 1) += xxw0 * yyw0 * norm;
    m(0, 2) += xxw0 * zzw0 * norm;
    m(1, 0) += yyw0 * xxw0 * norm;
    m(1, 1) += yyw0 * yyw0 * norm;
    m(1, 2) += yyw0 * zzw0 * norm;
    m(2, 0) += zzw0 * xxw0 * norm;
    m(2, 1) += zzw0 * yyw0 * norm;
    m(2, 2) += zzw0 * zzw0 * norm;
    nextValid = traj.NextValidPoint(nextValid + 1);
  }
  //
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i = 1; i < 3; ++i) {
    if (eigenval(i) > maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  }
  //
  pcdir = Vector_t(eigenvec(0, maxevalidx), eigenvec(1, maxevalidx), eigenvec(2, maxevalidx));
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir) < 0.) pcdir *= -1.;
  //
}

double TrajectoryMCSFitter::mcsLikelihood(double p,
                                          std::vector<float>& dthetaX,
                                          std::vector<float>& dthetaY,
                                          std::vector<float>& angleDriftFrac,
                                          std::vector<float>& seg_nradl,
                                          std::vector<float>& cumLen,
                                          bool fwd,
                                          int pid) const
{
  //
  const int beg = (fwd ? 0 : (dthetaX.size() - 1));
  const int end = (fwd ? dthetaX.size() : -1);
  const int incr = (fwd ? +1 : -1);
  //
  const double m = mass(pid);
  const double m2 = m * m;
  const double Etot = std::sqrt(p * p + m2);
  double Eij2 = 0.;
  //
  double result = 0.;
  for (int i = beg; i != end; i += incr) {
    if (dthetaX[i] < -900. || dthetaY[i] < -900.) {
      continue;
    }
    //
    const double Eij = GetE(Etot, cumLen[i], m);
    Eij2 = Eij * Eij;
    //
    if (Eij2 <= m2) {
      result = std::numeric_limits<double>::max();
      break;
    }
    const double pij = std::sqrt(Eij2 - m2);
    const double kineticEnergyMeV = 1000.0 * std::max(Eij - m, 0.0);
    const double beta = std::sqrt(1. - (m2 / (pij * pij + m2)));
    constexpr double HighlandSecondTerm = 0.038;
    const double tH0 = (HighlandFirstTerm(pij) / (pij * beta)) *
                       (1.0 + HighlandSecondTerm * std::log(seg_nradl[i])) *
                       std::sqrt(seg_nradl[i]);

    // Build the calibrated two-Gaussian widths for this candidate momentum.
    // The Highland term supplies the expected MCS scattering scale, while the
    // smooth calibration functions adjust that scale and add detector/reco
    // resolution floors.
    double sigma1X = 0.0;
    double sigma2X = 0.0;
    double sigma1Y = 0.0;
    double sigma2Y = 0.0;
    double fracX = 0.0;
    double fracY = 0.0;

    // Calibrated double-Gaussian PDF:
    //   PDF(theta) = r * G(theta; sigma1) + (1-r) * G(theta; sigma2)
    // sigma1 is the primary/core width, sigma2 is the broad tail width, and
    // r is the primary Gaussian area fraction.
    const double norm1D = 1.0 / std::sqrt(2.0 * M_PI);
    const bool useDirY =
      useDirectionDependentYZ_ &&
      static_cast<size_t>(i) < angleDriftFrac.size() &&
      angleDriftFrac[i] >= 0.0;
    const size_t yDirBin = (useDirY ? SmoothYDirectionBin(angleDriftFrac[i]) : 0);

    const auto smoothYPdfForDirBin = [&](const size_t dirBin) -> double {
      const double pivotY = DirectionValue(smoothPivotYByDir_, dirBin, smoothPivotY_);
      const double scale1Y =
        DirectionSmoothScale(pij, smoothScale1YByDir_, smoothScale1Y_, pivotY, dirBin);
      const double scale2Y =
        DirectionSmoothScale(pij, smoothScale2YByDir_, smoothScale2Y_, pivotY, dirBin);
      const double res1Y = DirectionValue(smoothRes1YByDir_, dirBin, smoothRes1Y_);
      const double res2Y = DirectionValue(smoothRes2YByDir_, dirBin, smoothRes2Y_);
      const double sigma1 = std::sqrt(std::pow(scale1Y * tH0, 2) + std::pow(res1Y, 2));
      const double sigma2 = std::sqrt(std::pow(scale2Y * tH0, 2) + std::pow(res2Y, 2));

      if (sigma1 <= 0.0 || sigma2 <= 0.0) {
        return std::numeric_limits<double>::min();
      }

      const double fracLowY = DirectionValue(smoothFracLowYByDir_, dirBin, smoothFracLowY_);
      const double fracHighY = DirectionValue(smoothFracHighYByDir_, dirBin, smoothFracHighY_);
      const double fracSlopeY = DirectionValue(smoothFracSlopeYByDir_, dirBin, smoothFracSlopeY_);
      const double fracMidY = DirectionValue(smoothFracMidYByDir_, dirBin, smoothFracMidY_);
      const double frac = SmoothFrac(pij, fracLowY, fracHighY, fracSlopeY, fracMidY);

      const double g1 = (norm1D / sigma1) * std::exp(-0.5 * std::pow(dthetaY[i] / sigma1, 2));
      const double g2 = (norm1D / sigma2) * std::exp(-0.5 * std::pow(dthetaY[i] / sigma2, 2));
      return std::max(frac * g1 + (1.0 - frac) * g2, std::numeric_limits<double>::min());
    };

    const double scale1X = SmoothScale(pij, smoothScale1X_, smoothPivotX_);
    const double scale2X = SmoothScale(pij, smoothScale2X_, smoothPivotX_);
    const double pivotY =
      (useDirY ? DirectionValue(smoothPivotYByDir_, yDirBin, smoothPivotY_) : smoothPivotY_);
    const double scale1Y =
      (useDirY ? DirectionSmoothScale(pij, smoothScale1YByDir_, smoothScale1Y_, pivotY, yDirBin) :
                 SmoothScale(pij, smoothScale1Y_, smoothPivotY_));
    const double scale2Y =
      (useDirY ? DirectionSmoothScale(pij, smoothScale2YByDir_, smoothScale2Y_, pivotY, yDirBin) :
                 SmoothScale(pij, smoothScale2Y_, smoothPivotY_));
    const double res1Y =
      (useDirY ? DirectionValue(smoothRes1YByDir_, yDirBin, smoothRes1Y_) : smoothRes1Y_);
    const double res2Y =
      (useDirY ? DirectionValue(smoothRes2YByDir_, yDirBin, smoothRes2Y_) : smoothRes2Y_);

    sigma1X = std::sqrt(std::pow(scale1X * tH0, 2) + std::pow(smoothRes1X_, 2));
    sigma2X = std::sqrt(std::pow(scale2X * tH0, 2) + std::pow(smoothRes2X_, 2));
    sigma1Y = std::sqrt(std::pow(scale1Y * tH0, 2) + std::pow(res1Y, 2));
    sigma2Y = std::sqrt(std::pow(scale2Y * tH0, 2) + std::pow(res2Y, 2));

    fracX = SmoothFrac(pij,
                       smoothFracLowX_,
                       smoothFracHighX_,
                       smoothFracSlopeX_,
                       smoothFracMidX_);
    const double fracLowY =
      (useDirY ? DirectionValue(smoothFracLowYByDir_, yDirBin, smoothFracLowY_) : smoothFracLowY_);
    const double fracHighY =
      (useDirY ? DirectionValue(smoothFracHighYByDir_, yDirBin, smoothFracHighY_) : smoothFracHighY_);
    const double fracSlopeY =
      (useDirY ? DirectionValue(smoothFracSlopeYByDir_, yDirBin, smoothFracSlopeY_) : smoothFracSlopeY_);
    const double fracMidY =
      (useDirY ? DirectionValue(smoothFracMidYByDir_, yDirBin, smoothFracMidY_) : smoothFracMidY_);
    fracY = SmoothFrac(pij,
                       fracLowY,
                       fracHighY,
                       fracSlopeY,
                       fracMidY);

    if (sigma1X <= 0.0 || sigma2X <= 0.0 || sigma1Y <= 0.0 || sigma2Y <= 0.0) {
      std::cout << "Error: RMS cannot be zero!" << std::endl;
      return std::numeric_limits<double>::max();
    }

    const double argX = dthetaX[i] / sigma1X;
    const double argY = dthetaY[i] / sigma1Y;
    const double g1X = (norm1D / sigma1X) * std::exp(-0.5 * argX * argX);
    const double g2X = (norm1D / sigma2X) * std::exp(-0.5 * std::pow(dthetaX[i] / sigma2X, 2));
    const double g1Y = (norm1D / sigma1Y) * std::exp(-0.5 * argY * argY);
    const double g2Y = (norm1D / sigma2Y) * std::exp(-0.5 * std::pow(dthetaY[i] / sigma2Y, 2));

    const double pdfX =
      std::max(fracX * g1X + (1.0 - fracX) * g2X, std::numeric_limits<double>::min());
    double pdfY =
      std::max(fracY * g1Y + (1.0 - fracY) * g2Y, std::numeric_limits<double>::min());

    if (useDirectionDependentYZ_ && useYZHighVxFallback_ && useDirY) {
      if (yDirBin == 3) {
        const double wHigh = YZBlendWeight(kineticEnergyMeV, yzBlendKEEdgesMeV_[2]);
        pdfY = (1.0 - wHigh) * smoothYPdfForDirBin(3) + wHigh * smoothYPdfForDirBin(2);
      }
      else if (yDirBin == 4) {
        if (kineticEnergyMeV < yzBlendKEEdgesMeV_[1]) {
          const double wLow = YZBlendWeight(kineticEnergyMeV, yzBlendKEEdgesMeV_[0]);
          pdfY = (1.0 - wLow) * smoothYPdfForDirBin(4) + wLow * smoothYPdfForDirBin(3);
        }
        else {
          const double wHigh = YZBlendWeight(kineticEnergyMeV, yzBlendKEEdgesMeV_[2]);
          pdfY = (1.0 - wHigh) * smoothYPdfForDirBin(3) + wHigh * smoothYPdfForDirBin(2);
        }
      }
      else {
        pdfY = smoothYPdfForDirBin(yDirBin);
      }

      pdfY = std::max(pdfY, std::numeric_limits<double>::min());
    }

    result += -std::log(pdfX);
    result += -std::log(pdfY);
  }
  return result;
}

double TrajectoryMCSFitter::energyLossLandau(const double mass2,
                                             const double e2,
                                             const double x) const
{
  // Most-probable energy loss over a step x.  This is the default propagation
  // mode used by the historical MicroBooNE MCS fitter.
  if (x <= 0.) return 0.;
  constexpr double Iinv2 = 1. / (188.E-6 * 188.E-6);
  constexpr double matConst = 1.4 * 18. / 40.;
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  constexpr double j = 0.200;
  //
  const double beta2 = (e2 - mass2) / e2;
  const double gamma2 = 1. / (1.0 - beta2);
  const double epsilon = 0.5 * kappa * x * matConst / beta2;
  //
  return 0.001 * epsilon * (log(2. * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2);
}

double TrajectoryMCSFitter::energyLossBetheBloch(const double mass, const double e2) const
{
  // Mean Bethe-Bloch energy loss.  This remains available through eLossMode=2.
  constexpr double Iinv = 1. / 188.E-6;
  constexpr double matConst = 1.4 * 18. / 40.;
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  //
  const double beta2 = (e2 - mass * mass) / e2;
  const double gamma2 = 1. / (1.0 - beta2);
  const double massRatio = me / mass;
  const double argument = (2. * me * gamma2 * beta2 * Iinv) *
                          std::sqrt(1 + 2 * std::sqrt(gamma2) * massRatio + massRatio * massRatio);
  //
  double dedx = kappa * matConst / beta2;
  //
  if (mass == 0.0) return 0.0;
  if (argument <= exp(beta2)) {
    dedx = 0.;
  }
  else {
    dedx *= (log(argument) - beta2) * 1.E-3;
    if (dedx < 0.) dedx = 0.;
  }
  return dedx;
}

double TrajectoryMCSFitter::GetE(const double initial_E,
                                 const double length_travelled,
                                 const double m) const
{
  // Propagate the candidate total energy from the assumed track start to the
  // segment being evaluated.  The likelihood uses this segment-local momentum
  // in the Highland scattering term.
  if (eLossMode_ == 1) {
    constexpr double kcal = 0.002105;
    return (initial_E - kcal * length_travelled);
  }
  //
  const double step_size = length_travelled / nElossSteps_;
  //
  double current_E = initial_E;
  const double m2 = m * m;
  //
  for (auto i = 0; i < nElossSteps_; ++i) {
    if (eLossMode_ == 2) {
      double dedx = energyLossBetheBloch(m, current_E);
      current_E -= (dedx * step_size);
    }
    else {
      current_E -= energyLossLandau(m2, current_E * current_E, step_size);
    }
    if (current_E <= m) {
      return 0.;
    }
  }
  return current_E;
}


