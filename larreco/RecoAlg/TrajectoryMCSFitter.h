#ifndef TRAJECTORYMCSFITTER_H
#define TRAJECTORYMCSFITTER_H

// Framework includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "lardata/RecoObjects/TrackState.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace trkf {
  /**
   * @file  larreco/RecoAlg/TrajectoryMCSFitter.h
   * @class trkf::TrajectoryMCSFitter
   *
   * @brief Class for Maximum Likelihood fit of Multiple Coulomb Scattering
   *        angles between segments within a Track or Trajectory
   *
   * Main inputs are a reconstructed Track/Trajectory and the fitter
   * configuration: particle hypothesis, segmenting choices, momentum scan
   * range, energy-loss model, Highland scattering term, and calibrated
   * double-Gaussian angle-PDF parameters.
   *
   * Outputs are: a recob::MCSFitResult, containing:
   *   resulting momentum, momentum uncertainty, and best likelihood value
   *   (both for fwd and bwd fit);
   *   vector of cumulative segment (radiation) lengths, vector of scattering
   *   angles, and PID hypothesis used in the fit.
   *
   * The likelihood is evaluated using projected angles in the local x'-z' and
   * y'-z' planes.  The local z' axis is the previous segment direction.  The
   * local y' axis is built from z' cross detector-drift axis, and x' completes
   * the right-handed local frame.
   *
   * Each projected angle is evaluated with a calibrated two-Gaussian PDF:
   *   PDF(theta) = frac(p) * G(theta; sigma1(p))
   *              + (1 - frac(p)) * G(theta; sigma2(p)).
   * The widths combine a scaled Highland term with detector/reconstruction
   * resolution floors:
   *   sigma_i(p) = sqrt((scale_i(p) * sigma_Highland(p))^2 + res_i^2).
   *
   * For configuration options see TrajectoryMCSFitter::Config
   */
  class TrajectoryMCSFitter {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      // Basic fit controls: particle hypothesis, segmentation, energy-loss
      // propagation, and momentum scan range.
      fhicl::Atom<int> pIdHypothesis{
        Name("pIdHypothesis"),
        Comment("Default particle Id hypothesis to be used in the fit when not specified."),
        13};

      fhicl::Atom<int> minNumSegments{
        Name("minNumSegments"),
        Comment("Minimum number of segments the track is split into."),
        3};

      fhicl::Atom<double> segmentLength{
        Name("segmentLength"),
        Comment("Nominal length of track segments used in the fit."),
        14.};

      fhicl::Atom<int> minHitsPerSegment{
        Name("minHitsPerSegment"),
        Comment("Exclude segments with fewer hits than this value."),
        2};

      fhicl::Atom<int> nElossSteps{
        Name("nElossSteps"),
        Comment("Number of steps for computing energy loss upstream to current segment."),
        10};

      fhicl::Atom<int> eLossMode{
        Name("eLossMode"),
        Comment("Default is MPV Landau. Choose 1 for MIP (constant); 2 for Bethe-Bloch."),
        0};

      fhicl::Atom<double> pMin{
        Name("pMin"),
        Comment("Minimum momentum value in likelihood scan."),
        0.01};

      fhicl::Atom<double> pMax{
        Name("pMax"),
        Comment("Maximum momentum value in likelihood scan."),
        7.50};

      fhicl::Atom<double> pStepCoarse{
        Name("pStepCoarse"),
        Comment("Step in momentum value in initial coarse likelihood scan."),
        0.01};

      fhicl::Atom<double> pStep{
        Name("pStep"),
        Comment("Step in momentum value in fine grained likelihood scan."),
        0.01};

      fhicl::Atom<double> fineScanWindow{
        Name("fineScanWindow"),
        Comment("Window size for fine grained likelihood scan around result of coarse scan."),
        0.01};

      fhicl::Sequence<double, 5> hlParams{
        Name("hlParams"),
        Comment(
          "Parameters for tuning of Highland formula. Formula is "
          "hlParams[0]/(p*p) + hlParams[1]/p + hlParams[2] + "
          "hlParams[3]*p + hlParams[4]*p*p."),
        {0., 0., 13.6, 0., 0.}};

      fhicl::Atom<double> segLenTolerance{
        Name("segLenTolerance"),
        Comment("Tolerance in actual segment length (lower bound)."),
        1.0};

      fhicl::Atom<bool> applySCEcorr{
        Name("applySCEcorr"),
        Comment("Flag to turn the Space Charge Effect correction on/off."),
        false};

      // -----------------------------------------------------------------------
      // Calibrated double-Gaussian angle model.
      //
      // The fitter keeps the standard segmenting, projected-angle definitions,
      // and likelihood scan.  Only the per-angle PDF is calibrated.  These
      // parameters are always the active likelihood model in this cleanup copy:
      //
      //   sigma1 = sqrt((scale1(p) * Highland)^2 + res1^2)
      //   sigma2 = sqrt((scale2(p) * Highland)^2 + res2^2)
      //
      // The scale functions are quadratic polynomials in log(p/pivot), exponentiated
      // so the scales stay positive:
      //
      //   scale(p) = exp(c0 + c1*log(p/pivot) + c2*log(p/pivot)^2)
      //
      // The fraction is a smooth transition in log(momentum):
      //
      //   frac(p) = high + (low - high)/(1 + exp(slope*(log(p)-log(mid))))
      //
      // The optional direction-dependent y'-z' block below can replace the
      // global y'-z' parameters using bins in |v_x|, the drift-direction
      // component of the local segment direction.
      // -----------------------------------------------------------------------
      fhicl::Atom<double> smoothPivotX{
        Name("smoothPivotX"),
        Comment("Pivot momentum [GeV/c] used in the x'-z' scale polynomials."),
        1.76013};

      fhicl::Atom<double> smoothPivotY{
        Name("smoothPivotY"),
        Comment("Pivot momentum [GeV/c] used in the y'-z' scale polynomials."),
        1.76013};

      fhicl::Sequence<double> smoothScale1X{
        Name("smoothScale1X"),
        Comment("x'-z' primary-width scale polynomial coefficients in log(p/pivot)."),
        {-0.2606437005830629, -0.6097282621750665, -0.32643875622580065}};

      fhicl::Sequence<double> smoothScale2X{
        Name("smoothScale2X"),
        Comment("x'-z' tail-width scale polynomial coefficients in log(p/pivot)."),
        {0.28830565314289575, -2.4268569828240127, -1.2656618807916975}};

      fhicl::Sequence<double> smoothScale1Y{
        Name("smoothScale1Y"),
        Comment("y'-z' primary-width scale polynomial coefficients in log(p/pivot)."),
        {-1.701443481233165, -2.8354591664952666, -1.1096648802852689}};

      fhicl::Sequence<double> smoothScale2Y{
        Name("smoothScale2Y"),
        Comment("y'-z' tail-width scale polynomial coefficients in log(p/pivot)."),
        {2.3046181068318963, -0.014085500410843253, -0.38723571558075676}};

      fhicl::Atom<double> smoothRes1X{
        Name("smoothRes1X"),
        Comment("x'-z' core detector/reco resolution floor [mrad]."),
        4.60589};

      fhicl::Atom<double> smoothRes2X{
        Name("smoothRes2X"),
        Comment("x'-z' tail detector/reco resolution floor [mrad]."),
        17.7963};

      fhicl::Atom<double> smoothRes1Y{
        Name("smoothRes1Y"),
        Comment("y'-z' core detector/reco resolution floor [mrad]."),
        12.8517};

      fhicl::Atom<double> smoothRes2Y{
        Name("smoothRes2Y"),
        Comment("y'-z' tail detector/reco resolution floor [mrad]."),
        104.754};

      fhicl::Atom<double> smoothFracLowX{
        Name("smoothFracLowX"),
        Comment("Low-momentum x'-z' primary Gaussian fraction."),
        0.732874};

      fhicl::Atom<double> smoothFracHighX{
        Name("smoothFracHighX"),
        Comment("High-momentum x'-z' primary Gaussian fraction."),
        0.914986};

      fhicl::Atom<double> smoothFracSlopeX{
        Name("smoothFracSlopeX"),
        Comment("x'-z' log-momentum transition slope for the primary fraction."),
        34.6052};

      fhicl::Atom<double> smoothFracMidX{
        Name("smoothFracMidX"),
        Comment("x'-z' midpoint momentum [GeV/c] for the primary-fraction transition."),
        0.697278};

      fhicl::Atom<double> smoothFracLowY{
        Name("smoothFracLowY"),
        Comment("Low-momentum y'-z' primary Gaussian fraction."),
        0.87595};

      fhicl::Atom<double> smoothFracHighY{
        Name("smoothFracHighY"),
        Comment("High-momentum y'-z' primary Gaussian fraction."),
        0.800436};

      fhicl::Atom<double> smoothFracSlopeY{
        Name("smoothFracSlopeY"),
        Comment("y'-z' log-momentum transition slope for the primary fraction."),
        50.0};

      fhicl::Atom<double> smoothFracMidY{
        Name("smoothFracMidY"),
        Comment("y'-z' midpoint momentum [GeV/c] for the primary-fraction transition."),
        2.06062};

      fhicl::Atom<bool> useDirectionDependentYZ{
        Name("useDirectionDependentYZ"),
        Comment(
          "If true, choose the y'-z' smooth double-Gaussian parameters from bins in "
          "|v_x|, the absolute drift-direction component of the average adjoining "
          "segment direction. Defaults are copied from the current momentum-only "
          "YZ calibration, so enabling this flag changes physics only after "
          "direction-binned YZ parameters are supplied."),
        false};

      fhicl::Sequence<double> smoothYDirBinEdges{
        Name("smoothYDirBinEdges"),
        Comment("|v_x| bin edges for direction-dependent y'-z' tuning."),
        {0.0, 0.1, 0.2, 0.35, 0.75, 1.0}};

      fhicl::Sequence<double> smoothPivotYByDir{
        Name("smoothPivotYByDir"),
        Comment(
          "Direction-binned y'-z' pivot momenta [GeV/c]. These pivots must match "
          "the pivots used when fitting the direction-binned scale coefficients."),
        {1.76013, 1.76013, 1.76013, 1.76013, 1.76013}};

      fhicl::Sequence<double> smoothScale1YByDir{
        Name("smoothScale1YByDir"),
        Comment(
          "Flattened direction-binned y'-z' core scale coefficients. With five |v_x| "
          "bins and quadratic log-polynomials, this has 5*3 values. Defaults repeat "
          "the current momentum-only YZ coefficients in each direction bin."),
        {-1.701443481233165, -2.8354591664952666, -1.1096648802852689,
         -1.701443481233165, -2.8354591664952666, -1.1096648802852689,
         -1.701443481233165, -2.8354591664952666, -1.1096648802852689,
         -1.701443481233165, -2.8354591664952666, -1.1096648802852689,
         -1.701443481233165, -2.8354591664952666, -1.1096648802852689}};

      fhicl::Sequence<double> smoothScale2YByDir{
        Name("smoothScale2YByDir"),
        Comment(
          "Flattened direction-binned y'-z' tail scale coefficients. With five |v_x| "
          "bins and quadratic log-polynomials, this has 5*3 values. Defaults repeat "
          "the current momentum-only YZ coefficients in each direction bin."),
        {2.3046181068318963, -0.014085500410843253, -0.38723571558075676,
         2.3046181068318963, -0.014085500410843253, -0.38723571558075676,
         2.3046181068318963, -0.014085500410843253, -0.38723571558075676,
         2.3046181068318963, -0.014085500410843253, -0.38723571558075676,
         2.3046181068318963, -0.014085500410843253, -0.38723571558075676}};

      fhicl::Sequence<double> smoothRes1YByDir{
        Name("smoothRes1YByDir"),
        Comment("Direction-binned y'-z' core detector/reco resolution floors [mrad]."),
        {12.8517, 12.8517, 12.8517, 12.8517, 12.8517}};

      fhicl::Sequence<double> smoothRes2YByDir{
        Name("smoothRes2YByDir"),
        Comment("Direction-binned y'-z' tail detector/reco resolution floors [mrad]."),
        {104.754, 104.754, 104.754, 104.754, 104.754}};

      fhicl::Sequence<double> smoothFracLowYByDir{
        Name("smoothFracLowYByDir"),
        Comment("Direction-binned low-momentum y'-z' primary Gaussian fractions."),
        {0.87595, 0.87595, 0.87595, 0.87595, 0.87595}};

      fhicl::Sequence<double> smoothFracHighYByDir{
        Name("smoothFracHighYByDir"),
        Comment("Direction-binned high-momentum y'-z' primary Gaussian fractions."),
        {0.800436, 0.800436, 0.800436, 0.800436, 0.800436}};

      fhicl::Sequence<double> smoothFracSlopeYByDir{
        Name("smoothFracSlopeYByDir"),
        Comment("Direction-binned y'-z' log-momentum fraction transition slopes."),
        {50.0, 50.0, 50.0, 50.0, 50.0}};

      fhicl::Sequence<double> smoothFracMidYByDir{
        Name("smoothFracMidYByDir"),
        Comment("Direction-binned y'-z' fraction midpoint momenta [GeV/c]."),
        {2.06062, 2.06062, 2.06062, 2.06062, 2.06062}};

      fhicl::Atom<bool> useYZHighVxFallback{
        Name("useYZHighVxFallback"),
        Comment(
          "If true and direction-dependent y'-z' is enabled, blend sparse high-|v_x| "
          "y'-z' bins toward the nearest lower-|v_x| tune at high kinetic energy. "
          "This keeps the low-energy direction bins separate while avoiding "
          "unstable sparse-statistics behavior at high momentum."),
        false};

      fhicl::Sequence<double, 3> yzBlendKEEdgesMeV{
        Name("yzBlendKEEdgesMeV"),
        Comment(
          "Kinetic-energy edges [MeV] for the high-|v_x| y'-z' blending. "
          "Default is low=600, middle=950, high=1300 MeV."),
        {600.0, 950.0, 1300.0}};

      fhicl::Atom<double> yzBlendWidthMeV{
        Name("yzBlendWidthMeV"),
        Comment("Sigmoid width [MeV] for the high-|v_x| y'-z' blending."),
        50.0};
    };

    using Parameters = fhicl::Table<Config>;

    TrajectoryMCSFitter(int pIdHyp,
                        int minNSegs,
                        double segLen,
                        int minHitsPerSegment,
                        int nElossSteps,
                        int eLossMode,
                        double pMin,
                        double pMax,
                        double pStepCoarse,
                        double pStep,
                        double fineScanWindow,
                        const std::array<double, 5>& hlParams,
                        double segLenTolerance,
                        bool applySCEcorr,
                        double smoothPivotX = 1.76013,
                        double smoothPivotY = 1.76013,
                        std::vector<double> smoothScale1X = std::vector<double>{},
                        std::vector<double> smoothScale2X = std::vector<double>{},
                        std::vector<double> smoothScale1Y = std::vector<double>{},
                        std::vector<double> smoothScale2Y = std::vector<double>{},
                        double smoothRes1X = 4.60589,
                        double smoothRes2X = 17.7963,
                        double smoothRes1Y = 12.8517,
                        double smoothRes2Y = 104.754,
                        double smoothFracLowX = 0.732874,
                        double smoothFracHighX = 0.914986,
                        double smoothFracSlopeX = 34.6052,
                        double smoothFracMidX = 0.697278,
                        double smoothFracLowY = 0.87595,
                        double smoothFracHighY = 0.800436,
                        double smoothFracSlopeY = 50.0,
                        double smoothFracMidY = 2.06062,
                        bool useDirectionDependentYZ = false,
                        std::vector<double> smoothYDirBinEdges = std::vector<double>{},
                        std::vector<double> smoothPivotYByDir = std::vector<double>{},
                        std::vector<double> smoothScale1YByDir = std::vector<double>{},
                        std::vector<double> smoothScale2YByDir = std::vector<double>{},
                        std::vector<double> smoothRes1YByDir = std::vector<double>{},
                        std::vector<double> smoothRes2YByDir = std::vector<double>{},
                        std::vector<double> smoothFracLowYByDir = std::vector<double>{},
                        std::vector<double> smoothFracHighYByDir = std::vector<double>{},
                        std::vector<double> smoothFracSlopeYByDir = std::vector<double>{},
                        std::vector<double> smoothFracMidYByDir = std::vector<double>{},
                        bool useYZHighVxFallback = false,
                        const std::array<double, 3>& yzBlendKEEdgesMeV =
                          std::array<double, 3>{{600.0, 950.0, 1300.0}},
                        double yzBlendWidthMeV = 50.0)
    {
      pIdHyp_ = pIdHyp;
      minNSegs_ = minNSegs;
      segLen_ = segLen;
      minHitsPerSegment_ = minHitsPerSegment;
      nElossSteps_ = nElossSteps;
      eLossMode_ = eLossMode;
      pMin_ = pMin;
      pMax_ = pMax;
      pStepCoarse_ = pStepCoarse;
      pStep_ = pStep;
      fineScanWindow_ = fineScanWindow;
      hlParams_ = hlParams;
      segLenTolerance_ = segLenTolerance;
      applySCEcorr_ = applySCEcorr;
      smoothPivotX_ = smoothPivotX;
      smoothPivotY_ = smoothPivotY;
      smoothScale1X_ = std::move(smoothScale1X);
      smoothScale2X_ = std::move(smoothScale2X);
      smoothScale1Y_ = std::move(smoothScale1Y);
      smoothScale2Y_ = std::move(smoothScale2Y);
      smoothRes1X_ = smoothRes1X;
      smoothRes2X_ = smoothRes2X;
      smoothRes1Y_ = smoothRes1Y;
      smoothRes2Y_ = smoothRes2Y;
      smoothFracLowX_ = smoothFracLowX;
      smoothFracHighX_ = smoothFracHighX;
      smoothFracSlopeX_ = smoothFracSlopeX;
      smoothFracMidX_ = smoothFracMidX;
      smoothFracLowY_ = smoothFracLowY;
      smoothFracHighY_ = smoothFracHighY;
      smoothFracSlopeY_ = smoothFracSlopeY;
      smoothFracMidY_ = smoothFracMidY;
      useDirectionDependentYZ_ = useDirectionDependentYZ;
      smoothYDirBinEdges_ = std::move(smoothYDirBinEdges);
      smoothPivotYByDir_ = std::move(smoothPivotYByDir);
      smoothScale1YByDir_ = std::move(smoothScale1YByDir);
      smoothScale2YByDir_ = std::move(smoothScale2YByDir);
      smoothRes1YByDir_ = std::move(smoothRes1YByDir);
      smoothRes2YByDir_ = std::move(smoothRes2YByDir);
      smoothFracLowYByDir_ = std::move(smoothFracLowYByDir);
      smoothFracHighYByDir_ = std::move(smoothFracHighYByDir);
      smoothFracSlopeYByDir_ = std::move(smoothFracSlopeYByDir);
      smoothFracMidYByDir_ = std::move(smoothFracMidYByDir);
      useYZHighVxFallback_ = useYZHighVxFallback;
      yzBlendKEEdgesMeV_ = yzBlendKEEdgesMeV;
      yzBlendWidthMeV_ = yzBlendWidthMeV;
    }

    explicit TrajectoryMCSFitter(const Parameters& p)
      : TrajectoryMCSFitter(p().pIdHypothesis(),
                            p().minNumSegments(),
                            p().segmentLength(),
                            p().minHitsPerSegment(),
                            p().nElossSteps(),
                            p().eLossMode(),
                            p().pMin(),
                            p().pMax(),
                            p().pStepCoarse(),
                            p().pStep(),
                            p().fineScanWindow(),
                            p().hlParams(),
                            p().segLenTolerance(),
                            p().applySCEcorr(),
                            p().smoothPivotX(),
                            p().smoothPivotY(),
                            p().smoothScale1X(),
                            p().smoothScale2X(),
                            p().smoothScale1Y(),
                            p().smoothScale2Y(),
                            p().smoothRes1X(),
                            p().smoothRes2X(),
                            p().smoothRes1Y(),
                            p().smoothRes2Y(),
                            p().smoothFracLowX(),
                            p().smoothFracHighX(),
                            p().smoothFracSlopeX(),
                            p().smoothFracMidX(),
                            p().smoothFracLowY(),
                            p().smoothFracHighY(),
                            p().smoothFracSlopeY(),
                            p().smoothFracMidY(),
                            p().useDirectionDependentYZ(),
                            p().smoothYDirBinEdges(),
                            p().smoothPivotYByDir(),
                            p().smoothScale1YByDir(),
                            p().smoothScale2YByDir(),
                            p().smoothRes1YByDir(),
                            p().smoothRes2YByDir(),
                            p().smoothFracLowYByDir(),
                            p().smoothFracHighYByDir(),
                            p().smoothFracSlopeYByDir(),
                            p().smoothFracMidYByDir(),
                            p().useYZHighVxFallback(),
                            p().yzBlendKEEdgesMeV(),
                            p().yzBlendWidthMeV())
    {}

    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj) const
    {
      return fitMcs(traj, pIdHyp_);
    }

    recob::MCSFitResult fitMcs(const recob::Track& track) const { return fitMcs(track, pIdHyp_); }

    recob::MCSFitResult fitMcs(const recob::Trajectory& traj) const
    {
      return fitMcs(traj, pIdHyp_);
    }

    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj, int pid) const;

    recob::MCSFitResult fitMcs(const recob::Track& track, int pid) const
    {
      return fitMcs(track.Trajectory(), pid);
    }

    recob::MCSFitResult fitMcs(const recob::Trajectory& traj, int pid) const
    {
      recob::TrackTrajectory::Flags_t flags(traj.NPoints());
      const recob::TrackTrajectory tt(traj, std::move(flags));
      return fitMcs(tt, pid);
    }

    void breakTrajInSegments(const recob::TrackTrajectory& traj,
                             std::vector<size_t>& breakpoints,
                             std::vector<float>& segradlengths,
                             std::vector<float>& cumseglens) const;

    void linearRegression(const recob::TrackTrajectory& traj,
                          const size_t firstPoint,
                          const size_t lastPoint,
                          recob::tracking::Vector_t& pcdir) const;

    double mcsLikelihood(double p,
                         std::vector<float>& dthetaX,
                         std::vector<float>& dthetaY,
                         std::vector<float>& angleDriftFrac,
                         std::vector<float>& seg_nradl,
                         std::vector<float>& cumLen,
                         bool fwd,
                         int pid) const;

    struct ScanResult {
      ScanResult(double ap, double apUnc, double alogL) : p(ap), pUnc(apUnc), logL(alogL) {}
      double p, pUnc, logL;
    };

    // Fine + coarse scan: uses class-level pMin_, pMax_, pStepCoarse_, pStep_
    const ScanResult doLikelihoodScan(std::vector<float>& dthetaX,
                                      std::vector<float>& dthetaY,
                                      std::vector<float>& angleDriftFrac,
                                      std::vector<float>& seg_nradlengths,
                                      std::vector<float>& cumLen,
                                      bool fwdFit,
                                      int pid) const;

    // Raw scan over [pmin, pmax] with given pstep
    const ScanResult doLikelihoodScan(std::vector<float>& dthetaX,
                                      std::vector<float>& dthetaY,
                                      std::vector<float>& angleDriftFrac,
                                      std::vector<float>& seg_nradlengths,
                                      std::vector<float>& cumLen,
                                      bool fwdFit,
                                      int pid,
                                      float pmin,
                                      float pmax,
                                      float pstep) const;

    inline double HighlandFirstTerm(const double p) const
    {
      return hlParams_[0] / (p * p) + hlParams_[1] / p + hlParams_[2] + hlParams_[3] * p +
             hlParams_[4] * p * p;
    }

    inline double SmoothScale(const double p,
                              const std::vector<double>& coeffs,
                              const double pivot) const
    {
      // Positive scale multiplier for the Highland width.  The polynomial is
      // evaluated in log(p/pivot), then exponentiated so the scale cannot go
      // negative during the likelihood scan.
      if (coeffs.empty()) return 1.0;

      const double safePivot = std::max(pivot, 1e-6);
      const double x = std::log(std::max(p, 1e-6) / safePivot);
      double poly = 0.0;
      double xpow = 1.0;
      for (const double c : coeffs) {
        poly += c * xpow;
        xpow *= x;
      }
      return std::exp(poly);
    }

    inline double SmoothFrac(const double p,
                             const double low,
                             const double high,
                             const double slope,
                             const double mid) const
    {
      // Smooth transition for the primary Gaussian area fraction.  The low and
      // high parameters are the asymptotic fractions; slope and mid control
      // where and how quickly the transition happens in log(momentum).
      const double safeMid = std::max(mid, 1e-6);
      const double arg = slope * (std::log(std::max(p, 1e-6)) - std::log(safeMid));
      const double frac = high + (low - high) / (1.0 + std::exp(arg));
      return std::clamp(frac, 1e-6, 1.0 - 1e-6);
    }

    inline size_t SmoothYDirectionBin(const double driftFrac) const
    {
      // Choose which direction-dependent y'-z' tune to use from |v_x|, the
      // absolute detector-drift component of the local segment-pair direction.
      if (smoothYDirBinEdges_.size() < 2) return 0;

      const double v = std::clamp(std::abs(driftFrac), 0.0, 1.0);
      for (size_t i = 0; i + 1 < smoothYDirBinEdges_.size(); ++i) {
        if (v >= smoothYDirBinEdges_[i] && v < smoothYDirBinEdges_[i + 1]) return i;
      }
      return smoothYDirBinEdges_.size() - 2;
    }

    inline double DirectionValue(const std::vector<double>& values,
                                 const size_t bin,
                                 const double fallback) const
    {
      // Direction-binned scalar lookup with a global fallback if the configured
      // vector is missing a bin.
      return (bin < values.size() ? values[bin] : fallback);
    }

    inline double DirectionSmoothScale(const double p,
                                       const std::vector<double>& flatCoeffs,
                                       const std::vector<double>& fallbackCoeffs,
                                       const double pivot,
                                       const size_t bin) const
    {
      // Direction-binned version of SmoothScale.  The coefficients are stored
      // as a flat vector: [bin0 c0,c1,c2, bin1 c0,c1,c2, ...].
      const size_t nBins = (smoothYDirBinEdges_.size() >= 2 ? smoothYDirBinEdges_.size() - 1 : 0);
      if (nBins == 0 || flatCoeffs.empty() || flatCoeffs.size() % nBins != 0) {
        return SmoothScale(p, fallbackCoeffs, pivot);
      }

      const size_t coeffsPerBin = flatCoeffs.size() / nBins;
      if (coeffsPerBin == 0 || bin >= nBins) return SmoothScale(p, fallbackCoeffs, pivot);

      const double safePivot = std::max(pivot, 1e-6);
      const double x = std::log(std::max(p, 1e-6) / safePivot);
      double poly = 0.0;
      double xpow = 1.0;
      const size_t offset = bin * coeffsPerBin;
      for (size_t j = 0; j < coeffsPerBin; ++j) {
        poly += flatCoeffs[offset + j] * xpow;
        xpow *= x;
      }
      return std::exp(poly);
    }

    inline double YZBlendWeight(const double kineticEnergyMeV,
                                const double edgeMeV) const
    {
      // Smooth 0-to-1 blend used when high-|v_x| y'-z' bins are too sparse at
      // high energy.  The width keeps the transition continuous in the PDF.
      const double width = std::max(yzBlendWidthMeV_, 1e-6);
      const double arg = std::clamp((kineticEnergyMeV - edgeMeV) / width, -60.0, 60.0);
      return 1.0 / (1.0 + std::exp(-arg));
    }

    double mass(int pid) const
    {
      if (abs(pid) == 13) return mumass;
      if (abs(pid) == 211) return pimass;
      if (abs(pid) == 321) return kmass;
      if (abs(pid) == 2212) return pmass;
      return util::kBogusD;
    }

    double energyLossBetheBloch(const double mass, const double e2) const;
    double energyLossLandau(const double mass2, const double E2, const double x) const;
    double GetE(const double initial_E, const double length_travelled, const double mass) const;

    int minNSegs() const { return minNSegs_; }
    double segLen() const { return segLen_; }
    double segLenTolerance() const { return segLenTolerance_; }

  private:
    int pIdHyp_;
    int minNSegs_;
    double segLen_;
    int minHitsPerSegment_;
    int nElossSteps_;
    int eLossMode_;
    double pMin_;
    double pMax_;
    double pStepCoarse_;
    double pStep_;
    double fineScanWindow_;
    std::array<double, 5> hlParams_;
    double segLenTolerance_;
    bool applySCEcorr_;
    // Global calibrated two-Gaussian parameters.
    double smoothPivotX_;
    double smoothPivotY_;
    std::vector<double> smoothScale1X_;
    std::vector<double> smoothScale2X_;
    std::vector<double> smoothScale1Y_;
    std::vector<double> smoothScale2Y_;
    double smoothRes1X_;
    double smoothRes2X_;
    double smoothRes1Y_;
    double smoothRes2Y_;
    double smoothFracLowX_;
    double smoothFracHighX_;
    double smoothFracSlopeX_;
    double smoothFracMidX_;
    double smoothFracLowY_;
    double smoothFracHighY_;
    double smoothFracSlopeY_;
    double smoothFracMidY_;
    bool useDirectionDependentYZ_;
    std::vector<double> smoothYDirBinEdges_;
    std::vector<double> smoothPivotYByDir_;
    std::vector<double> smoothScale1YByDir_;
    std::vector<double> smoothScale2YByDir_;
    std::vector<double> smoothRes1YByDir_;
    std::vector<double> smoothRes2YByDir_;
    std::vector<double> smoothFracLowYByDir_;
    std::vector<double> smoothFracHighYByDir_;
    std::vector<double> smoothFracSlopeYByDir_;
    std::vector<double> smoothFracMidYByDir_;
    bool useYZHighVxFallback_;
    std::array<double, 3> yzBlendKEEdgesMeV_;
    double yzBlendWidthMeV_;
  };
}

#endif

