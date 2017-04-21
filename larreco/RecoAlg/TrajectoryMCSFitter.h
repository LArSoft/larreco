#ifndef TRAJECTORYMCSFITTER_H
#define TRAJECTORYMCSFITTER_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"

namespace trkf {
  namespace {
    constexpr double mumass2 = 0.105658367*0.105658367;// Muon
    constexpr double pimass2 = 0.13957*0.13957;        // Charged pion
    constexpr double kmass2  = 0.493677*0.493677;      // Charged kaon
    constexpr double pmass2  = 0.938272*0.938272;      // Proton
  }
  /**
   * @brief something
   *
   * Something
   *
   * Inputs are: 
   * Output are: 
   *
   */
  class TrajectoryMCSFitter {
    // 
  public:
    //
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> pIdHypothesis {
        Name("pIdHypothesis"),
	Comment("Particle Id Hypothesis to be used in the fit.")
      };
      fhicl::Atom<double> segmentLength {
        Name("segmentLength"),
	Comment("Target length of track segments used in the fit."),
	14.
      };
      fhicl::Atom<double> pMin {
        Name("pMin"),
	Comment("Minimum momentum value in likelihood scan."),
	0.01
      };
      fhicl::Atom<double> pMax {
        Name("pMax"),
	Comment("Maximum momentum value in likelihood scan."),
	7.50  
      };
      fhicl::Atom<double> pStep {
        Name("pStep"),
	Comment("Step in momentum value in likelihood scan."),
	0.01
      };
      fhicl::Atom<double> angResol {
        Name("angResol"),
	Comment("Angular resolution parameter used in modified Highland formula. Unit is mrad."),
	3.0
      };
    };
    using Parameters = fhicl::Table<Config>;
    //
    TrajectoryMCSFitter(int pIdHyp, double segLen, double pMin, double pMax, double pStep, double angResol){
      pIdHyp_ = pIdHyp;
      segLen_ = segLen;
      pMin_ = pMin;
      pMax_ = pMax;
      pStep_ = pStep;
      angResol_ = angResol;
    }
    explicit TrajectoryMCSFitter(const Parameters & p)
      : TrajectoryMCSFitter(p().pIdHypothesis(),p().segmentLength(),p().pMin(),p().pMax(),p().pStep(),p().angResol()) {}
    //
    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj) const;
    recob::MCSFitResult fitMcs(const recob::Track& track) const { return fitMcs(track.Trajectory()); }
    //
    void linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, recob::tracking::Vector_t& pcdir) const;
    double mcsLikelihood(double p, double theta0x, std::vector<double>& dthetaij, std::vector<double>& seg_nradl, std::vector<double>& eLoss2, bool fwd, bool momDepConst) const;
    struct ScanResult {
    public:
    ScanResult(double ap, double apUnc, double alogL) : p(ap), pUnc(apUnc), logL(alogL) {}
      double p, pUnc, logL;
    };
    const ScanResult doLikelihoodScan(std::vector<double>& dtheta, std::vector<double>& seg_nradlengths, std::vector<double>& energyLoss2, bool fwdFit) const;
    //
    inline double MomentumDependentConstant(const double p) const {
      constexpr double a = 0.1049;
      constexpr double c = 11.0038;
      return (a/(p*p)) + c;
    }
    double mass2() const {
      if (abs(pIdHyp_)==13)   { return mumass2; } 
      if (abs(pIdHyp_)==211)  { return pimass2; }
      if (abs(pIdHyp_)==321)  { return kmass2;  } 
      if (abs(pIdHyp_)==2212) { return pmass2;  }
      return util::kBogusD;
    }
    //
  private:
    int pIdHyp_;
    double segLen_;
    double pMin_;
    double pMax_;
    double pStep_;
    double angResol_;
  };
}

#endif
