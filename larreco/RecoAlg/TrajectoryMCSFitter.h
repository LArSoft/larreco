#ifndef TRAJECTORYMCSFITTER_H
#define TRAJECTORYMCSFITTER_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/RecoObjects/TrackState.h"

namespace trkf {
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
	Comment("Particle Id Hypothesis to be used in the fit."),
	13
      };
      fhicl::Atom<int> minNumSegments {
        Name("minNumSegments"),
	Comment("Minimum number of segments the track is split into."),
	3
      };
      fhicl::Atom<double> segmentLength {
        Name("segmentLength"),
	Comment("Nominal length of track segments used in the fit."),
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
    TrajectoryMCSFitter(int pIdHyp, int minNSegs, double segLen, double pMin, double pMax, double pStep, double angResol){
      pIdHyp_ = pIdHyp;
      minNSegs_ = minNSegs;
      segLen_ = segLen;
      pMin_ = pMin;
      pMax_ = pMax;
      pStep_ = pStep;
      angResol_ = angResol;
    }
    explicit TrajectoryMCSFitter(const Parameters & p)
      : TrajectoryMCSFitter(p().pIdHypothesis(),p().minNumSegments(),p().segmentLength(),p().pMin(),p().pMax(),p().pStep(),p().angResol()) {}
    //
    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj, bool momDepConst = true) const;
    recob::MCSFitResult fitMcs(const recob::Track& track, bool momDepConst = true) const { return fitMcs(track.Trajectory(),momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Trajectory& traj, bool momDepConst = true) const {
      recob::TrackTrajectory::Flags_t flags(traj.NPoints());
      const recob::TrackTrajectory tt(traj,std::move(flags));
      return fitMcs(tt,momDepConst);
    }
    //
    void linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, recob::tracking::Vector_t& pcdir) const;
    double mcsLikelihood(double p, double theta0x, std::vector<double>& dthetaij, std::vector<double>& seg_nradl, std::vector<double>& cumLen, bool fwd, bool momDepConst) const;
    struct ScanResult {
    public:
    ScanResult(double ap, double apUnc, double alogL) : p(ap), pUnc(apUnc), logL(alogL) {}
      double p, pUnc, logL;
    };
    const ScanResult doLikelihoodScan(std::vector<double>& dtheta, std::vector<double>& seg_nradlengths, std::vector<double>& cumLen, bool fwdFit, bool momDepConst) const;
    //
    inline double MomentumDependentConstant(const double p) const {
      constexpr double a = 0.1049;
      constexpr double c = 11.0038;
      return (a/(p*p)) + c;
    }
    double mass() const {
      if (abs(pIdHyp_)==13)   { return mumass; }
      if (abs(pIdHyp_)==211)  { return pimass; }
      if (abs(pIdHyp_)==321)  { return kmass;  }
      if (abs(pIdHyp_)==2212) { return pmass;  }
      return util::kBogusD;
    }
    double energyLossBetheBloch(const double mass,const double p) const;//fixme: remove if worse than landau
    double energyLossLandau(const double mass,const double p, const double x) const;
    //
  private:
    int pIdHyp_;
    int    minNSegs_;
    double segLen_;
    double pMin_;
    double pMax_;
    double pStep_;
    double angResol_;
  };
}

#endif
