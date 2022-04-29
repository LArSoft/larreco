#ifndef TRACKCREATIONBOOKKEEPER_H
#define TRACKCREATIONBOOKKEEPER_H

#include "larreco/RecoAlg/TrackTrajectoryCreationBookKeeper.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "larreco/TrackFinder/TrackMaker.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

namespace trkmkr {

  /**
   * @file  larreco/RecoAlg/TrackCreationBookKeeper.h
   * @class trkmkr::TrackCreationBookKeeper
   *
   * @brief Helper class to aid the creation of a recob::Track, keeping data vectors in sync.
   *
   * Helper class to aid the creation of a recob::Track, keeping output data (vectors of recob::tracking::Point_t, recob::tracking::Vector_t, recob::TrackTrajectory::PointFlags_t, recob::Hit, and trkmkr::OptionalOutputs struct) in sync.
   * It internally stores and uses a trkmkr::TrackTrajectoryCreationBookKeeper object. Elements of those vectors are added sequentially using the addPoint functions.
   * Once all points have been added a call to the function finalizeTrack, builds the track moving the content of the vectors.
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */

  using Point_t = recob::tracking::Point_t;
  using Vector_t = recob::tracking::Vector_t;
  using PointFlags_t = recob::TrackTrajectory::PointFlags_t;

  class TrackCreationBookKeeper {
  public:
    /// Constructor: needs reference to output hit vector, optional outputs struct, and other parameters needed when creating the track object
  TrackCreationBookKeeper(std::vector<art::Ptr<recob::Hit> >& outhits, OptionalOutputs& optionals, int tkID, int pdgHyp, bool hasMomenta, int nfitpars = 4)
    : ttcbk_(outhits, hasMomenta), tkID_(tkID), pdgHyp_(pdgHyp), totChi2_(0), opts(&optionals), nfittedpars(nfitpars)
      {
	opts->reset();
      }
    //
    //@{
    /// Avoid copies of this object
    TrackCreationBookKeeper(const TrackCreationBookKeeper&) = delete;
    TrackCreationBookKeeper(TrackCreationBookKeeper&&) = delete;
    TrackCreationBookKeeper& operator=(const TrackCreationBookKeeper&) = delete;
    TrackCreationBookKeeper& operator=(TrackCreationBookKeeper&& ) = delete;
    //@}
    //
    //@{
    /// Add a single point; different version of the functions are provided using const references or rvalue references, with and without an OptionalPointElement argument.
    void addPoint(const Point_t& point, const Vector_t& vect, art::Ptr<recob::Hit> hit, const PointFlags_t& flag, double chi2) {
      ttcbk_.addPoint(point, vect, hit, flag);
      if (chi2>=0) {
        chi2v.push_back(chi2);
        totChi2_+=chi2;
      }
    }
    void addPoint(const Point_t& point, const Vector_t& vect, art::Ptr<recob::Hit> hit, const PointFlags_t& flag, double chi2, OptionalPointElement& ope) {
      addPoint(point, vect, hit, flag, chi2);
      opts->addPoint(ope);
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<recob::Hit> hit, PointFlags_t&& flag, double chi2) {
      ttcbk_.addPoint(std::move(point), std::move(vect), hit, std::move(flag));
      if (chi2>=0) {
        chi2v.push_back(chi2);
        totChi2_+=chi2;
      }
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<recob::Hit> hit, PointFlags_t&& flag, double chi2, OptionalPointElement& ope) {
      addPoint(std::move(point), std::move(vect), hit, std::move(flag), chi2);
      opts->addPoint(ope);
    }
    //@}
    //
    /// Set the total chi2 value
    void setTotChi2(double totChi2) { totChi2_ = totChi2; }
    //
    //@{
    /// Get the finalized recob::Track; needs the start and end covariance matrices.
    recob::Track finalizeTrack(const recob::tracking::SMatrixSym55& covStart, const recob::tracking::SMatrixSym55& covEnd) {
      return recob::Track(ttcbk_.finalizeTrackTrajectory(),pdgHyp_,totChi2_,int(chi2v.size())-nfittedpars,
		   recob::tracking::SMatrixSym55(covStart),recob::tracking::SMatrixSym55(covEnd),tkID_);
    }
    recob::Track finalizeTrack(recob::tracking::SMatrixSym55&& covStart, recob::tracking::SMatrixSym55&& covEnd) {
      return recob::Track(ttcbk_.finalizeTrackTrajectory(),pdgHyp_,totChi2_,int(chi2v.size())-nfittedpars,
		   std::move(covStart),std::move(covEnd),tkID_);
    }
    //@}
  private:
    trkmkr::TrackTrajectoryCreationBookKeeper ttcbk_;
    int tkID_;
    int pdgHyp_;
    double totChi2_;
    OptionalOutputs*             opts;
    std::vector<double>          chi2v;
    int nfittedpars; // hits are 1D measurement, i.e. each hit is one d.o.f.; no B field: 4 fitted parameters by default
    //
  };

}
#endif
