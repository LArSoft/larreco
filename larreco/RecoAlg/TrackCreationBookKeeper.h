#ifndef TRACKCREATIONBOOKKEEPER_H
#define TRACKCREATIONBOOKKEEPER_H

////////////////////////////////////////////////////////////////////////
// Class:       TrackCreationBookKeeper
// File:        TrackCreationBookKeeper.h
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrackTrajectoryCreationBookKeeper.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/TrackFinder/TrackMaker.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

namespace trkmkr {

  /**
   * @brief Helper class to aid the creation of a recob::Track, keeping data vectors in sync.
   *
   * Helper class to aid the creation of a recob::Track, keeping data vectors (Point_t, Vector_t, PointFlags_t, Hit, TrackFitHitInfo) in sync.
   * Elements of those vectors are added sequentially using the addPoint functions.
   * Once all points have been added a call to the function finalizeTrack, builds the track moving the content of the vectors.
   */

  using namespace recob;
  using Point_t = tracking::Point_t;
  using Vector_t = tracking::Vector_t;
  using PointFlags_t = TrackTrajectory::PointFlags_t;

  class TrackCreationBookKeeper {
  public:
  TrackCreationBookKeeper(std::vector<art::Ptr<Hit> >& outhits, OptionalOutputs& optionals, int tkID, int pdgHyp, bool hasMomenta, int nfitpars = 4)
    : ttcbk_(outhits, hasMomenta), tkID_(tkID), pdgHyp_(pdgHyp), totChi2_(0), opts(&optionals), nfittedpars(nfitpars)
      {
	opts->reset();
      }
    //
    // avoid copies of this object
    TrackCreationBookKeeper(const TrackCreationBookKeeper&) = delete;
    TrackCreationBookKeeper(TrackCreationBookKeeper&&) = delete;
    TrackCreationBookKeeper& operator=(const TrackCreationBookKeeper&) = delete;
    TrackCreationBookKeeper& operator=(TrackCreationBookKeeper&& ) = delete;
    //
    void addPoint(const Point_t& point, const Vector_t& vect, art::Ptr<Hit> hit, const PointFlags_t& flag, double chi2) {
      ttcbk_.addPoint(point, vect, hit, flag);
      chi2v.push_back(chi2);
      totChi2_+=chi2;
    }
    void addPoint(const Point_t& point, const Vector_t& vect, art::Ptr<Hit> hit, const PointFlags_t& flag, double chi2, OptionalPointElement& ope) {
      addPoint(point, vect, hit, flag, chi2);
      opts->addPoint(ope);
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<Hit> hit, PointFlags_t&& flag, double chi2) {
      ttcbk_.addPoint(std::move(point), std::move(vect), hit, std::move(flag));
      chi2v.push_back(chi2);
      totChi2_+=chi2;
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<Hit> hit, PointFlags_t&& flag, double chi2, OptionalPointElement& ope) {
      addPoint(std::move(point), std::move(vect), hit, std::move(flag), chi2);
      opts->addPoint(ope);
    }
    //
    void setTotChi2(double totChi2) { totChi2_ = totChi2; }
    //
    Track finalizeTrack(const tracking::SMatrixSym55& covStart, const tracking::SMatrixSym55& covEnd) {
      return Track(ttcbk_.finalizeTrackTrajectory(),pdgHyp_,totChi2_,int(chi2v.size())-nfittedpars,
		   tracking::SMatrixSym55(covStart),tracking::SMatrixSym55(covEnd),tkID_);
    }
    Track finalizeTrack(tracking::SMatrixSym55&& covStart, tracking::SMatrixSym55&& covEnd) {
      return Track(ttcbk_.finalizeTrackTrajectory(),pdgHyp_,totChi2_,int(chi2v.size())-nfittedpars,
		   std::move(covStart),std::move(covEnd),tkID_);
    }
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
