#ifndef TRACKCREATIONBOOKKEEPER_H
#define TRACKCREATIONBOOKKEEPER_H

////////////////////////////////////////////////////////////////////////
// Class:       TrackCreationBookKeeper
// File:        TrackCreationBookKeeper.h
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
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
  TrackCreationBookKeeper(Track& track, std::vector<art::Ptr<Hit> >& outhits, OptionalOutputs& optionals, int tkID, int pdgHyp, bool hasMomenta)
    : tkID_(tkID), pdgHyp_(pdgHyp), hasMomenta_(hasMomenta), totChi2_(0), outTrack(&track), hits(&outhits), opts(&optionals)
      {
	hits->clear();
	if (opts->isTrackFitInfosInit()) opts->trackFitHitInfos()->clear();
      }
    //
    void addPoint(const Point_t& point, const Vector_t& vect, art::Ptr<Hit> hit, const PointFlags_t& flag, double chi2) {
	positions.push_back(std::move(point));
	momenta.push_back(std::move(vect));
	hits->push_back(hit);
	flags.push_back(std::move(flag));
	chi2v.push_back(chi2);
	totChi2_+=chi2;
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<Hit> hit, PointFlags_t&& flag, double chi2) {
	positions.push_back(std::move(point));
	momenta.push_back(std::move(vect));
	hits->push_back(hit);
	flags.push_back(std::move(flag));
	chi2v.push_back(chi2);
	totChi2_+=chi2;
    }
    void addPoint(Point_t&& point, Vector_t&& vect, art::Ptr<Hit> hit, PointFlags_t&& flag, double chi2, TrackFitHitInfo&& tfhi) {
      addPoint(std::move(point), std::move(vect), hit, std::move(flag), chi2);
      if (opts->isTrackFitInfosInit()) opts->trackFitHitInfos()->push_back(tfhi);
      else {
	throw cet::exception("TrackCreationBookKeeper")
	  << "Passing TrackFitHitInfo elements to addPoint, but TrackFitHitInfo collection not initialized.\n";
      }
    }
    //
    void setTotChi2(double totChi2) { totChi2_ = totChi2; }
    //
    bool hasZeroMomenta() {
      for (const auto& mom : momenta) {
	if (mom.Mag2() <= 1.0e-9) return true;
      }
      return false;
    }
    //
    void finalizeTrack(tracking::SMatrixSym55& covStart, tracking::SMatrixSym55& covEnd) {
      *outTrack = Track(std::move(positions),std::move(momenta),std::move(flags),
			hasMomenta_,pdgHyp_,totChi2_,hits->size()-fittedpars,
			std::move(covStart),std::move(covEnd),tkID_);
    }
    void finalizeTrack(tracking::SMatrixSym55&& covStart, tracking::SMatrixSym55&& covEnd) {
	*outTrack = Track(std::move(positions),std::move(momenta),std::move(flags),
			  hasMomenta_,pdgHyp_,totChi2_,hits->size()-fittedpars,
			  std::move(covStart),std::move(covEnd),tkID_);
    }
  private:
    int tkID_;
    int pdgHyp_;
    bool hasMomenta_;
    double totChi2_;
    Track* outTrack;
    std::vector<art::Ptr<Hit> >* hits;
    OptionalOutputs*             opts;
    std::vector<Point_t>         positions;
    std::vector<Vector_t>        momenta;
    std::vector<PointFlags_t>    flags;
    std::vector<double>          chi2v;
    constexpr static int fittedpars = 4;//hits are 1D measurement, i.e. each hit is one d.o.f.; no B field: 4 fitted parameters
  };

}
#endif
