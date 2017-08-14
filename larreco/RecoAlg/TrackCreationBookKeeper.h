#ifndef TRACKCREATIONBOOKKEEPER_H
#define TRACKCREATIONBOOKKEEPER_H

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "larreco/TrackFinder/TrackMaker.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"

namespace trkmkr {

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
      //else fixme: error message?
    }
    //
    void setTotChi2(double totChi2) { totChi2_ = totChi2; }
    //
    void rejectPoint(unsigned int pos) {
      if (pos>=positions.size()) return;//fixme: error message?
      // move elements to the end of their vectors
      std::rotate( (positions.begin()+pos),(positions.begin()+pos+1),positions.end() );
      std::rotate( (momenta.begin()+pos),(momenta.begin()+pos+1),momenta.end() );
      std::rotate( (flags.begin()+pos),(flags.begin()+pos+1),flags.end() );
      std::rotate( (hits->begin()+pos),(hits->begin()+pos+1),hits->end() );
      // overwrite values
      positions.back() = Point_t(util::kBogusD,util::kBogusD,util::kBogusD);
      momenta.back() = Vector_t(util::kBogusD,util::kBogusD,util::kBogusD);
      auto mask = flags.back().mask();
      auto fhit = flags.back().fromHit();
      mask.set(recob::TrajectoryPointFlagTraits::HitIgnored,recob::TrajectoryPointFlagTraits::NoPoint);
      if (mask.isSet(recob::TrajectoryPointFlagTraits::Rejected)==0) mask.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
      flags.back() = recob::TrajectoryPointFlags(fhit,mask);
      // do the same for optional outputs
      if (opts->isTrackFitInfosInit()) {
	std::rotate( (opts->trackFitHitInfos()->begin()+pos),(opts->trackFitHitInfos()->begin()+pos+1),opts->trackFitHitInfos()->end() );
	SVector5 fakePar5(util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD);
	SMatrixSym55 fakeCov55;
	for (int i=0;i<5;i++) for (int j=i;j<5;j++) fakeCov55(i,j) = util::kBogusD;
	opts->trackFitHitInfos()->back() = recob::TrackFitHitInfo(opts->trackFitHitInfos()->back().hitMeas(),
								  opts->trackFitHitInfos()->back().hitMeasErr2(),
								  fakePar5,fakeCov55,opts->trackFitHitInfos()->back().WireId());
      }
    }
    Point_t  posAt(unsigned int pos) const { return positions[pos]; }
    Vector_t momAt(unsigned int pos) const { return momenta[pos]; }
    unsigned int size() const { return positions.size(); }
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
