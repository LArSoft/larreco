///////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterAlg.cxx
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/KalmanFilterAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"

// Local functions.

namespace {

  void update_momentum(const trkf::KVector<2>::type& defl,
		       const trkf::KSymMatrix<2>::type& errc,
		       const trkf::KSymMatrix<2>::type& errn,
		       double mass, double& invp, double& var_invp)
  // Momentum updater.
  //
  // Arguments: defl - Deflection (2D slope residual between two
  //                   sampled points on track).
  //            errc - Error matrix of deflection residual 
  //                   exclusive of multiple scattering.
  //            errn - Part of deflection noise error matrix
  //                   proportional to 1/(beta*momentum).
  //            mass - Mass of particle.
  //            invp - Modified 1/momentum track parameter.
  //            var_invp - Modified 1/momentum variance.
  //
  // Returns: True if success.
  //
  // The total deflection residual error matrix is 
  //
  // err = errc + [1/(beta*momentum)] * errn.
  //
  // The inverse momentum and error are updated using using log-likelihood
  //
  // log(L) = -0.5 * log(det(err)) - 0.5 * defl^T * err^(-1) * defl.
  // log(L) = -0.5 * tr(log(err)) - 0.5 * defl^T * err^(-1) * defl.
  //
  {
    // Calculate original k = 1./(beta*p), and original variance of k.

    double invp2 = invp*invp;
    double invp3 = invp*invp2;
    double invp4 = invp2*invp2;
    double mass2 = mass*mass;
    double k = std::sqrt(invp2 + mass2 * invp4);
    double dkdinvp = (invp + 2.*mass2*invp3) / k;
    double vark = var_invp * dkdinvp*dkdinvp;

    // First, find current inverse error matrix using momentum hypothesis.

    trkf::KSymMatrix<2>::type inverr = errc + k * errn;
    trkf::syminvert(inverr);

    // Find the first and second derivatives of the log likelihood
    // with respact to k.

    trkf::KMatrix<2, 2>::type temp1 = prod(inverr, errn);
    trkf::KMatrix<2, 2>::type temp2 = prod(temp1, temp1);

    trkf::KVector<2>::type vtemp1 = prod(inverr, defl);
    trkf::KVector<2>::type vtemp2 = prod(temp1, vtemp1);
    trkf::KVector<2>::type vtemp3 = prod(temp1, vtemp2);
    double derivk1 = -0.5 * trkf::trace(temp1) + 0.5 * inner_prod(defl, vtemp2);
    double derivk2 = 0.5 * trkf::trace(temp2) - inner_prod(defl, vtemp3);

    // We expect the log-likelihood to be most nearly Gaussian
    // with respect to variable q = k^(-1/2) = std::sqrt(beta*p).
    // Therefore, transform the original variables and log-likelihood
    // derivatives to q-space.

    double q = 1./std::sqrt(k);
    double varq = vark / (4.*k*k*k);
    double derivq1 = (-2.*k/q) * derivk1;
    double derivq2 = 6.*k*k * derivk1 + 4.*k*k*k * derivk2;

    if(derivq2 < 0.) {

      // Estimate the measurement in q-space.

      double q1 = q - derivq1 / derivq2;
      double varq1 = -1./derivq2;

      // Get updated estimated q, combining original estimate and
      // update from the current measurement.

      double newvarq = 1. / (1./varq + 1./varq1);
      double newq = newvarq * (q/varq + q1/varq1);
      q = newq;
      varq = newvarq;
      
      // Calculate updated c = 1./p and variance.

      double q2 = q*q;
      double q4 = q2*q2;
      double c2 = 2. / (q4 + q2 * std::sqrt(q4 + 4.*mass2));
      double c = std::sqrt(c2);
      double dcdq = -2. * (c/q) * (1. + mass2*c2) / (1. + 2.*mass2*c2);
      double varc = varq * dcdq*dcdq;

      // Update result.

      invp = c;
      var_invp = varc;
    }
  }
}

/// Constructor.
  
trkf::KalmanFilterAlg::KalmanFilterAlg(const fhicl::ParameterSet& pset) :
  fTrace(false),
  fMaxPErr(0.),
  fGoodPErr(0.),
  fMaxIncChisq(0.),
  fMaxEndChisq(0.),
  fMinLHits(0),
  fMaxLDist(0.),
  fMaxPredDist(0.),
  fMaxPropDist(0.),
  fMinSortDist(0.),
  fMaxSortDist(0.),
  fMaxSamePlane(0),
  fGapDist(0.),
  fMaxNoiseHits(0),
  fMinSampleDist(0.),
  fFitMomRange(false),
  fFitMomMS(false),
  fPlane(-1)
{
  mf::LogInfo("KalmanFilterAlg") << "KalmanFilterAlg instantiated.";

  // Load fcl parameters.

  reconfigure(pset);
}

/// Destructor.
trkf::KalmanFilterAlg::~KalmanFilterAlg()
{}

/// Reconfigure method.
void trkf::KalmanFilterAlg::reconfigure(const fhicl::ParameterSet& pset)
{
  fTrace = pset.get<bool>("Trace");
  fMaxPErr = pset.get<double>("MaxPErr");
  fGoodPErr = pset.get<double>("GoodPErr");
  fMaxIncChisq = pset.get<double>("MaxIncChisq");
  fMaxEndChisq = pset.get<double>("MaxEndChisq");
  fMinLHits = pset.get<int>("MinLHits");
  fMaxLDist = pset.get<double>("MaxLDist");
  fMaxPredDist = pset.get<double>("MaxPredDist");
  fMaxPropDist = pset.get<double>("MaxPropDist");
  fMinSortDist = pset.get<double>("MinSortDist");
  fMaxSortDist = pset.get<double>("MaxSortDist");
  fMaxSamePlane = pset.get<int>("MaxSamePlane");
  fGapDist = pset.get<double>("GapDist");
  fMaxNoiseHits = pset.get<double>("MaxNoiseHits");
  fMinSampleDist = pset.get<double>("MinSampleDist");
  fFitMomRange = pset.get<bool>("FitMomRange");
  fFitMomMS = pset.get<bool>("FitMomMS");
}

/// Add hits to track.
///
/// Arguments:
///
/// trk      - Starting track.
/// trg      - Result global track.
/// prop     - Propagator.
/// dir      - Direction.
/// hits     - Candidate hits.
///
/// Returns: True if success.
///
/// This method makes a unidirectional Kalman fit in the specified
/// direction, visiting each surface of the passed KHitContainer and
/// updating the track.  In case of multiple measurements on the same
/// surface, keep (at most) the one with the smallest incremental
/// chisquare.  Any measurements that fail the incremental chisquare
/// cut are rejected.  Resolved hits (accepted or rejected) are moved
/// to the unused list in KHitContainer.  The container is sorted at
/// the start of the method, and may be resorted during the progress
/// of the fit.
///
bool trkf::KalmanFilterAlg::buildTrack(const KTrack& trk,
				       KGTrack& trg,
				       const Propagator* prop,
				       const Propagator::PropDirection dir,
				       KHitContainer& hits) const
{
  assert(prop != 0);

  // Direction must be forward or backward (unknown is not allowed).

  if(dir != Propagator::FORWARD && dir != Propagator::BACKWARD)
    throw cet::exception("KalmanFilterAlg") 
	<< "No direction for Kalman fit.\n";

  // Sort container using this seed track.

  hits.sort(trk, true, prop, Propagator::UNKNOWN);

  // Loop over measurements (KHitGroup) from sorted list.

  double tchisq = 0.;        // Cumulative chisquare.
  double path = 0.;          // Cumulative path distance.
  int step = 0;              // Step count.
  int nsame = 0;             // Number of consecutive measurements in same plane.
  int last_plane = -1;       // Last plane.

  // Make a copy of the starting track, in the form of a KFitTrack,
  // which we will update as we go.

  TrackError err;
  trk.getSurface()->getStartingError(err);
  KETrack tre(trk, err);
  KFitTrack trf(tre, path, tchisq);

  // Also use the starting track as the reference track for linearized
  // propagation, until the track is established with reasonably small
  // errors.

  KTrack ref(trk);
  KTrack* pref = &ref;

  mf::LogInfo log("KalmanFilterAlg");

  // Loop over measurement groups (KHitGroups).

  while(hits.getSorted().size() > 0) {
    ++step;
    if(fTrace) {
      log << "Build Step " << step << "\n";
      log << "KGTrack has " << trg.numHits() << " hits.\n";
      log << trf;
    }

    // Get an iterator for the next KHitGroup.

    std::list<KHitGroup>::iterator it;
    if(dir == Propagator::FORWARD)
      it = hits.getSorted().begin();
    else {
      assert(dir == Propagator::BACKWARD);
      it = hits.getSorted().end();
      --it;
    }
    const KHitGroup& gr = *it;

    if(fTrace) {
      double path_est = gr.getPath();
      log << "Next surface: " << *(gr.getSurface()) << "\n";
      log << "  Estimated distance = " << path_est << "\n";
    }

    // Get the next prediction surface.  If this KHitGroup is on the
    // preferred plane, use that as the prediction surface.
    // Otherwise, use the current track surface as the prediction
    // surface.

    std::shared_ptr<const Surface> psurf = trf.getSurface();
    assert(gr.getPlane() >= 0);
    if(fPlane < 0 || gr.getPlane() < 0 || fPlane == gr.getPlane())
      psurf = gr.getSurface();

    // Propagate track to the prediction surface.

    boost::optional<double> dist = prop->noise_prop(trf, psurf, Propagator::UNKNOWN, true, pref);
    if(!!dist && std::abs(*dist) > fMaxPropDist)
      dist = boost::optional<double>(false, 0.);
    double ds = 0.;

    if(!!dist) {

      // Propagation succeeded.
      // Update cumulative path distance and track status.

      ds = *dist;
      path += ds;
      trf.setPath(path);
      if(dir == Propagator::FORWARD)
	trf.setStat(KFitTrack::FORWARD_PREDICTED);
      else {
	assert(dir == Propagator::BACKWARD);
	trf.setStat(KFitTrack::BACKWARD_PREDICTED);
      }
      if(fTrace) {
	log << "After propagation\n";
	log << "  Incremental distance = " << ds << "\n";
	log << "  Actual distance = " << path << "\n";
	log << "KGTrack has " << trg.numHits() << " hits.\n";
	log << trf;
      }

      // Loop over measurements in this group.

      const std::vector<std::shared_ptr<const KHitBase> >& hits = gr.getHits();
      double best_chisq = 0.;
      std::shared_ptr<const KHitBase> best_hit;
      for(std::vector<std::shared_ptr<const KHitBase> >::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {
	const KHitBase& hit = **ihit;

	// Update predction using current track hypothesis and get
	// incremental chisquare.

	bool ok = hit.predict(trf, prop);
	if(ok) {
	  double chisq = hit.getChisq();
	  double preddist = hit.getPredDistance();
	  if((pref != 0 || best_chisq < fMaxIncChisq) && abs(preddist) < fMaxPredDist &&
	     (best_hit.get() == 0 || chisq < best_chisq) ) {
	    best_hit = *ihit;
	    best_chisq = chisq;
	  }
	}
      }
      if(fTrace) {
	if(best_hit.get() != 0) {
	  log << "Hit after prediction\n";
	  log << *best_hit;
	}
      }

      // If we found a best measurement, and if the incremental
      // chisquare passes the cut, add it to the track and update 
      // fit information.

      if(best_hit.get() != 0) {
	ds += best_hit->getPredDistance();
	best_hit->update(trf);
	tchisq += best_chisq;
	trf.setChisq(tchisq);
	if(dir == Propagator::FORWARD)
	  trf.setStat(KFitTrack::FORWARD);
	else {
	  assert(dir == Propagator::BACKWARD);
	  trf.setStat(KFitTrack::BACKWARD);
	}

	// If the pointing error got too large, and there is no
	// reference track, quit.

	if(pref == 0 && trf.PointingError() > fMaxPErr) {
	  if(fTrace)
	    log << "Quitting because pointing error got too large.\n";
	  break;
	}
	  
	// Test number of consecutive measurements in the same plane.

	if(gr.getPlane() >= 0) {
	  if(gr.getPlane() == last_plane)
	    ++nsame;
	  else {
	    nsame = 1;
	    last_plane = gr.getPlane();
	  }
	}
	else {
	  nsame = 0;
	  last_plane = -1;
	}
	if(nsame <= fMaxSamePlane) {

	  // Make a KHitTrack and add it to the KGTrack.

	  KHitTrack trh(trf, best_hit);
	  trg.addTrack(trh);

	  // Decide if we want to kill the reference track.

	  if(pref != 0 && int(trg.numHits()) >= fMinLHits &&
	     (trf.PointingError() < fGoodPErr || path > fMaxLDist)) {
	    pref = 0;
	    if(fTrace)
	      log << "Killing reference track.\n";
	  }

	  if(fTrace) {
	    log << "After update\n";
	    log << "KGTrack has " << trg.numHits() << " hits.\n";
	    log << trf;
	  }
	}
      }
    }

    // The current KHitGroup is now resolved.
    // Move it to unused list.

    hits.getUnused().splice(hits.getUnused().end(), hits.getSorted(), it);

    // If the propagation distance was the wrong direction, resort the measurements.

    if(pref == 0 && !!dist && 
       ((dir == Propagator::FORWARD && (ds < fMinSortDist || ds > fMaxSortDist)) ||
	(dir == Propagator::BACKWARD && (-ds < fMinSortDist || -ds > fMaxSortDist)))) {
      if(fTrace)
	log << "Resorting measurements.\n";
      hits.sort(trf, true, prop, dir);
    }
  }

  // Clean track.

  cleanTrack(trg);

  // Set the fit status of the last added KHitTrack to optimal and get
  // the final chisquare and path length.

  double fchisq = 0.;
  path = 0.;
  if(trg.isValid()) {
    path = trg.endTrack().getPath() - trg.startTrack().getPath();
    if(dir == Propagator::FORWARD) {
      trg.endTrack().setStat(KFitTrack::OPTIMAL);
      fchisq = trg.endTrack().getChisq();
    }
    else {
      assert(dir == Propagator::BACKWARD);
      trg.startTrack().setStat(KFitTrack::OPTIMAL);
      fchisq = trg.startTrack().getChisq();
    }
  }

  // Summary.

  log << "KalmanFilterAlg build track summary.\n"
      << "Build direction = " << (dir == Propagator::FORWARD ? "FORWARD" : "BACKWARD") << "\n"
      << "Track has " << trg.numHits() << " hits.\n"
      << "Track length = " << path << "\n"
      << "Track chisquare = " << fchisq << "\n";

  // Done.  Return success if we added at least one measurement.

  return trg.numHits() > 0;
}

/// Smooth track.
///
/// Arguments:
///
/// trg  - Track to be smoothed.
/// trg1 - Track to receive result of unidirectional fit.
/// prop - Propagator.
///
/// Returns: True if success.
///
/// The starting track should be a global track that has been fit in
/// one direction.  Fit status should be optimal at (at least) one
/// end.  It is an error if the fit status is not optimal at either
/// end.  If the fit status is optimal at both ends, do nothing, but
/// return success.
///
/// If the second argument is non-null, save the result of the
/// unidirectional track fit produced as a byproduct of the smoothing
/// operation.  This track can be smoothed in order to iterate the
/// Kalman fit, etc.
///
/// The Kalman smoothing algorithm starts at the optimal end and fits
/// the track in the reverse direction, calculating optimal track
/// parameters at each measurement surface.
///
/// All measurements are included in the reverse fit.  No incremental
/// chisquare cut is applied.
///
/// If any measurement surface can not be reached because of a
/// measurement error, the entire smoothing operation is considered a
/// failure.  In that case, false is returned and the track is left in
/// an undefined state.
///
bool trkf::KalmanFilterAlg::smoothTrack(KGTrack& trg,
					KGTrack* trg1,
					const Propagator* prop) const
{
  assert(prop != 0);

  // Default result failure.

  bool result = false;

  // It is an error if the KGTrack is not valid.

  if(trg.isValid()) {

    // Examine the track endpoints and figure out which end of the track
    // to start from.  The fit always starts at the optimal end.  It is
    // an error if neither end point is optimal.  Do nothing and return
    // success if both end points are optimal.

    const KHitTrack& trh0 = trg.startTrack();
    const KHitTrack& trh1 = trg.endTrack();
    KFitTrack::FitStatus stat0 = trh0.getStat();
    KFitTrack::FitStatus stat1 = trh1.getStat();
    bool dofit = false;

    // Remember starting direction, track, and distance.

    Propagator::PropDirection dir = Propagator::UNKNOWN;
    const KTrack* trk = 0;
    double path = 0.;

    if(stat0 == KFitTrack::OPTIMAL) {
      if(stat1 == KFitTrack::OPTIMAL) {

	// Both ends optimal (do nothing, return success).

	dofit = false;
	result = true;

      }
      else {

	// Start optimal.

	dofit = true;
	dir = Propagator::FORWARD;
	trk = &trh0;
	path = 0.;
      }
    }
    else {
      if(stat1 == KFitTrack::OPTIMAL) {

	// End optimal.

	dofit = true;
	dir = Propagator::BACKWARD;
	trk = &trh1;
	path = trh1.getPath();
      }
      else {

	// Neither end optimal (do nothing, return failure).

	dofit = false;
	result = false;
      }
    }
    if(dofit) {
      assert(dir == Propagator::FORWARD || dir == Propagator::BACKWARD);
      assert(trk != 0);

      // Cumulative chisquare.

      double tchisq = 0.;

      // Construct starting KFitTrack with track information and distance
      // taken from the optimal end, but with "infinite" errors.

      TrackError err;
      trk->getSurface()->getStartingError(err);
      KETrack tre(*trk, err);
      KFitTrack trf(tre, path, tchisq);

      // Make initial reference track to be same as initial fit track.

      KTrack ref(trf);

      // Loop over KHitTracks contained in KGTrack.

      std::multimap<double, KHitTrack>::iterator it;
      std::multimap<double, KHitTrack>::iterator itend;
      if(dir == Propagator::FORWARD) {
	it = trg.getTrackMap().begin();
	itend = trg.getTrackMap().end();
      }
      else {
	assert(dir == Propagator::BACKWARD);
	it = trg.getTrackMap().end();
	itend = trg.getTrackMap().begin();
      }

      mf::LogInfo log("KalmanFilterAlg");

      // Loop starts here.

      result = true;             // Result success unless we find an error.
      int step = 0;              // Step count.
      while(dofit && it != itend) {
	++step;
	if(fTrace) {
	  log << "Smooth Step " << step << "\n";
	  log << "Reverse fit track:\n";
	  log << trf;
	}

	// For backward fit, decrement iterator at start of loop.

	if(dir == Propagator::BACKWARD)
	  --it;

	KHitTrack& trh = (*it).second;
	if(fTrace) {
	  log << "Forward track:\n";
	  log << trh;
	}

	// Extract measurement.

	const KHitBase& hit = *(trh.getHit());

	// Propagate KFitTrack to the next track surface.
	  
	std::shared_ptr<const Surface> psurf = trh.getSurface();
	boost::optional<double> dist = prop->noise_prop(trf, psurf, Propagator::UNKNOWN,
							true, &ref);

	// Check if propagation succeeded.  If propagation fails, this
	// measurement will be dropped from the unidirectional fit
	// track.  This measurement will still be in the original
	// track, but with a status other than optimal.

	if(!!dist) {

	  // Propagation succeeded.
	  // Update cumulative path distance and track status.

	  double ds = *dist;
	  path += ds;
	  trf.setPath(path);
	  if(dir == Propagator::FORWARD)
	    trf.setStat(KFitTrack::FORWARD_PREDICTED);
	  else {
	    assert(dir == Propagator::BACKWARD);
	    trf.setStat(KFitTrack::BACKWARD_PREDICTED);
	  }
	  if(fTrace) {
	    log << "Reverse fit track after propagation:\n";
	    log << "  Propagation distance = " << ds << "\n";
	    log << trf;
	  }

	  // See if we have the proper information to calculate an optimal track
	  // at this surface (should normally be possible).

	  KFitTrack::FitStatus stat = trh.getStat();
	  KFitTrack::FitStatus newstat = trf.getStat();

	  if((newstat == KFitTrack::FORWARD_PREDICTED && stat == KFitTrack::BACKWARD) ||
	     (newstat == KFitTrack::BACKWARD_PREDICTED && stat == KFitTrack::FORWARD)) {

	    // Update stored KHitTrack to be optimal.

	    bool ok = trh.combineFit(trf);

	    // Update the stored path distance to be from the currently fitting track.

	    trh.setPath(trf.getPath());

	    // Update reference track.

	    ref = trh;

	    // If combination failed, abandon the fit and return failure.

	    if(!ok) {
	      dofit = false;
	      result = false;
	      break;
	    }
	    if(fTrace) {
	      log << "Combined track:\n";
	      log << trh;
	    }
	  }

	  // Update measurement predction using current track hypothesis.

	  bool ok = hit.predict(trf, prop, &ref);
	  if(!ok) {

	    // If prediction failed, abandon the fit and return failure.

	    dofit = false;
	    result = false;
	    break;
	  }
	  else {

	    // Prediction succeeded.  Get incremental chisquare.  If
	    // this hit fails incremental chisquare cut, this hit will
	    // be dropped from the unidirecitonal Kalman fit track,
	    // but may still be in the smoothed track.

	    double chisq = hit.getChisq();
	    if(chisq < fMaxIncChisq) {
	      tchisq += chisq;
	      trf.setChisq(tchisq);

	      // Update the reverse fitting track using the current measurement
	      // (both track parameters and status).

	      hit.update(trf);
	      if(dir == Propagator::FORWARD)
		trf.setStat(KFitTrack::FORWARD);
	      else {
		assert(dir == Propagator::BACKWARD);
		trf.setStat(KFitTrack::BACKWARD);
	      }
	      if(fTrace) {
		log << "Reverse fit track after update:\n";
		log << trf;
	      }

	      // If unidirectional track pointer is not null, make a
	      // KHitTrack and save it in the unidirectional track.

	      if(trg1 != 0) {
		KHitTrack trh1(trf, trh.getHit());
		trg1->addTrack(trh1);
	      }
	    }
	  }
	}

	// For forward fit, increment iterator at end of loop.

	if(dir == Propagator::FORWARD)
	  ++it;

      }    // Loop over KHitTracks.

      // If fit was successful and the unidirectional track pointer
      // is not null and the track is valid, set the fit status of
      // the last added KHitTrack to optimal.

      if(result && trg1 != 0 && trg1->isValid()) {
	if(dir == Propagator::FORWARD)
	  trg1->endTrack().setStat(KFitTrack::OPTIMAL);
	else {
	  assert(dir == Propagator::BACKWARD);
	  trg1->startTrack().setStat(KFitTrack::OPTIMAL);
	}
      }

      // Recalibrate track map.

      trg.recalibrate();

    }      // Do fit.

    // Get the final chisquare.

    double fchisq = 0.5 * (trg.startTrack().getChisq() + trg.endTrack().getChisq());

    // Summary.

    mf::LogInfo log("KalmanFilterAlg");
    log << "KalmanFilterAlg smooth track summary.\n"
	<< "Smooth direction = " << (dir == Propagator::FORWARD ? "FORWARD" : "BACKWARD") << "\n"
	<< "Track has " << trg.numHits() << " hits.\n"
	<< "Track length = " << trg.endTrack().getPath() - trg.startTrack().getPath() << "\n"
	<< "Track chisquare = " << fchisq << "\n";
  }

  // Done.

  return result;
}

/// Add hits to existing track.
///
/// Arguments:
///
/// trg - Track to extend.
/// prop - Propagator.
/// hits - Hit collection to choose hits from.
///
/// This method extends a KGTrack by adding hits.  The KGTrack must
/// have previously been produced by a unidirectional Kalman fit (it
/// should be optimal at one end).  This method finds the optimal end
/// and extends the track in that direction.  If any hits are added,
/// the originally optimal end has its status reset to forward or
/// backward, and the new endpoint is optimal.  In any case, the final
/// result is unidirectionally fit KGTrack.
///
bool trkf::KalmanFilterAlg::extendTrack(KGTrack& trg,
					const Propagator* prop,
					KHitContainer& hits) const
{
  assert(prop != 0);

  // Default result failure.

  bool result = false;

  // Remember the original number of measurement.

  unsigned int nhits0 = trg.numHits();

  // It is an error if the KGTrack is not valid.

  if(trg.isValid()) {
    mf::LogInfo log("KalmanFilterAlg");

    // Examine the track endpoints and figure out which end of the
    // track to extend.  The track is always extended from the optimal
    // end.  It is an error if neither end point is optimal, or both
    // endoints are optimal.  Reset the status of the optimal, and
    // make a copy of the starting track fit.  Also get starting path
    // and chisquare.

    KHitTrack& trh0 = trg.startTrack();
    KHitTrack& trh1 = trg.endTrack();
    KFitTrack::FitStatus stat0 = trh0.getStat();
    KFitTrack::FitStatus stat1 = trh1.getStat();
    bool dofit = false;
    Propagator::PropDirection dir = Propagator::UNKNOWN;
    KFitTrack trf;
    double path = 0.;
    double tchisq = 0.;

    if(stat0 == KFitTrack::OPTIMAL) {
      if(stat1 == KFitTrack::OPTIMAL) {

	// Both ends optimal (do nothing, return failure).

	dofit = false;
	result = false;
	return result;
      }
      else {

	// Start optimal.  Extend backward.

	dofit = true;
	dir = Propagator::BACKWARD;
	trh0.setStat(KFitTrack::BACKWARD);
	trf = trh0;
	path = trh0.getPath();
	tchisq = trh0.getChisq();
      }
    }
    else {
      if(stat1 == KFitTrack::OPTIMAL) {

	// End optimal.  Extend forward.

	dofit = true;
	dir = Propagator::FORWARD;
	trh1.setStat(KFitTrack::FORWARD);
	trf = trh1;
	path = trh1.getPath();
	tchisq = trh1.getChisq();

	// Make sure forward extend track momentum is over some
	// minimum value.

	if(trf.getVector()(4) > 2.) {
	  trf.getVector()(4) = 2.;
	  trf.getError()(4,4) = 2.;
	}
      }
      else {

	// Neither end optimal (do nothing, return failure).

	dofit = false;
	result = false;
	return result;
      }
    }
    if(dofit) {

      // Sort hit container using starting track.

      hits.sort(trf, true, prop, dir);

      // Extend loop starts here.

      int step = 0;
      int nsame = 0;
      int last_plane = -1;
      while(hits.getSorted().size() > 0) {
	++step;
	if(fTrace) {
	  log << "Extend Step " << step << "\n";
	  log << "KGTrack has " << trg.numHits() << " hits.\n";
	  log << trf;
	}

	// Get an iterator for the next KHitGroup.

	std::list<KHitGroup>::iterator it;
	if(dir == Propagator::FORWARD)
	  it = hits.getSorted().begin();
	else {
	  assert(dir == Propagator::BACKWARD);
	  it = hits.getSorted().end();
	  --it;
	}
	const KHitGroup& gr = *it;

	if(fTrace) {
	  double path_est = gr.getPath();
	  log << "Next surface: " << *(gr.getSurface()) << "\n";
	  log << "  Estimated distance = " << path_est << "\n";
	}

	// Get the next prediction surface.  If this KHitGroup is on the
	// preferred plane, use that as the prediction surface.
	// Otherwise, use the current track surface as the prediction
	// surface.

	std::shared_ptr<const Surface> psurf = trf.getSurface();
	assert(gr.getPlane() >= 0);
	if(fPlane < 0 || gr.getPlane() < 0 || fPlane == gr.getPlane())
	  psurf = gr.getSurface();

	// Propagate track to the prediction surface.

	boost::optional<double> dist = prop->noise_prop(trf, psurf, Propagator::UNKNOWN, true);
	if(!!dist && std::abs(*dist) > fMaxPropDist)
	  dist = boost::optional<double>(false, 0.);
	double ds = 0.;

	if(!!dist) {

	  // Propagation succeeded.
	  // Update cumulative path distance and track status.

	  ds = *dist;
	  path += ds;
	  trf.setPath(path);
	  if(dir == Propagator::FORWARD)
	    trf.setStat(KFitTrack::FORWARD_PREDICTED);
	  else {
	    assert(dir == Propagator::BACKWARD);
	    trf.setStat(KFitTrack::BACKWARD_PREDICTED);
	  }
	  if(fTrace) {
	    log << "After propagation\n";
	    log << "  Incremental distance = " << ds << "\n";
	    log << "  Actual distance = " << path << "\n";
	    log << "KGTrack has " << trg.numHits() << " hits.\n";
	    log << trf;
	  }

	  // Loop over measurements in this group.

	  const std::vector<std::shared_ptr<const KHitBase> >& hits = gr.getHits();
	  double best_chisq = 0.;
	  std::shared_ptr<const KHitBase> best_hit;
	  for(std::vector<std::shared_ptr<const KHitBase> >::const_iterator ihit = hits.begin();
	      ihit != hits.end(); ++ihit) {
	    const KHitBase& hit = **ihit;

	    // Update predction using current track hypothesis and get
	    // incremental chisquare.

	    bool ok = hit.predict(trf, prop);
	    if(ok) {
	      double chisq = hit.getChisq();
	      double preddist = hit.getPredDistance();
	      if(chisq < fMaxIncChisq && abs(preddist) < fMaxPredDist &&
		 (best_hit.get() == 0 || chisq < best_chisq)) {
		best_hit = *ihit;
		best_chisq = chisq;
	      }
	    }
	  }
	  if(fTrace) {
	    if(best_hit.get() != 0) {
	      log << "Hit after prediction\n";
	      log << *best_hit;
	    }
	  }

	  // If we found a best measurement, and if the incremental
	  // chisquare passes the cut, add it to the track and update 
	  // fit information.

	  if(best_hit.get() != 0) {
	    ds += best_hit->getPredDistance();
	    best_hit->update(trf);
	    tchisq += best_chisq;
	    trf.setChisq(tchisq);
	    if(dir == Propagator::FORWARD)
	      trf.setStat(KFitTrack::FORWARD);
	    else {
	      assert(dir == Propagator::BACKWARD);
	      trf.setStat(KFitTrack::BACKWARD);
	    }

	    // If the pointing error got too large, quit.

	    if(trf.PointingError() > fMaxPErr) {
	      if(fTrace)
		log << "Quitting because pointing error got too large.\n";
	      break;
	    }

	    // Test number of consecutive measurements in the same plane.

	    if(gr.getPlane() >= 0) {
	      if(gr.getPlane() == last_plane)
		++nsame;
	      else {
		nsame = 1;
		last_plane = gr.getPlane();
	      }
	    }
	    else {
	      nsame = 0;
	      last_plane = -1;
	    }
	    if(nsame <= fMaxSamePlane) {

	      // Make a KHitTrack and add it to the KGTrack.

	      KHitTrack trh(trf, best_hit);
	      trg.addTrack(trh);

	      if(fTrace) {
		log << "After update\n";
		log << "KGTrack has " << trg.numHits() << " hits.\n";
		log << trf;
	      }
	    }
	  }
	}

	// The current KHitGroup is now resolved.
	// Move it to unused list.

	hits.getUnused().splice(hits.getUnused().end(), hits.getSorted(), it);

	// If the propagation distance was the wrong direction, resort the measurements.

	if(!!dist &&
	   ((dir == Propagator::FORWARD && (ds < fMinSortDist || ds > fMaxSortDist)) ||
	    (dir == Propagator::BACKWARD && (-ds < fMinSortDist || -ds > fMaxSortDist)))) {
	  if(fTrace)
	    log << "Resorting measurements.\n";
	  hits.sort(trf, true, prop, dir);
	}
      }
    }

    // Clean track.

    cleanTrack(trg);

    // Set the fit status of the last added KHitTrack to optimal and
    // get the final chisquare and path length.

    double fchisq = 0.;
    path = 0.;
    if(trg.isValid()) {
      path = trg.endTrack().getPath() - trg.startTrack().getPath();
      if(dir == Propagator::FORWARD) {
	trg.endTrack().setStat(KFitTrack::OPTIMAL);
	fchisq = trg.endTrack().getChisq();
      }
      else {
	assert(dir == Propagator::BACKWARD);
	trg.startTrack().setStat(KFitTrack::OPTIMAL);
	fchisq = trg.startTrack().getChisq();
      }
    }

    // Summary.

    log << "KalmanFilterAlg extend track summary.\n"
	<< "Extend direction = " << (dir == Propagator::FORWARD ? "FORWARD" : "BACKWARD") << "\n"
	<< "Track has " << trg.numHits() << " hits.\n"
	<< "Track length = " << path << "\n"
	<< "Track chisquare = " << fchisq << "\n";
  }

  // Done.

  result = (trg.numHits() > nhits0);
  return result;
}

/// Estimate track momentum using range.
///
/// Arguments:
///
/// trg    - Global track.
/// prop   - Propagator.
/// tremom - Track with momentum estimate.
///
/// Returns: True if success.
///
/// This method generates a momentum-estimating track by extracting
/// the last track from a global track, and setting its momentum to
/// some small value.
///
bool trkf::KalmanFilterAlg::fitMomentumRange(const KGTrack& trg,
					     const Propagator* prop,
					     KETrack& tremom) const
{
  if(!trg.isValid())
    return false;

  // Extract track with lowest momentum.

  const KHitTrack& trh = trg.endTrack();
  assert(trh.getStat() != KFitTrack::INVALID);
  tremom = trh;

  // Set track momentum to a small value.

  tremom.getVector()(4) = 100.;
  tremom.getError()(4,0) = 0.;
  tremom.getError()(4,1) = 0.;
  tremom.getError()(4,2) = 0.;
  tremom.getError()(4,3) = 0.;
  tremom.getError()(4,4) = 10000.;

  // Done.

  return true;
}

/// Estimate track momentum using multiple scattering.
///
/// Arguments:
///
/// trg    - Global track.
/// prop   - Propagator.
/// tremom - Track containing momentum estimate.
///
/// Returns: True if success.
///
/// This method estimates the momentum of the specified track using
/// multiple scattering.  As a result of calling this method, the
/// original global track is not updated, but a KETrack is produced
/// near the starting surface that has the estimated momentum.
///
/// The global track passed as argument should have been smoothed
/// prior to calling this method so that all or most fits are optimal.
/// If either the first or last fit is not optimal, return false.
///
/// This method assumes that track parameter four is 1/p.  This sort
/// of momentum estimation only makes sense if the momentum track
/// parameter is uncorrelated with any other track parameter.  The
/// error matrix of the first and last fit is checked for this.  If it
/// is found that either error matrix has correlated track parameter
/// four with any other track parameter, this method returns without
/// doing anything (return false).
///
bool trkf::KalmanFilterAlg::fitMomentumMS(const KGTrack& trg,
					  const Propagator* prop,
					  KETrack& tremom) const
{
  // Get iterators pointing to the first and last tracks.

  const std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();
  if(trackmap.size() < 2)
    return false;
  std::multimap<double, KHitTrack>::const_iterator itend[2];
  itend[0] = trackmap.begin();
  itend[1] = trackmap.end();
  --itend[1];

  // Check the fit status and error matrix of the first and last
  // track.

  bool result = true;

  for(int i=0; result && i<2; ++i) {
    const KHitTrack& trh = itend[i]->second;
    KFitTrack::FitStatus stat = trh.getStat();
    if(stat != KFitTrack::OPTIMAL)
      result = false;
    const TrackError& err = trh.getError();
    for(int j=0; j<4; ++j) {
      if(err(4,j) != 0.)
	result = false;
    }
  }
  if(!result)
    return result;

  // We will periodically sample the track trajectory.  At each sample
  // point, collect the following information.
  //
  // 1.  The path distance at the sample point.
  // 2.  The original momentum of the track at the sample point and its error.
  // 3.  One copy of the track (KETrack) that will be propagated without noise
  //     (infinite momentum track).
  // 3.  A second copy of the track (KETrack) that will be propagated with the
  //     minimum allowed momentum (range out track).
  // 4.  A third copy of the track (KETrack) that will propagated with noise
  //     with some intermediate momentum (noise track).
  //
  // Collect the first sample from the maximum path distance track.

  double s_sample = itend[1]->first;
  const KETrack& tre = itend[1]->second;
  KETrack tre_inf(tre);
  KTrack trk_range(tre);
  KETrack tre_noise(tre);
  tre_inf.getVector()(4) = 0.;
  tre_inf.getError()(4,4) = 0.;
  trk_range.getVector()(4) = 100.;
  tre_noise.getError()(4,4) = 0.;
  tre_noise.getVector()(4) = 1.;
  tre_noise.getError()(4,4) = 10.;
  double invp0 = tre_noise.getVector()(4);
  double var_invp0 = tre_noise.getError()(4,4);

  // Loop over fits, starting at the high path distance (low momentum)
  // end.  

  for(std::multimap<double, KHitTrack>::const_reverse_iterator it = trackmap.rbegin();
      it != trackmap.rend(); ++it) {
    double s = it->first;
    const KHitTrack& trh = it->second;

    // Ignore non-optimal fits.

    KFitTrack::FitStatus stat = trh.getStat();
    if(stat != KFitTrack::OPTIMAL)
      continue;

    // See if this track is far enough from the previous sample to
    // make a new sample point.

    if(std::abs(s - s_sample) > fMinSampleDist) {

      // Propagate tracks to the current track surface.

      std::shared_ptr<const Surface> psurf = trh.getSurface();
      boost::optional<double> dist_inf = prop->err_prop(tre_inf, psurf,
							Propagator::UNKNOWN, false);
      boost::optional<double> dist_range = prop->vec_prop(trk_range, psurf,
							  Propagator::UNKNOWN, false);
      boost::optional<double> dist_noise = prop->noise_prop(tre_noise, psurf,
							    Propagator::UNKNOWN, true);

      // All propagations should normally succeed.  If they don't,
      // ignore this sample for the purpose of updating the momentum.

      bool momentum_updated = false;
      if(!!dist_inf && !!dist_range && !!dist_noise) {

	// Extract the momentum at the new sample point.

	double invp1 = tre_noise.getVector()(4);
	double var_invp1 = tre_noise.getError()(4,4);

	// Get the average momentum and error for this pair of
	// sample points, and other data.

	double invp = 0.5 * (invp0 + invp1);
	double var_invp = 0.5 * (var_invp0 + var_invp1);
	double mass = tre_inf.Mass();
	double beta = std::sqrt(1. + mass*mass*invp*invp);
	double invbp = invp / beta;

	// Extract slope subvectors and sub-error-matrices.
	// We have the following variables.
	//
	// slope0 - Predicted slope vector (from infinite momentum track).
	// slope1 - Measured slope vector (from new sample point).
	// defl - Deflection (slope residual = difference between measured 
	//        and predicted slope vector).
	// err0 - Slope error matrix of prediction.
	// err1 - Slope error matrix of measurement
	// errc - Slope residual error matrix = err0 + err1.
	// errn - Noise slope error matrix divided by (1/beta*momentum).

	KVector<2>::type slope0 = project(tre_inf.getVector(), ublas::range(2, 4));
	KVector<2>::type slope1 = project(trh.getVector(), ublas::range(2, 4));
	KVector<2>::type defl = slope1 - slope0;
	KSymMatrix<2>::type err0 =
	  project(tre_inf.getError(), ublas::range(2, 4), ublas::range(2, 4));
	KSymMatrix<2>::type err1 =
	  project(trh.getError(), ublas::range(2, 4), ublas::range(2, 4));
	KSymMatrix<2>::type errc = err0 + err1;
	KSymMatrix<2>::type errn =
	  project(tre_noise.getError(), ublas::range(2, 4), ublas::range(2, 4));
	errn -= err0;
	errn /= invbp;

	// Calculate updated average momentum and error.

	double new_invp = invp;
	double new_var_invp = var_invp;
	update_momentum(defl, errc, errn, mass, new_invp, new_var_invp);

	// Calculate updated momentum and error at the second sample
	// point.

	double dp = 1./new_invp - 1./invp;
	invp0 = 1./(1./invp1 + dp);
	var_invp0 = new_var_invp;
	momentum_updated = true;

	// Make sure that updated momentum is not less than minimum
	// allowed momentum.

	double invp_range = trk_range.getVector()(4);
	if(invp0 > invp_range)
	  invp0 = invp_range;
      }

      // Update sample.

      if(momentum_updated) {
	s_sample = s;
	tre_inf = trh;
	tre_inf.getVector()(4) = 0.;
	tre_inf.getError()(4,4) = 0.;
	double invp_range = trk_range.getVector()(4);
	trk_range = trh;
	trk_range.getVector()(4) = invp_range;
	tre_noise = trh;
	tre_noise.getVector()(4) = invp0;
	tre_noise.getError()(4,4) = var_invp0;
      }
    }
  }

  // Propagate noise track to starting (high momentum) track surface
  // to get final starting momentum.  This propagation should normally
  // always succeed, but if it doesn't, don't update the track.

  const KHitTrack& trh0 = itend[0]->second;
  std::shared_ptr<const Surface> psurf = trh0.getSurface();
  boost::optional<double> dist_noise = prop->noise_prop(tre_noise, psurf,
							Propagator::UNKNOWN, true);
  result = !!dist_noise;

  // Update momentum-estimating track.

  mf::LogInfo log("KalmanFilterAlg");
  if(result)
    tremom = tre_noise;

  // Done.

  return result;
}

/// Estimate track momentum using either range or multiple scattering.
///
/// Arguments:
///
/// trg    - Global track whose momentum is to be updated.
/// prop   - Propagator.
/// tremom - Track with momentum estimate.
///
/// Returns: True if success.
///
bool trkf::KalmanFilterAlg::fitMomentum(const KGTrack& trg,
					const Propagator* prop,
					KETrack& tremom) const
{
  mf::LogInfo log("KalmanFilterAlg");
  double invp_range = 0.;
  double invp_ms = 0.;

  // Get multiple scattering momentum estimate.

  KETrack tremom_ms;
  bool ok_ms = false;
  if(fFitMomMS) {
    ok_ms = fitMomentumMS(trg, prop, tremom_ms);
    if(ok_ms) {
      KGTrack trg_ms(trg);
      ok_ms = updateMomentum(tremom_ms, prop, trg_ms);
      if(ok_ms) {
	invp_ms = trg_ms.startTrack().getVector()(4);
	double var_invp = trg_ms.startTrack().getError()(4,4);
	double p = 0.;
	if(invp_ms != 0.)
	  p = 1./invp_ms;
	double err_p = p*p * std::sqrt(var_invp);
	log << "Multiple scattering momentum estimate = " << p << "+-" << err_p << "\n";
      }
    }
  }

  // Get range momentum estimate.

  KETrack tremom_range;
  bool ok_range = false;
  if(fFitMomRange) {
    ok_range = fitMomentumRange(trg, prop, tremom_range);
    if(ok_range) {
      KGTrack trg_range(trg);
      ok_range = updateMomentum(tremom_range, prop, trg_range);
      if(ok_range) {
	invp_range = trg_range.startTrack().getVector()(4);
	double var_invp = trg_range.startTrack().getError()(4,4);
	double p = 0.;
	if(invp_range != 0.)
	  p = 1./invp_range;
	double err_p = p*p * std::sqrt(var_invp);
	log << "Range momentum estimate               = " << p << "+-" << err_p << "\n";
      }
    }
  }

  bool result = false;
  if(ok_range) {
    tremom = tremom_range;
    result = ok_range;
  }
  else if(ok_ms) {
    tremom = tremom_ms;
    result = ok_ms;
  }
  return result;
}

/// Set track momentum at each track surface.
///
/// Arguments:
///
/// tremom - Track containing momentum estimate.
/// prop   - Propagator.
/// trg    - Global track to be updated.
///
/// The track containing the momentum estimate is propagated to the
/// first or last track fit of the global track (whichever is closer),
/// then the momentum estimate is transfered to that track fit.  In
/// similar fashion, the momentum estimate is successively transfered
/// from that track fit to each track fit of the global track.
///
/// Only momentum track parameters of the global track fits are
/// updated.  Other track parameters and their errors are unmodified.
/// Unreachable track fits are deleted from the global track.  Overall
/// failure will occur if the momentum-estimating track can't be
/// propagated to the initial track fit, or if the final global track
/// has no valid track fits.
///
bool trkf::KalmanFilterAlg::updateMomentum(const KETrack& tremom,
					   const Propagator* prop,
					   KGTrack& trg) const
{
  // Get modifiable track map.

  std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();

  // If track map is empty, immediately return failure.

  if(trackmap.size() == 0)
    return false;

  // Make trial propagations to the first and last track fit to figure
  // out which track fit is closer to the momentum estimating track.

  // Find distance to first track fit.
  
  KETrack tre0(tremom);
  boost::optional<double> dist0 = prop->vec_prop(tre0, trackmap.begin()->second.getSurface(),
						 Propagator::UNKNOWN, false, 0, 0);
  // Find distance to last track fit.
  
  KETrack tre1(tremom);
  boost::optional<double> dist1 = prop->vec_prop(tre1, trackmap.rbegin()->second.getSurface(),
						 Propagator::UNKNOWN, false, 0, 0);

  // Based on distances, make starting iterator and direction flag.

  bool forward = true;
  std::multimap<double, KHitTrack>::iterator it = trackmap.begin();
  if(!!dist0) {

    // Propagation to first track succeeded.

    if(!!dist1) {

      // Propagation to both ends succeeded.  If the last track is
      // closer, initialize iterator and direction flag for reverse
      // direction.

      if(std::abs(*dist0) > std::abs(*dist1)) {
	it = trackmap.end();
	--it;
	forward = false;
      }
    }
  }
  else {

    // Propagation to first track failed.  Initialize iterator and
    // direction flag for reverse direction, provided that the
    // propagation to the last track succeeded.

    if(!!dist1) {
      it = trackmap.end();
      --it;
      forward = false;
    }
    else {

      // Propagation to both ends failed.  Return failure.

      return false;
    }
  }

  // Loop over track fits in global track.
  
  KETrack tre(tremom);
  bool done = false;
  while(!done) {
    KHitTrack& trh = it->second;
    assert(trh.getStat() != KFitTrack::INVALID);

    // Propagate momentum-estimating track to current track surface
    // and update momentum.

    boost::optional<double> dist = prop->noise_prop(tre, trh.getSurface(),
						    Propagator::UNKNOWN, true);

    // Copy momentum to global track.

    std::multimap<double, KHitTrack>::iterator erase_it = trackmap.end();
    if(!!dist) {
      trh.getVector()(4) = tre.getVector()(4);
      trh.getError()(4,0) = 0.;
      trh.getError()(4,1) = 0.;
      trh.getError()(4,2) = 0.;
      trh.getError()(4,3) = 0.;
      trh.getError()(4,4) = tre.getError()(4,4);
    }
    else {

      // If the propagation failed, remember that we are supposed to
      // erase this track from the global track.

      erase_it = it;
    }

    // Advance the iterator and set the done flag.

    if(forward) {
      ++it;
      done = (it == trackmap.end());
    }
    else {
      done = (it == trackmap.begin());
      if(!done)
	--it;
    }

    // Update momentum-estimating track from just-updated global track
    // fit, or erase global track fit.

    if(erase_it == trackmap.end())
      tre = trh;
    else
      trackmap.erase(erase_it);
  }

  bool result = (trackmap.size() > 0);

  // Print value of momentum at start of track.

  //if(result) {
  //  mf::LogInfo log("KalmanFilterAlg");
  //  double invp = trg.startTrack().getVector()(4);
  //  double var_invp = trg.startTrack().getError()(4,4);
  //  double p = 0.;
  //  if(invp != 0.)
  //    p = 1./invp;
  //  double err_p = p*p * std::sqrt(var_invp);
  //  log << "Updated momentum estimate = " << p << "+-" << err_p;
  //}

  return result;
}

/// Iteratively smooth a track.
///
/// Arguments:
///
/// nsmooth - Number of iterations.
/// trg     - Track to be smoothed.
/// prop    - Propagator.
///
/// Returns: True if success.
///
/// The initial track should have been unidirectionally fit.
///
bool trkf::KalmanFilterAlg::smoothTrackIter(int nsmooth,
					    KGTrack& trg,
					    const Propagator* prop) const
{
  bool ok = true;

  // Call smoothTrack method in a loop.

  for(int ismooth = 0; ok && ismooth < nsmooth-1; ++ismooth) {
    KGTrack trg1;
    ok = smoothTrack(trg, &trg1, prop);
    if(ok)
      trg = trg1;
  }

  // Do one final smooth without generating a new unidirectional
  // track.

  if(ok)
    ok = smoothTrack(trg, 0, prop);

  // Done.

  return ok;
}

/// Clean track by removing noise hits near endpoints.
///
/// Arguments:
///
/// trg - Track.
///
void trkf::KalmanFilterAlg::cleanTrack(KGTrack& trg) const
{
  // Get hold of a modifiable track map from trg.

  std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();

  // Do an indefinite loop until we no longer find any dirt.

  bool done = false;
  while(!done) {

    // If track map has fewer than fMaxNoiseHits tracks, then this is a
    // noise track.  Clear the map, making the whole track invalid.

    int ntrack = trackmap.size();
    if(ntrack <= fMaxNoiseHits) {
      trackmap.clear();
      done = true;
      break;
    }

    // Make sure the first two and last two tracks belong to different
    // views.  If not, remove the first or last track.

    if(ntrack >= 2) {

      // Check start.

      std::multimap<double, KHitTrack>::iterator it = trackmap.begin();
      const KHitTrack& trh1a = (*it).second;
      const KHitBase& hit1a = *(trh1a.getHit());
      int plane1 = hit1a.getMeasPlane();
      double chisq1 = hit1a.getChisq();
      ++it;
      const KHitTrack& trh2a = (*it).second;
      const KHitBase& hit2a = *(trh2a.getHit());
      int plane2 = hit2a.getMeasPlane();
      double chisq2 = hit2a.getChisq();
      if((plane1 >= 0 && plane2 >= 0 && plane1 == plane2) ||
	 chisq1 > fMaxEndChisq || chisq2 > fMaxEndChisq) {
	trackmap.erase(trackmap.begin(), it);
	done = false;
	continue;
      }

      // Check end.

      it = trackmap.end();
      --it;
      const KHitTrack& trh1b = (*it).second;
      const KHitBase& hit1b = *(trh1b.getHit());
      plane1 = hit1b.getMeasPlane();
      chisq1 = hit1b.getChisq();
      --it;
      const KHitTrack& trh2b = (*it).second;
      const KHitBase& hit2b = *(trh2b.getHit());
      plane2 = hit2b.getMeasPlane();
      chisq2 = hit2b.getChisq();
      if((plane1 >= 0 && plane2 >= 0 && plane1 == plane2) ||
	 chisq1 > fMaxEndChisq || chisq2 > fMaxEndChisq) {
	++it;
	trackmap.erase(it, trackmap.end());
	done = false;
	continue;
      }
    }

    // Loop over successive pairs of elements of track map.  Look for
    // adjacent elements with distance separation greater than
    // fGapDist.

    std::multimap<double, KHitTrack>::iterator it = trackmap.begin();
    std::multimap<double, KHitTrack>::iterator jt = trackmap.end();
    int nb = 0;                // Number of elements from begin to jt.
    int ne = ntrack;           // Number of elements it to end.
    bool found_noise = false;
    for(; it != trackmap.end(); ++it, ++nb, --ne) {
      if(jt == trackmap.end())
	jt = trackmap.begin();
      else {
	assert(nb >= 1);
	assert(ne >= 1);
	double disti = (*it).first;
	double distj = (*jt).first;
	double sep = disti - distj;
	assert(sep >= 0.);
	if(sep > fGapDist) {

	  // Found a gap.  See if we want to trim track.

	  if(nb <= fMaxNoiseHits) {

	    // Trim front.

	    found_noise = true;
	    trackmap.erase(trackmap.begin(), it);
	    break;
	  }
	  if(ne <= fMaxNoiseHits) {

	    // Trim back.

	    found_noise = true;
	    trackmap.erase(it, trackmap.end());
	    break;
	  }
	}
	++jt;
      }
    }

    // Decide if we are done.

    done = !found_noise;
  }
}
