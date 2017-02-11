#include "TrackKalmanFitter.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larreco/TrackFinder/TrackingPlaneHelper.h"
#include "lardata/RecoObjects/KFTrackState.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

/*

  TO DO LIST

  - fitTrack should work also with Trajectory and TrackTrajectory as input
  + check on input if there are bad flags
  + fitTrack should not reject hits but should rather set the flags
  - check what is the meaning of hittpindex
  + Update propagator so that copies are minimized (propagate should return the object)
  + Create a method for propagation in place
  + Make sure bools are returned for failures (propagation, state update, state combine)
  - Find a way to make propagator smarter in terms of the number of steps to do in one propagation (possibly remove masStep)
  + Update name of propagator class
  + Revisit TrackFitMeasurement to minimize copies and avoid duplicate information (plane? make sure hit and track state are on the same plane)
  + Revisit KFTrackState so that it makes sense both for propagated and updates states (add bool?)
  + Revisit fitTrack to minimize copies
  + Retrieve services once for all, and store as members of the fitter class
  + Create an object to be saved for hit-on-track information
  + TrackKalmanFitter should produce hit-on-track objects
  - check various fixme in the code for more todo...
  - add print/dump functions
  - document stuff!
  - remove tracking types that are not used
 */

bool trkf::TrackKalmanFitter::fitTrack(const recob::Track& track, const std::vector<art::Ptr<recob::Hit> >& hits,
				       const double pval, const int pdgid, const bool flipDirection,
				       recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits,
				       std::vector<recob::TrackFitHitInfo>& trackFitHitInfos) {

  auto position = track.Vertex();
  auto direction = track.VertexDirection();

  if (flipDirection) {
    position = track.End();
    direction = -track.EndDirection();
  }

  SVector5 trackStatePar(0.,0.,0.,0.,1./pval);
  auto trackStateCov = (flipDirection ? track.EndCovarianceLocal5D() : track.VertexCovarianceLocal5D() );
  if (track.NumberCovariance()==0) {
    trackStateCov(0, 0) = 1000.;
    trackStateCov(1, 1) = 1000.;
    trackStateCov(2, 2) = 0.25;
    trackStateCov(3, 3) = 0.25;
    trackStateCov(4, 4) = 10.;
  }

  Point_t  op(position(0),position(1),position(2));
  Vector_t od(direction(0),direction(1),direction(2));

  // setup the KFTrackState we'll use throughout the fit
  KFTrackState trackState(trackStatePar, trackStateCov, Plane(op,od), true, pdgid);//along direction by definition

  // figure out hit sorting based on minimum distance to first or last hit
  // (to be removed once there is a clear convention for hit sorting)
  geo::WireGeo const& wgeomF = geom->WireIDToWireGeo(hits.front()->WireID());
  geo::WireGeo const& wgeomB = geom->WireIDToWireGeo(hits.back()->WireID());
  Plane pF = recob::tracking::makePlane(wgeomF);
  Plane pB = recob::tracking::makePlane(wgeomB);
  bool success = true;
  double distF = propToPlane->distanceToPlane(success, op, od, pF);
  double distB = propToPlane->distanceToPlane(success, op, od, pB);
  if (!success) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }
  const bool reverseHits = distB<distF;

  // setup vector of HitStates and flags, with either same or inverse order as input hit vector
  // this is what we'll loop over during the fit
  std::vector<HitState>                            hitstatev;
  std::vector<recob::TrajectoryPointFlags::Mask_t> hitflagsv;
  unsigned int nplanes = 0;
  const int beg = (reverseHits ? hits.size()-1 : 0);
  const int end = (reverseHits ? -1 : hits.size());
  for (int ihit = beg; ihit!=end; (reverseHits ? ihit-- : ihit++)) {
    const auto& hit = hits[ihit];
    double t = hit->PeakTime();
    double terr = (useRMS_ ? hit->RMS() : hit->SigmaPeakTime() );
    double x = detprop->ConvertTicksToX(t, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
    double xerr = terr * detprop->GetXTicksCoefficient();
    hitstatev.push_back( std::move( HitState(x,hitErrScaleFact_*xerr*xerr,hit->WireID(),geom->WireIDToWireGeo(hit->WireID())) ) );
    hitflagsv.push_back( track.FlagsAtPoint(ihit).mask() );
    if ((hit->WireID().Plane+1)>nplanes) nplanes = hit->WireID().Plane+1;
  }

  // setup the track vectors we use to store the fit results
  // these three vectors are aligned
  std::vector<KFTrackState> fwdPrdTkState;
  std::vector<KFTrackState> fwdUpdTkState;
  std::vector<unsigned int> hitstateidx;
  std::vector<unsigned int> rejectedhsidx;
  fwdPrdTkState.reserve(hits.size());
  fwdUpdTkState.reserve(hits.size());
  hitstateidx.reserve(hits.size());

  if (sortHitsByPlane_) {
    //array of hit indices in planes, keeping the original sorting by plane
    std::vector<unsigned int> hitsInPlanes[nplanes] = { };
    for (unsigned int ihit = 0; ihit<hitstatev.size(); ihit++) {
      hitsInPlanes[hitstatev[ihit].wireId().Plane].push_back(ihit);
    }
    //array of indices, where iterHitsInPlanes[i] is the iterator over hitsInPlanes[i]
    unsigned int iterHitsInPlanes[nplanes] = {0};
    for (unsigned int p = 0; p<hitstatev.size(); ++p) {
      int min_plane = -1;
      double min_dist = DBL_MAX;
      //find the closest hit according to the sorting in each plane
      for (unsigned int iplane = 0; iplane<nplanes; ++iplane) {
	//note: ih is a reference, so when 'continue' we increment iterHitsInPlanes[iplane] and the hit is effectively rejected
	for (unsigned int& ih = iterHitsInPlanes[iplane]; ih<hitsInPlanes[iplane].size(); ++ih) {
	  unsigned int ihit = hitsInPlanes[iplane][ih];
	  const auto& hitstate = hitstatev[ihit];
	  const auto& hitflags = hitflagsv[ihit];
	  if (hitflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint        ) ||
	      hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ||
	      hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected       )) {
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  //get distance to measurement surface
	  bool success = true;
	  const double dist = propToPlane->distanceToPlane(success, trackState.trackState(), hitstate.plane());
	  if (!success) {
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  if (skipNegProp_ && dist<0.) {
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  if (dist<min_dist) {
	    min_plane = iplane;
	    min_dist = dist;
	  }
	  break;
	}
      }
      //now we know which is the closest wire: get the hitstate and increment the iterator
      if (min_plane<0) continue;
      const unsigned int ihit = hitsInPlanes[min_plane][iterHitsInPlanes[min_plane]];
      const auto& hitstate = hitstatev[ihit];
      auto& hitflags = hitflagsv[ihit];
      iterHitsInPlanes[min_plane]++;
      //propagate to measurement surface
      bool propok = true;
      trackState = propToPlane->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true);
      if (propok) {
	hitstateidx.push_back(ihit);
	fwdPrdTkState.push_back(trackState);
	//
	if (hitflags.isSet(recob::TrajectoryPointFlagTraits::HitIgnored   ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::Merged       ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::Shared       ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::DeltaRay     ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::DetectorIssue) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::Suspicious   )) {
	  //do not update the hit, mark as excluded from the fit
	  fwdUpdTkState.push_back(trackState);
	  hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
	  continue;
	}
	//
	//now update the forward fitted track
	bool upok = trackState.updateWithHitState(hitstate);
	if (upok==0) {
	  mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	  return false;
	}
	fwdUpdTkState.push_back(trackState);
      } else {
	//mf::LogWarning("TrackKalmanFitter") << "WARNING: forward propagation failed. Skip this hit...";
	rejectedhsidx.push_back(ihit);
	continue;
      }
    }
  } else {
    for (unsigned int ihit=0; ihit<hitstatev.size(); ++ihit) {
      const auto& hitstate = hitstatev[ihit];
      auto& hitflags = hitflagsv[ihit];
      if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ||
	  hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected)) {
	rejectedhsidx.push_back(ihit);
	continue;
      }
      if (skipNegProp_) {
	bool success = true;
	const double dist = propToPlane->distanceToPlane(success, trackState.trackState(), hitstate.plane());
	if (dist<0. || success==false) {
	  //mf::LogWarning("TrackKalmanFitter") << "WARNING: negative propagation distance. Skip this hit...";
	  rejectedhsidx.push_back(ihit);
	  continue;
	}
      }
      //propagate to measurement surface
      bool propok = true;
      trackState = propToPlane->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true);
      if (propok) {
	hitstateidx.push_back(ihit);
	fwdPrdTkState.push_back(trackState);
	//
	if (hitflags.isSet(recob::TrajectoryPointFlagTraits::HitIgnored   ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::Merged       ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::Shared       ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::DeltaRay     ) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::DetectorIssue) ||
	    hitflags.isSet(recob::TrajectoryPointFlagTraits::Suspicious   )) {
	  //do not update the hit, mark as excluded from the fit
	  fwdUpdTkState.push_back(trackState);
	  hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
	  continue;
	}
	//
	//now update the forward fitted track
	bool upok = trackState.updateWithHitState(hitstate);
	if (upok==0) {
	  mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	  return false;
	}
	fwdUpdTkState.push_back(trackState);
      } else {
	//mf::LogWarning("TrackKalmanFitter") << "WARNING: forward propagation failed. Skip this hit...";
	rejectedhsidx.push_back(ihit);
	continue;
      }
    }//for (auto hitstate : hitstatev)
  }

  assert( rejectedhsidx.size()+hitstateidx.size() == hitstatev.size());
  //std::cout << "TRACK AFTER FWD" << std::endl;
  //std::cout << "trackState parameters="   << trackState.parameters() << std::endl;
  //std::cout << "trackState covariance=\n" << trackState.covariance() << std::endl;

  //reinitialize trf for backward fit, scale the error to avoid biasing the backward fit
  trackState.setCovariance(100.*trackState.covariance());

  //backward loop over track states and hits in fwdUpdTracks: use hits for backward fit and fwd track states for smoothing
  float totChi2 = 0.;
  for (int itk = fwdPrdTkState.size()-1; itk>=0; itk--) {
    auto& fwdPrdTrackState = fwdPrdTkState[itk];
    auto& fwdUpdTrackState = fwdUpdTkState[itk];
    const auto& hitstate = hitstatev[hitstateidx[itk]];
    auto& hitflags = hitflagsv[hitstateidx[itk]];
    bool propok = true;
    trackState = propToPlane->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true);
    if (propok) {
      //combine forward predicted and backward predicted, add it to residuals (fixme)
      bool pcombok = fwdPrdTrackState.combineWithTrackState(trackState.trackState());
      if (pcombok==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	return false;
      }
      //combine forward updated and backward predicted, add it to output trajectory
      bool ucombok = fwdUpdTrackState.combineWithTrackState(trackState.trackState());
      if (ucombok==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	return false;
      }
      //update backward predicted, only if the hit was not excluded
      if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit)==0) {
	bool upok = trackState.updateWithHitState(hitstate);
	if (upok==0) {
	  mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	  return false;
	}
	//compute the chi2 between the combined predicted and the hit
	totChi2+=fwdPrdTrackState.chi2(hitstate);
      }
    } else {
      // ok, if the backward propagation failed we exclude this point from the rest of the fit,
      // but we can still use its position from the forward fit, so jsut mark it as ExcludedFromFit
      hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
      //mf::LogWarning("TrackKalmanFitter") << "WARNING: backward propagation failed. Skip this hit...";
      continue;
    }
  }//for (unsigned int itk = fwdPrdTkState.size()-1; itk>=0; itk--) {

  if (fwdUpdTkState.size()<2) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  //fill output trjectory objects with smoothed track and its hits
  std::vector<Point_t>                     positions;
  std::vector<Vector_t>                    momenta;
  std::vector<recob::TrajectoryPointFlags> flags;
  std::vector<unsigned int> hittpindex;
  if (sortOutputHitsMinLength_) {
    //try to sort fixing wires order on planes and picking the closest next plane
    std::vector<std::vector<unsigned int> > tracksInPlanes(nplanes);
    for (unsigned int p = 0; p<hitstateidx.size(); ++p) {
      const auto& hitstate = hitstatev[hitstateidx[p]];
      tracksInPlanes[hitstate.wireId().Plane].push_back(p);
    }
    //this assumes that the first hit/state is a good one, may want to check if that's the case
    std::vector<unsigned int> iterTracksInPlanes;
    for (auto it : tracksInPlanes) iterTracksInPlanes.push_back(0);
    auto& pos = fwdUpdTkState.front().position();
    auto& dir = fwdUpdTkState.front().momentum();
    for (unsigned int p = 0; p<fwdUpdTkState.size(); ++p) {
      int min_plane = -1;
      double min_dotp = DBL_MAX;
      for (unsigned int iplane = 0; iplane<iterTracksInPlanes.size(); ++iplane) {
	for (unsigned int& itk = iterTracksInPlanes[iplane]; itk<tracksInPlanes[iplane].size(); ++itk) {
	  auto& trackstate = fwdUpdTkState[tracksInPlanes[iplane][iterTracksInPlanes[iplane]]];
	  auto& tmppos = trackstate.position();
	  const double dotp = dir.Dot(tmppos-pos);
	  if (dotp<min_dotp) {
	    min_plane = iplane;
	    min_dotp = dotp;
	  }
	  break;
	}
      }
      if (min_plane<0) continue;
      const unsigned int ihit = tracksInPlanes[min_plane][iterTracksInPlanes[min_plane]];
      const auto& trackstate = fwdUpdTkState[ihit];
      const auto& hitflags = hitflagsv[hitstateidx[ihit]];
      positions.push_back(trackstate.position());
      momenta.push_back(trackstate.momentum());
      const unsigned int originalPos = (reverseHits ? hitstatev.size()-hitstateidx[ihit]-1 : hitstateidx[ihit]);
      //
      assert(originalPos>=0 && originalPos<hitstatev.size());
      //assert(hitstate.hitPtr().key()==hits[originalPos].key());
      //
      flags.push_back(recob::TrajectoryPointFlags(originalPos,hitflags));
      outHits.push_back(hits[originalPos]);
      hittpindex.push_back(p);
      iterTracksInPlanes[min_plane]++;
      //
      const auto& prdtrack = fwdPrdTkState[ihit];
      const auto& hitstate = hitstatev[hitstateidx[ihit]];
      assert(hitstate.wireId().Plane == hits[originalPos]->WireID().Plane);
      trackFitHitInfos.push_back( recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),prdtrack.parameters(),prdtrack.covariance(),hitstate.wireId()) );
      //
    }
  } else {
    for (unsigned int p = 0; p<fwdUpdTkState.size(); ++p) {
      const auto& trackstate = fwdUpdTkState[p];
      const auto& hitflags   = hitflagsv[hitstateidx[p]];
      const unsigned int originalPos = (reverseHits ? hitstatev.size()-hitstateidx[p]-1 : hitstateidx[p]);
      assert(originalPos>=0 && originalPos<hitstatev.size());
      positions.push_back(trackstate.position());
      momenta.push_back(trackstate.momentum());
      flags.push_back(recob::TrajectoryPointFlags(originalPos,hitflags));
      outHits.push_back(hits[originalPos]);
      //
      const auto& prdtrack = fwdPrdTkState[p];
      const auto& hitstate = hitstatev[hitstateidx[p]];
      assert(hitstate.wireId().Plane == hits[originalPos]->WireID().Plane);
      trackFitHitInfos.push_back( recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),prdtrack.parameters(),prdtrack.covariance(),hitstate.wireId()) );
      //
    }
  }

  for (unsigned int rejidx = 0; rejidx<rejectedhsidx.size(); ++rejidx) {
    positions.push_back(Point_t(util::kBogusD,util::kBogusD,util::kBogusD));
    momenta.push_back(Vector_t(util::kBogusD,util::kBogusD,util::kBogusD));
    const unsigned int originalPos = (reverseHits ? hitstatev.size()-rejectedhsidx[rejidx]-1 : rejectedhsidx[rejidx]);
    auto& mask = hitflagsv[rejectedhsidx[rejidx]];
    mask.set(recob::TrajectoryPointFlagTraits::HitIgnored,recob::TrajectoryPointFlagTraits::NoPoint,recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    flags.push_back(recob::TrajectoryPointFlags(originalPos,mask));
    outHits.push_back(hits[originalPos]);
    //
    const auto& hitstate = hitstatev[rejectedhsidx[rejidx]];
    SVector5 fakePar5(util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD);
    SMatrixSym55 fakeCov55;
    for (int i=0;i<5;i++) for (int j=i;j<5;j++) fakeCov55(i,j) = util::kBogusD;
    assert(hitstate.wireId().Plane == hits[originalPos]->WireID().Plane);
    trackFitHitInfos.push_back( recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),fakePar5,fakeCov55,hitstate.wireId()) );
    //
    //assert(hitstatev[rejectedhsidx[rejidx]].hitPtr().key()==hits[originalPos].key());
    //
  }

  assert(outHits.size()==hits.size());

  if (positions.size()<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  bool zeromom = false;
  for (const auto& mom : momenta) {
    if (mom.Mag2() == 0.) zeromom = true;
  }
  if (zeromom) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  bool propok = true;
  KFTrackState resultF = propToPlane->rotateToPlane(propok, fwdUpdTkState.front().trackState(),
						    Plane(fwdUpdTkState.front().position(),fwdUpdTkState.front().momentum()));
  KFTrackState resultB = propToPlane->rotateToPlane(propok, fwdUpdTkState.back().trackState(),
						    Plane(fwdUpdTkState.back().position(),fwdUpdTkState.back().momentum()));

  int ndof = positions.size()-4;//hits are 1D measurement, i.e. each hit is one d.o.f.; no B field: 4 fitted parameters

  outTrack = recob::Track(recob::TrackTrajectory(std::move(positions),std::move(momenta),std::move(flags),true),
			  pdgid,totChi2,ndof,std::move(resultF.covariance()),std::move(resultB.covariance()),track.ID());

  return true;

}
