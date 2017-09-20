#include "TrackKalmanFitter.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "larreco/RecoAlg/TrackCreationBookKeeper.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

bool trkf::TrackKalmanFitter::fitTrack(const recob::TrackTrajectory& traj, const int tkID, const SMatrixSym55& covVtx, const SMatrixSym55& covEnd,
				       const std::vector<art::Ptr<recob::Hit> >& hits, const double pval, const int pdgid, const bool flipDirection,
				       recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, trkmkr::OptionalOutputs& optionals) const {
  auto position = traj.Vertex();
  auto direction = traj.VertexDirection();

  if (flipDirection) {
    position = traj.End();
    direction = -traj.EndDirection();
  }

  auto trackStateCov = (flipDirection ? covEnd : covVtx );

  return fitTrack(position, direction, trackStateCov, hits, traj.Flags(), tkID, pval, pdgid, outTrack, outHits, optionals);
}

bool trkf::TrackKalmanFitter::fitTrack(const Point_t& position, const Vector_t& direction, SMatrixSym55& trackStateCov, 
				       const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags,
				       const int tkID, const double pval, const int pdgid,
				       recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, trkmkr::OptionalOutputs& optionals) const {
  if (hits.size()<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  // setup the KFTrackState we'll use throughout the fit
  KFTrackState trackState = setupInitialTrackState(position, direction, trackStateCov, pval, pdgid);

  // setup vector of HitStates and flags, with either same or inverse order as input hit vector
  // this is what we'll loop over during the fit
  bool reverseHits = false;
  std::vector<HitState>                            hitstatev;
  std::vector<recob::TrajectoryPointFlags::Mask_t> hitflagsv;
  bool inputok = setupInputStates(hits, flags, trackState, reverseHits, hitstatev, hitflagsv);
  if (!inputok) return false;

  // track and index vectors we use to store the fit results
  std::vector<KFTrackState> fwdPrdTkState;
  std::vector<KFTrackState> fwdUpdTkState;
  std::vector<unsigned int> hitstateidx;
  std::vector<unsigned int> rejectedhsidx;
  std::vector<unsigned int> sortedtksidx;

  // do the actual fit
  bool fitok = doFitWork(trackState, hitstatev, hitflagsv, fwdPrdTkState, fwdUpdTkState, hitstateidx, rejectedhsidx, sortedtksidx);
  if (!fitok && (skipNegProp_ || cleanZigzag_) && tryNoSkipWhenFails_) {
    mf::LogWarning("TrackKalmanFitter") << "Trying to recover with skipNegProp = false and cleanZigzag = false\n";
    fitok = doFitWork(trackState, hitstatev, hitflagsv, fwdPrdTkState, fwdUpdTkState, hitstateidx, rejectedhsidx, sortedtksidx, false);
  }
  if (!fitok) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failed for track with ID=" << tkID << "\n";
    return false;
  }

  // fill the track, the output hit vector, and the optional output products
  bool fillok = fillResult(hits, tkID, pdgid, reverseHits, hitstatev, hitflagsv, fwdPrdTkState, fwdUpdTkState, hitstateidx, rejectedhsidx, sortedtksidx, outTrack, outHits, optionals);
  return fillok;
}

trkf::KFTrackState trkf::TrackKalmanFitter::setupInitialTrackState(const Point_t& position, const Vector_t& direction, SMatrixSym55& trackStateCov, 
								   const double pval, const int pdgid) const {
  //start from large enough covariance matrix so that the fit is not biased
  if (trackStateCov==SMatrixSym55()) {
    trackStateCov(0, 0) = 1000.;
    trackStateCov(1, 1) = 1000.;
    trackStateCov(2, 2) = 0.25;
    trackStateCov(3, 3) = 0.25;
    trackStateCov(4, 4) = 10.;
  } else trackStateCov*=100.;
  // build vector of parameters on plane with point on the track and direction normal to the plane parallel to the track (so the first four parameters are zero by construction)
  SVector5 trackStatePar(0.,0.,0.,0.,1./pval);
  return KFTrackState(trackStatePar, trackStateCov, Plane(position,direction), true, pdgid);//along direction by definition
}

bool trkf::TrackKalmanFitter::setupInputStates(const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags, 
					       const KFTrackState& trackState, bool& reverseHits,
					       std::vector<HitState>& hitstatev, std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv) const {
  // figure out hit sorting based on minimum distance to first or last hit
  // (to be removed once there is a clear convention for hit sorting)
  geo::WireGeo const& wgeomF = geom->WireIDToWireGeo(hits.front()->WireID());
  geo::WireGeo const& wgeomB = geom->WireIDToWireGeo(hits.back()->WireID());
  Plane pF = recob::tracking::makePlane(wgeomF);
  Plane pB = recob::tracking::makePlane(wgeomB);
  bool success = true;
  double distF = propagator->distanceToPlane(success, trackState.trackState(), pF);
  double distB = propagator->distanceToPlane(success, trackState.trackState(), pB);
  if (dumpLevel_>1) std::cout << "distF=" << distF << " distB=" << distB << std::endl;
  if (!success) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }
  reverseHits = distB<distF;

  if (dumpLevel_>1) std::cout << "flags.size()=" << flags.size() << " hits.size()=" << hits.size() << " reverseHits=" << reverseHits << std::endl;

  // setup vector of HitStates and flags, with either same or inverse order as input hit vector
  // this is what we'll loop over during the fit
  const int fsize = flags.size();
  const int beg = (reverseHits ? hits.size()-1 : 0);
  const int end = (reverseHits ? -1 : hits.size());
  for (int ihit = beg; ihit!=end; (reverseHits ? ihit-- : ihit++)) {
    const auto& hit = hits[ihit];
    double t = hit->PeakTime();
    double terr = (useRMS_ ? hit->RMS() : hit->SigmaPeakTime() );
    double x = detprop->ConvertTicksToX(t, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
    double xerr = terr * detprop->GetXTicksCoefficient();
    hitstatev.push_back( std::move( HitState(x,hitErr2ScaleFact_*xerr*xerr,hit->WireID(),geom->WireIDToWireGeo(hit->WireID())) ) );
    //
    if (fsize>0 && ihit<fsize) hitflagsv.push_back( flags[ihit].mask() );
    else hitflagsv.push_back(recob::TrajectoryPointFlags::DefaultFlagsMask());
    //
    if (dumpLevel_>2) std::cout << "pushed flag mask=" << hitflagsv.back() << std::endl;
    //
    if (rejectHighMultHits_ && hit->Multiplicity()>1)   {
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::Merged);
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    }
    if (rejectHitsNegativeGOF_ && hit->GoodnessOfFit()<0) {
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::Suspicious);
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    }
  }
  assert(hits.size()==hitstatev.size());
  return true;
}

bool trkf::TrackKalmanFitter::doFitWork(KFTrackState& trackState, std::vector<HitState>& hitstatev, std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
					std::vector<KFTrackState>& fwdPrdTkState, std::vector<KFTrackState>& fwdUpdTkState, 
					std::vector<unsigned int>& hitstateidx, std::vector<unsigned int>& rejectedhsidx, std::vector<unsigned int>& sortedtksidx, 
					bool applySkipClean) const {

  fwdPrdTkState.clear();
  fwdUpdTkState.clear();
  hitstateidx.clear();
  rejectedhsidx.clear();
  sortedtksidx.clear();
  // these three vectors are aligned
  fwdPrdTkState.reserve(hitstatev.size());
  fwdUpdTkState.reserve(hitstatev.size());
  hitstateidx.reserve(hitstatev.size());

  if (sortHitsByPlane_) {
    //array of hit indices in planes, keeping the original sorting by plane
    const unsigned int nplanes = geom->MaxPlanes();
    std::vector<std::vector<unsigned int> > hitsInPlanes(nplanes);
    for (unsigned int ihit = 0; ihit<hitstatev.size(); ihit++) {
      hitsInPlanes[hitstatev[ihit].wireId().Plane].push_back(ihit);
    }
    //array of indices, where iterHitsInPlanes[i] is the iterator over hitsInPlanes[i]
    std::vector<unsigned int> iterHitsInPlanes(nplanes,0);
    for (unsigned int p = 0; p<hitstatev.size(); ++p) {
      if (dumpLevel_>1) std::cout << std::endl << "processing hit #" << p << std::endl;
      if (dumpLevel_>1) {
	std::cout << "compute distance from state=" << std::endl; trackState.dump();
      }
      int min_plane = -1;
      double min_dist = DBL_MAX;
      //find the closest hit according to the sorting in each plane
      for (unsigned int iplane = 0; iplane<nplanes; ++iplane) {
	//note: ih is a reference, so when 'continue' we increment iterHitsInPlanes[iplane] and the hit is effectively rejected
	for (unsigned int& ih = iterHitsInPlanes[iplane]; ih<hitsInPlanes[iplane].size(); ++ih) {
	  if (dumpLevel_>1) std::cout << "iplane=" << iplane << " nplanes=" << nplanes << " iterHitsInPlanes[iplane]=" << iterHitsInPlanes[iplane] << " hitsInPlanes[iplane].size()=" << hitsInPlanes[iplane].size() << std::endl;
 	  unsigned int ihit = hitsInPlanes[iplane][ih];
	  const auto& hitstate = hitstatev[ihit];
	  const auto& hitflags = hitflagsv[ihit];
	  if (hitflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint        ) ||
	      hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ||
	      hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected       )) {
	    if (dumpLevel_>1) std::cout << "rejecting hit flags=" << hitflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint        ) << ", " << hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) << ", " << hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected       ) << std::endl;
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  //get distance to measurement surface
	  bool success = true;
	  const double dist = propagator->distanceToPlane(success, trackState.trackState(), hitstate.plane());
	  if (dumpLevel_>1) std::cout << "distance to plane " << iplane << " wire=" << hitstate.wireId().Wire << " = " << dist << ", min_dist=" << min_dist << " min_plane=" << min_plane << " success=" << success << std::endl;
	  if (!success) {
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  if (applySkipClean && skipNegProp_ && dist<0.) {
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
      if (dumpLevel_>1) std::cout << "pick min_dist=" << min_dist << " min_plane=" << min_plane << std::endl;
      //now we know which is the closest wire: get the hitstate and increment the iterator
      if (min_plane<0) continue;
      const unsigned int ihit = hitsInPlanes[min_plane][iterHitsInPlanes[min_plane]];
      const auto& hitstate = hitstatev[ihit];
      if (dumpLevel_>1) hitstate.dump();
      auto& hitflags = hitflagsv[ihit];
      iterHitsInPlanes[min_plane]++;
      //propagate to measurement surface
      bool propok = true;
      trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::FORWARD);
      if (!propok && !(applySkipClean && skipNegProp_)) trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::BACKWARD);
      if (dumpLevel_>1) {
	std::cout << "hit state " << std::endl; hitstate.dump();
	std::cout << "propagation result=" << propok << std::endl;
	std::cout << "propagated state " << std::endl; trackState.dump();
	std::cout << "propagated planarity=" << hitstate.plane().direction().Dot(hitstate.plane().position()-trackState.position()) << std::endl;
      }
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
	//now update the forward fitted track
	bool upok = trackState.updateWithHitState(hitstate);
	if (upok==0) {
	  mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	  return false;
	}
	if (dumpLevel_>1) {
	  std::cout << "updated state " << std::endl; trackState.dump();
	}
	fwdUpdTkState.push_back(trackState);
      } else {
	if (dumpLevel_>0) std::cout << "WARNING: forward propagation failed. Skip this hit..." << std::endl;
	rejectedhsidx.push_back(ihit);
	continue;
      }
    }
  } else {
    for (unsigned int ihit=0; ihit<hitstatev.size(); ++ihit) {
      const auto& hitstate = hitstatev[ihit];
      if (dumpLevel_>1) {
	std::cout << std::endl << "processing hit #" << ihit << std::endl;
	hitstate.dump();
      }
      auto& hitflags = hitflagsv[ihit];
      if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ||
	  hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected)) {
	rejectedhsidx.push_back(ihit);
	continue;
      }
      if (applySkipClean && skipNegProp_) {
	bool success = true;
	const double dist = propagator->distanceToPlane(success, trackState.trackState(), hitstate.plane());
	if (dist<0. || success==false) {
	  if (dumpLevel_>0) std::cout << "WARNING: negative propagation distance. Skip this hit..." << std::endl;;
	  rejectedhsidx.push_back(ihit);
	  continue;
	}
      }
      //propagate to measurement surface
      bool propok = true;
      trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::FORWARD);
      if (!propok && !(applySkipClean && skipNegProp_)) trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::BACKWARD);
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
	bool upok = trackState.updateWithHitState(hitstate);
	if (upok==0) {
	  mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	  return false;
	}
	fwdUpdTkState.push_back(trackState);
      } else {
	if (dumpLevel_>0) std::cout << "WARNING: forward propagation failed. Skip this hit..." << std::endl;;
	rejectedhsidx.push_back(ihit);
	continue;
      }
    }//for (auto hitstate : hitstatev)
  }

  assert( rejectedhsidx.size()+hitstateidx.size() == hitstatev.size() );
  if (dumpLevel_>0) {
    std::cout << "TRACK AFTER FWD" << std::endl;
    trackState.dump();
  }

  //reinitialize trf for backward fit, scale the error to avoid biasing the backward fit
  trackState.setCovariance(100.*trackState.covariance());

  //backward loop over track states and hits in fwdUpdTracks: use hits for backward fit and fwd track states for smoothing
  int nvalidhits = 0;
  for (int itk = fwdPrdTkState.size()-1; itk>=0; itk--) {
    auto& fwdPrdTrackState = fwdPrdTkState[itk];
    auto& fwdUpdTrackState = fwdUpdTkState[itk];
    const auto& hitstate = hitstatev[hitstateidx[itk]];
    auto& hitflags = hitflagsv[hitstateidx[itk]];
    if (dumpLevel_>1) {
      std::cout << std::endl << "processing backward hit #" << itk << std::endl;
      hitstate.dump();
    }
    bool propok = true;
    trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::BACKWARD);
    if (!propok) trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::FORWARD);//do we want this?
    //
    if (dumpLevel_>1) {
      std::cout << "propagation result=" << propok << std::endl;
      std::cout << "propagated state " << std::endl; trackState.dump();
      std::cout << "propagated planarity=" << hitstate.plane().direction().Dot(hitstate.plane().position()-trackState.position()) << std::endl;
    }
    if (propok) {
      //combine forward predicted and backward predicted
      bool pcombok = fwdPrdTrackState.combineWithTrackState(trackState.trackState());
      if (pcombok==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	return false;
      }
      if (dumpLevel_>1) {
	std::cout << "combined prd state " << std::endl; fwdPrdTrackState.dump();
      }
      //combine forward updated and backward predicted
      bool ucombok = fwdUpdTrackState.combineWithTrackState(trackState.trackState());
      if (ucombok==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	return false;
      }
      if (dumpLevel_>1) {
	std::cout << "combined upd state " << std::endl; fwdUpdTrackState.dump();
      }
      //update backward predicted, only if the hit was not excluded
      if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit)==0) {
	bool upok = trackState.updateWithHitState(hitstate);
	if (upok==0) {
	  mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
	  return false;
	}
	if (dumpLevel_>1) {
	  std::cout << "updated state " << std::endl; trackState.dump();
	}
	//
	nvalidhits++;
      }
    } else {
      // ok, if the backward propagation failed we exclude this point from the rest of the fit,
      // but we can still use its position from the forward fit, so just mark it as ExcludedFromFit
      hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
      if (dumpLevel_>0) std::cout << "WARNING: backward propagation failed. Skip this hit..." << std::endl;;
      continue;
    }
  }//for (unsigned int itk = fwdPrdTkState.size()-1; itk>=0; itk--) {

  if (nvalidhits<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  // sort output states
  sortOutput(hitstatev, fwdUpdTkState, hitstateidx, rejectedhsidx, sortedtksidx, applySkipClean);
  if (sortedtksidx.size()<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  assert( rejectedhsidx.size()+sortedtksidx.size() == hitstatev.size() );
  return true;
}

void trkf::TrackKalmanFitter::sortOutput(std::vector<HitState>& hitstatev, std::vector<KFTrackState>& fwdUpdTkState, 
					 std::vector<unsigned int>& hitstateidx, std::vector<unsigned int>& rejectedhsidx,
					 std::vector<unsigned int>& sortedtksidx, bool applySkipClean) const {
  //
  if (sortOutputHitsMinLength_) {
    //try to sort fixing wires order on planes and picking the closest next plane
    const unsigned int nplanes = geom->MaxPlanes();
    std::vector<std::vector<unsigned int> > tracksInPlanes(nplanes);
    for (unsigned int p = 0; p<hitstateidx.size(); ++p) {
      const auto& hitstate = hitstatev[hitstateidx[p]];
      tracksInPlanes[hitstate.wireId().Plane].push_back(p);
    }
    //this assumes that the first hit/state is a good one, may want to check if that's the case
    std::vector<unsigned int> iterTracksInPlanes;
    for (auto it : tracksInPlanes) iterTracksInPlanes.push_back(0);
    auto pos = fwdUpdTkState.front().position();
    auto dir = fwdUpdTkState.front().momentum();
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
      if (applySkipClean && skipNegProp_ && min_dotp<0.) {
	rejectedhsidx.push_back(hitstateidx[ihit]);
	iterTracksInPlanes[min_plane]++;
	continue;
      }
      auto& trackstate = fwdUpdTkState[ihit];
      pos = trackstate.position();
      dir = trackstate.momentum();
      //
      sortedtksidx.push_back(ihit);
      iterTracksInPlanes[min_plane]++;
    }
  } else {
    for (unsigned int p = 0; p<fwdUpdTkState.size(); ++p) {
      sortedtksidx.push_back(p);
    }
  }
  //
  if (applySkipClean && cleanZigzag_) {
    std::vector<unsigned int> itoerase;
    bool clean = false;
    while (!clean) {
      bool broken = false;
      auto pos0 = fwdUpdTkState[sortedtksidx[0]].position();
      unsigned int i=1;
      unsigned int end=sortedtksidx.size()-1;
      for (;i<end;++i) {
	auto dir0 = fwdUpdTkState[sortedtksidx[i]].position()-pos0;
	auto dir2 = fwdUpdTkState[sortedtksidx[i+1]].position()-fwdUpdTkState[sortedtksidx[i]].position();
	dir0/=dir0.R();
	dir2/=dir2.R();
	if (dir2.Dot(dir0)<0.) {
	  broken = true;
	  end--;
	  break;
	} else pos0 = fwdUpdTkState[sortedtksidx[i]].position();
      }
      if (!broken) {
	clean = true;
      } else {
	rejectedhsidx.push_back(hitstateidx[sortedtksidx[i]]);
	sortedtksidx.erase(sortedtksidx.begin()+i);
      }
    }
  }
  //
}


bool trkf::TrackKalmanFitter::fillResult(const std::vector<art::Ptr<recob::Hit> >& inHits, const int tkID, const int pdgid, const bool reverseHits,
					 std::vector<HitState>& hitstatev, std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
					 std::vector<KFTrackState>& fwdPrdTkState, std::vector<KFTrackState>& fwdUpdTkState,
					 std::vector<unsigned int>& hitstateidx, std::vector<unsigned int>& rejectedhsidx, std::vector<unsigned int>& sortedtksidx,
					 recob::Track& outTrack, std::vector<art::Ptr<recob::Hit> >& outHits, trkmkr::OptionalOutputs& optionals) const {
  // fill output trajectory objects with smoothed track and its hits
  trkmkr::TrackCreationBookKeeper tcbk(outHits, optionals, tkID, pdgid, true);
  for (unsigned int p : sortedtksidx) {
    const auto& trackstate = fwdUpdTkState[p];
    const auto& hitflags   = hitflagsv[hitstateidx[p]];
    const unsigned int originalPos = (reverseHits ? hitstatev.size()-hitstateidx[p]-1 : hitstateidx[p]);
    assert(originalPos>=0 && originalPos<hitstatev.size());
    //
    const auto& prdtrack = fwdPrdTkState[p];
    const auto& hitstate = hitstatev[hitstateidx[p]];
    assert(hitstate.wireId().Plane == inHits[originalPos]->WireID().Plane);
    //
    trkmkr::OptionalPointElement ope;
    if (optionals.isTrackFitInfosInit()) {
      ope.setTrackFitHitInfo(recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),prdtrack.parameters(),prdtrack.covariance(),hitstate.wireId()));
    }
    tcbk.addPoint(trackstate.position(), trackstate.momentum(), inHits[originalPos],
		  recob::TrajectoryPointFlags(originalPos,hitflags), prdtrack.chi2(hitstate), ope);
  }

  // fill also with rejected hits information
  SMatrixSym55 fakeCov55;
  for (int i=0;i<5;i++) for (int j=i;j<5;j++) fakeCov55(i,j) = util::kBogusD;
  for (unsigned int rejidx = 0; rejidx<rejectedhsidx.size(); ++rejidx) {
    const unsigned int originalPos = (reverseHits ? hitstatev.size()-rejectedhsidx[rejidx]-1 : rejectedhsidx[rejidx]);
    auto& mask = hitflagsv[rejectedhsidx[rejidx]];
    mask.set(recob::TrajectoryPointFlagTraits::HitIgnored,recob::TrajectoryPointFlagTraits::NoPoint);
    if (mask.isSet(recob::TrajectoryPointFlagTraits::Rejected)==0) mask.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    //
    const auto& hitstate = hitstatev[rejectedhsidx[rejidx]];
    assert(hitstate.wireId().Plane == inHits[originalPos]->WireID().Plane);
    trkmkr::OptionalPointElement ope;
    if (optionals.isTrackFitInfosInit()) {
      ope.setTrackFitHitInfo( recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),
						     SVector5(util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD),fakeCov55,hitstate.wireId()) );
    }
    tcbk.addPoint(Point_t(util::kBogusD,util::kBogusD,util::kBogusD), Vector_t(util::kBogusD,util::kBogusD,util::kBogusD), inHits[originalPos],
		  recob::TrajectoryPointFlags(originalPos,mask), 0, ope);
  }

  if (dumpLevel_>1) std::cout << "outHits.size()=" << outHits.size() << " inHits.size()=" << inHits.size() << std::endl;
  assert(outHits.size()==inHits.size());

  bool propok = true;
  KFTrackState resultF = propagator->rotateToPlane(propok, fwdUpdTkState[sortedtksidx.front()].trackState(),
						   Plane(fwdUpdTkState[sortedtksidx.front()].position(),fwdUpdTkState[sortedtksidx.front()].momentum()));
  KFTrackState resultB = propagator->rotateToPlane(propok, fwdUpdTkState[sortedtksidx.back()].trackState(),
						   Plane(fwdUpdTkState[sortedtksidx.back()].position(),fwdUpdTkState[sortedtksidx.back()].momentum()));

  outTrack = tcbk.finalizeTrack(SMatrixSym55(resultF.covariance()),SMatrixSym55(resultB.covariance()));

  //if there are points with zero momentum return false
  size_t point = 0;
  while (outTrack.HasValidPoint(point)) {
    if (outTrack.MomentumAtPoint( outTrack.NextValidPoint(point++) ) <= 1.0e-9) return false;
  }

  if (dumpLevel_>0) {
    std::cout << "outTrack vertex=" << outTrack.Start()
	      << "\ndir=" << outTrack.StartDirection()
	      << "\ncov=\n" << outTrack.StartCovariance()
	      << "\nlength=" << outTrack.Length() //<< " inLenght=" << track.Length()
	      << std::endl;
  }

  return true;
}
