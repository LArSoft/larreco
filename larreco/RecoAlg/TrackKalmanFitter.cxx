#include "TrackKalmanFitter.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/RecoObjects/KFTrackState.h"
#include "lardata/RecoObjects/TrackingPlaneHelper.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

bool trkf::TrackKalmanFitter::fitTrack(const recob::Trajectory& track, int tkID,
				       const SMatrixSym55& covVtx, const SMatrixSym55& covEnd,
				       const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags,
				       const double pval, const int pdgid, const bool flipDirection,
				       recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits,
				       std::vector<recob::TrackFitHitInfo>& trackFitHitInfos) {
  outHits.clear();
  if (hits.size()<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  auto position = track.Vertex();
  auto direction = track.VertexDirection();

  if (flipDirection) {
    position = track.End();
    direction = -track.EndDirection();
  }

  auto trackStateCov = (flipDirection ? covEnd : covVtx );

  return fitTrack(position, direction, trackStateCov, tkID, hits, flags, pval, pdgid, outTrack, outHits, trackFitHitInfos);

}

bool trkf::TrackKalmanFitter::fitTrack(const Point_t& position, const Vector_t& direction,
				       SMatrixSym55& trackStateCov, int tkID,
				       const std::vector<art::Ptr<recob::Hit> >& hits, const std::vector<recob::TrajectoryPointFlags>& flags,
				       const double pval, const int pdgid,
				       recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits,
				       std::vector<recob::TrackFitHitInfo>& trackFitHitInfos) {
  outHits.clear();
  if (hits.size()<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  if (trackStateCov==SMatrixSym55()) {
    trackStateCov(0, 0) = 1000.;
    trackStateCov(1, 1) = 1000.;
    trackStateCov(2, 2) = 0.25;
    trackStateCov(3, 3) = 0.25;
    trackStateCov(4, 4) = 10.;
  }
  SVector5 trackStatePar(0.,0.,0.,0.,1./pval);

  // setup the KFTrackState we'll use throughout the fit
  KFTrackState trackState(trackStatePar, trackStateCov, Plane(position,direction), true, pdgid);//along direction by definition

  // figure out hit sorting based on minimum distance to first or last hit
  // (to be removed once there is a clear convention for hit sorting)
  geo::WireGeo const& wgeomF = geom->WireIDToWireGeo(hits.front()->WireID());
  geo::WireGeo const& wgeomB = geom->WireIDToWireGeo(hits.back()->WireID());
  Plane pF = recob::tracking::makePlane(wgeomF);
  Plane pB = recob::tracking::makePlane(wgeomB);
  bool success = true;
  double distF = propagator->distanceToPlane(success, position, direction, pF);
  double distB = propagator->distanceToPlane(success, position, direction, pB);
  if (dumpLevel_>1) std::cout << "distF=" << distF << " distB=" << distB << std::endl;
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
    hitstatev.push_back( std::move( HitState(x,hitErr2ScaleFact_*xerr*xerr,hit->WireID(),geom->WireIDToWireGeo(hit->WireID())) ) );
    //
    if (flags.size()>0) hitflagsv.push_back( flags[ihit].mask() );
    else hitflagsv.push_back(recob::TrajectoryPointFlags::makeMask());
    //
    if (rejectHighMultHits_ && hit->Multiplicity()>1)   {
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::Merged);
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    }
    if (rejectHitsNegativeGOF_ && hit->GoodnessOfFit()<0) {
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::Suspicious);
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    }
    //
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
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  //get distance to measurement surface
	  bool success = true;
	  const double dist = propagator->distanceToPlane(success, trackState.trackState(), hitstate.plane());
	  if (!success) {
	    rejectedhsidx.push_back(ihit);
	    continue;
	  }
	  if (dumpLevel_>1) std::cout << "distance to plane " << iplane << " wire=" << hitstate.wireId().Wire << " = " << dist << ", min_dist=" << min_dist << " min_plane=" << min_plane << std::endl;
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
      if (!propok && !skipNegProp_) trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::BACKWARD);
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
      if (skipNegProp_) {
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
      if (!propok && !skipNegProp_) trackState = propagator->propagateToPlane(propok, trackState.trackState(), hitstate.plane(), true, true, TrackStatePropagator::BACKWARD);
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

  assert( rejectedhsidx.size()+hitstateidx.size() == hitstatev.size());
  if (dumpLevel_>0) {
    std::cout << "TRACK AFTER FWD" << std::endl;
    trackState.dump();
  }

  //reinitialize trf for backward fit, scale the error to avoid biasing the backward fit
  trackState.setCovariance(100.*trackState.covariance());

  //backward loop over track states and hits in fwdUpdTracks: use hits for backward fit and fwd track states for smoothing
  float totChi2 = 0.;
  int nchi2 = 0;
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
	//compute the chi2 between the combined predicted and the hit
	totChi2+=fwdPrdTrackState.chi2(hitstate);
	nchi2++;
      }
    } else {
      // ok, if the backward propagation failed we exclude this point from the rest of the fit,
      // but we can still use its position from the forward fit, so just mark it as ExcludedFromFit
      hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
      if (dumpLevel_>0) std::cout << "WARNING: backward propagation failed. Skip this hit..." << std::endl;;
      continue;
    }
  }//for (unsigned int itk = fwdPrdTkState.size()-1; itk>=0; itk--) {

  if (nchi2<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  //fill output trajectory objects with smoothed track and its hits
  std::vector<Point_t>                     positions;
  std::vector<Vector_t>                    momenta;
  std::vector<recob::TrajectoryPointFlags> outFlags;
  std::vector<unsigned int> hittpindex;
  std::vector<unsigned int> updstatesindex;
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
      if (skipNegProp_ && min_dotp<0.) {
	rejectedhsidx.push_back(hitstateidx[ihit]);
	iterTracksInPlanes[min_plane]++;
	continue;
      }
      const auto& trackstate = fwdUpdTkState[ihit];
      const auto& hitflags = hitflagsv[hitstateidx[ihit]];
      positions.push_back(trackstate.position());
      momenta.push_back(trackstate.momentum());
      pos = trackstate.position();
      dir = trackstate.momentum();
      const unsigned int originalPos = (reverseHits ? hitstatev.size()-hitstateidx[ihit]-1 : hitstateidx[ihit]);
      //
      assert(originalPos>=0 && originalPos<hitstatev.size());
      //
      outFlags.push_back(recob::TrajectoryPointFlags(originalPos,hitflags));
      outHits.push_back(hits[originalPos]);
      hittpindex.push_back(p);
      updstatesindex.push_back(ihit);
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
      outFlags.push_back(recob::TrajectoryPointFlags(originalPos,hitflags));
      outHits.push_back(hits[originalPos]);
      updstatesindex.push_back(p);
      //
      const auto& prdtrack = fwdPrdTkState[p];
      const auto& hitstate = hitstatev[hitstateidx[p]];
      assert(hitstate.wireId().Plane == hits[originalPos]->WireID().Plane);
      trackFitHitInfos.push_back( recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),prdtrack.parameters(),prdtrack.covariance(),hitstate.wireId()) );
      //
    }
  }

  unsigned int nrej = 0;
  if (cleanZigzag_) {
    std::vector<unsigned int> itoerase;
    bool clean = false;
    while (!clean) {
      bool broken = false;
      auto pos0 = positions[0];
      unsigned int i=1;
      unsigned int end=positions.size()-1;
      for (;i<end;++i) {
	auto dir0 = positions[i]-pos0;
	auto dir2 = positions[i+1]-positions[i];
	dir0/=dir0.R();
	dir2/=dir2.R();
	if (dir2.Dot(dir0)<0.) {
	  broken = true;
	  end--;
	  break;
	} else pos0 = positions[i];
      }
      if (!broken) {
	clean = true;
      } else {
	nrej++;
	auto pos = i;
	positions.erase(positions.begin()+pos);
	momenta.erase(momenta.begin()+pos);
	auto mask = outFlags[pos].mask();
	auto fhit = outFlags[pos].fromHit();
	outFlags.erase(outFlags.begin()+pos);
	auto hit = outHits[pos];
	outHits.erase(outHits.begin()+pos);
	auto info = trackFitHitInfos[pos];
	trackFitHitInfos.erase(trackFitHitInfos.begin()+pos);
	//
	positions.push_back(Point_t(util::kBogusD,util::kBogusD,util::kBogusD));
	momenta.push_back(Vector_t(util::kBogusD,util::kBogusD,util::kBogusD));
	mask.set(recob::TrajectoryPointFlagTraits::HitIgnored,recob::TrajectoryPointFlagTraits::NoPoint);
	if (mask.isSet(recob::TrajectoryPointFlagTraits::Rejected)==0) mask.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
	outFlags.push_back(recob::TrajectoryPointFlags(fhit,mask));
	outHits.push_back(hit);
	SVector5 fakePar5(util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD);
	SMatrixSym55 fakeCov55;
	for (int i=0;i<5;i++) for (int j=i;j<5;j++) fakeCov55(i,j) = util::kBogusD;
	trackFitHitInfos.push_back(recob::TrackFitHitInfo(info.hitMeas(),info.hitMeasErr2(),fakePar5,fakeCov55,info.WireId()));
      }
    }
  }

  if (positions.size()-nrej<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  for (unsigned int rejidx = 0; rejidx<rejectedhsidx.size(); ++rejidx) {
    positions.push_back(Point_t(util::kBogusD,util::kBogusD,util::kBogusD));
    momenta.push_back(Vector_t(util::kBogusD,util::kBogusD,util::kBogusD));
    const unsigned int originalPos = (reverseHits ? hitstatev.size()-rejectedhsidx[rejidx]-1 : rejectedhsidx[rejidx]);
    auto& mask = hitflagsv[rejectedhsidx[rejidx]];
    mask.set(recob::TrajectoryPointFlagTraits::HitIgnored,recob::TrajectoryPointFlagTraits::NoPoint);
    if (mask.isSet(recob::TrajectoryPointFlagTraits::Rejected)==0) mask.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    outFlags.push_back(recob::TrajectoryPointFlags(originalPos,mask));
    outHits.push_back(hits[originalPos]);
    //
    const auto& hitstate = hitstatev[rejectedhsidx[rejidx]];
    SVector5 fakePar5(util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD,util::kBogusD);
    SMatrixSym55 fakeCov55;
    for (int i=0;i<5;i++) for (int j=i;j<5;j++) fakeCov55(i,j) = util::kBogusD;
    assert(hitstate.wireId().Plane == hits[originalPos]->WireID().Plane);
    trackFitHitInfos.push_back( recob::TrackFitHitInfo(hitstate.hitMeas(),hitstate.hitMeasErr2(),fakePar5,fakeCov55,hitstate.wireId()) );
  }

  if (dumpLevel_>1) std::cout << "outHits.size()=" << outHits.size() << " hits.size()=" << hits.size() << std::endl;
  assert(outHits.size()==hits.size());

  bool zeromom = false;
  for (const auto& mom : momenta) {
    if (mom.Mag2() == 0.) zeromom = true;
  }
  if (zeromom) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  bool propok = true;
  KFTrackState resultF = propagator->rotateToPlane(propok, fwdUpdTkState[updstatesindex.front()].trackState(),
						   Plane(fwdUpdTkState[updstatesindex.front()].position(),fwdUpdTkState[updstatesindex.front()].momentum()));
  KFTrackState resultB = propagator->rotateToPlane(propok, fwdUpdTkState[updstatesindex.back()].trackState(),
						   Plane(fwdUpdTkState[updstatesindex.back()].position(),fwdUpdTkState[updstatesindex.back()].momentum()));

  int ndof = nchi2-4;//hits are 1D measurement, i.e. each hit is one d.o.f.; no B field: 4 fitted parameters

  outTrack = recob::Track(recob::TrackTrajectory(std::move(positions),std::move(momenta),std::move(outFlags),true),
			  pdgid,totChi2,ndof,std::move(resultF.covariance()),std::move(resultB.covariance()),tkID);

  if (dumpLevel_>0) {
    std::cout << "outTrack vertex=" << outTrack.Start()
	      << "\ndir=" << outTrack.StartDirection()
	      << "\ncov=\n" << outTrack.StartCovariance()
	      << "\nlength=" << outTrack.Length() //<< " inLenght=" << track.Length()
	      << std::endl;
  }

  return true;

}
