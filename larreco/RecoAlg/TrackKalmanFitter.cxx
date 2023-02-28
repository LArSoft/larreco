#include "TrackKalmanFitter.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <float.h>
#include <iostream>
#include <stddef.h>
#include <string>

#include "Math/BinaryOperators.h"
#include "Math/Expression.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/MatrixRepresentationsStatic.h"
#include "Math/SMatrix.h"

#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/Utilities/BitMask.h"
#include "larreco/RecoAlg/TrackCreationBookKeeper.h"
#include "larreco/TrackFinder/TrackMaker.h"

bool trkf::TrackKalmanFitter::fitTrack(detinfo::DetectorPropertiesData const& detProp,
                                       const recob::TrackTrajectory& traj,
                                       const int tkID,
                                       const SMatrixSym55& covVtx,
                                       const SMatrixSym55& covEnd,
                                       const std::vector<art::Ptr<recob::Hit>>& hits,
                                       const double pval,
                                       const int pdgid,
                                       const bool flipDirection,
                                       recob::Track& outTrack,
                                       std::vector<art::Ptr<recob::Hit>>& outHits,
                                       trkmkr::OptionalOutputs& optionals) const
{
  auto position = traj.Vertex();
  auto direction = traj.VertexDirection();

  if (tryBothDirs_) {
    recob::Track fwdTrack;
    std::vector<art::Ptr<recob::Hit>> fwdHits;
    trkmkr::OptionalOutputs fwdoptionals;
    SMatrixSym55 fwdcov = covVtx;
    bool okfwd = fitTrack(detProp,
                          position,
                          direction,
                          fwdcov,
                          hits,
                          traj.Flags(),
                          tkID,
                          pval,
                          pdgid,
                          fwdTrack,
                          fwdHits,
                          fwdoptionals);

    recob::Track bwdTrack;
    std::vector<art::Ptr<recob::Hit>> bwdHits;
    trkmkr::OptionalOutputs bwdoptionals;
    SMatrixSym55 bwdcov = covEnd;
    bool okbwd = fitTrack(detProp,
                          position,
                          -direction,
                          bwdcov,
                          hits,
                          traj.Flags(),
                          tkID,
                          pval,
                          pdgid,
                          bwdTrack,
                          bwdHits,
                          bwdoptionals);

    if (okfwd == false && okbwd == false) { return false; }
    else if (okfwd == true && okbwd == true) {
      if ((fwdTrack.CountValidPoints() / (fwdTrack.Length() * fwdTrack.Chi2PerNdof())) >=
          (bwdTrack.CountValidPoints() / (bwdTrack.Length() * bwdTrack.Chi2PerNdof()))) {
        outTrack = fwdTrack;
        outHits = fwdHits;
        optionals = std::move(fwdoptionals);
      }
      else {
        outTrack = bwdTrack;
        outHits = bwdHits;
        optionals = std::move(bwdoptionals);
      }
    }
    else if (okfwd == true) {
      outTrack = fwdTrack;
      outHits = fwdHits;
      optionals = std::move(fwdoptionals);
    }
    else {
      outTrack = bwdTrack;
      outHits = bwdHits;
      optionals = std::move(bwdoptionals);
    }
    return true;
  }
  else {
    if (flipDirection) {
      position = traj.End();
      direction = -traj.EndDirection();
    }

    auto trackStateCov = (flipDirection ? covEnd : covVtx);

    return fitTrack(detProp,
                    position,
                    direction,
                    trackStateCov,
                    hits,
                    traj.Flags(),
                    tkID,
                    pval,
                    pdgid,
                    outTrack,
                    outHits,
                    optionals);
  }
}

bool trkf::TrackKalmanFitter::fitTrack(detinfo::DetectorPropertiesData const& detProp,
                                       const Point_t& position,
                                       const Vector_t& direction,
                                       SMatrixSym55& trackStateCov,
                                       const std::vector<art::Ptr<recob::Hit>>& hits,
                                       const std::vector<recob::TrajectoryPointFlags>& flags,
                                       const int tkID,
                                       const double pval,
                                       const int pdgid,
                                       recob::Track& outTrack,
                                       std::vector<art::Ptr<recob::Hit>>& outHits,
                                       trkmkr::OptionalOutputs& optionals) const
{
  if (dumpLevel_ > 1)
    std::cout << "Fitting track with tkID=" << tkID << " start pos=" << position
              << " dir=" << direction << " nHits=" << hits.size() << " mom=" << pval
              << " pdg=" << pdgid << std::endl;
  if (hits.size() < 4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  // setup the KFTrackState we'll use throughout the fit
  KFTrackState trackState = setupInitialTrackState(position, direction, trackStateCov, pval, pdgid);

  // setup vector of HitStates and flags, with either same or inverse order as input hit vector
  // this is what we'll loop over during the fit
  std::vector<HitState> hitstatev;
  std::vector<recob::TrajectoryPointFlags::Mask_t> hitflagsv;
  bool inputok = setupInputStates(detProp, hits, flags, trackState, hitstatev, hitflagsv);
  if (!inputok) return false;

  // track and index vectors we use to store the fit results
  std::vector<KFTrackState> fwdPrdTkState;
  std::vector<KFTrackState> fwdUpdTkState;
  std::vector<unsigned int> hitstateidx;
  std::vector<unsigned int> rejectedhsidx;
  std::vector<unsigned int> sortedtksidx;

  // do the actual fit
  bool fitok = doFitWork(trackState,
                         detProp,
                         hitstatev,
                         hitflagsv,
                         fwdPrdTkState,
                         fwdUpdTkState,
                         hitstateidx,
                         rejectedhsidx,
                         sortedtksidx);
  if (!fitok && (skipNegProp_ || cleanZigzag_) && tryNoSkipWhenFails_) {
    mf::LogWarning("TrackKalmanFitter")
      << "Trying to recover with skipNegProp = false and cleanZigzag = false\n";
    fitok = doFitWork(trackState,
                      detProp,
                      hitstatev,
                      hitflagsv,
                      fwdPrdTkState,
                      fwdUpdTkState,
                      hitstateidx,
                      rejectedhsidx,
                      sortedtksidx,
                      false);
  }
  if (!fitok) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failed for track with ID=" << tkID << "\n";
    return false;
  }

  // fill the track, the output hit vector, and the optional output products
  bool fillok = fillResult(hits,
                           tkID,
                           pdgid,
                           hitstatev,
                           hitflagsv,
                           fwdPrdTkState,
                           fwdUpdTkState,
                           hitstateidx,
                           rejectedhsidx,
                           sortedtksidx,
                           outTrack,
                           outHits,
                           optionals);
  return fillok;
}

trkf::KFTrackState trkf::TrackKalmanFitter::setupInitialTrackState(const Point_t& position,
                                                                   const Vector_t& direction,
                                                                   SMatrixSym55& trackStateCov,
                                                                   const double pval,
                                                                   const int pdgid) const
{
  //start from large enough covariance matrix so that the fit is not biased
  if (trackStateCov == SMatrixSym55()) {
    trackStateCov(0, 0) = 1000.;
    trackStateCov(1, 1) = 1000.;
    trackStateCov(2, 2) = 0.25;
    trackStateCov(3, 3) = 0.25;
    trackStateCov(4, 4) = 10.;
  }
  else
    trackStateCov *= 100.;
  // build vector of parameters on plane with point on the track and direction normal to the plane parallel to the track (so the first four parameters are zero by construction)
  SVector5 trackStatePar(0., 0., 0., 0., 1. / pval);
  return KFTrackState(trackStatePar,
                      trackStateCov,
                      Plane(position, direction),
                      true,
                      pdgid); //along direction by definition
}

bool trkf::TrackKalmanFitter::setupInputStates(
  detinfo::DetectorPropertiesData const& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<recob::TrajectoryPointFlags>& flags,
  const KFTrackState& trackState,
  std::vector<HitState>& hitstatev,
  std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv) const
{
  if (dumpLevel_ > 1)
    std::cout << "flags.size()=" << flags.size() << " hits.size()=" << hits.size() << std::endl;

  // setup vector of HitStates and flags, with either same or inverse order as input hit vector
  // this is what we'll loop over during the fit
  const size_t fsize = flags.size();
  for (size_t ihit = 0; ihit != hits.size(); ihit++) {
    const auto& hit = hits[ihit];
    double t = hit->PeakTime();
    double terr = (useRMS_ ? hit->RMS() : hit->SigmaPeakTime());
    double x =
      detProp.ConvertTicksToX(t, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
    double xerr = terr * detProp.GetXTicksCoefficient();
    hitstatev.emplace_back(
      x, hitErr2ScaleFact_ * xerr * xerr, hit->WireID(), geom->WireIDToWireGeo(hit->WireID()));
    //
    if (fsize > 0 && ihit < fsize)
      hitflagsv.push_back(flags[ihit].mask());
    else
      hitflagsv.push_back(recob::TrajectoryPointFlags::DefaultFlagsMask());
    //
    if (rejectHighMultHits_ && hit->Multiplicity() > 1) {
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::Merged);
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    }
    if (rejectHitsNegativeGOF_ && hit->GoodnessOfFit() < 0) {
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::Suspicious);
      hitflagsv.back().set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    }
    //
    if (dumpLevel_ > 2)
      std::cout << "pushed flag mask=" << hitflagsv.back()
                << " merged=" << hitflagsv.back().isSet(recob::TrajectoryPointFlagTraits::Merged)
                << " suspicious="
                << hitflagsv.back().isSet(recob::TrajectoryPointFlagTraits::Suspicious)
                << " nopoint=" << hitflagsv.back().isSet(recob::TrajectoryPointFlagTraits::NoPoint)
                << std::endl;
  }
  if (dumpLevel_ > 2) assert(hits.size() == hitstatev.size());
  return true;
}

bool trkf::TrackKalmanFitter::doFitWork(KFTrackState& trackState,
                                        detinfo::DetectorPropertiesData const& detProp,
                                        std::vector<HitState>& hitstatev,
                                        std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
                                        std::vector<KFTrackState>& fwdPrdTkState,
                                        std::vector<KFTrackState>& fwdUpdTkState,
                                        std::vector<unsigned int>& hitstateidx,
                                        std::vector<unsigned int>& rejectedhsidx,
                                        std::vector<unsigned int>& sortedtksidx,
                                        bool applySkipClean) const
{
  fwdPrdTkState.clear();
  fwdUpdTkState.clear();
  hitstateidx.clear();
  rejectedhsidx.clear();
  sortedtksidx.clear();
  // these three vectors are aligned
  fwdPrdTkState.reserve(hitstatev.size());
  fwdUpdTkState.reserve(hitstatev.size());
  hitstateidx.reserve(hitstatev.size());

  // keep a copy in case first propagation fails
  KFTrackState startState = trackState;

  if (sortHitsByPlane_) {
    //array of hit indices in planes, keeping the original sorting by plane
    const unsigned int nplanes = geom->MaxPlanes();
    std::vector<std::vector<unsigned int>> hitsInPlanes(nplanes);
    for (unsigned int ihit = 0; ihit < hitstatev.size(); ihit++) {
      hitsInPlanes[hitstatev[ihit].wireId().Plane].push_back(ihit);
    }
    if (sortHitsByWire_) {
      for (unsigned int iplane = 0; iplane < nplanes; ++iplane) {
        if (geom->Plane(iplane).GetIncreasingWireDirection<Vector_t>().Dot(trackState.momentum()) >
            0) {
          std::sort(hitsInPlanes[iplane].begin(),
                    hitsInPlanes[iplane].end(),
                    [hitstatev](const unsigned int& a, const unsigned int& b) -> bool {
                      return hitstatev[a].wireId().Wire < hitstatev[b].wireId().Wire;
                    });
        }
        else {
          std::sort(hitsInPlanes[iplane].begin(),
                    hitsInPlanes[iplane].end(),
                    [hitstatev](const unsigned int& a, const unsigned int& b) -> bool {
                      return hitstatev[a].wireId().Wire > hitstatev[b].wireId().Wire;
                    });
        }
      }
    }

    //dump hits sorted in each plane
    if (dumpLevel_ > 1) {
      int ch = 0;
      for (auto p : hitsInPlanes) {
        for (auto h : p) {
          std::cout << "hit #/Plane/Wire/x/mask: " << ch++ << " " << hitstatev[h].wireId().Plane
                    << " " << hitstatev[h].wireId().Wire << " " << hitstatev[h].hitMeas() << " "
                    << hitflagsv[h] << std::endl;
        }
      }
    }

    //array of indices, where iterHitsInPlanes[i] is the iterator over hitsInPlanes[i]
    std::vector<unsigned int> iterHitsInPlanes(nplanes, 0);
    for (unsigned int p = 0; p < hitstatev.size(); ++p) {
      if (dumpLevel_ > 1) std::cout << std::endl << "processing hit #" << p << std::endl;
      if (dumpLevel_ > 1)
        std::cout << "hit sizes: rej=" << rejectedhsidx.size() << " good=" << hitstateidx.size()
                  << " input=" << hitstatev.size() << std::endl;
      if (dumpLevel_ > 1) {
        std::cout << "compute distance from state=" << std::endl;
        trackState.dump();
      }
      int min_plane = -1;
      double min_dist = DBL_MAX;
      //find the closest hit according to the sorting in each plane
      for (unsigned int iplane = 0; iplane < nplanes; ++iplane) {
        //note: ih is a reference, so when 'continue' we increment iterHitsInPlanes[iplane] and the hit is effectively rejected
        if (dumpLevel_ > 1)
          std::cout << "iplane=" << iplane << " nplanes=" << nplanes
                    << " iterHitsInPlanes[iplane]=" << iterHitsInPlanes[iplane]
                    << " hitsInPlanes[iplane].size()=" << hitsInPlanes[iplane].size() << std::endl;
        for (unsigned int& ih = iterHitsInPlanes[iplane]; ih < hitsInPlanes[iplane].size(); ++ih) {
          //
          unsigned int ihit = hitsInPlanes[iplane][ih];
          if (dumpLevel_ > 1)
            std::cout << "ih=" << ih << " ihit=" << ihit << " size=" << hitsInPlanes[iplane].size()
                      << std::endl;
          const auto& hitstate = hitstatev[ihit];
          const auto& hitflags = hitflagsv[ihit];
          if (hitflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint) ||
              hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ||
              hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected)) {
            if (dumpLevel_ > 1)
              std::cout << "rejecting hit flags="
                        << hitflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint) << ", "
                        << hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) << ", "
                        << hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected) << " "
                        << hitflags << std::endl;
            rejectedhsidx.push_back(ihit);
            continue;
          }
          //get distance to measurement surface
          bool success = true;
          const double dist =
            propagator->distanceToPlane(success, trackState.trackState(), hitstate.plane());
          if (dumpLevel_ > 1)
            std::cout << "distance to plane " << iplane << " wire=" << hitstate.wireId().Wire
                      << " = " << dist << ", min_dist=" << std::min(min_dist, 999.)
                      << " min_plane=" << min_plane << " success=" << success
                      << " wirepo=" << hitstate.plane().position() << std::endl;
          if (!success) {
            rejectedhsidx.push_back(ihit);
            continue;
          }
          if (applySkipClean && skipNegProp_ && fwdUpdTkState.size() > 0 &&
              dist < negDistTolerance_) {
            rejectedhsidx.push_back(ihit);
            continue;
          }
          if (dist < min_dist) {
            min_plane = iplane;
            min_dist = dist;
          }
          break;
        }
      }
      if (dumpLevel_ > 1)
        std::cout
          << "pick min_dist=" << min_dist << " min_plane=" << min_plane << " wire="
          << (min_plane < 0 ?
                -1 :
                hitstatev[hitsInPlanes[min_plane][iterHitsInPlanes[min_plane]]].wireId().Wire)
          << std::endl;
      //
      //now we know which is the closest wire: get the hitstate and increment the iterator
      if (min_plane < 0) {
        if (rejectedhsidx.size() + hitstateidx.size() == hitstatev.size())
          break;
        else
          continue;
      }
      unsigned int ihit = hitsInPlanes[min_plane][iterHitsInPlanes[min_plane]];
      const auto* hitstate = &hitstatev[ihit];
      iterHitsInPlanes[min_plane]++;
      //
      //propagate to measurement surface
      bool propok = true;
      trackState = propagator->propagateToPlane(propok,
                                                detProp,
                                                trackState.trackState(),
                                                hitstate->plane(),
                                                true,
                                                true,
                                                TrackStatePropagator::FORWARD);
      if (!propok && !(applySkipClean && fwdUpdTkState.size() > 0 && skipNegProp_)) {
        if (dumpLevel_ > 1) std::cout << "attempt backward prop" << std::endl;
        trackState = propagator->propagateToPlane(propok,
                                                  detProp,
                                                  trackState.trackState(),
                                                  hitstate->plane(),
                                                  true,
                                                  true,
                                                  TrackStatePropagator::BACKWARD);
      }
      if (dumpLevel_ > 1) {
        std::cout << "hit state " << std::endl;
        hitstate->dump();
        std::cout << "propagation result=" << propok << std::endl;
        std::cout << "propagated state " << std::endl;
        trackState.dump();
        std::cout << "propagated planarity="
                  << hitstate->plane().direction().Dot(hitstate->plane().position() -
                                                       trackState.position())
                  << std::endl;
        std::cout << "residual=" << trackState.residual(*hitstate)
                  << " chi2=" << trackState.chi2(*hitstate) << std::endl;
      }
      if (propok) {

        //reject the first hit if its residue is too large (due to e.g. spurious hit or bad hit sorting)
        if (applySkipClean && fwdUpdTkState.size() == 0 &&
            std::abs(trackState.residual(*hitstate)) > maxResidueFirstHit_) {
          if (dumpLevel_ > 1)
            std::cout << "rejecting first hit with residual=" << trackState.residual(*hitstate)
                      << std::endl;
          rejectedhsidx.push_back(ihit);
          trackState = startState;
          continue;
        }

        //if there is >1 consecutive hit on the same wire, choose the one with best chi2 and exclude the others (should this be optional?)
        if (pickBestHitOnWire_) {
          auto wire = hitstate->wireId().Wire;
          float min_chi2_wire = trackState.chi2(*hitstate);
          for (unsigned int jh = iterHitsInPlanes[min_plane]; jh < hitsInPlanes[min_plane].size();
               ++jh) {
            const unsigned int jhit = hitsInPlanes[min_plane][jh];
            const auto& jhitstate = hitstatev[jhit];
            if (jhitstate.wireId().Wire != wire) break;
            if (ihit != jhit) iterHitsInPlanes[min_plane]++;
            float chi2 = trackState.chi2(jhitstate);
            if (dumpLevel_ > 1 && ihit != jhit)
              std::cout << "additional hit on the same wire with jhit=" << jhit << " chi2=" << chi2
                        << " ihit=" << ihit << " min_chi2_wire=" << min_chi2_wire << std::endl;
            if (chi2 < min_chi2_wire) {
              min_chi2_wire = chi2;
              if (dumpLevel_ > 1 && ihit != jhit) {
                std::cout << "\tnow using ";
                jhitstate.dump();
              }
              if (ihit != jhit) rejectedhsidx.push_back(ihit);
              ihit = jhit;
            }
            else {
              rejectedhsidx.push_back(jhit);
            }
          }
        }

        hitstate = &hitstatev[ihit];
        auto& hitflags = hitflagsv[ihit];

        //reject hits failing maxResidue, maxChi2, or maxDist cuts
        if (fwdUpdTkState.size() > 0 && applySkipClean &&
            (std::abs(trackState.residual(*hitstate)) > maxResidue_ ||
             trackState.chi2(*hitstate) > maxChi2_ || min_dist > maxDist_)) {
          //
          if (dumpLevel_ > 1)
            std::cout << "rejecting hit with res=" << std::abs(trackState.residual(*hitstate))
                      << " chi2=" << trackState.chi2(*hitstate) << " dist=" << min_dist
                      << std::endl;
          // reset the current state, do not update the hit, and mark as excluded from the fit
          if (fwdUpdTkState.size() > 0)
            trackState = fwdUpdTkState.back();
          else
            trackState = startState;
          rejectedhsidx.push_back(ihit);
          hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
          hitflags.set(
            recob::TrajectoryPointFlagTraits::
              NoPoint); //fixme: this is only for event display, also needed by current definition of ValidPoint
          continue;
        }

        hitstateidx.push_back(ihit);
        fwdPrdTkState.push_back(trackState);

        //hits with problematic flags are kept but not used for update and flagged as excluded from fit
        if (hitflags.isSet(recob::TrajectoryPointFlagTraits::HitIgnored) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::Merged) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::Shared) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::DeltaRay) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::DetectorIssue) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::Suspicious)) {
          //
          // reset the current state, do not update the hit, and mark as excluded from the fit
          if (fwdUpdTkState.size() > 0)
            trackState = fwdUpdTkState.back();
          else
            trackState = startState;
          fwdUpdTkState.push_back(trackState);
          hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
          hitflags.set(
            recob::TrajectoryPointFlagTraits::
              NoPoint); //fixme: this is only for event display, also needed by current definition of ValidPoint
          continue;
        }
        //now update the forward fitted track
        bool upok = trackState.updateWithHitState(*hitstate);
        if (upok == 0) {
          mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
          return false;
        }
        if (dumpLevel_ > 1) {
          std::cout << "updated state " << std::endl;
          trackState.dump();
        }
        fwdUpdTkState.push_back(trackState);
      }
      else {
        if (dumpLevel_ > 0)
          std::cout << "WARNING: forward propagation failed. Skip this hit..." << std::endl;
        // restore the last successful propagation
        if (fwdPrdTkState.size() > 0)
          trackState = fwdPrdTkState.back();
        else
          trackState = startState;
        rejectedhsidx.push_back(ihit);
        continue;
      }
      if (rejectedhsidx.size() + hitstateidx.size() == hitstatev.size()) break;
    }
  }
  else {
    for (unsigned int ihit = 0; ihit < hitstatev.size(); ++ihit) {
      const auto& hitstate = hitstatev[ihit];
      if (dumpLevel_ > 1) {
        std::cout << std::endl << "processing hit #" << ihit << std::endl;
        hitstate.dump();
      }
      auto& hitflags = hitflagsv[ihit];
      if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ||
          hitflags.isSet(recob::TrajectoryPointFlagTraits::Rejected)) {
        rejectedhsidx.push_back(ihit);
        continue;
      }
      bool success = true;
      const double dist =
        propagator->distanceToPlane(success, trackState.trackState(), hitstate.plane());
      if (applySkipClean && fwdUpdTkState.size() > 0 && skipNegProp_) {
        if (dist < 0. || success == false) {
          if (dumpLevel_ > 0)
            std::cout << "WARNING: negative propagation distance. Skip this hit..." << std::endl;
          ;
          rejectedhsidx.push_back(ihit);
          continue;
        }
      }
      //propagate to measurement surface
      bool propok = true;
      trackState = propagator->propagateToPlane(propok,
                                                detProp,
                                                trackState.trackState(),
                                                hitstate.plane(),
                                                true,
                                                true,
                                                TrackStatePropagator::FORWARD);
      if (!propok && !(applySkipClean && skipNegProp_))
        trackState = propagator->propagateToPlane(propok,
                                                  detProp,
                                                  trackState.trackState(),
                                                  hitstate.plane(),
                                                  true,
                                                  true,
                                                  TrackStatePropagator::BACKWARD);
      if (propok) {
        hitstateidx.push_back(ihit);
        fwdPrdTkState.push_back(trackState);
        //
        if (hitflags.isSet(recob::TrajectoryPointFlagTraits::HitIgnored) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::Merged) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::Shared) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::DeltaRay) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::DetectorIssue) ||
            hitflags.isSet(recob::TrajectoryPointFlagTraits::Suspicious) ||
            (fwdUpdTkState.size() > 0 && applySkipClean &&
             (std::abs(trackState.residual(hitstate)) > maxResidue_ ||
              trackState.chi2(hitstate) > maxChi2_ || dist > maxDist_))) {
          if (dumpLevel_ > 1)
            std::cout << "rejecting hit with res=" << std::abs(trackState.residual(hitstate))
                      << " chi2=" << trackState.chi2(hitstate) << " dist=" << dist << std::endl;
          // reset the current state, do not update the hit, mark as excluded from the fit
          if (fwdUpdTkState.size() > 0)
            trackState = fwdUpdTkState.back();
          else
            trackState = startState;
          fwdUpdTkState.push_back(trackState);
          hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
          hitflags.set(
            recob::TrajectoryPointFlagTraits::
              NoPoint); //fixme: this is only for event display, also needed by current definition of ValidPoint
          continue;
        }
        //
        bool upok = trackState.updateWithHitState(hitstate);
        if (upok == 0) {
          mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
          return false;
        }
        fwdUpdTkState.push_back(trackState);
      }
      else {
        if (dumpLevel_ > 0)
          std::cout << "WARNING: forward propagation failed. Skip this hit..." << std::endl;
        ;
        // restore the last successful propagation
        if (fwdPrdTkState.size() > 0)
          trackState = fwdPrdTkState.back();
        else
          trackState = startState;
        rejectedhsidx.push_back(ihit);
        continue;
      }
    } //for (auto hitstate : hitstatev)
  }

  if (dumpLevel_ > 2) assert(rejectedhsidx.size() + hitstateidx.size() == hitstatev.size());
  if (dumpLevel_ > 0) {
    std::cout << "TRACK AFTER FWD" << std::endl;
    std::cout << "hit sizes=" << rejectedhsidx.size() << " " << hitstateidx.size() << " "
              << hitstatev.size() << std::endl;
    trackState.dump();
  }

  //reinitialize trf for backward fit, scale the error to avoid biasing the backward fit
  trackState.setCovariance(100. * trackState.covariance());

  startState = trackState;

  //backward loop over track states and hits in fwdUpdTracks: use hits for backward fit and fwd track states for smoothing
  int nvalidhits = 0;
  for (int itk = fwdPrdTkState.size() - 1; itk >= 0; itk--) {
    auto& fwdPrdTrackState = fwdPrdTkState[itk];
    auto& fwdUpdTrackState = fwdUpdTkState[itk];
    const auto& hitstate = hitstatev[hitstateidx[itk]];
    auto& hitflags = hitflagsv[hitstateidx[itk]];
    if (dumpLevel_ > 1) {
      std::cout << std::endl << "processing backward hit #" << itk << std::endl;
      hitstate.dump();
    }
    bool propok = true;
    trackState = propagator->propagateToPlane(propok,
                                              detProp,
                                              trackState.trackState(),
                                              hitstate.plane(),
                                              true,
                                              true,
                                              TrackStatePropagator::BACKWARD);
    if (!propok)
      trackState = propagator->propagateToPlane(propok,
                                                detProp,
                                                trackState.trackState(),
                                                hitstate.plane(),
                                                true,
                                                true,
                                                TrackStatePropagator::FORWARD); //do we want this?
    //
    if (dumpLevel_ > 1) {
      std::cout << "propagation result=" << propok << std::endl;
      std::cout << "propagated state " << std::endl;
      trackState.dump();
      std::cout << "propagated planarity="
                << hitstate.plane().direction().Dot(hitstate.plane().position() -
                                                    trackState.position())
                << std::endl;
    }
    if (propok) {
      //combine forward predicted and backward predicted
      bool pcombok = fwdPrdTrackState.combineWithTrackState(trackState.trackState());
      if (pcombok == 0) {
        mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
        return false;
      }
      if (dumpLevel_ > 1) {
        std::cout << "combined prd state " << std::endl;
        fwdPrdTrackState.dump();
      }
      // combine forward updated and backward predicted and update backward predicted, only if the hit was not excluded
      if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) == 0) {
        bool ucombok = fwdUpdTrackState.combineWithTrackState(trackState.trackState());
        if (ucombok == 0) {
          mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
          return false;
        }
        if (dumpLevel_ > 1) {
          std::cout << "combined upd state " << std::endl;
          fwdUpdTrackState.dump();
        }
        bool upok = trackState.updateWithHitState(hitstate);
        if (upok == 0) {
          mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
          return false;
        }
        if (dumpLevel_ > 1) {
          std::cout << "updated state " << std::endl;
          trackState.dump();
        }
        // Keep a copy in case a future propagation fails
        startState = trackState;
        //
        nvalidhits++;
      }
      else {
        fwdUpdTrackState = fwdPrdTrackState;
      }
    }
    else {
      // ok, if the backward propagation failed we exclude this point from the rest of the fit,
      // but we can still use its position from the forward fit, so just mark it as ExcludedFromFit
      hitflags.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
      hitflags.set(
        recob::TrajectoryPointFlagTraits::NoPoint); // fixme: this is only for event display, also
                                                    // needed by current definition of ValidPoint
      if (dumpLevel_ > 0)
        std::cout << "WARNING: backward propagation failed. Skip this hit... hitstateidx[itk]="
                  << hitstateidx[itk] << " itk=" << itk << std::endl;
      ;
      // restore the last successful propagation
      trackState = startState;
      continue;
    }
  } //for (unsigned int itk = fwdPrdTkState.size()-1; itk>=0; itk--) {

  if (nvalidhits < 4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  if (dumpLevel_ > 1) std::cout << "sort output with nvalidhits=" << nvalidhits << std::endl;

  // sort output states
  sortOutput(
    hitstatev, fwdUpdTkState, hitstateidx, rejectedhsidx, sortedtksidx, hitflagsv, applySkipClean);
  size_t nsortvalid = 0;
  for (auto& idx : sortedtksidx) {
    auto& hitflags = hitflagsv[hitstateidx[idx]];
    if (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) == 0) nsortvalid++;
  }
  if (nsortvalid < 4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  if (dumpLevel_ > 2) assert(rejectedhsidx.size() + sortedtksidx.size() == hitstatev.size());
  return true;
}

void trkf::TrackKalmanFitter::sortOutput(
  std::vector<HitState>& hitstatev,
  std::vector<KFTrackState>& fwdUpdTkState,
  std::vector<unsigned int>& hitstateidx,
  std::vector<unsigned int>& rejectedhsidx,
  std::vector<unsigned int>& sortedtksidx,
  std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
  bool applySkipClean) const
{
  //
  if (sortOutputHitsMinLength_) {
    //sort hits keeping fixed the order on planes and picking the closest next plane
    const unsigned int nplanes = geom->MaxPlanes();
    std::vector<std::vector<unsigned int>> tracksInPlanes(nplanes);
    for (unsigned int p = 0; p < hitstateidx.size(); ++p) {
      const auto& hitstate = hitstatev[hitstateidx[p]];
      tracksInPlanes[hitstate.wireId().Plane].push_back(p);
    }
    if (dumpLevel_ > 2) {
      for (auto s : fwdUpdTkState) {
        std::cout << "state pos=" << s.position() << std::endl;
      }
    }
    //find good starting point
    std::vector<unsigned int> iterTracksInPlanes;
    for (auto it : tracksInPlanes)
      iterTracksInPlanes.push_back(0);
    auto pos = fwdUpdTkState.front().position();
    auto dir = fwdUpdTkState.front().momentum();
    unsigned int p = 0;
    for (; p < fwdUpdTkState.size(); ++p) {
      if (hitflagsv[hitstateidx[p]].isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) == 0) {
        pos = fwdUpdTkState[p].position();
        dir = fwdUpdTkState[p].momentum();
        if (dumpLevel_ > 2)
          std::cout << "sort output found point not excluded with p=" << p
                    << " hitstateidx[p]=" << hitstateidx[p] << " pos=" << pos << std::endl;
        break;
      }
      else {
        rejectedhsidx.push_back(hitstateidx[p]);
      }
    }
    if (dumpLevel_ > 1)
      std::cout << "sort output init with pos=" << pos << " dir=" << dir << std::endl;
    //pick hits based on minimum dot product with respect to position and direction at previous hit
    for (; p < fwdUpdTkState.size(); ++p) {
      int min_plane = -1;
      double min_dotp = DBL_MAX;
      for (unsigned int iplane = 0; iplane < iterTracksInPlanes.size(); ++iplane) {
        for (unsigned int& itk = iterTracksInPlanes[iplane]; itk < tracksInPlanes[iplane].size();
             ++itk) {
          auto& trackstate = fwdUpdTkState[tracksInPlanes[iplane][iterTracksInPlanes[iplane]]];
          auto& tmppos = trackstate.position();
          const double dotp = dir.Dot(tmppos - pos);
          if (dumpLevel_ > 2)
            std::cout << "iplane=" << iplane << " tmppos=" << tmppos << " tmpdir=" << tmppos - pos
                      << " dotp=" << dotp << std::endl;
          if (dotp < min_dotp) {
            min_plane = iplane;
            min_dotp = dotp;
          }
          break;
        }
      }
      if (min_plane < 0) continue;
      const unsigned int ihit = tracksInPlanes[min_plane][iterTracksInPlanes[min_plane]];
      if (applySkipClean && skipNegProp_ && min_dotp < negDistTolerance_ &&
          sortedtksidx.size() > 0) {
        if (dumpLevel_ > 2)
          std::cout << "sort output rejecting hit #" << ihit << " plane=" << min_plane
                    << " with min_dotp=" << min_dotp << std::endl;
        rejectedhsidx.push_back(hitstateidx[ihit]);
        iterTracksInPlanes[min_plane]++;
        continue;
      }
      if (dumpLevel_ > 2)
        std::cout << "sort output picking hit #" << ihit << " plane=" << min_plane
                  << " with min_dotp=" << min_dotp << std::endl;
      auto& trackstate = fwdUpdTkState[ihit];
      pos = trackstate.position();
      dir = trackstate.momentum();
      //
      sortedtksidx.push_back(ihit);
      iterTracksInPlanes[min_plane]++;
    }
  }
  else {
    for (unsigned int p = 0; p < fwdUpdTkState.size(); ++p) {
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
      unsigned int i = 1;
      unsigned int end = sortedtksidx.size() - 1;
      for (; i < end; ++i) {
        auto dir0 = fwdUpdTkState[sortedtksidx[i]].position() - pos0;
        auto dir2 =
          fwdUpdTkState[sortedtksidx[i + 1]].position() - fwdUpdTkState[sortedtksidx[i]].position();
        dir0 /= dir0.R();
        dir2 /= dir2.R();
        if (dir2.Dot(dir0) < 0.) {
          broken = true;
          end--;
          break;
        }
        else
          pos0 = fwdUpdTkState[sortedtksidx[i]].position();
      }
      if (!broken) { clean = true; }
      else {
        rejectedhsidx.push_back(hitstateidx[sortedtksidx[i]]);
        sortedtksidx.erase(sortedtksidx.begin() + i);
      }
    }
  }
  //
}

bool trkf::TrackKalmanFitter::fillResult(
  const std::vector<art::Ptr<recob::Hit>>& inHits,
  const int tkID,
  const int pdgid,
  std::vector<HitState>& hitstatev,
  std::vector<recob::TrajectoryPointFlags::Mask_t>& hitflagsv,
  std::vector<KFTrackState>& fwdPrdTkState,
  std::vector<KFTrackState>& fwdUpdTkState,
  std::vector<unsigned int>& hitstateidx,
  std::vector<unsigned int>& rejectedhsidx,
  std::vector<unsigned int>& sortedtksidx,
  recob::Track& outTrack,
  std::vector<art::Ptr<recob::Hit>>& outHits,
  trkmkr::OptionalOutputs& optionals) const
{
  // fill output trajectory objects with smoothed track and its hits
  int nvalidhits = 0;
  trkmkr::TrackCreationBookKeeper tcbk(outHits, optionals, tkID, pdgid, true);
  for (unsigned int p : sortedtksidx) {
    const auto& trackstate = fwdUpdTkState[p];
    const auto& hitflags = hitflagsv[hitstateidx[p]];
    const unsigned int originalPos = hitstateidx[p];
    if (dumpLevel_ > 2) assert(originalPos >= 0 && originalPos < hitstatev.size());
    //
    const auto& prdtrack = fwdPrdTkState[p];
    const auto& hitstate = hitstatev[hitstateidx[p]];
    if (dumpLevel_ > 2) assert(hitstate.wireId().Plane == inHits[originalPos]->WireID().Plane);
    //
    float chi2 =
      (hitflags.isSet(recob::TrajectoryPointFlagTraits::ExcludedFromFit) ? -1. :
                                                                           prdtrack.chi2(hitstate));
    if (hitflags.isSet(recob::TrajectoryPointFlagTraits::NoPoint) == false) nvalidhits++;
    //
    trkmkr::OptionalPointElement ope;
    if (optionals.isTrackFitInfosInit()) {
      ope.setTrackFitHitInfo(recob::TrackFitHitInfo(hitstate.hitMeas(),
                                                    hitstate.hitMeasErr2(),
                                                    prdtrack.parameters(),
                                                    prdtrack.covariance(),
                                                    hitstate.wireId()));
    }
    tcbk.addPoint(trackstate.position(),
                  trackstate.momentum(),
                  inHits[originalPos],
                  recob::TrajectoryPointFlags(originalPos, hitflags),
                  chi2,
                  ope);
  }
  if (dumpLevel_ > 0) std::cout << "fillResult nvalidhits=" << nvalidhits << std::endl;

  // fill also with rejected hits information
  SMatrixSym55 fakeCov55;
  for (int i = 0; i < 5; i++)
    for (int j = i; j < 5; j++)
      fakeCov55(i, j) = util::kBogusD;
  for (unsigned int rejidx = 0; rejidx < rejectedhsidx.size(); ++rejidx) {
    const unsigned int originalPos = rejectedhsidx[rejidx];
    auto& mask = hitflagsv[rejectedhsidx[rejidx]];
    mask.set(recob::TrajectoryPointFlagTraits::HitIgnored,
             recob::TrajectoryPointFlagTraits::NoPoint);
    if (mask.isSet(recob::TrajectoryPointFlagTraits::Rejected) == 0)
      mask.set(recob::TrajectoryPointFlagTraits::ExcludedFromFit);
    //
    const auto& hitstate = hitstatev[rejectedhsidx[rejidx]];
    if (dumpLevel_ > 2) assert(hitstate.wireId().Plane == inHits[originalPos]->WireID().Plane);
    trkmkr::OptionalPointElement ope;
    if (optionals.isTrackFitInfosInit()) {
      ope.setTrackFitHitInfo(recob::TrackFitHitInfo(
        hitstate.hitMeas(),
        hitstate.hitMeasErr2(),
        SVector5(util::kBogusD, util::kBogusD, util::kBogusD, util::kBogusD, util::kBogusD),
        fakeCov55,
        hitstate.wireId()));
    }
    tcbk.addPoint(Point_t(util::kBogusD, util::kBogusD, util::kBogusD),
                  Vector_t(util::kBogusD, util::kBogusD, util::kBogusD),
                  inHits[originalPos],
                  recob::TrajectoryPointFlags(originalPos, mask),
                  -1.,
                  ope);
  }

  if (dumpLevel_ > 0)
    std::cout << "outHits.size()=" << outHits.size() << " inHits.size()=" << inHits.size()
              << std::endl;
  if (dumpLevel_ > 2) assert(outHits.size() == inHits.size());

  bool propok = true;
  KFTrackState resultF =
    propagator->rotateToPlane(propok,
                              fwdUpdTkState[sortedtksidx.front()].trackState(),
                              Plane(fwdUpdTkState[sortedtksidx.front()].position(),
                                    fwdUpdTkState[sortedtksidx.front()].momentum()));
  KFTrackState resultB =
    propagator->rotateToPlane(propok,
                              fwdUpdTkState[sortedtksidx.back()].trackState(),
                              Plane(fwdUpdTkState[sortedtksidx.back()].position(),
                                    fwdUpdTkState[sortedtksidx.back()].momentum()));

  outTrack =
    tcbk.finalizeTrack(SMatrixSym55(resultF.covariance()), SMatrixSym55(resultB.covariance()));

  //if there are points with zero momentum return false
  size_t point = 0;
  while (outTrack.HasValidPoint(point)) {
    if (outTrack.MomentumAtPoint(outTrack.NextValidPoint(point++)) <= 1.0e-9) {
      mf::LogWarning("TrackKalmanFitter") << "found point with zero momentum!" << std::endl;
      return false;
    }
  }

  if (dumpLevel_ > 0) {
    std::cout << "outTrack vertex=" << outTrack.Start() << "\ndir=" << outTrack.StartDirection()
              << "\ncov=\n"
              << outTrack.StartCovariance() << "\nlength=" << outTrack.Length()
              << "\nchi2/ndof=" << outTrack.Chi2PerNdof() << std::endl;
  }

  return true;
}
