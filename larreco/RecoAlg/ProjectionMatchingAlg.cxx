////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()
#include "lardata/ArtDataHelper/ToElement.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larreco/RecoAlg/ImagePatternAlgs/DataProvider/DataProviderAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaSegment3D.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1F.h"

#include "range/v3/view.hpp"

using geo::vect::toPoint;
using lar::to_element;
using namespace ranges;

namespace {
  constexpr double step{0.3};
}

pma::ProjectionMatchingAlg::ProjectionMatchingAlg(const pma::ProjectionMatchingAlg::Config& config)
  : fOptimizationEps{config.OptimizationEps()}
  , fFineTuningEps{config.FineTuningEps()}
  , fTrkValidationDist2D{config.TrkValidationDist2D()}
  , fHitTestingDist2D{config.HitTestingDist2D()}
  , fMinTwoViewFraction{config.MinTwoViewFraction()}
  , fGeom{lar::providerFrom<geo::Geometry>()}
  , fChannelMapAlg{art::ServiceHandle<geo::ExptGeoHelperInterface const>()->ChannelMapAlgPtr()}
{
  pma::Node3D::SetMargin(config.NodeMargin3D());

  pma::Element3D::SetOptFactor(geo::kU, config.HitWeightU());
  pma::Element3D::SetOptFactor(geo::kV, config.HitWeightV());
  pma::Element3D::SetOptFactor(geo::kZ, config.HitWeightZ());
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate_on_adc(
  const detinfo::DetectorPropertiesData& detProp,
  const lariov::ChannelStatusProvider& channelStatus,
  const pma::Track3D& trk,
  const img::DataProviderAlg& adcImage,
  float const thr) const
{
  unsigned int nAll = 0, nPassed = 0;
  unsigned int testPlane = adcImage.Plane();

  std::vector<unsigned int> trkTPCs = trk.TPCs();
  std::vector<unsigned int> trkCryos = trk.Cryos();

  // check how pixels with a high signal are distributed along the track
  // namely: are there track sections crossing empty spaces, except dead wires?
  pma::Vector3D p(
    trk.front()->Point3D().X(), trk.front()->Point3D().Y(), trk.front()->Point3D().Z());
  for (auto const* seg : trk.Segments()) {
    if (seg->TPC() < 0) // skip segments between tpc's, look only at those contained in tpc
    {
      p = seg->End();
      continue;
    }
    pma::Vector3D p0 = seg->Start();
    pma::Vector3D p1 = seg->End();

    pma::Node3D* node = static_cast<pma::Node3D*>(seg->Prev());

    unsigned tpc = seg->TPC();
    unsigned cryo = seg->Cryo();

    pma::Vector3D dc = step * seg->GetDirection3D();

    double f = pma::GetSegmentProjVector(p, p0, p1);
    while ((f < 1.0) && node->SameTPC(p)) {
      geo::PlaneID const planeID{cryo, tpc, testPlane};
      pma::Vector2D const p2d(fGeom->Plane(planeID).WireCoordinate(toPoint(p)), p.X());
      geo::WireID const wireID(planeID, (int)p2d.X());

      int widx = (int)p2d.X();
      int didx = (int)detProp.ConvertXToTicks(p2d.Y(), testPlane, tpc, cryo);

      if (fGeom->HasWire(wireID)) {
        raw::ChannelID_t ch = fChannelMapAlg->PlaneWireToChannel(wireID);
        if (channelStatus.IsGood(ch)) {
          float max_adc = adcImage.poolMax(widx, didx, 2); // +/- 2 wires, can be parameterized
          if (max_adc > thr) nPassed++;

          nAll++;
        }
      }

      p += dc;
      f = pma::GetSegmentProjVector(p, p0, p1);
    }

    p = seg->End(); // need to have it at the end due to the p in the first iter
                    // set to the hit position, not segment start
  }

  if (nAll > 0) {
    double v = nPassed / (double)nAll;
    mf::LogVerbatim("ProjectionMatchingAlg") << "  trk fraction ok: " << v;
    return v;
  }

  return 1.0;
}

namespace {
  struct hits_on_plane {
    bool operator()(recob::Hit const& hit) const { return hit.WireID().Plane == plane_; }
    unsigned int const plane_;
  };
}

// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate_on_adc_test(
  const detinfo::DetectorPropertiesData& detProp,
  const lariov::ChannelStatusProvider& channelStatus,
  const pma::Track3D& trk,
  const img::DataProviderAlg& adcImage,
  const std::vector<art::Ptr<recob::Hit>>& hits,
  TH1F* histoPassing,
  TH1F* histoRejected) const
{
  double max_d = fTrkValidationDist2D;
  double const max_d2 = max_d * max_d;
  unsigned int nAll = 0, nPassed = 0;
  unsigned int testPlane = adcImage.Plane();

  std::vector<unsigned int> trkTPCs = trk.TPCs();
  std::vector<unsigned int> trkCryos = trk.Cryos();
  std::map<std::pair<unsigned int, unsigned int>, std::pair<TVector2, TVector2>> ranges;
  std::map<std::pair<unsigned int, unsigned int>, double> wirePitch;
  for (auto const& [t, c] : views::cartesian_product(trkTPCs, trkCryos)) {
    ranges[{t, c}] = trk.WireDriftRange(detProp, testPlane, t, c);
    wirePitch[{t, c}] = fGeom->TPC(geo::TPCID(c, t)).Plane(testPlane).WirePitch();
  }

  unsigned int tpc, cryo;
  std::map<std::pair<unsigned int, unsigned int>, std::vector<pma::Vector2D>> all_close_points;

  for (auto const& h :
       hits | views::transform(to_element) | views::filter(hits_on_plane{testPlane})) {
    tpc = h.WireID().TPC;
    cryo = h.WireID().Cryostat;
    std::pair<unsigned int, unsigned int> tpc_cryo(tpc, cryo);
    std::pair<TVector2, TVector2> rect = ranges[tpc_cryo];

    if ((h.WireID().Wire > rect.first.X() - 10) &&  // check only hits in the rectangle around
        (h.WireID().Wire < rect.second.X() + 10) && // the track projection, it is faster than
        (h.PeakTime() > rect.first.Y() - 100) &&    // calculation of trk.Dist2(p2d, testPlane)
        (h.PeakTime() < rect.second.Y() + 100)) {
      TVector2 p2d(wirePitch[tpc_cryo] * h.WireID().Wire,
                   detProp.ConvertTicksToX(h.PeakTime(), testPlane, tpc, cryo));

      double const d2 = trk.Dist2(p2d, testPlane, tpc, cryo);

      if (d2 < max_d2) { all_close_points[tpc_cryo].emplace_back(p2d.X(), p2d.Y()); }
    }
  }

  // then check how points-close-to-the-track-projection are distributed along
  // the track, namely: are there track sections crossing empty spaces, except
  // dead wires?
  pma::Vector3D p(
    trk.front()->Point3D().X(), trk.front()->Point3D().Y(), trk.front()->Point3D().Z());
  for (auto const* seg : trk.Segments()) {
    if (seg->TPC() < 0) // skip segments between tpc's, look only at those contained in tpc
    {
      p = seg->End();
      continue;
    }
    pma::Vector3D p0 = seg->Start();
    pma::Vector3D p1 = seg->End();

    pma::Node3D* node = static_cast<pma::Node3D*>(seg->Prev());

    tpc = seg->TPC();
    cryo = seg->Cryo();

    pma::Vector3D dc = step * seg->GetDirection3D();

    auto const& points = all_close_points[{tpc, cryo}];

    double f = pma::GetSegmentProjVector(p, p0, p1);

    geo::PlaneID const planeID{cryo, tpc, testPlane};
    double wirepitch = fGeom->Plane(planeID).WirePitch();
    while ((f < 1.0) && node->SameTPC(p)) {
      pma::Vector2D p2d(fGeom->Plane(planeID).WireCoordinate(toPoint(p)), p.X());
      geo::WireID const wireID{planeID, static_cast<unsigned int>(p2d.X())};

      int widx = (int)p2d.X();
      int didx = (int)detProp.ConvertXToTicks(p2d.Y(), planeID);

      if (fGeom->HasWire(wireID)) {
        raw::ChannelID_t ch = fChannelMapAlg->PlaneWireToChannel(wireID);
        if (channelStatus.IsGood(ch)) {
          bool is_close = false;
          float max_adc = adcImage.poolMax(widx, didx, 2);

          if (points.size()) {
            p2d.SetX(wirepitch * p2d.X());
            for (const auto& h : points) {
              if (pma::Dist2(p2d, h) < max_d2) {
                is_close = true;
                nPassed++;
                break;
              }
            }
          }
          nAll++;

          // now fill the calibration histograms
          if (is_close) {
            if (histoPassing) histoPassing->Fill(max_adc);
          }
          else {
            if (histoRejected) histoRejected->Fill(max_adc);
          }
        }
      }

      p += dc;
      f = pma::GetSegmentProjVector(p, p0, p1);
    }
    p = seg->End();
  }

  if (nAll > 0) {
    double v = nPassed / (double)nAll;
    mf::LogVerbatim("ProjectionMatchingAlg") << "  trk fraction ok: " << v;
    return v;
  }

  return 1.0;
}

// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(const detinfo::DetectorPropertiesData& detProp,
                                            const lariov::ChannelStatusProvider& channelStatus,
                                            const pma::Track3D& trk,
                                            const std::vector<art::Ptr<recob::Hit>>& hits) const
{
  if (hits.empty()) { return 0; }

  double max_d = fTrkValidationDist2D;
  double const max_d2 = max_d * max_d;
  unsigned int nAll = 0, nPassed = 0;
  unsigned int testPlane = hits.front()->WireID().Plane;

  std::vector<unsigned int> trkTPCs = trk.TPCs();
  std::vector<unsigned int> trkCryos = trk.Cryos();
  std::map<std::pair<unsigned int, unsigned int>, std::pair<TVector2, TVector2>> ranges;
  std::map<std::pair<unsigned int, unsigned int>, double> wirePitch;
  for (auto const& [t, c] : views::cartesian_product(trkTPCs, trkCryos)) {
    ranges[{t, c}] = trk.WireDriftRange(detProp, testPlane, t, c);
    wirePitch[{t, c}] = fGeom->TPC(geo::TPCID(c, t)).Plane(testPlane).WirePitch();
  }

  unsigned int tpc, cryo;
  std::map<std::pair<unsigned int, unsigned int>, std::vector<pma::Vector2D>> all_close_points;

  for (auto const& h :
       hits | views::transform(to_element) | views::filter(hits_on_plane{testPlane})) {
    tpc = h.WireID().TPC;
    cryo = h.WireID().Cryostat;
    std::pair<unsigned int, unsigned int> tpc_cryo(tpc, cryo);
    std::pair<TVector2, TVector2> rect = ranges[tpc_cryo];

    if ((h.WireID().Wire > rect.first.X() - 10) &&  // chceck only hits in the rectangle around
        (h.WireID().Wire < rect.second.X() + 10) && // the track projection, it is faster than
        (h.PeakTime() > rect.first.Y() - 100) &&    // calculation of trk.Dist2(p2d, testPlane)
        (h.PeakTime() < rect.second.Y() + 100)) {
      TVector2 p2d(wirePitch[tpc_cryo] * h.WireID().Wire,
                   detProp.ConvertTicksToX(h.PeakTime(), testPlane, tpc, cryo));

      double const d2 = trk.Dist2(p2d, testPlane, tpc, cryo);

      if (d2 < max_d2) all_close_points[tpc_cryo].emplace_back(p2d.X(), p2d.Y());
    }
  }

  // then check how points-close-to-the-track-projection are distributed along
  // the track, namely: are there track sections crossing empty spaces, except
  // dead wires?
  pma::Vector3D p(
    trk.front()->Point3D().X(), trk.front()->Point3D().Y(), trk.front()->Point3D().Z());
  for (auto const* seg : trk.Segments()) {
    if (seg->TPC() < 0) // skip segments between tpc's, look only at those contained in tpc
    {
      p = seg->End();
      continue;
    }
    pma::Vector3D p0 = seg->Start();
    pma::Vector3D p1 = seg->End();

    pma::Node3D* node = static_cast<pma::Node3D*>(seg->Prev());

    tpc = seg->TPC();
    cryo = seg->Cryo();

    pma::Vector3D dc = step * seg->GetDirection3D();

    auto const& points = all_close_points[{tpc, cryo}];

    double f = pma::GetSegmentProjVector(p, p0, p1);

    geo::PlaneID const planeID{cryo, tpc, testPlane};
    auto const& plane = fGeom->Plane(planeID);
    double const wirepitch = plane.WirePitch();
    while ((f < 1.0) && node->SameTPC(p)) {
      pma::Vector2D p2d(plane.WireCoordinate(toPoint(p)), p.X());
      geo::WireID const wireID{planeID, static_cast<unsigned int>(p2d.X())};
      if (fGeom->HasWire(wireID)) {
        raw::ChannelID_t ch = fChannelMapAlg->PlaneWireToChannel(wireID);
        if (channelStatus.IsGood(ch)) {
          if (points.size()) {
            p2d.SetX(wirepitch * p2d.X());
            for (const auto& h : points) {
              if (pma::Dist2(p2d, h) < max_d2) {
                nPassed++;
                break;
              }
            }
          }
          nAll++;
        }
      }

      p += dc;
      f = pma::GetSegmentProjVector(p, p0, p1);
    }
    p = seg->End();
  }

  if (nAll > 0) {
    double v = nPassed / (double)nAll;
    mf::LogVerbatim("ProjectionMatchingAlg") << "  trk fraction ok: " << v;
    return v;
  }

  return 1.0;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(detinfo::DetectorPropertiesData const& detProp,
                                            const lariov::ChannelStatusProvider& channelStatus,
                                            const TVector3& p0,
                                            const TVector3& p1,
                                            const std::vector<art::Ptr<recob::Hit>>& hits,
                                            unsigned int testPlane,
                                            unsigned int tpc,
                                            unsigned int cryo) const
{
  double max_d = fTrkValidationDist2D;
  double d2, max_d2 = max_d * max_d;
  unsigned int nAll = 0, nPassed = 0;

  TVector3 p(p0);
  TVector3 dc(p1);
  dc -= p;
  dc *= step / dc.Mag();

  double f = pma::GetSegmentProjVector(p, p0, p1);
  geo::PlaneID const planeID{cryo, tpc, testPlane};
  double const wirepitch = fGeom->Plane(planeID).WirePitch();
  while (f < 1.0) {
    TVector2 p2d(fGeom->Plane(planeID).WireCoordinate(toPoint(p)), p.X());
    geo::WireID wireID(planeID, (int)p2d.X());
    if (fGeom->HasWire(wireID)) {
      raw::ChannelID_t ch = fChannelMapAlg->PlaneWireToChannel(wireID);
      if (channelStatus.IsGood(ch)) {
        p2d.Set(wirepitch * p2d.X(), p2d.Y());
        for (const auto& h : hits)
          if ((h->WireID().Plane == testPlane) && (h->WireID().TPC == tpc) &&
              (h->WireID().Cryostat == cryo)) {
            d2 = pma::Dist2(
              p2d,
              pma::WireDriftToCm(detProp, h->WireID().Wire, h->PeakTime(), testPlane, tpc, cryo));
            if (d2 < max_d2) {
              nPassed++;
              break;
            }
          }
        nAll++;
      }
    }

    p += dc;
    f = pma::GetSegmentProjVector(p, p0, p1);
  }

  if (nAll > 3) // validate actually only if 2D projection in testPlane has some minimum length
  {
    double v = nPassed / (double)nAll;
    mf::LogVerbatim("ProjectionMatchingAlg") << "  segment fraction ok: " << v;
    return v;
  }
  else
    return 1.0;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::twoViewFraction(pma::Track3D& trk) const
{
  trk.SelectHits();
  trk.DisableSingleViewEnds();

  size_t idx = 0;
  while ((idx < trk.size() - 1) && !trk[idx]->IsEnabled())
    idx++;
  double l0 = trk.Length(0, idx + 1);

  idx = trk.size() - 1;
  while ((idx > 1) && !trk[idx]->IsEnabled())
    idx--;
  double l1 = trk.Length(idx - 1, trk.size() - 1);

  trk.SelectHits();

  return 1.0 - (l0 + l1) / trk.Length();
}
// ------------------------------------------------------

size_t pma::ProjectionMatchingAlg::getSegCount_(size_t const trk_size)
{
  int const nSegments = 0.8 * trk_size / sqrt(trk_size);
  return std::max(size_t{1}, static_cast<size_t>(nSegments));
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildTrack(
  const detinfo::DetectorPropertiesData& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits_1,
  const std::vector<art::Ptr<recob::Hit>>& hits_2) const
{
  pma::Track3D* trk = new pma::Track3D(); // track candidate
  trk->AddHits(detProp, hits_1);
  trk->AddHits(detProp, hits_2);

  mf::LogVerbatim("ProjectionMatchingAlg") << "track size: " << trk->size();
  std::vector<unsigned int> tpcs = trk->TPCs();
  for (size_t t = 0; t < tpcs.size(); ++t) {
    mf::LogVerbatim("ProjectionMatchingAlg") << "  tpc:" << tpcs[t];
  }
  mf::LogVerbatim("ProjectionMatchingAlg")
    << "  #coll:" << trk->NHits(geo::kZ) << "  #ind2:" << trk->NHits(geo::kV)
    << "  #ind1:" << trk->NHits(geo::kU);

  size_t nSegments = getSegCount_(trk->size());
  size_t nNodes = (size_t)(nSegments - 1); // n nodes to add

  mf::LogVerbatim("ProjectionMatchingAlg") << "  initialize trk";
  if (!trk->Initialize(detProp)) {
    mf::LogWarning("ProjectionMatchingAlg") << "  initialization failed, delete trk";
    delete trk;
    return 0;
  }

  double f = twoViewFraction(*trk);
  if (f > fMinTwoViewFraction) {
    double g = 0.0;
    mf::LogVerbatim("ProjectionMatchingAlg") << "  optimize trk (" << nSegments << " seg)";
    if (nNodes) {
      g = trk->Optimize(detProp, nNodes, fOptimizationEps, false, true, 25,
                        10); // build nodes
      mf::LogVerbatim("ProjectionMatchingAlg") << "  nodes done, g = " << g;
      trk->SelectAllHits();
    }
    g = trk->Optimize(detProp, 0, fFineTuningEps); // final tuning
    mf::LogVerbatim("ProjectionMatchingAlg") << "  tune done, g = " << g;

    trk->SortHits();
    return trk;
  }
  else {
    mf::LogVerbatim("ProjectionMatchingAlg") << "  clusters do not match, f = " << f;
    delete trk;
    return 0;
  }
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildMultiTPCTrack(
  const detinfo::DetectorPropertiesData& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits) const
{
  std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>> hits_by_tpc;
  for (auto const& h : hits) {
    hits_by_tpc[h->WireID().TPC].push_back(h);
  }

  std::vector<pma::Track3D*> tracks;
  for (auto const& hsel : hits_by_tpc) {
    pma::Track3D* trk = buildSegment(detProp, hsel.second);
    if (trk) tracks.push_back(trk);
  }

  bool need_reopt = false;
  while (tracks.size() > 1) {
    need_reopt = true;

    pma::Track3D* first = tracks.front();
    pma::Track3D* best = 0;
    double d, dmin = 1.0e12;
    size_t t_best = 0, cfg = 0;
    for (size_t t = 1; t < tracks.size(); ++t) {
      pma::Track3D* second = tracks[t];

      d = pma::Dist2(first->front()->Point3D(), second->front()->Point3D());
      if (d < dmin) {
        dmin = d;
        best = second;
        t_best = t;
        cfg = 0;
      }

      d = pma::Dist2(first->front()->Point3D(), second->back()->Point3D());
      if (d < dmin) {
        dmin = d;
        best = second;
        t_best = t;
        cfg = 1;
      }

      d = pma::Dist2(first->back()->Point3D(), second->front()->Point3D());
      if (d < dmin) {
        dmin = d;
        best = second;
        t_best = t;
        cfg = 2;
      }

      d = pma::Dist2(first->back()->Point3D(), second->back()->Point3D());
      if (d < dmin) {
        dmin = d;
        best = second;
        t_best = t;
        cfg = 3;
      }
    }
    if (best) {
      switch (cfg) {
      default:
      case 0:
      case 1:
        mergeTracks(detProp, *best, *first, false);
        tracks[0] = best;
        delete first;
        break;

      case 2:
      case 3:
        mergeTracks(detProp, *first, *best, false);
        delete best;
        break;
      }
      tracks.erase(tracks.begin() + t_best);
    }
    else
      break; // should not happen
  }

  pma::Track3D* trk = 0;
  if (!tracks.empty()) {
    trk = tracks.front();
    if (need_reopt) {
      double g = trk->Optimize(detProp, 0, fOptimizationEps);
      mf::LogVerbatim("ProjectionMatchingAlg")
        << "  reopt after merging initial tpc segments: done, g = " << g;
    }

    int nSegments = getSegCount_(trk->size()) - trk->Segments().size() + 1;
    if (nSegments > 0) // need to add segments
    {
      double g = 0.0;
      size_t nNodes = (size_t)(nSegments - 1); // n nodes to add
      if (nNodes) {
        mf::LogVerbatim("ProjectionMatchingAlg") << "  optimize trk (add " << nSegments << " seg)";

        g = trk->Optimize(detProp, nNodes, fOptimizationEps, false, true, 25,
                          10); // build nodes
        mf::LogVerbatim("ProjectionMatchingAlg") << "  nodes done, g = " << g;
        trk->SelectAllHits();
      }
      g = trk->Optimize(detProp, 0, fFineTuningEps); // final tuning
      mf::LogVerbatim("ProjectionMatchingAlg") << "  tune done, g = " << g;
    }
    trk->SortHits();
  }
  return trk;
}

// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildShowerSeg(
  const detinfo::DetectorPropertiesData& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const pma::Vector3D& vtx) const
{
  geo::Point_t const vtxpoint{vtx.X(), vtx.Y(), vtx.Z()};

  if (!fGeom->HasTPC(fGeom->FindTPCAtPosition(vtxpoint))) return 0;

  TVector3 vtxv3(vtx.X(), vtx.Y(), vtx.Z());

  const size_t tpc = fGeom->FindTPCAtPosition(vtxpoint).TPC;
  const size_t cryo = fGeom->PositionToCryostatID(vtxpoint).Cryostat;
  const geo::TPCGeo& tpcgeom = fGeom->TPC(geo::TPCID(cryo, tpc));

  // use only hits from tpc where the vtx is
  std::vector<art::Ptr<recob::Hit>> hitstpc;
  for (size_t h = 0; h < hits.size(); ++h)
    if (hits[h]->WireID().TPC == tpc) hitstpc.push_back(hits[h]);

  if (!hitstpc.size()) return 0;

  std::vector<art::Ptr<recob::Hit>> hitstrk;
  size_t view = 0;
  size_t countviews = 0;
  while (view < 3) {
    mf::LogVerbatim("ProjectionMatchinAlg") << "collecting hits from view: " << view;
    if (!tpcgeom.HasPlane(view)) {
      ++view;
      continue;
    }

    // select hits only for a single view
    std::vector<art::Ptr<recob::Hit>> hitsview;
    for (size_t h = 0; h < hitstpc.size(); ++h)
      if (hitstpc[h]->WireID().Plane == view) hitsview.push_back(hitstpc[h]);
    if (!hitsview.size()) {
      ++view;
      continue;
    }

    // filter our small groups of hits, far from main shower
    std::vector<art::Ptr<recob::Hit>> hitsfilter;
    TVector2 proj_pr = pma::GetProjectionToPlane(vtxv3, view, tpc, cryo);

    mf::LogVerbatim("ProjectionMatchinAlg") << "start filter out: ";
    FilterOutSmallParts(detProp, 2.0, hitsview, hitsfilter, proj_pr);
    mf::LogVerbatim("ProjectionMatchingAlg") << "after filter out";

    for (size_t h = 0; h < hitsfilter.size(); ++h)
      hitstrk.push_back(hitsfilter[h]);

    if (hitsfilter.size() >= 5) countviews++;

    ++view;
  }

  if (!hitstrk.size() || (countviews < 2)) {
    mf::LogWarning("ProjectionMatchinAlg") << "too few hits, segment not built";
    return 0;
  }

  // hits are prepared, finally segment is built

  pma::Track3D* trk = new pma::Track3D();
  trk = buildSegment(detProp, hitstrk, vtxv3);

  // make shorter segment to estimate direction more precise
  ShortenSeg_(detProp, *trk, tpcgeom);

  if (trk->size() < 10) {
    mf::LogWarning("ProjectionMatchingAlg") << "too short segment, delete segment";
    delete trk;
    return 0;
  }

  return trk;
}

// ------------------------------------------------------

void pma::ProjectionMatchingAlg::FilterOutSmallParts(
  const detinfo::DetectorPropertiesData& detProp,
  double r2d,
  const std::vector<art::Ptr<recob::Hit>>& hits_in,
  std::vector<art::Ptr<recob::Hit>>& hits_out,
  const TVector2& vtx2d) const
{
  size_t min_size = hits_in.size() / 5;
  if (min_size < 3) min_size = 3;

  std::vector<size_t> used;
  std::vector<std::vector<art::Ptr<recob::Hit>>> sub_groups;
  std::vector<art::Ptr<recob::Hit>> close_hits;

  float mindist2 = 99999.99;
  size_t id_sub_groups_save = 0;
  size_t id = 0;
  while (GetCloseHits_(detProp, r2d, hits_in, used, close_hits)) {

    sub_groups.emplace_back(close_hits);

    for (size_t h = 0; h < close_hits.size(); ++h) {
      TVector2 hi_cm = pma::WireDriftToCm(detProp,
                                          close_hits[h]->WireID().Wire,
                                          close_hits[h]->PeakTime(),
                                          close_hits[h]->WireID().Plane,
                                          close_hits[h]->WireID().TPC,
                                          close_hits[h]->WireID().Cryostat);

      float dist2 = pma::Dist2(hi_cm, vtx2d);
      if (dist2 < mindist2) {
        id_sub_groups_save = id;
        mindist2 = dist2;
      }
    }

    id++;
  }

  for (size_t i = 0; i < sub_groups.size(); ++i) {
    if (i == id_sub_groups_save) {
      for (auto h : sub_groups[i])
        hits_out.push_back(h);
    }
  }

  for (size_t i = 0; i < sub_groups.size(); ++i) {
    if ((i != id_sub_groups_save) && (hits_out.size() < 10) && (sub_groups[i].size() > min_size))
      for (auto h : sub_groups[i])
        hits_out.push_back(h);
  }
}

// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::GetCloseHits_(const detinfo::DetectorPropertiesData& detProp,
                                               double r2d,
                                               const std::vector<art::Ptr<recob::Hit>>& hits_in,
                                               std::vector<size_t>& used,
                                               std::vector<art::Ptr<recob::Hit>>& hits_out) const
{
  hits_out.clear();

  const double gapMargin = 5.0; // can be changed to f(id_tpc1, id_tpc2)
  size_t idx = 0;

  while ((idx < hits_in.size()) && Has_(used, idx))
    idx++;

  if (idx < hits_in.size()) {
    hits_out.push_back(hits_in[idx]);
    used.push_back(idx);

    unsigned int tpc = hits_in[idx]->WireID().TPC;
    unsigned int cryo = hits_in[idx]->WireID().Cryostat;
    unsigned int view = hits_in[idx]->WireID().Plane;
    double wirePitch = fGeom->Plane(geo::PlaneID(cryo, tpc, view)).WirePitch();
    double driftPitch = detProp.GetXTicksCoefficient(tpc, cryo);

    double r2d2 = r2d * r2d;
    double gapMargin2 = sqrt(2 * gapMargin * gapMargin);
    gapMargin2 = (gapMargin2 + r2d) * (gapMargin2 + r2d);

    bool collect = true;
    while (collect) {
      collect = false;
      for (size_t i = 0; i < hits_in.size(); i++)
        if (!Has_(used, i)) {
          art::Ptr<recob::Hit> const& hi = hits_in[i];
          TVector2 hi_cm(wirePitch * hi->WireID().Wire, driftPitch * hi->PeakTime());

          bool accept = false;

          for (auto const& ho : hits_out) {

            TVector2 ho_cm(wirePitch * ho->WireID().Wire, driftPitch * ho->PeakTime());
            double d2 = pma::Dist2(hi_cm, ho_cm);

            if (d2 < r2d2) {
              accept = true;
              break;
            }
          }
          if (accept) {
            collect = true;
            hits_out.push_back(hi);
            used.push_back(i);
          }
        }
    }
    return true;
  }
  else
    return false;
}

// ------------------------------------------------------

void pma::ProjectionMatchingAlg::ShortenSeg_(const detinfo::DetectorPropertiesData& detProp,
                                             pma::Track3D& trk,
                                             const geo::TPCGeo& tpcgeom) const
{
  double mse = trk.GetMse();
  mf::LogWarning("ProjectionMatchingAlg") << "initial value of mse: " << mse;
  while ((mse > 0.5) && TestTrk_(trk, tpcgeom)) {
    mse = trk.GetMse();
    for (size_t h = 0; h < trk.size(); ++h)
      if (std::sqrt(pma::Dist2(trk[h]->Point3D(), trk.front()->Point3D())) > 0.8 * trk.Length())
        trk[h]->SetEnabled(false);

    RemoveNotEnabledHits(trk);

    // trk.Optimize(0.0001, false); // BUG: first argument missing; tentatively:
    trk.Optimize(detProp, 0, 0.0001, false);
    trk.SortHits();

    mf::LogWarning("ProjectionMatchingAlg") << " mse: " << mse;
    if (mse == trk.GetMse()) break;

    mse = trk.GetMse();
  }

  RemoveNotEnabledHits(trk);
}

// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::TestTrk_(pma::Track3D& trk, const geo::TPCGeo& tpcgeom) const
{
  bool test = false;

  if (tpcgeom.HasPlane(0) && tpcgeom.HasPlane(1)) {
    if ((trk.NEnabledHits(0) > 5) && (trk.NEnabledHits(1) > 5)) test = true;
  }
  else if (tpcgeom.HasPlane(0) && tpcgeom.HasPlane(2)) {
    if ((trk.NEnabledHits(0) > 5) && (trk.NEnabledHits(2) > 5)) test = true;
  }
  else if (tpcgeom.HasPlane(1) && tpcgeom.HasPlane(2)) {
    if ((trk.NEnabledHits(1) > 5) && (trk.NEnabledHits(2) > 5)) test = true;
  }

  double length = 0.0;
  for (size_t h = 0; h < trk.size(); ++h) {
    if (!trk[h]->IsEnabled()) break;
    length = std::sqrt(pma::Dist2(trk.front()->Point3D(), trk[h]->Point3D()));
  }

  mf::LogWarning("ProjectionMatchingAlg") << "length of segment: " << length;
  if (length < 3.0) test = false; // cm

  return test;
}

// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::Has_(const std::vector<size_t>& v, size_t idx) const
{
  for (auto c : v)
    if (c == idx) return true;
  return false;
}

// ------------------------------------------------------

void pma::ProjectionMatchingAlg::RemoveNotEnabledHits(pma::Track3D& trk) const
{
  size_t h = 0;
  while (h < trk.size()) {
    if (trk[h]->IsEnabled())
      ++h;
    else
      (trk.release_at(h));
  }
}

// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
  const detinfo::DetectorPropertiesData& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits_1,
  const std::vector<art::Ptr<recob::Hit>>& hits_2) const
{
  pma::Track3D* trk = new pma::Track3D();
  trk->SetEndSegWeight(0.001F);
  trk->AddHits(detProp, hits_1);
  trk->AddHits(detProp, hits_2);

  if (trk->HasTwoViews() && (trk->TPCs().size() == 1)) // now only in single tpc
  {
    if (!trk->Initialize(detProp, 0.001F)) {
      mf::LogWarning("ProjectionMatchingAlg") << "initialization failed, delete segment";
      delete trk;
      return 0;
    }
    trk->Optimize(detProp, 0, fFineTuningEps);

    trk->SortHits();
    return trk;
  }
  else {
    mf::LogWarning("ProjectionMatchingAlg") << "need at least two views in single tpc";
    delete trk;
    return 0;
  }
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
  const detinfo::DetectorPropertiesData& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits_1,
  const std::vector<art::Ptr<recob::Hit>>& hits_2,
  const geo::Point_t& point) const
{
  pma::Track3D* trk = buildSegment(detProp, hits_1, hits_2);

  if (trk) {
    double dfront = pma::Dist2(trk->front()->Point3D(), point);
    double dback = pma::Dist2(trk->back()->Point3D(), point);
    if (dfront > dback) trk->Flip();

    trk->Nodes().front()->SetPoint3D({point.X(), point.Y(), point.Z()});
    trk->Nodes().front()->SetFrozen(true);
    trk->Optimize(detProp, 0, fFineTuningEps);

    trk->SortHits();
  }
  return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
  const detinfo::DetectorPropertiesData& detProp,
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const TVector3& point) const
{
  pma::Track3D* trk = buildSegment(detProp, hits);

  if (trk) {
    double dfront = pma::Dist2(trk->front()->Point3D(), point);
    double dback = pma::Dist2(trk->back()->Point3D(), point);
    if (dfront > dback) trk->Flip(); // detProp);

    trk->Nodes().front()->SetPoint3D(point);
    trk->Nodes().front()->SetFrozen(true);
    trk->Optimize(detProp, 0, fFineTuningEps);

    trk->SortHits();
  }
  return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::extendTrack(
  const detinfo::DetectorPropertiesData& detProp,
  const pma::Track3D& trk,
  const std::vector<art::Ptr<recob::Hit>>& hits,
  bool add_nodes) const
{
  pma::Track3D* copy = new pma::Track3D(trk);
  copy->AddHits(detProp, hits);

  mf::LogVerbatim("ProjectionMatchingAlg")
    << "ext. track size: " << copy->size() << "  #coll:" << copy->NHits(geo::kZ)
    << "  #ind2:" << copy->NHits(geo::kV) << "  #ind1:" << copy->NHits(geo::kU);

  if (add_nodes) {
    size_t nSegments = getSegCount_(copy->size());
    int nNodes = nSegments - copy->Nodes().size() + 1; // n nodes to add
    if (nNodes < 0) nNodes = 0;

    if (nNodes) {
      mf::LogVerbatim("ProjectionMatchingAlg") << "  add " << nNodes << " nodes";
      copy->Optimize(detProp, nNodes, fOptimizationEps);
    }
  }
  double g = copy->Optimize(detProp, 0, fFineTuningEps);
  mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt done, g = " << g;

  return copy;
}
// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::chkEndpointHits_(
  const detinfo::DetectorPropertiesData& detProp,
  int wire,
  int wdir,
  double drift_x,
  int view,
  unsigned int tpc,
  unsigned int cryo,
  const pma::Track3D& trk,
  const std::vector<art::Ptr<recob::Hit>>& hits) const
{
  size_t nCloseHits = 0;
  int forwWires = 3, backWires = -1;
  double xMargin = 0.4;
  for (auto h : hits)
    if ((view == (int)h->WireID().Plane) && (tpc == h->WireID().TPC) &&
        (cryo == h->WireID().Cryostat)) {
      bool found = false;
      for (size_t ht = 0; ht < trk.size(); ht++)
        if (trk[ht]->Hit2DPtr().key() == h.key()) {
          found = true;
          break;
        }
      if (found) continue;

      int dw = wdir * (h->WireID().Wire - wire);
      if ((dw <= forwWires) && (dw >= backWires)) {
        double x = detProp.ConvertTicksToX(h->PeakTime(), view, tpc, cryo);
        if (fabs(x - drift_x) < xMargin) nCloseHits++;
      }
    }
  if (nCloseHits > 1)
    return false;
  else
    return true;
}

bool pma::ProjectionMatchingAlg::addEndpointRef_(
  const detinfo::DetectorPropertiesData& detProp,
  pma::Track3D& trk,
  const std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>>& hits,
  std::pair<int, int> const* wires,
  double const* xPos,
  unsigned int tpc,
  unsigned int cryo) const
{
  double x = 0.0;
  std::vector<std::pair<int, unsigned int>> wire_view;
  for (unsigned int i = 0; i < 3; i++)
    if (wires[i].first >= 0) {
      const auto hiter = hits.find(i);
      if (hiter != hits.end()) {
        if (chkEndpointHits_(detProp,
                             wires[i].first,
                             wires[i].second,
                             xPos[i],
                             i,
                             tpc,
                             cryo,
                             trk,
                             hiter->second)) {
          x += xPos[i];
          wire_view.emplace_back(wires[i].first, i);
        }
      }
    }
  if (wire_view.size() > 1) {
    x /= wire_view.size();
    auto const [wire0, plane0] = wire_view[0];
    auto const [wire1, plane1] = wire_view[1];
    auto const [y, z, _] = fGeom
                             ->WireIDsIntersect(geo::WireID(cryo, tpc, plane0, wire0),
                                                geo::WireID(cryo, tpc, plane1, wire1))
                             .value_or(geo::WireIDIntersection::invalid());

    trk.AddRefPoint(x, y, z);
    mf::LogVerbatim("ProjectionMatchingAlg")
      << "trk tpc:" << tpc << " size:" << trk.size() << " add ref.point (" << x << "; " << y << "; "
      << z << ")";
    return true;
  }

  mf::LogVerbatim("ProjectionMatchingAlg")
    << "trk tpc:" << tpc << " size:" << trk.size()
    << " wire-plane-parallel track, but need two clean views of endpoint";
  return false;
}

void pma::ProjectionMatchingAlg::guideEndpoints(
  const detinfo::DetectorPropertiesData& detProp,
  pma::Track3D& trk,
  const std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>>& hits) const
{
  unsigned int tpc = trk.FrontTPC(), cryo = trk.FrontCryo();
  if ((tpc != trk.BackTPC()) || (cryo != trk.BackCryo())) {
    mf::LogWarning("ProjectionMatchingAlg") << "Please, apply before TPC stitching.";
    return;
  }

  const double maxCosXZ = 0.992546; // 7 deg

  pma::Segment3D* segFront = trk.Segments().front();
  if (trk.Segments().size() > 2) {
    pma::Segment3D* segFront1 = trk.Segments()[1];
    if ((segFront->Length() < 0.8) && (segFront1->Length() > 5.0)) segFront = segFront1;
  }
  pma::Vector3D dirFront = segFront->GetDirection3D();
  pma::Vector3D dirFrontXZ(dirFront.X(), 0., dirFront.Z());
  dirFrontXZ *= 1.0 / dirFrontXZ.R();

  pma::Segment3D* segBack = trk.Segments().back();
  if (trk.Segments().size() > 2) {
    pma::Segment3D* segBack1 = trk.Segments()[trk.Segments().size() - 2];
    if ((segBack->Length() < 0.8) && (segBack1->Length() > 5.0)) segBack = segBack1;
  }
  pma::Vector3D dirBack = segBack->GetDirection3D();
  pma::Vector3D dirBackXZ(dirBack.X(), 0., dirBack.Z());
  dirBackXZ *= 1.0 / dirBackXZ.R();

  if ((fabs(dirFrontXZ.Z()) < maxCosXZ) && (fabs(dirBackXZ.Z()) < maxCosXZ)) {
    return; // front & back are not parallel to wire planes => exit
  }

  unsigned int nPlanesFront = 0, nPlanesBack = 0;
  std::pair<int, int> wiresFront[3], wiresBack[3]; // wire index; index direction
  double xFront[3], xBack[3];

  for (unsigned int i = 0; i < 3; i++) {
    bool frontPresent = false, backPresent = false;
    if (fGeom->TPC(geo::TPCID(cryo, tpc)).HasPlane(i)) {
      int idxFront0 = trk.NextHit(-1, i);
      int idxBack0 = trk.PrevHit(trk.size(), i);
      if ((idxFront0 >= 0) && (idxFront0 < (int)trk.size()) && (idxBack0 >= 0) &&
          (idxBack0 < (int)trk.size())) {
        int idxFront1 = trk.NextHit(idxFront0, i);
        int idxBack1 = trk.PrevHit(idxBack0, i);
        if ((idxFront1 >= 0) && (idxFront1 < (int)trk.size()) && (idxBack1 >= 0) &&
            (idxBack1 < (int)trk.size())) {
          int wFront0 = trk[idxFront0]->Wire();
          int wBack0 = trk[idxBack0]->Wire();

          int wFront1 = trk[idxFront1]->Wire();
          int wBack1 = trk[idxBack1]->Wire();

          wiresFront[i].first = wFront0;
          wiresFront[i].second = wFront0 - wFront1;
          xFront[i] = detProp.ConvertTicksToX(trk[idxFront0]->PeakTime(), i, tpc, cryo);

          wiresBack[i].first = wBack0;
          wiresBack[i].second = wBack0 - wBack1;
          xBack[i] = detProp.ConvertTicksToX(trk[idxBack0]->PeakTime(), i, tpc, cryo);

          if (wiresFront[i].second) {
            if (wiresFront[i].second > 0)
              wiresFront[i].second = 1;
            else
              wiresFront[i].second = -1;

            frontPresent = true;
            nPlanesFront++;
          }

          if (wiresBack[i].second) {
            if (wiresBack[i].second > 0)
              wiresBack[i].second = 1;
            else
              wiresBack[i].second = -1;

            backPresent = true;
            nPlanesBack++;
          }
        }
      }
    }
    if (!frontPresent) { wiresFront[i].first = -1; }
    if (!backPresent) { wiresBack[i].first = -1; }
  }

  bool refAdded = false;
  if ((nPlanesFront > 1) && (fabs(dirFrontXZ.Z()) >= maxCosXZ)) {
    refAdded |= addEndpointRef_(detProp, trk, hits, wiresFront, xFront, tpc, cryo);
  }

  if ((nPlanesBack > 1) && (fabs(dirBackXZ.Z()) >= maxCosXZ)) {
    refAdded |= addEndpointRef_(detProp, trk, hits, wiresBack, xBack, tpc, cryo);
  }
  if (refAdded) {
    mf::LogVerbatim("ProjectionMatchingAlg") << "guide wire-plane-parallel track endpoints";
    double g = trk.Optimize(detProp, 0, 0.1 * fFineTuningEps);
    mf::LogVerbatim("ProjectionMatchingAlg") << "  done, g = " << g;
  }
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::guideEndpoints(
  const detinfo::DetectorPropertiesData& detProp,
  pma::Track3D& trk,
  pma::Track3D::ETrackEnd endpoint,
  const std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>>& hits) const
{
  const double maxCosXZ = 0.992546; // 7 deg

  unsigned int tpc, cryo;
  pma::Segment3D* seg0 = 0;
  pma::Segment3D* seg1 = 0;

  if (endpoint == pma::Track3D::kBegin) {
    tpc = trk.FrontTPC(), cryo = trk.FrontCryo();
    seg0 = trk.Segments().front();
    if (trk.Segments().size() > 2) { seg1 = trk.Segments()[1]; }
  }
  else {
    tpc = trk.BackTPC(), cryo = trk.BackCryo();
    seg0 = trk.Segments().back();
    if (trk.Segments().size() > 2) { seg1 = trk.Segments()[trk.Segments().size() - 2]; }
  }
  if (seg1 && (seg0->Length() < 0.8) && (seg1->Length() > 5.0)) { seg0 = seg1; }
  pma::Vector3D dir0 = seg0->GetDirection3D();
  pma::Vector3D dir0XZ(dir0.X(), 0., dir0.Z());
  dir0XZ *= 1.0 / dir0XZ.R();

  if (fabs(dir0XZ.Z()) < maxCosXZ) { return; } // not parallel to wire planes => exit

  unsigned int nPlanes = 0;
  std::pair<int, int> wires[3]; // wire index; index direction
  double x0[3];

  for (unsigned int i = 0; i < 3; i++) {
    bool present = false;
    if (fGeom->TPC(geo::TPCID(cryo, tpc)).HasPlane(i)) {
      int idx0 = -1, idx1 = -1;
      if (endpoint == pma::Track3D::kBegin) { idx0 = trk.NextHit(-1, i); }
      else {
        idx0 = trk.PrevHit(trk.size(), i);
      }

      if ((idx0 >= 0) && (idx0 < (int)trk.size()) && (trk[idx0]->TPC() == tpc) &&
          (trk[idx0]->Cryo() == cryo)) {
        if (endpoint == pma::Track3D::kBegin) { idx1 = trk.NextHit(idx0, i); }
        else {
          idx1 = trk.PrevHit(idx0, i);
        }

        if ((idx1 >= 0) && (idx1 < (int)trk.size()) && (trk[idx1]->TPC() == tpc) &&
            (trk[idx1]->Cryo() == cryo)) {
          int w0 = trk[idx0]->Wire();
          int w1 = trk[idx1]->Wire();

          wires[i].first = w0;
          wires[i].second = w0 - w1;
          x0[i] = detProp.ConvertTicksToX(trk[idx0]->PeakTime(), i, tpc, cryo);

          if (wires[i].second) {
            if (wires[i].second > 0)
              wires[i].second = 1;
            else
              wires[i].second = -1;

            present = true;
            nPlanes++;
          }
        }
      }
    }
    if (!present) { wires[i].first = -1; }
  }

  if ((nPlanes > 1) && (fabs(dir0XZ.Z()) >= maxCosXZ) &&
      addEndpointRef_(detProp, trk, hits, wires, x0, tpc, cryo)) {
    mf::LogVerbatim("ProjectionMatchingAlg") << "guide wire-plane-parallel track endpoint";
    double g = trk.Optimize(detProp, 0, 0.1 * fFineTuningEps);
    mf::LogVerbatim("ProjectionMatchingAlg") << "  done, g = " << g;
  }
}
// ------------------------------------------------------

std::vector<pma::Hit3D*> pma::ProjectionMatchingAlg::trimTrackToVolume(pma::Track3D&,
                                                                       TVector3,
                                                                       TVector3) const
{
  return {};
}
// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::alignTracks(pma::Track3D& first, pma::Track3D& second) const
{
  unsigned int k = 0;
  double const distFF = pma::Dist2(first.front()->Point3D(), second.front()->Point3D());
  double dist = distFF;

  double const distFB = pma::Dist2(first.front()->Point3D(), second.back()->Point3D());
  if (distFB < dist) {
    k = 1;
    dist = distFB;
  }

  double const distBF = pma::Dist2(first.back()->Point3D(), second.front()->Point3D());
  if (distBF < dist) {
    k = 2;
    dist = distBF;
  }

  double distBB = pma::Dist2(first.back()->Point3D(), second.back()->Point3D());
  if (distBB < dist) {
    k = 3;
    dist = distBB;
  }

  switch (k) // flip to get dst end before src start, do not merge if track's order reversed
  {
  case 0:
    first.Flip(); // detProp);
    break;
  case 1: mf::LogError("PMAlgTrackMaker") << "Tracks in reversed order."; return false;
  case 2: break;
  case 3:
    second.Flip(); // detProp);
    break;
  default:
    throw cet::exception("pma::ProjectionMatchingAlg")
      << "alignTracks: should never happen." << std::endl;
  }
  return true;
}

void pma::ProjectionMatchingAlg::mergeTracks(detinfo::DetectorPropertiesData const& detProp,
                                             pma::Track3D& dst,
                                             pma::Track3D& src,
                                             bool const reopt) const
{
  if (!alignTracks(dst, src)) return;

  unsigned int tpc = src.FrontTPC();
  unsigned int cryo = src.FrontCryo();
  double lmean = dst.Length() / (dst.Nodes().size() - 1);
  if ((pma::Dist2(dst.Nodes().back()->Point3D(), src.Nodes().front()->Point3D()) > 0.5 * lmean) ||
      (tpc != dst.BackTPC()) || (cryo != dst.BackCryo())) {
    dst.AddNode(detProp, src.Nodes().front()->Point3D(), tpc, cryo);
    if (src.Nodes().front()->IsFrozen()) dst.Nodes().back()->SetFrozen(true);
  }
  for (size_t n = 1; n < src.Nodes().size(); n++) {
    pma::Node3D* node = src.Nodes()[n];

    dst.AddNode(detProp, src.Nodes()[n]->Point3D(), node->TPC(), node->Cryo());

    if (node->IsFrozen()) dst.Nodes().back()->SetFrozen(true);
  }
  for (size_t h = 0; h < src.size(); h++) {
    dst.push_back(detProp, src[h]->Hit2DPtr());
  }
  if (reopt) {
    double g = dst.Optimize(detProp, 0, fFineTuningEps);
    mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt after merging done, g = " << g;
  }
  else {
    dst.MakeProjection();
  }

  dst.SortHits();
  dst.ShiftEndsToHits();

  dst.MakeProjection();
  dst.SortHits();
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::selectInitialHits(pma::Track3D& trk,
                                                     unsigned int view,
                                                     unsigned int* nused) const
{
  for (size_t i = 0; i < trk.size(); i++) {
    pma::Hit3D* hit = trk[i];
    if (hit->View2D() == view) {
      if ((hit->GetDistToProj() > 0.5) || // more than 0.5cm away away from the segment
          (hit->GetSegFraction() < -1.0)) // projects before segment start (to check!!!)
        hit->TagOutlier(true);
      else
        hit->TagOutlier(false);
    }
  }

  unsigned int nhits = 0;
  double last_x, dx = 0.0, last_q, dq = 0.0, dqdx = 0.0;
  int ih = trk.NextHit(-1, view);

  pma::Hit3D* hit = trk[ih];
  pma::Hit3D* lastHit = hit;

  if ((ih >= 0) && (ih < (int)trk.size())) {
    hit->TagOutlier(true);

    ih = trk.NextHit(ih, view);
    while ((dx < 2.5) && (ih >= 0) && (ih < (int)trk.size())) {
      hit = trk[ih];

      if (lar::util::absDiff(hit->Wire(), lastHit->Wire()) > 2)
        break; // break on gap in wire direction

      last_x = trk.HitDxByView(ih, view);
      last_q = hit->SummedADC();
      if (dx + last_x < 3.0) {
        dq += last_q;
        dx += last_x;
        nhits++;
      }
      else
        break;

      lastHit = hit;
      ih = trk.NextHit(ih, view);
    }
    while ((ih >= 0) && (ih < (int)trk.size())) {
      hit = trk[ih];
      hit->TagOutlier(true);
      ih = trk.NextHit(ih, view);
    }
  }
  else {
    mf::LogError("ProjectionMatchingAlg") << "Initial part selection failed.";
  }

  if (!nhits) {
    mf::LogError("ProjectionMatchingAlg") << "Initial part too short to select useful hits.";
  }

  if (dx > 0.0) dqdx = dq / dx;

  if (nused) (*nused) = nhits;

  return dqdx;
}
// ------------------------------------------------------
