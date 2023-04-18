/**
 *  @file   PmaTrack3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 *
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Build 3D segments and whole tracks by simultaneous matching hits in 2D projections.
 *          See PmaTrack3D.h file for details.
 */

#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/PmaSegment3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <array>

#include "range/v3/algorithm.hpp"
#include "range/v3/view.hpp"

pma::Track3D::Track3D() = default;

pma::Track3D::Track3D(const Track3D& src)
  : fMaxHitsPerSeg(src.fMaxHitsPerSeg)
  , fPenaltyFactor(src.fPenaltyFactor)
  , fMaxSegStopFactor(src.fMaxSegStopFactor)
  , fSegStopValue(src.fSegStopValue)
  , fMinSegStop(src.fMinSegStop)
  , fMaxSegStop(src.fMaxSegStop)
  , fSegStopFactor(src.fSegStopFactor)
  , fPenaltyValue(src.fPenaltyValue)
  , fEndSegWeight(src.fEndSegWeight)
  , fHitsRadius(src.fHitsRadius)
  , fT0(src.fT0)
  , fT0Flag(src.fT0Flag)
  , fTag(src.fTag)
{
  fHits.reserve(src.fHits.size());
  for (auto const* hit : src.fHits) {
    pma::Hit3D* h3d = new pma::Hit3D(*hit);
    h3d->fParent = this;
    fHits.push_back(h3d);
  }

  fNodes.reserve(src.fNodes.size());
  for (auto const* node : src.fNodes)
    fNodes.push_back(new pma::Node3D(*node));

  for (auto const* point : src.fAssignedPoints)
    fAssignedPoints.push_back(new TVector3(*point));

  RebuildSegments();
  MakeProjection();
}

pma::Track3D::~Track3D()
{
  for (size_t i = 0; i < fHits.size(); i++)
    delete fHits[i];
  for (size_t i = 0; i < fAssignedPoints.size(); i++)
    delete fAssignedPoints[i];

  for (size_t i = 0; i < fSegments.size(); i++)
    delete fSegments[i];
  for (size_t i = 0; i < fNodes.size(); i++)
    if (!fNodes[i]->NextCount() && !fNodes[i]->Prev()) delete fNodes[i];
}

bool pma::Track3D::Initialize(detinfo::DetectorPropertiesData const& detProp, float initEndSegW)
{
  if (!HasTwoViews(2)) {
    mf::LogError("pma::Track3D") << "Need min. 2 hits per view, at least two views.";
    return false;
  }

  auto cryos = Cryos();
  if (cryos.size() > 1) {
    mf::LogError("pma::Track3D") << "Only one cryostat for now, please.";
    return false;
  }
  int cryo = cryos.front();

  auto tpcs = TPCs();
  if (tpcs.size() > 1) {
    mf::LogError("pma::Track3D") << "Only one TPC, please.";
    return false;
  }
  // single tpc, many tpc's are ok, but need to be handled from ProjectionMatchingAlg::buildMultiTPCTrack()
  int tpc = tpcs.front();

  if (InitFromRefPoints(detProp, tpc, cryo))
    mf::LogVerbatim("pma::Track3D") << "Track initialized with 3D reference points.";
  else if (InitFromHits(detProp, tpc, cryo, initEndSegW))
    mf::LogVerbatim("pma::Track3D") << "Track initialized with hit positions.";
  else {
    InitFromMiddle(detProp, tpc, cryo);
    mf::LogVerbatim("pma::Track3D") << "Track initialized in the module center.";
  }

  UpdateHitsRadius();
  return true;
}

void pma::Track3D::ClearNodes()
{
  for (size_t i = 0; i < fNodes.size(); i++)
    delete fNodes[i];
  fNodes.clear();
}

bool pma::Track3D::InitFromHits(detinfo::DetectorPropertiesData const& detProp,
                                int tpc,
                                int cryo,
                                float initEndSegW)
{
  art::ServiceHandle<geo::Geometry const> geom;
  geo::TPCID const tpcid(cryo, tpc);

  float wtmp = fEndSegWeight;
  fEndSegWeight = initEndSegW;

  // endpoints for the first combination:
  TVector3 v3d_1(0., 0., 0.), v3d_2(0., 0., 0.);

  assert(!fHits.empty());

  pma::Hit3D* hit0_a = fHits.front();
  pma::Hit3D* hit1_a = fHits.front();

  geo::PlaneID const hit0_a_planeid{tpcid, hit0_a->View2D()};
  geo::PlaneID const hit1_a_planeid{tpcid, hit1_a->View2D()};

  float minX = detProp.ConvertTicksToX(hit0_a->PeakTime(), hit0_a_planeid);
  float maxX = detProp.ConvertTicksToX(hit1_a->PeakTime(), hit1_a_planeid);
  for (auto hit : fHits) {
    double const x = detProp.ConvertTicksToX(hit->PeakTime(), geo::PlaneID{tpcid, hit->View2D()});
    if (x < minX) {
      minX = x;
      hit0_a = hit;
    }
    if (x > maxX) {
      maxX = x;
      hit1_a = hit;
    }
  }

  pma::Hit3D* hit0_b = nullptr;
  pma::Hit3D* hit1_b = nullptr;

  float minDiff0 = 5000, minDiff1 = 5000;
  for (auto hit : fHits) {
    double const x = detProp.ConvertTicksToX(hit->PeakTime(), geo::PlaneID{tpcid, hit->View2D()});
    double diff = fabs(x - minX);
    if ((diff < minDiff0) && (hit->View2D() != hit0_a->View2D())) {
      minDiff0 = diff;
      hit0_b = hit;
    }
    diff = fabs(x - maxX);
    if ((diff < minDiff1) && (hit->View2D() != hit1_a->View2D())) {
      minDiff1 = diff;
      hit1_b = hit;
    }
  }

  if (hit0_a && hit0_b && hit1_a && hit1_b) {
    geo::PlaneID const hit0_b_planeid{tpcid, hit0_b->View2D()};
    geo::PlaneID const hit1_b_planeid{tpcid, hit1_b->View2D()};
    double x = 0.5 * (detProp.ConvertTicksToX(hit0_a->PeakTime(), hit0_a_planeid) +
                      detProp.ConvertTicksToX(hit0_b->PeakTime(), hit0_b_planeid));
    double y, z;
    geom->IntersectionPoint(geo::WireID{hit0_a_planeid, hit0_a->Wire()},
                            geo::WireID{hit0_b_planeid, hit0_b->Wire()},
                            y,
                            z);
    v3d_1.SetXYZ(x, y, z);

    x = 0.5 * (detProp.ConvertTicksToX(hit1_a->PeakTime(), hit1_a_planeid) +
               detProp.ConvertTicksToX(hit1_b->PeakTime(), hit1_b_planeid));
    geom->IntersectionPoint(geo::WireID{hit1_a_planeid, hit1_a->Wire()},
                            geo::WireID{hit1_b_planeid, hit1_b->Wire()},
                            y,
                            z);
    v3d_2.SetXYZ(x, y, z);

    ClearNodes();
    AddNode(detProp, v3d_1, tpc, cryo);
    AddNode(detProp, v3d_2, tpc, cryo);

    MakeProjection();
    UpdateHitsRadius();
    Optimize(detProp, 0, 0.01F, false, true, 100);
    SelectAllHits();
  }
  else {
    mf::LogVerbatim("pma::Track3D") << "Good hits not found.";
    fEndSegWeight = wtmp;
    return false;
  }

  if (pma::Dist2(fNodes.front()->Point3D(), fNodes.back()->Point3D()) < (0.3 * 0.3)) {
    mf::LogVerbatim("pma::Track3D") << "Short initial segment.";
    fEndSegWeight = wtmp;
    return false;
  }

  fEndSegWeight = wtmp;
  return true;
}

bool pma::Track3D::InitFromRefPoints(detinfo::DetectorPropertiesData const& detProp,
                                     int tpc,
                                     int cryo)
{
  if (fAssignedPoints.size() < 2) return false;

  ClearNodes();

  TVector3 mean(0., 0., 0.), stdev(0., 0., 0.), p(0., 0., 0.);
  for (size_t i = 0; i < fAssignedPoints.size(); i++) {
    p = *(fAssignedPoints[i]);
    mean += p;
    p.SetXYZ(p.X() * p.X(), p.Y() * p.Y(), p.Z() * p.Z());
    stdev += p;
  }
  stdev *= 1.0 / fAssignedPoints.size();
  mean *= 1.0 / fAssignedPoints.size();
  p = mean;
  p.SetXYZ(p.X() * p.X(), p.Y() * p.Y(), p.Z() * p.Z());
  stdev -= p;

  double sx = stdev.X(), sy = stdev.Y(), sz = stdev.Z();
  if (sx >= 0.0)
    sx = sqrt(sx);
  else
    sx = 0.0;
  if (sy >= 0.0)
    sy = sqrt(sy);
  else
    sy = 0.0;
  if (sz >= 0.0)
    sz = sqrt(sz);
  else
    sz = 0.0;
  stdev.SetXYZ(sx, sy, sz);

  double scale = 2.0 * stdev.Mag();
  double iscale = 1.0 / scale;

  size_t max_index = 0;
  double norm2, max_norm2 = 0.0;
  std::vector<TVector3> data;
  for (size_t i = 0; i < fAssignedPoints.size(); i++) {
    p = *(fAssignedPoints[i]);
    p -= mean;
    p *= iscale;
    norm2 = p.Mag2();
    if (norm2 > max_norm2) {
      max_norm2 = norm2;
      max_index = i;
    }
    data.push_back(p);
  }

  double y = 0.0, kappa = 1.0, prev_kappa, kchg = 1.0;
  TVector3 w(data[max_index]);

  while (kchg > 0.0001)
    for (size_t i = 0; i < data.size(); i++) {
      y = (data[i] * w);
      w += (y / kappa) * (data[i] - y * w);

      prev_kappa = kappa;
      kappa += y * y;
      kchg = fabs((kappa - prev_kappa) / prev_kappa);
    }
  w *= 1.0 / w.Mag();

  TVector3 v1(w), v2(w);
  v1 *= scale;
  v1 += mean;
  v2 *= -scale;
  v2 += mean;
  std::sort(fAssignedPoints.begin(), fAssignedPoints.end(), pma::bSegmentProjLess(v1, v2));
  for (size_t i = 0; i < fAssignedPoints.size(); i++) {
    AddNode(detProp, *(fAssignedPoints[i]), tpc, cryo);
  }

  RebuildSegments();
  MakeProjection();

  if (size()) UpdateHitsRadius();

  Optimize(detProp, 0, 0.01F, false, true, 100);
  SelectAllHits();

  return true;
}

void pma::Track3D::InitFromMiddle(detinfo::DetectorPropertiesData const& detProp, int tpc, int cryo)
{
  art::ServiceHandle<geo::Geometry const> geom;

  const auto& tpcGeo = geom->TPC(geo::TPCID(cryo, tpc));

  double minX = tpcGeo.MinX(), maxX = tpcGeo.MaxX();
  double minY = tpcGeo.MinY(), maxY = tpcGeo.MaxY();
  double minZ = tpcGeo.MinZ(), maxZ = tpcGeo.MaxZ();

  TVector3 v3d_1(0.5 * (minX + maxX), 0.5 * (minY + maxY), 0.5 * (minZ + maxZ));
  TVector3 v3d_2(v3d_1);

  TVector3 shift(5.0, 5.0, 5.0);
  v3d_1 += shift;
  v3d_2 -= shift;

  ClearNodes();
  AddNode(detProp, v3d_1, tpc, cryo);
  AddNode(detProp, v3d_2, tpc, cryo);

  MakeProjection();
  UpdateHitsRadius();

  Optimize(detProp, 0, 0.01F);
}

int pma::Track3D::index_of(const pma::Node3D* n) const
{
  for (size_t i = 0; i < fNodes.size(); ++i)
    if (fNodes[i] == n) return (int)i;
  return -1;
}

int pma::Track3D::index_of(const pma::Hit3D* hit) const
{
  for (size_t i = 0; i < size(); i++)
    if (fHits[i] == hit) return (int)i;
  return -1;
}

pma::Hit3D* pma::Track3D::release_at(size_t index)
{
  pma::Hit3D* h3d = fHits[index];
  fHits.erase(fHits.begin() + index);
  return h3d;
}

bool pma::Track3D::push_back(const detinfo::DetectorPropertiesData& detProp,
                             const art::Ptr<recob::Hit>& hit)
{
  for (auto const& trk_hit : fHits) {
    if (trk_hit->fHit == hit) return false;
  }
  pma::Hit3D* h3d = new pma::Hit3D(detProp, hit);
  h3d->fParent = this;
  fHits.push_back(h3d);
  return true;
}

bool pma::Track3D::erase(const art::Ptr<recob::Hit>& hit)
{
  for (size_t i = 0; i < size(); i++) {
    if (hit == fHits[i]->fHit) {
      pma::Hit3D* h3d = fHits[i];
      fHits.erase(fHits.begin() + i);
      delete h3d;
      return true;
    }
  }
  return false;
}

pma::Vector3D pma::Track3D::GetDirection3D(size_t index) const
{
  pma::Hit3D* h = fHits[index];

  for (auto s : fSegments) {
    if (s->HasHit(h)) return s->GetDirection3D();
  }
  for (auto n : fNodes) {
    if (n->HasHit(h)) return n->GetDirection3D();
  }

  auto pe = GetNearestElement(h->Point2D(), h->View2D(), h->TPC());
  if (pe) {
    mf::LogWarning("pma::Track3D")
      << "GetDirection3D(): had to update hit assignment to segment/node.";
    pe->AddHit(h);
    return pe->GetDirection3D();
  }
  else {
    throw cet::exception("pma::Track3D") << "GetDirection3D(): direction of a not assigned hit "
                                         << index << " (size: " << fHits.size() << ")" << std::endl;
  }
}

void pma::Track3D::AddHits(detinfo::DetectorPropertiesData const& detProp,
                           const std::vector<art::Ptr<recob::Hit>>& hits)
{
  fHits.reserve(fHits.size() + hits.size());
  for (auto const& hit : hits)
    push_back(detProp, hit);
}

void pma::Track3D::RemoveHits(const std::vector<art::Ptr<recob::Hit>>& hits)
{
  unsigned int n = 0;
  for (auto const& hit : hits) {
    if (erase(hit)) n++;
  }
  if (n) MakeProjection();
}

unsigned int pma::Track3D::NHits(unsigned int view) const
{
  unsigned int n = 0;
  for (auto hit : fHits) {
    if (hit->View2D() == view) n++;
  }
  return n;
}

unsigned int pma::Track3D::NEnabledHits(unsigned int view) const
{
  unsigned int n = 0;
  for (auto hit : fHits) {
    if (hit->IsEnabled() && ((view == geo::kUnknown) || (view == hit->View2D()))) n++;
  }
  return n;
}

bool pma::Track3D::HasTwoViews(size_t const nmin) const
{
  unsigned int nviews = 0;
  if (NHits(geo::kU) >= nmin) nviews++;
  if (NHits(geo::kV) >= nmin) nviews++;
  if (NHits(geo::kZ) >= nmin) nviews++;
  return nviews > 1;
}

std::vector<unsigned int> pma::Track3D::TPCs() const
{
  std::vector<unsigned int> tpc_idxs;
  for (auto hit : fHits) {
    unsigned int tpc = hit->TPC();

    bool found = false;
    for (unsigned int const tpc_idx : tpc_idxs)
      if (tpc_idx == tpc) {
        found = true;
        break;
      }

    if (!found) tpc_idxs.push_back(tpc);
  }
  return tpc_idxs;
}

std::vector<unsigned int> pma::Track3D::Cryos() const
{
  std::vector<unsigned int> cryo_idxs;
  for (auto hit : fHits) {
    unsigned int cryo = hit->Cryo();

    bool found = false;
    for (size_t j = 0; j < cryo_idxs.size(); j++)
      if (cryo_idxs[j] == cryo) {
        found = true;
        break;
      }

    if (!found) cryo_idxs.push_back(cryo);
  }
  return cryo_idxs;
}

std::pair<TVector2, TVector2> pma::Track3D::WireDriftRange(
  detinfo::DetectorPropertiesData const& detProp,
  unsigned int view,
  unsigned int tpc,
  unsigned int cryo) const
{
  std::pair<TVector2, TVector2> range(TVector2(0., 0.), TVector2(0., 0.));

  size_t n0 = 0;
  while ((n0 < fNodes.size()) && (fNodes[n0]->TPC() != (int)tpc))
    n0++;

  if (n0 < fNodes.size()) {
    size_t n1 = n0;
    while ((n1 < fNodes.size()) && (fNodes[n1]->TPC() == (int)tpc))
      n1++;

    if (n0 > 0) n0--;
    if (n1 == fNodes.size()) n1--;

    TVector2 p0 = fNodes[n0]->Projection2D(view);
    p0 = pma::CmToWireDrift(detProp, p0.X(), p0.Y(), view, tpc, cryo);

    TVector2 p1 = fNodes[n1]->Projection2D(view);
    p1 = pma::CmToWireDrift(detProp, p1.X(), p1.Y(), view, tpc, cryo);

    if (p0.X() > p1.X()) {
      double tmp = p0.X();
      p0.Set(p1.X(), p0.Y());
      p1.Set(tmp, p1.Y());
    }
    if (p0.Y() > p1.Y()) {
      double tmp = p0.Y();
      p0.Set(p0.X(), p1.Y());
      p1.Set(p1.X(), tmp);
    }

    range.first = p0;
    range.second = p1;
  }
  return range;
}

bool pma::Track3D::Flip(detinfo::DetectorPropertiesData const& detProp,
                        std::vector<pma::Track3D*>& allTracks)
{
  if (fNodes.size() < 2) { return true; }

  std::vector<pma::Track3D*> toSort;

  pma::Node3D* n = fNodes.front();
  if (n->Prev()) {
    pma::Segment3D* s = static_cast<pma::Segment3D*>(n->Prev());
    pma::Track3D* t = s->Parent();

    if (t->NextSegment(n)) // starts from middle of another track: need to break that one first
    {
      int idx = t->index_of(n);
      if (idx >= 0) {
        pma::Track3D* u = t->Split(detProp, idx, false); // u is in front of t
        if (u) {
          allTracks.push_back(u);
          if (u->Flip(detProp, allTracks)) {
            InternalFlip(toSort);
            toSort.push_back(this);
          }
          else {
            mf::LogWarning("pma::Track3D")
              << "Flip(): Could not flip after split (but new track is preserved!!).";
            return false;
          }
        }
        else {
          return false;
        } // could not flip/break associated tracks, so give up on this one
      }
      else {
        throw cet::exception("pma::Track3D") << "Node not found." << std::endl;
      }
    }
    else // flip root
    {
      if (t->Flip(detProp, allTracks)) // all ok, so can flip this track
      {
        InternalFlip(toSort);
        toSort.push_back(this);
      }
      else {
        return false;
      } // could not flip/break associated tracks, so give up on this one
    }
  }
  else // simply flip
  {
    InternalFlip(toSort);
    toSort.push_back(this);
  }

  for (size_t t = 0; t < toSort.size(); t++) {
    bool sorted = false;
    for (size_t u = 0; u < t; u++)
      if (toSort[u] == toSort[t]) {
        sorted = true;
        break;
      }
    if (!sorted) {
      toSort[t]->MakeProjection();
      toSort[t]->SortHits();
    }
  }
  return true;
}

void pma::Track3D::InternalFlip(std::vector<pma::Track3D*>& toSort)
{
  for (size_t i = 0; i < fNodes.size() - 1; i++)
    if (fNodes[i]->NextCount() > 1) {
      for (size_t j = 0; j < fNodes[i]->NextCount(); j++) {
        pma::Segment3D* s = static_cast<pma::Segment3D*>(fNodes[i]->Next(j));
        if (s->Parent() != this) toSort.push_back(s->Parent());
      }
    }

  if (fNodes.back()->NextCount())
    for (size_t j = 0; j < fNodes.back()->NextCount(); j++) {
      pma::Segment3D* s = static_cast<pma::Segment3D*>(fNodes.back()->Next(j));
      toSort.push_back(s->Parent());
    }

  if (fNodes.front()->Prev()) {
    pma::Segment3D* s = static_cast<pma::Segment3D*>(fNodes.front()->Prev());
    toSort.push_back(s->Parent());
    s->Parent()->InternalFlip(toSort);
  }

  std::reverse(fNodes.begin(), fNodes.end());
  toSort.push_back(this);
  RebuildSegments();
}

void pma::Track3D::Flip()
{
  std::vector<pma::Track3D*> toSort;
  InternalFlip(toSort);
  toSort.push_back(this);

  for (size_t t = 0; t < toSort.size(); t++) {
    bool sorted = false;
    for (size_t u = 0; u < t; u++)
      if (toSort[u] == toSort[t]) {
        sorted = true;
        break;
      }
    if (!sorted) {
      toSort[t]->MakeProjection();
      toSort[t]->SortHits();
    }
  }
}

bool pma::Track3D::CanFlip() const
{
  if (!fNodes.size()) { return false; }

  pma::Node3D* n = fNodes.front();
  if (n->Prev()) {
    pma::Segment3D* s = static_cast<pma::Segment3D*>(n->Prev());
    pma::Track3D* t = s->Parent();
    if (t->NextSegment(n)) { return false; } // cannot flip if starts from middle of another track
    else {
      return t->CanFlip();
    } // check root
  }
  else {
    return true;
  }
}

void pma::Track3D::AutoFlip(pma::Track3D::EDirection dir, double thr, unsigned int n)
{
  unsigned int nViews = 3;
  pma::dedx_map dedxInViews[3];
  for (unsigned int i = 0; i < nViews; i++) {
    GetRawdEdxSequence(dedxInViews[i], i, 1);
  }
  unsigned int bestView = 2;
  if (dedxInViews[0].size() > 2 * dedxInViews[2].size()) bestView = 0;
  if (dedxInViews[1].size() > 2 * dedxInViews[2].size()) bestView = 1;

  std::vector<std::vector<double>> dedx;
  for (size_t i = 0; i < size(); i++) {
    auto it = dedxInViews[bestView].find(i);
    if (it != dedxInViews[bestView].end()) { dedx.push_back(it->second); }
  }
  if (!dedx.empty()) dedx.pop_back();

  float dEdxStart = 0.0F, dEdxStop = 0.0F;
  float dEStart = 0.0F, dxStart = 0.0F;
  float dEStop = 0.0F, dxStop = 0.0F;
  if (dedx.size() > 4) {
    if (!n) // use default options
    {
      if (dedx.size() > 30)
        n = 12;
      else if (dedx.size() > 20)
        n = 8;
      else if (dedx.size() > 10)
        n = 4;
      else
        n = 3;
    }

    size_t k = (dedx.size() - 2) >> 1;
    if (n > k) n = k;

    for (size_t i = 1, j = 0; j < n; i++, j++) {
      dEStart += dedx[i][5];
      dxStart += dedx[i][6];
    }
    if (dxStart > 0.0F) dEdxStart = dEStart / dxStart;

    for (size_t i = dedx.size() - 2, j = 0; j < n; i--, j++) {
      dEStop += dedx[i][5];
      dxStop += dedx[i][6];
    }
    if (dxStop > 0.0F) dEdxStop = dEStop / dxStop;
  }
  else if (dedx.size() == 4) {
    dEStart = dedx[0][5] + dedx[1][5];
    dxStart = dedx[0][6] + dedx[1][6];
    dEStop = dedx[2][5] + dedx[3][5];
    dxStop = dedx[2][6] + dedx[3][6];
    if (dxStart > 0.0F) dEdxStart = dEStart / dxStart;
    if (dxStop > 0.0F) dEdxStop = dEStop / dxStop;
  }
  else if (dedx.size() > 1) {
    if (dedx.front()[2] > 0.0F) dEdxStart = dedx.front()[5] / dedx.front()[6];
    if (dedx.back()[2] > 0.0F) dEdxStop = dedx.back()[5] / dedx.back()[6];
  }
  else
    return;

  if ((dir == pma::Track3D::kForward) && ((1.0 + thr) * dEdxStop < dEdxStart)) {
    mf::LogVerbatim("pma::Track3D")
      << "Auto-flip fired (1), thr: " << (1.0 + thr) << ", value: " << dEdxStart / dEdxStop;
    Flip(); // particle stops at the end of the track
  }
  if ((dir == pma::Track3D::kBackward) && (dEdxStop > (1.0 + thr) * dEdxStart)) {
    mf::LogVerbatim("pma::Track3D")
      << "Auto-flip fired (2), thr: " << (1.0 + thr) << ", value: " << dEdxStop / dEdxStart;
    Flip(); // particle stops at the front of the track
  }
}

bool pma::Track3D::AutoFlip(detinfo::DetectorPropertiesData const& detProp,
                            std::vector<pma::Track3D*>& allTracks,
                            pma::Track3D::EDirection dir,
                            double thr,
                            unsigned int n)
{
  unsigned int nViews = 3;
  pma::dedx_map dedxInViews[3];
  for (unsigned int i = 0; i < nViews; i++) {
    GetRawdEdxSequence(dedxInViews[i], i, 1);
  }
  unsigned int bestView = 2;
  if (dedxInViews[0].size() > 2 * dedxInViews[2].size()) bestView = 0;
  if (dedxInViews[1].size() > 2 * dedxInViews[2].size()) bestView = 1;

  std::vector<std::vector<double>> dedx;
  for (size_t i = 0; i < size(); i++) {
    auto it = dedxInViews[bestView].find(i);
    if (it != dedxInViews[bestView].end()) { dedx.push_back(it->second); }
  }
  if (!dedx.empty()) dedx.pop_back();

  float dEdxStart = 0.0F, dEdxStop = 0.0F;
  float dEStart = 0.0F, dxStart = 0.0F;
  float dEStop = 0.0F, dxStop = 0.0F;
  if (dedx.size() > 4) {
    if (!n) // use default options
    {
      if (dedx.size() > 30)
        n = 12;
      else if (dedx.size() > 20)
        n = 8;
      else if (dedx.size() > 10)
        n = 4;
      else
        n = 3;
    }

    size_t k = (dedx.size() - 2) >> 1;
    if (n > k) n = k;

    for (size_t i = 1, j = 0; j < n; i++, j++) {
      dEStart += dedx[i][5];
      dxStart += dedx[i][6];
    }
    if (dxStart > 0.0F) dEdxStart = dEStart / dxStart;

    for (size_t i = dedx.size() - 2, j = 0; j < n; i--, j++) {
      dEStop += dedx[i][5];
      dxStop += dedx[i][6];
    }
    if (dxStop > 0.0F) dEdxStop = dEStop / dxStop;
  }
  else if (dedx.size() == 4) {
    dEStart = dedx[0][5] + dedx[1][5];
    dxStart = dedx[0][6] + dedx[1][6];
    dEStop = dedx[2][5] + dedx[3][5];
    dxStop = dedx[2][6] + dedx[3][6];
    if (dxStart > 0.0F) dEdxStart = dEStart / dxStart;
    if (dxStop > 0.0F) dEdxStop = dEStop / dxStop;
  }
  else if (dedx.size() > 1) {
    if (dedx.front()[2] > 0.0F) dEdxStart = dedx.front()[5] / dedx.front()[6];
    if (dedx.back()[2] > 0.0F) dEdxStop = dedx.back()[5] / dedx.back()[6];
  }
  else
    return false;

  bool done = false;
  if ((dir == pma::Track3D::kForward) && ((1.0 + thr) * dEdxStop < dEdxStart)) {
    mf::LogVerbatim("pma::Track3D")
      << "Auto-flip fired (1), thr: " << (1.0 + thr) << ", value: " << dEdxStart / dEdxStop;
    done = Flip(detProp, allTracks); // particle stops at the end of the track
  }
  if ((dir == pma::Track3D::kBackward) && (dEdxStop > (1.0 + thr) * dEdxStart)) {
    mf::LogVerbatim("pma::Track3D")
      << "Auto-flip fired (2), thr: " << (1.0 + thr) << ", value: " << dEdxStop / dEdxStart;
    done = Flip(detProp, allTracks); // particle stops at the front of the track
  }
  return done;
}

double pma::Track3D::TestHitsMse(const detinfo::DetectorPropertiesData& detProp,
                                 const std::vector<art::Ptr<recob::Hit>>& hits,
                                 bool normalized) const
{
  if (hits.empty()) {
    mf::LogWarning("pma::Track3D") << "TestHitsMse(): Empty cluster.";
    return -1.0;
  }

  double mse = 0.0;
  for (auto const& h : hits) {
    unsigned int cryo = h->WireID().Cryostat;
    unsigned int tpc = h->WireID().TPC;
    unsigned int view = h->WireID().Plane;
    unsigned int wire = h->WireID().Wire;
    float drift = h->PeakTime();

    mse += Dist2(pma::WireDriftToCm(detProp, wire, drift, view, tpc, cryo), view, tpc, cryo);
  }
  if (normalized)
    return mse / hits.size();
  else
    return mse;
}

unsigned int pma::Track3D::TestHits(detinfo::DetectorPropertiesData const& detProp,
                                    const std::vector<art::Ptr<recob::Hit>>& hits,
                                    double dist) const
{
  if (hits.empty()) {
    mf::LogWarning("pma::Track3D") << "TestHits(): Empty cluster.";
    return 0;
  }

  TVector3 p3d;
  double tst, d2 = dist * dist;
  unsigned int nhits = 0;
  for (auto const& h : hits)
    if (GetUnconstrainedProj3D(detProp, h, p3d, tst) && tst < d2) ++nhits;

  return nhits;
}

double pma::Track3D::Length(size_t start, size_t stop, size_t step) const
{
  auto const nHits = size();
  if (nHits < 2) return 0.0;

  if (start > stop) {
    size_t tmp = stop;
    stop = start;
    start = tmp;
  }
  if (start >= nHits - 1) return 0.0;
  if (stop >= nHits) stop = nHits - 1;

  double result = 0.0;
  pma::Hit3D *h1 = nullptr, *h0 = fHits[start];

  if (!step) step = 1;
  size_t i = start + step;
  while (i <= stop) {
    h1 = fHits[i];
    result += sqrt(pma::Dist2(h1->Point3D(), h0->Point3D()));
    h0 = h1;
    i += step;
  }
  if (i - step < stop) // last step jumped beyond the stop point
  {                    // so need to add the last short segment
    result += sqrt(pma::Dist2(h0->Point3D(), back()->Point3D()));
  }
  return result;
}

int pma::Track3D::NextHit(int index, unsigned int view, bool inclDisabled) const
{
  pma::Hit3D* hit = nullptr;
  if (index < -1) index = -1;
  while (++index < (int)size()) // look for the next index of hit from the view
  {
    hit = fHits[index];
    if (hit->View2D() == view) {
      if (inclDisabled)
        break;
      else if (hit->IsEnabled())
        break;
    }
  }
  return index;
}

int pma::Track3D::PrevHit(int index, unsigned int view, bool inclDisabled) const
{
  pma::Hit3D* hit = nullptr;
  if (index > (int)size()) index = (int)size();
  while (--index >= 0) // look for the prev index of hit from the view
  {
    hit = fHits[index];
    if (hit->View2D() == view) {
      if (inclDisabled)
        break;
      else if (hit->IsEnabled())
        break;
    }
  }
  return index;
}

double pma::Track3D::HitDxByView(size_t index,
                                 unsigned int view,
                                 pma::Track3D::EDirection dir,
                                 bool secondDir) const
{
  pma::Hit3D* nexthit = nullptr;
  pma::Hit3D* hit = fHits[index];

  if (hit->View2D() != view) {
    mf::LogWarning("pma::Track3D") << "Function used with the hit not matching specified view.";
  }

  double dx = 0.0; // [cm]
  bool hitFound = false;
  int i = index;
  switch (dir) {
  case pma::Track3D::kForward:
    while (!hitFound && (++i < (int)size())) {
      nexthit = fHits[i];
      dx += sqrt(pma::Dist2(hit->Point3D(), nexthit->Point3D()));

      if (nexthit->View2D() == view)
        hitFound = true;
      else
        hitFound = false;

      hit = nexthit;
    }
    if (!hitFound) {
      if (!secondDir)
        dx = 0.5 * HitDxByView(index, view, pma::Track3D::kBackward, true);
      else {
        dx = Length();
        mf::LogWarning("pma::Track3D") << "Single hit in this view.";
      }
    }
    break;

  case pma::Track3D::kBackward:
    while (!hitFound && (--i >= 0)) {
      nexthit = fHits[i];
      dx += sqrt(pma::Dist2(hit->Point3D(), nexthit->Point3D()));

      if (nexthit->View2D() == view)
        hitFound = true;
      else
        hitFound = false;

      hit = nexthit;
    }
    if (!hitFound) {
      if (!secondDir)
        dx = 0.5 * HitDxByView(index, view, pma::Track3D::kForward, true);
      else {
        dx = Length();
        mf::LogWarning("pma::Track3D") << "Single hit in this view.";
      }
    }
    break;

  default: mf::LogError("pma::Track3D") << "Direction undefined."; break;
  }
  return dx;
}

double pma::Track3D::HitDxByView(size_t index, unsigned int view) const
{
  if (index < size()) {
    return 0.5 * (HitDxByView(index, view, pma::Track3D::kForward) +
                  HitDxByView(index, view, pma::Track3D::kBackward));
  }
  else {
    mf::LogError("pma::Track3D") << "Hit index out of range.";
    return 0.0;
  }
}

pma::Segment3D* pma::Track3D::NextSegment(pma::Node3D* vtx) const
{
  pma::Segment3D* seg = nullptr;
  unsigned int nCount = vtx->NextCount();
  unsigned int k = 0;
  while (k < nCount) {
    seg = static_cast<pma::Segment3D*>(vtx->Next(k));
    if (seg && (seg->Parent() == this)) return seg;
    k++;
  }
  return 0;
}

pma::Segment3D* pma::Track3D::PrevSegment(pma::Node3D* vtx) const
{
  if (vtx->Prev()) {
    auto seg = static_cast<pma::Segment3D*>(vtx->Prev());
    if (seg->Parent() == this) return seg;
  }
  return nullptr;
}

double pma::Track3D::GetRawdEdxSequence(std::map<size_t, std::vector<double>>& dedx,
                                        unsigned int view,
                                        unsigned int skip,
                                        bool inclDisabled) const
{
  dedx.clear();

  if (!size()) return 0.0;

  size_t step = 1;

  pma::Hit3D* hit = nullptr;

  double rv, dr, dR, dq, dEq, qSkipped = 0.0;

  size_t j = NextHit(-1, view, inclDisabled), s = skip;
  if (j >= size()) return 0.0F; // no charged hits at all
  while (j < size())            // look for the first hit index
  {
    hit = fHits[j];
    dq = hit->SummedADC();
    if (s) {
      qSkipped += dq;
      s--;
    }
    else
      break;

    j = NextHit(j, view, inclDisabled);
  }

  size_t jmax = PrevHit(size(), view, inclDisabled);

  std::vector<size_t> indexes;
  TVector3 p0(0., 0., 0.), p1(0., 0., 0.);
  TVector2 c0(0., 0.), c1(0., 0.);
  while (j <= jmax) {
    indexes.clear(); // prepare to collect hit indexes used for this dE/dx entry

    indexes.push_back(j);
    hit = fHits[j];

    p0 = hit->Point3D();
    p1 = hit->Point3D();

    c0.Set(hit->Wire(), hit->PeakTime());
    c1.Set(hit->Wire(), hit->PeakTime());

    dEq = hit->SummedADC(); // [now it is ADC sum]

    dr = HitDxByView(j,
                     view,
                     pma::Track3D::kForward); // protection against hits on the same position
    rv = HitDxByView(j, view);                // dx seen by j-th hit
    fHits[j]->fDx = rv;
    dR = rv;

    size_t m = 1; // number of hits with charge > 0
    while (((m < step) || (dR < 0.1) || (dr == 0.0)) && (j <= jmax)) {
      j = NextHit(j, view); // just next, even if tagged as outlier
      if (j > jmax) break;  // no more hits in this view

      hit = fHits[j];
      if (!inclDisabled && !hit->IsEnabled()) {
        if (dr == 0.0)
          continue;
        else
          break;
      }
      indexes.push_back(j);

      p1 = hit->Point3D();

      c1.Set(hit->Wire(), hit->PeakTime());

      dq = hit->SummedADC();

      dEq += dq;

      dr = HitDxByView(j, view, pma::Track3D::kForward);
      rv = HitDxByView(j, view);
      fHits[j]->fDx = rv;
      dR += rv;
      m++;
    }
    p0 += p1;
    p0 *= 0.5;
    c0 += c1;
    c0 *= 0.5;

    double range = Length(0, j);

    std::vector<double> trk_section;
    trk_section.push_back(c0.X());
    trk_section.push_back(c0.Y());
    trk_section.push_back(p0.X());
    trk_section.push_back(p0.Y());
    trk_section.push_back(p0.Z());
    trk_section.push_back(dEq);
    trk_section.push_back(dR);
    trk_section.push_back(range);

    for (auto const idx : indexes)
      dedx[idx] = trk_section;

    j = NextHit(j, view, inclDisabled);
  }

  return qSkipped;
}

std::vector<float> pma::Track3D::DriftsOfWireIntersection(
  detinfo::DetectorPropertiesData const& detProp,
  unsigned int wire,
  unsigned int view) const
{
  std::vector<float> drifts;
  for (size_t i = 0; i < fNodes.size() - 1; i++) {
    int tpc = fNodes[i]->TPC(), cryo = fNodes[i]->Cryo();
    if ((tpc != fNodes[i + 1]->TPC()) || (cryo != fNodes[i + 1]->Cryo())) continue;

    TVector2 p0 = pma::CmToWireDrift(detProp,
                                     fNodes[i]->Projection2D(view).X(),
                                     fNodes[i]->Projection2D(view).Y(),
                                     view,
                                     fNodes[i]->TPC(),
                                     fNodes[i]->Cryo());
    TVector2 p1 = pma::CmToWireDrift(detProp,
                                     fNodes[i + 1]->Projection2D(view).X(),
                                     fNodes[i + 1]->Projection2D(view).Y(),
                                     view,
                                     fNodes[i + 1]->TPC(),
                                     fNodes[i + 1]->Cryo());

    if ((p0.X() - wire) * (p1.X() - wire) <= 0.0) {
      double frac_w = (wire - p0.X()) / (p1.X() - p0.X());
      double d = p0.Y() + frac_w * (p1.Y() - p0.Y());
      drifts.push_back((float)d);
    }
  }
  return drifts;
}

size_t pma::Track3D::CompleteMissingWires(detinfo::DetectorPropertiesData const& detProp,
                                          unsigned int view)
{
  int dPrev, dw, w, wx, wPrev, i = NextHit(-1, view);

  pma::Hit3D* hitPrev = nullptr;
  pma::Hit3D* hit = nullptr;

  std::vector<pma::Hit3D*> missHits;

  while (i < (int)size()) {
    hitPrev = hit;
    hit = fHits[i];

    if (hit->View2D() == view) {
      w = hit->Wire();
      if (hitPrev) {
        wPrev = hitPrev->Wire();
        dPrev = (int)hitPrev->PeakTime();
        if (abs(w - wPrev) > 1) {
          if (w > wPrev)
            dw = 1;
          else
            dw = -1;
          wx = wPrev + dw;
          int k = 1;
          while (wx != w) {
            int peakTime = 0;
            unsigned int iWire = wx;
            std::vector<float> drifts = DriftsOfWireIntersection(detProp, wx, view);
            if (drifts.size()) {
              peakTime = drifts[0];
              for (size_t d = 1; d < drifts.size(); d++)
                if (fabs(drifts[d] - dPrev) < fabs(peakTime - dPrev)) peakTime = drifts[d];
            }
            else {
              mf::LogVerbatim("pma::Track3D") << "Track does not intersect with the wire " << wx;
              break;
            }
            pma::Hit3D* hmiss =
              new pma::Hit3D(detProp, iWire, view, hit->TPC(), hit->Cryo(), peakTime, 1.0, 1.0);
            missHits.push_back(hmiss);
            wx += dw;
            k++;
          }
        }
      }
    }
    else
      hit = hitPrev;

    i = NextHit(i, view, true);
    while ((i < (int)size()) && (hit->TPC() != fHits[i]->TPC())) {
      hitPrev = hit;
      hit = fHits[i];
      i = NextHit(i, view, true);
    }
  }

  if (missHits.size()) {
    for (size_t hi = 0; hi < missHits.size(); hi++) {
      push_back(missHits[hi]);
    }

    MakeProjection();
    SortHits();
  }

  return missHits.size();
}

void pma::Track3D::AddNode(pma::Node3D* node)
{
  fNodes.push_back(node);
  if (fNodes.size() > 1) RebuildSegments();
}

bool pma::Track3D::AddNode(detinfo::DetectorPropertiesData const& detProp)
{
  pma::Segment3D* seg;
  pma::Segment3D* maxSeg = nullptr;

  size_t si = 0;
  while (si < fSegments.size()) {
    if (!fSegments[si]->IsFrozen()) {
      maxSeg = fSegments[si];
      break;
    }
    else
      si++;
  }
  if (!maxSeg) return false;

  unsigned int nHitsByView[3];
  unsigned int nHits, maxHits = 0;
  unsigned int vIndex = 0, segHits, maxSegHits = 0;
  float segLength, maxLength = maxSeg->Length();
  for (unsigned int i = si + 1; i < fNodes.size(); i++) {
    seg = static_cast<pma::Segment3D*>(fNodes[i]->Prev());
    if (seg->IsFrozen()) continue;

    nHitsByView[0] = seg->NEnabledHits(geo::kU);
    nHitsByView[1] = seg->NEnabledHits(geo::kV);
    nHitsByView[2] = seg->NEnabledHits(geo::kZ);
    segHits = nHitsByView[0] + nHitsByView[1] + nHitsByView[2];

    if (segHits < 15) {
      if ((nHitsByView[0] == 0) && ((nHitsByView[1] < 4) || (nHitsByView[2] < 4))) continue;
      if ((nHitsByView[1] == 0) && ((nHitsByView[0] < 4) || (nHitsByView[2] < 4))) continue;
      if ((nHitsByView[2] == 0) && ((nHitsByView[0] < 4) || (nHitsByView[1] < 4))) continue;
    }

    nHits = fNodes[i]->NEnabledHits() + seg->NEnabledHits() + fNodes[i - 1]->NEnabledHits();

    if (nHits > maxHits) {
      maxHits = nHits;
      maxLength = seg->Length();
      maxSegHits = segHits;
      maxSeg = seg;
      vIndex = i;
    }
    else if (nHits == maxHits) {
      segLength = seg->Length();
      if (segLength > maxLength) {
        maxLength = segLength;
        maxSegHits = segHits;
        maxSeg = seg;
        vIndex = i;
      }
    }
  }

  if (maxSegHits > 1) {
    maxSeg->SortHits();

    nHitsByView[0] = maxSeg->NEnabledHits(geo::kU);
    nHitsByView[1] = maxSeg->NEnabledHits(geo::kV);
    nHitsByView[2] = maxSeg->NEnabledHits(geo::kZ);

    unsigned int maxViewIdx = 2, midViewIdx = 2;
    if ((nHitsByView[2] >= nHitsByView[1]) && (nHitsByView[1] >= nHitsByView[0])) {
      maxViewIdx = 2;
      midViewIdx = 1;
    }
    else if ((nHitsByView[1] >= nHitsByView[2]) && (nHitsByView[2] >= nHitsByView[0])) {
      maxViewIdx = 1;
      midViewIdx = 2;
    }
    else if ((nHitsByView[0] >= nHitsByView[2]) && (nHitsByView[2] >= nHitsByView[1])) {
      maxViewIdx = 0;
      midViewIdx = 2;
    }
    else if ((nHitsByView[2] >= nHitsByView[0]) && (nHitsByView[0] >= nHitsByView[1])) {
      maxViewIdx = 2;
      midViewIdx = 0;
    }
    else if ((nHitsByView[0] >= nHitsByView[1]) && (nHitsByView[1] >= nHitsByView[2])) {
      maxViewIdx = 0;
      midViewIdx = 1;
    }
    else if ((nHitsByView[1] >= nHitsByView[0]) && (nHitsByView[0] >= nHitsByView[2])) {
      maxViewIdx = 1;
      midViewIdx = 0;
    }
    if (nHitsByView[midViewIdx] < 2) midViewIdx = maxViewIdx;

    if (nHitsByView[midViewIdx] < 2) {
      mf::LogVerbatim("pma::Track3D") << "AddNode(): too few hits.";
      return false;
    }

    unsigned int mHits[3] = {0, 0, 0};
    unsigned int halfIndex = (nHitsByView[midViewIdx] >> 1) - 1;
    unsigned int n = 0, i = 0, i0 = 0, i1 = 0;
    while (i < maxSeg->NHits() - 1) {
      if (maxSeg->Hit(i).IsEnabled()) {
        mHits[maxSeg->Hit(i).View2D()]++;
        if (maxSeg->Hit(i).View2D() == midViewIdx) {
          if (n == halfIndex) break;
          n++;
        }
      }
      i++;
    }

    i0 = i;
    i++;
    while ((i < maxSeg->NHits()) &&
           !((maxSeg->Hit(i).View2D() == midViewIdx) && maxSeg->Hit(i).IsEnabled())) {
      i++;
    }
    i1 = i;

    if (!nHitsByView[0]) {
      if (nHitsByView[1] && (mHits[1] < 2)) {
        mf::LogVerbatim("pma::Track3D") << "AddNode(): low Ind2 hits.";
        return false;
      }
      if (nHitsByView[2] && (mHits[2] < 2)) {
        mf::LogVerbatim("pma::Track3D") << "AddNode(): low Coll hits.";
        return false;
      }
    }

    maxSeg->SetProjection(maxSeg->Hit(i0));
    maxSeg->SetProjection(maxSeg->Hit(i1));

    unsigned int tpc = maxSeg->Hit(i0).TPC();
    unsigned int cryo = maxSeg->Hit(i0).Cryo();

    pma::Node3D* p = new pma::Node3D(detProp,
                                     (maxSeg->Hit(i0).Point3D() + maxSeg->Hit(i1).Point3D()) * 0.5,
                                     tpc,
                                     cryo,
                                     false,
                                     fNodes[vIndex]->GetDriftShift());

    fNodes.insert(fNodes.begin() + vIndex, p);

    maxSeg->AddNext(fNodes[vIndex]);

    seg = new pma::Segment3D(this, fNodes[vIndex], fNodes[vIndex + 1]);
    fSegments.insert(fSegments.begin() + vIndex, seg);

    return true;
  }
  else
    return false;
}

void pma::Track3D::InsertNode(detinfo::DetectorPropertiesData const& detProp,
                              TVector3 const& p3d,
                              size_t at_idx,
                              unsigned int tpc,
                              unsigned int cryo)
{
  pma::Node3D* vtx =
    new pma::Node3D(detProp, p3d, tpc, cryo, false, fNodes[at_idx]->GetDriftShift());
  fNodes.insert(fNodes.begin() + at_idx, vtx);

  if (fNodes.size() > 1) RebuildSegments();
}

bool pma::Track3D::RemoveNode(size_t idx)
{
  if ((fNodes.size() > 1) && (idx < fNodes.size())) {
    pma::Node3D* vtx = fNodes[idx];
    fNodes.erase(fNodes.begin() + idx);
    RebuildSegments();

    if (!vtx->NextCount() && !vtx->Prev()) { delete vtx; }

    return true;
  }
  else
    return false;
}

pma::Track3D* pma::Track3D::Split(detinfo::DetectorPropertiesData const& detProp,
                                  size_t idx,
                                  bool try_start_at_idx)
{
  if (!idx || (idx + 1 >= fNodes.size())) return 0;

  pma::Node3D* n = nullptr;
  pma::Track3D* t0 = new pma::Track3D();
  t0->fT0 = fT0;
  t0->fT0Flag = fT0Flag;
  t0->fTag = fTag;

  for (size_t i = 0; i < idx; ++i) {
    n = fNodes.front();
    n->ClearAssigned();

    pma::Segment3D* s = static_cast<pma::Segment3D*>(n->Prev());
    if (s && (s->Parent() == this)) s->RemoveNext(n);

    size_t k = 0;
    while (k < n->NextCount()) {
      s = static_cast<pma::Segment3D*>(n->Next(k));
      if (s->Parent() == this)
        n->RemoveNext(s);
      else
        k++;
    }

    fNodes.erase(fNodes.begin());
    t0->fNodes.push_back(n);
  }

  n = fNodes.front();
  t0->fNodes.push_back(
    new pma::Node3D(detProp, n->Point3D(), n->TPC(), n->Cryo(), false, n->GetDriftShift()));
  t0->RebuildSegments();
  RebuildSegments();

  size_t h = 0;
  while (h < size()) {
    pma::Hit3D* h3d = fHits[h];
    unsigned int view = h3d->View2D(), tpc = h3d->TPC(), cryo = h3d->Cryo();
    double dist2D_old = Dist2(h3d->Point2D(), view, tpc, cryo);
    double dist2D_new = t0->Dist2(h3d->Point2D(), view, tpc, cryo);

    if ((dist2D_new < dist2D_old) && t0->HasTPC(tpc))
      t0->push_back(release_at(h));
    else if (!HasTPC(tpc) && t0->HasTPC(tpc))
      t0->push_back(release_at(h));
    else
      h++;
  }

  bool passed = true;
  if (HasTwoViews() && t0->HasTwoViews()) {
    if (try_start_at_idx && t0->CanFlip()) {
      mf::LogVerbatim("pma::Track3D") << "  attach new t0 and this trk to a common start node";
      t0->Flip();
      passed = t0->AttachTo(fNodes.front());
    }
    else {
      mf::LogVerbatim("pma::Track3D") << "  attach this trk to the new t0 end node";
      passed = AttachTo(t0->fNodes.back());
    }
  }
  else {
    mf::LogVerbatim("pma::Track3D") << "  single-view track";
    passed = false;
  }

  if (!passed) {
    mf::LogVerbatim("pma::Track3D") << "  undo split";
    while (t0->size())
      push_back(t0->release_at(0));

    for (size_t i = 0; i < idx; ++i) {
      fNodes.insert(fNodes.begin() + i, t0->fNodes.front());
      t0->fNodes.erase(t0->fNodes.begin());
    }

    RebuildSegments();
    delete t0;
    t0 = 0;
  }

  return t0;
}

bool pma::Track3D::AttachTo(pma::Node3D* vStart, bool noFlip)
{
  pma::Node3D* vtx = fNodes.front();

  if (vtx == vStart) return true; // already connected!

  pma::Track3D* rootThis = GetRoot();
  if (vStart->Prev()) {
    pma::Track3D* rootPrev = static_cast<pma::Segment3D*>(vStart->Prev())->Parent()->GetRoot();
    if (rootPrev->IsAttachedTo(rootThis)) { return false; }
  }
  else if (vStart->NextCount()) {
    pma::Track3D* rootNext = static_cast<pma::Segment3D*>(vStart->Next(0))->Parent()->GetRoot();
    if (rootNext->IsAttachedTo(rootThis)) { return false; }
  }
  else {
    return false;
  }

  for (auto n : fNodes)
    if (n == vStart) {
      mf::LogError("pma::Track3D") << "Don't create loops!";
      return false;
    }

  if (!noFlip && CanFlip() && (vStart->TPC() == fNodes.back()->TPC()) &&
      (pma::Dist2(vtx->Point3D(), vStart->Point3D()) >
       pma::Dist2(fNodes.back()->Point3D(), vStart->Point3D())) &&
      (fNodes.back()->NextCount() == 0)) {
    mf::LogError("pma::Track3D") << "Flip, endpoint closer to vStart.";
    Flip();
  }

  if (vStart->TPC() == vtx->TPC())
    return AttachToSameTPC(vStart);
  else
    return AttachToOtherTPC(vStart);
}

bool pma::Track3D::AttachToOtherTPC(pma::Node3D* vStart)
{
  if (fNodes.front()->Prev()) return false;

  mf::LogVerbatim("pma::Track3D") << "Attach to track in another TPC.";

  fNodes.insert(fNodes.begin(), vStart);

  fSegments.insert(fSegments.begin(), new pma::Segment3D(this, fNodes[0], fNodes[1]));

  return true;
}

bool pma::Track3D::AttachToSameTPC(pma::Node3D* vStart)
{
  pma::Node3D* vtx = fNodes.front();

  if (vtx->Prev()) {
    pma::Segment3D* segPrev = static_cast<pma::Segment3D*>(vtx->Prev());
    pma::Track3D* tpPrev = segPrev->Parent();
    if (tpPrev->NextSegment(vtx)) { return false; }
    else if (tpPrev->CanFlip()) {
      tpPrev->Flip();
    } // flip in local vtx, no problem
    else {
      if (vStart->Prev()) {
        pma::Segment3D* segNew = static_cast<pma::Segment3D*>(vStart->Prev());
        pma::Track3D* tpNew = segNew->Parent();
        if (tpNew->NextSegment(vStart)) { return false; }
        else if (tpNew->CanFlip()) // flip in remote vStart, no problem
        {
          tpNew->Flip();

          mf::LogVerbatim("pma::Track3D") << "Reconnect prev to vStart.";
          tpPrev->fNodes[tpPrev->fNodes.size() - 1] = vStart;
          segPrev->AddNext(vStart);
        }
        else {
          mf::LogError("pma::Track3D") << "Flips in vtx and vStart not possible, cannot attach.";
          return false;
        }
      }
      else {
        mf::LogVerbatim("pma::Track3D") << "Reconnect prev to vStart.";
        tpPrev->fNodes[tpPrev->fNodes.size() - 1] = vStart;
        segPrev->AddNext(vStart);
      }
    }
  }

  while (vtx->NextCount()) // reconnect nexts to vStart
  {
    pma::Segment3D* seg = static_cast<pma::Segment3D*>(vtx->Next(0));
    pma::Track3D* trk = seg->Parent();

    vtx->RemoveNext(seg);
    trk->fNodes[0] = vStart;
    vStart->AddNext(seg);
  }

  if (vtx->NextCount() || vtx->Prev()) // better throw here
  {
    throw cet::exception("pma::Track3D") << "Something is still using disconnected vertex.";
  }
  else
    delete vtx; // ok
  return true;
}

bool pma::Track3D::AttachBackTo(pma::Node3D* vStart)
{
  pma::Node3D* vtx = fNodes.back();

  if (vtx == vStart) return true; // already connected!

  pma::Track3D* rootThis = GetRoot();
  if (vStart->Prev()) {
    pma::Track3D* rootPrev = static_cast<pma::Segment3D*>(vStart->Prev())->Parent()->GetRoot();
    if (rootPrev->IsAttachedTo(rootThis)) { return false; }
  }
  else if (vStart->NextCount()) {
    pma::Track3D* rootNext = static_cast<pma::Segment3D*>(vStart->Next(0))->Parent()->GetRoot();
    if (rootNext->IsAttachedTo(rootThis)) { return false; }
  }
  else {
    return false;
  }

  for (auto n : fNodes)
    if (n == vStart) {
      mf::LogError("pma::Track3D") << "Don't create loops!";
      return false;
    }

  if (vStart->TPC() == vtx->TPC())
    return AttachBackToSameTPC(vStart);
  else
    return AttachBackToOtherTPC(vStart);
}

bool pma::Track3D::AttachBackToOtherTPC(pma::Node3D* vStart)
{
  if (vStart->Prev()) return false;

  mf::LogVerbatim("pma::Track3D") << "Attach to track in another TPC.";

  fNodes.push_back(vStart);

  size_t idx = fNodes.size() - 1;
  fSegments.push_back(new pma::Segment3D(this, fNodes[idx - 1], fNodes[idx]));

  return true;
}

bool pma::Track3D::AttachBackToSameTPC(pma::Node3D* vStart)
{
  pma::Node3D* vtx = fNodes.back();

  if (vStart->Prev()) {
    pma::Segment3D* seg = static_cast<pma::Segment3D*>(vStart->Prev());
    pma::Track3D* tp = seg->Parent();
    if (tp->NextSegment(vStart)) {
      mf::LogError("pma::Track3D") << "Cannot attach back to inner node of other track.";
      return false;
    }

    if (tp->CanFlip())
      tp->Flip(); // flip in remote vStart, no problem
    else {
      mf::LogError("pma::Track3D") << "Flip not possible, cannot attach.";
      return false;
    }
  }
  fNodes[fNodes.size() - 1] = vStart;
  fSegments[fSegments.size() - 1]->AddNext(vStart);

  while (vtx->NextCount()) // reconnect nexts to vStart
  {
    pma::Segment3D* seg = static_cast<pma::Segment3D*>(vtx->Next(0));
    pma::Track3D* trk = seg->Parent();

    vtx->RemoveNext(seg);
    trk->fNodes[0] = vStart;
    vStart->AddNext(seg);
  }

  if (vtx->NextCount() || vtx->Prev()) {
    throw cet::exception("pma::Track3D")
      << "Something is still using disconnected vertex." << std::endl;
  }
  else
    delete vtx; // ok

  return true;
}

void pma::Track3D::ExtendWith(pma::Track3D* src)
{
  if (src->fNodes.front()->Prev()) {
    throw cet::exception("pma::Track3D")
      << "Cant extend with a track being a daughter of another track." << std::endl;
  }

  src->DeleteSegments(); // disconnect from vertices and delete
  for (size_t i = 0; i < src->fNodes.size(); ++i) {
    fNodes.push_back(src->fNodes[i]);

    size_t idx = fNodes.size() - 1;
    fSegments.push_back(new pma::Segment3D(this, fNodes[idx - 1], fNodes[idx]));
  }
  src->fNodes.clear(); // just empty the node collection

  while (src->size()) {
    push_back(src->release_at(0));
  } // move hits

  for (auto p : src->fAssignedPoints) {
    fAssignedPoints.push_back(p);
  } // move 3D reference points
  src->fAssignedPoints.clear();

  MakeProjection();
  SortHits();

  delete src;
}

pma::Track3D* pma::Track3D::GetRoot()
{
  pma::Track3D* trk = nullptr;

  if (fNodes.size()) {
    pma::Segment3D* seg = static_cast<pma::Segment3D*>(fNodes.front()->Prev());
    if (seg) trk = seg->Parent()->GetRoot();
    if (!trk) trk = this;
  }
  else
    throw cet::exception("pma::Track3D") << "Broken track, no nodes.";

  return trk;
}

bool pma::Track3D::GetBranches(std::vector<pma::Track3D const*>& branches, bool skipFirst) const
{
  for (auto trk : branches)
    if (trk == this) { throw cet::exception("pma::Track3D") << "Track tree has loop."; }

  branches.push_back(this);

  size_t i0 = 0;
  if (skipFirst) i0 = 1;

  for (size_t i = i0; i < fNodes.size(); ++i)
    for (size_t n = 0; n < fNodes[i]->NextCount(); n++) {
      pma::Segment3D* seg = static_cast<pma::Segment3D*>(fNodes[i]->Next(n));
      if (seg && (seg->Parent() != this)) seg->Parent()->GetBranches(branches, true);
    }

  return true;
}

bool pma::Track3D::IsAttachedTo(pma::Track3D const* trk) const
{
  if (trk == this) return true;

  std::vector<pma::Track3D const*> branchesThis, branchesTrk;

  if (!trk->GetBranches(branchesTrk)) return true; // has loop - better check the reason
  if (!GetBranches(branchesThis)) return true;     // has loop - better check the reason

  for (auto bThis : branchesThis)
    for (auto bTrk : branchesTrk)
      if (bThis == bTrk) return true;
  return false;
}

bool pma::Track3D::HasRefPoint(TVector3* p) const
{
  for (auto point : fAssignedPoints)
    if (point == p) return true;
  return false;
}

double pma::Track3D::GetMse(unsigned int view) const
{
  double sumMse = 0.0;
  unsigned int nEnabledHits = 0;
  for (auto n : fNodes) {
    sumMse += n->SumDist2(view);
    nEnabledHits += n->NEnabledHits(view);
  }
  for (auto s : fSegments) {
    sumMse += s->SumDist2(view);
    nEnabledHits += s->NEnabledHits(view);
  }

  if (nEnabledHits)
    return sumMse / nEnabledHits;
  else
    return 0.0;
}

double pma::Track3D::GetObjFunction(float penaltyFactor) const
{
  double sum = 0.0;
  float p = penaltyFactor * fPenaltyValue;
  for (size_t i = 0; i < fNodes.size(); i++) {
    sum += fNodes[i]->GetObjFunction(p, fEndSegWeight);
  }
  return sum / fNodes.size();
}

double pma::Track3D::Optimize(const detinfo::DetectorPropertiesData& detProp,
                              int nNodes,
                              double eps,
                              bool selAllHits,
                              bool setAllNodes,
                              size_t selSegHits,
                              size_t selVtxHits)
{
  if (!fNodes.size()) {
    mf::LogError("pma::Track3D") << "Track not initialized.";
    return 0.0;
  }

  if (!UpdateParams()) {
    mf::LogError("pma::Track3D") << "Track empty.";
    return 1.0e10;
  }
  double g0 = GetObjFunction(), g1 = 0.0;
  if (g0 == 0.0) return g0;

  // set branching flag only at the beginning, optimization is not changing
  // that and new nodes are not branching
  for (auto n : fNodes)
    n->SetVertexToBranching(setAllNodes);

  bool stop = false;
  fMinSegStop = fSegments.size();
  fMaxSegStop = (int)(size() / fMaxSegStopFactor) + 1;
  do {
    if (selSegHits || selVtxHits) SelectRndHits(selSegHits, selVtxHits);

    bool stepDone = true;
    unsigned int stepIter = 0;
    do {
      double gstep = 1.0;
      unsigned int iter = 0;
      while ((gstep > eps) && (iter < 1000)) {
        if ((fNodes.size() < 4) || (iter % 10 == 0))
          MakeProjection();
        else
          MakeFastProjection();

        if (!UpdateParams()) {
          mf::LogError("pma::Track3D") << "Track empty.";
          return 0.0;
        }

        for (auto n : fNodes)
          n->Optimize(fPenaltyValue, fEndSegWeight);

        g1 = g0;
        g0 = GetObjFunction();

        if (g0 == 0.0F) {
          MakeProjection();
          break;
        }
        gstep = fabs(g0 - g1) / g0;
        iter++;
      }

      stepIter++;
      if (fNodes.size() > 2) {
        stepDone = !(CheckEndSegment(pma::Track3D::kEnd) || CheckEndSegment(pma::Track3D::kBegin));
      }
    } while (!stepDone && (stepIter < 5));

    if (selAllHits) {
      selAllHits = false;
      selSegHits = 0;
      selVtxHits = 0;
      if (SelectAllHits()) continue;
    }

    switch (nNodes) {
    case 0: stop = true; break; // just optimize existing vertices

    case -1: // grow and optimize until automatic stop condition
      mf::LogVerbatim("pma::Track3D") << "optimized segments: " << fSegments.size();
      if ((fSegments.size() >= fSegStopValue) || (fSegments.size() >= fMaxSegStop)) { stop = true; }
      else {
        if (!AddNode(detProp)) stop = true;
      }
      break;

    default: // grow and optimize until fixed number of vertices is added
      if (nNodes > 12) {
        if (AddNode(detProp)) {
          MakeProjection();
          nNodes--;
        }
        else {
          mf::LogVerbatim("pma::Track3D") << "stop (3)";
          stop = true;
          break;
        }

        if (AddNode(detProp)) {
          MakeProjection();
          nNodes--;
          if (AddNode(detProp)) nNodes--;
        }
      }
      else if (nNodes > 4) {
        if (AddNode(detProp)) {
          MakeProjection();
          nNodes--;
        }
        else {
          mf::LogVerbatim("pma::Track3D") << "stop (2)";
          stop = true;
          break;
        }

        if (AddNode(detProp)) nNodes--;
      }
      else {
        if (AddNode(detProp)) { nNodes--; }
        else {
          mf::LogVerbatim("pma::Track3D") << "stop (1)";
          stop = true;
          break;
        }
      }
      break;
    }

  } while (!stop);

  MakeProjection();
  return GetObjFunction();
}

bool pma::Track3D::UpdateParamsInTree(bool skipFirst, size_t& depth)
{
  constexpr size_t maxTreeDepth = 100; // really big tree...

  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());

    if (++depth > maxTreeDepth) { mf::LogError("pma::Track3D") << "Tree depth above allowed max."; }
  }

  bool isOK = true;

  while (vtx) {
    auto segThis = NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) { isOK &= seg->Parent()->UpdateParamsInTree(true, depth); }
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  if (!UpdateParams()) {
    mf::LogError("pma::Track3D") << "Track empty.";
    isOK = false;
  }

  --depth;

  return isOK;
}

double pma::Track3D::TuneSinglePass(bool skipFirst)
{
  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  double g = 0.0;
  while (vtx) {
    vtx->Optimize(fPenaltyValue, fEndSegWeight);
    auto segThis = NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) g += seg->Parent()->TuneSinglePass(true);
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  return g + GetObjFunction();
}

pma::Track3D* pma::Track3D::GetNearestTrkInTree(const TVector3& p3d_cm,
                                                double& dist,
                                                bool skipFirst)
{
  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  pma::Track3D* result = this;
  dist = sqrt(Dist2(p3d_cm));

  pma::Track3D* candidate = nullptr;
  while (vtx) {
    auto segThis = NextSegment(vtx);

    double d;
    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) {
        candidate = seg->Parent()->GetNearestTrkInTree(p3d_cm, d, true);
        if (d < dist) {
          dist = d;
          result = candidate;
        }
      }
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  return result;
}

pma::Track3D* pma::Track3D::GetNearestTrkInTree(const TVector2& p2d_cm,
                                                unsigned view,
                                                unsigned int tpc,
                                                unsigned int cryo,
                                                double& dist,
                                                bool skipFirst)
{
  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  pma::Track3D* result = this;
  dist = Dist2(p2d_cm, view, tpc, cryo);

  pma::Track3D* candidate = nullptr;
  while (vtx) {
    auto segThis = NextSegment(vtx);

    double d;
    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) {
        candidate = seg->Parent()->GetNearestTrkInTree(p2d_cm, view, tpc, cryo, d, true);
        if (d < dist) {
          dist = d;
          result = candidate;
        }
      }
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  dist = sqrt(dist);
  return result;
}

void pma::Track3D::ReassignHitsInTree(pma::Track3D* trkRoot)
{
  bool skipFirst;
  if (trkRoot)
    skipFirst = true;
  else {
    trkRoot = this;
    skipFirst = false;
  }

  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  while (vtx) {
    auto segThis = NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) seg->Parent()->ReassignHitsInTree(trkRoot);
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  double d0, dmin;
  pma::Hit3D* hit = nullptr;
  pma::Track3D* nearestTrk = nullptr;
  size_t i = 0;
  while (HasTwoViews(2) && (i < size())) {
    hit = fHits[i];
    d0 = hit->GetDistToProj();

    nearestTrk = GetNearestTrkInTree(hit->Point2D(), hit->View2D(), hit->TPC(), hit->Cryo(), dmin);
    if ((nearestTrk != this) && (dmin < 0.5 * d0)) {
      nearestTrk->push_back(release_at(i));
      mf::LogVerbatim("pma::Track3D") << "*** hit moved to another track ***";
    }
    else
      i++;
  }
  if (!size()) { mf::LogError("pma::Track3D") << "ALL hits moved to other tracks."; }
}

void pma::Track3D::MakeProjectionInTree(bool skipFirst)
{
  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  while (vtx) {
    auto segThis = NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) seg->Parent()->MakeProjectionInTree(true);
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  MakeProjection();
}

void pma::Track3D::SortHitsInTree(bool skipFirst)
{
  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  while (vtx) {
    auto segThis = NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) seg->Parent()->SortHitsInTree(true);
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  SortHits();
}

double pma::Track3D::GetObjFnInTree(bool skipFirst)
{
  pma::Node3D* vtx = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(vtx)) vtx = static_cast<pma::Node3D*>(segThis->Next());
  }

  double g = 0.0;
  while (vtx) {
    auto segThis = NextSegment(vtx);

    for (size_t i = 0; i < vtx->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(vtx->Next(i));
      if (seg != segThis) g += seg->Parent()->GetObjFnInTree(true);
    }

    if (segThis)
      vtx = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  return g + GetObjFunction();
}

double pma::Track3D::TuneFullTree(double eps, double gmax)
{
  if (size_t depth = 1; !UpdateParamsInTree(false, depth)) {
    mf::LogError("pma::Track3D") << "TuneFullTree failed.";
    return -2; // negetive to tag destroyed tree
  }

  double g0 = GetObjFnInTree(), g1 = 0.0;
  if (!std::isfinite(g0)) {
    mf::LogError("pma::Track3D") << "Infinifte value of g.";
    return -2; // negetive to tag destroyed tree
  }
  if (g0 > gmax) {
    mf::LogWarning("pma::Track3D") << "Too high value of g: " << g0;
    return -1; // negetive to tag bad tree
  }
  if (g0 == 0.0) return g0;

  mf::LogVerbatim("pma::Track3D") << "Tune tree, g = " << g0;
  unsigned int stepIter = 0;
  do {
    double gstep = 1.0;
    unsigned int iter = 0;
    while ((gstep > eps) && (iter < 60)) {
      g1 = g0;
      g0 = TuneSinglePass();

      MakeProjectionInTree();

      if (size_t depth = 1; !UpdateParamsInTree(false, depth)) {
        g0 = -2;
        break;
      } // negetive to tag destroyed tree

      if (g0 == 0.0F) break;

      gstep = fabs(g0 - g1) / g0;
      iter++;
    }

    stepIter++;

  } while (stepIter < 5);

  MakeProjectionInTree();
  SortHitsInTree();

  if (g0 >= 0) { mf::LogVerbatim("pma::Track3D") << "  done, g = " << g0; }
  else {
    mf::LogError("pma::Track3D") << "TuneFullTree failed.";
  }
  return g0;
}

void pma::Track3D::ApplyDriftShiftInTree(const detinfo::DetectorClocksData& clockData,
                                         const detinfo::DetectorPropertiesData& detProp,
                                         double dx,
                                         bool skipFirst)
{
  pma::Node3D* node = fNodes.front();

  if (skipFirst) {
    if (auto segThis = NextSegment(node)) node = static_cast<pma::Node3D*>(segThis->Next());
  }

  while (node) {
    node->ApplyDriftShift(dx);

    auto segThis = NextSegment(node);
    for (size_t i = 0; i < node->NextCount(); i++) {
      auto seg = static_cast<pma::Segment3D*>(node->Next(i));
      if (seg != segThis) seg->Parent()->ApplyDriftShiftInTree(clockData, detProp, dx, true);
    }

    if (segThis)
      node = static_cast<pma::Node3D*>(segThis->Next());
    else
      break;
  }

  for (auto h : fHits) {
    h->fPoint3D[0] += dx;
  }

  for (auto p : fAssignedPoints) {
    (*p)[0] += dx;
  }

  // For T0 we need to make sure we use the total shift, not just this current
  // one in case of multiple shifts
  double newdx = fNodes.front()->GetDriftShift();

  // Now convert this newdx into T0 and store in fT0
  SetT0FromDx(clockData, detProp, newdx);
}

void pma::Track3D::SetT0FromDx(detinfo::DetectorClocksData const& clockData,
                               detinfo::DetectorPropertiesData const& detProp,
                               double dx)
{
  auto const* geom = lar::providerFrom<geo::Geometry>();
  const geo::TPCGeo& tpcGeo = geom->TPC(geo::TPCID(fNodes.front()->Cryo(), fNodes.front()->TPC()));

  // GetXTicksCoefficient has a sign that we don't care about. We need to decide
  // the sign for ourselves based on the drift direction of the TPC
  double correctedSign = 0;
  if (tpcGeo.DetectDriftDirection() > 0) {
    if (dx > 0) { correctedSign = 1.0; }
    else {
      correctedSign = -1.0;
    }
  }
  else {
    if (dx > 0) { correctedSign = -1.0; }
    else {
      correctedSign = 1.0;
    }
  }

  // The magnitude of x in ticks is fine
  double dxInTicks = fabs(dx / detProp.GetXTicksCoefficient());
  // Now correct the sign
  dxInTicks = dxInTicks * correctedSign;
  // At this stage, we have xInTicks relative to the trigger time
  double dxInTime = dxInTicks * clockData.TPCClock().TickPeriod();
  // Reconstructed times are relative to the trigger (t=0), so this is our T0
  fT0 = dxInTime;

  mf::LogDebug("pma::Track3D") << dx << ", " << dxInTicks << ", " << correctedSign << ", " << fT0
                               << ", " << tpcGeo.DetectDriftDirection()
                               << " :: " << clockData.TriggerTime() << ", "
                               << clockData.TriggerOffsetTPC() << std::endl;

  // Mark this track as having a measured T0
  fT0Flag = true;
}

void pma::Track3D::DeleteSegments()
{
  for (size_t i = 0; i < fSegments.size(); i++)
    delete fSegments[i];
  fSegments.clear();
}

void pma::Track3D::RebuildSegments()
{
  DeleteSegments();

  fSegments.reserve(fNodes.size() - 1);
  for (size_t i = 1; i < fNodes.size(); i++) {
    fSegments.push_back(new pma::Segment3D(this, fNodes[i - 1], fNodes[i]));
  }
}

void pma::Track3D::CleanupTails()
{
  unsigned int nhits = 0;
  while (!nhits && (fNodes.size() > 2) && !fNodes.front()->IsBranching()) {
    pma::Node3D* vtx = fNodes.front();
    nhits = vtx->NHits();

    pma::Segment3D* seg = NextSegment(vtx);
    nhits += seg->NHits();

    if (!nhits) {
      if (vtx->NextCount() < 2) delete vtx;
      fNodes.erase(fNodes.begin());

      delete seg;
      fSegments.erase(fSegments.begin());

      MakeProjection();
    }
  }

  nhits = 0;
  while (!nhits && (fNodes.size() > 2) && !fNodes.back()->IsBranching()) {
    pma::Node3D* vtx = fNodes.back();
    nhits = vtx->NHits();

    pma::Segment3D* seg = PrevSegment(vtx);
    nhits += seg->NHits();

    if (!nhits) {
      if (vtx->NextCount() < 1) delete fNodes.back();
      fNodes.pop_back();

      delete seg;
      fSegments.pop_back();

      MakeProjection();
    }
  }
}

bool pma::Track3D::ShiftEndsToHits()
{
  pma::Element3D* el;
  pma::Node3D* vtx;

  if (!(fNodes.front()->Prev())) {
    el = GetNearestElement(front()->Point3D());
    vtx = dynamic_cast<pma::Node3D*>(el);
    if (vtx) {
      if (vtx == fNodes.front())
        fNodes.front()->SetPoint3D(front()->Point3D());
      else {
        mf::LogWarning("pma::Track3D") << "First hit is projected to inner node.";
        return false;
      }
    }
    else {
      pma::Segment3D* seg = dynamic_cast<pma::Segment3D*>(el);
      if (seg) {
        if (seg->Prev() == fNodes.front()) {
          double l0 = seg->Length();
          fNodes.front()->SetPoint3D(front()->Point3D());
          if ((seg->Length() < 0.2 * l0) && (fNodes.size() > 2)) {
            pma::Node3D* toRemove = fNodes[1];
            if (toRemove->NextCount() == 1) {
              mf::LogWarning("pma::Track3D")
                << "ShiftEndsToHits(): Short start segment, node removed.";
              fNodes.erase(fNodes.begin() + 1);

              RebuildSegments();
              MakeProjection();

              delete toRemove;
            }
            else {
              mf::LogWarning("pma::Track3D") << "Branching node, not removed.";
              return false;
            }
          }
        }
        else {
          mf::LogWarning("pma::Track3D") << "First hit is projected to inner segment.";
          return false;
        }
      }
    }
  }

  if (!(fNodes.back()->NextCount())) {
    el = GetNearestElement(back()->Point3D());
    vtx = dynamic_cast<pma::Node3D*>(el);
    if (vtx) {
      if (vtx == fNodes.back())
        fNodes.back()->SetPoint3D(back()->Point3D());
      else {
        mf::LogWarning("pma::Track3D") << "First hit is projected to inner node.";
        return false;
      }
    }
    else {
      pma::Segment3D* seg = dynamic_cast<pma::Segment3D*>(el);
      if (seg) {
        if (seg->Next() == fNodes.back()) {
          double l0 = seg->Length();
          fNodes.back()->SetPoint3D(back()->Point3D());
          if ((seg->Length() < 0.2 * l0) && (fNodes.size() > 2)) {
            size_t idx = fNodes.size() - 2;
            pma::Node3D* toRemove = fNodes[idx];
            if (toRemove->NextCount() == 1) {
              mf::LogWarning("pma::Track3D")
                << "ShiftEndsToHits(): Short end segment, node removed.";
              fNodes.erase(fNodes.begin() + idx);

              RebuildSegments();
              MakeProjection();

              delete toRemove;
            }
            else {
              mf::LogWarning("pma::Track3D") << "Branching node, not removed.";
              return false;
            }
          }
        }
        else {
          mf::LogWarning("pma::Track3D") << "First hit is projected to inner segment.";
          return false;
        }
      }
    }
  }

  return true;
}

double pma::Track3D::Dist2(const TVector2& p2d,
                           unsigned int view,
                           unsigned int tpc,
                           unsigned int cryo) const
{
  double dist, min_dist = 1.0e12;

  int t = (int)tpc, c = (int)cryo;
  auto n0 = fNodes.front();
  if ((n0->TPC() == t) && (n0->Cryo() == c)) {
    dist = n0->GetDistance2To(p2d, view);
    if (dist < min_dist) min_dist = dist;
  }
  auto n1 = fNodes.back();
  if ((n1->TPC() == t) && (n1->Cryo() == c)) {
    dist = n1->GetDistance2To(p2d, view);
    if (dist < min_dist) min_dist = dist;
  }

  for (auto s : fSegments)
    if ((s->TPC() == t) && (s->Cryo() == c)) {
      dist = s->GetDistance2To(p2d, view);
      if (dist < min_dist) min_dist = dist;
    }
  return min_dist;
}

double pma::Track3D::Dist2(const TVector3& p3d) const
{
  using namespace ranges;
  auto to_distance2 = [&p3d](auto seg) { return seg->GetDistance2To(p3d); };
  return min(fSegments | views::transform(to_distance2));
}

pma::Element3D* pma::Track3D::GetNearestElement(const TVector2& p2d,
                                                unsigned int view,
                                                int tpc,
                                                bool skipFrontVtx,
                                                bool skipBackVtx) const
{
  if (fSegments.front()->TPC() < 0) skipFrontVtx = false;
  if (fSegments.back()->TPC() < 0) skipBackVtx = false;

  if (skipFrontVtx && skipBackVtx && (fSegments.size() == 1))
    return fSegments.front(); // no need for searching...

  size_t v0 = 0, v1 = fNodes.size();
  if (skipFrontVtx) v0 = 1;
  if (skipBackVtx) --v1;

  pma::Element3D* pe_min = nullptr;
  auto min_dist = std::numeric_limits<double>::max();
  for (size_t i = v0; i < v1; i++)
    if (fNodes[i]->TPC() == tpc) {
      double dist = fNodes[i]->GetDistance2To(p2d, view);
      if (dist < min_dist) {
        min_dist = dist;
        pe_min = fNodes[i];
      }
    }
  for (auto segment : fSegments)
    if (segment->TPC() == tpc) {
      if (segment->TPC() < 0) continue; // segment between TPC's

      double const dist = segment->GetDistance2To(p2d, view);
      if (dist < min_dist) {
        min_dist = dist;
        pe_min = segment;
      }
    }
  if (!pe_min) throw cet::exception("pma::Track3D") << "Nearest element not found." << std::endl;
  return pe_min;
}

pma::Element3D* pma::Track3D::GetNearestElement(const TVector3& p3d) const
{
  pma::Element3D* pe_min = fNodes.front();
  double dist, min_dist = pe_min->GetDistance2To(p3d);
  for (size_t i = 1; i < fNodes.size(); i++) {
    dist = fNodes[i]->GetDistance2To(p3d);
    if (dist < min_dist) {
      min_dist = dist;
      pe_min = fNodes[i];
    }
  }
  for (size_t i = 0; i < fSegments.size(); i++) {
    dist = fSegments[i]->GetDistance2To(p3d);
    if (dist < min_dist) {
      min_dist = dist;
      pe_min = fSegments[i];
    }
  }
  return pe_min;
}

bool pma::Track3D::GetUnconstrainedProj3D(detinfo::DetectorPropertiesData const& detProp,
                                          art::Ptr<recob::Hit> hit,
                                          TVector3& p3d,
                                          double& dist2) const
{
  TVector2 p2d = pma::WireDriftToCm(detProp,
                                    hit->WireID().Wire,
                                    hit->PeakTime(),
                                    hit->WireID().Plane,
                                    hit->WireID().TPC,
                                    hit->WireID().Cryostat);

  pma::Segment3D* seg = nullptr;
  double d2, min_d2 = 1.0e100;
  int tpc = (int)hit->WireID().TPC;
  int view = hit->WireID().Plane;
  for (size_t i = 0; i < fSegments.size(); i++) {
    if (fSegments[i]->TPC() != tpc) continue;

    d2 = fSegments[i]->GetDistance2To(p2d, view);
    if (d2 < min_d2) {
      min_d2 = d2;
      seg = fSegments[i];
    }
  }
  if (seg) {
    p3d = seg->GetUnconstrainedProj3D(p2d, view);
    dist2 = min_d2;

    pma::Node3D* prev = static_cast<pma::Node3D*>(seg->Prev());
    return prev->SameTPC(p3d); // 3D can be beyond the segment endpoints => in other TPC
  }

  return false;
}

void pma::Track3D::SortHits()
{
  std::vector<pma::Hit3D*> hits_tmp;
  hits_tmp.reserve(size());

  pma::Node3D* vtx = fNodes.front();
  pma::Segment3D* seg = NextSegment(vtx);
  while (vtx) {
    vtx->SortHits();
    for (size_t i = 0; i < vtx->NHits(); i++) {
      pma::Hit3D* h3d = &(vtx->Hit(i));
      if (h3d->fParent == this) hits_tmp.push_back(h3d);
    }

    if (seg) {
      seg->SortHits();
      for (size_t i = 0; i < seg->NHits(); i++) {
        pma::Hit3D* h3d = &(seg->Hit(i));
        if (h3d->fParent == this) hits_tmp.push_back(h3d);
      }

      vtx = static_cast<pma::Node3D*>(seg->Next());
      seg = NextSegment(vtx);
    }
    else
      break;
  }

  if (size() == hits_tmp.size()) {
    for (size_t i = 0; i < size(); i++) {
      fHits[i] = hits_tmp[i];
    }
  }
  else
    mf::LogError("pma::Track3D") << "Hit sorting problem.";
}

unsigned int pma::Track3D::DisableSingleViewEnds()
{
  SortHits();

  unsigned int nDisabled = 0;

  std::array<int, 4> hasHits{{}};

  pma::Hit3D* nextHit = nullptr;
  int hitIndex = -1;

  bool rebuild = false, stop = false;
  int nViews = 0;

  do {
    pma::Node3D* vtx = fNodes.front();
    pma::Segment3D* seg = NextSegment(vtx);
    if (!seg) break;

    if (vtx->NPoints() + seg->NPoints() > 0) hasHits[3] = 1;

    for (size_t i = 0; i < vtx->NHits(); i++) {
      hitIndex = index_of(&(vtx->Hit(i)));
      if ((hitIndex >= 0) && (hitIndex + 1 < (int)size()))
        nextHit = fHits[hitIndex + 1];
      else
        nextHit = 0;

      if (vtx->Hit(i).IsEnabled()) hasHits[vtx->Hit(i).View2D()] = 1;
      if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
      nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
      if (nViews < 2) {
        if (vtx->Hit(i).IsEnabled()) {
          vtx->Hit(i).SetEnabled(false);
          nDisabled++;
        }
      }
    }
    for (size_t i = 0; i < seg->NHits(); i++) {
      hitIndex = index_of(&(seg->Hit(i)));
      if ((hitIndex >= 0) && (hitIndex + 1 < (int)size()))
        nextHit = fHits[hitIndex + 1];
      else
        nextHit = 0;

      if (seg->Hit(i).IsEnabled()) hasHits[seg->Hit(i).View2D()] = 1;
      if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
      nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
      if (nViews < 2) {
        if (seg->Hit(i).IsEnabled()) {
          seg->Hit(i).SetEnabled(false);
          nDisabled++;
        }
      }
    }

    if (fNodes.size() < 3) break;

    nViews = hasHits[1] + hasHits[2] + hasHits[3];
    if (hasHits[0] || (nViews > 1))
      stop = true;
    else if (!fNodes.front()->IsBranching()) {
      pma::Node3D* vtx_front = fNodes.front();
      fNodes.erase(fNodes.begin());
      delete vtx_front;
      rebuild = true;
    }

  } while (!stop);

  stop = false;
  nViews = 0;
  hasHits.fill(0);

  do {
    pma::Node3D* vtx = fNodes.back();
    pma::Segment3D* seg = PrevSegment(vtx);
    if (!seg) break;

    if (vtx->NPoints() || seg->NPoints()) hasHits[3] = 1;

    for (int i = vtx->NHits() - 1; i >= 0; i--) {
      hitIndex = index_of(&(vtx->Hit(i)));
      if ((hitIndex >= 0) && (hitIndex - 1 >= 0))
        nextHit = fHits[hitIndex - 1];
      else
        nextHit = 0;

      if (vtx->Hit(i).IsEnabled()) hasHits[vtx->Hit(i).View2D()] = 1;
      if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
      nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
      if (nViews < 2) {
        if (vtx->Hit(i).IsEnabled()) {
          vtx->Hit(i).SetEnabled(false);
          nDisabled++;
        }
      }
    }
    for (int i = seg->NHits() - 1; i >= 0; i--) {
      hitIndex = index_of(&(seg->Hit(i)));
      if ((hitIndex >= 0) && (hitIndex - 1 >= 0))
        nextHit = fHits[hitIndex - 1];
      else
        nextHit = 0;

      if (seg->Hit(i).IsEnabled()) hasHits[seg->Hit(i).View2D()] = 1;
      if (nextHit && nextHit->IsEnabled()) hasHits[nextHit->View2D()] = 1;
      nViews = hasHits[0] + hasHits[1] + hasHits[2] + hasHits[3];
      if (nViews < 2) {
        if (seg->Hit(i).IsEnabled()) {
          seg->Hit(i).SetEnabled(false);
          nDisabled++;
        }
      }
    }

    if (fNodes.size() < 3) break;

    nViews = hasHits[1] + hasHits[2] + hasHits[3];
    if (hasHits[0] || (nViews > 1))
      stop = true;
    else if (!fNodes.back()->IsBranching()) {
      pma::Node3D* vtx_back = fNodes.back();
      fNodes.pop_back();
      delete vtx_back;
      rebuild = true;
    }

  } while (!stop);

  if (rebuild) {
    RebuildSegments();
    MakeProjection();
  }

  return nDisabled;
}

bool pma::Track3D::SelectHits(float fraction)
{
  if (fraction < 0.0F) fraction = 0.0F;
  if (fraction > 1.0F) fraction = 1.0F;
  if (fraction < 1.0F) std::sort(fHits.begin(), fHits.end(), pma::bTrajectory3DDistLess());

  size_t nHitsColl = (size_t)(fraction * NHits(geo::kZ));
  size_t nHitsInd2 = (size_t)(fraction * NHits(geo::kV));
  size_t nHitsInd1 = (size_t)(fraction * NHits(geo::kU));
  size_t coll = 0, ind2 = 0, ind1 = 0;

  bool b, changed = false;
  for (auto h : fHits) {
    b = h->IsEnabled();
    if (fraction < 1.0F) {
      h->SetEnabled(false);
      switch (h->View2D()) {
      case geo::kZ:
        if (coll++ < nHitsColl) h->SetEnabled(true);
        break;
      case geo::kV:
        if (ind2++ < nHitsInd2) h->SetEnabled(true);
        break;
      case geo::kU:
        if (ind1++ < nHitsInd1) h->SetEnabled(true);
        break;
      }
    }
    else
      h->SetEnabled(true);

    changed |= (b != h->IsEnabled());
  }

  if (fraction < 1.0F) {
    FirstElement()->SelectAllHits();
    LastElement()->SelectAllHits();
  }
  return changed;
}

bool pma::Track3D::SelectRndHits(size_t segmax, size_t vtxmax)
{
  bool changed = false;
  for (auto n : fNodes)
    changed |= n->SelectRndHits(vtxmax);
  for (auto s : fSegments)
    changed |= s->SelectRndHits(segmax);
  return changed;
}

bool pma::Track3D::SelectAllHits()
{
  bool changed = false;
  for (auto h : fHits) {
    changed |= !(h->IsEnabled());
    h->SetEnabled(true);
  }
  return changed;
}

void pma::Track3D::MakeProjection()
{
  for (auto n : fNodes)
    n->ClearAssigned(this);
  for (auto s : fSegments)
    s->ClearAssigned(this);

  pma::Element3D* pe = nullptr;

  bool skipFrontVtx = false, skipBackVtx = false;

  // skip outermost vertices if not branching
  if (!(fNodes.front()->IsFrozen()) && !(fNodes.front()->Prev()) &&
      (fNodes.front()->NextCount() == 1)) {
    skipFrontVtx = true;
  }
  if (!(fNodes.front()->IsFrozen()) && (fNodes.back()->NextCount() == 0)) { skipBackVtx = true; }

  for (auto h : fHits) // assign hits to nodes/segments
  {
    pe = GetNearestElement(h->Point2D(), h->View2D(), h->TPC(), skipFrontVtx, skipBackVtx);
    pe->AddHit(h);
  }

  for (auto p : fAssignedPoints) // assign ref points to nodes/segments
  {
    pe = GetNearestElement(*p);
    pe->AddPoint(p);
  }

  for (auto n : fNodes)
    n->UpdateHitParams();
  for (auto s : fSegments)
    s->UpdateHitParams();
}

void pma::Track3D::MakeFastProjection()
{
  std::vector<std::pair<pma::Hit3D*, pma::Element3D*>> assignments;
  assignments.reserve(fHits.size());

  for (auto hi : fHits) {
    pma::Element3D* pe = nullptr;

    for (auto s : fSegments) {
      for (size_t j = 0; j < s->NHits(); ++j)
        if (hi == s->Hits()[j]) // look at next/prev vtx,seg,vtx
        {
          pe = s;
          double min_d2 = s->GetDistance2To(hi->Point2D(), hi->View2D());
          int const tpc = hi->TPC();

          pma::Node3D* nnext = static_cast<pma::Node3D*>(s->Next());
          if (nnext->TPC() == tpc) {
            double const d2 = nnext->GetDistance2To(hi->Point2D(), hi->View2D());
            if (d2 < min_d2) {
              min_d2 = d2;
              pe = nnext;
            }

            pma::Segment3D* snext = NextSegment(nnext);
            if (snext && (snext->TPC() == tpc)) {
              double const d2 = snext->GetDistance2To(hi->Point2D(), hi->View2D());
              if (d2 < min_d2) {
                min_d2 = d2;
                pe = snext;
              }

              nnext = static_cast<pma::Node3D*>(snext->Next());
              if (nnext->TPC() == tpc) {
                double const d2 = nnext->GetDistance2To(hi->Point2D(), hi->View2D());
                if (d2 < min_d2) {
                  min_d2 = d2;
                  pe = nnext;
                }
              }
            }
          }

          pma::Node3D* nprev = static_cast<pma::Node3D*>(s->Prev());
          if (nprev->TPC() == tpc) {
            double const d2 = nprev->GetDistance2To(hi->Point2D(), hi->View2D());
            if (d2 < min_d2) {
              min_d2 = d2;
              pe = nprev;
            }

            pma::Segment3D* sprev = PrevSegment(nprev);
            if (sprev && (sprev->TPC() == tpc)) {
              double const d2 = sprev->GetDistance2To(hi->Point2D(), hi->View2D());
              if (d2 < min_d2) {
                min_d2 = d2;
                pe = sprev;
              }

              nprev = static_cast<pma::Node3D*>(sprev->Prev());
              if (nprev->TPC() == tpc) {
                double const d2 = nprev->GetDistance2To(hi->Point2D(), hi->View2D());
                if (d2 < min_d2) {
                  min_d2 = d2;
                  pe = nprev;
                }
              }
            }
          }

          s->RemoveHitAt(j);
          break;
        }
      if (pe) break;
    }

    if (!pe)
      for (auto n : fNodes) {
        for (size_t j = 0; j < n->NHits(); ++j)
          if (hi == n->Hits()[j]) // look at next/prev seg,vtx,seg
          {
            pe = n;
            double d2, min_d2 = n->GetDistance2To(hi->Point2D(), hi->View2D());
            int tpc = hi->TPC();

            pma::Segment3D* snext = NextSegment(n);
            if (snext && (snext->TPC() == tpc)) {
              d2 = snext->GetDistance2To(hi->Point2D(), hi->View2D());
              if (d2 < min_d2) {
                min_d2 = d2;
                pe = snext;
              }

              pma::Node3D* nnext = static_cast<pma::Node3D*>(snext->Next());
              if (nnext->TPC() == tpc) {
                d2 = nnext->GetDistance2To(hi->Point2D(), hi->View2D());
                if (d2 < min_d2) {
                  min_d2 = d2;
                  pe = nnext;
                }

                snext = NextSegment(nnext);
                if (snext && (snext->TPC() == tpc)) {
                  d2 = snext->GetDistance2To(hi->Point2D(), hi->View2D());
                  if (d2 < min_d2) {
                    min_d2 = d2;
                    pe = snext;
                  }
                }
              }
            }

            pma::Segment3D* sprev = PrevSegment(n);
            if (sprev && (sprev->TPC() == tpc)) {
              d2 = sprev->GetDistance2To(hi->Point2D(), hi->View2D());
              if (d2 < min_d2) {
                min_d2 = d2;
                pe = sprev;
              }

              pma::Node3D* nprev = static_cast<pma::Node3D*>(sprev->Prev());
              if (nprev->TPC() == tpc) {
                d2 = nprev->GetDistance2To(hi->Point2D(), hi->View2D());
                if (d2 < min_d2) {
                  min_d2 = d2;
                  pe = nprev;
                }

                sprev = PrevSegment(nprev);
                if (sprev && (sprev->TPC() == tpc)) {
                  d2 = sprev->GetDistance2To(hi->Point2D(), hi->View2D());
                  if (d2 < min_d2) {
                    min_d2 = d2;
                    pe = sprev;
                  }
                }
              }
            }

            n->RemoveHitAt(j);
            break;
          }
        if (pe) break;
      }

    if (pe)
      assignments.emplace_back(hi, pe);
    else
      mf::LogWarning("pma::Track3D") << "Hit was not assigned to any element.";
  }

  for (auto const& a : assignments)
    a.second->AddHit(a.first);

  for (auto n : fNodes)
    n->UpdateHitParams();
  for (auto s : fSegments)
    s->UpdateHitParams();
}

void pma::Track3D::UpdateProjection()
{
  for (auto n : fNodes)
    n->UpdateProjection();
  for (auto s : fSegments)
    s->UpdateProjection();
}

double pma::Track3D::AverageDist2() const
{
  double sum = 0.0;
  unsigned int count = 0;

  for (auto n : fNodes) {
    sum += n->SumDist2();
    count += n->NEnabledHits();
  }

  for (auto s : fSegments) {
    sum += s->SumDist2();
    count += s->NEnabledHits();
  }

  if (count) { return sum / count; }
  else {
    mf::LogError("pma::Track3D") << "0 enabled hits in AverageDist2 calculation.";
    return 0;
  }
}

bool pma::Track3D::UpdateParams()
{
  size_t n = size();
  if (n == 0) {
    fPenaltyValue = 1;
    fSegStopValue = 1;
    return false;
  }

  float nCubeRoot = pow((double)n, 1.0 / 3.0);
  float avgDist2Root = sqrt(AverageDist2());
  if (avgDist2Root > 0) {
    fPenaltyValue = fPenaltyFactor * pow((double)fSegments.size(), 1.8) * avgDist2Root /
                    (fHitsRadius * nCubeRoot);

    fSegStopValue = (int)(fSegStopFactor * nCubeRoot * fHitsRadius / avgDist2Root);
    if (fSegStopValue < fMinSegStop) fSegStopValue = fMinSegStop;
    return true;
  }
  else {
    fPenaltyValue = 1;
    fSegStopValue = 1;
    return false;
  }
}

bool pma::Track3D::SwapVertices(size_t v0, size_t v1)
{
  if (v0 == v1) return false;

  if (v0 > v1) {
    size_t vx = v0;
    v0 = v1;
    v1 = vx;
  }

  pma::Node3D* vtmp;
  if (v1 - v0 == 1) {
    pma::Segment3D* midSeg = NextSegment(fNodes[v0]);
    pma::Segment3D* prevSeg = PrevSegment(fNodes[v0]);
    pma::Segment3D* nextSeg = NextSegment(fNodes[v1]);

    fNodes[v1]->RemoveNext(nextSeg);
    midSeg->Disconnect();

    vtmp = fNodes[v0];
    fNodes[v0] = fNodes[v1];
    fNodes[v1] = vtmp;

    if (prevSeg) prevSeg->AddNext(fNodes[v0]);
    fNodes[v0]->AddNext(midSeg);
    midSeg->AddNext(fNodes[v1]);
    if (nextSeg) fNodes[v1]->AddNext(nextSeg);

    return false;
  }
  else {
    vtmp = fNodes[v0];
    fNodes[v0] = fNodes[v1];
    fNodes[v1] = vtmp;
    return true;
  }
}

bool pma::Track3D::CheckEndSegment(pma::Track3D::ETrackEnd endCode)
{
  unsigned int v1, v2;
  switch (endCode) {
  case pma::Track3D::kBegin:
    if (fSegments.front()->IsFrozen()) return false;
    if (fNodes.front()->NextCount() > 1) return false;
    v1 = 0;
    v2 = 1;
    break;
  case pma::Track3D::kEnd:
    if (fSegments.back()->IsFrozen()) return false;
    if (fNodes.back()->NextCount() > 1) return false;
    v1 = fNodes.size() - 1;
    v2 = fNodes.size() - 2;
    break;
  default: return false;
  }
  if (fNodes[v1]->TPC() != fNodes[v2]->TPC()) return false;

  double g1, g0 = GetObjFunction();

  if (SwapVertices(v1, v2)) RebuildSegments();
  MakeProjection();
  g1 = GetObjFunction();

  if (g1 >= g0) {
    if (SwapVertices(v1, v2)) RebuildSegments();
    MakeProjection();
    return false;
  }
  else
    return true;
}

void pma::Track3D::UpdateHitsRadius()
{
  std::vector<pma::Hit3D*> hitsColl, hitsInd1, hitsInd2;
  for (auto hit : fHits) {
    switch (hit->View2D()) {
    case geo::kZ: hitsColl.push_back(hit); break;
    case geo::kV: hitsInd2.push_back(hit); break;
    case geo::kU: hitsInd1.push_back(hit); break;
    }
  }
  fHitsRadius = std::max({pma::GetHitsRadius2D(hitsColl, true),
                          pma::GetHitsRadius2D(hitsInd2, true),
                          pma::GetHitsRadius2D(hitsInd1, true)});
}
