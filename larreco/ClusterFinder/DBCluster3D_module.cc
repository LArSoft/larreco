////////////////////////////////////////////////////////////////////////
// Class:       DBCluster3D
// Plugin Type: producer (art v2_11_02)
// File:        DBCluster3D_module.cc
//
// Generated at Sat Jun 16 07:43:32 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/DBScan3DAlg.h"

#include <memory>

namespace cluster {
  class DBCluster3D;
}

class cluster::DBCluster3D : public art::EDProducer {

  using Point_t = recob::tracking::Point_t;
  using Vector_t = recob::tracking::Vector_t;

public:
  explicit DBCluster3D(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DBCluster3D(DBCluster3D const&) = delete;
  DBCluster3D(DBCluster3D&&) = delete;
  DBCluster3D& operator=(DBCluster3D const&) = delete;
  DBCluster3D& operator=(DBCluster3D&&) = delete;

private:
  // Required functions.
  void produce(art::Event& e) override;

  const art::InputTag fHitModuleLabel;
  const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fSPHitAssnLabel;

  DBScan3DAlg fDBScan;

  geo::GeometryCore const* fGeom;

  double tickToDist;
  double fMinHitDis;
};

cluster::DBCluster3D::DBCluster3D(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fHitModuleLabel(p.get<art::InputTag>("HitModuleLabel"))
  , fSpacePointModuleLabel(p.get<art::InputTag>("SpacePointModuleLabel"))
  , fSPHitAssnLabel(p.get<art::InputTag>("SPHitAssnLabel"))
  , fDBScan(p.get<fhicl::ParameterSet>("DBScan3DAlg"))
  , fMinHitDis(p.get<double>("MinHitDis"))
{
  produces<std::vector<recob::Slice>>();
  produces<art::Assns<recob::Slice, recob::Hit>>();
  produces<art::Assns<recob::Slice, recob::SpacePoint>>();

  fGeom = art::ServiceHandle<geo::Geometry const>().get();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const det_prop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clock_data);

  tickToDist = det_prop.DriftVelocity(det_prop.Efield(), det_prop.Temperature());
  tickToDist *= 1.e-3 * sampling_rate(clock_data); // 1e-3 is conversion of 1/us to 1/ns
  fMinHitDis *= fMinHitDis;
}

void cluster::DBCluster3D::produce(art::Event& evt)
{
  auto scol = std::make_unique<std::vector<recob::Slice>>();
  auto slc_hit_assn = std::make_unique<art::Assns<recob::Slice, recob::Hit>>();
  auto slc_sps_assn = std::make_unique<art::Assns<recob::Slice, recob::SpacePoint>>();

  auto hitsHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);
  auto spsHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointModuleLabel);

  // all hits in the collection
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hitsHandle);

  // all space points in the collection
  std::vector<art::Ptr<recob::SpacePoint>> sps;
  art::fill_ptr_vector(sps, spsHandle);

  art::FindManyP<recob::SpacePoint> spFromHit(hitsHandle, evt, fSPHitAssnLabel);
  if (!spFromHit.isValid()) {
    std::cout << "spFromHit is invalid\n";
    return;
  }
  // Find the first Hit - SpacePoint assn and check consistency on the first event
  static bool first = true;
  if (first) {
    bool success = false;
    bool foundsps = false;
    for (auto& hit : hits) {
      auto& sps = spFromHit.at(hit.key());
      if (sps.empty()) continue;
      success = (sps[0].id() == spsHandle.id());
      foundsps = true;
      break;
    } // hit
    if ((!success) && foundsps)
      throw cet::exception("DBCluster3D")
        << "HitModuleLabel, SpacePointModuleLabel and SPHitAssnLabel are inconsistent\n";
    first = false;
  } // first

  art::FindManyP<recob::Hit> hitFromSp(spsHandle, evt, fSPHitAssnLabel);
  if (!hitFromSp.isValid()) {
    std::cout << "hitFromSp is invalid\n";
    return;
  }
  fDBScan.init(sps, hitFromSp);
  fDBScan.dbscan();

  //Find number of slices
  int maxid = INT_MIN;
  for (size_t i = 0; i < fDBScan.points.size(); ++i) {
    if (fDBScan.points[i].cluster_id > maxid) maxid = fDBScan.points[i].cluster_id;
  }
  size_t nslc = 0;
  if (maxid >= 0) nslc = maxid + 1;

  //Save hits associated with each slice
  std::vector<std::vector<art::Ptr<recob::Hit>>> slcHits(nslc);
  //Save hits on each PlaneID with pfparticle index
  std::map<geo::PlaneID, std::vector<std::pair<art::Ptr<recob::Hit>, unsigned int>>> hitmap;
  for (auto& hit : hits) {
    auto& sps = spFromHit.at(hit.key());
    if (sps.size()) { //found associated space point
      if (fDBScan.points[sps[0].key()].cluster_id >= 0) {
        slcHits[fDBScan.points[sps[0].key()].cluster_id].push_back(hit);
        hitmap[geo::PlaneID(hit->WireID())].push_back(
          std::make_pair(hit, fDBScan.points[sps[0].key()].cluster_id));
      }
    } // sps.size()
  }   // hit

  //Save hits not associated with any spacepoints
  for (auto& hit : hits) {
    bool found = false;
    for (size_t i = 0; i < slcHits.size(); ++i) {
      if (std::find(slcHits[i].begin(), slcHits[i].end(), hit) != slcHits[i].end()) {
        found = true;
        break;
      }
    }
    if (!found) {
      double wirePitch = fGeom->WirePitch(hit->WireID());
      double UnitsPerTick = tickToDist / wirePitch;
      double x0 = hit->WireID().Wire;
      double y0 = hit->PeakTime() * UnitsPerTick;
      double mindis = DBL_MAX;
      unsigned slcIndex = UINT_MAX;
      for (auto& hit2 : hitmap[geo::PlaneID(hit->WireID())]) {
        double dx = hit2.first->WireID().Wire - x0;
        double dy = hit2.first->PeakTime() * UnitsPerTick - y0;
        double dis = dx * dx + dy * dy;
        if (dis < mindis) {
          mindis = dis;
          slcIndex = hit2.second;
        }
      }
      if (slcIndex != UINT_MAX && mindis < fMinHitDis) { slcHits[slcIndex].push_back(hit); }
    }
  }

  //Save spacepoints for each slice
  std::vector<std::vector<art::Ptr<recob::SpacePoint>>> sps_in_slc(nslc);
  for (size_t i = 0; i < fDBScan.points.size(); ++i) {
    if (fDBScan.points[i].cluster_id >= 0) {
      sps_in_slc[fDBScan.points[i].cluster_id].push_back(sps[i]);
    }
  } // i

  // calculate the properties of the slice
  for (size_t isl = 0; isl < nslc; ++isl) {
    double sum = sps_in_slc[isl].size();
    // find the center
    double center[3] = {0.};
    // TODO: calculate the charge using recob::PointCharge instead of just counting hits
    float charge = slcHits[isl].size();
    for (auto& spt : sps_in_slc[isl]) {
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        center[xyz] += spt->XYZ()[xyz];
    } // spt
    for (unsigned short xyz = 0; xyz < 3; ++xyz)
      center[xyz] /= sum;
    // do a linear fit
    double sumx = 0, sumy = 0., sumz = 0., sumx2 = 0, sumy2 = 0, sumz2 = 0;
    double sumxy = 0, sumxz = 0;
    for (auto& spt : sps_in_slc[isl]) {
      double xx = spt->XYZ()[0] - center[0];
      double yy = spt->XYZ()[1] - center[1];
      double zz = spt->XYZ()[2] - center[2];
      sumx += xx;
      sumy += yy;
      sumz += zz;
      sumx2 += xx * xx;
      sumy2 += yy * yy;
      sumz2 += zz * zz;
      sumxy += xx * yy;
      sumxz += xx * zz;
    } // spt
    double delta = sum * sumx2 - sumx * sumx;
    if (delta <= 0) continue;
    // calculate the slopes
    double dydx = (sumxy * sum - sumx * sumy) / delta;
    double dzdx = (sumxz * sum - sumx * sumz) / delta;
    // and convert to direction vector
    double direction[3];
    double norm = std::sqrt(1 + dydx * dydx + dzdx * dzdx);
    direction[0] = 1 / norm;
    direction[1] = dydx / norm;
    direction[2] = dzdx / norm;
    // Find the points that are farthest away from the center along the slice axis
    unsigned int imax = 0, imin = 0;
    float minAlong = 1E6, maxAlong = -1E6;
    double tmpVec[3];
    for (unsigned int ipt = 0; ipt < sps_in_slc[isl].size(); ++ipt) {
      auto& spt = sps_in_slc[isl][ipt];
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        tmpVec[xyz] = spt->XYZ()[xyz] - center[xyz];
      double dotp = 0;
      for (unsigned short xyz = 0; xyz < 3; ++xyz)
        dotp += tmpVec[xyz] * direction[xyz];
      if (dotp < minAlong) {
        minAlong = dotp;
        imin = ipt;
      }
      if (dotp > maxAlong) {
        maxAlong = dotp;
        imax = ipt;
      }
    } // spt
    // Find the aspect ratio
    float cnt = 0.;
    float aspectRatio = 0;
    double arg = sum * sumy2 - sumy * sumy;
    if (arg > 0) {
      aspectRatio += std::abs(sum * sumxy - sumx * sumy) / sqrt(delta * arg);
      ++cnt;
    }
    arg = sum * sumz2 - sumz * sumz;
    if (arg > 0) {
      aspectRatio += std::abs(sum * sumxz - sumx * sumz) / sqrt(delta * arg);
      ++cnt;
    }
    if (cnt > 1) aspectRatio /= cnt;
    int id = isl + 1;
    Point_t ctr(center[0], center[1], center[2]);
    Vector_t dir(direction[0], direction[1], direction[2]);
    auto pos0 = sps_in_slc[isl][imin]->XYZ();
    Point_t ep0(pos0[0], pos0[1], pos0[2]);
    auto pos1 = sps_in_slc[isl][imax]->XYZ();
    Point_t ep1(pos1[0], pos1[1], pos1[2]);
    scol->emplace_back(id, ctr, dir, ep0, ep1, aspectRatio, charge);
    util::CreateAssn(evt, *scol, slcHits[isl], *slc_hit_assn);
    util::CreateAssn(evt, *scol, sps_in_slc[isl], *slc_sps_assn);
  } // isl

  evt.put(std::move(scol));
  evt.put(std::move(slc_hit_assn));
  evt.put(std::move(slc_sps_assn));
}

DEFINE_ART_MODULE(cluster::DBCluster3D)
