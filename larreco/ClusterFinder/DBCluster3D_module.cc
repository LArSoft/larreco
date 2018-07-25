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
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/DBScan3DAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <memory>

namespace cluster {
  class DBCluster3D;
}


class cluster::DBCluster3D : public art::EDProducer {
public:
  explicit DBCluster3D(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DBCluster3D(DBCluster3D const &) = delete;
  DBCluster3D(DBCluster3D &&) = delete;
  DBCluster3D & operator = (DBCluster3D const &) = delete;
  DBCluster3D & operator = (DBCluster3D &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  const art::InputTag fHitModuleLabel;
  const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fSPHitAssnLabel;

  DBScan3DAlg fDBScan;

  geo::GeometryCore const* fGeom;
  detinfo::DetectorProperties const* fDetProp;

  double tickToDist;
  double fMinHitDis;

};


cluster::DBCluster3D::DBCluster3D(fhicl::ParameterSet const & p)
  : fHitModuleLabel(p.get< art::InputTag >("HitModuleLabel"))
  , fSpacePointModuleLabel(p.get< art::InputTag >("SpacePointModuleLabel"))
  , fSPHitAssnLabel(p.get< art::InputTag >("SPHitAssnLabel"))
  , fDBScan(p.get< fhicl::ParameterSet >("DBScan3DAlg"))
  , fMinHitDis(p.get< double >("MinHitDis"))
{
  produces< std::vector<recob::Slice> >();
  produces< art::Assns<recob::Slice, recob::Hit> >();

  fGeom = &*(art::ServiceHandle<geo::Geometry>());
  fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  tickToDist = fDetProp->DriftVelocity(fDetProp->Efield(),fDetProp->Temperature());
  tickToDist *= 1.e-3 * fDetProp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns     
}

void cluster::DBCluster3D::produce(art::Event & evt)
{

  std::unique_ptr<std::vector<recob::Slice>> 
    slcCol(new std::vector<recob::Slice>);

  std::unique_ptr<art::Assns<recob::Slice, recob::Hit>> 
    slc_hit_assn(new art::Assns<recob::Slice, recob::Hit>);

  auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
  auto spsHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpacePointModuleLabel);

  // all hits in the collection
  std::vector< art::Ptr<recob::Hit> > hits;
  art::fill_ptr_vector(hits, hitsHandle);

  // all space points in the collection
  std::vector< art::Ptr<recob::SpacePoint> > sps;
  art::fill_ptr_vector(sps, spsHandle);

  art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSPHitAssnLabel);

  fDBScan.init(sps);
  fDBScan.dbscan();

  //Find number of slices
  int maxid = INT_MIN;
  for (size_t i = 0; i<fDBScan.points.size(); ++i){
//    std::cout<<"Space point index "<<i<<" "<<fDBScan.points[i].sp.key()<<" "<<fDBScan.points[i].cluster_id<<std::endl;
    if (fDBScan.points[i].cluster_id > maxid) maxid = fDBScan.points[i].cluster_id;
  }
  size_t nslc = 0;
  if (maxid>=0) nslc = maxid + 1;

  //Save hits associated with each slice
  std::vector<std::vector<art::Ptr<recob::Hit>>> slcHits(nslc);
  //Save hits on each PlaneID with pfparticle index
  std::map<geo::PlaneID, std::vector<std::pair<art::Ptr<recob::Hit>, unsigned int>>> hitmap;
  for (auto &hit : hits){
    auto &sps = spFromHit.at(hit.key());
    if (sps.size()){//found associated space point
      if (fDBScan.points[sps[0].key()].cluster_id>=0){
        slcHits[fDBScan.points[sps[0].key()].cluster_id].push_back(hit);
        hitmap[geo::PlaneID(hit->WireID())].push_back(std::make_pair(hit, fDBScan.points[sps[0].key()].cluster_id));
      }
    }
  }

  //Save hits not associated with any spacepoints
  for (auto &hit : hits){
    bool found = false;
    for (size_t i = 0; i<slcHits.size(); ++i){
      if (std::find(slcHits[i].begin(), slcHits[i].end(), hit) != slcHits[i].end()){
        found = true;
        break;
      }
    }
    if (!found){
      double wirePitch = fGeom->WirePitch(hit->WireID());
      double UnitsPerTick = tickToDist / wirePitch;
      double x0 = hit->WireID().Wire;
      double y0 = hit->PeakTime() * UnitsPerTick;
      double mindis = DBL_MAX;
      unsigned slcIndex = UINT_MAX;
      for (auto &hit2 : hitmap[geo::PlaneID(hit->WireID())]) {
        double dx = hit2.first->WireID().Wire - x0;
        double dy = hit2.first->PeakTime() * UnitsPerTick - y0;
        double dis = dx*dx + dy*dy;
        if (dis<mindis){
          mindis = dis;
          slcIndex = hit2.second;
        }
      }
      if (slcIndex != UINT_MAX && mindis < fMinHitDis){
        slcHits[slcIndex].push_back(hit);
      }
    }
  }

  //Save spacepoints for each slice
  std::vector<std::vector<art::Ptr<recob::SpacePoint>>> sps_in_slc(nslc);
  for (size_t i = 0; i<fDBScan.points.size(); ++i){
    if (fDBScan.points[i].cluster_id>=0){
      sps_in_slc[fDBScan.points[i].cluster_id].push_back(sps[i]);
    }
  } // i

  // calculate the slice aspect ratios
  std::vector<float> aspectRatio(nslc, 0);
  for(size_t isl = 0; isl < nslc; ++isl) {
    double sumx = 0, sumy = 0., sumz = 0., sumx2 = 0, sumy2 = 0, sumz2 = 0;
    double sumxy = 0, sumxz = 0;
    for(auto& spt : sps_in_slc[isl]) {
      double xx = spt->XYZ()[0];
      double yy = spt->XYZ()[1];
      double zz = spt->XYZ()[2];
      sumx += xx;
      sumy += yy;
      sumz += zz;
      sumx2 += xx * xx;
      sumy2 += yy * yy;
      sumz2 += zz * zz;
      sumxy += xx * yy;
      sumxz += xx * zz;
    } // spt
    double sum = sps_in_slc[isl].size();
    double delta = sum * sumx2 - sumx * sumx;
    if(delta <= 0) continue;
    float cnt = 0.;
    double arg = sum * sumy2 - sumy * sumy;
    if(arg > 0) {
      aspectRatio[isl] += std::abs(sum * sumxy - sumx * sumy) / sqrt(delta * arg);
      ++cnt;
    }
    arg = sum * sumz2 - sumz * sumz;
    if(arg > 0) {
      aspectRatio[isl] += std::abs(sum * sumxz - sumx * sumz) / sqrt(delta * arg);
      ++cnt;
    }
    if(cnt > 1) aspectRatio[isl] /= cnt;
  } // isl
  
  for(size_t isl = 0; isl < nslc; ++isl) {
    int id = isl + 1;
    slcCol->push_back(recob::Slice(id, aspectRatio[isl]));
    util::CreateAssn(*this, evt, *slcCol, slcHits[isl], *slc_hit_assn);
  } // i
  evt.put(std::move(slcCol));
  evt.put(std::move(slc_hit_assn));
}

DEFINE_ART_MODULE(cluster::DBCluster3D)
