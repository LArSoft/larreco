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
  , fSPHitAssnLabel(p.get< art::InputTag >("fSPHitAssnLabel"))
  , fDBScan(p.get< fhicl::ParameterSet >("DBScan3DAlg"))
  , fMinHitDis(p.get< double >("MinHitDis"))
{
  produces< std::vector<recob::Cluster> >();
  produces< std::vector<recob::PFParticle> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
  produces< art::Assns<recob::PFParticle, recob::Cluster> >();
  produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();

  fGeom = &*(art::ServiceHandle<geo::Geometry>());
  fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  tickToDist = fDetProp->DriftVelocity(fDetProp->Efield(),fDetProp->Temperature());
  tickToDist *= 1.e-3 * fDetProp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns     
}

void cluster::DBCluster3D::produce(art::Event & evt)
{
  std::unique_ptr<std::vector<recob::Cluster>> 
    clucol(new std::vector<recob::Cluster>);
  std::unique_ptr<std::vector<recob::PFParticle>> 
    pfpcol(new std::vector<recob::PFParticle>);

  std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>> 
    cls_hit_assn(new art::Assns<recob::Cluster, recob::Hit>);

  std::unique_ptr<art::Assns<recob::PFParticle, recob::Cluster>>
    pfp_cls_assn(new art::Assns<recob::PFParticle, recob::Cluster>);

  std::unique_ptr<art::Assns<recob::PFParticle, recob::SpacePoint>>
    pfp_sps_assn(new art::Assns<recob::PFParticle, recob::SpacePoint>);

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

  //Find number of pfparticles
  int maxid = INT_MIN;
  for (size_t i = 0; i<fDBScan.points.size(); ++i){
//    std::cout<<"Space point index "<<i<<" "<<fDBScan.points[i].sp.key()<<" "<<fDBScan.points[i].cluster_id<<std::endl;
    if (fDBScan.points[i].cluster_id > maxid) maxid = fDBScan.points[i].cluster_id;
  }
  size_t npfp = 0;
  if (maxid>=0) npfp = maxid + 1;

  //Save hits associated with each pfparticle
  std::vector<std::vector<art::Ptr<recob::Hit>>> cluhits(npfp);
  //Save hits on each PlaneID with pfparticle index
  std::map<geo::PlaneID, std::vector<std::pair<art::Ptr<recob::Hit>, unsigned int>>> hitmap;
  for (auto &hit : hits){
    auto &sps = spFromHit.at(hit.key());
    if (sps.size()){//found associated space point
      if (fDBScan.points[sps[0].key()].cluster_id>=0){
        cluhits[fDBScan.points[sps[0].key()].cluster_id].push_back(hit);
        hitmap[geo::PlaneID(hit->WireID())].push_back(std::make_pair(hit, fDBScan.points[sps[0].key()].cluster_id));
      }
    }
  }

  //Save hits not associated with any spacepoints
  for (auto &hit : hits){
    bool found = false;
    for (size_t i = 0; i<cluhits.size(); ++i){
      if (std::find(cluhits[i].begin(), cluhits[i].end(), hit) != cluhits[i].end()){
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
      unsigned pfpindex = UINT_MAX;
      for (auto &hit2 : hitmap[geo::PlaneID(hit->WireID())]){
        double x1 = hit2.first->WireID().Wire;
        double y1 = hit2.first->PeakTime() * UnitsPerTick;
        double dis = sqrt(pow(x1-x0,2)+pow(y1-y0,2));
        if (dis<mindis){
          mindis = dis;
          pfpindex = hit2.second;
        }
      }
      if (pfpindex != UINT_MAX && mindis < fMinHitDis){
        cluhits[pfpindex].push_back(hit);
      }
    }
  }

  //Save spacepoints for each pfparticle
  std::vector<std::vector<art::Ptr<recob::SpacePoint>>> sps_in_pfp(npfp);
  for (size_t i = 0; i<fDBScan.points.size(); ++i){
    if (fDBScan.points[i].cluster_id>=0){
      sps_in_pfp[fDBScan.points[i].cluster_id].push_back(sps[i]);
    }
  }

  for (size_t i = 0; i < npfp; ++i){
    //Save pfparticle
    int pdgCode = 0;
    size_t self = pfpcol->size();
    size_t parent = 0;
    std::vector<size_t> daughters;
    pfpcol->push_back(recob::PFParticle(pdgCode, self, parent, daughters));
    size_t cluStart = clucol->size();
    //Find hits in each planeid
    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit>>> hitsview;
    for (auto &hit: cluhits[i]){
      hitsview[geo::PlaneID(hit->WireID())].push_back(hit);
    }
    for (auto & vhits : hitsview){
      if (vhits.second.size()){
        //Save cluster
        float start_wire = 0; float sigma_start_wire = 0; float start_tick = 0; float sigma_start_tick = 0; float start_charge = 0; float start_angle = 0; float start_opening = 0; float end_wire = 0; float sigma_end_wire = 0; float end_tick = 0; float sigma_end_tick = 0; float end_charge = 0; float end_angle = 0; float end_opening = 0; float integral = 0; float integral_stddev = 0; float summedADC = 0; float summedADC_stddev = 0; unsigned int n_hits = vhits.second.size(); float multiple_hit_density = 0; float width = 0; size_t ID = clucol->size(); geo::View_t view = vhits.second[0]->View(); geo::PlaneID plane = vhits.first;
        clucol->push_back(recob::Cluster(start_wire, sigma_start_wire, start_tick, sigma_start_tick, start_charge, start_angle, start_opening, end_wire, sigma_end_wire, end_tick, sigma_end_tick, end_charge, end_angle, end_opening, integral, integral_stddev, summedADC, summedADC_stddev, n_hits, multiple_hit_density, width, ID, view, plane));
        //Save hit-cluster association
        util::CreateAssn(*this, evt, *clucol, vhits.second, *cls_hit_assn);
      }
    }
    size_t cluEnd = clucol->size();
    //Save cluster-pfparticle association
    util::CreateAssn(*this, evt, *pfpcol, *clucol, *pfp_cls_assn, cluStart, cluEnd);
    //Save sps-pfparticle association
    util::CreateAssn(*this, evt, *pfpcol, sps_in_pfp[i], *pfp_sps_assn);
  }
    
  evt.put(std::move(pfpcol));
  evt.put(std::move(clucol));
  evt.put(std::move(cls_hit_assn));
  evt.put(std::move(pfp_cls_assn));
  evt.put(std::move(pfp_sps_assn));
}

DEFINE_ART_MODULE(cluster::DBCluster3D)
