////////////////////////////////////////////////////////////////////////
// Class:       TCShower
// Plugin Type: producer (art v2_11_02)
// File:        TCShower_module.cc
//
// Generated at Fri Jun  8 14:55:04 2018 by Rory Fitzpatrick using cetskelgen
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
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

#include <memory>

namespace shower {
  class TCShower;
}


class shower::TCShower : public art::EDProducer {
public:
  explicit TCShower(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TCShower(TCShower const &) = delete;
  TCShower(TCShower &&) = delete;
  TCShower & operator = (TCShower const &) = delete;
  TCShower & operator = (TCShower &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  std::string fClusterModuleLabel;
  std::string fTrackModuleLabel;

  // Declare member data here.

};


shower::TCShower::TCShower(fhicl::ParameterSet const & pset) : 
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel", "trajcluster2" ) ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel", "pmtrack" ) )
{

  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  //  produces<art::Assns<recob::Shower, recob::Cluster> >();
  //  produces<art::Assns<recob::Shower, recob::Track> >();
}

void shower::TCShower::produce(art::Event & evt)
{
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  //  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitShowerAssociations(new art::Assns<recob::Shower, recob::Hit>);
  //  std::unique_ptr<art::Assns<recob::Shower, recob::Track> > trackAssociations(new art::Assns<recob::Shower, recob::Track>);

  // clusters
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  // tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

  // detector properties
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  art::FindManyP<recob::Hit> cls_fm(clusterListHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Hit> trk_fm(trackListHandle, evt, fTrackModuleLabel);

  int tolerance = 20;

  double maxDist = 10;

  for (size_t i = 0; i < tracklist.size(); ++i) {

    std::vector< art::Ptr<recob::Hit> > trk_hitlist = trk_fm.at(i);
    std::vector< art::Ptr<recob::Hit> > showerHits;
    
    for (size_t ii = 0; ii < trk_hitlist.size(); ++ii) {
      showerHits.push_back(trk_hitlist[ii]);
      std::cout<<trk_hitlist[ii]->WireID().Plane<<" "<<trk_hitlist[ii]->WireID().Wire<<" "<<trk_hitlist[ii]->PeakTime()<<std::endl;
    } // loop over track hits


    //TrajectoryPoint(size_t i)

    int nShowerHits = 0;
    bool showerCandidate = false;

    TVector3 trkStart = tracklist[i]->Vertex();
    TVector3 trkDir = tracklist[i]->VertexDirection();

    TVector3 trkPt2 = trkStart + trkDir;

    std::vector<double> trk_tick1(2); 
    std::vector<double> trk_wire1(2);

    std::vector<double> trk_tick2(2); 
    std::vector<double> trk_wire2(2);

    trk_tick1[0] = detprop->ConvertXToTicks(trkStart[0], 0, 0, 0);
    trk_tick1[1] = detprop->ConvertXToTicks(trkStart[0], 1, 0, 0); // second argument is plane

    trk_wire1[0] = geom->WireCoordinate(trkStart[1], trkStart[2], geo::PlaneID(0, 0, 0));
    trk_wire1[1] = geom->WireCoordinate(trkStart[1], trkStart[2], geo::PlaneID(0, 0, 1)); // last argument is plane

    trk_tick2[0] = detprop->ConvertXToTicks(trkPt2[0], 0, 0, 0);
    trk_tick2[1] = detprop->ConvertXToTicks(trkPt2[0], 1, 0, 0); // second argument is plane

    trk_wire2[0] = geom->WireCoordinate(trkPt2[1], trkPt2[2], geo::PlaneID(0, 0, 0));
    trk_wire2[1] = geom->WireCoordinate(trkPt2[1], trkPt2[2], geo::PlaneID(0, 0, 1)); // last argument is plane
    
    for (size_t j = 0; j < clusterlist.size(); ++j) {
      if (clusterlist[j]->ID() > 0) continue;

      std::vector< art::Ptr<recob::Hit> > cls_hitlist = cls_fm.at(j); 
      
      for (size_t k = 0; k < cls_hitlist.size(); ++k) {
	int planeNum = cls_hitlist[k]->WireID().Plane;

	double wirePitch = geom->WirePitch(planeNum);         
	double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());                    
	tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns 
	double UnitsPerTick = tickToDist / wirePitch; 

	double x0 = cls_hitlist[k]->WireID().Wire; 
	double y0 = cls_hitlist[k]->PeakTime() * UnitsPerTick;

	double x1 = trk_wire1[planeNum];	
	double y1 = trk_tick1[planeNum] * UnitsPerTick;

	double x2 = trk_wire2[planeNum];	
	double y2 = trk_tick2[planeNum] * UnitsPerTick;

	double dist = std::abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/std::sqrt( pow((y2-y1), 2) + pow((x2-x1), 2) );
	
	if (dist<maxDist) {
	  nShowerHits++;
	  showerHits.push_back(cls_hitlist[k]);
	}


      } // loop over hits in cluster

    } // loop over clusters

    if (nShowerHits > tolerance) showerCandidate = true;

    TVector3 dcosVtxErr;
    TVector3 xyzErr;

    std::vector<double> totalEnergy(2);
    std::vector<double> totalEnergyErr(2);

    std::vector<double> dEdx(2);
    std::vector<double> dEdxErr(2);

    if (showerCandidate) {

      std::cout << "track ID " << tracklist[i]->ID() << std::endl;
      showers->push_back(recob::Shower(trkDir, dcosVtxErr, trkStart, xyzErr, totalEnergy, totalEnergyErr, dEdx, dEdx, 0, 0));

      showers->back().set_id(showers->size()-1);
      
      util::CreateAssn(*this, evt, *(showers.get()), showerHits, *(hitShowerAssociations.get()) );

      break;

    }


  } // loop over tracks

  evt.put(std::move(showers));
  evt.put(std::move(hitShowerAssociations));

} // produce

void shower::TCShower::beginJob()
{

  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(shower::TCShower)
