////////////////////////////////////////////////////////////////////////
// $Id: ClusterCheater_module.cc Exp $
//
// ClusterCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTER_CLUSTERCHEATER_H
#define CLUSTER_CLUSTERCHEATER_H
#include <string>
#include <vector>
#include <algorithm>

// ROOT includes
#include "TStopwatch.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/AssociationUtil.h"
#include "Simulation/EmEveIdCalculator.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"


// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

namespace cluster {
  class ClusterCheater : public art::EDProducer {
  public:
    explicit ClusterCheater(fhicl::ParameterSet const& pset);
    virtual ~ClusterCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:
    
    std::string  fMCGeneratorLabel;  ///< label for module to get MC truth information
    std::string  fHitModuleLabel;    ///< label for module creating recob::Hit objects	   
    std::string  fG4ModuleLabel;     ///< label for module running G4 and making particles, etc
    unsigned int fMinHits;           ///< minimum number of hits to make a cluster

  };
}

namespace cluster{

  struct eveLoc{

    eveLoc(int id, geo::PlaneID plnID)
      : eveID(id)
      , planeID(plnID)
    {}

    friend bool operator < (eveLoc const& a, eveLoc const& b)
    { 
      if(a.eveID    != b.eveID)    
	return a.eveID    < b.eveID;
      
      return a.planeID < b.planeID;
    }
    
    int          eveID;
    geo::PlaneID planeID;
  };

  bool sortHitsByWire(art::Ptr<recob::Hit> a, art::Ptr<recob::Hit> b)
  {
    return a->WireID().Wire < b->WireID().Wire;
  }

  //--------------------------------------------------------------------
  ClusterCheater::ClusterCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  }

  //--------------------------------------------------------------------
  ClusterCheater::~ClusterCheater()
  {
  }

  //--------------------------------------------------------------------
  void ClusterCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMCGeneratorLabel  = pset.get< std::string  >("MCGeneratorLabel",  "generator");
    fHitModuleLabel    = pset.get< std::string  >("HitModuleLabel",    "hit"      );
    fG4ModuleLabel     = pset.get< std::string  >("G4ModuleLabel",     "largeant" );
    fMinHits           = pset.get< unsigned int >("MinHits",           1          );

    return;
  }

  //--------------------------------------------------------------------
  void ClusterCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<geo::Geometry>      geo;
    art::ServiceHandle<cheat::BackTracker> bt;
    
    // grab the hits that have been reconstructed
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fHitModuleLabel, hitcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitcol);
    
    // adopt an EmEveIdCalculator to find the eve ID.  
    // will return a primary particle if it doesn't find 
    // a responsible particle for an EM process
    bt->SetEveIdCalculator(new sim::EmEveIdCalculator);

    LOG_DEBUG("ClusterCheater") << bt->ParticleList();

    // make a map of vectors of art::Ptrs keyed by eveID values and 
    // location in cryostat, TPC, plane coordinates of where the hit originated
    std::map< eveLoc, std::vector< art::Ptr<recob::Hit> > > eveHitMap;

    // loop over all hits and fill in the map
    for( auto const& itr : hits ){

      std::vector<cheat::TrackIDE> eveides = bt->HitToEveID(itr);

      // loop over all eveides for this hit
      for(size_t e = 0; e < eveides.size(); ++e){

	// don't worry about eve particles that contribute less than 10% of the
	// energy in the current hit
	if( eveides[e].energyFrac < 0.1) continue;

	eveLoc el(eveides[e].trackID, 
		  itr->WireID().planeID());

	eveHitMap[el].push_back(itr);

      } // end loop over eve IDs for this hit

    }// end loop over hits

    // loop over the map and make clusters
    std::unique_ptr< std::vector<recob::Cluster> > clustercol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

    for(auto hitMapItr : eveHitMap){

      // ================================================================================
      // === Only keeping clusters with fMinHits 
      // ================================================================================
      if(hitMapItr.second.size() < fMinHits) continue;
	            
      // get the center of this plane in world coordinates
      double xyz[3]   = {0.};
      double xyz2[3]  = {0.};
      double local[3] = {0.};
      unsigned int cryostat = hitMapItr.first.planeID.Cryostat;
      unsigned int tpc      = hitMapItr.first.planeID.TPC;
      unsigned int plane    = hitMapItr.first.planeID.Plane;
      geo->Cryostat(cryostat).TPC(tpc).Plane(plane).LocalToWorld(local, xyz);

      LOG_DEBUG("ClusterCheater") << "make cluster for eveID: " << hitMapItr.first.eveID
				  << " in cryostat: "           << cryostat
				  << " tpc: "         	        << tpc     
				  << " plane: "       	        << plane
				  << " view: "                  << hitMapItr.second.at(0)->View();

      // get the direction of this particle in the current cryostat, tpc and plane
      const simb::MCParticle *part = bt->TrackIDToParticle(hitMapItr.first.eveID);

      // now set the y and z coordinates of xyz to be the first point on the particle
      // trajectory and use the initial directions to determine the dT/dW
      // multiply the direction cosine by 10 to give a decent lever arm for determining 
      // dW
      xyz[1]  = part->Vy();
      xyz[2]  = part->Vz();
      xyz2[0] = xyz[0];
      xyz2[1] = xyz[1] + 10.*part->Py()/part->P();
      xyz2[2] = xyz[2] + 10.*part->Pz()/part->P();

      // convert positions to wire and time
      unsigned int w1 = 0;
      unsigned int w2 = 0;
			
      try{
	w1 = geo->NearestWire(xyz, plane, tpc, cryostat); 
      }
      catch(cet::exception& e){
	w1 = atoi(e.explain_self().substr(e.explain_self().find("#")+1,5).c_str());
      }
      try{
	w2 = geo->NearestWire(xyz2, plane, tpc, cryostat); 
      }
      catch(cet::exception& e){
	w2 = atoi(e.explain_self().substr(e.explain_self().find("#")+1,5).c_str());
      }

      // sort the vector of hits with respect to the directionality of the wires determined by 
      if(w2 < w1) 
	std::sort(hitMapItr.second.rbegin(), hitMapItr.second.rend(), sortHitsByWire);
      else
	std::sort(hitMapItr.second.begin(),  hitMapItr.second.end(),  sortHitsByWire);

      // set the start and end wires and times
      double startWire = hitMapItr.second.front()->WireID().Wire;
      double startTime = hitMapItr.second.front()->StartTime();
      double endWire   = hitMapItr.second.back()->WireID().Wire;
      double endTime   = hitMapItr.second.back()->EndTime();
      double totalQ    =  0.;
      double dTdW      =  1.e6;
      double dQdW      =  1.e6;
      
      if(startWire != endWire){
	dTdW = (endTime - startTime)/(endWire - startWire);
	///\todo now figure out the dQdW      
      }

      for(auto const& h : hitMapItr.second) totalQ += h->Charge();

      // add a cluster to the collection.  Make the ID be the eve particle
      // trackID*1000 + plane number*100 + tpc*10 + cryostat that the current hits are from
      ///\todo: The above encoding of the ID probably won't work for LBNE and should be revisited
      
      clustercol->push_back(recob::Cluster(startWire, 0.,
					   startTime, 0.,
					   endWire,   0.,
					   endTime,   0.,
					   dTdW,      0.,
					   dQdW,      0.,
					   totalQ,
					   hitMapItr.second.at(0)->View(),
					   (hitMapItr.first.eveID*1000 + 
					    hitMapItr.first.planeID.Plane*100  + 
					    hitMapItr.first.planeID.TPC*10     + 
					    hitMapItr.first.planeID.Cryostat)
					   )
			    );
      
      // association the hits to this cluster
      util::CreateAssn(*this, evt, *clustercol, hitMapItr.second, *assn);
      
      mf::LogInfo("ClusterCheater") << "adding cluster: \n" 
				    << clustercol->back()
				    << "\nto collection.";
      
    } // end loop over the map

    evt.put(std::move(clustercol));
    evt.put(std::move(assn));

    return;

  } // end produce

} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(ClusterCheater)

}

#endif


