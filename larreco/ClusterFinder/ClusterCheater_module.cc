////////////////////////////////////////////////////////////////////////
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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "lardata/Utilities/GeometryUtilities.h"


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
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"

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
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    
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
    pi_serv->SetEveIdCalculator(new sim::EmEveIdCalculator);

    LOG_DEBUG("ClusterCheater") << pi_serv->ParticleList();

    // make a map of vectors of art::Ptrs keyed by eveID values and 
    // location in cryostat, TPC, plane coordinates of where the hit originated
    std::map< eveLoc, std::vector< art::Ptr<recob::Hit> > > eveHitMap;

    // loop over all hits and fill in the map
    for( auto const& itr : hits ){ 

      std::vector<sim::TrackIDE> eveides = bt_serv->HitToEveTrackIDEs(itr);

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

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;
    
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
      const simb::MCParticle *part = pi_serv->TrackIdToParticle_P(hitMapItr.first.eveID);
      if(!part) continue;

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
      double startTime = hitMapItr.second.front()->PeakTimeMinusRMS();
      double endWire   = hitMapItr.second.back()->WireID().Wire;
      double endTime   = hitMapItr.second.back()->PeakTimePlusRMS();

      // add a cluster to the collection.  Make the ID be the eve particle
      // trackID*1000 + plane number*100 + tpc*10 + cryostat that the current hits are from
      ///\todo: The above encoding of the ID probably won't work for DUNE and should be revisited
      const geo::PlaneID& planeID = hitMapItr.first.planeID;
      recob::Cluster::ID_t clusterID = (((
               hitMapItr.first.eveID
        )*10 + planeID.Plane
        )*10 + planeID.TPC // 10 is weird choice for DUNE FD... should be 1000! FIXME
        )*10 + planeID.Cryostat
        ;
      
      // feed the algorithm with all the cluster hits
      ClusterParamAlgo.ImportHits(hitMapItr.second);
      
      // create the recob::Cluster directly in the vector
      ClusterCreator cluster(
        ClusterParamAlgo,               // algo
        startWire,                      // start_wire
        0.,                             // sigma_start_wire
        startTime,                      // start_tick
        0.,                             // sigma_start_tick
        endWire,                        // end_wire 
        0.,                             // sigma_end_wire
        endTime,                        // end_tick
        0.,                             // sigma_end_tick
        clusterID,                      // ID
        hitMapItr.second.at(0)->View(), // view
        planeID,                        // plane
        recob::Cluster::Sentry          // sentry
        );
      
      clustercol->emplace_back(cluster.move());
      
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


