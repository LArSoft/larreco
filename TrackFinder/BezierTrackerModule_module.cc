#include "art/Persistency/Common/PtrVector.h"

#ifndef BEZIERTRACKERMOD_H
#define BEZIERTRACKERMOD_H

//
// Name: BezierTrackerModule.h
//
// Purpose: Header file for module BezierTrackerModule.  This modules makes
//          bezier tracks out of seed collections, or hits, or clusters
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SeedFinderAlgorithm.h"
#include "Geometry/Geometry.h"


namespace recob
{
  class Seed;
  class Track;
  class Hit;
}


namespace trkf {

  class BezierTrack;
  class BezierTrackerAlgorithm;
  class SpacePointAlg;

  class BezierTrackerModule : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackerModule(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackerModule();

    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    
    

  private:

    // Fcl Attributes.

    std::string fClusterModuleLabel;

    
    trkf::BezierTrackerAlgorithm * fBTrackAlg;
    

    void GetHitsFromClusters(std::string ClusterModuleLabel, art::Event& evt,     std::vector< std::vector<art::PtrVector<recob::Hit> > >& ReturnVec, std::vector<std::map<int, int> >& ClusterMap);
    
  };
}

#endif 

#include "art/Framework/Core/ModuleMacros.h" 


namespace trkf {
  DEFINE_ART_MODULE(BezierTrackerModule)
}

//
// Name: BezierTrackerModule.cxx
//
// Purpose: Implementation file for module BezierTrackerModule.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoObjects/BezierTrack.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/BezierTrackerAlgorithm.h"

namespace trkf {

  BezierTrackerModule::BezierTrackerModule(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Track> >();
    produces< art::Assns<recob::Track, recob::Vertex> >();
      
  }

  BezierTrackerModule::~BezierTrackerModule()
  {
  }

  void BezierTrackerModule::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel= pset.get<std::string>("ClusterModuleLabel");
    
    fBTrackAlg = new trkf::BezierTrackerAlgorithm(pset.get<fhicl::ParameterSet>("BezierTrackerAlgorithm"));
      
}

  void BezierTrackerModule::beginJob()
  {}


  void BezierTrackerModule::produce(art::Event& evt)
  {
 
    // Declare products to store

    std::unique_ptr< std::vector<recob::Track > > btracks ( new std::vector<recob::Track>);
    std::unique_ptr< std::vector<recob::Vertex > > vertices ( new std::vector<recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Track, recob::Hit > >  assnhit( new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Track, recob::Vertex > > assnvtx( new art::Assns<recob::Track, recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Track > > assnclus( new art::Assns<recob::Cluster, recob::Track>);
   

    std::vector<trkf::BezierTrack >           BTracks;
    
    std::vector<art::PtrVector<recob::Hit> >  HitsForAssns;
    
    std::vector<std::map<int, int> > ClusterMap;
    std::vector< std::vector<art::PtrVector<recob::Hit> > > SortedHits;
    // Produce appropriately organized hit object
    GetHitsFromClusters(fClusterModuleLabel, evt, SortedHits,ClusterMap);
 
    // Produce bezier tracks
    std::vector<std::vector<std::pair<int, int> > > ClustersUsed;
    BTracks = fBTrackAlg->MakeTracks(SortedHits, HitsForAssns, ClustersUsed);
    
    // Attempt to mitigate clustering imperfections
    fBTrackAlg->FilterOverlapTracks( BTracks, HitsForAssns, ClustersUsed );
    fBTrackAlg->SortTracksByLength(  BTracks, HitsForAssns, ClustersUsed );
    fBTrackAlg->MakeOverlapJoins(    BTracks, HitsForAssns, ClustersUsed ); 
    fBTrackAlg->SortTracksByLength(  BTracks, HitsForAssns, ClustersUsed );  
    fBTrackAlg->MakeDirectJoins(     BTracks, HitsForAssns, ClustersUsed );
    fBTrackAlg->FilterOverlapTracks( BTracks, HitsForAssns, ClustersUsed );    
    
    // Perform bezier vertexing
    std::vector<recob::Vertex> Vertices;
    std::vector<std::vector<int> > VertexMapping;
    fBTrackAlg->MakeVertexJoins(BTracks, Vertices, VertexMapping);
      
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	
	std::unique_ptr<recob::Track>  ToStore = BTracks.at(i).GetBaseTrack();
	btracks->push_back(*ToStore);
	util::CreateAssn(*this, evt, *(btracks.get()), HitsForAssns.at(i), *(assnhit.get()));
      }

    // Store vertex assns
    for(size_t v=0; v!=Vertices.size(); ++v)
      {
	vertices->push_back(Vertices.at(v));
	for(size_t t=0; t!=VertexMapping.at(v).size(); ++t)
	  {
	    util::CreateAssn(*this, evt, *(btracks), *(vertices), *(assnvtx.get()), v, v, VertexMapping[v][t]);
	  }
      }
    
    // Store cluster assns
    std::vector<art::Ptr<recob::Cluster> > Clusters;
    
    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);
    
    if(clusterh.isValid()) {
      art::fill_ptr_vector(Clusters, clusterh);
      

      for(size_t t=0; t!=btracks->size(); ++t)
	{
	  for(size_t i=0; i!=ClustersUsed[t].size(); ++i)
	    {
	      int View = ClustersUsed[t][i].first;
	      int Index = ClustersUsed[t][i].second;
	      int ClusterID = ClusterMap[View][Index];
	      util::CreateAssn(*this, evt, *(btracks), Clusters.at(ClusterID), *(assnclus.get()), t);
	    }
	}
    }

    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Storing in evt - check"<<std::endl;
    evt.put(std::move(btracks));
    evt.put(std::move(vertices));
    evt.put(std::move(assnhit));
    evt.put(std::move(assnvtx));
    evt.put(std::move(assnclus));

    // Now tidy up
    trkf::SpacePointAlg *Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();
    Sptalg->clearHitMap();
    for(size_t i=0; i!=HitsForAssns.size(); ++i)
      HitsForAssns.at(i).clear();
    HitsForAssns.clear();
  
  }
  //-----------------------------------------

  void BezierTrackerModule::GetHitsFromClusters(std::string ClusterModuleLabel, art::Event& evt,  std::vector<std::vector<art::PtrVector<recob::Hit> > > & ReturnVec, std::vector<std::map<int, int> >& ClusterMap )
  {
    ReturnVec.clear();
    ReturnVec.resize(3);
    ClusterMap.clear();
    ClusterMap.resize(3);

    std::vector<art::Ptr<recob::Cluster> > Clusters;
    
    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);
    
    if(clusterh.isValid()) {
      art::fill_ptr_vector(Clusters, clusterh);
    }
    
    art::FindManyP<recob::Hit> fm(clusterh, evt, ClusterModuleLabel);
      
    for(size_t iclus = 0; iclus < Clusters.size(); ++iclus) {
      art::Ptr<recob::Cluster> ThisCluster = Clusters.at(iclus);
      
      std::vector< art::Ptr<recob::Hit> > ihits = fm.at(iclus);
      
      art::PtrVector<recob::Hit> HitsThisCluster;    
      for(std::vector< art::Ptr<recob::Hit> >::const_iterator i = ihits.begin();
	  i != ihits.end(); ++i)
	HitsThisCluster.push_back(*i);

      int View = ThisCluster->View();
      ClusterMap[View][ReturnVec[View].size()] = iclus;      
      ReturnVec[View].push_back(HitsThisCluster);
    }
  }
  


  //----------------------------------------------------------------------
  void BezierTrackerModule::endJob()
  {

  }
}
