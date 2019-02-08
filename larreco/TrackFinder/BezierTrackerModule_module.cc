#include "canvas/Persistency/Common/PtrVector.h"

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
#include "canvas/Persistency/Common/FindManyP.h"
#include "larreco/RecoAlg/SeedFinderAlgorithm.h"
#include "larcore/Geometry/Geometry.h"


namespace recob
{
  class Seed;
  class Trajectory;
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
    double      fTrajPtPitch;
    
    std::unique_ptr<trkf::BezierTrackerAlgorithm> fBTrackAlg;
    

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
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/BezierTrackerAlgorithm.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace

namespace trkf {

  BezierTrackerModule::BezierTrackerModule(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
  {
    reconfigure(pset);
    produces< std::vector<recob::Trajectory> >("bezierformat");
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Trajectory, recob::Hit> >("bezierformat");
    produces< art::Assns<recob::Cluster, recob::Trajectory> >("bezierformat");
    produces< art::Assns<recob::Trajectory, recob::Vertex> >("bezierformat");
      
    produces< std::vector<recob::Trajectory> >();
    produces< art::Assns<recob::Trajectory, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Trajectory> >();
    produces< art::Assns<recob::Trajectory, recob::Vertex> >();

  }

  BezierTrackerModule::~BezierTrackerModule()
  {
  }

  void BezierTrackerModule::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fTrajPtPitch        = pset.get<double>("TrajPtPitch",0.1);
    
    fBTrackAlg.reset( new trkf::BezierTrackerAlgorithm(pset.get<fhicl::ParameterSet>("BezierTrackerAlgorithm")) );
      
}

  void BezierTrackerModule::beginJob()
  {}


  void BezierTrackerModule::produce(art::Event& evt)
  {
 
    // Declare products to store

    std::unique_ptr< std::vector<recob::Trajectory > > btracks ( new std::vector<recob::Trajectory>);
    std::unique_ptr< std::vector<recob::Vertex > > vertices ( new std::vector<recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Trajectory, recob::Hit > >  assnhit( new art::Assns<recob::Trajectory, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Trajectory, recob::Vertex > > assnvtx( new art::Assns<recob::Trajectory, recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Trajectory> > assnclus( new art::Assns<recob::Cluster, recob::Trajectory>);
   
    std::unique_ptr< std::vector<recob::Trajectory > > ub_tracks ( new std::vector<recob::Trajectory>);
    std::unique_ptr< art::Assns<recob::Trajectory, recob::Hit > >  ub_assnhit( new art::Assns<recob::Trajectory, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Trajectory, recob::Vertex > > ub_assnvtx( new art::Assns<recob::Trajectory, recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Trajectory > > ub_assnclus( new art::Assns<recob::Cluster, recob::Trajectory>);

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
	
	btracks->push_back(BTracks.at(i).GetTrajectory());

	std::vector<TVector3> track_xyz,track_dir;
	BTracks.at(i).FillTrackVectors(track_xyz,track_dir,fTrajPtPitch);

	ub_tracks->emplace_back(geo::vect::convertCollToPoint(track_xyz), geo::vect::convertCollToVector(track_dir), false);

	util::CreateAssn(*this, evt, *(btracks.get()), HitsForAssns.at(i), *(assnhit.get()));
	util::CreateAssn(*this, evt, *(ub_tracks), HitsForAssns.at(i), *(assnhit.get()));
      }

    // Store vertex assns
    for(size_t v=0; v!=Vertices.size(); ++v)
      {
	vertices->push_back(Vertices.at(v));
	for(size_t t=0; t!=VertexMapping.at(v).size(); ++t)
	  {
	    util::CreateAssn(*this, evt, *(btracks), *(vertices), *(assnvtx.get()), v, v, VertexMapping[v][t]);
	    util::CreateAssn(*this, evt, *(ub_tracks), *(vertices), *(assnvtx.get()), v, v, VertexMapping[v][t]);
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
	      util::CreateAssn(*this, evt, *(ub_tracks), Clusters.at(ClusterID), *(assnclus.get()), t);
	    }
	}
    }

    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Storing in evt - check"<<std::endl;
    evt.put(std::move(btracks),"bezierformat");
    evt.put(std::move(vertices));
    evt.put(std::move(assnhit),"bezierformat");
    evt.put(std::move(assnvtx),"bezierformat");
    evt.put(std::move(assnclus),"bezierformat");

    evt.put(std::move(ub_tracks));
    evt.put(std::move(ub_assnhit));
    evt.put(std::move(ub_assnvtx));
    evt.put(std::move(ub_assnclus));

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
