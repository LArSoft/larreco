#ifndef SEEDFINDER_H
#define SEEDFINDER_H

//
// Name: SeedFinderModule.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SeedFinderAlgorithm.h"
#include "Utilities/AssociationUtil.h"


namespace recob
{
  class SpacePoint;
  class Cluster;
  class Seed;
  class Hit;
}

namespace trkf {

  class SeedFinderModule : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SeedFinderModule(fhicl::ParameterSet const& pset);
    virtual ~SeedFinderModule();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    
    void endJob();
  
    
    art::PtrVector<recob::Hit>       GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt);
    void                             GetSortedHitsFromClusters(std::string ClusterModuleLabel, art::Event& evt,std::vector< std::vector< art::PtrVector<recob::Hit> > > & SortedHits );



  private:

   
    // Fcl Attributes.

    SeedFinderAlgorithm     fSeedAlg;                  // Algorithm object
    std::string             fInputModuleLabel;         // Where to find hits, if we need them
    int                     fInputSource;              // 1: Use Clusters
                                                      // 2: Use Hits

  };
  
}

#endif // SEEDFINDER_H



#include "art/Framework/Core/ModuleMacros.h" 


namespace trkf {
  DEFINE_ART_MODULE(SeedFinderModule)
}

//
// Name: SeedFinderModule.cxx
//
// Purpose: Implementation file for module SeedFinderModule.
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
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace trkf {

  //----------------------------------------------------------------------------
  SeedFinderModule::SeedFinderModule(const fhicl::ParameterSet& pset) :
    fSeedAlg(pset.get<fhicl::ParameterSet>("SeedAlg"))
  {
    reconfigure(pset);
    produces<std::vector<recob::Seed> >();
  
  
  }

  //----------------------------------------------------------------------------
  SeedFinderModule::~SeedFinderModule()
  {
  }

  //----------------------------------------------------------------------------
  void SeedFinderModule::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSeedAlg.reconfigure ( pset.get<fhicl::ParameterSet>("SeedAlg") );
    fInputModuleLabel      = pset.get<std::string>("InputModuleLabel");
    fInputSource           = pset.get<int>("InputSource");
  
  }

  //----------------------------------------------------------------------------
  void SeedFinderModule::beginJob()
  {}




  //----------------------------------------------------------------------------
  void SeedFinderModule::produce(art::Event& evt)
  {
    
    std::unique_ptr<std::vector<recob::Seed> > seeds(new std::vector<recob::Seed>);

    std::vector<std::vector<recob::SpacePoint> > SpacePointsWithSeeds;
    std::vector<recob::Seed> SeedVector;

    // Make seeds from clusters    
    if(fInputSource==1)
      {
	std::vector<std::vector<art::PtrVector<recob::Hit> > > HitsPerSeed;
	
	std::vector< std::vector< art::PtrVector<recob::Hit> > > SortedHits;
        GetSortedHitsFromClusters(fInputModuleLabel, evt, SortedHits);

	std::vector<std::vector<recob::Seed> > Seeds = fSeedAlg.GetSeedsFromSortedHits(SortedHits, HitsPerSeed);

	for(size_t i=0; i!=Seeds.size(); ++i)
	  for(size_t j=0; j!=Seeds.at(i).size(); ++j)
	    SeedVector.push_back(Seeds.at(i).at(j));
      }
    
    // Make seeds from unsorted hits
    else if(fInputSource==0)
      {
	
	art::PtrVector<recob::Hit> Hits = GetHitsFromEvent(fInputModuleLabel, evt);
	std::vector<art::PtrVector<recob::Hit> > HitCatalogue;
  	SeedVector = fSeedAlg.GetSeedsFromUnSortedHits(Hits, HitCatalogue);
      }
    else	  
      {
	throw cet::exception("SeedFinderModule") << 
	  "Unkown source mode " << fInputSource<<"\n";
      }

    if(SeedVector.size()>0)
      {
	for(size_t i=0; i!=SeedVector.size(); ++i)
	  {
	    seeds->push_back(SeedVector.at(i));
	  }
	
      }
    else
      mf::LogInfo("SeedFinder")<<"Seed finder made no seeds : no space points in event"<<std::endl;
    
    evt.put(std::move(seeds));
    


  }




  //----------------------------------------------------------------------------
  // Get the hits associated with stored clusters
  //

  void SeedFinderModule::GetSortedHitsFromClusters(std::string ClusterModuleLabel, art::Event& evt,std::vector< std::vector< art::PtrVector<recob::Hit> > > & SortedHits )
  {
    
    SortedHits.clear();
    SortedHits.resize(3);
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

      art::PtrVector<recob::Hit>  HitsThisCluster;
      for(std::vector< art::Ptr<recob::Hit> >::const_iterator i = ihits.begin();
          i != ihits.end(); ++i)
        HitsThisCluster.push_back(*i);

      SortedHits[ThisCluster->View()].push_back(HitsThisCluster);
      
    }
  }






  
  //----------------------------------------------------------------------------
  // Extract vector of hits from event
  //

  art::PtrVector<recob::Hit> SeedFinderModule::GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt)
  {
    
    art::PtrVector<recob::Hit> TheHits;
    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(HitModuleLabel, hith);
    for(unsigned int i=0; i<hith->size(); ++i)
      {
	art::Ptr<recob::Hit> hit(hith, i);
	TheHits.push_back(hit);
      }
    
    return TheHits;
  }

  




  //----------------------------------------------------------------------------
  void SeedFinderModule::endJob()
  {

  }
}
