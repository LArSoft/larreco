////////////////////////////////////////////////////////////////////////
//
// \file SmallClusterFinder.cxx
//
// \author corey.adams@yale.edu
//
// This algorithm is designed to find small clusters that could correspond to gammas
// or low energy electrons.
//
/*	There are two parameters that matter from the fcl file:
		fNHitsInClust is the number of hits that should be in these small clusters
			^-- Gamma events seem to rarely have more than 4 hits in the cluster
			^-- SN events are unclear.  Should this even be used for SN?
		fRadiusSizePar is the distance (in cm) between the small clusters and any other hits.
	
	This algorithm sorts the hits by plane, and then looks at each hit individually.  If
	there is a hit within RadiusSizePar, it gets added to a local list.  All other hits
	are ignored.  Then, if the number of hits that got added to the local list is greater 
	then NHitsInClust, the original hit is ignored.  If it's less, the original hit is 
	presumed to be part of a very small (or single hit) cluster.  So its added to the list
	of hits in the small cluster.
	
	All of the small clusters are then split apart into groups in the way you would expect.
	Each cluster is assigned an ID number to distinguish it, and the hits that aren't 
	identified as small clusters all end up in the "leftover" cluster.  The numbering scheme
	is ID = 100*iPlane + Cluster on that plane, and the leftover hits are the first (0th)
	cluster written out.
	
	All of the heavy lifting is done is the SmallClusterFinderAlg.	
	This module remains responsible for interacting with the framework, getting hits,
	writing clusters, etc.
	
	Thanks to Andrzej for the basic alg.

	-Corey
*/
// 
///////////////////////////////////////////////////////////////////////

#include <iostream>

// include the proper bit of the framework
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/ModuleMacros.h"

//All the larsoft goodies:
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/SmallClusterFinderAlg.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/AssociationUtil.h"


namespace cluster {

  class SmallClusterFinder : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit SmallClusterFinder(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~SmallClusterFinder();                               /**Destructor*/
    void beginJob();                                     
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Routine that finds the cluster and sets the dTdW of the 2D shower*/
   
  private:
    
    art::ServiceHandle<geo::Geometry> geom;
    
    //input parameters
    std::string fHitFinderModuleLabel; 
    bool verbose;
  //  double fRadiusSizePar;		//Determines the max radius of the cluster, must be separated
  //  double fNHitsInClust;		//Forces cluster to have a max number of hits
    							// Remember, this is the *small* cluster finder
    
  	SmallClusterFinderAlg fSmallClusterFinderAlg; //object that actually finds and sorts gammas

    unsigned int fNPlanes; // number of planes  
    
    int GetPlaneAndTPC(	art::Ptr<recob::Hit> a,
    					unsigned int &p,
    					unsigned int &cs,
    					unsigned int &t,
    					unsigned int &w);
    
  }; // class SmallAngleFinder

}
namespace cluster{
	SmallClusterFinder::SmallClusterFinder(fhicl::ParameterSet const& pset)
	{
		this->reconfigure(pset);
		produces< std::vector<recob::Cluster> >();				//This code makes clusters
		produces< art::Assns<recob::Cluster, recob::Hit>  >();  //Matches clusters with hits
	}

	void SmallClusterFinder::reconfigure(fhicl::ParameterSet const& pset) 
	{
		fHitFinderModuleLabel 	=pset.get< 	std::string > ("HitFinderModuleLabel"); 
		verbose					=pset.get<	bool		> ("Verbose");
	
		//Let the clusterAlg have access to the pset too:
		fSmallClusterFinderAlg.reconfigure(pset.get<fhicl::ParameterSet> ("smallClustAlg") );
	 }

	// ***************** //
	SmallClusterFinder::~SmallClusterFinder()
	{
		//Nothing to do in the destructor
	}

	//____________________________________________________________________________
	void SmallClusterFinder::beginRun(art::Run& run)
	  {
		//nothing to do at beginRun()
		return;
	  }

	//-----------------------------------------------

	// ***************** //
	void SmallClusterFinder::beginJob()
	{
		// this will not change on a run per run basis.
		fNPlanes = geom->Nplanes(); 				//get the number of planes in the TPC
		/**Get TFileService and define output Histograms*/
		art::ServiceHandle<art::TFileService> tfs;
		return;
	}
  
  
  
	// ***************** //
	// This method actually makes the clusters.
	void SmallClusterFinder::produce(art::Event& evt)
	{ 
	  /**Get Clusters*/
	
		//Get the hits for this event:
		art::Handle< std::vector<recob::Hit> > HitListHandle;
		evt.getByLabel(fHitFinderModuleLabel,HitListHandle);  

		//A vector to hold hits, not yet filled:
		std::vector< art::Ptr<recob::Hit> > hitlist;
	
		//How many hits in this event?  Tell user:
		if (verbose) std::cout << " ++++ Hitsreceived received " << HitListHandle->size() << " +++++ " << std::endl;
		//Catch the case were there are no hits in the event:
		if(HitListHandle->size() ==0 )
		{
			if (verbose) std::cout << "No hits received! Exiting." << std::endl;
			return;
		}
		hitlist.resize(HitListHandle->size());

		//wrap the hits in art::Ptrs to pass to the Alg
		for (unsigned int iHit = 0; iHit < hitlist.size(); iHit++){
			hitlist[iHit] = art::Ptr< recob::Hit>(HitListHandle, iHit);
		}
	
		//std::cout << "Passing " << hitlist.size() << " hits to the alg." << std::endl;

  
		//Now run the alg to find the gammas:
		fSmallClusterFinderAlg.FindSmallClusters(hitlist);
	
  
 
		// make an art::PtrVector of the clusters
		std::unique_ptr< std::vector<recob::Cluster> > SmallClusterFinder(new std::vector<recob::Cluster>);
		std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
	
  
		for(unsigned int iplane=0;iplane<fNPlanes;iplane++){
		
			//Get the leftover hits for this plane:
			std::vector<art::Ptr<recob::Hit> > leftoverHits = fSmallClusterFinderAlg.GetLeftoversByPlane(iplane);    
	
			//write the leftover hits as a cluster:
			if (leftoverHits.size() != 0){
				if (verbose) std::cout << "Writing leftover hits to cluster ID: " << iplane*100 << std::endl;
				recob::Cluster leftover( 0,  0, 	//start wire, error in start wire
								 0,  0,		//start time, error in start time
								 0., 0.,	//end   wire, error in end   wire
								 0., 0.,  	//end   time, error in end   wire
								 0., 0.,	//dTdW, error in dTdW ??
								 0., 0.,	//dQdW, error in dQdW ??
								 5.,		//Total Q
								 geom->Plane(iplane,0,0).View(), 	//View (geo::View_t)
								 iplane*100);	//id for cluster, given as plane here
	
				SmallClusterFinder->push_back(leftover);
				util::CreateAssn(*this, evt, *SmallClusterFinder, leftoverHits, *assn);  
			} //leftovers are written for this plane, if they exist.
		
			//Now get the small clusters;
			std::vector< std::vector<art::Ptr<recob::Hit> > > smallClusters;
			smallClusters = fSmallClusterFinderAlg.GetSmallClustersByPlane(iplane);
		
			for (unsigned int i = 0; i < smallClusters.size(); i++){
				recob::Cluster clust( 0,  0, 	//start wire, error in start wire
								 0,  0,		//start time, error in start time
								 0., 0.,	//end   wire, error in end   wire
								 0., 0.,  	//end   time, error in end   wire
								 0., 0.,	//dTdW, error in dTdW ??
								 0., 0.,	//dQdW, error in dQdW ??
								 5.,		//Total Q
								 geom->Plane(iplane,0,0).View(), 	//View (geo::View_t)
								 iplane*100 + i + 1);	//id for cluster, given as plane here

	
				SmallClusterFinder->push_back(clust);
				// associate the hits to this cluster
				util::CreateAssn(*this, evt, *SmallClusterFinder, smallClusters[i], *assn);
			}
		
			//Just in case there is nothing for this event, make an empty cluster
			//so that the clusters are defined.  Tag with ID -1.
		/*	if (smallClusters.size() == 0 && leftoverHits.size() == 0){
				std::vector< art::Ptr < recob::Hit> > emptyVector;
				recob::Cluster clust( 0,  0, 	//start wire, error in start wire
								 0,  0,		//start time, error in start time
								 0., 0.,	//end   wire, error in end   wire
								 0., 0.,  	//end   time, error in end   wire
								 0., 0.,	//dTdW, error in dTdW ??
								 0., 0.,	//dQdW, error in dQdW ??
								 0.,		//Total Q
								 geom->Plane(iplane,0,0).View(), 	//View (geo::View_t)
								 -iplane);	//id for cluster, given as plane here

	
				SmallClusterFinder->push_back(clust);
				// associate the hits to this cluster
				util::CreateAssn(*this, evt, *SmallClusterFinder, emptyVector, *assn);
			}
		*/
			
		}
	
		//Finish up:
		evt.put(std::move(SmallClusterFinder));
		evt.put(std::move(assn));
 
		return;
	} //end produce

} // end namespace cluster

namespace cluster {

  DEFINE_ART_MODULE(SmallClusterFinder)

}


