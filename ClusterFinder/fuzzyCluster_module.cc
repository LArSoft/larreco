////////////////////////////////////////////////////////////////////////
// $Id: DBSCANfinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// \file fuzzyCluster_module.cc
//
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TGeoManager.h"
#include "TH1.h"

//Framework includes:
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/JamesRandom.h"

// LArSoft includes
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "RecoAlg/fuzzyClusterAlg.h"
#include "Filters/ChannelFilter.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/SeedCreator.h"


class TH1F;

namespace cluster{
   
  //--------------------------------------------------------------- 
  class fuzzyCluster : public art::EDProducer
  {
  public:
    explicit fuzzyCluster(fhicl::ParameterSet const& pset); 
    ~fuzzyCluster();
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
  
    std::string fhitsModuleLabel;
    long int fHoughSeed;
   
    fuzzyClusterAlg ffuzzyCluster; ///< object that implements the fuzzy cluster algorithm
  };

}

namespace cluster{


  //-------------------------------------------------
  fuzzyCluster::fuzzyCluster(fhicl::ParameterSet const& pset) :
    ffuzzyCluster(pset.get< fhicl::ParameterSet >("fuzzyClusterAlg"))
  {  
    this->reconfigure(pset);
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  
    // Create random number engine needed for PPHT
    createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");
  }
  
  //-------------------------------------------------
  fuzzyCluster::~fuzzyCluster()
  {
  }
  
  //-------------------------------------------------
  void fuzzyCluster::reconfigure(fhicl::ParameterSet const& p)
  {
    fhitsModuleLabel  = p.get< std::string >("HitsModuleLabel");
    fHoughSeed = p.get< long int >("HoughSeed");
    ffuzzyCluster.reconfigure(p.get< fhicl::ParameterSet >("fuzzyClusterAlg"));
  }
  
  //-------------------------------------------------
  void fuzzyCluster::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
  
    fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
    fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
    fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
      
  }
  
  //-----------------------------------------------------------------
  void fuzzyCluster::produce(art::Event& evt)
  {
     
    //get a collection of clusters   
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
  
    art::ServiceHandle<geo::Geometry> geom;
  
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fhitsModuleLabel,hitcol);
   
    // loop over all hits in the event and look for clusters (for each plane)
    std::vector<art::Ptr<recob::Hit> > allhits;

    // If a nonzero random number seed has been provided, 
    // overwrite the seed already initialized
    if(fHoughSeed != 0){
      art::ServiceHandle<art::RandomNumberGenerator> rng;
      CLHEP::HepRandomEngine &engine = rng->getEngine();
      engine.setSeed(fHoughSeed,0);
    } 

    // get the ChannelFilter
    filter::ChannelFilter chanFilt;
      
    // make a map of the geo::PlaneID to vectors of art::Ptr<recob::Hit>
    std::map<geo::PlaneID, std::vector< art::Ptr<recob::Hit> > > planeIDToHits;
    for(size_t i = 0; i < hitcol->size(); ++i)
      planeIDToHits[hitcol->at(i).WireID().planeID()].push_back(art::Ptr<recob::Hit>(hitcol, i));


    for(auto & itr : planeIDToHits){
      
      geo::SigType_t sigType = geom->SignalType(itr.first);
      allhits.resize(itr.second.size());
      allhits.swap(itr.second);
      
      //Begin clustering with fuzzy
      
      ffuzzyCluster.InitFuzzy(allhits, chanFilt.SetOfBadChannels());
      
      //----------------------------------------------------------------
      for(unsigned int j = 0; j < ffuzzyCluster.fps.size(); ++j){
  	
	if(allhits.size() != ffuzzyCluster.fps.size()) break;
  	
	fhitwidth->Fill(ffuzzyCluster.fps[j][2]);
  	
	if(sigType == geo::kInduction)  fhitwidth_ind_test->Fill(ffuzzyCluster.fps[j][2]);
	if(sigType == geo::kCollection) fhitwidth_coll_test->Fill(ffuzzyCluster.fps[j][2]);
      }
      
      //*******************************************************************
      ffuzzyCluster.run_fuzzy_cluster(allhits);
      
      //End clustering with fuzzy
      
      
      for(size_t i = 0; i < ffuzzyCluster.fclusters.size(); ++i){
	std::vector<art::Ptr<recob::Hit> > clusterHits;
	double totalQ = 0.;
  	
	for(size_t j = 0; j < ffuzzyCluster.fpointId_to_clusterId.size(); ++j){
	  if(ffuzzyCluster.fpointId_to_clusterId[j]==i){ 
	    clusterHits.push_back(allhits[j]);
	    totalQ += clusterHits.back()->Charge();
	  }
	} 
	
	
	////////
	if (clusterHits.size()>0){
	  /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	  const geo::WireID& wireID = clusterHits.front()->WireID();
	  unsigned int sw = wireID.Wire;
	  unsigned int ew = clusterHits.back()->WireID().Wire;
	  
	  // create the recob::Cluster directly in the vector
	  ccol->emplace_back(sw*1., 0.,
				 clusterHits.front()->PeakTime(), clusterHits.front()->SigmaPeakTime(),
				 ew*1., 0.,
				 clusterHits.back()->PeakTime(), clusterHits.back()->SigmaPeakTime(),
				 -999., 0., 
				 -999., 0.,
				 totalQ,
				 clusterHits.front()->View(),
				 ccol->size(),
				 wireID.planeID());
	  
	  // associate the hits to this cluster
	  util::CreateAssn(*this, evt, *ccol, clusterHits, *assn);
  	  
	  clusterHits.clear();
  	  
	}//end if clusterHits has at least one hit
     
      }//end loop over fclusters
  	
      allhits.clear();
    } // end loop over map
  
    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "fuzzyCluster Summary:";
    for(size_t i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;
  
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
  
    return;
  } // end produce
  
} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(fuzzyCluster)
  
} 







