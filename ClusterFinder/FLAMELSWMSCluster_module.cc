////////////////////////////////////////////////////////////////////////
//
// \file FLAMELSWMSCluster_module.cc
//
// \author bcarls@fnal.gov
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
#include "RecoAlg/FLAMEClusterAlg.h"
#include "RecoAlg/CornerFinderAlg.h"
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
  class FLAMELSWMSCluster : public art::EDProducer
  {
  public:
    explicit FLAMELSWMSCluster(fhicl::ParameterSet const& pset); 
    ~FLAMELSWMSCluster();
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
  
    std::string fhitsModuleLabel;
    std::string fCalDataModuleLabel;
   
    FLAMEClusterAlg fFLAMECluster; ///< object that implements the FLAMELSWMS cluster algorithm
  };

}

namespace cluster{


  //-------------------------------------------------
  FLAMELSWMSCluster::FLAMELSWMSCluster(fhicl::ParameterSet const& pset) :
    fFLAMECluster(pset.get< fhicl::ParameterSet >("FLAMEClusterAlg"))
  {  
    this->reconfigure(pset);
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  
  }
  
  //-------------------------------------------------
  FLAMELSWMSCluster::~FLAMELSWMSCluster()
  {
  }
  
  //-------------------------------------------------
  void FLAMELSWMSCluster::reconfigure(fhicl::ParameterSet const& p)
  {
    fhitsModuleLabel  = p.get< std::string >("HitsModuleLabel");
    fCalDataModuleLabel  = p.get< std::string >("CalDataModuleLabel");
    fFLAMECluster.reconfigure(p.get< fhicl::ParameterSet >("FLAMEClusterAlg"));
  }
  
  //-------------------------------------------------
  void FLAMELSWMSCluster::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
  
    fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
    fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
    fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
      
  }
  
  //-----------------------------------------------------------------
  void FLAMELSWMSCluster::produce(art::Event& evt)
  {
     
    //get a collection of clusters   
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
  
    art::ServiceHandle<geo::Geometry> geom;
  
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fhitsModuleLabel,hitcol);

    art::Handle< std::vector<recob::Wire> > wirecol;
    evt.getByLabel(fCalDataModuleLabel,wirecol);
   
    // loop over all hits in the event and look for clusters (for each plane)
    std::vector<art::Ptr<recob::Hit> > allhits;

    // loop over all end points in the event to help look for clusters (for each plane)
    //std::vector<art::Ptr<recob::EndPoint2D> > allends;
    std::vector<recob::EndPoint2D> allends;

    std::vector<recob::EndPoint2D> endcol; 

    // Pass information into CornerFinder
    //fcornerfinder.TakeInRaw(evt);
    //fcornerfinder.get_feature_points_LineIntegralScore(endcol);

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

      //double maxStrength = 0;
      //int iMaxStrength = -1;
      //for(size_t i = 0; i< endcol.size(); ++i){
      //if(endcol[i].WireID().Plane    == plane && 
      //endcol[i].WireID().TPC      == tpc   && 
      //endcol[i].WireID().Cryostat == cstat){allends.push_back(endcol[i]);  
      //std::cout << "wire: " << endcol[i].WireID().Wire << " time: " << endcol[i].DriftTime() << " strength: " << endcol[i].Strength() << std::endl;
      //if(maxStrength < endcol[i].Strength()){
      //maxStrength = endcol[i].Strength();
      //iMaxStrength = i;
      //}
      //}
      //}  

      //std::cout << "Strongest vertex" << std::endl;
      //std::cout << "wire: " << endcol[iMaxStrength].WireID().Wire << " time: " << endcol[iMaxStrength].DriftTime() << " strength: " << endcol[iMaxStrength].Strength() << std::endl;
      











      
      //Begin clustering with FLAME
      
      
      fFLAMECluster.InitFLAME(allhits, chanFilt.SetOfBadChannels());

      //----------------------------------------------------------------
      for(unsigned int j = 0; j < fFLAMECluster.fps.size(); ++j){
  	
	if(allhits.size() != fFLAMECluster.fps.size()) break;
  	
	fhitwidth->Fill(fFLAMECluster.fps[j][2]);
  	
	if(sigType == geo::kInduction)  fhitwidth_ind_test->Fill(fFLAMECluster.fps[j][2]);
	if(sigType == geo::kCollection) fhitwidth_coll_test->Fill(fFLAMECluster.fps[j][2]);
      }



      fFLAMECluster.run_FLAME_cluster(allhits);
      
      //End clustering with FLAMELSWMS
      
     









      for(size_t i = 0; i < fFLAMECluster.fclusters.size(); ++i){
	std::vector<art::Ptr<recob::Hit> > clusterHits;
	double totalQ = 0.;
  	
	for(size_t j = 0; j < fFLAMECluster.fpointId_to_clusterId.size(); ++j){
	  if(fFLAMECluster.fpointId_to_clusterId[j]==i){ 
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
    mf::LogVerbatim("Summary") << "FLAMELSWMSCluster Summary:";
    for(size_t i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;
  
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
  
    return;
  } // end produce
  
} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(FLAMELSWMSCluster)
  
} 







