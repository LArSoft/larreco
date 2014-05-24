////////////////////////////////////////////////////////////////////////
//
//  
////////////////////////////////////////////////////////////////////////

#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

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

//LArSoft includes:
#include "Geometry/Geometry.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
//#include "RecoAlg/ClusterParamsAlg.h"
#include "RecoAlg/ClusterMergeAlg.h"



#ifndef CLUSTERCRAWLERSHOWER_H
#define CLUSTERCRAWLERSHOWER_H


namespace cluster {
   
  class ClusterCrawlerShower : public art::EDProducer {
    
  public:
    
    explicit ClusterCrawlerShower(fhicl::ParameterSet const& pset); 
    ~ClusterCrawlerShower();
    
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
        

    unsigned int  fMinHitListSize;	 
    double	  fMaxDistance;	 


    //ClusterParamsAlg fCParAlg;
    ClusterMergeAlg fCMergeAlg;
   
    std::string     fClusterCrawlerModuleLabel;
  
  
  protected: 
    
  }; // class ClusterCrawlerShower

}

#endif // CLUSTERCRAWLERSHOWER_H



//namespace cluster{

  //-------------------------------------------------
cluster::ClusterCrawlerShower::ClusterCrawlerShower(fhicl::ParameterSet const& pset) 
//: fCParAlg(pset.get<fhicl::ParameterSet >("ClusterParamsAlg"),pset.get< std::string > ("module_type"))
  : fCMergeAlg(pset.get<fhicl::ParameterSet >("ClusterMergeAlg"))
  , fClusterCrawlerModuleLabel(pset.get<std::string>("ClusterCrawlerModuleLabel"))
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
}


  void cluster::ClusterCrawlerShower::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fClusterCrawlerModuleLabel	=pset.get< std::string >("ClusterCrawlerModuleLabel");
    //fCParAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"));
    fMinHitListSize	  	=pset.get<unsigned int >("MinHitListSize");
    fMaxDistance		=pset.get<double>("MaxDistance");
  }

  //-------------------------------------------------
  cluster::ClusterCrawlerShower::~ClusterCrawlerShower()
  {
  }

  //-------------------------------------------------
  void cluster::ClusterCrawlerShower::beginJob()
  {
    //this doesn't do anything now, but it might someday
  }
    
  //------------------------------------------------------------------------------------//
  void cluster::ClusterCrawlerShower::produce(art::Event& evt)
  { 
    // Get a Handle for the input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterCrawlerModuleLabel,clusterVecHandle);

    // Get a Handle for the input hits, assuming for the time being that the hits label will be that of Cluster Crawler
    art::Handle< std::vector<recob::Hit> > hitVecHandle;
    evt.getByLabel(fClusterCrawlerModuleLabel,hitVecHandle);

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    art::FindManyP<recob::Hit> fmh(clusterVecHandle, evt, fClusterCrawlerModuleLabel);

    
    // make a map of the geo::View to vectors of art::Ptr<recob::Hit>
    std::map<geo::View_t, std::vector< art::Ptr<recob::Hit> > > viewToHits;
    for(size_t i = 0; i < hitVecHandle->size(); ++i)
      viewToHits[hitVecHandle->at(i).View()].push_back(art::Ptr<recob::Hit>(hitVecHandle, i));


    fCMergeAlg.ClearEventInfo();
    std::vector<recob::Cluster> inputShowerClusters;
    std::vector<std::vector<art::Ptr<recob::Hit> > > inputShowerClusterHits;
    for(unsigned int iClust = 0; iClust < clusterVecHandle->size(); iClust++){

      art::Ptr<recob::Cluster> cl(clusterVecHandle, iClust);
      // Hits that are in the cluster
      std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
      std::vector< art::Ptr<recob::Hit> > showerHitList;


      double wirePitch = geom->WirePitch(cl->View());
      double xyScale  = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
      xyScale        *= detprop->SamplingRate()/wirePitch;

      if(hitlist.size()<=fMinHitListSize)
        continue;


      // Need to find the hits in and around the cluster crawler cluster
      //for(size_t i = 0; i < viewToHits[cl->View()].size(); i++){
      for(auto hitsItr = viewToHits[cl->View()].begin(); hitsItr < viewToHits[cl->View()].end(); hitsItr++){

	//std::cout << (*hitsItr)->StartTime() << std::endl;

	double clIntercept = cl->StartPos().back() - cl->dTdW()*cl->StartPos().front();
	double distance = (TMath::Abs((*hitsItr)->PeakTime()-cl->dTdW()*(double)((*hitsItr)->WireID().Wire)-clIntercept)/(std::sqrt(pow(xyScale*cl->dTdW(),2)+1)));

	// Sum up background hits, use smart distance
	double peakTimePerpMin=-(1/cl->dTdW())*(double)((*hitsItr)->WireID().Wire)+cl->StartPos().back()+(1/cl->dTdW())*(cl->StartPos().front());
	double peakTimePerpMax=-(1/cl->dTdW())*(double)((*hitsItr)->WireID().Wire)+cl->EndPos().back()+(1/cl->dTdW())*(cl->EndPos().front());


	if(distance > 1*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.))
         && distance < 25*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.))){
	  if((cl->dTdW() < 0 && (*hitsItr)->PeakTime() < peakTimePerpMin && (*hitsItr)->PeakTime() > peakTimePerpMax)
            || (cl->dTdW() > 0 && (*hitsItr)->PeakTime() > peakTimePerpMin && (*hitsItr)->PeakTime() < peakTimePerpMax)){
	  showerHitList.push_back(*hitsItr);
	  }
	}
      }

      /*
      double wire_start = cl->StartPos().front();
      double wire_end = cl->EndPos().front();
      double time_start = cl->StartPos().back();
      double time_end = cl->EndPos().back();

      if(!fCParAlg.isShower(cl->dTdW(),
			    wire_start,
			    time_start,
			    wire_end,
			    time_end,
			    showerHitList))
	continue;
      */
      
      //// Store this processed shower-like cluster
      //inputShowerClusters.push_back(std::move(temp));
      //inputShowerClusterHits.push_back(hitlist);

      //// Feed this cluster & hitlist into Merge algorithm
      //fCMergeAlg.AppendClusterInfo(inputShowerClusters.back(),hitlist);

    } // End loop on clusters.



    //evt.put(std::move(SuperClusters));
    //evt.put(std::move(assn));

    return;

  }


//} // end namespace 



namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawlerShower)
  
} // end namespace 

