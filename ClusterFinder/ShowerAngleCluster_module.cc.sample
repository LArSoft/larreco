////////////////////////////////////////////////////////////////////////
/// \file  ShowerAngleCluster_module.cc
/// \brief Create a Cluster with an angle 
///
/// \version $Id: AngleCluster.cxx,v 0.1 19/07/2011 12:45:16 PM  andrzejs $
/// \author andrzej.szelc@yale.edu
///  
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// Create a cluster saving its slope angle

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

#ifndef SHOWERANGLECLUSTER_H
#define SHOWERANGLECLUSTER_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include <utility> // std::move()
#include <vector>
#include <string>



#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"
//#include "RecoAlg/HoughBaseAlg.h"
#include "RecoAlg/ClusterParamsAlg.h"
#include "RecoAlg/ClusterMatchAlg.h"
#include "RecoAlg/ClusterMergeAlg.h"

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>


// Framework includes



#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

//#include "TMatrixD.h"
//#include "TVectorD.h"
//#include "TDecompSVD.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"



// LArSoft includes
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "SimulationBase/MCTruth.h"
#include "Geometry/geo.h"
#include "RecoBase/Cluster.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/PlaneGeo.h"

//#include "SimulationBase/simbase.h"
//#include "RawData/RawDigit.h"
//#include "SummaryData/summary.h"
#include "CLHEP/Random/JamesRandom.h"
#include "Utilities/SeedCreator.h"




namespace cluster {

  class ShowerAngleCluster : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerAngleCluster(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~ShowerAngleCluster();                               /**Destructor*/
    void beginJob();                                     
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Routine that finds the cluster and sets the dTdW of the 2D shower*/
   
      

  private:

    double Get2DAngleForHit( art::Ptr<recob::Hit> starthit,std::vector < art::Ptr < recob::Hit> > hitlist);
    double Get2DAngleForHit( unsigned int wire, double time,std::vector < art::Ptr < recob::Hit> > hitlist);
    void ClearandResizeVectors(unsigned int nClusters);

 
    //HoughBaseAlg fHBAlg; 
    ClusterParamsAlg fCParAlg;
    ClusterMergeAlg fCMergeAlg;
    ClusterMatchAlg fCMatchAlg;
    bool fExternalStartPoints;
 //   double fWiretoCm,fTimetoCm,fWireTimetoCmCm;
    
    std::vector< unsigned int >fNWires;
    double fNTimes;
    
    recob::Cluster MainClusterLoop(art::Ptr<recob::Cluster> inCluster,
				   std::vector< art::Ptr<recob::Hit> > hitlist, 
				   unsigned int iClustInput, 
				   unsigned int iClustOutput); 

    recob::Cluster MergeClusterLoop(std::vector< art::Ptr<recob::Hit> > &hitlist,
                                    unsigned int iClustOutput);
    
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geo;
    util::GeometryUtilities gser;
   
    float fTimeTick; // time sample in us
    float fPresamplings;
    float fDriftVelocity;
     
    int fRun,fEvent,fSubRun;
    bool fForceRightGoing;

    //input parameter labels:
 
    std::string fClusterModuleLabel;
  
    std::vector <double> xangle;       // in cm, cm
   
    std::vector <double> lineslopetest;   // in wire, time
    std::vector <double>  lineinterctest;   
    std::vector<double>  fWireVertex,fTimeVertex;
    std::vector<double>  fWireEnd,fTimeEnd;
    std::vector<double> fVerticalness;
    std::vector< double > fErrors;
    unsigned int fNPlanes; // number of planes  

    TH1F *  fh_omega_single;
    TTree* ftree_cluster;
    bool matchflag;
    unsigned int  fMinHitListSize;	 
  }; // class ShowerAngleCluster

}

#endif // SHOWERANGLECLUSTER_H


////////////////////////////////////////////////// Module definition

// ***************** //

//------------------------------------------------------------------------------
cluster::ShowerAngleCluster::ShowerAngleCluster(fhicl::ParameterSet const& pset)
 :// fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg")),
  fCParAlg(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"),pset.get< std::string >("module_type")),
  fCMergeAlg(pset.get< fhicl::ParameterSet >("ClusterMergeAlg")),
  fCMatchAlg(pset.get< fhicl::ParameterSet >("ClusterMatchAlg"))
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit>  >(); 
  // produces< art::Assns<recob::Cluster, recob::Hit>  >(); 
  //produces< art::Assns<recob::Cluster, recob::Cluster>  >(); 
  produces< std::vector < art::PtrVector <recob::Cluster> >  >();

  // Store std::vector< art::PtrVector<recob::SpacePoint> > > if requested by ClusterMatchAlg
  if(fCMatchAlg.StoreSpacePoints()) {
    produces< std::vector <recob::SpacePoint> >();
    produces< std::vector < art::PtrVector <recob::SpacePoint> > >();
  }
  // Create random number engine needed for PPHT
  createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");

  
}


void cluster::ShowerAngleCluster::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel 		=pset.get< std::string >("ClusterModuleLabel");
  fCParAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"));
  fExternalStartPoints		=pset.get< bool >("ExternalStartPoints");
  fMinHitListSize	  	=pset.get<unsigned int >("MinHitListSize");
 }

// ***************** //
cluster::ShowerAngleCluster::~ShowerAngleCluster()
{
}

//____________________________________________________________________________
 void cluster::ShowerAngleCluster::beginRun(art::Run& /*run*/)
  {
    
  //  //std::cout << "------------- In SHowANgle preBeginRun"<<  larp->Efield() << std::endl;
  //  


    return;
  }




//-----------------------------------------------
// namespace cluster {
// struct SortByWire 
// {
//   bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
//   { return 
//       h1.Wire()->RawDigit()->Channel() < h2.Wire()->RawDigit()->Channel() ;
//   }
// };
// }

// ***************** //
void cluster::ShowerAngleCluster::beginJob()
{

    
  // this will not change on a run per run basis.
  fNPlanes = geo->Nplanes();
 // fWirePitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  fTimeTick=detp->SamplingRate()/1000.; 
 
  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  
   // fNTimes=geo->DetHalfWidth(tpc)*2/(fTimetoCm);
    fNWires.resize(fNPlanes);
    
  for(unsigned int i=0;i<fNPlanes;++i){
   
     fNWires[i]=geo->Nwires(i); //Plane(i,tpc).Nwires();
  }
  

  
  
  
  
  ftree_cluster =tfs->make<TTree>("ShowerAngleCluster","Results");/**All-knowing tree with reconstruction information*/
   fh_omega_single= tfs->make<TH1F>("fh_omega_single","Theta distribution Hit",720,-180., 180.) ;
    
    ftree_cluster->Branch("run",&fRun,"run/I");
    ftree_cluster->Branch("subrun",&fSubRun,"subrun/I");
    ftree_cluster->Branch("event",&fEvent,"event/I");
    ftree_cluster->Branch("nplanes",&fNPlanes,"nplanes/I");
  
    ftree_cluster->Branch("wire_vertex","std::vector<double>", &fWireVertex);
    ftree_cluster->Branch("time_vertex","std::vector<double>", &fTimeVertex);
      
    ftree_cluster->Branch("wire_last","std::vector<double>", &fWireEnd);
    ftree_cluster->Branch("time_last","std::vector<double>", &fTimeEnd);
    
    
}

// ************************************* //
void cluster::ShowerAngleCluster::ClearandResizeVectors(unsigned int nClusters) {
  
  ///////////////

  
   fVerticalness.clear();
    
  
  
  // startflag.clear();
   lineslopetest.clear();
   lineinterctest.clear();
  
   xangle.clear();
   fWireVertex.clear();
   fTimeVertex.clear();
   fWireEnd.clear();
   fTimeEnd.clear();
   xangle.resize(nClusters); 


   fWireVertex.resize(nClusters);
   fTimeVertex.resize(nClusters);
   fWireEnd.resize(nClusters);
   fTimeEnd.resize(nClusters);
    
   fVerticalness.resize(nClusters);
   lineslopetest.resize(nClusters);
   lineinterctest.resize(nClusters);
  // startflag.resize(nClusters);
   
   fErrors.resize(nClusters);
    
  
  
}
  


  
  
  
// ***************** //
void cluster::ShowerAngleCluster::produce(art::Event& evt)
{ 
  

  // make an art::PtrVector of the clusters
  std::unique_ptr<std::vector<recob::Cluster> > ShowerAngleCluster(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
  std::unique_ptr< std::vector < art::PtrVector < recob::Cluster > > > classn(new std::vector < art::PtrVector < recob::Cluster > >);

  std::unique_ptr<std::vector<recob::SpacePoint> > MatchSpacePoint(new std::vector<recob::SpacePoint>);
  std::unique_ptr< std::vector< art::PtrVector < recob::SpacePoint > > > sps_assn(new std::vector < art::PtrVector < recob::SpacePoint > >);

 mf::LogWarning("ShowerAngleCluster") << "In produce module " ; 

  //%%%%% this goes into ana module.
  //Find run, subrun and event number:
  fRun = evt.id().run();
  fSubRun = evt.id().subRun();
  fEvent = evt.id().event();

  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
 

  
  art::Handle< std::vector<art::PtrVector < recob::Cluster> > > clusterCollectionHandle;
  //matchflag=true;  
 // try{
//     art::FindManyP<recob::Cluster>  fmclust= art::FindManyP<recob::Cluster>(clusterListHandle, evt, fClusterModuleLabel);
   
  evt.getByLabel(fClusterModuleLabel,clusterCollectionHandle);

 //    //std::cout << " CLusters associated " << fmclust.at(0).size() << " " <<std::endl; 
    //  }
  //catch(cet::exception e) {
 //     mf::LogWarning("ShowerAngleCluster") << "caught exception \n"
     //                              << e;
   //   matchflag=false;		
     
  //    }
  
  
 
 // art::PtrVector<recob::Cluster> clusters;
 if( clusterCollectionHandle.isValid() ) // received matched clusters
   {
     //std::cout << " cluster Assoc Handle size: " << clusterCollectionHandle->size() << std::endl;
     unsigned int nCollections= clusterCollectionHandle->size();
     std::vector < art::PtrVector<recob::Cluster> >::const_iterator clusterSet = clusterCollectionHandle->begin();
     
     for(unsigned int iCol=0;iCol<nCollections;iCol++)
       {
	 
	 const art::PtrVector<recob::Cluster> pvcluster(*(clusterSet++));
	 auto it = pvcluster.begin();
	 int nClusts = pvcluster.size();
	 ClearandResizeVectors(nClusts);   // need to do this smarter (coutn at start and have an overall index?)
	 for(int iClust = 0; iClust < nClusts; ++iClust) {
	   const art::Ptr<recob::Cluster> pclust(*(it++)); 
	   auto pcoll { pclust };
	   art::FindManyP<recob::Hit> fs( pcoll, evt, fClusterModuleLabel);
	   std::vector< art::Ptr<recob::Hit> > hitlist = fs.at(0);
	   //std::cout << "---------- going into pclust " << pclust->StartPos()[0] << " " << pclust->StartPos()[1] << std::endl;
	   ShowerAngleCluster->push_back(std::move(MainClusterLoop( pclust,hitlist, iClust, ShowerAngleCluster->size())));
	   util::CreateAssn(*this, evt, *(ShowerAngleCluster.get()), hitlist, *(assn.get()));
	 }
	 
	 art::PtrVector < recob::Cluster > cvec;
	 //std::cout << "---------- going into aptr outside: " << (*ShowerAngleCluster).size() << " nClusts " << nClusts << std::endl;	
	 
	 // for(unsigned int ip=0;ip<fNPlanes;ip++)  {
	 for(int ip=(int)((*ShowerAngleCluster).size()-1);ip>=(int)((*ShowerAngleCluster).size()-nClusts) && ip>=0 ;ip--)  {
	   art::ProductID aid = this->getProductID< std::vector < recob::Cluster > >(evt);
	   //std::cout << " ip is equal to: " << ip  << std::endl;
	   art::Ptr< recob::Cluster > aptr(aid, ip, evt.productGetter(aid));
	   ////std::cout << "---------- going into aptr " << aptr->StartPos()[0] << " " << aptr->StartPos()[1] << std::endl;
	   cvec.push_back(aptr);
	 }
	 if((*ShowerAngleCluster).size()>0)   //make sure to make associations to something that exists.
	   classn->push_back(cvec);   
       }
     
   }
 else // no matching was done. Need to do it ourselves, but first run through all of the clusters
   {
     
     evt.getByLabel(fClusterModuleLabel,clusterListHandle);
     art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);
     
     //std::cout << " ++++ Clusters received " << clusterListHandle->size() << " +++++ " << std::endl;
     if(clusterListHandle->size() ==0 )
       {
	 //std::cout << " no clusters received! exiting " << std::endl;
	 return;
       }
     ClearandResizeVectors(clusterListHandle->size());
     // resizing once cluster size is known.
     
     //endflag=false;
     
     
     /** 
	 Loop over input clusters, process & re-evaluate input clusters' start/end/angle 
	 by using ClusterParamsAlg. The results are filtered to only keep shower-like
	 clusters. The output of this loop is then used for merging algorithm.
	 -- Kazu Dec. 27 2013
     */
     fCMergeAlg.ClearEventInfo();
     std::vector<recob::Cluster> inputShowerClusters;
     std::vector<std::vector<art::Ptr<recob::Hit> > > inputShowerClusterHits;
     for(unsigned int iClust = 0; iClust < clusterListHandle->size(); iClust++){
       
       art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
       std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
       
       if(hitlist.size()<=fMinHitListSize )
	 continue;

       recob::Cluster temp = MainClusterLoop( cl, hitlist, iClust, inputShowerClusters.size());
       if(!fCParAlg.isShower(lineslopetest[iClust],
			     fWireVertex[iClust],
			     fTimeVertex[iClust],
			     fWireEnd[iClust],
			     fTimeEnd[iClust], 
			     hitlist))
	 continue;
       
       // Store this processed shower-like cluster
       inputShowerClusters.push_back(std::move(temp));
       inputShowerClusterHits.push_back(hitlist);

       // Feed this cluster & hitlist into Merge algorithm
       fCMergeAlg.AppendClusterInfo(inputShowerClusters.back(),hitlist);

     } // End loop on clusters.
     
     // Run merging algorithm & retrieve the result
     fCMergeAlg.ProcessMergeAlg();

     std::vector<std::vector<unsigned int> > mergedClusterSets = fCMergeAlg.GetClusterSets();
     
     /**
	Loop over the merging algorithm result. For clusters to be merged,
	re-compute cluster parameters (start,end,angle). For clusters that
	are not to be merged, copy already-computed parameters and only
	modify the cluster ID to an appropriate value.
     */
     std::vector<std::vector<art::Ptr<recob::Hit> > >outputShowerClusterHits;
     outputShowerClusterHits.reserve(mergedClusterSets.size());
     ShowerAngleCluster->reserve(mergedClusterSets.size());
     for(auto const& cluster_set : mergedClusterSets) {

       int outClusterID = -1;

       // If length == 1, no merging required
       if(cluster_set.empty())

	 mf::LogError(__FUNCTION__) << "Encountered 0-length merged cluster sets (logic error)!";
       
       else if(cluster_set.size()==1) {

	 unsigned int iClust = cluster_set[0];
	 outClusterID = (int)(ShowerAngleCluster->size());
	 ShowerAngleCluster->push_back(recob::Cluster(inputShowerClusters[iClust].StartPos()[0],
						      inputShowerClusters[iClust].SigmaStartPos()[0],
						      inputShowerClusters[iClust].StartPos()[1],
						      inputShowerClusters[iClust].SigmaStartPos()[1],
						      inputShowerClusters[iClust].EndPos()[0],
						      inputShowerClusters[iClust].SigmaEndPos()[0],
						      inputShowerClusters[iClust].EndPos()[1],
						      inputShowerClusters[iClust].SigmaEndPos()[1],
						      inputShowerClusters[iClust].dTdW(),
						      inputShowerClusters[iClust].SigmadTdW(),
						      inputShowerClusters[iClust].dQdW(),
						      inputShowerClusters[iClust].SigmadQdW(),
						      inputShowerClusters[iClust].Charge(),
						      inputShowerClusters[iClust].View(),
						      outClusterID,
						      inputShowerClusters[iClust].Plane()));
	 outputShowerClusterHits.push_back(inputShowerClusterHits.at(iClust));
       }
       // else we need to merge hits & re-evaluate cluster parameters.
       else{

	 std::vector<art::Ptr<recob::Hit> > hitlist;
	 for(auto const& iClust : cluster_set)

	   for(auto const& iHit : inputShowerClusterHits[iClust])

	     hitlist.push_back(iHit);
	 outClusterID = (int)(ShowerAngleCluster->size());
	 ShowerAngleCluster->push_back(std::move(MergeClusterLoop(hitlist,outClusterID)));
	 outputShowerClusterHits.push_back(hitlist);
       }

       util::CreateAssn(*this, evt, *(ShowerAngleCluster.get()), outputShowerClusterHits.at(outClusterID), *(assn.get()));

     }

     if(!matchflag)
       {
	 //matching code here.
	 //// declare,  matching table
	 std::vector< std::vector < unsigned int > > ClusterSets;
	 
	 //
	 // Following code segment added by Kazu to use a separate module
	 // to find possible cluster combination pairs and fill ClusterSets variable.
	 // -- Kazu on Nov. 30 2013
	 //
	 // Update:
	 // Modified Matching code to work with computed ShowerAngleCluster output recob::Cluster.
	 // -- Kazu on Dec. 27 2013
	 //
         fCMatchAlg.ClearEventInfo();
         fCMatchAlg.SetMCTruthModName("generator"); // Won't crash even if data product do not exist.
         fCMatchAlg.FillMCInfo(evt);
	 for(size_t i=0; i<ShowerAngleCluster->size(); ++i )
	  
	   fCMatchAlg.AppendClusterInfo(ShowerAngleCluster->at(i),
					outputShowerClusterHits.at(i));

	 if(fNPlanes==2)
	   fCMatchAlg.MatchTwoPlanes();
	 else if(fNPlanes==3)
	   fCMatchAlg.MatchThreePlanes();
	 else
	   mf::LogError("ShowerAngleCluster")<<"Matching with NPlanes <2 or >3 not supported!";
	 
	 // Retrieve the result
	 ClusterSets = fCMatchAlg.GetMatchedClusters();

	 // Retrieve computed SpacePoint (and store if configured to do so)
	 std::vector<std::vector<recob::SpacePoint> > matchSPS_v = fCMatchAlg.GetMatchedSpacePoints();
	 
	 if(fCMatchAlg.StoreSpacePoints()){
	   
	   for(auto const& spsVector : matchSPS_v)
	     
	     for(auto const& sps : spsVector)
	       
	       MatchSpacePoint->push_back(sps);
	   
	 }
	 
	 
	 //
	 // End of Kazu's added code on Nov. 30 2013
	 //
	 
	 // 	std::vector < int > maximums;
	 // 	std::vector < unsigned int > maxclusts;
	 // 	maximums.resize(fNPlanes);
	 // 	maxclusts.resize(fNPlanes,0);
	 // 	
	 // 	for(unsigned int iClust = 0; iClust < clusterListHandle->size(); iClust++){
	 // 	  
	 // 	  art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
	 // 	  std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
	 // 	  if( (int)hitlist.size() > maximums[cl->View()])  
	 // 	   {
	 // 	     maximums[cl->View()]=hitlist.size();
	 // 	     maxclusts[cl->View()]=iClust;
	 // 	   }
	 // 	}
	 
	 //
	 // Now that ClusterSets are filled, I comment out the following section
	 // -- Kazu on Nov. 30 2013
	 //
	 /*
	 //sanity check:
	 for(unsigned int ii=0;ii<maxclusts.size();ii++)
	 //std::cout << "ii" << ii << " " << maxclusts[ii] << std::endl;
	 ClusterSets.push_back(maxclusts);  
	 
	 */ //End of commented out block by Kazu on Nov. 30 2013
	 
	 //  for(unsigned int ii=0;ii<clusterListHandle->size();ii++){
	 //for(unsigned int jj=ii+1;jj<clusterListHandle->size();jj++){
	 //  // don't try to match clusters in the same view.
	 //  if(ShowerAngleCluster[ii]->View()==ShowerAngleCluster[jj]->View())
	 //   continue;
	 
	 // std::vector<double > xpos1=ShowerAngleCluster[ii]->StartPos();
	 // std::vector<double > xpos2=ShowerAngleCluster[jj]->StartPos();
	 
	 
	 //} //end jj loop
	 // } // end ii loop
	 
	 //
	 
	 //
	 // Tip: ClusterSets is std::vector<std::vector<unsigned int> > of length J x K
	 // where J is the # of wire planes and K is the length of matched vector in each plane.
	 // ShowerReco code prefers a different index convention: a 1-D vector of consecutive
	 // Cluster index ID, which totally makes sense for a conventional purpose. 
	 // So we do a little bit unusual looping here over ClusterSets.
	 // -- Kazu 03 Dec 2013
	 //
	 art::PtrVector < recob::Cluster > cvec;
	 art::PtrVector < recob::SpacePoint > spsvec;
	 for(unsigned int iSet=0; iSet<ClusterSets[geo::kU].size(); iSet++) {
	   cvec.clear();
	   spsvec.clear();
	   for(unsigned int ip=0; ip<fNPlanes; ++ip) {

	     // Obtain a poitner to the output cluster through an index
	     art::ProductID aid = this->getProductID< std::vector < recob::Cluster > >(evt);
	     art::Ptr< recob::Cluster > aptr(aid, ClusterSets[ip][iSet], evt.productGetter(aid));
	     cvec.push_back(aptr);
	   }
	   // Do the same for SpacePoint if MatchAlg requests to store it
	   if(cvec.size() && fCMatchAlg.StoreSpacePoints()){
	     size_t index_offset = 0;
	     for(size_t i=0; i<iSet; ++i)
	       index_offset += matchSPS_v[i].size();
	     for(size_t i=0; i<matchSPS_v[i].size(); ++i) {
	       art::ProductID aid = this->getProductID< std::vector < recob::SpacePoint > >(evt);
	       art::Ptr< recob::SpacePoint > aptr(aid, (unsigned int)(index_offset+i), evt.productGetter(aid));
	       spsvec.push_back(aptr);
	     }
	   }
	   
	   //make sure to make associations to something that exists.  
	   if(ClusterSets[geo::kU].size()>0 && cvec.size()) {
	     classn->push_back(cvec);		
	     if(spsvec.size())
	       sps_assn->push_back(spsvec);
	   }
	 } // end of loop on sets
	 //At least one of the first two planes 	
       } // end if matchflag is false
     /**Fill the output tree with all information */
     ftree_cluster->Fill();
     
   } // end else if collection is valid
 //  //std::cout << "--- classn size, saving: " << classn.size() << std::endl << std::endl;
 
 evt.put(std::move(ShowerAngleCluster));
 evt.put(std::move(assn));
 evt.put(std::move(classn));
 if(fCMatchAlg.StoreSpacePoints()) {
   evt.put(std::move(MatchSpacePoint));
   evt.put(std::move(sps_assn));
 }
 
}




recob::Cluster cluster::ShowerAngleCluster::MainClusterLoop( art::Ptr<recob::Cluster> inCluster,
							     std::vector< art::Ptr<recob::Hit> > hitlist, 
							     unsigned int iClustInput,
							     unsigned int iClustOutput) {

      std::vector< double > spos=inCluster->StartPos();
      std::vector< double > sposerr=inCluster->SigmaStartPos();
      
      std::vector< double > epos=inCluster->EndPos();
      std::vector< double > eposerr=inCluster->SigmaEndPos();
      
     // Start positions are determined elsewhere and accepted here
     
	
      ///! Change for accepting DBCluster and cheatcluster, so that it doesn't get fooled.
      //startflag.push_back(false); 
      
    //std::cout << " hitlist size: " << hitlist.size() << std::endl;
   

      
    double lineslope, lineintercept,goodness,wire_start,time_start,wire_end,time_end;
    int nofshowerclusters=0;
   // //std::cout << "++++ hitlist size " << hitlist.size() << std::endl;
    
    //////////////////////////////////
   // fCParAlg.Find2DAxisRoughHighCharge(lineslope,lineintercept,goodness,hitlist);
   //  //std::cout << "%%%%%%%% lineslope, intercept " << lineslope << " "<< lineintercept << std::endl;

    
    
    fCParAlg.Find2DAxisRough(lineslope,lineintercept,goodness,hitlist);
    fVerticalness[iClustInput]=goodness;
    //if(hitlist_high.size()<=3 )
	//continue;
    
    fCParAlg.Find2DStartPointsHighCharge( hitlist,wire_start,time_start,wire_end,time_end);
    
    fWireVertex[iClustInput]=wire_start;
    fTimeVertex[iClustInput]=time_start;
   

    
    double wstn=0,tstn=0,wendn=0,tendn=0;
    fCParAlg.FindTrunk(hitlist,wstn,tstn,wendn,tendn,lineslope,lineintercept);
    int fDirection = (wstn<wendn)  ? 1 : -1 ;     // if current starting point is less then end point then direction is 1.
    
    
    double HiBin,LowBin,invHiBin,invLowBin;
    fCParAlg.FindDirectionWeights(lineslope,wstn,tstn,wendn,tendn,hitlist,HiBin,LowBin,invHiBin,invLowBin); 
    
    if(invHiBin+invLowBin> 1000)
      nofshowerclusters++;
    
    fDirection=fCParAlg.DecideClusterDirection(hitlist,lineslope,wstn,tstn,wendn,tendn);
     //std::cout << "%%%%%%%% direction start points: (" << wstn<<","<<tstn<<"), end: ( "<<wendn << ","<<tendn <<")" << "Direction: " << fDirection << std::endl;
    wire_start=wstn;
    time_start=tstn;
    wire_end=wendn;
    time_end=tendn;
    fCParAlg.RefineStartPointsHough(hitlist, wire_start,time_start,wire_end,time_end,fDirection); 
     //std::cout << "%%%%%%%% Hough line refine start points: (" << wire_start<<","<<time_start<<"), end: ( "<<wire_end << ","<<time_end <<")" << "Direction: " << fDirection << std::endl; 
    if(fExternalStartPoints==false)  // start points have been handed over from external
    {
      fWireVertex[iClustInput]=wire_start;
      fTimeVertex[iClustInput]=time_start;
      fWireEnd[iClustInput]=wire_end;
      fTimeEnd[iClustInput]=time_end; 
    }
    else
    {
     fWireVertex[iClustInput]=spos[0];
     fTimeVertex[iClustInput]=spos[1]; 
     fWireEnd[iClustInput]=epos[0];
     fTimeEnd[iClustInput]=epos[1];
    }
    
    
    
    lineslopetest[iClustInput]=lineslope; 
    lineinterctest[iClustInput]=lineintercept;
      
	  
    xangle[iClustInput]=Get2DAngleForHit( fWireVertex[iClustInput],fTimeVertex[iClustInput], hitlist);
  
  
    //std::cout << " in save loop, \"iplane\" " << iClustInput <<  std::endl;
    //std::cout << " hitlist size: " << hitlist.size() <<  std::endl;
     
  
    if(fVerticalness[iClustInput]<1) 
      fErrors[iClustInput]=0.1;
    else
      fErrors[iClustInput]=10; 
  
      
    // Following 2 lines commented out: we have an access to the plane through viewfix variable
    //unsigned int xplane;
    //xplane=hitlist[0]->WireID().Plane;

    geo::View_t viewfix = hitlist[0]->View();
    
    //std::cout << " wire " << fWireVertex[iClustInput] << std::endl;
    //std::cout << " fErrors " << fErrors[iClustInput] << std::endl;
    //std::cout << " time " << fTimeVertex[iClustInput] << std::endl;
    //std::cout << " wirelast " << fWireEnd[iClustInput] << std::endl;
    //std::cout << " timelast " <<  fTimeEnd[iClustInput] << std::endl; 
    //std::cout << " xangle " <<  xangle[iClustInput] << std::endl; 
    //std::cout <<"lineslope " <<  lineslopetest[iClustInput] << std::endl; 
    //std::cout <<"lineinterc " << lineinterctest[iClustInput] <<std::endl;	
    //std::cout <<"plane  " << xplane <<std::endl;	

    return recob::Cluster /* outCluster */ (
			      fWireVertex[iClustInput], fErrors[iClustInput],
			      fTimeVertex[iClustInput], fErrors[iClustInput],
			      fWireEnd[iClustInput], fWireEnd[iClustInput]*0.05,
			      fTimeEnd[iClustInput], fTimeEnd[iClustInput]*0.05,  
			      xangle[iClustInput], xangle[iClustInput]*0.05, lineslopetest[iClustInput],lineinterctest[iClustInput],5.,
			      viewfix,
			      iClustOutput,
			      hitlist.front()->WireID().planeID());

//     ShowerAngleCluster->push_back(temp);
//   
//     util::CreateAssn(*this, evt, *(ShowerAngleCluster.get()), hitlist, *(assn.get()));
//     return outCluster;
    
} // cluster::ShowerAngleCluster::MainClusterLoop()

recob::Cluster cluster::ShowerAngleCluster::MergeClusterLoop( std::vector<art::Ptr<recob::Hit> > &hitlist,
                                                              unsigned int iClustOutput) {

  double lineslope, lineintercept,goodness,wire_start,time_start,wire_end,time_end;

  // initialize
  lineslope = lineintercept = goodness = wire_start = time_start = wire_end = time_end = 0;

  fCParAlg.Find2DAxisRough(lineslope,lineintercept,goodness,hitlist);

  fCParAlg.Find2DStartPointsHighCharge( hitlist, wire_start, time_start, wire_end, time_end);

  fCParAlg.FindTrunk(hitlist,wire_start,time_start,wire_end,time_end,lineslope,lineintercept);
  int fDirection = (wire_start<wire_end)  ? 1 : -1 ;     // if current starting point is less then end point then direction is 1.

  double HiBin,LowBin,invHiBin,invLowBin;
  fCParAlg.FindDirectionWeights(lineslope,
                                wire_start, time_start,
                                wire_end, time_end,
                                hitlist,
                                HiBin,LowBin,invHiBin,invLowBin);

  fDirection=fCParAlg.DecideClusterDirection(hitlist,
                                             lineslope,
                                             wire_start, time_start,
                                             wire_end,   time_end);

  fCParAlg.RefineStartPointsHough(hitlist, wire_start,time_start,wire_end,time_end,fDirection);

  if(fExternalStartPoints)  // start points have been handed over from external
    {
      // For merging, "external" start/end points has to be re-generated here because
      // the input clusters' hit list has been merged unless merging algorithm produces
      // a relevant start/end point on the merged cluster
      mf::LogError("ShowerAngleCluster")<<"External start points not implemented for merging stage!"<<std::endl;
    }

  double angle_2d = Get2DAngleForHit(wire_start, time_start, hitlist);

  double vtx_error = ( (goodness < 1) ? 0.1 : 10);

  geo::View_t viewfix = hitlist[0]->View();

  recob::Cluster outCluster(wire_start, vtx_error,
                            time_start, vtx_error,
                            wire_end, wire_end*0.05,
                            time_end, time_end*0.05,
                            angle_2d, angle_2d*0.05, lineslope, lineintercept, 5.,
                            viewfix,
                            iClustOutput,
                            hitlist.front()->WireID().planeID()
                            );
  return outCluster;
}

////////////////////////////////////////////////////////////////////////////////
// Method to get the 2D angle ogf a Cluster based on its starting wire and time.
////////////////////////////////////////////////////////////////////////////////

double cluster::ShowerAngleCluster::Get2DAngleForHit( unsigned int swire,double stime,std::vector < art::Ptr < recob::Hit> > hitlist) {
  
  fh_omega_single->Reset();
  
  unsigned int wire;
  // this should changed on the loop on the cluster of the shower
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
     // art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = (*hitIter)->PeakTime();  
      wire=(*hitIter)->WireID().Wire; 
      double omx=gser.Get2Dangle((double)wire,(double)swire,time,stime);
      fh_omega_single->Fill(180*omx/TMath::Pi(),(*hitIter)->Charge());
     }
    
  double omega = fh_omega_single->GetBinCenter(fh_omega_single->GetMaximumBin());// Mean value of the fit
   
  return omega; // in degrees.
}






namespace cluster {

  DEFINE_ART_MODULE(ShowerAngleCluster)

}
