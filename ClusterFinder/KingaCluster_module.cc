////////////////////////////////////////////////////////////////////////
// $Id: DBSCANfinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// \file HoughLineFinder_module.cc
//
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////


extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <dirent.h>
#include "TMath.h"
#include <string>

// ROOT includes
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

// Framework includes
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes 
#include "RawData/RawDigit.h"
#include "Geometry/Geometry.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"



//#ifndef KingaCluster_H
//#define KingaCluster_H


class TH1F;
class TH2F;

namespace recob { class Hit; }
namespace cluster {
   
  class KingaCluster : public art::EDProducer {
    
  public:
    
    explicit KingaCluster(fhicl::ParameterSet const& pset); 
    KingaCluster();
    ~KingaCluster();
         
    void produce(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void AngularDistribution(unsigned int cstat, 
			     unsigned int tpc, 
			     unsigned int plane);
    void FindMax(unsigned int cstat, 
		 unsigned int tpc, 
		 unsigned int plane);
    //void FinalPeaks();
    void FitAngularDistributions();
    void FindClusters(unsigned int cstat, 
		      unsigned int tpc, 
		      unsigned int plane);
    void ReassignHitID(unsigned int cstat, 
		       unsigned int tpc, 
		       unsigned int plane,
		       unsigned int HitPosition,
		       unsigned int WrongPeakNo);
 
  private:
  
    std::string fDBScanModuleLabel;  
    //std::string fGenieGenModuleLabel;  
    std::string fEndPoint2DModuleLabel;
   
    TH1F *fh_theta_ind;
    TH1F *fh_theta_coll;
    TH1F *fh_theta_ind_2D;
    TH1F *fh_theta_coll_2D;
    TH1F *fh_theta_coll_Area;
    TH1F *fh_theta_ind_Area;
    TH1F *Hit_Area_Ind;
    TH1F *Hit_Area_Coll;
   
    art::ServiceHandle<geo::Geometry> fGeom;

    std::vector< art::Ptr<recob::Hit> > allhits;
    std::vector<int>           maxBin;             ///< stores bin # of local maximum
    std::vector<int>           MaxStartPoint;      ///< bin no of the starting point of a peak
    std::vector<int>           MaxEndPoint;        ///< bin no of the end point of a peak
    std::vector<int>           MaxStartPointTheta; ///< theta value of the starting point of a peak
    std::vector<int>           MaxEndPointTheta;   ///< theta value of the end point of a peak
    std::vector<int>           fwire_vertex;
    std::vector<int>           fwire_vertex_reco;
    std::vector<double>        ftime_vertex;
    std::vector<double>        ftime_vertex_reco;
    std::vector<unsigned int>  HitsWithClusterID;
    double                     ftimetick;          ///< get from parameterset
    double                     fdriftvelocity;     ///< get from paramtereset 9either k and V)
    double                     fpi;
    //int                        fMC;
    //double                     MCvertex [3];

    std::vector<double>        maxBinValues;
    std::vector<double>        OriginalmaxBinValues;
    std::vector<int>           SortedMaxBin;
    std::vector<int>           FinalPeaks;
    int                        fpeaks_found;            ///< flag to determine whether the program should continue or not
    int                        fpeak_problem;           ///< flag to determine whether the program should continue or not
    bool                       need_to_reassign_hitsIDs;
    bool                       go_ahead_at_reassign;
  protected:

    
  };
  

}

//#endif // KingaCluster_H

namespace cluster{

  //-------------------------------------------------
  KingaCluster::KingaCluster(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  }
  
  //-------------------------------------------------
  KingaCluster::~KingaCluster()
  {
  }
  
  //-------------------------------------------------
  void KingaCluster::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDBScanModuleLabel     = pset.get< std::string >("DBScanModuleLabel");
    fEndPoint2DModuleLabel = pset.get< std::string >("EndPoint2DModuleLabel");
  
    return;
  }
  
  //-------------------------------------------------
  void KingaCluster::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
   
    fh_theta_ind       = tfs->make<TH1F>("fh_theta_ind","theta angle in degrees, Induction Plane", 180,-180 ,180  );
    fh_theta_coll      = tfs->make<TH1F>("fh_theta_coll","theta angle in degrees, Collection Plane", 180,-180 ,180  );
    fh_theta_ind_2D    = tfs->make<TH1F>("fh_theta_ind_2D","theta angle in degrees, Induction Plane", 180,-180 ,180  );
    fh_theta_coll_2D   = tfs->make<TH1F>("fh_theta_coll_2D","theta angle in degrees, Collection Plane", 180,-180 ,180  );
    fh_theta_ind_Area  = tfs->make<TH1F>("fh_theta_ind_Area","Hit Area vs theta angle in degrees, Induction Plane", 180,-180 ,180  );
    fh_theta_coll_Area = tfs->make<TH1F>("fh_theta_coll_Area","Hit Area vs theta angle in degrees, Collection Plane", 180,-180 ,180  );
  
    Hit_Area_Ind       = tfs->make<TH1F>("Hit_Area_Ind","Hit Area, Induction Plane", 100,0 ,1  );
    Hit_Area_Coll      = tfs->make<TH1F>("Hit_Area_Coll","Hit Area, Collection Plane", 100,0 ,1  );
  
  }
  
  //-----------------------------------------------------------------
  namespace cluster {
    struct SortByWire {
      bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
      { 
        return h1.Wire()->Channel() < h2.Wire()->Channel() ;
      }
    };
  }
  
  //-----------------------------------------------------------------
  void KingaCluster::produce(art::Event& evt)
  {
    
    mf::LogInfo("KingaCluster") <<"In KingaCluster::produce(art::Event& evt)";
    mf::LogInfo("KingaCluster") <<" Working on " << "Run: " << evt.run() << " Event: " << evt.id().event();
  
    fpeak_problem=0;
    fpeaks_found=1;
  
    ftime_vertex_reco.clear();
    fwire_vertex_reco.clear();
  
    //////////////////////////////////////////////////////
    // here is how to get a collection of objects out of the file
    // and connect it to an art::Handle
    //////////////////////////////////////////////////////
    // Read in the clusterList object(s).
     
    art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
    evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
    art::PtrVector<recob::EndPoint2D> endpointlist;
   
    for (size_t i = 0; i < endpointListHandle->size(); ++i){
      art::Ptr<recob::EndPoint2D> endpointHolder(endpointListHandle,i);
      endpointlist.push_back(endpointHolder);
    }
      
    mf::LogInfo("KingaCluster")<<"SHOULD BE GETTING RECO VERTEX, endpointlist.size()= "<<endpointlist.size();
    
    if(endpointlist.size()==0){
      mf::LogInfo("KingaCluster")<<"ATTENTION: NO VERTEX FOUND, KINGACLUSTER WILL EXIT";
      return;
    }
  
    for (size_t j = 0; j < endpointlist.size(); ++j){
  
      ftime_vertex_reco.push_back(endpointlist[j]->DriftTime());
      fwire_vertex_reco.push_back(endpointlist[j]->WireID().Wire);
            
      mf::LogInfo("KingaCluster") << "j=" << j 
  				<< " vtx2d_w=" <<endpointlist[j]->WireID().Wire
  				<< " vtx2d_t=" <<endpointlist[j]->DriftTime();
            
    }
        
       
    //********************************************************************     
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  
    //Point to a collection of clusters to output.
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
  
    art::PtrVector<recob::Cluster> clusIn;
   
    std::vector< art::Ptr<recob::Hit> > hits;
    art::PtrVector<recob::Hit> clusterHits;
    
    for(size_t ii = 0; ii < clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
  
    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fDBScanModuleLabel);
    mf::LogInfo("KingaCluster")<<"No of DBSCAN clusters= "<<clusIn.size();
  
        // make a map of the geo::PlaneID to vectors of art::Ptr<recob::Hit>
    std::map<geo::PlaneID, std::vector< art::Ptr<recob::Hit> > > planeIDToHits;

    for(size_t j = 0; j < clusIn.size(); ++j) {
      
      hits = fmh.at(j);
      
      for(size_t i = 0; i < hits.size(); ++i)
	planeIDToHits[hits.at(i)->WireID().planeID()].push_back(hits.at(i));
    }


    for(auto & itr : planeIDToHits){
      geo::PlaneID planeID = itr.first;
      need_to_reassign_hitsIDs=0;
      go_ahead_at_reassign=0;
      
      allhits.resize(itr.second.size());
      allhits.swap(itr.second);

      //mf::LogInfo("KingaCluster")<<"allhits.size()="<<allhits.size();
      
      //Now we have hits for the plane that we are on right now, so let's do some work:
      //maxBin.clear();
      mf::LogInfo("KingaCluster") << "ATTENTION, STARTING WORK ON PLANE# " << planeID.Plane;
  
      AngularDistribution(planeID.Cryostat, planeID.TPC, planeID.Plane);
  
      FindMax(planeID.Cryostat, planeID.TPC, planeID.Plane);
  
      if(fpeak_problem == 1){
	mf::LogInfo("KingaCluster") << "bye bye";
  	  
	allhits.clear();
	maxBin.clear();
	maxBinValues.clear();
	SortedMaxBin.clear();
	MaxStartPoint.clear();
	MaxEndPoint.clear();
	MaxStartPointTheta.clear();
	MaxEndPointTheta.clear();
	HitsWithClusterID.clear();
	FinalPeaks.clear();
	OriginalmaxBinValues.clear();
  	
	clusterHits.clear();
  	
	for(int bin = 0; bin < fh_theta_ind_2D->GetNbinsX(); ++bin){
	  
	  fh_theta_ind_2D->SetBinContent(bin,0);
	  fh_theta_coll_2D->SetBinContent(bin,0);
	  fh_theta_ind->SetBinContent(bin,0);
	  fh_theta_coll->SetBinContent(bin,0);
	  fh_theta_coll_Area->SetBinContent(bin,0);
	  fh_theta_ind_Area->SetBinContent(bin,0);
	}
	return;
      }// end if peak problem
  
      if(fpeaks_found == 0){
	mf::LogInfo("KingaCluster") << "KingaClusters FAILED on this event because "
				    << "no peaks were found. Perhaps your threshold "
				    << "for peak's height is too big. Goodbye! ";
	allhits.clear();
	maxBin.clear();
	maxBinValues.clear();
	SortedMaxBin.clear();
	MaxStartPoint.clear();
	MaxEndPoint.clear();
	MaxStartPointTheta.clear();
	MaxEndPointTheta.clear();
	HitsWithClusterID.clear();
	FinalPeaks.clear();
	OriginalmaxBinValues.clear();
	for(int bin = 0; bin< fh_theta_ind_2D->GetNbinsX(); ++bin){
	  
	  fh_theta_ind_2D->SetBinContent(bin,0);
	  fh_theta_coll_2D->SetBinContent(bin,0);
	  fh_theta_ind->SetBinContent(bin,0);
	  fh_theta_coll->SetBinContent(bin,0);
	  fh_theta_coll_Area->SetBinContent(bin,0);
	  fh_theta_ind_Area->SetBinContent(bin,0);
	}
  	
	return;
      }
      
      //FinalPeaks();
      
      FindClusters(planeID.Cryostat, planeID.TPC, planeID.Plane);
      if(need_to_reassign_hitsIDs==1){
	mf::LogInfo("KingaCluster") <<"***************************************************************\n"
				    <<"***************  ATTENTION   ***********************\n"
				    <<" WILL NEED TO REASSIGN HIT IDs\n"
				    <<"***************************************************************\n";
	FindClusters(planeID.Cryostat, planeID.TPC, planeID.Plane);
      }
      mf::LogInfo("KingaCluster")<<"HitsWithClusterID.size()= "<< HitsWithClusterID.size()
				 << "compare with allhits.size()= "<< allhits.size();
  	
  	
      //********************************************************************         
      for(size_t ClusterNo = 0; ClusterNo < MaxStartPoint.size(); ++ClusterNo) {
	
	double totalQ = 0.;
	for(size_t j = 0; j < HitsWithClusterID.size(); ++j){
	  
	  if(HitsWithClusterID[j] == (ClusterNo+1)){
	    clusterHits.push_back(allhits[j]);
	    totalQ += clusterHits.back()->Charge();
	  } //if
  	  
	} //loop over HitsWithClusterID
  	
  	  // let's look at the clusters produced:
	mf::LogInfo("KingaCluster")<<"For Cluster # "<<ClusterNo<<" we have "
				   <<clusterHits.size()<<" hits :";
	
	//.................................
	if (clusterHits.size() > 0){
	  
	  /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	  unsigned int p = 0; 
	  unsigned int t = 0; 
	  unsigned int c = 0; 
	  unsigned int sw = clusterHits[0]->WireID().Wire;
	  unsigned int ew = clusterHits[clusterHits.size()-1]->WireID().Wire;
  	  
	  clusterHits.sort(cluster::SortByWire());
  	  
	  recob::Cluster cluster(sw*1., 0.,
				 clusterHits[0]->PeakTime(), 
				 clusterHits[0]->SigmaPeakTime(),
				 ew*1., 0.,
				 clusterHits[clusterHits.size()-1]->PeakTime(), 
				 clusterHits[clusterHits.size()-1]->SigmaPeakTime(),
				 -999., 0., 
				 -999., 0.,
				 totalQ,
				 fGeom->Cryostat(c).TPC(t).Plane(p).View(),
				 ClusterNo);
	  
	  ccol->push_back(cluster);
  	  
	  // associate the hits to this cluster
	  util::CreateAssn(*this, evt, *(ccol.get()), clusterHits, *(assn.get()));
  	  
	  mf::LogInfo("KingaCluster") << "Produced Cluster #" << ClusterNo;
	  //mf::LogInfo("KingaCluster")<<"no of hits for this cluster is "<<clusterHits.size()<<std::endl;
	  // mf::LogInfo("KingaCluster")<<cluster.StartPos()[0]<<", "<<cluster.StartPos()[1]<<" --> "
	  //<<cluster.EndPos()[0]<<", "<<cluster.EndPos()[1]<<std::endl;
	  clusterHits.clear();
	  //////
	} // end if hits in current cluster
  	
      } //clusters
      
      allhits.clear();
      maxBin.clear();
      maxBinValues.clear();
      SortedMaxBin.clear();
      MaxStartPoint.clear();
      MaxEndPoint.clear();
      MaxStartPointTheta.clear();
      MaxEndPointTheta.clear();
      HitsWithClusterID.clear();
      FinalPeaks.clear();
      OriginalmaxBinValues.clear();
      mf::LogInfo("KingaCluster")<<"Should be starting to work on the other plane now";
    }//PlaneIDs
   
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
      
    return;
  }
      
  //..............................................................  
  
  void KingaCluster::AngularDistribution(unsigned int cstat, 
					 unsigned int tpc, 
					 unsigned int plane){   
    
    // get the signal type for the passed in plane
    geo::SigType_t sigType = fGeom->Plane(plane, tpc, cstat).SignalType();
  
    if(sigType == geo::kInduction){
      for(int bin = 0; bin <= fh_theta_ind->GetNbinsX(); ++bin){
        
        fh_theta_ind_2D->SetBinContent(bin,0);
        fh_theta_ind->SetBinContent(bin,0);
        fh_theta_ind_Area->SetBinContent(bin,0);
      }
    }
     
    if(sigType == geo::kCollection){
      for(int bin = 0; bin <= fh_theta_ind_Area->GetNbinsX(); ++bin){
        
        fh_theta_coll_2D->SetBinContent(bin,0);
        fh_theta_coll->SetBinContent(bin,0);
        fh_theta_coll_Area->SetBinContent(bin,0);
        
      }
    }
     
     
    //std::vector<unsigned int> fwire_vertex,ftime_vertex;
    double a_polar, b_polar,theta_polar;
    ftimetick      =  0.198; //get from parameterset
    fdriftvelocity =  0.157;  //get from paramtereset 9either k and V)
    fpi            = TMath::Pi();
        
    //   mf::LogInfo("KingaCluster")<<"for PLANE "<<plane
    // 			     <<" fwire_vertex= "<<fwire_vertex[plane]
    // 			     <<" ftime_vertex= "<<ftime_vertex[plane];
     
     
    mf::LogInfo("KingaCluster")<<"No of HITS for plane "<<plane<<" is: "<<allhits.size();
     
    for(size_t i = 0; i < allhits.size(); ++i){
       
      // if(i==0){
      //         fwire_vertex=allhits[i]->Wire()->RawDigit()->Channel();
      //         ftime_vertex=allhits[i]->PeakTime();
      //     mf::LogInfo("KingaCluster")<<"for PLANE "<<plane
      // 			       <<" fwire_vertex= "<<fwire_vertex
      // 			       <<" ftime_vertex= "<<ftime_vertex;
      //           }   
       
      //mf::LogInfo("KingaCluster")<<"................................";
      //     mf::LogInfo("KingaCluster")<<" For hit# "<<i
      // 			       <<" w= "<<w
      // 			       <<" fwire_vertex[plane]= "
      // 			       <<fwire_vertex[plane]
      // 			       <<" ftime_vertex[plane]= "<<ftime_vertex[plane];
       
      int diff_w = allhits[i]->WireID().Wire - fwire_vertex_reco[plane];
      b_polar = diff_w*0.4; /**in cm*/
       
      //mf::LogInfo("KingaCluster")<<" diff_w= "<<diff_w;
      //mf::LogInfo("KingaCluster")<<" b_polar= "<<b_polar;
      a_polar = (allhits[i]->PeakTime() - ftime_vertex_reco[plane])* ftimetick *fdriftvelocity; /** in cm*/
       
       
      theta_polar = std::abs(std::asin(a_polar/std::sqrt(pow(a_polar,2)+pow(b_polar,2)))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/
       
       
      // We have 4 cases depending on which quater a hit is (origin being defined on a vertex):
       
      if     (b_polar == 0 && a_polar == 0) theta_polar =  90;/** in deg*/
      else if(b_polar == 0 && a_polar <  0) theta_polar = 180;/** in deg*/
      else if(b_polar == 0 && a_polar >  0) theta_polar =   0;/** in deg*/
      else if(b_polar >  0 && a_polar == 0) theta_polar =  90;/** in deg*/
      else if(b_polar <  0 && a_polar == 0) theta_polar = -90;/** in deg*/
      else if(b_polar >  0 && a_polar >  0) theta_polar =  90 - theta_polar;/** in deg*/
      else if(b_polar >  0 && a_polar <  0) theta_polar =  90 + theta_polar;/** in deg*/
      else if(b_polar <  0 && a_polar >  0) theta_polar = -(90 - theta_polar);/** in deg*/
      else if(b_polar <  0 && a_polar <  0) theta_polar = -(90 + theta_polar);/** in deg*/
        
      //fh_theta[plane]->Fill(theta_polar,allhits[i]->Charge());
      //mf::LogInfo("KingaCluster")<<"**** theta_polar= "<<theta_polar;
      // mf::LogInfo("KingaCluster")<<" a_polar= "<<a_polar<<" b_polar= "<<b_polar;
       
      //mf::LogInfo("KingaCluster")<<"w= "<<w<<" t= "<<std::setprecision(10)<<allhits[i]->PeakTime();
      //   mf::LogInfo("KingaCluster")<<" theta_polar= "<<theta_polar
      //		       <<" hit area= "
      //		       <<(allhits[i]->EndTime()-allhits[i]->StartTime())*ftimetick*fdriftvelocity*0.4;
  
      if (sigType == geo::kInduction ) {
         
        //       mf::LogInfo("KingaCluster")<<"plane= "<<plane
        // 				 <<" theta= "<<theta_polar
        // 				 <<"  channel= "<<allhits[i]->Wire()->RawDigit()->Channel()
        // 				 <<" time= "<<allhits[i]->PeakTime()
        // 				 <<" fwire_vertex[plane]= "<<fwire_vertex[plane]
        // 				 <<" ftime_vertex[plane]= "<<ftime_vertex[plane];
         
         
        fh_theta_ind->Fill(theta_polar);
        fh_theta_ind_2D->Fill(theta_polar,allhits[i]->Charge());
        fh_theta_ind_Area->Fill(theta_polar,(allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);
         
        Hit_Area_Ind->Fill((allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);
      }
      if (sigType == geo::kCollection ) {
         
        //       mf::LogInfo("KingaCluster")<<"plane= "
        // 				 <<plane<<"  w= "
        // 				 <<w<<" time= "
        // 				 <<allhits[i]->PeakTime()
        // 				 <<" theta= "<<theta_polar
        // 				 <<" fwire_vertex[plane]= "<<fwire_vertex[plane]
        // 				 <<" ftime_vertex[plane]= "<<ftime_vertex[plane];
         
         
        fh_theta_coll->Fill(theta_polar);
        fh_theta_coll_2D->Fill(theta_polar,allhits[i]->Charge());
        fh_theta_coll_Area->Fill(theta_polar,(allhits[i]->EndTime()-allhits[i]->StartTime())*ftimetick*fdriftvelocity*0.4);
        Hit_Area_Coll->Fill((allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);
         
      }// end if kCollection
    } // end loop over all hits
     
    return;
  }
  
      
      
  //..............................................................   
  void KingaCluster::FindMax(unsigned int cstat, 
			     unsigned int tpc, 
			     unsigned int plane){  
  
    // mf::LogInfo("KingaCluster")<<"No of bins= "<<fh_theta_ind->GetNbinsX();
    //   mf::LogInfo("KingaCluster")<<" Bincontent= "<<fh_theta_ind->GetBinContent(48);
    //   mf::LogInfo("KingaCluster")<<" Bincontent= "<<fh_theta_ind->GetBinContent(49);
    //   mf::LogInfo("KingaCluster")<<" Bincontent= "<<fh_theta_ind->GetBinContent(50);
    // std::vector<int> PossibleFinalMax, FinalMax;
  
    std::vector<int> startTimes;  //stores time of 1st local minimum
    // std::vector<int> maxBin;    //stores time of local maximum
    std::vector<int> endTimes;    //stores time of 2nd local minimum
    bool maxFound = false; //Flag for whether a value>threshold has been found
    
    int minTimeHolder;
   
    startTimes.clear();
    maxBin.clear();
    endTimes.clear();
  
    //mf::LogInfo("KingaCluster")<<"We have "<<allhits.size()<<" hits for plane "<<plane;
  
    //  double threshold=76*allhits.size()/2;
    //   double MinThreshold=1000;
    //  double threshold=6;
    //   double MinThreshold=4;
  
    // worked for pizero and nu_e :
    // double threshold=80000;
    //   double MinThreshold=60000;
    //............................
      
    // double threshold=30000;
    //   double MinThreshold=20000;
  
    //for lines:
    // double threshold=20000;
    //   double MinThreshold=10000;
  
    // double threshold=600;
    //   double MinThreshold=100;
  
    double threshold=0.2;
    //double MinHitThreshold=1; // Making sure that the peak found in fh_theta_coll_Area 
    // and _ind corresponds to more than 1 hit. 	       
    // Sometimes you can have a very large hit area that 
    // will produce a peak in that distribution but it   
    // only corresponds to 1 hit.                        
    double MinThreshold=0.1;
    
    
    int time=1;
    int ValidPeak=0;
    mf::LogInfo("KingaCluster")<< "Threshold that a peak must have in order to be considered a peak = "
			       << threshold
			       << ". For inflection points we must have the value to drop to "
			       <<MinThreshold;
  
    // get the signal type for the passed in plane
    geo::SigType_t sigType = fGeom->Plane(plane, tpc, cstat).SignalType();
    
    //......................................................... 
    //collection plane:
    if (sigType == geo::kCollection){
   
      for(int bin = 1; bin < fh_theta_coll_Area->GetNbinsX()+1; ++bin){
   
        if(fh_theta_coll_Area->GetBinContent(bin)   > fh_theta_coll_Area->GetBinContent(bin+1) && 
	   fh_theta_coll_Area->GetBinContent(bin+1) < fh_theta_coll_Area->GetBinContent(bin+2)) {
	  //only add points if we've already found a local max above threshold.
	  if(maxFound) {
	    endTimes.push_back(time+1);
	    maxFound = false;
	    //keep these in case new hit starts right away
	    minTimeHolder = time+2;
	  }
	  else minTimeHolder = time+1; 
        }
  
        //if not a minimum,-> test if we are at a local maximum
        //if so and the max value is above threshold add it and proceed.
        else if(fh_theta_coll_Area->GetBinContent(bin)   < fh_theta_coll_Area->GetBinContent(bin+1) &&
		fh_theta_coll_Area->GetBinContent(bin+1) > fh_theta_coll_Area->GetBinContent(bin+2) &&
		fh_theta_coll_Area->GetBinContent(bin+1) > threshold ) {
	  maxFound = true;
	  maxBin.push_back(time+1);
	  startTimes.push_back(minTimeHolder);         
        }
  
        time++; 
       
      }
    
      if(maxBin.size()==0){
    
        mf::LogInfo("KingaCluster")<<"COLLECTION PLANE: ";
        mf::LogInfo("KingaCluster")<<" COULDN'T FIND ANY MAXIMA IN YOUR THETA DISTRIBUTION!!!!"
				   << "  PROGRAM WILL NOT PROCEED!!!";
        fpeaks_found=0;
      }
      if(maxBin.size()>0){     
       
        // Lets make sure that the first bin in the maxBin corresponds to the highest peak:
  
        //std::vector<double> maxBinValues;
        for(unsigned int i=0;i<maxBin.size();i++){
	  maxBinValues.push_back(fh_theta_coll_Area->GetBinContent(maxBin[i]));
	  OriginalmaxBinValues.push_back(fh_theta_coll_Area->GetBinContent(maxBin[i]));
        }
        
        //       mf::LogInfo("KingaCluster")<<"The largest is at position:  "
        // 				 <<std::distance(maxBinValues.begin(),
        // 						 std::max_element(maxBinValues.begin(),
        // 								  maxBinValues.end()))
        // 				 <<" which corresponds to bin #   "
        // 				 <<maxBin[std::distance(maxBinValues.begin(),
        // 							std::max_element(maxBinValues.begin(),
        // 									 maxBinValues.end()))]
        // 				 <<" and its value= "
        // 				 <<std::max_element(maxBinValues.begin(),maxBinValues.end());
    
        // sort values from the largest to the smallest in maxBinValues, 
        // then find the corresponding bin numbers to create SortedMaxBin:
    
        sort(maxBinValues.begin(), maxBinValues.end());
  
        reverse(maxBinValues.begin(), maxBinValues.end());
  
        for(size_t i = 0; i < maxBinValues.size(); ++i){
  
	  std::vector<double>::iterator pos = std::find(OriginalmaxBinValues.begin(), 
							OriginalmaxBinValues.end(),
							maxBinValues[i]);
	  SortedMaxBin.push_back(maxBin[pos-OriginalmaxBinValues.begin()]);
  
        }
  
  
        //create SortedMaxBin vector whose first element is the global max:
    
        //SortedMaxBin.push_back(maxBin[std::distance(maxBinValues.begin(),
        // 						  std::max_element(maxBinValues.begin(),maxBinValues.end()))]);
        //   for(int i=0; i<maxBin.size();i++){
        //    if(i!=std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end())))
        //   SortedMaxBin.push_back(maxBin[i]);
    
        //   }
    
        //  mf::LogInfo("KingaCluster")<<"SortexMaxBin elements are: "<<std::endl;
        // for(unsigned int i=0; i<SortedMaxBin.size(); i++)
        // {
        // 
        // mf::LogInfo("KingaCluster")<<SortedMaxBin[i]<<std::endl;
        // 
        // }
        // int ValidPeak=0;
        // loop over maxima and find where they start on the left and right side, form cluster for each:
    
        for(size_t maxNo = 0; maxNo < SortedMaxBin.size(); ++maxNo){
    
	  //loop over the ranges and make sure that your peaks don't fall into already formed clusters
	  //mf::LogInfo("KingaCluster")<<"Right now we have "<<MaxStartPoint.size()<<" ranges";
	  if(MaxStartPoint.size() == 0) ValidPeak = 1;
	  for(size_t NoRange = 0; NoRange < MaxStartPoint.size(); ++NoRange){
	    //mf::LogInfo("KingaCluster")<<"Checking peak "<<SortedMaxBin[maxNo];
	    if(SortedMaxBin[maxNo]>MaxStartPoint[NoRange] && SortedMaxBin[maxNo]<MaxEndPoint[NoRange]){
	      //maxNo++;
	      //mf::LogInfo("KingaCluster")<<"this peak is out of the picture! --> "<<SortedMaxBin[maxNo]<<std::endl;
	      ValidPeak=0;
	      break;
	    }
	    else{
	      ValidPeak=1;
	      //mf::LogInfo("KingaCluster")<<"passed"<<std::endl;
	    }
	  }
	  if(ValidPeak==1){
     
	    // mf::LogInfo("KingaCluster")<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
	    FinalPeaks.push_back(SortedMaxBin[maxNo]);
	    mf::LogInfo("KingaCluster")<<"We are working on peak at bin #"<<SortedMaxBin[maxNo];
  
	    //start at the peak and go left
	    for(int LeftBin = SortedMaxBin[maxNo]-1; LeftBin > SortedMaxBin[maxNo]-30; --LeftBin){
	      if(fh_theta_coll_Area->GetBinContent(LeftBin)   < MinThreshold && 
		 fh_theta_coll_Area->GetBinContent(LeftBin-1) > fh_theta_coll_Area->GetBinContent(LeftBin)) {
		MaxStartPoint.push_back(LeftBin);
		// mf::LogInfo("KingaCluster")<<"picked option 1, startin point @bin "<<LeftBin;
		mf::LogInfo("KingaCluster")<<"For peak at bin= "<<SortedMaxBin[maxNo]
					   <<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"
					   <<" LeftBin="<<LeftBin<<"("<<-180+2*LeftBin<<" degrees)"
					   <<" RightBin= ";
		break;
	      }
	      else if(fh_theta_coll_Area->GetBinContent(LeftBin) < MinThreshold && 
		      fh_theta_coll_Area->GetBinContent(LeftBin-1) == 0){
		MaxStartPoint.push_back(LeftBin-1);
		//mf::LogInfo("KingaCluster")<<"picked option 2, startin point @bin "<<LeftBin-1;
		mf::LogInfo("KingaCluster")<<"For peak at bin= "<<SortedMaxBin[maxNo]
					   <<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"
					   <<" LeftBin="<<LeftBin-1<<"("<<-180+2*(LeftBin-1)<<" degrees)"
					   <<" RightBin= ";
		break;
             
	      }
	      else if (LeftBin==SortedMaxBin[maxNo]-29){
		mf::LogInfo("KingaCluster")<<" cannot find starting point of the peak!!!!!";
		fpeak_problem=1;
		return;
	      }
           
	    }// end loop over sorted max
       
	    for(int RightBin = SortedMaxBin[maxNo]+1; RightBin < SortedMaxBin[maxNo]+30; ++RightBin){
  
	      if(fh_theta_coll_Area->GetBinContent(RightBin)   < MinThreshold && 
		 fh_theta_coll_Area->GetBinContent(RightBin+1) > fh_theta_coll_Area->GetBinContent(RightBin)) {
		MaxEndPoint.push_back(RightBin);
		mf::LogInfo("KingaCluster")<<RightBin<<"("<<-180+2*RightBin<<" degrees)";
		break;
	      }
	      else if(fh_theta_coll_Area->GetBinContent(RightBin) < MinThreshold && 
		      fh_theta_coll_Area->GetBinContent(RightBin+1) == 0){
		MaxEndPoint.push_back(RightBin+1);
		mf::LogInfo("KingaCluster")<<RightBin+1<<"("<<-180+2*(RightBin+1)<<" degrees)";
		break;
	      }
	      else if(RightBin == SortedMaxBin[maxNo]+29){
		mf::LogInfo("KingaCluster")<<" cannot find end point of the peak!!!!!";
		fpeak_problem=1;
		return;
	      }
           
	    } // end loop over sorted max
        
	  } //valid peak
    
	  ValidPeak=0;
    
        }//peaks
  
      } //if maxBin.size()>0
    
      startTimes.clear();
      maxBin.clear();
      endTimes.clear();
    
    } //plane 1
   
    time=1;
    ValidPeak=0;
    
    //......................................................... 
   
    //induction plane
    if(sigType == geo::kInduction){
   
   
      //mf::LogInfo("KingaCluster")<<"No of bins= "<<fh_theta_ind_Area->GetNbinsX()<<std::endl;
      for(int bin = 1; bin < fh_theta_ind_Area->GetNbinsX()+1; ++bin){
   
        if(fh_theta_ind_Area->GetBinContent(bin)   > fh_theta_ind_Area->GetBinContent(bin+1) && 
	   fh_theta_ind_Area->GetBinContent(bin+1) < fh_theta_ind_Area->GetBinContent(bin+2)) {
	  //only add points if we've already found a local max above threshold.
	  if(maxFound) {
	    endTimes.push_back(time+1);
	    maxFound = false;
	    //keep these in case new hit starts right away
	    minTimeHolder = time+2;
	  }
	  else {
	    minTimeHolder = time+1; 
	  }
        }
        //if not a minimum test if we are at a local maximum
        //if so and the max value is above threshold add it and proceed.
        else if(fh_theta_ind_Area->GetBinContent(bin)   < fh_theta_ind_Area->GetBinContent(bin+1) &&
		fh_theta_ind_Area->GetBinContent(bin+1) > fh_theta_ind_Area->GetBinContent(bin+2) &&
		fh_theta_ind_Area->GetBinContent(bin+1) > threshold) {
	  maxFound = true;
	  maxBin.push_back(time+1);
	  startTimes.push_back(minTimeHolder); 
        
	  // && fh_theta_ind->GetBinContent(bin+1) > MinHitThreshold
        }
        time++;
   
      }
      mf::LogInfo("KingaCluster")<<"INDUCTION PLANE: ";
  
      if(maxBin.size()==0){
        mf::LogInfo("KingaCluster")<< " COULDN'T FIND ANY MAXIMA IN YOUR THETA DISTRIBUTION!!!! "
				   << " PROGRAM WILL NOT PROCEED!!!"<<std::endl;
        fpeaks_found=0;
      }
      
      if(maxBin.size()>0){   
        mf::LogInfo("KingaCluster")<<"maxBin.size()= "<<maxBin.size();
  
        for(size_t i = 0; i < maxBin.size(); ++i){
  
	  mf::LogInfo("KingaCluster")<<"maxTime is at bin = "
				     <<maxBin[i]<<" ("<<-180+2*maxBin[i]<<" degrees)"
				     <<" and its value is "
				     <<fh_theta_ind_Area->GetBinContent(maxBin[i])
				     <<"....................................";
        }//maxBin
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
  
        // Lets make sure that the first bin in the maxBin corresponds to the highest peak:
  
        //std::vector<double> maxBinValues;
        for(size_t i = 0; i < maxBin.size(); ++i){
	  maxBinValues.push_back(fh_theta_ind_Area->GetBinContent(maxBin[i]));
	  OriginalmaxBinValues.push_back(fh_theta_ind_Area->GetBinContent(maxBin[i]));
        }
  
        //       mf::LogInfo("KingaCluster")<<"The largest is at position:  "
        // 				 <<std::distance(maxBinValues.begin(),
        // 						 std::max_element(maxBinValues.begin(),maxBinValues.end()))
        // 				 <<" which corresponds to bin #   "
        // 				 <<maxBin[std::distance(maxBinValues.begin(),
        // 							std::max_element(maxBinValues.begin(),maxBinValues.end()))]
        // 				 <<" and its value= "
        // 				 <<std::max_element(maxBinValues.begin(),maxBinValues.end());
    
        // sort values from the largest to the smallest in maxBinValues, 
        // then find the corresponding bin numbers to create SortedMaxBin:
    
        sort(maxBinValues.begin(), maxBinValues.end());
        reverse (maxBinValues.begin(), maxBinValues.end());
  
        for(size_t i = 0; i < maxBinValues.size(); ++i){
	  std::vector<double>::iterator pos = std::find( OriginalmaxBinValues.begin(), 
							 OriginalmaxBinValues.end(),
							 maxBinValues[i]);
	  SortedMaxBin.push_back(maxBin[pos-OriginalmaxBinValues.begin()]);
        }
    
        //create SortedMaxBin vector whose first element is the global max:
    
    
        // std::vector<int> SortedMaxBin;
        //SortedMaxBin.push_back(maxBin[std::distance(maxBinValues.begin(),
        //					   std::max_element(maxBinValues.begin(),maxBinValues.end()))]);
        //   for(int i=0; i<maxBin.size();i++){
        //    if(i!=std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end())))
        //   SortedMaxBin.push_back(maxBin[i]);
        //   
        //   }
    
        //  mf::LogInfo("KingaCluster")<<"SortexMaxBin elements are: "<<std::endl;
        // for(unsigned int i=0; i<SortedMaxBin.size(); i++)
        // {
        // 
        // mf::LogInfo("KingaCluster")<<SortedMaxBin[i]<<std::endl;
        // 
        // }
  
        // loop over maxima and find where they start on the left and right side, form cluster for each:
    
        for(size_t maxNo = 0; maxNo < SortedMaxBin.size(); ++maxNo){
    
	  //loop over the ranges and make sure that your peaks don't fall into already formed clusters
	  if(MaxStartPoint.size()==0){ValidPeak=1;}
	  for(size_t NoRange = 0; NoRange < MaxStartPoint.size(); ++NoRange){
	    if(SortedMaxBin[maxNo]>MaxStartPoint[NoRange] && SortedMaxBin[maxNo]<MaxEndPoint[NoRange]){
	      ValidPeak=0;
	      break;
	    }
	    else {ValidPeak=1;}
	  }
    
	  if(ValidPeak==1){
     
	    //mf::LogInfo("KingaCluster")<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
	    FinalPeaks.push_back(SortedMaxBin[maxNo]);
	    //start at the peak and go left
	    for(int LeftBin = SortedMaxBin[maxNo]-1; LeftBin > SortedMaxBin[maxNo]-30; --LeftBin){
  
	      if(fh_theta_ind_Area->GetBinContent(LeftBin)   < MinThreshold && 
		 fh_theta_ind_Area->GetBinContent(LeftBin-1) > fh_theta_ind_Area->GetBinContent(LeftBin)) {
		MaxStartPoint.push_back(LeftBin);
              
		mf::LogInfo("KingaCluster")<<"For peak at bin= "<<SortedMaxBin[maxNo]
					   <<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"
					   <<" LeftBin="<<LeftBin<<"("<<-180+2*LeftBin<<" degrees)"
					   <<" RightBin= ";
          
		break;
	      }
	      else if(fh_theta_ind_Area->GetBinContent(LeftBin) < MinThreshold && 
		      fh_theta_ind_Area->GetBinContent(LeftBin-1) == 0){
		MaxStartPoint.push_back(LeftBin-1);
		//mf::LogInfo("KingaCluster")<<"picked option 2, startin point @bin "<<LeftBin-1;
		mf::LogInfo("KingaCluster")<<"For peak at bin= "<<SortedMaxBin[maxNo]
					   <<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"
					   <<" LeftBin="<<LeftBin-1<<"("<<-180+2*(LeftBin-1)<<" degrees)"
					   <<" RightBin= ";
		break;
             
	      }
	      else if (LeftBin==SortedMaxBin[maxNo]-30){
		mf::LogInfo("KingaCluster")<<"cannot find starting point of the peak!!!!!";
		fpeak_problem=1;
		return;
	      }         
	    }
    
    
	    for(int RightBin = SortedMaxBin[maxNo]+1; RightBin < SortedMaxBin[maxNo]+30; ++RightBin){
	      if(fh_theta_ind_Area->GetBinContent(RightBin)   < MinThreshold && 
		 fh_theta_ind_Area->GetBinContent(RightBin+1) > fh_theta_ind_Area->GetBinContent(RightBin)) {
		MaxEndPoint.push_back(RightBin);
		mf::LogInfo("KingaCluster")<<RightBin<<"("<<-180+2*RightBin<<" degrees)";
		break;
	      }
	      else if(fh_theta_ind_Area->GetBinContent(RightBin) < MinThreshold && 
		      fh_theta_ind_Area->GetBinContent(RightBin+1) == 0){
		MaxEndPoint.push_back(RightBin+1);
		mf::LogInfo("KingaCluster")<<RightBin+1<<"("<<-180+2*(RightBin+1)<<" degrees)";
		break;
	      }
	      else if(RightBin == SortedMaxBin[maxNo]+20){
		mf::LogInfo("KingaCluster")<<"cannot find end point of the peak!!!!!";
		fpeak_problem=1;
		return;
	      }
           
	    }
    
	  }
	  ValidPeak=0;
    
        }//peaks
      }// if maxBin.size()>0 this means that we can find peaks in the theta distribution 
    } //induction plane
  
    //......................................................... 
    //Now clear the histograms:
   
    // for(int bin=0; bin< fh_theta_ind_2D->GetNbinsX(); bin++){
    //    
    //    fh_theta_ind_2D->SetBinContent(bin,0);
    //    fh_theta_coll_2D->SetBinContent(bin,0);
    //    fh_theta_ind->SetBinContent(bin,0);
    //    fh_theta_coll->SetBinContent(bin,0);
    //    fh_theta_coll_Area->SetBinContent(bin,0);
    //    fh_theta_ind_Area->SetBinContent(bin,0);
    //  }
   
  }   
  
  //..............................................................   
  // void KingaCluster::FinalPeaks(){  
  // mf::LogInfo("KingaCluster")<<"In FinalPeaks()"<<std::endl;
  // // for(int i=0;i<maxBin.size();i++){
  // //   mf::LogInfo("KingaCluster")<<"maxTime is at bin = "<<maxBin[i]<<" and its value is "<<fh_theta_ind_Area->GetBinContent(maxBin[i])<<std::endl;
  // //   }
  // 
  // 
  // 
  // }
  
  //..............................................................   
  void KingaCluster::FindClusters(unsigned int cstat, 
				  unsigned int tpc, 
				  unsigned int plane){ 
  
    // get the signal type for the passed in plane
    geo::SigType_t sigType = fGeom->Plane(plane, tpc, cstat).SignalType();
  
    // First order check: make sure ranges of the peaks don't overlap, 
    // if they do you need to correct for that. The highest peaks should 
    // be left alone and we should work with the ones of the smallest peak 
    // value first. This means start from the end of the sorted peaks.
  
    std::vector<int> peak_with_wrong_range;
    std::vector<int> peak_with_which_it_ovelaps;
  
    for( int pk = FinalPeaks.size()-1; pk >= 0; --pk){
      mf::LogInfo("KingaCluster")<<"pk= "<<pk;
  
      for(size_t pk2 = 0; pk2 < FinalPeaks.size(); ++pk2){
        if(pk != (int) pk2 && 
	   ((MaxStartPoint[pk]<MaxEndPoint[pk2] && MaxStartPoint[pk]>MaxStartPoint[pk2]) ||
	    ( MaxEndPoint[pk]>MaxStartPoint[pk2] && MaxEndPoint[pk]<MaxEndPoint[pk2] ))
	   ){
	  mf::LogInfo("KingaCluster")<<"WRONG RANGE, NEED TO FIX IT FOR PEAK AT BIN #"<<FinalPeaks[pk];
      
	  peak_with_wrong_range.push_back(pk); //this gives peak#, NOT a bin#
	  peak_with_which_it_ovelaps.push_back(pk2); //this gives peak#, NOT a bin#
        }
      }
    }// end loop over pk
  
    int diff=0;
  
    for(size_t i = 0; i < peak_with_wrong_range.size(); ++i){
  
      //if front of a range overlaps with the back of the range already in place
      if(MaxStartPoint[peak_with_wrong_range[i]] < MaxEndPoint[peak_with_which_it_ovelaps[i]] && 
         MaxStartPoint[peak_with_wrong_range[i]] > MaxStartPoint[peak_with_which_it_ovelaps[i]]){
   
        diff=MaxEndPoint[peak_with_which_it_ovelaps[i]]-MaxStartPoint[peak_with_wrong_range[i]];
        MaxStartPoint[peak_with_wrong_range[i]]=MaxStartPoint[peak_with_wrong_range[i]]+ diff;
        mf::LogInfo("KingaCluster")<<"changing startpoint to bin #"<<MaxStartPoint[peak_with_wrong_range[i]];
  
      }
  
      //if back of a range overlaps with the front of the range already in place
      if(MaxEndPoint[peak_with_wrong_range[i]] > MaxStartPoint[peak_with_which_it_ovelaps[i]] &&
         MaxEndPoint[peak_with_wrong_range[i]] < MaxEndPoint[peak_with_which_it_ovelaps[i]]){
   
        diff=MaxEndPoint[peak_with_wrong_range[i]]-MaxStartPoint[peak_with_which_it_ovelaps[i]];
        MaxEndPoint[peak_with_wrong_range[i]]=MaxEndPoint[peak_with_wrong_range[i]]- diff;
    
        mf::LogInfo("KingaCluster")<<"changing endpoint to bin #"<<MaxEndPoint[peak_with_wrong_range[i]];
      }
      diff=0;
    }
  
    // First let's make sure that each range of peaks contains 
    // more than some specified number of hits. You can do it 
    // by knowing the ranges, going into histograms and counting 
    // the entries from start of the range to the end. If some range 
    // contains less than the specified number you need to remove the peak.
  
    mf::LogInfo("KingaCluster")<<" NO OF FINALPEAKS BEFORE EVALUATION IS: "<<FinalPeaks.size();
  
    int MinHitsInRange=2; //later make it a parameter. DEFINED BELOW, LOOK->>>
    double no_hits_in_range=0;
    std::vector<int> TempFinalPeaks;
    std::vector<int> TempMaxStartPoint;
    std::vector<int> TempMaxEndPoint;
  
    std::vector<int> positive_diff_end_minus_start;
    std::vector<int> positive_diff_start_minus_end;
    std::vector<int> diff_end_minus_start;
    std::vector<int> diff_start_minus_end;
  
    double closest_range_right_side=0;
    double closest_range_left_side=0;
    int this_is_the_first_range=0;
    int this_is_the_last_range=0;
    int well_separated=0;
  
    //---------------START BASIC EVALUATION FIRST---------------
  
    // Now I will eliminate all the peaks that have only 1 
    // hit in it's RANGE(not in actual peak):
    // However, if there is a group of peaks that have just 
    // one hit in their range and are all close to each other 
    // then we need to form a cluster out of them, call it "hand_made_peak". 
    // But this peak has to be well seprated from the rest!!!! so need 
    // to determine it's range also in order to check the separation
  
    std::vector<int> one_hit_peaks;
    std::vector<int> one_hit_peaks_start_point;
    std::vector<int> one_hit_peaks_end_point;
    std::vector<int> double_hit_peaks;
    std::vector<int> double_hit_peaks_start_point;
    std::vector<int> double_hit_peaks_end_point;
  
    for(size_t peak = 0; peak < MaxStartPoint.size(); ++peak){
      for(int bin=MaxStartPoint[peak]; bin<MaxEndPoint[peak];bin++){
  
        //mf::LogInfo("KingaCluster")<<"bin= "<<bin<<std::endl;
        if(sigType == geo::kInduction){
	  no_hits_in_range+=fh_theta_ind->GetBinContent(bin);
	  //mf::LogInfo("KingaCluster")<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range;
        }
        if(sigType == geo::kCollection){
	  no_hits_in_range+=fh_theta_coll->GetBinContent(bin);
	  //mf::LogInfo("KingaCluster")<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range;
        }
   
      }//loop thru bins for each peak
    
      if(no_hits_in_range == 1){
        one_hit_peaks.push_back(FinalPeaks[peak]);
        one_hit_peaks_start_point.push_back(MaxStartPoint[peak]);
        one_hit_peaks_end_point.push_back(MaxEndPoint[peak]);
      }
    
      if(no_hits_in_range == 2){
        double_hit_peaks.push_back(FinalPeaks[peak]);
        double_hit_peaks_start_point.push_back(MaxStartPoint[peak]);
        double_hit_peaks_end_point.push_back(MaxEndPoint[peak]);
      }
    
    
      mf::LogInfo("KingaCluster")<<"no_hits_in_range= "<<no_hits_in_range
				 <<" for peak at bin# "<<FinalPeaks[peak];
      if(no_hits_in_range > 2){
        TempFinalPeaks.push_back(FinalPeaks[peak]);
        TempMaxStartPoint.push_back(MaxStartPoint[peak]);
        TempMaxEndPoint.push_back(MaxEndPoint[peak]);
      }
      no_hits_in_range = 0;
    } //loop thru all peaks
  
    FinalPeaks=TempFinalPeaks;
    MaxStartPoint=TempMaxStartPoint;
    MaxEndPoint=TempMaxEndPoint;
  
    TempFinalPeaks.clear();
    TempMaxStartPoint.clear();
    TempMaxEndPoint.clear();
    no_hits_in_range=0;
  
    mf::LogInfo("KingaCluster")<<" So now the size of FinalPeaks.size()="<<FinalPeaks.size();
  
    std::vector<int> grouped_double_hit_peaks;
  
    // First, let's look at 2-hit peaks and a case in which 
    // there are TWO 2-hit peaks right next to each other. 
    // If this is the case make one peak out of them by 
    // increstd::asing the range of the peak, make the peak be 
    // in the middle of the range.
  
    //.......2-hit peaks START................................
  
    // std::vector<int> used_double_hit_peaks;
    // used_double_hit_peaks.clear();
    // 
    //   for(unsigned int i=0;i<double_hit_peaks.size(); i++){
    //     for(unsigned int j=0;j<double_hit_peaks.size(); j++){
    // 
    //   if(i!=j && double_hit_peaks_end_point[i]==double_hit_peaks_start_point[j]){
    //   
    //   //create one peak with range of both peaks:
    //   //Also make sure you do the combination only once (you can have i=1 and j=2 && //i=2 and j=1...you should just do one since they are the same)
    //   
    //   if(std::find(FinalPeaks.begin(),FinalPeaks.end(),(double_hit_peaks[i]+double_hit_peaks[j])*0.5)==FinalPeaks.end()){
    //   FinalPeaks.push_back((double_hit_peaks[i]+double_hit_peaks[j])*0.5);
    //   MaxStartPoint.push_back(double_hit_peaks_start_point[i]);
    //   MaxEndPoint.push_back(double_hit_peaks_end_point[j]);
    //   
    //   used_double_hit_peaks.push_back(i);
    //   used_double_hit_peaks.push_back(j);
    //   mf::LogInfo("KingaCluster")<<"Creating one peak from TWO 2-hit peaks at peak at bin #"<<(double_hit_peaks[i]+double_hit_peaks[j])*0.5<<" and range in bins: ["<<double_hit_peaks_start_point[i]<<", "<<double_hit_peaks_end_point[j]<<"]"<<std::endl;
    //   
    //   }//if this peak doesn't exist yet, make it
    //   
    //   }
    // 
    // 
    // 
    // 
    //     }
    //    }
  
    //Now fill FinalPeaks with the rest of double_hit_peaks, unless you can form a bigger cluster out of very close-by 2-hit-clusters:
  
  
    //............. GROUP 2-HIT-CLUSTERS HERE, THE SAME WAY AS 1-HIT-CLUSTERS( see work on 1-hit clusters below)
  
    // std::vector<int> temp_double_hit_peaks;
    // temp_double_hit_peaks.clear();
    // std::vector<int> temp_double_hit_peaks_start_point;
    // temp_double_hit_peaks_start_point.clear();
    // std::vector<int> temp_double_hit_peaks_end_point;
    // temp_double_hit_peaks_end_point.clear();
  
  
    // for(unsigned int k=0;k<double_hit_peaks.size();k++){
    // 
    // if(std::find(used_double_hit_peaks.begin(),used_double_hit_peaks.end(),double_hit_peaks[k])==used_double_hit_peaks.end()){
    // 
    // temp_double_hit_peaks.push_back(double_hit_peaks[k]);
    // temp_double_hit_peaks_start_point.push_back(double_hit_peaks_start_point[k]);
    // temp_double_hit_peaks_end_point.push_back(double_hit_peaks_end_point[k]);
    // }
    // }
  
    // double_hit_peaks.clear();
    // double_hit_peaks_start_point.clear();
    // double_hit_peaks_end_point.clear();
    // double_hit_peaks=temp_double_hit_peaks;
    // double_hit_peaks_start_point=temp_double_hit_peaks_start_point;
    // double_hit_peaks_end_point=temp_double_hit_peaks_end_point;
  
    int diff_between_double_hit_peaks=7;
    std::vector<int> used_peak;
    used_peak.clear();
    if(double_hit_peaks.size() >= 2){
  
      for(size_t i = 0; i < double_hit_peaks.size(); ++i){
        for(size_t j = 0; j < double_hit_peaks.size(); ++j){
      
	  if(grouped_double_hit_peaks.size() == 0){
        
	    // loop thru existing ranges and make sure that your 
	    // 2-hit peaks are close to each other but also well 
	    // separated from the existing ranges
	    for(size_t r = 0; r < FinalPeaks.size(); ++r){ 
	      mf::LogInfo("KingaCluster")<<"### FinalPeaks.size()= "<<FinalPeaks.size();
	      mf::LogInfo("KingaCluster")<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]
					 <<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r];
  
	      if(i!=j && 
		 abs(double_hit_peaks[i]-double_hit_peaks[j])<=diff_between_double_hit_peaks ){ 
         
		// 	     mf::LogInfo("KingaCluster")<<"abs(one_hit_peaks[i]-MaxStartPoint[r])= "
		// 					<<abs(one_hit_peaks[i]-MaxStartPoint[r])<<" :: "
		// 					<<one_hit_peaks[i]<<" - "<<MaxStartPoint[r];
         
		// 	      mf::LogInfo("KingaCluster")<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]
		// 					 <<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r];
         
		// 	      mf::LogInfo("KingaCluster")<<"one_hit_peaks[i]= "<<one_hit_peaks[i]
		// 					 <<" one_hit_peaks[j]= "<<one_hit_peaks[j];
         
         
		//if you can't find this peak in the group, then add it:
		if(std::find(grouped_double_hit_peaks.begin(),
			     grouped_double_hit_peaks.end(),
			     double_hit_peaks[i])==grouped_double_hit_peaks.end()){
  
		  grouped_double_hit_peaks.push_back(double_hit_peaks[i]);
  
		  mf::LogInfo("KingaCluster")<<"ADDING 2-HIT PEAK #"<<double_hit_peaks[i]<<" TO THE GROUP";
        
		  used_peak.push_back(double_hit_peaks[i]);
		}
		if(std::find(grouped_double_hit_peaks.begin(),
			     grouped_double_hit_peaks.end(),
			     double_hit_peaks[j])==grouped_double_hit_peaks.end()){
		  grouped_double_hit_peaks.push_back(double_hit_peaks[j]);
		  mf::LogInfo("KingaCluster")<<"ADDING 2-HIT PEAK #"<<double_hit_peaks[j]<<" TO THE GROUP";
       
		  used_peak.push_back(double_hit_peaks[j]);
		}
	      } //if difference between peak bin#s is small 
	    }// loop thru already existing peaks 
	  } //if size=0 
        } //first loop thru each 1-hit peak
      } //second loop thru each 1-hit peak
    
      int failed=0;
      //Now go thru all the peaks again to find the third matching peak:
      if(grouped_double_hit_peaks.size()>0){
        for(size_t y = 0; y < double_hit_peaks.size(); ++y){
	  //now every added peak must be very close to EACH already added peak in the group
	  for(size_t q = 0; q < grouped_double_hit_peaks.size(); ++q){
         
	    if(abs(double_hit_peaks[y]-grouped_double_hit_peaks[q]) > 5){
	      failed=1;
	      break;
	    }
         
	  } //for
         
         
	  if(failed==0 && 
	     std::find(grouped_double_hit_peaks.begin(),
		       grouped_double_hit_peaks.end(),
		       double_hit_peaks[y]) == grouped_double_hit_peaks.end() ){
         
	    grouped_double_hit_peaks.push_back(double_hit_peaks[y]);
         
	    used_peak.push_back(double_hit_peaks[y]);
         
	  }
	  failed=0;
        } //loop thru all 2-hit peaks
      }// if size >0
    }
  
  
    //so now we should have a group of 2-hit peaks that are close to each other. 
    if(grouped_double_hit_peaks.size()>=2){
  
      //form a peak by taking the average of them and picking the one that's closest to that value:
      int sum=0;
      int made_peak=200;
   
      for(size_t k=0; k < grouped_double_hit_peaks.size(); ++k){
      
        sum+=grouped_double_hit_peaks[k];
        //       mf::LogInfo("KingaCluster")<<"grouped_double_hit_peaks.size()= "
        // 				 <<grouped_double_hit_peaks.size()
        // 				 <<" sum= "<<sum;
      }
     
      made_peak=sum/grouped_double_hit_peaks.size();
      mf::LogInfo("KingaCluster")<<"---------MADE A HOME_MADE_PEAK FROM GROUPS OF ***2***"
				 << "-HIT CLUSTERS--------THIS PEAK IS AT BIN #"<<made_peak;
  
      FinalPeaks.push_back(made_peak);
      //Now work on the range. For now just start at the smallest peak bin-1 and end at the largest peak bin+1:    
      std::sort(grouped_double_hit_peaks.begin(),grouped_double_hit_peaks.end());
      
      MaxStartPoint.push_back(grouped_double_hit_peaks[0]-1);
      MaxEndPoint.push_back(grouped_double_hit_peaks[grouped_double_hit_peaks.size()-1]+1);
      
      mf::LogInfo("KingaCluster")<<" its range is ["<<grouped_double_hit_peaks[0]-1
				 <<", "<<grouped_double_hit_peaks[grouped_double_hit_peaks.size()-1]+1<<"]";
      
    }
  
  
    //If we weren't able to group 2-hit-clusters together then just add them to FinalPeaks the way they are:
    for(size_t f = 0; f < double_hit_peaks.size(); ++f){
  
      if(std::find(used_peak.begin(),used_peak.end(),double_hit_peaks[f])==used_peak.end()){
        mf::LogInfo("KingaCluster")<<"adding peak at bin #"<<double_hit_peaks[f];
        FinalPeaks.push_back(double_hit_peaks[f]);
        MaxStartPoint.push_back(double_hit_peaks_start_point[f]);
        MaxEndPoint.push_back(double_hit_peaks_end_point[f]);
      }
    }
    //..............END OF GROUPING 2-HIT-CLUSTERS
  
    mf::LogInfo("KingaCluster")<<"---------------------*******-----------------------------------";
    mf::LogInfo("KingaCluster")<<"No of FinalPeaks after 2-hit evaluation = "<<FinalPeaks.size();
  
    for(size_t peak = 0; peak < FinalPeaks.size(); ++peak){
      mf::LogInfo("KingaCluster")<<"peak at bin # "<<FinalPeaks[peak]
				 <<" ("<<-180+2*FinalPeaks[peak]
				 <<" degrees). Its range is ["<<-180+2*MaxStartPoint[peak]
				 <<", "<<-180+2*MaxEndPoint[peak]<<" ]";
    }
  
  
    mf::LogInfo("KingaCluster")<<"-----------------------*******---------------------------------";
  
    //.............. 2-hit peaks END............................
  
    //Before I start actuall work on these 1-hit peaks I will 
    // get rid of all the ones that are too close to the already existing ranges:
  
    std::vector<int> temp_one_hit_peaks;
    int bad_one_hit_peak=0;
  
    for(size_t i = 0; i < one_hit_peaks.size(); ++i){
      for(size_t g = 0; g < FinalPeaks.size(); ++g){
  
        if(abs(one_hit_peaks[i]-MaxStartPoint[g])<6 || abs(one_hit_peaks[i]-MaxEndPoint[g])<6){
     
	  bad_one_hit_peak=1;
	  break;
     
        }
      }
   
      if (bad_one_hit_peak==0) temp_one_hit_peaks.push_back(one_hit_peaks[i]);
      bad_one_hit_peak=0;
    }
  
    one_hit_peaks.clear();
    one_hit_peaks=temp_one_hit_peaks;
  
  
    mf::LogInfo("KingaCluster")<<" ### NOW THE UPDATED NO OF one_hit_peaks= "
			       <<one_hit_peaks.size()<<" they are at bins : ";
  
    for(size_t u = 0; u < one_hit_peaks.size(); ++u){
  
      mf::LogInfo("KingaCluster")<<one_hit_peaks[u];
  
    }
  
    //-----------WORK ON 1-HIT PEAKS----------------------------
  
    std::vector<int> grouped_one_hit_peaks;
    int hand_made_peak=0;
    int very_well_separated_two_hit_peak=0;
    int two_hits_only=0;
    int diff_between_one_hit_peaks=0;
    int event_is_clean=0;
  
  
    /// \todo where do the 4 and 7 come from, should be tunable parameters
    if(one_hit_peaks.size() <= 4){
      event_is_clean=1;
      mf::LogInfo("KingaCluster")<<"event_is_clean";
    }
  
    if(one_hit_peaks.size() >= 4) diff_between_one_hit_peaks=4;
    else diff_between_one_hit_peaks = 7;
  
    mf::LogInfo("KingaCluster")<<"one_hit_peaks size= "<<one_hit_peaks.size();
  
    if(one_hit_peaks.size() >= 2){
  
      for(size_t i = 0; i < one_hit_peaks.size(); ++i){
        for(size_t j = 0; j < one_hit_peaks.size(); ++j){
      
	  if(grouped_one_hit_peaks.size() == 0){
        
	    // loop thru existing ranges and make sure that your 1-hit 
	    // peaks are not close to each other but also well separated 
	    // from the existing ranges
	    for(size_t r = 0; r < FinalPeaks.size(); ++r){ 
	      mf::LogInfo("KingaCluster")<<"### FinalPeaks.size()= "<<FinalPeaks.size();
	      mf::LogInfo("KingaCluster")<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]
					 <<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r];
	      if(i!=j && abs(one_hit_peaks[i]-one_hit_peaks[j])<=diff_between_one_hit_peaks ){ 
         
		//if you can't find this peak in the group, then add it:
		if(std::find(grouped_one_hit_peaks.begin(),
			     grouped_one_hit_peaks.end(),
			     one_hit_peaks[i]) == grouped_one_hit_peaks.end()){
		  grouped_one_hit_peaks.push_back(one_hit_peaks[i]);
		  mf::LogInfo("KingaCluster")<<"ADDING PEAK #"<<one_hit_peaks[i]<<" TO THE GROUP";
		}
		if(std::find(grouped_one_hit_peaks.begin(),
			     grouped_one_hit_peaks.end(),
			     one_hit_peaks[j]) == grouped_one_hit_peaks.end()){
		  grouped_one_hit_peaks.push_back(one_hit_peaks[j]);
		  mf::LogInfo("KingaCluster")<<"ADDING PEAK #"<<one_hit_peaks[j]<<" TO THE GROUP";
		}
	      } //if difference between peak bin#s is small 
	    }// loop thru already existing peaks 
	  } //if size=0 
        } //first loop thru each 1-hit peak
      } //second loop thru each 1-hit peak
    
      int failed=0;
      //Now go thru all the peaks again to find the third matching peak:
      if(grouped_one_hit_peaks.size() > 0){
        for(size_t y = 0;y < one_hit_peaks.size(); ++y){
	  //now every added peak must be very close to EACH already added peak in the group
	  for(size_t q = 0; q < grouped_one_hit_peaks.size(); ++q){
         
	    if(abs(one_hit_peaks[y]-grouped_one_hit_peaks[q])>5){
	      failed=1;
	      break;
	    }
         
	  } //for
         
	  if(failed == 0 && 
	     std::find(grouped_one_hit_peaks.begin(),
		       grouped_one_hit_peaks.end(),
		       one_hit_peaks[y])==grouped_one_hit_peaks.end() ){
         
	    grouped_one_hit_peaks.push_back(one_hit_peaks[y]);
         
	  }
	  failed=0;
        } //loop thru all 1-hit peaks
      }// if size >0
    } //if we have more than 3 1-hit peaks
  
    //so now we should have a group of 1-hit peaks that are close to each other.
    /// \todo 3 and 2 are magic numbers, these should either be parameters or explained
    if(grouped_one_hit_peaks.size() >= 3 || 
       (grouped_one_hit_peaks.size() == 2 && event_is_clean==1)){
  
      //form a peak by taking the average of them and picking the one that's closest to that value:
      int sum=0;
      int made_peak=200;
   
      for(size_t k = 0; k < grouped_one_hit_peaks.size(); ++k){
      
        sum += grouped_one_hit_peaks[k];
      }
     
      made_peak=sum/grouped_one_hit_peaks.size();
      mf::LogInfo("KingaCluster")<<"---------MADE A HOME_MADE_PEAK FROM GROUPS OF 1-HIT CLUSTERS"
				 << "--------THIS PEAK IS AT BIN #"<<made_peak;
      FinalPeaks.push_back(made_peak);
      // Now work on the range. For now just start at the smallest 
      // peak bin-1 and end at the largest peak bin+1:
      
      std::sort(grouped_one_hit_peaks.begin(),grouped_one_hit_peaks.end());
      
      MaxStartPoint.push_back(grouped_one_hit_peaks[0]-1);
      MaxEndPoint.push_back(grouped_one_hit_peaks[grouped_one_hit_peaks.size()-1]+1);
      
      mf::LogInfo("KingaCluster")<<" its range is ["<<grouped_one_hit_peaks[0]-1
				 <<", "<<grouped_one_hit_peaks[grouped_one_hit_peaks.size()-1]+1<<"]";
      
      if(grouped_one_hit_peaks.size()>=3){hand_made_peak=1;}
      if(grouped_one_hit_peaks.size()==2){two_hits_only=1;}
      
    }
  
  
  
    //-----------END OF WORK ON 1-HIT PEAKS---------------------
  
  
    mf::LogInfo("KingaCluster")<<"After BASIC EVALUATION we now have "<<FinalPeaks.size()<<" peaks";
    mf::LogInfo("KingaCluster")<<"FinalPeaks are at bin(s):  ";
    for(size_t i = 0; i < FinalPeaks.size(); ++i){
      mf::LogInfo("KingaCluster")<<FinalPeaks[i]<<" which corresponds to angle=  "<<-180+2*FinalPeaks[i];
    }
  
    //---------------------END BASIC EVALUATION-----------------------------
  
    for(size_t peak = 0; peak < MaxStartPoint.size(); ++peak){
  
  
      // let's calculate how many bins away is the closest cluster 
      // to the one in question. You need to look to the left and 
      // to the right of each range to determine it. This is needed 
      // to pick the right MinHitsInRange value. Motivation: for very 
      // separated clusters in histos we want to be more lenient even 
      // though their peaks are not high. For very crowded environment 
      // want to be more strict.
   
      //loop thru all the other peaks to figure out the bin distance:
      for(size_t peak2 = 0; peak2 < MaxStartPoint.size(); ++peak2){
    
        if(peak != peak2){
	  diff_end_minus_start.push_back(MaxStartPoint[peak2]-MaxEndPoint[peak]);
	  diff_start_minus_end.push_back(MaxStartPoint[peak]-MaxEndPoint[peak2]);
        }
      }
    
      for(size_t diff=0; diff < diff_end_minus_start.size(); ++diff){
      
        if(diff_end_minus_start[diff] >= 0){
	  positive_diff_end_minus_start.push_back(diff_end_minus_start[diff]); 
        }
        if(diff_start_minus_end[diff] >= 0){
	  positive_diff_start_minus_end.push_back(diff_start_minus_end[diff]); 
        }
      
      }
  
      //now take the minimum and this is your closest range in bin numbers:
      if(positive_diff_end_minus_start.size() > 0){
        closest_range_right_side = *std::min_element(positive_diff_end_minus_start.begin(),
						     positive_diff_end_minus_start.end());
      }
      else if(positive_diff_end_minus_start.size() == 0) this_is_the_last_range = 1;
   
      if(positive_diff_start_minus_end.size() > 0){
        closest_range_left_side = *std::min_element(positive_diff_start_minus_end.begin(),
						    positive_diff_start_minus_end.end());
      }
      else if(positive_diff_start_minus_end.size() == 0) this_is_the_first_range=1;
   
   
      //if we only have 2 hits in the peak, it better be very well separated, otherwise not a valid peak!
      /// \todo the 16 appears to be a magic number and should be explained or a parameter
      if((two_hits_only == 1) && 
         (closest_range_right_side >= 16 || this_is_the_last_range  == 1 ) && 
         (closest_range_left_side  >= 16 || this_is_the_first_range == 1)){
        MinHitsInRange = 2;
        very_well_separated_two_hit_peak=1;
      }
    
      //if range is well separated (or first or last):
      /// \todo the 8 appears to be a magic number and should be explained or a parameter
      if((closest_range_right_side >= 8 || this_is_the_last_range  == 1 ) && 
         (closest_range_left_side  >= 8 || this_is_the_first_range == 1)){
        MinHitsInRange=2;
        well_separated=1;
        mf::LogInfo("KingaCluster")<<"Peak at bin #"<<FinalPeaks[peak]<<" is ^^^ WELL SEPARATED ^^^";
   
        if(this_is_the_last_range  == 1) mf::LogInfo("KingaCluster")<<" because it's the last range and ";
        if(this_is_the_first_range == 1) mf::LogInfo("KingaCluster")<<" because it's the first range and ";
        if(closest_range_right_side >= 8){
	  mf::LogInfo("KingaCluster")<<"The closest range from right side is "
				     <<closest_range_right_side<<" bins away";
        }
        if(closest_range_left_side >= 8){
	  mf::LogInfo("KingaCluster")<<"The closest range from left side is "
				     <<closest_range_left_side<<" bins away";
        }
   
      }
      else MinHitsInRange=3; 
   
   
      //........................
      for(int bin = MaxStartPoint[peak]; bin < MaxEndPoint[peak]; ++bin){
   
        //mf::LogInfo("KingaCluster")<<"bin= "<<bin<<std::endl;
        if(sigType == geo::kInduction){
	  no_hits_in_range+=fh_theta_ind->GetBinContent(bin);
	  //mf::LogInfo("KingaCluster")<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range;
        }
        if(sigType == geo::kCollection){
	  no_hits_in_range+=fh_theta_coll->GetBinContent(bin);
	  //mf::LogInfo("KingaCluster")<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range;
        }
   
      }
   
      mf::LogInfo("KingaCluster")<<"no_hits_in_range= "<<no_hits_in_range
				 <<" for peak at bin # "<<FinalPeaks[peak]
				 <<" ("<<-180+2*FinalPeaks[peak]<<" degrees). Its range is ["
				 <<-180+2*MaxStartPoint[peak]<<", "<<-180+2*MaxEndPoint[peak]<<" ]";
   
      if(sigType == geo::kInduction){
        if((fh_theta_ind_Area->GetBinContent(FinalPeaks[peak]) > 0.4 && 
	    no_hits_in_range>=MinHitsInRange) || 
	   hand_made_peak==1 || 
	   very_well_separated_two_hit_peak==1){
	  TempFinalPeaks.push_back(FinalPeaks[peak]);
	  TempMaxStartPoint.push_back(MaxStartPoint[peak]);
	  TempMaxEndPoint.push_back(MaxEndPoint[peak]);
        }
        else if(fh_theta_ind_Area->GetBinContent(FinalPeaks[peak]) <= 0.4 && 
		no_hits_in_range >= MinHitsInRange && 
		well_separated == 1){
	  TempFinalPeaks.push_back(FinalPeaks[peak]);
	  TempMaxStartPoint.push_back(MaxStartPoint[peak]);
	  TempMaxEndPoint.push_back(MaxEndPoint[peak]);
        }
   
      }//inductionplane 
   
      if(sigType == geo::kCollection){
        if((fh_theta_coll_Area->GetBinContent(FinalPeaks[peak]) > 0.4 && 
	    no_hits_in_range>=MinHitsInRange) || 
	   hand_made_peak==1){
	  TempFinalPeaks.push_back(FinalPeaks[peak]);
	  TempMaxStartPoint.push_back(MaxStartPoint[peak]);
	  TempMaxEndPoint.push_back(MaxEndPoint[peak]);
        }
        else if(fh_theta_coll_Area->GetBinContent(FinalPeaks[peak]) < 0.4 && 
		no_hits_in_range>=MinHitsInRange && 
		well_separated == 1){
	  TempFinalPeaks.push_back(FinalPeaks[peak]);
	  TempMaxStartPoint.push_back(MaxStartPoint[peak]);
	  TempMaxEndPoint.push_back(MaxEndPoint[peak]);
   
        }
   
      }//collection plane
   
   
      well_separated=0;
      no_hits_in_range=0;
      positive_diff_end_minus_start.clear();
      positive_diff_start_minus_end.clear();
      diff_end_minus_start.clear();
      diff_start_minus_end.clear();
      this_is_the_first_range=0;
      this_is_the_last_range=0;
   
    } //for each peak
  
  
    FinalPeaks=TempFinalPeaks;
    MaxStartPoint=TempMaxStartPoint;
    MaxEndPoint=TempMaxEndPoint;
  
    mf::LogInfo("KingaCluster")<<" NO OF FINALPEAKS ***AFTER*** EVALUATION IS: "<<FinalPeaks.size();
  
    // If you just get 1 peak (which most likely corresponds to the 
    // muon track in case of CCQE) then let's check if there is a 
    // possibility of forming another cluster. Allow this formation 
    // if the event is CLEAN (small number of 1-hit peaks, and overall 
    // small no of all found peaks in the histos). Look for a peak 
    // that is VERY FAR AWAY from the already formed peak(s). But make 
    // sure that we get at least 2-hit in this newely formed cluster. 
    // Allow to run this test for both 1-peak-hits only since we
    //already checked for 2-hit-peak possibilities above.
  
    if(FinalPeaks.size() == 1 && event_is_clean == 1){
      mf::LogInfo("KingaCluster")<<"$$$$  will try to come up with another peak since the event is clean";
      int separation=0;
      std::vector<int> separated_one_hit_peaks;
      std::vector<int> close_and_separated_one_h_pk;
  
      // now take all the 1 hit peaks and see if you can find two which 
      // are very far from the already formed peak:
      for(size_t onepk = 0; onepk < one_hit_peaks.size(); ++onepk){
        for(size_t fpk = 0; fpk < FinalPeaks.size(); ++fpk){
  	
	  //first pick the ones which are very far from the already existing peak
	  /// \todo 17 is a magic number and should be made a parameter
	  if(abs(one_hit_peaks[onepk]-FinalPeaks[fpk]) > 17){
   
	    separated_one_hit_peaks.push_back(one_hit_peaks[onepk]);
	    mf::LogInfo("KingaCluster")<<"adding "<<one_hit_peaks[onepk]<<" to separated_one_hit_peaks";
   
	  }
        }
      }
  
      //Now check if you can find 2 1-hit-peaks which are also not so far from each other:
      for(size_t i = 0; i < separated_one_hit_peaks.size(); ++i){
        for(size_t j = 0; j < separated_one_hit_peaks.size(); ++j){
   
	  /// \todo 11 is a magic number and should be a parameter
	  if(i != j && abs(separated_one_hit_peaks[i]-separated_one_hit_peaks[j]) < 11){
	    separation=abs(separated_one_hit_peaks[i]-separated_one_hit_peaks[j]);
	    close_and_separated_one_h_pk.push_back(separated_one_hit_peaks[i]);
	    close_and_separated_one_h_pk.push_back(separated_one_hit_peaks[j]);
	    break;
	  }
   
        }
        if(close_and_separated_one_h_pk.size() == 2) break;
      }
  
      //form a peak out of 2 1-hit-peaks:
      if(close_and_separated_one_h_pk.size() == 2){
  
        int sum=0; 
  
        for(size_t y = 0; y < close_and_separated_one_h_pk.size(); ++y){
  
	  sum+=close_and_separated_one_h_pk[y];
	  // and the range should be the start point of the first hit, 
	  // the end point should be the end point of the second hit, 
	  // so you could sort them. But, actually the range doesnt 
	  // really matter here, just make it a few bins to the right and left!
  
        }
  
        int made_peak = sum/close_and_separated_one_h_pk.size();
        mf::LogInfo("KingaCluster")<<"---------------------------------------------------------"
				   <<" WARNING:: WILL MAKE AN EXCEPTION FOR THIS EVENT "
				   << " (event is clean, only 1 peak originally found) "
				   << " AND CREATE A NEW PEAK AT BIN #"<<made_peak
				   <<"---------------------------------------------------------";
  
        FinalPeaks.push_back(made_peak);
        MaxStartPoint.push_back(made_peak-0.5*separation);
        MaxEndPoint.push_back(made_peak+0.5*separation);
      }//if we have 2 1-hit-clusters
    }
  
    mf::LogInfo("KingaCluster")<<" FinalPeaks.size()="<<FinalPeaks.size();
  
    // One last check to Make sure that the peaks I have now are not too 
    // close to each other, otherwise throw the smallest which is very close. 
    // Right now they are ordered based on their area height, so we should 
    // be dropping them from the back, ie the smallest. (This piece of code 
    // might be move more up if needed, perhaps at the top of the basic 
    // evaluation section)
  
    TempFinalPeaks.clear();
    TempMaxStartPoint.clear();
    TempMaxEndPoint.clear();
  
    // outside loop is an int instead of size_t to prevent --pk 
    // trying to take it less than 0
    /// \todo should really use iterators for these loops
    for(int pk = FinalPeaks.size()-1; pk >= 0; --pk){
      for(size_t pk2 = 0; pk2 < FinalPeaks.size(); ++pk2){
  
        /// \todo 3 is a magic number
        if(pk != (int)pk2 && abs(FinalPeaks[pk]-FinalPeaks[pk2]) < 3){
    
	  for(size_t i = 0; i < FinalPeaks.size(); ++i){
	    if(int(i) != pk){ 
	      TempFinalPeaks.push_back(FinalPeaks[i]);
	      TempMaxStartPoint.push_back(MaxStartPoint[i]);
	      TempMaxEndPoint.push_back(MaxEndPoint[i]);
	    }
	  }
    
	  FinalPeaks.clear();
	  MaxStartPoint.clear();
	  MaxEndPoint.clear();
	  FinalPeaks=TempFinalPeaks;
	  MaxStartPoint=TempMaxStartPoint;
	  MaxEndPoint=TempMaxEndPoint;
	  TempFinalPeaks.clear();
	  TempMaxStartPoint.clear();
	  TempMaxEndPoint.clear();
	  break;
   
        }
      }
    }
  
    mf::LogInfo("KingaCluster")<<" After making sure that each peak is more than "
			       << "3 bins away from the previous one, FinalPeaks.size()="
			       <<FinalPeaks.size();
  
    // Lastly, let's make sure that if there are ranges of peaks 
    // that lay right next to each other, their peak signal is 
    // actually different and not too small:
    // I already set that the minimum peak must be greater than 
    // 0.4 so here the check is for small peaks of that order up to ~0.6
  
    std::vector<int> bad_small_peak,marked;
    bad_small_peak.clear();
    marked.clear();
  
    for(size_t peak = 0; peak < MaxStartPoint.size(); ++peak){
      for(size_t peak2 = 0; peak2 < MaxStartPoint.size(); ++peak2){
  
        if((sigType == geo::kInduction && 
	    peak != peak2 && 
	    abs(MaxStartPoint[peak]-MaxEndPoint[peak2]) <= 1 && 
	    (fh_theta_ind_Area->GetBinContent(FinalPeaks[peak])  < 0.6 || 
	     fh_theta_ind_Area->GetBinContent(FinalPeaks[peak2]) < 0.6)) || 
	   (sigType == geo::kCollection && 
	    peak != peak2 && 
	    abs(MaxStartPoint[peak]-MaxEndPoint[peak2]) <= 1 && 
	    (fh_theta_coll_Area->GetBinContent(FinalPeaks[peak]) < 0.6 || 
	     fh_theta_coll_Area->GetBinContent(FinalPeaks[peak2]) < 0.6))
	   ) {
  
	  //get rid of one of them
	  if(std::find(marked.begin(), marked.end(),peak)  == marked.end() && 
	     std::find(marked.begin(), marked.end(),peak2) == marked.end()){
	    bad_small_peak.push_back(FinalPeaks[peak]);
	    marked.push_back(FinalPeaks[peak]);
	    marked.push_back(FinalPeaks[peak2]);
	  }//if not analyzed already
  
        }
      }
    }
  
    //now copy the right peaks, if we found any small peak ranges right next to each other:
    if(bad_small_peak.size()>0){
      mf::LogInfo("KingaCluster")<<" ATTENTION: WILL NEED TO DELETE PEAKS AT THE FOLLOWING BIN #s"
				 << ", b/c its range is right next to some other range and the "
				 << "peak signal is < 0.6 ";
  
      for(size_t bin = 0; bin < bad_small_peak.size(); ++bin){
        mf::LogInfo("KingaCluster") <<bin;
      }
  
      TempFinalPeaks.clear();
      TempMaxStartPoint.clear();
      TempMaxEndPoint.clear();
  
      for(size_t i = 0; i < FinalPeaks.size(); ++i){
        if(std::find(bad_small_peak.begin(),bad_small_peak.end(),FinalPeaks[i])==bad_small_peak.end()){
	  TempFinalPeaks.push_back(FinalPeaks[i]);
	  TempMaxStartPoint.push_back(MaxStartPoint[i]);
	  TempMaxEndPoint.push_back(MaxEndPoint[i]);
        }
      }
  
      bad_small_peak.clear();
      marked.clear();
  
      FinalPeaks.clear();
      MaxStartPoint.clear();
      MaxEndPoint.clear();
      FinalPeaks=TempFinalPeaks;
      MaxStartPoint=TempMaxStartPoint;
      MaxEndPoint=TempMaxEndPoint;
      TempFinalPeaks.clear();
      TempMaxStartPoint.clear();
      TempMaxEndPoint.clear();
    }//if need to make corrections
   
    need_to_reassign_hitsIDs=0;
   
    mf::LogInfo("KingaCluster")<<"FORMING CLUSTERS NOW :) ";
    mf::LogInfo("KingaCluster")<<"FinalPeaks are at bin(s):  ";
    for(size_t i = 0; i < FinalPeaks.size(); ++i){
      mf::LogInfo("KingaCluster")<<FinalPeaks[i]<<" which corresponds to angle=  "
				 <<-180+2*FinalPeaks[i];
    }
  
  
    //int no_noise_hits=0;
    mf::LogInfo("KingaCluster")<<"In FindClusters(int plane), we should be producing "
			       <<MaxStartPoint.size()<<" clusters";
    double a_polar, b_polar,theta_polar;
   
    std::vector<double> DiffAngles;
  
    for(size_t i = 0; i < allhits.size(); ++i){
     
     
      int diff_w= allhits[i]->WireID().Wire - fwire_vertex_reco[plane];
      b_polar = diff_w*0.4; /**in cm*/
      //b_polar = (w - fwire_vertex[plane])* 0.4; /**in cm*/
      a_polar = (allhits[i]->PeakTime() - ftime_vertex_reco[plane])* ftimetick *fdriftvelocity; /** in cm*/
      theta_polar = std::abs(std::asin(a_polar/std::sqrt(pow(a_polar,2)+pow(b_polar,2)))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/
     
      if     (b_polar == 0 && a_polar == 0) theta_polar =  90; /** in deg*/
      else if(b_polar == 0 && a_polar <  0) theta_polar = 180; /** in deg*/
      else if(b_polar == 0 && a_polar >  0) theta_polar =   0; /** in deg*/
      else if(b_polar >  0 && a_polar == 0) theta_polar =  90; /** in deg*/
      else if(b_polar <  0 && a_polar == 0) theta_polar = -90; /** in deg*/
      else if(b_polar >  0 && a_polar >  0) theta_polar = 90 - theta_polar; /** in deg*/
      else if(b_polar >  0 && a_polar <  0) theta_polar = 90 + theta_polar; /** in deg*/
      else if(b_polar <  0 && a_polar >  0) theta_polar = -(90 - theta_polar); /** in deg*/
      else if(b_polar <  0 && a_polar <  0) theta_polar = -(90 + theta_polar); /** in deg*/
        
      for(size_t ClusterNo = 0; ClusterNo < MaxStartPoint.size(); ++ClusterNo){
  
        if(theta_polar >= (-180+2*MaxStartPoint[ClusterNo]) && 
	   theta_polar <= (-180+2*MaxEndPoint[ClusterNo])){
	  //want to start counting from 1, O is reserved for hits that will be marked as noise
	  HitsWithClusterID.push_back(ClusterNo+1);
	  break;
        }
        else if(ClusterNo == MaxStartPoint.size()-1){
	  //decide where noise hits go
       
	  for(unsigned int peakNo=0;peakNo<FinalPeaks.size();peakNo++){
	    DiffAngles.push_back(std::abs(-180+2*FinalPeaks[peakNo]-theta_polar));
	  }
	  // now take minimum of DiffAngles and find at which position 
	  // it is at, this position corresponds to clusterNo +1 , 
	  // because we don't want to mark hits with zero 
       
	  int position = std::distance(DiffAngles.begin(),std::min_element(DiffAngles.begin(),DiffAngles.end()));
     
	  HitsWithClusterID.push_back(position+1);
	  // if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
	  //    mf::LogInfo("KingaCluster")<<"  This hit is closest to cluster # "<<position+1<<std::endl;
	  //    }
	  //no_noise_hits++;
	  DiffAngles.clear();
        }
      } //loop over all ranges
    } //allhits
  
  
    //mf::LogInfo("KingaCluster")<<"In FindClusters(int plane), marked "<<no_noise_hits<<" hits as NOISE"<<std::endl;
    //mf::LogInfo("KingaCluster")<<"HitsWithClusterID contains the following clusterIDs:"<<std::endl;
   
    //for(int i=0; i<HitsWithClusterID.size();i++){
  
    //mf::LogInfo("KingaCluster")<<HitsWithClusterID[i]<<"  ";
  
  
    //}
  
    //mf::LogInfo("KingaCluster")<<std::endl;
  
    //Now let's pick a minimum number of hits that we require in a cluster, 
    // MinHitsInCluster. If just formed cluster containes less than the desired 
    // number than it need to be assigned to the nearest cluster according to 
    // distance. This is done by removing the peak that corresponded to that 
    // cluster and rerunning assignement of hits again. We want to remove all 
    // the wrong peaks at once so reassignement is done only once. Put all the 
    // wrong peaks into a vector:
  
    std::vector<unsigned int> WrongPeakNo;
    std::vector<int> WireNo;
    int span=0;
  
    /// \todo this variable should be a parameter
    int MinHitsInCluster=3; 
  
    for(size_t NClus = 0; NClus < MaxStartPoint.size(); ++NClus){
      //search for clusters with too little hits (ie 1 or less than your desired parameter):
  
      int NoHitsInCluster = std::count(HitsWithClusterID.begin(),HitsWithClusterID.end(),NClus+1);
      mf::LogInfo("KingaCluster") <<"*** No of Hits for cluster # "
				  <<NClus+1<<" = "<<NoHitsInCluster;
  
      //........Check the span for small clusters................
      /// \todo 5 is now a magic number
      if(NoHitsInCluster < 5 && NoHitsInCluster > 0){
      
        WireNo.clear();
        for(size_t h = 0; h < HitsWithClusterID.size(); ++h) {
      
          if(HitsWithClusterID[h] == NClus+1){
	    WireNo.push_back(allhits[h]->Channel());
          }
      
        }
  
        //now order WireNo and subtract first and last element to get the span in wire number:
        std::sort(WireNo.begin(),WireNo.end());
  
        span = WireNo[WireNo.size()-1]-WireNo[0];
  
        /// \todo and now 0.66 is magic
        if(span > (NoHitsInCluster+0.66*NoHitsInCluster)){
	  need_to_reassign_hitsIDs=1;
	  WrongPeakNo.push_back(NClus);
        }
      }//if clusters containing less than 5 hits
  
      //......end of checking span for small clusters...............
      //
      // let's give a chance for 2-hit-clusters IF it is very well 
      // separated from the rest of the peaks, say about 60 degrees 
      // which is 30 bins then do NOT reassign hits but leave this 
      // small 2-hit cluster as a valid one
  
      if(NoHitsInCluster < MinHitsInCluster && NoHitsInCluster != 2){
        need_to_reassign_hitsIDs=1;
        WrongPeakNo.push_back(NClus);
      }
    
    
      if(NoHitsInCluster == 2){
        for(size_t pk = 0; pk < FinalPeaks.size(); ++pk){
	  if(NClus != pk && abs(FinalPeaks[NClus]-FinalPeaks[pk]) < 30){ 
	    go_ahead_at_reassign=1;
	    break;
	  }
        }
     
        if(go_ahead_at_reassign == 1) WrongPeakNo.push_back(NClus);
        if(go_ahead_at_reassign == 0){
     
	  mf::LogInfo("KingaCluster")<<"FOUND A 2-HIT-CLUSTER WHICH IS VERY WELL "
				     << "SEPARATED FROM THE REST, I will leave it "
				     << "(This is peak at bin #"<<FinalPeaks[NClus];
        } 
     
      } //2-hit-clusters
      go_ahead_at_reassign=0;
  
    }//loop thru all the clusters
  
    if(need_to_reassign_hitsIDs == 1){
      //Now we need to change all the vectors containing info about peaks and subtract all the fake peaks:
      HitsWithClusterID.clear();
      std::vector<int> FinalPeaksTemporary;
      std::vector<int> MaxStartPointTemporary;
      std::vector<int> MaxEndPointTemporary;
  
      std::vector<unsigned int>::iterator iter;
      size_t initialNoPeaks = FinalPeaks.size();
      for(size_t f = 0; f < initialNoPeaks; ++f){
        iter = find(WrongPeakNo.begin(),WrongPeakNo.end(),f);
        if(iter == WrongPeakNo.end()){
	  FinalPeaksTemporary.push_back(FinalPeaks[f]);
	  MaxStartPointTemporary.push_back(MaxStartPoint[f]);
	  MaxEndPointTemporary.push_back(MaxEndPoint[f]);
        }
      }
      FinalPeaks.clear();
      MaxStartPoint.clear();
      MaxEndPoint.clear();
  
      FinalPeaks=FinalPeaksTemporary;
      MaxStartPoint=MaxStartPointTemporary;
      MaxEndPoint=MaxEndPointTemporary;
    }//if need to reassign
  }


} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(KingaCluster)
  
} 

