////////////////////////////////////////////////////////////////////////
//
// LineMerger class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// biagio.rossi@lhep.unibe.ch
// msoderbe@syr.edu
// joshua.spitz@yale.edu
//
// This algorithm is designed to merge 2D lines with similar slope and endpoints 
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



//#ifndef LINEMERGER_H
//#define LINEMERGER_H


namespace cluster {
   
  class LineMerger : public art::EDProducer {
    
  public:
    
    explicit LineMerger(fhicl::ParameterSet const& pset); 
    ~LineMerger();
    
    void produce(art::Event& evt);
    void beginJob();
    
  private:
        
    std::string     fClusterModuleLabel;
    double          fSlope; // tolerance for matching angles between two lines (in units of radians) 
    double          fEndpointWindow; // tolerance for matching endpoints (in units of time samples) 
   
    bool SlopeCompatibility(double slope1,double slope2);
    int  EndpointCompatibility(std::vector<double> sclstart, std::vector<double> sclend,std::vector<double> cl2start, std::vector<double> cl2end);
    
  protected: 
    
  }; // class LineMerger

}

//#endif // LINEMERGER_H



namespace cluster{

  //-------------------------------------------------
  LineMerger::LineMerger(fhicl::ParameterSet const& pset) 
    : fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel"))
    , fSlope             (pset.get<double     >("Slope"))
    , fEndpointWindow    (pset.get<double     >("EndpointWindow"))
  {
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  }

  //-------------------------------------------------
  LineMerger::~LineMerger()
  {
  }

  //-------------------------------------------------
  void LineMerger::beginJob()
  {
    //this doesn't do anything now, but it might someday
  }
    
  //------------------------------------------------------------------------------------//
  void LineMerger::produce(art::Event& evt)
  { 
    // Get a Handle for the input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel,clusterVecHandle);

    art::ServiceHandle<geo::Geometry> geo;
    int nplanes = geo->Nplanes();

    //one PtrVector for each plane in the geometry
    std::vector< art::PtrVector<recob::Cluster> > Cls(nplanes);

    //vector with indicators for whether a cluster has been merged already
    std::vector< std::vector<int> > Cls_matches(nplanes);

    // loop over the input Clusters
    for(size_t i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a art::Ptr to each Cluster
      art::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      
      switch(cl->View()){
      case geo::kU :
	Cls[0].push_back(cl);
	Cls_matches[0].push_back(0);
	break;
      case geo::kV :
	Cls[1].push_back(cl);
	Cls_matches[1].push_back(0);
	break;
      case geo::kZ :
	Cls[2].push_back(cl);
	Cls_matches[2].push_back(0);
	break;
      default :
	break;
      }// end switch on view
    }// end loop over input clusters

    std::unique_ptr<std::vector<recob::Cluster> >             SuperClusters(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

    for(int i = 0; i < nplanes; ++i){

      int clustersfound = 0; // how many merged clusters found in each plane
      int clsnum1       = 0;

      for(size_t c = 0; c < Cls[i].size(); ++c){
	if(Cls_matches[i][clsnum1] == 1){
	  ++clsnum1;
	  continue;
	}

	art::FindManyP<recob::Hit> fmh(Cls[i], evt, fClusterModuleLabel);
	
	// make a new cluster to put into the SuperClusters collection 
	// because we want to be able to adjust it later
	recob::Cluster cl1( *( Cls[i][c].get() ) ); 

	// find the hits associated with the current cluster
	std::vector< art::Ptr<recob::Hit> > ptrvs = fmh.at(c);

	Cls_matches[i][clsnum1] = 1; 
	++clustersfound;
	
	int clsnum2 = 0;
	for(size_t c2 = 0; c2 < Cls[i].size(); ++c2){

	  recob::Cluster cl2( *(Cls[i][c2].get()) );

	  if(Cls_matches[i][clsnum2] == 1){
	    ++clsnum2;
	    continue;
	  }

	  // find the hits associated with this second cluster
	  std::vector< art::Ptr<recob::Hit> > ptrvs2 = fmh.at(c2);
	  
	  // check that the slopes are the same
	  // added 13.5 ticks/wirelength in ArgoNeuT. 
	  // \todo need to make this detector agnostic
	  // would be nice to have a LArProperties function that returns ticks/wire.
	  bool sameSlope = SlopeCompatibility(cl1.dTdW()*(1./13.5),
					      cl2.dTdW()*(1./13.5));  
	  
	  // check that the endpoints fall within a circular window of each other 
	  // done in place of intercept matching
	  int sameEndpoint = EndpointCompatibility(cl1.StartPos(), cl1.EndPos(),
						   cl2.StartPos(), cl2.EndPos());
	  
	  // if the slopes and end points are the same, combine the clusters
	  // note that after 1 combination cl1 is no longer what we started 
	  // with
	  if(sameSlope && sameEndpoint){
	    cl1 = cl1 + cl2;
	    Cls_matches[i][clsnum2] = 1;       

	    // combine the hit collections
	    // take into account order when merging hits from two clusters: doc-1776
	    if (sameEndpoint == 1) for(size_t h = 0; h < ptrvs2.size(); ++h) ptrvs.push_back(ptrvs2[h]);
	    else if (sameEndpoint == -1){
	      for(size_t h = 0; h < ptrvs.size(); ++h) ptrvs2.push_back(ptrvs[h]);
	      ptrvs = ptrvs2;
	    }
	  }
	  
	  ++clsnum2;
	}// end loop over second cluster iterator

	// now add the final version of cl1 to the collection of SuperClusters
	// and create the association between the super cluster and the hits
	SuperClusters->push_back(cl1);
	util::CreateAssn(*this, evt, *(SuperClusters.get()), ptrvs, *(assn.get()));	
	++clsnum1;

      }// end loop over first cluster iterator
    }// end loop over planes

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "LineMerger Summary:";
    for(size_t i = 0; i < SuperClusters->size(); ++i) 
      mf::LogVerbatim("Summary") << SuperClusters->at(i);

    evt.put(std::move(SuperClusters));
    evt.put(std::move(assn));

    return;

  }

  //------------------------------------------------------------------------------------//
  //checks the difference between angles of the two lines
  bool LineMerger::SlopeCompatibility(double slope1, double slope2)
  { 
    double sl1 = atan(slope1);
    double sl2 = atan(slope2);

    //the units of fSlope are radians
    bool comp  = std::abs(sl1-sl2) < fSlope ? true : false;

    return comp;
  }
  //------------------------------------------------------------------------------------//
  int LineMerger::EndpointCompatibility(std::vector<double> sclstart, 
					 std::vector<double> sclend,
					 std::vector<double> cl2start, 
					 std::vector<double> cl2end)
  { 
    double sclstartwire = sclstart[0];
    double sclstarttime = sclstart[1];
    double sclendwire   = sclend[0];
    double sclendtime   = sclend[1];
    
    double cl2startwire = cl2start[0];
    double cl2starttime = cl2start[1];
    double cl2endwire   = cl2end[0];
    double cl2endtime   = cl2end[1];

    // \todo 13.5 ticks/wire. need to make this detector agnostic--spitz
    double distance = std::sqrt((pow(sclendwire-cl2startwire,2)*13.5) + pow(sclendtime-cl2starttime,2));

    //not sure if this line is necessary--spitz
    double distance2 = std::sqrt((pow(sclstartwire-cl2endwire,2)*13.5) + pow(sclstarttime-cl2endtime,2));
    
//    bool comp = (distance  < fEndpointWindow ||
//		 distance2 < fEndpointWindow) ? true : false;

    //determine which way the two clusters should be merged. TY
    int comp = 0;
    if (distance < fEndpointWindow) 
      comp = 1;
    else if (distance2 < fEndpointWindow)
      comp = -1;
    return comp;
  }



} // end namespace





namespace cluster{

  DEFINE_ART_MODULE(LineMerger)
  
} // end namespace 

