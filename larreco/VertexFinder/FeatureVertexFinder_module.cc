////////////////////////////////////////////////////////////////////////
//
// FeatureVertexFinder class
//
// jasaadi@fnal.gov
//
// This algorithm is designed to reconstruct the vertices using the
// 2D cluster information and CornerFinderAlg to find 2-d and 3-d verticies
//
// 06/05/14
// This is a major rewrite of the original code and will be factored to take
// advantage of improvements in LArSoft. The code now goes something like this:
//
// 1)Take up EndPoint2d's from Cluster Crawler and keep ones that make sense in 3d
// (this algorithm wants ClusterCrawler to have been run, but should totally work
//  even if it hasn't)
//
// 2)Take up all EndPoint2d's from the image finding CornerFinder algorithm
// (this algortithm wants CornerFinder to have been run, but should totally work
//  even if it hasn't)
//
// 3)Take up an additional clustering algorithm (dBcluster by default but can work 
//   with any clustering algorithm) and use slopes and endpoints to find intersections
//
// 4)Now merge together any duplicates and sort them by their weight (which needs its units still
//   decided somehow)
//
// 5)Return a list of verticies based on what mode you run this algorithm in
//	a) Primary mode: This mode will return only a single 3d vertex and an appropriate
//	   		 number of EndPoint2d's corresponding to the most likely neutrino vertex
//	b) All mode: This mode returns all 3d/2d verticies that this algorithm has found but
//		     sorts the list by the most likely candidate
// 
//
//
// This is Preliminary Work and needs modifications
//
// ////////////////////////////////////////////////////////////////////////

// ##########################
// ### Basic C++ Includes ###
// ##########################
#include <string>
#include <iostream>
#include <iomanip>
#include <ios>
#include <fstream>
#include <cmath>
#include <algorithm>

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
//#include "RecoAlg/ClusterParamsAlg.h"

// #####################
// ### ROOT Includes ###
// #####################
#include "TH1D.h"
#include "TVectorD.h"
#include "TGeoManager.h"
#include "TGraph.h"
#include "TF1.h"

// =====================================================================================================
// =====================================================================================================

class TH1D;

///vertex reconstruction
namespace vertex {
   
 class FeatureVertexFinder :  public art::EDProducer {
    
  public:
    
    explicit FeatureVertexFinder(fhicl::ParameterSet const& pset); 

  private:
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt);
   
   // ### This function will take in and EndPoint2d from either cluster crawler
   // ### or corner finder and only save points that make a 3d-candidate 
   void Get3dVertexCandidates(std::vector< art::Ptr<recob::EndPoint2D> > EndPoints, bool PlaneDet);
   
   // ### This function will take in 2d Clusters (by default dBCluster)
   // ### and return 2d vertex candidates for sorting later
   void Find2dClusterVertexCandidates(art::PtrVector<recob::Cluster> RawClusters, art::FindManyP<recob::Hit> fmhit);
   
   // ### This function takes in 2d Vertex candidates (found by Find2dClusterVertexCandidates), 
   // ### does some basic merging and then finds 3d-candidates for use later
   void Find3dVtxFrom2dClusterVtxCand(std::vector<double> Wire_2dvtx, std::vector<double> Time_2dvtx, std::vector<double> Plane_2dvtx);
   
   // ### This function merges and sorts all the 3d vertex candidates, it merges 
   // ### them if they are within 2.0 cm in X, Y, and Z simultaneously and then sorts
   // ### the list by vertex strength and Z location (lowest Z is first on the list
   void MergeAndSort3dVtxCandidate(std::vector<double> merge_vtxX, std::vector<double> merge_vtxY, std::vector<double> merge_vtxZ, std::vector<double> merge_vtxStgth);
   
  
   std::string fClusterModuleLabel;
   std::string fHitModuleLabel;
   std::string fCornerFinderModuleLabel;
   std::string fCCrawlerEndPoint2dModuleLabel;
   //cluster::ClusterParamsAlg fClParAlg;
   Double_t fRunningMode;
   
   bool GT2PlaneDetector = false;
   
   // ########################################################
   // ### Unsorted Raw list of 3d-Vertex candidates filled ###
   // ########################################################
   std::vector<double> candidate_x = {0.};
   std::vector<double> candidate_y = {0.};
   std::vector<double> candidate_z = {0.};
   std::vector<double> candidate_strength = {0.};
   
   
   // #############################################################
   // ### Merged and sorted list of 3d-Vertex candidates filled ###
   // #############################################################
   std::vector<double> MergeSort3dVtx_xpos = {0.};
   std::vector<double> MergeSort3dVtx_ypos = {0.};
   std::vector<double> MergeSort3dVtx_zpos = {0.};
   std::vector<double> MergeSort3dVtx_strength = {0.};
   
   
   // ##################################
   // ### 2d Cluster based Variables ###
   // ################################## 
   std::vector<std::vector<int> > Cls;		 	//<---- Index to clusters in each view
   std::vector<double> dtdwstart		= {0.};	//<----Slope (delta Time Tick vs delta Wire) 
   std::vector<double> dtdwstartError 		= {0.};	//<---Error on the slope 
   std::vector<double> Clu_Plane 		= {0.};	//<---Plane of the current cluster
   std::vector<double> Clu_StartPos_Wire 	= {0.};	//<---Starting wire number of cluster
   std::vector<double> Clu_StartPos_TimeTick 	= {0.};	//<---Starting TDC value of the cluster
  
   std::vector<double> Clu_EndPos_Wire 		= {0.};	//<---Ending wire number of cluster
   std::vector<double> Clu_EndPos_TimeTick 	= {0.}; //<---Ending TDC value of the cluster
  
   std::vector<double> Clu_Slope 		= {0.};	//<---Calculated Slope of the cluster (TDC/Wire)
   std::vector<double> Clu_Yintercept 		= {0.};	//<---Clusters Y Intercept using start positions
   std::vector<double> Clu_Yintercept2 		= {0.};	//<---Clusters Y Intercept using end positions
   std::vector<double> Clu_Length 		= {0.}; //<---Calculated Length of the cluster
    
   
   // ##################################
   // ### 2d Cluster based verticies ###
   // ##################################
   std::vector<double> TwoDvtx_wire = {0.};
   std::vector<double> TwoDvtx_time = {0.};
   std::vector<double> TwoDvtx_plane = {0.};
  };
    
}//<---End namespace vertex

// =====================================================================================================
// =====================================================================================================







// =====================================================================================================
// =====================================================================================================
namespace vertex{

//-----------------------------------------------------------------------------
// fhicl::ParameterSet
  FeatureVertexFinder::FeatureVertexFinder(fhicl::ParameterSet const& pset) :
    EDProducer{pset}
  //fClParAlg(pset.get<fhicl::ParameterSet>("ClusterParamsAlg"), pset.get< std::string >("module_type"))
  {  
    /*this->*/reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::EndPoint2D, recob::Hit> >();
    
    // === Don't think I'll actually have any of these associations ===
    // ===         should consider removing in the future
    produces< art::Assns<recob::Vertex, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
    
    art::ServiceHandle<geo::Geometry const> geom;
    Cls.resize(geom->Nplanes(),std::vector<int>());
  }

//---------------------------------------------------------------------------
void FeatureVertexFinder::reconfigure(fhicl::ParameterSet const& p) 
  {
    fCornerFinderModuleLabel  	   = p.get< std::string >("CornerFinderModuleLabel");
    fClusterModuleLabel       	   = p.get< std::string >("ClusterModuleLabel");
    fHitModuleLabel	      	   = p.get< std::string >("HitModuleLabel");
    fCCrawlerEndPoint2dModuleLabel = p.get< std::string >("CCrawlerEndPoint2dModuleLabel");
    fRunningMode                   = p.get< double       >("RunningMode");
  }

// -----------------------------------------------------------------------------
// Produce
void FeatureVertexFinder::produce(art::Event& evt)
   {

   // #########################
   // ### Geometry Services ###
   // #########################
   art::ServiceHandle<geo::Geometry const> geom;
    
   // ####################################
   // ### Detector Properties Services ###
   // ####################################
   const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   
   // ######################################################
   // ### Figuring out if I have a 2 or 3 plane detector ###
   // ######################################################
   GT2PlaneDetector = false;
   // ##############################
   // ### Looping over cryostats ###
   // ##############################
   for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
   	{
    	// ##########################
      	// ### Looping over TPC's ###
      	// ##########################
      	for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
	   {
     	   if (geom->Cryostat(cstat).TPC(tpc).Nplanes() > 2){GT2PlaneDetector = true;}
	   }//<---End tpc loop
	}//<---End cstat loop
		
   if(GT2PlaneDetector){std::cout<<"yeah"<<std::endl;}
   
   // #######################################################
   // ### These are the things I want to put on the event ###
   // #######################################################
   std::unique_ptr<std::vector<recob::Vertex> >                 vcol(new std::vector<recob::Vertex>);          //3D vertex
   std::unique_ptr<std::vector<recob::EndPoint2D> >             epcol(new std::vector<recob::EndPoint2D>);  //2D vertex
   std::unique_ptr< art::Assns<recob::EndPoint2D, recob::Hit> > assnep(new art::Assns<recob::EndPoint2D, recob::Hit>);
   std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> >  assnsh(new art::Assns<recob::Vertex, recob::Shower>);
   std::unique_ptr< art::Assns<recob::Vertex, recob::Track> >   assntr(new art::Assns<recob::Vertex, recob::Track>);
   std::unique_ptr< art::Assns<recob::Vertex, recob::Hit> >     assnh(new art::Assns<recob::Vertex, recob::Hit>);  
   
   //std::cout<<"Setup of outputs"<<std::endl;
   
//-------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------- ClusterCrawler EndPoint2d's -------------------------------------------------------------  
//-------------------------------------------------------------------------------------------------------------------------------------------------
   // ####################################################
   // ### Getting the EndPoint2d's for cluster crawler ###
   // ####################################################
   try
      {
      art::Handle< std::vector<recob::EndPoint2D> > ccrawlerFinderHandle;
      evt.getByLabel(fCCrawlerEndPoint2dModuleLabel,ccrawlerFinderHandle);
      std::vector< art::Ptr<recob::EndPoint2D> > ccrawlerEndPoints;
      art::fill_ptr_vector(ccrawlerEndPoints, ccrawlerFinderHandle);
   
      //std::cout<<"Getting the EndPoint2d's for cluster crawler"<<std::endl;
      //std::cout<<"Length of ccrawlerEndPoints vector = "<<ccrawlerEndPoints.size()<<std::endl;
      // ########################################################
      // ### Passing in the EndPoint2d's from Cluster Crawler ###
      // ########################################################
      Get3dVertexCandidates(ccrawlerEndPoints, GT2PlaneDetector);
      //std::cout<<"Get3dVertexCandidates Cluster Crawler"<<std::endl;
      //std::cout<<"Number of candidate verticies after cluster crawler = "<<candidate_x.size()<<std::endl;
      }
   catch(...)
      {mf::LogWarning("FeatureVertexFinder") << "Failed to get EndPoint2d's from Cluster Crawler";}
//-------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------- CornerFinder EndPoint2d's --------------------------------------------------------------  
//-------------------------------------------------------------------------------------------------------------------------------------------------   
   // ##################################################
   // ### Getting the EndPoint2d's for Corner Finder ###
   // ##################################################
   try
      {
      art::Handle< std::vector<recob::EndPoint2D> > CornerFinderHandle;
      evt.getByLabel(fCornerFinderModuleLabel,CornerFinderHandle);
      std::vector< art::Ptr<recob::EndPoint2D> > cornerEndPoints;
      art::fill_ptr_vector(cornerEndPoints, CornerFinderHandle);
   
      //std::cout<<"Getting the EndPoint2d's for Corner Finder"<<std::endl;
      //std::cout<<"Length of cornerEndPoints vector = "<<cornerEndPoints.size()<<std::endl;
   
      // ######################################################
      // ### Passing in the EndPoint2d's from Corner Finder ###
      // ######################################################
      Get3dVertexCandidates(cornerEndPoints, GT2PlaneDetector);
      //std::cout<<"Get3dVertexCandidates Corner Finder"<<std::endl;
      //std::cout<<"Number of candidate verticies after corner finder = "<<candidate_x.size()<<std::endl;
      
      }
   catch(...)
      {mf::LogWarning("FeatureVertexFinder") << "Failed to get EndPoint2d's from Corner Finder";}


//---------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------- Making Cluster Slope EndPoint2d's --------------------------------------------------------------  
//--------------------------------------------------------------------------------------------------------------------------------------------------------- 
   // ###################################################
   // ### Retreiving the Cluster Module for the event ###
   // ###################################################
   try
      {
      art::Handle< std::vector<recob::Cluster> > clusterListHandle;
      evt.getByLabel(fClusterModuleLabel,clusterListHandle);
   
      //std::cout<<"Retreiving the Cluster Module for the event"<<std::endl;
      
   
      // #################################################
      // ### Finding hits associated with the clusters ###
      // #################################################
      art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);
      //std::cout<<"Finding hits associated with the clusters"<<std::endl;
   
      // ##################################
      // ### Filling the Cluster Vector ###
      // ##################################
      art::PtrVector<recob::Cluster> clusters;
      for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
     	art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
     	clusters.push_back(clusterHolder);
     	}//<---End ii loop
      
      //std::cout<<"Number of clusters found = "<<clusters.size()<<std::endl;
  
      // ####################################################
      // ### Passing in the clusters to find 2d Verticies ###
      // ####################################################
      Find2dClusterVertexCandidates(clusters, fmh);
      //std::cout<<"Made 2d vertex candidates"<<std::endl;
      //std::cout<<"Number of 2d cluster vertex candidates found = "<<TwoDvtx_wire.size()<<std::endl;
      // ###############################################################
      // ### Finding 3d Candidates from 2d cluster vertex candidates ###
      // ###############################################################
      Find3dVtxFrom2dClusterVtxCand(TwoDvtx_wire, TwoDvtx_time, TwoDvtx_plane);
      //std::cout<<"Made 3d vertex candidates from 2d cluster candidates"<<std::endl;
      //std::cout<<"Number of candidate verticies after cluster step = "<<candidate_x.size()<<std::endl;
      
      }
   catch(...)
      {mf::LogWarning("FeatureVertexFinder") << "Failed to get Cluster from default cluster module";}
  
  
   // ################################################
   // ### Merging and sorting 3d vertex candidates ###
   // ################################################
   MergeAndSort3dVtxCandidate(candidate_x, candidate_y, candidate_z, candidate_strength);  
   /*std::cout<<std::endl;
   std::cout<<"Merged and sorted the cadidates"<<std::endl;
   std::cout<<"Number of merged and sorted verticies = "<<MergeSort3dVtx_xpos.size()<<std::endl;
   std::cout<<"MergeSort3dVtx_xpos[0] = "<<MergeSort3dVtx_xpos[0]<<std::endl;
   std::cout<<"MergeSort3dVtx_ypos[0] = "<<MergeSort3dVtx_ypos[0]<<std::endl;
   std::cout<<"MergeSort3dVtx_zpos[0] = "<<MergeSort3dVtx_zpos[0]<<std::endl;*/


//------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------- Putting Verticies on the event --------------------------------------------------------------  
//------------------------------------------------------------------------------------------------------------------------------------------------------  
   
   // ########################################################################################################
   // ###           Now that I have a list of 3d vertex candidates I will return 3d/2d verticies           ###
   // ###                           based on which option the user has chosen                              ###
   // ### fRunningMode == 0 (this returns a full list of all 3d/2d vertex candidates)                      ###
   // ### fRunningMode == 1 (this returns only one vertex which is established as the most likely primary) ###
   // ########################################################################################################
   
   // =======================================
   // === Returning all vertex candidates ===
   // =======================================
   if(fRunningMode == 0)
      {
      // ######################################
      // ### Looping over Primary Verticies ###
      // ######################################
      for(size_t pri = 0; pri < MergeSort3dVtx_xpos.size(); pri++)
	{
	// ###############################################
	// ### Push each primary vertex onto the event ###
	// ###############################################
	double tempxyz[3] = {MergeSort3dVtx_xpos[pri], MergeSort3dVtx_ypos[pri], MergeSort3dVtx_zpos[pri]};
	// ######################################
	// ### Skipping a vertex that is zero ###
	// ######################################
	if(tempxyz[0] == 0 && tempxyz[1] == 0 && tempxyz[2] == 0){continue;}
	recob::Vertex the3Dvertex(tempxyz, vcol->size());
	vcol->push_back(the3Dvertex);
	// ---------------------------------------------------------------------
	// --- Now go make the 2DEndPoints that correspond to each 3d vertex ---
	// ---------------------------------------------------------------------
		
	// ##############################
    	// ### Looping over cryostats ###
    	// ##############################
    	for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
    	   {
      	   // ### Looping over TPC's ###
      	   for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
	      {
	      // ### Loop over the wire planes ###
	      for (size_t plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane)
	      	{
		double temp2dXYZ[3] = {MergeSort3dVtx_xpos[pri], MergeSort3dVtx_ypos[pri], MergeSort3dVtx_zpos[pri]};
		double temp2dStrength = MergeSort3dVtx_strength[pri];
		// ######################################
		// ### Skipping a vertex that is zero ###
		// ######################################
		if(temp2dXYZ[0] == 0 && temp2dXYZ[1] == 0 && temp2dXYZ[2] == 0){continue;}
		
		// ######################################################################
		// ### Converting the 3d vertex into 2d time ticks, wire, and channel ###
		// ######################################################################
		double EndPoint2d_TimeTick = detprop->ConvertXToTicks(temp2dXYZ[0],plane, tpc, cstat);
		int EndPoint2d_Wire = 0;
		int EndPoint2d_Channel = 0;
		// ### Putting in protection in case NearestWire Fails ###
		try
		   {EndPoint2d_Wire = geom->NearestWire(temp2dXYZ , plane, tpc, cstat);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      	   continue;}
		// ### Putting in protection in case NearestChannel Fails ###
		try
		   {EndPoint2d_Channel     = geom->NearestChannel(temp2dXYZ, plane, tpc, cstat);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      	    continue;}
		    
		 // ### Making geo::WireID and getting the current View number ###	
		geo::View_t View = geom->View(EndPoint2d_Channel);
		geo::WireID wireID(cstat,tpc,plane,EndPoint2d_Wire);
		
		// ################################################
		// ### Putting the 2d Vertex found on the event ###
		// ################################################
		recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
					  wireID ,		//<---geo::WireID
					  temp2dStrength ,	//<---Vtx strength (JA: ?)
					  epcol->size() ,	//<---Vtx ID (JA: ?)
					  View ,		//<---Vtx View 	
					  1 );			//<---Total Charge (JA: Need to figure this one?)
		epcol->push_back(vertex);
		}//<---End Plane loop
      	      }//<---End TPC loop
	   }//<---End cstat loop
	}//<---End pri loop
      }//<---End fRunningMode == 0
      
   
   // ================================================
   // === Returning only primary vertex candidates ===
   // ================================================
   if(fRunningMode != 0)
      {
      int position = 0;
      int bail = 0;
      // ######################################
      // ### Looping over Primary Verticies ###
      // ######################################
      for(size_t pri = 0; pri < MergeSort3dVtx_xpos.size(); pri++)
	{
	// ###############################################
	// ### Push each primary vertex onto the event ###
	// ###############################################
	double tempxyz[3] = {MergeSort3dVtx_xpos[pri], MergeSort3dVtx_ypos[pri], MergeSort3dVtx_zpos[pri]};
	// ######################################
	// ### Skipping a vertex that is zero ###
	// ######################################
	if(bail > 0){continue;}
	if(tempxyz[0] == 0 && tempxyz[1] == 0 && tempxyz[2] == 0){continue;}
	position = pri;
	bail++;
	recob::Vertex the3Dvertex(tempxyz, vcol->size());
	vcol->push_back(the3Dvertex);
	
	
	}
	
	// ---------------------------------------------------------------------
	// --- Now go make the 2DEndPoints that correspond to each 3d vertex ---
	// ---------------------------------------------------------------------
		
	// ##############################
    	// ### Looping over cryostats ###
    	// ##############################
    	for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
    	   {
      	   // ### Looping over TPC's ###
      	   for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
	      {
	      // ### Loop over the wire planes ###
	      for (size_t plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane)
	      	{
		double temp2dXYZ[3] = {MergeSort3dVtx_xpos[position], MergeSort3dVtx_ypos[position], MergeSort3dVtx_zpos[position]};
		double temp2dStrength = MergeSort3dVtx_strength[position];
		
		// ######################################################################
		// ### Converting the 3d vertex into 2d time ticks, wire, and channel ###
		// ######################################################################
		double EndPoint2d_TimeTick = detprop->ConvertXToTicks(temp2dXYZ[0],plane, tpc, cstat);
		int EndPoint2d_Wire = 0;
		int EndPoint2d_Channel = 0;
		// ### Putting in protection in case NearestWire Fails ###
		try
		   {EndPoint2d_Wire = geom->NearestWire(temp2dXYZ , plane, tpc, cstat);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      	   continue;}
		// ### Putting in protection in case NearestChannel Fails ###
		try
		   {EndPoint2d_Channel     = geom->NearestChannel(temp2dXYZ, plane, tpc, cstat);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      	    continue;}
		    
		 // ### Making geo::WireID and getting the current View number ###	
		geo::View_t View = geom->View(EndPoint2d_Channel);
		geo::WireID wireID(cstat,tpc,plane,EndPoint2d_Wire);
		
		// ################################################
		// ### Putting the 2d Vertex found on the event ###
		// ################################################
		recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
					  wireID ,		//<---geo::WireID
					  temp2dStrength ,	//<---Vtx strength (JA: ?)
					  epcol->size() ,	//<---Vtx ID (JA: ?)
					  View ,		//<---Vtx View 	
					  1 );			//<---Total Charge (JA: Need to figure this one?)
		epcol->push_back(vertex);
		}//<---End Plane loop
      	      }//<---End TPC loop
	   }//<---End cstat loop
      }//<---End fRunningMode == 1   
   
   
    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "FeatureVertexFinder Summary:";
    for(size_t i = 0; i<epcol->size(); ++i) mf::LogVerbatim("Summary") << epcol->at(i) ;
    for(size_t i = 0; i<vcol->size(); ++i) mf::LogVerbatim("Summary") << vcol->at(i) ;
    
    
    /*for(size_t j = 0; j<epcol->size(); ++j) std::cout<<" EndPoint2d = " << epcol->at(j) ;
    std::cout<<std::endl;
    for(size_t j = 0; j<vcol->size(); ++j) {std::cout<< " Vertex 3d = " << vcol->at(j) << std::endl;}
    std::cout<<std::endl;*/
    
    
    
    evt.put(std::move(epcol));
    evt.put(std::move(vcol));
    evt.put(std::move(assnep));
    evt.put(std::move(assntr));
    evt.put(std::move(assnsh));
    evt.put(std::move(assnh));
    
    // ################################################
    // ### Clearing vectors at the end of the event ###
    // ################################################
    vcol.reset();
    epcol.reset();
    candidate_x.clear();
    candidate_y.clear();
    candidate_z.clear();
    candidate_strength.clear();
    MergeSort3dVtx_xpos.clear();
    MergeSort3dVtx_ypos.clear();
    MergeSort3dVtx_zpos.clear();
    MergeSort3dVtx_strength.clear();
    dtdwstart.clear(); 
    dtdwstartError.clear(); 
    Clu_Plane.clear();
    Clu_StartPos_Wire.clear();
    Clu_StartPos_TimeTick.clear();
    Clu_EndPos_Wire.clear();
    Clu_EndPos_TimeTick.clear();
    Clu_Slope.clear();
    Clu_Yintercept.clear();
    Clu_Yintercept2.clear();
    Clu_Length.clear();
    TwoDvtx_wire.clear();
    TwoDvtx_time.clear();
    TwoDvtx_plane.clear();
    
  }//<--End FeatureVertexFinder::produce





//==============================================================================================================================================================================
//==============================================================================================================================================================================
//==============================================================================================================================================================================
//==============================================================================================================================================================================





// -----------------------------------------------------------------------------
// Get 3d Vertex Candidates
// -----------------------------------------------------------------------------
void vertex::FeatureVertexFinder::Get3dVertexCandidates(std::vector< art::Ptr<recob::EndPoint2D> > EndPoints, bool PlaneDet)
   {
   // #########################
   // ### Geometry Services ###
   // #########################
   art::ServiceHandle<geo::Geometry const> geom;
   
   // ####################################
   // ### Detector Properties Services ###
   // ####################################
   const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   
   double y = 0., z = 0.;
   double yy = 0., zz = 0.;
   double yy2 = 0., zz2 = 0.;
   double yy3 = 0., zz3 = 0.;
   // ##############################
   // ### Looping over cryostats ###
   // ##############################
   for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
      {
       // ##########################
       // ### Looping over TPC's ###
       // ##########################
       for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
	{
	// ###############################################################################
    	// ### Now loop over those Endpoints and see if they match up in between views ###
    	// ###############################################################################
	for(size_t endpt1 = 0; endpt1 < EndPoints.size(); endpt1++)
	   {
	   // #########################################
	   // ### Looping over the rest of the list ###
	   // #########################################
	   for(size_t endpt2 = endpt1+1; endpt2 < EndPoints.size(); endpt2++)
	   	{
		// ##########################################################################
		// ### Check to make sure we are comparing features from different planes ###
		// ##########################################################################
		if(EndPoints.at(endpt1)->WireID().Plane != EndPoints.at(endpt2)->WireID().Plane)
		   {
		   // #############################################################################
		   // ### Get the appropriate time offset for the two planes we are considering ###
		   // #############################################################################
		   float tempXFeature1 = detprop->ConvertTicksToX(EndPoints.at(endpt1)->DriftTime(), EndPoints.at(endpt1)->WireID().Plane, tpc, cstat);		
		   float tempXFeature2 = detprop->ConvertTicksToX(EndPoints.at(endpt2)->DriftTime(), EndPoints.at(endpt2)->WireID().Plane, tpc, cstat);
		   // ####################################################################
		   // ### Checking to see if these features have intersecting channels ###
		   // ###             and are within 0.5 cm in projected X             ###
		   // ####################################################################
		   if( geom->ChannelsIntersect( geom->PlaneWireToChannel(EndPoints.at(endpt2)->WireID().Plane, EndPoints.at(endpt2)->WireID().Wire, tpc, cstat), 
		                                geom->PlaneWireToChannel(EndPoints.at(endpt1)->WireID().Plane, EndPoints.at(endpt1)->WireID().Wire, tpc, cstat), 
						yy, zz) &&
					        std::abs(tempXFeature1 - tempXFeature2) < 0.5)
		      {
		      // #####################################################################################
		      // ### Use this fill if we are in a detector with fewer than 3 plane (e.g. ArgoNeuT) ###
		      // #####################################################################################
		      if(!PlaneDet)
		      	{
		      	candidate_x.push_back(tempXFeature1);
		      	candidate_y.push_back(yy);
		      	candidate_z.push_back(zz);
		      	candidate_strength.push_back(EndPoints.at(endpt1)->Strength() + EndPoints.at(endpt2)->Strength());
			}//<---End fill for 2 plane detector
		      // ##############################################################################
		      // ### Adding a check to see if I am in a 3-plane detector and therefore need ###
		      // ###           to check for a match across more than 2 planes               ###
		      // ##############################################################################
		      if(PlaneDet)
		      	{
			// #########################################
	   		// ### Looping over the rest of the list ###
	   		// #########################################
	   		for(size_t endpt3 = endpt2+1; endpt3 < EndPoints.size(); endpt3++)
	   		   {
			   // ##########################################################################
			   // ### Check to make sure we are comparing features from different planes ###
			   // ##########################################################################
			   if(EndPoints.at(endpt3)->WireID().Plane != EndPoints.at(endpt2)->WireID().Plane && 
			      EndPoints.at(endpt3)->WireID().Plane != EndPoints.at(endpt1)->WireID().Plane && 
			      EndPoints.at(endpt1)->WireID().Plane != EndPoints.at(endpt2)->WireID().Plane )
			      {
			      float tempXFeature3 = detprop->ConvertTicksToX(EndPoints.at(endpt3)->DriftTime(), EndPoints.at(endpt3)->WireID().Plane, tpc, cstat);
			      // ####################################################################################
			      // ### Checking to make sure our third feature has an intersecting channel with our ###
			      // ###         other two channels and is within 1.0 cm projected in X               ###
			      // ####################################################################################
			      if(geom->ChannelsIntersect( geom->PlaneWireToChannel(EndPoints.at(endpt3)->WireID().Plane, EndPoints.at(endpt3)->WireID().Wire, tpc, cstat), 
		                                          geom->PlaneWireToChannel(EndPoints.at(endpt1)->WireID().Plane, EndPoints.at(endpt1)->WireID().Wire, tpc, cstat), 
						          yy3, zz3) &&
				 geom->ChannelsIntersect( geom->PlaneWireToChannel(EndPoints.at(endpt3)->WireID().Plane, EndPoints.at(endpt3)->WireID().Wire, tpc, cstat), 
		                                          geom->PlaneWireToChannel(EndPoints.at(endpt2)->WireID().Plane, EndPoints.at(endpt2)->WireID().Wire, tpc, cstat), 
						          yy2, zz2) &&
				 geom->ChannelsIntersect( geom->PlaneWireToChannel(EndPoints.at(endpt2)->WireID().Plane, EndPoints.at(endpt2)->WireID().Wire, tpc, cstat), 
		                                          geom->PlaneWireToChannel(EndPoints.at(endpt1)->WireID().Plane, EndPoints.at(endpt1)->WireID().Wire, tpc, cstat), 
						          yy, zz) &&
				 std::abs(tempXFeature3 - tempXFeature2) < 1.0 && std::abs(tempXFeature3 - tempXFeature1) < 1.0 &&
				 std::abs(tempXFeature1 - tempXFeature2) < 1.0 )
			      	{
				candidate_x.push_back(detprop->ConvertTicksToX(EndPoints.at(endpt1)->DriftTime(), EndPoints.at(endpt1)->WireID().Plane, tpc, cstat));
				
				// ###################################
				// ### Finding intersection points ###
				// ###################################
				geom->IntersectionPoint(EndPoints.at(endpt1)->WireID().Wire, EndPoints.at(endpt2)->WireID().Wire, 
				                        EndPoints.at(endpt1)->WireID().Plane, EndPoints.at(endpt2)->WireID().Plane, cstat, tpc, y, z);
							
				candidate_y.push_back(y);
				candidate_z.push_back(z);
				candidate_strength.push_back(EndPoints.at(endpt1)->Strength()+EndPoints.at(endpt2)->Strength()+EndPoints.at(endpt3)->Strength());
				
				// ### Note: If I've made it here I have a matched triplet...since I don't want to use any of these ###
				// ###    features again I am going to iterate each of the counters so we move to the next one      ###
				if(endpt1 < EndPoints.size())
				   {endpt1++;}
				if(endpt2 < EndPoints.size())
				   {endpt2++;}
				if(endpt3 < EndPoints.size())
				   {endpt3++;}
				}//<---End finding 3d point across all three planes
			      
			      }//<---End checking for all different planes
			   }//<---End endpt3
			
			}//<---End fill for 3 plane detector
		      
		      }//<---End intersecting channels
		      
	   	   }//<---End making sure we are looking across planes
		}//<---End endpt2 loop
	   }//<---End endpt1 loop
	}//<---End TPC loop
      }//<---End cstat
   
   
   
   }//<---End Get3dVertexCandidates




// -----------------------------------------------------------------------------   
// Get 2d Vertex Candidates from clusters
// -----------------------------------------------------------------------------   
void vertex::FeatureVertexFinder::Find2dClusterVertexCandidates(art::PtrVector<recob::Cluster> RawClusters, art::FindManyP<recob::Hit> fmhit)
   {
   // #########################
   // ### Geometry Services ###
   // #########################
   art::ServiceHandle<geo::Geometry const> geom;
   
   // ####################################
   // ### Detector Properties Services ###
   // ####################################
   const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   
   int nClustersFound = 0;
   
   // Initialize Cls
   for(auto& c : Cls) c.clear();
      
   for(size_t iclu = 0; iclu < RawClusters.size(); ++iclu)
      {
      // ### Gathering the hits associated with the current cluster ###
      std::vector< art::Ptr<recob::Hit> > hit = fmhit.at(iclu);
   
      // ### I only want to consider this cluster if it has a sufficient number of hits ###
      if(hit.size() < 15){continue;}
         
      // ######################################################
      // ### Determine which view the current cluster is in ###
      // ######################################################
      const geo::View_t view = RawClusters.at(iclu)->View();
      if(view >= Cls.size()) {
          std::cerr << __PRETTY_FUNCTION__
	  	    <<"\033[93m Logic error in my code ... view " << view << " not supported ! \033[00m"
		    << std::endl;
          throw std::exception();      
      }
      
      Cls.at(RawClusters.at(iclu)->View()).push_back(iclu);
      /*
      switch(RawClusters[iclu]->View())
	{
	case geo::kU :
	   Cls[0].push_back(iclu);
	   break;
	case geo::kV :
	   Cls[1].push_back(iclu);
	   break;
	case geo::kZ :
	   Cls[2].push_back(iclu);
	   break;
	default :
	  break;
	}  
	*/
   // #############################################################
   // ### Filling wires and times into a TGraph for the cluster ###
   // #############################################################	  
   std::vector<double> wires;
   std::vector<double> times;
   
   // ### Counting the number of hits in the current cluster (n) ###
   int n = 0;
   // ############################################
   // ### Looping over the hits in the cluster ###
   // ############################################
   for(size_t i = 0; i < hit.size(); ++i)
      {
      wires.push_back(hit[i]->WireID().Wire);
      times.push_back(hit[i]->PeakTime());
      ++n;
      }//<---End loop over hits (i)
   // ################################################################
   // ### If there are 2 or more hits in the cluster fill a TGraph ###
   // ###         and fit a from a polynomial or order 1           ###
   // ################################################################
   if(n>=15)
      {
      // ###################################
      // ### Push the hits into a TGraph ###
      // ###################################
      TGraph *the2Dtrack = new TGraph(n,&wires[0],&times[0]);  
      // === Try to fit the TGraph with a 1st order polynomial ===
      try
	{
	the2Dtrack->Fit("pol1","Q");
	TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
	double par[2];
	double parerror[2];
	pol1->GetParameters(par);
	parerror[1] = pol1->GetParError(1);
	double FitChi2 = pol1->GetChisquare();
	double FitNDF = pol1->GetNDF();
			
	double fitGoodness = FitChi2/FitNDF;
			
	// ######################################################
	// ### Skipping the fitted slope if the Chi2/NDF < 5  ###
	// ######################################################
	if( fitGoodness > 10)
	   {
	   dtdwstart.push_back(std::tan(RawClusters[iclu]->StartAngle()));
	   continue;
				
	   }//<---End check on chi2/ndf fit

	// #######################################################################
	// ### Take change in time tick vs change in wire (dT/dW) from the fit ###
	// #######################################################################
	dtdwstart.push_back(par[1]);
	dtdwstartError.push_back(parerror[1]);
	}//<---End try to fit with a polynomial order 1
			

	// ### If the fitter fails just take dT/dW from the cluster ###
	catch(...)
	   {
	   mf::LogWarning("FeatureVertexFinder") << "Fitter failed, using the clusters default dTdW()";
	   delete the2Dtrack;
	   dtdwstart.push_back(std::tan(RawClusters[iclu]->StartAngle()));
	   continue;
	   }
	   
      delete the2Dtrack;
      }//<---End if the cluster has 2 or more hits
   // #################################################
   // ### If the cluster has fewer than 2 hits just ### 
   // ###      take the dT/dW from the cluster      ###
   // #################################################
   else {dtdwstart.push_back(std::tan(RawClusters[iclu]->StartAngle()));}
   }//<---End loop over clusters iclu
   
   
   
   // ##########################################################################################
   // ##########################################################################################
   // ### Now that I have slopes for all the clusters move on to finding intersection points ###
   // ##########################################################################################
   // ##########################################################################################
   
   
   // ##############################
   // ### Looping over cryostats ###
   // ##############################
   for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
      {
      // ##########################
      // ### Looping over TPC's ###
      // ##########################
      for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
	{
    	// #################################
	// ### Loop over the wire planes ###
	// #################################
	for (unsigned int i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
	   {
	   // ##############################################
	   // ### If there is at least one cluster found ###
	   // ##############################################
	   if (Cls[i].size() >= 1)
	      {
	      // ##############################
	      // ### Loop over each cluster ###
	      // ##############################
	      for (unsigned int j = 0; j<Cls[i].size(); ++j)
	      	{
		// === Current Clusters Plane ===
	        Clu_Plane.push_back(RawClusters.at(Cls.at(i).at(j))->View());
		// === Current Clusters StartPos ===
	        Clu_StartPos_Wire.push_back(RawClusters.at(Cls.at(i).at(j))->StartWire());
		Clu_StartPos_TimeTick.push_back(RawClusters.at(Cls.at(i).at(j))->StartTick());
	        // === Current Clusters EndPos ===
		Clu_EndPos_Wire.push_back(RawClusters.at(Cls.at(i).at(j))->EndWire());
	        Clu_EndPos_TimeTick.push_back(RawClusters.at(Cls.at(i).at(j))->EndTick());
	        // === Current Clusters Slope (In Wire and Time Tick)
		Clu_Slope.push_back(dtdwstart[Cls[i][j]]);
		Clu_Length.push_back(std::sqrt(pow((RawClusters.at(Cls.at(i).at(j))->StartWire()-RawClusters.at(Cls.at(i).at(j))->EndWire())*13.5,2) + 
		                              pow(RawClusters.at(Cls.at(i).at(j))->StartTick()-RawClusters.at(Cls.at(i).at(j))->EndTick(),2)));
		// ######################################################
		// ### Given a slope and a point find the y-intercept ###
		// ###                   c = y-mx                     ###
		// ######################################################
		Clu_Yintercept.push_back(RawClusters.at(Cls.at(i).at(j))->StartTick() - (dtdwstart[Cls[i][j]] * RawClusters.at(Cls.at(i).at(j))->StartWire()));
		// #################################################################
		// ###     Also calculating the y-intercept but using the        ###
		// ###   end time of the cluster correct for the possibility     ###
		// ### that the clustering didn't get start and end points right ###
		// #################################################################
		Clu_Yintercept2.push_back(RawClusters.at(Cls.at(i).at(j))->EndTick() - (dtdwstart[Cls[i][j]] * RawClusters.at(Cls.at(i).at(j))->EndWire()));
		
		// #######################################################
		// ### Iterating on the total number of clusters found ###
		// #######################################################
		
		nClustersFound++;
	      	}//<---End loop over all clusters
		
	      }//<---End check if we have at least one cluster
   	   // ################################################################
	   // ## If no clusters were found then put in dummy vertex values ###
	   // ################################################################
	   else
	      {
	      TwoDvtx_wire.push_back(-1);
	      TwoDvtx_time.push_back(-1);
	      TwoDvtx_plane.push_back(-1);
	      }//<---End no clusters found else statement
	   
	   }//<---End loop over planes (i)
	}//<---End loop over tpc's
      }//<---End loop over cryostats
   
   // ##########################################################################
   // ###  Now loop over all the clusters found and establish a preliminary  ###
   // ### set of 2d-verticies based on the slope/intercept of those clusters ###
   // ##########################################################################
   for(unsigned int n = nClustersFound; n > 0; n--)
      {
      // #######################################################
      // ###   Looping over the clusters starting from the   ###
      // ### first cluster and checking against the nCluster ###
      // #######################################################
      for (unsigned int m = 0; m < n; m++)
      	{
      	// ###########################################################
	// ### Checking to make sure clusters are in the same view ###
	// ###########################################################
	if(Clu_Plane[n] == Clu_Plane[m])
	   {
	   // --- Skip the vertex if the lines slope don't intercept ---
	   if(Clu_Slope[m] - Clu_Slope[n] == 0){break;}
	   // ============================================================
	   // === X intersection = (yInt2 - yInt1) / (slope1 - slope2) ===
	   float intersection_X = (Clu_Yintercept[n] - Clu_Yintercept[m]) / (Clu_Slope[m] - Clu_Slope[n]);
	   // ================================================
	   // === Y intersection = (slope1 * XInt) + yInt1 ===
	   float intersection_Y = (Clu_Slope[m] * intersection_X) + Clu_Yintercept[m];
	   // ============================================================
	   // === X intersection = (yInt2 - yInt1) / (slope1 - slope2) ===
	   float intersection_X2 = (Clu_Yintercept2[n] - Clu_Yintercept2[m]) / (Clu_Slope[m] - Clu_Slope[n]);
	   // ================================================
	   // === Y intersection = (slope1 * XInt) + yInt1 ===
	   float intersection_Y2 = (Clu_Slope[m] * intersection_X2) + Clu_Yintercept2[m];
	   
	   // #########################################
	   // ### Skipping crap intersection points ###
	   // #########################################
	   if(intersection_X2 < 1){intersection_X2 = -999;}
	   if(intersection_X2 > geom->Nwires(Clu_Plane[m],0,0)){intersection_X2 = -999;}
	   if(intersection_Y2 < 0){intersection_Y2 = -999;}
	   if(intersection_Y2 > detprop->NumberTimeSamples() ){intersection_Y2 = -999;}			
	   if(intersection_X < 1){intersection_X = -999;}
	   if(intersection_X > geom->Nwires(Clu_Plane[m],0,0)){intersection_X = -999;}
	   if(intersection_Y < 0){intersection_Y = -999;}
	   if(intersection_Y > detprop->NumberTimeSamples() ){intersection_Y = -999;}
	   
	   // ############################################################
	   // ### Putting in a protection for the findManyHit function ###
	   // ############################################################
	   try
	      {
	      // ### Gathering the hits associated with the current cluster ###
	      std::vector< art::Ptr<recob::Hit> > hitClu1 = fmhit.at(m);
	      std::vector< art::Ptr<recob::Hit> > hitClu2 = fmhit.at(n);
	   
	      // ### If the intersection point is 80 or more wires away from either cluster
	      // ### and one of the clusters has fewer than 8 hits the intersection
	      // ### is likely a crap one and we won't save this point
	      if((abs( Clu_EndPos_Wire[m] - intersection_X2 ) > 80 && hitClu1.size() < 8) ||
	         (abs( Clu_EndPos_Wire[n] - intersection_X2 ) > 80 && hitClu2.size() < 8) )
	         {intersection_X2 = -999;intersection_Y2 = -999;}
	      
	      if((abs( Clu_StartPos_Wire[m] - intersection_X ) > 80 && hitClu1.size() < 8) ||
	         (abs( Clu_StartPos_Wire[n] - intersection_X ) > 80 && hitClu2.size() < 8) )
	         {intersection_X = -999;intersection_Y = -999;}
	      
	      // ### If the intersection point is 50 or more wires away from either cluster
	      // ### and the one of the clusters has fewer than 3 hits the intersection
	      // ### is likely a crap one and we won't save this point 
	      if((abs( Clu_EndPos_Wire[m] - intersection_X2 ) > 50 && hitClu1.size() < 4) ||
	         (abs( Clu_EndPos_Wire[n] - intersection_X2 ) > 50 && hitClu2.size() < 4) )
	         {intersection_X2 = -999;intersection_Y2 = -999;}
						    
	      if((abs( Clu_StartPos_Wire[m] - intersection_X ) > 50 && hitClu1.size() < 4) ||
	         (abs( Clu_StartPos_Wire[n] - intersection_X ) > 50 && hitClu2.size() < 4) )
	         {intersection_X = -999;intersection_Y = -999;}
	      
	      }
	   catch(...) 
	      {
	      mf::LogWarning("FeatureVertexFinder") << "FindManyHit Function faild";
	      intersection_X = -999;intersection_Y = -999;
	      intersection_X2 = -999;intersection_Y2 = -999;
	      continue;
	      }
	   
	   // ##########################################################################
	   // ### Push back a candidate 2dClusterVertex if it is inside the detector ###
	   // ##########################################################################
	   if( intersection_X2 > 1 && intersection_Y2 > 0 && 
	     ( intersection_X2 < geom->Nwires(Clu_Plane[m],0,0) ) &&
	     ( intersection_Y2 < detprop->NumberTimeSamples() ) )
	      {
	      
	      TwoDvtx_wire.push_back(intersection_X2);
	      TwoDvtx_time.push_back(intersection_Y2);
	      TwoDvtx_plane.push_back(Clu_Plane[m]);				
	      }//<---End saving a "good 2d vertex" candidate
	      
	   // ##########################################################################
	   // ### Push back a candidate 2dClusterVertex if it is inside the detector ###
	   // ##########################################################################
	   if( intersection_X > 1 && intersection_Y > 0 && 
	     ( intersection_X < geom->Nwires(Clu_Plane[m],0,0) ) &&
	     ( intersection_Y < detprop->NumberTimeSamples() ) )
	      {
	      TwoDvtx_wire.push_back(intersection_X);
	      TwoDvtx_time.push_back(intersection_Y);
	      TwoDvtx_plane.push_back(Clu_Plane[m]);				
	      }//<---End saving a "good 2d vertex" candidate

	   }//<---End check that they are in differnt planes
      	}//<---End m loop
      }//<---End n loop

   }//<---End 2dClusterVertexCandidates












// -----------------------------------------------------------------------------   
// Get 3d Vertex Candidates from clusters 2d Vertex candidates
// ----------------------------------------------------------------------------- 
void vertex::FeatureVertexFinder::Find3dVtxFrom2dClusterVtxCand(std::vector<double> Wire_2dvtx, std::vector<double> Time_2dvtx, std::vector<double> Plane_2dvtx)
   {
   
   //std::cout<<"####################################################################"<<std::endl;
   //std::cout<<"Get 3d Vertex Candidates from clusters 2d Vertex candidates Function"<<std::endl;
   //std::cout<<std::endl;
   std::vector<double> vtx_wire_merged;
   std::vector<double> vtx_time_merged;
   std::vector<double> vtx_plane_merged;
   
   double y_coord = 0, z_coord = 0;
   
   bool merged = false;

   // #########################
   // ### Geometry Services ###
   // #########################
   art::ServiceHandle<geo::Geometry const> geom;
   
   // ####################################
   // ### Detector Properties Services ###
   // ####################################
   const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   
   // ---------------------- MERGING THE LONG LIST OF 2D CANDIDATES ---------------------------
   
   // #########################################
   // ### Looping over 2d-verticies (loop1) ###
   // #########################################
   for(size_t vtxloop1 = 0 ; vtxloop1 < Wire_2dvtx.size(); vtxloop1++)
      {
      if(Wire_2dvtx[vtxloop1] < 0){continue;}
	
      merged = false;
      // #########################################
      // ### Looping over 2d-verticies (loop2) ###
      // #########################################
      for(size_t vtxloop2 = vtxloop1+1; vtxloop2 < Wire_2dvtx.size(); vtxloop2++)
	{
	if(Wire_2dvtx[vtxloop2] < 0){continue;} 
		
	// ########################################################
	// ### Make sure the 2d-Verticies are in the same plane ###
	// ########################################################
	if(Plane_2dvtx[vtxloop1] == Plane_2dvtx[vtxloop2])
	   {
	   // ###############################################			
	   // ### Considering merging 2d vertices if they ###
	   // ###    are within 3 wires of each other     ###
	   // ###############################################
	   if( fabs(Wire_2dvtx[vtxloop1] - Wire_2dvtx[vtxloop2]) < 4 )
	      {
	      // ############################################################
	      // ### Merge the verticies if they are within 10 time ticks ###
	      // ############################################################
	      if( fabs(Time_2dvtx[vtxloop1] - Time_2dvtx[vtxloop2]) < 10 )
	      	{
		vtx_wire_merged.push_back( ((Wire_2dvtx[vtxloop2] + Wire_2dvtx[vtxloop1])/ 2) );
		vtx_time_merged.push_back( ((Time_2dvtx[vtxloop2] + Time_2dvtx[vtxloop1])/ 2) ) ;
		vtx_plane_merged.push_back( Plane_2dvtx[vtxloop1] );
					
		merged = true;		
		if(vtxloop2<Wire_2dvtx.size()){vtxloop2++;}
		if(vtxloop1<Wire_2dvtx.size()){vtxloop1++;}				
		}//<---End the check within 10 time ticks for merging
	      }//<---Looking at vertices that are within 3 wires of each other
	   }//<---Only looking at vertices that are in the same plane	
	}//<---End vtxloop2
	if(!merged)
		{
		vtx_wire_merged.push_back( Wire_2dvtx[vtxloop1] );
		vtx_time_merged.push_back( Time_2dvtx[vtxloop1] );
		vtx_plane_merged.push_back( Plane_2dvtx[vtxloop1] );
		}//<---end saving unmerged verticies
      }//<---End vtxloop1
   
   // #####################################
   // ### Variables for channel numbers ###
   // #####################################
   uint32_t     vtx1_channel = 0;
   uint32_t     vtx2_channel = 0;
   // --------------------------------------------------------------------------
   // ---   Having now found a very long list of potential 2-d end points    ---
   // --- we need to check if any of them match between planes and only keep ---
   // ---                       those that have matches                      ---
   // --------------------------------------------------------------------------
   // ##############################
   // ### Looping over cryostats ###
   // ##############################
   for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
      {
       // ##########################
       // ### Looping over TPC's ###
       // ##########################
       for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
      	{
    	for(unsigned int vtx = vtx_wire_merged.size(); vtx > 0; vtx--)
    	   {
	   for (unsigned int vtx1 = 0; vtx1 < vtx; vtx1++)
	      {
	      // ###########################################################################
	      // ### Check to make sure we are comparing verticies from different planes ###
	      // ###########################################################################
	      if(vtx_plane_merged[vtx1] != vtx_plane_merged[vtx])
	      	{
		// === To figure out if these two verticies are from a common point 
		// === we need to check if the channels intersect and if they are 
		// === close in time ticks as well...to do this we have to do some 
		// === converting to use geom->PlaneWireToChannel(PlaneNo, Wire, tpc, cstat)
		bool match = false;
			
		unsigned int vtx1_plane   = vtx_plane_merged[vtx1];
		unsigned int vtx1_wire    = vtx_wire_merged[vtx1];
		try
		   {vtx1_channel = geom->PlaneWireToChannel(vtx1_plane, vtx1_wire, tpc, cstat);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "PlaneWireToChannel Failed";
		   match = false;
	            continue;}
		
		unsigned int vtx2_plane   = vtx_plane_merged[vtx];
		unsigned int vtx2_wire    = vtx_wire_merged[vtx];
		try
		   {vtx2_channel = geom->PlaneWireToChannel(vtx2_plane, vtx2_wire, tpc, cstat);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "PlaneWireToChannel Failed";
		   match = false;
	            continue;}
			
		// ##############################################################################
		// ### Check to see if the channels intersect and save the y and z coordinate ###
		// ##############################################################################
		
		try
		   {match = geom->ChannelsIntersect( vtx1_channel, vtx2_channel, y_coord, z_coord);}
		catch(...)
		   {mf::LogWarning("FeatureVertexFinder") << "match failed for some reason";
	            match = false;
	      	    continue;}
		// ####################################################################
		// ### If the channels intersect establish if they are close in "X" ###
		// ####################################################################
		if( match )
		   {
		   float tempXCluster1 = detprop->ConvertTicksToX(vtx_time_merged[vtx1], vtx1_plane, tpc, cstat);
		   float tempXCluster2 = detprop->ConvertTicksToX(vtx_time_merged[vtx], vtx2_plane, tpc, cstat);
		   // ###############################################################################
		   // ### Now check if the matched channels are within 0.5 cm when projected in X ###
		   // ###         and that we have less than 100 of these candidates...           ###
		   // ###                 because more than that seems silly                      ###
		   // ###############################################################################
		   if(std::abs(tempXCluster1 - tempXCluster2) < 0.5 && candidate_x.size() < 101)
		      {
		      //        detprop->ConvertTicksToX(ticks, plane, tpc, cryostat)
		      candidate_x.push_back( detprop->ConvertTicksToX(vtx_time_merged[vtx1], vtx_plane_merged[vtx1], tpc, cstat) );
		      candidate_y.push_back( y_coord );
		      candidate_z.push_back( z_coord );
		      candidate_strength.push_back( 10 ); //<--For cluster verticies I give it a strength of "10" arbitrarily for now
					
		      }//<---End Checking if the vertices agree "well enough" in time tick
		   }//<---End Checking if verticies intersect

		}//<--- End checking we are in different planes
	      }//<---end vtx1 for loop
	   }//<---End vtx for loop
      	}//<---End loop over TPC's 
      }//<---End loop over cryostats

   }//<---End Find3dVtxFrom2dClusterVtxCand









// -----------------------------------------------------------------------------   
// Get 3d Vertex Candidates from clusters 2d Vertex candidates
// ----------------------------------------------------------------------------- 
void vertex::FeatureVertexFinder::MergeAndSort3dVtxCandidate(std::vector<double> merge_vtxX, std::vector<double> merge_vtxY, std::vector<double> merge_vtxZ, std::vector<double> merge_vtxStgth)
   {
   
   std::vector<double> x_3dVertex_dupRemoved = {0.};
   std::vector<double> y_3dVertex_dupRemoved = {0.};
   std::vector<double> z_3dVertex_dupRemoved = {0.};
   std::vector<double> strength_dupRemoved   = {0.};
   
   // #####################################################
   // ### Looping over the 3d candidates found thus far ###
   // #####################################################
   for(size_t dup = 0; dup < merge_vtxX.size(); dup ++)
      {
      // ### Temperary storing the current vertex ###
      float tempX_dup = merge_vtxX[dup];
      float tempY_dup = merge_vtxY[dup];
      float tempZ_dup = merge_vtxZ[dup];
      float tempStgth = merge_vtxStgth[dup];
      
      // #######################################################
      // ### Setting a boolian to see if this is a duplicate ###
      // #######################################################
      bool duplicate_found = false;
      // #############################################################
      // ### Looping over the rest of the list for duplicate check ###
      // #############################################################
      for(size_t check = dup+1; check < merge_vtxX.size(); check++)
	{
	// #############################################################################
	// ### I am going to call a duplicate vertex one that matches in x, y, and z ###
	// ###        within 0.1 cm for all 3 coordinates simultaneously             ###
	// #############################################################################
	if(std::abs( merge_vtxX[check] - tempX_dup ) < 0.1 && std::abs( merge_vtxY[check] - tempY_dup ) < 0.1 &&
	   std::abs( merge_vtxZ[check] - tempZ_dup ) < 0.1 )
	   {duplicate_found = true;}//<---End checking to see if this is a duplicate vertex

	}//<---End check for loop
	
	// ######################################################################
	// ### If we didn't find a duplicate then lets save this 3d vertex as ###
	// ###            a real candidate for consideration                  ###
	// ######################################################################
	if(!duplicate_found && tempX_dup > 0)
	   {
	   x_3dVertex_dupRemoved.push_back(tempX_dup);
	   y_3dVertex_dupRemoved.push_back(tempY_dup);
	   z_3dVertex_dupRemoved.push_back(tempZ_dup);
	   strength_dupRemoved.push_back(tempStgth);
	   }//<---End storing only non-duplicates	

	}//<---End dup for loop
   
   // ########################################################################################
   // ### Sorting the verticies I have found such that the first in the list is the vertex ###
   // ###         with the highest vertex strength and the lowest z location               ###
   // ########################################################################################

   int flag = 1;
   double tempX, tempY, tempZ, tempS;
   
   // #####################################################
   // ### Looping over all duplicate removed candidates ###
   // #####################################################
   for(size_t npri = 0; (npri < x_3dVertex_dupRemoved.size()) && flag; npri++)
      {
      flag = 0;
      for(size_t mpri = 0; mpri < x_3dVertex_dupRemoved.size() -1; mpri++)
	{
	// Swap the order of the two elements
	if(strength_dupRemoved[mpri+1] > strength_dupRemoved[mpri] || 
	  (strength_dupRemoved[mpri+1] == strength_dupRemoved[mpri] && z_3dVertex_dupRemoved[mpri] > z_3dVertex_dupRemoved[mpri+1]) )
	   {
	   tempX = x_3dVertex_dupRemoved[mpri];
	   x_3dVertex_dupRemoved[mpri] = x_3dVertex_dupRemoved[mpri+1];
	   x_3dVertex_dupRemoved[mpri+1] = tempX;
			
	   tempY = y_3dVertex_dupRemoved[mpri];
	   y_3dVertex_dupRemoved[mpri] = y_3dVertex_dupRemoved[mpri+1];
	   y_3dVertex_dupRemoved[mpri+1] = tempY;
			
	   tempZ = z_3dVertex_dupRemoved[mpri];
	   z_3dVertex_dupRemoved[mpri] = z_3dVertex_dupRemoved[mpri+1];
	   z_3dVertex_dupRemoved[mpri+1] = tempZ;
			
	   tempS = strength_dupRemoved[mpri];
	   strength_dupRemoved[mpri] = strength_dupRemoved[mpri+1];
	   strength_dupRemoved[mpri+1] = tempS;
			
	   flag = 1;
			
	   }//<---Inside swap loop
	}//<---End mpri	
      }//<---End npri loop
      
   // ############################################################
   // ### Pushing into a vector of merged and sorted verticies ###
   // ############################################################
   for(size_t count = 0; count < x_3dVertex_dupRemoved.size(); count++)
      {
      MergeSort3dVtx_xpos.push_back(x_3dVertex_dupRemoved[count]);
      MergeSort3dVtx_ypos.push_back(y_3dVertex_dupRemoved[count]);
      MergeSort3dVtx_zpos.push_back(z_3dVertex_dupRemoved[count]);
      MergeSort3dVtx_strength.push_back(strength_dupRemoved[count]);
      
      
      
      }//<---End count loop    
   
   
   }// End MergeAndSort3dVtxCandidate

  DEFINE_ART_MODULE(FeatureVertexFinder)   
   
// -----------------------------------------------------------------------------

}//<---End namespace vertex
