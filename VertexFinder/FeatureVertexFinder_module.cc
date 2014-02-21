////////////////////////////////////////////////////////////////////////
//
// FeatureVertexFinder class
//
// jasaadi@fnal.gov
//
// This algorithm is designed to reconstruct the vertices using the
// 2D cluster information and CornerFinderAlg to find 2-d and 3-d verticies
//
// We will utilize many of the methods found in VeretxFinder2d but modify them 
// to return vertex candidates and only a point that has matching in at least 
// 2 views will be used
//
//
//
// This is Preliminary Work and needs modifications
//
// ////////////////////////////////////////////////////////////////////////

// ##########################
// ### Basic C++ Includes ###
// ##########################
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Vertex.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/ClusterParamsAlg.h"

// #####################
// ### ROOT Includes ###
// #####################
#include "TMath.h"
#include "TH1D.h"
#include "TVectorD.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TVector3.h"

// =====================================================================================================
// =====================================================================================================
#ifndef FeatureVertexFinder_H
#define FeatureVertexFinder_H

class TH1D;

///vertex reconstruction
namespace vertex {
   
 class FeatureVertexFinder :  public art::EDProducer {
    
  public:
    
    explicit FeatureVertexFinder(fhicl::ParameterSet const& pset); 
    virtual ~FeatureVertexFinder();        
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);

  private:
    
    TH1D *dtIC;
  
   std::string fClusterModuleLabel;
   std::string fHitModuleLabel;
   std::string fCornerFinderModuleLabel;
   cluster::ClusterParamsAlg fClParAlg;
  };
    
}//<---End namespace vertex

#endif // FeatureVertexFinder_H
// =====================================================================================================
// =====================================================================================================








// =====================================================================================================
// =====================================================================================================
namespace vertex{

//-----------------------------------------------------------------------------
  FeatureVertexFinder::FeatureVertexFinder(fhicl::ParameterSet const& pset):
  fClParAlg(pset.get<fhicl::ParameterSet>("ClusterParamsAlg"), pset.get< std::string >("module_type"))
  {  
    /*this->*/reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::EndPoint2D, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
  }
//-----------------------------------------------------------------------------
  FeatureVertexFinder::~FeatureVertexFinder()
  {
  }

  //---------------------------------------------------------------------------
  void FeatureVertexFinder::reconfigure(fhicl::ParameterSet const& p) 
  {
    fCornerFinderModuleLabel  = p.get< std::string >("CornerFinderModuleLabel");
    fClusterModuleLabel       = p.get< std::string >("ClusterModuleLabel");
    fHitModuleLabel	      = p.get< std::string >("HitModuleLabel");

    return;
  }
  //-------------------------------------------------------------------------
  void FeatureVertexFinder::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    dtIC = tfs->make<TH1D>("dtIC","It0-Ct0",100,-5,5);
    dtIC->Sumw2();
  }

// //-----------------------------------------------------------------------------
  void FeatureVertexFinder::produce(art::Event& evt)
  {
  
  
    //std::cout<<"Top of FeatureVertexFinder Produce"<<std::endl;  
    // #########################
    // ### Geometry Services ###
    // #########################
    art::ServiceHandle<geo::Geometry> geom;
    
    // ###############################
    // ### LAr Properties Services ###
    // ###############################
    art::ServiceHandle<util::LArProperties> larprop;
    
    // ####################################
    // ### Detector Properties Services ###
    // ####################################
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    // ########################################
    // ### cout Run Number and Event Number ###
    // ########################################
    /*std::cout<< std::endl;
    std::cout << " ======================================================" << std::endl;
    std::cout << "Run    : " << evt.run() <<" Event  : "<<evt.id().event() << std::endl;
    std::cout << " ======================================================" << std::endl;
    std::cout << std::endl;*/

    
    
    // ###############################
    // ### Defining TPC Parameters ###
    // ###############################
    // === TPC Name === (JA: Not sure if I need this)
    TString tpcName = geom->GetLArTPCVolumeName();
    // === Length of the detector ===
    double DetectorLength =  (geom->DetHalfHeight())*2.;
    // === wire angle with respect to the vertical direction ===
    double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; 
    
    // ==================================================
    // === Calculating the Timetick to CM conversion  ===
    // ==================================================
    double TimeTick = detprop->SamplingRate()/1000.; //<---To get units of microsecond...not nanosec
    double Efield_drift = larprop->Efield();      // Electric Field in the drift region in kV/cm
    double Temperature  = larprop->Temperature(); // LAr Temperature in K
    // === Drift velocity in the drift region (cm/us) ===
    double DriftVelocity = larprop->DriftVelocity(Efield_drift,Temperature); 
    // === Time sample (cm) === 
    double TimetoCm = DriftVelocity*TimeTick; 
    
    
    double presamplings = detprop->TriggerOffset(); //trigger offset
    
    std::vector<double> vtx_wire = {0.};
    std::vector<double> vtx_time = {0.};
    std::vector<double> vtx_plane = {0.};
    
    
    
 //-------------------------------------------------------------------------------------------------------------------------------------------------    
    
    //std::cout<<"Making new collections for output"<<std::endl;
    //Point to a collection of vertices to output.
    std::unique_ptr<std::vector<recob::Vertex> >                 vcol(new std::vector<recob::Vertex>);          //3D vertex
    std::unique_ptr<std::vector<recob::EndPoint2D> >             epcol(new std::vector<recob::EndPoint2D>);  //2D vertex
    std::unique_ptr< art::Assns<recob::EndPoint2D, recob::Hit> > assnep(new art::Assns<recob::EndPoint2D, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> >  assnsh(new art::Assns<recob::Vertex, recob::Shower>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track> >   assntr(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit> >     assnh(new art::Assns<recob::Vertex, recob::Hit>);
    

    // ###################################################
    // ### Retreiving the Cluster Module for the event ###
    // ###################################################
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);
    
    // #################################################
    // ### Finding hits associated with the clusters ###
    // #################################################
    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);
    
   
    
    // ##################################
    // ### Filling the Cluster Vector ###
    // ##################################
    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
      
    }//<---End ii loop
    
    // ###############################################
    // ### Retreiving the Hit Module for the event ###
    // ###############################################
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    evt.getByLabel(fHitModuleLabel,hitListHandle);
    
    // #########################################
    // ### Putting Hits into a vector (hits) ###
    // #########################################
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitListHandle);   
    
//-------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------ CORNERFINDER ALG -------------------------------------------------------------------  
//-------------------------------------------------------------------------------------------------------------------------------------------------
 
    // ####################################################################
    // ### Utilizing the CornerFinderAlg to get a list of EndPoint2d's  ###
    // ####################################################################
    // Change from Wes (2/21/14): Take CornerFinderModule as input to algorithm here.
    
    std::vector<unsigned int>   feature_wire = {0};
    std::vector<double>   feature_time = {0.};
    std::vector<unsigned int>   feature_plane = {0};
    std::vector<double>	  feature_strength = {0.};
    std::vector<uint32_t> feature_channel =  {0};
    int nFeatures = 0;
    
    float x_feature[1000] = {0.}, y_feature[1000] = {0.}, z_feature[1000] = {0.}, strength_feature[1000] = {0.};
    double y = 0., z = 0.;
    double yy = 0., zz = 0.;
    double yy2 = 0., zz2 = 0.;
    double yy3 = 0., zz3 = 0.;
    int n3dFeatures = 0;
    
    
    art::Handle< std::vector<recob::EndPoint2D> > CornerFinderHandle;
    evt.getByLabel(fCornerFinderModuleLabel,CornerFinderHandle);
    std::vector< art::Ptr<recob::EndPoint2D> > EndPoints;
    art::fill_ptr_vector(EndPoints, CornerFinderHandle);   
        
    // ########################################################
    // ### Loop over all the features and record their info ###
    // ########################################################
    for(size_t i=0; i!=EndPoints.size(); ++i)
    	{
	
	feature_wire.push_back(EndPoints.at(i)->WireID().Wire);
	feature_time.push_back(EndPoints.at(i)->DriftTime());
	feature_plane.push_back(EndPoints.at(i)->WireID().Plane);
	feature_channel.push_back(geom->PlaneWireToChannel(EndPoints.at(i)->WireID().Plane, EndPoints.at(i)->WireID().Wire, 0, 0));
	feature_strength.push_back(EndPoints.at(i)->Strength());
	nFeatures++;
	
	
	
       	}//<---End i loop finding 2d Features
    
    
    bool GT2PlaneDetector = false;
    
    
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
    		// ##############################################################################
    		// ### Now loop over those features and see if they match up in between views ###
    		// ##############################################################################
    		for(int feature1 = 0; feature1 < nFeatures; feature1 ++)
    			{
			for(int feature2 = feature1+1; feature2 < nFeatures; feature2++)
				{
				// ##########################################################################
				// ### Check to make sure we are comparing features from different planes ###
				// ##########################################################################
				if(feature_plane[feature1] != feature_plane[feature2])
					{
					
					// #############################################################################
					// ### Get the appropriate time offset for the two planes we are considering ###
					// #############################################################################
					float tempXFeature1 = detprop->ConvertTicksToX(feature_time[feature1], feature_plane[feature1], tpc, cstat);
					
					float tempXFeature2 = detprop->ConvertTicksToX(feature_time[feature2], feature_plane[feature2], tpc, cstat);
					
					
					// ####################################################################
					// ### Checking to see if these features have intersecting channels ###
					// ###               and are within 1 cm in projected X             ###
					// ####################################################################
					if( geom->ChannelsIntersect( feature_channel[feature2], feature_channel[feature1], yy, zz) &&
					    std::abs(tempXFeature1 - tempXFeature2) < 0.5)
					    	{
						// #####################################################################################
						// ### Use this fill if we are in a detector with fewer than 3 plane (e.g. ArgoNeuT) ###
						// #####################################################################################
						if(!GT2PlaneDetector)
							{
							x_feature[n3dFeatures] = detprop->ConvertTicksToX(feature_time[feature1], feature_plane[feature1], tpc, cstat); 
							y_feature[n3dFeatures] = yy;
							z_feature[n3dFeatures] = zz;
							strength_feature[n3dFeatures] = feature_strength[feature1] + feature_strength[feature1];
				
							n3dFeatures++;
							}//<---End fill for 2 plane detector
						
						
						// ##############################################################################
						// ### Adding a check to see if I am in a 3-plane detector and therefore need ###
						// ###           to check for a match across more than 2 planes               ###
						// ##############################################################################
						if(GT2PlaneDetector)
							{
							
							for(int feature3 = feature2+1; feature3 < nFeatures; feature3++)
								{
								// ##########################################################################
								// ### Check to make sure we are comparing features from different planes ###
								// ##########################################################################
								if(feature_plane[feature3] != feature_plane[feature2] && feature_plane[feature3] != feature_plane[feature1] && 
								   feature_plane[feature1] != feature_plane[feature2] )
									{
									float tempXFeature3 = detprop->ConvertTicksToX(feature_time[feature3], feature_plane[feature3], tpc, cstat);
									
									// ####################################################################################
									// ### Checking to make sure our third feature has an intersecting channel with our ###
									// ###         other two channels and is within 1 cm projected in X                 ###
									// ####################################################################################
									if( geom->ChannelsIntersect( feature_channel[feature3], feature_channel[feature1], yy3, zz3) &&
									    geom->ChannelsIntersect( feature_channel[feature3], feature_channel[feature2], yy2, zz2) &&
									    geom->ChannelsIntersect( feature_channel[feature2], feature_channel[feature1], yy, zz) &&
					    				     std::abs(tempXFeature3 - tempXFeature2) < 0.5 && std::abs(tempXFeature3 - tempXFeature1) < 0.5 &&
									     std::abs(tempXFeature1 - tempXFeature2) < 0.5 )
										{
										x_feature[n3dFeatures] = detprop->ConvertTicksToX(feature_time[feature1], feature_plane[feature1], tpc, cstat); 
										geom->IntersectionPoint(feature_wire[feature1], feature_wire[feature2], feature_plane[feature1], feature_plane[feature2], cstat, tpc, y, z);
										//std::cout<<"y = "<<y<<" , z = "<<z<<std::endl;
										y_feature[n3dFeatures] = y;
										z_feature[n3dFeatures] = z;
										strength_feature[n3dFeatures] = feature_strength[feature1] + feature_strength[feature2] + feature_strength[feature3];

										n3dFeatures++;
										
										
										// ### Note: If I've made it here I have a matched triplet...since I don't want to use any of these ###
										// ###    features again I am going to iterate each of the counters so we move to the next one      ###
										if(feature1 < nFeatures)
											{feature1++;}
										if(feature2 < nFeatures)
											{feature2++;}
										if(feature3 < nFeatures)
											{feature3++;}
									
									
										}//<---End Finding a 3d point that matches in all three planes
									}//<----End checking that we have different planes
								}//<---End feature3 loop
							}//<---End Greater than 2 plane detector check
	
						}//<---End Checking if the feature matches in 3d
					}//<---End Checking features are in different planes
				}//<---End feature2 loop
			}//<---End feature1 loop
		}//<---End looping over TPC's
	}//<---End looping over cryostats
 	

 // === Skipping 3-d 	
 int nGood3dFeatures = 0;
 float temp_x[1000] = {0.}, temp_y[1000] = {0.}, temp_z[1000] = {0.}, temp_s[1000];
 
 for(int check = n3dFeatures; check > 0; check--)
 	{
	if(strength_feature[check] < 100){continue;}

	if(std::abs(x_feature[check] - x_feature[check-1]) < 0.1 && std::abs(y_feature[check] - y_feature[check-1]) < 0.1 && 
	   std::abs(z_feature[check] - z_feature[check-1]) < 0.1 )
	   	{continue;}
		
	temp_x[nGood3dFeatures] = x_feature[check];
	temp_y[nGood3dFeatures] = y_feature[check];
	temp_z[nGood3dFeatures] = z_feature[check];
	temp_s[nGood3dFeatures] = strength_feature[check];
	nGood3dFeatures++;

	}
 
 

// #########################################################################
// ### Zeroing the main counter since we are looping over a trimmed list ###
// #########################################################################
n3dFeatures = 0;

bool Plane0HitMatch = false;
bool Plane1HitMatch = false;
bool Plane2HitMatch = false;	
for(int checkz = nGood3dFeatures; checkz > 0; checkz--)
	{
	if(temp_x[checkz] == temp_x[checkz-1] && temp_y[checkz] == temp_y[checkz-1] && 
	   temp_z[checkz] == temp_z[checkz-1] )
	   	{continue;}
	
	
	
	
	
	double temp_xyz[3] = {0.};
	temp_xyz[0] = temp_x[checkz];
	temp_xyz[1] = temp_y[checkz];
	temp_xyz[2] = temp_z[checkz];
	
	
	// =============================================================
	// === Looping over hits to see if any are near this feature ===
	for(size_t nFeatureHits = 0; nFeatureHits < hitListHandle->size(); nFeatureHits++)
		{
							
		// === Finding Wire/Plane associated with the hit ===
       		art::Ptr<recob::Hit> FeatureHit(hitListHandle, nFeatureHits);
		
		// ==========================
		// === Looking in Plane 0 ===
		// ==========================
		if(FeatureHit->WireID().Plane == 0)
			{
			//std::cout<<"Checking Plane 0"<<std::endl;
			// ===========================================================================================
			// === If there is a hit within 5 wires and 20 time ticks in this plane set boolian = true ===
			if(std::abs(FeatureHit->WireID().Wire - geom->NearestWire(temp_xyz , 0, 0, 0)) < 6 &&
			   std::abs(FeatureHit->PeakTime() - detprop->ConvertXToTicks(temp_xyz[0],0, 0, 0)) < 20)
			   	{
				Plane0HitMatch = true;}
			
			
			}//<---End Plane 0
		// ==========================
		// === Looking in Plane 1 ===
		// ==========================
		if(FeatureHit->WireID().Plane == 1)
			{
			//std::cout<<"Checking Plane 1"<<std::endl;
			// ===========================================================================================
			// === If there is a hit within 5 wires and 20 time ticks in this plane set boolian = true ===
			if(std::abs(FeatureHit->WireID().Wire - geom->NearestWire(temp_xyz , 1, 0, 0)) < 6 &&
			   std::abs(FeatureHit->PeakTime() - detprop->ConvertXToTicks(temp_xyz[0],1, 0, 0)) < 20)
			   	{
				Plane1HitMatch = true;}
			
			
			}//<---End Plane 1
			
		// ==========================
		// === Looking in Plane 2 ===
		// ==========================
		if(GT2PlaneDetector && FeatureHit->WireID().Plane == 2)
			{
			//std::cout<<"Checking Plane 2 "<<std::endl;
			// ===========================================================================================
			// === If there is a hit within 5 wires and 20 time ticks in this plane set boolian = true ===
			if(std::abs(FeatureHit->WireID().Wire - geom->NearestWire(temp_xyz , 2, 0, 0)) < 6 &&
			   std::abs(FeatureHit->PeakTime() - detprop->ConvertXToTicks(temp_xyz[0],2, 0, 0)) < 20)
			   	{
				Plane2HitMatch = true;}
			
			
			}//<---End Plane 1
	
		}//<---End nFeatureHit loop 
	
	// ==========================================================================================
	// === Requiring that the feature can be matched to a nearby hit in all possible views to === 
	// ===              be considered as a candidate 3d fearture vertex point                 ===
	// ==========================================================================================
	if(Plane0HitMatch && Plane1HitMatch && Plane2HitMatch)
		{
		/*std::cout<<"X = "<<temp_x[checkz] <<" , Y = "<<temp_y[checkz]<<" , Z = "<<temp_z[checkz]<<" strength = "<<temp_s[checkz]<<std::endl;
		std::cout<<std::endl;
		std::cout<<"Plane = 0 , Wire = "<<geom->NearestWire(temp_xyz , 0, 0, 0)<<" Time = "<<detprop->ConvertXToTicks(temp_xyz[0],0, 0, 0)<<std::endl;
		std::cout<<"Plane = 1 , Wire = "<<geom->NearestWire(temp_xyz , 1, 0, 0)<<" Time = "<<detprop->ConvertXToTicks(temp_xyz[0],1, 0, 0)<<std::endl;
		std::cout<<"Plane = 2 , Wire = "<<geom->NearestWire(temp_xyz , 2, 0, 0)<<" Time = "<<detprop->ConvertXToTicks(temp_xyz[0],2, 0, 0)<<std::endl;
		std::cout<<std::endl;*/
	
		x_feature[n3dFeatures] = temp_x[checkz];
		y_feature[n3dFeatures] = temp_y[checkz];
		z_feature[n3dFeatures] = temp_z[checkz];
		strength_feature[n3dFeatures] = temp_s[checkz];
		n3dFeatures++;
		
		}
	}


//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------- CLUSTER INTERCEPT INFORMATION -------------------------------------------------------- 
//-----------------------------------------------------------------------------------------------------------------------------------------------



    // ==============================================================================
    // ===  Loop over all the clusters that are found in the cluster list and fit ===
    // ==  the hits in those clusters to establish a vector of slopes (dtdwstart) ===
    // === that has the same number and order as our cluster PtrVector (clusters) ===
    // ==============================================================================
    std::vector<int> Cls[3]; //<---- Index to clusters in each view
    std::vector<double> dtdwstart; //<----Slope (delta Time Tick vs delta Wire) 
    std::vector<double> dtdwstartError; //<---Error on the slope    
    
    for(size_t iclu = 0; iclu < clusters.size(); ++iclu)
    	{
	// ######################################################
	// ### Determine which view the current cluster is in ###
	// ######################################################
	switch(clusters[iclu]->View())
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
	  
	// #############################################################
	// ### Filling wires and times into a TGraph for the cluster ###
	// #############################################################	  
	std::vector<double> wires;
	std::vector<double> times;
	  
	// ### Gathering the hits associated with the current cluster ###
	std::vector< art::Ptr<recob::Hit> > hit = fmh.at(iclu);
	
	// ### Counting the number of hits in the current cluster (n) ###
	int n = 0;
	// ############################################
	// ### Looping over the hits in the cluster ###
	// ############################################
	for(size_t i = 0; i < hit.size(); ++i)
		{
		// ### JA: Should I convert wire and time to cm? Will that help the fitting? ###
		
	    	wires.push_back(hit[i]->WireID().Wire);
	    	times.push_back(hit[i]->PeakTime());
	    	++n;
	  	}//<---End loop over hits (i)
	// ################################################################
	// ### If there are 2 or more hits in the cluster fill a TGraph ###
	// ###         and fit a from a polynomial or order 1           ###
	// ################################################################
	if(n>=2)
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
			
			// ############################################################
			// ### Skipping the fitted slope if the 0.5 < Chi2/NDF < 5  ###
			// ############################################################
			if( fitGoodness > 5)
				{
				mf::LogWarning("FeatureVertexFinder") << "Fitter returned poor Chi2/NDF, using the clusters default dTdW()";
				dtdwstart.push_back(clusters[iclu]->dTdW());

				continue;
				
				}//<---End check on chi2/ndf fit

	      		// #######################################################################
	      		// ### Take change in time tick vs change in wire (dT/dW) from the fit ###
	      		// #######################################################################
	      		dtdwstart.push_back(par[1]);
			dtdwstartError.push_back(parerror[1]);
	    		}//<---End try to fit with a polynomial order 1
			
	    	// ############################################################
	    	// ### If the fitter fails just take dT/dW from the cluster ###
	    	// ############################################################
	    	catch(...)
			{
	      		mf::LogWarning("FeatureVertexFinder") << "Fitter failed, using the clusters default dTdW()";
	      		delete the2Dtrack;
	      		dtdwstart.push_back(clusters[iclu]->dTdW());
	      		continue;
	    		}
	    	delete the2Dtrack;
	  	}//<---End if the cluster has 2 or more hits
	// #################################################
	// ### If the cluster has fewer than 2 hits just ### 
	// ###      take the dT/dW from the cluster      ###
	// #################################################
	else {dtdwstart.push_back(clusters[iclu]->dTdW());}
	}//<---End loop over clusters iclu
	
	
//-------------------------------------------------------------------------------------------------------------------------------------------------    
    
    
    
    
    // ##########################################
    // ### Variables to save for each cluster ###
    // ##########################################
    int   nClustersFound = 0;				//<---Number of clusters to be evaluated in a single plane 
    							//    (this gets zeroed after looping over all clusters in a plane)
							
    int   n2dVertexCandidates = 0;			//<---Number of candidate 2d Vertices found
    float Clu_Plane[10000] = {0};         	     	//<---Plane of the current cluster
    float Clu_StartPos_Wire[10000]= {0};     		//<---Starting wire number of cluster
    float Clu_StartPos_TimeTick[10000]= {0};      	//<---Starting TDC value of the cluster
  
    float Clu_EndPos_Wire[10000]= {0};       		//<---Ending wire number of cluster
    float Clu_EndPos_TimeTick[10000]= {0};        	//<---Ending TDC value of the cluster
  
    float Clu_Slope[10000]= {0};	  	   		//<---Calculated Slope of the cluster (TDC/Wire)
    float Clu_Yintercept[10000]= {0};			//<---Clusters Y Intercept using start positions
    float Clu_Yintercept2[10000]= {0};			//<---Clusters Y Intercept using end positions
    float Clu_Length[10000]= {0};			//<---Calculated Length of the cluster
    
    int   AllCluster = 0;
    float AllCluster_Plane[10000] = {0.};		//<---Storing the plane # for all clusters
    float AllCluster_Length[10000] = {0.};		//<---Storing the length for all clusters
    float AllCluster_StartWire[10000] = {0.};		//<---Storing the Start Wire for all clusters
    float AllCluster_StartTime[10000] = {0.};		//<---Storing the Start Time for all clusters
    float AllCluster_EndWire[10000] = {0.};		//<---Storing the End Wire for all clusters
    float AllCluster_EndTime[10000] = {0.};		//<---Storing the End Time for all clusters
    
    
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
		// unsigned int nplanes = geom->Cryostat(cstat).TPC(tpc).Nplanes();
		
		
    		
    		// #################################
		// ### Loop over the wire planes ###
		// #################################
		for (unsigned int i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
			{
    			//           geom->WirePitch(Wire1, Wire2, Plane#, TPC#, Cyro#);
    			// double wire_pitch = geom->WirePitch(0,1,i,tpc,cstat); // currently unused
    			
			// ##############################################
	  		// ### If there is at least one cluster found ###
	  		// ##############################################
	  		if (Cls[i].size() >= 1)
				{
				// ##############################
	    			// ### Loop over each cluster ###
	    			// ##############################
	    			for (unsigned j = 0; j<Cls[i].size(); ++j)
					{
					
					// === Current Clusters Plane ===
					Clu_Plane[nClustersFound]  		  = clusters[Cls[i][j]]->View();
					AllCluster_Plane[AllCluster]		  = Clu_Plane[nClustersFound];
					
					// === Current Clusters StartPos ===
					Clu_StartPos_Wire[nClustersFound]	  = clusters[Cls[i][j]]->StartPos()[0];
					AllCluster_StartWire[AllCluster]	  = Clu_StartPos_Wire[nClustersFound];
					Clu_StartPos_TimeTick[nClustersFound]	  = clusters[Cls[i][j]]->StartPos()[1];
					AllCluster_StartTime[AllCluster]	  = Clu_StartPos_TimeTick[nClustersFound];
					
					// === Current Clusters EndPos ===
					Clu_EndPos_Wire[nClustersFound]		  = clusters[Cls[i][j]]->EndPos()[0];
					AllCluster_EndWire[AllCluster]		  = Clu_EndPos_Wire[nClustersFound];
					Clu_EndPos_TimeTick[nClustersFound]	  = clusters[Cls[i][j]]->EndPos()[1];
					AllCluster_EndTime[AllCluster]		  = Clu_EndPos_TimeTick[nClustersFound];
					
					// === Current Clusters Slope (In Wire and Time Tick)
					Clu_Slope[nClustersFound] 		  = dtdwstart[Cls[i][j]];
					
					
					Clu_Length[nClustersFound] 		  = std::sqrt(pow((clusters[Cls[i][j]]->StartPos()[0]-clusters[Cls[i][j]]->EndPos()[0])*13.5,2) //<---JA: Note this 13.5 is a hard coded
					 						    +pow(clusters[Cls[i][j]]->StartPos()[1]-clusters[Cls[i][j]]->EndPos()[1],2));	// number that came from ArgoNeuT MC...fix me!
					
					AllCluster_Length[AllCluster]		  = Clu_Length[nClustersFound];
					// ######################################################
					// ### Given a slope and a point find the y-intercept ###
					// ###                   c = y-mx                     ###
					// ######################################################
					Clu_Yintercept[nClustersFound] = Clu_StartPos_TimeTick[nClustersFound] - (Clu_Slope[nClustersFound] * Clu_StartPos_Wire[nClustersFound]);
					
					// #################################################################
					// ###     Also calculating the y-intercept but using the        ###
					// ###   end time of the cluster correct for the possibility     ###
					// ### that the clustering didn't get start and end points right ###
					// #################################################################
					Clu_Yintercept2[nClustersFound] = Clu_EndPos_TimeTick[nClustersFound] - (Clu_Slope[nClustersFound] * Clu_EndPos_Wire[nClustersFound]);
					
					
					nClustersFound++;
					AllCluster++;
    					}//<---End looping over clusters
					
					
    				}//<---End checking that we have at least one cluster found
			// ################################################################
			// ## If no clusters were found then put in dummy vertex values ###
			// ################################################################
			else
				{
	    			vtx_wire.push_back(-1);
	    			vtx_time.push_back(-1);
				vtx_plane.push_back(-1);
	  			}//<---End no clusters found else statement
			
			
			
			
			// ################################################################################
			// ### Now we try to find a 2-d vertex in the plane we are currently looking in ###
			// ################################################################################
			for (unsigned int n = nClustersFound; n > 0; n--)
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
						
						
							
						
						// #####################################################
						// ### Skipping calculating for end points to reduce ###
						// ###    the number of candidate verticies found    ###
						// #####################################################
						
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
						if(intersection_X2 > geom->Nwires(Clu_Plane[m],tpc,cstat)){intersection_X2 = -999;}
						if(intersection_Y2 < 0){intersection_Y2 = -999;}
						if(intersection_Y2 > detprop->NumberTimeSamples() ){intersection_Y2 = -999;}
						
						if(intersection_X < 1){intersection_X = -999;}
						if(intersection_X > geom->Nwires(Clu_Plane[m],tpc,cstat)){intersection_X = -999;}
						if(intersection_Y < 0){intersection_Y = -999;}
						if(intersection_Y > detprop->NumberTimeSamples() ){intersection_Y = -999;}
						
						
						// ### If the intersection point is 80 or more wires away from either cluster
						// ### and one of the clusters has fewer than 8 hits the intersection
						// ### is likely a crap one and we won't save this point 
						
						// ### Gathering the hits associated with the current cluster ###
						std::vector< art::Ptr<recob::Hit> > hitClu1 = fmh.at(m);
						std::vector< art::Ptr<recob::Hit> > hitClu2 = fmh.at(n);
						
						if( (abs( Clu_EndPos_Wire[m] - intersection_X2 ) > 80 && hitClu1.size() < 8) ||
						    (abs( Clu_EndPos_Wire[n] - intersection_X2 ) > 80 && hitClu2.size() < 8) )
						    {
						    intersection_X2 = -999;
						    intersection_Y2 = -999;
						    }
						    
						if( (abs( Clu_StartPos_Wire[m] - intersection_X ) > 80 && hitClu1.size() < 8) ||
						    (abs( Clu_StartPos_Wire[n] - intersection_X ) > 80 && hitClu2.size() < 8) )
						    {
						    intersection_X = -999;
						    intersection_Y = -999;
						    }
						
						// ### If the intersection point is 50 or more wires away from either cluster
						// ### and the one of the clusters has fewer than 3 hits the intersection
						// ### is likely a crap one and we won't save this point 
						if( (abs( Clu_EndPos_Wire[m] - intersection_X2 ) > 50 && hitClu1.size() < 4) ||
						    (abs( Clu_EndPos_Wire[n] - intersection_X2 ) > 50 && hitClu2.size() < 4) )
						    {
						    intersection_X2 = -999;
						    intersection_Y2 = -999;
						    }
						    
						if( (abs( Clu_StartPos_Wire[m] - intersection_X ) > 50 && hitClu1.size() < 4) ||
						    (abs( Clu_StartPos_Wire[n] - intersection_X ) > 50 && hitClu2.size() < 4) )
						    {
						    intersection_X = -999;
						    intersection_Y = -999;
						    }     
						    
						
						// ##################################################
						// ### Now checking if this intersection point is ###
						// ###       near 1 (???) or more hits            ###
						// ##################################################
						
						int nOverlapHitsEndPoint   = 0;
						int nOverlapHitsStartPoint = 0;
						// =========================
						// === Looping over hits ===
						for(size_t nHits = 0; nHits < hitListHandle->size(); nHits++)
							{
							
							// === Finding Wire/Plane associated with the hit ===
       							art::Ptr<recob::Hit> hit(hitListHandle, nHits);
       							
       							//w = hit->WireID().Wire; 
							
							// ===============================================================
							// === Skipping hits that aren't in the plane with the cluster ===
							// ===============================================================
							if(Clu_Plane[m] != hit->WireID().Plane)	{continue;}
							
							double hit_wire = hit->WireID().Wire;
							double hit_time = hit->PeakTime();
							
							if( fabs(intersection_X2 - hit_wire) < 4 && fabs(intersection_Y2 - hit_time) < 20)
								{
								nOverlapHitsEndPoint++;

								
								}//<--End checking the hit overlap with the intersection
								
							if( fabs(intersection_X - hit_wire) < 4 && fabs(intersection_Y - hit_time) < 20)
								{
								nOverlapHitsStartPoint++;
								}//<--End checking the hit overlap with the intersection
							
							
							}//<---End nHits loop
						
						
						// ### If there was no hit found near this intercept point this is likely ###
						// ###                          not a 2-d vertex                          ###
						if(nOverlapHitsEndPoint <1)
							{
							intersection_X2 = -999;
							intersection_Y2 = -999;
							}
							
						if(nOverlapHitsStartPoint < 1)
							{
							intersection_X = -999;
							intersection_Y = -999;
							}

						
						// ##########################################################
						// ### Filling the vector of Vertex Wire, Time, and Plane ###
						// ##########################################################
						
						// -----------------------------------------------------------------------------
						// --- Skip this vertex if the X  and Y intersection is outside the detector ---
						// --- using geom->Nwires(plane,tpc,cyrostat) & detprop->NumberTimeSamples() ---
						// -----------------------------------------------------------------------------
						if( intersection_X2 > 1 && 
						    intersection_Y2 > 0 && 
						  ( intersection_X2 < geom->Nwires(Clu_Plane[m],tpc,cstat) ) &&
						  ( intersection_Y2 < detprop->NumberTimeSamples() ) )
							{
							
							
							vtx_wire.push_back(intersection_X2);
							vtx_time.push_back(intersection_Y2);
							vtx_plane.push_back(Clu_Plane[m]);
							n2dVertexCandidates++;
							}//<---End saving a "good 2d vertex" candidate	
							
						// -----------------------------------------------------------------------------
						// --- Skip this vertex if the X  and Y intersection is outside the detector ---
						// --- using geom->Nwires(plane,tpc,cyrostat) & detprop->NumberTimeSamples() ---
						// -----------------------------------------------------------------------------
						if( intersection_X > 1 && 
						    intersection_Y > 0 && 
						  ( intersection_X < geom->Nwires(Clu_Plane[m],tpc,cstat) ) &&
						  ( intersection_Y < detprop->NumberTimeSamples() ) )
							{
							
							
							vtx_wire.push_back(intersection_X);
							vtx_time.push_back(intersection_Y);
							vtx_plane.push_back(Clu_Plane[m]);
							n2dVertexCandidates++;
							}//<---End saving a "good 2d vertex" candidate
						
						
						}//<---End making sure we are in the same plane
					}//<---End m ++ loop
				}//<--- End n-- loop
			
			// === Zeroing the clusters for the current view
			nClustersFound = 0;
    			}//<---End loop over wireplanes
    		}//<---End looping over TPC's
    	}//<---End looping over cryostats


    // ########################################################################
    // ### Introducing a merge step for the candidate 2-d vertex candidates ###
    // ########################################################################
    double Wire[100000]  = {0.};
    double Time[100000]  = {0.};
    double Plane[100000] = {0.};
    
    for(int vtxloop = 0 ; vtxloop < n2dVertexCandidates; vtxloop++)
	{

	Wire[vtxloop]  = vtx_wire[vtxloop];
	Time[vtxloop]  = vtx_time[vtxloop];
	Plane[vtxloop] = vtx_plane[vtxloop];
	
	}
    
    double vtx_wire_merged[100000]  = {0.};
    double vtx_time_merged[100000]  = {0.};
    double vtx_plane_merged[100000] = {0.};
    
    bool merged = false;
    
    int n2dMergedVertices = 0;
    

    // #########################################
    // ### Looping over 2d-verticies (loop1) ###
    // #########################################
    for(int vtxloop1 = 0 ; vtxloop1 < n2dVertexCandidates; vtxloop1++)
	{
	if(Wire[vtxloop1] < 0){continue;}
	
	merged = false;
	// #########################################
	// ### Looping over 2d-verticies (loop2) ###
	// #########################################
	for(int vtxloop2 = vtxloop1+1; vtxloop2 < n2dVertexCandidates; vtxloop2++)
    		{
		if(Wire[vtxloop2] < 0){continue;} 
		
		// ########################################################
		// ### Make sure the 2d-Verticies are in the same plane ###
		// ########################################################
		if(Plane[vtxloop1] == Plane[vtxloop2])
			{
			// ###############################################			
			// ### Considering merging 2d vertices if they ###
			// ###    are within 3 wires of each other     ###
			// ###############################################
			if( fabs(Wire[vtxloop1] - Wire[vtxloop2]) < 4 )
				{
				// ############################################################
				// ### Merge the verticies if they are within 10 time ticks ###
				// ############################################################
				if( fabs(Time[vtxloop1] - Time[vtxloop2]) < 10 )
					{
					vtx_wire_merged[n2dMergedVertices] = ((Wire[vtxloop2] + Wire[vtxloop1])/ 2) ;
					vtx_time_merged[n2dMergedVertices] = ((Time[vtxloop2] + Time[vtxloop1])/ 2) ;
					vtx_plane_merged[n2dMergedVertices] = Plane[vtxloop1];
					
					merged = true;
					n2dMergedVertices++;
					
					if(vtxloop2<n2dVertexCandidates)
						{vtxloop2++;}
					if(vtxloop1<n2dVertexCandidates)
						{vtxloop1++;}	
					
					
					}//<---End the check within 10 time ticks for merging
				}//<---Looking at vertices that are within 3 wires of each other
			}//<---Only looking at vertices that are in the same plane
		
		}//<---End vtxloop2
	if(!merged)
		{
		vtx_wire_merged[n2dMergedVertices]  = Wire[vtxloop1];
		vtx_time_merged[n2dMergedVertices]  = Time[vtxloop1] ;
		vtx_plane_merged[n2dMergedVertices] = Plane[vtxloop1];
		n2dMergedVertices++;
		}//<---end saving unmerged verticies
			
	
	}//<---End vtxloop1
	


    
    double y_coord = 0, z_coord = 0;
    
    double x_3dVertex[100000] = {0.}, y_3dVertex[100000] = {0.}, z_3dVertex[100000] = {0.};
    int n3dVertex = 0;
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
    		for(unsigned int vtx = n2dMergedVertices; vtx > 0; vtx--)
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
					// === JA: Need to include vtx tpc, and cstat to make detector agnositc
			
					unsigned int vtx1_plane   = vtx_plane_merged[vtx1];
					unsigned int vtx1_wire    = vtx_wire_merged[vtx1];
					try
						{
						vtx1_channel = geom->PlaneWireToChannel(vtx1_plane, vtx1_wire, tpc, cstat);
				
						}
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "PlaneWireToChannel Failed";
	      
	      					continue;
	    					}
				
				
			
					unsigned int vtx2_plane   = vtx_plane_merged[vtx];
					unsigned int vtx2_wire    = vtx_wire_merged[vtx];
			
					try
						{
						vtx2_channel = geom->PlaneWireToChannel(vtx2_plane, vtx2_wire, tpc, cstat);
				
						}
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "PlaneWireToChannel Failed";
	      
	      					continue;
				
						}
			
					// ##############################################################################
					// ### Check to see if the channels intersect and save the y and z coordinate ###
					// ##############################################################################
					bool match = false;
			
					try
						{
						match = geom->ChannelsIntersect( vtx1_channel, vtx2_channel, y_coord, z_coord);
				
				
						}
					catch(...)
						{
	      					mf::LogWarning("FeatureVertexFinder") << "match failed for some reason";
	      
	      					match = false;
	      					continue;
	    					}
					if( match )
						{
				
						// #############################################################################
						// ### Now check if the matched channels are within 1 cm when projected in X ###
						// #############################################################################
						float tempXCluster1 = detprop->ConvertTicksToX(vtx_time_merged[vtx1], vtx1_plane, tpc, cstat);
						float tempXCluster2 = detprop->ConvertTicksToX(vtx_time_merged[vtx], vtx2_plane, tpc, cstat);
						
						if(std::abs(tempXCluster1 - tempXCluster2) < 0.5)
							{
							//        detprop->ConvertTicksToX(ticks, plane, tpc, cryostat)
							x_3dVertex[n3dVertex] = detprop->ConvertTicksToX(vtx_time_merged[vtx1], vtx_plane_merged[vtx1], tpc, cstat); 
							y_3dVertex[n3dVertex] = y_coord;
							z_3dVertex[n3dVertex] = z_coord;
							n3dVertex++;
					
							}//<---End Checking if the vertices agree "well enough" in time tick
				
						}//<---End Checking if verticies intersect

					}//<--- End checking we are in different planes
				}//<---end vtx1 for loop
			}//<---End vtx for loop
		}//<---End loop over TPC's 
	}//<---End loop over cryostats




//-----------------------------------------------------------------------------------------------------------------------------

// ###################################################################
// ### Now we need to make sure that we remove duplicate verticies ###
// ###    just in case we have any in our list of n3dVertex        ###
// ###################################################################

double x_3dVertex_dupRemoved[100000] = {0.}, y_3dVertex_dupRemoved[100000] = {0.}, z_3dVertex_dupRemoved[100000] = {0.};
double strength_3dVertex_dubRemoved[10000] = {0.};
int n3dVertex_dupRemoved = 0;

for(int dup = 0; dup < n3dVertex; dup ++)
	{
	float tempX_dup = x_3dVertex[dup];
	float tempY_dup = y_3dVertex[dup];
	float tempZ_dup = z_3dVertex[dup];
	
	
	bool duplicate_found = false;
	
	for(int check = dup+1; check < n3dVertex; check++)
		{
		
		// #############################################################################
		// ### I am going to call a duplicate vertex one that matches in x, y, and z ###
		// ###        within 0.1 cm for all 3 coordinates simultaneously             ###
		// #############################################################################
		if(std::abs( x_3dVertex[check] - tempX_dup ) < 0.1 && std::abs( y_3dVertex[check] - tempY_dup ) < 0.1 &&
		   std::abs( z_3dVertex[check] - tempZ_dup ) < 0.1 )
		   	{
			duplicate_found = true;
			
			}//<---End checking to see if this is a duplicate vertex
	
	
	
		}//<---End check for loop
	
	// ######################################################################
	// ### If we didn't find a duplicate then lets save this 3d vertex as ###
	// ###            a real candidate for consideration                  ###
	// ######################################################################
	if(!duplicate_found)
		{
		x_3dVertex_dupRemoved[n3dVertex_dupRemoved] = tempX_dup;
		y_3dVertex_dupRemoved[n3dVertex_dupRemoved] = tempY_dup;
		z_3dVertex_dupRemoved[n3dVertex_dupRemoved] = tempZ_dup;
		
		n3dVertex_dupRemoved++;
		}	

	}//<---End dup for loop



//-----------------------------------------------------------------------------------------------------------------------------
// ---- Adding a section where I pair together 3d feature points with 3d cluster vertex candidates and see how many I have ----

double pVtx_X[1000] = {0.}, pVtx_Y[1000] = {0.}, pVtx_Z[1000] = {0.}, pVtx_Stgth[1000] = {0.};

int nprimaries = 0;
for(int nCV = 0; nCV < n3dVertex_dupRemoved; nCV++)
	{
	for(int nFV = 0; nFV <  n3dFeatures; nFV++)
		{
		if(std::abs(x_feature[nFV] - x_3dVertex_dupRemoved[nCV]) < 1 &&
		   std::abs(y_feature[nFV] - y_3dVertex_dupRemoved[nCV]) < 1 &&
		   std::abs(z_feature[nFV] - z_3dVertex_dupRemoved[nCV]) < 1 )
		   	{
		  	
			pVtx_X[nprimaries]     = x_3dVertex_dupRemoved[nCV];
			pVtx_Y[nprimaries]     = y_3dVertex_dupRemoved[nCV];
			pVtx_Z[nprimaries]     = z_3dVertex_dupRemoved[nCV];
			pVtx_Stgth[nprimaries] = strength_feature[nFV];
			nprimaries++;
			

			}
	
		}//<---End 
	}//<---End nCV loop 



// ##################################################################################################
// ### Sorting the "primary" verticies I have found such that the first in the list is the vertex ###
// ###              with the highest vertex strength and the lowest z location                    ###
// ##################################################################################################

int flag = 1;
double tempX, tempY, tempZ, tempS;

for(int npri = 0; (npri < nprimaries) && flag; npri++)
	{
	flag = 0;
	for(int mpri = 0; mpri < nprimaries -1; mpri++)
		{
		// Swap the order of the two elements
		if(pVtx_Stgth[mpri+1] > pVtx_Stgth[mpri] || (pVtx_Stgth[mpri+1] == pVtx_Stgth[mpri] && pVtx_Z[mpri] > pVtx_Z[mpri+1]) )
			{
			tempX = pVtx_X[mpri];
			pVtx_X[mpri] = pVtx_X[mpri+1];
			pVtx_X[mpri+1] = tempX;
			
			tempY = pVtx_Y[mpri];
			pVtx_Y[mpri] = pVtx_Y[mpri+1];
			pVtx_Y[mpri+1] = tempY;
			
			tempZ = pVtx_Z[mpri];
			pVtx_Z[mpri] = pVtx_Z[mpri+1];
			pVtx_Z[mpri+1] = tempZ;
			
			tempS = pVtx_Stgth[mpri];
			pVtx_Stgth[mpri] = pVtx_Stgth[mpri+1];
			pVtx_Stgth[mpri+1] = tempS;
			
			flag = 1;
			
			}//<---Inside swap loop
		}//<---End mpri	
	}//<---End npri loop
	


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------

// ########################################################################################################
// ### Case 1) IF WE HAVE "PRIMARY VERTICIES" (Cluster verticies which are matched across all planes and 
// ###         paired with a 3d feature candidate and sorted into an descending order with the first 
// ###         vertex in the list likely the primary vertex
// ########################################################################################################
if(nprimaries > 0)
	{
	// ######################################
	// ### Looping over Primary Verticies ###
	// ######################################
	for(int pri = 0; pri < nprimaries; pri++)
		{
		
		// ###############################################
		// ### Push each primary vertex onto the event ###
		// ###############################################
		double tempxyz[3] = {pVtx_X[pri], pVtx_Y[pri], pVtx_Z[pri]};
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
    			// ##########################
      			// ### Looping over TPC's ###
      			// ##########################
      			for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
				{
    				// #################################
				// ### Loop over the wire planes ###
				// #################################
				for (size_t plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane)
					{
					double temp2dXYZ[3] = {pVtx_X[pri], pVtx_Y[pri], pVtx_Z[pri]};
					double temp2dStrength = pVtx_Stgth[pri];
					
					// ######################################################################
					// ### Converting the 3d vertex into 2d time ticks, wire, and channel ###
					// ######################################################################
					double EndPoint2d_TimeTick = detprop->ConvertXToTicks(temp2dXYZ[0],plane, tpc, cstat);
					int EndPoint2d_Wire = 0;
					int EndPoint2d_Channel     = 0;
					
					// ### Putting in protection in case NearestWire Fails ###
					try
					
						{EndPoint2d_Wire 	   = geom->NearestWire(temp2dXYZ , plane, tpc, cstat);}
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      					continue;
						
						}
					// ### Putting in protection in case NearestChannel Fails ###
					try
						{EndPoint2d_Channel     = geom->NearestChannel(temp2dXYZ, plane, tpc, cstat);}
						
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      					continue;
						
						}
					// ### Making geo::WireID and getting the current View number ###	
					geo::View_t View	   = geom->View(EndPoint2d_Channel);
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
					
					
					}//<---End plane loop
				}//<---End tpc loop
			}//<---End cstat loop
		}//<---End primary loop	
	
	}//<---End nprimaries > 0 check
       

      

//----------------------------------------------------------------------------------------------------------------------------- 
// ##################################################################################################################################################
// ### Case 2) Here we have at least 1 3d candidate cluster vertex found but no matches to a 3d feature point (and thus no primary vertex found)
// ###         In this case we will try to wittle down the number of cluster verticies found by merging them together...once merged we will
// ### 	       store them onto the event (right now in any old order....I will eventually sort these in a smarter fashion
// ##################################################################################################################################################
int totalGood = 0;
int nMerges = 0;
double x_good[10000] = {0.}, y_good[10000] = {0.}, z_good[10000] = {0.}, strength_good[10000] = {0.};
	
int totalGood2 = 0;
double x_good2[10000] = {0.}, y_good2[10000] = {0.}, z_good2[10000] = {0.};

int match2planes = 0;
int TwoDvertexStrength = 0;
	
if (n3dVertex_dupRemoved > 0 && nprimaries == 0)
	{
	// ######################################################
	// ###     Setting a limit to the number of merges    ###
	// ### to be 3 times the number of 3d verticies found ###
	// ######################################################
	int LimitMerge = 3 * n3dVertex_dupRemoved;
	// ### Trying to merge nearby verticies found ###
	for( int merge1 = 0; merge1 < n3dVertex_dupRemoved; merge1++)
		{
		TwoDvertexStrength = 1;  
		for( int merge2 = n3dVertex ; merge2 > merge1; merge2--)
			{
			double temp1_x = x_3dVertex[merge1];
			double temp1_y = y_3dVertex[merge1];
			double temp1_z = z_3dVertex[merge1];
			
			double temp2_x = x_3dVertex[merge2];
			double temp2_y = y_3dVertex[merge2];
			double temp2_z = z_3dVertex[merge2];
				
			if( temp1_x == 0 || temp1_y == 0|| temp1_z == 0 ||
			    temp2_x == 0 || temp2_y == 0 || temp2_z == 0) {continue;}
				    
				
			// ### Merge the verticies if they are within 1.0 cm of each other ###
			if ( (std::abs( temp1_x - temp2_x ) < 2.0 && temp1_x != 0 && temp2_x !=0) &&
			     (std::abs( temp1_y - temp2_y ) < 2.0 && temp1_y != 0 && temp2_y !=0) &&
			     (std::abs( temp1_z - temp2_z ) < 2.0 && temp1_z != 0 && temp2_z !=0) &&
			      nMerges < LimitMerge)
			    	{
				nMerges++;
				// ### Zero the vertex that I am merging ###
				x_3dVertex[merge2] = 0.0;
				y_3dVertex[merge2] = 0.0;
				z_3dVertex[merge2] = 0.0;
					
				n3dVertex_dupRemoved++;
					
				// ### Add the merged vertex to the end of the vector ###
				x_3dVertex[n3dVertex_dupRemoved] = (temp1_x + temp2_x)/2;
				y_3dVertex[n3dVertex_dupRemoved] = (temp1_y + temp2_y)/2;
				z_3dVertex[n3dVertex_dupRemoved] = (temp1_z + temp2_z)/2;
					
				// ######################################################################
				// ### If we merged the verticies then increase its relative strength ###
				// ######################################################################
				TwoDvertexStrength++;
				strength_3dVertex_dubRemoved[n3dVertex_dupRemoved] = TwoDvertexStrength;
					
					
				}//<---End merging verticies
			}//<---End merge2 loop
		}//<---End merge1 loop
			
	// ######################################################################
	// ### Now loop over the new list and only save the non-zero vertices ###
	// ######################################################################
	for (int goodvtx = 0; goodvtx < n3dVertex; goodvtx++)
		{
		if( x_3dVertex[goodvtx] != 0.0 && 
		    y_3dVertex[goodvtx] != 0.0 &&
		    z_3dVertex[goodvtx] != 0.0 )
		    	{
			bool duplicate = false;
			// ###############################################
			// ### Check to make sure this isn't a copy of ###
			// ###        a previously found vertex        ###
			// ###############################################
			for (int check = goodvtx+1; check > 0; check--)
				{
				// #########################################################################
				// ### Adding a protection that we don't check the vertex against itself ###
				// #########################################################################
				if (check != goodvtx)
					{
					// ### check if this vertex exists in the list ###
					if (x_3dVertex[goodvtx] == x_3dVertex[check] &&
					    y_3dVertex[goodvtx] == y_3dVertex[check] &&
					    z_3dVertex[goodvtx] == z_3dVertex[check] )
					    {
					    duplicate = true;
					    
					    }//<---End duplicate 
					}//<---Don't check the vertex against itself
					
				}//<---end check loop
			if(!duplicate)
				{	
				x_good[totalGood] =  x_3dVertex[goodvtx];
				y_good[totalGood] =  y_3dVertex[goodvtx];
				z_good[totalGood] =  z_3dVertex[goodvtx];
				strength_good[totalGood] = strength_3dVertex_dubRemoved[goodvtx];
				totalGood++;
				}//<---End removing duplicates
			
			}//<---End making sure the vetex isn't equal to zero
		}//<---End goodvtx loop
		
		
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ###################################################
	// ### Loop over good verticies to only save those ###
	// ###       that match some number of hits        ###
	// ###################################################
		
	// ##############################
	// ### Looping over verticies ###
	// ##############################
	for(int ll = 0; ll < totalGood; ll++)
		{
			
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
				for (size_t i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
					{
					double xyz[3] = {x_good[ll], y_good[ll], z_good[ll]};
					double EndPoint2d_TimeTick = detprop->ConvertXToTicks(x_good[ll],i, tpc, cstat);
					int EndPoint2d_Wire = 0;
					try
				
					{EndPoint2d_Wire 	   = geom->NearestWire(xyz , i, tpc, cstat);}
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      					continue;
						
						}

					// #############################################################
					// ### Only going to save the vertex if it is within 3 wires ###
					// ###       and 20 time ticks of a hit in the event         ###
					// #############################################################
						
					int matchedVtxToHit = 0;
					// =========================
					// === Looping over hits ===
					for(size_t nHits = 0; nHits < hitListHandle->size(); nHits++)
						{
							
						// === Finding Wire/Plane associated with the hit ===
       						art::Ptr<recob::Hit> hit(hitListHandle, nHits);
       							
							
						// ===============================================================
						// === Skipping hits that aren't in the plane with the cluster ===
						// ===============================================================
						if(i != hit->WireID().Plane)	{continue;}
							
						double hit_wire = hit->WireID().Wire;
						double hit_time = hit->PeakTime();
							
						if( fabs(EndPoint2d_Wire - hit_wire) < 4 && fabs(EndPoint2d_TimeTick - hit_time) < 20)
							{
							matchedVtxToHit++;
								
								
							}
								
						}//<---End loop over hits
						
					// #####################################################
					// ### Saving the 2d Vertex found if matched to more ###
					// ### hits than what you would expect along a track ###
					// #####################################################
					if(matchedVtxToHit > 5)
						{
						match2planes++;
						}	

					}//<---End loop over Planes
				}//<---End loop over tpc's
			}//<---End loop over cryostats
		// #########################################
		// ### Saving the 3d vertex if there was ###
		// ###    a match in at least 2 planes   ###
		// #########################################
		if(match2planes > 1)
			{
			x_good2[totalGood2] = x_good[ll];
			y_good2[totalGood2] = y_good[ll];
			z_good2[totalGood2] = z_good[ll];
			totalGood2++;
			
			}	

		match2planes = 0;
		
		}//<---End ll loop		
		
	// ##############################
	// ### Looping over verticies ###
	// ##############################
	for(int l = 0; l < totalGood2; l++)
		{
		// #############################
		// ### Looping over features ###
		// #############################
		for (int f = 0; f < n3dFeatures; f++)
			{
			// ###########################################################
			// ### Looking for features and verticies within 3cm in 3d ###
			// ###########################################################
			if(std::abs(x_good2[l] - x_feature[f]) <= 3 &&  
			   std::abs(y_good2[l] - y_feature[f]) <= 3 &&
			   std::abs(z_good2[l] - z_feature[f]) <= 3)
			   	{
					
				// ##########################################################
				// ### If there is also a 3d feature near this point then ###
				// ###      increase the vertex strength again            ###
				// ##########################################################
					
				strength_good[l] = strength_good[l]+1;
					
				}//<---End finding a match between features and verticies 
			}//<--- End feature for loop
				
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
				for (size_t i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
					{
					double xyz[3] = {x_good2[l], y_good2[l], z_good2[l]};
					double EndPoint2d_TimeTick = detprop->ConvertXToTicks(x_good2[l],i, tpc, cstat);
					int EndPoint2d_Wire = 0;
					int EndPoint2d_Channel = 0;
					try
				
						{EndPoint2d_Wire 	   = geom->NearestWire(xyz , i, tpc, cstat);}
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	    
	      					continue;				
						}
							
					try
						{EndPoint2d_Channel     = geom->NearestChannel(xyz, i, tpc, cstat);}
					catch(...)
						{
						mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
      						continue;
						
						}
							
					geo::View_t View	   = geom->View(EndPoint2d_Channel);
					geo::WireID wireID(cstat,tpc,i,EndPoint2d_Wire);
						
					int EndPoint2dstrength = strength_good[l];
					
					recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
								  wireID ,		//<---geo::WireID
								  EndPoint2dstrength ,	//<---Vtx strength (JA: ?)
								  epcol->size() ,	//<---Vtx ID (JA: ?)
								  View ,		//<---Vtx View 	
								  1 );			//<---Vtx Total Charge (JA: Need to figure this one?)
					epcol->push_back(vertex);
							

					}//<---End loop over Planes
				}//<---End loop over tpc's
			}//<---End loop over cryostats
		// #########################################
		// ### Saving the 3d vertex if there was ###
		// ###    a match in at least 2 planes   ###
		// #########################################
			
		double xyz2[3] = {x_good2[l], y_good2[l], z_good2[l]};
		recob::Vertex the3Dvertex(xyz2, vcol->size());
		vcol->push_back(the3Dvertex);
				
		}//<--End l for loop
		
	}//<---End Case 2, many 3d Verticies found
		
		
		
//--------------------------------------------------------------------------------------------------------------------------------------------- 

// ##################################################################################################################	
// ### Case3) If we have no Primary Verticies and no Cluster Verticies found we are going to try a bunch of 
// ###	      bail out techniques including:
// ###		1) Seeing if there is a feature matched to the longest cluster in the planes start point
// ###		2) Just taking the longest cluster in the plane's Start Point and making a 3d point out of that
// ###		3) Finally if all that fails just take the strongest 3d feature point and save that to the event
// ###	THIS SECTION COULD PROBABLY USE THE MOST WORK!!!!
// ##################################################################################################################
	      	
	double x_featClusMatch = 0, y_featClusMatch = 0, z_featClusMatch = 0;
	
	double LongestClusterStartWire[100] = {0};
	double LongestClusterStartTime[100] = {0};
	if (n3dVertex_dupRemoved == 0)
		{
		bool NothingFoundYet = true;
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
				for (size_t plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane)
					{
					// ############################################################################
					// ### We reinitialize the strength for each plane since we haven't checked ###
					// ###  at the top of this loop if this has been matched to anything else   ###
					// ############################################################################
					TwoDvertexStrength = 0;
					
					// === Defining a bunch of temp variables ===
					float Length_LongestClusterInPlane = 0;
					float StartWire_LongestClusterInPlane = 0;
					float StartTime_LongestClusterInPlane = 0;
					float EndWire_LongestClusterInPlane = 0;
					float EndTime_LongestClusterInPlane = 0;
					//float Plane_LongestClusterInPlane = 0;
					// #######################################################################
					// ### Looping over clusters to find the longest cluster in each plane ###
					// #######################################################################
					for (int clus = 0; clus < AllCluster; clus++)
						{
						// ### Checking that the current cluster is in the current plane ###
						if(AllCluster_Plane[clus] == plane)
							{
							
							// Storing the current cluster length
							float tempClusterLength = AllCluster_Length[clus];
							

							
							// Is the current cluster the longest one found?
							if( tempClusterLength > Length_LongestClusterInPlane)
								{
								// Save the longest clusters information in this view
								Length_LongestClusterInPlane    = tempClusterLength;
								StartWire_LongestClusterInPlane = AllCluster_StartWire[clus];
								StartTime_LongestClusterInPlane = AllCluster_StartTime[clus];
								EndWire_LongestClusterInPlane   = AllCluster_EndWire[clus];
								EndTime_LongestClusterInPlane   = AllCluster_EndTime[clus];
								//Plane_LongestClusterInPlane     = AllCluster_Plane[clus];
								}//<---End saving the longest cluster
							}//<---Current cluster in current plane	
						}///<---End clus for loop

					
					// #####################################################################################
					// ###       Now that we have the longest cluster in this plane lets loop over       ###
					// ### the 3d-Feature points found and see if any are consistant with this 2-d point ###
					// #####################################################################################
					
					NothingFoundYet = true;
					
					// ### Set the score for this 2d vertex in this plane to be = 0 ###
					// ### This will increase for each feature it finds that matches ###
					// ### down to it
					
					/// #### Get the 3d feature point, project it back down to the current plane...check to see if consistant
					float temp_feat_score = 0;
					
					for (int feat = 0; feat < n3dFeatures; feat++)
						{
						// #########################################################
						// ### Finding the nearest wire for this 3dfeature point ###
						// #########################################################
						float xyzfeat[3]   = {x_feature[feat],y_feature[feat],z_feature[feat]};
						float xyzfeat_score = strength_feature[feat];
						double feat_2dWire = -1;
						double feat_TimeT = -1;
						
						try
							{
							feat_2dWire = geom->NearestWire(xyzfeat, plane, tpc, cstat);
							feat_TimeT  = detprop->ConvertXToTicks(x_feature[feat], plane, tpc, cstat);
							}
						catch(...)
							{
							mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
							continue;
							}

						// ######################################################################################
						// ###           Now check to see if this feature point is within 2 wires and         ###
						// ###   5 time ticks of the start time of the longest cluster found in this plane    ###
						// ###                and is the strongest xyz feature in the list                    ###
						// ######################################################################################
						if( std::abs(feat_2dWire - StartWire_LongestClusterInPlane) <= 2 && std::abs(StartTime_LongestClusterInPlane - feat_TimeT ) <=5 &&
						    xyzfeat_score >= temp_feat_score  )
							{
							temp_feat_score = xyzfeat_score;
							// ###########################################################
							// ### If a feature matches the 2d point increase strength ###
							// ###########################################################
							TwoDvertexStrength++;

							
							if(NothingFoundYet)
								{
								
								x_featClusMatch = xyzfeat[0];
								y_featClusMatch = xyzfeat[1];
								z_featClusMatch = xyzfeat[2];
								
								double xyz[3] = {xyzfeat[0], xyzfeat[1], xyzfeat[2]};
								double EndPoint2d_TimeTick = detprop->ConvertXToTicks(xyzfeat[0],plane, tpc, cstat);
						
								int EndPoint2d_Wire 	   = 0;
								int EndPoint2d_Channel     = 0;
								
								try
					
									{EndPoint2d_Wire 	   = geom->NearestWire(xyz , plane, tpc, cstat);}
								catch(...)
									{
									mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      								continue;
						
									}
							
								try
									{EndPoint2d_Channel     = geom->NearestChannel(xyz, plane, tpc, cstat);}
								catch(...)
									{
									mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	    
	      								continue;
						
									}
								geo::View_t View	   = geom->View(EndPoint2d_Channel);
								geo::WireID wireID(cstat,tpc,plane,EndPoint2d_Wire);
								
								// ### Saving the 2d Vertex found ###
								recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
										  wireID ,		//<---geo::WireID
										  TwoDvertexStrength ,			//<---Vtx strength (JA: ?)
										  epcol->size() ,	//<---Vtx ID (JA: ?)
										  View ,		//<---Vtx View 	
										  1 );			//<---Vtx charge (JA: Need to figure this one?)
								epcol->push_back(vertex);
								
								
								}//<---End saving this vertex found
							
							NothingFoundYet = false;
							}//<---End finding a match between the start time and the feature
						}//<---End feat for loop
					
					// ##################################################################################
					// ### After checking all the start times against the feature points, if we still ###
					// ###    don't have a vertex found check the end times of the longest cluster    ###
					// ###           (This should protect {some} against the clustering               ###
					// ###               putting the start and end in the wrong spot)                 ###
					// ##################################################################################
					if (NothingFoundYet )
						{
						for (int feat2 = 0; feat2 < n3dFeatures; feat2++)
							{
							// #########################################################
							// ### Finding the nearest wire for this 3dfeature point ###
							// #########################################################
							float xyzfeat2[3]   = {x_feature[feat2],y_feature[feat2],z_feature[feat2]};
							float xyzfeat_score2 = strength_feature[feat2];
							
							double feat_2dWire = -1;
							double feat_TimeT = -1;
							try
								{
								feat_2dWire = geom->NearestWire(xyzfeat2, plane, tpc, cstat);
								feat_TimeT  = detprop->ConvertXToTicks(x_feature[feat2], plane, tpc, cstat);
								}
							catch(...)
								{
								mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
								continue;
								}

							// ######################################################################################
							// ###           Now check to see if this feature point is within 2 wires and         ###
							// ###   5 time ticks of the start time of the longest cluster found in this plane    ###
							// ######################################################################################
							if( std::abs(feat_2dWire - EndWire_LongestClusterInPlane) <= 2 && std::abs(EndTime_LongestClusterInPlane - feat_TimeT ) <=5 && 
							    xyzfeat_score2 >= temp_feat_score )
								{
								
								temp_feat_score = xyzfeat_score2;
								// ###########################################################
								// ### If a feature matches the 2d point increase strength ###
								// ###########################################################
								TwoDvertexStrength++;
								
								if(NothingFoundYet)
									{
								
									x_featClusMatch = xyzfeat2[0];
									y_featClusMatch = xyzfeat2[1];
									z_featClusMatch = xyzfeat2[2];
									
									double xyz[3] = {xyzfeat2[0], xyzfeat2[1], xyzfeat2[2]};
									double EndPoint2d_TimeTick = detprop->ConvertXToTicks(xyzfeat2[0],plane, tpc, cstat);
									
									try
										{
										int EndPoint2d_Wire 	   = geom->NearestWire(xyz , plane, tpc, cstat);
										int EndPoint2d_Channel     = geom->NearestChannel(xyz, plane, tpc, cstat);
										geo::View_t View	   = geom->View(EndPoint2d_Channel);
										geo::WireID wireID(cstat,tpc,plane,EndPoint2d_Wire);
									
										// ### Saving the 2d Vertex found ###
										recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
											  wireID ,		//<---geo::WireID
											  TwoDvertexStrength ,	//<---Vtx strength (JA: ?)
											  epcol->size() ,	//<---Vtx ID (JA: ?)
											  View ,		//<---Vtx View 	
											  1 );			//<---Vtx charge (JA: Need to figure this one?)
										epcol->push_back(vertex);
										}
									catch(...)
										{
										mf::LogWarning("FeatureVertexFinder") << "2dWire failed";
	      									continue;
										
										}
								
								
									}//<---End saving this vertex found
								
								
								NothingFoundYet = false;
								}//<---End finding a match between the 
							}//<---End feat for loop
					
						}//<---End Nothing Found Yet
					
					// ####################################################################################	
					// ### After all that...if we still don't have a 2-d / 3d Vertex then we just take  ###
					// ### the starting point of the longest cluster as out 2d vertex and in each plane ###
					// ### For the 3d point lets take the Feature with the greates strength as the 3-d  ###
					// ### and for now (JA: come back to) we won't enforce that these points match      ###
					// ####################################################################################
					if (NothingFoundYet)
						{
						
						double EndPoint2d_TimeTick = StartTime_LongestClusterInPlane;
						int EndPoint2d_Wire 	   = StartWire_LongestClusterInPlane;
						int EndPoint2d_Channel     = geom->PlaneWireToChannel(plane,EndPoint2d_Wire,tpc,cstat);
						geo::View_t View	   = geom->View(EndPoint2d_Channel);
						geo::WireID wireID(cstat,tpc,plane,EndPoint2d_Wire);
						
						// ### Saving the 2d Vertex found ###
						recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
									  wireID ,		//<---geo::WireID
									  TwoDvertexStrength ,	//<---Vtx strength (JA: ?)
									  epcol->size() ,	//<---Vtx ID (JA: ?)
									  View ,		//<---Vtx View 	
									  1 );			//<---Vtx Strength (JA: Need to figure this one?)
						epcol->push_back(vertex);
						
						
						
						
						
						}//<---End last check for something found
						
					LongestClusterStartWire[plane] = StartWire_LongestClusterInPlane;
					LongestClusterStartTime[plane] = StartTime_LongestClusterInPlane;
					}//<--End looping over planes
				}//<--End looping over TPC's
			}//<---End looping over cryostats
		
		// ########################################################
		// ### Last Ditch attempt to get a 3d point (if needed) ###'
		// ########################################################
		if (NothingFoundYet)
			{
			
			// ### Taking the highest scoring feature point ###
			float tempFeatureScore = 0;
			// ### Looping over features ###
			for (int feat3 = 0; feat3 < n3dFeatures; feat3++)
				{
							
				if( strength_feature[feat3]  >  tempFeatureScore)
					{
					tempFeatureScore = strength_feature[feat3];
					x_featClusMatch = x_feature[feat3];
					y_featClusMatch = y_feature[feat3];
					z_featClusMatch = z_feature[feat3];
								
								
								
					}//<---End checking if the current feature is the strongest feature
							
				}//<---End feat3 loop
			
			
			if (x_featClusMatch != 0) {NothingFoundYet = false;}
			
			}//<---Last ditch attempt to get a 3d vertex
			
		// ### In a last, last ditch effort we will take the start point
		// ### of the two longest clusters we found in each plane and just compute
		// ### an XYZ coordinate for them...regardless of how well they match
		if (NothingFoundYet)
			{
			double temp_wire_pitch0 = geom->WirePitch(0,1,0,0,0);
			double temp_wire_pitch1 = geom->WirePitch(0,1,1,0,0);
			double Iw0 = (LongestClusterStartWire[0]+3.95)*temp_wire_pitch0;
		        double Cw0 = (LongestClusterStartWire[0]+1.84)*temp_wire_pitch1;
			double It0 = LongestClusterStartTime[0] - presamplings;
			It0 *= TimetoCm;
			double Ct0 = LongestClusterStartTime[1] - presamplings ;
			Ct0 *= TimetoCm;
			x_featClusMatch = Ct0;
			y_featClusMatch = (Cw0-Iw0)/(2.*TMath::Sin(Angle));
			z_featClusMatch = (Cw0+Iw0)/(2.*TMath::Cos(Angle))-DetectorLength/2.*TMath::Tan(Angle);
			
			}
		
		double xyz2[3] = {x_featClusMatch, y_featClusMatch, z_featClusMatch};
			recob::Vertex the3Dvertex(xyz2, vcol->size());
			vcol->push_back(the3Dvertex);	

		}//<---End case where there was no good match	
       	
    
//----------------------------------------------------------------------------------------------------------------------------- 
    

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

  } // end of produce
}


// =====================================================================================================
// =====================================================================================================

// -------------------------------------------------------------------------------------------------------------------------

// #######################################################################
// ###   Trying to incorporate ClusterParmsAlg into here to utilize    ###
// ### the finding axis of a cluster instead of doing simple pol1 fits ###
// #######################################################################

/*    // ##########################################
    // ### Start by looping over the clusters ###
    // ##########################################
    for(unsigned int iClust = 0; iClust < clusterListHandle->size(); iClust++)
    	{

    	art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
	
   	std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
	
	double lineslope, lineintercept,goodness,wire_start,time_start,wire_end,time_end;

	fClParAlg.Find2DAxisRough(lineslope,lineintercept,goodness,hitlist);
	// hmmm...seems that my hitlist is always too small....that's strange....
	
	std::cout<<std::endl;
	std::cout<<"iClust = "<<iClust<<" , StartWire = "<<cl->StartPos()[0]<<" , StartTime = "<<cl->StartPos()[1]<<std::endl;
	std::cout<<"iClust = "<<iClust<<" , EndWire   = "<<cl->EndPos()[0]<<" , EndTime   = "<<cl->EndPos()[1]<<std::endl;
	std::cout<<"iClust = "<<iClust<<" , Slope     = "<<lineslope<<" , intercept = "<<lineintercept<<" , goodness = "<<goodness<<std::endl;
	
	
	fClParAlg.Find2DStartPointsHighCharge( hitlist,wire_start,time_start,wire_end,time_end);
	
	std::cout<<std::endl;
	std::cout<<"Recalculated StartWire = "<<wire_start<<" , Recalculated StartTime = "<<time_start<<std::endl;
	std::cout<<"Recalculated EndWire   = "<<wire_end<<" , Recalculated EndTime = "<<time_end<<std::endl;


	}//<---End Cluster list

*/





// =====================================================================================================
// =====================================================================================================
namespace vertex{

  DEFINE_ART_MODULE(FeatureVertexFinder)

} // end of vertex namespace
// =====================================================================================================
// =====================================================================================================

