////////////////////////////////////////////////////////////////////////
//
// \file Track3Dreco.cxx
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port and edits: soderber,echurch@fnal.gov
//  This algorithm is designed to reconstruct 3D tracks through a simple 
//  2D-track matching algorithm
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDProducer.h"

#include <vector>
#include <string>
#include <iomanip>

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

// Framework includes
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

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

namespace trkf {
   
  class Track3Dreco : public art::EDProducer {
    
  public:
    
    explicit Track3Dreco(fhicl::ParameterSet const& pset);
    ~Track3Dreco();
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:
        
    int             ftmatch;             ///< tolerance for time matching (in time samples) 
    double          fchi2dof;            ///< tolerance for chi2/dof of cluster fit to function
    std::string     fClusterModuleLabel; ///< label for input cluster collection

  
  }; // class Track3Dreco

}

namespace trkf {

//-------------------------------------------------
Track3Dreco::Track3Dreco(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Track>                        >();
  produces< std::vector<recob::SpacePoint>                   >();
  produces< art::Assns<recob::Track,      recob::Cluster>    >();
  produces< art::Assns<recob::Track,      recob::SpacePoint> >();
  produces< art::Assns<recob::SpacePoint, recob::Hit>        >();
  produces< art::Assns<recob::Track,      recob::Hit>        >();
}

//-------------------------------------------------
Track3Dreco::~Track3Dreco()
{
}

void Track3Dreco::reconfigure(fhicl::ParameterSet const& pset)
{
  fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
  ftmatch                 = pset.get< int    >("TMatch");
  fchi2dof                = pset.get< double >("Chi2DOFmax");
}

//-------------------------------------------------
void Track3Dreco::beginJob()
{
  

}

void Track3Dreco::endJob()
{
}

//------------------------------------------------------------------------------------//
void Track3Dreco::produce(art::Event& evt)
{ 
   // get services
   art::ServiceHandle<geo::Geometry> geom;
   art::ServiceHandle<util::LArProperties> larprop;

   std::unique_ptr<std::vector<recob::Track>                    > tcol(new std::vector<recob::Track>);
   std::unique_ptr<std::vector<recob::SpacePoint>               > spacepoints(new std::vector<recob::SpacePoint>);
   std::unique_ptr< art::Assns<recob::Track, recob::Cluster>    > cassn (new art::Assns<recob::Track,      recob::Cluster>);
   std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint> > sassn (new art::Assns<recob::Track,      recob::SpacePoint>);
   std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit>   > shassn(new art::Assns<recob::SpacePoint, recob::Hit>);
   std::unique_ptr< art::Assns<recob::Track, recob::Hit>        > hassn (new art::Assns<recob::Track,      recob::Hit>);

   // define TPC parameters
   TString tpcName = geom->GetLArTPCVolumeName();

   //  double YC =  (m_TPCHalfZ-5.)*2.; // TPC height in cm
   double YC =  (geom->DetHalfHeight())*2.; // *ArgoNeuT* TPC active-volume height in cm
   double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
   // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
   double timetick = 0.198;    //time sample in us
   double presamplings = 60.;
   const double wireShift=50.; // half the number of wires from the Induction(Collection) plane intersecting with a wire from the Collection(Induction) plane.
   double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
   double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
   double Efield_drift = 0.5;  // Electric Field in the drift region in kV/cm
   double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
   double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
   double Temperature = 87.6;  // LAr Temperature in K

   double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature); //drift velocity in the drift 
                                                                            //region (cm/us)
   double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature); //drift velocity between shield 
                                                                            //and induction (cm/us)
   double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature); //drift velocity between induction 
                                                                            //and collection (cm/us)
   double timepitch = driftvelocity*timetick;                               //time sample (cm) 
   double tSI = plane_pitch/driftvelocity_SI/timetick;                      //drift time between Shield and 
                                                                            //Collection planes (time samples)
   double tIC = plane_pitch/driftvelocity_IC/timetick;                      //drift time between Induction and 
                                                                            //Collection planes (time samples)


   // get input Cluster object(s).
   art::Handle< std::vector<recob::Cluster> > clusterListHandle;
   evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  
   // Declare some vectors..
   // Induction
   std::vector<double> Iwirefirsts;       // in cm
   std::vector<double> Iwirelasts;        // in cm
   std::vector<double> Itimefirsts;       // in cm
   std::vector<double> Itimelasts;        // in cm
   std::vector<double> Itimefirsts_line;  // in cm
   std::vector<double> Itimelasts_line;   // in cm
   std::vector < std::vector<art::Ptr<recob::Hit> > > IclusHitlists;
   std::vector<unsigned int> Icluster_count; 

   // Collection
   std::vector<double> Cwirefirsts;       // in cm
   std::vector<double> Cwirelasts;        // in cm
   std::vector<double> Ctimefirsts;       // in cm
   std::vector<double> Ctimelasts;        // in cm
   std::vector<double> Ctimefirsts_line;  // in cm
   std::vector<double> Ctimelasts_line;   // in cm
   std::vector< std::vector< art::Ptr<recob::Hit> > > CclusHitlists;
   std::vector<unsigned int> Ccluster_count; 

   // Some variables for the hit
   float time;            //hit time at maximum
   unsigned int wire;     //hit wire number
   unsigned int plane;    //hit plane number

   size_t startSPIndex = spacepoints->size(); //index for knowing which spacepoints are with which cluster
   size_t endSPIndex = spacepoints->size(); //index for knowing which spacepoints are with which cluster

   art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);

   for(size_t ii = 0; ii < clusterListHandle->size(); ++ii){

     art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
     
     /////////////////////////
     //////////// 2D track FIT
     /////////////////////////
     
     // Gaaaaaah! Change me soon!!! But, for now, 
     // let's just chuck one plane's worth of info. EC, 30-Mar-2011.
     ///\todo: This is very bad practice and should be changed ASAP
     if (cl->View() == geo::kZ) continue;      
     
     std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(ii);

     if(hitlist.size() == 1) continue;//only one Hit in this Cluster...will cause TGraph fit to fail.

     // sort the hit list to be sure it is in the correct order 
     // using the Hit < operator
     std::sort(hitlist.begin(), hitlist.end());
     
     TGraph *the2Dtrack = new TGraph(hitlist.size());
     
     std::vector<double> wires;
     std::vector<double> times;
     
      int np = 0;
      //loop over cluster hits
      for(art::PtrVector<recob::Hit>::const_iterator theHit = hitlist.begin(); theHit != hitlist.end();  theHit++){
	//recover the Hit
	//      recob::Hit* theHit = (recob::Hit*)(*hitIter);
	time = (*theHit)->PeakTime() ;
	
	time -= presamplings;
	
	plane = (*theHit)->WireID().Plane;
	wire = (*theHit)->WireID().Wire;
	
	//correct for the distance between wire planes
	if(plane == 1) time -= tIC;   // Collection
	
	//transform hit wire and time into cm
	double wire_cm;
         if(plane == 0)
	   wire_cm = (double)((wire+3.95) * wire_pitch);
         else
	   wire_cm = (double)((wire+1.84) * wire_pitch);
         
         double time_cm;
         if(time > tSI) time_cm = (double)( (time-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
         else time_cm = time*driftvelocity_SI*timetick;
         
         wires.push_back(wire_cm);
         times.push_back(time_cm);
	
         the2Dtrack->SetPoint(np,wire_cm,time_cm);
         np++;
      }//end of loop over cluster hits
      
      // fit the 2Dtrack and get some info to store
      try{
	the2Dtrack->Fit("pol1","Q");
      }
      catch(...){
	mf::LogWarning("Track3Dreco") << "The 2D track fit failed";
	continue;
      }
      
      TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
      double par[2];
      pol1->GetParameters(par);
      double intercept = par[0];
      double slope = par[1];
      
      double w0 = wires.front();               // first hit wire (cm)
      double w1 = wires.back();       	       // last hit wire (cm)
      double t0 = times.front();      	       // first hit time (cm)
      double t1 = times.back();       	       // last hit time (cm)
      double t0_line = intercept + (w0)*slope; // time coordinate at wire w0 on the fit line (cm)  
      double t1_line = intercept + (w1)*slope; // time coordinate at wire w1 on the fit line (cm)

      // actually store the 2Dtrack info
      switch(plane){
         case 0:
            Iwirefirsts.push_back(w0);
            Iwirelasts.push_back(w1);
            Itimefirsts.push_back(t0);
            Itimelasts.push_back(t1); 
            Itimefirsts_line.push_back(t0_line);
            Itimelasts_line.push_back(t1_line);    
            IclusHitlists.push_back(hitlist);
            Icluster_count.push_back(ii);
            break;
         case 1:
            Cwirefirsts.push_back(w0);
            Cwirelasts.push_back(w1);
            Ctimefirsts.push_back(t0);
            Ctimelasts.push_back(t1);
            Ctimefirsts_line.push_back(t0_line);
            Ctimelasts_line.push_back(t1_line);
            CclusHitlists.push_back(hitlist);
            Ccluster_count.push_back(ii);
            break;   
      }
      //delete the2Dtrack;
      delete pol1;
   }// end of loop over all input clusters
  
   /////////////////////////////////////////////////////
   /////// 2D Track Matching and 3D Track Reconstruction
   /////////////////////////////////////////////////////

   //loop over Collection view 2D tracks
   for(size_t collectionIter = 0; collectionIter < CclusHitlists.size(); ++collectionIter){  
     // Recover previously stored info
     double Cw0 = Cwirefirsts[collectionIter];
     double Cw1 = Cwirelasts[collectionIter];
     //double Ct0 = Ctimefirsts[collectionIter];
     //double Ct1 = Ctimelasts[collectionIter];
     double Ct0_line = Ctimefirsts_line[collectionIter];
     double Ct1_line = Ctimelasts_line[collectionIter];
     std::vector< art::Ptr<recob::Hit> > hitsCtrk = CclusHitlists[collectionIter];
     
     double collLength = TMath::Sqrt( TMath::Power(Ct1_line - Ct0_line,2) + TMath::Power(Cw1 - Cw0,2));
     
     //loop over Induction view 2D tracks
     for(size_t inductionIter = 0; inductionIter < IclusHitlists.size(); ++inductionIter){   
       // Recover previously stored info
       double Iw0 = Iwirefirsts[inductionIter];
       double Iw1 = Iwirelasts[inductionIter];
       //double It0 = Itimefirsts[inductionIter];
       //double It1 = Itimelasts[inductionIter];
       double It0_line = Itimefirsts_line[inductionIter];
       double It1_line = Itimelasts_line[inductionIter];
       std::vector< art::Ptr<recob::Hit> > hitsItrk = IclusHitlists[inductionIter];
       
       double indLength = TMath::Sqrt( TMath::Power(It1_line - It0_line,2) + TMath::Power(Iw1 - Iw0,2)); 
       
       bool forward_match = ((std::abs(Ct0_line-It0_line)<ftmatch*timepitch) && 
			     (std::abs(Ct1_line-It1_line)<ftmatch*timepitch));
       bool backward_match = ((std::abs(Ct0_line-It1_line)<ftmatch*timepitch) && 
			      (std::abs(Ct1_line-It0_line)<ftmatch*timepitch));
       
         

       // match 2D tracks
       if(forward_match || backward_match ){ 	
	 
	 // Reconstruct the 3D track
	 TVector3 XYZ0, XYZ1;  // track endpoints
	 if(forward_match){
	   XYZ0.SetXYZ(Ct0_line,(Cw0-Iw0)/(2.*TMath::Sin(Angle)),
		       (Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
	   XYZ1.SetXYZ(Ct1_line,(Cw1-Iw1)/(2.*TMath::Sin(Angle)),
		       (Cw1+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
	 }
	 else{
	   XYZ0.SetXYZ(Ct0_line,(Cw0-Iw1)/(2.*TMath::Sin(Angle)),
		       (Cw0+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
	   XYZ1.SetXYZ(Ct1_line,(Cw1-Iw0)/(2.*TMath::Sin(Angle)),
		       (Cw1+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
	 }

	 //compute track direction in Local co-ordinate system
	 //WARNING:  There is an ambiguity introduced here for the case of backwards-going tracks.  
	 //If available, vertex info. could sort this out.
	 TVector3 startpointVec,endpointVec;
	 TVector2 collVtx, indVtx;
	 if(XYZ0.Z() <= XYZ1.Z()){
	   startpointVec.SetXYZ(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
	   endpointVec.SetXYZ(XYZ1.X(),XYZ1.Y(),XYZ1.Z());
	   if(forward_match){
	     collVtx.Set(Ct0_line,Cw0);
	     indVtx.Set(It0_line,Iw0);
	   }
	   else{
	     collVtx.Set(Ct0_line,Cw0);
	     indVtx.Set(It1_line,Iw1);
	   }
	 }
	 else{
	   startpointVec.SetXYZ(XYZ1.X(),XYZ1.Y(),XYZ1.Z());
	   endpointVec.SetXYZ(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
	   if(forward_match){
	     collVtx.Set(Ct1_line,Cw1);
	     indVtx.Set(It1_line,Iw1);
	   }
	   else{
	     collVtx.Set(Ct1_line,Cw1);
	     indVtx.Set(It0_line,Iw0);
	   }
	 }
         
	 //compute track (normalized) cosine directions in the TPC co-ordinate system
	 TVector3 DirCos = endpointVec - startpointVec;
	 
	 //SetMag casues a crash if the magnitude of the vector is zero
	 try{
	   DirCos.SetMag(1.0);//normalize vector
	 }
	 catch(...){
	   mf::LogWarning("Track3Dreco") <<"The Spacepoint is infinitely small";
	   continue;
	 }
	 
	 art::Ptr <recob::Cluster> cl1(clusterListHandle,Icluster_count[inductionIter]);
	 art::Ptr <recob::Cluster> cl2(clusterListHandle,Ccluster_count[collectionIter]);
	 art::PtrVector<recob::Cluster> clustersPerTrack;
	 clustersPerTrack.push_back(cl1);
	 clustersPerTrack.push_back(cl2);
	 
	 
	 /////////////////////////////
	 // Match hits
	 ////////////////////////////
	        
	 std::vector< art::Ptr<recob::Hit> > minhits = hitsCtrk.size() <= hitsItrk.size() ? hitsCtrk : hitsItrk;
	 std::vector< art::Ptr<recob::Hit> > maxhits = hitsItrk.size() <= hitsCtrk.size() ? hitsCtrk : hitsItrk;
	 
	 
	 std::vector<bool> maxhitsMatch(maxhits.size());
	 for(size_t it = 0; it < maxhits.size(); ++it) maxhitsMatch[it] = false;
	 
	 std::vector<recob::Hit*> hits3Dmatched;
	 // For the matching start from the view where the track projection presents less hits
	 unsigned int imaximum = 0;
	 //loop over hits

	 startSPIndex = spacepoints->size();

	 for(size_t imin = 0; imin < minhits.size(); ++imin){ 
	   //get wire - time coordinate of the hit
	   unsigned int wire,plane1,plane2;
	   wire = minhits[imin]->WireID().Wire;
	   plane1 = minhits[imin]->WireID().Plane;

	   // get the wire-time co-ordinates of the hit to be matched
	   double w1;
	   if(plane1 == 0)
	     w1 = (double)((wire+3.95)*wire_pitch);
	   else
	     w1 = (double)((wire+1.84) * wire_pitch);
	   double temptime1 = minhits[imin]->PeakTime()-presamplings;
	   if(plane1 == 1) temptime1 -= tIC;
	   double t1;
	   if(temptime1>tSI) t1 = (double)( (temptime1-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
	   else t1 = temptime1*driftvelocity_SI*timetick;
	   
	   //get the track origin co-ordinates in the two views
	   TVector2 minVtx2D;
	   plane1==1 ? minVtx2D.Set(collVtx.X(),collVtx.Y()): minVtx2D.Set(indVtx.X(),indVtx.Y());
	   TVector2 maxVtx2D;
	   plane1==1 ? maxVtx2D.Set(indVtx.X(),indVtx.Y()): maxVtx2D.Set(collVtx.X(),collVtx.Y());
	   
	   double ratio = (collLength>indLength) ? collLength/indLength : indLength/collLength;
           
	   //compute the distance of the hit (imin) from the relative track origin
	   double minDistance = ratio*TMath::Sqrt(TMath::Power(t1-minVtx2D.X(),2) 
						  + TMath::Power(w1-minVtx2D.Y(),2));
	   
	   //core matching algorithm
	   double difference = 9999999.;	
  
	   //loop over hits of the other view	   
	   for(size_t imax = 0; imax < maxhits.size(); ++imax){ 
	     if(!maxhitsMatch[imax]){
	       //get wire - time coordinate of the hit
	       wire = maxhits[imax]->WireID().Wire;
	       plane2 = maxhits[imax]->WireID().Plane;
	       
	       double w2;
	       if(plane2 == 0)
		 w2 = (double)((wire+3.95)*wire_pitch);
	       else
		 w2 = (double)((wire+1.84)*wire_pitch);
	       double temptime2 = maxhits[imax]->PeakTime()-presamplings;
	       if(plane2 == 1) temptime2 -= tIC;
	       double t2;
	       if(temptime2 > tSI) t2 = (double)( (temptime2-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
	       else t2 = temptime2*driftvelocity_SI*timetick;
	      
	       
	       bool timematch = (std::abs(t1-t2)<ftmatch*timepitch);
	       bool wirematch = (std::abs(w1-w2)<wireShift*wire_pitch);
	       
	       double maxDistance = TMath::Sqrt(TMath::Power(t2-maxVtx2D.X(),2)+TMath::Power(w2-maxVtx2D.Y(),2));
	       if (wirematch && timematch && std::abs(maxDistance-minDistance)<difference) {
		 difference = std::abs(maxDistance-minDistance);
		 imaximum = imax;
	       }
	     }
	   }// end loop over max hits
	   maxhitsMatch[imaximum]=true;
	  
	   art::PtrVector<recob::Hit> sp_hits;
	   if(difference!= 9999999.){
	     sp_hits.push_back(minhits[imin]);
	     sp_hits.push_back(maxhits[imaximum]);
	   }
	   
	   // Get the time-wire co-ordinates of the matched hit
	   wire = maxhits[imaximum]->WireID().Wire;
	   plane2 = maxhits[imaximum]->WireID().Plane;
	   
	   double w1_match;
	   if(plane2 == 0)
	     w1_match = (double)((wire+3.95)*wire_pitch);
	   else
	     w1_match = (double)((wire+1.84)*wire_pitch);
	   double temptime3 = maxhits[imaximum]->PeakTime()-presamplings;
	   if(plane2 == 1) temptime3 -= tIC;
	   double t1_match;
	   if(temptime3 > tSI) t1_match = (double)( (temptime3-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
	   else t1_match = temptime3*driftvelocity_SI*timetick;
	   
	   
	   // create the 3D hit, compute its co-ordinates and add it to the 3D hits list	  
	   double Ct = plane1==1?t1:t1_match;
	   double Cw = plane1==1?w1:w1_match;
	   double Iw = plane1==1?w1_match:w1;
	   
	   const TVector3 hit3d(Ct,(Cw-Iw)/(2.*TMath::Sin(Angle)),(Cw+Iw)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle)); 
	   Double_t hitcoord[3];       
	   hitcoord[0] = hit3d.X();
	   hitcoord[1] = hit3d.Y();
	   hitcoord[2] = hit3d.Z();           

	   Double_t hitcoord_errs[3];
	   for (int i=0; i<3; i++) hitcoord_errs[i]=-1.000;

	   //3d point at end of track
	   recob::SpacePoint mysp(hitcoord, hitcoord_errs, -1., spacepoints->size());

	   spacepoints->push_back(mysp);

	   // associate the hits to the space point
	   util::CreateAssn(*this, evt, *spacepoints, sp_hits, *shassn);
	   
	 }//loop over min-hits
      
	 endSPIndex = spacepoints->size();

	 // Add the 3D track to the vector of the reconstructed tracks
	 if(spacepoints->size() > startSPIndex || clustersPerTrack.size()>0){
	   
	   std::vector<TVector3>               xyz;
	   xyz.push_back(startpointVec);
	   xyz.push_back(endpointVec);
	   std::vector<TVector3>               dir_xyz;
	   dir_xyz.push_back(DirCos);
	   dir_xyz.push_back(DirCos);
	   std::vector< std::vector <double> > dQdx = std::vector< std::vector<double> >(0);
	   std::vector<double>                 fitMomentum = std::vector<double>(2, util::kBogusD);
	  
	   recob::Track  the3DTrack(xyz,dir_xyz,dQdx, fitMomentum,tcol->size());
	   tcol->push_back(the3DTrack);

	   // associate the track with its spacepoints
	   util::CreateAssn(*this, evt, *tcol, *spacepoints, *sassn, startSPIndex, endSPIndex);

	   // associate the track with its clusters
	   util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *cassn);
	   
	   art::FindManyP<recob::Hit> fmhc(clustersPerTrack, evt, fClusterModuleLabel);

	   // get the hits associated with each cluster and associate those with the track
	   for(size_t p = 0; p < clustersPerTrack.size(); ++p){
	     std::vector< art::Ptr<recob::Hit> > hits = fmhc.at(p);
	     util::CreateAssn(*this, evt, *tcol, hits, *hassn);
	   }
	   
	 }
       
       } //close match 2D tracks
       
       
     }//close loop over Induction view 2D tracks
    
   }//close loop over Collection view 2D tracks

   mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
   mf::LogVerbatim("Summary") << "Track3Dreco Summary:";
   for(unsigned int i = 0; i<tcol->size(); ++i){
      mf::LogVerbatim("Summary") << tcol->at(i) ;
   } 
 
   evt.put(std::move(tcol));
   evt.put(std::move(spacepoints));
   evt.put(std::move(cassn));
   evt.put(std::move(sassn));
   evt.put(std::move(hassn));
   evt.put(std::move(shassn));

   return;
}

  DEFINE_ART_MODULE(Track3Dreco)

} // namespace
