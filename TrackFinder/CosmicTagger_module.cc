////////////////////////////////////////////////////////////////////////
// Class:       CosmicTagger
// Module Type: producer
// File:        CosmicTagger_module.cc
//
// Generated at Tue Apr  9 10:24:15 2013 by Sarah Lockwitz using artmod
// from art v1_02_06.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicTagger
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "Geometry/Geometry.h"
#include "Geometry/geo.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
//\todo Reconstruction Producers should never include Simulation headers
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
//#include "TTree.h"
//#include "TH1.h"

namespace simb{
  class MCTruth;
  class MCParticle;
}

namespace trkf {
  class CosmicTagger;
  class SpacePoint;
  class Track;
}

class trkf::CosmicTagger : public art::EDProducer {
public:
  explicit CosmicTagger(fhicl::ParameterSet const & p);
  virtual ~CosmicTagger();

  void produce(art::Event & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void endJob() override;

private:

  // Declare member data here.

  SpacePointAlg fSptalg;

  std::vector<recob::SpacePoint> spts;
  std::vector<float> all_x, all_y, all_z; 

  int   fReadOutWindowSize;
  float fSamplingRate;

  // stuff to set in the fcl file
  float fTotalBoundaryLimit; // 15
  float f3DSpillDistance;    // 12
  int   fSpillVetoCtr;       // 2
  int   fdTLimit;            // 8
  int   fdWLimit;            // 8
};


trkf::CosmicTagger::CosmicTagger(fhicl::ParameterSet const & p)
// :
 : fSptalg( p.get< fhicl::ParameterSet >("SpacePointAlg") )
// Initialize member data here.
{

  this->reconfigure(p);

  // Call appropriate Produces<>() functions here.
  produces< std::vector<recob::Track> >();
  produces< std::vector<recob::Hit> >();
  produces< std::vector<recob::SpacePoint> >();
  produces< art::Assns<recob::Track, recob::Hit> >();
  produces< art::Assns<recob::Track, recob::SpacePoint> >();
  produces< art::Assns<recob::Track, recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
  produces< std::vector<recob::Cluster> >();

}

trkf::CosmicTagger::~CosmicTagger()
{
  // Clean up dynamic memory and other resources here.
}

void trkf::CosmicTagger::produce(art::Event & e)
{
  // Implementation of required member function here.




  std::unique_ptr< std::vector<recob::Track > >      outTracks     ( new std::vector<recob::Track>);
  std::unique_ptr< std::vector<recob::Hit > >        outHits       ( new std::vector<recob::Hit>);
  std::unique_ptr< std::vector<recob::SpacePoint > > outSpacePoints( new std::vector<recob::SpacePoint>);
  std::unique_ptr< std::vector<recob::Cluster> >     outClusters   ( new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Track, recob::Hit > >        assnOutTracksHits   ( new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr< art::Assns<recob::Track, recob::SpacePoint > > assnOutTracksSpts   ( new art::Assns<recob::Track, recob::SpacePoint>);
  std::unique_ptr< art::Assns<recob::Track, recob::Cluster > >    assnOutTracksCluster( new art::Assns<recob::Track, recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Cluster, recob::Hit > >      assnOutClustersHits ( new art::Assns<recob::Cluster, recob::Hit>);

  std::vector<recob::SpacePoint> badSpts;
  std::vector<int> badSptsIDs;
  std::vector<float> vx1, vx2,vy1,vy2,vz1,vz2;


  std::string trackModuleLabel          = "track";
  std::string trackModuleLabelPreSpill  = trackModuleLabel;
  std::string trackModuleLabelPostSpill = trackModuleLabel;
  trackModuleLabelPreSpill             +="PreSpill";
  trackModuleLabelPostSpill            +="PostSpill";

  art::Handle<std::vector<recob::Track> > Trk_h;
  art::Handle<std::vector<recob::Track> > TrkPre_h;
  art::Handle<std::vector<recob::Track> > TrkPost_h;
  e.getByLabel(trackModuleLabel         ,Trk_h    );
  e.getByLabel(trackModuleLabelPreSpill ,TrkPre_h );
  e.getByLabel(trackModuleLabelPostSpill,TrkPost_h);
  art::PtrVector<recob::Track> TrkVec;
  art::PtrVector<recob::Track> TrkVecPre;
  art::PtrVector<recob::Track> TrkVecPost;

  std::vector<double> preSpillEndPtsX, preSpillEndPtsY, preSpillEndPtsZ;
  std::vector<double> postSpillEndPtsX, postSpillEndPtsY, postSpillEndPtsZ;
  std::vector<double> preSpillAngles, postSpillAngles;
  preSpillEndPtsX.clear(); preSpillEndPtsY.clear(); preSpillEndPtsZ.clear();
  postSpillEndPtsX.clear(); postSpillEndPtsY.clear(); postSpillEndPtsZ.clear();
  preSpillAngles.clear(); postSpillAngles.clear();


  for(unsigned int i=0; i < Trk_h->size(); i++) {
    art::Ptr<recob::Track> track(Trk_h,i);
    TrkVec.push_back(track);
  }
  for(unsigned int i=0; i < TrkPre_h->size(); i++) {
    art::Ptr<recob::Track> track(TrkPre_h,i);
    TrkVecPre.push_back(track);
  }
  for(unsigned int i=0; i < TrkPost_h->size(); i++) {
    art::Ptr<recob::Track> track(TrkPost_h,i);
    TrkVecPost.push_back(track);
  }


  //////////////////////////////////////////
  // Checking Pre- & Post-Spill Clusters
  //////////////////////////////////////////
  std::string clusterModuleLabel          = "cluster";
  std::string clusterModuleLabelPreSpill  = "clusterPreSpill";
  std::string clusterModuleLabelPostSpill = "clusterPostSpill";
  art::Handle<std::vector<recob::Cluster> > Cluster_h;
  art::Handle<std::vector<recob::Cluster> > ClusterPreSpill_h;
  art::Handle<std::vector<recob::Cluster> > ClusterPostSpill_h;
  e.getByLabel(clusterModuleLabel         , Cluster_h         );
  e.getByLabel(clusterModuleLabelPreSpill , ClusterPreSpill_h );
  e.getByLabel(clusterModuleLabelPostSpill, ClusterPostSpill_h);

  std::vector<double> dTdW_PreSpill[3], dQdW_PreSpill[3];
  std::vector<double> endPtsT_PreSpill[3], endPtsW_PreSpill[3];

  std::vector<double> dTdW_PostSpill[3], dQdW_PostSpill[3];
  std::vector<double> endPtsT_PostSpill[3], endPtsW_PostSpill[3];


  for( size_t i = 0; i < ClusterPreSpill_h->size(); i++ ) {
    art::Ptr<recob::Cluster> cc(ClusterPreSpill_h, i);
    int view = cc->View();
    
    dTdW_PreSpill[view].push_back( cc->dTdW() );
    // crappy work around:    
    //dTdW_PreSpill[view].push_back( (cc->EndPos()[1]-cc->StartPos()[1])/(cc->EndPos()[0]-cc->StartPos()[0]) ); 

    dQdW_PreSpill[view].push_back( cc->dQdW() );
    double tt = fabs(cc->StartPos()[1]-fReadOutWindowSize) < fabs(cc->EndPos()[1]-fReadOutWindowSize) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
    double ww = fabs(cc->StartPos()[1]-fReadOutWindowSize) < fabs(cc->EndPos()[1]-fReadOutWindowSize) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];

    mf::LogInfo("CosmicTagger Info") << "View: "<< view <<" FILLING PRE SPILL T'S & W'S.  t: " << tt << " w: " << ww 
				     << " startT:" << cc->StartPos()[1]
				     << " startW:" << cc->StartPos()[0]
				     << " endT:" << cc->EndPos()[1]
				     << " endW:" << cc->EndPos()[0] ;

    endPtsT_PreSpill[view].push_back( tt );
    endPtsW_PreSpill[view].push_back( ww );

  }

  for( size_t i = 0; i < ClusterPostSpill_h->size(); i++ ) {
    art::Ptr<recob::Cluster> cc(ClusterPostSpill_h, i);
    int view = cc->View();

    dTdW_PostSpill[view].push_back( cc->dTdW() );
    //dTdW_PostSpill[view].push_back( (cc->EndPos()[1]-cc->StartPos()[1])/(cc->EndPos()[0]-cc->StartPos()[0]) ); // crappy work around    

    dQdW_PostSpill[view].push_back( cc->dQdW() );
    double tt = fabs(cc->StartPos()[1]-2*fReadOutWindowSize) < fabs(cc->EndPos()[1]-2*fReadOutWindowSize) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
    double ww = fabs(cc->StartPos()[1]-2*fReadOutWindowSize) < fabs(cc->EndPos()[1]-2*fReadOutWindowSize) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];
    mf::LogInfo("CosmicTagger Info") << "View: "<< view << "FILLING POST SPILL T'S & W'S.  t: " 
				     << tt << " w: " << ww 
				     << " startT:" << cc->StartPos()[1]
				     << " startW:" << cc->StartPos()[0]
				     << " endT:" << cc->EndPos()[1]
				     << " endW:" << cc->EndPos()[0];
    endPtsT_PostSpill[view].push_back( tt );
    endPtsW_PostSpill[view].push_back( ww );

  }
  // End of Pre- & Post-Spill Clusters
  /////////////////////////////////////////////


  //////////////////////////////
  // Filling PreSpill End Pts 
  //////////////////////////////


  art::FindMany<recob::SpacePoint> sptsPre(TrkPre_h, e, trackModuleLabelPreSpill);


  for(unsigned int i=0; i<TrkPre_h->size(); i++) {

    std::vector<const recob::SpacePoint* > spts = sptsPre.at(i);

    if(spts.size()<=0) { continue;}
    std::vector<double> tempX, tempY, tempZ;

    for(unsigned int sp=0; sp<spts.size(); sp++) {      
      tempX.push_back(spts[sp]->XYZ()[0]);
      tempY.push_back(spts[sp]->XYZ()[1]);
      tempZ.push_back(spts[sp]->XYZ()[2]);
    }

    if(tempX.size()<2) continue;


    float maxDist=0;
    int endPtA=-1, endPtB=-1;
    for(unsigned int aa=0; aa<tempX.size()-1;aa++) {
      for(unsigned int bb=aa+1; bb<tempX.size(); bb++) {
	float dist = sqrt( pow(tempX[aa]-tempX[bb],2) + pow(tempY[aa]-tempY[bb],2) + pow(tempZ[aa]-tempZ[bb],2) );
	if(dist>maxDist) {maxDist=dist; endPtA=aa; endPtB=bb;}
      }
    }

    //    art::Ptr<recob::Track> tTrackPre = TrkVecPre.at(i);
    //std::cerr << tempX[endPtA] << " "<< tempY[endPtA] << " "<< tempZ[endPtA] << std::endl;
    //std::cerr << tempX[endPtB] << " "<< tempY[endPtB] << " "<< tempZ[endPtB] << std::endl;


    if( tempY[endPtA] < tempY[endPtB] ) { int tempXX=endPtB; endPtB=endPtA; endPtA=tempXX;}

    preSpillEndPtsX.push_back(tempX[endPtA]);
    preSpillEndPtsY.push_back(tempY[endPtA]);
    preSpillEndPtsZ.push_back(tempZ[endPtA]);
    
    preSpillEndPtsX.push_back(tempX[endPtB]);
    preSpillEndPtsY.push_back(tempY[endPtB]);
    preSpillEndPtsZ.push_back(tempZ[endPtB]);
    

    double dX = tempX[endPtA] - tempX[endPtB];
    double dY = tempY[endPtA] - tempY[endPtB];
    double dZ = tempZ[endPtA] - tempZ[endPtB];
    preSpillAngles.push_back(atan(sqrt(pow(dX,2)+pow(dZ,2))/dY));

    spts.clear();
  }

  /////////////////////////////////
  // Filling Post Spill End Pts 
  /////////////////////////////////
  art::FindManyP<recob::SpacePoint> sptsPost(TrkPost_h, e, trackModuleLabelPostSpill);
  for(unsigned int i=0; i<TrkPost_h->size(); i++) {

    std::vector<art::Ptr<recob::SpacePoint> > spts = sptsPost.at(i);

    if(spts.size()<=0) { continue;}
    std::vector<double> tempX, tempY, tempZ;
    for(unsigned int sp=0; sp<spts.size(); sp++) {

      tempX.push_back(spts[sp]->XYZ()[0]);
      tempY.push_back(spts[sp]->XYZ()[1]);
      tempZ.push_back(spts[sp]->XYZ()[2]);
    }

    if(tempX.size()<2) continue;

    float maxDist=0;
    int endPtA=-1, endPtB=-1;
    for(unsigned int aa=0; aa<tempX.size()-1;aa++) {
      for(unsigned int bb=aa+1; bb<tempX.size(); bb++) {
	float dist = sqrt( pow(tempX[aa]-tempX[bb],2) + pow(tempY[aa]-tempY[bb],2) + pow(tempZ[aa]-tempZ[bb],2) );
	if(dist>maxDist) {maxDist=dist; endPtA=aa; endPtB=bb;}
      }
    }

    //    art::Ptr<recob::Track> tTrackPost = TrkVecPost.at(i);

    if(tempY[endPtA]<tempY[endPtB]) { int tempXX=endPtB; endPtB=endPtA; endPtA=tempXX;}
    //    std::cerr << tempX[endPtA] << " "<< tempY[endPtA] << " "<< tempZ[endPtA] << std::endl;
    //    std::cerr << tempX[endPtB] << " "<< tempY[endPtB] << " "<< tempZ[endPtB] << std::endl;
    postSpillEndPtsX.push_back(tempX[endPtA]);
    postSpillEndPtsY.push_back(tempY[endPtA]);
    postSpillEndPtsZ.push_back(tempZ[endPtA]);
    
    postSpillEndPtsX.push_back(tempX[endPtB]);
    postSpillEndPtsY.push_back(tempY[endPtB]);
    postSpillEndPtsZ.push_back(tempZ[endPtB]);
    
    double dX = tempX[endPtA] - tempX[endPtB];
    double dY = tempY[endPtA] - tempY[endPtB];
    double dZ = tempZ[endPtA] - tempZ[endPtB];
    preSpillAngles.push_back(atan(sqrt(pow(dX,2)+pow(dZ,2))/dY));
    
    spts.clear();
  }



  /////////////////////////////////
  // LOOPING OVER INSPILL TRACKS
  /////////////////////////////////
  art::FindManyP<recob::SpacePoint> sptsSpill   (Trk_h, e, trackModuleLabel);
  art::FindManyP<recob::Hit>        hitsSpill   (Trk_h, e, trackModuleLabel);
  art::FindManyP<recob::Cluster>    ClusterSpill(Trk_h, e, trackModuleLabel);

  std::vector<int> badIDs;
  std::vector<int> badClusterIDs;
  badIDs.clear();
  badClusterIDs.clear();


  for( unsigned int iTrack=0; iTrack<Trk_h->size(); iTrack++ ) {
    art::ServiceHandle<geo::Geometry> geo;
    int origin      = -1;
    int isCosmic    =  0;

    int face1 = -1; 
    int face2 = -1; 


    art::Ptr<recob::Track> tTrack = TrkVec.at(iTrack);

    spts.clear();
    std::vector<art::Ptr<recob::SpacePoint> > spts   = sptsSpill.at(iTrack);
    std::vector<art::Ptr<recob::Hit> >       HitVec  = hitsSpill.at(iTrack);
    std::vector< art::Ptr< recob::Cluster> > ClusterVect = ClusterSpill.at(iTrack);
    art::PtrVector<recob::Cluster> ptrvs;

    if(spts.size()<3) continue;

    /////////////////////
    // cluster check
    std::vector<double> tdDTDW_post(3,-999), tdDQDW_post(3,-999), tdT_post(3,-999), tdW_post(3,-999);
    std::vector<double> tdDTDW_pre(3,-999), tdDQDW_pre(3,-999), tdT_pre(3,-999), tdW_pre(3,-999);
    std::vector<double> tdDist_window_post(3,-999), tdDist_window_pre(3,-999);
    std::vector<double> tDTDW(3,-999), tDTDW_pre(3,-999), tDTDW_post(3,-999);



    int preSpillVetoCtr=0;
    int postSpillVetoCtr=0;

    for( unsigned int k = 0; k < ClusterVect.size(); k++ ) {

      ptrvs.push_back( ClusterVect[k] );

      int clID = ClusterVect[k]->ID();
      int view = ClusterVect[k]->View();
      double dtdw = ClusterVect[k]->dTdW();
      double dqdw = ClusterVect[k]->dQdW();


      double t0 = ClusterVect[k]->StartPos()[1] < ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[1] : ClusterVect[k]->EndPos()[1];
      double t1 = ClusterVect[k]->StartPos()[1] > ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[1] : ClusterVect[k]->EndPos()[1];
      double w0 = ClusterVect[k]->StartPos()[1] < ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[0] : ClusterVect[k]->EndPos()[0];
      double w1 = ClusterVect[k]->StartPos()[1] > ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[0] : ClusterVect[k]->EndPos()[0];


      dtdw = ClusterVect[k]->dTdW(); 

      tDTDW[view]=dtdw;
      mf::LogInfo("CosmicTagger") << "Track number:" 
				  << iTrack << " Looking for " << clID<< ", T:"
				  << ClusterVect[k]->StartPos()[1]
				  <<","<< ClusterVect[k]->EndPos()[1] 
				  << " t0:"<< t0 << " t1:" << t1 <<", W:"
				  << ClusterVect[k]->StartPos()[0]<<","
				  <<ClusterVect[k]->EndPos()[0] <<", view:"<< view << ", dtdw:"<< dtdw;


      // JUST LOOP OVER PRE & POST SPILL CLUSTERS HERE
      int timeLimit=1000;

      
      // check against prespill window
      if( fabs(t0-fReadOutWindowSize) < timeLimit ) { // is it near the prespill boundary?
	int indx=-1;
	double dist=9999;
	for(unsigned int j=0; j<endPtsW_PreSpill[view].size(); j++) {
	  double thisDist = sqrt( pow(t0-endPtsT_PreSpill[view].at(j) ,2) + pow(w0-endPtsW_PreSpill[view].at(j) ,2));
	  if(thisDist<dist) { indx=j; dist=thisDist;}
	}
	if(indx<0) continue;


	tDTDW_pre[view] = dTdW_PreSpill[view] .at(indx);

	tdDTDW_pre[view]        = fabs( dtdw - dTdW_PreSpill[view] .at(indx) );
	tdDQDW_pre[view]        = fabs( dqdw - dQdW_PreSpill[view] .at(indx) );
	tdT_pre   [view]        = t0 - endPtsT_PreSpill[view].at(indx) ;
	tdW_pre   [view]        = fabs( w0 - endPtsW_PreSpill[view].at(indx) );
	tdDist_window_pre[view] = dist;


	mf::LogInfo("CosmicTagger") << "closest prespill:  dtdw:" << dTdW_PreSpill[view] .at(indx) << " T:" << endPtsT_PreSpill[view].at(indx) 
				    << " W:" << endPtsW_PreSpill[view].at(indx)
				    << "  Saving dist_window_pre: " << tdDist_window_pre[view] 
				    << " dT_pre:"<< tdT_pre[view] 
				    << " dW_pre:" << tdW_pre[view];

	if(tdT_pre[view]<=fdTLimit && tdW_pre[view]<=fdWLimit) preSpillVetoCtr++;

      } // end of prespill check

      if( fabs(t1-2*fReadOutWindowSize) < timeLimit ) { // is it near the postspill boundary?
	int indx=-1;
	double dist=9999;
	for(unsigned int j=0; j<endPtsW_PostSpill[view].size(); j++) {
	  double thisDist = sqrt( pow(t1-endPtsT_PostSpill[view].at(j) ,2) + pow(w1-endPtsW_PostSpill[view].at(j) ,2));
	  if(thisDist<dist) { indx=j; dist=thisDist;}
	}
	if(indx<0) continue;

	mf::LogInfo("CosmicTagger") << "Post spill closest postspill:  dtdw:" << dTdW_PostSpill[view] .at(indx) 
				    << " T:" << endPtsT_PostSpill[view].at(indx) 
				    << " W:" << endPtsW_PostSpill[view].at(indx);

	tDTDW_post[view] = dTdW_PostSpill[view] .at(indx);

	tdDTDW_post[view]        = fabs( dtdw - dTdW_PostSpill[view] .at(indx) );
	tdDQDW_post[view]        = fabs( dqdw - dQdW_PostSpill[view] .at(indx) );
	tdT_post   [view]        = endPtsT_PostSpill[view].at(indx) - t1 ;
	tdW_post   [view]        = fabs( w1 - endPtsW_PostSpill[view].at(indx) );
	tdDist_window_post[view] = dist;

	if(tdT_post[view]<=fdTLimit && tdW_post[view]<=fdWLimit) postSpillVetoCtr++;

      } // end of postspill check
     
    } // End of Cluster Check
    /////////////////////////////////



    all_x.clear();
    all_y.clear();
    all_z.clear();



    // The two methods below yield slightly different results
    // geo->Det... seems to have a smaller good region than PositionToTPC

//    for(unsigned int j=0; j<spts.size(); j++) {
//	
//      if( spts[j]->XYZ()[0]> 2.0*geo->DetHalfWidth() || spts[j]->XYZ()[0]<0 ) {
//	std::cerr << "FAILING Z!!  " << spts[j]->XYZ()[0] << std::endl;
//	continue;
//      }
//      if( fabs( spts[j]->XYZ()[1] ) > geo->DetHalfHeight() ) {
//	std::cerr << "FAILING Y!!  " << spts[j]->XYZ()[1] << std::endl;
//	continue; 
//      }
//      if( spts[j]->XYZ()[2] < 0 || spts[j]->XYZ()[2] > geo->DetLength() ) {
//	std::cerr << "FAILING Z!!  " << spts[j]->XYZ()[2] << std::endl;
//	continue; 
//      }
//
//      all_x.push_back( spts[j]->XYZ()[0] );
//      all_y.push_back( spts[j]->XYZ()[1] );
//      all_z.push_back( spts[j]->XYZ()[2] );
//    }



    double testXYZ[3]={0.};
    unsigned int testTPC=0, testCRYO=0;
    for(unsigned int j=0; j<spts.size(); j++) {
      try{
	testXYZ[0]=spts[j]->XYZ()[0];
	testXYZ[1]=spts[j]->XYZ()[1];
	testXYZ[2]=spts[j]->XYZ()[2];
	geo->PositionToTPC(testXYZ, testTPC, testCRYO);
      }
      catch(cet::exception &e) {
	continue;
      }
      all_x.push_back( spts[j]->XYZ()[0] );
      all_y.push_back( spts[j]->XYZ()[1] );
      all_z.push_back( spts[j]->XYZ()[2] );
    }




    // TRUTH CHECK
    art::ServiceHandle<cheat::BackTracker> bt;
    const art::Ptr<simb::MCTruth> mcpart = bt->TrackIDToMCTruth( tTrack->ID() );
    origin = mcpart->Origin();
    //    const simb::MCParticle *simpart = bt->TrackIDToParticle( tTrack->ID() );




    // This is a way to see how "wide" a track is
    TVectorD xV;
    float centroidX=0, centroidY=0, centroidZ=0;
    TMatrixD mCluster(all_x.size(), 3);

    centroidX = accumulate(all_x.begin(),all_x.end(),0);
    centroidY = accumulate(all_y.begin(),all_y.end(),0);
    centroidZ = accumulate(all_z.begin(),all_z.end(),0);
    centroidX /= 1.0*all_x.size();
    centroidY /= 1.0*all_y.size();
    centroidZ /= 1.0*all_z.size();


    for(unsigned int ii=0; ii<all_x.size(); ii++) {
      double vals[3];
      vals[0] = all_x[ii] - centroidX;
      vals[1] = all_y[ii] - centroidY;
      vals[2] = all_z[ii] - centroidZ;
      xV.Use(3, vals);
      TMatrixDRow(mCluster, ii) = xV;
    }


    TMatrixD Aw = mCluster;
    if(Aw.GetNrows()<3) continue;
    TDecompSVD svd(Aw);
    TMatrixD m_svdV = svd.GetV();
    TVectorD c_svd = svd.GetSig();



    double xyz[3], dxyz[3];//, xyzout[6];
    xyz[0] = centroidX;
    xyz[1] = centroidY;
    xyz[2] = centroidZ;
    dxyz[0] = m_svdV(0,0);
    dxyz[1] = m_svdV(1,0);
    dxyz[2] = m_svdV(2,0);






    double xyzout1[3], xyzout2[3];
    geo::ProjectToBoxEdge( xyz, dxyz, 0, 2.*geo->DetHalfWidth(), -1.*geo->DetHalfHeight(), geo->DetHalfHeight(), 0, geo->DetLength(), xyzout1 );
    double ndxyz[3]; 
    ndxyz[0] = -dxyz[0];
    ndxyz[1] = -dxyz[1];
    ndxyz[2] = -dxyz[2];
    geo::ProjectToBoxEdge( xyz, ndxyz, 0, 2.*geo->DetHalfWidth(), -1.*geo->DetHalfHeight(), geo->DetHalfHeight(), 0, geo->DetLength(), xyzout2 );

    if( xyzout1[0]<2 ) face1=2;
    else if( fabs( 2.*geo->DetHalfWidth()-xyzout1[0] ) < 2 ) face1=5;
    else if( fabs( geo->DetHalfHeight() - xyzout1[1] ) < 2 ) face1=3;
    else if( fabs( xyzout1[1] +geo->DetHalfHeight() ) < 2 )  face1=4;
    else if( fabs( xyzout1[2] ) < 2 )  face1=1;
    else if( fabs( geo->DetLength() - xyzout1[2] ) < 2 )  face1=6;

    if( xyzout2[0]<2 ) face2=2;
    else if( fabs( 2.*geo->DetHalfWidth()-xyzout2[0] ) < 2 ) face2=5;
    else if( fabs( geo->DetHalfHeight() - xyzout2[1] ) < 2 ) face2=3;
    else if( fabs( xyzout2[1] + geo->DetHalfHeight() ) < 2 ) face2=4;
    else if( fabs( xyzout2[2] ) < 2 ) face2=1;
    else if( fabs( geo->DetLength() - xyzout2[2] ) < 2 ) face2=6;

    if( face1<0 || face2<0 ) {
      mf::LogInfo("CosmicTagger")<< "Problem projecting to box edge.  Pt 1 Values: " <<
	" " << xyzout1[0]                                << 
	" " << fabs( 2.*geo->DetHalfWidth()-xyzout1[0] ) << 
	" " << fabs( geo->DetHalfHeight() - xyzout1[1] ) << 
	" " << fabs( xyzout1[1] + geo->DetHalfHeight() ) << 
	" " << fabs( xyzout1[2] )			       << 
	" " << fabs( geo->DetLength() - xyzout1[2] );
      
      mf::LogInfo("CosmicTagger")<< "Problem projecting to box edge. Pt 2 Values: " <<
	" "<< xyzout2[0]                                << 
	" "<< fabs( 2.*geo->DetHalfWidth()-xyzout2[0] ) << 
	" "<< fabs( geo->DetHalfHeight() - xyzout2[1] ) << 
	" "<< fabs( xyzout2[1] + geo->DetHalfHeight() ) << 
	" "<< fabs( xyzout2[2] )			      << 
	" "<< fabs( geo->DetLength() - xyzout2[2] );


    }




    // POOR MAN'S WAY OF FINDING END PTS :
    float pt2[3];
    pt2[0]=xyz[0]+10*dxyz[0];
    pt2[1]=xyz[1]+10*dxyz[1];
    pt2[2]=xyz[2]+10*dxyz[2];
    float aaa = sqrt( pow(xyz[0]-pt2[0],2) + pow(xyz[1]-pt2[1],2) + pow(xyz[2]-pt2[2],2) );

    //double yhi =  1.0*geo->DetHalfHeight();
    float maxVal=0;
    int endPt1=-1, endPt2=-1;
    for(unsigned int aa=0; aa+1<all_x.size(); aa++) {

      float bbb = sqrt( pow(xyz[0]-all_x[aa],2) + pow(xyz[1]-all_y[aa],2) + pow(xyz[2]-all_z[aa],2) );
      float ccc = sqrt( pow(pt2[0]-all_x[aa],2) + pow(pt2[1]-all_y[aa],2) + pow(pt2[2]-all_z[aa],2) );
      float aDist = sqrt( ccc*ccc - pow((ccc*ccc+aaa*aaa-bbb*bbb)/(2*aaa) ,2) );
      if(aDist>5) continue;

      for(unsigned int bb=aa+1; bb<all_x.size(); bb++) {

	bbb = sqrt( pow(xyz[0]-all_x[bb],2) + pow(xyz[1]-all_y[bb],2) + pow(xyz[2]-all_z[bb],2) );
	ccc = sqrt( pow(pt2[0]-all_x[bb],2) + pow(pt2[1]-all_y[bb],2) + pow(pt2[2]-all_z[bb],2) );
	aDist = sqrt( ccc*ccc - pow((ccc*ccc+aaa*aaa-bbb*bbb)/(2*aaa) ,2) );
	if(aDist>5) continue;

	float hereVal = sqrt( pow(all_x[aa]-all_x[bb],2) + pow(all_y[aa]-all_y[bb],2) + pow(all_z[aa]-all_z[bb],2) );
	if(hereVal>maxVal) {maxVal=hereVal; endPt1=aa; endPt2=bb;}
      }
    }

    if( all_y[endPt1]<all_y[endPt2] ) {
      float tempIndx = endPt1;
      endPt1 = endPt2;
      endPt2 = tempIndx;
    }







    //////////////////////////////////////////////////////////////////////
    //Let's check Cluster connections to pre & post spill planes first
    //////////////////////////////////////////////////////////////////////

    if( postSpillVetoCtr >= fSpillVetoCtr || preSpillVetoCtr >= fSpillVetoCtr ) isCosmic = 4;
    


    float preSpillDistance=1000, postSpillDistance=1000;
    if( isCosmic == 0 ) {
      for(unsigned int psEp=0; psEp<preSpillEndPtsY.size(); psEp++) {
	float dist =99;	
	float dist1 = sqrt( pow(all_x[endPt1]-(preSpillEndPtsX[psEp]-2*geo->DetHalfWidth()),2) + 
			    pow(all_y[endPt1]-preSpillEndPtsY[psEp],2) + 
			    pow(all_z[endPt1]-preSpillEndPtsZ[psEp],2) );
	float dist2 = sqrt( pow(all_x[endPt2]-(preSpillEndPtsX[psEp]-2*geo->DetHalfWidth()),2) + 
			    pow(all_y[endPt2]-preSpillEndPtsY[psEp],2) + 
			    pow(all_z[endPt2]-preSpillEndPtsZ[psEp],2) );
	dist = std::min(dist1,dist2);
	if( dist < preSpillDistance ) preSpillDistance = dist;
	if( dist < f3DSpillDistance ) isCosmic=2;
      }
      
      
      for(unsigned int psEp=0; psEp<postSpillEndPtsY.size(); psEp++) {
	float dist =99;
	float dist1 = sqrt( pow(all_x[endPt1]-2.0*geo->DetHalfWidth()-postSpillEndPtsX[psEp],2) + 
			    pow(all_y[endPt1]-postSpillEndPtsY[psEp],2) + 
			    pow(all_z[endPt1]-postSpillEndPtsZ[psEp],2) );
	float dist2 = sqrt( pow(all_x[endPt2]-2.0*geo->DetHalfWidth()-postSpillEndPtsX[psEp],2) + 
			    pow(all_y[endPt2]-postSpillEndPtsY[psEp],2) + 
			    pow(all_z[endPt2]-postSpillEndPtsZ[psEp],2) );
	dist = std::min(dist1,dist2);
	if( dist < postSpillDistance ) postSpillDistance = dist;
	if( dist < f3DSpillDistance ) isCosmic=2;
      }
    }


    
    
//    // then check y & z boundaries
//    int nBd=0;
//    float bndDist = 5;
//    if(fabs(all_y[endPt1] - geo->DetHalfHeight())<bndDist ) nBd++;
//    if(fabs(all_y[endPt2] + geo->DetHalfHeight())<bndDist ) nBd++;
//    if(fabs(all_z[endPt1] - geo->DetLength())<bndDist || fabs(all_z[endPt2] - geo->DetLength())<bndDist ) nBd++;
//    if(fabs(all_z[endPt1] )<bndDist || fabs(all_z[endPt2] )<bndDist ) nBd++;
//    if(isCosmic==0 && nBd>1) isCosmic=1;



    /////////////////////////////////////////
    // check modified crossed two boundaries
    /////////////////////////////////////////

    float totalBoundaryDistance=0;

    if(face1==1 || face2==1) {
      totalBoundaryDistance += std::min(all_z[endPt1], all_z[endPt2]);
    }
    if(face1==6 || face2==6 ) {
      totalBoundaryDistance += std::min( geo->DetLength() - all_z[endPt1], geo->DetLength() - all_z[endPt2] ) ;
    }
    if( face1==2 || face2==2 ) {
      totalBoundaryDistance += std::min( all_x[endPt1], all_x[endPt2]);
    }
    if( face1==5 || face2==5 ) {
      totalBoundaryDistance += std::min(2.*geo->DetHalfWidth() - all_x[endPt1], 2.*geo->DetHalfWidth() - all_x[endPt2]);
    }
    if( face1==3 || face2==3) {
      totalBoundaryDistance += std::min( geo->DetHalfHeight() - all_y[endPt1], geo->DetHalfHeight() - all_y[endPt2] );
    }
    if( face1==4 || face2==4 ) {
      totalBoundaryDistance += std::min( geo->DetHalfHeight() + all_y[endPt1], geo->DetHalfHeight() + all_y[endPt2] );
    }
    if ( isCosmic==0 && totalBoundaryDistance < fTotalBoundaryLimit ) isCosmic=5;

    mf::LogInfo("CosmicTagger Info") << "Boundaries: "<< face1 << ", " << face2 << ", totalBoundaryDistance   " << totalBoundaryDistance;

    
    
       
    mf::LogInfo("CosmicTagger Results") << "The IsCosmic value is "<< isCosmic << " origin: " << origin
					<< " | preSpillDistance: " << preSpillDistance
					<<  " postSpillDistance " << postSpillDistance << " | " 
					<< all_x[endPt1]<<","<< all_y[endPt1] << "," << all_z[endPt1]<< " | | " 
					<< all_x[endPt2]<< ","<< all_y[endPt2] <<"," << all_z[endPt2];


    ////////////////////////
    // Save Cosmic Stuff
    ////////////////////////    

    if( !isCosmic ) { //saving clusters that are not cosmics
      badIDs.push_back(tTrack->ID());
      
      for( unsigned int k = 0; k < ClusterVect.size(); k++ ) {
	badClusterIDs.push_back(ClusterVect[k]->ID() );
      }
    }
    else {
      outTracks->push_back( *tTrack );
      
      for (unsigned int p=0; p< HitVec.size(); p++) {
	outHits->push_back(*HitVec[p]);
      }
      
      util::CreateAssn(*this, e, *outTracks, HitVec, *assnOutTracksHits);
      util::CreateAssn(*this, e, *outTracks, spts, *assnOutTracksSpts);
      util::CreateAssn(*this, e, *outTracks, ptrvs, *assnOutTracksCluster );
    }
    
    

    HitVec.clear();



  }
  // END OF LOOPING OVER INSPILL TRACKS

  ///////////////////////////////////////////////////
  // looping over clusters
  art::PtrVector<recob::Cluster> tClusterVec;
  for(unsigned int i=0;i<Cluster_h->size(); i++) {
    art::Ptr<recob::Cluster> cltr(Cluster_h, i);
    tClusterVec.push_back( cltr );
  }


  art::FindManyP<recob::Hit> tHits   (Cluster_h, e, clusterModuleLabel);
  for (unsigned int i =0; i<tClusterVec.size(); i++) {
    art::Ptr<recob::Cluster> tCl = tClusterVec.at(i);
    std::vector< art::Ptr<recob::Hit> > tHitVec = tHits.at(i);

    int tempID = tCl->ID();
    //    if(tHitVec.size()<2) continue;
    if( find(badClusterIDs.begin(), badClusterIDs.end(), tempID) == badClusterIDs.end() ) {
      outClusters->push_back(*tCl);
      util::CreateAssn(*this, e, *(outClusters.get()), tHitVec, *(assnOutClustersHits.get()));
    }
  }
  // end looping over clusters
  ///////////////////////////////////////////////////



  e.put(std::move(outTracks));
  e.put(std::move(assnOutTracksCluster));
  e.put(std::move(assnOutTracksHits));
  e.put(std::move(assnOutTracksSpts));
  e.put(std::move(outClusters));
  e.put(std::move(outHits));
  e.put(std::move(outSpacePoints));
  e.put(std::move(assnOutClustersHits));



  TrkVec.clear();
  TrkVecPre.clear();
  TrkVecPost.clear();


} // end of produce

void trkf::CosmicTagger::beginJob()
{
  // Implementation of optional member function here.
}

void trkf::CosmicTagger::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  art::ServiceHandle<util::DetectorProperties> detp;
  fReadOutWindowSize = detp->ReadOutWindowSize();
  fSamplingRate = detp->SamplingRate();

  fTotalBoundaryLimit = p.get<float>("TPCBoundaryLimit", 15);
  f3DSpillDistance    = p.get<float>("SpillDistance",12); 

  fSpillVetoCtr = p.get<int>("SpillVetoCounter", 2);
  fdTLimit      = p.get<int>("dTLimit",8);
  fdWLimit      = p.get<int>("dWLimit",8);
}

void trkf::CosmicTagger::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(trkf::CosmicTagger)
