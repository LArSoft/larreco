////////////////////////////////////////////////////////////////////////
// Class:       CosmicTagger
// Module Type: producer
// File:        CosmicTagger_module.cc
//
// Generated at Mon Sep 24 18:21:00 2012 by Sarah Lockwitz using artmod
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CosmicTagger
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "Geometry/Geometry.h"
#include "Geometry/geo.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Track.h"

#include "AnalysisBase/CosmicTag.h"

//#include "SimulationBase/MCTruth.h"
//#include "SimulationBase/MCParticle.h"


#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TTree.h"
#include "TH1.h"




class TTree;
class TH1;

//namespace simb{
//  class MCTruth;
//  class MCParticle;
//}

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
  void fillOutOfSpillClusters();
  void doClusterCheck(  art::FindManyP<recob::Hit> tHits, art::FindManyP<recob::SpacePoint> tSpacePoints   );// override;
  void doTrackClusterCheck( std::vector< art::Ptr< recob::Cluster> > ClusterVect, std::vector<double> &t1Times, 
			    std::vector<double> &t2Times, std::vector<int> &fail );
  void doSomeSpacePointStuff(art::FindManyP<recob::SpacePoint> sptsSpill , int iTrack, int &fa, int &fb);
  //  void truthCheck(int trackID, int &origin, int &trackPdg, double &trackTime, double &endTime, float &trackP);

private:

  // Declare member data here.

  int   fReadOutWindowSize;
  float fSamplingRate;

  // stuff to set in the fcl file
  float fTotalBoundaryLimit; // 15
  float f3DSpillDistance;    // 12
  int   fSpillVetoCtr;       // 2
  int   fdTLimit;            // 8
  int   fdWLimit;            // 8
  std::string fTrackModuleLabel;
  std::string fClusterModuleLabel;
  int fDoTrackCheck;
  int fDoClusterCheck;
  int fClusterAssociatedToTracks;





  std::vector<double> tdDist_window_post, tdDist_window_pre;

  int failClusterCheck;




  art::Handle<std::vector<recob::Cluster> > Cluster_h;

  std::unique_ptr< std::vector<recob::Cluster> >     outClusters   ;

  std::unique_ptr< art::Assns<anab::CosmicTag, recob::Cluster > >    assnOutCosmicTagCluster;


  art::ServiceHandle<cheat::BackTracker> bt;


  std::unique_ptr< std::vector< anab::CosmicTag > > cosmicTagVector;
  std::unique_ptr< std::vector< anab::CosmicTag > > cosmicTagTrackVector;

};


int nEvent;

std::vector<float> all_x, all_y, all_z; 

TTree *tree;
typedef struct {
  int event;
  int trackNum;
  int org;
  int isCosmicVal;
  int clusterCheck;
  int pdg;
  float svd0;
  float svd1;
  float svd2;
  float x1;
  float x2;
  float x1gen;
  float x2gen;
  float y1; 
  float y2;
  float z1;
  float z2;
  float length;
  float yDistance;
  float preSpillDistance;
  float postSpillDistance;
  float time;
  float enterP;
  int face1;
  int face2;
  float totalBoundaryDistance;
  int type;


  std::vector<double> dist_window_post;
  std::vector<double> dist_window_pre;

  // would like to get a generator time in here
} cfTrack_t;
cfTrack_t cTrack;





trkf::CosmicTagger::CosmicTagger(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  this->reconfigure(p);

  // Call appropriate Produces<>() functions here.

  produces< std::vector<recob::Track> >();
  produces< std::vector<recob::Cluster> >();
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();


}

trkf::CosmicTagger::~CosmicTagger() {
  // Clean up dynamic memory and other resources here.
}

void trkf::CosmicTagger::produce(art::Event & e) {
  // Implementation of required member function here.




  std::unique_ptr< std::vector<recob::Track > >      outTracksForTags( new std::vector<recob::Track>);

  outClusters    = std::unique_ptr<std::vector<recob::Cluster> > ( new std::vector<recob::Cluster> );

  cosmicTagTrackVector = std::unique_ptr< std::vector< anab::CosmicTag > >( new std::vector<anab::CosmicTag> );

  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >    assnOutCosmicTagTrack( new art::Assns<recob::Track, anab::CosmicTag>);










  art::Handle<std::vector<recob::Track> > Trk_h;

  if(fDoTrackCheck) e.getByLabel(fTrackModuleLabel        ,Trk_h    );


  /*
art::Handle< std::vector<recob::Track> > trackh;
evt.getByLabel(fTrackModuleLabel, trackh);
std::vector<art::Ptr<recob::Track> > tracklist;
art::fill_ptr_vector(tracklist,trackh);


art::FindManyP<anab::CosmicTag> cosmictag( tracklist ,evt,fCosmicTagAssocLabel);
   */


  art::PtrVector<recob::Track> TrkVec;

  if(fDoTrackCheck) {
    for(unsigned int i=0; i < Trk_h->size(); i++) {
      art::Ptr<recob::Track> track(Trk_h,i);
      TrkVec.push_back(track);
    }
  }






  //////////////////////////////////////////////////////////////////////
  ////////////////

  

  failClusterCheck = 0;
  tdDist_window_post.clear();
  tdDist_window_pre .clear();

  for( int hh=0; hh<3; hh++ ) {
    tdDist_window_post.push_back(-999);
    tdDist_window_pre .push_back(-999);
  }


  if( fDoClusterCheck ) {
    e.getByLabel(fClusterModuleLabel        , Cluster_h         );
  }








  /////////////////////////////////
  // LOOPING OVER INSPILL TRACKS
  /////////////////////////////////
  

  // BEZIER TRACKING CURRENTLY DOES NOT SAVE THE SPACEPOINTS
  // CHANGE THIS BACK WHEN IT DOES
  //  if(fDoClusterCheck) art::FindManyP<recob::SpacePoint> sptsSpill   (Trk_h, e, fTrackModuleLabel);

  std::cout << "fDoTrackCheck " << fDoTrackCheck <<std::endl;//<< " number of tracks " << Trk_h->size() << std::endl;
  if( fDoTrackCheck ) {

    art::FindManyP<recob::Hit>        hitsSpill   (Trk_h, e, fTrackModuleLabel);

    for( unsigned int iTrack=0; iTrack<Trk_h->size(); iTrack++ ) {

      art::ServiceHandle<geo::Geometry> geo;
      int origin      = -1;
      int isCosmic    =  0;
      int trackPdg    = -999;
      int face1       = -1; 
      int face2       = -1; 
      int trackType   = -1;

      std::cerr << "WHAT TRACK ARE WE ON " << iTrack << std::endl;
      art::Ptr<recob::Track>              tTrack       = TrkVec.at(iTrack);

      std::vector<anab::CosmicTag> damnAssns;

      std::vector<art::Ptr<recob::Hit> >  HitVec  = hitsSpill.at(iTrack);

      std::vector< art::Ptr< recob::Cluster> > ClusterVect;

      if(fClusterAssociatedToTracks ) {
	art::FindManyP<recob::Cluster>    ClusterSpill(Trk_h, e, fTrackModuleLabel);
	ClusterVect = ClusterSpill.at(iTrack);
      }

      //spts.clear();
    
      std::vector<art::Ptr<recob::SpacePoint> > spts;

    


      art::PtrVector<recob::Hit> HitsToProcess;
      for( size_t i=0; i!=HitVec.size(); ++i) HitsToProcess.push_back(HitVec.at(i));
      art::PtrVector<recob::Cluster> ptrvs;




      ////////////////////////////////////////////////////
      // CHECKING TRUTH STUFF
      float trackP     = -999;
      double trackTime = -999;
      //double endTime   = -999;
      //      truthCheck( tTrack->ID(), origin, trackPdg, trackTime, endTime, trackP);


      ////////////////////////////////////////////////////
      // CHECKING SOME SPACEPOINT STUFF
      if( fDoClusterCheck && 0) {
	art::FindManyP<recob::SpacePoint> sptsSpill   (Trk_h, e, fTrackModuleLabel);
	doSomeSpacePointStuff(sptsSpill, iTrack, face1, face2);
      }

	



      // A BETTER WAY OF FINDING END POINTS:
      TVector3 tVector1 = tTrack->Vertex();
      TVector3 tVector2 = tTrack->End();


      float trackEndPt1_X = tVector1[0]; //all_x[endPt1];
      float trackEndPt1_Y = tVector1[1]; //all_y[endPt1];
      float trackEndPt1_Z = tVector1[2]; //all_z[endPt1];
      float trackEndPt2_X = tVector2[0]; //all_x[endPt2];
      float trackEndPt2_Y = tVector2[1]; //all_y[endPt2];
      float trackEndPt2_Z = tVector2[2]; //all_z[endPt2];

      if(trackEndPt1_X != trackEndPt1_X ) {
	std::cerr << "!!! FOUND A PROBLEM... the length is: " << tTrack->Length() << 
	  " np: " << tTrack->NumberTrajectoryPoints() << " id: " << tTrack->ID() << " " << tTrack << std::endl;
	for( size_t hh=0; hh<tTrack->NumberTrajectoryPoints(); hh++) {
	  std::cerr << hh << " " << tTrack->LocationAtPoint(hh)[0] << ", " <<
	    tTrack->LocationAtPoint(hh)[1] << ", " <<
	    tTrack->LocationAtPoint(hh)[2] << std::endl;
	}
	continue; // I don't want to deal with these bullshit "tracks"
      }


      // Let's define some sort of track isolation
      // isolation near end points is more important






    


     float yDist = std::min( 1.0*geo->DetHalfHeight() - trackEndPt1_Y, 1.0*geo->DetHalfHeight() - trackEndPt2_Y );


     mf::LogInfo("CosmicTagger Info") << "X values at endpts: " << trackEndPt1_X << " " << trackEndPt2_X ;
     mf::LogInfo("CosmicTagger Info") << "Y values at endpts: " << trackEndPt1_Y << " " << trackEndPt2_Y << " yhi:"<< 1.0*geo->DetHalfHeight();
     mf::LogInfo("CosmicTagger Info") << "Z values at endpts: " << trackEndPt1_Z << " " << trackEndPt2_Z ;
     mf::LogInfo("CosmicTagger Info") << "Trk length: " << tTrack->Length() ;

      //      float minDistPt1 = 99;
      //float minDistPt2 = 99;
      //      float dy1 = std::min(geo->DetHalfHeight() - trackEndPt1_Y, trackEndPt1_Y + geo->DetHalfHeight());
      //      float dy2 = std::min(geo->DetHalfHeight() - trackEndPt2_Y, trackEndPt2_Y + geo->DetHalfHeight());
      //      float dx1 = std::min(2.0*geo->DetHalfWidth() - trackEndPt1_X, 1.0*trackEndPt1_X );
      //      float dx2 = std::min(2.0*geo->DetHalfWidth() - trackEndPt2_X, 1.0*trackEndPt2_X );
      //      float dz1 = std::min(geo->DetLength() - trackEndPt1_Z, 1.0*trackEndPt1_Z );
      //      float dz2 = std::min(geo->DetLength() - trackEndPt2_Z, 1.0*trackEndPt2_Z );
      
      //      minDistPt1 = std::min(dx1,dy1);
      //      minDistPt1 = std::min(minDistPt1, dz1);
      //      minDistPt2 = std::min(dx2,dy2);
      //      minDistPt2 = std::min(minDistPt2, dz2);

      mf::LogInfo("CosmicTagger Info") << geo->DetHalfHeight() << " y's"<<trackEndPt1_Y << ": " << geo->DetHalfHeight() - trackEndPt1_Y 
				       << " " <<  trackEndPt1_Y + geo->DetHalfHeight();
      mf::LogInfo("CosmicTagger Info") << geo->DetHalfHeight() << " y's"<<trackEndPt2_Y << ": " << geo->DetHalfHeight() - trackEndPt2_Y 
				       << " " << trackEndPt2_Y + geo->DetHalfHeight() ;
      mf::LogInfo("CosmicTagger Info") << geo->DetHalfWidth()  << " x's"<<trackEndPt1_X << ": " << 2.0*geo->DetHalfWidth() - trackEndPt1_X 
				       << " " << 1.0*trackEndPt1_X 		;
      mf::LogInfo("CosmicTagger Info") << geo->DetHalfWidth()  << " x's"<<trackEndPt2_X << ": " << 2.0*geo->DetHalfWidth() - trackEndPt2_X 
				       << " " <<  1.0*trackEndPt2_X 		 ;
      mf::LogInfo("CosmicTagger Info") << geo->DetLength()     << " z's"<<trackEndPt1_Z << ": " << geo->DetLength() - trackEndPt1_Z<< " " 
				       <<  1.0*trackEndPt1_Z 			 ;
      mf::LogInfo("CosmicTagger Info") << geo->DetLength()     << " z's"<<trackEndPt2_Z << ": " << geo->DetLength() - trackEndPt2_Z<< " " 
				       <<  1.0*trackEndPt2_Z                        ;
      

      //      std::cerr << "distance to detector boundary for endpt1: " <<  minDistPt1 << std::endl;
      //      std::cerr << "distance to detector boundary for endpt2: " <<  minDistPt2 << std::endl;



      //////////////////////////////////////////////////////////////////////
      //Let's check Cluster connections to pre & post spill planes first

      // THIS REQUIRES THAT WE CHECK THE CLUSTERS WITHIN THIS TRACK LOOP
      // SO DO SOMETHING LIKE:
      std::vector <double> t1Times, t2Times;
      std::vector <int> fail;
      if(fDoClusterCheck) doTrackClusterCheck( ClusterVect, t1Times, t2Times, fail );


      if( count( fail.begin(), fail.end(), 1 ) > 0 ) {
	mf::LogInfo("CosmicTaggerAna Info") << "Track is failing cluster check!!! " ;
	isCosmic = 4;
      }




    
      /////////////////////////////////
      // Now check Y & Z boundaries:
      /////////////////////////////////
      int nBd = 0;
      float bndDist = 5;
      if(fabs(trackEndPt1_Y - geo->DetHalfHeight())<bndDist ) nBd++;
      if(fabs(trackEndPt2_Y + geo->DetHalfHeight())<bndDist ) nBd++;
      if(fabs(trackEndPt1_Z - geo->DetLength())<bndDist || fabs(trackEndPt2_Z - geo->DetLength())<bndDist ) nBd++;
      if(fabs(trackEndPt1_Z )<bndDist || fabs(trackEndPt2_Z )<bndDist ) nBd++;
      if(isCosmic==0 && nBd>1) isCosmic=1;



      /////////////////////////////////////////
      // check modified crossed two boundaries
      /////////////////////////////////////////

      float totalBoundaryDistance=0;
      /*
      std::cerr << "isCosmic: " << isCosmic << " " << face1 << " " << face2  << " " << totalBoundaryDistance  << std::endl;
      mf::LogInfo("CosmicTagger Info") << "face: " << face1 << " " << face2 ;
      if(face1==1 || face2==1) {
	totalBoundaryDistance += std::min(trackEndPt1_Z, trackEndPt2_Z);
      }
      if(face1==6 || face2==6 ) {
	totalBoundaryDistance += std::min( geo->DetLength() - trackEndPt1_Z, geo->DetLength() - trackEndPt2_Z ) ;
      }
      if( face1==2 || face2==2 ) {
	totalBoundaryDistance += std::min( trackEndPt1_X, trackEndPt2_X);
      }
      if( face1==5 || face2==5 ) {
	totalBoundaryDistance += std::min(2.*geo->DetHalfWidth() - trackEndPt1_X, 2.*geo->DetHalfWidth() - trackEndPt2_X);
      }
      if( face1==3 || face2==3) {
	totalBoundaryDistance += std::min( geo->DetHalfHeight() - trackEndPt1_Y, geo->DetHalfHeight() - trackEndPt2_Y );
      }
      if( face1==4 || face2==4 ) {
	totalBoundaryDistance += std::min( geo->DetHalfHeight() + trackEndPt1_Y, geo->DetHalfHeight() + trackEndPt2_Y );
      }
      if ( isCosmic==0 && totalBoundaryDistance > fTotalBoundaryLimit ) isCosmic=5;
      std::cerr << "isCosmic: " << isCosmic << " " << face1 << " " << face2  << " " << totalBoundaryDistance  << std::endl;

      // Check boundaries considering when it passed through
      // this assumes the track comes from above:
      double xGen1 = -999;
      double xGen2 = -999;
      if( trackTime > -990 ) {
	xGen1 = trackEndPt1_X+ (geo->DetHalfWidth()*2.0/fReadOutWindowSize)*(fReadOutWindowSize - trackTime);
	xGen2 = trackEndPt2_X+ (geo->DetHalfWidth()*2.0/fReadOutWindowSize)*(fReadOutWindowSize - endTime);
      }
      mf::LogInfo("CosmicTaggerAna Info") << "test: " <<  trackEndPt1_X << " " << xGen1 << " t1: " << trackTime
					  << ", 2nd: " << trackEndPt2_X << " " << xGen2 << " t2 " << endTime ;
      if(isCosmic==0) {
	if(xGen1<bndDist || (2.*geo->DetHalfWidth() - xGen1) < bndDist ) nBd++;
	else if(xGen2<bndDist || (2.*geo->DetHalfWidth() - xGen2) < bndDist ) nBd++;
	if(nBd>1) isCosmic = 3;
      }
      */







      //mf::LogInfo("CosmicTagger Info") << "Boundaries: "<< face1 << ", " << face2 << ", totalBoundaryDistance   " << totalBoundaryDistance;

    
    
       



      float tick1 =  9999;
      float tick2 = -9999;
      
      for (unsigned int p=0; p< HitVec.size(); p++) {
	//	  outHits->push_back(*HitVec[p]);
	if( HitVec[p]->PeakTime() < tick1 ) tick1 =  HitVec[p]->PeakTime();
	if( HitVec[p]->PeakTime() > tick2 ) tick2 =  HitVec[p]->PeakTime();
      }
      
      
      
      
      
      std::vector<float> endPt1;
      endPt1.push_back( trackEndPt1_X );
      endPt1.push_back( trackEndPt1_Y );
      endPt1.push_back( trackEndPt1_Z );
      std::vector<float> endPt2;
      endPt2.push_back( trackEndPt2_X );
      endPt2.push_back( trackEndPt2_Y );
      endPt2.push_back( trackEndPt2_Z );
      
      
      anab::CosmicTag cctt = anab::CosmicTag(endPt1, endPt2, 0, isCosmic );
      
      int nXBd =0;
      float xBnd1 = cctt.getXInteraction(endPt1[0], 2.0*geo->DetHalfWidth(), fReadOutWindowSize, trackTime, std::floor(tick1) );
      float xBnd2 = cctt.getXInteraction(endPt1[0], 2.0*geo->DetHalfWidth(), fReadOutWindowSize, trackTime, std::floor(tick2) );
      if(xBnd1 < bndDist || xBnd2 < bndDist) nXBd++;
      if( ( 2.*geo->DetHalfWidth() - xBnd1 < bndDist ) || ( 2.*geo->DetHalfWidth() - xBnd1 < bndDist ) ) nXBd++;
      
      if(tick1 < fReadOutWindowSize || tick2 > 2*fReadOutWindowSize ) isCosmic = 1;
      if( nBd>1 ) isCosmic = 2;
      if( isCosmic==0 && nXBd+nBd>1 ) isCosmic = 3;
      
      float cosmicScore = isCosmic > 0 ? 1 : 0;
    
    
    cosmicTagTrackVector->push_back( anab::CosmicTag(endPt1,
						     endPt2,
						     cosmicScore,
						     isCosmic
						     ) );

    damnAssns.push_back( anab::CosmicTag(endPt1,
					  endPt2,
					  cosmicScore,
					  isCosmic
					  ) );
    
    mf::LogInfo("CosmicTagger Results") << "The IsCosmic value is "<< isCosmic << " origin: " << origin
					<< trackEndPt1_X<<","<< trackEndPt1_Y << "," << trackEndPt1_Z<< " | | " 
					<< trackEndPt2_X<< ","<< trackEndPt2_Y <<"," << trackEndPt2_Z;
    
    /*
      anab::CosmicTag cctt = anab::CosmicTag(endPt1, endPt2, -999, cosmicScore, isCosmic );
      
      const simb::MCParticle *simpart = bt->TrackIDToParticle( tTrack->ID() );
      double xx1 = -999;
      double xx2 = -999;
      if( simpart ) {
	size_t numtraj = simpart->NumberTrajectoryPoints();
	size_t firstPt = -999;
	size_t lastPt = -999;
	for(size_t t = 0; t < numtraj; ++t){
          try{
            // check if the particle is inside a TPC              
            double pos[3] = {simpart->Vx(t), simpart->Vy(t), simpart->Vz(t)};
            unsigned int tpc   = 0;
            unsigned int cstat = 0;
            geo->PositionToTPC(pos, tpc, cstat);
          }
          catch(cet::exception &e){
            continue;
          }
	  if(firstPt==-999) firstPt=t;
	  else lastPt = t;
	}
	if(firstPt>-1) xx1 = simpart->Vx(firstPt); 
	if(lastPt>-1) xx2 = simpart->Vx(lastPt); 
      }

      int trackTime1 = endPt1[0]/(2.0*geo->DetHalfWidth()/fReadOutWindowSize);
      int trackTime2 = endPt2[0]/(2.0*geo->DetHalfWidth()/fReadOutWindowSize);

      std::cerr << "What the hell is going on " 
		<< cctt.getXInteraction(endPt1[0], 2.0*geo->DetHalfWidth(), fReadOutWindowSize, trackTime, trackTime1 )
		<< " " 
		<< cctt.getXInteraction(endPt2[0], 2.0*geo->DetHalfWidth(), fReadOutWindowSize, trackTime, trackTime2 ) 
		<< " | " <<xx1 << " " << xx2  <<std::endl;
      */




      outTracksForTags->push_back( *tTrack );
      art::PtrVector<recob::Track> tV;
      tV.clear();
      tV.push_back( tTrack );



      //util::CreateAssn(*this, e, *cosmicTagTrackVector, tV, *assnOutCosmicTagTrack );
      std::cerr << "creating (?) association  " << cosmicTagTrackVector->back().CosmicScore() << std::endl;
      util::CreateAssn(*this, e, *cosmicTagTrackVector, tTrack, *assnOutCosmicTagTrack );


      


      cTrack.event             = e.event();
      cTrack.trackNum          = iTrack;
      cTrack.org               = origin;
      cTrack.isCosmicVal       = isCosmic;
      cTrack.clusterCheck = failClusterCheck;
      cTrack.pdg               = trackPdg;
      cTrack.svd0              = -1;//c_svd(0);
      cTrack.svd1              = -1;//c_svd(1);
      cTrack.svd2              = -1;//c_svd(2);
      cTrack.x1                = trackEndPt1_X;
      cTrack.x2                = trackEndPt2_X;
      cTrack.x1gen = -1;//xGen1;
      cTrack.x2gen = -1;//xGen2;
      cTrack.y1                = trackEndPt1_Y;
      cTrack.y2                = trackEndPt2_Y;
      cTrack.z1                = trackEndPt1_Z;
      cTrack.z2                = trackEndPt2_Z;
      cTrack.length            = tTrack->Length();
      cTrack.yDistance         = yDist;
      cTrack.preSpillDistance  = -999;///preSpillDistance;
      cTrack.postSpillDistance = -999;///postSpillDistance;
      cTrack.time              = trackTime;
      cTrack.enterP = trackP;
      cTrack.face1 = face1;
      cTrack.face2 = face2;
      cTrack.totalBoundaryDistance = totalBoundaryDistance;
      cTrack.type = trackType;

      cTrack.dist_window_post.clear();
      cTrack.dist_window_pre.clear();


      cTrack.dist_window_post  = tdDist_window_post;  
      cTrack.dist_window_pre   = tdDist_window_pre;

      tree->Fill();






    }
    // END OF LOOPING OVER INSPILL TRACKS
  }





  e.put( std::move(outTracksForTags) );
  e.put( std::move(outClusters) );
  e.put( std::move(cosmicTagTrackVector) );
  e.put( std::move(assnOutCosmicTagTrack) );




  TrkVec.clear();



} // end of produce


///void trkf::CosmicTagger::fillOutOfSpillClusters() {
///
///
///  for( size_t i = 0; i < ClusterPreSpill_h->size(); i++ ) {
///
///    art::Ptr<recob::Cluster> cc(ClusterPreSpill_h, i);
///    int view = cc->View();
///    
///    dTdW_PreSpill[view].push_back( cc->dTdW() );
///
///    dQdW_PreSpill[view].push_back( cc->dQdW() );
///    //double tt = fabs(cc->StartPos()[1]-fReadOutWindowSize) < fabs(cc->EndPos()[1]-fReadOutWindowSize) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
///    //double ww = fabs(cc->StartPos()[1]-fReadOutWindowSize) < fabs(cc->EndPos()[1]-fReadOutWindowSize) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];
///
///    double tt = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
///    double ww = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];
///
///    std::cerr << "View: "<< view <<" FILLING PRE SPILL T'S & W'S.  t: " << tt << " w: " << ww 
///	      << " startT:" << cc->StartPos()[1]
///	      << " startW:" << cc->StartPos()[0]
///	      << " endT:" << cc->EndPos()[1]
///	      << " endW:" << cc->EndPos()[0] << std::endl;
///    tt = fReadOutWindowSize - tt;
///    endPtsT_PreSpill[view].push_back( tt );
///    endPtsW_PreSpill[view].push_back( ww );
///
///  }
///
///  for( size_t i = 0; i < ClusterPostSpill_h->size(); i++ ) {
///    art::Ptr<recob::Cluster> cc(ClusterPostSpill_h, i);
///    int view = cc->View();
///
///    dTdW_PostSpill[view].push_back( cc->dTdW() );
///    dQdW_PostSpill[view].push_back( cc->dQdW() );
///    //double tt = fabs(cc->StartPos()[1]-2*fReadOutWindowSize) < fabs(cc->EndPos()[1]-2*fReadOutWindowSize) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
///    //double ww = fabs(cc->StartPos()[1]-2*fReadOutWindowSize) < fabs(cc->EndPos()[1]-2*fReadOutWindowSize) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];
///
///    double tt = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
///    double ww = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];
///
///    std::cerr << "View: "<< view << "FILLING POST SPILL T'S & W'S.  t: " << tt << " w: " << ww 
///	      << " startT:" << cc->StartPos()[1]
///	      << " startW:" << cc->StartPos()[0]
///	      << " endT:" << cc->EndPos()[1]
///	      << " endW:" << cc->EndPos()[0] << std::endl;
///    tt = fReadOutWindowSize + tt;
///    endPtsT_PostSpill[view].push_back( tt );
///    endPtsW_PostSpill[view].push_back( ww );
///
///  }
///  // End of Pre- & Post-Spill Clusters
///  /////////////////////////////////////////////
///
///}


/////////////////////////////////////////////////////
// Do Cluster Check
// Does this cluster connect to another window?
/////////////////////////////////////////////////////
///void trkf::CosmicTagger::doClusterCheck(  art::FindManyP<recob::Hit> tHits, art::FindManyP<recob::SpacePoint> tSpacePoints  ) {
///  // MAYBE DO A CLUSTER CHECK HERE:
///  //  std::string clusterModuleLabel          = "cluster";
///
///
///
///  badClusterIDs.clear();
///  int timeLimit = 1000;
///
///  //  std::vector<int> endPointsWires[3];
///  //   std::vector<int> beginPointsWires[3];
///  //   std::vector<float> endPointsTicks[3];
///  //   std::vector<float> beginPointsTicks[3];
///  //   std::vector<int> IDs[3];
///
///  /////////////////////////////////
///  // NOW CHECK IN-SPILL CLUSTERS FOR 
///  // CONNECTIONS TO OTHER WINDOWS
///  ///////////////////////////////////
///  
///
///  
///  for( size_t i = 0; i < Cluster_h->size(); i++ ) {
///
///    art::Ptr<recob::Cluster> cc(Cluster_h, i);
///    int view = cc->View();
///    
///    double dtdw = cc->dTdW();
///    double dqdw = cc->dQdW();
///    
///    int failClusterCheck =0;
///    
///    // Check Low Time First
///    double t0 = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->StartPos()[1] : (double)cc->EndPos()[1];
///    double w0 = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->StartPos()[0] : (double)cc->EndPos()[0];
///    double t1 = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->EndPos()[1] : (double)cc->StartPos()[1];
///    double w1 = fabs(cc->StartPos()[1]) < fabs(cc->EndPos()[1]) ? (double)cc->EndPos()[0] : (double)cc->StartPos()[0];
///    std::cerr << "View: "<< view << "IN SPILL T'S & W'S.  t: " << t0 << " w: " << w0 
///	      << " startT:" << cc->StartPos()[1]
///	      << " startW:" << cc->StartPos()[0]
///	      << " endT:" << cc->EndPos()[1]
///	      << " endW:" << cc->EndPos()[0] << std::endl;
///    
/// 
///    // check against prespill window
///    //    if( fabs(t0-fReadOutWindowSize) < timeLimit ) { // is it near the prespill boundary?
///    if( fabs(t0) < timeLimit ) { // is it near the prespill boundary?
///      int indx=-1;
///      double dist=9999;
///      for(unsigned int j=0; j<endPtsW_PreSpill[view].size(); j++) {
///	double thisDist = sqrt( pow(t0-endPtsT_PreSpill[view].at(j) ,2) + pow(w0-endPtsW_PreSpill[view].at(j) ,2));
///	if(thisDist<dist) { indx=j; dist=thisDist;}
///      }
///      if(indx<0) continue;
///      
///      std::cerr << "closest prespill:  dtdw:" << dTdW_PreSpill[view] .at(indx) << " T:" << endPtsT_PreSpill[view].at(indx) 
///		<< " W:" << endPtsW_PreSpill[view].at(indx) << std::endl;
///      tDTDW_pre[view] = dTdW_PreSpill[view] .at(indx);
///      
///      tdDTDW_pre[view]        = fabs( dtdw - dTdW_PreSpill[view] .at(indx) );
///      tdDQDW_pre[view]        = fabs( dqdw - dQdW_PreSpill[view] .at(indx) );
///      tdT_pre   [view]        = t0 - endPtsT_PreSpill[view].at(indx) ;
///      tdW_pre   [view]        = fabs( w0 - endPtsW_PreSpill[view].at(indx) );
///      tdDist_window_pre[view] = dist;
///      
///      std::cerr << "!!!!!   Saving dist_window_pre: " << tdDist_window_pre[view] << " dT_pre:"<< tdT_pre[view] 
///		<< " dW_pre:" << tdW_pre[view] << std::endl;
///      
///      if(tdT_pre[view] <= fdTLimit && tdW_pre[view] <= fdWLimit) {
///	failClusterCheck++;
///	preSpillVetoCtr++;
///	badClusterIDs.push_back( cc->ID() );
///      }
///    } // end of prespill check
///    
///    
///    //    if( fabs(t1-2*fReadOutWindowSize) < timeLimit ) { // is it near the postspill boundary?
///    if( fabs(t1-fReadOutWindowSize) < timeLimit ) { // is it near the postspill boundary?
///      int indx=-1;
///      double dist=9999;
///      for(unsigned int j=0; j<endPtsW_PostSpill[view].size(); j++) {
///	std::cerr << "Post spill Checking T: " << endPtsT_PostSpill[view].at(j) << " W: " << endPtsW_PostSpill[view].at(j) << std::endl;
///	double thisDist = sqrt( pow(t1-endPtsT_PostSpill[view].at(j) ,2) + pow(w1-endPtsW_PostSpill[view].at(j) ,2));
///	if(thisDist<dist) { indx=j; dist=thisDist;}
///      }
///      if(indx<0) continue;
///      
///      std::cerr << "closest postspill:  dtdw:" << dTdW_PostSpill[view] .at(indx) << " T:" << endPtsT_PostSpill[view].at(indx) 
///		<< " W:" << endPtsW_PostSpill[view].at(indx) << std::endl;
///      
///      tDTDW_post[view] = dTdW_PostSpill[view] .at(indx);
///      
///      tdDTDW_post[view]        = fabs( dtdw - dTdW_PostSpill[view] .at(indx) );
///      tdDQDW_post[view]        = fabs( dqdw - dQdW_PostSpill[view] .at(indx) );
///      tdT_post   [view]        = endPtsT_PostSpill[view].at(indx) - t1 ;
///      tdW_post   [view]        = fabs( w1 - endPtsW_PostSpill[view].at(indx) );
///      tdDist_window_post[view] = dist;
///      
///      if(tdT_post[view] <= fdTLimit && tdW_post[view] <= fdWLimit) {
///	postSpillVetoCtr++;
///	failClusterCheck++;
///	badClusterIDs.push_back( cc->ID() );
///      }
///      
///    } // end of postspill check
///
///
///    float enCosmic=0, enNotCosmic=0;
///    // doing a truth check
///    int origin = -1;
///    int truthType1 =0;
///    int lessThanPre = 0, greaterThanPre = 0;
///    int lessThanPost = 0, greaterThanPost = 0;
///    std::vector< art::Ptr<recob::Hit> > tHitVec = tHits.at(i);
///    std::cerr << "=======  "  << tHitVec.size() << std::endl;
///    for ( auto const& itr : tHitVec ) {
///      std::vector<cheat::TrackIDE> eveides = bt->HitToEveID(itr);
///      double time = -1;
///      std::cerr << "eveides size() " << eveides.size() << " hit charge: " << itr->Charge() << std::endl;
///      for( size_t ee = 0; ee < eveides.size(); ee++ ) {
///	if( eveides[ee].energyFrac < 0.1 ) continue;
///	const art::Ptr<simb::MCTruth> mcpart = bt->TrackIDToMCTruth( eveides[ee].trackID );
///	// WE DON'T WANT THIS TO JUST TAKE THE LAST ORIGIN
///	// THINK ABOUT THIS A BIT MORE
///	origin = mcpart->Origin();
///	if(origin==1) enNotCosmic += (eveides[ee].energyFrac)*(eveides[ee].energy);
///	if(origin==2) enCosmic += (eveides[ee].energyFrac)*(eveides[ee].energy);
///      
///	const simb::MCParticle *simpart = bt->TrackIDToParticle( eveides[ee].trackID );
///	double trackTime = simpart->T();
///	double endTime = simpart->EndT();
///	time = trackTime;
///
///
///	//	std::cerr << "TEST times " << trackTime << " " << endTime << std::endl;
///      }
///
///      std::vector<double> xyz ;
///      try{
///	xyz = bt->HitToXYZ( itr );
///	double xWindow1 = 256 - (256/1.6e6)*time;
///	double xWindow2 = 2*256 - (256/1.6e6)*time;
///	if(xyz[0]<xWindow1) lessThanPre=1;
///	if(xyz[0]>xWindow1) greaterThanPre=1;
///	if(xyz[0]<xWindow2) lessThanPost=1;
///	if(xyz[0]>xWindow2) greaterThanPost=1;
///	
///	std::cerr << "some debugging " << xyz[0] << " " << time << " " << xWindow1 << " " << xWindow2 << " | " << origin << std::endl;
///      }
///      catch(cet::exception &e) {
///	continue;
///      }
///
///    } // over hit vector
///
///    if( ( (lessThanPre>0 && greaterThanPre>0) || (lessThanPre>0 && greaterThanPre>0) ) && origin==2 ) truthType1 =1;
///    std::cerr <<"Is this really a type 1 cosmic? " << truthType1 << std::endl;
///
///    // Is it possible to kinda check cluster endpoints in this routine?
///    //for (int j=0; j<cc->StartPos().size(); j++ ) std::cerr << "Cluster endpoints: " << cc->StartPos()[j] << " " << cc->EndPos()[j] << std::endl;
///    //     if( cc->StartPos()[1] < cc->EndPos()[1] ) {
///    //       beginPointsWires[view].push_back( cc->StartPos()[0] );
///    //       endPointsWires[view].push_back( cc->EndPos()[0] );
///    //       beginPointsTicks[view].push_back( cc->StartPos()[1] );
///    //       endPointsTicks[view].push_back( cc->EndPos()[1] );
///    //     }
///    //     else {
///    //       beginPointsWires[view].push_back( cc->EndPos()[0] );
///    //       endPointsWires[view].push_back( cc->StartPos()[0] );
///    //       beginPointsTicks[view].push_back( cc->EndPos()[1] );
///    //       endPointsTicks[view].push_back( cc->StartPos()[1] );
///    //     }
///
///    std::vector< art::Ptr<recob::SpacePoint> > tSptsVec = tSpacePoints.at(i);
///    double testXYZ[3]={0.};
///    unsigned int testTPC=0, testCRYO=0;
///    std::vector<float> test_all_x, test_all_y, test_all_z;
///    art::ServiceHandle<geo::Geometry> geo;
///
///
///    for(unsigned int j=0; j<tSptsVec.size(); j++) {
///      try{
///	testXYZ[0]= tSptsVec[j]->XYZ()[0];
///	testXYZ[1]= tSptsVec[j]->XYZ()[1];
///	testXYZ[2]= tSptsVec[j]->XYZ()[2];
///	geo->PositionToTPC(testXYZ, testTPC, testCRYO);
///      }
///      catch(cet::exception &e) {
///	continue;
///      }
///      test_all_x.push_back( tSptsVec[j]->XYZ()[0] );
///      test_all_y.push_back( tSptsVec[j]->XYZ()[1] );
///      test_all_z.push_back( tSptsVec[j]->XYZ()[2] );
///    }
///
///    if(test_all_x.size()>2) {
///      TVectorD xV;
///      float centroidX=0, centroidY=0, centroidZ=0;
///      TMatrixD mCluster(test_all_x.size(), 3);
///      
///      centroidX = accumulate(test_all_x.begin(), test_all_x.end(),0);
///      centroidY = accumulate(test_all_y.begin(), test_all_y.end(),0);
///      centroidZ = accumulate(test_all_z.begin(), test_all_z.end(),0);
///      centroidX /= 1.0*test_all_x.size();
///      centroidY /= 1.0*test_all_y.size();
///      centroidZ /= 1.0*test_all_z.size();
///      for(unsigned int ii=0; ii<test_all_x.size(); ii++) {
///	double vals[3];
///	vals[0] = test_all_x[ii] - centroidX;
///	vals[1] = test_all_y[ii] - centroidY;
///	vals[2] = test_all_z[ii] - centroidZ;
///	xV.Use(3, vals);
///	TMatrixDRow(mCluster, ii) = xV;
///      }
///     
///      TMatrixD Aw = mCluster;
///      std::cerr << "nrows: " << Aw.GetNrows() << " ncols: " << Aw.GetNcols() << std::endl;
///      if(Aw.GetNrows()<3) break;
///      TDecompSVD svd(Aw);
///      TMatrixD m_svdV = svd.GetV();
///      TVectorD c_svd = svd.GetSig();
///
///      double xyz[3], dxyz[3];//, xyzout[6];
///      xyz[0] = centroidX;
///      xyz[1] = centroidY;
///      xyz[2] = centroidZ;
///      dxyz[0] = m_svdV(0,0);
///      dxyz[1] = m_svdV(1,0);
///      dxyz[2] = m_svdV(2,0);
///
///
///      int endPt1 = -1;
///      int endPt2 = -1;
///      float maxDist = -1;
///      float aa = sqrt( pow( dxyz[0]*10,2) + pow( dxyz[1]*10, 2) + pow( dxyz[2]*10, 2) );
///      for( unsigned int sp=0; sp<test_all_x.size(); sp++) {
/// 	float dist = sqrt( pow( test_all_x[sp] - centroidX, 2) + pow( test_all_y[sp] - centroidY, 2) + pow( test_all_z[sp] - centroidZ, 2) );
/// 	float ee =  sqrt( pow(test_all_x[sp]- (xyz[0]+dxyz[0]*10),2) + pow(test_all_y[sp]- (xyz[1]+dxyz[1]*10),2) + pow(test_all_z[sp]- (xyz[2]+dxyz[2]*10),2) );
/// 	float dd =  sqrt( pow(test_all_x[sp]- (xyz[0]),2) + pow(test_all_y[sp]- (xyz[1]),2) + pow(test_all_z[sp]- (xyz[2]),2) );
/// 	float lineDist = sqrt( ee*ee - pow( (aa*aa-dd*dd)/(2*aa),2) );
///	lineDist=0;
/// 	if( dist > maxDist && lineDist< 3 ) {
/// 	  maxDist = dist;
/// 	  endPt1 = sp;
/// 	}
///      }
///
///      maxDist=-1;
///      for( unsigned int sp=0; sp<test_all_x.size(); sp++) {
/// 	float dist = sqrt( pow( test_all_x[sp] - test_all_x[endPt1], 2) + pow( test_all_y[sp] - test_all_y[endPt1], 2) + pow( test_all_z[sp] - test_all_z[endPt1], 2) );
/// 	float ee =  sqrt( pow(test_all_x[sp]- (xyz[0]+dxyz[0]*10),2) + pow(test_all_y[sp]- (xyz[1]+dxyz[1]*10),2) + pow(test_all_z[sp]- (xyz[2]+dxyz[2]*10),2) );
/// 	float dd =  sqrt( pow(test_all_x[sp]- (xyz[0]),2) + pow(test_all_y[sp]- (xyz[1]),2) + pow(test_all_z[sp]- (xyz[2]),2) );
/// 	float lineDist = sqrt( ee*ee - pow( (aa*aa-dd*dd)/(2*aa),2) );
///	lineDist=0;
/// 	if( dist > maxDist && lineDist< 3 ) {
/// 	  maxDist = dist;
/// 	  endPt2 = sp;
/// 	}
///      }
///
///      float length = sqrt( pow( test_all_x[endPt1] - test_all_x[endPt2], 2) + pow( test_all_y[endPt1] - test_all_y[endPt2], 2) + pow( test_all_z[endPt1] - test_all_z[endPt2], 2) );
///      std::cerr << test_all_x.size() << " Length:" << length << " The end points found are: " << test_all_x[endPt1] << "," << test_all_y[endPt1] << "," <<test_all_z[endPt1] << " and " << test_all_x[endPt2] << "," << test_all_y[endPt2] << "," << test_all_z[endPt2] << std::endl;
///
///
///
///    } // if more than 2 spts
///
///
///
///
///
///
///
///
///    // dtwindow, dwwindow, drwindow, distboundary1, distboundary2
///    float dtwindow = -1;
///    float dwwindow = -1;
///    float drwindow = -1;
///    std::cerr << " ** 1 ** " << std::endl;
///    if( tdDist_window_pre[view] < tdDist_window_post[view] ) {
///      dtwindow = tdT_pre[view];
///      dwwindow = tdW_pre[view];
///      drwindow = tdDist_window_pre[view];
///    }
///    else {
///      dtwindow = tdT_post[view];
///      dwwindow = tdW_post[view];
///      drwindow = tdDist_window_post[view];
///    }
///    std::cerr << " ** 2 ** " << std::endl;
///    //    art::Ptr<anab::CosmicTag> tempCT = art::Ptr<anab::CosmicTag>(dtwindow, dwwindow, drwindow, -1, -1);
///    //    cosmicTagVector.push_back( tempCT );
///    // SEL:
///    std::cerr << " ** 3 ** " << std::endl;
///    //    cosmicTagVector->push_back( anab::CosmicTag(dtwindow, dwwindow, drwindow, -1, -1) );
///    //cosmicTagVector->push_back( anab::CosmicTag() );
///      std::vector<float> endPt1;
///      endPt1.push_back( -1 );
///      endPt1.push_back( -1 );
///      endPt1.push_back( -1 );
///      std::vector<float> endPt2;
///      endPt2.push_back( -1 );
///      endPt2.push_back( -1 );
///      endPt2.push_back( -1 );
///      std::vector<int> dT;
///      //dT.assign(3,-1); // at the moment, I can't do this because tracks are not associated to clusters
///      std::vector<int> dW;
///      //dW.assign(3,-1); // "
///      std::vector<float> ddTdW;
///      // ddTdW.assign(3,-1); // "
///      std::vector<int> bitPass;
///      //bitPass.set(); // all are 1 now 
///       cosmicTagTrackVector->push_back( anab::CosmicTag(-1, -1, -1, -1,-1,
/// 						       endPt1,
/// 						       endPt2,
/// 						       dT,
/// 						       dW,
/// 						       ddTdW,
/// 						       -1, -1, -1, -1,
/// 						       bitPass
/// 						       ) );
///    std::cerr << " ** 4 ** " << std::endl;
///
///  } // loop over inspill clusters
///
///
///
///
///  /*
///  // ENDPOINTS CHECK
///  // "Begin" is just the lower wire number
///  for( int a=0; a<endPointsWires[0].size(); a++ ) {
///  float beginTickA = beginPointsTicks[0].at(a);
///  float endTickA = endPointsTicks[0].at(a);
///
///  for( int b=0; b<endPointsWires[1].size(); b++ ) {
///  float beginTickB = beginPointsTicks[1].at(b);
///  float endTickB = endPointsTicks[1].at(b);
///
///  if( fabs(beginTickA-beginTickB) < 5 && fabs(endTickA-endTickB) < 5) {
///  for( int c=0; c<endPointsWires[2].size(); c++ ) {
///  float beginTickC = beginPointsTicks[2].at(c);
///  float endTickC = endPointsTicks[2].at(c);
///  if( fabs(beginTickC-beginTickB) < 5 && fabs(endTickC-endTickB) < 5) std::cerr << "found match" << std::endl;
///  }
///  }
///  }
///  }
///  */
///
///  // CHECK ASSOCIATED SPACEPOINTS
///
///
///
///  // End of in-spill cluster check
///  /////////////////////////////////
///
///}
///////////// END OF CLUSTER CHECK


void trkf::CosmicTagger::doTrackClusterCheck( std::vector< art::Ptr< recob::Cluster> > ClusterVect, 
					      std::vector<double> &t1Times, std::vector<double> &t2Times, 
					      std::vector<int> &fail  ) {
  // Let's have this return the earliest and latest t0 times
  // Also some sort of boolean for if these are outside of the
  // spill window

  t1Times.clear();
  t2Times.clear();
  fail.clear();
  int timeLimit = 5;

  for( unsigned int k = 0; k < ClusterVect.size(); k++ ) {
    fail.push_back(0);

    double t0 = ClusterVect[k]->StartPos()[1] < ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[1] : ClusterVect[k]->EndPos()[1];
    double t1 = ClusterVect[k]->StartPos()[1] > ClusterVect[k]->EndPos()[1] ? ClusterVect[k]->StartPos()[1] : ClusterVect[k]->EndPos()[1]; 

    t1Times.push_back(t0);
    t1Times.push_back(t1);

    if( t0+timeLimit < fReadOutWindowSize ) { // This is into the pre-spill window
      fail.at(k)=1;
    }
    if( t0-timeLimit > 2*fReadOutWindowSize ) { // This is into the post-spill window
      fail.at(k)=1;
    }

    std::cerr << "----------- Times are: " << t0 << ", " << t1 << " fail? " << fail.at(k) << " " << fReadOutWindowSize << std::endl;

  }

  return;
}






void trkf::CosmicTagger::doSomeSpacePointStuff(art::FindManyP<recob::SpacePoint> sptsSpill , int iTrack, int &face1, int &face2) {


  std::vector<art::Ptr<recob::SpacePoint> > spts = sptsSpill.at(iTrack);
  
  std::cerr << " SPACE POINTS SIZE IS : " << spts.size() << std::endl;
  //      if(spts.size()<3 ) continue;
  
  art::ServiceHandle<geo::Geometry> geo;
  
  
  if(spts.size()>3) {
    
    all_x.clear();
    all_y.clear();
    all_z.clear();
    
    
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
    std::cerr << "nrows: " << Aw.GetNrows() << " ncols: " << Aw.GetNcols() << std::endl;
    if(Aw.GetNrows()<3) return;
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
    
    

    /*
    
    // see if the track should pass two y planes
    bool shouldPassTwoYPlanes=true;
    float t = (geo->DetHalfHeight()-xyz[1])/dxyz[1];
    float testX = xyz[0]+dxyz[0]*t;
    float testZ = xyz[2]+dxyz[2]*t;
    if( testX<0 || testX>2.*geo->DetHalfWidth() ) shouldPassTwoYPlanes=false;
    if( testZ<0 || testZ>geo->DetLength() )       shouldPassTwoYPlanes=false;
    t = (geo->DetHalfHeight()-xyz[1])/dxyz[1];
    testX = xyz[0]+dxyz[0]*t;
    testZ = xyz[2]+dxyz[2]*t;
    if( testX<0 || testX>2.*geo->DetHalfWidth() ) shouldPassTwoYPlanes=false;
    if( testZ<0 || testZ>geo->DetLength() )       shouldPassTwoYPlanes=false;
    */
    
    
    double xyzout1[3], xyzout2[3];
    geo::ProjectToBoxEdge( xyz, dxyz, 0, 2.*geo->DetHalfWidth(), -1.*geo->DetHalfHeight(), geo->DetHalfHeight(), 0, geo->DetLength(), xyzout1 );
    double ndxyz[3]; 
    ndxyz[0] = -dxyz[0];
    ndxyz[1] = -dxyz[1];
    ndxyz[2] = -dxyz[2];
    geo::ProjectToBoxEdge( xyz, ndxyz, 0, 2.*geo->DetHalfWidth(), -1.*geo->DetHalfHeight(), geo->DetHalfHeight(), 0, geo->DetLength(), xyzout2 );
    mf::LogInfo("CosmicTaggerAna Info") << "1 project to box edges: " << xyzout1[0] << " " << xyzout1[1] << " " << xyzout1[2] ;
    mf::LogInfo("CosmicTaggerAna Info") << "2 project to box edges: " << xyzout2[0] << " " << xyzout2[1] << " " << xyzout2[2] ;
    if( xyzout1[0]<2 ) face1=2;
    else if( fabs( 2.*geo->DetHalfWidth()-xyzout1[0] ) < 2 ) face1=5;
    else if( fabs( geo->DetHalfHeight() - xyzout1[1] ) < 2 ) face1=3;
    else if( fabs( xyzout1[1] +geo->DetHalfHeight() ) < 2 ) face1=4;
    else if( fabs( xyzout1[2] ) < 2 ) face1=1;
    else if( fabs( geo->DetLength() - xyzout1[2] ) < 2 ) face1=6;
    
    if( xyzout2[0]<2 ) face2=2;
    else if( fabs( 2.*geo->DetHalfWidth()-xyzout2[0] ) < 2 ) face2=5;
    else if( fabs( geo->DetHalfHeight() - xyzout2[1] ) < 2 ) face2=3;
    else if( fabs( xyzout2[1] + geo->DetHalfHeight() ) < 2 ) face2=4;
    else if( fabs( xyzout2[2] ) < 2 ) face2=1;
    else if( fabs( geo->DetLength() - xyzout2[2] ) < 2 ) face2=6;
    
    
    
    if( face1<0 || face2<0 ) {
      mf::LogInfo("CosmicTagger Info") << "!!!!!!!!!!!! Problem projecting to box edge!!!! " ;
      mf::LogInfo("CosmicTagger Info") << xyzout1[0]                                ;
      mf::LogInfo("CosmicTagger Info") << fabs( 2.*geo->DetHalfWidth()-xyzout1[0] ) ;
      mf::LogInfo("CosmicTagger Info") << fabs( geo->DetHalfHeight() - xyzout1[1] ) ;
      mf::LogInfo("CosmicTagger Info") << fabs( xyzout1[1] + geo->DetHalfHeight() ) ;
      mf::LogInfo("CosmicTagger Info") << fabs( xyzout1[2] )			    ;
      mf::LogInfo("CosmicTagger Info") << fabs( geo->DetLength() - xyzout1[2] )     ;
      
      mf::LogInfo("CosmicTagger Info") << "---- " ;
      mf::LogInfo("CosmicTagger Info") << xyzout2[0]                                ;
      mf::LogInfo("CosmicTagger Info") << fabs( 2.*geo->DetHalfWidth()-xyzout2[0] ) ;
      mf::LogInfo("CosmicTagger Info") << fabs( geo->DetHalfHeight() - xyzout2[1] ) ;
      mf::LogInfo("CosmicTagger Info") << fabs( xyzout2[1] + geo->DetHalfHeight() ) ;
      mf::LogInfo("CosmicTagger Info") << fabs( xyzout2[2] )			    ;
      mf::LogInfo("CosmicTagger Info") << fabs( geo->DetLength() - xyzout2[2] )     ;
      
      
    }
  }// spts > 3 so we can to the matrix stuff
  
  
  
} // doSomeSpacePointStuff



// void trkf::CosmicTagger::truthCheck(int trackID, int &origin, int &trackPdg, double &trackTime, double &endTime, float &trackP) {
  
//   std::cerr << "Doing Truth Check..." << std::endl;
//   const art::Ptr<simb::MCTruth> mcpart = bt->TrackIDToMCTruth( trackID );
//   origin = mcpart->Origin();
//   const simb::MCParticle *simpart = bt->TrackIDToParticle( trackID );


//   if( simpart ) {
//     trackPdg = simpart->PdgCode();
//     trackTime = simpart->T();
//     endTime = simpart->EndT();
//     trackTime = trackTime*1.0/fSamplingRate;
//     endTime = endTime*1.0/fSamplingRate;
//     trackP = simpart->Momentum().P();
//     trackP = simpart->E();
//     mf::LogInfo("CosmicTaggerAna Info") << "energy is " << trackP ;
//   }
//   else {
//     std::cerr << "Backtracker is failing to get the MCParticle..." << std::endl;
//   }
  
//}// end of truthCheck




void trkf::CosmicTagger::beginJob() {

  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<geo::Geometry> geo;

  tree = tfs->make<TTree>("CosmicTree", "CosmicTree");
  tree->Branch("event",  &cTrack.event          , "event/I"); 
  tree->Branch("trackNum", &cTrack.trackNum     , "trackNum/I")             ;
  tree->Branch("org",  &cTrack.org              , "org/I");                
  tree->Branch("isCosmicVal",  &cTrack.isCosmicVal    , "isCosmicVal/I");        
  tree->Branch("clusterCheck",  &cTrack.clusterCheck  , "clusterCheck/I");        
  tree->Branch("pdg", &cTrack.pdg         , "pdg/I");
  tree->Branch("svd0",  &cTrack.svd0      , "svd0/F");             
  tree->Branch("svd1",  &cTrack.svd1      , "svd1/F");             
  tree->Branch("svd2",  &cTrack.svd2      , "svd2/F");             
  tree->Branch("x1",  &cTrack.x1          , "x1/F");               
  tree->Branch("x2",  &cTrack.x2          , "x2/F");               
  tree->Branch("x1gen",  &cTrack.x1gen    , "x1gen/F");               
  tree->Branch("x2gen",  &cTrack.x2gen    , "x2gen/F");               
  tree->Branch("y1",  &cTrack.y1               , "y1/F");               
  tree->Branch("y2",  &cTrack.y2               , "y2/F");               
  tree->Branch("z1",  &cTrack.z1               , "z1/F");               
  tree->Branch("z2",  &cTrack.z2               , "z2/F");               
  tree->Branch("length",  &cTrack.length       , "length/F");           
  tree->Branch("yDistance",  &cTrack.yDistance , "yDistance/F");        
  tree->Branch("preSpillDistance",  &cTrack.preSpillDistance  , "preSpillDistance/F"); 
  tree->Branch("postSpillDistance",  &cTrack.postSpillDistance, "postSpillDistance/F");
  tree->Branch("time",  &cTrack.time      , "time/F");
  tree->Branch("enterP",  &cTrack.enterP  , "enterP/F");
  tree->Branch("face1",  &cTrack.face1    , "face1/I");
  tree->Branch("face2",  &cTrack.face2    , "face2/I");
  tree->Branch("totalBoundaryDistance",  &cTrack.totalBoundaryDistance, "totalBoundaryDistance/F");
  tree->Branch("type",  &cTrack.type      , "type/I");

  tree->Branch("dist_window_post", "std::vector<double>", &cTrack.dist_window_post    );
  tree->Branch("dist_window_pre" , "std::vector<double>", &cTrack.dist_window_pre    );



}

void trkf::CosmicTagger::reconfigure(fhicl::ParameterSet const & p) {
  // Implementation of optional member function here.
  
  ////////  fSptalg                = new trkf::SpacePointAlg(p.get<fhicl::ParameterSet>("SpacePointAlg"));


  art::ServiceHandle<util::DetectorProperties> detp;
  fReadOutWindowSize = detp->ReadOutWindowSize();
  fSamplingRate = detp->SamplingRate();
  fTotalBoundaryLimit = p.get<float>("TPCBoundaryLimit", 15);
  f3DSpillDistance    = p.get<float>("SpillDistance",12); 

  fSpillVetoCtr = p.get<int>("SpillVetoCounter", 2);
  fdTLimit      = p.get<int>("dTLimit",8);
  fdWLimit      = p.get<int>("dWLimit",8);

  fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel", "cluster");
  fTrackModuleLabel   = p.get< std::string >("TrackModuleLabel", "track");
  fDoTrackCheck = p.get< int >("DoTrackCheck", 0);
  fDoClusterCheck = p.get< int >("DoClusterCheck", 1);
  fClusterAssociatedToTracks = p.get< int >("ClustersAssociatedToTracks",1);
}

void trkf::CosmicTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(trkf::CosmicTagger)
