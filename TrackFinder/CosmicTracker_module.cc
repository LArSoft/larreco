////////////////////////////////////////////////////////////////////////
//
//  CosmicTracker
//
//  Tracker to reconstruct cosmic ray muons and neutrino interactions
// 
//  tjyang@fnal.gov
// 
//
//  ** Modified by Muhammad Elnimr to check multiple TPCs for the 35t prototype.
//   
//  mmelnimr@as.ua.edu
//  April 2014
//
//  Major modifications to separate algorithms from the module and 
//  use Bruce's TrackTrajectoryAlg to get trajectory points. T. Yang
//  June 2015
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/FindManyP.h"
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
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/ClusterMatchTQ.h"
#include "RecoAlg/CosmicTrackerAlg.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1D.h"
#include "TVirtualFitter.h"

struct trkPoint{
public:
  TVector3 pos;
  TVector3 dir; //direction cosines
  art::Ptr<recob::Hit> hit;
};

bool AnglesConsistent(const TVector3& p1, const TVector3& p2, const TVector3& a1, const TVector3& a2, double angcut){
  double angle1 = a1.Angle(a2);
  if (angle1>TMath::PiOver2()) angle1 = TMath::Pi() - angle1;
  double angle2 = a1.Angle(p1-p2);
  if (angle2>TMath::PiOver2()) angle2 = TMath::Pi() - angle2;
  if (angle1<angcut&&angle2<angcut) return true;
  else return false;
}

// compare end points and directions
bool MatchTrack(const std::vector<trkPoint>& trkpts1, const std::vector<trkPoint>& trkpts2, double discut, double angcut){
  bool match = false;
  if (!trkpts1.size()) return match;
  if (!trkpts2.size()) return match;
  if ((trkpts1[0].hit)->WireID().Cryostat == (trkpts2[0].hit)->WireID().Cryostat &&
      (trkpts1[0].hit)->WireID().TPC == (trkpts2[0].hit)->WireID().TPC) return match; 
//  art::ServiceHandle<geo::Geometry> geom;
//  const geo::TPCGeo &thetpc1 = geom->TPC((trkpts1[0].hit)->WireID().TPC, (trkpts1[0].hit)->WireID().Cryostat);
//  const geo::TPCGeo &thetpc2 = geom->TPC((trkpts2[0].hit)->WireID().TPC, (trkpts2[0].hit)->WireID().Cryostat);

  //std::cout<<trkpts1[0].pos.Y()<<" "<<trkpts1.back().pos.Y()<<" "<<trkpts1[0].pos.Z()<<" "<<trkpts1.back().pos.Z()<<std::endl;
  //std::cout<<trkpts2[0].pos.Y()<<" "<<trkpts2.back().pos.Y()<<" "<<trkpts2[0].pos.Z()<<" "<<trkpts2.back().pos.Z()<<std::endl;

//  if (thetpc1.InFiducialY(trkpts1[0].pos.Y(),5)&&
//      thetpc1.InFiducialY(trkpts1.back().pos.Y(),5)&&
//      thetpc1.InFiducialZ(trkpts1[0].pos.Z(),5)&&
//      thetpc1.InFiducialZ(trkpts1.back().pos.Z(),5)) return match;
//  //std::cout<<"pass 1"<<std::endl;
//  if (thetpc2.InFiducialY(trkpts2[0].pos.Y(),5)&&
//      thetpc2.InFiducialY(trkpts2.back().pos.Y(),5)&&
//      thetpc2.InFiducialZ(trkpts2[0].pos.Z(),5)&&
//      thetpc2.InFiducialZ(trkpts2.back().pos.Z(),5)) return match;
//  //std::cout<<"pass 2"<<std::endl;

  if (AnglesConsistent(trkpts1[0].pos, trkpts2[0].pos, 
		       trkpts1[0].dir, trkpts2[0].dir, angcut)) match = true;
  if (AnglesConsistent(trkpts1.back().pos, trkpts2[0].pos, 
		       trkpts1.back().dir, trkpts2[0].dir, angcut)) match = true;
  if (AnglesConsistent(trkpts1[0].pos, trkpts2.back().pos, 
		       trkpts1[0].dir, trkpts2.back().dir, angcut)) match = true;
  if (AnglesConsistent(trkpts1.back().pos, trkpts2.back().pos, 
		       trkpts1.back().dir, trkpts2.back().dir, angcut)) match = true;

  return match;
}

bool SortByWire (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) { 
  return h1->WireID().Wire < h2->WireID().Wire;
}

bool sp_sort_x0(const trkPoint &tp1, const trkPoint &tp2)
{
  return tp1.pos.X() < tp2.pos.X();
}

bool sp_sort_x1(const trkPoint &tp1, const trkPoint &tp2)
{
  return tp1.pos.X() > tp2.pos.X();
}

bool sp_sort_y0(const trkPoint &tp1, const trkPoint &tp2)
{
  return tp1.pos.Y() < tp2.pos.Y();
}

bool sp_sort_y1(const trkPoint &tp1, const trkPoint &tp2)
{
  return tp1.pos.Y() > tp2.pos.Y();
}

bool sp_sort_z0(const trkPoint &tp1, const trkPoint &tp2)
{
  return tp1.pos.Z() < tp2.pos.Z();
}

bool sp_sort_z1(const trkPoint &tp1, const trkPoint &tp2)
{
  return tp1.pos.Z() > tp2.pos.Z();
}

bool spt_sort_x0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] < xyz2[0];
}

bool spt_sort_x1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] > xyz2[0];
}

bool spt_sort_y0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[1] < xyz2[1];
}

bool spt_sort_y1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[1] > xyz2[1];
}

bool spt_sort_z0(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] < xyz2[2];
}

bool spt_sort_z1(const recob::SpacePoint h1, const recob::SpacePoint h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] > xyz2[2];
}


namespace trkf {

  class CosmicTracker : public art::EDProducer {
    
  public:
    
    explicit CosmicTracker(fhicl::ParameterSet const& pset);
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:

    cluster::ClusterMatchTQ  fClusterMatch;
    trkf::CosmicTrackerAlg   fCTAlg;

    std::string     fClusterModuleLabel; ///< label for input cluster collection
    
    std::string     fSortDir;            ///< sort space points 

    bool            fStitchTracks;       ///< Stitch tracks from different TPCs
    double          fDisCut;             ///< Distance cut for track merging
    double          fAngCut;             ///< Angle cut for track merging

    bool            fTrajOnly;           ///< Only use trajectory points from TrackTrajectoryAlg for debugging
  

  }; // class CosmicTracker

}

namespace trkf {

  //-------------------------------------------------
  CosmicTracker::CosmicTracker(fhicl::ParameterSet const& pset) :
    fClusterMatch(pset.get< fhicl::ParameterSet >("ClusterMatch")),
    fCTAlg(pset.get< fhicl::ParameterSet >("CTAlg"))
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
  void CosmicTracker::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
    fSortDir                = pset.get< std::string >("SortDirection","+z");
    fStitchTracks           = pset.get< bool   >("StitchTracks");
    fDisCut                 = pset.get< double >("DisCut");
    fAngCut                 = pset.get< double >("AngCut");
    fTrajOnly               = pset.get< bool   >("TrajOnly");
  }

  //-------------------------------------------------
  void CosmicTracker::beginJob()
  {
  }

  //-------------------------------------------------
  void CosmicTracker::endJob()
  {
  }

  //------------------------------------------------------------------------------------//
  void CosmicTracker::produce(art::Event& evt){
  
    // get services
    art::ServiceHandle<geo::Geometry> geom;

    std::unique_ptr<std::vector<recob::Track>      >              tcol (new std::vector<recob::Track>);           
    std::unique_ptr<std::vector<recob::SpacePoint> >                 spcol(new std::vector<recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
    std::unique_ptr<art::Assns<recob::Track, recob::Cluster> >    tcassn(new art::Assns<recob::Track, recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit> >        thassn(new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> >   shassn(new art::Assns<recob::SpacePoint, recob::Hit>);

    //double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
    //double presamplings = detprop->TriggerOffset(); // presamplings in ticks  
    //double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
    //double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
    //double Efield_drift = larprop->Efield(0);  // Electric Field in the drift region in kV/cm
    //double Temperature = larprop->Temperature();  // LAr Temperature in K

    //double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
    //double timepitch = driftvelocity*timetick;                         //time sample (cm) 


    // get input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    std::vector<art::Ptr<recob::Cluster> > clusterlist;
    if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);

    art::FindManyP<recob::Hit> fm(clusterListHandle, evt, fClusterModuleLabel);

    // find matched clusters
    fClusterMatch.ClusterMatch(clusterlist,fm);
    std::vector<std::vector<unsigned int> > &matchedclusters = fClusterMatch.matchedclusters;

    // get track space points
    std::vector<std::vector<trkPoint>> trkpts(matchedclusters.size());
    for (size_t itrk = 0; itrk<matchedclusters.size(); ++itrk){//loop over tracks

      std::vector<art::Ptr<recob::Hit> > hitlist;
      for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){//loop over clusters

        std::vector< art::Ptr<recob::Hit> > hits = fm.at(matchedclusters[itrk][iclu]);
	for (size_t ihit = 0; ihit<hits.size(); ++ihit){
	  hitlist.push_back(hits[ihit]);
	}
      }
      //reconstruct space points and directions
      fCTAlg.SPTReco(hitlist);
      if (!fTrajOnly){
	if (!fCTAlg.trkPos.size()) continue;
	for (size_t i = 0; i<hitlist.size(); ++i){
	  trkPoint trkpt;
	  trkpt.pos = fCTAlg.trkPos[i];
	  trkpt.dir = fCTAlg.trkDir[i];
	  trkpt.hit = hitlist[i];
	  trkpts[itrk].push_back(trkpt);
	}
	if (fSortDir=="+x") std::sort(trkpts[itrk].begin(),trkpts[itrk].end(),sp_sort_x0);
	if (fSortDir=="-x") std::sort(trkpts[itrk].begin(),trkpts[itrk].end(),sp_sort_x1);
	if (fSortDir=="+y") std::sort(trkpts[itrk].begin(),trkpts[itrk].end(),sp_sort_y0);
	if (fSortDir=="-y") std::sort(trkpts[itrk].begin(),trkpts[itrk].end(),sp_sort_y1);
	if (fSortDir=="+z") std::sort(trkpts[itrk].begin(),trkpts[itrk].end(),sp_sort_z0);
	if (fSortDir=="-z") std::sort(trkpts[itrk].begin(),trkpts[itrk].end(),sp_sort_z1);
      }

      if (fTrajOnly){//debug only
	if (!fCTAlg.trajPos.size()) continue;
	size_t spStart = spcol->size();
	std::vector<recob::SpacePoint> spacepoints;
	//      for (size_t ihit = 0; ihit<hitlist.size(); ++ihit){
	//	if (fCTAlg.usehit[ihit] == 1){
	for (size_t ipt = 0; ipt<fCTAlg.trajPos.size(); ++ipt){
	  art::PtrVector<recob::Hit> sp_hits;
	  //sp_hits.push_back(hitlist[ihit]);
	  double hitcoord[3];
	  hitcoord[0] = fCTAlg.trajPos[ipt].X();
	  hitcoord[1] = fCTAlg.trajPos[ipt].Y();
	  hitcoord[2] = fCTAlg.trajPos[ipt].Z();
//	  if (itrk==1){
//	    std::cout<<"hitcoord "<<hitcoord[0]<<" "<<hitcoord[1]<<" "<<hitcoord[2]<<std::endl;
//	  }
	  double err[6] = {util::kBogusD};
	  recob::SpacePoint mysp(hitcoord, 
				 err, 
				 util::kBogusD, 
				 spStart + spacepoints.size());//3d point at end of track
	  spacepoints.push_back(mysp);
	  spcol->push_back(mysp);        
	  util::CreateAssn(*this, evt, *spcol, fCTAlg.trajHit[ipt], *shassn);
	  //}//
	}//ihit
	size_t spEnd = spcol->size();
	//sort in z direction
	//std::sort(spacepoints.begin(),spacepoints.end(),sp_sort_z0);
	//std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,sp_sort_z0);
	if (fSortDir == "+x"){
	  std::sort(spacepoints.begin(),spacepoints.end(),spt_sort_x0);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,spt_sort_x0);
	}
	if (fSortDir == "-x"){
	  std::sort(spacepoints.begin(),spacepoints.end(),spt_sort_x1);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,spt_sort_x1);
	}
	if (fSortDir == "+y"){
	  std::sort(spacepoints.begin(),spacepoints.end(),spt_sort_y0);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,spt_sort_y0);
	}
	if (fSortDir == "-y"){
	  std::sort(spacepoints.begin(),spacepoints.end(),spt_sort_y1);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,spt_sort_y1);
	}
	if (fSortDir == "+z"){
	  std::sort(spacepoints.begin(),spacepoints.end(),spt_sort_z0);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,spt_sort_z0);
	}
	if (fSortDir == "-z"){
	  std::sort(spacepoints.begin(),spacepoints.end(),spt_sort_z1);
	  std::sort(spcol->begin()+spStart,spcol->begin()+spEnd,spt_sort_z1);
	}

	if(spacepoints.size()>0){
	  
	  // make a vector of the trajectory points along the track
	  std::vector<TVector3> xyz(spacepoints.size());
	  for(size_t s = 0; s < spacepoints.size(); ++s){
	    xyz[s] = TVector3(spacepoints[s].XYZ());
	  }        
	  //Calculate track direction cosines 
	  TVector3 startpointVec,endpointVec, DirCos;
	  startpointVec = xyz[0];
	  endpointVec = xyz.back();
	  DirCos = endpointVec - startpointVec;
	  //SetMag casues a crash if the magnitude of the vector is zero
	  try
	    {
	      DirCos.SetMag(1.0);//normalize vector
	    }
	  catch(...){std::cout<<"The Spacepoint is infinitely small"<<std::endl;
	    continue;
	  }
	  //std::cout<<DirCos.x()<<" "<<DirCos.y()<<" "<<DirCos.z()<<std::endl;
	  std::vector<TVector3> dircos(spacepoints.size(), DirCos);
	  
	  std::vector< std::vector<double> > dQdx;
	  std::vector<double> mom(2, util::kBogusD);
	  tcol->push_back(recob::Track(xyz, dircos, dQdx, mom, tcol->size()));
	  
	  // make associations between the track and space points
	  util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);
	  
	  // now the track and clusters
	  //util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *tcassn);
	  
	  // and the hits and track
	  std::vector<art::Ptr<recob::Hit> > trkhits;       
	  for (size_t ihit = 0; ihit<hitlist.size(); ++ihit){
	    //if (fCTAlg.usehit[ihit] == 1){
	    trkhits.push_back(hitlist[ihit]);
	    //}
	  }
	  util::CreateAssn(*this, evt, *tcol, trkhits, *thassn);
	}
      }
    }//itrk

    // by default creates a spacepoint and direction for each hit
    if (!fTrajOnly){
      std::vector<std::vector<unsigned int>> trkidx;
      //stitch tracks from different TPCs
      if (fStitchTracks){
	for (size_t itrk1 = 0; itrk1<trkpts.size(); ++itrk1){
	  for (size_t itrk2 = itrk1+1; itrk2<trkpts.size(); ++itrk2){
	    if (MatchTrack(trkpts[itrk1], trkpts[itrk2], fDisCut, fAngCut)){
	      int found1 = -1;
	      int found2 = -1;
	      for (size_t i = 0; i<trkidx.size(); ++i){
		for (size_t j = 0; j<trkidx[i].size(); ++j){
		  if (trkidx[i][j]==itrk1) found1 = i;
		  if (trkidx[i][j]==itrk2) found2 = i;
		}
	      }
	      //std::cout<<itrk1<<" "<<itrk2<<" "<<found1<<" "<<found2<<std::endl;
	      if (found1==-1&&found2==-1){
		std::vector<unsigned int> tmp;
		tmp.push_back(itrk1);
		tmp.push_back(itrk2);
		trkidx.push_back(tmp);
	      }
	      else if(found1==-1&&found2!=-1){
		trkidx[found2].push_back(itrk1);
	      }
	      else if(found1!=-1&&found2==-1){
		trkidx[found1].push_back(itrk2);
	      }
	      else if (found1!=found2){//merge two vectors
		trkidx[found1].insert(trkidx[found1].end(),
				      trkidx[found2].begin(),
				      trkidx[found2].end());
		trkidx.erase(trkidx.begin()+found2);
	      }
	    }//found match
	  }//itrk2
	}//itrk1
	for (size_t itrk = 0; itrk<trkpts.size(); ++itrk){
	  bool found = false;
	  for (size_t i = 0; i<trkidx.size(); ++i){
	    for (size_t j = 0; j<trkidx[i].size(); ++j){
	      if (trkidx[i][j]==itrk) found = true;
	    }
	  }
	  if (!found) {
	    std::vector<unsigned int> tmp;
	    tmp.push_back(itrk);
	    trkidx.push_back(tmp);
	  }
	}
      }//stitch
      else{
	trkidx.resize(trkpts.size());
	for (size_t i = 0; i<trkpts.size(); ++i){
	  trkidx[i].push_back(i);
	}
      }
      
      //make recob::track and associations
      for (size_t i = 0; i<trkidx.size(); ++i){
	//all track points
	std::vector<trkPoint> finaltrkpts;
	//all the clusters associated with the current track
	std::vector<art::Ptr<recob::Cluster>> clustersPerTrack;
	//all hits
	std::vector<art::Ptr<recob::Hit> > hitlist;
	for(size_t j = 0; j<trkidx[i].size(); ++j){
	  for (size_t k = 0; k<trkpts[trkidx[i][j]].size(); ++k){
	    finaltrkpts.push_back(trkpts[trkidx[i][j]][k]);
	    hitlist.push_back(trkpts[trkidx[i][j]][k].hit);
	  }//k
	  for (size_t iclu = 0; iclu<matchedclusters[trkidx[i][j]].size(); ++iclu){
	    art::Ptr <recob::Cluster> cluster(clusterListHandle,matchedclusters[trkidx[i][j]][iclu]);
	    clustersPerTrack.push_back(cluster);
	  }
	  
	}//j
	if (fStitchTracks){
	  if (fSortDir=="+x") std::sort(finaltrkpts.begin(),finaltrkpts.end(),sp_sort_x0);
	  if (fSortDir=="-x") std::sort(finaltrkpts.begin(),finaltrkpts.end(),sp_sort_x1);
	  if (fSortDir=="+y") std::sort(finaltrkpts.begin(),finaltrkpts.end(),sp_sort_y0);
	  if (fSortDir=="-y") std::sort(finaltrkpts.begin(),finaltrkpts.end(),sp_sort_y1);
	  if (fSortDir=="+z") std::sort(finaltrkpts.begin(),finaltrkpts.end(),sp_sort_z0);
	  if (fSortDir=="-z") std::sort(finaltrkpts.begin(),finaltrkpts.end(),sp_sort_z1);
	}
	size_t spStart = spcol->size();
	std::vector<recob::SpacePoint> spacepoints;
	for (size_t ipt = 0; ipt<finaltrkpts.size(); ++ipt){
	  art::PtrVector<recob::Hit> sp_hits;
	  sp_hits.push_back(finaltrkpts[ipt].hit);
	  double hitcoord[3];
	  hitcoord[0] = finaltrkpts[ipt].pos.X();
	  hitcoord[1] = finaltrkpts[ipt].pos.Y();
	  hitcoord[2] = finaltrkpts[ipt].pos.Z();
	  //std::cout<<"hitcoord "<<hitcoord[0]<<" "<<hitcoord[1]<<" "<<hitcoord[2]<<std::endl;
	  double err[6] = {util::kBogusD};
	  recob::SpacePoint mysp(hitcoord, 
				 err, 
				 util::kBogusD, 
				 spStart + spacepoints.size());//3d point at end of track
	  spacepoints.push_back(mysp);
	  spcol->push_back(mysp);   
	  util::CreateAssn(*this, evt, *spcol, sp_hits, *shassn);
	}//ipt
	size_t spEnd = spcol->size();
	if(spacepoints.size()>0){
	  // make a vector of the trajectory points along the track
	  std::vector<TVector3> xyz(spacepoints.size());
	  std::vector<TVector3> dircos(spacepoints.size());
	  for(size_t s = 0; s < spacepoints.size(); ++s){
	    xyz[s] = TVector3(spacepoints[s].XYZ());
	    dircos[s] = finaltrkpts[s].dir;
	    //flip direction if needed.
	    if (spacepoints.size()>1){
	      if (s==0){
		TVector3 xyz1 = TVector3(spacepoints[s+1].XYZ());
		TVector3 dir = xyz1-xyz[s];
		if (dir.Angle(dircos[s])>0.8*TMath::Pi()){
		  dircos[s] = -dircos[s];
		}
	      }
	      else{
		TVector3 dir = xyz[s]-xyz[s-1];
		if (dir.Angle(dircos[s])>0.8*TMath::Pi()){
		  dircos[s] = -dircos[s];
		}
	      }
	    }
	    //std::cout<<s<<" "<<xyz[s].X()<<" "<<xyz[s].Y()<<" "<<xyz[s].Z()<<" "<<dircos[s].X()<<" "<<dircos[s].Y()<<" "<<dircos[s].Z()<<std::endl;
	  }        
	  std::vector< std::vector<double> > dQdx;
	  std::vector<double> mom(2, util::kBogusD);
	  tcol->push_back(recob::Track(xyz, dircos, dQdx, mom, tcol->size()));
	  
	  // make associations between the track and space points
	  util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);
	  
	  // now the track and clusters
	  util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *tcassn);
	  
	  // and the hits and track
	  std::vector<art::Ptr<recob::Hit> > trkhits;       
	  for (size_t ihit = 0; ihit<hitlist.size(); ++ihit){
	    //if (fCTAlg.usehit[ihit] == 1){
	    trkhits.push_back(hitlist[ihit]);
	    //}
	  }
	  util::CreateAssn(*this, evt, *tcol, trkhits, *thassn);
	}
	
      }//i
    }
    mf::LogVerbatim("Summary") << std::setfill('-') 
                               << std::setw(175) 
                               << "-" 
                               << std::setfill(' ');
    mf::LogVerbatim("Summary") << "CosmicTracker Summary:";
    for(unsigned int i = 0; i<tcol->size(); ++i) mf::LogVerbatim("Summary") << tcol->at(i) ;
    mf::LogVerbatim("Summary") << "CosmicTracker Summary End:";
    
    evt.put(std::move(tcol));
    evt.put(std::move(spcol));
    evt.put(std::move(tspassn));
    evt.put(std::move(tcassn));
    evt.put(std::move(thassn));
    evt.put(std::move(shassn));

    return;
  }
  
  DEFINE_ART_MODULE(CosmicTracker)

} // namespace
