////////////////////////////////////////////////////////////////////////
// CosmicTrackerAlg.cxx
//
// tjyang@fnal.gov
//
// CosmicTrackerAlg class
//
/////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "CosmicTrackerAlg.h"

#include "TH1D.h"
#include "TFile.h"

#include <iostream>

bool SortByMultiplet(art::Ptr<recob::Hit> const& a, 
		     art::Ptr<recob::Hit> const& b){
  // compare the wire IDs first:
  if (a->WireID().Plane != b->WireID().Plane) 
    return a->WireID().Plane < b->WireID().Plane;
  if (a->WireID().Wire != b->WireID().Wire)
    return a->WireID().Wire > b->WireID().Wire;
  if (a->StartTick() != b->StartTick()) return a->StartTick() < b->StartTick();
  // if still undecided, resolve by local index
  return a->LocalIndex() < b->LocalIndex(); // if still unresolved, it's a bug!
}   

namespace trkf{

  CosmicTrackerAlg::CosmicTrackerAlg(fhicl::ParameterSet const& pset){
    this->reconfigure(pset);
    hitindex.resize(3);
    w0.resize(3);
    w1.resize(3);
    t0.resize(3);
    t1.resize(3);
  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::reconfigure(fhicl::ParameterSet const& pset){
    fMinNumCluHits = pset.get< unsigned int> ("MinNumCluHits");
    fSegmentSize   = pset.get< unsigned int >("SegmentSize");
  }
  
  //---------------------------------------------------------------------
  void CosmicTrackerAlg::SPTReco(std::vector<art::Ptr<recob::Hit> >&fHits){

    std::reverse(fHits.begin(), fHits.end());

    trkPos.clear();
    trkDir.clear();

    // Track hit X and WireIDs in each plane
    std::array<std::vector<std::pair<double, geo::WireID>>,3> trkXW;
    // Track hit charge ...
    std::array<std::vector<double>,3> trkChg;
    std::array< std::vector<art::Ptr<recob::Hit>>,3> trkHits;
    for (size_t i = 0; i<fHits.size(); ++i){
      trkHits[fHits[i]->WireID().Plane].push_back(fHits[i]);
    }
    // make the track trajectory
    for(size_t ipl = 0; ipl < 3; ++ipl) {
      trkXW[ipl].resize(trkHits[ipl].size());
      trkChg[ipl].resize(trkHits[ipl].size());
      for(size_t iht = 0; iht < trkHits[ipl].size(); ++iht) {
        double xx = detprop->ConvertTicksToX(trkHits[ipl][iht]->PeakTime(), 
					     ipl, trkHits[ipl][iht]->WireID().TPC, 
					     trkHits[ipl][iht]->WireID().Cryostat);
        trkXW[ipl][iht] = std::make_pair(xx, trkHits[ipl][iht]->WireID());
        trkChg[ipl][iht] = trkHits[ipl][iht]->Integral();
      } // iht
    } // ip
    fTrackTrajectoryAlg.TrackTrajectory(trkXW, trkPos, trkDir, trkChg);
    for (size_t i = 0; i< trkPos.size(); ++i){
      std::cout<<trkPos[i].X()<<" "<<trkPos[i].Y()<<" "<<trkPos[i].Z()<<std::endl;
    }

    /*
    usehit.assign(fHits.size(),0);
    spx.assign(fHits.size(),-9999);
    spy.assign(fHits.size(),-9999);
    spz.assign(fHits.size(),-9999);

    //sort input hits by wire ID number, 
    //then by start of the region of interest in time, then by the multiplet
    std::reverse(fHits.begin(), fHits.end());
    //std::sort(fHits.begin(), fHits.end(), &SortByMultiplet);
    for (auto &ihit: fHits){
      std::cout<<ihit->WireID().Plane<<" "<<ihit->WireID().Wire<<" "<<ihit->PeakTime()<<std::endl;
    }

    for (size_t i = 0; i<fHits.size(); ++i){
      hitindex[fHits[i]->WireID().Plane].push_back(i);
    }
    for (size_t i = 0; i<3; ++i) std::cout<<i<<" "<<hitindex[i].size()<<std::endl;

    GetClusterSegment(fHits);

    int ipl = -1;
    int jpl = -1;
    int kpl = -1;
    unsigned int maxnhits = 0;
    unsigned int secondmaxnhits = 0;
    for (size_t i = 0; i<3; ++i){
      if (hitindex[i].size()>maxnhits){
	jpl = ipl;
	secondmaxnhits = maxnhits;
	maxnhits = hitindex[i].size();
	ipl = i;
      }
      else if (hitindex[i].size()>secondmaxnhits){
	secondmaxnhits = hitindex[i].size();
	jpl = i;
      }
    }
    if (geom->Nplanes() == 3){
      if ((ipl==0&&jpl==1)||(ipl==1&&jpl==0)) kpl = 2;
      if ((ipl==0&&jpl==2)||(ipl==2&&jpl==0)) kpl = 1;
      if ((ipl==1&&jpl==2)||(ipl==2&&jpl==1)) kpl = 0;
    }

    std::cout<<ipl<<" "<<jpl<<" "<<kpl<<std::endl;
    
    if (ipl==-1||jpl==-1){
      mf::LogWarning("ClusterCrawlerAlg") << "Need at least two views to reconstruct a track";
      return;
    }
    if (hitindex[ipl].size() < fMinNumCluHits||
	hitindex[jpl].size() < fMinNumCluHits){
      mf::LogWarning("ClusterCrawlerAlg") << "Too few hits to recontruct tracks.";
      return;
    }

    GetSPs(fHits,ipl,jpl);
    GetSPs(fHits,jpl,ipl);
    if (kpl!=-1&& hitindex[kpl].size() >= fMinNumCluHits){
      GetSPs(fHits,kpl,ipl);
    }
    */
  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::GetClusterSegment(std::vector<art::Ptr<recob::Hit> >&fHits){


    for (size_t i = 0; i<3; ++i){
      if (hitindex[i].size() < fMinNumCluHits) continue;
      w0[i].push_back(fHits[hitindex[i][0]]->WireID().Wire);
      t0[i].push_back(fHits[hitindex[i][0]]->PeakTime()
		      - detprop->GetXTicksOffset(fHits[hitindex[i][0]]->WireID().Plane,
						 fHits[hitindex[i][0]]->WireID().TPC,
						 fHits[hitindex[i][0]]->WireID().Cryostat));
      unsigned int nhits = 0;
      for (size_t j = 1; j < hitindex[i].size(); ++j){
	if (nhits<fSegmentSize){
	  ++nhits;
	}
	else{
	  w1[i].push_back(fHits[hitindex[i][j]]->WireID().Wire);
	  t1[i].push_back(fHits[hitindex[i][j]]->PeakTime()
			  - detprop->GetXTicksOffset(fHits[hitindex[i][j]]->WireID().Plane,
						     fHits[hitindex[i][j]]->WireID().TPC,
						     fHits[hitindex[i][j]]->WireID().Cryostat));
	  w0[i].push_back(fHits[hitindex[i][j]]->WireID().Wire);
	  t0[i].push_back(fHits[hitindex[i][j]]->PeakTime()
			  - detprop->GetXTicksOffset(fHits[hitindex[i][j]]->WireID().Plane,
						     fHits[hitindex[i][j]]->WireID().TPC,
						     fHits[hitindex[i][j]]->WireID().Cryostat));
	  nhits = 0;
	}
      }
      w1[i].push_back(fHits[hitindex[i][hitindex[i].size()-1]]->WireID().Wire);
      t1[i].push_back(fHits[hitindex[i][hitindex[i].size()-1]]->PeakTime()
		      - detprop->GetXTicksOffset(fHits[hitindex[i][hitindex[i].size()-1]]->WireID().Plane,
						 fHits[hitindex[i][hitindex[i].size()-1]]->WireID().TPC,
						 fHits[hitindex[i][hitindex[i].size()-1]]->WireID().Cryostat));
    }
    
    for (size_t i = 0; i<3; ++i){
      for (size_t j = 0; j<w0[i].size(); ++j){
	std::cout<<i<<" "<<w0[i][j]<<" "<<w1[i][j]<<" "<<t0[i][j]<<" "<<t1[i][j]<<std::endl;
      }
    }
  }
 
  //---------------------------------------------------------------------
  void CosmicTrackerAlg::GetSPs(std::vector<art::Ptr<recob::Hit> >&fHits, 
				unsigned int ipl0,
				unsigned int ipl1){


    for (size_t hit = 0; hit<hitindex[ipl0].size(); ++hit){
      geo::WireID wireid = fHits[hitindex[ipl0][hit]]->WireID();
      double t = fHits[hitindex[ipl0][hit]]->PeakTime()
	- detprop->GetXTicksOffset(fHits[hitindex[ipl0][hit]]->WireID().Plane,
				   fHits[hitindex[ipl0][hit]]->WireID().TPC,
				   fHits[hitindex[ipl0][hit]]->WireID().Cryostat);
      size_t iseg = 0;
      for (; iseg<t0[ipl1].size(); ++iseg){
	if ((t>=t0[ipl1][iseg]&&t<=t1[ipl1][iseg])||
	    (t>=t1[ipl1][iseg]&&t<=t0[ipl1][iseg])){
	  if (t0[ipl1][iseg]-t1[ipl1][iseg]){
	    double ww = (t-t1[ipl1][iseg])/(t0[ipl1][iseg]-t1[ipl1][iseg])*(w0[ipl1][iseg]-w1[ipl1][iseg])+w1[ipl1][iseg];
	    int nwires = geom->Nwires(fHits[hitindex[ipl1][0]]->WireID().Plane,
				      fHits[hitindex[ipl1][0]]->WireID().TPC,
				      fHits[hitindex[ipl1][0]]->WireID().Cryostat);
	    int wire1 = int(ww);
	    int wire2 = wire1+1;
	    bool wire1ok = false;
	    bool wire2ok = false;
	    double y1 = 0;
	    double y2 = 0;
	    double z1 = 0;
	    double z2 = 0;
	    if (wire1>=0&&wire1<nwires){
	      geo::WireID wireid1(fHits[hitindex[ipl1][0]]->WireID().Cryostat,
				  fHits[hitindex[ipl1][0]]->WireID().TPC,
				  fHits[hitindex[ipl1][0]]->WireID().Plane, wire1);
	      geo::WireIDIntersection widIntersect;
	      if (geom->WireIDsIntersect(wireid,wireid1,widIntersect)){
		wire1ok = true;
		y1 = widIntersect.y;
		z1 = widIntersect.z;
	      }
	    }
	    if (wire2>=0&&wire2<nwires){
	      geo::WireID wireid2(fHits[hitindex[ipl1][0]]->WireID().Cryostat,
				  fHits[hitindex[ipl1][0]]->WireID().TPC,
				  fHits[hitindex[ipl1][0]]->WireID().Plane, wire2);
	      geo::WireIDIntersection widIntersect;
	      if (geom->WireIDsIntersect(wireid,wireid2,widIntersect)){
		wire2ok = true;
		y2 = widIntersect.y;
		z2 = widIntersect.z;
	      }
	    }
	    if (wire1ok||wire2ok){
	      usehit[hitindex[ipl0][hit]] = 1;
	      spx[hitindex[ipl0][hit]] = detprop->ConvertTicksToX(fHits[hitindex[ipl0][hit]]->PeakTime(),
								  fHits[hitindex[ipl0][hit]]->WireID().Plane,
								  fHits[hitindex[ipl0][hit]]->WireID().TPC,
								  fHits[hitindex[ipl0][hit]]->WireID().Cryostat);
	      if (wire1ok&&wire2ok){
		spy[hitindex[ipl0][hit]] = (ww - wire1)*(y2-y1)+y1;
		spz[hitindex[ipl0][hit]] = (ww - wire1)*(z2-z1)+z1;
	      }
	      else if (wire1ok){
		spy[hitindex[ipl0][hit]] = y1;
		spz[hitindex[ipl0][hit]] = z1;
	      }
	      else if (wire2ok){
		spy[hitindex[ipl0][hit]] = y2;
		spz[hitindex[ipl0][hit]] = z2;
	      }
	      break;
	    }
//	    if (ww0>=0&&ww0<int(geom->Nwires(fHits[hitindex[ipl1][0]]->WireID().Plane,
//					     fHits[hitindex[ipl1][0]]->WireID().TPC,
//					     fHits[hitindex[ipl1][0]]->WireID().Cryostat))){
//	      geo::WireID wire0(fHits[hitindex[ipl1][0]]->WireID().Cryostat,
//				fHits[hitindex[ipl1][0]]->WireID().TPC,
//				fHits[hitindex[ipl1][0]]->WireID().Plane, ww0);
//	      geo::WireIDIntersection widIntersect;
//	      //std::cout<<wire.Wire<<" "<<wire0.Wire<<std::endl;
//	      if (geom->WireIDsIntersect(wire,wire0,widIntersect)){
//		std::cout<<widIntersect.y<<" "<<widIntersect.z<<std::endl;
//		spx[hitindex[ipl0][hit]] = detprop->ConvertTicksToX(fHits[hitindex[ipl0][hit]]->PeakTime(),
//								    fHits[hitindex[ipl0][hit]]->WireID().Plane,
//								    fHits[hitindex[ipl0][hit]]->WireID().TPC,
//								    fHits[hitindex[ipl0][hit]]->WireID().Cryostat);
//		spy[hitindex[ipl0][hit]] = widIntersect.y;
//		spz[hitindex[ipl0][hit]] = widIntersect.z;
//		usehit[hitindex[ipl0][hit]] = 1;
//		
//		break;
//	      }
//	    }

	  }
	}
      }//iseg
      if (iseg==t0[ipl1].size()) iseg = 0;
    }//hit
  }

}//namespace trkf
