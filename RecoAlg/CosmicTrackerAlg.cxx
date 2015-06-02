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

/*
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
*/

namespace trkf{

  CosmicTrackerAlg::CosmicTrackerAlg(fhicl::ParameterSet const& pset){
    this->reconfigure(pset);
  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::reconfigure(fhicl::ParameterSet const& pset){
  }
  
  //---------------------------------------------------------------------
  void CosmicTrackerAlg::SPTReco(std::vector<art::Ptr<recob::Hit> >&fHits){

    //this sorting is for clustercrawler input so that the hit list starts from downstream
    std::reverse(fHits.begin(), fHits.end());

    trajPos.clear();
    trajDir.clear();

    trkPos.clear();
    trkDir.clear();

    // Track hit X and WireIDs in each plane
    std::array<std::vector<std::pair<double, geo::WireID>>,3> trajXW;
    // Track hit charge ...
    std::array<std::vector<double>,3> trajChg;
    std::array< std::vector<art::Ptr<recob::Hit>>,3> trajHits;
    for (size_t i = 0; i<fHits.size(); ++i){
      trajHits[fHits[i]->WireID().Plane].push_back(fHits[i]);
    }
    // make the track trajectory
    for(size_t ipl = 0; ipl < 3; ++ipl) {
      trajXW[ipl].resize(trajHits[ipl].size());
      trajChg[ipl].resize(trajHits[ipl].size());
      for(size_t iht = 0; iht < trajHits[ipl].size(); ++iht) {
        double xx = detprop->ConvertTicksToX(trajHits[ipl][iht]->PeakTime(), 
					     ipl, trajHits[ipl][iht]->WireID().TPC, 
					     trajHits[ipl][iht]->WireID().Cryostat);
        trajXW[ipl][iht] = std::make_pair(xx, trajHits[ipl][iht]->WireID());
        trajChg[ipl][iht] = trajHits[ipl][iht]->Integral();
      } // iht
    } // ip
    fTrackTrajectoryAlg.TrackTrajectory(trajXW, trajPos, trajDir, trajChg);
    vw.clear();
    vt.clear();
    vtraj.clear();
    vw.resize(geom->Ncryostats());
    vt.resize(geom->Ncryostats());
    vtraj.resize(geom->Ncryostats());
    for (size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      vw[cstat].resize(geom->Cryostat(cstat).NTPC());
      vt[cstat].resize(geom->Cryostat(cstat).NTPC());
      vtraj[cstat].resize(geom->Cryostat(cstat).NTPC());
      for (size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
	const geo::TPCGeo& tpcgeom = geom->Cryostat(cstat).TPC(tpc);
	vw[cstat][tpc].resize(tpcgeom.Nplanes());
	vt[cstat][tpc].resize(tpcgeom.Nplanes());
	vtraj[cstat][tpc].resize(tpcgeom.Nplanes());
	for (size_t plane = 0; plane < tpcgeom.Nplanes(); ++plane){
	  for (size_t i = 0; i< trajPos.size(); ++i){
	    double wirecord = geom->WireCoordinate(trajPos[i].Y(),
						   trajPos[i].Z(),
						   plane,tpc,cstat);
	    double tick = detprop->ConvertXToTicks(trajPos[i].X(),
						   plane,tpc,cstat);
	    //if (int(wirecord)>=0&&int(wirecord)<int(geom->Nwires(plane,tpc,cstat))
	    //&&int(tick)>=0&&int(tick)<int(detprop->ReadOutWindowSize())){
	    vw[cstat][tpc][plane].push_back(wirecord);
	    vt[cstat][tpc][plane].push_back(tick);
	    vtraj[cstat][tpc][plane].push_back(i);
	    //}//if wire and tick make sense
	  }//trajPos.size()
	}//plane
      }//tpc
    }//cstat

    MakeSPT(fHits);

  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::MakeSPT(std::vector<art::Ptr<recob::Hit> >&fHits){
    
    double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
    double Efield_drift = larprop->Efield(0);  // Electric Field in the drift region in kV/cm
    double Temperature = larprop->Temperature();  // LAr Temperature in K
    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
    double time_pitch = driftvelocity*timetick;   //time sample (cm) 

    for (size_t i = 0; i<fHits.size(); ++i){
      int ip1 = -1;
      int ip2 = -1;
      int ip2p = -1;
      int ih1 = -1;
      int ih2 = -1;
      int ih2p = -1;
      double mindis1 = 1e9;
      double mindis2 = 1e9;
      double mindis2p = 1e9;
      geo::WireID wireid = fHits[i]->WireID();
      unsigned int wire = wireid.Wire;
      unsigned int plane = wireid.Plane;
      unsigned int tpc = wireid.TPC;
      unsigned int cstat = wireid.Cryostat;
      double wire_pitch = geom->WirePitch(0,1,plane,tpc,cstat);    //wire pitch in cm
      //find the two trajectory points that enclose the hit
      //find the nearest trajectory point first
      for (size_t j = 0; j<vw[cstat][tpc][plane].size(); ++j){
	double dis = sqrt(pow(wire_pitch*(wire-vw[cstat][tpc][plane][j]),2)+
			  pow(time_pitch*(fHits[i]->PeakTime()-vt[cstat][tpc][plane][j]),2));
	if (dis<mindis1){
	  mindis1 = dis;
	  ip1 = vtraj[cstat][tpc][plane][j];
	  ih1 = j;
	}
      }
      //find the second nearest trajectory point, prefer the point on the other side
      for (size_t j = 0; j<vw[cstat][tpc][plane].size(); ++j){
	if (int(j)==ih1) continue;
	double dis = sqrt(pow(wire_pitch*(wire-vw[cstat][tpc][plane][j]),2)+
			  pow(time_pitch*(fHits[i]->PeakTime()-vt[cstat][tpc][plane][j]),2));
	if (dis<mindis2){
	  mindis2 = dis;
	  ip2 = vtraj[cstat][tpc][plane][j];
	  ih2 = j;
	}
	TVector3 v1(wire_pitch*(wire-vw[cstat][tpc][plane][ih1]),
		    time_pitch*(fHits[i]->PeakTime()-vt[cstat][tpc][plane][j]),
		    0);
	TVector3 v2(wire_pitch*(wire-vw[cstat][tpc][plane][j]),
		    time_pitch*(fHits[i]->PeakTime()-vt[cstat][tpc][plane][j]),
		    0);
	if (v1.Dot(v2)<0){
	  if (dis<mindis2p){
	    mindis2p = dis;
	    ip2p = vtraj[cstat][tpc][plane][j];
	    ih2p = j;
	  }
	}
      }
      if (ip2p != -1){
	ip2 = ip2p;
	ih2 = ih2p;
      }
      //std::cout<<i<<" "<<ip1<<" "<<ip2<<std::endl;
      if (ip1==-1||ip2==-1){
	throw cet::exception("CosmicTrackerAlg")<<"Cannot find two nearest trajectory points.";
      }
      TVector3 v1(wire_pitch*(wire-vw[cstat][tpc][plane][ih1]),
		  time_pitch*(fHits[i]->PeakTime()-vt[cstat][tpc][plane][ih1]),
		  0);
      TVector3 v2(wire_pitch*(vw[cstat][tpc][plane][ih2]-vw[cstat][tpc][plane][ih1]),
		  time_pitch*(vt[cstat][tpc][plane][ih2]-vt[cstat][tpc][plane][ih1]),
		  0);
      if (!v2.Mag()){
	throw cet::exception("CosmicTrackerAlg")<<"Two identical trajectory points.";
      }	
      double ratio = v1.Dot(v2)/v2.Mag2();
      TVector3 pos = trajPos[ip1]+ratio*(trajPos[ip2]-trajPos[ip1]);
      trkPos.push_back(pos);
      if (trajDir[ip1].Mag()){
	trkDir.push_back(trajDir[ip1]);
      }
      else if (trajDir[ip2].Mag()){
	trkDir.push_back(trajDir[ip2]);
      }
      else{//only got two trajectory points?
	trkDir.push_back((trajPos[ip2]-trajPos[ip1]).Unit());
      }
    }//i
  }

}//namespace trkf
