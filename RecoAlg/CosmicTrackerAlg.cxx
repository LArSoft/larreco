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
struct PlnLen{
  unsigned short index;
  float length;
};

bool greaterThan1 (PlnLen p1, PlnLen p2) { return (p1.length > p2.length);}

namespace trkf{

  CosmicTrackerAlg::CosmicTrackerAlg(fhicl::ParameterSet const& pset){
    this->reconfigure(pset);
  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::reconfigure(fhicl::ParameterSet const& pset){
  }
  
  //---------------------------------------------------------------------
  void CosmicTrackerAlg::SPTReco(std::vector<art::Ptr<recob::Hit> >&fHits){

    trajPos.clear();
    trajDir.clear();

    trkPos.clear();
    trkDir.clear();

    // Track hit vectors for fitting the trajectory
    std::array<std::vector<geo::WireID>,3> trkWID;
    std::array<std::vector<double>,3> trkX;
    std::array<std::vector<double>,3> trkXErr;

    std::array< std::vector<art::Ptr<recob::Hit>>,3> trajHits;

    for (size_t i = 0; i<fHits.size(); ++i){
      trajHits[fHits[i]->WireID().Plane].push_back(fHits[i]);
    }

    std::vector<PlnLen> spl;
    PlnLen plnlen;
    for(unsigned int ipl = 0; ipl < 3; ++ipl) {
      plnlen.index = ipl;
      plnlen.length = trajHits[ipl].size();
      //if(plnlen.length > 0)
      spl.push_back(plnlen);
    }
    std::sort (spl.begin(),spl.end(), greaterThan1);
    // spl[0] has the most hits and spl.back() has the least

    for (size_t ipl = 0; ipl <3; ++ipl){
      if (!trajHits[ipl].size()) continue;
      //if (ipl == spl[0].index) continue;
      double xbeg0 = detprop->ConvertTicksToX(trajHits[spl[0].index][0]->PeakTime(),
					      trajHits[spl[0].index][0]->WireID().Plane,
					      trajHits[spl[0].index][0]->WireID().TPC,
					      trajHits[spl[0].index][0]->WireID().Cryostat);
      double xend0 = detprop->ConvertTicksToX(trajHits[spl[0].index].back()->PeakTime(),
					      trajHits[spl[0].index].back()->WireID().Plane,
					      trajHits[spl[0].index].back()->WireID().TPC,
					      trajHits[spl[0].index].back()->WireID().Cryostat);
      double xbeg1 = detprop->ConvertTicksToX(trajHits[ipl][0]->PeakTime(),
					      trajHits[ipl][0]->WireID().Plane,
					      trajHits[ipl][0]->WireID().TPC,
					      trajHits[ipl][0]->WireID().Cryostat);
      double xend1 = detprop->ConvertTicksToX(trajHits[ipl].back()->PeakTime(),
					      trajHits[ipl].back()->WireID().Plane,
					      trajHits[ipl].back()->WireID().TPC,
					      trajHits[ipl].back()->WireID().Cryostat);
      double dx1 = std::abs(xbeg0-xbeg1)+std::abs(xend0-xend1);
      double dx2 = std::abs(xbeg0-xend1)+std::abs(xend0-xbeg1);
      if (std::abs(xbeg1-xend1)>0.2 // this is to make sure the track is not completely flat in t, different by ~2.5 ticks
	  &&dx2<dx1){
	std::reverse(trajHits[ipl].begin(),trajHits[ipl].end());
      }
    }	

    // make the track trajectory
    for(size_t ipl = 0; ipl < 3; ++ipl) {
      //if (ipl == spl.back().index) continue;
      trkWID[ipl].resize(trajHits[ipl].size());
      trkX[ipl].resize(trajHits[ipl].size());
      trkXErr[ipl].resize(trajHits[ipl].size());
      for(size_t iht = 0; iht < trajHits[ipl].size(); ++iht) {
        trkWID[ipl][iht] = trajHits[ipl][iht]->WireID();
        trkX[ipl][iht] = detprop->ConvertTicksToX(trajHits[ipl][iht]->PeakTime(),ipl, trajHits[ipl][iht]->WireID().TPC, trajHits[ipl][iht]->WireID().Cryostat);
        trkXErr[ipl][iht] = 0.2 * trajHits[ipl][iht]->RMS() * trajHits[ipl][iht]->Multiplicity();
      } // iht
    } // ip
    fTrackTrajectoryAlg.TrackTrajectory(trkWID, trkX, trkXErr, trajPos, trajDir);
    //remove duplicated points
    for (auto ipos = trajPos.begin(), idir = trajDir.begin(); ipos != trajPos.end() - 1; ){
      auto ipos1 = ipos +1;
      if (*ipos1 == *ipos){
	ipos = trajPos.erase(ipos);
	idir = trajDir.erase(idir);
      }
      else{
	++ipos;
	++idir;
      }
    }
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
      //std::cout<<v1.X()<<" "<<v1.Y()<<" "<<v2.X()<<" "<<v2.Y()<<" "<<v1.Mag()<<" "<<v2.Mag()<<std::endl;
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
