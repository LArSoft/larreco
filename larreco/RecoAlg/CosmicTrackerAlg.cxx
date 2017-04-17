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
namespace {
  template <typename T>
  inline T sqr(T v) { return v*v; }
} // local namespace

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
    fSPTAlg   = pset.get< int    >("SPTAlg",   0    );
    fTrajOnly = pset.get< bool   >("TrajOnly", false);
    ftmatch   = pset.get< double >("TMatch");
    fsmatch   = pset.get< double >("SMatch");
  }
  
  //---------------------------------------------------------------------
  void CosmicTrackerAlg::SPTReco(std::vector<art::Ptr<recob::Hit> >&fHits){

    trajPos.clear();
    trajDir.clear();
    trajHit.clear();

    trkPos.clear();
    trkDir.clear();

    if (fSPTAlg==0){
      TrackTrajectory(fHits);
    }
    else if (fSPTAlg==1){
      Track3D(fHits);
    }
    else{
      throw cet::exception("CosmicTrackerAlg")<<"Unknown SPTAlg "<<fSPTAlg<<", needs to be 0 or 1"; 
    }

    if (!fTrajOnly&&trajPos.size()) MakeSPT(fHits);
    else{
      trkPos = trajPos;
      trkDir = trajDir;
    }
  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::TrackTrajectory(std::vector<art::Ptr<recob::Hit> >&fHits){
    
    detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

/*
    // Track hit X and WireIDs in each plane
    std::array<std::vector<std::pair<double, geo::WireID>>,3> trajXW;
    // Track hit charge ...
    std::array<std::vector<double>,3> trajChg;
    std::array< std::vector<art::Ptr<recob::Hit>>,3> trajHits;
    for (size_t i = 0; i<fHits.size(); ++i){
      trajHits[fHits[i]->WireID().Plane].push_back(fHits[i]);
    }
*/
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
/*
    // make the track trajectory
    for(size_t ipl = 0; ipl < 3; ++ipl) {
      //if (ipl == spl.back().index) continue;
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
*/
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
    if (!trajPos.size()||!trajDir.size()){
      mf::LogWarning("CosmicTrackerAlg")<<"Failed to reconstruct trajectory points.";
      return;
    }
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

  }

  //---------------------------------------------------------------------
  void CosmicTrackerAlg::Track3D(std::vector<art::Ptr<recob::Hit> >&fHits){

    larprop = lar::providerFrom<detinfo::LArPropertiesService>();
    detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    //save time/hit information along track trajectory
    std::vector<std::map<int,double> > vtimemap(3);
    std::vector<std::map<int,art::Ptr<recob::Hit> > > vhitmap(3);

    //find hit on each wire along the fitted line
    for (size_t ihit = 0; ihit<fHits.size(); ++ihit){//loop over hits
      geo::WireID hitWireID = fHits[ihit]->WireID();
      unsigned int w = hitWireID.Wire;
      unsigned int pl = hitWireID.Plane;
      double time = fHits[ihit]->PeakTime();
      time -= detprop->GetXTicksOffset(fHits[ihit]->WireID().Plane,
				       fHits[ihit]->WireID().TPC,
				       fHits[ihit]->WireID().Cryostat);
      vtimemap[pl][w] = time;
      vhitmap[pl][w] = fHits[ihit];
    }

    // Find two clusters with the most numbers of hits, and time ranges
    int iclu1 = -1;
    int iclu2 = -1;
    int iclu3 = -1;
    unsigned maxnumhits0 = 0;
    unsigned maxnumhits1 = 0;
    
    std::vector<double> tmin(vtimemap.size());
    std::vector<double> tmax(vtimemap.size());
    for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
      tmin[iclu] = 1e9;
      tmax[iclu] = -1e9;
    }
    
    for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
      for (auto itime = vtimemap[iclu].begin(); itime!=vtimemap[iclu].end(); ++itime){
	if (itime->second>tmax[iclu]){
	  tmax[iclu] = itime->second;
	}
	if (itime->second<tmin[iclu]){
	  tmin[iclu] = itime->second;
	}
      }
      if (vtimemap[iclu].size()>maxnumhits0){
	if (iclu1!=-1){
	  iclu2 = iclu1;
	  maxnumhits1 = maxnumhits0;
	}
	iclu1 = iclu;
	maxnumhits0 = vtimemap[iclu].size();
      }
      else if (vtimemap[iclu].size()>maxnumhits1){
	iclu2 = iclu;
	maxnumhits1 = vtimemap[iclu].size();
      }
    }
    
    std::swap(iclu1,iclu2); //now iclu1 has fewer hits than iclu2
    
    //find iclu3
    for (int iclu = 0; iclu<(int)vtimemap.size(); ++iclu){
      if (iclu!=iclu1&&iclu!=iclu2) iclu3 = iclu;
    }
    
    if (iclu1!=-1&&iclu2!=-1){//at least two good clusters
      //select hits in a common time range
      auto ihit = vhitmap[iclu1].begin();
      auto itime = vtimemap[iclu1].begin();
      while (itime!=vtimemap[iclu1].end()){
	if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
	    itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
	  vtimemap[iclu1].erase(itime++);
	  vhitmap[iclu1].erase(ihit++);
	}
	else{
	  ++itime;
	  ++ihit;
	}
      }
      
      ihit = vhitmap[iclu2].begin();
      itime = vtimemap[iclu2].begin();
      while (itime!=vtimemap[iclu2].end()){
	if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
	    itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
	  vtimemap[iclu2].erase(itime++);
	  vhitmap[iclu2].erase(ihit++);
	}
	else{
	  ++itime;
	  ++ihit;
	}
      }
      
      //if one cluster is empty, replace it with iclu3
      if (!vtimemap[iclu1].size()){
	if (iclu3!=-1){
	  std::swap(iclu3,iclu1);
	}
      }
      if (!vtimemap[iclu2].size()){
	if (iclu3!=-1){
	  std::swap(iclu3,iclu2);
	  std::swap(iclu1,iclu2);
	}
      }
      if ((!vtimemap[iclu1].size())||(!vtimemap[iclu2].size())) return;
      
      bool rev = false;
      auto times1 = vtimemap[iclu1].begin();
      auto timee1 = vtimemap[iclu1].end();
      --timee1;
      auto times2 = vtimemap[iclu2].begin();
      auto timee2 = vtimemap[iclu2].end();
      --timee2;
      
      double ts1 = times1->second;
      double te1 = timee1->second;
      double ts2 = times2->second;
      double te2 = timee2->second;
      
      //find out if we need to flip ends
      if (std::abs(ts1-ts2)+std::abs(te1-te2)>std::abs(ts1-te2)+std::abs(te1-ts2)){
	rev = true;
      }

      double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
      double Efield_drift = detprop->Efield(0);  // Electric Field in the drift region in kV/cm
      double Temperature = detprop->Temperature();  // LAr Temperature in K
      
      double driftvelocity = detprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
      double timepitch = driftvelocity*timetick;         

      double wire_pitch = geom->WirePitch(0,1,
					  vhitmap[0].begin()->second->WireID().Plane,
					  vhitmap[0].begin()->second->WireID().TPC,
					  vhitmap[0].begin()->second->WireID().Cryostat);    //wire pitch in cm
      
      //find out 2d track length for all clusters associated with track candidate
      std::vector<double> vtracklength;      
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
	
        double tracklength = 0;
        if (vtimemap[iclu].size()==1){
          tracklength = wire_pitch;
        }
        else if (!vtimemap[iclu].empty()) {
          std::map<int,double>::const_iterator iw = vtimemap[iclu].cbegin(),
            wend = vtimemap[iclu].cend();
          double t0 = iw->first, w0 = iw->second;
          while (++iw != wend) {
            tracklength += std::sqrt(sqr((iw->first-w0)*wire_pitch)+sqr((iw->second-t0)*timepitch));
            w0 = iw->first;
            t0 = iw->second;
          } // while
        }
	vtracklength.push_back(tracklength);
      }
      
      std::map<int,int> maxhitsMatch;
      
      auto ihit1 = vhitmap[iclu1].begin();
      for (auto itime1 = vtimemap[iclu1].begin(); 
	   itime1 != vtimemap[iclu1].end(); 
	   ++itime1, ++ihit1){//loop over min-hits
	std::vector<art::Ptr<recob::Hit>> sp_hits;
	sp_hits.push_back(ihit1->second);
	double hitcoord[3] = {0.,0.,0.};
	double length1 = 0;
	if (vtimemap[iclu1].size()==1){
	  length1 = wire_pitch;
	}
	else{
	  for (auto iw1 = vtimemap[iclu1].begin(); iw1!=itime1; ++iw1){
	    auto iw2 = iw1;
	    ++iw2;
	    length1 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)
				 +std::pow((iw1->second-iw2->second)*timepitch,2));
	  }
	}
	double difference = 1e10; //distance between two matched hits
	auto matchedtime = vtimemap[iclu2].end();
	auto matchedhit  = vhitmap[iclu2].end();
	
	auto ihit2 = vhitmap[iclu2].begin();
	for (auto itime2 = vtimemap[iclu2].begin(); 
	     itime2!=vtimemap[iclu2].end(); 
	     ++itime2, ++ihit2){//loop over max-hits
	  if (maxhitsMatch[itime2->first]) continue;
	  double length2 = 0;
	  if (vtimemap[iclu2].size()==1){
	    length2 = wire_pitch;
	  }
	  else{
	    for (auto iw1 = vtimemap[iclu2].begin(); iw1!=itime2; ++iw1){
	      auto iw2 = iw1;
	      ++iw2;
	      length2 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)+std::pow((iw1->second-iw2->second)*timepitch,2));
	    }
	  }
	  if (rev) length2 = vtracklength[iclu2] - length2;
	  length2 = vtracklength[iclu1]/vtracklength[iclu2]*length2;
	  bool timematch = std::abs(itime1->second-itime2->second)<ftmatch;
	  if (timematch &&std::abs(length2-length1)<difference){
	    difference = std::abs(length2-length1);
	    matchedtime = itime2;
	    matchedhit = ihit2;
	  }
	}//loop over hits2
	if (difference<fsmatch){
	  hitcoord[0] = detprop->ConvertTicksToX(matchedtime->second+
						 detprop->GetXTicksOffset((matchedhit->second)->WireID().Plane,
									  (matchedhit->second)->WireID().TPC,
									  (matchedhit->second)->WireID().Cryostat),
						 (matchedhit->second)->WireID().Plane,
						 (matchedhit->second)->WireID().TPC,
						 (matchedhit->second)->WireID().Cryostat);

	  hitcoord[1] = -1e10;
	  hitcoord[2] = -1e10;

	  //WireID is the exact segment of the wire where the hit is on (1 out of 3 for the 35t)

	  geo::WireID c1=(ihit1->second)->WireID();
	  geo::WireID c2=(matchedhit->second)->WireID();
	  
	  geo::WireIDIntersection tmpWIDI;            
	  bool sameTpcOrNot=geom->WireIDsIntersect(c1,c2, tmpWIDI);
	  
	  if(sameTpcOrNot){
	    hitcoord[1]=tmpWIDI.y;
	    hitcoord[2]=tmpWIDI.z;
	  }
	  
	  if (hitcoord[1]>-1e9&&hitcoord[2]>-1e9){
	    maxhitsMatch[matchedtime->first] = 1;
	    sp_hits.push_back(matchedhit->second);
	  }
	}//if (difference<fsmatch)
	if (sp_hits.size()>1){
	  trajPos.push_back(TVector3(hitcoord));
	  trajHit.push_back(sp_hits);
	}
      }//loop over hits1
    }//if (iclu1!=-1&&iclu2!=-1){//at least two good clusters

  }


  //---------------------------------------------------------------------
  void CosmicTrackerAlg::MakeSPT(std::vector<art::Ptr<recob::Hit> >&fHits){
    
    larprop = lar::providerFrom<detinfo::LArPropertiesService>();
    detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
    double Efield_drift = detprop->Efield(0);  // Electric Field in the drift region in kV/cm
    double Temperature = detprop->Temperature();  // LAr Temperature in K
    double driftvelocity = detprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
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
	return;
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
