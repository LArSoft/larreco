#include "TCShowerAlg.h"

struct pfpStuff {
  art::Ptr<recob::PFParticle> pfp;
  art::Ptr<recob::Track> trk;
  art::Ptr<recob::Vertex> vtx;
  std::vector<art::Ptr<recob::Hit> > hits;

  std::vector<int> clsIDs;

  double score;
};

bool comparePFP(const pfpStuff& l, const pfpStuff& r) {

  art::Ptr<recob::Vertex> lvtx = l.vtx;
  art::Ptr<recob::Vertex> rvtx = r.vtx;

  double lz = l.hits.size();
  double rz = r.hits.size();
  
  // RSF: USED TO BE 50
  int hitthres = 80; // TODO: ADJUST THIS THRESHOLD 

  if (lz > hitthres && rz <= hitthres) return false;
  else if (lz <= hitthres && rz > hitthres) return true;
  return lvtx->position().Z() > rvtx->position().Z();

  //  return l.score > r.score;
}

bool compareHit(const art::Ptr<recob::Hit>& l, const art::Ptr<recob::Hit>& r) {

  int lwire = l->WireID().asWireID().Wire;
  int rwire = r->WireID().asWireID().Wire;
  /*
  double ltick = l->PeakTime();
  double rtick = r->PeakTime();
  */
  return lwire > rwire;
}

namespace shower {

  TCShowerAlg::TCShowerAlg(fhicl::ParameterSet const& pset) :
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg") ),
    fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")){
  }

  int TCShowerAlg::makeShowers(std::vector<art::Ptr<recob::PFParticle> > pfplist, std::vector<art::Ptr<recob::Vertex> > vertexlist, std::vector<art::Ptr<recob::Cluster> > clusterlist, std::vector<art::Ptr<recob::Hit> > hitlist, art::FindManyP<recob::Hit> cls_fm, art::FindManyP<recob::Cluster> clspfp_fm, art::FindManyP<recob::Vertex> vtxpfp_fm, art::FindManyP<recob::PFParticle> hit_fm, art::FindManyP<recob::Cluster> hitcls_fm, art::FindManyP<recob::Track> trkpfp_fm, art::FindManyP<anab::Calorimetry> fmcal) {

    totalEnergy.resize(2);
    totalEnergyErr.resize(2);
    dEdx.resize(2);
    dEdxErr.resize(2);   

    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<pfpStuff> allpfps;

    // put together pfparticle information
    for (size_t i = 0; i < pfplist.size(); ++i) {
      pfpStuff thispfp;
      thispfp.hits.clear();
      thispfp.clsIDs.clear();
      thispfp.pfp = pfplist[i];
      
      std::vector<art::Ptr<recob::Vertex> > thisvtxlist =  vtxpfp_fm.at(pfplist[i].key());
      if (thisvtxlist.size()) thispfp.vtx = thisvtxlist[0];

      std::vector<art::Ptr<recob::Track> > thistrklist = trkpfp_fm.at(pfplist[i].key());
      if (thistrklist.size()) thispfp.trk = thistrklist[0];

      std::vector<art::Ptr<recob::Cluster> > thisclusterlist = clspfp_fm.at(pfplist[i].key());
      std::vector<int> clustersize;

      for (size_t j = 0; j < thisclusterlist.size(); ++j) {

	thispfp.clsIDs.push_back(thisclusterlist[j]->ID());

	std::vector<art::Ptr<recob::Hit> > thishitlist = cls_fm.at(thisclusterlist[j].key());
	clustersize.push_back((int)thishitlist.size());

	for (size_t k = 0; k < thishitlist.size(); ++k) {
	  thispfp.hits.push_back(thishitlist[k]);
	} // loop through hits

      } // loop through clusters

      if (clustersize.size() == 3) {	
	if (!thispfp.vtx) continue;
	if (!thispfp.trk) continue;

	allpfps.push_back(thispfp);
	
	double tick = detprop->ConvertXToTicks(thispfp.vtx->position().X(), geo::PlaneID(0,0,2) );
        int wire = geom->WireCoordinate(thispfp.vtx->position().Y(), thispfp.vtx->position().Z(), geo::PlaneID(0,0,2));
	
	std::cout << "pfp " << thispfp.pfp->Self() + 1 << " cluster sizes " << clustersize[0] << ":" << clustersize[1] << ":" << clustersize[2] << " vertex " << thispfp.vtx->ID() << " " << tick << ":" << wire << " z " << thispfp.vtx->position().Z() << std::endl;

      } // add pfp to list

    } // loop through pfparticles

    std::sort(allpfps.begin(), allpfps.end(), comparePFP);
    std::reverse(allpfps.begin(), allpfps.end());

    std::cout << "sorted pfps: ";
    for (size_t i = 0; i < allpfps.size(); ++i)
      std::cout << allpfps[i].pfp->Self() + 1 << " ";
    std::cout << std::endl;

    bool showerCandidate = false;

    for (size_t i = 0; i < allpfps.size(); ++i) {

      showerHits.clear();

      art::Ptr<recob::Vertex> pfpvtx = allpfps[i].vtx;
      art::Ptr<recob::Track> pfptrk = allpfps[i].trk;
      std::vector<art::Ptr<recob::Hit> > pfphits = allpfps[i].hits;
      std::vector<art::Ptr<recob::Cluster> > pfpcls = clspfp_fm.at(allpfps[i].pfp.key());

      //      if (pfphits.size() == 0) continue;

      std::cout << "pfp " << allpfps[i].pfp->Self() + 1 << " hits " << pfphits.size() << std::endl;

      TVector3 vtx;
      vtx[0] = pfpvtx->position().X();
      vtx[1] = pfpvtx->position().Y();
      vtx[2] = pfpvtx->position().Z();

      if (pfptrk->Vertex()[2] < pfptrk->End()[2]) {
	shwvtx = pfptrk->Vertex();
	shwDir = pfptrk->VertexDirection();
      }
      else {
	shwvtx = pfptrk->End();
	shwDir = -pfptrk->EndDirection();
      }

      //      int tolerance = 100; // how many shower like cluster you need to define a shower              
      int tolerance = 60; // how many shower like cluster you need to define a shower              
      double pullTolerance = 0.7; // hits should be evenly distributed around the track
      double maxDist = 20; // how far a shower like cluster can be from the track
      double minDistVert = 15; // exclude tracks near the vertex

      // TODO: count number of shower-like hits hits "behind" vertex
      // if high, pick an earlier cluster and refit hits
      // this is going to take some restructuring

      if (pfphits.size() < 30) continue;
      //      if (pfphits.size() < 15) continue;
      if (pfphits.size() > 500) continue;  
      // adjust tolerances for short tracks
      if (pfphits.size() < 90) {
	tolerance = 50;
	pullTolerance = 0.9;
      }
      if (pfphits.size() > 400) tolerance = 200;
      else if (pfphits.size() > 100) tolerance = 100; // RSF added used to be 100, 100
 
      // add pfp hits to shower
      for (size_t ii = 0; ii < pfphits.size(); ++ii) {
	if ( addShowerHit(pfphits[ii], showerHits) ) showerHits.push_back(pfphits[ii]);
      } // loop over pfphits

      int nShowerHits = 0;
      double showerHitPull = 0;

      TVector3 pfpStart = shwvtx;
      TVector3 pfpPt2 = shwvtx+shwDir; // a second point along the track

      // track vertex
      std::map<geo::PlaneID, double> trk_tick1;
      std::map<geo::PlaneID, double> trk_wire1;

      // second track point       
      std::map<geo::PlaneID, double> trk_tick2;
      std::map<geo::PlaneID, double> trk_wire2;

      for (auto iPlane = geom->begin_plane_id(); iPlane != geom->end_plane_id(); ++iPlane){
        trk_tick1[*iPlane] = detprop->ConvertXToTicks(pfpStart[0], *iPlane);
        trk_wire1[*iPlane] = geom->WireCoordinate(pfpStart[1], pfpStart[2], *iPlane);
        trk_tick2[*iPlane] = detprop->ConvertXToTicks(pfpPt2[0], *iPlane);
        trk_wire2[*iPlane] = geom->WireCoordinate(pfpPt2[1], pfpPt2[2], *iPlane);
      }

      for (size_t j = 0; j < clusterlist.size(); ++j) {
	std::vector< art::Ptr<recob::Hit> > cls_hitlist = cls_fm.at(clusterlist[j].key());

	if (clusterlist[j]->ID() > 0 && cls_hitlist.size() > 10) continue;
	if (cls_hitlist.size() > 50) continue;

	bool isGoodCluster = false; // true if the hit belongs to a cluster that should be added to the shower

	bool skipit = false; // skip clusters already in the pfp
	for (size_t k = 0; k < allpfps[i].clsIDs.size(); ++k) {
	  if (allpfps[i].clsIDs[k] == clusterlist[j]->ID()) skipit = true;
	}
	if (skipit) continue;

	for (size_t jj = 0; jj < cls_hitlist.size(); ++jj) {
	  int isGoodHit = goodHit(cls_hitlist[jj], maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2);

	  if (isGoodHit == -1) {
	    isGoodCluster = false;
	    break;
	  }
	  else if (isGoodHit == 1) {
	    isGoodCluster = true;
	  }

	} // loop over hits in cluster 

	// add hits to shower  
	if (isGoodCluster) {
	  for (size_t jj = 0; jj < cls_hitlist.size(); ++jj) {
	    nShowerHits++;

	    int showerHitPullAdd = 0;
	    goodHit(cls_hitlist[jj], maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2, showerHitPullAdd);
	    showerHitPull += showerHitPullAdd;

	    if ( addShowerHit(cls_hitlist[jj], showerHits) ) showerHits.push_back(cls_hitlist[jj]);

	  } // loop over hits in cluster
	} // cluster contains hit close to track 

      } // loop through cluserlist

      showerHitPull /= nShowerHits; 

      std::cout << "shower hits " << showerHits.size() << " " << nShowerHits << " shower pull " << showerHitPull << std::endl; 

      if (nShowerHits > tolerance && std::abs(showerHitPull) < pullTolerance) {
	showerCandidate = true;
	std::cout << "SHOWER CANDIDATE" << std::endl;
	// loop over hits to find those that aren't associated with any clusters
	if (nShowerHits > 400) maxDist *= 2; // TODO: optimize this threshold              
	for (size_t k = 0; k < hitlist.size(); ++k) {
	  std::vector< art::Ptr<recob::Cluster> > hit_clslist = hitcls_fm.at(hitlist[k].key());
	  if (hit_clslist.size()) continue;

	  int isGoodHit = goodHit(hitlist[k], maxDist*2, minDistVert*2, trk_wire1, trk_tick1, trk_wire2, trk_tick2);
	  if (isGoodHit == 1 && addShowerHit(hitlist[k], showerHits) ) showerHits.push_back(hitlist[k]);
	} // loop over hits

	// loop over clusters to see if any fall within the shower
	for (size_t k = 0; k < clusterlist.size(); ++k) {
	  std::vector< art::Ptr<recob::Hit> > cls_hitlist = cls_fm.at(clusterlist[k].key());
	  //	  if (clusterlist[k]->ID() < 0) continue;
	  if (clusterlist[k]->ID() > 0 && cls_hitlist.size() > 50) continue; 

	  double thisDist = maxDist;
	  double thisMin = minDistVert;

	  if (cls_hitlist.size() < 10) {
	    thisDist *= 4;
	    thisMin *= 4;
	  }
	  else if (cls_hitlist.size() < 30) thisDist *= 2;

	  int nhits = 0;
          int ngoodhits = 0;

	  // are the cluster hits close?
          for (size_t kk = 0; kk < cls_hitlist.size(); ++kk) {
            nhits++;
            int isGoodHit = goodHit(cls_hitlist[kk], thisDist, thisMin, trk_wire1, trk_tick1, trk_wire2, trk_tick2);
            if (isGoodHit == -1){
              ngoodhits = 0;
              break;
            }
            else if (isGoodHit == 1) {
              ngoodhits++;
            }
          } // loop over cluster hits 

	  double fracGood = (double)ngoodhits/nhits;  

	  bool isGoodTrack = fracGood > 0.4;

          if (isGoodTrack) {
            for (size_t kk = 0; kk < cls_hitlist.size(); ++kk) {
	      if ( addShowerHit(cls_hitlist[kk], showerHits) ) showerHits.push_back(cls_hitlist[kk]);
            } // loop over hits to add them to showe
	  }

	} // loop through clusterlist

      } // decide if shower

      /*
      // get dE/dx

      //      auto vhit = fmthm.at(tracklist[i].key());
      //      auto vmeta = fmthm.data(tracklist[i].key());

      TVector3 dir = trkPt2-trkStart;
      dir = dir.Unit();

      unsigned int bestplanetemp = 0;
      //      double minpitch = 999;

      for (unsigned int plane = 0; plane < geom->MaxPlanes(); ++plane) {
	std::vector<float> vQ;
	double pitch = 0;
	double totQ = 0;
	double avgT = 0;
	int nhits = 0;

	for (size_t h = 0; h < vhit.size(); ++h) {
	  unsigned int thisplane = vhit[h]->WireID().planeID().Plane;
	  if (thisplane != plane) continue;

	  if (!pitch) { // find pitch if it hasn't been calculated 
	    double wirePitch = geom->WirePitch(vhit[h]->WireID().planeID());
	    double angleToVert = geom->WireAngleToVertical(geom->Plane(vhit[h]->WireID().planeID()).View(), vhit[h]->WireID().planeID()) - 0.5 * ::util::pi<>();

	    double cosgamma = std::abs(std::sin(angleToVert) * dir[1] + std::cos(angleToVert) * dir[2] );
	    if (cosgamma > 0) pitch = wirePitch/cosgamma;

	    if (pitch < minpitch) {
	      minpitch = pitch;
	      bestplanetemp = plane;
	    }
	  } // calculate pit h
	  double x = tracklist[i]->TrajectoryPoint(vmeta[h]->Index()).position.X();
	  double y = tracklist[i]->TrajectoryPoint(vmeta[h]->Index()).position.Y();
	  double z = tracklist[i]->TrajectoryPoint(vmeta[h]->Index()).position.Z();

	  double x0 = tracklist[i]->Vertex().X();
	  double y0 = tracklist[i]->Vertex().Y();
	  double z0 = tracklist[i]->Vertex().Z();

	  double dist = sqrt( pow(x-x0,2) + pow(y-y0,2) + pow(z-z0,2) );
	  if (dist > 5) continue;

	  vQ.push_back(vhit[h]->Integral());
	  totQ += vhit[h]->Integral();
	  avgT += vhit[h]->PeakTime();
	  ++nhits;

	} // loop through hits
	if (totQ) {
	  double dQdx = TMath::Median(vQ.size(), &vQ[0])/pitch;
	  dEdx[plane] = fCalorimetryAlg.dEdx_AREA(dQdx, avgT/nhits, plane);
	}

      } // loop through planes

      if (showerCandidate) {
	shwDir = (trkPt2-trkStart).Unit();
	shwvtx = tracklist[i]->Vertex();
	bestplane = int(bestplanetemp);
	break;
      }
      */
      if (showerCandidate) {
	std::cout << "THIS IS THE SHOWER PFP: " << allpfps[i].pfp->Self() + 1 << std::endl;
	break;
      }

    } // loop through allpfps

    if (showerCandidate) return 1;

    return 0;
  } // makeShowers

  // ----------------------------------------------------- 
  // return -1 if hit is too close to track vertex or has a wide opening angle
  // return 1 if hit is close to the shower axis
  // return 0 otherwise                                     

  int TCShowerAlg::goodHit(art::Ptr<recob::Hit> hit, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2){

    int pull = 0;
    return goodHit(hit, maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2, pull);

  } // goodHit

  // -----------------------------------------------------   
  // return -1 if hit is too close to track vertex or has a wide opening angle
  // return 1 if hit is close to the shower axis
  // return 0 otherwise   
  
  int TCShowerAlg::goodHit(art::Ptr<recob::Hit> hit, double maxDist, double minDistVert, std::map<geo::PlaneID, double> trk_wire1, std::map<geo::PlaneID, double> trk_tick1, std::map<geo::PlaneID, double> trk_wire2, std::map<geo::PlaneID, double> trk_tick2, int& pull){

    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<geo::Geometry> geom;

    double wirePitch = geom->WirePitch(hit->WireID());
    double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    double UnitsPerTick = tickToDist / wirePitch;

    double x0 = hit->WireID().Wire;
    double y0 = hit->PeakTime() * UnitsPerTick;

    double x1 = trk_wire1[hit->WireID()];
    double y1 = trk_tick1[hit->WireID()] * UnitsPerTick;

    double x2 = trk_wire2[hit->WireID()];
    double y2 = trk_tick2[hit->WireID()] * UnitsPerTick;

    double distToVert = std::sqrt( pow(x0 - x1, 2) + pow(y0 - y1, 2) );
    if (distToVert < minDistVert) return -1;

    if (x2 > x1) {
      if (distToVert < 50) maxDist /= 4;
      else if (distToVert < 100) maxDist /= 2; // trying to exclude photon showers
    }

    // exclude cluster if it's "behind" the vertex
    double a = std::sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) );
    double b = std::sqrt( pow(x0 - x1, 2) + pow(y0 - y1, 2) );
    double c = std::sqrt( pow(x0 - x2, 2) + pow(y0 - y2, 2) );
    double costheta = -( pow(c,2) - pow(a,2) - pow(b,2) ) / (2 * a * b);
    if (costheta < 0) return -1;

    double dist = std::abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/std::sqrt( pow((y2-y1), 2) + pow((x2-x1), 2) );

    if (dist < maxDist) {
      if ( ( (y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) > 0) pull = 1;
      else pull = -1;
      return 1;
    }

    return 0;

  } // goodHit  

  // -----------------------------------------------------
  
  bool TCShowerAlg::addShowerHit(art::Ptr<recob::Hit> hit, std::vector< art::Ptr<recob::Hit> > showerhits) {

    for (size_t i = 0; i < showerhits.size(); ++i) {
      if ( hit.key() == showerhits[i].key() ) return false;
    }
    
    return true;
    
  } // addShowerHit  

  // -----------------------------------------------------

} // namespace shower

