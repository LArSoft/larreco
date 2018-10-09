#include "TCShowerAlg.h"

struct pfpStuff {
  art::Ptr<recob::PFParticle> pfp;
  art::Ptr<recob::Vertex> vtx;
  std::vector<art::Ptr<recob::EndPoint2D> > vx2;
  std::vector<art::Ptr<recob::Hit> > hits;
  double score;
};

bool comparePFP(const pfpStuff& l, const pfpStuff& r) {

  art::Ptr<recob::Vertex> lvtx = l.vtx;
  art::Ptr<recob::Vertex> rvtx = r.vtx;

  double lz = l.hits.size();
  double rz = r.hits.size();
   
  int hitthres = 50; // TODO: ADJUST THIS THRESHOLD 

  if (lz > hitthres && rz <= hitthres) return false;
  else if (lz <= hitthres && rz > hitthres) return true;
  return lvtx->position().Z() > rvtx->position().Z();

  //  return l.score > r.score;
}

bool compareHit(const art::Ptr<recob::Hit>& l, const art::Ptr<recob::Hit>& r) {

  // TODO: TRY FLIPPING HIT SORTING

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

  int TCShowerAlg::makeShowers(std::vector<art::Ptr<recob::PFParticle> > pfplist, std::vector<art::Ptr<recob::Vertex> > vertexlist, std::vector<art::Ptr<recob::Cluster> > clusterlist, std::vector<art::Ptr<recob::Hit> > hitlist, std::vector<art::Ptr<recob::EndPoint2D> > vx2list, art::FindManyP<recob::Hit> cls_fm, art::FindManyP<recob::Cluster> clspfp_fm, art::FindManyP<recob::Vertex> vtxpfp_fm, art::FindManyP<recob::PFParticle> hit_fm, art::FindManyP<recob::Cluster> hitcls_fm, art::FindManyP<anab::Calorimetry> fmcal) {

    totalEnergy.resize(2);
    totalEnergyErr.resize(2);
    dEdx.resize(2);
    dEdxErr.resize(2);   

    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<pfpStuff> allpfps;

    /*
    for (size_t i = 0; i < vx2list.size(); ++i) {
      std::cout << "vx2 " << vx2list[i]->ID() << " strength " << vx2list[i]->Strength() << "\t" << vx2list[i]->View() << " " << vx2list[i]->WireID().asWireID().Wire << ":" << vx2list[i]->DriftTime() << std::endl;
    } // loop through 2d vertices
    */
    // put together pfparticle information
    for (size_t i = 0; i < pfplist.size(); ++i) {
      pfpStuff thispfp;
      thispfp.hits.clear();
      thispfp.pfp = pfplist[i];
      
      thispfp.vtx = vtxpfp_fm.at(pfplist[i].key())[0];
      std::vector<art::Ptr<recob::Cluster> > thisclusterlist = clspfp_fm.at(pfplist[i].key());

      std::vector<int> clustersize;
      std::vector<double> clusterscore;

      for (size_t j = 0; j < thisclusterlist.size(); ++j) {
	double thisclusterscore = 999;
	int vx2Index = -1;
	
	std::vector<art::Ptr<recob::Hit> > thishitlist = cls_fm.at(thisclusterlist[j].key());

	clustersize.push_back((int)thishitlist.size());

	for (size_t k = 0; k < thishitlist.size(); ++k) {
	  thispfp.hits.push_back(thishitlist[k]);
	  //	  if (thishitlist.size() > 10) thispfp.hits.push_back(thishitlist[k]); // NEED TO TEST

	  // loop through vx2s to see which one is closest to the hit
	  // see if doca is smaller than the current best

	  for (size_t l = 0; l < vx2list.size(); ++l) {
	    if (thishitlist[k]->View() != vx2list[l]->View() ) continue;

	    double wirePitch = geom->WirePitch(thishitlist[k]->WireID());
	    double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
	    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
	    double UnitsPerTick = tickToDist / wirePitch;
	    
	    double x1 = thishitlist[k]->WireID().Wire;
	    double y1 = thishitlist[k]->PeakTime() * UnitsPerTick;

	    double x2 = vx2list[l]->WireID().Wire;
	    double y2 = vx2list[l]->DriftTime() * UnitsPerTick;

	    double dist = std::sqrt(pow(x1-x2, 2) + pow(y1-y2, 2) );

	    if (dist < thisclusterscore) {
	      thisclusterscore = dist;
	      vx2Index = l;
	      // std::cout << "dist " << dist << " vtx " << vx2list[l]->ID() << " pfp " << pfplist[i]->Self()+1 << std::endl;
	    }

	  }

	} // loop through hits

	clusterscore.push_back(thisclusterscore);
	if (vx2Index != -1) thispfp.vx2.push_back(vx2list[vx2Index]);

      } // loop through clusters

      //      std::cout << "pfp " << thispfp.pfp->Self() + 1 << " nclusters " << clustersize.size() << " nvx2 " << thispfp.vx2.size() << std::endl;  

      //      if (clustersize.size() == 3 && thispfp.vx2.size() == 3) {
      if (clustersize.size() == 3) {
	
	thispfp.score = clusterscore[0] + clusterscore[1] + clusterscore[2];

	allpfps.push_back(thispfp);
	
	double tick = detprop->ConvertXToTicks(thispfp.vtx->position().X(), geo::PlaneID(0,0,2) );
        int wire = geom->WireCoordinate(thispfp.vtx->position().Y(), thispfp.vtx->position().Z(), geo::PlaneID(0,0,2));

	std::cout << "pfp " << thispfp.pfp->Self() + 1 << " cluster sizes " << clustersize[0] << ":" << clustersize[1] << ":" << clustersize[2] << " pfp score " << thispfp.score << " vertex " << thispfp.vtx->ID() << " " << tick << ":" << wire << " z " << thispfp.vtx->position().Z() << std::endl;

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
      std::vector<art::Ptr<recob::Hit> > pfphits = allpfps[i].hits;
      std::vector<art::Ptr<recob::Cluster> > pfpcls = clspfp_fm.at(allpfps[i].pfp.key());

      //      if (pfphits.size() == 0) continue;

      std::cout << "pfp " << allpfps[i].pfp->Self() + 1 << " hits " << pfphits.size() << std::endl;

      pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(pfphits);
      std::vector<TVector3> spts;
      for (size_t j = 0; j < pmatrack->size(); ++j) {
	if ( (*pmatrack)[j]->IsEnabled() ) {
	  TVector3 p3d = (*pmatrack)[j]->Point3D();
	  spts.push_back(p3d);
	  //	  std::cout << "track point " << p3d[0] << " " << p3d[1] << " " << p3d[2] << std::endl;
	}
      } // pmatrack points

      TVector3 vtx;
      vtx[0] = pfpvtx->position().X();
      vtx[1] = pfpvtx->position().Y();
      vtx[2] = pfpvtx->position().Z();

      if (spts.size() >= 2) {
	//	if ( (spts[0]-vtx).Mag() < (spts.back()-vtx).Mag() ) { // choose direction based on pfparticle vertex
	if ( (spts[0])[2] < spts.back()[2] ) { // choose direction based on z
	  shwvtx = spts[0];
	  size_t sp = 5;
	  if (spts.size() - 1 < sp) sp = spts.size() - 1;
	  shwDir = spts[sp] - spts[0];
	  shwDir = shwDir.Unit();
	}
	else {
	  shwvtx = spts.back();
	  size_t sp = 0;
	  if (spts.size() > 6) sp = spts.size() - 6;
	  shwDir = spts[sp] - spts[spts.size() - 1];
	  shwDir = shwDir.Unit();
	}
      }

      //      int tolerance = 100; // how many shower like cluster you need to define a shower              
      int tolerance = 60; // how many shower like cluster you need to define a shower              
      double pullTolerance = 0.7; // hits should be evenly distributed around the track
      double maxDist = 10; // how far a shower like cluster can be from the track
      double minDistVert = 15; // exclude tracks near the vertex

      if (pfphits.size() < 30) continue;
      //      if (pfphits.size() < 15) continue;
      if (pfphits.size() > 500) continue;  
      // adjust tolerances for short tracks
      if (pfphits.size() < 90) {
	tolerance = 50;
	pullTolerance = 0.9;
      }

      // add pfp hits to shower
      for (size_t ii = 0; ii < pfphits.size(); ++ii) {
	if ( addShowerHit(pfphits[ii], showerHits) ) showerHits.push_back(pfphits[ii]);
      } // loop over pfphits

      int nShowerHits = 0;
      double showerHitPull = 0;

      TVector3 pfpStart = shwvtx;
      TVector3 pfpPt2 = shwvtx+shwDir; // a second point along the track

      /*
      pfpStart[0] = pfpvtx->position().X();
      pfpStart[1] = pfpvtx->position().Y();
      pfpStart[2] = pfpvtx->position().Z();
      */

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
      /*
      bool clusterTooSmall = false;

      for (size_t ii = 0; ii < pfpcls.size(); ++ii) {
	std::vector<art::Ptr<recob::Hit> > clshitlist = cls_fm.at(pfpcls[ii].key() );

	int firstIndex = 2;
	int secondIndex = 9;

	if (clshitlist.size() < 11) {
	  //	  clusterTooSmall = true;
	  //	  break;
	  secondIndex = clshitlist.size() - 1;
	}
	if (clshitlist.size() < 4) {
	  firstIndex = 0;
	  secondIndex = clshitlist.size() - 1;
	}

	std::sort(clshitlist.begin(), clshitlist.end(), compareHit);
	std::reverse(clshitlist.begin(), clshitlist.end());

	// need to decide a scheme for flipping hits on a plane
	//	bool flipHits = false;

	// flip according to 2D vertices
	for (size_t jj = 0; jj < allpfps[i].vx2.size(); ++jj) {

	  art::Ptr<recob::EndPoint2D> thisvx2 = allpfps[i].vx2[jj];

	  if (thisvx2->View() != pfpcls[ii]->View() ) continue;

	  int w1 = thisvx2->WireID().Wire;
	  int w2 = clshitlist[0]->WireID().Wire;
	  int w3 = clshitlist[clshitlist.size()-2]->WireID().Wire;

	  std::cout << thisvx2->View() << " "  << w1 << " " << w2 << " " << w3 << std::endl;

	  int d1 = abs(w1-w2);
	  int d2 = abs(w1-w3);
	  if (d2 < d1) flipHits = true;

	} // pfp.vx2

	for (auto iPlane = geom->begin_plane_id(); iPlane != geom->end_plane_id(); ++iPlane) {
	  if (pfpcls[ii]->Plane() != *iPlane) continue;

	  double wirePitch = geom->WirePitch(clshitlist[0]->WireID());
	  double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
	  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns                            
	  double UnitsPerTick = tickToDist / wirePitch;

	  //	  double vtxtick = detprop->ConvertXToTicks(allpfps[i].vtx->position().X(), *iPlane );
	  //	  int vtxwire = geom->WireCoordinate(allpfps[i].vtx->position().Y(), allpfps[i].vtx->position().Z(), *iPlane);
	  double vtxtick = detprop->ConvertXToTicks(shwvtx[0], *iPlane );
	  int vtxwire = geom->WireCoordinate(shwvtx[1], shwvtx[2], *iPlane);

	  double x0 = (double) vtxwire;
	  double y0 = vtxtick * UnitsPerTick;

          double x1 = clshitlist[0]->WireID().Wire;
	  double y1 = clshitlist[0]->PeakTime() * UnitsPerTick;

          double x2 = clshitlist[clshitlist.size()-1]->WireID().Wire;
	  double y2 = clshitlist[clshitlist.size()-1]->PeakTime() * UnitsPerTick;

	  double d1 =  sqrt( pow(x0-x1,2) + pow(y0-y1, 2) );
	  double d2 =  sqrt( pow(x0-x2,2) + pow(y0-y2, 2) );

	  std::cout << "d1 " << d1 << " d2 " << d2 << std::endl;

	  if (d1 > d2) flipHits = true;

	}
	std::cout << "FLIP HITS " << flipHits << std::endl; 
	//	if (flipHits) std::reverse(clshitlist.begin(), clshitlist.end());

	// check if need to reverse shower!
	auto iPlane = pfpcls[ii]->Plane();

	trk_tick1[iPlane] = clshitlist[firstIndex]->PeakTime();
	trk_wire1[iPlane] = clshitlist[firstIndex]->WireID().asWireID().Wire;
	trk_tick2[iPlane] = clshitlist[secondIndex]->PeakTime();
	trk_wire2[iPlane] = clshitlist[secondIndex]->WireID().asWireID().Wire;


      } // pfpcls

      if (clusterTooSmall) continue;
      */

      for (size_t j = 0; j < clusterlist.size(); ++j) {
	std::vector< art::Ptr<recob::Hit> > cls_hitlist = cls_fm.at(clusterlist[j].key());

	if (clusterlist[j]->ID() > 0 && cls_hitlist.size() > 10) continue;

	bool isGoodCluster = false; // true if the hit belongs to a cluster that should be added to the shower

	for (size_t jj = 0; jj < cls_hitlist.size(); ++jj) {
	  // don't count hits associated with the track
	  // TODO: get PFParticle from hit and check if ID matches parent
	  /*
	  std::vector< art::Ptr<recob::Track> > hit_trklist = hit_fm.at(cls_hitlist[jj].key());
	  if (hit_trklist.size()) {
	    if (hit_trklist[0]->ID() == tracklist[i]->ID()) continue;
	  }
	  */
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

	// loop over hits to find those that aren't associated with any clusters
	if (nShowerHits > 400) maxDist *= 2; // TODO: optimize this threshold              
	for (size_t k = 0; k < hitlist.size(); ++k) {
	  std::vector< art::Ptr<recob::Cluster> > hit_clslist = hitcls_fm.at(hitlist[k].key());
	  if (hit_clslist.size()) continue;
	  int isGoodHit = goodHit(hitlist[k], maxDist*2, minDistVert*2, trk_wire1, trk_tick1, trk_wire2, trk_tick2);
	  if (isGoodHit == 1 && addShowerHit(hitlist[k], showerHits) ) showerHits.push_back(hitlist[k]);
	} // loop over hits

	/*
	TVector3 parentDir = trkPt2 - trkStart;

	// loop over tracks to see if any fall within the shower
	for (size_t k = 0; k < tracklist.size(); ++k) {
	  if (k == i) continue; // don't check current track  

	  TVector3 trkStartOther = tracklist[k]->Vertex();
	  TVector3 trkPt2Other; // a second point along the track

	  recob::Track::Point_t trkPt2temp  = tracklist[k]->TrajectoryPoint(5).position;
	  trkPt2Other[0] = trkPt2temp.X();
	  trkPt2Other[1] = trkPt2temp.Y();
	  trkPt2Other[2] = trkPt2temp.Z();

	  TVector3 otherDir = trkPt2Other - trkStartOther;

	  double ang = parentDir.Angle(otherDir);

	  std::vector< art::Ptr<recob::Hit> > trk_hitlistOther = trk_fm.at(tracklist[k].key());

	  int nhits = 0;
	  int ngoodhits = 0;

	  maxDist = 10;
	  if (tracklist[k]->Length() > 30) maxDist /= 2;

	  // are the track hits close?
	  for (size_t kk = 0; kk < trk_hitlistOther.size(); ++kk) {
	    nhits++;
	    int isGoodHit = goodHit(trk_hitlist[kk], maxDist, minDistVert, trk_wire1, trk_tick1, trk_wire2, trk_tick2);

	    if (isGoodHit == -1){
	      ngoodhits = 0;
	      break;
	    }
	    else if (isGoodHit == 1) {
	      ngoodhits++;
	    }
	  } // loop over track hits  

	  double fracGood = (double)ngoodhits/nhits;

	  double x1 = tracklist[i]->End().X();
	  double y1 = tracklist[i]->End().Y();
	  double z1 = tracklist[i]->End().Z();

	  double x2 = tracklist[k]->End().X();
	  double y2 = tracklist[k]->End().Y();
	  double z2 = tracklist[k]->End().Z();

	  double dist = sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2) );

	  if (dist > 20)
	    if (z2 - z1 < 0) continue;

	  bool isGoodTrack = (ang < 0.5 && fracGood > 0.4 && nhits >= 50);
	  if (ang < 1 && fracGood > 0.5 && nhits < 50 && nhits > 20 ) isGoodTrack = true;
	  if (ang < 1 && fracGood > 0.6 && nhits <= 20) isGoodTrack = true;
	  if (tracklist[k]->Length() > 30) isGoodTrack = (ang < 0.5 && fracGood > 0.7);

	  if (nhits > 80 && ang > 0.2) isGoodTrack = false;

	  if (isGoodTrack) {

	    for (size_t kk = 0; kk < trk_hitlistOther.size(); ++kk) {
	      if ( addShowerHit(trk_hitlistOther[kk], showerHits) ) showerHits.push_back(trk_hitlistOther[kk]);
	    } // loop over hits to add them to shower
	  } // good track
	} // loop over tracks  
	*/
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

