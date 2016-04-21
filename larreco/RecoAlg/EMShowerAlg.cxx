///////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/EMShowerAlg.h"

shower::EMShowerAlg::EMShowerAlg(fhicl::ParameterSet const& pset) : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
								    fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
								    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
								    fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")) {

  fMinTrackLength  = pset.get<double>("MinTrackLength");
  fdEdxTrackLength = pset.get<double>("dEdxTrackLength");
  fNfitpass        = pset.get<unsigned int>("Nfitpass");
  fNfithits        = pset.get<std::vector<unsigned int> >("Nfithits");
  fToler           = pset.get<std::vector<double> >("Toler");
  if (fNfitpass!=fNfithits.size() ||
      fNfitpass!=fToler.size()) {
    throw art::Exception(art::errors::Configuration)
      << "EMShowerAlg: fNfithits and fToler need to have size fNfitpass";
  }
  fDebug = pset.get<int>("Debug",0);
  fDetector = pset.get<std::string>("Detector","dune35t");

  hTrueDirection = tfs->make<TH1I>("trueDir","",2,0,2);

  //this->MakePicture();

}

void shower::EMShowerAlg::AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
						     art::FindManyP<recob::Hit> const& fmh,
						     art::FindManyP<recob::Track> const& fmt,
						     std::map<int,std::vector<int> >& clusterToTracks,
						     std::map<int,std::vector<int> >& trackToClusters) {

  std::vector<int> clustersToIgnore = {-999};
  this->AssociateClustersAndTracks(clusters, fmh, fmt, clustersToIgnore, clusterToTracks, trackToClusters);

  return;

}

void shower::EMShowerAlg::AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
						     art::FindManyP<recob::Hit> const& fmh,
						     art::FindManyP<recob::Track> const& fmt,
						     std::vector<int> const& clustersToIgnore,
						     std::map<int,std::vector<int> >& clusterToTracks,
						     std::map<int,std::vector<int> >& trackToClusters) {

  // Look through all the clusters
  for (std::vector<art::Ptr<recob::Cluster> >::const_iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt) {

    // Get the hits in this cluster
    std::vector<art::Ptr<recob::Hit> > clusterHits = fmh.at(clusterIt->key());

    // Look at all these hits and find the associated tracks
    for (std::vector<art::Ptr<recob::Hit> >::iterator clusterHitIt = clusterHits.begin(); clusterHitIt != clusterHits.end(); ++clusterHitIt) {

      // Get the tracks associated with this hit
      std::vector<art::Ptr<recob::Track> > clusterHitTracks = fmt.at(clusterHitIt->key());
      if (clusterHitTracks.size() > 1) { std::cout << "More than one track associated with this hit!" << std::endl; continue; }
      if (clusterHitTracks.size() < 1) continue;
      if (clusterHitTracks.at(0)->Length() < fMinTrackLength) {
	if (fDebug > 1)
	  std::cout << "Track " << clusterHitTracks.at(0)->ID() << " is too short! (" << clusterHitTracks.at(0)->Length() << ")" << std::endl;
	continue;
      }

      // Add this cluster to the track map
      int track = clusterHitTracks.at(0).key();
      //int trackID = clusterHitTracks.at(0)->ID();
      int cluster = (*clusterIt).key();
      if (std::find(clustersToIgnore.begin(), clustersToIgnore.end(), cluster) != clustersToIgnore.end())
	continue;
      if (std::find(trackToClusters[track].begin(), trackToClusters[track].end(), cluster) == trackToClusters[track].end())
	trackToClusters[track].push_back(cluster);
      if (std::find(clusterToTracks[cluster].begin(), clusterToTracks[cluster].end(), track) == clusterToTracks[cluster].end())
	clusterToTracks[cluster].push_back(track);

    }

  }

  return;

}

void shower::EMShowerAlg::CheckIsolatedHits(std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap) {

  std::map<int,std::vector<int> > firstTPC;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    firstTPC[showerHitsIt->second.at(0)->WireID().TPC].push_back(showerHitsIt->first);

  // If all in the same TPC then that's great!
  if (firstTPC.size() == 1)
    return;

  // If they are in more than two TPCs, not much we can do
  else if (firstTPC.size() > 2)
    return;

  // If we get to this point, there should be something we can do!

  // Find the problem plane
  int problemPlane = -1;
  for (std::map<int,std::vector<int> >::iterator firstTPCIt = firstTPC.begin(); firstTPCIt != firstTPC.end(); ++firstTPCIt)
    if (firstTPCIt->second.size() == 1)
      problemPlane = firstTPCIt->second.at(0);

  // Require three hits
  if (showerHitsMap.at(problemPlane).size() < 3)
    return;

  // and get the other planes with at least three hits
  std::vector<int> otherPlanes;
  for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane)
    if (plane != problemPlane and
	showerHitsMap.count(plane) and
	showerHitsMap.at(plane).size() >= 3)
      otherPlanes.push_back(plane);

  if (otherPlanes.size() == 0)
    return;

  // Look at the hits after the first one
  if (showerHitsMap.at(problemPlane).at(0)->WireID().TPC == showerHitsMap.at(problemPlane).at(1)->WireID().TPC)
    return;

  std::map<int,int> tpcCount;
  for (std::vector<int>::iterator otherPlaneIt = otherPlanes.begin(); otherPlaneIt != otherPlanes.end(); ++otherPlaneIt)
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = std::next(showerHitsMap.at(*otherPlaneIt).begin());
	 hitIt != showerHitsMap.at(*otherPlaneIt).end() and
	   hitIt != std::next(showerHitsMap.at(*otherPlaneIt).begin(),2);
	 ++hitIt)
      ++tpcCount[(*hitIt)->WireID().TPC];

  // Remove the first hit if it is in the wrong TPC
  if (tpcCount.size() == 1 and tpcCount.begin()->first == (int)showerHitsMap.at(problemPlane).at(1)->WireID().TPC) {
    art::Ptr<recob::Hit> naughty_hit = showerHitsMap.at(problemPlane).at(0);
    showerHitsMap.at(problemPlane).erase(showerHitsMap.at(problemPlane).begin());
    showerHitsMap.at(problemPlane).push_back(naughty_hit);
  }

  return;

}

bool shower::EMShowerAlg::CheckShowerHits(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) {

  bool consistencyCheck = true;

  if (showerHitsMap.size() < 2)
    consistencyCheck = true;

  else if (showerHitsMap.size() == 2) {

    // With two views, we can check:
    //  -- timing between views is consistent
    //  -- the 3D start point makes sense when projected back onto the individual planes

    std::vector<art::Ptr<recob::Hit> > startHits;
    std::vector<geo::PlaneID> planes;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      startHits.push_back(showerHitsIt->second.front());
      planes.push_back(showerHitsIt->second.front()->WireID().planeID());
    }

    TVector3 showerStartPos = Construct3DPoint(startHits.at(0), startHits.at(1));
    TVector2 proj1 = Project3DPointOntoPlane(showerStartPos, planes.at(0));
    TVector2 proj2 = Project3DPointOntoPlane(showerStartPos, planes.at(1));

    double timingDifference = TMath::Abs( startHits.at(0)->PeakTime() - startHits.at(1)->PeakTime() );
    double projectionDifference = ( (HitPosition(startHits.at(0)) - proj1).Mod() + (HitPosition(startHits.at(1)) - proj2).Mod() ) / (double)2;

    if (timingDifference > 40 or
	projectionDifference > 1 or
	showerStartPos.X() == -9999 or showerStartPos.Y() == -9999 or showerStartPos.Z() == -9999)
      consistencyCheck = false;

    if (fDebug > 0)
      std::cout << "Timing difference is " << timingDifference << " and projection distance is " << projectionDifference << " (start is (" << showerStartPos.X() << ", " << showerStartPos.Y() << ", " << showerStartPos.Z() << ")" << std::endl;

  }

  else if (showerHitsMap.size() == 3) {

    // With three views, we can check:
    //  -- the timing between views is consistent
    //  -- the 3D start point formed by two views and projected back into the third is close to the start point in that view

    std::map<int,art::Ptr<recob::Hit> > start2DMap;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitIt = showerHitsMap.begin(); showerHitIt != showerHitsMap.end(); ++showerHitIt)
      start2DMap[showerHitIt->first] = showerHitIt->second.front();

    std::map<int,double> projDiff;
    std::map<int,double> timingDiff;

    for (int plane = 0; plane < 3; ++plane) {

      std::vector<int> otherPlanes;
      for (int otherPlane = 0; otherPlane < 3; ++otherPlane)
	if (otherPlane != plane)
	  otherPlanes.push_back(otherPlane);

      TVector3 showerStartPos = Construct3DPoint(start2DMap.at(otherPlanes.at(0)), start2DMap.at(otherPlanes.at(1)));
      TVector2 showerStartProj = Project3DPointOntoPlane(showerStartPos, start2DMap.at(plane)->WireID().planeID());

      if (fDebug > 1) {
	std::cout << "Plane... " << plane << std::endl;
	std::cout << "Start position in this plane is " << HitPosition(start2DMap.at(plane)).X() << ", " << HitPosition(start2DMap.at(plane)).Y() << ")" << std::endl;
	std::cout << "Shower start from other two planes is (" << showerStartPos.X() << ", " << showerStartPos.Y() << ", " << showerStartPos.Z() << ")" << std::endl;
	std::cout << "Projecting the other two planes gives position (" << showerStartProj.X() << ", " << showerStartProj.Y() << ")" << std::endl;
      }

      double projDiff = TMath::Abs((showerStartProj-HitPosition(start2DMap.at(plane))).Mod());
      double timeDiff = TMath::Max(TMath::Abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(0))->PeakTime()),
				   TMath::Abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(1))->PeakTime()));

      if (fDebug > 0)
	std::cout << "Plane " << plane << " has projDiff " << projDiff << " and timeDiff " << timeDiff << std::endl;

      if (projDiff > 1 or timeDiff > 40)
	consistencyCheck = false;

    }

  }

  if (fDebug > 1)
    std::cout << "Consistency check is " << consistencyCheck << std::endl;

  return consistencyCheck;

}

std::vector<int> shower::EMShowerAlg::CheckShowerPlanes(std::vector<std::vector<int> > const& initialShowers,
							std::vector<art::Ptr<recob::Cluster> > const& clusters,
							art::FindManyP<recob::Hit> const& fmh) {

  std::vector<int> clustersToIgnore;

  // // Look at each shower
  // for (std::vector<std::vector<int> >::const_iterator initialShowerIt = initialShowers.begin(); initialShowerIt != initialShowers.end(); ++initialShowerIt) {

  //   // Map the clusters and cluster hits to each view
  //   std::map<int,std::vector<art::Ptr<recob::Cluster> > > planeClusters;
  //   std::map<int,std::vector<art::Ptr<recob::Hit> > > planeClusterHits;
  //   for (std::vector<int>::const_iterator clusterIt = initialShowerIt->begin(); clusterIt != initialShowerIt->end(); ++clusterIt) {
  //     art::Ptr<recob::Cluster> cluster = clusters.at(*clusterIt);
  //     std::vector<art::Ptr<recob::Hit> > hits = fmh.at(cluster.key());
  //     planeClusters[cluster->Plane().Plane].push_back(cluster);
  //     for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
  // 	planeClusterHits[cluster.key()].push_back(*hitIt);
  //   }

  //   // Can't do much with fewer than three views
  //   if (planeClusters.size() < 3)
  //     continue;

  //   // Look at the average RMS of clusters in each view
  //   std::map<int,double> avRMS;
  //   for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeClusterIt = planeClusters.begin(); planeClusterIt != planeClusters.end(); ++planeClusterIt) {
  //     std::cout << "Plane " << planeClusterIt->first << std::endl;
  //     for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeClusterIt->second.begin(); clusterIt != planeClusterIt->second.end(); ++clusterIt) {
  // 	double rms = ShowerHitRMS(planeClusterHits.at(clusterIt->key()));
  // 	std::cout << "Cluster " << clusterIt->key() << " has RMS " << rms << std::endl;
  //     }
  //   }
  // }


  // Look at each shower
  for (std::vector<std::vector<int> >::const_iterator initialShowerIt = initialShowers.begin(); initialShowerIt != initialShowers.end(); ++initialShowerIt) {

    if (std::distance(initialShowers.begin(),initialShowerIt) > 0)
      continue;

    // Map the clusters and cluster hits to each view
    std::map<int,std::vector<art::Ptr<recob::Cluster> > > planeClusters;
    std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHits;
    for (std::vector<int>::const_iterator clusterIt = initialShowerIt->begin(); clusterIt != initialShowerIt->end(); ++clusterIt) {
      art::Ptr<recob::Cluster> cluster = clusters.at(*clusterIt);
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(cluster.key());
      planeClusters[cluster->Plane().Plane].push_back(cluster);
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
  	planeHits[(*hitIt)->WireID().Plane].push_back(*hitIt);
    }

    TFile* outFile = new TFile("chargeDistributions.root","RECREATE");
    std::map<int,TH1D*> chargeDist;
    for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeIt = planeClusters.begin(); planeIt != planeClusters.end(); ++planeIt) {
      for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeIt->second.begin(); clusterIt != planeIt->second.end(); ++clusterIt) {
	chargeDist[planeIt->first] = new TH1D(std::string("chargeDist_Plane"+std::to_string(planeIt->first)+"_Cluster"+std::to_string(clusterIt->key())).c_str(),"",150,0,1000);
	std::vector<art::Ptr<recob::Hit> > hits = fmh.at(clusterIt->key());
	for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
	  chargeDist[planeIt->first]->Fill((*hitIt)->Integral());
	outFile->cd();
	chargeDist[planeIt->first]->Write();
      }
    }
    outFile->Close();
    delete outFile;

    // Can't do much with fewer than three views
    if (planeClusters.size() < 3)
      continue;

    // Look at how many clusters each plane has, and the proportion of hits each one uses
    std::map<int,std::vector<double> > planeClusterSizes;
    for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeClustersIt = planeClusters.begin(); planeClustersIt != planeClusters.end(); ++planeClustersIt) {
      for (std::vector<art::Ptr<recob::Cluster> >::iterator planeClusterIt = planeClustersIt->second.begin(); planeClusterIt != planeClustersIt->second.end(); ++planeClusterIt) {
  	std::vector<art::Ptr<recob::Hit> > hits = fmh.at(planeClusterIt->key());
        planeClusterSizes[planeClustersIt->first].push_back((double)hits.size()/(double)planeHits.at(planeClustersIt->first).size());
      }
    }

    // Find the average hit fraction across all clusters in the plane
    std::map<int,double> planeClustersAvSizes;
    for (std::map<int,std::vector<double> >::iterator planeClusterSizesIt = planeClusterSizes.begin(); planeClusterSizesIt != planeClusterSizes.end(); ++planeClusterSizesIt) {
      double average = 0;
      for (std::vector<double>::iterator planeClusterSizeIt = planeClusterSizesIt->second.begin(); planeClusterSizeIt != planeClusterSizesIt->second.end(); ++planeClusterSizeIt)
  	average += *planeClusterSizeIt;
      average /= planeClusterSizesIt->second.size();
      planeClustersAvSizes[planeClusterSizesIt->first] = average;
    }

    // Now decide if there is one plane which is ruining the reconstruction
    // If two planes have a low average cluster fraction and one high, this plane likely merges two particle deposits together
    int badPlane = -1;
    for (std::map<int,double>::iterator clusterAvSizeIt = planeClustersAvSizes.begin(); clusterAvSizeIt != planeClustersAvSizes.end(); ++clusterAvSizeIt) {

      // Get averages from other planes and add in quadrature
      std::vector<double> otherAverages;
      for (std::map<int,double>::iterator otherClustersAvSizeIt = planeClustersAvSizes.begin(); otherClustersAvSizeIt != planeClustersAvSizes.end(); ++otherClustersAvSizeIt)
    	if (clusterAvSizeIt->first != otherClustersAvSizeIt->first)
    	  otherAverages.push_back(otherClustersAvSizeIt->second);
      double quadOtherPlanes = 0;
      for (std::vector<double>::iterator otherAvIt = otherAverages.begin(); otherAvIt != otherAverages.end(); ++otherAvIt)
    	quadOtherPlanes += TMath::Power(*otherAvIt,2);
      quadOtherPlanes = TMath::Sqrt(quadOtherPlanes);

      // If this plane has an average higher than the quadratic sum of the others, it may be bad
      if (clusterAvSizeIt->second >= quadOtherPlanes)
    	badPlane = clusterAvSizeIt->first;

    }

    if (badPlane != -1) {
      if (fDebug > 0)
  	std::cout << "Bad plane is " << badPlane << std::endl;
      for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeClusters.at(badPlane).begin(); clusterIt != planeClusters.at(badPlane).end(); ++clusterIt)
  	clustersToIgnore.push_back(clusterIt->key());
    }

  }

  return clustersToIgnore;

}

TVector3 shower::EMShowerAlg::Construct3DPoint(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2) {

  // x is average of the two x's
  double x = (fDetProp->ConvertTicksToX(hit1->PeakTime(), hit1->WireID().planeID()) + fDetProp->ConvertTicksToX(hit2->PeakTime(), hit2->WireID().planeID())) / (double)2;

  // y and z got from the wire interections
  geo::WireIDIntersection intersection;
  fGeom->WireIDsIntersect(hit1->WireID(), hit2->WireID(), intersection);

  return TVector3(x, intersection.y, intersection.z);

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& hits1,
								  std::vector<art::Ptr<recob::Hit> > const& hits2,
								  std::map<geo::PlaneID,TVector2> const& showerCentreMap) {

  std::unique_ptr<recob::Track> track;

  std::vector<art::Ptr<recob::Hit> > track1, track2;

  // Check the TPCs
  if ((*hits1.begin())->WireID().TPC != (*hits2.begin())->WireID().TPC) {
    mf::LogWarning("EMShowerAlg") << "Warning: attempting to construct a track from two different TPCs.  Returning a null track.";
    return track;
  }

  // Check for tracks crossing TPC boundaries
  std::map<int,int> tpcMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits1.begin(); hitIt != hits1.end(); ++hitIt)
    ++tpcMap[(*hitIt)->WireID().TPC];
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits2.begin(); hitIt != hits2.end(); ++hitIt)
    ++tpcMap[(*hitIt)->WireID().TPC];
  if (tpcMap.size() > 1) {
    mf::LogWarning("EMShowerAlg") << "Warning: attempting to construct a track which crosses more than one TPC -- PMTrack can't handle this right now.  Returning a track made just from hits in the first TPC it traverses.";
    unsigned int firstTPC1 = hits1.at(0)->WireID().TPC, firstTPC2 = hits2.at(0)->WireID().TPC;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits1.begin(); hitIt != hits1.end(); ++hitIt)
      if ((*hitIt)->WireID().TPC == firstTPC1) track1.push_back(*hitIt);
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits2.begin(); hitIt != hits2.end(); ++hitIt)
      if ((*hitIt)->WireID().TPC == firstTPC2) track2.push_back(*hitIt);    
  }
  else {
    track1 = hits1;
    track2 = hits2;
  }

  if (fDebug > 0) {
    std::cout << "About to make me a track from these 'ere 'its... " << std::endl;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit1 = track1.begin(); hit1 != track1.end(); ++hit1)
      std::cout << "Hit (" << HitCoordinates(*hit1).X() << ", " << HitCoordinates(*hit1).Y() << ") (real wire " << (*hit1)->WireID().Wire << ") in TPC " << (*hit1)->WireID().TPC << std::endl;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit2 = track2.begin(); hit2 != track2.end(); ++hit2)
      std::cout << "Hit (" << HitCoordinates(*hit2).X() << ", " << HitCoordinates(*hit2).Y() << ") (real wire " << (*hit2)->WireID().Wire << ") in TPC " << (*hit2)->WireID().TPC << std::endl;
  }

  TVector3 trackStart = Construct3DPoint(track1.at(0), track2.at(0));
  pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(track1, track2, trackStart);

  if (!pmatrack) {
    mf::LogInfo("EMShowerAlg") << "Skipping this event because not enough hits in two views";
    return track;
  }

  std::vector<TVector3> xyz, dircos;
  std::vector<std::vector<double> > dEdx; // Right now, not finding the dE/dx for these tracks.  Can extend if needed.

  for (unsigned int i = 0; i < pmatrack->size(); i++) {

    xyz.push_back((*pmatrack)[i]->Point3D());

    if (i < pmatrack->size()-1) {
      size_t j = i+1;
      double mag = 0.0;
      TVector3 dc(0., 0., 0.);
      while ((mag == 0.0) and (j < pmatrack->size())) {
      	dc = (*pmatrack)[j]->Point3D();
      	dc -= (*pmatrack)[i]->Point3D();
      	mag = dc.Mag();
      	++j;
      }
      if (mag > 0.0) dc *= 1.0 / mag;
      else if (!dircos.empty()) dc = dircos.back();
      // TVector3 dc((*pmatrack)[i+1]->Point3D());
      // dc -= (*pmatrack)[i]->Point3D();
      // dc *= 1.0 / dc.Mag();
      dircos.push_back(dc);
    }
    else dircos.push_back(dircos.back());

  }

  // Orient the track correctly
  std::map<geo::PlaneID,double> distanceToVertex, distanceToEnd;
  TVector3 vertex = *xyz.begin(), end = *xyz.rbegin();

  // Loop over all the planes and find the distance from the vertex and end projections to the centre in each plane
  for (std::map<geo::PlaneID,TVector2>::const_iterator showerCentreIt = showerCentreMap.begin(); showerCentreIt != showerCentreMap.end(); ++showerCentreIt) {

    // Project the vertex and the end point onto this plane
    TVector2 vertexProj = Project3DPointOntoPlane(vertex, showerCentreIt->first);
    TVector2 endProj    = Project3DPointOntoPlane(end, showerCentreIt->first);

    // Find the distance of each to the centre of the cluster
    distanceToVertex[showerCentreIt->first] = (vertexProj - showerCentreIt->second).Mod();
    distanceToEnd[showerCentreIt->first] = (endProj - showerCentreIt->second).Mod();

  }

  // Find the average distance to the vertex and the end across the planes
  double avDistanceToVertex = 0, avDistanceToEnd = 0;
  for (std::map<geo::PlaneID,double>::iterator distanceToVertexIt = distanceToVertex.begin(); distanceToVertexIt != distanceToVertex.end(); ++distanceToVertexIt)
    avDistanceToVertex += distanceToVertexIt->second;
  avDistanceToVertex /= distanceToVertex.size();

  for (std::map<geo::PlaneID,double>::iterator distanceToEndIt = distanceToEnd.begin(); distanceToEndIt != distanceToEnd.end(); ++distanceToEndIt)
    avDistanceToEnd += distanceToEndIt->second;
  avDistanceToEnd /= distanceToEnd.size();

  if (fDebug > 1)
    std::cout << "Distance to vertex is " << avDistanceToVertex << " and distance to end is " << avDistanceToEnd << std::endl;

  // Change order if necessary
  if (avDistanceToEnd > avDistanceToVertex) {
    std::reverse(xyz.begin(), xyz.end());
    std::transform(dircos.begin(), dircos.end(), dircos.begin(), [](TVector3 const& vec){return -1*vec;});
  }

  if (xyz.size() != dircos.size())
    mf::LogError("EMShowerAlg") << "Problem converting pma::Track3D to recob::Track";

  track = std::make_unique<recob::Track>(xyz, dircos, dEdx);

  return track;

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
								  std::vector<art::Ptr<recob::Hit> > const& track2) {

  std::map<geo::PlaneID,TVector2> showerCentreMap;

  return this->ConstructTrack(track1, track2, showerCentreMap);

}

double shower::EMShowerAlg::FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, std::unique_ptr<recob::Track> const& track) {

  double totalCharge = 0, totalDistance = 0, avHitTime = 0;
  unsigned int nHits = 0;

  if (!track)
    return -999;

  // Get the pitch
  double pitch = 0;
  try { pitch = lar::utils::TrackPitchInView(*track, trackHits.at(0)->View()); }
  catch(...) { pitch = 0; }

  // Deal with large pitches
  if (pitch > fdEdxTrackLength) {
    double dEdx = fCalorimetryAlg.dEdx_AREA(*trackHits.begin(), pitch);
    return dEdx;
  }

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
    if (totalDistance + pitch < fdEdxTrackLength) {
      totalDistance += pitch;
      totalCharge += (*trackHitIt)->Integral();
      avHitTime += (*trackHitIt)->PeakTime();
      ++nHits;
    }
  }

  avHitTime /= (double)nHits;

  double dQdx = totalCharge / totalDistance;
  double dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avHitTime, trackHits.at(0)->WireID().Plane);

  return dEdx;

}

void shower::EMShowerAlg::FindInitialTrack(art::PtrVector<recob::Hit> const& hits,
					   std::unique_ptr<recob::Track>& initialTrack,
					   std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits, int plane) {

  /// Finding the initial track requires three stages:
  ///  -- put the hits in the correct order in each view
  ///  -- find the initial track-like hits in each view
  ///  -- use these to construct a track

  // First, order the hits into the correct shower order in each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap = OrderShowerHits(hits, plane);

  // Now find the hits belonging to the track
  initialTrackHits = FindShowerStart(showerHitsMap);

  if (fDebug > 0) {
    std::cout << "Here are the initial shower hits... " << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator initialTrackHitsIt = initialTrackHits.begin(); initialTrackHitsIt != initialTrackHits.end(); ++initialTrackHitsIt) {
      std::cout << "  Plane " << initialTrackHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator initialTrackHitIt = initialTrackHitsIt->second.begin(); initialTrackHitIt != initialTrackHitsIt->second.end(); ++initialTrackHitIt)
	std::cout << "    Hit is (" << HitCoordinates(*initialTrackHitIt).X() << " (real hit " << (*initialTrackHitIt)->WireID() << "), " << HitCoordinates(*initialTrackHitIt).Y() << ")" << std::endl;
    }
  }

  // Now we have the track hits -- can make a track!
  initialTrack = MakeInitialTrack(initialTrackHits);

  if (initialTrack and fDebug > 0) {
    std::cout << "The track start is (" << initialTrack->Vertex().X() << ", " << initialTrack->Vertex().Y() << ", " << initialTrack->Vertex().Z() << ")" << std::endl;
    std::cout << "The track direction is (" << initialTrack->VertexDirection().X() << ", " << initialTrack->VertexDirection().Y() << ", " << initialTrack->VertexDirection().Z() << ")" << std::endl;
  }

  // // Fill correct or incorrect direction histogram
  // std::map<int,int> trackHits;
  // for (art::PtrVector<recob::Hit>::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
  //   ++trackHits[FindTrackID(*hitIt)];
  // int trueTrack = -9999;
  // for (std::map<int,int>::iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt)
  //   if (trackHitIt->second/(double)hits.size() > 0.8)
  //     trueTrack = trackHitIt->first;
  // if (trueTrack != -9999) {
  //   const simb::MCParticle* trueParticle = backtracker->TrackIDToParticle(trueTrack);
  //   TVector3 trueStart = trueParticle->Position().Vect();
  //   if (initialTrack) {
  //     TVector3 recoStart = initialTrack->Vertex();
  //     if ((recoStart-trueStart).Mag() < 5)
  // 	hTrueDirection->Fill(1);
  //     else
  // 	hTrueDirection->Fill(0);
  //   }
  // }
  // else
  //   hTrueDirection->Fill(0);

  return;

}

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindOrderOfHits(std::vector<art::Ptr<recob::Hit> > const& hits, bool perpendicular) {

  // Find the charge-weighted centre (in [cm]) of this shower
  TVector2 centre = ShowerCentre(hits);

  // Find a rough shower 'direction'
  TVector2 direction = ShowerDirection(hits);

  if (perpendicular)
    direction = direction.Rotate(TMath::Pi()/2);

  // Find how far each hit (projected onto this axis) is from the centre
  TVector2 pos;
  std::map<double,art::Ptr<recob::Hit> > hitProjection;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    pos = HitPosition(*hit) - centre;
    hitProjection[direction*pos] = *hit;
  }

  // Get a vector of hits in order of the shower
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(hitProjection.begin(), hitProjection.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

  // TGraph* graph = new TGraph();
  // for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = showerHits.begin(); hitIt != showerHits.end(); ++hitIt) {
  //   //std::cout << "Hit at wire " << (*hitIt)->WireID() << " and tick " << (*hitIt)->PeakTime() << " is pos (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
  //   graph->SetPoint(graph->GetN(), HitPosition(*hitIt).X(), HitPosition(*hitIt).Y());
  // }
  // graph->SetMarkerStyle(8);
  // graph->SetMarkerSize(2);
  // TCanvas* canvas = new TCanvas();
  // graph->Draw("AP");
  // TLine line;
  // line.SetLineColor(2);
  // line.DrawLine(centre.X()-1000*direction.X(),centre.Y()-1000*direction.Y(),centre.X()+1000*direction.X(),centre.Y()+1000*direction.Y());
  // canvas->SaveAs("thisCanvas.png");
  // delete canvas; delete graph;

  return showerHits;

}

std::vector<std::vector<int> > shower::EMShowerAlg::FindShowers(std::map<int,std::vector<int> > const& trackToClusters) {

  // Showers are vectors of clusters
  std::vector<std::vector<int> > showers;

  // Loop over all tracks 
  for (std::map<int,std::vector<int> >::const_iterator trackToClusterIt = trackToClusters.begin(); trackToClusterIt != trackToClusters.end(); ++ trackToClusterIt) {

    // Find which showers already made are associated with this track
    std::vector<int> matchingShowers;
    for (unsigned int shower = 0; shower < showers.size(); ++shower)
      for (std::vector<int>::const_iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if ( (std::find(showers.at(shower).begin(), showers.at(shower).end(), *cluster) != showers.at(shower).end()) and
	     (std::find(matchingShowers.begin(), matchingShowers.end(), shower)) == matchingShowers.end() )
	  matchingShowers.push_back(shower);

    // THINK THERE PROBABLY CAN BE MORE THAN ONE!
    // IN FACT, THIS WOULD BE A SUCCESS OF THE SHOWERING METHOD!
    // // Shouldn't be more than one
    // if (matchingShowers.size() > 1)
    //   mf::LogInfo("EMShowerAlg") << "Number of showers this track matches is " << matchingShowers.size() << std::endl;

    // New shower
    if (matchingShowers.size() < 1)
      showers.push_back(trackToClusterIt->second);

    // Add to existing shower
    else {
      for (std::vector<int>::const_iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if (std::find(showers.at(matchingShowers.at(0)).begin(), showers.at(matchingShowers.at(0)).end(), *cluster) == showers.at(matchingShowers.at(0)).end())
	  showers.at(matchingShowers.at(0)).push_back(*cluster);
    }
  }

  return showers;

}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::FindShowerStart(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap) {

  std::map<int,std::vector<art::Ptr<recob::Hit> > > initialHitsMap;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator orderedShowerIt = orderedShowerMap.begin(); orderedShowerIt != orderedShowerMap.end(); ++orderedShowerIt) {

    std::vector<art::Ptr<recob::Hit> > initialHits;
    const std::vector<art::Ptr<recob::Hit> > orderedShower = orderedShowerIt->second;

    // Find if the shower is travelling along ticks or wires
    bool wireDirection = true;
    std::vector<int> wires;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt)
      wires.push_back(std::round(HitCoordinates(*hitIt).X()));
    std::sort(wires.begin(), wires.end());
    if (TMath::Abs(*wires.begin()-std::round(HitCoordinates(*orderedShower.begin()).X())) > 3 and
    	TMath::Abs(*wires.rbegin()-std::round(HitCoordinates(*orderedShower.begin()).X())) > 3)
      wireDirection = false;

    // Deal with showers travelling along wires
    if (wireDirection) {
      bool increasing = HitCoordinates(*orderedShower.rbegin()).X() > HitCoordinates(*orderedShower.begin()).X();
      std::map<int,std::vector<art::Ptr<recob::Hit> > > wireHitMap;
      int multipleWires = 0;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt)
	wireHitMap[std::round(HitCoordinates(*hitIt).X())].push_back(*hitIt);
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt) {
	int wire = std::round(HitCoordinates(*hitIt).X());
	if (wireHitMap[wire].size() > 1) {
	  ++multipleWires;
	  if (multipleWires > 5) break;
	  continue;
	}
	else if ( (increasing and wireHitMap[wire+1].size() > 1 and (wireHitMap[wire+2].size() > 1 or wireHitMap[wire+3].size() > 1)) or
		  (!increasing and wireHitMap[wire-1].size() > 1 and (wireHitMap[wire-2].size() > 1 or wireHitMap[wire-3].size() > 1)) ) {
	  if ( (increasing and (wireHitMap[wire+4].size() < 2 and wireHitMap[wire+5].size() < 2 and wireHitMap[wire+6].size() < 2 and wireHitMap[wire+9].size() > 1)) or
	       (!increasing and (wireHitMap[wire-4].size() < 2 and wireHitMap[wire-5].size() < 2 and wireHitMap[wire-6].size() < 2) and wireHitMap[wire-9].size() > 1) )
	    initialHits.push_back(*hitIt);
	  else
	    break;
	}
	else
	  initialHits.push_back(*hitIt);
      }
      if (!initialHits.size()) initialHits.push_back(*orderedShower.begin());
    }

    // Deal with showers travelling along ticks
    else {
      bool increasing = HitCoordinates(*orderedShower.rbegin()).Y() > HitCoordinates(*orderedShower.begin()).Y();
      std::map<int,std::vector<art::Ptr<recob::Hit> > > tickHitMap;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt)
	tickHitMap[std::round(HitCoordinates(*hitIt).Y())].push_back(*hitIt);
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt) {
	int tick = std::round(HitCoordinates(*hitIt).Y());
	if ( (increasing and (tickHitMap[tick+1].size() or tickHitMap[tick+2].size() or tickHitMap[tick+3].size() or tickHitMap[tick+4].size() or tickHitMap[tick+5].size())) or
	     (!increasing and (tickHitMap[tick-1].size() or tickHitMap[tick-2].size() or tickHitMap[tick-3].size() or tickHitMap[tick-4].size() or tickHitMap[tick-5].size())) )
	  break;
	else
	  initialHits.push_back(*hitIt);
      }
      if (!initialHits.size()) initialHits.push_back(*orderedShower.begin());
    }

    // Need at least two hits to make a track
    if (initialHits.size() == 1 and orderedShower.size() > 2)
      initialHits.push_back(orderedShower.at(1));

    // Quality check -- make sure there isn't a single hit in a different TPC (artefact of tracking failure)
    std::vector<art::Ptr<recob::Hit> > newInitialHits;
    std::map<int,int> tpcHitMap;
    std::vector<int> singleHitTPCs;
    for (std::vector<art::Ptr<recob::Hit> >::iterator initialHitIt = initialHits.begin(); initialHitIt != initialHits.end(); ++initialHitIt)
      ++tpcHitMap[(*initialHitIt)->WireID().TPC];
    for (std::map<int,int>::iterator tpcIt = tpcHitMap.begin(); tpcIt != tpcHitMap.end(); ++tpcIt)
      if (tpcIt->second == 1) singleHitTPCs.push_back(tpcIt->first);
    if (singleHitTPCs.size()) {
      if (fDebug > 1)
	for (std::vector<int>::iterator tpcIt = singleHitTPCs.begin(); tpcIt != singleHitTPCs.end(); ++tpcIt)
	  std::cout << "Removed hits in TPC " << *tpcIt << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator initialHitIt = initialHits.begin(); initialHitIt != initialHits.end(); ++initialHitIt)
	if (std::find(singleHitTPCs.begin(), singleHitTPCs.end(), (*initialHitIt)->WireID().TPC) == singleHitTPCs.end())
	  newInitialHits.push_back(*initialHitIt);
      if (!newInitialHits.size()) newInitialHits.push_back(*orderedShower.begin());
    }
    else
      newInitialHits = initialHits;

    initialHitsMap[orderedShowerIt->first] = newInitialHits;

  }

  return initialHitsMap;

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::MakeInitialTrack(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap) {

  std::unique_ptr<recob::Track> track;

  // Can't do much with just one view
  if (initialHitsMap.size() < 2) {
    mf::LogInfo("EMShowerAlg") << "Only one useful view for this shower; nothing can be done" << std::endl;
    return track;
  }

  std::vector<std::pair<int,int> > initialHitsSize;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator initialHitIt = initialHitsMap.begin(); initialHitIt != initialHitsMap.end(); ++initialHitIt)
    initialHitsSize.push_back(std::make_pair(initialHitIt->first, initialHitIt->second.size()));

  // Sort the planes by number of hits
  std::sort(initialHitsSize.begin(), initialHitsSize.end(), [](std::pair<int,int> const& size1, std::pair<int,int> const& size2) { return size1.second > size2.second; });

  return ConstructTrack(initialHitsMap.at(initialHitsSize.at(0).first), initialHitsMap.at(initialHitsSize.at(1).first));

  // std::map<int,std::vector<art::Ptr<recob::Hit> > > trackLengthMap;
  // for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator initialHitsIt = initialHitsMap.begin(); initialHitsIt != initialHitsMap.end(); ++initialHitsIt)
  //   trackLengthMap[initialHitsIt->second.size()] = initialHitsIt->second;

  // std::vector<std::vector<art::Ptr<recob::Hit> > > longestTracks;
  // std::transform(trackLengthMap.rbegin(), trackLengthMap.rend(), std::back_inserter(longestTracks),
  // 		 [](std::pair<int,std::vector<art::Ptr<recob::Hit> > > const& initialTrackHits) { return initialTrackHits.second; });

  // return ConstructTrack(longestTracks.at(0), longestTracks.at(1));

}

recob::Shower shower::EMShowerAlg::MakeShower(art::PtrVector<recob::Hit> const& hits,
					      std::unique_ptr<recob::Track> const& initialTrack,
					      std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap) {

  //return recob::Shower();

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  int bestPlane = -1;
  unsigned int highestNumberOfHits = 0;
  std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;

  // Look at each of the planes
  for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {

    // If there's hits on this plane...
    if (planeHitsMap.count(plane) != 0) {
      dEdx.push_back(FinddEdx(initialHitsMap.at(plane), initialTrack));
      totalEnergy.push_back(fShowerEnergyAlg.ShowerEnergy(planeHitsMap.at(plane), plane));
      if (planeHitsMap.at(plane).size() > highestNumberOfHits) {
	bestPlane = plane;
	highestNumberOfHits = planeHitsMap.at(plane).size();
      }
    }

    // If not...
    else {
      dEdx.push_back(0);
      totalEnergy.push_back(0);
    }

  }

  TVector3 direction, directionError, showerStart, showerStartError;
  if (initialTrack) {
    direction = initialTrack->VertexDirection();
    showerStart = initialTrack->Vertex();
  }

  if (fDebug > 0) {
    std::cout << "Best plane is " << bestPlane << std::endl;
    std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
    std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
    std::cout << "The shower start is (" << showerStart.X() << ", " << showerStart.Y() << ", " << showerStart.Z() << ")" << std::endl;
    std::cout << "The shower direction is (" << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << std::endl;
  }

  return recob::Shower(direction, directionError, showerStart, showerStartError, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);

}

recob::Shower shower::EMShowerAlg::MakeShower(art::PtrVector<recob::Hit> const& hits,
					      art::Ptr<recob::Vertex> const& vertex,
					      int & iok) {
  
  iok = 1;

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->WireID().Plane].push_back(*hit);

  std::vector<std::vector<art::Ptr<recob::Hit> > > initialTrackHits(3);

  int pl0 = -1;
  int pl1 = -1;
  unsigned maxhits0 = 0;
  unsigned maxhits1 = 0;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {
    
    std::vector<art::Ptr<recob::Hit> > showerHits;
    OrderShowerHits(planeHits->second, showerHits, vertex);
    //if (!isCleanShower(showerHits)) continue;
    FindInitialTrackHits(showerHits, vertex, initialTrackHits[planeHits->first]);
    if ((planeHits->second).size()>maxhits0){
      if (pl0!=-1){
	maxhits1 = maxhits0;
	pl1 = pl0;
      }
      pl0 = planeHits->first;
      maxhits0 = (planeHits->second).size();
    }
    else if ((planeHits->second).size()>maxhits1){
      pl1 = planeHits->first;
      maxhits1 = (planeHits->second).size();
    }

  }
  //std::cout<<pl0<<" "<<pl1<<std::endl;
//  if (pl0!=-1&&pl1!=-1) {
//    pl0 = 1;
//    pl1 = 2;
//  }
  if (pl0!=-1&&pl1!=-1
      &&initialTrackHits[pl0].size()>=2
      &&initialTrackHits[pl1].size()>=2
      &&initialTrackHits[pl0][0]->WireID().TPC==
      initialTrackHits[pl1][0]->WireID().TPC){
    double xyz[3];
    vertex->XYZ(xyz);
    TVector3 vtx(xyz);
//    std::vector<art::Ptr<recob::Hit>> alltrackhits;
//    for (size_t i = 0; i<3; ++i){
//      for (auto const&hit : initialTrackHits[i]){
//	alltrackhits.push_back(hit);
//      }
//    }
    //std::cout<<"vertex "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
    //for (auto const&hit : initialTrackHits[pl0]) std::cout<<*hit<<std::endl;
    //for (auto const&hit : initialTrackHits[pl1]) std::cout<<*hit<<std::endl;
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(initialTrackHits[pl0], initialTrackHits[pl1]);
    //std::cout<<pmatrack->size()<<std::endl;
    //pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(alltrackhits);
    std::vector<TVector3> spts;
    double xshift = pmatrack->GetXShift();
    bool has_shift = (xshift != 0.0);
    for (size_t i = 0; i<pmatrack->size(); ++i){
      if ((*pmatrack)[i]->IsEnabled()){
	TVector3 p3d = (*pmatrack)[i]->Point3D();
	if (has_shift) p3d.SetX(p3d.X() + xshift);
	//std::cout<<p3d.X()<<" "<<p3d.Y()<<" "<<p3d.Z()<<std::endl;
	spts.push_back(p3d);
      }
    }
    if (spts.size()>=2){ //at least two space points
      TVector3 shwxyz, shwxyzerr;
      TVector3 shwdir, shwdirerr;
      std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;
      int bestPlane = pl0;
      double minpitch = 1000;
      std::vector<TVector3> dirs;
      if ((spts[0]-vtx).Mag()<(spts.back()-vtx).Mag()){
	shwxyz = spts[0];
	size_t i = 5;
	if (spts.size()-1<5) i = spts.size()-1;
	shwdir = spts[i] - spts[0];
	shwdir = shwdir.Unit();
      }
      else{
	shwxyz = spts.back();
	size_t i = 0;
	if (spts.size()>6) i = spts.size() - 6;
	shwdir = spts[i] - spts[spts.size()-1];
	shwdir = shwdir.Unit();
      }
      //std::cout<<shwxyz.X()<<" "<<shwxyz.Y()<<" "<<shwxyz.Z()<<std::endl;
      //std::cout<<shwdir.X()<<" "<<shwdir.Y()<<" "<<shwdir.Z()<<std::endl;
      for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {
	if (planeHitsMap.find(plane)!=planeHitsMap.end()){
	  totalEnergy.push_back(fShowerEnergyAlg.ShowerEnergy(planeHitsMap[plane], plane));
	}
	else{
	  totalEnergy.push_back(0);
	}
	if (initialTrackHits[plane].size()){
	  double fdEdx = 0;
	  double totQ = 0;
	  double avgT = 0;
	  double pitch = 0;
	  double wirepitch = fGeom->WirePitch(initialTrackHits[plane][0]->WireID().planeID());
	  double angleToVert = fGeom->WireAngleToVertical(fGeom->Plane(plane).View(),initialTrackHits[plane][0]->WireID().planeID()) - 0.5*TMath::Pi();
	  double cosgamma = std::abs(sin(angleToVert)*shwdir.Y()+
				     cos(angleToVert)*shwdir.Z());
	  if (cosgamma>0) pitch = wirepitch/cosgamma;
	  if (pitch){
	    if (pitch<minpitch){
	      minpitch = pitch;
	      bestPlane = plane;
	    }
	    int nhits = 0;
	    //std::cout<<"pitch = "<<pitch<<std::endl;
	    for (auto const& hit: initialTrackHits[plane]){
	      //std::cout<<hit->WireID()<<" "<<hit->PeakTime()<<" "<<std::abs((hit->WireID().Wire-initialTrackHits[plane][0]->WireID().Wire)*pitch)<<" "<<fdEdxTrackLength<<std::endl;
	      int w1 = hit->WireID().Wire;
	      int w0 = initialTrackHits[plane][0]->WireID().Wire;
	      if (std::abs((w1-w0)*pitch)<fdEdxTrackLength){
		totQ += hit->Integral();
		avgT+= hit->PeakTime();
		++nhits;
		//std::cout<<hit->WireID()<<" "<<hit->PeakTime()<<" "<<hit->Integral()<<" "<<totQ<<" "<<avgT<<std::endl;
	      }
	    }
	    if (totQ) {
	      double dQdx = totQ/(nhits*pitch);
	      fdEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avgT/nhits, initialTrackHits[plane][0]->WireID().Plane);
	    }
	  }
	  dEdx.push_back(fdEdx);
	}
	else{
	  dEdx.push_back(0);
	}
      }
      iok = 0;
      if (fDebug > 0) {
	std::cout << "Best plane is " << bestPlane << std::endl;
	std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
	std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
	std::cout << "The shower start is (" << shwxyz.X() << ", " << shwxyz.Y() << ", " << shwxyz.Z() << ")" << std::endl;
	shwxyz.Print();
      }

      return recob::Shower(shwdir, shwdirerr, shwxyz, shwxyzerr, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);
    }
  }
  return recob::Shower();
}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::OrderShowerHits(art::PtrVector<recob::Hit> const& shower, int plane) {

  /// Ordering the shower hits requires three stages:
  ///  -- putting all the hits in a given plane in some kind of order
  ///  -- use the properties of the hits in all three planes to check this order
  ///  -- orient the hits correctly using properties of the shower

  // ------------- Put hits in order ------------

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
    showerHitsMap[(*hit)->WireID().Plane].push_back(*hit);

  // Order the hits, get the RMS and the RMS gradient for the hits in this plane
  std::map<int,double> planeRMSGradients, planeRMS;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    std::vector<art::Ptr<recob::Hit> > orderedHits = FindOrderOfHits(showerHitsIt->second);
    planeRMS[showerHitsIt->first] = ShowerHitRMS(orderedHits);
    planeRMSGradients[showerHitsIt->first] = ShowerHitRMSGradient(orderedHits);
    showerHitsMap[showerHitsIt->first] = orderedHits;
  }

  if (fDebug > 0)
    for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt)
      std::cout << "Plane " << planeRMSIt->first << " has RMS " << planeRMSIt->second << " and RMS gradient " << planeRMSGradients.at(planeRMSIt->first) << std::endl;

  // ------------- Check between the views to ensure consistency of ordering -------------

  // Check between the views to make sure there isn't a poorly formed shower in just one view
  // First, determine the average RMS and RMS gradient across the other planes
  std::map<int,double> planeOtherRMS, planeOtherRMSGradients;
  for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt) {
    planeOtherRMS[planeRMSIt->first] = 0;
    planeOtherRMSGradients[planeRMSIt->first] = 0;
    int nOtherPlanes = 0;
    for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane) {
      if (plane != planeRMSIt->first and planeRMS.count(plane)) {
	planeOtherRMS[planeRMSIt->first] += planeRMS.at(plane);
	planeOtherRMSGradients[planeRMSIt->first] += planeRMSGradients.at(plane);
	++nOtherPlanes;
      }
    }
    planeOtherRMS[planeRMSIt->first] /= (double)nOtherPlanes;
    planeOtherRMSGradients[planeRMSIt->first] /= (double)nOtherPlanes;
  }

  // Look to see if one plane has a particularly high RMS (compared to the others) whilst having a similar gradient
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    if (planeRMS.at(showerHitsIt->first) > planeOtherRMS.at(showerHitsIt->first) * 2 and
	TMath::Abs(planeRMSGradients.at(showerHitsIt->first) / planeOtherRMSGradients.at(showerHitsIt->first)) < 0.1) {
      if (fDebug > 0)
	std::cout << "Plane " << showerHitsIt->first << " was perpendicular... recalculating" << std::endl;
      std::vector<art::Ptr<recob::Hit> > orderedHits = this->FindOrderOfHits(showerHitsIt->second, true);
      showerHitsMap[showerHitsIt->first] = orderedHits;
      planeRMSGradients[showerHitsIt->first] = this->ShowerHitRMSGradient(orderedHits);
    }
  }

  // ------------- Orient the shower correctly ---------------

  // If there is only one view then not much cross-checking we can do!
  //  -- reverse the shower if the RMS gradient was negative
  if (planeRMSGradients.size() < 2) {
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      if (planeRMSGradients.at(showerHitsIt->first) < 0)
	std::reverse(showerHitsIt->second.begin(), showerHitsIt->second.end());
    }
  }

  // Can do a bit more if there are two views:
  //  -- check all the gradients
  //  -- if any are significant negative gradients reverse the shower
  //  -- check the shower views
  //     -- if inconsistent, reverse the order of the smallest RMS-bin gradient view
  if (planeRMSGradients.size() == 2) {

    std::map<double,int> gradientMap;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      double gradient = planeRMSGradients.at(showerHitsIt->first);
      gradientMap[TMath::Abs(gradient)] = showerHitsIt->first;
      if (TMath::Abs(gradient) < 0.01)
	continue;
      if (gradient < 0)
	std::reverse(showerHitsIt->second.begin(), showerHitsIt->second.end());
    }

    CheckIsolatedHits(showerHitsMap);

    if (!CheckShowerHits(showerHitsMap)) {
      int planeToReverse = gradientMap.begin()->second;
      std::reverse(showerHitsMap.at(planeToReverse).begin(), showerHitsMap.at(planeToReverse).end());
    }

  }

  // Can do the most checks if we have three available views:
  //  -- check all the gradients
  //  -- if any are significant negative gradients reverse the shower
  //  -- check the shower views
  //     -- if inconsistent:
  //        -- if insignificant RMS-bin gradient shower views, reverse the shower order
  //        -- if none, reverse the smallest gradient
  if (planeRMSGradients.size() == 3) {

    if (fDebug > 0)
      for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = showerHitsMap.begin(); planeIt != showerHitsMap.end(); ++planeIt)
	std::cout << "Before any reversing: Plane " << planeIt->first << ": start is (" << HitCoordinates(planeIt->second.front()).X() << ", " << HitCoordinates(planeIt->second.front()).Y() << ")" << std::endl;

    std::map<double,int> gradientMap;
    std::vector<int> ignoredPlanes;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
      double gradient = planeRMSGradients.at(showerHitsIt->first);
      gradientMap[TMath::Abs(gradient)] = showerHitsIt->first;
      if (TMath::Abs(gradient) < 0.01) {
	ignoredPlanes.push_back(showerHitsIt->first);
	continue;
      }
      if (gradient < 0)
	std::reverse(showerHitsIt->second.begin(), showerHitsIt->second.end());
    }

    if (fDebug > 0)
      for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = showerHitsMap.begin(); planeIt != showerHitsMap.end(); ++planeIt)
	std::cout << "After reversing: Plane " << planeIt->first << ": start is (" << HitCoordinates(planeIt->second.front()).X() << ", " << HitCoordinates(planeIt->second.front()).Y() << ")" << std::endl;

    CheckIsolatedHits(showerHitsMap);

    if (fDebug > 0)
      for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = showerHitsMap.begin(); planeIt != showerHitsMap.end(); ++planeIt)
	std::cout << "After checking isolated hits: Plane " << planeIt->first << ": start is (" << HitCoordinates(planeIt->second.front()).X() << ", " << HitCoordinates(planeIt->second.front()).Y() << ")" << std::endl;

    if (!CheckShowerHits(showerHitsMap)) {
      int planeToReverse;
      if (ignoredPlanes.size())
	planeToReverse = ignoredPlanes.at(0);
      else
	planeToReverse = gradientMap.begin()->second;
      if (fDebug > 0)
	std::cout << "Plane to reverse is " << planeToReverse << std::endl;
      std::reverse(showerHitsMap.at(planeToReverse).begin(), showerHitsMap.at(planeToReverse).end());
    }

  }

  if (fDebug > 1) {
    std::cout << "End of OrderShowerHits: here are the order of hits:" << std::endl;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHitsIt = showerHitsMap.begin(); planeHitsIt != showerHitsMap.end(); ++planeHitsIt) {
      std::cout << "  Plane " << planeHitsIt->first << std::endl;
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = planeHitsIt->second.begin(); hitIt != planeHitsIt->second.end(); ++hitIt)
	std::cout << "    Hit (" << HitCoordinates(*hitIt).X() << " (real wire " << (*hitIt)->WireID() << "), " << HitCoordinates(*hitIt).Y() << ") -- pos (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
    }
  }

  return showerHitsMap;

}

void shower::EMShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
					    std::vector<art::Ptr<recob::Hit> >& showerHits,
					    art::Ptr<recob::Vertex> const& vertex){

  showerHits = FindOrderOfHits(shower);

  // Find TPC for the vertex
  double xyz[3];
  vertex->XYZ(xyz);
  geo::TPCID tpc = fGeom->FindTPCAtPosition(xyz);
  if (!tpc.isValid&&showerHits.size()) tpc = geo::TPCID(showerHits[0]->WireID());
  //std::cout<<tpc<<std::endl;
  // Find hits in the same TPC
  art::Ptr<recob::Hit> hit0, hit1;
  for (auto &hit: showerHits){
    if (hit->WireID().TPC==tpc.TPC){
      if (hit0.isNull()){
	hit0 = hit;
      }
      hit1 = hit;
    }
  }
  if (hit0.isNull()||hit1.isNull()) return;
  TVector2 coord0 = TVector2(hit0->WireID().Wire, hit0->PeakTime());
  TVector2 coord1 = TVector2(hit1->WireID().Wire, hit1->PeakTime());
  TVector2 coordvtx = TVector2(fGeom->WireCoordinate(xyz[1], xyz[2], hit0->WireID().planeID()),
			       fDetProp->ConvertXToTicks(xyz[0],  hit0->WireID().planeID()));
//  std::cout<<coord0.X()<<" "<<coord0.Y()<<std::endl;
//  std::cout<<coord1.X()<<" "<<coord1.Y()<<std::endl;
//  std::cout<<coordvtx.X()<<" "<<coordvtx.Y()<<std::endl;
//  std::cout<<hit0->WireID()<<" "<<hit1->WireID()<<std::endl;
  if ((coord1-coordvtx).Mod()<(coord0-coordvtx).Mod()){
    std::reverse(showerHits.begin(), showerHits.end());
  }
  //std::cout<<showerHits[0]->WireID()<<" "<<showerHits.back()->WireID()<<std::endl;
}

void shower::EMShowerAlg::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >const& showerHits,
					       art::Ptr<recob::Vertex> const& vertex,
					       std::vector<art::Ptr<recob::Hit> >& trackHits){

  // Find TPC for the vertex
  //std::cout<<"here"<<std::endl;
  double xyz[3];
  vertex->XYZ(xyz);
  //std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
  geo::TPCID tpc = fGeom->FindTPCAtPosition(xyz);
  //std::cout<<tpc<<std::endl;
  //vertex cannot be projected into a TPC, find the TPC that has the most hits
  if (!tpc.isValid){
    std::map<geo::TPCID, unsigned int> tpcmap;
    unsigned maxhits = 0;
    for (auto const&hit : showerHits){
      ++tpcmap[geo::TPCID(hit->WireID())];
    }
    for (auto const&t : tpcmap){
      if (t.second > maxhits){
	maxhits = t.second;
	tpc = t.first;
      }
    }
  }
  //std::cout<<tpc<<std::endl;
    //if (!tpc.isValid&&showerHits.size()) tpc = geo::TPCID(showerHits[0]->WireID());
  if (!tpc.isValid) return;
  //std::cout<<"here 1"<<std::endl;

  double parm[2];
  int fitok = 0;
  std::vector<double> wfit;
  std::vector<double> tfit;
  std::vector<double> cfit;
    
  for (size_t i = 0; i<fNfitpass; ++i){

    // Fit a straight line through hits
    unsigned int nhits = 0;
    for (auto &hit: showerHits){
      //std::cout<<i<<" "<<hit->WireID()<<" "<<tpc<<std::endl;
      if (hit->WireID().TPC==tpc.TPC){
	TVector2 coord = HitCoordinates(hit);
	//std::cout<<i<<" "<<hit->WireID()<<" "<<hit->PeakTime()<<std::endl;
	if (i==0||(std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<fToler[i-1])||fitok==1){
	  ++nhits;
	  if (nhits==fNfithits[i]+1) break;
	  wfit.push_back(coord.X());
	  tfit.push_back(coord.Y());
	//cfit.push_back(hit->Integral());
	  cfit.push_back(1.);
	  if (i==fNfitpass-1) {
	    trackHits.push_back(hit);
	  }
	//std::cout<<*hit<<std::endl;
//
//<<hit->PeakTime()<<" "<<std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<<std::endl;
	}
      }
    }
  
    if (i<fNfitpass-1&&wfit.size()){
      fitok = WeightedFit(wfit.size(), &wfit[0], &tfit[0], &cfit[0], &parm[0]);
    }
    wfit.clear();
    tfit.clear();
    cfit.clear();
  }

}


TVector2 shower::EMShowerAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) {

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());

}

TVector2 shower::EMShowerAlg::HitPosition(art::Ptr<recob::Hit> const& hit) {

  geo::PlaneID planeID = hit->WireID().planeID();

  return HitPosition(HitCoordinates(hit), planeID);

}

TVector2 shower::EMShowerAlg::HitPosition(TVector2 const& pos, geo::PlaneID planeID) {

  return TVector2(pos.X() * fGeom->WirePitch(planeID),
		  fDetProp->ConvertTicksToX(pos.Y(), planeID));

}

double shower::EMShowerAlg::GlobalWire(const geo::WireID& wireID) {

  double wireCentre[3];
  fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);

  double globalWire = -999;
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }
  else {
    // FOR COLLECTION WIRES, HARD CODE THE GEOMETRY FOR GIVEN DETECTORS
    // THIS _SHOULD_ BE TEMPORARY. GLOBAL WIRE SUPPORT IS BEING ADDED TO THE LARSOFT GEOMETRY AND SHOULD BE AVAILABLE SOON
    if (fDetector == "dune35t") {
      std::cout << "Detector is 35t" <<std::endl;
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
      else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
      else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
      else mf::LogError("BlurredClusterAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC << " (geometry" << fDetector << ")";
    }
    else if (fDetector == "dune10kt") {
      std::cout << "Detector is 10kt!"<< std::endl;
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      // Detector geometry has four TPCs, two on top of each other, repeated along z...
      int block = wireID.TPC / 4;
      globalWire = (nwires*block) + wireID.Wire;
    }
    else {
      std::cout << "Detector is other" << std::endl;
      if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
      else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
    }
  }

  return globalWire;

}

TVector2 shower::EMShowerAlg::ShowerDirection(const std::vector<art::Ptr<recob::Hit> >& showerHits) {

  TVector2 pos;
  int nhits = 0;
  double sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = showerHits.begin(); hit != showerHits.end(); ++hit) {
    ++nhits;
    pos = HitPosition(*hit);
    sumx += pos.X();
    sumy += pos.Y();
    sumx2 += pos.X() * pos.X();
    sumxy += pos.X() * pos.Y();
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  TVector2 direction = TVector2(1,gradient).Unit();

  return direction;

}

TVector2 shower::EMShowerAlg::ShowerCentre(const std::vector<art::Ptr<recob::Hit> >& showerHits) {

  TVector2 pos, chargePoint = TVector2(0,0);
  double totalCharge = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = showerHits.begin(); hit != showerHits.end(); ++hit) {
    pos = HitPosition(*hit);
    chargePoint += (*hit)->Integral() * pos;
    totalCharge += (*hit)->Integral();
  }
  TVector2 centre = chargePoint / totalCharge;

  return centre;

}

double shower::EMShowerAlg::ShowerHitRMS(const std::vector<art::Ptr<recob::Hit> >& showerHits) {

  TVector2 direction = ShowerDirection(showerHits);
  TVector2 centre = ShowerCentre(showerHits);

  std::vector<double> distanceToAxis;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitsIt = showerHits.begin(); showerHitsIt != showerHits.end(); ++showerHitsIt) {
    TVector2 proj = (HitPosition(*showerHitsIt) - centre).Proj(direction) + centre;
    distanceToAxis.push_back((HitPosition(*showerHitsIt) - proj).Mod());
  }
  double RMS = TMath::RMS(distanceToAxis.begin(), distanceToAxis.end());

  return RMS;

}

double shower::EMShowerAlg::ShowerHitRMSGradient(const std::vector<art::Ptr<recob::Hit> >& showerHits) {

  // Don't forget to clean up the header file!
  bool makeDirectionPlot = false;

  // Find a rough shower 'direction' and centre
  TVector2 direction = ShowerDirection(showerHits);

  // Bin the hits into discreet chunks
  int nShowerSegments = 5;
  double lengthOfShower = (HitPosition(showerHits.back()) - HitPosition(showerHits.front())).Mod();
  double lengthOfSegment = lengthOfShower / (double)nShowerSegments;
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerSegments;
  std::map<int,double> segmentCharge;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitsIt = showerHits.begin(); showerHitsIt != showerHits.end(); ++showerHitsIt) {
    showerSegments[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment].push_back(*showerHitsIt);
    segmentCharge[(int)(HitPosition(*showerHitsIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment] += (*showerHitsIt)->Integral();
  }

  TGraph* graph = new TGraph();
  std::vector<std::pair<int,double> > binVsRMS;

  // Loop over the bins to find the distribution of hits as the shower progresses
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerSegmentIt = showerSegments.begin(); showerSegmentIt != showerSegments.end(); ++showerSegmentIt) {

    // Get the mean position of the hits in this bin
    TVector2 meanPosition(0,0);
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt)
      meanPosition += HitPosition(*hitInSegmentIt);
    meanPosition /= (double)showerSegmentIt->second.size();

    // Get the RMS of this bin
    std::vector<double> distanceToAxisBin;
    for (std::vector<art::Ptr<recob::Hit> >::iterator hitInSegmentIt = showerSegmentIt->second.begin(); hitInSegmentIt != showerSegmentIt->second.end(); ++hitInSegmentIt) {
      TVector2 proj = (HitPosition(*hitInSegmentIt) - meanPosition).Proj(direction) + meanPosition;
      distanceToAxisBin.push_back((HitPosition(*hitInSegmentIt) - proj).Mod());
    }

    double RMSBin = TMath::RMS(distanceToAxisBin.begin(), distanceToAxisBin.end());
    if (makeDirectionPlot)
      graph->SetPoint(graph->GetN(), showerSegmentIt->first, RMSBin);//*segmentCharge.at(showerSegmentIt->first));
    binVsRMS.push_back(std::make_pair(showerSegmentIt->first, RMSBin));

  }

  // Get the gradient of the RMS-bin plot
  int nhits = 0;
  double sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (std::vector<std::pair<int,double> >::iterator binVsRMSIt = binVsRMS.begin(); binVsRMSIt != binVsRMS.end(); ++binVsRMSIt) {
    ++nhits;
    sumx += binVsRMSIt->first;
    sumy += binVsRMSIt->second;
    sumx2 += binVsRMSIt->first * binVsRMSIt->first;
    sumxy += binVsRMSIt->first * binVsRMSIt->second;
  }
  double RMSgradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);

  if (makeDirectionPlot) {
    TCanvas* canv = new TCanvas();
    graph->Fit("pol1");
    TF1* fit = graph->GetFunction("pol1");
    Double_t graphGradient = fit->GetParameter(1);
    graph->Draw();
    canv->SaveAs("direction.png");
    if (fDebug > 0)
      std::cout << "Gradient from graph is " << graphGradient << " and from vector is " << RMSgradient << std::endl;
  }
  delete graph;

  return RMSgradient;

}

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, geo::PlaneID planeID) {

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  geo::TPCID tpcID = fGeom->FindTPCAtPosition(pointPosition);
  int tpc = 0;
  if (tpcID.isValid)
    tpc = tpcID.TPC;
  else
    tpc = 0;

  TVector2 wireTickPos = TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), planeID.Plane, tpc % 2, 0),
				  fDetProp->ConvertXToTicks(point.X(), planeID.Plane, tpc % 2, 0));

  //return wireTickPos;
  return HitPosition(wireTickPos, planeID);

}

Int_t shower::EMShowerAlg::WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm){

  Double_t sumx=0.;
  Double_t sumx2=0.;
  Double_t sumy=0.;
  Double_t sumy2=0.;
  Double_t sumxy=0.;
  Double_t sumw=0.;
  Double_t eparm[2];
    
  parm[0]  = 0.;
  parm[1]  = 0.;
  eparm[0] = 0.;
  eparm[1] = 0.;
    
  for (Int_t i=0; i<n; i++) {
    sumx += x[i]*w[i];
    sumx2 += x[i]*x[i]*w[i];
    sumy += y[i]*w[i]; 
    sumy2 += y[i]*y[i]*w[i];
    sumxy += x[i]*y[i]*w[i];
    sumw += w[i];
  }
    
  if (sumx2*sumw-sumx*sumx==0.) return 1;
  if (sumx2-sumx*sumx/sumw==0.) return 1;
    
  parm[0] = (sumy*sumx2-sumx*sumxy)/(sumx2*sumw-sumx*sumx);
  parm[1] = (sumxy-sumx*sumy/sumw)/(sumx2-sumx*sumx/sumw);
    
  eparm[0] = sumx2*(sumx2*sumw-sumx*sumx);
  eparm[1] = (sumx2-sumx*sumx/sumw);
    
  if (eparm[0]<0. || eparm[1]<0.) return 1;
    
  eparm[0] = sqrt(eparm[0])/(sumx2*sumw-sumx*sumx);
  eparm[1] = sqrt(eparm[1])/(sumx2-sumx*sumx/sumw);
    
  return 0;
    
}

bool shower::EMShowerAlg::isCleanShower(std::vector<art::Ptr<recob::Hit> > const& hits){

  if (!hits.size()) return false;
  if (hits.size()>2000) return true;
  if (hits.size()<20) return true;
  std::map<int, int> hitmap;
  unsigned nhits = 0;
  for (auto const&hit : hits){
    ++nhits;
    if (nhits>2)
      ++hitmap[hit->WireID().Wire];
    if (nhits==20) break;
  }
  //std::cout<<hits.size()<<" "<<float(nhits-2)/hitmap.size()<<std::endl;
  if (float(nhits-2)/hitmap.size()>1.4) return false;
  else return true;
}

shower::HitPosition::HitPosition()
  : fGeom(lar::providerFrom<geo::Geometry>())
  , fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {}


// // Code to make wire maps showing the global wire coordinates
// struct TPCWire {
//   int fPlane, fTPC, fWire, fGlobalWire;
//   TVector3 fStart, fEnd;
//   TPCWire(int plane, int tpc, int wire, int globalWire, TVector3 start, TVector3 end) {
//     fPlane = plane;
//     fTPC = tpc;
//     fWire = wire;
//     fGlobalWire = globalWire;
//     fStart = start;
//     fEnd = end;
//   }
//   TPCWire() { }
//   void SetProps(int plane, int tpc, int wire, int globalWire, TVector3 start, TVector3 end) {
//     fPlane = plane;
//     fTPC = tpc;
//     fWire = wire;
//     fGlobalWire = globalWire;
//     fStart = start;
//     fEnd = end;
//   }
// };

// void shower::EMShowerAlg::MakePicture() {

//   std::vector<TPCWire> allWires;

//   for (geo::WireID const& wireID : fGeom->IterateWireIDs()) {

//     if (wireID.TPC % 2 == 0)
//       continue;

//     double xyzStart[3], xyzEnd[3];
//     fGeom->WireEndPoints(wireID, xyzStart, xyzEnd);
//     int globalWire = GlobalWire(wireID);

//     allWires.emplace_back(wireID.Plane, wireID.TPC, wireID.Wire, globalWire, TVector3(xyzStart[0],xyzStart[1],xyzStart[2]), TVector3(xyzEnd[0],xyzEnd[1],xyzEnd[2]));

//   } // for all wires

//   TCanvas* uplane = new TCanvas();
//   TCanvas* vplane = new TCanvas();
//   TCanvas* zplane = new TCanvas();
//   TText* number = new TText();
//   number->SetTextSize(0.002);
//   TLine* line = new TLine(0,100,0,100);
//   line->SetLineWidth(0.5);
//   uplane->Range(-5,-100,160,140);
//   vplane->Range(-5,-100,160,140);
//   zplane->Range(-5,-100,160,140);

//   for (std::vector<TPCWire>::iterator wireIt = allWires.begin(); wireIt != allWires.end(); ++wireIt) {
//     if (wireIt->fPlane == 0)
//       uplane->cd();
//     else if (wireIt->fPlane == 1)
//       vplane->cd();
//     else if (wireIt->fPlane == 2)
//       zplane->cd();
//     line->DrawLine(wireIt->fStart.Z(), wireIt->fStart.Y(), wireIt->fEnd.Z(), wireIt->fEnd.Y());
//     number->DrawText(wireIt->fStart.Z(), wireIt->fStart.Y(), (std::to_string(wireIt->fGlobalWire)+" ("+std::to_string(wireIt->fWire)+")").c_str());
//     number->DrawText(wireIt->fEnd.Z(), wireIt->fEnd.Y(), (std::to_string(wireIt->fGlobalWire)+" ("+std::to_string(wireIt->fWire)+")").c_str());
//   }

//   uplane->SaveAs("UPlane.pdf");
//   vplane->SaveAs("VPlane.pdf");
//   zplane->SaveAs("ZPlane.pdf");

// }

//
// // Timing code...
//
// auto start_time = std::chrono::high_resolution_clock::now();
// // Put stuff here!
// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count();
// std::cout << "Duration is " << duration/1000000.0 << " s " << std::endl;;
//
