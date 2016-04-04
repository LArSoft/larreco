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

shower::EMShowerAlg::EMShowerAlg(fhicl::ParameterSet const& pset)
  : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  , fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg"))
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  , fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg"))
{
  fMinTrackLength  = pset.get<double>("MinTrackLength");
  fdEdxTrackLength = pset.get<double>("dEdxTrackLength");
  fNfitpass        = pset.get<unsigned int>("Nfitpass");
  fNfithits        = pset.get<std::vector<unsigned int> >("Nfithits");
  fToler           = pset.get<std::vector<double> >("Toler");
  if (fNfitpass!=fNfithits.size()||
      fNfitpass!=fToler.size()){
    throw art::Exception(art::errors::Configuration)
      <<"EMShowerAlg: fNfithits and fToler need to have size fNfitpass";
  }

  hTrueDirection = tfs->make<TH1I>("trueDir","",2,0,2);
}

void shower::EMShowerAlg::AssociateClustersAndTracks(std::vector<art::Ptr<recob::Cluster> > const& clusters,
						     art::FindManyP<recob::Hit> const& fmh,
						     art::FindManyP<recob::Track> const& fmt,
						     std::map<int,std::vector<int> >& clusterToTracks,
						     std::map<int,std::vector<int> >& trackToClusters) {

  /// Map associated tracks and clusters together given their associated hits

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

  /// Map associated tracks and clusters together given their associated hits, whilst ignoring certain clusters

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
	//std::cout << "Track " << clusterHitTracks.at(0)->ID() << " is too short! (" << clusterHitTracks.at(0)->Length() << ")" << std::endl;
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
      if (std::find(clusterToTracks[cluster].begin(), clusterToTracks[cluster].end(), track) == clusterToTracks[cluster].end()) {
	// if (trackID == 65550 or trackID == 65549 or trackID == 65548 or trackID == 65536)
	//   std::cout << "Track " << trackID << " is associated with cluster " << cluster << std::endl;
	clusterToTracks[cluster].push_back(track);
      }

    }

  }

  return;

}

bool shower::EMShowerAlg::CheckShowerHits(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) {

  /// Takes the shower hits in all views and ensure the ordering is consistent
  /// Returns bool, indicating whether or not everything makes sense!

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

      // std::cout << "Plane... " << plane << std::endl;
      // std::cout << "Start position in this plane is " << HitPosition(start2DMap.at(plane)).X() << ", " << HitPosition(start2DMap.at(plane)).Y() << ")" << std::endl;
      // std::cout << "Shower start from other two planes is (" << showerStartPos.X() << ", " << showerStartPos.Y() << ", " << showerStartPos.Z() << ")" << std::endl;
      // std::cout << "Projecting the other two planes gives position (" << showerStartProj.X() << ", " << showerStartProj.Y() << ")" << std::endl;

      double projDiff = TMath::Abs((showerStartProj-HitPosition(start2DMap.at(plane))).Mod());
      double timeDiff = TMath::Max(TMath::Abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(0))->PeakTime()),
				   TMath::Abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(1))->PeakTime()));

      std::cout << "Plane " << plane << " has projDiff " << projDiff << " and timeDiff " << timeDiff << std::endl;
      if (projDiff > 1 or timeDiff > 40)
	consistencyCheck = false;

    }

  }

  return consistencyCheck;

}

void shower::EMShowerAlg::CheckShowerPlanes(std::vector<std::vector<int> > const& initialShowers,
					    std::vector<int>& clustersToIgnore,
					    std::vector<art::Ptr<recob::Cluster> > const& clusters,
					    art::FindManyP<recob::Hit> const& fmh) {

  /// Takes the initial showers found and tries to resolve issues where one bad view ruins the event

  for (std::vector<std::vector<int> >::const_iterator initialShowerIt = initialShowers.begin(); initialShowerIt != initialShowers.end(); ++initialShowerIt) {

    // Make maps of all clusters and cluster hits in each view
    std::map<int,std::vector<art::Ptr<recob::Cluster> > > planeClusters;
    std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHits;

    // Loop over the clusters comprising this shower
    for (std::vector<int>::const_iterator clusterIt = initialShowerIt->begin(); clusterIt != initialShowerIt->end(); ++clusterIt) {
      art::Ptr<recob::Cluster> cluster = clusters.at(*clusterIt);
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(cluster.key());
      planeClusters[cluster->Plane().Plane].push_back(cluster);
      for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
	planeHits[(*hitIt)->WireID().Plane].push_back(*hitIt);
    }

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

    std::cout << "Looking at bad plane: Shower " << std::distance(initialShowers.begin(),initialShowerIt) << std::endl;
    for (std::map<int,double>::iterator planeClustersAvSizesIt = planeClustersAvSizes.begin(); planeClustersAvSizesIt != planeClustersAvSizes.end(); ++planeClustersAvSizesIt)
      std::cout << "Plane " << planeClustersAvSizesIt->first << " has average hit thing " << planeClustersAvSizesIt->second << std::endl;

    // Now decide if there is one plane which is ruining the reconstruction
    // If two planes have a low average cluster fraction and one high, this plane likely merges two particle deposits together
    int badPlane = -1;
    for (std::map<int,double>::iterator clusterAvSizeIt = planeClustersAvSizes.begin(); clusterAvSizeIt != planeClustersAvSizes.end(); ++clusterAvSizeIt) {
      double avClusterSizeThisPlane = clusterAvSizeIt->second;
      double sumClusterSizeOtherPlanes = 0;
      for (std::map<int,double>::iterator otherClustersAvSizeIt = planeClustersAvSizes.begin(); otherClustersAvSizeIt != planeClustersAvSizes.end(); ++otherClustersAvSizeIt)
	if (clusterAvSizeIt->first != otherClustersAvSizeIt->first)
    	  sumClusterSizeOtherPlanes += otherClustersAvSizeIt->second;
      if (avClusterSizeThisPlane >= sumClusterSizeOtherPlanes)
    	badPlane = clusterAvSizeIt->first;
    }

    // // Now, decide if there is one plane which is ruining the reconstruction
    // // If there are two planes with a low average cluster fraction and one with a high one, this plane likely merges two particle deposits together
    // std::vector<int> highAverage, lowAverage;
    // for (std::map<int,double>::iterator planeClusterAverageSizeIt = planeClustersAvSizes.begin(); planeClusterAverageSizeIt != planeClustersAvSizes.end(); ++planeClusterAverageSizeIt) {
    //   if (planeClusterAverageSizeIt->second > 0.9) highAverage.push_back(planeClusterAverageSizeIt->first);
    //   else lowAverage.push_back(planeClusterAverageSizeIt->first);
    // }
    // int badPlane = -1;
    // if (highAverage.size() == 1 and highAverage.size() < lowAverage.size())
    //   badPlane = highAverage.at(0);

    if (badPlane != -1) {
      std::cout << "Bad plane is " << badPlane << std::endl;
      for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeClusters.at(badPlane).begin(); clusterIt != planeClusters.at(badPlane).end(); ++clusterIt)
	clustersToIgnore.push_back(clusterIt->key());
    }

  }

  return;

}

TVector3 shower::EMShowerAlg::Construct3DPoint(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2) {

  /// Constructs a 3D point (in [cm]) to represent the hits given in two views

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

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All PMA methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).
  /// This implementation also orients the track in the correct direction if a map of shower centres (units [cm]) in each view is provided.

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
    unsigned int firstTPC1 = (*hits1.begin())->WireID().TPC, firstTPC2 = (*hits2.begin())->WireID().TPC;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits1.begin(); hitIt != hits1.end(); ++hitIt)
      if ((*hitIt)->WireID().TPC == firstTPC1) track1.push_back(*hitIt);
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits2.begin(); hitIt != hits2.end(); ++hitIt)
      if ((*hitIt)->WireID().TPC == firstTPC2) track2.push_back(*hitIt);    
  }
  else {
    track1 = hits1;
    track2 = hits2;
  }

  if (debug) {
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
      TVector3 dc((*pmatrack)[i+1]->Point3D());
      dc -= (*pmatrack)[i]->Point3D();
      dc *= 1.0 / dc.Mag();
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

  if (debug)
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

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).

  std::map<geo::PlaneID,TVector2> showerCentreMap;

  return this->ConstructTrack(track1, track2, showerCentreMap);

}

double shower::EMShowerAlg::FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, std::unique_ptr<recob::Track> const& track) {

  /// Finds dE/dx for the track given a set of hits

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
      // double dEdx = fCalorimetryAlg.dEdx_AREA(shower.at(*trackHitIt), pitch);
      // std::cout << "After " << totalDistance << ", the dE/dx is " << dEdx << std::endl;
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

  /// Finds the initial track-like part of the shower and the hits in all views associated with it

  // First, order the hits into the correct shower order in each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap = OrderShowerHits(hits, plane);

  // Now find the hits belonging to the track
  initialTrackHits = FindShowerStart(showerHitsMap);

  std::cout << "Here are the initial shower hits... " << std::endl;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator initialTrackHitsIt = initialTrackHits.begin(); initialTrackHitsIt != initialTrackHits.end(); ++initialTrackHitsIt) {
    std::cout << "  Plane " << initialTrackHitsIt->first << std::endl;
    for (std::vector<art::Ptr<recob::Hit> >::iterator initialTrackHitIt = initialTrackHitsIt->second.begin(); initialTrackHitIt != initialTrackHitsIt->second.end(); ++initialTrackHitIt)
      std::cout << "    Hit is (" << HitCoordinates(*initialTrackHitIt).X() << " (real hit " << (*initialTrackHitIt)->WireID() << "), " << HitCoordinates(*initialTrackHitIt).Y() << ")" << std::endl;
  }

  // Now we have the track hits -- can make a track!
  initialTrack = MakeInitialTrack(initialTrackHits);

  if (initialTrack) {
    std::cout << "The track start is " << std::cout;
    initialTrack->Vertex().Print();
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

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindOrderOfHits(std::vector<art::Ptr<recob::Hit> > const& hits) {

  /// Orders hits along the best fit line through the charge weighted centre of the hits

  // Find the charge-weighted centre (in [cm]) of this shower
  TVector2 centre = ShowerCentre(hits);

  // Find a rough shower 'direction'
  TVector2 direction = ShowerDirection(hits);

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
  //   std::cout << "Hit at wire " << (*hitIt)->WireID() << " and tick " << (*hitIt)->PeakTime() << " is pos (" << HitPosition(*hitIt).X() << ", " << HitPosition(*hitIt).Y() << ")" << std::endl;
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

void shower::EMShowerAlg::FindShowers(std::map<int,std::vector<int> > const& trackToClusters,
				      std::vector<std::vector<int> >& showers) {

  /// Makes showers given a map between tracks and all clusters associated with them

  // Loop over all tracks 
  for (std::map<int,std::vector<int> >::const_iterator trackToClusterIt = trackToClusters.begin(); trackToClusterIt != trackToClusters.end(); ++ trackToClusterIt) {

    // Find which showers already made are associated with this track
    std::vector<int> matchingShowers;
    for (unsigned int shower = 0; shower < showers.size(); ++shower)
      for (std::vector<int>::const_iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if ( (std::find(showers.at(shower).begin(), showers.at(shower).end(), *cluster) != showers.at(shower).end()) and
	     (std::find(matchingShowers.begin(), matchingShowers.end(), shower)) == matchingShowers.end() )
	  matchingShowers.push_back(shower);

    // Shouldn't be more than one
    if (matchingShowers.size() > 1)
      mf::LogInfo("EMShowerAlg") << "Number of showers this track matches is " << matchingShowers.size() << std::endl;

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

  return;

}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::FindShowerStart(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap) {

  /// Takes a map of the shower hits on each plane (ordered from what has been decided to be the start)
  /// Returns a map of the initial track-like part of the shower on each plane

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
      if (debug)
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

  /// Takes initial track hits from multiple views and forms a track object which best represents the start of the shower

  std::unique_ptr<recob::Track> track;

  std::map<int,std::vector<art::Ptr<recob::Hit> > > trackLengthMap;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator initialHitsIt = initialHitsMap.begin(); initialHitsIt != initialHitsMap.end(); ++initialHitsIt)
    trackLengthMap[initialHitsIt->second.size()] = initialHitsIt->second;

  // Can't do much with just one view
  if (trackLengthMap.size() < 2) {
    mf::LogInfo("EMShowerAlg") << "Only one useful view for this shower; nothing can be done" << std::endl;
    return track;
  }

  std::vector<std::vector<art::Ptr<recob::Hit> > > longestTracks;
  std::transform(trackLengthMap.rbegin(), trackLengthMap.rend(), std::back_inserter(longestTracks),
		 [](std::pair<int,std::vector<art::Ptr<recob::Hit> > > const& initialTrackHits) { return initialTrackHits.second; });

  return ConstructTrack(longestTracks.at(0), longestTracks.at(1));

}

recob::Shower shower::EMShowerAlg::MakeShower(art::PtrVector<recob::Hit> const& hits,
					      std::unique_ptr<recob::Track> const& initialTrack,
					      std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap) {

  /// Makes a recob::Shower object given the hits in the shower and the initial track-like part

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

  std::cout << "Best plane is " << bestPlane << std::endl;
  std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
  std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
  std::cout << "The shower start is " << std::endl;
  showerStart.Print();

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
      std::cout << "Best plane is " << bestPlane << std::endl;
      std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
      std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
      std::cout << "The shower start is " << std::endl;
      shwxyz.Print();

      return recob::Shower(shwdir, shwdirerr, shwxyz, shwxyzerr, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);
    }
  }
  return recob::Shower();
}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::OrderShowerHits(art::PtrVector<recob::Hit> const& shower, int plane) {

  /// Takes the hits associated with a shower and orders them so they follow the direction of the shower

  // Don't forget to clean up the header file!

  // Save RMS, and the gradient of the RMS for each plane
  std::map<int,double> planeRMSGradients;
  std::map<int,double> planeRMS;

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
    showerHitsMap[(*hit)->WireID().Plane].push_back(*hit);

  // Find if the shower isn't well-formed along a central axis in one view
  int badPlane = -1;
  std::map<int,double> planeOtherRMS;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt)
    planeRMS[showerHitsIt->first] = ShowerHitRMS(showerHitsIt->second);
  if (planeRMS.size() == 2) {
    std::vector<int> planes;
    for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt)
      planes.push_back(planeRMSIt->first);
    for (std::vector<int>::iterator plane1It = planes.begin(); plane1It != planes.end(); ++plane1It)
      for (std::vector<int>::iterator plane2It = planes.begin(); plane2It != planes.end(); ++plane2It)
  	if (*plane1It != *plane2It)
  	  planeOtherRMS[*plane1It] = planeRMS.at(*plane2It);
  }
  else if (planeRMS.size() > 2) {
    for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt) {
      std::vector<double> otherRMSs;
      for (int plane = 0; plane < (int)fGeom->MaxPlanes(); ++plane)
  	if (plane != planeRMSIt->first)
  	  otherRMSs.push_back(planeRMS.at(plane));
      double avOtherRMSs = 0;
      for (std::vector<double>::iterator otherRMSIt = otherRMSs.begin(); otherRMSIt != otherRMSs.end(); ++otherRMSIt)
  	avOtherRMSs += *otherRMSIt;
      avOtherRMSs /= (double)otherRMSs.size();
      planeOtherRMS[planeRMSIt->first] = avOtherRMSs;
    }
  }
  for (std::map<int,double>::iterator planeOtherRMSIt = planeOtherRMS.begin(); planeOtherRMSIt != planeOtherRMS.end(); ++planeOtherRMSIt)
    if (planeRMS.at(planeOtherRMSIt->first) > planeOtherRMSIt->second * 2.5) {
      badPlane = planeOtherRMSIt->first;
      std::cout << "Too high: " << badPlane << std::endl;
      mf::LogInfo("EMShowerAlg") << "Ommitting view " << badPlane << " for this shower; hits do not appear to lie along an axis" << std::endl;
    }

  std::cout << "Bad plane is " << badPlane << std::endl;

  // Order the hits if they appear to be well defined along a shower 'axis'
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {
    if (showerHitsIt->first != plane and plane != -1) continue;
    if (showerHitsIt->first == badPlane)
      continue;
    showerHitsMap[showerHitsIt->first] = this->FindOrderOfHits(showerHitsIt->second);
  }

  // Order the hits in each plane first
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitsIt = showerHitsMap.begin(); showerHitsIt != showerHitsMap.end(); ++showerHitsIt) {

    // std::cout << "Hits in order for plane " << showerHitsIt->first << ":" << std::endl;
    // for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitIt = showerHitsIt->second.begin(); showerHitIt != showerHitsIt->second.end(); ++showerHitIt)
    //   std::cout << "Hit at position (" << HitPosition(*showerHitIt).X() << ", " << HitPosition(*showerHitIt).Y() << ") has real wire " << (*showerHitIt)->WireID() << std::endl;

    std::cout << "Plane " << showerHitsIt->first << " has start (" << HitCoordinates(showerHitsIt->second.front()).X() << " (real wire " << showerHitsIt->second.front()->WireID() << "), " << HitCoordinates(showerHitsIt->second.front()).Y() << ") and end (" << HitCoordinates(showerHitsIt->second.back()).X() << " (real wire " << showerHitsIt->second.back()->WireID() << "), " << HitCoordinates(showerHitsIt->second.back()).Y() << ")" << std::endl;

    // First, order the hits along the shower
    // Then we need to see if this is correct or if we need to swap the order
    std::vector<art::Ptr<recob::Hit> > showerHits = showerHitsIt->second;

    // Find a rough shower 'direction' and centre
    TVector2 direction = ShowerDirection(showerHits);

    // Bin the hits into discreet chunks
    int nShowerSegments = 5;
    double lengthOfShower = (HitPosition(showerHits.back()) - HitPosition(showerHits.front())).Mod();
    double lengthOfSegment = lengthOfShower / (double)nShowerSegments;
    std::map<int,std::vector<art::Ptr<recob::Hit> > > showerSegments;
    std::map<int,double> segmentCharge;
    for (std::vector<art::Ptr<recob::Hit> >::iterator showerHitIt = showerHits.begin(); showerHitIt != showerHits.end(); ++showerHitIt) {
      showerSegments[(int)(HitPosition(*showerHitIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment].push_back(*showerHitIt);
      segmentCharge[(int)(HitPosition(*showerHitIt)-HitPosition(showerHits.front())).Mod() / lengthOfSegment] += (*showerHitIt)->Integral();
    }

    //TGraph* graph = new TGraph();
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
      //graph->SetPoint(graph->GetN(), showerSegmentIt->first, RMSBin);//*segmentCharge.at(showerSegmentIt->first));
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

    // TCanvas* canv = new TCanvas();
    // graph->Fit("pol1");
    // TF1* fit = graph->GetFunction("pol1");
    // Double_t graphGradient = fit->GetParameter(1);
    // graph->Draw();
    // canv->SaveAs("direction.png");
    // std::cout << "Gradient from graph is " << graphGradient << " and from vector is " << RMSgradient << std::endl;
    //delete graph;

    planeRMSGradients[showerHitsIt->first] = RMSgradient;

  }

  for (std::map<int,double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end(); ++planeRMSIt)
    std::cout << "Plane " << planeRMSIt->first << " has RMS " << planeRMSIt->second << " and RMS gradient " << planeRMSGradients.at(planeRMSIt->first) << std::endl;

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

    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = showerHitsMap.begin(); planeIt != showerHitsMap.end(); ++planeIt)
    std::cout << "After reversing: Plane " << planeIt->first << ": start is (" << HitCoordinates(planeIt->second.front()).X() << ", " << HitCoordinates(planeIt->second.front()).Y() << ")" << std::endl;

    if (!CheckShowerHits(showerHitsMap)) {
      int planeToReverse;
      if (ignoredPlanes.size())
	planeToReverse = ignoredPlanes.at(0);
      else
	planeToReverse = gradientMap.begin()->second;
      std::cout << "Plane to reverse is " << planeToReverse << std::endl;
      std::reverse(showerHitsMap.at(planeToReverse).begin(), showerHitsMap.at(planeToReverse).end());
    }

  }

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = showerHitsMap.begin(); planeIt != showerHitsMap.end(); ++planeIt)
    std::cout << "End of OrderShowerHits: Plane " << planeIt->first << ": start is (" << HitCoordinates(planeIt->second.front()).X() << ", " << HitCoordinates(planeIt->second.front()).Y() << ")" << std::endl;

  return showerHitsMap;

}

void shower::EMShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
					    std::vector<art::Ptr<recob::Hit> >& showerHits,
					    art::Ptr<recob::Vertex> const& vertex){

  /// Takes the hits associated with a shower and orders then so they follow the direction of the shower

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

  /// Return the coordinates of this hit in global wire/tick space

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());

}

TVector2 shower::EMShowerAlg::HitPosition(art::Ptr<recob::Hit> const& hit) {

  /// Return the coordinates of this hit in units of cm

  geo::PlaneID planeID = hit->WireID().planeID();

  return HitPosition(HitCoordinates(hit), planeID);

}

TVector2 shower::EMShowerAlg::HitPosition(TVector2 const& pos, geo::PlaneID planeID) {

  /// Return the coordinates of this hit in units of cm

  return TVector2(pos.X() * fGeom->WirePitch(planeID),
		  fDetProp->ConvertTicksToX(pos.Y(), planeID));

}

double shower::EMShowerAlg::GlobalWire(geo::WireID wireID) {

  /// Find the global wire position

  double wireCentre[3];
  fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);

  double globalWire = -999;
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }
  else {
    unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
    if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
    else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
    else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
    else mf::LogError("EMShowerAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC;
  }

  return globalWire;

}

TVector2 shower::EMShowerAlg::ShowerDirection(const std::vector<art::Ptr<recob::Hit> >& showerHits) {

  /// Returns a rough shower 'direction' given the hits in the shower

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

  /// Returns the charge-weighted shower centre

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

  /// Returns the RMS of the hits from the central shower 'axis' along the length of the shower

  TVector2 direction = ShowerDirection(showerHits);
  TVector2 centre = ShowerCentre(showerHits);

  std::vector<double> distanceToAxis;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitIt = showerHits.begin(); showerHitIt != showerHits.end(); ++showerHitIt) {
    TVector2 proj = (HitPosition(*showerHitIt) - centre).Proj(direction) + centre;
    distanceToAxis.push_back((HitPosition(*showerHitIt) - proj).Mod());
  }
  double RMS = TMath::RMS(distanceToAxis.begin(), distanceToAxis.end());

  return RMS;

}

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, geo::PlaneID planeID) {

  /// Projects a 3D point (units [cm]) onto a 2D plane
  /// Returns 2D point (units [cm])

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
