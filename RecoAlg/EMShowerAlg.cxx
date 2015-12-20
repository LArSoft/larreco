////////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#include "RecoAlg/EMShowerAlg.h"

shower::EMShowerAlg::EMShowerAlg(fhicl::ParameterSet const& pset) : fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
                                                                    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
								    fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")) {
  fMinTrackLength  = pset.get<double>("MinTrackLength");
  fdEdxTrackLength = pset.get<double>("dEdxTrackLength");
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
      if (clusterHitTracks.at(0)->Length() < fMinTrackLength) continue;

      // Add this cluster to the track map
      int track = clusterHitTracks.at(0).key();
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

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::CheckShowerHits(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap,
										       std::map<int,double> const& goodnessOfOrderMap) {

  /// Takes the shower hits in all views and ensures the ordering is consistent between views, fixing things if not
  /// This method makes sure the two views are consistent and sorts things out if not.

  std::map<int,std::vector<art::Ptr<recob::Hit> > > newShowerHitsMap;
  std::vector<int> planesToReverse;

  if (orderedShowerMap.size() == 1)
    newShowerHitsMap = orderedShowerMap;

  else if (orderedShowerMap.size() == 2) {

    art::Ptr<recob::Hit> start1 = *orderedShowerMap.begin()->second.begin();
    art::Ptr<recob::Hit> start2 = *orderedShowerMap.rbegin()->second.begin();

    // Find the start point of the two views
    TVector3 start = Construct3DPoint(start1, start2);
    std::cout << "The start is 3D is (" << start.X() << ", " << start.Y() << ", " << start.Z() << ")" << std::endl;

    // Project this back onto each plane
    TVector2 proj1 = Project3DPointOntoPlane(start, start1->WireID().planeID());
    TVector2 proj2 = Project3DPointOntoPlane(start, start2->WireID().planeID());
    //std::cout << "Projection 1 lands (" << proj1.X() << ", " << proj1.Y() << ") and projection2 lands (" << proj2.X() << ", " << proj2.Y() << ")" << std::endl;

    // Get the average distance from this projection to the actual start
    std::cout << "Distance from start " << (HitPosition(start1) - proj1).Mod() << " and from end " << (HitPosition(start2) - proj2).Mod() << std::endl;
    double projDiff = ( (HitPosition(start1) - proj1).Mod() + (HitPosition(start2) - proj2).Mod() ) / (double)2;

    std::cout << "Proj diff is " << projDiff << std::endl;

    // Do we need to swap the order in one of the views?
    bool flip = false;
    if (TMath::Abs(start1->PeakTime() - start2->PeakTime()) > 40)
      flip = true;
    else if (start.X() == -9999 or start.Y() == -9999 or start.Z() == -9999)
      flip = true;
    else if (projDiff > 10)
      flip = true;

    std::cout << "Here's the goodness of order... " << std::endl;
    for (std::map<int,double>::const_iterator orderIt = goodnessOfOrderMap.begin(); orderIt != goodnessOfOrderMap.end(); ++orderIt)
      std::cout << "Plane " << orderIt->first << " has goodness of order " << orderIt->second << std::endl;

    // Which view?
    if (flip) {
      if (goodnessOfOrderMap.begin()->second > goodnessOfOrderMap.rbegin()->second)
	planesToReverse.push_back(goodnessOfOrderMap.rbegin()->first);
      else
	planesToReverse.push_back(goodnessOfOrderMap.begin()->first);
    }

    if (planesToReverse.size() == 1)
      std::cout << "Flipping plane " << planesToReverse.at(0) << std::endl;

  }

  else if (orderedShowerMap.size() == 3) {

    planesToReverse = IdentifyBadPlanes(orderedShowerMap);

  }

  // Make the new shower hits map
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerPlaneIt = orderedShowerMap.begin(); showerPlaneIt != orderedShowerMap.end(); ++showerPlaneIt) {
    if (std::find(planesToReverse.begin(), planesToReverse.end(), showerPlaneIt->first) == planesToReverse.end()) {
      newShowerHitsMap[showerPlaneIt->first] = showerPlaneIt->second;
    }
    else {
      std::vector<art::Ptr<recob::Hit> > tmpShowerHits = showerPlaneIt->second;
      std::reverse(tmpShowerHits.begin(), tmpShowerHits.end());
      newShowerHitsMap[showerPlaneIt->first] = tmpShowerHits;
    }
  }

  return newShowerHitsMap;

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

    // Look at how many clusters each plane has, and the proportion of hits each of uses
    std::map<int,std::vector<double> > planeClusterSizes;
    for (std::map<int,std::vector<art::Ptr<recob::Cluster> > >::iterator planeClustersIt = planeClusters.begin(); planeClustersIt != planeClusters.end(); ++planeClustersIt) {
      for (std::vector<art::Ptr<recob::Cluster> >::iterator planeClusterIt = planeClustersIt->second.begin(); planeClusterIt != planeClustersIt->second.end(); ++planeClusterIt) {
	std::vector<art::Ptr<recob::Hit> > hits = fmh.at(planeClusterIt->key());
        planeClusterSizes[planeClustersIt->first].push_back((double)hits.size()/(double)planeHits.at(planeClustersIt->first).size());
      }
    }

    // Find the average hit fraction across all clusters in the plane
    std::map<int,double> planeClusterAverageSizes;
    for (std::map<int,std::vector<double> >::iterator planeClusterSizesIt = planeClusterSizes.begin(); planeClusterSizesIt != planeClusterSizes.end(); ++planeClusterSizesIt) {
      double average = 0;
      for (std::vector<double>::iterator planeClusterSizeIt = planeClusterSizesIt->second.begin(); planeClusterSizeIt != planeClusterSizesIt->second.end(); ++planeClusterSizeIt)
	average += *planeClusterSizeIt;
      average /= planeClusterSizesIt->second.size();
      planeClusterAverageSizes[planeClusterSizesIt->first] = average;
    }

    // Now, decide if there is one plane which is ruining the reconstruction
    // If there are two planes with a low average cluster fraction and one with a high one, this plane likely merges two particle deposits together
    std::vector<int> highAverage, lowAverage;
    for (std::map<int,double>::iterator planeClusterAverageSizeIt = planeClusterAverageSizes.begin(); planeClusterAverageSizeIt != planeClusterAverageSizes.end(); ++planeClusterAverageSizeIt) {
      if (planeClusterAverageSizeIt->second > 0.9) highAverage.push_back(planeClusterAverageSizeIt->first);
      else lowAverage.push_back(planeClusterAverageSizeIt->first);
    }
    int badPlane = -1;
    if (highAverage.size() == 1 and highAverage.size() < lowAverage.size())
      badPlane = highAverage.at(0);

    if (badPlane != -1) 
      for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = planeClusters.at(badPlane).begin(); clusterIt != planeClusters.at(badPlane).end(); ++clusterIt)
	clustersToIgnore.push_back(clusterIt->key());

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

  std::cout << "About to make me a track from these 'ere 'its... " << std::endl;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit1 = track1.begin(); hit1 != track1.end(); ++hit1)
    std::cout << "Hit (" << HitCoordinates(*hit1).X() << ", " << HitCoordinates(*hit1).Y() << ") (real wire " << (*hit1)->WireID().Wire << ") in TPC " << (*hit1)->WireID().TPC << std::endl;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit2 = track2.begin(); hit2 != track2.end(); ++hit2)
    std::cout << "Hit (" << HitCoordinates(*hit2).X() << ", " << HitCoordinates(*hit2).Y() << ") (real wire " << (*hit2)->WireID().Wire << ") in TPC " << (*hit2)->WireID().TPC << std::endl;

  pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(track1, track2);

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
  try { pitch = track->PitchInView(trackHits.at(0)->View()); }
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

// void shower::EMShowerAlg::FindInitialTrack(art::PtrVector<recob::Hit> const& hits,
// 					   std::unique_ptr<recob::Track>& initialTrack,
// 					   std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits) {

//   /// Finds the initial track-like part of the shower and the hits in all views associated with it

//   art::FindManyP<recob::Cluster> fmc;

//   return FindInitialTrack(hits, initialTrack, initialTrackHits, fmc, 1);

// }

void shower::EMShowerAlg::FindInitialTrack(art::PtrVector<recob::Hit> const& hits,
					   std::unique_ptr<recob::Track>& initialTrack,
					   std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits,
					   art::FindManyP<recob::Cluster> const& fmc, int plane) {

  /// Finds the initial track-like part of the shower and the hits in all views associated with it

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->WireID().Plane].push_back(*hit);

  std::map<int,std::vector<art::Ptr<recob::Hit> > > orderedShowerMap;
  std::map<int,double> goodnessOfOrderMap;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::cout << std::endl << "Plane " << planeHits->first << std::endl;
    if (planeHits->first != plane and plane != -1) continue;

    std::vector<art::Ptr<recob::Hit> > showerHits;
    double goodness = OrderShowerHits(planeHits->second, showerHits, fmc);
    orderedShowerMap[planeHits->first] = showerHits;
    goodnessOfOrderMap[planeHits->first] = goodness;

    std::cout << "After ordering: Plane " << planeHits->first << ": start is determined to be (" << HitCoordinates(*showerHits.begin()).X() << ", " << HitCoordinates(*showerHits.begin()).Y() << ") (TPC " << (*showerHits.begin())->WireID().TPC << ", wire " << (*showerHits.begin())->WireID().Wire << ") and the end is determined to be (" << HitCoordinates(*showerHits.rbegin()).X() << ", " << HitCoordinates(*showerHits.rbegin()).Y() << ") (TPC " << (*showerHits.rbegin())->WireID().TPC << ", wire " << (*showerHits.rbegin())->WireID().Wire << ")" << std::endl;

  }

  std::cout << std::endl;

  //return;

  // Now we want to ensure that all views are consistent with each other
  std::map<int,std::vector<art::Ptr<recob::Hit> > > showerHitsMap = CheckShowerHits(orderedShowerMap, goodnessOfOrderMap);

  // Now find the hits belonging to the track
  initialTrackHits = FindShowerStart(showerHitsMap);

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator trackHitIt = initialTrackHits.begin(); trackHitIt != initialTrackHits.end(); ++trackHitIt)
    std::cout << "Plane " << trackHitIt->first << " has track start (" << HitCoordinates(*trackHitIt->second.begin()).X() << ", " << HitCoordinates(*trackHitIt->second.begin()).Y() << ") and end (" << HitCoordinates(*trackHitIt->second.rbegin()).X() << ", " << HitCoordinates(*trackHitIt->second.rbegin()).Y() << ")" << std::endl;

  // Now we have the track hits -- can make a track!
  initialTrack = MakeInitialTrack(initialTrackHits, showerHitsMap);

  if (initialTrack)
    std::cout << "Found track! Start is (" << initialTrack->Vertex().X() << ", " << initialTrack->Vertex().Y() << ", " << initialTrack->Vertex().Z() << ") and end is (" << initialTrack->End().X() << ", " << initialTrack->End().Y() << ", " << initialTrack->End().Z() << "); direction is (" << initialTrack->VertexDirection().X() << ", " << initialTrack->VertexDirection().Y() << ", " << initialTrack->VertexDirection().Z() << ")" << std::endl;

  return;

}

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindOrderOfHits(std::vector<art::Ptr<recob::Hit> > const& hits) {

  /// Orders hits along the best fit line through the charge weighted centre of the hits

  // Find the charge-weighted centre (in [cm]) of this shower
  TVector2 pos, chargePoint = TVector2(0,0);
  double totalCharge = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    pos = HitPosition(*hit);
    chargePoint += (*hit)->Integral() * pos;
    totalCharge += (*hit)->Integral();
  }
  TVector2 centre = chargePoint / totalCharge;

  // Find a rough shower 'direction'
  int nhits = 0;
  double sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    ++nhits;
    pos = HitPosition(*hit);
    sumx += pos.X();
    sumy += pos.Y();
    sumx2 += pos.X() * pos.X();
    sumxy += pos.X() * pos.Y();
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  TVector2 direction = TVector2(1,gradient).Unit();

  // Find how far each hit (projected onto this axis) is from the centre
  std::map<double,art::Ptr<recob::Hit> > hitProjection;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    pos = HitPosition(*hit) - centre;
    hitProjection[direction*pos] = *hit;
  }

  // Get a vector of hits in order of the shower
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(hitProjection.begin(), hitProjection.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

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
      std::cout << "WARNING! Number of showers this track matches is " << matchingShowers.size() << std::endl;

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

    // Quality check -- make sure there isn't a single hit in a different TPC (artefact of tracking failure)
    std::vector<art::Ptr<recob::Hit> > newInitialHits;
    if (initialHits.size() > 1) {
      std::map<int,int> tpcHitMap;
      std::vector<int> singleHitTPCs;
      for (std::vector<art::Ptr<recob::Hit> >::iterator initialHitIt = initialHits.begin(); initialHitIt != initialHits.end(); ++initialHitIt)
	++tpcHitMap[(*initialHitIt)->WireID().TPC];
      for (std::map<int,int>::iterator tpcIt = tpcHitMap.begin(); tpcIt != tpcHitMap.end(); ++tpcIt)
	if (tpcIt->second == 1) singleHitTPCs.push_back(tpcIt->first);
      if (singleHitTPCs.size()) {
	for (std::vector<int>::iterator tpcIt = singleHitTPCs.begin(); tpcIt != singleHitTPCs.end(); ++tpcIt)
	  std::cout << *tpcIt << std::endl;
	for (std::vector<art::Ptr<recob::Hit> >::iterator initialHitIt = initialHits.begin(); initialHitIt != initialHits.end(); ++initialHitIt)
	  if (std::find(singleHitTPCs.begin(), singleHitTPCs.end(), (*initialHitIt)->WireID().TPC) == singleHitTPCs.end())
	    newInitialHits.push_back(*initialHitIt);
      }
      else
	newInitialHits = initialHits;
    }
    else newInitialHits = initialHits;

    initialHitsMap[orderedShowerIt->first] = newInitialHits;

  }

  return initialHitsMap;

}

std::vector<double> shower::EMShowerAlg::GetShowerDirectionProperties(std::vector<art::Ptr<recob::Hit> > const& showerHits, TVector2 const& direction, std::string end) {

  /// Returns some ints which can be used to determine the likelihood that the shower propagates along the hits as ordered in showerHits

  TCanvas* canv = new TCanvas();
  TGraph* graph = new TGraph();

  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(2);
  canv->cd();
  graph->Draw("AP");

  // // Don't forget to clean up the .h file!

  // Count how often the perpendicular distance from the shower axis increases
  int spreadingHits = 0;
  double previousDistanceFromAxisPos = 0, previousDistanceFromAxisNeg = 0;
  int hitsAfterOffshoot = 0;
  int hitsAbove = 0, hitsBelow = 0;
  TVector2 pos, projPos, start = HitPosition(*showerHits.begin());
  //std::cout << "Start is (" << start.X() << ", " << start.Y() << ")" << std::endl;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = showerHits.begin(); hitIt != showerHits.end(); ++hitIt) {
    pos = HitPosition(*hitIt) - start;
    graph->SetPoint(graph->GetN(),pos.X(),pos.Y());
    double angle = pos.DeltaPhi(direction);
    projPos = (pos).Proj(direction);
    double distanceFromAxis = TMath::Sqrt(TMath::Power(pos.Mod(),2) - TMath::Power(projPos.Mod(),2));
    //std::cout << "X is " << pos.X() << " and distance from axis is " << distanceFromAxis << std::endl;
    //std::cout << "Hit (" << HitCoordinates(*hitIt).X() << ", " << HitCoordinates(*hitIt).Y() << ") has (relative position is (" << pos.X() << ", " << pos.Y() << ")); angle " << angle << " and projected distance " << distanceFromAxis << std::endl;
    if (distanceFromAxis > 1.2) {
      if (angle > 0) ++hitsAbove;
      if (angle < 0) ++hitsBelow;
      if ( (angle > 0 and distanceFromAxis > previousDistanceFromAxisPos) or
	   (angle < 0 and distanceFromAxis > previousDistanceFromAxisNeg) ) {
	++spreadingHits;
	hitsAfterOffshoot = 0;
	if (distanceFromAxis > 0) previousDistanceFromAxisPos = distanceFromAxis;
	else previousDistanceFromAxisNeg = distanceFromAxis;
      }
      else {
	++hitsAfterOffshoot;
	--spreadingHits;
      }
      if (hitsAfterOffshoot > 5) {
	hitsAfterOffshoot = 0;
	previousDistanceFromAxisPos = 0;
	previousDistanceFromAxisNeg = 0;
      }
    }
    //std::cout << "End of hit; spreadingHits is " << spreadingHits << " and offshoot is " << hitsAfterOffshoot << std::endl;
  }

  //std::cout << "Total hits away from centre is " << hitsAbove+hitsBelow << std::endl;
  std::vector<double> directionProperties = {(double)spreadingHits, hitsAbove/(double)showerHits.size(), hitsBelow/(double)showerHits.size()};

  std::cout << "Direction properties above " << directionProperties.at(1) << ", below " << directionProperties.at(2) << std::endl;

  TLine line;
  line.SetLineColor(2);
  canv->cd();
  line.DrawLine(-1000*direction.X(),-1000*direction.Y(),1000*direction.X(),1000*direction.Y());

  canv->SaveAs(TString("ShowerPlane")+TString(end)+TString(".png"));

  return directionProperties;

}

std::vector<int> shower::EMShowerAlg::IdentifyBadPlanes(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) {

  /// Takes the shower hits in all views and decides if a particular plane(s) is inconsistent with the other two.
  /// Normally this involves the hits being ordered incorrectly.

  std::vector<int> badPlanes;

  // Put the start hits of the showers in a map
  std::map<int,art::Ptr<recob::Hit> > start2DMap;
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerHitIt = showerHitsMap.begin(); showerHitIt != showerHitsMap.end(); ++showerHitIt)
    start2DMap[showerHitIt->first] = *showerHitIt->second.begin();

  std::map<int,double> projDiff;
  std::map<int,bool> isolated;

  for (int plane = 0; plane < 3; ++plane) {

    // Looking for two things -- a difference in the average projection of a hit made from two planes in the third,
    // and an isolated projected hit
    projDiff[plane] = 0;
    isolated[plane] = false;

    // Get the other two planes
    std::vector<int> otherPlanes;
    for (int otherPlane = 0; otherPlane < 3; ++otherPlane)
      if (plane != otherPlane)
	otherPlanes.push_back(otherPlane);

    // Find the start point for the other two views
    TVector3 start3D = Construct3DPoint(start2DMap.at(otherPlanes.at(0)), start2DMap.at(otherPlanes.at(1)));

    // Project this back onto each plane
    for (int allPlanes = 0; allPlanes < 3; ++allPlanes) {
      TVector2 proj = Project3DPointOntoPlane(start3D, start2DMap.at(allPlanes)->WireID().planeID());
      projDiff[plane] += (proj - HitPosition(start2DMap.at(allPlanes))).Mod();
      double minDist = 1e9;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = showerHitsMap.at(allPlanes).begin(); hitIt != showerHitsMap.at(allPlanes).end(); ++hitIt) {
	double separation = (HitPosition(*hitIt) - proj).Mod();
	if (separation < minDist) minDist = separation;
      }
      std::cout << "Plane " << plane << " has a minDist " << minDist << std::endl;
      if (minDist > 15) isolated[plane] = true;
    }
    projDiff[plane] /= 3;

  }

  for (std::map<int,double>::iterator projDiffIt = projDiff.begin(); projDiffIt != projDiff.end(); ++projDiffIt)
    std::cout << "For combo " << projDiffIt->first << ", the average difference after projecting is " << projDiffIt->second << std::endl;

  std::vector<int> highAverage, lowAverage, isolations;
  for (std::map<int,double>::iterator projDiffIt = projDiff.begin(); projDiffIt != projDiff.end(); ++projDiffIt) {
    if (projDiffIt->second > 15) highAverage.push_back(projDiffIt->first);
    else lowAverage.push_back(projDiffIt->first);
  }
  for (std::map<int,bool>::iterator isolationIt = isolated.begin(); isolationIt != isolated.end(); ++isolationIt)
    if (isolationIt->second) isolations.push_back(isolationIt->first);

  // If the projection onto a particular plane was far from all the shower hits something is wrong!
  if (isolations.size() == 1) {
    std::vector<int> otherPlanes;
    for (int plane = 0; plane < 3; ++plane)
      if (std::find(isolations.begin(), isolations.end(), plane) == isolations.end()) otherPlanes.push_back(plane);
    if (projDiff.at(otherPlanes.at(0)) < projDiff.at(otherPlanes.at(1)))
      badPlanes.push_back(otherPlanes.at(0));
    else
      badPlanes.push_back(otherPlanes.at(1));
  }

  // If no isolations but one plane has a much lower average projection than the other two, this is also an indication something is wrong
  else if (highAverage.size() == 2 and lowAverage.size() == 1)
    badPlanes.push_back(*lowAverage.begin());

  return badPlanes;

}

std::unique_ptr<recob::Track> shower::EMShowerAlg::MakeInitialTrack(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap,
								    std::map<int,std::vector<art::Ptr<recob::Hit> > > const& showerHitsMap) {

  /// Takes initial track hits from multiple views and forms a track object which best represents the start of the shower
  /// Returns false if this isn't possible.

  // Don't want to make any assumptions on the number of planes, or the number of planes with hits on
  // Try to resolve any conflicts as we go...

  std::unique_ptr<recob::Track> track;

  // Look at planes
  std::vector<int> planes, emptyPlanes;
  std::map<geo::PlaneID,TVector2> showerCentreMap;
  for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {
    if (showerHitsMap.count(plane) == 0) {
      emptyPlanes.push_back(plane);
      continue;
    }
    planes.push_back(plane);
    // Find the charge-weighted centre (in [wire/tick]) of this shower
    TVector2 pos, chargePoint = TVector2(0,0);
    double totalCharge = 0;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = showerHitsMap.at(plane).begin(); hit != showerHitsMap.at(plane).end(); ++hit) {
      pos = HitPosition(*hit);
      chargePoint += (*hit)->Integral() * pos;
      totalCharge += (*hit)->Integral();
    }
    TVector2 centre = chargePoint / totalCharge;
    showerCentreMap[showerHitsMap.at(plane).at(0)->WireID().planeID()] = centre;
  } 

  // If there's only one plane, there's not much we can do!
  if (planes.size() == 1)
    return track;

  // If there are two planes then we can make the object
  else if (planes.size() == 2) {
    if (TMath::Abs(HitCoordinates(initialHitsMap.at(planes.at(0)).at(0)).Y() - HitCoordinates(initialHitsMap.at(planes.at(1)).at(0)).Y()) <= 15)
      track = ConstructTrack(initialHitsMap.at(planes.at(0)), initialHitsMap.at(planes.at(1)), showerCentreMap);
  }

  // If we have three planes to consider then we can do more checks!
  else if (planes.size() == 3) {
    std::vector<std::vector<art::Ptr<recob::Hit> > > showerPlanes;
    // See if there's a plane we should avoid
    std::vector<int> planesToAvoid = IdentifyBadPlanes(showerHitsMap);
    if (planesToAvoid.size() == 1) {
      for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator initialHitsIt = initialHitsMap.begin(); initialHitsIt != initialHitsMap.end(); ++initialHitsIt)
	if (std::find(planesToAvoid.begin(), planesToAvoid.end(), initialHitsIt->first) == planesToAvoid.end())
	  showerPlanes.push_back(initialHitsIt->second);
    }
    else {
      // Take the longest two initial tracks
      std::map<int,int> trackSizeMap;
      for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerPlaneIt = initialHitsMap.begin(); showerPlaneIt != initialHitsMap.end(); ++showerPlaneIt)
	trackSizeMap[showerPlaneIt->second.size()] = showerPlaneIt->first;
      if (trackSizeMap.size() == 1)
	return track;
      for (std::map<int,int>::reverse_iterator trackSizeIt = trackSizeMap.rbegin(); trackSizeIt != trackSizeMap.rend(); ++trackSizeIt) {
	if (std::distance(trackSizeMap.rbegin(), trackSizeIt) > 1)
	  break;
	showerPlanes.push_back(initialHitsMap.at(trackSizeIt->second));
	std::cout << "Using plane " << trackSizeIt->second << std::endl;
      }
    }
    if (TMath::Abs(HitCoordinates(showerPlanes.at(0).at(0)).Y() - HitCoordinates(showerPlanes.at(1).at(0)).Y()) <= 15)
      track = ConstructTrack(showerPlanes.at(0), showerPlanes.at(1), showerCentreMap);
  }

  return track;

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
  std::vector<double> totalEnergy(3), totalEnergyError(3), dEdx(3), dEdxError(3);

  // Look at each of the planes
  for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {

    // If hits exist on this plane...
    if (planeHitsMap.count(plane) != 0) {
      dEdx.push_back(FinddEdx(initialHitsMap.at(plane), initialTrack));
      totalEnergy.push_back(fShowerEnergyAlg.ShowerEnergy(planeHitsMap.at(plane), plane));
      if (planeHitsMap.at(plane).size() > highestNumberOfHits) {
	bestPlane = plane;
	highestNumberOfHits = planeHitsMap.at(plane).size();
      }
    }

    // If there are no hits on this plane
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

double shower::EMShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower,
					    std::vector<art::Ptr<recob::Hit> >& showerHits,
					    art::FindManyP<recob::Cluster> const& fmc) {

  /// Takes the hits associated with a shower and orders then so they follow the direction of the shower

  showerHits = FindOrderOfHits(shower);

  // Get the clusters associated with this shower in this plane
  std::vector<art::Ptr<recob::Cluster> > lineClusters;
  std::vector<int> clustersInVector;
  std::map<int,std::vector<art::Ptr<recob::Hit> > > lineClusterHits;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    std::vector<art::Ptr<recob::Cluster> > associatedClusters = fmc.at(hit->key());
    for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = associatedClusters.begin(); clusterIt != associatedClusters.end(); ++clusterIt) {
      lineClusterHits[clusterIt->key()].push_back(*hit);
      if (std::find(clustersInVector.begin(), clustersInVector.end(), clusterIt->key()) == clustersInVector.end()) {
	lineClusters.push_back(*clusterIt);
	clustersInVector.push_back(clusterIt->key());
      }
    }
  }

  std::cout << "Number of line clusters is " << lineClusters.size() << std::endl;

  // Find the longest line cluster and get its direction
  std::map<int,art::Ptr<recob::Cluster> > clusterSizes;
  for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = lineClusters.begin(); clusterIt != lineClusters.end(); ++clusterIt)
    clusterSizes[(*clusterIt)->NHits()] = *clusterIt;
  art::Ptr<recob::Cluster> longestLineCluster = clusterSizes.rbegin()->second;
  std::vector<art::Ptr<recob::Hit> > longestLineClusterHits = lineClusterHits.at(longestLineCluster.key());//FindOrderOfHits(lineClusterHits.at(longestLineCluster.key()));

  // Try instead using the direction of shower through the charge weighted centre (classic Wallbank-recon)!
  // Find the charge-weighted centre (in [cm]) of this shower
  TVector2 pos, chargePoint = TVector2(0,0);
  double totalCharge = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    pos = HitPosition(*hit);
    chargePoint += (*hit)->Integral() * pos;
    totalCharge += (*hit)->Integral();
  }
  TVector2 showerCentre = chargePoint / totalCharge;

  // Find a rough shower 'direction'
  int nhits = 0;
  double sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    ++nhits;
    pos = HitPosition(*hit);
    sumx += pos.X();
    sumy += pos.Y();
    sumx2 += pos.X() * pos.X();
    sumxy += pos.X() * pos.Y();
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  TVector2 showerDirection = TVector2(1,gradient).Unit();

  // Now find the start and end for each of the other showers (assume they start by the central axis and move outwards)
  std::map<int,std::pair<TVector2,TVector2> > clusterEnds;
  for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = lineClusters.begin(); clusterIt != lineClusters.end(); ++clusterIt) {

    if (clusterIt->key() == longestLineCluster.key())
      continue;

    // Order the cluster hits along their central axis
    std::vector<art::Ptr<recob::Hit> > clusterHits = FindOrderOfHits(lineClusterHits.at(clusterIt->key()));

    // Get the two ends
    TVector2 clusterStart = HitPosition(*clusterHits.begin()), clusterEnd = HitPosition(*clusterHits.rbegin());
    //std::cout << "Cluster " << clusterIt->key() << " start is (" << HitCoordinates(*clusterHits.begin()).X() << ", " << HitCoordinates(*clusterHits.begin()).Y() << ") and end is (" << HitCoordinates(*clusterHits.rbegin()).X() << ", " << HitCoordinates(*clusterHits.rbegin()).Y() << ") (size is " << lineClusterHits.at(clusterIt->key()).size() << ")" << std::endl;

    // Now we have the two ends, can project both onto the longest line cluster and see which is closest
    double distance1 = (((clusterStart-showerCentre).Proj(showerDirection))+showerCentre - clusterStart).Mod();
    double distance2 = (((clusterEnd-showerCentre).Proj(showerDirection))+showerCentre - clusterEnd).Mod();

    //std::cout << "Distance from start is " << distance1 << " and from end is " << distance2 << std::endl;

    // And make a pair of the ends
    if (distance1 < distance2)
      clusterEnds[clusterIt->key()] = std::make_pair(clusterStart, clusterEnd);
    else {
      //std::cout << "Flipped" << std::endl;
      clusterEnds[clusterIt->key()] = std::make_pair(clusterEnd, clusterStart);
    }

  }

  // Find the cluster(s) associated with the two ends of the shower
  std::vector<art::Ptr<recob::Hit> >::iterator beginIt = showerHits.begin(), beginSecondIt = ++beginIt;
  std::vector<art::Ptr<recob::Hit> >::reverse_iterator endIt = showerHits.rbegin(), endSecondIt = ++endIt;
  std::vector<art::Ptr<recob::Cluster> > startClusters = fmc.at(beginIt->key()), endClusters = fmc.at(endIt->key());
  if (startClusters.size() == 0) startClusters = fmc.at(beginSecondIt->key());
  if (endClusters.size() == 0) endClusters = fmc.at(endSecondIt->key());
  std::vector<int> startClusterNums, endClusterNums;
  std::transform(startClusters.begin(), startClusters.end(), std::back_inserter(startClusterNums), [](art::Ptr<recob::Cluster> const& cluster){return cluster.key();});
  std::transform(endClusters.begin(), endClusters.end(), std::back_inserter(endClusterNums), [](art::Ptr<recob::Cluster> const& cluster){return cluster.key();});

  // If there aren't many clusters, see if the longest cluster is associated with either of the ends
  if (lineClusters.size() <= 3) {
    bool startIsStart = false, endIsStart = false;
    for (std::vector<art::Ptr<recob::Cluster> >::iterator startClusterIt = startClusters.begin(); startClusterIt != startClusters.end(); ++startClusterIt)
      if (startClusterIt->key() == longestLineCluster.key()) startIsStart = true;
    for (std::vector<art::Ptr<recob::Cluster> >::iterator endClusterIt = endClusters.begin(); endClusterIt != endClusters.end(); ++endClusterIt)
      if (endClusterIt->key() == longestLineCluster.key()) endIsStart = true;
    if (startIsStart and !endIsStart) {
      std::cout << "Returning start is start" << std::endl;
      return 0.55;
    }
    else if (endIsStart and !startIsStart) {
      std::reverse(showerHits.begin(), showerHits.end());
      std::cout << "Returning end is start!" << std::endl;
      return 0.45;
    }
  }

  // Get the average angle between the line clusters and the shower direction
  double averageAngleStart = 0, averageAngleEnd = 0;
  int clustersCountedFromStart = 0, clustersCountedFromEnd = 0;
  TVector2 showerStart = HitPosition(*showerHits.begin()), showerEnd = HitPosition(*showerHits.rbegin());
  for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = lineClusters.begin(); clusterIt != lineClusters.end(); ++clusterIt) {

    if (clusterIt->key() == longestLineCluster.key())
      continue;
    // if (lineClusterHits.at(clusterIt->key()).size() < 5)
    //   continue;

    // The 'ends' of the cluster
    std::pair<TVector2,TVector2> thisClusterEnds = clusterEnds.at(clusterIt->key());
    TVector2 clusterDirection = (thisClusterEnds.second - thisClusterEnds.first).Unit();

    // From start
    if (std::find(startClusterNums.begin(), startClusterNums.end(), clusterIt->key()) == startClusterNums.end()) {
      ++clustersCountedFromStart;

      // Angle the shower from the start
      TVector2 direction;
      if ( (thisClusterEnds.first-showerStart)*showerDirection > 1 )
	direction = showerDirection;
      else
	direction = -1 * showerDirection;

      averageAngleStart += TMath::Abs(direction.DeltaPhi(clusterDirection));

      std::cout << "Cluster " << clusterIt->key() << " has angle from start " << TMath::Abs(direction.DeltaPhi(clusterDirection));

    }

    // From end
    if (std::find(endClusterNums.begin(), endClusterNums.end(), clusterIt->key()) == endClusterNums.end()) {
      ++clustersCountedFromEnd;

      // Angle the shower from the end
      TVector2 direction;
      if ( (thisClusterEnds.first-showerEnd)*showerDirection > 1 )
	direction = showerDirection;
      else
	direction = -1 * showerDirection;

      averageAngleEnd += TMath::Abs(direction.DeltaPhi(clusterDirection));

      std::cout << ", and angle from end " << TMath::Abs(direction.DeltaPhi(clusterDirection)) << std::endl;

    }

  }

  averageAngleStart /= clustersCountedFromStart;
  averageAngleEnd /= clustersCountedFromEnd;

  std::cout << "Average angle from start is " << averageAngleStart << " and from end is " << averageAngleEnd << std::endl;

  double goodness;
  if (averageAngleEnd < averageAngleStart) {
    std::reverse(showerHits.begin(), showerHits.end());
    goodness = 1 - (averageAngleEnd / (double)(averageAngleEnd+averageAngleStart));
    std::cout << "Reversing end and start" << std::endl;
  }
  else
    goodness = 1 - (averageAngleStart / (double)(averageAngleStart+averageAngleEnd));

  std::cout << "About to return goodness of " << goodness << std::endl;

  return goodness;

}

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower) {

  /// Takes the hits associated with a shower and orders them so they follow the direction of the shower

  std::vector<art::Ptr<recob::Hit> > showerHits = FindOrderOfHits(shower);

  TVector2 direction;

  // Find the shower direction properties
  std::vector<double> propertiesFromBeginning = GetShowerDirectionProperties(showerHits, direction, "Beginning");
  std::reverse(showerHits.begin(), showerHits.end());
  std::vector<double> propertiesFromEnd = GetShowerDirectionProperties(showerHits, direction, "End");
  
  // Decide which is the end
  //double goodnessOfOrder = 0;
  double asymmetryDetermination = 0.5;
  bool fromBeginning = true;
  if ( (propertiesFromBeginning.at(1) > asymmetryDetermination or propertiesFromBeginning.at(2) > asymmetryDetermination) and
       (propertiesFromEnd.at(1) > asymmetryDetermination or propertiesFromEnd.at(2) > asymmetryDetermination) ) {
    int hitAsymmetryBeginning = TMath::Max(propertiesFromBeginning.at(0), propertiesFromBeginning.at(1));
    int hitAsymmetryEnd = TMath::Max(propertiesFromEnd.at(0), propertiesFromEnd.at(1));
    std::cout << "Both beginning and end have missquewed directions.  Asymmetry from beginning is " << hitAsymmetryBeginning << " and end is " << hitAsymmetryEnd << std::endl;
    if (hitAsymmetryEnd < hitAsymmetryBeginning) fromBeginning = false;
  }
  else if ( (propertiesFromBeginning.at(1) > asymmetryDetermination or propertiesFromBeginning.at(2) > asymmetryDetermination) and
	    propertiesFromEnd.at(1) <= asymmetryDetermination and propertiesFromEnd.at(2) <= asymmetryDetermination ) {
    std::cout << "Only looking from the beginning has a lot of asymmetry.  Above is " << propertiesFromBeginning.at(1) << " and below is " << propertiesFromBeginning.at(2) << ".  Pick the end" << std::endl;
    fromBeginning = false;
  }
  else if ( (propertiesFromEnd.at(1) > asymmetryDetermination or propertiesFromEnd.at(2) > asymmetryDetermination) and
	    propertiesFromBeginning.at(1) <= asymmetryDetermination and propertiesFromBeginning.at(2) <= asymmetryDetermination ) {
    std::cout << "Only looking from the end has a lot of asymmetry.  Above is " << propertiesFromEnd.at(1) << " and below is " << propertiesFromEnd.at(2) << ".  Pick the end" << std::endl;
    fromBeginning = true;
  }
  else if ( propertiesFromBeginning.at(1) <= asymmetryDetermination and propertiesFromBeginning.at(2) <= asymmetryDetermination and
	    propertiesFromEnd.at(1) <= asymmetryDetermination and propertiesFromEnd.at(2) <= asymmetryDetermination ) {
    std::cout << "From beginning is " << propertiesFromBeginning.at(0) << " and from end is " << propertiesFromEnd.at(0) << std::endl;
    if (propertiesFromBeginning.at(0) < propertiesFromEnd.at(0)) fromBeginning = false;
  }
  
  if (fromBeginning)
    std::reverse(showerHits.begin(), showerHits.end());

  return showerHits;

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

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, geo::PlaneID planeID) {

  /// Projects a 3D point (units [cm]) onto a 2D plane
  /// Returns 2D point (units [cm])

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  TVector2 wireTickPos = TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), planeID.Plane, fGeom->FindTPCAtPosition(pointPosition).TPC % 2, 0),
				  fDetProp->ConvertXToTicks(point.X(), planeID.Plane, fGeom->FindTPCAtPosition(pointPosition).TPC % 2, 0));

  //return wireTickPos;
  return HitPosition(wireTickPos, planeID);

}
