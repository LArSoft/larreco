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
                                                                    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")) {
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

  /// Map associated tracks and clusters together given their associated hits

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

void shower::EMShowerAlg::MakeShowers(std::map<int,std::vector<int> > const& trackToClusters,
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

void shower::EMShowerAlg::FindShowerProperties(art::PtrVector<recob::Hit> const& hits, art::FindManyP<recob::Track> const& fmt,
					       TVector3& direction, TVector3& directionError, TVector3& vertex, TVector3& vertexError,
					       std::vector<double>& totalEnergy, std::vector<double>& totalEnergyError, std::vector<double>& dEdx, std::vector<double>& dEdxError,
					       int& bestPlane) {

  /// Finds the properties of the shower from the hits in it

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  // Save information about the hits in each of the views
  std::map<int,art::Ptr<recob::Hit> > vertexMap;
  std::map<int,TVector2> showerCentreMap;
  std::map<int,art::Ptr<recob::Track> > trackMap;
  std::map<int,std::vector<art::Ptr<recob::Hit> > > trackHitsMap;
  std::map<int,std::pair<std::vector<art::Ptr<recob::Hit> >, std::vector<art::Ptr<recob::Hit> > > > trackHitsBothEndsMap;

  unsigned int highestNumberOfHits = 0;

  // Consider each plane separately to determine shower properties
  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::cout << "Plane " << planeHits->first << std::endl;

    // Find best plane
    if (planeHits->second.size() > highestNumberOfHits) {
      highestNumberOfHits = planeHits->second.size();
      bestPlane = planeHits->first;
    }

    // Find the charge-weighted centre of this shower
    TVector2 pos, chargePoint = TVector2(0,0);
    double totalCharge = 0;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator planeHit = planeHits->second.begin(); planeHit != planeHits->second.end(); ++planeHit) {
      pos = HitCoordinates(*planeHit);
      chargePoint += (*planeHit)->Integral() * pos;
      totalCharge += (*planeHit)->Integral();
    }
    TVector2 centre = chargePoint / totalCharge;
    showerCentreMap[planeHits->first] = centre;

    // Find rough 'ends' of the shower!
    art::Ptr<recob::Hit> end1, end2;
    FindShowerEnds(planeHits->second, centre, end1, end2);

    // // Decide which end is the vertex
    // art::Ptr<recob::Hit> vertex = FindVertex(planeHits->second, end1, end2);
    // art::Ptr<recob::Hit> showerEnd = vertex.key() == end1.key() ? end2 : end1;
    // vertexMap.push_back(vertex);

    // // Get the initial track
    // std::vector<int> trackHits = FindTrack(planeHits->second, HitCoordinates(vertex), HitCoordinates(showerEnd));
    // art::Ptr<recob::Track> track = fmt.at(hits.at(trackHits.at(0)).key()).at(0);

    // Get the initial track
    std::vector<art::Ptr<recob::Hit> > trackHits1 = FindTrack(planeHits->second, HitCoordinates(end1), HitCoordinates(end2));
    std::vector<art::Ptr<recob::Hit> > trackHits2 = FindTrack(planeHits->second, HitCoordinates(end2), HitCoordinates(end1));
    trackHitsBothEndsMap[planeHits->first] = std::make_pair(trackHits1, trackHits2);
    std::vector<art::Ptr<recob::Hit> > trackHits = trackHits1.size() > trackHits2.size() ? trackHits1 : trackHits2;
    trackHitsMap[planeHits->first] = trackHits;

    // Find a track object
    art::Ptr<recob::Track> track;
    for (std::vector<art::Ptr<recob::Hit> >::iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
      std::vector<art::Ptr<recob::Track> > tracks = fmt.at(trackHitIt->key());
      if (tracks.size() == 0)
	continue;
      track = tracks.at(0);
      break;
    }
    trackMap[planeHits->first] = track;

    // Vertex
    art::Ptr<recob::Hit> planevertex = *trackHits.begin();
    vertexMap[planeHits->first] = planevertex;

  }

  // Find 3D vertex and direction
  art::Ptr<recob::Track> vertexTrack = FindVertexTrack(vertexMap, trackMap, trackHitsMap);
  if (vertexTrack.isNull()) return;
  FindShowerStartDirection(vertexTrack, showerCentreMap, vertex, direction);

  // Use the 3D vertex to resolve any problems in 2D
  ProjectVertexIn2D(vertex, trackHitsMap, trackHitsBothEndsMap);

  // Find energy and dE/dx
  for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {
    if (planeHitsMap.count(plane) != 0) {
      dEdx.push_back(FinddEdx(trackHitsMap.at(plane), vertexTrack, plane, vertexMap.at(plane)->View()));
      totalEnergy.push_back(fShowerEnergyAlg.ShowerEnergy(planeHitsMap.at(plane), plane));
    }
    else {
      dEdx.push_back(0);
      totalEnergy.push_back(0);
    }
  }

  std::cout << "Best plane is " << bestPlane << std::endl;
  std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
  std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
  std::cout << "The shower start is " << std::endl;
  vertex.Print();

  return;

}

void shower::EMShowerAlg::ProjectVertexIn2D(TVector3 const& vertex,
					    std::map<int,std::vector<art::Ptr<recob::Hit> > >& trackHitsMap,
					    std::map<int,std::pair<std::vector<art::Ptr<recob::Hit> >,std::vector<art::Ptr<recob::Hit> > > > const& trackHitsBothEndsMap) {

  /// Projects the 3D direction into all the 2D views to make sure the correct end of the track is selected in each view

  for (std::map<int,std::pair<std::vector<art::Ptr<recob::Hit> >, std::vector<art::Ptr<recob::Hit> > > >::const_iterator trackHitsIt = trackHitsBothEndsMap.begin();
       trackHitsIt != trackHitsBothEndsMap.end();
       ++trackHitsIt) {

    TVector2 end1 = HitCoordinates(*trackHitsIt->second.first.begin());
    TVector2 end2 = HitCoordinates(*trackHitsIt->second.second.begin());

    TVector2 vertexProj = Project3DPointOntoPlane(vertex, trackHitsIt->first);

    if ( (vertexProj - end1).Mod() < (vertexProj - end2).Mod() )
      trackHitsMap[trackHitsIt->first] = trackHitsIt->second.first;
    else
      trackHitsMap[trackHitsIt->first] = trackHitsIt->second.second;

  }

}

void shower::EMShowerAlg::FindShowerStartDirection(art::Ptr<recob::Track> const& vertexTrack,
						   std::map<int,TVector2> const& showerCentreMap,
						   TVector3& showerVertex,
						   TVector3& showerDirection) {

  /// Finds the start of the shower, and its direction

  TVector3 vertex = vertexTrack->Vertex();
  TVector3 end = vertexTrack->End();

  std::map<int,double> distanceToVertex;
  std::map<int,double> distanceToEnd;

  // Loop over all the planes and find the distance from the vertex and end projections to the centre in each plane
  for (std::map<int,TVector2>::const_iterator showerCentreIt = showerCentreMap.begin(); showerCentreIt != showerCentreMap.end(); ++showerCentreIt) {

    // Project the vertex and the end point onto this plane
    TVector2 vertexProj = Project3DPointOntoPlane(vertex, showerCentreIt->first);
    TVector2 endProj    = Project3DPointOntoPlane(end, showerCentreIt->first);

    // Find the distance of each to the centre of the cluster
    distanceToVertex[showerCentreIt->first] = (vertexProj - showerCentreIt->second).Mod();
    distanceToEnd[showerCentreIt->first] = (endProj - showerCentreIt->second).Mod();

  }

  // Find the average distance to the vertex and the end across the planes
  double avDistanceToVertex = 0, avDistanceToEnd = 0;
  for (std::map<int,double>::iterator distanceToVertexIt = distanceToVertex.begin(); distanceToVertexIt != distanceToVertex.end(); ++distanceToVertexIt)
    avDistanceToVertex += distanceToVertexIt->second;
  avDistanceToVertex /= distanceToVertex.size();

  for (std::map<int,double>::iterator distanceToEndIt = distanceToEnd.begin(); distanceToEndIt != distanceToEnd.end(); ++distanceToEndIt)
    avDistanceToEnd += distanceToEndIt->second;
  avDistanceToEnd /= distanceToEnd.size();

  // Set the vertex and directions for this shower
  if (avDistanceToVertex > avDistanceToEnd) {
    showerVertex = vertex;
    showerDirection = vertexTrack->VertexDirection();
  }
  else {
    showerVertex = end;
    showerDirection = (-1) * vertexTrack->VertexDirection();
  }

  return;

}

// void shower::EMShowerAlg::FindShowerCentre(std::map<int,TVector2> const& showerCentreMap, TVector3& showerCentre) {

//   /// Takes a map of shower centres from all the views and finds the rough 3D point

//   // Find all the x and wire positions
//   std::vector<double> xPos, wires;
//   for (std::map<int,TVector2>::const_iterator showerCentreIt = showerCentreMap.begin(); showerCentreIt != showerCentreMap.end(); ++showerCentreIt) {
//     xPos.push_back(fDetProp->ConvertTicksToX(showerCentreIt->second.Y()));
//     wires.push_back(showerCentreIt->first);
//   }

//   // Find the average x position
//   double xAv;
//   for (std::vector<double>::iterator xIt = xPos.begin(); xIt != xPos.end(); ++xIt)
//     xAv += *xIt;
//   xAv /= xPos.size();

//   // Find the y and z coordinates
//   if (wires.size() == 2) {
    

// }

art::Ptr<recob::Track> shower::EMShowerAlg::FindVertexTrack(std::map<int,art::Ptr<recob::Hit> > const& vertexMap,
							    std::map<int,art::Ptr<recob::Track> > const& trackMap,
							    std::map<int,std::vector<art::Ptr<recob::Hit> > > const& trackHitsMap) {

  /// Finds the 3D vertex and direction given the tracks associated with the 2D vertex

  art::Ptr<recob::Track> vertexTrack;

  // First, find out if two views agree on a track (normally they will)
  std::map<int,std::vector<int> > trackIDToPlanes;
  std::vector<int> planesWithTrack;
  for (std::map<int,art::Ptr<recob::Track> >::const_iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
    if (!trackIt->second.isNull()) {
      trackIDToPlanes[trackIt->second->ID()].push_back(trackIt->first);
      planesWithTrack.push_back(trackIt->first);
    }
  for (std::map<int,std::vector<int> >::iterator trackVertexIt = trackIDToPlanes.begin(); trackVertexIt != trackIDToPlanes.end(); ++trackVertexIt)
    if (trackVertexIt->second.size() > 1)
      vertexTrack = trackMap.at(*trackVertexIt->second.begin());

  // If there are no planes with a track reconstructed, don't carry on right now
  if (planesWithTrack.size() == 0)
    return vertexTrack;

  // If they don't, try to use the third view (it exsits) to pick the correct track
  if (vertexTrack.isNull() and trackMap.size() > 2) {
    std::map<int,double> distanceOfVertexFromVertex;
    for (std::map<int,art::Ptr<recob::Hit> >::const_iterator vertexIt = vertexMap.begin(); vertexIt != vertexMap.end(); ++vertexIt) {
      if (std::find(planesWithTrack.begin(), planesWithTrack.end(), vertexIt->first) != planesWithTrack.end())
	continue;
      for (std::vector<int>::iterator planesWithTrackIt = planesWithTrack.begin(); planesWithTrackIt != planesWithTrack.end(); ++planesWithTrackIt)
	distanceOfVertexFromVertex[*planesWithTrackIt] = ( TMath::Abs( HitCoordinates(vertexIt->second).Y() - HitCoordinates(vertexMap.at(*planesWithTrackIt)).Y() ) ) / HitCoordinates(vertexMap.at(*planesWithTrackIt)).Y();
    }
    std::vector<int> thirdViewClose;
    for (std::map<int,double>::iterator distanceOfVertexFromVertexIt = distanceOfVertexFromVertex.begin(); distanceOfVertexFromVertexIt != distanceOfVertexFromVertex.end(); ++distanceOfVertexFromVertexIt)
      if (distanceOfVertexFromVertexIt->second < 0.1)
	thirdViewClose.push_back(distanceOfVertexFromVertexIt->first);
    if (thirdViewClose.size() == 1)
      vertexTrack = trackMap.at(thirdViewClose.at(0));
  }

  // Finally, if all else fails, just pick the view with the longest initial track reconstructed in 2D
  if (vertexTrack.isNull()) {
    std::map<int,int> lengthOfTrackToPlane;
    for (std::vector<int>::iterator planesWithTrackIt = planesWithTrack.begin(); planesWithTrackIt != planesWithTrack.end(); ++planesWithTrackIt)
      lengthOfTrackToPlane[trackHitsMap.at(*planesWithTrackIt).size()] = *planesWithTrackIt;
    vertexTrack = trackMap.at(lengthOfTrackToPlane.rbegin()->second);
  }

  return vertexTrack;

}

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindTrack(std::vector<art::Ptr<recob::Hit> > const& shower, TVector2 const& start, TVector2 const& end) {

  /// Finds the track from the start of the shower

  std::vector<art::Ptr<recob::Hit> > trackHits;

  // Map of hit on each wire
  std::map<int,std::vector<art::Ptr<recob::Hit> > > hitWires;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
    hitWires[(int)HitCoordinates(*hit).X()].push_back(*hit);

  // Call the track ended when there is more than one hit on a wire twice in three wires
  std::vector<int> lastTwoWires = {0, 0};

  // Find out which way to look!
  int startWire = (int)start.X(), endWire = (int)end.X();
  bool increasing;
  if (startWire < endWire) increasing = true;
  else increasing = false;

  // Look through the hits from the start
  if (increasing) {
    for (int wire = startWire; wire < endWire; ++wire) {
      int numberOfHitsOnThisWire;
      if (hitWires.find(wire) != hitWires.end()) numberOfHitsOnThisWire = hitWires.at(wire).size();
      else numberOfHitsOnThisWire = 0;
      if (numberOfHitsOnThisWire >= 2 and (lastTwoWires.at(0) >= 2 or lastTwoWires.at(1) >= 2))
	break;
      else {
	if (numberOfHitsOnThisWire)
	  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitWireIt = hitWires.at(wire).begin(); hitWireIt != hitWires.at(wire).end(); ++hitWireIt)
	    trackHits.push_back(*hitWireIt);
	lastTwoWires[0] = lastTwoWires[1];
	lastTwoWires[1] = numberOfHitsOnThisWire;
      }
    }
  }

  else {
    for (int wire = startWire; wire > endWire; --wire) {
      int numberOfHitsOnThisWire;
      if (hitWires.find(wire) != hitWires.end()) numberOfHitsOnThisWire = hitWires.at(wire).size();
      else numberOfHitsOnThisWire = 0;
      if (numberOfHitsOnThisWire >= 2 and (lastTwoWires.at(0) >= 2 or lastTwoWires.at(1) >= 2))
	break;
      else {
	if (numberOfHitsOnThisWire)
	  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitWireIt = hitWires.at(wire).begin(); hitWireIt != hitWires.at(wire).end(); ++hitWireIt)
	    trackHits.push_back(*hitWireIt);
	lastTwoWires[0] = lastTwoWires[1];
	lastTwoWires[1] = numberOfHitsOnThisWire;
      }
    }
  }

  if (trackHits.size() == 0)
    trackHits.push_back(hitWires.begin()->second.at(0));

  return trackHits;

}

void shower::EMShowerAlg::FindShowerEnds(std::vector<art::Ptr<recob::Hit> > const& shower,
					 TVector2 const& centre2,
					 art::Ptr<recob::Hit>& end1,
					 art::Ptr<recob::Hit>& end2) {

  /// Roughly finds the two 'ends' of the shower (one is the vertex, one isn't well defined!)

  // Find the charge-weighted centre (in cm) of this shower
  TVector2 pos, chargePoint = TVector2(0,0);
  double totalCharge = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    pos = HitPosition(*hit);
    chargePoint += (*hit)->Integral() * pos;
    totalCharge += (*hit)->Integral();
  }
  TVector2 centre = chargePoint / totalCharge;

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
  TVector2 direction = TVector2(1,gradient).Unit();

  // Find how far each hit (projected onto this axis) is from the centre
  std::map<double,art::Ptr<recob::Hit> > hitProjection;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    pos = HitPosition(*hit) - centre;
    hitProjection[direction*pos] = *hit;
  }

  end1 = hitProjection.begin()->second;
  end2 = hitProjection.rbegin()->second;

  return;

}

double shower::EMShowerAlg::FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, art::Ptr<recob::Track> const& track, unsigned int plane, geo::View_t const& view) {

  /// Finds dE/dx for the track given a set of hits

  double totalCharge = 0, totalDistance = 0, avHitTime = 0;
  unsigned int nHits = 0;

  // Get the pitch
  double pitch = 0;
  try { pitch = track->PitchInView(view); }
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
  double dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avHitTime, plane);

  return dEdx;

}

// void shower::EMShowerAlg::FindVertex(art::PtrVector<recob::Hit> const& shower, TVector2 const& end1, TVector2 const& end2, std::vector<int>& trackHits) {

//   /// Decides which 'end' of the shower is the true vertex and finds the initial track

//   // Map of hit on each wire
//   std::map<int,std::vector<int> > hitWires;
//   for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
//     hitWires[(int)HitCoordinates(*hit).X()].push_back(std::distance(shower.begin(), hit));

//   std::vector<int> trackHits1, trackHits2;
//   FindTrack(end1, end2, hitWires, trackHits1);
//   FindTrack(end2, end1, hitWires, trackHits2);

//   trackHits = trackHits1.size() > trackHits2.size() ? trackHits1 : trackHits2;

//   std::cout << "Vertex is (" << HitCoordinates(shower.at(*trackHits.begin())).X() << ", " << HitCoordinates(shower.at(*trackHits.begin())).Y() << ") and the end of the track is (" << HitCoordinates(shower.at(*trackHits.rbegin())).X() << ", " << HitCoordinates(shower.at(*trackHits.rbegin())).Y() << ")" << std::endl;

//   return;

// }

// void shower::EMShowerAlg::FindTrack(art::PtrVector<recob::Hit> const& shower, std::map<double,int> const& hitToEnd, std::vector<int>& trackHits) {

//   /// Contructs a track from the 'end' of a shower

//   TVector2 end = HitCoordinates(shower.at(hitToEnd.begin()->second));
//   std::cout << "End coordinates are " << std::endl;
//   end.Print();

//   // Look through the hits from the end
//   for (std::map<double,int>::const_iterator hitToEndIt = hitToEnd.begin(); hitToEndIt != hitToEnd.end(); ++hitToEndIt) {

//     // Put first hit in vector
//     if (hitToEnd.begin() == hitToEndIt) {
//       trackHits.push_back(hitToEndIt->second);
//       continue;
//     }

//     // Find the direction from all hits up to this one
//     double nhits, sumx, sumy, sumx2, sumxy;
//     for (std::map<double,int>::const_iterator hitsInTrackIt = hitToEnd.begin(); hitsInTrackIt != std::next(hitToEndIt,1); ++hitsInTrackIt) {
//       ++nhits;
//       TVector2 hitpos = HitCoordinates(shower.at(hitsInTrackIt->second));
//       sumx += hitpos.X();
//       sumy += hitpos.Y();
//       sumx2 += hitpos.X() * hitpos.X();
//       sumxy += hitpos.X() * hitpos.Y();
//     }
//     double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
//     TVector2 direction = TVector2(1,gradient).Unit();

//     // Project all the hits onto this direction and find average displacement
//     double totalDisplacement = 0;
//     for (std::map<double,int>::const_iterator hitsInTrackIt = hitToEnd.begin(); hitsInTrackIt != std::next(hitToEndIt,1); ++hitsInTrackIt) {
//       TVector2 hitpos = HitCoordinates(shower.at(hitsInTrackIt->second)) - end;
//       totalDisplacement += TMath::Sqrt(TMath::Power(hitpos.Mod(),2) - TMath::Power((hitpos.Proj(direction)).Mod(),2));
//     }
//     double averageDisplacement = totalDisplacement / (std::distance(hitToEnd.begin(), hitToEndIt) + 1);

//     std::cout << "Hit at (" << HitCoordinates(shower.at(hitToEndIt->second)).X() << ", " << HitCoordinates(shower.at(hitToEndIt->second)).Y() << ") [distance of " << hitToEndIt->first << "] leads to av displacement of " << averageDisplacement << std::endl;

//     // if (averageDisplacement < 1)
//     //   trackHits.push_back(hitToEndIt->second);
//     // else
//     //   break;

//   }

//   return;

// }

// void shower::EMShowerAlg::FindVertex(art::PtrVector<recob::Hit> const& shower, TVector2 const& end1, TVector2 const& end2, std::vector<int>& trackHits) {

//   /// Decides which 'end' of the shower is the true vertex and finds the initial track

//   // Distance of each hit from either 'end'
//   std::map<double,int> hitToEnd1, hitToEnd2;
//   TVector2 hitpos;
//   for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
//     hitpos = HitCoordinates(*hit);
//     hitToEnd1[(hitpos - end1).Mod()] = std::distance(shower.begin(), hit);
//     hitToEnd2[(hitpos - end2).Mod()] = std::distance(shower.begin(), hit);
//   }

//   // Try to form a straight line from each end to see which is the vertex
//   std::vector<int> trackHits1, trackHits2;
//   FindTrack(shower, hitToEnd1, trackHits1);
//   FindTrack(shower, hitToEnd2, trackHits2);

//   // The vertex is the end with the longest initial straight track
//   trackHits = trackHits1.size() > trackHits2.size() ? trackHits1 : trackHits2;

//   return;

// }

// void shower::EMShowerAlg::FindTrack(art::PtrVector<recob::Hit> const& shower, std::map<double,int> const& hitToEnd, std::vector<int>& trackHits) {

//   /// Contructs a track from the 'end' of a shower

//   TVector2 end = HitCoordinates(shower.at(hitToEnd.begin()->second));
//   std::cout << "End coordinates are " << std::endl;
//   end.Print();

//   // Look through the hits from the end
//   for (std::map<double,int>::const_iterator hitToEndIt = hitToEnd.begin(); hitToEndIt != hitToEnd.end(); ++hitToEndIt) {

//     // Put first hit in vector
//     if (hitToEnd.begin() == hitToEndIt) {
//       trackHits.push_back(hitToEndIt->second);
//       continue;
//     }

//     // Find the direction from all hits up to this one
//     double nhits, sumx, sumy, sumx2, sumxy;
//     for (std::map<double,int>::const_iterator hitsInTrackIt = hitToEnd.begin(); hitsInTrackIt != std::next(hitToEndIt,1); ++hitsInTrackIt) {
//       ++nhits;
//       TVector2 hitpos = HitCoordinates(shower.at(hitsInTrackIt->second));
//       sumx += hitpos.X();
//       sumy += hitpos.Y();
//       sumx2 += hitpos.X() * hitpos.X();
//       sumxy += hitpos.X() * hitpos.Y();
//     }
//     double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
//     TVector2 direction = TVector2(1,gradient).Unit();

//     // Project all the hits onto this direction and find average displacement
//     double totalDisplacement = 0;
//     for (std::map<double,int>::const_iterator hitsInTrackIt = hitToEnd.begin(); hitsInTrackIt != std::next(hitToEndIt,1); ++hitsInTrackIt) {
//       TVector2 hitpos = HitCoordinates(shower.at(hitsInTrackIt->second)) - end;
//       totalDisplacement += TMath::Sqrt(TMath::Power(hitpos.Mod(),2) - TMath::Power((hitpos.Proj(direction)).Mod(),2));
//     }
//     double averageDisplacement = totalDisplacement / (std::distance(hitToEnd.begin(), hitToEndIt) + 1);

//     std::cout << "Hit at (" << HitCoordinates(shower.at(hitToEndIt->second)).X() << ", " << HitCoordinates(shower.at(hitToEndIt->second)).Y() << ") [distance of " << hitToEndIt->first << "] leads to av displacement of " << averageDisplacement << std::endl;

//     // if (averageDisplacement < 1)
//     //   trackHits.push_back(hitToEndIt->second);
//     // else
//     //   break;

//   }

//   return;

// }

// art::Ptr<recob::Hit> shower::EMShowerAlg::FindVertex(art::PtrVector<recob::Hit> const& shower, art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2) {

//   /// Decides which 'end' of the shower is the true vertex

//   TVector2 end1 = HitCoordinates(hit1);
//   TVector2 end2 = HitCoordinates(hit2);

//   // Find the geometrical centre of the shower and the direction vector
//   TVector2 centre = (end1 + end2) / 2;
//   TVector2 direction = (end1 - end2).Unit();

//   // Sum all deposited charge for each end of the shower
//   double chargeEnd1 = 0, chargeEnd2 = 0;

//   // Project all hits onto this vector to determine which end each is closer to
//   TVector2 proj;
//   double distanceEnd1, distanceEnd2;
//   for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
//     proj = (HitCoordinates(*hit)-centre).Proj(direction);
//     distanceEnd1 = (end1 - centre - proj).Mod();
//     distanceEnd2 = (end2 - centre - proj).Mod();
//     if (distanceEnd1 < distanceEnd2)
//       chargeEnd1 += (*hit)->Integral();
//     else
//       chargeEnd2 += (*hit)->Integral();
//   }

//   std::cout << "Charge nearest end1 is " << chargeEnd1 << " and end2 is " << chargeEnd2 << std::endl;

//   // The half of the shower nearest the vertex will deposit less energy
//   return chargeEnd1 < chargeEnd2 ? hit1 : hit2;

// }

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

  double globalWire;
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }
  else {
    if (wireID.TPC % 2 == 0) globalWire = wireID.Wire + ((wireID.TPC/2) * fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat));
    else globalWire = wireID.Wire + ((int)(wireID.TPC/2) * fGeom->Nwires(wireID.Plane, 1, wireID.Cryostat));
  }

  return globalWire;

}

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, unsigned int plane) {

  /// Projects a 3D point onto a 2D plane

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  // temp fix for tpc for lariat!
  geo::TPCID tpcID = fGeom->FindTPCAtPosition(pointPosition);
  int tpc = 0;
  if (tpcID.isValid)
    tpc = tpcID.TPC;
  else
    tpc = 0;

  return TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), plane, tpc % 2, 0),
		  fDetProp->ConvertXToTicks(point.X(), plane, tpc % 2, 0));

}

// int shower::EMShowerAlg::FindTrackID(art::Ptr<recob::Hit> const& hit) {
//   double particleEnergy = 0;
//   int likelyTrackID = 0;
//   std::vector<sim::TrackIDE> trackIDs = bt->HitToTrackID(hit);
//   for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
//     if (trackIDs.at(idIt).energy > particleEnergy) {
//       particleEnergy = trackIDs.at(idIt).energy;
//       likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
//     }
//   }
//   return likelyTrackID;
// }

// int shower::EMShowerAlg::FindTrueTrack(std::vector<art::Ptr<recob::Hit> > const& showerHits) {
//   std::map<int,double> trackMap;
//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator showerHitIt = showerHits.begin(); showerHitIt != showerHits.end(); ++showerHitIt) {
//     art::Ptr<recob::Hit> hit = *showerHitIt;
//     int trackID = FindTrackID(hit);
//     trackMap[trackID] += hit->Integral();
//   }
//   //return std::max_element(trackMap.begin(), trackMap.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) {return p1.second < p2.second;} )->first;
//   double highestCharge = 0;
//   int clusterTrack = 0;
//   for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
//     if (trackIt->second > highestCharge) {
//       highestCharge = trackIt->second;
//       clusterTrack  = trackIt->first;
//     }
//   return clusterTrack;
// }

