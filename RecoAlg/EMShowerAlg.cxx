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

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
								  std::vector<art::Ptr<recob::Hit> > const& track2,
								  std::map<int,TVector2> const& showerCentreMap) {

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).
  /// This implementation also orients the track in the correct direction if a map of shower centres (units [wire/tick]) in each view is provided.

  std::unique_ptr<recob::Track> track;

  // Check for tracks crossing TPC boundaries
  std::map<int,int> tpcMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = track1.begin(); hitIt != track1.end(); ++hitIt)
    ++tpcMap[(*hitIt)->WireID().TPC];
  if (tpcMap.size() > 1) {
    mf::LogWarning("EMShowerAlg") << "Warning: attempting to construct a track which crosses more than one TPC -- PMTrack can't handle this right now.  Returning a null track.";
    return track;
  }

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
  std::map<int,double> distanceToVertex, distanceToEnd;
  TVector3 vertex = *xyz.begin(), end = *xyz.rbegin();

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

  std::map<int,TVector2> showerCentreMap;

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

void shower::EMShowerAlg::FindInitialTrack(art::PtrVector<recob::Hit> const& hits,
					   std::unique_ptr<recob::Track>& initialTrack,
					   std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits, int plane) {

  /// Finds the initial track-like part of the shower and the hits in all views associated with it

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  std::map<int,std::vector<art::Ptr<recob::Hit> > > orderedShowerMap;
  std::map<int,TVector2> showerCentreMap;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::cout << std::endl << "Plane " << planeHits->first << std::endl;
    if (planeHits->first != plane and plane != -1) continue;

    // Find the charge-weighted centre (in [wire/tick]) of this shower
    TVector2 pos, chargePoint = TVector2(0,0);
    double totalCharge = 0;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = planeHits->second.begin(); hit != planeHits->second.end(); ++hit) {
      pos = HitCoordinates(*hit);
      chargePoint += (*hit)->Integral() * pos;
      totalCharge += (*hit)->Integral();
    }
    TVector2 centre = chargePoint / totalCharge;
    showerCentreMap[planeHits->first] = centre;

    std::vector<art::Ptr<recob::Hit> > showerHits = OrderShowerHits(planeHits->second);
    orderedShowerMap[planeHits->first] = showerHits;

    std::cout << "After ordering: Plane " << planeHits->first << ": start is determined to be (" << HitCoordinates(*showerHits.begin()).X() << ", " << HitCoordinates(*showerHits.begin()).Y() << ") and the end is determined to be (" << HitCoordinates(*showerHits.rbegin()).X() << ", " << HitCoordinates(*showerHits.rbegin()).Y() << ")" << std::endl;

  }

  std::cout << std::endl;

  // Now find the hits belonging to the track
  initialTrackHits = FindShowerStart(orderedShowerMap);

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator trackHitIt = initialTrackHits.begin(); trackHitIt != initialTrackHits.end(); ++trackHitIt)
    std::cout << "Plane " << trackHitIt->first << " has track start (" << HitCoordinates(*trackHitIt->second.begin()).X() << ", " << HitCoordinates(*trackHitIt->second.begin()).Y() << ") and end (" << HitCoordinates(*trackHitIt->second.rbegin()).X() << ", " << HitCoordinates(*trackHitIt->second.rbegin()).Y() << ")" << std::endl;

  // Now we have the track hits -- can make a track!
  initialTrack = MakeInitialTrack(initialTrackHits, showerCentreMap);

  if (initialTrack)
    std::cout << "Found track! Start is (" << initialTrack->Vertex().X() << ", " << initialTrack->Vertex().Y() << ", " << initialTrack->Vertex().Z() << ") and end is (" << initialTrack->End().X() << ", " << initialTrack->End().Y() << ", " << initialTrack->End().Z() << "); direction is (" << initialTrack->VertexDirection().X() << ", " << initialTrack->VertexDirection().Y() << ", " << initialTrack->VertexDirection().Z() << ")" << std::endl;

  return;

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
      std::map<int,std::vector<art::Ptr<recob::Hit> > > wireHitMap;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt) {
	int wire = std::round(HitCoordinates(*hitIt).X());
	wireHitMap[wire].push_back(*hitIt);
	if ( wireHitMap[wire].size() > 1 and ((wireHitMap[wire-1].size() > 1 or wireHitMap[wire-2].size() > 1)) )
	  break;
	else
	  initialHits.push_back(*hitIt);
      }
    }

    // Deal with showers travelling along ticks
    else {
      std::map<int,std::vector<art::Ptr<recob::Hit> > > tickHitMap;
      for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShower.begin(); hitIt != orderedShower.end(); ++hitIt) {
	int tick = std::round(HitCoordinates(*hitIt).Y());
	tickHitMap[tick].push_back(*hitIt);
	if ( (tickHitMap[tick].size() > 1) and ((tickHitMap[tick-1].size() > 1 or tickHitMap[tick-2].size() > 1)) )
	  break;
	else
	  initialHits.push_back(*hitIt);
      }
    }

    initialHitsMap[orderedShowerIt->first] = initialHits;

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
      else
	++hitsAfterOffshoot;
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

std::unique_ptr<recob::Track> shower::EMShowerAlg::MakeInitialTrack(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& initialHitsMap,
								    std::map<int,TVector2> const& showerCentreMap) {

  /// Takes initial track hits from multiple views and forms a track object which best represents the start of the shower
  /// Returns false if this isn't possible.

  // Don't want to make any assumptions on the number of planes, or the number of planes with hits on
  // Try to resolve any conflicts as we go...

  std::unique_ptr<recob::Track> track;

  // If there's only one plane, there's not much we can do!
  if (initialHitsMap.size() == 1)
    return track;

  // If there are two planes then we can make the object
  else if (initialHitsMap.size() == 2) {
    std::vector<std::vector<art::Ptr<recob::Hit> > > showerPlanes;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerPlaneIt = initialHitsMap.begin(); showerPlaneIt != initialHitsMap.end(); ++showerPlaneIt)
      showerPlanes.push_back(showerPlaneIt->second);
    track = ConstructTrack(showerPlanes.at(0), showerPlanes.at(1), showerCentreMap);
  }

  // If we have three planes to consider then we can do more checks!
  else if (initialHitsMap.size() == 3) {
    std::map<int,std::vector<art::Ptr<recob::Hit> > > trackSizeMap;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator showerPlaneIt = initialHitsMap.begin(); showerPlaneIt != initialHitsMap.end(); ++showerPlaneIt)
      trackSizeMap[showerPlaneIt->second.size()] = showerPlaneIt->second;
    std::vector<std::vector<art::Ptr<recob::Hit> > > showerPlanes;
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::reverse_iterator trackSizeIt = trackSizeMap.rbegin(); trackSizeIt != trackSizeMap.rend(); ++trackSizeIt) {
      if (std::distance(trackSizeMap.rbegin(), trackSizeIt) > 1)
	break;
      showerPlanes.push_back(trackSizeIt->second);
    }
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

std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> > const& shower) {

  /// Takes the hits associated with a shower and orders them so they follow the direction of the shower

  // Find the charge-weighted centre (in [cm]) of this shower
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

  // Get a vector of hits in order of the shower
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(hitProjection.begin(), hitProjection.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });
  
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

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, unsigned int plane) {

  /// Projects a 3D point (units [cm]) onto a 2D plane
  /// Returns 2D point (units [wire/tick])

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  return TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), plane, fGeom->FindTPCAtPosition(pointPosition).TPC % 2, 0),
		  fDetProp->ConvertXToTicks(point.X(), plane, fGeom->FindTPCAtPosition(pointPosition).TPC % 2, 0));

}
