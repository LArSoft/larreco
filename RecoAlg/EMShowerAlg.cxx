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


void shower::EMShowerAlg::FindInitialTrack(art::PtrVector<recob::Hit> const& hits, recob::Track& initialTrack, std::vector<art::Ptr<recob::Hit> >& initialTrackHits) {

  /// Finds the initial track-like part of the shower and the hits in all views associated with it

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  std::map<int,std::vector<art::Ptr<recob::Hit> > > orderedShowerMap;
  std::map<int,double> goodnessOfOrderMap;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::cout << std::endl << "Plane " << planeHits->first << std::endl;
    if (planeHits->first != 0) continue;

    // Find the charge-weighted centre (in cm) of this shower
    TVector2 pos, chargePoint = TVector2(0,0);
    double totalCharge = 0;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = planeHits->second.begin(); hit != planeHits->second.end(); ++hit) {
      pos = HitPosition(*hit);
      chargePoint += (*hit)->Integral() * pos;
      totalCharge += (*hit)->Integral();
    }
    TVector2 centre = chargePoint / totalCharge;

    // Find a rough shower 'direction'
    int nhits = 0;
    double sumx=0., sumy=0., sumx2=0., sumxy=0.;
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = planeHits->second.begin(); hit != planeHits->second.end(); ++hit) {
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
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = planeHits->second.begin(); hit != planeHits->second.end(); ++hit) {
      pos = HitPosition(*hit) - centre;
      hitProjection[direction*pos] = *hit;
    }

    // Get a vector of hits in order of the shower
    std::vector<art::Ptr<recob::Hit> > showerHits;
    std::transform(hitProjection.begin(), hitProjection.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

    // Find the shower direction properties
    std::vector<double> propertiesFromBeginning = GetShowerDirectionProperties(showerHits, direction);
    std::reverse(showerHits.begin(), showerHits.end());
    //std::vector<double> propertiesFromEnd = GetShowerDirectionProperties(showerHits, direction);
    std::vector<double> propertiesFromEnd = {0.,0.,0.};

    // Decide which is the end
    double goodnessOfOrder = 0;
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

    orderedShowerMap[planeHits->first] = showerHits;
    goodnessOfOrderMap[planeHits->first] = goodnessOfOrder;

    std::cout << "After ordering: Plane " << planeHits->first << ": start is determined to be (" << HitCoordinates(*showerHits.begin()).X() << ", " << HitCoordinates(*showerHits.begin()).Y() << ") and the end is determined to be (" << HitCoordinates(*showerHits.rbegin()).X() << ", " << HitCoordinates(*showerHits.rbegin()).Y() << ") -- goodness is " << goodnessOfOrder << std::endl;

  }

  std::cout << std::endl;

  // // Make the object using all this information
  // recob::Track initialTrack;
  // std::vector<art::Ptr<recob::Hit> > trackHits;
  // MakeVertexTrack(initialTrack, trackHits, orderedShowerMap, goodnessOfOrderMap);

  // Now find the hits belonging to the track
  std::map<int,std::vector<art::Ptr<recob::Hit> > > trackHits = FindShowerStart(orderedShowerMap);

  return;

  // Don't want to make any assumptions on the number of planes, or the number of planes with hits on
  // Try to resolve any conflicts as we go...

  // If there's only one plane, there's not much we can do!
  if (planeHitsMap.size() == 1)
    return;

  // If there are two planes then we can try

}

// void shower::EMShowerAlg::MakeVertexTrack(recob::Track& initialTrack,
// 					  std::vector<art::Ptr<recob::Hit> >& trackHits,
// 					  std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap,
// 					  std::map<int,double> const& goodnessOfOrderMap) {

//   /// Makes a track to represent the start of the shower given the ordered shower hits

//   // Don't want to make any assumptions on the number of planes present

// }

std::vector<double> shower::EMShowerAlg::GetShowerDirectionProperties(std::vector<art::Ptr<recob::Hit> > const& showerHits, TVector2 const& direction) {

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
    //pos = HitCoordinates(*hitIt) - start;
    graph->SetPoint(graph->GetN(),pos.X(),pos.Y());
    double angle = pos.DeltaPhi(direction);
    projPos = (pos).Proj(direction);
    double distanceFromAxis = TMath::Sqrt(TMath::Power(pos.Mod(),2) - TMath::Power(projPos.Mod(),2));
    //std::cout << "X is " << pos.X() << " and distance from axis is " << distanceFromAxis << std::endl;
    //std::cout << "Hit (" << HitCoordinates(*hitIt).X() << ", " << HitCoordinates(*hitIt).Y() << ") has (relative position is (" << pos.X() << ", " << pos.Y() << ")); angle " << angle << " and projected distance " << distanceFromAxis << std::endl;
    if (angle > 0 and distanceFromAxis > 1) ++hitsAbove;
    if (angle < 0 and distanceFromAxis > 1) ++hitsBelow;
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
    //std::cout << "End of hit; spreadingHits is " << spreadingHits << " and offshoot is " << hitsAfterOffshoot << std::endl;
  }

  //std::cout << "Total hits away from centre is " << hitsAbove+hitsBelow << std::endl;
  std::vector<double> directionProperties = {(double)spreadingHits, hitsAbove/(double)showerHits.size(), hitsBelow/(double)showerHits.size()};

  std::cout << "Direction properties above " << directionProperties.at(1) << ", below " << directionProperties.at(2) << std::endl;

  TLine line;
  line.SetLineColor(2);
  canv->cd();
  line.DrawLine(-1000*direction.X(),-1000*direction.Y(),1000*direction.X(),1000*direction.Y());

  canv->SaveAs("ShowerPlane.png");

  return directionProperties;

}

std::map<int,std::vector<art::Ptr<recob::Hit> > > shower::EMShowerAlg::FindShowerStart(std::map<int,std::vector<art::Ptr<recob::Hit> > > const& orderedShowerMap) {

  /// Takes a map of the shower hits on each plane (ordered from what has been decided to be the start) along with some metric of goodness
  /// Returns a map of the initial track-like part of the shower on each plane

  std::map<int,std::vector<art::Ptr<recob::Hit> > > initialHitsMap;

  for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator orderedShowerIt = orderedShowerMap.begin(); orderedShowerIt != orderedShowerMap.end(); ++orderedShowerIt) {

    std::vector<art::Ptr<recob::Hit> > initialHits;

    double previousGradient = 0;

    // Look at the shower from the start and decide which hits belong to the initial shower
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = orderedShowerIt->second.begin(); hitIt != orderedShowerIt->second.end(); ++hitIt) {
      if (std::distance(orderedShowerIt->second.begin(),hitIt) == 0) {
	initialHits.push_back(*hitIt);
	continue;
      }
      if (std::distance(orderedShowerIt->second.begin(),hitIt) == 1) {
	initialHits.push_back(*hitIt);
	previousGradient = (HitPosition(initialHits.at(0))-HitPosition(initialHits.at(1))).Y() / (double)(HitPosition(initialHits.at(0))-HitPosition(initialHits.at(1))).X();
	continue;
      }
      double currentGradient = (HitPosition(*initialHits.rbegin())-HitPosition(*hitIt)).Y() / (double)(HitPosition(*initialHits.rbegin())-HitPosition(*hitIt)).X();
      double gradientFractionalChange = TMath::Abs(1 - currentGradient/previousGradient);
      if (gradientFractionalChange < 0.1) {
	initialHits.push_back(*hitIt);
	previousGradient = currentGradient;
      }
      else
	break;
    }

    initialHitsMap[orderedShowerIt->first] = initialHits;

  }

  return initialHitsMap;

}


//     // Find the track-like part of the shower
//     std::vector<art::Ptr<recob::Hit> > trackHits = FindShowerStart(planeHits->second);
//     //trackMap[planeHits->first] = track;
//     trackHitsMap[planeHits->first] = trackHits;

//     // Find a track object
//     art::Ptr<recob::Track> track;
//     for (std::vector<art::Ptr<recob::Hit> >::iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
//       std::vector<art::Ptr<recob::Track> > tracks = fmt.at(trackHitIt->key());
//       if (tracks.size() == 0)
//     	continue;
//       track = tracks.at(0);
//       break;
//     }
//     trackMap[planeHits->first] = track;

//     // Vertex
//     art::Ptr<recob::Hit> planevertex = *trackHits.begin();
//     vertexMap[planeHits->first] = planevertex;
//     std::cout << "Vertex is (" << HitCoordinates(planevertex).X() << ", " << HitCoordinates(planevertex).Y() << ")" << std::endl;

//   }

//   // Made a 3D track object for the start of the shower
//   trackHitsMap = CheckTrackHits(trackHitsMap);
//   art::Ptr<recob::Track> vertexTrack = MakeVertexTrack(trackHitsMap);

//   // // Find 3D vertex and direction
//   // art::Ptr<recob::Track> vertexTrack = FindVertexTrack(vertexMap, trackMap, trackHitsMap);
//   // if (vertexTrack.isNull()) return;
//   FindShowerStartDirection(vertexTrack, showerCentreMap, vertex, direction);

//   // Use the 3D vertex to resolve any problems in 2D
//   ProjectVertexIn2D(vertex, trackHitsMap, trackHitsBothEndsMap);

//   return;

// }

// std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindInitialTrack(std::vector<art::Ptr<recob::Hit> > const& shower) {

//   /// Finds the initial track like part of the shower in a particular view

//   // Find the charge-weighted centre (in cm) of this shower
//   TVector2 pos, chargePoint = TVector2(0,0);
//   double totalCharge = 0;
//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
//     pos = HitPosition(*hit);
//     chargePoint += (*hit)->Integral() * pos;
//     totalCharge += (*hit)->Integral();
//   }
//   TVector2 centre = chargePoint / totalCharge;

//   // Find a rough shower 'direction'
//   int nhits = 0;
//   double sumx=0., sumy=0., sumx2=0., sumxy=0.;
//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
//     ++nhits;
//     pos = HitPosition(*hit);
//     sumx += pos.X();
//     sumy += pos.Y();
//     sumx2 += pos.X() * pos.X();
//     sumxy += pos.X() * pos.Y();
//   }
//   double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
//   TVector2 direction = TVector2(1,gradient).Unit();

//   // Find how far each hit (projected onto this axis) is from the centre
//   std::map<double,art::Ptr<recob::Hit> > hitProjection;
//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
//     pos = HitPosition(*hit) - centre;
//     hitProjection[direction*pos] = *hit;
//   }

//   // Construct a line from the initial part of the shower

//   // forwards
//   nhits = 0;
//   sumx=0., sumy=0., sumx2=0., sumxy=0.;
//   for (std::map<double,art::Ptr<recob::Hit> >::iterator hitIt = hitProjection.begin(); hitIt != hitProjection.end(); ++hitIt) {
//     if (std::distance(hitProjection.begin(),hitIt) > 5) break;
//     ++nhits;
//     pos = HitPosition(hitIt->second);
//     sumx += pos.X();
//     sumy += pos.Y();
//     sumx2 += pos.X() * pos.X();
//     sumxy += pos.X() * pos.Y();
//   }
//   gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
//   TVector2 directionFromStart = TVector2(1,gradient).Unit();

//   // backwards
//   nhits = 0;
//   sumx=0., sumy=0., sumx2=0., sumxy=0.;
//   for (std::map<double,art::Ptr<recob::Hit> >::reverse_iterator hitIt = hitProjection.rbegin(); hitIt != hitProjection.rend(); ++hitIt) {
//     if (std::distance(hitProjection.rbegin(),hitIt) > 5) break;
//     ++nhits;
//     pos = HitPosition(hitIt->second);
//     sumx += pos.X();
//     sumy += pos.Y();
//     sumx2 += pos.X() * pos.X();
//     sumxy += pos.X() * pos.Y();
//   }
//   gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
//   TVector2 directionFromEnd = TVector2(1,gradient).Unit();

//   // Count how often the distance increases
//   int fromBeginning = 0, fromEnd = 0;
//   double previousDistanceFromAxis = 0;
//   int hitsAfterOffshoot = 0;
//   TVector2 start = HitPosition(hitProjection.begin()->second);
//   for (std::map<double,art::Ptr<recob::Hit> >::iterator hitIt = hitProjection.begin(); hitIt != hitProjection.end(); ++hitIt) {
//     TVector2 pos = HitPosition(hitIt->second);
//     TVector2 projPos = pos.Proj(directionFromStart);
//     double distanceFromAxis = TMath::Sqrt(TMath::Power((projPos-start).Mod(),2) + TMath::Power(pos*directionFromStart,2));
//     if (distanceFromAxis > previousDistanceFromAxis) {
//       ++fromBeginning;
//       hitsAfterOffshoot = 0;
//       previousDistanceFromAxis = distanceFromAxis;
//     }
//     else
//       ++hitsAfterOffshoot;
//     if (hitsAfterOffshoot > 5) {
//       hitsAfterOffshoot = 0;
//       previousDistanceFromAxis = 0;
//     }
//   }

//   previousDistanceFromAxis = 0;
//   hitsAfterOffshoot = 0;
//   TVector2 end = HitPosition(hitProjection.rbegin()->second);
//   for (std::map<double,art::Ptr<recob::Hit> >::reverse_iterator hitIt = hitProjection.rbegin(); hitIt != hitProjection.rend(); ++hitIt) {
//     TVector2 pos = HitPosition(hitIt->second);
//     TVector2 projPos = pos.Proj(directionFromEnd);
//     double distanceFromAxis = TMath::Sqrt(TMath::Power((projPos-start).Mod(),2) + TMath::Power(pos*directionFromStart,2));
//     if (distanceFromAxis > previousDistanceFromAxis) {
//       ++fromEnd;
//       hitsAfterOffshoot = 0;
//       previousDistanceFromAxis = distanceFromAxis;
//     }
//     else
//       ++hitsAfterOffshoot;
//     if (hitsAfterOffshoot > 5) {
//       hitsAfterOffshoot = 0;
//       previousDistanceFromAxis = 0;
//     }
//   }

//   int startWire = fromBeginning > fromEnd ? HitCoordinates(hitProjection.begin()->second).X() : HitCoordinates(hitProjection.rbegin()->second).X();

//   std::vector<art::Ptr<recob::Hit> > trackHits = Construct2DShowerStart(shower, startWire);

//   return trackHits;

// }

// void shower::EMShowerAlg::ProjectVertexIn2D(TVector3 const& vertex,
// 					    std::map<int,std::vector<art::Ptr<recob::Hit> > >& trackHitsMap,
// 					    std::map<int,std::pair<std::vector<art::Ptr<recob::Hit> >,std::vector<art::Ptr<recob::Hit> > > > const& trackHitsBothEndsMap) {

//   /// Projects the 3D direction into all the 2D views to make sure the correct end of the track is selected in each view

//   for (std::map<int,std::pair<std::vector<art::Ptr<recob::Hit> >, std::vector<art::Ptr<recob::Hit> > > >::const_iterator trackHitsIt = trackHitsBothEndsMap.begin();
//        trackHitsIt != trackHitsBothEndsMap.end();
//        ++trackHitsIt) {

//     TVector2 end1 = HitCoordinates(*trackHitsIt->second.first.begin());
//     TVector2 end2 = HitCoordinates(*trackHitsIt->second.second.begin());

//     TVector2 vertexProj = Project3DPointOntoPlane(vertex, trackHitsIt->first);

//     if ( (vertexProj - end1).Mod() < (vertexProj - end2).Mod() )
//       trackHitsMap[trackHitsIt->first] = trackHitsIt->second.first;
//     else
//       trackHitsMap[trackHitsIt->first] = trackHitsIt->second.second;

//   }

// }

// void shower::EMShowerAlg::FindShowerStartDirection(art::Ptr<recob::Track> const& vertexTrack,
// 						   std::map<int,TVector2> const& showerCentreMap,
// 						   TVector3& showerVertex,
// 						   TVector3& showerDirection) {

//   /// Finds the start of the shower, and its direction

//   TVector3 vertex = vertexTrack->Vertex();
//   TVector3 end = vertexTrack->End();

//   std::map<int,double> distanceToVertex;
//   std::map<int,double> distanceToEnd;

//   // Loop over all the planes and find the distance from the vertex and end projections to the centre in each plane
//   for (std::map<int,TVector2>::const_iterator showerCentreIt = showerCentreMap.begin(); showerCentreIt != showerCentreMap.end(); ++showerCentreIt) {

//     // Project the vertex and the end point onto this plane
//     TVector2 vertexProj = Project3DPointOntoPlane(vertex, showerCentreIt->first);
//     TVector2 endProj    = Project3DPointOntoPlane(end, showerCentreIt->first);

//     // Find the distance of each to the centre of the cluster
//     distanceToVertex[showerCentreIt->first] = (vertexProj - showerCentreIt->second).Mod();
//     distanceToEnd[showerCentreIt->first] = (endProj - showerCentreIt->second).Mod();

//   }

//   // Find the average distance to the vertex and the end across the planes
//   double avDistanceToVertex = 0, avDistanceToEnd = 0;
//   for (std::map<int,double>::iterator distanceToVertexIt = distanceToVertex.begin(); distanceToVertexIt != distanceToVertex.end(); ++distanceToVertexIt)
//     avDistanceToVertex += distanceToVertexIt->second;
//   avDistanceToVertex /= distanceToVertex.size();

//   for (std::map<int,double>::iterator distanceToEndIt = distanceToEnd.begin(); distanceToEndIt != distanceToEnd.end(); ++distanceToEndIt)
//     avDistanceToEnd += distanceToEndIt->second;
//   avDistanceToEnd /= distanceToEnd.size();

//   // Set the vertex and directions for this shower
//   if (avDistanceToVertex > avDistanceToEnd) {
//     showerVertex = vertex;
//     showerDirection = vertexTrack->VertexDirection();
//   }
//   else {
//     showerVertex = end;
//     showerDirection = (-1) * vertexTrack->VertexDirection();
//   }

//   return;

// }

// art::Ptr<recob::Track> shower::EMShowerAlg::FindVertexTrack(std::map<int,art::Ptr<recob::Hit> > const& vertexMap,
// 							    std::map<int,art::Ptr<recob::Track> > const& trackMap,
// 							    std::map<int,std::vector<art::Ptr<recob::Hit> > > const& trackHitsMap) {

//   /// Finds the 3D vertex and direction given the tracks associated with the 2D vertex

//   art::Ptr<recob::Track> vertexTrack;

//   // First, find out if two views agree on a track (normally they will)
//   std::map<int,std::vector<int> > trackIDToPlanes;
//   std::vector<int> planesWithTrack;
//   for (std::map<int,art::Ptr<recob::Track> >::const_iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
//     if (!trackIt->second.isNull()) {
//       trackIDToPlanes[trackIt->second->ID()].push_back(trackIt->first);
//       planesWithTrack.push_back(trackIt->first);
//     }
//   for (std::map<int,std::vector<int> >::iterator trackVertexIt = trackIDToPlanes.begin(); trackVertexIt != trackIDToPlanes.end(); ++trackVertexIt)
//     if (trackVertexIt->second.size() > 1)
//       vertexTrack = trackMap.at(*trackVertexIt->second.begin());

//   // If there are no planes with a track reconstructed, don't carry on right now
//   if (planesWithTrack.size() == 0)
//     return vertexTrack;

//   // If they don't, try to use the third view (it exsits) to pick the correct track
//   if (vertexTrack.isNull() and trackMap.size() > 2) {
//     std::map<int,double> distanceOfVertexFromVertex;
//     for (std::map<int,art::Ptr<recob::Hit> >::const_iterator vertexIt = vertexMap.begin(); vertexIt != vertexMap.end(); ++vertexIt) {
//       if (std::find(planesWithTrack.begin(), planesWithTrack.end(), vertexIt->first) != planesWithTrack.end())
// 	continue;
//       for (std::vector<int>::iterator planesWithTrackIt = planesWithTrack.begin(); planesWithTrackIt != planesWithTrack.end(); ++planesWithTrackIt)
// 	distanceOfVertexFromVertex[*planesWithTrackIt] = ( TMath::Abs( HitCoordinates(vertexIt->second).Y() - HitCoordinates(vertexMap.at(*planesWithTrackIt)).Y() ) ) / HitCoordinates(vertexMap.at(*planesWithTrackIt)).Y();
//     }
//     std::vector<int> thirdViewClose;
//     for (std::map<int,double>::iterator distanceOfVertexFromVertexIt = distanceOfVertexFromVertex.begin(); distanceOfVertexFromVertexIt != distanceOfVertexFromVertex.end(); ++distanceOfVertexFromVertexIt)
//       if (distanceOfVertexFromVertexIt->second < 0.1)
// 	thirdViewClose.push_back(distanceOfVertexFromVertexIt->first);
//     if (thirdViewClose.size() == 1)
//       vertexTrack = trackMap.at(thirdViewClose.at(0));
//   }

//   // Finally, if all else fails, just pick the view with the longest initial track reconstructed in 2D
//   if (vertexTrack.isNull()) {
//     std::map<int,int> lengthOfTrackToPlane;
//     for (std::vector<int>::iterator planesWithTrackIt = planesWithTrack.begin(); planesWithTrackIt != planesWithTrack.end(); ++planesWithTrackIt)
//       lengthOfTrackToPlane[trackHitsMap.at(*planesWithTrackIt).size()] = *planesWithTrackIt;
//     vertexTrack = trackMap.at(lengthOfTrackToPlane.rbegin()->second);
//   }

//   return vertexTrack;

// }

// std::vector<art::Ptr<recob::Hit> > shower::EMShowerAlg::FindTrack(std::vector<art::Ptr<recob::Hit> > const& shower, TVector2 const& start, TVector2 const& end) {

//   /// Finds the track from the start of the shower

//   std::vector<art::Ptr<recob::Hit> > trackHits;

//   // Map of hit on each wire
//   std::map<int,std::vector<art::Ptr<recob::Hit> > > hitWires;
//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
//     hitWires[(int)HitCoordinates(*hit).X()].push_back(*hit);

//   // Call the track ended when there is more than one hit on a wire twice in three wires
//   std::vector<int> lastTwoWires = {0, 0};

//   // Find out which way to look!
//   int startWire = (int)start.X(), endWire = (int)end.X();
//   bool increasing;
//   if (startWire < endWire) increasing = true;
//   else increasing = false;

//   // Look through the hits from the start
//   if (increasing) {
//     for (int wire = startWire; wire < endWire; ++wire) {
//       int numberOfHitsOnThisWire;
//       if (hitWires.find(wire) != hitWires.end()) numberOfHitsOnThisWire = hitWires.at(wire).size();
//       else numberOfHitsOnThisWire = 0;
//       if (numberOfHitsOnThisWire >= 2 and (lastTwoWires.at(0) >= 2 or lastTwoWires.at(1) >= 2))
// 	break;
//       else {
// 	if (numberOfHitsOnThisWire)
// 	  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitWireIt = hitWires.at(wire).begin(); hitWireIt != hitWires.at(wire).end(); ++hitWireIt)
// 	    trackHits.push_back(*hitWireIt);
// 	lastTwoWires[0] = lastTwoWires[1];
// 	lastTwoWires[1] = numberOfHitsOnThisWire;
//       }
//     }
//   }

//   else {
//     for (int wire = startWire; wire > endWire; --wire) {
//       int numberOfHitsOnThisWire;
//       if (hitWires.find(wire) != hitWires.end()) numberOfHitsOnThisWire = hitWires.at(wire).size();
//       else numberOfHitsOnThisWire = 0;
//       if (numberOfHitsOnThisWire >= 2 and (lastTwoWires.at(0) >= 2 or lastTwoWires.at(1) >= 2))
// 	break;
//       else {
// 	if (numberOfHitsOnThisWire)
// 	  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitWireIt = hitWires.at(wire).begin(); hitWireIt != hitWires.at(wire).end(); ++hitWireIt)
// 	    trackHits.push_back(*hitWireIt);
// 	lastTwoWires[0] = lastTwoWires[1];
// 	lastTwoWires[1] = numberOfHitsOnThisWire;
//       }
//     }
//   }

//   if (trackHits.size() == 0)
//     trackHits.push_back(hitWires.begin()->second.at(0));

//   return trackHits;

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

  //std::cout << "Looking at wire " << wireID.Wire << ", in TPC " << wireID.TPC << " and on plane " << wireID.Plane << std::endl;

  double globalWire = -999;
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }
  else {
    //std::cout << "It's collection!" << std::endl;
    unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
    //std::cout << "Number of wires is " << nwires << std::endl;
    if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
    else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
    else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
    else mf::LogError("EMShowerAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC;
    //std::cout << "Global wire is " << globalWire << std::endl;
  }

  return globalWire;

}

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane(TVector3 const& point, unsigned int plane) {

  /// Projects a 3D point onto a 2D plane

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  return TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), plane, fGeom->FindTPCAtPosition(pointPosition).TPC % 2, 0),
		  fDetProp->ConvertXToTicks(point.X(), plane, fGeom->FindTPCAtPosition(pointPosition).TPC % 2, 0));

}

recob::Shower shower::EMShowerAlg::MakeShower(art::PtrVector<recob::Hit> const& hits,
					      recob::Track const& initialTrack,
					      std::vector<art::Ptr<recob::Hit> > const& initialTrackHits) {

  /// Makes a recob::Shower object given the hits in the shower and the initial track-like part

  return recob::Shower();

  // Find the shower hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  // Find the initial track hits on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > initialHitsMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator initialHitIt = initialTrackHits.begin(); initialHitIt != initialTrackHits.end(); ++initialHitIt)
    initialHitsMap[(*initialHitIt)->View()].push_back(*initialHitIt);

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

  TVector3 direction = initialTrack.VertexDirection(), directionError;
  TVector3 showerStart = initialTrack.Vertex(), showerStartError;

  std::cout << "Best plane is " << bestPlane << std::endl;
  std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
  std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;
  std::cout << "The shower start is " << std::endl;
  showerStart.Print();

  return recob::Shower(direction, directionError, showerStart, showerStartError, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);

}

double shower::EMShowerAlg::FinddEdx(std::vector<art::Ptr<recob::Hit> > const& trackHits, recob::Track const& track) {

  /// Finds dE/dx for the track given a set of hits

  double totalCharge = 0, totalDistance = 0, avHitTime = 0;
  unsigned int nHits = 0;

  // Get the pitch
  double pitch = 0;
  try { pitch = track.PitchInView(trackHits.at(0)->View()); }
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

recob::Track shower::EMShowerAlg::ConstructTrack(std::vector<art::Ptr<recob::Hit> > const& track1,
						 std::vector<art::Ptr<recob::Hit> > const& track2) {

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).

  pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(track1, track2);

  std::vector<TVector3> xyz(pmatrack->size()), dircos(pmatrack->size());
  std::vector<std::vector<double> > dst_dQdx(3, std::vector<double>());

  unsigned int cryo = pmatrack->FrontCryo();
  unsigned int tpc = pmatrack->FrontTPC();

  typedef std::map<size_t,std::vector<double> > dedx_map;
  std::map<unsigned int,dedx_map> pmatrack_dQdx;
  if (fGeom->TPC(tpc, cryo).HasPlane(geo::kU)) {
    pmatrack_dQdx[geo::kU] = dedx_map();
    pmatrack->GetRawdEdxSequence(pmatrack_dQdx[geo::kU], geo::kU);
  }
  if (fGeom->TPC(tpc, cryo).HasPlane(geo::kV)) {
    pmatrack_dQdx[geo::kV] = dedx_map();
    pmatrack->GetRawdEdxSequence(pmatrack_dQdx[geo::kV], geo::kV);
  }
  if (fGeom->TPC(tpc, cryo).HasPlane(geo::kZ)) {
    pmatrack_dQdx[geo::kZ] = dedx_map();
    pmatrack->GetRawdEdxSequence(pmatrack_dQdx[geo::kZ], geo::kZ);
  }
  
  TVector3 p3d;
  double xshift = pmatrack->GetXShift();
  bool has_shift = (xshift != 0.0);
  for (size_t i = 0; i < pmatrack->size(); i++)
    if ((*pmatrack)[i]->IsEnabled()) {
      p3d = (*pmatrack)[i]->Point3D();
      if (has_shift) p3d.SetX(p3d.X() + xshift);
      xyz.push_back(p3d);
      
      if (i < pmatrack->size() - 1) {
	TVector3 dc((*pmatrack)[i+1]->Point3D());
	dc -= (*pmatrack)[i]->Point3D();
	dc *= 1.0 / dc.Mag();
	dircos.push_back(dc);
      }
      else dircos.push_back(dircos.back());
	
      double dQ = 0., dx = 0.;
      dst_dQdx[geo::kU].push_back(0.);
      dst_dQdx[geo::kV].push_back(0.);
      dst_dQdx[geo::kZ].push_back(0.);
	
      double dQdx;
      for (auto const& m : pmatrack_dQdx) {
	auto it = m.second.find(i);
	if (it != m.second.end()) {
	  dQ = it->second[5];
	  dx = it->second[6];
	  if (dx > 0.) dQdx = dQ/dx;
	  else dQdx = 0.;
	    
	  size_t backIdx = dst_dQdx[m.first].size() - 1;
	  dst_dQdx[m.first][backIdx] = dQdx;
	  
	  break;
	}
      }
    }

  if (xyz.size() != dircos.size())
    mf::LogError("EMShowerAlg") << "Problem converting pma::Track3D to recob::Track";
  
  return recob::Track(xyz, dircos, dst_dQdx);

}
