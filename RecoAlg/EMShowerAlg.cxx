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

void shower::EMShowerAlg::MakeShowers(std::map<int,std::vector<int> > const& trackToClusters, std::vector<std::vector<int> >& showers) {

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

void shower::EMShowerAlg::FindShowerProperties(art::PtrVector<recob::Hit> const& hits, art::FindManyP<recob::Track> const& fmt, calo::CalorimetryAlg const& calo,
					       TVector3& direction, TVector3& directionError, TVector3& vertex, TVector3& vertexError,
					       std::vector<double>& totalEnergy, std::vector<double>& totalEnergyError, std::vector<double>& dEdx, std::vector<double>& dEdxError,
					       int& bestPlane) {

  /// Finds the properties of the shower from the hits in it

  // Consider each plane separately
  std::map<int,art::PtrVector<recob::Hit> > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  std::map<int,art::Ptr<recob::Hit> > vertexMap;
  std::map<int,art::Ptr<recob::Track> > trackMap;
  std::map<int,std::vector<int> > trackHitsMap;

  unsigned int highestNumberOfHits = 0;

  // Loop through planes to find vertices in each view
  for (std::map<int,art::PtrVector<recob::Hit> >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {
    
    std::cout << "Plane " << planeHits->first << std::endl;

    // Find best plane
    if (planeHits->second.size() > highestNumberOfHits) {
      highestNumberOfHits = planeHits->second.size();
      bestPlane = planeHits->first;
    }

    // Find rough 'ends' of the shower!
    art::Ptr<recob::Hit> end1, end2;
    FindShowerEnds(planeHits->second, end1, end2);

    // // Decide which end is the vertex
    // art::Ptr<recob::Hit> vertex = FindVertex(planeHits->second, end1, end2);
    // art::Ptr<recob::Hit> showerEnd = vertex.key() == end1.key() ? end2 : end1;
    // vertexMap.push_back(vertex);

    // // Get the initial track
    // std::vector<int> trackHits = FindTrack(planeHits->second, HitCoordinates(vertex), HitCoordinates(showerEnd));
    // art::Ptr<recob::Track> track = fmt.at(hits.at(trackHits.at(0)).key()).at(0);

    // Get the initial track
    std::vector<int> trackHits1 = FindTrack(planeHits->second, HitCoordinates(end1), HitCoordinates(end2));
    std::vector<int> trackHits2 = FindTrack(planeHits->second, HitCoordinates(end2), HitCoordinates(end1));
    std::vector<int> trackHits = trackHits1.size() > trackHits2.size() ? trackHits1 : trackHits2;
    trackHitsMap[planeHits->first] = trackHits;

    // Find a track object
    art::Ptr<recob::Track> track;
    for (std::vector<int>::iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
      std::vector<art::Ptr<recob::Track> > tracks = fmt.at(planeHits->second.at(*trackHitIt).key());
      if (tracks.size() == 0)
	continue;
      track = tracks.at(0);
      break;
    }
    trackMap[planeHits->first] = track;

    // Vertex
    art::Ptr<recob::Hit> planevertex = planeHits->second.at(*trackHits.begin());
    vertexMap[planeHits->first] = planevertex;
    std::cout << "The vertex is " << std::endl;
    HitCoordinates(planevertex).Print();

  }

  // Find 3D vertex and direction
  art::Ptr<recob::Track> vertexTrack;
  FindVertexTrack(vertexTrack, vertexMap, trackMap, trackHitsMap);
  vertex = vertexTrack->Vertex();
  direction = vertexTrack->VertexDirection();

  // Find energy and dE/dx
  for (unsigned int plane = 0; plane < fGeom->MaxPlanes(); ++plane) {
    if (planeHitsMap.count(plane) != 0) {
      dEdx.push_back(FinddEdx(planeHitsMap.at(plane), vertexTrack, calo, vertexMap.at(plane)->View(), trackHitsMap.at(plane)));
      totalEnergy.push_back(FindTotalEnergy(planeHitsMap.at(plane), plane));
    }
    else {
      dEdx.push_back(0);
      totalEnergy.push_back(0);
    }
  }

  std::cout << "Best plane is " << bestPlane << std::endl;
  std::cout << "dE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2] << std::endl;
  std::cout << "Total energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1] << " and " << totalEnergy[2] << std::endl;

  return;

}

void shower::EMShowerAlg::FindVertexTrack(art::Ptr<recob::Track>& vertexTrack, std::map<int,art::Ptr<recob::Hit> > const& vertexMap, std::map<int,art::Ptr<recob::Track> > const& trackMap, std::map<int,std::vector<int> > const& trackHitsMap) {

  /// Finds the 3D vertex and direction given the tracks associated with the 2D vertex

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
    return;

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

  return;

}

double shower::EMShowerAlg::FindTotalEnergy(art::PtrVector<recob::Hit> const& hits, int plane) {

  /// Finds the total energy deposited by the shower in this view

  double totalCharge = 0, totalEnergy = 0;

  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    totalCharge += (*hit)->Integral();

  double Uintercept = -1519.33, Ugradient = 148867;
  double Vintercept = -1234.91, Vgradient = 149458;
  double Zintercept = -1089.73, Zgradient = 145372;

  switch (plane) {
  case 0:
    totalEnergy = (double)(totalCharge - Uintercept)/(double)Ugradient;
    break;
  case 1:
    totalEnergy = (double)(totalCharge - Vintercept)/(double)Vgradient;
    break;
  case 2:
    totalEnergy = (double)(totalCharge - Zintercept)/(double)Zgradient;
    break;
  }

  return totalEnergy;

}

double shower::EMShowerAlg::FinddEdx(art::PtrVector<recob::Hit> const& shower, art::Ptr<recob::Track> const& track, calo::CalorimetryAlg const& calo, geo::View_t const& view, std::vector<int> const& trackHits) {

  /// Finds dE/dx for the track given a set of hits

  std::vector<double> dEdx;

  for (std::vector<int>::const_iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt)
    dEdx.push_back(calo.dEdx_AREA(shower.at(*trackHitIt), track->PitchInView(view)));

  double avdEdx = 0;
  for (std::vector<double>::iterator dEdxIt = dEdx.begin(); dEdxIt != dEdx.end(); ++dEdxIt)
    avdEdx += *dEdxIt;

  avdEdx /= dEdx.size();

  return avdEdx;

}

std::vector<int> shower::EMShowerAlg::FindTrack(art::PtrVector<recob::Hit> const& shower, TVector2 const& start, TVector2 const& end) {

  /// Finds the track from the start of the shower

  std::vector<int> trackHits;

  // Map of hit on each wire
  std::map<int,std::vector<int> > hitWires;
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
    hitWires[(int)HitCoordinates(*hit).X()].push_back(std::distance(shower.begin(), hit));

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
	  for (std::vector<int>::const_iterator hitWireIt = hitWires.at(wire).begin(); hitWireIt != hitWires.at(wire).end(); ++hitWireIt)
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
	  for (std::vector<int>::const_iterator hitWireIt = hitWires.at(wire).begin(); hitWireIt != hitWires.at(wire).end(); ++hitWireIt)
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

void shower::EMShowerAlg::FindShowerEnds(art::PtrVector<recob::Hit> const& shower, art::Ptr<recob::Hit>& end1, art::Ptr<recob::Hit>& end2) {

  /// Roughly finds the two 'ends' of the shower (one is the vertex, one isn't well defined!)

  // Find the charge-weighted centre of this shower
  TVector2 pos, chargePoint = TVector2(0,0);
  double totalCharge = 0;
  for (art::PtrVector<recob::Hit>::const_iterator planeHit = shower.begin(); planeHit != shower.end(); ++planeHit) {
    pos = HitCoordinates(*planeHit);
    chargePoint += (*planeHit)->Integral() * pos;
    totalCharge += (*planeHit)->Integral();
  }
  TVector2 centre = chargePoint / totalCharge;

  // First need to find each end of the shower
  std::map<double,int> hitDistances;
  std::map<int,TVector2> hitCentreVector;

  // Find the distance of each hit from the shower centre
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    TVector2 pos = HitCoordinates(*hit);
    double distanceFromCentre = (pos - centre).Mod();
    hitDistances[distanceFromCentre] = std::distance(shower.begin(),hit);
    hitCentreVector[std::distance(shower.begin(),hit)] = pos - centre;
  }

  // Use this to find the end points
  int oneHit = 0, otherHit = 0;
  for (std::map<double,int>::reverse_iterator hitDistance = hitDistances.rbegin(); hitDistance != hitDistances.rend(); ++hitDistance) {
    if (hitDistance == hitDistances.rbegin())
      oneHit = hitDistance->second;
    else
      if (TMath::Abs(hitCentreVector.at(oneHit).DeltaPhi(hitCentreVector.at(hitDistance->second))) > TMath::Pi()/2) {
	otherHit = hitDistance->second;
	break;
      }
  }

  end1 = shower.at(oneHit);
  end2 = shower.at(otherHit);

  return;

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
