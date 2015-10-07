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

void shower::EMShowerAlg::FindShowerProperties(art::PtrVector<recob::Hit> const& hits, art::FindManyP<recob::Track> const& fmt,
					       TVector3& direction, TVector3& directionError, TVector3& vertex, TVector3& vertexError,
					       std::vector<double>& totalEnergy, std::vector<double>& totalEnergyError, std::vector<double>& dEdx, std::vector<double>& dEdxError,
					       int bestPlane) {

  /// Finds the properties of the shower from the hits in it

  // Consider each plane separately
  std::map<int,art::PtrVector<recob::Hit> > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  std::map<int,TVector2> vertexMap;
  std::map<int,art::Ptr<recob::Hit> > vertexHitMap;

  // Loop through planes
  for (std::map<int,art::PtrVector<recob::Hit> >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::cout << "Plane " << planeHits->first << std::endl;

    // Find rough 'ends' of the shower!
    art::Ptr<recob::Hit> end1, end2;
    FindShowerEnds(planeHits->second, end1, end2);

    // Decide which end is the vertex and get the initial track
    art::Ptr<recob::Hit> vertex = FindVertex(planeHits->second, end1, end2);
    art::Ptr<recob::Hit> showerEnd = vertex.key() == end1.key() ? end2 : end1;
    std::vector<int> track = FindTrack(planeHits->second, HitCoordinates(vertex), HitCoordinates(showerEnd));

    std::cout << "The vertex is " << std::endl;
    HitCoordinates(vertex).Print();
    std::cout << "and the other end of the shower is " << std::endl;
    HitCoordinates(showerEnd).Print();

  }

  //art::Ptr<recob::Track> initialTrack = fmt.at(vertex.key()).at(0);
  if (vertexMap.size() < 2) return;

  // std::cout << "Vertex x is " << fDetProp->ConvertTicksToX(vertexHitMap.at(0)->PeakTime(), vertexHitMap.at(0)->WireID().planeID()) << ", " << fDetProp->ConvertTicksToX(vertexHitMap.at(1)->PeakTime(), vertexHitMap.at(1)->WireID().planeID()) << " and " << fDetProp->ConvertTicksToX(vertexHitMap.at(2)->PeakTime(), vertexHitMap.at(2)->WireID().planeID()) << std::endl;

  // geo::WireIDIntersection yz;
  // fGeom->WireIDsIntersect(vertexHitMap.at(0)->WireID(), vertexHitMap.at(1)->WireID(), yz);
  // std::cout << "y and z: From 0 and 1... " << yz.y << " and " << yz.z << std::endl;
  // fGeom->WireIDsIntersect(vertexHitMap.at(0)->WireID(), vertexHitMap.at(2)->WireID(), yz);
  // std::cout << "y and z: From 0 and 2... " << yz.y << " and " << yz.z << std::endl;
  // fGeom->WireIDsIntersect(vertexHitMap.at(1)->WireID(), vertexHitMap.at(2)->WireID(), yz);
  // std::cout << "y and z: From 1 and 2... " << yz.y << " and " << yz.z << std::endl;

  return;

}

art::Ptr<recob::Hit> shower::EMShowerAlg::FindVertex(art::PtrVector<recob::Hit> const& shower, art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2) {

  /// Decides which 'end' of the shower is the true vertex

  TVector2 end1 = HitCoordinates(hit1);
  TVector2 end2 = HitCoordinates(hit2);

  // Find the geometrical centre of the shower and the direction vector
  TVector2 centre = (end1 + end2) / 2;
  TVector2 direction = (end1 - end2).Unit();

  // Sum all deposited charge for each end of the shower
  double chargeEnd1 = 0, chargeEnd2 = 0;

  // Project all hits onto this vector to determine which end each is closer to
  TVector2 proj;
  double distanceEnd1, distanceEnd2;
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit) {
    proj = (HitCoordinates(*hit)-centre).Proj(direction);
    distanceEnd1 = (end1 - centre - proj).Mod();
    distanceEnd2 = (end2 - centre - proj).Mod();
    if (distanceEnd1 < distanceEnd2)
      chargeEnd1 += (*hit)->Integral();
    else
      chargeEnd2 += (*hit)->Integral();
  }

  std::cout << "Charge nearest end1 is " << chargeEnd1 << " and end2 is " << chargeEnd2 << std::endl;

  // The half of the shower nearest the vertex will deposit less energy
  return chargeEnd1 < chargeEnd2 ? hit1 : hit2;

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
  int oneHit, otherHit;
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
