////////////////////////////////////////////////////////////////////////////
// Class:       EMShower
// Module Type: producer
// File:        EMShower_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
//
////////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"
#include "MCCheater/BackTracker.h"

// LArSoft includes
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Shower.h"

// ROOT includes
#include "TPrincipal.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TF2.h"
#include "TTree.h"

//#include "RecoAlg/EMShowerAlg.hxx"

namespace shower {
  class EMShower;
}

class shower::EMShower : public art::EDProducer {
public:

  EMShower(fhicl::ParameterSet const& pset);

  void produce(art::Event &evt);
  void reconfigure(fhicl::ParameterSet const& p);
  void MakeShowers(std::map<int,std::vector<int> > const& trackToClusters, std::vector<std::vector<int> >& showers);
  void FindInitialTrack(art::PtrVector<recob::Hit> const& hits);
  void FindShowerEnds(art::PtrVector<recob::Hit> const& shower, TVector2 const& centre, TVector2& end1, TVector2& end2);
  void FindVertex(art::PtrVector<recob::Hit> const& shower, TVector2 const& end1, TVector2 const& end2, std::vector<int>& trackHits);
  void FindTrack(art::PtrVector<recob::Hit> const& shower, std::map<double,int> const& hitToEnd, std::vector<int>& trackHits);
  void FindTrack(TVector2 const& start, TVector2 const& end, std::map<int,std::vector<int> > const& hitWires, std::vector<int>& trackHits);
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);
  double GlobalWire(geo::WireID wireID);

private:

  std::string fHitsModuleLabel, fClusterModuleLabel, fTrackModuleLabel;

  //EMShowerAlg fEMShowerAlg;

  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;


};

shower::EMShower::EMShower(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::Track> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
}

void shower::EMShower::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel    = p.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel = p.get<std::string>("ClusterModuleLabel");
  fTrackModuleLabel   = p.get<std::string>("TrackModuleLabel");
}

void shower::EMShower::produce(art::Event &evt) {

  // Output -- showers and associations with hits and clusters
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitAssociations(new art::Assns<recob::Shower, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Track> > trackAssociations(new art::Assns<recob::Shower, recob::Track>);
  std::unique_ptr<art::Assns<recob::Shower, recob::SpacePoint> > spacePointAssociations(new art::Assns<recob::Shower, recob::SpacePoint>);

  // Event has hits, tracks and clusters found already

  // Hits
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Tracks
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel(fTrackModuleLabel,trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);

  // Clusters
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > clusters;
  if (evt.getByLabel(fClusterModuleLabel,clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Track> fmt(hitHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmsp(trackHandle, evt, fTrackModuleLabel);

  // Map between tracks and clusters
  std::map<int,std::vector<int> > trackToClusters;
  std::map<int,std::vector<int> > clusterToTracks;

  // Look through all the clusters
  for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt) {

    // Get the hits in this cluster
    std::vector<art::Ptr<recob::Hit> > clusterHits = fmh.at(clusterIt->key());

    // Look at all these hits and find the associated tracks
    for (std::vector<art::Ptr<recob::Hit> >::iterator clusterHitIt = clusterHits.begin(); clusterHitIt != clusterHits.end(); ++clusterHitIt) {

      // Get the tracks associated with this hit
      std::vector<art::Ptr<recob::Track> > clusterHitTracks = fmt.at(clusterHitIt->key());
      if (clusterHitTracks.size() > 1) { std::cout << "More than one track associated with this hit!" << std::endl; continue; }
      if (clusterHitTracks.size() < 1) continue;

      // Add this cluster to the track map
      int track = clusterHitTracks.at(0).key();
      int cluster = (*clusterIt).key();
      if (std::find(trackToClusters[track].begin(), trackToClusters[track].end(), cluster) == trackToClusters[track].end())
	trackToClusters[track].push_back(cluster);
      if (std::find(clusterToTracks[cluster].begin(), clusterToTracks[cluster].end(), track) == clusterToTracks[cluster].end())
	clusterToTracks[cluster].push_back(track);

    }     

  }

  // Make showers
  std::vector<std::vector<int> > newShowers;
  MakeShowers(trackToClusters, newShowers);

  // Make output larsoft products -----------------------

  int showerNum = 0;
  for (std::vector<std::vector<int> >::iterator newShower = newShowers.begin(); newShower != newShowers.end(); ++newShower, ++showerNum) {

    // New shower
    std::cout << std::endl << "Start shower " << showerNum << std::endl;

    // New associations
    art::PtrVector<recob::Hit> showerHits;
    art::PtrVector<recob::Cluster> showerClusters;
    art::PtrVector<recob::Track> showerTracks;
    art::PtrVector<recob::SpacePoint> showerSpacePoints;

    std::vector<int> associatedTracks;

    // Make showers and associations
    for (std::vector<int>::iterator showerCluster = (*newShower).begin(); showerCluster != (*newShower).end(); ++showerCluster) {

      // Clusters
      art::Ptr<recob::Cluster> cluster = clusters.at(*showerCluster);
      showerClusters.push_back(cluster);

      // Hits
      std::vector<art::Ptr<recob::Hit> > showerClusterHits = fmh.at(cluster.key());
      for (std::vector<art::Ptr<recob::Hit> >::iterator showerClusterHit = showerClusterHits.begin(); showerClusterHit != showerClusterHits.end(); ++showerClusterHit)
	showerHits.push_back(*showerClusterHit);

      // Tracks
      std::vector<int> clusterTracks = clusterToTracks.at(*showerCluster);
      for (std::vector<int>::iterator clusterTracksIt = clusterTracks.begin(); clusterTracksIt != clusterTracks.end(); ++clusterTracksIt)
	if (std::find(associatedTracks.begin(), associatedTracks.end(), *clusterTracksIt) == associatedTracks.end())
	  associatedTracks.push_back(*clusterTracksIt);

    }

    // Tracks and space points
    for (std::vector<int>::iterator associatedTracksIt = associatedTracks.begin(); associatedTracksIt != associatedTracks.end(); ++associatedTracksIt) {
      art::Ptr<recob::Track> showerTrack = tracks.at(*associatedTracksIt);
      showerTracks.push_back(showerTrack);
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmsp.at(showerTrack.key());
      for (std::vector<art::Ptr<recob::SpacePoint> >::iterator spacePointsIt = spacePoints.begin(); spacePointsIt != spacePoints.end(); ++spacePointsIt)
	showerSpacePoints.push_back(*spacePointsIt);
    }

    // Find the initial track-like part of the shower
    FindInitialTrack(showerHits);

    // Make shower object and associations
    recob::Shower shower = recob::Shower();
    shower.set_id(showerNum);
    showers->push_back(shower);

    util::CreateAssn(*this, evt, *(showers.get()), showerHits,        *(hitAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerClusters,    *(clusterAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerTracks,      *(trackAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerSpacePoints, *(spacePointAssociations.get()));

  }

  // Put in event
  evt.put(std::move(showers));
  evt.put(std::move(hitAssociations));
  evt.put(std::move(clusterAssociations));
  evt.put(std::move(trackAssociations));
  evt.put(std::move(spacePointAssociations));

}

void shower::EMShower::MakeShowers(std::map<int,std::vector<int> > const& trackToClusters, std::vector<std::vector<int> >& showers) {

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

void shower::EMShower::FindInitialTrack(art::PtrVector<recob::Hit> const& hits) {

  /// Finds the initial track-like part of the shower given all the hits in it

  // Consider each plane separately
  std::map<int,art::PtrVector<recob::Hit> > planeHitsMap;
  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit)
    planeHitsMap[(*hit)->View()].push_back(*hit);

  std::map<int,TVector2> vertexMap;
  std::map<int,art::Ptr<recob::Hit> > vertexHitMap;

  // Loop through planes
  for (std::map<int,art::PtrVector<recob::Hit> >::iterator planeHits = planeHitsMap.begin(); planeHits != planeHitsMap.end(); ++planeHits) {

    std::cout << "Plane " << planeHits->first << std::endl;

    // Find the charge-weighted centre of this shower
    TVector2 pos, chargePoint = TVector2(0,0);
    double totalCharge = 0;
    for (std::vector<art::Ptr<recob::Hit> >::iterator planeHit = planeHits->second.begin(); planeHit != planeHits->second.end(); ++planeHit) {
      pos = HitCoordinates(*planeHit);
      chargePoint += (*planeHit)->Integral() * pos;
      totalCharge += (*planeHit)->Integral();
    }
    TVector2 centre = chargePoint / totalCharge;

    // Find rough 'ends' of the shower!
    TVector2 end1, end2;
    FindShowerEnds(planeHits->second, centre, end1, end2);

    // Decide which end is the vertex and get the initial track
    std::vector<int> trackHits;
    FindVertex(planeHits->second, end1, end2, trackHits);
    // TVector2 vertex = HitCoordinates(planeHits->second.at(*trackHits.begin()));
    // art::Ptr<recob::Hit> vertexHit = planeHits->second.at(*trackHits.begin());
    // TVector2 trackEnd = HitCoordinates(planeHits->second.at(*trackHits.rbegin()));

    // vertexMap[planeHits->first] = vertex;
    // vertexHitMap[planeHits->first] = vertexHit;

    // std::cout << "Vertex in this plane is " << std::endl;
    // vertex.Print();
    // std::cout << "End of track in this plane is " << std::endl;
    // trackEnd.Print();

  }

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

void shower::EMShower::FindShowerEnds(art::PtrVector<recob::Hit> const& shower, TVector2 const& centre, TVector2& end1, TVector2& end2) {

  /// Roughly finds the two 'ends' of the shower (one is the vertex, one isn't well defined!)

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

  end1 = HitCoordinates(shower.at(oneHit));
  end2 = HitCoordinates(shower.at(otherHit));

  return;

}

void shower::EMShower::FindVertex(art::PtrVector<recob::Hit> const& shower, TVector2 const& end1, TVector2 const& end2, std::vector<int>& trackHits) {

  /// Decides which 'end' of the shower is the true vertex and finds the initial track

  // Map of hit on each wire
  std::map<int,std::vector<int> > hitWires;
  for (art::PtrVector<recob::Hit>::const_iterator hit = shower.begin(); hit != shower.end(); ++hit)
    hitWires[(int)HitCoordinates(*hit).X()].push_back(std::distance(shower.begin(), hit));

  std::vector<int> trackHits1, trackHits2;
  FindTrack(end1, end2, hitWires, trackHits1);
  FindTrack(end2, end1, hitWires, trackHits2);

  trackHits = trackHits1.size() > trackHits2.size() ? trackHits1 : trackHits2;

  std::cout << "Vertex is (" << HitCoordinates(shower.at(*trackHits.begin())).X() << ", " << HitCoordinates(shower.at(*trackHits.begin())).Y() << ") and the end of the track is (" << HitCoordinates(shower.at(*trackHits.rbegin())).X() << ", " << HitCoordinates(shower.at(*trackHits.rbegin())).Y() << ")" << std::endl;

  return;

}

void shower::EMShower::FindTrack(TVector2 const& start, TVector2 const& end, std::map<int,std::vector<int> > const& hitWires, std::vector<int>& trackHits) {

  /// Finds the track from the start of the shower

  // Call the track ended when there is more than one hit on a wire twice in three wires
  std::vector<int> lastTwoWires = {0, 0};

  // Find out which way too look!
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

}

void shower::EMShower::FindTrack(art::PtrVector<recob::Hit> const& shower, std::map<double,int> const& hitToEnd, std::vector<int>& trackHits) {

  /// Contructs a track from the 'end' of a shower

  TVector2 end = HitCoordinates(shower.at(hitToEnd.begin()->second));
  std::cout << "End coordinates are " << std::endl;
  end.Print();

  // Look through the hits from the end
  for (std::map<double,int>::const_iterator hitToEndIt = hitToEnd.begin(); hitToEndIt != hitToEnd.end(); ++hitToEndIt) {

    // Put first hit in vector
    if (hitToEnd.begin() == hitToEndIt) {
      trackHits.push_back(hitToEndIt->second);
      continue;
    }

    // Find the direction from all hits up to this one
    double nhits, sumx, sumy, sumx2, sumxy;
    for (std::map<double,int>::const_iterator hitsInTrackIt = hitToEnd.begin(); hitsInTrackIt != std::next(hitToEndIt,1); ++hitsInTrackIt) {
      ++nhits;
      TVector2 hitpos = HitCoordinates(shower.at(hitsInTrackIt->second));
      sumx += hitpos.X();
      sumy += hitpos.Y();
      sumx2 += hitpos.X() * hitpos.X();
      sumxy += hitpos.X() * hitpos.Y();
    }
    double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
    TVector2 direction = TVector2(1,gradient).Unit();

    // Project all the hits onto this direction and find average displacement
    double totalDisplacement = 0;
    for (std::map<double,int>::const_iterator hitsInTrackIt = hitToEnd.begin(); hitsInTrackIt != std::next(hitToEndIt,1); ++hitsInTrackIt) {
      TVector2 hitpos = HitCoordinates(shower.at(hitsInTrackIt->second)) - end;
      totalDisplacement += TMath::Sqrt(TMath::Power(hitpos.Mod(),2) - TMath::Power((hitpos.Proj(direction)).Mod(),2));
    }
    double averageDisplacement = totalDisplacement / (std::distance(hitToEnd.begin(), hitToEndIt) + 1);

    std::cout << "Hit at (" << HitCoordinates(shower.at(hitToEndIt->second)).X() << ", " << HitCoordinates(shower.at(hitToEndIt->second)).Y() << ") [distance of " << hitToEndIt->first << "] leads to av displacement of " << averageDisplacement << std::endl;

    // if (averageDisplacement < 1)
    //   trackHits.push_back(hitToEndIt->second);
    // else
    //   break;

  }

  return;

}

// void shower::EMShower::FindVertex(art::PtrVector<recob::Hit> const& shower, TVector2 const& end1, TVector2 const& end2, std::vector<int>& trackHits) {

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

// void shower::EMShower::FindTrack(art::PtrVector<recob::Hit> const& shower, std::map<double,int> const& hitToEnd, std::vector<int>& trackHits) {

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

TVector2 shower::EMShower::HitCoordinates(art::Ptr<recob::Hit> const& hit) {

  /// Return the coordinates of this hit in global wire/tick space

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());

}

double shower::EMShower::GlobalWire(geo::WireID wireID) {

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

DEFINE_ART_MODULE(shower::EMShower)
