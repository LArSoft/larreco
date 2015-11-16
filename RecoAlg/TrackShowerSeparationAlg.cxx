////////////////////////////////////////////////////////////////////////
// Class: TrackShowerSeparationAlg
// File:  TrackShowerSeparationAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/TrackShowerSeparationAlg.h"

shower::TrackShowerSeparationAlg::TrackShowerSeparationAlg(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  // tree
  ftree = tfs->make<TTree>("tree","tree");
  ftree->Branch("Run",&Run);
  ftree->Branch("Event",&Event);
  ftree->Branch("Distance",&Distance);
  ftree->Branch("Angle",&Angle);
  ftree->Branch("Length",&Length);
  ftree->Branch("TrackID",&TrackID);
  ftree->Branch("PDG",&pdg);
  ftree->Branch("NSpacePoints",&NSpacePoints);
  ftree->Branch("AvDistance",&AvDistance);
}

void shower::TrackShowerSeparationAlg::reconfigure(fhicl::ParameterSet const& pset) {
  fAngleCut = 5;
  fDistanceCut = 50;
  fVertexProximityCut = 10;
}

void shower::TrackShowerSeparationAlg::IdentifyTracksFromEventCentre(const std::vector<art::Ptr<recob::Track> >& tracks,
								     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
								     const art::FindManyP<recob::Track>& fmtsp) {

  /// Attempt to identify tracks in the event by using the centre of the event

  // Find the charge weighted centre of the entire event!
  // ---- except space points don't have charge (could look at associated hits, but cba!)
  TVector3 centre = TVector3(0,0,0);
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
    centre += TVector3((*spacePointIt)->XYZ()[0], (*spacePointIt)->XYZ()[1], (*spacePointIt)->XYZ()[2]);
  centre *= 1/spacePoints.size();

  TVector3 trackVertex, trackEnd, trackDirection;

  // Look at all tracks
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {

    trackVertex = ( ((*trackIt)->Vertex()-centre).Mag() > ((*trackIt)->End()-centre).Mag() ) ? (*trackIt)->Vertex() : (*trackIt)->End();
    trackEnd = ( ((*trackIt)->Vertex()-centre).Mag() > ((*trackIt)->End()-centre).Mag() ) ? (*trackIt)->End() : (*trackIt)->Vertex();
    trackDirection = ( ((*trackIt)->Vertex()-centre).Mag() > ((*trackIt)->End()-centre).Mag() ) ? (*trackIt)->VertexDirection() : (-1)*(*trackIt)->EndDirection();

    // Get the space points not associated with the current track
    std::vector<art::Ptr<recob::SpacePoint> > surroundingSpacePoints;
    GetSurroundingSpacePoints(spacePoints, surroundingSpacePoints, fmtsp, (*trackIt)->ID());

    // // TRUE -- use truth to find which particle this track is associated with --------
    // std::vector<art::Ptr<recob::Hit> > trackHits = fmh.at(trackIt->key());
    // int trueTrackID = FindTrueTrack(trackHits);
    // const simb::MCParticle* trueParticle = backtracker->TrackIDToParticle(trueTrackID);
    // pdg = trueParticle->PdgCode();
    // Length = (*trackIt)->Length();
    // TrackID = (*trackIt)->ID();
    // // -------------------------------------------------------------------------------

    bool showerLike = IdentifyShowerLikeTrack(trackEnd, trackDirection, surroundingSpacePoints);

    if (showerLike)
      continue;

    else
      fTrackLikeIDs.push_back((*trackIt)->ID());

  }

  return;

}

void shower::TrackShowerSeparationAlg::IdentifyTracksNearVertex(art::Ptr<recob::Vertex> const& vertex,
								std::vector<art::Ptr<recob::Track> > const& tracks,
								std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints,
								art::FindManyP<recob::Track> const& fmtsp) {

  /// Identifies hadron-like tracks which originate from near the interaction vertex

  double xyz[3];
  vertex->XYZ(xyz);
  TVector3 vertexPos = TVector3(xyz[0], xyz[1], xyz[2]);

  // Look at each track
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {

    TVector3 end, direction;

    // See if either end is close to the vertex
    if ( ((*trackIt)->Vertex() - vertexPos).Mag() < fVertexProximityCut or
	 ((*trackIt)->End() - vertexPos).Mag() < fVertexProximityCut) {
      end = ((*trackIt)->Vertex() - vertexPos).Mag() < ((*trackIt)->End() - vertexPos).Mag() ? (*trackIt)->End() : (*trackIt)->VertexDirection();
      direction = ((*trackIt)->Vertex() - vertexPos).Mag() < ((*trackIt)->End() - vertexPos).Mag() ? (*trackIt)->VertexDirection() : (-1)*(*trackIt)->VertexDirection();
    }

    else
      continue;

    // Get the space points not associated with the current track
    std::vector<art::Ptr<recob::SpacePoint> > surroundingSpacePoints;
    GetSurroundingSpacePoints(spacePoints, surroundingSpacePoints, fmtsp, (*trackIt)->ID());

    // Make sure this track start isn't a shower start
    bool showerLike = IdentifyShowerLikeTrack(end, direction, surroundingSpacePoints);

    if (showerLike)
      continue;

    // This track originates from near the interaction vertex and is not shower-like
    fTrackLikeIDs.push_back((*trackIt)->ID());

  }

  return;

}

bool shower::TrackShowerSeparationAlg::IdentifyShowerLikeTrack(TVector3 const& end,
							       TVector3 const& direction,
							       std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints) {

  /// Determines whether or not a track is actually the start of a shower

  std::vector<art::Ptr<recob::SpacePoint> > spacePointsInCone;
  GetSpacePointsInCone(spacePoints, spacePointsInCone, end, direction);

  if (spacePointsInCone.size() < 2)
    return false;

  // Get the average spread of these space points
  double spread = SpacePointSpread(spacePointsInCone);

  if (spread < 1)
    return false;

  return true;

}

void shower::TrackShowerSeparationAlg::FillHitsToCluster(const std::vector<art::Ptr<recob::Hit> >& initialHits,
							 std::vector<art::Ptr<recob::Hit> >& hitsToCluster,
							 const art::FindManyP<recob::Track>& fmt) {

  /// Fill the output container with all the hits not associated to track-like objects

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator initialHit = initialHits.begin(); initialHit != initialHits.end(); ++initialHit) {
    std::vector<art::Ptr<recob::Track> > showerTracks = fmt.at(initialHit->key());
    if ( showerTracks.size() and (std::find(fTrackLikeIDs.begin(), fTrackLikeIDs.end(), showerTracks.at(0)->ID()) == fTrackLikeIDs.end()) )
      hitsToCluster.push_back(*initialHit);
  }

  return;

}

int shower::TrackShowerSeparationAlg::FindTrackID(art::Ptr<recob::Hit> &hit) {
  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackID(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }
  return likelyTrackID;
}

int shower::TrackShowerSeparationAlg::FindTrueTrack(std::vector<art::Ptr<recob::Hit> > &trackHits) {
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
    art::Ptr<recob::Hit> hit = *trackHitIt;
    int trackID = FindTrackID(hit);
    trackMap[trackID] += hit->Integral();
  }
  //return std::max_element(trackMap.begin(), trackMap.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) {return p1.second < p2.second;} )->first;
  double highestCharge = 0;
  int clusterTrack = 0;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      clusterTrack  = trackIt->first;
    }
  return clusterTrack;
}

void shower::TrackShowerSeparationAlg::GetSurroundingSpacePoints(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
								 std::vector<art::Ptr<recob::SpacePoint> >& surroundingSpacePoints,
								 const art::FindManyP<recob::Track>& fmt,
								 unsigned int trackID) {

  /// Finds the space points surrounding the track but not part of it

  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {

    bool spacePointIsInCurrentTrack = false;
  
    std::vector<art::Ptr<recob::Track> > spacePointTracks = fmt.at(spacePointIt->key());
    for (std::vector<art::Ptr<recob::Track> >::iterator spacePointTrackIt = spacePointTracks.begin(); spacePointTrackIt != spacePointTracks.end(); ++spacePointTrackIt)
      if (spacePointTrackIt->key() == trackID) spacePointIsInCurrentTrack = true;
    
    if (!spacePointIsInCurrentTrack)
      surroundingSpacePoints.push_back(*spacePointIt);

  }

  return;

}

void shower::TrackShowerSeparationAlg::GetSpacePointsInCone(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
							    std::vector<art::Ptr<recob::SpacePoint> >& spacePointsInCone,
							    const TVector3& trackEnd,
							    const TVector3& trackDirection) {

  /// Look for space points near the track and within a narrow cone

  TVector3 spacePointPos, spacePointProj;
  double displacementFromAxis, distanceFromTrackEnd, projDistanceFromTrackEnd, angleFromTrackEnd;

  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {

    // Get the properties of this space point
    const double* xyz = (*spacePointIt)->XYZ();
    spacePointPos = TVector3(xyz[0], xyz[1], xyz[2]);
    spacePointProj = ((spacePointPos-trackEnd).Dot(trackDirection))*trackDirection + trackEnd;
    displacementFromAxis = (spacePointProj - spacePointPos).Mag();
    distanceFromTrackEnd = (spacePointPos - trackEnd).Mag();
    projDistanceFromTrackEnd = (spacePointProj - trackEnd).Mag();
    angleFromTrackEnd = TMath::ASin(displacementFromAxis/distanceFromTrackEnd) * 180 / TMath::Pi();

    if (projDistanceFromTrackEnd < fDistanceCut and angleFromTrackEnd < fAngleCut)
      spacePointsInCone.push_back(*spacePointIt);

  }

  return;

}

void shower::TrackShowerSeparationAlg::RemoveTrackHits(std::vector<art::Ptr<recob::Hit> > const& initialHits,
						       std::vector<art::Ptr<recob::Track> > const& tracks,
						       std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints,
						       std::vector<art::Ptr<recob::Vertex> > const& vertices,
						       art::FindManyP<recob::Track> const& fmth,
						       art::FindManyP<recob::Track> const& fmtsp,
						       art::FindManyP<recob::Hit> const& fmh,
						       std::vector<art::Ptr<recob::Hit> >& hitsToCluster,
						       int event,
						       int run) {

  Event = event;
  Run = run;

  /// Takes specific previously reconstructed quantites and removes hits which are considered track-like

  if (spacePoints.size() == 0) {
    hitsToCluster = initialHits;
    return;
  }

  // Find the vertex furthest upstream (if it exists)
  art::Ptr<recob::Vertex> vertex;
  for (std::vector<art::Ptr<recob::Vertex> >::const_iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt) {
    double xyzNew[3], xyzOld[3];
    (*vertexIt)->XYZ(xyzNew);
    if (vertex.isNull())
      vertex = *vertexIt;
    else {
      vertex->XYZ(xyzOld);
      if (xyzNew[2] < xyzOld[2])
	vertex = *vertexIt;
    }
  }

  // If we have a vertex then, for the time being, use it to remove tracks which are too close!
  if (!vertex.isNull())
    this->IdentifyTracksNearVertex(vertex, tracks, spacePoints, fmtsp);

  // If no vertex, things are harder
  else
    this->IdentifyTracksFromEventCentre(tracks, spacePoints, fmtsp);

  this->FillHitsToCluster(initialHits, hitsToCluster, fmth);

  return;

}

double shower::TrackShowerSeparationAlg::SpacePointSpread(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) {

  /// Finds the spread of a set of space points about their central axis

  // Find the centre of the space points
  TVector3 centre = TVector3(0,0,0);
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
    centre += TVector3((*spacePointIt)->XYZ()[0], (*spacePointIt)->XYZ()[1], (*spacePointIt)->XYZ()[2]);
  centre *= 1/spacePoints.size();

  // Find the central axis of the space points
  TPrincipal* pca = new TPrincipal(3,"");
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
    pca->AddRow((*spacePointIt)->XYZ());
  pca->MakePrincipals();
  const TMatrixD* eigenvectors = pca->GetEigenVectors();
  TVector3 spacePointDir = TVector3((*eigenvectors)[0][0], (*eigenvectors)[1][0], (*eigenvectors)[2][0]);
  delete pca;

  // See if the space points form something which may resemble a track (i.e. straight line) or a shower
  double avDistanceFromCentralAxis = 0;
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {
    const double* xyz = (*spacePointIt)->XYZ();
    TVector3 spacePointPos = TVector3(xyz[0], xyz[1], xyz[2]);
    TVector3 projectionOntoCentralAxis = ((spacePointPos-centre).Dot(spacePointDir))*spacePointDir + centre;
    avDistanceFromCentralAxis += (spacePointPos - projectionOntoCentralAxis).Mag();
  }
  avDistanceFromCentralAxis /= spacePoints.size();

  return avDistanceFromCentralAxis;

}
