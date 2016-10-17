////////////////////////////////////////////////////////////////////////
// Class: TrackShowerSeparationAlg
// File:  TrackShowerSeparationAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Track/shower separation class.
// Provides methods for removing hits associated with track-like
// objects.
// To be run after track reconstruction, before shower reconstruction.
////////////////////////////////////////////////////////////////////////

#ifndef TrackShowerSeparationAlg_hxx
#define TrackShowerSeparationAlg_hxx

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// larsoft
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"

// ROOT
#include "TTree.h"
#include "TMath.h"
#include "TPrincipal.h"

namespace shower {
  class TrackShowerSeparationAlg;
}

class shower::TrackShowerSeparationAlg {
 public:

  TrackShowerSeparationAlg(fhicl::ParameterSet const& pset);

  /// Read in configurable parameters from provided parameter set
  void reconfigure(fhicl::ParameterSet const& pset);

  std::vector<art::Ptr<recob::Hit> > RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& hits,
						     const std::vector<art::Ptr<recob::Track> >& tracks,
						     const art::FindManyP<recob::Hit>& fmh);

  // --------------------------- OLD (late 2015) -------------------------------

  /// Takes specific previously reconstructed quantites and removes hits which are considered track-like
  /// Returns a vector of hits which are determined to be shower-like
  std::vector<art::Ptr<recob::Hit> > RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& initialHits,
						     const std::vector<art::Ptr<recob::Track> >& tracks,
						     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
						     const std::vector<art::Ptr<recob::Vertex> >& vertices,
						     const art::FindManyP<recob::Track>& fmth,
						     const art::FindManyP<recob::Track>& fmtsp,
						     const art::FindManyP<recob::Hit>& fmh,
						     int event,
						     int run);

  /// Uses information from Pandora reconstruction to separate track-like and shower-like hits
  /// Returns a vector of hits which are determined to be shower-like
  std::vector<art::Ptr<recob::Hit> > RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& hits,
						     const std::vector<art::Ptr<recob::PFParticle> > pfParticles,
						     const art::FindManyP<recob::Cluster>& fmc,
						     const art::FindManyP<recob::Hit>& fmh);

 private:

  int InitialTrackLikeSegment(const std::vector<art::Ptr<recob::Track> >& tracks, std::vector<int>& showerTracks, std::vector<int>& trackTracks);

  // --------------------------- OLD (late 2015) -------------------------------

  /// Fill the output container with all the hits not associated with track-like objects
  std::vector<art::Ptr<recob::Hit> > FillHitsToCluster(const std::vector<art::Ptr<recob::Hit> >& initialHits,
						       const art::FindManyP<recob::Track>& fmt);

  /// Find the true track most likely associated with this hit
  int FindTrackID(const art::Ptr<recob::Hit>& hit);

  /// Find the true track most likely associated with this set of hits
  int FindTrueTrack(const std::vector<art::Ptr<recob::Hit> >& trackHits);

  /// Finds the space points surrounding the track but not part of it
  std::vector<art::Ptr<recob::SpacePoint> > GetSurroundingSpacePoints(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
								      const art::FindManyP<recob::Track>& fmt,
								      unsigned int trackID);

  /// Look for space points near the track and within a narrow cone
  std::vector<art::Ptr<recob::SpacePoint> > GetSpacePointsInCone(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
								 const TVector3& trackEnd,
								 const TVector3& trackDirection);

  /// Determines whether or not a track is actually the start of a shower
  bool IdentifyShowerLikeTrack(const TVector3& end,
			       const TVector3& direction,
			       const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints);

  /// Attempt to identify tracks in the event by using the centre of the event
  void IdentifyTracksFromEventCentre(const std::vector<art::Ptr<recob::Track> >& tracks,
				     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
				     const art::FindManyP<recob::Track>& fmtsp);

  /// Identifies tracks which start just after previously identified tracks end
  void IdentifyTracksNearTracks(const std::vector<art::Ptr<recob::Track> >& tracks);

  /// Identifies hadron-like tracks which originate from near the interaction vertex
  void IdentifyTracksNearVertex(const art::Ptr<recob::Vertex>& vertex,
				const std::vector<art::Ptr<recob::Track> >& tracks,
				const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
				const art::FindManyP<recob::Track>& fmtsp);

  /// Finds the spread of a set of space points about their central axis
  double SpacePointSpread(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints);

  std::vector<int> fTrackLikeIDs, fShowerLikeIDs;

  // Configurable parameters
  double fAngleCut, fDistanceCut, fVertexProximityCut, fTrackProximityCut, fAvTrackHitDistance;

  art::ServiceHandle<cheat::BackTracker> backtracker;
  art::ServiceHandle<art::TFileService> tfs;

  TTree* ftree;
  double Distance, Angle, Length, AvDistance;
  int Event, Run, TrackID, pdg, NSpacePoints;

};

#endif
