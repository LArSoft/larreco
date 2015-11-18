////////////////////////////////////////////////////////////////////////
// Class: TrackShowerSeparationAlg
// File:  TrackShowerSeparationAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Track/shower discrimination class
////////////////////////////////////////////////////////////////////////

#ifndef TrackShowerSeparationAlg_hxx
#define TrackShowerSeparationAlg_hxx

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Core/FindManyP.h"

// larsoft
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Vertex.h"

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

  void reconfigure(fhicl::ParameterSet const& pset);
  void IdentifyTracksFromEventCentre(const std::vector<art::Ptr<recob::Track> >& tracks, const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, art::FindManyP<recob::Track> const& fmt);
  void IdentifyTracksNearTracks(std::vector<art::Ptr<recob::Track> > const& tracks);
  void IdentifyTracksNearVertex(art::Ptr<recob::Vertex> const& vertex, std::vector<art::Ptr<recob::Track> > const& tracks, std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints, art::FindManyP<recob::Track> const& fmt);
  bool IdentifyShowerLikeTrack(TVector3 const& end, TVector3 const& direction, std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints);
  void FillHitsToCluster(std::vector<art::Ptr<recob::Hit> > const& initialHits, std::vector<art::Ptr<recob::Hit> >& hitsToCluster, art::FindManyP<recob::Track> const& fmt);
  int FindTrackID(art::Ptr<recob::Hit> &hit);
  int FindTrueTrack(std::vector<art::Ptr<recob::Hit> > &trackHits);
  void GetSurroundingSpacePoints(const std::vector<art::Ptr<recob::SpacePoint> >& allsp, std::vector<art::Ptr<recob::SpacePoint> >& sp, const art::FindManyP<recob::Track>& fmt, unsigned int trackID);
  void GetSpacePointsInCone(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, std::vector<art::Ptr<recob::SpacePoint> >& spacePointsInCone, const TVector3& trackEnd, const TVector3& trackDirection);
  void RemoveTrackHits(std::vector<art::Ptr<recob::Hit> > const& initialHits,
		       std::vector<art::Ptr<recob::Track> > const& tracks,
		       std::vector<art::Ptr<recob::SpacePoint> > const& spacePoints,
		       std::vector<art::Ptr<recob::Vertex> > const& vertices,
		       art::FindManyP<recob::Track> const& fmth,
		       art::FindManyP<recob::Track> const& fmtsp,
		       art::FindManyP<recob::Hit> const& fmh,
		       std::vector<art::Ptr<recob::Hit> >& hitsToCluster,
		       int event,
		       int run);
  double SpacePointSpread(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints);

 private:

  std::vector<int> fTrackLikeIDs, fShowerLikeIDs;

  double fAngleCut, fDistanceCut, fVertexProximityCut;

  art::ServiceHandle<cheat::BackTracker> backtracker;
  art::ServiceHandle<art::TFileService> tfs;

  TTree* ftree;
  double Distance, Angle, Length, AvDistance;
  int Event, Run, TrackID, pdg, NSpacePoints;

};

#endif
