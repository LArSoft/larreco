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
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
namespace fhicl { class ParameterSet; }

// larsoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// ROOT
#include "TMath.h"
#include "TVector3.h"

namespace shower {
  class TrackShowerSeparationAlg;
  class ReconTrack;
}

#include <algorithm>
#include <iterator>
#include <map>
#include <vector>

class shower::ReconTrack {
 public:

  ReconTrack(int id) {
    fID = id;
    fTrack = false;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;
  }

  // Setters
  void SetVertex(TVector3 vertex) { fVertex = vertex; }
  void SetEnd(TVector3 end) { fEnd = end; }
  void SetLength(double length) { fLength = length; }
  void SetVertexDir(TVector3 vertexDir) { fVertexDir = vertexDir; }
  void SetDirection(TVector3 direction) { fDirection = direction; }
  void SetHits(std::vector<art::Ptr<recob::Hit> > hits) { fHits = hits; }
  void SetSpacePoints(std::vector<art::Ptr<recob::SpacePoint> > spacePoints) { fSpacePoints = spacePoints; }

  void AddForwardTrack(int track) { if (std::find(fForwardConeTracks.begin(), fForwardConeTracks.end(), track) == fForwardConeTracks.end()) fForwardConeTracks.push_back(track); }
  void AddBackwardTrack(int track) { if (std::find(fBackwardConeTracks.begin(), fBackwardConeTracks.end(), track) == fBackwardConeTracks.end()) fBackwardConeTracks.push_back(track); }
  void AddShowerTrack(int track) { fShowerTracks.push_back(track); }

  void AddForwardSpacePoint(int spacePoint) { fForwardSpacePoints.push_back(spacePoint); }
  void AddBackwardSpacePoint(int spacePoint) { fBackwardSpacePoints.push_back(spacePoint); }
  void AddCylinderSpacePoint(int spacePoint) { fCylinderSpacePoints.push_back(spacePoint); }
  void AddSphereSpacePoint(int spacePoint) { fSphereSpacePoints.push_back(spacePoint); }
  void AddIsolationSpacePoint(int spacePoint, double distance) { fIsolationSpacePoints[spacePoint] = distance; }

  // Getters
  int ID() const { return fID; }
  TVector3 Vertex() const { return fVertex; }
  TVector3 End() const { return fEnd; }
  double Length() const { return fLength; }
  TVector3 VertexDirection() const { return fVertexDir; }
  TVector3 Direction() const { return fDirection; }
  const std::vector<art::Ptr<recob::Hit> >& Hits() const { return fHits; }
  const std::vector<art::Ptr<recob::SpacePoint> >& SpacePoints() const { return fSpacePoints; }

  void FlipTrack() {
    TVector3 tmp = fEnd;
    fEnd = fVertex;
    fVertex = tmp;
    fDirection *= -1;
  }

  void MakeShower() {
    if (fTrack)
      this->MakeShowerTrack();
    else
      this->MakeShowerCone();
  }
  void MakeShowerTrack() {
    fShower = true;
    fShowerTrack = true;
    fShowerCone = false;
    fTrack = false;
  }
  void MakeShowerCone() {
    fShower = true;
    fShowerCone = true;
    fShowerTrack = false;
    fTrack = false;
  }
  void MakeTrack() {
    fTrack = true;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;
  }

  bool IsShower() const { return fShower; }
  bool IsShowerTrack() const { return fShowerTrack; }
  bool IsShowerCone() const { return fShowerCone; }
  bool IsTrack() const { return fTrack; }
  bool IsUndetermined() const { return !fTrack and !fShower; }

  int TrackConeSize() const { return (int)fForwardConeTracks.size() - (int)fBackwardConeTracks.size(); }
  bool ShowerTrackCandidate() const { return TrackConeSize() > 5; }
  const std::vector<int>& ShowerTracks() const { return fShowerTracks; }
  const std::vector<int>& ForwardConeTracks() const { return fForwardConeTracks; }

  int ConeSize() const { return (int)fForwardSpacePoints.size() - (int)fBackwardSpacePoints.size(); }
  int ForwardSpacePoints() const { return fForwardSpacePoints.size(); }
  int NumCylinderSpacePoints() const { return fCylinderSpacePoints.size(); }
  double CylinderSpacePointRatio() const { return (double)fCylinderSpacePoints.size()/(double)fSpacePoints.size(); }
  int NumSphereSpacePoints() const { return fSphereSpacePoints.size(); }
  double SphereSpacePointDensity(double scale) const { return (double)fSphereSpacePoints.size()/(4*TMath::Pi()*TMath::Power((scale*fLength/2.),3)/3.); }
  double IsolationSpacePointDistance() const { std::vector<double> distances;
    std::transform(fIsolationSpacePoints.begin(), fIsolationSpacePoints.end(), std::back_inserter(distances), [](const std::pair<int,double>& p){return p.second;});
    return TMath::Mean(distances.begin(), distances.end()); }

 private:

  int fID;
  TVector3 fVertex;
  TVector3 fEnd;
  double fLength;
  TVector3 fVertexDir;
  TVector3 fDirection;
  std::vector<art::Ptr<recob::Hit> > fHits;
  std::vector<art::Ptr<recob::SpacePoint> > fSpacePoints;

  std::vector<int> fForwardConeTracks;
  std::vector<int> fBackwardConeTracks;
  std::vector<int> fShowerTracks;

  std::vector<int> fForwardSpacePoints;
  std::vector<int> fBackwardSpacePoints;
  std::vector<int> fCylinderSpacePoints;
  std::vector<int> fSphereSpacePoints;
  std::map<int,double> fIsolationSpacePoints;

  bool fShower;
  bool fShowerTrack;
  bool fShowerCone;
  bool fTrack;

};

class shower::TrackShowerSeparationAlg {
 public:

  TrackShowerSeparationAlg(fhicl::ParameterSet const& pset);

  /// Read in configurable parameters from provided parameter set
  void reconfigure(fhicl::ParameterSet const& pset);

  std::vector<art::Ptr<recob::Hit> > SelectShowerHits(int event,
						      const std::vector<art::Ptr<recob::Hit> >& hits,
						      const std::vector<art::Ptr<recob::Track> >& tracks,
						      const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
						      const art::FindManyP<recob::Hit>& fmht,
						      const art::FindManyP<recob::Track>& fmth,
						      const art::FindManyP<recob::SpacePoint>& fmspt,
                                                      const art::FindManyP<recob::Track>& fmtsp) const;

 private:

  ///
  std::vector<int> InitialTrackLikeSegment(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks) const;

  ///
  TVector3 Gradient(const std::vector<TVector3>& points, const std::unique_ptr<TVector3>& dir) const;

  ///
  TVector3 Gradient(const art::Ptr<recob::Track>& track) const;

  ///
  TVector3 Gradient(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) const;

  /// Projects the 3D point given along the line defined by the specified direction
  /// Coordinate system is assumed to be centred at the origin unless a difference point is specified
  TVector3 ProjPoint(const TVector3& point, const TVector3& direction, const TVector3& origin = TVector3(0,0,0)) const;

  /// Return 3D point of this space point
  TVector3 SpacePointPos(const art::Ptr<recob::SpacePoint>& spacePoint) const;

  ///
  double SpacePointsRMS(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) const;

  // Parameters
  int fDebug;

  // Properties
  double fConeAngle;
  double fCylinderRadius;

  // Cuts
  double fTrackVertexCut;
  double fCylinderCut;
  double fShowerConeCut;

};

#endif
