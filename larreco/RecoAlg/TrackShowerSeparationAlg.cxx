////////////////////////////////////////////////////////////////////////
// Class: TrackShowerSeparationAlg
// File:  TrackShowerSeparationAlg.cxx
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Track/shower separation class.
// Provides methods for removing hits associated with track-like
// objects.
// To be run after track reconstruction, before shower reconstruction.
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrackShowerSeparationAlg.h"
#include "fhiclcpp/ParameterSet.h"

#include "TMathBase.h"
#include "TVector2.h"

shower::TrackShowerSeparationAlg::TrackShowerSeparationAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

void shower::TrackShowerSeparationAlg::reconfigure(fhicl::ParameterSet const& pset)
{
  fConeAngle = pset.get<double>("ConeAngle");
  fCylinderRadius = pset.get<double>("CylinderRadius");
  fTrackVertexCut = pset.get<double>("TrackVertexCut");
  fCylinderCut = pset.get<double>("CylinderCut");
  fShowerConeCut = pset.get<double>("ShowerConeCut");

  fDebug = pset.get<int>("Debug", 0);
}

std::vector<art::Ptr<recob::Hit>> shower::TrackShowerSeparationAlg::SelectShowerHits(
  int event,
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<art::Ptr<recob::Track>>& tracks,
  const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
  const art::FindManyP<recob::Hit>& fmht,
  const art::FindManyP<recob::Track>& fmth,
  const art::FindManyP<recob::SpacePoint>& fmspt,
  const art::FindManyP<recob::Track>& fmtsp) const
{

  // Ok, here we are again
  // Playing the game in which no one wins, but everyone loses
  // Trying to separate tracks from showers
  //    Ode to track/shower separation
  //    M Wallbank, Oct 2016

  // Showers are comprised of two topologically different parts:
  //  -- the 'shower track' (the initial part of the shower)
  //  -- the 'shower cone' (the part of the shower after the initial track when showering happens)

  //                                            -
  //                                          --
  //                 -             --       --
  //               --         ---         ---
  //             --      ----       --
  //  -----------   -----       -  --        ---
  //                    ---             ---
  //                       --               ---
  //  {=========}{==============================}
  // shower track          shower cone

  std::map<int, std::unique_ptr<ReconTrack>> reconTracks;
  for (std::vector<art::Ptr<recob::Track>>::const_iterator trackIt = tracks.begin();
       trackIt != tracks.end();
       ++trackIt) {
    std::unique_ptr<ReconTrack> track = std::make_unique<ReconTrack>(trackIt->key());
    track->SetVertex((*trackIt)->Vertex<TVector3>());
    track->SetEnd((*trackIt)->End<TVector3>());
    track->SetVertexDir((*trackIt)->VertexDirection<TVector3>());
    track->SetLength((*trackIt)->Length());
    track->SetDirection(Gradient(*trackIt));
    track->SetHits(fmht.at(trackIt->key()));
    track->SetSpacePoints(fmspt.at(trackIt->key()));
    const std::vector<art::Ptr<recob::SpacePoint>> spsss = fmspt.at(trackIt->key());
    reconTracks[trackIt->key()] = std::move(track);
  }

  // Consider the space point cylinder situation
  double avCylinderSpacePoints = 0;
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    // Get the 3D properties of the track
    TVector3 point = trackIt->second->Vertex();
    TVector3 direction = trackIt->second->Direction();

    // Count space points in the volume around the track
    for (std::vector<art::Ptr<recob::SpacePoint>>::const_iterator spacePointIt =
           spacePoints.begin();
         spacePointIt != spacePoints.end();
         ++spacePointIt) {
      const std::vector<art::Ptr<recob::Track>> spTracks = fmtsp.at(spacePointIt->key());
      if (find_if(spTracks.begin(), spTracks.end(), [&trackIt](const art::Ptr<recob::Track>& t) {
            return (int)t.key() == trackIt->first;
          }) != spTracks.end())
        continue;
      // Get the properties of this space point
      TVector3 pos = SpacePointPos(*spacePointIt);
      TVector3 proj = ProjPoint(pos, direction, point);
      if ((pos - proj).Mag() < fCylinderRadius)
        trackIt->second->AddCylinderSpacePoint(spacePointIt->key());
    }
    avCylinderSpacePoints += trackIt->second->CylinderSpacePointRatio();
  }
  avCylinderSpacePoints /= (double)reconTracks.size();

  if (fDebug > 1) {
    std::cout << std::endl << "Cylinder space point ratios:" << std::endl;
    for (std::map<int, std::unique_ptr<ReconTrack>>::const_iterator trackIt = reconTracks.begin();
         trackIt != reconTracks.end();
         ++trackIt)
      std::cout << "  Track " << trackIt->first << " has cylinder space point ratio "
                << trackIt->second->CylinderSpacePointRatio() << " (event average "
                << avCylinderSpacePoints << ")" << std::endl;
  }

  // Identify tracks
  if (fDebug > 0) std::cout << std::endl << "Identifying tracks:" << std::endl;
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt)
    if (trackIt->second->CylinderSpacePointRatio() / avCylinderSpacePoints < fCylinderCut) {
      if (fDebug > 0)
        std::cout << "  Making track " << trackIt->first << " a track (Type I)" << std::endl;
      trackIt->second->MakeTrack();
    }

  // Identify further tracks
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    if (trackIt->second->IsTrack()) continue;
    for (std::map<int, std::unique_ptr<ReconTrack>>::const_iterator otherTrackIt =
           reconTracks.begin();
         otherTrackIt != reconTracks.end();
         ++otherTrackIt) {
      if (trackIt->first == otherTrackIt->first or !otherTrackIt->second->IsTrack()) continue;
      if ((trackIt->second->Vertex() - otherTrackIt->second->Vertex()).Mag() < fTrackVertexCut or
          (trackIt->second->Vertex() - otherTrackIt->second->End()).Mag() < fTrackVertexCut or
          (trackIt->second->End() - otherTrackIt->second->Vertex()).Mag() < fTrackVertexCut or
          (trackIt->second->End() - otherTrackIt->second->End()).Mag() < fTrackVertexCut) {
        if (fDebug > 0)
          std::cout << "  Making track " << trackIt->first << " a track (Type II)" << std::endl;
        trackIt->second->MakeTrack();
      }
    }
  }

  // Consider removing false tracks by looking at their closest approach to any other track

  // Consider the space point cone situation
  std::vector<art::Ptr<recob::SpacePoint>> showerSpacePoints;
  for (std::vector<art::Ptr<recob::SpacePoint>>::const_iterator spacePointIt = spacePoints.begin();
       spacePointIt != spacePoints.end();
       ++spacePointIt) {
    bool showerSpacePoint = true;
    const std::vector<art::Ptr<recob::Track>> spacePointTracks = fmtsp.at(spacePointIt->key());
    for (std::vector<art::Ptr<recob::Track>>::const_iterator trackIt = spacePointTracks.begin();
         trackIt != spacePointTracks.end();
         ++trackIt)
      if (reconTracks[trackIt->key()]->IsTrack()) showerSpacePoint = false;
    if (showerSpacePoint) showerSpacePoints.push_back(*spacePointIt);
  }

  // Identify tracks which slipped through and shower tracks
  // For the moment, until the track tagging gets better at least, don't try to identify tracks from this
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    for (std::vector<art::Ptr<recob::SpacePoint>>::const_iterator spacePointIt =
           showerSpacePoints.begin();
         spacePointIt != showerSpacePoints.end();
         ++spacePointIt) {
      bool associatedSpacePoint = false;
      const std::vector<art::Ptr<recob::Track>> spTracks = fmtsp.at(spacePointIt->key());
      for (std::vector<art::Ptr<recob::Track>>::const_iterator spTrackIt = spTracks.begin();
           spTrackIt != spTracks.end();
           ++spTrackIt)
        if ((int)spTrackIt->key() == trackIt->first) associatedSpacePoint = true;
      if (associatedSpacePoint) continue;
      if ((SpacePointPos(*spacePointIt) - trackIt->second->Vertex())
            .Angle(trackIt->second->Direction()) < fConeAngle * TMath::Pi() / 180) {
        trackIt->second->AddForwardSpacePoint(spacePointIt->key());
        trackIt->second->AddForwardTrack(spTracks.at(0).key());
      }
      if ((SpacePointPos(*spacePointIt) - trackIt->second->Vertex())
            .Angle(-1 * trackIt->second->Direction()) < fConeAngle * TMath::Pi() / 180) {
        trackIt->second->AddBackwardSpacePoint(spacePointIt->key());
        trackIt->second->AddBackwardTrack(spTracks.at(0).key());
      }
    }
  }
  if (fDebug > 0) std::cout << std::endl << "Identifying showers:" << std::endl;
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    if (fDebug > 1)
      std::cout << "    Track " << trackIt->first << " has space point cone size "
                << trackIt->second->ConeSize() << " and track cone size "
                << trackIt->second->TrackConeSize() << std::endl;
    if (TMath::Abs(trackIt->second->ConeSize()) > 30 and
        TMath::Abs(trackIt->second->TrackConeSize()) > 3) {
      trackIt->second->MakeShower();
      if (fDebug > 0)
        std::cout << "  Making track " << trackIt->first << " a shower (Type I)" << std::endl;
      if (trackIt->second->ConeSize() < 0) trackIt->second->FlipTrack();
    }
  }

  // Look for shower cones
  std::cout << std::endl << "  Shower cones:" << std::endl;
  for (std::map<int, std::unique_ptr<ReconTrack>>::const_iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    if (trackIt->second->IsShower()) {
      if (fDebug > 1) std::cout << "    Track " << trackIt->first << std::endl;
      for (std::map<int, std::unique_ptr<ReconTrack>>::iterator otherTrackIt = reconTracks.begin();
           otherTrackIt != reconTracks.end();
           ++otherTrackIt) {
        if (trackIt->first == otherTrackIt->first or !otherTrackIt->second->IsUndetermined())
          continue;
        if ((otherTrackIt->second->Vertex() - trackIt->second->Vertex())
                .Angle(trackIt->second->Direction()) < fConeAngle * TMath::Pi() / 180 or
            (otherTrackIt->second->End() - trackIt->second->Vertex())
                .Angle(trackIt->second->Direction()) < fConeAngle * TMath::Pi() / 180) {
          std::cout << "      " << otherTrackIt->first << std::endl;
          otherTrackIt->second->MakeShower();
          if (fDebug > 0)
            std::cout << "      Making track " << otherTrackIt->first << " a shower (Type II)"
                      << std::endl;
        }
      }
    }
  }

  // Look at remaining undetermined tracks
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    if (trackIt->second->IsUndetermined()) { trackIt->second->MakeShower(); }
  }

  std::cout << std::endl << "Event " << event << " track shower separation:" << std::endl;
  std::cout << "Shower initial tracks are:" << std::endl;
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt)
    if (trackIt->second->IsShowerTrack()) std::cout << "  " << trackIt->first << std::endl;

  std::cout << "Track tracks are:" << std::endl;
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt)
    if (trackIt->second->IsTrack()) std::cout << "  " << trackIt->first << std::endl;

  // Shower hits -- select all hits which aren't associated with a determined track
  std::vector<art::Ptr<recob::Hit>> showerHits;
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator hitIt = hits.begin(); hitIt != hits.end();
       ++hitIt) {
    bool showerHit = true;
    const std::vector<art::Ptr<recob::Track>> hitTracks = fmth.at(hitIt->key());
    for (std::vector<art::Ptr<recob::Track>>::const_iterator hitTrackIt = hitTracks.begin();
         hitTrackIt != hitTracks.end();
         ++hitTrackIt)
      if (reconTracks[hitTrackIt->key()]->IsTrack()) showerHit = false;
    if (showerHit) showerHits.push_back(*hitIt);
  }

  return showerHits;
}

std::vector<int> shower::TrackShowerSeparationAlg::InitialTrackLikeSegment(
  std::map<int, std::unique_ptr<ReconTrack>>& reconTracks) const
{

  // Consider the cones for each track
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {

    for (std::map<int, std::unique_ptr<ReconTrack>>::iterator otherTrackIt = reconTracks.begin();
         otherTrackIt != reconTracks.end();
         ++otherTrackIt) {

      if (trackIt->first == otherTrackIt->first) continue;
      if ((otherTrackIt->second->Vertex() - trackIt->second->Vertex())
              .Angle(trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180 or
          (otherTrackIt->second->End() - trackIt->second->Vertex())
              .Angle(trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180) {
        trackIt->second->AddForwardTrack(otherTrackIt->first);
        otherTrackIt->second->AddShowerTrack(trackIt->first);
      }
      if ((otherTrackIt->second->Vertex() - trackIt->second->Vertex())
              .Angle(-1 * trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180 or
          (otherTrackIt->second->End() - trackIt->second->Vertex())
              .Angle(-1 * trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180)
        trackIt->second->AddBackwardTrack(otherTrackIt->first);
    }
  }

  // Determine if any of these tracks are actually shower tracks
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {

    if (!trackIt->second->ShowerTrackCandidate()) continue;

    std::cout << "Track " << trackIt->first
              << " is a candidate, with shower cone tracks:" << std::endl;
    const std::vector<int>& showerConeTracks = trackIt->second->ForwardConeTracks();
    for (std::vector<int>::const_iterator showerConeTrackIt = showerConeTracks.begin();
         showerConeTrackIt != showerConeTracks.end();
         ++showerConeTrackIt)
      std::cout << "  " << *showerConeTrackIt << std::endl;

    bool isBestCandidate = true;
    const std::vector<int>& showerTracks = trackIt->second->ShowerTracks();
    for (std::vector<int>::const_iterator showerTrackIt = showerTracks.begin();
         showerTrackIt != showerTracks.end();
         ++showerTrackIt) {
      if (!reconTracks[*showerTrackIt]->ShowerTrackCandidate()) continue;
      if (std::find(showerConeTracks.begin(), showerConeTracks.end(), *showerTrackIt) ==
          showerConeTracks.end())
        continue;
      if (reconTracks[*showerTrackIt]->IsShowerTrack()) continue;
      if (trackIt->second->TrackConeSize() < reconTracks[*showerTrackIt]->TrackConeSize())
        isBestCandidate = false;
    }

    if (isBestCandidate) trackIt->second->MakeShowerTrack();
  }

  // Determine which tracks are shower cones
  std::vector<int> showerTracks;
  for (std::map<int, std::unique_ptr<ReconTrack>>::iterator trackIt = reconTracks.begin();
       trackIt != reconTracks.end();
       ++trackIt) {
    if (trackIt->second->IsShowerTrack()) {
      showerTracks.push_back(trackIt->first);
      const std::vector<int>& coneTracks = trackIt->second->ForwardConeTracks();
      for (std::vector<int>::const_iterator coneTrackIt = coneTracks.begin();
           coneTrackIt != coneTracks.end();
           ++coneTrackIt)
        reconTracks[*coneTrackIt]->MakeShowerCone();
    }
  }

  return showerTracks;
}

TVector3 shower::TrackShowerSeparationAlg::Gradient(const std::vector<TVector3>& points,
                                                    const std::unique_ptr<TVector3>& dir) const
{

  int nhits = 0;
  double sumx = 0., sumy = 0., sumz = 0., sumx2 = 0., sumxy = 0., sumxz = 0.;
  for (std::vector<TVector3>::const_iterator pointIt = points.begin(); pointIt != points.end();
       ++pointIt) {
    ++nhits;
    sumx += pointIt->X();
    sumy += pointIt->Y();
    sumz += pointIt->Z();
    sumx2 += pointIt->X() * pointIt->X();
    sumxy += pointIt->X() * pointIt->Y();
    sumxz += pointIt->X() * pointIt->Z();
  }

  double dydx = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  double yint = (sumy * sumx2 - sumx * sumxy) / (nhits * sumx2 - sumx * sumx);
  double dzdx = (nhits * sumxz - sumx * sumz) / (nhits * sumx2 - sumx * sumx);
  double zint = (sumz * sumx2 - sumx * sumxz) / (nhits * sumx2 - sumx * sumx);
  TVector2 directionXY = TVector2(1, dydx).Unit(), directionXZ = TVector2(1, dzdx).Unit();
  TVector3 direction = TVector3(1, dydx, dzdx).Unit();
  TVector3 intercept = TVector3(0, yint, zint);

  // Make sure the best fit direction is pointing correctly
  if (dir and TMath::Abs(direction.Angle(*dir)) > TMath::Pi() / 2.) direction *= -1;

  return direction;
}

TVector3 shower::TrackShowerSeparationAlg::Gradient(const art::Ptr<recob::Track>& track) const
{

  std::vector<TVector3> points;
  std::unique_ptr<TVector3> dir;

  for (unsigned int traj = 0; traj < track->NumberTrajectoryPoints(); ++traj)
    points.push_back(track->LocationAtPoint<TVector3>(traj));
  dir = std::make_unique<TVector3>(track->VertexDirection<TVector3>());

  return Gradient(points, dir);
}

TVector3 shower::TrackShowerSeparationAlg::Gradient(
  const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints) const
{

  std::vector<TVector3> points;
  std::unique_ptr<TVector3> dir;

  for (std::vector<art::Ptr<recob::SpacePoint>>::const_iterator spacePointIt = spacePoints.begin();
       spacePointIt != spacePoints.end();
       ++spacePointIt)
    points.push_back(SpacePointPos(*spacePointIt));

  return Gradient(points, dir);
}

TVector3 shower::TrackShowerSeparationAlg::ProjPoint(const TVector3& point,
                                                     const TVector3& direction,
                                                     const TVector3& origin) const
{
  return (point - origin).Dot(direction) * direction + origin;
}

TVector3 shower::TrackShowerSeparationAlg::SpacePointPos(
  const art::Ptr<recob::SpacePoint>& spacePoint) const
{
  const double* xyz = spacePoint->XYZ();
  return TVector3(xyz[0], xyz[1], xyz[2]);
}

double shower::TrackShowerSeparationAlg::SpacePointsRMS(
  const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints) const
{

  TVector3 point = SpacePointPos(spacePoints.at(0));
  TVector3 direction = Gradient(spacePoints);

  std::vector<double> distances;
  for (std::vector<art::Ptr<recob::SpacePoint>>::const_iterator spacePointIt = spacePoints.begin();
       spacePointIt != spacePoints.end();
       ++spacePointIt) {
    TVector3 pos = SpacePointPos(*spacePointIt);
    TVector3 proj = ProjPoint(pos, direction, point);
    distances.push_back((pos - proj).Mag());
  }

  double rms = TMath::RMS(distances.begin(), distances.end());

  return rms;
}
