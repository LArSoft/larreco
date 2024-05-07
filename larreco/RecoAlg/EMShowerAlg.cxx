///////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#include "cetlib/container_algorithms.h"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcorealg/Geometry/Exceptions.h"
#include "lardata/ArtDataHelper/ToElement.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "larreco/RecoAlg/EMShowerAlg.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
#include "TMath.h"
#include "TMathBase.h"
#include "TMultiGraph.h"
#include "TProfile.h"

#include "range/v3/numeric.hpp"
#include "range/v3/view.hpp"

#include <limits>

using lar::to_element;
using namespace ranges;

shower::EMShowerAlg::EMShowerAlg(fhicl::ParameterSet const& pset, int const debug)
  : fDebug{debug}
  , fMinTrackLength{pset.get<double>("MinTrackLength")}
  , fdEdxTrackLength{pset.get<double>("dEdxTrackLength")}
  , fSpacePointSize{pset.get<double>("SpacePointSize")}
  , fNfitpass{pset.get<unsigned int>("Nfitpass")}
  , fNfithits{pset.get<std::vector<unsigned int>>("Nfithits")}
  , fToler{pset.get<std::vector<double>>("Toler")}
  , fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg"))
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  , fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg"))
  , fDetector{pset.get<std::string>("Detector", "dune35t")} // tmp
  , fMakeGradientPlot{pset.get<bool>("MakeGradientPlot", false)}
  , fMakeRMSGradientPlot{pset.get<bool>("MakeRMSGradientPlot", false)}
  , fNumShowerSegments{pset.get<int>("NumShowerSegments", 4)}
{
  if (fNfitpass != fNfithits.size() || fNfitpass != fToler.size()) {
    throw art::Exception(art::errors::Configuration)
      << "EMShowerAlg: fNfithits and fToler need to have size fNfitpass";
  }
}

void shower::EMShowerAlg::AssociateClustersAndTracks(
  std::vector<art::Ptr<recob::Cluster>> const& clusters,
  art::FindManyP<recob::Hit> const& fmh,
  art::FindManyP<recob::Track> const& fmt,
  std::map<int, std::vector<int>>& clusterToTracks,
  std::map<int, std::vector<int>>& trackToClusters) const
{
  std::vector<int> clustersToIgnore = {-999};
  AssociateClustersAndTracks(
    clusters, fmh, fmt, clustersToIgnore, clusterToTracks, trackToClusters);
}

void shower::EMShowerAlg::AssociateClustersAndTracks(
  std::vector<art::Ptr<recob::Cluster>> const& clusters,
  art::FindManyP<recob::Hit> const& fmh,
  art::FindManyP<recob::Track> const& fmt,
  std::vector<int> const& clustersToIgnore,
  std::map<int, std::vector<int>>& clusterToTracks,
  std::map<int, std::vector<int>>& trackToClusters) const
{
  // Look through all the clusters
  for (auto const& clusterPtr : clusters) {

    // Get the hits in this cluster
    auto const& clusterHits = fmh.at(clusterPtr.key());

    // Look at all these hits and find the associated tracks
    for (auto const& clusterHitPtr : clusterHits) {
      // Get the tracks associated with this hit
      auto const& clusterHitTracks = fmt.at(clusterHitPtr.key());
      if (clusterHitTracks.size() > 1) {
        std::cout << "More than one track associated with this hit!\n";
        continue;
      }

      if (clusterHitTracks.size() < 1) continue;

      auto const& clusterHitTrackPtr = clusterHitTracks[0];
      if (clusterHitTrackPtr->Length() < fMinTrackLength) {
        if (fDebug > 2)
          std::cout << "Track " << clusterHitTrackPtr->ID() << " is too short! ("
                    << clusterHitTrackPtr->Length() << ")\n";
        continue;
      }

      // Add this cluster to the track map
      int track = clusterHitTrackPtr.key();
      int cluster = clusterPtr.key();
      if (cet::search_all(clustersToIgnore, cluster)) continue;
      if (not cet::search_all(trackToClusters[track], cluster))
        trackToClusters[track].push_back(cluster);
      if (not cet::search_all(clusterToTracks[cluster], track))
        clusterToTracks[cluster].push_back(track);
    }
  }
}

void shower::EMShowerAlg::CheckIsolatedHits_(
  std::map<int, std::vector<art::Ptr<recob::Hit>>>& showerHitsMap) const
{
  std::map<int, std::vector<int>> firstTPC;
  for (auto const& [plane, hits] : showerHitsMap)
    firstTPC[hits.at(0)->WireID().TPC].push_back(plane);

  // If all in the same TPC then that's great!
  if (firstTPC.size() != 2) return;

  // If they are in more than two TPCs, not much we can do
  else if (firstTPC.size() > 2)
    return;

  // If we get to this point, there should be something we can do!

  // Find the problem plane
  int problemPlane = -1;
  for (auto const& planes : firstTPC | views::values)
    if (planes.size() == 1) problemPlane = planes[0];

  // Require three hits
  if (showerHitsMap.at(problemPlane).size() < 3) return;

  // and get the other planes with at least three hits
  std::vector<int> otherPlanes;
  for (int plane = 0; plane < (int)fWireReadoutGeom->MaxPlanes(); ++plane)
    if (plane != problemPlane and showerHitsMap.count(plane) and
        showerHitsMap.at(plane).size() >= 3)
      otherPlanes.push_back(plane);

  if (otherPlanes.size() == 0) return;

  // Look at the hits after the first one
  unsigned int wrongTPC = showerHitsMap.at(problemPlane).at(0)->WireID().TPC;
  unsigned int nHits = 0;
  for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt = showerHitsMap.at(problemPlane).begin();
       hitIt != showerHitsMap.at(problemPlane).end() and (*hitIt)->WireID().TPC == wrongTPC;
       ++hitIt)
    ++nHits;

  // If there are more than two hits in the 'wrong TPC', we can't be sure it is indeed wrong
  if (nHits > 2) return;

  // See if at least the next four times as many hits are in a different TPC
  std::map<int, int> otherTPCs;
  for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt =
         std::next(showerHitsMap.at(problemPlane).begin(), nHits);
       hitIt != showerHitsMap.at(problemPlane).end() and
       std::distance(std::next(showerHitsMap.at(problemPlane).begin(), nHits), hitIt) < 4 * nHits;
       ++hitIt)
    ++otherTPCs[(*hitIt)->WireID().TPC];

  if (otherTPCs.size() > 1) return;

  // If we get this far, we can move the problem hits from the front of the shower to the back
  std::map<int, int> tpcCount;
  for (int const otherPlane : otherPlanes)
    for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt =
           std::next(showerHitsMap.at(otherPlane).begin());
         hitIt != showerHitsMap.at(otherPlane).end() and
         hitIt != std::next(showerHitsMap.at(otherPlane).begin(), 2);
         ++hitIt)
      ++tpcCount[(*hitIt)->WireID().TPC];

  // Remove the first hit if it is in the wrong TPC
  if (tpcCount.size() == 1 and
      tpcCount.begin()->first ==
        (int)(*std::next(showerHitsMap.at(problemPlane).begin(), nHits))->WireID().TPC) {
    std::vector<art::Ptr<recob::Hit>> naughty_hits;
    for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt = showerHitsMap.at(problemPlane).begin();
         hitIt != std::next(showerHitsMap.at(problemPlane).begin(), nHits);
         ++hitIt) {
      naughty_hits.push_back(*hitIt);
      showerHitsMap.at(problemPlane).erase(hitIt);
    }
    for (auto const& naughty_hit : naughty_hits)
      showerHitsMap.at(problemPlane).push_back(naughty_hit);
  }
}

bool shower::EMShowerAlg::CheckShowerHits_(
  detinfo::DetectorPropertiesData const& detProp,
  std::map<int, std::vector<art::Ptr<recob::Hit>>> const& showerHitsMap) const
{
  bool consistencyCheck = true;

  if (showerHitsMap.size() < 2) { consistencyCheck = true; }
  else if (showerHitsMap.size() == 2) {

    // With two views, we can check:
    //  -- timing between views is consistent
    //  -- the 3D start point makes sense when projected back onto the individual planes

    std::vector<art::Ptr<recob::Hit>> startHits;
    std::vector<int> planes;
    for (auto const& [plane, hits] : showerHitsMap) {
      startHits.push_back(hits.front());
      planes.push_back(plane);
    }

    auto const showerStartPos = Construct3DPoint_(detProp, startHits.at(0), startHits.at(1));
    TVector2 proj1 = Project3DPointOntoPlane_(detProp, showerStartPos, planes.at(0));
    TVector2 proj2 = Project3DPointOntoPlane_(detProp, showerStartPos, planes.at(1));

    double timingDifference = std::abs(startHits.at(0)->PeakTime() - startHits.at(1)->PeakTime());
    double projectionDifference = ((HitPosition_(detProp, *startHits.at(0)) - proj1).Mod() +
                                   (HitPosition_(detProp, *startHits.at(1)) - proj2).Mod()) /
                                  (double)2;

    if (timingDifference > 40 or projectionDifference > 5 or showerStartPos.X() == -9999 or
        showerStartPos.Y() == -9999 or showerStartPos.Z() == -9999)
      consistencyCheck = false;

    if (fDebug > 1)
      std::cout << "Timing difference is " << timingDifference << " and projection distance is "
                << projectionDifference << " (start is (" << showerStartPos.X() << ", "
                << showerStartPos.Y() << ", " << showerStartPos.Z() << ")\n";
  }
  else if (showerHitsMap.size() == 3) {

    // With three views, we can check:
    //  -- the timing between views is consistent
    //  -- the 3D start point formed by two views and projected back into the third is close to the start point in that view

    std::map<int, art::Ptr<recob::Hit>> start2DMap;
    for (auto const& [plane, hits] : showerHitsMap) {
      start2DMap[plane] = hits.front();
    }

    std::map<int, double> projDiff;
    std::map<int, double> timingDiff;

    for (int plane = 0; plane < 3; ++plane) {

      std::vector<int> otherPlanes;
      for (int otherPlane = 0; otherPlane < 3; ++otherPlane)
        if (otherPlane != plane) otherPlanes.push_back(otherPlane);

      auto const showerStartPos = Construct3DPoint_(
        detProp, start2DMap.at(otherPlanes.at(0)), start2DMap.at(otherPlanes.at(1)));
      TVector2 showerStartProj = Project3DPointOntoPlane_(detProp, showerStartPos, plane);

      if (fDebug > 2) {
        std::cout << "Plane... " << plane;
        std::cout << "\nStart position in this plane is "
                  << HitPosition_(detProp, *start2DMap.at(plane)).X() << ", "
                  << HitPosition_(detProp, *start2DMap.at(plane)).Y() << ")\n";
        std::cout << "Shower start from other two planes is (" << showerStartPos.X() << ", "
                  << showerStartPos.Y() << ", " << showerStartPos.Z() << ")\n";
        std::cout << "Projecting the other two planes gives position (" << showerStartProj.X()
                  << ", " << showerStartProj.Y() << ")\n";
      }

      double projDiff =
        std::abs((showerStartProj - HitPosition_(detProp, *start2DMap.at(plane))).Mod());
      double timeDiff = TMath::Max(
        std::abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(0))->PeakTime()),
        std::abs(start2DMap.at(plane)->PeakTime() - start2DMap.at(otherPlanes.at(1))->PeakTime()));

      if (fDebug > 1)
        std::cout << "Plane " << plane << " has projDiff " << projDiff << " and timeDiff "
                  << timeDiff << '\n';

      if (projDiff > 5 or timeDiff > 40) consistencyCheck = false;
    }
  }

  if (fDebug > 1) std::cout << "Consistency check is " << consistencyCheck << '\n';

  return consistencyCheck;
}

std::vector<int> shower::EMShowerAlg::CheckShowerPlanes(
  std::vector<std::vector<int>> const& initialShowers,
  std::vector<art::Ptr<recob::Cluster>> const& clusters,
  art::FindManyP<recob::Hit> const& fmh) const
{
  std::vector<int> clustersToIgnore;

  // Look at each shower
  for (auto initialShowerIt = initialShowers.cbegin(); initialShowerIt != initialShowers.cend();
       ++initialShowerIt) {

    if (std::distance(initialShowers.cbegin(), initialShowerIt) > 0) continue;

    // Map the clusters and cluster hits to each view
    std::map<int, std::vector<art::Ptr<recob::Cluster>>> planeClusters;
    std::map<int, std::vector<art::Ptr<recob::Hit>>> planeHits;
    for (int const clusterIndex : *initialShowerIt) {
      art::Ptr<recob::Cluster> const& cluster = clusters.at(clusterIndex);
      planeClusters[cluster->Plane().Plane].push_back(cluster);
      for (auto const& hit : fmh.at(cluster.key()))
        planeHits[hit->WireID().Plane].push_back(hit);
    }

    TFile* outFile = new TFile("chargeDistributions.root", "RECREATE");
    std::map<int, TH1D*> chargeDist;
    for (auto const& [plane, clusterPtrs] : planeClusters) {
      for (auto const& clusterPtr : clusterPtrs) {
        chargeDist[plane] = new TH1D(std::string("chargeDist_Plane" + std::to_string(plane) +
                                                 "_Cluster" + std::to_string(clusterPtr.key()))
                                       .c_str(),
                                     "",
                                     150,
                                     0,
                                     1000);
        auto const& hits = fmh.at(clusterPtr.key());
        for (auto const& hit : hits | views::transform(to_element)) {
          chargeDist[plane]->Fill(hit.Integral());
        }
        outFile->cd();
        chargeDist[plane]->Write();
      }
    }
    outFile->Close();
    delete outFile;

    // Can't do much with fewer than three views
    if (planeClusters.size() < 3) continue;

    // Look at how many clusters each plane has, and the proportion of hits each one uses
    std::map<int, std::vector<double>> planeClusterSizes;
    for (std::map<int, std::vector<art::Ptr<recob::Cluster>>>::iterator planeClustersIt =
           planeClusters.begin();
         planeClustersIt != planeClusters.end();
         ++planeClustersIt) {
      for (std::vector<art::Ptr<recob::Cluster>>::iterator planeClusterIt =
             planeClustersIt->second.begin();
           planeClusterIt != planeClustersIt->second.end();
           ++planeClusterIt) {
        std::vector<art::Ptr<recob::Hit>> hits = fmh.at(planeClusterIt->key());
        planeClusterSizes[planeClustersIt->first].push_back(
          (double)hits.size() / (double)planeHits.at(planeClustersIt->first).size());
      }
    }

    // Find the average hit fraction across all clusters in the plane
    std::map<int, double> planeClustersAvSizes;
    for (auto const& [plane, cluster_sizes] : planeClusterSizes) {
      double const average = accumulate(cluster_sizes, 0.) / cluster_sizes.size();
      planeClustersAvSizes[plane] = average;
    }

    // Now decide if there is one plane which is ruining the reconstruction
    // If two planes have a low average cluster fraction and one high, this plane likely merges two particle deposits together
    int badPlane = -1;
    for (auto const [plane, avg] : planeClustersAvSizes) {

      // Get averages from other planes and add in quadrature
      std::vector<double> otherAverages;
      for (auto const [other_plane, other_avg] : planeClustersAvSizes)
        if (plane != other_plane) otherAverages.push_back(other_avg);

      double const sumSquareAvgsOtherPlanes = accumulate(
        otherAverages, 0., [](double sum, double avg) { return sum + cet::square(avg); });
      double const quadOtherPlanes = std::sqrt(sumSquareAvgsOtherPlanes);

      // If this plane has an average higher than the quadratic sum of the
      // others, it may be bad
      if (avg >= quadOtherPlanes) badPlane = plane;
    }

    if (badPlane != -1) {
      if (fDebug > 1) std::cout << "Bad plane is " << badPlane << '\n';
      for (auto const& cluster : planeClusters.at(badPlane))
        clustersToIgnore.push_back(cluster.key());
    }
  }

  return clustersToIgnore;
}

geo::Point_t shower::EMShowerAlg::Construct3DPoint_(detinfo::DetectorPropertiesData const& detProp,
                                                    art::Ptr<recob::Hit> const& hit1,
                                                    art::Ptr<recob::Hit> const& hit2) const
{
  // x is average of the two x's
  double x = (detProp.ConvertTicksToX(hit1->PeakTime(), hit1->WireID().planeID()) +
              detProp.ConvertTicksToX(hit2->PeakTime(), hit2->WireID().planeID())) /
             (double)2;

  // y and z got from the wire interections
  auto intersection = fWireReadoutGeom->WireIDsIntersect(hit1->WireID(), hit2->WireID())
                        .value_or(geo::WireIDIntersection::invalid());
  return {x, intersection.y, intersection.z};
}

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> const& hits1,
  std::vector<art::Ptr<recob::Hit>> const& hits2,
  std::map<int, TVector2> const& showerCentreMap) const
{
  std::unique_ptr<recob::Track> track;

  std::vector<art::Ptr<recob::Hit>> track1, track2;

  // Check the TPCs
  if ((*hits1.begin())->WireID().TPC != (*hits2.begin())->WireID().TPC) {
    mf::LogWarning("EMShowerAlg") << "Warning: attempting to construct a track from two different "
                                     "TPCs.  Returning a null track.";
    return track;
  }

  // Check for tracks crossing TPC boundaries
  std::map<int, int> tpcMap;
  for (auto const& hit : hits1)
    ++tpcMap[hit->WireID().TPC];
  for (auto const& hit : hits2)
    ++tpcMap[hit->WireID().TPC];
  if (tpcMap.size() > 1) {
    mf::LogWarning("EMShowerAlg")
      << "Warning: attempting to construct a track which crosses more than one TPC -- PMTrack "
         "can't handle this right now.  Returning a track made just from hits in the first TPC it "
         "traverses.";
    unsigned int firstTPC1 = hits1.at(0)->WireID().TPC, firstTPC2 = hits2.at(0)->WireID().TPC;
    for (auto const& hit : hits1)
      if (hit->WireID().TPC == firstTPC1) track1.push_back(hit);
    for (auto const& hit : hits2)
      if (hit->WireID().TPC == firstTPC2) track2.push_back(hit);
  }
  else {
    track1 = hits1;
    track2 = hits2;
  }

  if (fDebug > 1) {
    std::cout << "About to make a track from these hits:\n";
    auto print_hits = [this](auto const& track) {
      for (auto const& hit : track | views::transform(to_element)) {
        std::cout << "Hit (" << HitCoordinates_(hit).X() << ", " << HitCoordinates_(hit).Y()
                  << ") (real wire " << hit.WireID().Wire << ") in TPC " << hit.WireID().TPC
                  << '\n';
      }
    };
    print_hits(track1);
    print_hits(track2);
  }

  auto const trackStart = Construct3DPoint_(detProp, track1.at(0), track2.at(0));
  pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(detProp, track1, track2, trackStart);

  if (!pmatrack) {
    mf::LogInfo("EMShowerAlg") << "Skipping this event because not enough hits in two views";
    return track;
  }

  std::vector<TVector3> xyz, dircos;

  for (unsigned int i = 0; i < pmatrack->size(); i++) {

    xyz.push_back((*pmatrack)[i]->Point3D());

    if (i < pmatrack->size() - 1) {
      size_t j = i + 1;
      double mag = 0.0;
      TVector3 dc(0., 0., 0.);
      while ((mag == 0.0) and (j < pmatrack->size())) {
        dc = (*pmatrack)[j]->Point3D();
        dc -= (*pmatrack)[i]->Point3D();
        mag = dc.Mag();
        ++j;
      }
      if (mag > 0.0)
        dc *= 1.0 / mag;
      else if (!dircos.empty())
        dc = dircos.back();
      dircos.push_back(dc);
    }
    else
      dircos.push_back(dircos.back());
  }

  // Orient the track correctly
  std::map<int, double> distanceToVertex, distanceToEnd;
  using geo::vect::toPoint;
  geo::Point_t const vertex = toPoint(*xyz.begin());
  geo::Point_t const end = toPoint(*xyz.rbegin());

  // Loop over all the planes and find the distance from the vertex and end
  // projections to the centre in each plane
  for (std::map<int, TVector2>::const_iterator showerCentreIt = showerCentreMap.begin();
       showerCentreIt != showerCentreMap.end();
       ++showerCentreIt) {

    // Project the vertex and the end point onto this plane
    TVector2 vertexProj = Project3DPointOntoPlane_(detProp, vertex, showerCentreIt->first);
    TVector2 endProj = Project3DPointOntoPlane_(detProp, end, showerCentreIt->first);

    // Find the distance of each to the centre of the cluster
    distanceToVertex[showerCentreIt->first] = (vertexProj - showerCentreIt->second).Mod();
    distanceToEnd[showerCentreIt->first] = (endProj - showerCentreIt->second).Mod();
  }

  // Find the average distance to the vertex and the end across the planes
  double avDistanceToVertex = 0, avDistanceToEnd = 0;
  for (std::map<int, double>::iterator distanceToVertexIt = distanceToVertex.begin();
       distanceToVertexIt != distanceToVertex.end();
       ++distanceToVertexIt)
    avDistanceToVertex += distanceToVertexIt->second;
  avDistanceToVertex /= distanceToVertex.size();

  for (std::map<int, double>::iterator distanceToEndIt = distanceToEnd.begin();
       distanceToEndIt != distanceToEnd.end();
       ++distanceToEndIt)
    avDistanceToEnd += distanceToEndIt->second;
  if (distanceToEnd.size() != 0) avDistanceToEnd /= distanceToEnd.size();

  if (fDebug > 2)
    std::cout << "Distance to vertex is " << avDistanceToVertex << " and distance to end is "
              << avDistanceToEnd << '\n';

  // Change order if necessary
  if (avDistanceToEnd > avDistanceToVertex) {
    std::reverse(xyz.begin(), xyz.end());
    std::transform(
      dircos.begin(), dircos.end(), dircos.begin(), [](TVector3 const& vec) { return -1 * vec; });
  }

  if (xyz.size() != dircos.size())
    mf::LogError("EMShowerAlg") << "Problem converting pma::Track3D to recob::Track";

  track = std::make_unique<recob::Track>(
    recob::TrackTrajectory(recob::tracking::convertCollToPoint(xyz),
                           recob::tracking::convertCollToVector(dircos),
                           recob::Track::Flags_t(xyz.size()),
                           false),
    0,
    -1.,
    0,
    recob::tracking::SMatrixSym55(),
    recob::tracking::SMatrixSym55(),
    -1);

  return track;
}

std::unique_ptr<recob::Track> shower::EMShowerAlg::ConstructTrack(
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> const& track1,
  std::vector<art::Ptr<recob::Hit>> const& track2) const
{
  std::map<int, TVector2> showerCentreMap;
  return ConstructTrack(detProp, track1, track2, showerCentreMap);
}

double shower::EMShowerAlg::FinddEdx_(detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      std::vector<art::Ptr<recob::Hit>> const& trackHits,
                                      std::unique_ptr<recob::Track> const& track) const
{
  assert(not empty(trackHits));
  if (!track) return -999;

  recob::Hit const& firstHit = *trackHits.front();

  // Get the pitch
  double pitch = 0;
  try {
    pitch = lar::util::TrackPitchInView(*track, firstHit.View());
  }
  catch (...) {
  }

  // Deal with large pitches
  if (pitch > fdEdxTrackLength) {
    return fCalorimetryAlg.dEdx_AREA(clockData, detProp, firstHit, pitch);
  }

  double totalCharge = 0, totalDistance = 0, avHitTime = 0;
  unsigned int nHits = 0;

  for (auto const& hit : trackHits | views::transform(to_element)) {
    if (totalDistance + pitch < fdEdxTrackLength) {
      totalDistance += pitch;
      totalCharge += hit.Integral();
      avHitTime += hit.PeakTime();
      ++nHits;
    }
  }

  avHitTime /= (double)nHits;

  double const dQdx = totalCharge / totalDistance;
  return fCalorimetryAlg.dEdx_AREA(clockData, detProp, dQdx, avHitTime, firstHit.WireID().Plane);
}

void shower::EMShowerAlg::FindInitialTrack(
  detinfo::DetectorPropertiesData const& detProp,
  const std::map<int, std::vector<art::Ptr<recob::Hit>>>& showerHitsMap,
  std::unique_ptr<recob::Track>& initialTrack,
  std::map<int, std::vector<art::Ptr<recob::Hit>>>& initialTrackHits) const
{

  /// Finding the initial track requires three stages:
  ///  -- find the initial track-like hits in each view
  ///  -- use these to construct a track

  // Now find the hits belonging to the track
  if (fDebug > 1)
    std::cout << " ------------------ Finding initial track hits "
                 "-------------------- \n";
  initialTrackHits = FindShowerStart_(showerHitsMap);
  if (fDebug > 1) {
    std::cout << "Here are the initial shower hits... \n";
    for (auto const& [plane, hitPtrs] : initialTrackHits) {
      std::cout << "  Plane " << plane << '\n';
      for (auto const& hit : hitPtrs | views::transform(to_element)) {
        std::cout << "    Hit is (" << HitCoordinates_(hit).X() << " (real hit " << hit.WireID()
                  << "), " << HitCoordinates_(hit).Y() << ")\n";
      }
    }
  }
  if (fDebug > 1)
    std::cout << " ------------------ End finding initial track hits "
                 "-------------------- \n";

  // Now we have the track hits -- can make a track!
  if (fDebug > 1) std::cout << " ------------------ Finding initial track -------------------- \n";
  initialTrack = MakeInitialTrack_(detProp, initialTrackHits, showerHitsMap);
  if (initialTrack and fDebug > 1) {
    std::cout << "The track start is (" << initialTrack->Vertex().X() << ", "
              << initialTrack->Vertex().Y() << ", " << initialTrack->Vertex().Z() << ")\n";
    std::cout << "The track direction is (" << initialTrack->VertexDirection().X() << ", "
              << initialTrack->VertexDirection().Y() << ", " << initialTrack->VertexDirection().Z()
              << ")\n";
  }
  if (fDebug > 1)
    std::cout << " ------------------ End finding initial track "
                 "-------------------- \n";
}

std::vector<art::Ptr<recob::Hit>> shower::EMShowerAlg::FindOrderOfHits_(
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> const& hits,
  bool perpendicular) const
{
  // Find the charge-weighted centre (in [cm]) of this shower
  TVector2 centre = ShowerCenter_(detProp, hits);

  // Find a rough shower 'direction'
  TVector2 direction = ShowerDirection_(detProp, hits);

  if (perpendicular) direction = direction.Rotate(TMath::Pi() / 2);

  // Find how far each hit (projected onto this axis) is from the centre
  TVector2 pos;
  std::map<double, art::Ptr<recob::Hit>> hitProjection;
  for (auto const& hitPtr : hits) {
    pos = HitPosition_(detProp, *hitPtr) - centre;
    hitProjection[direction * pos] = hitPtr;
  }

  // Get a vector of hits in order of the shower
  std::vector<art::Ptr<recob::Hit>> showerHits;
  cet::transform_all(
    hitProjection, std::back_inserter(showerHits), [](auto const& pr) { return pr.second; });

  // Make gradient plot
  if (fMakeGradientPlot) {
    std::map<int, TGraph*> graphs;
    for (auto const& hit : showerHits | views::transform(to_element)) {
      int tpc = hit.WireID().TPC;
      if (graphs[tpc] == nullptr) graphs[tpc] = new TGraph();
      graphs[tpc]->SetPoint(
        graphs[tpc]->GetN(), HitPosition_(detProp, hit).X(), HitPosition_(detProp, hit).Y());
    }
    TMultiGraph* multigraph = new TMultiGraph();
    for (auto const [color, graph] : graphs) {
      graph->SetMarkerColor(color);
      graph->SetMarkerStyle(8);
      graph->SetMarkerSize(2);
      multigraph->Add(graph);
    }
    TCanvas* canvas = new TCanvas();
    multigraph->Draw("AP");
    TLine line;
    line.SetLineColor(2);
    line.DrawLine(centre.X() - 1000 * direction.X(),
                  centre.Y() - 1000 * direction.Y(),
                  centre.X() + 1000 * direction.X(),
                  centre.Y() + 1000 * direction.Y());
    canvas->SaveAs("Gradient.png");
    delete canvas;
    delete multigraph;
  }

  return showerHits;
}

std::vector<std::vector<int>> shower::EMShowerAlg::FindShowers(
  std::map<int, std::vector<int>> const& trackToClusters) const
{
  // Showers are vectors of clusters
  std::vector<std::vector<int>> showers;

  // Loop over all tracks
  for (auto const& clusters : trackToClusters | views::values) {

    // Find which showers already made are associated with this track
    std::vector<int> matchingShowers;
    for (unsigned int shower = 0; shower < showers.size(); ++shower)
      for (int const cluster : clusters) {
        if (cet::search_all(showers[shower], cluster) and
            not cet::search_all(matchingShowers, shower)) {
          matchingShowers.push_back(shower);
        }
      }

    // THINK THERE PROBABLY CAN BE MORE THAN ONE!
    // IN FACT, THIS WOULD BE A SUCCESS OF THE SHOWERING METHOD!
    // // Shouldn't be more than one
    // if (matchingShowers.size() > 1)
    //   mf::LogInfo("EMShowerAlg") << "Number of showers this track matches is " << matchingShowers.size() << std::endl;

    // New shower
    if (matchingShowers.size() < 1) showers.push_back(clusters);

    // Add to existing shower
    else {
      for (int const cluster : clusters) {
        if (not cet::search_all(showers.at(matchingShowers[0]), cluster))
          showers.at(matchingShowers.at(0)).push_back(cluster);
      }
    }
  }

  return showers;
}

std::map<int, std::vector<art::Ptr<recob::Hit>>> shower::EMShowerAlg::FindShowerStart_(
  std::map<int, std::vector<art::Ptr<recob::Hit>>> const& orderedShowerMap) const
{

  std::map<int, std::vector<art::Ptr<recob::Hit>>> initialHitsMap;

  for (auto const& [plane, orderedShower] : orderedShowerMap) {
    std::vector<art::Ptr<recob::Hit>> initialHits;

    // Find if the shower is traveling along ticks or wires
    bool wireDirection = true;
    std::vector<int> wires;
    for (auto const& hit : orderedShower | views::transform(to_element))
      wires.push_back(std::round(HitCoordinates_(hit).X()));

    cet::sort_all(wires);
    if (std::abs(*wires.begin() - std::round(HitCoordinates_(**orderedShower.begin()).X())) > 3 and
        std::abs(*wires.rbegin() - std::round(HitCoordinates_(**orderedShower.begin()).X())) > 3)
      wireDirection = false;

    // Deal with showers traveling along wires
    if (wireDirection) {
      bool increasing = HitCoordinates_(**orderedShower.rbegin()).X() >
                        HitCoordinates_(**orderedShower.begin()).X();
      std::map<int, std::vector<art::Ptr<recob::Hit>>> wireHitMap;
      int multipleWires = 0;
      for (auto const& hitPtr : orderedShower)
        wireHitMap[std::round(HitCoordinates_(*hitPtr).X())].push_back(hitPtr);

      for (auto const& hitPtr : orderedShower) {
        int wire = std::round(HitCoordinates_(*hitPtr).X());
        if (wireHitMap[wire].size() > 1) {
          ++multipleWires;
          if (multipleWires > 5) break;
          continue;
        }
        else if ((increasing and wireHitMap[wire + 1].size() > 1 and
                  (wireHitMap[wire + 2].size() > 1 or wireHitMap[wire + 3].size() > 1)) or
                 (!increasing and wireHitMap[wire - 1].size() > 1 and
                  (wireHitMap[wire - 2].size() > 1 or wireHitMap[wire - 3].size() > 1))) {
          if ((increasing and
               (wireHitMap[wire + 4].size() < 2 and wireHitMap[wire + 5].size() < 2 and
                wireHitMap[wire + 6].size() < 2 and wireHitMap[wire + 9].size() > 1)) or
              (!increasing and
               (wireHitMap[wire - 4].size() < 2 and wireHitMap[wire - 5].size() < 2 and
                wireHitMap[wire - 6].size() < 2) and
               wireHitMap[wire - 9].size() > 1))
            initialHits.push_back(hitPtr);
          else
            break;
        }
        else
          initialHits.push_back(hitPtr);
      }
      if (!initialHits.size()) initialHits.push_back(*orderedShower.begin());
    }

    // Deal with showers travelling along ticks
    else {
      bool increasing = HitCoordinates_(**orderedShower.rbegin()).Y() >
                        HitCoordinates_(**orderedShower.begin()).Y();
      std::map<int, std::vector<art::Ptr<recob::Hit>>> tickHitMap;
      for (std::vector<art::Ptr<recob::Hit>>::const_iterator hitIt = orderedShower.begin();
           hitIt != orderedShower.end();
           ++hitIt)
        tickHitMap[std::round(HitCoordinates_(**hitIt).Y())].push_back(*hitIt);

      for (auto const& hitPtr : orderedShower) {
        int const tick = std::round(HitCoordinates_(*hitPtr).Y());
        if ((increasing and (tickHitMap[tick + 1].size() or tickHitMap[tick + 2].size() or
                             tickHitMap[tick + 3].size() or tickHitMap[tick + 4].size() or
                             tickHitMap[tick + 5].size())) or
            (!increasing and (tickHitMap[tick - 1].size() or tickHitMap[tick - 2].size() or
                              tickHitMap[tick - 3].size() or tickHitMap[tick - 4].size() or
                              tickHitMap[tick - 5].size())))
          break;
        else
          initialHits.push_back(hitPtr);
      }
      if (initialHits.empty()) initialHits.push_back(*orderedShower.begin());
    }

    // Need at least two hits to make a track
    if (initialHits.size() == 1 and orderedShower.size() > 2)
      initialHits.push_back(orderedShower[1]);

    // Quality check -- make sure there isn't a single hit in a different TPC (artefact of tracking failure)
    std::vector<art::Ptr<recob::Hit>> newInitialHits;
    std::map<int, int> tpcHitMap;
    std::vector<int> singleHitTPCs;
    for (auto const& hit : initialHits | views::transform(to_element))
      ++tpcHitMap[hit.WireID().TPC];

    for (auto const [tpc, count] : tpcHitMap)
      if (count == 1) singleHitTPCs.push_back(tpc);

    if (singleHitTPCs.size()) {
      if (fDebug > 2)
        for (int const tpc : singleHitTPCs)
          std::cout << "Removed hits in TPC " << tpc << '\n';

      for (auto const& hitPtr : initialHits)
        if (not cet::search_all(singleHitTPCs, hitPtr->WireID().TPC))
          newInitialHits.push_back(hitPtr);
      if (!newInitialHits.size()) newInitialHits.push_back(*orderedShower.begin());
    }
    else
      newInitialHits = initialHits;

    initialHitsMap[plane] = newInitialHits;
  }

  return initialHitsMap;
}

std::map<int, std::map<int, bool>> shower::EMShowerAlg::GetPlanePermutations_(
  const detinfo::DetectorPropertiesData& detProp,
  const std::map<int, std::vector<art::Ptr<recob::Hit>>>& showerHitsMap) const
{

  // The map to return
  std::map<int, std::map<int, bool>> permutations;

  // Get the properties of the shower hits across the planes which will be used to
  // determine the likelihood of a particular reorientation permutation
  //   -- relative width in the wire direction (if showers travel along the wire
  //      direction in a particular plane)
  //   -- the RMS gradients (how likely it is the RMS of the hit positions from
  //      central axis increases along a particular orientation)

  // Find the RMS, RMS gradient and wire widths
  std::map<int, double> planeRMSGradients, planeRMS;
  for (auto const& [plane, hitPtrs] : showerHitsMap) {
    planeRMS[plane] = ShowerHitRMS_(detProp, hitPtrs);
    planeRMSGradients[plane] = ShowerHitRMSGradient_(detProp, hitPtrs);
  }

  // Order these backwards so they can be used to discriminate between planes
  std::map<double, int> gradientMap;
  for (int const plane : showerHitsMap | views::keys)
    gradientMap[std::abs(planeRMSGradients.at(plane))] = plane;

  std::map<double, int> wireWidthMap = RelativeWireWidth_(showerHitsMap);

  if (fDebug > 1)
    for (auto const [gradient, plane] : wireWidthMap)
      std::cout << "Plane " << plane << " has relative width in wire of " << gradient
                << "; and an RMS gradient of " << planeRMSGradients[plane] << '\n';

  // Find the correct permutations
  int perm = 0;
  std::vector<std::map<int, bool>> usedPermutations;

  // Most likely is to not change anything
  for (int const plane : showerHitsMap | views::keys)
    permutations[perm][plane] = 0;
  ++perm;

  // Use properties of the shower to determine the middle cases
  for (int const plane : wireWidthMap | views::values) {
    std::map<int, bool> permutation;
    permutation[plane] = true;
    for (int const plane2 : wireWidthMap | views::values)
      if (plane != plane2) permutation[plane2] = false;

    if (not cet::search_all(usedPermutations, permutation)) {
      permutations[perm] = permutation;
      usedPermutations.push_back(permutation);
      ++perm;
    }
  }
  for (int const plane : wireWidthMap | views::reverse | views::values) {
    std::map<int, bool> permutation;
    permutation[plane] = false;
    for (int const plane2 : wireWidthMap | views::values)
      if (plane != plane2) permutation[plane2] = true;

    if (not cet::search_all(usedPermutations, permutation)) {
      permutations[perm] = permutation;
      usedPermutations.push_back(permutation);
      ++perm;
    }
  }

  // Least likely is to change everything
  for (int const plane : showerHitsMap | views::keys)
    permutations[perm][plane] = 1;
  ++perm;

  if (fDebug > 1) {
    std::cout << "Here are the permutations!\n";
    for (auto const& [index, permutation] : permutations) {
      std::cout << "Permutation " << index << '\n';
      for (auto const [plane, value] : permutation)
        std::cout << "  Plane " << plane << " has value " << value << '\n';
    }
  }

  return permutations;
}

std::unique_ptr<recob::Track> shower::EMShowerAlg::MakeInitialTrack_(
  detinfo::DetectorPropertiesData const& detProp,
  std::map<int, std::vector<art::Ptr<recob::Hit>>> const& initialHitsMap,
  std::map<int, std::vector<art::Ptr<recob::Hit>>> const& showerHitsMap) const
{
  // Can't do much with just one view
  if (initialHitsMap.size() < 2) {
    mf::LogInfo("EMShowerAlg") << "Only one useful view for this shower; nothing can be done\n";
    return std::unique_ptr<recob::Track>();
  }

  std::vector<std::pair<int, int>> initialHitsSize;
  for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::const_iterator initialHitIt =
         initialHitsMap.begin();
       initialHitIt != initialHitsMap.end();
       ++initialHitIt)
    initialHitsSize.push_back(std::make_pair(initialHitIt->first, initialHitIt->second.size()));

  // Sort the planes by number of hits
  std::sort(initialHitsSize.begin(),
            initialHitsSize.end(),
            [](std::pair<int, int> const& size1, std::pair<int, int> const& size2) {
              return size1.second > size2.second;
            });

  // Pick the two planes to use to make the track
  //   -- if more than two planes, can choose the two views
  //      more accurately
  //   -- if not, just use the two views available

  std::vector<int> planes;

  if (initialHitsSize.size() > 2 and !CheckShowerHits_(detProp, showerHitsMap)) {
    int planeToIgnore = WorstPlane_(showerHitsMap);
    if (fDebug > 1)
      std::cout << "Igoring plane " << planeToIgnore << " in creation of initial track\n";
    for (std::vector<std::pair<int, int>>::iterator hitsSizeIt = initialHitsSize.begin();
         hitsSizeIt != initialHitsSize.end() and planes.size() < 2;
         ++hitsSizeIt) {
      if (hitsSizeIt->first == planeToIgnore) continue;
      planes.push_back(hitsSizeIt->first);
    }
  }
  else
    planes = {initialHitsSize.at(0).first, initialHitsSize.at(1).first};

  return ConstructTrack(detProp, initialHitsMap.at(planes.at(0)), initialHitsMap.at(planes.at(1)));
}

recob::Shower shower::EMShowerAlg::MakeShower(
  detinfo::DetectorClocksData const& clockData,
  detinfo::DetectorPropertiesData const& detProp,
  art::PtrVector<recob::Hit> const& hits,
  std::unique_ptr<recob::Track> const& initialTrack,
  std::map<int, std::vector<art::Ptr<recob::Hit>>> const& initialHitsMap) const
{

  // Find the shower hits on each plane
  std::map<int, std::vector<art::Ptr<recob::Hit>>> planeHitsMap;
  for (auto const& hitPtr : hits)
    planeHitsMap[hitPtr->View()].push_back(hitPtr);

  int bestPlane = -1;
  unsigned int highestNumberOfHits = 0;
  std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;

  // Look at each of the planes
  for (unsigned int plane = 0; plane < fWireReadoutGeom->MaxPlanes(); ++plane) {

    // If there's hits on this plane...
    if (planeHitsMap.count(plane)) {
      dEdx.push_back(FinddEdx_(clockData, detProp, initialHitsMap.at(plane), initialTrack));
      totalEnergy.push_back(
        fShowerEnergyAlg.ShowerEnergy(clockData, detProp, planeHitsMap.at(plane), plane));
      if (planeHitsMap.at(plane).size() > highestNumberOfHits and initialHitsMap.count(plane)) {
        bestPlane = plane;
        highestNumberOfHits = planeHitsMap.at(plane).size();
      }
    }

    // If not...
    else {
      dEdx.push_back(0);
      totalEnergy.push_back(0);
    }
  }

  TVector3 direction, directionError, showerStart, showerStartError;
  if (initialTrack) {
    direction = initialTrack->VertexDirection<TVector3>();
    showerStart = initialTrack->Vertex<TVector3>();
  }

  if (fDebug > 0) {
    std::cout << "Best plane is " << bestPlane;
    std::cout << "\ndE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and " << dEdx[2];
    std::cout << "\nTotal energy for each plane is: " << totalEnergy[0] << ", " << totalEnergy[1]
              << " and " << totalEnergy[2];
    std::cout << "\nThe shower start is (" << showerStart.X() << ", " << showerStart.Y() << ", "
              << showerStart.Z() << ")\n";
    std::cout << "The shower direction is (" << direction.X() << ", " << direction.Y() << ", "
              << direction.Z() << ")\n";
  }

  return recob::Shower(direction,
                       directionError,
                       showerStart,
                       showerStartError,
                       totalEnergy,
                       totalEnergyError,
                       dEdx,
                       dEdxError,
                       bestPlane);
}

recob::Shower shower::EMShowerAlg::MakeShower(detinfo::DetectorClocksData const& clockData,
                                              detinfo::DetectorPropertiesData const& detProp,
                                              art::PtrVector<recob::Hit> const& hits,
                                              art::Ptr<recob::Vertex> const& vertex,
                                              int& iok) const
{
  iok = 1;

  // Find the shower hits on each plane
  std::map<int, std::vector<art::Ptr<recob::Hit>>> planeHitsMap;
  for (auto const& hitPtr : hits)
    planeHitsMap[hitPtr->WireID().Plane].push_back(hitPtr);

  std::vector<std::vector<art::Ptr<recob::Hit>>> initialTrackHits(3);

  int pl0 = -1;
  int pl1 = -1;
  unsigned maxhits0 = 0;
  unsigned maxhits1 = 0;

  for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::iterator planeHits = planeHitsMap.begin();
       planeHits != planeHitsMap.end();
       ++planeHits) {

    std::vector<art::Ptr<recob::Hit>> showerHits;
    OrderShowerHits_(detProp, planeHits->second, showerHits, vertex);
    FindInitialTrackHits(showerHits, vertex, initialTrackHits[planeHits->first]);
    if ((planeHits->second).size() > maxhits0) {
      if (pl0 != -1) {
        maxhits1 = maxhits0;
        pl1 = pl0;
      }
      pl0 = planeHits->first;
      maxhits0 = (planeHits->second).size();
    }
    else if ((planeHits->second).size() > maxhits1) {
      pl1 = planeHits->first;
      maxhits1 = (planeHits->second).size();
    }
  }
  if (pl0 != -1 && pl1 != -1 && initialTrackHits[pl0].size() >= 2 &&
      initialTrackHits[pl1].size() >= 2 &&
      initialTrackHits[pl0][0]->WireID().TPC == initialTrackHits[pl1][0]->WireID().TPC) {
    double xyz[3];
    vertex->XYZ(xyz);
    TVector3 vtx(xyz);
    pma::Track3D* pmatrack =
      fProjectionMatchingAlg.buildSegment(detProp, initialTrackHits[pl0], initialTrackHits[pl1]);
    std::vector<TVector3> spts;

    for (size_t i = 0; i < pmatrack->size(); ++i) {
      if ((*pmatrack)[i]->IsEnabled()) {
        TVector3 p3d = (*pmatrack)[i]->Point3D();
        spts.push_back(p3d);
      }
    }
    if (spts.size() >= 2) { // at least two space points
      TVector3 shwxyz, shwxyzerr;
      TVector3 shwdir, shwdirerr;
      std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;
      int bestPlane = pl0;
      double minpitch = 1000;
      std::vector<TVector3> dirs;
      if ((spts[0] - vtx).Mag() < (spts.back() - vtx).Mag()) {
        shwxyz = spts[0];
        size_t i = 5;
        if (spts.size() - 1 < 5) i = spts.size() - 1;
        shwdir = spts[i] - spts[0];
        shwdir = shwdir.Unit();
      }
      else {
        shwxyz = spts.back();
        size_t i = 0;
        if (spts.size() > 6) i = spts.size() - 6;
        shwdir = spts[i] - spts[spts.size() - 1];
        shwdir = shwdir.Unit();
      }
      for (unsigned int plane = 0; plane < fWireReadoutGeom->MaxPlanes(); ++plane) {
        if (planeHitsMap.find(plane) != planeHitsMap.end()) {
          totalEnergy.push_back(
            fShowerEnergyAlg.ShowerEnergy(clockData, detProp, planeHitsMap[plane], plane));
        }
        else {
          totalEnergy.push_back(0);
        }
        if (initialTrackHits[plane].size()) {
          double fdEdx = 0;
          double totQ = 0;
          double avgT = 0;
          double pitch = 0;
          double wirepitch =
            fWireReadoutGeom->Plane(initialTrackHits[plane][0]->WireID().planeID()).WirePitch();
          double angleToVert = fWireReadoutGeom->WireAngleToVertical(
                                 fWireReadoutGeom->Plane(geo::PlaneID{0, 0, plane}).View(),
                                 initialTrackHits[plane][0]->WireID().planeID()) -
                               0.5 * TMath::Pi();
          double cosgamma = std::abs(sin(angleToVert) * shwdir.Y() + cos(angleToVert) * shwdir.Z());
          if (cosgamma > 0) pitch = wirepitch / cosgamma;
          if (pitch) {
            if (pitch < minpitch) {
              minpitch = pitch;
              bestPlane = plane;
            }
            int nhits = 0;
            std::vector<float> vQ;
            for (auto const& hit : initialTrackHits[plane]) {
              int w1 = hit->WireID().Wire;
              int w0 = initialTrackHits[plane][0]->WireID().Wire;
              if (std::abs((w1 - w0) * pitch) < fdEdxTrackLength) {
                vQ.push_back(hit->Integral());
                totQ += hit->Integral();
                avgT += hit->PeakTime();
                ++nhits;
              }
            }
            if (totQ) {
              double dQdx = TMath::Median(vQ.size(), &vQ[0]) / pitch;
              fdEdx = fCalorimetryAlg.dEdx_AREA(
                clockData, detProp, dQdx, avgT / nhits, initialTrackHits[plane][0]->WireID().Plane);
            }
          }
          dEdx.push_back(fdEdx);
        }
        else {
          dEdx.push_back(0);
        }
      }
      iok = 0;
      if (fDebug > 1) {
        std::cout << "Best plane is " << bestPlane;
        std::cout << "\ndE/dx for each plane is: " << dEdx[0] << ", " << dEdx[1] << " and "
                  << dEdx[2];
        std::cout << "\nTotal energy for each plane is: " << totalEnergy[0] << ", "
                  << totalEnergy[1] << " and " << totalEnergy[2];
        std::cout << "\nThe shower start is (" << shwxyz.X() << ", " << shwxyz.Y() << ", "
                  << shwxyz.Z() << ")\n";
        shwxyz.Print();
      }

      return recob::Shower(shwdir,
                           shwdirerr,
                           shwxyz,
                           shwxyzerr,
                           totalEnergy,
                           totalEnergyError,
                           dEdx,
                           dEdxError,
                           bestPlane);
    }
  }
  return recob::Shower();
}

std::vector<recob::SpacePoint> shower::EMShowerAlg::MakeSpacePoints(
  detinfo::DetectorPropertiesData const& detProp,
  std::map<int, std::vector<art::Ptr<recob::Hit>>> const& showerHits,
  std::vector<std::vector<art::Ptr<recob::Hit>>>& hitAssns) const
{
  // Space points to return
  std::vector<recob::SpacePoint> spacePoints;

  // Make space points
  // Use the following procedure:
  //  -- Consider hits plane by plane
  //  -- For each hit on the first plane, consider the 3D point made by combining with each hit from the second plane
  //  -- Project this 3D point back into the two planes
  //  -- Determine how close to a the original hits this point lies
  //  -- If close enough, make a 3D space point from this point
  //  -- Discard these used hits in future iterations, along with hits in the
  //     third plane (if exists) close to the projection of the point into this
  //     plane

  // Container to hold used hits
  std::vector<int> usedHits;

  // Look through plane by plane
  for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::const_iterator showerHitIt =
         showerHits.begin();
       showerHitIt != showerHits.end();
       ++showerHitIt) {

    // Find the other planes with hits
    std::vector<int> otherPlanes;
    for (unsigned int otherPlane = 0; otherPlane < fWireReadoutGeom->MaxPlanes(); ++otherPlane)
      if ((int)otherPlane != showerHitIt->first and showerHits.count(otherPlane))
        otherPlanes.push_back(otherPlane);

    // Can't make space points if we only have one view
    if (otherPlanes.size() == 0) return spacePoints;

    // Look at all hits on this plane
    for (std::vector<art::Ptr<recob::Hit>>::const_iterator planeHitIt = showerHitIt->second.begin();
         planeHitIt != showerHitIt->second.end();
         ++planeHitIt) {

      if (std::find(usedHits.begin(), usedHits.end(), planeHitIt->key()) != usedHits.end())
        continue;

      // Make a 3D point with every hit on the second plane
      const std::vector<art::Ptr<recob::Hit>> otherPlaneHits = showerHits.at(otherPlanes.at(0));
      for (std::vector<art::Ptr<recob::Hit>>::const_iterator otherPlaneHitIt =
             otherPlaneHits.begin();
           otherPlaneHitIt != otherPlaneHits.end() and
           std::find(usedHits.begin(), usedHits.end(), planeHitIt->key()) == usedHits.end();
           ++otherPlaneHitIt) {

        if ((*otherPlaneHitIt)->WireID().TPC != (*planeHitIt)->WireID().TPC or
            std::find(usedHits.begin(), usedHits.end(), otherPlaneHitIt->key()) != usedHits.end())
          continue;

        auto const point = Construct3DPoint_(detProp, *planeHitIt, *otherPlaneHitIt);
        std::vector<art::Ptr<recob::Hit>> pointHits;
        bool truePoint = false;

        if (otherPlanes.size() > 1) {

          TVector2 projThirdPlane = Project3DPointOntoPlane_(detProp, point, otherPlanes.at(1));
          const std::vector<art::Ptr<recob::Hit>> otherOtherPlaneHits =
            showerHits.at(otherPlanes.at(1));

          for (std::vector<art::Ptr<recob::Hit>>::const_iterator otherOtherPlaneHitIt =
                 otherOtherPlaneHits.begin();
               otherOtherPlaneHitIt != otherOtherPlaneHits.end() and !truePoint;
               ++otherOtherPlaneHitIt) {

            if ((*otherOtherPlaneHitIt)->WireID().TPC == (*planeHitIt)->WireID().TPC and
                (projThirdPlane - HitPosition_(detProp, **otherOtherPlaneHitIt)).Mod() <
                  fSpacePointSize) {

              truePoint = true;

              // Remove hits used to make the point
              usedHits.push_back(planeHitIt->key());
              usedHits.push_back(otherPlaneHitIt->key());
              usedHits.push_back(otherOtherPlaneHitIt->key());

              pointHits.push_back(*planeHitIt);
              pointHits.push_back(*otherPlaneHitIt);
              pointHits.push_back(*otherOtherPlaneHitIt);
            }
          }
        }

        else if ((Project3DPointOntoPlane_(detProp, point, (*planeHitIt)->WireID().Plane) -
                  HitPosition_(detProp, **planeHitIt))
                     .Mod() < fSpacePointSize and
                 (Project3DPointOntoPlane_(detProp, point, (*otherPlaneHitIt)->WireID().Plane) -
                  HitPosition_(detProp, **otherPlaneHitIt))
                     .Mod() < fSpacePointSize) {

          truePoint = true;

          usedHits.push_back(planeHitIt->key());
          usedHits.push_back(otherPlaneHitIt->key());

          pointHits.push_back(*planeHitIt);
          pointHits.push_back(*otherPlaneHitIt);
        }

        // Make space point
        if (truePoint) {
          double xyz[3] = {point.X(), point.Y(), point.Z()};
          double xyzerr[6] = {fSpacePointSize,
                              fSpacePointSize,
                              fSpacePointSize,
                              fSpacePointSize,
                              fSpacePointSize,
                              fSpacePointSize};
          double chisq = 0.;
          spacePoints.emplace_back(xyz, xyzerr, chisq);
          hitAssns.push_back(pointHits);
        }

      } // end loop over second plane hits

    } // end loop over first plane hits

  } // end loop over planes

  if (fDebug > 0) {
    std::cout << "-------------------- Space points -------------------\n";
    std::cout << "There are " << spacePoints.size() << " space points:\n";
    if (fDebug > 1)
      for (std::vector<recob::SpacePoint>::const_iterator spacePointIt = spacePoints.begin();
           spacePointIt != spacePoints.end();
           ++spacePointIt) {
        const double* xyz = spacePointIt->XYZ();
        std::cout << "  Space point (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")\n";
      }
  }

  return spacePoints;
}

std::map<int, std::vector<art::Ptr<recob::Hit>>> shower::EMShowerAlg::OrderShowerHits(
  detinfo::DetectorPropertiesData const& detProp,
  art::PtrVector<recob::Hit> const& shower,
  int desired_plane) const
{
  /// Ordering the shower hits requires three stages:
  ///  -- putting all the hits in a given plane in some kind of order
  ///  -- use the properties of the hits in all three planes to check this order
  ///  -- orient the hits correctly using properties of the shower

  // ------------- Put hits in order ------------

  // Find the shower hits on each plane
  std::map<int, std::vector<art::Ptr<recob::Hit>>> showerHitsMap;
  for (auto const& hitPtr : shower) {
    showerHitsMap[hitPtr->WireID().Plane].push_back(hitPtr);
  }

  // Order the hits, get the RMS and the RMS gradient for the hits in this plane
  std::map<int, double> planeRMSGradients, planeRMS;
  for (auto const& [plane, hits] : showerHitsMap) {
    if (desired_plane != plane and desired_plane != -1) continue;
    std::vector<art::Ptr<recob::Hit>> orderedHits = FindOrderOfHits_(detProp, hits);
    planeRMS[plane] = ShowerHitRMS_(detProp, orderedHits);
    planeRMSGradients[plane] = ShowerHitRMSGradient_(detProp, orderedHits);
    showerHitsMap[plane] = orderedHits;
  }

  if (fDebug > 1) {
    for (auto const [plane, shower_hit_rms] : planeRMS) {
      std::cout << "Plane " << plane << " has RMS " << shower_hit_rms << " and RMS gradient "
                << planeRMSGradients.at(plane) << '\n';
    }
  }

  if (fDebug > 2) {
    std::cout << "\nHits in order; after ordering, before reversing...\n";
    for (auto const& [plane, hitPtrs] : showerHitsMap) {
      std::cout << "  Plane " << plane << '\n';
      for (auto const& hit : hitPtrs | views::transform(to_element)) {
        std::cout << "    Hit at (" << HitCoordinates_(hit).X() << ", " << HitCoordinates_(hit).Y()
                  << ") -- real wire " << hit.WireID() << ", hit position ("
                  << HitPosition_(detProp, hit).X() << ", " << HitPosition_(detProp, hit).Y()
                  << ")\n";
      }
    }
  }

  // ------------- Check between the views to ensure consistency of ordering -------------

  // Check between the views to make sure there isn't a poorly formed shower in just one view
  // First, determine the average RMS and RMS gradient across the other planes
  std::map<int, double> planeOtherRMS, planeOtherRMSGradients;
  for (std::map<int, double>::iterator planeRMSIt = planeRMS.begin(); planeRMSIt != planeRMS.end();
       ++planeRMSIt) {
    planeOtherRMS[planeRMSIt->first] = 0;
    planeOtherRMSGradients[planeRMSIt->first] = 0;
    int nOtherPlanes = 0;
    for (int plane = 0; plane < (int)fWireReadoutGeom->MaxPlanes(); ++plane) {
      if (plane != planeRMSIt->first and planeRMS.count(plane)) {
        planeOtherRMS[planeRMSIt->first] += planeRMS.at(plane);
        planeOtherRMSGradients[planeRMSIt->first] += planeRMSGradients.at(plane);
        ++nOtherPlanes;
      }
    }
    planeOtherRMS[planeRMSIt->first] /= (double)nOtherPlanes;
    planeOtherRMSGradients[planeRMSIt->first] /= (double)nOtherPlanes;
  }

  // Look to see if one plane has a particularly high RMS (compared to the
  // others) whilst having a similar gradient
  for (auto const& [plane, hitPtrs] : showerHitsMap) {
    if (planeRMS.at(plane) > planeOtherRMS.at(plane) * 2) {
      if (fDebug > 1) std::cout << "Plane " << plane << " was perpendicular... recalculating\n";
      std::vector<art::Ptr<recob::Hit>> orderedHits =
        this->FindOrderOfHits_(detProp, hitPtrs, true);
      showerHitsMap[plane] = orderedHits;
      planeRMSGradients[plane] = this->ShowerHitRMSGradient_(detProp, orderedHits);
    }
  }

  // ------------- Orient the shower correctly ---------------

  if (fDebug > 1) {
    std::cout << "Before reversing, here are the start and end points...\n";
    for (auto const& [plane, hitPtrs] : showerHitsMap) {
      std::cout << "  Plane " << plane << " has start (" << HitCoordinates_(*hitPtrs.front()).X()
                << ", " << HitCoordinates_(*hitPtrs.front()).Y() << ") (real wire "
                << hitPtrs.front()->WireID() << ") and end ("
                << HitCoordinates_(*hitPtrs.back()).X() << ", "
                << HitCoordinates_(*hitPtrs.back()).Y() << ") (real wire "
                << hitPtrs.back()->WireID() << ")\n";
    }
  }

  // Use the RMS gradient information to get an initial ordering
  for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::iterator showerHitsIt =
         showerHitsMap.begin();
       showerHitsIt != showerHitsMap.end();
       ++showerHitsIt) {
    double gradient = planeRMSGradients.at(showerHitsIt->first);
    if (gradient < 0) std::reverse(showerHitsIt->second.begin(), showerHitsIt->second.end());
  }

  if (fDebug > 2) {
    std::cout << "\nHits in order; after reversing, before checking isolated hits...\n";
    for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::iterator showerHitsIt =
           showerHitsMap.begin();
         showerHitsIt != showerHitsMap.end();
         ++showerHitsIt) {
      std::cout << "  Plane " << showerHitsIt->first << '\n';
      for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt = showerHitsIt->second.begin();
           hitIt != showerHitsIt->second.end();
           ++hitIt)
        std::cout << "    Hit at (" << HitCoordinates_(**hitIt).X() << ", "
                  << HitCoordinates_(**hitIt).Y() << ") -- real wire " << (*hitIt)->WireID()
                  << ", hit position (" << HitPosition_(detProp, **hitIt).X() << ", "
                  << HitPosition_(detProp, **hitIt).Y() << ")\n";
    }
  }

  CheckIsolatedHits_(showerHitsMap);

  if (fDebug > 2) {
    std::cout << "\nHits in order; after checking isolated hits...\n";
    for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::iterator showerHitsIt =
           showerHitsMap.begin();
         showerHitsIt != showerHitsMap.end();
         ++showerHitsIt) {
      std::cout << "  Plane " << showerHitsIt->first << '\n';
      for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt = showerHitsIt->second.begin();
           hitIt != showerHitsIt->second.end();
           ++hitIt)
        std::cout << "    Hit at (" << HitCoordinates_(**hitIt).X() << ", "
                  << HitCoordinates_(**hitIt).Y() << ") -- real wire " << (*hitIt)->WireID()
                  << ", hit position (" << HitPosition_(detProp, **hitIt).X() << ", "
                  << HitPosition_(detProp, **hitIt).Y() << ")\n";
    }
  }

  if (fDebug > 1) {
    std::cout << "After reversing and checking isolated hits, here are the "
                 "start and end points...\n";
    for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::iterator showerHitsIt =
           showerHitsMap.begin();
         showerHitsIt != showerHitsMap.end();
         ++showerHitsIt)
      std::cout << "  Plane " << showerHitsIt->first << " has start ("
                << HitCoordinates_(*showerHitsIt->second.front()).X() << ", "
                << HitCoordinates_(*showerHitsIt->second.front()).Y() << ") (real wire "
                << showerHitsIt->second.front()->WireID() << ") and end ("
                << HitCoordinates_(*showerHitsIt->second.back()).X() << ", "
                << HitCoordinates_(*showerHitsIt->second.back()).Y() << ")\n";
  }

  // Check for views in which the shower travels almost along the wire planes
  // (shown by a small relative wire width)
  std::map<double, int> wireWidths = RelativeWireWidth_(showerHitsMap);
  std::vector<int> badPlanes;
  if (fDebug > 1) std::cout << "Here are the relative wire widths... \n";
  for (auto const [relative_wire_width, plane] : wireWidths) {
    if (fDebug > 1)
      std::cout << "Plane " << plane << " has relative wire width " << relative_wire_width << '\n';
    if (relative_wire_width < 1 / (double)wireWidths.size()) badPlanes.push_back(plane);
  }

  std::map<int, std::vector<art::Ptr<recob::Hit>>> ignoredPlanes;
  if (badPlanes.size() == 1) {
    int const badPlane = badPlanes[0];
    if (fDebug > 1) std::cout << "Ignoring plane " << badPlane << " when orientating\n";
    ignoredPlanes[badPlane] = showerHitsMap.at(badPlane);
    showerHitsMap.erase(badPlane);
  }

  // Consider all possible permutations of planes (0/1, oriented
  // correctly/incorrectly)
  std::map<int, std::map<int, bool>> permutations = GetPlanePermutations_(detProp, showerHitsMap);

  // Go through all permutations until we have a satisfactory orientation
  auto const originalShowerHitsMap = showerHitsMap;

  int n = 0;
  while (!CheckShowerHits_(detProp, showerHitsMap) and
         n < TMath::Power(2, (int)showerHitsMap.size())) {
    if (fDebug > 1) std::cout << "Permutation " << n << '\n';
    for (int const plane : showerHitsMap | views::keys) {
      auto hits = originalShowerHitsMap.at(plane);
      if (permutations[n][plane] == 1) { std::reverse(hits.begin(), hits.end()); }
      showerHitsMap[plane] = hits;
    }
    ++n;
  }

  // Go back to original if still no consistency
  if (!CheckShowerHits_(detProp, showerHitsMap)) showerHitsMap = originalShowerHitsMap;

  if (fDebug > 2) {
    std::cout << "End of OrderShowerHits: here are the order of hits:\n";
    for (auto const& [plane, hitPtrs] : showerHitsMap) {
      std::cout << "  Plane " << plane << '\n';
      for (auto const& hit : hitPtrs | views::transform(to_element)) {
        std::cout << "    Hit (" << HitCoordinates_(hit).X() << " (real wire " << hit.WireID()
                  << "), " << HitCoordinates_(hit).Y() << ") -- pos ("
                  << HitPosition_(detProp, hit).X() << ", " << HitPosition_(detProp, hit).Y()
                  << ")\n";
      }
    }
  }

  if (ignoredPlanes.size())
    showerHitsMap[ignoredPlanes.begin()->first] = ignoredPlanes.begin()->second;

  return showerHitsMap;
}

void shower::EMShowerAlg::OrderShowerHits_(detinfo::DetectorPropertiesData const& detProp,
                                           std::vector<art::Ptr<recob::Hit>> const& shower,
                                           std::vector<art::Ptr<recob::Hit>>& showerHits,
                                           art::Ptr<recob::Vertex> const& vertex) const
{
  showerHits = FindOrderOfHits_(detProp, shower);

  // Find TPC for the vertex
  auto const& xyz = vertex->position();
  geo::TPCID tpc = fGeom->FindTPCAtPosition(xyz);
  if (!tpc.isValid && showerHits.size()) tpc = geo::TPCID(showerHits[0]->WireID());

  // Find hits in the same TPC
  art::Ptr<recob::Hit> hit0, hit1;
  for (auto& hit : showerHits) {
    if (hit->WireID().TPC == tpc.TPC) {
      if (hit0.isNull()) { hit0 = hit; }
      hit1 = hit;
    }
  }
  if (hit0.isNull() || hit1.isNull()) return;
  TVector2 coord0 = TVector2(hit0->WireID().Wire, hit0->PeakTime());
  TVector2 coord1 = TVector2(hit1->WireID().Wire, hit1->PeakTime());
  auto const& planeID = hit0->WireID().parentID();
  TVector2 coordvtx = TVector2(fWireReadoutGeom->Plane(planeID).WireCoordinate(xyz),
                               detProp.ConvertXToTicks(xyz.X(), planeID));
  if ((coord1 - coordvtx).Mod() < (coord0 - coordvtx).Mod()) {
    std::reverse(showerHits.begin(), showerHits.end());
  }
}

void shower::EMShowerAlg::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit>> const& showerHits,
                                               art::Ptr<recob::Vertex> const& vertex,
                                               std::vector<art::Ptr<recob::Hit>>& trackHits) const
{
  // Find TPC for the vertex
  auto const& xyz = vertex->position();
  geo::TPCID tpc = fGeom->FindTPCAtPosition(xyz);

  // vertex cannot be projected into a TPC, find the TPC that has the most hits
  if (!tpc.isValid) {
    std::map<geo::TPCID, unsigned int> tpcmap;
    unsigned maxhits = 0;
    for (auto const& hit : showerHits) {
      ++tpcmap[geo::TPCID(hit->WireID())];
    }
    for (auto const& t : tpcmap) {
      if (t.second > maxhits) {
        maxhits = t.second;
        tpc = t.first;
      }
    }
  }
  if (!tpc.isValid) return;

  double parm[2];
  int fitok = 0;
  std::vector<double> wfit;
  std::vector<double> tfit;
  std::vector<double> cfit;

  for (size_t i = 0; i < fNfitpass; ++i) {

    // Fit a straight line through hits
    unsigned int nhits = 0;
    for (auto& hit : showerHits) {
      if (hit->WireID().TPC == tpc.TPC) {
        TVector2 coord = HitCoordinates_(*hit);
        if (i == 0 ||
            (std::abs((coord.Y() - (parm[0] + coord.X() * parm[1])) * cos(atan(parm[1]))) <
             fToler[i - 1]) ||
            fitok == 1) {
          ++nhits;
          if (nhits == fNfithits[i] + 1) break;
          wfit.push_back(coord.X());
          tfit.push_back(coord.Y());
          cfit.push_back(1.);
          if (i == fNfitpass - 1) { trackHits.push_back(hit); }
        }
      }
    }

    if (i < fNfitpass - 1 && wfit.size()) {
      fitok = WeightedFit(wfit.size(), &wfit[0], &tfit[0], &cfit[0], &parm[0]);
    }
    wfit.clear();
    tfit.clear();
    cfit.clear();
  }
}

TVector2 shower::EMShowerAlg::HitCoordinates_(recob::Hit const& hit) const
{
  return TVector2(GlobalWire_(hit.WireID()), hit.PeakTime());
}

TVector2 shower::EMShowerAlg::HitPosition_(detinfo::DetectorPropertiesData const& detProp,
                                           recob::Hit const& hit) const
{
  geo::PlaneID planeID = hit.WireID().planeID();
  return HitPosition_(detProp, HitCoordinates_(hit), planeID);
}

TVector2 shower::EMShowerAlg::HitPosition_(detinfo::DetectorPropertiesData const& detProp,
                                           TVector2 const& pos,
                                           geo::PlaneID planeID) const
{
  return {pos.X() * fWireReadoutGeom->Plane(planeID).WirePitch(),
          detProp.ConvertTicksToX(pos.Y(), planeID)};
}

double shower::EMShowerAlg::GlobalWire_(const geo::WireID& wireID) const
{
  double globalWire = -999;

  // Induction
  if (fWireReadoutGeom->SignalType(wireID) == geo::kInduction) {
    auto const wireCenter = fWireReadoutGeom->Wire(wireID).GetCenter();
    globalWire = fWireReadoutGeom
                   ->Plane({wireID.Cryostat,
                            wireID.TPC % 2, // 0 or 1
                            wireID.Plane})
                   .WireCoordinate(wireCenter);
  }

  // Collection
  else {
    // FOR COLLECTION WIRES, HARD CODE THE GEOMETRY FOR GIVEN DETECTORS
    // THIS _SHOULD_ BE TEMPORARY. GLOBAL WIRE SUPPORT IS BEING ADDED TO THE LARSOFT GEOMETRY AND SHOULD BE AVAILABLE SOON
    if (fDetector == "dune35t") {
      unsigned int nwires =
        fWireReadoutGeom->Nwires(geo::PlaneID{wireID.Cryostat, 0, wireID.Plane});
      if (wireID.TPC == 0 or wireID.TPC == 1)
        globalWire = wireID.Wire;
      else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5)
        globalWire = nwires + wireID.Wire;
      else if (wireID.TPC == 6 or wireID.TPC == 7)
        globalWire = (2 * nwires) + wireID.Wire;
      else
        mf::LogError("BlurredClusterAlg")
          << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC
          << " (geometry" << fDetector << ")";
    }
    else if (fDetector == "dune10kt") {
      unsigned int nwires =
        fWireReadoutGeom->Nwires(geo::PlaneID{wireID.Cryostat, 0, wireID.Plane});
      // Detector geometry has four TPCs, two on top of each other, repeated along z...
      int block = wireID.TPC / 4;
      globalWire = (nwires * block) + wireID.Wire;
    }
    else {
      auto const wireCenter = fWireReadoutGeom->Wire(wireID).GetCenter();
      globalWire = fWireReadoutGeom
                     ->Plane({wireID.Cryostat,
                              wireID.TPC % 2, // 0 or 1
                              wireID.Plane})
                     .WireCoordinate(wireCenter);
    }
  }

  return globalWire;
}

std::map<double, int> shower::EMShowerAlg::RelativeWireWidth_(
  const std::map<int, std::vector<art::Ptr<recob::Hit>>>& showerHitsMap) const
{

  // Get the wire widths
  std::map<int, int> planeWireLength;
  for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::const_iterator showerHitsIt =
         showerHitsMap.begin();
       showerHitsIt != showerHitsMap.end();
       ++showerHitsIt)
    planeWireLength[showerHitsIt->first] =
      std::abs(HitCoordinates_(*showerHitsIt->second.front()).X() -
               HitCoordinates_(*showerHitsIt->second.back()).X());

  // Find the relative wire width for each plane with respect to the others
  std::map<int, double> planeOtherWireLengths;
  for (std::map<int, int>::iterator planeWireLengthIt = planeWireLength.begin();
       planeWireLengthIt != planeWireLength.end();
       ++planeWireLengthIt) {
    double quad = 0.;
    for (int plane = 0; plane < (int)fWireReadoutGeom->MaxPlanes(); ++plane)
      if (plane != planeWireLengthIt->first and planeWireLength.count(plane))
        quad += cet::square(planeWireLength[plane]);
    quad = std::sqrt(quad);
    planeOtherWireLengths[planeWireLengthIt->first] =
      planeWireLength[planeWireLengthIt->first] / (double)quad;
  }

  // Order these backwards
  std::map<double, int> wireWidthMap;
  for (std::map<int, std::vector<art::Ptr<recob::Hit>>>::const_iterator showerHitsIt =
         showerHitsMap.begin();
       showerHitsIt != showerHitsMap.end();
       ++showerHitsIt)
    wireWidthMap[planeOtherWireLengths.at(showerHitsIt->first)] = showerHitsIt->first;

  return wireWidthMap;
}

TVector2 shower::EMShowerAlg::ShowerDirection_(
  detinfo::DetectorPropertiesData const& detProp,
  const std::vector<art::Ptr<recob::Hit>>& showerHits) const
{

  TVector2 pos;
  double weight = 1;
  double sumx = 0., sumy = 0., sumx2 = 0., sumxy = 0., sumweight = 0.;
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator hit = showerHits.begin();
       hit != showerHits.end();
       ++hit) {
    //++nhits;
    pos = HitPosition_(detProp, **hit);
    weight = cet::square((*hit)->Integral());
    sumweight += weight;
    sumx += weight * pos.X();
    sumy += weight * pos.Y();
    sumx2 += weight * pos.X() * pos.X();
    sumxy += weight * pos.X() * pos.Y();
  }
  double gradient = (sumweight * sumxy - sumx * sumy) / (sumweight * sumx2 - sumx * sumx);
  TVector2 direction = TVector2(1, gradient).Unit();

  return direction;
}

TVector2 shower::EMShowerAlg::ShowerCenter_(
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::Hit>> const& showerHits) const
{

  TVector2 pos, chargePoint = TVector2(0, 0);
  double totalCharge = 0;
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator hit = showerHits.begin();
       hit != showerHits.end();
       ++hit) {
    pos = HitPosition_(detProp, **hit);
    chargePoint += (*hit)->Integral() * pos;
    totalCharge += (*hit)->Integral();
  }
  TVector2 centre = chargePoint / totalCharge;

  return centre;
}

double shower::EMShowerAlg::ShowerHitRMS_(detinfo::DetectorPropertiesData const& detProp,
                                          const std::vector<art::Ptr<recob::Hit>>& showerHits) const
{

  TVector2 direction = ShowerDirection_(detProp, showerHits);
  TVector2 centre = ShowerCenter_(detProp, showerHits);

  std::vector<double> distanceToAxis;
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator showerHitsIt = showerHits.begin();
       showerHitsIt != showerHits.end();
       ++showerHitsIt) {
    TVector2 proj = (HitPosition_(detProp, **showerHitsIt) - centre).Proj(direction) + centre;
    distanceToAxis.push_back((HitPosition_(detProp, **showerHitsIt) - proj).Mod());
  }
  double RMS = TMath::RMS(distanceToAxis.begin(), distanceToAxis.end());

  return RMS;
}

double shower::EMShowerAlg::ShowerHitRMSGradient_(
  detinfo::DetectorPropertiesData const& detProp,
  const std::vector<art::Ptr<recob::Hit>>& showerHits) const
{
  // Find a rough shower 'direction' and center
  TVector2 direction = ShowerDirection_(detProp, showerHits);

  // Bin the hits into discreet chunks
  int nShowerSegments = fNumShowerSegments;
  double lengthOfShower =
    (HitPosition_(detProp, *showerHits.back()) - HitPosition_(detProp, *showerHits.front())).Mod();
  double lengthOfSegment = lengthOfShower / (double)nShowerSegments;
  std::map<int, std::vector<art::Ptr<recob::Hit>>> showerSegments;
  std::map<int, double> segmentCharge;
  for (auto const& hitPtr : showerHits) {
    auto const& hit = *hitPtr;
    int const segment =
      (HitPosition_(detProp, hit) - HitPosition_(detProp, *showerHits.front())).Mod() /
      lengthOfSegment;
    showerSegments[segment].push_back(hitPtr);
    segmentCharge[segment] += hit.Integral();
  }

  TGraph* graph = new TGraph();
  std::vector<std::pair<int, double>> binVsRMS;

  // Loop over the bins to find the distribution of hits as the shower
  // progresses
  for (auto const& [segment, hitPtrs] : showerSegments) {

    // Get the mean position of the hits in this bin
    TVector2 meanPosition(0, 0);
    for (auto const& hit : hitPtrs | views::transform(to_element))
      meanPosition += HitPosition_(detProp, hit);
    meanPosition /= (double)hitPtrs.size();

    // Get the RMS of this bin
    std::vector<double> distanceToAxisBin;
    for (auto const& hit : hitPtrs | views::transform(to_element)) {
      TVector2 proj = (HitPosition_(detProp, hit) - meanPosition).Proj(direction) + meanPosition;
      distanceToAxisBin.push_back((HitPosition_(detProp, hit) - proj).Mod());
    }

    double RMSBin = TMath::RMS(distanceToAxisBin.begin(), distanceToAxisBin.end());
    binVsRMS.emplace_back(segment, RMSBin);
    if (fMakeRMSGradientPlot) graph->SetPoint(graph->GetN(), segment, RMSBin);
  }

  // Get the gradient of the RMS-bin plot
  double sumx = 0., sumy = 0., sumx2 = 0., sumxy = 0., sumweight = 0.;
  for (auto const& [bin, RMSBin] : binVsRMS) {
    double weight = segmentCharge.at(bin);
    sumweight += weight;
    sumx += weight * bin;
    sumy += weight * RMSBin;
    sumx2 += weight * bin * bin;
    sumxy += weight * bin * RMSBin;
  }
  double RMSgradient = (sumweight * sumxy - sumx * sumy) / (sumweight * sumx2 - sumx * sumx);

  if (fMakeRMSGradientPlot) {
    TVector2 direction = TVector2(1, RMSgradient).Unit();
    TCanvas* canv = new TCanvas();
    graph->Draw();
    graph->GetXaxis()->SetTitle("Shower segment");
    graph->GetYaxis()->SetTitle("RMS of hit distribution");
    TVector2 centre = TVector2(graph->GetMean(1), graph->GetMean(2));
    TLine line;
    line.SetLineColor(2);
    line.DrawLine(centre.X() - 1000 * direction.X(),
                  centre.Y() - 1000 * direction.Y(),
                  centre.X() + 1000 * direction.X(),
                  centre.Y() + 1000 * direction.Y());
    canv->SaveAs("RMSGradient.png");
    delete canv;
  }
  delete graph;

  return RMSgradient;
}

TVector2 shower::EMShowerAlg::Project3DPointOntoPlane_(
  detinfo::DetectorPropertiesData const& detProp,
  geo::Point_t const& point,
  int plane,
  int cryostat) const
{
  TVector2 wireTickPos{-999., -999.};

  geo::TPCID tpcID = fGeom->FindTPCAtPosition(point);
  int tpc = 0;
  if (tpcID.isValid)
    tpc = tpcID.TPC;
  else
    return wireTickPos;

  // Construct wire ID for this point projected onto the plane
  geo::PlaneID const planeID(cryostat, tpc, plane);
  geo::WireID wireID;
  try {
    wireID = fWireReadoutGeom->Plane(planeID).NearestWireID(point);
  }
  catch (geo::InvalidWireError const& e) {
    wireID = e.suggestedWireID(); // pick the closest valid wire
  }

  wireTickPos = TVector2(GlobalWire_(wireID), detProp.ConvertXToTicks(point.X(), planeID));

  return HitPosition_(detProp, wireTickPos, planeID);
}

int shower::EMShowerAlg::WorstPlane_(
  const std::map<int, std::vector<art::Ptr<recob::Hit>>>& showerHitsMap) const
{
  // Get the width of the shower in wire coordinate
  std::map<int, int> planeWireLength;
  std::map<int, double> planeOtherWireLengths;
  for (auto const& [plane, hits] : showerHitsMap)
    planeWireLength[plane] =
      std::abs(HitCoordinates_(*hits.front()).X() - HitCoordinates_(*hits.back()).X());
  for (std::map<int, int>::iterator planeWireLengthIt = planeWireLength.begin();
       planeWireLengthIt != planeWireLength.end();
       ++planeWireLengthIt) {
    double quad = 0.;
    for (int plane = 0; plane < (int)fWireReadoutGeom->MaxPlanes(); ++plane)
      if (plane != planeWireLengthIt->first and planeWireLength.count(plane))
        quad += cet::square(planeWireLength[plane]);
    quad = std::sqrt(quad);
    planeOtherWireLengths[planeWireLengthIt->first] =
      planeWireLength[planeWireLengthIt->first] / (double)quad;
  }

  if (fDebug > 1)
    for (auto const [plane, relative_width] : planeOtherWireLengths) {
      std::cout << "Plane " << plane << " has " << planeWireLength[plane]
                << " wire width and therefore has relative width in wire of " << relative_width
                << '\n';
    }

  std::map<double, int> wireWidthMap;
  for (int const plane : showerHitsMap | views::keys) {
    double wireWidth = planeWireLength.at(plane);
    wireWidthMap[wireWidth] = plane;
  }

  return wireWidthMap.begin()->second;
}

Int_t shower::EMShowerAlg::WeightedFit(const Int_t n,
                                       const Double_t* x,
                                       const Double_t* y,
                                       const Double_t* w,
                                       Double_t* parm) const
{

  Double_t sumx = 0.;
  Double_t sumx2 = 0.;
  Double_t sumy = 0.;
  Double_t sumxy = 0.;
  Double_t sumw = 0.;
  Double_t eparm[2];

  parm[0] = 0.;
  parm[1] = 0.;
  eparm[0] = 0.;
  eparm[1] = 0.;

  for (Int_t i = 0; i < n; i++) {
    sumx += x[i] * w[i];
    sumx2 += x[i] * x[i] * w[i];
    sumy += y[i] * w[i];
    sumxy += x[i] * y[i] * w[i];
    sumw += w[i];
  }

  if (sumx2 * sumw - sumx * sumx == 0.) return 1;
  if (sumx2 - sumx * sumx / sumw == 0.) return 1;

  parm[0] = (sumy * sumx2 - sumx * sumxy) / (sumx2 * sumw - sumx * sumx);
  parm[1] = (sumxy - sumx * sumy / sumw) / (sumx2 - sumx * sumx / sumw);

  eparm[0] = sumx2 * (sumx2 * sumw - sumx * sumx);
  eparm[1] = (sumx2 - sumx * sumx / sumw);

  if (eparm[0] < 0. || eparm[1] < 0.) return 1;

  eparm[0] = sqrt(eparm[0]) / (sumx2 * sumw - sumx * sumx);
  eparm[1] = sqrt(eparm[1]) / (sumx2 - sumx * sumx / sumw);

  return 0;
}

bool shower::EMShowerAlg::isCleanShower(std::vector<art::Ptr<recob::Hit>> const& hits) const
{
  if (hits.empty()) return false;
  if (hits.size() > 2000) return true;
  if (hits.size() < 20) return true;

  std::map<int, int> hitmap;
  unsigned nhits = 0;
  for (auto const& hit : hits | views::transform(to_element)) {
    ++nhits;
    if (nhits > 2) ++hitmap[hit.WireID().Wire];
    if (nhits == 20) break;
  }
  if (float(nhits - 2) / hitmap.size() > 1.4)
    return false;
  else
    return true;
}
