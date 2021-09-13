////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgCosmicTagger
// Author:      L. Whitehead (leigh.howard.whitehead@cern.ch),
//              R. Sulej (robert.sulej@cern.ch) March 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgCosmicTagger.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larreco/RecoAlg/PMAlg/PmaNode3D.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include <numeric> // std::accumulate

void
pma::PMAlgCosmicTagger::tag(detinfo::DetectorClocksData const& clockData,
                            pma::TrkCandidateColl& tracks)
{
  // Get the detector dimensions
  GetDimensions();

  mf::LogInfo("pma::PMAlgCosmicTagger")
    << "Passed " << tracks.size() << " tracks for tagging cosmics.";

  size_t n = 0;

  if (fTagOutOfDriftTracks) n += outOfDriftWindow(tracks);
  if (fTagFullHeightTracks) n += fullHeightCrossing(tracks);
  if (fTagFullWidthTracks) n += fullWidthCrossing(tracks);
  if (fTagFullLengthTracks) n += fullLengthCrossing(tracks);
  if (fTagNonBeamT0Tracks) n += nonBeamT0Tag(clockData, tracks);
  if (fTagTopFrontBack) n += tagTopFrontBack(tracks);
  if (fTagApparentStopper) n += tagApparentStopper(tracks);

  mf::LogInfo("pma::PMAlgCosmicTagger") << "...tagged " << n << " cosmic-like tracks.";
}

size_t
pma::PMAlgCosmicTagger::outOfDriftWindow(pma::TrkCandidateColl& tracks) const
{
  mf::LogInfo("pma::PMAlgCosmicTagger") << "   - tag tracks out of 1 drift window;";
  size_t n = 0;

  auto const* geom = lar::providerFrom<geo::Geometry>();

  for (auto& t : tracks.tracks()) {

    pma::Track3D& trk = *(t.Track());

    double min, max, p;
    bool node_out_of_drift_min = false;
    bool node_out_of_drift_max = false;
    for (size_t nidx = 0; nidx < trk.Nodes().size(); ++nidx) {
      auto const& node = *(trk.Nodes()[nidx]);
      auto const& tpcGeo = geom->TPC(node.TPC(), node.Cryo());
      // DetectDriftDirection returns a short int, but switch requires an int
      int driftDir = abs(tpcGeo.DetectDriftDirection());
      p = node.Point3D()[driftDir - 1];
      switch (driftDir) {
      case 1:
        min = tpcGeo.MinX();
        max = tpcGeo.MaxX();
        break;
      case 2:
        min = tpcGeo.MinY();
        max = tpcGeo.MaxY();
        break;
      case 3:
        min = tpcGeo.MinZ();
        max = tpcGeo.MaxZ();
        break;
      default:
        throw cet::exception("PMAlgCosmicTagger")
          << "Drift direction unknown: " << driftDir << std::endl;
      }
      if (p < min - fOutOfDriftMargin) { node_out_of_drift_min = true; }
      if (p > max + fOutOfDriftMargin) { node_out_of_drift_max = true; }
    }

    if (node_out_of_drift_min && node_out_of_drift_max) {
      trk.SetTagFlag(pma::Track3D::kCosmic);
      trk.SetTagFlag(pma::Track3D::kOutsideDrift_Complete);
      ++n;
    }
    else if (node_out_of_drift_min || node_out_of_drift_max) {
      trk.SetTagFlag(pma::Track3D::kCosmic);
      trk.SetTagFlag(pma::Track3D::kOutsideDrift_Partial);
      ++n;
    }
  }

  mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " tracks out of 1 drift window.";

  return n;
}

// Leigh: Make use of the fact that our cathode and anode crossing tracks have a reconstructed T0.
// Check to see if this time is consistent with the beam
size_t
pma::PMAlgCosmicTagger::nonBeamT0Tag(detinfo::DetectorClocksData const& clockData,
                                     pma::TrkCandidateColl& tracks) const
{
  size_t n = 0;

  // Search through all of the tracks
  for (auto& t : tracks.tracks()) {

    // Non zero T0 means we reconstructed it
    if (t.Track()->GetT0() != 0.0) {
      mf::LogInfo("pma::PMAlgCosmicTagger") << " - track with T0 = " << t.Track()->GetT0();

      if (fabs(t.Track()->GetT0() - trigger_offset(clockData)) > fNonBeamT0Margin) {
        ++n;
        t.Track()->SetTagFlag(pma::Track3D::kCosmic);
        t.Track()->SetTagFlag(pma::Track3D::kBeamIncompatible);
      }
    }
  }

  mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " non-beam T0 tracks.";
  return n;
}

size_t
pma::PMAlgCosmicTagger::tagTopFrontBack(pma::TrkCandidateColl& tracks) const
{

  size_t n = 0;

  auto const* geom = lar::providerFrom<geo::Geometry>();

  short int hIdx = ConvertDirToInt(geom->TPC(0, 0).HeightDir());
  short int lIdx = ConvertDirToInt(geom->TPC(0, 0).LengthDir());

  // Loop over the tracks
  for (auto& t : tracks.tracks()) {

    // Get the first and last positions from the track.
    auto const& node0 = *(t.Track()->Nodes()[0]);
    auto const& node1 = *(t.Track()->Nodes()[t.Track()->Nodes().size() - 1]);

    // Check which end is the vertex (assume the largest height)
    TVector3 vtx =
      (node0.Point3D()[hIdx] > node1.Point3D()[hIdx]) ? node0.Point3D() : node1.Point3D();
    TVector3 end =
      (node0.Point3D()[hIdx] <= node1.Point3D()[hIdx]) ? node0.Point3D() : node1.Point3D();

    // Check we have a track starting at the top of the detector
    bool top = isTopVertex(vtx, fTopFrontBackMargin, hIdx);

    // Check the track ends at the front or back of the detector
    bool frontBack = isFrontBackVertex(end, fTopFrontBackMargin, lIdx);

    // Check we path both criteria but without letting either the start or end of the track fulfill both
    if (top && frontBack) {
      ++n;
      t.Track()->SetTagFlag(pma::Track3D::kCosmic);
      t.Track()->SetTagFlag(pma::Track3D::kGeometry_YZ);
    }
  }

  mf::LogInfo("pma::PMAlgCosmicTagger")
    << " - Tagged " << n << " tracks crossing from top to front/back." << std::endl;

  return n;
}

size_t
pma::PMAlgCosmicTagger::tagApparentStopper(pma::TrkCandidateColl& tracks) const
{

  size_t n = 0;

  // Tracks entering from the top of the detector that stop in the fiducial volume
  // are very likely to be cosmics that leave through the APA, but have their
  // drift coordinate incorrectly set due to lack of T0
  auto const* geom = lar::providerFrom<geo::Geometry>();

  short int hIdx = ConvertDirToInt(geom->TPC(0, 0).HeightDir());

  // Loop over the tracks
  for (auto& t : tracks.tracks()) {

    // Get the first and last positions from the track.
    auto const& node0 = *(t.Track()->Nodes()[0]);
    auto const& node1 = *(t.Track()->Nodes()[t.Track()->Nodes().size() - 1]);

    // Check which end is the vertex (assume the largest height)
    TVector3 vtx =
      (node0.Point3D()[hIdx] > node1.Point3D()[hIdx]) ? node0.Point3D() : node1.Point3D();
    TVector3 end =
      (node0.Point3D()[hIdx] <= node1.Point3D()[hIdx]) ? node0.Point3D() : node1.Point3D();

    if (fabs(vtx[hIdx] - fDimensionsMax[hIdx]) < fApparentStopperMargin) {
      // Check the other element to see if it ends away from the bottom of the detector
      if (fabs(end[hIdx] - fDimensionsMin[hIdx]) > 5. * fApparentStopperMargin) {

        // We now need to loop over all of the tracks to see if any start within fStopperBuffer of our end point.
        bool foundTrack = false;
        for (auto const& tt : tracks.tracks()) {
          // Don't match with itself!
          if ((&tt) == (&t)) continue;

          // Compare this track with our main track
          TVector3 trkVtx = (tt.Track()->Nodes()[0])->Point3D();
          TVector3 trkEnd = (tt.Track()->Nodes()[tt.Track()->Nodes().size() - 1])->Point3D();

          if ((end - trkVtx).Mag() < fStopperBuffer || (end - trkEnd).Mag() < fStopperBuffer) {
            foundTrack = true;
            break;
          }
        }
        if (foundTrack) {
          // This isn't really a stopping particle, so move on
          continue;
        }

        // If we don't mind about tagging all stopping particles then this satisfies our requirements
        if (!fVetoActualStopper) {
          ++n;
          t.Track()->SetTagFlag(pma::Track3D::kCosmic);
          t.Track()->SetTagFlag(pma::Track3D::kGeometry_Y);
          continue;
        }

        // If we want to actually ignore the stopping particles, use de/dx...
        // Store the number of sigma from the mean for the final dedx point in each view
        std::vector<float> nSigmaPerView;

        // Loop over the views
        for (auto const view : geom->Views()) {

          // Get the dedx for this track and view
          std::map<size_t, std::vector<double>> dedx;
          t.Track()->GetRawdEdxSequence(dedx, view);

          std::vector<double> trk_dedx;

          for (int h = t.Track()->NextHit(-1, view); h != -1; h = t.Track()->NextHit(h, view)) {
            // If this is the last hit then this method won't work
            if (h > t.Track()->PrevHit(t.Track()->size(), view)) break;
            // Make sure we have a reasonable value
            if (dedx[h][5] / dedx[h][6] <= 0 || dedx[h][5] / dedx[h][6] > 1e6) continue;
            trk_dedx.push_back(dedx[h][5] / dedx[h][6]);
          }

          if (trk_dedx.size() == 0) {
            mf::LogInfo("pma::PMAlgCosmicTagger")
              << "View " << view << " has no hits." << std::endl;
            continue;
          }

          double sum = std::accumulate(std::begin(trk_dedx), std::end(trk_dedx), 0.0);
          double mean = sum / static_cast<double>(trk_dedx.size());
          double accum = 0.0;
          std::for_each(std::begin(trk_dedx), std::end(trk_dedx), [&](const double d) {
            accum += (d - mean) * (d - mean);
          });
          double stdev = sqrt(accum / static_cast<double>(trk_dedx.size() - 1));

          mf::LogInfo("pma::PMAlgCosmicTagger")
            << " View " << view << " has average dedx " << mean << " +/- " << stdev
            << " and final dedx " << trk_dedx[trk_dedx.size() - 1] << std::endl;

          nSigmaPerView.push_back(fabs((trk_dedx[trk_dedx.size() - 1] - mean) / stdev));
        }

        bool notStopper = true;
        short unsigned int n2Sigma = 0;
        short unsigned int n3Sigma = 0;
        for (auto const nSigma : nSigmaPerView) {
          if (nSigma >= 2.0) ++n2Sigma;
          if (nSigma >= 3.0) ++n3Sigma;
        }

        if (n3Sigma > 0) notStopper = false;
        if (n2Sigma == nSigmaPerView.size()) notStopper = false;

        if (notStopper) {
          ++n;
          t.Track()->SetTagFlag(pma::Track3D::kCosmic);
          t.Track()->SetTagFlag(pma::Track3D::kGeometry_Y);
        }
      } // Check on bottom position
    }   // Check on top position
  }     // End loop over tracks

  mf::LogInfo("pma::PMAlgCosmicTagger")
    << " - Tagged " << n << " tracks stopping in the detector after starting at the top."
    << std::endl;

  return n;
}

size_t
pma::PMAlgCosmicTagger::fullHeightCrossing(pma::TrkCandidateColl& tracks) const
{

  // Just use the first tpc to get the height dimension (instead of assuming y).
  auto const* geom = lar::providerFrom<geo::Geometry>();
  TVector3 dir = geom->TPC(0, 0).HeightDir();

  size_t n = fullCrossingTagger(tracks, ConvertDirToInt(dir));

  mf::LogInfo("pma::PMAlgCosmicTagger")
    << " - Tagged " << n << " tracks crossing the full detector height";
  return n;
}

size_t
pma::PMAlgCosmicTagger::fullWidthCrossing(pma::TrkCandidateColl& tracks) const
{

  // Just use the first tpc to get the width dimension (instead of assuming x).
  auto const* geom = lar::providerFrom<geo::Geometry>();
  TVector3 dir = geom->TPC(0, 0).WidthDir();

  size_t n = fullCrossingTagger(tracks, ConvertDirToInt(dir));

  mf::LogInfo("pma::PMAlgCosmicTagger")
    << " - Tagged " << n << " tracks crossing the full detector width";
  return n;
}

size_t
pma::PMAlgCosmicTagger::fullLengthCrossing(pma::TrkCandidateColl& tracks) const
{

  // Just use the first tpc to get the length dimension (instead of assuming z).
  auto const* geom = lar::providerFrom<geo::Geometry>();
  TVector3 dir = geom->TPC(0, 0).LengthDir();

  size_t n = fullCrossingTagger(tracks, ConvertDirToInt(dir));

  mf::LogInfo("pma::PMAlgCosmicTagger")
    << " - Tagged " << n << " tracks crossing the full detector length";
  return n;
}

size_t
pma::PMAlgCosmicTagger::fullCrossingTagger(pma::TrkCandidateColl& tracks, int direction) const
{

  if (direction == -1) {
    mf::LogWarning("pma::PMAlgCosmicTagger")
      << " - Could not recognise direction, not attempting to perform fullCrossingTagger.";
    return 0;
  }

  size_t n = 0;

  double detDim = fDimensionsMax[direction] - fDimensionsMin[direction];

  pma::Track3D::ETag dirTag = pma::Track3D::kNotTagged;
  switch (direction) {
  case 0: dirTag = pma::Track3D::kGeometry_XX; break;
  case 1: dirTag = pma::Track3D::kGeometry_YY; break;
  case 2: dirTag = pma::Track3D::kGeometry_ZZ; break;
  default: dirTag = pma::Track3D::kNotTagged; break;
  }

  // Loop over the tracks
  for (auto& t : tracks.tracks()) {

    // Get the first and last positions from the track.
    auto const& node0 = *(t.Track()->Nodes()[0]);
    auto const& node1 = *(t.Track()->Nodes()[t.Track()->Nodes().size() - 1]);

    // Get the length of the track in the requested direction
    double trkDim = fabs(node0.Point3D()[direction] - node1.Point3D()[direction]);

    if ((detDim - trkDim) < fFullCrossingMargin) {
      ++n;
      t.Track()->SetTagFlag(pma::Track3D::kCosmic);
      t.Track()->SetTagFlag(dirTag);
      mf::LogInfo("pma::PMAlgCosmicTagger") << " -- track tagged in direction " << direction
                                            << " with " << trkDim << " (c.f. " << detDim << ")";
    }
  }

  return n;
}

bool
pma::PMAlgCosmicTagger::isTopVertex(const TVector3& pos, double tolerance, short int dirIndx) const
{

  return (fabs(pos[dirIndx] - fDimensionsMax[dirIndx]) < tolerance);
}

bool
pma::PMAlgCosmicTagger::isFrontBackVertex(const TVector3& pos,
                                          double tolerance,
                                          short int dirIndx) const
{

  bool front = (fabs(pos[dirIndx] - fDimensionsMin[dirIndx]) < tolerance);
  bool back = (fabs(pos[dirIndx] - fDimensionsMax[dirIndx]) < tolerance);

  return front || back;
}

void
pma::PMAlgCosmicTagger::GetDimensions()
{

  // Need to find the minimum and maximum height values from the geometry.
  double minX = 1.e6;
  double maxX = -1.e6;
  double minY = 1.e6;
  double maxY = -1.e6;
  double minZ = 1.e6;
  double maxZ = -1.e6;

  auto const* geom = lar::providerFrom<geo::Geometry>();

  // Since we can stack TPCs, we can't just use geom::TPCGeom::Height()
  for (geo::TPCID const& tID : geom->IterateTPCIDs()) {
    geo::TPCGeo const& TPC = geom->TPC(tID);

    // We need to make sure we only include the real TPCs
    // We have dummy TPCs in the protoDUNE and DUNE geometries
    // The dummy ones have a drift distance of only ~13 cm.
    if (TPC.DriftDistance() < 25.0) { continue; }

    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5 * TPC.Length()};

    if (center[0] - tpcDim[0] < minX) minX = center[0] - tpcDim[0];
    if (center[0] + tpcDim[0] > maxX) maxX = center[0] + tpcDim[0];
    if (center[1] - tpcDim[1] < minY) minY = center[1] - tpcDim[1];
    if (center[1] + tpcDim[1] > maxY) maxY = center[1] + tpcDim[1];
    if (center[2] - tpcDim[2] < minZ) minZ = center[2] - tpcDim[2];
    if (center[2] + tpcDim[2] > maxZ) maxZ = center[2] + tpcDim[2];
  } // for all TPC

  fDimensionsMin.clear();
  fDimensionsMax.clear();
  fDimensionsMin.push_back(minX);
  fDimensionsMin.push_back(minY);
  fDimensionsMin.push_back(minZ);
  fDimensionsMax.push_back(maxX);
  fDimensionsMax.push_back(maxY);
  fDimensionsMax.push_back(maxZ);
}

short int
pma::PMAlgCosmicTagger::ConvertDirToInt(const TVector3& dir) const
{

  if (dir.X() > 0.99) return 0;
  if (dir.Y() > 0.99) return 1;
  if (dir.Z() > 0.99)
    return 2;

  else
    return -1;
}
