///////////////////////////////////////////////////////////////////////
///
/// \file   SpacePointAlg.cxx
///
/// \brief  Algorithm for generating space points from hits.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/SpacePointAlg.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
//\todo Remove include of BackTrackerService.h once this algorithm is stripped of test for MC
#include "lardata/RecoObjects/KHitTrack.h"
#include "lardata/RecoObjects/KHitWireX.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larsim/MCCheater/BackTrackerService.h"

//----------------------------------------------------------------------
// Constructor.
//
namespace trkf {

  SpacePointAlg::SpacePointAlg(const fhicl::ParameterSet& pset)
    : fMaxDT{pset.get<double>("MaxDT")}
    , fMaxS{pset.get<double>("MaxS")}
    , fMinViews{pset.get<int>("MinViews")}
    , fEnableU{pset.get<bool>("EnableU")}
    , fEnableV{pset.get<bool>("EnableV")}
    , fEnableW{pset.get<bool>("EnableW")}
    , fFilter{pset.get<bool>("Filter")}
    , fMerge{pset.get<bool>("Merge")}
    , fPreferColl{pset.get<bool>("PreferColl")}
    , fTickOffsetU{pset.get<double>("TickOffsetU", 0.)}
    , fTickOffsetV{pset.get<double>("TickOffsetV", 0.)}
    , fTickOffsetW{pset.get<double>("TickOffsetW", 0.)}
  {
    // Only allow one of fFilter and fMerge to be true.

    if (fFilter && fMerge)
      throw cet::exception("SpacePointAlg") << "Filter and Merge flags are both true.\n";

    // Report.

    std::cout << "SpacePointAlg configured with the following parameters:\n"
              << "  MaxDT = " << fMaxDT << "\n"
              << "  MaxS = " << fMaxS << "\n"
              << "  MinViews = " << fMinViews << "\n"
              << "  EnableU = " << fEnableU << "\n"
              << "  EnableV = " << fEnableV << "\n"
              << "  EnableW = " << fEnableW << "\n"
              << "  Filter = " << fFilter << "\n"
              << "  Merge = " << fMerge << "\n"
              << "  PreferColl = " << fPreferColl << "\n"
              << "  TickOffsetU = " << fTickOffsetU << "\n"
              << "  TickOffsetV = " << fTickOffsetV << "\n"
              << "  TickOffsetW = " << fTickOffsetW << std::endl;
  }

  //----------------------------------------------------------------------
  // Print geometry and properties constants.
  //
  void SpacePointAlg::update(detinfo::DetectorPropertiesData const& detProp) const
  {
    // Generate info report on first call only.

    static bool first = true;
    bool report = first;
    first = false;

    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;

    // Calculate and print geometry information.

    if (report) mf::LogInfo("SpacePointAlg") << "Updating geometry constants.\n";

    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      // Fill view-dependent quantities.

      geo::View_t view = plane.View();
      std::string viewname = "?";
      if (view == geo::kU) { viewname = "U"; }
      else if (view == geo::kV) {
        viewname = "V";
      }
      else if (view == geo::kZ) {
        viewname = "Z";
      }
      else
        throw cet::exception("SpacePointAlg") << "Bad view = " << view << "\n";

      std::string sigtypename = "?";
      geo::SigType_t sigtype = geom->SignalType(plane.ID());
      if (sigtype == geo::kInduction)
        sigtypename = "Induction";
      else if (sigtype == geo::kCollection)
        sigtypename = "Collection";
      else
        throw cet::exception("SpacePointAlg") << "Bad signal type = " << sigtype << "\n";

      std::string orientname = "?";
      geo::Orient_t orient = plane.Orientation();
      if (orient == geo::kVertical)
        orientname = "Vertical";
      else if (orient == geo::kHorizontal)
        orientname = "Horizontal";
      else
        throw cet::exception("SpacePointAlg") << "Bad orientation = " << orient << "\n";

      if (report) {
        auto const xyz = plane.GetCenter();
        auto const& tpcgeom = geom->TPC(plane.ID());
        mf::LogInfo("SpacePointAlg")
          << '\n'
          << plane.ID() << '\n'
          << "  View: " << viewname << "\n"
          << "  SignalType: " << sigtypename << "\n"
          << "  Orientation: " << orientname << "\n"
          << "  Plane location: " << xyz.X() << "\n"
          << "  Plane pitch: " << tpcgeom.Plane0Pitch(plane.ID().Plane) << "\n"
          << "  Wire angle: " << plane.Wire(0).ThetaZ() << "\n"
          << "  Wire pitch: " << tpcgeom.WirePitch() << "\n"
          << "  Time offset: " << detProp.GetXTicksOffset(plane.ID()) << "\n";
      }

      if (orient != geo::kVertical)
        throw cet::exception("SpacePointAlg") << "Horizontal wire geometry not implemented.\n";
    } // end loop over planes
  }

  //----------------------------------------------------------------------
  // Get corrected time for the specified hit.
  double SpacePointAlg::correctedTime(detinfo::DetectorPropertiesData const& detProp,
                                      const recob::Hit& hit) const
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;

    // Correct time for trigger offset and plane-dependent time offsets.

    double t = hit.PeakTime() - detProp.GetXTicksOffset(hit.WireID());
    if (hit.View() == geo::kU)
      t -= fTickOffsetU;
    else if (hit.View() == geo::kV)
      t -= fTickOffsetV;
    else if (hit.View() == geo::kW)
      t -= fTickOffsetW;

    return t;
  }

  //----------------------------------------------------------------------
  // Spatial separation of hits (zero if two or fewer).
  double SpacePointAlg::separation(const art::PtrVector<recob::Hit>& hits) const
  {
    // Get geometry service.

    art::ServiceHandle<geo::Geometry const> geom;

    // Trivial case - fewer than three hits.

    if (hits.size() < 3) return 0.;

    // Error case - more than three hits.

    if (hits.size() > 3) {
      mf::LogError("SpacePointAlg") << "Method separation called with more than three htis.";
      return 0.;
    }

    // Got exactly three hits.

    // Calculate angles and distance of each hit from origin.

    double dist[3] = {0., 0., 0.};
    double sinth[3] = {0., 0., 0.};
    double costh[3] = {0., 0., 0.};
    unsigned int cstats[3];
    unsigned int tpcs[3];
    unsigned int planes[3];

    for (int i = 0; i < 3; ++i) {

      // Get tpc, plane, wire.

      const recob::Hit& hit = *(hits[i]);
      const geo::WireGeo& wgeom = geom->WireIDToWireGeo(hit.WireID());
      cstats[i] = hit.WireID().Cryostat;
      tpcs[i] = hit.WireID().TPC;
      planes[i] = hit.WireID().Plane;

      // Check tpc and plane errors.

      for (int j = 0; j < i; ++j) {

        if (cstats[j] != hit.WireID().Cryostat) {
          mf::LogError("SpacePointAlg")
            << "Method separation called with hits from multiple cryostats..";
          return 0.;
        }

        if (tpcs[j] != hit.WireID().TPC) {
          mf::LogError("SpacePointAlg")
            << "Method separation called with hits from multiple tpcs..";
          return 0.;
        }

        if (planes[j] == hit.WireID().Plane) {
          mf::LogError("SpacePointAlg")
            << "Method separation called with hits from the same plane..";
          return 0.;
        }
      }

      // Get angles and distance of wire.

      double const hl = wgeom.HalfL();
      auto const xyz = wgeom.GetCenter();
      auto const xyz1 = wgeom.GetEnd();
      double s = (xyz1.Y() - xyz.Y()) / hl;
      double c = (xyz1.Z() - xyz.Z()) / hl;
      sinth[hit.WireID().Plane] = s;
      costh[hit.WireID().Plane] = c;
      dist[hit.WireID().Plane] = xyz.Z() * s - xyz.Y() * c;
    }

    double S = ((sinth[1] * costh[2] - costh[1] * sinth[2]) * dist[0] +
                (sinth[2] * costh[0] - costh[2] * sinth[0]) * dist[1] +
                (sinth[0] * costh[1] - costh[0] * sinth[1]) * dist[2]);
    return S;
  }

  //----------------------------------------------------------------------
  // Check hits for compatibility.
  // Check hits pairwise for different views and maximum time difference.
  // Check three hits for spatial compatibility.
  bool SpacePointAlg::compatible(detinfo::DetectorPropertiesData const& detProp,
                                 const art::PtrVector<recob::Hit>& hits,
                                 bool useMC) const
  {
    art::ServiceHandle<geo::Geometry const> geom;

    int nhits = hits.size();

    // Fewer than two or more than three hits can never be compatible.

    bool result = nhits >= 2 && nhits <= 3;
    bool mc_ok = true;
    unsigned int tpc = 0;
    unsigned int cstat = 0;

    if (result) {

      // First do pairwise tests.
      // Do double loop over hits.

      for (int ihit1 = 0; result && ihit1 < nhits - 1; ++ihit1) {
        const recob::Hit& hit1 = *(hits[ihit1]);
        geo::WireID hit1WireID = hit1.WireID();
        geo::View_t view1 = hit1.View();

        double t1 = hit1.PeakTime() -
                    detProp.GetXTicksOffset(hit1WireID.Plane, hit1WireID.TPC, hit1WireID.Cryostat);

        // If using mc information, get a collection of track ids for hit 1.
        // If not using mc information, this section of code will trigger the
        // insertion of a single invalid HitMCInfo object into fHitMCMap.

        const HitMCInfo& mcinfo1 = fHitMCMap[(useMC ? &hit1 : 0)];
        const std::vector<int>& tid1 = mcinfo1.trackIDs;
        bool only_neg1 = tid1.size() > 0 && tid1.back() < 0;

        // Loop over second hit.

        for (int ihit2 = ihit1 + 1; result && ihit2 < nhits; ++ihit2) {
          const recob::Hit& hit2 = *(hits[ihit2]);
          geo::WireID hit2WireID = hit2.WireID();
          geo::View_t view2 = hit2.View();

          // Test for same tpc and different views.

          result = result && hit1WireID.TPC == hit2WireID.TPC && view1 != view2 &&
                   hit1WireID.Cryostat == hit2WireID.Cryostat;
          if (result) {

            // Remember which tpc and cryostat we are in.

            tpc = hit1WireID.TPC;
            cstat = hit1WireID.Cryostat;

            double t2 = hit2.PeakTime() - detProp.GetXTicksOffset(
                                            hit2WireID.Plane, hit2WireID.TPC, hit2WireID.Cryostat);

            // Test maximum time difference.

            result = result && std::abs(t1 - t2) <= fMaxDT;

            // Test mc truth.

            if (result && useMC) {

              // Test whether hits have a common parent track id.

              const HitMCInfo& mcinfo2 = fHitMCMap[&hit2];
              std::vector<int> tid2 = mcinfo2.trackIDs;
              bool only_neg2 = tid2.size() > 0 && tid2.back() < 0;
              std::vector<int>::iterator it = std::set_intersection(
                tid1.begin(), tid1.end(), tid2.begin(), tid2.end(), tid2.begin());
              tid2.resize(it - tid2.begin());

              // Hits are compatible if they have parents in common.
              // If the only parent id in common is negative (-999),
              // then hits are compatible only if both hits have only
              // negative parent tracks.

              bool only_neg3 = tid2.size() > 0 && tid2.back() < 0;
              mc_ok = tid2.size() > 0 && (!only_neg3 || (only_neg1 && only_neg2));
              result = result && mc_ok;

              // If we are still OK, check that either hit is
              // the nearest neighbor of the other.

              if (result) {
                result = mcinfo1.pchit[hit2WireID.Plane] == &hit2 ||
                         mcinfo2.pchit[hit1WireID.Plane] == &hit1;
              }
            }
          }
        }
      }

      // If there are exactly three hits, and they pass pairwise tests, check
      // for spatial compatibility.

      if (result && nhits == 3) {

        // Loop over hits.

        double dist[3] = {0., 0., 0.};
        double sinth[3] = {0., 0., 0.};
        double costh[3] = {0., 0., 0.};

        for (int i = 0; i < 3; ++i) {

          // Get tpc, plane, wire.

          const recob::Hit& hit = *(hits[i]);
          geo::WireID hitWireID = hit.WireID();

          const geo::WireGeo& wgeom = geom->WireIDToWireGeo(hit.WireID());
          if ((hitWireID.TPC != tpc) || (hitWireID.Cryostat != cstat))
            throw cet::exception("SpacePointAlg") << "compatible(): geometry mismatch\n";

          // Get angles and distance of wire.

          double const hl = wgeom.HalfL();
          auto const xyz = wgeom.GetCenter();
          auto const xyz1 = wgeom.GetEnd();
          double s = (xyz1.Y() - xyz.Y()) / hl;
          double c = (xyz1.Z() - xyz.Z()) / hl;
          sinth[hit.WireID().Plane] = s;
          costh[hit.WireID().Plane] = c;
          dist[hit.WireID().Plane] = xyz.Z() * s - xyz.Y() * c;
        }

        // Do space cut.

        double S = ((sinth[1] * costh[2] - costh[1] * sinth[2]) * dist[0] +
                    (sinth[2] * costh[0] - costh[2] * sinth[0]) * dist[1] +
                    (sinth[0] * costh[1] - costh[0] * sinth[1]) * dist[2]);

        result = result && std::abs(S) < fMaxS;
      }
    }

    // Done.

    return result;
  }

  //----------------------------------------------------------------------
  // Fill one space point using a colleciton of hits.
  // Assume points have already been tested for compatibility.
  //
  void SpacePointAlg::fillSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                     const art::PtrVector<recob::Hit>& hits,
                                     std::vector<recob::SpacePoint>& sptv,
                                     int sptid) const
  {
    art::ServiceHandle<geo::Geometry const> geom;

    double timePitch = detProp.GetXTicksCoefficient();

    int nhits = hits.size();

    // Remember associated hits internally.

    if (fSptHitMap.find(sptid) != fSptHitMap.end())
      throw cet::exception("SpacePointAlg") << "fillSpacePoint(): hit already present!\n";
    fSptHitMap[sptid] = hits;

    // Calculate position and error matrix.

    double xyz[3] = {0., 0., 0.};
    double errxyz[6] = {0., 0., 0., 0., 0., 0.};

    // Calculate x using drift times.
    // Loop over all hits and calculate the weighted average drift time.
    // Also calculate time variance and chisquare.

    double sumt2w = 0.;
    double sumtw = 0.;
    double sumw = 0.;

    for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
         ++ihit) {

      const recob::Hit& hit = **ihit;
      geo::WireID hitWireID = hit.WireID();

      // Correct time for trigger offset and view-dependent time offsets.

      double t0 = detProp.GetXTicksOffset(hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);
      double t = hit.PeakTime() - t0;
      double et = hit.SigmaPeakTime();
      double w = 1. / (et * et);

      sumt2w += w * t * t;
      sumtw += w * t;
      sumw += w;
    }

    double drift_time = 0.;
    double var_time = 0.;
    double chisq = 0.;
    if (sumw != 0.) {
      drift_time = sumtw / sumw;
      //var_time = sumt2w / sumw - drift_time * drift_time;
      var_time = 1. / sumw;
      if (var_time < 0.) var_time = 0.;
      chisq = sumt2w - sumtw * drift_time;
      if (chisq < 0.) chisq = 0.;
    }
    xyz[0] = drift_time * timePitch;
    errxyz[0] = var_time * timePitch * timePitch;

    // Calculate y, z using wires (need at least two hits).

    if (nhits >= 2) {

      // Calculate y and z by chisquare minimization of wire coordinates.

      double sus = 0.; // sum w_i u_i sin_th_i
      double suc = 0.; // sum w_i u_i cos_th_i
      double sc2 = 0.; // sum w_i cos2_th_i
      double ss2 = 0.; // sum w_i sin2_th_i
      double ssc = 0.; // sum w_i sin_th_i cos_th_i

      // Loop over points.

      for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
           ++ihit) {

        const recob::Hit& hit = **ihit;
        geo::WireID hitWireID = hit.WireID();
        const geo::WireGeo& wgeom = geom->WireIDToWireGeo(hit.WireID());

        // Calculate angle and wire coordinate in this view.

        double const hl = wgeom.HalfL();
        auto const cen = wgeom.GetCenter();
        auto const cen1 = wgeom.GetEnd();
        double s = (cen1.Y() - cen.Y()) / hl;
        double c = (cen1.Z() - cen.Z()) / hl;
        double u = cen.Z() * s - cen.Y() * c;
        double eu = geom->WirePitch(hitWireID.asPlaneID()) / std::sqrt(12.);
        double w = 1. / (eu * eu);

        // Summations

        sus += w * u * s;
        suc += w * u * c;
        sc2 += w * c * c;
        ss2 += w * s * s;
        ssc += w * s * c;
      }

      // Calculate y,z

      double denom = sc2 * ss2 - ssc * ssc;
      if (denom != 0.) {
        xyz[1] = (-suc * ss2 + sus * ssc) / denom;
        xyz[2] = (sus * sc2 - suc * ssc) / denom;
        errxyz[2] = ss2 / denom;
        errxyz[4] = ssc / denom;
        errxyz[5] = sc2 / denom;
      }

      // Make space point.

      recob::SpacePoint spt(xyz, errxyz, chisq, sptid);
      sptv.push_back(spt);
    }
    return;
  }

  /// Fill a collection of space points.
  ///
  /// Arguments:
  ///
  /// spts   - Collection of space points to fill.
  /// sptalg - Space point algorithm object.
  ///
  /// This method uses the hits contained in this track to construct
  /// space points.
  ///
  /// This method does not have any knowledge of what constitutes a
  /// good space point, except that Hits are required to be
  /// consecutive when sorted by path distance, and space points are
  /// required to pass compatibility tests used by the space point
  /// algorithm object.  This method will make space points from
  /// either two or three Hits (even for three-plane detectors), if
  /// the space point algorithm is configured to allow it.
  ///
  void SpacePointAlg::fillSpacePoints(detinfo::DetectorPropertiesData const& detProp,
                                      std::vector<recob::SpacePoint>& spts,
                                      std::multimap<double, KHitTrack> const& trackMap) const
  {
    // Loop over KHitTracks.

    art::PtrVector<recob::Hit> hits;
    art::PtrVector<recob::Hit> compatible_hits;
    for (std::multimap<double, KHitTrack>::const_iterator it = trackMap.begin();
         it != trackMap.end();
         ++it) {
      const KHitTrack& track = (*it).second;

      // Extrack Hit from track.

      const std::shared_ptr<const KHitBase>& hit = track.getHit();
      const KHitWireX* phit = dynamic_cast<const KHitWireX*>(&*hit);
      if (phit != 0) {
        const art::Ptr<recob::Hit> prhit = phit->getHit();

        // Test this hit for compatibility.

        hits.push_back(prhit);
        bool ok = this->compatible(detProp, hits);
        if (!ok) {

          // The new hit is not compatible.  Make a space point out of
          // the last known compatible hits, provided there are at least
          // two.

          if (compatible_hits.size() >= 2) {
            this->fillSpacePoint(detProp, compatible_hits, spts, this->numHitMap());
            compatible_hits.clear();
          }

          // Forget about any previous hits.

          hits.clear();
          hits.push_back(prhit);
        }

        // Update the list of known compatible hits.

        compatible_hits = hits;
      }
    }

    // Maybe make one final space point.

    if (compatible_hits.size() >= 2) {
      this->fillSpacePoint(detProp, compatible_hits, spts, this->numHitMap());
    }
  }

  //----------------------------------------------------------------------
  // Fill one space point using a colleciton of hits.
  // Assume points have already been tested for compatibility.
  // This version assumes there can be multiple hits per view,
  // and gives unequal weight to different hits.
  //
  void SpacePointAlg::fillComplexSpacePoint(detinfo::DetectorPropertiesData const& detProp,
                                            const art::PtrVector<recob::Hit>& hits,
                                            std::vector<recob::SpacePoint>& sptv,
                                            int sptid) const
  {
    art::ServiceHandle<geo::Geometry const> geom;

    // Calculate time pitch.

    double timePitch = detProp.GetXTicksCoefficient(); // cm / tick

    // Figure out which tpc we are in.

    unsigned int tpc0 = 0;
    unsigned int cstat0 = 0;
    int nhits = hits.size();
    if (nhits > 0) {
      tpc0 = hits.front()->WireID().TPC;
      cstat0 = hits.front()->WireID().Cryostat;
    }

    // Remember associated hits internally.

    if (fSptHitMap.count(sptid) != 0)
      throw cet::exception("SpacePointAlg") << "fillComplexSpacePoint(): hit already present!\n";
    fSptHitMap[sptid] = hits;

    // Do a preliminary scan of hits.
    // Determine weight given to hits in each view.

    unsigned int nplanes = geom->TPC(geo::TPCID(cstat0, tpc0)).Nplanes();
    std::vector<int> numhits(nplanes, 0);
    std::vector<double> weight(nplanes, 0.);

    for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
         ++ihit) {

      const recob::Hit& hit = **ihit;
      geo::WireID hitWireID = hit.WireID();
      // kept as assertions for performance reasons
      assert(hitWireID.Cryostat == cstat0);
      assert(hitWireID.TPC == tpc0);
      assert(hitWireID.Plane < nplanes);
      ++numhits[hitWireID.Plane];
    }

    for (unsigned int plane = 0; plane < nplanes; ++plane) {
      double np = numhits[plane];
      if (np > 0.) weight[plane] = 1. / (np * np * np);
    }

    // Calculate position and error matrix.

    double xyz[3] = {0., 0., 0.};
    double errxyz[6] = {0., 0., 0., 0., 0., 0.};

    // Calculate x using drift times.
    // Loop over all hits and calculate the weighted average drift time.
    // Also calculate time variance and chisquare.

    double sumt2w = 0.;
    double sumtw = 0.;
    double sumw = 0.;

    for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
         ++ihit) {

      const recob::Hit& hit = **ihit;
      geo::WireID hitWireID = hit.WireID();

      // Correct time for trigger offset and view-dependent time offsets.

      double t0 = detProp.GetXTicksOffset(hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);
      double t = hit.PeakTime() - t0;
      double et = hit.SigmaPeakTime();
      double w = weight[hitWireID.Plane] / (et * et);

      sumt2w += w * t * t;
      sumtw += w * t;
      sumw += w;
    }

    double drift_time = 0.;
    double var_time = 0.;
    double chisq = 0.;
    if (sumw != 0.) {
      drift_time = sumtw / sumw;
      var_time = sumt2w / sumw - drift_time * drift_time;
      if (var_time < 0.) var_time = 0.;
      chisq = sumt2w - sumtw * drift_time;
      if (chisq < 0.) chisq = 0.;
    }
    xyz[0] = drift_time * timePitch;
    errxyz[0] = var_time * timePitch * timePitch;

    // Calculate y, z using wires (need at least two hits).

    if (nhits >= 2) {

      // Calculate y and z by chisquare minimization of wire coordinates.

      double sus = 0.; // sum w_i u_i sin_th_i
      double suc = 0.; // sum w_i u_i cos_th_i
      double sc2 = 0.; // sum w_i cos2_th_i
      double ss2 = 0.; // sum w_i sin2_th_i
      double ssc = 0.; // sum w_i sin_th_i cos_th_i

      // Loop over points.

      for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
           ++ihit) {

        const recob::Hit& hit = **ihit;
        geo::WireID hitWireID = hit.WireID();
        const geo::WireGeo& wgeom = geom->WireIDToWireGeo(hit.WireID());

        // Calculate angle and wire coordinate in this view.

        double const hl = wgeom.HalfL();
        auto const cen = wgeom.GetCenter();
        auto const cen1 = wgeom.GetEnd();
        double s = (cen1.Y() - cen.Y()) / hl;
        double c = (cen1.Z() - cen.Z()) / hl;
        double u = cen.Z() * s - cen.Y() * c;
        double eu = geom->WirePitch(hitWireID.asPlaneID()) / std::sqrt(12.);
        double w = weight[hitWireID.Plane] / (eu * eu);

        // Summations

        sus += w * u * s;
        suc += w * u * c;
        sc2 += w * c * c;
        ss2 += w * s * s;
        ssc += w * s * c;
      }

      // Calculate y,z

      double denom = sc2 * ss2 - ssc * ssc;
      if (denom != 0.) {
        xyz[1] = (-suc * ss2 + sus * ssc) / denom;
        xyz[2] = (sus * sc2 - suc * ssc) / denom;
        errxyz[2] = ss2 / denom;
        errxyz[4] = ssc / denom;
        errxyz[5] = sc2 / denom;
      }

      // Make space point.

      recob::SpacePoint spt(xyz, errxyz, chisq, sptid);
      sptv.push_back(spt);
    }
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (non-mc-truth version).
  //
  void SpacePointAlg::makeSpacePoints(detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      const art::PtrVector<recob::Hit>& hits,
                                      std::vector<recob::SpacePoint>& spts) const
  {
    makeSpacePoints(clockData, detProp, hits, spts, false);
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (mc-truth version).
  //
  void SpacePointAlg::makeMCTruthSpacePoints(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             const art::PtrVector<recob::Hit>& hits,
                                             std::vector<recob::SpacePoint>& spts) const
  {
    makeSpacePoints(clockData, detProp, hits, spts, true);
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (general version).
  //
  void SpacePointAlg::makeSpacePoints(detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      const art::PtrVector<recob::Hit>& hits,
                                      std::vector<recob::SpacePoint>& spts,
                                      bool useMC) const
  {
    art::ServiceHandle<geo::Geometry const> geom;

    fSptHitMap.clear();

    // Print diagnostic information.

    update(detProp);

    // First make sure result vector is empty.

    spts.erase(spts.begin(), spts.end());

    // Statistics.

    int n2 = 0;     // Number of two-hit space points.
    int n3 = 0;     // Number of three-hit space points.
    int n2filt = 0; // Number of two-hit space points after filtering/merging.
    int n3filt = 0; // Number of three-hit space pointe after filtering/merging.

    // Sort hits into maps indexed by [cryostat][tpc][plane][wire].
    // If using mc information, also generate maps of sim::IDEs and mc
    // position indexed by hit.

    std::vector<std::vector<std::vector<std::multimap<unsigned int, art::Ptr<recob::Hit>>>>> hitmap;
    fHitMCMap.clear();

    unsigned int ncstat = geom->Ncryostats();
    hitmap.resize(ncstat);
    for (auto const& cryoid : geom->Iterate<geo::CryostatID>()) {
      unsigned int cstat = cryoid.Cryostat;
      unsigned int ntpc = geom->Cryostat(cryoid).NTPC();
      hitmap[cstat].resize(ntpc);
      for (unsigned int tpc = 0; tpc < ntpc; ++tpc) {
        int nplane = geom->Cryostat(cryoid).TPC(tpc).Nplanes();
        hitmap[cstat][tpc].resize(nplane);
      }
    }

    for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
         ++ihit) {
      const art::Ptr<recob::Hit>& phit = *ihit;
      geo::View_t view = phit->View();
      if ((view == geo::kU && fEnableU) || (view == geo::kV && fEnableV) ||
          (view == geo::kZ && fEnableW)) {
        geo::WireID phitWireID = phit->WireID();
        hitmap[phitWireID.Cryostat][phitWireID.TPC][phitWireID.Plane].insert(
          std::make_pair(phitWireID.Wire, phit));
      }
    }

    // Fill mc information, including IDEs and closest neighbors
    // of each hit.
    ///\todo Why are we still checking on whether this is MC or not?
    ///\todo Such checks should not be in reconstruction code.
    if (useMC) {
      art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

      // First loop over hits and fill track ids and mc position.
      for (auto const& id : geom->Iterate<geo::PlaneID>()) {
        auto const [cstat, tpc, plane] = std::make_tuple(id.Cryostat, id.TPC, id.Plane);
        int nplane = geom->TPC(id).Nplanes();
        for (std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator ihit =
               hitmap[cstat][tpc][plane].begin();
             ihit != hitmap[cstat][tpc][plane].end();
             ++ihit) {
          const art::Ptr<recob::Hit>& phit = ihit->second;
          const recob::Hit& hit = *phit;
          HitMCInfo& mcinfo = fHitMCMap[&hit]; // Default HitMCInfo.

          // Fill default nearest neighbor information (i.e. none).

          mcinfo.pchit.resize(nplane, 0);
          mcinfo.dist2.resize(nplane, 1.e20);

          // Get sim::IDEs for this hit.

          std::vector<sim::IDE> ides = bt_serv->HitToAvgSimIDEs(clockData, phit);

          // Get sorted track ids. for this hit.

          mcinfo.trackIDs.reserve(ides.size());
          for (std::vector<sim::IDE>::const_iterator i = ides.begin(); i != ides.end(); ++i)
            mcinfo.trackIDs.push_back(i->trackID);
          sort(mcinfo.trackIDs.begin(), mcinfo.trackIDs.end());

          // Get position of ionization for this hit.

          try {
            mcinfo.xyz = bt_serv->SimIDEsToXYZ(ides);
          }
          catch (cet::exception& x) {
            mcinfo.xyz.clear();
          }
        } // end loop over ihit
      }   // end loop oer planes

      // Loop over hits again and fill nearest neighbor information for real.
      for (auto const& id : geom->Iterate<geo::PlaneID>()) {
        auto const [cstat, tpc, plane] = std::make_tuple(id.Cryostat, id.TPC, id.Plane);
        int nplane = geom->TPC(id).Nplanes();
        for (std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator ihit =
               hitmap[cstat][tpc][plane].begin();
             ihit != hitmap[cstat][tpc][plane].end();
             ++ihit) {
          const art::Ptr<recob::Hit>& phit = ihit->second;
          const recob::Hit& hit = *phit;
          HitMCInfo& mcinfo = fHitMCMap[&hit];
          if (mcinfo.xyz.size() != 0) {
            assert(mcinfo.xyz.size() == 3);

            // Fill nearest neighbor information for this hit.

            for (int plane2 = 0; plane2 < nplane; ++plane2) {
              for (std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator jhit =
                     hitmap[cstat][tpc][plane2].begin();
                   jhit != hitmap[cstat][tpc][plane2].end();
                   ++jhit) {
                const art::Ptr<recob::Hit>& phit2 = jhit->second;
                const recob::Hit& hit2 = *phit2;
                const HitMCInfo& mcinfo2 = fHitMCMap[&hit2];

                if (mcinfo2.xyz.size() != 0) {
                  assert(mcinfo2.xyz.size() == 3);
                  double dx = mcinfo.xyz[0] - mcinfo2.xyz[0];
                  double dy = mcinfo.xyz[1] - mcinfo2.xyz[1];
                  double dz = mcinfo.xyz[2] - mcinfo2.xyz[2];
                  double dist2 = dx * dx + dy * dy + dz * dz;
                  if (dist2 < mcinfo.dist2[plane2]) {
                    mcinfo.dist2[plane2] = dist2;
                    mcinfo.pchit[plane2] = &hit2;
                  }
                } // end if mcinfo2.xyz valid
              }   // end loop over jhit
            }     // end loop over plane2
          }       // end if mcinfo.xyz valid.
        }         // end loop over ihit
      }           // end loop over planes
    }             // end if MC

    // use mf::LogDebug instead of MF_LOG_DEBUG because we reuse it in many lines
    // insertions are protected by mf::isDebugEnabled()
    mf::LogDebug debug("SpacePointAlg");
    if (mf::isDebugEnabled()) {
      debug << "Total hits = " << hits.size() << "\n\n";
      for (auto const& id : geom->Iterate<geo::TPCID>()) {
        auto const [cstat, tpc] = std::make_tuple(id.Cryostat, id.TPC);
        int nplane = hitmap[cstat][tpc].size();
        for (int plane = 0; plane < nplane; ++plane) {
          debug << "TPC, Plane: " << tpc << ", " << plane
                << ", hits = " << hitmap[cstat][tpc][plane].size() << "\n";
        }
      } // end loop over TPCs
    }   // if debug

    // Make empty multimap from hit pointer on preferred
    // (most-populated or collection) plane to space points that
    // include that hit (used for sorting, filtering, and
    // merging).

    typedef const recob::Hit* sptkey_type;
    std::multimap<sptkey_type, recob::SpacePoint> sptmap;
    std::set<sptkey_type> sptkeys; // Keys of multimap.

    // Loop over TPCs.
    for (auto const& tpcid : geom->Iterate<geo::TPCID>()) {

      auto const [cstat, tpc] = std::make_tuple(tpcid.Cryostat, tpcid.TPC);
      // Sort maps in increasing order of number of hits.
      // This is so that we can do the outer loops over hits
      // over the views with fewer hits.
      //
      // If config parameter PreferColl is true, treat the colleciton
      // plane as if it had the most hits, regardless of how many
      // hits it actually has.  This will force space points to be
      // filtered and merged with respect to the collection plane
      // wires.  It will also force space points to be sorted by
      // collection plane wire.

      int nplane = hitmap[cstat][tpc].size();
      std::vector<int> index(nplane);

      for (int i = 0; i < nplane; ++i)
        index[i] = i;

      for (int i = 0; i < nplane - 1; ++i) {

        for (int j = i + 1; j < nplane; ++j) {
          bool icoll =
            fPreferColl && geom->SignalType(geo::PlaneID(tpcid, index[i])) == geo::kCollection;
          bool jcoll =
            fPreferColl && geom->SignalType(geo::PlaneID(tpcid, index[j])) == geo::kCollection;
          if ((hitmap[cstat][tpc][index[i]].size() > hitmap[cstat][tpc][index[j]].size() &&
               !jcoll) ||
              icoll) {
            int temp = index[i];
            index[i] = index[j];
            index[j] = temp;
          }
        }
      } // end loop over i

      // how many views with hits?
      // This will allow for the special case where we might have only 2 planes of information and
      // still want space points even if a three plane TPC
      std::vector<std::multimap<unsigned int, art::Ptr<recob::Hit>>>& hitsByPlaneVec =
        hitmap[cstat][tpc];
      int nViewsWithHits(0);

      for (int i = 0; i < nplane; i++) {
        if (hitsByPlaneVec[index[i]].size() > 0) nViewsWithHits++;
      }

      // If two-view space points are allowed, make a double loop
      // over hits and produce space points for compatible hit-pairs.

      if ((nViewsWithHits == 2 || nplane == 2) && fMinViews <= 2) {

        // Loop over pairs of views.
        for (int i = 0; i < nplane - 1; ++i) {
          unsigned int plane1 = index[i];

          if (hitmap[cstat][tpc][plane1].empty()) continue;

          for (int j = i + 1; j < nplane; ++j) {
            unsigned int plane2 = index[j];

            if (hitmap[cstat][tpc][plane2].empty()) continue;

            // Get angle, pitch, and offset of plane2 wires.
            geo::PlaneID const plane2_id{tpcid, plane2};
            const geo::WireGeo& wgeo2 = geom->Plane(plane2_id).Wire(0);
            double const hl2 = wgeo2.HalfL();
            auto const xyz21 = wgeo2.GetStart();
            auto const xyz22 = wgeo2.GetEnd();
            double s2 = (xyz22.Y() - xyz21.Y()) / (2. * hl2);
            double c2 = (xyz22.Z() - xyz21.Z()) / (2. * hl2);
            double dist2 = -xyz21.Y() * c2 + xyz21.Z() * s2;
            double pitch2 = geom->WirePitch(plane2_id);

            if (!fPreferColl &&
                hitmap[cstat][tpc][plane1].size() > hitmap[cstat][tpc][plane2].size())
              throw cet::exception("SpacePointAlg")
                << "makeSpacePoints(): hitmaps with incompatible size\n";

            // Loop over pairs of hits.

            art::PtrVector<recob::Hit> hitvec;
            hitvec.reserve(2);

            for (std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator ihit1 =
                   hitmap[cstat][tpc][plane1].begin();
                 ihit1 != hitmap[cstat][tpc][plane1].end();
                 ++ihit1) {

              const art::Ptr<recob::Hit>& phit1 = ihit1->second;
              geo::WireID phit1WireID = phit1->WireID();
              const geo::WireGeo& wgeo = geom->WireIDToWireGeo(phit1WireID);

              // Get endpoint coordinates of this wire.
              // (kept as assertions for performance reasons)
              assert(phit1WireID.Cryostat == cstat);
              assert(phit1WireID.TPC == tpc);
              assert(phit1WireID.Plane == plane1);
              auto const xyz1 = wgeo.GetStart();
              auto const xyz2 = wgeo.GetEnd();

              // Find the plane2 wire numbers corresponding to the endpoints.

              double wire21 = (-xyz1.Y() * c2 + xyz1.Z() * s2 - dist2) / pitch2;
              double wire22 = (-xyz2.Y() * c2 + xyz2.Z() * s2 - dist2) / pitch2;

              int wmin = std::max(0., std::min(wire21, wire22));
              int wmax = std::max(0., std::max(wire21, wire22) + 1.);

              std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator
                ihit2 = hitmap[cstat][tpc][plane2].lower_bound(wmin),
                ihit2end = hitmap[cstat][tpc][plane2].upper_bound(wmax);

              for (; ihit2 != ihit2end; ++ihit2) {

                const art::Ptr<recob::Hit>& phit2 = ihit2->second;

                // Check current pair of hits for compatibility.
                // By construction, hits should always have compatible views
                // and times, but may not have compatible mc information.

                hitvec.clear();
                hitvec.push_back(phit1);
                hitvec.push_back(phit2);
                bool ok = compatible(detProp, hitvec, useMC);
                if (ok) {

                  // Add a space point.

                  ++n2;

                  // make a dummy vector of recob::SpacePoints
                  // as we are filtering or merging and don't want to
                  // add the created SpacePoint to the final collection just yet
                  // This dummy vector will hold just one recob::SpacePoint,
                  // which will go into the multimap and then the vector
                  // will go out of scope.

                  std::vector<recob::SpacePoint> sptv;
                  fillSpacePoint(detProp, hitvec, sptv, sptmap.size());
                  sptkey_type key = &*phit2;
                  sptmap.insert(std::pair<sptkey_type, recob::SpacePoint>(key, sptv.back()));
                  sptkeys.insert(key);
                }
              }
            }
          }
        }
      } // end if fMinViews <= 2

      // If three-view space points are allowed, make a triple loop
      // over hits and produce space points for compatible triplets.

      if (nplane >= 3 && fMinViews <= 3) {

        // Loop over triplets of hits.

        art::PtrVector<recob::Hit> hitvec;
        hitvec.reserve(3);

        unsigned int plane1 = index[0];
        unsigned int plane2 = index[1];
        unsigned int plane3 = index[2];

        // Get angle, pitch, and offset of plane1 wires.

        geo::PlaneID const plane1_id{tpcid, plane1};
        const geo::WireGeo& wgeo1 = geom->Plane(plane1_id).Wire(0);
        double const hl1 = wgeo1.HalfL();
        auto const xyz11 = wgeo1.GetStart();
        auto const xyz12 = wgeo1.GetEnd();
        double s1 = (xyz12.Y() - xyz11.Y()) / (2. * hl1);
        double c1 = (xyz12.Z() - xyz11.Z()) / (2. * hl1);
        double dist1 = -xyz11.Y() * c1 + xyz11.Z() * s1;
        double pitch1 = geom->WirePitch(plane1_id);
        const double TicksOffset1 = detProp.GetXTicksOffset(plane1_id);

        // Get angle, pitch, and offset of plane2 wires.

        geo::PlaneID const plane2_id{tpcid, plane2};
        const geo::WireGeo& wgeo2 = geom->Plane(plane2_id).Wire(0);
        double const hl2 = wgeo2.HalfL();
        auto const xyz21 = wgeo2.GetStart();
        auto const xyz22 = wgeo2.GetEnd();
        double s2 = (xyz22.Y() - xyz21.Y()) / (2. * hl2);
        double c2 = (xyz22.Z() - xyz21.Z()) / (2. * hl2);
        double dist2 = -xyz21.Y() * c2 + xyz21.Z() * s2;
        double pitch2 = geom->WirePitch(plane2_id);
        const double TicksOffset2 = detProp.GetXTicksOffset(plane2_id);

        // Get angle, pitch, and offset of plane3 wires.

        geo::PlaneID const plane3_id{tpcid, plane3};
        const geo::WireGeo& wgeo3 = geom->Plane(plane3_id).Wire(0);
        double const hl3 = wgeo3.HalfL();
        auto const xyz31 = wgeo3.GetStart();
        auto const xyz32 = wgeo3.GetEnd();
        double s3 = (xyz32.Y() - xyz31.Y()) / (2. * hl3);
        double c3 = (xyz32.Z() - xyz31.Z()) / (2. * hl3);
        double dist3 = -xyz31.Y() * c3 + xyz31.Z() * s3;
        double pitch3 = geom->WirePitch(plane3_id);
        const double TicksOffset3 = detProp.GetXTicksOffset(plane3_id);

        // Get sine of angle differences.

        double s12 = s1 * c2 - s2 * c1; // sin(theta1 - theta2).
        double s23 = s2 * c3 - s3 * c2; // sin(theta2 - theta3).
        double s31 = s3 * c1 - s1 * c3; // sin(theta3 - theta1).

        // Loop over hits in plane1.

        std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator
          ihit1 = hitmap[cstat][tpc][plane1].begin(),
          ihit1end = hitmap[cstat][tpc][plane1].end();
        for (; ihit1 != ihit1end; ++ihit1) {

          unsigned int wire1 = ihit1->first;
          const art::Ptr<recob::Hit>& phit1 = ihit1->second;
          geo::WireID phit1WireID = phit1->WireID();
          const geo::WireGeo& wgeo = geom->WireIDToWireGeo(phit1WireID);

          // Get endpoint coordinates of this wire from plane1.
          // (kept as assertions for performance reasons)
          assert(phit1WireID.Cryostat == cstat);
          assert(phit1WireID.TPC == tpc);
          assert(phit1WireID.Plane == plane1);
          assert(phit1WireID.Wire == wire1);
          auto const xyz1 = wgeo.GetStart();
          auto const xyz2 = wgeo.GetEnd();

          // Get corrected time and oblique coordinate of first hit.

          double t1 = phit1->PeakTime() - TicksOffset1;
          double u1 = wire1 * pitch1 + dist1;

          // Find the plane2 wire numbers corresponding to the endpoints.

          double wire21 = (-xyz1.Y() * c2 + xyz1.Z() * s2 - dist2) / pitch2;
          double wire22 = (-xyz2.Y() * c2 + xyz2.Z() * s2 - dist2) / pitch2;

          int wmin = std::max(0., std::min(wire21, wire22));
          int wmax = std::max(0., std::max(wire21, wire22) + 1.);

          std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator
            ihit2 = hitmap[cstat][tpc][plane2].lower_bound(wmin),
            ihit2end = hitmap[cstat][tpc][plane2].upper_bound(wmax);

          for (; ihit2 != ihit2end; ++ihit2) {

            int wire2 = ihit2->first;
            const art::Ptr<recob::Hit>& phit2 = ihit2->second;

            // Get corrected time of second hit.

            double t2 = phit2->PeakTime() - TicksOffset2;

            // Check maximum time difference with first hit.

            bool dt12ok = std::abs(t1 - t2) <= fMaxDT;
            if (dt12ok) {

              // Test first two hits for compatibility before looping
              // over third hit.

              hitvec.clear();
              hitvec.push_back(phit1);
              hitvec.push_back(phit2);
              bool h12ok = compatible(detProp, hitvec, useMC);
              if (h12ok) {

                // Get oblique coordinate of second hit.

                double u2 = wire2 * pitch2 + dist2;

                // Predict plane3 oblique coordinate and wire number.

                double u3pred = (-u1 * s23 - u2 * s31) / s12;
                double w3pred = (u3pred - dist3) / pitch3;
                double w3delta = std::abs(fMaxS / (s12 * pitch3));
                int w3min = std::max(0., std::ceil(w3pred - w3delta));
                int w3max = std::max(0., std::floor(w3pred + w3delta));

                std::map<unsigned int, art::Ptr<recob::Hit>>::const_iterator
                  ihit3 = hitmap[cstat][tpc][plane3].lower_bound(w3min),
                  ihit3end = hitmap[cstat][tpc][plane3].upper_bound(w3max);

                for (; ihit3 != ihit3end; ++ihit3) {

                  int wire3 = ihit3->first;
                  const art::Ptr<recob::Hit>& phit3 = ihit3->second;

                  // Get corrected time of third hit.

                  double t3 = phit3->PeakTime() - TicksOffset3;

                  // Check time difference of third hit compared to first two hits.

                  bool dt123ok = std::abs(t1 - t3) <= fMaxDT && std::abs(t2 - t3) <= fMaxDT;
                  if (dt123ok) {

                    // Get oblique coordinate of third hit and check spatial separation.

                    double u3 = wire3 * pitch3 + dist3;
                    double S = s23 * u1 + s31 * u2 + s12 * u3;
                    bool sok = std::abs(S) <= fMaxS;
                    if (sok) {

                      // Test triplet for compatibility.

                      hitvec.clear();
                      hitvec.push_back(phit1);
                      hitvec.push_back(phit2);
                      hitvec.push_back(phit3);
                      bool h123ok = compatible(detProp, hitvec, useMC);
                      if (h123ok) {

                        // Add a space point.

                        ++n3;

                        // make a dummy vector of recob::SpacePoints
                        // as we are filtering or merging and don't want to
                        // add the created SpacePoint to the final collection just yet
                        // This dummy vector will hold just one recob::SpacePoint,
                        // which will go into the multimap and then the vector
                        // will go out of scope.

                        std::vector<recob::SpacePoint> sptv;
                        fillSpacePoint(detProp, hitvec, sptv, sptmap.size() - 1);
                        sptkey_type key = &*phit3;
                        sptmap.insert(std::pair<sptkey_type, recob::SpacePoint>(key, sptv.back()));
                        sptkeys.insert(key);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } // end if fMinViews <= 3

      // Do Filtering.

      if (fFilter) {

        // Transfer (some) space points from sptmap to spts.

        spts.reserve(spts.size() + sptkeys.size());

        // Loop over keys of space point map.
        // Space points that have the same key are candidates for filtering.

        for (std::set<sptkey_type>::const_iterator i = sptkeys.begin(); i != sptkeys.end(); ++i) {
          sptkey_type key = *i;

          // Loop over space points corresponding to the current key.
          // Choose the single best space point from among this group.

          double best_chisq = 0.;
          const recob::SpacePoint* best_spt = 0;

          for (std::multimap<sptkey_type, recob::SpacePoint>::const_iterator j =
                 sptmap.lower_bound(key);
               j != sptmap.upper_bound(key);
               ++j) {
            const recob::SpacePoint& spt = j->second;
            if (best_spt == 0 || spt.Chisq() < best_chisq) {
              best_spt = &spt;
              best_chisq = spt.Chisq();
            }
          }

          // Transfer best filtered space point to result vector.

          if (!best_spt)
            throw cet::exception("SpacePointAlg") << "makeSpacePoints(): no best point\n";
          spts.push_back(*best_spt);
          if (fMinViews <= 2)
            ++n2filt;
          else
            ++n3filt;
        }
      } // end if filtering

      // Do merging.

      else if (fMerge) {

        // Transfer merged space points from sptmap to spts.

        spts.reserve(spts.size() + sptkeys.size());

        // Loop over keys of space point map.
        // Space points that have the same key are candidates for merging.

        for (std::set<sptkey_type>::const_iterator i = sptkeys.begin(); i != sptkeys.end(); ++i) {
          sptkey_type key = *i;

          // Loop over space points corresponding to the current key.
          // Make a collection of hits that is the union of the hits
          // from each candidate space point.

          std::multimap<sptkey_type, recob::SpacePoint>::const_iterator jSPT =
                                                                          sptmap.lower_bound(key),
                                                                        jSPTend =
                                                                          sptmap.upper_bound(key);

          art::PtrVector<recob::Hit> merged_hits;
          for (; jSPT != jSPTend; ++jSPT) {
            const recob::SpacePoint& spt = jSPT->second;

            // Loop over hits from this space points.
            // Add each hit to the collection of all hits.

            const art::PtrVector<recob::Hit>& spt_hits = getAssociatedHits(spt);
            merged_hits.reserve(merged_hits.size() +
                                spt_hits.size()); // better than nothing, but not ideal
            for (art::PtrVector<recob::Hit>::const_iterator k = spt_hits.begin();
                 k != spt_hits.end();
                 ++k) {
              const art::Ptr<recob::Hit>& hit = *k;
              merged_hits.push_back(hit);
            }
          }

          // Remove duplicates.

          std::sort(merged_hits.begin(), merged_hits.end());
          art::PtrVector<recob::Hit>::iterator it =
            std::unique(merged_hits.begin(), merged_hits.end());
          merged_hits.erase(it, merged_hits.end());

          // Construct a complex space points using merged hits.

          fillComplexSpacePoint(detProp, merged_hits, spts, sptmap.size() + spts.size() - 1);

          if (fMinViews <= 2)
            ++n2filt;
          else
            ++n3filt;
        }
      } // end if merging

      // No filter, no merge.

      else {

        // Transfer all space points from sptmap to spts.

        spts.reserve(spts.size() + sptkeys.size());

        // Loop over space points.

        for (std::multimap<sptkey_type, recob::SpacePoint>::const_iterator j = sptmap.begin();
             j != sptmap.end();
             ++j) {
          const recob::SpacePoint& spt = j->second;
          spts.push_back(spt);
        }

        // Update statistics.

        n2filt = n2;
        n3filt = n3;
      }
    } // end loop over tpcs

    if (mf::isDebugEnabled()) {
      debug << "\n2-hit space points = " << n2 << "\n"
            << "3-hit space points = " << n3 << "\n"
            << "2-hit filtered/merged space points = " << n2filt << "\n"
            << "3-hit filtered/merged space points = " << n3filt;
    } // if debug
  }

  //----------------------------------------------------------------------
  // Get hits associated with a particular space point, based on most recent
  // call of any make*SpacePoints method.
  const art::PtrVector<recob::Hit>& SpacePointAlg::getAssociatedHits(
    const recob::SpacePoint& spt) const
  {
    // It is an error if no hits are associated with this space point (throw exception).

    std::map<int, art::PtrVector<recob::Hit>>::const_iterator it = fSptHitMap.find(spt.ID());
    if (it == fSptHitMap.end()) {
      mf::LogWarning("SpacePointAlg")
        << "Looking for ID " << spt.ID() << " from " << fSptHitMap.size() << std::endl;
      throw cet::exception("SpacePointAlg") << "No Hits associated with space point.\n";
    }
    return (*it).second;
  }

} // end namespace trkf
