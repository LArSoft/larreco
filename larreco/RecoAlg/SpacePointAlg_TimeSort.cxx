/*!
 * Title:   SpacePointAlg_TimeSort class
 * Author:  wketchum@lanl.gov
 * Inputs:  std::vector<recob::Hit> (one for each plane)
 * Outputs: std::vector<recob::SpacePoint>
 *
 * Description:
 * This space point algorithm tries to improve speed by
 * (1) keeping hit collections distinct among planes;
 * (2) sorting hit collections by time; and,
 * (3) having a lookup table for (y,z) coordinate positions.
 * The original use case for this code was with the TTHitFinder,
 * which produces an incredibly large number of hits per plane,
 * making some sorted space point alg more attractive.
 *
 * This code is totally microboone specific, btw.
 */
#include "larreco/RecoAlg/SpacePointAlg_TimeSort.h"

#include <math.h>

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

//boost includes
#include "boost/multi_array.hpp"

namespace {
  inline bool HitTimeComparison(art::Ptr<recob::Hit> const& a, art::Ptr<recob::Hit> const& b)
  {
    return a->PeakTime() < b->PeakTime();
  }
}

namespace sppt {

  //-------------------------------------------------
  SpacePointAlg_TimeSort::SpacePointAlg_TimeSort(fhicl::ParameterSet const& pset)
  {
    fTimeDiffMax = pset.get<float>("TimeDiffMax");
    fZDiffMax = pset.get<float>("ZDiffMax");
    fYDiffMax = pset.get<float>("YDiffMax");

    //enforce a minimum time diff
    if (fTimeDiffMax < 0) {
      throw cet::exception("SpacePointAlg_TimeSort")
        << "Time difference must be greater than zero.";
    }
    if (fZDiffMax < 0) {
      throw cet::exception("SpacePointAlg_TimeSort")
        << "Z-coordinate difference must be greater than zero.";
    }
    if (fYDiffMax < 0) {
      throw cet::exception("SpacePointAlg_TimeSort")
        << "Y-coordinate difference must be greater than zero.";
    }
  }

  //-------------------------------------------------
  void SpacePointAlg_TimeSort::setTimeOffsets(detinfo::DetectorPropertiesData const& detProp)
  {
    TIME_OFFSET_U = -1 * detProp.GetXTicksOffset(geo::View_t::kU, 0, 0);
    TIME_OFFSET_V = -1 * detProp.GetXTicksOffset(geo::View_t::kV, 0, 0);
    TIME_OFFSET_Y = -1 * detProp.GetXTicksOffset(geo::View_t::kZ, 0, 0);
    TICKS_TO_X = detProp.GetXTicksCoefficient();

    TIME_OFFSET_SET = true;
  }

  //-------------------------------------------------
  void SpacePointAlg_TimeSort::fillCoordinatesArrays()
  {
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    constexpr geo::TPCID tpcid{0, 0};
    constexpr geo::PlaneID uplane_id{tpcid, geo::View_t::kU};
    constexpr geo::PlaneID vplane_id{tpcid, geo::View_t::kV};
    constexpr geo::PlaneID zplane_id{tpcid, geo::View_t::kZ};

    unsigned int nwires_u = wireReadoutGeom.Nwires(uplane_id);
    unsigned int nwires_v = wireReadoutGeom.Nwires(vplane_id);
    unsigned int nwires_y = wireReadoutGeom.Nwires(zplane_id);

    coordinates_UV_y.resize(boost::extents[nwires_v][nwires_u]);
    coordinates_UV_z.resize(boost::extents[nwires_v][nwires_u]);
    coordinates_UY_y.resize(boost::extents[nwires_y][nwires_u]);
    coordinates_UY_z.resize(boost::extents[nwires_y][nwires_u]);
    for (unsigned int iu = 0; iu < nwires_u; iu++) {
      geo::WireID const uplane_wire_id{uplane_id, iu};
      for (unsigned int iv = 0; iv < nwires_v; iv++) {
        auto intersection =
          wireReadoutGeom.WireIDsIntersect(uplane_wire_id, geo::WireID{vplane_id, iv})
            .value_or(geo::WireIDIntersection::invalid());
        coordinates_UV_y[iv][iu] = intersection.y;
        coordinates_UV_z[iv][iu] = intersection.z;
      }
      for (unsigned int iy = 0; iy < nwires_y; iy++) {
        auto intersection =
          wireReadoutGeom.WireIDsIntersect(uplane_wire_id, geo::WireID{zplane_id, iy})
            .value_or(geo::WireIDIntersection::invalid());
        coordinates_UY_y[iy][iu] = intersection.y;
        coordinates_UY_z[iy][iu] = intersection.z;
      }
    }

    COORDINATES_FILLED = true;
  }

  //-------------------------------------------------
  void SpacePointAlg_TimeSort::createSpacePoints(
    detinfo::DetectorPropertiesData const& detProp,
    std::vector<art::Ptr<recob::Hit>>& hitVec_U,
    std::vector<art::Ptr<recob::Hit>>& hitVec_V,
    std::vector<art::Ptr<recob::Hit>>& hitVec_Y,
    std::unique_ptr<std::vector<recob::SpacePoint>>& spptCollection,
    std::unique_ptr<std::vector<std::vector<art::Ptr<recob::Hit>>>>& spptAssociatedHits)
  {

    if (!TIME_OFFSET_SET) {
      mf::LogWarning("SpacePointAlg_TimeSort")
        << "Time offsets not set before createSpacePoints call!"
        << "\nYou should call SpacePointAlg_TimeSort::setTimeOffsets() in beginRun()!"
        << "\nWill be set now, but you should modify your code!";
      setTimeOffsets(detProp);
    }
    if (!COORDINATES_FILLED) {
      mf::LogWarning("SpacePointAlg_TimeSort")
        << "Coordinate arrays not filled before createSpacePoints call!"
        << "\nYou should call SpacePointAlg_TimeSort::fillCoordinateArrays() in beginRun()!"
        << "\nWill be filled now, but you should modify your code!";
      fillCoordinatesArrays();
    }

    //sort the hits by the time
    sortHitsByTime(hitVec_U);
    sortHitsByTime(hitVec_V);
    sortHitsByTime(hitVec_Y);

    MF_LOG_DEBUG("SpacePointAlg_TimeSort")
      << "Sorted " << hitVec_U.size() << " u hits, " << hitVec_V.size() << " v hits, "
      << hitVec_Y.size() << " y hits.";

    //now, do the loop to search for like-timed hits across the three planes
    std::vector<art::Ptr<recob::Hit>>::iterator ihitu = hitVec_U.begin();
    std::vector<art::Ptr<recob::Hit>>::iterator ihitv = hitVec_V.begin();
    std::vector<art::Ptr<recob::Hit>>::iterator ihity = hitVec_Y.begin();
    std::vector<art::Ptr<recob::Hit>>::iterator ihitv_inner, ihity_inner;
    double time_hitu = (*ihitu)->PeakTime() + TIME_OFFSET_U;
    double time_hitv = (*ihitv)->PeakTime() + TIME_OFFSET_V;
    double time_hity = (*ihity)->PeakTime() + TIME_OFFSET_Y;
    double time_hitv_inner, time_hity_inner;
    while (ihitu != hitVec_U.end()) {
      time_hitu = (*ihitu)->PeakTime() + TIME_OFFSET_U;

      mf::LogInfo("SpacePointAlg_TimeSort")
        << "Hit times (u,v,y)=(" << time_hitu << "," << time_hitv << "," << time_hity << ")";

      //if time_hitu is too much bigger than time_hitv, need to advance hitv iterator
      while ((time_hitu - time_hitv) > fTimeDiffMax) {
        ihitv++;
        if (ihitv == hitVec_V.end()) break;
        time_hitv = (*ihitv)->PeakTime() + TIME_OFFSET_V;
      }
      if (ihitv == hitVec_V.end()) break;

      //same thing with time_hitu and time_hity
      while ((time_hitu - time_hity) > fTimeDiffMax) {
        ihity++;
        if (ihity == hitVec_Y.end()) break;
        time_hity = (*ihity)->PeakTime() + TIME_OFFSET_Y;
      }
      if (ihity == hitVec_Y.end()) break;

      //OK, now we know time_hitu <= time_hitv and time_hitu <= time_hity.
      //Next, check if time_hitu is near time_hitv and time_hit y. If not,
      //we have to increment ihitu.
      if (std::abs(time_hitu - time_hitv) > fTimeDiffMax ||
          std::abs(time_hitu - time_hity) > fTimeDiffMax) {
        ihitu++;
        continue;
      }

      //OK! Note we KNOW that these three match in time:
      //  -- time_hitu is within fTimeDiffMax of both time_hitv and time_hity; and
      //  -- time_hitu <= time_hitv AND time_hitu <=time_hity, so time_hitv and time_hity are near too

      mf::LogInfo("SpacePointAlg_TimeSort") << "Matching hit times (u,v,y)=(" << time_hitu << ","
                                            << time_hitv << "," << time_hity << ")";

      //Next thing to do, we need to loop over all possible 3-hit matches for our given u-hit.
      //We need new iterators in v and y at this location, and will loop over those
      ihitv_inner = ihitv;
      time_hitv_inner = (*ihitv_inner)->PeakTime() + TIME_OFFSET_V;
      ihity_inner = ihity;
      time_hity_inner = (*ihity_inner)->PeakTime() + TIME_OFFSET_Y;

      while (std::abs(time_hitu - time_hitv_inner) < fTimeDiffMax &&
             std::abs(time_hitu - time_hity_inner) < fTimeDiffMax) {

        unsigned int uwire = (*ihitu)->WireID().Wire;
        unsigned int vwire = (*ihitv_inner)->WireID().Wire;
        unsigned int ywire = (*ihity_inner)->WireID().Wire;

        mf::LogInfo("SpacePointAlg_TimeSort")
          << "(y,z) coordinate for uv/uy: (" << coordinates_UV_y[vwire][uwire] << ","
          << coordinates_UV_z[vwire][uwire] << ")/(" << coordinates_UY_y[ywire][uwire] << ","
          << coordinates_UY_z[ywire][uwire] << ")";

        if (std::abs(coordinates_UV_y[vwire][uwire] - coordinates_UY_y[ywire][uwire]) < fYDiffMax &&
            std::abs(coordinates_UV_z[vwire][uwire] - coordinates_UY_z[ywire][uwire]) < fZDiffMax) {

          double xyz[3];
          double xyz_err[6];
          //triangular error matrix:
          // | 0.         |
          // | 0.  0.     |
          // | 0.  0.  0. |

          //get average y and z, with errors
          xyz[1] = (coordinates_UV_y[vwire][uwire] + coordinates_UY_y[ywire][uwire]) * 0.5;
          xyz_err[2] = std::abs(coordinates_UV_y[vwire][uwire] - xyz[1]);
          xyz[2] = (coordinates_UV_z[vwire][uwire] + coordinates_UY_z[ywire][uwire]) * 0.5;
          xyz_err[5] = std::abs(coordinates_UV_z[vwire][uwire] - xyz[2]);

          double t_val = (time_hitu + time_hitv_inner + time_hity_inner) / 3.;
          double t_err = 0.5 * std::sqrt((time_hitu - t_val) * (time_hitu - t_val) +
                                         (time_hitv_inner - t_val) * (time_hitv_inner - t_val) +
                                         (time_hity_inner - t_val) * (time_hity_inner - t_val));
          xyz[0] = TICKS_TO_X * t_val;
          xyz_err[0] = TICKS_TO_X * t_err;

          //make space point to put on event
          recob::SpacePoint spt(xyz, xyz_err, 0., spptCollection->size());
          spptCollection->push_back(spt);

          //make association with hits
          std::vector<art::Ptr<recob::Hit>> myhits = {*ihitu, *ihitv_inner, *ihity_inner};
          spptAssociatedHits->push_back(myhits);
        }

        //now increment the v or y hit, whichever is smalles (closest to u hit) in time
        if (time_hitv_inner <= time_hity_inner) {
          ihitv_inner++;
          if (ihitv_inner == hitVec_V.end()) break;
          time_hitv_inner = (*ihitv_inner)->PeakTime() + TIME_OFFSET_V;
        }
        else {
          ihity_inner++;
          if (ihity_inner == hitVec_Y.end()) break;
          time_hity_inner = (*ihity_inner)->PeakTime() + TIME_OFFSET_Y;
        }
      }

      ihitu++;
    } // end while looping over u hits

    MF_LOG_DEBUG("SpacePointAlg_TimeSort")
      << "Finished with " << spptCollection->size() << " spacepoints.";

  } //end createSpacePoints

  //-------------------------------------------------
  void SpacePointAlg_TimeSort::sortHitsByTime(std::vector<art::Ptr<recob::Hit>>& hitVec) const
  {
    std::sort(hitVec.begin(), hitVec.end(), HitTimeComparison);
  }

} //end sppt namespace
