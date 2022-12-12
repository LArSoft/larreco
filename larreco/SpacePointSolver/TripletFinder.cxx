#include "larreco/SpacePointSolver/TripletFinder.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TVector3.h"

#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

namespace reco3d {
  // -------------------------------------------------------------------------
  TripletFinder::TripletFinder(const detinfo::DetectorPropertiesData& detProp,
                               const std::vector<art::Ptr<recob::Hit>>& xhits,
                               const std::vector<art::Ptr<recob::Hit>>& uhits,
                               const std::vector<art::Ptr<recob::Hit>>& vhits,
                               const std::vector<raw::ChannelID_t>& xbad,
                               const std::vector<raw::ChannelID_t>& ubad,
                               const std::vector<raw::ChannelID_t>& vbad,
                               double distThresh,
                               double distThreshDrift,
                               double xhitOffset)
    : geom(art::ServiceHandle<geo::Geometry const>()->provider())
    , channelMapAlg{art::ServiceHandle<geo::ExptGeoHelperInterface const>()->ChannelMapAlgPtr()}
    , fDistThresh(distThresh)
    , fDistThreshDrift(distThreshDrift)
    , fXHitOffset(xhitOffset)
  {
    FillHitMap(detProp, xhits, fX_by_tpc);
    FillHitMap(detProp, uhits, fU_by_tpc);
    FillHitMap(detProp, vhits, fV_by_tpc);

    FillBadMap(xbad, fXbad_by_tpc);
    FillBadMap(ubad, fUbad_by_tpc);
    FillBadMap(vbad, fVbad_by_tpc);
  }

  // -------------------------------------------------------------------------
  void TripletFinder::FillHitMap(const detinfo::DetectorPropertiesData& detProp,
                                 const std::vector<art::Ptr<recob::Hit>>& hits,
                                 std::map<geo::TPCID, std::vector<HitOrChan>>& out)
  {
    for (const art::Ptr<recob::Hit>& hit : hits) {
      for (geo::TPCID tpc : channelMapAlg->ROPtoTPCs(channelMapAlg->ChannelToROP(hit->Channel()))) {
        double xpos = 0;
        for (geo::WireID wire : channelMapAlg->ChannelToWire(hit->Channel())) {
          if (geo::TPCID(wire) == tpc) {
            xpos = detProp.ConvertTicksToX(hit->PeakTime(), wire);
            if (channelMapAlg->SignalType(wire) == geo::kCollection) xpos += fXHitOffset;
          }
        }

        out[tpc].emplace_back(hit.get(), xpos);
      }
    }
    for (auto& it : out)
      std::sort(it.second.begin(), it.second.end(), [](auto a, auto b) { return a.xpos < b.xpos; });
  }

  // -------------------------------------------------------------------------
  void TripletFinder::FillBadMap(const std::vector<raw::ChannelID_t>& bads,
                                 std::map<geo::TPCID, std::vector<raw::ChannelID_t>>& out)
  {
    for (raw::ChannelID_t chan : bads) {
      for (geo::TPCID tpc : channelMapAlg->ROPtoTPCs(channelMapAlg->ChannelToROP(chan))) {
        out[tpc].push_back(chan);
      }
    }
  }

  // -------------------------------------------------------------------------
  class IntersectionCache {
  public:
    IntersectionCache(geo::GeometryCore const* geo, geo::ChannelMapAlg const* cmAlg, geo::TPCID tpc)
      : geom(geo), channelMapAlg{cmAlg}, fTPC(tpc)
    {}

    std::optional<geo::WireIDIntersection> const& operator()(raw::ChannelID_t a, raw::ChannelID_t b)
    {
      auto const key = std::make_pair(a, b);

      auto it = fPtMap.find(key);
      if (it != fPtMap.end()) { return it->second; }

      return fPtMap.try_emplace(key, ISect(a, b)).first->second;
    }

  private:
    std::optional<geo::WireIDIntersection> ISect(raw::ChannelID_t chanA,
                                                 raw::ChannelID_t chanB) const
    {
      for (geo::WireID awire : channelMapAlg->ChannelToWire(chanA)) {
        if (geo::TPCID(awire) != fTPC) continue;
        for (geo::WireID bwire : channelMapAlg->ChannelToWire(chanB)) {
          if (geo::TPCID(bwire) != fTPC) continue;

          if (auto pt = geom->WireIDsIntersect(awire, bwire)) return pt;
        }
      }

      return std::nullopt;
    }

    const geo::GeometryCore* geom;
    const geo::ChannelMapAlg* channelMapAlg;

    using key_t = std::pair<raw::ChannelID_t, raw::ChannelID_t>;
    std::map<key_t, std::optional<geo::WireIDIntersection>> fPtMap;

    geo::TPCID fTPC;
  };

  // -------------------------------------------------------------------------
  bool TripletFinder::CloseDrift(double xa, double xb) const
  {
    return fabs(xa - xb) < fDistThreshDrift;
  }

  // -------------------------------------------------------------------------
  bool TripletFinder::CloseSpace(geo::WireIDIntersection ra, geo::WireIDIntersection rb) const
  {
    const TVector3 pa(ra.y, ra.z, 0);
    const TVector3 pb(rb.y, rb.z, 0);

    return (pa - pb).Mag() < fDistThresh;
  }

  bool LessThanXHit(const ChannelDoublet& a, const ChannelDoublet& b)
  {
    // Make sure the bad hits get sorted too
    if (a.a.hit == 0 && b.a.hit == 0) return a.a.chan < b.a.chan;
    // But mostly just order the real hits (in some arbitrary order)
    return a.a.hit < b.a.hit;
  }

  bool SameXHit(const ChannelDoublet& a, const ChannelDoublet& b)
  {
    if (a.a.hit == 0 && b.a.hit == 0) return a.a.chan == b.a.chan;
    return a.a.hit == b.a.hit;
  }

  // -------------------------------------------------------------------------
  std::vector<HitTriplet> TripletFinder::Triplets()
  {
    std::vector<HitTriplet> ret;

    for (const auto& it : fX_by_tpc) {
      const geo::TPCID& tpc = it.first;

      std::vector<ChannelDoublet> xus = DoubletsXU(tpc);
      std::vector<ChannelDoublet> xvs = DoubletsXV(tpc);

      // Cache to prevent repeating the same questions
      IntersectionCache isectUV(geom, channelMapAlg, tpc);

      // For the efficient looping below to work we need to sort the doublet
      // lists so the X hits occur in the same order.
      std::sort(xus.begin(), xus.end(), LessThanXHit);
      std::sort(xvs.begin(), xvs.end(), LessThanXHit);

      auto xvit_begin = xvs.begin();

      int nxuv = 0;
      for (const ChannelDoublet& xu : xus) {
        const HitOrChan& x = xu.a;
        const HitOrChan& u = xu.b;

        // Catch up until we're looking at the same X hit in XV
        while (xvit_begin != xvs.end() && LessThanXHit(*xvit_begin, xu))
          ++xvit_begin;

        // Loop through all those matching hits
        for (auto xvit = xvit_begin; xvit != xvs.end() && SameXHit(*xvit, xu); ++xvit) {
          const HitOrChan& v = xvit->b;

          // Only allow one bad channel per triplet
          if (!x.hit && !u.hit) continue;
          if (!x.hit && !v.hit) continue;
          if (!u.hit && !v.hit) continue;

          if (u.hit && v.hit && !CloseDrift(u.xpos, v.xpos)) continue;

          auto maybe_ptUV = isectUV(u.chan, v.chan);
          if (!maybe_ptUV) continue;

          auto const& ptUV = *maybe_ptUV;
          if (!CloseSpace(xu.pt, xvit->pt) || !CloseSpace(xu.pt, ptUV) ||
              !CloseSpace(xvit->pt, ptUV))
            continue;

          double xavg = 0;
          int nx = 0;
          if (x.hit) {
            xavg += x.xpos;
            ++nx;
          }
          if (u.hit) {
            xavg += u.xpos;
            ++nx;
          }
          if (v.hit) {
            xavg += v.xpos;
            ++nx;
          }
          xavg /= nx;

          const XYZ pt{
            xavg, (xu.pt.y + xvit->pt.y + ptUV.y) / 3, (xu.pt.z + xvit->pt.z + ptUV.z) / 3};

          ret.emplace_back(HitTriplet{x.hit, u.hit, v.hit, pt});
          ++nxuv;
        } // end for xv
      }   // end for xu

      std::cout << tpc << " " << xus.size() << " XUs and " << xvs.size() << " XVs -> " << nxuv
                << " XUVs" << std::endl;

    } // end for tpc

    std::cout << ret.size() << " XUVs total" << std::endl;

    return ret;
  }

  // -------------------------------------------------------------------------
  std::vector<HitTriplet> TripletFinder::TripletsTwoView()
  {
    std::vector<HitTriplet> ret;

    for (const auto& it : fX_by_tpc) {
      const geo::TPCID& tpc = it.first;

      std::vector<ChannelDoublet> xus = DoubletsXU(tpc);

      for (const ChannelDoublet& xu : xus) {
        const HitOrChan& x = xu.a;
        const HitOrChan& u = xu.b;

        double xavg = x.xpos;
        int nx = 1;
        if (u.hit) {
          xavg += u.xpos;
          ++nx;
        }
        xavg /= nx;

        const XYZ pt{xavg, xu.pt.y, xu.pt.z};

        ret.emplace_back(HitTriplet{x.hit, u.hit, 0, pt});
      } // end for xu
    }   // end for tpc

    std::cout << ret.size() << " XUs total" << std::endl;

    return ret;
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelDoublet> TripletFinder::DoubletsXU(geo::TPCID tpc)
  {
    std::vector<ChannelDoublet> ret =
      DoubletHelper(tpc, fX_by_tpc[tpc], fU_by_tpc[tpc], fUbad_by_tpc[tpc]);

    // Find X(bad)+U(good) doublets, have to flip them for the final result
    for (auto it : DoubletHelper(tpc, fU_by_tpc[tpc], {}, fXbad_by_tpc[tpc])) {
      ret.push_back({it.b, it.a, it.pt});
    }

    return ret;
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelDoublet> TripletFinder::DoubletsXV(geo::TPCID tpc)
  {
    std::vector<ChannelDoublet> ret =
      DoubletHelper(tpc, fX_by_tpc[tpc], fV_by_tpc[tpc], fVbad_by_tpc[tpc]);

    // Find X(bad)+V(good) doublets, have to flip them for the final result
    for (auto it : DoubletHelper(tpc, fV_by_tpc[tpc], {}, fXbad_by_tpc[tpc])) {
      ret.push_back({it.b, it.a, it.pt});
    }

    return ret;
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelDoublet> TripletFinder::DoubletHelper(
    geo::TPCID tpc,
    const std::vector<HitOrChan>& ahits,
    const std::vector<HitOrChan>& bhits,
    const std::vector<raw::ChannelID_t>& bbads) const
  {
    std::vector<ChannelDoublet> ret;

    IntersectionCache isect(geom, channelMapAlg, tpc);

    auto b_begin = bhits.begin();

    for (const HitOrChan& a : ahits) {
      // Bad channels are easy because there's no timing constraint
      for (raw::ChannelID_t b : bbads) {
        if (auto pt = isect(a.chan, b)) { ret.emplace_back(a, b, *pt); }
      }

      while (b_begin != bhits.end() && b_begin->xpos < a.xpos && !CloseDrift(b_begin->xpos, a.xpos))
        ++b_begin;

      for (auto bit = b_begin; bit != bhits.end(); ++bit) {
        const HitOrChan& b = *bit;

        if (b.xpos > a.xpos && !CloseDrift(b.xpos, a.xpos)) break;

        auto pt = isect(a.chan, b.chan);
        if (!pt) continue;

        ret.emplace_back(a, b, *pt);
      } // end for b
    }   // end for a

    return ret;
  }
}
