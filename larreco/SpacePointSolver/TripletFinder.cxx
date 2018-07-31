#include "larreco/SpacePointSolver/TripletFinder.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/Hit.h"

namespace reco3d
{
  // -------------------------------------------------------------------------
  TripletFinder::TripletFinder(const std::vector<art::Ptr<recob::Hit>>& xhits,
                               const std::vector<art::Ptr<recob::Hit>>& uhits,
                               const std::vector<art::Ptr<recob::Hit>>& vhits,
                               const std::vector<raw::ChannelID_t>& ubad,
                               const std::vector<raw::ChannelID_t>& vbad,
                               double distThresh, double distThreshDrift)
    : geom(art::ServiceHandle<geo::Geometry>()->provider()),
      detprop(art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider()),
      fDistThresh(distThresh),
      fDistThreshDrift(distThreshDrift)
  {
    std::vector<HitOrChan> xs, us, vs;
    for(const auto& x: xhits) xs.emplace_back(x->Channel(), x.get());
    for(const auto& u: uhits) us.emplace_back(u->Channel(), u.get());
    for(const auto& v: vhits) vs.emplace_back(v->Channel(), v.get());

    FillHitMap(xs, fX_by_tpc);
    FillHitMap(us, fU_by_tpc);
    FillHitMap(vs, fV_by_tpc);

    FillBadMap(ubad, fUbad_by_tpc);
    FillBadMap(vbad, fVbad_by_tpc);
  }

  // -------------------------------------------------------------------------
  bool sortByHitTime(const HitOrChan& a, const HitOrChan& b)
  {
    return a.hit->PeakTime() < b.hit->PeakTime();
  }

  // -------------------------------------------------------------------------
  void TripletFinder::
  FillHitMap(const std::vector<HitOrChan>& hits,
             std::map<geo::TPCID, std::vector<HitOrChan>>& out)
  {
    for(const HitOrChan& hit: hits){
      for(geo::TPCID tpc: geom->ROPtoTPCs(geom->ChannelToROP(hit.chan))){
        out[tpc].push_back(hit);
      }
    }
    for(auto& it: out) std::sort(it.second.begin(), it.second.end(), sortByHitTime);
  }

  // -------------------------------------------------------------------------
  void TripletFinder::
  FillBadMap(const std::vector<raw::ChannelID_t>& bads,
             std::map<geo::TPCID, std::vector<raw::ChannelID_t>>& out)
  {
    for(raw::ChannelID_t chan: bads){
      for(geo::TPCID tpc: geom->ROPtoTPCs(geom->ChannelToROP(chan))){
        out[tpc].push_back(chan);
      }
    }
  }

  // -------------------------------------------------------------------------
  bool sortByXHit(const ChannelDoublet& a, const ChannelDoublet& b)
  {
    return a.a.hit < b.a.hit;
  }

  // -------------------------------------------------------------------------
  class IntersectionCache
  {
  public:
    IntersectionCache(geo::TPCID tpc)
      : geom(art::ServiceHandle<geo::Geometry>()->provider()),
        fTPC(tpc)
    {
    }

    bool operator()(raw::ChannelID_t a, raw::ChannelID_t b,
                    geo::WireIDIntersection& pt)
    {
      const auto key = std::make_pair(a, b);

      auto it = fMap.find(key);
      if(it != fMap.end()){
        pt = fPtMap[key];
        return it->second;
      }

      const bool res = ISect(a, b, pt);
      fMap.insert({key, res});
      fPtMap.insert({key, pt});
      return res;
    }

  protected:
    bool ISect(raw::ChannelID_t chanA, raw::ChannelID_t chanB,
               geo::WireIDIntersection& pt) const
    {
      for(geo::WireID awire: geom->ChannelToWire(chanA)){
        if(geo::TPCID(awire) != fTPC) continue;
        for(geo::WireID bwire: geom->ChannelToWire(chanB)){
          if(geo::TPCID(bwire) != fTPC) continue;

          if(geom->WireIDsIntersect(awire, bwire, pt)) return true;
        }
      }

      return false;
    }

    const geo::GeometryCore* geom;

    std::map<std::pair<raw::ChannelID_t, raw::ChannelID_t>, bool> fMap;
    std::map<std::pair<raw::ChannelID_t, raw::ChannelID_t>, geo::WireIDIntersection> fPtMap;

    geo::TPCID fTPC;
  };

  // -------------------------------------------------------------------------
  double TripletFinder::HitToXPos(raw::ChannelID_t chan,
                                  double t,
                                  geo::TPCID tpc) const
  {
    for(geo::WireID w: geom->ChannelToWire(chan))
      if(geo::TPCID(w) == tpc)
        return detprop->ConvertTicksToX(t, w);
    // Not reached
    abort();
  }

  // -------------------------------------------------------------------------
  bool TripletFinder::CloseDrift(double xa, double xb) const
  {
    return fabs(xa-xb) < fDistThreshDrift;
  }

  // -------------------------------------------------------------------------
  bool TripletFinder::CloseTime(geo::TPCID tpc,
                                raw::ChannelID_t c1, raw::ChannelID_t c2,
                                double t1, double t2) const
  {
    return CloseDrift(HitToXPos(c1, t1, tpc),
                      HitToXPos(c2, t2, tpc));
  }

  // -------------------------------------------------------------------------
  bool TripletFinder::CloseSpace(geo::WireIDIntersection ra,
                                 geo::WireIDIntersection rb) const
  {
    const TVector3 pa(ra.y, ra.z, 0);
    const TVector3 pb(rb.y, rb.z, 0);

    return (pa-pb).Mag() < fDistThresh;
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelTriplet> TripletFinder::Triplets()
  {
    std::vector<ChannelTriplet> ret;

    for(const auto& it: fX_by_tpc){
      const geo::TPCID& tpc = it.first;

      std::vector<ChannelDoublet> xus = DoubletsXU(tpc);
      std::vector<ChannelDoublet> xvs = DoubletsXV(tpc);

      // Cache to prevent repeating the same questions
      IntersectionCache isectUV(tpc);

      // For the efficient looping below to work we need to sort the doublet
      // lists so the X hits occur in the same order.
      std::sort(xus.begin(), xus.end(), sortByXHit);
      std::sort(xvs.begin(), xvs.end(), sortByXHit);

      auto xvit_begin = xvs.begin();

      int nxuv = 0;
      for(const ChannelDoublet& xu: xus){
        const HitOrChan& x = xu.a;
        const HitOrChan& u = xu.b;

        // Catch up until we're looking at the same X hit in XV
        while(xvit_begin != xvs.end() && xvit_begin->a.hit < x.hit) ++xvit_begin;

        // Loop through all those matching hits
        for(auto xvit = xvit_begin; xvit != xvs.end() && xvit->a.hit == x.hit; ++xvit){
          const HitOrChan& v = xvit->b;

          // Only allow one bad channel per triplet
          if(!u.hit && !v.hit) continue;

          if(u.hit && v.hit && !CloseTime(tpc,
                                          u.chan, v.chan,
                                          u.hit->PeakTime(),
                                          v.hit->PeakTime())) continue;

          geo::WireIDIntersection ptUV;
          if(!isectUV(u.chan, v.chan, ptUV)) continue;

          if(!CloseSpace(xu.pt, xvit->pt) ||
             !CloseSpace(xu.pt, ptUV) ||
             !CloseSpace(xvit->pt, ptUV)) continue;

          double xavg = HitToXPos(x.chan, x.hit->PeakTime(), tpc);
          int nx = 1;
          if(u.hit){xavg += HitToXPos(u.chan, u.hit->PeakTime(), tpc); ++nx;}
          if(v.hit){xavg += HitToXPos(v.chan, v.hit->PeakTime(), tpc); ++nx;}
          xavg /= nx;

          const XYZ pt{xavg,
              (xu.pt.y + xvit->pt.y + ptUV.y)/3,
              (xu.pt.z + xvit->pt.z + ptUV.z)/3};

          ret.emplace_back(x, u, v, pt);
          ++nxuv;
        } // end for xv
      } // end for xu

      std::cout << tpc << " " << xus.size() << " XUs and " << xvs.size() << " XVs -> " << nxuv << " XUVs" << std::endl;

    } // end for tpc

    std::cout << ret.size() << " XUVs total" << std::endl;

    return ret;
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelTriplet> TripletFinder::TripletsTwoView()
  {
    std::vector<ChannelTriplet> ret;

    for(const auto& it: fX_by_tpc){
      const geo::TPCID& tpc = it.first;

      std::vector<ChannelDoublet> xus = DoubletsXU(tpc);

      for(const ChannelDoublet& xu: xus){
        const HitOrChan& x = xu.a;
        const HitOrChan& u = xu.b;

        double xavg = HitToXPos(x.chan, x.hit->PeakTime(), tpc);
        int nx = 1;
        if(u.hit){xavg += HitToXPos(u.chan, u.hit->PeakTime(), tpc); ++nx;}
        xavg /= nx;

        const XYZ pt{xavg, xu.pt.y, xu.pt.z};

        ret.emplace_back(x, u, HitOrChan{0, 0}, pt);
      } // end for xu
    } // end for tpc

    std::cout << ret.size() << " XUs total" << std::endl;

    return ret;
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelDoublet> TripletFinder::DoubletsXU(geo::TPCID tpc)
  {
    return DoubletHelper(tpc, fX_by_tpc[tpc], fU_by_tpc[tpc], fUbad_by_tpc[tpc]);
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelDoublet> TripletFinder::DoubletsXV(geo::TPCID tpc)
  {
    return DoubletHelper(tpc, fX_by_tpc[tpc], fV_by_tpc[tpc], fVbad_by_tpc[tpc]);
  }

  // -------------------------------------------------------------------------
  std::vector<ChannelDoublet> TripletFinder::
  DoubletHelper(geo::TPCID tpc,
                const std::vector<HitOrChan>& ahits,
                const std::vector<HitOrChan>& bhits,
                const std::vector<raw::ChannelID_t>& bbads) const
  {
    std::vector<ChannelDoublet> ret;

    IntersectionCache isect(tpc);

    for(const HitOrChan& a: ahits){
      // Bad channels are easy because there's no timing constraint
      for(raw::ChannelID_t b: bbads){
        geo::WireIDIntersection pt;
        if(isect(a.chan, b, pt)){
          ret.emplace_back(a, HitOrChan{b, 0}, pt);
        }
      }

      auto b_begin = bhits.begin();
      while(b_begin != bhits.end() &&
            b_begin->hit->PeakTime() < a.hit->PeakTime() &&
            !CloseTime(tpc, b_begin->chan, a.chan,
                       b_begin->hit->PeakTime(),
                       a.hit->PeakTime())) ++b_begin;

      for(auto bit = b_begin; bit != bhits.end(); ++bit){
        const HitOrChan& b = *bit;

        if(b.hit->PeakTime() > a.hit->PeakTime() &&
           !CloseTime(tpc, b.chan, a.chan,
                      b.hit->PeakTime(),
                      a.hit->PeakTime())) break;

        geo::WireIDIntersection pt;
        if(!isect(a.chan, b.chan, pt)) continue;

        ret.emplace_back(a, b, pt);
      } // end for b
    } // end for a

    return ret;
  }
}
