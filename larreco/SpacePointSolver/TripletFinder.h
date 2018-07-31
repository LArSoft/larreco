// Christopher Backhouse - bckhouse@fnal.gov

#ifndef RECO3D_TRIPLETFINDER_H
#define RECO3D_TRIPLETFINDER_H

#include "larcore/Geometry/Geometry.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <map>
#include <vector>

namespace recob{class Hit;}

namespace reco3d
{
  struct HitOrChan
  {
    HitOrChan(raw::ChannelID_t c, const recob::Hit* h) : chan(c), hit(h) {}

    raw::ChannelID_t chan;
    const recob::Hit* hit; // null for bad channel
  };

  struct ChannelDoublet
  {
    ChannelDoublet(HitOrChan a_, HitOrChan b_, geo::WireIDIntersection p)
      : a(a_), b(b_), pt(p)
    {
    }

    HitOrChan a, b;
    geo::WireIDIntersection pt;
  };

  struct XYZ
  {
    double x, y, z;
  };

  struct HitTriplet
  {
    const recob::Hit *x, *u, *v;
    XYZ pt;
  };

  class TripletFinder
  {
  public:
    TripletFinder(const std::vector<art::Ptr<recob::Hit>>& xhits,
                  const std::vector<art::Ptr<recob::Hit>>& uhits,
                  const std::vector<art::Ptr<recob::Hit>>& vhits,
                  const std::vector<raw::ChannelID_t>& ubad,
                  const std::vector<raw::ChannelID_t>& vbad,
                  double distThresh, double distThreshDrift);

    std::vector<HitTriplet> Triplets();
    /// Only search for XU intersections
    std::vector<HitTriplet> TripletsTwoView();

  protected:
    const geo::GeometryCore* geom;
    const detinfo::DetectorProperties* detprop;

    /// Helper for constructor
    void FillHitMap(const std::vector<HitOrChan>& hits,
                    std::map<geo::TPCID, std::vector<HitOrChan>>& out);
    /// Helper for constructor
    void FillBadMap(const std::vector<raw::ChannelID_t>& bads,
                    std::map<geo::TPCID, std::vector<raw::ChannelID_t>>& out);

    double HitToXPos(raw::ChannelID_t chan,
                     double t,
                     geo::TPCID tpc) const;

    bool CloseDrift(double xa, double xb) const;
    bool CloseTime(geo::TPCID tpc,
                   raw::ChannelID_t c1, raw::ChannelID_t c2,
                   double t1, double t2) const;
    bool CloseSpace(geo::WireIDIntersection ra,
                    geo::WireIDIntersection rb) const;

    std::vector<ChannelDoublet> DoubletsXU(geo::TPCID tpc);
    std::vector<ChannelDoublet> DoubletsXV(geo::TPCID tpc);

    std::vector<ChannelDoublet>
    DoubletHelper(geo::TPCID tpc,
                  const std::vector<HitOrChan>& ahits,
                  const std::vector<HitOrChan>& bhits,
                  const std::vector<raw::ChannelID_t>& bbads) const;

    double fDistThresh;
    double fDistThreshDrift;

    std::map<geo::TPCID, std::vector<HitOrChan>> fX_by_tpc;
    std::map<geo::TPCID, std::vector<HitOrChan>> fU_by_tpc;
    std::map<geo::TPCID, std::vector<HitOrChan>> fV_by_tpc;

    std::map<geo::TPCID, std::vector<raw::ChannelID_t>> fUbad_by_tpc;
    std::map<geo::TPCID, std::vector<raw::ChannelID_t>> fVbad_by_tpc;
  };
}

#endif
