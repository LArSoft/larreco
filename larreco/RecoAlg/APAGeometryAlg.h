/////////////////////////////////////////////////////////////////
//  \file APAGeometryAlg.h
//  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////
#ifndef APAGeometryALG_H
#define APAGeometryALG_H

#include <stdint.h>
#include <vector>

#include "TVector3.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace fhicl {
  class ParameterSet;
}

namespace apa {

  // each APA has 4 separate views
  typedef enum _apa_plane_proj {
    kU,  ///< U view on both sides of the APA
    kV,  ///< V view on both sides of the APA
    kZ0, ///< Z view on the smaller-x side of the APA
    kZ1, ///< Z view on the larger-x side of the APA
    kUnknown
  } APAView_t;

  //---------------------------------------------------------------
  class APAGeometryAlg {
  public:
    APAGeometryAlg(fhicl::ParameterSet const& pset);
    APAGeometryAlg();

    void reconfigure(fhicl::ParameterSet const& p);

    void Init(); ///< Initialize some chanel numbers to speed up other methods

    bool APAChannelsIntersect(uint32_t chan1,
                              uint32_t chan2,
                              std::vector<geo::WireIDIntersection>& IntersectVector) const;
    ///< If the channels intersect, get all intersections

    bool LineSegChanIntersect(TVector3 xyzStart,
                              TVector3 xyzEnd,
                              uint32_t chan,
                              std::vector<geo::WireID>& widsCrossed,
                              bool ExtendLine) const;
    ///< If a line given by start/end points intersects a channel

    std::vector<geo::WireID> ChanSegsPerSide(uint32_t chan, unsigned int side) const;
    std::vector<geo::WireID> ChanSegsPerSide(std::vector<geo::WireID> wids,
                                             unsigned int side) const;

    std::vector<double> ThreeChanPos(uint32_t u, uint32_t v, uint32_t z) const;
    ///< Find the center of the 3 intersections, choose best if multiple

    geo::WireID NearestWireIDOnChan(const double WorldLoc[3],
                                    uint32_t chan,
                                    unsigned int const plane,
                                    unsigned int const tpc = 0,
                                    unsigned int const cstat = 0) const;

    unsigned int ChannelToAPA(
      uint32_t chan) const; ///< Get number of the APA containing the given channel
    void ChannelToAPA(uint32_t chan, unsigned int& apa, unsigned int& cryo) const;
    APAView_t APAView(uint32_t chan) const; ///< Get which of the 4 APA views the channel is in
    unsigned int ChannelsInView(geo::View_t geoview) const;
    uint32_t FirstChannelInView(geo::View_t geoview, unsigned int apa, unsigned int cryo) const;
    uint32_t FirstChannelInView(geo::View_t geoview, uint32_t chan) const;
    uint32_t FirstChannelInView(uint32_t chan) const;
    unsigned int ChannelsInAPAView(APAView_t apaview) const;
    unsigned int ChannelsPerAPA() const { return fChannelsPerAPA; };

  private:
    art::ServiceHandle<geo::Geometry const> fGeom; // handle to geometry service

    unsigned int fChannelsPerAPA; ///< All APAs have this same number of channels
    unsigned int fAPAsPerCryo;

    // channel boundaries to avoid calling ChannelToWire repetitively
    uint32_t fFirstU;
    uint32_t fLastU;
    uint32_t fFirstV;
    uint32_t fLastV;
    uint32_t fFirstZ0;
    uint32_t fLastZ0;
    uint32_t fFirstZ1;
    uint32_t fLastZ1;

    double fChannelRange[2]; // for each induction view: U=0, V=1

  }; // class APAGeometryAlg

} // namespace apa

#endif // ifndef APAGeometryALG_H
