#ifndef TRACKCALORIMETRYALG_H
#define TRACKCALORIMETRYALG_H
/*!
 * Title:   Track Calorimetry Algorithim Class
 * Author:  Wes Ketchum (wketchum@lanl.gov), based on code the Calorimetry_module
 *
 * Description: Algorithm that produces a calorimetry object given a track
 * Input:       recob::Track, Assn<recob::Spacepoint,recob::Track>, Assn<recob::Hit,recob::Track>
 * Output:      anab::Calorimetry, (and Assn<anab::Calorimetry,recob::Track>)
*/
#include <iostream>

#include "larcorealg/CoreUtils/ProviderPack.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TVector3.h"
#include <set>

namespace anab {
  class Calorimetry;
}
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
  class LArProperties;
}
namespace fhicl {
  class ParameterSet;
}
namespace geo {
  struct PlaneID;
  struct WireID;
}
namespace recob {
  class Hit;
  class Track;
}

namespace calo {
  class TrackCalorimetryAlg;
}

class calo::TrackCalorimetryAlg {
public:
  using Providers_t = lar::ProviderPack<geo::GeometryCore, detinfo::LArProperties>;

  TrackCalorimetryAlg(fhicl::ParameterSet const& p);

  void ExtractCalorimetry(detinfo::DetectorClocksData const& clock_data,
                          detinfo::DetectorPropertiesData const& det_prop,
                          std::vector<recob::Track> const&,
                          std::vector<recob::Hit> const&,
                          std::vector<std::vector<size_t>> const&,
                          std::vector<anab::Calorimetry>&,
                          std::vector<size_t>&,
                          Providers_t providers);

private:
  CalorimetryAlg caloAlg;
  unsigned int fNHitsToDetermineStart;

  struct HitProperties {
    HitProperties() = default;
    HitProperties(float q, float dqdx, float dedx, float p, TVector3 pos, float pf)
      : charge(q), dQdx(dqdx), dEdx(dedx), pitch(p), xyz(pos), path_fraction(pf)
    {}
    float charge;
    float dQdx;
    float dEdx;
    float pitch;
    TVector3 xyz;
    float path_fraction;
    void
    Print() const
    {
      std::cout << "\tCharge " << charge << "  dQdx " << dQdx << "  dEdx " << dEdx << "  pitch "
                << pitch << "  (x,y,z) (" << xyz.X() << "," << xyz.Y() << "," << xyz.Z() << ")"
                << " path_fraction " << path_fraction << std::endl;
    }
  };
  struct HitPropertySorter {
    bool
    operator()(HitProperties const& i, HitProperties const& j) const
    {
      return i.path_fraction < j.path_fraction;
    }
  };

  using HitPropertiesMultiset_t = std::multiset<HitProperties, HitPropertySorter>;

  void
  ClearInternalVectors()
  {}
  void
  ReserveInternalVectors(size_t s)
  {}

  std::vector<float> CreatePathLengthFractionVector(recob::Track const& track);

  void AnalyzeHit(detinfo::DetectorClocksData const&,
                  detinfo::DetectorPropertiesData const&,
                  recob::Hit const&,
                  recob::Track const&,
                  std::vector<std::pair<geo::WireID, float>> const&,
                  std::vector<float> const&,
                  HitPropertiesMultiset_t&,
                  geo::GeometryCore const&);

  bool IsInvertedTrack(HitPropertiesMultiset_t const&);

  void MakeCalorimetryObject(HitPropertiesMultiset_t const& hpm,
                             recob::Track const& track,
                             size_t const& i_track,
                             std::vector<anab::Calorimetry>& caloVector,
                             std::vector<size_t>& assnTrackCaloVector,
                             geo::PlaneID const& planeID);

  void PrintHitPropertiesMultiset(HitPropertiesMultiset_t const& hpm);
};

#endif
