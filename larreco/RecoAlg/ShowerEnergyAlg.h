////////////////////////////////////////////////////////////////////////
// Class: ShowerEnergyAlg
// File:  ShowerEnergyAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#ifndef ShowerEnergyAlg_hxx
#define ShowerEnergyAlg_hxx

#include <array>

// Framework
#include "canvas/Persistency/Common/Ptr.h"
namespace fhicl {
  class ParameterSet;
}

// larsoft
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace shower {
  class ShowerEnergyAlg;
}

class shower::ShowerEnergyAlg {
public:
  explicit ShowerEnergyAlg(fhicl::ParameterSet const& pset);

  /// This overload is preferred as it does not rely on the cached
  /// DetectorProperties data member.
  double ShowerEnergy(detinfo::DetectorClocksData const& clockData,
                      detinfo::DetectorPropertiesData const& detProp,
                      std::vector<art::Ptr<recob::Hit>> const& hits,
                      geo::PlaneID::PlaneID_t plane) const;

private:
  struct LinearFunction {
    double gradient;
    double intercept;
    double energy_from(double const charge) const noexcept { return charge * gradient + intercept; }
  };
  std::array<LinearFunction, 3> const fLinearFunctions;
};

#endif
